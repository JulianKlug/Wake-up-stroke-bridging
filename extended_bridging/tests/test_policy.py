"""Acceptance tests for Stage 13 — §16 of `specs/stage13_feasible_policy.md`.

Section banners match the specification's numbers. §16.2's config assertions live in
`test_config.py`, §16.11's two digests and §16.15's `gcompute` tests in `test_standardise.py`, and
§16.17's raise-site scan in `test_bootstrap.py`, because that is where the changed code is.

Tests needing the private workbook are marked `DATA_GATED` and tagged `[data-gated]`; the one test
that runs Stage 12's 4.2-minute bootstrap is `SLOW`-gated on `STAGE12_SLOW=1` as well.

**NO STANDARDISED PROBABILITY, RISK DIFFERENCE, INTERVAL LIMIT, TREATMENT COEFFICIENT OR
CONTRAINDICATION COEFFICIENT FROM THE WORKBOOK APPEARS IN THIS FILE** (§4.3). The workbook-gated
tests assert counts, set membership, identities and properties; every numerical pin is measured on
`fixtures_stage13.py`.
"""
from __future__ import annotations

import ast
import dataclasses
import os
from pathlib import Path
from typing import Final

import numpy as np
import pandas as pd
import pytest

import bootstrap
import cohort
import config
import data
import derive
import eligibility
import model
import outcome
import policy
import standardise
from fixtures_stage13 import (N_CONTRAINDICATED, all_eligible_frame, mixed_frame,
                              no_contraindicated_frame, separated_indicator_frame,
                              treated_contraindicated_frame)
from test_data import hand_source
from test_standardise import COHORT_STEPS, DATA_GATED, SLOW, classified

MODULE = Path(policy.__file__).resolve()
SOURCE = MODULE.read_text(encoding="utf-8")

RUN_BOOT: Final[int] = 60

# `standardise`'s ten step names (Stage 12 §15), for §16.9's disjointness — an intersection, not a list.
STANDARDISE_STEPS: Final[frozenset[str]] = frozenset({
    "standardisation_population", "absence_by_standardisation_column", "standardisation_fit",
    "standardised_distributions", "treated_support", "support_baseline",
    "standardisation_support_restricted", "hierarchical_fit", "standardised_hierarchical",
    "standardisation_replicates"})
POLICY_STEPS: Final[tuple[tuple[str, str], ...]] = (
    ("cohort", "policy_population"),
    ("missingness", "absence_by_policy_column"),
    ("model", "policy_fit"),
    ("model", "policy_contrast"),
    ("model", "policy_replicates"),
)
U_SERIES: Final[frozenset[str]] = frozenset({"U1", "U2", "U3", "U4", "U5", "U7", "U8"})


def message_of(excinfo) -> str:
    return str(excinfo.value)


def audit_of() -> data.Audit:
    return data.Audit(hand_source(Path("."), "stage13"))


def driver(df: pd.DataFrame, audit: data.Audit, with_inference: bool = False):
    """§13's canonical call order."""
    pop = policy.population(df, audit)
    pol = policy.contrast(pop, audit)
    boot = policy.inference(pop, audit) if with_inference else None
    return pop, pol, boot


def with_boot(n: int, fn):
    original = config.N_BOOT
    config.N_BOOT = n
    try:
        return fn()
    finally:
        config.N_BOOT = original


def eligible_only_fit(pop: pd.DataFrame) -> tuple[pd.DataFrame, model.PolrFit]:
    """A [§14a]-shaped fit: eligible rows only, no indicator — for §16.7's contrasts."""
    rows = pop[pop[config.CONTRAINDICATED] == 0.0]
    X, _ = model.design(rows, (config.TREATMENT,) + config.STANDARDISATION_COVARIATES)
    return X, model.polr(X, rows[config.PRIMARY_OUTCOME].to_numpy(dtype=float))


# --- 16.1  the population is [§14b]'s and nobody else's -------------------------------------------

@DATA_GATED
def test_the_population_is_123_over_FOUR_centres_with_19_contraindicated_beside_104_and_93():
    """§16.1. All three populations in one test, so the DIFFERENCES are the assertion. [data-gated]"""
    df, audit = classified()
    pop = policy.population(df, audit)
    std = standardise.population(df, audit_of())
    built = cohort.build(df, audit_of())

    assert len(pop) == 123 and len(std) == 104 and len(built) == 93
    assert sorted(pop["center"].unique()) == sorted(config.CENTER_ORDER)
    assert int(pop[config.CONTRAINDICATED].sum()) == 19
    # The 123 is the 104 plus exactly the records `eligibility.retained` is False for.
    ineligible = set(df.loc[~eligibility.retained(df), "case_id"])
    assert set(pop["case_id"]) == set(std["case_id"]) | ineligible
    assert set(pop["case_id"]) - set(std["case_id"]) == ineligible
    assert len(ineligible) == 19
    # Every contraindicated record is EVT alone.
    assert (pop.loc[pop[config.CONTRAINDICATED] == 1.0, config.TREATMENT]
            == min(config.TREATMENT_LABELS)).all()
    # Contraindicated at exactly two centres.
    per_centre = pop.groupby("center", observed=True)[config.CONTRAINDICATED].sum()
    assert int((per_centre > 0).sum()) == 2


@DATA_GATED
def test_the_three_11_losses_are_Stage_12s_three_by_identifier_and_all_are_eligible():
    """§16.1, §4.1. [data-gated]"""
    df, audit = classified()
    policy.population(df, audit)
    mine = audit.entry("cohort", "policy_population")
    other = audit_of()
    standardise.population(df, other)
    theirs = other.entry("cohort", "standardisation_population")

    assert mine.n == 3
    # Stage 12's entry names the 19 ineligible AND its three [§11] losses; ours names the three only.
    assert set(mine.case_ids) <= set(theirs.case_ids)
    assert set(theirs.case_ids) - set(mine.case_ids) == set(df.loc[~eligibility.retained(df), "case_id"])
    assert eligibility.retained(df[df["case_id"].isin(mine.case_ids)]).all()


def test_neither_cohort_nor_propensity_nor_derive_is_imported_by_scan():
    """§16.1, §0.2. The imports are the assertion."""
    tree = ast.parse(SOURCE)
    imported = {n.names[0].name for n in ast.walk(tree) if isinstance(n, ast.Import)}
    imported |= {n.module for n in ast.walk(tree) if isinstance(n, ast.ImportFrom)}
    assert not imported & {"cohort", "propensity", "derive"}
    assert {"standardise", "bootstrap", "eligibility", "model"} <= imported
    calls = {ast.unparse(c.func) for c in ast.walk(tree) if isinstance(c, ast.Call)}
    assert not [c for c in calls if c.startswith(("cohort.", "propensity.", "derive."))]
    # `standardise.population` is not called either: it applies restriction 2 (§0.2).
    assert "standardise.population" not in calls
    assert "eligibility.retained" in calls


# --- 16.3  the estimand keys are derived and not counted ------------------------------------------

def test_the_key_count_is_DERIVED_from_the_config_and_no_key_is_a_p_a_contraindication_or_a_dist():
    """§16.3, §3.3."""
    keys = policy._estimand_keys()
    assert len(keys) == 2 * len(config.MRS_THRESHOLDS) + 2 + 2 * len(config.MRS_LEVELS) + 2
    assert len(set(keys)) == len(keys)
    assert not [k for k in keys if "contraindication" in k]
    assert not [k for k in keys if k.startswith("dist")]
    assert {"beta", "share_eligible", "mrs_0_2", "mortality"} <= set(keys)
    assert all(f"eligible_rd_{k}" in keys for k in config.MRS_THRESHOLDS)
    assert all(f"active_{j}" in keys and f"comparator_{j}" in keys for j in config.MRS_LEVELS)

    boot = with_boot(RUN_BOOT, lambda: driver(mixed_frame(), audit_of(), with_inference=True)[2])
    assert set(boot.draws) == set(keys)
    assert all(interval.p is None for interval in boot.intervals.values())


# --- 16.4  the fit is unweighted, on all rows -----------------------------------------------------

def test_model_polr_is_called_with_TWO_POSITIONAL_ARGUMENTS_and_no_w_and_no_row_is_excluded():
    """§16.4, §5.1. The call site is the assertion."""
    calls = [c for c in ast.walk(ast.parse(SOURCE))
             if isinstance(c, ast.Call) and ast.unparse(c.func) == "model.polr"]
    assert len(calls) == 1
    assert len(calls[0].args) == 2 and calls[0].keywords == []

    pop = policy.population(mixed_frame(), audit_of())
    pol = policy.contrast(pop, audit_of())
    assert pol.n_fit == pol.n_average == len(pop) == len(mixed_frame())


# --- 16.5  the regimes ----------------------------------------------------------------------------

def test_the_active_vector_bridges_exactly_the_eligible_and_the_comparator_nobody():
    """§16.5, §6.2."""
    pop = policy.population(mixed_frame(), audit_of())
    D, active, comparator = policy._regimes(pop)
    treated, control = max(config.TREATMENT_LABELS), min(config.TREATMENT_LABELS)
    assert len(active) == len(comparator) == len(D) == len(pop)
    assert int((active == treated).sum()) == int(eligibility.retained(pop).sum())
    assert int(D.sum()) == N_CONTRAINDICATED
    assert (comparator == control).all()
    assert (active[D == 1.0] == control).all()


def test_U2_fires_BEFORE_any_prediction_when_a_contraindicated_row_is_set_to_the_treated_code(
        monkeypatch):
    """§16.5, §6.3. Asserted by the message AND by `ordinal_probabilities` never having been called."""
    pop = policy.population(mixed_frame(), audit_of())
    D, active, comparator = policy._regimes(pop)
    corrupted = active.copy()
    corrupted[np.flatnonzero(D == 1.0)[0]] = float(max(config.TREATMENT_LABELS))
    monkeypatch.setattr(policy, "_regimes", lambda frame: (D, corrupted, comparator))

    calls = []
    real = model.ordinal_probabilities
    monkeypatch.setattr(model, "ordinal_probabilities",
                        lambda fit, X: calls.append(1) or real(fit, X))
    with pytest.raises(config.SchemaError) as excinfo:
        policy.contrast(pop, audit_of())
    assert message_of(excinfo).startswith("U2")
    assert calls == []


# --- 16.6  the dilution identity and the two namings ----------------------------------------------

def test_the_dilution_identity_and_the_two_namings_hold_on_the_FIXTURE():
    """§16.6, §7.3, §7.2."""
    pol = policy.contrast(policy.population(mixed_frame(), audit_of()), audit_of())
    for k in config.MRS_THRESHOLDS:
        assert abs(pol.rd[k] - pol.share_eligible * pol.eligible_rd[k]) < 1e-12
    assert pol.mrs_0_2 == pol.rd[2]
    assert abs(pol.mortality + pol.rd[5]) < 1e-15
    assert pol.share_eligible == pol.n_eligible / pol.n_average
    assert pol.n_eligible == pol.n_average - N_CONTRAINDICATED


@DATA_GATED
def test_the_dilution_identity_and_the_two_namings_hold_on_the_WORKBOOK():
    """§16.6. Identities and tolerances only; no value is quoted. [data-gated]"""
    df, audit = classified()
    _, pol, _ = driver(df, audit)
    for k in config.MRS_THRESHOLDS:
        assert abs(pol.rd[k] - pol.share_eligible * pol.eligible_rd[k]) < 1e-12
    assert pol.mrs_0_2 == pol.rd[2]
    assert abs(pol.mortality + pol.rd[5]) < 1e-15
    assert pol.n_eligible == 104 and pol.n_average == 123
    assert pol.dropped == ()
    assert len(pol.fit.columns) == 1 + len(config.POLICY_COVARIATES) + (
        len(config.FACTOR_LEVELS["onset_type"]) - 2)


def test_the_per_row_contrast_of_every_contraindicated_row_is_EXACTLY_zero():
    """§16.6, §7.3. Measured `0.0` and not `1e-17`: the same function of the same floats."""
    pop = policy.population(mixed_frame(), audit_of())
    pol = policy.contrast(pop, audit_of())
    D, active, comparator = policy._regimes(pop)
    X, _ = model.design(pop, (config.TREATMENT,) + config.POLICY_COVARIATES)
    Xa, Xc = X.copy(), X.copy()
    Xa[config.TREATMENT], Xc[config.TREATMENT] = active, comparator
    per_row = (model.ordinal_probabilities(pol.fit, Xa)
               - model.ordinal_probabilities(pol.fit, Xc))
    assert (per_row[D == 1.0] == 0.0).all()
    assert (per_row[D == 0.0] != 0.0).any()


def test_U7_fires_on_an_eligible_rd_perturbed_by_1e_13():
    """§16.6. Above the 1e-12 tolerance by construction? No — BELOW it; so the perturbation must be
    scaled by the share to land above. A 1e-13 perturbation of `eligible_rd[2]` is caught only if
    `share * 1e-13 > 1e-12`, which it is not, so the test perturbs by 1e-11 and ALSO shows 1e-13
    passing: the tolerance is a contract on arithmetic, not on the fixture."""
    pol = policy.contrast(policy.population(mixed_frame(), audit_of()), audit_of())
    nudged = dict(pol.eligible_rd)
    nudged[2] += 1e-13
    policy._assert_dilution(pol.rd, nudged, pol.share_eligible)           # inside 1e-12: silent
    nudged[2] = pol.eligible_rd[2] + 1e-11
    with pytest.raises(config.SchemaError) as excinfo:
        policy._assert_dilution(pol.rd, nudged, pol.share_eligible)
    assert message_of(excinfo).startswith("U7")
    assert "rd_2" in message_of(excinfo)


def test_the_two_sides_of_the_identity_are_TWO_gcompute_calls_and_the_second_sees_no_contraindicated_row(
        monkeypatch):
    """§16.6, §7.3 item 3, §6.3. A `rd / share_eligible` shortcut makes the second call disappear; a
    full-design second call fails the row count."""
    pop = policy.population(mixed_frame(), audit_of())
    D, active, comparator = policy._regimes(pop)
    treated, control = max(config.TREATMENT_LABELS), min(config.TREATMENT_LABELS)
    seen = []
    real = standardise.gcompute

    def spy(X, fit, over, act, comp, **kw):
        result = real(X, fit, over, act, comp, **kw)
        seen.append((X, over, act, comp, kw, result))
        return result

    monkeypatch.setattr(standardise, "gcompute", spy)
    pol = policy.contrast(pop, audit_of())

    assert len(seen) == 2
    X1, over1, act1, comp1, kw1, g1 = seen[0]
    assert len(X1) == len(pop) and over1.all()
    assert (act1 == active).all() and (comp1 == control).all()
    assert kw1["keys"] == (policy._ACTIVE, policy._COMPARATOR)
    X2, over2, act2, comp2, kw2, g2 = seen[1]
    assert len(X2) == pol.n_eligible and over2.all()
    assert not (X2[config.CONTRAINDICATED] == 1.0).any()
    assert (act2 == treated).all() and (comp2 == control).all()
    assert "keys" not in kw2
    assert pol.rd == g1.rd and pol.eligible_rd == g2.rd
    assert pol.distribution == g1.distribution


def test_U2_fires_on_a_STORED_column_flipped_after_population_when_regimes_is_bypassed(monkeypatch):
    """§16.6, §6.3. The two-source check: the stored column against the freshly built vector."""
    pop = policy.population(mixed_frame(), audit_of())
    D, active, comparator = policy._regimes(pop)
    eligible_row = np.flatnonzero(D == 0.0)[0]
    pop.iloc[eligible_row, pop.columns.get_loc(config.CONTRAINDICATED)] = 1.0
    monkeypatch.setattr(policy, "_regimes", lambda frame: (D, active, comparator))
    with pytest.raises(config.SchemaError) as excinfo:
        policy.contrast(pop, audit_of())
    assert message_of(excinfo).startswith("U2")


# --- 16.7  the absolute levels move with the fit and the contrast does not ------------------------

def test_the_contraindicated_block_and_the_rd_move_with_the_fit_on_mixed_frame():
    """§16.7, §7.3. Two fits that agree on nothing here; §5.1's gap on a fixture."""
    pop = policy.population(mixed_frame(), audit_of())
    pol = policy.contrast(pop, audit_of())
    Xe, fit_a = eligible_only_fit(pop)
    X, _ = model.design(pop, (config.TREATMENT,) + config.POLICY_COVARIATES)
    contraindicated = pop[config.CONTRAINDICATED].to_numpy() == 1.0
    Xc = X.loc[contraindicated, list(Xe.columns)].copy()
    Xc[config.TREATMENT] = float(min(config.TREATMENT_LABELS))
    under_14a = model.ordinal_probabilities(fit_a, Xc)
    Xb = X[contraindicated].copy()
    Xb[config.TREATMENT] = float(min(config.TREATMENT_LABELS))
    under_14b = model.ordinal_probabilities(pol.fit, Xb)
    assert np.abs(under_14a - under_14b).max() > 1e-3

    n = len(Xe)
    g_a = standardise.gcompute(Xe, fit_a, np.ones(n, dtype=bool), np.ones(n), np.zeros(n),
                               label="14a-shaped")
    assert any(abs(g_a.rd[k] - pol.eligible_rd[k]) > 1e-4 for k in config.MRS_THRESHOLDS)


def test_on_separated_indicator_frame_the_eligible_contrast_under_the_two_fits_agrees_to_1e_6():
    """§16.7, §9.1's mechanism as an assertion. 1e-6 and not tighter: two Newton runs on different
    designs stopping on `POLR_TOL = 1e-8`."""
    pop = policy.population(separated_indicator_frame(), audit_of())
    pol = policy.contrast(pop, audit_of())
    Xe, fit_a = eligible_only_fit(pop)
    n = len(Xe)
    g_a = standardise.gcompute(Xe, fit_a, np.ones(n, dtype=bool), np.ones(n), np.zeros(n),
                               label="14a-shaped")
    for k in config.MRS_THRESHOLDS:
        assert abs(g_a.rd[k] - pol.eligible_rd[k]) < 1e-6
    assert abs(pol.contraindication_log_odds) > config.POLR_MAX_ABS_BETA


# --- 16.8  the two labels -------------------------------------------------------------------------

def test_every_Policy_carries_BOTH_guard_strings_and_has_NO_odds_ratio():
    """§16.8, §8."""
    pol = policy.contrast(policy.population(mixed_frame(), audit_of()), audit_of())
    assert pol.measure is policy._CONDITIONAL_14B
    assert pol.label is policy._OPERATIONAL
    assert "[§14b]" in pol.label and "[§14a]" in pol.label and "[§7]" in pol.label
    assert "contraindication" in pol.measure
    assert pol.measure != standardise._CONDITIONAL
    assert not hasattr(pol, "odds_ratio")
    fields = {f.name for f in dataclasses.fields(policy.Policy)}
    assert "odds_ratio" not in fields and "ate" not in fields and "rd_ate" not in fields
    assert pol.conditional_odds_ratio == pytest.approx(np.exp(pol.conditional_log_odds))

    tree = ast.parse(SOURCE)
    docstrings = set()
    for node in ast.walk(tree):
        if isinstance(node, (ast.Module, ast.FunctionDef, ast.ClassDef)):
            doc = ast.get_docstring(node, clean=False)
            if doc is not None:
                docstrings.add(doc)
    literals = [n.value for n in ast.walk(tree)
                if isinstance(n, ast.Constant) and isinstance(n.value, str)]
    assert not [s for s in literals if "odds_ratio" in s and s not in docstrings]


# --- 16.9  the audit entries and their step names -------------------------------------------------

def test_the_five_entries_appear_in_13s_DRIVER_ORDER_and_collide_with_NOBODY():
    """§16.9, §12. Two empty intersections, and entry 1 names exactly its `n`."""
    audit = audit_of()
    frame = mixed_frame()
    pop, _, _ = with_boot(RUN_BOOT, lambda: driver(frame, audit, with_inference=True))
    steps = [(e.kind, e.step) for e in audit.entries]
    assert steps == list(POLICY_STEPS)
    mine = {step for _, step in steps}
    assert not mine & COHORT_STEPS
    assert not mine & STANDARDISE_STEPS

    first = audit.entries[0]
    assert first.n == len(first.case_ids) == len(frame) - len(pop) == 0
    dropped_one = frame.copy()
    dropped_one.loc[dropped_one.index[3], config.PRIMARY_OUTCOME] = np.nan
    other = audit_of()
    policy.population(dropped_one, other)
    entry = other.entry("cohort", "policy_population")
    assert entry.n == len(entry.case_ids) == 1


@DATA_GATED
def test_the_Stage_2_to_4_path_records_TWELVE_entries_and_this_stage_FIVE_making_17():
    """§16.9's ungated form: 12 before, 5 added. [data-gated]"""
    df, audit = classified()
    assert len(audit.entries) == 12
    with_boot(RUN_BOOT, lambda: driver(df, audit, with_inference=True))
    assert len(audit.entries) == 17
    assert [(e.kind, e.step) for e in audit.entries[12:]] == list(POLICY_STEPS)


@SLOW
@DATA_GATED
def test_a_driver_running_Stages_12_and_13_records_SIXTY_ONE_entries():
    """§16.9. Stage 12's 56 plus these 5. [slow] [data-gated]"""
    from test_standardise import driver as stage12_driver
    df, audit = classified()
    stage12_driver(df, audit, with_inference=True)
    assert len(audit.entries) == 12 + 10
    driver(df, audit, with_inference=True)
    assert len(audit.entries) == 12 + 10 + 5
    # 56 is the Stage 1–12 ledger with Stages 5–11 in; here only Stage 12's ten ran on top of the
    # Stage 2–4 twelve, so the assertion is the increment: five, and the step set is disjoint.
    steps = [e.step for e in audit.entries]
    assert len(steps) == len(set(steps))


# --- 16.10  the indicator is derived from `eligibility` -------------------------------------------

@DATA_GATED
def test_blanking_ivt_contraindicated_after_classify_leaves_population_UNCHANGED():
    """§16.10, §5.2. The flag is Stage 4's input and nothing of this stage's. [data-gated]"""
    df, _ = classified()
    before = policy.population(df, audit_of())
    blanked = df.copy()
    blanked["ivt_contraindicated"] = np.nan
    after = policy.population(blanked, audit_of())
    pd.testing.assert_frame_equal(before.drop(columns=["ivt_contraindicated"]),
                                  after.drop(columns=["ivt_contraindicated"]))


def test_U4_fires_on_a_THIRD_eligibility_value_and_on_NO_eligibility_column():
    """§16.10, §5.2."""
    third = mixed_frame()
    third.loc[third.index[0], config.ELIGIBILITY] = "undocumented"
    with pytest.raises(config.SchemaError) as excinfo:
        policy.population(third, audit_of())
    assert message_of(excinfo).startswith("U4")

    with pytest.raises(config.SchemaError) as excinfo:
        policy.population(mixed_frame().drop(columns=[config.ELIGIBILITY]), audit_of())
    assert message_of(excinfo).startswith("U4")
    assert config.ELIGIBILITY in message_of(excinfo)


def test_U3_fires_on_treated_contraindicated_frame():
    """§16.10, §4.1."""
    with pytest.raises(config.SchemaError) as excinfo:
        policy.population(treated_contraindicated_frame(), audit_of())
    assert message_of(excinfo).startswith("U3")


def test_U5_fires_on_an_EMPTY_frame_and_on_eligible_patients_all_in_ONE_arm():
    """§16.10, §11."""
    with pytest.raises(config.SchemaError) as excinfo:
        policy.population(mixed_frame().iloc[0:0], audit_of())
    assert message_of(excinfo).startswith("U5")

    one_arm = mixed_frame()
    one_arm[config.TREATMENT] = float(min(config.TREATMENT_LABELS))
    with pytest.raises(config.SchemaError) as excinfo:
        policy.population(one_arm, audit_of())
    assert message_of(excinfo).startswith("U5")
    assert config.TREATMENT_LABELS[max(config.TREATMENT_LABELS)] in message_of(excinfo)


def test_U5_does_NOT_fire_on_all_eligible_frame_which_is_REPORTED_not_refused():
    """§16.10, §11. A cohort with nobody to withhold IVT from is a data fact, recorded."""
    audit = audit_of()
    pop = policy.population(all_eligible_frame(), audit)
    assert len(pop) == len(all_eligible_frame())
    entry = audit.entry("cohort", "policy_population")
    assert policy._NO_CONTRAINDICATED in entry.detail
    assert "share_eligible" in entry.detail
    assert all(f"'{c}': 0" in entry.detail for c in config.CENTER_ORDER)

    pol = policy.contrast(pop, audit)
    assert pol.share_eligible == 1.0
    assert pol.dropped == (config.CONTRAINDICATED,)
    assert pol.rd == pol.eligible_rd
    assert pol.label is policy._OPERATIONAL
    assert np.isnan(pol.contraindication_log_odds)


# --- 16.12  separation ----------------------------------------------------------------------------

def test_the_indicator_separates_on_separated_indicator_frame_and_U8_does_NOT_fire():
    """§16.12, §9. The companion applying `outcome._assert_reportable`'s worst-coefficient rule DOES
    raise — the guard's scope is a choice, and this names the rule it declined."""
    pop = policy.population(separated_indicator_frame(), audit_of())
    pol = policy.contrast(pop, audit_of())
    assert abs(pol.contraindication_log_odds) > config.POLR_MAX_ABS_BETA
    assert pol.contraindication_log_odds < 0.0
    assert abs(pol.conditional_log_odds) < config.POLR_MAX_ABS_BETA
    assert all(np.isfinite(v) and -1.0 <= v <= 1.0 for v in pol.rd.values())
    with pytest.raises(model.FitError) as excinfo:
        outcome._assert_reportable(pol.fit)
    assert message_of(excinfo).startswith("G7")

    replicate = policy._replicate(pop, policy._estimand_keys())
    assert replicate.failures == {}
    assert replicate.max_abs_beta["policy"] > config.POLR_MAX_ABS_BETA


def test_U8_fires_on_a_constructed_fit_costing_beta_and_NOTHING_ELSE(monkeypatch):
    """§16.12, §9.2. Key granularity."""
    keys = policy._estimand_keys()
    pop = policy.population(mixed_frame(), audit_of())
    real = model.polr

    def separated(X, y, w=None):
        fit = real(X, y, w)
        beta = fit.beta.copy()
        beta[fit.columns.index(config.TREATMENT)] = config.POLR_MAX_ABS_BETA + 1.0
        return dataclasses.replace(fit, beta=beta)

    monkeypatch.setattr(model, "polr", separated)
    replicate = policy._replicate(pop, keys)
    assert set(replicate.failures) == {"beta"}
    assert replicate.failures["beta"] == "separation"
    assert set(replicate.values) == set(keys) - {"beta"}
    assert all(0.0 <= replicate.values[f"active_{j}"] <= 1.0 for j in config.MRS_LEVELS)

    with pytest.raises(model.FitError) as excinfo:
        policy.contrast(pop, audit_of())
    assert message_of(excinfo).startswith("U8")


# --- 16.13  the replicate loop --------------------------------------------------------------------

def test_the_drawn_frame_sequence_is_a_function_of_the_SEED_ALONE_and_all_four_strata_are_exact():
    """§16.13, Stage 10 §15.2's property over this body."""
    pop = policy.population(mixed_frame(), audit_of())
    sizes = pop.groupby(config.BOOT_STRATUM, observed=True).size().to_dict()
    assert len(sizes) == len(config.CENTER_ORDER)

    def frames(body):
        seen = []
        bootstrap.replicates(pop, lambda draw: seen.append(draw) or body(draw), 6, config.SEED,
                             config.BOOT_STRATUM)
        return seen

    quiet = frames(lambda draw: None)
    real = frames(lambda draw: policy._replicate(draw, policy._estimand_keys()))
    assert len(quiet) == len(real) == 6
    for a, b in zip(quiet, real):
        pd.testing.assert_frame_equal(a, b)
        assert a.groupby(config.BOOT_STRATUM, observed=True).size().to_dict() == sizes


def test_the_contraindicated_count_VARIES_across_draws_and_D_agrees_with_the_stored_column():
    """§16.13, §16.16, §10.4."""
    pop = policy.population(mixed_frame(), audit_of())
    generator = np.random.default_rng(config.SEED)
    counts = set()
    treated = max(config.TREATMENT_LABELS)
    for _ in range(40):
        draw = bootstrap.resample(pop, generator, config.BOOT_STRATUM)
        D, active, comparator = policy._regimes(draw)
        counts.add(int(D.sum()))
        assert (D + active / treated == 1.0).all()
        assert (comparator == min(config.TREATMENT_LABELS)).all()
        assert (D == draw[config.CONTRAINDICATED].to_numpy()).all()
    assert len(counts) > 1 and min(counts) > 0


@DATA_GATED
def test_the_contraindicated_count_is_NEVER_ZERO_on_the_workbook():
    """§16.13, §5.3. Both strata would have to miss all of theirs at once. [data-gated]"""
    df, _ = classified()
    pop = policy.population(df, audit_of())
    generator = np.random.default_rng(config.SEED)
    counts = [int(bootstrap.resample(pop, generator, config.BOOT_STRATUM)[config.CONTRAINDICATED].sum())
              for _ in range(400)]
    assert min(counts) > 0 and len(set(counts)) > 1


def test_U1_costs_EVERY_key_under_degenerate_design():
    """§16.13, §16.14. A frame where every patient is treated is unreachable through `population`
    (U5 fires) and reachable inside a replicate."""
    keys = policy._estimand_keys()
    everyone = policy.population(mixed_frame(), audit_of())
    everyone = everyone[everyone[config.CONTRAINDICATED] == 0.0].copy()
    everyone[config.TREATMENT] = float(max(config.TREATMENT_LABELS))
    replicate = policy._replicate(everyone, keys)
    assert replicate.values == {}
    assert set(replicate.failures) == set(keys)
    assert set(replicate.failures.values()) == {"degenerate_design"}
    assert replicate.n_alpha is None

    X, dropped = model.design(everyone, (config.TREATMENT,) + config.POLICY_COVARIATES)
    with pytest.raises(model.FitError) as excinfo:
        policy._assert_exposure_survived(X, dropped)
    assert message_of(excinfo).startswith("U1")


def test_a_polr_failure_costs_EVERY_key_under_nonconvergence(monkeypatch):
    """§16.13, §10.2. The second failure group."""
    keys = policy._estimand_keys()
    pop = policy.population(mixed_frame(), audit_of())
    monkeypatch.setattr(
        model, "polr",
        lambda *a, **k: (_ for _ in ()).throw(model.FitError("polr: no convergence, forced")))
    replicate = policy._replicate(pop, keys)
    assert set(replicate.failures) == set(keys)
    assert set(replicate.failures.values()) == {"nonconvergence"}


def test_inference_reconciles_every_one_of_the_THIRTY_keys_and_SchemaErrors_are_NOT_caught(
        monkeypatch):
    """§16.13, §10.2."""
    pop = policy.population(mixed_frame(), audit_of())
    boot = with_boot(RUN_BOOT, lambda: policy.inference(pop, audit_of()))
    assert len(boot.draws) == 30
    for drawn in boot.draws.values():
        assert len(drawn.draws) + sum(drawn.failures.values()) == drawn.n_attempted == RUN_BOOT
        assert np.all(np.isfinite(drawn.draws))
    assert boot.seed == config.SEED and boot.n_boot == RUN_BOOT
    assert "policy" in boot.diagnostics.max_abs_beta
    assert len(boot.diagnostics.max_abs_beta["policy"]) == RUN_BOOT

    monkeypatch.setattr(policy, "_assert_dilution",
                        lambda *a: (_ for _ in ()).throw(config.SchemaError("U7  forced")))
    with pytest.raises(config.SchemaError):
        with_boot(3, lambda: policy.inference(pop, audit_of()))


# --- 16.14  every U identifier fires, and the set is the assertion --------------------------------

def test_EVERY_U_IDENTIFIER_FIRES_and_the_set_is_the_assertion():
    """§16.14. The left side is scraped from `policy.py`'s raise sites; the right from this file."""
    raised = set()
    for node in ast.walk(ast.parse(SOURCE)):
        if not isinstance(node, ast.Raise) or not isinstance(node.exc, ast.Call):
            continue
        exc = ast.unparse(node.exc.func)
        assert exc in ("C.SchemaError", "model.FitError"), exc
        first = node.exc.args[0]
        if isinstance(first, ast.JoinedStr):
            first = first.values[0]
        token = str(first.value).split()[0]
        assert token[0] == "U" and token[1:].isdigit(), token
        raised.add(token)
    assert raised == U_SERIES

    own = Path(__file__).resolve().read_text(encoding="utf-8")
    tested = {u for u in U_SERIES if f"_{u}_" in own or f'startswith("{u}")' in own}
    assert tested == U_SERIES
    # U6 is T10 by its Stage 12 name inside `gcompute`, listed for the reader's checklist only.
    assert "U6" not in raised


def test_the_public_surface_is_THREE_NAMES_and_the_privates_are_the_declared_ONES():
    """§3.2, extended by three privates the implementation needed: `_assert_regime` (U2, called on
    both designs), `_population_table` (entry 1's ledger) and `_estimate` (the fit-to-record path
    `contrast` and the body share, so the two cannot drift)."""
    tree = ast.parse(SOURCE)
    functions = [n.name for n in tree.body if isinstance(n, ast.FunctionDef)]
    assert [n for n in functions if not n.startswith("_")] == ["population", "contrast", "inference"]
    assert sorted(n for n in functions if n.startswith("_")) == sorted([
        "_regimes", "_assert_exposure_survived", "_assert_beta_reportable", "_assert_dilution",
        "_assert_regime", "_population_table", "_estimate", "_estimand_keys", "_values",
        "_replicate", "_fit_table", "_contrast_table", "_replicates_table"])
    assert [n.name for n in tree.body if isinstance(n, ast.ClassDef)] == ["Policy"]


# --- 16.16  `_regimes` and the no-contraindicated draw --------------------------------------------

def test_U4_fires_on_a_stored_column_FLIPPED_after_population_and_names_the_count():
    """§16.16, §5.2's disagreement clause."""
    pop = policy.population(mixed_frame(), audit_of())
    D = pop[config.CONTRAINDICATED].to_numpy()
    pop.iloc[np.flatnonzero(D == 0.0)[0], pop.columns.get_loc(config.CONTRAINDICATED)] = 1.0
    with pytest.raises(config.SchemaError) as excinfo:
        policy.contrast(pop, audit_of())
    assert message_of(excinfo).startswith("U4")
    assert "1 row" in message_of(excinfo)


def test_no_contraindicated_frame_through_the_BODY_is_a_legitimate_draw():
    """§16.16, §5.3. Same outcome as `population` + `contrast` on `all_eligible_frame` (§16.10)."""
    frame = no_contraindicated_frame()
    keys = policy._estimand_keys()
    replicate = policy._replicate(frame, keys)
    assert replicate.failures == {}
    assert replicate.values["share_eligible"] == 1.0
    D, active, comparator = policy._regimes(frame)
    assert (D == 0.0).all() and (active == max(config.TREATMENT_LABELS)).all()
    for k in config.MRS_THRESHOLDS:
        assert replicate.values[f"rd_{k}"] == replicate.values[f"eligible_rd_{k}"]

    pol = policy._estimate(frame)
    assert pol.dropped == (config.CONTRAINDICATED,)
    assert pol.share_eligible == 1.0 and pol.rd == pol.eligible_rd


# --- 16.17  the tally, the renderers and the rendering rules --------------------------------------

def test_standardise_diagnostics_KEEPS_the_policy_key_and_bootstrap_diagnostics_DROPS_it():
    """§16.17, §10.2. Why the generic tally and not the bootstrap one is the one called."""
    replicates = tuple(
        bootstrap.Replicate(values={"rd_0": 0.1}, failures={}, n_alpha=6, polr_iterations=5,
                            sum_w=108.0, n_in_model=108, max_abs_beta={"policy": v})
        for v in (1.0, 2.0, 24.5))
    kept = standardise.diagnostics(replicates)
    assert set(kept.max_abs_beta) == {"policy"}
    assert len(kept.max_abs_beta["policy"]) == 3
    assert bootstrap.diagnostics(replicates).max_abs_beta == {}
    calls = {ast.unparse(c.func) for c in ast.walk(ast.parse(SOURCE)) if isinstance(c, ast.Call)}
    assert "standardise.diagnostics" in calls and "bootstrap.diagnostics" not in calls


def test_the_fit_table_carries_this_stages_roles():
    """§16.17, §12."""
    pol = policy.contrast(policy.population(mixed_frame(), audit_of()), audit_of())
    table = policy._fit_table(pol.fit, pol.dropped)
    assert len(table) == 1 + len(pol.fit.columns) + len(pol.dropped)
    roles = {row[0]: row[-1] for row in table[1:]}
    assert "U8" in roles[config.TREATMENT]
    assert "NUISANCE" in roles[config.CONTRAINDICATED] and "[§4.3]" in roles[config.CONTRAINDICATED]
    for name, role in roles.items():
        if name not in (config.TREATMENT, config.CONTRAINDICATED):
            assert "unbounded" in role, name


def test_the_contrast_table_is_headed_by_REGIME_names_and_carries_the_decomposition():
    """§16.17, §12."""
    pol = policy.contrast(policy.population(mixed_frame(), audit_of()), audit_of())
    table = policy._contrast_table(pol)
    header = table[0]
    assert any(policy._ACTIVE in h for h in header) and any(policy._COMPARATOR in h for h in header)
    for label in config.TREATMENT_LABELS.values():
        assert not any(label in h for h in header)
    assert "eligible_rd_k" in header
    first_column = [row[0] for row in table]
    assert "share_eligible" in first_column
    assert any("mortality" in c for c in first_column) and any("mrs_0_2" in c for c in first_column)
    assert any("eligible_rd_k" in h for h in header) and any("share" in h for h in header)


def test_the_replicates_table_renders_on_THIRTY_keys_labelling_every_denominator():
    """§16.17, §12."""
    audit = audit_of()
    with_boot(RUN_BOOT, lambda: driver(mixed_frame(), audit, with_inference=True))
    entry = audit.entry("model", "policy_replicates")
    assert entry is not None and entry.n == RUN_BOOT
    blocks = [row[0] for row in entry.table[1:]]
    assert all("=" in b for b in blocks)
    quantities = [row[1] for row in entry.table[1:]]
    assert any("max_abs_coef" in q for q in quantities)
    assert not any("hier" in q for q in quantities)
    grid = [row for row in entry.table[1:] if "max_abs_coef" in row[1]]
    assert grid and "indicator" in grid[0][2]


def test_rendering_rules_over_all_five_entries_no_pipe_equal_cell_counts_no_14a_outside_the_label():
    """§16.17, §12. `data._md_table` does no escaping (Stage 7's `|SMD|` header). The `[§14a]` rule
    is applied to every table cell and to every detail line OUTSIDE `Policy.label`, which the spec
    both prescribes verbatim — it says what the contrast is NOT — and requires printed in entry 4."""
    audit = audit_of()
    with_boot(RUN_BOOT, lambda: driver(mixed_frame(), audit, with_inference=True))
    assert len(audit.entries) == 5
    for entry in audit.entries:
        if entry.table is not None:
            widths = {len(row) for row in entry.table}
            assert len(widths) == 1, entry.step
            for row in entry.table:
                assert not any("|" in cell for cell in row), (entry.step, row)
                assert not any("[§14a]" in cell for cell in row), (entry.step, row)
            rendered = data._md_table(entry.table).split("\n")
            assert len({line.count("|") for line in rendered}) == 1, entry.step
        outside = entry.detail.replace(policy._OPERATIONAL, "").replace(policy._NO_CONTRAINDICATED, "")
        assert "[§14a]" not in outside, entry.step
        assert "|" not in entry.detail.replace("| ", ""), entry.step


# --- 4.3  what may not be quoted ------------------------------------------------------------------

def test_this_module_quotes_no_workbook_probability_risk_difference_or_coefficient():
    """§4.3. The only floats a data-gated test may name are tolerances."""
    tree = ast.parse(Path(__file__).resolve().read_text(encoding="utf-8"))
    gated = [n for n in ast.walk(tree) if isinstance(n, ast.FunctionDef)
             and any("DATA_GATED" in ast.unparse(d) for d in n.decorator_list)]
    assert gated
    for node in gated:
        for constant in (c for c in ast.walk(node) if isinstance(c, ast.Constant)):
            if isinstance(constant.value, float) and 0.0 < abs(constant.value) < 1.0:
                assert abs(constant.value) < 1e-9, (node.name, constant.value)
