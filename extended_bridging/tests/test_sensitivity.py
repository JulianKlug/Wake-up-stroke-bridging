"""Stage 11 acceptance tests — §15 of `specs/stage11_multiplicity_subgroups_evalue.md`.

Banner comments are the `### 15.x` headings of that document, lower-cased. Sections whose subject
lives in another module go in that module's test file and are tagged there: R9 and the three
newly-public `bootstrap` names in `test_bootstrap.py` (§15.2, §15.11), the registry in
`test_config.py` (§15.10), the two entry points in `test_propensity.py`, the role column in
`test_balance.py` (§15.3), and `outcome.estimation_population` in `test_outcome.py` (§15.7, §15.12).

Tests needing the private workbook are marked `DATA_GATED` and tagged `[data-gated]` in their
docstrings. Everything else runs on `fixtures_stage11.py`'s constructed frames.

**No number this file produces appears in any document under `specs/`** — no adjusted p, no E-value,
no subgroup odds ratio and no arm limit (§1). What it asserts are shapes, counts, identities and
guards.

`C.N_BOOT` is moved to `RUN_BOOT` for every test that drives a replicate loop, which is
`test_bootstrap.py`'s own device: at the prespecified 2000 this file would be a nine-minute suite,
and every property §15 asserts is a property of the loop rather than of its length.
"""
from __future__ import annotations

import ast
import inspect
import math
from pathlib import Path
from typing import Final

import numpy as np
import pandas as pd
import pytest

import balance
import bootstrap
import config
import data
import model
import outcome
import propensity
import sensitivity
from fixtures_stage11 import (SUBGROUP, family_fixtures, interval_fixtures,
                              near_separated_subgroup_frame, separable_subgroup_frame,
                              subgroup_frame, tail_count_50)

DATA_GATED = pytest.mark.skipif(
    not config.DATA_XLSX.exists(),
    reason="the private workbook is gitignored and absent from this checkout")

MODULE = Path(sensitivity.__file__).resolve()
SOURCE = MODULE.read_text(encoding="utf-8")
THIS_FILE = Path(__file__).resolve().read_text(encoding="utf-8")

RUN_BOOT: Final[int] = 60


def message_of(excinfo) -> str:
    return str(excinfo.value)


@pytest.fixture
def boot_n(monkeypatch):
    """`C.N_BOOT` at RUN_BOOT for the duration of one test."""
    monkeypatch.setattr(config, "N_BOOT", RUN_BOOT)
    return RUN_BOOT


def fitted(df: pd.DataFrame):
    """(df, ps, est, audit) — Stages 6 and 8 over one frame, against ONE Audit."""
    audit = data.Audit(data.WORKBOOK)
    ps = propensity.fit(df, audit)
    return df, ps, outcome.primary(df, ps, audit), audit


@pytest.fixture(scope="module")
def stage11_run():
    """(df, ps, est, sec, audit, boot) from ONE Stage 6-10 drive at RUN_BOOT replicates.

    Module-scoped because it is the expensive thing in this file. `monkeypatch` is function-scoped,
    so the constant is moved and restored by hand — `test_bootstrap.py`'s own pattern for its own
    reason.
    """
    df = subgroup_frame()
    original = config.N_BOOT
    config.N_BOOT = RUN_BOOT
    try:
        audit = data.Audit(data.WORKBOOK)
        ps = propensity.fit(df, audit)
        bal = balance.assess(df, ps, audit)
        est = outcome.primary(df, ps, audit)
        sec = outcome.secondary(df, ps, audit)
        yield df, ps, bal, est, sec, audit, bootstrap.run(df, ps, est, sec, audit)
    finally:
        config.N_BOOT = original


# --- 15.0  the frames and fixtures this stage is tested on --------------------------------------------

def test_the_subgroup_frame_carries_WHAT_STAGE_TEN_DID_NOT_and_nothing_else_moves():
    df = subgroup_frame()
    for name in (*config.NEGATIVE_CONTROLS, *config.SUBGROUPS):
        assert name in df.columns
    for name in config.SUBGROUPS:
        assert str(df[name].dtype) == "Int64"          # derive.py's dtype, `<NA>` available
        assert set(df[name].dropna().unique()) <= {0, 1}
    # the four risk factors VARY, so C.PROPENSITY_FULL is not C.PROPENSITY_PRIMARY with four names
    # in `dropped` — a frame on which the arm cannot be told from the primary passes anything.
    for name in config.NEGATIVE_CONTROLS:
        assert df[name].nunique(dropna=False) == 2


@pytest.mark.parametrize("kind,token", [("empty_cell", "O6"), ("constant_product", "G8"),
                                        ("one_level", "G8")])
def test_each_degenerate_frame_reaches_the_route_it_was_built_for(kind, token):
    """§15.0's own claim, asserted rather than described. §8.3's two routes are DIFFERENT mechanisms
    and a fixture reaching the wrong one would test the wrong guard."""
    df, ps, est, _ = fitted(separable_subgroup_frame(kind))
    with pytest.raises(model.FitError) as e:
        sensitivity._subgroup_fit(df, ps, est.in_estimate, SUBGROUP)
    assert message_of(e).split()[0] == token


def test_the_near_separated_frame_makes_a_fit_that_SUCCEEDS_and_is_rejected():
    """*This is the fixture the whole of §8.4 rests on.* A frame on which the fit FAILED instead
    would pass a wrong implementation that has no G9 at all, because O1-O6 would already have
    raised — so the test asserts the fit reaches `polr`, converges, and is rejected after."""
    df, ps, est, _ = fitted(near_separated_subgroup_frame())
    mask = est.in_estimate & df[SUBGROUP].notna()
    sub = df.loc[mask].copy()
    sub[f"{config.TREATMENT}_x_{SUBGROUP}"] = (
        sub[config.TREATMENT].astype(float) * sub[SUBGROUP].astype(float))
    X, dropped = model.design(sub, (config.TREATMENT, SUBGROUP,
                                    f"{config.TREATMENT}_x_{SUBGROUP}"))
    assert dropped == () and tuple(X.columns) == (config.TREATMENT, SUBGROUP,
                                                  f"{config.TREATMENT}_x_{SUBGROUP}")
    fit = model.polr(X, sub[config.PRIMARY_OUTCOME].to_numpy(dtype=float),
                     ps.w.loc[mask].to_numpy(dtype=float))
    assert np.all(np.isfinite(fit.beta)) and fit.iterations > 0       # it CONVERGED

    with pytest.raises(model.FitError) as e:
        sensitivity._subgroup_fit(df, ps, est.in_estimate, SUBGROUP)
    assert message_of(e).startswith("G9")

    # ...and it crossed on a REPORTED quantity, never on `delta`, which is §8.4's decision.
    coefficients = dict(zip(X.columns, (float(v) for v in fit.beta)))
    beta = coefficients[config.TREATMENT]
    gamma = coefficients[f"{config.TREATMENT}_x_{SUBGROUP}"]
    assert max(abs(beta), abs(gamma), abs(beta + gamma)) >= config.POLR_MAX_ABS_BETA


def test_the_family_and_interval_fixtures_carry_EVERY_branch_section_six_and_seven_have():
    assert set(family_fixtures()) == {"complete", "thinned", "extra_keys", "none_p", "floor",
                                      "no_beta"}
    assert set(interval_fixtures()) == {"spanning", "lo_at_null", "hi_at_null", "both_at_null",
                                        "wide", "asymmetric", "outside"}
    assert len(tail_count_50()) == config.N_BOOT


# --- 15.1  the preconditions, one frame per branch ----------------------------------------------------
#
# Ten H identifiers and R9, and the count comes from §9.1 rather than from a list here. Each is
# asserted on a frame that reaches it and nothing else, and each `C.SchemaError` message is matched by
# its IDENTIFIER — `match="H4"` and not on prose — so a reworded message does not break the suite and
# a renumbered guard does.

def test_H1_fires_when_the_ARMS_OWN_propensity_is_passed_as_the_reference_half(boot_n):
    df = subgroup_frame()
    audit = data.Audit(data.WORKBOOK)
    arm_ps = propensity.fit_full(df, audit)
    with pytest.raises(config.SchemaError, match="H1"):
        sensitivity.full_covariate(df, arm_ps, audit)


def test_H2_fires_on_a_MISALIGNED_propensity(boot_n):
    df, ps, _, audit = fitted(subgroup_frame())
    bent = ps.__class__(e=ps.e[::-1], w=ps.w[::-1], in_model=ps.in_model, ess=ps.ess, fit=ps.fit,
                        dropped=ps.dropped, spec=ps.spec)
    with pytest.raises(config.SchemaError, match="H2"):
        sensitivity.full_covariate(df, bent, audit)


def test_H3_fires_on_an_Interval_whose_p_is_None():
    sec, boot = family_fixtures()["none_p"]
    with pytest.raises(config.SchemaError, match="H3"):
        sensitivity.multiplicity(sec, boot, data.Audit(data.WORKBOOK))


def test_H3_fires_on_a_family_member_with_no_DRAWS_key_at_all():
    sec, boot = family_fixtures()["complete"]
    starved = f"{config.BINARY_OUTCOMES[0]}.rd"
    thinned = bootstrap.Bootstrap(
        boot.seed, boot.n_boot, {k: v for k, v in boot.draws.items() if k != starved},
        {k: v for k, v in boot.intervals.items() if k != starved}, boot.diagnostics)
    with pytest.raises(config.SchemaError, match="H3"):
        sensitivity.multiplicity(sec, thinned, data.Audit(data.WORKBOOK))


def test_H5_fires_when_there_is_no_beta_interval_at_all():
    """The primary is UNCORRECTED, not UNREPORTED. Without this check an absent key is a bare
    `KeyError` from a reported field's initialiser (§6.1)."""
    sec, boot = family_fixtures()["no_beta"]
    with pytest.raises(config.SchemaError, match="H5"):
        sensitivity.multiplicity(sec, boot, data.Audit(data.WORKBOOK))


def test_H6_fires_on_a_NON_FINITE_odds_ratio():
    df, ps, est, audit = fitted(subgroup_frame())
    bal = balance.assess(df, ps, audit)
    boot = _bootstrap_over(interval_fixtures()["spanning"], est)
    broken = _replace_primary(est, odds_ratio=float("nan"))
    with pytest.raises(config.SchemaError, match="H6"):
        sensitivity.e_value_primary(broken, bal, boot, audit)


def test_H7_fires_when_the_point_estimate_is_OUTSIDE_ITS_OWN_INTERVAL():
    """Stage 10 §8.1: a percentile interval is not a function of the point estimate, so this is
    arithmetically reachable. The pilot reports it as "the interval spans the null" — a different and
    benign condition — and that is the sharpest single reason this guard raises (§7.3)."""
    df, ps, est, audit = fitted(subgroup_frame())
    bal = balance.assess(df, ps, audit)
    outside = interval_fixtures()["outside"]
    boot = _bootstrap_over(outside, est)
    broken = _replace_primary(est, beta=outside.lo - 1.0)
    with pytest.raises(config.SchemaError, match="H7"):
        sensitivity.e_value_primary(broken, bal, boot, audit)


@pytest.mark.parametrize("kind", ["empty_cell", "constant_product", "one_level"])
def test_H8_fires_when_a_POINT_ESTIMATE_subgroup_fit_raises_and_NAMES_THE_TOKEN(kind, boot_n):
    """§8.5, and the two doors into `_subgroup_fit` behave DIFFERENTLY on one frame: inside the
    replicate loop this same token is dropped and counted, and here it raises. Omitting the subgroup
    from `estimates` is `pilots/analysis.py:811-813`'s uncounted `continue`."""
    df, ps, est, audit = fitted(separable_subgroup_frame(kind))
    with pytest.raises(config.SchemaError) as e:
        sensitivity.subgroups(df, ps, est, audit)
    assert "H8" in message_of(e)
    assert SUBGROUP in message_of(e)
    assert ("O6" if kind == "empty_cell" else "G8") in message_of(e)


def test_H8_fires_on_the_NEAR_SEPARATED_frame_too_and_the_fit_had_CONVERGED(boot_n):
    df, ps, est, audit = fitted(near_separated_subgroup_frame())
    with pytest.raises(config.SchemaError) as e:
        sensitivity.subgroups(df, ps, est, audit)
    assert "H8" in message_of(e) and "G9" in message_of(e)


def test_H9_fires_when_in_estimate_is_not_a_SUBSET_of_in_model(boot_n):
    df, ps, est, audit = fitted(subgroup_frame())
    widened = est.in_estimate | ~ps.in_model
    with pytest.raises(config.SchemaError, match="H9"):
        sensitivity.subgroups(df, ps, _replace_primary(est, in_estimate=widened), audit)


def test_H9_fires_on_a_MISALIGNED_index(boot_n):
    df, ps, est, audit = fitted(subgroup_frame())
    moved = est.in_estimate.copy()
    moved.index = moved.index[::-1]
    with pytest.raises(config.SchemaError, match="H9"):
        sensitivity.subgroups(df, ps, _replace_primary(est, in_estimate=moved), audit)


def test_H10_fires_on_EITHER_HALF_of_the_pair_being_wrong(boot_n):
    """`Arm.__post_init__` is the check from the other side, on the RECORD rather than on the call,
    and it is the one a caller reaches by building the record by hand (§3.1, §4.5)."""
    df, ps, _, audit = fitted(subgroup_frame())
    arm = sensitivity.full_covariate(df, ps, audit)
    fields = dict(ps=arm.ps, balance=arm.balance, estimate=arm.estimate, reference=arm.reference,
                  seed=arm.seed, n_boot=arm.n_boot, draws=arm.draws, intervals=arm.intervals,
                  diagnostics=arm.diagnostics,
                  differs_only_in_specification=arm.differs_only_in_specification)
    with pytest.raises(config.SchemaError, match="H10"):
        sensitivity.Arm(**{**fields, "ps": arm.reference})
    with pytest.raises(config.SchemaError, match="H10"):
        sensitivity.Arm(**{**fields, "reference": arm.ps})


def test_every_H_identifier_is_matched_by_its_IDENTIFIER_and_never_by_prose():
    """§15.1's rule as a scan over this file: `match="H4"` and not on a sentence, so a reworded
    message does not break the suite and a renumbered guard does."""
    tree = ast.parse(THIS_FILE)
    matched = {
        keyword.value.value
        for node in ast.walk(tree)
        if isinstance(node, ast.Call) and isinstance(node.func, ast.Attribute)
        and node.func.attr == "raises"
        for keyword in node.keywords
        if keyword.arg == "match" and isinstance(keyword.value, ast.Constant)}
    assert {f"H{i}" for i in (1, 2, 3, 5, 6, 7, 9, 10)} <= matched
    assert all(len(pattern) <= 3 for pattern in matched)


def _replace_primary(est: outcome.Primary, **changes) -> outcome.Primary:
    import dataclasses
    return dataclasses.replace(est, **changes)


def _bootstrap_over(interval: bootstrap.Interval, est: outcome.Primary) -> bootstrap.Bootstrap:
    """A `Bootstrap` carrying ONE hand-built `beta` interval and draws that reconcile with it."""
    draws = bootstrap.Draws(quantity="beta",
                            draws=np.linspace(interval.lo, interval.hi, interval.n_draws),
                            n_attempted=interval.n_draws, failures={})
    return bootstrap.Bootstrap(config.SEED, config.N_BOOT, {"beta": draws}, {"beta": interval},
                               bootstrap.diagnostics(()))


# --- 15.3  the arm reuses every landed stage and adds no estimator ------------------------------------

def attribute_calls(source: str) -> set[str]:
    """`module.name` for every attribute access in `source`. `test_balance.py:1425`'s instrument."""
    return {f"{node.value.id}.{node.attr}"
            for node in ast.walk(ast.parse(source))
            if isinstance(node, ast.Attribute) and isinstance(node.value, ast.Name)}


def function_named(name: str) -> ast.FunctionDef:
    return next(n for n in ast.parse(SOURCE).body
                if isinstance(n, ast.FunctionDef) and n.name == name)


def test_each_propensity_entry_point_is_named_in_EXACTLY_THE_RIGHT_FUNCTIONS():
    """The one scan that catches the two bodies being wired to each other's specification (§15.3).

    A subgroup body calling `fit_full` estimates the [§13] subgroups under a sensitivity
    specification no section asks for; an arm body calling `fit` produces intervals whose replicates
    are the primary's while the point estimates are the arm's, which is R9's own failure one module
    over — and neither raises anything.

    `fit_full` is named TWICE and that is correct: once in the arm's replicate body, and once in
    `full_covariate` for the point estimate the replicates are the distribution of. `fit` is named
    ONCE, in the subgroup body — `subgroups` itself never refits, because it takes the [§7]
    `Propensity` as an argument.
    """
    named = {n.name: attribute_calls(ast.unparse(n)) for n in ast.parse(SOURCE).body
             if isinstance(n, ast.FunctionDef)}
    assert {f for f, calls in named.items() if "propensity.fit" in calls} == {"_subgroup_replicate"}
    assert {f for f, calls in named.items() if "propensity.fit_full" in calls} == {
        "_arm_replicate", "full_covariate"}
    assert "propensity.fit_full" not in named["_subgroup_replicate"]
    assert "propensity.fit" not in named["_arm_replicate"]


def test_that_the_attribute_scan_FIRES():
    assert "propensity.fit_full" in attribute_calls("ps = propensity.fit_full(df, audit)\n")


def test_the_arm_carries_the_full_specification_and_the_reference_carries_the_primary(boot_n):
    df, ps, _, audit = fitted(subgroup_frame())
    arm = sensitivity.full_covariate(df, ps, audit)
    assert arm.ps.spec is config.PROPENSITY_FULL
    assert arm.reference.spec is config.PROPENSITY_PRIMARY
    assert arm.reference is ps                      # READ, never refitted here
    assert arm.balance.spec is config.PROPENSITY_FULL
    assert arm.seed == config.SEED and arm.n_boot == config.N_BOOT


def test_NO_p_value_on_ANY_of_the_arms_seven_keys(boot_n):
    """[§13] prescribes sensitivity analyses that REPORT and prescribes no test of any of them. A p
    on the arm's `beta` would be an eighth test [§13]'s correction is not told about (§5.5)."""
    df, ps, _, audit = fitted(subgroup_frame())
    arm = sensitivity.full_covariate(df, ps, audit)
    assert len(arm.draws) == 7
    assert set(arm.draws) == {"beta", *(f"rd_{k}" for k in config.MRS_THRESHOLDS)}
    assert all(interval.p is None for interval in arm.intervals.values())


def test_the_lambda_FALSE_is_LOAD_BEARING_and_not_decorative(boot_n):
    """§15.3's companion: the SAME helper with `_tested_subgroup` produces p-values on the same
    draws, so the argument is what enforces the rule rather than the draws happening to lack one.

    `bootstrap._tested` is deliberately not reused for either caller, and this is why: `_tested`
    returns True for `"beta"`, so an implementation reaching for it emits a p on the arm's `beta`.
    """
    df, ps, _, audit = fitted(subgroup_frame())
    arm = sensitivity.full_covariate(df, ps, audit)
    with_a_rule = sensitivity._intervals(arm.draws, bootstrap._tested)
    assert with_a_rule["beta"].p is not None
    assert arm.intervals["beta"].p is None
    assert with_a_rule["beta"].lo == arm.intervals["beta"].lo     # the same draws, the same limits


def test_the_arms_diagnostics_reconcile_and_the_augmented_fields_are_EMPTY(boot_n):
    """`max_abs_beta` and `or_corrected` are empty because the arm produces no binary estimate and no
    `m_a(X)` — `Diagnostics`' own "AUGMENTED ONLY" rule satisfied by there being nothing augmented
    (§3.1). `n_alpha` sums to the replicates that reached `outcome.primary`."""
    df, ps, _, audit = fitted(subgroup_frame())
    arm = sensitivity.full_covariate(df, ps, audit)
    assert arm.diagnostics.max_abs_beta == {}
    assert set(arm.diagnostics.or_corrected.values()) == {0}
    reached = sum(arm.diagnostics.n_alpha.values())
    assert reached == len(arm.draws["beta"].draws)
    assert sum(arm.diagnostics.polr_iterations.values()) == reached
    assert arm.diagnostics.sum_w.size == sum(arm.diagnostics.n_in_model.values())


def test_the_arms_beta_draws_DIFFER_from_the_primarys_on_the_SAME_drawn_frames(boot_n):
    """The behavioural falsifier, with the companion that makes it attributable: the k-th drawn frame
    is IDENTICAL in the two runs, so the difference is the specification and nothing else."""
    df, ps, est, audit = fitted(subgroup_frame())
    sec = outcome.secondary(df, ps, audit)
    boot = bootstrap.run(df, ps, est, sec, audit)
    arm = sensitivity.full_covariate(df, ps, audit)
    assert len(arm.draws["beta"].draws) == len(boot.draws["beta"].draws)
    assert not np.allclose(arm.draws["beta"].draws, boot.draws["beta"].draws)


# --- 15.4  benjamini-hochberg, and the vector that catches both defects -------------------------------

def test_the_vector_that_catches_BOTH_defects_to_exact_float_equality():
    """§6.3. Registry-order (0.04, 0.01, 0.03) at m = 3 gives (0.04, 0.03, 0.04).

    Three distinct wrong answers from one input: an implementation missing the running min returns
    0.045 in the second position, one missing the unsort returns (0.03, 0.04, 0.04), and one missing
    both returns (0.03, 0.045, 0.04).
    """
    got = sensitivity.benjamini_hochberg(np.array([0.04, 0.01, 0.03]))
    # to EXACT float equality on the products, and never to a tolerance
    assert list(got) == [0.04 * 3 / 3, 0.01 * 3 / 1, 0.04 * 3 / 3]
    assert got[1] != 0.03 * 3 / 2                                 # the running min fired
    assert list(got) != [0.03, 0.04, 0.04]                        # the unsort fired


def test_the_COMPANION_case_a_wrong_implementation_ALSO_gets_right():
    """§6.3 labels it as such, so that a suite consisting only of it cannot look like coverage."""
    assert list(sensitivity.benjamini_hochberg(np.array([0.01, 0.02, 0.03]))) == [0.03, 0.03, 0.03]


def test_a_family_of_ONE_is_bit_identical_to_its_input():
    p = np.array([0.037])
    assert list(sensitivity.benjamini_hochberg(p)) == list(p)


def test_the_cap_at_one_is_specified_and_INERT():
    """q_(m) = min(1, m*p_(m)/m) = p_(m) for any valid p, and every other q_(i) is at most that, so
    the cap can bind only on an input §6.1 already rejects (§6.2)."""
    def uncapped(p: np.ndarray) -> np.ndarray:
        """The same step-up with the `np.minimum(1.0, ...)` REMOVED, and nothing else changed."""
        m = p.size
        order = np.argsort(p, kind="stable")
        q = p[order] * m / np.arange(1, m + 1)
        q = np.minimum.accumulate(q[::-1])[::-1]
        out = np.empty_like(q)
        out[order] = q
        return out

    g = np.random.default_rng(config.SEED)
    grid = np.concatenate([[1.0 / (config.N_BOOT - 1)],
                           2.0 * np.arange(1, config.N_BOOT // 2) / (config.N_BOOT - 2)])
    for m in (3, 4):
        for _ in range(2000):
            p = g.choice(grid, m)
            capped = sensitivity.benjamini_hochberg(p)
            assert np.all(capped <= 1.0)
            # the cap BINDS in zero of them: with it removed the answer is bit-identical
            assert list(capped) == list(uncapped(p))


def test_permutation_invariance_and_validity_over_the_ACHIEVABLE_grid():
    """The grid is `{1/(B+1)} u {2k/B}` and NOT `U(0, 1)`, because TIES are the reachable condition
    and a uniform draw never produces one (§6.4)."""
    g = np.random.default_rng(config.SEED)
    grid = np.concatenate([[1.0 / (config.N_BOOT - 1)],
                           2.0 * np.arange(1, config.N_BOOT // 2) / (config.N_BOOT - 2)])
    ties = 0
    for _ in range(2000):
        m = int(g.integers(3, 5))
        p = g.choice(grid, m)
        ties += int(len(set(p)) < m)
        adjusted = sensitivity.benjamini_hochberg(p)
        assert np.all(adjusted >= p - 1e-15) and np.all(adjusted <= 1.0)
        order = g.permutation(m)
        permuted = sensitivity.benjamini_hochberg(p[order])
        assert np.allclose(permuted, adjusted[order], rtol=1e-12, atol=0.0)
    assert ties > 0, "the grid produced no tie, so the stable sort was never exercised"


def test_the_oracle_is_statsmodels_and_it_agrees_over_the_achievable_grid():
    """§6.4: `multipletests(..., method="fdr_bh")` is the oracle and MAY NOT be imported into a
    shipped module — an adjusted p-value is a REPORTED ESTIMATE, which is what `pyproject.toml`'s
    test-only policy exists to keep statsmodels out of."""
    multipletests = pytest.importorskip("statsmodels.stats.multitest").multipletests
    g = np.random.default_rng(config.SEED)
    grid = np.concatenate([[1.0 / (config.N_BOOT - 1)],
                           2.0 * np.arange(1, config.N_BOOT // 2) / (config.N_BOOT - 2)])
    hand = [np.array(v) for v in ([0.04, 0.01, 0.03], [0.01, 0.02, 0.03], [0.037],
                                  [0.05, 0.05, 0.05], [1.0, 1.0, 1.0, 1.0],
                                  [0.001, 0.5, 0.5, 0.9])]
    for p in [*hand, *(g.choice(grid, int(g.integers(3, 5))) for _ in range(3000))]:
        assert np.allclose(sensitivity.benjamini_hochberg(p),
                           multipletests(p, method="fdr_bh")[1], rtol=1e-12, atol=0.0)


def test_the_family_sizes_come_from_by_family_and_NOT_from_counting_intervals():
    """Stage 10 emits 26 keys of which 8 carry a p, so `len(boot.intervals)` is neither 3 nor 4
    (§6.1). Asserted on a `Bootstrap` carrying EXTRA keys, so an implementation counting intervals
    gets a different m and fails."""
    sec, complete = family_fixtures()["complete"]
    _, extra = family_fixtures()["extra_keys"]
    assert len(extra.intervals) > len(complete.intervals)
    plain = sensitivity.multiplicity(sec, complete, data.Audit(data.WORKBOOK))
    widened = sensitivity.multiplicity(sec, extra, data.Audit(data.WORKBOOK))
    assert {f: (c.m_declared, c.m_used) for f, c in plain.families.items()} == {
        "secondary": (3, 3), "safety": (4, 4)}
    assert ({f: c.adjusted for f, c in plain.families.items()}
            == {f: c.adjusted for f, c in widened.families.items()})


def test_H4s_denominator_is_THAT_KEYS_n_draws_and_a_check_against_N_BOOT_is_WEAKER():
    """§15.4, and it is a test of the CHECK and not of the procedure.

    Two `Interval`s on one family with different `n_draws` — 1982 and 1998 — and a raw p at exactly
    `1/1983`. On the 1982-draw key it is the floor and passes; the same value on the 1998-draw key is
    BELOW `1/1999` and raises H4. A check written against `C.N_BOOT` floors at `1/2001` and passes
    BOTH — so it accepts a p below the smallest value `bootstrap_p` could have returned for that key,
    which is the one thing this precondition exists to catch.
    """
    sec, boot = family_fixtures()["floor"]
    keys = [f"{k}.rd" for k in config.BINARY_OUTCOMES]
    passing, failing = boot.intervals[keys[0]], boot.intervals[keys[1]]
    assert passing.n_draws > failing.n_draws              # so passing's floor is the SMALLER one
    assert passing.p == 1.0 / (passing.n_draws + 1)       # exactly ON its own floor
    assert failing.p < 1.0 / (failing.n_draws + 1)        # BELOW its own
    assert failing.p > 1.0 / (config.N_BOOT + 1)          # the WEAKER check would accept it

    with pytest.raises(config.SchemaError, match="H4") as e:
        sensitivity.multiplicity(sec, boot, data.Audit(data.WORKBOOK))
    assert keys[1] in message_of(e) and keys[0] not in message_of(e)


# --- 15.5  the absent p, and the family denominator `[unreachable on v7]` -----------------------------

def test_a_member_with_no_INTERVAL_is_disclosed_as_a_row_and_is_not_an_error():
    """§6.5. `m` is the number of tests PERFORMED, so the thinned family's scaling factor is
    `m_used`; the alternative — holding `m` at the declared size — is also valid and is conservative,
    and it is rejected rather than overlooked, because it inflates every adjusted p in the family to
    pay for a test nobody ran."""
    sec, boot = family_fixtures()["thinned"]
    audit = data.Audit(data.WORKBOOK)
    result = sensitivity.multiplicity(sec, boot, audit)
    starved = config.BINARY_OUTCOMES[0]
    family = result.families[config.OUTCOMES[starved].family]

    assert family.m_used == family.m_declared - 1
    assert starved not in family.raw and starved not in family.adjusted
    assert [name for name, _ in family.absent] == [starved]
    reason = dict(family.absent)[starved]
    assert reason.startswith("no interval:")
    assert str(config.ci_min_draws()) in reason
    assert str(len(boot.draws[f"{starved}.rd"].draws)) in reason

    # ...and the scaling factor really is m_used, checkable by hand from the two printed counts.
    hand = sensitivity.benjamini_hochberg(np.array(list(family.raw.values())))
    assert list(family.adjusted.values()) == list(hand)

    # ...and it is a ROW of the rendered table, never a blank and never absent.
    entry = audit.entry("model", "multiplicity")
    rows = {row[0]: row[1:] for row in entry.table[1:]}
    assert starved in rows
    assert reason in rows[starved]


def test_m_used_zero_is_empty_dicts_plus_a_FULL_absent_and_BH_is_not_called_at_all():
    """§6.5's "says so", and there is no sentinel and no sentence field: a reader given
    `m_declared = 4`, `m_used = 0` and four `absent` rows has been told everything."""
    sec, boot = family_fixtures()["complete"]
    emptied = bootstrap.Bootstrap(
        boot.seed, boot.n_boot, boot.draws,
        {k: v for k, v in boot.intervals.items()
         if k == "beta" or config.OUTCOMES[k.split(".")[0]].family != "safety"},
        boot.diagnostics)
    result = sensitivity.multiplicity(sec, emptied, data.Audit(data.WORKBOOK))
    safety = result.families["safety"]
    assert safety.m_declared == 4 and safety.m_used == 0
    assert safety.raw == {} and safety.adjusted == {}
    assert len(safety.absent) == 4
    assert result.families["secondary"].m_used == 3          # the other family is untouched


def test_the_primary_is_a_FIELD_and_the_safety_family_carries_the_descriptive_only_LABEL():
    sec, boot = family_fixtures()["complete"]
    result = sensitivity.multiplicity(sec, boot, data.Audit(data.WORKBOOK))
    assert result.primary_uncorrected == boot.intervals["beta"].p
    assert result.families["safety"].descriptive_only is True
    assert result.families["secondary"].descriptive_only is False
    # no verdict of any kind lives on either record
    assert not hasattr(result, "significant")
    assert not any(f.startswith("significant") for f in vars(result.families["safety"]))


# --- 15.6  the e-value ---------------------------------------------------------------------------------

def test_the_two_closed_form_values_bit_exact():
    assert sensitivity.e_value(1.0) == 1.0
    assert sensitivity.e_value(2.0) == 2.0 + math.sqrt(2.0)


def test_E_is_symmetric_under_inversion_and_the_pin_is_ROOT_FIRST_THEN_INVERT():
    """`sqrt(1/OR)` and `1/sqrt(OR)` are mathematically equal and NOT bit-identical for every input,
    so the order is pinned. Asserted to rtol=1e-15 because exactness is NOT claimed (§7.1)."""
    for odds_ratio in (0.25, 0.3, 0.7, 0.9, 0.99):
        assert sensitivity.e_value(math.sqrt(odds_ratio)) == pytest.approx(
            sensitivity.e_value(math.sqrt(1.0 / odds_ratio)), rel=1e-15, abs=0.0)


def test_E_is_strictly_monotone_away_from_the_null():
    grid = np.geomspace(1.0 + 1e-9, 1e6, 400)
    values = [sensitivity.e_value(rr) for rr in grid]
    assert all(b > a for a, b in zip(values, values[1:]))


def test_E_is_NOT_SMOOTH_at_the_null_and_that_is_pinned_as_a_PROPERTY():
    """A change of 1e-12 in the odds ratio moves E by about 7e-07 — six orders of magnitude of
    amplification. Consequence for reporting: two decimals, and the statement that an E-value of 1.05
    and one of 1.15 differ by about a hundredfold in the odds ratio behind them (§7.4)."""
    assert sensitivity.e_value(math.sqrt(1.0 + 1e-12)) - 1.0 > 1e-7


def test_the_E_value_is_BOUNDED_because_G7_bounds_the_coefficient():
    """G7 is what makes every numerical edge case in this deliverable unreachable, and that is worth
    asserting because a reader will otherwise assume the E-value is unbounded (§7.4)."""
    biggest = math.exp(config.POLR_MAX_ABS_BETA)
    assert math.isfinite(sensitivity.e_value(math.sqrt(biggest)))
    assert sensitivity.e_value(math.sqrt(biggest)) < 1e4


@pytest.mark.parametrize("kind", ["spanning", "lo_at_null", "hi_at_null", "both_at_null", "wide"])
def test_the_limits_E_value_is_EXACTLY_ONE_whenever_the_interval_spans_the_null(kind):
    """Not approximately, not nan, and never the formula applied to a null-crossing limit — which is
    minimised at 1 and increases in BOTH directions once inverted, so a limit on the far side of the
    null scores as strong evidence (§7.3). Ties AT the null count as spanning."""
    df, ps, est, audit = fitted(subgroup_frame())
    bal = balance.assess(df, ps, audit)
    interval = interval_fixtures()[kind]
    est = _replace_primary(est, beta=0.0, odds_ratio=1.0)
    ev = sensitivity.e_value_primary(est, bal, _bootstrap_over(interval, est), audit)
    assert ev.spans_null is True
    assert ev.e_limit == 1.0
    assert ev.limit is None
    assert ev.n_draws == interval.n_draws


def test_the_E_value_taken_from_a_NULL_CROSSING_LIMIT_would_be_LARGE(boot_n):
    """Definition of done item 5: the 1.0 rule's value is a MEASURED fact in the suite rather than a
    sentence. The `wide` fixture reaches to -4.0, and the formula applied to it returns a number that
    reads as robustness while meaning the opposite."""
    wide = interval_fixtures()["wide"]
    wrong = sensitivity.e_value(math.sqrt(math.exp(wide.lo)))
    assert wrong > 10.0
    assert sensitivity.e_value(1.0) == 1.0


def test_the_nearest_limit_is_selected_BY_SIGN_and_never_by_the_smaller_absolute_value():
    """A `min(|lo|, |hi|)` implementation picks the wrong limit ONLY in the spanning case, where the
    1.0 rule then masks it — so the defect is invisible until the first analysis whose interval
    excludes the null, which is the analysis nobody wants to find it in (§7.3)."""
    df, ps, est, audit = fitted(subgroup_frame())
    bal = balance.assess(df, ps, audit)
    interval = interval_fixtures()["asymmetric"]
    est = _replace_primary(est, beta=1.0, odds_ratio=math.e)
    ev = sensitivity.e_value_primary(est, bal, _bootstrap_over(interval, est), audit)
    assert ev.spans_null is False
    assert ev.limit == interval.lo
    # taking `hi` instead gives a VISIBLY different number, which is what makes the fixture work
    assert sensitivity.e_value(math.sqrt(math.exp(interval.hi))) > 2.0 * ev.e_limit


def test_all_three_routes_to_null_spanning_AGREE_on_the_tail_count_fifty_fixture():
    """§7.2's three routes, on the ONE fixture where route 3 can be told apart from routes 1 and 2:
    Stage 10 §9.4 measured agreement holding there under `inverted_cdf` and FAILING under `linear`,
    so it is the frame on which a wrong percentile method is visible."""
    draws = tail_count_50()
    lo, hi = bootstrap.percentile_ci(draws, config.CI_LEVEL)
    p = bootstrap.bootstrap_p(draws)

    # The two-sided alpha is SNAPPED exactly as `percentile_ci` snaps its quantile, and for the same
    # reason: `1.0 - 0.95` is 0.050000000000000044, `bootstrap_p` returns a multiple of 2/B, and at
    # tail 50 the p is 0.05 EXACTLY — so an unsnapped `p >= 1.0 - C.CI_LEVEL` reads False and route 3
    # contradicts routes 1 and 2 at precisely the boundary route 3 exists to decide (§7.2).
    alpha = round(1.0 - config.CI_LEVEL, 9)
    assert alpha == 0.05 and (1.0 - config.CI_LEVEL) != 0.05

    route_one = lo <= 0.0 <= hi
    route_two = math.exp(lo) <= 1.0 <= math.exp(hi)
    route_three = p >= alpha
    assert route_one is True                       # Stage 10 §9.4 measured this under inverted_cdf
    assert route_one == route_two == route_three

    # ...and it is the fixture where a WRONG percentile method is visible: under `linear` the same
    # draws give an interval EXCLUDING zero against the same p, so routes 2 and 3 part company.
    linear_lo, linear_hi = np.percentile(draws, [2.5, 97.5], method="linear")
    assert not (linear_lo <= 0.0 <= linear_hi)


def test_route_three_uses_GREATER_OR_EQUAL_and_the_boundary_is_the_case_it_decides():
    """`bootstrap_p` counts a draw of exactly 0.0 in both tails, so a p landing exactly on
    `1 - C.CI_LEVEL` is the arithmetic image of a limit landing exactly on 0 — which the tie rule
    calls spanning. With `>` the two rules disagree at exactly that point (§7.2)."""
    boundary = 1.0 - config.CI_LEVEL
    assert (boundary >= 1.0 - config.CI_LEVEL) is True
    assert (boundary > 1.0 - config.CI_LEVEL) is False
    df, ps, est, audit = fitted(subgroup_frame())
    bal = balance.assess(df, ps, audit)
    for kind in ("lo_at_null", "hi_at_null"):
        interval = interval_fixtures()[kind]
        assert interval.p == boundary
        est_at_null = _replace_primary(est, beta=0.0, odds_ratio=1.0)
        ev = sensitivity.e_value_primary(est_at_null, bal,
                                         _bootstrap_over(interval, est_at_null), audit)
        assert ev.spans_null is True


@pytest.mark.parametrize("n", [41, 97, 400, 1998, 2000])
def test_exp_of_a_percentile_limit_IS_the_percentile_limit_of_exp_under_the_pinned_method(n):
    """§7.2 route two, and it is a SECOND consequence of Stage 10 §8.2's pin that §8.2 does not draw:
    it is what licenses reporting an odds-ratio interval as `exp` of a `beta` interval at all.
    `inverted_cdf` selects an ORDER STATISTIC and a strictly monotone transform commutes with
    order-statistic selection exactly; `linear` interpolates and does not."""
    # DRAWN and not evenly spaced: on a uniform grid `linear` interpolates onto a grid point at
    # some sizes and the second half of this test passes for the wrong reason. Measured.
    draws = np.sort(np.random.default_rng(config.SEED).normal(0.5, 1.0, n))
    lo, hi = bootstrap.percentile_ci(draws, config.CI_LEVEL)
    q_lo = round(50.0 * (1.0 - config.CI_LEVEL), 9)
    q_hi = round(100.0 - q_lo, 9)
    # `np.exp` and NOT `math.exp`: the two differ by an ULP on some inputs, and the property under
    # test is that the ORDER STATISTIC is the same element — not that two libm implementations agree.
    pinned = np.percentile(np.exp(draws), [q_lo, q_hi], method=config.PERCENTILE_METHOD)
    assert list(np.exp(np.array([lo, hi]))) == list(pinned)
    # The companion — `linear` does NOT commute — holds only where `linear` actually interpolates.
    # At n = 41 both quantile positions land on integers (0.025*40 = 1.0 and 0.975*40 = 39.0), so it
    # selects the same order statistics and the two methods coincide for a structural reason rather
    # than a numerical one. Measured, and asserted as the condition rather than assumed away.
    positions = [q / 100.0 * (n - 1) for q in (q_lo, q_hi)]
    interpolated = np.percentile(np.exp(draws), [q_lo, q_hi], method="linear")
    if any(position != int(position) for position in positions):
        assert list(np.exp(np.array([lo, hi]))) != list(interpolated)
    else:
        assert list(np.exp(np.array([lo, hi]))) == list(interpolated)


def test_n_draws_is_populated_when_there_is_NO_INTERVAL_AT_ALL():
    """§7.4: the `None` must say "too few surviving replicates, and here is how few" rather than
    standing alone. Read off `boot.draws["beta"]`, which Stage 10 keeps when it withholds the
    interval."""
    df, ps, est, audit = fitted(subgroup_frame())
    bal = balance.assess(df, ps, audit)
    starved = bootstrap.Draws(quantity="beta", draws=np.linspace(-1.0, 1.0, 7), n_attempted=7,
                              failures={})
    boot = bootstrap.Bootstrap(config.SEED, config.N_BOOT, {"beta": starved}, {},
                               bootstrap.diagnostics(()))
    ev = sensitivity.e_value_primary(est, bal, boot, audit)
    assert ev.e_limit is None and ev.limit is None and ev.spans_null is False
    assert ev.n_draws == 7
    assert math.isfinite(ev.e_point)                # the point estimate is still reported


def test_the_approximation_is_READ_FROM_THE_FIELD_and_carries_VanderWeeles_citation():
    """The roadmap requires the conversion STATED wherever the E-value is reported, and §19c requires
    the citation in the string rather than only in a spec section: the conversion is VanderWeele's
    and attributing it to this plan's roadmap overstates whose choice it is."""
    df, ps, est, audit = fitted(subgroup_frame())
    bal = balance.assess(df, ps, audit)
    ev = sensitivity.e_value_primary(est, bal, _bootstrap_over(interval_fixtures()["spanning"], est),
                                     audit)
    assert "VanderWeele 2017" in ev.approximation and "28(6):e58" in ev.approximation
    assert "15%" in ev.approximation
    entry = audit.entry("model", "e_value_primary")
    assert ev.approximation in entry.detail or any(
        ev.approximation in cell for row in entry.table for cell in row)
    assert ev.measure == "common odds ratio [§8]"


def test_the_E_value_entry_carries_the_WORST_RESIDUAL_SMD_beside_it():
    """§7.6, and the rule is STRUCTURAL rather than editorial. The E-value bounds UNMEASURED
    confounding and this analysis has measured confounding it has not fixed, so an E-value printed
    without the worst residual |SMD| invites the reading "residual confounding would have to be
    strong" when a covariate that is IN THE MODEL is visibly not balanced."""
    df, ps, est, audit = fitted(subgroup_frame())
    bal = balance.assess(df, ps, audit)
    sensitivity.e_value_primary(est, bal, _bootstrap_over(interval_fixtures()["spanning"], est),
                                audit)
    entry = audit.entry("model", "e_value_primary")
    rows = {row[0]: row[1:] for row in entry.table[1:]}
    worst, undefined = bal.worst()
    assert rows["worst residual abs SMD"][0] == data._fmt(abs(worst.weighted))
    assert rows["rows no SMD could be computed on"][0] == str(len(undefined))
    assert worst.covariate in entry.detail
    # `bal` is a PARAMETER, so the pairing cannot be dropped by a reporting layer
    assert list(inspect.signature(sensitivity.e_value_primary).parameters) == [
        "est", "bal", "boot", "audit"]


def test_the_E_value_entry_names_EVERY_THRESHOLDS_baseline_risk_against_the_fifteen_percent_band():
    """§12.2 item 4a: Stage 14 must name the thresholds outside the conversion's documented range,
    and it cannot do that from a log that does not print them — so the entry carries one row per
    `C.MRS_THRESHOLDS` entry, read off `est.cumulative` and equal to it row for row (§15.13)."""
    df, ps, est, audit = fitted(subgroup_frame())
    bal = balance.assess(df, ps, audit)
    sensitivity.e_value_primary(est, bal, _bootstrap_over(interval_fixtures()["spanning"], est),
                                audit)
    rows = {row[0]: row[1:] for row in audit.entry("model", "e_value_primary").table[1:]}
    control = min(config.TREATMENT_LABELS)
    for k in config.MRS_THRESHOLDS:
        cell = rows[f"control P(Y <= {k})"]
        assert cell[0] == data._fmt(est.cumulative[k][control])
        above = est.cumulative[k][control] > 0.15
        assert ("NOT licensed" in cell[1]) is not above


# --- 15.7  the subgroup model, and what it may not contain --------------------------------------------

def test_the_design_is_THREE_COLUMNS_and_no_confounder_enters_it():
    """[§8] makes every estimate marginal in the overlap population: adding the [§6] confounders here
    would make `exp(beta)` a CONDITIONAL odds ratio, which is precisely what the roadmap forbids
    Stage 12 from mislabelling (§8.2). Asserted by inspecting the design the fit RAN on."""
    df, ps, est, _ = fitted(subgroup_frame())
    _, _, fit, _ = sensitivity._subgroup_fit(df, ps, est.in_estimate, SUBGROUP)
    assert fit.columns == (config.TREATMENT, SUBGROUP, f"{config.TREATMENT}_x_{SUBGROUP}")
    assert not set(fit.columns) & set(config.PS_COVARIATES) - {config.TREATMENT}


def test_G8s_predicate_IS_15_7s_assertion_and_not_a_second_spelling_of_it():
    """§8.3: a weaker predicate over set membership would pass a future `design` that REORDERED
    columns, while `beta` and `gamma` are read back BY NAME — so the two would silently stop
    describing the same fit."""
    body = ast.unparse(function_named("_subgroup_fit"))
    assert "tuple(X.columns) != declared" in body
    assert "declared = (C.TREATMENT, s, ix)" in body


def test_the_subgroup_mask_equals_in_estimate_on_this_frame_and_is_written_ANYWAY():
    """A no-op here, because both subgroup columns derive from [§6] covariates and `complete_cases`
    has already removed any row missing either. Without the mask an `Int64` `<NA>` becomes `nan`,
    reaches `polr` and raises O1 — droppable and COUNTED — so the failure would be absorbed as a
    sparse replicate rather than reported as a fact about the frame (§8.2). Asserted POSITIVELY."""
    df, ps, est, _ = fitted(subgroup_frame())
    for name in config.SUBGROUPS:
        _, _, _, mask = sensitivity._subgroup_fit(df, ps, est.in_estimate, name)
        assert mask.equals(est.in_estimate & df[name].notna())
        assert mask.equals(est.in_estimate)


def test_both_level_odds_ratios_come_from_ONE_fit_and_cannot_disagree_about_the_model(boot_n):
    df, ps, est, audit = fitted(subgroup_frame())
    subs = sensitivity.subgroups(df, ps, est, audit)
    for name, estimate in subs.estimates.items():
        beta, gamma, _, _ = sensitivity._subgroup_fit(df, ps, est.in_estimate, name)
        assert estimate.or_level[0] == float(np.exp(beta))
        assert estimate.or_level[1] == float(np.exp(beta + gamma))
        assert estimate.or_ratio == float(np.exp(gamma))
        assert estimate.gamma == gamma
        assert estimate.or_level[1] == pytest.approx(estimate.or_level[0] * estimate.or_ratio)


def test_the_level_coding_is_ARITHMETIC_and_gammas_SIGN_does_not_move_with_a_LABEL():
    """The pilot takes its reference level from `pd.get_dummies(series.astype(str), drop_first=True)`
    (`pilots/analysis.py:827`), which makes the SIGN of `gamma` a function of string sort order.
    Here level 1 is the named condition and the column is arithmetic, so relabelling cannot move it
    (§8.2). Asserted by building a frame whose labels sort the other way."""
    df, ps, est, _ = fitted(subgroup_frame())
    _, gamma, _, _ = sensitivity._subgroup_fit(df, ps, est.in_estimate, SUBGROUP)

    relabelled = df.copy()
    relabelled["onset_type"] = relabelled["onset_type"].map(
        {"witnessed": "zzz_witnessed", "unwitnessed": "aaa_unwitnessed",
         "wake_up": "aaa_wake_up"}).astype("string")
    # the SUBGROUP column is unchanged — it is arithmetic and does not read the label at all
    assert relabelled[SUBGROUP].equals(df[SUBGROUP])
    _, again, _, _ = sensitivity._subgroup_fit(relabelled, ps, est.in_estimate, SUBGROUP)
    assert again == gamma


def test_NOTHING_asserts_that_the_two_level_odds_ratios_BRACKET_the_primarys():
    """§8.7's prohibition. A weighted proportional-odds model is NOT collapsible, so the primary
    `exp(beta)` is not a weighted average of the two and they need not bracket it. A reviewer's
    instinct is to add `assert min(or_level) <= primary_or <= max(or_level)` as a sanity check, and
    it would fail on correct output — so the module and this file are scanned for the comparison."""
    # By AST, and over BOTH files: any chained comparison — `a <= b <= c` — in which a subgroup odds
    # ratio and a primary one both appear is the assertion this prohibition is about. A literal scan
    # cannot be written here, because the literal would appear in the scan itself.
    for source in (SOURCE, THIS_FILE):
        for node in ast.walk(ast.parse(source)):
            if not isinstance(node, ast.Compare) or len(node.ops) < 2:
                continue
            spelled = ast.unparse(node)
            assert not ("or_level" in spelled and "odds_ratio" in spelled), spelled
    # ...and there is no field inviting it, which is where the comparison would otherwise live.
    assert "primary_odds_ratio" not in {f for f in vars(sensitivity.SubgroupEstimate)}
    assert "odds_ratio" not in sensitivity.SubgroupEstimate.__annotations__


def test_the_six_cumulative_RD_k_are_NOT_reported_per_level_and_the_DISTRIBUTIONS_are(boot_n):
    """§8.6. [§13] does not ask for them; twelve more intervals from a hypothesis-generating analysis
    is twelve more numbers a reader will read as findings. What IS reported is the pair of weighted
    cumulative distributions per level, POINT VALUES ONLY, no interval and no p."""
    df, ps, est, audit = fitted(subgroup_frame())
    subs = sensitivity.subgroups(df, ps, est, audit)
    for estimate in subs.estimates.values():
        assert set(estimate.cumulative) == {0, 1}
        for level, distribution in estimate.cumulative.items():
            assert set(distribution) == set(config.MRS_THRESHOLDS)
            for by_arm in distribution.values():
                assert set(by_arm) == set(config.TREATMENT_LABELS)
    assert not any(key.startswith("rd_") or ".rd" in key for key in subs.draws)
    assert not any(key.startswith("rd_") or ".rd" in key for key in subs.intervals)


def test_the_six_estimand_keys_are_COMPUTED_from_the_registry_and_never_listed(boot_n):
    """§3.3, and it is Stage 9 §4.1's move two stages on: a third subgroup restored to the registry
    gains its three keys, its interval and its p in the same edit, and cannot fail to."""
    df, ps, est, audit = fitted(subgroup_frame())
    subs = sensitivity.subgroups(df, ps, est, audit)
    assert set(subs.draws) == {f"{s}.{q}" for s in config.SUBGROUPS
                               for q in ("beta", "gamma", "beta_plus_gamma")}
    assert len(subs.draws) == 3 * len(config.SUBGROUPS)
    for estimate in subs.estimates.values():
        assert estimate.n_interaction_tests == len(config.SUBGROUPS)
        assert estimate.hypothesis_generating is True
    assert "n_interaction_tests=len(C.SUBGROUPS)" in SOURCE.replace(" ", "")


def test_a_p_value_is_reported_on_GAMMA_and_on_NOTHING_ELSE(boot_n):
    """[§13] asks for AN interaction test; two level-specific p-values would be two tests of one
    contrast, which is [§8]'s own objection to six threshold-wise tests one stage on (§8.5)."""
    df, ps, est, audit = fitted(subgroup_frame())
    subs = sensitivity.subgroups(df, ps, est, audit)
    tested = {key for key, interval in subs.intervals.items() if interval.p is not None}
    assert tested == {f"{s}.gamma" for s in config.SUBGROUPS}
    assert all(sensitivity._tested_subgroup(key) is (key in tested) for key in subs.draws)
    # a FUNCTION and not a frozen set, so a third subgroup gains its p in the same edit
    assert sensitivity._tested_subgroup("a_third_subgroup.gamma") is True
    assert sensitivity._tested_subgroup("a_third_subgroup.beta") is False


def test_subgroups_carries_NO_diagnostics_field_and_that_is_the_decision(boot_n):
    """`Replicate.n_alpha` is a SCALAR and this body runs TWO `polr` fits, so any single value would
    silently describe one subgroup while being tallied as the replicate's (§3.1, §8.5). Every counter
    `subgroup_replicates` prints therefore comes off `Draws`."""
    df, ps, est, audit = fitted(subgroup_frame())
    subs = sensitivity.subgroups(df, ps, est, audit)
    assert not hasattr(subs, "diagnostics")
    assert "bootstrap.diagnostics" not in ast.unparse(function_named("subgroups"))
    body = ast.unparse(function_named("_subgroup_replicate"))
    assert "n_alpha=None" in body and "polr_iterations=None" in body
    # ...and `sum_w` and `n_in_model` ARE populated, because there is exactly one propensity fit
    assert "sum_w=float(ps.w[ps.in_model].sum())" in body


# --- 15.8  the two degeneracies, and the bucket map is scanned not trusted ----------------------------

@pytest.mark.parametrize("kind,token,bucket", [("empty_cell", "O6", "degenerate_design"),
                                               ("constant_product", "G8", "degenerate_design")])
def test_each_degeneracy_is_DROPPED_AND_COUNTED_inside_the_loop(kind, token, bucket, boot_n):
    """[roadmap invariant 5]: never substituted with a different estimator, and never silently
    omitted. The token classifies into the bucket §9 declares, and `bucket` RAISES on an
    unrecognised one rather than defaulting — so a token added without its map entry is a crash and
    not a counter reading zero."""
    df = separable_subgroup_frame(kind)
    ps, = (propensity.fit(df, data.Audit(data.WORKBOOK)),)
    replicate = sensitivity._subgroup_replicate(
        df, tuple(f"{s}.{q}" for s in config.SUBGROUPS
                  for q in ("beta", "gamma", "beta_plus_gamma")), data.WORKBOOK)
    assert config.FAILURE_BUCKETS[token] == bucket
    assert all(v == bucket for k, v in replicate.failures.items() if k.startswith(f"{SUBGROUP}."))


def test_G9_is_DROPPED_AND_COUNTED_as_separation_and_not_as_nonconvergence(boot_n):
    """Stage 8 §11 requires the separation count reported SEPARATELY from the convergence count,
    because it measured that the second never fires on this estimator and a single counter reading
    zero therefore says nothing. G9 is G7's guard on a different set of coefficients (§9)."""
    df = near_separated_subgroup_frame()
    replicate = sensitivity._subgroup_replicate(
        df, tuple(f"{s}.{q}" for s in config.SUBGROUPS
                  for q in ("beta", "gamma", "beta_plus_gamma")), data.WORKBOOK)
    assert config.FAILURE_BUCKETS["G9"] == "separation"
    assert replicate.failures[f"{SUBGROUP}.gamma"] == "separation"


def test_THE_ATOMICITY_one_subgroup_fails_and_the_other_carries_DRAWS_in_ONE_replicate():
    """§8.5's granularity rule as an assertion rather than a docstring: a subgroup's own failure
    costs that subgroup's THREE keys atomically and costs the other subgroup NOTHING."""
    df = separable_subgroup_frame("empty_cell")           # breaks `unknown_onset` only
    keys = tuple(f"{s}.{q}" for s in config.SUBGROUPS
                 for q in ("beta", "gamma", "beta_plus_gamma"))
    replicate = sensitivity._subgroup_replicate(df, keys, data.WORKBOOK)
    other = next(s for s in config.SUBGROUPS if s != SUBGROUP)

    assert set(replicate.failures) == {f"{SUBGROUP}.{q}"
                                       for q in ("beta", "gamma", "beta_plus_gamma")}
    assert set(replicate.values) == {f"{other}.{q}"
                                     for q in ("beta", "gamma", "beta_plus_gamma")}
    assert not set(replicate.values) & set(replicate.failures)     # Replicate's own partition


def test_a_PROPENSITY_failure_fails_the_WHOLE_replicate_and_is_counted_against_EVERY_key():
    """§7.1, inherited: nothing downstream runs without `e`, `w` and `in_model`, so there is no
    estimand the failure does not affect — and a reading in which those replicates are never
    ATTEMPTED is precisely the silent drop `Draws.__post_init__` exists to make impossible."""
    df = subgroup_frame()
    keys = tuple(f"{s}.{q}" for s in config.SUBGROUPS
                 for q in ("beta", "gamma", "beta_plus_gamma"))
    one_armed = df[df[config.TREATMENT] == 1.0]
    with pytest.raises(config.SchemaError):
        sensitivity._subgroup_replicate(one_armed, keys, data.WORKBOOK)   # F2 is a SchemaError

    # Two [§6] covariates made exactly collinear, so `model.firth` raises F3 on a rank deficiency —
    # `model.FitError`, which IS droppable, and NOT a `SchemaError`, which is not.
    collinear = df.copy()
    collinear["tmax6_ml"] = collinear["core_ml"]
    with pytest.raises(model.FitError, match="F3"):
        propensity.fit(collinear, data.Audit(data.WORKBOOK))
    replicate = sensitivity._subgroup_replicate(collinear, keys, data.WORKBOOK)
    assert set(replicate.failures) == set(keys) and replicate.values == {}
    assert replicate.n_alpha is None and replicate.sum_w is None and replicate.n_in_model is None


def test_no_H_identifier_is_a_FAILURE_BUCKETS_key_and_every_FitError_token_is():
    """§9.1's separation, from BOTH sides. `bucket` raises on an unrecognised leading token rather
    than defaulting, so an H token arriving there would present as a taxonomy gap — it cannot,
    because `C.SchemaError` is never caught anywhere in this module and the only `except` in either
    replicate body is on `model.FitError`."""
    tokens = {
        str(node.exc.args[0].values[0].value if isinstance(node.exc.args[0], ast.JoinedStr)
            else node.exc.args[0].value).split()[0]
        for node in ast.walk(ast.parse(SOURCE))
        if isinstance(node, ast.Raise) and isinstance(node.exc, ast.Call)
        and getattr(node.exc.func, "attr", getattr(node.exc.func, "id", "")) == "FitError"}
    assert tokens == {"G8", "G9"}
    assert tokens <= set(config.FAILURE_BUCKETS)

    identifiers = {f"H{i}" for i in range(1, 11)}
    assert not identifiers & set(config.FAILURE_BUCKETS)
    handlers = [ast.unparse(node.type) for node in ast.walk(ast.parse(SOURCE))
                if isinstance(node, ast.ExceptHandler) and node.type is not None]
    assert set(handlers) == {"model.FitError"}


# --- 15.9  the three runs draw the same frames, and the two dicts are never merged --------------------

def test_all_three_runs_draw_the_SAME_SEQUENCE_of_frames_and_a_HAND_DRIVEN_loop_matches(boot_n):
    """Same seed, same frame, same stratum. `replicates` creates ONE `default_rng(seed)` and calls
    `resample` OUTSIDE `body` with no `try`, so the stream is a function of the seed alone and not of
    what the body does or raises (§5.4, §8.5). The pairing is ASSERTED AND NEVER USED: nothing
    computes a contrast between an arm draw and a primary draw, so no false precision arises."""
    df, ps, est, audit = fitted(subgroup_frame())
    sec = outcome.secondary(df, ps, audit)
    subgroup_keys = tuple(f"{s}.{q}" for s in config.SUBGROUPS
                          for q in ("beta", "gamma", "beta_plus_gamma"))
    arm_keys = ("beta", *(f"rd_{k}" for k in est.rd))
    paths = {key: e.augmented_path for key, e in sec.estimates.items()}
    stage_ten_keys = bootstrap._estimand_keys(est, sec)

    def fingerprints(body):
        seen: list[tuple[str, ...]] = []

        def spy(draw):
            seen.append(tuple(draw["case_id"]))
            return body(draw)

        bootstrap.replicates(df, spy, config.N_BOOT, config.SEED, config.BOOT_STRATUM)
        return seen

    hand = []
    rng = np.random.default_rng(config.SEED)
    for _ in range(config.N_BOOT):
        hand.append(tuple(bootstrap.resample(df, rng, config.BOOT_STRATUM)["case_id"]))

    stage_ten = fingerprints(
        lambda draw: bootstrap._replicate(draw, stage_ten_keys, paths, data.WORKBOOK))
    arm = fingerprints(lambda draw: sensitivity._arm_replicate(draw, arm_keys, data.WORKBOOK))
    subgroups = fingerprints(
        lambda draw: sensitivity._subgroup_replicate(draw, subgroup_keys, data.WORKBOOK))

    assert stage_ten == arm == subgroups == hand
    assert len(set(hand)) > 1                        # the frames really do differ from one another


def test_the_TWO_DRAWS_DICTS_ARE_NEVER_MERGED_and_the_test_MEASURES_what_the_merge_costs(boot_n):
    """§3.3's prohibition, as a measured fact in the suite rather than a warning in a docstring.

    The arm's keys are DELIBERATELY the primary's, so the two tables are compared key for key — which
    is what "reported beside the primary ones" means when the reader is a person and not a formatter.
    The cost is that `{**boot.draws, **arm.draws}` is a SILENT OVERWRITE: seven keys collide, nothing
    in `Draws` notices, and the result is a dict of 26 keys in which seven describe a different
    specification.
    """
    df, ps, est, audit = fitted(subgroup_frame())
    sec = outcome.secondary(df, ps, audit)
    boot = bootstrap.run(df, ps, est, sec, audit)
    arm = sensitivity.full_covariate(df, ps, audit)

    assert set(arm.draws) < set(boot.draws)                       # every arm key is a primary key
    merged = {**boot.draws, **arm.draws}
    assert len(merged) == len(boot.draws) == 26                   # 26 and NOT 33
    assert len(boot.draws) + len(arm.draws) == 33
    for key in arm.draws:
        assert merged[key] is arm.draws[key]                      # ...and seven now describe [§13]
    assert "{**boot.draws" not in SOURCE and "**arm.draws" not in SOURCE


# --- 15.11  the module boundary, and the public surface is six names ----------------------------------

def test_the_public_surface_is_SIX_NAMES():
    """§3.2's count, re-derived rather than repeated."""
    functions = [n.name for n in ast.parse(SOURCE).body if isinstance(n, ast.FunctionDef)]
    public = [n for n in functions if not n.startswith("_")]
    assert public == ["benjamini_hochberg", "e_value", "multiplicity", "e_value_primary",
                      "subgroups", "full_covariate"]


def test_the_six_privates_that_matter_are_all_present_and_there_is_NO_second_bucket_or_collect():
    """§5.4 makes Stage 10's three public rather than writing second ones here: a second classifier
    would be a second definition of the separation-versus-convergence split, a second `Draws` loop a
    second implementation of the reconciliation property, and a second tally a second answer to "how
    many replicates fitted five cutpoints"."""
    functions = {n.name for n in ast.parse(SOURCE).body if isinstance(n, ast.FunctionDef)}
    assert {"_arm_replicate", "_subgroup_fit", "_subgroup_replicate", "_assert_arm_inputs",
            "_tested_subgroup", "_intervals"} <= functions
    assert not {"_bucket", "_collect", "_diagnostics"} & functions
    calls = attribute_calls(SOURCE)
    assert {"bootstrap.bucket", "bootstrap.collect", "bootstrap.diagnostics"} <= calls


@pytest.mark.parametrize("forbidden", ["statsmodels", "scipy", "sklearn"])
def test_NO_SHIPPED_MODULE_IMPORTS_THE_ORACLE(forbidden):
    """Stage 6 §12.7's scan, EXTENDED rather than duplicated. `pyproject.toml` declares statsmodels
    *"retained for unpenalised cross-checks in tests, not for any reported estimate"* — and an
    adjusted p-value IS a reported estimate, so importing `multipletests` into a shipped module to
    produce one is the thing that policy exists to forbid (§6.4)."""
    root = Path(config.__file__).resolve().parent
    for path in sorted(root.glob("*.py")):
        imported: set[str] = set()
        for node in ast.walk(ast.parse(path.read_text(encoding="utf-8"))):
            if isinstance(node, ast.Import):
                imported.update(alias.name.split(".")[0] for alias in node.names)
            elif isinstance(node, ast.ImportFrom) and node.module:
                imported.add(node.module.split(".")[0])
        assert forbidden not in imported, f"{path.name} imports {forbidden}"


def test_sensitivity_names_NO_outcome_private_and_NEITHER_derive_entry_point():
    """`derive_cohort` and `_core_above_median` are the [§13] median subgroup's own definition, and
    the whole point of freezing it is that a replicate RESAMPLES the column rather than recomputing
    it — a per-replicate cut-point would give every replicate its own subgroup, silently (§8.2)."""
    # By AST and never by substring, which is Stage 10 §15.11's instrument: this module DISCUSSES
    # both names in a comment — the argument for why the median is resampled rather than recomputed
    # is exactly where they belong — and what must not happen is that either is CALLED.
    named: set[str] = set()
    for node in ast.walk(ast.parse(SOURCE)):
        if isinstance(node, ast.Name):
            named.add(node.id)
        elif isinstance(node, ast.Attribute):
            named.add(node.attr)
    assert not {c for c in attribute_calls(SOURCE) if c.startswith("outcome._")}
    assert "derive" not in {n.split(".")[0] for n in attribute_calls(SOURCE)}
    for forbidden in ("derive_cohort", "_core_above_median"):
        assert forbidden not in named


@pytest.mark.parametrize("function", ["benjamini_hochberg", "e_value"])
def test_the_TWO_DEFINITIONS_name_NO_outcome_NO_family_NO_covariate_and_NO_centre(function):
    """§3.2: each is a DEFINITION rather than a step, and a definition with one caller is still a
    definition somebody will want to check against an oracle without constructing a frame."""
    names: set[str] = set()
    for child in ast.walk(function_named(function)):
        if isinstance(child, ast.Name):
            names.add(child.id)
        elif isinstance(child, ast.Attribute):
            names.add(child.attr)
        elif isinstance(child, ast.Constant) and isinstance(child.value, str):
            names.add(child.value)
    forbidden = (set(config.PS_COVARIATES) | set(config.OUTCOMES) | set(config.CENTER_ORDER)
                 | set(config.SUBGROUPS) | {config.TREATMENT, "BINARY_OUTCOMES", "SUBGROUPS",
                                            "PRIMARY_OUTCOME", "MRS_THRESHOLDS"})
    assert not forbidden & names, f"{function} names {sorted(forbidden & names)}"


def test_the_four_returned_records_carry_NO_VERDICT_OF_ANY_KIND():
    """Stage 8 §3's rule binds harder here than anywhere before it, because three of the four records
    are the ones a reader most wants a verdict from: an adjusted p, an E-value and an interaction
    test (§3.1)."""
    for record in (sensitivity.FamilyCorrection, sensitivity.Multiplicity, sensitivity.EValue,
                   sensitivity.SubgroupEstimate, sensitivity.Subgroups, sensitivity.Arm):
        for name in dir(record):
            assert not name.startswith("favours")
            assert name not in {"significant", "verdict", "is_significant", "reject"}


# --- 15.13  the fifteen audit entries, their step names, and the ledger -------------------------------

STEP_NAMES: Final[tuple[str, ...]] = (
    "multiplicity",
    "e_value_primary",
    *(f"subgroup_fit_{s}" for s in config.SUBGROUPS),
    "subgroup_replicates",
    "covariate_completeness_full_covariate",
    "design_matrix_full_covariate",
    "propensity_fit_full_covariate",
    "overlap_weights_full_covariate",
    "balance_smd_full_covariate",
    "overlap_by_centre_full_covariate",
    "outcome_completeness_full_covariate",
    "primary_fit_full_covariate",
    "cumulative_rd_full_covariate",
    "sensitivity_arm",
)


@pytest.fixture(scope="module")
def staged(stage11_run):
    """The four entry points driven once, in §14's order, against the SAME Audit."""
    df, ps, bal, est, sec, audit, boot = stage11_run
    original = config.N_BOOT
    config.N_BOOT = RUN_BOOT
    try:
        before = len(audit.entries)
        mult = sensitivity.multiplicity(sec, boot, audit)
        ev = sensitivity.e_value_primary(est, bal, boot, audit)
        subs = sensitivity.subgroups(df, ps, est, audit)
        arm = sensitivity.full_covariate(df, ps, audit)
        yield df, audit, before, mult, ev, subs, arm
    finally:
        config.N_BOOT = original


def test_the_stage_appends_FIFTEEN_entries_with_EXACTLY_these_step_names(staged):
    _, audit, before, *_ = staged
    added = audit.entries[before:]
    assert len(added) == 15
    assert [entry.step for entry in added] == list(STEP_NAMES)


def test_all_fifteen_are_the_MODEL_kind_and_data_KINDS_stays_NINE(staged):
    """No new `data.KINDS` entry: all fifteen are claims about what was fitted to a population, which
    is the kind `data.py` declares for exactly that. KINDS has been nine for four stages (§10)."""
    _, audit, before, *_ = staged
    assert {entry.kind for entry in audit.entries[before:]} == {"model"}
    assert len(data.KINDS) == 9


def test_case_ids_is_EMPTY_on_all_fifteen(staged):
    """Nine `model` entries name cases and these remove no patient; the arm's own excluded record is
    named by `propensity._record_exclusion` under the SUFFIXED step, which is that function's rule
    applying to the arm unchanged (§10)."""
    _, audit, before, *_ = staged
    named = {entry.step: entry.case_ids for entry in audit.entries[before:]
             if entry.step != "covariate_completeness_full_covariate"}
    assert all(ids == () for ids in named.values()), named


def test_the_step_name_COLLISION_has_been_seen_and_the_six_suffixed_are_DISJOINT(staged):
    """Definition of done item 3. One `Audit`, both specifications, and
    `audit.entry("model", "propensity_fit")` returning the PRIMARY's."""
    df, audit, *_ = staged
    suffixed = {step for step in STEP_NAMES if step.endswith("_full_covariate")}
    unsuffixed = {step[:-len("_full_covariate")] for step in suffixed}
    # NINE and not §15.13's six: `propensity._fit` records four, `balance.assess` two and
    # `outcome.primary` three, and §10's own list enumerates all nine. §15.13's "six" predates the
    # count and §21 is normative — so the number here is derived from §10's list rather than quoted.
    assert len(suffixed) == 9
    assert not suffixed & unsuffixed
    for step in unsuffixed:
        assert audit.entry("model", step) is not None            # the primary's, and only its
    primary = audit.entry("model", "propensity_fit")
    arm = audit.entry("model", "propensity_fit_full_covariate")
    assert primary.n != arm.n or primary.detail != arm.detail
    assert config.PROPENSITY_PRIMARY.label in primary.detail
    assert config.PROPENSITY_FULL.label in arm.detail


def test_the_arm_entry_carries_BOTH_specifications_and_its_counters_come_off_the_RECORD(staged):
    """§10's first decision: the Accept-when is about a PAIR, and an entry showing one half makes the
    comparison the reader's arithmetic."""
    _, audit, _, _, _, _, arm = staged
    entry = audit.entry("model", "sensitivity_arm")
    assert entry.n == int(arm.ps.in_model.sum())
    assert entry.table[0][1] == config.PROPENSITY_PRIMARY.label
    assert entry.table[0][2] == config.PROPENSITY_FULL.label
    rows = {row[0]: row[1:] for row in entry.table[1:]}
    assert rows["in_model"][:2] == (str(int(arm.reference.in_model.sum())),
                                    str(int(arm.ps.in_model.sum())))
    for code, label in config.TREATMENT_LABELS.items():
        assert rows[f"ESS ({label})"][:2] == (data._fmt(arm.reference.ess[code]),
                                              data._fmt(arm.ps.ess[code]))
    for key, d in arm.draws.items():
        assert rows[f"replicates surviving: {key}"][1] == f"{len(d.draws)} of {d.n_attempted}"
    # the clause is SELECTED from the computed field and never written as prose (§5.2)
    assert rows["the two in_model masks are IDENTICAL"][0] == (
        "yes" if arm.differs_only_in_specification else "no")


def test_the_arm_entry_NAMES_THE_SPENT_NEGATIVE_CONTROLS_and_prints_BOTH_DIRECTIONS(staged):
    """§5.3: a lower worst residual |SMD| is NOT uniformly better balance, and the summary statistic
    alone reads as though it were."""
    _, audit, *_ = staged
    entry = audit.entry("model", "sensitivity_arm")
    rows = {row[0]: row[1:] for row in entry.table[1:]}
    for name in config.NEGATIVE_CONTROLS:
        assert f"abs SMD after weighting: {name}" in rows
    assert any(key.startswith("abs SMD after weighting: center = ") for key in rows)
    assert "worst residual abs SMD" in rows
    assert "spends the negative controls" in entry.detail.lower() or "spends" in entry.detail
    assert "gets WORSE" in entry.detail


def test_the_subgroup_entries_print_the_FOUR_CELLS_and_the_treated_count_per_level(staged):
    """§10: the number that governs §8.3 and §8.4 both is visible to a reader of the log WITHOUT the
    specification."""
    _, audit, _, _, _, subs, _ = staged
    for name, estimate in subs.estimates.items():
        entry = audit.entry("model", f"subgroup_fit_{name}")
        rows = {row[0]: row[1:] for row in entry.table[1:]}
        for cell, count in estimate.n_by_level_arm.items():
            assert rows[f"n {cell}"][0] == str(count)
        assert "hypothesis-generating" in entry.detail.lower()
        assert str(estimate.n_interaction_tests) in rows["interaction tests performed"][0]


def test_the_subgroup_replicates_entry_reports_off_DRAWS_and_names_the_TRUNCATION(staged):
    """§10: per key, `n_attempted`, the surviving count and the bucket counts — and the sentence the
    drop rate needs, which is NOT the same sentence as a count of lost replicates (§8.4)."""
    _, audit, _, _, _, subs, _ = staged
    entry = audit.entry("model", "subgroup_replicates")
    rows = {row[0]: row[1:] for row in entry.table[1:]}
    for key, d in subs.draws.items():
        assert rows[key][0] == str(d.n_attempted)
        assert rows[key][1] == str(len(d.draws))
    assert "in_model per replicate" in rows
    assert "SELECTION ON THE ESTIMATE" in entry.detail
    assert "narrower" in entry.detail


def test_the_multiplicity_entry_prints_the_PRIMARY_UNCORRECTED_beside_the_families(staged):
    _, audit, _, mult, *_ = staged
    entry = audit.entry("model", "multiplicity")
    rows = {row[0]: row[1:] for row in entry.table[1:]}
    assert rows[config.PRIMARY_OUTCOME][0] == "primary"
    assert rows[config.PRIMARY_OUTCOME][3] == data._fmt(mult.primary_uncorrected)
    assert "uncorrected" in rows[config.PRIMARY_OUTCOME][4]
    for family, correction in mult.families.items():
        for name in correction.raw:
            assert rows[name][0] == family
            assert rows[name][1] == str(correction.m_declared)
            assert rows[name][2] == str(correction.m_used)
            assert rows[name][3] == data._fmt(correction.raw[name])
            assert rows[name][4] == data._fmt(correction.adjusted[name])
            assert ("descriptive only" in rows[name][5]) is correction.descriptive_only


def test_NO_cell_of_any_new_grid_contains_a_PIPE_and_every_row_has_the_SAME_WIDTH(staged):
    """`data._md_table` computes its column widths from `rows[0]`, so a short row raises IndexError
    and a pipe inside a cell breaks the rendering — `propensity._padded`'s own finding."""
    _, audit, before, *_ = staged
    for entry in audit.entries[before:]:
        assert entry.table is not None
        widths = {len(row) for row in entry.table}
        assert len(widths) == 1, (entry.step, sorted(widths))
        assert not any("|" in cell for row in entry.table for cell in row), entry.step
    rendered = audit.to_markdown()
    assert "sensitivity_arm" in rendered and "subgroup_replicates" in rendered


# --- 15.12  stage 11 adds nothing, refits nothing, and breaks nothing earlier -------------------------

def test_the_landed_pipeline_is_BYTE_IDENTICAL_with_and_without_this_stage():
    """§15.12, and it is how `outcome.py`'s estimator is asserted unchanged — BY NUMBER rather than
    by diff, because four things in it do change: three audit step names and the extraction of
    `estimation_population`. The extraction is the one that could move a number and does not: if it
    were not equivalent, every downstream number would move at once and none of them does.
    """
    def landed():
        df = subgroup_frame()
        audit = data.Audit(data.WORKBOOK)
        ps = propensity.fit(df, audit)
        bal = balance.assess(df, ps, audit)
        est = outcome.primary(df, ps, audit)
        sec = outcome.secondary(df, ps, audit)
        return df, ps, bal, est, sec, audit

    original = config.N_BOOT
    config.N_BOOT = RUN_BOOT
    try:
        df_a, ps_a, bal_a, est_a, sec_a, audit_a = landed()
        boot_a = bootstrap.run(df_a, ps_a, est_a, sec_a, audit_a)
        alone = audit_a.to_markdown()

        df_b, ps_b, bal_b, est_b, sec_b, audit_b = landed()
        boot_b = bootstrap.run(df_b, ps_b, est_b, sec_b, audit_b)
        ledger = len(audit_b.entries)
        sensitivity.multiplicity(sec_b, boot_b, audit_b)
        sensitivity.e_value_primary(est_b, bal_b, boot_b, audit_b)
        sensitivity.subgroups(df_b, ps_b, est_b, audit_b)
        sensitivity.full_covariate(df_b, ps_b, audit_b)
    finally:
        config.N_BOOT = original

    # every landed number
    assert est_a.beta == est_b.beta and est_a.rd == est_b.rd
    assert list(ps_a.w.dropna()) == list(ps_b.w.dropna())
    assert [row.weighted for row in bal_a.covariates] == [row.weighted for row in bal_b.covariates]
    assert {k: list(v.draws) for k, v in boot_a.draws.items()} == {
        k: list(v.draws) for k, v in boot_b.draws.items()}
    # ...and the entries the landed stages recorded, byte for byte
    assert [(e.kind, e.step, e.n, e.detail, e.case_ids, e.table)
            for e in audit_a.entries] == [(e.kind, e.step, e.n, e.detail, e.case_ids, e.table)
                                          for e in audit_b.entries[:ledger]]
    assert alone == data.Audit(data.WORKBOOK).__class__.to_markdown(audit_a)


# --- 15.14  `[data-gated]` — the workbook, end to end -------------------------------------------------

@pytest.fixture(scope="module")
def workbook():
    """The whole Stage 1-10 pipeline over the private workbook, at the PRESPECIFIED replicate count
    for the point estimates and at RUN_BOOT for the loops.

    The replicate count is moved for the reason `test_bootstrap.py` moves it: at 2000 this fixture is
    a nine-minute drive, and every property §15.14 asserts is a property of the shape rather than of
    the length. §14's own numbers are produced by the shipped pipeline and live in the gitignored
    log, never here (§1).
    """
    import cohort
    import derive
    import eligibility
    original = config.N_BOOT
    config.N_BOOT = RUN_BOOT
    try:
        df, audit = data.load(data.WORKBOOK)
        df = derive.derive(df, audit)
        df = eligibility.classify(df, audit)
        df = cohort.build(df, audit)
        ps = propensity.fit(df, audit)
        bal = balance.assess(df, ps, audit)
        est = outcome.primary(df, ps, audit)
        sec = outcome.secondary(df, ps, audit)
        yield df, ps, bal, est, sec, audit, bootstrap.run(df, ps, est, sec, audit)
    finally:
        config.N_BOOT = original


@pytest.fixture(scope="module")
def workbook_staged(workbook):
    """The four entry points driven ONCE over the workbook, in §14's order, against its own Audit.

    Module-scoped and shared, because every `[data-gated]` assertion below is about the SAME drive —
    a second one would double a two-minute fixture to assert the same shapes, and a per-test drive
    over a shared `Audit` would make the ledger count depend on test ORDER.
    """
    df, ps, bal, est, sec, audit, boot = workbook
    original = config.N_BOOT
    config.N_BOOT = RUN_BOOT
    try:
        before = len(audit.entries)
        mult = sensitivity.multiplicity(sec, boot, audit)
        ev = sensitivity.e_value_primary(est, bal, boot, audit)
        subs = sensitivity.subgroups(df, ps, est, audit)
        arm = sensitivity.full_covariate(df, ps, audit)
        yield before, mult, ev, subs, arm
    finally:
        config.N_BOOT = original


@DATA_GATED
def test_all_four_entry_points_run_on_the_workbook_with_ZERO_SchemaError(workbook, workbook_staged):
    """`[data-gated]`. Definition of done item 7, and no number this test produces appears in this
    file, in the specification, or anywhere under version control (§1)."""
    df, ps, bal, est, sec, audit, boot = workbook
    _, mult, ev, subs, arm = workbook_staged

    # the family sizes are three and four, off `by_family()` and never off `intervals`
    assert {f: c.m_declared for f, c in mult.families.items()} == {"secondary": 3, "safety": 4}
    assert all(c.absent == () for c in mult.families.values())

    # every surviving count clears the floor, so every estimand carries an interval
    for record in (subs, arm):
        assert set(record.intervals) == set(record.draws)
        for d in record.draws.values():
            assert len(d.draws) >= config.ci_min_draws()

    # the four re-roled rows are EXACTLY C.NEGATIVE_CONTROLS, and the table does not move
    primary_rows = {row.covariate: row.role for row in bal.covariates}
    arm_rows = {row.covariate: row.role for row in arm.balance.covariates}
    assert list(primary_rows) == list(arm_rows)
    assert {name for name in primary_rows if primary_rows[name] != arm_rows[name]} == set(
        config.NEGATIVE_CONTROLS)

    # the two in_model masks are equal, so the roadmap's Accept-when holds ON THIS WORKBOOK
    assert ps.in_model.equals(arm.ps.in_model)
    assert arm.differs_only_in_specification is True

    # the E-value is finite and the interaction tests are counted rather than corrected
    assert math.isfinite(ev.e_point)
    assert all(e.n_interaction_tests == len(config.SUBGROUPS) for e in subs.estimates.values())
    assert set(subs.estimates) == set(config.SUBGROUPS)


@DATA_GATED
def test_the_ledger_reaches_FORTY_SIX_on_the_workbook(workbook, workbook_staged):
    """`[data-gated]`. 31 + 15. A ledger count derived by arithmetic rather than by running is
    exactly the kind Stage 10 §10.1 required be asserted (§10)."""
    _, _, _, _, _, audit, _ = workbook
    before, *_ = workbook_staged
    assert before == 31
    assert len(audit.entries) == 46
    assert [e.step for e in audit.entries[before:]] == list(STEP_NAMES)
    assert {e.kind for e in audit.entries[before:]} == {"model"}


@DATA_GATED
def test_the_five_treated_at_witnessed_onset_is_a_MEASURED_FACT_of_this_cohort(workbook_staged):
    """`[data-gated]`. §16 item 2: [§13]'s amendment of 2026-08-10 justified keeping this subgroup on
    the cohort split, and the number that governs ESTIMABILITY is a different one — the treated count
    in the ATO population at the level. It is asserted to be SMALL rather than to a value, because
    §1 keeps this stage's numbers out of version control and because the PI question it raises is
    about the shape and not about the digit."""
    _, _, _, subs, _ = workbook_staged
    cells = subs.estimates["unknown_onset"].n_by_level_arm
    treated_at_witnessed = cells[f"S=0,A={max(config.TREATMENT_LABELS)}"]
    assert 0 < treated_at_witnessed < 10
    assert sum(cells.values()) == subs.estimates["unknown_onset"].n


# --- the Definition of done's mutation companions -----------------------------------------------------

def test_WITHOUT_THE_RUNNING_MIN_and_WITHOUT_THE_UNSORT_the_step_up_IS_SEEN_TO_FAIL():
    """Definition of done item 4. Both wrong implementations are WRITTEN OUT and driven on §6.3's
    vector, so the two properties are measured facts in the suite rather than assertions about a
    function nobody wrote.

    One input, three distinct wrong answers — and each is silent in a different way: without the
    running min the middle outcome reports a LARGER adjusted p than the one below it and two of three
    secondary outcomes swap rank; without the unsort every p is reattached to the wrong outcome.
    """
    def without_running_min(p):
        m, order = p.size, np.argsort(p, kind="stable")
        q = np.minimum(1.0, p[order] * m / np.arange(1, m + 1))
        out = np.empty_like(q)
        out[order] = q
        return out

    def without_unsort(p):
        m, order = p.size, np.argsort(p, kind="stable")
        q = np.minimum(1.0, p[order] * m / np.arange(1, m + 1))
        return np.minimum.accumulate(q[::-1])[::-1]

    def without_both(p):
        m, order = p.size, np.argsort(p, kind="stable")
        return np.minimum(1.0, p[order] * m / np.arange(1, m + 1))

    p = np.array([0.04, 0.01, 0.03])
    right = list(sensitivity.benjamini_hochberg(p))
    assert right == [0.04, 0.03, 0.04]
    assert list(without_running_min(p)) == [0.04, 0.03, 0.045]      # NOT monotone in rank
    assert list(without_unsort(p)) == [0.03, 0.04, 0.04]            # the family, SORTED
    assert list(without_both(p)) == [0.03, 0.045, 0.04]
    for wrong in (without_running_min, without_unsort, without_both):
        assert list(wrong(p)) != right

    # ...and the companion case is the one ALL THREE also get right, which is why a suite consisting
    # only of it cannot look like coverage (§6.3).
    companion = np.array([0.01, 0.02, 0.03])
    assert list(without_running_min(companion)) == list(sensitivity.benjamini_hochberg(companion))


def test_THE_TWO_DOORS_INTO_subgroup_fit_behave_DIFFERENTLY_ON_ONE_FRAME(boot_n):
    """Definition of done item 12, on ONE frame rather than on two.

    Through the replicate loop the same degeneracy is `model.FitError`, is DROPPED and is COUNTED
    into `C.FAILURE_BUCKETS`' bucket; through `subgroups` it is `C.SchemaError` H8 and stops the
    stage. The droppable-and-counted argument is about the replicate loop, and the point estimate is
    not in one — there is nothing for it to be dropped into (§8.5).
    """
    df = near_separated_subgroup_frame()
    keys = tuple(f"{s}.{q}" for s in config.SUBGROUPS
                 for q in ("beta", "gamma", "beta_plus_gamma"))

    # door one — the replicate loop: dropped, counted, and never substituted
    replicate = sensitivity._subgroup_replicate(df, keys, data.WORKBOOK)
    assert replicate.failures[f"{SUBGROUP}.gamma"] == "separation"
    assert f"{SUBGROUP}.gamma" not in replicate.values

    # door two — the point estimate: `C.SchemaError` H8, naming the subgroup and the token
    _, ps, est, audit = fitted(df)
    with pytest.raises(config.SchemaError) as e:
        sensitivity.subgroups(df, ps, est, audit)
    assert "H8" in message_of(e) and "G9" in message_of(e) and SUBGROUP in message_of(e)
    # ...and it is NOT a FitError, which is the whole of the difference
    assert not isinstance(e.value, model.FitError)
