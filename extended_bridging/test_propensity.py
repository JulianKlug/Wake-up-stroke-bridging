"""Acceptance tests for Stage 6b — §12.8-12.14 of `specs/stage6_propensity_and_weights.md`.

This is the first test module that takes the **whole pipeline end to end with no `data/`**, which no
stage before it could do past Stage 5's P2. `test_cohort.cohort_frame()`'s five-record cohort looks
unfittable — 9 [§6] covariates over 5 patients — and works, because the constant-column rule strips
nine of the twelve candidate columns, leaving `core_ml`, `tmax6_ml`, `center_CHUV` and an intercept:
four parameters over five records, full rank, converging in 11 iterations. **It is also the best
available test of §4.3**, because nine dropped columns is a bigger signal than the workbook's one.

Fixtures are imported from the modules that declare them, never re-declared (Stage 3 §12's rule).
Tests needing the private workbook are marked `skipif(not DATA_XLSX.exists())` and take this file's
**own** module-scoped `workbook` fixture.

**Bare frame or normalised frame** — Stage 5 §12.0.1's distinction carries forward unchanged. Every
frame below is the normalised one: `run(...)`, then `derive`, then `classify`, then `build`. A test
meaning to reach `design` and passing a bare frame gets D3, not the check it intended, and the failure
looks like a passing `pytest.raises` whose message is never read.

This file is **not** exempt from the Stage 1 §7 raw-name scan and must not become exempt.
"""
from __future__ import annotations

import ast
import inspect
import subprocess
import sys
import textwrap
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

import cohort
import config
import data
import derive
import eligibility
import model
import propensity
from test_cohort import built, cohort_frame, classified, set_cell

DATA_GATED = pytest.mark.skipif(
    not config.DATA_XLSX.exists(),
    reason="the private workbook is gitignored and absent from this checkout")

MODULE = Path(propensity.__file__).resolve()
MODULE_DIR = MODULE.parent
SOURCE = MODULE.read_text(encoding="utf-8")

ARM_CODES = tuple(sorted(config.TREATMENT_LABELS))


def message_of(excinfo) -> str:
    return str(excinfo.value)


def rows_of(entry: data.AuditEntry) -> dict[str, tuple[str, ...]]:
    """A table's body keyed by its first cell, so a test names the row it is asserting."""
    return {row[0]: row[1:] for row in entry.table[1:]}


def fitted(df=None):
    """The whole pipeline through Stage 6, on the hand-built cohort. No `data/` needed."""
    frame, audit = built(df)
    return frame, audit, propensity.fit(frame, audit)


# --- 12.0.2  the golden vector -----------------------------------------------------------------------
#
# Pinned as literals from the SYNTHETIC frame, so that no patient-derived number enters git (§7.3): the
# fitted coefficients live in the audit log, which is gitignored, and the workbook fit is asserted on
# properties instead (§12.14).
#
# ONE TRAP, and it is why the assertion is by identifier. HAND-2's weight and COHORT-3's weight AND
# propensity are all 0.288470, because the fixture's two CHUV records happen to have linear predictors
# of equal magnitude and opposite sign. A test asserting weights by value alone would pass while
# confusing the two records, so every assertion below names its record and none compares against a
# sorted list of values.
#
# Asserted to 1e-6, which is looser than the agreement §16b measures against an external kernel and
# deliberately so: this is a regression pin against an edit to OUR algorithm, and pinning it at machine
# precision would make it fail on a numpy patch release rather than on a mistake.

GOLDEN = {
    #             e         w
    "HAND-1":   (0.929038, 0.070962),   # HUG,  bridging
    "HAND-2":   (0.711530, 0.288470),   # CHUV, bridging
    "HAND-5":   (0.702237, 0.297763),   # HUG,  bridging
    "COHORT-1": (0.257067, 0.257067),   # HUG,  control  — w == e in the control arm
    "COHORT-3": (0.288470, 0.288470),   # CHUV, control  — and this is the 0.288470 collision
}


@pytest.mark.parametrize("case_id", sorted(GOLDEN))
def test_the_golden_vector_by_identifier(case_id):
    df, _, ps = fitted()
    where = df["case_id"] == case_id
    expected_e, expected_w = GOLDEN[case_id]
    assert float(ps.e[where].iloc[0]) == pytest.approx(expected_e, abs=1e-6)
    assert float(ps.w[where].iloc[0]) == pytest.approx(expected_w, abs=1e-6)


def test_the_golden_fit_and_its_dropped_columns():
    df, _, ps = fitted()
    assert ps.fit.iterations == 11
    assert ps.fit.converged_on == "likelihood"
    assert list(ps.fit.columns) == ["core_ml", "tmax6_ml", "center_CHUV"]
    assert len(ps.dropped) == 9                     # nine of twelve candidates — §4.3's best signal
    assert ps.ess[1] == pytest.approx(2.4413, abs=1e-4)
    assert ps.ess[0] == pytest.approx(1.9934, abs=1e-4)


# --- 12.8  `fit` returns aligned Series with the missingness reimposed --------------------------------

def test_e_w_and_in_model_carry_the_cohorts_index_EQUAL_not_merely_the_same_length():
    df, _, ps = fitted()
    for series in (ps.e, ps.w, ps.in_model):
        assert series.index.equals(df.index)


def test_in_model_is_total_boolean_and_never_missing():
    _, _, ps = fitted()
    assert ps.in_model.dtype == bool
    assert not ps.in_model.isna().any()


def test_the_sentinel_is_np_nan_and_pd_NA_would_not_fit_in_the_column():
    # Measured: a float64 Series has exactly ONE missing value and pd.NA is not it. The masked Float64
    # dtype would take it and is declined, because every downstream .to_numpy(dtype=float) converts it
    # back to nan silently — §3.1's first hazard, on the one column most likely to be summed.
    _, _, ps = fitted()
    assert ps.e.dtype == "float64" and ps.w.dtype == "float64"
    with pytest.raises(TypeError):
        pd.Series(pd.NA, index=ps.e.index, dtype="float64")


def _blanked():
    """The built cohort with HAND-1's `core_ml` blanked, and a fresh Audit. §12.8's frame.

    Blanked on the BUILT frame rather than upstream, so that the cohort is the same five records and
    the only difference is the covariate — blanking before `build` would also change `penumbra_ml`,
    which Stage 2 recomputes from it, and the [§13] median the cohort freezes.
    """
    frame, audit = built()
    return set_cell(frame, "HAND-1", "core_ml", np.nan), audit


def test_blanking_a_covariate_excludes_exactly_that_record():
    frame, audit = _blanked()
    ps = propensity.fit(frame, audit)
    where = frame["case_id"] == "HAND-1"
    assert not bool(ps.in_model[where].iloc[0])
    assert bool(ps.e[where].isna().iloc[0]) and bool(ps.w[where].isna().iloc[0])
    assert int(ps.in_model.sum()) == len(frame) - 1


def test_the_excluded_record_LEFT_the_design_rather_than_entering_it_as_a_zero():
    """The real form of the property, and an earlier draft of §12.8 got it wrong.

    That draft asserted the four survivors match the golden vector to 1e-12, "which is what proves the
    excluded record left the design rather than entering it as a zero". **That is false, and measured:
    blanking HAND-1's `core_ml` moves the four survivors by up to 0.0478** — a logistic fit over four
    records is a different fit from one over five, so the change is the CORRECT behaviour and the
    assertion would have failed. Worse, the natural repair — loosening the tolerance until it passes —
    tests nothing at all.

    What actually distinguishes "left the design" from "entered as a zero" is that refitting the
    four-record frame DIRECTLY, with no blanking, gives identical probabilities. A zero row would
    change the fit, and this is what would catch it.
    """
    frame_b, audit_b = _blanked()
    ps_b = propensity.fit(frame_b, audit_b)

    frame_full, audit_d = built()
    frame_d = frame_full[frame_full["case_id"] != "HAND-1"]
    ps_d = propensity.fit(frame_d, audit_d)

    for case_id in ("HAND-2", "HAND-5", "COHORT-1", "COHORT-3"):
        a = float(ps_b.e[frame_b["case_id"] == case_id].iloc[0])
        b = float(ps_d.e[frame_d["case_id"] == case_id].iloc[0])
        assert a == pytest.approx(b, abs=1e-12)
    # and the survivors really did move relative to the five-record fit, which is why the
    # earlier draft's assertion could not have passed
    moved = max(abs(float(ps_b.e[frame_b["case_id"] == c].iloc[0]) - GOLDEN[c][0])
                for c in ("HAND-2", "HAND-5", "COHORT-1", "COHORT-3"))
    assert moved > 1e-3


def test_the_design_handed_to_firth_has_as_many_rows_as_in_model():
    frame, audit = _blanked()
    ps = propensity.fit(frame, audit)
    assert len(ps.fit.p) == int(ps.in_model.sum()) == len(frame) - 1


def test_F5_asserts_strict_interiority_rather_than_clipping():
    _, _, ps = fitted()
    inside = ps.e[ps.in_model]
    assert ((inside > 0.0) & (inside < 1.0)).all()
    # and the assertion is reachable: a Fit whose probabilities hit the boundary raises rather than
    # being clipped into range, which is what would turn a degenerate fit into a weight of 0.0.
    degenerate = model.Fit(np.zeros(2), np.array([1.0, 0.5]), 3, "likelihood", ("x",), 0.0, 0, 0)
    with pytest.raises(model.FitError) as e:
        propensity._assert_probabilities(degenerate, pd.Series(["A", "B"]))
    assert "F5" in message_of(e) and "not strictly in (0, 1)" in message_of(e)


# --- 12.9  weights and ESS ------------------------------------------------------------------------------

@pytest.mark.parametrize("code", ARM_CODES)
def test_the_weight_is_one_minus_e_treated_and_e_control(code):
    # Asserted through TREATMENT_LABELS' codes rather than literal 0/1.
    df, _, ps = fitted()
    arm = ps.in_model & (df[config.TREATMENT] == code)
    expected = (1.0 - ps.e[arm]) if code == 1 else ps.e[arm]
    assert np.allclose(ps.w[arm].to_numpy(float), expected.to_numpy(float), atol=0, rtol=0)


def test_every_weight_is_in_the_unit_interval_and_nothing_is_trimmed():
    _, _, ps = fitted()
    inside = ps.w[ps.in_model]
    assert ((inside >= 0.0) & (inside <= 1.0)).all()


def test_ess_is_hand_computable():
    # §6.2's two hand computations. A weighted sample of four with one weight of three is worth three
    # observations: 6² / 12 = 3.
    assert propensity.ess([1] * 10) == 10.0
    assert propensity.ess([1, 1, 1, 3]) == 3.0


def test_ess_raises_on_an_all_zero_arm_and_the_pilots_zero_is_declined():
    with pytest.raises(model.FitError) as e:
        propensity.ess([0.0, 0.0, 0.0])
    assert "every weight in this arm is zero" in message_of(e)
    # the companion: `pilots/analysis.py:140` returns 0.0 here, and 0.0 is a number every downstream
    # ratio accepts — a structural non-positivity rendered as a finite effective sample size.
    assert 10.0 / (0.0 + 1e-12) > 0        # a plausible finite ratio, which is the whole problem


def test_ess_on_an_EMPTY_arm_raises_a_DIFFERENT_message():
    # They share a branch by arithmetic — measured, np.sum(np.asarray([]) ** 2) is 0.0 — and they are
    # different diagnoses: all-zero means the fit put every patient on the boundary; empty means the
    # mask lost an arm F2 had already established was present. Stage 7 calls `ess` per centre per arm
    # [§9], where the empty case is the one that actually occurs.
    assert float(np.sum(np.asarray([]) ** 2)) == 0.0
    with pytest.raises(model.FitError) as empty:
        propensity.ess([])
    with pytest.raises(model.FitError) as zeros:
        propensity.ess([0.0, 0.0])
    assert "no weighted patient at all" in message_of(empty)
    assert "every weight in this arm is zero" in message_of(zeros)
    assert message_of(empty) != message_of(zeros)


def test_ess_drops_missing_before_summing():
    # So passing the full cohort-length `w` and the masked `w` give the same number. A Kish sum over a
    # nan is nan, and the tempting "fix" for that is the fillna(0) this module forbids.
    df, _, ps = fitted()
    arm = df[config.TREATMENT] == 1
    assert propensity.ess(ps.w[arm]) == propensity.ess(ps.w[arm & ps.in_model])


# --- 12.10  the weights do not balance exactly under Firth, and do under an MLE score -------------------

def test_the_arm_weight_sums_are_NOT_equal_under_firth():
    """Asserted as an INEQUALITY rather than tolerated (§6.3).

    Overlap weights have an exact-balance property under an UNPENALISED logistic MLE, because that is
    the score equation. Firth's modified score is a different equation, so the property holds only
    approximately. A test asserting equality to a loose tolerance would pass under Firth today and pass
    under a silently reinstated MLE tomorrow — which is the substitution [§7] forbids.

    This half needs no `data/`, so the discrimination §6.3 depends on survives a plain checkout even
    though its MLE companion does not.
    """
    df, _, ps = fitted()
    sums = [float(ps.w[ps.in_model & (df[config.TREATMENT] == code)].sum()) for code in ARM_CODES]
    assert abs(sums[0] - sums[1]) > 1e-6


@DATA_GATED
def test_the_arm_weight_sums_ARE_equal_under_an_mle_score(workbook):
    """The companion, DATA-GATED onto the workbook design — and an earlier draft ran it on the hand
    cohort, where it cannot hold.

    `cohort_frame()`'s cohort is 5 records over 3 columns and COMPLETELY SEPARATED, so there is no
    maximum-likelihood estimate to compare against: measured, `statsmodels.Logit` Newton reports
    converged=False with max|coef| = 130 and a weight-sum difference of −4.4e-11, which "passes" a
    1e-10 assertion only because both sums have underflowed, while BFGS converges to a different point
    with max|coef| = 168 and a difference of 3.9e-07, which fails it. The assertion's verdict was a
    function of the optimiser. On the workbook design the MLE converges and the property holds
    properly.
    """
    import statsmodels.api as sm
    df, _ = workbook
    mask = model.complete_cases(df, config.PS_COVARIATES)
    X, _ = model.design(df.loc[mask], config.PS_COVARIATES)
    a = df.loc[mask, config.TREATMENT].to_numpy(dtype=float)
    mle = sm.Logit(a, sm.add_constant(X.to_numpy(dtype=float))).fit(disp=0)
    e = np.asarray(mle.predict(), dtype=float)
    w = np.where(a == 1.0, 1.0 - e, e)
    assert abs(float(w[a == 1].sum()) - float(w[a == 0].sum())) < 1e-8


# --- 12.11  the four audit entries, and the log ----------------------------------------------------------

STEPS = ("covariate_completeness", "design_matrix", "propensity_fit", "overlap_weights")


def test_the_four_entries_appear_in_order_at_positions_captured_BEFORE_the_call():
    # Never a tail slice of the form [-n:], which is the repair Stage 5 §10 landed in test_derive.py
    # and which this stage must not reintroduce — it inserts four entries where it would bite.
    frame, audit = built()
    before = len(audit.entries)
    propensity.fit(frame, audit)
    added = audit.entries[before:]
    assert tuple(e.step for e in added) == STEPS
    assert {e.kind for e in added} == {"model"}


def test_the_model_entries_render_under_fitted_models_between_cohort_and_structural():
    _, audit, _ = fitted()
    rendered = audit.to_markdown()
    assert "## Fitted models" in rendered
    order = [rendered.index(f"## {data._HEADINGS[k]}")
             for k in ("cohort", "model", "structural")]
    assert order == sorted(order)
    for step in STEPS:
        assert f"- **{step}**" in rendered


def test_data_kinds_is_the_declared_nine_with_model_in_position():
    assert data.KINDS == (
        "provenance", "contract", "correction", "observation",
        "derivation", "cohort", "model", "structural", "missingness")
    assert data._HEADINGS["model"] == "Fitted models"


def test_covariate_completeness_names_as_many_identifiers_as_its_n():
    frame, audit = _blanked()
    propensity.fit(frame, audit)
    entry = audit.entry("model", "covariate_completeness")
    assert entry.n == 1
    assert entry.case_ids == ("HAND-1",)


def test_the_exclusion_check_FIRES_on_a_duplicated_case_id_and_records_NOTHING():
    """Modelled on test_cohort.py's working precedent, and the second assertion is what the
    raise-before-record ordering buys: `Audit` has no removal path, so recording first would leave the
    very inconsistency the raise exists to prevent.
    """
    frame, audit = built()
    frame = frame.copy()
    frame.loc[frame.index[:2], "case_id"] = "DUPLICATE"
    in_model = pd.Series(False, index=frame.index)
    before = len(audit.entries)
    with pytest.raises(config.SchemaError) as e:
        propensity._record_exclusion(frame, in_model, audit)
    assert "can name" in message_of(e)
    assert len(audit.entries) == before        # nothing was recorded


def test_a_plain_sorted_list_could_NEVER_fail_that_check():
    """The companion: the necessity of the `set` and the `notna` is MEASURED, not commented.

    `Audit.record` normalises with `tuple(sorted(str(c) for c in case_ids))` (data.py:255) — it sorts
    WITHOUT deduplicating, and `str(pd.NA)` is the four characters `<NA>`. So a `_record_exclusion`
    collecting a plain `sorted(...)` list has `len(case_ids) == n` by construction, and §12.11's
    acceptance test could not have passed.
    """
    excluded = pd.Series(["DUPLICATE", "DUPLICATE", pd.NA], dtype="object")
    plain = sorted(str(c) for c in excluded)
    assert len(plain) == 3                                     # the check is dead code
    guarded = tuple(sorted(set(excluded[excluded.notna()])))
    assert len(guarded) == 1                                   # and alive with both mechanisms


@DATA_GATED
def test_completeness_agrees_with_absence_by_cohort_column_for_every_covariate(workbook):
    """What replaces the per-centre columns §7.4 removed: the numbers stay in ONE place, and a drift
    between the two tables is a red suite."""
    _, audit = workbook
    ours = rows_of(audit.entry("model", "covariate_completeness"))
    theirs = rows_of(audit.entry("missingness", "absence_by_cohort_column"))
    for covariate in config.PS_COVARIATES:
        assert ours[covariate][0] == theirs[covariate][2]      # n_absent, read from both tables
    excluded_by = {c for c in config.PS_COVARIATES if ours[c][1] == "yes"}
    assert excluded_by == {"core_ml", "tmax6_ml"}


def test_the_design_matrix_table_has_one_row_per_LEVEL_and_marks_the_references():
    _, audit, ps = fitted()
    entry = audit.entry("model", "design_matrix")
    rows = rows_of(entry)
    # one row per linear covariate plus one per declared level of every factor — NOT one per column
    expected = len([c for c in config.PS_COVARIATES if c not in config.CATEGORICAL]) + sum(
        len(config.FACTOR_LEVELS[c]) for c in config.PS_COVARIATES if c in config.CATEGORICAL)
    assert len(rows) == expected
    assert rows["center_HUG"][-1] == "reference (baseline)"
    assert rows["onset_type_witnessed"][-1] == "reference (baseline)"
    assert rows["center_USZ"][-1] == "dropped: constant"
    # n is the PARAMETER count including the intercept, and is NOT the row count. The two being
    # different is the point: a level that contributes no parameter still gets a row.
    assert entry.n == len(ps.fit.columns) + 1
    assert entry.n != len(rows)


def test_the_reference_rows_are_in_NEITHER_the_matrix_nor_dropped():
    # Which is why `_design_table` ranges over the declaration: a table built from `X.columns +
    # dropped` would name every column of the model and none of its baselines.
    df, _, ps = fitted()
    X, dropped = model.design(df.loc[ps.in_model], config.PS_COVARIATES)
    for reference in ("center_HUG", "onset_type_witnessed"):
        assert reference not in X.columns and reference not in dropped


def test_the_coefficient_table_is_in_Fit_columns_order():
    _, audit, ps = fitted()
    rows = [row[0] for row in audit.entry("model", "propensity_fit").table[1:]]
    assert rows == ["intercept", *ps.fit.columns]


def test_the_overlap_weights_table_renders_a_row_per_declared_arm_plus_all():
    _, audit, _ = fitted()
    rows = rows_of(audit.entry("model", "overlap_weights"))
    for label in config.TREATMENT_LABELS.values():
        assert label in rows
    assert "all (pooled)" in rows


def test_the_overlap_weights_detail_carries_the_conditionality_statement():
    _, audit, _ = fitted()
    detail = audit.entry("model", "overlap_weights").detail
    assert "CHANGING THE PROPENSITY SPECIFICATION CHANGES THE POPULATION" in detail
    assert "h(X) = e(X){1 − e(X)}" in detail


def test_the_weighted_population_table_reconciles_its_unweighted_column_with_the_frame():
    """§7.6, and the unweighted column IS the reconciliation: it is asserted against the cohort's own
    means computed directly, so the description of the weighted population cannot drift from the
    population it describes.
    """
    df, audit, ps = fitted()
    rows = rows_of(audit.entry("model", "overlap_weights"))
    sub = df.loc[ps.in_model]
    for code, label in config.TREATMENT_LABELS.items():
        column = list(config.TREATMENT_LABELS.values()).index(label)
        arm = sub[sub[config.TREATMENT] == code]
        assert rows["core_ml"][column] == data._fmt(arm["core_ml"].mean())
        assert rows["center = CHUV"][column] == data._fmt((arm["center"] == "CHUV").mean())


def test_a_factor_contributes_one_row_per_LEVEL_to_the_weighted_population_table():
    _, audit, _ = fitted()
    rows = rows_of(audit.entry("model", "overlap_weights"))
    for level in config.FACTOR_LEVELS["center"]:
        assert f"center = {level}" in rows
    for level in config.FACTOR_LEVELS["onset_type"]:
        assert f"onset_type = {level}" in rows


def test_the_ATO_weighted_column_DIFFERS_from_the_unweighted_one():
    # A weighted population identical to the unweighted one would mean the weights did nothing and the
    # table would be decoration.
    _, audit, _ = fitted()
    rows = rows_of(audit.entry("model", "overlap_weights"))
    assert any(rows[key][:2] != rows[key][2:4]
               for key in rows if key.startswith(("core_ml", "tmax6_ml", "center = ")))


def test_no_entry_contains_a_numpy_repr():
    # Every float reaches the log through data._fmt, which is what makes byte identity hold.
    _, audit, _ = fitted()
    for entry in (e for e in audit.entries if e.kind == "model"):
        cells = [c for row in (entry.table or ()) for c in row]
        assert not [c for c in cells if "np." in c or "array(" in c or "float64" in c]


def test_two_runs_render_identical_markdown():
    assert fitted()[1].to_markdown() == fitted()[1].to_markdown()


def test_the_log_is_identical_across_interpreters_with_different_hash_seeds(tmp_path):
    """The two-seed driver of Stage 5 §12.9, extended one stage.

    The driver is written out rather than sketched, for Stage 3 §12.12's reason: an implementer
    choosing it freely can choose one that fits nothing at all and still see two identical outputs.
    """
    driver = textwrap.dedent("""
        import sys
        sys.path.insert(0, %r)
        import data, propensity
        from test_propensity import fitted
        _, audit, _ = fitted()
        sys.stdout.write(audit.to_markdown())
    """) % str(MODULE_DIR)
    script = tmp_path / "driver.py"
    script.write_text(driver, encoding="utf-8")
    outputs = []
    for seed in ("0", "1"):
        done = subprocess.run([sys.executable, str(script)], capture_output=True, text=True,
                              env={"PYTHONHASHSEED": seed, "PATH": "/usr/bin:/bin"},
                              cwd=str(MODULE_DIR))
        assert done.returncode == 0, done.stderr
        outputs.append(done.stdout)
    assert outputs[0] == outputs[1]
    assert "## Fitted models" in outputs[0]


# --- 12.12  the module boundary holds ---------------------------------------------------------------------

def test_fit_takes_exactly_df_and_audit():
    # No covariate parameter, no `scheme`, no centre list (§6.4, §16). Asserted by `inspect`, so adding
    # a defaulted keyword fails the test rather than passing silently.
    assert list(inspect.signature(propensity.fit).parameters) == ["df", "audit"]
    assert all(p.default is inspect.Parameter.empty
               for p in inspect.signature(propensity.fit).parameters.values())


def test_ess_is_public_and_lives_here_rather_than_in_model():
    assert callable(propensity.ess)
    assert not hasattr(model, "ess")


# --- 12.13  Stage 6 adds nothing to the frame and needs Stage 5 --------------------------------------------

def test_the_frame_comes_back_unchanged():
    frame, audit = built()
    before = frame.copy(deep=True)
    propensity.fit(frame, audit)
    pd.testing.assert_frame_equal(frame, before)
    assert list(frame.columns) == list(before.columns)      # no column added (§0.2)


def test_fit_on_the_UNRESTRICTED_frame_is_legitimate_and_keeps_the_ineligible():
    """Not asserted to raise: it is a legitimate call, and Stage 12 makes something like it.

    What is asserted is the property Stage 12 needs and roadmap INVARIANT 2 forbids for [§7]: the
    ineligible patients are IN this fit. The guarantee that they are not in the [§7] population is
    Stage 5's, not this stage's.
    """
    frame, audit = classified()
    ps = propensity.fit(frame, audit)
    ineligible = frame[config.ELIGIBILITY] == config.INELIGIBLE
    assert int(ineligible.sum()) >= 1
    assert bool(ps.in_model[ineligible].all())              # in the fit, deliberately
    assert len(frame) == 9                                  # the hand frame, not the workbook


def test_fit_raises_F2_on_a_single_armed_frame():
    frame, audit = built()
    one_arm = frame[frame[config.TREATMENT] == 1]
    with pytest.raises(config.SchemaError) as e:
        propensity.fit(one_arm, audit)
    assert "F2" in message_of(e) and "both arms" in message_of(e)


def test_fit_raises_F1_on_a_frame_whose_exposure_is_missing():
    frame, audit = built()
    frame = set_cell(frame, "HAND-1", config.TREATMENT, pd.NA)
    with pytest.raises(config.SchemaError) as e:
        propensity.fit(frame, audit)
    assert "F1" in message_of(e) and "HAND-1" in message_of(e)


def test_fit_works_on_a_frame_with_no_subgroup_column():
    # Stage 6 never reads a subgroup.
    frame, audit = built()
    frame = frame.drop(columns=list(config.COHORT_DEPENDENT_SUBGROUPS))
    assert propensity.fit(frame, audit).fit.iterations >= 1


# --- 12.14  structural facts from the workbook [data-gated] ------------------------------------------------

@pytest.fixture(scope="module")
def workbook():
    """This module's OWN fixture, as Stage 5's test file has its own. Never shared across modules."""
    df, audit = data.load(data.WORKBOOK)
    df = derive.derive(df, audit)
    df = eligibility.classify(df, audit)
    df = cohort.build(df, audit)
    propensity.fit(df, audit)
    return df, audit


@pytest.fixture(scope="module")
def workbook_ps(workbook):
    df, audit = workbook
    return df, propensity.fit(df, data.Audit(data.WORKBOOK))


@DATA_GATED
def test_the_cohort_is_93_and_in_model_is_92(workbook_ps):
    df, ps = workbook_ps
    assert len(df) == 93
    assert int(ps.in_model.sum()) == 92
    excluded = df.loc[~ps.in_model]
    assert len(excluded) == 1
    assert excluded["center"].iloc[0] == "Lugano"
    assert int(excluded[config.TREATMENT].iloc[0]) == 0            # control arm
    assert bool(excluded[["core_ml", "tmax6_ml"]].isna().all(axis=None))


@DATA_GATED
def test_the_design_is_eleven_columns_and_3_25_treated_per_parameter(workbook_ps):
    df, ps = workbook_ps
    assert ps.dropped == ("center_USZ",)
    assert len(ps.fit.columns) == 11
    parameters = len(ps.fit.columns) + 1
    assert parameters == 12
    treated = int((df.loc[ps.in_model, config.TREATMENT] == 1).sum())
    assert treated == 39
    # Asserted as a NUMBER, so a covariate added without a [§13] amendment moves it.
    assert treated / parameters == pytest.approx(3.25, abs=0.005)


@DATA_GATED
def test_the_workbook_fit_converges_in_seven_iterations(workbook_ps):
    _, ps = workbook_ps
    assert ps.fit.iterations == 7
    assert ps.fit.converged_on == "likelihood"
    inside = ps.e[ps.in_model]
    assert float(inside.min()) > 0.02 and float(inside.max()) < 0.94


@DATA_GATED
def test_the_workbook_ess_is_30_34_treated_and_29_20_control(workbook_ps):
    df, ps = workbook_ps
    assert ps.ess[1] == pytest.approx(30.34, abs=0.005)
    assert ps.ess[0] == pytest.approx(29.20, abs=0.005)
    # each against its own denominator, which is why [§11] requires it beside the estimate: the
    # control denominator is 53 rather than 54 because of the one excluded record
    assert int((df.loc[ps.in_model, config.TREATMENT] == 1).sum()) == 39
    assert int((df.loc[ps.in_model, config.TREATMENT] == 0).sum()) == 53


@DATA_GATED
def test_the_largest_workbook_weight_is_below_one_and_no_weight_is_zero(workbook_ps):
    _, ps = workbook_ps
    inside = ps.w[ps.in_model]
    assert 0.9 < float(inside.max()) < 1.0            # measured 0.965
    assert float(inside.min()) > 0.0


@DATA_GATED
def test_constant_covariates_is_empty_while_dropped_is_not(workbook_ps):
    # §10's detect-versus-drop distinction, asserted rather than described: `constant_covariates`
    # DETECTS and logs at Stage 5, `design` DROPS at Stage 6, and two stages cannot both own it.
    df, ps = workbook_ps
    assert derive.constant_covariates(df, config.PS_COVARIATES) == ()
    assert ps.dropped != ()


@DATA_GATED
def test_the_workbook_audit_ledger_is_seven_four_one_six_four(workbook):
    _, audit = workbook
    counts = {kind: len([e for e in audit.entries if e.kind == kind]) for kind in data.KINDS}
    assert counts["model"] == 4
    assert len(audit.entries) == 22        # load 7 / derive 4 / classify 1 / build 6 / fit 4


# --- the entry-position rule, verified by grep (DoD-13) -----------------------------------------------------

@pytest.mark.parametrize("name", ["test_model.py", "test_propensity.py"])
def test_no_tail_slice_of_the_entries_appears_in_either_new_test_file(name):
    """Stage 5 §10's repair, and this stage inserts four entries where it would bite.

    An assertion on `audit.entries[-4:]` passes whether or not the entries it names are the ones the
    call added, because any four trailing entries satisfy it.
    """
    source = (MODULE_DIR / name).read_text(encoding="utf-8")
    offenders = [
        node.lineno for node in ast.walk(ast.parse(source))
        if isinstance(node, ast.Subscript) and isinstance(node.slice, ast.Slice)
        and isinstance(node.slice.lower, ast.UnaryOp)
        and isinstance(node.slice.lower.op, ast.USub)]
    assert offenders == []
