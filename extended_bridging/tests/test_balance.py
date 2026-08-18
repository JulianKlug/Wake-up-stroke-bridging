"""Acceptance tests for Stage 7 — §12 of `specs/stage7_balance_and_overlap.md`.

Three kinds of input, and the split follows Stage 6 §12.0's. **Hand-built vectors** carry `smd`'s
arithmetic and all four of §5.3's undefined branches, and need no frame. **`cohort_frame()`'s
cohort**, taken through `derive → classify → build → fit`, exercises more of this stage than the
workbook does: 14 of its 19 rows read exactly 0.0, `nihss_baseline` is the row §5.3a is about, every
centre has exactly one control so every within-centre SMD is undefined while the ESS and propensity
range are not, and two of its four declared centres carry no cohort record at all. **The workbook**
is data-gated and reaches §12.12 only.

Fixtures are imported from the modules that declare them, never re-declared (Stage 3 §12's rule) —
`cohort_frame` from `test_cohort.py` and `hand_source`/`run` from `test_data.py`, which are two
modules where the natural guess is one.

**Bare frame or normalised frame** — Stage 5 §12.0.1's distinction carries forward unchanged.
`cohort_frame()` returns raw centre codes and `run(...)` maps them to `CENTER_ORDER`'s labels; every
frame below is the normalised one. A test meaning to reach `assess` and passing a bare frame never
gets there, because Stage 5's C4 raises first, and the failure looks like a passing `pytest.raises`
whose message is never read.

This file is **not** exempt from the Stage 1 §7 raw-name scan and must not become exempt.
"""
from __future__ import annotations

import ast
import inspect
import subprocess
import sys
import textwrap
import types
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

import balance
import cohort
import config
import data
import derive
import eligibility
import model
import propensity
from test_cohort import built, cohort_frame, set_cell

DATA_GATED = pytest.mark.skipif(
    not config.DATA_XLSX.exists(),
    reason="the private workbook is gitignored and absent from this checkout")

MODULE = Path(balance.__file__).resolve()
MODULE_DIR = MODULE.parent                       # the flat module root
TESTS_DIR = Path(__file__).resolve().parent      # this file's own directory
SOURCE = MODULE.read_text(encoding="utf-8")

ARM_CODES = tuple(sorted(config.TREATMENT_LABELS))


def message_of(excinfo) -> str:
    return str(excinfo.value)


def rows_of(entry: data.AuditEntry) -> dict[str, tuple[str, ...]]:
    """A table's body keyed by its first cell, so a test names the row it is asserting."""
    return {row[0]: row[1:] for row in entry.table[1:]}


def assessed(df: pd.DataFrame | None = None):
    """The whole pipeline through Stage 7, on the hand-built cohort. No `data/` needed."""
    frame, audit = built(df)
    ps = propensity.fit(frame, audit)
    return frame, ps, audit, balance.assess(frame, ps, audit)


def fitted(df: pd.DataFrame | None = None):
    """Through Stage 6 only, for the tests that drive `assess` themselves."""
    frame, audit = built(df)
    return frame, propensity.fit(frame, audit), audit


# --- 12.0  the hand-built vectors ------------------------------------------------------------------

def hand_smd():
    """§12.3 — the roadmap's hand-computed value. sd = sqrt((0.5 + 0.5)/2), (3.5 - 1.5)/sd."""
    return (np.array([1.0, 2.0, 3.0, 4.0]), np.array([0.0, 0.0, 1.0, 1.0]), np.ones(4))


def separating():
    """§12.4a — constant within each arm, different between: the branch pilots returns 0.0 for."""
    return (np.array([0.0, 0.0, 1.0, 1.0]), np.array([0.0, 0.0, 1.0, 1.0]), np.ones(4))


def constant_everywhere():
    """§12.4a — the SAME constant in both arms, and NOT an indicator, so §3.1's 1-ulp fact bites.

    14.0 rather than 1.0 deliberately: np.average of a vector of exact 1.0s returns exactly 1.0
    under any weights, so an indicator cannot exercise §5.3a and 17 of 19 fixture rows are
    indicators.

    **The control arm's weights are (0.2, 0.8) and that is not free either.** What makes the two
    weighted means differ is whether an arm's Σ(w·x)/Σw division ROUNDS, which is a property of the
    weights and not of the arm size — `test_only_SOME_weight_vectors_round` is that measurement. An
    earlier draft of this fixture used (0.1, 0.9) against (0.3, 0.7), two vectors that are both
    exact, so both means were 14.0 to every bit and the rejected-form companion below returned 0.0
    rather than nan: the fixture passed while exercising nothing.
    """
    return (np.full(4, 14.0), np.array([0.0, 0.0, 1.0, 1.0]),
            np.array([0.2, 0.8, 0.3, 0.7]))


def one_armed():
    """§12.4 — a covariate observed in only one arm, after the caller's own mask. Branch 1 and 3."""
    return (np.array([1.0, 2.0, 3.0]), np.array([0.0, 1.0, 1.0]), np.ones(3))


def pilots_smd(x, d, w) -> float:
    """`pilots/analysis.py:194-208`, verbatim in form, as the companion for §12.4 and §12.4a.

    Written out here rather than imported: `pilots/` is a separate project with its own config, and
    what is being compared is the *form* of the estimator — its guards and its fall-through — not a
    call into a package this repository does not depend on.
    """
    x = np.asarray(x, dtype=float)
    t, c = d == 1, d == 0
    mt = np.average(x[t], weights=w[t]) if w[t].sum() else np.nan
    mc = np.average(x[c], weights=w[c]) if w[c].sum() else np.nan
    sd = np.sqrt((x[t].var(ddof=1) + x[c].var(ddof=1)) / 2)
    if not np.isfinite(sd):
        return np.nan
    return float((mt - mc) / sd) if sd > 0 else 0.0


# --- 12.3  `smd` reproduces hand-computed values [roadmap] -----------------------------------------

def test_smd_reproduces_the_roadmaps_hand_computed_value():
    # sd = sqrt((0.5 + 0.5)/2) = 0.70711, and (3.5 - 1.5)/sd = 2*sqrt(2). Computed rather than
    # pinned, so the criterion is the algebra and not a transcription of it.
    x, a, w = hand_smd()
    assert balance.smd(x, a, w) == pytest.approx(2 * np.sqrt(2), abs=1e-12)


def test_the_numerator_is_weighted_and_the_hand_value_says_so():
    # Non-unit weights, so the numerator's weighting is exercised independently of the denominator.
    # Treated {3, 4} weighted 1:3 has mean 3.75; control {1, 2} weighted 3:1 has mean 1.25; the
    # pooled SD is unchanged at sqrt(0.5), so the ratio is 2.5/0.70711.
    x, a, _ = hand_smd()
    w = np.array([3.0, 1.0, 1.0, 3.0])
    assert balance.smd(x, a, w) == pytest.approx(2.5 / np.sqrt(0.5), abs=1e-12)


def test_the_second_hand_value_with_weights_that_reverse_the_sign():
    # The same vector with the weights swapped between arms: treated mean 3.25, control mean 1.75,
    # so the numerator halves and the sign is unchanged. Two hand values rather than one, because a
    # single one is satisfied by an estimator that ignores `w` in exactly the right place.
    x, a, _ = hand_smd()
    w = np.array([1.0, 3.0, 3.0, 1.0])
    assert balance.smd(x, a, w) == pytest.approx(1.5 / np.sqrt(0.5), abs=1e-12)


def test_the_yardstick_does_not_move_when_the_weights_do():
    """[§9]'s sentence, asserted on the number itself and not on a consequence of it (§3, §5.5)."""
    x, a, _ = hand_smd()
    yardstick = balance._pooled_sd(x, a)
    for w in (np.ones(4), np.array([3.0, 1.0, 1.0, 3.0]), np.array([0.1, 0.9, 0.9, 0.1])):
        assert balance._pooled_sd(x, a) == yardstick
        # and the ratio really is that number's reciprocal times the weighted mean difference
        t, c = a == 1.0, a == 0.0
        numerator = np.average(x[t], weights=w[t]) - np.average(x[c], weights=w[c])
        assert balance._ratio(x, a, w, yardstick) == pytest.approx(numerator / yardstick, abs=1e-15)


def test_the_pooled_sd_uses_ddof_one():
    # §16b's oracle comparison turns on this: a ddof = 0 convention differs by sqrt(n/(n-1)), which
    # at n = 92 is 0.5% — visible at six significant figures and invisible to the eye.
    x, a, _ = hand_smd()
    assert balance._pooled_sd(x, a) == pytest.approx(np.sqrt((0.5 + 0.5) / 2), abs=1e-15)
    ddof_zero = np.sqrt((x[a == 1.0].var(ddof=0) + x[a == 0.0].var(ddof=0)) / 2)
    assert balance._pooled_sd(x, a) != pytest.approx(ddof_zero, abs=1e-6)


def test_smd_is_exactly_the_composition_of_the_two_privates():
    # §5.5's split cannot drift: the public function is `_ratio` over `_pooled_sd` and nothing else.
    for x, a, w in (hand_smd(), constant_everywhere(),
                    (np.array([1.0, 5.0, 2.0, 9.0]), np.array([0.0, 0.0, 1.0, 1.0]),
                     np.array([0.2, 0.8, 0.5, 0.5]))):
        expected = balance._ratio(x, a, w, balance._pooled_sd(x, a))
        got = balance.smd(x, a, w)
        assert (got == expected) or (np.isnan(got) and np.isnan(expected))


def test_smd_accepts_series_as_well_as_arrays():
    # `smd` is public and Stage 12 calls it [§14a]; its docstring promises either.
    x, a, w = hand_smd()
    assert balance.smd(pd.Series(x), pd.Series(a), pd.Series(w)) == pytest.approx(2 * np.sqrt(2))


# --- 12.4  the four undefined branches, each distinguished from the others -------------------------

UNDEFINED = {
    "branch1_arm_of_one": (np.array([1.0, 2.0, 3.0]), np.array([0.0, 0.0, 1.0]), np.ones(3)),
    "branch2_zero_total_weight": (np.array([1.0, 2.0, 3.0, 4.0]), np.array([0.0, 0.0, 1.0, 1.0]),
                                  np.array([1.0, 1.0, 0.0, 0.0])),
    "branch3_observed_in_one_arm_only": one_armed(),
    "branch4_separates_the_arms": separating(),
}


@pytest.mark.parametrize("name", sorted(UNDEFINED))
def test_each_undefined_branch_returns_missing_and_never_zero(name):
    x, a, w = UNDEFINED[name]
    assert np.isnan(balance.smd(x, a, w))


def test_one_armed_returns_nan_and_not_zero_and_the_pilots_agree_here():
    """The roadmap's second acceptance criterion, and the branch the pilots already get right.

    Honest about which of the two criteria `pilots/analysis.py` met: this one it did, through the
    same `if not np.isfinite(sd)` guard §16 lifts. §12.4a is the one it did not.
    """
    x, a, w = one_armed()
    assert np.isnan(balance.smd(x, a, w))
    assert np.isnan(pilots_smd(x, a, w))


def test_the_zero_weight_guard_precedes_np_average_which_RAISES(recwarn):
    """§3.1's first fact, measured rather than asserted in a comment.

    `np.average` raises `ZeroDivisionError` on weights summing to zero rather than returning `nan`,
    so an arm whose weights all vanish must be checked BEFORE the mean. A `ZeroDivisionError`
    escaping `assess` would stop a diagnostic mid-table, and [§9] wants the table.
    """
    x, a, w = UNDEFINED["branch2_zero_total_weight"]
    assert np.isnan(balance.smd(x, a, w))

    t = a == 1.0
    with pytest.raises(ZeroDivisionError):                    # the guard removed
        np.average(x[t], weights=w[t])


def test_the_arm_size_branch_is_about_the_ARM_and_not_about_the_data():
    # An arm of one is undefined; the same data with one record added to that arm is a number.
    x, a, w = np.array([1.0, 2.0, 3.0]), np.array([0.0, 0.0, 1.0]), np.ones(3)
    assert np.isnan(balance.smd(x, a, w))
    assert np.isfinite(balance.smd(np.append(x, 5.0), np.append(a, 1.0), np.ones(4)))


def test_a_nan_in_x_or_in_w_returns_nan_by_two_DIFFERENT_routes():
    """§5.3's two recorded routes. Neither is a branch and neither returns a number.

    `smd` does not mask — its docstring says the caller has already masked all three vectors — and
    `_table` does exactly that per row. This is what a caller who forgets gets: the undefined answer
    rather than a plausible one, which is the property [§9] needs.
    """
    x, a, w = hand_smd()
    assert np.isnan(balance.smd(np.array([1.0, np.nan, 3.0, 4.0]), a, w))   # via branch 3
    assert np.isnan(balance.smd(x, a, np.array([1.0, np.nan, 1.0, 1.0])))   # by propagation
    # and the second really does NOT fire branch 2 — `not nan` is False, so np.average propagates
    assert not (not np.array([1.0, np.nan])[np.array([True, True])].sum())


# --- 12.4a  the zero-variance split, both halves ---------------------------------------------------

def test_a_covariate_that_perfectly_separates_the_arms_is_undefined_here_and_zero_in_the_pilots():
    """The strongest imbalance expressible, reported by the pilots as the best-balanced row.

    `pilots/analysis.py`'s only guard is `sd > 0` and its fall-through is `0.0`, so a covariate on
    which every treated patient scores 1 and every control 0 passes `abs(smd) < SMD_THRESHOLD` and
    appears in the balance table as perfectly balanced. §5.3 branch 4 returns `nan` instead.
    """
    x, a, w = separating()
    assert np.isnan(balance.smd(x, a, w))
    assert pilots_smd(x, a, w) == 0.0


def rejected_zero_variance_form(x, a, w) -> float:
    """§5.3a's rejected branch 4: compare the two WEIGHTED means for equality.

    The natural way to write it, and wrong. The shipped form compares the raw arm constants, which
    is exact by construction because both arm variances are zero whenever the pooled SD is.
    """
    t, c = a == 1.0, a == 0.0
    return 0.0 if np.average(x[t], weights=w[t]) == np.average(x[c], weights=w[c]) else np.nan


def test_the_same_constant_in_both_arms_reads_zero():
    x, a, w = constant_everywhere()
    assert balance.smd(x, a, w) == 0.0
    assert np.isfinite(balance.smd(x, a, w))


def test_the_rejected_weighted_mean_form_calls_a_constant_covariate_UNDEFINED():
    """§5.3a's second half. Both halves are required — the first alone passes for the wrong
    estimator and the second alone passes for one that never returns 0.0 at all.

    The rejected form returns `nan` for a covariate that is 14.0 on every record of both arms: no
    imbalance at all, reported as undefined, rendering `missing` in the log and landing the row in
    `worst()`'s undefined list.
    """
    x, a, w = constant_everywhere()
    assert balance.smd(x, a, w) == 0.0
    assert np.isnan(rejected_zero_variance_form(x, a, w))


def test_the_rejected_form_on_the_REAL_fixture_cohort():
    """Where the defect was actually found (§18): `nihss_baseline`, 14.0 on every record.

    The real arms and the real [§7] weights, so this is the measurement rather than a reconstruction
    of it. `cohort_frame()`'s cohort is 3 bridging and 2 control, and the two weighted means come
    back as 14.000000000000002 and 14.0 — a difference of 1.776e-15.
    """
    df, ps, _ = fitted()
    x = df.loc[ps.in_model, "nihss_baseline"].to_numpy(dtype=float)
    a = df.loc[ps.in_model, config.TREATMENT].to_numpy(dtype=float)
    w = ps.w.loc[ps.in_model].to_numpy(dtype=float)
    assert set(np.unique(x)) == {14.0}                       # the premise, asserted not assumed
    assert balance.smd(x, a, w) == 0.0
    assert np.isnan(rejected_zero_variance_form(x, a, w))


def test_nihss_baseline_reads_zero_and_is_DEFINED_on_the_fixture_cohort():
    """§12.4a's real-frame witness, running with no `data/`. The row §5.3a is about."""
    _, _, _, bal = assessed()
    row = next(r for r in bal.covariates if r.covariate == "nihss_baseline")
    assert row.weighted == 0.0 and row.unweighted == 0.0
    assert np.isfinite(row.weighted)
    assert "nihss_baseline" not in bal.worst()[1]          # not among the undefined


def test_only_SOME_weight_vectors_round_and_that_is_what_decides_the_branch():
    """What makes the two weighted means differ is the WEIGHTS, not the arm size (§18).

    This test exists because the fixture above was originally declared with two weight vectors that
    both divide exactly, so it exercised nothing and passed. Arm size is not the mechanism: two
    records per arm bite or do not bite depending on which weights they carry.
    """
    exact = (np.array([0.1, 0.9]), np.array([0.3, 0.7]))
    rounds = (np.array([0.2, 0.8]), np.array([0.3, 0.7, 0.11]))
    for w in exact:
        assert np.average(np.full(len(w), 14.0), weights=w) == 14.0
    for w in rounds:
        assert np.average(np.full(len(w), 14.0), weights=w) != 14.0

    # both arms of two records, and the pair still differs — so it is not about the arm size
    x, a, w = constant_everywhere()
    assert int((a == 1.0).sum()) == int((a == 0.0).sum()) == 2
    assert np.average(x[a == 1.0], weights=w[a == 1.0]) != np.average(x[a == 0.0],
                                                                     weights=w[a == 0.0])


def test_the_one_ulp_fact_itself_and_why_indicators_hid_it():
    """§3.1, asserted directly — this is the test that explains the other two.

    `np.average` divides a weighted sum by a weight sum, and the rounding of that division depends
    on the weights. A vector of exact 1.0s or 0.0s is immune, which is why 17 of the fixture's 19
    rows could never have shown this and the one integer-valued constant covariate did.
    """
    w1, w0 = np.array([0.3, 0.7, 0.11]), np.array([0.1, 0.9])
    assert np.average(np.full(3, 14.0), weights=w1) != np.average(np.full(2, 14.0), weights=w0)
    assert np.average(np.full(3, 1.0), weights=w1) == 1.0
    assert np.average(np.full(2, 1.0), weights=w0) == 1.0
    assert np.average(np.full(3, 0.0), weights=w1) == 0.0


# --- 12.1  the balance set and the role column ------------------------------------------------------

def declared_level_count() -> int:
    """The row count the table must have, computed from the declarations rather than pinned at 19."""
    return sum(len(config.FACTOR_LEVELS[c]) if c in config.CATEGORICAL else 1
               for c in config.BALANCE_SET)


def test_the_balance_set_is_the_fourteen_declared_names_in_declaration_order():
    assert config.BALANCE_SET == config.PS_COVARIATES + config.BALANCE_ONLY
    assert len(config.BALANCE_SET) == 14


def test_the_table_has_one_row_per_declared_level_and_the_decomposition_is_fourteen_plus_five():
    """19 rows: 14 from PS_COVARIATES' expansion — 7 linear + 3 onset_type + 4 center — and 5 from
    BALANCE_ONLY. Asserted against the declarations, never against a literal 19.

    An earlier draft of the spec wrote the decomposition as 12 + 5, which is 17: twelve is
    `design_matrix`'s parameter count `X.shape[1] + 1` (`propensity.py:477`) and not any table's row
    count (§4.2, §20.2).
    """
    df, ps, _ = fitted()
    sub = df.loc[ps.in_model]
    rows = balance._table(sub, sub[config.TREATMENT].to_numpy(dtype=float),
                          ps.w.loc[ps.in_model].to_numpy(dtype=float), config.BALANCE_SET)
    assert len(rows) == declared_level_count()

    from_ps = sum(len(config.FACTOR_LEVELS[c]) if c in config.CATEGORICAL else 1
                  for c in config.PS_COVARIATES)
    from_balance_only = len(config.BALANCE_ONLY)          # none of the five is a factor
    assert (from_ps, from_balance_only) == (14, 5)
    assert from_ps + from_balance_only == len(rows)
    assert [row.covariate for row in rows][:2] == ["age", "sex"]     # BALANCE_SET's own order


def test_the_three_roles_partition_the_table_with_no_row_unlabelled():
    df, ps, _ = fitted()
    sub = df.loc[ps.in_model]
    rows = balance._table(sub, sub[config.TREATMENT].to_numpy(dtype=float),
                          ps.w.loc[ps.in_model].to_numpy(dtype=float), config.BALANCE_SET)
    counts = {role: sum(1 for row in rows if row.role == role)
              for role in ("propensity model", "negative control", "excluded [§6]")}
    assert counts == {"propensity model": 14, "negative control": 4, "excluded [§6]": 1}
    assert sum(counts.values()) == len(rows)


def test_penumbra_ml_is_in_the_table_and_is_NOT_a_negative_control():
    """[§6] excludes `penumbra_ml` because it is a deterministic function of `core_ml` and
    `tmax6_ml` — including all three makes the design exactly singular — and gives it a balance-table
    role of its own. A covariate excluded for exact collinearity is not a covariate weighting was
    ever expected to fix, so reporting it as a negative control would misstate what the row means.
    This is the roadmap correction of §19 in test form.
    """
    assert "penumbra_ml" in config.BALANCE_ONLY
    assert "penumbra_ml" not in config.NEGATIVE_CONTROLS
    assert balance._role("penumbra_ml") == "excluded [§6]"
    assert set(config.BALANCE_ONLY) - set(config.NEGATIVE_CONTROLS) == {"penumbra_ml"}


@pytest.mark.parametrize("name", config.NEGATIVE_CONTROLS)
def test_every_vascular_risk_factor_is_a_negative_control(name):
    assert balance._role(name) == "negative control"


def test_a_fifth_risk_factor_becomes_a_negative_control_by_being_declared():
    """§4.3's "cannot fail to", measured rather than described.

    `config.py` is executed from its own source with one literal patched, so what is exercised is the
    shipped declaration and not a re-spelling of it in this test. Standard library only, so the exec
    costs nothing and touches no filesystem.
    """
    source = Path(config.__file__).read_text(encoding="utf-8")
    needle = '"hypertension", "hyperlipidemia", "diabetes", "smoking")                    # [§13]'
    assert source.count(needle) == 1, "the PS_COVARIATES_FULL declaration moved; repoint this patch"
    patched = source.replace(needle, needle.replace(
        '"smoking")', '"smoking", "vasculitis")'))

    # A real module object registered in sys.modules, because `config.py` declares dataclasses and
    # `@dataclass` resolves its own module by name — a bare namespace dict raises inside dataclasses.
    patched_module = types.ModuleType("config_with_a_fifth_risk_factor")
    patched_module.__file__ = config.__file__
    sys.modules[patched_module.__name__] = patched_module
    try:
        exec(compile(patched, config.__file__, "exec"), patched_module.__dict__)
    finally:
        del sys.modules[patched_module.__name__]

    assert "vasculitis" in patched_module.NEGATIVE_CONTROLS
    assert patched_module.NEGATIVE_CONTROLS == config.NEGATIVE_CONTROLS + ("vasculitis",)
    # and BALANCE_SET does not move, because a [§13] addition is not a BALANCE_ONLY name
    assert patched_module.BALANCE_SET == config.BALANCE_SET


# --- 12.2  declared levels, including the reference and the absent one -------------------------------

@pytest.mark.parametrize("factor", config.CATEGORICAL)
def test_every_declared_level_of_every_factor_has_a_row(factor):
    # Parametrised over CATEGORICAL, so a third factor is covered by this test as it stands.
    df, ps, _ = fitted()
    labels = [label for label, _ in balance._levels(df.loc[ps.in_model], factor)]
    assert labels == [f"{factor} = {level}" for level in config.FACTOR_LEVELS[factor]]


def test_the_reference_levels_model_design_drops_are_PRESENT_here():
    """`center = HUG` and `onset_type = witnessed` are the two levels `design` drops by name.

    `center = HUG` carries the workbook's worst residual imbalance, so a table without it reports the
    worst among the covariates it kept as though it were the worst overall (§4.2).
    """
    df, ps, _ = fitted()
    sub = df.loc[ps.in_model]
    labels = [row.covariate for row in balance._table(
        sub, sub[config.TREATMENT].to_numpy(dtype=float),
        ps.w.loc[ps.in_model].to_numpy(dtype=float), config.BALANCE_SET)]
    for factor, level in config.REFERENCE_LEVELS.items():
        assert f"{factor} = {level}" in labels
    assert "center = HUG" in labels and "onset_type = witnessed" in labels


def test_a_declared_level_NOBODY_has_is_a_row_reading_zero():
    """A level that vanishes is a fact nobody sees. On `cohort_frame()`'s cohort two of the four
    declared centres have no record at all, and both are rows reading 0.0 in both columns — which is
    what §5.3's legitimate-zero branch is for, and what makes the [§3] restriction readable off the
    table.
    """
    df, ps, _, bal = assessed()
    present = set(df.loc[ps.in_model, "center"])
    absent = [c for c in config.CENTER_ORDER if c not in present]
    assert absent, "the fixture cohort is supposed to leave declared centres empty"
    for centre in absent:
        row = next(r for r in bal.covariates if r.covariate == f"center = {centre}")
        assert row.unweighted == 0.0 and row.weighted == 0.0


def design_column_names(sub: pd.DataFrame, names) -> list[str]:
    """The shortcut §4.2 forbids: a balance table built over `model.design`'s columns."""
    X, _ = model.design(sub, names)
    return list(X.columns)


@DATA_GATED
def test_the_design_shortcut_is_SIXTEEN_rows_and_omits_exactly_the_three_indicators(workbook):
    """§4.2 as a measurement rather than an argument (§12.2, DoD-6).

    Built over `BALANCE_SET` — the same fourteen declared names — so the entire difference is the
    ENCODING and nothing is attributable to the covariate list. Nothing raises: the shortcut produces
    a shorter, complete-looking table.
    """
    df, ps, _ = workbook
    sub = df.loc[ps.in_model]
    rows = balance._table(sub, sub[config.TREATMENT].to_numpy(dtype=float),
                          ps.w.loc[ps.in_model].to_numpy(dtype=float), config.BALANCE_SET)
    columns = design_column_names(sub, config.BALANCE_SET)
    assert len(rows) == 19 and len(columns) == 16

    # the table's labels in `design`'s own naming, so the comparison needs no second declaration
    missing = [row.covariate.replace(" = ", "_") for row in rows
               if row.covariate.replace(" = ", "_") not in columns]
    assert sorted(missing) == ["center_HUG", "center_USZ", "onset_type_witnessed"]


@DATA_GATED
def test_the_row_the_shortcut_misses_carries_the_workbook_worst_residual_imbalance(
        workbook, workbook_balance):
    df, ps, _ = workbook
    in_model_rows = [r for r in workbook_balance.covariates if r.role == "propensity model"]
    worst = max(in_model_rows, key=lambda r: abs(r.weighted))
    assert worst.covariate == "center = HUG"
    assert worst.covariate.replace(" = ", "_") not in design_column_names(
        df.loc[ps.in_model], config.BALANCE_SET)


@DATA_GATED
def test_the_form_the_pilots_actually_use_is_ELEVEN_rows_missing_eight_of_nineteen(workbook):
    """`pilots/analysis.py:249-265` builds its table over `design(df, C.PS_COVARIATES)`.

    Worse than the `BALANCE_SET` form by five further rows — the whole of `BALANCE_ONLY`. Both
    numbers are asserted, because the first isolates the encoding and the second measures the
    shortcut that exists in code today. Neither is 12, which is `design_matrix`'s parameter count.
    """
    df, ps, _ = workbook
    columns = design_column_names(df.loc[ps.in_model], config.PS_COVARIATES)
    assert len(columns) == 11
    assert 19 - len(columns) == 8


def test_the_levels_missingness_mask_is_UNREACHABLE_from_assess():
    """§21.3's finding: every `CATEGORICAL` name is a PS covariate, so a record with an absent factor
    value is complete-cased out of `in_model` before `_levels` ever sees it.

    The mask stays — it is correct, and Stage 12 will need it over a population `complete_cases` did
    not build — and what is asserted is the property that makes it unreachable, which fails if a
    third factor is ever declared outside the propensity model.
    """
    assert [c for c in config.CATEGORICAL if c not in config.PS_COVARIATES] == []

    # and the mask itself does what it claims, asserted directly on `_levels` rather than through a
    # frame no pipeline can produce
    frame = pd.DataFrame({"center": pd.Series(["HUG", None], dtype="string")})
    values = dict(balance._levels(frame, "center"))["center = HUG"]
    assert values.tolist()[0] == 1.0 and np.isnan(values.tolist()[1])


def test_the_three_pandas_facts_levels_rests_on():
    # Measured on pandas 2.3.3 and asserted here so no implementer re-derives them (§12.2).
    series = pd.Series(["HUG", "CHUV", None], dtype="string")
    assert (series == "HUG").dtype == "boolean"                       # nullable, not bool
    assert np.isnan((series == "HUG").astype(float).mask(series.isna()).to_numpy()[2])
    frame = pd.DataFrame({"x": [1, 2, 3]}, index=[0, 1, 2])
    assert len(frame.loc[pd.Series([True, False, pd.NA], dtype="boolean")]) == 1   # <NA> is False


# --- 12.5a  `_worst` and `_over_threshold`, over a hand-built list ------------------------------------

def rows_for_verdicts() -> tuple[balance.CovariateBalance, ...]:
    """Four rows, one of them undefined: the shape both verdicts are about."""
    return tuple(
        balance.CovariateBalance(covariate=name, role="propensity model", n=10, sd=1.0,
                                 unweighted=0.0, weighted=value)
        for name, value in (("age", 0.05), ("sex", -0.20), ("core_ml", 0.11),
                            ("center = USZ", np.nan)))


def test_worst_returns_the_largest_DEFINED_row_and_the_undefined_names_from_ONE_call():
    worst, undefined = balance._worst(rows_for_verdicts())
    assert worst.covariate == "sex" and worst.weighted == -0.20
    assert undefined == ("center = USZ",)


def test_the_two_wrong_answers_worst_exists_to_make_unavailable():
    """§3.1: pandas SKIPS the missing value and numpy POISONS the result — the same expression, two
    opposite wrong answers, and `worst()` returns the number AND the name so neither can be reported
    alone.
    """
    column = [row.weighted for row in rows_for_verdicts()]
    assert pd.Series(column).abs().max() == 0.20              # skips the undefined row silently
    assert np.isnan(np.abs(np.array(column)).max())           # poisons the whole column
    worst, undefined = balance._worst(rows_for_verdicts())
    assert abs(worst.weighted) == 0.20 and undefined


def test_over_threshold_excludes_the_undefined_rows_and_the_two_collections_partition_the_table():
    rows = rows_for_verdicts()
    over = balance._over_threshold(rows)
    _, undefined = balance._worst(rows)
    assert over == ("sex", "core_ml")                # 0.11 reaches 0.10; 0.05 does not
    assert set(over) & set(undefined) == set()
    balanced = [r.covariate for r in rows
                if np.isfinite(r.weighted) and abs(r.weighted) < config.SMD_THRESHOLD]
    assert set(over) | set(undefined) | set(balanced) == {r.covariate for r in rows}


def test_worst_is_empty_handed_on_a_table_with_no_defined_row():
    rows = (balance.CovariateBalance("age", "propensity model", 3, np.nan, np.nan, np.nan),)
    assert balance._worst(rows) == (None, ("age",))
    assert balance._over_threshold(rows) == ()


# --- 12.6/12.7  the centre rows, asserted directly on `_centre` and `_status` -------------------------

def pooled_yardstick(df: pd.DataFrame, ps: propensity.Propensity) -> dict[str, float]:
    """§5.5's one yardstick, harvested from the table over `in_model` exactly as §8 harvests it."""
    sub = df.loc[ps.in_model]
    rows = balance._table(sub, sub[config.TREATMENT].to_numpy(dtype=float),
                          ps.w.loc[ps.in_model].to_numpy(dtype=float), config.BALANCE_SET)
    return {row.covariate: row.sd for row in rows}


def one_armed_centre():
    """A centre whose only control is complete-cased out of `in_model` — the `_NO_CONTRAST` case.

    Constructed, because no frame reaches it: Stage 5's P2 guarantees both arms at every retained
    centre, and on both the workbook and the fixture `in_model` empties none of them (§6.3). Blanking
    a [§6] covariate on CHUV's only control is the one-line edit that does.
    """
    frame, audit = built()
    frame = set_cell(frame, "COHORT-3", "core_ml", np.nan)
    return frame, propensity.fit(frame, audit), audit


def test_status_returns_each_of_the_three_and_there_is_no_fourth():
    index = pd.RangeIndex(3)
    at_none = pd.Series([False, False, False], index=index)
    at_all = pd.Series([True, True, True], index=index)
    both_arms = {0: np.array([True, False, False]), 1: np.array([False, True, True])}
    one_arm = {0: np.array([False, False, False]), 1: np.array([True, True, True])}

    assert balance._status(at_none, both_arms) == balance._STRUCTURAL
    assert balance._status(at_all, one_arm) == balance._NO_CONTRAST
    assert balance._status(at_all, both_arms) == balance._REPORTED
    # the three module constants are the whole of the column's vocabulary
    assert len({balance._REPORTED, balance._STRUCTURAL, balance._NO_CONTRAST}) == 3


def test_the_structural_branch_beats_the_empty_arm_branch_and_the_order_is_the_specification():
    # A centre with no cohort record has no arm either, so both conditions hold; [§3] removed it and
    # that is the fact the column must carry, not complete-casing's consequence.
    at_none = pd.Series([False, False], index=pd.RangeIndex(2))
    empty = {0: np.array([False, False]), 1: np.array([False, False])}
    assert balance._status(at_none, empty) == balance._STRUCTURAL


def test_centre_builds_one_row_per_declared_centre_including_the_ones_with_no_record():
    df, ps, _ = fitted()
    sds = pooled_yardstick(df, ps)
    rows = [balance._centre(df, ps, centre, sds) for centre in config.CENTER_ORDER]
    assert [row.centre for row in rows] == list(config.CENTER_ORDER)

    present = set(df["center"])
    for row in rows:
        if row.centre in present:
            assert row.status == balance._REPORTED and row.in_cohort > 0
        else:
            assert row.status == balance._STRUCTURAL
            assert row.in_cohort == 0 and row.weighted == 0
            assert np.isnan(row.max_weight) and np.isnan(row.weight_share)
            assert all(np.isnan(row.ess[code]) for code in config.TREATMENT_LABELS)
            assert np.isnan(row.worst)


def test_an_arm_of_ONE_is_not_an_empty_arm_and_the_two_land_in_different_columns():
    """§6.3, on the fixture cohort where every centre has exactly one control.

    The ESS and the propensity range are computable — Kish of a single weight is 1.0 — while every
    within-centre SMD is undefined, because §5.3's first branch needs two records per arm. So the
    row reads `overlap reported` with a `missing` worst: the centre's overlap IS reported and its
    balance is not judgeable, which are two different facts in two different columns.
    """
    df, ps, _ = fitted()
    sds = pooled_yardstick(df, ps)
    reported = [balance._centre(df, ps, c, sds) for c in config.CENTER_ORDER
                if int((df["center"] == c).sum())]
    assert reported, "the fixture cohort is supposed to retain at least one centre"
    for row in reported:
        assert row.n[0] == 1                                   # exactly one control
        assert row.ess[0] == pytest.approx(1.0)                # Kish of a single weight
        assert all(np.isfinite(v) for v in row.e_range[0])
        assert np.isnan(row.worst)                             # and the balance is not judgeable
        assert row.status == balance._REPORTED


def test_the_pooled_row_comes_from_the_same_function_and_carries_no_status_of_its_own():
    df, ps, _ = fitted()
    sds = pooled_yardstick(df, ps)
    pooled = balance._centre(df, ps, None, sds)
    assert pooled.centre == "all (pooled)"
    assert pooled.status == balance._REPORTED           # no fourth literal (§6.2, §7.4)
    assert pooled.in_cohort == len(df)
    assert pooled.weighted == int(ps.in_model.sum())
    assert pooled.weight_share == 1.0                   # §8 divides a sum by itself
    for code in config.TREATMENT_LABELS:
        assert pooled.ess[code] == pytest.approx(ps.ess[code])


def test_a_centre_rows_sd_is_the_POOLED_one_and_never_the_centres_own():
    """§5.5's one yardstick: the numerator is the centre's and the denominator is the whole weighted
    set's. Asserted by name and exactly, because a centre row computed against its own variances is
    finite, plausible and up to 46% out.
    """
    df, ps, _ = fitted()
    sds = pooled_yardstick(df, ps)
    centre = next(c for c in config.CENTER_ORDER if int((df["center"] == c).sum()))
    sub = df.loc[(df["center"] == centre) & ps.in_model]
    rows = balance._table(
        sub, sub[config.TREATMENT].to_numpy(dtype=float),
        ps.w.loc[(df["center"] == centre) & ps.in_model].to_numpy(dtype=float),
        [c for c in config.BALANCE_SET if c != "center"], sds)
    for row in rows:
        assert (row.sd == sds[row.covariate]) or (np.isnan(row.sd) and np.isnan(sds[row.covariate]))
        # and every label really is a key of `sds`, so `_table`'s .get fallback is measured inert
        assert row.covariate in sds


def test_the_centre_indicators_are_dropped_INSIDE_a_centre_and_judged_pooled():
    df, ps, _ = fitted()
    sds = pooled_yardstick(df, ps)
    centre = next(c for c in config.CENTER_ORDER if int((df["center"] == c).sum()))
    at = (df["center"] == centre) & ps.in_model
    names = [c for c in config.BALANCE_SET if c != "center"]
    labels = [row.covariate for row in balance._table(
        df.loc[at], df.loc[at, config.TREATMENT].to_numpy(dtype=float),
        ps.w.loc[at].to_numpy(dtype=float), names, sds)]
    assert not any(label.startswith("center = ") for label in labels)
    # and pooled they are judged: the whole declared set, centre levels included
    pooled_labels = [row.covariate for row in balance._table(
        df.loc[ps.in_model], df.loc[ps.in_model, config.TREATMENT].to_numpy(dtype=float),
        ps.w.loc[ps.in_model].to_numpy(dtype=float), config.BALANCE_SET, sds)]
    assert "center = HUG" in pooled_labels


def test_ess_is_NOT_called_on_an_empty_arm_and_stage_sixes_raise_is_intact_behind_the_guard():
    """§6.3's two halves. The guard is shown to be doing the work, and the raise is shown to be there.

    Without the guard the diagnostic aborts mid-table on a centre complete-casing emptied — which is
    a finding to report, not an exception to propagate.
    """
    df, ps, _ = one_armed_centre()
    sds = pooled_yardstick(df, ps)
    row = balance._centre(df, ps, "CHUV", sds)

    assert row.status == balance._NO_CONTRAST
    assert row.n[0] == 0 and row.n[1] > 0
    assert np.isnan(row.ess[0]) and np.isfinite(row.ess[1])
    assert all(np.isnan(v) for v in row.e_range[0])
    assert np.isnan(row.worst)
    assert row.in_cohort > row.weighted            # the control is in the cohort and not weighted

    # the companion: Stage 6's raise, on that same arm
    emptied = ps.w[(df["center"] == "CHUV") & ps.in_model & (df[config.TREATMENT] == 0)]
    with pytest.raises(model.FitError) as excinfo:
        propensity.ess(emptied)
    assert "in_model mask lost one" in message_of(excinfo)


def test_in_cohort_and_weighted_are_read_from_at_and_at_and_in_model():
    df, ps, _ = one_armed_centre()
    sds = pooled_yardstick(df, ps)
    row = balance._centre(df, ps, "CHUV", sds)
    at = df["center"] == "CHUV"
    assert row.in_cohort == int(at.sum())
    assert row.weighted == int((at & ps.in_model).sum())
    assert row.in_cohort != row.weighted           # the two columns are not one column


def test_the_defined_weight_shares_sum_to_one_and_a_structural_centres_share_is_nan_not_zero():
    """Both qualifications are load-bearing (§12.6): a structurally non-positive centre's share is
    `nan` — `len(w)` is 0, so §8 returns nan and `_fmt` renders `missing` — and a `no weighted
    contrast` centre carries a positive share while not being a REPORTED centre, so "over the
    reported centres" would be the wrong index set.
    """
    df, ps, _ = fitted()
    sds = pooled_yardstick(df, ps)
    rows = [balance._centre(df, ps, c, sds) for c in config.CENTER_ORDER]
    defined = [row.weight_share for row in rows if np.isfinite(row.weight_share)]
    structural = [row for row in rows if row.status == balance._STRUCTURAL]
    assert structural and all(np.isnan(row.weight_share) for row in structural)
    assert sum(defined) == pytest.approx(1.0, abs=1e-12)
    assert balance._centre(df, ps, None, sds).weight_share == 1.0


# --- 12.0.2  the golden balance vector ---------------------------------------------------------------
#
# Pinned as literals from the SYNTHETIC frame, so that no patient-derived number enters git (Stage 6
# §7.3). Asserted to 1e-6, which is a regression pin against an edit to our own arithmetic: pinning
# at machine precision would make it fail on a numpy patch release rather than on a mistake.
#
# The assertion is BY COVARIATE NAME, one row at a time. `center = HUG` and `center = CHUV` differ
# only in sign, so an assertion against a sorted list of values would pass while confusing them.
#
# Every figure is `_fmt`'s six significant figures, correctly ROUNDED — which is what makes the 1e-6
# tolerance the right one. An earlier draft of §12.0.2 wrote `penumbra_ml` as 1.630960, truncated at
# the seventh digit and 1.27e-6 from the true 1.6309612667875129, so a test pinning that literal at
# the stated tolerance failed. §18 carries the measurement.

GOLDEN_BALANCE = {
    #                    unweighted    weighted
    "core_ml":         (0.723747,    0.635113),
    "tmax6_ml":        (2.248299,    1.839089),
    "center = HUG":    (0.258199,    0.139181),
    "center = CHUV":   (-0.258199,  -0.139181),
    "penumbra_ml":     (1.630961,    1.318924),
}


@pytest.mark.parametrize("covariate", sorted(GOLDEN_BALANCE))
def test_the_golden_balance_vector_by_covariate_name(covariate):
    _, _, _, bal = assessed()
    row = next(r for r in bal.covariates if r.covariate == covariate)
    unweighted, weighted = GOLDEN_BALANCE[covariate]
    assert row.unweighted == pytest.approx(unweighted, abs=1e-6)
    assert row.weighted == pytest.approx(weighted, abs=1e-6)


def test_the_fixture_table_is_nineteen_rows_fourteen_of_them_exactly_zero_and_none_undefined():
    _, _, _, bal = assessed()
    assert len(bal.covariates) == declared_level_count()
    assert sum(1 for r in bal.covariates
               if r.unweighted == 0.0 and r.weighted == 0.0) == 14
    worst, undefined = bal.worst()
    assert undefined == ()
    assert worst.covariate == "tmax6_ml"
    assert abs(worst.weighted) == pytest.approx(1.839089, abs=1e-6)


# --- 12.0.3  the workbook fixture returns three things from ONE audit ----------------------------------

@pytest.fixture(scope="module")
def workbook():
    """(df, ps, audit) from ONE linear run against ONE Audit. §12.0.3.

    `test_propensity.py`'s `workbook_ps` deliberately does the opposite — it refits against a FRESH
    Audit, because that module's tests are about `fit`'s own four entries and a shared audit would
    make their positions depend on what ran before. Copying that idiom here would put
    `overlap_weights` in one Audit and `balance_smd` in another, and §12.8's reconciliation would
    have nothing to read. The two fixtures look alike and differ in the one way that matters.
    """
    df, audit = data.load(data.WORKBOOK)
    df = derive.derive(df, audit)
    df = eligibility.classify(df, audit)
    df = cohort.build(df, audit)
    return df, propensity.fit(df, audit), audit


@pytest.fixture(scope="module")
def workbook_balance(workbook):
    """`assess` run ONCE over the fixture above, so the reading tests share one call."""
    df, ps, audit = workbook
    return balance.assess(df, ps, audit)


# --- 12.2a  the per-row denominator [§11] ---------------------------------------------------------------
#
# The column §4.4 argues for, tested on a CONSTRUCTED frame because no available frame varies it:
# every workbook row reads 92 and every fixture row reads 5. Without this section the `n` column is a
# constant everywhere in the suite, and an implementation masking globally instead of per row is
# green.

def blanked_penumbra():
    """One `in_model` record's `penumbra_ml` blanked — a BALANCE_ONLY covariate, so `in_model` cannot
    move. That non-movement is the second half of the test and the contrast that makes the first half
    mean something."""
    frame, audit = built()
    frame = set_cell(frame, "HAND-1", "penumbra_ml", np.nan)
    return frame, propensity.fit(frame, audit), audit


def test_a_blanked_balance_only_covariate_gives_ITS_row_a_smaller_n_and_leaves_the_others(recwarn):
    df, ps, audit = blanked_penumbra()
    bal = balance.assess(df, ps, audit)
    rows = {row.covariate: row for row in bal.covariates}
    assert rows["penumbra_ml"].n == 4
    assert all(row.n == 5 for name, row in rows.items() if name != "penumbra_ml")


def test_in_model_does_NOT_move_when_a_balance_only_covariate_is_blanked():
    # Half the test: without it, an implementation complete-casing on BALANCE_SET would pass the
    # first assertion by dropping the record from every row.
    df, ps, _ = blanked_penumbra()
    assert int(ps.in_model.sum()) == 5
    assert ps.in_model.all()


def test_the_row_is_computed_over_ITS_OWN_records_and_its_own_pooled_sd():
    df, ps, audit = blanked_penumbra()
    bal = balance.assess(df, ps, audit)
    row = next(r for r in bal.covariates if r.covariate == "penumbra_ml")

    values = df.loc[ps.in_model, "penumbra_ml"].to_numpy(dtype=float)
    present = np.isfinite(values)
    x = values[present]
    a = df.loc[ps.in_model, config.TREATMENT].to_numpy(dtype=float)[present]
    w = ps.w.loc[ps.in_model].to_numpy(dtype=float)[present]
    assert row.weighted == balance.smd(x, a, w)          # the same number, not merely close
    assert row.sd == balance._pooled_sd(x, a)            # a different denominator from every other row
    assert row.sd != next(r.sd for r in bal.covariates if r.covariate == "tmax6_ml")


def test_blanking_a_PS_COVARIATE_instead_moves_in_model_AND_every_rows_n_together():
    # The contrast. This is what makes the three assertions above mean something rather than pass.
    frame, audit = built()
    frame = set_cell(frame, "HAND-1", "core_ml", np.nan)
    ps = propensity.fit(frame, audit)
    bal = balance.assess(frame, ps, audit)
    assert int(ps.in_model.sum()) == 4                   # in_model DID move
    assert {row.n for row in bal.covariates} == {4}      # and every row's n moved with it


# --- 12.5  the threshold is a verdict and `assess` does not raise on it ---------------------------------

def test_a_grossly_imbalanced_frame_returns_normally_and_names_its_rows():
    """§5.4. Residual imbalance is a finding for [§16] and, above the threshold, a question for the
    PI under [§13] — not a malformed input. `assess` raises only on B1-B5, which are all statements
    about the FRAME rather than about the DATA.

    The fixture cohort is grossly imbalanced already: its worst weighted |SMD| is 1.84, eighteen
    times [§9]'s threshold, and five of its nineteen rows are at or above it.
    """
    _, _, _, bal = assessed()                            # no pytest.raises anywhere on this path
    assert bal.unbalanced() == ("core_ml", "tmax6_ml", "center = HUG", "center = CHUV",
                                "penumbra_ml")
    assert abs(bal.worst()[0].weighted) > 10 * config.SMD_THRESHOLD


def balance_with(rows):
    """A `Balance` over hand-built covariate rows, with a real `CentreOverlap` for `pooled`."""
    df, ps, _ = fitted()
    pooled = balance._centre(df, ps, None, pooled_yardstick(df, ps))
    return balance.Balance(covariates=tuple(rows), centres=(), pooled=pooled)


def test_unbalanced_excludes_the_undefined_rows_and_worst_returns_them_instead():
    bal = balance_with(rows_for_verdicts())
    assert "center = USZ" not in bal.unbalanced()
    assert bal.worst()[1] == ("center = USZ",)
    assert bal.worst()[0].covariate == "sex"
    # the two collections partition the table with the defined-and-balanced rows
    assert set(bal.unbalanced()) | set(bal.worst()[1]) | {"age"} == {
        row.covariate for row in bal.covariates}


def test_the_balance_methods_delegate_to_the_module_level_verdicts():
    rows = rows_for_verdicts()
    bal = balance_with(rows)
    assert bal.unbalanced() == balance._over_threshold(rows)
    assert bal.worst()[1] == balance._worst(rows)[1]


def test_the_detail_names_EXACTLY_the_rows_unbalanced_returns():
    """One threshold predicate, two callers (§3). The failure this invites is a `detail` string
    naming four rows above a table flagging five."""
    _, _, audit, bal = assessed()
    detail = audit.entry("model", "balance_smd").detail
    assert f"{len(bal.unbalanced())} row(s) reach" in detail
    for name in bal.unbalanced():
        assert name in detail
    # and the table agrees with both: one `no` cell per named row
    table = audit.entry("model", "balance_smd").table
    assert sum(1 for row in table[1:] if row[-1] == "no") == len(bal.unbalanced())


def test_the_threshold_is_read_from_config_and_never_from_a_literal(monkeypatch):
    _, _, audit, bal = assessed()
    before = {row[0]: row[-1] for row in audit.entry("model", "balance_smd").table[1:]}
    assert before["center = HUG"] == "no"                # |SMD| 0.139 against 0.10

    monkeypatch.setattr(config, "SMD_THRESHOLD", 0.5)
    _, _, audit_patched, bal_patched = assessed()
    after = {row[0]: row[-1] for row in audit_patched.entry("model", "balance_smd").table[1:]}
    assert after["center = HUG"] == "yes"                # the verdict moved with the constant
    assert "center = HUG" not in bal_patched.unbalanced()
    assert audit_patched.entry("model", "balance_smd").table[0][-1] == "abs SMD < 0.5"


# --- 12.6  the within-centre table, through `assess` ---------------------------------------------------

def test_centres_is_exactly_CENTER_ORDER_and_pooled_is_a_field_beside_it():
    _, _, _, bal = assessed()
    assert [row.centre for row in bal.centres] == list(config.CENTER_ORDER)
    assert len(bal.centres) == len(config.CENTER_ORDER)
    assert isinstance(bal.pooled, balance.CentreOverlap)
    assert bal.pooled not in bal.centres
    assert "all (pooled)" not in [row.centre for row in bal.centres]


def test_the_RENDERED_table_is_those_rows_plus_the_pooled_one_in_that_order():
    # The return value and the table have different shapes on purpose, and both are pinned: the
    # log's reader wants one grid with a labelled row, a Python caller wants `centres` to mean the
    # declared centres and nothing else (§3, §6.1).
    _, _, audit, _ = assessed()
    rendered = [row[0] for row in audit.entry("model", "overlap_by_centre").table[1:]]
    assert rendered == [*config.CENTER_ORDER, "all (pooled)"]


@DATA_GATED
def test_in_cohort_and_weighted_DIFFER_for_lugano_on_the_workbook(workbook, workbook_balance):
    df, ps, _ = workbook
    lugano = next(row for row in workbook_balance.centres if row.centre == "Lugano")
    assert (lugano.in_cohort, lugano.weighted) == (31, 30)
    assert lugano.in_cohort == int((df["center"] == "Lugano").sum())
    assert lugano.weighted == int(((df["center"] == "Lugano") & ps.in_model).sum())


@DATA_GATED
def test_every_arms_e_range_lies_inside_the_pooled_range_and_the_pooled_one_is_stage_sixes(
        workbook, workbook_balance):
    pooled = workbook_balance.pooled
    low = min(pooled.e_range[code][0] for code in config.TREATMENT_LABELS)
    high = max(pooled.e_range[code][1] for code in config.TREATMENT_LABELS)
    for row in workbook_balance.centres:
        for code in config.TREATMENT_LABELS:
            lo, hi = row.e_range[code]
            if np.isfinite(lo):
                assert low <= lo <= hi <= high

    # one number in one place (Stage 6 §20.2's rule): the range Stage 6 already published
    _, _, audit = workbook
    detail = audit.entry("model", "propensity_fit").detail
    assert f"range {data._fmt(low)} to {data._fmt(high)}" in detail


@DATA_GATED
def test_the_within_centre_worst_is_on_the_POOLED_yardstick_and_the_local_one_would_differ(
        workbook, workbook_balance):
    """§5.5 and DoD-14, pinned by a test rather than by a comment.

    Under each centre's own variances no two rows of that column share a yardstick, none shares one
    with `SMD after` in the entry above it, and the numbers are up to 46% out. What does NOT change
    is the covariate identified as worst at each centre, or the order of the three centres — which
    is what makes the change one of scale and not of finding.
    """
    df, ps, _ = workbook
    names = [c for c in config.BALANCE_SET if c != "center"]
    expected = {"HUG": 0.4956, "CHUV": 1.4100, "Lugano": -1.1087}
    expected_local = {"HUG": 0.5385, "CHUV": 2.0608, "Lugano": -1.5275}
    expected_covariate = {"HUG": "diabetes", "CHUV": "onset_type = unwitnessed",
                          "Lugano": "hyperlipidemia"}
    sds = pooled_yardstick(df, ps)

    for row in workbook_balance.centres:
        if row.status != balance._REPORTED:
            continue
        at = (df["center"] == row.centre) & ps.in_model
        arms = df.loc[at, config.TREATMENT].to_numpy(dtype=float)
        weights = ps.w.loc[at].to_numpy(dtype=float)
        pooled_rows = balance._table(df.loc[at], arms, weights, names, sds)
        local_rows = balance._table(df.loc[at], arms, weights, names)      # sds=None — the defect

        assert row.worst == pytest.approx(expected[row.centre], abs=5e-4)
        worst_pooled = max((r for r in pooled_rows if np.isfinite(r.weighted)),
                           key=lambda r: abs(r.weighted))
        worst_local = max((r for r in local_rows if np.isfinite(r.weighted)),
                          key=lambda r: abs(r.weighted))
        assert worst_pooled.weighted == pytest.approx(row.worst, abs=1e-12)
        assert worst_local.weighted == pytest.approx(expected_local[row.centre], abs=5e-4)
        # the same covariate is identified either way; only the scale moves
        assert worst_pooled.covariate == worst_local.covariate == expected_covariate[row.centre]
        assert abs(worst_local.weighted) > abs(worst_pooled.weighted)


@DATA_GATED
def test_a_centre_rows_sd_equals_the_pooled_rows_for_the_same_covariate_on_the_workbook(workbook):
    df, ps, _ = workbook
    sds = pooled_yardstick(df, ps)
    at = (df["center"] == "HUG") & ps.in_model
    rows = balance._table(df.loc[at], df.loc[at, config.TREATMENT].to_numpy(dtype=float),
                          ps.w.loc[at].to_numpy(dtype=float),
                          [c for c in config.BALANCE_SET if c != "center"], sds)
    for row in rows:
        assert row.sd == sds[row.covariate]


@DATA_GATED
def test_dropping_the_centre_indicators_inside_a_centre_changes_no_centres_worst(
        workbook, workbook_balance):
    # Measured inert on this data, and kept for its reason rather than its effect: a centre's own
    # indicators are constant there by construction, so those rows can only ever read 0.0.
    df, ps, _ = workbook
    sds = pooled_yardstick(df, ps)
    for row in workbook_balance.centres:
        if row.status != balance._REPORTED:
            continue
        at = (df["center"] == row.centre) & ps.in_model
        with_centre = balance._table(
            df.loc[at], df.loc[at, config.TREATMENT].to_numpy(dtype=float),
            ps.w.loc[at].to_numpy(dtype=float), config.BALANCE_SET, sds)
        defined = [r.weighted for r in with_centre if np.isfinite(r.weighted)]
        assert max(defined, key=abs) == pytest.approx(row.worst, abs=1e-12)


# --- 12.7  structural non-positivity, and `ess` is never allowed to raise --------------------------------

@DATA_GATED
def test_USZ_has_a_row_whose_status_is_the_structural_one_and_whose_cells_are_all_missing(
        workbook_balance, workbook):
    usz = next(row for row in workbook_balance.centres if row.centre == "USZ")
    assert usz.status == balance._STRUCTURAL            # the STRING, because Stage 14 reads it
    assert usz.in_cohort == 0 and usz.weighted == 0

    _, _, audit = workbook
    rendered = rows_of(audit.entry("model", "overlap_by_centre"))["USZ"]
    assert rendered[-1] == balance._STRUCTURAL
    assert rendered[:4] == ("0", "0", "0", "0")         # in cohort, weighted, and both arms
    assert all(cell == "missing" for cell in rendered[4:-1])   # every overlap cell, worst included


def test_the_status_values_over_the_whole_table_are_a_subset_of_the_declared_three():
    _, _, _, bal = assessed()
    declared = {balance._REPORTED, balance._STRUCTURAL, balance._NO_CONTRAST}
    assert {row.status for row in (*bal.centres, bal.pooled)} <= declared
    assert bal.pooled.status == balance._REPORTED       # what makes the claim true of the TABLE


def test_two_fixture_centres_are_structural_and_two_report_with_a_MISSING_worst():
    df, _, audit, bal = assessed()
    structural = [row for row in bal.centres if row.status == balance._STRUCTURAL]
    reported = [row for row in bal.centres if row.status == balance._REPORTED]
    assert len(structural) == 2 and len(reported) == 2
    for row in reported:
        assert np.isnan(row.worst)                      # an arm of one: not judgeable
        assert np.isfinite(row.ess[0]) and np.isfinite(row.ess[1])
    # and it renders as `missing`, not as 0.0 — asserting only the status would pass on an
    # implementation reporting a worst of zero
    rendered = rows_of(audit.entry("model", "overlap_by_centre"))
    for row in reported:
        assert rendered[row.centre][-2] == "missing"


def test_assess_does_NOT_raise_when_complete_casing_empties_a_centres_arm():
    df, ps, audit = one_armed_centre()
    bal = balance.assess(df, ps, audit)                 # no FitError mid-table
    chuv = next(row for row in bal.centres if row.centre == "CHUV")
    assert chuv.status == balance._NO_CONTRAST
    assert np.isnan(chuv.ess[0]) and np.isnan(chuv.worst)
    assert "CHUV" in [row.centre for row in bal.centres]          # still present in the table
    assert "CHUV" in rows_of(audit.entry("model", "overlap_by_centre"))
    # and the entry's own count excludes it: `n` is the centres that contributed a contrast
    assert audit.entry("model", "overlap_by_centre").n == sum(
        1 for row in bal.centres if row.status == balance._REPORTED)


def kish_divisions(source: str) -> list[int]:
    """Line numbers where a squared sum is divided by something — a Kish sum written by hand."""
    return [node.lineno for node in ast.walk(ast.parse(source))
            if isinstance(node, ast.BinOp) and isinstance(node.op, ast.Div)
            and any(isinstance(inner, ast.BinOp) and isinstance(inner.op, ast.Pow)
                    and isinstance(inner.right, ast.Constant) and inner.right.value == 2
                    for inner in ast.walk(node))]


def test_no_second_effective_sample_size_is_computed_in_this_module():
    # `propensity.ess` is public precisely so this stage does not write a second Kish sum
    # (Stage 6 §3), and `assess` calls it per centre per arm.
    assert kish_divisions(SOURCE) == []
    assert "propensity.ess" in SOURCE


def test_the_kish_scan_fires_on_a_pasted_kish_sum():
    # A scan that silently matches nothing would otherwise pass as a green test.
    assert kish_divisions("x = float(np.sum(v) ** 2 / np.sum(v ** 2))\n") == [1]


# --- 12.8  the pooled row reconciles with Stage 6 ---------------------------------------------------------
#
# Every assertion here reads `Balance.pooled` and, on the other side, the `overlap_weights` entry
# pulled out of the SAME Audit (§12.0.3) rather than recomputed.
#
# **Every comparison goes through `_fmt`, and that is not a formality.** An `AuditEntry`'s table holds
# STRINGS — Stage 6 built those cells with `_fmt` (`propensity.py:333-337`) — so comparing a
# full-precision float against six significant figures fails on the seventh digit. What is exact is
# the RENDERING: two entries of the same log agreeing to the precision the log prints is the property
# a reader can check, and the strongest one available without recomputing what Stage 6 published.

@DATA_GATED
def test_the_pooled_ess_equals_propensity_ess_and_stage_sixes_rendered_cells(
        workbook, workbook_balance):
    _, ps, audit = workbook
    weights_table = rows_of(audit.entry("model", "overlap_weights"))
    for code, label in config.TREATMENT_LABELS.items():
        assert workbook_balance.pooled.ess[code] == ps.ess[code]          # exactly, not approx
        n_cell, _sum_w, ess_cell = weights_table[label][:3]
        assert data._fmt(workbook_balance.pooled.ess[code]) == ess_cell
        assert str(workbook_balance.pooled.n[code]) == n_cell


@DATA_GATED
def test_the_pooled_max_weight_equals_stage_sixes_pooled_cell(workbook, workbook_balance):
    _, _, audit = workbook
    pooled_row = rows_of(audit.entry("model", "overlap_weights"))["all (pooled)"]
    assert data._fmt(workbook_balance.pooled.max_weight) == pooled_row[3]   # `max w`
    assert data._fmt(workbook_balance.pooled.weight_share) == pooled_row[4]


@DATA_GATED
def test_the_pooled_worst_is_the_balance_tables_worst_computed_by_the_same_call(workbook_balance):
    # One quantity computed twice, in two places, by the same `_table` call with the same arguments
    # (§8). The duplication is the price of the one-code-path property, and this assertion is what
    # makes it evidence rather than an assumption — a shared computation would agree by construction.
    assert workbook_balance.pooled.worst == pytest.approx(
        workbook_balance.worst()[0].weighted, abs=1e-12)


# --- 12.9  the two audit entries, and the log --------------------------------------------------------------

STEPS = ("balance_smd", "overlap_by_centre")


def test_the_two_entries_appear_in_order_at_positions_captured_BEFORE_the_call():
    # Never a tail slice of the form [-n:], which is the repair Stage 5 §10 landed and which every
    # stage since has carried: any two trailing entries satisfy such an assertion.
    df, ps, audit = fitted()
    before = len(audit.entries)
    balance.assess(df, ps, audit)
    added = audit.entries[before:]
    assert tuple(entry.step for entry in added) == STEPS
    assert {entry.kind for entry in added} == {"model"}


def test_data_kinds_is_STILL_the_declared_nine_and_nothing_new_renders():
    # §7.1's claim in test form, and the assertion that fails if someone adds a `balance` kind.
    assert data.KINDS == (
        "provenance", "contract", "correction", "observation",
        "derivation", "cohort", "model", "structural", "missingness")
    assert len(data.KINDS) == 9
    assert set(data._HEADINGS) == set(data.KINDS)
    assert data._HEADINGS["model"] == "Fitted models"


def test_both_entries_render_under_fitted_models_after_overlap_weights():
    _, _, audit, _ = assessed()
    rendered = audit.to_markdown()
    assert "## Fitted models" in rendered
    positions = [rendered.index(f"- **{step}**")
                 for step in ("overlap_weights", *STEPS)]
    assert positions == sorted(positions)


def test_the_entry_n_is_not_the_row_count_in_either_entry():
    _, ps, audit, bal = assessed()
    smd_entry = audit.entry("model", "balance_smd")
    overlap_entry = audit.entry("model", "overlap_by_centre")

    assert smd_entry.n == int(ps.in_model.sum())                 # records judged
    assert smd_entry.n != len(smd_entry.table) - 1               # and not the 19 rows
    assert overlap_entry.n == sum(1 for row in bal.centres if row.status == balance._REPORTED)
    assert overlap_entry.n != len(overlap_entry.table) - 1       # 2 reported against 5 rows


@DATA_GATED
def test_the_workbook_entry_n_matches_propensity_fits_n(workbook, workbook_balance):
    _, _, audit = workbook
    assert audit.entry("model", "balance_smd").n == audit.entry("model", "propensity_fit").n
    assert audit.entry("model", "overlap_by_centre").n == 3      # three centres contributed


def test_the_verdict_column_holds_only_yes_no_and_undefined():
    df, ps, audit = one_armed_centre()
    balance.assess(df, ps, audit)
    verdicts = {row[-1] for row in audit.entry("model", "balance_smd").table[1:]}
    assert verdicts <= {"yes", "no", "undefined"}
    for forbidden in ("True", "False", "1", "0", "nan"):
        assert forbidden not in verdicts


def test_an_undefined_smd_renders_as_the_literal_missing_and_never_as_nan():
    rows = (balance.CovariateBalance("age", "propensity model", 5, np.nan, np.nan, np.nan),)
    table = balance._smd_table(rows)
    assert table[1] == ("age", "propensity model", "5", "missing", "missing", "undefined")
    assert balance._range((np.nan, np.nan)) == "missing"         # one `missing`, not two
    assert balance._range((0.257067, 0.257067)) == "0.257067–0.257067"


def test_no_cell_of_either_table_contains_a_numpy_repr():
    _, _, audit, _ = assessed()
    for step in STEPS:
        cells = [cell for row in audit.entry("model", step).table for cell in row]
        assert not [c for c in cells if "np." in c or "array(" in c or "float64" in c]


@pytest.mark.parametrize("step", STEPS)
def test_no_cell_holds_a_pipe_and_the_rendered_grid_is_rectangular(step):
    """DoD-13. `data._md_table` does no escaping and sizes the separator from `len(rows[0])`
    (data.py:242-245), so one pipe inside a header cell gives a header row with more markdown cells
    than its own separator and body. Measured before the two columns were renamed: `balance_smd`
    rendered 8 cells against a 6-cell separator and `overlap_by_centre` 15 against 13, **in every log
    this stage writes**, and byte-identity across hash seeds passed the whole time — a table
    identically broken under both seeds is still identical. Look at the log, not only at its hash.
    """
    _, _, audit, _ = assessed()
    table = audit.entry("model", step).table
    assert not [cell for row in table for cell in row if "|" in cell]

    lines = data._md_table(table).splitlines()
    counts = {len(line.split("|")) for line in lines}
    assert len(counts) == 1, f"{step}: header, separator and body disagree on cell count"


def test_that_grid_assertion_FIRES_on_a_pasted_pipe():
    # The companion, because byte-identity does not catch this and neither would an eye.
    broken = (("covariate", "|SMD| < 0.1"), ("age", "yes"))
    lines = data._md_table(broken).splitlines()
    assert len({len(line.split("|")) for line in lines}) == 2


def test_two_runs_render_identical_markdown():
    assert assessed()[2].to_markdown() == assessed()[2].to_markdown()


def test_the_log_is_identical_across_interpreters_with_different_hash_seeds(tmp_path):
    """The two-seed driver of Stage 5 §12.9, extended one stage further.

    Written out rather than sketched, for Stage 3 §12.12's reason: an implementer choosing it freely
    can choose one that fits nothing at all and still see two identical outputs.
    """
    driver = textwrap.dedent("""
        import sys
        sys.path.insert(0, %r)
        sys.path.insert(0, %r)
        import balance, data
        from test_balance import assessed
        _, _, audit, _ = assessed()
        sys.stdout.write(audit.to_markdown())
    """) % (str(MODULE_DIR), str(TESTS_DIR))
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
    assert "- **balance_smd**" in outputs[0] and "- **overlap_by_centre**" in outputs[0]


def test_the_audit_ledger_is_twenty_two_before_and_twenty_four_after_on_the_fixture():
    df, ps, audit = fitted()
    before = len(audit.entries)
    balance.assess(df, ps, audit)
    assert len(audit.entries) - before == 2


@DATA_GATED
def test_the_workbook_ledger_reaches_twenty_four(workbook, workbook_balance):
    _, _, audit = workbook
    # load 7 / derive 4 / classify 1 / build 6 / fit 4 / assess 2
    assert len(audit.entries) == 24


# --- 12.10  Stage 7 adds nothing, refits nothing, and imports no fitter --------------------------------------

def test_the_frame_comes_back_unchanged_cell_for_cell():
    frame, audit = built()
    ps = propensity.fit(frame, audit)
    before = frame.copy(deep=True)
    balance.assess(frame, ps, audit)
    pd.testing.assert_frame_equal(frame, before)
    assert list(frame.columns) == list(before.columns)
    assert frame.index.equals(before.index)
    assert frame.dtypes.equals(before.dtypes)


def test_the_propensity_comes_back_unchanged_and_no_attribute_is_added_to_the_frame():
    frame, audit = built()
    ps = propensity.fit(frame, audit)
    e, w, in_model = ps.e.copy(), ps.w.copy(), ps.in_model.copy()
    columns = list(frame.columns)
    balance.assess(frame, ps, audit)
    pd.testing.assert_series_equal(ps.e, e)
    pd.testing.assert_series_equal(ps.w, w)
    pd.testing.assert_series_equal(ps.in_model, in_model)
    assert list(frame.columns) == columns


def imported_modules(source: str) -> set[str]:
    names: set[str] = set()
    for node in ast.walk(ast.parse(source)):
        if isinstance(node, ast.Import):
            names.update(alias.name.split(".")[0] for alias in node.names)
        elif isinstance(node, ast.ImportFrom) and node.module:
            names.add(node.module.split(".")[0])
    return names


def attribute_calls(source: str) -> set[str]:
    return {f"{node.value.id}.{node.attr}"
            for node in ast.walk(ast.parse(source))
            if isinstance(node, ast.Attribute) and isinstance(node.value, ast.Name)}


@pytest.mark.parametrize("forbidden", ["model", "statsmodels", "sklearn", "scipy"])
def test_balance_imports_no_fitter_and_above_all_NOT_model(forbidden):
    """The `model` half is the one that matters and is the reason this scan exists (§4.2).

    Importing it would not fail. It would silently produce a 16-row table over `BALANCE_SET` — or 11
    over `PS_COVARIATES` — omitting `center = HUG`, which carries this cohort's worst residual
    imbalance, and reporting the worst among the rows it kept as though it were the worst overall.
    """
    assert forbidden not in imported_modules(SOURCE)


def test_that_import_scan_FIRES_on_a_pasted_import():
    assert "model" in imported_modules("import config as C\nimport propensity, model\n")
    assert "model" in imported_modules("from model import design\n")


def test_balance_names_neither_model_firth_nor_propensity_fit():
    calls = attribute_calls(SOURCE)
    assert "model.firth" not in calls
    assert "propensity.fit" not in calls
    assert "propensity.ess" in calls              # the one thing it does call over there
    assert "model.design" not in calls


def test_that_attribute_scan_FIRES():
    assert "propensity.fit" in attribute_calls("ps = propensity.fit(df, audit)\n")


def test_assess_takes_exactly_df_ps_and_audit():
    # No covariate list, no threshold, no centre list (§4.1, §5.4, §15). Asserted by `inspect`, so a
    # defaulted keyword fails this test rather than passing silently: a keyword makes a prespecified
    # choice look like an option.
    parameters = inspect.signature(balance.assess).parameters
    assert list(parameters) == ["df", "ps", "audit"]
    assert all(p.default is inspect.Parameter.empty for p in parameters.values())


def test_smds_signature_is_unchanged_by_the_pooled_sd_split():
    # A `sd=` keyword here would make [§9]'s prescribed denominator look like a caller's option.
    assert list(inspect.signature(balance.smd).parameters) == ["x", "a", "w"]


def test_balance_contains_no_bare_assert():
    offenders = [node.lineno for node in ast.walk(ast.parse(SOURCE))
                 if isinstance(node, ast.Assert)]
    assert offenders == []


def raw_names_in(source: str) -> list[tuple[int, str]]:
    contract = set(config.COLUMN_CONTRACT)
    return [(node.lineno, node.value) for node in ast.walk(ast.parse(source))
            if isinstance(node, ast.Constant) and isinstance(node.value, str)
            and node.value in contract]


def test_neither_new_file_names_a_raw_header_nor_is_exempt_from_the_scan():
    assert raw_names_in(SOURCE) == []
    assert raw_names_in(Path(__file__).read_text(encoding="utf-8")) == []
    assert "balance.py" not in config.EXEMPT_FROM_RAW_NAME_SCAN
    assert "test_balance.py" not in config.EXEMPT_FROM_RAW_NAME_SCAN


def test_that_raw_name_scan_FIRES():
    # The header is ASSEMBLED rather than written, and that is not a flourish: `test_config.py`'s
    # repo-wide scan reads this file too, and a companion holding the literal would make the shipped
    # scan fail on the test that proves it works. Concatenated fragments are the residual hole that
    # scan's own docstring declares, which is what makes this safe to rely on here.
    header = "Prestroke" + "mRS "
    assert raw_names_in(f"x = frame[{header!r}]\n") == [(1, header)]


def test_no_tail_slice_of_the_entries_appears_in_this_file():
    """Stage 5 §10's repair. `test_propensity.py` carries the same scan over its own two files and
    stays unedited (DoD-3), so this module carries its own.

    An assertion on `audit.entries[-2:]` passes whether or not the entries it names are the ones the
    call added, because any two trailing entries satisfy it.
    """
    source = Path(__file__).read_text(encoding="utf-8")
    offenders = [node.lineno for node in ast.walk(ast.parse(source))
                 if isinstance(node, ast.Subscript) and isinstance(node.slice, ast.Slice)
                 and isinstance(node.slice.lower, ast.UnaryOp)
                 and isinstance(node.slice.lower.op, ast.USub)]
    assert offenders == []


# --- 12.11  the preconditions ----------------------------------------------------------------------------

def shifted(series: pd.Series) -> pd.Series:
    """The same values on an index sharing no label with the frame's — equal LENGTH, wrong rows."""
    return pd.Series(series.to_numpy(), index=series.index + 1000, name=series.name)


def misaligned_propensity(ps: propensity.Propensity, which: str) -> propensity.Propensity:
    return propensity.Propensity(
        e=shifted(ps.e) if which == "e" else ps.e,
        w=shifted(ps.w) if which == "w" else ps.w,
        in_model=shifted(ps.in_model) if which == "in_model" else ps.in_model,
        ess=ps.ess, fit=ps.fit, dropped=ps.dropped)


@pytest.mark.parametrize("which", ["e", "w", "in_model"])
def test_B1_fires_on_each_of_the_three_series_and_NAMES_which(which):
    """The `e`-alone case is the one an earlier draft missed: it reached `_centre`'s
    `ps.e.loc[weighted]` and raised pandas' `IndexingError` AFTER `balance_smd` had been recorded —
    a log describing half a diagnostic, from an input B1 was meant to reject.
    """
    df, ps, audit = fitted()
    before = len(audit.entries)
    with pytest.raises(config.SchemaError) as excinfo:
        balance.assess(df, misaligned_propensity(ps, which), audit)

    message = message_of(excinfo)
    assert "B1" in message and f"{which} carries a different index" in message
    assert "index label(s)" in message                    # not the two lengths, which are both 5
    assert len(audit.entries) == before                   # and NOTHING was recorded


def test_B1s_message_reports_the_shared_label_count_because_the_lengths_are_equal():
    df, ps, _ = fitted()
    bad = misaligned_propensity(ps, "in_model")
    assert len(bad.in_model) == len(df)                   # a message quoting lengths says nothing
    with pytest.raises(config.SchemaError) as excinfo:
        balance._assert_balance_inputs(df, bad)
    message = message_of(excinfo)
    assert f"{len(bad.in_model)} mask row(s) against {len(df)} frame row(s)" in message
    assert "sharing 0 index label(s)" in message          # the number that says something


def test_WITHOUT_the_phase_one_raise_the_same_input_dies_with_a_BARE_AssertionError():
    """§4.5a, measured. B1's paragraph is assembled into `bad` and then discarded by an exception
    pandas raises two checks later — with the empty string as its message, on the one failure mode
    B1 exists to describe.
    """
    df, ps, _ = fitted()
    bad = misaligned_propensity(ps, "in_model")
    with pytest.raises(AssertionError) as excinfo:        # what B4/B5 would hit
        df.loc[bad.in_model, config.TREATMENT]
    assert message_of(excinfo) == ""                      # no name, no explanation, nothing


def test_WITHOUT_B1_ENTIRELY_the_positional_read_produces_a_COMPLETE_FINITE_table():
    """The failure mode B1 exists for, and the reason it is an index check and not a length check.

    Two Series of equal length on different indexes align to nothing under `.loc` and to the WRONG
    ROWS under `.to_numpy()` — and the second is what this module does. Nothing raises and every
    cell is finite: a balance table for a population that does not exist.
    """
    df, ps, _ = fitted()
    sub = df.loc[ps.in_model]
    a = sub[config.TREATMENT].to_numpy(dtype=float)
    wrong = ps.w.loc[ps.in_model].to_numpy(dtype=float)[::-1]      # the same values, wrong records
    rows = balance._table(sub, a, wrong, config.BALANCE_SET)
    assert len(rows) == declared_level_count()
    assert all(np.isfinite(row.weighted) for row in rows)
    # and it disagrees with the truth, silently
    right = {row.covariate: row.weighted for row in balance._table(
        sub, a, ps.w.loc[ps.in_model].to_numpy(dtype=float), config.BALANCE_SET)}
    assert any(row.weighted != right[row.covariate] for row in rows)


def three_valued_mask(ps: propensity.Propensity) -> propensity.Propensity:
    mask = ps.in_model.astype("object")
    mask.iloc[0] = pd.NA
    return propensity.Propensity(e=ps.e, w=ps.w, in_model=mask, ess=ps.ess, fit=ps.fit,
                                 dropped=ps.dropped)


def test_B2_fires_on_a_three_valued_mask_and_records_nothing():
    df, ps, audit = fitted()
    before = len(audit.entries)
    with pytest.raises(config.SchemaError) as excinfo:
        balance.assess(df, three_valued_mask(ps), audit)
    assert "B2" in message_of(excinfo)
    assert len(audit.entries) == before


def test_WITHOUT_the_phase_one_raise_a_three_valued_mask_dies_in_pandas_instead():
    df, ps, _ = fitted()
    with pytest.raises(ValueError) as excinfo:
        df.loc[three_valued_mask(ps).in_model, config.TREATMENT]
    assert "Cannot mask with non-boolean array containing NA / NaN values" in message_of(excinfo)


@pytest.mark.parametrize("column", config.BALANCE_SET)
def test_B3_fires_on_a_frame_missing_any_balance_set_column(column):
    df, ps, audit = fitted()
    before = len(audit.entries)
    with pytest.raises(config.SchemaError) as excinfo:
        balance.assess(df.drop(columns=[column]), ps, audit)
    assert "B3" in message_of(excinfo) and column in message_of(excinfo)
    assert len(audit.entries) == before


def test_B3_fires_on_a_frame_missing_the_TREATMENT_column_TOO():
    # `ivt` is in neither BALANCE_SET nor CENTER_ORDER and is read by B5, by §8's `a`, and by every
    # arm mask in §5 and §6.
    df, ps, audit = fitted()
    with pytest.raises(config.SchemaError) as excinfo:
        balance.assess(df.drop(columns=[config.TREATMENT]), ps, audit)
    assert "B3" in message_of(excinfo) and config.TREATMENT in message_of(excinfo)


def test_WITH_B3_IN_PHASE_TWO_the_TREATMENT_case_would_die_at_B5_with_a_KeyError():
    """§21.1's finding: §20's E5 named `TREATMENT` in B3 and left B5 reading it two lines later, so
    the same defect survived the repair one check further along. This companion is what would have
    caught it.
    """
    df, ps, _ = fitted()
    with pytest.raises(KeyError) as excinfo:
        df.drop(columns=[config.TREATMENT]).loc[ps.in_model, config.TREATMENT]
    assert config.TREATMENT in str(excinfo.value)


def test_a_frame_failing_B1_and_B4_reports_only_B1_and_says_B4_was_not_run(monkeypatch):
    """The alternative reading — two of five assertions named on a frame that fails three — is a
    clean bill of health for the ones nobody checked (§4.5a).
    """
    df, ps, _ = fitted()
    monkeypatch.setattr(config, "BALANCE_SET", config.BALANCE_SET + ("onset_to_groin_min",))
    with pytest.raises(config.SchemaError) as excinfo:
        balance._assert_balance_inputs(df, misaligned_propensity(ps, "in_model"))
    message = message_of(excinfo)
    # "B4  " with the check's own two-space prefix — the summary line says "B4 and B5 were not run",
    # which is the sentence being asserted rather than a fired assertion.
    assert "B1  " in message and "B4  " not in message
    assert "B4 and B5 were not run" in message


def test_B4_fires_on_a_post_time_zero_member_of_the_balance_set(monkeypatch):
    # Unreachable on the declared set — `test_config.py` asserts BALANCE_SET and POST_TIME_ZERO are
    # disjoint — and recorded rather than treated as a reason to omit the check, exactly as Stage 5's
    # C1-C4 and Stage 6's D1-D4 are.
    assert set(config.BALANCE_SET) & config.POST_TIME_ZERO == set()

    df, ps, audit = fitted()
    before = len(audit.entries)
    monkeypatch.setattr(config, "BALANCE_SET", ("age", "mrs_90d"))
    with pytest.raises(config.SchemaError) as excinfo:
        balance.assess(df, ps, audit)
    assert "B4" in message_of(excinfo) and "mrs_90d" in message_of(excinfo)
    assert len(audit.entries) == before


def test_B5_fires_when_in_model_leaves_one_arm_and_names_both_arms_and_their_counts():
    df, ps, audit = fitted()
    before = len(audit.entries)
    # `.astype(bool)` because `ps.in_model & (df[TREATMENT] == 1)` over an Int64 column is nullable
    # `boolean`, which B2 catches first — a different assertion from the one under test.
    single = propensity.Propensity(
        e=ps.e, w=ps.w,
        in_model=(ps.in_model & (df[config.TREATMENT] == 1)).astype(bool),
        ess=ps.ess, fit=ps.fit, dropped=ps.dropped)
    with pytest.raises(config.SchemaError) as excinfo:
        balance.assess(df, single, audit)
    message = message_of(excinfo)
    assert "B5" in message
    for label in config.TREATMENT_LABELS.values():
        assert label in message
    assert "an arm missing HERE was lost by in_model" in message
    assert len(audit.entries) == before


def test_several_failures_WITHIN_a_phase_produce_one_error_naming_each():
    df, ps, audit = fitted()
    # phase 1: B1 and B2 together
    both = three_valued_mask(misaligned_propensity(ps, "e"))
    with pytest.raises(config.SchemaError) as excinfo:
        balance.assess(df, both, audit)
    assert "B1" in message_of(excinfo) and "B2" in message_of(excinfo)
    assert "2 balance assertion(s) failed" in message_of(excinfo)


def test_phase_two_collects_B4_and_B5_together(monkeypatch):
    df, ps, audit = fitted()
    monkeypatch.setattr(config, "BALANCE_SET", ("age", "mrs_90d"))
    # `.astype(bool)` because `ps.in_model & (df[TREATMENT] == 1)` over an Int64 column is nullable
    # `boolean`, which B2 catches first — a different assertion from the one under test.
    single = propensity.Propensity(
        e=ps.e, w=ps.w,
        in_model=(ps.in_model & (df[config.TREATMENT] == 1)).astype(bool),
        ess=ps.ess, fit=ps.fit, dropped=ps.dropped)
    with pytest.raises(config.SchemaError) as excinfo:
        balance.assess(df, single, audit)
    message = message_of(excinfo)
    assert "B4" in message and "B5" in message
    assert "2 balance assertion(s) failed" in message
    assert "of them weighted" in message                  # phase 2's summary line, not phase 1's


# --- 12.12  structural facts from the workbook [data-gated] ------------------------------------------------

@DATA_GATED
def test_the_workbook_table_is_nineteen_rows_every_n_is_92_and_none_is_undefined(
        workbook, workbook_balance):
    _, ps, _ = workbook
    assert len(workbook_balance.covariates) == 19
    assert {row.n for row in workbook_balance.covariates} == {int(ps.in_model.sum())} == {92}
    assert workbook_balance.worst()[1] == ()


@DATA_GATED
def test_the_worst_before_is_center_HUG_and_the_worst_after_is_a_NEGATIVE_CONTROL(workbook_balance):
    # Asserted to two decimals, so a change in the estimator moves them.
    before = max(workbook_balance.covariates, key=lambda row: abs(row.unweighted))
    after, undefined = workbook_balance.worst()
    assert before.covariate == "center = HUG"
    assert abs(before.unweighted) == pytest.approx(1.34, abs=5e-3)
    assert after.covariate == "diabetes" and after.role == "negative control"
    assert abs(after.weighted) == pytest.approx(0.38, abs=5e-3)
    assert undefined == ()


@DATA_GATED
def test_exactly_five_rows_reach_the_threshold_and_exactly_TWO_are_in_the_propensity_model(
        workbook_balance):
    """§6.4's finding, pinned so a future workbook or specification cannot quietly change it.

    Two of the five are IN the propensity model, and that is the finding: overlap weights have an
    exact-balance property and the [§7] Firth estimator does not deliver it. It is reported here and
    handed to the PI under [§9] and [§13]; it is not a licence to change the estimator.
    """
    unbalanced = workbook_balance.unbalanced()
    assert unbalanced == ("center = HUG", "center = Lugano", "hypertension",
                          "hyperlipidemia", "diabetes")
    in_model = [row.covariate for row in workbook_balance.covariates
                if row.covariate in unbalanced and row.role == "propensity model"]
    assert in_model == ["center = HUG", "center = Lugano"]
    assert len(unbalanced) == 5 and len(in_model) == 2

    # eleven of nineteen were above the threshold BEFORE weighting
    assert sum(1 for row in workbook_balance.covariates
               if abs(row.unweighted) >= config.SMD_THRESHOLD) == 11


@DATA_GATED
def test_the_in_model_worst_is_above_the_threshold_at_0_211_on_center_HUG(workbook_balance):
    in_model = [row for row in workbook_balance.covariates if row.role == "propensity model"]
    worst = max(in_model, key=lambda row: abs(row.weighted))
    assert worst.covariate == "center = HUG"
    assert worst.weighted == pytest.approx(0.211121, abs=1e-5)
    assert abs(worst.weighted) > config.SMD_THRESHOLD          # [§9]'s threshold, exceeded


@DATA_GATED
def test_center_USZ_reads_zero_in_both_columns_and_its_centre_row_is_structural(workbook_balance):
    row = next(r for r in workbook_balance.covariates if r.covariate == "center = USZ")
    assert (row.unweighted, row.weighted) == (0.0, 0.0)
    usz = next(r for r in workbook_balance.centres if r.centre == "USZ")
    assert usz.status == balance._STRUCTURAL


WORKBOOK_CENTRES = {
    #          ESS control  ESS treated   max w      weight share
    "HUG":    (10.6761,     25.3441,      0.831072,  0.566719),
    "CHUV":   (12.2544,     6.8295,       0.630395,  0.275094),
    "Lugano": (19.8179,     1.9765,       0.964955,  0.158186),
}


@DATA_GATED
@pytest.mark.parametrize("centre", sorted(WORKBOOK_CENTRES))
def test_the_three_reported_centres_carry_the_measured_overlap_statistics(centre, workbook_balance):
    row = next(r for r in workbook_balance.centres if r.centre == centre)
    ess_control, ess_treated, max_w, share = WORKBOOK_CENTRES[centre]
    assert row.status == balance._REPORTED
    assert row.ess[0] == pytest.approx(ess_control, abs=1e-4)
    assert row.ess[1] == pytest.approx(ess_treated, abs=1e-4)
    assert row.max_weight == pytest.approx(max_w, abs=1e-4)
    assert row.weight_share == pytest.approx(share, abs=1e-4)


@DATA_GATED
def test_the_defined_weight_shares_sum_to_one_with_USZs_reading_missing(workbook, workbook_balance):
    defined = [row.weight_share for row in workbook_balance.centres
               if np.isfinite(row.weight_share)]
    assert len(defined) == 3
    assert sum(defined) == pytest.approx(1.0, abs=1e-12)

    _, _, audit = workbook
    assert rows_of(audit.entry("model", "overlap_by_centre"))["USZ"][-3] == "missing"


@DATA_GATED
def test_statsmodels_is_a_declared_dependency_so_the_companion_below_needs_no_importorskip():
    # Checked rather than assumed: `statsmodels` is in `pyproject.toml`'s main dependency list
    # ("retained for unpenalised cross-checks in tests") and not in the gated `reference` group, so
    # the companion below runs on any `uv sync`. Stage 6's oracle is the one that needs the guard,
    # and its package is deliberately not named here — `test_model.py` owns that name and a scan
    # over there asserts it appears nowhere else.
    text = (MODULE_DIR / "pyproject.toml").read_text(encoding="utf-8")
    main, groups = text.split("[dependency-groups]", 1)
    assert "statsmodels" in main and "statsmodels" not in groups
    assert "reference = [" in groups             # the gated group exists and is not where we are


@DATA_GATED
def test_an_unpenalised_MLE_score_balances_every_in_model_row_and_the_FIRTH_score_does_not(workbook):
    """§6.4's attribution, and the pair is the assertion.

    Overlap weights solve the balance equation for every covariate in the score BECAUSE THAT IS THE
    SCORE EQUATION; Firth's modified score is `X'(y − p + h(0.5 − p)) = 0`, a different equation, so
    the property holds only approximately. The first half alone would pass on a silently reinstated
    MLE — the substitution Stage 6 §12.10 exists to catch — and the second half is what makes this a
    test about Firth rather than about arithmetic.

    The negative controls are imbalanced under BOTH scores, because exact balance is a property of
    the covariates in the score and of no others. That is what [§6] retains them for.

    This is an ORACLE, run in a test. Swapping it into the pipeline would exchange a documented
    residual imbalance for the mixture-of-two-estimators bootstrap [§7] was written against.
    """
    import statsmodels.api as sm

    df, ps, _ = workbook
    X, _ = model.design(df.loc[ps.in_model], config.PS_COVARIATES)
    a = df.loc[ps.in_model, config.TREATMENT].to_numpy(dtype=float)
    mle = sm.Logit(a, sm.add_constant(X.to_numpy(dtype=float))).fit(disp=0)

    e = pd.Series(np.nan, index=df.index, dtype="float64")
    e.loc[ps.in_model] = np.asarray(mle.predict(), dtype=float)
    w = pd.Series(np.nan, index=df.index, dtype="float64")
    w.loc[ps.in_model] = np.where(a == 1.0, 1.0 - e.loc[ps.in_model], e.loc[ps.in_model])
    unpenalised = propensity.Propensity(e=e, w=w, in_model=ps.in_model, ess=ps.ess, fit=ps.fit,
                                        dropped=ps.dropped)

    # a FRESH Audit: this companion is an oracle and must not append to the pipeline's log
    mle_balance = balance.assess(df, unpenalised, data.Audit(data.WORKBOOK))
    in_model_rows = [row for row in mle_balance.covariates if row.role == "propensity model"]
    controls = [row for row in mle_balance.covariates if row.role == "negative control"]

    assert max(abs(row.weighted) for row in in_model_rows) < 1e-12       # exactly balanced
    assert max(abs(row.weighted) for row in controls) > 0.1             # and never were

    # and under the [§7] Firth score the same quantity is above 0.2
    firth_rows = [row for row in balance.assess(df, ps, data.Audit(data.WORKBOOK)).covariates
                  if row.role == "propensity model"]
    assert max(abs(row.weighted) for row in firth_rows) > 0.2


# --- 12.9a  the two `detail` strings, read as PROSE under every one of their branches ------------------
#
# §21.5 named this as where a third review round should start: two rounds rendered the tables and read
# their cells, and neither checked that a `detail` sentence's counts match the table beside it or that
# its conditional clauses read as English at 0, 1 and many. Reading them found one defect — the
# undefined-rows clause ran into the next sentence for want of a full stop — which is why these
# assertions are on the SENTENCE and not only on the count.

def sentences_of(detail: str) -> list[str]:
    return [s.strip() for s in detail.split(". ") if s.strip()]


def test_the_smd_detail_reads_as_english_when_NO_row_is_undefined():
    _, _, audit, bal = assessed()
    detail = audit.entry("model", "balance_smd").detail
    assert bal.worst()[1] == ()
    assert ". No row is undefined. A row's `n` is its own denominator" in detail
    assert f"{len(bal.unbalanced())} row(s) reach |SMD| = 0.1 after weighting: " in detail
    assert f"worst {data._fmt(abs(bal.worst()[0].weighted))} on {bal.worst()[0].covariate}" in detail


def test_the_smd_detail_reads_as_english_when_EVERY_row_is_undefined():
    df, ps, audit = one_armed_centre()
    bal = balance.assess(df, ps, audit)
    detail = audit.entry("model", "balance_smd").detail

    assert len(bal.worst()[1]) == len(bal.covariates) and bal.unbalanced() == ()
    assert "0 row(s) reach |SMD| = 0.1 after weighting. " in detail        # no colon, no list
    assert f"{len(bal.worst()[1])} row(s) are undefined" in detail
    # the clause ENDS, rather than running into the sentence after it
    assert "penumbra_ml. A row's `n` is its own denominator" in detail
    assert not [s for s in sentences_of(detail) if s.startswith("A row's") and "penumbra" in s]


def test_the_overlap_detail_reads_as_english_at_zero_one_and_many():
    # one structural centre and no thin one, on the workbook's shape; two structural on the fixture
    _, _, audit, bal = assessed()
    detail = audit.entry("model", "overlap_by_centre").detail
    structural = [row.centre for row in bal.centres if row.status == balance._STRUCTURAL]
    assert f"{len(structural)} centre(s) contribute no cohort record" in detail
    assert ": " + ", ".join(structural) + "." in detail
    assert "lost an arm to complete-casing" not in detail          # none here, so no clause

    # and with a thin centre, the clause appears and closes
    df, ps, audit_thin = one_armed_centre()
    balance.assess(df, ps, audit_thin)
    thin_detail = audit_thin.entry("model", "overlap_by_centre").detail
    assert "1 centre(s) lost an arm to complete-casing and carry no weighted contrast: CHUV." \
        in thin_detail
    assert "CHUV. `restrict_centres` above" in thin_detail


def test_every_count_in_both_details_is_interpolated_and_never_written_out():
    # §21.4's finding 7: `_smd_detail` hard-coded "the four vascular risk factors" against §4.3's
    # computed-set promise and against the very patch §12.1 applies.
    source = MODULE.read_text(encoding="utf-8")
    detail_source = source[source.index("def _smd_detail"):source.index("def _overlap_detail")]
    for written_out in ("four vascular", "19 declared", "14 [§6]", "0.10 after"):
        assert written_out not in detail_source
    assert "len(C.NEGATIVE_CONTROLS)" in detail_source
    assert "len(C.BALANCE_SET)" in detail_source
