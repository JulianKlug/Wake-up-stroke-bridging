"""Acceptance tests for Stage 8b — §14 of `specs/stage8_primary_outcome_estimator.md`.

`polr`'s own sections — §14.3, §14.4, §14.7's fit half and §14.8 — live in `test_model.py`, because
that is where the function is. Everything here reads a `Primary` or an `Audit`.

Tests needing the private workbook are marked `DATA_GATED` and take this file's **own** module-scoped
`workbook` fixture. Everything else runs on a plain checkout with no `data/`.

**§14.12 ASSERTS PROPERTIES AND NOT VALUES, AND §4.3 IS WHY.** Stage 7 pinned 1.343 and 0.380 into
`test_balance.py` and was right to: a balance diagnostic is evidence *about* the analysis. Stage 8's
output IS the analysis, and every stage of this pipeline records its decisions as taken before the
outcome was examined by arm. So `beta`, `exp(beta)` and the six `RD_k` **on the workbook** are not in
this file, are not in any file under `specs/`, and are produced at run time into the gitignored log.
The regression pin is §14.0.2's golden vector instead, from a fixture-derived cohort with no patient
data in it — and a meta-assertion at the foot of §14.12 scans that section's own source for a float
literal outside a declared allowlist, which is what keeps this true as the file is edited.

This file is **not** exempt from the Stage 1 §7 raw-name scan and must not become exempt.

**The three silent failures this stage exists to make loud**, each tested with a companion showing
what happens *without* the guard::

    a separated fit        G7   without it: exp(beta) = 6.5e15, converged in 17 iterations   §14.7
                                on the score criterion, every counter 0, nothing raised
    a vanished exposure    G6   without it: a width-0 design FITS, and the caller gets       §14.11
                                either IndexError on beta[0] or an odds ratio of exactly 1.0
    a misaligned score     G1   without it: odds_ratio 58.01 against the correct 5.76,       §14.14
                                because w arrives in the Propensity's order
"""
from __future__ import annotations

import ast
import inspect
import re
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
import outcome
import propensity
from test_cohort import built, cohort_frame, set_cell       # tests/test_cohort.py:116, :147, :152
from test_data import hand_source, run                      # tests/test_data.py:98, :107
from test_model import (deaths_both_arms_frame, hand_ordinal, orientation_frame,
                        rare_category_frame, reference_frame, separated_frame)

DATA_GATED = pytest.mark.skipif(
    not config.DATA_XLSX.exists(),
    reason="the private workbook is gitignored and absent from this checkout")

MODULE = Path(outcome.__file__).resolve()
MODULE_DIR = MODULE.parent
TESTS_DIR = Path(__file__).resolve().parent
SOURCE = MODULE.read_text(encoding="utf-8")
THIS_FILE = Path(__file__).resolve().read_text(encoding="utf-8")


def message_of(excinfo) -> str:
    return str(excinfo.value)


def rows_of(entry: data.AuditEntry) -> dict[str, tuple[str, ...]]:
    """A table's body keyed by its first cell, so a test names the row it is asserting."""
    return {row[0]: row[1:] for row in entry.table[1:]}


# --- 14.0  the frames and fixtures this stage is tested on --------------------------------------------
#
# THE FIRST OF THEM IS A PROBLEM NO EARLIER STAGE HAD. The committed fixture cohort cannot reach this
# stage at all: `cohort_frame()` builds its three cohort records from the `HAND-4` template, so every
# analysis column not overridden is constant across the frame — and `mrs_90d` is one of them, at 2.
# `polr` on it raises O5 with one weighted category. Stage 6 fitted a propensity model on that frame
# because the EXPOSURE varies there by construction; Stage 8 is the first stage whose RESPONSE is a
# column the fixture holds constant.
#
# So this stage declares its frames as OVERRIDES on the committed fixture rather than as new frames,
# for Stage 3 §12's reason: a second hand-built cohort would be a second place the pipeline's fixture
# data lives, and the two would diverge on the next Stage 1 amendment.


def ordinal_cohort() -> pd.DataFrame:
    """`cohort_frame()` with mrs_90d spread across the declared levels. §14.0.

    The spread is a `pd.Series` over the frame's index and covers 0 through 6 in order, so which of
    the nine records survive both restrictions decides which levels the cohort carries — measured: the
    five survivors take mrs_90d 0, 1, 4, 6 and 3, so the cohort is missing levels 2 and 5 and §5.3's
    collapse fires on it. That is the point: the collapse branch has a real-frame witness that runs
    with no `data/`, and the workbook has all seven levels and cannot provide one.

    The `[: len(base)]` slice is a ONE-DIRECTIONAL guard and stays one, which is a decision rather
    than an oversight: it is total if `cohort_frame()` ever shrinks and raises a bare pandas length
    mismatch if it GROWS. Cycling MRS_LEVELS over the index instead would be total in both directions
    — and would give a DIFFERENT spread in the surviving positions, `[0, 1, 4, 6, 1]` against
    `[0, 1, 4, 6, 3]`, so `categories` would become `(0, 1, 4, 6)` and every number in §14.0.2's
    golden vector would move. A totality repair that changes the pin the fixture exists to carry is
    not a repair. The hazard is recorded rather than removed, and the survivors assertion below is
    what names the cause when the mismatch fires.
    """
    base = cohort_frame()
    spread = pd.Series([0, 1, 2, 3, 4, 5, 6, 2, 3][: len(base)],
                       index=base.index, dtype="Int64")
    return cohort_frame(mrs_90d=spread)


def truncated_cohort() -> pd.DataFrame:
    """`ordinal_cohort()` with the worst mRS removed, so the top FITTED level is not the top DECLARED
    one. §9.2, §14.9.

    `ordinal_cohort()`'s survivors take mrs_90d 0, 1, 4, 6 and 3, so its top fitted category is 6 —
    which is `MRS_LEVELS[-1]`, by coincidence and not by design. The workbook occupies all seven
    levels, so its top fitted category is 6 as well. **Every other frame in this suite therefore has
    `fit.categories[-1] == MRS_LEVELS[-1]`, and §9.2's status branch is untested without this one.**

    Replacing the 6 with a 2 gives categories (0, 1, 2, 3, 4): five occupied levels, four cutpoints,
    mRS 5 and 6 both absent, and the top fitted level is 4. That is the frame on which the draft §9.2
    corrects rendered `alpha at mRS <= 4` as `missing` under status `fitted`, and `alpha at mRS <= 6`
    as the highest declared level with no mention that nobody is in it.

    Which element to change is derived from a measurement rather than guessed: the survivors carry
    `[0, 1, 4, 6, 3]` against the spread `[0, 1, 2, 3, 4, 5, 6, 2, 3]`, so the surviving index
    positions are 0, 1, 4, 6 and 8 — the seventh element is one of them and the sixth is not. Which
    records survive does not depend on `mrs_90d` at all: the [§3] restrictions are on centre, arm and
    eligibility, so the override moves the level set and never the population.

    A [§10] replicate that draws no death produces exactly this shape, which is why the branch is
    worth a fixture rather than a note.
    """
    base = cohort_frame()
    spread = pd.Series([0, 1, 2, 3, 4, 5, 2, 2, 3][: len(base)],
                       index=base.index, dtype="Int64")
    return cohort_frame(mrs_90d=spread)


def separated_cohort() -> tuple[pd.DataFrame, propensity.Propensity, data.Audit]:
    """`ordinal_cohort()`'s cohort with the outcome made to separate the arms perfectly. §14.7.

    The Propensity is the REAL one, fitted before the outcome is overwritten, so the weights are a
    genuine [§7] score and only the response is constructed — which is what makes `primary` raise G7
    for the reason §6 gives rather than because the frame is unreachable.

    It is `test_outcome.py`'s and not `test_model.py`'s because "`primary` on the same frame raises
    G7" needs a cohort frame and a `Propensity`, which `separated_frame()`'s three arrays are not.
    """
    df, audit = built(ordinal_cohort())
    ps = propensity.fit(df, audit)
    df = df.copy()
    df.loc[df[config.TREATMENT] == 1, config.PRIMARY_OUTCOME] = config.MRS_LEVELS[0]
    df.loc[df[config.TREATMENT] == 0, config.PRIMARY_OUTCOME] = config.MRS_LEVELS[-1]
    return df, ps, audit


def estimated(frame=None):
    """(df, ps, est, audit) — the whole stage run once on a fixture-derived cohort."""
    df, audit = built(ordinal_cohort() if frame is None else frame)
    ps = propensity.fit(df, audit)
    return df, ps, outcome.primary(df, ps, audit), audit


FIXTURE_COHORTS = {"ordinal_cohort": ordinal_cohort, "truncated_cohort": truncated_cohort}


@pytest.fixture(scope="module")
def workbook():
    """(df, ps, audit) from ONE linear run against ONE Audit. §14.0.3.

    `test_propensity.py`'s `workbook_ps` deliberately does the opposite, and Stage 7 §12.0.3 records
    why. The shape is specified rather than left to the pattern because §14.10 reads
    `overlap_weights` out of the audit and compares it with Stage 8's own entries, so both must be in
    the SAME `Audit`. `balance.assess` is not called by it, so the ledger reaching 27 is the
    full-pipeline count of §11 and not this fixture's.
    """
    df, audit = data.load(data.WORKBOOK)
    df = derive.derive(df, audit)
    df = eligibility.classify(df, audit)
    df = cohort.build(df, audit)
    return df, propensity.fit(df, audit), audit


@pytest.fixture(scope="module")
def workbook_primary(workbook):
    """`primary` run ONCE over the fixture above, so the reading tests share one call."""
    df, ps, audit = workbook
    return outcome.primary(df, ps, audit)


# --- 14.0.1  what `ordinal_cohort()` exercises ---------------------------------------------------------

def test_the_COMMITTED_fixture_cohort_cannot_reach_this_stage_at_all():
    """Measured, and it is why §14.0 declares an override at all.

    `cohort_frame()`'s five-record cohort has `mrs_90d` [2, 2, 2, 2, 2] — one category — so `polr`
    raises O5 and no fit exists. A test that meant to reach `primary` and passed the committed frame
    would see a passing `pytest.raises` whose message is never read.
    """
    df, audit = built()
    assert list(df[config.PRIMARY_OUTCOME]) == [2, 2, 2, 2, 2]
    ps = propensity.fit(df, audit)
    with pytest.raises(model.FitError) as e:
        outcome.primary(df, ps, audit)
    assert "O5" in message_of(e) and "1 response category" in message_of(e)


def test_ordinal_cohort_builds_the_same_five_records_and_spreads_the_outcome():
    """If the override changed which records survive, every number in §14.0.2 would be describing a
    different population — so the survivors are asserted against the five `case_id`s
    `test_cohort.py` already pins, and the assertion is what names the cause if the fixture moves.
    """
    df, _ = built(ordinal_cohort())                      # asserted to raise nothing
    assert set(df["case_id"]) == {"HAND-1", "HAND-2", "HAND-5", "COHORT-1", "COHORT-3"}
    assert list(df[config.PRIMARY_OUTCOME]) == [0, 1, 4, 6, 3]
    assert int((df[config.TREATMENT] == 1).sum()) == 3
    assert int((df[config.TREATMENT] == 0).sum()) == 2


def test_the_derived_dichotomies_come_back_RECOMPUTED_from_the_ordinal_source():
    # Stage 3's contract, and this override is the first thing in the repository to exercise it on a
    # VARYING mrs_90d. A carried column would give the template's constant 2 everywhere.
    df, _ = built(ordinal_cohort())
    assert list(df["mrs_0_2_90d"]) == [1, 1, 0, 0, 0]
    assert list(df["mrs_0_1_90d"]) == [1, 1, 0, 0, 0]
    assert list(df["death_90d"]) == [0, 0, 0, 1, 0]
    assert list(df["mrs_5_6_90d"]) == [0, 0, 0, 1, 0]


def test_the_BARE_versus_NORMALISED_rule_still_applies():
    # Stage 5 §12.0.1's distinction carries forward: `cohort_frame()` returns raw centre codes and
    # `run(...)` maps them. A test that means to reach `primary` and passes a BARE frame never gets
    # there — Stage 5's C4 raises first — and the failure looks like a passing `pytest.raises` whose
    # message is never read.
    with pytest.raises(config.SchemaError) as e:
        cohort.build(ordinal_cohort(), data.Audit(hand_source(Path("."), "bare")))
    assert "C4" in message_of(e) or "C1" in message_of(e)


# --- 14.0.2  the golden vector -------------------------------------------------------------------------
#
# Pinned as literals, from `ordinal_cohort()` so that NO PATIENT-DERIVED NUMBER ENTERS GIT (Stage 6
# §7.3) — and, at this stage, so that no ESTIMATE enters it either (§4.3). Asserted to 1e-6, for
# Stage 6 §12.0.2's and Stage 7 §12.0.2's reason: this is a regression pin against an edit to our own
# arithmetic, and pinning it at machine precision would make it fail on a numpy patch release rather
# than on a mistake.
#
# THE ASSERTION IS BY NAME AND BY THRESHOLD, ONE AT A TIME — alpha against cut_levels and RD_k against
# k — because the alpha vector is monotone and a set comparison would pass while permuting it.

GOLDEN_ALPHA = {0: -4.0128444951, 1: -1.8393869358, 3: -0.7580715177, 4: 0.4483367817}
GOLDEN_RD = {0: +0.1079772317, 1: +0.5469186443, 2: +0.5469186443,
             3: +0.0181363724, 4: +0.4712177281, 5: +0.4712177281}
GOLDEN_CUMULATIVE = {
    0: {0: 0.0000000000, 1: 0.1079772317},
    1: {0: 0.0000000000, 1: 0.5469186443},
    2: {0: 0.0000000000, 1: 0.5469186443},
    3: {0: 0.5287822719, 1: 0.5469186443},
    4: {0: 0.5287822719, 1: 1.0000000000},
    5: {0: 0.5287822719, 1: 1.0000000000},
}


def test_the_golden_vector_reproduces_the_fit():
    _, _, est, _ = estimated()
    assert est.fit.categories == (0, 1, 3, 4, 6)
    assert est.cut_levels == (0, 1, 3, 4)
    assert est.beta == pytest.approx(1.7516026138, abs=1e-6)
    assert est.odds_ratio == pytest.approx(5.7638324754, abs=1e-6)
    assert (est.fit.iterations, est.fit.converged_on) == (4, "likelihood")
    assert (est.fit.rescales, est.fit.halvings) == (0, 0)


@pytest.mark.parametrize("level", sorted(GOLDEN_ALPHA))
def test_the_golden_vector_reproduces_each_cutpoint_BY_ITS_DECLARED_LEVEL(level):
    """`cut_levels` is `(0, 1, 3, 4)` and NOT `(0, 1, 2, 3)`: the third cutpoint sits above mRS 3
    because level 2 is absent, and labelling it `MRS_THRESHOLDS[2]` would name a cutpoint with a
    threshold it is not."""
    _, _, est, _ = estimated()
    assert est.alpha[est.cut_levels.index(level)] == pytest.approx(GOLDEN_ALPHA[level], abs=1e-6)


@pytest.mark.parametrize("k", sorted(GOLDEN_RD))
def test_the_golden_vector_reproduces_each_RD_and_each_cumulative_probability(k):
    _, _, est, _ = estimated()
    assert est.rd[k] == pytest.approx(GOLDEN_RD[k], abs=1e-6)
    for code in config.TREATMENT_LABELS:
        assert est.cumulative[k][code] == pytest.approx(GOLDEN_CUMULATIVE[k][code], abs=1e-6)


def test_the_golden_vectors_three_LOAD_BEARING_features_each_assert_separately():
    """So that a failure names its own cause rather than reading as "the pin moved"."""
    _, _, est, _ = estimated()
    # k = 1 and k = 2 are identical, and so are k = 4 and k = 5: levels 2 and 5 are unoccupied, so
    # the cumulative probability cannot move across them. NON-DECREASING, never strictly increasing.
    assert est.rd[1] == est.rd[2] and est.rd[4] == est.rd[5]
    assert 2 not in est.fit.categories and 5 not in est.fit.categories
    # the third cutpoint is labelled by the level it cuts at or below, not by a threshold index
    assert est.cut_levels == est.fit.categories[:-1]
    # P(Y <= 3) in the control arm is one control record's whole weight arriving at once on a
    # two-record arm, which is what makes this vector reproducible by hand — and a weak witness for
    # anything about the workbook (§4.3, §15).
    assert est.cumulative[2][0] == 0.0 and est.cumulative[3][0] > 0.5


# --- 14.1  the population and the collapse -------------------------------------------------------------

def test_in_estimate_is_in_model_AND_the_outcome_present_boolean_and_TOTAL():
    df, ps, est, _ = estimated()
    expected = ps.in_model & df[config.PRIMARY_OUTCOME].notna()
    assert est.in_estimate.dtype == bool
    assert not est.in_estimate.isna().any()
    assert est.in_estimate.index.equals(df.index)
    assert len(est.in_estimate) == len(df)
    assert list(est.in_estimate) == list(expected)


def test_on_ordinal_cohort_THE_COLLAPSE_FIRES_and_that_is_the_real_frame_witness():
    # The workbook occupies all seven levels and cannot provide one; this runs with no `data/`.
    _, _, est, _ = estimated()
    assert est.fit.categories == (0, 1, 3, 4, 6)
    assert len(est.fit.alpha) == 4
    assert tuple(v for v in config.MRS_LEVELS if v not in est.fit.categories) == (2, 5)


@pytest.mark.parametrize("constant", list(config.MRS_LEVELS))
def test_a_constant_response_raises_O5_naming_the_category_count(constant):
    df, audit = built(ordinal_cohort())
    ps = propensity.fit(df, audit)
    df = set_cell(df, "HAND-1", config.PRIMARY_OUTCOME, constant)
    for case_id in ("HAND-2", "HAND-5", "COHORT-1", "COHORT-3"):
        df = set_cell(df, case_id, config.PRIMARY_OUTCOME, constant)
    with pytest.raises(model.FitError) as e:
        outcome.primary(df, ps, audit)
    assert "O5" in message_of(e) and "1 response category" in message_of(e)


# --- 14.2  the [§11] denominator is per estimate, and no available frame varies it ----------------------
#
# The column §4.1 argues for, tested on a CONSTRUCTED frame because nothing available varies it: every
# workbook record and every `ordinal_cohort()` record has its outcome present. Without this section the
# mask is a copy of `in_model` on every frame in the suite, and an implementation that never built it
# would be green.


def _outcome_blanked():
    df, audit = built(ordinal_cohort())
    ps = propensity.fit(df, audit)
    return set_cell(df, "HAND-1", config.PRIMARY_OUTCOME, pd.NA), ps, audit


def test_blanking_ONE_outcome_moves_in_estimate_and_LEAVES_in_model_ALONE():
    df, ps, audit = _outcome_blanked()
    est = outcome.primary(df, ps, audit)
    assert int(ps.in_model.sum()) == 5              # HALF THE TEST: unchanged
    assert int(est.in_estimate.sum()) == 4
    entry = audit.entry("model", "outcome_completeness")
    assert entry.n == 1
    assert entry.case_ids == ("HAND-1",)            # and it NAMES the record


def test_WITHOUT_the_in_model_half_an_outcome_complete_caser_would_pass_by_removing_it_twice():
    """The assertion above is only meaningful beside this one.

    An implementation that complete-cased the propensity model on the OUTCOME too would satisfy
    "in_estimate is 4" by removing the record from both masks, and `in_model` would read 4 as well.
    `propensity.fit` is asserted not to know the outcome exists.
    """
    df, ps, _ = _outcome_blanked()
    refitted = propensity.fit(df, data.Audit(hand_source(Path("."), "refit")))
    assert int(refitted.in_model.sum()) == 5
    assert list(refitted.in_model) == list(ps.in_model)


def test_the_fit_is_computed_over_the_FOUR_with_their_own_weights():
    # The estimate shown to be over its OWN records and its OWN weights rather than over the table's.
    df, ps, audit = _outcome_blanked()
    est = outcome.primary(df, ps, audit)
    sub = df.loc[est.in_estimate]
    X, dropped = model.design(sub, (config.TREATMENT,))
    direct = model.polr(X, sub[config.PRIMARY_OUTCOME].to_numpy(dtype=float),
                        ps.w.loc[est.in_estimate].to_numpy(dtype=float))
    assert dropped == ()
    assert np.array_equal(est.fit.beta, direct.beta)
    assert np.array_equal(est.fit.alpha, direct.alpha)
    assert est.fit.categories == direct.categories


def test_THE_CONTRAST_blanking_a_PS_COVARIATE_moves_BOTH_masks_together():
    """What makes the three assertions above mean something rather than merely pass.

    And it is the second half of §4.1's rule: a record whose OUTCOME is present but whose COVARIATES
    are not is in NEITHER mask, so `in_estimate` is the intersection and not the outcome mask.
    """
    df, audit = built(ordinal_cohort())
    df = set_cell(df, "HAND-1", "age", pd.NA)
    ps = propensity.fit(df, audit)
    est = outcome.primary(df, ps, audit)
    assert int(ps.in_model.sum()) == 4                     # in_model moved
    assert int(est.in_estimate.sum()) == 4                 # and in_estimate moved with it
    assert df.loc[df["case_id"] == "HAND-1", config.PRIMARY_OUTCOME].notna().all()
    assert not bool(est.in_estimate[df["case_id"] == "HAND-1"].iloc[0])
    assert audit.entry("model", "outcome_completeness").n == 0   # the outcome mask removed nobody


# --- 14.5  the orientation, and the sign every reference implementation has [roadmap] -------------------

def test_THE_ORIENTATION_CRITERION_on_a_construction_whose_truth_is_KNOWN():
    """The construction's own truth is asserted FIRST, so the test states what "downward" means on
    its own data before it asks the fit about it."""
    X, y, w = orientation_frame()
    a = X[config.TREATMENT].to_numpy(dtype=float)
    treated_mean = float(np.mean(y[a == 1.0]))
    control_mean = float(np.mean(y[a == 0.0]))
    assert treated_mean < control_mean                     # measured 1.04 against 3.77
    fit = model.polr(X, y, w)
    assert float(np.exp(fit.beta[0])) > 1.0                # measured 58.35


def test_THE_COMPANION_the_other_convention_reports_the_OPPOSITE_direction():
    """What makes the criterion above a measurement rather than a restatement.

    A test that only checks `exp(beta) > 1` on data where the effect is large passes for an
    implementation that returns `exp(|beta|)`. The `alpha_k - beta*A` reading of the SAME fit is
    asserted strictly below 1 **and** equal to `1 / exp(beta)` — the reciprocal identity is what shows
    the companion is reading the same fit rather than a second one.
    """
    fit = model.polr(*orientation_frame())
    ours = float(np.exp(fit.beta[0]))
    theirs = float(np.exp(-fit.beta[0]))
    assert theirs < 1.0 < ours
    assert theirs == pytest.approx(1.0 / ours, rel=1e-12)


def test_THE_PILOTS_NEGATION_applied_to_OUR_beta_INVERTS_THE_PRIMARY_RESULT():
    """The single most consequential do-not-lift in this stage (§7.2, §18).

    `pilots/analysis.py:475` writes `float(np.exp(-res.params.iloc[0]))` with the comment "A positive
    coefficient means a shift towards higher (worse) mRS; invert so the common OR is oriented like
    every other outcome". That is CORRECT for `statsmodels`' `theta_k - x'beta` and WRONG for §5.1's
    `alpha_k + beta*A`, which is [§14a]'s own form. The negation is not a convention the pilots chose;
    it is the conversion between two parametrisations, and lifting the line without lifting the
    parametrisation inverts the primary result of the study.
    """
    fit = model.polr(*orientation_frame())
    lifted = float(np.exp(-fit.beta[0]))
    assert (lifted > 1.0) is not (float(np.exp(fit.beta[0])) > 1.0)


def test_the_unweighted_fit_against_STATSMODELS_with_the_ASYMMETRY_pinned_by_name():
    """§18b's second oracle, and BOTH halves are the point.

    **The coefficients come back NEGATED and the cutpoints do NOT.** A test asserting only that the
    two fits "agree" would pass on a comparison of the cutpoints alone, which is exactly how the
    mistake survives a review: a reader sees four cutpoints agreeing to five decimals and concludes
    the implementations agree.

    The third assertion is what makes the first two able to fail: `ours.beta - sm.beta` is asserted
    LARGER than 1e-2, so a frame on which both coefficients happened to be near zero cannot pass the
    sum test by having nothing to negate. Measured 1.425e-05, 8.947e-06 and 1.52.

    The tolerance is `statsmodels`' `bfgs` and not ours, and the test says so — §14.3's replication
    oracle is the one that holds at machine precision.
    """
    from statsmodels.miscmodels.ordinal_model import OrderedModel
    X, y = reference_frame()
    ours = model.polr(X, y)
    specification = OrderedModel(y, X.to_numpy(dtype=float), distr="logit")
    theirs = specification.fit(method="bfgs", disp=0)
    k = X.shape[1]
    their_beta = np.asarray(theirs.params[:k], dtype=float)
    their_thresholds = specification.transform_threshold_params(
        np.asarray(theirs.params[k:], dtype=float))[1:-1]
    assert float(np.max(np.abs(ours.beta + their_beta))) < 1e-4          # NEGATED
    assert float(np.max(np.abs(ours.alpha - their_thresholds))) < 1e-4   # and NOT negated
    assert float(np.max(np.abs(ours.beta - their_beta))) > 1e-2          # something to negate


def test_statsmodels_needs_no_importorskip_because_it_is_a_declared_dependency():
    # Checked rather than assumed (§2), as Stage 7 §12.12 checked the same thing. It is a project
    # dependency and not a `reference`-group one, so a plain checkout has it.
    import statsmodels                                     # noqa: F401  — the import IS the test
    assert statsmodels.__version__


@pytest.mark.parametrize("field,value", [("higher_is_better", True), ("kind", "binary")])
def test_orientation_RAISES_on_a_patched_registry(monkeypatch, field, value):
    """§7.1's assertion shown to FIRE rather than assumed to.

    If an ordinal outcome for which higher is better were ever registered as primary, a hardcoded
    sentence would print `> 1 favours bridging` beside a coefficient meaning the opposite, and nothing
    anywhere would disagree with it.
    """
    from dataclasses import replace
    registry = dict(config.OUTCOMES)
    registry[config.PRIMARY_OUTCOME] = replace(registry[config.PRIMARY_OUTCOME], **{field: value})
    monkeypatch.setattr(config, "OUTCOMES", registry)
    with pytest.raises(config.SchemaError) as e:
        outcome._orientation()
    assert config.PRIMARY_OUTCOME in message_of(e)
    assert "favours bridging" in message_of(e)


def test_the_orientation_sentence_is_COMPUTED_from_the_registry():
    sentence = outcome._orientation()
    assert config.TREATMENT_LABELS[1] in sentence
    assert config.OUTCOMES[config.PRIMARY_OUTCOME].label in sentence
    assert "alpha_k + beta*A" in sentence


# --- 14.6  the cumulative risk differences [roadmap, amended] -------------------------------------------
#
# `hand_ordinal()`'s weights are integers summing to 21 in each arm, so every weighted cumulative
# probability below is a sum of at most six of them over a denominator of 21 — which is what makes
# these the one set of VALUE assertions in this section, and it can be because that generator holds no
# seed.

HAND_RD = {0: +1 / 21, 1: -3 / 21, 2: -5 / 21, 3: -5 / 21, 4: -3 / 21, 5: +1 / 21}


@pytest.mark.parametrize("k", sorted(HAND_RD))
def test_RD_k_against_a_HAND_COMPUTED_weighted_proportion_difference(k):
    """Threshold by threshold, each against a literal computed by hand from that generator's twelve
    integer weights.

    Treated arm: y = 0,1,2,3,4,5 with w = 1,2,3,4,5,6 (Sigma 21). Control arm: y = 1,2,3,4,5,6 with
    w = 6,5,4,3,2,1 (Sigma 21). So P(Y <= 0) is 1/21 against 0/21, P(Y <= 1) is 3/21 against 6/21,
    and so on — every one a ratio of small integers.
    """
    y, a, w = hand_ordinal()
    rd, cumulative = outcome.cumulative_rd(y, a, w)
    assert rd[k] == pytest.approx(HAND_RD[k], abs=1e-12)
    treated = float(np.sum(w[(a == 1.0) & (y <= k)])) / 21.0
    control = float(np.sum(w[(a == 0.0) & (y <= k)])) / 21.0
    assert cumulative[k][1] == pytest.approx(treated, abs=1e-12)
    assert cumulative[k][0] == pytest.approx(control, abs=1e-12)


@pytest.mark.parametrize("name,triple", [
    ("hand_ordinal", hand_ordinal),
    ("deaths_both_arms", lambda: deaths_both_arms_frame()[1:]),
])
def test_the_cumulative_probabilities_are_NON_DECREASING_BY_CONSTRUCTION(name, triple):
    """**Satisfied by construction, and therefore unable to detect the forbidden route.**

    `P_w(Y <= k | arm)` is a cumulative sum of non-negative weights over a denominator that does not
    depend on `k`, so it is non-decreasing arithmetically. Measured: 2000 samples of a 40-record
    seven-category frame gave 0 non-monotone results by this route — and the threshold-model route
    gave a strict decrease in only 8 of 1990, by less than 5e-07. So this assertion holds of the code,
    of any correct alternative implementation, AND of a wrong one that reached the same numbers a
    different way, which is why §21's roadmap amendment 3 replaces the criterion rather than
    restating it. The companion below is the one that can fail.
    """
    y, a, w = triple()
    _, cumulative = outcome.cumulative_rd(y, a, w)
    for code in config.TREATMENT_LABELS:
        series = [cumulative[k][code] for k in config.MRS_THRESHOLDS]
        assert all(later >= earlier for earlier, later in zip(series, series[1:]))


def test_THE_COMPANION_the_six_threshold_model_route_gives_DIFFERENT_numbers():
    """The assertion that CAN fail: the route taken is asserted to be the empirical one.

    Different, not crossing — and that is the honest form. The crossing [§8] warns about is real and
    is rare: over 2000 samples of a 40-record seven-category frame the threshold route produced a
    strict decrease in 8 of the 1990 that fitted, 0.40%, and the largest decrease measured was below
    5e-07, which is float noise rather than a crossing anyone would see. So the argument for the
    empirical route rests on [§8]'s reasoning and on construction, not on a measured failure rate.
    """
    y, a, w = hand_ordinal()
    empirical, _ = outcome.cumulative_rd(y, a, w)
    X = pd.DataFrame({config.TREATMENT: a})
    through_models = {}
    for k in config.MRS_THRESHOLDS:
        indicator = (y <= float(k)).astype(float)
        if len(np.unique(indicator)) < 2:
            continue
        fit = model.polr(X, indicator, w)
        # P(indicator = 1 | arm) under the threshold model, which is P(Y <= k | arm)
        through_models[k] = tuple(
            float(1.0 / (1.0 + np.exp(-(fit.alpha[0] + fit.beta[0] * float(code)))))
            for code in config.TREATMENT_LABELS)
    assert through_models, "no threshold model fitted — the companion has no subject"
    differences = [abs(empirical[k] - (through_models[k][1] - through_models[k][0]))
                   for k in through_models]
    assert max(differences) > 1e-6


def test_RD_6_is_structurally_zero_and_is_NOT_REPORTED_while_RD_5_is_not():
    """The pair is the assertion: without the second, a reader concludes the last reported threshold
    is the trivial one. `RD_5` is the SURVIVAL contrast."""
    _, y, a, w = deaths_both_arms_frame()
    rd, cumulative = outcome.cumulative_rd(y, a, w)
    assert sorted(rd) == list(config.MRS_THRESHOLDS)
    assert len(rd) == 6 and 6 not in rd
    at_six = outcome.weighted_proportion((y <= 6.0).astype(float), a, w)
    assert at_six[1] == 1.0 and at_six[0] == 1.0
    assert at_six[1] - at_six[0] == 0.0                            # structurally zero
    assert rd[5] != 0.0                                            # and the survival contrast is not
    assert cumulative[5][1] < 1.0 and cumulative[5][0] < 1.0       # a death in each arm


def test_TIES_are_asserted_where_a_declared_level_is_unoccupied():
    # An implementation asserting strict monotonicity fails here and is wrong to: levels 2 and 5 are
    # unoccupied on `ordinal_cohort()`, so the cumulative probability cannot move across them.
    _, _, est, _ = estimated()
    assert 2 not in est.fit.categories and 5 not in est.fit.categories
    assert est.rd[1] == est.rd[2]
    assert est.rd[4] == est.rd[5]


def test_weighted_proportion_returns_NAN_and_not_ZERO_for_a_zero_weight_arm():
    x = np.array([1.0, 0.0, 1.0, 0.0])
    a = np.array([1.0, 1.0, 0.0, 0.0])
    share = outcome.weighted_proportion(x, a, np.array([0.5, 0.5, 0.0, 0.0]))
    assert np.isnan(share[0])
    assert share[1] == pytest.approx(0.5, abs=1e-12)
    assert set(share) == set(config.TREATMENT_LABELS)               # keyed by arm, BOTH arms


def test_WITHOUT_the_guard_np_average_RAISES_rather_than_returning_nan():
    # Stage 7 §3.1's fact, one stage on, on the same expression — which is why the check comes BEFORE
    # the mean rather than after it.
    with pytest.raises(ZeroDivisionError):
        np.average(np.array([1.0, 0.0]), weights=np.array([0.0, 0.0]))


def test_cumulative_rd_takes_EXACTLY_y_a_w_and_no_threshold_list():
    # MRS_THRESHOLDS is prespecified at Stage 1 for this stage; a defaulted keyword would make a
    # prespecified choice look like an option, and would fail this test rather than pass silently.
    signature = inspect.signature(outcome.cumulative_rd)
    assert list(signature.parameters) == ["y", "a", "w"]
    assert all(p.default is inspect.Parameter.empty for p in signature.parameters.values())


def test_the_arm_codes_are_read_from_the_REGISTRY_and_not_written_as_1_and_0():
    """§8.1. The sign of `RD_k` is the one thing in this module a reader cannot check by reading, so
    a third declared arm code must fail here rather than silently make `RD_k` a contrast of the
    extremes."""
    assert config.TREATMENT_LABELS[outcome._TREATED] == config.TREATMENT_LABELS[1]
    assert config.TREATMENT_LABELS[outcome._COMPARATOR] == config.TREATMENT_LABELS[0]
    assert outcome._TREATED != outcome._COMPARATOR
    assert len(config.TREATMENT_LABELS) == 2


# --- 14.7  the separation guard, and G7 is what turns a converged fit into a failure --------------------
#
# `polr`'s half of this section is `test_model.py`'s — that a perfectly separated fit CONVERGES, and
# the three detectors that do not work. This is the guard's half.

def test_primary_on_a_SEPARATED_cohort_raises_G7_naming_the_odds_ratio_it_rejected():
    """So a reader of a Stage 10 failure log can tell a separated replicate from a non-converged one.

    The `Propensity` is the REAL one, fitted before the outcome was overwritten, so the weights are a
    genuine [§7] score and only the response is constructed — which is what makes this raise for the
    reason §6 gives rather than because the frame is unreachable.
    """
    df, ps, audit = separated_cohort()
    with pytest.raises(model.FitError) as e:
        outcome.primary(df, ps, audit)
    text = message_of(e)
    assert "G7" in text
    assert "separation" in text and "does NOT present as non-convergence" in text
    assert "the common odds ratio is" in text                 # the rejected value, in the message
    assert f"bound of {config.POLR_MAX_ABS_BETA:g}" in text


def test_G7_raises_FitError_and_NOT_SchemaError_because_that_is_what_stage_10_catches():
    # [§10] catches `FitError` to drop and count a replicate and may catch nothing else: a
    # SchemaError here is a bug in the resampler, not a sparse replicate (propensity.py:468-469).
    df, ps, audit = separated_cohort()
    with pytest.raises(model.FitError):
        outcome.primary(df, ps, audit)
    assert not issubclass(model.FitError, config.SchemaError)


def test_WITHOUT_G7_the_separated_fit_RETURNS_an_odds_ratio_no_manuscript_could_print(monkeypatch):
    """DoD-5, and it is the only finding in this stage that was found by running the code.

    Nothing looks wrong: the fit converges on the score criterion with every safeguard counter at
    zero and every fitted probability finite, and `exp(beta)` comes back as a finite float that every
    downstream table will accept.
    """
    monkeypatch.setattr(outcome, "_assert_reportable", lambda fit: None)
    df, ps, audit = separated_cohort()
    est = outcome.primary(df, ps, audit)
    assert est.fit.converged_on == "score"
    assert (est.fit.rescales, est.fit.halvings) == (0, 0)
    assert np.isfinite(est.odds_ratio) and est.odds_ratio > 1e12
    assert abs(est.beta) > config.POLR_MAX_ABS_BETA


def test_the_bound_is_read_from_CONFIG_and_never_from_a_literal(monkeypatch):
    # Asserted by patching the constant low and watching a HEALTHY fit start raising. It lives in
    # `config.py` and not in a default argument for the reason FIRTH_TOL does: [§10] refits this
    # model in every one of N_BOOT replicates, so a bound is a property of the sampling distribution
    # and not a runtime knob, and changing it changes which replicates are dropped.
    frame, audit = built(ordinal_cohort())
    ps = propensity.fit(frame, audit)
    outcome.primary(frame, ps, audit)                       # healthy, does not raise
    monkeypatch.setattr(config, "POLR_MAX_ABS_BETA", 1.0)
    with pytest.raises(model.FitError) as e:
        outcome.primary(frame, ps, data.Audit(audit.source))
    assert "G7" in message_of(e) and "bound of 1" in message_of(e)


def test_the_guard_bounds_BETA_and_NOT_ALPHA(monkeypatch):
    """Asserted by PATCHING THE BOUND rather than by finding a frame whose `alpha` exceeds 14.0 —
    because no such frame is constructible.

    `alpha` for a rare category grows like `log n`: measured, one rare record in 61 gives
    `max|alpha|` 4.16, in 601 gives 6.40, in 6001 gives 8.70 and in 60001 gives 11.00, so reaching
    14.0 needs order 1e7 records. An earlier draft asserted that the test "constructs one", and it
    cannot.

    The constructible form: on `rare_category_frame()` — `max|alpha|` 4.1599, `|beta|` 0.126998 — the
    guard PASSES with the bound patched to 3.0, where a bound applied to `alpha` would raise, and
    raises only when the bound goes below `|beta|`. That asserts the field the guard reads.
    """
    fit = model.polr(*rare_category_frame())
    assert float(np.max(np.abs(fit.alpha))) > 3.0           # an alpha-bound at 3.0 would reject it
    assert float(np.max(np.abs(fit.beta))) < 3.0

    monkeypatch.setattr(config, "POLR_MAX_ABS_BETA", 3.0)
    outcome._assert_reportable(fit)                          # and the real guard passes

    monkeypatch.setattr(config, "POLR_MAX_ABS_BETA", 0.1)
    with pytest.raises(model.FitError) as e:
        outcome._assert_reportable(fit)
    assert "G7" in message_of(e)


# --- 14.9  the three audit entries, and the log --------------------------------------------------------

STEPS = ("outcome_completeness", "primary_fit", "cumulative_rd")


def test_the_three_entries_appear_in_the_declared_order_at_CAPTURED_positions():
    # Never a tail slice, which is the repair Stage 5 §10 landed and Stages 6 and 7 carried.
    df, audit = built(ordinal_cohort())
    ps = propensity.fit(df, audit)
    before = len(audit.entries)
    outcome.primary(df, ps, audit)
    assert len(audit.entries) - before == 3
    assert [e.step for e in audit.entries[before:before + 3]] == list(STEPS)
    assert [e.kind for e in audit.entries[before:before + 3]] == ["model"] * 3


def test_KINDS_stays_the_declared_NINE_and_no_heading_is_added():
    """§9.1's claim in test form, and the assertion that fails if someone adds an `estimate` kind.

    `data.py:151-152` declares that Stages 7-13 all render under `model` or under an existing kind.
    Stage 7 tested that claim once; this stage tests it twice, and `data.py` is untouched — which is
    why none of the literal pins in the six existing test modules moves.
    """
    assert len(data.KINDS) == 9
    assert data.KINDS == ("provenance", "contract", "correction", "observation", "derivation",
                          "cohort", "model", "structural", "missingness")
    assert set(data._HEADINGS) == set(data.KINDS)
    assert data._HEADINGS["model"] == "Fitted models"


def test_all_three_render_under_Fitted_models_AFTER_overlap_weights():
    _, _, _, audit = estimated()
    rendered = audit.to_markdown()
    section = rendered.split("## Fitted models")[1].split("\n## ")[0]
    positions = [section.index(f"- **{step}**") for step in ("overlap_weights", *STEPS)]
    assert positions == sorted(positions)


def test_each_entrys_n_is_a_DIFFERENT_QUANTITY():
    """Entry 1's is records removed by the outcome mask, entry 2's is records estimated, entry 3's is
    `len(MRS_THRESHOLDS)`. Asserted on a frame where the three differ, because that is the property a
    reader needs."""
    df, ps, audit = _outcome_blanked()
    outcome.primary(df, ps, audit)
    ns = [audit.entry("model", step).n for step in STEPS]
    assert ns == [1, 4, 6]
    assert len(set(ns)) == 3


def test_entry_1_is_recorded_even_when_its_n_is_ZERO():
    """The v7 case. An entry that appears only when it has something to say is an entry whose absence
    a reader has to interpret, and [§11] asks for the denominator of EVERY estimate."""
    _, _, _, audit = estimated()
    entry = audit.entry("model", "outcome_completeness")
    assert entry is not None
    assert entry.n == 0 and entry.case_ids == ()
    assert entry.table is not None and len(entry.table) == 5     # header + four population rows


def test_entry_2_has_one_row_per_DECLARED_mRS_level_including_the_absent_ones():
    _, _, est, audit = estimated()
    rows = rows_of(audit.entry("model", "primary_fit"))
    for level in config.MRS_LEVELS:
        assert f"alpha at mRS <= {level}" in rows
    for absent in (2, 5):
        level, estimate, status = rows[f"alpha at mRS <= {absent}"]
        assert estimate == "missing"
        assert status.startswith("absent from this population")


def test_the_TOP_FITTED_LEVEL_status_is_asserted_BY_NAME_on_truncated_cohort():
    """§14.0.5's frame, and the branch §9.2 corrects.

    **The previous test passes on the broken table**: the row that rendered `missing` under status
    `fitted` is not an absent row, so an assertion about absent rows never looks at it. `categories`
    is asserted FIRST, because without it a change to which records survive would silently turn this
    fixture back into one whose top fitted level is the top declared one and the branch would go back
    to being untested with every test still green.
    """
    _, _, est, audit = estimated(truncated_cohort())
    assert est.fit.categories == (0, 1, 2, 3, 4)               # NOT MRS_LEVELS[-1] at the top
    rows = rows_of(audit.entry("model", "primary_fit"))

    _, estimate, status = rows["alpha at mRS <= 4"]
    assert estimate == "missing"
    assert status == "highest FITTED level — no cutpoint, P(Y <= it) = 1"
    assert status != "fitted"

    _, estimate, status = rows["alpha at mRS <= 6"]
    assert estimate == "missing"
    assert status.startswith("absent from this population")
    assert "highest declared level" in status

    _, estimate, status = rows["alpha at mRS <= 5"]
    assert estimate == "missing" and status.startswith("absent from this population")
    assert "highest declared level" not in status


NO_CUTPOINT_STATUSES = (
    "absent from this population — no cutpoint [§5.3]",
    "absent from this population — no cutpoint [§5.3], and it is the highest declared level",
    "highest FITTED level — no cutpoint, P(Y <= it) = 1",
)
FIT_TABLE_STATUSES = NO_CUTPOINT_STATUSES + (
    "fitted", "treatment coefficient [§8]", "common odds ratio — see the orientation above")


@pytest.mark.parametrize("name,frame", list(FIXTURE_COHORTS.items()))
def test_the_estimate_column_reads_MISSING_IF_AND_ONLY_IF_there_is_no_cutpoint(name, frame):
    """The invariant that generalises §9.2's correction, and catches every miswiring of that branch
    including ones nobody anticipated. Verified on both fixture frames, and VIOLATED on
    `truncated_cohort()` under the draft's `level == MRS_LEVELS[-1]` branch."""
    _, _, _, audit = estimated(frame())
    for row in audit.entry("model", "primary_fit").table[1:]:
        _, _, estimate, status = row
        if not row[0].startswith("alpha at mRS <= "):
            continue
        assert (estimate == "missing") is (status in NO_CUTPOINT_STATUSES), row


@pytest.mark.parametrize("name,frame", list(FIXTURE_COHORTS.items()))
def test_the_status_column_holds_only_the_SIX_declared_strings(name, frame):
    # SIX, not four: an earlier draft of this bullet said four while `_fit_table` emitted five, and
    # §9.2's correction makes it six. Asserted against a literal tuple, so a seventh status added
    # without a test fails here — and never `True`, `1` or `nan`.
    _, _, _, audit = estimated(frame())
    statuses = {row[3] for row in audit.entry("model", "primary_fit").table[1:]}
    assert statuses <= set(FIT_TABLE_STATUSES)
    assert all(isinstance(s, str) and s for s in statuses)
    assert len(FIT_TABLE_STATUSES) == 6


@pytest.mark.parametrize("step", STEPS)
def test_no_cell_holds_a_pipe_and_the_rendered_grid_is_RECTANGULAR(step):
    """`data._md_table` does no escaping and sizes its separator from `len(rows[0])`, so one pipe
    inside a header cell gives a header row with more markdown cells than its own separator and body
    — and byte-identity across hash seeds passes the whole time, because a table identically broken
    under both seeds is still identical (Stage 7 §21.1). Look at the log, not only at its hash.

    The column that would have carried one is entry 3's cumulative probability, whose natural header
    is `P(Y <= k | bridging)`; it is written `P(Y <= k), bridging` and the mathematical bar is
    deliberately not used.
    """
    _, _, _, audit = estimated()
    table = audit.entry("model", step).table
    assert not [cell for row in table for cell in row if "|" in cell]
    lines = data._md_table(table).splitlines()
    assert len({len(line.split("|")) for line in lines}) == 1


def test_that_grid_assertion_FIRES_on_a_pasted_pipe():
    broken = (("threshold", "P(Y <= k | bridging)"), ("mRS <= 0", "0.1"))
    assert len({len(line.split("|")) for line in data._md_table(broken).splitlines()}) == 2


def test_entry_2s_detail_reads_as_PROSE_under_its_ABSENT_LEVELS_branch():
    """Stage 7 §21.5a found an unterminated clause in exactly this position that three fence
    assemblies and two review rounds walked past, because the branch never rendered. This stage
    asserts the branch's prose from the start, on the frame that reaches it."""
    _, _, est, audit = estimated()
    detail = audit.entry("model", "primary_fit").detail
    assert ("so 4 cutpoint(s) were estimated; 2 declared level(s) are absent from this population "
            "and are marked as such below: 2, 5. ") in detail
    assert detail.endswith(".")
    assert "  " not in detail                                   # no doubled space at a join
    assert "below: 2, 5. exp(beta) is the common odds ratio" in detail


def test_entry_2s_detail_reads_as_PROSE_under_its_NO_ABSENT_LEVELS_branch():
    # The other branch, on a frame that occupies every declared level, so the clause is absent and
    # the sentence still terminates. `deaths_both_arms_frame()` has all seven in both arms.
    _, y, a, w = deaths_both_arms_frame()
    fit = model.polr(pd.DataFrame({config.TREATMENT: a}), y, w)
    assert fit.categories == tuple(config.MRS_LEVELS)
    detail = outcome._fit_detail(fit, 14)
    assert "so 6 cutpoint(s) were estimated. exp(beta) is the common odds ratio" in detail
    assert "declared level(s) are absent" not in detail
    assert detail.endswith(".")


def test_entry_3s_detail_reads_as_PROSE_under_its_PER_ARM_TOTALS_clause():
    _, _, _, audit = estimated()
    detail = audit.entry("model", "cumulative_rd").detail
    assert "denominators [§11]: EVT alone " in detail
    assert ", bridging " in detail
    assert "mRS <= 6 is not reported" in detail and "mRS <= 5 is the survival contrast" in detail
    assert detail.endswith(".")


@pytest.mark.parametrize("step", STEPS)
def test_each_details_interpolated_counts_MATCH_the_table_beside_it(step):
    _, _, est, audit = estimated()
    entry = audit.entry("model", step)
    if step == "primary_fit":
        fitted_rows = [r for r in entry.table[1:] if r[3] == "fitted"]
        assert f"{len(fitted_rows)} cutpoint(s) were estimated" in entry.detail
        assert len(fitted_rows) == len(est.alpha)
    if step == "cumulative_rd":
        assert f"{len(entry.table) - 1} declared threshold(s)" in entry.detail
        assert len(entry.table) - 1 == len(config.MRS_THRESHOLDS)
    if step == "outcome_completeness":
        assert f"{int(est.in_estimate.sum())} of " in entry.detail


def test_two_runs_render_identical_markdown():
    assert estimated()[3].to_markdown() == estimated()[3].to_markdown()


def test_the_log_is_identical_across_interpreters_with_different_hash_seeds(tmp_path):
    """The two-seed driver of Stage 5 §12.9, extended one stage.

    Written out rather than sketched, for Stage 3 §12.12's reason: an implementer choosing it freely
    can choose one that fits nothing at all and still see two identical outputs.
    """
    driver = textwrap.dedent("""
        import sys
        sys.path.insert(0, %r)
        sys.path.insert(0, %r)
        from test_outcome import estimated
        _, _, _, audit = estimated()
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
    for step in STEPS:
        assert f"- **{step}**" in outputs[0]


# --- 14.10  the estimate and the risk differences are ONE contrast --------------------------------------

def test_the_fit_and_the_RD_are_computed_over_THE_SAME_MASK_bound_once():
    """§4.2's "bound once" as an assertion rather than a code comment.

    Two presentations of one contrast computed over two populations are two contrasts, and both
    numbers would be finite and each correct for its own population.
    """
    df, ps, est, _ = estimated()
    sub = df.loc[est.in_estimate]
    X, _ = model.design(sub, (config.TREATMENT,))
    y = sub[config.PRIMARY_OUTCOME].to_numpy(dtype=float)
    a = sub[config.TREATMENT].to_numpy(dtype=float)
    w = ps.w.loc[est.in_estimate].to_numpy(dtype=float)

    direct_fit = model.polr(X, y, w)
    assert np.array_equal(est.fit.beta, direct_fit.beta)
    assert np.array_equal(est.fit.alpha, direct_fit.alpha)

    direct_rd, direct_cumulative = outcome.cumulative_rd(y, a, w)
    assert est.rd == direct_rd
    assert est.cumulative == direct_cumulative


@pytest.mark.parametrize("name,frame", list(FIXTURE_COHORTS.items()))
def test_the_orientation_agrees_between_the_two_scales_THROUGH_THE_WEIGHTED_MEAN(name, frame):
    """Under §5.1's parametrisation a positive `beta` moves mass to LOW mRS in the treated arm, so the
    treated arm's weighted mean mRS must be the lower one. A sign error in either scale alone breaks
    this, and neither alone would be caught by §14.5 or §14.6.

    **THE PER-THRESHOLD FORM OF THIS ASSERTION IS FALSE**, and an earlier draft asserted it: "on a
    frame where `beta > 0`, the `RD_k` are non-negative at every `k` where the two arms' cumulative
    probabilities differ" is not a property of a proportional-odds fit. `beta` is a single summary of
    a shift the data need not exhibit at every cut, which is exactly the assumption §15 files as
    untested — so where proportional odds fails, `beta` and an individual `RD_k` legitimately
    disagree in sign. **The suite's own fixtures are counter-examples**: `truncated_cohort()` has
    `beta = +0.3448` while `RD_3 = -0.4531`, because both control records sit at mRS <= 3 so
    `P(Y <= 3 | comparator)` is exactly 1.0. The weighted mean is a monotone functional of the same
    distributions and is the form that holds.
    """
    df, ps, est, _ = estimated(frame())
    sub = df.loc[est.in_estimate]
    y = sub[config.PRIMARY_OUTCOME].to_numpy(dtype=float)
    a = sub[config.TREATMENT].to_numpy(dtype=float)
    w = ps.w.loc[est.in_estimate].to_numpy(dtype=float)
    means = outcome.weighted_proportion(y, a, w)
    assert np.sign(est.beta) == -np.sign(means[outcome._TREATED] - means[outcome._COMPARATOR])


def test_the_PER_THRESHOLD_form_is_FALSE_on_this_suites_own_fixtures():
    # Asserted so the weaker-looking form above is not mistaken for timidity: this is the measurement
    # that made the strong form unwritable, and it would have failed at T6 on correct code.
    _, _, est, _ = estimated(truncated_cohort())
    assert est.beta > 0.0
    assert est.rd[3] < 0.0
    assert est.cumulative[3][outcome._COMPARATOR] == 1.0

    y, a, w = hand_ordinal()
    rd, _ = outcome.cumulative_rd(y, a, w)
    fit = model.polr(pd.DataFrame({config.TREATMENT: a}), y, w)
    assert float(fit.beta[0]) < 0.0
    assert rd[0] > 0.0 and rd[5] > 0.0


def test_the_STRICT_per_k_form_is_KEPT_where_proportional_odds_holds_BY_DESIGN():
    # `orientation_frame()` is the frame on which the strong form means something: the treated arm's
    # mass is at mRS 0-3 and the control's at 2-6 by construction, so every RD_k carries beta's sign.
    X, y, w = orientation_frame()
    a = X[config.TREATMENT].to_numpy(dtype=float)
    fit = model.polr(X, y, w)
    rd, _ = outcome.cumulative_rd(y, a, w)
    assert float(fit.beta[0]) > 0.0
    assert all(rd[k] > 0.0 for k in config.MRS_THRESHOLDS)


def test_entry_3s_per_arm_totals_reconcile_with_stage_6s_overlap_weights_THROUGH_FMT():
    """Stage 7 §12.8's rule: an `AuditEntry.table` holds strings and a full-precision float will not
    equal one, so the comparison goes through `_fmt` on both sides. Read out of the SAME `Audit`."""
    _, _, _, audit = estimated()
    overlap = rows_of(audit.entry("model", "overlap_weights"))
    detail = audit.entry("model", "cumulative_rd").detail
    for code, label in config.TREATMENT_LABELS.items():
        sum_w = overlap[label][1]                       # ("n", "sum w", "ESS", …)
        assert f"{label} {sum_w}" in detail


# --- 14.11  Stage 8 adds nothing, refits nothing, and reads no diagnostic --------------------------------

def test_the_frame_comes_back_unchanged_cell_for_cell():
    df, audit = built(ordinal_cohort())
    ps = propensity.fit(df, audit)
    before = df.copy(deep=True)
    outcome.primary(df, ps, audit)
    pd.testing.assert_frame_equal(df, before)
    assert list(df.columns) == list(before.columns)
    assert df.index.equals(before.index)
    assert df.dtypes.equals(before.dtypes)


def test_the_propensity_comes_back_unchanged_and_no_attribute_is_added_to_the_frame():
    df, audit = built(ordinal_cohort())
    ps = propensity.fit(df, audit)
    e, w, in_model = ps.e.copy(), ps.w.copy(), ps.in_model.copy()
    columns = list(df.columns)
    outcome.primary(df, ps, audit)
    pd.testing.assert_series_equal(ps.e, e)
    pd.testing.assert_series_equal(ps.w, w)
    pd.testing.assert_series_equal(ps.in_model, in_model)
    assert list(df.columns) == columns


def _imported_roots(source: str) -> set[str]:
    roots: set[str] = set()
    for node in ast.walk(ast.parse(source)):
        if isinstance(node, ast.Import):
            roots.update(alias.name.split(".")[0] for alias in node.names)
        if isinstance(node, ast.ImportFrom) and node.module:
            roots.add(node.module.split(".")[0])
    return roots


def test_outcome_imports_NEITHER_BALANCE_nor_a_second_estimator():
    """**The `balance` half is the one that matters**: importing it would not fail, and Stage 7 §11
    states that Stage 8 does not read the `Balance`. No estimate in this SAP is conditioned on a
    balance diagnostic, and a pipeline in which a threshold silently selects an estimator is the
    mixture [§7] forbids one level up."""
    roots = _imported_roots(SOURCE)
    assert "balance" not in roots
    assert roots & {"statsmodels", "sklearn", "scipy"} == set()
    assert roots >= {"config", "model", "propensity", "data", "numpy", "pandas"}


def test_that_import_scan_FIRES():
    pasted = "import balance\nimport statsmodels.api as sm\nfrom scipy import optimize\n"
    assert _imported_roots(pasted) >= {"balance", "statsmodels", "scipy"}


def test_outcome_never_names_propensity_fit():
    # It READS a Propensity; it does not fit one. A refit here would be a second score for the same
    # population, and [§10] refits the whole pipeline per replicate.
    calls = {f"{node.func.value.id}.{node.func.attr}"
             for node in ast.walk(ast.parse(SOURCE))
             if isinstance(node, ast.Call) and isinstance(node.func, ast.Attribute)
             and isinstance(node.func.value, ast.Name)}
    assert "propensity.fit" not in calls
    assert "model.polr" in calls and "model.design" in calls


def test_primary_takes_EXACTLY_df_ps_audit():
    # No outcome name, no covariate list, no threshold list, no bound — [§8] fixes all four, and a
    # defaulted keyword fails this test rather than passing silently.
    signature = inspect.signature(outcome.primary)
    assert list(signature.parameters) == ["df", "ps", "audit"]
    assert all(p.default is inspect.Parameter.empty for p in signature.parameters.values())


PRIMARY_FIELDS = ("beta", "odds_ratio", "alpha", "cut_levels", "rd", "cumulative",
                  "in_estimate", "fit")
POLR_FIT_FIELDS = ("beta", "alpha", "categories", "columns", "iterations", "converged_on",
                   "first_step_norm", "rescales", "halvings")


@pytest.mark.parametrize("cls,fields", [(outcome.Primary, PRIMARY_FIELDS),
                                        (model.PolrFit, POLR_FIT_FIELDS)])
def test_the_result_types_carry_EXACTLY_their_declared_fields(cls, fields):
    import dataclasses
    assert tuple(f.name for f in dataclasses.fields(cls)) == fields


@pytest.mark.parametrize("cls", [outcome.Primary, model.PolrFit])
def test_NO_FIELD_ANYWHERE_carries_a_standard_error_or_an_interval(cls):
    """§5.6 is the argument and this is the assertion that makes it survive a later edit.

    The naive weighted-likelihood covariance is WRONG here, not merely out of scope: `-H^-1` at the
    optimum is the inverse observed information of a likelihood in which `w_i` counts observations,
    and the overlap weights are a tilting function of an ESTIMATED propensity score, so the sampling
    variability of `e` enters the estimate and `-H^-1` omits it entirely. A `se` attribute on the
    object holding `beta` would be in a table within two stages.
    """
    import dataclasses
    forbidden = ("se", "cov", "ci", "pval", "interval", "stderr")
    offenders = [f.name for f in dataclasses.fields(cls)
                 if any(word in f.name.lower() for word in forbidden)]
    assert offenders == []


def test_Primary_has_NO_METHOD_RETURNING_A_VERDICT():
    """[§8] specifies no threshold on `beta`: the null is tested at Stage 10 from the bootstrap
    distribution, and a method here returning `odds_ratio > 1.0` would be a significance test with no
    interval behind it, one attribute access away from a reporting layer. Stage 7's `Balance` carries
    `worst()` and `unbalanced()` because [§9] specifies a threshold; this one carries nothing."""
    methods = [name for name, value in vars(outcome.Primary).items()
               if callable(value) and not name.startswith("__")]
    assert methods == []
    for forbidden in ("significant", "favours_bridging", "unbalanced", "worst"):
        assert not hasattr(outcome.Primary, forbidden)


def test_WITHOUT_G6_a_one_armed_frame_gives_EITHER_an_IndexError_OR_an_odds_ratio_of_ONE():
    """DoD-9, and the point is that NEITHER is a failure the caller can see.

    `design`'s constant-column rule is right for a covariate (Stage 6 §4.3) and catastrophic for the
    SOLE predictor: the design comes back with no treatment column, `polr` fits an intercept-only
    ordinal model over the cutpoints alone, CONVERGES, and returns a `beta` of length zero. Neither
    function raises. Whether the caller then gets an `IndexError` on `beta[0]` or a silent `1.0` from
    `exp(beta.sum())` depends on how the reporting layer happens to be written, and both are worse
    than a raise — which is why `primary` reads `fit.beta[0]` and G6 is what makes the index safe.
    """
    df, audit = built(ordinal_cohort())
    ps = propensity.fit(df, audit)
    one_armed = df[df[config.TREATMENT] == 1]

    X, dropped = model.design(one_armed, (config.TREATMENT,))
    assert X.shape[1] == 0 and dropped == (config.TREATMENT,)        # `design` does not raise

    fit = model.polr(X, one_armed[config.PRIMARY_OUTCOME].to_numpy(dtype=float),
                     ps.w.loc[one_armed.index].to_numpy(dtype=float))
    assert len(fit.beta) == 0 and fit.iterations >= 1                # and `polr` CONVERGES

    with pytest.raises(IndexError):
        float(fit.beta[0])                                           # reading 1
    assert float(np.exp(fit.beta.sum())) == 1.0                      # reading 2 — "no effect"

    # and G6 is what stops either from happening
    with pytest.raises(model.FitError) as e:
        outcome._assert_exposure_survived(X, dropped)
    assert "G6" in message_of(e)


@pytest.mark.parametrize("name", ["outcome.py", "model.py"])
def test_neither_module_contains_a_bare_assert(name):
    # Python strips `assert` under -O, so under a flag nobody remembers setting, a module whose whole
    # purpose is to fail loudly would succeed silently.
    source = (MODULE_DIR / name).read_text(encoding="utf-8")
    assert [n.lineno for n in ast.walk(ast.parse(source)) if isinstance(n, ast.Assert)] == []


@pytest.mark.parametrize("name", ["outcome.py", "model.py", "test_outcome.py"])
def test_no_module_of_this_stage_is_exempt_from_the_raw_name_scan(name):
    assert name not in config.EXEMPT_FROM_RAW_NAME_SCAN


@pytest.mark.parametrize("name", ["outcome.py", "model.py"])
def test_neither_module_names_a_RAW_HEADER(name):
    source = (MODULE_DIR / name).read_text(encoding="utf-8")
    contract = set(config.COLUMN_CONTRACT)
    offenders = [node.value for node in ast.walk(ast.parse(source))
                 if isinstance(node, ast.Constant) and isinstance(node.value, str)
                 and node.value in contract]
    assert offenders == []


def test_that_raw_name_scan_FIRES():
    """The header is BUILT FROM THE CONTRACT and never written as a literal, and that is not
    fastidiousness: `test_config.py`'s repository-wide raw-name scan walks every `.py` under the tree
    that is not in `EXEMPT_FROM_RAW_NAME_SCAN`, and DoD-3 requires `test_outcome.py` never to be
    exempt. A pasted header here would make this file its own only offender, and the repair would be
    to gut the scan (Stage 6 §12.7's rule).
    """
    contract = set(config.COLUMN_CONTRACT)
    header = next(raw for raw, column in config.COLUMN_CONTRACT.items()
                  if column.name == config.PRIMARY_OUTCOME)
    pasted = f'x = df[{header!r}]\n'
    assert [n.value for n in ast.walk(ast.parse(pasted))
            if isinstance(n, ast.Constant) and n.value in contract] == [header]


# --- 14.12  properties from the workbook [data-gated] ---------------------------------------------------
#
# THIS SECTION ASSERTS PROPERTIES AND NOT VALUES, AND §4.3 IS WHY. Stage 7 §12.12 pinned 1.343 and
# 0.380; the analogous numbers here are the study's result. Every assertion below can fail, and none
# of them quotes the answer — no `beta`, no `exp(beta)`, no `RD_k`, and no mRS distribution split by
# arm. The counts are pinned deliberately: 93 / 92 / 92, seven levels and six cutpoints are counts and
# not estimates. §14.12a's meta-assertion is what keeps that true as this file is edited.

@DATA_GATED
def test_workbook_the_three_populations_and_an_entry_that_removes_nobody(workbook, workbook_primary):
    df, ps, audit = workbook
    assert len(df) == 93
    assert int(ps.in_model.sum()) == 92
    assert int(workbook_primary.in_estimate.sum()) == 92
    assert audit.entry("model", "outcome_completeness").n == 0
    # On v7 the third population EQUALS the second, which is exactly the condition under which an
    # implementation that never built the mask would be green — §14.2 is the constructed test.
    assert list(workbook_primary.in_estimate) == list(ps.in_model)


@DATA_GATED
def test_workbook_all_seven_declared_levels_are_occupied_so_the_collapse_is_INERT(workbook_primary):
    # Recorded as inert on this workbook rather than assumed to be. `ordinal_cohort()` is the frame
    # that exercises the branch, and it runs with no `data/`.
    assert workbook_primary.fit.categories == tuple(config.MRS_LEVELS)
    assert len(workbook_primary.fit.alpha) == 6
    assert len(workbook_primary.alpha) == 6
    assert workbook_primary.cut_levels == config.MRS_THRESHOLDS


@DATA_GATED
def test_workbook_the_fit_CONVERGED_on_a_declared_route(workbook_primary):
    fit = workbook_primary.fit
    assert fit.converged_on in ("likelihood", "score")
    assert 1 <= fit.iterations < config.POLR_MAX_ITER


@DATA_GATED
def test_workbook_alpha_is_STRICTLY_ASCENDING(workbook_primary):
    alpha = np.asarray(workbook_primary.alpha, dtype=float)
    assert len(alpha) == 6
    assert bool(np.all(np.diff(alpha) > 0.0))


@DATA_GATED
def test_workbook_G7_DOES_NOT_FIRE_on_v7(workbook_primary):
    """Asserted, because a guard that fired on the real data would mean the primary analysis has no
    estimate — and that is a finding for the PI, not a passing test.

    The separation guard is measured to fire on constructions and has never fired on data; this
    assertion is the first thing that will tell anyone, which is the right place for it and is not the
    same as knowing.
    """
    assert abs(workbook_primary.beta) < config.POLR_MAX_ABS_BETA
    assert float(np.max(np.abs(workbook_primary.fit.beta))) < config.POLR_MAX_ABS_BETA


@DATA_GATED
def test_workbook_the_odds_ratio_is_FINITE_POSITIVE_and_is_exp_of_beta(workbook_primary):
    assert np.isfinite(workbook_primary.odds_ratio)
    assert workbook_primary.odds_ratio > 0.0
    assert workbook_primary.odds_ratio == pytest.approx(
        float(np.exp(workbook_primary.beta)), rel=1e-12)


@DATA_GATED
def test_workbook_the_cumulative_probabilities_are_NON_DECREASING_in_k_within_each_arm(
        workbook_primary):
    for code in config.TREATMENT_LABELS:
        series = [workbook_primary.cumulative[k][code] for k in config.MRS_THRESHOLDS]
        assert all(later >= earlier for earlier, later in zip(series, series[1:]))
        assert all(0.0 <= value <= 1.0 for value in series)
    assert len(workbook_primary.rd) == 6


@DATA_GATED
def test_workbook_each_arms_P_of_mRS_at_most_5_matches_whether_that_arm_HAS_A_DEATH(
        workbook, workbook_primary):
    """§8.3's statement tied to the workbook without quoting a difference: `P(Y <= 5)` is strictly
    below 1 exactly where the arm contains a death, and exactly 1 where it does not."""
    df, _, _ = workbook
    estimated_rows = df.loc[workbook_primary.in_estimate]
    for code in config.TREATMENT_LABELS:
        arm = estimated_rows[config.TREATMENT] == code
        has_death = bool((arm & (estimated_rows[config.PRIMARY_OUTCOME]
                                 == config.MRS_LEVELS[-1])).any())
        share = workbook_primary.cumulative[config.MRS_THRESHOLDS[-1]][code]
        assert (share < 1.0) is has_death
        assert (share == 1.0) is not has_death


@DATA_GATED
def test_workbook_the_orientation_agrees_between_the_scales_THROUGH_THE_WEIGHTED_MEAN(
        workbook, workbook_primary):
    """The strongest property available without quoting a magnitude, and the one that catches a sign
    error introduced downstream of §14.5's synthetic construction.

    **It reports no magnitude and no distribution split by arm**: a sign is not an estimate.

    **NOT the per-threshold form**, and that matters more here than in §14.10 because this section
    runs on the workbook. A draft asserted that every `RD_k` carries `beta`'s sign wherever the arms
    differ; that is false wherever proportional odds fails, which §15 already declares untested — so
    asserted on v7 it could go RED FOR A LEGITIMATE REASON, with a message saying "sign error", on the
    one stage where §4.3 makes that impossible to distinguish in advance. A data-gated test that can
    fail for a reason the analysis already declares untested converts a known limitation into an alarm
    at the exact moment the PI is reading the estimate for the first time.
    """
    df, ps, _ = workbook
    sub = df.loc[workbook_primary.in_estimate]
    means = outcome.weighted_proportion(
        sub[config.PRIMARY_OUTCOME].to_numpy(dtype=float),
        sub[config.TREATMENT].to_numpy(dtype=float),
        ps.w.loc[workbook_primary.in_estimate].to_numpy(dtype=float))
    assert np.sign(workbook_primary.beta) == -np.sign(
        means[outcome._TREATED] - means[outcome._COMPARATOR])


@DATA_GATED
def test_workbook_entry_3s_per_arm_totals_reconcile_with_overlap_weights(workbook, workbook_primary):
    _, _, audit = workbook
    overlap = rows_of(audit.entry("model", "overlap_weights"))
    detail = audit.entry("model", "cumulative_rd").detail
    for code, label in config.TREATMENT_LABELS.items():
        assert f"{label} {overlap[label][1]}" in detail


@DATA_GATED
def test_workbook_the_three_entries_render_and_their_counts_match_their_tables(
        workbook, workbook_primary):
    _, _, audit = workbook
    assert len(audit.entries) == 25            # 22 after fit, and this fixture does not call assess
    fit_entry = audit.entry("model", "primary_fit")
    assert fit_entry.n == 92
    assert len([r for r in fit_entry.table[1:] if r[3] == "fitted"]) == 6
    assert "6 cutpoint(s) were estimated" in fit_entry.detail
    assert "declared level(s) are absent" not in fit_entry.detail
    rd_entry = audit.entry("model", "cumulative_rd")
    assert rd_entry.n == 6
    assert len(rd_entry.table) - 1 == 6
    assert "6 declared threshold(s)" in rd_entry.detail
    for step in STEPS:
        assert f"- **{step}**" in audit.to_markdown()


# --- 14.12a  the float-literal meta-assertion that keeps §4.3 true -------------------------------------
#
# Its boundary and its allowlist are specified rather than left to the implementer, because a scan
# whose scope is "this section" is either unwritable or vacuous.

SECTION_BANNER = "# --- 14.12  "
# The NEXT `# ---` banner at column 0, in FILE order and not in section number: `test_balance.py`'s
# banners are not in numeric order, so a scan keyed on `14.13` would silently run to the end of the
# file. The convention already exists in test_cohort.py, test_model.py, test_propensity.py and
# test_balance.py, so the delimiter is a fact about the repository rather than something this invents.
NEXT_BANNER = "\n# --- "


def _section_source(source: str, banner: str) -> str:
    start = source.index(banner)
    end = source.index(NEXT_BANNER, start + len(banner))
    return source[start:end]


_SECTION_REFERENCE = re.compile(r"§\s*\d+(?:\.\d+)*[a-z]?")
_MAGNITUDE = re.compile(r"\d+\.\d+(?:e[-+]?\d+)?|\d+e[-+]?\d+", re.IGNORECASE)


def _unallowlisted_float_literals(section: str) -> list[object]:
    """Every magnitude in `section` outside the declared allowlist — as a float literal OR in a
    docstring.

    ALLOWED, and it is short by design: (a) the value of a `rel=` or `abs=` keyword, or the second
    positional argument of `pytest.approx`; and (b) `0.0` or `1.0`, which are the two boundary values
    §14.12 asserts `odds_ratio` and the cumulative probabilities against and which carry no magnitude.
    Everything else fails. INTEGERS ARE UNRESTRICTED: 93 / 92 / 92, the six cutpoints and the seven
    levels are counts, and §14.12 pins them deliberately.

    **THE DOCSTRING HALF IS NOT OPTIONAL AND IT IS WHY THIS IS NOT A ONE-LINE AST WALK.** The natural
    place to paste `exp(beta)` while reading the log is the docstring of the assertion it justifies,
    where it is an `ast.Constant` of type `str` and no float-literal scan can see it. So string
    constants are searched too, with `§`-prefixed numbers stripped first — every numeric-looking
    substring §14.12 legitimately contains is a section reference (`§4.3`, `§8.3`, `§14.5`, `§14.10`),
    measured, and nothing else.

    **A float in a `#` comment IS invisible to this**, because comments are not in the AST. That is a
    real hole and a small one, stated rather than papered over: §15's first bullet and DoD-17's
    read-by-eye stand behind it, and a comment carrying the estimate does not make an assertion depend
    on it.
    """
    tree = ast.parse(textwrap.dedent(section))
    allowed: set[int] = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.keyword) and node.arg in ("rel", "abs"):
            allowed.add(id(node.value))
        if (isinstance(node, ast.Call) and isinstance(node.func, ast.Attribute)
                and node.func.attr == "approx" and len(node.args) > 1):
            allowed.add(id(node.args[1]))

    offenders: list[object] = []
    for node in ast.walk(tree):
        if not isinstance(node, ast.Constant):
            continue
        if isinstance(node.value, float):
            if node.value not in (0.0, 1.0) and id(node) not in allowed:
                offenders.append(node.value)
        elif isinstance(node.value, str):
            offenders += _MAGNITUDE.findall(_SECTION_REFERENCE.sub("", node.value))
    return offenders


def test_the_section_boundary_is_FOUND_AT_BOTH_ENDS_and_the_range_is_NON_EMPTY():
    """Asserted BEFORE anything is asserted about what is in it.

    A scan that silently matched nothing is the green-test failure `test_the_raw_name_scan_actually_
    fires` (Stage 1 §7) exists to prevent, three modules earlier.
    """
    assert SECTION_BANNER in THIS_FILE
    section = _section_source(THIS_FILE, SECTION_BANNER)
    assert section.startswith(SECTION_BANNER)
    assert len(section.splitlines()) > 20
    assert "def test_workbook_" in section
    assert "# --- 14.12a" not in section                 # the terminator is not inside the range


def test_NO_FLOAT_LITERAL_outside_the_allowlist_appears_in_section_14_12():
    # This is what keeps §4.3 true as the file is edited: after this stage there is no pinned
    # regression number on the primary effect anywhere in git, and a magnitude pasted here while
    # reading the log would put one there quietly.
    assert _unallowlisted_float_literals(_section_source(THIS_FILE, SECTION_BANNER)) == []


def test_that_float_scan_FIRES_on_a_pasted_magnitude_AND_on_one_in_a_DOCSTRING():
    pasted = 'def test_x():\n    assert est.beta == pytest.approx(0.7412, abs=1e-6)\n'
    assert _unallowlisted_float_literals(pasted) == [0.7412]

    # DoD-17's second half, and the one the scan exists for: the natural place to paste exp(beta)
    # while reading the log is the docstring of the assertion it justifies.
    in_docstring = 'def test_x():\n    """The workbook gives exp(beta) = 2.0984."""\n    assert 1\n'
    assert _unallowlisted_float_literals(in_docstring) == ["2.0984"]

    both = ('def test_x():\n    """exp(beta) is 2.0984."""\n'
            '    assert est.beta == pytest.approx(0.7412, abs=1e-6)\n')
    assert sorted(str(v) for v in _unallowlisted_float_literals(both)) == ["0.7412", "2.0984"]


def test_that_float_scan_does_NOT_fire_on_the_forms_the_allowlist_declares():
    allowed = 'def test_x():\n    assert a == pytest.approx(b, rel=1e-12)\n    assert c > 0.0\n'
    assert _unallowlisted_float_literals(allowed) == []

    # section references are stripped before the string search, or every docstring in this file
    # would be an offender — measured: §4.3, §8.3, §14.5 and §14.10 are the only numeric-looking
    # substrings §14.12 legitimately contains
    citations = 'def test_x():\n    """See §4.3, §8.3, §14.5 and §14.10."""\n    assert 1\n'
    assert _unallowlisted_float_literals(citations) == []

    # and integers are unrestricted: 93 / 92 / 92, six cutpoints and seven levels are COUNTS
    counts = 'def test_x():\n    """93 records, 92 in_model, 6 cutpoints."""\n    assert n == 93\n'
    assert _unallowlisted_float_literals(counts) == []


# --- 14.13  the module boundary -------------------------------------------------------------------------

def _config_attributes(source: str) -> set[str]:
    return {node.attr for node in ast.walk(ast.parse(source))
            if isinstance(node, ast.Attribute) and isinstance(node.value, ast.Name)
            and node.value.id in {"C", "config"}}


def test_model_names_neither_the_exposure_nor_any_outcome_after_this_stages_amendment():
    """Stage 6 §12.12's scan, extended to the constants this stage adds — `POLR_*` passes it and an
    `OUTCOME_MAX_ABS_BETA` would not, which is §12's reason for the prefix."""
    attributes = _config_attributes((MODULE_DIR / "model.py").read_text(encoding="utf-8"))
    assert [name for name in attributes
            if name.startswith(("PS_", "OUTCOME_", "PRIMARY_"))] == []
    assert {"TREATMENT", "PRIMARY_OUTCOME", "MRS_LEVELS", "MRS_THRESHOLDS"} & attributes == set()
    assert {"POLR_MAX_ITER", "POLR_TOL", "POLR_SCORE_TOL", "POLR_MAX_HALVINGS",
            "POLR_MAX_STEP", "POLR_ETA_CLIP"} <= attributes
    # and POLR_MAX_ABS_BETA is NOT among them: G7 is the caller's, for Stage 6 F5's reason — "a
    # boundary value is a degenerate fit and not a confident one" is a statement about the
    # specification and not about the arithmetic (§6.3).
    assert "POLR_MAX_ABS_BETA" not in attributes
    assert "POLR_MAX_ABS_BETA" in _config_attributes(SOURCE)


def test_model_names_no_ivt_or_mrs_STRING_LITERAL():
    source = (MODULE_DIR / "model.py").read_text(encoding="utf-8")
    literals = [node.value for node in ast.walk(ast.parse(source))
                if isinstance(node, ast.Constant) and isinstance(node.value, str)]
    # scoped past the module docstring, which legitimately names `propensity.py` and the SAP
    body = literals[1:]
    assert [s for s in body if "ivt" in s or "mrs_" in s] == []


def test_outcome_is_the_ONLY_SHIPPED_MODULE_naming_PRIMARY_OUTCOME():
    shipped = sorted(p for p in MODULE_DIR.glob("*.py"))
    naming = [p.name for p in shipped
              if "PRIMARY_OUTCOME" in p.read_text(encoding="utf-8") and p.name != "config.py"]
    assert naming == ["outcome.py"]
    assert "PRIMARY_OUTCOME" not in (MODULE_DIR / "model.py").read_text(encoding="utf-8")


# --- 14.14  the preconditions ---------------------------------------------------------------------------

def _reindexed(ps: propensity.Propensity, which: str, index) -> propensity.Propensity:
    """`ps` with ONE of its three Series carried on a different index. §14.14."""
    fields = {"e": ps.e, "w": ps.w, "in_model": ps.in_model}
    fields[which] = fields[which].reindex(index)
    return propensity.Propensity(**fields, ess=ps.ess, fit=ps.fit, dropped=ps.dropped)


MISALIGNMENTS = {
    "superset": lambda idx: pd.Index(list(idx) + [10 ** 6]),
    "subset": lambda idx: idx[:-1],
    "permutation": lambda idx: idx[::-1],
    "dtype": lambda idx: idx.astype(str),
}


@pytest.mark.parametrize("series", ["e", "w", "in_model"])
@pytest.mark.parametrize("kind", sorted(MISALIGNMENTS))
def test_G1_fires_on_every_series_and_every_kind_of_misalignment(series, kind):
    """Parametrised over the three Series AND the four kinds, so G1 is COMPLETE for §10 rather than
    assumed to be. Stage 7 §21.4 found that an `e`-only misalignment reached `_centre` and raised
    pandas' `IndexingError` after an entry had been recorded; the lesson is inherited rather than
    re-learned.

    `e` alone is the one case that CANNOT produce a wrong estimate, because `primary` never reads `e`
    — G1 is stricter than §10 strictly needs, which is the right side to err on.
    """
    df, audit = built(ordinal_cohort())
    ps = propensity.fit(df, audit)
    bad = _reindexed(ps, series, MISALIGNMENTS[kind](df.index))
    assert not getattr(bad, series).index.equals(df.index)      # all four DO fail index.equals
    with pytest.raises(config.SchemaError) as e:
        outcome.primary(df, bad, audit)
    text = message_of(e)
    assert "G1" in text and series in text
    assert "G4 and G5 were not run" in text


def test_WITHOUT_G1_a_PERMUTED_propensity_gives_a_plausible_estimate_for_MIS_WEIGHTED_records(
        monkeypatch):
    """What G1 PREVENTS, which no earlier draft measured.

    Every other guard in this stage has a with-it-removed witness producing the wrong number, and G1
    had none while its own message claims it is the difference between an estimate and a finite common
    odds ratio for a population that does not exist. **The witness is a permutation, not `e`**:
    `primary` never reads `e`, so an `e`-only misalignment cannot produce a wrong estimate at all.

    Measured: with G1 removed, a `Propensity` carrying the frame's labels in reversed order gives an
    odds ratio of 58.01 against the correct 5.76, with nothing raising, because `w` arrives ordered by
    the `Propensity`'s index while the response and the arm vector come from the frame's.
    """
    df, audit = built(ordinal_cohort())
    ps = propensity.fit(df, audit)
    correct = outcome.primary(df, ps, audit)

    monkeypatch.setattr(outcome, "_assert_primary_inputs", lambda frame, score: None)
    permuted = _reindexed(_reindexed(_reindexed(ps, "e", df.index[::-1]),
                                     "w", df.index[::-1]), "in_model", df.index[::-1])
    wrong = outcome.primary(df, permuted, data.Audit(audit.source))
    assert wrong.odds_ratio == pytest.approx(58.0097, abs=1e-3)
    assert correct.odds_ratio == pytest.approx(5.7638, abs=1e-3)
    assert abs(wrong.odds_ratio - correct.odds_ratio) > 50.0


def test_WITHOUT_G1_a_SUBSET_index_leaves_the_denominator_in_SILENCE():
    """The quieter second witness: the missing labels resolve to False inside the `&`, so the record
    leaves the [§11] denominator with nothing raising and nothing to notice."""
    df, audit = built(ordinal_cohort())
    ps = propensity.fit(df, audit)
    subset = ps.in_model.reindex(df.index[:-1])
    silent = subset & df[config.PRIMARY_OUTCOME].notna()
    assert silent.dtype == bool
    assert int(silent.sum()) == int(ps.in_model.sum()) - 1
    assert not silent.isna().any()


def test_WITHOUT_G1_a_SUPERSET_index_is_BENIGN_which_is_why_G1_is_stricter_than_needed(monkeypatch):
    df, audit = built(ordinal_cohort())
    ps = propensity.fit(df, audit)
    correct = outcome.primary(df, ps, audit)
    monkeypatch.setattr(outcome, "_assert_primary_inputs", lambda frame, score: None)
    extra = pd.Index(list(df.index) + [10 ** 6])
    benign = _reindexed(_reindexed(_reindexed(ps, "e", extra), "w", extra), "in_model", extra)
    same = outcome.primary(df, benign, data.Audit(audit.source))
    assert same.odds_ratio == pytest.approx(correct.odds_ratio, rel=1e-12)


@pytest.mark.parametrize("blank", [None, pd.NA, np.nan], ids=["None", "pd.NA", "np.nan"])
def test_G2_fires_on_a_three_valued_mask(blank):
    df, audit = built(ordinal_cohort())
    ps = propensity.fit(df, audit)
    mask = ps.in_model.astype(object)
    mask.iloc[0] = blank
    bad = propensity.Propensity(e=ps.e, w=ps.w, in_model=mask, ess=ps.ess, fit=ps.fit,
                                dropped=ps.dropped)
    with pytest.raises(config.SchemaError) as e:
        outcome.primary(df, bad, audit)
    assert "G2" in message_of(e) and "boolean and TOTAL" in message_of(e)


@pytest.mark.parametrize("blank", [None, pd.NA, np.nan], ids=["None", "pd.NA", "np.nan"])
def test_WITHOUT_G2_the_record_leaves_the_denominator_and_NOTHING_RAISES(blank):
    """The companion asserts the SILENT exclusion, and an earlier draft said the input "dies with
    pandas' ValueError". It does not: `object & bool` COERCES the missing entry to False inside the
    `&`, so the mask comes back a clean boolean Series of the frame's length with 4 records estimated
    against `in_model`'s 5 and no exception anywhere.

    **That makes G2 stronger, not weaker** — a guard whose absence is silent is a guard worth having,
    and one whose absence raises anyway is nearly redundant.
    """
    df, audit = built(ordinal_cohort())
    ps = propensity.fit(df, audit)
    mask = ps.in_model.astype(object)
    mask.iloc[0] = blank
    silent = mask & df[config.PRIMARY_OUTCOME].notna()
    assert silent.dtype == bool
    assert not silent.isna().any()
    assert len(silent) == len(df)
    assert int(silent.sum()) == int(ps.in_model.sum()) - 1


@pytest.mark.parametrize("column", ["PRIMARY_OUTCOME", "TREATMENT", "case_id"])
def test_G3_fires_on_a_frame_missing_a_column_it_reads(column):
    name = getattr(config, column, column)
    df, audit = built(ordinal_cohort())
    ps = propensity.fit(df, audit)
    with pytest.raises(config.SchemaError) as e:
        outcome.primary(df.drop(columns=[name]), ps, audit)
    assert "G3" in message_of(e) and name in message_of(e)


def test_WITH_G3_IN_PHASE_2_the_outcome_case_dies_with_a_bare_KeyError_instead():
    """§4.4a's fourth failure, and the companion that would have caught Stage 7's incomplete first
    repair — applied here before the fact.

    Phase 2's first line builds the mask from `df[C.PRIMARY_OUTCOME]`, so naming a column in an
    absent-column check is necessary and is NOT sufficient unless that check raises before the column
    is read.
    """
    df, audit = built(ordinal_cohort())
    ps = propensity.fit(df, audit)
    dropped = df.drop(columns=[config.PRIMARY_OUTCOME])
    with pytest.raises(KeyError) as e:
        _ = ps.in_model & dropped[config.PRIMARY_OUTCOME].notna()
    assert config.PRIMARY_OUTCOME in str(e.value)
    # and with G3 where it belongs, the same frame gets a message instead
    with pytest.raises(config.SchemaError) as schema:
        outcome.primary(dropped, ps, audit)
    assert "G3" in message_of(schema)


def test_G4_fires_on_a_one_armed_estimation_population_naming_BOTH_arms_and_their_counts():
    df, audit = built(ordinal_cohort())
    ps = propensity.fit(df, audit)
    one_armed = df[df[config.TREATMENT] == 1]
    ps_one = propensity.Propensity(
        e=ps.e.loc[one_armed.index], w=ps.w.loc[one_armed.index],
        in_model=ps.in_model.loc[one_armed.index], ess=ps.ess, fit=ps.fit, dropped=ps.dropped)
    with pytest.raises(config.SchemaError) as e:
        outcome.primary(one_armed, ps_one, audit)
    text = message_of(e)
    assert "G4" in text
    for label in config.TREATMENT_LABELS.values():
        assert label in text
    assert "EVT alone: 0" in text and "bridging: 3" in text


def test_G5_fires_on_a_HAND_BUILT_propensity_with_a_nan_weight_inside_in_estimate():
    """Unreachable from `propensity.fit` — `w` is nan off `in_model` by construction — so it is
    constructed, in the posture Stage 7 §12.7 declares.

    It is not defensive. It is the only check that names the DELETION mechanism, and the mechanism is
    the one thing about this stage a reader will not guess: a nan weight here does not produce a nan
    estimate, it produces a different, finite, plausible estimate over a coarser outcome scale.
    """
    df, audit = built(ordinal_cohort())
    ps = propensity.fit(df, audit)
    w = ps.w.copy()
    w.iloc[0] = np.nan
    bad = propensity.Propensity(e=ps.e, w=w, in_model=ps.in_model, ess=ps.ess, fit=ps.fit,
                                dropped=ps.dropped)
    with pytest.raises(config.SchemaError) as e:
        outcome.primary(df, bad, audit)
    text = message_of(e)
    assert "G5" in text and "1 record(s)" in text
    assert "DELETES" in text and "Never fill it" in text


@pytest.mark.parametrize("dropped", [(config.TREATMENT,), ()])
def test_G6_fires_on_a_width_zero_design_and_raises_FitError_NOT_SchemaError(dropped):
    with pytest.raises(model.FitError) as e:
        outcome._assert_exposure_survived(pd.DataFrame(index=range(3)), dropped)
    assert "G6" in message_of(e)
    assert not isinstance(e.value, config.SchemaError)


def test_G6_fires_on_a_design_that_grew_an_extra_column():
    # `dropped` is computed from the MATRIX and not from the arm counts, so any future change to
    # `design`'s rule reaches here first.
    X = pd.DataFrame({config.TREATMENT: [1.0, 0.0], "age": [70.0, 71.0]})
    with pytest.raises(model.FitError) as e:
        outcome._assert_exposure_survived(X, ())
    assert "G6" in message_of(e) and "SOLE predictor" in message_of(e)


def test_G6_PASSES_on_the_design_primary_actually_builds():
    df, audit = built(ordinal_cohort())
    X, dropped = model.design(df, (config.TREATMENT,))
    assert tuple(X.columns) == (config.TREATMENT,) and dropped == ()
    outcome._assert_exposure_survived(X, dropped)                 # does not raise


G1_TO_G5 = ("G1", "G2", "G3", "G4", "G5")


def _broken(kind: str):
    """(df, ps) failing exactly one of G1-G5, for the records-nothing assertion."""
    df, audit = built(ordinal_cohort())
    ps = propensity.fit(df, audit)
    if kind == "G1":
        return df, _reindexed(ps, "w", df.index[::-1]), audit
    if kind == "G2":
        mask = ps.in_model.astype(object)
        mask.iloc[0] = pd.NA
        return df, propensity.Propensity(e=ps.e, w=ps.w, in_model=mask, ess=ps.ess, fit=ps.fit,
                                         dropped=ps.dropped), audit
    if kind == "G3":
        return df.drop(columns=[config.TREATMENT]), ps, audit
    if kind == "G4":
        one = df[df[config.TREATMENT] == 1]
        return one, propensity.Propensity(
            e=ps.e.loc[one.index], w=ps.w.loc[one.index], in_model=ps.in_model.loc[one.index],
            ess=ps.ess, fit=ps.fit, dropped=ps.dropped), audit
    w = ps.w.copy()
    w.iloc[0] = np.nan
    return df, propensity.Propensity(e=ps.e, w=w, in_model=ps.in_model, ess=ps.ess, fit=ps.fit,
                                     dropped=ps.dropped), audit


@pytest.mark.parametrize("kind", G1_TO_G5)
def test_G1_to_G5_append_NOTHING_AT_ALL_to_the_audit(kind):
    # They are preconditions and they run before any `audit.record`. Delta is exactly 0.
    df, ps, audit = _broken(kind)
    before = len(audit.entries)
    with pytest.raises(config.SchemaError):
        outcome.primary(df, ps, audit)
    assert len(audit.entries) - before == 0


@pytest.mark.parametrize("guard,make", [
    ("G7", separated_cohort),
])
def test_G6_and_G7_append_EXACTLY_ENTRY_ONE_and_never_entries_two_or_three(guard, make):
    """The one place this stage departs from Stage 7's "records nothing", and the departure IS the
    point rather than a leak.

    Both fire after entry 1 would naturally be written, and `primary` orders
    `_record_outcome_completeness` before the design deliberately: the [§11] denominator is a fact
    about the frame established before any modelling choice, and a log that names it is useful
    precisely when what follows fails (propensity.py:440-443).

    And a log carrying the coefficient of a fit that was then REJECTED would render under
    `## Fitted models` beside a `missing` estimate and read as a partial run — which is why
    `_assert_reportable` sits between the fit and `audit.record`.
    """
    df, ps, audit = make()
    before = len(audit.entries)
    with pytest.raises(model.FitError) as e:
        outcome.primary(df, ps, audit)
    assert guard in message_of(e)
    assert len(audit.entries) - before == 1
    assert audit.entries[before].step == "outcome_completeness"
    assert audit.entry("model", "primary_fit") is None
    assert audit.entry("model", "cumulative_rd") is None


def test_a_frame_failing_G1_AND_G4_reports_ONLY_G1_and_says_the_others_were_not_run():
    df, audit = built(ordinal_cohort())
    ps = propensity.fit(df, audit)
    one = df[df[config.TREATMENT] == 1]
    bad = propensity.Propensity(
        e=ps.e, w=ps.w, in_model=ps.in_model, ess=ps.ess, fit=ps.fit, dropped=ps.dropped)
    with pytest.raises(config.SchemaError) as e:
        outcome.primary(one, bad, audit)
    text = message_of(e)
    assert "G1" in text
    # `G4  ` with its own two-space prefix is the CHECK; the bare "G4" below is the sentence saying
    # it was never run, which is the half that has to be there.
    assert "G4  " not in text
    assert "G4 and G5 were not run" in text


def test_G1_and_G2_together_give_ONE_error_naming_TWO():
    df, audit = built(ordinal_cohort())
    ps = propensity.fit(df, audit)
    mask = ps.in_model.astype(object)
    mask.iloc[0] = pd.NA
    bad = propensity.Propensity(e=ps.e, w=ps.w.reindex(df.index[::-1]), in_model=mask,
                                ess=ps.ess, fit=ps.fit, dropped=ps.dropped)
    with pytest.raises(config.SchemaError) as e:
        outcome.primary(df, bad, audit)
    text = message_of(e)
    assert "G1" in text and "G2" in text
    assert "2 primary-estimate assertion(s) failed" in text


def test_G4_and_G5_together_give_ONE_error_naming_TWO():
    df, audit = built(ordinal_cohort())
    ps = propensity.fit(df, audit)
    one = df[df[config.TREATMENT] == 1]
    w = ps.w.loc[one.index].copy()
    w.iloc[0] = np.nan
    bad = propensity.Propensity(
        e=ps.e.loc[one.index], w=w, in_model=ps.in_model.loc[one.index],
        ess=ps.ess, fit=ps.fit, dropped=ps.dropped)
    with pytest.raises(config.SchemaError) as e:
        outcome.primary(one, bad, audit)
    text = message_of(e)
    assert "G4" in text and "G5" in text
    assert "2 primary-estimate assertion(s) failed" in text
