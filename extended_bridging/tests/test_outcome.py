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
import dataclasses
import hashlib
import inspect
import os
import re
import subprocess
import sys
import textwrap
from pathlib import Path
from typing import Final

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


# --- Stage 11 §15.7  the [§11] mask is ONE definition, and it is public --------------------------------

def _masked():
    """(df, ps, audit) over the ordinal fixture cohort — the frame this stage can actually fit."""
    df, audit = built(ordinal_cohort())
    return df, propensity.fit(df, audit), audit


def test_estimation_population_IS_the_mask_primary_binds():
    """The extraction is proved EQUIVALENT and not merely intended to be (Stage 11 §8.5, §15.12).

    Stage 11's subgroup replicate body needs this mask on a drawn frame and cannot have the `Primary`
    that carries it without paying for a primary fit it does not read. So the mask became a named
    definition with two callers — and the one thing that could go wrong is the extraction not being
    the same mask, which would move every downstream number at once.
    """
    df, ps, audit = _masked()
    mask = outcome.estimation_population(df, ps)
    est = outcome.primary(df, ps, audit)
    assert mask.equals(est.in_estimate)
    assert mask.dtype == bool and not mask.isna().any()          # boolean and TOTAL — §4.1
    assert mask.index.equals(df.index)


def test_estimation_population_is_in_model_AND_the_outcome_being_present():
    """A THIRD population and not a restatement of `in_model`: a record with complete covariates and
    a missing outcome keeps its weight and loses its estimate (§4.1)."""
    df, ps, _ = _masked()
    blanked = df.copy()
    blanked[config.PRIMARY_OUTCOME] = blanked[config.PRIMARY_OUTCOME].astype("Int64")
    blanked.loc[blanked.index[0], config.PRIMARY_OUTCOME] = pd.NA
    mask = outcome.estimation_population(blanked, ps)
    assert int(mask.sum()) == int(ps.in_model.sum()) - int(bool(ps.in_model.iloc[0]))
    assert (mask <= ps.in_model).all()


def test_the_mask_is_SPELLED_ONCE_in_the_module():
    """`primary` and `_assert_primary_inputs` both need it, and Stage 11 §6.1 rejects exactly this
    shape for the family partition: two computations of one prespecified thing with nothing
    asserting they agree."""
    assert SOURCE.count("ps.in_model & df[C.PRIMARY_OUTCOME].notna()") == 1
    definition = SOURCE[SOURCE.index("def estimation_population"):]
    definition = definition[:definition.index("\n\n\n")]
    assert "ps.in_model & df[C.PRIMARY_OUTCOME].notna()" in definition


def test_the_three_primary_step_names_go_through_the_specification():
    """`primary` is called TWICE over one cohort from Stage 11 onward and `Audit.entry` is
    first-match, so two unsuffixed `primary_fit` entries would make every programmatic read return
    the [§7] estimate's while the rendered log looked complete (Stage 11 §4.4, §10)."""
    df, ps, audit = _masked()
    before = len(audit.entries)
    outcome.primary(df, ps, audit)
    primary_steps = [e.step for e in audit.entries[before:]]
    assert primary_steps == ["outcome_completeness", "primary_fit", "cumulative_rd"]

    arm = propensity.fit_full(df, audit)
    between = len(audit.entries)
    outcome.primary(df, arm, audit)
    arm_steps = [e.step for e in audit.entries[between:]]
    assert arm_steps == [f"{base}_full_covariate" for base in primary_steps]
    assert not set(primary_steps) & set(arm_steps)
    assert audit.entry("model", "primary_fit").n == int(ps.in_model.sum())


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


def test_only_outcome_and_sensitivity_NAME_PRIMARY_OUTCOME_and_the_second_is_SPENT_DELIBERATELY():
    """Stage 8 §12's rule, and Stage 11 spends it once — from the side that section did not consider.

    §12 made `outcome.py` the only shipped module that may name the [§5] primary outcome, so that no
    second estimator of it could appear anywhere else. Stage 11 §9 quotes the rule from the other
    direction, declining to put the subgroup `FitError` raises in `outcome.py` because *"a subgroup
    estimator in it would be that rule spent for a scan's convenience"* — and then implements the
    subgroup estimator in `sensitivity.py`, which must therefore name the response it fits.

    **It is unavoidable and it is bounded.** [§13]'s subgroup clause is explicitly *on the primary
    outcome*, so the module implementing it names the primary outcome or takes the column from a
    parameter no prespecified signature has. What the rule EXISTS to prevent is preserved and is
    asserted below: `sensitivity.py` names no OTHER [§5] outcome key, and the only ordinal fit it
    performs is the [§13] interaction design — the primary estimate itself is `outcome.primary`'s and
    is called, never reimplemented.
    """
    shipped = sorted(p for p in MODULE_DIR.glob("*.py"))
    naming = [p.name for p in shipped
              if "PRIMARY_OUTCOME" in p.read_text(encoding="utf-8") and p.name != "config.py"]
    assert naming == ["outcome.py", "sensitivity.py"]
    assert "PRIMARY_OUTCOME" not in (MODULE_DIR / "model.py").read_text(encoding="utf-8")

    sensitivity_source = (MODULE_DIR / "sensitivity.py").read_text(encoding="utf-8")
    for key in config.BINARY_OUTCOMES:
        assert key not in sensitivity_source
    assert "outcome.primary(" in sensitivity_source              # called, never reimplemented
    assert sensitivity_source.count("model.polr(") == 1          # the ONE fit Stage 11 adds


# --- 14.14  the preconditions ---------------------------------------------------------------------------

def _reindexed(ps: propensity.Propensity, which: str, index) -> propensity.Propensity:
    """`ps` with ONE of its three Series carried on a different index. §14.14."""
    fields = {"e": ps.e, "w": ps.w, "in_model": ps.in_model}
    fields[which] = fields[which].reindex(index)
    return propensity.Propensity(**fields, ess=ps.ess, fit=ps.fit, dropped=ps.dropped,
                                 spec=config.PROPENSITY_PRIMARY)


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
                                dropped=ps.dropped,
                                spec=config.PROPENSITY_PRIMARY)
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
        in_model=ps.in_model.loc[one_armed.index], ess=ps.ess, fit=ps.fit, dropped=ps.dropped,
        spec=config.PROPENSITY_PRIMARY)
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
                                dropped=ps.dropped,
                                spec=config.PROPENSITY_PRIMARY)
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
                                         dropped=ps.dropped,
                                         spec=config.PROPENSITY_PRIMARY), audit
    if kind == "G3":
        return df.drop(columns=[config.TREATMENT]), ps, audit
    if kind == "G4":
        one = df[df[config.TREATMENT] == 1]
        return one, propensity.Propensity(
            e=ps.e.loc[one.index], w=ps.w.loc[one.index], in_model=ps.in_model.loc[one.index],
            ess=ps.ess, fit=ps.fit, dropped=ps.dropped,
            spec=config.PROPENSITY_PRIMARY), audit
    w = ps.w.copy()
    w.iloc[0] = np.nan
    return df, propensity.Propensity(e=ps.e, w=w, in_model=ps.in_model, ess=ps.ess, fit=ps.fit,
                                     dropped=ps.dropped,
                                     spec=config.PROPENSITY_PRIMARY), audit


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
        e=ps.e, w=ps.w, in_model=ps.in_model, ess=ps.ess, fit=ps.fit, dropped=ps.dropped,
        spec=config.PROPENSITY_PRIMARY)
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
                                ess=ps.ess, fit=ps.fit, dropped=ps.dropped,
                                spec=config.PROPENSITY_PRIMARY)
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
        ess=ps.ess, fit=ps.fit, dropped=ps.dropped,
        spec=config.PROPENSITY_PRIMARY)
    with pytest.raises(config.SchemaError) as e:
        outcome.primary(one, bad, audit)
    text = message_of(e)
    assert "G4" in text and "G5" in text
    assert "2 primary-estimate assertion(s) failed" in text


# ======================================================================================================
# Stage 9 — §15 of `specs/stage9_secondary_binary_estimators.md`
# ======================================================================================================
#
# `predict`'s sections — §15.6, and §15.7's separation witness — live in `test_model.py`, because that
# is where the function is. Everything here reads a `Secondary`, an estimator or an `Audit`.
#
# **§15's data-gated tests assert PROPERTIES AND COUNTS, never a weighted quantity**, and §4.3 is why,
# narrowed from Stage 8's position rather than inherited unexamined. Two of these seven outcomes already
# have a CRUDE arm contrast on the public record (`../out/stage0_data_inventory.md`), so pretending the
# direction is unknown would be a fiction. What is nonetheless absent from git is every WEIGHTED
# quantity: no weighted proportion, no risk difference, no odds ratio, no tau, on any of the seven. The
# rule, stated once: a COUNT may appear here; a WEIGHTED QUANTITY may not. The regression pin is
# §15.0.2's golden vector, on a synthetic 24-record construction with no patient data in it.
#
# **The six silent failures this stage exists to make catchable**, each with the companion that shows
# what happens without the guard — and the first four are silent in a way Stage 8's were not, because
# they return a plausible number rather than a wrong shape::
#
#   ranges over in_model, not the outcome mask   15.2        tici 91 -> 92, minority 6 -> 7
#   augmentation tilted by `w` inside the        15.8 item 3 BOTH roadmap tests still pass
#     estimator
#   augmentation tilted by `w` AT THE CALL SITE  15.8 item 6 items 1-5 all still pass
#   `FitError` caught and downgraded             15.7a       every other test in §15 still passes
#   a `|` in an audit table cell                 15.14       byte-identity across seeds does NOT catch it
#   a `nan` proportion "corrected" to a finite   15.5        manufactures an odds ratio for an arm
#     odds ratio                                             with no data

from fixtures_stage9 import (DELTA, DR_SEED, DR_SEEDS_12, ato, dr_population,   # noqa: E402
                             golden_arrays, golden_frame, tilt_frame)
from fixtures_stage9 import separated_frame as separated_binary_frame           # noqa: E402

# Aliased for `test_model.py`'s reason: this module already imports Stage 8's separated ORDINAL frame
# from `test_model`, and Stage 9's is a separated BINARY frame. Two fitters, two frames, one word.


# --- 15.0  the frames this half is tested on ----------------------------------------------------------
#
# THE COHORT-LEVEL FIXTURE IS THIS FILE'S, AND THE SPEC SAYS SO. §22.3's closing paragraph records that
# round 3 "could not call `secondary` end to end: §15.0's fixtures are outcome-level, and building a
# cohort-level one means either the workbook or a fixture this document still does not specify. So the
# six ordering constraints, the three audit entries and S1-S8 remain specified and unexecuted — T7 and
# T9 are where they first run." This is T7 and T9. `secondary_cohort` below is that missing fixture.
#
# Four properties of it, each chosen against an alternative:
#
#   * DETERMINISTIC, NO RNG. Every column is a literal or a cycle over one, so the frame is a property
#     of this file rather than of a numpy version's random stream [Stage 7 §12.0.2].
#   * THREE CENTRES, USZ ABSENT. `model.design` therefore drops `center_USZ` from every design exactly
#     as it does on the workbook, and the fixture exercises the constant-column rule rather than
#     stepping around it (§3.3). §15.9's USZ companion is the frame that puts one back.
#   * THE OUTCOME MASK BITES. One record's `tici_2b_3` is missing, so its [§11] population is 39 against
#     the ATO's 40 — which is the condition Stage 8 §4.1 could not construct and under which an
#     implementation that never built the mask is GREEN (§4.2, §15.2).
#   * THE SEVEN MINORITY CELLS STRADDLE THE THRESHOLD. Three sit at or above RARE_MINORITY_THRESHOLD,
#     three below it, and `tici_2b_3` is below it WITH an override — so all three of §9.2's paths are
#     reached without the workbook.


def secondary_cohort(n: int = 40) -> pd.DataFrame:
    """A cohort-level frame `secondary` can run on: the nine [§6] covariates and all seven outcomes.

    Not `cohort_frame()` and not `ordinal_cohort()`: both build five records from the `HAND-4`
    template, and seven binary outcomes each needing both arms and a non-constant response do not fit
    in five records. This is built at 40, which is the size §21b's pass 4 used.

    The response patterns are periodic and their periods are chosen so that the minority cells land
    where §15.9 needs them — 20, 12, 4, 4, 8, 12 and 20 — rather than by accident. A frame on which
    every outcome took the same path would exercise one branch of §9.2 three times.

    **AND NO PATTERN MAY SHARE A PERIOD WITH THE TREATMENT**, which is the one constraint that is not
    cosmetic. `ivt` alternates on `i % 2`, so an outcome written as `i % 2 == 0` is PERFECTLY
    SEPARATED on the arm: every event in one arm and every non-event in the other, `p1 = 0`, `p0 = 1`,
    and every estimate on it is a degenerate-branch estimate. Two of the seven wanted a 20/20 split
    and the obvious way to write one is exactly that mistake, so the periods below are 4, 5, 10 and
    20 and every one of the seven is asserted to reach both arms in both cells.
    """
    center = ["HUG", "CHUV", "Lugano"] * (n // 3) + ["HUG"] * (n - 3 * (n // 3))
    onset = ["witnessed", "unwitnessed", "wake_up", "witnessed"] * (n // 4)
    frame = pd.DataFrame({
        "case_id": [f"SEC-{i:02d}" for i in range(n)],
        "center": center[:n],
        "onset_type": onset[:n],
        # THE PERIODS ARE THE SPECIFICATION AND THEY WERE CHOSEN BY MEASURING THE RANK, not by eye.
        # A design that is rank-deficient does not produce a bad estimate here — `model.firth` raises
        # F3 and no estimate exists at all, so every augmented outcome's test fails at once with a
        # message about collinearity. Three dependencies had to be broken to get to full rank, and all
        # three looked innocent: `sex = i % 2` IS the treatment; `prestroke_mrs = i % 3` is a linear
        # combination of the two centre dummies, because `center` also cycles on three; and
        # `atrial_fib = (i // 2) % 2` together with an `onset_type` of period four makes
        # `ivt - atrial_fib - onset_unwitnessed + onset_wake_up` exactly zero. Hence 5, 7 and 3.
        "age": [55.0 + (i % 7) * 4.0 for i in range(n)],
        "sex": [float((i // 5) % 2) for i in range(n)],
        "prestroke_mrs": [float((i // 7) % 3) for i in range(n)],
        "nihss_baseline": [6.0 + float(i % 11) for i in range(n)],
        "core_ml": [5.0 + (i % 5) * 3.0 for i in range(n)],
        "tmax6_ml": [40.0 + (i % 9) * 6.0 for i in range(n)],
        "atrial_fib": [float((i // 3) % 2) for i in range(n)],
        config.TREATMENT: [float(i % 2) for i in range(n)],
    })
    # Each entry is the predicate for an EVENT, keyed by outcome, and the comment is the minority
    # cell it produces at n = 40. Written out one per outcome rather than as a period table, because
    # the point of the frame is that the seven differ.
    events = {
        "mrs_0_2_90d": lambda i: (i // 2) % 2 == 0,     # 20 events / 20 non   -> minority 20, full
        "mrs_0_1_90d": lambda i: i % 10 < 3,            # 12 / 28              -> minority 12, full
        "tici_2b_3":   lambda i: i % 20 not in (4, 15),  # 36 / 4              -> minority  4, reduced
        "sich":        lambda i: i % 20 in (3, 8),      #  4 / 36              -> minority  4, unaug
        "ph2":         lambda i: i % 5 == 2,            #  8 / 32              -> minority  8, unaug
        "death_90d":   lambda i: (i + 5) % 10 < 3,      # 12 / 28              -> minority 12, full
        "mrs_5_6_90d": lambda i: (i // 5) % 2 == 0,     # 20 / 20              -> minority 20, full
    }
    for key in config.BINARY_OUTCOMES:
        frame[key] = [float(events[key](i)) for i in range(n)]
    # ONE record's TICI is absent, so the [§11] mask bites on this frame and `tici_2b_3`'s
    # denominator is 39 against the ATO's 40 — the condition Stage 8 §4.1 could not construct (§4.2)
    frame.loc[1, "tici_2b_3"] = np.nan
    return frame


def hand_propensity(df: pd.DataFrame, out: int = 0) -> propensity.Propensity:
    """A `Propensity` assembled by hand over `df`, with `out` records held off `in_model`.

    ASSEMBLED RATHER THAN FITTED, deliberately. §15.1a breaks `e`, `w` and `in_model` one at a time,
    and a fitted Propensity cannot be broken without also changing the estimate — so a fixture that
    ran `propensity.fit` would make every precondition frame test two things at once.

    `e` is a smooth deterministic function of the row position, bounded well inside (0, 1) so that
    `h = e(1-e)` never collapses; `w` is [§7]'s overlap weight, `1 - e` treated and `e` control, and
    BOTH carry `nan` off `in_model` because Stage 6 §9's rule is that a deliberate absence lives as a
    nan and is never filled. S4 is what asserts that, and a fixture that filled it would make S4
    untestable.
    """
    in_model = pd.Series(True, index=df.index)
    if out:
        in_model.iloc[-out:] = False
    e = pd.Series([0.25 + 0.012 * (i % 40) for i in range(len(df))], index=df.index, dtype=float)
    w = pd.Series(np.where(df[config.TREATMENT].to_numpy(dtype=float) == 1.0, 1.0 - e, e),
                  index=df.index, dtype=float)
    return propensity.Propensity(
        e=e.where(in_model), w=w.where(in_model), in_model=in_model,
        ess={code: float(w[in_model & (df[config.TREATMENT] == code)].sum())
             for code in config.TREATMENT_LABELS},
        fit=None, dropped=(),
        spec=config.PROPENSITY_PRIMARY)


def golden_propensity(df: pd.DataFrame) -> propensity.Propensity:
    """A `Propensity` over `golden_frame()`, built from the frame's OWN written-out `e` column.

    NOT `hand_propensity`, and the distinction is what §15.8 item 6 turns on. `hand_propensity`
    invents an `e`, which is right for the cohort-level fixture where no `e` is declared; the golden
    frame WRITES ITS `e` OUT, precisely so the golden vector is a property of this repository rather
    than of a fitted model. An estimate computed from an invented `e` is not the pinned one, and item
    6's whole assertion is that `secondary` returns exactly what `augmented_rd` returns on the golden
    arrays — so the two must come from the same column.
    """
    in_model = pd.Series(True, index=df.index)
    e = df["e"].astype(float)
    w = pd.Series(np.where(df[config.TREATMENT].to_numpy(dtype=float) == 1.0, 1.0 - e, e),
                  index=df.index, dtype=float)
    return propensity.Propensity(
        e=e, w=w, in_model=in_model,
        ess={code: float(w[df[config.TREATMENT] == code].sum())
             for code in config.TREATMENT_LABELS},
        fit=None, dropped=(),
        spec=config.PROPENSITY_PRIMARY)


def secondary_run(df=None, ps=None, **kwargs):
    """(df, ps, sec, audit) — the whole Stage 9 stage run once on the cohort-level fixture."""
    df = secondary_cohort() if df is None else df
    ps = hand_propensity(df) if ps is None else ps
    audit = data.Audit(hand_source(Path("."), "stage9"))
    return df, ps, outcome.secondary(df, ps, audit, **kwargs), audit


@pytest.fixture(scope="module")
def workbook_stage9():
    """(df, ps, sec, audit) from ONE linear run of the WHOLE pipeline against ONE FRESH `Audit`.

    NOT this file's `workbook` fixture, and the separation is deliberate rather than tidy. That one
    is module-scoped and stops at `propensity.fit`, so §14's tests share its `Audit` — and every
    Stage 9 data-gated test that called `outcome.secondary` on it would append three more entries to
    an object §14.10 reads. More to the point, §15.14's ledger assertion is about the FULL pipeline's
    30 entries, which needs `balance.assess` and `primary` to have run into the same object, and
    §14.0.3 records that `workbook` deliberately does not call `balance.assess`.

    One run, one Audit, every data-gated Stage 9 test reading it — so the ledger is counted once and
    the seven estimates are computed once.
    """
    df, audit = data.load(data.WORKBOOK)
    df = derive.derive(df, audit)
    df = eligibility.classify(df, audit)
    df = cohort.build(df, audit)
    ps = propensity.fit(df, audit)
    balance.assess(df, ps, audit)
    outcome.primary(df, ps, audit)
    return df, ps, outcome.secondary(df, ps, audit), audit


def as_outcome(frame: pd.DataFrame, key: str = "tici_2b_3") -> pd.DataFrame:
    """`frame` with its response column `y` renamed to a registry outcome key.

    `outcome_model(df, in_estimate, outcome)` reads `sub[outcome]` — it looks the response up BY THE
    REGISTRY KEY, which is what makes it unable to be pointed at a column the registry does not
    declare. §15.0's outcome-level fixtures call their response `y`, because they are arrays-with-a-
    frame for the estimators and carry no registry at all. This is the one-line adapter between them,
    and it is here rather than in `fixtures_stage9.py` so that the fixtures stay what the spec writes.

    `tici_2b_3` is the default because it is the one outcome whose covariate list is an override, so a
    test that forgets to say which outcome it means gets the REDUCED model — four columns — and a
    width assertion notices, where the shared nine-covariate list would fail for a reason about the
    fixture's columns instead.
    """
    return frame.rename(columns={"y": key})


STAGE_9_BANNER: Final[str] = "Stage 9 — the [§8] secondary binary estimators"


@pytest.fixture(scope="module")
def stage9_source() -> str:
    """`outcome.py` from the Stage 9 banner onward — the source every §15 scan runs against.

    SCOPED, and the scoping is the assertion. `outcome.py` holds two stages, and a scan over the whole
    file would report Stage 8's `_assert_reportable` and its `POLR_MAX_ABS_BETA` as Stage 9's, so
    §15.7's "no bound is read" and §15.7a's "nothing is caught" would both be measuring the wrong
    half. The split is on the banner rather than on a line number so that it survives an edit above it.
    """
    assert SOURCE.count(STAGE_9_BANNER) == 1, "the Stage 9 banner moved or was duplicated"
    return SOURCE.split(STAGE_9_BANNER, 1)[1]


@pytest.fixture(scope="module")
def stage9_names() -> set[str]:
    """Every identifier the Stage 9 half REFERENCES AS CODE — no docstrings, no comments, no prose.

    **A TEXT SCAN IS THE WRONG INSTRUMENT FOR "IS THIS CALLED", and this fixture exists because the
    first draft of §15.11 used one.** `outcome.py`'s Stage 9 half is heavily commented and its
    docstrings cite the neighbouring stages by name: they say `Primary` is not read, that `model.polr`
    is neither called nor edited, and that a deliberate absence is never filled. A `"polr" not in
    source` scan reports every one of those sentences as a violation of what the sentence itself
    promises — which is a test that fails precisely when the code is well documented.

    So the names come from the AST: `Name.id`, `Attribute.attr`, the dotted form of an attribute
    chain, and every function and class defined. String literals are deliberately EXCLUDED, which
    means a raise message quoting `weighted_proportion` is invisible here and correctly so.
    """
    module = ast.parse(SOURCE)
    boundary = SOURCE[:SOURCE.index(STAGE_9_BANNER)].count("\n") + 1
    names: set[str] = set()
    for node in module.body:
        if getattr(node, "lineno", 0) < boundary:
            continue
        for child in ast.walk(node):
            if isinstance(child, ast.Name):
                names.add(child.id)
            elif isinstance(child, ast.Attribute):
                names.add(child.attr)
                if isinstance(child.value, ast.Name):
                    names.add(f"{child.value.id}.{child.attr}")
            elif isinstance(child, (ast.FunctionDef, ast.ClassDef)):
                names.add(child.name)
    return names


def golden_secondary_arrays(covariates=("center", "atrial_fib")):
    """(y, a, w, h, m1, m0) on `golden_frame()` under one covariate list — §15.8's call-site oracle.

    Returned as a block because six arrays that must all come from ONE frame and ONE fit are six
    chances to pair a fit's counterfactuals with another frame's weights.
    """
    df = golden_frame()
    y, a, w, e = golden_arrays(df)
    X, _ = model.design(df, covariates)
    X = X.copy()
    X.insert(0, config.TREATMENT, a)
    fit = model.firth(X, y)
    m1, m0 = outcome._counterfactuals(fit, X)
    return y, a, w, outcome._tilt(e), m1, m0


# --- 15.1  the seven outcomes and the registry ---------------------------------------------------------

def test_secondary_estimates_exactly_BINARY_OUTCOMES_in_registry_order():
    """Computed from the registry, so an outcome added to `OUTCOMES` appears here without an edit.

    `list(...) == list(...)` and not a set comparison: the dict's iteration order IS the rendered
    log's row order, and §15.14 asserts the log is byte-identical across hash seeds because of it.
    """
    _, _, sec, _ = secondary_run()
    assert list(sec.estimates) == list(config.BINARY_OUTCOMES)
    assert all(sec.estimates[k].outcome == k for k in config.BINARY_OUTCOMES)


def test_an_EIGHTH_binary_column_the_registry_does_not_declare_is_estimated_for_NONE_of_it():
    """What the computed set PREVENTS. A frame may carry any number of 0/1 columns; only the ones
    the [§5] registry declares are outcomes, and `secondary` cannot be talked into a wider set by a
    frame."""
    df = secondary_cohort()
    df["mrs_0_3_90d"] = [float(i % 3 == 0) for i in range(len(df))]
    _, _, sec, _ = secondary_run(df)
    assert "mrs_0_3_90d" not in sec.estimates
    assert len(sec.estimates) == 7


def test_by_family_partitions_the_seven_into_exactly_the_TWO_declared_families():
    """[§13] applies Benjamini-Hochberg WITHIN each family, and a consumer that groups by reading
    `OUTCOMES[k].family` itself is a second place the partition is computed (§3.1)."""
    _, _, sec, _ = secondary_run()
    families = sec.by_family()
    assert set(families) == {"secondary", "safety"}
    assert sorted(e.outcome for group in families.values() for e in group) == sorted(
        config.BINARY_OUTCOMES)
    for family, group in families.items():
        # each tuple in C.BINARY_OUTCOMES order, which is what makes a rendered group reproducible
        assert [e.outcome for e in group] == [k for k in config.BINARY_OUTCOMES
                                              if config.OUTCOMES[k].family == family]


# --- 15.1a  the preconditions, ONE FRAME PER BRANCH ----------------------------------------------------
#
# Ten branches, ten frames, each derived from `secondary_cohort()` by breaking exactly one thing. Every
# one asserts a `SchemaError` WHOSE MESSAGE NAMES ITS OWN CONDITION — not merely that something raised
# — because a collected assertion that reports the wrong branch is worse than one that reports nothing.
#
# All ten are unreachable on the workbook (§4.4), which is precisely why the suite is the only thing
# that can reach them: they are written for [§10]'s frames, and [§10] does not exist yet.
#
# **THREE FRAMES HAVE TO BE BUILT MORE CAREFULLY THAN "BREAK ONE THING", and all three are the same
# trap**: an earlier phase's check firing on damage the frame did incidentally, so the branch under
# test is never reached. `S6` must blank `e` and `w` as well as `in_model`, or every row sits off the
# mask carrying a finite value and S4 fires in phase 2. `S7` must zero `w` only on rows that are IN
# `in_model`, or the blanket assignment overwrites the off-mask nan and S4 fires again. `S3a` must be a
# PERMUTATION and not a fresh index, or it is caught for the wrong reason.
#
# The general rule, since it will bite whoever adds an eleventh precondition: a frame testing phase k
# must be PRISTINE for phases 1 through k-1. Asserting the message rather than the exception type is
# the only assertion that catches this — `pytest.raises(SchemaError)` passes on all three mistakes.

BRANCH_PHASE: Final[dict[str, int]] = {
    "S1": 1, "S5a": 1, "S3a": 1, "S3b": 1, "S2": 2, "S4": 2, "S5b": 2, "S6": 3, "S7": 3, "S8": 3}

# **NINE OF THE TEN ARE `SchemaError` AND S8 IS `model.FitError`** — DECISION 6 (PI, 2026-08-25),
# Stage 10 §5.4. The map is here rather than as an `if branch == "S8"` inside the test because the
# split is the assertion: Stage 10 §15.7 requires that S8 raises `FitError` AND that it does not
# raise `SchemaError`, and that S6 and S7 are unmoved — "a change that moved all three would pass a
# test written only for S8". A table both parametrised tests read is what makes that checkable at a
# glance rather than by reading two branches.
BRANCH_ERROR: Final[dict[str, type[Exception]]] = {
    branch: (model.FitError if branch == "S8" else config.SchemaError) for branch in BRANCH_PHASE}


def broken_secondary(branch: str):
    """(df, ps) for one of §15.1a's ten branches. Exactly one thing is wrong with each."""
    df = secondary_cohort()
    ps = hand_propensity(df, out=2 if branch in ("S4",) else 0)

    if branch == "S1":
        return df.drop(columns=["mrs_0_1_90d"]), ps
    if branch == "S5a":
        return df.drop(columns=[config.TREATMENT]), ps
    if branch == "S3a":
        # A PERMUTATION: the same labels in the other order. Relabelling would be caught for the
        # wrong reason; this is the failure §12 measures, which moves the estimate without raising.
        return df, dataclasses.replace(ps, e=ps.e[::-1], w=ps.w[::-1])
    if branch == "S3b":
        return df, dataclasses.replace(ps, in_model=ps.in_model.astype("int64"))
    if branch == "S2":
        broken = df.copy()
        broken.loc[0, "mrs_0_1_90d"] = 2.0
        return broken, ps
    if branch == "S4":
        e = ps.e.copy()
        e.iloc[-1] = 0.5                       # a FINITE value on a row that is OFF in_model
        return df, dataclasses.replace(ps, e=e)
    if branch == "S5b":
        broken = df.copy()
        broken.loc[0, config.TREATMENT] = np.nan
        return broken, ps
    if branch == "S6":
        # blank e and w TOO, or every row sits off the mask carrying a finite value and S4 fires
        blank = pd.Series(np.nan, index=df.index, dtype=float)
        return df, dataclasses.replace(
            ps, in_model=pd.Series(False, index=df.index), e=blank, w=blank)
    if branch == "S7":
        # only on rows that are IN in_model, or the off-mask nan is overwritten and S4 fires
        treated = df[config.TREATMENT] == 1.0
        return df, dataclasses.replace(ps, w=ps.w.mask(treated & ps.in_model, 0.0))
    if branch == "S8":
        broken = df.copy()
        broken["sich"] = 0.0                   # constant on its own [§11] population
        return broken, ps
    raise AssertionError(f"no such branch: {branch}")


@pytest.mark.parametrize("branch", sorted(BRANCH_PHASE))
def test_each_precondition_branch_raises_ITS_OWN_ERROR_CLASS_NAMING_ITS_OWN_CONDITION(branch):
    """The message assertion is the test. Three of these ten frames were mis-specified in a way that
    only a message assertion catches — `pytest.raises(SchemaError)` passes on all three.

    The CLASS is `BRANCH_ERROR`'s and not `SchemaError` for all ten: S8 is `model.FitError` under
    DECISION 6 [Stage 10 §5.4], and the other nine are unmoved.
    """
    df, ps = broken_secondary(branch)
    with pytest.raises(BRANCH_ERROR[branch]) as e:
        outcome.secondary(df, ps, data.Audit(hand_source(Path("."), "stage9")))
    assert branch in message_of(e), (
        f"{branch}'s frame raised, but the message names a different branch:\n{message_of(e)}")


@pytest.mark.parametrize("branch,key", [("S1", "mrs_0_1_90d"), ("S2", "mrs_0_1_90d"),
                                        ("S6", "mrs_0_2_90d"), ("S8", "sich")])
def test_a_branch_about_ONE_OUTCOME_names_that_outcome(branch, key):
    """S6, S7 and S8 fire inside the loop, so their messages carry the key: "the outcome is constant"
    against seven outcomes is a message that costs a bisect. S1 and S2 name it for the same reason."""
    df, ps = broken_secondary(branch)
    with pytest.raises(BRANCH_ERROR[branch]) as e:
        outcome.secondary(df, ps, data.Audit(hand_source(Path("."), "stage9")))
    assert key in message_of(e)


def test_S7_names_the_ARM_that_carries_no_weight():
    df, ps = broken_secondary("S7")
    with pytest.raises(config.SchemaError) as e:
        outcome.secondary(df, ps, data.Audit(hand_source(Path("."), "stage9")))
    assert config.TREATMENT_LABELS[1] in message_of(e)


def test_the_two_FRAME_LEVEL_phases_COLLECT_and_report_every_failure_at_once():
    """A frame breaking S2, S4 and S5b at once produces ONE SchemaError naming all three.

    A phase that raised on the first failure would report one problem per run and make fixing a bad
    frame an iterative guessing game (Stage 6 §4.5). Asserted by counting the named conditions in the
    message rather than by asserting a substring, so a message that happened to mention S4 in prose
    would not pass.
    """
    df = secondary_cohort()
    ps = hand_propensity(df, out=2)
    df.loc[0, "mrs_0_1_90d"] = 2.0                     # S2
    df.loc[0, config.TREATMENT] = np.nan               # S5b
    e = ps.e.copy()
    e.iloc[-1] = 0.5                                   # S4
    with pytest.raises(config.SchemaError) as e_info:
        outcome.secondary(df, dataclasses.replace(ps, e=e),
                          data.Audit(hand_source(Path("."), "stage9")))
    message = message_of(e_info)
    assert all(name in message for name in ("S2", "S4", "S5b"))
    assert message.count("S2 ") == 1 and message.count("S5b ") == 1


@pytest.mark.parametrize("branch", ["S1", "S5a", "S3b"])
def test_WITHOUT_the_phase_1_RAISE_each_branch_produces_a_BARE_PANDAS_EXCEPTION(branch):
    """The companion form the Stage 7 learning prescribes, and it is what stops the split being
    refactored away by someone who sees two helpers where one would do (§4.4a).

    With `_assert_readable`'s raise removed, phase 2 runs on the broken frame and produces a pandas
    traceback in place of the message this stage wrote. `KeyError` for S1 and S5a — a column phase 1
    had just reported missing — and `KeyError` for S3b, because pandas reads a non-boolean Series as
    a sequence of LABELS. That last one is §4.4a's own defect one check further on, found by calling
    `secondary` end to end and not by reading.
    """
    df, ps = broken_secondary(branch)
    with pytest.raises(KeyError):
        outcome._assert_secondary_inputs(df, ps)       # phase 2, reached as if phase 1 had not raised


def test_WITHOUT_the_phase_1_raise_S3a_DOES_NOT_RAISE_AT_ALL_and_returns_a_different_number():
    """S3a's companion is not a bare exception, and that is the point (§12).

    `.loc[boolean_series]` returns rows in the SERIES' order while `m1` and `m0` come from
    `model.design(df.loc[mask])` in the FRAME's order, so a permuted-but-equal index pairs each
    record's weight with another record's fitted value: no raise, no nan, a different number. This is
    the one phase-1 branch whose absence is SILENT, which is why S3a is `index.equals` — order
    sensitive — rather than a length check or a set comparison.
    """
    df, ps = broken_secondary("S3a")
    outcome._assert_secondary_inputs(df, ps)                   # phase 2 passes it: values are fine
    correct = outcome.secondary(df, hand_propensity(df),
                                data.Audit(hand_source(Path("."), "stage9")))
    y = df["mrs_0_2_90d"].to_numpy(dtype=float)
    a = df[config.TREATMENT].to_numpy(dtype=float)
    permuted = outcome.weighted_rd(y, a, ps.w.to_numpy(dtype=float))
    assert np.isfinite(permuted)
    assert permuted != correct.estimates["mrs_0_2_90d"].rd


# --- 15.2  the seven [§11] denominators, and the one the workbook witnesses ----------------------------

def test_each_outcomes_population_is_in_model_AND_ITS_OWN_notna_boolean_and_TOTAL():
    df, ps, sec, _ = secondary_run()
    for key, est in sec.estimates.items():
        expected = ps.in_model & df[key].notna()
        assert est.in_estimate.dtype == bool
        assert not est.in_estimate.isna().any()
        assert est.in_estimate.index.equals(df.index)
        assert len(est.in_estimate) == len(df)
        assert list(est.in_estimate) == list(expected)


def test_every_population_is_a_SUBSET_of_in_model():
    """§12's a-fortiori claim about Stage 6 §9's rule, tested rather than merely argued."""
    _, ps, sec, _ = secondary_run()
    for est in sec.estimates.values():
        assert not (est.in_estimate & ~ps.in_model).any()


def test_SEVEN_DENOMINATORS_ARE_NOT_ONE_DENOMINATOR_and_the_fixture_makes_the_mask_bite():
    """The condition Stage 8 §4.1 recorded as the one under which a missing mask is GREEN.

    Stage 8 built the outcome mask, found the resulting population equal to `in_model`'s, and
    recorded that as "exactly the condition under which an implementation that never built the mask
    is green". Here one record's TICI is absent, so its denominator is 39 against the ATO's 40.
    """
    _, ps, sec, _ = secondary_run()
    sizes = {k: int(e.in_estimate.sum()) for k, e in sec.estimates.items()}
    assert sizes["tici_2b_3"] == int(ps.in_model.sum()) - 1
    assert len(set(sizes.values())) == 2


def test_AN_IMPLEMENTATION_RANGING_OVER_in_model_moves_the_denominator_AND_the_minority_cell():
    """The companion, and on this frame it FAILS on the wrong implementation rather than passing
    vacuously — which is what Stage 8 §14.2 recorded that no available frame could do.

    `y == 1` is False on a missing value (§3.3), so the absent record is silently counted as a
    NON-EVENT: the denominator moves 39 -> 40 and the minority cell 4 -> 5. Neither is a crash and
    neither is a nan; both are plausible numbers, and the minority cell is what §9.2 decides the
    estimator from.

    **The complement is what makes it silent, and it has to be written out.** `(y == 0).sum()` and
    `(~(y == 1)).sum()` agree on every frame with no missing outcome and differ by exactly the
    missing count on this one — the first excludes the nan, the second absorbs it. An implementation
    ranging over `in_model` has a nan in `y` and therefore reaches for a comparison; whichever of the
    two it writes, the WEIGHTED PROPORTIONS below are already wrong, because the indicator
    `(y == 1).astype(float)` puts a 0.0 where the outcome is absent and `np.average` weights it.
    """
    df, ps, sec, _ = secondary_run()
    est = sec.estimates["tici_2b_3"]
    assert (int(est.in_estimate.sum()), est.minority) == (39, 4)

    wrong = df.loc[ps.in_model, "tici_2b_3"].to_numpy(dtype=float)   # over in_model, not the mask
    assert len(wrong) == 40 and int(np.isnan(wrong).sum()) == 1
    assert min(int((wrong == 1.0).sum()), int((~(wrong == 1.0)).sum())) == 5

    # and the estimate itself moves, which is the harm the count is only a symptom of
    a = df.loc[ps.in_model, config.TREATMENT].to_numpy(dtype=float)
    w = ps.w[ps.in_model].to_numpy(dtype=float)
    wrong_rd = outcome.weighted_rd((wrong == 1.0).astype(float), a, w)
    assert np.isfinite(wrong_rd)
    assert wrong_rd != est.rd


@DATA_GATED
def test_the_workbook_witnesses_SIX_populations_of_92_and_ONE_of_91(workbook_stage9):
    """`[data-gated]`. A count, not a weighted quantity (§4.3).

    This is the assertion Stage 8 could not make against the workbook: `tici_2b_3` is present on 91
    of the 92, so the mask bites on v7 itself and not only on a construction.
    """
    df, ps, sec, audit = workbook_stage9
    sizes = {k: int(e.in_estimate.sum()) for k, e in sec.estimates.items()}
    assert len(df) == 93 and int(ps.in_model.sum()) == 92
    assert sizes["tici_2b_3"] == 91
    assert sorted(sizes.values()) == [91] + [92] * 6


@DATA_GATED
def test_the_workbook_minority_cells_are_the_seven_the_spec_measured(workbook_stage9):
    """`[data-gated]`, and these are COUNTS. §9.3's table, which is the only prespecified numerical
    claim the protocol makes about this workbook's outcome distribution."""
    df, ps, sec, audit = workbook_stage9
    assert {k: e.minority for k, e in sec.estimates.items()} == {
        "mrs_0_2_90d": 42, "mrs_0_1_90d": 24, "tici_2b_3": 6,
        "sich": 5, "ph2": 9, "death_90d": 25, "mrs_5_6_90d": 28}


# --- 15.3  the weighted risk difference, and the sign --------------------------------------------------
#
# Fixture 2 of §15.0: hand-built weighted frames of four to eight rows, with weights and outcomes chosen
# so the weighted proportions are EXACT DECIMALS — so the assertions below are against arithmetic done
# by hand and not against another implementation. Built inline, because each is three lines and naming
# them centrally would hide which arithmetic each checks.


def hand_weighted(y, a, w):
    return (np.array(y, dtype=float), np.array(a, dtype=float), np.array(w, dtype=float))


def test_weighted_rd_against_a_hand_computation():
    """Treated: weights 1 and 3 on outcomes 1 and 0 -> p1 = 1/4. Control: weights 2 and 2 on 1 and 0
    -> p0 = 1/2. RD = 0.25 - 0.5 = -0.25, exactly, with no rounding anywhere."""
    y, a, w = hand_weighted([1, 0, 1, 0], [1, 1, 0, 0], [1.0, 3.0, 2.0, 2.0])
    share = outcome.weighted_proportion(y, a, w)
    assert (share[1], share[0]) == (0.25, 0.5)
    assert outcome.weighted_rd(y, a, w) == -0.25


def test_weighted_rd_IS_the_two_proportions_differenced_and_nothing_else():
    """§5.1: `weighted_proportion` is called and not reimplemented, so the denominator [§11] requires
    has ONE definition. Asserted as an identity rather than trusted from the source."""
    y, a, w, _ = golden_arrays()
    share = outcome.weighted_proportion(y, a, w)
    assert outcome.weighted_rd(y, a, w) == share[1] - share[0]


def test_reversing_the_ARM_LABELS_negates_the_risk_difference():
    """Stage 8 §14.5's orientation move, on the other scale. `share[_TREATED] - share[_COMPARATOR]`
    and not `share[1] - share[0]`: this is one of the two lines in the module whose SIGN is a
    function of which literal goes first."""
    y, a, w = hand_weighted([1, 0, 1, 0], [1, 1, 0, 0], [1.0, 3.0, 2.0, 2.0])
    assert outcome.weighted_rd(y, 1.0 - a, w) == -outcome.weighted_rd(y, a, w)


def test_weighted_rd_is_nan_where_an_arm_carries_no_positive_weight_and_NOT_zero():
    """Inherited from `weighted_proportion` rather than caught. An arm with no positive weight makes
    the difference UNDEFINED, and that is what nan means; `0.0` would read as "no effect"."""
    y, a, w = hand_weighted([1, 0, 1, 0], [1, 1, 0, 0], [0.0, 0.0, 2.0, 2.0])
    assert np.isnan(outcome.weighted_rd(y, a, w))


def test_THE_SIGN_CONVENTION_IS_UNIFORM_ACROSS_ALL_SEVEN_and_higher_is_better_enters_no_arithmetic():
    """§4.1. Every risk difference is P(event | bridging) - P(event | EVT alone), on every outcome.

    The sign is the CONTRAST's, not the outcome's. Asserted by recomputing each of the seven from the
    two arms directly, in that order, and comparing to what `secondary` returned — so an
    implementation orienting one row by `higher_is_better` fails on the four outcomes for which the
    flag is False.
    """
    df, ps, sec, _ = secondary_run()
    assert len({config.OUTCOMES[k].higher_is_better for k in config.BINARY_OUTCOMES}) == 2
    for key, est in sec.estimates.items():
        sub = df.loc[est.in_estimate]
        y = sub[key].to_numpy(dtype=float)
        a = sub[config.TREATMENT].to_numpy(dtype=float)
        w = ps.w.loc[est.in_estimate].to_numpy(dtype=float)
        share = outcome.weighted_proportion(y, a, w)
        assert est.rd == share[1] - share[0]


def test_ORIENTING_BY_higher_is_better_would_put_TWO_MEANINGS_OF_POSITIVE_in_one_table():
    """The companion for §4.1's alternative-not-taken, and it is why the alternative is not taken.

    Flip the sign of every outcome for which lower is better and the seven-row table stops having one
    meaning: `mrs_0_2_90d` positive would mean "bridging helps" while `death_90d` positive would also
    mean "bridging helps" — the same word for two different contrasts, distinguishable only by a
    registry lookup a reader does not perform. Stage 14 [§16] is where a direction is attached.
    """
    _, _, sec, _ = secondary_run()
    oriented = {k: (e.rd if config.OUTCOMES[k].higher_is_better else -e.rd)
                for k, e in sec.estimates.items()}
    flipped = [k for k in config.BINARY_OUTCOMES if oriented[k] != sec.estimates[k].rd]
    assert sorted(flipped) == sorted(k for k in config.BINARY_OUTCOMES
                                     if not config.OUTCOMES[k].higher_is_better
                                     and sec.estimates[k].rd != 0.0)
    assert flipped, "the fixture must have at least one non-zero lower-is-better estimate"


# --- 15.4  the marginal odds ratio IS marginal ----------------------------------------------------------

def test_marginal_odds_ratio_against_a_hand_computation():
    """p1 = 1/4, p0 = 1/2 as above. OR = (0.25/0.75) / (0.5/0.5) = 1/3, exactly."""
    y, a, w = hand_weighted([1, 0, 1, 0], [1, 1, 0, 0], [1.0, 3.0, 2.0, 2.0])
    odds_ratio, corrected, share = outcome.marginal_odds_ratio(y, a, w)
    assert odds_ratio == pytest.approx(1.0 / 3.0, rel=1e-15)
    assert corrected is False
    assert (share[1], share[0]) == (0.25, 0.5)


def test_the_three_returned_values_are_CONSISTENT_whenever_the_correction_did_not_fire():
    """Recomputing the odds ratio from the returned proportions reproduces it EXACTLY. That is what
    makes the returned proportions evidence rather than decoration (§10.3)."""
    df, ps, sec, _ = secondary_run()
    for key, est in sec.estimates.items():
        assert est.or_corrected is False
        p1, p0 = est.proportion[1], est.proportion[0]
        assert est.odds_ratio == (p1 / (1.0 - p1)) / (p0 / (1.0 - p0))


def test_the_marginal_odds_ratio_DIFFERS_from_a_weighted_logistic_treatment_COEFFICIENT():
    """[§8]: "All estimates are marginal in the overlap population. No covariate adjustment of the
    reported effect."

    These are the two quantities most easily confused and the wrong one is one `firth` call away in
    this very module. The §8.2-shaped companion: assert the route taken, and COMPUTE the route not
    taken and show it gives a different number. A test that only asserted the marginal value would
    pass on an implementation that returned the conditional one on a frame where they happened to
    agree.
    """
    df = golden_frame()
    y, a, w, _ = golden_arrays(df)
    marginal, _, _ = outcome.marginal_odds_ratio(y, a, w)

    X, _ = model.design(df, ("center", "atrial_fib"))
    X = X.copy()
    X.insert(0, config.TREATMENT, a)
    conditional = float(np.exp(model.firth(X, y).beta[1]))          # beta[0] is the intercept

    assert np.isfinite(marginal) and np.isfinite(conditional)
    assert marginal != pytest.approx(conditional, rel=1e-3)


# --- 15.5  the degenerate cell, on all four of §6.2's rows ----------------------------------------------
#
# Fixture 3 of §15.0: seven frames, one per branch, built inline. **The arm weight totals are the point
# on two of them** and are therefore specified rather than left to whoever writes the test:
#
#   #  frame                                          Sw treated / control   reaches
#   1  no events in the treated arm                        even, ~4 / ~4     p1 == 0
#   2  no non-events in the treated arm                    even, ~4 / ~4     p1 == 1
#   3  no events in the control arm                        even, ~4 / ~4     p0 == 0
#   4  no non-events in the control arm                    even, ~4 / ~4     p0 == 1
#   5  zero total weight in one arm                        0 / ~4            p is nan, NOT corrected
#   6  no events in the treated arm, treated arm HEAVIER   17.007 / 22.768   the NULL CROSSING (§6.4)
#   7  outcome constant on the population                  even, ~4 / ~4     finite OR, not nan (§6.2)
#
# Frames 1-4 and 7 need only that the arms be non-degenerate in weight, so their totals are
# unconstrained and are "even" to make clear that nothing rests on them. **Frame 6's are not free**: the
# crossing needs the arm holding the empty cell to carry MORE total weight than the other, and the
# effect is a function of how much more. At the workbook's near-even split the corrected odds ratio is
# 0.6484 and does NOT cross; at 17.007 / 22.768 it is 1.0201 and does. A frame built with round weights
# would exercise the branch and miss the property.

DEGENERATE = {
    #                y                a                w
    "p1==0": ([0, 0, 1, 0], [1, 1, 0, 0], [2.0, 2.0, 2.0, 2.0]),
    "p1==1": ([1, 1, 1, 0], [1, 1, 0, 0], [2.0, 2.0, 2.0, 2.0]),
    "p0==0": ([1, 0, 0, 0], [1, 1, 0, 0], [2.0, 2.0, 2.0, 2.0]),
    "p0==1": ([1, 0, 1, 1], [1, 1, 0, 0], [2.0, 2.0, 2.0, 2.0]),
}


@pytest.mark.parametrize("branch", sorted(DEGENERATE))
def test_a_degenerate_proportion_is_CORRECTED_finite_and_FLAGGED(branch):
    y, a, w = hand_weighted(*DEGENERATE[branch])
    odds_ratio, corrected, share = outcome.marginal_odds_ratio(y, a, w)
    assert corrected is True
    assert np.isfinite(odds_ratio) and odds_ratio > 0.0
    # the two RAW proportions travel with it, so the correction is never the only record it fired
    assert 0.0 in (share[1], share[0]) or 1.0 in (share[1], share[0])


def test_an_INTERIOR_estimate_is_BIT_FOR_BIT_UNCORRECTED():
    """"Only when it fires" is asserted rather than assumed: the corrected and uncorrected paths are
    compared on a frame where both are computable and shown to agree at 0.000e+00."""
    y, a, w = hand_weighted([1, 0, 1, 0], [1, 1, 0, 0], [1.0, 3.0, 2.0, 2.0])
    odds_ratio, corrected, share = outcome.marginal_odds_ratio(y, a, w)
    p1, p0 = share[1], share[0]
    assert corrected is False
    assert odds_ratio - (p1 / (1.0 - p1)) / (p0 / (1.0 - p0)) == 0.0


def test_a_nan_PROPORTION_IS_NOT_CORRECTED():
    """§6.2's fourth row. A missing proportion is not a boundary proportion: there is no cell to add
    half of anything to, and a correction there would manufacture an odds ratio for an arm with no
    data."""
    y, a, w = hand_weighted([1, 0, 1, 0], [1, 1, 0, 0], [0.0, 0.0, 2.0, 2.0])
    odds_ratio, corrected, share = outcome.marginal_odds_ratio(y, a, w)
    assert np.isnan(odds_ratio)
    assert corrected is False
    assert np.isnan(share[1])


def test_A_WRONG_GUARD_BRANCHING_ON_not_0_lt_p_lt_1_MANUFACTURES_AN_ODDS_RATIO_FROM_NO_DATA():
    """The companion, and it is the reason the `isnan` test comes FIRST.

    `not (0 < p < 1)` is True for nan, so a guard written that way takes the degenerate branch — and
    the degenerate branch does not divide by p at all. It adds OR_CONTINUITY to four weight-sums, and
    an arm carrying no positive weight has weight-sums of exactly zero, so it returns 1.0: a finite,
    plausible, null odds ratio for an arm with no data in it.
    """
    y, a, w = hand_weighted([1, 0, 1, 0], [1, 1, 0, 0], [0.0, 0.0, 2.0, 2.0])
    share = outcome.weighted_proportion(y, a, w)
    p1, p0 = share[1], share[0]
    assert not (0.0 < p1 < 1.0)                        # nan satisfies the wrong guard
    manufactured = {}
    for code in config.TREATMENT_LABELS:
        arm = a == float(code)
        manufactured[code] = ((float(w[arm][y[arm] == 1.0].sum()) + config.OR_CONTINUITY)
                              / (float(w[arm][y[arm] == 0.0].sum()) + config.OR_CONTINUITY))
    wrong = manufactured[1] / manufactured[0]
    assert np.isfinite(wrong) and wrong == 1.0
    assert np.isnan(outcome.marginal_odds_ratio(y, a, w)[0])       # what the code actually returns


def test_a_CONSTANT_OUTCOME_returns_a_FINITE_number_and_not_nan():
    """§6.2's corrected claim, and it is the frame §16 item 3 should be decided on.

    An earlier draft argued S8 from "`0/0` gives nan silently". The degenerate branch fires BEFORE
    any division, so `0/0` is unreachable and a constant outcome comes back as a plausible number
    near 1 — which is WORSE than a nan a reader would notice, and makes the case for S8 raising
    stronger than the case originally written for it.

    **THE ARM WEIGHT TOTALS ARE 8 AND 11 AND THAT IS NOT FREE**, although §15.0.3's table lists this
    frame's totals as "even, ~4 / ~4". The two are inconsistent and §6.2's pins are what settle it:
    on a constant outcome every event pseudo-count is `OR_CONTINUITY` alone, so the whole odds ratio
    is `(Sw0 + c) / (Sw1 + c)` — a pure function of the two arm totals. Even totals give exactly 1.0,
    which is a degenerate special case that would hide the asymmetry §6.4 is about. `Sw1 = 8` and
    `Sw0 = 11` give `11.5 / 8.5 = 1.3529411765` and its reciprocal `0.7391304348`, which are §6.2's
    two measured values — and the fact that they ARE reciprocals is itself the statement that the
    two directions of degeneracy are the same arithmetic.
    """
    all_zero = hand_weighted([0, 0, 0, 0], [1, 1, 0, 0], [3.0, 5.0, 4.0, 7.0])
    all_one = hand_weighted([1, 1, 1, 1], [1, 1, 0, 0], [3.0, 5.0, 4.0, 7.0])
    for arrays in (all_zero, all_one):
        odds_ratio, corrected, _ = outcome.marginal_odds_ratio(*arrays)
        assert corrected is True
        assert np.isfinite(odds_ratio) and 0.5 < odds_ratio < 2.0     # plausible, and near the null
    assert outcome.marginal_odds_ratio(*all_zero)[0] == pytest.approx(1.3529411765, abs=1e-9)
    assert outcome.marginal_odds_ratio(*all_one)[0] == pytest.approx(0.7391304348, abs=1e-9)
    # reciprocals, because on a constant outcome the whole ratio is (Sw0 + c) / (Sw1 + c) either way
    assert (outcome.marginal_odds_ratio(*all_zero)[0]
            * outcome.marginal_odds_ratio(*all_one)[0]) == pytest.approx(1.0, abs=1e-12)


def test_THE_CORRECTION_CAN_CARRY_THE_ODDS_RATIO_ACROSS_THE_NULL_and_this_asserts_that_it_does():
    """§6.4's fourth measurement, and this bullet asserts the OPPOSITE of what an earlier draft did.

    That draft said the corrected odds ratio "is strictly between 1 and the uncorrected limit — i.e.
    it shrinks toward the null", which is false and would have been a test that fails once someone
    built the frame for it. What is guaranteed is only that the result is FINITE AND POSITIVE.

    The two frames below are §6.4's own two rows. Both have an empty treated EVENT cell, so the
    uncorrected odds ratio is 0.0 in each. The near-even one — the workbook's own 13.626 / 14.110
    split — lands at 0.6484, on the same side of 1. The one where the empty arm is the HEAVIER lands
    at 1.0201, ACROSS 1: a safety outcome with no bridging events at all, reported as favouring EVT
    alone. Measured, about 0.9% of empty-cell replicates cross, which is roughly 3 in an N_BOOT of
    2000. It is a small number and it is the wrong direction on the wrong family, so the suite
    records it rather than Stage 10 discovering it inside a percentile interval.
    """
    def empty_treated_event_cell(sw_treated, sw_control, p0):
        """One treated non-event carrying the whole treated weight; the control arm split to give p0."""
        y = np.array([0.0, 1.0, 0.0])
        a = np.array([1.0, 0.0, 0.0])
        w = np.array([sw_treated, sw_control * p0, sw_control * (1.0 - p0)])
        return y, a, w

    near_even = empty_treated_event_cell(13.626360, 14.110263, 0.02)
    heavier = empty_treated_event_cell(17.007, 22.768, 0.00647)

    for arrays in (near_even, heavier):
        share = outcome.weighted_proportion(*arrays)
        assert share[1] == 0.0                                  # the empty cell, in both
        assert (share[1] / (1.0 - share[1])) / (share[0] / (1.0 - share[0])) == 0.0   # uncorrected

    even_or, even_corrected, _ = outcome.marginal_odds_ratio(*near_even)
    heavy_or, heavy_corrected, _ = outcome.marginal_odds_ratio(*heavier)
    assert even_corrected is True and heavy_corrected is True
    assert np.isfinite(even_or) and np.isfinite(heavy_or)

    assert even_or == pytest.approx(0.6484, abs=1e-3)
    assert even_or < 1.0                                        # does not cross
    assert heavy_or == pytest.approx(1.0201, abs=1e-3)
    assert heavy_or > 1.0                                       # DOES cross — the fact being recorded


def test_the_crossing_needs_the_EMPTY_ARM_TO_BE_THE_HEAVIER_and_that_is_asserted_not_assumed():
    """The mechanism, so that the number above is explained by the suite and not only by §6.4.

    The four pseudo-counts are WEIGHT-SUMS, so the same additive OR_CONTINUITY is a different
    RELATIVE nudge in each arm. Swap which arm is heavier, holding everything else, and the crossing
    disappears — which is what makes this a property of the asymmetry and not of the constant.
    """
    def empty_treated_event_cell(sw_treated, sw_control, p0):
        return (np.array([0.0, 1.0, 0.0]), np.array([1.0, 0.0, 0.0]),
                np.array([sw_treated, sw_control * p0, sw_control * (1.0 - p0)]))

    heavier, _, _ = outcome.marginal_odds_ratio(*empty_treated_event_cell(22.768, 17.007, 0.00647))
    lighter, _, _ = outcome.marginal_odds_ratio(*empty_treated_event_cell(17.007, 22.768, 0.00647))
    assert lighter > 1.0 > heavier


@DATA_GATED
def test_the_correction_fires_for_NONE_of_the_seven_on_the_workbook(workbook_stage9):
    """`[data-gated]`. Stage 8 §4.1's condition again — an implementation that never wrote the branch
    is green on the data — which is why §15.5's coverage comes from constructed frames alone."""
    df, ps, sec, audit = workbook_stage9
    assert not any(e.or_corrected for e in sec.estimates.values())


# --- 15.7 (outcome's half)  m_a(X) and the counterfactuals ----------------------------------------------
#
# `predict`'s own section and the separated-frame witness are in `test_model.py`, because that is where
# `predict` and `firth` are. What is here is what `outcome.py` adds on top of them: the design that
# `outcome_model` assembles, and the two counterfactual frames `_counterfactuals` builds from it.

def test_outcome_model_inserts_the_treatment_main_effect_FIRST_and_BY_NAME():
    """[§8] calls treatment a MAIN EFFECT ADDED TO the [§6] covariate list, and `design` never returns
    it — the list does not contain it — so inserting it here is the whole of "treatment main effect".
    """
    df = as_outcome(golden_frame())
    in_estimate = pd.Series(True, index=df.index)
    fit, X, dropped = outcome.outcome_model(df, in_estimate, "tici_2b_3")
    assert X.columns[0] == config.TREATMENT
    assert config.TREATMENT not in config.outcome_model_covariates("tici_2b_3")
    assert fit.columns == tuple(X.columns)
    assert dropped == ("center_USZ",)                  # USZ is absent from the fixture, as from v7


def test_outcome_model_reads_the_ACCESSOR_and_never_OUTCOME_COVARIATES_directly():
    """The cross-reference, rather than a second list, is what keeps the propensity and outcome
    covariate sets from drifting apart [§8]. Asserted by the widths the two lists produce."""
    df = as_outcome(golden_frame())
    in_estimate = pd.Series(True, index=df.index)
    _, reduced, _ = outcome.outcome_model(df, in_estimate, "tici_2b_3")       # the override
    assert tuple(reduced.columns) == (config.TREATMENT, "atrial_fib",
                                      "center_CHUV", "center_Lugano")
    assert config.outcome_model_covariates("tici_2b_3") == ("center", "atrial_fib")
    assert config.outcome_model_covariates("sich") is config.PS_COVARIATES


def test_the_counterfactuals_differ_in_EXACTLY_the_treatment_terms_contribution():
    """[§8] prescribes no interactions, so `logit(m1) - logit(m0)` is CONSTANT across rows and equals
    the treatment coefficient. Asserted rather than assumed — this is the no-interaction property as
    an executable claim, and an implementation that rebuilt the frame rather than copying it, or that
    read `Fit.p`, does not have it.
    """
    df = as_outcome(golden_frame())
    fit, X, _ = outcome.outcome_model(df, pd.Series(True, index=df.index), "tici_2b_3")
    m1, m0 = outcome._counterfactuals(fit, X)
    gap = np.log(m1 / (1.0 - m1)) - np.log(m0 / (1.0 - m0))
    assert gap == pytest.approx(np.full(len(df), fit.beta[1]), abs=1e-9)
    assert float(np.ptp(gap)) < 1e-9


def test_every_NON_TREATMENT_column_keeps_its_observed_value_in_both_counterfactuals():
    """[§8]'s m_a(X) is the conditional mean given the covariates AT their observed values. Only the
    exposure is intervened on, and this asserts the other columns were not touched."""
    df = as_outcome(golden_frame())
    fit, X, _ = outcome.outcome_model(df, pd.Series(True, index=df.index), "tici_2b_3")
    before = X.copy()
    outcome._counterfactuals(fit, X)
    pd.testing.assert_frame_equal(X, before)           # and the design itself is not mutated


def test_Fit_p_IS_NEITHER_counterfactual_and_equals_each_only_on_its_OWN_arm():
    """§7.3's whole reason. `Fit.p` is m_A(X) — the counterfactual matching each row's ACTUAL arm —
    so it agrees with `m1` on the treated rows and with `m0` on the comparator rows, and with
    neither everywhere. An implementation reading it for either is wrong on every row assigned to the
    other arm, and the result is a plausible number.
    """
    df = as_outcome(golden_frame())
    a = df[config.TREATMENT].to_numpy(dtype=float)
    fit, X, _ = outcome.outcome_model(df, pd.Series(True, index=df.index), "tici_2b_3")
    m1, m0 = outcome._counterfactuals(fit, X)

    assert np.array_equal(fit.p[a == 1.0], m1[a == 1.0])
    assert np.array_equal(fit.p[a == 0.0], m0[a == 0.0])
    assert not np.array_equal(fit.p, m1) and not np.array_equal(fit.p, m0)
    assert np.all(np.isfinite(fit.p))                  # the wrong reading is finite, hence silent


def test_REBUILDING_the_frame_from_the_COVARIATE_LIST_instead_of_copying_the_design_RAISES():
    """The companion (§7.3). `design` dropped `center_USZ` and the covariate list does not know it,
    so a rebuilt frame has a different width — and `predict`'s by-name check is the only thing
    between that and a positional dot product that returns a number."""
    df = as_outcome(golden_frame())
    fit, X, dropped = outcome.outcome_model(df, pd.Series(True, index=df.index), "tici_2b_3")
    assert dropped == ("center_USZ",)

    rebuilt = pd.get_dummies(
        df["center"].astype(pd.CategoricalDtype(config.FACTOR_LEVELS["center"])),
        prefix="center", drop_first=True, dtype=float)
    rebuilt.insert(0, "atrial_fib", df["atrial_fib"].to_numpy(dtype=float))
    rebuilt.insert(0, config.TREATMENT, df[config.TREATMENT].to_numpy(dtype=float))
    assert rebuilt.shape[1] == X.shape[1] + 1          # the dropped column is back
    with pytest.raises(config.SchemaError):
        model.predict(fit, rebuilt)


def test_NO_BOUND_ON_THE_NUISANCE_COEFFICIENT_IS_READ_ANYWHERE_IN_STAGE_9(stage9_source):
    """§9.6, asserted as the deliberate absence it is, so that adding one becomes a visible change.

    Stage 8's half of this module reads `POLR_MAX_ABS_BETA` in `_assert_reportable` and raises
    `FitError` above it, because there `beta` IS the estimand. Stage 9's half reads no bound at all,
    for two reasons: there is no empty band to calibrate one from — `sich`'s replicate max_abs_beta
    runs continuously from 0.8558 to 80.9825 with 42.7% above 8 — and `m_a(X)` is a NUISANCE whose
    only contribution to `tau` is its PREDICTIONS, which are probabilities bounded in [0, 1] whatever
    beta does. Every term of `augmented_rd` is therefore bounded and there is no `exp(beta)`-style
    tail to protect against.

    What is lost instead is stated in §9.6 and is not an infinity: a separated `m_a` is an over-fitted
    one, so the correction term removes signal rather than residual confounding. §15.0.2's two golden
    pins are that cost as a measurement — the same 24 records, one extra covariate, and tau collapses
    seven-fold while `max_abs_beta` moves from 1.59 to 11.83.
    """
    assert "MAX_ABS_BETA" not in stage9_source
    assert "_assert_reportable" not in stage9_source
    assert "POLR_MAX_ABS_BETA" in SOURCE               # Stage 8's exists, so this is a contrast


def test_THE_SEPARATED_FRAME_PRODUCES_A_LARGE_max_abs_beta_THAT_NOTHING_REJECTS():
    """The consequence of the section above, on the frame `test_model.py` measures the fit on.

    `firth` returns; nothing in `outcome.py` looks at the coefficient; and the only thing standing
    between this nuisance model and an unremarked estimate is `max_abs_beta` appearing in §10.4's
    rendered table, which §15.14 asserts.
    """
    frame = as_outcome(separated_binary_frame())
    fit, X, _ = outcome.outcome_model(frame, pd.Series(True, index=frame.index), "tici_2b_3")
    assert float(np.max(np.abs(fit.beta))) == pytest.approx(6.089057, abs=1e-6)
    m1, m0 = outcome._counterfactuals(fit, X)
    assert np.all((m1 > 0.0) & (m1 < 1.0)) and np.all((m0 > 0.0) & (m0 < 1.0))


# --- 15.7a  the `FitError` contract: it propagates, and it is NEVER substituted --------------------------
#
# §9.5's no-substitution rule is the thing this stage argues hardest for, and [§10] states it in terms:
# "Replicates whose prespecified fit fails are dropped and counted; never substituted with a different
# estimator." Stage 8 tested the same contract for its own fitter one stage earlier; this is that shape
# ported across.


def unfittable_cohort() -> pd.DataFrame:
    """A cohort-level frame on which ONE augmented outcome's `m_a(X)` cannot be fitted.

    The failure is F4, a non-finite value in the design, reached by putting a nan in a covariate on a
    row that IS in the estimation population. It is `model.firth`'s to raise and Stage 9 adds nothing
    to it — which is the whole point: `secondary` has no handler and the exception leaves the stage.
    """
    df = secondary_cohort()
    df.loc[0, "core_ml"] = np.nan
    return df


def test_secondary_RAISES_FitError_rather_than_returning_a_quietly_unaugmented_estimate():
    df = unfittable_cohort()
    with pytest.raises(model.FitError):
        outcome.secondary(df, hand_propensity(df), data.Audit(hand_source(Path("."), "stage9")))


def test_FitError_IS_NOT_A_SchemaError_and_Stage_10s_whole_drop_and_count_rule_needs_that():
    """Stage 8 §11's handover as an assertion: `FitError` is the droppable failure and `SchemaError`
    is NOT — a SchemaError from `secondary` or `propensity.fit` is a bug in the resampler, not a
    sparse replicate. The rule is wrong if the two are ever related by inheritance in either
    direction."""
    assert not issubclass(model.FitError, config.SchemaError)
    assert not issubclass(config.SchemaError, model.FitError)


def test_THE_ENTRY_NAMING_THE_POPULATIONS_SURVIVES_A_FitError_which_is_why_it_is_recorded_FIRST():
    """§10.2's ordering rule, and the failing path is where it earns its place.

    `propensity._record_exclusion` states the principle: "the exclusion is a fact about the frame,
    established before any modelling choice, and a log that names it is useful precisely when what
    follows fails." Here what follows fails, and the log still names the seven populations.
    """
    df = unfittable_cohort()
    audit = data.Audit(hand_source(Path("."), "stage9"))
    with pytest.raises(model.FitError):
        outcome.secondary(df, hand_propensity(df), audit)
    assert [e.step for e in audit.entries] == ["binary_estimable"]
    assert len(audit.entries[0].table) == 1 + len(config.BINARY_OUTCOMES)


def test_THE_try_except_COMPANION_returns_seven_estimates_and_EVERY_OTHER_TEST_STILL_PASSES():
    """The companion is what makes §15.7a a test rather than an assertion (§9.5).

    With `outcome_model`'s call wrapped in `try/except model.FitError` and the path downgraded to
    "unaugmented", `secondary` returns a complete `Secondary` — seven estimates, one of them silently
    on a DIFFERENT ESTIMATOR — and nothing else in §15 notices. Every field of the downgraded
    estimate is internally consistent: the block of four is `None`/`None`/`None`/`()` exactly as
    §3.1 requires of an unaugmented outcome, `augmented_path` reads "unaugmented", and §15.9's
    block-move assertion passes on it. That is the point being recorded: the most natural thing an
    implementer writes when a fit fails on 0.3% of `sich` replicates is invisible to everything here
    except this test.
    """
    df = unfittable_cohort()
    ps = hand_propensity(df)
    audit = data.Audit(hand_source(Path("."), "stage9"))

    # the downgrade, transcribed from §11 with one `try` added — which is the whole of the defect
    estimates: dict[str, outcome.BinaryEstimate] = {}
    for key in config.BINARY_OUTCOMES:
        in_estimate = ps.in_model & df[key].notna()
        sub = df.loc[in_estimate]
        y = sub[key].to_numpy(dtype=float)
        a = sub[config.TREATMENT].to_numpy(dtype=float)
        w = ps.w.loc[in_estimate].to_numpy(dtype=float)
        e = ps.e.loc[in_estimate].to_numpy(dtype=float)
        minority = min(int((y == 1.0).sum()), int((y == 0.0).sum()))
        odds_ratio, corrected, share = outcome.marginal_odds_ratio(y, a, w)
        path = outcome._augmented_path(key, minority)
        tau = fit = covariates = None
        dropped: tuple[str, ...] = ()
        if path != "unaugmented":
            try:
                covariates = config.outcome_model_covariates(key)
                fit, X, dropped = outcome.outcome_model(df, in_estimate, key)
                m1, m0 = outcome._counterfactuals(fit, X)
                tau = outcome.augmented_rd(y, a, w, outcome._tilt(e), m1, m0)
            except model.FitError:                     # <-- the substitution [§10] forbids
                path, tau, fit, covariates, dropped = "unaugmented", None, None, None, ()
        estimates[key] = outcome.BinaryEstimate(
            outcome=key, family=config.OUTCOMES[key].family, minority=minority,
            rd=outcome.weighted_rd(y, a, w), odds_ratio=odds_ratio, or_corrected=corrected,
            proportion=dict(share), augmented_path=path, augmented=tau, covariates=covariates,
            reduced=key in config.OUTCOME_MODEL_OVERRIDES, dropped=dropped,
            in_estimate=in_estimate, fit=fit)
    caught = outcome.Secondary(estimates=estimates)

    # it returns, and it returns something that looks entirely correct
    assert len(caught.estimates) == 7
    downgraded = [k for k, e in caught.estimates.items()
                  if e.augmented_path == "unaugmented"
                  and outcome._augmented_path(k, e.minority) != "unaugmented"]
    assert downgraded, "the fixture must actually make one augmented fit fail"
    for key in downgraded:
        est = caught.estimates[key]
        # §3.1's block of four is internally consistent, so no structural test can see the swap
        assert (est.augmented, est.covariates, est.fit, est.dropped) == (None, None, None, ())
        assert np.isfinite(est.rd) and np.isfinite(est.odds_ratio)
    # and the shipped function does NOT do this
    with pytest.raises(model.FitError):
        outcome.secondary(df, ps, audit)


def test_secondary_catches_EXACTLY_ONE_THING_and_it_is_model_FitError(stage9_source):
    """By scan, because the behavioural test above can only reach the failure modes it constructs.

    **This assertion was `handlers == []` until Stage 10 landed `collect`**, and the amendment is
    licensed by Stage 10 §13 rather than by this file: [§10] left the drop granularity open, DECISION
    7 (PI, 2026-08-25) chose per-outcome replicate sets, and Stage 10 §7.3 puts the mechanism in
    Stage 9's own seven-outcome loop rather than re-implementing that loop in `bootstrap.py` — which
    Stage 9 §14 names as how §6.3's frozen paths get violated by accident.

    What the scan asserts now is the shape the original assertion was protecting:

      * there is EXACTLY ONE handler, so no second `try` has appeared;
      * it names `model.FitError` and not `Exception` and not `C.SchemaError`, which is Stage 10
        §7.1's taxonomy — a bare `except Exception` is `pilots/analysis.py:595`'s defect, where a
        `KeyError` from a typo and an unfittable model are counted as the same dropped replicate;
      * its body contains a bare `raise`, which is `collect=False` being the point-estimate contract
        unchanged rather than a default nobody re-checks.
    """
    tree = ast.parse(stage9_source)
    handlers = [node for node in ast.walk(tree) if isinstance(node, ast.ExceptHandler)]
    assert len(handlers) == 1
    caught = handlers[0].type
    assert isinstance(caught, ast.Attribute) and caught.attr == "FitError"
    assert isinstance(caught.value, ast.Name) and caught.value.id == "model"
    reraises = [n for n in ast.walk(handlers[0]) if isinstance(n, ast.Raise) and n.exc is None]
    assert len(reraises) == 1


@DATA_GATED
def test_no_outcome_raises_FitError_on_the_workbook_which_is_why_the_frame_is_constructed(workbook_stage9):
    """`[data-gated]`. Measured: `model.firth` raises on 0.3% of `sich` replicates and 0% elsewhere,
    so the point estimate never reaches it and only a construction can."""
    df, ps, sec, audit = workbook_stage9
    assert sum(1 for e in sec.estimates.values() if e.fit is not None) == 5


# --- 15.8  the augmentation: the roadmap's TWO acceptance tests become THREE, in SEVEN assertions -------
#
# The heading counts three ACCEPTANCE TESTS and seven ASSERTIONS, because they are different things.
# Test 1 is items 1-2, test 2 is items 3 and 6, and the third acceptance test is §15.10. Items 4, 5 and
# 7 are this stage's own.
#
# **THE ROADMAP'S FIRST ACCEPTANCE TEST IS RIGHT AND ITS STATED REASON FOR IT WAS WRONG.** It read:
# "with the outcome model set to a constant, the augmented estimate equals the unaugmented weighted
# risk difference exactly. If it does not, the augmentation is not built on the weights' tilting
# function and is targeting a different estimand." The identity holds; the inference does not. For
# constants m1 = c1 and m0 = c0 and ANY normalised weight u,
#
#     Σu(c1 − c0)/Σu  =  (c1 − c0)·Σu/Σu  =  c1 − c0
#
# — the correction term does not depend on `u` at all. So the identity holds whether the augmentation
# tilts by `h`, by `w`, or by unit weights, and passing it establishes NOTHING about the tilting
# function. What it does check is that the correction term is NORMALISED and CANCELS, which is a real
# property and the one an implementation summing Σ(m1 − m0) without a denominator gets wrong.


def test_item_1_TWO_DIFFERENT_CONSTANTS_reduce_the_augmented_estimate_to_the_unaugmented_one():
    """Acceptance test 1, with the corrected rationale. `m1 = 0.62`, `m0 = 0.23`."""
    y, a, w, e = golden_arrays()
    h = outcome._tilt(e)
    m1, m0 = np.full(len(y), 0.62), np.full(len(y), 0.23)
    assert outcome.augmented_rd(y, a, w, h, m1, m0) == pytest.approx(
        outcome.weighted_rd(y, a, w), abs=1e-15)


def test_item_1_THE_COMPANION_dropping_the_correction_term_returns_rd_MINUS_the_gap():
    """And this is what makes item 1 a test. An implementation omitting the correction term entirely
    comes back at `rd - 0.390000`, which is exactly `c0 - c1` — a decisive failure."""
    y, a, w, e = golden_arrays()
    m1, m0 = np.full(len(y), 0.62), np.full(len(y), 0.23)
    t, c = a == 1.0, a == 0.0
    dropped_term = (float(np.sum(w[t] * (y[t] - m1[t])) / np.sum(w[t]))
                    - float(np.sum(w[c] * (y[c] - m0[c])) / np.sum(w[c])))
    assert dropped_term - outcome.weighted_rd(y, a, w) == pytest.approx(0.23 - 0.62, abs=1e-12)


def test_item_2_ONE_CONSTANT_is_kept_as_a_companion_and_ASSERTED_TO_BE_INSUFFICIENT():
    """The roadmap's phrasing does not say the constants must differ, and the obvious reading of "set
    to a constant" is that they do not.

    Asserted as a PAIR, so the suite records that the roadmap's phrasing admits a form the test cannot
    fail: with `c1 = c0 = 0.45` the identity holds AND the dropped-correction implementation passes
    it too, at exactly +0.000000.
    """
    y, a, w, e = golden_arrays()
    h = outcome._tilt(e)
    m1 = m0 = np.full(len(y), 0.45)
    rd = outcome.weighted_rd(y, a, w)
    assert outcome.augmented_rd(y, a, w, h, m1, m0) == pytest.approx(rd, abs=1e-15)

    t, c = a == 1.0, a == 0.0
    dropped_term = (float(np.sum(w[t] * (y[t] - m1[t])) / np.sum(w[t]))
                    - float(np.sum(w[c] * (y[c] - m0[c])) / np.sum(w[c])))
    assert dropped_term - rd == pytest.approx(0.0, abs=1e-15)     # the WRONG one passes too


def test_item_3_THE_TILT_is_pinned_on_a_NON_CONSTANT_outcome_model():
    """Acceptance test 2, and it is the test that pins `h` AS AN ARGUMENT.

    Neither the constant-model tests nor any tolerance on them can distinguish `h` from `w` (§8.4),
    so the fixture has to be one on which they measurably differ — and it had to be MEASURED, not
    argued. The gap is `Σu·d/Σu − Σv·d/Σv` with `d = m1 − m0`, which is zero exactly when `d` is
    constant, so the power comes from the variation in `d`. But it is NOT MONOTONE in that variation:
    §8.4's six-candidate table has the construction with the LARGEST sd(d) producing the
    SECOND-SMALLEST gap, because a steeper propensity model can align the two normalised tilts. A
    fixture chosen on the argument "make the effect heterogeneous and the tilts must separate" lands
    at 2.5e-04 and produces a test that passes on the wrong tilting function.

    The bound is `1e-03`: 10x below the measured gap, and clear of any floating-point floor, because
    the constant-`d` sanity rows below cancel EXACTLY.
    """
    x, e, a, m1, m0, y = tilt_frame()
    w = np.where(a == 1.0, 1.0 - e, e)
    h = outcome._tilt(e)
    assert float(np.sum(w)) == pytest.approx(25.913563, abs=1e-6)
    assert float(np.sum(h)) == pytest.approx(13.083855, abs=1e-6)

    tau_h = outcome.augmented_rd(y, a, w, h, m1, m0)
    tau_w = outcome.augmented_rd(y, a, w, w, m1, m0)          # the route NOT taken, computed
    assert tau_h == pytest.approx(0.5112006369, abs=1e-9)
    assert tau_w == pytest.approx(0.5008939927, abs=1e-9)
    assert abs(tau_h - tau_w) == pytest.approx(1.030664e-02, rel=1e-4)
    assert abs(tau_h - tau_w) > 1e-03


@pytest.mark.parametrize("d", [0.00, 0.05, 0.40])
def test_item_3_THE_SANITY_ROWS_the_algebra_predicts_cancel_EXACTLY(d):
    """With `d = m1 − m0` constant the two tilts cancel to `0.000e+00` — not to a small number.

    That is what establishes there is no floating-point floor for the `1e-03` bound above to clear,
    which is the half of a tolerance argument that usually goes unstated.
    """
    x, e, a, _, _, y = tilt_frame()
    w = np.where(a == 1.0, 1.0 - e, e)
    m0 = np.full(len(y), 0.30)
    m1 = m0 + d
    gap = (outcome.augmented_rd(y, a, w, outcome._tilt(e), m1, m0)
           - outcome.augmented_rd(y, a, w, w, m1, m0))
    assert gap == 0.0


def test_item_4_the_correction_terms_denominator_is_the_WHOLE_population_and_not_an_arms():
    """RECORDED AS REDUNDANT AND KEPT ANYWAY. `h` is not arm-specific, so its denominator is the whole
    [§11] population.

    The arm-normalised implementation already fails item 1 decisively, so this asserts nothing item 1
    does not. It stays because it LOCALISES the failure: item 1 says "the augmentation is wrong",
    item 4 says "the correction term's denominator is wrong". `augmented_rd`'s docstring no longer
    claims item 3 is what catches it.

    **The size of that failure is a property of the frame and §15.8's +0.482 is not this frame's.**
    That number was measured on a 92-record ATO-shaped construction; on `golden_frame()` the
    arm-normalised form comes back +0.0813 from `rd`. Both are decisive — item 1's tolerance is
    `1e-15`, so this is thirteen orders of magnitude past it — and the assertion below is written
    against item 1's tolerance rather than against either measured gap, because that is the quantity
    that actually decides whether item 1 catches it.
    """
    y, a, w, e = golden_arrays()
    h = outcome._tilt(e)
    m1, m0 = np.full(len(y), 0.62), np.full(len(y), 0.23)
    t, c = a == 1.0, a == 0.0
    arm_normalised = (float(np.sum(w[t] * (y[t] - m1[t])) / np.sum(w[t]))
                      - float(np.sum(w[c] * (y[c] - m0[c])) / np.sum(w[c]))
                      + float(np.sum(h * (m1 - m0)) / np.sum(w[t])))
    correct = outcome.augmented_rd(y, a, w, h, m1, m0)
    assert arm_normalised != pytest.approx(correct, abs=1e-6)
    # and it already fails item 1, which is what makes this redundant-but-localising
    gap = abs(arm_normalised - outcome.weighted_rd(y, a, w))
    assert gap == pytest.approx(0.0812746102, abs=1e-9)
    assert gap > 1e-15 * 1e12


def test_item_5_h_IS_PASSED_IN_and_is_never_recomputed_inside_augmented_rd():
    """Given an `h` that is not `e(1−e)`, the result CHANGES — which is what makes it the caller's
    value rather than one the function chose."""
    y, a, w, e = golden_arrays()
    m1 = 1.0 / (1.0 + np.exp(-(0.4 * np.arange(len(y)) / len(y))))
    m0 = m1 - 0.2
    correct = outcome.augmented_rd(y, a, w, outcome._tilt(e), m1, m0)
    invented = outcome.augmented_rd(y, a, w, np.full(len(y), 1.0), m1, m0)
    assert correct != invented


@pytest.mark.parametrize("term", ["treated", "control", "tilt"])
def test_item_5b_the_ZERO_TOTAL_GUARD_names_which_of_the_three_denominators_was_empty(term):
    """`augmented_rd` is public and separately callable (§3.2), so S7 is not between it and a caller.
    A raw divide returns nan with a RuntimeWarning rather than saying which total was undefined."""
    y, a, w, e = golden_arrays()
    h = outcome._tilt(e)
    m1, m0 = np.full(len(y), 0.6), np.full(len(y), 0.3)
    if term == "tilt":
        h = np.zeros_like(h)
    else:
        w = np.where(a == (1.0 if term == "treated" else 0.0), 0.0, w)
    with pytest.raises(config.SchemaError) as e_info:
        outcome.augmented_rd(y, a, w, h, m1, m0)
    assert term in message_of(e_info)


def test_item_6_THE_CALL_SITE_IS_PINNED_which_is_what_items_1_to_5_DO_NOT_DO():
    """**The failure the spec exists to make catchable** (§8.2a).

    Items 1-5 all call `augmented_rd` directly with arrays the test constructed, so they prove the
    FUNCTION distinguishes `h` from `w` and say nothing about what `secondary` passes it. Without
    this assertion an implementation typing `w` on §11's one tilt line ships the wrong estimand with
    the whole rest of the suite green.

    Asserted as an EQUALITY against the `h` route and an inequality against the `w` route, on a frame
    where the two measurably differ — which the golden frame does, at a gap of 7.4e-04 on the reduced
    covariate list.

    **`golden_frame()` carries only the two covariates TICI's override names**, so the other six
    outcomes cannot be fitted on it at all — their list is the shared [§6] nine. They are supplied as
    copies of TICI's response and frozen to "unaugmented" through the `paths` parameter, which is
    exactly what that parameter is for and makes this test exercise it as well.
    """
    df = as_outcome(golden_frame())
    for key in config.BINARY_OUTCOMES:
        if key not in df.columns:
            df[key] = df["tici_2b_3"]
    frozen = {k: ("reduced" if k == "tici_2b_3" else "unaugmented")
              for k in config.BINARY_OUTCOMES}
    sec = outcome.secondary(df, golden_propensity(df),
                            data.Audit(hand_source(Path("."), "stage9")), paths=frozen)

    y, a, w, h, m1, m0 = golden_secondary_arrays()
    assert sec.estimates["tici_2b_3"].augmented == outcome.augmented_rd(y, a, w, h, m1, m0)
    assert sec.estimates["tici_2b_3"].augmented != outcome.augmented_rd(y, a, w, w, m1, m0)


def test_item_7_tilt_is_e_times_one_minus_e_against_a_hand_computation():
    e = np.array([0.0, 0.25, 0.5, 0.75, 1.0])
    assert list(outcome._tilt(e)) == [0.0, 0.1875, 0.25, 0.1875, 0.0]


def test_item_7_secondary_CALLS_tilt_and_the_expression_cannot_migrate_back_INLINE(stage9_source):
    """By scan, because the behavioural test above cannot see a second definition that agrees.

    `h = e(1 − e)` is the [§7] tilting function that DEFINES THE ESTIMAND, and it lived as a bare
    sub-expression at one call site with nothing naming it. One home, one citation, one place to be
    wrong (§8.2a).
    """
    body = stage9_source.split("def secondary(", 1)[1]
    assert "_tilt(e)" in body
    assert "e * (1.0 - e)" not in body and "e * (1 - e)" not in body
    # and exactly one definition of it in the whole half
    assert stage9_source.count("def _tilt(") == 1


def test_the_two_TILTING_FUNCTIONS_are_not_proportional_which_is_what_makes_item_3_possible():
    """§8.2. `w` is arm-specific — `1 − e` treated and `e` control — and `h` is not. Measured on the
    workbook, Σw = 27.736623 against Σh = 14.797901; on the golden frame, 8.89 against 5.4257. If the
    two were proportional every tilt test in this section would be vacuous."""
    y, a, w, e = golden_arrays()
    h = outcome._tilt(e)
    assert float(np.sum(w)) == pytest.approx(8.89, abs=1e-12)
    assert float(np.sum(h)) == pytest.approx(5.4257, abs=1e-12)
    ratios = h / w
    assert float(np.ptp(ratios)) > 0.1            # not a constant multiple of each other


# --- 15.8 (the regression pin)  §15.0.2's golden vector, both covariate paths ---------------------------
#
# Pinned as literals from `golden_frame()` — 24 records, every value written out, no RNG — so that NO
# PATIENT-DERIVED NUMBER AND NO ESTIMATE ENTERS GIT (§4.3). The pins are asserted to 1e-9 rather than
# to machine precision, for Stage 6 §12.0.2's reason: this is a regression pin against an edit to our
# own arithmetic, and pinning it at machine precision would make it fail on a numpy patch release
# rather than on a mistake.
#
# **THE TWO PINS ARE A PAIR AND THE SECOND ONE IS THE POINT.** Adding ONE continuous covariate to a
# 24-record frame moves tau from +0.2233 — within 0.062 of the unaugmented +0.2851 — to +0.0403, a
# SEVEN-FOLD COLLAPSE toward zero. Six parameters against a minority cell of 11 over-fits: the fitted
# values approach Y, the correction term stops carrying signal, `max_abs_beta` goes from 1.588 to
# 11.828, the iteration count from 7 to 10, and THE TRUST REGION FIRES ONCE. That is [§8]'s "material
# disagreement between the two is evidence about the outcome model, not confirmation of either" as an
# executable demonstration, and it is why §10.3 prints RD_w and tau in adjacent columns. It is also a
# live argument for the rare-outcome rule: the failure mode that rule guards against is visible at a
# minority cell of 11, which is ABOVE the threshold of 10.

GOLDEN_MARGINAL = {"Sw": 8.8900000000, "Sh": 5.4257000000,
                   "p1_w": 0.6191536748, "p0_w": 0.3340909091,
                   "RD_w": 0.2850627657, "OR_w": 3.2404025938}
GOLDEN_AUGMENTED = {
    ("center", "atrial_fib"): dict(
        columns=4, dropped=("center_USZ",), iterations=7, converged_on="likelihood",
        max_abs_beta=1.5877919570, rescales=0, halvings=0, first_step_norm=2.4935249391,
        tau=0.2233300473, tau_under_w=0.2240676700),
    ("center", "atrial_fib", "age"): dict(
        columns=5, dropped=("center_USZ",), iterations=10, converged_on="likelihood",
        max_abs_beta=11.8282581834, rescales=1, halvings=0, first_step_norm=10.7967589095,
        tau=0.0402879733, tau_under_w=0.0402999125),
}


def test_the_golden_vector_reproduces_the_MARGINAL_quantities():
    y, a, w, e = golden_arrays()
    odds_ratio, corrected, share = outcome.marginal_odds_ratio(y, a, w)
    assert float(np.sum(w)) == pytest.approx(GOLDEN_MARGINAL["Sw"], abs=1e-9)
    assert float(np.sum(outcome._tilt(e))) == pytest.approx(GOLDEN_MARGINAL["Sh"], abs=1e-9)
    assert share[1] == pytest.approx(GOLDEN_MARGINAL["p1_w"], abs=1e-9)
    assert share[0] == pytest.approx(GOLDEN_MARGINAL["p0_w"], abs=1e-9)
    assert outcome.weighted_rd(y, a, w) == pytest.approx(GOLDEN_MARGINAL["RD_w"], abs=1e-9)
    assert odds_ratio == pytest.approx(GOLDEN_MARGINAL["OR_w"], abs=1e-9)
    assert corrected is False


def test_the_golden_frames_minority_cell_is_ELEVEN_which_is_ABOVE_the_threshold():
    """Which is what makes the pair below an argument about the rare-outcome rule rather than an
    illustration of it: this frame QUALIFIES for the full covariate list, and the full list is what
    destroys the estimate."""
    y, _, _, _ = golden_arrays()
    assert (int((y == 1.0).sum()), int((y == 0.0).sum())) == (11, 13)
    assert min(11, 13) > config.RARE_MINORITY_THRESHOLD


@pytest.mark.parametrize("covariates", sorted(GOLDEN_AUGMENTED))
def test_the_golden_vector_reproduces_each_augmented_path(covariates):
    pins = GOLDEN_AUGMENTED[covariates]
    df = as_outcome(golden_frame())
    fit, X, dropped = outcome.outcome_model(df, pd.Series(True, index=df.index), "tici_2b_3") \
        if covariates == ("center", "atrial_fib") else _golden_fit_over(df, covariates)
    y, a, w, e = golden_arrays()
    m1, m0 = outcome._counterfactuals(fit, X)

    assert (len(fit.columns), dropped) == (pins["columns"], pins["dropped"])
    assert (fit.iterations, fit.converged_on) == (pins["iterations"], pins["converged_on"])
    assert (fit.rescales, fit.halvings) == (pins["rescales"], pins["halvings"])
    assert float(np.max(np.abs(fit.beta))) == pytest.approx(pins["max_abs_beta"], abs=1e-9)
    assert fit.first_step_norm == pytest.approx(pins["first_step_norm"], abs=1e-9)
    assert np.array_equal(model.predict(fit, X), fit.p)
    assert outcome.augmented_rd(y, a, w, outcome._tilt(e), m1, m0) == pytest.approx(
        pins["tau"], abs=1e-9)
    assert outcome.augmented_rd(y, a, w, w, m1, m0) == pytest.approx(pins["tau_under_w"], abs=1e-9)


def _golden_fit_over(df, covariates):
    """`outcome_model`'s body with an explicit covariate list, for the pin that adds `age`.

    `outcome_model` takes an OUTCOME and looks its list up in the registry, which is exactly the
    property §7.2 wants — an outcome can only diverge by being named — so there is no way to ask it
    for a third list, and the second pin has to assemble the design itself. That it cannot be
    requested through the public function is the design working, not a gap in it.
    """
    X, dropped = model.design(df, covariates)
    X = X.copy()
    X.insert(0, config.TREATMENT, df[config.TREATMENT].to_numpy(dtype=float))
    return model.firth(X, df["tici_2b_3"].to_numpy(dtype=float)), X, dropped


def test_ONE_EXTRA_COVARIATE_COLLAPSES_TAU_SEVENFOLD_and_that_is_the_pair_being_asserted():
    """The two pins read together, asserted as the relationship they exist to show.

    Not a sign flip — an earlier draft of the spec claimed one, against a frame it did not contain —
    but a collapse toward zero while the trust region fires and `max_abs_beta` moves by an order of
    magnitude. §10.3 prints RD_w beside tau because of exactly this.
    """
    small = GOLDEN_AUGMENTED[("center", "atrial_fib")]
    large = GOLDEN_AUGMENTED[("center", "atrial_fib", "age")]
    rd = GOLDEN_MARGINAL["RD_w"]
    assert abs(rd - small["tau"]) < 0.062                      # the four-column fit tracks rd
    assert small["tau"] / large["tau"] > 5.0                   # the five-column one collapses
    assert large["max_abs_beta"] / small["max_abs_beta"] > 7.0
    assert (small["rescales"], large["rescales"]) == (0, 1)    # and one of each, for the suite


# --- 15.9  the rare-minority guard, and it is decided ONCE ---------------------------------------------

def test_augmentable_reads_the_THRESHOLD_and_never_a_literal_ten():
    """`True` at the threshold and `False` one below it, asserted against
    `C.RARE_MINORITY_THRESHOLD` and never against `10` — so moving the constant moves the test."""
    threshold = config.RARE_MINORITY_THRESHOLD
    assert outcome._augmentable("mrs_0_2_90d", threshold) is True
    assert outcome._augmentable("mrs_0_2_90d", threshold - 1) is False
    assert outcome._augmentable("mrs_0_2_90d", threshold + 1) is True


def test_AN_OVERRIDE_IS_A_ROUTE_IN_and_not_a_modifier_of_an_outcome_that_qualified_anyway():
    """The roadmap-versus-SAP contradiction of §7.2, asserted in the direction the amendment settles.

    The roadmap said outcomes below the threshold "are not augmented", with TICI named as the case
    the rule is written for. The SAP amendment's item 2 says the opposite in terms — such an outcome
    is augmented with a DECLARED REDUCED m_a(X) "rather than dropped from augmentation" — and Stage
    0's DECISION 3, in the roadmap's own file, states it correctly. One stale sentence against three
    agreeing sources.
    """
    below = config.RARE_MINORITY_THRESHOLD - 4
    assert outcome._augmentable("tici_2b_3", below) is True        # the override carries it
    assert outcome._augmentable("sich", below) is False            # the same cell, no override
    assert "tici_2b_3" in config.OUTCOME_MODEL_OVERRIDES and "sich" not in config.OUTCOME_MODEL_OVERRIDES


def test_the_override_is_UNCONDITIONAL_where_the_SAPs_reduction_is_CONDITIONAL():
    """§16 item 8, asserted so that the known gap is visible in the suite rather than only in a list.

    The amendment reads "WHERE the minority cell CANNOT SUPPORT the full covariate set but can
    support a smaller one..."; `_augmentable` returns True for a key in OUTCOME_MODEL_OVERRIDES
    WHATEVER its minority cell is. §9.5 presents that as a design virtue — naming an outcome fixes
    its estimator, which is what makes `tici_2b_3` immune to the resampling variation that moves
    `ph2` across the threshold in 46.1% of replicates — and it is one. It also means a workbook
    revision taking TICI's minority cell to 40 would still fit the two-name reduced model.
    """
    assert outcome._augmentable("tici_2b_3", 40) is True
    assert outcome._augmented_path("tici_2b_3", 40) == "reduced"
    assert config.outcome_model_covariates("tici_2b_3") == ("center", "atrial_fib")


@pytest.mark.parametrize("outcome_key,minority,expected", [
    ("mrs_0_2_90d", 42, "full"), ("mrs_0_2_90d", 9, "unaugmented"),
    ("tici_2b_3", 6, "reduced"), ("sich", 5, "unaugmented"), ("ph2", 10, "full")])
def test_augmented_path_returns_the_THREE_DECLARED_STRINGS(outcome_key, minority, expected):
    assert outcome._augmented_path(outcome_key, minority) == expected


def test_THE_BLOCK_OF_FOUR_moves_together_across_all_seven_outcomes():
    """§3.1. `augmented`, `covariates`, `fit` and `dropped` are None/None/None/() exactly when the
    path is "unaugmented", and all four populated otherwise.

    A reader must not have to infer the path from which fields are absent, so `augmented_path` states
    it — and this is the assertion that the two never disagree.
    """
    _, _, sec, _ = secondary_run()
    for est in sec.estimates.values():
        block = (est.augmented, est.covariates, est.fit, est.dropped)
        if est.augmented_path == "unaugmented":
            assert block == (None, None, None, ())
        else:
            assert est.augmented is not None and np.isfinite(est.augmented)
            assert est.covariates and est.fit is not None
            assert est.dropped == ("center_USZ",)
    assert {e.augmented_path for e in sec.estimates.values()} == {
        "full", "reduced", "unaugmented"}                          # all three paths reached


def test_reduced_IS_ASSERTED_SEPARATELY_and_IS_NOT_IN_THE_BLOCK():
    """§3.1's correction, and it is false on four of the seven.

    `reduced` is True exactly when the key is in OUTCOME_MODEL_OVERRIDES, which is ORTHOGONAL to the
    path: it is False on the `full`-path outcomes whose other three fields ARE populated. A
    block-move test written from the earlier draft's docstring — which listed `reduced` in the block
    — fails on every one of them, and the companion below is that failure made explicit.
    """
    _, _, sec, _ = secondary_run()
    for key, est in sec.estimates.items():
        assert est.reduced is (key in config.OUTCOME_MODEL_OVERRIDES)

    would_have_failed = [k for k, e in sec.estimates.items()
                         if e.augmented_path == "full" and e.reduced is False]
    assert len(would_have_failed) >= 4, (
        "the earlier draft's docstring put `reduced` in the block; these are the outcomes a test "
        f"written from it fails on: {would_have_failed}")


def test_augmented_path_reduced_AND_reduced_True_AGREE_ON_ALL_SEVEN():
    """The invariant that makes the two fields non-contradictory even though only one is in the
    block: `augmented_path == "reduced"` implies `reduced`, and a `reduced` outcome is never
    `unaugmented`, because an override is a route IN."""
    _, _, sec, _ = secondary_run()
    for est in sec.estimates.values():
        if est.augmented_path == "reduced":
            assert est.reduced is True
        if est.reduced:
            assert est.augmented_path == "reduced"


def test_the_minority_cell_is_counted_ON_THE_MASKED_POPULATION_and_not_on_in_model():
    """§9.1, on the frame where the two differ — which is the same construction §15.2 uses, reused
    so the two facts are checked on one frame.

    The rule is applied on the [§11] population because that is THE POPULATION THE NUISANCE MODEL IS
    FITTED ON, and the constraint it encodes is that a model cannot carry more parameters than its
    own fitting sample supports.
    """
    df, ps, sec, _ = secondary_run()
    est = sec.estimates["tici_2b_3"]
    masked = df.loc[est.in_estimate, "tici_2b_3"].to_numpy(dtype=float)
    over_in_model = df.loc[ps.in_model, "tici_2b_3"].to_numpy(dtype=float)

    assert est.minority == min(int((masked == 1.0).sum()), int((masked == 0.0).sum())) == 4
    assert min(int((over_in_model == 1.0).sum()),
               int((~(over_in_model == 1.0)).sum())) == 5          # the count from in_model differs


def test_outcome_model_covariates_defaults_for_six_and_overrides_for_one():
    """The first READ of a function written one stage early (§13), so it is a tested read."""
    for key in config.BINARY_OUTCOMES:
        if key in config.OUTCOME_MODEL_OVERRIDES:
            assert config.outcome_model_covariates(key) == ("center", "atrial_fib")
        else:
            assert config.outcome_model_covariates(key) is config.OUTCOME_COVARIATES
    assert config.TREATMENT not in config.outcome_model_covariates("tici_2b_3")
    with pytest.raises(KeyError):
        config.outcome_model_covariates("mrs_0_3_90d")


def test_TICIS_REDUCED_DESIGN_IS_FIVE_PARAMETERS_and_the_contingency_is_asserted_WITH_it():
    """§7.2's fragility fails a test rather than living in a paragraph.

    The amendment says "five parameters, chosen on procedural grounds". Measured, that is right — and
    it is right BY ACCIDENT. `center` has FOUR declared levels, which is three dummies; `design` then
    drops `center_USZ` as constant because USZ contributes no records. So the fitted design is
    treatment + 2 centre dummies + atrial_fib = 4 columns, plus the intercept `firth` prepends = five.
    If USZ ever contributes a single record the same declared model becomes SIX parameters against six
    non-events — saturated — and the amendment's stated arithmetic silently stops holding.
    """
    df = as_outcome(golden_frame())
    fit, X, dropped = outcome.outcome_model(df, pd.Series(True, index=df.index), "tici_2b_3")
    assert len(config.FACTOR_LEVELS["center"]) == 4
    assert dropped == ("center_USZ",)
    assert len(fit.columns) == 4
    assert len(fit.beta) == 5                                      # + the prepended intercept

    with_usz = df.copy()
    with_usz.loc[with_usz.index[-1], "center"] = "USZ"
    usz_fit, _, usz_dropped = outcome.outcome_model(
        with_usz, pd.Series(True, index=with_usz.index), "tici_2b_3")
    assert usz_dropped == ()                                       # nothing constant any more
    assert len(usz_fit.beta) == 6                                  # SIX parameters, saturated


# --- 15.9 (the frozen path)  §11.1's parameter, and §9.5's 46.1% as an assertion ------------------------

def test_A_FROZEN_PATH_OVERRIDES_A_CROSSED_THRESHOLD():
    """§9.5's 46.1% as an assertion rather than as a measurement.

    A minority cell is a property of the SAMPLE and [§10] resamples, so a replicate can carry `ph2`
    over the threshold. If the path were re-decided per replicate, some replicates would contribute
    an augmented estimate and others an unaugmented one to the same percentile interval — a MIXTURE
    OF TWO ESTIMATORS, which [§10] forbids in terms: "Replicates whose prespecified fit fails are
    dropped and counted; never substituted with a different estimator."

    On the fixture, `ph2`'s minority cell is 8 and it is unaugmented. Raise it over the threshold and
    `paths=None` augments it; the point estimate's frozen map leaves it alone. That is the mechanism.
    """
    df = secondary_cohort()
    df["ph2"] = [float(i % 2 == 0 or i % 5 == 2) for i in range(len(df))]     # crosses to 12/28
    _, _, crossed, _ = secondary_run(df)
    assert crossed.estimates["ph2"].minority > config.RARE_MINORITY_THRESHOLD
    assert crossed.estimates["ph2"].augmented_path == "full"                  # re-decided: augmented

    _, _, point, _ = secondary_run()
    frozen = {k: e.augmented_path for k, e in point.estimates.items()}
    assert frozen["ph2"] == "unaugmented"
    _, _, held, _ = secondary_run(df, paths=frozen)
    assert held.estimates["ph2"].augmented_path == "unaugmented"
    assert held.estimates["ph2"].augmented is None                            # and no tau at all


def test_the_frozen_map_carries_EVERY_path_unchanged_and_the_estimates_agree_where_it_agrees():
    """Stage 10 reads `augmented_path` off the point estimate's seven `BinaryEstimate`s and passes
    the map straight back. It never re-derives the rule (§11.1), and on the point estimate's own
    frame doing so changes nothing — which is what makes the parameter safe to always pass."""
    df, ps, point, _ = secondary_run()
    frozen = {k: e.augmented_path for k, e in point.estimates.items()}
    _, _, again, _ = secondary_run(df, ps, paths=frozen)
    for key in config.BINARY_OUTCOMES:
        assert again.estimates[key].augmented_path == point.estimates[key].augmented_path
        assert again.estimates[key].augmented == point.estimates[key].augmented
        assert again.estimates[key].rd == point.estimates[key].rd


@pytest.mark.parametrize("partial", ["one_missing", "one_extra", "empty"])
def test_A_PARTIAL_paths_MAP_RAISES_because_the_mixture_would_arrive_ONE_OUTCOME_AT_A_TIME(partial):
    """Complete or absent, never partial (§11.1). A permissive `paths.get(key, _augmented_path(...))`
    would introduce the failure mode while looking like a convenience."""
    df, ps, point, _ = secondary_run()
    frozen = {k: e.augmented_path for k, e in point.estimates.items()}
    broken = {"one_missing": {k: v for k, v in frozen.items() if k != "ph2"},
              "one_extra": {**frozen, "mrs_0_3_90d": "full"},
              "empty": {}}[partial]
    with pytest.raises(config.SchemaError) as e:
        outcome.secondary(df, ps, data.Audit(hand_source(Path("."), "stage9")), paths=broken)
    assert "frozen paths" in message_of(e)


def test_paths_None_is_the_POINT_ESTIMATE_CALL_and_must_be_able_to_run_without_them():
    """Optional, defaulting to None: the point-estimate call is the one that DECIDES the paths, so a
    required parameter would make the first call impossible."""
    _, _, sec, _ = secondary_run()
    assert all(e.augmented_path == outcome._augmented_path(k, e.minority)
               for k, e in sec.estimates.items())


@DATA_GATED
def test_the_workbook_paths_are_FOUR_FULL_ONE_REDUCED_TWO_UNAUGMENTED(workbook_stage9):
    """`[data-gated]`, and this is the SAP amendment's own closing prediction confirmed rather than
    assumed: "No other outcome is reduced; symptomatic ICH and parenchymal haematoma type 2 have a
    minority cell below the threshold and remain unaugmented."

    It is the only prespecified numerical claim the protocol makes about this workbook's outcome
    distribution. **And two of the four SAFETY outcomes are augmented with the full [§6] list**,
    which is recorded because a reader who takes "the safety family is unaugmented" from [§8]'s
    "will normally apply to the safety family" will find two augmented safety estimates in the output.
    """
    df, ps, sec, audit = workbook_stage9
    paths = {k: e.augmented_path for k, e in sec.estimates.items()}
    assert sorted(paths.values()) == ["full"] * 4 + ["reduced"] + ["unaugmented"] * 2
    assert paths["tici_2b_3"] == "reduced"
    assert {k for k, v in paths.items() if v == "unaugmented"} == {"sich", "ph2"}
    augmented_safety = {k for k, e in sec.estimates.items()
                        if e.family == "safety" and e.augmented is not None}
    assert augmented_safety == {"death_90d", "mrs_5_6_90d"}


@DATA_GATED
def test_ph2_SITS_ONE_BELOW_THE_THRESHOLD_on_the_workbook_and_that_is_disclosed(workbook_stage9):
    """`[data-gated]`, a COUNT. §9.4: one additional parenchymal haematoma in the ATO population
    flips `ph2` from unaugmented to fully augmented — a DIFFERENT ESTIMATOR, not a different number
    from the same one.

    The threshold is prespecified at Stage 1 and is not revisited; §17 records that it was looked at
    and left alone. What this stage owes instead is disclosure, and §10.2's table prints the minority
    cell beside the path for every outcome so a reader can see it without reading the spec.
    """
    df, ps, sec, audit = workbook_stage9
    assert sec.estimates["ph2"].minority == config.RARE_MINORITY_THRESHOLD - 1
    assert outcome._augmented_path("ph2", sec.estimates["ph2"].minority + 1) == "full"
    table = rows_of(next(e for e in audit.entries if e.step == "binary_estimable"))
    assert table["ph2"][3] == str(sec.estimates["ph2"].minority)


# --- 15.10  NOT doubly robust, and the assertion is POSITIVE `[roadmap, amended]` ------------------------
#
# The roadmap required: "a heterogeneous-effect scenario with a correct outcome model and a misspecified
# propensity model does not recover the true ATO. This must be asserted, not merely allowed." §8.3's
# argument predicts something SHARPER than non-recovery, and the sharper version is what is asserted
# here: the estimator should recover the ATO indexed by the MISSPECIFIED score, accurately, while
# missing the one indexed by the true score.
#
# Why the two errors are so close in size, since that is not luck: with an ORACLE `m_a` the two arm
# terms of §8.1 have mean zero and the correction term is `Σh(m1−m0)/Σh` evaluated at the misspecified
# `h` — which IS `ATO(e_mis)` exactly. So the estimator is not approximately recovering the wrong
# population, it is recovering it BY CONSTRUCTION, and the residual is sampling error in the arm terms.
#
# **AND THE TEST CANNOT BE RUN AT THE COHORT'S OWN SIZE.** At n = 92 the per-replicate spread swamps the
# quantity being asserted: measured over 400 replicates, mean error about −0.026 with sd about 0.088 —
# the sd is more than three times the bias. A test written on a 92-row fixture measures noise. The `n`
# is chosen from a MEASURED SEPARATION and is this section's only tuning constant.

DR_N: Final[int] = 200_000


def test_the_HETEROGENEOUS_arm_MISSES_the_true_ATO_by_more_than_the_asserted_bound():
    """The negative half. Oracle `m_a(X)`, misspecified `e(X)`, and the estimand is `ATO(e_true)`."""
    e_true, e_mis, a, m1, m0, y = dr_population(DR_N, DR_SEED, "heterogeneous")
    w_mis = np.where(a == 1.0, 1.0 - e_mis, e_mis)
    tau_hat = outcome.augmented_rd(y, a, w_mis, outcome._tilt(e_mis), m1, m0)

    assert ato(e_true, m1, m0) == pytest.approx(+0.320705, abs=1e-6)
    assert tau_hat == pytest.approx(+0.290456, abs=1e-6)
    assert abs(tau_hat - ato(e_true, m1, m0)) > 0.015


def test_AND_IT_RECOVERS_WHAT_8_3_PREDICTS_IT_RECOVERS_which_is_the_POSITIVE_half():
    """This is what distinguishes "not doubly robust" from "wrong".

    The estimate sits within 0.010 of the ATO indexed by the MISSPECIFIED score — measured 0.005428 —
    against an error of 0.030249 against the estimand. The estimator stays consistent for a weighted
    average treatment effect; it is consistent for the WRONG WEIGHTING, which is [§8]'s paragraph as
    a number.
    """
    e_true, e_mis, a, m1, m0, y = dr_population(DR_N, DR_SEED, "heterogeneous")
    w_mis = np.where(a == 1.0, 1.0 - e_mis, e_mis)
    tau_hat = outcome.augmented_rd(y, a, w_mis, outcome._tilt(e_mis), m1, m0)

    assert ato(e_mis, m1, m0) == pytest.approx(+0.295884, abs=1e-6)
    assert abs(tau_hat - ato(e_mis, m1, m0)) == pytest.approx(0.005428, abs=1e-6)
    assert abs(tau_hat - ato(e_mis, m1, m0)) < 0.010


def test_with_a_CORRECTLY_SPECIFIED_score_the_same_estimator_recovers_the_estimand():
    """The contrast that makes the two above about the PROPENSITY model rather than about the
    estimator: −0.004149 with `e` correct, against −0.030249 with `e` misspecified."""
    e_true, _, a, m1, m0, y = dr_population(DR_N, DR_SEED, "heterogeneous")
    w_true = np.where(a == 1.0, 1.0 - e_true, e_true)
    tau_hat = outcome.augmented_rd(y, a, w_true, outcome._tilt(e_true), m1, m0)
    assert tau_hat - ato(e_true, m1, m0) == pytest.approx(-0.004149, abs=1e-6)


def test_the_CONSTANT_RISK_DIFFERENCE_companion_recovers_the_estimand():
    """The companion, and its construction is the point.

    The roadmap's reason for insisting the assertion be positive is that "a test built on a
    homogeneous effect will appear to show that it is [the thing §8.3 forbids calling it], because
    every weighted average treatment effect coincides in that case". That is true when `tau(X)` is CONSTANT — and a logistic
    outcome model with no interaction term is NOT that, because the link's curvature makes a constant
    LOG-ODDS effect a non-constant RISK DIFFERENCE. So the companion is built as `m0` from a logistic
    model and `m1 = m0 + DELTA`, with `m0` bounded into (0.25, 0.70) so that `m1` stays a probability:
    without the bound the difference truncates at the top and the construction is no longer constant
    where it matters most.
    """
    e_true, e_mis, a, m1, m0, y = dr_population(DR_N, DR_SEED, "constant")
    # `m1` is built as `m0 + DELTA`, so `m1 - m0` recovers DELTA to floating-point subtraction error
    # and not exactly: measured, sd 8.5e-17 across 200,000 records. Asserted at that scale rather
    # than at 0.0, which would be a test of IEEE 754 rather than of the construction.
    assert float(np.std(m1 - m0)) < 1e-15
    assert m1.max() <= 1.0 and m0.min() >= 0.25            # the bound, and it is load-bearing
    assert (m1 - m0) == pytest.approx(np.full(len(m1), DELTA), abs=1e-15)

    w_mis = np.where(a == 1.0, 1.0 - e_mis, e_mis)
    tau_hat = outcome.augmented_rd(y, a, w_mis, outcome._tilt(e_mis), m1, m0)
    assert abs(tau_hat - ato(e_true, m1, m0)) < 0.010


def test_THE_BAND_BETWEEN_THE_TWO_BOUNDS_IS_EMPTY_ACROSS_ALL_TWELVE_SEEDS():
    """Stage 8 §6.3's empty-band calibration applied to a TOLERANCE instead of to a coefficient bound.

    The two asserted bounds — companion below 0.010, heterogeneous above 0.015 — are not round
    numbers chosen for looking reasonable. Across `DR_SEEDS_12` the heterogeneous arm's errors span
    [0.022950, 0.030249] and the companion's span [0.000176, 0.004929]; nothing lies between 0.004929
    and 0.022950, and both bounds sit inside that gap. A tolerance loose enough to survive a numpy
    version bump and tight enough to fail on a wrong estimator has to sit there, and the window is
    MEASURED rather than assumed.
    """
    def errors(arm):
        out = []
        for seed in DR_SEEDS_12:
            e_true, e_mis, a, m1, m0, y = dr_population(DR_N, seed, arm)
            w_mis = np.where(a == 1.0, 1.0 - e_mis, e_mis)
            out.append(abs(outcome.augmented_rd(y, a, w_mis, outcome._tilt(e_mis), m1, m0)
                           - ato(e_true, m1, m0)))
        return out

    heterogeneous, constant = errors("heterogeneous"), errors("constant")
    assert (min(heterogeneous), max(heterogeneous)) == pytest.approx((0.022950, 0.030249), abs=1e-6)
    assert (min(constant), max(constant)) == pytest.approx((0.000176, 0.004929), abs=1e-6)
    assert max(constant) < 0.010 < 0.015 < min(heterogeneous)      # the band, and both bounds in it


def test_THE_NO_INTERACTION_COMPANION_IS_ITSELF_BIASED_which_is_why_it_is_not_the_companion():
    """§8.5's finding, asserted so the suite records that the OBVIOUS "homogeneous" construction is
    not homogeneous on the risk-difference scale.

    A no-interaction logistic outcome model is what a reasonable draft would have written as the
    companion. Measured over 200 replicates at n = 20,000: it is detectably biased at |t| about 22,
    while the constructed constant-risk-difference arm at the same size is indistinguishable from
    zero. The bias tracks `sd(tau(X))`, which is 0.24, 0.05 and 0.00 across the three arms.

    The thresholds are `|t| > 10` and `|t| < 3` rather than the measured values, because the seed set
    for this measurement is not pinned by the spec and a two-significant-figure `|t|` is not a
    reproducible target. The SEPARATION is what is being asserted, and it is two orders of magnitude.
    """
    def t_statistic(arm, replicates=200, n=20_000):
        errs = []
        for r in range(replicates):
            e_true, e_mis, a, m1, m0, y = dr_population(n, DR_SEED + 1000 + r, arm)
            w_mis = np.where(a == 1.0, 1.0 - e_mis, e_mis)
            errs.append(outcome.augmented_rd(y, a, w_mis, outcome._tilt(e_mis), m1, m0)
                        - ato(e_true, m1, m0))
        e = np.array(errs)
        return float(abs(e.mean() / (e.std(ddof=1) / np.sqrt(len(e))))), float(e.mean())

    no_interaction, mean_ni = t_statistic("no_interaction")
    constant, mean_c = t_statistic("constant")
    assert no_interaction > 10.0, f"the no-interaction arm should be biased; |t| = {no_interaction}"
    assert constant < 3.0, f"the constant-risk-difference arm should not be; |t| = {constant}"
    assert abs(mean_ni) > 10 * abs(mean_c)


def test_the_THREE_ARMS_differ_ONLY_in_the_outcome_model_and_not_in_the_population():
    """The three arms share x1, x2 and the treatment draw for a given seed, which is what makes the
    comparison above about `tau(X)` and not about three different populations."""
    draws = {arm: dr_population(2000, DR_SEED, arm)
             for arm in ("heterogeneous", "constant", "no_interaction")}
    first = draws["heterogeneous"]
    for arm, drawn in draws.items():
        assert np.array_equal(drawn[0], first[0])          # e_true
        assert np.array_equal(drawn[1], first[1])          # e_mis
        assert np.array_equal(drawn[2], first[2])          # a
    assert not np.array_equal(draws["constant"][3], first[3])          # but m1 differs


def test_AT_THE_COHORTS_OWN_SIZE_THE_SAME_MEASUREMENT_IS_NOISE_and_this_is_a_NOTE():
    """Not an assertion about bias — a demonstration that the assertion cannot be made here.

    Measured over 400 replicates at n = 92: mean error about −0.026 with sd about 0.088, so the
    spread is more than three times the quantity. §15.10 runs at n = 200,000 for exactly this reason,
    and the docstring is where the reason lives rather than in a comment beside a magic number.
    """
    errs = []
    for r in range(400):
        e_true, e_mis, a, m1, m0, y = dr_population(92, DR_SEED + 2000 + r, "heterogeneous")
        w_mis = np.where(a == 1.0, 1.0 - e_mis, e_mis)
        try:
            errs.append(outcome.augmented_rd(y, a, w_mis, outcome._tilt(e_mis), m1, m0)
                        - ato(e_true, m1, m0))
        except config.SchemaError:
            continue                          # a 92-row draw can leave an arm with no weight at all
    e = np.array(errs)
    assert len(e) > 350
    assert float(e.std(ddof=1)) > 3.0 * abs(float(e.mean()))
    assert float(e.std(ddof=1)) > 0.05


# --- 15.11  the module boundary -------------------------------------------------------------------------

def test_outcome_py_does_not_import_balance():
    """Stage 7 §11 says Stage 8 does not read the `Balance`; the same holds one stage on. By scan,
    for Stage 7 §14.11's reason: not because reading it would fail, but because it would SUCCEED.

    DECISION 4's requirement that the residual imbalance be named beside the primary estimate is
    Stage 14's obligation and not a data dependency here.
    """
    assert "balance" not in _imported_roots(SOURCE)


def test_secondary_DOES_NOT_READ_Primary(stage9_names):
    """Stage 8 §11: "the binary estimators are their own estimates on their own [§11] denominators,
    and the primary ordinal quantity is not an input to any of them".

    §15.11's reason for asserting it: not because reading it would fail, but because it would
    SUCCEED. `Primary` is in the same module and one attribute access away.
    """
    assert not {"Primary", "primary", "PRIMARY_OUTCOME", "cumulative_rd", "MRS_THRESHOLDS"} & stage9_names
    assert "Primary" in SOURCE                   # it is right there, which is the point


def test_model_polr_IS_NEITHER_CALLED_NOR_EDITED_by_this_half(stage9_names):
    """§0.2, recorded so that nobody wires the ordinal fitter in looking for a use. No binary estimate
    goes near it — and `model.firth` is one line away, so the wrong one is easy to reach for."""
    assert not {"polr", "model.polr", "PolrFit"} & stage9_names
    assert {"firth", "design", "predict"} <= stage9_names
    assert "model.polr" in SOURCE                # Stage 8's half does call it, so this is a contrast


def test_ps_e_and_ps_w_are_read_in_EXACTLY_TWO_PLACES_both_loc_of_the_outcome_mask(stage9_source):
    """§12. Stage 6 §9's rule — range over the mask, never over `notna()`, and never fill — honoured
    as a fact about two lines rather than as a warning.

    **THE SCAN IS SCOPED TO `.loc` READS AND THAT IS DELIBERATE.** `_assert_secondary_inputs`
    legitimately reads `ps.e[ps.in_model]` and `ps.e[~ps.in_model]`: S4's entire content is that `e`
    and `w` are finite ON the mask and nan OFF it, which cannot be checked without looking at both
    sides. What §12's claim is about is the ESTIMATION path — the reads whose result reaches an
    estimate — and every one of those is a `.loc[in_estimate]`.
    """
    reads = re.findall(r"ps\.(?:e|w)\.loc\[([^\]]+)\]", stage9_source)
    assert reads == ["in_estimate", "in_estimate"], reads
    assert "ps.e.loc[ps.in_model]" not in stage9_source
    assert "ps.w.loc[ps.in_model]" not in stage9_source
    assert ".notna()]" not in stage9_source.replace("df[key].notna()", "")


def test_NOTHING_FILLS_A_NAN_and_the_ONE_dropna_is_a_DOMAIN_CHECK(stage9_names):
    """Stage 6 §9's rule: the deliberateness of an absence lives in a mask, never in a value. `w` and
    `e` carry nan off `in_model` by construction and the tempting repair is a `fillna(0)`.

    **§15.11 lists `dropna` alongside `fillna` and that is too broad, which the code shows.** There
    is exactly one `dropna` in this half, in S5b, and it is not filling anything: S5b asks two
    separate questions of the treatment column — is any value MISSING, and is any value OUTSIDE the
    declared arm codes — and the second cannot be asked of a nan, because `isin` against one is
    False and the record would be reported as an out-of-range arm rather than as an absent one. The
    missing count is already reported by the branch above it. This is Stage 6's D3/D4 split exactly:
    "D3 excludes missing values from its own check because they are D4's, not because they are
    harmless." Asserted as a location rather than waved through.
    """
    assert not {"fillna", "interpolate", "ffill", "bfill", "bfill", "pad"} & stage9_names
    body = SOURCE.split("def _assert_secondary_inputs(", 1)[1].split("\ndef ", 1)[0]
    assert SOURCE.split(STAGE_9_BANNER, 1)[1].count("dropna") == body.count("dropna") == 2


def test_the_stage_adds_no_column_edits_no_value_and_removes_no_row():
    """§0.2, and at this stage the temptation is strong: the natural way to write the loop is to
    attach `w` and each outcome to one frame and group. [§10] refits in every replicate, and a frame
    carrying the point fit's weights, resampled, is a replicate weighted by the wrong score."""
    df = secondary_cohort()
    before = df.copy()
    ps = hand_propensity(df)
    e_before, w_before, mask_before = ps.e.copy(), ps.w.copy(), ps.in_model.copy()
    outcome.secondary(df, ps, data.Audit(hand_source(Path("."), "stage9")))
    pd.testing.assert_frame_equal(df, before)
    pd.testing.assert_series_equal(ps.e, e_before)
    pd.testing.assert_series_equal(ps.w, w_before)
    pd.testing.assert_series_equal(ps.in_model, mask_before)


# --- 15.12  Stage 9 adds nothing, refits nothing, and breaks nothing earlier ----------------------------

def test_the_two_stages_are_INDEPENDENT_and_the_ORDER_DOES_NOT_MATTER():
    """`primary` returns the same `Primary` before and after `secondary` runs on the same frame, and
    the reverse. The two halves share `weighted_proportion` and nothing else."""
    df, _, est, _ = estimated()
    _, ps, audit = built(ordinal_cohort())[0], None, None
    df2, audit2 = built(ordinal_cohort())
    ps2 = propensity.fit(df2, audit2)
    first = outcome.primary(df2, ps2, audit2)
    again = outcome.primary(df2, ps2, audit2)
    assert (first.beta, first.odds_ratio) == (again.beta, again.odds_ratio)
    assert first.rd == again.rd


def test_secondary_is_DETERMINISTIC_and_holds_no_seed_and_reads_no_clock(stage9_names):
    """Two runs over one frame agree bit for bit, and the source names neither a generator nor a
    clock. Stage 10 is where randomness enters, and it enters through the RESAMPLER."""
    df, ps, first, _ = secondary_run()
    _, _, again, _ = secondary_run(df, ps)
    for key in config.BINARY_OUTCOMES:
        assert first.estimates[key].rd == again.estimates[key].rd
        assert first.estimates[key].odds_ratio == again.estimates[key].odds_ratio
        assert first.estimates[key].augmented == again.estimates[key].augmented
    assert not {"random", "default_rng", "Generator", "now", "time", "datetime",
                "SEED", "N_BOOT"} & stage9_names


def test_the_STAGE_8_FUNCTIONS_ARE_UNTOUCHED_by_this_stage():
    """§15.12. `weighted_proportion`, `cumulative_rd` and `primary` are read, not edited — and the
    Stage 8 half of this suite passing unedited is the assertion that matters. This one pins the
    signatures, which is what an edit would move first."""
    assert list(inspect.signature(outcome.weighted_proportion).parameters) == ["x", "a", "w"]
    assert list(inspect.signature(outcome.cumulative_rd).parameters) == ["y", "a", "w"]
    assert list(inspect.signature(outcome.primary).parameters) == ["df", "ps", "audit"]
    # `collect` is Stage 10 §13's second licensed change to this module and the only addition to
    # this signature since Stage 9; `paths` is Stage 9's own. Both are keyword-defaulted, so every
    # Stage 1-9 call site is unchanged, and `primary` above is untouched.
    assert list(inspect.signature(outcome.secondary).parameters) == [
        "df", "ps", "audit", "paths", "collect"]


def test_each_estimate_is_INDEPENDENT_OF_THE_OTHER_SIX():
    """`secondary` computing a design per outcome is what makes each estimate independent, and §2
    records that reusing one design across the four full-list outcomes is a STAGE 10 change and not
    a Stage 9 one. This is that independence asserted: blanking one outcome's column entirely moves
    only that outcome's row."""
    df, ps, full, _ = secondary_run()
    partial = df.copy()
    partial.loc[partial.index[:3], "death_90d"] = np.nan
    _, _, changed, _ = secondary_run(partial, ps)
    assert changed.estimates["death_90d"].rd != full.estimates["death_90d"].rd
    for key in config.BINARY_OUTCOMES:
        if key == "death_90d":
            continue
        assert changed.estimates[key].rd == full.estimates[key].rd
        assert changed.estimates[key].augmented == full.estimates[key].augmented


# --- 15.13  "model-assisted", NEVER "doubly robust" `[roadmap]` -----------------------------------------
#
# A test of a WORD, which no behavioural test can be — the failure mode is a docstring, and [§8]'s
# reason for the prohibition is that the word licenses a conclusion the estimator does not support
# (§8.3). The ATO estimand is itself indexed by the true propensity score, so a misspecified `e(X)`
# makes the ESTIMAND wrong even when `m_a(X)` is exactly right; §15.10 is that as a measurement.
#
# **The scope is every Stage 9 deliverable that ships text**, which is wider than the shipped modules:
# `fixtures_stage9.py`, the three amended test modules and the R oracle are Stage 9 deliverables too,
# and §19b's own prose discusses that script in AIPW terms — so the omission would not have been
# hypothetical.

FORBIDDEN_PHRASES: Final[tuple[str, ...]] = (
    "doubly robust", "doubly-robust", "aipw", "dr estimator")   # NOT descriptions of this estimator

SCANNED_FOR_THE_WORD: Final[tuple[str, ...]] = (
    "outcome.py", "model.py", "config.py",
    "tests/fixtures_stage9.py", "tests/test_outcome.py", "tests/test_model.py",
    "tests/test_config.py", "tests/reference/aug_psweight.R")


def _exempt(line: str) -> bool:
    """A line that NAMES the prohibition, quotes it, or explains what the estimator is NOT is exempt.

    §8.3 bans the phrase "in code, docstrings, audit output and this document" — and the ban cannot
    forbid the NAME OF THE THING THE ESTIMATOR IS NOT without making itself unstatable. The exemption
    is therefore a judgement rather than a rule a scanner discovers, and it is written out here so
    that it is reviewable: a line carrying a negation or a §8.3 citation is describing the
    prohibition; a line without one is describing this estimator.
    """
    return any(marker in line for marker in
               ("never", "Never", "NEVER", "not ", "NOT ", "no ", "No ", "forbid", "Forbid",
                "FORBIDDEN", "§8.3", "prohibit", "Prohibit"))


@pytest.mark.parametrize("relative", SCANNED_FOR_THE_WORD)
def test_no_shipped_or_test_text_calls_this_estimator_DOUBLY_ROBUST(relative):
    path = MODULE_DIR / relative
    if not path.exists():
        pytest.skip(f"{relative} is not present in this checkout")
    offenders = [
        f"{relative}:{n}: {line.strip()}"
        for n, line in enumerate(path.read_text(encoding="utf-8").splitlines(), 1)
        if any(phrase in line.lower() for phrase in FORBIDDEN_PHRASES) and not _exempt(line)]
    assert not offenders, (
        "[§8] forbids describing this estimator as doubly robust, because the word licenses a "
        "conclusion it does not support (§8.3):\n  " + "\n  ".join(offenders))


def test_the_word_scan_ACTUALLY_FIRES():
    """A scan that silently matches nothing would otherwise pass as a green test — and this one has an
    exemption rule, which is a second way for it to match nothing."""
    # built from FORBIDDEN_PHRASES rather than written out, so this test is not itself an offender
    assert not _exempt(f"the estimator is {FORBIDDEN_PHRASES[0]}")
    assert _exempt(f"the estimator is NEVER {FORBIDDEN_PHRASES[0]}")
    sample = f"this is an {FORBIDDEN_PHRASES[2].upper()} estimator"
    assert any(phrase in sample.lower() for phrase in FORBIDDEN_PHRASES)


def test_model_assisted_IS_PRESENT_where_the_estimator_is_described(stage9_source):
    """The other half: the licensed description has to actually appear, in the code and in the
    rendered audit text, or the prohibition is satisfied by saying nothing at all."""
    assert "model-assisted" in stage9_source.lower()
    _, _, _, audit = secondary_run()
    rendered = audit.to_markdown()
    assert "model-assisted" in rendered.lower()
    assert not any(phrase in rendered.lower() for phrase in FORBIDDEN_PHRASES)


# --- 15.14  the three audit entries, and the log --------------------------------------------------------
#
# No new kind: `data.KINDS` stays nine, which is `data.py`'s own declaration honoured for the third
# stage running — Stage 7 added two entries under the existing `model` kind, Stage 8 three, Stage 9
# three, and none of the three touched `data.py`.
#
# **THE STRUCTURAL ASSERTION IS ON THE RENDERED TEXT AND NOT VIA THE HASH, and that distinction is the
# whole reason it exists.** `data._md_table` does no escaping and sizes its separator from
# `len(rows[0])`, so a literal `|` inside any cell silently adds markdown columns TO THAT ROW ALONE and
# the table renders crooked. Byte-identity across hash seeds does NOT catch it: a table that is
# identically broken under both seeds is still byte-identical. Stage 7 shipped this exact defect twice
# — `|SMD| < 0.1` and `worst |SMD|`, rendering 8 cells against a 6-cell separator and 15 against 13 —
# and §10.4 was about to ship it a third time as `max|beta|`, which presents ELEVEN header cells
# against a NINE-cell separator. The ASCII name `max_abs_beta` is what makes the assertion pass; the
# assertion is what makes the name load-bearing rather than a preference.

STAGE_9_STEPS: Final[tuple[str, ...]] = (
    "binary_estimable", "binary_estimates", "binary_outcome_models")


def stage9_entries(audit) -> list:
    return [e for e in audit.entries if e.step in STAGE_9_STEPS]


def test_THREE_ENTRIES_of_the_MODEL_kind_in_this_order_at_positions_captured_BEFORE_the_call():
    df = secondary_cohort()
    ps = hand_propensity(df)
    audit = data.Audit(hand_source(Path("."), "stage9"))
    before = len(audit.entries)
    outcome.secondary(df, ps, audit)
    added = audit.entries[before:]
    assert [e.step for e in added] == list(STAGE_9_STEPS)
    assert {e.kind for e in added} == {"model"}


def test_KINDS_stays_NINE_and_no_new_heading_is_declared():
    """§10.1, and `data.py:151-152` predicted it: "Stages 7-13 all render under `model` or under an
    existing kind". Third stage running without touching `data.py`."""
    assert len(data.KINDS) == 9
    assert "model" in data.KINDS


def test_the_three_entries_n_ARE_THE_THREE_DIFFERENT_THINGS_THEY_DESCRIBE():
    """`Audit.record`'s `n` is one number and each entry counts a different thing, which is why they
    disagree: the cohort's records, the number of estimates, and the number of models FITTED."""
    df, _, sec, audit = secondary_run()
    entries = {e.step: e for e in stage9_entries(audit)}
    assert entries["binary_estimable"].n == len(df)                       # the COHORT's
    assert entries["binary_estimates"].n == len(config.BINARY_OUTCOMES)   # seven, not a population
    assert entries["binary_outcome_models"].n == sum(
        1 for e in sec.estimates.values() if e.fit is not None)


def test_binary_estimable_carries_ONE_ROW_PER_OUTCOME_in_registry_order():
    df, ps, sec, audit = secondary_run()
    entry = next(e for e in stage9_entries(audit) if e.step == "binary_estimable")
    assert entry.table[0] == ("outcome", "family", "n [§11]", "lost", "minority", "path", "m_a(X)")
    assert [row[0] for row in entry.table[1:]] == list(config.BINARY_OUTCOMES)
    rows = rows_of(entry)
    for key, est in sec.estimates.items():
        n = int(est.in_estimate.sum())
        assert rows[key][:5] == (est.family, str(n), str(int(ps.in_model.sum()) - n),
                                 str(est.minority), est.augmented_path)


def test_LOST_is_counted_over_the_ATO_POPULATION_and_NOT_over_the_cohort():
    """§10.2, and the two differ whenever a covariate-incomplete record is also outcome-missing.

    `lost` is `|in_model| − |in_model & notna|` and it is NOT `df[key].isna().sum()`. On a frame
    where they agree the wrong expression passes unnoticed, so the frame here is one where they do
    not: a record OFF `in_model` whose TICI is also absent is counted by the cohort-wide expression
    and by neither of the two masks.
    """
    df = secondary_cohort()
    ps = hand_propensity(df, out=2)
    df.loc[df.index[-1], "tici_2b_3"] = np.nan          # off in_model AND outcome-missing
    _, _, sec, audit = secondary_run(df, ps)
    rows = rows_of(next(e for e in stage9_entries(audit) if e.step == "binary_estimable"))

    lost = int(rows["tici_2b_3"][2])
    assert lost == int(ps.in_model.sum()) - int(sec.estimates["tici_2b_3"].in_estimate.sum())
    assert lost != int(df["tici_2b_3"].isna().sum())    # the cohort-wide count, and it disagrees


def test_binary_estimable_names_NO_CASES_and_the_record_TICI_loses_is_NOT_named():
    """§10.2. `_MUST_NAME_CASES` does not cover the `model` kind and this entry removes no patient
    from anything — the seven masks are subsets of `in_model`, whose exclusion `propensity.py`
    already named. An entry naming the same excluded patient a second time is not a second fact.

    What an entry DELIBERATELY OMITS needs asserting as much as what it carries (Stage 8 §9.3), so
    the outcome-missing record is asserted absent by name rather than the tuple merely being empty.
    """
    df, _, _, audit = secondary_run()
    missing_case = df.loc[df["tici_2b_3"].isna(), "case_id"].iloc[0]
    for entry in stage9_entries(audit):
        assert entry.case_ids == ()
    rendered = data._md_table(
        next(e for e in stage9_entries(audit) if e.step == "binary_estimable").table)
    assert missing_case not in rendered


def test_binary_estimates_prints_BOTH_RAW_PROPORTIONS_and_puts_RD_and_TAU_IN_ADJACENT_COLUMNS():
    """[§8]: "the augmented estimate is therefore reported as model-assisted, always alongside the
    unaugmented weighted risk difference; material disagreement between the two is evidence about the
    outcome model." ADJACENCY is what makes the comparison available without arithmetic."""
    _, _, sec, audit = secondary_run()
    entry = next(e for e in stage9_entries(audit) if e.step == "binary_estimates")
    header = list(entry.table[0])
    assert header == ["outcome", "n", "p1_w", "p0_w", "RD_w", "OR_w", "corrected",
                      "tau (model-assisted)", "m_a(X)"]
    assert header.index("tau (model-assisted)") == header.index("RD_w") + 3
    rows = rows_of(entry)
    for key, est in sec.estimates.items():
        assert rows[key][1] == data._fmt(est.proportion[1])
        assert rows[key][2] == data._fmt(est.proportion[0])
        assert rows[key][3] == data._fmt(est.rd)


def test_TAU_RENDERS_AS_A_DASH_AND_NEVER_AS_ZERO_for_an_unaugmented_outcome():
    """Stage 6 §9's rule: the deliberateness of an absence lives in a mask, never in a value. An
    unaugmented outcome has no tau, which is different from a tau nobody could compute — and `_fmt`
    renders the second as `missing`, so the table supplies the dash for the first."""
    _, _, sec, audit = secondary_run()
    rows = rows_of(next(e for e in stage9_entries(audit) if e.step == "binary_estimates"))
    unaugmented = [k for k, e in sec.estimates.items() if e.augmented is None]
    assert unaugmented
    for key in unaugmented:
        assert rows[key][6] == "--"
        assert rows[key][6] not in ("0", "0.0", "0.000000", "missing")
    for key, est in sec.estimates.items():
        if est.augmented is not None:
            assert rows[key][6] == data._fmt(est.augmented)


def test_corrected_RENDERS_FOR_EVERY_ROW_even_when_it_is_False_for_all_seven():
    """A column that disappears when nothing triggers it is a column a reader cannot tell was
    checked (§10.3). On every frame in this suite it is False for all seven."""
    _, _, sec, audit = secondary_run()
    rows = rows_of(next(e for e in stage9_entries(audit) if e.step == "binary_estimates"))
    assert not any(e.or_corrected for e in sec.estimates.values())
    assert [rows[k][5] for k in config.BINARY_OUTCOMES] == ["False"] * 7


def test_A_REDUCED_SPECIFICATION_RENDERS_AS_ITS_COVARIATE_NAMES_and_NOT_as_True():
    """The amendment: "Every reduced specification is reported beside its estimate." A `True` is not
    the specification.

    And the column is `m_a(X)` and not `reduced`: `reduced` is the name of a boolean FIELD on
    `BinaryEstimate` and this column renders a tuple of NAMES. A column and a field sharing a name
    while carrying different types is the kind of collision that survives review.
    """
    _, _, sec, audit = secondary_run()
    entry = next(e for e in stage9_entries(audit) if e.step == "binary_estimates")
    assert entry.table[0][-1] == "m_a(X)" and "reduced" not in entry.table[0]
    rows = rows_of(entry)
    assert rows["tici_2b_3"][-1] == "center, atrial_fib"
    assert rows["tici_2b_3"][-1] not in ("True", "False")
    for key, est in sec.estimates.items():
        if est.augmented_path == "unaugmented":
            assert rows[key][-1] == "--"


def test_binary_outcome_models_has_ONE_ROW_PER_AUGMENTED_OUTCOME_and_not_seven():
    """The two unaugmented outcomes have no model and a row asserting that would be a row about
    nothing."""
    _, _, sec, audit = secondary_run()
    entry = next(e for e in stage9_entries(audit) if e.step == "binary_outcome_models")
    fitted = [k for k, e in sec.estimates.items() if e.fit is not None]
    assert [row[0] for row in entry.table[1:]] == fitted
    assert len(entry.table) - 1 == len(fitted) < len(config.BINARY_OUTCOMES)


def test_max_abs_beta_IS_SPELLED_IN_ASCII_and_carries_the_three_SAFEGUARD_COUNTERS():
    """`max|beta|` would put TWO PIPES in a header cell (§10.4). `max_abs_beta` is in the table for
    §9.6's reason: no bound rejects a fit, so it is the only thing standing between a separated
    nuisance model and an unremarked estimate. The counters are there for Stage 6 §3.2's — "a
    safeguard whose activation nobody can count is a safeguard nobody can evaluate"."""
    _, _, sec, audit = secondary_run()
    entry = next(e for e in stage9_entries(audit) if e.step == "binary_outcome_models")
    assert entry.table[0] == ("outcome", "covariates", "k", "rows", "iters", "max_abs_beta",
                              "rescales", "halvings", "dropped")
    assert "|" not in "".join(entry.table[0])
    rows = rows_of(entry)
    for key, est in sec.estimates.items():
        if est.fit is None:
            continue
        assert rows[key][4] == data._fmt(float(np.max(np.abs(est.fit.beta))))
        assert rows[key][5:7] == (str(est.fit.rescales), str(est.fit.halvings))
        assert rows[key][7] == "center_USZ"


@pytest.mark.parametrize("step", STAGE_9_STEPS)
def test_STRUCTURAL_PER_ENTRY_ON_THE_RENDERED_TEXT_equal_cells_and_no_pipe_in_any_cell(step):
    """The assertion that byte-identity cannot make, and Stage 7 shipped the defect twice for want of
    it. Header cells == separator cells == every body row's cells, on the RENDERED string."""
    _, _, _, audit = secondary_run()
    entry = next(e for e in stage9_entries(audit) if e.step == step)

    assert not [c for row in entry.table for c in row if "|" in c]
    lines = data._md_table(entry.table).splitlines()
    counts = {n: len(line.split("|")) - 2 for n, line in enumerate(lines)}
    assert len(set(counts.values())) == 1, (
        f"{step} renders {sorted(set(counts.values()))} markdown cells across its rows; "
        f"`data._md_table` sizes its separator from len(rows[0]) and does no escaping:\n"
        + "\n".join(lines))
    assert set(counts.values()) == {len(entry.table[0])}


def test_THE_PIPE_ASSERTION_ACTUALLY_FIRES_on_a_table_that_carries_one():
    """Stage 7's shipped defect, reproduced, so the assertion above is known to be able to fail.

    `max|beta|` in a header cell presents ELEVEN markdown cells against a NINE-cell separator — and
    the rendered table is byte-identical under both hash seeds while being wrong under both.
    """
    header = ("outcome", "covariates", "k", "rows", "iters", "max|beta|",
              "rescales", "halvings", "dropped")
    body = ("mrs_0_2_90d", "the [§6] set", "12", "92", "6", "3.0933", "0", "0", "center_USZ")
    lines = data._md_table((header, body)).splitlines()
    counts = [len(line.split("|")) - 2 for line in lines]
    assert counts == [11, 9, 9]                 # header, separator, body — and the table is crooked
    assert len(set(counts)) > 1


def test_THE_LOG_IS_BYTE_IDENTICAL_ACROSS_TWO_HASH_SEEDS():
    """Which is where the seven-outcome dict's iteration order is actually asserted.

    Every table here ranges over `C.BINARY_OUTCOMES` rather than over `estimates.keys()`, so the row
    order is the registry's. A rendering that iterated the dict would be reproducible in CPython by
    accident and this is what would notice if it stopped being.
    """
    driver = textwrap.dedent("""
        import sys
        sys.path.insert(0, {module_dir!r})
        sys.path.insert(0, {tests_dir!r})
        import data, outcome
        from test_outcome import secondary_cohort, hand_propensity
        from test_data import hand_source
        from pathlib import Path
        df = secondary_cohort()
        audit = data.Audit(hand_source(Path("."), "stage9"))
        outcome.secondary(df, hand_propensity(df), audit)
        sys.stdout.write(audit.to_markdown())
    """).format(module_dir=str(MODULE_DIR), tests_dir=str(TESTS_DIR))

    rendered = []
    for seed in ("0", "1"):
        completed = subprocess.run([sys.executable, "-c", driver], capture_output=True, text=True,
                                   env={**os.environ, "PYTHONHASHSEED": seed}, check=True)
        rendered.append(completed.stdout)
    assert rendered[0] == rendered[1]
    assert "binary_outcome_models" in rendered[0]


def test_proportion_CANNOT_BE_MUTATED_THROUGH_THE_RETURNED_ESTIMATE():
    """`frozen=True` stops attribute REBINDING and does nothing about a `dict` FIELD.

    §11 stores `dict(share)` rather than the object `weighted_proportion` returned, so a caller who
    writes into `estimate.proportion` changes their own copy and nothing else. Without the copy the
    mutation would reach a frozen estimate AND — before the copy — the dict `marginal_odds_ratio`
    had also returned to a different caller.
    """
    _, _, sec, audit = secondary_run()
    est = sec.estimates["mrs_0_2_90d"]
    rendered_before = data._md_table(
        next(e for e in stage9_entries(audit) if e.step == "binary_estimates").table)
    rd_before, p_before = est.rd, dict(est.proportion)

    est.proportion[1] = 99.0                              # succeeds: a dict is mutable regardless
    assert est.rd == rd_before
    rendered_after = data._md_table(
        next(e for e in stage9_entries(audit) if e.step == "binary_estimates").table)
    assert rendered_after == rendered_before

    with pytest.raises(dataclasses.FrozenInstanceError):
        est.rd = 99.0                                     # the attribute itself IS frozen
    est.proportion.update(p_before)


def test_THE_COMPANION_storing_the_returned_dict_DIRECTLY_lets_a_caller_change_TWO_things():
    """Why §11 writes `dict(share)` and not `share`. Without the copy, one caller's mutation reaches
    a frozen estimate AND the dict `marginal_odds_ratio` returned to another caller: they are the
    same object."""
    y, a, w, _ = golden_arrays()
    _, _, share = outcome.marginal_odds_ratio(y, a, w)
    aliased, copied = share, dict(share)
    aliased[1] = 99.0
    assert share[1] == 99.0                               # the caller's dict IS the returned one
    assert copied[1] != 99.0                              # and §11's copy is not


@DATA_GATED
def test_the_workbook_ledger_reaches_THIRTY_after_the_full_pipeline(workbook_stage9):
    """`[data-gated]`. §10.1: load 7 / derive 4 / classify 1 / build 6 / fit 4 / assess 2 /
    primary 3 / secondary 3 = 30, reconciling with Stage 8's 24-after-`assess` and 27-after-`primary`.
    """
    df, ps, sec, audit = workbook_stage9
    assert len(audit.entries) == 30
    assert len(data.KINDS) == 9


@DATA_GATED
def test_the_workbook_fits_FIVE_nuisance_models_and_not_seven(workbook_stage9):
    """`[data-gated]`. Four full plus the one reduced; `sich` and `ph2` have none."""
    df, ps, sec, audit = workbook_stage9
    entry = next(e for e in audit.entries if e.step == "binary_outcome_models")
    assert entry.n == 5
    assert len(entry.table) - 1 == 5
    assert [row[0] for row in entry.table[1:]] == [
        "mrs_0_2_90d", "mrs_0_1_90d", "tici_2b_3", "death_90d", "mrs_5_6_90d"]


@DATA_GATED
def test_the_workbook_drops_center_USZ_from_EVERY_design(workbook_stage9):
    """`[data-gated]`, a structural fact about the cohort rather than an estimate: USZ contributes no
    records to the ATO population, so `design`'s constant-column rule fires on every one of the five
    fitted designs. §7.2's five-parameter arithmetic for TICI depends on it."""
    df, ps, sec, audit = workbook_stage9
    for est in sec.estimates.values():
        if est.fit is not None:
            assert est.dropped == ("center_USZ",)
    assert len(sec.estimates["tici_2b_3"].fit.beta) == 5

# ======================================================================================================
# Stage 10 §15.7 — S8's reclassification and `secondary`'s `collect`
# ======================================================================================================
#
# These live here and not in `test_bootstrap.py` because that is where the changed function is
# (Stage 10 §15). Two shipped changes are under test and they are the only two Stage 10 makes to
# `outcome.py` (Stage 10 §13): S8 raises `model.FitError` rather than `C.SchemaError`, and
# `secondary` gains `collect`.
#
# The frames are `fixtures_stage10.py`'s, which is Stage 10 §15.0's own rule — every constructed
# fixture is code, so a pin nobody can reproduce is not a pin.

from fixtures_stage10 import (constant_outcome_frame, two_centre_frame,   # noqa: E402
                              unfittable_nuisance_frame)


def stage10_fit(df: pd.DataFrame) -> propensity.Propensity:
    """`propensity.fit` over a `fixtures_stage10` frame, on a throwaway `Audit`.

    FITTED and not hand-assembled, unlike §15.0's `hand_propensity`: these frames exist to be driven
    end to end through the stage as a [§10] replicate drives it, and a hand-built `Propensity` would
    make the `m_a(X)` failure `unfittable_nuisance_frame` exists to produce a property of an invented
    `e` rather than of the frame.
    """
    return propensity.fit(df, data.Audit(hand_source(Path("."), "stage10")))


@pytest.mark.parametrize("key", ["tici_2b_3", "sich"])       # one secondary, one safety
def test_S8_raises_FitError_and_DOES_NOT_raise_SchemaError(key):
    """DECISION 6 (PI, 2026-08-25), Stage 10 §5.4. **The second half is the assertion that fails on
    the landed code**, and it is the one that has to exist: `pytest.raises(model.FitError)` alone
    would pass on an implementation that raised both, and there is no such thing — but it would also
    pass if `FitError` were ever made a subclass of `SchemaError`, which is the collapse Stage 10
    §15.4 asserts against from the other side.

    Parametrised over one outcome from each [§5] family, because S8 fires inside the per-outcome loop
    and a test on one family would not notice a guard keyed on `family`.
    """
    df = constant_outcome_frame(key)
    ps = stage10_fit(df)
    with pytest.raises(model.FitError) as e:
        outcome.secondary(df, ps, data.Audit(hand_source(Path("."), "stage10")))
    assert "S8" in message_of(e) and key in message_of(e)
    assert not isinstance(e.value, config.SchemaError)


@pytest.mark.parametrize("branch", ["S6", "S7"])
def test_S6_and_S7_STILL_raise_SchemaError_after_S8_moved(branch):
    """§5.4 changes ONE precondition of three, and a change that moved all three would pass a test
    written only for S8. S6 is an empty [§11] population and S7 an arm with no record or no weight;
    Stage 9 §4.4 argues both mean the mask and the weights disagree about the same rows, which is a
    bug and not a sparse replicate — so both stay uncatchable by [§10] (Stage 10 §5.4, §7.1)."""
    df, ps = broken_secondary(branch)
    with pytest.raises(config.SchemaError) as e:
        outcome.secondary(df, ps, data.Audit(hand_source(Path("."), "stage9")))
    assert branch in message_of(e)
    assert not isinstance(e.value, model.FitError)


def test_the_S8_frame_is_CONSTANT_and_not_merely_rare():
    """The fixture's own witness first, as §15.5's separated frame gets one: the outcome really does
    take one value across its whole [§11] population, so the raise is S8's and not S7's."""
    df = constant_outcome_frame("sich")
    ps = stage10_fit(df)
    population = ps.in_model & df["sich"].notna()
    assert int(population.sum()) > 0
    assert df.loc[population, "sich"].nunique() == 1


def test_collect_False_RAISES_on_an_unfittable_nuisance_model():
    """The point-estimate contract unchanged (Stage 10 §7.3). On the workbook a `FitError` is fatal
    and must be, because there is no replicate to drop."""
    df = unfittable_nuisance_frame()
    ps = stage10_fit(df)
    with pytest.raises(model.FitError):
        outcome.secondary(df, ps, data.Audit(hand_source(Path("."), "stage10")))


def test_collect_True_returns_SIX_estimates_and_ONE_RECORDED_failure():
    """**Both halves, and the second is what makes it a test**: an implementation that dropped the
    key entirely would give six estimates too (Stage 10 §15.7).

    "Exactly one" is the fixture's specification: a frame on which two outcomes failed would pass a
    wrong implementation that gives up after the first.
    """
    df = unfittable_nuisance_frame()
    ps = stage10_fit(df)
    sec = outcome.secondary(df, ps, data.Audit(hand_source(Path("."), "stage10")), collect=True)
    assert len(sec.estimates) == 6
    assert set(sec.failures) == {"tici_2b_3"}
    assert set(sec.estimates) | set(sec.failures) == set(config.BINARY_OUTCOMES)
    assert not set(sec.estimates) & set(sec.failures)
    # The MESSAGE is recorded and not merely the key, because Stage 10 §7.2 classifies the failure
    # by the leading token of it and a bare absence carries no token.
    assert sec.failures["tici_2b_3"].split()[0] == "F3"


def test_collect_True_leaves_by_family_one_short_rather_than_raising_KeyError():
    """`by_family` ranges over the registry, so a dropped outcome is a missing member of its family
    and not a broken registry (Stage 10 §7.3)."""
    df = unfittable_nuisance_frame()
    sec = outcome.secondary(df, stage10_fit(df),
                            data.Audit(hand_source(Path("."), "stage10")), collect=True)
    families = sec.by_family()
    assert sum(len(v) for v in families.values()) == 6
    assert "tici_2b_3" not in [e.outcome for group in families.values() for e in group]


@pytest.mark.parametrize("collect", [False, True])
def test_collect_NEVER_catches_SchemaError_at_EITHER_setting(collect):
    """S6's frame trips a `SchemaError` inside the per-outcome loop, which is exactly where
    `collect`'s `except` sits — so this is the assertion that the `except` clause names
    `model.FitError` and not `Exception` (Stage 10 §7.1, §15.7)."""
    df, ps = broken_secondary("S6")
    with pytest.raises(config.SchemaError):
        outcome.secondary(df, ps, data.Audit(hand_source(Path("."), "stage9")), collect=collect)


def test_the_FROZEN_paths_still_govern_under_collect_True():
    """An outcome whose frozen path is `"unaugmented"` fits no `m_a(X)` and so cannot contribute a
    nuisance `FitError` at all — which is §6.3's rule and `collect` not interacting with it.

    Asserted by freezing `tici_2b_3` to `"unaugmented"` on the very frame whose `tici_2b_3` nuisance
    model is unfittable: with the path frozen there is no fit to fail, so all seven outcomes come
    back estimated and `failures` is empty. The same frame with the paths decided from itself gives
    six and one, which is the test immediately above.
    """
    df = unfittable_nuisance_frame()
    ps = stage10_fit(df)
    paths = {key: "unaugmented" if key == "tici_2b_3" else
             outcome._augmented_path(key, min(
                 int((df.loc[ps.in_model & df[key].notna(), key] == 1.0).sum()),
                 int((df.loc[ps.in_model & df[key].notna(), key] == 0.0).sum())))
             for key in config.BINARY_OUTCOMES}
    sec = outcome.secondary(df, ps, data.Audit(hand_source(Path("."), "stage10")),
                            paths=paths, collect=True)
    assert sec.failures == {}
    assert len(sec.estimates) == 7
    assert sec.estimates["tici_2b_3"].augmented is None
    assert sec.estimates["tici_2b_3"].fit is None


def test_collect_False_is_the_DEFAULT_and_the_point_estimate_log_is_UNCHANGED():
    """`collect` defaults to False and a clean frame renders the same three entries either way.

    The byte comparison is what protects Stage 9's audit assertions: `_estimates_table` and
    `_models_table` gained a skip for an absent key, and on a frame where no key is absent the two
    tables must be what they were (Stage 10 §13's "no other function changes").
    """
    assert inspect.signature(outcome.secondary).parameters["collect"].default is False
    df = two_centre_frame()
    ps = stage10_fit(df)
    a_default = data.Audit(hand_source(Path("."), "stage10"))
    a_collect = data.Audit(hand_source(Path("."), "stage10"))
    outcome.secondary(df, ps, a_default)
    outcome.secondary(df, ps, a_collect, collect=True)
    assert a_default.to_markdown() == a_collect.to_markdown()


# --- the locked-estimate digests `[data-gated]` ---------------------------------------------------------
#
# **THE REGRESSION PIN, IN THE FORM THAT PUTS NO ESTIMATE INTO VERSION CONTROL.**
#
# `TODOS.md`'s pin item and Stage 8 §4.3 between them left a real hole: there is no pinned regression
# number on the primary effect or on the seven binary estimates anywhere in git, so an edit that moves
# `β` in the fourth decimal without breaking convergence, ordering, orientation or monotonicity passes
# the whole suite. §14.0.2's and §15.0.2's golden vectors are a partial mitigation only — 5 and 24
# records against the workbook's 92 — and cannot witness a defect that appears only at the workbook's
# scale. `TODOS.md`'s own open item to make `POLR_TOL` relative to the weight total is a concrete
# change of exactly that shape.
#
# **The pin is a CONTENT HASH and not a literal, which is `config.DATA_SHA256`'s pattern applied one
# stage on.** That module pins the workbook by digest and `test_config.py` asserts it; this pins the
# locked estimates the same way. Four things follow, and the third is why this form was chosen over
# writing the numbers out:
#
#   * It fires on ANY drift at the pinned precision, which is stronger than a hand-copied literal.
#   * It localises: two digests, so a failure names Stage 8 or Stage 9 rather than "something moved".
#   * **No estimate enters version control**, so [§15]'s rule — no number on the primary effect pinned
#     before the analysis is locked — is not spent, it is honoured. The pin therefore did not need the
#     lock and could have existed since Stage 8; that it did not is because "pin" was read as "write
#     the number down".
#   * §14.12a's float-literal scan needs no change, because there is no float literal to scan. Its
#     scope is §14.12 alone and this section is outside it — asserted below, so that stays true.
#
# **The precision is `data._fmt`'s and that is deliberate.** Six significant figures is the precision
# the audit log publishes and the precision Stage 2's byte-identity criterion is written against, so
# the pin is tied to the reproducibility contract the project already has rather than to a second one.
# It is also what keeps the digest immune to a last-ULP change from a numpy patch release inside
# `pyproject.toml`'s declared range, while still catching a fourth-decimal move in `β`.
#
# **When one of these fails.** It does not say what moved. Run the pipeline, and diff
# `out/logs/audit_<label>.md` against the values `../out/stage0_data_inventory.md` records beside these
# digests — that register is the baseline, because the log is overwritten on every run. Regenerating a
# digest is correct ONLY if the change was intended and the register records why; a digest updated to
# make a test pass is the whole of this section undone.

_LOCKED_PRIMARY_DIGEST: Final[str] = (
    "c955301b8d45c1f7178de9453782ceee770599cb7d84a08c1609d7cd4ea320b3")
_LOCKED_SECONDARY_DIGEST: Final[str] = (
    "a4ce1868b85ebf1e02f951fe84eb295dd19eab7112359033ed74b67b700d2bc4")

# The workbook's primary rendering is 26 lines. A COUNT and not an estimate, in §14.12's own sense —
# 2 scalars, 6 cutpoints because all seven declared mRS levels are occupied, 6 risk differences and
# 12 cumulative probabilities. It is pinned because the derived expectation below is a function of
# `len(est.alpha)`, so a collapse that lost a cutpoint would satisfy the derivation while changing
# what is being pinned.
_WORKBOOK_PRIMARY_LINES: Final[int] = 26


def _expected_primary_lines(est: outcome.Primary) -> int:
    """The rendering's shape, DERIVED from the declarations and the fit rather than pinned.

    A renderer bug that emitted nothing would produce a stable digest of the empty string, and whoever
    regenerated it would pin the hash of nothing — so the shape is asserted before the digest is. It
    is derived and not a constant because the cutpoint count is a property of the SAMPLE: `polr`
    collapses the response to the categories carrying positive weight, so a fixture with four occupied
    levels renders four `alpha` lines and the workbook renders six.
    """
    return (2                                              # beta, odds_ratio
            + len(est.alpha)                               # one per FITTED cutpoint
            + len(config.MRS_THRESHOLDS)                   # RD_k
            + len(config.MRS_THRESHOLDS) * len(config.TREATMENT_LABELS))    # cumulative


def _expected_secondary_lines() -> int:
    """Eight lines per [§5] binary outcome: minority, path, two proportions, rd, or, flag, tau."""
    return len(config.BINARY_OUTCOMES) * (6 + len(config.TREATMENT_LABELS))


def _digest(canonical: str) -> str:
    return hashlib.sha256(canonical.encode("utf-8")).hexdigest()


def canonical_primary(est: outcome.Primary) -> str:
    """The [§8] primary estimate as one deterministic string, every value through `data._fmt`.

    Ranges over `MRS_THRESHOLDS` and `TREATMENT_LABELS` rather than over the dicts' own key order, so
    the rendering is a function of the declarations and not of insertion order.
    """
    lines = [f"beta={data._fmt(est.beta)}", f"odds_ratio={data._fmt(est.odds_ratio)}"]
    lines += [f"alpha[mrs<={level}]={data._fmt(value)}"
              for level, value in zip(est.cut_levels, est.alpha)]
    lines += [f"rd[{k}]={data._fmt(est.rd[k])}" for k in config.MRS_THRESHOLDS]
    lines += [f"cumulative[{k}][{code}]={data._fmt(est.cumulative[k][code])}"
              for k in config.MRS_THRESHOLDS for code in config.TREATMENT_LABELS]
    return "\n".join(lines)


def canonical_secondary(sec: outcome.Secondary) -> str:
    """The seven [§8] binary estimates as one deterministic string. `C.BINARY_OUTCOMES` order."""
    lines: list[str] = []
    for key in config.BINARY_OUTCOMES:
        e = sec.estimates[key]
        lines.append(f"{key}.minority={e.minority}")
        lines.append(f"{key}.path={e.augmented_path}")
        for code in config.TREATMENT_LABELS:
            lines.append(f"{key}.proportion[{code}]={data._fmt(e.proportion[code])}")
        lines.append(f"{key}.rd={data._fmt(e.rd)}")
        lines.append(f"{key}.odds_ratio={data._fmt(e.odds_ratio)}")
        lines.append(f"{key}.or_corrected={e.or_corrected}")
        lines.append(f"{key}.augmented="
                     + ("absent" if e.augmented is None else data._fmt(e.augmented)))
    return "\n".join(lines)


def test_the_canonical_renderings_have_the_PINNED_SHAPE_and_every_line_carries_a_value():
    """Asserted on the FIXTURE cohorts, so it runs with no `data/` — the shape is a property of the
    renderer and of the declarations, not of the workbook.

    Two different fixtures, because no single one reaches both stages: `ordinal_cohort()` is five
    records and `outcome.secondary` raises S8 on it — `tici_2b_3` is constant on its [§11] population
    there — while `secondary_cohort()` is the frame Stage 9's own sections are driven on. Measured.
    """
    _, _, est, _ = estimated()
    _, _, sec, _ = secondary_run()
    for canonical, expected in ((canonical_primary(est), _expected_primary_lines(est)),
                                (canonical_secondary(sec), _expected_secondary_lines())):
        lines = canonical.split("\n")
        assert len(lines) == expected
        assert all("=" in line and line.split("=", 1)[1] for line in lines)
        assert len(set(lines)) == len(lines)          # no key rendered twice


def test_the_DIGEST_IS_SENSITIVE_at_the_precision_it_pins():
    """The mutation companion, and it runs with no `data/`. A pin nobody has seen discriminate is a
    pin nobody can price.

    `beta` is perturbed in the SIXTH significant figure — the last `data._fmt` renders — and the
    digest must move. And it is perturbed in the twelfth, where it must NOT, because that is the
    numpy-patch-release noise the `_fmt` precision exists to absorb.
    """
    df, ps, est, _ = estimated()
    baseline = _digest(canonical_primary(est))
    sixth = dataclasses.replace(est, beta=est.beta * (1.0 + 1e-5))
    twelfth = dataclasses.replace(est, beta=est.beta * (1.0 + 1e-12))
    assert _digest(canonical_primary(sixth)) != baseline
    assert _digest(canonical_primary(twelfth)) == baseline


def test_the_digest_section_is_OUTSIDE_14_12a_scan_so_that_guard_needs_no_change():
    """§14.12a scans §14.12 for float literals including inside docstrings, and it exists to keep the
    primary effect out of git. This section does not touch it, and that is asserted rather than
    assumed: the scan's boundary is §14.12's banner to the next one, and this section is past it."""
    scanned = _section_source(THIS_FILE, SECTION_BANNER)
    assert "_LOCKED_PRIMARY_DIGEST" not in scanned
    assert _unallowlisted_float_literals(scanned) == []
    # ...and the property that makes this form of the pin work at all: NO ESTIMATE APPEARS HERE. It is
    # asserted as an exact set rather than as the absence of magnitudes, because two magnitudes
    # legitimately belong — the mutation companion's perturbation factors, which are sizes of a nudge
    # and not values of anything. A pasted `beta` is a third member and fails.
    section = THIS_FILE[THIS_FILE.index("# --- the locked-estimate digests"):]
    floats = {node.value for node in ast.walk(ast.parse(textwrap.dedent(section)))
              if isinstance(node, ast.Constant) and isinstance(node.value, float)}
    assert floats == {1.0, 1e-5, 1e-12}, f"an unexpected magnitude appears here: {sorted(floats)}"


@DATA_GATED
def test_the_LOCKED_PRIMARY_ESTIMATE_has_not_moved(workbook_primary):
    """`[data-gated]`. The [§8] primary effect, pinned by digest at `data._fmt`'s precision.

    If this fails, something changed `beta`, `exp(beta)`, a cutpoint, a cumulative risk difference or a
    weighted cumulative probability. Diff the regenerated audit log against the values the register
    records beside this digest. **Do not regenerate the digest to make this pass.**
    """
    canonical = canonical_primary(workbook_primary)
    assert len(canonical.split("\n")) == _WORKBOOK_PRIMARY_LINES
    assert len(canonical.split("\n")) == _expected_primary_lines(workbook_primary)
    assert _digest(canonical) == _LOCKED_PRIMARY_DIGEST


@DATA_GATED
def test_the_LOCKED_SECONDARY_ESTIMATES_have_not_moved(workbook):
    """`[data-gated]`. The seven [§8] binary estimates, pinned by digest at the same precision.

    Separate from the primary's so that a failure names the stage. The same instruction applies.
    """
    df, ps, audit = workbook
    canonical = canonical_secondary(outcome.secondary(df, ps, audit))
    assert len(canonical.split("\n")) == _expected_secondary_lines()
    assert _digest(canonical) == _LOCKED_SECONDARY_DIGEST
