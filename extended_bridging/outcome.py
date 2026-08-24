"""Stage 8b — the [§8] primary outcome estimate.

The [§8] specification: a weighted proportional-odds model with treatment as the sole predictor,
oriented so that ``exp(beta) > 1`` favours bridging, and the weighted empirical cumulative risk
differences ``RD_k`` that present the same contrast on an absolute scale. The numerics it stands on
are ``model.py``'s, because [§14a] prescribes the same estimator with a wider design and unit weights
and Stage 12 calls ``model.polr`` directly [Stage 8 §0.1].

**SEPARATION IN THIS ESTIMATOR DOES NOT PRESENT AS NON-CONVERGENCE, AND THAT IS WHY G7 EXISTS.** The
roadmap's *Accept when* reads as though a sparse ordinal fit failed loudly, as an unpenalised
propensity fit does. It does not. Measured: a perfectly separated 20-against-20 frame **converges in
17 iterations on the score criterion**, with every safeguard counter at zero, every fitted probability
finite and nothing anywhere out of range — and returns ``exp(beta) = 6.5e15``. Four hundred sparse
replicates of this cohort's shape produced **zero** convergence failures, because there are none to
produce. So [§10]'s "replicates whose prespecified fit fails are dropped and counted" has nothing to
count, and a percentile interval's upper limit becomes a quantile of a distribution with a tail at
1e+09 on the odds-ratio scale while the p-value barely moves — the interval and the p-value
disagreeing is the one outcome [§10] wrote its p-value definition to prevent. ``_assert_reportable``
is the only thing that turns a degenerate fit into a failure [§10] can count. **It is not a magic
number to be deleted.**

**No standard error, no interval and no p-value leaves this stage, and that is a correctness claim
rather than a division of labour** [Stage 8 §5.6]. The naive weighted-likelihood covariance treats the
overlap weights as frequencies, which they are not: they are a tilting function of an **estimated**
propensity score, so the sampling variability of estimating it enters the estimate and ``-H^-1`` omits
it entirely. [§7] says so directly. [§10]'s percentile bootstrap is the prespecified interval and the
source of the primary p-value.

**This stage adds no column to the cohort frame, edits no value, removes no row and re-weights
nothing.** Stage 6 §0.2 and Stage 7 §0.2 made the same promise; here the temptation is the strongest
yet, because the natural way to write this is to attach ``w`` and the outcome to one frame and group.
[§10] refits in every one of ``N_BOOT`` replicates, and a frame carrying the point fit's weights,
resampled into a replicate, is a replicate weighted by the wrong score.

**Nothing here reads the ``Balance``.** Stage 7 §11 states it and Stage 7 §16 records why: no estimate
in this SAP is conditioned on a balance diagnostic, and a pipeline in which a threshold silently
selects an estimator is the mixture [§7] forbids one level up.

Section references in brackets are to ``statistical_analysis_plan.md``. The specification for this
module is ``specs/stage8_primary_outcome_estimator.md``; nothing here is invented outside it.

::

    cohort.build(df, audit)        →  93 rows x 33 columns              [Stage 5 §11]
    propensity.fit(cohort, audit)  →  Propensity(e, w, in_model, …)     [Stage 6 §11]
    balance.assess(...)            →  Balance — NOT read here           [Stage 7 §11]
                       │
                       ▼
    ┌─────────────────────────────────────────────────────────────────────────────┐
    │  STAGE 8b — outcome.py       reads config, data.Audit, model, propensity    │
    │                                                                             │
    │   weighted_proportion(x, a, w)   the one weighted per-arm mean      (§8.1)  │
    │   cumulative_rd(y, a, w)         RD_k over MRS_THRESHOLDS           (§8)    │
    │                                                                             │
    │   primary(df, ps, audit) -> Primary                                         │
    │     ├─ _assert_primary_inputs    G1,G2,G3 THEN G4,G5 → SchemaError  (§4.4)  │
    │     ├─ in_estimate = in_model & outcome present     the [§11] third (§4.1)  │
    │     ├─ _record_outcome_completeness                 entry 1        (§9.2)  │
    │     ├─ design(sub, (TREATMENT,))  → G6: the exposure SURVIVED       (§4.5)  │
    │     ├─ model.polr(X, y, w)                                          (§5)    │
    │     ├─ _assert_reportable        G7: |beta| bounded → FitError      (§6)    │
    │     ├─ cumulative_rd from the weighted distributions                (§8)    │
    │     └─ three `model` audit entries                                  (§9)    │
    └─────────────────────────────────────────────────────────────────────────────┘
                       │
                       ▼
     Primary(beta, odds_ratio, alpha, cut_levels, rd, cumulative, in_estimate, fit)
                       │
                       ▼
     Stage 10 [§10] refits this in every replicate    Stage 14 [§16] reports it

``primary`` appends three entries of the existing ``model`` kind. It writes no file, as no stage
before it does — and **on this stage the audit log is the only place the primary estimate exists**
[Stage 8 §4.3].

This module is **not** exempt from the Stage 1 §7 raw-name scan and must never become exempt.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Final

import numpy as np
import pandas as pd

import config as C
import model
import propensity
# `_fmt` is private to data.py and is imported anyway, for the reason derive.py, propensity.py and
# balance.py give: it is the pipeline's *one* float formatter, and a second one is a second way for
# two runs to disagree.
from data import Audit, _fmt


# --- the arm codes, named once. §8.1 ---------------------------------------------------------------
#
# [§8]'s RD_k is P(. | bridging) - P(. | EVT alone) and THE ORDER OF THE SUBTRACTION IS THE SIGN,
# which §7 is four hundred words about. `share[1] - share[0]` is correct and is two bare integers in
# the one expression where a reader checking the orientation most needs a name to check against.
# Derived from TREATMENT_LABELS rather than written as 1 and 0, for PRIMARY_OUTCOME's reason: the
# registry already says which code is the treated arm, and a literal here is a second place for it to
# be said. `config.py`'s TREATMENT comment documents the coding as `1 = IVT before EVT (bridging),
# 0 = EVT alone`, so the treated arm is the higher code and that is the property being read.

_COMPARATOR: Final[int] = min(C.TREATMENT_LABELS)
_TREATED: Final[int] = max(C.TREATMENT_LABELS)

# The three audit step names. Named as constants rather than written at the three `audit.record` call
# sites because the acceptance suite asserts the entries appear in this order at positions captured
# before the call, and a test comparing against a string literal it also writes is a test of nothing.
_STEP_COMPLETENESS: Final[str] = "outcome_completeness"
_STEP_FIT: Final[str] = "primary_fit"
_STEP_RD: Final[str] = "cumulative_rd"


# --- what the stage returns -------------------------------------------------------------------------

@dataclass(frozen=True)
class Primary:
    """The [§8] primary estimate. Point estimates only: no interval, no p-value, no standard error.

    `in_estimate` is boolean and TOTAL, exactly as `Propensity.in_model` is, and for Stage 6 §9's
    reason: the deliberateness of an absence lives in a mask, never in a value. Stage 8 adds no `nan`
    column of its own and so has nothing to fill; what it has is a THIRD population, and a caller that
    wants the [§11] denominator reads the mask rather than counting `notna()` on something.

    **There is no `favours_bridging()` method and no verdict of any kind.** Stage 7's `Balance`
    carries `worst()` and `unbalanced()` because [§9] specifies a threshold and a verdict is what a
    threshold is for. [§8] specifies no threshold on `beta`: the null is tested at Stage 10 from the
    bootstrap distribution, and a method here returning `odds_ratio > 1.0` would be a significance
    test with no interval behind it, one attribute access away from a reporting layer.
    """

    beta: float                              # the treatment coefficient — §5.1's parametrisation
    odds_ratio: float                        # exp(beta), oriented so > 1 favours bridging — §7
    alpha: tuple[float, ...]                 # the fitted cutpoints, ascending
    cut_levels: tuple[int, ...]              # the declared mRS level each alpha cuts at or below
    rd: dict[int, float]                     # RD_k, keyed by MRS_THRESHOLDS — §8
    cumulative: dict[int, dict[int, float]]  # P(Y <= k | arm), keyed threshold then arm code — §8
    in_estimate: pd.Series                   # boolean, TOTAL — the [§11] population — §4.1
    fit: model.PolrFit


# --- the weighted per-arm proportion, and RD_k [§8] --------------------------------------------------

def weighted_proportion(x: np.ndarray, a: np.ndarray, w: np.ndarray) -> dict[int, float]:
    """The weighted mean of `x` within each declared arm. Keyed by arm code, as TREATMENT_LABELS is.

    Public, and it is Stage 9's as much as this stage's — the [§8] weighted risk difference is this
    function differenced and the weighted marginal odds ratio is these two proportions — which is why
    it lives here rather than inside `cumulative_rd`. A second implementation in a later stage would
    be a second definition of the denominator, and [§11] requires the denominator. The same move
    Stage 6 §3 made with `ess` and Stage 7 §0.1 made with `smd`.

    `nan` where an arm carries no positive weight, never 0.0. **The check comes BEFORE the mean, and
    the order is load-bearing**: `np.average` with weights summing to zero RAISES `ZeroDivisionError`
    rather than returning nan [Stage 7 §3.1] — the same fact, in the same place, one stage on.
    """
    out: dict[int, float] = {}
    for code in C.TREATMENT_LABELS:
        arm = a == float(code)
        total = float(w[arm].sum())
        out[code] = float(np.average(x[arm], weights=w[arm])) if total > 0.0 else np.nan
    return out


def cumulative_rd(y: np.ndarray, a: np.ndarray,
                  w: np.ndarray) -> tuple[dict[int, float], dict[int, dict[int, float]]]:
    """RD_k and the two weighted cumulative distributions, over MRS_THRESHOLDS. [§8]

    Takes no threshold list: MRS_THRESHOLDS is prespecified at Stage 1 for this stage, and a keyword
    would make a prespecified choice look like an option [Stage 7 §15, Stage 6 §6.4].

    Returns BOTH the differences and the two distributions they came from, from one call, for
    Stage 7 §3's reason: `RD_k` alone cannot be checked against anything, and a reader given only the
    difference cannot see that P(Y <= k) reached 1 three thresholds early. §9.3's table prints all
    three columns because of this signature, not despite it.

    **Taken from the weighted DISTRIBUTIONS and never from six threshold models** [§8]. The cumulative
    probabilities are then non-decreasing in `k` by construction — a cumulative sum of non-negative
    weights over a denominator that does not depend on `k` — rather than by luck, which is [§8]'s own
    reason: augmenting each threshold separately would replace one primary number with six and can
    produce crossing cumulative probabilities.

    Each threshold's indicator carries the response's own missingness — a comparison against a missing
    value is False in numpy, which would silently count an absent outcome as a non-event — and there
    is no missingness to carry on this path, because §4.1's mask removed it. The mask is the guard;
    this is the statement that it is.

    `share[_TREATED] - share[_COMPARATOR]` and not `share[1] - share[0]`: the arm codes are named once
    at module scope from TREATMENT_LABELS and read here, because this is the one line in the stage
    whose SIGN is a function of which literal goes first, and §7 is four hundred words about that
    sign. Two bare integers in a subtraction are the one place a reader checking the orientation has
    nothing to check against.
    """
    rd: dict[int, float] = {}
    cumulative: dict[int, dict[int, float]] = {}
    for k in C.MRS_THRESHOLDS:
        share = weighted_proportion((y <= float(k)).astype(float), a, w)
        cumulative[k] = share
        rd[k] = share[_TREATED] - share[_COMPARATOR]
    return rd, cumulative


# --- the preconditions [§4.4] ------------------------------------------------------------------------

def _assert_primary_inputs(df: pd.DataFrame, ps: propensity.Propensity) -> None:
    """G1-G5. Phase 1 is `can this be read`; phase 2 is `is the data estimable`. §4.4a.

    Stage 7 §4.5a is the precedent and the argument is unchanged: G4 and G5 and everything in
    `primary` read through the frame, the columns and the mask that G1, G2 and G3 are about, and
    pandas does not return a wrong answer for a broken mask — it raises, from inside the collection,
    before the SchemaError is ever assembled.

    This stage adds a fourth failure of that shape, and it is the reason G3 names the OUTCOME column
    and not only the exposure: phase 2's first line builds the mask from `df[C.PRIMARY_OUTCOME]`, so
    an absent outcome column raises `KeyError('mrs_90d')` there — before G4 and before any message is
    returned. Naming a column in an absent-column check is necessary and is **not** sufficient unless
    that check raises before the column is read.
    """
    # PHASE 1 — can this frame, these columns and this mask be READ?
    bad: list[str] = []

    misaligned = [name for name, s in (("e", ps.e), ("w", ps.w), ("in_model", ps.in_model))
                  if not s.index.equals(df.index)]
    if misaligned:
        bad.append(
            f"G1  the Propensity is not aligned to this frame: {', '.join(misaligned)} "
            f"carr{'ies' if len(misaligned) == 1 else 'y'} a different index. "
            f"{len(ps.in_model)} mask row(s) against {len(df)} frame row(s), sharing "
            f"{len(ps.in_model.index.intersection(df.index))} index label(s). e, w and in_model are "
            "Series on the COHORT's index [Stage 6 §0.2]. This function's caller reads them by LABEL "
            "throughout, so the failure is not a positional read: it is that `w` arrives ordered by "
            "the Propensity's index while the response and the arm vector come from the frame's, so "
            "each record is weighted by another record's score. Measured on the fixture cohort — a "
            "reversed index gives a common odds ratio of 58.01 against the correct 5.76, with "
            "nothing raising. A subset index is worse still: the missing labels resolve to False and "
            "leave the denominator in silence. This is Stage 7's B1 one stage on, and here the wrong "
            "number is the estimate itself.")

    if ps.in_model.dtype != bool or ps.in_model.isna().any():
        bad.append(
            f"G2  in_model is {ps.in_model.dtype} and carries "
            f"{int(ps.in_model.isna().sum())} missing value(s). It is boolean and TOTAL by "
            "construction [Stage 6 §4.4]. A missing entry does NOT raise downstream and does not "
            "propagate: `object & bool` COERCES it to False inside the `&`, so the estimation mask "
            "comes back a clean boolean Series of the right length with that record silently outside "
            "the [§11] denominator. Measured for None, pd.NA and np.nan alike — 4 records estimated "
            "where 5 are in_model, nothing raised.")

    absent = [c for c in (C.PRIMARY_OUTCOME, C.TREATMENT, "case_id") if c not in df.columns]
    if absent:
        bad.append(
            f"G3  {', '.join(absent)}: not a column of the frame. {C.PRIMARY_OUTCOME} is the [§5] "
            f"primary outcome, {C.TREATMENT} is the exposure, and case_id is what entry 1 names its "
            "removed records by; G4, the design and every arm mask read the first two, and "
            "`_record_outcome_completeness` dereferences the third. Stage 7's B3 had to move into "
            "phase 1 for exactly this reason [Stage 7 §4.5a, §21.1] and this check is written there "
            "from the start. case_id is in this list because without it `primary` raises a bare "
            "KeyError('case_id') from inside entry 1 — measured — which is neither of the two types "
            "its docstring declares, and `model.py:164-170` is the landed precedent for exactly this "
            "defect one module over.")

    if bad:
        raise C.SchemaError(
            "\n  " + "\n  ".join(bad)
            + f"\n  {len(bad)} primary-estimate assertion(s) failed over {len(df)} records. The "
            "frame, a column or the mask cannot be read, so G4 and G5 were not run: they read "
            "through all three (§4.4a).")

    # PHASE 2 — readable. These two describe the DATA, over the population §4.1 defines.
    bad = []
    in_estimate = ps.in_model & df[C.PRIMARY_OUTCOME].notna()

    present = [code for code in C.TREATMENT_LABELS
               if int((df.loc[in_estimate, C.TREATMENT] == code).sum())]
    if len(present) < 2:
        bad.append(
            "G4  " + ", ".join(
                f"{C.TREATMENT_LABELS[c]}: "
                f"{int((df.loc[in_estimate, C.TREATMENT] == c).sum())}"
                for c in C.TREATMENT_LABELS)
            + ". A treatment contrast needs both arms. With one, `design` drops the sole predictor "
            "as a constant column and `polr` fits an intercept-only ordinal model and RETURNS a beta "
            "of length zero — measured, and nothing raises on either side (§3.1, §4.5). Stage 6's F2 "
            "and Stage 7's B5 check the same property over their own populations; an arm missing "
            "HERE was lost by in_model or by the outcome mask.")

    nan_w = int((in_estimate & ps.w.isna()).sum())
    if nan_w:
        bad.append(
            f"G5  {nan_w} record(s) are in the estimation population and carry no weight. `w` is "
            "nan off in_model by construction [Stage 6 §9], so this is unreachable while the "
            "estimation mask is built from in_model — and it is checked because a nan weight does "
            "not propagate here, it DELETES: `nan > 0.0` is False, so the positive-weight rule drops "
            "the whole outcome category that record sits in. Measured: one nan weight took a "
            "four-category fit to three and moved beta by 1.48, with nothing raising (§3.1). Never "
            "fill it [Stage 6 §9]; the frame is wrong, not the value.")

    if bad:
        raise C.SchemaError(
            "\n  " + "\n  ".join(bad)
            + f"\n  {len(bad)} primary-estimate assertion(s) failed over {len(df)} records, "
            f"{int(in_estimate.sum())} of them estimable.")


def _assert_exposure_survived(X: pd.DataFrame, dropped: tuple[str, ...]) -> None:
    """The design handed to `polr` carries the exposure, and exactly it. §4.5.

    Not folded into G4. G4 is about the DATA — both arms are present — and this is about the MATRIX,
    and the two can come apart: `design`'s constant-column rule is keyed on `nunique(dropna=False)`
    over the built column, so a future rule change, a covariate list that grew, or an exposure column
    of one distinct value for a reason G4 does not model all arrive here. A design of width 0 is a
    model with no treatment term whose fit CONVERGES (§3.1, measured), so this is the last point at
    which the estimand is still checkable.

    **It raises `model.FitError` and not `SchemaError`**, and that placement is deliberate: [§10]
    catches `FitError` to drop and count a replicate and *"may catch nothing else: a SchemaError here
    is a bug in the resampler, not a sparse replicate"* (`propensity.py:468-469`). A replicate whose
    stratified resample happened to draw one arm at every centre is a sparse replicate, not a bug —
    so it must be droppable, which means it must raise the type Stage 10 catches.
    """
    if tuple(X.columns) != (C.TREATMENT,) or dropped:
        raise model.FitError(
            f"G6  the [§8] design is {tuple(X.columns)} with {dropped or 'nothing'} dropped, and "
            f"[§8] prescribes treatment as the SOLE predictor: exactly ({C.TREATMENT!r},) and "
            "nothing dropped. A width-0 design fits an intercept-only ordinal model, converges, and "
            "returns a beta of length zero — the estimand vanishes and neither `design` nor `polr` "
            "raises (§4.5, measured). `design` drops constant columns by declaration "
            "[Stage 6 §4.3], which is right for a covariate and wrong for the exposure.")


def _assert_reportable(fit: model.PolrFit) -> None:
    """G7 — the fitted treatment coefficient is inside the declared bound. §6.

    Stage 6's F5 is the precedent, structurally and in placement: F5 asserts the fitted probabilities
    are strictly interior and lives in `propensity.py` rather than in `model.py`, because "a boundary
    value is a degenerate fit and not a confident one" is a statement about the specification and not
    about the arithmetic. G7 is the same statement for this estimator, and it is on the CALLER's side
    for the same reason.

    **It is a bound on `beta` and NOT on `alpha`.** A cutpoint's magnitude legitimately grows when a
    category is rare — the separated construction's alpha is -18.2, and there is nothing wrong with a
    rare-category cutpoint of that size on a frame that is not separated — while `beta` is the
    reported effect and `exp(beta)` is the number a manuscript prints. Measured: `alpha` for a rare
    category grows like `log n`, so one rare record in 61 gives `max|alpha|` 4.16 and in 60001 gives
    11.00; reaching 14.0 needs order 1e7 records, and no such frame is constructible.
    """
    worst = float(np.max(np.abs(fit.beta)))
    if worst >= C.POLR_MAX_ABS_BETA:
        raise model.FitError(
            f"G7  |coefficient| reached {worst:.4g} against a bound of {C.POLR_MAX_ABS_BETA:g}, so "
            f"the common odds ratio is {float(np.exp(worst)):.4g}. This is separation, and it does "
            "NOT present as non-convergence: measured, a perfectly separated frame converges in 17 "
            "iterations on the score criterion with every safeguard counter at zero and every fitted "
            "probability finite [§6.1]. Over 4800 fits at twelve true effect sizes, no "
            "non-degenerate fit exceeded 11.04 at any cutpoint count and no degenerate one came "
            "below 18.98; nothing at all lies between them, and every bound from 12 to 20 drops the "
            "same replicates [§6.3]. The likelihood flattens as the coefficient grows and the loop "
            "returns normally. [§10] drops and counts a replicate whose fit fails; without this "
            "raise there is no failure to count, and the percentile interval is a quantile of a "
            f"distribution with a tail at exp({worst:.1f}) [§6.2]. There is no second estimator "
            "[§8, invariant 5].")


# --- the orientation [§7] ----------------------------------------------------------------------------

def _orientation() -> str:
    """The [§8] orientation sentence, computed from the outcome registry. §7.1.

    Not a literal. If an ordinal outcome for which higher is better were ever registered as primary, a
    hardcoded sentence would print `> 1 favours bridging` beside a coefficient meaning the opposite,
    and nothing anywhere would disagree with it. The assertion is the point; the returned string is a
    convenience.
    """
    outcome = C.OUTCOMES[C.PRIMARY_OUTCOME]
    if outcome.kind != "ordinal" or outcome.higher_is_better:
        raise C.SchemaError(
            f"the [§5] primary outcome {C.PRIMARY_OUTCOME!r} is kind={outcome.kind!r} with "
            f"higher_is_better={outcome.higher_is_better}. §5.1's parametrisation gives "
            "`exp(beta) > 1 favours bridging` only for an ordinal outcome on which LOWER is better; "
            "for one on which higher is better the same coefficient favours the comparator, and the "
            "orientation statement [§8] requires would be exactly inverted with every number in the "
            "table still finite and still plausible.")
    return (f"exp(beta) is the common odds ratio and is oriented so that > 1 favours "
            f"{C.TREATMENT_LABELS[_TREATED]}: under logit P(Y <= k) = alpha_k + beta*A a positive "
            f"beta raises P(Y <= k) in the treated arm at every k, and lower {outcome.label} is "
            "better.")


# --- the three `model` audit entries [§9] --------------------------------------------------------------
#
# No new kind, and that is `data.py`'s own declaration honoured: `data.py:143-152` declares `model` for
# "what was fitted to the population `cohort` describes" and adds that "Stages 7-13 all render under
# `model` or under an existing kind". A weighted proportional-odds fit is a fitted model, its
# coefficient is the estimate, and the cumulative risk differences are that estimate on an absolute
# scale. So KINDS stays nine and none of the literal pins in the six existing test modules moves.
#
# Four rules about the tables, all inherited:
#   * every DECLARED thing is rendered whether or not the data fills it — a category that vanishes
#     from a table is a category nobody knows was considered (Stage 6 §7.5, Stage 7 §4.2)
#   * every cell is computed from the object it describes, never by subtracting another row: a table
#     that reconciles with itself by construction has stopped being evidence (Stage 5 §7.3)
#   * `_fmt` is the only float formatter and `missing` the only rendering of an absent value
#   * NO CELL OF ANY TABLE CONTAINS A `|`. `data._md_table` does no escaping and sizes its separator
#     from `len(rows[0])`, so one pipe in a header cell gives a header row with more markdown cells
#     than its body — and byte-identity across hash seeds does not catch it, because a table
#     identically broken under both seeds is still identical (Stage 7 §21.1).


def _completeness_table(df: pd.DataFrame, ps: propensity.Propensity,
                        in_estimate: pd.Series) -> tuple[tuple[str, ...], ...]:
    """The [§11] denominator chain, one row per population step. §4.1.

    Four rows, and THEY ARE NOT A CUMULATIVE CHAIN — row 3 is a count over the whole cohort, as its
    own third cell says, so adjacent differences are not restriction costs. On the workbook the four
    read 93, 92, 93, 92 and the differences are -1, +1, -1: a "restriction" that adds a record.
    `cohort_flow` (Stage 5 §7.3) IS a nested chain; this is a chain of two with a parallel denominator
    printed beside it, because [§11] asks for the outcome's own completeness as well as the estimate's.

    Row 3 is kept and not made nested, because the quantity a reader wants when the last row is short
    is "how many outcomes are present at all", and deriving it by subtraction from a nested row is the
    reconciles-with-itself table Stage 5 §7.3 forbids. The last row is what the estimate below it is
    computed on, and it is the number [§11] asks to be reported for every estimate.
    """
    present = df[C.PRIMARY_OUTCOME].notna()
    header = ("population", "n", "what it is")
    return (
        header,
        ("[§3] primary cohort", str(len(df)), "both restrictions applied [Stage 5 §11]"),
        ("in_model", str(int(ps.in_model.sum())),
         "covariate-complete; the ATO population [Stage 6 §4.4]"),
        (f"{C.PRIMARY_OUTCOME} present (whole cohort)", str(int(present.sum())),
         "NOT a step in the chain — the outcome's own completeness, before the mask below"),
        ("in_estimate", str(int(in_estimate.sum())),
         "in_model AND outcome present — this estimate's [§11] denominator"),
    )


def _record_outcome_completeness(df: pd.DataFrame, ps: propensity.Propensity,
                                 in_estimate: pd.Series, audit: Audit) -> None:
    """Entry 1, and it must name as many patients as it says it removed.

    `propensity._record_exclusion` (propensity.py:400-432) line for line, both mechanisms included,
    because both are what make the comparison able to fail at all: the identifiers go through a `set`,
    so a duplicated case_id collapses to one name, and a missing case_id is dropped by `notna()`
    rather than stringified. Without either, `len(case_ids)` equals `n` by construction and the check
    is dead code.

    The raise comes BEFORE `audit.record`, again as Stages 5 and 6 do it: recording first and reading
    the entry back leaves an unreconcilable entry in an `Audit` that has no removal path.

    On v7 this entry's `n` is 0 and it names nobody (§4.1, measured), and it is recorded anyway. An
    entry that appears only when it has something to say is an entry whose absence a reader has to
    interpret, and [§11] asks for the denominator of every estimate rather than for the denominators
    that happen to differ.
    """
    removed = ps.in_model & ~in_estimate
    excluded = df.loc[removed & df["case_id"].notna(), "case_id"]
    case_ids = tuple(sorted(set(excluded)))
    n_removed = int(removed.sum())
    if len(case_ids) != n_removed:
        raise C.SchemaError(
            f"{_STEP_COMPLETENESS} removes {n_removed} record(s) from the ATO population and can "
            f"name {len(case_ids)}. An estimate whose denominator the log cannot reconstruct is an "
            "estimate nobody can check [§11].")
    audit.record("model", _STEP_COMPLETENESS, n_removed,
                 f"[§11] complete-case per estimate. {int(in_estimate.sum())} of "
                 f"{int(ps.in_model.sum())} weighted record(s) carry the [§5] primary outcome and "
                 f"are estimated; {n_removed} do not, and every one of them is named above. Such a "
                 "record keeps its weight and loses its estimate — it is in the ATO population "
                 "[Stage 6 §4.4] and not in this estimate's denominator, which is why this is a "
                 "THIRD population and not a restatement of that one [§4.1]. The primary estimate "
                 "and the cumulative risk differences below are computed on the same mask, bound "
                 "once, because [§8] presents them as one contrast on two scales.",
                 case_ids=case_ids, table=_completeness_table(df, ps, in_estimate))


def _fit_table(fit: model.PolrFit) -> tuple[tuple[str, ...], ...]:
    """One row per DECLARED mRS level, then the coefficient and the odds ratio. §9.2.

    Ranges over MRS_LEVELS and never over `fit.categories`, for `propensity._design_table`'s reason
    (propensity.py:242-261): a level absent from the fitted set is a fact about the population, and a
    row that vanishes is a fact nobody sees. The `status` column is where the collapse of §5.3 becomes
    visible in the log rather than inferable from a parameter count.
    """
    fitted = {level: value for level, value in zip(fit.categories[:-1], fit.alpha)}
    header = ("parameter", "mRS level", "estimate", "status")
    rows: list[tuple[str, ...]] = []
    for level in C.MRS_LEVELS:
        # THE NO-CUTPOINT TEST KEYS ON `fit.categories[-1]` AND NEVER ON `MRS_LEVELS[-1]`. There are
        # two different reasons a declared level has no cutpoint — it is absent from the population,
        # or it is the TOP OF THE FITTED SCALE where P(Y <= it) = 1 — and only the second is a
        # property of the fit. The two tests coincide on every frame this suite has (the workbook
        # occupies all seven levels and `ordinal_cohort()`'s top category is 6 by coincidence) and
        # come apart in any [§10] replicate that draws no death, where the declared form renders
        # `missing` under status `fitted` with nothing raising. Do not simplify it back.
        if level not in fit.categories:
            status = ("absent from this population — no cutpoint [§5.3]"
                      + (", and it is the highest declared level"
                         if level == C.MRS_LEVELS[-1] else ""))
        elif level == fit.categories[-1]:
            status = "highest FITTED level — no cutpoint, P(Y <= it) = 1"
        else:
            status = "fitted"
        rows.append((f"alpha at mRS <= {level}", str(level), _fmt(fitted.get(level)), status))
    for name, value in zip(fit.columns, fit.beta):
        rows.append((name, "—", _fmt(value), "treatment coefficient [§8]"))
        rows.append((f"exp({name})", "—", _fmt(float(np.exp(value))),
                     "common odds ratio — see the orientation above"))
    return (header, *rows)


def _fit_detail(fit: model.PolrFit, n: int) -> str:
    absent = tuple(level for level in C.MRS_LEVELS if level not in fit.categories)
    return (
        f"[§8] weighted proportional-odds model with treatment as the sole predictor, over {n} "
        f"record(s). Converged in {fit.iterations} iteration(s) on the {fit.converged_on} "
        f"criterion; {fit.rescales} step(s) were shortened by the trust region and {fit.halvings} "
        f"halving(s) taken. {len(fit.categories)} outcome category/ies carried positive weight, so "
        f"{len(fit.alpha)} cutpoint(s) were estimated"
        + (f"; {len(absent)} declared level(s) are absent from this population and are marked as "
           f"such below: {', '.join(str(level) for level in absent)}. " if absent else ". ")
        + _orientation()
        + " NO STANDARD ERROR AND NO INTERVAL IS PRODUCED HERE, and that is a property of the "
        "estimator rather than a division of labour: the observation weights are a tilting function "
        "of an ESTIMATED propensity score and not frequencies, so the inverse observed information "
        "omits the variability of estimating it [§7]. [§10]'s percentile bootstrap over N_BOOT "
        "replicates is the prespecified interval and the source of the primary p-value, which is a "
        "test of this coefficient and NOT a risk-difference-scale test [§10]. The estimate is "
        "marginal in the overlap population [§8] and is conditional on the [§7] propensity "
        "specification, which defines that population [§7].")


def _rd_table(rd: dict[int, float],
              cumulative: dict[int, dict[int, float]]) -> tuple[tuple[str, ...], ...]:
    """One row per DECLARED threshold: the two weighted cumulative probabilities and RD_k. §8

    No cell and no header contains a `|` — the columns are named `P(Y <= k | …)` in mathematics and
    `P(Y <= k), <arm>` here, because `data._md_table` does no escaping and a pipe in a header cell
    gives a header row with more markdown cells than its body [Stage 7 §21.1].

    The RD column is read from `rd` and not computed from the two printed cells: a table that
    reconciles with itself by construction has stopped being evidence [Stage 5 §7.3].
    """
    header = ("threshold", *(f"P(Y <= k), {label}" for label in C.TREATMENT_LABELS.values()),
              "RD_k")
    rows = tuple(
        (f"mRS <= {k}", *(_fmt(cumulative[k][code]) for code in C.TREATMENT_LABELS), _fmt(rd[k]))
        for k in C.MRS_THRESHOLDS)
    return (header, *rows)


def _rd_detail(totals: dict[int, float], n: int) -> str:
    return (
        f"[§8] weighted empirical cumulative risk differences, RD_k = P(mRS <= k | "
        f"{C.TREATMENT_LABELS[_TREATED]}) - P(mRS <= k | {C.TREATMENT_LABELS[_COMPARATOR]}), over "
        f"{len(C.MRS_THRESHOLDS)} declared threshold(s) and {n} record(s) — the same records and the "
        "same weights as the coefficient above [§4.2]. Taken from the weighted DISTRIBUTIONS and "
        "not from threshold-specific models, so the cumulative probabilities are non-decreasing in k "
        "by construction rather than by luck [§8]. Ties between adjacent thresholds are expected "
        "wherever a declared level is unoccupied and are not a defect [§8.3]. Per-arm weight totals, "
        "which are these proportions' denominators [§11]: "
        + ", ".join(f"{C.TREATMENT_LABELS[code]} {_fmt(total)}" for code, total in totals.items())
        + f". mRS <= {C.MRS_LEVELS[-1]} is not reported because P(Y <= {C.MRS_LEVELS[-1]}) is 1 in "
        f"both arms by definition, so its risk difference is structurally zero; mRS <= "
        f"{C.MRS_THRESHOLDS[-1]} is the survival contrast and is not [§8.3]. THESE CARRY CONFIDENCE "
        "INTERVALS BUT NO p-VALUES [§8]: they are the absolute-scale presentation of one contrast, "
        "and six threshold-wise tests of it would only invite selection of the most favourable cut. "
        "The primary outcome has exactly one test [§8, §10].")


# --- the stage [§8] ------------------------------------------------------------------------------------

def primary(df: pd.DataFrame, ps: propensity.Propensity, audit: Audit) -> Primary:
    """The [§8] primary outcome estimate over the [§3] cohort and the [§7] weights.

    Takes no outcome name: [§8] rests on there being ONE primary quantity with one test, and
    PRIMARY_OUTCOME is computed from the [§5] registry. Takes no covariate list: treatment is the sole
    predictor by specification. Takes no threshold list: MRS_THRESHOLDS is Stage 1's. A keyword for
    any of the three would make a prespecified choice look like an option [Stage 6 §6.4, Stage 7 §15].

    Returns a Primary; adds no column to `df` and edits nothing (§0.2). Appends three `model` entries.

    Raises SchemaError on G1-G5 and on D1-D4, and model.FitError on G6, G7, O1-O6 or a fit that does
    not converge. [§10] catches the second to drop and count a replicate, and may catch nothing else:
    a SchemaError here is a bug in the resampler, not a sparse replicate (propensity.py:468-469).

    Six things about the order below, each of which is a failure if moved:

      * `_assert_primary_inputs` comes FIRST, in two phases, before the mask is built for anything
        but its own phase-2 checks (§4.4a).
      * `in_estimate` is bound ONCE and every subsequent line reads it — the design, the response,
        the weights, the arm vector and the cumulative distributions. Two separate reads would be two
        chances to compute the coefficient and the risk differences over different populations, which
        [§8] presents as one contrast (§4.2).
      * `_record_outcome_completeness` comes BEFORE the design, so the log names the [§11] denominator
        even if `design` then raises on D2 — the exclusion is a fact about the frame, established
        before any modelling choice, and a log that names it is useful precisely when what follows
        fails (propensity.py:440-443).
      * `_assert_exposure_survived` sits BETWEEN the design and the fit, because a width-0 design FITS
        and the estimand's absence stops being visible one line later (§4.5).
      * `_assert_reportable` sits BETWEEN the fit and everything that reads it, because a separated
        fit's `exp(beta)` is a finite float that every downstream table will accept (§6).
      * `cumulative_rd` runs AFTER the fit rather than before, so that a frame which cannot support
        the primary quantity does not produce the absolute-scale presentation of a contrast that has
        no estimate. [§8] makes RD_k a presentation of `beta`, not an alternative to it.

    `fit.beta[0]` and not `fit.beta.sum()`, and G6 is what makes the index safe. A width-0 design
    gives `beta` of length zero, so `beta[0]` would be an `IndexError` and `beta.sum()` would be
    `0.0` — an odds ratio of exactly 1.0, which reads as "no effect" and is the estimand having
    vanished (§4.5, measured).

    `cut_levels` is `fit.categories[:-1]` and that is the identity `alpha` needs: `alpha[j]` is the
    cutpoint between `categories[j]` and `categories[j+1]`, so the level it cuts at or below is
    `categories[j]`. On a collapsed fit those are NOT `MRS_THRESHOLDS` — a cohort missing mRS 2 has a
    cutpoint at `mRS <= 1` whose next category is 3, and reporting it against `MRS_THRESHOLDS[1]`
    would label a cutpoint with a threshold it is not.
    """
    _assert_primary_inputs(df, ps)                                     # G1-G3 | G4-G5

    in_estimate = ps.in_model & df[C.PRIMARY_OUTCOME].notna()          # §4.1 — bound ONCE
    _record_outcome_completeness(df, ps, in_estimate, audit)           # entry 1

    sub = df.loc[in_estimate]
    X, dropped = model.design(sub, (C.TREATMENT,))                     # D1-D4 inherited
    _assert_exposure_survived(X, dropped)                              # G6, §4.5

    y = sub[C.PRIMARY_OUTCOME].to_numpy(dtype=float)
    a = sub[C.TREATMENT].to_numpy(dtype=float)
    w = ps.w.loc[in_estimate].to_numpy(dtype=float)

    fit = model.polr(X, y, w)                                          # raises FitError; §5
    _assert_reportable(fit)                                            # G7, §6

    n = int(in_estimate.sum())
    audit.record("model", _STEP_FIT, n, _fit_detail(fit, n), table=_fit_table(fit))

    rd, cumulative = cumulative_rd(y, a, w)                            # §8
    totals = {code: float(w[a == float(code)].sum()) for code in C.TREATMENT_LABELS}
    audit.record("model", _STEP_RD, len(C.MRS_THRESHOLDS), _rd_detail(totals, n),
                 table=_rd_table(rd, cumulative))

    beta = float(fit.beta[0])
    return Primary(
        beta=beta,
        odds_ratio=float(np.exp(beta)),
        alpha=tuple(float(v) for v in fit.alpha),
        cut_levels=tuple(int(level) for level in fit.categories[:-1]),
        rd=rd,
        cumulative=cumulative,
        in_estimate=in_estimate,
        fit=fit)
