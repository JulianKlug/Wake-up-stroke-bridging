"""Stages 8b and 9 — the [§8] outcome estimates, primary and secondary.

**Two stages, one module, and Stage 8 §0.1 is why**: *"Stage 9's binary estimators extend that module
so the weighted per-arm proportion has one home."* [§11] requires the denominator of every estimate
and ``weighted_proportion`` **is** that denominator — the [§8] weighted risk difference is it
differenced and the weighted marginal odds ratio is those two proportions. A second implementation in
a second module would be a second definition of it. The file therefore has two halves, banner-divided:
everything down to ``primary`` is Stage 8's ordinal quantity, and everything from ``BinaryEstimate``
onward is Stage 9's seven binary ones. The second half reads nothing the first half produces —
``Primary`` is not an input to any binary estimate, and §15.11 asserts that by scan for Stage 7
§14.11's reason: not because reading it would fail, but because it would succeed.

Stage 9's own specification is ``specs/stage9_secondary_binary_estimators.md``; the banner at the head
of the second half states what it adds and what it deliberately leaves alone.

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

from dataclasses import dataclass, field
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

# The three audit step BASES. Named as constants rather than written at the three `audit.record` call
# sites because the acceptance suite asserts the entries appear in this order at positions captured
# before the call, and a test comparing against a string literal it also writes is a test of nothing.
#
# **Every one reaches `audit.record` through `ps.spec.step`** [Stage 11 §4.4, §10]. `primary` is
# called TWICE over one cohort from Stage 11 onward — once under [§7] and once under [§13]'s
# full-covariate specification — and `Audit.entry` is first-match (data.py:273-275), so two
# unsuffixed `primary_fit` entries would make every programmatic read return the [§7] estimate's
# while the rendered log looked complete. The primary specification's suffix is the empty string, so
# the three names below are byte-identical on a run without the arm and no landed assertion moves.
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


# --- the [§11] estimation population [§4.1] ---------------------------------------------------------

def estimation_population(df: pd.DataFrame, ps: propensity.Propensity) -> pd.Series:
    """The [§11] population the primary estimate is computed on: weighted AND outcome present.

    Boolean and TOTAL on `df`'s index, exactly as `Propensity.in_model` is and for Stage 6 §9's
    reason: the deliberateness of an absence lives in a mask, never in a value.

    **Public, and it is a DEFINITION rather than a step** — `weighted_proportion`'s own argument
    applied to the quantity that argument names: *"a second implementation in a later stage would be
    a second definition of the denominator, and [§11] requires the denominator."* `primary` binds it
    once and every line of that function reads the result; `_assert_primary_inputs`' phase 2 reads
    the same definition rather than spelling it again; and Stage 11's subgroup replicate body needs
    this mask on a DRAWN frame without paying for a primary fit it does not read, which is why the
    definition is a name and not three words inside one function (Stage 11 §8.5).

    It is a THIRD population and not a restatement of `in_model`: a record whose covariates are
    complete and whose outcome is missing keeps its weight and loses its estimate (§4.1). Nothing
    here is applied — the mask is returned for a caller to report and to log, which is
    `model.complete_cases`' rule one stage on.
    """
    return ps.in_model & df[C.PRIMARY_OUTCOME].notna()


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
    in_estimate = estimation_population(df, ps)             # the ONE definition — Stage 11 §8.5

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
            f"{ps.spec.step(_STEP_COMPLETENESS)} removes {n_removed} record(s) from the ATO "
            f"population and can "
            f"name {len(case_ids)}. An estimate whose denominator the log cannot reconstruct is an "
            "estimate nobody can check [§11].")
    audit.record("model", ps.spec.step(_STEP_COMPLETENESS), n_removed,
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

    in_estimate = estimation_population(df, ps)                        # §4.1 — bound ONCE
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
    audit.record("model", ps.spec.step(_STEP_FIT), n, _fit_detail(fit, n), table=_fit_table(fit))

    rd, cumulative = cumulative_rd(y, a, w)                            # §8
    totals = {code: float(w[a == float(code)].sum()) for code in C.TREATMENT_LABELS}
    audit.record("model", ps.spec.step(_STEP_RD), len(C.MRS_THRESHOLDS), _rd_detail(totals, n),
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


# ======================================================================================================
# Stage 9 — the [§8] secondary binary estimators
# ======================================================================================================
#
# Seven binary outcomes — three secondary and four safety — each with a weighted risk difference, a
# weighted marginal odds ratio, and, where the minority cell supports a nuisance model, a
# MODEL-ASSISTED augmented risk difference tilted by h = e(1 - e). Each on its own [§11] denominator.
# Point estimates only: every interval and every p-value is [§10]'s.
#
# **Why this is here and not in a module of its own** [Stage 9 §0.1]. `weighted_proportion` above IS
# the [§11] denominator, and the weighted risk difference is that function differenced while the
# weighted marginal odds ratio is those two proportions. A second implementation in a second module
# would be a second definition of the denominator — the move Stage 6 §3 declined for `ess` and
# Stage 7 §0.1 for `smd`. The split that IS made is the one Stages 6 and 8 both made: the numerics go
# in `model.py` — Stage 9 adds exactly one function there, `predict` — and the [§8] specification
# stays here, where the exposure is already named.
#
# **`Primary` is not read, and `Balance` is neither read nor imported** [Stage 9 §0.2]. The binary
# estimators are their own estimates on their own denominators, and the primary ordinal quantity is
# not an input to any of them. `model.polr` is neither called nor edited: no binary estimate goes near
# the ordinal fitter, and that is recorded so nobody wires it in looking for a use.
#
# **"MODEL-ASSISTED", NEVER "DOUBLY ROBUST"** [§8, roadmap]. The ATO estimand is itself indexed by the
# true propensity score — the target population is defined by h(X) = e(X){1 - e(X)} and corresponds to
# no observable subset of the data — so a misspecified e(X) makes the ESTIMAND wrong even when m_a(X)
# is exactly right. The estimator stays consistent for a weighted average treatment effect, for the
# wrong weighting. The word licenses a conclusion this estimator does not support, and the failure
# mode is a word in a docstring that no behavioural test can catch, so the suite scans for it instead.


# --- what the stage returns [§3.1] ---------------------------------------------------------------------

@dataclass(frozen=True)
class BinaryEstimate:
    """One [§8] binary outcome's three estimates, on its own [§11] population.

    **Four fields move as a block and `reduced` is NOT one of them.** `augmented`, `covariates`,
    `fit` and `dropped` are `None`/`None`/`None`/`()` exactly when `augmented_path` is
    `"unaugmented"`, and all four are populated otherwise: an outcome is either augmented or it is
    not, and §9.2 is the only thing that decides which. A reader must not have to infer the path from
    which fields are absent, so `augmented_path` states it as one of three declared strings and
    §15.9 asserts the block.

    `reduced` is **orthogonal to the path**, and an earlier draft of this docstring listed it in the
    block, which is false on four of the seven outcomes: `reduced` is `True` exactly when the key is
    in `OUTCOME_MODEL_OVERRIDES`, so on v7 it is `True` for one outcome alone and `False` for the four
    `full`-path outcomes *whose other three fields are populated*. A block-move test written from the
    earlier sentence fails on four of the seven. §15.9 asserts the block and `reduced` separately.

    **There is no verdict and no `favours_bridging()`**, for Stage 8 §3's reason: [§8] specifies no
    threshold on a risk difference, the null is tested at Stage 10 from the bootstrap distribution,
    and a method here returning `rd > 0.0` would be a significance test with no interval behind it.

    `or_corrected` means A WEIGHTED PROPORTION REACHED 0 OR 1, which is very nearly but not exactly
    "a cell was empty". The guard tests the float, so a cell carrying a weight small enough that the
    proportion rounds to exactly 1.0 sets the flag while the cell is not empty -- measured, a
    non-event weight of 1e-18 gives p1 == 1.0 and `or_corrected` True with an uncorrected odds ratio
    of +inf. One-directional: an empty cell always gives exactly 0.0 or 1.0, so there are no false
    negatives, and the correction is right in both cases because +inf is what it exists to avoid.
    §10.3's column therefore reads "corrected" and not "empty cell", and §15.5 asserts the
    distinction rather than leaving a reader to infer a cell count from a flag.

    `or_corrected` is a FLAG AND NOT A MODE. [§8] gives no degenerate-cell rule and the roadmap says
    only "kept finite"; §6.3 prescribes the correction and requires that it be printed wherever the
    estimate is. The two uncorrected proportions travel with it in `proportion` so that a reader can
    see which cell was empty — the correction is never the only record that it fired.
    """

    outcome: str                             # the C.OUTCOMES key
    family: str                              # "secondary" | "safety" — Stage 11 groups on this
    minority: int                            # min(events, non-events) on in_estimate — §9.1
    rd: float                                # the weighted risk difference — §5.2
    odds_ratio: float                        # the weighted marginal odds ratio — §6
    or_corrected: bool                       # §6.3's correction fired
    proportion: dict[int, float]             # the two weighted proportions, UNCORRECTED — §5.1
    augmented_path: str                      # "full" | "reduced" | "unaugmented" — §9.2
    augmented: float | None                  # tau — §8; None when unaugmented
    covariates: tuple[str, ...] | None       # m_a(X)'s declared list, treatment excluded — §7.1
    reduced: bool                            # the list is an OUTCOME_MODEL_OVERRIDES entry — §7.2
    dropped: tuple[str, ...]                 # design columns dropped as constant — §7.1
    in_estimate: pd.Series                   # boolean, TOTAL — the [§11] population — §4.2
    fit: model.Fit | None                    # the m_a(X) fit; None when unaugmented


@dataclass(frozen=True)
class Secondary:
    """The seven [§8] binary estimates, keyed by outcome in C.BINARY_OUTCOMES order.

    A wrapper rather than a bare dict for one reason: `by_family` is Stage 11's grouping and [§13]
    applies Benjamini-Hochberg WITHIN each family. A consumer that groups by reading
    `C.OUTCOMES[k].family` itself is a second place the family partition is computed, and [§13]'s
    correction is wrong if the two disagree. This is not a verdict method — it partitions, it does
    not judge.

    **`estimates` and `failures` PARTITION `C.BINARY_OUTCOMES` and `failures` is empty at the
    point-estimate contract** [Stage 10 §7.3]. `secondary(..., collect=True)` is Stage 10's call and
    the only one that can fill it: an outcome whose `m_a(X)` cannot be fitted is recorded here by
    name and message rather than raised past the other six, which is [§10]'s "dropped and counted" at
    the granularity [§10] left open. With `collect=False` — every caller before Stage 10 — a
    `FitError` propagates and `failures` is `{}`, so nothing about the point estimate moves.

    The key is RECORDED and not simply absent, because "we tried this outcome and could not fit it"
    and "this outcome was never attempted" are the two things [§10]'s counters have to tell apart
    (Stage 10 §3.1, §15.7).
    """

    estimates: dict[str, BinaryEstimate]
    failures: dict[str, str] = field(default_factory=dict)   # outcome -> the FitError message

    def by_family(self) -> dict[str, tuple[BinaryEstimate, ...]]:
        """The estimates grouped by [§5] family, each tuple in C.BINARY_OUTCOMES order.

        Ranges over the registry and SKIPS an outcome that is not in `estimates`, which under the
        point-estimate contract is none of them. Under `collect=True` a family can come back one
        short, and a `KeyError` there would be this method reporting a Stage 10 replicate's dropped
        outcome as a broken registry [Stage 10 §7.3].
        """
        out: dict[str, list[BinaryEstimate]] = {}
        for key in C.BINARY_OUTCOMES:
            if key in self.estimates:
                out.setdefault(self.estimates[key].family, []).append(self.estimates[key])
        return {fam: tuple(v) for fam, v in out.items()}


# --- the preconditions, in THREE phases [§4.4] ---------------------------------------------------------
#
# Three collections, not two, and the third split is Stage 7's own defect repeating. Stage 7's learnings
# entry states the rule: "the boundary is can-this-be-READ vs is-the-DATA-judgeable — NOT mask-checks vs
# the-rest. A column-presence check is a precondition of every read that follows it, exactly as a mask
# check is." A collection is only a collection if every check in it can RUN: with S1 and S2 in one
# phase, a frame missing an outcome column has S1 record its failure and S2 then raise a bare KeyError
# before the collected SchemaError is ever built — so the caller sees a pandas traceback instead of the
# message written here.
#
# Stage 6 solved it the other way and it is recorded because the two are not equally safe:
# `model._assert_design_inputs` keeps one phase and re-filters inside it (model.py:152-154). That is
# correct BY VIGILANCE — every future check must remember to re-filter. The three-phase split is correct
# BY STRUCTURE: a check placed in phase 2 cannot see a column phase 1 rejected, because phase 1 has
# already raised.
#
# Every one of S1-S8 is unreachable on the workbook. That is recorded rather than treated as a reason to
# skip one — Stage 6 §4.5's position — and [§10]'s frames are the population they are written for, which
# is precisely why the suite is the only thing that can reach them.


def _assert_readable(df: pd.DataFrame, ps: propensity.Propensity) -> None:
    """S1, S5a, S3a, S3b: what must hold BEFORE ANY VALUE IS READ. Raises; never returns a verdict.

    Phase 1 of three (§4.4a). Every check here is a precondition of a READ that a later phase
    performs, so this function raising is what makes the later phases' collection honest.

    `.index.equals(df.index)` and NOT a length check or a set comparison, which is S3a's whole
    content: `.loc[boolean_series]` returns rows in the SERIES' order, while `m1` and `m0` come from
    `model.design(df.loc[mask])` in the FRAME's order. A permuted-but-equal index therefore pairs
    each record's weight with a different record's fitted value -- no raise, no nan, a different
    number. `propensity.fit` builds e and w from df so the orders agree there, and nothing in S4
    constrains a Propensity a caller assembles by hand [Stage 9 §12].

    S3b IS HERE AND NOT IN PHASE 2, because S4's `ps.e[ps.in_model]` is a masked read: an integer
    in_model is interpreted as LABELS and raises KeyError before any collected SchemaError exists
    (§4.4). Measured, not reasoned -- it came out of calling `secondary` end to end.
    """
    failures: list[str] = []
    if ps.in_model.dtype != bool:
        failures.append(
            f"S3b ps.in_model has dtype {ps.in_model.dtype} and not bool. Every read below masks "
            "with it, and pandas reads a non-boolean Series as LABELS -- so this raises KeyError "
            "rather than reporting a mask problem [Stage 9 §4.4a].")
    if len(ps.in_model) != len(df):
        failures.append(
            f"S3b ps.in_model has {len(ps.in_model)} entries against the frame's {len(df)}. The "
            "mask is TOTAL [Stage 6 §9]: absence lives in a False, never in a missing row.")
    absent = [k for k in C.BINARY_OUTCOMES if k not in df.columns]
    if absent:
        failures.append(
            f"S1  {', '.join(absent)}: not a column of the frame. Every C.BINARY_OUTCOMES key is "
            "estimated by this stage, so a registry entry with no column is a frame that cannot "
            "answer what [§8] asks of it.")
    if C.TREATMENT not in df.columns:
        failures.append(
            f"S5a {C.TREATMENT}: not a column of the frame. Every estimate here is a contrast "
            "between the two declared arms [Stage 4 §4.2].")
    for name, series in (("in_model", ps.in_model), ("e", ps.e), ("w", ps.w)):
        if not series.index.equals(df.index):
            failures.append(
                f"S3a ps.{name} is not indexed like the frame. `.loc` returns the SERIES' order, so "
                "a permuted index pairs each record's weight with another record's fitted value and "
                "returns a different number without raising [Stage 9 §4.4b].")
    if failures:
        raise C.SchemaError("\n".join(["outcome.secondary cannot read this frame:", *failures]))


def _assert_secondary_inputs(df: pd.DataFrame, ps: propensity.Propensity) -> None:
    """S2, S4, S5b: whether the DATA is judgeable. Collected into one SchemaError.

    Phase 2 of three. Runs only on what `_assert_readable` passed, which is what lets it read
    df[key], df[C.TREATMENT], ps.e and ps.w without re-filtering -- Stage 6's
    `_assert_design_inputs` re-filters instead (model.py:152-154) and §4.4a is why this one does not.

    S2 is re-asserted on the four DERIVED_DICHOTOMIES, which C.BINARY_COLUMNS does not cover: those
    columns are Stage 3's and their 0/1-ness is a property of the derivation rather than of the
    workbook's contract.
    """
    failures: list[str] = []
    for key in C.BINARY_OUTCOMES:
        offside = df.loc[df[key].notna(), key]
        bad = sorted(set(offside[~offside.isin((0.0, 1.0))].tolist()))
        if bad:
            failures.append(
                f"S2  {key}: carries {bad} outside {{0, 1}} and missing. `y == 1` is False on any "
                "other value, so a third level would be counted as a non-event silently "
                "[Stage 9 §3.3].")
    for name, series in (("e", ps.e), ("w", ps.w)):
        inside, outside = series[ps.in_model], series[~ps.in_model]
        if not np.isfinite(inside.to_numpy(dtype=float)).all():
            failures.append(
                f"S4  ps.{name} is not finite on in_model. Every weighted quantity divides by a sum "
                "of these [Stage 6 §9].")
        if not outside.isna().all():
            failures.append(
                f"S4  ps.{name} carries {int(outside.notna().sum())} value(s) OFF in_model. Stage 6 "
                "§9's rule is that a deliberate absence lives as nan and is never filled; a finite "
                "value off the mask is a value some consumer will range over.")
    arm = df[C.TREATMENT]
    if arm.isna().any():
        failures.append(
            f"S5b {C.TREATMENT}: {int(arm.isna().sum())} record(s) carry no arm. Treatment is "
            "assigned for every cohort member by Stage 4 [§4.2]; a missing arm is a bug upstream.")
    bad_arm = sorted(set(arm.dropna()[~arm.dropna().isin([float(c) for c in C.TREATMENT_LABELS])]))
    if bad_arm:
        failures.append(
            f"S5b {C.TREATMENT}: carries {bad_arm} outside the declared arm codes "
            f"{sorted(C.TREATMENT_LABELS)}.")
    if failures:
        raise C.SchemaError("\n".join(["outcome.secondary cannot estimate on this frame:", *failures]))


def _assert_estimable(key: str, y: np.ndarray, a: np.ndarray, w: np.ndarray) -> None:
    """S6, S7, S8, on ONE outcome's [§11] population. Names the outcome in every message.

    Phase 3 of three, and the only one that runs inside the loop -- so its messages carry `key`,
    because "the outcome is constant" against seven outcomes costs a bisect (§15.1a).

    S7 RAISES rather than tolerating the nan `weighted_proportion` would return, and §4.4 is the
    argument: that function returns nan when an arm has no positive weight, so S7 firing means the
    MASK and the WEIGHTS disagree about the same rows, which is a bug and not a sparse replicate.

    **S8 IS `model.FitError` AND S6 AND S7 ARE `C.SchemaError`, AND THE SPLIT IS THE WHOLE OF
    DECISION 6** (PI, 2026-08-25) [Stage 10 §5.4]. Stage 9 wrote S8 as `SchemaError` and named the
    classification as the open question of §16 item 3; `TODOS.md` set the trigger -- the first
    Stage 10 replicate that hits it, and a non-negligible rate -- and the trigger fired. Measured
    over 1998 stratified replicates: `sich` 16 (0.8%), `tici_2b_3` 3 (0.2%), 19 replicates and 1.0%
    in total. `sich` has five events in ninety-two, so a resample drawing none of them is a sparse
    replicate and not a bug, and Stage 8 §11 forbids [§10] catching a `SchemaError` -- so as landed
    this raised something no correct implementation could catch and no correct resampler prevent.

    Stage 9's own reason for the other choice is what settles it. The objection was to RETURNING a
    plausible number -- a constant outcome does not yield nan, it yields a finite 1.3529411765
    (§6.2) for an outcome nobody observed an event in -- and `FitError` returns nothing: it drops
    the replicate and counts it in a named bucket, which is [§10]'s mechanism for exactly this. S6
    is an empty population and S7 is an arm with no record or no weight; both mean the mask and the
    weights disagree about the same rows, which is a bug, so both stay `SchemaError`. Neither fired
    in 2000 replicates [Stage 10 §5.4, §7.1].

    **S7 and S8 raise on the FIRST failure rather than collecting**, unlike the two frame-level
    phases, and the asymmetry is deliberate: a frame-level phase reports on a frame the caller can
    fix in one edit, where a per-outcome failure means THIS outcome is unestimable and the remaining
    six are a separate question. Collecting across the loop would report seven problems for one bad
    mask.
    """
    if y.size == 0:
        raise C.SchemaError(
            f"S6  {key}: its [§11] population is empty. `in_model & notna()` selected no record, so "
            "there is no denominator and every estimate below would divide by zero.")
    for code in C.TREATMENT_LABELS:
        arm = a == float(code)
        if not arm.any():
            raise C.SchemaError(
                f"S7  {key}: arm {code} ({C.TREATMENT_LABELS[code]}) has no record on its [§11] "
                "population. Every estimate here is a contrast and a contrast needs both arms.")
        if not float(w[arm].sum()) > 0.0:
            raise C.SchemaError(
                f"S7  {key}: arm {code} ({C.TREATMENT_LABELS[code]}) carries zero total weight over "
                f"{int(arm.sum())} record(s). weighted_proportion returns nan rather than raising "
                "here (outcome.py:155-157), so reaching this means the mask and the weights "
                "disagree about the same rows -- a bug, not a sparse replicate [Stage 9 §4.4].")
    if len(set(y.tolist())) < 2:
        raise model.FitError(
            f"S8  {key}: constant at {y[0]!r} on its [§11] population of {y.size}. The risk "
            "difference is 0 and both weighted proportions are degenerate in the SAME direction, so "
            "§6.3's correction returns a finite odds ratio near 1 for an outcome with no contrast in "
            "it -- a plausible number rather than a legible failure [Stage 9 §6.2]. It is FitError "
            "and not SchemaError -- DECISION 6, PI 2026-08-25 -- because an outcome with five events "
            "in ninety-two having none in a resample is a sparse replicate rather than a bug: [§10] "
            "drops and counts it, which is what the objection to a plausible number asks for "
            "[Stage 10 §5.4].")


# --- the weighted risk difference [§8, §5] --------------------------------------------------------------

def weighted_rd(y: np.ndarray, a: np.ndarray, w: np.ndarray) -> float:
    """The [§8] weighted risk difference: P_w(Y = 1 | bridging) - P_w(Y = 1 | EVT alone).

    `share[_TREATED] - share[_COMPARATOR]` and not `share[1] - share[0]`, for the reason
    `cumulative_rd` gives at outcome.py:190-194: the arm codes are named once at module scope from
    TREATMENT_LABELS, and this is one of the two lines in the module whose SIGN is a function of
    which literal goes first. Two bare integers in a subtraction are the one place a reader checking
    the orientation has nothing to check against.

    `higher_is_better` enters NO arithmetic here and the sign convention is uniform across all seven
    outcomes (§4.1): the sign is the CONTRAST's, not the outcome's, and orienting it per outcome would
    make a table of seven numbers in which the meaning of "positive" varies by row. Stage 14 [§16] is
    where a direction is attached to a number.

    `nan` propagates from `weighted_proportion` rather than being caught. An arm with no positive
    weight makes the difference undefined and that is what nan means; S7 is what decides whether
    such a frame reaches here at all (§4.4).

    Public, and it is two lines, for §3.2's reason: §15.8's identity test calls it directly against
    `augmented_rd`, and a test that can only reach it through a seven-outcome loop over a workbook is
    a test of the loop.
    """
    share = weighted_proportion(y, a, w)
    return share[_TREATED] - share[_COMPARATOR]


# --- the weighted marginal odds ratio, and the cell that makes it infinite [§8, §6] ----------------------

def marginal_odds_ratio(y: np.ndarray, a: np.ndarray,
                        w: np.ndarray) -> tuple[float, bool, dict[int, float]]:
    """The [§8] weighted marginal odds ratio, its correction flag, and the raw proportions.

    MARGINAL, NOT CONDITIONAL. [§8]'s first line is "All estimates are marginal in the overlap
    population. No covariate adjustment of the reported effect", and this is a contrast of two
    marginal proportions and never a coefficient from a model — which is one `firth` call away in this
    very module, so §15.4 asserts the two differ rather than leaving it to be noticed.

    Returns THREE things from one call, for Stage 8 §8.1's reason: an odds ratio alone cannot be
    checked against anything, and a reader given only the corrected value cannot see which cell was
    empty. §10.3's table prints all three because of this signature, not despite it.

    The correction is Haldane-Anscombe: OR_CONTINUITY added to all four weighted pseudo-counts --
    the per-arm weight-sums of events and non-events. It fires ONLY when a proportion is degenerate,
    so an interior estimate is bit-for-bit what §6.1 specifies and is not perturbed to buy
    smoothness it does not need.

    THE nan ROW IS NOT CORRECTED AND MUST NOT BE. A missing proportion is not a boundary proportion:
    there is no cell to add half of anything to, and a correction there would manufacture an odds
    ratio for an arm with no data. The branch tests `isnan` FIRST for that reason, and a wrong
    implementation branching on `not (0 < p < 1)` — which nan satisfies — is what §15.5's companion
    shows manufacturing one.

    **IT DOES NOT NECESSARILY SHRINK TOWARD THE NULL, AND AN EARLIER FORM OF THIS DOCSTRING SAID IT
    DOES.** The four pseudo-counts are WEIGHT-SUMS, and the two arms' weight totals are not equal
    across replicates, so the same additive 0.5 is a different relative nudge in each arm. When the
    arm holding the empty cell is the HEAVIER one, the corrected odds ratio can land on the far side
    of 1 from the uncorrected limit -- measured, an empty treated event cell with Sw1 = 17.007
    against Sw0 = 22.768 and p0 = 0.00647 returns 1.0201 where the uncorrected value is 0.0. §6.4
    carries the measurements and §17 files what it implies for a decision that is PI-reversible.
    What IS guaranteed is only that the result is finite and positive.

    `0/0` and `inf/inf` are UNREACHABLE: the degenerate branch is taken before any division, so a
    constant outcome returns a finite number — measured, 1.3529411765 all-zero and 0.7391304348
    all-one — and not a nan. That is S8's real justification and it is the stronger one (§6.2).
    """
    share = weighted_proportion(y, a, w)
    p1, p0 = share[_TREATED], share[_COMPARATOR]
    if np.isnan(p1) or np.isnan(p0):
        return np.nan, False, share
    if min(p1, p0) > 0.0 and max(p1, p0) < 1.0:
        return (p1 / (1.0 - p1)) / (p0 / (1.0 - p0)), False, share
    c = C.OR_CONTINUITY
    odds = {}
    for code in C.TREATMENT_LABELS:
        arm = a == float(code)
        ev = float(w[arm][y[arm] == 1.0].sum()) + c
        non = float(w[arm][y[arm] == 0.0].sum()) + c
        odds[code] = ev / non
    return odds[_TREATED] / odds[_COMPARATOR], True, share


# --- the outcome regression m_a(X) [§8, §7] --------------------------------------------------------------

def outcome_model(df: pd.DataFrame, in_estimate: pd.Series,
                  outcome: str) -> tuple[model.Fit, pd.DataFrame, tuple[str, ...]]:
    """m_a(X) for `outcome`: a Firth logistic fit of treatment plus its declared covariate list.

    [§8]: "treatment main effect plus the §6 covariate set, linear terms, fitted by Firth logistic. No
    interactions, no splines, no cross-fitting — at this treated-arm size flexible nuisance models add
    variance rather than robustness, and cross-fitting folds would be too small to serve their
    purpose."

    Returns the fit, the design it was fitted on, and the columns `design` dropped as constant. The
    design comes back because §7.3's counterfactuals are built by COPYING it and overwriting one
    column -- never by rebuilding from the covariate list, which would not know what was dropped.

    The covariate list comes from `C.outcome_model_covariates(outcome)` and never from
    C.OUTCOME_COVARIATES directly (config.py:489-503). That accessor defaults to the shared [§6]
    entry and returns an override only for an outcome that declares one, which is what makes the two
    nuisance models unable to drift apart: [§8] says "the cross-reference, rather than a second list,
    keeps the two sets from drifting apart in later revisions".

    Treatment is inserted FIRST and by name from C.TREATMENT. `design` never returns it -- the
    covariate list does not contain it, and [§8] calls it a main effect added to that list -- so
    inserting it here is the whole of "treatment main effect". `firth` prepends the intercept, so a
    design matrix is never ambiguous about whether it has one [Stage 6 §5.2].

    `model.firth` raises `FitError` and never returns a fallback (model.py:393). Stage 9 inherits that
    and adds nothing: an m_a(X) that cannot be fitted is a FitError, and [§10] drops and counts the
    replicate rather than substituting the unaugmented estimator (§9.5).
    """
    sub = df.loc[in_estimate]
    X, dropped = model.design(sub, C.outcome_model_covariates(outcome))
    X = X.copy()
    X.insert(0, C.TREATMENT, sub[C.TREATMENT].to_numpy(dtype=float))
    return model.firth(X, sub[outcome].to_numpy(dtype=float)), X, dropped


def _counterfactuals(fit: model.Fit, X: pd.DataFrame) -> tuple[np.ndarray, np.ndarray]:
    """m_1 and m_0 on every row of `X`: the design with the exposure set to each declared arm.

    A COPY with one column overwritten, never a rebuild from the covariate list, because `design`
    dropped constant columns and the list does not know which (§7.1). Every other column keeps its
    observed value: [§8]'s m_a(X) is the conditional mean given the covariates AT their observed
    values, and only the exposure is intervened on.

    `Fit.p` is NOT either of these and an implementation reading it for either is wrong on every row
    assigned to the other arm: `p` is the fitted probability at the OBSERVED treatment, because the X
    handed to `firth` carried the observed treatment column (model.py:104, §7.3).
    """
    out = []
    for code in (_TREATED, _COMPARATOR):
        frame = X.copy()
        frame[C.TREATMENT] = float(code)
        out.append(model.predict(fit, frame))
    return out[0], out[1]


# --- the augmentation [§8, §8] -----------------------------------------------------------------------
#
#   tau  =  Σ₁ w(Y − m₁)/Σ₁ w  −  Σ₀ w(Y − m₀)/Σ₀ w  +  Σ h(m₁ − m₀)/Σ h        with h = e(1 − e)


def _tilt(e: np.ndarray) -> np.ndarray:
    """The [§7] ATO tilting function, h = e(1 - e). One home, one citation, one place to be wrong.

    §0.1's rule applied to the tilting function: `weighted_proportion` is here rather than inside
    `cumulative_rd` because "a second implementation in a later stage would be a second definition
    of the denominator", and Stage 6 §3 and Stage 7 §0.1 made the same move for `ess` and `smd`.
    `h` is the one quantity in [§8] that was left as a sub-expression, and it is the quantity that
    decides WHICH weighted average treatment effect is being estimated -- so it is the one where a
    second definition costs most (§8.3).

    **It exists because it was inline and untested.** Written as a bare `e * (1.0 - e)` at
    `secondary`'s one tilt line, nothing named it and nothing tested the call site: every test that
    distinguishes h from w calls `augmented_rd` directly with arrays the test constructed, so an
    implementer typing `w` on that one line shipped the wrong estimand with the whole suite green.

    `w` is arm-specific — `1 - e` treated and `e` control [§7] — and `h` is not; `h` is the ATO's own
    tilting function and the target population [§7] defines. Measured on the workbook, Σw = 27.736623
    against Σh = 14.797901: the two differ by roughly a factor of two and are NOT proportional, which
    is what makes §15.8 item 3 able to separate them at all.

    NOT a method on Propensity and not a field of it: `propensity.py` is not amended by this stage
    (§13), and h is read by exactly one consumer. If a second consumer appears, moving it to
    Propensity is the right change and §16 records the trigger.

    Takes `e` and not a Propensity, so it is callable on the raw array §11 already sliced and on a
    synthetic one a test constructs.
    """
    return e * (1.0 - e)


def augmented_rd(y: np.ndarray, a: np.ndarray, w: np.ndarray,
                 h: np.ndarray, m1: np.ndarray, m0: np.ndarray) -> float:
    """The [§8] MODEL-ASSISTED augmented risk difference. Never "doubly robust" (§8.3).

    `h` is passed in and is not computed here, because a tilting function is not this function's to
    choose: `augmented_rd` is arithmetic over six arrays and §8.2a's `_tilt` is the one place
    h = e(1-e) is written. An earlier draft of this docstring said h "is `Propensity`'s", which is
    false -- `Propensity` carries e, w and in_model and no h -- and the falsehood mattered, because
    it let the expression live inline at one call site with nothing asserting it (§8.2a).

    The three terms are summed in the order [§8] writes them. Two arm terms, each normalised by its
    OWN arm's weight total, then one correction term normalised over BOTH arms -- h is not
    arm-specific, so its denominator is the whole [§11] population and not either arm's. An
    implementation normalising the correction by an arm total is caught by §15.8 item 1, not by item
    4: measured, the arm-normalised form returns tau - rd = +0.482 on the two-constant test, which is
    a decisive failure. An earlier draft of this docstring cited "§15.8's third test" for this, which
    is the tilt test, and §15.8 item 4 turns out to test nothing item 1 does not (§22.3 item 13).

    THE THREE DENOMINATORS ARE CHECKED BEFORE THEY DIVIDE, for `weighted_proportion`'s reason
    (§3.3): this function is public and separately callable (§3.2), so S7 is not between it and a
    caller, and a raw divide returns nan with a RuntimeWarning rather than saying which total was
    empty. An earlier draft divided raw and §5.1's claim that Stage 9 "inherits the check" was true
    of `weighted_rd` and false of this function (§22.3 item 12).
    """
    t, c = a == float(_TREATED), a == float(_COMPARATOR)
    totals = {"treated": float(np.sum(w[t])), "control": float(np.sum(w[c])),
              "tilt": float(np.sum(h))}
    empty = sorted(k for k, v in totals.items() if not v > 0.0)
    if empty:
        raise C.SchemaError(
            f"augmented_rd was given zero total weight for: {', '.join(empty)}. "
            f"Totals are {totals}. np.sum over an empty selection is 0.0 and the division would "
            "return nan with a warning rather than naming which term was undefined "
            "[Stage 9 §8.1]. S7 is the guard on the `secondary` path; this is the guard on the "
            "public one.")
    return (float(np.sum(w[t] * (y[t] - m1[t])) / totals["treated"])
            - float(np.sum(w[c] * (y[c] - m0[c])) / totals["control"])
            + float(np.sum(h * (m1 - m0)) / totals["tilt"]))


# --- the rare-minority guard, decided ONCE [§9] ---------------------------------------------------------

def _augmentable(outcome: str, minority: int) -> bool:
    """Whether `outcome` is augmented [§8, as amended 2026-08-07].

    Two ways in, and the second is not a special case of the first:

      * the minority cell supports the shared [§6] covariate list, or
      * the outcome DECLARES a reduced m_a(X) in C.OUTCOME_MODEL_OVERRIDES.

    The amendment's item 2 is explicit that a declared reduction is augmented "rather than dropped
    from augmentation", so an override is a route IN and not a modifier of an outcome that qualified
    anyway. The roadmap said such outcomes "are not augmented"; that sentence was stale and §22
    item 1 amends it (§7.2).

    `min(events, non-events)` and NOT the event count, on the outcome's OWN [§11] population — which
    has to be said because the two available answers differ, and the coincidence that they agree on
    the non-event count is exactly what would let a wrong choice go unnoticed (§9.1). The rule is
    applied on the [§11] population because that is the population the nuisance model is FITTED on,
    and the constraint it encodes is that a model cannot carry more parameters than its own fitting
    sample supports.

    `C.RARE_MINORITY_THRESHOLD` is read and never compared against a literal.

    Reading the registry rather than a second list of names is what makes this unable to drift: an
    outcome cannot become reduced except by being named in config.py, and cannot be named there
    without becoming augmentable here.
    """
    return (minority >= C.RARE_MINORITY_THRESHOLD
            or outcome in C.OUTCOME_MODEL_OVERRIDES)


def _augmented_path(outcome: str, minority: int) -> str:
    """"full" | "reduced" | "unaugmented" -- the three declared paths, named once."""
    if not _augmentable(outcome, minority):
        return "unaugmented"
    return "reduced" if outcome in C.OUTCOME_MODEL_OVERRIDES else "full"


# --- the three `model` audit entries [§10] ---------------------------------------------------------------
#
# No new kind: `data.KINDS` stays nine, which is `data.py:151-152`'s own prediction honoured for the
# third stage running. Stage 7 added two entries under the existing `model` kind, Stage 8 three, and
# Stage 9 adds three — taking the ledger to
#
#   load 7 / derive 4 / classify 1 / build 6 / fit 4 / assess 2 / primary 3 / secondary 3  =  30
#
# NO CELL OF ANY OF THE THREE TABLES CONTAINS A `|`, and that is why §10.4's column is `max_abs_beta`
# and not `max|beta|`. `data._md_table` (data.py:242-245) does no escaping and sizes its separator from
# `len(rows[0])`, so two pipes in a header cell present eleven markdown cells against a nine-cell
# separator and the whole table renders misaligned. Stage 7 shipped this exact defect twice — `|SMD| <
# 0.1` and `worst |SMD|`, rendering 8/6 and 15/13 — and it is INVISIBLE to the byte-identity criterion,
# because a table identically broken under both hash seeds is still byte-identical. §15.14 asserts the
# structural property on the rendered text instead.
#
# What these entries deliberately do NOT carry: no coefficient vector — m_a(X) is a nuisance and [§8]'s
# first line forbids reporting adjusted effects, so printing its treatment coefficient invites exactly
# the misreading that a conditional odds ratio is the estimate — and no standard error, no interval and
# no p-value anywhere. [§10]'s percentile bootstrap is the prespecified interval and there is no second.


def _estimable_table(df: pd.DataFrame, ps: propensity.Propensity,
                     populations: dict[str, pd.Series],
                     minorities: dict[str, int]) -> tuple[tuple[str, ...], ...]:
    """§10.2's seven rows. Takes `ps` because `lost` is not derivable without in_model (§10.2).

    `lost` is |in_model| - |in_model & notna|, and it is NOT df[key].isna().sum(): that counts over
    the cohort rather than over the ATO population, and the two differ whenever a
    covariate-incomplete record is also outcome-missing. On v7 they happen to agree for six of the
    seven, which is exactly the coincidence that would let the wrong expression pass unnoticed.

    Ranges over C.BINARY_OUTCOMES so the row order is the registry's, which is what makes the
    rendered log byte-identical across hash seeds (§15.14).

    `df` is in the signature and unread, matching `_estimable_detail` beside it: the two are one
    table and one caption over the same four objects, and a helper pair whose signatures differ by a
    parameter invites a caller to pass them different things.
    """
    in_model = int(ps.in_model.sum())
    header = ("outcome", "family", "n [§11]", "lost", "minority", "path", "m_a(X)")
    rows: list[tuple[str, ...]] = [header]
    for key in C.BINARY_OUTCOMES:
        n = int(populations[key].sum())
        path = _augmented_path(key, minorities[key])
        covariates = ("--" if path == "unaugmented"
                      else ", ".join(C.outcome_model_covariates(key)))
        rows.append((key, C.OUTCOMES[key].family, str(n), str(in_model - n),
                     str(minorities[key]), path, covariates))
    return tuple(rows)


def _estimable_detail(df: pd.DataFrame, ps: propensity.Propensity,
                      populations: dict[str, pd.Series]) -> str:
    """The three populations in the Stage 8 §4.1 form: the cohort's, the ATO's, then per outcome.

    Stated as a RANGE over the seven rather than seven numbers, because the table already carries
    them per outcome and a detail string that repeats a table is a second place to disagree with it.
    """
    sizes = sorted({int(populations[k].sum()) for k in C.BINARY_OUTCOMES})
    span = str(sizes[0]) if len(sizes) == 1 else f"{sizes[0]}-{sizes[-1]}"
    return (f"{len(df)} cohort record(s); {int(ps.in_model.sum())} in the [§7] overlap population; "
            f"{span} on the seven [§8] outcome populations, one per outcome in the table. Each "
            f"estimate uses its OWN denominator [§11] and no outcome is estimated on another's.")


def _estimates_table(estimates: dict[str, BinaryEstimate]) -> tuple[tuple[str, ...], ...]:
    """§10.3's seven rows: both proportions, RD, OR, the correction flag, tau, the m_a(X) list.

    The unaugmented risk difference and tau are ADJACENT columns, and that is a reporting requirement
    rather than a layout choice. [§8]: "the augmented estimate is therefore reported as
    model-assisted, always alongside the unaugmented weighted risk difference; material disagreement
    between the two is evidence about the outcome model." Adjacency is what makes the comparison the
    protocol requires available without arithmetic.

    `corrected` is rendered for every row even when it is False for all seven, which on v7 it is: a
    column that disappears when nothing triggers it is a column a reader cannot tell was checked
    (§10.3). `tau` renders `--` and never 0.0 for an unaugmented outcome, because the deliberateness
    of an absence lives in a mask and never in a value [Stage 6 §9].

    The last column is `m_a(X)` and NOT `reduced`: `reduced` is the name of a boolean FIELD on
    BinaryEstimate and this column renders a tuple of covariate NAMES. A column and a field sharing a
    name while carrying different types is the kind of collision that survives review and then costs
    an afternoon. The amendment requires the specification and not a `True`: "Every reduced
    specification is reported beside its estimate."

    NO CELL MAY CONTAIN A PIPE and none does: `data._md_table` does no escaping (§10.3, §15.14).

    An outcome absent from `estimates` is SKIPPED, which under the point-estimate contract is none of
    them: only `secondary(..., collect=True)` can produce a partial dict, and that is Stage 10's
    replicate call whose `Audit` is a throwaway nobody writes [Stage 10 §6.2, §7.3]. Rendering a row
    of dashes for it would put a row about nothing in a log nobody reads, and rendering the message
    would put a `FitError`'s prose — which is not pipe-checked — in a markdown cell.
    """
    header = ("outcome", "n", "p1_w", "p0_w", "RD_w", "OR_w", "corrected",
              "tau (model-assisted)", "m_a(X)")
    rows: list[tuple[str, ...]] = [header]
    for key in C.BINARY_OUTCOMES:
        if key not in estimates:
            continue
        est = estimates[key]
        rows.append((
            key, str(int(est.in_estimate.sum())),
            _fmt(est.proportion[_TREATED]), _fmt(est.proportion[_COMPARATOR]),
            _fmt(est.rd), _fmt(est.odds_ratio), str(est.or_corrected),
            "--" if est.augmented is None else _fmt(est.augmented),
            "--" if est.covariates is None else ", ".join(est.covariates)))
    return tuple(rows)


def _models_table(estimates: dict[str, BinaryEstimate]) -> tuple[tuple[str, ...], ...]:
    """§10.4's one row per AUGMENTED outcome -- five on v7, not seven.

    The two unaugmented outcomes have no model and a row asserting that would be a row about
    nothing. `max_abs_beta` is spelled in ASCII: `max|beta|` would put two pipes in a header cell and
    break the rendered table in a way byte-identity does not catch (§10.4, §15.14).

    `max_abs_beta` is in the table for §9.6's reason: no bound rejects a fit here, so the only thing
    standing between a separated nuisance model and an unremarked estimate is a number in the log
    that a reader can see. The three safeguard counters are here for Stage 6 §3.2's reason -- "a
    safeguard whose activation nobody can count is a safeguard nobody can evaluate" -- and `dropped`
    because `center_USZ` vanishing from every design is a fact about the cohort nobody should
    rediscover.
    """
    header = ("outcome", "covariates", "k", "rows", "iters", "max_abs_beta",
              "rescales", "halvings", "dropped")
    rows: list[tuple[str, ...]] = [header]
    for key in C.BINARY_OUTCOMES:
        if key not in estimates:            # collect=True dropped it — `_estimates_table`'s reason
            continue
        est = estimates[key]
        if est.fit is None:
            continue
        rows.append((
            key, ", ".join(est.covariates or ()), str(len(est.fit.columns)),
            str(int(est.in_estimate.sum())), str(est.fit.iterations),
            _fmt(float(np.max(np.abs(est.fit.beta)))),
            str(est.fit.rescales), str(est.fit.halvings),
            ", ".join(est.dropped) if est.dropped else "--"))
    return tuple(rows)


def _record_estimable(df: pd.DataFrame, ps: propensity.Propensity,
                      populations: dict[str, pd.Series], minorities: dict[str, int],
                      audit: Audit) -> None:
    """§10.2's entry. `n` is the COHORT's, because Audit.record takes one number and there are seven.

    Recorded BEFORE any estimate is computed, so that a failure in the estimator leaves a log that
    still names the populations. This is `propensity._record_exclusion`'s ordering rule
    (propensity.py:441-443) applied one stage on: "the exclusion is a fact about the frame,
    established before any modelling choice, and a log that names it is useful precisely when what
    follows fails."

    No `case_ids`: `_MUST_NAME_CASES` does not cover the `model` kind (data.py:196, 215) and this
    entry removes no patient from anything -- the seven masks are subsets of in_model, whose
    exclusion `propensity._record_exclusion` already named. An entry naming the same excluded patient
    a second time is not a second fact (§10.2).
    """
    audit.record("model", "binary_estimable", len(df),
                 _estimable_detail(df, ps, populations),
                 table=_estimable_table(df, ps, populations, minorities))


def _record_estimates(estimates: dict[str, BinaryEstimate], audit: Audit) -> None:
    """§10.3's entry. `n` is the number of [§8] binary estimates, which is seven and not a population.

    The detail names what the table cannot: that the unaugmented risk difference and tau sit in
    ADJACENT columns because [§8] requires the comparison, and that "model-assisted" is the licensed
    description (§8.3). §15.13 scans this string for the forbidden phrase.
    """
    augmented = sum(1 for e in estimates.values() if e.augmented is not None)
    corrected = sum(1 for e in estimates.values() if e.or_corrected)
    audit.record(
        "model", "binary_estimates", len(estimates),
        f"{len(estimates)} [§8] binary outcome(s), each on its own [§11] denominator. "
        f"{augmented} carry a model-assisted augmented risk difference beside the unaugmented one, "
        "in adjacent columns, because [§8] requires the comparison and calls material disagreement "
        "between them evidence about the outcome model rather than confirmation of either. "
        f"§6.3's continuity correction fired on {corrected} of them. No interval and no p-value: "
        "[§10] owns those.",
        table=_estimates_table(estimates))


def _record_models(estimates: dict[str, BinaryEstimate], audit: Audit) -> None:
    """§10.4's entry. `n` is the number of nuisance models FITTED -- five on v7, not seven.

    No coefficient vector, for §10.4's reason: `m_a(X)` is a nuisance and [§8]'s first line forbids
    reporting adjusted effects, so printing its treatment coefficient invites the misreading that a
    conditional odds ratio is the estimate. No standard error and no interval anywhere.
    """
    fitted = [e for e in estimates.values() if e.fit is not None]
    audit.record(
        "model", "binary_outcome_models", len(fitted),
        f"{len(fitted)} Firth nuisance model(s) for the augmented outcomes; "
        f"{len(estimates) - len(fitted)} outcome(s) are unaugmented under §9.2 and have none. "
        "max_abs_beta is reported because no bound rejects a fit [§9.6], so it is the only thing "
        "between a separated nuisance model and an unremarked estimate. No coefficient vector, no "
        "standard error, no interval: only the predictions of these models enter any estimate.",
        table=_models_table(estimates))


# --- the stage [§8, §11] ---------------------------------------------------------------------------------

def secondary(df: pd.DataFrame, ps: propensity.Propensity, audit: Audit,
              paths: dict[str, str] | None = None, collect: bool = False) -> Secondary:
    """The [§8] binary estimates: seven outcomes, each on its own [§11] denominator.

    Takes no outcome list: `C.BINARY_OUTCOMES` is computed from the [§5] registry, so a ninth outcome
    added to `OUTCOMES` is estimated here in the same edit and cannot fail to be (§4.1).

    `paths` is Stage 10's frozen point-estimate decision (§9.5, §11.1). None means "decide from this
    frame", which is correct on the point estimate and WRONG in a replicate: `ph2` crosses the
    threshold in 46.1% of replicates and [§10] forbids mixing estimators across one interval. An
    earlier draft omitted the parameter entirely, which made §12.1's and §14's "the paths must be
    passed IN" unsatisfiable through this function (§22.3 item 11). It is optional because the
    point-estimate call is the one that DECIDES the paths, and complete-or-absent because a partial
    map is the two-estimator mixture arrived at one outcome at a time.

    `collect` turns a per-outcome `model.FitError` into a recorded failure instead of a raise
    [Stage 10 §7.3]. Default False, which is the point-estimate contract unchanged: on the workbook a
    `FitError` is fatal and must be, because there is no replicate to drop. Stage 10 passes True, and
    then a single outcome's unfittable `m_a(X)` leaves the other six estimated — which is [§10]'s
    "dropped and counted" at the granularity [§10] left open and Stage 10 §7.3 chose. The key lands
    in `Secondary.failures` with its message, because Stage 10 classifies the failure by the leading
    token of that message (Stage 10 §7.2) and a bare absence carries no token.

    **`collect` NEVER catches `C.SchemaError`, at either setting.** S6 and S7 stay fatal, which is
    Stage 10 §5.4's whole split: an empty [§11] population or an arm with no weight means the mask
    and the weights disagree about the same rows, and that is a bug rather than a sparse replicate.
    S8 alone was reclassified, and it is a `FitError` now, so it is what `collect` collects.

    **The scope of the `except` is ONE OUTCOME'S WHOLE ESTIMATE and not its augmentation**
    [Stage 10 §7.3]. An `m_a(X)` failure costs that outcome's `rd`, `odds_ratio` AND `augmented`
    together, because Stage 9 §14 forbids taking percentiles of `tau` and of `rd` from different
    replicate sets — [§10.3] reports the two as a comparison. Catching around `outcome_model` alone
    would leave `rd` estimated and `tau` absent, which is exactly that prohibition violated.

    Six things about the order below, each of which is a failure if moved:

      * `_assert_readable` runs FIRST, raising before any value is read, so a frame missing an
        outcome column reports S1 rather than a bare pandas KeyError from S2 (§4.4a).
      * `_assert_secondary_inputs` runs SECOND and over the whole frame, before any outcome is
        touched, so a frame that fails S2-S5 reports that rather than a message about one outcome.
      * `_record_estimable` runs BEFORE the estimator loop, so a raise inside an estimate still
        leaves a log naming the seven populations (§10.2).
      * `in_estimate` is built per outcome and `e`, `w` are sliced to IT, never to `in_model`. This is
        the one place either is read, which is how Stage 6 §9's "range over the mask, never over
        notna()" is honoured -- as a fact about two lines rather than as a warning (§12).
      * The minority cell is counted on `in_estimate` and not on the cohort (§9.1).
      * `_augmentable` is consulted BEFORE `outcome_model` is called, so an unaugmented outcome never
        fits a nuisance model. Reversing them would make `sich`'s 17-iteration fit (§7.4) happen on
        every call and be discarded, and would put a FitError on a path [§8] says has no model.

    NO FitError IS SWALLOWED HERE, and that is the whole of §9.5's no-substitution rule as code. A
    `model.firth` failure on a declared nuisance model is never downgraded to the unaugmented form:
    at `collect=False` it propagates, and at `collect=True` it removes the outcome's whole estimate
    and is recorded. Neither route produces a number from a fallback estimator, which is roadmap
    invariant 5. A try/except around `outcome_model` that let `rd` survive is the most natural thing
    an implementer writes when a fit fails on 0.3% of `sich` replicates and it is the thing that
    turns this stage into a two-estimator mixture. §15.7a asserts it and shows the caught version
    passing every other test.

    Deterministic, and writes no file. No seed is held, no clock is read, and the only mutable object
    touched is the `Audit` passed in. Stage 10 must pass a throwaway one per replicate, for Stage 8
    §11's reason and more so: three entries per replicate at N_BOOT = 2000 is 6000 entries, and one
    of them carries a seven-row rendered table.

    Adds no column to `df`, edits no value, removes no row and re-weights nothing (§0.2). Reads
    neither `Primary` nor `Balance`.
    """
    _assert_readable(df, ps)
    _assert_secondary_inputs(df, ps)
    if paths is not None and set(paths) != set(C.BINARY_OUTCOMES):
        raise C.SchemaError(
            f"secondary was given frozen paths for {sorted(paths)} against "
            f"{sorted(C.BINARY_OUTCOMES)}. A partial map would let some outcomes carry the "
            "point estimate's path and others re-decide from the replicate, which is the "
            "two-estimator mixture §9.5 forbids arrived at one outcome at a time.")

    populations: dict[str, pd.Series] = {}
    minorities: dict[str, int] = {}
    for key in C.BINARY_OUTCOMES:
        populations[key] = ps.in_model & df[key].notna()
        y = df.loc[populations[key], key].to_numpy(dtype=float)
        minorities[key] = min(int((y == 1.0).sum()), int((y == 0.0).sum()))
    _record_estimable(df, ps, populations, minorities, audit)

    estimates: dict[str, BinaryEstimate] = {}
    failures: dict[str, str] = {}
    for key in C.BINARY_OUTCOMES:
        try:
            in_estimate = populations[key]
            sub = df.loc[in_estimate]
            y = sub[key].to_numpy(dtype=float)
            a = sub[C.TREATMENT].to_numpy(dtype=float)
            w = ps.w.loc[in_estimate].to_numpy(dtype=float)
            e = ps.e.loc[in_estimate].to_numpy(dtype=float)
            _assert_estimable(key, y, a, w)          # S6, S7 SchemaError; S8 FitError [§5.4]

            odds_ratio, corrected, share = marginal_odds_ratio(y, a, w)
            path = paths[key] if paths is not None else _augmented_path(key, minorities[key])
            tau: float | None = None
            fit: model.Fit | None = None
            covariates: tuple[str, ...] | None = None
            dropped: tuple[str, ...] = ()
            if path != "unaugmented":
                covariates = C.outcome_model_covariates(key)
                fit, X, dropped = outcome_model(df, in_estimate, key)
                m1, m0 = _counterfactuals(fit, X)
                tau = augmented_rd(y, a, w, _tilt(e), m1, m0)

            estimates[key] = BinaryEstimate(
                outcome=key, family=C.OUTCOMES[key].family, minority=minorities[key],
                rd=weighted_rd(y, a, w), odds_ratio=odds_ratio, or_corrected=corrected,
                proportion=dict(share), augmented_path=path, augmented=tau, covariates=covariates,
                reduced=key in C.OUTCOME_MODEL_OVERRIDES, dropped=dropped,
                in_estimate=in_estimate, fit=fit)
        except model.FitError as failure:
            # `collect=False` re-raises, which is the point-estimate contract unchanged. Only
            # `model.FitError` reaches here: `C.SchemaError` is not a subclass of it, so S6 and S7
            # propagate at both settings and no `except` in this module can catch them
            # [Stage 10 §7.1, §15.7].
            if not collect:
                raise
            failures[key] = str(failure)

    _record_estimates(estimates, audit)
    _record_models(estimates, audit)
    return Secondary(estimates=estimates, failures=failures)
