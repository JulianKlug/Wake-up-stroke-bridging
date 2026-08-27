"""Stage 6a — the pipeline's numerics [§6, §7].

The design matrix and the Firth fit, and nothing that knows what the exposure is. This module is
**outcome-agnostic by construction**: it reads ``config`` and nothing else, takes no ``Audit``, and
names neither the treatment nor any covariate list. Its caller logs; it computes.

**``predict`` does not bend that rule, and this sentence is here because it looks as though it might**
[Stage 9 §7.3]. Stage 9 builds its counterfactuals ``m_1`` and ``m_0`` by copying a design, overwriting
one column and calling ``predict`` — so a reader who greps for the function finds it being used to
intervene on an exposure. ``predict`` cannot tell an exposure from a covariate: it takes a ``Fit`` and
a frame, checks the column names agree as a sequence, and evaluates. The knowledge that one of those
columns is a treatment lives in ``outcome.py``, where the exposure is already named.

**Stages 8, 9 and 12 come into this module directly and never through ``propensity.py``.** Stage 9's
outcome regression ``m_a(X)`` is a Firth logistic fit on a design matrix over
``outcome_model_covariates(outcome)`` — the same two functions, a different covariate list and a
different response — and Stage 12's proportional-odds model is a design matrix over
``STANDARDISATION_COVARIATES``. If ``design`` and ``firth`` lived in ``propensity.py``, Stage 9 would
import a module named for the exposure in order to fit an **outcome** model, and the first reader to
tidy that import would be right to. That sentence is the one thing standing between this module and
somebody moving it there for tidiness [Stage 6 §0.1].

Section references in brackets are to ``statistical_analysis_plan.md``. The specification for this
module is ``specs/stage6_propensity_and_weights.md``; nothing here is invented outside it.

Where it sits::

    cohort.build(df, audit)  →  93 rows x 33 columns                    [Stage 5 §11]
                       │
                       ▼
    ┌──────────────────────────────────────────────────────────────────────────┐
    │  STAGE 6a — model.py                              reads config.py only   │
    │                                                                          │
    │   design(df, covariates)        [§6] reference-coded design matrix       │
    │     ├─ _assert_design_inputs    D1…D4 → SchemaError                      │
    │     ├─ declared FACTOR_LEVELS → Categorical → reference-coded    (§4.2)  │
    │     ├─ constant columns dropped, AFTER dummying                  (§4.3)  │
    │     └─ returns (X, dropped)                                              │
    │                                                                          │
    │   complete_cases(df, covariates)   the [§11] mask, RETURNED not applied   │
    │                                                                          │
    │   firth(X, y)                   [§7] penalised logistic          (§5)    │
    │     ├─ _assert_fittable         F3 rank / F4 X / F6 y → FitError (§5.4)  │
    │     ├─ l(b) + ½log|I(b)|,  score X'(y − p + h(0.5 − p))                  │
    │     ├─ step-halving on the PENALISED likelihood                          │
    │     ├─ converge on Δpll  OR a flat modified score                (§5.3)  │
    │     └─ raise FitError — never a second estimator      [invariant 5]      │
    └──────────────────────────────────────────────────────────────────────────┘
         │                                       ▲                 ▲
         │                         Stage 9 [§8] m_a(X)   Stage 12 [§14a]
         ▼                         enter HERE and never through propensity.py
    Stage 6b — propensity.py

**Roadmap invariant 4 becomes enforceable here**, and nowhere before: "no post-time-zero variable
appears in any model's design matrix — assert against an explicit denylist, not by inspection".
``design`` is the one function every model's covariates pass through, D2 is that assertion, and
Stages 8, 9 and 12 inherit it by *calling* ``design`` rather than by remembering.

**An observation-weight vector does not breach the no-exposure rule, and this sentence is here so the
next reader does not have to re-derive it** [Stage 8 §0.1, §12]. ``propensity.ess``'s docstring calls
the rule a *"no-exposure-no-weights rule"*, and read against this module's own docstring that clause is
about ``ess`` — a Kish sum, which fits nothing — and not about a fitter. ``firth(X, y)`` already takes
the **exposure** as an opaque ``y`` and is outcome-agnostic because it does not know what ``y`` *is*; a
non-negative observation-weight vector is opaque in exactly the same way. The rule being honoured is
"this module names neither the treatment nor any covariate list", and ``polr(X, y, w)`` names none of
the three.

This module is **not** exempt from the Stage 1 §7 raw-name scan and must never become exempt. It names
no raw header, and it names no exposure — §12.12 scans for both.
"""
from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass

import numpy as np
import pandas as pd

import config as C


# --- the one exception type [invariant 5] --------------------------------------------------------
#
# The evidence for the no-fallback rule is in this repository rather than in principle:
# `pilots/analysis.py`'s docstrings record that SAP v1.0's scikit-learn fallback fired in 16.6% of
# bootstrap replicates, so one interval in six was a mixture of two estimators' sampling distributions
# rather than one estimator's. [§7] is written against that measurement.

class FitError(RuntimeError):
    """The prespecified fit did not converge on this sample. [§7]: there is no second estimator."""


# --- what one fit is ------------------------------------------------------------------------------

@dataclass(frozen=True)
class Fit:
    """One Firth fit. `p` is aligned to the rows of the X it was given, positionally.

    The last three fields describe how the fit was *reached* rather than what it is, which is a
    different kind of thing from `beta` and `p` and earns its own justification [Stage 6 §3.2]. The
    precedent is `converged_on`. §5.2's trust region and step-halving are the safeguards that fire
    when a bootstrap replicate is hard, [§10] refits N_BOOT times, and a safeguard whose activation
    nobody can count is a safeguard nobody can evaluate. They are counters, not state: the loop
    computes all three anyway, no downstream estimate depends on them, and without them the four
    safeguard tests of §12.4a cannot be written at all — they are properties of the *iteration*, and
    testing a safeguard by side effect ("does it converge?") passes a broken trust region that happens
    to converge.
    """

    beta: np.ndarray            # (k + 1,), intercept first
    p: np.ndarray               # (n,), fitted probabilities
    iterations: int
    converged_on: str           # "likelihood" | "score"  — §5.3
    columns: tuple[str, ...]    # the design's column names, in beta's order after the intercept
    first_step_norm: float      # ‖step‖ at iteration 1, BEFORE any rescale   ┐
    rescales: int               # steps shortened by the trust region [§5.2a] ├─ §3.2
    halvings: int               # total step-halvings across all iterations   ┘


# --- the design matrix [§6] -----------------------------------------------------------------------
#
#   Reference-coded dummies for the declared factors, linear terms for everything else,
#   constant columns dropped [roadmap Stage 6].
#
#   Factors are CATEGORICAL. Their level sets are FACTOR_LEVELS and their baselines are
#   REFERENCE_LEVELS. No interaction, no spline, no transform — [§6] says "linear terms only".


def _assert_design_inputs(df: pd.DataFrame, covariates: Sequence[str]) -> None:
    """D1-D4, collected: one SchemaError reports every one that failed [Stage 6 §4.5].

    Following Stage 2's `_assert_schema`, Stage 4's `_assert_classifier_inputs` and Stage 5's
    `_assert_cohort_inputs`. Every branch is unreachable on the workbook — every centre is in
    CENTER_ORDER, every onset_type is in its declared levels, and the [§6] covariates are disjoint
    from POST_TIME_ZERO, which test_config.py already asserts — and that is recorded rather than
    treated as a reason to skip one. They are written for the workbook that has not arrived yet, and
    for the covariate list a later stage passes in.
    """
    bad: list[str] = []

    absent = [c for c in covariates if c not in df.columns]
    if absent:
        bad.append(
            f"D1  {', '.join(absent)}: not a column of the frame. A design matrix over a covariate "
            "the frame does not carry is a model over a different specification than the one "
            "requested. Stage 3's derived columns exist only after derive(); the [§13] subgroups "
            "only after derive_cohort().")

    denied = [c for c in covariates if c in C.POST_TIME_ZERO]
    if denied:
        bad.append(
            f"D2  {', '.join(denied)}: post-time-zero [invariant 4, §12]. This is the assertion "
            "against the explicit denylist that invariant 4 requires, and this function is where it "
            "can be made: every model's covariates pass through here. onset_to_groin_min is "
            "reported by arm and never adjusted for [§12]; an outcome in a design matrix is the "
            "outcome conditioning on itself.")

    # D1 fires before D3 can KeyError on a covariate, so the generator re-filters on `c in df.columns`
    # — the checks are collected, so D1 must not prevent D3 and D4 from running. D4 is emitted before
    # D3 although it is numbered after; the order is inert because the messages are collected.
    for c in (c for c in covariates if c in C.CATEGORICAL and c in df.columns):
        absent_values = df[c].isna()
        if absent_values.any():
            bad.append(
                f"D4  {c}: {int(absent_values.sum())} record(s) carry no value. A missing factor "
                "value is NOT caught downstream and this is the only place it can be: the Categorical "
                "conversion maps it to NaN, get_dummies emits an ALL-ZERO row, and that row is "
                "byte-identical to a record genuinely at the reference level — measured "
                "[Stage 6 §4.5]. There is no nan left for F4 to find and D3 excludes missing values "
                "from its own check by design. Use model.complete_cases() [§4.4]; every caller owns "
                "its own mask [§9].")

        # D3 excludes missing values from its own check because they are D4's, not because they are
        # harmless. Without the notna clause D3 would fire on every incomplete record and the two
        # failures — a level outside the declared set, and no level at all — would report as one, with
        # one message describing the wrong mechanism.
        offside = ~df[c].isin(C.FACTOR_LEVELS[c]) & df[c].notna()
        if not offside.any():
            continue
        # `case_id` is used to NAME the offenders and is not a precondition of building a design.
        # Indexing it unconditionally made `design` raise a bare KeyError on any frame without it —
        # including the single-centre and factor-only frames §12.1 describes, and anything Stage 9
        # or Stage 12 assembles. The count is the assertion; the names are a courtesy, and a courtesy
        # may not impose a schema.
        named = (", ".join(sorted(df.loc[offside, "case_id"])) if "case_id" in df.columns
                 else f"index {sorted(df.index[offside])} — the frame carries no case_id")
        off = df.index[offside]
        if len(off):
            bad.append(
                f"D3  {c}: {len(off)} record(s) carry a level outside FACTOR_LEVELS[{c!r}]: "
                f"{named}. The Categorical conversion turns such a value into NaN "
                f"and get_dummies then encodes it as an all-zero row — which is the encoding of the "
                f"reference level, {C.REFERENCE_LEVELS[c]!r}. The record would be modelled as though "
                "it were at the reference, with no missing value anywhere and nothing to notice "
                "[Stage 6 §3.1].")

    if bad:
        raise C.SchemaError(
            "\n  " + "\n  ".join(bad)
            + f"\n  {len(bad)} design assertion(s) failed over {len(df)} records and "
            f"{len(covariates)} covariate(s).")


def design(df: pd.DataFrame, covariates: Sequence[str]) -> tuple[pd.DataFrame, tuple[str, ...]]:
    """The [§6] design matrix over `covariates`, reference-coded on REFERENCE_LEVELS.

    No intercept column: `firth` prepends its own, so a design matrix is never ambiguous about
    whether it has one. Constant columns are dropped and named in the return value rather than
    logged here, because model.py takes no Audit [Stage 6 §0.1] and its caller is the step.

    Four properties, each load-bearing, and each one a thing `pilots/analysis.py:25-41` does
    differently [Stage 6 §4.2]:

    * **The level set is declared, not observed.** `pd.Categorical(..., categories=FACTOR_LEVELS[c])`
      makes the column set a function of config.py and not of the rows. Without it, get_dummies on the
      bare `string` column emits only the levels present — so an absent level's all-zero column never
      exists, the constant-column rule has nothing to do, and the width becomes a property of the
      sample.
    * **The reference column is dropped by name, never by position — and the reason is not the obvious
      one.** On a *declared* Categorical, `drop_first=True` drops the first *declared* level, which IS
      the reference, so the two forms are equivalent and neither rebaselines; measured. What by-name
      buys is independence from a coupling nothing else states — `drop_first` is correct only while
      `REFERENCE_LEVELS[c] == FACTOR_LEVELS[c][0]`, which test_config.py now asserts. A replicate that
      loses the reference level entirely fails as a rank deficiency (F3), not as a KeyError: the
      declared Categorical guarantees the column exists.
    * **`nunique(dropna=False)`, not `nunique()`.** A column that is 1.0 on one row and nan on the
      rest has `nunique() == 1` and would be dropped as constant while carrying a nan into nothing.
      Moot here because `complete_cases` runs first; written this way so it stays moot if a caller
      ever forgets.
    * **The drop happens after dummying**, so both an all-zero dummy (an absent level) and an all-one
      dummy (a single-level factor in a subset) are removed. A factor with one level present
      contributes no column at all, which is correct: it is collinear with the intercept.

    **`remove_unused_categories` is not called and must not be.** Stage 2 §6 chose `string` for
    `center` precisely so no dead level survives a restriction, and this function creates its
    Categorical fresh from the declared levels on every call. Calling it would delete the all-zero
    absent-level column *before* the constant-column rule could name it, and `dropped` — which the
    audit log publishes — would be silently empty.

    Nothing here mutates a covariate list: the dropping happens to the matrix, and `dropped` reports
    it [Stage 3 §7.2].
    """
    _assert_design_inputs(df, covariates)
    factors = [c for c in covariates if c in C.CATEGORICAL]
    sub = df[list(covariates)].copy()
    for c in factors:
        sub[c] = pd.Categorical(sub[c], categories=C.FACTOR_LEVELS[c])
    X = pd.get_dummies(sub, columns=factors, dtype=float)
    X = X.drop(columns=[f"{c}_{C.REFERENCE_LEVELS[c]}" for c in factors]).astype(float)
    dropped = tuple(c for c in X.columns if X[c].nunique(dropna=False) <= 1)
    return X.drop(columns=list(dropped)), dropped


def complete_cases(df: pd.DataFrame, covariates: Sequence[str]) -> pd.Series:
    """Boolean, total, on `df`'s index: True where every covariate in `covariates` is present.

    [§11] is complete-case *per estimate*, with the denominator reported for each. So this is
    returned for a caller to report and to log, never applied here: `design(df.dropna(...))` would
    make the exclusion a property of a call nobody can see, and the excluded patient would leave
    the analysis with every table still reconciling.

    **Not folded into `design`**, because `design` is called by Stages 9 and 12 over different
    covariate lists and the complete-case set is a property of the list. An outcome model over a
    reduced list may exclude nobody, and a `design` that silently dropped rows would give its caller a
    design matrix and a response of different lengths, aligned by luck.
    """
    return df[list(covariates)].notna().all(axis=1)


# --- the Firth fit [§7] ---------------------------------------------------------------------------
#
# [§7] prescribes one estimator: Firth-penalised logistic regression, refitted identically in every
# bootstrap replicate. The penalised log-likelihood is
#
#     l*(b)  =  l(b)  +  ½ log|I(b)|,        I(b) = X' W X,  W = diag(p(1 − p))
#
# whose stationary point solves the modified score equation
#
#     X'(y − p + h(0.5 − p))  =  0,          h = diag(H),  H = W^½ X (X' W X)⁻¹ X' W^½
#
# The ½log|I(b)| term is the Jeffreys prior, and h(0.5 − p) is what pulls the fit back from a boundary
# an unpenalised likelihood would run to. That is the entire reason [§7] names it: an unpenalised MLE
# fails to converge in a minority of sparse replicates, and [§10] refits 2000 times.
#
# It is implemented here rather than imported because there is no maintained Firth implementation for
# Python — statsmodels has none, the one working package was last released in 2022 and requires a
# dependency Stage 1 excluded, the other is a two-line stub, and R cannot be a dependency of a
# 2000-replicate bootstrap. Investigated rather than assumed; §16b is the record, and what replaces the
# dependency is four oracles rather than a comment. The package names live in the spec and in the test
# that uses one of them, never here: §12.7 scans the repository to keep that true.


def _penalised_loglik(Xc: np.ndarray, y: np.ndarray, beta: np.ndarray) -> float:
    """l(b) + ½log|I(b)|, or -inf where the information matrix is singular.

    -inf rather than a raise, so that step-halving can *reject* a step and try a shorter one; a raise
    there would turn a recoverable step into a failed fit. slogdet returns -inf for a singular matrix
    rather than raising [§3.1], so the branch is about a non-finite logdet and not an exception.

    THE BRANCH DOES NOT FIRE, and that is measured rather than hoped [§3.1, §18]. FIRTH_WEIGHT_FLOOR
    floors w, so X'WX >= FIRTH_WEIGHT_FLOOR * X'X, which is nonsingular whenever X is full rank — and
    F3 established that before the loop began. A design scaled small enough to underflow the
    determinant loses rank under matrix_rank first, so F3 raises before this function is ever called.
    It is kept because it is the contract between the floor and the halving loop: remove the floor, or
    weaken F3 to a warning, and this branch becomes the thing that stops a crash. Do not delete it as
    dead code without deleting the reason it is dead.

    `np.logaddexp(0.0, eta)`, never `np.log(1 + np.exp(eta))`: the second overflows at eta ~ 710 and
    loses precision long before, while the clip keeps exp in range and logaddexp is exact at the ends.
    """
    eta = np.clip(Xc @ beta, -C.FIRTH_ETA_CLIP, C.FIRTH_ETA_CLIP)
    ll = float(np.sum(y * eta - np.logaddexp(0.0, eta)))
    w = np.clip(1.0 / (1.0 + np.exp(-eta)) * (1.0 - 1.0 / (1.0 + np.exp(-eta))),
                C.FIRTH_WEIGHT_FLOOR, None)
    _, logdet = np.linalg.slogdet((Xc * w[:, None]).T @ Xc)
    return ll + 0.5 * logdet if np.isfinite(logdet) else -np.inf


def _probabilities(Xc: np.ndarray, beta: np.ndarray) -> np.ndarray:
    """The fitted probabilities, on the SAME clipped linear predictor the loop used.

    Written out rather than inlined because it is called on both convergence routes and F5 asserts
    strict interiority on its output (§5.5): if this applied a different clip — or none — from the one
    inside the loop, F5 would be testing a quantity the weights are not computed from, and a fit whose
    linear predictor hit FIRTH_ETA_CLIP could pass F5 while producing a weight of exactly 0.0.
    """
    return 1.0 / (1.0 + np.exp(-np.clip(Xc @ beta, -C.FIRTH_ETA_CLIP, C.FIRTH_ETA_CLIP)))


#   THE THREE SILENT FAILURES THIS FUNCTION IS THE ONLY GUARD AGAINST [Stage 6 §3.1, §5.4]
#
#   Int64 + pd.NA  --to_numpy(float)-->  nan   SILENTLY, no warning
#   one nan in X   --firth-->  all-nan beta, all-nan p, and the loop RETURNS NORMALLY
#                              with a plausible iteration count.        [Stage 6 §3.1]
#   one nan in y   --firth-->  raises, but with the step-halving message: a sentence
#                              about [§7]'s single estimator for a missing outcome.
#   y = an ordinal --firth-->  FITS. Converges. Returns coefficients for a logistic
#                              model of a seven-level ordinal, with nothing missing
#                              and nothing to notice.                   [Stage 6 §5.4]
#   So this function is not defensive. F4 is the only thing between a missing covariate
#   and a [nan, nan] confidence interval in Stage 10; F6 is the only thing between an
#   ordinal response and a table of confident nonsense. Stages 8, 9 and 12 call firth
#   directly [§0.1, §9] — there is no second checkpoint downstream of this one.

def _assert_fittable(Xc: np.ndarray, y: np.ndarray, columns: tuple[str, ...]) -> None:
    """F3, F4 and F6 — the three things that make `firth` return, or blame, the wrong answer [§3.1].

    Ordered design → response → rank: F4 and F6 are O(n) and their failures are the common ones, and
    F3's matrix_rank is an SVD, so it goes last. Unlike D1-D4 these raise at the first failure rather
    than collecting, because each one makes the next check meaningless — a rank computed over a design
    containing nan is not a rank.

    F6's finiteness branch runs BEFORE its domain branch, and the order is load-bearing: `np.isin`
    against a nan is False, so without it a missing response would be reported as a domain error.

    The messages name neither the exposure nor any outcome, and that is a constraint rather than a
    stylistic choice. The natural message for the first branch names a specific caller's mask and for
    the second a specific ordinal column — and both would put the name of a response inside the module
    §0.1 defines as outcome-agnostic. A module that says that in a raise is one grep away from being a
    module that knows what the exposure is, and the next reader to add a check will follow the
    precedent. The concrete examples live in the spec instead, where they are more useful; §12.12
    scans this module's raise strings for exactly that.
    """
    if not np.all(np.isfinite(Xc)):
        bad = [columns[j - 1] if j else "intercept"
               for j in np.unique(np.argwhere(~np.isfinite(Xc))[:, 1])]
        raise FitError(
            f"F4  the design carries a non-finite value in: {', '.join(bad)}. An Int64 column with "
            "pd.NA converts to nan SILENTLY [§3.1], the fit then runs to completion and returns "
            "all-nan coefficients with a plausible iteration count, and the first thing anyone sees "
            "is a [nan, nan] interval in Stage 10. Use model.complete_cases() [§4.4].")

    if not np.all(np.isfinite(y)):
        raise FitError(
            f"F6  the response carries {int((~np.isfinite(y)).sum())} non-finite value(s) of "
            f"{len(y)}. Every candidate step is then nan, no step is ever accepted, and the fit "
            "raises step-halving-exhausted at iteration 1 — a message about the estimator for a "
            "failure of the data [§3.1]. The response is complete-case per estimate [§11] and the "
            "mask is the caller's: this module does not know what the response is [Stage 6 §4.4, §9].")

    off = np.unique(y[~np.isin(y, (0.0, 1.0))])
    if off.size:
        raise FitError(
            f"F6  the response takes {off.size} value(s) outside {{0, 1}}: "
            f"{', '.join(f'{v:g}' for v in off[:10])}. This is a logistic fit; a response on any "
            "other scale produces confident coefficients for a model of something else, with no "
            "missing value and nothing to notice. An ordinal response belongs to a "
            "proportional-odds fit and does not come through here [Stage 6 §5.4, §9].")

    rank = int(np.linalg.matrix_rank(Xc))
    if rank < Xc.shape[1]:
        raise FitError(
            f"F3  the design has rank {rank} of {Xc.shape[1]} columns: "
            f"{', '.join(('intercept',) + columns)}. The information matrix is singular at every "
            "beta, so no penalised optimum is identified. The constant-column rule [§4.3] removes an "
            "absent factor level; exact collinearity between two covariates is a [§6] "
            "specification error — penumbra_ml with core_ml and tmax6_ml is the one [§6] names.")


def firth(X: pd.DataFrame, y: np.ndarray) -> Fit:
    """Firth-penalised logistic regression [§7]. Raises FitError; never returns a fallback.

    `X` carries no intercept — one is prepended here, so a design matrix is never ambiguous about
    whether it has one, and `beta[0]` is always the intercept.

    Six numerical details, all lifted from `pilots/analysis.py:346-393` and all load-bearing:

    * `beta` starts at zero, so the first `p` is 0.5 everywhere and the first information matrix is
      X'X/4 — the best-conditioned it will ever be, and the iteration where a rank deficiency shows.
    * The hat diagonal is computed on the square-root-weighted design via `einsum`, not as
      `np.diag(H)`: the full H is n x n and only its diagonal is used.
    * The step is **rescaled, not rejected**, when it exceeds the trust radius. A trust region: the
      direction is a Newton direction and is kept, only its length is bounded, so a first iteration
      far from the optimum cannot leap into a region where the likelihood is -inf and burn all the
      halvings.
    * Step-halving is on the **penalised** likelihood, not the unpenalised one. Halving on l(b) would
      drive the fit toward the boundary the penalty exists to hold it back from — while converging.
    * `ll_old` is updated only when the loop continues, so the acceptance test at the next iteration
      compares against the last *accepted* value. The pilots' asymmetry here is deliberate and kept.
    * `np.linalg.inv` is deliberately NOT wrapped. `pilots/analysis.py:366-369` catches LinAlgError
      and falls back to `np.linalg.pinv` — a fallback not to a different estimator but to a different
      *estimand*: the pseudo-inverse silently picks the minimum-norm solution among infinitely many,
      so the fit returns coefficients for a design that identifies none. Measured: `inv` does raise on
      an exactly singular information matrix, so the pinv branch is what hides it.

    **The trust radius is RELATIVE to the iterate** — `‖step‖ / max(1, ‖beta‖)` — and that is a
    correction, not a detail [Stage 6 §5.2a]. A Firth estimate is equivariant under rescaling a
    covariate; an absolute bound is not, because FIRTH_MAX_STEP x FIRTH_MAX_ITER caps how far the
    coefficient vector can travel in total. Measured: a four-record design whose optimum sits at
    ‖beta‖ ~ 2008 was reported as a failure under the absolute form and converges in 10 iterations
    under this one, and the same data scaled by 1000 converged either way. If failure depends on
    covariate units, *which* replicates [§10] drops is partly an artifact of how the data was
    recorded — selection on the replicate.

    **One failure mode is not guarded, and it is stated plainly because it is not** [Stage 6 §5.2c].
    `l(b)` is concave; `½log|I(b)|` is not, and their sum is not. So this loop finds *a* stationary
    point of the penalised likelihood where [§7] prescribes *the maximiser*, and on some designs those
    differ — measured at 1 replicate in 2000 on this cohort and 4 of 600 synthetic separated designs,
    with `e` interior, the design full rank, and every counter below reading zero. There is no
    symptom. Four detectors were designed and measured to fire 0 of 1, because the basins are narrow
    and design-dependent so only probes that **vary per fit** find them — and this function fits one
    design, holds no seed, and cannot draw varying probes without either breaking determinism or
    inventing a seed from the data's own bytes. Stage 10 owns the detection: it has the replicate loop
    and C.SEED. Do not add a fixed-seed multistart here; it was tried, and it catches nothing while
    moving every pinned number.
    """
    Xc = np.column_stack([np.ones(len(y)), X.to_numpy(dtype=float)])
    y = np.asarray(y, dtype=float)
    _assert_fittable(Xc, y, tuple(X.columns))            # §5.4
    beta = np.zeros(Xc.shape[1])
    ll_old = _penalised_loglik(Xc, y, beta)
    first_step_norm, rescales, halvings = 0.0, 0, 0      # §3.2; §12.4a asserts on all three

    for iteration in range(1, C.FIRTH_MAX_ITER + 1):
        eta = np.clip(Xc @ beta, -C.FIRTH_ETA_CLIP, C.FIRTH_ETA_CLIP)
        p = 1.0 / (1.0 + np.exp(-eta))
        w = np.clip(p * (1.0 - p), C.FIRTH_WEIGHT_FLOOR, None)
        Xw = Xc * np.sqrt(w)[:, None]
        info_inv = np.linalg.inv(Xw.T @ Xw)              # LinAlgError is NOT caught — §5.4
        h = np.clip(np.einsum("ij,jk,ik->i", Xw, info_inv, Xw), 0.0, 1.0)
        score = Xc.T @ (y - p + h * (0.5 - p))
        step = info_inv @ score

        norm = float(np.linalg.norm(step))
        if iteration == 1:
            first_step_norm = norm
        # RELATIVE to the current iterate, never absolute — §5.2a. A trust region, not a rejection:
        # the Newton direction is kept and only its length is bounded. `max(1.0, ...)` and not
        # `‖beta‖`, which is unbounded at beta = 0.
        size = norm / max(1.0, float(np.linalg.norm(beta)))
        if size > C.FIRTH_MAX_STEP:
            step = step * (C.FIRTH_MAX_STEP / size)
            rescales += 1

        for _ in range(C.FIRTH_MAX_HALVINGS):
            ll_new = _penalised_loglik(Xc, y, beta + step)
            if ll_new >= ll_old:
                break
            step = step / 2.0
            halvings += 1
        else:
            raise FitError(
                f"Firth: step-halving exhausted {C.FIRTH_MAX_HALVINGS} halvings at iteration "
                f"{iteration} without increasing the penalised likelihood. [§7] prescribes one "
                "estimator; there is no second one to try.")

        beta = beta + step
        # TWO CONVERGENCE ROUTES, either one sufficient, and the Fit records which fired [§5.3]:
        #
        #     |Δ l*(b)| < FIRTH_TOL                  →  converged_on = "likelihood"
        #     max |modified score| < FIRTH_SCORE_TOL →  converged_on = "score"
        #
        # THE STEP NORM IS DELIBERATELY NOT A THIRD. Under near-separation the penalised surface is
        # genuinely flat: the fit is finite and correct, and the Newton step keeps drifting along a
        # ridge without changing the objective. Requiring ‖step‖ → 0 as well would report that fit as
        # a failure — and [§10] drops failed replicates, so the ones dropped would be exactly the
        # sparse ones, which is selection on the replicate.
        #
        # The pilots' comment attributed this rule to the R reference implementation. **That
        # attribution is struck**: its source has been read (Stage 6 §18g) and it tests the
        # likelihood change, the score AND the coefficient step CONJUNCTIVELY, so it requires the
        # very criterion excluded here and requires all three at once rather than either of two. The
        # argument above is unaffected, because it never depended on the attribution — which is
        # exactly why the spec called it a sentence to strike rather than a decision to revisit.
        #
        # Two measured facts about the score route, neither known when the criteria were written:
        #
        #   * `score` is bound at the TOP of the iteration, from the PRE-step beta, and both tests run
        #     after `beta = beta + step`. So a Fit reporting "score" reports a flat score one step
        #     behind its own coefficients. Harmless — when the score is below the tolerance the step is
        #     negligible — but stated because §12 asserts converged_on as a property of the return.
        #   * The score route fires exactly where it was designed to: on an ILL-CONDITIONED
        #     information matrix. Near the optimum Δl* ≈ ½·sᵀI⁻¹s ~ ½·s²·‖I⁻¹‖. An earlier reading
        #     assumed ‖I⁻¹‖ ~ 1 and concluded this route was dead; at ‖I⁻¹‖ = 1e6, |s| = 1e-6 gives
        #     Δl* ~ 5e-7, comfortably ABOVE FIRTH_TOL — so the likelihood is still moving when the
        #     score has gone flat. A nearly-singular information matrix and a flat penalised surface
        #     are the same fact stated twice, which is the condition [§7] chose Firth for. A replicate
        #     that comes back "score" is telling you its design is close to singular, not that
        #     something is broken.
        #
        # The ORDER the two run in is part of the estimator: [§10] refits in every replicate, so
        # reordering them changes the sampling distribution. That is a [§13] amendment, not a repair.
        moved = abs(ll_new - ll_old)                     # BEFORE ll_old is rebound — §5.2b
        if moved < C.FIRTH_TOL:
            return Fit(beta, _probabilities(Xc, beta), iteration, "likelihood", tuple(X.columns),
                       first_step_norm, rescales, halvings)
        if float(np.max(np.abs(score))) < C.FIRTH_SCORE_TOL:
            return Fit(beta, _probabilities(Xc, beta), iteration, "score", tuple(X.columns),
                       first_step_norm, rescales, halvings)
        ll_old = ll_new

    raise FitError(
        f"Firth: no convergence in {C.FIRTH_MAX_ITER} iterations. The penalised likelihood moved "
        f"{moved:.3g} on the last step against a tolerance of {C.FIRTH_TOL:g}, and the largest "
        f"modified score component was {float(np.max(np.abs(score))):.3g} against "
        f"{C.FIRTH_SCORE_TOL:g}. {rescales} step(s) were shortened by the trust region and {halvings} "
        f"halving(s) were taken [§5.2a]. [§7] prescribes one estimator, and [§10] drops and counts "
        "the replicate rather than substituting another.")


def predict(fit: Fit, X: pd.DataFrame) -> np.ndarray:
    """The fitted probabilities of `fit` on `X`. `X` carries no intercept; one is prepended.

    Named neither for treatment nor for any covariate: this is `Fit` evaluated somewhere, and the
    somewhere is the caller's business [Stage 6 §0.1]. `outcome.py` is what knows that one of the
    columns is an exposure and that setting it to a constant makes the result a counterfactual.

    **Why it exists at all**: `Fit.p` is the fitted probability at the OBSERVED treatment, because the
    `X` handed to `firth` carries the observed treatment column. It is `m_A(X)` — the counterfactual
    matching each row's actual arm — and it is therefore neither `m_1` nor `m_0`. A caller reading it
    for either is wrong on every row assigned to the other arm [Stage 9 §7.3].

    The column check is BY NAME and is not a length check. `design` drops constant columns and
    returns their names, so a frame built from a covariate list has a different width from the one a
    fit saw -- and a positional dot product against a mismatched design returns a number rather than
    an error. Comparing `fit.columns` to `X.columns` as a SEQUENCE catches a reordering too, which a
    set comparison would not: `firth`'s beta is positional after the intercept.

    **THIS IS `_probabilities` EVALUATED SOMEWHERE ELSE, AND IT IS THE SAME EXPRESSION ON PURPOSE.**
    Two details, and `test_model.py` asserts the consequence as an EQUALITY against `fit.p` rather
    than as a tolerance:

      * The intercept is PREPENDED AS A COLUMN, never added as a scalar. `firth` computes
        `eta = Xc @ beta` over an intercept-prepended `Xc`; `beta[0] + X @ beta[1:]` is the same
        mathematics in a different summation order and does not reproduce `Fit.p` -- measured,
        max|difference| 1.11e-16 on 10 of 60 rows [Stage 9 §3.3].
      * The linear predictor IS CLIPPED, at FIRTH_ETA_CLIP, because `_probabilities` clips at
        FIRTH_ETA_CLIP (model.py:320) and model.py:317-318 gives the reason: a reader of `Fit` that
        applies "a different clip -- or none -- from the one inside the loop" is reading a quantity
        the fit was not computed from. 1/(1+exp(-eta)) is better conditioned than
        exp(eta)/(1+exp(eta)) at the POSITIVE tail only; at eta = -800 it overflows and returns
        exactly 0.0 where the clipped form returns 7.124576406741285e-218 [Stage 9 §3.3]. The clip
        for this form is FIRTH_ETA_CLIP and has nothing to do with Stage 8's POLR_ETA_CLIP.
    """
    if tuple(X.columns) != fit.columns:
        raise C.SchemaError(
            f"predict was given columns {tuple(X.columns)} for a fit on {fit.columns}. "
            "The coefficient vector is positional after the intercept, so a mismatched or "
            "reordered design returns a number instead of an error [Stage 9 §7.3].")
    Xc = np.column_stack([np.ones(len(X)), X.to_numpy(dtype=float)])
    eta = np.clip(Xc @ fit.beta, -C.FIRTH_ETA_CLIP, C.FIRTH_ETA_CLIP)
    return 1.0 / (1.0 + np.exp(-eta))


# --- Stage 8a — the weighted proportional-odds fit [§8] --------------------------------------------
#
#   polr(X, y, w)  is a SECOND fitter beside `firth`, not a generalisation of it. The two share only
#   `FitError`, and they differ in the one place a reader will assume they agree: `firth` PREPENDS an
#   intercept and `polr` must not be given one, because the K cutpoints ARE the intercepts and a
#   column of ones is exactly collinear with their sum. O6 is what turns that mistake into a raise.
#
#   It is general in its covariates from the first line and that is not speculative generality:
#   [§14a] prescribes `logit P(Y <= k | A, X) = alpha_k + beta*A + gamma'X` over
#   STANDARDISATION_COVARIATES, which is this estimator with a wider design and unit weights. One
#   implementation, two callers [Stage 8 §0.1].


@dataclass(frozen=True)
class PolrFit:
    """One weighted proportional-odds fit, in the `logit P(Y <= k) = alpha_k + x'beta` form.

    `alpha` is ascending and `categories` says which response values were FITTED, so a caller can
    name the declared level each cutpoint belongs to WITHOUT re-deriving the collapse: `alpha[j]` is
    the cutpoint at or below `categories[j]` [Stage 8 §5.3].

    The last four fields describe how the fit was *reached* rather than what it is, which is `Fit`'s
    precedent and earns its place for `Fit`'s reason (Stage 6 §3.2): [§10] refits N_BOOT times and a
    safeguard whose activation nobody can count is a safeguard nobody can evaluate.

    **`categories` is a field and not a length**, and it earns its place three times: the collapse
    makes the number of cutpoints a property of the SAMPLE rather than of the declared level set; the
    audit table renders one row per *declared* level and marks the absent ones, which it cannot do
    from a count; and [§10] hands Stage 10 a `beta` whose comparability across replicates depends on
    which categories each replicate had. An `n_categories: int` would carry the shape and lose the
    identity, and the identity is what the log has to print.

    **THERE IS NO STANDARD ERROR FIELD, and that is a correctness claim rather than an omission**
    [Stage 8 §5.6]. `-H^-1` at the optimum is the inverse observed information of a likelihood in
    which `w_i` counts observations; observation weights supplied by a caller need not be frequencies,
    and where they are a tilting function of an ESTIMATED nuisance parameter the sampling variability
    of estimating it is omitted entirely. [§10]'s percentile bootstrap is the prespecified interval.
    """

    beta: np.ndarray             # (m,), one per design column — NO intercept, and no alpha
    alpha: np.ndarray            # (K,), ascending; K = len(categories) - 1
    categories: tuple[int, ...]  # the response values FITTED, ascending
    columns: tuple[str, ...]     # the design's column names, in beta's order
    iterations: int
    converged_on: str            # "likelihood" | "score"
    first_step_norm: float       # ‖step‖ at iteration 1, BEFORE any rescale   ┐
    rescales: int                # steps shortened by the trust region         ├─ Stage 8 §5.2
    halvings: int                # total step-halvings across all iterations   ┘


def _weighted_categories(y: np.ndarray, w: np.ndarray) -> tuple[int, ...]:
    """The response levels carrying POSITIVE total weight, ascending. [Stage 8 §5.3]

    One function and two callers — `_assert_polr_fittable`'s O5 counts what this returns and `polr`
    collapses to it — because the expression is not the obvious one and O5's ordering turns on a
    subtlety inside it: `np.nan > 0.0` is False, so a single nan weight makes its level's total
    weight nan and DELETES the level rather than propagating. Written twice, the two copies are one
    edit from disagreeing about which levels a fit has, and the disagreement would be silent in both
    directions: O5 would count a different number of categories from the number `polr` then fits, and
    both numbers would be finite and plausible (measured).

    `int(v)` is safe because O4 has already established that every response value is integral, and O4
    precedes O5 for that reason as well. This function does NOT re-run O3: it is called from inside
    O5, so a nan weight has already raised by the time `polr` calls it a second time.
    """
    return tuple(int(v) for v in np.unique(y) if w[y == v].sum() > 0.0)


def _ord_pieces(Xb: np.ndarray, alpha: np.ndarray, upper: np.ndarray, lower: np.ndarray,
                has_u: np.ndarray, has_l: np.ndarray) -> tuple[np.ndarray, ...]:
    """The two cumulative probabilities bracketing each observation's category, and two derivatives.

    `upper`/`lower` are cutpoint indices and `has_u`/`has_l` say whether each exists: an observation
    in the LOWEST category has no lower cutpoint (P(Y <= -1) = 0) and one in the HIGHEST has no upper
    one (P(Y <= K) = 1). Those two boundary conventions are what make the first and last categories
    contribute a single-sided derivative rather than a special case — so `gu` defaults to 1.0 and `gl`
    to 0.0, and NOT to the clipped expit of a clipped index.

    The clip on the linear predictor is the one `firth` uses and for the same reason: it keeps `exp`
    in range. Measured: the largest |linear predictor| this estimator has produced is 18.2, in the
    perfectly separated construction, so POLR_ETA_CLIP is never approached on any real input —
    recorded rather than relied on.
    """
    zu = np.where(has_u, alpha[np.clip(upper, 0, len(alpha) - 1)] + Xb, 0.0)
    zl = np.where(has_l, alpha[np.clip(lower, 0, len(alpha) - 1)] + Xb, 0.0)
    gu = np.where(has_u, 1.0 / (1.0 + np.exp(-np.clip(zu, -C.POLR_ETA_CLIP, C.POLR_ETA_CLIP))), 1.0)
    gl = np.where(has_l, 1.0 / (1.0 + np.exp(-np.clip(zl, -C.POLR_ETA_CLIP, C.POLR_ETA_CLIP))), 0.0)
    du = np.where(has_u, gu * (1.0 - gu), 0.0)              # d gamma_u / d z
    dl = np.where(has_l, gl * (1.0 - gl), 0.0)
    ddu = np.where(has_u, du * (1.0 - 2.0 * gu), 0.0)       # d2 gamma_u / d z2
    ddl = np.where(has_l, dl * (1.0 - 2.0 * gl), 0.0)
    return gu, gl, du, dl, ddu, ddl


def _ord_loglik(Xn: np.ndarray, y_idx: np.ndarray, w: np.ndarray, alpha: np.ndarray,
                beta: np.ndarray, K: int) -> float:
    """The weighted log-likelihood, or -inf where any category probability is non-positive.

    -inf rather than a raise, so that step-halving can REJECT a step and try a shorter one; a raise
    there would turn a recoverable step into a failed fit. This is `_penalised_loglik`'s posture
    (Stage 6 §5.2), and the branch is what makes a crossed pair of cutpoints unreachable rather than
    caught: a crossing makes some observation's category probability non-positive, so the halving
    rejects the crossing rather than an assertion detecting it after the step was taken.

    THE BRANCH IS MEASURED NOT TO FIRE on any input this estimator has been given: 3000 attempts at
    reaching a crossing with pure Newton from a crowded start never lost the ordering. It is kept
    because it is the CONTRACT between the collapse and the halving loop — remove the collapse, or
    start from an unordered alpha, and this branch becomes the thing that stops a crash. Do not
    delete it as dead code without deleting the reason it is dead.
    """
    upper, lower = y_idx, y_idx - 1
    has_u, has_l = y_idx <= K - 1, y_idx >= 1
    gu, gl = _ord_pieces(Xn @ beta, alpha, upper, lower, has_u, has_l)[:2]
    p = gu - gl
    if np.any(p <= 0.0) or not np.all(np.isfinite(p)):
        return -np.inf
    return float(np.sum(w * np.log(p)))


def _ord_score_hess(Xn: np.ndarray, y_idx: np.ndarray, w: np.ndarray, alpha: np.ndarray,
                    beta: np.ndarray, K: int) -> tuple[np.ndarray, np.ndarray]:
    """The analytic weighted score and Hessian over (alpha, beta), in that order.

    Analytic and not differenced: a finite-difference gradient is 2*(K+m) extra likelihood
    evaluations per Newton step, and [§10] refits N_BOOT times. Verified against central differences
    at three parameter values on two frames, with the Hessian negative definite at every one of them
    — which is the concavity claim as a measurement rather than a citation.

    The derivation, once, because a wrong sign here is a plausible number: with u = P(Y <= j),
    v = P(Y <= j-1) and p = u - v, the per-observation log-likelihood is log p and

        d/dzu  = u'/p = du/p          d2/dzu2      = ddu/p - (du/p)^2
        d/dzl  = -v'/p = -dl/p        d2/dzl2      = -ddl/p - (dl/p)^2
                                      d2/dzu dzl   = du*dl/p^2

    and zu = alpha_ju + x'beta, zl = alpha_jl + x'beta, so d z_k / d alpha_m = 1{k = m} and
    d z_k / d beta = x. The beta block therefore collects Huu + 2*Hul + Hll, which is the sum of
    both single-sided terms PLUS TWICE the cross term, and the alpha-beta blocks collect one
    single-sided term plus the cross term each.
    """
    n, m = Xn.shape
    upper, lower = y_idx, y_idx - 1
    has_u, has_l = y_idx <= K - 1, y_idx >= 1
    gu, gl, du, dl, ddu, ddl = _ord_pieces(Xn @ beta, alpha, upper, lower, has_u, has_l)
    p = gu - gl

    A, B = du / p, -dl / p
    Huu = ddu / p - (du / p) ** 2
    Hll = -ddl / p - (dl / p) ** 2
    Hul = du * dl / p ** 2

    g = np.zeros(K + m)
    H = np.zeros((K + m, K + m))
    np.add.at(g, upper[has_u], (w * A)[has_u])
    np.add.at(g, lower[has_l], (w * B)[has_l])
    g[K:] = Xn.T @ (w * (A + B))

    np.add.at(H, (upper[has_u], upper[has_u]), (w * Huu)[has_u])
    np.add.at(H, (lower[has_l], lower[has_l]), (w * Hll)[has_l])
    both = has_u & has_l
    np.add.at(H, (upper[both], lower[both]), (w * Hul)[both])
    np.add.at(H, (lower[both], upper[both]), (w * Hul)[both])

    cu = np.where(both, Huu + Hul, np.where(has_u, Huu, 0.0))
    cl = np.where(both, Hll + Hul, np.where(has_l, Hll, 0.0))
    for r in range(m):
        column = np.zeros(K)
        np.add.at(column, upper[has_u], (w * cu * Xn[:, r])[has_u])
        np.add.at(column, lower[has_l], (w * cl * Xn[:, r])[has_l])
        H[:K, K + r] += column
        H[K + r, :K] += column
    H[K:, K:] = (Xn * (w * (Huu + 2.0 * Hul + Hll))[:, None]).T @ Xn
    return g, H


#   THE TWO SILENT FAILURES THIS FUNCTION IS THE ONLY GUARD AGAINST [Stage 8 §5.5, §3.1]
#
#   one nan weight   --polr-->  its response level's total weight is nan, `nan > 0.0` is
#                               False, and the POSITIVE-WEIGHT RULE DELETES THE LEVEL. It
#                               does not propagate. Measured: a four-category fit became a
#                               three-category one, beta moved from -1.775 to -1.205,
#                               iterations 4, converged_on "likelihood", nothing raised.
#   one nan response --polr-->  np.unique returns nan as a level, `w[y == nan]` is empty so
#                               its total weight is 0.0, and the RECORD IS DROPPED by the
#                               same rule. Measured: 80 records fitted as 79, with nothing
#                               missing anywhere in the output.
#   So O3 MUST PRECEDE O5. With O3 first the nan is a raise; with O5 first it is a different
#   model, over a coarser scale, reported with a plausible iteration count.

def _assert_polr_fittable(Xn: np.ndarray, y: np.ndarray, w: np.ndarray,
                          columns: tuple[str, ...]) -> None:
    """O1-O6 — the six things that make `polr` return, or blame, the wrong answer.

    Ordered cheap-to-expensive and raising at the FIRST failure rather than collecting, as
    `_assert_fittable` does (Stage 6 §5.4) and for its reason: each one makes the next meaningless.
    A rank computed over a design containing nan is not a rank.

    **O4 MUST precede O5** for the same shape of reason O3 does, and `int(v)` inside
    `_weighted_categories` is safe only because it does.

    The messages name neither the exposure nor any outcome, which is a constraint and not a stylistic
    choice: this module is outcome-agnostic (Stage 6 §0.1) and the acceptance suite scans its raise
    strings for both. The concrete examples live in the spec instead.
    """
    if not np.all(np.isfinite(Xn)):
        bad = [columns[j] for j in np.unique(np.argwhere(~np.isfinite(Xn))[:, 1])]
        raise FitError(
            f"O1  the design carries a non-finite value in: {', '.join(bad)}. An Int64 column with "
            "pd.NA converts to nan SILENTLY, and every candidate step is then nan, so no step is "
            "accepted and the fit blames step-halving for a failure of the data. Use "
            "model.complete_cases().")

    if not np.all(np.isfinite(y)):
        raise FitError(
            f"O2  the response carries {int((~np.isfinite(y)).sum())} non-finite value(s) of "
            f"{len(y)}. It does NOT propagate here — np.unique returns nan as a level, its total "
            "weight is 0, and the record is DROPPED by the positive-weight rule with nothing "
            "missing anywhere in the result. Measured: 80 records fitted as 79. The response is "
            "complete-case per estimate and the mask is the caller's.")

    # O3 BEFORE O5, and the order is the finding rather than a preference: `np.nan > 0.0` is False,
    # so without this a nan weight DELETES the response category it sits in instead of raising, and
    # the fit comes back finite, plausible, and over a coarser scale (measured, above).
    if not np.all(np.isfinite(w)) or np.any(w < 0.0):
        raise FitError(
            f"O3  the weights carry {int((~np.isfinite(w)).sum())} non-finite and "
            f"{int(np.sum(w < 0.0))} negative value(s). A non-finite weight DELETES an outcome "
            "category rather than propagating (see O5's ordering); a negative one makes the "
            "objective non-concave, so the Hessian may be indefinite and a Newton step may ascend "
            "away from the maximiser while every counter reads clean.")

    off = y[np.not_equal(np.mod(y, 1.0), 0.0)]
    if off.size:
        raise FitError(
            f"O4  the response takes {off.size} non-integer value(s): "
            f"{', '.join(f'{v:g}' for v in np.unique(off)[:10])}. This is an ORDINAL fit and the "
            "response is a category code: one cutpoint is estimated per distinct value, so a "
            "continuous response silently makes the parameter count the number of distinct values "
            "and every category a singleton. It fits. A continuous response belongs to no model in "
            "this pipeline.")

    categories = _weighted_categories(y, w)                                # one definition
    if len(categories) < 2:
        raise FitError(
            f"O5  {len(categories)} response category/ies carry positive weight: "
            f"{tuple(float(v) for v in categories)}. A proportional-odds fit needs at least two: "
            "with one, there is no cutpoint to estimate and the likelihood is constant in every "
            "parameter.")

    rank = int(np.linalg.matrix_rank(Xn)) if Xn.shape[1] else 0
    if Xn.shape[1] and rank < Xn.shape[1]:
        raise FitError(
            f"O6  the design has rank {rank} of {Xn.shape[1]} columns: {', '.join(columns)}. The "
            "cutpoints ARE the intercepts, so a column of ones is exactly collinear with their sum "
            "and a design carrying its own intercept fails here — which is the opposite of the "
            "penalised logistic fit in this module, and the one place a reader will assume the two "
            "agree.")


def polr(X: pd.DataFrame, y: np.ndarray, w: np.ndarray | None = None) -> PolrFit:
    """Weighted proportional-odds regression [§8]. Raises FitError; never returns a fallback.

    THE PARAMETRISATION, spelled out, because every reference implementation uses the other sign::

        logit P(Y <= k | x)  =  alpha_k  +  x'beta ,        k = 0 … K-1
        alpha_1 < alpha_2 < … < alpha_K                     the cutpoints, ascending

        l(theta)  =  SUM_i  w_i * log( P(Y = y_i | x_i) )   the log-likelihood is WEIGHTED
                     with P(Y = j) = P(Y <= j) - P(Y <= j-1)

    This is [§14a]'s own form, so it is the SAP's parametrisation and not a choice made here.
    `statsmodels`' `OrderedModel` and R's `polr` and `clm` all write `theta_k - x'beta` instead, and
    measured against both: **the coefficients come back NEGATED and the cutpoints do NOT.** That
    asymmetry is what makes the mistake survive review — a reader comparing two fits sees four
    cutpoints agreeing to five decimals and concludes the implementations agree. Do not "fix" the sign
    to match a textbook; a caller reading this coefficient under the other convention inverts the
    direction of its own result with every number still finite and still plausible.

    `X` carries NO intercept and must not: the K cutpoints are the intercepts, and adding a column of
    ones makes the design exactly singular against their sum. That is the opposite of `firth`, which
    prepends its own — so the two functions differ in the one place a reader will assume they agree,
    and O6's rank check is what turns the mistake into a raise rather than into a pseudo-inverse.

    `w` is a non-negative observation-weight vector or None. None is the UNWEIGHTED fit and is not a
    synonym for ones: a caller who passes nothing has said something different from a caller who
    passes ones, and the acceptance suite turns on the two being indistinguishable in the result and
    visible in the call.

    Four numerical details, each load-bearing:

    * **The response is collapsed to the categories carrying positive weight**, before anything else
      numeric happens, and `categories` reports what was fitted.
    * **The start values are beta = 0 and alpha_k = logit of the weighted cumulative share.** At
      beta = 0 that is the exact maximiser of the intercept-only model, so iteration 1 starts at the
      solution of a nested model rather than at an arbitrary point — measured 4 to 5 iterations on
      every frame. And they are FINITE BY CONSTRUCTION because of the collapse: every kept category
      carries positive weight, so every cumulative share is strictly inside (0, 1).
    * **The step is rescaled, not rejected**, when it exceeds the trust radius, and the radius is
      RELATIVE to the iterate — `firth`'s correction (Stage 6 §5.2a), inherited for its reason: a
      proportional-odds estimate is equivariant under rescaling a covariate and an absolute bound is
      not, so with one, WHICH replicates [§10] drops would depend on how the data was recorded.
    * **`np.linalg.solve`, not `inv`, and it is deliberately NOT wrapped.** Stage 6 §5.2's argument
      applies unchanged: a `pinv` fallback silently picks the minimum-norm solution among infinitely
      many, so the fit returns coefficients for a design that identifies none. O6 establishes full
      rank before the loop begins.

    Unlike `firth` there is no penalty and no ½log|I| term, so **the objective is concave** — the
    weighted sum of concave per-observation log-likelihoods, verified negative definite at the optimum
    and away from it. Stage 6 §5.2c's one unguarded failure mode, "this loop finds *a* stationary
    point where [§7] prescribes *the* maximiser", therefore does not recur here: for this estimator
    they are the same point.

    **What is far more dangerous instead is that A SEPARATED FIT CONVERGES AND RETURNS.** Measured on
    20 against 20 with no overlap: beta 36.4058, exp(beta) 6.5e15, 17 iterations on the score
    criterion, zero rescales, zero halvings, every fitted probability finite and nothing anywhere out
    of range. The likelihood flattens as beta grows and the loop returns normally, so there is no
    failure for [§10] to drop and count. This function does not guard it — a bound on a coefficient is
    a statement about the specification and not about the arithmetic, so it lives on the caller's
    side, exactly as Stage 6 put F5 in `propensity.py` and not here.
    """
    Xn = np.asarray(X.to_numpy(dtype=float))
    y = np.asarray(y, dtype=float)
    w = np.ones(len(y)) if w is None else np.asarray(w, dtype=float)
    _assert_polr_fittable(Xn, y, w, tuple(X.columns))                      # O1-O6

    # POSITIVELY-WEIGHTED and not merely OBSERVED, and the two reasons are one decision seen twice:
    # it is what makes every cumulative share below strictly interior, so the start values are finite
    # by construction; and it is what makes a crossed pair of cutpoints unreachable, because a
    # crossing only drives a category probability non-positive FOR OBSERVATIONS IN THAT CATEGORY, and
    # if they all carry zero weight the `-inf` never fires. The rule is one function with two callers
    # — O5 counts what this returns — because the nan subtlety above lives inside the expression.
    categories = _weighted_categories(y, w)
    keep = np.isin(y, np.asarray(categories, dtype=float))
    Xn, y, w = Xn[keep], y[keep], w[keep]
    K, m = len(categories) - 1, Xn.shape[1]
    y_idx = np.searchsorted(np.asarray(categories, dtype=float), y)

    share = np.array([w[y_idx <= k].sum() / w.sum() for k in range(K)])    # strictly interior
    par = np.concatenate([np.log(share / (1.0 - share)), np.zeros(m)])
    ll_old = _ord_loglik(Xn, y_idx, w, par[:K], par[K:], K)
    first_step_norm, rescales, halvings = 0.0, 0, 0

    for iteration in range(1, C.POLR_MAX_ITER + 1):
        g, H = _ord_score_hess(Xn, y_idx, w, par[:K], par[K:], K)
        step = -np.linalg.solve(H, g)                       # LinAlgError is NOT caught
        norm = float(np.linalg.norm(step))
        if iteration == 1:
            first_step_norm = norm
        size = norm / max(1.0, float(np.linalg.norm(par)))
        if size > C.POLR_MAX_STEP:
            step = step * (C.POLR_MAX_STEP / size)
            rescales += 1

        for _ in range(C.POLR_MAX_HALVINGS):
            ll_new = _ord_loglik(Xn, y_idx, w, (par + step)[:K], (par + step)[K:], K)
            if ll_new >= ll_old:
                break
            step = step / 2.0
            halvings += 1
        else:
            raise FitError(
                f"polr: step-halving exhausted {C.POLR_MAX_HALVINGS} halvings at iteration "
                f"{iteration} without increasing the weighted log-likelihood. [§8] prescribes one "
                "estimator; there is no second one to try [invariant 5].")

        par = par + step
        # TWO CONVERGENCE ROUTES, either one sufficient, and the PolrFit records which fired —
        # `firth`'s rule (Stage 6 §5.3), with the step norm deliberately not a third for its reason.
        # The ORDER they run in is part of the estimator: [§10] refits in every replicate, so
        # reordering them changes the sampling distribution. That is a [§13] amendment, not a repair.
        #
        # `g` is bound at the TOP of the iteration and both tests run after `par = par + step`, so a
        # PolrFit reporting "score" reports a flat score one iterate behind its own coefficients.
        # Harmless — when the score is below the tolerance the step is negligible — but stated
        # because the acceptance suite asserts `converged_on` as a property of the return, and an
        # implementer who recomputed the score after the step to "fix" this would change which
        # iteration the loop exits on.
        #
        # AND THE SCORE ROUTE IS THE ONE A SEPARATED FIT RETURNS ON. Stage 6 §5.3 established that a
        # fit coming back "score" is telling you its design is close to singular; here it is telling
        # you something stronger and worse.
        moved = abs(ll_new - ll_old)
        if moved < C.POLR_TOL:
            return PolrFit(par[K:], par[:K], categories, tuple(X.columns), iteration,
                           "likelihood", first_step_norm, rescales, halvings)
        if float(np.max(np.abs(g))) < C.POLR_SCORE_TOL:
            return PolrFit(par[K:], par[:K], categories, tuple(X.columns), iteration,
                           "score", first_step_norm, rescales, halvings)
        ll_old = ll_new

    raise FitError(
        f"polr: no convergence in {C.POLR_MAX_ITER} iterations. The weighted log-likelihood moved "
        f"{moved:.3g} against a tolerance of {C.POLR_TOL:g}, and the largest score component was "
        f"{float(np.max(np.abs(g))):.3g} against {C.POLR_SCORE_TOL:g}. {rescales} step(s) were "
        f"shortened by the trust region and {halvings} halving(s) taken. [§8] prescribes one "
        "estimator, and [§10] drops and counts the replicate rather than substituting another.")


def ordinal_probabilities(fit: PolrFit, X: pd.DataFrame) -> np.ndarray:
    """P(Y = c | x) for every c in `fit.categories`: (len(X), len(fit.categories)). [Stage 12 §6]

    The counterpart of `predict` for the ordinal fitter, and it inherits three of its four rules:
    the column check is BY NAME against `fit.columns` and a reorder raises `C.SchemaError`; the
    linear predictor is clipped at POLR_ETA_CLIP; and the function names neither treatment nor
    covariate, so a caller supplies the counterfactual meaning by overwriting a column
    [Stage 12 §7.1].

    The fourth rule is where it differs and the difference is load-bearing. `predict` prepends its
    own intercept because `firth` fits one; this must not, because THE CUTPOINTS ARE THE INTERCEPTS
    and O6 has already established that `X` carries none::

        P(Y <= k | x) = expit(alpha_k + x'beta),   with P(Y <= -1) = 0 and P(Y <= K) = 1 EXACTLY
        P(Y = c_j | x) = P(Y <= j) - P(Y <= j-1)

    The two structural bounds are written as literal 0.0 and 1.0 columns and never as an expit of a
    large number: `expit(POLR_ETA_CLIP)` is 1.0 in float64 but that is an accident of the clip, and
    a row's probabilities must sum to exactly 1.0 by CONSTRUCTION for Stage 12 §20.5's assertion to
    be about the arithmetic rather than about the clip.

    **Nothing else in this pipeline computes this, and the gap is exact** [Stage 12 §6.1]. `predict`
    evaluates a Firth binary `Fit` and returns one probability per row. `_ord_pieces` computes the
    two cumulative probabilities bracketing ONE category per row -- the row's own observed one -- and
    is private and used only inside the log-likelihood, the score and the Hessian. Neither is
    `P(Y = j | x)` for every `j`, which is what g-computation averages.

    **THE RETURNED LEVEL SET IS THE FITTED ONE AND NOT THE DECLARED ONE.** `polr` collapses the
    response to the categories carrying positive weight before it fits anything, so on a frame
    missing an mRS level this returns six columns where `C.MRS_LEVELS` declares seven -- measured at
    4.15% of Stage 12's replicates, always mRS 5 [Stage 12 §6.2]. Re-expressing the result on a
    declared level set, with an unoccupied level as a structural zero, is the CALLER's
    (Stage 12 §6.3): this function returns what the fit is about, and a level for which the fit
    provides no cutpoint has no probability this function could invent.
    """
    if tuple(X.columns) != fit.columns:
        raise C.SchemaError(
            f"T5  ordinal_probabilities was given columns {tuple(X.columns)} for a fit on "
            f"{fit.columns}. The coefficient vector is positional, so a mismatched or reordered "
            "design returns a number instead of an error [Stage 9 §7.3, Stage 12 §6.1]. The "
            f"difference is: {tuple(c for c in X.columns if c not in fit.columns) or 'order only'}.")
    eta = np.clip(X.to_numpy(dtype=float) @ fit.beta, -C.POLR_ETA_CLIP, C.POLR_ETA_CLIP)
    # The cumulative sequence, bracketed by EXACT 0.0 and EXACT 1.0. `alpha` is ascending and the
    # clip is monotone, so the differences below are non-negative by construction.
    z = np.clip(fit.alpha[None, :] + eta[:, None], -C.POLR_ETA_CLIP, C.POLR_ETA_CLIP)
    cumulative = np.column_stack([
        np.zeros(len(X)),                        # P(Y <= -1) = 0, EXACTLY -- never an expit
        1.0 / (1.0 + np.exp(-z)),
        np.ones(len(X)),                         # P(Y <= K)  = 1, EXACTLY -- never an expit
    ])
    return np.diff(cumulative, axis=1)


# --- Stage 12 — the random-centre-intercept proportional-odds fit [§14a] ----------------------------
#
#   polr_ri(X, y, groups)  is a THIRD fitter beside `firth` and `polr`, and it is the only estimator
#   Stage 12 adds. [§14a]'s sensitivity 2 prescribes
#
#       logit P(Y <= k | x, b_c)  =  alpha_k + x'beta + b_c ,      b_c ~ N(0, sigma^2)
#
#   with ONE common treatment effect -- no A x centre interaction and no random slope, both forbidden
#   by [§14a] and [§15] at four centres.
#
#   NOTHING IN THE ENVIRONMENT CAN FIT IT [Stage 12 §12]. `statsmodels` has `MixedLM` for the linear
#   case and Bayesian mixed GLMs for binomial and Poisson, and NO ordinal mixed model at any API;
#   `scipy` is test-only by policy [Stage 6 §2]. So it is written here, and the two things that make
#   it affordable and correct are both MEASUREMENTS rather than conventions:
#
#     * THE QUADRATURE IS ADAPTIVE.  Non-adaptive Gauss-Hermite is off by 0.48 log-likelihood units
#       at 31 nodes at sigma = 3 and is NON-MONOTONE in the node count; adaptive quadrature is off by
#       8.9e-4 at THREE nodes [Stage 12 §12.2]. 24% of replicates fit sigma_hat above 1, squarely in
#       the region where the non-adaptive rule is not converged, so this is not a refinement -- it is
#       what makes the arm computable at all.
#     * THE GRADIENT IS ANALYTIC.  The same fit costs 39.0 s under Newton with a central-difference
#       Hessian and 0.76 s under BFGS with an analytic gradient -- 21.7 hours against 17.5 minutes at
#       N_BOOT [Stage 12 §12.4, §12.8]. The gradient is part of the specification, not an
#       optimisation, and test_model.py asserts it against central differences.
#
#   THE ORDINAL DERIVATIVES ARE `_ord_pieces`' AND ARE NOT WRITTEN A SECOND TIME. `b` enters the
#   linear predictor exactly as an intercept shift does, so d/db of an observation's log-probability
#   IS the `A + B` that `_ord_score_hess` already assembles for d/dbeta, and d2/db2 IS its
#   `Huu + 2*Hul + Hll`. A second derivation here would be a second definition of the same
#   arithmetic, verified against nothing.


@dataclass(frozen=True)
class RIFit:
    """A random-intercept proportional-odds fit. No standard error, for `PolrFit`'s reason.

    **`nodes` is a field and `PolrFit` has no counterpart, because the node count CHANGES THE ANSWER**
    [Stage 12 §12.3] while `POLR_TOL` changes only how precisely the same answer is found. Measured:
    non-adaptive quadrature returns sigma_hat 0.720 / 0.803 / 0.832 at 5 / 7 / 9 nodes against
    0.542426850 at 31. A record that omitted it would let two fits with different numerical content
    compare equal on every field.

    `b` and `b_sd` are the conditional modes and the curvature at them -- the empirical-Bayes
    intercepts [§14a]'s standardisation is conditioned on (Stage 12 §12.6). They are FREE: the
    adaptive quadrature computes them at every objective evaluation, so carrying them costs nothing
    and recomputing them in the caller would be a second mode search.

    `at_floor` records that `sigma` reached `POLR_RI_SIGMA_FLOOR`, and **such a fit is NOT a failure**:
    [§14a] names sigma^2_C = 0 as a legitimate answer -- *"if outcomes truly do not differ by centre
    given X then sigma^2_C = 0 and it collapses to the pooled model"* -- and dropping those replicates
    would select the bootstrap on the value of the very parameter the arm exists to examine
    [Stage 12 §12.5].
    """

    beta: np.ndarray             # (m,), one per design column — NO intercept
    alpha: np.ndarray            # (K,), ascending
    sigma: float                 # the between-group SD, > 0
    groups: tuple[str, ...]      # the group labels, in the order `b` is keyed
    b: np.ndarray                # (G,), the conditional modes — Stage 12 §12.6
    b_sd: np.ndarray             # (G,), the curvature at each mode
    categories: tuple[int, ...]
    columns: tuple[str, ...]
    nodes: int                   # the quadrature node count this fit used — Stage 12 §12.3
    iterations: int
    converged_on: str            # "likelihood" | "score"
    first_step_norm: float
    rescales: int
    halvings: int
    at_floor: bool               # sigma reached POLR_RI_SIGMA_FLOOR — Stage 12 §12.5


def _ord_b_derivatives(Xb: np.ndarray, alpha: np.ndarray, y_idx: np.ndarray,
                       K: int) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Per-observation P(Y = y_i), and the pieces of its first two derivatives in an INTERCEPT SHIFT.

    `Xb` is the linear predictor already carrying whatever shift is being differentiated at, and may
    be any shape broadcastable against `y_idx[..., None]`-style indexing — the adaptive quadrature
    evaluates an (n, Q) grid of one shift per observation per node in one call.

    THE FOUR RETURNS ARE `_ord_score_hess`' OWN QUANTITIES AND ARE NOT A SECOND DERIVATION. `b`
    enters `eta` exactly as a coefficient on a column of ones does, so::

        d log p / d b    =  A + B                 which is d/dbeta's multiplier
        d2 log p / d b2  =  Huu + 2*Hul + Hll     which is the beta block's multiplier

    `A` and `B` are returned SEPARATELY and not pre-summed, because the cutpoint gradient needs them
    apart: `A` accumulates at the observation's UPPER cutpoint index and `B` at its lower, which is
    that function's two `np.add.at` calls. Everything else needs only the sum. Both are lifted from
    its derivation verbatim, and its docstring is where the algebra is written out once.
    """
    upper, lower = y_idx, y_idx - 1
    has_u, has_l = y_idx <= K - 1, y_idx >= 1
    gu, gl, du, dl, ddu, ddl = _ord_pieces(Xb, alpha, upper, lower, has_u, has_l)
    p = np.maximum(gu - gl, _RI_P_FLOOR)          # see _RI_P_FLOOR -- a far node, not a bad step
    A, B = du / p, -dl / p
    Huu = ddu / p - (du / p) ** 2
    Hll = -ddl / p - (dl / p) ** 2
    Hul = du * dl / p ** 2
    return p, A, B, Huu + 2.0 * Hul + Hll


# The conditional-mode search's own two numbers. They are NOT config constants and the distinction is
# the one Stage 12 §12.3 draws: `POLR_RI_NODES` changes the ANSWER, so it is prespecified; these
# govern an inner solve whose only job is to place the quadrature nodes, and their accuracy is
# verified END TO END by §20.12's gradient test against central differences on the true objective —
# measured 1e-8 relative at every sigma the bootstrap reaches. A looser inner solve shows up there as
# a worse gradient, which is a test failure rather than a silently different estimate.
_RI_MODE_TOL: Final[float] = 1e-12
_RI_MODE_ITER: Final[int] = 100
_SQRT2: Final[float] = float(np.sqrt(2.0))

# A FLOOR ON THE PER-OBSERVATION CATEGORY PROBABILITY, and it is a correctness fix rather than a
# guard. At a quadrature node the posterior has moved far away from, `gu` and `gl` are both 1.0 in
# float64 and `p = gu - gl` is EXACTLY ZERO by cancellation -- so `log p` is -inf and `du/p`,
# `(du/p)**2` and `du*dl/p**2` are all `nan`. Measured: those nans reach the gradient, the trial step
# is rejected by the `isfinite` test in the halving loop, and the fit converges to the right answer
# through a line search that is doing the wrong job. A node where the integrand underflows is
# NEGLIGIBLE, not invalid, and a quadrature rule's whole business is to weight it at zero.
#
# `polr`'s convention is the opposite -- `_ord_loglik` returns -inf and lets step-halving reject the
# step -- and it is right THERE, where a non-positive probability means a crossed pair of cutpoints
# at the iterate. Here it means a far node, which is not a property of the iterate at all.
#
# THE EXPONENT IS CHOSEN AGAINST `p**2` AND NOT AGAINST `p`. `Hul = du*dl/p**2`, so a floor below
# about 1.5e-154 makes `p**2` itself underflow to zero and reintroduces the 0/0 this fixes. 1e-150
# leaves `p**2` at 1e-300, comfortably normal, while `log(1e-150)` is -345 -- so a floored node sits
# 345 log-units below the peak and the log-sum-exp weights it at exp(-345), which is zero to every
# decimal anyone will read.
_RI_P_FLOOR: Final[float] = 1e-150


def _ri_modes(Xb: np.ndarray, alpha: np.ndarray, y_idx: np.ndarray, K: int,
              gidx: np.ndarray, n_groups: int, sigma: float) -> tuple[np.ndarray, np.ndarray]:
    """Each group's conditional mode `b_hat_c` and posterior SD `tau_c`. Newton, from b = 0.

    This is what makes the quadrature ADAPTIVE (Stage 12 §12.2): the nodes go at
    `b_hat_c + sqrt(2)*tau_c*t_q` rather than at `sigma*sqrt(2)*t_q`, which is where the mass
    actually is once the data pull a centre's posterior away from the prior.

    The objective per group is concave — `d2/db2` of the ordinal log-likelihood is negative and the
    prior contributes a further `-1/sigma^2` — so Newton ascends from any start and `b = 0` is the
    prior mode. **`1/sigma^2` is why `POLR_RI_SIGMA_FLOOR` exists**: under a `log sigma`
    parametrisation the boundary is at -inf, where this expression overflows, and Stage 12 §12.5
    observed exactly that as an overflow in the curvature and a `-inf` in a quadrature weight.
    """
    b = np.zeros(n_groups)
    prior_curvature = 1.0 / sigma ** 2
    for _ in range(_RI_MODE_ITER):
        _, A, B, d2 = _ord_b_derivatives(Xb + b[gidx], alpha, y_idx, K)
        g1 = np.zeros(n_groups)
        np.add.at(g1, gidx, A + B)
        g2 = np.zeros(n_groups)
        np.add.at(g2, gidx, d2)
        step = -(g1 - b * prior_curvature) / (_curvature(g2, prior_curvature))
        b = b + step
        if float(np.max(np.abs(step))) < _RI_MODE_TOL:
            break
    _, _, _, d2 = _ord_b_derivatives(Xb + b[gidx], alpha, y_idx, K)
    g2 = np.zeros(n_groups)
    np.add.at(g2, gidx, d2)
    return b, np.sqrt(-1.0 / _curvature(g2, prior_curvature))


def _curvature(g2: np.ndarray, prior_curvature: float) -> np.ndarray:
    """The per-group log-posterior curvature, held to the sign the mathematics guarantees.

    **THIS CLAMP ENFORCES A PROPERTY RATHER THAN HIDING A VIOLATION OF ONE, and the distinction is
    the whole reason it is written this way.** The per-group log-posterior is the sum of concave
    ordinal log-likelihood terms and a Gaussian log-prior, so `d2/db2` is negative EXACTLY, and the
    prior alone contributes `-1/sigma^2`. The total therefore cannot exceed `-1/sigma^2`.

    In float64 it sometimes does. `Huu = ddu/p - (du/p)**2` is a difference of two large quantities in
    the tails, and cancellation can return a small POSITIVE value for a term that is mathematically
    negative. When the sum crosses zero, `-1/(g2 - 1/sigma^2)` is negative, `sqrt` returns `nan`, the
    quadrature nodes are `nan`, the objective is `nan`, and the line search rejects the step — so the
    fit *still returns the right answer*, through a halving loop doing the wrong job.
    **Measured: 2 of 2000 replicates exhausted the halvings that way and raised T2**, which is a
    replicate dropped for a floating-point artefact rather than for anything about the data.

    So the ordinal contribution is clamped at its own guaranteed sign — `min(g2, 0)` — and the prior's
    is exact. The result is at most `-1/sigma^2 < 0` and the bound is TIGHT: it is attained whenever
    the data contribute no curvature, which is precisely the far-tail case that produces the noise.
    No tolerance is invented and no magnitude is chosen.
    """
    return np.minimum(g2, 0.0) - prior_curvature


def _ri_objective(par: np.ndarray, Xn: np.ndarray, y_idx: np.ndarray, K: int, m: int,
                  gidx: np.ndarray, n_groups: int, t: np.ndarray, logw: np.ndarray,
                  adaptive: bool) -> tuple[float, np.ndarray, np.ndarray, np.ndarray]:
    """The marginal log-likelihood at `par`, its analytic gradient, and the modes it placed nodes at.

    Returns `(log_likelihood, gradient, b_hat, tau)`. `par` is `(alpha, beta, log sigma)` and
    **`log sigma` and not `sigma`, so the positivity constraint is structural** (Stage 12 §12.1).

    THE INTEGRAL, and both quadrature rules, because the difference between them is 5 orders of
    magnitude and the specification is the second one::

        L_c  =  INTEGRAL  phi(b; 0, sigma^2) * PROD_i P(Y_i = y_i | x_i, b)  db

        NON-ADAPTIVE   b_q = sigma*sqrt(2)*t_q
                       log L_c = -log(pi)/2 + LSE_q[ log w_q + s_cq ]
        ADAPTIVE       b_q = b_hat_c + sqrt(2)*tau_c*t_q
                       log L_c = log(sqrt(2)*tau_c)
                                 + LSE_q[ log w_q + t_q^2 + log phi(b_q) + s_cq ]

    with `s_cq = SUM_{i in c} log P(Y_i | x_i, b_q)` and `t_q, w_q` the Gauss-Hermite nodes and
    weights for the weight `exp(-t^2)` — `numpy.polynomial.hermite.hermgauss`, which is why no
    `scipy` dependency is added (Stage 12 §2).

    **LSE is log-sum-exp and it is not a micro-optimisation.** `s_cq` is a sum of up to 41 log
    probabilities and the exponentials underflow to exactly zero at nodes the posterior has moved
    away from; a naive `log(sum(exp(...)))` returns `-inf` for a group whose every node underflows,
    and the fit then reports a finite objective for a model it never evaluated.

    **THE GRADIENT IS EXACT FOR THE FROZEN-NODE QUADRATURE SUM AND NEGLECTS THE DERIVATIVE OF THE
    NODE POSITIONS** (Stage 12 §12.4, departure 2). `b_hat_c` and `tau_c` are recomputed at the
    current parameters and then held fixed for that evaluation's gradient. Measured against central
    differences on the TRUE objective: relative error 5.0e-8 / 9.3e-9 / 9.5e-9 / 6.0e-8 at
    sigma = 0.20 / 0.54 / 1.50 / 3.00 — six to eight orders below `POLR_SCORE_TOL`, at every sigma
    the bootstrap reaches. The approximation is named rather than hidden, and §20.12 pins the
    measurement so a change to the mode search that made it worse fails a test.

    Differentiating `log L_c = c_c + LSE_q[term_cq]` gives `SUM_q pi_cq * d term_cq / d theta` with
    `pi_cq` the softmax over the nodes, so every derivative below is a posterior-weighted average of
    a per-node one. `alpha` and `beta` come through `s_cq` and are `_ord_score_hess`' own; the
    `log sigma` derivative is where the two rules differ, because non-adaptive nodes MOVE with sigma
    and adaptive ones are frozen.
    """
    alpha, beta = par[:K], par[K:K + m]
    sigma = float(np.exp(par[-1]))
    Xb = Xn @ beta
    n_nodes = len(t)

    if adaptive:
        b_hat, tau = _ri_modes(Xb, alpha, y_idx, K, gidx, n_groups, sigma)
        bq = b_hat[:, None] + _SQRT2 * tau[:, None] * t[None, :]          # (G, Q)
        const = np.log(_SQRT2 * tau)
        prior = -0.5 * np.log(2.0 * np.pi * sigma ** 2) - bq ** 2 / (2.0 * sigma ** 2)
        extra = logw[None, :] + t[None, :] ** 2 + prior
    else:
        b_hat, tau = np.zeros(n_groups), np.full(n_groups, sigma)
        bq = np.broadcast_to(sigma * _SQRT2 * t[None, :], (n_groups, n_nodes)).copy()
        const = np.full(n_groups, -0.5 * np.log(np.pi))
        extra = np.broadcast_to(logw[None, :], (n_groups, n_nodes)).copy()

    # One (n, Q) evaluation for the whole sample: every observation at every one of ITS group's
    # nodes. `y_idx[:, None]` is what makes `_ord_b_derivatives` broadcast over the node axis.
    p, A, B, _ = _ord_b_derivatives(Xb[:, None] + bq[gidx], alpha, y_idx[:, None], K)
    s = np.zeros((n_groups, n_nodes))
    np.add.at(s, gidx, np.log(p))

    term = extra + s
    peak = term.max(axis=1, keepdims=True)
    lse = peak[:, 0] + np.log(np.exp(term - peak).sum(axis=1))
    weights = np.exp(term - lse[:, None])                                 # (G, Q), rows sum to 1

    # The posterior weight of each observation's own node, so every per-observation derivative below
    # is already averaged over the nodes before it reaches the cutpoint accumulation.
    w_iq = weights[gidx]
    ra = (w_iq * A).sum(axis=1)
    rb = (w_iq * B).sum(axis=1)

    grad = np.zeros(K + m + 1)
    has_u, has_l = y_idx <= K - 1, y_idx >= 1
    np.add.at(grad, y_idx[has_u], ra[has_u])                              # the UPPER cutpoint
    np.add.at(grad, (y_idx - 1)[has_l], rb[has_l])                        # the lower one
    grad[K:K + m] = Xn.T @ (ra + rb)
    if adaptive:
        # d/d log_sigma of  -log sigma - b_q^2/(2 sigma^2), the nodes being frozen.
        grad[-1] = float(np.sum(weights * (-1.0 + bq ** 2 / sigma ** 2)))
    else:
        # The nodes MOVE: d b_q / d log_sigma = b_q, so the chain rule runs through `s`.
        d1 = np.zeros((n_groups, n_nodes))
        np.add.at(d1, gidx, A + B)
        grad[-1] = float(np.sum(weights * d1 * bq))

    return float(np.sum(const + lse)), grad, b_hat, tau


def _polr_ri(X: pd.DataFrame, y: np.ndarray, groups: pd.Series | np.ndarray,
             nodes: int, adaptive: bool) -> RIFit:
    """`polr_ri`'s body, with the two things the SPECIFICATION fixes exposed as parameters.

    `polr_ri` below passes `C.POLR_RI_NODES` and `adaptive=True` and is the only caller any shipped
    module has. This private form exists for ONE reason: Stage 12 §20.12 requires the measurement
    §12.2 rests on asserted as a test — that sigma_hat at 9, 11, 15 and 31 ADAPTIVE nodes agrees to
    1e-7 while sigma_hat at 5, 7 and 9 NON-ADAPTIVE nodes does not — and that assertion cannot be
    made against a function with no way to ask for the rule the specification rejected.

    It is private, and `adaptive=False` appears nowhere outside `test_model.py`, so the rejected rule
    is reachable by a test and not by a caller.
    """
    Xn = np.asarray(X.to_numpy(dtype=float))
    y = np.asarray(y, dtype=float)
    w = np.ones(len(y))
    _assert_polr_fittable(Xn, y, w, tuple(X.columns))                      # O1-O6, unweighted

    labels = np.asarray([str(g) for g in np.asarray(groups)])
    if len(labels) != len(y):
        raise C.SchemaError(
            f"T4  the grouping vector holds {len(labels)} label(s) for {len(y)} response value(s). "
            "The two are aligned by POSITION, so a mismatch groups the wrong patients together and "
            "returns a between-group variance for a grouping nobody specified.")

    categories = _weighted_categories(y, w)                                # one definition
    keep = np.isin(y, np.asarray(categories, dtype=float))
    Xn, y, labels = Xn[keep], y[keep], labels[keep]
    K, m = len(categories) - 1, Xn.shape[1]
    y_idx = np.searchsorted(np.asarray(categories, dtype=float), y)

    # THE GROUP SET IS THE ONE THAT SURVIVES THE COLLAPSE, and it is sorted for `resample`'s reason
    # (Stage 10 §4.1): the order is the label's own and never the frame's, so `b` is keyed
    # deterministically regardless of how the rows arrived.
    unique = tuple(sorted(set(labels.tolist())))
    if len(unique) < 2:
        raise FitError(
            f"T4  the grouping variable has {len(unique)} group(s) with records: {unique}. A "
            "between-group variance needs at least two: with one, `b_c` is exactly collinear with "
            "the cutpoints and sigma is not identified — the likelihood is flat in it, so the fit "
            "returns whatever the start value was, finite and plausible.")
    gidx = np.searchsorted(np.asarray(unique), labels)
    n_groups = len(unique)

    # THE START VALUES ARE THE POOLED FIT'S, and Stage 12 §12.8 is why that is load-bearing rather
    # than convenient: `(alpha, beta)` from `polr` is a nested model's EXACT maximiser in every
    # coordinate but one, which is why 25 BFGS iterations suffice and why the arm costs 17.5 minutes
    # rather than hours. It is also why a pooled `FitError` costs the hierarchical arm too
    # (Stage 12 §13.2) — there is no hierarchical arm without a pooled fit, and this call is where
    # that dependency is structural rather than merely stated.
    start = polr(pd.DataFrame(Xn, columns=list(X.columns)), y)
    log_floor = float(np.log(C.POLR_RI_SIGMA_FLOOR))
    par = np.concatenate([start.alpha, start.beta, [np.log(0.5)]])

    t, quad_w = np.polynomial.hermite.hermgauss(nodes)
    logw = np.log(quad_w)

    def objective(theta: np.ndarray) -> tuple[float, np.ndarray, np.ndarray, np.ndarray]:
        return _ri_objective(theta, Xn, y_idx, K, m, gidx, n_groups, t, logw, adaptive)

    def project(theta: np.ndarray) -> np.ndarray:
        """`log sigma` is held at or above the floor rather than clipped inside the objective.

        A clip inside would make the objective FLAT below the floor, where the true derivative is
        zero and a gradient reported from the clipped point traps the iterate there permanently. A
        projection keeps every evaluation at a feasible point, so the gradient is always the real one
        and a fit sitting at the floor is sitting there because the data put it there
        [Stage 12 §12.5].
        """
        out = theta.copy()
        out[-1] = max(out[-1], log_floor)
        return out

    par = project(par)
    ll_old, grad, b_hat, tau = objective(par)
    inverse_hessian = np.eye(K + m + 1)
    first_step_norm, rescales, halvings = 0.0, 0, 0
    moved = np.inf

    for iteration in range(1, C.POLR_RI_MAX_ITER + 1):
        step = inverse_hessian @ grad                    # ASCENT: B approximates (-H)^-1
        norm = float(np.linalg.norm(step))
        if iteration == 1:
            first_step_norm = norm
        size = norm / max(1.0, float(np.linalg.norm(par)))
        if size > C.POLR_MAX_STEP:
            step = step * (C.POLR_MAX_STEP / size)
            rescales += 1

        # Step-halving PER ITERATION, `polr`'s discipline transplanted. Stage 12 §12.4 declares the
        # one departure: the TOTAL across the fit is much larger — measured 85 over 25 iterations —
        # because BFGS proposes long steps early and the line search shortens them. That is the
        # algorithm working, and the total is recorded so a change that makes it grow is visible.
        for _ in range(C.POLR_MAX_HALVINGS):
            trial = project(par + step)
            ll_new, grad_new, b_new, tau_new = objective(trial)
            if np.isfinite(ll_new) and ll_new >= ll_old:
                break
            step = step / 2.0
            halvings += 1
        else:
            raise FitError(
                f"polr_ri: step-halving exhausted {C.POLR_MAX_HALVINGS} halvings at iteration "
                f"{iteration} without increasing the marginal log-likelihood. [§14a] prescribes one "
                "estimator; there is no second one to try [invariant 5].")

        # The BFGS update is on the MINIMISATION problem, so `y_curv` is the change in the gradient
        # of the NEGATIVE log-likelihood. The curvature condition is checked rather than assumed: a
        # non-positive `y_curv . s` makes the update indefinite, and the safeguard is to keep the
        # previous approximation rather than to take a step along an ascent direction that is not one.
        s_vec = trial - par
        y_curv = -(grad_new - grad)
        curvature = float(y_curv @ s_vec)
        if curvature > 0.0:
            rho = 1.0 / curvature
            left = np.eye(K + m + 1) - rho * np.outer(s_vec, y_curv)
            inverse_hessian = left @ inverse_hessian @ left.T + rho * np.outer(s_vec, s_vec)

        par, grad, b_hat, tau = trial, grad_new, b_new, tau_new
        moved = abs(ll_new - ll_old)
        ll_old = ll_new

        # TWO CONVERGENCE ROUTES, either one sufficient, in `polr`'s order and for its reason: [§10]
        # refits in every replicate, so reordering them changes the sampling distribution. Measured
        # (Stage 12 §12.5): 200 of 200 replicates converge, every one on "likelihood" — which is what
        # a projected iterate at the floor looks like, the step in `log sigma` being zero there.
        if moved < C.POLR_TOL:
            return _ri_fit(par, K, m, categories, tuple(X.columns), unique, b_hat, tau, nodes,
                           iteration, "likelihood", first_step_norm, rescales, halvings, log_floor)
        if float(np.max(np.abs(grad))) < C.POLR_SCORE_TOL:
            return _ri_fit(par, K, m, categories, tuple(X.columns), unique, b_hat, tau, nodes,
                           iteration, "score", first_step_norm, rescales, halvings, log_floor)

    raise FitError(
        f"polr_ri: no convergence in {C.POLR_RI_MAX_ITER} iterations. The marginal log-likelihood "
        f"moved {moved:.3g} against a tolerance of {C.POLR_TOL:g}, and the largest gradient "
        f"component was {float(np.max(np.abs(grad))):.3g} against {C.POLR_SCORE_TOL:g}. "
        f"{rescales} step(s) were shortened by the trust region and {halvings} halving(s) taken. "
        "[§14a] prescribes one estimator, and [§10] drops and counts the replicate rather than "
        "substituting another.")


def _ri_fit(par: np.ndarray, K: int, m: int, categories: tuple[int, ...],
            columns: tuple[str, ...], groups: tuple[str, ...], b_hat: np.ndarray,
            tau: np.ndarray, nodes: int, iterations: int, converged_on: str,
            first_step_norm: float, rescales: int, halvings: int, log_floor: float) -> RIFit:
    """`RIFit` from the optimiser's state. One constructor, two call sites, no repeated field list."""
    # AT THE FLOOR, `sigma` IS THE FLOOR EXACTLY, and the snap is not cosmetic:
    # `exp(log(1e-4))` is 0.00010000000000000009, so a fit that reported it would make
    # `sigma == C.POLR_RI_SIGMA_FLOOR` False on a fit whose `at_floor` is True -- two fields
    # disagreeing about one fact, which is the shape of defect Stage 12 §3.1 declines a discriminator
    # field to avoid. §20.12 asserts the equality, so this is what makes it assertable.
    at_floor = bool(par[-1] <= log_floor + 1e-12)
    sigma = C.POLR_RI_SIGMA_FLOOR if at_floor else float(np.exp(par[-1]))
    return RIFit(
        beta=par[K:K + m].copy(), alpha=par[:K].copy(), sigma=sigma, groups=groups,
        b=b_hat.copy(), b_sd=tau.copy(), categories=categories, columns=columns, nodes=nodes,
        iterations=iterations, converged_on=converged_on, first_step_norm=first_step_norm,
        rescales=rescales, halvings=halvings, at_floor=at_floor)


def polr_ri(X: pd.DataFrame, y: np.ndarray, groups: pd.Series | np.ndarray) -> RIFit:
    """Random-intercept proportional-odds regression [§14a, Stage 12 §12]. Raises FitError.

        logit P(Y <= k | x, b_c)  =  alpha_k + x'beta + b_c ,      b_c ~ N(0, sigma^2)

    ONE common `beta` across groups — no `A x group` term and no random slope, both forbidden by
    [§14a] and [§15] at four centres, neither being identifiable in any useful sense there.

    `X` is the same design `polr` fits and carries the exposure column; **the grouping variable is
    NOT in it and could not be** — a fixed group effect and a random one are the same parameter
    twice. `groups` is aligned by POSITION with `y`, as `smd`'s arguments are, because the caller has
    already masked both to the same records.

    THE PARAMETRISATION IS `polr`'s AND THEREFORE [§14a]'s, `alpha_k + x'beta`, and `ordinal::clmm`
    writes `zeta_k - x'beta` instead. Measured against it the asymmetry is THREE-WAY: the
    coefficients come back negated, the thresholds do NOT, and **`sigma` does not either** because it
    is a scale. `tests/reference/polr_ri_clmm.R` asserts all three rather than that "the fits agree",
    which would pass on the thresholds alone (Stage 12 §12.7).

    Reads `C.POLR_RI_NODES` and takes no node argument, for the reason `polr` takes no tolerance:
    the node count is part of the definition of the objective and [§10] refits this model in every
    one of `N_BOOT` replicates, so it is prespecified rather than a runtime knob (Stage 12 §12.3).

    **A FIT AT `POLR_RI_SIGMA_FLOOR` IS NOT A FAILURE AND DOES NOT RAISE.** [§14a] names
    sigma^2_C = 0 as a legitimate answer, `RIFit.at_floor` records it, and Stage 12 §12.5 measured it
    in 44 of 200 replicates. Dropping those would select the bootstrap on the value of the parameter
    the arm exists to examine.

    **There is no standard error, for `PolrFit`'s reason** [Stage 8 §5.6], and one more specific to
    this fit: at four clusters a variance parameter's asymptotic standard error is the least
    trustworthy number in the output. [§10]'s percentile bootstrap is the prespecified interval.
    """
    return _polr_ri(X, y, groups, C.POLR_RI_NODES, adaptive=True)
