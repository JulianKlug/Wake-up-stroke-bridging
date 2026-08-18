"""Stage 6a — the pipeline's numerics [§6, §7].

The design matrix and the Firth fit, and nothing that knows what the exposure is. This module is
**outcome-agnostic by construction**: it reads ``config`` and nothing else, takes no ``Audit``, and
names neither the treatment nor any covariate list. Its caller logs; it computes.

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
