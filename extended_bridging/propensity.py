"""Stage 6b — the [§7] propensity score and overlap weights.

The [§7] specification, and the only module in this pipeline that knows what the exposure is. The
numerics it stands on — the design matrix and the Firth fit — are ``model.py``'s, because Stage 9's
outcome regression and Stage 12's standardisation model need the same two functions with a different
covariate list and a different response [Stage 6 §0.1].

**What ``fit`` returns, and the one rule every later stage owes.** ``e`` and ``w`` are Series over the
**whole** cohort, carrying ``np.nan`` wherever ``in_model`` is False — not over the fitted subset. The
sentinel is ``np.nan`` and not ``pd.NA``: a ``float64`` Series has exactly one missing value and
``pd.NA`` is not it, and the masked ``Float64`` dtype is declined because every downstream
``.to_numpy(dtype=float)`` would convert it back to ``nan`` silently, on the one column most likely to
be summed. So **the deliberateness of the absence lives in ``in_model``**, which is boolean, total and
never missing — not in the value.

That distinction is the whole of this stage's second decision. **A weight of zero is a patient who was
weighed and found irrelevant; an absent weight is a patient who was never weighed.** The two are the
same in every sum and different in every denominator, and [§11] requires the denominator. ``nan`` is
the safe direction — a stage that forgets the mask gets ``nan`` rather than a plausible wrong number —
but that holds **only while nothing fills it**. A ``fillna(0)`` anywhere downstream converts "never
weighed" into "weighed and found irrelevant", silently, in the direction that changes an estimate
rather than breaking it, and it is the natural thing to reach for when a weighted mean comes back
``nan``. So: **range over ``in_model``, never over ``notna()``, and never fill.**

This stage **adds no column to the cohort frame and edits no value.** A ``propensity`` column would be
a second place the score lives, and [§10] refits in every one of ``N_BOOT`` replicates: a frame
carrying the point fit's column, resampled into a replicate, is a replicate silently weighted by the
wrong score [Stage 6 §0.2].

Section references in brackets are to ``statistical_analysis_plan.md``. The specification for this
module is ``specs/stage6_propensity_and_weights.md``; nothing here is invented outside it.

::

    cohort.build(df, audit)  →  93 rows x 33 columns                     [Stage 5 §11]
                       │
                       ▼
    ┌──────────────────────────────────────────────────────────────────────────┐
    │  STAGE 6b — propensity.py        reads config.py, data.Audit, model      │
    │                                                                          │
    │   fit(df, audit)      -> _fit(df, PROPENSITY_PRIMARY, audit)   [§7]      │
    │   fit_full(df, audit) -> _fit(df, PROPENSITY_FULL, audit)      [§13]     │
    │                                                                          │
    │   _fit(df, spec, audit) -> Propensity                                    │
    │     ├─ _assert_fit_inputs(df)                F1, F2 → SchemaError        │
    │     ├─ complete_cases(df, spec.covariates)   92 of 93        (§4.4)      │
    │     ├─ _record_exclusion(...)                entry 1, spec.step          │
    │     ├─ design(df[in_model], spec.covariates) 11 or 16 columns (§4)       │
    │     ├─ firth(X, a)                           7 iterations    (§5)        │
    │     ├─ _assert_probabilities(...)             F5 → FitError   (§5.5)     │
    │     ├─ weights, per-arm ESS                                  (§6)        │
    │     └─ four `model` audit entries, spec.step(base)           (§7)        │
    └──────────────────────────────────────────────────────────────────────────┘
                       │
                       ▼
     Propensity(e, w, in_model, ess, fit, dropped, spec)  →  Stages 7-11

**There are TWO named entry points and no covariate keyword** [Stage 11 §4]. Stage 6 §6.4 declined a
``covariates=`` parameter because the [§13] propensity-specification rows are *different target
populations* and a parameter makes them look like options, and it left a promissory note — that when
one of them is un-deferred it arrives *"as a named specification with a [§13] amendment, not as a
keyword"*. DECISION 4 un-deferred exactly one, so the note is paid here: ``C.Specification`` is a
frozen record, the registry holds exactly two instances, and each has its own public entry point
taking ``(df, audit)`` and nothing else. **The record and not a bare covariate tuple**, because the
seam has four jobs and a tuple does one: the audit step names must not collide, and ``Audit.entry``
is first-match — so two fits recording ``propensity_fit`` would make every programmatic read return
the primary's while the rendered log looked complete [Stage 11 §4.2, §4.4].

``_fit`` appends four entries of the ``model`` kind to the ``Audit`` that ``load()`` created, each
named ``spec.step(base)``. The primary's suffix is the empty string, so every step name Stages 6 and
7 already record is byte-identical. It writes no file, as no stage before it does — and the audit log is the **only** place in this repository the
fitted coefficients ever appear [Stage 6 §7.3].

This module is **not** exempt from the Stage 1 §7 raw-name scan and must never become exempt.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd

import config as C
import model
# `_fmt` is private to data.py and is imported anyway, for the reason derive.py gives: it is the
# pipeline's *one* float formatter, and a second one is a second way for two runs to disagree. Every
# coefficient, ESS and weighted mean below reaches the log through it, at six significant figures,
# which is what makes §12.11's byte-identity criterion hold.
from data import Audit, _fmt


# --- what the stage returns ------------------------------------------------------------------------

@dataclass(frozen=True)
class Propensity:
    """The [§7] propensity score and overlap weights, over the frame `fit` was given.

    A return value and never a column on the frame [Stage 6 §0.2], so a caller that wants the point
    score and a replicate's score at once has to hold two objects and name them.

    **``spec`` is appended and carries NO DEFAULT** [Stage 11 §4.3]. A default of
    ``C.PROPENSITY_PRIMARY`` would keep every direct construction site green, and that is exactly
    the hazard: an arm-derived ``Propensity`` that forgot the field would role-label in
    ``balance._role`` as the primary's, and every number downstream would read as the [§7]
    specification's. Without a default those sites raise ``TypeError`` — loud, and each a one-line
    fix. ``dataclasses.replace`` carries the field correctly and stays the safe way to bend a
    ``Propensity``.
    """

    e: pd.Series                # float64, the cohort's index, np.nan where not in_model (§8)
    w: pd.Series                # float64, same index and same missingness
    in_model: pd.Series         # boolean, TOTAL — never missing; this is where the intent lives
    ess: dict[int, float]       # keyed by arm code, as TREATMENT_LABELS is
    fit: model.Fit
    dropped: tuple[str, ...]    # design columns dropped as constant, in design order
    spec: C.Specification       # which prespecified specification this is — Stage 11 §4.2


# --- the effective sample size [§7] ----------------------------------------------------------------

def ess(w: pd.Series | np.ndarray) -> float:
    """Kish effective sample size, (Σw)² / Σw². Raises on an empty arm and on a zero-weight arm.

    Reported per arm [§7]. `pilots/analysis.py:140` returns 0.0 when no weight is positive; that is
    declined. An arm whose weights all vanish is an arm in which the fit assigned every patient a
    propensity of exactly 0 or 1 — a structural non-positivity the [§3] restrictions were supposed
    to have removed — and 0.0 is a number every downstream ratio will accept [Stage 6 §3.1].

    Public, and it lives here rather than in `model.py` for two reasons. Stage 7's within-centre
    overlap table needs it per centre per arm [§9], and a second implementation there would be a
    second definition of the effective sample size; and a Kish sum has nothing to do with fitting
    anything, so it would breach `model.py`'s no-exposure-no-weights rule [Stage 6 §0.1].

    `FitError` is qualified as `model.FitError`: it is declared in model.py (§5.5) and this module
    imports the module, not the name. Stage 7 calls this per centre per arm, where the empty branch is
    the one that fires — a centre with no patient in one arm is exactly what `treating_centres`
    excludes, and this raise is what stops Stage 7 plotting it as a zero.

    `dropna()` first, so passing the full 93-row `w` and the 92-row masked `w` give the same number: a
    Kish sum over a nan is nan, and the tempting "fix" for that is a `fillna(0)` this module's
    docstring forbids.
    """
    v = np.asarray(pd.Series(w).dropna(), dtype=float)
    if not v.size:
        raise model.FitError(
            "ESS: this arm has no weighted patient at all, so the effective sample size is 0/0 over "
            "an empty sum. This is NOT the all-zero-weight case below and must not report as it: an "
            "empty arm means the in_model mask lost one, since F2 established both arms were present "
            "before the fit [Stage 6 §8]. Measured: np.sum([] ** 2) is 0.0, so the two reach the "
            "same branch and only a separate check tells them apart [§3.1].")
    total = float(np.sum(v ** 2))
    if not total:
        raise model.FitError(
            "ESS: every weight in this arm is zero, so the effective sample size is 0/0. The fit "
            "put every patient in the arm at a propensity of exactly 0 or 1 — structural "
            "non-positivity, which [§3]'s restrictions remove by design [Stage 5 §5.2 P2].")
    return float(np.sum(v) ** 2 / total)


# --- the preconditions ------------------------------------------------------------------------------

def _assert_fit_inputs(df: pd.DataFrame) -> None:
    """F1 and F2, collected: one SchemaError reports both. Before the mask and before the design."""
    bad: list[str] = []

    off_arm = df.loc[~df[C.TREATMENT].isin((0, 1)), "case_id"]
    if len(off_arm):
        bad.append(
            f"F1  {C.TREATMENT}: {len(off_arm)} record(s) whose exposure is missing or outside "
            f"{{0, 1}}: {', '.join(sorted(off_arm))}. The response of the propensity model is the "
            "exposure; a missing one is a row with no y, and Int64's pd.NA becomes nan silently "
            "[§3.1]. Stage 4's E5 and Stage 5's C3 check the same property upstream.")

    present = [code for code in C.TREATMENT_LABELS if int((df[C.TREATMENT] == code).sum())]
    if len(present) < 2:
        bad.append(
            "F2  " + ", ".join(f"{C.TREATMENT_LABELS[c]}: "
                               f"{int((df[C.TREATMENT] == c).sum())}" for c in C.TREATMENT_LABELS)
            + ". A propensity model needs both arms: with one, the fit is a model of a constant and "
            "every overlap weight is 0 or 1. Stage 5's P2 guarantees both arms at every retained "
            "centre, so this can only be a caller subsetting after build().")

    if bad:
        raise C.SchemaError(
            "\n  " + "\n  ".join(bad)
            + f"\n  {len(bad)} propensity assertion(s) failed against {len(df)} records.")


def _assert_probabilities(fitted: model.Fit, case_ids: pd.Series) -> None:
    """F5 — the fitted probabilities are finite and strictly inside (0, 1).

    Not a clip. `pilots/analysis.py:82` clips into [1e-8, 1-1e-8]; a probability of exactly 1.0
    means the linear predictor hit FIRTH_ETA_CLIP [§3.1], which under a penalised likelihood means the
    fit is degenerate, and a clip turns that into a weight of 0.0 that every sum accepts [§5.5].

    It sits between the fit and the weights, because a degenerate `e` becomes a weight of exactly 0.0
    one line later and stops being visible.
    """
    bad = ~np.isfinite(fitted.p) | (fitted.p <= 0.0) | (fitted.p >= 1.0)
    if bad.any():
        raise model.FitError(
            f"F5  {int(bad.sum())} fitted probability/ies are not strictly in (0, 1): "
            f"{', '.join(sorted(case_ids.to_numpy()[bad]))}. Firth's penalty exists to keep the fit "
            "off the boundary, so a boundary value is a degenerate fit and not a confident one. "
            "[§10] drops and counts the replicate; it does not clip it into range.")


# --- the four audit entries [§7] ---------------------------------------------------------------------
#
# **Every one of the eight helpers below takes the `C.Specification`** [Stage 11 §4.3], because each
# either read `C.PS_COVARIATES` or asserted "[§6]", and three of their `detail` strings are FALSE
# under the [§13] specification — one of them with no referent at all: `_weights_detail` said
# "conditional on this specification", and with two fits in one log "this" points at nothing.
#
# Five rules about the tables, three of them inherited from Stages 4 and 5:
#
#   * Every declared thing is rendered whether or not the data fills it. `design_matrix` has a row for
#     an absent centre reading `dropped: constant`, because a column that vanishes from a table is a
#     column nobody knows was considered; `overlap_weights` has a row per TREATMENT_LABELS entry.
#     `pd.crosstab` is used nowhere in this repository, for that reason.
#   * Every cell is computed from the object it describes, never by subtracting another row: a table
#     that reconciles with itself by construction has stopped being evidence [Stage 5 §7.3].
#   * `_fmt` is the only float formatter and `missing` the only rendering of an absent value.
#   * The `all` row of `overlap_weights` carries an ESS computed over both arms POOLED, and says so —
#     not the sum of the two arms', which is not a Kish sum of anything.
#   * Every step name is `spec.step(base)` and never a literal, because `Audit.entry` is first-match
#     and the arm records the same four bases over the same cohort [Stage 11 §4.4].


def _completeness_detail(df: pd.DataFrame, in_model: pd.Series, spec: C.Specification) -> str:
    n_excluded = int((~in_model).sum())
    return (
        f"[§11] complete-case on the {spec.label} propensity specification's "
        f"{len(spec.covariates)} covariate(s) {spec.sap}: {int(in_model.sum())} of {len(df)} cohort "
        f"record(s) carry every one and are fitted; {n_excluded} record(s) do not, and every one of "
        "them is named above. Such a record stays in the cohort and loses its weight — `e` and `w` "
        "are absent for it and "
        "`in_model` is False — because a covariate that was never recorded is not a [§3] restriction. "
        "So the ATO target population is this complete-case set rather than the [§3] cohort, which is "
        "stated here rather than left implicit [Stage 6 §4.4]. It is a property of THIS "
        "specification's covariate list and not of the cohort: a specification adding a covariate "
        "that carries its own missingness excludes more records, so two specifications' denominators "
        "have to be compared rather than assumed equal [Stage 11 §5.2]. The per-centre breakdown of "
        "every "
        "absence is in `absence_by_cohort_column`, under Missingness and denominators; this entry "
        "carries only what it uniquely owns, which is what the missingness cost THIS fit.")


def _completeness_table(df: pd.DataFrame, in_model: pd.Series,
                        spec: C.Specification) -> tuple[tuple[str, ...], ...]:
    """One row per covariate of `spec`, in its declared order — every one, complete or not.

    `excluded_by` is the column that earns this entry. `absence_by_cohort_column` answers "what is
    missing"; only this answers "what did the missingness cost this fit", which is [§11]'s
    per-estimate denominator [Stage 6 §7.4].
    """
    header = ("covariate", "n_absent", "excluded_by")
    rows = tuple(
        (c, str(int(df[c].isna().sum())),
         "yes" if int((df[c].isna() & ~in_model).sum()) else "no")
        for c in spec.covariates)
    return (header, *rows)


def _design_detail(X: pd.DataFrame, dropped: tuple[str, ...], spec: C.Specification) -> str:
    return (
        f"{spec.label} design matrix {spec.sap}: {X.shape[1]} column(s) plus an intercept over "
        f"{X.shape[0]} fitted "
        f"record(s), reference-coded on REFERENCE_LEVELS. {len(dropped)} column(s) dropped as "
        f"constant{': ' + ', '.join(dropped) if dropped else ''}. Factor levels are DECLARED "
        "(FACTOR_LEVELS), not observed, so the width is a property of config.py and not of the "
        "sample, and a level absent from this frame appears below as `dropped: constant` rather than "
        "vanishing. The reference level of each factor contributes no parameter and is marked as such "
        f"— every other coefficient is relative to it. n is the PARAMETER count including the "
        f"intercept, {X.shape[1] + 1}, which is what the [§13] degrees-of-freedom budget is computed "
        "against; it is not the row count of this table.")


def _design_table(X: pd.DataFrame, dropped: tuple[str, ...],
                  spec: C.Specification) -> tuple[tuple[str, ...], ...]:
    """One row per level of every declared factor, plus one per linear covariate. §7.5.

    Ranges over the DECLARATION, never over the matrix: a level absent from both X and `dropped` is a
    reference level, and a table built from `X.columns + dropped` would omit exactly those — which is
    a worse omission here than anywhere else, because the audit log is the only place the fitted
    coefficients appear. `center_CHUV = 0.41` means *relative to HUG*; without HUG on the page the
    number is uninterpretable.

    It takes the SPECIFICATION and reads `spec.covariates` off it, and that is not incidental: from
    `X.columns` and `dropped` alone the
    reference rows can only be *inferred*, by looking for `{factor}_{level}` name prefixes among the
    surviving columns and guessing which declared level is missing — which happens to work and is the
    wrong shape, since a linear covariate sharing a prefix would be misread as a dummy and
    REFERENCE_LEVELS has to be consulted regardless.

    Row order is `spec.covariates`' order, which is NOT the design matrix's column order (get_dummies
    appends the dummies). That is deliberate: the matrix's order exists to map a coefficient to a name,
    and this table's order exists to be read by a human against the declared covariate list.
    """
    header = ("column", "term", "factor", "level", "status")
    rows: list[tuple[str, ...]] = []
    for c in spec.covariates:
        if c not in C.CATEGORICAL:
            rows.append((c, "linear", "—", "—",
                         "dropped: constant" if c in dropped else "kept"))
            continue
        for level in C.FACTOR_LEVELS[c]:
            name = f"{c}_{level}"
            status = ("reference (baseline)" if level == C.REFERENCE_LEVELS[c]
                      else "dropped: constant" if name in dropped
                      else "kept")
            rows.append((name, "dummy", c, level, status))
    return (header, *rows)


def _fit_detail(fitted: model.Fit, e: pd.Series, spec: C.Specification) -> str:
    return (
        f"[§7] Firth-penalised logistic regression under the {spec.label} specification "
        f"{spec.sap}, converged in {fitted.iterations} iteration(s) on "
        f"the {fitted.converged_on} criterion. Fitted probabilities range "
        f"{_fmt(e.min())} to {_fmt(e.max())}, all finite and strictly interior (F5). "
        f"{fitted.rescales} step(s) were shortened by the trust region and {fitted.halvings} "
        "halving(s) taken [§5.2a]. The penalty is what keeps the fit off a boundary an unpenalised "
        "likelihood runs to under near-separation, which is why [§7] names it; there is no second "
        "estimator on any path, and a fit that does not converge raises rather than falling back. "
        "The coefficients below are the only place in this repository they appear, and they are "
        "relative to each factor's reference level as the design table marks it.")


def _fit_table(fitted: model.Fit) -> tuple[tuple[str, ...], ...]:
    """One row per parameter, in `Fit.columns`' order — which IS `beta`'s order after the intercept.

    A wrong order here is correctly-computed numbers under wrong labels, which is the one error a
    reader of the log cannot detect.
    """
    header = ("parameter", "coefficient")
    names = ("intercept",) + fitted.columns
    return (header, *((name, _fmt(value)) for name, value in zip(names, fitted.beta)))


def _weights_detail(spec: C.Specification) -> str:
    return (
        "[§7] overlap weights: w = 1 − e for the treated arm and w = e for the control arm, the ATO "
        "tilt h(X) = e(X){1 − e(X)} split between the arms. Bounded in [0, 1] by construction and "
        "largest at e = 0.5 — patients near equipoise. No trimming, no truncation, no stabilisation: "
        "there is nothing to trim, because the pathology trimming exists for (an inverse-probability "
        "weight exploding as e → 0) cannot arise. "
        "The ATO target population is h(X) = e(X){1 − e(X)}, so CHANGING THE PROPENSITY "
        "SPECIFICATION CHANGES THE POPULATION, not merely the precision of the estimate — and both "
        "the effective sample sizes and the weighted-population description below are conditional "
        f"on the {spec.label} specification {spec.sap} and on no other. **The specification is NAMED "
        "here rather than called `this`**, because more than one propensity model is fitted over this "
        "cohort and `this` points at nothing in a log carrying two of them [§7, Stage 6 §6.4, "
        "Stage 11 §4.3]. The second table describes who that weighted "
        "population comprises: overlap weighting up-weights patients near equipoise and down-weights "
        "those whose treatment was nearly determined, so it is not the unweighted cohort. It is a "
        "description and not a balance assessment — [§9]'s standardised mean differences are Stage "
        "7's, and this stage does not judge the balance it reports.")


def _weights_table(df: pd.DataFrame, in_model: pd.Series, w: pd.Series,
                   per_arm: dict[int, float]) -> tuple[tuple[str, ...], ...]:
    """One row per declared arm plus `all`, whether or not the data fills it.

    The `all` row's ESS is computed over both arms POOLED and is labelled, because the sum of two Kish
    sums is not a Kish sum. `weight share` is the arm's share of the total weight, which is the number
    that says how much of the estimate each arm is carrying.
    """
    header = ("arm", "n", "sum w", "ESS", "max w", "weight share")
    total = float(w[in_model].sum())
    rows: list[tuple[str, ...]] = []
    for code, label in C.TREATMENT_LABELS.items():
        arm = in_model & (df[C.TREATMENT] == code)
        values = w[arm]
        rows.append((label, str(int(arm.sum())), _fmt(values.sum()), _fmt(per_arm[code]),
                     _fmt(values.max()), _fmt(values.sum() / total)))
    pooled = w[in_model]
    rows.append(("all (pooled)", str(int(in_model.sum())), _fmt(pooled.sum()),
                 _fmt(ess(pooled)), _fmt(pooled.max()), _fmt(1.0)))
    return (header, *rows)


def _weighted_population_table(df: pd.DataFrame, in_model: pd.Series, w: pd.Series,
                               spec: C.Specification) -> tuple[tuple[str, ...], ...]:
    """One row per covariate of `spec` — per LEVEL for a factor — unweighted and ATO-weighted. §7.6.

    [§7] asks for two things and the ESS discharges only one: "report the effective sample size per
    arm and DESCRIBE THE WEIGHTED POPULATION. Because the ATO population is defined statistically,
    state explicitly who it comprises." A sentence saying the population depends on the specification
    is a caveat *about* the population, not a description *of* it.

    Ranges over the declaration for the same reason `_design_table` does: a level nobody in the
    weighted population has is a fact about the population, and a row that vanishes is a fact nobody
    sees. A continuous covariate contributes its mean, a factor level its proportion.

    `np.average(..., weights=...)` is the weighted mean WITHIN an arm — what "who does this arm
    comprise after weighting" asks. It is not a between-arm contrast, and nothing here is a
    standardised mean difference or a negative control: those are Stage 7's [§9].
    """
    header = ("covariate", *(f"{lab} ({kind})"
                             for kind in ("unweighted", "ATO-weighted")
                             for lab in C.TREATMENT_LABELS.values()))
    sub, ww = df.loc[in_model], w.loc[in_model].to_numpy(dtype=float)
    arms = {code: (sub[C.TREATMENT] == code).to_numpy() for code in C.TREATMENT_LABELS}

    def cells(v: np.ndarray) -> tuple[str, ...]:
        return (*(_fmt(v[a].mean()) for a in arms.values()),
                *(_fmt(np.average(v[a], weights=ww[a])) for a in arms.values()))

    rows: list[tuple[str, ...]] = []
    for c in spec.covariates:
        if c in C.CATEGORICAL:
            for level in C.FACTOR_LEVELS[c]:
                rows.append((f"{c} = {level}",
                             *cells((sub[c] == level).to_numpy(dtype=float))))
        else:
            rows.append((c, *cells(sub[c].to_numpy(dtype=float))))
    return (header, *rows)


def _padded(table: tuple[tuple[str, ...], ...], width: int) -> tuple[tuple[str, ...], ...]:
    """Widen every row to `width` with the not-applicable dash `_design_table` already uses.

    **This exists because §7.2 asks one audit entry to carry TWO tables and `AuditEntry` holds one.**
    `overlap_weights` is a single entry by declaration — §12.11 asserts exactly four `model` entries in
    `fit`'s order, so a fifth is not available — and its two tables are 6 and 5 columns wide. §8's
    literal line, `_weights_table(...) + _weighted_population_table(...)[1:]`, therefore does two
    things it did not intend: `data._md_table` computes its column widths from `rows[0]` and raises
    IndexError on the narrower rows, and the `[1:]` discards the population table's header, which is
    what names its four arm-by-weighting columns. Without that header the second block would render
    under `arm | n | sum w | ESS | max w | weight share` — six labels for five columns of something
    else. §18c parsed every code fence in the spec but states it did not run them as modules with
    their audit calls, which is exactly where this surfaces.

    So the two tables are concatenated into one grid, the population header is KEPT as a labelled row
    inside the body, and the short rows are padded. Nothing is dropped and nothing is unlabelled; the
    rendering stays deterministic, so §12.11's byte-identity criterion is unaffected.
    """
    return tuple(row + ("—",) * (width - len(row)) for row in table)


def _record_exclusion(df: pd.DataFrame, in_model: pd.Series, spec: C.Specification,
                      audit: Audit) -> None:
    """The complete-case entry, which must name as many patients as it says it excluded.

    This is the rule `_MUST_NAME_CASES` cannot express for this kind [Stage 6 §7.1] — `model` carries
    four entries and only this one removes patients — stated here where it can be stated exactly. It
    is Stage 5's `_record_removal` (cohort.py:185-206) line for line, both mechanisms included,
    because both are what make the comparison able to fail at all:

      * the identifiers go through a `set`, so a duplicated case_id collapses to one name;
      * a missing case_id is dropped by `notna()` rather than stringified.

    Without either, `len(case_ids)` equals `n` by construction and the check is dead code. That is not
    hypothetical — `Audit.record` normalises with `tuple(sorted(str(c) for c in case_ids))`
    (data.py:271), which sorts but does not deduplicate, and `str(pd.NA)` is the four characters
    `<NA>`. A draft collecting a plain `sorted(...)` list could never raise.

    The raise comes BEFORE `audit.record`, again as Stage 5 does it. Recording first and reading the
    entry back leaves an unreconcilable entry in the Audit when the check fires, and `Audit` has no
    removal path — the log would carry the very inconsistency the raise exists to prevent.
    """
    excluded = df.loc[~in_model & df["case_id"].notna(), "case_id"]
    case_ids = tuple(sorted(set(excluded)))
    n_excluded = int((~in_model).sum())
    if len(case_ids) != n_excluded:
        raise C.SchemaError(
            f"{spec.step('covariate_completeness')} excludes {n_excluded} record(s) and can name "
            f"{len(case_ids)}. "
            "An estimate whose denominator the log cannot reconstruct is an estimate nobody can "
            "check [§11]. A duplicated or missing case_id is forbidden by Stage 2's A2, so this is a "
            "belt over those braces — written here because this is the first stage whose denominator "
            "differs from the cohort's [Stage 6 §4.4].")
    audit.record("model", spec.step("covariate_completeness"), n_excluded,
                 _completeness_detail(df, in_model, spec), case_ids=case_ids,
                 table=_completeness_table(df, in_model, spec))


# --- the private fit and the two public entry points -------------------------------------------------
#
# Six things about the order below, each of which is a failure if moved [Stage 6 §8]:
#
#   * `_assert_fit_inputs` comes FIRST, before the mask and before the design. A frame failing F1 or F2
#     otherwise produces a fit that raises later with a message about rank.
#   * `_record_exclusion` comes BEFORE the design, so the log names the excluded patients even if
#     `design` then raises on D2 or D3. The exclusion is a fact about the frame, established before any
#     modelling choice, and a log that names it is useful precisely when what follows fails.
#   * `e` and `w` are built as all-nan float64 Series and FILLED, never by `reindex` on a shorter
#     Series — which would produce the same nan by accident rather than by construction, and would
#     silently accept a mismatched index.
#   * `a` is bound ONCE, from the complete-case subset, and used for both the fit and the weights. Two
#     separate reads would be two chances to align the weights to a different subset than the fit.
#   * `per_arm` ranges over TREATMENT_LABELS, not over the arms observed. An arm absent here is a
#     FitError from `ess`, which is right: F2 established both arms are present, so an empty one means
#     the mask lost it.
#   * `_assert_probabilities` sits BETWEEN the fit and the weights, because a degenerate `e` becomes a
#     weight of exactly 0.0 one line later and stops being visible.
#
# And one about the shape [Stage 11 §4.3]: `_fit` is private and `fit`/`fit_full` are one line each,
# so BOTH public entry points keep `(df, audit)` with no defaults — `test_propensity.py`'s signature
# guard survives verbatim, parametrised over the two names rather than rewritten, and the door it is
# keeps its hinges.


def _fit(df: pd.DataFrame, spec: C.Specification, audit: Audit) -> Propensity:
    """The propensity score and overlap weights of ONE declared specification, over the [§3] cohort.

    **Private, and the two public names below are its only callers** [Stage 11 §4.1, §4.3]. It takes
    a `C.Specification` and not a covariate list, and it is not part of the public surface, so a
    caller cannot reach a third specification without one being declared in `config.py` with a [§13]
    amendment behind it. A `covariates=` keyword on a public entry point was declined twice — once at
    Stage 6 §6.4 and once here — and the second reason is the decisive one: under
    `covars or C.PS_COVARIATES` a call site that forgets the keyword produces the PRIMARY
    specification, and every number downstream reads as the sensitivity arm's.

    Every audit step name goes through `spec.step`, because `Audit.entry` is first-match
    (data.py:273-275) and two fits over one cohort otherwise record the same four names. With the
    primary's empty suffix nothing landed moves; with the arm's, six entries are disjoint from six
    [Stage 11 §4.4].

    Returns a Propensity carrying `spec`; adds no column to `df` and edits nothing. Appends four
    `model` entries.

    Raises SchemaError on F1/F2 and on D1-D4, and model.FitError on F3/F4/F5/F6 or a fit that does not
    converge. Stage 10 catches the second to drop and count a replicate, and may catch nothing else: a
    SchemaError here is a bug in the resampler, not a sparse replicate.
    """
    _assert_fit_inputs(df)                                     # F1, F2

    in_model = model.complete_cases(df, spec.covariates)
    _record_exclusion(df, in_model, spec, audit)               # §7.2, entry 1

    X, dropped = model.design(df.loc[in_model], spec.covariates)
    audit.record("model", spec.step("design_matrix"), X.shape[1] + 1,
                 _design_detail(X, dropped, spec),
                 table=_design_table(X, dropped, spec))

    a = df.loc[in_model, C.TREATMENT].to_numpy(dtype=float)
    fitted = model.firth(X, a)                                 # raises FitError; no fallback [§5.5]
    _assert_probabilities(fitted, df.loc[in_model, "case_id"])  # F5

    e = pd.Series(np.nan, index=df.index, dtype="float64")   # np.nan, NOT pd.NA — §3.1
    e.loc[in_model] = fitted.p
    w = pd.Series(np.nan, index=df.index, dtype="float64")
    w.loc[in_model] = np.where(a == 1.0, 1.0 - fitted.p, fitted.p)

    audit.record("model", spec.step("propensity_fit"), int(in_model.sum()),
                 _fit_detail(fitted, e, spec), table=_fit_table(fitted))

    per_arm = {code: ess(w[in_model & (df[C.TREATMENT] == code)]) for code in C.TREATMENT_LABELS}
    weights = _weights_table(df, in_model, w, per_arm)
    population = _weighted_population_table(df, in_model, w, spec)
    audit.record("model", spec.step("overlap_weights"), int(in_model.sum()),
                 _weights_detail(spec),
                 table=weights + _padded(population, len(weights[0])))

    return Propensity(e=e, w=w, in_model=in_model, ess=per_arm, fit=fitted, dropped=dropped,
                      spec=spec)


def fit(df: pd.DataFrame, audit: Audit) -> Propensity:
    """The [§7] propensity score and overlap weights over the [§3] cohort. THE estimand's fit.

    Takes no covariate list: the propensity specification is `C.PROPENSITY_PRIMARY` and the estimand
    is indexed by it [§7, Stage 6 §6.4]. `pilots/analysis.py:953` takes `covars or C.PS_COVARIATES`;
    that was declined, because the [§13] propensity-specification rows are *different target
    populations* and a parameter makes them look like options.

    **That refusal carried a promissory note and Stage 11 paid it** — *"when they are un-deferred
    they arrive as a named specification with a [§13] amendment, not as a keyword"*. DECISION 4 (PI,
    2026-08-24) un-deferred exactly one of the six [§13] rows, so there are now TWO named entry
    points over one private `_fit`, each taking `(df, audit)` and neither taking a default. The other
    is `fit_full` below. A third row arrives as a third `C.Specification`, never as an argument here
    [Stage 11 §4].

    Returns a Propensity whose `spec` is `C.PROPENSITY_PRIMARY`; adds no column to `df` and edits
    nothing. Appends four `model` entries under their UNSUFFIXED names, so every landed ledger
    assertion stays byte-identical.

    Raises SchemaError on F1/F2 and on D1-D4, and model.FitError on F3/F4/F5/F6 or a fit that does not
    converge. Stage 10 catches the second to drop and count a replicate, and may catch nothing else: a
    SchemaError here is a bug in the resampler, not a sparse replicate.
    """
    return _fit(df, C.PROPENSITY_PRIMARY, audit)


def fit_full(df: pd.DataFrame, audit: Audit) -> Propensity:
    """The [§13, DECISION 4] full-covariate propensity score and overlap weights. A SENSITIVITY fit.

    `C.PROPENSITY_FULL` adds the four vascular risk factors [§6] declares as balance negative
    controls, which is exactly `C.NEGATIVE_CONTROLS` — so this specification SPENDS them, and the
    balance table's `role` column says so without one row moving [Stage 7 §4.1, Stage 11 §4.4].

    **This is a different target population and not a better-adjusted version of the same one.** The
    ATO estimand is indexed by the propensity score itself: h(X) = e(X){1 − e(X)} moves when the
    specification does, so the number this fit weights toward is not the [§7] estimand measured more
    carefully. `Propensity.spec` is what keeps the two apart downstream, and `bootstrap.run`'s R9
    raises rather than accepting this object — `run`'s replicate body calls `fit` unconditionally, so
    a `run` over this `Propensity` would return intervals whose point estimates are the arm's and
    whose replicates are the primary's specification, with every number finite [Stage 11 §4.5].

    Takes no covariate list, for `fit`'s reason. Appends the same four `model` entries under names
    suffixed `_full_covariate`, so `Audit.entry`'s first-match lookup cannot return the primary's.
    """
    return _fit(df, C.PROPENSITY_FULL, audit)
