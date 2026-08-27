"""Stage 7 — the [§9] balance and overlap diagnostics.

The first stage that returns a **verdict**. Stages 1-5 classify and restrict and Stage 6 fits; every
one of them produces something a later stage consumes. What this stage produces only a *reader*
consumes: [§9]'s threshold, applied to [§6]'s covariates, reported for [§16]. Nothing downstream of
it computes with its output, which is why it is a module and not a function on ``propensity.py``.

**Balance is judged per declared factor LEVEL, not per design-matrix column, and that is why this
module does not import ``model``.** ``model.design`` drops each factor's reference dummy by name and
then drops constant columns [Stage 6 §4.2, §4.3], so a table built from its columns silently omits
``center_HUG``, ``onset_type_witnessed`` and ``center_USZ`` — three of the seven indicators the [§6]
factors expand to on this cohort. ``center = HUG`` carries this cohort's **worst residual
imbalance**, so the shortcut produces a shorter, complete-looking table whose reported worst is
wrong, with nothing raising. What a balance table needs instead — indicators over ``FACTOR_LEVELS``
— is four lines, and ``levels`` is them [Stage 7 §4.2].

This stage **adds no column to the cohort frame, edits no value, refits nothing, re-weights nothing
and computes no second effective sample size** [Stage 6 §11, Stage 7 §0.2]. The temptation here is
different from Stage 6's: the natural way to write a balance table is to attach ``w`` to the frame
and group, and a frame carrying the point fit's weights, resampled into a [§10] replicate, is a
replicate weighted by the wrong score. ``propensity.ess`` is called per centre per arm and a Kish sum
appears nowhere below.

Section references in brackets are to ``statistical_analysis_plan.md``. The specification for this
module is ``specs/stage7_balance_and_overlap.md``; nothing here is invented outside it.

::

    cohort.build(df, audit)         →  93 rows x 33 columns             [Stage 5 §11]
    propensity.fit(cohort, audit)   →  Propensity(e, w, in_model, …)    [Stage 6 §11]
                       │
                       ▼
    ┌────────────────────────────────────────────────────────────────────────────┐
    │  STAGE 7 — balance.py     reads config, data.Audit, propensity              │
    │                           and NOT model — §4.2                              │
    │                                                                             │
    │   smd(x, a, w)                  [§9] one standardised mean difference       │
    │     ├─ weighted means, UNWEIGHTED pooled SD          (§5.1)                 │
    │     └─ four undefined branches, each its own         (§5.3)                 │
    │                                                                             │
    │   assess(df, ps, audit) -> Balance                                          │
    │     ├─ _assert_balance_inputs   B1,B2,B3 THEN B4,B5 → SchemaError  (§4.5)   │
    │     ├─ the balance set, expanded to DECLARED LEVELS  (§4.2)                 │
    │     ├─ SMD before and after, over in_model's 92      (§4.4, §5)             │
    │     ├─ the within-centre overlap table + pooled row  (§6)                   │
    │     └─ two `model` audit entries                     (§7)                   │
    └────────────────────────────────────────────────────────────────────────────┘
                       │
                       ▼
     Balance(covariates, centres, pooled)  →  Stage 11 [§13], Stage 14 [§16]

``assess`` appends two entries of the existing ``model`` kind — no new kind and no new heading, which
is ``data.py:144-149``'s declaration honoured [Stage 7 §7.1]. It writes no file, as no stage before
it does.

``smd`` is public and is Stage 12's as much as this stage's: [§14a]'s support check is a standardised
mean difference over a different population with unit weights, and a second implementation there
would be a second definition of the yardstick [Stage 7 §9].

This module is **not** exempt from the Stage 1 §7 raw-name scan and must never become exempt.
"""
from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass
from typing import Final

import numpy as np
import pandas as pd

import config as C
import propensity
# `_fmt` is private to data.py and is imported anyway, for the reason derive.py and propensity.py
# give: it is the pipeline's *one* float formatter, and a second one is a second way for two runs to
# disagree. Every SMD, range, ESS and weight share below reaches the log through it.
from data import Audit, _fmt


# --- what the stage returns ---------------------------------------------------------------------

@dataclass(frozen=True)
class CovariateBalance:
    """One row of the [§9] balance table: one declared level, or one linear covariate."""

    covariate: str              # "age", or "center = HUG" for a declared factor level
    role: str                   # §4.3 — computed, never declared
    n: int                      # records this covariate is present on, within in_model [§11]
    sd: float                   # the UNWEIGHTED pooled SD both columns are divided by (§5.1)
    unweighted: float           # SMD before weighting; nan where undefined (§5.3)
    weighted: float             # SMD after weighting;  nan where undefined (§5.3)


@dataclass(frozen=True)
class CentreOverlap:
    """One row of the [§9] within-centre overlap table, or the pooled row.

    `status` holds one of §6.2's THREE and nothing else — `_REPORTED`, `_STRUCTURAL` or
    `_NO_CONTRAST`, written here as the constant names rather than as prose so that this comment
    cannot drift from §7.4's literals. Stage 14 reads this column BY VALUE before it draws anything
    (§6.2, §9), so a docstring naming a string the code never produces tells the next reader to
    match on the wrong one. **The pooled row is not a fourth status**: it carries `_REPORTED` like
    any other both-arm row, and what makes it the pooled row is that it lives in `Balance.pooled`
    rather than in `Balance.centres`.
    """

    centre: str
    in_cohort: int              # records at this centre in the [§3] cohort
    weighted: int               # of those, records inside in_model — NOT the same number (§6.1)
    n: dict[int, int]           # keyed by arm code, as TREATMENT_LABELS is
    e_range: dict[int, tuple[float, float]]
    ess: dict[int, float]       # propensity.ess per arm; nan where an arm is empty (§6.3)
    max_weight: float
    weight_share: float
    worst: float                # worst within-centre |SMD|; nan where none is defined
    status: str                 # §6.2's three, verbatim


# The two verdicts are module-level functions over a row sequence, and the dataclass delegates.
# §7's `_smd_detail` needs both, and it has only the rows — the centres do not exist yet when the
# first audit entry is recorded (§8). An earlier draft reached them by constructing a throwaway
# `Balance(tuple(rows), ())`, which broke the moment `Balance` gained a third field and — worse
# while it worked — wrote the threshold predicate a SECOND time, inline, beside `unbalanced()`'s.
# Two spellings of "which rows reach 0.10" in one module is one edit away from a `detail` string
# that names four rows above a table that flags five. One definition, two callers, no temporary.


def _worst(rows: Sequence[CovariateBalance]) -> tuple[CovariateBalance | None, tuple[str, ...]]:
    """The largest |weighted| SMD among the DEFINED rows, AND the undefined ones' names.

    Both, from one call, and that is the whole design of this pair [§3.1]. `max` over a column
    containing `nan` skips it in pandas and poisons it in numpy, so a caller asking only for the
    worst gets either a number that ignores the rows nobody could judge or a `nan` that hides the
    rows anybody could. Returning the pair makes it impossible to report the first without being
    handed the second.
    """
    defined = [row for row in rows if np.isfinite(row.weighted)]
    undefined = tuple(row.covariate for row in rows if not np.isfinite(row.weighted))
    return (max(defined, key=lambda row: abs(row.weighted)) if defined else None), undefined


def _over_threshold(rows: Sequence[CovariateBalance]) -> tuple[str, ...]:
    """Covariates whose |weighted| SMD reaches SMD_THRESHOLD. Undefined rows are NOT here.

    The undefined rows are `_worst`'s second element, so the two collections partition the table
    with the defined-and-balanced ones and no row is in neither by accident.
    """
    return tuple(row.covariate for row in rows
                 if np.isfinite(row.weighted) and abs(row.weighted) >= C.SMD_THRESHOLD)


@dataclass(frozen=True)
class Balance:
    """What Stage 7 returns. The frame and the Propensity come back untouched [§0.2].

    `pooled` is a field and not the last element of `centres`, and that is a decision. An earlier
    draft returned `centres + (pooled,)` as one tuple, which cost three things at once:
    `len(centres)` stopped being `len(CENTER_ORDER)`, so no consumer could iterate the declared
    centres without knowing to drop a row; the pooled row needed a status of its own — a fourth
    literal §6.2's three-status table never declared and Stage 14 would have had to match on; and a
    caller counting reported centres off `centres` would have counted the whole cohort as one of
    them. The split costs nothing, because the RENDERED table is still one table: §8 passes
    `centres + (pooled,)` to `_overlap_table`, so [§9]'s "within each centre as well as pooled" is
    one grid computed by one function, and only the return value distinguishes the two kinds of row.
    """

    covariates: tuple[CovariateBalance, ...]
    centres: tuple[CentreOverlap, ...]   # exactly CENTER_ORDER, in CENTER_ORDER's order
    pooled: CentreOverlap                # the `all (pooled)` row — its own field, never a centre

    def worst(self) -> tuple[CovariateBalance | None, tuple[str, ...]]:
        """§3's pair, over this table's rows. The public name; `_worst` is the definition."""
        return _worst(self.covariates)

    def unbalanced(self) -> tuple[str, ...]:
        """The rows reaching SMD_THRESHOLD, over this table's rows. Undefined rows are not here."""
        return _over_threshold(self.covariates)


# --- the standardised mean difference [§9] -----------------------------------------------------
#
#                        mean_w(x | treated)  −  mean_w(x | control)
#     SMD(x, a, w)  =  ───────────────────────────────────────────────
#                         √( ( var(x | treated) + var(x | control) ) / 2 )
#
# The NUMERATOR is weighted; the DENOMINATOR is not, and it is the same denominator in both columns
# [§9]. Sample variances, ddof = 1. The unweighted column is this expression with w ≡ 1 and the
# weighted column is the same expression with the [§7] overlap weights — ONE function, called twice,
# with only the weights differing, because the alternative is two functions that can disagree about
# the denominator, which is precisely what [§9]'s "so the yardstick does not move" forbids.
#
# Three functions rather than one, and the split is §5.5's: `_pooled_sd` is the yardstick, `_ratio`
# applies a GIVEN yardstick to a numerator, and `smd` is [§9]'s definition — the two over the same
# records. §6.1's within-centre column has to divide a centre's means by the WHOLE weighted set's SD,
# so the separation exists whether or not a function name admits it; hiding it inside `smd` did not
# prevent it, it only made it the one thing in this stage nobody could assert.


def _pooled_sd(x: np.ndarray, a: np.ndarray) -> float:
    """[§9]'s UNWEIGHTED pooled SD over the records given. nan where there is none. §5.1.

    The one definition of the yardstick. Sample variances, ddof = 1 — which §16b's oracle asserts
    rather than assumes, because a ddof = 0 convention differs by sqrt(n/(n-1)) and at n = 92 that
    is 0.5%, visible at six significant figures and invisible to the eye.

    An arm of fewer than two records has no sample variance, so there is no pooled SD: that is
    §5.3's branch 1 and it is a fact about the DENOMINATOR's population, which is why it lives here
    rather than in `_ratio`.
    """
    t, c = a == 1.0, a == 0.0
    if int(t.sum()) < 2 or int(c.sum()) < 2:
        return np.nan
    sd = float(np.sqrt((x[t].var(ddof=1) + x[c].var(ddof=1)) / 2.0))
    return sd if np.isfinite(sd) else np.nan


def _ratio(x: np.ndarray, a: np.ndarray, w: np.ndarray, sd: float) -> float:
    """Weighted mean difference over a GIVEN yardstick. §5.1, and §6.1's foreign-SD case.

    Returns nan — never 0.0 — wherever the quantity is undefined, and §5.3 is the taxonomy.
    `pilots/analysis.py:194-208` returns 0.0 in one of those cases, and that case is a covariate
    which PERFECTLY SEPARATES the arms reported as perfectly balanced [§5.3a].

    **The two-record floor is a floor and no longer a variance requirement.** With `sd` supplied
    from another population, a weighted mean over one record is arithmetically fine — and a [§9]
    diagnostic standardising a single observation is not a diagnostic. So the guard stays, for a
    reason that changed: §5.3 branch 1 was about `var(ddof=1)`, this is about what a number means.
    Measured (§18): it is what keeps every within-centre row of the fixture cohort — one control per
    centre — reported as `missing` rather than as a number nobody should read.
    """
    t, c = a == 1.0, a == 0.0

    if int(t.sum()) < 2 or int(c.sum()) < 2:
        return np.nan                                  # §5.3, first — the FLOOR
    if not w[t].sum() or not w[c].sum():
        return np.nan                                  # §5.3, second — BEFORE np.average [§3.1]
    if not np.isfinite(sd):
        return np.nan                                  # §5.3, third
    if sd == 0.0:
        # §5.3a — the ARM CONSTANTS, never the weighted means. Both arm variances are zero whenever
        # the pooled SD is (they are non-negative and they sum to zero), so each arm is constant and
        # its constant is the raw datum x[t][0], with no summation, no weights and no rounding in it.
        # Comparing the two WEIGHTED means instead fails at 1 ulp — measured on the committed
        # fixture, where `nihss_baseline` is 14.0 on every record of both arms and np.average gives
        # 14.0 against 14.000000000000002 — so a covariate with no imbalance at all was reported as
        # undefined. A tolerance is not the repair: it would need a scale, and it would make this
        # branch fire on covariates that are merely NEARLY constant, which is its opposite.
        return 0.0 if x[t][0] == x[c][0] else np.nan   # §5.3, fourth

    return float((np.average(x[t], weights=w[t]) - np.average(x[c], weights=w[c])) / sd)


def smd(x: pd.Series | np.ndarray, a: pd.Series | np.ndarray,
        w: pd.Series | np.ndarray) -> float:
    """The [§9] standardised mean difference: weighted means over an UNWEIGHTED pooled SD.

    Numerator and denominator over the SAME records — [§9]'s definition, and the only form Stage 12
    needs [§9, §14a]. Its signature is unchanged by §5.5's split and must stay unchanged: a `sd=`
    keyword here would make [§9]'s prescribed denominator look like a caller's option, which is
    Stage 6 §6.4's argument and §15's.

    Positional, not by keyword, and aligned by position rather than by index: the caller has
    already masked all three to the same records, and B1 is what makes that safe.
    """
    x = np.asarray(x, dtype=float)
    a = np.asarray(a, dtype=float)
    w = np.asarray(w, dtype=float)
    return _ratio(x, a, w, _pooled_sd(x, a))


# --- the balance set [§6, §9] -------------------------------------------------------------------
#
# Balance is judged against the FULL [§6] confounder set, not only the covariates a given
# specification put in its propensity model [§9]. The set is BALANCE_SET — 14 declared names —
# expanded to one row per DECLARED LEVEL of each factor and one per linear covariate: 19 rows.
#
# The set does not change when the specification does, and that is the whole purpose of the `role`
# column: under [§13]'s full-covariate propensity model the four vascular risk factors move from
# `negative control` to `propensity model` and NOT ONE ROW moves in or out of the table, so the two
# specifications are read against the same fourteen names.


def levels(df: pd.DataFrame, name: str) -> list[tuple[str, pd.Series]]:
    """One (label, indicator) pair per DECLARED level of a factor; one pair for a linear covariate.

    **Public, and not for symmetry: [§14a]'s support-check baseline table needs it and may not own a
    second copy** [Stage 12 §10.3]. Stage 7 §5.5 made `smd` public on the same argument -- a second
    implementation of `smd` would be a second definition of the yardstick, and a second
    implementation of this would be a second definition of WHAT A ROW IS.

    **And the failure of a second copy would be invisible.** Stage 12 §10.1 measured that all three
    of `onset_type`'s declared levels are present in [§14a]'s treated arm, so a copy that ranged over
    OBSERVED levels would produce a byte-identical table on v7 and silently drop a row on the next
    workbook. The rename to a public name is the whole of Stage 12's change to this module.

    Ranges over FACTOR_LEVELS, never over the frame and never over a design matrix. A level nobody
    in the cohort has is a row reading 0, which is a fact about the population; a level that
    vanishes is a fact nobody sees [Stage 6 §7.5].

    The indicator carries the factor's own missingness — `mask(isna())` — rather than encoding an
    absent factor value as a zero in every level, which is Stage 6 §4.5 D4's failure in a table
    instead of in a design. That mask is **unreachable from `assess`** while every `CATEGORICAL`
    name is a PS covariate, because such a record is complete-cased out of `in_model` first
    (measured, §12.2); it stays because it is correct, because Stage 12 will need it over a
    population `complete_cases` did not build, and because a guard costing one method call is not
    worth removing to make a test honest.
    """
    if name not in C.CATEGORICAL:
        return [(name, df[name].astype("Float64").astype(float))]
    return [(f"{name} = {level}",
             (df[name] == level).astype(float).mask(df[name].isna()))
            for level in C.FACTOR_LEVELS[name]]


def _role(name: str) -> str:
    """What this covariate is TO the specification being diagnosed. §4.3.

    Reads PS_COVARIATES for the same reason `propensity.fit` does [Stage 6 §6.4]: the estimand is
    indexed by the propensity model, so "in the model" is a property of the one prespecified
    specification and not an argument a caller may vary.

    `BALANCE_ONLY` is NOT the negative-control set: it holds five names and [§6] names four vascular
    risk factors as negative controls. The fifth is `penumbra_ml`, excluded because it is a
    deterministic function of `core_ml` and `tmax6_ml`, and a table equating the two sets asserts
    that a collinear covariate is a covariate weighting was never expected to fix.
    """
    if name in C.PS_COVARIATES:
        return "propensity model"
    if name in C.NEGATIVE_CONTROLS:
        return "negative control"
    return "excluded [§6]"


def _table(sub: pd.DataFrame, a: np.ndarray, w: np.ndarray, names: Sequence[str],
           sds: dict[str, float] | None = None) -> tuple[CovariateBalance, ...]:
    """One CovariateBalance per declared level, over the records `sub` already restricts to.

    The per-row mask is the covariate's own: [§11] is complete-case PER ESTIMATE, and a row of this
    table is an estimate. `n` is what survives it, which is why the column exists (§4.4), and
    §12.2a is what tells this form from one mask over the frame — on every frame this stage has, the
    column is a constant, so a global mask would be green.

    `sds` is §5.5's yardstick. None means "compute each row's own", which is right for the table
    over `in_model` — that IS the population the yardstick is defined on. A dict means "use these",
    which is every within-centre call: the numerator is the centre's and the denominator is not. It
    is a parameter with a default rather than two functions, because the two calls differ in one
    argument and nothing else; and it is NOT on `smd`, whose signature stays [§9]'s (§5.2).

    `sds.get(label, np.nan)` and not `sds[label]`: a label the pooled table did not produce has no
    yardstick, and the honest answer for its row is `missing` rather than a `KeyError` inside a
    diagnostic. It cannot happen today — every within-centre call passes a subset of `BALANCE_SET` —
    and §12.6 asserts that, so the `get` stays inert while the assertion is what would notice.
    """
    rows: list[CovariateBalance] = []
    for name in names:
        for label, indicator in levels(sub, name):
            values = indicator.to_numpy(dtype=float)
            present = np.isfinite(values)
            x, arm, weight = values[present], a[present], w[present]
            sd = _pooled_sd(x, arm) if sds is None else sds.get(label, np.nan)
            rows.append(CovariateBalance(
                covariate=label, role=_role(name), n=int(present.sum()), sd=sd,
                unweighted=_ratio(x, arm, np.ones_like(weight), sd),
                weighted=_ratio(x, arm, weight, sd)))
    return tuple(rows)


# --- overlap [§9] ---------------------------------------------------------------------------------
#
# [§9] asks for overlap "within each centre as well as pooled", and the pooled report is the
# `all (pooled)` ROW of that same table rather than a second table: Stage 6's `overlap_weights` entry
# already carries the pooled arm cells, so a third rendering would be the duplication Stage 6 §7.4
# declined, and a pooled row computed by the SAME CODE as the centre rows cannot disagree with them
# about what a column means. In the return value it is a field of its own (§3), because the log's
# reader wants one table with a labelled row and a Python caller wants `centres` to mean the declared
# centres and nothing else.

_REPORTED: Final[str] = "overlap reported"
_STRUCTURAL: Final[str] = "structural non-positivity [§3]"
_NO_CONTRAST: Final[str] = "no weighted contrast"


def _status(at: pd.Series, arms: dict[int, np.ndarray]) -> str:
    """§6.2's three, and the order is the specification. Called for the pooled row too.

    `at` is COHORT membership and `arms` is what survived `in_model`, so the two branches are
    different facts: a declared centre with no cohort record was removed by [§3] restriction 1,
    and a centre with records but an empty arm was emptied by complete-casing. Reading the second
    as the first would report a data-handling consequence as a design one.

    There is no pooled branch and no fourth constant. The pooled row's `at` is True everywhere and
    B5 has already established both arms, so it reaches `_REPORTED` — which is true of it. What
    marks it as pooled is `Balance.pooled` and its `centre` cell, never this column (§3, §6.2).
    """
    if not int(at.sum()):
        return _STRUCTURAL
    if any(not int(mask.sum()) for mask in arms.values()):
        return _NO_CONTRAST
    return _REPORTED


def _centre(df: pd.DataFrame, ps: propensity.Propensity, centre: str | None,
            sds: dict[str, float]) -> CentreOverlap:
    """One row of §6.1's table. `centre=None` builds the `all (pooled)` row from the same code.

    **`propensity.ess` is called only where the arm is non-empty**, and that is a decision rather
    than defensive padding (§6.3). `ess` raises `model.FitError` on an empty arm with a message
    reading "the in_model mask lost one" (Stage 6 §6.2) — correct for the [§7] fit, where an empty
    arm means something went wrong, and wrong here, where an empty arm within one centre is a
    FINDING TO REPORT. Stage 5's P2 guarantees both arms at every retained centre; nothing
    guarantees both arms at every retained centre AFTER complete-casing, and on a workbook where one
    centre's only control is missing a CTP volume the diagnostic would abort mid-table rather than
    saying so. Stage 6's raise stays the backstop for a caller that passes an empty arm anyway.

    `sds` is §5.5's one yardstick and is REQUIRED, with no default: a centre row computed against
    its own variances is the defect §5.5 records — up to 46% out, every cell finite and plausible —
    and a defaulted `None` here is the one line that would reintroduce it silently.
    """
    at = pd.Series(True, index=df.index) if centre is None else df["center"] == centre
    weighted = at & ps.in_model
    sub = df.loc[weighted]
    e, w = ps.e.loc[weighted], ps.w.loc[weighted]
    arms = {code: (sub[C.TREATMENT] == code).to_numpy() for code in C.TREATMENT_LABELS}

    worst, defined = np.nan, [
        row.weighted for row in _table(
            sub, sub[C.TREATMENT].to_numpy(dtype=float), w.to_numpy(dtype=float),
            # `center` is dropped INSIDE a centre: its indicators are constant there by
            # construction, so those rows can only ever read 0.0 and a row that cannot vary is
            # not evidence. Measured: dropping them changes no centre's worst (§18).
            [c for c in C.BALANCE_SET if c != "center"] if centre is not None else C.BALANCE_SET,
            sds)                                    # §5.5 — the pooled yardstick, not the centre's
        if np.isfinite(row.weighted)]
    if defined:
        worst = float(max(defined, key=abs))

    return CentreOverlap(
        centre="all (pooled)" if centre is None else centre,
        in_cohort=int(at.sum()), weighted=int(weighted.sum()),
        n={code: int(mask.sum()) for code, mask in arms.items()},
        e_range={code: ((float(e[mask].min()), float(e[mask].max())) if mask.any()
                        else (np.nan, np.nan)) for code, mask in arms.items()},
        ess={code: (propensity.ess(w[mask]) if mask.any() else np.nan)
             for code, mask in arms.items()},
        max_weight=float(w.max()) if len(w) else np.nan,
        weight_share=float(w.sum() / ps.w[ps.in_model].sum()) if len(w) else np.nan,
        worst=worst,
        status=_status(at, arms))


# --- the two `model` audit entries [§7] -------------------------------------------------------------
#
# No new kind, and that is `data.py:144-149`'s declaration honoured: "Stages 7-13 all render under
# `model` or under an existing kind" [Stage 7 §7.1]. A balance table is not obviously a "fitted
# model", so the argument is made rather than assumed: `model`'s entries are CLAIMS ABOUT A FIT — the
# population it ran on, the matrix it used, the coefficients it produced, the weights it implies —
# and a standardised mean difference is a claim about the same fit, namely how well its weights did
# what they were for. The log then reads: this is what was fitted, these are the weights, this is
# what the weights achieved. A separate heading would put the verdict in a different section from
# the thing it judges.
#
# Four rules about the tables, all inherited:
#
#   * Every declared thing is rendered whether or not the data fills it — every level of every
#     factor including the reference and one nobody has, every centre in CENTER_ORDER including one
#     the [§3] restriction removed. `pd.crosstab` is used nowhere in this repository for that reason.
#   * Every cell is computed from the object it describes, never by subtracting another row [Stage 5
#     §7.3]. The pooled row is computed over `in_model`, not by combining the centre rows.
#   * `_fmt` is the only float formatter and `missing` the only rendering of an absent value, so an
#     undefined SMD renders as `missing` and never as `nan`, `0`, or a numpy repr.
#   * The verdict column is a STRING, never a bool: `_fmt(True)` is `float(True)` is "1", which would
#     read as a number beside an SMD.
#
# And one rule this stage adds, measured: **no cell of either table may contain a `|`.**
# `data._md_table` does no escaping (data.py:242-245) — it joins cells with " | " and sizes the
# separator from `len(rows[0])` — so one pipe inside a header cell gives a header row with more
# markdown cells than its own separator and body. Measured before the columns were renamed:
# `balance_smd` rendered 8 cells against a 6-cell separator and `overlap_by_centre` 15 against 13, in
# every log this stage writes, and byte-identity across hash seeds passed the whole time. Hence
# `abs SMD < 0.1` and `worst abs SMD`.


def _verdict(value: float) -> str:
    """`yes`, `no` or `undefined` — a STRING, because `_fmt(True)` is `float(True)` is "1" [§7.2].

    Three-valued and not two: a row whose SMD is undefined has neither passed the threshold nor
    failed it, and rendering it `no` would put a covariate nobody could judge among the findings.
    """
    if not np.isfinite(value):
        return "undefined"
    return "yes" if abs(value) < C.SMD_THRESHOLD else "no"


def _smd_table(rows: Sequence[CovariateBalance]) -> tuple[tuple[str, ...], ...]:
    """One row per declared level, in BALANCE_SET order. §7.2.

    The threshold is interpolated into the header from `C.SMD_THRESHOLD` rather than written as
    `0.10`, so the column cannot claim a threshold the verdict was not computed against. It renders
    as `0.1`, because `_fmt` is six significant figures — the header says what the code compares.
    """
    header = ("covariate", "role", "n", "SMD before", "SMD after",
              f"abs SMD < {_fmt(C.SMD_THRESHOLD)}")
    return (header, *(
        (row.covariate, row.role, str(row.n), _fmt(row.unweighted), _fmt(row.weighted),
         _verdict(row.weighted))
        for row in rows))


def _range(bounds: tuple[float, float]) -> str:
    """The propensity range as ONE cell, both numbers through `_fmt`. §7.2.

    One cell rather than two columns per arm, which would take this table to sixteen. `missing`
    rather than `missing–missing` when the arm is empty, so an absent range reads as the single
    fact it is.

    **An arm of ONE renders as `0.257067–0.257067`, and that is deliberate.** Measured on the fixture
    cohort, where every centre has exactly one control. Collapsing it to a single number would make
    the cell's width carry information the `n` column already carries, and would make an arm of one
    indistinguishable from an arm of many whose propensities happen to coincide — which is a real
    state and a different one. A degenerate range that LOOKS degenerate is the honest rendering; the
    reader learns the arm size from the column two cells to the left [§9].
    """
    low, high = bounds
    return "missing" if not np.isfinite(low) else f"{_fmt(low)}–{_fmt(high)}"


def _overlap_table(rows: Sequence[CentreOverlap]) -> tuple[tuple[str, ...], ...]:
    """One row per declared centre, then `all (pooled)`. §6.1, §7.2.

    `worst abs SMD` renders the absolute value the column names, while `CentreOverlap.worst` keeps
    the sign for a caller that wants to know which arm it favours. It is standardised by the POOLED
    in_model SD, not by the centre's own — §5.5, and the column name says nothing about that because
    a column name cannot; `_overlap_detail` does.
    """
    labels = C.TREATMENT_LABELS
    header = ("centre", "in cohort", "weighted",
              *(f"n ({label})" for label in labels.values()),
              *(f"e ({label})" for label in labels.values()),
              *(f"ESS ({label})" for label in labels.values()),
              "max w", "weight share", "worst abs SMD", "status")
    return (header, *(
        (row.centre, str(row.in_cohort), str(row.weighted),
         *(str(row.n[code]) for code in labels),
         *(_range(row.e_range[code]) for code in labels),
         *(_fmt(row.ess[code]) for code in labels),
         _fmt(row.max_weight), _fmt(row.weight_share), _fmt(abs(row.worst)), row.status)
        for row in rows))


def _smd_detail(df: pd.DataFrame, ps: propensity.Propensity,
                rows: Sequence[CovariateBalance]) -> str:
    worst, undefined = _worst(rows)              # §3 — never a throwaway Balance
    over = _over_threshold(rows)                 # §3 — the SAME predicate `unbalanced()` uses
    return (
        f"[§9] standardised mean differences over {len(rows)} declared level(s) of "
        f"{len(C.BALANCE_SET)} [§6] covariate(s), before and after the [§7] overlap weights, on "
        f"{int(ps.in_model.sum())} of {len(df)} cohort record(s) — the weighted set, in BOTH "
        "columns, so the pooled standard deviation is one yardstick over one population "
        "[Stage 7 §4.4]. The denominator is the UNWEIGHTED pooled SD in both, which is what makes "
        "the two columns comparable [§9]. Balance is judged against the FULL [§6] confounder set "
        f"regardless of what this specification fitted, and the {len(C.NEGATIVE_CONTROLS)} vascular "
        "risk factors are "
        "negative controls: residual imbalance on them shows what weighting does not fix [§6]. "
        f"{len(over)} row(s) reach |SMD| = {_fmt(C.SMD_THRESHOLD)} after weighting"
        + (": " + ", ".join(over) if over else "")
        + (f"; worst {_fmt(abs(worst.weighted))} on {worst.covariate}" if worst else "")
        # The trailing full stop on the undefined branch is load-bearing and the spec's own fence
        # omitted it: without it the sentence runs into the next one — "…, penumbra_ml A row's `n` is
        # its own denominator" — which is only visible by reading a rendered `detail` as prose, on a
        # frame that HAS an undefined row. Neither review round had one (§21.5).
        + (f". {len(undefined)} row(s) are undefined and are reported as such rather than as zero: "
           + ", ".join(undefined) + "." if undefined else ". No row is undefined.")
        + " A row's `n` is its own denominator [§11]; the record(s) outside the fit are named in "
        "`covariate_completeness` above. This entry reports balance and does not judge the "
        "estimator: [§7] prescribes one propensity model and there is no second one to try.")


def _overlap_detail(centres: Sequence[CentreOverlap]) -> str:
    absent = [c.centre for c in centres if c.status == _STRUCTURAL]
    thin = [c.centre for c in centres if c.status == _NO_CONTRAST]
    return (
        f"[§9] overlap within each centre as well as pooled, over {len(C.CENTER_ORDER)} declared "
        "centre(s). A pooled distribution can look acceptable while treatment is nearly determined "
        "by centre, and only the within-centre view distinguishes patient-level equipoise from "
        f"that. {len(absent)} centre(s) contribute no cohort record and are reported as "
        "structurally non-positive rather than given an overlap plot [§9]"
        + (": " + ", ".join(absent) if absent else "")
        + (f". {len(thin)} centre(s) lost an arm to complete-casing and carry no weighted contrast: "
           + ", ".join(thin) if thin else "")
        + ". `restrict_centres` above names which patients each [§3] restriction removed. The "
        "pooled row is computed by the same code as the centre rows and reconciles with "
        "`overlap_weights`; no second effective sample size is computed anywhere in this stage "
        "[Stage 6 §11]. Every `worst abs SMD` in this table is standardised by the pooled unweighted "
        "standard deviation over the WEIGHTED SET, not by the centre's own, so the column is "
        "comparable across these rows and against `balance_smd` above [Stage 7 §5.5]. It is still "
        "not a ranking: an abs SMD over an arm of two carries almost no information, which is why "
        "the arm sizes are in the same row [Stage 7 §13].")


# --- the preconditions [§4.5] -----------------------------------------------------------------------
#
# Five checks in TWO PHASES, each phase collected and raised together — and that departure from every
# other assertion helper in this repository is load-bearing rather than stylistic (§4.5a).
#
# Stages 2, 5 and 6 collect their failures and raise once, which is right wherever the checks are
# independent. Here they are not: B4, B5 and everything in `assess` read through the frame and the
# mask that B1, B2 and B3 are about, and pandas does not return a wrong answer for a broken mask or
# an absent column — it RAISES, from inside the collection, before the SchemaError is ever assembled.
# Measured on pandas 2.3.3, all three landing on B5's `df.loc[ps.in_model, C.TREATMENT]`:
#
#     Propensity reindexed to a shifted index  →  AssertionError('')   BARE, no message at all
#     in_model cast to object with one <NA>    →  ValueError("Cannot mask with non-boolean array…")
#     the TREATMENT column dropped             →  KeyError('ivt')
#
# In every case the assembled message — B1's paragraph about a finite table for a population that
# does not exist, B2's about a three-valued mask deciding the denominator, B3's naming the column —
# was thrown away by an exception raised a check or two later. So the boundary is not "mask checks"
# against "the rest": it is CAN THIS BE READ against IS THE DATA JUDGEABLE, and a column-presence
# check belongs to the first for exactly the reason a mask check does. A reader who tidies this back
# into one collection reintroduces the defect; §12.11 carries a companion for each of the three.


def _assert_balance_inputs(df: pd.DataFrame, ps: propensity.Propensity) -> None:
    # PHASE 1 — can this frame and this mask be READ? B1, B2, B3.
    bad: list[str] = []

    misaligned = [name for name, s in (("e", ps.e), ("w", ps.w), ("in_model", ps.in_model))
                  if not s.index.equals(df.index)]
    if misaligned:
        bad.append(
            f"B1  the Propensity is not aligned to this frame: {', '.join(misaligned)} "
            f"carr{'ies' if len(misaligned) == 1 else 'y'} a different index. "
            f"{len(ps.in_model)} mask row(s) against {len(df)} frame row(s), sharing "
            f"{len(ps.in_model.index.intersection(df.index))} index label(s). e, w and in_model are "
            "Series on the COHORT's index [Stage 6 §0.2], and a positional read of a misaligned pair "
            "produces a balance table for a population that does not exist, with every cell finite.")

    if ps.in_model.dtype != bool or ps.in_model.isna().any():
        bad.append(
            f"B2  in_model is {ps.in_model.dtype} and carries "
            f"{int(ps.in_model.isna().sum())} missing value(s). It is boolean and TOTAL by "
            "construction [Stage 6 §4.4]; a three-valued mask resolves by branch order at the "
            "first `if`, and here that decides who is in the denominator.")

    absent = [c for c in (*C.BALANCE_SET, C.TREATMENT) if c not in df.columns]
    if absent:
        bad.append(
            f"B3  {', '.join(absent)}: not a column of the frame. [§9] judges balance against the "
            "full [§6] confounder set, so a name this frame cannot supply is a balance table that "
            f"silently judges less than it claims — and {C.TREATMENT} is checked with them because "
            "B5 and every row of §5 read the arm from it.")

    if bad:
        raise C.SchemaError(
            "\n  " + "\n  ".join(bad)
            + f"\n  {len(bad)} balance assertion(s) failed over {len(df)} records. The frame or the "
            "mask cannot be read, so B4 and B5 were not run: they read through both (§4.5a).")

    # PHASE 2 — the frame and the mask are readable. These two describe the DATA, and they are
    # collected as Stages 2, 5 and 6 collect theirs.
    bad = []

    denied = [c for c in C.BALANCE_SET if c in C.POST_TIME_ZERO]
    if denied:
        bad.append(
            f"B4  {', '.join(denied)}: post-time-zero [invariant 4, §12]. Balance is a property of "
            "confounders measured at or before time zero; a post-exposure variable judged against "
            "[§9]'s threshold reads as a confounder that weighting failed to fix, when it is a "
            "variable no weighting should touch. onset_to_groin_min is reported by arm [§12].")

    present = [code for code in C.TREATMENT_LABELS
               if int((df.loc[ps.in_model, C.TREATMENT] == code).sum())]
    if len(present) < 2:
        bad.append(
            "B5  " + ", ".join(
                f"{C.TREATMENT_LABELS[c]}: "
                f"{int((df.loc[ps.in_model, C.TREATMENT] == c).sum())}"
                for c in C.TREATMENT_LABELS)
            + ". A standardised mean difference needs both arms, so with one every row of this "
            "table is undefined and the table is not a diagnostic. Stage 6's F2 established both "
            "arms were present before the fit; an arm missing HERE was lost by in_model.")

    if bad:
        raise C.SchemaError(
            "\n  " + "\n  ".join(bad)
            + f"\n  {len(bad)} balance assertion(s) failed over {len(df)} records, "
            f"{int(ps.in_model.sum())} of them weighted.")


# --- the public entry point ---------------------------------------------------------------------------
#
# Five things about the order below, each of which is a failure if moved [§8]:
#
#   * `_assert_balance_inputs` comes FIRST, before the mask and before any arithmetic — Stage 5 §4.4's
#     rule. A misaligned Propensity otherwise produces a complete table for a population that never
#     existed, and §4.5's phase 1 is what makes that a message rather than a bare AssertionError.
#   * `sub`, `a` and `w` are bound ONCE, from `in_model`, and every row of every table is computed
#     from those three. Two separate `.loc` reads would be two chances to align the weights to a
#     different subset than the means. This is also where Stage 6 §9's rule is honoured — range over
#     `in_model`, never over `notna()`, and never fill — and nothing else below reads `ps.w`.
#   * `_table` is called by both entries with the same signature and the same arguments, so the
#     pooled worst |SMD| and the balance table's worst are the same number from the same call. It is
#     computed TWICE — once here, once inside `_centre(df, ps, None, sds)` — and that is the price of
#     the one-code-path property, priced at 7.4 ms. Not "computed once": sharing the result would
#     make §12.8's reconciliation unfalsifiable, and a table that agrees with itself by construction
#     has stopped being evidence [Stage 5 §7.3].
#   * `_centre(df, ps, None, sds)` builds the pooled row, so the pooled and per-centre cells cannot
#     mean different things (§6.1). It goes into `Balance.pooled` rather than into `centres`, so
#     `centres` is exactly CENTER_ORDER (§3).
#   * `sds` is harvested from `covariates`, after the first entry is recorded and before any centre
#     row. That ordering is the whole of §5.5's implementation: the yardstick exists as a value only
#     because `CovariateBalance` carries `sd`, and it is read from the table computed over `in_model`
#     — never recomputed, so the pooled row and the centre rows cannot disagree about it.


def assess(df: pd.DataFrame, ps: propensity.Propensity, audit: Audit) -> Balance:
    """The [§9] balance and overlap diagnostics over the [§3] cohort and its [§7] weights.

    Takes no covariate list and no threshold: the balance set is BALANCE_SET and the threshold is
    SMD_THRESHOLD, both prespecified [§9]. Adds no column to `df`, edits nothing, refits nothing
    and computes no second effective sample size [Stage 6 §11].

    Raises SchemaError on B1-B5 and on nothing else. Imbalance is a finding, not a failure [§5.4]:
    residual imbalance is a question for the PI under [§13] and there is nothing for a caller to
    catch. This is the one place a reader might expect Stage 6's posture and not get it.
    """
    _assert_balance_inputs(df, ps)                             # B1, B2, B3 then B4, B5 — §4.5

    sub = df.loc[ps.in_model]
    a = sub[C.TREATMENT].to_numpy(dtype=float)
    w = ps.w.loc[ps.in_model].to_numpy(dtype=float)

    covariates = _table(sub, a, w, C.BALANCE_SET)
    audit.record("model", "balance_smd", int(ps.in_model.sum()),
                 _smd_detail(df, ps, covariates), table=_smd_table(covariates))

    # §5.5 — the ONE yardstick, harvested from the table computed over in_model and handed to every
    # centre row. A centre's numerator is its own; its denominator is the whole weighted set's.
    sds = {row.covariate: row.sd for row in covariates}
    centres = tuple(_centre(df, ps, centre, sds) for centre in C.CENTER_ORDER)
    pooled = _centre(df, ps, None, sds)                        # the `all (pooled)` row — §6.1
    audit.record("model", "overlap_by_centre",
                 sum(1 for c in centres if c.status == _REPORTED),
                 _overlap_detail(centres), table=_overlap_table(centres + (pooled,)))

    return Balance(covariates=covariates, centres=centres, pooled=pooled)
