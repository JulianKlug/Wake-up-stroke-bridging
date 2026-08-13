"""Stage 3 — derived variables [§5, §6, §13].

Every variable the analysis needs that the workbook does not contain, built from the registries in
``config.py`` rather than from literals, with missingness reimposed everywhere a comparison would
otherwise manufacture a zero.

This stage **adds columns**. It drops no row, drops no column, and edits no value that Stage 2
delivered. If a value looks like it needs fixing, the fix belongs in Stage 2's ``_correct`` with its
own audit entry — corrections are that stage's, and its log is the record of them.

Section references in brackets are to ``statistical_analysis_plan.md``. DECISION *n* refers to the
decisions recorded in ``../out/stage0_data_inventory.md``. The specification for this module is
``specs/stage3_derived_variables.md``; nothing here is invented outside it.

Two entry points, called at two different points of the pipeline, because the five things this stage
builds are not all the same kind of thing::

    data.load()  →  (df, audit)          25 analysis columns, 126 rows   [Stage 2 §11]
                         │
                         ▼
    ┌──────────────────────────────────────────────────────────────────────┐
    │  derive(df, audit)               row-wise: each record's own values   │
    │    ├─ _assert_onset_flags(df)    both flags positive → SchemaError    │
    │    ├─ onset_type                 three levels, <NA> if a flag is <NA> │
    │    ├─ unknown_onset              [§13], from onset_type               │
    │    ├─ _dichotomies(df)           four, through OUTCOMES and OPS       │
    │    └─ absence_by_column(…)       one row per derived column           │
    └──────────────────────────────────────────────────────────────────────┘
                         │
                         ▼
                31 columns  →  Stage 4 eligibility  →  Stage 5 cohort restrictions
                                                                │
                                                                ▼
    ┌──────────────────────────────────────────────────────────────────────┐
    │  derive_cohort(df, audit)        cohort-wise: needs the whole frame   │
    │    ├─ core_above_median          [§13]; median OF THIS FRAME, frozen  │
    │    └─ constant_covariates(…)     detected and logged, never dropped   │
    └──────────────────────────────────────────────────────────────────────┘
                         │
                         ▼
                33 columns  →  Stage 6 propensity  →  Stage 10 bootstrap
                                                      resamples core_above_median,
                                                      never recomputes it

Three of the five are properties of a **record** and can be computed the moment the frame is read.
Two are properties of a **cohort** — a median, and whether a covariate varies — and a cohort does not
exist until Stage 5 has applied the [§3] restrictions. The numbers say why that matters: the core
volume median is 6.0 mL over all 126 records and 5.0 mL over the 93-record primary cohort, so one
function computing both halves at Stage 3's natural position would have split the cohort on the wrong
number, with nothing to notice because both are plausible. **Stage 5 is ``derive_cohort``'s only
caller.**

Both entry points ``df.copy()`` before writing, so a caller's frame is never mutated underneath it,
and both append to the ``Audit`` that ``load()`` created rather than opening a second one. Neither
writes a file.

Every derivation here is a comparison against a column that may be missing, so the whole module rests
on what pandas does with ``<NA>``. The obvious statement of it — "a comparison against a missing value
returns false" — is true of numpy and **false** of the nullable dtypes this pipeline reads into, and
which of the two applies is a Stage 1 decision rather than a Stage 3 one. Verified on pandas 2.3.3::

    Int64    s <= 2                       → <NA>   propagates
    string   s != "witnessed"             → <NA>   propagates
    float64  s > 5.0   (NaN > 5.0)        → False  ABSORBS
    object   s <= 2                       → False  ABSORBS

    Series.mask(cond) where cond carries <NA>:  <NA> is treated as True — the row is masked
                                                and takes the replacement value.

``mrs_90d`` is declared ``Int64`` and propagates; ``core_ml`` declares no dtype and arrives as
``float64`` from the workbook and ``int64`` from the fixture, so it absorbs. **Every mask below stays
regardless**, and none of them is justified by "the comparison would return false" — they are
justified by making the result independent of a dtype no derivation controls. Each site says which of
the two it is today, because a reader told that a redundant call is load-bearing stops believing the
ones that are.

This module is **not** exempt from the Stage 1 §7 raw-name scan and must never become exempt. It
names no raw header, no outcome threshold and no factor level; every one of those is read from
``config.py``.
"""
from __future__ import annotations

from collections.abc import Sequence

import pandas as pd

import config as C
# `_fmt` is private to data.py and is imported anyway, deliberately: it is the pipeline's *one* float
# formatter, and a second one is a second way for two runs to disagree. It is also what renders a
# missing median as the literal `missing` rather than as a pandas repr, which is how an all-missing
# `core_ml` stays visible in the log instead of degrading silently.
from data import Audit, _fmt, absence_by_column


# --- the onset factor [§5, §6] --------------------------------------------------------------------

def _assert_onset_flags(df: pd.DataFrame) -> None:
    """Raise if any record is positive on more than one onset flag.

    Before any level is assigned, because a record positive on both belongs to none of [§5]'s three
    levels and every way of giving it one is a clinical statement nobody has made. A resolution —
    prefer wake-up, prefer unwitnessed, call it witnessed — would apply silently to every future
    workbook, and `onset_type` is a [§6] covariate that enters the propensity model, the outcome
    regression and the balance table. Stage 2 took the same posture toward the 999 groin time and
    DECISION 2 toward the contradictory death flags: a derivation rule where one is available, a loud
    stop where none is.

    v7 has no such record, so this branch is unreachable today. It is written for the workbook that
    has not arrived yet.

    The flags are read from ONSET_TYPE_FROM_FLAG rather than written, for two reasons. The first is
    the ordinary one: a third onset flag is then covered by this assertion by existing, in the same
    edit that gives it a level. The second is particular to this factor — two of the flag *columns*
    are named identically to the *levels* they produce, so writing them here would put a string equal
    to a declared level into this module, which is exactly what §12.4's scan forbids and what would
    let a renamed level drift away from its second declaration.
    """
    flags = [flag for flag, _level in C.ONSET_TYPE_FROM_FLAG]
    # `.fillna(False)` is belt-and-braces here, and the honest reason is not the obvious one.
    # `(flag == 1)` on an Int64 column yields <NA> where the flag is missing; the accumulator is
    # `Int64` too, so without the fill that <NA> would propagate into the count, `<NA> > 1` would be
    # <NA>, and summing a BooleanDtype skips <NA> — the count below would be the same number. It
    # stays for two things it does buy: it makes each mask total, so the expression means the same
    # if a flag ever loses its Int64 declaration and starts absorbing; and it keeps `conflicting`
    # free of <NA>, which `.loc` would otherwise treat as False and silently pass over — verified on
    # pandas 2.3.3, and the opposite of `.mask`, which treats <NA> as True and applies the
    # replacement. `conflicting` carries no <NA> by construction, so the fill is belt-and-braces
    # here; Stage 4 §3.1 is where that asymmetry is load-bearing rather than tidy. A record with one
    # missing flag is not evidence that two are positive — that record is handled in _onset_type,
    # which gives it <NA>.
    #
    # The accumulator must stay `Int64` and not `int64`: `astype("int64")` on a BooleanDtype
    # carrying <NA> raises, which would turn the fill from a guard into a load-bearing call and make
    # the paragraph above false.
    positive = pd.Series(0, index=df.index, dtype="Int64")
    for flag in flags:
        positive = positive + (df[flag] == 1).fillna(False).astype("Int64")

    conflicting = positive > 1
    if int(conflicting.sum()):
        raise C.SchemaError(
            f"{' and '.join(flags)} are positive together on {int(conflicting.sum())} record(s): "
            f"{', '.join(sorted(df.loc[conflicting, 'case_id']))}. [§5]'s onset factor has "
            f"{len(C.FACTOR_LEVELS['onset_type'])} levels and a record positive on more than one "
            "flag belongs to none of them. Correct the workbook or amend [§5]; do not choose a "
            "level here.")


def _onset_type(df: pd.DataFrame) -> pd.Series:
    """The three-level [§5] onset factor, as a `string` — never a Categorical.

    ::

        wake_up  unwitnessed        onset_type
        ───────  ───────────        ──────────────────────────────────────────
           0          0             REFERENCE_LEVELS["onset_type"]  (the baseline)
           1          0             the level ONSET_TYPE_FROM_FLAG gives wake_up
           0          1             the level ONSET_TYPE_FROM_FLAG gives unwitnessed
           1          1             impossible — _assert_onset_flags has already raised
          <NA>        *             <NA>, and NEVER the baseline
           *         <NA>           <NA>, and NEVER the baseline

    The last two rows are the ones an implementer "tidying the masks" deletes, and they are the
    dangerous ones: `onset_type` is a [§6] **covariate**, so a record silently labelled with the
    baseline enters the propensity model, the outcome regression and the balance table carrying a
    fabricated value. `<NA>` drops it by complete-case [§11] instead, which is what [§11] prescribes
    and what the reported denominator will then say. Nothing upstream forbids a missing flag —
    Stage 2's A5 allows missing throughout BINARY_COLUMNS — and v7 has none, so **no test running
    against the real workbook can reach this branch**. `pilots/pilot_data.py` gets it wrong, with an
    `np.select(..., "witnessed")` whose default absorbs a missing flag into the baseline.

    No level name is written here. FACTOR_LEVELS declares the level set and REFERENCE_LEVELS declares
    the baseline; a literal in this module would give each level a second declaration, and the two
    would be free to drift on the one edit — a renamed level — that nothing else would catch.

    A `string` and not a Categorical, for Stage 2 §6's reason about `center`: Stage 5 restricts the
    cohort, and a categorical that keeps a dead level makes every subsequent `groupby(observed=False)`
    resurrect it as an all-missing row — in the balance table, the [§9] within-centre overlap table
    and every subgroup table. Stage 6 builds the categorical from FACTOR_LEVELS at the point of use.
    """
    onset = pd.Series(C.REFERENCE_LEVELS["onset_type"], index=df.index, dtype="string")
    unknown = pd.Series(False, index=df.index)
    for flag, level in C.ONSET_TYPE_FROM_FLAG:
        # Two independent guards, and each one alone is wrong differently. On a record whose flag is
        # <NA>:  both → <NA> (correct);  fillna but no unknown mask → the baseline;  unknown mask but
        # no fillna → <NA>, correct only by accident, because the <NA> condition masks the row to the
        # LEVEL and the second mask then overwrites it;  neither → a fabricated *non-baseline* level,
        # which reads more plausibly in a table than the fabricated baseline does.
        onset = onset.mask((df[flag] == 1).fillna(False), level)
        # Accumulated over the same tuple, so a third onset flag contributes both its level and its
        # missingness in one edit. Written out as `df["wake_up"].isna() | df["unwitnessed"].isna()`
        # it would contribute only the level.
        unknown = unknown | df[flag].isna()
    # The loop's order is not a precedence rule: both flags positive is impossible by
    # _assert_onset_flags. If that assertion were ever deleted the last entry would win silently,
    # which is precisely why it is not optional.
    return onset.mask(unknown)


def _unknown_onset(onset: pd.Series) -> pd.Series:
    """[§13] subgroup: onset was not witnessed. The complement of the baseline over the three levels.

    Built from `onset_type` and **never from the two flags again**. Two derivations of one concept
    drift, and this pair would drift silently: both are correct on v7, and they differ only where a
    flag is missing — where the flag route yields 0 and this route yields `<NA>`.

    `<NA>` propagates and is never 1. "We do not know the onset type" is not "the onset was unknown",
    and the distinction would be easy to lose in the count: 93 of 126 records are already unwitnessed
    or wake-up, so one more would not look wrong.
    """
    return (onset != C.REFERENCE_LEVELS["onset_type"]).astype("Int64").mask(onset.isna())


# --- the outcome dichotomies [§5] -----------------------------------------------------------------

def _dichotomies(df: pd.DataFrame) -> pd.DataFrame:
    """Every [§5] outcome with an ordinal source, built from OUTCOMES and OPS.

    Four columns result, all from `mrs_90d`. `tici_2b_3`, `sich` and `ph2` are read directly from the
    workbook and are already columns of the frame; the `continue` is what leaves them alone, and
    `Outcome.source is None` is the only test of which is which.

    **No threshold, operator or outcome name is written here.** A threshold typed twice can have the
    manuscript print `mRS 0-1` while the analysis computes `mRS 0-2`, with no test failing. The
    iteration order is OUTCOMES' insertion order, which is [§5]'s table order — any order computes
    the same four columns, but a *declared* order is what makes the audit table reproducible.

    `Int64` and not `boolean`: every consumer treats these as 0/1 numerics — weighted means, risk
    differences, the Firth outcome model, the augmented estimator's `Y`. Both keep `<NA>`; only
    `Int64` keeps the arithmetic. They are deliberately **not** added to BINARY_COLUMNS, which is
    Stage 2's domain check over the frame as read: adding them would have Stage 2 assert columns that
    do not exist when it runs.

    **This writes into the frame it is given, and that is deliberate.** It is the one function in this
    module that does; `derive` has already copied, so the frame it hands over is its own and an
    in-place write is free. The consequence is a contract on the caller rather than on the function,
    so it is written down instead of inferred: *a test calling `_dichotomies` directly passes a frame
    it owns and must not assert anything about that frame afterwards.*
    """
    for key, o in C.OUTCOMES.items():
        if o.source is None:
            continue                                    # read directly; nothing to derive
        src = df[o.source]
        # `.mask(src.isna())` is the whole of invariant 6 — but not for the reason it looks like.
        # `mrs_90d` is declared Int64, so the comparison already propagates <NA> and this call is a
        # no-op on every frame this repository builds today. What it buys is that invariant 6 holds
        # for a dtype no derivation controls: on a float64 or object source the comparison ABSORBS,
        # and the records with no 90-day mRS become non-events in all four columns at once — not
        # good-outcome, not poor-outcome, and alive. That frame is not hypothetical: it is what any
        # read whose `dtype=` argument was dropped delivers. The failure is invisible — the column is
        # full, the counts are plausible, every downstream table reconciles — and [§11] names both
        # halves of it: the denominator moves from 124 to 126, and a patient whose outcome is unknown
        # is counted as having had a specific one.
        df[key] = C.OPS[o.op](src, o.threshold).astype("Int64").mask(src.isna())
    return df


# --- the median subgroup [§13] --------------------------------------------------------------------

def _core_above_median(df: pd.DataFrame) -> tuple[pd.Series, float]:
    """[§13] subgroup: core volume above the median **of the frame this is given**.

    ::

        median = core_ml.median()      skipna, of THIS frame — 6.0 mL over all 126 records,
                                       5.0 mL over the 93-record primary cohort
              core > median   →  1
              core == median  →  0     ties sit BELOW. Not a formality: 3 cohort records sit
                                       exactly on 5.0 mL and 39 have a core of exactly 0
              core is missing →  <NA>
                    │
                    ▼
        computed ONCE, on the analysis cohort, and FROZEN. Stage 10 resamples this column
        across its 2000 replicates; it does not recompute it. A replicate that recomputed
        would get its own cut-point and therefore its own subgroup, and no test would fail.

    The mask is **the one mask in Stage 3 that is load-bearing today**. `core_ml` declares no dtype in
    the contract, so it arrives as `float64` from the workbook and `NaN > 5.0` is `False`, not `<NA>`
    — the comparison absorbs. Without it, the one cohort record with no core volume is silently
    assigned to the below-median group, in whichever arm it happens to sit, with no `<NA>` anywhere to
    notice. Contrast `_dichotomies`, where the identical-looking call on an `Int64` source is
    redundant. Complete-case [§11] then drops the record from the subgroup analysis, which is correct.

    An all-missing `core_ml` degrades quietly rather than wrongly: `median()` returns `NaN`, the
    comparison is `False` everywhere, and the mask sets every row to `<NA>` — so the subgroup is empty
    rather than wrong. It is not silent either, because `data._fmt` renders the median as the literal
    `missing` in the audit entry. No guard is added for it; a cohort all-missing on a [§6] covariate
    cannot get past Stage 2's A6 and the balance table without being noticed first.
    """
    core = df["core_ml"]
    median = float(core.median())                       # skipna, the pandas default
    return (core > median).astype("Int64").mask(core.isna()), median


# --- zero-variance covariates ---------------------------------------------------------------------

def constant_covariates(
        df: pd.DataFrame, covariates: Sequence[str] = C.PS_COVARIATES_FULL) -> tuple[str, ...]:
    """Names in `covariates` with at most one distinct non-missing value in `df`, in that order.

    `nunique(dropna=True) <= 1` catches three cases in one expression: a constant column, an
    all-missing column (0 distinct), and a factor left with a single level after a restriction. It
    works unchanged on `onset_type` and `center`, which are strings, and on the numeric covariates.

    Ranges over the caller's sequence, so the order of the result is the caller's and the log is
    deterministic — never over a set.

    **The default argument carries a precondition.** `PS_COVARIATES_FULL` contains `onset_type`, which
    the contract cannot supply, so `constant_covariates(df)` raises `KeyError` on any frame `derive`
    has not run on. That is the correct failure: a caller asking which covariates are constant before
    the covariates exist has a sequencing bug, and silently skipping the missing name would answer a
    question about a different covariate set. `derive_cohort` is downstream of `derive` by
    construction and Stage 6 passes its own design-frame column list, so neither reaches it.

    Returns names and changes nothing — not the frame, and not `covariates`. Detection is Stage 3's;
    the *dropping* is Stage 6's, at the point of use, which is also the only place that can see a
    bootstrap replicate having emptied a factor level. Two stages cannot both own the deletion, and
    this is the division that leaves invariant 3 intact: `OUTCOME_COVARIATES` **is** `PS_COVARIATES`
    by object identity, so a Stage 3 that removed a name from the lists would break the identity
    invariant 3 is asserted on, and one that removed a column from the frame would leave those lists
    naming a column that no longer exists.

    A constant covariate is legitimate — [§6] anticipates exactly one, `first_image_mri`, which the
    contract has already dropped by kind so it cannot reach here at all. What is not legitimate is a
    *silent* drop, which is why this is logged rather than raised.
    """
    return tuple(c for c in covariates if df[c].nunique(dropna=True) <= 1)


# --- the entry points -------------------------------------------------------------------------

def derive(df: pd.DataFrame, audit: Audit) -> pd.DataFrame:
    """The row-wise derivations: the onset factor, its [§13] subgroup, and the four dichotomies.

    Returns a new frame of 31 columns — the 25 Stage 2 delivered, unchanged, plus 6 — and appends
    three `derivation` entries and one `missingness` entry to `audit`. Writes no file.

    Raises SchemaError if any record is positive on both onset flags.
    """
    df = df.copy()
    _assert_onset_flags(df)

    onset = _onset_type(df)
    df["onset_type"] = onset
    counts = onset.value_counts(dropna=False)
    n_missing = int(onset.isna().sum())
    audit.record(
        "derivation", "onset_type", int(onset.notna().sum()),
        f"three levels from wake_up and unwitnessed [§5]; the baseline is "
        f"REFERENCE_LEVELS['onset_type'] = {C.REFERENCE_LEVELS['onset_type']!r}. {n_missing} "
        "record(s) carry <NA> because a source flag is missing — never the baseline [§4.3].",
        table=(("level", "n"),
               *((level, str(int(counts.get(level, 0))))
                 for level in C.FACTOR_LEVELS["onset_type"]),
               ("<NA>", str(n_missing))))

    unknown = _unknown_onset(onset)
    df["unknown_onset"] = unknown
    audit.record(
        "derivation", "unknown_onset", int((unknown == 1).sum()),
        f"[§13] subgroup: onset_type != {C.REFERENCE_LEVELS['onset_type']!r}. "
        f"{int((unknown == 1).sum())} of {len(df)} record(s); {int(unknown.isna().sum())} carry "
        "<NA>, inherited from onset_type and never read as 0.")

    df = _dichotomies(df)
    audit.record(
        "derivation", "dichotomies", len(C.DERIVED_DICHOTOMIES),
        f"{len(C.DERIVED_DICHOTOMIES)} dichotomies derived from their ordinal source [§5], each "
        "carrying exactly its source's missingness. No threshold is written outside OUTCOMES.",
        # `rule` is computed by Outcome from the same fields the loop applies, so the printed
        # derivation and the applied one cannot disagree.
        table=(("outcome", "rule", "n", "events", "missing"),
               *((key, C.OUTCOMES[key].rule, str(int(df[key].notna().sum())),
                  str(int((df[key] == 1).sum())), str(int(df[key].isna().sum())))
                 for key in C.DERIVED_DICHOTOMIES)))

    # Not DERIVED_NAMES: `core_above_median` does not exist yet, and is accounted for in
    # derive_cohort's own entry. Every absence a derived column carries is inherited, so every row of
    # this table is `missing` or `complete` — `structural` means "the value does not exist for this
    # patient", and it is the one label that tells a reader an absence is not data loss.
    absence_by_column(
        df, audit, sorted(C.ROW_WISE_DERIVED), "absence_by_derived_column",
        "one row per derived column; kind is missing or complete. No derived column is structurally "
        "non-applicable: every absence here is inherited from a source column [§8.3].")
    return df


def derive_cohort(df: pd.DataFrame, audit: Audit) -> pd.DataFrame:
    """The cohort-wise derivations. **Stage 5 is the only caller**, after both [§3] restrictions.

    Returns a new frame of 33 columns and appends two `derivation` entries to `audit`. Raises
    SchemaError if called twice on the same frame — see below.
    """
    # Structural, not documentary. The natural shape of a bootstrap replicate is "resample, then
    # re-run the pipeline", and under that shape a derive_cohort that simply recomputed would give
    # every replicate its own cut-point and therefore its own subgroup: membership would move between
    # replicates, the estimand would differ between them, and no test would fail — the estimates would
    # stay plausible and the interval would simply be of a slightly different quantity than the point
    # estimate. Stage 10 must resample the frame that already carries the column.
    already = [c for c in C.COHORT_DEPENDENT_SUBGROUPS if c in df.columns]
    if already:
        raise C.SchemaError(
            f"{', '.join(already)} already present. The [§13] median subgroup is computed once, on "
            "the analysis cohort, and frozen [Stage 3 §6.4]: a bootstrap replicate resamples this "
            "column, it does not recompute it. Calling derive_cohort twice would give the replicate "
            "its own cut-point and its own subgroup, silently.")

    df = df.copy()
    core = df["core_ml"]
    above, median = _core_above_median(df)
    df["core_above_median"] = above
    audit.record(
        "derivation", "core_above_median", int((above == 1).sum()),
        "[§13] subgroup: core_ml > the median of the frame this was computed on. The median is a "
        "property of the cohort [Stage 1 §6] and is frozen here — Stage 10 resamples this column, it "
        "does not recompute it [§6.4]. Ties sit below.",
        table=(("statistic", "value"),
               ("median (mL)", _fmt(median)),
               ("records", str(len(df))),
               ("above", str(int((above == 1).sum()))),
               ("at the median", str(int((core == median).sum()))),
               ("below", str(int((core < median).sum()))),
               ("missing", str(int(core.isna().sum())))))

    constant = constant_covariates(df)
    audit.record(
        "derivation", "constant_covariates", len(constant),
        f"{len(constant)} of {len(C.PS_COVARIATES_FULL)} covariate(s) in PS_COVARIATES_FULL have at "
        "most one distinct non-missing value in this frame. Nothing is dropped here: Stage 6's "
        "design matrix drops constant columns at the point of use [§7.1].",
        table=((("covariate", "distinct values"),
                *((c, str(int(df[c].nunique(dropna=True)))) for c in constant))
               if constant else None))
    return df
