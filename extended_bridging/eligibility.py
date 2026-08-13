"""Stage 4 — eligibility classification [§3, §11].

Every patient carries one of three eligibility classes, decided from ``ivt_contraindicated`` and the
presence of a reason under DECISION 1, with the one combination that would silently absorb a
contradiction raised on instead of resolved, and with the retained set declared once so that
``== "eligible"`` and ``!= "ineligible"`` cannot mean different things in different stages.

This stage **adds one column**. It drops no row, drops no column, and edits no value that Stage 2 or
Stage 3 delivered. ``classify`` calls ``df.copy()`` before writing, so a caller's frame is never
mutated underneath it.

Section references in brackets are to ``statistical_analysis_plan.md``. DECISION *n* refers to the
decisions recorded in ``../out/stage0_data_inventory.md``. The specification for this module is
``specs/stage4_eligibility_classification.md``; nothing here is invented outside it.

Where the stage sits, and why it is a module of its own rather than part of Stage 5's cohort
construction — the two arrows out are the whole of the answer::

    derive(df, audit)  →  31 columns                              [Stage 3 §11]
                         │
                         ▼
    ┌──────────────────────────────────────────────────────────────────────┐
    │  classify(df, audit)                                                  │
    │    ├─ _assert_classifier_inputs(df)   E1…E5 → SchemaError             │
    │    ├─ _classify(df)                   the four cases, every mask      │
    │    │                                  filled, labels from config      │
    │    └─ _crosstab(df)                   centre x arm x class, every     │
    │                                       cell rendered; the cells must   │
    │                                       reconcile to len(df) or raise   │
    │                                                                       │
    │  retained(df)                         != ineligible, from             │
    │                                       ELIGIBILITY_RETAINED   [§5.1]   │
    └──────────────────────────────────────────────────────────────────────┘
                         │
                         ▼
                32 columns  →  Stage 5 restrictions  →  derive_cohort  →  33 columns
                         │                                                    │
                         │  the UNRESTRICTED frame keeps the column,           ▼
                         └─→  Stage 13 [§14b]                           Stages 6-11
                              the only analysis whose population
                              includes ineligible patients

[§14b] is "the only §14 analysis whose population includes contraindicated patients", and its active
regime assigns treatment *as a function of eligibility*. So the classification has to exist on the
unrestricted 126-record frame and survive independently of cohort construction: a classifier living
inside a Stage 5 cohort function would compute the one thing [§14b] needs and then discard the rows
that need it.

**Stage 4 has no dependency on Stage 3.** It reads ``ivt_contraindicated``,
``contraindication_reason``, ``ivt``, ``center`` and ``case_id``, all of which Stage 2 delivers, and
reads nothing ``derive`` produces. It is placed after ``derive`` because Stage 3's diagram places it
there and one pipeline order is better than two, not because it needs anything from there.

This module is one long chain of masks over nullable columns, and the pandas facts it turns on point
in **opposite** directions. Verified on pandas 2.3.3, the version ``uv.lock`` pins::

    (flag == 1)  on Int64 carrying <NA>          → <NA>                    propagates
    (flag == 0) & reason.notna()  with <NA>      → <NA> where reason is present
                                                   False where it is absent

    Series.mask(cond)      with <NA> in cond     → <NA> is treated as TRUE
                                                   the row IS masked, and takes the replacement
    df.loc[cond] / s[cond] with <NA> in cond     → <NA> is treated as FALSE
                                                   the row is NOT selected

    Series.isin(values)    with <NA> in the      → False, NEVER <NA>
                           series                  so ~isin SELECTS the missing row, and no fill
                                                   is needed — or wanted — on an isin mask

The third fact is why ``~…isin(…)`` is the right shape for E4, E5 and ``retained``, and why none of
them carries a ``.fillna(False)``. The direction of the damage is worth being blunt about, because it
is the opposite of the harmless one: under the mask chain of ``_classify`` with the fills removed, a
record whose flag is missing is classified ``ineligible`` — and an ineligible record is one Stage 5
deletes. The failure is not a fabricated label in a table. It is a patient silently removed from the
analysis, with a full column, plausible counts and every downstream denominator reconciling.

This module is **not** exempt from the Stage 1 §7 raw-name scan and must never become exempt. It names
no raw header — ``IVT_contraindicated_binary`` and ``Contraindications_to_IVT`` are read by their
analysis names — and it writes no class label: all three come from ``config.py``.
"""
from __future__ import annotations

import pandas as pd

import config as C
from data import Audit


def _assert_classifier_inputs(df: pd.DataFrame) -> None:
    """The five preconditions of the classifier, reported together rather than one at a time.

    Three decide whether the mask chain in `_classify` can be trusted at all (E1, E2, E5), one is the
    contradiction the *rule* would otherwise absorb (E3), and one is the contradiction the *table*
    would otherwise absorb (E4). Collected rather than raised singly, following Stage 2's
    `_assert_schema`: five independent properties of four columns are worth reporting together, so a
    corrected workbook is diagnosed in one run instead of five.

    Four of the five have a Stage 2 counterpart — E1 is A4b, E2 is A5, E5 is A4 and A5 together, E4 is
    close to A3 — and they are re-checked here deliberately. `classify` is called on frames Stage 2
    never saw: Stage 12 and Stage 13 take the unrestricted frame, every acceptance test builds its own,
    and a notebook gets no A4b at all. An assertion that only runs when the caller happened to come
    through `load()` is a guarantee about a code path, not about a function.

    Every branch is unreachable on v7 — the flag is 0/1 on all 126 records, never missing and never 1
    on a treated patient; the exposure likewise; every centre is declared — and that is recorded rather
    than treated as a reason to skip one. Like Stage 3's `_assert_onset_flags`, these are written for
    the workbook that has not arrived yet.
    """
    flag, treated = df["ivt_contraindicated"], df[C.TREATMENT] == 1
    bad: list[str] = []

    missing = df.loc[flag.isna(), "case_id"]
    if len(missing):
        bad.append(
            f"E1  ivt_contraindicated: {len(missing)} record(s) carry no classifier value: "
            f"{', '.join(sorted(missing))}. [§3] forbids reading an absent contraindication flag "
            "as 'no contraindication', and this stage will not invent a class for it.")

    outside = df.loc[flag.notna() & ~flag.isin((0, 1)), "case_id"]
    if len(outside):
        bad.append(
            f"E2  ivt_contraindicated: {len(outside)} record(s) outside {{0, 1}}: "
            f"{', '.join(sorted(outside))}. The classifier tests equality against both values; a "
            "third value would fall through to indeterminate and read as an undocumented reason.")

    conflicting = df.loc[(treated & (flag == 1)).fillna(False), "case_id"]
    if len(conflicting):
        bad.append(
            f"E3  {C.TREATMENT} and ivt_contraindicated: {len(conflicting)} treated record(s) are "
            f"flagged contraindicated: {', '.join(sorted(conflicting))}. One of the two is wrong "
            "and the revealed-fact rule would absorb the contradiction without a trace. Resolve it "
            "with the data owner; do not choose a class here.")

    # E4 and E5 carry NO `.fillna(False)`, four lines from three masks in `_classify` that must all
    # carry one. `isin` answers a membership question and a missing value is not a member, so it
    # returns False rather than propagating: `~isin` *catches* the missing row, and the mask carries no
    # <NA> at all — numpy `bool` on the `string` centre, pandas `boolean` on the Int64 exposure, and
    # neither nullable in practice. Verified on pandas 2.3.3; §12.3 pins both halves.
    #
    # So a fill here would be a no-op rather than a bug — but it is still wrong to write one, and the
    # reason is not the one it is tempting to give. It would assert that this mask can carry <NA> when
    # it cannot, which is what makes a reader stop believing the three fills in `_classify` that are
    # load-bearing. And it would become a real bug the moment either check is rewritten as a comparison
    # chain, where <NA> does propagate and a fill flips "never recorded" from caught to ignored.
    off_centre = df.loc[~df["center"].isin(C.CENTER_ORDER), "case_id"]
    if len(off_centre):
        bad.append(
            f"E4  center: {len(off_centre)} record(s) carry a centre outside CENTER_ORDER: "
            f"{', '.join(sorted(off_centre))}. _crosstab renders one row per declared centre, so "
            "such a record would be classified and then vanish from the cross-tabulation, leaving a "
            "table whose empty cells no longer mean what [§3] restriction 1 reads them as meaning.")

    off_arm = df.loc[~df[C.TREATMENT].isin((0, 1)), "case_id"]
    if len(off_arm):
        bad.append(
            f"E5  {C.TREATMENT}: {len(off_arm)} record(s) whose exposure is missing or outside "
            f"{{0, 1}}: {', '.join(sorted(off_arm))}. The classifier reads the exposure twice — E3 "
            "and the revealed-fact mask — and both read an absent value as 'not treated', so a "
            "treated patient with no recorded exposure would be classified from the flag alone.")

    if bad:
        raise C.SchemaError(
            "\n  " + "\n  ".join(bad)
            + f"\n  {len(bad)} eligibility assertion(s) failed against {len(df)} records. [§3] "
            "classifies from time-zero information applied symmetrically across arms; none of "
            "these is a reason to relax that.")


# DECISION 1's four cases, and the <NA> row that the fills below are about:
#
#   treated                                        →  eligible       by revealed fact
#   ivt_contraindicated == 1                       →  ineligible     [§3] restriction 2
#   ivt_contraindicated == 0, a reason recorded     →  eligible       documented: no absolute
#                                                                    contraindication
#   ivt_contraindicated == 0, none recorded         →  indeterminate  [§3]: a blank is NEVER read
#                                                                    as "no contraindication"
#   ivt_contraindicated missing                     →  reaches no case, and WITHOUT the fills takes
#                                                     the last one it appears in — `ineligible`,
#                                                     which Stage 5 deletes. That is the line an
#                                                     implementer "tidying the masks" deletes.
#
# The free text is never read. `contraindication_reason` contributes exactly one bit — whether a reason
# was recorded at all — so `.notna()` is the whole of it: no strip, no case fold, no membership test.
# Stage 2's A9 already raises on a whitespace-only cell, and it is the one Stage 2 assertion this stage
# depends on: if A9 is ever relaxed, this classifier must be revisited in the same commit.
def _classify(df: pd.DataFrame) -> pd.Series:
    flag, reason = df["ivt_contraindicated"], df["contraindication_reason"]
    out = pd.Series(C.INDETERMINATE, index=df.index, dtype="string")
    out = out.mask(((flag == 0) & reason.notna()).fillna(False), C.ELIGIBLE)
    out = out.mask((flag == 1).fillna(False), C.INELIGIBLE)
    out = out.mask((df[C.TREATMENT] == 1).fillna(False), C.ELIGIBLE)
    return out


def _crosstab(df: pd.DataFrame) -> tuple[tuple[str, ...], ...]:
    """Eligibility by centre and arm, every declared cell rendered whether or not the data fills it.

    Never `pd.crosstab`, which drops absent combinations. The empty cells are the information: on v7
    the `USZ / bridging` row is all zeros, and that row **is** the [§3] restriction 1 finding — the
    never-IVT centre, visible one stage before Stage 5 acts on it.

    Takes the frame, not the frame and the Series: `classify` assigns the column and then tabulates
    from the column it is about to return, so the table is provably about the shipped column rather
    than about a Series computed alongside it.

    Raises SchemaError if the cells do not reconcile to `len(df)`. That is a second guard independent
    of E4 rather than belt-and-braces — E4 catches an *undeclared centre value*, this catches any
    reason a patient falls in no cell: a TREATMENT_LABELS that stopped covering the arm coding, a
    CENTER_ORDER shortened by an edit, a future third dimension. The `all / both` row cannot
    substitute, because it is computed from the frame rather than from the cells, so a vanished patient
    leaves it correct.
    """
    eligibility = df[C.ELIGIBILITY]
    header = ("centre", "arm", *C.ELIGIBILITY_ORDER, "n")
    rows: list[tuple[str, ...]] = []
    rendered = 0
    for centre in C.CENTER_ORDER:
        # sorted, so the arm order is 0 then 1 — control before treated — in every interpreter. A dict
        # iteration would be equally deterministic today and would stop being a *declared* order.
        for code in sorted(C.TREATMENT_LABELS):
            cell = ((df["center"] == centre) & (df[C.TREATMENT] == code)).fillna(False)
            rendered += int(cell.sum())
            rows.append((
                centre, C.TREATMENT_LABELS[code],
                *(str(int((cell & (eligibility == cls)).sum())) for cls in C.ELIGIBILITY_ORDER),
                str(int(cell.sum()))))
    if rendered != len(df):
        raise C.SchemaError(
            f"\n  the cross-tabulation renders {rendered} of {len(df)} record(s). Every patient "
            "must fall in exactly one CENTER_ORDER x TREATMENT_LABELS cell. The absent ones would "
            "be classified and then missing from the table, and the `all / both` row would still "
            "reconcile against the frame, so nothing downstream would disagree.")
    rows.append((
        "all", "both",
        *(str(int((eligibility == cls).sum())) for cls in C.ELIGIBILITY_ORDER),
        str(len(df))))
    return (header, *rows)


def classify(df: pd.DataFrame, audit: Audit) -> pd.DataFrame:
    """[§3] eligibility under DECISION 1. Appends one column and one audit entry; drops nothing.

    Read `ivt_contraindicated`, `contraindication_reason`, `ivt`, `center` and `case_id` — all of
    them Stage 2's — and nothing `derive` produces. That independence is what Stage 12 and Stage 13
    rely on when they classify the unrestricted frame.

    The assertion runs **before** the copy, so a frame that fails is never copied and E1-E5 are
    reported about the caller's own frame. `_assert_classifier_inputs` only reads, so ordering it
    first is safe as well as cheaper.

    The audit entry's kind is `derivation` — an entry describing a column built from other columns,
    which is what this is. It names no patient: the one place this stage names patients is the
    assertion, which raises, so such a frame never reaches a log.
    """
    _assert_classifier_inputs(df)
    df = df.copy()
    df[C.ELIGIBILITY] = _classify(df)

    # From the produced column, via one pass — never from a second traversal of the flag and the
    # reason, which would give the log and the frame two sources for one number. `.get(cls, 0)` and
    # not `counts[cls]`: a frame in which no patient is ineligible — the fixture is one — has no
    # `ineligible` entry, and indexing it would raise on the cleanest possible input.
    counts = df[C.ELIGIBILITY].value_counts()
    n_e, n_i, n_x = (int(counts.get(cls, 0))
                     for cls in (C.ELIGIBLE, C.INDETERMINATE, C.INELIGIBLE))
    audit.record(
        "derivation", C.ELIGIBILITY, len(df),
        f"[§3] eligibility under DECISION 1: treated → {C.ELIGIBLE}; ivt_contraindicated = 1 → "
        f"{C.INELIGIBLE}; flag = 0 with a reason recorded → {C.ELIGIBLE}; flag = 0 with none "
        f"recorded → {C.INDETERMINATE}. The free text is never read — the reason column contributes "
        f"one bit, whether a reason exists. {n_e} eligible, {n_i} indeterminate, {n_x} ineligible "
        f'of {len(df)} record(s). A blank is never read as "no contraindication" [§3]; the '
        "indeterminate group is retained and Stage 5 drops only the ineligible.",
        table=_crosstab(df))
    return df


def retained(df: pd.DataFrame) -> pd.Series:
    """[§3]'s retained set: every class but `ineligible`. Boolean, and total — never <NA>.

    The one place the difference between `== eligible` and `!= ineligible` is decided. It is 43 of
    126 records — 80% of the control arm in the primary cohort — so two stages resolving it
    differently would not look like a bug in either of them. Stage 5 restricts on this, [§14a]
    draws its population with it, and [§14b] uses the three-level column instead.

    `isin` over the declared tuple, never `!= C.INELIGIBLE`. The two agree today and would diverge on
    the one edit — a fourth class — that nothing else would catch, and the tuple is where a fourth
    class would be declared.

    Raises KeyError on a frame `classify` has not run on. That is the correct failure: a caller
    asking who is retained before anyone has been classified has a sequencing bug, and a predicate
    that answered anyway would answer about a different population.
    """
    return df[C.ELIGIBILITY].isin(C.ELIGIBILITY_RETAINED)
