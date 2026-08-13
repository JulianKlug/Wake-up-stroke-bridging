"""Stage 5 — cohort construction [§2, §3].

The [§3] primary cohort: two restrictions applied in a declared order, each naming the patients it
removed, followed by Stage 3's cohort-wise derivations so that the [§13] median subgroup is frozen on
the cohort rather than on the frame it came from.

This stage **removes rows and adds one column**. It edits no value that Stages 2, 3 or 4 delivered,
and ``build`` calls ``df.copy()`` before writing, so a caller's frame is never mutated underneath it.
The index is preserved and never reset: a cohort record keeps the label it carried on the unrestricted
frame, which is what makes it traceable back to the population Stages 12 and 13 analyse.

Section references in brackets are to ``statistical_analysis_plan.md``. DECISION *n* refers to the
decisions recorded in ``../out/stage0_data_inventory.md``. The specification for this module is
``specs/stage5_cohort_construction.md``; nothing here is invented outside it.

Where the stage sits — and the two arrows out are the whole of why it is a module of its own::

    classify(df, audit)  →  32 columns, 126 rows                     [Stage 4 §11]
                       │
                       ├──────────────────────────────→  Stage 12 [§14a], Stage 13 [§14b]
                       │                                 the UNRESTRICTED frame, 126 rows
                       ▼
    ┌────────────────────────────────────────────────────────────────────────┐
    │  STAGE 5 — cohort.py           reads config.py, data.Audit, eligibility │
    │                                                                        │
    │   build(df, audit)                                                     │
    │     ├─ _assert_cohort_inputs(df)    C1…C4 → SchemaError                │
    │     ├─ restriction 1, inline        [§3] centres          126 → 104    │
    │     ├─ restriction 2, inline        [§3] eligibility      104 →  93    │
    │     ├─ _flow_table(…)               one row per step, + the            │
    │     │                               no-reason-on-file line             │
    │     ├─ derive.derive_cohort(…)      Stage 3's second entry point,      │
    │     │                               HERE and nowhere else  32 → 33     │
    │     ├─ absence_by_column(…)         the cohort's own denominators [§11]│
    │     └─ _assert_cohort(out)          P1…P4 → SchemaError                │
    │                                                                        │
    │   treating_centres(df)              the [§3] restriction-1 predicate,  │
    │                                     computed, never declared           │
    └────────────────────────────────────────────────────────────────────────┘
                       │
                       ▼
              33 columns, 93 rows  →  Stages 6-11

The second arrow is why this module removes rows from a *copy* rather than from the pipeline: [§14a]
and [§14b] take the unrestricted classified frame, and a restriction living inside a Stage 6 propensity
function would be recomputed by every caller that wanted a cohort and by none that wanted the frame.

``build`` appends four entries to the ``Audit`` that ``load()`` created — three of the ``cohort`` kind
and one ``missingness`` — plus the two ``derivation`` entries ``derive_cohort`` records. It writes no
file.

**This module imports ``derive`` and ``eligibility``, which is a first for this pipeline.** It is
deliberate: ``build`` calls ``derive.derive_cohort`` so that nothing else has to, which is what makes
the [§13] median a property of the cohort rather than of whichever frame a driver happened to pass;
and it calls ``eligibility.retained`` so that restriction 2 cannot be written as a comparison of this
stage's own.

The pandas facts this stage turns on, verified on pandas 2.3.3 as pinned by ``uv.lock``::

    df[mask]               with <NA> in mask       → the row is NOT selected
                                                     (Stage 4 §3.1 records the same fact for .loc)

    Series.isin(values)    with <NA> in the series → False, NEVER <NA>
                                                     so ~isin SELECTS the missing row

    df[mask] preserves the index                   → cohort labels are the frame's labels
    df[mask] on a `string` column leaves no        → no dead level; `center` is not a Categorical,
        dead level                                   per Stage 3 §4.4

The first fact is why the preconditions come before the masks and not after. A record whose exposure or
centre is missing is *silently excluded* by ``df[mask]``: it is not selected by restriction 1's
keep-mask, so it leaves the cohort with no entry anywhere saying it did, the flow table still
reconciles, and the identifier appears among restriction 1's removals as though the patient had been at
a never-treating centre. C3 and C4 make that unreachable. The fourth fact is why nothing here calls
``remove_unused_categories``: ``center`` is a ``string`` by Stage 2 §6's decision, taken for exactly
this reason.

This module is **not** exempt from the Stage 1 §7 raw-name scan and must never become exempt. It names
no raw header and writes no eligibility class label: the restriction reads ``eligibility.retained``,
which reads ``ELIGIBILITY_RETAINED``.
"""
from __future__ import annotations

from typing import Final

import pandas as pd

import config as C
import derive
import eligibility
from data import Audit, absence_by_column


# --- restriction 1's predicate --------------------------------------------------------------------

def treating_centres(df: pd.DataFrame) -> tuple[str, ...]:
    """The centres in `df` that contributed at least one bridging patient, in CENTER_ORDER.

    [§3] restriction 1 keeps these and drops every other declared centre. Computed from the frame
    rather than declared as a constant: a centre's treatment availability is a property of the data,
    and the one edit that would matter — a workbook in which USZ starts administering IVT, or in
    which Lugano's two bridging patients turn out to be miscoded — must change the cohort rather
    than require someone to remember to change a list. `pilots/pilot_config.py:78` declares
    `NEVER_IVT_CENTERS = ["USZ"]`; Stage 5 §16 is where that is declined.

    Ranges over CENTER_ORDER and never over the frame's observed values, so the result is a declared
    order and the log is deterministic. A declared centre contributing no rows at all is therefore
    not returned — vacuously correct, since it contributed no bridging patient, and visible in
    §7.1's table as an all-zero row rather than as an absence.

    Public, and not for symmetry: **Stage 12's support check needs its complement.** [§14a] requires
    the proportion of never-IVT-centre patients outside the treated support, and a baseline table
    comparing them with the IVT-treated — both range over the centres this does *not* return. A second
    implementation over there would be a second definition of the restriction-1 population, and the two
    would agree on v7 and be free to drift on the workbook that follows it. The *dropped* set is
    `tuple(c for c in C.CENTER_ORDER if c not in keep)`, computed at the one place that needs it.

    The count is over the bridging arm only, because [§3] restriction 1 names one direction — the
    direction v7 showed. A centre with no *control* arm is a positivity violation too; P2 raises on it
    rather than silently extending [§3]. The asymmetry is [§3]'s and this stage does not resolve it.
    """
    return tuple(c for c in C.CENTER_ORDER
                 if int(((df["center"] == c) & (df[C.TREATMENT] == 1)).sum()) > 0)


# --- the preconditions [§3] -----------------------------------------------------------------------
#
# Before any row is removed, and before the frame is copied, following Stage 2's `_assert_schema` and
# Stage 4's `_assert_classifier_inputs`. C4 is emitted before C3 although it is numbered after: the
# numbering follows Stage 4's E-numbering of the same two properties, so a reader moving between the
# files is not made to re-learn which is which, while the emission order follows the columns' order in
# the frame. The order is inert, because the messages are collected.
#
# Every branch is unreachable on v7 — every centre is declared and `ivt` is 0/1 and never missing —
# and that is recorded rather than treated as a reason to skip one. Like Stage 3's
# `_assert_onset_flags` and Stage 4's E1-E5, they are written for the workbook that has not arrived
# yet.

def _assert_cohort_inputs(df: pd.DataFrame) -> None:
    """The four preconditions, collected: one SchemaError reports every one that failed."""
    bad: list[str] = []

    if C.ELIGIBILITY not in df.columns:
        bad.append(
            f"C1  {C.ELIGIBILITY}: absent. Stage 4's classify() has not run on this frame. [§3]'s "
            "second restriction is a function of eligibility and this stage will not classify: "
            "eligibility.classify() is one call and doing it here would give the column two "
            "producers and the log two entries.")

    already = [c for c in C.COHORT_DEPENDENT_SUBGROUPS if c in df.columns]
    if already:
        bad.append(
            f"C2  {', '.join(already)}: already present. The [§13] median subgroup is frozen on the "
            "COHORT [Stage 3 §6.4], so a frame that already carries it was restricted once already. "
            "derive_cohort raises on this too — but only after both restrictions have been applied "
            "and logged, leaving an audit log describing a cohort no caller received.")

    off_centre = df.loc[~df["center"].isin(C.CENTER_ORDER), "case_id"]
    if len(off_centre):
        bad.append(
            f"C4  center: {len(off_centre)} record(s) carry a centre outside CENTER_ORDER: "
            f"{', '.join(sorted(off_centre))}. treating_centres ranges over CENTER_ORDER, so such a "
            "record is kept by no centre and dropped by restriction 1 without ever having been at a "
            "never-treating centre. The flow table would reconcile and the identifier would appear "
            "in restriction 1's removals.")

    off_arm = df.loc[~df[C.TREATMENT].isin((0, 1)), "case_id"]
    if len(off_arm):
        bad.append(
            f"C3  {C.TREATMENT}: {len(off_arm)} record(s) whose exposure is missing or outside "
            f"{{0, 1}}: {', '.join(sorted(off_arm))}. restriction 1's predicate counts bridging "
            "patients per centre, so a treated patient with no recorded exposure could take their "
            "centre's count to zero and remove the whole centre.")

    if bad:
        raise C.SchemaError(
            "\n  " + "\n  ".join(bad)
            + f"\n  {len(bad)} cohort assertion(s) failed against {len(df)} records. [§3]'s two "
            "restrictions are applied by design, before any weighting; none of these is a reason "
            "to relax that.")


# --- the one way a restriction reaches the log ----------------------------------------------------

def _record_removal(audit: Audit, step: str, removed: pd.DataFrame, detail: str,
                    table: tuple[tuple[str, ...], ...]) -> None:
    """The one way a [§3] restriction reaches the log. Names every record it removed, or raises.

    Stronger than data.py's kind-keyed rule, which asks only that *some* case be named: this asserts
    the entry names EXACTLY as many patients as it says it removed. A restriction that reported 11
    and named 3 would satisfy _MUST_NAME_CASES and would be a cohort nobody could reconstruct
    [Stage 5 §6.3].

    The identifiers are collected through a `set`, and a missing one is dropped rather than sorted
    against a string. That is what makes the comparison able to fire at all: a duplicated identifier
    collapses to one name and a missing one contributes none, so either leaves the entry naming fewer
    patients than it removed. Both are forbidden by Stage 2's A2, so this is a belt over Stage 2's
    braces — written here because this is the stage whose output is unreconstructable without it.
    """
    named = removed.loc[removed["case_id"].notna(), "case_id"]
    case_ids = tuple(sorted(set(named)))
    if len(case_ids) != len(removed):
        raise C.SchemaError(
            f"{step} removed {len(removed)} record(s) and names {len(case_ids)}. A [§3] restriction "
            "that cannot name what it removed makes the cohort unreconstructable from the log.")
    audit.record("cohort", step, len(removed), detail, case_ids=case_ids, table=table)


# --- the four audit entries [§7] -------------------------------------------------------------------
#
# `contraindication_reason` is read by exactly three of the helpers below, for COUNTING and never for
# classifying, and that distinction is the whole of what survives DECISION 1a. The column may reach a
# table or a `detail` string; it must never reach a mask that decides a class, and §12.6's
# per-function AST scan holds the module to the three names. So the arm-restricted mask is written out
# at each of the three sites rather than factored into a shared helper: a helper *returning* a mask
# over that column is a mask that can escape into a restriction, and the scan — which sees only where
# the column is named — would not notice.
#
# Arm-restricted, and the restriction is not cosmetic. Over all arms the same count is 82 on v7 rather
# than 43, because no treated patient carries a reason either — and a count over all arms would be
# measuring nothing, since a treated patient's eligibility is established by the fact that they
# received IVT. The assumption DECISION 1a introduces is about CONTROLS with no documented reason.
#
# The arm columns come from TREATMENT_LABELS in sorted() key order — control before treated, in every
# interpreter, and identical to Stage 4 §7.1's crosstab so the two tables can be read against each
# other. The class columns come from ELIGIBILITY_ORDER. `pd.crosstab` is used nowhere: it drops absent
# combinations, and the empty cells are the information.

_ARM_CODES: Final[tuple[int, ...]] = tuple(sorted(C.TREATMENT_LABELS))
_ARM_LABELS: Final[tuple[str, ...]] = tuple(C.TREATMENT_LABELS[code] for code in _ARM_CODES)


def _arm_cells(df: pd.DataFrame) -> tuple[str, ...]:
    return tuple(str(int((df[C.TREATMENT] == code).sum())) for code in _ARM_CODES)


def _class_cells(df: pd.DataFrame) -> tuple[str, ...]:
    return tuple(str(int((df[C.ELIGIBILITY] == cls).sum())) for cls in C.ELIGIBILITY_ORDER)


def _centres_detail(keep: tuple[str, ...], dropped: pd.DataFrame) -> str:
    """§7.1's detail. Takes no frame: every number in it is a count of centres or of removed records.

    §8 writes this call as `_centres_detail(df, keep, dropped)`; the frame argument would be unused, so
    it is dropped here. Nothing else about §8's twelve lines moves.
    """
    return (
        f"[§3] restriction 1: {len(keep)} of {len(C.CENTER_ORDER)} declared centre(s) contributed at "
        f"least one {C.TREATMENT_LABELS[1]} patient and are retained; "
        f"{len(C.CENTER_ORDER) - len(keep)} contributed none and are dropped with their "
        f"{len(dropped)} patient(s). P(IVT = 1 | centre, X) = 0 there structurally, not by chance — "
        "no weighting recovers a contrast that was never available, and a penalised score would "
        "report shrinkage as though it were treatment availability.")


def _centres_table(df: pd.DataFrame, keep: tuple[str, ...],
                   dropped: pd.DataFrame) -> tuple[tuple[str, ...], ...]:
    """One row per declared centre, rendered whether or not the data fills it, then `all`.

    The empty cells are the information: this is the table in which the [§3] restriction-1 finding is
    the all-zero `bridging` cell. The `status` cell is derived from membership of `keep`, never
    recomputed, so the table and the restriction cannot disagree.

    **The `all` row is not a reconciliation.** It is computed from the frame, so a record at an
    undeclared centre would leave it correct while vanishing from the four centre rows. That failure
    is C4's, which raises before this table is built, and P4's, which is the stronger of the two
    because it sums against the input frame rather than against a row of this table.
    """
    header = ("centre", *_ARM_LABELS, "n", "status")
    rows = tuple(
        (centre, *_arm_cells(df[df["center"] == centre]),
         str(int((df["center"] == centre).sum())),
         "retained" if centre in keep else "dropped")
        for centre in C.CENTER_ORDER)
    return (header, *rows,
            ("all", *_arm_cells(df), str(len(df)), f"{len(dropped)} removed"))


def _eligibility_detail(after_1: pd.DataFrame, ineligible: pd.DataFrame) -> str:
    reason = after_1["contraindication_reason"]
    undocumented = int((reason.isna() & (after_1[C.TREATMENT] == 0)
                        & eligibility.retained(after_1)).sum())
    # In CENTER_ORDER, and a centre contributing no rows does not count: `notna().sum()` is vacuously
    # zero there, and reporting the dropped centre as one that never collected a reason would describe
    # an absence of patients as an absence of documentation.
    centres = [centre for centre in C.CENTER_ORDER
               if int((after_1["center"] == centre).sum())
               and not int(reason[after_1["center"] == centre].notna().sum())]
    return (
        f"[§3] restriction 2: {len(ineligible)} patient(s) carry ivt_contraindicated = 1 and are "
        "removed — they were never candidates for bridging, and retaining them makes \"no IVT\" a "
        "marker of contraindications and their prognosis. Eligibility is classified from the flag "
        f"alone [DECISION 1a]; the {undocumented} retained {C.TREATMENT_LABELS[0]} patient(s) with "
        "no documented contraindication reason are eligible, and the free text is not read. The flag "
        f"was recorded on every record, including at the {len(centres)} centre(s) that never "
        "collected a reason — that it means the same thing at all of them is an assumption, and [§3] "
        "requires it stated as a limitation.")


def _eligibility_table(after_1: pd.DataFrame,
                       keep: tuple[str, ...]) -> tuple[tuple[str, ...], ...]:
    """One row per **retained** centre, over the frame restriction 1 produced, then `all`.

    Ranging over `keep` and not over CENTER_ORDER is deliberate: the never-treating centre is gone by
    the time this restriction runs, and rendering an all-zero row for it would suggest its patients
    were considered for eligibility and found eligible. `keep` is declared, so the order is still not
    the data's.

    The `removed` column is the `ineligible` column by construction. Printed anyway, because it is the
    column a reader checks the entry's `n` against, and because if the two ever differ the restriction
    has stopped being "drop the ineligible". The `all / eligible` cell is the cohort's size — an
    identity rather than a coincidence, and the one place the two restrictions' arithmetic meets before
    the flow table states it.
    """
    header = ("centre", *C.ELIGIBILITY_ORDER, "n", "removed")
    rows = tuple(
        (centre, *_class_cells(after_1[after_1["center"] == centre]),
         str(int((after_1["center"] == centre).sum())),
         str(int(((after_1["center"] == centre)
                  & ~eligibility.retained(after_1)).sum())))
        for centre in keep)
    return (header, *rows,
            ("all", *_class_cells(after_1), str(len(after_1)),
             str(int((~eligibility.retained(after_1)).sum()))))


def _flow_detail(source: pd.DataFrame, cohort: pd.DataFrame) -> str:
    treated, control = (int((cohort[C.TREATMENT] == code).sum()) for code in (1, 0))
    undocumented = int((cohort["contraindication_reason"].isna()
                        & (cohort[C.TREATMENT] == 0)).sum())
    return (
        f"[§3] cohort flow. The primary cohort is {len(cohort)} of {len(source)} record(s) at "
        f"{int(cohort['center'].nunique())} centre(s): {treated} {C.TREATMENT_LABELS[1]} and "
        f"{control} {C.TREATMENT_LABELS[0]}. Restrictions are applied in the declared order and each "
        "record is counted once, under the first restriction that removed it — the retained set does "
        "not depend on that order but the attribution does [Stage 5 §4.5]. "
        f"{undocumented} retained {C.TREATMENT_LABELS[0]} patient(s) have no documented "
        "contraindication reason and are eligible from the flag alone [DECISION 1a].")


def _flow_table(source: pd.DataFrame, after_1: pd.DataFrame,
                after_2: pd.DataFrame) -> tuple[tuple[str, ...], ...]:
    """Four rows, on every frame, each computed from the frame at that step.

    Never by subtracting the row above: subtraction would make the table reconcile with itself by
    construction and stop being evidence.

    Every `centres` cell is `nunique()` over the row's own frame — the centres contributing at least
    one record — never `len(CENTER_ORDER)` and never the retained tuple's length. The three agree on
    v7 by coincidence and disagree on the fixture, so the rule is declared here rather than left to be
    inferred from the table.

    The last row survives DECISION 1a and no longer reports a *class*: it reports the same patients as
    a **provenance** count, the retained controls whose eligibility rests on the flag with no
    documented reason behind it. That is what a reader needs in order to size the assumption DECISION
    1a makes, and it is the group a future [§13] sensitivity analysis would be sized against. Its
    `centres` cell is the finding — the undocumented group is centre-driven, not patient-driven.
    """
    no_reason = after_2[after_2["contraindication_reason"].isna()
                        & (after_2[C.TREATMENT] == 0)]
    steps = (
        ("as classified", source),
        ("after restriction 1", after_1),
        ("after restriction 2", after_2),
        (f"of which {C.TREATMENT_LABELS[0]}, no reason on file", no_reason),
    )
    header = ("step", "centres", "records", *_ARM_LABELS, *C.ELIGIBILITY_ORDER)
    return (header, *(
        (label, str(int(frame["center"].nunique())), str(len(frame)),
         *_arm_cells(frame), *_class_cells(frame))
        for label, frame in steps))


# The detail of the fourth entry is a module constant rather than a literal at the call site, so
# `build` stays twelve readable lines; the other three are computed per call and stay in their helpers.
_ABSENCE_DETAIL: Final[str] = (
    "one row per analysis and derived column, over the PRIMARY COHORT rather than over the frame "
    "as read [§11]. Every estimate carries its own denominator; this is where the denominators "
    "of the [§7] analyses are read off. Nothing here is imputed.")


# --- the postconditions ----------------------------------------------------------------------------

def _assert_cohort(df: pd.DataFrame, source: pd.DataFrame, removed: int) -> None:
    """The four postconditions, collected. `source` is the frame `build` was given, bound before the
    copy; `removed` is restriction 1's count plus restriction 2's. Raises SchemaError or returns."""
    bad: list[str] = []

    remaining = df.loc[df[C.ELIGIBILITY] == C.INELIGIBLE, "case_id"]
    if len(remaining):
        bad.append(
            f"P1  {C.ELIGIBILITY}: {len(remaining)} ineligible patient(s) remain in the cohort: "
            f"{', '.join(sorted(remaining))}. [§3] restriction 2 removes exactly these, and roadmap "
            "invariant 2 forbids one reaching [§7]. This is what fails if the mask is inverted.")

    empty = [(centre, C.TREATMENT_LABELS[code])
             for centre in treating_centres(source)
             for code in sorted(C.TREATMENT_LABELS)
             if not int(((df["center"] == centre) & (df[C.TREATMENT] == code)).sum())]
    if empty:
        bad.append(
            "P2  " + "; ".join(f"{centre} contributes no {arm} patient" for centre, arm in empty)
            + ". Restriction 1 guarantees a bridging arm at every retained centre, so this can only "
            "be restriction 2 having removed a centre's whole control arm: P(IVT = 0 | centre, X) = 0 "
            "there, the mirror image of the [§3] restriction-1 violation, and [§3] names only the "
            "bridging direction. Dropping the centre would extend [§3]; keeping it puts a "
            "structurally non-positive stratum in the propensity model, and `center` is a [§6] "
            "covariate. Amend [§3] with the PI; do not choose here.")

    if not len(df) or not treating_centres(df):
        bad.append(
            f"P3  the cohort is {len(df)} record(s) at {len(treating_centres(df))} centre(s). An "
            "empty cohort is a bug, not a result: every stage after this one divides by zero or fits "
            "over no data, and the first place it surfaces is a LinAlgError in Stage 6 with nothing "
            "pointing back here.")

    if removed + len(df) != len(source):
        bad.append(
            f"P4  the flow does not reconcile: {removed} removed + {len(df)} retained != "
            f"{len(source)} given. A record selected by NEITHER restriction's keep-mask leaves every "
            "table internally consistent [§3.1]; only a sum against the input frame notices.")

    if bad:
        raise C.SchemaError(
            "\n  " + "\n  ".join(bad)
            + f"\n  {len(bad)} cohort postcondition(s) failed on {len(df)} retained of "
            f"{len(source)} record(s).")


# --- the public entry point ------------------------------------------------------------------------
#
# The restrictions are applied 1 then 2, as [§3] lists them and as roadmap Stage 5 requires.
#
#   1 then 2   (the declared order)      126 → 104 → 93      22 removed, then 11
#   2 then 1                             126 → 107 → 93      19 removed, then 14
#
# The resulting COHORT does not depend on the order — both are row masks over independent predicates,
# so the retained set is an intersection and intersection is commutative; verified on v7 rather than
# reasoned, identifier for identifier. The ATTRIBUTION does depend on it, and that is what the flow
# table publishes: a reader of it is reading a decomposition, not a derivation. Eleven patients are
# both ineligible and at a never-treating centre, and they are counted once, under restriction 1.
#
# That is the line an implementer "simplifying the two masks into one" deletes.

def build(df: pd.DataFrame, audit: Audit) -> pd.DataFrame:
    """The [§3] primary cohort: two restrictions in the declared order, then Stage 3's cohort-wise
    derivations. Returns a new frame of 33 columns; appends three `cohort` entries, one
    `missingness` entry and — through derive_cohort — two `derivation` entries.

    **The only caller of derive.derive_cohort in the repository**, which is what makes the [§13]
    median a property of the cohort rather than of whichever frame a driver happened to pass. The
    median is 5.0 mL here and 6.0 mL on the frame this receives [§4.5].

    Raises SchemaError on C1-C4 before any row is removed, and on P1-P4 after the cohort is built.
    """
    _assert_cohort_inputs(df)
    source = df
    df = df.copy()

    keep = treating_centres(df)
    dropped = df[~df["center"].isin(keep)]
    _record_removal(audit, "restrict_centres", dropped, _centres_detail(keep, dropped),
                    _centres_table(df, keep, dropped))
    after_1 = df[df["center"].isin(keep)]

    ineligible = after_1[~eligibility.retained(after_1)]
    _record_removal(audit, "restrict_eligibility", ineligible,
                    _eligibility_detail(after_1, ineligible), _eligibility_table(after_1, keep))
    after_2 = after_1[eligibility.retained(after_1)]

    audit.record("cohort", "cohort_flow", len(after_2),
                 _flow_detail(source, after_2), table=_flow_table(source, after_1, after_2))

    #   126 records, as classified   6.0 mL
    #   104, after restriction 1     5.0 mL   <-- equal to the cohort's on v7, so the median
    #    93, the cohort              5.0 mL       CANNOT tell this position from the right one.
    #                                             The audit entry's `records` cell can: 93. [§4.5]
    out = derive.derive_cohort(after_2, audit)
    absence_by_column(out, audit, [*sorted(C.ANALYSIS_NAMES), *sorted(C.DERIVED_NAMES)],
                      "absence_by_cohort_column", _ABSENCE_DETAIL)
    _assert_cohort(out, source, len(dropped) + len(ineligible))
    return out
