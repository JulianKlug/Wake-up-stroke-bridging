"""Acceptance tests for Stage 5 — §12 of `specs/stage5_cohort_construction.md`.

Tests that need the private workbook are marked `skipif(not DATA_XLSX.exists())` and take this file's
**own** module-scoped `workbook` fixture; every other test runs on any checkout, against
`tests/fixture_schema.xlsx` or a hand-built frame.

This file is **not** exempt from the Stage 1 §7 raw-name scan and must not become exempt. Every frame
below is built with analysis names.

**Neither existing hand-built input can reach the end of this stage**, and that is measured rather than
anticipated (§18)::

                        records   restriction 1     restriction 2   per retained centre (t, c)
  fixture_schema.xlsx      2      keeps HUG only          1         HUG (1, 0)      → P2 RAISES
  hand_frame()             6      keeps HUG, CHUV         3         HUG (2, 0)      → P2 RAISES
                                                                    CHUV (1, 0)

Both collapse to a cohort with no control arm, because both were built for stages that never removed a
row. Neither is wrong; neither is a cohort. So this file declares `cohort_frame()` — `hand_frame()`
plus three records — and the two existing frames become *positive* tests of the both-arms guard
(§12.4) rather than fixtures anything runs on. Its verified classification and flow (§12.0, §18)::

    HAND-1   HUG     bridging  flag 0  no reason   → eligible     removed by nothing
    HAND-2   CHUV    bridging  flag 0  no reason   → eligible     removed by nothing
    HAND-3   Lugano  control   flag 0  reason      → eligible     removed by restriction 1
    HAND-4   USZ     control   flag 0  reason      → eligible     removed by restriction 1
    HAND-5   HUG     bridging  flag 0  no reason   → eligible     removed by nothing
    HAND-6   Lugano  control   flag 0  no reason   → eligible     removed by restriction 1
    COHORT-1 HUG     control   flag 0  reason      → eligible     RETAINED
    COHORT-2 HUG     control   flag 1  reason      → ineligible   removed by restriction 2
    COHORT-3 CHUV    control   flag 0  no reason   → eligible     RETAINED, undocumented

    9 records → restriction 1 keeps HUG and CHUV, removing 3 → restriction 2 removes 1 → 5
    cohort: HUG (2 bridging, 1 control), CHUV (1 bridging, 1 control)
    core_ml median: 10.0 over the 9, 10.0 after restriction 1, 8.0 over the cohort

**Bare frame or normalised frame — the distinction every test below turns on** (§12.0.1)::

    cohort_frame()                       raw codes    → C4 fires: no centre is in CENTER_ORDER
    run(cohort_frame(), hand_source(…))  labels       → the frame the restrictions are written for

Unless a test says otherwise, every frame below is the normalised one: `run(...)`, then `derive`, then
`classify`. The deliberate exceptions are §12.3's C1+C4 case, which uses the bare frame precisely
because it is one, and §12.11's, which stops after `classify`. A test that means to reach P2 and passes
a bare frame gets a C4 raise instead and proves nothing — the failure looks like a passing
`pytest.raises(SchemaError)` and the message is never read.

Two of the four preconditions guard against inputs the earlier stages refuse to produce, so the frames
that reach them are corrupted **after `classify`** rather than built by `run()`::

    C1  eligibility column absent      unreachable through classify(); drop the column after it
    C2  core_above_median present      reachable: call build twice
    C3  ivt = 2 / ivt missing          Stage 2's A5 and A4 raise first, and so does Stage 4's E5
    C4  centre = "Bern" / missing      Stage 2's A3 raises first, and so does Stage 4's E4
    C1 + C4 together                   the BARE cohort_frame(): raw codes, never classified
"""
from __future__ import annotations

import ast
import os
import subprocess
import sys
import warnings
from pathlib import Path

import pandas as pd
import pytest

import cohort
import config
import data
import derive
import eligibility
from test_data import ANALYSIS_DTYPES, CENTER_CODES, hand_frame, hand_source, run

DATA_GATED = pytest.mark.skipif(
    not config.DATA_XLSX.exists(),
    reason="the private workbook is gitignored and absent from this checkout")

MODULE = Path(cohort.__file__).resolve()
MODULE_DIR = MODULE.parent          # the flat module root
TESTS_DIR = Path(__file__).resolve().parent   # this file's own directory

CLASSES = config.ELIGIBILITY_ORDER
ARMS = tuple(config.TREATMENT_LABELS[code] for code in sorted(config.TREATMENT_LABELS))

# The analysis name of the column DECISION 1a stopped reading. Named once, so §12.6's scan and the
# frames below cannot drift apart, and so this file states which column the rule is about.
REASON = "contraindication_reason"

COHORT_IDS = {"HAND-1", "HAND-2", "HAND-5", "COHORT-1", "COHORT-3"}


# --- 12.0  the frame this stage is tested on -------------------------------------------------------

_COHORT_RECORDS: list[dict[str, object]] = [
    # a control at a treating centre with a DOCUMENTED reason — without this one, no centre has
    # both arms
    dict(case_id="COHORT-1", center=CENTER_CODES[0], ivt=0, ivt_contraindicated=0,
         contraindication_reason="Clinician decision", core_ml=8.0, tmax6_ml=50.0,
         penumbra_ml=42.0, onset_to_ivt_min=None, onset_to_groin_min=330.0),
    # the only ineligible record anywhere in the file: restriction 2 removes exactly this one
    dict(case_id="COHORT-2", center=CENTER_CODES[0], ivt=0, ivt_contraindicated=1,
         contraindication_reason="Anticoagulation", core_ml=12.0, tmax6_ml=70.0,
         penumbra_ml=58.0, onset_to_ivt_min=None, onset_to_groin_min=280.0),
    # a control with NO documented reason that SURVIVES both restrictions. Under DECISION 1a this
    # is the record that proves an undocumented patient is eligible rather than a third class, and
    # it is the one §7.3's last row counts. It is deliberately at a DIFFERENT centre from COHORT-1,
    # so the flow table's `of which` row reads 1 centre against 2 in the row above it
    dict(case_id="COHORT-3", center=CENTER_CODES[1], ivt=0, ivt_contraindicated=0,
         contraindication_reason=None, core_ml=3.0, tmax6_ml=40.0,
         penumbra_ml=37.0, onset_to_ivt_min=None, onset_to_groin_min=420.0),
]


def cohort_frame(**overrides: object) -> pd.DataFrame:
    """`hand_frame()` plus three records, so that a cohort with both arms exists. §12.0.

    Built as one DataFrame from records rather than by `pd.concat`, which emits a FutureWarning on
    pandas 2.3.3 when an operand carries an all-NA column — `onset_to_ivt_min` is None on all three
    records above — and whose dtype resolution for that case is documented as changing.

    `HAND-4` is the template because it is the one record whose volumes are internally consistent and
    which carries no correction: `HAND-1` carries the copy-paste penumbra that Stage 2 rewrites, and
    `HAND-2` the 999 groin time that raises an observation.
    """
    base = hand_frame()
    template = base[base["case_id"] == "HAND-4"].iloc[0].to_dict()
    records = [*base.to_dict("records"), *({**template, **r} for r in _COHORT_RECORDS)]
    df = pd.DataFrame(records)
    for column, dtype in ANALYSIS_DTYPES.items():
        df[column] = df[column].astype(dtype)
    for column, value in overrides.items():
        df[column] = value
    return df[sorted(config.ANALYSIS_NAMES)]


_SOURCE = hand_source(Path("."), "hand_built_stage5")


def classified(df: pd.DataFrame | None = None) -> tuple[pd.DataFrame, data.Audit]:
    """A **bare** frame put through Stage 2, Stage 3 and Stage 4 — the frame `build` is written for."""
    frame, audit = run(cohort_frame() if df is None else df, _SOURCE)
    return eligibility.classify(derive.derive(frame, audit), audit), audit


def built(df: pd.DataFrame | None = None) -> tuple[pd.DataFrame, data.Audit]:
    frame, audit = classified(df)
    return cohort.build(frame, audit), audit


def set_cell(df: pd.DataFrame, case_id: str, column: str, value: object) -> pd.DataFrame:
    """One cell replaced on a frame that has already been through a stage. Matched on `case_id`."""
    df = df.copy()
    df.loc[df["case_id"] == case_id, column] = value
    return df


def ids(df: pd.DataFrame) -> set[str]:
    return set(df["case_id"])


def message_of(excinfo) -> str:
    return str(excinfo.value)


def rows_of(entry: data.AuditEntry) -> dict[str, tuple[str, ...]]:
    """A table's body keyed by its first cell, so a test names the row it is asserting."""
    return {row[0]: row[1:] for row in entry.table[1:]}


def entry_of(audit: data.Audit, step: str) -> data.AuditEntry:
    found = audit.entry("cohort", step) or audit.entry("derivation", step) or audit.entry(
        "missingness", step)
    assert found is not None, step
    return found


def test_the_frame_this_file_declares_passes_every_stage_2_assertion():
    # A1-A9 green, so any SchemaError a test below sees comes from what that test did.
    frame, _ = run(cohort_frame(), _SOURCE)
    assert len(frame) == 9
    assert list(frame.columns) == sorted(config.ANALYSIS_NAMES)


def test_its_dtypes_are_the_hand_frames_and_construction_is_warning_free():
    # The `pd.concat` form is not: it emits a FutureWarning about all-NA columns on pandas 2.3.3, and
    # the dtype resolution it warns about is what `onset_to_ivt_min` would be resolved by.
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        frame = cohort_frame()
    assert frame.dtypes.to_dict() == hand_frame().dtypes.to_dict()


# --- 12.1  the two restrictions -------------------------------------------------------------------

def test_the_cohort_is_the_five_retained_records_by_identifier():
    # By identifier, never by count: a cohort that removed COHORT-1 and kept COHORT-2 — the inversion
    # §12.6 exists to prevent — would leave the count at 5.
    out, _ = built()
    assert ids(out) == COHORT_IDS


def test_treating_centres_is_the_bridging_arms_centres_in_center_order():
    frame, _ = classified()
    assert cohort.treating_centres(frame) == ("HUG", "CHUV")


def test_the_order_is_center_orders_and_not_the_frames():
    # Reversing the rows reverses the order the centres first appear in, and must not move the result.
    frame, _ = classified()
    assert cohort.treating_centres(frame.iloc[::-1]) == ("HUG", "CHUV")


def test_a_centre_contributing_no_bridging_patient_is_not_returned():
    frame, _ = classified()
    # Lugano contributes two controls and no treated patient; USZ one control. Both are dropped, and
    # the difference between "contributed none" and "contributed no rows" is §7.1's all-zero row.
    assert "Lugano" not in cohort.treating_centres(frame)
    assert "USZ" not in cohort.treating_centres(frame)
    assert int((frame["center"] == "Lugano").sum()) == 2


def test_restriction_1_removes_records_at_more_than_one_centre():
    _, audit = built()
    entry = entry_of(audit, "restrict_centres")
    assert entry.n == 3
    assert set(entry.case_ids) == {"HAND-3", "HAND-4", "HAND-6"}


def test_restriction_2_removes_exactly_the_flagged_control():
    _, audit = built()
    entry = entry_of(audit, "restrict_eligibility")
    assert entry.n == 1 and entry.case_ids == ("COHORT-2",)


def test_no_ineligible_patient_remains():
    # Roadmap invariant 2, read off the returned frame rather than off the log.
    out, _ = built()
    assert int((out[config.ELIGIBILITY] == config.INELIGIBLE).sum()) == 0


# --- 12.2  an undocumented patient is eligible, and the reason column is not read [DECISION 1a] ----

def test_an_undocumented_control_is_eligible_and_in_the_cohort():
    out, _ = built()
    undocumented = out.loc[out["case_id"] == "COHORT-3"].iloc[0]
    assert pd.isna(undocumented[REASON])
    assert undocumented[config.ELIGIBILITY] == config.ELIGIBLE
    # And the record that differs from it only in carrying a reason is in the cohort too. Both halves,
    # because either alone passes against a broken restriction — the first against one that retains
    # everything, the second against one that still reads the reason and drops the blank.
    assert "COHORT-1" in ids(out)


def test_blanking_the_reason_column_on_every_record_leaves_the_cohort_unmoved():
    """The whole of DECISION 1a in one assertion. Under DECISION 1 this moved three records into
    `indeterminate`; under DECISION 1a the column contributes nothing to any class, so the cohort
    cannot move — identifier for identifier, not merely in size."""
    baseline, _ = built()
    blanked, _ = built(cohort_frame(contraindication_reason=None))
    assert ids(blanked) == ids(baseline) == COHORT_IDS


# --- 12.3  C1-C4 fire, naming their cases ---------------------------------------------------------

def test_c1_fires_on_a_frame_classify_has_not_run_on_and_leaves_the_audit_untouched():
    frame, audit = classified()
    before = len(audit.entries)
    with pytest.raises(config.SchemaError) as e:
        cohort.build(frame.drop(columns=[config.ELIGIBILITY]), audit)
    assert "C1" in message_of(e) and config.ELIGIBILITY in message_of(e)
    # The half C1 exists for: without it, restriction 1 has already removed three patients and
    # recorded an entry by the time `retained` raises its KeyError.
    assert len(audit.entries) == before


def test_c2_fires_when_build_is_called_twice_and_leaves_the_audit_untouched():
    # Also the shape a bootstrap replicate would take, which is why it is C2 and not derive_cohort's.
    out, audit = built()
    before = len(audit.entries)
    with pytest.raises(config.SchemaError) as e:
        cohort.build(out, audit)
    assert "C2" in message_of(e) and config.COHORT_DEPENDENT_SUBGROUPS[0] in message_of(e)
    assert len(audit.entries) == before


@pytest.mark.parametrize("value", [2, pd.NA])
def test_c3_fires_on_an_exposure_outside_zero_and_one(value):
    # Corrupted after `classify`: Stage 2's A5 and A4 raise on both values, and so does Stage 4's E5.
    frame, audit = classified()
    with pytest.raises(config.SchemaError) as e:
        cohort.build(set_cell(frame, "HAND-1", config.TREATMENT, value), audit)
    assert "C3" in message_of(e) and "HAND-1" in message_of(e)


@pytest.mark.parametrize("value", ["Bern", pd.NA])
def test_c4_fires_on_a_centre_outside_center_order(value):
    """Injected after `classify`, and the reason is stronger than C3's: two earlier checks fire on the
    same cell first, in two different modules — Stage 2's A3 rejects any centre outside
    CENTER_RECODE.values() and Stage 4's E4 any centre outside CENTER_ORDER. A test that corrupts
    early passes while asserting nothing about this stage.

    The missing case is the one that pins §3.1's second fact: `isin` returns False for <NA>, so
    `~isin` catches the row rather than letting it through.
    """
    frame, audit = classified()
    with pytest.raises(config.SchemaError) as e:
        cohort.build(set_cell(frame, "HAND-1", "center", value), audit)
    assert "C4" in message_of(e) and "HAND-1" in message_of(e)


def test_the_isin_mask_catches_a_missing_centre_rather_than_propagating():
    # The cell of §3.1's second fact C4 rests on, asserted rather than assumed.
    frame, _ = classified()
    frame = set_cell(frame, "HAND-1", "center", pd.NA)
    caught = ~frame["center"].isin(config.CENTER_ORDER)
    assert int(caught.sum()) == 1 and int(caught.isna().sum()) == 0


def test_a_frame_tripping_two_checks_reports_both_in_one_message():
    # The BARE frame is the cheapest such frame: no `eligibility` column, and no centre in
    # CENTER_ORDER, so C1 and C4 fire together.
    audit = data.Audit(_SOURCE)
    with pytest.raises(config.SchemaError) as e:
        cohort.build(cohort_frame(), audit)
    message = message_of(e)
    assert "C1" in message and "C4" in message
    assert "2 cohort assertion(s) failed" in message
    assert len(audit.entries) == 0


# --- 12.4  P1-P4 fire, and the two existing frames are what fire P2 -------------------------------

def _fixture_through_stage_4() -> tuple[pd.DataFrame, data.Audit]:
    df, audit = data.load(data.FIXTURE)
    return eligibility.classify(derive.derive(df, audit), audit), audit


def test_p2_fires_on_the_schema_fixture_naming_the_centre_and_the_arm():
    # Data-gate-free: the fixture is on every checkout. It restricts to one bridging patient at HUG.
    frame, audit = _fixture_through_stage_4()
    with pytest.raises(config.SchemaError) as e:
        cohort.build(frame, audit)
    message = message_of(e)
    assert "P2" in message and "HUG" in message and config.TREATMENT_LABELS[0] in message


def test_p2_on_the_hand_frame_names_both_centres_rather_than_the_first():
    # The NORMALISED hand frame (§12.0.1): the bare one raises C4 long before P2 and the test would
    # pass on the wrong exception. Its cohort is HUG (2, 0) and CHUV (1, 0), so the message enumerates.
    frame, audit = classified(hand_frame())
    with pytest.raises(config.SchemaError) as e:
        cohort.build(frame, audit)
    message = message_of(e)
    assert "P2" in message and "HUG" in message and "CHUV" in message


def test_p2_fires_when_restriction_2_removes_a_centres_whole_control_arm():
    """The failure §5.2 is written for, on a frame that is otherwise a valid cohort: HUG's only
    remaining control is the one restriction 2 removes."""
    frame = cohort_frame()
    frame.loc[frame["case_id"] == "COHORT-1", "ivt_contraindicated"] = 1
    stage_4, audit = classified(frame)
    with pytest.raises(config.SchemaError) as e:
        cohort.build(stage_4, audit)
    message = message_of(e)
    assert "P2" in message and "HUG" in message and config.TREATMENT_LABELS[0] in message
    # It says which restriction can have caused it, because an implementer debugging it will
    # otherwise look at the wrong one.
    assert "restriction 2" in message
    assert "Amend [§3] with the PI; do not choose here." in message


def test_p1_fires_when_restriction_2_is_inverted():
    """Driven by a locally reimplemented restriction that keeps the ineligible — never by
    monkeypatching `cohort.py`, per Stage 3 §12.3. A patched module would test the patch."""
    frame, _ = classified()
    keep = cohort.treating_centres(frame)
    without_restriction_2 = frame[frame["center"].isin(keep)]
    with pytest.raises(config.SchemaError) as e:
        cohort._assert_cohort(without_restriction_2, frame,
                              len(frame) - len(without_restriction_2))
    message = message_of(e)
    assert "P1" in message and "COHORT-2" in message


def test_p3_fires_when_restriction_1_keeps_no_centre():
    # No treated patient anywhere, so restriction 1 keeps nothing. The empty frame WALKS to P3:
    # derive_cohort returns 0 rows with its median rendering `missing`, constant_covariates reports
    # every covariate constant, and absence_by_column guards its own division (§5.1, §18).
    frame, audit = classified(cohort_frame(ivt=0, onset_to_ivt_min=None))
    with pytest.raises(config.SchemaError) as e:
        cohort.build(frame, audit)
    message = message_of(e)
    assert "P3" in message and "0 record(s) at 0 centre(s)" in message
    assert "P2" not in message                       # no retained centre for P2 to range over


def test_p4_is_the_only_thing_that_notices_a_record_selected_by_neither_mask():
    """§3.1's first fact, and the reason P4 sums against the input frame.

    The reimplementation drops by *centre membership* rather than by `~keep`, so a record whose centre
    is missing is selected by neither mask. It leaves the analysis, the flow table reconciles row by
    row, and its identifier appears nowhere. The missing centre is injected after `classify`, for
    §12.3's reason.
    """
    frame, _ = classified()
    frame = set_cell(frame, "HAND-3", "center", pd.NA)

    keep = cohort.treating_centres(frame)
    dropped_centres = tuple(c for c in config.CENTER_ORDER if c not in keep)
    dropped = frame[frame["center"].isin(dropped_centres)]
    after_1 = frame[frame["center"].isin(keep)]
    after_2 = after_1[eligibility.retained(after_1)]

    # Every row of the flow table is internally consistent, which is the whole problem.
    table = cohort._flow_table(frame, after_1, after_2)
    for row in table[1:]:
        arms, classes = row[3:3 + len(ARMS)], row[3 + len(ARMS):]
        assert sum(int(c) for c in arms) == int(row[2]), row
        assert sum(int(c) for c in classes) == int(row[2]), row

    with pytest.raises(config.SchemaError) as e:
        cohort._assert_cohort(after_2, frame, len(dropped) + len(after_1) - len(after_2))
    message = message_of(e)
    assert "P4" in message and "does not reconcile" in message
    # 2 dropped by centre membership + 1 ineligible + 5 retained is 8, against the 9 given: HAND-3 is
    # in no term of that sum, and no other check has a term to notice it with.
    assert "3 removed + 5 retained != 9 given" in message
    for other in ("P1", "P2", "P3"):
        assert other not in message, other


# --- 12.5  the order does not change the cohort, and restriction 1 is invariant to restriction 2 ---

def _forward(df: pd.DataFrame) -> pd.DataFrame:
    after_1 = df[df["center"].isin(cohort.treating_centres(df))]
    return after_1[eligibility.retained(after_1)]


def _reverse(df: pd.DataFrame) -> pd.DataFrame:
    after_2 = df[eligibility.retained(df)]
    return after_2[after_2["center"].isin(cohort.treating_centres(after_2))]


def test_the_reverse_order_gives_the_same_identifiers():
    # The same identifiers, not the same count: on v7 both orders retain the same 93 records while
    # removing 22 then 11 and 19 then 14 respectively, so a count-only assertion says nothing.
    frame, _ = classified()
    assert ids(_forward(frame)) == ids(_reverse(frame)) == COHORT_IDS


def test_restriction_1_is_invariant_to_restriction_2():
    """The property §4.5 derives from DECISION 1 and which *explains* the test above: restriction 2 can
    never remove a treated patient, so the predicate is unmoved by it. Asserted separately, because
    the two would come apart on exactly the frame Stage 4's E3 forbids and a test checking only the
    sets would not say which property had broken."""
    frame, _ = classified()
    assert cohort.treating_centres(frame) == cohort.treating_centres(
        frame[eligibility.retained(frame)])


# --- 12.6  the class labels and the reason column — three AST scans --------------------------------

def _label_literals_in_source(source: str) -> list[tuple[int, str]]:
    """Every string literal in `source` that is exactly a declared class, with its line number.

    Stage 1 9.4's technique, as Stage 4 §12.6 used it. A label typed in the module would give the class
    a second declaration, and the two would be free to drift on the one edit nothing else would catch.
    """
    labels = set(CLASSES)
    return [(node.lineno, node.value)
            for node in ast.walk(ast.parse(source))
            if isinstance(node, ast.Constant) and isinstance(node.value, str)
            and node.value in labels]


def _order_subscripts_in_source(source: str) -> list[int]:
    """Every integer subscript of a declared class tuple. Indexing a *display* order would couple the
    restriction to it, so reordering a table's columns would silently rewrite the rule."""
    registries = {"ELIGIBILITY_ORDER", "ELIGIBILITY_RETAINED"}
    found: list[int] = []
    for node in ast.walk(ast.parse(source)):
        if not isinstance(node, ast.Subscript):
            continue
        target = node.value
        name = (target.attr if isinstance(target, ast.Attribute)
                else target.id if isinstance(target, ast.Name) else None)
        if name in registries and isinstance(node.slice, ast.Constant) and isinstance(
                node.slice.value, int):
            found.append(node.lineno)
    return found


def _reason_readers(source: str) -> set[str]:
    """The functions in `source` that name `contraindication_reason`, walked function by function.

    `"<module>"` stands for a reference outside every function, which would escape a per-function
    allowlist. Comments are not in the AST and are therefore invisible here, which is correct: the rule
    is about what the code reads, and the module docstring discusses the column at length.
    """
    tree = ast.parse(source)
    def hits(node: ast.AST) -> set[int]:
        return {id(n) for n in ast.walk(node)
                if isinstance(n, ast.Constant) and n.value == REASON}

    everywhere, inside, readers = hits(tree), set(), set()
    for node in ast.walk(tree):
        if isinstance(node, ast.FunctionDef) and hits(node):
            readers.add(node.name)
            inside |= hits(node)
    if everywhere - inside:
        readers.add("<module>")
    return readers


def test_cohort_py_writes_no_class_label():
    offenders = _label_literals_in_source(MODULE.read_text(encoding="utf-8"))
    assert not offenders, f"cohort.py names a class at {offenders}"


def test_the_label_scan_actually_fires():
    snippet = f'keep = df[df[C.ELIGIBILITY] != {config.INELIGIBLE!r}]\n'
    assert _label_literals_in_source(snippet) == [(1, config.INELIGIBLE)]


def test_cohort_py_never_indexes_the_display_order():
    offenders = _order_subscripts_in_source(MODULE.read_text(encoding="utf-8"))
    assert not offenders, f"cohort.py indexes a class tuple at lines {offenders}"


def test_the_subscript_scan_actually_fires():
    assert _order_subscripts_in_source("keep = retained(df, C.ELIGIBILITY_RETAINED[0])\n") == [1]
    # and it does not fire on the legitimate keyed lookups this module makes
    assert _order_subscripts_in_source("label = C.TREATMENT_LABELS[code]\n") == []


def test_the_reason_column_reaches_only_the_three_rendering_helpers():
    """Under DECISION 1a the column may reach a *table* or a *`detail` string* and must never reach a
    mask that decides a class. Three functions, not one: §7.2's detail interpolates the undocumented
    count and the centre count from it, so an allowlist of `_flow_table` alone would fail against the
    first honest implementation.

    This is the one place Stage 5 can catch a partial revert of DECISION 1a — an implementer who
    half-remembers "a blank is never read as no contraindication" and puts `.notna()` back into the
    restriction. That edit takes the cohort from 93 to 50 on the workbook with every table still
    reconciling.
    """
    assert _reason_readers(MODULE.read_text(encoding="utf-8")) == {
        "_eligibility_detail", "_flow_detail", "_flow_table"}


def test_the_reason_scan_fires_on_a_reference_pasted_into_build():
    # `build` is the function the restrictions live in, and therefore the one place the column must
    # never appear.
    snippet = (f"def build(df, audit):\n"
               f"    keep = eligibility.retained(df) & df[{REASON!r}].notna()\n"
               f"    return df[keep]\n")
    assert _reason_readers(snippet) == {"build"}


def test_the_reason_scan_sees_a_reference_outside_every_function():
    assert _reason_readers(f"_MASK = {REASON!r}\n") == {"<module>"}


# --- 12.7  derive_cohort is called on the cohort ---------------------------------------------------

def _statistic(entry: data.AuditEntry, name: str) -> str:
    return rows_of(entry)[name][0]


def test_the_records_cell_of_the_frozen_median_entry_is_the_cohorts_size():
    """The assertion that holds on every frame, including the one where the median does not move.

    A `derive_cohort` moved to the top of `build` writes 9 here and 126 on v7; moved *between* the two
    restrictions — the likelier slip, and the position a reader tidying the function would choose — it
    writes 6 here and 104 on v7, where the median cannot tell the difference [§4.5, §21 R1].
    """
    out, audit = built()
    entry = entry_of(audit, "core_above_median")
    assert _statistic(entry, "records") == str(len(out)) == "5"


def test_the_frozen_median_is_the_cohorts_and_not_either_frames_above_it():
    out, audit = built()
    # Both cells are strings — `data._fmt` renders through %.6g, so this is "8" and never "8.0".
    assert _statistic(entry_of(audit, "core_above_median"), "median (mL)") == "8"
    assert float(out["core_ml"].median()) == 8.0
    # What it is discriminating against, stated rather than left implicit.
    frame, _ = classified()
    after_1 = frame[frame["center"].isin(cohort.treating_centres(frame))]
    assert float(frame["core_ml"].median()) == 10.0
    assert float(after_1["core_ml"].median()) == 10.0


def _derive_cohort_callers() -> set[str]:
    """Every file under `extended_bridging/` whose AST calls `derive_cohort`."""
    found: set[str] = set()
    for dirpath, dirnames, filenames in os.walk(MODULE_DIR):
        dirnames[:] = sorted(d for d in dirnames if d not in {".venv", "__pycache__"})
        for name in sorted(f for f in filenames if f.endswith(".py")):
            path = Path(dirpath) / name
            tree = ast.parse(path.read_text(encoding="utf-8"))
            for node in ast.walk(tree):
                if not isinstance(node, ast.Call):
                    continue
                target = node.func
                called = (target.attr if isinstance(target, ast.Attribute)
                          else target.id if isinstance(target, ast.Name) else None)
                if called == "derive_cohort":
                    found.add(path.name)
    return found


def test_build_is_derive_cohorts_only_caller_outside_derive_and_the_tests():
    # `derive.py`'s "Stage 5 is derive_cohort's only caller" is a sentence until something checks it.
    # `derive.py` itself only defines it, and this file reaches it exclusively through `build`, so the
    # one production call site in the repository is `cohort.build` and the one test driver is Stage 3's.
    assert _derive_cohort_callers() == {"cohort.py", "test_derive.py"}


# --- 12.8  the three tables -----------------------------------------------------------------------

def test_the_centres_table_renders_every_declared_centre():
    _, audit = built()
    entry = entry_of(audit, "restrict_centres")
    assert entry.table[0] == ("centre", *ARMS, "n", "status")
    assert len(entry.table) == len(config.CENTER_ORDER) + 2      # + header + `all`
    assert [row[0] for row in entry.table[1:]] == [*config.CENTER_ORDER, "all"]


def test_the_centres_table_is_asserted_cell_for_cell():
    _, audit = built()
    rows = rows_of(entry_of(audit, "restrict_centres"))
    assert rows["HUG"] == ("2", "2", "4", "retained")
    assert rows["CHUV"] == ("1", "1", "2", "retained")
    assert rows["Lugano"] == ("2", "0", "2", "dropped")          # contributed none — the all-zero cell
    assert rows["USZ"] == ("1", "0", "1", "dropped")
    assert rows["all"] == ("6", "3", "9", "3 removed")


def test_the_status_cells_partition_and_agree_with_treating_centres():
    frame, audit = classified()
    cohort.build(frame, audit)
    keep = cohort.treating_centres(frame)
    rows = rows_of(entry_of(audit, "restrict_centres"))
    for centre in config.CENTER_ORDER:                 # the `all` row is a total, not a status
        assert rows[centre][-1] == ("retained" if centre in keep else "dropped"), centre


def test_the_eligibility_table_ranges_over_the_retained_centres_only():
    _, audit = built()
    entry = entry_of(audit, "restrict_eligibility")
    assert entry.table[0] == ("centre", *CLASSES, "n", "removed")
    assert [row[0] for row in entry.table[1:]] == ["HUG", "CHUV", "all"]
    rows = rows_of(entry)
    assert rows["HUG"] == ("3", "1", "4", "1")
    assert rows["CHUV"] == ("2", "0", "2", "0")
    assert rows["all"] == ("5", "1", "6", "1")


def test_the_eligibility_tables_removed_column_is_its_ineligible_column():
    out, audit = built()
    entry = entry_of(audit, "restrict_eligibility")
    for label, cells in rows_of(entry).items():
        assert cells[CLASSES.index(config.INELIGIBLE)] == cells[-1], label
    # and the `all` row's counts are the entry's n and the cohort's size — the identity §7.2 notes
    assert rows_of(entry)["all"][-1] == str(entry.n)
    assert rows_of(entry)["all"][CLASSES.index(config.ELIGIBLE)] == str(len(out))


def test_the_flow_table_is_four_rows_asserted_cell_for_cell():
    _, audit = built()
    entry = entry_of(audit, "cohort_flow")
    assert entry.table[0] == ("step", "centres", "records", *ARMS, *CLASSES)
    assert len(entry.table) == 5
    rows = rows_of(entry)
    assert rows["as classified"] == ("4", "9", "6", "3", "8", "1")
    assert rows["after restriction 1"] == ("2", "6", "3", "3", "5", "1")
    assert rows["after restriction 2"] == ("2", "5", "2", "3", "5", "0")
    assert rows[f"of which {config.TREATMENT_LABELS[0]}, no reason on file"] == (
        "1", "1", "1", "0", "1", "0")


def test_every_flow_row_reconciles_against_its_own_frame():
    _, audit = built()
    for label, cells in rows_of(entry_of(audit, "cohort_flow")).items():
        arms, classes = cells[2:2 + len(ARMS)], cells[2 + len(ARMS):]
        assert sum(int(c) for c in arms) == int(cells[1]), label
        assert sum(int(c) for c in classes) == int(cells[1]), label


def test_the_last_flow_row_is_arm_restricted_and_its_centres_cell_is_the_finding():
    """Two numbers, both named, because that difference is §7.3's finding: the undocumented group is
    centre-driven rather than patient-driven, and dropping the arm restriction would inflate a cell
    nobody checks rather than fail."""
    out, audit = built()
    rows = rows_of(entry_of(audit, "cohort_flow"))
    of_which = rows[f"of which {config.TREATMENT_LABELS[0]}, no reason on file"]
    assert of_which[0] == "1" and rows["after restriction 2"][0] == "2"
    # over ALL arms the same count is 4 at 2 centres, and before restriction 1 it is 5 at 3
    all_arms = out[out[REASON].isna()]
    assert len(all_arms) == 4 and int(all_arms["center"].nunique()) == 2
    frame, _ = classified()
    unrestricted = frame[frame[REASON].isna()]
    assert len(unrestricted) == 5 and int(unrestricted["center"].nunique()) == 3


def test_the_centres_cell_is_the_frames_own_nunique_at_that_step():
    # Declared, not inferred from the table: the three candidate rules — nunique, len(CENTER_ORDER)
    # and the retained tuple's length — agree on v7 and disagree here.
    frame, audit = classified()
    out = cohort.build(frame, audit)
    rows = rows_of(entry_of(audit, "cohort_flow"))
    assert rows["as classified"][0] == str(int(frame["center"].nunique())) == "4"
    assert rows["after restriction 2"][0] == str(int(out["center"].nunique())) == "2"
    assert len(config.CENTER_ORDER) == 4 and len(cohort.treating_centres(out)) == 2


def test_the_tables_are_not_built_with_crosstab():
    # An AST scan and not a text scan, because the module's own comments explain why it is not used:
    # `pd.crosstab` drops absent combinations, and the empty cells are the information.
    called = {node.func.attr for node in ast.walk(ast.parse(MODULE.read_text(encoding="utf-8")))
              if isinstance(node, ast.Call) and isinstance(node.func, ast.Attribute)}
    assert "crosstab" not in called


# --- 12.9  the audit inventory --------------------------------------------------------------------

_STAGE_5_INVENTORY: list[tuple[str, str]] = [
    ("cohort", "restrict_centres"),
    ("cohort", "restrict_eligibility"),
    ("cohort", "cohort_flow"),
    ("derivation", "core_above_median"),
    ("derivation", "constant_covariates"),
    ("missingness", "absence_by_cohort_column"),
]


def test_the_six_entries_appear_in_the_declared_order():
    # Position against an index captured immediately before the call, never against a remembered
    # total: that is the trap TODOS.md records for test_derive.py's tail slice.
    frame, audit = classified()
    before = len(audit.entries)
    cohort.build(frame, audit)
    assert [(e.kind, e.step) for e in audit.entries[before:]] == _STAGE_5_INVENTORY


def test_the_two_removals_name_their_patients_and_the_flow_names_nobody():
    _, audit = built()
    for step in ("restrict_centres", "restrict_eligibility"):
        entry = entry_of(audit, step)
        assert len(entry.case_ids) == entry.n > 0, step
    flow = entry_of(audit, "cohort_flow")
    assert flow.n == 5 and flow.case_ids == ()


def test_record_removal_raises_when_it_cannot_name_what_it_removed():
    # Fired directly, on a duplicated identifier: A2 forbids one, so this is the only way to reach it.
    frame, audit = classified()
    duplicated = set_cell(frame, "HAND-3", "case_id", "HAND-6")
    before = len(audit.entries)
    with pytest.raises(config.SchemaError) as e:
        cohort._record_removal(audit, "restrict_centres",
                               duplicated[duplicated["center"] == "Lugano"], "detail",
                               (("centre",), ("Lugano",)))
    assert "removed 2 record(s) and names 1" in message_of(e)
    assert len(audit.entries) == before


def test_the_entries_render_under_the_cohort_construction_heading():
    _, audit = built()
    rendered = audit.to_markdown()
    assert f"## {data._HEADINGS['cohort']}" in rendered
    assert "## Cohort construction" in rendered
    for step in ("restrict_centres", "restrict_eligibility", "cohort_flow"):
        assert f"- **{step}**" in rendered


def test_data_kinds_is_the_declared_nine_in_order():
    # So a future implementer cannot quietly add a tenth. `cohort` sits between `derivation` and
    # `model` [§6.2], which is what puts "who is in the analysis" before what was fitted to them and
    # before its denominators. `model` is Stage 6's [Stage 6 §7.1].
    assert data.KINDS == (
        "provenance", "contract", "correction", "observation",
        "derivation", "cohort", "model", "structural", "missingness")
    assert data._HEADINGS["cohort"] == "Cohort construction"


def test_two_runs_render_identical_markdown():
    assert built()[1].to_markdown() == built()[1].to_markdown()


def test_the_log_is_identical_across_interpreters_with_different_hash_seeds():
    """It imports a test module in a subprocess, and that is the price of §12.0: the fixture cannot
    reach the end of this stage, so the frame has to come from Python.

    `hand_source`'s path need not exist — nothing after the read touches the filesystem, and the log
    header prints `path.name` and the label rather than opening anything. The `build` call stands on
    its own line so that a SchemaError surfaces as a traceback rather than being swallowed into an
    empty string that then compares equal across both seeds.
    """
    # Both directories, for the reason test_propensity.py's driver gives: the shipped modules are at
    # the project root and the test modules one level down in `tests/`, and a subprocess running plain
    # `python` gets neither from pytest.
    script = (f"import sys; sys.path[:0] = [{str(MODULE_DIR)!r}, {str(TESTS_DIR)!r}]\n"
              "import pathlib\n"
              "import derive, eligibility, cohort, test_data, test_cohort\n"
              "df, audit = test_data.run(test_cohort.cohort_frame(),\n"
              "                          test_data.hand_source(pathlib.Path('/tmp')))\n"
              "df = eligibility.classify(derive.derive(df, audit), audit)\n"
              "cohort.build(df, audit)\n"
              "sys.stdout.write(audit.to_markdown())\n")
    rendered = []
    for seed in ("0", "1"):
        result = subprocess.run(
            [sys.executable, "-c", script],
            cwd=MODULE_DIR, capture_output=True, text=True,
            env={**os.environ, "PYTHONHASHSEED": seed})
        assert result.returncode == 0, result.stderr
        rendered.append(result.stdout)
    assert rendered[0] == rendered[1]
    # and the driver really did build a cohort, so the comparison is of something
    assert "## Cohort construction" in rendered[0]
    assert "restrict_centres" in rendered[0] and "cohort_flow" in rendered[0]


# --- 12.10  the cohort's denominators -------------------------------------------------------------

def test_the_cohort_missingness_entry_is_over_the_cohort_and_not_the_frame():
    out, audit = built()
    entry = entry_of(audit, "absence_by_cohort_column")
    assert entry.kind == "missingness"
    assert entry.n == len(out) == 5                   # not 9, which is the frame `build` was given
    assert [row[0] for row in entry.table[1:]] == [
        *sorted(config.ANALYSIS_NAMES), *sorted(config.DERIVED_NAMES)]
    # core_above_median could only be produced after derive_cohort, so this row pins the call order a
    # second way.
    assert config.COHORT_DEPENDENT_SUBGROUPS[0] in [row[0] for row in entry.table[1:]]


def test_its_per_centre_columns_come_from_center_order():
    # So the dropped centre renders as an all-zero column rather than vanishing.
    _, audit = built()
    entry = entry_of(audit, "absence_by_cohort_column")
    assert entry.table[0] == ("column", "kind", "n", "n_absent", "pct", *config.CENTER_ORDER)


# --- 12.11  Stage 5 needs Stage 4 and nothing else new --------------------------------------------

def test_build_on_a_classified_but_not_derived_frame_raises_key_error_on_onset_type():
    """It gets all the way to `derive_cohort` and raises there — and what it raises is
    `KeyError('onset_type')`, not `SchemaError`: `_core_above_median` succeeds because `core_ml` is a
    Stage 2 column, and it is `constant_covariates`, ranging over `PS_COVARIATES_FULL`, that fails.

    A test written against a `SchemaError` and a helpful message would never go green, and the
    temptation would be to "fix" `build` with a precondition this stage does not want: it genuinely
    does not need `derive`, and the boundary is documented rather than enforced.
    """
    stage_2, audit = run(cohort_frame(), _SOURCE)
    classified_only = eligibility.classify(stage_2, audit)
    before = len(audit.entries)
    with pytest.raises(KeyError) as e:
        cohort.build(classified_only, audit)
    assert "onset_type" in str(e.value)
    # The raise happens late, unlike C1 and C2 — the ordinary consequence of a postcondition failing
    # after the log has been written, and the difference is why those two are placed early.
    assert [(kind, step) for kind, step in
            ((e.kind, e.step) for e in audit.entries[before:])] == _STAGE_5_INVENTORY[:4]


# --- 12.12  no bare assert, no raw header, not exempt ---------------------------------------------

def _assert_lines(source: str) -> list[int]:
    return [node.lineno for node in ast.walk(ast.parse(source)) if isinstance(node, ast.Assert)]


def test_cohort_py_contains_no_bare_assert():
    # A check written as `assert` disappears under -O. Every check in the module raises SchemaError.
    assert _assert_lines(MODULE.read_text(encoding="utf-8")) == []


def test_the_assert_scan_actually_fires():
    assert _assert_lines("def f():\n    assert True\n") == [2]


def test_neither_file_is_exempt_from_the_raw_name_scan():
    # test_config.py's 9.4 sweep covers both modules the moment they exist — asserted here too, so the
    # guarantee is visible from the stage that has to keep it, and cannot be granted quietly.
    assert "cohort.py" not in config.EXEMPT_FROM_RAW_NAME_SCAN
    assert "test_cohort.py" not in config.EXEMPT_FROM_RAW_NAME_SCAN


# --- 12.13  structural facts from the workbook ----------------------------------------------------

@pytest.fixture(scope="module")
def workbook() -> tuple[pd.DataFrame, data.Audit]:
    """`load → derive → classify` over v7. This file's own fixture: `test_data.py`'s and
    `test_eligibility.py`'s do not cross files."""
    df, audit = data.load()
    return eligibility.classify(derive.derive(df, audit), audit), audit


@DATA_GATED
def test_the_workbook_retains_three_centres(workbook):
    frame, _ = workbook
    assert cohort.treating_centres(frame) == ("HUG", "CHUV", "Lugano")


@DATA_GATED
def test_the_flow_is_126_to_104_to_93(workbook):
    frame, audit = workbook
    out = cohort.build(frame, audit)
    assert len(frame) == 126 and len(out) == 93
    assert entry_of(audit, "restrict_centres").n == 22
    assert entry_of(audit, "restrict_eligibility").n == 11
    assert int((out[config.TREATMENT] == 1).sum()) == 39
    assert int((out[config.TREATMENT] == 0).sum()) == 54


@DATA_GATED
def test_the_classification_is_107_to_19_and_93_to_0_in_the_cohort(workbook):
    frame, audit = workbook
    out = cohort.build(frame, audit)
    counts = frame[config.ELIGIBILITY].value_counts()
    assert (int(counts[config.ELIGIBLE]), int(counts[config.INELIGIBLE])) == (107, 19)
    assert int((out[config.ELIGIBILITY] == config.ELIGIBLE).sum()) == 93
    assert int((out[config.ELIGIBILITY] == config.INELIGIBLE).sum()) == 0


@DATA_GATED
def test_the_retained_set_is_what_decision_1s_three_class_rule_produced(workbook):
    """The assertion that pins DECISION 1a as a relabelling of the cohort rather than a change to it.

    DECISION 1 retained a patient unless they were `ineligible`, and the 43 it labelled
    `indeterminate` — flag 0 with no documented reason — were retained. Reimplemented locally, over
    both readings of the reason column, the primary cohort is the same 93 patients (§18, §20).
    """
    frame, audit = workbook
    out = cohort.build(frame, audit)
    under_decision_1 = frame["ivt_contraindicated"] != 1        # every class but ineligible
    after_1 = frame[frame["center"].isin(cohort.treating_centres(frame))]
    reproduced = after_1[under_decision_1[after_1.index]]
    assert set(reproduced["case_id"]) == set(out["case_id"])
    assert len(reproduced) == 93


@DATA_GATED
def test_43_retained_controls_carry_no_documented_reason_at_two_centres(workbook):
    frame, audit = workbook
    out = cohort.build(frame, audit)
    controls = out[out[config.TREATMENT] == 0]
    undocumented = controls[controls[REASON].isna()]
    assert len(undocumented) == 43 and len(controls) == 54
    assert int(undocumented["center"].nunique()) == 2
    # over ALL arms the same count is 82 at 3 centres, because no treated patient carries a reason
    all_arms = out[out[REASON].isna()]
    assert len(all_arms) == 82 and int(all_arms["center"].nunique()) == 3


@DATA_GATED
def test_the_three_tables_reproduce_the_specification_cell_for_cell(workbook):
    frame, audit = workbook
    cohort.build(frame, audit)
    centres = rows_of(entry_of(audit, "restrict_centres"))
    assert centres["HUG"] == ("22", "30", "52", "retained")
    assert centres["CHUV"] == ("14", "7", "21", "retained")
    assert centres["Lugano"] == ("29", "2", "31", "retained")
    assert centres["USZ"] == ("22", "0", "22", "dropped")       # the [§3] restriction-1 finding
    assert centres["all"] == ("87", "39", "126", "22 removed")

    eligible = rows_of(entry_of(audit, "restrict_eligibility"))
    assert eligible["HUG"] == ("41", "11", "52", "11")
    assert eligible["CHUV"] == ("21", "0", "21", "0")
    assert eligible["Lugano"] == ("31", "0", "31", "0")
    assert eligible["all"] == ("93", "11", "104", "11")
    assert "USZ" not in eligible

    flow = rows_of(entry_of(audit, "cohort_flow"))
    assert flow["as classified"] == ("4", "126", "87", "39", "107", "19")
    assert flow["after restriction 1"] == ("3", "104", "65", "39", "93", "11")
    assert flow["after restriction 2"] == ("3", "93", "54", "39", "93", "0")
    assert flow[f"of which {config.TREATMENT_LABELS[0]}, no reason on file"] == (
        "2", "43", "43", "0", "43", "0")


@DATA_GATED
def test_the_frozen_median_is_5_and_the_intermediate_frames_is_also_5(workbook):
    """The coincidence §4.5 is written about, asserted so that a future workbook in which it stops
    being true is noticed rather than assumed — and so that the `records` cell is seen to be the guard
    that does not rest on it."""
    frame, audit = workbook
    out = cohort.build(frame, audit)
    after_1 = frame[frame["center"].isin(cohort.treating_centres(frame))]
    assert float(frame["core_ml"].median()) == 6.0
    assert float(after_1["core_ml"].median()) == 5.0            # identical to the cohort's
    assert float(out["core_ml"].median()) == 5.0
    entry = entry_of(audit, "core_above_median")
    assert _statistic(entry, "median (mL)") == "5"
    assert _statistic(entry, "records") == "93"


@DATA_GATED
def test_no_covariate_is_constant_over_the_cohort(workbook):
    # The premise Stage 6's design-matrix builder rests on: no [§6] covariate is lost to the
    # restriction. Not a property of every frame in this file — cohort_frame()'s five records leave 10
    # of the 13 constant, which is why the workbook is where this can be asserted at all.
    frame, audit = workbook
    out = cohort.build(frame, audit)
    assert derive.constant_covariates(out) == ()
    assert entry_of(audit, "constant_covariates").n == 0
    assert len(derive.constant_covariates(built()[0])) == 10


@DATA_GATED
def test_the_order_invariance_properties_hold_on_the_workbook(workbook):
    frame, _ = workbook
    forward, reverse = _forward(frame), _reverse(frame)
    assert set(forward["case_id"]) == set(reverse["case_id"])
    assert len(forward) == 93
    # the same 93 while removing 22 then 11 one way and 19 then 14 the other — which is why the
    # assertion above is on identifiers rather than on counts
    assert len(frame) - len(frame[frame["center"].isin(cohort.treating_centres(frame))]) == 22
    assert len(frame) - len(frame[eligibility.retained(frame)]) == 19
    assert cohort.treating_centres(frame) == cohort.treating_centres(
        frame[eligibility.retained(frame)])


@DATA_GATED
def test_every_centre_in_the_workbooks_cohort_is_both_armed(workbook):
    # So P2 is unreachable on v7: HUG (30, 11), CHUV (7, 14), Lugano (2, 29).
    frame, audit = workbook
    out = cohort.build(frame, audit)
    for centre, expected in (("HUG", (30, 11)), ("CHUV", (7, 14)), ("Lugano", (2, 29))):
        at = out[out["center"] == centre]
        assert (int((at[config.TREATMENT] == 1).sum()),
                int((at[config.TREATMENT] == 0).sum())) == expected, centre


@DATA_GATED
def test_every_precondition_is_unreachable_on_the_workbook(workbook):
    # C3's and C4's premises, which is why §4.4's branches are written for a workbook that has not
    # arrived rather than for this one.
    frame, _ = workbook
    assert bool(frame["center"].isin(config.CENTER_ORDER).all())
    assert bool(frame[config.TREATMENT].isin((0, 1)).all())


# --- 12.14  Stage 5 removes rows, adds one column, and changes nothing else -----------------------

def _assert_only_rows_removed(df: pd.DataFrame, audit: data.Audit) -> None:
    before = df.copy()
    out = cohort.build(df, audit)

    assert list(out.columns) == [*df.columns, *config.COHORT_DEPENDENT_SUBGROUPS]
    # preserved and never reset: a cohort record keeps the label it carried on the unrestricted frame,
    # so it can be traced back to the population Stages 12 and 13 analyse.
    assert out.index.is_unique and bool(out.index.isin(df.index).all())
    # the line that fails if any value is edited or any dtype "tidied" on the way past
    assert out.drop(columns=list(config.COHORT_DEPENDENT_SUBGROUPS)).equals(df.loc[out.index])
    # and the line that fails if df.copy() is moved below a write, which nothing else would notice
    assert df.equals(before)


def test_build_removes_rows_and_adds_one_column():
    frame, audit = classified()
    _assert_only_rows_removed(frame, audit)


@DATA_GATED
def test_build_removes_rows_and_adds_one_column_on_the_workbook(workbook):
    frame, audit = workbook
    _assert_only_rows_removed(frame.copy(), data.Audit(audit.source))
