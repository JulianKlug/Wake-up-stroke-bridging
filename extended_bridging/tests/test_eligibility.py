"""Acceptance tests for Stage 4 — §12 of `specs/stage4_eligibility_classification.md`.

Tests that need the private workbook are marked `skipif(not DATA_XLSX.exists())`; every other test
runs on any checkout, against `tests/fixture_schema.xlsx` or a hand-built frame.

This file is **not** exempt from the Stage 1 §7 raw-name scan and must not become exempt. Every frame
below is built with analysis names.

The hand-built frame is `test_data.py`'s, imported rather than re-declared, for Stage 3 §12's reason.
It needs no extending: it already carries both arms, both reason states and all four centres. Its
classification falls out of the records as they stand and is not hand-tuned::

    HAND-1  HUG     treated  flag 0  no reason        → eligible
    HAND-2  CHUV    treated  flag 0  no reason        → eligible
    HAND-3  Lugano  control  flag 0  "Anticoagul…"    → eligible
    HAND-4  USZ     control  flag 0  "Clinician d…"   → eligible
    HAND-5  HUG     treated  flag 0  no reason        → eligible
    HAND-6  Lugano  control  flag 0  no reason        → eligible   ← DECISION 1a: no reason on
                                                                    file is ELIGIBLE, not a class

    6 eligible, 0 ineligible — so every ineligible case below is made by
    corrupt("ivt_contraindicated", 1, where=…), which is also the only way to reach E3. Under
    DECISION 1 this frame gave 5 eligible and 1 indeterminate; only HAND-6's label moved.

**But it carries `center` as a raw code, not as a label**, and every per-centre assertion here turns
on that. `CENTER_CODES = tuple(config.CENTER_RECODE)` is `("1", "Lausanne", "Lugano", "USZ")`, and
Stage 2's `centre_recode` is what maps them to `CENTER_ORDER`. Two of the four spell their own labels
and two do not, which is why the **bare** frame fails E4 on exactly three of six records and renders
exactly three of six in the table — on two *different* subsets, which is the evidence that E4 and
`_crosstab`'s reconciliation are two guards and not one written twice (§12.8). Every classification
above is therefore stated for the frame put **through Stage 2**.

Four of the five assertions guard against inputs Stage 2 refuses to produce, so the frames that reach
them cannot be built by `run()` and are corrupted *after* it::

    E1  ivt_contraindicated missing    A4b raises first
    E2  ivt_contraindicated = 2        A5  raises first
    E5  ivt missing                    A4  raises first
    E5  ivt = 2                        A5  raises first

    E3  treated and flagged            reachable through run() — corrupt() alone suffices
    E4  centre outside CENTER_ORDER    reachable, and the BARE hand_frame() already is one
"""
from __future__ import annotations

import ast
import os
import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest

import config
import data
import derive
import eligibility
from test_data import corrupt, hand_frame, hand_source, run

DATA_GATED = pytest.mark.skipif(
    not config.DATA_XLSX.exists(),
    reason="the private workbook is gitignored and absent from this checkout")

MODULE = Path(eligibility.__file__).resolve()
MODULE_DIR = MODULE.parent

CLASSES = config.ELIGIBILITY_ORDER
N_CELLS = len(config.CENTER_ORDER) * len(config.TREATMENT_LABELS)


# --- helpers -----------------------------------------------------------------------------------
#
# `_stage2` puts a hand-built frame through the pipeline Stage 4 receives its input from, which
# matters for one thing in particular: `center` arrives as a raw code and Stage 2 recodes it, and
# every cell of §7.1's table is keyed on a CENTER_ORDER label.

_SOURCE = hand_source(Path("."), "hand_built_stage4")


def _stage2(df: pd.DataFrame) -> tuple[pd.DataFrame, data.Audit]:
    return run(df, _SOURCE)


def classified(df: pd.DataFrame | None = None) -> tuple[pd.DataFrame, data.Audit]:
    """The hand frame (or `df`) through Stage 2 and then through `classify`."""
    frame, audit = _stage2(hand_frame() if df is None else df)
    return eligibility.classify(frame, audit), audit


def after_stage_2(df: pd.DataFrame | None = None) -> pd.DataFrame:
    """The frame Stage 2 delivers, so a test can corrupt a cell Stage 2 would have refused."""
    return _stage2(hand_frame() if df is None else df)[0]


def where(df: pd.DataFrame, case_id: str, column: str):
    """One cell, addressed by identifier — never by row position."""
    return df.loc[df["case_id"] == case_id, column].iloc[0]


def eligibility_of(df: pd.DataFrame, case_id: str) -> str:
    return where(df, case_id, config.ELIGIBILITY)


def message_of(excinfo) -> str:
    return str(excinfo.value)


def set_cell(df: pd.DataFrame, case_id: str, column: str, value: object) -> pd.DataFrame:
    """`corrupt`, but on a frame that has already been through Stage 2. Matched on `case_id`."""
    df = df.copy()
    df.loc[df["case_id"] == case_id, column] = value
    return df


# The §4.3 mask, reimplemented locally — `eligibility.py` is never monkey-patched, per Stage 3 §12.3.
# A patched module would test the patch; these test what the shipped mask buys.

def _chain_without_the_fill(df: pd.DataFrame) -> pd.Series:
    flag = df["ivt_contraindicated"]
    out = pd.Series(config.ELIGIBLE, index=df.index, dtype="string")
    return out.mask(flag == 1, config.INELIGIBLE)


def _chain_with_the_deleted_revealed_fact_mask(df: pd.DataFrame) -> pd.Series:
    """DECISION 1's last mask, kept here only to prove its deletion changed nothing.

    Under DECISION 1 it forced every treated patient to `eligible` and moved 39 of 126 records,
    because no reason is recorded for any of them. Under DECISION 1a `flag == 0` already classifies
    all 39, so the mask is a no-op — and E3 is what keeps it one.
    """
    return eligibility._classify(df).mask(
        (df[config.TREATMENT] == 1).fillna(False), config.ELIGIBLE)


# --- 12.1  the two cases of the rule -------------------------------------------------------------

def test_every_treated_record_is_eligible_from_the_flag_alone():
    out, _ = classified()
    for case_id in ("HAND-1", "HAND-2", "HAND-5"):
        assert eligibility_of(out, case_id) == config.ELIGIBLE, case_id


def test_a_control_with_a_recorded_reason_is_eligible():
    # Documented: a reason was recorded and it is not an absolute contraindication.
    out, _ = classified()
    assert eligibility_of(out, "HAND-3") == config.ELIGIBLE
    assert eligibility_of(out, "HAND-4") == config.ELIGIBLE


def test_the_flag_overrides_a_recorded_reason():
    # The mask order of §4.4: flagged-and-documented is the ordinary case, not an edge one — on v7
    # all 19 flagged records carry a reason.
    out, _ = classified(corrupt("ivt_contraindicated", 1, where="HAND-3"))
    assert eligibility_of(out, "HAND-3") == config.INELIGIBLE
    assert eligibility_of(out, "HAND-4") == config.ELIGIBLE      # only the corrupted record moved


def test_the_hand_frame_classifies_six_eligible_and_none_ineligible():
    out, _ = classified()
    assert out[config.ELIGIBILITY].value_counts().to_dict() == {config.ELIGIBLE: 6}


# --- 12.2  no documented reason is ELIGIBLE, and the reason column is not read [DECISION 1a] -------

def test_a_record_with_no_documented_reason_is_eligible():
    # HAND-6 has flag 0 and no reason. Under DECISION 1 it was `indeterminate`; under DECISION 1a it
    # is eligible, and [§3]'s amendment of 2026-08-13 is the protocol change that says so.
    out, _ = classified()
    assert eligibility_of(out, "HAND-6") == config.ELIGIBLE


def test_the_reason_column_feeds_no_class_at_all():
    """The whole of DECISION 1a in one assertion: blank the reason column on EVERY record and the
    classification does not move.

    Under DECISION 1 this moved three records into `indeterminate`. It is the test that fails if an
    implementer half-remembers "a blank is never read as no contraindication" and puts `.notna()`
    back into the classifier — an edit that would silently relabel 43 of 126 records on v7.
    """
    baseline, _ = classified()
    blanked, _ = classified(hand_frame(contraindication_reason=None))
    pd.testing.assert_series_equal(
        baseline[config.ELIGIBILITY], blanked[config.ELIGIBILITY], check_dtype=False)


# --- 12.3  E2, E3, E4 and E5 fire, naming their cases ----------------------------------------------

def test_a_treated_record_flagged_contraindicated_raises_naming_it():
    with pytest.raises(config.SchemaError) as e:
        classified(corrupt("ivt_contraindicated", 1, where="HAND-1"))
    message = message_of(e)
    assert "E3" in message and "HAND-1" in message
    # It refuses to choose rather than resolving, as Stage 3's onset assertion does.
    assert "do not choose a class here" in message


def test_the_message_enumerates_the_offenders_rather_than_counting_them():
    # Every treated record flagged at once. A message reporting only a count would pass the test
    # above and fail this one.
    with pytest.raises(config.SchemaError) as e:
        classified(hand_frame(ivt_contraindicated=1))
    message = message_of(e)
    for case_id in ("HAND-1", "HAND-2", "HAND-5"):
        assert case_id in message, case_id


def test_a_flagged_control_does_not_raise():
    # E3 is about the contradiction, not about the flag. A flagged control is ordinary [§3].
    assert eligibility._assert_classifier_inputs(
        after_stage_2(corrupt("ivt_contraindicated", 1, where="HAND-3"))) is None


def test_a_third_value_in_the_flag_raises_e2():
    # Corrupted after Stage 2: A5 raises on a flag of 2 before _classify is reached.
    df = set_cell(after_stage_2(), "HAND-4", "ivt_contraindicated", 2)
    with pytest.raises(config.SchemaError) as e:
        eligibility._assert_classifier_inputs(df)
    message = message_of(e)
    assert "E2" in message and "HAND-4" in message


def test_a_centre_outside_center_order_raises_e4():
    df = set_cell(after_stage_2(), "HAND-4", "center", "Bern")
    with pytest.raises(config.SchemaError) as e:
        eligibility._assert_classifier_inputs(df)
    message = message_of(e)
    assert "E4" in message and "HAND-4" in message


def test_a_missing_centre_is_caught_too_because_isin_returns_false_for_na():
    # §3.1's third fact, asserted rather than assumed: `isin` answers a membership question and a
    # missing value is not a member, so it returns False rather than propagating — and `~isin`
    # therefore *catches* the missing row. This is why E4 carries no `.fillna(False)`.
    df = set_cell(after_stage_2(), "HAND-2", "center", pd.NA)
    with pytest.raises(config.SchemaError) as e:
        eligibility._assert_classifier_inputs(df)
    assert "E4" in message_of(e) and "HAND-2" in message_of(e)


def test_the_isin_masks_carry_no_na_on_either_dtype_they_are_used_over():
    """The property E4's and E5's correctness rests on, pinned rather than left in a comment.

    It is why neither carries a `.fillna(False)` and why writing one would be a no-op rather than the
    bug the spec's §4.2 claims: `~isin` never propagates. The dtype differs by column and neither is
    nullable in practice — numpy `bool` over the `string` centre, pandas `boolean` over the Int64
    exposure. Re-check on any pandas major bump: if `isin` ever starts propagating, E4 and E5 stop
    catching a missing value and this test is what says so.
    """
    df = set_cell(after_stage_2(), "HAND-5", config.TREATMENT, pd.NA)
    df = set_cell(df, "HAND-2", "center", pd.NA)

    arm = ~df[config.TREATMENT].isin((0, 1))
    centre = ~df["center"].isin(config.CENTER_ORDER)
    assert not pd.Series(arm).isna().any()
    assert not pd.Series(centre).isna().any()
    # and the missing rows are the ones selected, not the ones passed over
    assert set(df.loc[arm, "case_id"]) == {"HAND-5"}
    assert set(df.loc[centre, "case_id"]) == {"HAND-2"}


def test_an_exposure_outside_zero_one_raises_e5():
    df = set_cell(after_stage_2(), "HAND-5", config.TREATMENT, 2)
    with pytest.raises(config.SchemaError) as e:
        eligibility._assert_classifier_inputs(df)
    assert "E5" in message_of(e) and "HAND-5" in message_of(e)


def test_a_missing_exposure_raises_e5():
    # The one a `.fillna(False)` on E5's mask would silently absorb. Both readers of the exposure —
    # E3 and the revealed-fact mask — treat absence as "not treated", so nothing else would notice.
    df = set_cell(after_stage_2(), "HAND-5", config.TREATMENT, pd.NA)
    with pytest.raises(config.SchemaError) as e:
        eligibility._assert_classifier_inputs(df)
    assert "E5" in message_of(e) and "HAND-5" in message_of(e)


def test_several_failures_are_reported_in_one_message():
    # What the collected raise is for: a corrected workbook is diagnosed in one run, not five.
    df = set_cell(after_stage_2(), "HAND-4", "ivt_contraindicated", 2)
    df = set_cell(df, "HAND-2", "center", "Bern")
    with pytest.raises(config.SchemaError) as e:
        eligibility._assert_classifier_inputs(df)
    message = message_of(e)
    assert "E2" in message and "HAND-4" in message
    assert "E4" in message and "HAND-2" in message
    assert "2 eligibility assertion(s) failed" in message


def test_the_hand_frame_through_stage_2_trips_no_assertion():
    assert eligibility._assert_classifier_inputs(after_stage_2()) is None


# --- 12.4  E1 fires, and the fabrication square ---------------------------------------------------

def _missing_flag_frame() -> pd.DataFrame:
    """HAND-6's flag removed, after Stage 2 — A4b raises first otherwise."""
    return set_cell(after_stage_2(), "HAND-6", "ivt_contraindicated", pd.NA)


def test_a_missing_flag_raises_e1_naming_it():
    with pytest.raises(config.SchemaError) as e:
        eligibility.classify(_missing_flag_frame(), data.Audit(_SOURCE))
    message = message_of(e)
    assert "E1" in message and "HAND-6" in message
    assert "will not invent a class for it" in message


def test_without_the_assertion_and_without_the_fills_a_missing_flag_is_fabricated_as_ineligible():
    """The square's second cell, and the one Definition of done 5 is watched against.

    `ineligible` is the class Stage 5 deletes, so the failure is not a fabricated label in a table —
    it is a patient silently removed from the analysis, with a full column, plausible counts and
    every downstream denominator reconciling.
    """
    df = _missing_flag_frame()
    fabricated = _chain_without_the_fill(df)
    assert fabricated[df["case_id"] == "HAND-6"].iloc[0] == config.INELIGIBLE


def test_without_the_assertion_but_with_the_fill_it_is_eligible_and_therefore_retained():
    # The square's third cell. Still fabricated, but retained rather than deleted — so the two cells
    # fail differently, and the shipped code guards both ways: E1 raises, and the fill makes the
    # `ineligible` outcome unreachable a second time, independently of whether E1 is ever weakened.
    # Under DECISION 1 this cell read `indeterminate`; the fabricated label changed with the class
    # set, the *shape* of the square did not.
    df = _missing_flag_frame()
    with_fill = eligibility._classify(df)
    assert with_fill[df["case_id"] == "HAND-6"].iloc[0] == config.ELIGIBLE


def test_only_the_record_with_the_missing_flag_differs_between_the_two_chains():
    df = _missing_flag_frame()
    differ = eligibility._classify(df) != _chain_without_the_fill(df)
    assert list(df.loc[differ, "case_id"]) == ["HAND-6"]


# --- 12.5  the classifier is total, and its values are the declared ones ---------------------------

def _every_classifiable_frame() -> list[tuple[str, pd.DataFrame]]:
    """Every frame in this file that classifies rather than raising."""
    return [
        ("hand", classified()[0]),
        ("flagged control", classified(corrupt("ivt_contraindicated", 1, where="HAND-3"))[0]),
        ("reason added", classified(corrupt("contraindication_reason", "x", where="HAND-6"))[0]),
        ("fixture", eligibility.classify(*data.load(data.FIXTURE))),
    ]


@pytest.mark.parametrize("label, df", _every_classifiable_frame())
def test_the_column_is_never_missing_and_holds_only_declared_classes(label, df):
    column = df[config.ELIGIBILITY]
    assert int(column.isna().sum()) == 0, label
    assert set(column) <= set(CLASSES), label


@pytest.mark.parametrize("label, df", _every_classifiable_frame())
def test_the_column_is_a_string_and_not_a_categorical(label, df):
    # Stage 5 removes a whole class from the frame, and a categorical that keeps a dead level makes
    # every subsequent groupby(observed=False) resurrect `ineligible` as an all-missing row.
    assert str(df[config.ELIGIBILITY].dtype) == "string", label
    assert not isinstance(df[config.ELIGIBILITY].dtype, pd.CategoricalDtype), label


# --- 12.6  the labels come from the configuration -------------------------------------------------

def _label_literals_in_source(source: str) -> list[tuple[int, str]]:
    """Every string literal in `source` that is exactly a declared class, with its line number.

    Stage 1 9.4's technique. A label typed in the module would give the class a second declaration,
    and the two would be free to drift on the one edit — a renamed class — that nothing else in the
    repository would catch.
    """
    labels = set(CLASSES)
    return [(node.lineno, node.value)
            for node in ast.walk(ast.parse(source))
            if isinstance(node, ast.Constant) and isinstance(node.value, str)
            and node.value in labels]


def _order_subscripts_in_source(source: str) -> list[int]:
    """Every integer subscript of a declared class tuple, with its line number.

    The subtler mistake, and worse than a literal: indexing into a *display* order couples the
    classification rule to it, so reordering the columns of a table would silently rewrite the rule.
    """
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


def test_eligibility_py_writes_no_class_label():
    offenders = _label_literals_in_source(MODULE.read_text(encoding="utf-8"))
    assert not offenders, f"eligibility.py names a class at {offenders}"


def test_the_label_scan_actually_fires():
    # A scan that silently matches nothing would otherwise pass as a green test.
    snippet = f'out = pd.Series({config.ELIGIBLE!r}, index=df.index)\n'
    assert _label_literals_in_source(snippet) == [(1, config.ELIGIBLE)]


def test_eligibility_py_never_indexes_the_display_order():
    offenders = _order_subscripts_in_source(MODULE.read_text(encoding="utf-8"))
    assert not offenders, f"eligibility.py indexes a class tuple at lines {offenders}"


def test_the_subscript_scan_actually_fires():
    assert _order_subscripts_in_source("out = out.mask(cond, C.ELIGIBILITY_ORDER[2])\n") == [1]
    assert _order_subscripts_in_source("out = ELIGIBILITY_ORDER[0]\n") == [1]
    # and it does not fire on a legitimate keyed lookup by a variable
    assert _order_subscripts_in_source("label = C.TREATMENT_LABELS[code]\n") == []


# --- 12.7  `retained` is the declared set ---------------------------------------------------------

def _mixed_frame() -> pd.DataFrame:
    """Both classes present: HAND-3 flagged and therefore ineligible, the other five eligible."""
    return classified(corrupt("ivt_contraindicated", 1, where="HAND-3"))[0]


def test_retained_equals_not_ineligible():
    df = _mixed_frame()
    pd.testing.assert_series_equal(
        eligibility.retained(df), df[config.ELIGIBILITY] != config.INELIGIBLE,
        check_names=False, check_dtype=False)
    # Under DECISION 1a the two readings coincide, so this asserts an identity rather than a
    # decision. It stays because it is the identity a third class would break, and because
    # `retained` is what Stages 5, 12 and 13 all call.
    assert int(eligibility.retained(df).sum()) == 5      # every record but the flagged one


def test_retained_moves_under_a_patched_registry(monkeypatch):
    """What proves the predicate is registry-driven rather than a comparison that happens to agree.

    Under DECISION 1 the patch was `(ELIGIBLE,)` — the point where `!= ineligible` and
    `isin(RETAINED)` stopped agreeing. That IS the shipped tuple now, so the patch has to invert the
    registry instead: retained becomes `(INELIGIBLE,)` and the predicate must follow it exactly.
    """
    df = _mixed_frame()
    monkeypatch.setattr(config, "ELIGIBILITY_RETAINED", (config.INELIGIBLE,))
    patched = eligibility.retained(df)
    assert int(patched.sum()) == 1                       # only the flagged record, inverted
    assert patched[df["case_id"] == "HAND-3"].iloc[0]
    assert not patched[df["case_id"] == "HAND-6"].iloc[0]

    monkeypatch.undo()
    assert int(eligibility.retained(df).sum()) == 5
    assert eligibility.retained(df)[df["case_id"] == "HAND-6"].iloc[0]


def test_retained_is_bool_and_total():
    predicate = eligibility.retained(_mixed_frame())
    assert predicate.dtype == bool
    assert int(predicate.isna().sum()) == 0


def test_retained_raises_key_error_on_a_frame_classify_has_not_run_on():
    # The correct failure: a caller asking who is retained before anyone has been classified has a
    # sequencing bug, and a predicate that answered anyway would answer about a different population.
    with pytest.raises(KeyError):
        eligibility.retained(after_stage_2())


# --- 12.8  the cross-tabulation -------------------------------------------------------------------

def test_the_header_comes_from_the_registries():
    table = eligibility._crosstab(classified()[0])
    assert table[0] == ("centre", "arm", *CLASSES, "n")


def test_every_declared_cell_is_rendered_even_on_a_two_record_frame():
    # The fixture's 2 records occupy 2 of the 8 declared cells — one bridging at HUG, one EVT alone at
    # CHUV — so the other 6 render as zeros. The empty ones are the information: on v7 the all-zero
    # `USZ / bridging` row IS the [§3] restriction 1 finding, and pd.crosstab would omit the row and
    # the finding with it.
    table = eligibility._crosstab(eligibility.classify(*data.load(data.FIXTURE)))
    assert len(table) == N_CELLS + 1 + 1 == 10
    cells = table[1:-1]
    assert len(cells) == N_CELLS == 8
    assert sum(1 for row in cells if row[-1] == "0") == 6


def test_the_arm_order_is_control_then_treated_in_every_interpreter():
    table = eligibility._crosstab(classified()[0])
    arms = [row[1] for row in table[1:1 + len(config.TREATMENT_LABELS)]]
    assert arms == [config.TREATMENT_LABELS[0], config.TREATMENT_LABELS[1]]


def test_each_row_reconciles_to_its_own_n():
    for label, df in _every_classifiable_frame():
        for row in eligibility._crosstab(df)[1:]:
            assert sum(int(c) for c in row[2:-1]) == int(row[-1]), f"{label}: {row}"


def test_each_column_reconciles_to_the_all_both_row():
    for label, df in _every_classifiable_frame():
        table = eligibility._crosstab(df)
        cells, total = table[1:-1], table[-1]
        for j in range(2, 2 + len(CLASSES) + 1):
            assert sum(int(row[j]) for row in cells) == int(total[j]), f"{label}: column {j}"
        assert int(total[-1]) == len(df), label


def test_the_nine_rows_are_the_hand_frames_own_counts():
    # Asserted on the frame through Stage 2 — never on the bare one, whose raw centre codes would
    # render three of six records.
    table = eligibility._crosstab(classified()[0])
    assert table[1:] == (
        ("HUG",    config.TREATMENT_LABELS[0], "0", "0", "0"),
        ("HUG",    config.TREATMENT_LABELS[1], "2", "0", "2"),
        ("CHUV",   config.TREATMENT_LABELS[0], "0", "0", "0"),
        ("CHUV",   config.TREATMENT_LABELS[1], "1", "0", "1"),
        ("Lugano", config.TREATMENT_LABELS[0], "2", "0", "2"),
        ("Lugano", config.TREATMENT_LABELS[1], "0", "0", "0"),
        ("USZ",    config.TREATMENT_LABELS[0], "1", "0", "1"),
        ("USZ",    config.TREATMENT_LABELS[1], "0", "0", "0"),
        ("all",    "both",                     "6", "0", "6"))


def test_the_bare_hand_frame_renders_three_of_six_and_refuses_to_produce_the_table():
    """A positive test of §7.1's reconciliation guard, using the bare frame as a fixture for the
    failure rather than as a source of expected counts.

    Note what the `all / both` row cannot do on its own: it is computed from the frame rather than
    from the cells, so a vanished patient leaves it *correct*. That is precisely why the
    reconciliation is a check inside `_crosstab` and not merely a property a reader could notice.
    """
    bare = hand_frame()
    bare[config.ELIGIBILITY] = eligibility._classify(bare)
    with pytest.raises(config.SchemaError) as e:
        eligibility._crosstab(bare)
    assert "renders 3 of 6 record(s)" in message_of(e)


def test_e4_and_the_reconciliation_catch_different_subsets():
    """The evidence that the two guards are not one written twice, and the reason §7.1 keeps both.

    E4 catches the records whose raw centre code is not a declared label — `"1"` and `"Lausanne"`.
    The reconciliation catches the records that fall in no cell, which is the complement: the three
    whose raw codes happen to spell their own labels are the three that render. If removing E4 ever
    changes nothing, the two guards have been collapsed into one.
    """
    bare = hand_frame()
    with pytest.raises(config.SchemaError) as e:
        eligibility._assert_classifier_inputs(bare)
    named = {c for c in ("HAND-1", "HAND-2", "HAND-3", "HAND-4", "HAND-5", "HAND-6")
             if c in message_of(e)}
    assert named == {"HAND-1", "HAND-2", "HAND-5"}

    bare[config.ELIGIBILITY] = eligibility._classify(bare)
    rendered = ((bare["center"].isin(config.CENTER_ORDER))
                & (bare[config.TREATMENT].isin(tuple(config.TREATMENT_LABELS))))
    assert set(bare.loc[rendered, "case_id"]) == {"HAND-3", "HAND-4", "HAND-6"}
    assert named & set(bare.loc[rendered, "case_id"]) == set()


# --- 12.9  the audit inventory and reproduction ---------------------------------------------------

_STAGE_4_INVENTORY: list[tuple[str, str]] = [("derivation", "eligibility")]


def _pipeline() -> tuple[pd.DataFrame, data.Audit]:
    df, audit = data.load(data.FIXTURE)
    return eligibility.classify(derive.derive(df, audit), audit), audit


def test_exactly_one_entry_is_appended_and_it_is_the_declared_one():
    """Its position is asserted against the driver, never against a remembered total.

    On the fixture the index happens to be 10 — six from `load()` and four from `derive()` — but the
    count differs on the workbook, and `derive_cohort`'s two entries are Stage 5's and arrive *after*
    this one.
    """
    df, audit = data.load(data.FIXTURE)
    df = derive.derive(df, audit)
    before = len(audit.entries)
    out = eligibility.classify(df, audit)
    assert [(e.kind, e.step) for e in audit.entries[before:]] == _STAGE_4_INVENTORY

    entry = audit.entries[before]
    assert entry.n == len(out)              # the row count, not one class's count
    assert entry.case_ids == ()             # a derivation describes a column, not a patient
    assert len(entry.table) == N_CELLS + 1 + 1


def test_the_step_is_the_column_constant():
    # So the column and its log entry cannot be renamed apart.
    _, audit = _pipeline()
    assert audit.entry("derivation", config.ELIGIBILITY) is not None


def test_a_frame_with_no_ineligible_patient_does_not_raise():
    # `counts.get(cls, 0)`, never `counts[cls]`: the fixture's two records are both eligible, so
    # `value_counts()` has no `ineligible` key and indexing it would raise on the cleanest input.
    out, audit = _pipeline()
    assert set(out[config.ELIGIBILITY]) == {config.ELIGIBLE}
    assert f"0 {config.INELIGIBLE}" in audit.entry("derivation", config.ELIGIBILITY).detail


def test_the_detail_names_both_classes_the_counts_and_the_undocumented_group():
    _, audit = _pipeline()
    detail = audit.entry("derivation", config.ELIGIBILITY).detail
    for cls in CLASSES:
        assert cls in detail
    assert "DECISION 1a" in detail
    assert "neither its text nor its presence" in detail
    # The log must carry the group whose eligibility rests on the flag alone. It is the only thing
    # the amendment leaves saying so at this stage [§3] amendment, 2026-08-13].
    assert "carry no documented reason" in detail
    assert config.TREATMENT_LABELS[0] in detail


def test_it_renders_under_the_existing_derivations_heading():
    _, audit = _pipeline()
    rendered = audit.to_markdown()
    assert f"## {data._HEADINGS['derivation']}" in rendered
    assert f"- **{config.ELIGIBILITY}** (n=2)" in rendered


def test_this_stage_adds_no_kind_and_stage_5_added_exactly_cohort():
    # Asserted directly, so a future implementer cannot quietly add a kind. Stage 4 added none: its
    # entry is a `derivation`, because it describes a column. Stage 5 added exactly one — `cohort`,
    # for rows *removed*, which no other kind is about — and left `_MUST_NAME_CASES` alone, because a
    # kind-keyed rule cannot say "the two entries that remove patients must name them, the flow table
    # that removes nobody need not"; `cohort._record_removal` states it exactly instead [Stage 5 §6].
    assert data.KINDS == (
        "provenance", "contract", "correction", "observation",
        "derivation", "cohort", "model", "structural", "missingness")
    assert data._MUST_NAME_CASES == frozenset({"correction", "observation"})


def test_two_runs_render_identical_markdown():
    assert _pipeline()[1].to_markdown() == _pipeline()[1].to_markdown()


def test_the_log_is_identical_across_interpreters_with_different_hash_seeds():
    """The driver is written out rather than sketched, for Stage 3 §12.12's reason: an implementer
    choosing it freely can choose one that renders no classification at all and still see two
    identical outputs.

    FIXTURE and not WORKBOOK, so the test runs on a checkout with no `data/`.
    """
    script = ("import sys, data, derive, eligibility\n"
              "df, audit = data.load(data.FIXTURE)\n"
              "df = derive.derive(df, audit)\n"
              "df = eligibility.classify(df, audit)\n"
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
    # and the driver really did render the classification, so the comparison is of something
    assert "## Derivations" in rendered[0]
    assert config.ELIGIBILITY in rendered[0]
    for centre in config.CENTER_ORDER:
        assert centre in rendered[0]


# --- 12.10  no bare assert in eligibility.py ------------------------------------------------------

def _assert_lines(source: str) -> list[int]:
    return [node.lineno for node in ast.walk(ast.parse(source)) if isinstance(node, ast.Assert)]


def test_eligibility_py_contains_no_bare_assert():
    # A check written as `assert` disappears under -O. Every check in the module raises SchemaError.
    assert _assert_lines(MODULE.read_text(encoding="utf-8")) == []


def test_the_assert_scan_actually_fires():
    assert _assert_lines("def f():\n    assert True\n") == [2]


def test_neither_file_is_exempt_from_the_raw_name_scan():
    # test_config.py's 9.4 sweep covers both modules the moment they exist — asserted here too, so
    # the guarantee is visible from the stage that has to keep it.
    assert "eligibility.py" not in config.EXEMPT_FROM_RAW_NAME_SCAN
    assert "test_eligibility.py" not in config.EXEMPT_FROM_RAW_NAME_SCAN


# --- 12.11  Stage 4 does not depend on Stage 3 ----------------------------------------------------

def test_classify_runs_on_a_stage_2_frame_and_produces_the_identical_column():
    """What Stages 12 and 13 rely on when they classify the unrestricted frame — cheap to assert and
    expensive to rediscover."""
    stage_2, audit = _stage2(hand_frame())
    from_stage_2 = eligibility.classify(stage_2, data.Audit(_SOURCE))
    from_derived = eligibility.classify(derive.derive(stage_2, audit), audit)

    assert len(stage_2.columns) == 25 and len(from_stage_2.columns) == 26
    assert len(from_derived.columns) == 32
    pd.testing.assert_series_equal(
        from_stage_2[config.ELIGIBILITY], from_derived[config.ELIGIBILITY])


# --- 12.12  the Stage 1 amendments ----------------------------------------------------------------
#
# In test_config.py's 9.3c rather than here: the constants are Stage 1's, and their acceptance
# criteria belong with the file that owns them.


# --- 12.13  structural facts from the workbook ----------------------------------------------------
#
# Data-gated, against v7 and the spec's §18 record. This file declares its own `workbook` fixture:
# test_data.py's is module-scoped, does not cross files, and produces an unclassified frame in any
# case. Two of the tests below moved here from test_data.py §12.14 — both assert Stage 4's
# preconditions, and leaving them there would mean the property was tested against one workbook and
# not against the classifier.

# Under DECISION 1a. The `indeterminate` column is gone and its counts have moved into `eligible`;
# every `n` and the [§3] restriction-1 row are unchanged, which is the table-level evidence that the
# amendment relabelled records rather than moving them.
_V7_CROSSTAB: tuple[tuple[str, ...], ...] = (
    ("HUG",    config.TREATMENT_LABELS[0], "11",  "11", "22"),
    ("HUG",    config.TREATMENT_LABELS[1], "30",  "0",  "30"),
    ("CHUV",   config.TREATMENT_LABELS[0], "14",  "0",  "14"),
    ("CHUV",   config.TREATMENT_LABELS[1], "7",   "0",  "7"),
    ("Lugano", config.TREATMENT_LABELS[0], "29",  "0",  "29"),
    ("Lugano", config.TREATMENT_LABELS[1], "2",   "0",  "2"),
    ("USZ",    config.TREATMENT_LABELS[0], "14",  "8",  "22"),
    ("USZ",    config.TREATMENT_LABELS[1], "0",   "0",  "0"),   # ← [§3] restriction 1
    ("all",    "both",                     "107", "19", "126"))


@pytest.fixture(scope="module")
def workbook() -> tuple[pd.DataFrame, data.Audit]:
    df, audit = data.load()
    return eligibility.classify(derive.derive(df, audit), audit), audit


@DATA_GATED
def test_no_treated_patient_is_flagged_as_contraindicated(workbook):
    # E3's premise, and the assertion this stage exists to add. Moved here from test_data.py §12.14:
    # here it is a property of the classifier as well as a data-gated fact about one workbook.
    df, _ = workbook
    flagged = df[(df[config.TREATMENT] == 1) & (df["ivt_contraindicated"] == 1)]
    assert len(flagged) == 0


@DATA_GATED
def test_the_eligibility_classifier_is_never_missing(workbook):
    # E1's and E2's premises. Moved here from test_data.py §12.14 for the same reason.
    df, _ = workbook
    assert int(df["ivt_contraindicated"].isna().sum()) == 0
    assert set(df["ivt_contraindicated"]) == {0, 1}


@DATA_GATED
def test_every_centre_is_declared(workbook):
    # E4's premise, and the reason the cells of the table below total 126.
    df, _ = workbook
    assert set(df["center"]) <= set(config.CENTER_ORDER)


@DATA_GATED
def test_the_exposure_is_zero_one_and_never_missing(workbook):
    # E5's premise. The classifier reads it twice and both readers treat absence as "not treated".
    df, _ = workbook
    assert int(df[config.TREATMENT].isna().sum()) == 0
    assert set(df[config.TREATMENT]) == {0, 1}


@DATA_GATED
def test_the_classification_is_107_eligible_and_19_ineligible(workbook):
    df, _ = workbook
    assert df[config.ELIGIBILITY].value_counts().to_dict() == {
        config.ELIGIBLE: 107, config.INELIGIBLE: 19}
    # 107 = DECISION 1's 64 eligible + 43 indeterminate. The retained set is what the amendment
    # leaves untouched, and it is the number Stage 5 restricts from.
    assert int(eligibility.retained(df).sum()) == 107


@DATA_GATED
def test_the_cross_tab_reproduces_the_nine_rows_including_the_all_zero_bridging_row(workbook):
    _, audit = workbook
    table = audit.entry("derivation", config.ELIGIBILITY).table
    assert table[0] == ("centre", "arm", *CLASSES, "n")
    assert table[1:] == _V7_CROSSTAB


@DATA_GATED
def test_all_nineteen_ineligible_records_carry_a_reason_and_no_treated_patient_does(workbook):
    df, _ = workbook
    ineligible = df[config.ELIGIBILITY] == config.INELIGIBLE
    assert int(ineligible.sum()) == 19
    assert int(df.loc[ineligible, "contraindication_reason"].notna().sum()) == 19
    treated = df[config.TREATMENT] == 1
    assert int(df.loc[treated, "contraindication_reason"].notna().sum()) == 0


@DATA_GATED
def test_the_deleted_revealed_fact_mask_would_change_nothing(workbook):
    """DECISION 1a deletes the revealed-fact mask, and this is the test that it was safe to delete.

    Under DECISION 1 that one line moved 39 of 126 records. Under DECISION 1a `flag == 0` already
    classifies all 39 treated patients as eligible, so putting the mask back changes nothing — and
    **E3 is what keeps that true**. A workbook in which one treated patient carried the flag would
    make the deletion consequential, and E3 raises on it rather than letting this test discover it.
    """
    df, _ = workbook
    with_mask = _chain_with_the_deleted_revealed_fact_mask(df)
    pd.testing.assert_series_equal(
        df[config.ELIGIBILITY], with_mask, check_names=False, check_dtype=False)
    assert int((df[config.TREATMENT] == 1).sum()) == 39
    assert set(df.loc[df[config.TREATMENT] == 1, config.ELIGIBILITY]) == {config.ELIGIBLE}


@DATA_GATED
def test_43_retained_controls_carry_no_documented_reason(workbook):
    """[§3]'s amendment of 2026-08-13 makes this group indistinguishable in the frame, so it is
    counted here and in Stage 5's cohort-flow table or nowhere.

    It is the same 43 records DECISION 1 labelled `indeterminate`: every one a control, at the two
    centres that never collected a reason. Their eligibility now rests on `ivt_contraindicated = 0`
    alone, and nothing in the data can test whether that 0 was an assessment or an unfilled default.
    """
    df, _ = workbook
    undocumented = df["contraindication_reason"].isna() & (df[config.TREATMENT] == 0)
    assert int(undocumented.sum()) == 43
    assert set(df.loc[undocumented, "center"]) == {"CHUV", "Lugano"}
    assert set(df.loc[undocumented, config.ELIGIBILITY]) == {config.ELIGIBLE}
    # Over all arms it is 82, because no treated patient carries a reason either — which is why the
    # count that matters is arm-restricted [Stage 5 §7.3].
    assert int(df["contraindication_reason"].isna().sum()) == 82


@DATA_GATED
def test_classify_takes_derives_frame_from_31_to_32_columns(workbook):
    df, _ = workbook
    assert len(df.columns) == 32


# --- 12.14  Stage 4 adds one column and changes nothing else --------------------------------------

def _promise_frames() -> list[tuple[str, pd.DataFrame, data.Audit]]:
    hand, hand_audit = _stage2(hand_frame())
    fixture, fixture_audit = data.load(data.FIXTURE)
    return [("hand through Stage 2", hand, hand_audit),
            ("fixture through derive", derive.derive(fixture, fixture_audit), fixture_audit)]


@pytest.mark.parametrize("label, df, audit", _promise_frames())
def test_one_column_is_appended_and_nothing_else_moves(label, df, audit):
    """§0.2's promise, asserted rather than rested on the `copy()` call.

    The last assertion is the one that fails if `df.copy()` is ever moved below the assignment or
    turned into a slice, and nothing else in §12 would notice. The third is the one that fails if a
    future edit "tidies" a dtype on the way past.
    """
    before = df.copy()
    out = eligibility.classify(df, audit)

    assert list(out.columns) == [*df.columns, config.ELIGIBILITY], label   # appended, and last
    assert len(out) == len(df), label                                     # no row dropped
    assert out.drop(columns=[config.ELIGIBILITY]).equals(before), label   # no value edited
    assert df.equals(before), label                                       # caller's frame untouched
