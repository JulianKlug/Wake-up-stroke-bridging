"""Acceptance tests for Stage 3 — §12 of `specs/stage3_derived_variables.md`.

Tests that need the private workbook are marked `skipif(not DATA_XLSX.exists())`; every other test
runs on any checkout, against `tests/fixture_schema.xlsx` or a hand-built frame.

This file is **not** exempt from the Stage 1 §7 raw-name scan and must not become exempt. Every frame
below is built with analysis names.

The hand-built frame is `test_data.py`'s, imported rather than re-declared: two hand frames drift,
and this one is already contract-valid over six records, so any error a test sees comes from the
corruption that test applied. It needs no extending — every column Stage 3 reads is already in it —
and it already encodes the missing-volume, tie and both-arm cases. What it needs is per-record
*variation*, and `hand_frame(**overrides)` cannot give it: that assigns a scalar to a whole column.
`corrupt(column, value, where="HAND-N")` is the per-record tool, and it matches on `case_id`, never
on row position.

The frame's derived expectations fall out of the records as they stand, and are not hand-tuned::

    core_ml over the six records         (13, 20, None, 5, 0, 30)
      median (skipna, over 0, 5, 13, 20, 30)   13.0    ← equals HAND-1's own value
      above  (> 13.0)                          2       HAND-2 (20), HAND-6 (30)
      at     (== 13.0)                         1       HAND-1  → core_above_median = 0  [ties below]
      below  (< 13.0)                          2       HAND-5 (0), HAND-4 (5)
      missing                                  1       HAND-3  → core_above_median = <NA>

    mrs_90d = 2 on every record          the four dichotomies are NOT uniformly events:
      mrs_0_2_90d   2 <= 2   →  1            one event column and three non-event columns, which
      mrs_0_1_90d   2 <= 1   →  0            is what makes the frame usable for 12.5 at all — a
      death_90d     2 == 6   →  0            frame where every dichotomy read alike could not
      mrs_5_6_90d   2 >= 5   →  0            distinguish a registry-driven loop from a constant
"""
from __future__ import annotations

import ast
import dataclasses
import os
import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest

import config
import data
import derive
from test_data import corrupt, hand_frame, hand_source, run

DATA_GATED = pytest.mark.skipif(
    not config.DATA_XLSX.exists(),
    reason="the private workbook is gitignored and absent from this checkout")

MODULE = Path(derive.__file__).resolve()
MODULE_DIR = MODULE.parent

ONSET_LEVELS = config.FACTOR_LEVELS["onset_type"]
BASELINE = config.REFERENCE_LEVELS["onset_type"]


# --- helpers ----------------------------------------------------------------------------------
#
# `_stage2` puts a hand-built frame through the pipeline Stage 3 actually receives its input from,
# which matters for one thing in particular: `center` arrives as a raw code and Stage 2 recodes it,
# and the derived missingness table counts absences per CENTER_ORDER label. A frame handed straight
# to `derive` would render four columns of zeros and 12.9's per-centre reconciliation would be
# asserting nothing.

_SOURCE = hand_source(Path("."), "hand_built_stage3")


def _stage2(df: pd.DataFrame) -> tuple[pd.DataFrame, data.Audit]:
    return run(df, _SOURCE)


def derived(df: pd.DataFrame | None = None) -> tuple[pd.DataFrame, data.Audit]:
    """The hand frame (or `df`) through Stage 2 and then through `derive`."""
    frame, audit = _stage2(hand_frame() if df is None else df)
    return derive.derive(frame, audit), audit


def where(df: pd.DataFrame, case_id: str, column: str):
    """One cell, addressed by identifier — never by row position."""
    return df.loc[df["case_id"] == case_id, column].iloc[0]


def message_of(excinfo) -> str:
    return str(excinfo.value)


# --- 12.1  the onset assertion fires ---------------------------------------------------------------

def test_a_record_positive_on_both_onset_flags_raises_naming_it():
    df = corrupt("wake_up", 1, where="HAND-2")
    df.loc[df["case_id"] == "HAND-2", "unwitnessed"] = 1
    with pytest.raises(config.SchemaError) as e:
        derive.derive(*_stage2(df))
    message = message_of(e)
    assert "HAND-2" in message
    assert "wake_up" in message and "unwitnessed" in message
    # It refuses to choose rather than resolving: a resolution would apply silently to every future
    # workbook, and onset_type is a [§6] covariate.
    assert "do not choose a level here" in message


def test_each_flag_positive_on_a_different_record_does_not_raise():
    df = corrupt("wake_up", 1, where="HAND-2")
    df.loc[df["case_id"] == "HAND-3", "unwitnessed"] = 1
    out, _ = derived(df)
    assert where(out, "HAND-2", "onset_type") == "wake_up"
    assert where(out, "HAND-3", "onset_type") == "unwitnessed"


def test_one_missing_flag_is_not_evidence_that_both_are_positive():
    # A missing flag is ordinary missingness [§13], handled by _onset_type, not by the assertion.
    df = corrupt("wake_up", pd.NA, where="HAND-4")
    df.loc[df["case_id"] == "HAND-4", "unwitnessed"] = 1
    assert derive._assert_onset_flags(derived(df)[0]) is None


def test_the_hand_frame_itself_trips_no_assertion():
    assert len(derived()[0]) == 6


# --- 12.2  the three levels partition the frame ------------------------------------------------------

def three_level_frame() -> pd.DataFrame:
    """One record per level. Built with `corrupt`, not with overrides, which assign a whole column."""
    df = corrupt("wake_up", 1, where="HAND-2")
    df.loc[df["case_id"] == "HAND-3", "unwitnessed"] = 1
    return df


def test_every_record_carries_exactly_one_declared_level():
    out, _ = derived(three_level_frame())
    onset = out["onset_type"]
    assert set(onset.dropna()) <= set(ONSET_LEVELS)
    assert set(onset.dropna()) == set(ONSET_LEVELS), "the frame must cover all three levels"
    counts = [int((onset == level).sum()) for level in ONSET_LEVELS]
    assert sum(counts) == len(out) - int(onset.isna().sum())


def test_no_fourth_value_appears():
    out, _ = derived(three_level_frame())
    assert set(out["onset_type"].dropna()) - set(ONSET_LEVELS) == set()


def test_onset_type_is_a_string_and_not_a_categorical():
    # Stage 5 restricts the cohort; a categorical that keeps a dead level makes every later
    # groupby(observed=False) resurrect it as an all-missing row. Stage 6 builds the categorical
    # from FACTOR_LEVELS, at the point of use.
    out, _ = derived(three_level_frame())
    assert not isinstance(out["onset_type"].dtype, pd.CategoricalDtype)
    assert out["onset_type"].dtype == "string"


def test_the_baseline_is_the_declared_reference_level():
    out, _ = derived()
    # _HAND_CONSTANT has both flags 0 on every record, so the whole frame sits at the baseline.
    assert set(out["onset_type"]) == {BASELINE}


# --- 12.3  a missing flag yields <NA>, never the baseline and never a level ----------------------------
#
# The §4.3 branch that no data-gated test can reach: Stage 0 records both flags as non-missing on all
# 126 records, so an implementation that omitted the mask would be green on v7 and wrong on v8.
#
# The two guards are pinned separately, because each fails differently and the combined assertion
# sees neither. The variants below are reimplemented here, locally and in two lines — derive.py is
# never monkey-patched — so the shipped function keeps both guards while the test still documents
# what each one buys.


def _missing_flag_frame() -> pd.DataFrame:
    """HAND-4 has `wake_up` <NA> and `unwitnessed` 0. Nothing upstream forbids it: Stage 2's A5
    asserts BINARY_COLUMNS take values in {0, 1} *or missing*, and both flags are in that list."""
    return corrupt("wake_up", pd.NA, where="HAND-4")


def _onset_without_the_unknown_mask(df: pd.DataFrame) -> pd.Series:
    """§4.2's square, cell 2: `.fillna(False)` kept, the `unknown` mask dropped."""
    onset = pd.Series(BASELINE, index=df.index, dtype="string")
    for flag, level in config.ONSET_TYPE_FROM_FLAG:
        onset = onset.mask((df[flag] == 1).fillna(False), level)
    return onset


def _onset_without_the_fillna(df: pd.DataFrame) -> pd.Series:
    """Cell 3: the `unknown` mask kept, `.fillna(False)` dropped. Correct only by accident."""
    onset = pd.Series(BASELINE, index=df.index, dtype="string")
    unknown = pd.Series(False, index=df.index)
    for flag, level in config.ONSET_TYPE_FROM_FLAG:
        onset = onset.mask(df[flag] == 1, level)
        unknown = unknown | df[flag].isna()
    return onset.mask(unknown)


def _onset_without_either_guard(df: pd.DataFrame) -> pd.Series:
    """Cell 4: the plausible "tidy the masks" edit."""
    onset = pd.Series(BASELINE, index=df.index, dtype="string")
    for flag, level in config.ONSET_TYPE_FROM_FLAG:
        onset = onset.mask(df[flag] == 1, level)
    return onset


def test_a_missing_flag_yields_na_for_onset_type():
    out, _ = derived(_missing_flag_frame())
    assert pd.isna(where(out, "HAND-4", "onset_type"))
    assert where(out, "HAND-4", "onset_type") is not BASELINE


def test_a_missing_flag_yields_na_for_unknown_onset_too():
    # "We do not know the onset type" is not "the onset was unknown", and never 0 either.
    out, _ = derived(_missing_flag_frame())
    assert pd.isna(where(out, "HAND-4", "unknown_onset"))


def test_only_the_record_with_the_missing_flag_is_affected():
    out, _ = derived(_missing_flag_frame())
    assert int(out["onset_type"].isna().sum()) == 1
    assert set(out["onset_type"].dropna()) == {BASELINE}


def test_dropping_the_unknown_mask_fabricates_the_baseline():
    # Cell 2. The fabricated covariate value §4.3 warns about: the record would enter the propensity
    # model, the outcome regression and the balance table labelled `witnessed`.
    frame, _ = _stage2(_missing_flag_frame())
    variant = _onset_without_the_unknown_mask(frame)
    assert variant[frame["case_id"] == "HAND-4"].iloc[0] == BASELINE
    assert pd.isna(derive._onset_type(frame)[frame["case_id"] == "HAND-4"].iloc[0])


def test_dropping_the_fillna_alone_is_invisible():
    # Cell 3 — correct, but by accident: the <NA> condition masks the row to the LEVEL, and the
    # unknown mask then overwrites it. Pinned so that the square is complete and cell 4 is not read
    # as evidence that `.fillna(False)` is doing this work on its own.
    frame, _ = _stage2(_missing_flag_frame())
    pd.testing.assert_series_equal(
        _onset_without_the_fillna(frame), derive._onset_type(frame), check_names=False)


def test_dropping_both_guards_fabricates_a_non_baseline_level():
    # Cell 4, and the worst of the four: a fabricated `wake_up` reads more plausibly in a table than
    # a fabricated `witnessed` does. This is the failure the "a mask carrying <NA> is not a mask"
    # rationale hid — <NA> in a condition IS treated as True, and takes the replacement.
    frame, _ = _stage2(_missing_flag_frame())
    variant = _onset_without_either_guard(frame)
    fabricated = variant[frame["case_id"] == "HAND-4"].iloc[0]
    assert fabricated == "wake_up"
    assert fabricated in ONSET_LEVELS and fabricated != BASELINE


# --- 12.4  the levels come from the configuration ------------------------------------------------------

def _level_literals_in_source(source: str) -> list[tuple[int, str]]:
    """Every string literal in `source` that is exactly a declared onset level, with its line number.

    The same technique as Stage 1's 9.4 raw-name scan and Stage 2's 12.10 assert scan. A level typed
    here would give it a second declaration, and the two would be free to drift on the one edit — a
    renamed level — that nothing else in the repository would catch.
    """
    levels = set(ONSET_LEVELS)
    return [(node.lineno, node.value)
            for node in ast.walk(ast.parse(source))
            if isinstance(node, ast.Constant) and isinstance(node.value, str)
            and node.value in levels]


def test_derive_py_writes_no_onset_level():
    offenders = _level_literals_in_source(MODULE.read_text(encoding="utf-8"))
    assert not offenders, f"derive.py names an onset level at {offenders}"


def test_the_level_scan_actually_fires():
    # A scan that silently matches nothing would otherwise pass as a green test.
    snippet = f'onset = np.select(conditions, choices, {BASELINE!r})\n'
    assert _level_literals_in_source(snippet) == [(1, BASELINE)]


def test_derive_py_names_no_raw_workbook_header():
    # test_config.py's 9.4 sweep covers this module the moment it exists — asserted here too, so the
    # guarantee is visible from the stage that has to keep it.
    assert "derive.py" not in config.EXEMPT_FROM_RAW_NAME_SCAN
    assert "test_derive.py" not in config.EXEMPT_FROM_RAW_NAME_SCAN


# --- 12.5  the dichotomies -------------------------------------------------------------------------

_READ_DIRECTLY = tuple(k for k, o in config.OUTCOMES.items() if o.source is None and k != "mrs_90d")


def test_exactly_one_column_per_registry_entry_with_a_source():
    before, _ = _stage2(hand_frame())
    existing = set(before.columns)
    after, _ = derived()
    new = set(after.columns) - existing
    assert new == {"onset_type", "unknown_onset", *config.DERIVED_DICHOTOMIES}
    assert set(config.DERIVED_DICHOTOMIES) == {
        k for k, o in config.OUTCOMES.items() if o.source is not None}


def test_the_outcomes_read_directly_are_left_alone():
    # `Outcome.source is None` is the only test of which is which, so this is what the `continue`
    # buys. The three are already columns of the frame; a loop that rebuilt them would need a source
    # that does not exist.
    frame, _ = _stage2(hand_frame())
    before = {key: frame[key].tolist() for key in _READ_DIRECTLY}
    after = derive._dichotomies(frame)              # writes in place; this test owns the frame
    for key in _READ_DIRECTLY:
        assert after[key].tolist() == before[key], key


@pytest.mark.parametrize("key", config.DERIVED_DICHOTOMIES)
def test_each_dichotomy_takes_values_in_zero_one_or_missing(key):
    out, _ = derived()
    assert str(out[key].dtype) == "Int64"           # Int64, not boolean: the consumers do arithmetic
    assert set(out[key].dropna()) <= {0, 1}


@pytest.mark.parametrize("key", config.DERIVED_DICHOTOMIES)
def test_each_dichotomy_equals_its_registry_rule_wherever_the_source_is_present(key):
    frame, _ = _stage2(corrupt("mrs_90d", pd.NA, where="HAND-3"))
    outcome = config.OUTCOMES[key]
    src = frame[outcome.source]
    expected = config.OPS[outcome.op](src, outcome.threshold)
    out = derive._dichotomies(frame)                # this test owns the frame it passes
    present = src.notna()
    assert (out.loc[present, key] == expected[present].astype("Int64")).all()


def test_the_derived_dichotomies_are_not_added_to_the_binary_domain_check():
    # BINARY_COLUMNS is Stage 2's A5 domain check over the frame *as read*. Adding these names to it
    # would have Stage 2 assert columns that do not exist when it runs.
    assert set(config.DERIVED_DICHOTOMIES) & set(config.BINARY_COLUMNS) == set()


# --- 12.6  registry-driven, not literal ---------------------------------------------------------------

def test_a_patched_threshold_moves_the_produced_column(monkeypatch):
    """Run the loop under a registry whose mRS 0-2 threshold is 3. If the column does not move, a
    threshold is written in derive.py — and a threshold typed twice can have the manuscript print
    `mRS 0-1` while the analysis computes `mRS 0-2`, with no test failing."""
    key = "mrs_0_2_90d"
    patched = dict(config.OUTCOMES)
    patched[key] = dataclasses.replace(config.OUTCOMES[key], threshold=3)
    monkeypatch.setattr(config, "OUTCOMES", patched)

    # mrs_90d = 3 on every record: an event under the patched threshold, a non-event under the real
    # one, so the two registries cannot agree by accident.
    frame, _ = _stage2(hand_frame(mrs_90d=3))
    assert set(derive._dichotomies(frame)[key]) == {1}

    monkeypatch.undo()
    frame, _ = _stage2(hand_frame(mrs_90d=3))
    assert set(derive._dichotomies(frame)[key]) == {0}


def test_the_printed_rule_is_computed_from_the_same_fields_the_loop_applies():
    _, audit = derived()
    table = audit.entry("derivation", "dichotomies").table
    printed = {row[0]: row[1] for row in table[1:]}
    assert printed == {key: config.OUTCOMES[key].rule for key in config.DERIVED_DICHOTOMIES}


# --- 12.7  missingness is reimposed  [roadmap Stage 3's acceptance criterion, invariant 6] --------------
#
# Parametrised over every derived dichotomy AND over both source dtypes, and the second dtype is what
# makes this criterion testable at all.
#
# `mrs_90d` is declared Int64, so the comparison already propagates <NA> and `.mask(src.isna())` is a
# no-op on every frame this repository builds. The Int64 half of this parametrisation therefore
# PASSES WITH THE MASK DELETED. The float64 half is the one that fails — and float64 is the shape any
# frame built without READ_DTYPES has: a read whose `dtype=` argument was dropped, a frame assembled
# in a notebook, a fixture built from a dict and not cast.
#
# Definition of done 5 is watched against the float64 half. A future reader who runs only the Int64
# case will otherwise conclude the mask is dead code and delete it.

_SOURCE_DTYPES = ["Int64", "float64"]


def _frame_with_a_missing_mrs(dtype: str) -> pd.DataFrame:
    frame, _ = _stage2(corrupt("mrs_90d", pd.NA, where="HAND-3"))
    frame["mrs_90d"] = frame["mrs_90d"].astype(dtype)
    return frame


@pytest.mark.parametrize("dtype", _SOURCE_DTYPES)
@pytest.mark.parametrize("key", config.DERIVED_DICHOTOMIES)
def test_each_dichotomy_carries_exactly_its_sources_missingness(key, dtype):
    frame = _frame_with_a_missing_mrs(dtype)
    expected = frame["mrs_90d"].isna()
    out = derive._dichotomies(frame)                # this test owns the frame it passes
    pd.testing.assert_series_equal(out[key].isna(), expected, check_names=False)


@pytest.mark.parametrize("dtype", _SOURCE_DTYPES)
def test_a_record_with_no_mrs_is_missing_in_all_four_and_is_zero_in_none(dtype):
    """The trap, in its own test. Without the mask on a float64 source the record becomes a non-event
    in every one of the four — not good-outcome, not poor-outcome, and alive — and two things break
    at once that [§11] names: the denominator moves, and a patient whose outcome is unknown is
    counted as having had a specific one. The column is full, the counts are plausible, and every
    downstream table reconciles, which is why it is asserted rather than trusted."""
    frame = _frame_with_a_missing_mrs(dtype)
    out = derive._dichotomies(frame)
    row = out["case_id"] == "HAND-3"
    for key in config.DERIVED_DICHOTOMIES:
        value = out.loc[row, key].iloc[0]
        assert pd.isna(value), f"{key} is {value!r} on a record with no 90-day mRS"
        # and it is not a non-event either. `<NA> == 0` is <NA>, which sums to 0 in Int64, so this
        # counts the record only if the mask went and the comparison manufactured a zero.
        assert int((out.loc[row, key] == 0).sum()) == 0, key


def test_the_float64_half_is_the_one_that_can_fail():
    """Pins the asymmetry itself, so that the parametrisation above is not read as two equal halves.

    On the declared Int64 source the comparison propagates <NA> and the mask changes nothing; on a
    float64 source it absorbs, and the mask is the only thing standing between a missing outcome and
    a fabricated non-event.
    """
    frame = _frame_with_a_missing_mrs("float64")
    o = config.OUTCOMES["mrs_0_2_90d"]
    unmasked = config.OPS[o.op](frame["mrs_90d"], o.threshold).astype("Int64")
    assert unmasked[frame["case_id"] == "HAND-3"].iloc[0] == 0        # fabricated non-event

    frame = _frame_with_a_missing_mrs("Int64")
    unmasked = config.OPS[o.op](frame["mrs_90d"], o.threshold).astype("Int64")
    assert pd.isna(unmasked[frame["case_id"] == "HAND-3"].iloc[0])    # propagates; mask redundant


# --- 12.8  unknown_onset follows onset_type ------------------------------------------------------------

def test_unknown_onset_is_the_complement_of_the_baseline_level():
    out, _ = derived(three_level_frame())
    expected = (out["onset_type"] != BASELINE).astype("Int64").mask(out["onset_type"].isna())
    pd.testing.assert_series_equal(out["unknown_onset"], expected, check_names=False)
    assert str(out["unknown_onset"].dtype) == "Int64"


def test_unknown_onset_carries_missingness_and_is_never_one_where_onset_type_is_missing():
    out, _ = derived(_missing_flag_frame())
    pd.testing.assert_series_equal(
        out["unknown_onset"].isna(), out["onset_type"].isna(), check_names=False)


def test_unknown_onset_is_derived_from_onset_type_and_not_from_the_flags():
    """Overwrite `onset_type` and confirm `unknown_onset` moves with it.

    Two derivations of one concept drift, and this pair would drift silently: both are correct on v7,
    and they differ only where a flag is missing — where the flag route yields 0 and the onset_type
    route yields <NA>.
    """
    out, _ = derived()
    assert set(out["unknown_onset"]) == {0}          # every record is at the baseline
    overwritten = pd.Series(ONSET_LEVELS[-1], index=out.index, dtype="string")
    assert set(derive._unknown_onset(overwritten)) == {1}
    # and the flags are untouched, so a flag-driven implementation would still read 0
    assert set(out["wake_up"]) == {0} and set(out["unwitnessed"]) == {0}


# --- 12.9  the derived missingness table [§11] -----------------------------------------------------------

def _derived_rows(audit: data.Audit) -> dict[str, tuple[str, ...]]:
    table = audit.entry("missingness", "absence_by_derived_column").table
    return {row[0]: row for row in table[1:]}


def test_the_derived_table_has_one_row_per_row_wise_derived_column():
    # ROW_WISE_DERIVED, not DERIVED_NAMES: core_above_median does not exist at this point, and is
    # accounted for in derive_cohort's own entry. The two counts differ by one by design.
    _, audit = derived()
    rows = _derived_rows(audit)
    assert sorted(rows) == sorted(config.ROW_WISE_DERIVED)
    assert "core_above_median" not in rows
    assert len(rows) == len(config.DERIVED_NAMES) - 1


@pytest.mark.parametrize("column", sorted(config.ROW_WISE_DERIVED))
def test_no_derived_column_is_structural_or_informative_absence(column):
    """Every absence a derived column carries is inherited — `onset_type` and `unknown_onset` from
    the onset flags, the four dichotomies from `mrs_90d`. `structural` means "the value does not
    exist for this patient", and it is the one label that tells a reader an absence is not data
    loss; a derived column given it to make a table look better would misrepresent exactly what the
    distinction was built to represent."""
    _, audit = derived(corrupt("mrs_90d", pd.NA, where="HAND-3"))
    assert _derived_rows(audit)[column][1] in {"missing", "complete"}


def test_every_derived_row_reconciles_to_the_frame():
    out, audit = derived(corrupt("mrs_90d", pd.NA, where="HAND-3"))
    for column, row in _derived_rows(audit).items():
        assert int(row[2]) + int(row[3]) == len(out), column
        per_centre = [int(cell) for cell in row[5:]]
        assert sum(per_centre) == int(row[3]), f"{column}: per-centre counts do not sum"


def test_the_derived_table_uses_the_same_header_as_stage_2s():
    # One classification, called twice: the per-centre columns, the CENTER_ORDER header and the
    # reconciliation all come along unchanged rather than being written a second time.
    _, audit = derived()
    assert audit.entry("missingness", "absence_by_derived_column").table[0] == (
        audit.entry("missingness", "absence_by_column").table[0])


def test_a_derived_column_with_a_missing_input_is_counted_as_missing():
    _, audit = derived(corrupt("mrs_90d", pd.NA, where="HAND-3"))
    rows = _derived_rows(audit)
    for key in config.DERIVED_DICHOTOMIES:
        assert rows[key][1] == "missing"
        assert int(rows[key][3]) == 1
    assert rows["onset_type"][1] == "complete"


# --- 12.10  core_above_median ------------------------------------------------------------------------

_HAND_MEDIAN = 13.0                    # over (13, 20, None, 5, 0, 30), skipna — HAND-1's own value


def test_derive_alone_does_not_create_the_median_subgroup():
    out, _ = derived()
    assert "core_above_median" not in out.columns


def test_derive_cohort_adds_it():
    out, audit = derived()
    cohort = derive.derive_cohort(out, audit)
    assert "core_above_median" in cohort.columns
    assert str(cohort["core_above_median"].dtype) == "Int64"
    assert len(cohort.columns) == len(out.columns) + 1


def test_a_second_call_raises_naming_the_column_and_the_frozen_median_rule():
    """The frozen-median guard. The natural shape of a bootstrap replicate is "resample, then re-run
    the pipeline", and under that shape a derive_cohort that recomputed would give every replicate
    its own cut-point and therefore its own subgroup — with no test failing, because the estimates
    would stay plausible and only the estimand would move."""
    out, audit = derived()
    cohort = derive.derive_cohort(out, audit)
    with pytest.raises(config.SchemaError) as e:
        derive.derive_cohort(cohort, audit)
    message = message_of(e)
    assert "core_above_median" in message
    assert "resamples this column" in message and "does not recompute" in message


def test_the_whole_split_matches_the_frames_own_numbers():
    out, audit = derived()
    cohort = derive.derive_cohort(out, audit)
    core, above = cohort["core_ml"], cohort["core_above_median"]
    assert float(core.median()) == _HAND_MEDIAN
    assert int((above == 1).sum()) == 2                         # HAND-2 (20), HAND-6 (30)
    assert int((core == _HAND_MEDIAN).sum()) == 1               # HAND-1
    assert int((core < _HAND_MEDIAN).sum()) == 2                # HAND-4 (5), HAND-5 (0)
    assert int(above.isna().sum()) == 1                         # HAND-3, no core volume


def test_ties_sit_below():
    # Not a formality: 3 records sit exactly on the primary cohort's median of 5.0 mL and 39 have a
    # core of exactly 0. `>` and `>=` differ on those 3 patients, and nothing else in the repository
    # would catch which was used.
    out, audit = derived()
    cohort = derive.derive_cohort(out, audit)
    assert where(cohort, "HAND-1", "core_ml") == _HAND_MEDIAN
    assert where(cohort, "HAND-1", "core_above_median") == 0


def test_a_missing_core_volume_yields_na():
    out, audit = derived()
    cohort = derive.derive_cohort(out, audit)
    assert pd.isna(where(cohort, "HAND-3", "core_ml"))
    assert pd.isna(where(cohort, "HAND-3", "core_above_median"))


def test_without_the_mask_the_missing_core_joins_the_below_median_group():
    """The one mask in Stage 3 that is load-bearing today, and unlike 12.7's Int64 half this fails on
    the first run: `core_ml` declares no dtype in the contract, so it arrives as float64 and
    `NaN > 13.0` is False — the comparison ABSORBS. The record would be silently assigned to the
    below-median group, in whichever arm it happens to sit, with no <NA> anywhere to notice."""
    out, _ = derived()
    assert str(out["core_ml"].dtype) == "float64", "the premise of this test"
    core = out["core_ml"]
    unmasked = (core > float(core.median())).astype("Int64")
    assert unmasked[out["case_id"] == "HAND-3"].iloc[0] == 0     # swept into below-median
    assert pd.isna(derive._core_above_median(out)[0][out["case_id"] == "HAND-3"].iloc[0])


def test_the_median_is_of_the_frame_it_is_given():
    """Which is why this lives in derive_cohort and not in derive. Over all 126 records the workbook's
    core median is 6.0 mL; over the 93-record primary cohort it is 5.0. One function computing both
    halves at Stage 3's natural position would have split the cohort on the wrong number, and both
    numbers are plausible."""
    out, _ = derived()
    subset = out[out["case_id"] != "HAND-6"]                     # drops the 30 mL record
    whole_median = derive._core_above_median(out)[1]
    subset_median = derive._core_above_median(subset)[1]
    assert whole_median == _HAND_MEDIAN
    assert subset_median == 9.0                                 # median of (0, 5, 13, 20)
    assert whole_median != subset_median
    # and the membership really moves with it: HAND-1 is at the median of the whole and above the
    # median of the subset
    assert derive._core_above_median(out)[0][out["case_id"] == "HAND-1"].iloc[0] == 0
    assert derive._core_above_median(subset)[0][subset["case_id"] == "HAND-1"].iloc[0] == 1


def test_the_audit_entry_prints_the_median_and_all_four_counts():
    out, audit = derived()
    derive.derive_cohort(out, audit)
    values = dict(audit.entry("derivation", "core_above_median").table[1:])
    assert values == {
        "median (mL)": "13", "records": "6",
        "above": "2", "at the median": "1", "below": "2", "missing": "1"}


# --- 12.11  constant_covariates -----------------------------------------------------------------------
#
# On the hand frame this returns TEN of the thirteen PS_COVARIATES_FULL entries, not (). _HAND_CONSTANT
# gives each of them a single value, and the expectation is pinned as a literal tuple rather than
# described — in PS_COVARIATES_FULL order, which is also what proves the caller's-order clause.
#
# Do not "fix" the frame to make this return (). It cannot be done through hand_frame(**overrides),
# which assigns a scalar and therefore makes a column *more* constant, and doing it through ten
# corrupt() calls would build a second hand frame by the back door. Asserting the ten is the stronger
# test anyway: it is the only assertion in §12 that would catch constant_covariates ranging over a set
# instead of over the caller's sequence.

_CONSTANT_ON_THE_HAND_FRAME = (
    "age", "sex", "prestroke_mrs", "nihss_baseline", "onset_type",
    "atrial_fib", "hypertension", "hyperlipidemia", "diabetes", "smoking")


def test_it_finds_the_ten_single_valued_covariates_in_the_callers_order():
    out, _ = derived()
    assert derive.constant_covariates(out) == _CONSTANT_ON_THE_HAND_FRAME
    # the three that are absent are the three that vary
    assert set(config.PS_COVARIATES_FULL) - set(_CONSTANT_ON_THE_HAND_FRAME) == {
        "core_ml", "tmax6_ml", "center"}


def test_the_result_is_in_the_callers_sequence_order():
    out, _ = derived()
    found = derive.constant_covariates(out)
    assert found == tuple(c for c in config.PS_COVARIATES_FULL if c in set(found))
    reversed_order = tuple(reversed(config.PS_COVARIATES_FULL))
    assert derive.constant_covariates(out, reversed_order) == tuple(reversed(found))


def test_a_one_level_factor_is_found_and_it_is_a_string_column():
    # Free, and covered by the baseline assertion itself: every record in the hand frame is at the
    # baseline, so `onset_type` is among the ten with no construction at all.
    out, _ = derived()
    assert out["onset_type"].dtype == "string"
    assert "onset_type" in derive.constant_covariates(out)


def test_a_numeric_column_made_constant_joins_the_ten():
    out, _ = derived(hand_frame(core_ml=7.0))
    found = derive.constant_covariates(out)
    assert set(found) == set(_CONSTANT_ON_THE_HAND_FRAME) | {"core_ml"}
    # and it joins them *in PS_COVARIATES_FULL order*, where core_ml actually sits — not appended to
    # the end, which is what a function ranging over a set and sorting would produce.
    assert found == tuple(c for c in config.PS_COVARIATES_FULL if c in set(found))
    assert found.index("core_ml") < found.index("atrial_fib")


def test_an_all_missing_column_joins_them_too():
    # nunique(dropna=True) is 0 there, so `<= 1` catches a constant column, an all-missing one and a
    # factor left with a single level, all in one expression.
    out, _ = derived(hand_frame(tmax6_ml=float("nan")))
    assert out["tmax6_ml"].nunique(dropna=True) == 0
    assert "tmax6_ml" in derive.constant_covariates(out)


def test_it_mutates_neither_the_frame_nor_the_covariate_list():
    out, _ = derived()
    before = out.copy()
    derive.constant_covariates(out)
    pd.testing.assert_frame_equal(out, before)
    # invariant 3 rests on these two resolving to one object; a function that filtered a list in
    # place would break it in a way that only shows up as a covariate quietly missing from one of
    # the two nuisance models.
    assert config.OUTCOME_COVARIATES is config.PS_COVARIATES
    assert config.PS_COVARIATES_FULL[:len(config.PS_COVARIATES)] == config.PS_COVARIATES


def test_it_asks_for_a_covariate_the_contract_cannot_supply_and_says_so():
    # The docstring's precondition: PS_COVARIATES_FULL contains onset_type, so calling this before
    # derive has run is a sequencing bug. The correct failure is loud — a silent skip would answer a
    # question about a different covariate set.
    frame, _ = _stage2(hand_frame())
    with pytest.raises(KeyError):
        derive.constant_covariates(frame)


# --- 12.12  the audit inventory and reproduction --------------------------------------------------------

_STAGE_3_INVENTORY: list[tuple[str, str]] = [
    ("derivation", "onset_type"),
    ("derivation", "unknown_onset"),
    ("derivation", "dichotomies"),
    ("missingness", "absence_by_derived_column"),
    ("derivation", "core_above_median"),
    ("derivation", "constant_covariates"),
]


def _full_pipeline() -> tuple[pd.DataFrame, data.Audit]:
    df, audit = data.load(data.FIXTURE)
    df = derive.derive(df, audit)
    return derive.derive_cohort(df, audit), audit


def test_the_six_entries_appear_with_the_declared_kinds_and_steps_in_order():
    _, audit = _full_pipeline()
    stage_3 = [(e.kind, e.step) for e in audit.entries][-len(_STAGE_3_INVENTORY):]
    assert stage_3 == _STAGE_3_INVENTORY


def test_derives_four_entries_precede_derive_cohorts_two():
    df, audit = data.load(data.FIXTURE)
    before = len(audit.entries)
    df = derive.derive(df, audit)
    assert [(e.kind, e.step) for e in audit.entries[before:]] == _STAGE_3_INVENTORY[:4]
    derive.derive_cohort(df, audit)
    assert [(e.kind, e.step) for e in audit.entries[before + 4:]] == _STAGE_3_INVENTORY[4:]


def test_the_rendered_log_carries_a_derivations_heading():
    _, audit = _full_pipeline()
    assert f"## {data._HEADINGS['derivation']}" in audit.to_markdown()
    assert "## Derivations" in audit.to_markdown()


def test_constant_covariates_is_recorded_even_when_nothing_is_constant():
    # Each of the fixture's 13 covariates has 2 distinct values over its 2 records. The entry is
    # still recorded, with n = 0 and no table — the shape Stage 2 gives penumbra_recomputed.
    _, audit = _full_pipeline()
    entry = audit.entry("derivation", "constant_covariates")
    assert entry is not None and entry.n == 0
    assert entry.table is None
    assert entry.case_ids == ()


def test_no_derivation_entry_names_a_patient():
    # A derivation describes a column, not a patient. The one place Stage 3 names patients is the
    # both-onset-flags assertion, which raises — such a frame never reaches a log.
    _, audit = _full_pipeline()
    assert all(e.case_ids == () for e in audit.entries if e.kind == "derivation")


def test_two_runs_render_identical_markdown():
    assert _full_pipeline()[1].to_markdown() == _full_pipeline()[1].to_markdown()


def test_the_log_is_identical_across_interpreters_with_different_hash_seeds():
    """Stage 3 adds two dict-ordered loops — OUTCOMES and SUBGROUPS — and one sorted one, and a
    same-process comparison cannot see iteration order at all.

    The driver is written out rather than sketched, for the reason Stage 2 §9.6 gives about pinning
    log contents: an implementer choosing it freely can choose one that renders no derivation at all
    and still see two identical outputs. FIXTURE and not WORKBOOK, because the workbook is gitignored
    and a data-gated reproduction test is no test on most checkouts.
    """
    script = ("import sys, data, derive\n"
              "df, audit = data.load(data.FIXTURE)\n"
              "df = derive.derive(df, audit)\n"
              "derive.derive_cohort(df, audit)\n"
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
    # and the driver really did render the derivations, so the comparison is of something
    assert "## Derivations" in rendered[0]
    for _, step in _STAGE_3_INVENTORY:
        assert step in rendered[0]


def test_the_fixture_reaches_every_entry():
    """Verified against tests/fixture_schema.xlsx as committed: 2 records, wake_up 1 and 0 so
    onset_type renders two of its three level rows and its <NA> row, mrs_90d <NA> on one record so
    the dichotomies render a missing count, core_ml 0 and 12 giving a median of 6.0 with one record
    either side, and both onset flags present so the assertion does not raise.

    Note the dtype asymmetry: the fixture's core_ml arrives as int64 because both its values are
    integers, where the workbook's is float64 — so this test exercises the int64 path and 12.10
    exercises float64. Neither replaces the other.
    """
    df, audit = _full_pipeline()
    assert len(df) == 2
    assert str(df["core_ml"].dtype) == "int64"
    assert dict(audit.entry("derivation", "core_above_median").table[1:])["median (mL)"] == "6"
    dichotomies = {row[0]: row for row in audit.entry("derivation", "dichotomies").table[1:]}
    assert all(row[4] == "1" for row in dichotomies.values())       # one missing mRS
    levels = dict(audit.entry("derivation", "onset_type").table[1:])
    assert levels[BASELINE] == "1" and levels["wake_up"] == "1"
    assert levels["<NA>"] == "0"


def test_neither_entry_point_mutates_its_caller_s_frame():
    df, audit = data.load(data.FIXTURE)
    before = df.copy()
    out = derive.derive(df, audit)
    pd.testing.assert_frame_equal(df, before)
    after_derive = out.copy()
    derive.derive_cohort(out, audit)
    pd.testing.assert_frame_equal(out, after_derive)


def test_stage_3_drops_no_row_and_no_column_and_edits_no_value():
    df, audit = data.load(data.FIXTURE)
    out = derive.derive(df, audit)
    cohort = derive.derive_cohort(out, audit)
    assert len(cohort) == len(df)
    assert set(df.columns) <= set(cohort.columns)
    for column in df.columns:
        pd.testing.assert_series_equal(cohort[column], df[column], obj=column)


# --- 12.14  no bare assert in derive.py -------------------------------------------------------------

def _assert_lines(source: str) -> list[int]:
    return [node.lineno for node in ast.walk(ast.parse(source)) if isinstance(node, ast.Assert)]


def test_derive_py_contains_no_bare_assert():
    # Python strips `assert` under -O, so under a flag nobody remembers setting a check whose whole
    # purpose is to fail loudly would succeed silently.
    offenders = _assert_lines(MODULE.read_text(encoding="utf-8"))
    assert not offenders, f"derive.py has assert statements at line(s) {offenders}"


def test_the_assert_scan_actually_fires():
    assert _assert_lines("def f(x):\n    assert x > 0\n") == [2]


# --- 12.15  structural facts from the workbook ------------------------------------------------------
#
# The facts §18 records, re-checked in code so that updating DATA_SHA256 also re-checks the premises.

@pytest.fixture(scope="module")
def workbook_derived() -> tuple[pd.DataFrame, data.Audit]:
    df, audit = data.load()
    return derive.derive(df, audit), audit


@pytest.fixture(scope="module")
def primary_cohort(workbook_derived) -> pd.DataFrame:
    """The [§3] restrictions, derived from the data rather than hardcoded — Stage 5 owns this for
    real; here it is only the population §18's median is quoted over."""
    df, _ = workbook_derived
    never_ivt = [c for c in config.CENTER_ORDER
                 if int(((df["center"] == c) & (df[config.TREATMENT] == 1)).sum()) == 0]
    assert tuple(never_ivt) == config.EXPECTED_NEVER_IVT
    return df[~df["center"].isin(never_ivt) & (df["ivt_contraindicated"] != 1)]


@DATA_GATED
def test_the_onset_levels_are_33_36_and_57(workbook_derived):
    df, _ = workbook_derived
    counts = df["onset_type"].value_counts()
    assert (int(counts["witnessed"]), int(counts["unwitnessed"]), int(counts["wake_up"])) == (
        33, 36, 57)
    assert int(counts.sum()) == config.N_RECORDS_EXPECTED == 126     # the three partition the frame


@DATA_GATED
def test_neither_onset_flag_is_ever_missing_and_no_record_is_positive_on_both(workbook_derived):
    df, _ = workbook_derived
    for flag, _level in config.ONSET_TYPE_FROM_FLAG:
        assert int(df[flag].isna().sum()) == 0
    assert int(((df["wake_up"] == 1) & (df["unwitnessed"] == 1)).sum()) == 0
    assert int(df["onset_type"].isna().sum()) == 0


@DATA_GATED
def test_unknown_onset_is_one_on_93_records(workbook_derived):
    df, _ = workbook_derived
    assert int((df["unknown_onset"] == 1).sum()) == 93


@DATA_GATED
@pytest.mark.parametrize("key, events", [
    ("mrs_0_2_90d", 53), ("mrs_0_1_90d", 31), ("death_90d", 37), ("mrs_5_6_90d", 42)])
def test_the_dichotomies_reproduce_stage_0s_full_frame_counts(key, events, workbook_derived):
    df, _ = workbook_derived
    assert int((df[key] == 1).sum()) == events
    assert int(df[key].isna().sum()) == 2            # the two records with no 90-day mRS


@DATA_GATED
def test_the_full_frame_median_is_6_ml(workbook_derived):
    df, _ = workbook_derived
    assert float(df["core_ml"].median()) == 6.0


@DATA_GATED
def test_the_primary_cohort_is_93_records_and_splits_at_5_ml(primary_cohort):
    # 5.0 over the cohort against 6.0 over all 126, which is the whole reason derive_cohort exists
    # as a second entry point.
    assert len(primary_cohort) == 93
    cohort = derive.derive_cohort(primary_cohort, data.Audit(data.WORKBOOK))
    core, above = cohort["core_ml"], cohort["core_above_median"]
    assert float(core.median()) == 5.0
    assert int((above == 1).sum()) == 44
    assert int((core == 5.0).sum()) == 3             # ties, which sit below
    assert int((core < 5.0).sum()) == 45
    assert int(above.isna().sum()) == 1


@DATA_GATED
def test_the_core_volume_is_exactly_zero_on_39_cohort_records(primary_cohort):
    assert int((primary_cohort["core_ml"] == 0).sum()) == 39
