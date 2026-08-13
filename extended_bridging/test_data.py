"""Acceptance tests for Stage 2 — §12 of `specs/stage2_load_and_clean.md`.

Tests that need the private workbook are marked `skipif(not DATA_XLSX.exists())`; every other test
runs on any checkout, against `tests/fixture_schema.xlsx` or a hand-built frame.

This file is **not** exempt from the Stage 1 §7 raw-name scan and must not become exempt — that scan
covers the two modules that handle the workbook, which is the whole of what it guards. So the frames
below are built with analysis names, `FIXTURE` carries the read path, and the two places that need a
raw header derive it from `COLUMN_CONTRACT` rather than writing one.
"""
from __future__ import annotations

import ast
import dataclasses
import os
import re
import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest

import config
import data

DATA_GATED = pytest.mark.skipif(
    not config.DATA_XLSX.exists(),
    reason="the private workbook is gitignored and absent from this checkout")

MODULE = Path(data.__file__).resolve()
MODULE_DIR = MODULE.parent

# Derived, never written out: the scan forbids a raw header as a literal here, and deriving it also
# means these tests cannot drift away from the contract.
RAW_OF = {analysis: raw for raw, analysis in config.RENAME.items()}
CENTER_CODES = tuple(config.CENTER_RECODE)              # ("1", "Lausanne", "Lugano", "USZ")
ANALYSIS_DTYPES = {c.name: c.dtype for c in config.COLUMN_CONTRACT.values()
                   if c.name is not None and c.dtype is not None}


# --- hand-built frames --------------------------------------------------------------------------
#
# Six records, exactly as §12.5 specifies: one carrying the copy-paste penumbra of §7.2, one
# carrying onset_to_groin_min = 999, one missing core_ml and tmax6_ml together, one control and one
# treated so A7 has both sides, and one ordinary record. The frame is contract-valid — it passes
# A1-A9 — so any SchemaError a test sees comes from the corruption that test applied.

_HAND_RECORDS: list[dict[str, object]] = [
    # the copy-paste: the core volume sits in the penumbra cell, so the stored 13 is 83 mL short of
    # the 109 - 13 = 96 that [§6]'s definition gives. The gap is deliberately large — the correction
    # has to be visible in the balance table, not lost inside _TOL
    dict(case_id="HAND-1", center=CENTER_CODES[0], ivt=1, onset_to_ivt_min=180.0,
         onset_to_groin_min=250.0, core_ml=13.0, tmax6_ml=109.0, penumbra_ml=13.0,
         contraindication_reason=None),
    # the 999 groin time, with an onset-to-IVT that makes the implied interval unremarkable
    dict(case_id="HAND-2", center=CENTER_CODES[1], ivt=1, onset_to_ivt_min=955.0,
         onset_to_groin_min=999.0, core_ml=20.0, tmax6_ml=60.0, penumbra_ml=40.0,
         contraindication_reason=None),
    # both volumes missing on the same record, which is the v7 pattern
    dict(case_id="HAND-3", center=CENTER_CODES[2], ivt=0, onset_to_ivt_min=None,
         onset_to_groin_min=400.0, core_ml=None, tmax6_ml=None, penumbra_ml=None,
         contraindication_reason="Anticoagulation"),
    dict(case_id="HAND-4", center=CENTER_CODES[3], ivt=0, onset_to_ivt_min=None,
         onset_to_groin_min=310.0, core_ml=5.0, tmax6_ml=45.0, penumbra_ml=40.0,
         contraindication_reason="Clinician decision"),
    dict(case_id="HAND-5", center=CENTER_CODES[0], ivt=1, onset_to_ivt_min=120.0,
         onset_to_groin_min=200.0, core_ml=0.0, tmax6_ml=88.0, penumbra_ml=88.0,
         contraindication_reason=None),
    dict(case_id="HAND-6", center=CENTER_CODES[2], ivt=0, onset_to_ivt_min=None,
         onset_to_groin_min=500.0, core_ml=30.0, tmax6_ml=95.0, penumbra_ml=65.0,
         contraindication_reason=None),
]

_HAND_CONSTANT: dict[str, object] = dict(
    age=71, sex=0, hypertension=1, hyperlipidemia=0, diabetes=0, smoking=1, atrial_fib=0,
    prestroke_mrs=0, wake_up=0, unwitnessed=0, nihss_baseline=14, mrs_90d=2,
    sich=0, ph2=0, tici_2b_3=1, ivt_contraindicated=0)


def hand_frame(**overrides: object) -> pd.DataFrame:
    """The six-record frame, typed as the reader would deliver it. `overrides` set a whole column."""
    df = pd.DataFrame([{**_HAND_CONSTANT, **record} for record in _HAND_RECORDS])
    for column, dtype in ANALYSIS_DTYPES.items():
        df[column] = df[column].astype(dtype)
    for column, value in overrides.items():
        df[column] = value
    return df[sorted(config.ANALYSIS_NAMES)]


def corrupt(column: str, value: object, where: str = "HAND-1") -> pd.DataFrame:
    """The hand frame with one cell replaced, matched on `case_id` — never by row position."""
    df = hand_frame()
    df.loc[df["case_id"] == where, column] = value
    return df


def hand_source(tmp_path: Path, label: str = "hand_built") -> data.Source:
    """A test-local Source. Legal by §3: SOURCES pins what data.py declares, not what Source permits.

    The path need not exist — nothing after the read touches the filesystem, and the log header
    prints `path.name` and `label`, never a path that is opened.
    """
    return data.Source(tmp_path / "hand.xlsx", None, None, label)


def run(df: pd.DataFrame, source: data.Source) -> tuple[pd.DataFrame, data.Audit]:
    audit = data.Audit(source)
    return data._pipeline_after_read(df, source, audit), audit


def message_of(excinfo) -> str:
    return str(excinfo.value)


# --- 12.1  version and provenance ----------------------------------------------------------------

def test_sources_are_exactly_the_workbook_and_the_fixture():
    # Asserted against the literal, so a third module-level source fails here rather than appearing
    # quietly, in the same way test_config.py pins OUTCOME_MODEL_OVERRIDES.
    assert data.SOURCES == (data.WORKBOOK, data.FIXTURE)
    assert len(data.SOURCES) == 2


def test_the_workbook_carries_both_guarantees_and_the_fixture_neither():
    assert data.WORKBOOK.sha256 == config.DATA_SHA256
    assert data.WORKBOOK.n_records == config.N_RECORDS_EXPECTED
    assert data.FIXTURE.sha256 is None
    assert data.FIXTURE.n_records is None


def test_source_rejects_attribute_assignment():
    with pytest.raises(dataclasses.FrozenInstanceError):
        data.FIXTURE.sha256 = "0" * 64


def test_a_wrong_hash_raises_naming_both_hashes_and_the_remedy():
    # Over the fixture's path, which exists on every checkout — so the version check is tested
    # without `data/` rather than only under DATA_GATED.
    wrong = data.Source(config.FIXTURE_XLSX, "0" * 64, None, "fixture_wrong_hash")
    with pytest.raises(config.DataVersionError) as e:
        data._verify_version(wrong)
    message = message_of(e)
    assert "0" * 8 in message                        # expected
    assert "fixture_wrong_hash" in message
    assert "DATA_SHA256" in message and "stage0_data_inventory.py" in message


def test_an_unpinned_source_skips_the_version_check():
    assert data._verify_version(data.FIXTURE) is None


@DATA_GATED
def test_the_workbook_matches_its_pinned_hash():
    assert data._verify_version(data.WORKBOOK) is None


# --- 12.2  read contract -------------------------------------------------------------------------

def test_the_fixture_reads_under_the_contract():
    raw = data._read(data.FIXTURE, data.Audit(data.FIXTURE))
    assert raw.shape[1] == config.N_COLUMNS_EXPECTED == 42
    assert raw[RAW_OF["case_id"]].dtype == "string"


def test_the_literal_na_sentinel_reads_as_missing_and_keeps_its_dtype():
    # §4.2: the workbook holds literal 'N/A' strings in the primary outcome and in tici_2b_3, so
    # NA_VALUES is load-bearing rather than a guard. The fixture carries one cell of it.
    raw = data._read(data.FIXTURE, data.Audit(data.FIXTURE))
    tici = raw[RAW_OF["tici_2b_3"]]
    assert str(tici.dtype) == "Int64"
    assert tici.isna().sum() == 1


def test_a_header_only_read_is_named_as_the_cause():
    header_only = pd.read_excel(config.FIXTURE_XLSX, sheet_name=config.SHEET, nrows=0)
    with pytest.raises(config.SchemaError) as e:
        config.assert_column_contract(header_only.columns)
    assert "header-only read" in message_of(e)


def test_a_row_count_mismatch_names_both_counts():
    wrong = data.Source(config.FIXTURE_XLSX, None, 126, "fixture_wrong_row_count")
    raw = data._read(wrong, data.Audit(wrong))
    with pytest.raises(config.SchemaError) as e:
        data._apply_contract(raw, wrong, data.Audit(wrong))
    message = message_of(e)
    assert "2" in message and "126" in message


def test_a_source_declaring_no_row_count_does_not_raise():
    audit = data.Audit(data.FIXTURE)
    raw = data._read(data.FIXTURE, audit)
    assert len(data._apply_contract(raw, data.FIXTURE, audit)) == 2


# --- 12.3  contract application ------------------------------------------------------------------

def test_rename_and_drop_leave_exactly_the_analysis_names():
    df, _ = data.load(data.FIXTURE)
    assert set(df.columns) == config.ANALYSIS_NAMES
    assert len(df.columns) == 25


def test_no_raw_header_and_no_dropped_column_survives():
    df, _ = data.load(data.FIXTURE)
    assert set(df.columns) & set(config.COLUMN_CONTRACT) == set()
    assert set(df.columns) & config.DROPPED == set()


# --- 12.4  normalisation -------------------------------------------------------------------------

def test_every_centre_is_recoded():
    df, _ = data.load(data.FIXTURE)
    assert set(df["center"]) <= set(config.CENTER_RECODE.values())


def test_centre_is_not_a_categorical():
    # Stage 5 drops the zero-bridging centre; a categorical that keeps a dead level makes every
    # later groupby(observed=False) resurrect it as an all-missing row. Stage 6 builds the
    # categorical from FACTOR_LEVELS, at the point of use.
    df, _ = data.load(data.FIXTURE)
    assert not isinstance(df["center"].dtype, pd.CategoricalDtype)
    assert df["center"].dtype == "string"


def test_an_unmapped_centre_code_raises_naming_the_value(tmp_path):
    with pytest.raises(config.SchemaError) as e:
        run(corrupt("center", "Zurich"), hand_source(tmp_path))
    message = message_of(e)
    assert "A3" in message and "center" in message and "Zurich" in message


def test_a_padded_identifier_raises_rather_than_being_repaired(tmp_path):
    with pytest.raises(config.SchemaError) as e:
        run(corrupt("case_id", " HAND-1 "), hand_source(tmp_path))
    message = message_of(e)
    assert "A2" in message and "HAND-1" in message


# --- 12.5  corrections are content-driven --------------------------------------------------------

def test_shuffling_the_input_changes_neither_the_frame_nor_the_log(tmp_path):
    """The Stage 2 counterpart of Stage 1's AST scan for raw names.

    The seam is the whole pipeline after the read, not `_correct` alone: `_normalise` and
    `_missingness` also write into the log, and a row-indexed edit in either would otherwise pass
    green. The log this produces carries five of the seven entries and not the two the read records,
    because the shuffle starts after the read — neither of those can depend on row order, one counts
    rows and the other counts columns.
    """
    source = hand_source(tmp_path)
    hand = hand_frame()
    base_frame, base_log = None, None
    for seed in (0, 1, 2, 3, 4):
        shuffled = hand.sample(frac=1, random_state=seed).reset_index(drop=True)
        out, audit = run(shuffled, source)
        out = out.sort_values("case_id").reset_index(drop=True)
        rendered = audit.to_markdown()
        if base_frame is None:
            base_frame, base_log = out, rendered
            continue
        pd.testing.assert_frame_equal(out, base_frame)
        assert rendered == base_log, f"the audit log changed under seed {seed}"


def test_the_returned_frame_keeps_a_default_range_index(tmp_path):
    out, _ = run(hand_frame(), hand_source(tmp_path))
    assert isinstance(out.index, pd.RangeIndex)
    assert list(out.index) == list(range(len(out)))


# --- 12.6  penumbra recomputation -----------------------------------------------------------------

def _penumbra_entry(audit: data.Audit) -> data.AuditEntry:
    return audit.entry("correction", "penumbra_recomputed")


def test_the_copy_pasted_penumbra_is_recomputed_and_its_case_named(tmp_path):
    out, audit = run(hand_frame(), hand_source(tmp_path))
    entry = _penumbra_entry(audit)
    assert entry.n == 1
    assert entry.case_ids == ("HAND-1",)
    corrected = out.loc[out["case_id"] == "HAND-1", "penumbra_ml"].iloc[0]
    assert corrected == pytest.approx(96.0)


def test_penumbra_is_recomputed_on_every_record_not_only_the_disagreeing_one(tmp_path):
    out, _ = run(hand_frame(), hand_source(tmp_path))
    expected = out["tmax6_ml"] - out["core_ml"]
    pd.testing.assert_series_equal(out["penumbra_ml"], expected, check_names=False)


def test_no_recomputed_penumbra_is_negative(tmp_path):
    out, _ = run(hand_frame(), hand_source(tmp_path))
    assert (out["penumbra_ml"].dropna() >= 0).all()


def test_source_missingness_is_preserved_exactly(tmp_path):
    out, _ = run(hand_frame(), hand_source(tmp_path))
    missing_volume = out["core_ml"].isna() | out["tmax6_ml"].isna()
    pd.testing.assert_series_equal(out["penumbra_ml"].isna(), missing_volume, check_names=False)
    assert int(missing_volume.sum()) == 1        # no record gained a value


def test_a_frame_that_already_agrees_records_no_identifiers(tmp_path):
    agreed = hand_frame()
    agreed["penumbra_ml"] = agreed["tmax6_ml"] - agreed["core_ml"]
    _, audit = run(agreed, hand_source(tmp_path))
    entry = _penumbra_entry(audit)
    assert entry.n == 0
    assert entry.case_ids == ()
    assert entry.table is None


def test_the_standing_query_names_its_case(tmp_path):
    # 999 is left standing, not set to missing: it is inside the observed range and the implied
    # IVT-to-groin interval is unremarkable. The pilot's correction is deliberately not carried.
    out, audit = run(hand_frame(), hand_source(tmp_path))
    entry = audit.entry("observation", "onset_to_groin_999")
    assert entry.n == 1
    assert entry.case_ids == ("HAND-2",)
    assert out.loc[out["case_id"] == "HAND-2", "onset_to_groin_min"].iloc[0] == 999


# --- 12.7  every correction names its cases -------------------------------------------------------

@pytest.mark.parametrize("kind", data.KINDS)
def test_an_entry_touching_records_without_naming_them(kind):
    must_name = kind in data._MUST_NAME_CASES
    if must_name:
        with pytest.raises(ValueError) as e:
            data.AuditEntry(kind, "step", 1, "detail")
        assert "names no case" in message_of(e)
    else:
        assert data.AuditEntry(kind, "step", 126, "detail").case_ids == ()


@pytest.mark.parametrize("kind", data.KINDS)
def test_an_entry_that_touched_nothing_needs_no_identifier(kind):
    assert data.AuditEntry(kind, "step", 0, "detail").n == 0


def test_the_two_kinds_that_must_name_cases():
    # Roadmap Stage 2 demands the rule of corrections; an observation nobody can look up is not a
    # query, so it is held to the same rule. Pinned against the literal.
    assert data._MUST_NAME_CASES == frozenset({"correction", "observation"})


def test_an_unknown_kind_raises():
    with pytest.raises(ValueError) as e:
        data.AuditEntry("cleanup", "step", 0, "detail")
    assert "not an audit kind" in message_of(e)


def test_every_kind_has_a_section_heading():
    assert tuple(data._HEADINGS) == data.KINDS


def test_the_headings_are_the_declared_eight_in_pipeline_order():
    # Pinned against the literal, like SOURCES and _MUST_NAME_CASES: `derivation` is Stage 3's and
    # sits between `observation` and `structural`, because a document that printed the derivations
    # after the missingness table describing them would read backwards.
    #
    # `cohort` is Stage 5's and sits between `derivation` and `structural` [Stage 5 §6.2]. Neither of
    # KINDS' two properties is perfect there, because Stage 3 has two entry points on opposite sides
    # of Stage 5 and both record `derivation`: `derive`'s four entries run before the restrictions and
    # `derive_cohort`'s two after them. The position is wrong about two entries rather than four, and
    # it puts "who is in the analysis" immediately before the section reporting its denominators.
    assert data.KINDS == (
        "provenance", "contract", "correction", "observation",
        "derivation", "cohort", "structural", "missingness")
    assert list(data._HEADINGS.values()) == [
        "Provenance", "Contract", "Corrections", "Observations",
        "Derivations", "Cohort construction", "Structural non-applicability",
        "Missingness and denominators"]


# --- 12.8  byte-identical reproduction -------------------------------------------------------------

def test_two_loads_render_identical_markdown():
    first = data.load(data.FIXTURE)[1].to_markdown()
    second = data.load(data.FIXTURE)[1].to_markdown()
    assert first == second


def test_two_writes_produce_identical_bytes(tmp_path):
    a = data.load(data.FIXTURE)[1].write(tmp_path / "a.md")
    b = data.load(data.FIXTURE)[1].write(tmp_path / "b.md")
    assert a.read_bytes() == b.read_bytes()


def test_the_log_carries_no_timestamp_and_no_clock_time():
    # So the property cannot be passing merely because two runs landed in the same second.
    rendered = data.load(data.FIXTURE)[1].to_markdown()
    assert not re.search(r"\d{4}-\d{2}-\d{2}", rendered)
    assert not re.search(r"\d{2}:\d{2}:\d{2}", rendered)


def test_the_log_carries_no_absolute_path():
    rendered = data.load(data.FIXTURE)[1].to_markdown()
    assert str(config.ROOT) not in rendered
    assert str(config.FIXTURE_XLSX) not in rendered
    assert config.FIXTURE_XLSX.name in rendered      # the name, not the path


def test_the_log_is_identical_across_interpreters_with_different_hash_seeds():
    """A same-process comparison cannot see frozenset ordering, which is the one way this has broken.

    `DROPPED` is a frozenset and CPython randomises str hashing per process, so a drop list built by
    iterating it renders the contract entry's table in a different order in every interpreter —
    while passing every same-process byte comparison.
    """
    script = ("import sys, data\n"
              "sys.stdout.write(data.load(data.FIXTURE)[1].to_markdown())\n")
    rendered = []
    for seed in ("0", "1"):
        result = subprocess.run(
            [sys.executable, "-c", script],
            cwd=MODULE_DIR, capture_output=True, text=True,
            env={**os.environ, "PYTHONHASHSEED": seed})
        assert result.returncode == 0, result.stderr
        rendered.append(result.stdout)
    assert rendered[0] == rendered[1]


def test_load_writes_no_file(tmp_path, monkeypatch):
    logs = tmp_path / "logs"
    monkeypatch.setattr(config, "LOGS", logs)
    data.load(data.FIXTURE)
    assert not logs.exists(), "load() created the log directory; the caller decides, not load()"


def test_the_two_declared_sources_cannot_collide():
    assert data._audit_path(data.WORKBOOK) != data._audit_path(data.FIXTURE)
    assert data.WORKBOOK.label in data._audit_path(data.WORKBOOK).name
    assert data.FIXTURE.label in data._audit_path(data.FIXTURE).name


def test_write_creates_the_log_directory_and_returns_the_path(tmp_path, monkeypatch):
    logs = tmp_path / "logs"
    monkeypatch.setattr(config, "LOGS", logs)
    _, audit = data.load(data.FIXTURE)
    written = audit.write()
    # `audit_`, not `stage2_audit_`: the Audit spans the pipeline, and Stage 3 appends to this same
    # object rather than opening a second log.
    assert written == logs / f"audit_{data.FIXTURE.label}.md"
    assert written.read_text(encoding="utf-8") == audit.to_markdown()


# --- 12.9  schema assertions fire -------------------------------------------------------------------

_CORRUPTIONS = [
    ("A6", "age", lambda: corrupt("age", 12), "12"),
    ("A6", "nihss_baseline", lambda: corrupt("nihss_baseline", 50), "50"),
    ("A6", "mrs_90d", lambda: corrupt("mrs_90d", 7), "7"),
    ("A6", "prestroke_mrs", lambda: corrupt("prestroke_mrs", 6), "6"),
    ("A6", "core_ml", lambda: corrupt("core_ml", -1.0), "-1"),
    ("A5", "ivt", lambda: corrupt("ivt", 2), "2"),
    ("A4", "ivt", lambda: corrupt("ivt", pd.NA), "missing the treatment"),
    ("A4b", "ivt_contraindicated", lambda: corrupt("ivt_contraindicated", pd.NA),
     "missing the eligibility classifier"),
    ("A2", "case_id", lambda: corrupt("case_id", "HAND-2"), "HAND-2"),
    ("A2", "case_id", lambda: corrupt("case_id", " HAND-1 "), "HAND-1"),
    ("A9", "contraindication_reason", lambda: corrupt("contraindication_reason", "  "),
     "whitespace-only"),
    ("A7", "onset_to_ivt_min", lambda: corrupt("onset_to_ivt_min", 200.0, where="HAND-4"),
     "HAND-4"),
    ("A7", "onset_to_ivt_min", lambda: corrupt("onset_to_ivt_min", None, where="HAND-5"),
     "HAND-5"),
    ("A8", "onset_to_ivt_min", lambda: corrupt("onset_to_ivt_min", 300.0, where="HAND-5"),
     "HAND-5"),
    ("A3", "center", lambda: corrupt("center", "Zurich"), "Zurich"),
]


@pytest.mark.parametrize(
    "check, column, build, fragment", _CORRUPTIONS,
    ids=[f"{check}-{column}-{fragment}" for check, column, _, fragment in _CORRUPTIONS])
def test_one_corruption_raises_naming_the_column_and_the_value(
        check, column, build, fragment, tmp_path):
    with pytest.raises(config.SchemaError) as e:
        run(build(), hand_source(tmp_path))
    message = message_of(e)
    assert check in message
    assert column in message
    assert fragment in message


def test_two_corruptions_produce_one_error_naming_both(tmp_path):
    # Without this, §8.1's collect-all rule is untested and a fail-fast implementation passes every
    # other case above. A first-failure exception turns "three columns drifted in v8" into three
    # read-parse-fail round trips, with the operator learning one fact per run.
    df = corrupt("mrs_90d", 7)
    df.loc[df["case_id"] == "HAND-3", "core_ml"] = -1.0
    with pytest.raises(config.SchemaError) as e:
        run(df, hand_source(tmp_path))
    message = message_of(e)
    assert "mrs_90d" in message
    assert "core_ml" in message
    assert "2 assertion(s) failed" in message


def test_the_hand_frame_itself_passes_every_assertion(tmp_path):
    # Every case above corrupts this frame, so if it did not pass on its own the tests would be
    # asserting on the wrong violation.
    assert len(run(hand_frame(), hand_source(tmp_path))[0]) == 6


# --- 12.10  no bare assert in data.py ----------------------------------------------------------------

def _assert_lines(source: str) -> list[int]:
    return [node.lineno for node in ast.walk(ast.parse(source)) if isinstance(node, ast.Assert)]


def test_data_py_contains_no_bare_assert():
    # Python strips `assert` under -O, so under a flag nobody remembers setting a module whose
    # entire purpose is to fail loudly would succeed silently.
    offenders = _assert_lines(MODULE.read_text(encoding="utf-8"))
    assert not offenders, f"data.py has assert statements at line(s) {offenders}"


def test_the_assert_scan_actually_fires():
    # A scan that silently matches nothing would otherwise pass as a green test.
    assert _assert_lines("def f(x):\n    assert x > 0\n") == [2]


# --- 12.11  structural non-applicability [§11] --------------------------------------------------------

def _missingness_rows(audit: data.Audit) -> dict[str, tuple[str, ...]]:
    table = audit.entry("missingness", "absence_by_column").table
    return {row[0]: row for row in table[1:]}


def _missingness_header(audit: data.Audit) -> tuple[str, ...]:
    return audit.entry("missingness", "absence_by_column").table[0]


def test_the_two_kinds_of_absence_are_labelled_not_counted_as_missing(tmp_path):
    _, audit = run(hand_frame(), hand_source(tmp_path))
    rows = _missingness_rows(audit)
    assert rows["onset_to_ivt_min"][1] == "structural"
    assert rows["contraindication_reason"][1] == "not recorded"


@pytest.mark.parametrize(
    "column", sorted(set(config.STRUCTURALLY_NON_APPLICABLE) | set(config.INFORMATIVE_ABSENCE)))
def test_every_labelled_column_is_an_analysis_name_and_appears_with_its_reason(column, tmp_path):
    assert column in config.ANALYSIS_NAMES
    _, audit = run(hand_frame(), hand_source(tmp_path))
    entry = audit.entry("missingness", "absence_by_column")
    reason = {**config.STRUCTURALLY_NON_APPLICABLE, **config.INFORMATIVE_ABSENCE}[column]
    assert column in _missingness_rows(audit)
    assert reason in entry.detail


def test_a_structural_entry_is_recorded_per_declared_column(tmp_path):
    _, audit = run(hand_frame(), hand_source(tmp_path))
    steps = [e.step for e in audit.entries if e.kind == "structural"]
    assert steps == list(config.STRUCTURALLY_NON_APPLICABLE)


@DATA_GATED
def test_the_structural_count_equals_the_number_of_controls():
    df, audit = data.load()
    entry = audit.entry("structural", "onset_to_ivt_min")
    assert entry.n == int((df[config.TREATMENT] == 0).sum())


# --- 12.12  no imputation ----------------------------------------------------------------------------

def _raw_after_contract(source: data.Source) -> pd.DataFrame:
    audit = data.Audit(source)
    return data._apply_contract(data._read(source, audit), source, audit)


def _assert_no_imputation(source: data.Source) -> None:
    """Per column, not in aggregate: an aggregate count can net a gain against a loss."""
    before = _raw_after_contract(source)
    after, _ = data.load(source)
    for column in sorted(config.ANALYSIS_NAMES):
        expected = (before["core_ml"].isna() | before["tmax6_ml"].isna()
                    if column == "penumbra_ml" else before[column].isna())
        pd.testing.assert_series_equal(
            after[column].isna(), expected, check_names=False,
            obj=f"missingness pattern of {column}")


def test_the_fixture_gains_no_value():
    _assert_no_imputation(data.FIXTURE)


@DATA_GATED
def test_the_workbook_gains_no_value():
    _assert_no_imputation(data.WORKBOOK)


# --- 12.13  denominators ------------------------------------------------------------------------------

def test_the_missingness_table_has_one_row_per_analysis_column(tmp_path):
    _, audit = run(hand_frame(), hand_source(tmp_path))
    rows = _missingness_rows(audit)
    assert sorted(rows) == sorted(config.ANALYSIS_NAMES)
    assert len(rows) == 25


def test_the_centre_columns_come_from_center_order(tmp_path):
    # Never from the data: a run in which one centre contributes no rows would otherwise render a
    # narrower table that still reconciles, and the missing centre is the information.
    _, audit = run(hand_frame(), hand_source(tmp_path))
    assert _missingness_header(audit) == (
        "column", "kind", "n", "n_absent", "pct", *config.CENTER_ORDER)


def test_every_row_reconciles_to_the_frame(tmp_path):
    frame, audit = run(hand_frame(), hand_source(tmp_path))
    for column, row in _missingness_rows(audit).items():
        assert int(row[2]) + int(row[3]) == len(frame), column
        per_centre = [int(cell) for cell in row[5:]]
        assert sum(per_centre) == int(row[3]), f"{column}: per-centre counts do not sum"


def test_a_centre_with_no_rows_renders_a_column_of_zeros(tmp_path):
    # The fixture has two centres of four. The honest rendering is zeros, not a narrower table.
    _, audit = data.load(data.FIXTURE)
    header = _missingness_header(audit)
    assert len(header) == 5 + len(config.CENTER_ORDER)


# --- 12.14  structural facts from the workbook ---------------------------------------------------------
#
# The facts §7 and §8 were written against. DATA_SHA256 makes them safe; these tests document the
# dependency, so updating the hash also re-checks the premises.

@pytest.fixture(scope="module")
def workbook() -> tuple[pd.DataFrame, data.Audit]:
    return data.load()


@DATA_GATED
def test_the_full_read_is_126_by_42_and_a_header_only_read_loses_a_column():
    raw = data._read(data.WORKBOOK, data.Audit(data.WORKBOOK))
    assert raw.shape == (config.N_RECORDS_EXPECTED, config.N_COLUMNS_EXPECTED)
    header_only = pd.read_excel(config.DATA_XLSX, sheet_name=config.SHEET, nrows=0)
    assert len(header_only.columns) == config.N_COLUMNS_EXPECTED - 1


@DATA_GATED
def test_the_arms_are_39_treated_and_87_control(workbook):
    df, _ = workbook
    counts = df[config.TREATMENT].value_counts()
    assert int(counts.get(1)) == 39
    assert int(counts.get(0)) == 87


@DATA_GATED
def test_exactly_one_record_is_corrected_and_its_discrepancy_is_83_ml(workbook):
    _, audit = workbook
    entry = audit.entry("correction", "penumbra_recomputed")
    assert entry.n == 1
    assert len(entry.case_ids) == 1
    stored, recomputed, difference = entry.table[1][1:]
    assert (float(stored), float(recomputed), float(difference)) == (13.0, 96.0, 83.0)


@DATA_GATED
def test_exactly_one_record_carries_the_999_groin_time(workbook):
    _, audit = workbook
    entry = audit.entry("observation", "onset_to_groin_999")
    assert entry.n == 1
    assert len(entry.case_ids) == 1


# The two eligibility premises that used to sit here — no treated patient carries the flag, and the
# flag is never missing — moved to `test_eligibility.py` §12.13 when Stage 4 landed. They assert the
# classifier's preconditions, and one of them is Stage 4's headline assertion (E3); leaving them here
# would have meant the property was tested against one workbook and not against the classifier.


@DATA_GATED
def test_every_assertion_passes_against_the_workbook(workbook):
    df, _ = workbook
    assert data._assert_schema(df, data.WORKBOOK) is None


# --- 12.16  the audit log's inventory -------------------------------------------------------------------
#
# §12.8 compares a run against itself, which holds for *any* log — including one whose contents an
# implementer chose freely. This pins the entries themselves, so an implementation cannot quietly
# drop an entry, rename a step or reorder a section.

_INVENTORY_ON_THE_FIXTURE: list[tuple[str, str]] = [
    ("provenance", "read"),
    ("contract", "rename_and_drop"),
    ("provenance", "centre_recode"),
    ("correction", "penumbra_recomputed"),
    ("structural", "onset_to_ivt_min"),
    ("missingness", "absence_by_column"),
]


def test_the_fixture_log_contains_exactly_the_declared_entries():
    _, audit = data.load(data.FIXTURE)
    assert [(e.kind, e.step) for e in audit.entries] == _INVENTORY_ON_THE_FIXTURE


def test_the_structural_entries_are_one_per_declared_column():
    _, audit = data.load(data.FIXTURE)
    structural = [e.step for e in audit.entries if e.kind == "structural"]
    assert structural == list(config.STRUCTURALLY_NON_APPLICABLE)


def test_the_fixture_triggers_no_observation():
    _, audit = data.load(data.FIXTURE)
    assert audit.entry("observation", "onset_to_groin_999") is None


def test_penumbra_is_recorded_even_when_nothing_disagreed():
    # One of the five unconditional entries. With n = 0 it carries no identifiers and no table.
    _, audit = data.load(data.FIXTURE)
    entry = audit.entry("correction", "penumbra_recomputed")
    assert entry is not None and entry.n == 0 and entry.case_ids == ()


def test_every_heading_appears_and_an_empty_section_says_none():
    """The empty section is `observation`, and its *successor* comes from KINDS rather than a name.

    This test used to split on `structural` by name, which made it the one test in the repository
    that any new `kind` breaks: Stage 3 inserted `derivation` between the two, so the slice became
    `_none_ ## Derivations _none_` and the assertion failed. Re-pointing it at `derivation` would
    have been the same defect one position along — it rots on the next kind inserted between them.
    Taking the neighbour from KINDS is the repair that cannot recur.
    """
    rendered = data.load(data.FIXTURE)[1].to_markdown()
    for kind in data.KINDS:
        assert f"## {data._HEADINGS[kind]}" in rendered
    successor = data.KINDS[data.KINDS.index("observation") + 1]
    observations = rendered.split(f"## {data._HEADINGS['observation']}")[1]
    assert observations.split(f"## {data._HEADINGS[successor]}")[0].strip() == "_none_"


def test_the_header_agrees_with_the_entries_it_is_rendered_from():
    # A header recomputed independently of the entries is the one way these can disagree.
    df, audit = data.load(data.FIXTURE)
    read = audit.entry("provenance", "read")
    contract = audit.entry("contract", "rename_and_drop")
    n_mapped = len(config.RENAME)
    header = data.Audit.to_markdown(audit).split("\n")
    line = next(l for l in header if l.startswith("Read:"))
    assert line == (f"Read: {read.n} rows x {config.N_COLUMNS_EXPECTED} columns; "
                    f"{n_mapped} mapped, {contract.n} dropped")
    assert read.n == len(df)


def test_the_sections_render_in_kinds_order():
    rendered = data.load(data.FIXTURE)[1].to_markdown()
    positions = [rendered.index(f"## {data._HEADINGS[kind]}") for kind in data.KINDS]
    assert positions == sorted(positions)
