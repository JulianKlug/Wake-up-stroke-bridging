"""Stage 14 acceptance criteria — `specs/stage14_outputs_and_guardrails.md` §12."""

from __future__ import annotations

import ast
import hashlib
import os
from pathlib import Path
from typing import Final

import numpy as np
import pytest

import bootstrap
import config
import data
import report
from fixtures_stage14 import synthetic_run

MODULE = Path(report.__file__).resolve()
SOURCE = MODULE.read_text(encoding="utf-8")

DATA_GATED = pytest.mark.skipif(
    not config.DATA_XLSX.exists(),
    reason="the private workbook is gitignored and absent from this checkout")
SLOW = pytest.mark.skipif(
    os.environ.get("STAGE14_SLOW") != "1",
    reason="§12.12 is a full 2000-replicate run of Stages 1–13; open the gate with STAGE14_SLOW=1.")

# Captured on 2026-08-29 from a full run AFTER the pipe-in-header fixes and BEFORE the row-builder
# extraction (§8); a hash quotes no number.
STAGE13_AUDIT_SHA256: Final[str] = "d04c1f699cf7a54fe5381818cb7368b80ddf6f255b41b75d0c80f705075df27b"

TABLE_IDS = [i for i in config.OUTPUT_IDS if i.startswith("T")]
FIGURE_IDS = [i for i in config.OUTPUT_IDS if i.startswith("F")]


@pytest.fixture(scope="module")
def mixed(tmp_path_factory):
    run = synthetic_run("mixed", separation=5)
    out = tmp_path_factory.mktemp("mixed")
    return run, report.write(run, out), out


@pytest.fixture(scope="module")
def one_signed(tmp_path_factory):
    run = synthetic_run("one", separation=0)
    out = tmp_path_factory.mktemp("one")
    return run, report.write(run, out), out


def table(out: Path, oid: str, suffix: str = ".md") -> str:
    return (out / "tables" / f"{oid}_{config.OUTPUT_IDS[oid]}{suffix}").read_text(encoding="utf-8")


def figure(out: Path, oid: str) -> str:
    return (out / "figures" / f"{oid}_{config.OUTPUT_IDS[oid]}.svg").read_text(encoding="utf-8")


# --- 12.1  surface ---------------------------------------------------------------------------------

def test_the_public_surface_is_THIRTEEN_NAMES_and_TWO_dataclasses_in_source_order():
    """`run`, `write` and `main` drive the stage; the other five are shared with the manuscript.

    The five were private until `figures_and_tables/` needed them, and they are public for
    `balance.levels`'s reason [Stage 7 §5.5]: a manuscript exhibit carrying a [§16] sentence it
    composed itself is a second definition of what the report says, and the failure would be a
    sentence that stops tracking the run. `baseline_rows` and `binary_rows` are promoted by
    EXTRACTION and `_t02`, `_t10` and `_t11` are still their only callers here, so neither can rot
    into a manuscript-only branch. `baseline_summary` and `outcome_summary` are the two names here this stage
    does not itself call: T02 keeps the arm MEANS its [§9] standardised mean differences summarise, and the clinical
    `n (%)` / `median (IQR)` presentation belongs to the manuscript. Both walk `C.BALANCE_SET` through
    `balance.levels`, and `tests/test_manuscript.py` pins that they produce the same rows. Nothing else was
    promoted — the row builders, the figure
    builders and the renderers stay private, so the manuscript cannot reach past the statements.
    """
    tree = ast.parse(SOURCE)
    functions = [n.name for n in tree.body if isinstance(n, ast.FunctionDef)]
    assert [n for n in functions if not n.startswith("_")] == [
        "run", "direction_clause", "exceedance_clause", "pmf", "interval_note", "baseline_rows",
        "baseline_summary", "outcome_summary", "binary_rows", "arm_clauses", "render_csv",
        "write", "main"]
    assert [n.name for n in tree.body if isinstance(n, ast.ClassDef)] == ["Run", "Manifest"]


def test_report_written_is_disjoint_from_every_other_step_and_reaches_only_data_privates():
    others = set()
    for path in MODULE.parent.glob("*.py"):
        if path.name in ("report.py", "stage0_data_inventory.py"):
            continue
        others |= {n.value for n in ast.walk(ast.parse(path.read_text(encoding="utf-8")))
                   if isinstance(n, ast.Constant) and isinstance(n.value, str)}
    assert report._STEP_WRITTEN not in others
    calls = {ast.unparse(c.func) for c in ast.walk(ast.parse(SOURCE)) if isinstance(c, ast.Call)}
    assert {c for c in calls if "._" in c} <= {"data._assert_no_pipe"}


# --- 12.2  coverage --------------------------------------------------------------------------------

def test_every_output_id_lands_once_and_the_manifest_names_every_file(mixed):
    _, manifest, out = mixed
    assert len(manifest.paths) == 2 * len(TABLE_IDS) + len(FIGURE_IDS) + 1
    for path in manifest.paths:
        assert (out / path).exists(), path
    for oid in TABLE_IDS:
        assert table(out, oid) and table(out, oid, ".csv")
    for oid in FIGURE_IDS:
        assert figure(out, oid).startswith("<?xml")


# --- 12.3  the direction clause tracks the data ----------------------------------------------------

def test_the_direction_clause_is_ABSENT_when_every_RD_k_shares_a_sign_and_PRESENT_otherwise(mixed, one_signed):
    _, _, out_mixed = mixed
    _, _, out_one = one_signed
    assert report.CONSTANT_SHIFT in table(out_one, "T07")
    assert report.DIRECTION not in table(out_one, "T07")
    assert report.CONSTANT_SHIFT in table(out_mixed, "T07")
    assert report.DIRECTION in table(out_mixed, "T07")


def test_the_direction_clause_is_computed_from_rd_and_zero_is_its_own_sign():
    assert report.direction_clause({0: 0.1, 1: 0.2}) == ""
    assert report.direction_clause({0: 0.1, 1: -0.2}) == report.DIRECTION
    assert report.direction_clause({0: 0.0, 1: 0.2}) == report.DIRECTION


# --- 12.4  the thinning clause tracks the counters ------------------------------------------------

def test_the_thinning_clause_names_surviving_counts_only_where_separation_removed_a_replicate(mixed, one_signed):
    assert report.THINNING in table(mixed[2], "T13")
    assert "55 of 60 survive" in table(mixed[2], "T13")
    assert report.THINNING not in table(one_signed[2], "T13")
    assert report.THINNING in table(mixed[2], "T19") and report.THINNING not in table(one_signed[2], "T19")


# --- 12.4a  the sigma boundary rate tracks the draws ----------------------------------------------

def test_the_sigma_floor_rate_counts_draws_AT_the_floor_and_never_merely_NEAR_it():
    floor = config.POLR_RI_SIGMA_FLOOR
    at = np.array([floor, floor, 2.0 * floor, 0.5])
    assert report._floor_rate(bootstrap.Draws("hier.sigma", at, 4, {})) == "2 of 4"
    assert report._floor_rate(bootstrap.Draws("hier.sigma", np.array([0.5]), 1, {})) == "0 of 1"
    assert report._floor_rate(None) == report._DASH


def test_T15_names_the_sigma_BOUNDARY_RATE_and_not_only_whether_the_POINT_estimate_is_at_the_floor(mixed):
    rendered = table(mixed[2], "T15")
    assert "replicates at the floor: 0 of 60" in rendered
    assert "is not a value for the between-centre SD" in rendered
    assert "Point estimate at the floor: no" in rendered


# --- 12.5  [§14] beside the primary ----------------------------------------------------------------

def test_T22_and_F05_carry_the_populations_statement_and_T07_does_not(mixed):
    run, _, out = mixed
    for text in (table(out, "T22"), figure(out, "F05")):
        assert report.POPULATIONS in text and report.TRANSPORTED in text
    assert report.POPULATIONS not in table(out, "T07") and report.TRANSPORTED not in table(out, "T07")
    header = table(out, "T22").split("\n")[0]
    assert "odds_ratio" not in header and "eligible_rd" not in header and "support.beta" not in header
    assert " p " not in table(out, "T18").split("\n")[0]
    assert run.pol.label in table(out, "T18") and run.pol.measure in table(out, "T18")
    assert "support.beta" not in table(out, "T15")


# --- 12.6  labels are read, not typed --------------------------------------------------------------

def test_augmented_rows_are_model_assisted_safety_rows_descriptive_subgroups_hypothesis_generating(mixed):
    run, _, out = mixed
    t10 = [line for line in table(out, "T10").split("\n") if line.startswith("| ") and "---" not in line][1:]
    for line in t10:
        assert report.MODEL_ASSISTED in table(out, "T10").split("\n")[0]
    for est in run.secondary.by_family()["safety"]:
        row = next(line for line in table(out, "T11").split("\n") if config.OUTCOMES[est.outcome].label in line)
        assert f"| {report.DESCRIPTIVE_SAFETY} " in row
    assert report.HYPOTHESIS_GENERATING in table(out, "T13") and report.HYPOTHESIS_GENERATING in figure(out, "F03")
    assert run.e_value.approximation in table(out, "T09")
    assert report.E_VALUE_HEURISTIC in table(out, "T09")
    assert "diabetes (+0.250)" in table(out, "T09") and "diabetes (+0.250)" in table(out, "T07")


def test_the_e_value_range_names_exactly_the_thresholds_outside_the_licence():
    cumulative = {k: {0: p, 1: p} for k, p in zip(config.MRS_THRESHOLDS, (0.05, 0.2, 0.5, 0.7, 0.9, 0.95))}
    clause = report._e_value_range(cumulative)
    assert "RD_0" in clause and "RD_4" in clause and "RD_5" in clause
    assert "RD_1" not in clause and "RD_2" not in clause and "RD_3" not in clause


# --- 12.7  denominators ----------------------------------------------------------------------------

def test_every_estimate_row_carries_a_non_empty_denominator(mixed):
    run, _, out = mixed
    n = str(int(run.primary.in_estimate.sum()))
    for oid in ("T07", "T10", "T11"):
        lines = table(out, oid).split("\n")
        header = [c.strip() for c in lines[0].strip("|").split("|")]
        col = header.index("denominator")
        body = [line for line in lines[2:] if line.startswith("| ")]
        assert body
        for line in body:
            cell = [c.strip() for c in line.strip("|").split("|")][col]
            assert cell == n or cell == "—", (oid, line)
    assert f"| {n} " in table(out, "T13")


# --- 12.8  crude table -----------------------------------------------------------------------------

def test_T03_carries_the_descriptive_label_and_no_difference_column(mixed):
    text = table(mixed[2], "T03")
    assert report.DESCRIPTIVE_CRUDE in text
    assert "difference" not in text.split("\n")[0].lower()
    assert report.DESCRIPTIVE_CRUDE in table(mixed[2], "T03", ".csv")


# --- 12.9  rendering -------------------------------------------------------------------------------

def test_no_cell_contains_a_pipe_every_row_has_the_header_width_and_csv_agrees(mixed):
    run, _, out = mixed
    for oid in TABLE_IDS:
        rows, clauses = report._TABLES[oid](run)
        assert len({len(r) for r in rows}) == 1, oid
        assert not any("|" in cell for row in rows for cell in row), oid
        md_body = [line for line in table(out, oid).split("\n") if line.startswith("| ")]
        assert len({line.count("|") for line in md_body}) == 1, oid
        csv_rows = [line for line in table(out, oid, ".csv").split("\n") if line and not line.startswith("# ")]
        assert len(csv_rows) == len(md_body) - 1 == len(rows), oid   # md has the separator row
        assert sum(1 for line in table(out, oid, ".csv").split("\n") if line.startswith("# ")) == len(clauses)


def test_V1_fires_on_a_pipe_in_a_cell():
    with pytest.raises(config.SchemaError) as excinfo:
        data.coefficient_rows(("a|b",), (1.0,), (), lambda name: "role")
    assert str(excinfo.value).startswith("V1")
    with pytest.raises(config.SchemaError):
        data.replicate_rows([("attempted = 1", [("q", "x|y")])])


# --- 12.10  checklist ------------------------------------------------------------------------------

def test_every_checklist_id_is_an_output_and_V3_fires_on_an_empty_item(monkeypatch, mixed):
    for _, ids in config.CHECKLIST:
        assert ids and set(ids) <= set(config.OUTPUT_IDS)
    covered = {i for _, ids in config.CHECKLIST for i in ids}
    assert covered == set(config.OUTPUT_IDS)
    monkeypatch.setattr(config, "CHECKLIST", (("orphan", ()),))
    with pytest.raises(config.SchemaError) as excinfo:
        report._t21(mixed[0])
    assert str(excinfo.value).startswith("V3")


# --- 12.11  determinism ----------------------------------------------------------------------------

def test_two_writes_of_one_run_yield_identical_digests_and_the_log_names_them(tmp_path):
    run = synthetic_run("mixed", separation=5)
    a = report.write(run, tmp_path / "a")
    b = report.write(run, tmp_path / "b")
    assert a.sha256 == b.sha256
    assert sum(1 for e in run.audit.entries if e.step == report._STEP_WRITTEN) == 1
    written = run.audit.entry("provenance", report._STEP_WRITTEN)
    assert written.n == len(a.paths) - 1                     # the log itself is not in its own table
    assert len(written.table) == written.n + 1
    for path in a.paths:
        text = (tmp_path / "a" / path).read_bytes()
        assert b"Date" not in text or path.suffix != ".svg"
        assert str(tmp_path).encode() not in text


def test_V4_fires_on_a_None_field():
    import dataclasses
    run = synthetic_run()
    with pytest.raises(config.SchemaError) as excinfo:
        dataclasses.replace(run, pol_boot=None)
    assert str(excinfo.value).startswith("V4")


# --- 12.13  nothing leaks --------------------------------------------------------------------------

def test_outputs_land_under_out_and_no_tracked_file_quotes_a_float(mixed):
    _, manifest, _ = mixed
    for path in manifest.paths:
        assert not path.is_absolute() and path.parts[0] in ("tables", "figures", "logs")
    import re
    spec = (MODULE.parent / "specs" / "stage14_outputs_and_guardrails.md").read_text(encoding="utf-8")
    for text, name in ((SOURCE, "report.py"), (spec, "spec")):
        assert not re.search(r"\d\.\d{4,}", text), name


# --- 12.14  no fitter ------------------------------------------------------------------------------

def test_report_imports_no_fitter_and_calls_only_the_entry_points_run_names():
    tree = ast.parse(SOURCE)
    imported = {a.name.split(".")[0] for n in ast.walk(tree) if isinstance(n, ast.Import) for a in n.names}
    imported |= {n.module.split(".")[0] for n in ast.walk(tree) if isinstance(n, ast.ImportFrom) and n.module}
    assert not imported & {"model", "scipy", "statsmodels"}
    calls = {ast.unparse(c.func) for c in ast.walk(tree) if isinstance(c, ast.Call)}
    allowed = {
        "propensity.fit", "outcome.primary", "outcome.secondary", "bootstrap.run",
        "sensitivity.multiplicity", "sensitivity.e_value_primary", "sensitivity.subgroups",
        "sensitivity.full_covariate", "standardise.population", "standardise.all_centre",
        "standardise.support", "standardise.hierarchical", "standardise.inference",
        "policy.population", "policy.contrast", "policy.inference"}
    stage_calls = {c for c in calls if c.split(".")[0] in {
        "propensity", "outcome", "bootstrap", "sensitivity", "standardise", "policy"}}
    assert stage_calls <= allowed, stage_calls - allowed


# --- 12.12 / 12.15  upstream log unchanged, the frame unmutated — [slow] [data-gated] -------------

@SLOW
@DATA_GATED
def test_STAGE13_AUDIT_SHA256_pins_the_full_log_across_the_row_builder_extraction(tmp_path):
    """The gate is a run: `STAGE14_SLOW=1 uv run pytest -k STAGE13_AUDIT -rs` must report passed."""
    run = report.run()
    log = run.audit.write(tmp_path / "audit.md")
    assert hashlib.sha256(log.read_bytes()).hexdigest() == STAGE13_AUDIT_SHA256
    assert len(run.audit) == 61


@DATA_GATED
def test_run_returns_the_classified_frame_unmutated_and_writes_nothing(tmp_path, monkeypatch):
    import derive
    import eligibility
    original = config.N_BOOT
    config.N_BOOT = 60
    try:
        run = report.run()
    finally:
        config.N_BOOT = original
    df, _ = data.load()
    expected = eligibility.classify(derive.derive(df, data.Audit(data.WORKBOOK)), data.Audit(data.WORKBOOK))
    import pandas.testing as pdt
    pdt.assert_frame_equal(run.df, expected)
    manifest = report.write(run, tmp_path)
    assert len(manifest.paths) == 2 * len(TABLE_IDS) + len(FIGURE_IDS) + 1
    assert len(run.audit) == 62
