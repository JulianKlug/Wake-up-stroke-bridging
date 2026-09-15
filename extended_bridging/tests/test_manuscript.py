"""Acceptance tests for the manuscript layer — `figures_and_tables/manuscript.py`.

Data-free: everything here runs on `fixtures_stage14.synthetic_run`, the same synthetic `Run` the
Stage 14 tests use. What is tested is not a number; it is that the second reporting surface cannot
drift from the first. Stage 14's guardrails were written down *"so they cannot be lost in editing"*
(`implementation_roadmap.md`, Stage 14), and a manuscript exhibit is exactly where they would be: a
figure caption composed by hand reads correctly on the day it is written and stops tracking the run
the day after.
"""

from __future__ import annotations

import ast
import dataclasses
from pathlib import Path

import matplotlib
import pytest

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

import config  # noqa: E402
import manuscript  # noqa: E402
import report  # noqa: E402
from fixtures_stage14 import synthetic_run  # noqa: E402

MODULE = Path(manuscript.__file__).resolve()
SOURCE = MODULE.read_text(encoding="utf-8")

# What `manuscript.py` may import at module level. Its neighbour below is `report.py`; an estimator
# imported here is this layer reaching past it, and `config`/`data` are the contract both families
# read. The plotting and frame libraries are how it draws and formats, and nothing else.
# `IPython` is the notebook's display, imported inside `show`; the result types are named THROUGH
# `report`, so describing one does not require importing it.
ALLOWED_IMPORTS = {
    "config", "data", "report", "matplotlib", "numpy", "pandas", "IPython",
    "__future__", "hashlib", "pickle", "sys", "dataclasses", "enum", "pathlib", "typing"}


@pytest.fixture(scope="module")
def run():
    return synthetic_run("mixed", separation=5)


@pytest.fixture(scope="module")
def one_signed():
    return synthetic_run("one", separation=0)


def bundle_of(run) -> manuscript.Bundle:
    """The synthetic `Run`'s Stage 1-11 half, as a `Bundle`. The field names are the same by design."""
    return manuscript.Bundle(
        manuscript.Provenance.current(),
        *(getattr(run, f.name) for f in dataclasses.fields(manuscript.Bundle)[1:]))


# --- the two families stay apart --------------------------------------------------------------------

def test_no_manuscript_id_reaches_the_report_output_contract():
    """`C.OUTPUT_IDS` is what `report.write` iterates, and it may not learn about this family.

    An id added there reaches `report.write`, which has no builder for it and raises KeyError, and
    `C.CHECKLIST` — which may only name ids from that dict — would be claiming a STROBE item was
    discharged by a file Stage 14 never wrote.
    """
    assert not set(manuscript.MANUSCRIPT_IDS) & set(config.OUTPUT_IDS)
    discharging = {oid for _, ids in config.CHECKLIST for oid in ids}
    assert not set(manuscript.MANUSCRIPT_IDS) & discharging


def test_the_manuscript_layer_imports_no_estimator():
    """Stage 6 §12.7's rule, applied one directory down: this layer talks to `report`, not past it."""
    imported = set()
    for node in ast.walk(ast.parse(SOURCE)):
        if isinstance(node, ast.Import):
            imported.update(alias.name.split(".")[0] for alias in node.names)
        elif isinstance(node, ast.ImportFrom) and node.module:
            imported.add(node.module.split(".")[0])

    assert imported <= ALLOWED_IMPORTS, f"unexpected imports: {sorted(imported - ALLOWED_IMPORTS)}"


def test_every_balance_covariate_has_a_human_label():
    """A covariate added to [§9]'s set and not named here would reach a tick as its column name."""
    assert set(manuscript.COVARIATE_LABELS) == set(config.BALANCE_SET)
    assert manuscript.label("center = HUG") == "Centre: HUG"
    assert manuscript.label("onset_type = wake_up") == "Onset type: wake-up"


def test_an_id_that_is_not_an_exhibit_is_refused():
    with pytest.raises(manuscript.ExhibitError):
        manuscript._stem("Figure_9")


def test_the_promoted_report_surface_is_what_this_layer_calls():
    """Every name the manuscript reads off `report` is public there, and none of them is a builder."""
    for name in ("direction_clause", "exceedance_clause", "interval_note", "pmf", "baseline_rows",
                 "binary_rows", "arm_clauses", "render_csv", "fmt", "CONSTANT_SHIFT",
                 "ATO_DESCRIPTION", "MODEL_ASSISTED", "DESCRIPTIVE_SAFETY"):
        assert hasattr(report, name), name
    for gone in ("_direction_clause", "_exceedance_clause", "_interval_note", "_pmf",
                 "_binary_rows", "_render_csv"):
        assert not hasattr(report, gone), f"{gone} survived the promotion under both spellings"


# --- the cache is keyed on what moves a number ---------------------------------------------------------

@pytest.mark.parametrize("field", ["schema", "data_sha256", "modules", "seed", "n_boot", "ci_level",
                                   "percentile_method", "boot_stratum", "python", "numpy", "pandas"])
def test_the_cache_key_moves_when_any_component_moves(field):
    """A 60-replicate cache must never be readable as a 2000-replicate one, and so on for all eleven.

    `tests/test_report.py` monkeypatches `N_BOOT = 60` for its slow gate, which is precisely the
    accident this prevents: a bundle built under that patch and left on disk.
    """
    current = manuscript.Provenance.current()
    moved = dataclasses.replace(current, **{field: "moved" if isinstance(
        getattr(current, field), str) else getattr(current, field) + 1})

    assert moved.key() != current.key()
    assert current.differences(moved) == (field,)


def test_a_bundle_missing_a_stage_is_refused(run):
    with pytest.raises(config.SchemaError):
        dataclasses.replace(bundle_of(run), arm=None)


def test_the_bundle_names_its_fields_as_the_run_does(run):
    """Field-name identity is what lets `report`'s row builders read a `Bundle` unchanged."""
    shared = {f.name for f in dataclasses.fields(manuscript.Bundle)} - {"provenance"}
    assert shared <= {f.name for f in dataclasses.fields(report.Run)}
    assert report.baseline_rows(run.cohort, run.ps) == report.baseline_rows(
        bundle_of(run).cohort, bundle_of(run).ps)


def test_the_two_baseline_tables_describe_the_same_rows(run):
    """T02's means and the manuscript's `n (%)` / `median (IQR)` are two readings of one row set.

    They are allowed to present a covariate differently — that is the point of having both — but not
    to disagree about which covariates there are, in what order, or over which declared levels. Both
    walk `C.BALANCE_SET` through `balance.levels`, and this is what keeps them doing so.
    """
    means = report.baseline_rows(run.cohort, run.ps)
    clinical = report.baseline_summary(run.cohort, run.ps)

    assert [row[0] for row in means[1:]] == [row[0] for row in clinical[1:]]
    assert means[0][0] == clinical[0][0] == "covariate"
    assert [row[0] for row in clinical[1:]][-1] == "effective sample size"


def test_every_clinical_row_says_which_summary_it_is(run):
    summaries = {row[1] for row in report.baseline_summary(run.cohort, run.ps)[1:]}

    assert summaries == {"n (%)", "median (IQR)", "ESS"}


def test_counted_rows_are_declared_and_not_inferred_from_the_sample(run):
    """A 0/1 covariate is counted because `config` says its values are in {0, 1}, not because they are.

    Inferring it would format by what this cohort happens to hold: a continuous covariate that came in
    all-zero on one workbook would print as `n (%)` there and `median (IQR)` on the next.
    """
    rows = {row[0]: row[1] for row in report.baseline_summary(run.cohort, run.ps)[1:]}

    for name in ("sex", "atrial_fib", "hypertension"):
        assert rows[name] == "n (%)", name
    for name in ("age", "nihss_baseline", "core_ml"):
        assert rows[name] == "median (IQR)", name
    for level in ("witnessed", "unwitnessed", "wake_up"):
        assert rows[f"onset_type = {level}"] == "n (%)", level


# --- the [§16] guardrails, carried onto the second surface -----------------------------------------------

def test_every_exhibit_carrying_exp_beta_carries_the_constant_shift_statement(run):
    assert report.CONSTANT_SHIFT in manuscript.primary_clauses(run)


def test_the_direction_statement_is_computed_and_not_printed_unconditionally(run, one_signed):
    """The mitigation [§8] prescribes only works if the reader is told what it is for [§16a]."""
    assert report.DIRECTION in manuscript.primary_clauses(run)
    assert report.DIRECTION not in manuscript.primary_clauses(one_signed)
    assert report.CONSTANT_SHIFT in manuscript.primary_clauses(one_signed)


def test_the_residual_imbalance_and_the_interval_note_travel_with_the_estimate(run):
    clauses = manuscript.primary_clauses(run)
    assert any(clause.startswith(("Residual imbalance", "No covariate exceeds")) for clause in clauses)
    assert any(clause.startswith("Intervals:") for clause in clauses)


def test_safety_is_labelled_descriptive_and_augmented_is_labelled_model_assisted(run):
    rows, clauses = report.binary_rows(run, "safety")
    assert all(row[1] == report.DESCRIPTIVE_SAFETY for row in rows[1:] if row[1] not in
               ("no p-value", "not estimable"))
    assert any(report.DESCRIPTIVE_SAFETY.upper() in clause for clause in clauses)

    header = report.binary_rows(run, "secondary")[0][0]
    assert report.MODEL_ASSISTED in header[11]


def test_every_binary_row_carries_its_own_denominator(run):
    for family in ("secondary", "safety"):
        rows = report.binary_rows(run, family)[0]
        assert rows[0][2] == "denominator"
        assert all(row[2] for row in rows[1:])


def test_the_provenance_clause_names_the_run_and_the_drawing_library(run):
    clause = manuscript.provenance_clause(bundle_of(run))
    assert str(config.SEED) in clause and matplotlib.__version__ in clause


# --- what lands on disk --------------------------------------------------------------------------------

def test_save_table_writes_the_rows_then_the_clauses(tmp_path, monkeypatch, run):
    monkeypatch.setattr(manuscript, "OUT", tmp_path)
    frame = _frame()

    path = manuscript.save_table("Table_1", frame, ("a caption", "", report.ATO_DESCRIPTION))
    text = path.read_text(encoding=manuscript.CSV_ENCODING)

    assert path.name == "Table_1_baseline_and_weighting.csv"
    assert path.read_bytes().startswith(b"\xef\xbb\xbf"), "Excel needs the BOM to read UTF-8"
    assert text.startswith("quantity,value\r\n") or text.startswith("quantity,value\n")
    assert "# a caption" in text and f"# {report.ATO_DESCRIPTION}" in text
    assert "# \n" not in text, "an empty clause was rendered as a blank comment"


def test_save_figure_writes_the_png_and_its_source_data(tmp_path, monkeypatch):
    monkeypatch.setattr(manuscript, "OUT", tmp_path)
    fig, ax = plt.subplots()
    ax.plot([0, 1], [0, 1])

    path = manuscript.save_figure("Figure_1", fig, _frame(), ("Figure 1. A caption.",))

    assert path.name == "Figure_1_cohort_flow.png"
    assert path.read_bytes()[:8] == b"\x89PNG\r\n\x1a\n"
    source = (tmp_path / f"Figure_1_cohort_flow{manuscript.SOURCE_DATA_SUFFIX}").read_text(
        manuscript.CSV_ENCODING)
    assert "# Figure 1. A caption." in source
    assert not plt.get_fignums(), "save_figure left the figure open"


def test_save_slide_writes_the_deck_source_beside_the_exhibit(tmp_path, monkeypatch):
    monkeypatch.setattr(manuscript, "OUT", tmp_path)

    path = manuscript.save_slide("Figure_1", "---\ntitle: t\n---\n")

    assert path.name == "Figure_1_cohort_flow.md"
    assert path.read_text("utf-8").startswith("---")


def test_save_legends_collects_the_captions_and_hides_the_provenance(tmp_path, monkeypatch):
    monkeypatch.setattr(manuscript, "OUT", tmp_path)
    fig, _ = plt.subplots()
    manuscript.save_figure("Figure_1", fig, _frame(),
                           ("Figure 1. A caption.", "A required statement [§11].",
                            f"{manuscript.PROVENANCE_PREFIX} DATA_SHA256 abc, seed 1, 2 replicates."))

    text = manuscript.save_legends().read_text("utf-8")

    assert "## Figure 1 — cohort_flow" in text
    assert "Figure 1. A caption." in text
    assert "A required statement [§11]." in text
    # The provenance says which run drew the figure; a printed legend does not carry it.
    assert f"<!-- {manuscript.PROVENANCE_PREFIX}" in text
    assert not any(line.startswith(manuscript.PROVENANCE_PREFIX) for line in text.splitlines())


def test_save_legends_skips_a_figure_that_has_not_been_written(tmp_path, monkeypatch):
    monkeypatch.setattr(manuscript, "OUT", tmp_path)
    fig, _ = plt.subplots()
    manuscript.save_figure("Figure_2", fig, _frame(), ("Figure 2. Only this one.",))

    text = manuscript.save_legends().read_text("utf-8")

    assert "## Figure 2 — design_diagnostics" in text
    assert "Figure_1" not in text


def test_the_manifest_names_every_file_of_every_exhibit(tmp_path, monkeypatch):
    monkeypatch.setattr(manuscript, "OUT", tmp_path)
    fig, _ = plt.subplots()
    manuscript.save_figure("Figure_1", fig, _frame(), ("Figure 1. A caption.",))
    manuscript.save_table("Table_1", _frame(), ("Table 1. A different caption.",))

    written = manuscript.written()
    assert set(written["exhibit"]) == {"Figure_1", "Table_1"}
    assert len(written) == 3        # png, its source data, and the table
    assert written["sha256"].nunique() == 3


def _frame():
    import pandas as pd

    return pd.DataFrame({"quantity": ["a", "b"], "value": [1.5, 2.0]})


# --- the exhibits carry the run's numbers, not numbers that look like them ---------------------------
#
# Gated on both output families being present, because the manuscript CSVs are built from the private
# workbook and `out/` is gitignored. This is the check the whole layer exists to pass: an exhibit is a
# view of a locked analysis output, and the way it would fail is by drifting one cell at a time.

_TABLES = config.OUT / "tables"
_MANUSCRIPT = config.OUT / "manuscript"

WRITTEN = pytest.mark.skipif(
    not (_TABLES / "T07_primary.csv").exists() or not (_MANUSCRIPT / "Table_S1_primary_and_sensitivity.csv").exists(),
    reason="needs both `python -m report` and the seven notebooks to have run against the workbook")


def _read(path: Path):
    import pandas as pd

    return pd.read_csv(path, comment="#", dtype=str)      # pandas strips the BOM itself


@WRITTEN
def test_table_S1_reproduces_T07_cell_for_cell():
    table = _read(_MANUSCRIPT / "Table_S1_primary_and_sensitivity.csv")
    primary = _read(_TABLES / "T07_primary.csv")
    effect = table[table["row"].str.startswith("primary") & (table["block"] == "effect")]

    assert len(effect) == len(primary)
    for mine, theirs in zip(effect.itertuples(), primary.itertuples()):
        assert (mine.quantity, mine.estimate, mine.draws) == (theirs.quantity, theirs.estimate,
                                                              theirs.draws)


@WRITTEN
def test_table_S1s_sensitivity_block_reproduces_T12():
    table = _read(_MANUSCRIPT / "Table_S1_primary_and_sensitivity.csv")
    arm = _read(_TABLES / "T12_sensitivity_arm.csv")
    arm = arm[arm["row"].str.startswith("full-covariate")].iloc[0]
    mine = table[table["row"].str.startswith("full-covariate")].set_index("quantity")["estimate"]

    assert mine["common odds ratio exp(beta)"] == arm["exp(beta)"]
    assert mine[mine.index.str.startswith("worst residual")].item() == arm["worst residual SMD"]


@WRITTEN
def test_table_1_covers_T02s_rows_and_reproduces_T05s_balance():
    """The cells are a different presentation of T02's rows; the SMD columns are T05's, rounded.

    Table 1 reports `n (%)` and `median (IQR)` where T02 reports arm means, so the cells are NOT
    comparable — that is the point of having both. What must still hold is that it describes the same
    covariates, over the same declared levels, in the same order, and that its balance columns are
    T05's numbers and not a second computation of them.
    """
    table = _read(_MANUSCRIPT / "Table_1_baseline_and_weighting.csv")
    baseline = _read(_TABLES / "T02_baseline_by_arm.csv")
    smd = _read(_TABLES / "T05_balance_smd.csv")

    assert table["covariate"].tolist() == baseline["covariate"].tolist()
    assert table["n"].tolist() == baseline["n"].tolist()

    rows = table[table["SMD weighted"] != "—"]
    assert rows["exceeds"].tolist() == smd["exceeds"].tolist()
    for column in ("SMD unweighted", "SMD weighted"):
        assert rows[column].tolist() == [f"{float(v):.3f}" for v in smd[column]], column


@WRITTEN
def test_table_S2_is_T10_and_T11_in_one_table():
    import pandas as pd

    table = _read(_MANUSCRIPT / "Table_S2_secondary_and_safety.csv")
    both = pd.concat([_read(_TABLES / "T10_secondary.csv"), _read(_TABLES / "T11_safety.csv")],
                     ignore_index=True)

    for column in ("outcome", "label", "denominator", "RD", "p raw", "p adjusted", "odds ratio"):
        assert table[column].tolist() == both[column].tolist(), column


@WRITTEN
def test_figure_3s_source_data_reproduces_T08_and_T07():
    source = _read(_MANUSCRIPT / "Figure_3_primary_mrs_source_data.csv")
    distribution = _read(_TABLES / "T08_mrs_distribution_by_arm.csv")
    primary = _read(_TABLES / "T07_primary.csv")

    for arm, column in (("EVT alone", "P(Y=j) EVT alone"), ("bridging", "P(Y=j) bridging")):
        plotted = source[(source["panel"] == "A") & (source["row"] == arm)]["value"].tolist()
        assert plotted == distribution[column].tolist(), arm

    plotted = source[(source["panel"] == "B") & (source["quantity"] == "RD_k")]["value"].tolist()
    assert plotted == primary[primary["quantity"].str.startswith("RD_")]["estimate"].tolist()


@WRITTEN
def test_figure_2s_source_data_reproduces_T05_and_T06():
    source = _read(_MANUSCRIPT / "Figure_2_design_diagnostics_source_data.csv")
    smd = _read(_TABLES / "T05_balance_smd.csv")
    overlap = _read(_TABLES / "T06_overlap_by_centre.csv").set_index("centre")

    # Panel A is labelled, and it shows the covariates weighting was expected to act on: every row is
    # T05's number for that covariate, and the collinear one T05 also carries is deliberately absent.
    plotted = source[(source["panel"] == "A") & (source["quantity"] == "SMD weighted")]
    by_label = dict(zip(smd["covariate"].map(manuscript.label), smd["SMD weighted"]))

    assert plotted["value"].tolist() == [by_label[row] for row in plotted["row"]]
    assert manuscript.label("penumbra_ml") not in set(plotted["row"])
    assert "penumbra_ml" in set(smd["covariate"]), "T05 still reports the full [§9] set"

    for centre in source[source["panel"] == "B"]["row"].unique():
        rows = source[(source["panel"] == "B") & (source["row"] == centre)].set_index("quantity")
        assert rows.loc["n bridging", "value"] == overlap.loc[centre, "n bridging"]
        assert rows.loc["ESS EVT alone", "value"] == overlap.loc[centre, "ESS EVT alone"]


@WRITTEN
def test_figure_1s_source_data_reproduces_T01():
    source = _read(_MANUSCRIPT / "Figure_1_cohort_flow_source_data.csv")
    flow = _read(_TABLES / "T01_cohort_flow.csv")

    assert source["records"].tolist()[:3] == flow["records"].tolist()[:3]
    assert source["records"].tolist()[3] == flow[flow["step"].str.contains("in_model")]["records"].item()


@WRITTEN
def test_figure_4s_source_data_reproduces_T10_and_T11():
    import pandas as pd

    source = _read(_MANUSCRIPT / "Figure_4_secondary_and_safety_source_data.csv")
    both = pd.concat([_read(_TABLES / "T10_secondary.csv"), _read(_TABLES / "T11_safety.csv")],
                     ignore_index=True)
    unaugmented = source[source["mark"] != manuscript.Mark.MODEL_ASSISTED.name]

    assert unaugmented["risk difference"].tolist() == both["RD"].tolist()


@WRITTEN
@pytest.mark.parametrize("mid", sorted(manuscript.MANUSCRIPT_IDS))
def test_every_exhibit_landed_with_its_clauses(mid):
    stem = manuscript._stem(mid)
    csv = _MANUSCRIPT / (f"{stem}{manuscript.SOURCE_DATA_SUFFIX}" if mid.startswith("Figure")
                         else f"{stem}.csv")
    assert csv.exists(), csv
    if mid.startswith("Figure"):
        assert (_MANUSCRIPT / f"{stem}.png").exists()

    clauses = [line for line in csv.read_text(manuscript.CSV_ENCODING).splitlines()
               if line.startswith("# ")]
    assert clauses, f"{mid} landed without a caption"
    assert clauses[0].startswith(f"# {mid.replace('_', ' ')}."), clauses[0]
    assert len(clauses) == len(set(clauses)), f"{mid} repeats a clause"


@WRITTEN
def test_the_slide_deck_carries_the_figures_counts_and_not_typed_ones():
    """The PowerPoint diagram is a second rendering of Figure 1, not a second set of numbers.

    Its boxes are editable, which is the point and also the risk: a co-author can retype a count in
    PowerPoint and nothing would notice. What this pins is the file the notebook writes — every box
    text carries the count the run produced, and there is one exclusion annotation per arrow.
    """
    import json

    deck = (_MANUSCRIPT / "Figure_1_cohort_flow.md").read_text("utf-8")
    spec = json.loads(deck.split("- diagram: ", 1)[1].strip())
    source = _read(_MANUSCRIPT / "Figure_1_cohort_flow_source_data.csv")

    assert len(spec["nodes"]) == len(source)
    assert len(spec["labels"]) == len(spec["edges"]) == len(source) - 1
    for node, row in zip(spec["nodes"], source.itertuples()):
        assert f"n = {row.records}   ({row._3} " in node["text"], node["text"]


@WRITTEN
def test_the_legends_file_carries_every_figures_own_caption():
    """One legend page, and every line of it read back off the figure it belongs to.

    A journal asks for the legends separately from the figures, and a legend retyped onto that page is
    the one piece of manuscript text with nothing checking it against the exhibit.
    """
    text = (_MANUSCRIPT / manuscript.LEGENDS).read_text("utf-8")
    figures = [mid for mid in manuscript.MANUSCRIPT_IDS if mid.startswith("Figure")]

    for mid in figures:
        csv = _MANUSCRIPT / f"{manuscript._stem(mid)}{manuscript.SOURCE_DATA_SUFFIX}"
        clauses = [line[2:] for line in csv.read_text(manuscript.CSV_ENCODING).splitlines()
                   if line.startswith("# ")]
        assert f"## {mid.replace('_', ' ')} — {manuscript.MANUSCRIPT_IDS[mid]}" in text
        for clause in clauses:
            assert clause in text, (mid, clause)

    assert text.count("## Figure") == len(figures)
    assert not any(line.startswith(manuscript.PROVENANCE_PREFIX) for line in text.splitlines())


@WRITTEN
@pytest.mark.parametrize("name", ["Table_1_baseline_and_weighting.csv",
                                  "Figure_3_primary_mrs_source_data.csv"])
def test_every_manuscript_csv_announces_its_encoding(name):
    """Three bytes, so that an en dash opened in a spreadsheet is an en dash.

    The Markdown files carry no BOM: a renderer would show it and a diff would carry it.
    """
    assert (_MANUSCRIPT / name).read_bytes().startswith(b"\xef\xbb\xbf")
    assert not (_MANUSCRIPT / manuscript.LEGENDS).read_bytes().startswith(b"\xef\xbb\xbf")


@WRITTEN
def test_the_prose_quotes_the_locked_primary_result():
    """The Results sentence and T07 are one number, not two that happen to agree.

    Prose is the one manuscript surface a regeneration does not touch: a figure gets redrawn when the
    analysis moves and a sentence gets left behind. `manuscript_text.ipynb` interpolates rather than
    types, and this is what says so out loud.
    """
    # The prose is re-flowed, so a quoted interval can straddle a line break: compare on one line.
    text = " ".join((_MANUSCRIPT / "methods_and_results.md").read_text("utf-8").split())
    primary = _read(_TABLES / "T07_primary.csv").set_index("quantity")
    arm = _read(_TABLES / "T12_sensitivity_arm.csv")
    arm = arm[arm["row"].str.startswith("full-covariate")].iloc[0]

    odds = primary.loc["common odds ratio exp(beta)"]
    assert f"{float(odds['estimate']):.2f} (95% CI, {float(odds['ci lo']):.2f} to " \
           f"{float(odds['ci hi']):.2f}; P={float(odds['p']):.3f})" in text
    assert f"{float(arm['exp(beta)']):.2f} (95% CI, {float(arm['ci lo']):.2f} to " \
           f"{float(arm['ci hi']):.2f})" in text
    assert f"{float(arm['worst residual SMD']):.3f}" in text

    rd2 = primary.loc["RD_2  P(mRS<=k) difference"]
    assert f"{float(rd2['estimate']):.3f} (95% CI, {float(rd2['ci lo']):.3f} to " \
           f"{float(rd2['ci hi']):.3f})" in text


@WRITTEN
def test_the_prose_still_marks_what_only_the_study_team_can_supply():
    text = (_MANUSCRIPT / "methods_and_results.md").read_text("utf-8")

    assert text.count("<!-- TO COMPLETE") == 2, "the inclusion criteria and the ethics statement"


@WRITTEN
def test_the_prose_quotes_the_locked_transported_results():
    """The [§14] sentences are read off T22 and T16, not recomputed and not retyped.

    Those analyses are not in the `Bundle` — they are the eighteen minutes of Stage 12 no main exhibit
    reads — so the notebook quotes the locked tables that already carry them. This is what says the
    quoting worked.
    """
    text = " ".join((_MANUSCRIPT / "methods_and_results.md").read_text("utf-8").split())
    comparison = _read(_TABLES / "T22_primary_vs_14.csv").set_index("quantity")
    support = _read(_TABLES / "T16_support_box.csv").set_index("quantity")

    for column in ("[§14a] all eligible", "[§14b] policy"):
        row = comparison.loc["mRS 0-2 difference"]
        lo, hi = row[f"{column} ci"].split(" to ")[0], row[f"{column} ci"].split(" to ")[1].split(" (")[0]
        assert f"{float(row[column]):.3f} (95% CI, {float(lo):.3f} to {float(hi):.3f})" in text, column
        assert f"({int(float(comparison.loc['population', column]))} patients)" in text, column

    assert f"of the {int(support.loc['never-IVT patients', 'value'])} patients" in text


@WRITTEN
def test_the_prose_does_not_call_a_transported_analysis_a_sensitivity_analysis():
    """[§14a] is explicit: it is an ATE in another population, NOT the ATO estimated with more centres.

    The two must not be read as one effect with and without extra centres, and the sentence that says
    so is the one a reader most easily loses.
    """
    text = " ".join((_MANUSCRIPT / "methods_and_results.md").read_text("utf-8").split())

    assert "Neither is a sensitivity analysis on the primary estimand" in text
    assert "not a like-for-like comparison" in text
    assert "transported" in text


@WRITTEN
def test_table_2s_weighted_mrs_classes_are_T08s():
    """The descriptive table and the primary figure print one weighted distribution, not two."""
    table = _read(_MANUSCRIPT / "Table_2_outcomes_by_arm.csv").set_index("outcome")
    distribution = _read(_TABLES / "T08_mrs_distribution_by_arm.csv").set_index("mRS")

    for level in distribution.index:
        row = table.loc[f"mRS {level} at 90 days"]
        for arm in ("EVT alone", "bridging"):
            share = float(distribution.loc[level, f"P(Y=j) {arm}"])
            assert row[f"weighted {arm}"] == f"{100 * share:.1f}%", (level, arm)


@WRITTEN
def test_table_2s_weighted_arms_differ_by_the_risk_difference_reported_in_table_S2():
    """A share is a description and a risk difference is an estimate; they are the same two numbers.

    The weighted arms are rounded to a tenth of a percentage point, so the two roundings allow a fifth
    of a point of slack and no more.
    """
    import pandas as pd

    table = _read(_MANUSCRIPT / "Table_2_outcomes_by_arm.csv").set_index("outcome")
    both = pd.concat([_read(_TABLES / "T10_secondary.csv"), _read(_TABLES / "T11_safety.csv")],
                     ignore_index=True).set_index("outcome")

    for outcome, rd in both["RD"].items():
        row = next(r for label, r in table.iterrows() if label.startswith(outcome))
        difference = float(row["weighted bridging"].rstrip("%")) - float(row["weighted EVT alone"].rstrip("%"))
        assert abs(difference - 100 * float(rd)) < 0.2, outcome


@WRITTEN
def test_table_2_carries_the_crude_rate_warning():
    """[§16b]: any table of raw event rates by arm says they are not unadjusted estimates."""
    text = (_MANUSCRIPT / "Table_2_outcomes_by_arm.csv").read_text(manuscript.CSV_ENCODING)

    assert report.DESCRIPTIVE_CRUDE in text
    assert report.ATO_DESCRIPTION in text


def test_the_model_results_are_supplementary_and_the_main_tables_are_not():
    """The main text carries who was compared and what happened; the estimates are the supplement.

    An `S` id is not decoration: it is how `_stem` names the file, how the caption opens, and how the
    supplement document finds the two tables it compiles.
    """
    main = [mid for mid in manuscript.MANUSCRIPT_IDS if not mid.startswith("Table_S")]
    supplementary = [mid for mid in manuscript.MANUSCRIPT_IDS if mid.startswith("Table_S")]

    assert main == ["Figure_1", "Figure_2", "Figure_3", "Figure_4", "Table_1", "Table_2"]
    assert supplementary == ["Table_S1", "Table_S2"]
    assert manuscript._stem("Table_S1").startswith("Table_S1_")


@WRITTEN
def test_the_supplement_compiles_both_estimate_tables_and_the_transported_ones():
    """Gated on the untracked notebook having run: the document is a deliverable, not an exhibit."""
    document = _MANUSCRIPT / "supplement_tables.md"
    if not document.exists():
        pytest.skip("the supplement notebook is untracked and has not been run here")

    text = document.read_text("utf-8")
    for heading in ("## Table S1.", "## Table S2.", "## Table S3."):
        assert heading in text, heading
    assert report.POPULATIONS in text
    assert report.TRANSPORTED in text
    assert (_MANUSCRIPT / "supplement_tables.docx").exists()
