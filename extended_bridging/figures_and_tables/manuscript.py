"""The manuscript exhibits — four figures, two tables and their supplement, drawn from one locked run.

`report.py` writes the analysis outputs: 22 tables and 5 figures under `OUT/tables` and `OUT/figures`,
in `OUTPUT_IDS` order, discharging [§16]'s STROBE/RECORD checklist. Those are outputs of the *analysis*.
A manuscript needs a different order, because it makes an argument:

::

    Figure 1   who was compared            cohort flow
    Figure 2   were they comparable        balance (Love) + propensity overlap by centre
    Figure 3   what happened               weighted mRS distribution + the six cumulative RD_k
    Figure 4   what else happened          secondary and safety risk differences

    Table 1    baseline and weighting      by arm, unweighted and weighted, with both SMD columns
    Table 2    outcomes by arm             mRS classes and binary outcomes, observed and weighted

    Table S1   the primary result          common OR + RD_k, then [§13]'s full-covariate arm
    Table S2   secondary and safety        RD, model-assisted RD, OR, raw and adjusted p

Everything else stays where it is and becomes supplementary material.

**This module computes nothing.** It loads result objects, formats them and draws them. No estimator
runs here, no interval is recomputed, no number is rounded by hand -- `report.fmt` is the one
formatter and `report`'s clause functions are the one source of [§16]'s sentences.

::

    config.py, report.py, the Stage 1-11 modules
                       │                      (this module's neighbour below is report.py;
                       ▼                       it reads result objects and calls report.*)
    ┌──────────────────────────────────────────────────────────────────────────────┐
    │  manuscript.py                                                               │
    │                                                                              │
    │   load(rebuild)          the cached Bundle — Stages 1-11, built once          │
    │     ├─ Provenance        what the numbers depend on; the cache key            │
    │     └─ Bundle            the nine result objects the seven exhibits read      │
    │                                                                              │
    │   stacked_mrs / forest / love / overlap / flow      the drawing primitives    │
    │   primary_clauses / provenance_clause               the [§16] sentences       │
    │   save_figure / save_table                          OUT/manuscript, PNG + CSV │
    └──────────────────────────────────────────────────────────────────────────────┘
                       │
                       ▼
    figures_and_tables/*.ipynb    one notebook per exhibit; each composes, none computes

Section references in brackets are to `../statistical_analysis_plan.md`. The specification for this
module is `spec.md`, beside it.
"""

from __future__ import annotations

import hashlib
import pickle
import sys
from dataclasses import dataclass, fields
from enum import Enum
from pathlib import Path
from typing import Final, Sequence

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.legend_handler import HandlerTuple

# The shipped modules are flat one directory up; this one is deliberately not among them.
_ROOT: Final[Path] = Path(__file__).resolve().parent.parent
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))

import config as C  # noqa: E402
import data  # noqa: E402
import report  # noqa: E402

# --- where the exhibits land ----------------------------------------------------------------------
#
# A SECOND family, and not an extension of `C.OUTPUT_IDS`. `report.write` iterates that dict and would
# raise KeyError on an id it has no builder for, and [§16]'s checklist would claim a STROBE item was
# discharged by a file that stage never wrote. Ids are file-name stems; slugs are names, never titles.

MANUSCRIPT_IDS: Final[dict[str, str]] = {
    "Figure_1": "cohort_flow",
    "Figure_2": "design_diagnostics",
    "Figure_3": "primary_mrs",
    "Figure_4": "secondary_and_safety",
    "Table_1": "baseline_and_weighting",
    "Table_2": "outcomes_by_arm",
    # The model results are supplementary: the main text carries who was compared and what happened to
    # them, and the estimates live in the supplement the `S` ids compile into.
    "Table_S1": "primary_and_sensitivity",
    "Table_S2": "secondary_and_safety",
}

OUT: Final[Path] = C.OUT / "manuscript"
CACHE: Final[Path] = OUT / ".cache"
SOURCE_DATA_SUFFIX: Final[str] = "_source_data.csv"


class ExhibitError(Exception):
    """An id that is not a manuscript exhibit, or a source-data frame that does not match its figure."""


def _stem(mid: str) -> str:
    if mid not in MANUSCRIPT_IDS:
        raise ExhibitError(f"{mid!r} is not a manuscript exhibit: {sorted(MANUSCRIPT_IDS)}")
    return f"{mid}_{MANUSCRIPT_IDS[mid]}"


# --- what the numbers depend on --------------------------------------------------------------------

SCHEMA: Final[int] = 1        # bump when `Bundle`'s fields change; invalidates every cache


class Rebuild(Enum):
    """Whether `load` may use a cache it finds."""

    IF_STALE = "if stale"     # reuse a cache whose provenance matches; otherwise rebuild
    ALWAYS = "always"         # rebuild and overwrite, whatever is on disk


@dataclass(frozen=True)
class Provenance:
    """Everything that can move a number, as one record. It is both the cache key and its receipt.

    The interpreter and the two library versions are in here because `pyproject.toml` says why they
    are pinned: numpy's NEP 50 promotion and pandas 3.0's string inference each move results, so a
    cache built under one and read under another is a cache of different numbers. `modules` is a
    digest of the shipped source, so editing an estimator invalidates the cache without anyone having
    to remember to. Editing this module or a notebook does not -- neither produces a number.
    """

    schema: int
    data_sha256: str
    modules: str
    seed: int
    n_boot: int
    ci_level: float
    percentile_method: str
    boot_stratum: str
    python: str
    numpy: str
    pandas: str

    @staticmethod
    def current() -> Provenance:
        return Provenance(
            SCHEMA, C.DATA_SHA256, _module_digest(), C.SEED, C.N_BOOT, C.CI_LEVEL,
            C.PERCENTILE_METHOD, C.BOOT_STRATUM, ".".join(str(v) for v in sys.version_info[:2]),
            np.__version__, pd.__version__)

    def key(self) -> str:
        parts = "|".join(f"{f.name}={getattr(self, f.name)}" for f in fields(self))
        return hashlib.sha256(parts.encode("utf-8")).hexdigest()[:12]

    def differences(self, other: Provenance) -> tuple[str, ...]:
        return tuple(f.name for f in fields(self)
                     if getattr(self, f.name) != getattr(other, f.name))


def _module_digest() -> str:
    """One digest over the shipped modules -- the files that turn the workbook into the estimates."""
    h = hashlib.sha256()
    for path in sorted(_ROOT.glob("*.py")):
        h.update(path.name.encode("utf-8"))
        h.update(path.read_bytes())
    return h.hexdigest()[:16]


# --- the cached run ---------------------------------------------------------------------------------

@dataclass(frozen=True)
class Bundle:
    """Stages 1-11: every result object the seven exhibits read, and no more.

    **Not a `report.Run`, deliberately.** A `Run` carries the [§14a] and [§14b] objects too, and its
    V4 guard exists so that a partial one cannot reach `report.write` and print T22 without its
    comparators. Filling those seven fields with placeholders to satisfy the guard would produce an
    object that lies about what it holds; leaving them `None` is what the guard forbids. So this is a
    different record with the SAME field names, which is what lets `report`'s row builders read it
    unchanged -- none of `_t01`, `_t02`, `_t05`, `_t07`, `_t12` or `_binary_rows` touches a Stage
    12/13 field.

    **Nothing unpickled may re-enter an estimator.** `C.Specification` is compared by identity in four
    shipped call sites, and `pickle` restores a copy, not the singleton -- so a `ps` off this record
    fed back to `sensitivity.full_covariate` would raise on a specification that is in fact correct.
    The bundle is input to formatting and drawing, and to nothing else.
    """

    provenance: Provenance
    source: data.Source
    audit: data.Audit
    df: pd.DataFrame
    cohort: pd.DataFrame
    ps: "report.propensity.Propensity"
    balance: "report.balance.Balance"
    primary: "report.outcome.Primary"
    secondary: "report.outcome.Secondary"
    boot: "report.bootstrap.Bootstrap"
    multiplicity: "report.sensitivity.Multiplicity"
    e_value: "report.sensitivity.EValue"
    subgroups: "report.sensitivity.Subgroups"
    arm: "report.sensitivity.Arm"

    def __post_init__(self) -> None:
        missing = [f.name for f in fields(self) if getattr(self, f.name) is None]
        if missing:
            raise C.SchemaError(f"Bundle is missing {missing}: a stage was skipped.")

    @staticmethod
    def build(source: data.Source = data.WORKBOOK) -> Bundle:
        """Stages 1-11 in `report.py`'s canonical order, which is where that order is written."""
        classified, df, audit, ps, bal, est, sec, boot = report._stages_1_10(source)
        mult, ev, sub, arm = report._stage_11(df, ps, bal, est, sec, boot, audit)
        return Bundle(Provenance.current(), source, audit, classified, df, ps, bal, est, sec, boot,
                      mult, ev, sub, arm)


def load(rebuild: Rebuild = Rebuild.IF_STALE, source: data.Source = data.WORKBOOK) -> Bundle:
    """The bundle, from cache when its provenance still holds, otherwise rebuilt and cached.

    Building costs a few minutes -- Stage 10's 2000 replicates, then Stage 11's subgroup and
    full-covariate refits. Seven notebooks each rebuilding it would cost that seven times over and,
    worse, would leave seven results that are only equal because the seed is fixed. One cache is one
    run, and every exhibit is a view of it.
    """
    provenance = Provenance.current()
    path = CACHE / f"bundle_{provenance.key()}.pkl"

    if rebuild is Rebuild.IF_STALE and path.exists():
        bundle = pickle.loads(path.read_bytes())
        bundle.__post_init__()                       # pickle restores state without running it
        moved = provenance.differences(bundle.provenance)
        if moved:
            raise C.SchemaError(f"cache {path.name} was written under different {moved}")
        return bundle

    bundle = Bundle.build(source)
    CACHE.mkdir(parents=True, exist_ok=True)
    path.write_bytes(pickle.dumps(bundle))
    return bundle


# --- house style -------------------------------------------------------------------------------------
#
# One look for the four figures. These are manuscript geometry and they live here rather than in
# `report.py`, whose figures are byte-locked by `C.SVG_RC` so that two runs produce identical bytes.
# Restyling a manuscript figure during review must not be able to move a locked analysis output.

DPI: Final[int] = 300
ARM_COLOUR: Final[dict[int, str]] = {0: "#8c8c8c", 1: "#1f6f8b"}      # EVT alone, bridging
ACCENT: Final[str] = "#b5443a"                                        # thresholds, nulls, exceedance
RULE: Final[str] = "#4d4d4d"
MRS_CMAP: Final[str] = "RdYlGn_r"

RC: Final[dict[str, object]] = {
    "figure.dpi": 110, "savefig.dpi": DPI, "savefig.bbox": "tight",
    "font.size": 8, "axes.titlesize": 9, "axes.labelsize": 8,
    "xtick.labelsize": 7, "ytick.labelsize": 7, "legend.fontsize": 7,
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.grid": False, "figure.facecolor": "white", "savefig.facecolor": "white"}

COLUMN: Final[tuple[float, float]] = (3.4, 2.6)      # single journal column, inches
PAGE: Final[tuple[float, float]] = (7.0, 4.0)        # full width


def use_house_style() -> None:
    """Apply `RC` to the session. Call it once per notebook, in the import cell."""
    plt.rcParams.update(RC)


# --- what a covariate is called ---------------------------------------------------------------------
#
# Analysis names are column names. Every exhibit and the prose read a covariate's clinical term from
# HERE, so the Love plot's ticks and the Methods' covariate list cannot come to call the same thing two
# things. Sentence case, because prose is the harder of the two renderings to fix up; `label` supplies
# the capital a tick wants. `tests/test_manuscript.py` pins that the set covers `C.BALANCE_SET`.

COVARIATE_LABELS: Final[dict[str, str]] = {
    "age": "age",
    "sex": "sex",
    "prestroke_mrs": "prestroke mRS",
    "nihss_baseline": "baseline NIHSS",
    "onset_type": "onset type",
    "core_ml": "ischaemic core volume",
    "tmax6_ml": "Tmax>6 s volume",
    "atrial_fib": "atrial fibrillation",
    "center": "centre",
    "hypertension": "hypertension",
    "hyperlipidemia": "hyperlipidaemia",
    "diabetes": "diabetes",
    "smoking": "smoking",
    "penumbra_ml": "penumbra volume",
}

LEVEL_LABELS: Final[dict[str, str]] = {"wake_up": "wake-up"}     # a level whose name is not its word


def label(covariate: str) -> str:
    """A balance row's name as a reader meets it: `center = HUG` becomes `Centre: HUG`."""
    name, _, level = covariate.partition(" = ")
    text = COVARIATE_LABELS[name]
    if level:
        text = f"{text}: {LEVEL_LABELS.get(level, level)}"
    return text[0].upper() + text[1:]


class Mark(Enum):
    """What a forest row is, which is what decides how its point is drawn.

    `key` is how the mark reads in prose, and it is a field for the reason the marker itself is: the
    legend has to name what the reader sees, and a marker changed here while the legend still says
    "hollow square" is a legend that lies. Nothing is drawn onto the figure to explain a mark.
    """

    ESTIMATE = ("o", True, "filled circles")          # the estimate being reported
    DESCRIPTIVE = ("s", False, "hollow squares")      # [§10, §13] descriptive — estimation, not testing
    MODEL_ASSISTED = ("^", False, "triangles")        # augmented, read beside the unaugmented row

    def __init__(self, marker: str, filled: bool, key: str) -> None:
        self.marker = marker
        self.filled = filled
        self.key = key


# --- the drawing primitives ---------------------------------------------------------------------------

def stacked_mrs(ax, rows: Sequence[tuple[str, Sequence[float]]], title: str = "") -> None:
    """Horizontal stacked bars of P(mRS = j), one bar per row, mRS 0 on the left."""
    cmap = plt.get_cmap(MRS_CMAP, len(C.MRS_LEVELS))
    left = np.zeros(len(rows))

    for j in C.MRS_LEVELS:
        share = np.array([shares[j] for _, shares in rows])
        ax.barh(range(len(rows)), share, left=left, color=cmap(j), edgecolor="white",
                linewidth=0.6, label=f"mRS {j}")
        for i, (value, start) in enumerate(zip(share, left)):
            if value > 0.05:
                ax.text(start + value / 2, i, str(j), ha="center", va="center", fontsize=6,
                        color="white" if j >= 4 else "black")
        left += share

    ax.set_yticks(range(len(rows)), [label for label, _ in rows])
    ax.set_xlim(0, 1)
    ax.set_xlabel("share of the weighted arm")
    ax.invert_yaxis()
    if title:
        ax.set_title(title, loc="left")


def forest(ax, rows: Sequence[tuple[str, float, float, float, Mark]], xlabel: str,
           null: float, log: bool = False) -> None:
    """One row per estimate: the point, its interval, and a reference line at the null."""
    for i, (_, point, lo, hi, mark) in enumerate(rows):
        ax.plot([lo, hi], [i, i], color=RULE, linewidth=1.0, solid_capstyle="butt")
        ax.plot([point], [i], mark.marker, color=RULE, markersize=4.5,
                markerfacecolor=RULE if mark.filled else "white")

    ax.axvline(null, color=ACCENT, linestyle=":", linewidth=1.0)
    ax.set_yticks(range(len(rows)), [label for label, *_ in rows])
    ax.set_xlabel(xlabel)
    ax.invert_yaxis()
    if log:
        ax.set_xscale("log")


def love(ax, covariates: Sequence["report.balance.CovariateBalance"],
         threshold: float = C.SMD_THRESHOLD) -> tuple[list, list[str]]:
    """Absolute SMD before and after weighting, one row per covariate, with [§9]'s threshold drawn.

    Rows are `balance.CovariateBalance`. They are ordered by residual imbalance, largest first, so the
    covariates [§9]'s amendment requires named beside the estimate are the ones the eye lands on.
    """
    ordered = sorted(covariates, key=lambda r: (np.isnan(r.weighted), -abs(r.weighted)))
    y = range(len(ordered))

    for i, row in enumerate(ordered):
        ax.plot([abs(row.unweighted), abs(row.weighted)], [i, i], color="#c8c8c8", linewidth=0.8,
                zorder=1)
    before = ax.scatter([abs(r.unweighted) for r in ordered], y, s=14, facecolor="white",
                        edgecolor=RULE, linewidth=0.8, zorder=2)

    # The weighted points are TWO scatters, not one scatter with a colour array, because the colour
    # carries a fact — red is at or above [§9]'s threshold — and a legend built from a colour array
    # shows whichever colour happened to come first. Two artists let the key show both.
    over = [i for i, r in enumerate(ordered) if abs(r.weighted) >= threshold]
    within = [i for i, r in enumerate(ordered) if i not in over]
    balanced = ax.scatter([abs(ordered[i].weighted) for i in within], within, s=16,
                          color=ARM_COLOUR[1], zorder=3)
    exceeds = ax.scatter([abs(ordered[i].weighted) for i in over], over, s=16, color=ACCENT, zorder=3)

    ax.axvline(threshold, color=ACCENT, linestyle=":", linewidth=1.0)
    ax.set_yticks(list(y), [label(r.covariate) for r in ordered])
    # The threshold is the dotted line and the caption says what it is; naming it again on the axis
    # makes the axis about the rule rather than about the quantity.
    ax.set_xlabel("|standardised mean difference|")
    ax.invert_yaxis()
    return ([before, (balanced, exceeds)], ["unweighted", "overlap-weighted"])


def overlap(ax, bundle: Bundle, centre: str, bins: int = 20) -> None:
    """Propensity scores at one centre, one histogram per arm. The panel is headed by the centre.

    The arm counts are NOT in the heading. They are the substance of this panel and they belong where
    a reader can take the exact number from — the caption and the source data — rather than crowding a
    heading whose job is to say which centre this is.
    """
    at = (bundle.cohort["center"] == centre) & bundle.ps.in_model
    edges = np.linspace(0, 1, bins + 1)

    for arm in sorted(C.TREATMENT_LABELS):
        e = bundle.ps.e[at & (bundle.cohort[C.TREATMENT] == arm)].dropna()
        ax.hist(e, bins=edges, alpha=0.65, color=ARM_COLOUR[arm], label=C.TREATMENT_LABELS[arm])

    ax.set_title(centre, loc="left")
    ax.set_xlim(0, 1)
    ax.set_xlabel("propensity score e(X)")


# A legend entry whose handle is a tuple of artists draws all of them, side by side.
PAIRED_HANDLER: Final[dict] = {tuple: HandlerTuple(ndivide=None)}


def panel_letters(fig, axes: dict, dy: float = 0.05) -> None:
    """The bold panel letter above each panel's top-left corner, clear of that panel's own heading."""
    for letter, ax in axes.items():
        box = ax.get_position()
        fig.text(box.x0, box.y1 + dy, letter, fontweight="bold", fontsize=10, va="bottom", ha="left")


def flow(ax, steps: Sequence[tuple[str, str]], drops: Sequence[str]) -> None:
    """A STROBE flow: one box per step, one right-hand annotation per exclusion between boxes."""
    ax.set_xlim(0, 11)
    ax.set_ylim(0, len(steps) * 2)
    ax.axis("off")

    for i, (label, count) in enumerate(steps):
        y = (len(steps) - i) * 2 - 1.0
        ax.add_patch(plt.Rectangle((0.2, y - 0.62), 5.4, 1.24, facecolor="#f2f4f5",
                                   edgecolor=RULE, linewidth=0.8))
        # Label above, count below: one line each, because a cohort label is long and a count is not.
        ax.text(0.45, y + 0.26, label, va="center", ha="left", fontsize=7.5)
        ax.text(0.45, y - 0.28, count, va="center", ha="left", fontsize=7.5, fontweight="bold")

        if i == len(steps) - 1:
            continue
        ax.annotate("", xy=(2.9, y - 1.40), xytext=(2.9, y - 0.68),
                    arrowprops={"arrowstyle": "-|>", "color": RULE, "linewidth": 0.8})
        if i < len(drops):
            ax.text(5.85, y - 1.04, drops[i], va="center", ha="left", fontsize=7, color=ACCENT)


# --- the [§16] sentences -------------------------------------------------------------------------------
#
# Every one of these comes off `report`, which computes it from the result object. A manuscript
# exhibit that composed its own would be a second definition of what the report says, and the failure
# would be silent: a sentence that stops tracking the run still reads correctly.

def primary_clauses(bundle: Bundle) -> tuple[str, ...]:
    """What [§16] requires wherever the common odds ratio appears."""
    return tuple(clause for clause in (
        report.CONSTANT_SHIFT,
        report.direction_clause(bundle.primary.rd),
        report.exceedance_clause(bundle.balance),
        report.interval_note(bundle.boot)) if clause)


PROVENANCE_PREFIX: Final[str] = "Drawn from the run pinned by"


def provenance_clause(bundle: Bundle) -> str:
    """Which run and which drawing library produced the exhibit.

    matplotlib carries no upper bound in `pyproject.toml`, because figure APPEARANCE has none of the
    reproducibility guarantees the numbers have. Recording the version is what makes a figure that
    redraws differently traceable to that rather than to the analysis.
    """
    p = bundle.provenance
    return (f"{PROVENANCE_PREFIX} DATA_SHA256 {p.data_sha256[:12]}, seed {p.seed}, "
            f"{p.n_boot} replicates; matplotlib {plt.matplotlib.__version__}.")


# --- writing -------------------------------------------------------------------------------------------

# A CSV is opened in a spreadsheet, and a spreadsheet does not sniff the encoding: Excel reads a
# BOM-less file as the system code page and turns the en dash in "74–85", the em dash standing for a
# missing cell and the section sign in "[§6]" into mojibake. The BOM is three bytes that make the file
# say what it is. It goes on the CSVs only — pandas strips it transparently — and never on the
# Markdown, where a renderer would show it and a diff would carry it.
CSV_ENCODING: Final[str] = "utf-8-sig"


def _write(stem: str, suffix: str, content: bytes) -> Path:
    OUT.mkdir(parents=True, exist_ok=True)
    path = OUT / f"{stem}{suffix}"
    path.write_bytes(content)
    return path


def _clauses(clauses: Sequence[str]) -> list[str]:
    """Non-empty, in order, each once.

    Two families of binary outcomes carry the same interval note, and an exhibit that shows both
    would otherwise print it twice — which reads as two different notes about two different
    bootstraps. Dropping the repeat is a presentation decision and belongs here, not in `report.py`,
    where each table carries its own clauses and there is nothing to deduplicate against.
    """
    return list(dict.fromkeys(clause for clause in clauses if clause))


def _rows_of(frame: pd.DataFrame) -> report.Rows:
    """A frame as `report.Rows` -- every cell a string, formatted by the one formatter."""
    header = tuple(str(c) for c in frame.columns)
    body = tuple(tuple(report.fmt(v) if isinstance(v, (int, float, np.floating)) else str(v)
                       for v in row) for row in frame.itertuples(index=False))
    return (header, *body)


def save_table(mid: str, frame: pd.DataFrame, clauses: Sequence[str]) -> Path:
    """`Table_n_<slug>.csv`, rendered by `report.render_csv` so both families share one convention."""
    content = report.render_csv(_rows_of(frame), _clauses(clauses))
    return _write(_stem(mid), ".csv", content.encode(CSV_ENCODING))


def save_figure(mid: str, fig, source_data: pd.DataFrame, clauses: Sequence[str]) -> Path:
    """`Figure_n_<slug>.png` and its source data, which is where the caption lives.

    A PNG carries no text, and [§16]'s statements are required *beside the estimate*. So the numbers
    the figure plots are written as a CSV with the statements beneath them, in `report.render_csv`'s
    order -- the same file the journal asks for as source data. The caption belongs first in
    `clauses`, so it heads that block.
    """
    stem = _stem(mid)
    png = OUT / f"{stem}.png"
    OUT.mkdir(parents=True, exist_ok=True)
    fig.savefig(png, dpi=DPI, bbox_inches="tight", facecolor="white")
    plt.close(fig)

    content = report.render_csv(_rows_of(source_data), _clauses(clauses))
    _write(stem, SOURCE_DATA_SUFFIX, content.encode(CSV_ENCODING))
    return png


def save_slide(mid: str, deck: str) -> Path:
    """`Figure_n_<slug>.md` — the exhibit as a slide deck source, for the powerpoint skill.

    A matplotlib figure is a picture of a diagram; a co-author cannot move a box in it. A cohort flow
    is the one exhibit here that is pure geometry and no data-derived shape, so it is also worth having
    as native, editable shapes. What this writes is the SPEC, with every count taken from the run; the
    rendering is a separate command and a separate tool, which is why this layer does not shell out to
    it. `spec.md` §7 carries the command.
    """
    return _write(_stem(mid), ".md", deck.encode("utf-8"))


LEGENDS: Final[str] = "figure_legends.md"


def save_legends() -> Path:
    """Every figure legend in one file, read back out of the figures themselves.

    Journals ask for the legends as a page of their own, and a legend retyped there is a legend that
    stops tracking the figure. So this reads the `# ` clause block each figure already wrote beside its
    PNG: the caption first, then the [§16] statements that must be printed with it. It rebuilds the
    whole file from whatever is on disk, so running one notebook leaves it complete rather than partial.

    The provenance sentence is kept as an HTML comment — it says which run drew the figure, which a
    reader of the manuscript does not need and a co-author checking the figure does.
    """
    blocks = [f"# Figure legends\n",
              "Generated by `manuscript.save_legends()` from the clause block beside each figure — "
              "edit the notebook that writes the figure, never this file.\n"]

    for mid, slug in MANUSCRIPT_IDS.items():
        if not mid.startswith("Figure"):
            continue
        csv = OUT / f"{_stem(mid)}{SOURCE_DATA_SUFFIX}"
        if not csv.exists():
            continue

        clauses = [line[2:] for line in csv.read_text(encoding=CSV_ENCODING).splitlines()
                   if line.startswith("# ")]
        if not clauses:
            continue

        printed = [c for c in clauses[1:] if not c.startswith(PROVENANCE_PREFIX)]
        hidden = [c for c in clauses if c.startswith(PROVENANCE_PREFIX)]
        blocks.append(f"## {mid.replace('_', ' ')} — {slug}\n")
        blocks.append(clauses[0] + "\n")
        if printed:
            blocks.append(" ".join(printed) + "\n")
        blocks.extend(f"<!-- {c} -->\n" for c in hidden)

    return _write(LEGENDS.removesuffix(".md"), ".md", "\n".join(blocks).encode("utf-8"))


def save_prose(stem: str, text: str) -> Path:
    """A section of manuscript prose, with every number interpolated from the run.

    Prose is where a result is most likely to be wrong and least likely to be caught: a figure is
    regenerated, a sentence is retyped. So the Methods and Results are written as a template the
    notebook fills from the same bundle the exhibits are drawn from, and no number in them is keyed
    by hand. What the template cannot supply — the source cohort's inclusion criteria, the ethics
    statement — is left as a marked placeholder rather than invented.
    """
    return _write(stem, ".md", text.encode("utf-8"))


def show(path: Path):
    """The written PNG, displayed in the notebook -- the artefact itself, not a second rendering.

    `report.py` sets the Agg backend at import [report.py §matplotlib], so whether an inline figure
    appears depends on an import order the reader cannot see. Displaying the file removes the
    question: what the notebook shows is the bytes that went to `OUT/manuscript`.
    """
    from IPython.display import Image

    return Image(filename=str(path))


def written() -> pd.DataFrame:
    """Every exhibit on disk with its size and digest -- the manuscript family's own manifest."""
    rows = []
    for mid in MANUSCRIPT_IDS:
        for path in sorted(OUT.glob(f"{_stem(mid)}*")):
            content = path.read_bytes()
            rows.append((mid, path.name, len(content), hashlib.sha256(content).hexdigest()[:16]))
    return pd.DataFrame(rows, columns=["exhibit", "file", "bytes", "sha256"])


__all__ = [
    "Bundle", "ExhibitError", "Mark", "Provenance", "Rebuild", "MANUSCRIPT_IDS", "OUT",
    "load", "use_house_style", "stacked_mrs", "forest", "love", "overlap", "flow",
    "COVARIATE_LABELS", "PAIRED_HANDLER", "label", "panel_letters", "primary_clauses", "provenance_clause", "save_table", "save_figure", "save_slide", "save_legends", "save_prose", "show",
    "written"]
