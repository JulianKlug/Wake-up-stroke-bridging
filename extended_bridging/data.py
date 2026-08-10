"""Stage 2 — load and clean [§11].

One reader that turns the private workbook into the analysis frame, and one audit log that makes
every difference between the two visible. Nothing is imputed, nothing is corrected by row position,
and no correction happens that the log does not name the patients for.

Section references in brackets are to ``statistical_analysis_plan.md``. DECISION *n* refers to the
decisions recorded in ``../out/stage0_data_inventory.md``. The specification for this module is
``specs/stage2_load_and_clean.md``; nothing here is invented outside it.

::

    data/…v7…with_abs_contra_indication.xlsx        (gitignored, patient data)
                       │
                       ▼
    ┌──────────────────────────────────────────────────────────────────────┐
    │  STAGE 2 — data.py                                   reads config.py │
    │                                                                      │
    │   load(source)                                                       │
    │     │                                                                │
    │     ├─ audit = Audit(source)     created first: every step below     │
    │     │                            records into this one object        │
    │     ├─ _verify_version(source)   sha256 vs DATA_SHA256               │
    │     │            ▼                              → DataVersionError   │
    │     ├─ _read(source, audit)      READ_DTYPES, NA_VALUES, full parse  │
    │     │            ▼               records provenance/read             │
    │     ├─ _apply_contract(raw, source, audit)  assert_column_contract,  │
    │     │            ▼               rows, RENAME, drop → SchemaError    │
    │     │                            records contract/rename_and_drop    │
    │     └─ _pipeline_after_read(df, source, audit)  ← the tests' seam    │
    │            │                                                         │
    │            ├─ _normalise(df, audit)      CENTER_RECODE               │
    │            │            ▼                                            │
    │            ├─ _correct(df, audit)        content-driven only,        │
    │            │            ▼                every case named            │
    │            ├─ _assert_schema(df, source) A1…A9 + A4b, every          │
    │            │            ▼                violation collected         │
    │            └─ _missingness(df, audit)    structural / not recorded / │
    │                         ▼                missing / complete          │
    └──────────────────────────────────────────────────────────────────────┘
         │
         ▼
    (df, Audit)          load() writes NOTHING. The caller decides:
         │                   audit.write()  →  out/logs/stage2_audit_<label>.md
         │                                     (gitignored: names patients)
         ├──────────────┬──────────────┬──────────────┐
         ▼              ▼              ▼              ▼
     Stage 3        Stage 4        Stage 5      Stage 14
     derive         eligibility    cohort       reporting (denominators)

``load()`` has no filesystem side effect, and the log path carries the source's label. Both follow
from one fact: the acceptance tests run the whole pipeline against ``FIXTURE`` on every ``pytest``
invocation, including on machines that have ``data/``. A ``load()`` that wrote a fixed path would
overwrite the workbook's audit log with the fixture's two-row log on every test run, silently — the
file is gitignored, so no diff would ever show it.

The ``Audit`` is created first and threaded through, not built at the end: two entries are recorded
inside ``_read`` and ``_apply_contract``, and the log's header is rendered *from the entries* rather
than from a second traversal of the source. One number, one origin.

Every arrow out of this module carries analysis names only. This module is **not** exempt from the
Stage 1 §7 raw-name scan and must never become exempt.
"""
from __future__ import annotations

import hashlib
from collections.abc import Iterable, Sequence
from dataclasses import dataclass, field
from pathlib import Path
from typing import Final

import pandas as pd

import config as C


# --- the source, and why it is a type ----------------------------------------------------------
#
# load() must run against two files: the private workbook, and the Stage 1 fixture, which carries
# the 42 verbatim headers and two invented rows. Without the second, every Stage 2 test needs
# `data/`, and `data/` is gitignored.
#
# The obvious shape, load(path=None, check_hash=True, check_rows=True), is the wrong one: it lets a
# caller switch off a guarantee with a keyword, and a keyword in a call site is invisible from the
# module that declares the guarantee. Instead the file and its guarantees travel together.

@dataclass(frozen=True)
class Source:
    """A file load() may be run against, carrying the guarantees that apply to it."""

    path: Path
    sha256: str | None       # None: no version check. Only the fixture may be None.
    n_records: int | None    # None: no row-count check. Only the fixture may be None.
    label: str               # for the audit log header; never a path


WORKBOOK: Final[Source] = Source(
    C.DATA_XLSX, C.DATA_SHA256, C.N_RECORDS_EXPECTED,
    "v7_july26_with_abs_contra_indication")

FIXTURE: Final[Source] = Source(
    C.FIXTURE_XLSX, None, None, "fixture_schema")

# Pinned against this literal by the acceptance tests, so a third module-level source cannot appear
# without a test failing. A *test-local* Source — over a tmp_path file, or over the fixture's path
# with a deliberately wrong hash — is legal and is what lets the version and shuffle tests run on a
# checkout with no `data/`. What this tuple pins is what data.py declares, not what Source permits.
SOURCES: Final[tuple[Source, ...]] = (WORKBOOK, FIXTURE)


# --- the audit log ------------------------------------------------------------------------------
#
# Roadmap Stage 2 accepts when "the audit log reproduces byte-identically on a re-run, and every
# correction in it names the cases it touched". Both clauses are properties of the *rendering*, so
# both live here rather than in a caller.
#
# Byte identity is not automatic. Every rule below is a way it has been lost before:
#   - no timestamp, no clock time, no run identifier: the log is a function of the data
#   - no absolute path: the header prints `source.path.name` and `source.label`
#   - no pandas or numpy repr: every string in the log is built here
#   - no set iteration: where a loop ranges over columns it ranges over a declared tuple
#   - one float formatter, `_fmt`, with missing rendered as the literal `missing`
#   - identifiers sorted, not in frame order — which is what makes the log invariant to input row
#     order, and so what lets the shuffle test detect a correction written by row position

KINDS: Final[tuple[str, ...]] = (
    "provenance", "contract", "correction", "observation", "structural", "missingness")

# The six section headings of the rendered document, in KINDS order. A section with no entries
# still prints its heading and the single line `_none_`, so a clean workbook's log reads "nothing
# was observed" rather than looking truncated, and the skeleton is fixed for the eye to scan.
_HEADINGS: Final[dict[str, str]] = {
    "provenance":  "Provenance",
    "contract":    "Contract",
    "correction":  "Corrections",
    "observation": "Observations",
    "structural":  "Structural non-applicability",
    "missingness": "Missingness and denominators",
}

# Roadmap Stage 2 demands the names-its-cases rule of corrections. `observation` is held to it too:
# the one observation this stage defines is a standing query with the data owner (§7.3), and a query
# that names no patient cannot be answered, which makes it exactly as useless as an unattributed
# correction. The other four kinds legitimately name none — provenance, contract and missingness
# describe columns rather than patients, and structural describes a whole arm.
_MUST_NAME_CASES: Final[frozenset[str]] = frozenset({"correction", "observation"})


@dataclass(frozen=True)
class AuditEntry:
    """One recorded step. The correction-names-its-cases rule lives here, not in a reviewer's eye."""

    kind: str
    step: str
    n: int
    detail: str
    case_ids: tuple[str, ...] = ()
    table: tuple[tuple[str, ...], ...] | None = None   # header row first

    def __post_init__(self) -> None:
        if self.kind not in KINDS:
            raise ValueError(
                f"{self.kind!r} is not an audit kind. Known kinds: {list(KINDS)}. "
                "An entry whose kind is not declared would render under no section heading.")
        if self.kind in _MUST_NAME_CASES and self.n > 0 and not self.case_ids:
            raise ValueError(
                f"{self.kind} {self.step!r} touched {self.n} record(s) but names no case. "
                "Roadmap Stage 2 requires every correction to name the cases it touched, "
                "and an observation nobody can look up is not a query.")


def _fmt(x: object) -> str:
    """The one float formatter. Six significant figures is exact for every value in this workbook.

    Declared once and used everywhere a number reaches the log, because a second formatter is a
    second way for two runs to disagree. Missing renders as the literal `missing`, never as a pandas
    repr, which is free to change between minor releases.
    """
    if x is None or (not isinstance(x, str) and pd.isna(x)):
        return "missing"
    return f"{float(x):.6g}"


def _md_table(rows: Sequence[Sequence[str]]) -> str:
    """First row is the header; every cell is already a string; columns padded to their widest cell.

    Fifteen lines rather than `tabulate`, deliberately: byte-identical reproduction is an acceptance
    criterion, and a formatting library's column padding, alignment rules and float repr are free to
    change between minor versions. A criterion that rests on a third party fails on an unrelated
    `uv sync`.
    """
    widths = [max(len(row[j]) for row in rows) for j in range(len(rows[0]))]
    body = ["| " + " | ".join(c.ljust(w) for c, w in zip(row, widths)) + " |" for row in rows]
    separator = "| " + " | ".join("---" for _ in widths) + " |"
    return "\n".join([body[0], separator, *body[1:]])


def _audit_path(source: Source) -> Path:
    """Derived from the label, so the fixture's log can never overwrite the workbook's."""
    return C.LOGS / f"stage2_audit_{source.label}.md"


class Audit:
    """The ordered collection of entries, and their renderer."""

    def __init__(self, source: Source) -> None:
        self.source = source
        self.entries: list[AuditEntry] = []

    def record(self, kind: str, step: str, n: int, detail: str,
               case_ids: Iterable[object] = (),
               table: tuple[tuple[str, ...], ...] | None = None) -> None:
        """Append an entry, normalising identifiers to a sorted tuple of strings on the way in."""
        self.entries.append(AuditEntry(
            kind, step, int(n), detail, tuple(sorted(str(c) for c in case_ids)), table))

    def entry(self, kind: str, step: str) -> AuditEntry | None:
        """The first entry with this (kind, step), or None. The header reads its numbers back."""
        return next((e for e in self.entries if e.kind == kind and e.step == step), None)

    def to_markdown(self) -> str:
        return "\n".join([*self._header(), "", *self._sections(), ""])

    def write(self, path: Path | None = None) -> Path:
        """Write the log and return the path. load() never calls this; a caller asks for it."""
        path = _audit_path(self.source) if path is None else path
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(self.to_markdown(), encoding="utf-8")
        return path

    # The header's numbers come from the entries, never from a second traversal of the source.
    # Recomputing them here would give the log two sources for one number, and the one that drifts
    # is always the one nobody tests. This is why the Audit is created at the top of load().
    def _header(self) -> list[str]:
        version = (f"`{self.source.sha256[:8]}…{self.source.sha256[-7:]}`"
                   if self.source.sha256 is not None else f"not pinned ({self.source.label})")
        return [
            "# Stage 2 — data audit",
            "",
            f"Source: `{self.source.path.name}`, sheet `{C.SHEET}`",
            f"Version: {version}",
            f"Read: {self._read_line()}",
            "Generated by `extended_bridging/data.py`. Regenerate rather than edit.",
        ]

    def _read_line(self) -> str:
        read, contract = self.entry("provenance", "read"), self.entry("contract", "rename_and_drop")
        if read is None or contract is None or contract.table is None:
            # The pipeline was driven from a frame that was already read — the acceptance tests do
            # exactly this. There is no read to report, and inventing one would be the second source
            # of truth this method exists to avoid.
            return "not recorded (this audit did not perform the read)"
        n_columns = len(contract.table) - 1          # one row per raw column, header first
        return (f"{read.n} rows x {n_columns} columns; "
                f"{n_columns - contract.n} mapped, {contract.n} dropped")

    def _sections(self) -> list[str]:
        out: list[str] = []
        for kind in KINDS:
            out.append(f"## {_HEADINGS[kind]}")
            out.append("")
            entries = [e for e in self.entries if e.kind == kind]
            if not entries:
                out.extend(["_none_", ""])
                continue
            for entry in entries:
                out.extend(_render_entry(entry))
        return out

    def __len__(self) -> int:
        return len(self.entries)


def _render_entry(entry: AuditEntry) -> list[str]:
    """The pilot's shape, which reads well: a bullet, then the cases, then any table beneath."""
    detail = entry.detail.split("\n")
    lines = [f"- **{entry.step}** (n={entry.n}): {detail[0]}"]
    lines.extend(f"  {line}" if line else "" for line in detail[1:])
    if entry.case_ids:
        lines[-1] = lines[-1] + "  "                       # a markdown hard break, not decoration
        lines.append("  Affected: `" + ", ".join(entry.case_ids) + "`")
    if entry.table is not None:
        lines.append("")
        lines.extend(f"  {row}" for row in _md_table(entry.table).split("\n"))
    lines.append("")
    return lines


# --- the read -----------------------------------------------------------------------------------

def _verify_version(source: Source) -> None:
    """sha256 the file before anything else: "which file is this" beats "column X is missing"."""
    if source.sha256 is None:
        return
    digest = hashlib.sha256(source.path.read_bytes()).hexdigest()
    if digest != source.sha256:
        raise C.DataVersionError(
            f"\n  expected {source.sha256[:8]}…  ({source.label})"
            f"\n  actual   {digest[:8]}…"
            "\n  The workbook changed. Re-run stage0_data_inventory.py, review"
            "\n  ../out/stage0_data_inventory.md, and update DATA_SHA256 in the same commit.")


def _read(source: Source, audit: Audit) -> pd.DataFrame:
    """A full parse, never nrows=0.

    The workbook's final column is headerless and pandas materialises it only once rows are read: a
    header-only read returns 41 columns, a full read 42. `assert_column_contract` has a dedicated
    branch that names this cause, and it can only help if this reader never triggers it deliberately.

    `keep_default_na` is passed explicitly so that it is a decision rather than a default. Setting it
    False — the plausible move for someone hardening the reader — would leave the workbook's literal
    'N/A' cells as strings in Int64 columns and crash the read. That failure is loud, so it is not
    dangerous; the danger is the repair, which would be to widen NA_VALUES blindly.
    """
    raw = pd.read_excel(
        source.path,
        sheet_name=C.SHEET,
        dtype=C.READ_DTYPES,
        na_values=list(C.NA_VALUES),
        keep_default_na=True,
    )
    state = (f"pinned and verified ({source.sha256:.8})" if source.sha256 is not None
             else f"not pinned ({source.label})")
    audit.record(
        "provenance", "read", len(raw),
        f"sheet {C.SHEET!r}; {len(raw)} rows x {raw.shape[1]} columns; version {state}")
    return raw


def _apply_contract(raw: pd.DataFrame, source: Source, audit: Audit) -> pd.DataFrame:
    """Check the columns, check the rows, then rename and drop. The order is the specification.

    The contract entry is recorded *after* both checks pass: a log that records what was dropped
    from a frame that then failed `assert_column_contract` describes a frame that never existed.
    """
    C.assert_column_contract(raw.columns)

    if source.n_records is not None and len(raw) != source.n_records:
        raise C.SchemaError(
            f"{source.path.name} has {len(raw)} rows; {source.label} declares "
            f"{source.n_records}. A filtered sheet or a partial read produces a smaller cohort "
            "that every downstream denominator then reports as fact.")

    # The drop list ranges over COLUMN_CONTRACT, never over DROPPED itself. DROPPED is a frozenset,
    # CPython randomises str hashing per process unless PYTHONHASHSEED is set, and so list(DROPPED)
    # is in a different order in every interpreter. For .drop() alone that is harmless — pandas
    # preserves the frame's column order regardless — but the table below is rendered into the log,
    # and byte-identical reproduction would then be lost *across processes* while still passing a
    # same-process comparison.
    dropped = [key for key in C.COLUMN_CONTRACT if key in C.DROPPED]
    df = raw.rename(columns=C.RENAME).drop(columns=dropped)

    table = (("raw column", "fate"),) + tuple(
        (key, column.name if column.name is not None else "dropped")
        for key, column in C.COLUMN_CONTRACT.items())
    audit.record(
        "contract", "rename_and_drop", len(dropped),
        f"{len(C.RENAME)} columns mapped to analysis names, {len(dropped)} dropped; every dropped "
        "column's reason is in config.COLUMN_CONTRACT",
        table=table)
    return df


# --- normalisation ------------------------------------------------------------------------------

def _normalise(df: pd.DataFrame, audit: Audit) -> pd.DataFrame:
    """Two operations, neither of which changes a value's meaning.

    `center` is recoded through CENTER_RECODE and stays a `string`. It is deliberately not made a
    Categorical here: Stage 5 drops the zero-bridging centre, and a categorical that keeps a dead
    level makes every subsequent groupby(observed=False) resurrect it as an all-missing row.
    Stage 6 builds the categorical from FACTOR_LEVELS at the point of use.

    `case_id` is asserted (A2), not stripped. Stripping it would be a silent correction of the kind
    §7 forbids; asserting it turns a future padded identifier into an error naming the value.
    `contraindication_reason` is not normalised at all: under DECISION 1 the free text is never
    read, so there is nothing to normalise, and normalising it would create the appearance that the
    text feeds a classifier.
    """
    recoded = df["center"].map(C.CENTER_RECODE)
    # An unmapped code keeps its original text rather than becoming <NA>, so that A3 can name the
    # offending value. A code silently mapped to <NA> would drop that centre from every by-centre
    # table and from the [§3] restriction, and the error would name nothing.
    #
    # The cast back to `string` is not cosmetic: `.map` returns object dtype, and [§11] promises
    # Stage 3 a `string` centre column. An object column would still compare and group correctly,
    # so nothing downstream would fail — it would simply stop being the column the contract
    # describes, which is the kind of drift this stage exists to make visible.
    df["center"] = recoded.fillna(df["center"]).astype("string")

    by_label = {label: code for code, label in C.CENTER_RECODE.items()}
    table = (("code", "label", "n"),) + tuple(
        (by_label[label], label, str(int((df["center"] == label).sum())))
        for label in C.CENTER_ORDER)
    audit.record(
        "provenance", "centre_recode", int(recoded.notna().sum()),
        f"CENTER_RECODE applied; {len(C.CENTER_RECODE)} codes mapped to "
        f"{len(C.CENTER_ORDER)} labels; no value altered, a code was given its name",
        table=table)
    return df


# --- corrections --------------------------------------------------------------------------------
#
# Two findings, two different fates. The asymmetry is the whole of §7's argument and is the thing
# most likely to be "tidied" by someone who has not read §7.3:
#
#   penumbra_ml disagrees with tmax6_ml - core_ml      onset_to_groin_min == 999
#   on 1 record (the core volume was copy-pasted)      on 1 record
#          │                                                  │
#          ▼                                                  ▼
#   RECOMPUTED on every record, because [§6] makes     LEFT STANDING. 999 is inside the observed
#   the definition authoritative and the stored        range, the implied IVT-to-groin interval is
#   column merely a copy of it                         44 min against a treated median of 70, and
#          │                                           neighbouring values are equally round. The
#          ▼                                           entire case for "placeholder" is the digits
#   logged as a `correction`, naming the case                 │
#                                                             ▼
#                                                      logged as an `observation`, naming the case;
#                                                      a standing query with the data owner
#
# Every correction is content-driven: a boolean mask over column values, never .iloc[i], never
# .loc[i] with an integer, never .index[…], never .head/.tail. The returned frame keeps the
# workbook's row order and a default RangeIndex, and that index is not an identifier.

_TOL: Final[float] = 1e-6            # absolute; the volumes are millilitres in the tens and hundreds
_GROIN_SENTINEL: Final[int] = 999


def _correct(df: pd.DataFrame, audit: Audit) -> pd.DataFrame:
    """Recompute `penumbra_ml` from its definition; flag the 999 groin time without touching it."""
    computed = df["tmax6_ml"] - df["core_ml"]
    disagree = (df["penumbra_ml"] - computed).abs() > _TOL

    # A record with a missing volume is never counted as disagreeing: (NA - x).abs() > _TOL is
    # false, not missing, so those records fall out without a guard. That is the behaviour we want,
    # and it is what `n` below means — a reader who assumed missing rows counted as disagreements
    # would read v7's n = 1 as evidence of a bug rather than as the answer.
    threshold = C.TARGET_MISMATCH["min_penumbra_ml"]
    crossing = disagree & ((df["penumbra_ml"] > threshold) != (computed > threshold))

    # Applied to *every* record, not only the disagreeing ones. A patch of the disagreeing rows
    # leaves the stored column authoritative wherever it happens to agree, and a second copy-paste
    # that is self-consistent by accident would survive it. Recomputing everywhere is content-driven
    # by construction: there is no mask to get wrong and no row to name. Missingness is preserved
    # because subtraction propagates it.
    stored = df["penumbra_ml"]
    df["penumbra_ml"] = computed

    n = int(disagree.sum())
    detail = (
        "[§6] defines penumbra as tmax6_ml - core_ml, so the definition is authoritative and the "
        f"stored column is a copy of it. Recomputed on all {len(df)} records; {n} disagreed beyond "
        f"{_fmt(_TOL)}.")
    table = None
    if n:
        detail += (
            f" {int(crossing.sum())} record(s) cross TARGET_MISMATCH['min_penumbra_ml'] = "
            f"{_fmt(threshold)}, so [§13] target-mismatch membership changes for them.")
        # The per-record values go in the table, not in the sentence: a sentence naming "stored X vs
        # recomputed Y" reads correctly for the one record v7 has and becomes ambiguous for the
        # second one a corrected workbook brings.
        rows = sorted(
            (str(case_id), _fmt(was), _fmt(now), _fmt(now - was))
            for case_id, was, now in zip(
                df.loc[disagree, "case_id"], stored[disagree], computed[disagree]))
        table = (("case_id", "stored", "recomputed", "difference"),) + tuple(rows)
    audit.record("correction", "penumbra_recomputed", n, detail,
                 df.loc[disagree, "case_id"], table)

    _observe_groin_sentinel(df, audit)
    return df


def _observe_groin_sentinel(df: pd.DataFrame, audit: Audit) -> None:
    """Log the 999 groin times as a standing query. The value stands; setting it missing would delete
    real data on the strength of a digit pattern.

    `pilots/data.py` sets these to missing, on the reasoning that 999 is a placeholder. Stage 2 does
    not carry that correction: the value is inside the observed range, its implied IVT-to-groin
    interval is unremarkable, and neighbouring observed values are equally round. This is the same
    posture DECISION 2 took toward the shipped death flags — a derivation rule where one is
    available, a logged standing query where none is, and never a silent resolution.
    """
    groin, ivt_time = df["onset_to_groin_min"], df["onset_to_ivt_min"]
    sentinel = groin == _GROIN_SENTINEL
    n = int(sentinel.sum())
    if not n:
        return

    treated = df[C.TREATMENT] == 1
    intervals = (groin - ivt_time)[treated & groin.notna() & ivt_time.notna()]
    gaps = ", ".join(_fmt(g) for g in (groin - ivt_time)[sentinel])
    audit.record(
        "observation", "onset_to_groin_999", n,
        f"onset_to_groin_min == {_GROIN_SENTINEL} on {n} record(s). Not corrected (§7.3): the value "
        f"is inside the observed range {_fmt(groin.min())}-{_fmt(groin.max())}, the implied "
        f"IVT-to-groin interval is {gaps} min against a treated-arm median of "
        f"{_fmt(intervals.median())}, and neighbouring observed values are equally round. "
        "Standing query with the data owner.",
        df.loc[sentinel, "case_id"])


# --- schema assertions --------------------------------------------------------------------------
#
# Run after the corrections, once, over the analysis frame. After, not before: penumbra_ml is
# recomputed in full, so asserting the stored column first would be asserting something Stage 2 is
# about to discard.
#
# Every check runs, and one SchemaError reports all the violations. Not fail-fast: this module's job
# is to meet a workbook it has never seen — the data owner will send a v8 — and a first-failure
# exception turns "three columns drifted" into three read-hash-parse-fail cycles, with the operator
# learning one fact per run. Each check guards its own preconditions, so no check depends on another
# having passed and the order below is presentation, not semantics.
#
# There are no bare `assert` statements here, or anywhere in this module. Python strips `assert`
# under -O, so under a flag nobody remembers setting, a module whose entire purpose is to fail
# loudly would succeed silently. The acceptance tests scan this file for ast.Assert.


def _values(series: pd.Series) -> str:
    """Sorted distinct values of an offending selection, for an error message."""
    return ", ".join(sorted({_fmt(v) if not isinstance(v, str) else repr(v) for v in series}))


def _assert_schema(df: pd.DataFrame, source: Source) -> None:
    """A1–A9 plus A4b. Every violation is collected; one SchemaError reports them all."""
    bad: list[str] = []

    # A1 — a post-condition, and unreachable today: nothing between the contract and here adds or
    # removes a row. It is kept anyway, at the cost of one comparison on 126 records, because it is
    # the check that catches a future Stage 2 edit which filters rows. Do not delete it as dead
    # code, and do not expect a test to reach it through load().
    if source.n_records is not None and len(df) != source.n_records:
        bad.append(f"A1  rows: {len(df)} record(s); {source.label} declares {source.n_records}")

    case_id = df["case_id"]
    if int(case_id.isna().sum()):
        bad.append(f"A2  case_id: {int(case_id.isna().sum())} record(s) have no identifier")
    duplicated = case_id[case_id.duplicated(keep=False) & case_id.notna()]
    if len(duplicated):
        bad.append(f"A2  case_id: duplicated identifier(s): {_values(duplicated)} — a duplicate "
                   "double-weights one patient through the propensity fit and the bootstrap")
    padded = case_id[case_id.notna() & (case_id != case_id.str.strip())]
    if len(padded):
        bad.append(f"A2  case_id: identifier(s) carrying surrounding whitespace: {_values(padded)}")

    center = df["center"]
    unmapped = center[~center.isin(tuple(C.CENTER_RECODE.values()))]
    if len(unmapped):
        bad.append(f"A3  center: value(s) outside CENTER_RECODE: {_values(unmapped)}")

    for check, column in (("A4", C.TREATMENT), ("A4b", "ivt_contraindicated")):
        n_missing = int(df[column].isna().sum())
        if n_missing:
            bad.append(f"{check}  {column}: {n_missing} record(s) missing the "
                       + ("treatment" if column == C.TREATMENT else "eligibility classifier"))

    for column in C.BINARY_COLUMNS:
        offending = df[column][df[column].notna() & ~df[column].isin((0, 1))]
        if len(offending):
            bad.append(f"A5  {column}: value(s) outside {{0, 1}}: {_values(offending)}")

    for column, (low, high) in C.PLAUSIBLE_RANGES.items():
        values = df[column][df[column].notna()]
        offending = values[(values < low) | (values > high if high is not None else False)]
        if len(offending):
            bad.append(f"A6  {column}: value(s) outside ({low}, {high}): {_values(offending)}")

    # A7 skips records whose exposure is missing and lets A4 report those: letting A7 compare
    # against a missing exposure produces a second confusing message about the same records.
    known = df[C.TREATMENT].notna()
    has_time = df["onset_to_ivt_min"].notna()
    treated_without = df.loc[known & (df[C.TREATMENT] == 1) & ~has_time, "case_id"]
    control_with = df.loc[known & (df[C.TREATMENT] == 0) & has_time, "case_id"]
    if len(treated_without):
        bad.append(f"A7  onset_to_ivt_min: absent for {len(treated_without)} treated patient(s): "
                   f"{_values(treated_without)} — that is not ordinary missingness [§11]")
    if len(control_with):
        bad.append(f"A7  onset_to_ivt_min: present for {len(control_with)} control(s): "
                   f"{_values(control_with)} — the arm label is wrong")

    both = df["onset_to_ivt_min"].notna() & df["onset_to_groin_min"].notna()
    out_of_order = df.loc[both & (df["onset_to_ivt_min"] > df["onset_to_groin_min"]), "case_id"]
    if len(out_of_order):
        bad.append(f"A8  onset_to_ivt_min: after groin puncture on {len(out_of_order)} record(s): "
                   f"{_values(out_of_order)} — the record is not bridging at all")

    reason = df["contraindication_reason"]
    blank = reason[reason.notna() & (reason.str.strip() == "")]
    if len(blank):
        bad.append(f"A9  contraindication_reason: {len(blank)} whitespace-only cell(s) — a blank "
                   "reads as 'a reason was recorded' and moves a patient from indeterminate to "
                   "eligible [DECISION 1, §3]")

    if bad:
        raise C.SchemaError(
            "\n  " + "\n  ".join(bad)
            + f"\n  {len(bad)} assertion(s) failed against {len(df)} records. Every one is a "
            "contract\n  violation; none of them is a reason to relax the contract.")


# --- missingness and denominators [§11] ----------------------------------------------------------
#
# Absence has three kinds, and a single "n missing" column misrepresents all three:
#
#   structural      a key of STRUCTURALLY_NON_APPLICABLE. onset_to_ivt_min is absent for all 87
#                   controls because they were never given IVT. That is structure, not data loss.
#   not recorded    a key of INFORMATIVE_ABSENCE. contraindication_reason is absent for 82 of 126
#                   records, and that absence *is* the signal: under DECISION 1 it is the one bit
#                   the column contributes, the thing that separates eligible from indeterminate.
#   missing         anything else with a non-zero absence count — genuine data loss
#   complete        no absences
#
# A reader who does not see the split reads 82/126 on contraindication_reason as two-thirds
# unusable, when Stage 4 uses every one of those 82. Stage 2 labels them; nothing imputes them.
#
# The per-centre columns are there because the pattern is centre-driven, not random: the
# contraindication reason was recorded for controls at HUG and USZ and nowhere else. A pooled count
# hides that shape entirely, and it is the shape that makes the indeterminate group most of the
# control arm.


def _missingness(df: pd.DataFrame, audit: Audit) -> None:
    """One structural entry per declared column, then one row per analysis column."""
    for column, reason in C.STRUCTURALLY_NON_APPLICABLE.items():
        n = int(df[column].isna().sum())
        audit.record("structural", column, n,
                     f"absent on {n} of {len(df)} records — structural, not data loss: {reason}")

    # The four centre columns come from CENTER_ORDER, never from the data and never from a literal.
    # From the data, a run in which one centre contributes no rows — the fixture, every future
    # subset, Stage 5's restricted cohort — silently renders a narrower table that still reconciles,
    # and the missing centre is the information. From a literal, CENTER_ORDER and the log drift
    # apart the first time a fifth centre joins.
    header = ("column", "kind", "n", "n_absent", "pct", *C.CENTER_ORDER)
    rows: list[tuple[str, ...]] = []
    labelled: list[str] = []
    for column in sorted(C.ANALYSIS_NAMES):
        absent = df[column].isna()
        n_absent = int(absent.sum())
        if column in C.STRUCTURALLY_NON_APPLICABLE:
            kind = "structural"
            labelled.append(f"{column} ({kind}): {C.STRUCTURALLY_NON_APPLICABLE[column]}")
        elif column in C.INFORMATIVE_ABSENCE:
            kind = "not recorded"
            labelled.append(f"{column} ({kind}): {C.INFORMATIVE_ABSENCE[column]}")
        else:
            kind = "missing" if n_absent else "complete"
        rows.append((
            column, kind, str(len(df) - n_absent), str(n_absent),
            _fmt(100 * n_absent / len(df)) if len(df) else _fmt(0),
            *(str(int((absent & (df["center"] == centre)).sum())) for centre in C.CENTER_ORDER)))

    detail = "\n".join([
        "one row per analysis column; kind is structural, not recorded, missing or complete. "
        "Nothing here is imputed.",
        *labelled])
    audit.record("missingness", "absence_by_column", len(df), detail,
                 table=(header, *rows))


# --- the pipeline -------------------------------------------------------------------------------

def _pipeline_after_read(df: pd.DataFrame, source: Source, audit: Audit) -> pd.DataFrame:
    """Everything after the read: normalise, correct, assert, report. The tests' seam.

    Named and exported to the tests because the shuffle test — the one that turns "never indexed by
    row" from a convention into a check — cannot go through load(), which takes a file. Driving
    _correct alone would leave _normalise and _missingness untested for row dependence, and both
    write into the log.
    """
    df = _normalise(df.copy(), audit)
    df = _correct(df, audit)
    _assert_schema(df, source)
    _missingness(df, audit)
    return df


def load(source: Source = WORKBOOK) -> tuple[pd.DataFrame, Audit]:
    """Read `source` into the analysis frame, with the audit log of everything that differed.

    Writes nothing. The caller decides whether the log is written, and where:

        df, audit = load()
        audit.write()          # → out/logs/stage2_audit_<label>.md, gitignored

    Raises DataVersionError if the file's bytes are not the pinned ones, and SchemaError if its
    columns, row count or values do not meet the contract.
    """
    audit = Audit(source)
    _verify_version(source)
    raw = _read(source, audit)
    df = _apply_contract(raw, source, audit)
    return _pipeline_after_read(df, source, audit), audit
