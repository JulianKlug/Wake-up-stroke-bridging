# Stage 2 spec — load and clean

Implements roadmap Stage 2 [§11]. Section references in brackets are to `statistical_analysis_plan.md`.
Numbers and decisions referenced as DECISION *n* are established in Stage 0 and recorded in
`../out/stage0_data_inventory.md`. Constants named in `SMALL_CAPS` are Stage 1's and live in
`config.py`; §1 of `stage1_config_and_data_contract.md` is their specification.

**Status.** Written 2026-08-08 against the workbook as verified in §18; revised 2026-08-10 after an
engineering review of this document. It is the **sole source for the Stage 2 implementation**:
everything the implementer needs is here, and anything not here is not to be invented.

Two prior statements are corrected, both marked **[REV]** with their evidence: the Stage 0 note's
claim that the workbook contains no `'N/A'` sentinels (§4.2), and the pilot's treatment of
`onset_to_groin_min = 999` as a placeholder (§7.3).

Eleven changes came out of the review, and eight further points that were settled in it are written
down here rather than left in a conversation nobody can read later. Each is marked **[ENG]** where it
appears, and §19 lists all nineteen with the defect or the omission each one closes — so a reader who
saw the 2026-08-08 draft can find what moved without re-reading the document.

**Nothing about Stage 2 is settled anywhere but here.** If a decision was made about this stage and
is not in this file, it is not a decision.

**Goal.** One reader that turns the private workbook into the analysis frame, and one audit log that
makes every difference between the two visible. Nothing is imputed, nothing is corrected by row
position, and no correction happens that the log does not name the patients for.

**Not in scope.** Deriving variables (Stage 3), classifying eligibility (Stage 4), restricting the
cohort (Stage 5). Stage 2 reads, renames, corrects, asserts and reports. It creates no new variable
except by recomputing one that [§6] already defines as a function of two others.

---

## 0. Where Stage 2 sits

```
  data/…v7…with_abs_contra_indication.xlsx        (gitignored, patient data)
                     │
                     ▼
  ┌────────────────────────────────────────────────────────────────────────┐
  │  STAGE 2 — data.py                                     reads config.py │
  │                                                                        │
  │   load(source)                                                         │
  │     │                                                                  │
  │     ├─ audit = Audit(source)     created first: every step below       │
  │     │                            records into this one object          │
  │     ├─ _verify_version(source)   sha256 vs DATA_SHA256                 │
  │     │            ▼                              → DataVersionError     │
  │     ├─ _read(source, audit)      READ_DTYPES, NA_VALUES, full parse    │
  │     │            ▼               records provenance/read               │
  │     ├─ _apply_contract(raw, source, audit)  assert_column_contract,    │
  │     │            ▼               rows, RENAME, drop → SchemaError      │
  │     │                            records contract/rename_and_drop      │
  │     └─ _pipeline_after_read(df, source, audit)  ← the tests' seam      │
  │            │                                                           │
  │            ├─ _normalise(df, audit)      CENTER_RECODE                 │
  │            │            ▼                                              │
  │            ├─ _correct(df, audit)        content-driven only,          │
  │            │            ▼                every case named              │
  │            ├─ _assert_schema(df, source) A1…A9 + A4b, every violation  │
  │            │            ▼                collected  → SchemaError      │
  │            └─ _missingness(df, audit)    structural / not recorded /   │
  │                         ▼                missing / complete            │
  └────────────────────────────────────────────────────────────────────────┘
       │
       ▼
  (df, Audit)          load() writes NOTHING. The caller decides:
       │                   audit.write()  →  out/logs/audit_<label>.md
       │                                     (gitignored: names patients)
       ├──────────────┬──────────────┬──────────────┐
       ▼              ▼              ▼              ▼
   Stage 3        Stage 4        Stage 5      Stage 14
   derive         eligibility    cohort       reporting (denominators)
```

**[ENG] `load()` has no filesystem side effect, and the log path carries the source's label.** Both
follow from one fact: acceptance §12 runs the whole pipeline against `FIXTURE` on every `pytest`
invocation, including on machines that have `data/`. A `load()` that wrote a fixed path would
overwrite the workbook's audit log with the fixture's two-row log on every test run, silently — the
file is gitignored, so no diff would ever show it. §9.5 gives the derived path; §17's definition of
done gives the two-line snippet that replaces the implicit write.

**The `Audit` is created first and threaded through, not built at the end. [ENG]** §9.6 puts two
entries — `provenance/read` and `contract/rename_and_drop` — inside `_read` and `_apply_contract`, so
those functions need the object before `_pipeline_after_read` exists. Threading it is also what lets
§9.5's header be rendered *from the entries* rather than from a second traversal of the source: the
row and column counts it prints are the ones `provenance/read` already carries, and the mapped and
dropped counts are the ones `contract/rename_and_drop` already carries. One number, one origin.

`load` is therefore, in full:

```python
def load(source: Source = WORKBOOK) -> tuple[pd.DataFrame, Audit]:
    audit = Audit(source)
    _verify_version(source)                                  # DataVersionError
    raw = _read(source, audit)                               # provenance/read
    df = _apply_contract(raw, source, audit)                 # contract/rename_and_drop
    df = _pipeline_after_read(df, source, audit)             # the remaining five entries
    return df, audit                                         # no file written
```

Every arrow out of Stage 2 carries analysis names only. `data.py` is **not** exempt from the Stage 1
§7 raw-name scan and must never become exempt; see §12.10.

### 0.1 Why three new types, and not fewer

Recorded because a reviewer will count them and ask. `data.py` introduces `Source`, `AuditEntry` and
`Audit`, and each earns its place by making a rule unrepresentable rather than merely documented:

- **`Source`** stops `load(check_hash=False)` from ever being writable at a call site (§3). A keyword
  that disables a guarantee is invisible from the module that declares the guarantee.
- **`AuditEntry`** carries the correction-names-its-cases rule in `__post_init__` (§9.2), so the rule
  holds for any code path that constructs an entry, not only for paths that go through `record`.
- **`Audit`** is the ordered collection and the renderer.

Collapsing `AuditEntry` into `Audit.record` was considered and rejected: it moves roadmap Stage 2's
acceptance criterion from a type into a method, and §12.7 would then test a call rather than a
constructor. Collapsing `Source` into module constants was considered and rejected: it reintroduces
the three parallel arguments — path, hash, row count — that §3 exists to bundle.

## 1. Deliverables

| Path | Contents |
|---|---|
| `extended_bridging/data.py` | the reader, the corrections, the assertions, the audit log |
| `extended_bridging/test_data.py` | the acceptance tests in §12 |
| `extended_bridging/config.py` | **amended** — one new constant, `INFORMATIVE_ABSENCE` (§10.2), and one corrected comment on `NA_VALUES` (§4.2) |
| `extended_bridging/test_config.py` | **amended [ENG]** — `INFORMATIVE_ABSENCE` joins 9.13's union, a disjointness test, one `'N/A'` cell in the fixture rows (§12.15) |
| `extended_bridging/tests/fixture_schema.xlsx` | **regenerated [ENG]** — from the amended `test_config.py`, via `--write-fixture` |
| `extended_bridging/specs/stage1_config_and_data_contract.md` | **amended [ENG]** — §5.1's "Stage 0 found no sentinels; this is a guard", corrected per §4.2 |
| `extended_bridging/pyproject.toml` | **amended** — `tabulate>=0.9` (§2.2), with `uv.lock` refreshed |

Six files, four of them one-line amendments. The three **[ENG]** rows are not scope growth: §4.2
proves a comment in `config.py` false and §12 requires a fixture behaviour the committed fixture does
not have, so the 2026-08-08 draft was under-counting its own deliverables rather than proposing less
work. §19 records why each was added.

`out/logs/audit_<label>.md` is an **output, not a deliverable**. It is written when a caller
asks for it, it is gitignored, and it must stay that way: it names patients by case identifier, which
is exactly what roadmap Stage 2's acceptance criterion demands of it. This is the same reason
`../out/stage0_data_inventory.md` is gitignored. Nothing under `specs/` may quote a case identifier.

## 2. Environment

### 2.1 Unchanged

`uv`, Python 3.12, dependencies as `stage1_config_and_data_contract.md` §2. Commands are run from
`extended_bridging/`; flat module layout, `import config as C`.

```bash
cd extended_bridging && uv sync && uv run pytest -v
```

### 2.2 `tabulate` is added, and Stage 2 does not use it  **[REV — new]**

`stage0_data_inventory.py` calls `DataFrame.to_markdown()` in five places. `to_markdown` requires the
optional dependency `tabulate`, which is in neither `pyproject.toml` nor `uv.lock`. **Stage 0
therefore cannot be regenerated in the committed environment** — it raises `ImportError` — while the
head of its own note says "Regenerate rather than edit" and roadmap Stage 0 claims re-running
reproduces the note byte-identically. Add `tabulate>=0.9` and re-lock, so that claim becomes true
again.

Stage 2's audit log nonetheless renders its own tables (§9.4). Byte-identical reproduction is an
acceptance criterion here, and a formatting library's column padding, alignment rules and float
repr are free to change between minor versions. A criterion that rests on a third party is a
criterion that fails on an unrelated `uv sync`.

### 2.3 What this costs to run  **[ENG — new]**

Recorded so that no later stage optimises something that was never slow, and so that §15's rejection
of a parsed-frame cache has a number behind it.

The workbook is 126 rows by 42 columns. `load()` reads the file twice — once as bytes for the sha256,
once through `read_excel` — and both are well under a second on the pinned openpyxl. Every correction
and every assertion is a single vectorised pass over at most 126 values. The audit log is a few
hundred formatted strings. §12.5's shuffle runs a six-record frame five times.

There is no N+1 pattern, no query, nothing to cache, and no memory concern: the whole cohort is
smaller than the module that describes it. If a future stage is slow it will be Stage 10's 2000
bootstrap refits, not this. **Do not add a cache, a lazy read, or a memoised `load()`** — §15 gives
the standing reason, which is that a cache is a second source of truth, not that it would be too
much work.

## 3. The source, and why it is a type

`load()` must run against two files: the private workbook, and `tests/fixture_schema.xlsx` — the
Stage 1 fixture, which carries the 42 verbatim headers and two invented rows. Without the second,
every Stage 2 test needs `data/`, and `data/` is gitignored.

The obvious shape, `load(path=None, check_hash=True, check_rows=True)`, is the wrong one: it lets a
caller switch off a guarantee with a keyword, and a keyword in a call site is invisible from the
module that declares the guarantee. Instead the file and its guarantees travel together.

```python
@dataclass(frozen=True)
class Source:
    path: Path
    sha256: str | None       # None: no version check. Only the fixture may be None.
    n_records: int | None    # None: no row-count check. Only the fixture may be None.
    label: str               # for the audit log header; never a path

WORKBOOK: Final[Source] = Source(
    C.DATA_XLSX, C.DATA_SHA256, C.N_RECORDS_EXPECTED,
    "v7_july26_with_abs_contra_indication")

FIXTURE: Final[Source] = Source(
    C.FIXTURE_XLSX, None, None, "fixture_schema")

SOURCES: Final[tuple[Source, ...]] = (WORKBOOK, FIXTURE)


def load(source: Source = WORKBOOK) -> tuple[pd.DataFrame, Audit]:
    ...
```

Two instances, and §12.1 asserts `SOURCES` against that literal — so a third source cannot appear
without a test failing, in the same way §9.5 of the Stage 1 spec pins `OUTCOME_MODEL_OVERRIDES`
against `{"tici_2b_3"}`. `WORKBOOK` carries both a hash and a row count; `FIXTURE` carries neither,
because it has 2 rows and no pinned bytes — the fixture generator in `test_config.py` may rewrite it.

Both sources use `C.SHEET` (`"Feuil1"`). The fixture uses that sheet name deliberately, so no `sheet`
field is needed and the read path is identical for both.

**[ENG] `SOURCES` pins what `data.py` declares, not what `Source` permits.** A test may construct its
own `Source` over a `tmp_path` file or over `FIXTURE`'s path with a deliberately wrong hash; §12.1's
literal comparison still fails if a third source is *declared at module level*, which is the drift it
exists to catch. Without this sentence an implementer reads "Only the fixture may be `None`" as a
prohibition on test-local instances and is left with no way to write §12.1's hash-mismatch case or
§12.5's shuffle case except against the private workbook — which would make both data-gated, and
Definition of done #2 requires the whole pipeline to run without `data/`. Three tests depend on this:

```python
# a wrong hash over a file that exists on every checkout — no data/ needed
WRONG = Source(C.FIXTURE_XLSX, "0" * 64, None, "fixture_wrong_hash")

# a hand-built frame's log needs a header; the label is what the header prints
HAND = Source(tmp_path / "hand.xlsx", None, None, "hand_built")
```

## 4. The read

### 4.1 The call

```python
def _read(source: Source, audit: Audit) -> pd.DataFrame:
    raw = pd.read_excel(
        source.path,
        sheet_name=C.SHEET,
        dtype=C.READ_DTYPES,        # from the contract; never a raw header written here
        na_values=list(C.NA_VALUES),
        keep_default_na=True,       # explicit, not defaulted. See §4.2.
    )
    audit.record("provenance", "read", n=len(raw), ...)     # §9.6
    return raw
```

**A full parse, never `nrows=0`.** The workbook's final column is headerless and pandas materialises
it only once rows are read: a header-only read returns 41 columns, a full read 42.
`assert_column_contract` has a dedicated branch that names this cause (Stage 1 §4), and it can only
help if Stage 2 never triggers it deliberately.

`READ_DTYPES` types 20 of the 42 columns — `string` where inference yields mixed types, `Int64`
wherever a column is a count, a score or a 0/1 flag, so that missingness survives as `<NA>` instead
of forcing the column to float. The volumes and times are left to inference. This reader is stricter
than `pilots/data.py`: a non-numeric value in a scored column raises at read time rather than being
absorbed as `NaN` (Stage 1 §4.4). That is intended. Do not widen `NA_VALUES` in response — add the
specific sentinel, or fix the workbook.

### 4.2 The workbook does contain `'N/A'`  **[REV]**

`../out/stage0_data_inventory.md` states, under Structural checks: "No `'N/A'` or empty-string
sentinels: missingness is already blank in every column." **That is not correct.** Read with
`dtype=object, keep_default_na=False`, the v7 workbook holds the literal three-character string
`N/A` in:

| Column | Literal `'N/A'` cells |
|---|---|
| `mRSscoreat90days` | 2 |
| `TICI_2b_3` | 3 |
| the six 24-hour NIHSS columns | 1 each |

They are invisible under a default read because `'N/A'` is already a member of pandas' default NA
list, so Stage 0 saw them as blanks and reported them as such. The counts line up exactly: the "2
missing" for `mrs_90d` and "3 missing" for `tici_2b_3` in the Stage 0 note *are* these strings.

Three consequences, all of which the implementation must respect:

- **`NA_VALUES` is load-bearing, not a guard.** Stage 1 §5.1 records it as "Stage 0 found no
  sentinels; this is a guard." It is doing real work on `mrs_90d`, the primary outcome, and on
  `tici_2b_3`.
- **`keep_default_na` must stay `True`, and is passed explicitly so that it is a decision rather than
  a default.** Setting it `False` — the plausible move for someone hardening the reader — would leave
  `'N/A'` as a string in an `Int64` column and crash the read. That failure is loud, so it is not
  dangerous; the danger is the repair, which would be to widen `NA_VALUES` blindly.
- **What `keep_default_na=True` costs**, stated so it is not discovered later: pandas also treats
  `NA`, `NULL`, `None`, `nan`, `-` and a dozen others as missing. Of the three `string` columns, none
  currently holds such a value — `contraindication_reason` has ten distinct texts, none NA-like, and
  `case_id` and `Center` are checked by §8. If a future workbook records a contraindication reason of
  literally `NA`, it would read as "no reason recorded" and move that patient from *eligible* to
  *indeterminate* [DECISION 1]. §8 A9 does not catch that case; nothing does. It is recorded in §13.

Also correct, while here: the Stage 0 note's header reads "126 rows, 43 columns". Its `load()` adds a
derived `center` column before the count is taken. The workbook has **42** columns, which is what
`N_COLUMNS_EXPECTED` says and what §12.14 asserts.

## 5. Applying the contract

`_apply_contract(raw, source, audit)`, in this order, and the order is the specification:

```python
def _apply_contract(raw: pd.DataFrame, source: Source, audit: Audit) -> pd.DataFrame:
    C.assert_column_contract(raw.columns)                  # SchemaError; Stage 1 §4
    if source.n_records is not None and len(raw) != source.n_records:
        raise SchemaError(...)                             # names both counts
    dropped = [key for key in C.COLUMN_CONTRACT if key in C.DROPPED]   # [ENG] declared order
    df = raw.rename(columns=C.RENAME).drop(columns=dropped)
    audit.record("contract", "rename_and_drop", n=len(dropped), ...)   # §9.6
    return df
```

The contract entry is recorded **after** both checks pass. A log that records what was dropped from a
frame that then failed `assert_column_contract` describes a frame that never existed.

Post-condition, asserted at §12.3: `set(df.columns) == C.ANALYSIS_NAMES`, 25 columns. Neither a raw
header nor a dropped column survives. That post-condition is what makes the Stage 1 §7 raw-name scan
worth anything downstream — the scan stops a module from *writing* `'Deathat90days'`, and this drop
stops it from *finding* it.

`RENAME` and `DROPPED` are derived views of `COLUMN_CONTRACT`; Stage 2 reads them and never rebuilds
either. The rename is applied before the drop so both operate on names the contract owns.

**[ENG] The drop list ranges over `COLUMN_CONTRACT`, never over `DROPPED` itself.** `DROPPED` is a
`frozenset` (Stage 1 §4.3), CPython randomises `str` hashing per process unless `PYTHONHASHSEED` is
set, and so `list(C.DROPPED)` is in a different order in every interpreter. For `.drop()` alone that
is harmless — pandas preserves the frame's column order regardless — but §9.6 renders the dropped
columns into the audit log's `contract` entry, and byte-identical reproduction (roadmap Stage 2's
acceptance criterion) is then lost *across processes* while still passing a same-process comparison.
§12.8 tests it across two interpreters with different seeds, because that is the only place it
breaks. This is the same rule as §9.3's "no set iteration", applied at the one site that predates it.

The version check precedes all of this (§7 of the Stage 1 spec's data-flow diagram): if the workbook
was replaced, "which file is this" is a more useful error than "column X is missing."

```
DataVersionError:
  expected 54934fbb…  (v7_july26_with_abs_contra_indication)
  actual   9a1c33be…
  The workbook changed. Re-run stage0_data_inventory.py, review
  ../out/stage0_data_inventory.md, and update DATA_SHA256 in the same commit.
```

## 6. Normalisation

Two operations, neither of which changes a value's meaning.

**`center`.** `df["center"] = df["center"].map(C.CENTER_RECODE)`, then assert totality (§8 A3). The
column is read as `string`, so HUG's integer code arrives as `"1"` and matches the contract's key.
The recode is logged as a `provenance` entry with the four counts, not as a correction: no datum is
altered, a code is given its name.

**`center` stays a `string`. It is not made a `Categorical` here.** Stage 5 drops the zero-bridging
centre, and a categorical that keeps a dead level makes every subsequent `groupby(observed=False)`
resurrect USZ as an all-missing row — in the balance table, the within-centre overlap table [§9] and
the cohort-flow table. Stage 1 §5.4 already assigns categorical construction to Stage 6, which builds
it from `FACTOR_LEVELS` at the point of use, where the level set is chosen deliberately.

**`case_id` is asserted, not stripped.** It arrives as `string` and is already clean (§18). Stripping
it would be a silent correction of the kind §7 forbids; asserting it (§8 A2) turns a future padded
identifier into an error naming the value.

**`contraindication_reason` is not normalised.** The pilot collapsed three spellings of "Clinician
decision" with a regex. Under DECISION 1 the free text is never read — only whether a reason is
present — so there is nothing to normalise, and normalising it would create the appearance that the
text feeds a classifier. §8 A9 protects the one bit that is read.

## 7. Corrections

### 7.1 The discipline

Roadmap Stage 2 and [§11] agree: every correction is **content-driven, matched on values, never
indexed by row**, and logged with the affected case identifiers.

Concretely, in `data.py`:

- No `.iloc[i]`, no `.loc[i]` with an integer, no `.index[…]`, no `.head`/`.tail` on the analysis
  frame. Every correction is a boolean mask over column values.
- The returned frame keeps the workbook's row order and a default `RangeIndex`. **That index is not
  an identifier and nothing may treat it as one.** `case_id` is the identifier.
- Every correction records an `AuditEntry` of kind `correction` carrying `df.loc[mask, "case_id"]`.
  An entry of kind `correction` with `n > 0` and no identifiers is rejected by the type (§9.2).

§12.5 turns this from a convention into a test: shuffling the input rows must produce an identical
frame after sorting by `case_id`, **and a byte-identical audit log**. Any use of row position breaks
one or the other. It is the Stage 2 counterpart of Stage 1's AST scan for raw names.

### 7.2 `penumbra_ml` is recomputed from its definition

[§6] excludes `penumbra_ml` from every model because it is "a deterministic function of core and
Tmax>6 s volumes; including all three makes the design exactly singular." That sentence makes the
definition authoritative and the stored column merely a copy of it, so Stage 2 recomputes it:

```python
computed = df["tmax6_ml"] - df["core_ml"]
disagree = (df["penumbra_ml"] - computed).abs() > _TOL      # _TOL = 1e-6
df["penumbra_ml"] = computed
```

Applied to **every** record, not only the disagreeing ones. A patch of the disagreeing rows leaves
the stored column authoritative wherever it happens to agree, and a second copy-paste that is
self-consistent by accident would survive it. Recomputing everywhere is content-driven by
construction: there is no mask to get wrong and no row to name.

Missingness is preserved without a `.mask()` because subtraction propagates it: `core_ml` and
`tmax6_ml` are missing on the same single record, which is also the record where `penumbra_ml` is
missing. §12.6 asserts the preservation rather than relying on the arithmetic.

**A record with a missing volume is never counted as disagreeing. [ENG]** `(NA - x).abs() > _TOL` is
false, not missing, under both numpy and pandas nullable semantics, so `disagree` excludes those
records without a guard — which is the behaviour we want and is why the guard is absent. It is stated
because it determines the correction entry's `n`, and `n` is what §9.6 renders and §12.6 asserts: a
reader who assumed missing rows counted as disagreements would read `n = 1` in v7 as evidence of a
bug rather than as the answer. `_TOL = 1e-6` is an absolute tolerance, which is right here because
the volumes are millilitres in the tens and hundreds, recorded to at most one decimal.

**What it corrects in v7.** Exactly one record — HUG, treated — stores `core_ml = 13`,
`tmax6_ml = 109`, `penumbra_ml = 13`. The penumbra cell holds the core value; the discrepancy is
83 mL. The pilot found the same defect and its audit entry read "one row had the core volume
copy-pasted into the penumbra cell."

**The consequence, stated because it is not cosmetic.** `penumbra_ml` feeds the balance table as a
`BALANCE_ONLY` variable [§6, §9], and an 83 mL error on one arm of a 126-record cohort is visible
there. A reader comparing the arms on penumbra volume would otherwise be comparing one number that
came from the definition against 125 that came from a copy-paste.

**[AMENDED 2026-08-10]** This section used to state a second consequence: that the correction moved
this patient across the `target_mismatch` subgroup's volume threshold. That subgroup is withdrawn by
the [§13] amendment, so the consequence no longer exists — and the two constants the sentence named
are deleted from `config.py`. The recomputation itself is untouched and its justification is
unchanged: it rests on [§6] defining penumbra as `tmax6_ml - core_ml`, which was never part of any
subgroup. See `stage3_derived_variables.md` §6.1.

The audit entry names the case and states the stored and recomputed values.

### 7.3 `onset_to_groin_min = 999` is flagged, not corrected  **[REV]**

`pilots/data.py:103` sets `onset_to_groin_min == 999` to missing, on the reasoning that 999 is a
placeholder. **Stage 2 does not carry that correction.** The value stands; the record is logged as an
`observation` naming the case, and it goes to the data owner as a standing query.

The evidence, all of it from v7:

| | |
|---|---|
| records with `onset_to_groin_min == 999` | 1 (HUG, treated) |
| that record's `onset_to_ivt_min` | 955 — so the IVT-to-groin interval is **44 minutes** |
| the same interval across all 39 treated patients | min 29, median 70, max 310 |
| observed `onset_to_groin_min` range | 253 – 1219 |
| other observed values in 900–1100 | 905, 905, 905, 950, **999**, 1035, 1045 |

A 44-minute interval between thrombolysis and groin puncture is unremarkable — it sits between the
minimum and the median. The value is inside the observed range, and neighbouring values are equally
"round". The entire case for calling it a placeholder is the digit pattern, and setting a plausible
value to missing on a digit pattern deletes real data.

This is the same posture DECISION 2 took toward the two records whose shipped death flag contradicts
their mRS: a derivation rule where one is available, a logged standing query where one is not, and
never a silent resolution. Here no derivation is available at all — nothing else in the workbook
determines a groin-puncture time — so the query is the whole of the response.

The stakes are low and should be stated as such: `onset_to_groin_min` is post-time-zero, is in
`POST_TIME_ZERO`, and [§12] confines it to a by-arm descriptive table. It enters no model under any
specification. If the data owner confirms it is a placeholder, the fix is one correction entry in
this section, and it changes one cell of one descriptive table.

### 7.4 Pilot corrections deliberately not carried

| `pilots/data.py` | Why not |
|---|---|
| `nihss_change_24h` recomputed from `nihss_24h − nihss_baseline` | all six 24-hour NIHSS columns are dropped by the contract; they are not in the [§5] registry, and restoring them requires amending [§5] first |
| `contraindication_ivt` regex normalisation | DECISION 1 — the free text is never read (§6) |
| `_rebuild` of the four dichotomies | Stage 3 owns them; Stage 1 §3.2 already gives the loop, driven by `OUTCOMES` and `OPS` |
| `_eligibility` | Stage 4, and it is built on `IVT_INELIGIBLE_REASONS` / `IVT_ELIGIBLE_REASONS`, which DECISION 1 retires |
| `center_hug`, `unknown_onset`, `mrs_shift_90d`, `onset_type` | Stage 3 — except `center_hug` and `mrs_shift_90d`, which are not in the [§5] registry, and `target_mismatch`, whose subgroup the [§13] amendment of 2026-08-10 withdrew |
| `hir` missingness entry | `HIR` has no analysis name in the contract at all, precisely so that no covariate list can name it [§6] |
| the `assert` statements in `_check` | the checks are lifted; the statement form is not (§8.3) |

## 8. Schema assertions

Run **after** the corrections, once, over the analysis frame. After, not before: `penumbra_ml` is
recomputed in full, so asserting the stored column first would be asserting something Stage 2 is
about to discard.

### 8.1 The checks

Each reads a Stage 1 constant where one exists, and each is justified by the specific silent failure
it prevents. Each names the column and the offending values.

**[ENG] Every check runs, and one `SchemaError` reports all the violations.** Not fail-fast. This
module's job is to meet a workbook it has never seen — the data owner will send a v8 — and a
first-failure exception turns "three columns drifted" into three read-hash-parse-fail cycles, one per
round trip, with the operator learning one fact per run. Each check appends zero or more lines to a
list; `_assert_schema` raises once at the end if the list is non-empty.

```
SchemaError:
  A4  ivt: 2 record(s) missing the treatment
  A6  mrs_90d: value(s) outside (0, 6): 7
  A6  core_ml: value(s) outside (0, None): -1
  3 assertion(s) failed against 126 records. Every one is a contract
  violation; none of them is a reason to relax the contract.
```

**Each check guards its own preconditions**, so no check depends on another having passed and the
table's order is presentation, not semantics. A7 in particular skips records where `ivt` is `<NA>`
and lets A4 report those — the alternative, letting A7 compare against a missing exposure, produces a
second confusing message about the same two records.

| | Check | Reads | The silent failure it prevents |
|---|---|---|---|
| A1 | row count equals `source.n_records`, when the source declares one | `N_RECORDS_EXPECTED` | a filtered sheet or a partial read producing a smaller cohort that every downstream denominator then reports as fact |
| A2 | `case_id` non-missing, unique, and equal to its own stripped form | — | a duplicated identifier double-weighting one patient through the propensity fit and the bootstrap; a padded identifier that fails to join to the audit log |
| A3 | every `center` value is in `CENTER_RECODE.values()` | `CENTER_RECODE` | a new or re-spelt centre code mapping to `<NA>` and quietly dropping that centre from every by-centre table and from the [§3] restriction |
| A4 | `TREATMENT` is never missing | `TREATMENT` | a missing exposure becoming a third arm, or being counted as a control |
| A4b **[ENG]** | `ivt_contraindicated` is never missing | — | a `<NA>` classifier input reaching Stage 4, which [DECISION 1] builds the whole *eligible* / *indeterminate* split on. §14 promises Stage 4 this guarantee; A5 does not deliver it, because A5 permits missing everywhere |
| A5 | every column in `BINARY_COLUMNS` takes values in `{0, 1}` or missing | `BINARY_COLUMNS` | `ivt = 2` flowing into every weight, estimate and arm label without an error; `ivt_contraindicated = 2` defeating the [§4] classifier [DECISION 1] |
| A6 | every column in `PLAUSIBLE_RANGES` lies inside its `(low, high)`, inclusive; `None` means unbounded above | `PLAUSIBLE_RANGES` | an mRS of 7 creating a seventh outcome category in the proportional-odds fit; a negative volume; a pre-stroke mRS of 6, which is death |
| A7 | `onset_to_ivt_min` is non-missing for exactly the patients with `ivt == 1` | `STRUCTURALLY_NON_APPLICABLE` | a treated patient with no IVT time being read as ordinary missingness [§11]; a control with an IVT time, which would mean the arm label is wrong |
| A8 | `onset_to_ivt_min <= onset_to_groin_min` wherever both are present | — | a time ordering under which thrombolysis follows puncture, i.e. the record is not bridging at all |
| A9 | no non-missing `contraindication_reason` is empty after stripping | — | a whitespace-only cell reading as "a reason was recorded" and moving a patient from *indeterminate* to *eligible* [DECISION 1, §3] |

A6 covers `penumbra_ml`, so after §7.2 it also asserts `tmax6_ml >= core_ml` — a negative recomputed
penumbra would mean the two volume columns are swapped.

**A1 is a post-condition and duplicates §5 deliberately. [ENG]** Nothing between §5 and here adds or
removes a row: the rename and drop are column-wise, the recode and the recomputation are value-wise.
So A1 cannot fire today, and it is kept anyway — it is the check that catches a future Stage 2 edit
which filters rows, at the cost of one comparison on 126 records. It is recorded as unreachable so
that nobody deletes it as dead code and nobody expects §12.9 to test it through `load()`.

All ten hold in v7 (§18). They are written for the workbook that has not arrived yet.

### 8.2 What Stage 2 does not assert

`wake_up` and `unwitnessed` are never both positive. That is roadmap Stage 3's acceptance criterion,
asserted where `onset_type` is built, and duplicating it here would give one invariant two owners.

**No treated patient carries `ivt_contraindicated == 1`. [ENG]** `config.py`'s contract entry for
that column assigns the assertion to Stage 4 — "Stage 4 asserts it is 0/1, never missing, and never 1
for a treated patient" — and Stage 4 is where the classifier lives, so that is where a contradiction
between the flag and the exposure has to be resolved. Stage 2 supplies the two preconditions Stage 4
needs (A5 for the domain, A4b for non-missingness) and nothing more. §12.14 records the v7 count as a
data-gated fact, not as a runtime check, for the reason §8.2 gives below about Stage 0's counts.

The arm sizes (39 treated, 87 control) and the other Stage 0 counts are not runtime assertions.
`DATA_SHA256` already pins the file, so a runtime check adds nothing; they appear instead as
data-gated tests (§12.14), following the precedent of Stage 1's acceptance 9.14 — the test documents
the dependency, so updating the hash also re-checks the premises.

### 8.3 No bare `assert`  **[REV]**

`pilots/data.py:177-193` states every one of its checks as `assert`. Python strips `assert` under
`-O`, so under a flag nobody remembers setting, a module whose entire purpose is to fail loudly
succeeds silently. Stage 2 raises `SchemaError` explicitly, everywhere, and §12.10 scans `data.py`
for `ast.Assert` and fails if it finds one.

`SchemaError` is reused rather than a new exception being minted. The contract is columns *and* the
domains those columns are declared to take — `PLAUSIBLE_RANGES` and `BINARY_COLUMNS` are as much a
part of it as `COLUMN_CONTRACT` — so a range violation is a contract violation. `DataVersionError`
covers the file's identity. Stage 1 declared both and nothing else is needed.

## 9. The audit log

### 9.1 What it has to be

Roadmap Stage 2 accepts when "the audit log reproduces byte-identically on a re-run, and every
correction in it names the cases it touched." Both clauses are properties of the *rendering*, so both
are specified here rather than left to the implementation.

### 9.2 The types

```python
KINDS: Final[tuple[str, ...]] = (
    "provenance", "contract", "correction", "observation",
    "derivation",                   # [AMENDED 2026-08-10] Stage 3's, inserted rather than appended
    "structural", "missingness")

@dataclass(frozen=True)
class AuditEntry:
    kind: str
    step: str
    n: int
    detail: str
    case_ids: tuple[str, ...] = ()
    table: tuple[tuple[str, ...], ...] | None = None   # header row first

    def __post_init__(self) -> None:
        if self.kind not in KINDS:
            raise ValueError(...)
        if self.kind in _MUST_NAME_CASES and self.n > 0 and not self.case_ids:
            raise ValueError(
                f"{self.kind} {self.step!r} touched {self.n} record(s) but names no case. "
                "Roadmap Stage 2 requires every correction to name the cases it touched, "
                "and an observation nobody can look up is not a query.")

_MUST_NAME_CASES: Final[frozenset[str]] = frozenset({"correction", "observation"})   # [ENG]
```

The rule lives in the type, not in a reviewer's attention. A correction that changed nothing carries
`n = 0` and no identifiers, which is legitimate and is what §12.6 checks on a clean frame.

**`observation` is held to the same rule as `correction`. [ENG]** §9.6 gives `onset_to_groin_999`
`case_ids`, and §7.3 makes it a standing query with the data owner — a query that names no patient
cannot be answered, which makes it exactly as useless as an unattributed correction. Roadmap Stage 2
demands the rule only of corrections, so this is a strengthening, and it is free: the one observation
this spec defines already carries its identifiers. The other four kinds legitimately carry none —
`provenance`, `contract` and `missingness` describe columns rather than patients, and `structural`
describes a whole arm (§9.6).

```python
class Audit:
    def __init__(self, source: Source) -> None
    def record(self, kind, step, n, detail, case_ids=(), table=None) -> None
    def to_markdown(self) -> str
    def write(self, path: Path | None = None) -> Path      # [ENG] default derived, see §9.5
```

`record` normalises identifiers on the way in: `tuple(sorted(str(c) for c in case_ids))`.

**Sorted, not in frame order.** Sorting is what makes the log invariant to input row order, which is
what lets §12.5 test "never indexed by row" by shuffling the input and comparing the rendered log
byte for byte. Frame order would be deterministic too, and would fail that test for a correct
implementation.

### 9.3 The rendering rules

Byte identity is not automatic. Every one of these is a way it has been lost before:

- **No timestamp, no clock time, no run identifier.** The log is a function of the data, not of when
  it was produced. §12.8 greps the rendered text for an ISO-8601 date and for `HH:MM:SS`, so the
  property cannot pass merely because two runs landed in the same second.
- **No absolute path.** The header prints `source.path.name` and `source.label`, never a path that
  varies by checkout.
- **No pandas or numpy repr.** Every string in the log is built by `data.py`. A pandas minor release
  is free to change column padding, float repr or `<NA>` rendering.
- **No set iteration and no dict insertion assumptions beyond config's own.** Where a loop ranges
  over columns it ranges over a declared tuple — `BINARY_COLUMNS`, `PLAUSIBLE_RANGES`,
  `ANALYSIS_NAMES` sorted — never over a `set`.
- **One float formatter, declared once:** `f"{float(x):.6g}"`, with missing rendered as the literal
  `missing`. Six significant figures is exact for every volume, time and score in this workbook and
  is stable across platforms.
- **Sorted identifiers** (§9.2), rendered as an inline code span, comma-separated.

### 9.4 The table helper

```python
def _md_table(rows: Sequence[Sequence[str]]) -> str
```

First row is the header; every cell is already a string; each column is padded to its widest cell;
the separator row is `---` per column. Roughly fifteen lines, no dependency, no version to drift.
§2.2 says why this is not `tabulate`.

### 9.5 The document

Sections in `KINDS` order, entries within a section in the order recorded:

**[AMENDED 2026-08-10]** The title and footer below name the pipeline rather than this stage, for the
same reason the path did: Stage 3 appends its `Derivations` section to this document. `KINDS` gains
`derivation` between `observation` and `structural`, so the section list below is seven rather than
six. Everything else about the header is unchanged — in particular its numbers are still read back
from the entries, never recomputed. See `stage3_derived_variables.md` §9.1 and §9.2.

```
# Analysis audit

Source: `Excel_bridging_…_v7_july26_with_abs_contra_indication.xlsx`, sheet `Feuil1`
Version: `54934fbb…f657371`            (or `not pinned (fixture_schema)`)
Read: 126 rows x 42 columns; 25 mapped, 17 dropped
Generated by the `extended_bridging` pipeline. Regenerate rather than edit.

## Provenance
## Contract
## Corrections
## Observations
## Derivations                      ← Stage 3's, inserted here rather than appended
## Structural non-applicability
## Missingness and denominators
```

Each entry renders as the pilot's shape, which reads well and is worth keeping:

```
- **<step>** (n=<n>): <detail>
  Affected: `<id>, <id>, …`
```

with the `Affected:` line omitted when there are no identifiers, and any table indented beneath. A
section whose entries are all absent still prints its heading, followed by the single line `_none_` —
so `Observations` on a clean workbook says "nothing was observed" rather than looking truncated, and
the headings are a fixed skeleton the eye can scan for. **[AMENDED 2026-08-10]** There are seven of
them since `KINDS` gained `derivation`; a test that treats two headings as adjacent must take the
neighbour from `KINDS` rather than naming it, or the next inserted kind breaks it.

**The header's numbers come from the entries, not from a second traversal. [ENG]** `Source` and
`Version` come from `self.source`; `Read: {n_raw} rows x {n_cols} columns` comes from the
`provenance/read` entry; `{n_mapped} mapped, {n_dropped} dropped` comes from the
`contract/rename_and_drop` entry. Recomputing them in the renderer would give the log two sources for
one number, and the one that drifts is always the one nobody tests. This is why the `Audit` is
created at the top of `load()` and threaded through the read (§0).

**The path is derived from the source's label, and `load()` never writes it. [ENG]**

```python
def _audit_path(source: Source) -> Path:
    # [AMENDED 2026-08-10] `audit_`, not `stage2_audit_`. The Audit this module creates is threaded
    # through the whole pipeline and Stage 3 appends to it, so the artefact is the analysis's rather
    # than this stage's [stage3_derived_variables.md §9.1]. The shape and the label are unchanged,
    # so the guarantee below holds exactly as before.
    return C.LOGS / f"audit_{source.label}.md"

# Audit.write(path=None) → path or _audit_path(self.source)
```

`write` creates `C.LOGS` if absent and returns the path. `out/` is gitignored; §17's definition of
done checks that with `git check-ignore`, because this file names patients.

A single fixed `AUDIT_MD` was the 2026-08-08 draft's design and it does not survive contact with §12:
acceptance 12.2, 12.4, 12.5 and 12.8 all run the pipeline against `FIXTURE`, so on any machine that
has `data/`, one `pytest` run would replace the workbook's audit log with the fixture's two-row log.
Nothing would report it — the file is gitignored, so it appears in no diff, and the next reader would
review an audit of two invented patients believing it described 126 real ones. Deriving the path from
`source.label` makes the collision unrepresentable rather than merely unlikely, and it costs one
f-string. §12.8 asserts the two paths differ; §12.8 also asserts `load()` creates no file at all.

### 9.6 The entry inventory  **[ENG — new]**

Byte-identical reproduction is roadmap Stage 2's acceptance criterion, and §12.8 tests it by
comparing a run against itself. That is necessary and not sufficient: it holds for *any* log,
including one whose contents an implementer chose freely. The 2026-08-08 draft specified four entries
across the declared `KINDS`, leaving the `Contract` and `Structural` headings of §9.5 with nothing to
put under them, and pinned no `step` string, no `detail` wording and no meaning for `n`. Two
implementers would have produced two different logs and both would have passed §12.

This is the log. Seven entries, in this order. `step` strings are exact.

| kind | step | recorded by | `n` counts | table | `case_ids` |
|---|---|---|---|---|---|
| `provenance` | `read` | `_read` | rows read | — | — |
| `provenance` | `centre_recode` | `_normalise` | rows recoded | 4 rows + header | — |
| `contract` | `rename_and_drop` | `_apply_contract` | columns dropped | 42 rows + header | — |
| `correction` | `penumbra_recomputed` | `_correct` | records whose stored value disagreed | `n` rows + header, only when `n > 0` | those records |
| `observation` | `onset_to_groin_999` | `_correct` | records matching | — | those records |
| `structural` | one per `STRUCTURALLY_NON_APPLICABLE` key, `step` = the column name | `_missingness` | records where the column is absent | — | — |
| `missingness` | `absence_by_column` | `_missingness`, through `absence_by_column` | rows in the frame | 25 rows + header | — |

`detail` templates, with `{}` filled from the data and every number through §9.3's one formatter:

- `read` — `sheet {sheet!r}; {n_raw} rows x {n_cols} columns; version {state}`, where `state` is
  `pinned and verified ({sha256:.8})` or `not pinned ({label})`.
- `centre_recode` — `CENTER_RECODE applied; {k} codes mapped to {k} labels; no value altered, a code
  was given its name`. The table is `code | label | n`, rows in `CENTER_ORDER`, one row per label.
- `rename_and_drop` — `{n_mapped} columns mapped to analysis names, {n_dropped} dropped; every
  dropped column's reason is in config.COLUMN_CONTRACT`. The table is `raw column | fate`, one row
  per key **in `COLUMN_CONTRACT` order** (§5's `[ENG]` note), `fate` being the analysis name or the
  literal `dropped`. The whole difference between workbook and frame, on one table, which is what §0
  means by "makes every difference between the two visible".
- `penumbra_recomputed` — `[§6] defines penumbra as tmax6_ml - core_ml, so the definition is
  authoritative and the stored column is a copy of it. Recomputed on all {n_rows} records; {n}
  disagreed beyond {tol}.` The per-record values go in the table, not in the sentence:
  `case_id | stored | recomputed | difference`, rows sorted by `case_id`, present only when `n > 0`.
  A sentence naming "stored X vs recomputed Y" reads correctly for the one record v7 has and becomes
  ambiguous for the second one a corrected workbook brings. With `n = 0` the table is omitted and
  `case_ids` is empty. **[AMENDED 2026-08-10]** A second sentence used to follow, counting the
  records whose recomputed penumbra crossed the withdrawn [§13] subgroup's volume threshold; it is
  deleted with the subgroup (§7.2). The table already carried stored, recomputed and difference per
  case, so the entry loses no information and no test asserted the sentence.
- `onset_to_groin_999` — `onset_to_groin_min == 999 on {n} record(s). Not corrected (§7.3): the value
  is inside the observed range {lo}-{hi}, the implied IVT-to-groin interval is {gap} min against a
  treated-arm median of {median}, and neighbouring observed values are equally round. Standing query
  with the data owner.`
- structural, per column — the reason string from `STRUCTURALLY_NON_APPLICABLE`, verbatim, prefixed
  `absent on {n} of {n_rows} records — structural, not data loss: `.
- `absence_by_column` — `one row per analysis column; kind is structural, not recorded, missing or
  complete. Nothing here is imputed.` The table is §10.3's.

Five entries are unconditional — `read`, `centre_recode`, `rename_and_drop`, `penumbra_recomputed`
and `absence_by_column` — with `penumbra_recomputed` carrying `n = 0` and no identifiers on a frame
where every stored value already agrees. `onset_to_groin_999` appears only when `n > 0`. The
structural entries are one per key of `STRUCTURALLY_NON_APPLICABLE`, so today exactly one. §12.16
asserts the inventory against this table, so an implementation cannot quietly drop an entry, rename a
step or reorder a section.

## 10. Missingness and denominators [§11]

### 10.1 The requirement

[§11] is five sentences and Stage 2 answers four of them: complete-case with the denominator reported
for every estimate; no imputation; structural non-applicability distinguished from missingness and
never imputed; every correction content-driven and logged. (The fifth, dichotomies rebuilt with
missingness reimposed, is Stage 3.)

**Stage 2 imputes nothing.** §12.12 asserts this mechanically rather than by inspection: for every
column, the cleaned frame's missingness pattern equals the raw frame's, with `penumbra_ml` the single
declared exception, where the pattern is inherited from `core_ml` and `tmax6_ml`.

### 10.2 Absence has three kinds, and one of them is new  **[REV to Stage 1]**

A single "n missing" column misrepresents this cohort. `onset_to_ivt_min` is absent for all 87
controls because they were never given IVT — that is structure, not data loss, and `config.py`
already says so in `STRUCTURALLY_NON_APPLICABLE`. `contraindication_reason` is absent for 82 of 126
records, and **that absence is itself the signal**: under DECISION 1 it is the one bit the column
contributes, the thing that separates *eligible* from *indeterminate* [§3]. Reporting 82/126 beside
genuine data loss invites a reader to see the field as two-thirds unusable, when Stage 4 uses every
one of those 82.

There is no constant for that third kind, so Stage 2 adds one to `config.py`:

```python
# Absence that carries information rather than data loss. `contraindication_reason` is missing for 82
# of 126 records, and that absence is exactly the bit [DECISION 1] reads: no reason was recorded,
# which [§3] forbids reading as "no contraindication" and Stage 4 maps to `indeterminate`. Listing
# those 82 beside genuine data loss misrepresents both. Stage 2 labels them; nothing imputes them.
INFORMATIVE_ABSENCE: Final[dict[str, str]] = {
    "contraindication_reason":
        "no reason was recorded — the one bit [DECISION 1] reads [§3]",
}
```

Declared in `config.py`, not in `data.py`, because that is where this project declares facts about
columns, and because Stage 14 will need the same distinction when it prints denominators. It is
recorded here as an amendment rather than added quietly; §12.11 covers it.

**[ENG] Two Stage 1 tests must move with it, and the 2026-08-08 draft said one already covered it.**
That draft's task T2 verified itself with "`pytest test_config.py` still green (9.13 covers the new
names)". It does not: 9.13 parametrises over `PS_COVARIATES_FULL | BALANCE_ONLY | BINARY_COLUMNS |
PLAUSIBLE_RANGES | STRUCTURALLY_NON_APPLICABLE | POST_TIME_ZERO`, and `INFORMATIVE_ABSENCE` is in
none of them, so the suite passes because nothing reads the new constant. A typo in the key — say
`contraindication_reasons` — would surface only in §10.3, as a column silently classified `missing`
instead of `not recorded`, which is the exact misreading §10.2 exists to prevent. So:

- `INFORMATIVE_ABSENCE` joins 9.13's union, and every key must resolve to an analysis name.
- A new Stage 1 test asserts `STRUCTURALLY_NON_APPLICABLE` and `INFORMATIVE_ABSENCE` are **disjoint**.
  §10.3's `kind` column has one branch per constant; an overlap would resolve by branch order rather
  than by declaration, and the column that carries the [DECISION 1] bit is the one at stake.

Both are §17's task T2 and are listed in §1's deliverables.

### 10.3 The table

One row per analysis column, in `sorted(ANALYSIS_NAMES)` order. The header is built as

```python
header = ["column", "kind", "n", "n_absent", "pct", *C.CENTER_ORDER]     # [ENG]
```

which renders, today, as

| column | kind | n | n_absent | pct | HUG | CHUV | Lugano | USZ |
|---|---|---|---|---|---|---|---|---|

**The four centre columns come from `CENTER_ORDER`, never from the data and never from a literal.
[ENG]** From the data, a run in which one centre contributes no rows — the fixture, every future
subset, and Stage 5's restricted cohort — silently renders a narrower table that still reconciles,
and the missing centre is the information. From a literal, `CENTER_ORDER` and the log drift apart the
first time a fifth centre joins. `CENTER_ORDER` is already the declared display order (Stage 1 §5.4)
and `test_config.py` already pins it equal to `CENTER_RECODE`'s values, so one expression inherits
both guarantees. A centre with no rows renders a column of zeros, which is the honest rendering.

`kind` is one of:

- `structural` — the column is a key of `STRUCTURALLY_NON_APPLICABLE`; the reason string is printed
  beneath the table
- `not recorded` — a key of `INFORMATIVE_ABSENCE`; likewise
- `missing` — anything else with a non-zero absence count
- `complete` — no absences

The per-centre columns are there because Stage 0 established that the pattern is centre-driven, not
random: the contraindication reason was recorded for controls at HUG and USZ and nowhere else, and
`HIR` — dropped by the contract — was missing for all of CHUV and Lugano. A pooled count hides that
shape entirely, and it is the shape that makes the *indeterminate* group 80% of the control arm.

Recorded as a `missingness` entry carrying the table. Its row count reconciles to the frame's row
count, which §12.13 checks.

## 11. Data flow into Stage 3

Stage 3 receives `(df, audit)` and may rely on exactly this:

```
  df
   ├─ 126 rows, workbook order, default RangeIndex — NOT an identifier
   ├─ exactly the 25 columns of ANALYSIS_NAMES, no raw header, no dropped column
   ├─ case_id      string, unique, non-missing, unpadded
   ├─ center       string, one of the four CENTER_RECODE values — not a Categorical
   ├─ typed columns per READ_DTYPES; Int64 where a score or flag, so <NA> survives
   ├─ penumbra_ml == tmax6_ml - core_ml exactly, everywhere
   ├─ ivt and ivt_contraindicated are 0/1 and never missing  (A4, A4b)
   ├─ missingness identical to the workbook's, penumbra_ml excepted (§10.1)
   └─ assertions A1-A9 (including A4b) have passed, all of them, in one pass

  audit
   └─ Audit, already carrying the seven §9.6 entries. Stage 3 appends; it does
      not re-open. No file has been written unless a caller asked for one.
```

Stage 3 then builds `onset_type`, the four dichotomies through `OUTCOMES` and `OPS`, and the [§13]
subgroups. It appends its own entries to the same `Audit`, so one log covers the pipeline.

## 12. Acceptance criteria

All tests live in `test_data.py`, with section banners matching these numbers, as `test_config.py`
does for Stage 1. Tests needing the private workbook reuse the `DATA_GATED` idiom
(`pytest.mark.skipif(not config.DATA_XLSX.exists())`); everything else runs against `FIXTURE` or a
hand-built frame.

> **`test_data.py` is not exempt from the Stage 1 §7 raw-name scan, and must not become exempt.**
> That scan walks `extended_bridging/**/*.py` and fails on any string literal equal to a contract
> key. So the tests build their frames with **analysis names**, use `FIXTURE` for the read path, and
> derive any raw header they need from `COLUMN_CONTRACT` keys rather than writing one. Adding
> `data.py` or `test_data.py` to `EXEMPT_FROM_RAW_NAME_SCAN` would disable the rule for the two
> modules that handle the workbook — which is the whole of what it guards.

1. **Version and provenance.** A file whose hash differs raises `DataVersionError` naming both hashes
   and the remedy — tested **[ENG]** with a test-local `Source` over `FIXTURE`'s path carrying a
   deliberately wrong hash, so the check runs on a checkout with no `data/` rather than only under
   `DATA_GATED`. `SOURCES == (WORKBOOK, FIXTURE)`, asserted against that literal so a third
   module-level source fails the test rather than appearing quietly. `WORKBOOK.sha256` and
   `WORKBOOK.n_records` are both set; `FIXTURE`'s are both `None`. `Source` rejects attribute
   assignment.
2. **Read contract.** Reading `FIXTURE` with `READ_DTYPES` and `NA_VALUES` succeeds, yields 42
   columns and a `string` dtype for `case_id`. A header-only read raises `SchemaError` through
   `assert_column_contract` with the header-only cause named. A frame whose row count differs from
   `source.n_records` raises `SchemaError` naming both counts; a source with `n_records = None` does
   not raise. The fixture's literal `'N/A'` cell in `tici_2b_3` reads as `<NA>` and the column keeps
   its `Int64` dtype **[ENG]** — the 2026-08-08 draft asserted this without the fixture having such a
   cell, so the criterion could only ever have run data-gated; §12.15 adds the cell.
3. **Contract application.** After rename and drop, `set(df.columns) == ANALYSIS_NAMES` and
   `len(df.columns) == 25`. No key of `COLUMN_CONTRACT` and no member of `DROPPED` survives.
4. **Normalisation.** Every `center` value is in `CENTER_RECODE.values()`. An unmapped code raises
   `SchemaError` naming the offending value. `center` is not a `CategoricalDtype`. `case_id` is not
   stripped — a padded identifier raises rather than being silently repaired.
5. **Corrections are content-driven.** Shuffling the input rows yields a frame identical after
   `sort_values("case_id").reset_index(drop=True)`, **and a byte-identical audit log**. Run on a
   hand-built frame that exercises both §7.2 and §7.3, so the property is tested where it could
   actually break.

   **The harness, stated because the 2026-08-08 draft left it open. [ENG]** `load()` takes a file, so
   this test cannot go through it; it drives `_pipeline_after_read(df, source, audit)` — the named
   seam in §0 covering normalise, correct, assert and missingness — with a test-local `Source` (§3),
   a fresh `Audit` per shuffle, and a contract-valid hand-built frame of six records: one carrying
   the copy-paste penumbra of §7.2, one carrying `onset_to_groin_min = 999`, one missing `core_ml`
   and `tmax6_ml` together, one control and one treated so A7 has both sides, and one ordinary
   record. Shuffled under five fixed seeds; frames compared after sorting, rendered logs compared
   byte for byte.

   ```python
   HAND = Source(tmp_path / "hand.xlsx", None, None, "hand_built")
   base, base_md = None, None
   for seed in (0, 1, 2, 3, 4):
       shuffled = hand.sample(frac=1, random_state=seed).reset_index(drop=True)
       out, audit = _pipeline_after_read(shuffled, HAND, Audit(HAND)), ...
       out = out.sort_values("case_id").reset_index(drop=True)
       ...                                    # compare against base / base_md
   ```

   `HAND.path` need not exist: nothing after the read touches the filesystem, and the header prints
   `source.path.name` and `source.label`, never a path that is opened. That is the same property §9.3
   relies on when it forbids absolute paths in the log.

   The log this produces carries five of §9.6's entries and not the two the read records
   (`provenance/read`, `contract/rename_and_drop`), because the shuffle starts after the read. That
   is correct and is not a coverage hole: those two are byte-compared by §12.8 and their content is
   asserted by §12.16, and neither can depend on row order — one counts rows, the other counts
   columns.

   The seam is the whole pipeline after the read, not `_correct` alone: `_normalise` and
   `_missingness` also write into the log, and a row-indexed edit in either would otherwise pass
   green while §7.1 claims the property covers them. Left at `_correct`, an implementer could equally
   have run this against `load(FIXTURE)` — whose two rows exercise neither §7.2 nor §7.3, so the test
   would pass against an implementation that got both corrections wrong.
6. **Penumbra recomputation.** On a hand-built frame with one copy-paste record: the value is
   corrected, the identifier is named, the entry's `n` is 1, and no value is negative. Source
   missingness is preserved exactly — a record missing `core_ml` is missing `penumbra_ml` after, and
   no record gains a value. On a frame where every stored value already agrees, the entry has `n = 0`
   and no identifiers.
7. **Every correction names its cases.** Constructing an `AuditEntry(kind="correction", n=1,
   case_ids=())` raises `ValueError`. The same with `n=0` does not. An unknown `kind` raises.
   **[ENG]** `kind="observation"` with `n=1` and no identifiers raises the same way (§9.2), and
   `kind="provenance"` with `n=126` and no identifiers does not — parametrised over `KINDS` so the
   rule's boundary is asserted for all six, not only for the two that must name cases.
8. **Byte-identical reproduction.** Two `load()` calls on the same source render identical markdown;
   two `write()` calls produce identical bytes. The rendered text matches neither `\d{4}-\d{2}-\d{2}`
   nor `\d{2}:\d{2}:\d{2}`, so the property cannot be passing because both runs fell in the same
   second. The rendered text contains no absolute path.

   Three additions **[ENG]**:
   - **Across processes, not just within one.** Two interpreters launched with `PYTHONHASHSEED=0` and
     `PYTHONHASHSEED=1` render the same log for the same source, byte for byte. A same-process
     comparison cannot see frozenset ordering (§5), which is the one way this has actually broken.
   - **`load()` writes nothing.** `load(FIXTURE)` leaves `C.LOGS` unchanged: if the directory did not
     exist it still does not, and if it did, no file in it was created or modified.
   - **The two sources cannot collide.** `_audit_path(WORKBOOK) != _audit_path(FIXTURE)`, and
     `write()` creates `C.LOGS` when absent and returns the path it wrote.
9. **Schema assertions fire.** Parametrised, one targeted corruption of a valid frame per case:
   `age = 12`; `nihss_baseline = 50`; `mrs_90d = 7`; `prestroke_mrs = 6`; `core_ml = -1`; `ivt = 2`;
   `ivt` missing; **`ivt_contraindicated` missing (A4b) [ENG]**; a duplicated `case_id`; a padded
   `case_id`; a `contraindication_reason` of `"  "`; `onset_to_ivt_min` present for a control;
   `onset_to_ivt_min` absent for a treated patient; `onset_to_ivt_min > onset_to_groin_min`; a
   `center` value outside `CENTER_RECODE`. Each raises `SchemaError`, and the message names the
   column and the offending value.

   **Plus one test that two corruptions produce one error naming both. [ENG]** A frame carrying both
   `mrs_90d = 7` and `core_ml = -1` raises a single `SchemaError` whose message names `mrs_90d`,
   `core_ml` and the count of failed assertions. Without it, §8.1's collect-all rule is untested and
   a fail-fast implementation passes every other case in this list.

   A1 is not in this list, and §8.1 says why: nothing between §5 and §8 changes the row count, so it
   is a post-condition with no reachable failure through `load()`.
10. **No bare `assert` in `data.py`.** An AST scan for `ast.Assert` over `data.py`, because `python
    -O` strips them (§8.3). A companion test parses a synthetic snippet containing an `assert` and
    confirms the scan fires, so a scan that silently matches nothing cannot pass as a green test.
11. **Structural non-applicability [§11].** The missingness table labels `onset_to_ivt_min`
    `structural` and `contraindication_reason` `not recorded`, and neither is counted as `missing`.
    Every key of `STRUCTURALLY_NON_APPLICABLE` and of `INFORMATIVE_ABSENCE` is an analysis name and
    appears in the table with its reason printed. Data-gated: the structural count for
    `onset_to_ivt_min` equals the number of controls.
12. **No imputation.** For every column of the returned frame, the missingness pattern equals that of
    the raw frame after rename and drop — `penumbra_ml` excepted, where it equals
    `core_ml.isna() | tmax6_ml.isna()`. Asserted per column, not in aggregate: an aggregate count can
    net a gain against a loss.
13. **Denominators.** The missingness table has one row per analysis column, its `n + n_absent`
    reconciles to the frame's row count on every row, and the per-centre absence counts sum to the
    overall count on every row.
14. **Structural facts from the workbook.** Data-gated, against the real file: full read is 126 × 42
    and a header-only read is 41; arms are 39 treated / 87 control; exactly one record is corrected
    by §7.2 and its discrepancy is 83 mL; exactly one record has `onset_to_groin_min == 999`; no
    treated patient carries `ivt_contraindicated == 1`; `ivt_contraindicated` is never missing (A4b);
    all of A1–A9 including A4b pass. These are the facts §7 and §8 were written against; the hash
    makes them safe, and the test documents the dependency so updating the hash re-checks the
    premises.
15. **Stage 1 amendments. [ENG]** These live in `test_config.py`, not `test_data.py`, because they
    test Stage 1 constants and the Stage 1 fixture. Every key of `INFORMATIVE_ABSENCE` resolves under
    9.13's union; `STRUCTURALLY_NON_APPLICABLE` and `INFORMATIVE_ABSENCE` are disjoint; the fixture
    holds the literal three-character string `N/A` in `TICI_2b_3` on its first row when read with
    `dtype=object, keep_default_na=False`, and reads as `<NA>` under `READ_DTYPES` and `NA_VALUES`.
    The last of these is what makes §12.2 runnable without `data/`, and it is why the fixture is
    regenerated: `test_config.py`'s own header calls the fixture's purpose reproducing "the two
    structural hazards of the real workbook", and §4.2 established a third.
16. **The audit log's inventory. [ENG]** The rendered log of a `load(FIXTURE)` contains exactly the
    §9.6 entries: the five unconditional ones, in the declared order, with the declared `kind` and
    `step` strings; no `onset_to_groin_999` entry, because the fixture has no such record; one
    structural entry per key of `STRUCTURALLY_NON_APPLICABLE`. `penumbra_recomputed` is present with
    `n = 0` and no identifiers. Asserted against the §9.6 table's `(kind, step)` pairs, so an entry
    cannot be dropped, renamed or reordered without this test failing.

    Three properties of the rendering go with it: every `KINDS` heading appears, the `Observations`
    section renders the single line `_none_` because the fixture triggers no observation (§9.5) —
    **[AMENDED 2026-08-10]** its extent runs to the heading of whichever kind follows `observation`
    in `KINDS`, which is now `derivation`, and the test reads that neighbour from `KINDS`; the
    header's `Read:` line agrees digit for digit with the `provenance/read` entry's `n` and the
    frame's raw column count; and the header's `mapped, dropped` counts agree with
    `contract/rename_and_drop`. A header recomputed independently of the entries is the one way these
    can disagree, and §9.5 forbids it.

### Coverage map

```
data.py                                            test_data.py
├── Source / WORKBOOK / FIXTURE / SOURCES          ├── 12.1  two sources, frozen, fields set
│   └── test-local Source is legal          [ENG]  └── 12.1  used by 12.1, 12.5
│
├── _verify_version(source)
│   ├── hash matches            → None             ├── 12.1
│   ├── hash differs            → DataVersionError ├── 12.1  wrong-hash Source over the
│   │                                              │         fixture — no data/ needed  [ENG]
│   └── sha256 is None          → skipped          └── 12.1  fixture path
│
├── _read(source, audit)
│   ├── full parse, READ_DTYPES → 42 cols          ├── 12.2
│   ├── 'N/A' cell              → <NA>, Int64 kept ├── 12.2  [REV] §4.2, needs 12.15  [ENG]
│   └── header-only read        → SchemaError      └── 12.2  cause named
│
├── _apply_contract(raw, source, audit)
│   ├── contract mismatch       → SchemaError      ├── 12.2
│   ├── row count mismatch      → SchemaError      ├── 12.2
│   ├── n_records is None       → no raise         ├── 12.2
│   ├── drop list in contract order         [ENG]  ├── 12.8  cross-process render
│   └── rename + drop           → ANALYSIS_NAMES   └── 12.3
│
├── _pipeline_after_read(df, source, audit) [ENG]  ├── 12.5  the seam the shuffle drives
│   │
│   ├── _normalise(df, audit)
│   │   ├── centre recode       → 4 labels         ├── 12.4
│   │   ├── unmapped code       → SchemaError      ├── 12.4
│   │   └── center is not Categorical              └── 12.4
│   │
│   ├── _correct(df, audit)
│   │   ├── penumbra recomputed → 1 case named     ├── 12.6
│   │   ├── penumbra already agrees → n=0, no ids  ├── 12.6
│   │   ├── groin 999           → observation      ├── 12.6 / 12.14
│   │   ├── row order irrelevant                   ├── 12.5  shuffle, 5 seeds, 6 records
│   │   └── missingness unchanged                  └── 12.12
│   │
│   ├── _assert_schema(df, source)  A1…A9 + A4b
│   │   ├── A1 row count → post-condition   [ENG]  ├── §8.1  unreachable by design,
│   │   │                                          │         deliberately not in 12.9
│   │   ├── A4b ivt_contraindicated missing [ENG]  ├── 12.9
│   │   ├── every other check              → raise ├── 12.9  one corruption each
│   │   ├── two violations → one error      [ENG]  ├── 12.9  collect-all is tested
│   │   └── no ast.Assert anywhere in data.py      └── 12.10 [REV] §8.3
│   │
│   └── _missingness(df, audit)
│       ├── structural entry per declared column   ├── 12.11
│       └── absence_by_column(…, ANALYSIS_NAMES)   │        [AMENDED 2026-08-10] public;
│           ├── kind: structural / not recorded /… ├── 12.11  Stage 3 calls it again on
│           ├── centre columns from CENTER_ORDER   ├── 12.13  its derived columns, so the
│           └── per-centre columns reconcile       └── 12.13  classification has one home
│
└── Audit
    ├── created before _read, threaded      [ENG]  ├── 12.16 header counts match entries
    ├── correction with n>0 and no ids → ValueError├── 12.7
    ├── observation with n>0 and no ids     [ENG]  ├── 12.7  → ValueError, same rule
    ├── provenance with n>0 and no ids      [ENG]  ├── 12.7  legal, no raise
    ├── unknown kind                   → ValueError├── 12.7
    ├── the §9.6 entry inventory            [ENG]  ├── 12.16
    ├── an empty section prints `_none_`    [ENG]  ├── 12.16 Observations, on FIXTURE
    ├── to_markdown deterministic, same process    ├── 12.8
    ├── to_markdown deterministic, across seeds    ├── 12.8  PYTHONHASHSEED 0 vs 1  [ENG]
    ├── load() writes no file               [ENG]  ├── 12.8
    ├── _audit_path differs per source      [ENG]  ├── 12.8
    ├── write() creates C.LOGS when absent  [ENG]  ├── 12.8
    └── write() bytes deterministic                └── 12.8

config.py / test_config.py                         test_config.py
├── INFORMATIVE_ABSENCE keys resolve        [ENG]  ├── 12.15 (9.13's union)
├── the two absence kinds are disjoint      [ENG]  ├── 12.15
├── fixture holds a literal 'N/A'           [ENG]  ├── 12.15
├── data.py is scanned for raw names               ├── 9.4   already green today
└── data.py is NOT in the exemption list            └── 9.4   already green today

Data-gated (skipif not DATA_XLSX.exists()):        └── 12.14 126x42 / 41, 39/87 arms,
                                                            1 penumbra fix at 83 mL,
                                                            1 groin 999, 0 treated flagged,
                                                            A4b holds
Every branch above has a test. No branch is untested.
```

**Two of the guards §12's blockquote asks for already exist and are green.** `test_config.py`'s 9.4
walks `extended_bridging/**/*.py` minus the exemptions, so `data.py` and `test_data.py` come under
the raw-name scan the moment they are created — no change needed. And 9.4's
`test_the_exemption_list_is_exactly_the_three_declared_files` asserts
`set(EXEMPT_FROM_RAW_NAME_SCAN) == {"config.py", "stage0_data_inventory.py", "test_config.py"}`, so
adding `data.py` to the exemptions fails a Stage 1 test rather than merely violating a convention.
Definition of done #3 is enforced, not requested.

### Failure modes

| Failure | Test | Error handling | Visible or silent |
|---|---|---|---|
| Workbook replaced | 12.1 (data-gated) | `DataVersionError` names both hashes and the fix | visible |
| Column added, removed or renamed | 12.2 | `SchemaError` via `assert_column_contract` | visible |
| Header-only read used by mistake | 12.2 | `SchemaError` names the cause | visible |
| Rows filtered upstream | 12.2 | `SchemaError` names both counts | visible |
| New centre code | 12.4 | `SchemaError` names the value | visible |
| Out-of-range score, volume or time | 12.9 | one `SchemaError` listing every violation | visible |
| Non-binary value in a 0/1 column | 12.9 | same, column and values named | visible |
| Treatment missing | 12.9 | same | visible |
| Eligibility classifier missing **[ENG]** | 12.9 | `SchemaError` from A4b, before Stage 4 sees `<NA>` | visible |
| Several columns drift at once **[ENG]** | 12.9 | one error names all of them; no fix-rerun cycle | visible |
| IVT time present for a control | 12.9 | `SchemaError` — the arm label is wrong | visible |
| Correction added without identifiers | 12.7 | `ValueError` from `AuditEntry` | visible |
| Standing query logged without identifiers **[ENG]** | 12.7 | `ValueError` — an observation nobody can look up is not a query | visible |
| Correction written by row position | 12.5 | shuffle test diverges over the whole post-read pipeline | visible |
| The log's header disagrees with its own entries **[ENG]** | 12.16 | unrepresentable: the header is rendered from the entries (§9.5) | visible |
| A check written as `assert`, run under `-O` | 12.10 | AST scan fails at test time | visible |
| Log made non-reproducible (timestamp, path) | 12.8 | byte comparison and pattern scan fail | visible |
| Log made non-reproducible across processes **[ENG]** | 12.8 | two-seed render comparison fails | visible |
| A test run overwrites the workbook's audit log **[ENG]** | 12.8 | unrepresentable: `load()` writes nothing and the path carries the label | visible |
| An audit entry dropped or renamed **[ENG]** | 12.16 | inventory comparison against §9.6 fails | visible |
| A value imputed | 12.12 | per-column missingness comparison fails | visible |
| A literal `NA` arrives as a contraindication reason | none | reads as "no reason recorded" | **silent** — §13 |

One silent failure mode remains, and it is recorded rather than papered over. Every other row above
is visible at the point of failure, and six of them were silent in the 2026-08-08 draft. Three are
listed as *unrepresentable* rather than as caught: the shape of the code denies them, which is a
stronger guarantee than a test and is the standard §3, §9.2 and §9.5 each hold themselves to.

## 13. Known gaps carried forward

- **A literal `NA` in a text column would read as missing.** `keep_default_na=True` is required by
  §4.2, and its NA list includes `NA`. In `contraindication_reason` that would move a patient from
  *eligible* to *indeterminate* [DECISION 1]. No check catches it, because a reason of `NA` and a
  blank cell are indistinguishable after the read. Accepted: the current ten values are all ordinary
  clinical text, and the alternative — reading that column with `keep_default_na=False` while the
  rest of the sheet keeps it — is a second read path and a worse trade.
- **`onset_to_groin_min = 999` is unresolved** (§7.3). A standing query with the data owner. The
  identifier is in the audit log, not here.
- **The two records whose shipped death flag contradicts their mRS** remain a standing query
  [DECISION 2]. Stage 2 has nothing to do about them: the contract drops both shipped columns, so no
  correction is available or needed.
- **The Stage 0 note is inaccurate in two places** (§4.2): the `'N/A'` claim, and the "43 columns"
  header. Re-running `stage0_data_inventory.py` will not fix either — the first is a claim its own
  probe cannot see, the second is a count taken after a derived column is added. Both are worth
  amending when Stage 0 is next touched; neither blocks Stage 2. **[ENG]** Folding the fix into T1
  was considered at review and declined: T1 restores Stage 0's *runnability*, which is what Stage 2
  needs from it, and correcting the note's content means editing an accepted stage's script from
  Stage 2's branch. This bullet is the record; there is deliberately no `TODOS.md` entry, because a
  second file restating this one is a second thing to keep true.
- **`sex` coding is undocumented**, carried from Stage 1 §10. `SEX_LABELS` stays `None`.

## 14. What Stage 2 deliberately does not decide

Recorded so a later stage does not look here for an answer that was never placed here.

- **How dichotomies are built.** Stage 3, from `OUTCOMES` and `OPS`, with missingness reimposed.
- **How `onset_type` is built**, and the assertion that wake-up and unwitnessed never coincide.
  Stage 3 (§8.2).
- **How the [§13] subgroups are built.** Stage 3 — `unknown_onset` from `onset_type`, and
  `core_above_median` from a median that is a property of the *cohort* and so is computed only after
  Stage 5's restrictions. Stage 2 guarantees the `penumbra_ml` the balance table reads is the one
  [§6] defines (§7.2), and nothing more. **[AMENDED 2026-08-10]** This line used to read "whether
  `target_mismatch` is met"; that subgroup is withdrawn.
- **Which patients are eligible.** Stage 4. Stage 2 guarantees `ivt_contraindicated` is 0/1 (A5) and
  never missing (**A4b [ENG]**), and that no reason cell is blank-but-present (A9); the
  classification, and the check that no treated patient is flagged, are Stage 4's (§8.2). The
  2026-08-08 draft credited the never-missing half to A5, which permits missing everywhere — so the
  guarantee this section promised Stage 4 was not one Stage 2 made. A4b is that guarantee.
- **Which centres and patients are dropped.** Stage 5.
- **How factors are coded.** Stage 6, from `FACTOR_LEVELS` and `REFERENCE_LEVELS` (§6).

## 15. NOT in scope for Stage 2

| Considered | Why deferred |
|---|---|
| Deriving any outcome or covariate | Stage 3. The one recomputation here restores a definition [§6] already gives, and creates no new column. |
| Imputation of any kind | Never [§11]. §12.12 asserts it mechanically. |
| Reading the contraindication free text | DECISION 1. Only presence is read, and A9 protects that bit. |
| Normalising `contraindication_reason` spellings | Same. The pilot's regex is not lifted (§6). |
| Making `center` a `Categorical` | Stage 6, at the point of use, after the Stage 5 restrictions (§6). |
| A per-column missingness *threshold* | No column is dropped for missingness in this plan; [§11] is complete-case with reported denominators, and the denominators are §10.3. |
| Caching the parsed frame | 126 rows. The read is not the cost, and a cache is a second source of truth. |
| A CLI entry point | Stage 14 is the single entry point [§16]. `data.py` is a library module. |

## 16. What already exists, and what to lift

`pilots/data.py` is gitignored but present, and is the closest prior art. Lift these, checking each
against §7 and §8:

| From the pilot | Status |
|---|---|
| the `Audit` shape — ordered log of `(step, n, detail, case_ids)` | lift the shape; add `kind`, the names-its-cases rule for corrections **and observations** (§9.2), sorted identifiers, the §9.3 determinism rules, the §9.6 inventory, and the label-derived path (§9.5) |
| the pilot's habit of building the log at the end of the run | **do not lift.** The `Audit` is constructed first and threaded through the read, so §9.6's `provenance/read` and `contract/rename_and_drop` entries exist where the facts do and the header is rendered from them (§0, §9.5) |
| `to_markdown`'s entry format | lift; it reads well |
| the penumbra recomputation | lift the correction; rewrite the justification from [§6] and apply it to every record (§7.2) |
| the `_check` battery — row count, unique ids, arms, ranges, non-negative volumes, IVT-before-groin, IVT time tracks the flag | lift the **checks**; not the `assert` statement form (§8.3), and read the constants rather than the literals |
| `missingness()` — overall plus per-centre | lift; extend with the three-way `kind` column (§10.3) |
| the `999` sentinel correction | **do not lift.** §7.3, with the evidence. |
| `_rebuild`, `_eligibility`, the derived columns, the reason regex, the `hir` entry | **do not lift.** §7.4. |
| bare `assert` | **do not lift.** §8.3. |
| `DROP_COLUMNS` / `DROP_AFTER_RENAME` | **do not lift.** `DROPPED` is a derived view of the contract; a second list is exactly the drift Stage 1 §4.3 exists to prevent. |

`stage0_data_inventory.py:97` has the correct `derive` idiom (`.mask(source.isna())`) — that belongs
to Stage 3, not here, and Stage 1 §3.2 already specifies it.

## 17. Implementation tasks

Ordered. Each independently verifiable; one lane, so no parallelisation applies.

- [x] **T1 (P1)** — `pyproject.toml`: add `tabulate>=0.9`; `uv lock`. Verify: `uv run python
      stage0_data_inventory.py` completes, and running it twice produces identical bytes.
- [x] **T2 (P1)** — `config.py`: add `INFORMATIVE_ABSENCE` per §10.2, and correct the `NA_VALUES`
      comment per §4.2. `test_config.py`: add `INFORMATIVE_ABSENCE` to 9.13's union and add the
      disjointness test. **[ENG]** Verify: acceptance 12.15's first two clauses — and see them fail
      first, by temporarily misspelling a key, because 9.13 did **not** cover the new constant before
      this task and a green suite would otherwise prove nothing.
- [x] **T2b (P1) [ENG]** — `test_config.py`: put the literal `'N/A'` in `_FIXTURE_ROWS[0]`'s
      `tici_2b_3`; regenerate with `uv run python test_config.py --write-fixture`; commit
      `tests/fixture_schema.xlsx`. Verify: acceptance 12.15's third clause, and `test_config.py` green
      — in particular `test_the_fixture_reads_under_the_declared_dtypes`, whose
      `mRSscoreat90days.isna().sum() == 1` is why the cell goes in `TICI_2b_3` and not in the mRS
      column.
- [x] **T2c (P2) [ENG]** — `specs/stage1_config_and_data_contract.md` §5.1: correct "Stage 0 found no
      sentinels; this is a guard" per §4.2, so the two specs agree. Verify: read.
- [x] **T3 (P1)** — `data.py` §3: `Source`, `WORKBOOK`, `FIXTURE`, `SOURCES`, `_verify_version`, and
      the test-local-`Source` allowance. Verify: acceptance 12.1, including the wrong-hash case with
      `data/` renamed away.
- [x] **T4 (P1)** — `AuditEntry` with `_MUST_NAME_CASES`, `Audit`, `_md_table`, `_audit_path` per §9,
      and the §9.6 inventory. `Audit(source)` must be constructible before anything is read — T5 and
      T9 both pass it in. Verify: acceptance 12.7, 12.8, 12.16.
- [x] **T5 (P1)** — `_read` and `_apply_contract` per §4 and §5, drop list in contract order. Verify:
      acceptance 12.2, 12.3. Depends on T2b for the `'N/A'` clause of 12.2.
- [x] **T6 (P1)** — `_normalise` per §6. Verify: acceptance 12.4.
- [x] **T7 (P1)** — `_correct` per §7, both entries. Verify: acceptance 12.6, 12.12.
- [x] **T8 (P1)** — `_assert_schema` per §8: A1–A9 plus A4b, every violation collected into one
      `SchemaError`, no bare `assert`. Verify: acceptance 12.9 including the two-violations case,
      and 12.10.
- [x] **T9 (P1)** — `_missingness` per §10, `_pipeline_after_read` as the §0 seam, and `load` wiring
      the pipeline. Verify: acceptance 12.5, 12.11, 12.13.
- [x] **T10 (P2)** — `test_data.py` remainder, including the data-gated tests behind `DATA_GATED` and
      the cross-process render of 12.8. Verify: `uv run pytest -v` green both with and without
      `data/` present.

T1, T2, T2b and T2c touch Stage 0 and Stage 1 files and nothing in `data.py`; T3–T10 touch `data.py`
and `test_data.py` and nothing else. So there are two lanes and one dependency across them: T5's
`'N/A'` acceptance clause needs T2b's regenerated fixture. Everything else in the second lane is
sequential on the first — the pipeline is one module — so the parallelism is worth at most one
worktree and is not worth taking.

### Diagrams that belong in the code, not only here  **[ENG — new]**

`config.py` carries its reasoning in comments beside the declarations, and `data.py` should do the
same. Three diagrams, and no more — a diagram nobody maintains is worse than none, because it is
believed:

- **`data.py`'s module docstring** — the §0 pipeline, from `load` down to the seven audit entries,
  including the "writes nothing" arrow. It is the map for every later stage that imports this module,
  and it is the first thing a reader sees.
- **Above `_correct`** — the two corrections and their fates side by side: `penumbra_ml` recomputed
  on every record and logged as a `correction`, `onset_to_groin_min == 999` left standing and logged
  as an `observation`. The asymmetry is §7's whole argument and is the thing most likely to be
  "tidied" by someone who has not read §7.3.
- **Above `_missingness`** — the three-way split of absence: `structural` (a key of
  `STRUCTURALLY_NON_APPLICABLE`), `not recorded` (a key of `INFORMATIVE_ABSENCE`), `missing`
  (anything else), `complete` (none). A reader who does not see the split reads 82/126 on
  `contraindication_reason` as two-thirds unusable, which is exactly the misreading §10.2 exists to
  stop. **[AMENDED 2026-08-10]** The classification now lives in the public `absence_by_column`,
  which Stage 3 calls a second time over its derived columns; the comment belongs with it, because a
  second implementation over there would drift on precisely the `structural`-versus-`missing` branch
  this note exists to protect.

**Keeping them true is part of any change that touches them.** If a later stage edits `_correct` or
`_missingness`, the diagram above it is edited in the same commit or the change is not finished. This
applies to §0's diagram in this file as well: it is the one place the pipeline's shape is written
down, and a stale one would send an implementer down the 2026-08-08 draft's path.

### Definition of done

Stage 2 is complete when all of the following hold, and not before:

1. `uv run pytest -v` is green **with** `data/` present — every test, including the data-gated ones.
2. `uv run pytest -v` is green **with `data/` temporarily renamed** — the data-gated tests skip and
   nothing else fails or errors at collection. `FIXTURE` carries the whole pipeline.
3. `test_config.py` is still green. The Stage 1 §7 raw-name scan now covers `data.py` and
   `test_data.py`, and neither has been added to `EXEMPT_FROM_RAW_NAME_SCAN`.
4. The workbook's audit log has been produced and reproduces byte-identically, **[ENG]** including
   across interpreters with different hash seeds:

   ```bash
   cd extended_bridging
   PYTHONHASHSEED=0 uv run python -c \
     "import data; _, a = data.load(); print(a.write())"          # out/logs/audit_v7_…md
   PYTHONHASHSEED=1 uv run python -c \
     "import data; _, a = data.load(); a.write(Path('/tmp/b.md'))"
   diff out/logs/audit_v7_july26_with_abs_contra_indication.md /tmp/b.md
   git check-ignore -v out/logs/audit_v7_july26_with_abs_contra_indication.md
   ```

   `diff` reports nothing and `check-ignore` reports a match — the log names patients and must never
   be committed. `load()` itself wrote no file; the caller asked for one.
5. Acceptance 12.10 has been **seen to fail** when a bare `assert True` is pasted into `data.py`, and
   to pass again when it is removed. A green scan that has never failed is not evidence.
6. Acceptance 12.5 has been **seen to fail** against a deliberately row-indexed correction
   (`df.loc[df.index[0], "penumbra_ml"] = …`) and to pass when it is rewritten as a mask. Same
   reasoning.
7. **[ENG]** Acceptance 12.15's disjointness and 9.13 clauses have been **seen to fail** against a
   misspelled `INFORMATIVE_ABSENCE` key, and acceptance 12.9's two-violations case against a
   fail-fast `_assert_schema`. Three tests, all of which pass trivially on the wrong implementation
   if they are never watched failing.
8. The audit log has been read end to end by a human, and every correction in it names its cases.
9. **[ENG]** `uv run pytest -v` was run once with `data/` present *after* a run with `data/` renamed
   away, and `out/logs/` contains no `audit_fixture_schema.md` — confirming the fixture path
   is derived and that no test wrote through `load()`.

## 18. Verification record

Checked against the workbook on 2026-08-08, before this spec was finalised, by a read-only probe
reporting aggregate counts only. Recorded so a reader can tell which numbers were verified rather
than carried from the Stage 0 note or the pilot, and so the checks are re-runnable after a workbook
update. Acceptance 12.14 re-checks the structural subset in code.

| Claim | Where used | Verified |
|---|---|---|
| `read_excel` with `READ_DTYPES` and `NA_VALUES` succeeds; full read 126 × 42, header-only 41 | §4.1, §5 | yes |
| Header-only column set equals the contract minus `Unnamed: 41` | §4.1 | yes |
| Literal `'N/A'` present: `mRSscoreat90days` 2, `TICI_2b_3` 3, six 24-h NIHSS columns 1 each | §4.2 | yes |
| `Center` takes exactly the four keys of `CENTER_RECODE` under `dtype="string"`; recode is total | §6, A3 | yes |
| `case_id` unique, non-missing, no surrounding whitespace, `string` dtype | §6, A2 | yes |
| `penumbra_ml` disagrees with `tmax6_ml − core_ml` on exactly 1 record; discrepancy 83 mL; the stored penumbra equals the stored core exactly | §7.2 | yes |
| ~~That record's criteria: core 13 < 70; ratio 109/13 = 8.4 > 1.8; penumbra 13 → 96 crosses `min_penumbra_ml` = 15~~ | ~~§7.2~~ | **[WITHDRAWN 2026-08-10]** — verified when checked, but the [§13] subgroup whose criteria these are no longer exists, so the row states a fact about nothing. §7.2's surviving claim is the row above it |
| `core_ml`, `tmax6_ml` and `penumbra_ml` are missing on the same single record | §7.2, §10.1 | yes |
| `onset_to_groin_min == 999` on exactly 1 record; its `onset_to_ivt_min` is 955, so the interval is 44 min | §7.3 | yes |
| IVT-to-groin interval across the 39 treated: min 29, median 70, max 310 | §7.3 | yes |
| `onset_to_groin_min` range 253–1219, with 905, 905, 905, 950, 999, 1035, 1045 in 900–1100 | §7.3 | yes |
| `onset_to_ivt_min` non-missing for exactly the 39 treated, missing for all 87 controls | A7, §10.2 | yes |
| `onset_to_ivt_min <= onset_to_groin_min` on all 39 | A8 | yes |
| All 13 `BINARY_COLUMNS` take values in {0, 1}; only `tici_2b_3` has missing (3) | A5 | yes |
| All 9 `PLAUSIBLE_RANGES` columns inside range: age 28–97, NIHSS 2–32, prestroke mRS 0–4, mRS 0–6, volumes and times ≥ 0 | A6 | yes |
| `contraindication_reason`: 44 non-missing, 10 distinct, none padded, none blank after stripping | §6, A9 | yes |
| `ivt_contraindicated` 0/1, never missing; 19 controls flagged; **0 treated flagged** | A5, §12.14 | yes |
| `wake_up` and `unwitnessed` never both 1, neither ever missing | §8.2 | yes |
| `tabulate` absent from `pyproject.toml` and `uv.lock`; `DataFrame.to_markdown()` raises `ImportError` | §2.2 | yes |
| `tests/fixture_schema.xlsx` uses sheet `Feuil1`, passes `assert_column_contract`, renames to exactly `ANALYSIS_NAMES`, and its two rows satisfy A5–A9 | §3, §12 | yes |
| The Stage 0 note's "43 columns" counts a derived `center` column added after the read | §4.2 | yes |

**[ENG] Claims added by the 2026-08-10 revision, and how each was checked.** These were verified
against the committed repository rather than against the workbook, so they are re-checkable by anyone
with a checkout — no `data/` needed.

| Claim | Where used | Verified |
|---|---|---|
| `test_config.py` 9.13 parametrises over `PS_COVARIATES_FULL`, `BALANCE_ONLY`, `BINARY_COLUMNS`, `PLAUSIBLE_RANGES`, `STRUCTURALLY_NON_APPLICABLE`, `POST_TIME_ZERO` — and nothing else | §10.2, §12.15 | yes, read |
| `tests/fixture_schema.xlsx` as committed contains no `'N/A'` cell: `_FIXTURE_ROWS` holds only ints, floats, `None`, `"FIX-001"`, `7001`, `1`, `"Lausanne"`, `"Clinician decision"` | §12.2, §12.15 | yes, read |
| `_FIXTURE_ROWS[0]["tici_2b_3"]` is not read by any existing Stage 1 test; `test_the_fixture_reads_under_the_declared_dtypes` asserts only on `CaseID`, `Age` and `mRSscoreat90days` | T2b | yes, read |
| `config.py`'s `NA_VALUES` line still carries the comment "Stage 0 found no sentinels; this is a guard" | §4.2, §1 | yes, read |
| `DROPPED` is a `frozenset`, and `CENTER_ORDER` a declared tuple equal to `CENTER_RECODE`'s values | §5, §10.3 | yes, read |
| `test_config.py` 9.4 walks `extended_bridging/**/*.py` minus the three exemptions, and pins the exemption set against that literal — so `data.py` is scanned, and exempting it fails a Stage 1 test | §12, DoD 3 | yes, read |
| `stage0_data_inventory.py` writes no timestamp, date or clock time — so T1's byte-identical claim is achievable once `tabulate` is installed | §2.2, T1 | yes, read |
| `A5` as written permits missing in every `BINARY_COLUMNS` column, so it does not deliver §14's never-missing guarantee | A4b, §14 | yes, read |

## 19. What the 2026-08-10 revision changed

Nineteen entries in two tables. The first eleven each close a specific defect in the 2026-08-08
draft. The remaining eight close a different kind of gap: reasoning that was settled during the
review and would otherwise have lived only in the conversation that produced it.

Recorded so that a reader who has the earlier version can find what moved, and so that none of this
is re-litigated as a fresh idea in Stage 3.

### Defects in the 2026-08-08 draft

| # | Change | Sections | The defect it closes |
|---|---|---|---|
| 1 | `load()` writes no file; the log path is derived from `source.label` | §0, §9.2, §9.5, §12.8, DoD 4, 9 | Every `pytest` run against `FIXTURE` overwrote the workbook's audit log, on any machine with `data/`, with no diff and no error — the file is gitignored |
| 2 | The audit log's entry inventory is pinned: seven entries, exact `kind`/`step`/`n`/`detail` | §9.6, §12.16 | Two of six declared `KINDS` had no entry specified and no `step`, `detail` or `n` was pinned, so byte-identity held per-implementation rather than per-specification |
| 3 | Every assertion runs; one `SchemaError` reports all violations; each check guards its own preconditions | §8.1, §12.9, failure modes | Fail-fast turned "three columns drifted in v8" into three read-parse-fail round trips, and A7 silently depended on A4 having passed |
| 4 | A4b: `ivt_contraindicated` is never missing | §8.1, §11, §12.9, §12.14, §14 | §14 promised Stage 4 a guarantee A5 does not make — A5 permits missing everywhere — and Stage 4's whole [DECISION 1] split rests on it |
| 5 | `config.py`'s `NA_VALUES` comment and Stage 1 spec §5.1 are corrected | §1, §4.2, T2, T2c | §4.2 proved the comment false and left it sitting on the line a future hardening pass would edit, having predicted that exact repair |
| 6 | The drop list and the missingness header range over declared order, not over a `frozenset` or a literal | §5, §10.3, §12.8 | `str` hashing is per-process randomised, so the `contract` entry's table would differ between runs while passing a same-process byte comparison |
| 7 | A1 is labelled an unreachable post-condition and kept | §8.1, §12.9 | It read as one of nine live checks while being untestable through `load()`, so its absence from §12.9 looked like a coverage gap |
| 8 | The fixture gains a literal `'N/A'` cell in `TICI_2b_3` | §1, §12.2, §12.15, T2b | §12.2 asserted a behaviour the committed fixture could not exercise, so the most consequential `[REV]` in this spec had no test that runs without `data/` |
| 9 | `INFORMATIVE_ABSENCE` joins 9.13's union, and the two absence kinds must be disjoint | §1, §10.2, §12.15, T2 | T2's verification step was a no-op: 9.13 does not reach the new constant, so a typo in the key would surface only as a misclassified column in §10.3 |
| 10 | §12.5's harness is specified: `_pipeline_after_read`, a test-local `Source`, six records, five seeds | §0, §3, §12.5 | The spec's strongest test named no function to call and no `Source` to build, and §3 read as forbidding a test-local one; the plausible readings either fail §12.1 or test nothing |
| 11 | Four editorial corrections: wrong-hash test runs without `data/`; `write()` creates `C.LOGS`; `_apply_contract`'s signature; §8.2 gains the treated/contraindicated invariant | §5, §8.2, §12.1, §12.8 | Small drifts between the prose, the coverage map and `config.py`'s statement of who owns which assertion |

### Settled in the same review, and written down here rather than left in conversation

These are not corrections to the draft. They are answers that existed only as reasoning during the
review, and each one is a question the next implementer would otherwise have had to re-answer —
differently, and with no record of why.

| # | Recorded | Sections | The question it stops being asked |
|---|---|---|---|
| 12 | The `Audit` is created at the top of `load()` and threaded through `_read` and `_apply_contract`; the header's counts are read back from those entries | §0, §4.1, §5, §9.5 | §9.6 puts two entries inside functions that ran before any `Audit` existed in the draft's diagram. Threading it also gives every number in the header exactly one origin |
| 13 | `observation` is held to the correction-names-its-cases rule; `_MUST_NAME_CASES` names both kinds | §9.2 | Roadmap Stage 2 demands it only of corrections, so an implementer would reasonably have left an unattributed observation — and §7.3's standing query is unanswerable without the identifier |
| 14 | Why three new types, and what collapsing each would cost | §0.1 | A reviewer counts `Source`, `AuditEntry`, `Audit` and asks. The answer was given once, in conversation, and would have been asked again at every future review |
| 15 | What Stage 2 costs to run, and the standing instruction not to cache | §2.3 | §15 rejects a parsed-frame cache on principle. The number behind the principle — 126 rows, two file reads, sub-second — was nowhere, so the rejection read as taste |
| 16 | A record with a missing volume is never counted as disagreeing; `_TOL` is absolute and why | §7.2 | `n` in the penumbra correction entry depends on it, and a reader who assumed missing counted as disagreement would read v7's `n = 1` as a bug |
| 17 | An empty section prints its heading and `_none_` | §9.5 | The headings are a fixed skeleton; without the rule, a clean workbook's log looks truncated and an implementer may drop the heading instead |
| 18 | Why the Stage 0 note fix is not folded into T1, and why there is no `TODOS.md` | §13 | The option was raised and declined at review. Without the record, the next reader sees a known-false note, an obviously adjacent task, and no reason not to |
| 19 | Which three ASCII diagrams go into `data.py`, and that keeping them true is part of any change touching them | §17 | `config.py` carries its reasoning beside its declarations; without saying so, `data.py`'s would live only in this file, which nobody has open while editing `_correct` |

Nothing in this revision changes a statistical decision, a column's fate, or a number. Every change
is about whether the specification says what it means, whether a test can tell, and whether the
reasoning survives the conversation it was produced in.
