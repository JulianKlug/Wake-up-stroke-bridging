# Stage 1 spec — configuration and data contract

Implements roadmap Stage 1. Section references in brackets are to `statistical_analysis_plan.md`.
Numbers and decisions referenced as DECISION *n* are established in Stage 0 and recorded in
`../out/stage0_data_inventory.md`.

**Status.** Revised 2026-08-07 after engineering review. This document is the **sole source for the
Stage 1 implementation**: everything the implementer needs is here, and anything not here is not to
be invented. Where the review changed a decision, the change is marked **[REV]** with its reasoning,
so a reader of the pilot or of the pre-review draft can see what moved and why.

**Goal.** One module that holds every fact about the source workbook and every prespecified constant,
so that no later module ever names a raw column, and so that a change to the workbook fails loudly
here rather than propagating silently into an estimate.

**Not in scope.** Reading or cleaning data (Stage 2), deriving variables (Stage 3), classifying
eligibility (Stage 4). This stage declares; it does not compute. Two exceptions, both pure functions
of their arguments and neither touching the filesystem: the column contract check, and the outcome
model covariate lookup. Both are called by later stages.

---

## 0. Where Stage 1 sits

```
  data/…v7…with_abs_contra_indication.xlsx        (gitignored, patient data)
                     │
                     │  read once, full parse
                     ▼
  ┌──────────────────────────────────────────────────────────────┐
  │  STAGE 1 — config.py        (stdlib only, no pandas import)  │
  │                                                              │
  │   COLUMN_CONTRACT ──┬─→ RENAME          (name is not None)   │
  │   42 verbatim keys  ├─→ READ_DTYPES     (dtype is not None)  │
  │                     └─→ DROPPED         (name is None)       │
  │                                                              │
  │   OUTCOMES ─────────┬─→ POST_TIME_ZERO  (denylist, inv. 4)   │
  │   8 entries         └─→ Outcome.rule    (derived from op/thr)│
  │                                                              │
  │   PS_COVARIATES ────┬─→ OUTCOME_COVARIATES   (same object)   │
  │   (tuple, frozen)   ├─→ STANDARDISATION_COVARIATES (derived) │
  │                     └─→ PS_COVARIATES_FULL       (derived)   │
  │                                                              │
  │   assert_column_contract()      outcome_model_covariates()   │
  │   DATA_SHA256, SEED, N_BOOT, thresholds, factor levels …     │
  └──────────────────────────────────────────────────────────────┘
       │            │            │            │            │
       ▼            ▼            ▼            ▼            ▼
   Stage 2      Stage 3      Stage 4      Stages 6-13   Stage 14
   read+clean   derive       eligibility  estimation    reporting
```

Every arrow out of Stage 1 carries analysis names only. No arrow carries a raw header. That property
is what §7 tests.

## 1. Deliverables

| Path | Contents |
|---|---|
| `extended_bridging/pyproject.toml` | project metadata, dependencies, pytest config |
| `extended_bridging/uv.lock` | committed lockfile |
| `extended_bridging/config.py` | the configuration module — everything below |
| `extended_bridging/test_config.py` | the acceptance tests in §9 |
| `extended_bridging/tests/fixture_schema.xlsx` | **[REV]** synthetic 42-header workbook, two invented rows, no patient data |

`.venv/` is gitignored and is not a deliverable. `data/` and `out/` are gitignored and are not
deliverables; nothing in this stage writes to either.

**[REV] Why the fixture exists.** `data/` is gitignored, so without a fixture the whole suite —
including the dozen tests that are pure logic — fails on any checkout that lacks the private
workbook. The fixture carries the 42 verbatim headers (including the three with trailing whitespace
and the headerless final column) and two rows of invented values, so the contract check, the
trailing-whitespace case and the header-only-read hazard are all testable anywhere. §9.1 pins a
second test that compares fixture headers against the real workbook whenever it is present, so the
fixture cannot drift unnoticed.

## 2. Environment

`uv`, Python 3.12 (`requires-python = ">=3.12,<3.13"`, matching the pilot so its outputs stay
comparable). Runtime dependencies: `pandas>=2.2`, `numpy>=1.26,<2.1`, `scipy>=1.13`,
`statsmodels>=0.14.2`, `matplotlib>=3.9`, `openpyxl>=3.1`. Dev group: `pytest>=8.0`.

**[REV] The numpy upper bound is not optional.** The pilot pins `numpy>=1.26,<2.1`. numpy 2.x changed
scalar type promotion (NEP 50), which can move floating-point results at the last digits — enough to
perturb a bootstrap percentile. "Matching the pilot so its outputs stay comparable" requires matching
this pin, not only the Python version.

No `scikit-learn`, no `seaborn`, no `causallib`/`causalinference`: the Firth fit (Stage 6), the
weighted proportional-odds fit (Stage 8) and the augmented estimator (Stage 9) are all implemented
directly, because none of those libraries supports the observation weights the plan requires.
`statsmodels` is retained for unpenalised cross-checks in tests, not for any reported estimate.

`pyproject.toml` carries no `[build-system]`: this is an analysis project, not a distributable
package, so `uv` treats it as a virtual project and `uv sync` installs dependencies without trying to
build anything. `[tool.pytest.ini_options] testpaths = ["."]`, matching the pilot.

```bash
cd extended_bridging && uv sync && uv run pytest -v
```

The machine's system Python is 3.9; `uv` downloads and pins 3.12 for this project on first `uv sync`.
Nothing here should be run with the system interpreter.

**Working directory contract.** All commands are run from `extended_bridging/`. Flat module layout
(`import config`), matching the pilot. This is a convenience for the test runner only — §5 anchors
every path to `__file__`, so no constant in `config.py` depends on the working directory.

**Check the fixture is committable.** The repository `.gitignore` excludes `data/`, `out/`, `pilots`,
`.venv/`, `__pycache__/` and `.pytest_cache/`. `extended_bridging/tests/fixture_schema.xlsx` matches
none of those and commits normally — confirm with `git check-ignore -v` before assuming it landed.
It contains no patient data by construction (§1); that is what makes committing it safe.

## 3. Module shape

Module-level constants in `config.py`, matching the pilot's idiom.

**`config.py` imports only the standard library** — `dataclasses`, `pathlib`, `operator`, `typing`.
No pandas. Two consequences worth having: importing the configuration is free, and the pure-logic
acceptance tests in §9 run without a scientific stack. Where a value exists to be handed to pandas
(`READ_DTYPES`, `FACTOR_LEVELS`), it is a plain string or tuple that the consuming stage converts.

### 3.1 Exceptions

```python
class SchemaError(Exception):
    """The workbook's columns do not match COLUMN_CONTRACT."""

class DataVersionError(Exception):
    """The workbook's content hash does not match DATA_SHA256."""
```

**[REV]** The pre-review draft named `SchemaError` in §4 without ever defining it, and had no
exception for a content change at all. Both are declared here and nowhere else.

### 3.2 Structured types

Both frozen dataclasses, because their fields are asserted against and a bare tuple would not survive
review.

```python
@dataclass(frozen=True)
class Column:
    name: str | None      # analysis name, or None if dropped
    reason: str           # why it is dropped, or what it is. Never empty, ever.
    dtype: str | None = None   # pandas dtype passed to read_excel; None means infer

@dataclass(frozen=True)
class Outcome:
    label: str            # for tables and figures
    kind: str             # "ordinal" | "binary"
    family: str           # "primary" | "secondary" | "safety"
    source: str | None    # analysis name of the ordinal source, or None if read directly
    op: str | None        # "<=" | ">=" | "==" ; None iff source is None
    threshold: int | None # None iff source is None
    higher_is_better: bool

    @property
    def rule(self) -> str | None:
        """Human-readable derivation, e.g. 'mrs_90d <= 2'. Derived, never stored."""
        if self.source is None:
            return None
        return f"{self.source} {self.op} {self.threshold}"

OPS = {"<=": operator.le, ">=": operator.ge, "==": operator.eq}
```

**[REV] `Column.dtype`.** The pre-review draft recorded "string, not numeric" for `CaseID` as prose
in a markdown Note column. Prose is not callable, so Stage 2 would have written
`dtype={"CaseID": str}` — naming a raw header, which §7 forbids and §9.4 would then flag. The
contract now carries the dtype, and `READ_DTYPES` (§5) is derived from it, so Stage 2 passes a
config-owned mapping to `read_excel` and never types a header.

**[REV] `Outcome.op` / `Outcome.threshold`, with `rule` derived.** The pre-review draft stored `rule`
as a human-readable string. A string cannot be applied, so Stage 3 would have re-typed each threshold
as code — a second copy, with nothing asserting the two agree. An edit to one and not the other would
have the manuscript print `mRS 0-1` while the analysis computed `mRS 0-2`, and no test would fail.
The threshold is now stored once, as data; Stage 3 applies it through `OPS` and the prose is computed
from the same fields. This is the same derived-not-duplicated discipline invariant 3 applies to
covariate lists, extended to outcome definitions, where the consequence of drift is larger.

Stage 3 derives all four dichotomies in one loop, with no threshold written anywhere but the registry:

```python
for key, o in OUTCOMES.items():
    if o.source is None:
        continue                                    # read directly, nothing to derive
    src = df[o.source]
    df[key] = OPS[o.op](src, o.threshold).astype("Int64").mask(src.isna())
```

The `.mask(src.isna())` is not optional and is not decoration: a comparison against a missing value
returns `False`, so without it a patient with no 90-day mRS silently becomes a non-event. Roadmap
Stage 3 asserts every derived dichotomy carries exactly the missingness of its source.

## 4. The column contract  [Stage 1 acceptance criterion]

`COLUMN_CONTRACT: dict[str, Column]` keyed by the **verbatim** raw header, including its whitespace.
Three headers carry trailing spaces (`'PrestrokemRS '`, `'mRS56at90days   '`, `'Status '`); the
contract must not strip them, and the loader must not strip headers before matching. Stripping either
side turns a schema change into a silent rename.

`assert_column_contract(raw_columns: Iterable[str]) -> None` raises `SchemaError` on any symmetric
difference between the contract's keys and the workbook's columns. It materialises its argument to a
list on entry, so a generator can be passed and the error message can still name the offenders.

```
assert_column_contract(raw_columns)
        │
        ├─ materialise to list ────────────────────────────────────┐
        │                                                          │
        ├─ len(cols) == 41 and set(cols) == CONTRACT - {'Unnamed: 41'}
        │       └─→ SchemaError: "41 columns — this is the signature of a
        │              header-only read (nrows=0). pandas materialises the
        │              headerless final column only once rows are read.
        │              Re-read with rows."                    [REV]
        │
        ├─ missing = CONTRACT - cols        (in contract, absent from workbook)
        ├─ extra   = cols - CONTRACT        (in workbook, absent from contract)
        │
        ├─ both empty ─→ return None
        │
        └─→ SchemaError naming every offender with repr(), and for each
            `extra`, any contract key equal to it after .strip() on both
            sides — "'PrestrokemRS' is not in the contract; did you mean
            'PrestrokemRS ' (note trailing space)?"
```

**[REV] The 41-column branch.** The pre-review draft documented the header-only-read hazard in prose
and left the check to produce a generic "missing `Unnamed: 41`" message. Anyone who hits it will read
that as a schema change and go looking in the workbook. Naming the actual cause in the error turns
the subtlest trap in this stage into a self-explaining failure. `N_COLUMNS_EXPECTED = 42` is exported
for the same reason.

> **The check must run against a full read, never a header-only read.** `parse(sheet, nrows=0)`
> returns 41 columns; a full parse returns 42. The workbook's last column is headerless
> (`Unnamed: 41`, 13 dates), and pandas only materialises it once rows are read. A contract test
> written against the header row would pass while blind to a column carrying data.

### 4.1 Mapped — 25 columns

`dtype` is the third `Column` field. A blank cell means `None`: pandas infers, which is correct for
the continuous volumes and times. `"Int64"` is the nullable integer dtype, chosen wherever a column
is a count, a score or a 0/1 flag so that missingness survives as `<NA>` rather than forcing the
column to float. `"string"` is chosen where inference would produce mixed types.

| Raw header | Analysis name | dtype | Note |
|---|---|---|---|
| `CaseID` | `case_id` | `string` | formats differ by centre (`SSR-HUG-…`, `L-…`, bare integers); inference yields mixed int/str |
| `Center` | `center` | `string` | HUG is the integer `1`, the others are names; recoded via `CENTER_RECODE` |
| `Age` | `age` | `Int64` | §6 |
| `Sex` | `sex` | `Int64` | §6; coding undocumented, see §10 |
| `MedHistHypertension` | `hypertension` | `Int64` | §6 negative control |
| `MedHistHyperlipidemia` | `hyperlipidemia` | `Int64` | §6 negative control |
| `MedHistDiabetes` | `diabetes` | `Int64` | §6 negative control |
| `MedHistSmoking` | `smoking` | `Int64` | §6 negative control |
| `MedHistAtrialFibr` | `atrial_fib` | `Int64` | §6 |
| `PrestrokemRS ` | `prestroke_mrs` | `Int64` | §6 |
| `Wakeupstroke` | `wake_up` | `Int64` | source of `onset_type` [§5] |
| `Unwitnessedstroke` | `unwitnessed` | `Int64` | source of `onset_type` [§5] |
| `NIHSSonadmission` | `nihss_baseline` | `Int64` | §6 |
| `IVTwithrtPA` | `ivt` | `Int64` | treatment |
| `TimefromONSETtoIVTmin` | `onset_to_ivt_min` | | descriptive only; structurally non-applicable in controls [§11]; denylisted |
| `TimefromONSETtogroinmin` | `onset_to_groin_min` | | post-exposure [§12]; reported by arm; denylisted |
| `HypoperfusedtissuevolumeTmax6sml` | `tmax6_ml` | | §6 |
| `IschemiccorevolumeCBF30ml` | `core_ml` | | §6 |
| `Penumbravolumeml` | `penumbra_ml` | | excluded from all models [§6]; balance-only and subgroup use |
| `mRSscoreat90days` | `mrs_90d` | `Int64` | primary outcome, and the source of four dichotomies |
| `Symptomaticintracranialhaemorrhage` | `sich` | `Int64` | safety outcome |
| `Parenchymalhaematomatype2` | `ph2` | `Int64` | safety outcome |
| `TICI_2b_3` | `tici_2b_3` | `Int64` | secondary outcome; no ordinal source exists |
| `Contraindications_to_IVT` | `contraindication_reason` | `string` | **presence only** — the text is never read [DECISION 1] |
| `IVT_contraindicated_binary` | `ivt_contraindicated` | `Int64` | the [§4] eligibility classifier [DECISION 1] |

### 4.2 Dropped, each with its reason — 17 columns

Dropped columns carry `name=None` and no dtype. They are still read — `read_excel` returns the whole
sheet — so the only thing standing between a later module and `df['Deathat90days']` is §7.

| Raw header | Reason recorded in the contract |
|---|---|
| `FirstbrainimageMRI` | constant 0 across all 126 records; a centre-pathway attribute, not a patient characteristic [§6] |
| `HIR` | excluded a priori from every model and sensitivity analysis [§6] |
| `NIHSSscoreat24h` | post-time-zero, and not in the [§5] outcome registry |
| `NIHSSscore02at24 h` | as above |
| `8pointsNIHSSreductionat24 h` | as above |
| `NIHSS01at24h` | as above |
| `EarlyNeurologicalRecovery` | as above |
| `ChangeinNIHSSscorefrombaselineto24h` | as above |
| `mRSscore01at90days` | shipped dichotomy; rebuilt from `mrs_90d` [§5] |
| `mRS02at90days` | shipped dichotomy; rebuilt from `mrs_90d` [§5] |
| `Deathat90days` | shipped dichotomy; rebuilt as `mrs_90d == 6` [DECISION 2]. Disagrees with the ordinal source on 2 records |
| `mRS56at90days   ` | shipped dichotomy; rebuilt as `mrs_90d >= 5` [DECISION 2]. Disagrees on 3 records |
| `Deathat7days` | not in the [§5] outcome registry |
| `Status ` | free-text administrative field (`DCD`, `CG`, dated contact notes). **Vital-status-adjacent and deliberately unread**: the 90-day mRS is the sole ground truth for vital status [DECISION 2] |
| `Unnamed: 38` | empty |
| `Unnamed: 39` | empty |
| `Unnamed: 41` | headerless HUG-only contact date; post-time-zero administrative |

The six 24-hour NIHSS columns and `Deathat7days` are analysed in the pilot but are **not** in this
plan's [§5] registry. Dropping them is a contract decision, not an oversight; adding any of them back
requires amending [§5] first.

`FirstbrainimageMRI` is dropped by kind, not by variance. [§6] excludes it as a centre-pathway
attribute largely collinear with `center`, and adds "dropped automatically if constant" as a second
argument for the same conclusion. It is therefore not an exception to the assert-don't-hardcode
principle that governs `EXPECTED_NEVER_IVT` (§5): the exclusion does not depend on the observed
constancy, so there is nothing to assert.

### 4.3 Derived views of the contract  **[REV — new]**

Every consumer reads one of these, never `COLUMN_CONTRACT` directly. All four are computed from the
contract; none may be written out as a second literal.

```python
RENAME = {raw: c.name for raw, c in COLUMN_CONTRACT.items() if c.name is not None}
DROPPED = frozenset(raw for raw, c in COLUMN_CONTRACT.items() if c.name is None)
ANALYSIS_NAMES = frozenset(RENAME.values())
READ_DTYPES = {raw: c.dtype for raw, c in COLUMN_CONTRACT.items() if c.dtype is not None}
```

**Analysis names must be unique.** If two raw headers mapped to the same analysis name, `df.rename`
would silently collapse them and one column's data would vanish with no error anywhere. Nothing in
the pre-review draft checked this. §9.2 asserts `len(RENAME) == len(ANALYSIS_NAMES) == 25`, and
additionally that no analysis name is also a contract key — a name that is both would make the §7
scan ambiguous.

**Duplicate raw headers.** If a corrected workbook ships the same header twice, pandas mangles the
second to `Header.1`, which `assert_column_contract` reports as one added and zero missing columns.
That is the correct outcome — a loud `SchemaError` — but the message will name `'Header.1'` rather
than saying "duplicate". Worth knowing when reading the error; no special handling is specified.

### 4.4 The `Int64` reader is stricter than the pilot's  **[REV — new]**

The `Int64` dtypes in §4.1 are a deliberate trade. They preserve missingness as `<NA>` instead of
forcing scored columns to float, which is what makes roadmap Stage 3's "every derived dichotomy has
exactly the missingness of its source" assertion clean. The cost is that a non-numeric value landing
in a scored column raises at read time rather than being absorbed as `NaN`.

`NA_VALUES = ("N/A", "")` (§5.1) covers the two sentinels roadmap Stage 2 anticipates, and pandas
applies `na_values` before dtype conversion, so those convert cleanly. Anything else — a stray
comment in a score cell, a `?`, a date — raises. That is the intended behaviour under this stage's
Goal, but it is a stricter reader than `pilots/data.py`, and the first symptom of a hand-edited
workbook will be a pandas conversion error rather than a missing value. Do not respond to that by
widening `NA_VALUES`: add the specific sentinel, or fix the workbook.

## 5. Registries and constants

### 5.1 Paths and data identity

Every path is anchored to the module file, never to the working directory. This matches
`stage0_data_inventory.py:20`, which already does it correctly.

```python
ROOT = Path(__file__).resolve().parents[1]          # the repository root
DATA_XLSX = ROOT / "data" / (
    "Excel_bridging_EXTEND_paper_HUG_CHUV_LUGANO_USZ_def_v7_july26_"
    "with_abs_contra_indication.xlsx")
SHEET = "Feuil1"
OUT = ROOT / "out"
TABLES, FIGURES, LOGS = OUT / "tables", OUT / "figures", OUT / "logs"

FIXTURE_XLSX = Path(__file__).resolve().parent / "tests" / "fixture_schema.xlsx"

N_RECORDS_EXPECTED = 126
N_COLUMNS_EXPECTED = 42
NA_VALUES = ("N/A", "")     # roadmap Stage 2; load-bearing, not a guard — see below

DATA_SHA256 = "54934fbb2ae22647a9c0a2cbaff7ac425ed8948bcaeabe36d69aca00df657371"
```

**[REV] Paths were bare relative strings** (`data/…`, `../out/`) in the pre-review draft, which
resolve against the working directory. Anchoring is not a style preference here: Stage 14 is a single
entry point that a reader will plausibly invoke from the repository root.

**[REV 2026-08-10] `NA_VALUES` is load-bearing, and this line's earlier comment — "Stage 0 found no
sentinels, this is a guard" — was false.** `specs/stage2_load_and_clean.md` §4.2 read the workbook
with `dtype=object, keep_default_na=False` and found the literal three-character string `N/A` in
`mRSscoreat90days` (2 cells), `TICI_2b_3` (3) and each of the six 24-hour NIHSS columns (1 each). The
Stage 0 note reports them as blanks because `'N/A'` is already in pandas' default NA list, and its
"2 missing" for `mrs_90d` and "3 missing" for `tici_2b_3` *are* these strings. So this tuple is doing
real work on the primary outcome. Two consequences carried into Stage 2: `keep_default_na` stays
`True` and is passed explicitly, and §4.4's instruction stands — a conversion error is repaired by
naming the specific sentinel, never by widening this tuple.

**[REV] `DATA_SHA256`.** The Goal above says a change to the workbook must fail loudly in Stage 1.
The column contract catches *schema* changes only. It catches nothing if the data owner returns a
corrected file with identical headers and two amended mRS values — and §10 records two open queries,
so exactly that is expected. `data/` is gitignored, so git cannot answer "which workbook produced
this estimate" either.

Stage 2 hashes `DATA_XLSX` and raises `DataVersionError` on mismatch:

```
DataVersionError:
  expected 54934fbb…  (v7_july26_with_abs_contra_indication)
  actual   9a1c33be…
  The workbook changed. Re-run stage0_data_inventory.py, review
  ../out/stage0_data_inventory.md, and update DATA_SHA256 in the same commit.
```

The point is not to forbid a new workbook. It is to make adopting one a deliberate commit that also
re-runs Stage 0, rather than a file drop that silently re-bases every estimate. Stage 14's run
summary records the hash alongside the seed.

### 5.2 Inference and thresholds

```python
SEED = 20260807          # recorded in the run summary
N_BOOT = 2000
SMD_THRESHOLD = 0.10     # §9
RARE_MINORITY_THRESHOLD = 10
```

`RARE_MINORITY_THRESHOLD` is named for the minority cell, `min(events, non-events)`, not for events
[§8 amendment, DECISION 3]. A constant named `RARE_EVENT_THRESHOLD` would invite the exact misreading
the amendment corrects — TICI 2b-3 has 114 events and 6 non-events, so it passes an event-count rule
while failing that rule's rationale.

### 5.3 Treatment and centres

```python
TREATMENT = "ivt"                                   # 1 = IVT before EVT (bridging), 0 = EVT alone
TREATMENT_LABELS = {0: "EVT alone", 1: "bridging"}

CENTER_RECODE = {"1": "HUG", "Lausanne": "CHUV", "Lugano": "Lugano", "USZ": "USZ"}
CENTER_ORDER = ("HUG", "CHUV", "Lugano", "USZ")     # HUG first (reference, n=52); rest per pilot
EXPECTED_NEVER_IVT = ("USZ",)
```

Observed centre sizes are HUG 52, CHUV 21, Lugano 31, USZ 22 — so `CENTER_ORDER` is *not* a
size ordering. HUG leads because it is the largest and the declared reference level (§5.4); the
remaining three keep the pilot's order so tables stay comparable.

**[REV] `TREATMENT` and `TREATMENT_LABELS` are new.** Without them every stage writes the literal
`"ivt"` and the literal `"bridging"`, in tables, figures and logs. `"ivt"` is an analysis name, so
§7's scan would not catch a typo in it. The pilot has `TREATMENT`; the pre-review draft dropped it.

**[REV] `CENTER_ORDER` had no value** in the pre-review draft — it was named and left blank. It is
the display and factor order, with HUG first as the largest centre and the reference level; the
remaining order matches the pilot.

`EXPECTED_NEVER_IVT` is an assertion target, never an operative rule. Stage 5 derives the
zero-bridging set from the data and asserts it equals this; hardcoding the exclusion instead would
keep excluding USZ if a corrected workbook gave it treated patients.

### 5.4 Factor levels and reference categories  **[REV — entirely new]**

```python
CATEGORICAL = ("onset_type", "center")

FACTOR_LEVELS = {
    "center":     CENTER_ORDER,
    "onset_type": ("witnessed", "unwitnessed", "wake_up"),
}
REFERENCE_LEVELS = {
    "center":     "HUG",         # largest centre, n = 52
    "onset_type": "witnessed",   # the clinical baseline
}
```

The pre-review draft named `CATEGORICAL` but declared neither the levels nor which level is the
reference. Roadmap Stage 6 then asks for "reference-coded dummies for factors, constant columns
dropped." Those two together are underdetermined, with two consequences:

1. `pd.get_dummies(drop_first=True)` drops the alphabetically first level, so the baseline centre
   would silently have been CHUV, and every coefficient in the output table would be relative to a
   baseline nobody chose.
2. Worse for [§10] inference. A bootstrap replicate that happens to contain no Lugano patients
   produces a design matrix with a *different width and different column meanings* — the Lugano
   dummy is absent rather than present-and-constant. Stage 10 stacks 2000 coefficient vectors; a
   subset would be misaligned.

Declaring the levels fixes both. Stage 6 builds each factor as a categorical with the declared level
set, so an absent level yields an all-zero column that the "constant columns dropped" rule then
removes — deterministically, in every replicate, with the drop recorded. The marginal estimates are
invariant to coding either way; the model-parameter outputs and the [§14a] conditional odds ratio are
not. `pilots/config.py:70` already carried `ONSET_ORDER`; the pre-review draft dropped it.

### 5.5 Covariates [§6]

```python
PS_COVARIATES: Final[tuple[str, ...]] = (
    "age", "sex", "prestroke_mrs", "nihss_baseline", "onset_type",
    "core_ml", "tmax6_ml", "atrial_fib", "center")

BALANCE_ONLY: Final[tuple[str, ...]] = (
    "hypertension", "hyperlipidemia", "diabetes", "smoking", "penumbra_ml")

OUTCOME_COVARIATES = PS_COVARIATES                                          # the same object
STANDARDISATION_COVARIATES = tuple(c for c in PS_COVARIATES if c != "center")   # §14a, derived
PS_COVARIATES_FULL = PS_COVARIATES + (
    "hypertension", "hyperlipidemia", "diabetes", "smoking")                # §13, derived

DERIVED_NAMES: Final[tuple[str, ...]] = ("onset_type",)   # produced by Stage 3, not by the contract
```

`OUTCOME_COVARIATES` binds the same object rather than repeating it, and the other two are computed
from it. None may be written out as a second literal: invariant 3 exists because two hand-maintained
lists drift.

**[REV] Tuples, not lists.** The aliasing above is right, and it is exactly what makes mutability
dangerous. A later module doing `covs = config.PS_COVARIATES; covs.append(...)` would rewrite the
specification for `PS_COVARIATES`, `OUTCOME_COVARIATES` and every module importing afterwards, for
the remainder of the process — including all 2000 refits in Stage 10. Every §9 test runs at import
time, before any of that could happen, so none would catch it. Tuples make the mutation an
`AttributeError` instead. Consumers write `df[list(PS_COVARIATES)]` where pandas wants a list; that
is the entire cost. `POST_TIME_ZERO` is a `frozenset` for the same reason.

**[REV] `DERIVED_NAMES` is new** so that §9.3 does not hardcode `"onset_type"` in a test. When Stage 3
gains a second derived covariate, the constant is the one place that changes.

### 5.6 Outcome-model overrides [§8 amendment, DECISION 3]

```python
OUTCOME_MODEL_OVERRIDES = {"tici_2b_3": ("center", "atrial_fib")}

def outcome_model_covariates(outcome: str) -> tuple[str, ...]:
    """m_a(X) covariates for `outcome`: the shared §6 set unless it declares an override.

    Treatment is NOT in the returned tuple. Every estimator adds the treatment main
    effect itself, so TICI's reduced model is intercept + treatment + center (2 df)
    + atrial_fib = 5 parameters against 6 non-events, as [§8 amendment] states.
    """
    if outcome not in OUTCOMES:
        raise KeyError(
            f"{outcome!r} is not in OUTCOMES. Known outcomes: {sorted(OUTCOMES)}. "
            "An outcome model may only be requested for a registered outcome.")
    return OUTCOME_MODEL_OVERRIDES.get(outcome, OUTCOME_COVARIATES)
```

Every consumer calls the function; nobody reads `OUTCOME_COVARIATES` directly for an outcome model.
A registry that defaults to the shared entry means an outcome can only diverge by being named. The
returned tuple must reach the reporting layer and print beside the TICI estimate.

**[REV] The unknown-key branch is new.** Without it a mistyped outcome name silently receives the
full §6 covariate set — a silent default of exactly the kind [§4] and [§8] forbid elsewhere.
**[REV] The treatment convention is new.** The pre-review draft never said whether treatment belonged
in the returned list, and the [§8 amendment]'s five-parameter count only reconciles if it does not.

### 5.7 Outcome registry [§5]

`OUTCOMES: dict[str, Outcome]`, exactly eight entries. **[REV]** The pre-review draft tabulated five
of the six `Outcome` fields, leaving `label` and `higher_is_better` to be invented; both are pinned
here. `higher_is_better` drives figure orientation only and never changes an estimate's sign.

| key | label | kind | family | source | op | threshold | higher_is_better |
|---|---|---|---|---|---|---|---|
| `mrs_90d` | mRS at 90 days | ordinal | primary | — | — | — | `False` |
| `mrs_0_2_90d` | mRS 0-2 at 90 days | binary | secondary | `mrs_90d` | `<=` | 2 | `True` |
| `mrs_0_1_90d` | mRS 0-1 at 90 days | binary | secondary | `mrs_90d` | `<=` | 1 | `True` |
| `tici_2b_3` | TICI 2b-3 | binary | secondary | — | — | — | `True` |
| `sich` | Symptomatic intracranial haemorrhage | binary | safety | — | — | — | `False` |
| `ph2` | Parenchymal haematoma type 2 | binary | safety | — | — | — | `False` |
| `death_90d` | Death at 90 days | binary | safety | `mrs_90d` | `==` | 6 | `False` |
| `mrs_5_6_90d` | mRS 5-6 at 90 days | binary | safety | `mrs_90d` | `>=` | 5 | `False` |

Exactly one entry has `family == "primary"`. The pilot had two, and [§8] rests on there being one
primary quantity with one test; §9.3 asserts the count.

### 5.8 Eligibility [§3, §4, DECISION 1]

```python
ELIGIBILITY_ORDER = ("eligible", "indeterminate", "ineligible")
```

The configuration holds **no list of contraindication reason strings**. Under DECISION 1 the
classifier is `ivt_contraindicated` plus whether `contraindication_reason` is present, so no string is
ever matched and there is nothing to enumerate. What Stage 4 asserts instead — recorded here as a
comment on the constant, since the [§4] "assert every reason is classified" requirement no longer has
anything to range over — is that `ivt_contraindicated` is 0/1 and never missing, and that **no treated
patient carries a 1**.

### 5.9 Post-time-zero denylist [invariant 4]

```python
POST_TIME_ZERO: Final[frozenset[str]] = frozenset(
    {"onset_to_ivt_min", "onset_to_groin_min", *OUTCOMES})
```

Built from the outcome registry so a new outcome is denylisted by existing, not by being remembered.
Defined after `OUTCOMES` in the module; the ordering is load-bearing.

### 5.10 Structural non-applicability and plausible ranges

```python
STRUCTURALLY_NON_APPLICABLE = {
    "onset_to_ivt_min": "control arm — no IVT was given, so the time does not exist [§11]",
}

PLAUSIBLE_RANGES = {                 # (low, high), inclusive; for Stage 2's schema assertions
    "age":             (18, 110),
    "nihss_baseline":  (0, 42),
    "prestroke_mrs":   (0, 5),       # 6 is death; impossible pre-stroke
    "mrs_90d":         (0, 6),
    "core_ml":         (0, None),
    "tmax6_ml":        (0, None),
    "penumbra_ml":     (0, None),
    "onset_to_ivt_min":   (0, None),
    "onset_to_groin_min": (0, None),
}

BINARY_COLUMNS: Final[tuple[str, ...]] = (
    "sex", "hypertension", "hyperlipidemia", "diabetes", "smoking", "atrial_fib",
    "wake_up", "unwitnessed", "ivt", "sich", "ph2", "tici_2b_3", "ivt_contraindicated")
```

**[REV] All three are new as structured constants.** The pre-review draft carried the ranges as prose
("volumes and times non-negative"), which Stage 2 cannot iterate over, and left structural
non-applicability entirely to Stage 3 to remember. `None` as an upper bound means unbounded above.

**[REV] `BINARY_COLUMNS` closes a hole the ranges left open.** Thirteen mapped columns are 0/1 flags
and none of them appeared in any assertion. Stage 2 asserts every value in these columns is in
`{0, 1}` or missing. Two of them matter more than the rest: `ivt` is the exposure, and a value of `2`
or `-1` there would flow into every weight, every estimate and every arm label without a single
error; `ivt_contraindicated` is the [§4] classifier, and [DECISION 1] already requires Stage 4 to
assert it is 0/1 and never missing — declaring the set here means Stage 4 asserts against a named
constant rather than a hardcoded pair of column names.

### 5.11 Labels pending a data query

```python
SEX_LABELS: dict[int, str] | None = None    # see §10; must stay None until the query is answered
```

**[REV] New.** §10 records that the `sex` coding is undocumented. Without a named constant, Stage 14
will guess when it builds the baseline table. Stage 14 must print the levels as `sex = 0` / `sex = 1`
while this is `None`, and must never invent labels.

## 6. Subgroups and target mismatch [§13]

```python
TARGET_MISMATCH = {"max_core_ml": 70, "min_ratio": 1.8, "min_penumbra_ml": 15}
MRS_THRESHOLDS = (0, 1, 2, 3, 4, 5)                   # cumulative RD_k [§8]
MISMATCH_RATIO_UNDEFINED_COUNTS_AS_MET = True
```

Core volume is exactly 0 in 50 of 125 non-missing records, so the mismatch *ratio* is undefined for
40% of the cohort. `MISMATCH_RATIO_UNDEFINED_COUNTS_AS_MET = True` means core = 0 with
penumbra > `min_penumbra_ml` satisfies the ratio criterion, so the choice appears in the
configuration rather than as an inline comparison in Stage 3.

**[REV] Missing core is not the same case as zero core.** One record has no core volume at all, and
one has no Tmax>6 s volume. The constant above covers `core == 0`; it does not cover `core is NA`.
Stage 3 must return **missing**, not "met" and not "not met", when any input to the criterion is
missing. A subgroup flag of `<NA>` drops that patient from the subgroup analysis by complete-case
[§11], which is correct; treating `<NA>` as 0 would silently sweep it into the mismatch-present group.

**[REV] Stage 3 must branch on the constant, not on a literal.** The whole reason
`MISMATCH_RATIO_UNDEFINED_COUNTS_AS_MET` exists is so the convention is visible in the configuration
rather than buried as an inline comparison. An implementation that reads
`if core == 0 and penumbra > 15:` has reproduced the convention while defeating its purpose, and the
constant becomes decorative. Stage 3 reads the constant and `TARGET_MISMATCH["min_penumbra_ml"]`.

The third subgroup, core volume above/below median [§13], has no constant here on purpose: the median
is a property of the cohort, computed by Stage 3 after the [§3] restrictions, and freezing it as a
number in the configuration would silently decouple it from the cohort it describes.

## 7. The no-raw-names rule

"Nothing outside this module may reference a raw column name" is testable, and is tested (§9.4)
rather than left as a convention.

**Where the rule comes from, and what it actually buys.** It began as a docstring aspiration in
`pilots/config.py`, became a build requirement in roadmap Stage 1, and is a tested criterion here.
Two things justify it, and the second is the load-bearing one:

- *Blast radius.* The raw headers are hostile: three carry trailing whitespace, two carry interior
  whitespace, one is invented by pandas. Scattered across ten modules, a v8 workbook breaks ten
  files in ten ways. Confined to the contract, it breaks once, in `assert_column_contract`.
- *It is the mechanical enforcement of DECISION 2 and DECISION 1.* `read_excel` returns the whole
  sheet, so `Deathat90days`, `mRS56at90days   `, `Status ` and the free text of
  `Contraindications_to_IVT` are all sitting in the DataFrame with no analysis name. Nothing but
  this rule stops a later module from reading the shipped death flag — which would move the primary
  safety outcome from 25 events to 26 — or from reading the free text that DECISION 1 rules out.

**[REV] The check is an AST string-literal scan, not a substring scan.** A raw header can only enter
code as a string literal, so that is the only thing worth checking. The pre-review draft scanned raw
file text for any contract key as a substring, which fails twice over: `Age` is a substring of
`average`, `Agent` and any comment reading "Age at baseline"; `HIR` is inside `THIRD`. And
`test_config.py` was not exempt, while §9.1 requires it to construct a trailing-whitespace variant —
so the check would have flagged its own test file on the first run. The predictable repair is to add
files to the exemption list until it passes, which disables the rule precisely for the modules that
handle raw headers. A rule exempted into uselessness is worse than no rule, because it reads as
enforcement.

```python
EXEMPT_FROM_RAW_NAME_SCAN = {
    "config.py":                "owns the contract",
    "stage0_data_inventory.py": "predates the configuration; reads raw headers by design",
    "test_config.py":           "must construct malformed headers to test the contract check",
}
```

The scan walks `extended_bridging/**/*.py` recursively, skipping `.venv/` and `__pycache__/`, parses
each file with `ast`, collects every `ast.Constant` whose value is a `str`, and fails on exact
equality with a contract key — naming the file, the line number and the header. Comments and
identifiers are not in the AST and are therefore invisible to it.

**The comparison is case-sensitive, and must stay that way.** Two contract keys differ from their own
analysis names only by case: `TICI_2b_3` → `tici_2b_3`, and `Center` → `center`. A case-insensitive
scan would flag every legitimate use of `tici_2b_3` and `center` in every downstream module — which
is most of them — and the repair would be to gut the scan. Exact, case-sensitive equality is the
specification, not an implementation detail.

The scan sees f-strings as `ast.JoinedStr` with `ast.Constant` parts, so `f"{df['Age']}"` would still
be caught; a header assembled at runtime from concatenated fragments would not. That residual hole is
accepted: it requires deliberate effort to construct, and no legitimate reason to.

## 8. Data flow into Stage 2

Stage 1 declares; Stage 2 is its first consumer. The handoff, so that the contract check and the
version check land in the right order:

```
  Stage 2 reader
      │
      ├─ hash DATA_XLSX ──── != DATA_SHA256 ──→ DataVersionError    (§5.1)
      │                                          "re-run Stage 0"
      ├─ read_excel(DATA_XLSX, sheet_name=SHEET,
      │             dtype=READ_DTYPES,          ← from the contract, §4.1
      │             na_values=list(NA_VALUES))  ← §5.1
      │             # full read. never nrows=0. see §4
      │
      ├─ assert_column_contract(df.columns) ──→ SchemaError         (§4)
      │
      ├─ len(df) != N_RECORDS_EXPECTED ──────→ SchemaError
      │
      ├─ rename via RENAME (contract, name is not None)
      ├─ drop columns where name is None
      ├─ recode center via CENTER_RECODE, assert no unmapped code
      ├─ content-driven corrections → audit log with case identifiers
      └─ range assertions from PLAUSIBLE_RANGES                     (§5.10)
```

The hash check precedes the schema check deliberately: if the workbook was replaced, "which file is
this" is the more useful error than "column X is missing."

## 9. Acceptance criteria

All tests live in `test_config.py`. Tests that need the private workbook are marked
`@pytest.mark.skipif(not DATA_XLSX.exists())`; every other test runs on any checkout, against
`FIXTURE_XLSX`.

1. **Contract exhaustiveness.** Against a **full** read of `FIXTURE_XLSX`, `assert_column_contract`
   accepts the 42 columns, and raises `SchemaError` on: an added column, a removed column, and a
   header differing only in trailing whitespace — with the strip-match hint present in the message.
   **[REV]** Two further cases: a header-only read (41 columns) raises with the header-only-read
   cause named, not a generic missing-column message; and an empty iterable raises. **[REV]** One
   data-gated test asserts the fixture's headers equal the real workbook's, so the fixture cannot
   drift.
2. **No silent passthrough.** **[REV]** Every `Column` has a non-empty `reason` — *every* column, not
   only dropped ones. The pre-review draft's "either a non-`None` name **or** a non-empty reason" is
   an OR, satisfied by a mapped column with a blank reason, which contradicts §3.2's "Never empty,
   ever." Concretely: `len(reason) >= 20` for every entry, and `reason.strip().lower()` is not in
   `{"", "-", "n/a", "tbd", "todo", "see above", "dropped"}`. The threshold is a number, not "more
   than a few characters", so the test is deterministic. **[REV]** Additionally the §4.3 derived
   views: `len(COLUMN_CONTRACT) == N_COLUMNS_EXPECTED == 42`; `len(RENAME) == 25` and
   `len(DROPPED) == 17`; `len(RENAME) == len(ANALYSIS_NAMES)`, so no two raw headers collapse onto
   one analysis name; and `ANALYSIS_NAMES & set(COLUMN_CONTRACT) == ∅`, so no analysis name is also
   a raw header.
3. **Registry integrity.** Every `Outcome.source` is either `None` or a mapped analysis name; every
   `PS_COVARIATES` entry is a mapped analysis name or a member of `DERIVED_NAMES`; the three families
   partition the registry and each is non-empty; `OUTCOMES` has exactly the eight [§5] entries.
   **[REV]** Additionally: exactly one entry has `family == "primary"`; `op` and `threshold` are both
   `None` iff `source` is `None`; every `op` appearing in the registry has an `OPS` entry; and
   `Outcome.rule` round-trips (`"mrs_90d <= 2"` for `mrs_0_2_90d`, `None` for `mrs_90d`).
4. **No raw names outside `config.py`.** The AST scan of §7 over `extended_bridging/**/*.py`, with
   `EXEMPT_FROM_RAW_NAME_SCAN` as the only exemptions, fails naming file, line and header.
   **[REV]** A companion test asserts the scan actually fires: it parses a synthetic snippet
   containing `'PrestrokemRS '` and confirms a failure, so a scan that silently matches nothing
   cannot pass as a green test.
5. **Invariant 3.** `outcome_model_covariates(o) is PS_COVARIATES` for every registered outcome
   except `tici_2b_3`; `outcome_model_covariates("tici_2b_3") == ("center", "atrial_fib")`; and
   `set(OUTCOME_MODEL_OVERRIDES) == {"tici_2b_3"}`, asserted against that literal so a second
   override fails the test rather than passing quietly. **[REV]** `outcome_model_covariates` raises
   `KeyError` on an unregistered name, and `TREATMENT` is in none of the returned tuples.
6. **Invariant 4.** `set(PS_COVARIATES) & POST_TIME_ZERO == set()`, and likewise for
   `PS_COVARIATES_FULL` and `STANDARDISATION_COVARIATES`. **[REV]** `POST_TIME_ZERO` contains every
   key of `OUTCOMES`, so adding an outcome cannot bypass the denylist.
7. **Derived-not-duplicated.** `OUTCOME_COVARIATES is PS_COVARIATES`; `"center" not in
   STANDARDISATION_COVARIATES` while every other §6 covariate is present; `PS_COVARIATES_FULL`
   starts with `PS_COVARIATES` and adds exactly the four vascular risk factors.
8. **[REV] Immutability.** `PS_COVARIATES`, `BALANCE_ONLY`, `CENTER_ORDER`, `ELIGIBILITY_ORDER`,
   `MRS_THRESHOLDS` and every `FACTOR_LEVELS` value are tuples; `POST_TIME_ZERO` is a frozenset;
   `Column` and `Outcome` reject attribute assignment.
9. **[REV] Read contract.** `READ_DTYPES` keys are a subset of `COLUMN_CONTRACT` keys and contain no
   dropped column; every dtype string is one of `{"string", "Int64", "Float64"}`; reading
   `FIXTURE_XLSX` with `dtype=READ_DTYPES` succeeds and yields a `string` dtype for `CaseID`.
10. **[REV] Factor declarations.** `set(FACTOR_LEVELS) == set(CATEGORICAL) == set(REFERENCE_LEVELS)`;
    every reference level is a member of its own level tuple; `CENTER_ORDER` equals
    `FACTOR_LEVELS["center"]` and its members are exactly the values of `CENTER_RECODE`;
    `EXPECTED_NEVER_IVT` is a subset of `CENTER_ORDER`.
11. **[REV] Constants sanity.** `DATA_SHA256` is 64 lowercase hex characters; `N_BOOT >= 1000`;
    `0 < SMD_THRESHOLD < 1`; `RARE_MINORITY_THRESHOLD >= 1`; `MRS_THRESHOLDS` equals
    `tuple(range(PLAUSIBLE_RANGES["mrs_90d"][1]))` so the cumulative thresholds cannot fall out of
    step with the mRS range; `SEED` is a positive int. Data-gated: the file at `DATA_XLSX` hashes to
    `DATA_SHA256`.
12. **[REV] Exclusions hold.** `penumbra_ml` and `hir` appear in no covariate tuple — not
    `PS_COVARIATES`, not `PS_COVARIATES_FULL`, not `STANDARDISATION_COVARIATES`, not any value of
    `OUTCOME_MODEL_OVERRIDES`. `BALANCE_ONLY ∩ PS_COVARIATES == ∅`. `hir` has no analysis name at
    all, so this is belt-and-braces against a future edit that maps it.
13. **[REV] Names resolve.** Every entry of `PS_COVARIATES_FULL`, `BALANCE_ONLY`, `BINARY_COLUMNS`,
    `PLAUSIBLE_RANGES`, `STRUCTURALLY_NON_APPLICABLE` and `POST_TIME_ZERO` is either in
    `ANALYSIS_NAMES`, in `DERIVED_NAMES`, or a key of `OUTCOMES`. A constant naming a variable that
    no stage can produce is a typo that would otherwise surface as a `KeyError` deep in Stage 6 or
    Stage 9, hours into a bootstrap.
14. **[REV] Structural claims from Stage 0 hold.** Data-gated, against the real workbook: 42 columns
    on a full read and 41 on a header-only read; 126 rows; `CaseID` unique; every `Center` value is a
    key of `CENTER_RECODE`; `FirstbrainimageMRI` has one distinct value. These are the facts §4 and
    §5 were written against, and a corrected workbook that satisfies `DATA_SHA256` cannot violate
    them — but the test documents the dependency, so updating the hash also re-checks the premises.

### Coverage map

```
config.py                                              test_config.py
├── module import (constants)                          ├── 9.2  reasons non-empty, >= 20 chars
│   ├── COLUMN_CONTRACT (42)                           ├── 9.1  exhaustive vs fixture + workbook
│   ├── RENAME / DROPPED / ANALYSIS_NAMES  (derived)   ├── 9.2  25/17 split, names unique
│   ├── READ_DTYPES        (derived)                   ├── 9.9  keys ⊆ contract, dtypes valid
│   ├── STANDARDISATION_COVARIATES (derived)           ├── 9.7  center excluded, rest present
│   ├── PS_COVARIATES_FULL (derived)                   ├── 9.7  prefix + exactly 4 added
│   ├── POST_TIME_ZERO     (derived, frozenset)        ├── 9.6  disjoint x3, covers OUTCOMES
│   ├── FACTOR_LEVELS / REFERENCE_LEVELS               ├── 9.10 refs ∈ levels, keys agree
│   ├── BINARY_COLUMNS / PLAUSIBLE_RANGES              ├── 9.13 every name resolves
│   ├── DATA_SHA256, SEED, N_BOOT, thresholds          ├── 9.11 format + range sanity
│   └── tuple/frozenset types                          └── 9.8  immutability
│
├── assert_column_contract(cols)
│   ├── exact 42 match          → None                 ├── 9.1
│   ├── added column            → SchemaError          ├── 9.1
│   ├── removed column          → SchemaError          ├── 9.1
│   ├── trailing-ws variant     → SchemaError + hint   ├── 9.1
│   ├── 41 cols (header-only)   → SchemaError + cause  ├── 9.1  [REV]
│   └── empty iterable          → SchemaError          └── 9.1  [REV]
│
├── outcome_model_covariates(o)
│   ├── "tici_2b_3"             → override tuple       ├── 9.5
│   ├── any other registered    → PS_COVARIATES (is)   ├── 9.5
│   └── unregistered            → KeyError             └── 9.5  [REV]
│
└── Outcome.rule
    ├── source is None          → None                 ├── 9.3  [REV]
    ├── source set              → "mrs_90d <= 2"       ├── 9.3  [REV]
    └── every op has an OPS entry                      └── 9.3  [REV]

Data-gated (skipif not DATA_XLSX.exists()):            ├── 9.1  fixture headers == workbook
                                                       ├── 9.11 file hashes to DATA_SHA256
                                                       └── 9.14 42/41 cols, 126 rows, CaseID
                                                                unique, centres map, MRI constant

Every branch above has a test. No branch is untested.
```

### Failure modes

| Failure | Test | Error handling | Visible or silent |
|---|---|---|---|
| Workbook replaced, same schema | 9.11 (data-gated) | `DataVersionError` names both hashes and the fix | visible |
| Column added or renamed | 9.1 | `SchemaError` names it with `repr()` | visible |
| Header whitespace stripped upstream | 9.1 | `SchemaError` + strip-match hint | visible |
| Header-only read used by mistake | 9.1 | `SchemaError` names the cause | visible |
| A later stage reads `Deathat90days` | 9.4 | AST scan fails, naming file and line | visible |
| A later stage mutates `PS_COVARIATES` | 9.8 | `AttributeError` at the mutation site | visible |
| Unregistered outcome requested | 9.5 | `KeyError` listing known outcomes | visible |
| Derived dichotomy loses missingness | roadmap Stage 3 | `.mask(src.isna())` in the §3.2 loop | visible via Stage 3 test |
| Subgroup criterion with missing core | roadmap Stage 3 (§6) | returns `<NA>`, complete-case drops it | visible in denominators |

No failure mode in this stage is both untested and silent.

## 10. Known gaps carried into later stages

- **`sex` coding is undocumented** (0/1, 74 vs 52). Adjustment is unaffected, but the baseline table
  cannot label the levels. A data query, not a blocker. `SEX_LABELS` (§5.11) stays `None` until it is
  answered, and Stage 14 prints `sex = 0` / `sex = 1` rather than guessing.
- **Two records where the shipped death flag contradicts the mRS** are read as alive [DECISION 2].
  Standing query with the data owner; identifiers are in the Stage 0 note. If a corrected workbook
  resolves them, `DATA_SHA256` forces the update to be a deliberate commit that re-runs Stage 0.
- **`onset_to_ivt_min` is structurally non-applicable in controls**, not missing [§11]. Recorded in
  `STRUCTURALLY_NON_APPLICABLE` (§5.10) and denylisted; Stage 3 marks the distinction.
- **One record has no core volume and one has no Tmax>6 s volume.** Stage 3's mismatch criterion must
  return `<NA>` for these rather than applying the core = 0 convention (§6).

## 11. What Stage 1 deliberately does not decide

Recorded so a later stage does not look here for an answer that was never placed here.

- **Which centres are dropped.** `EXPECTED_NEVER_IVT` is an assertion target. Stage 5 derives the set
  from the data and asserts equality (§5.3).
- **How eligibility is classified.** `ELIGIBILITY_ORDER` names the categories; Stage 4 owns the rule.
- **How `onset_type` is built.** `FACTOR_LEVELS` names the three levels; Stage 3 builds it from
  `wake_up` and `unwitnessed` and asserts they are never both positive.
- **Whether an outcome is augmented.** `RARE_MINORITY_THRESHOLD` is the number; Stage 9 applies it to
  `min(events, non-events)`.
- **How the design matrix is coded.** `FACTOR_LEVELS` and `REFERENCE_LEVELS` are the declarations;
  Stage 6 builds the matrix and drops constant columns.

## 12. NOT in scope for Stage 1

| Considered | Why deferred |
|---|---|
| Reading, cleaning, correcting data | Stage 2. Stage 1 declares only; the contract check is a pure function of a column list. |
| Deriving `onset_type` and the dichotomies | Stage 3, using `OPS` and the registry from §3.2. |
| Eligibility classification | Stage 4. Stage 1 holds only `ELIGIBILITY_ORDER`. |
| CI, packaging, distribution | This is internal analysis code with no published artifact and one runtime environment. Nothing to build or ship. |
| Full-covariate and no-centre sensitivity propensity models [§13] | Roadmap defers the [§13] sensitivity suite. `PS_COVARIATES_FULL` is declared now because it is a pure derivation; nothing consumes it yet. |
| A `dtype` for every mapped column | Only the columns where inference is wrong or dangerous are typed (§4.1). Over-typing turns a benign float into a read-time crash on a corrected workbook. |
| Validating `data/` exists at import | `config.py` must import without the private workbook, or the pure-logic tests cannot run anywhere. Existence is Stage 2's problem. |

## 13. What already exists, and what to lift

`pilots/config.py` is gitignored but present, and already carries correct versions of much of this.
Lift these rather than re-deriving them, checking each against §4 and §5 above:

| From the pilot | Status |
|---|---|
| `RENAME` (raw → analysis names) | lift the 25 mapped entries; the pilot maps several columns this contract drops |
| `CENTER_RECODE` | lift verbatim |
| `CENTER_ORDER`, `ONSET_ORDER` | lift; `ONSET_ORDER` becomes `FACTOR_LEVELS["onset_type"]` |
| `TARGET_MISMATCH` | lift verbatim |
| `MRS_THRESHOLDS` | lift, as a tuple |
| `PS_COVARIATES`, `BALANCE_ONLY`, `PS_COVARIATES_FULL` | lift, as tuples |
| `TREATMENT` | lift |
| `OUTCOMES` | **do not lift.** The pilot has 14 entries, two families marked primary, and outcomes this plan's [§5] registry excludes. Build from §5.7. |
| `OUTCOME_COVARIATES` | **do not lift.** The pilot uses a reduced set; [§8] requires the same object as `PS_COVARIATES`. |
| `IVT_INELIGIBLE_REASONS` / `IVT_ELIGIBLE_REASONS` | **do not lift.** DECISION 1 retires string matching entirely (§5.8). |
| `IPTW_TRIM`, `PS_COVARIATES_REDUCED`, `PS_COVARIATES_V1`, `NOT_REPRODUCIBLE` | **do not lift.** No [§13] or [§15] analysis in this plan uses them. |

`stage0_data_inventory.py:20` already has the correct path-anchoring idiom (`ROOT =
Path(__file__).resolve().parents[1]`); §5.1 reuses it exactly.

## 14. Implementation tasks

Ordered. Each is independently verifiable; there is one lane, so no parallelisation applies.

- [ ] **T1 (P1)** — `pyproject.toml` + `uv.lock`: dependencies per §2 including `numpy<2.1`, no
      `[build-system]`, `testpaths = ["."]`. Verify: `uv sync` succeeds.
- [ ] **T2 (P1)** — `config.py` §3: `SchemaError`, `DataVersionError`, `Column`, `Outcome`, `OPS`.
      Stdlib imports only. Verify: `python -c "import config"` with no scientific stack installed.
- [ ] **T3 (P1)** — `COLUMN_CONTRACT`, all 42 entries with dtypes and reasons per §4.1 and §4.2;
      `RENAME`, `DROPPED`, `ANALYSIS_NAMES`, `READ_DTYPES` derived per §4.3. Verify: acceptance
      9.2, 9.9.
- [ ] **T4 (P1)** — `assert_column_contract` per §4, including the 41-column branch and the
      strip-match hint. Verify: acceptance 9.1.
- [ ] **T5 (P1)** — `tests/fixture_schema.xlsx`: 42 verbatim headers, two invented rows, headerless
      final column carrying a value. Verify: acceptance 9.1 runs green without `data/` present.
- [ ] **T6 (P1)** — constants §5.1 through §5.11, in dependency order (`OUTCOMES` before
      `POST_TIME_ZERO`). Verify: acceptance 9.3, 9.6, 9.7, 9.10, 9.11, 9.12, 9.13.
- [ ] **T7 (P1)** — `outcome_model_covariates` per §5.6, with the unknown-key branch. Verify:
      acceptance 9.5.
- [ ] **T8 (P1)** — the AST scan and `EXEMPT_FROM_RAW_NAME_SCAN` per §7, case-sensitive, plus the
      test that proves the scan fires. Verify: acceptance 9.4.
- [ ] **T9 (P2)** — `test_config.py` remainder: acceptance 9.8, and the data-gated tests behind
      `skipif` (9.1 fixture-vs-workbook, 9.11 hash, 9.14 structural facts). Verify: `uv run pytest -v`
      green both with and without `data/` present.

### Definition of done

Stage 1 is complete when all of the following hold, and not before:

1. `uv run pytest -v` is green **with** `data/` present — every test, including the data-gated ones.
2. `uv run pytest -v` is green **with `data/` temporarily renamed** — the data-gated tests skip, and
   nothing else fails or errors at collection.
3. `git check-ignore -v extended_bridging/tests/fixture_schema.xlsx` reports nothing, and the fixture
   is committed alongside `uv.lock`.
4. Acceptance 9.4 fails when a raw header is deliberately pasted into a scratch `.py` file under
   `extended_bridging/`, and passes again when it is removed. A green scan that has never been seen
   to fail is not evidence.
5. No file other than `config.py`, `stage0_data_inventory.py` and `test_config.py` exists in
   `extended_bridging/` yet — Stages 2 onward are not started, so criterion 9.4 has nothing to find
   and criterion 4 above is the only proof it works.

## 15. Verification record

The facts in §4 and §5 were checked against the workbook on 2026-08-07, before this spec was
finalised. Recorded so that a reader can tell which numbers were verified rather than carried over
from the Stage 0 note or the pilot, and so the checks are re-runnable after a workbook update.

| Claim | Where used | Verified |
|---|---|---|
| Full read is 126 × 42; header-only read is 41 columns | §4, §5.1 | yes |
| `Unnamed: 41` absent from a header-only read, 13 non-null values on a full read | §4, §4.2 | yes |
| Trailing whitespace on exactly `'PrestrokemRS '`, `'mRS56at90days   '`, `'Status '` | §4 | yes |
| 25 mapped + 17 dropped covers all 42 with no overlap and no omission | §4.1, §4.2 | yes |
| 25 distinct analysis names | §4.3 | yes |
| `CaseID` has no duplicates and is object-typed (mixed formats) | §4.1 | yes |
| `Center` takes exactly the four keys of `CENTER_RECODE` | §5.3 | yes |
| `FirstbrainimageMRI` has one distinct value, `0` | §4.2 | yes |
| `sex` is 0/1 with counts 74 / 52 | §5.11, §10 | yes |
| `core_ml` is exactly 0 in 50 of 125 non-missing; 1 record missing | §6, §10 | yes |
| `tmax6_ml` has 1 missing; `mrs_90d` has 2 missing; `tici_2b_3` has 3 missing | §4.1, §10 | yes |
| Derived death (`mrs_90d == 6`) disagrees with `Deathat90days` on 2 records | §4.2 | yes |
| Derived `mrs_90d >= 5` disagrees with `mRS56at90days   ` on 3 records | §4.2 | yes |
| Derived `<= 2` and `<= 1` agree with their shipped columns on every record | §4.2 | yes |
| `Wakeupstroke` and `Unwitnessedstroke` are never both 1 | §5.4 | yes |
| `IVT_contraindicated_binary` is 0/1 with no missing (107 / 19) | §5.8 | yes |
| Observed ranges: age 28-97, NIHSS 2-32, prestroke mRS 0-4 — all inside `PLAUSIBLE_RANGES` | §5.10 | yes |
| SHA-256 of the workbook is `54934fbb…f657371` | §5.1 | yes |

Acceptance 9.14 re-checks the structural subset of this table in code, so the record does not go
stale silently.
