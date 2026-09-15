# Stage 14 spec — outputs and guardrails

Implements roadmap Stage 14 [§16]. Section references in brackets are to
`statistical_analysis_plan.md`. Decisions referenced as DECISION *n* are recorded in
`../out/stage0_data_inventory.md`. Constants in `SMALL_CAPS` live in `config.py`. Stage 12 §17.2 and
Stage 13 §14 are written as handovers to this document; both are answered in §0.3.

**Precondition.** Written against the landed Stages 1–13. This stage reads every result object they
return and computes no estimate, interval or p-value of its own. It fits nothing.

**Amendment, 2026-09-11 — eight names in this document are now public, and three are new.** The
manuscript exhibits in `figures_and_tables/` are a second reporting surface, and a [§16] sentence they
compose for themselves is a second definition of what the report says. So `_direction_clause`,
`_exceedance_clause`, `_interval_note`, `_pmf` and `_render_csv` lost their underscore, `_t02`'s and `_binary_rows`'
bodies were extracted as `baseline_rows` and `binary_rows`, `_t12`'s clause pair as
`arm_clauses` beside a new `AGREEMENT_NOT_REASSURANCE` constant, and §12.1's pinned public surface
grew from three names to eleven.
Nothing else moved: the row builders, the figure builders and `_render_md` stay private, `OUTPUT_IDS`
is unchanged, and `report.py` still knows nothing about the manuscript layer. Read every `_name` below
as its promoted spelling where the two differ.

**Status.** Written 2026-08-29; implemented the same day, and §16 records what was measured. Where
the implementation departed from the first draft the section says so in a **Landed** note. No DECISION is taken:
[§16]'s four amendments (a–d) already prescribe every statement this stage prints. DECISION 10's
reversible reporting choice — the eligible share and the undiluted eligible contrast as bootstrap keys
with intervals — stands as keys.

**Goal.** One entry point that runs Stages 1–13 in the canonical order, writes every table, figure and
log under `../out/`, and emits a run summary and the [§16] manuscript checklist. Every labelling rule
[§16] imposes is derived from a result object, so an edit to a template cannot drop it.

**Not in scope.** Manuscript prose. Any new estimator or estimand. Any change to a number Stages 1–13
produce. Figure styling beyond defaults.

---

## 0. Where Stage 14 sits

```
  report.run(source)                                       report.write(run, out)
  ┌──────────────────────────────────────────────┐        ┌─────────────────────────┐
  │ data.load → derive → eligibility → cohort    │        │ out/tables/T01..T22.md  │
  │ → propensity → balance → outcome.primary     │  Run   │              + .csv     │
  │ → outcome.secondary → bootstrap.run          │ ─────► │ out/figures/F01..F05.svg│
  │ → sensitivity.{multiplicity, e_value_primary,│        │ out/logs/audit_<label>.md│
  │   subgroups, full_covariate}                 │        │ Manifest (paths, sha256)│
  │ → standardise.{population, all_centre,       │        └─────────────────────────┘
  │   support, all_centre(over), hierarchical,   │
  │   inference}                                 │
  │ → policy.{population, contrast, inference}   │
  └──────────────────────────────────────────────┘
```

### 0.1 One module, `report.py`

It holds results, renders and writes. It imports no fitter (§12.14). The three test drivers that
carried the call order (`tests/test_sensitivity.py`, `tests/test_standardise.py`, `tests/test_policy.py`)
now delegate to `report._stages_1_10`, `report._stage_12`, `report._stage_13` — the private pieces
`run` is composed of — so the order exists once. **Landed:** the drivers keep their `with_inference`
flag, which is why the pieces and not `run` are what they call.

### 0.2 What Stage 14 touches upstream

`data.py` gains two row-builders (§8). `standardise.py` and `policy.py` become their callers, and
`policy.py` gains one public string alias (`NO_CONTRAINDICATED`). `config.py` gains the additions in §11. No public surface of any Stage 1–13 module changes;
the audit log is byte-identical before and after (§12.12).

### 0.3 The handovers, answered

| Owed by | Rule | Where discharged |
|---|---|---|
| Stage 12 §17.2 | different-populations statement computed from `n_average` and `in_estimate.sum()`, printed as numbers | T22, F05 — `_populations_clause` (§6) |
| Stage 12 §17.2 | `exp(beta)` conditional wherever it appears; never under `Primary.odds_ratio` | T15/T18 print `measure` beside `conditional_odds_ratio`; T22 has no OR column for [§14] |
| Stage 12 §17.2 | support check beside the estimate | T16 is rendered into T22's footer, not a supplement |
| Stage 12 §17.2 | cutpoint-collapse rate beside `RD_5` | T15/T18 footnote from `Draws.failures` and the `rd_5 == rd_4` count |
| Stage 12 §17.2 | `at_floor` rate and USZ posterior SD for the hierarchical arm | T15 hier block footer |
| Stage 12 §17.2 | one `beta` per fit; nothing for `support.` | T15 prints `beta`, `hier.beta`; §12.6 asserts no `support.beta` cell |
| Stage 13 §14 | label printed from `Policy.label` | T18 title cell |
| Stage 13 §14 | three populations printed as numbers | T22 header row |
| Stage 13 §14 | `share_eligible` beside every `rd_k`; `eligible_rd_k` printed as a factor, never in a column with Stage 12's `RD_k` | T18 layout; §12.5 asserts T22 has no `eligible_rd` column |
| Stage 13 §14 | `δ` nuisance, fit table only, no interval | T18 fit block inherits `policy._fit_table`'s role column |
| Stage 13 §14 | no [§14b] p-value, no [§13] family | T18 has no p column; §12.5 |
| Stage 13 §14 | a cohort without contraindicated patients is reported, not refused | `run` has no skip path (V4) |

---

## 1. Deliverables

| Path | Contents |
|---|---|
| `report.py` | **new.** `Run`, `Manifest`, `run`, `write`, `main` |
| `tests/test_report.py`, `tests/fixtures_stage14.py` | **new.** §12 |
| `data.py` | **amended.** `coefficient_rows`, `replicate_rows` (§8) |
| `standardise.py`, `policy.py` | **amended.** call the row-builders; five pipe-in-cell sites fixed (§8) |
| `config.py` | **amended.** §11 |
| `tests/test_sensitivity.py`, `tests/test_standardise.py`, `tests/test_policy.py` | **amended.** drivers delegate to `report._stage*` (§0.1) |
| `implementation_roadmap.md`, `TODOS.md` | **amended.** §13, §15 task 8 |
| `statistical_analysis_plan.md` | **unchanged.** Nothing here needs an amendment |

---

## 2. Environment

Pins unchanged. `matplotlib>=3.9` is already a dependency; `report.py` sets the `Agg` backend before
importing `pyplot`. Run: `cd extended_bridging && uv sync && uv run python -m report`, or `uv run pytest -v`.
**Landed:** no `[project.scripts]` entry — the project has no build system by design (`pyproject.toml`),
so a console script cannot be installed; `python -m report` is the entry point.
Cost: one full run is three bootstraps of `N_BOOT` plus Stage 11's; measured in §16.

---

## 3. Module shape

### 3.1 What the stage returns

```python
@dataclass(frozen=True)
class Run:
    source: data.Source
    audit: data.Audit
    df: pd.DataFrame                     # classified frame, unchanged by every later stage
    cohort: pd.DataFrame                 # the [§3] cohort Stages 6–11 ran on
    ps: propensity.Propensity
    balance: balance.Balance
    primary: outcome.Primary
    secondary: outcome.Secondary
    boot: bootstrap.Bootstrap
    multiplicity: sensitivity.Multiplicity
    e_value: sensitivity.EValue
    subgroups: sensitivity.Subgroups
    arm: sensitivity.Arm
    std_population: pd.DataFrame
    std: standardise.Standardisation     # all eligible
    support: standardise.Support
    std_support: standardise.Standardisation
    hier: standardise.Hierarchical
    std_boot: bootstrap.Bootstrap
    pol_population: pd.DataFrame
    pol: policy.Policy
    pol_boot: bootstrap.Bootstrap

@dataclass(frozen=True)
class Manifest:
    paths: tuple[Path, ...]              # sorted; relative to `out`
    sha256: dict[str, str]               # relative path -> digest of the written bytes
```

`Run.__post_init__` raises V4 if any field is `None`. Nothing is recomputed from `Run`; it is a record.

### 3.2 The public names — three

```python
def run(source: data.Source = data.WORKBOOK) -> Run: ...
def write(run: Run, out: Path = C.OUT) -> Manifest: ...
def main() -> None: ...                  # run() then write(); prints the manifest
```

### 3.3 The output inventory

Ids are fixed here and in `C.OUTPUT_IDS`. Every id is written exactly once (V2). `Guardrails` names
the §6 clause each carries.

| Id | Slug | Source | Guardrails |
|---|---|---|---|
| T01 | cohort_flow | `audit.entry("cohort","cohort_flow")`, `outcome.estimation_population` | denominator chain [§11] |
| T02 | baseline_by_arm | `df`, `ps.w`, `ps.ess` | ATO description [§7] |
| T03 | crude_event_rates | `df`, `C.OUTCOMES` | DESCRIPTIVE_CRUDE |
| T04 | onset_to_groin_by_arm | `df` | descriptive [§12] |
| T05 | balance_smd | `balance.covariates` | — |
| T06 | overlap_by_centre | `balance.centres`, `balance.pooled` | structural centres named [§9] |
| T07 | primary | `primary`, `boot.intervals`, `balance.unbalanced()` | CONSTANT_SHIFT, direction, exceedance, denominator |
| T08 | mrs_distribution_by_arm | `primary.cumulative` | — |
| T09 | e_value | `e_value`, `primary.cumulative`, `balance.unbalanced()` | E_VALUE_HEURISTIC, range, exceedance |
| T10 | secondary | `secondary.by_family()["secondary"]`, `boot`, `multiplicity` | MODEL_ASSISTED, reduced spec, denominator |
| T11 | safety | `secondary.by_family()["safety"]`, `boot`, `multiplicity` | DESCRIPTIVE_SAFETY, absent rows, denominator |
| T12 | sensitivity_arm | `arm` | ESS, worst SMD, deferred rows |
| T13 | subgroups | `subgroups` | HYPOTHESIS_GENERATING, SINGLE_SHIFT, NON_COLLAPSIBLE, thinning |
| T14 | subgroup_cumulative | `subgroups.estimates[*].cumulative` | — |
| T15 | standardised_14a | `std`, `std_support`, `hier`, `std_boot` | TRANSPORTED, measure, collapse, at_floor |
| T16 | support_box | `support` | — |
| T17 | support_baseline | `support.baseline` | — |
| T18 | policy_14b | `pol`, `pol_boot` | `Policy.label`, share beside rd, factor role |
| T19 | bootstrap_counters | all five `Bootstrap`/`Subgroups`/`Arm` draw sets | denominator per block |
| T20 | run_summary | `Run`, `C.SEED`, `C.DATA_SHA256`, `len(audit)` | — |
| T21 | strobe_record_checklist | `C.CHECKLIST` | V3 |
| T22 | primary_vs_14 | `primary`, `std`, `pol`, `support` | POPULATIONS, TRANSPORTED |
| F01 | mrs_bars_primary | `primary.cumulative` | — |
| F02 | forest_binary | T10, T11 rows | MODEL_ASSISTED, DESCRIPTIVE_SAFETY in legend |
| F03 | forest_subgroups | T13 rows | HYPOTHESIS_GENERATING in title |
| F04 | overlap_by_centre | `ps.e`, `df["center"]`, `balance.centres` | structural centres omitted, named in caption |
| F05 | mrs_bars_primary_vs_14 | `primary`, `std`, `pol` | POPULATIONS, TRANSPORTED in caption |

---

## 4. The canonical call order

Written once, in `report.run`, in §0's box order. Stage 12's `support` precedes the restricted
`all_centre` because `over=sup.inside` is its output. The classified frame is passed unchanged to
Stages 12 and 13; §12.15 asserts it is not mutated. `run` writes nothing.

---

## 5. Outputs

### 5.1 Tables

Path `out/tables/<id>_<slug>.md` and `.csv`. Both are rendered from one `(header, rows)` pair: the
`.md` through `data._md_table`, the `.csv` through `csv.writer(lineterminator="\n")`. Cells are
strings already passed through `data._fmt`; the CSV writes the same strings. Row counts agree (§12.9).
Clauses (§6) are written below the table as one paragraph per clause, in the `.md` only; the `.csv`
carries them as a trailing `# ` comment line each, so a CSV opened alone is still labelled.

Per-table rules, where not obvious from §3.3:

- **T01.** Rows: the four `cohort_flow` rows, then three denominator rows — cohort, covariate-complete
  (`ps.in_model.sum()`), outcome-present (`primary.in_estimate.sum()`) — then one row per
  `BinaryOutcome` denominator (`in_estimate.sum()`). Per-centre columns from `CENTER_ORDER`.
- **T02.** One row per `BALANCE_SET` covariate: unweighted mean/proportion by arm, weighted by arm,
  then `n` and ESS per arm from `ps.ess`. Header states it describes the ATO population.
- **T03.** One row per `C.OUTCOMES` entry: events/denominator/rate by arm, pooled and within centre.
  Title cell and footer carry `DESCRIPTIVE_CRUDE`. No difference column — a difference is an estimate.
- **T07.** Row 1 `exp(β)`: point, `Interval.lo/hi`, `p`, denominator. Rows 2–7 `RD_k`: point, lo, hi,
  `p` cell `—` (none exists), denominator. Footer: CONSTANT_SHIFT always; `_direction_clause` when
  it applies; `_exceedance_clause`. Never a column headed `odds_ratio` for any [§14] object.
- **T09.** `e_value.approximation` verbatim; `_e_value_range`; `_exceedance_clause`. Both required to
  appear together [§16d].
- **T10/T11.** One row per outcome: `rd` with interval, `odds_ratio` (`or_corrected` flagged), raw
  `p`, adjusted `p`, denominator. T10 adds `augmented` where `augmented_path != "unaugmented"` in a
  column headed `augmented (model-assisted)`, and the reduced covariate list where `reduced` is set.
  T11 adds `descriptive` to every row's `label` column and one row per `FamilyCorrection.absent`
  with its reason; header states `m_declared` and `m_used`.
- **T13.** Per subgroup: `n_by_level_arm`, `or_level[0]`, `or_level[1]`, `or_ratio` with interval and
  `p`, `n_interaction_tests`. Footer: HYPOTHESIS_GENERATING, SINGLE_SHIFT, NON_COLLAPSIBLE, and
  `_thinning_clause(subgroups.draws)` where it fires.
- **T15.** Three blocks (all eligible, treated support, random intercept), each: `n_average`,
  `n_fit`, distribution by arm, six `rd_k` with intervals, `mrs_0_2`, `mortality`, `measure`,
  `conditional_odds_ratio` for the two fits that have one. Footer: TRANSPORTED; `rd_5 == rd_4`
  count from `std_boot.draws`; `hier.at_floor` rate and `posterior_sd["USZ"]`.
- **T18.** One table, `row` column: per `k` — `rd_k`, interval, `share_eligible`, `eligible_rd_k`
  under the heading `factor: eligible contrast under this model`; then `mRS 0-2`, `mortality`,
  `n_average`, `n_eligible`; then one row per coefficient with its role (from `coefficient_rows`, the
  indicator labelled a nuisance, no interval); then the conditional odds ratio beside `measure`.
  Footer: `pol.label`, TRANSPORTED, the factor sentence, the collapse count, and
  `policy.NO_CONTRAINDICATED` when `share_eligible == 1` (a public alias added for this one read).
- **T19.** Five blocks via `replicate_rows`, each with its denominator column: `n_attempted`,
  surviving draws, one column per `FAILURE_BUCKETS` value. Block labels: `primary and secondary`,
  `subgroups`, `sensitivity arm`, `[§14a]`, `[§14b]`.
- **T20.** Key/value rows: source label, `DATA_SHA256`, `N_RECORDS_EXPECTED`, `SEED`, every
  `n_boot` (one row per bootstrap object), `CI_LEVEL`, `PERCENTILE_METHOD`, total failures per
  bucket, audit ledger length, output count.
- **T21.** One row per `CHECKLIST` item: STROBE/RECORD item, output ids. V3 if any id list is empty.
- **T22.** Header row prints the three populations as numbers: `primary.in_estimate.sum()` over
  `len(cohort.treating_centres)` centres; `std.n_average` over four; `pol.n_average` over four with
  `pol.n_average - pol.n_eligible` contraindicated. Rows `rd_k`, `mrs_0_2`, `mortality` with
  intervals for each of the three. No OR column. Footer: POPULATIONS, TRANSPORTED, support numbers
  (`support.n_never_ivt`, `n_never_ivt_outside`, largest non-grouping baseline |SMD|).

### 5.2 Figures

SVG only. `plt.rcParams.update(C.SVG_RC)`; `savefig(..., format="svg", metadata={"Date": None,
"Creator": None})`. Text is not converted to paths. Captions are `<title>` elements set from the
clause constants, so §12.5 can assert on the file bytes.

- **F01.** Stacked horizontal bars, one per arm, seven `MRS_LEVELS` segments from
  `primary.cumulative`; arm labels from `TREATMENT_LABELS`.
- **F02.** Forest of `rd` intervals for every binary outcome; secondary above safety; safety
  markers hollow, legend entry `descriptive`; augmented points as a second marker labelled
  `model-assisted`.
- **F03.** Forest of `or_level[0]`, `or_level[1]` per subgroup; title carries HYPOTHESIS_GENERATING.
- **F04.** One panel per centre in `CENTER_ORDER` whose `CentreOverlap.status` equals the reported
  literal (matched by value, `balance.py:94`); histograms of `ps.e` by arm; caption names omitted
  centres and their status.
- **F05.** F01's bars for the primary, `std`, `pol` stacked in three rows; caption POPULATIONS +
  TRANSPORTED.

### 5.3 The log

`write` calls `run.audit.write()` once, after recording `report_written` (§10). No other code path
writes the log.

---

## 6. Guardrails

Statements are module-level `Final[str]` constants. Conditional clauses are functions of a result
object and return `""` when they do not apply; a renderer concatenates, never branches on data.

| Constant / function | Source of truth | Text (shape) | Rendered in | Test |
|---|---|---|---|---|
| `CONSTANT_SHIFT` | [§16a] | common OR rests on an untested constant-shift assumption | T07 | 12.3 |
| `_direction_clause(rd)` | `Primary.rd` signs | present iff `len({sign(v) for v in rd.values()}) > 1` | T07 | 12.3 |
| `DESCRIPTIVE_CRUDE` | [§16b] | crude rates are descriptive, not unadjusted estimates | T03 | 12.8 |
| `_thinning_clause(draws)` | `Draws.failures["separation"]` | present iff any key's count > 0; states surviving count and selection on the estimate | T13, T19 | 12.4 |
| `E_VALUE_HEURISTIC` | [§16d] | `e_value.approximation` + heuristic sentence | T09 | 12.6 |
| `_e_value_range(cumulative)` | control-arm `P(Y<=k)` vs `E_VALUE_PREVALENCE_FLOOR` | names thresholds whose `min(p, 1-p)` is below the floor | T09 | 12.6 |
| `_exceedance_clause(names, covariates)` | `Balance.unbalanced()`, `CovariateBalance.weighted` | names each covariate with its magnitude | T07, T09 | 12.6 |
| `MODEL_ASSISTED` | `BinaryEstimate.augmented_path` | column header + legend | T10, F02 | 12.6 |
| `DESCRIPTIVE_SAFETY` | `Outcome.family == "safety"` | row label | T11, F02 | 12.6 |
| `HYPOTHESIS_GENERATING`, `SINGLE_SHIFT`, `NON_COLLAPSIBLE` | [§13], `SubgroupEstimate.hypothesis_generating` | footer | T13, F03 | 12.6 |
| `TRANSPORTED` | [§16], `Standardisation.population`, `Policy.label` | model-based transported results | T15, T18, T22, F05 | 12.5 |
| `_populations_clause(primary, std, pol)` | `in_estimate.sum()`, `n_average` | different populations, printed as numbers | T22, F05 | 12.5 |
| denominator column | `in_estimate.sum()` per object | every estimate row | T07, T10, T11, T13 | 12.7 |

`exp(β)` from `Standardisation` or `Policy` is `conditional_odds_ratio` and is printed only beside
`measure`. The attribute `odds_ratio` is deliberately absent there (`standardise.py:104`), so a
renderer reaching for it fails rather than mislabels.

---

## 7. Determinism

Inherited from `data.py:120-127`: no clock, no absolute path, no library repr, no set iteration, one
float formatter, sorted identifiers. Added: SVG metadata stripped, `svg.hashsalt` fixed in `SVG_RC`,
CSV line terminator fixed. Acceptance: `Manifest.sha256` identical across two `PYTHONHASHSEED`
values (§12.11).

---

## 8. Row-builder extraction

Closes the `TODOS.md` row-builder item on its stated trigger. Two functions in `data.py`:

```python
def coefficient_rows(fit, roles: Mapping[str, str], labels: Mapping[str, str]) -> tuple[tuple[str, ...], ...]
def replicate_rows(blocks: Sequence[tuple[str, str, Sequence[tuple[str, ...]]]]) -> tuple[tuple[str, ...], ...]
    # blocks: (block label, denominator label, rows); pads rows to a common width with "—"
```

Both raise V1 on a `|` in any cell. Callers: `standardise._fit_table`, `policy._fit_table`,
`standardise._replicates_table`, `policy._replicates_table`. **Landed:** `bootstrap.py` is NOT a
caller and `_padded` stays — `test_bootstrap.py` pins the module's private list, and Stage 10's grid
has a different shape (counters + padded diagnostics). Signatures as landed:
`coefficient_rows(columns, beta, dropped, role)` and `replicate_rows([(denominator, [(quantity, value)])])`.

**Order matters.** Stage 12's tables had `|` inside cells, so the log was malformed and byte-identity
would have frozen the defect. **Landed — five sites, all in `standardise.py`:** `_distribution_table`'s
header (`P | arm` → `P(Y=j) arm`) and its two `abs(...)` rows, `_fit_table`'s `|coef|` → `abs_coef`,
`_baseline_table`'s `|SMD| > 0.1` → `abs_SMD > 0.1`, `_replicates_table`'s `max |coefficient|` →
`max_abs_coef`. Task 2 was: fix these, run once at `N_BOOT`, capture the log's SHA-256
(`STAGE13_AUDIT_SHA256`), extract, run again, assert. `STAGE12_ENTRIES_SHA256` in `test_standardise.py`
moved with the header fixes and was re-captured from the same pre-extraction run;
`STAGE12_INTERVALS_SHA256` did not move.

---

## 9. Failure taxonomy — V series

| id | raised by | class | meaning |
|---|---|---|---|
| V1 | `data.coefficient_rows`, `data.replicate_rows` | `SchemaError` | a cell contains `\|` |
| V2 | `report.write` | `SchemaError` | an output id rendered twice |
| V3 | `report.write` | `SchemaError` | a `CHECKLIST` item maps to no output |
| V4 | `Run.__post_init__` | `SchemaError` | a field is `None` — a stage was skipped |
| V5 | `report.main` | `SchemaError` | manifest digests differ from a prior manifest passed via `--check` |

None is caught. No `FitError` exists here; nothing is fitted.

---

## 10. Audit entries — one

| kind | step | contents |
|---|---|---|
| `provenance` | `report_written` | n = output count; table of relative path and sha256 |

Ledger: Stage 13's 61 + 1 = **62**. `KINDS` unchanged. `report_written` is disjoint from every
existing step (§12.1).

---

## 11. Config additions

```python
# --- [§16] reporting ---
E_VALUE_PREVALENCE_FLOOR: Final[float] = 0.15   # [§16d] VanderWeele's documented licence
OUTPUT_IDS: Final[dict[str, str]] = {
    "T01": "cohort_flow", "T02": "baseline_by_arm", "T03": "crude_event_rates",
    "T04": "onset_to_groin_by_arm", "T05": "balance_smd", "T06": "overlap_by_centre",
    "T07": "primary", "T08": "mrs_distribution_by_arm", "T09": "e_value",
    "T10": "secondary", "T11": "safety", "T12": "sensitivity_arm", "T13": "subgroups",
    "T14": "subgroup_cumulative", "T15": "standardised_14a", "T16": "support_box",
    "T17": "support_baseline", "T18": "policy_14b", "T19": "bootstrap_counters",
    "T20": "run_summary", "T21": "strobe_record_checklist", "T22": "primary_vs_14",
    "F01": "mrs_bars_primary", "F02": "forest_binary", "F03": "forest_subgroups",
    "F04": "overlap_by_centre", "F05": "mrs_bars_primary_vs_14",
}
SVG_RC: Final[dict[str, object]] = {"svg.hashsalt": "stage14", "svg.fonttype": "none",
                                    "figure.dpi": 100, "path.simplify": False}
# STROBE + RECORD items that name an output. Item text is the checklist's own; ids are OUTPUT_IDS keys.
CHECKLIST: Final[tuple[tuple[str, tuple[str, ...]], ...]] = (
    ("STROBE 13 participants: numbers at each stage", ("T01",)),
    ("STROBE 14 descriptive data", ("T02", "T04")),
    ("STROBE 15 outcome data", ("T03",)),
    ("STROBE 16 main results", ("T07", "T08", "F01")),
    ("STROBE 16 other analyses: secondary, safety", ("T10", "T11", "F02")),
    ("STROBE 17 other analyses: subgroups, sensitivity", ("T12", "T13", "T14", "F03")),
    ("STROBE 12 statistical methods: balance and overlap", ("T05", "T06", "F04")),
    ("STROBE 19 limitations: unmeasured confounding", ("T09",)),
    ("RECORD 12.3 linkage and cleaning", ("T20",)),
    ("[§14] transported analyses", ("T15", "T16", "T17", "T18", "T22", "F05")),
    ("[§10] bootstrap and seed", ("T19", "T20")),
    ("[§16] manuscript checklist emitted as a file", ("T21",)),
)
```

`OUTPUT_IDS` is the one list; `write` iterates it, V2 guards duplicates, and §12.2 asserts every id
lands.

---

## 12. Acceptance criteria

### 12.0 `tests/fixtures_stage14.py`

A synthetic `Run` built from dataclasses without fitting: `Primary` with one-signed `rd` and a
second with mixed signs; `Bootstrap` draw sets with zero and non-zero `separation` counts; a
`Balance` with one exceedance; `Secondary` with one augmented and one unaugmented estimate and one
absent safety p. Real-data tests are `DATA_GATED`.

### 12.1 Surface
Three public names, two dataclasses, source order `run, write, main`. `report_written` not in any
other module's step set. `report.py` reaches into no other module's privates except `data._fmt`,
`data._md_table`.

### 12.2 Coverage
`write` produces every `OUTPUT_IDS` id once; every `Manifest.paths` entry exists; `.md` and `.csv`
counts equal `len(OUTPUT_IDS)` for `T*`, `.svg` for `F*`.

### 12.3 The direction clause tracks the data
On the one-signed fixture T07 contains `CONSTANT_SHIFT` and not `_direction_clause`'s text; on the
mixed fixture it contains both. The roadmap's accept-when.

### 12.4 The thinning clause tracks the counters
Absent when every `separation` count is zero; present with the surviving count otherwise.

### 12.5 [§14] beside the primary
T22 and F05 bytes contain `POPULATIONS` and `TRANSPORTED`; T07 contains neither. T22 has no
column named `odds_ratio`, `eligible_rd`, or `support.beta`. T18 has no `p` column.

### 12.6 Labels are read, not typed
Every T10 row with `augmented_path != "unaugmented"` carries `MODEL_ASSISTED`; every T11 row carries
`DESCRIPTIVE_SAFETY`; T13 and F03 carry `HYPOTHESIS_GENERATING`; T09 contains
`e_value.approximation` and every exceedance name with its magnitude, and names exactly the
thresholds whose control-arm baseline is below `E_VALUE_PREVALENCE_FLOOR`.

### 12.7 Denominators
Every row of T07, T10, T11, T13 has a non-empty `denominator` cell equal to the object's
`in_estimate.sum()` or `n`.

### 12.8 Crude table
T03 title and footer contain `DESCRIPTIVE_CRUDE`; T03 has no difference column.

### 12.9 Rendering
No cell of any table contains `|`; header, separator and every body row have equal cell counts;
CSV row count equals `.md` body row count + 1.

### 12.10 Checklist
Every `CHECKLIST` id is in `OUTPUT_IDS`; V3 fires on a fixture item with an empty id tuple.

### 12.11 Determinism
Two `write` calls under different `PYTHONHASHSEED` yield identical `Manifest.sha256`.

### 12.12 Upstream log unchanged
`STAGE13_AUDIT_SHA256`, captured after the header fix and before the extraction (§8), equals the
digest after. `DATA_GATED`; the gate is a run with zero skips.

### 12.13 Nothing leaks into version control
All `Manifest.paths` are under `C.OUT`; a grep of `specs/stage14_*.md` and `report.py` for a float
with four or more decimals returns nothing.

### 12.14 No fitter
AST scan: `report.py` imports none of `model`, `scipy`, `statsmodels`; calls nothing from
`propensity`, `outcome`, `bootstrap`, `sensitivity`, `standardise`, `policy` except the entry points
`run` names.

### 12.15 Upstream suites
Stages 1–13 test files unedited except the three driver replacements; all green; `run` returns
`df` equal to the classified frame passed to Stages 12 and 13.

---

## 13. Known gaps carried forward

1. **Figures have no independent oracle.** Assertions are on bytes and captions, not on geometry.
   **Trigger:** a reviewer disputing a figure.
2. **The checklist maps items to outputs; it writes no prose.** **Trigger:** manuscript drafting.
3. **Hybrid plug-in for contraindicated patients — closed negatively.** PI confirmation 2026-08-29:
   the model prediction stands; T18 prints `Policy.label` and `measure`. The `TODOS.md` item closes
   with this document.

## 14. NOT in scope

| Considered | Why declined |
|---|---|
| Word/LaTeX tables | new dependency; Markdown + CSV suffice for transcription |
| PNG figures | not byte-reproducible |
| A `--skip` flag for Stages 12–13 | a partial run would produce T22 without its comparators; V4 |
| A proportional-odds test | [§15] declines it; DECISION 5 |

## 15. Implementation tasks

In order, each landing green before the next begins.

1. `report.run` and `Run`; replace the three test drivers with it.
2. Fix `standardise._distribution_table` header; run; capture `STAGE13_AUDIT_SHA256`; add
   `data.coefficient_rows`, `data.replicate_rows`; move the six callers; assert the digest.
3. `config.py` additions (§11).
4. Tables T01–T22 with clause constants and functions (§6).
5. Figures F01–F05.
6. `fixtures_stage14.py`, `test_report.py` (§12).
7. `write`, `Manifest`, `report_written`, `main`, `pyproject.toml` script.
8. Documents: roadmap Stage 14 marked landed; `TODOS.md` — close row-builder, denominator-label and
   hybrid plug-in items; open the figure-oracle item.

## 16. Verification record

Measured 2026-08-29 on the workbook at `N_BOOT` = 2000: `run` 571 s; 27 outputs, 50 files plus the log;
ledger 61 before `report_written`, 62 after. `Manifest.sha256` identical across `PYTHONHASHSEED` 0 and 1
(60-replicate run). Two `write` calls of one `Run` are byte-identical (the entry is replaced, not
appended). Clauses on the workbook: `_direction_clause` FIRES (the six RD_k do not share a sign);
`_thinning_clause` FIRES on the `unknown_onset` subgroup and on no other draw set; `_exceedance_clause`
names five covariates; `_e_value_range` names thresholds at both tails. The pre-extraction log digest
equals the post-extraction one (§12.12, run twice at 2000). Not claimed: figure geometry (TODOS).
