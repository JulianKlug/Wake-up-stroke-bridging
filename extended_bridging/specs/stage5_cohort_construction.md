# Stage 5 spec — cohort construction

Implements roadmap Stage 5 [§2, §3]. Section references in brackets are to
`statistical_analysis_plan.md`. Numbers and decisions referenced as DECISION *n* are established in
Stage 0 and recorded in `../out/stage0_data_inventory.md`. Constants named in `SMALL_CAPS` are Stage
1's and live in `config.py`; `stage1_config_and_data_contract.md` is their specification. The frame
this stage receives is specified by `stage4_eligibility_classification.md` §11.

**Status.** Written 2026-08-13, against the frame and the workbook as verified in §18. Reviewed the
same day; §21 is the review record and every finding it raised is folded into the sections below. It
is the **sole source for the Stage 5 implementation**: everything the implementer needs is here, and
anything not here is not to be invented.

**Written against DECISION 1a, not DECISION 1.** On 2026-08-13 the PI amended the eligibility rule:
`ivt_contraindicated` is the classifier and the free text is not read at all, so **eligibility has
two classes — a patient whose contraindication reason was never documented is `eligible`**. The
`indeterminate` class is abolished.

**That amendment has landed. Stage 5 is unblocked.** All nine items §20 tracks are in the tree as of
2026-08-13 — `config.py`, `eligibility.py`, [§3], `../out/stage0_data_inventory.md`, Stage 4's spec
§21, the three test files and the roadmap — and `uv run pytest -v` is green at **490 tests**. §20 is
now the dated record of what the amendment required and where each item landed, not a list of work to
do. What the amendment did *not* change is the cohort: the retained set is the same 107 records and
the primary cohort the same 93 patients, identifier for identifier (§18). Every flow number in this
document is unaffected; only the class columns of §7's tables are.

**Goal.** The primary cohort exists as a frame, produced by two restrictions applied in a declared
order, each naming the patients it removed; the analysis population is asserted to be one no
estimator can quietly fail on; and the [§13] median subgroup is frozen on that cohort rather than on
the frame it came from.

**Not in scope.** Fitting anything (Stage 6), any balance or overlap diagnostic (Stage 7), the
all-centre population of [§14a] (Stage 12), and the policy population of [§14b] (Stage 13). Stage 5
**removes rows and adds one column**. It edits no value that Stages 2, 3 or 4 delivered.

**One thing this spec settles that no earlier stage could.** Stage 4 §9.1 deferred to Stage 5 the
question of whether a cohort-shaped audit `kind` is added to `data.py`. It is: §6 declares `cohort`,
in the position §6.2 argues for, and `_MUST_NAME_CASES` is deliberately left alone for the reason
§6.3 gives. That is the first `data.py` amendment since Stage 2.

**And one finding that reshapes the test suite before a line is written.** Neither
`tests/fixture_schema.xlsx` nor `test_data.hand_frame()` can reach the end of this stage: both
collapse to a cohort with no control arm, which §5.2's guard raises on. Measured, not anticipated
(§18). Stage 5 therefore declares its own frame, and the two existing ones become *positive* tests of
the guard rather than fixtures it runs on. §12.0 is that argument in full; skipping it produces a
test file that cannot be written.

---

## 0. Where Stage 5 sits

```
  data.load()  →  (df, audit)          25 analysis columns, 126 rows   [Stage 2 §11]
                       │
                       ▼
  derive(df, audit)                    31 columns                      [Stage 3 §11]
                       │
                       ▼
  classify(df, audit)                  32 columns                      [Stage 4 §11]
                       │
                       ├──────────────────────────────────→  Stage 12 [§14a], Stage 13 [§14b]
                       │                                     the UNRESTRICTED frame, 126 rows
                       ▼
  ┌────────────────────────────────────────────────────────────────────────┐
  │  STAGE 5 — cohort.py           reads config.py, data.Audit, eligibility │
  │                                                                        │
  │   build(df, audit)                                                     │
  │     ├─ _assert_cohort_inputs(df)    C1…C4 → SchemaError                │
  │     ├─ restriction 1, inline        [§3] centres          126 → 104    │
  │     ├─ restriction 2, inline        [§3] eligibility      104 →  93    │
  │     │                               both masks are written out in §8;  │
  │     │                               neither is a function of its own   │
  │     ├─ _flow_table(…)               one row per step, + the            │
  │     │                               no-reason-on-file line             │
  │     ├─ derive.derive_cohort(…)      Stage 3's second entry point,      │
  │     │                               HERE and nowhere else  32 → 33     │
  │     ├─ absence_by_column(…)         the cohort's own denominators [§11]│
  │     └─ _assert_cohort(out)          P1…P4 → SchemaError                │
  │                                                                        │
  │   treating_centres(df)              the [§3] restriction-1 predicate,  │
  │                                     computed, never declared           │
  └────────────────────────────────────────────────────────────────────────┘
                       │
                       ▼
              33 columns, 93 rows  →  Stages 6-11
```

`build` appends four entries to the `Audit` that `load()` created — three of the new `cohort` kind
and one `missingness` — plus the two `derivation` entries `derive_cohort` records. It writes no file.

### 0.1 Why a module of its own

- **`derive_cohort` needs exactly one caller and this is it.** `derive.py`'s docstring already says
  so — "**Stage 5 is `derive_cohort`'s only caller**" — and today that is a sentence rather than a
  structure. §4.5 makes it a structure: `build` is the only function that calls it, so the guarantee
  "the [§13] median is the cohort's" stops depending on a driver's discipline. The numbers are why it
  matters: the core-volume median is 6.0 mL over the 126-record frame and 5.0 mL over the 93-record
  cohort (§18), so a call placed one line too early splits the cohort on the wrong number, and both
  numbers are plausible.
- **Stage 12 and Stage 13 need the frame Stage 5 throws rows away from.** They take the unrestricted
  126-record classified frame, which is the second arrow out of the diagram above. A restriction
  living inside a Stage 6 propensity function would be recomputed by every caller that wanted a
  cohort and by none that wanted the frame.
- **The repository's unit is one module, one spec, one test file.** Stages 1-4 each hold that shape.

### 0.2 What Stage 5 does not touch

`build` calls `df.copy()` before writing, as `derive` and `classify` do, so a caller's frame is never
mutated underneath it. No existing column is edited and no value is corrected. What it *does* do,
and no earlier stage did, is **drop rows** — 33 of 126 — which is why §7's three tables and §12.14's
promise assertions carry more weight here than their Stage 4 equivalents.

**The index is preserved and never reset.** A cohort record keeps the label it carried on the
unrestricted frame, so it can be traced back to the frame Stages 12 and 13 analyse. `reset_index(drop=True)`
would destroy that link and nothing would fail — the cohort would be the same 93 patients under new
labels. `pilots/analysis.py:184` does exactly that, and §16 is where it is declined. §12.14 asserts
the index is preserved and unique.

## 1. Deliverables

| Path | Contents |
|---|---|
| `extended_bridging/cohort.py` | `build`, `treating_centres` |
| `extended_bridging/test_cohort.py` | the acceptance tests in §12, the `cohort_frame()` builder of §12.0, and its **own** module-scoped `workbook` fixture — `test_data.py`'s and `test_eligibility.py`'s do not cross files |
| `extended_bridging/data.py` | **amended** — `KINDS` and `_HEADINGS` gain `cohort`; `_MUST_NAME_CASES` deliberately unchanged (§6.3); and three prose corrections DECISION 1a's sweep missed — A9's runtime message, the `_missingness` comment block, and the `_MUST_NAME_CASES` comment's kind count (§10, §21 R7) |
| `extended_bridging/test_data.py` | **amended** — the two literal pins of `KINDS` and `_HEADINGS.values()`, and the name of the test that carries the word *seven* (§10) |
| `extended_bridging/test_eligibility.py` | **amended** — `test_data_kinds_is_unchanged_by_this_stage`, whose comment predicted this stage and whose literal is now wrong (§10) |
| `extended_bridging/test_derive.py` | **amended** — the tail-slice repair `TODOS.md` defers to "the commit that first builds the real pipeline order". That commit is this one (§10) |
| `extended_bridging/implementation_roadmap.md` | **amended, and landing with this document** — Stage 5 gains its `**Spec:**` line and an **Accept when** naming the both-arms guard and the frozen-median ordering (§19) |
| `TODOS.md` | **amended twice** — the `tests/fixture_cohort.xlsx` item was **added at review**, ahead of the implementation (§13, §19); the tail-slice item is struck by the implementation, because §10 lands it (§19) |

`out/logs/audit_<label>.md` is an **output, not a deliverable**, as in Stages 2, 3 and 4: gitignored,
written only when a caller asks, and it names patients. Nothing under `specs/` may quote a case
identifier, and nothing in this file does — every patient below is a `HAND-N` or `COHORT-N` fixture
record, or a count.

`config.py` is **not** amended. Stage 5 declares no new constant: the restriction predicates are
computed from `CENTER_ORDER`, `TREATMENT`, `TREATMENT_LABELS` and `ELIGIBILITY_RETAINED`, all of
which exist. §15 records the one constant that was considered and declined.

## 2. Environment

Unchanged from Stage 2 §2.1. `uv`, Python 3.12, no new dependency — `cohort.py` needs pandas and
nothing else. Commands run from `extended_bridging/`, flat module layout, `import config as C`.

```bash
cd extended_bridging && uv sync && uv run pytest -v
```

**What this costs to run.** Two boolean masks over at most 126 rows, four boolean reductions for the
preconditions, `len(CENTER_ORDER) × len(TREATMENT_LABELS)` counts for the first table, three
reductions per centre for the second, and one pass per flow row for the third. `derive_cohort` adds
one median and one `nunique` per covariate. There is nothing to cache and §15 says so as a standing
instruction.

## 3. Module shape

```python
# cohort.py
from __future__ import annotations

import pandas as pd

import config as C
import derive
import eligibility
from data import Audit, absence_by_column
```

**No new types**, for Stage 3 §3's reason: everything this stage produces is a frame, a column or an
audit entry, and all three already have types. The public surface is two functions:

```python
def build(df: pd.DataFrame, audit: Audit) -> pd.DataFrame
def treating_centres(df: pd.DataFrame) -> tuple[str, ...]
```

and the privates are `_assert_cohort_inputs`, `_record_removal`, `_assert_cohort`, and §7's three
detail/table helper pairs — `_centres_detail` / `_centres_table`, `_eligibility_detail` /
`_eligibility_table`, `_flow_detail` / `_flow_table` — plus the `_ABSENCE_DETAIL` constant of §7.4.

**There is no `_restrict_centres` and no `_restrict_eligibility`.** Both restrictions are plain masks
written inline in `build` (§8), and §8 is canonical: a reader traces both restrictions, both audit
entries and the `derive_cohort` call through twelve consecutive lines, which extracting them would
break for no gain — the extracted functions would have to thread `after_1` and `after_2` back out
through return values. The decision was taken at review on 2026-08-13 (§21 R2). Anything in this
document that reads as though those two functions exist is a defect; report it.

`treating_centres` is public for one concrete reason and not for symmetry: **Stage 12's support check
needs its complement.** [§14a] requires "the proportion of patients from never-IVT centres falling
outside [the treated support]" and "a baseline table comparing IVT-treated patients with
never-IVT-centre patients". Both range over the centres this function does *not* return. A second
implementation over there would be a second definition of the [§3] restriction-1 population, and the
two would agree on v7 and be free to drift on the workbook that follows it.

**`cohort.py` imports `derive` and `eligibility`, and that is a first for this pipeline** — Stages 2,
3 and 4 each import only `config` and `data`. It is deliberate and it is the whole of §0.1's first
bullet: `build` calls `derive.derive_cohort` so that nothing else has to, and calls
`eligibility.retained` so that restriction 2 cannot be written as a comparison. §12.6's AST scan is
what keeps the second promise.

`cohort.py` is **not** exempt from the Stage 1 §7 raw-name scan and must never become exempt.
`test_config.py` 9.4 walks `extended_bridging/**/*.py` minus a pinned three-file exemption set, so
this module and its test file are scanned from the moment they exist. Neither names a raw header.

### 3.1 The pandas facts this stage turns on

Declared once here, verified by running them on pandas 2.3.3 as pinned by `uv.lock` (§18).

```
  df[mask]               with <NA> in mask       → the row is NOT selected
                                                   (the same fact Stage 4 §3.1 records for .loc)

  Series.isin(values)    with <NA> in the series → False, NEVER <NA>
                                                   so ~isin SELECTS the missing row

  df[mask] preserves the index                   → cohort labels are the frame's labels
  df[mask] on a `string` column leaves no        → no dead level; `center` is not a Categorical,
      dead level                                   per Stage 3 §4.4
```

**The first fact is why the preconditions come before the masks and not after.** A record whose
exposure or centre is missing is *silently excluded* by `df[mask]` — it is not selected by
restriction 1's keep-mask, so it leaves the cohort with no entry anywhere saying it did. The removal
would be attributed to restriction 1 in the flow table, which would still reconcile, and the
identifier would appear in restriction 1's `case_ids` as though it had been at a never-treating
centre. §4.4's C3 and C4 are what make that unreachable, and they are the same two checks Stage 4
made for the same reason at E4 and E5.

**The fourth fact is why nothing here calls `remove_unused_categories`.** `center` and `onset_type`
are `string`, by Stage 2 §6's and Stage 3 §4.4's decisions, taken precisely so that this stage's
restriction leaves no dead level behind for `groupby(observed=False)` to resurrect as an all-missing
row in every downstream table. `pilots/analysis.py:185` needs that call because it made `center` a
`Categorical`; §16 is where that is declined.

## 4. The two restrictions [§3]

### 4.1 The rule

[§3], restated once so the code has something to be checked against:

```
  1.  Exclude centres that contributed no bridging patients.
      There, P(IVT = 1 | centre, X) = 0 STRUCTURALLY — not by chance. No weighting
      recovers a contrast that was never available, and a penalised score would report
      shrinkage as though it were treatment availability.

  2.  Exclude patients with an absolute contraindication to IVT.
      They were never candidates for bridging; retaining them makes "no IVT" a marker
      of contraindications and their prognosis.

      Eligibility is TWO classes, from `ivt_contraindicated` alone [DECISION 1a]:
      flag = 1 -> ineligible; flag = 0 -> eligible. A patient whose contraindication
      reason was never documented is ELIGIBLE, not a third class. The free text of
      `Contraindications_to_IVT` is not read, and — unlike under DECISION 1 — neither
      is its presence.
```

Both are applied **by design, before any weighting** — which is the whole of why they live in a
stage of their own, ahead of Stage 6, rather than as an argument to a propensity fit.

### 4.2 Restriction 1 is computed, never declared

```python
def treating_centres(df: pd.DataFrame) -> tuple[str, ...]:
    """The centres in `df` that contributed at least one bridging patient, in CENTER_ORDER.

    [§3] restriction 1 keeps these and drops every other declared centre. Computed from the frame
    rather than declared as a constant: a centre's treatment availability is a property of the data,
    and the one edit that would matter — a workbook in which USZ starts administering IVT, or in
    which Lugano's two bridging patients turn out to be miscoded — must change the cohort rather
    than require someone to remember to change a list. `pilots/pilot_config.py:78` declares
    `NEVER_IVT_CENTERS = ["USZ"]`; §16 is where that is declined.

    Ranges over CENTER_ORDER and never over the frame's observed values, so the result is a declared
    order and the log is deterministic. A declared centre contributing no rows at all is therefore
    not returned — vacuously correct, since it contributed no bridging patient, and visible in
    §7.1's table as an all-zero row rather than as an absence.
    """
    return tuple(c for c in C.CENTER_ORDER
                 if int(((df["center"] == c) & (df[C.TREATMENT] == 1)).sum()) > 0)
```

Three properties, each load-bearing:

- **`> 0`, over the bridging arm only.** [§3] restriction 1 names one direction, because that is the
  direction v7 showed. The other direction — a centre with no *control* — is a positivity violation
  too, and §5.2 raises on it rather than silently extending [§3]. The asymmetry is [§3]'s and Stage 5
  does not resolve it.
- **The arm code is `1`, from `C.TREATMENT`, and the count is over the frame given.** Not from
  `TREATMENT_LABELS`, which is a display registry; §7.1's table is where the labels are read.
- **It returns the *kept* centres, not the dropped ones.** The dropped set is
  `tuple(c for c in C.CENTER_ORDER if c not in keep)`, computed at the one place that needs it. The
  positive form is what Stage 12 wants the complement of, and one direction declared is one fewer
  thing to keep consistent.

### 4.3 Restriction 2 is `eligibility.retained`, and nothing else

```python
    keep = eligibility.retained(df)          # != ineligible, from ELIGIBILITY_RETAINED
```

Never a comparison written here. **Under DECISION 1a the two readings coincide** — with two classes,
`== C.ELIGIBLE` and `!= C.INELIGIBLE` select the same patients — so the 43-record hazard Stage 4 §5.1
was built against is gone, and this section is deliberately *weaker* than the draft it replaces. It
is not deleted, for two reasons:

- **`retained` remains the one seam a third class would come back through.** It reads
  `ELIGIBILITY_RETAINED`, which is where a class would be declared. If the [§13] sensitivity suite is
  ever un-deferred with an "undocumented reason" arm — the analysis DECISION 1a makes newly
  interesting, since the 43 are now inside the primary cohort rather than flagged in it — the two
  readings come apart again, at a stage that has already been written.
- **A literal is still a second declaration.** §12.6 keeps the AST scan for any string equal to a
  member of `ELIGIBILITY_ORDER`, and for a `Subscript` of `ELIGIBILITY_ORDER` by an integer constant,
  which would couple the restriction to a display order. Both have companion tests proving they fire.

**What is dropped is the `Compare` scan and its justification.** The earlier draft forbade
`df[df[C.ELIGIBILITY] != C.INELIGIBLE]` on the grounds that 43 of 93 cohort records turned on it.
They no longer do — that construct is now merely redundant, not dangerous — and a scan whose stated
reason has evaporated is worse than no scan, because the next reader believes the reason.

### 4.4 The preconditions come first [§3]

Before any row is removed, and before the frame is copied. Four checks, collected and raised
together, following Stage 2's `_assert_schema` and Stage 4's `_assert_classifier_inputs`:

```python
def _assert_cohort_inputs(df: pd.DataFrame) -> None:
    bad: list[str] = []

    if C.ELIGIBILITY not in df.columns:
        bad.append(
            f"C1  {C.ELIGIBILITY}: absent. Stage 4's classify() has not run on this frame. [§3]'s "
            "second restriction is a function of eligibility and this stage will not classify: "
            "eligibility.classify() is one call and doing it here would give the column two "
            "producers and the log two entries.")

    already = [c for c in C.COHORT_DEPENDENT_SUBGROUPS if c in df.columns]
    if already:
        bad.append(
            f"C2  {', '.join(already)}: already present. The [§13] median subgroup is frozen on the "
            "COHORT [Stage 3 §6.4], so a frame that already carries it was restricted once already. "
            "derive_cohort raises on this too — but only after both restrictions have been applied "
            "and logged, leaving an audit log describing a cohort no caller received.")

    off_centre = df.loc[~df["center"].isin(C.CENTER_ORDER), "case_id"]
    if len(off_centre):
        bad.append(
            f"C4  center: {len(off_centre)} record(s) carry a centre outside CENTER_ORDER: "
            f"{', '.join(sorted(off_centre))}. treating_centres ranges over CENTER_ORDER, so such a "
            "record is kept by no centre and dropped by restriction 1 without ever having been at a "
            "never-treating centre. The flow table would reconcile and the identifier would appear "
            "in restriction 1's removals.")

    off_arm = df.loc[~df[C.TREATMENT].isin((0, 1)), "case_id"]
    if len(off_arm):
        bad.append(
            f"C3  {C.TREATMENT}: {len(off_arm)} record(s) whose exposure is missing or outside "
            f"{{0, 1}}: {', '.join(sorted(off_arm))}. restriction 1's predicate counts bridging "
            "patients per centre, so a treated patient with no recorded exposure could take their "
            "centre's count to zero and remove the whole centre.")

    if bad:
        raise C.SchemaError(
            "\n  " + "\n  ".join(bad)
            + f"\n  {len(bad)} cohort assertion(s) failed against {len(df)} records. [§3]'s two "
            "restrictions are applied by design, before any weighting; none of these is a reason "
            "to relax that.")
```

- **C1 is checked rather than left to `retained`'s `KeyError`.** Stage 4 §5.1 declares that `KeyError`
  the correct failure for a caller asking who is retained before anyone is classified, and it stays
  correct — but by the time `build` reached it, restriction 1 would already have removed 22 patients
  and recorded an entry. C1 makes the failure happen before the log is touched.
- **C2 is checked here although `derive_cohort` checks it too**, and the duplication is the point:
  `derive_cohort`'s check fires at §4.5's call site, which is *after* both restrictions have been
  applied and all three `cohort` entries recorded. The audit would then describe a cohort that no
  caller ever received. One cheap membership test moves that failure to before the first removal.
- **C3 and C4 are Stage 4's E5 and E4, re-checked for Stage 4 §6.1's reason verbatim.** `build` is
  called on frames Stage 4 never saw — every acceptance test below builds its own, and a future
  caller subsetting between `classify` and `build` gets neither. An assertion that only runs when the
  caller happened to come through `classify` is a guarantee about a code path, not about a function.
  §3.1's first fact is what makes both silent rather than loud without them.
- **Every branch is unreachable on v7**, and that is recorded rather than treated as a reason to skip
  one (§18). Like Stage 3's `_assert_onset_flags` and Stage 4's E1-E5, they are written for the
  workbook that has not arrived yet.

**C4 is emitted before C3 although it is numbered after**, and the order is inert because the
messages are collected. The numbering follows Stage 4's E-numbering of the same two properties so a
reader moving between the files is not made to re-learn which is which; the emission order follows
the columns' order in the frame. If that trade reads badly to the implementer, renumber — nothing
depends on the numerals but §12's banners.

### 4.5 The order, and what does and does not turn on it

The restrictions are applied **1 then 2**, as [§3] lists them and as roadmap Stage 5 requires.

**The resulting cohort does not depend on the order.** Both are row masks over independent
predicates, so the retained set is an intersection, and intersection is commutative. Verified rather
than reasoned: applying 2 then 1 on v7 gives the same 93 records, identifier for identifier (§18).
§12.5 asserts it.

**The attribution does depend on the order, and that is what the flow table publishes.** On v7:

```
  1 then 2   (the declared order)      126 → 104 → 93      22 removed, then 11
  2 then 1                             126 → 107 → 93      19 removed, then 14
```

Both remove 33. A reader of the flow table is reading a decomposition, not a derivation, and §7.3's
row labels name the restriction each step applied so the decomposition cannot be misread as a claim
about which patients were "really" excluded for which reason. Eleven patients are both ineligible and
at a never-treating centre, and they are counted once, under restriction 1.

**Restriction 1's predicate is invariant to restriction 2, and that is a consequence of DECISION 1
rather than a coincidence.** Restriction 2 can never remove a treated patient: treated → `eligible`
by revealed fact, so no treated patient is `ineligible`, and Stage 4's E3 raises on the one frame
where that could fail. So `treating_centres` returns the same tuple before and after restriction 2 —
verified on v7 (§18) — which is why the two orders agree on the cohort and not merely on its size.
§12.5 asserts the invariance directly as well as the set equality, because the two would come apart
on exactly the frame E3 forbids, and a test that checked only the sets would not say which property
had broken.

Then, and only then:

```python
    df = derive.derive_cohort(df, audit)          # 32 → 33 columns; the median is the COHORT's
```

**This call is the reason `build` exists as a function rather than as three lines in a driver.**
`derive_cohort` freezes the [§13] core-volume median on the frame it is handed. A `derive_cohort`
called one line earlier produces a full, plausible `core_above_median` column describing a different
subgroup, and Stage 3 already records that no test would fail.

**The median alone does not catch every misplacement, and the reason is a measured fact about v7.**
There are three frames the call could land on, and on the workbook two of them share a median:

```
  frame                                  records   core_ml median
  as classified, before restriction 1      126        6.0 mL
  after restriction 1, before 2            104        5.0 mL      ← same as the cohort's
  the cohort                                93        5.0 mL
```

So a call moved to the top of `build` is caught by the median; a call moved **between the two
restrictions** — the likelier slip, and the position a reader tidying the function would choose — is
invisible to it on v7. On `cohort_frame()` the same three medians are 10.0 / 10.0 / 8.0, so the
hand-built frame does catch it, but only for as long as that frame keeps the property.

**What catches it on every frame is the entry's `records` cell.** `derive_cohort` writes
`("records", str(len(df)))` into its own `core_above_median` table (`derive.py`), and that cell reads
93 on v7 only if the call ran on the cohort — 126 or 104 otherwise. A row count cannot coincide the
way a median can. §12.7 therefore asserts **both**: the `records` cell equals `len(cohort)`, and the
`median (mL)` cell is the cohort's. Three cohort records sit exactly on 5.0 with ties placed below
(Stage 3 `_core_above_median`), which is why the median is worth asserting at all — it is the number
the subgroup is cut on — but the row count is what makes the guard total.

## 5. Postconditions

### 5.1 What is asserted about the result

Written out in full, as `_assert_cohort_inputs` and `build` are: this is the function the [§3]
guarantees are read off, and a reader reconstructing it from prose would be reconstructing the one
thing that must not be approximate.

```python
def _assert_cohort(df: pd.DataFrame, source: pd.DataFrame, removed: int) -> None:
    """The four postconditions, collected. `source` is the frame `build` was given, bound before the
    copy; `removed` is restriction 1's count plus restriction 2's. Raises SchemaError or returns."""
    bad: list[str] = []

    remaining = df.loc[df[C.ELIGIBILITY] == C.INELIGIBLE, "case_id"]
    if len(remaining):
        bad.append(
            f"P1  {C.ELIGIBILITY}: {len(remaining)} ineligible patient(s) remain in the cohort: "
            f"{', '.join(sorted(remaining))}. [§3] restriction 2 removes exactly these, and roadmap "
            "invariant 2 forbids one reaching [§7]. This is what fails if the mask is inverted.")

    empty = [(centre, C.TREATMENT_LABELS[code])
             for centre in treating_centres(source)
             for code in sorted(C.TREATMENT_LABELS)
             if not int(((df["center"] == centre) & (df[C.TREATMENT] == code)).sum())]
    if empty:
        bad.append(
            "P2  " + "; ".join(f"{centre} contributes no {arm} patient" for centre, arm in empty)
            + ". Restriction 1 guarantees a bridging arm at every retained centre, so this can only "
            "be restriction 2 having removed a centre's whole control arm: P(IVT = 0 | centre, X) = 0 "
            "there, the mirror image of the [§3] restriction-1 violation, and [§3] names only the "
            "bridging direction. Dropping the centre would extend [§3]; keeping it puts a "
            "structurally non-positive stratum in the propensity model, and `center` is a [§6] "
            "covariate. Amend [§3] with the PI; do not choose here.")

    if not len(df) or not treating_centres(df):
        bad.append(
            f"P3  the cohort is {len(df)} record(s) at {len(treating_centres(df))} centre(s). An "
            "empty cohort is a bug, not a result: every stage after this one divides by zero or fits "
            "over no data, and the first place it surfaces is a LinAlgError in Stage 6 with nothing "
            "pointing back here.")

    if removed + len(df) != len(source):
        bad.append(
            f"P4  the flow does not reconcile: {removed} removed + {len(df)} retained != "
            f"{len(source)} given. A record selected by NEITHER restriction's keep-mask leaves every "
            "table internally consistent [§3.1]; only a sum against the input frame notices.")

    if bad:
        raise C.SchemaError(
            "\n  " + "\n  ".join(bad)
            + f"\n  {len(bad)} cohort postcondition(s) failed on {len(df)} retained of "
            f"{len(source)} record(s).")
```

Four checks over the finished cohort, again collected:

- **P1 — no ineligible patient remains.** Roadmap invariant 2, asserted where it becomes true rather
  than where it is relied on. `(df[C.ELIGIBILITY] == C.INELIGIBLE).sum() == 0`. It is the one
  postcondition that restates the restriction that produced it, and it is cheap enough that the
  redundancy is worth having: it is what fails if restriction 2's mask is ever inverted.
- **P2 — every retained centre contributes both arms.** §5.2, which is most of this section.
- **P3 — the cohort is non-empty and retains at least one centre.** An empty cohort is a bug, not a
  result: every stage after this one would divide by zero or return a fit over no data, and the first
  place it would surface is a `LinAlgError` in Stage 6 with nothing pointing back here.

  **P3 is reached rather than crashed into, and that was measured rather than assumed** (§18). On a
  frame where restriction 1 keeps nothing, everything between the mask and P3 survives an empty
  frame: `derive_cohort` returns 0 rows and 33 columns with its median cell rendering the literal
  `missing` through `data._fmt`, `constant_covariates` reports every covariate constant, and
  `absence_by_column` guards its own `len(df)` division. So **do not add a defensive early return** —
  the empty frame walks to P3 and P3 is what names the problem. An early return would report the same
  failure without the flow table that shows which restriction emptied the cohort.
- **P4 — the flow reconciles.** `removed₁ + removed₂ + len(cohort) == len(source)`. This is the guard
  against the failure §3.1's first fact describes: a record silently unselected by a mask leaves
  every table internally consistent, and only a sum against the *input* frame notices. `_flow_table`
  computes its rows from the frames; P4 compares the total to the frame `build` was given.

### 5.2 The both-arms guard, which raises rather than choosing

Roadmap Stage 5's acceptance criterion is "every centre in the resulting cohort contains both arms".
It is a **runtime check that raises**, not a test-time assertion, and not a report. Its clause is the
`empty` comprehension in §5.1's body above; it is not repeated here, because two copies of a message
this specific is how the code and the specification come to say different things.

Four things about it:

- **It can only fire because of restriction 2.** Restriction 1 keeps exactly the centres with a
  bridging patient, so the treated arm is non-empty at every retained centre by construction. The
  message says so, because an implementer debugging it will otherwise look at the wrong restriction.
- **It is one workbook revision away, not hypothetical.** HUG contributes 22 controls of which 11 are
  ineligible (§18). A workbook in which the other 11 also carried a contraindication leaves HUG with
  30 bridging patients and no control arm, and HUG is 41 of the 93 cohort records.
- **Raising is the house posture and the alternatives are both worse.** Dropping the centre silently
  extends [§3] restriction 1 to a direction it does not name; keeping it hands Stage 6 a centre
  stratum in which the propensity score is 1 by construction, and `center` is a [§6] covariate that
  enters the design matrix. Stage 3 took this posture toward the both-onset-flags record and Stage 4
  toward the treated-and-flagged record: a derivation rule where one is available, a loud stop where
  none is.
- **It is what makes the fixture and the bare hand frame unusable as Stage 5 inputs**, which is §12.0
  and is the single largest consequence this spec has for the test suite.

## 6. The audit kind — the question Stage 4 deferred

### 6.1 `cohort` is added to `data.py`

Stage 4 §9.1 weighed a new kind and declined it *for that stage*, recording the reason it would
plausibly arrive here:

> A new `cohort` kind was weighed and declined **for this stage**, not in general. Stage 5 has to log
> rows *removed*, which no existing kind describes — `correction` is about values and `derivation` is
> about columns — so Stage 5 will plausibly insert one.

It does. `KINDS` gains `"cohort"` and `_HEADINGS` gains `"cohort": "Cohort construction"`. The
argument is exactly the one Stage 4 anticipated: the rule `KINDS` is organised on is *what the entry
is about*, and none of the seven existing kinds is about **rows removed**. Filing a restriction under
`derivation` would put "33 patients left the analysis" under a heading whose every other entry
describes a column, in a log whose §12.9-style inventory tests are the only thing standing between
this pipeline and an unrecorded exclusion.

Stage 3 §9.2's repair already made the insertion a one-line change: `test_data.py`'s
heading-adjacency test takes its neighbour from `KINDS` rather than naming it, so nothing about the
rendering breaks. What does break is three literal pins, and §10 is the ledger for them.

### 6.2 Its position: between `derivation` and `structural`

```
  provenance · contract · correction · observation · derivation · COHORT · structural · missingness
```

`KINDS` is documented as being "in pipeline order, which is also the order the document reads in".
Neither property can be had perfectly here, because Stage 3 has **two** entry points on opposite
sides of this stage: `derive`'s four entries run before the restrictions and `derive_cohort`'s two
run after them, and both are `derivation`. Whichever position `cohort` takes, one of those two groups
renders out of chronological order.

The position above is chosen because it is wrong about two entries rather than four, and because the
document then reads:

> what was read · what the contract did · what was corrected · what was observed · what was derived ·
> **who is in the analysis** · what is structurally absent · what is missing

which puts the cohort immediately before the section describing its denominators — and §7.4's entry
is a `missingness` table over the cohort, so the last two headings become "who is in, and what is
missing for them". Placing `cohort` before `derivation` would separate them and would be wrong about
`derive`'s four entries instead of `derive_cohort`'s two.

### 6.3 `_MUST_NAME_CASES` is deliberately unchanged

This looks like the natural amendment and it is the wrong one. `_MUST_NAME_CASES` is keyed by
**kind**, and it fires whenever `n > 0 and not case_ids`. Stage 5's kind carries two different sorts
of entry:

```
  restrict_centres        n = 22   removes patients   MUST name them
  restrict_eligibility    n = 11   removes patients   MUST name them
  cohort_flow             n = 93   removes nobody     has nobody to name
```

A kind-keyed rule cannot express that, and forcing it to would mean either giving `cohort_flow` an
`n` that misrepresents what it did, or printing the cohort's 93 identifiers under a summary table.

So the rule is enforced where it can be stated exactly — in `cohort.py`, by the single helper both
removals go through:

```python
def _record_removal(audit: Audit, step: str, removed: pd.DataFrame, detail: str,
                    table: tuple[tuple[str, ...], ...]) -> None:
    """The one way a [§3] restriction reaches the log. Names every record it removed, or raises.

    Stronger than data.py's kind-keyed rule, which asks only that *some* case be named: this asserts
    the entry names EXACTLY as many patients as it says it removed. A restriction that reported 11
    and named 3 would satisfy _MUST_NAME_CASES and would be a cohort nobody could reconstruct.
    """
    case_ids = tuple(sorted(removed["case_id"]))
    if len(case_ids) != len(removed):
        raise C.SchemaError(
            f"{step} removed {len(removed)} record(s) and names {len(case_ids)}. A [§3] restriction "
            "that cannot name what it removed makes the cohort unreconstructable from the log.")
    audit.record("cohort", step, len(removed), detail, case_ids=case_ids, table=table)
```

The `!=` can only fire on a duplicated or missing `case_id`, which Stage 2's A2 forbids — so this is
a guard on a property held one stage upstream, written here because this is the stage whose output
depends on it. §12.9 fires it against a hand-built frame.

`test_data.py`'s two kind-parametrised tests pick the new kind up automatically: the `else` branch of
`test_an_entry_touching_records_without_naming_them` asserts a `cohort` entry with `n = 126` and no
identifiers is legal, which is precisely the `cohort_flow` case, and no edit is needed there.

## 7. The four audit entries

### 7.1 `restrict_centres` — kind `cohort`

`n` is the number of records removed; `case_ids` names every one of them.

`detail`:

> `[§3] restriction 1: {k} of {len(CENTER_ORDER)} declared centre(s) contributed at least one
> {TREATMENT_LABELS[1]} patient and are retained; {d} contributed none and are dropped with their
> {n} patient(s). P(IVT = 1 | centre, X) = 0 there structurally, not by chance — no weighting
> recovers a contrast that was never available, and a penalised score would report shrinkage as
> though it were treatment availability.`

table — **one row per declared centre, rendered whether or not the data fills it**, for Stage 4
§7.1's reason: the empty cells are the information, and this is the table in which the [§3]
restriction-1 finding is the all-zero `bridging` cell.

```
  | centre | EVT alone | bridging | n   | status         |
  | HUG    | 22        | 30       | 52  | retained       |
  | CHUV   | 14        | 7        | 21  | retained       |
  | Lugano | 29        | 2        | 31  | retained       |
  | USZ    | 22        | 0        | 22  | dropped        |   ← [§3] restriction 1
  | all    | 87        | 39       | 126 | 22 removed     |
```

The arm columns come from `TREATMENT_LABELS` in `sorted()` key order — control before treated, in
every interpreter, and identical to Stage 4 §7.1's crosstab so the two tables can be read against
each other. The `status` cell is `retained` / `dropped` and is derived from membership of
`treating_centres(df)`, never recomputed, so the table and the restriction cannot disagree.

**The `all` row is not a reconciliation.** It is computed from the frame, so a record at an undeclared
centre would leave it correct while vanishing from the four centre rows. That failure is C4's, and C4
raises before this table is built. Stage 4 §7.1 needed both a precondition and an in-table
reconciliation because its table was the only thing rendering those cells; here the same job is done
by C4 plus P4, and P4 is the stronger of the two because it sums against the input frame rather than
against a row of the table.

### 7.2 `restrict_eligibility` — kind `cohort`

`n` is the number of records removed; `case_ids` names them.

`detail`:

> `[§3] restriction 2: {n} patient(s) carry ivt_contraindicated = 1 and are removed — they were
> never candidates for bridging, and retaining them makes "no IVT" a marker of contraindications
> and their prognosis. Eligibility is classified from the flag alone [DECISION 1a]; the {r}
> patient(s) with no documented contraindication reason are eligible, and the free text is not read.
> The flag was recorded on every record, including at the {c} centre(s) that never collected a
> reason — that it means the same thing at all of them is an assumption, and [§3] requires it stated
> as a limitation.`

table — one row per **retained** centre, over the frame restriction 1 produced:

```
  | centre | eligible | ineligible | n   | removed |
  | HUG    | 41       | 11         | 52  | 11      |
  | CHUV   | 21       | 0          | 21  | 0       |
  | Lugano | 31       | 0          | 31  | 0       |
  | all    | 93       | 11         | 104 | 11      |
```

The `all / eligible` cell is **93**, which is the cohort size — an identity rather than a
coincidence, since restriction 2 removes exactly the ineligible. It is worth noticing in the log
because it is the one place the two restrictions' arithmetic meets before §7.3 states it.

Ranging over `treating_centres`' result and not over `CENTER_ORDER` is deliberate: USZ is gone by the
time this restriction runs, and rendering an all-zero USZ row would suggest its patients were
considered for eligibility and found eligible. The tuple is declared, so the order is still not the
data's.

The class columns come from `ELIGIBILITY_ORDER`. The `removed` column is the `ineligible` column by
construction — printed anyway, because it is the column a reader checks the `n` against, and because
if the two ever differ the restriction has stopped being "drop the ineligible".

### 7.3 `cohort_flow` — kind `cohort`

`n` is `len(cohort)`. No `case_ids`: it removes nobody (§6.3).

`detail`:

> `[§3] cohort flow. The primary cohort is {n} of {N} record(s) at {k} centre(s): {t} bridging and
> {c} {TREATMENT_LABELS[0]}. Restrictions are applied in the declared order and each record is
> counted once, under the first restriction that removed it — the retained set does not depend on
> that order but the attribution does [Stage 5 §4.5]. {u} retained patient(s) have no documented
> contraindication reason and are eligible from the flag alone [DECISION 1a].`

table:

```
  | step                                  | centres | records | EVT alone | bridging | eligible | ineligible |
  | as classified                         | 4       | 126     | 87        | 39       | 107      | 19         |
  | after restriction 1                   | 3       | 104     | 65        | 39       | 93       | 11         |
  | after restriction 2                   | 3       | 93      | 54        | 39       | 93       | 0          |
  | of which EVT alone, no reason on file  | 2       | 43      | 43        | 0        | 43       | 0          |
```

- **The last row survives DECISION 1a, and it is the reason roadmap Stage 5's requirement survives
  too.** It no longer reports a *class* — there is no third class — but it reports the same 43
  patients, now as a **provenance** count: retained patients whose eligibility rests on the flag with
  no documented reason behind it. That is what a reader needs to size the assumption DECISION 1a
  makes, and it is the group a future [§13] sensitivity analysis would be sized against. Its
  `centres` cell reads **2**, not 3, and that difference is the finding: the undocumented group is
  centre-driven, not patient-driven — it is precisely the control arms of CHUV and Lugano. The 43 are
  80% of the cohort's 54 controls.
- **It is restricted to the `TREATMENT_LABELS[0]` arm, and that restriction is not cosmetic.** Over
  all arms the same count is **82**, because no treated patient carries a reason either (§18) — and a
  row reading 82 would be measuring nothing, since a treated patient's eligibility is established by
  the fact that they received IVT. The assumption DECISION 1a introduces is about **controls** with no
  documented reason, and this row must count exactly those. The row label carries the arm name from
  `TREATMENT_LABELS`, so it cannot say one thing and count another.
- **It is computed from `contraindication_reason.isna()`, read for *counting* and never for
  classifying.** The distinction is the whole of what survives DECISION 1a, and §12.6's scan is what
  keeps it: the column may reach a table or a `detail` string, and must never reach a mask that
  decides a class. Stage 5 is not the only reader — `eligibility.classify` counts the same
  undocumented controls for its own `detail` — so the rule is *how* the column is read, not *where*.
- **`ineligible` reads 0 in the cohort row**, which is P1 rendered rather than asserted. Both are
  kept: the column because a reader needs to see it fall to zero, P1 because a table is not a check.
- **Every column is computed from the frame at that step**, never by subtracting the row above.
  Subtraction would make the table reconcile with itself by construction and stop being evidence.

### 7.4 `absence_by_cohort_column` — kind `missingness`

One call to the existing `data.absence_by_column`, over the finished cohort, after `derive_cohort` so
that `core_above_median` exists:

```python
_ABSENCE_DETAIL: Final[str] = (
    "one row per analysis and derived column, over the PRIMARY COHORT rather than over the frame "
    "as read [§11]. Every estimate carries its own denominator; this is where the denominators "
    "of the [§7] analyses are read off. Nothing here is imputed.")

# in build(), on the returned frame and not on `after_2` — core_above_median must exist
    absence_by_column(out, audit, [*sorted(C.ANALYSIS_NAMES), *sorted(C.DERIVED_NAMES)],
                      "absence_by_cohort_column", _ABSENCE_DETAIL)
```

The detail is a module constant rather than a literal at the call site, so §8's `build` stays twelve
readable lines; the other three entries' details are computed per call and stay in their helpers.

**Why this entry exists at all.** The log currently reports missingness twice, both times over 126
records — Stage 2's `absence_by_column` and Stage 3's `absence_by_derived_column`. Every denominator
in [§7] and [§8] is a property of the 93-record cohort, and nothing in the log states one. [§11]
requires each estimate to carry its own denominator and [§16] requires them reported; a reader of the
audit log presently has to reconstruct them.

**Why it is Stage 5's and not Stage 14's.** Stage 14 is the reporting layer and writes tables and
figures; the audit log is the pipeline's own record of what it did to the data, and the cohort is
what this stage did. Deferring it would mean the one stage that creates the population is also the
one that says nothing about it.

**This is a decision, not a derivation, and it is the one item in this spec a reader may reasonably
want struck.** It adds a fourth entry and a `missingness` row per column to a log that already has
two such tables. If the PI prefers the denominators reported once, at Stage 14, striking it removes
§7.4, §12.10 and one line of §9 and changes nothing else — no other section depends on it. The
argument for keeping it is that the two existing tables are over a population no estimate uses.

The two column lists are concatenated rather than unioned, so the row order is declared twice over
and the log is invariant to `frozenset` iteration — `ANALYSIS_NAMES` is a `frozenset` and Stage 2
already sorts it for this reason. `eligibility` appears in neither list, deliberately: Stage 4 §8
keeps it out of `DERIVED_NAMES` so it can never become a covariate, and it is total by construction
so a missingness row over it could never say anything.

## 8. `build`, written out

The public function is short and every line is load-bearing, so it is given here rather than left to
be reconstructed from §0's diagram.

```python
def build(df: pd.DataFrame, audit: Audit) -> pd.DataFrame:
    """The [§3] primary cohort: two restrictions in the declared order, then Stage 3's cohort-wise
    derivations. Returns a new frame of 33 columns; appends three `cohort` entries, one
    `missingness` entry and — through derive_cohort — two `derivation` entries.

    **The only caller of derive.derive_cohort in the repository**, which is what makes the [§13]
    median a property of the cohort rather than of whichever frame a driver happened to pass. The
    median is 5.0 mL here and 6.0 mL on the frame this receives [§4.5].

    Raises SchemaError on C1-C4 before any row is removed, and on P1-P4 after the cohort is built.
    """
    _assert_cohort_inputs(df)
    source = df
    df = df.copy()

    keep = treating_centres(df)
    dropped = df[~df["center"].isin(keep)]
    _record_removal(audit, "restrict_centres", dropped, _centres_detail(df, keep, dropped),
                    _centres_table(df, keep, dropped))
    after_1 = df[df["center"].isin(keep)]

    ineligible = after_1[~eligibility.retained(after_1)]
    _record_removal(audit, "restrict_eligibility", ineligible,
                    _eligibility_detail(after_1, ineligible), _eligibility_table(after_1, keep))
    after_2 = after_1[eligibility.retained(after_1)]

    audit.record("cohort", "cohort_flow", len(after_2),
                 _flow_detail(source, after_2), table=_flow_table(source, after_1, after_2))

    out = derive.derive_cohort(after_2, audit)
    absence_by_column(out, audit, [*sorted(C.ANALYSIS_NAMES), *sorted(C.DERIVED_NAMES)],
                      "absence_by_cohort_column", _ABSENCE_DETAIL)
    _assert_cohort(out, source, len(dropped) + len(ineligible))
    return out
```

Seven things an implementer could reasonably get wrong:

- **The preconditions run before the copy**, so a frame that fails is never copied and C1-C4 are
  reported about the caller's own frame. `_assert_cohort_inputs` only reads.
- **`source` is bound to the caller's frame before the copy**, and P4 compares against it. Binding it
  after would compare the copy to itself and P4 would be a tautology.
- **`build` copies once, at the top, and the two restrictions are plain masks on the copy.** A
  `.copy()` inside each restriction would be two more full frames for no gain; `df[mask]` returns a
  new object either way, and nothing here writes into `after_1` or `after_2`.
- **`derive_cohort` is called on `after_2`, after both restrictions** — §4.5. Moving it above either
  one is the failure this stage exists to prevent and the only one that produces a full, plausible
  column.
- **`eligibility.retained` is called twice rather than stored**, once for the removed frame and once
  for the kept one, and the two must not be replaced by a stored mask and its negation *unless* the
  negation is written as `~mask`. `retained` is total (Stage 4 §5.1), so `~` is safe here — but the
  symmetry with `~df["center"].isin(keep)` above, where `isin` is also total, is what makes both
  readable. Either form is correct; what is not correct is `!= C.INELIGIBLE` (§4.3, §12.6).
- **`_record_removal` is used for both restrictions and `audit.record` directly for the flow**, which
  is §6.3's distinction made visible in the code: the two entries that remove patients go through the
  helper that forces them to be named, and the one that does not, does not.
- **P1-P4 run on `out`, after `derive_cohort`**, not on `after_2`. `derive_cohort` adds a column and
  drops no row, so the checks would pass either way — but asserting on the returned object is what
  makes them assertions about what the caller receives.

The `_*_detail` and `_*_table` helpers are §7's four entries; they are named here rather than inlined
because `build` is the function a reader traces the [§3] restrictions through and three f-strings
would bury it.

## 9. Data flow into Stages 6-14

```
  build(df, audit) returns
   ├─ the 32 columns it was given, unchanged — no value edited
   ├─ core_above_median   Int64, frozen on THIS cohort [§13, Stage 3 §6.4]
   ├─ 33 columns, 93 rows on v7
   └─ the index of the unrestricted frame, preserved and unique [§0.2]

  audit — the same object load() created. Entry positions are asserted against a captured
  index, never against a remembered total [Stage 4 §12.9]:

    load()            6 on the fixture, 7 on the workbook
    derive()          4
    classify()        1
    build()           3 cohort  →  restrict_centres, restrict_eligibility, cohort_flow
                      2 derivation (derive_cohort's, recorded from inside build)
                      1 missingness (absence_by_cohort_column)
```

What each later stage may rely on, and what each owes:

- **Stage 6 [§7]** — fits the propensity model on this frame and nothing else. Every retained centre
  has both arms (P2), so no centre stratum is structurally non-positive. `center` is a `string`, so
  the design-matrix builder makes its own `Categorical` from `FACTOR_LEVELS` at the point of use and
  a dropped centre leaves no empty dummy column.
- **Stage 7 [§9]** — the within-centre overlap table ranges over `treating_centres(cohort)`, which on
  this frame is every centre present. It reads §7.4's table for denominators.
- **Stage 10 [§10]** — resamples **this** frame, stratified by centre, and never calls `build` or
  `derive_cohort` inside a replicate. `core_above_median` is resampled along with the patient. Stage
  3's `derive_cohort` raises if a replicate tries, and C2 raises one call earlier.
- **Stage 12 [§14a]** — does **not** take this frame. It takes the unrestricted classified frame,
  draws its population with `eligibility.retained`, and calls `treating_centres` only to compute its
  complement for the support check (§3).
- **Stage 13 [§14b]** — likewise takes the unrestricted frame, and reads the `eligibility` column
  itself rather than the predicate. Under DECISION 1a that column is **binary**, which is what closed
  [§14b]'s open question about assigning treatment for a third value; there is no three-level column
  anywhere in the pipeline and a later stage must not look for one.
- **Stage 14 [§16]** — the cohort-flow table of §7.3 is the source for the manuscript's flow
  diagram, and §7.4's table for the per-estimate denominators [§11].

## 10. What Stage 5 amends in Stages 1-4

The full ledger, so that no amendment is discovered during implementation.

**Nine files is the cost of one audit kind, not scope creep, and the distinction is worth stating
because the file count is what a reviewer reacts to.** Two files are new. Of the seven amended, five
are two-or-three-line edits that `KINDS` gaining a member *forces*: three literal pins across
`test_data.py` and `test_eligibility.py`, plus two test renames. They cannot be deferred — a pin that
disagrees with `KINDS` is a red suite — and they cannot be avoided, because pinning against literals
is exactly what those tests are for. The two amendments that are genuinely Stage 5's own judgement
are `test_derive.py`'s tail slice (which `TODOS.md` already assigned to this commit) and `data.py`'s
prose. Checked at review rather than assumed: `test_data.py`'s two kind-parametrised tests range over
`data.KINDS` and need **no** edit, and neither does its heading-adjacency test, which takes its
neighbour from `KINDS` rather than naming it. T1's verification step is what turns that claim from a
belief into a check.

| File | Amendment | Why |
|---|---|---|
| `data.py` | `KINDS` gains `"cohort"` between `"derivation"` and `"structural"`; `_HEADINGS` gains `"cohort": "Cohort construction"` in the same position | §6.1, §6.2. The first `data.py` amendment since Stage 2, and the one Stage 4 §9.1 predicted |
| `data.py` | `_MUST_NAME_CASES` — **no change to the frozenset**, and the comment above it gains one sentence saying `cohort` is deliberately excluded and where the rule lives instead. Its "the other **five** kinds legitimately name none" becomes **six** | §6.3. Left undocumented, the omission reads as one; and a count that is off by one in the same comment is how a reader learns not to trust it |
| `data.py` | **A9's error message** loses its DECISION 1 clause. It currently reads "a blank reads as 'a reason was recorded' and moves a patient from indeterminate to eligible [DECISION 1, §3]" — a class that no longer exists, in **runtime text shown to the data owner**. A9 now guards data quality only, and must say so | §21 R7. DECISION 1a's sweep corrected the module that depended on A9 and left A9's own message describing the dependency it removed |
| `data.py` | the `_missingness` comment block above `absence_by_column` is rewritten for two classes. It presently says the reason column's absence "is the one bit the column contributes, the thing that separates eligible from indeterminate" and "makes the indeterminate group most of the control arm" | §21 R7. Both sentences describe a withdrawn rule. The absence is still centre-driven and still worth the per-centre columns — that argument survives; only the class does not |
| `test_data.py` | `test_the_headings_are_the_declared_seven_in_pipeline_order` — both literals gain the new kind and heading, and the test is **renamed** (`_seven_` → `_eight_`) | it pins `KINDS` and `list(_HEADINGS.values())` against literals, which is the point of it |
| `test_data.py` | the comment inside that test gains the `cohort` position argument of §6.2 | it currently explains only why `derivation` sits where it does |
| `test_eligibility.py` | `test_data_kinds_is_unchanged_by_this_stage` — the literal gains `"cohort"`; the test is **renamed** to say what it now pins, and its comment loses the prediction and gains the pointer to §6 | its comment says "Stage 5 […] will plausibly insert one". It has. A test whose name says *unchanged* and whose body pins a changed tuple is worse than either |
| `test_derive.py` | `test_the_six_entries_appear_with_the_declared_kinds_and_steps_in_order` — replace the tail slice `[-len(_STAGE_3_INVENTORY):]` with an index captured before the calls, as `test_derives_four_entries_precede_derive_cohorts_two` two lines below already does | **`TODOS.md`'s deferred item, whose "do it with Stage 5" this commit satisfies.** No longer hypothetical: `build` puts three `cohort` entries between `derive`'s four and `derive_cohort`'s two, so any test exercising the real order makes `[-6:]` span them |
| `implementation_roadmap.md` | Stage 5 gains its `**Spec:**` line; its **Accept when** gains the both-arms guard as a runtime raise and the frozen-median ordering | §19. Precedent: Stages 3 and 4 landed their spec lines with the document, not with the code |
| `TODOS.md` | the tail-slice item is struck | it is done by this commit; an open item that has been done is a trap in the other direction |

**Four things that look like they need amending and do not.**

- **`config.py`.** Stage 5 declares no constant (§1). `CENTER_ORDER`, `TREATMENT`,
  `TREATMENT_LABELS`, `COHORT_DEPENDENT_SUBGROUPS` and `ELIGIBILITY_RETAINED` all exist and all are
  read rather than written.
- **`derive.py`.** Its docstrings already carry the post-Stage-4 column counts — `31 → Stage 4 →
  Stage 5 restrictions`, and `derive_cohort` "returns a new frame of 33 columns" — because Stage 4's
  T6 landed them. Read, not assumed (§18). Its sentence "**Stage 5 is `derive_cohort`'s only
  caller**" becomes true rather than aspirational with this stage and needs no edit.
- **`eligibility.py`.** `retained` is called, not changed. Stage 4 §13's "retained is not enforced"
  is answered by §12.6's scan over `cohort.py`, not by an edit there.
- **`test_data.py`'s two kind-parametrised tests.** They range over `data.KINDS`, so `cohort` is
  covered the moment it is declared, and its `else` branch is the `cohort_flow` case (§6.3). Read,
  not assumed (§18).

## 11. Handover to Stage 6

```
  cohort  =  build(classify(derive(load()[0], audit), audit), audit)

    93 rows, 33 columns, v7
    every centre both-armed                      P2
    no ineligible patient                        P1, invariant 2
    no never-IVT centre                          restriction 1, invariant 1
    core_above_median frozen at 5.0 mL           [§13], Stage 3 §6.4
    index = the unrestricted frame's labels      §0.2
```

Stage 6 receives this and fits. It does not re-restrict, does not re-classify, and does not call
`derive_cohort`.

## 12. Acceptance criteria

All tests live in `test_cohort.py`, with section banners matching these numbers, as `test_derive.py`
and `test_eligibility.py` do. Tests needing the private workbook reuse the `DATA_GATED` idiom;
everything else runs against hand-built frames. `test_cohort.py` is not exempt from the raw-name scan
and builds its frames with analysis names.

### 12.0 The frame this stage is tested on, and why it is new

**Neither existing hand-built input can reach the end of this stage.** Measured (§18), not
anticipated:

```
                        records   restriction 1     restriction 2   per retained centre (t, c)
  fixture_schema.xlsx      2      keeps HUG only          1         HUG (1, 0)      → P2 RAISES
  hand_frame()             6      keeps HUG, CHUV         3         HUG (2, 0)      → P2 RAISES
                                                                    CHUV (1, 0)
```

Both collapse to a cohort with no control arm, because both were built for stages that never removed
a row: the fixture is a *schema* fixture — two records, enough to exercise the read contract — and
`hand_frame()`'s six records carry one control at each of the two centres that contribute no treated
patient. Neither is wrong; neither is a cohort.

Three routes were considered:

| Considered | Declined because |
|---|---|
| Extend `hand_frame()` with controls at HUG and CHUV | Stage 4 §12 states it "needs **no** extending" and every Stage 4 test is written against its exact six-record classification and crosstab totals. DECISION 1a already rewrites those numbers once (§20 item 7); rewriting them a second time, for a different reason, in a file this stage has no business touching, would make the two edits impossible to review apart |
| Add rows to `tests/fixture_schema.xlsx` | Stage 2, 3 and 4 tests assert its two records by count and by classification (`test_the_fixture_reaches_every_entry`, Stage 4 §12.9's driver). It is a committed binary with no generator script (§18), so the edit is both invasive and unreviewable in a diff |
| A second committed workbook, `tests/fixture_cohort.xlsx` | It would make the whole pipeline runnable from a file with no `data/`, which Stages 6-13 will want. But it is a new binary artifact and a new `data.Source` to justify a frame that ten lines of Python produce. Recorded in §13 as the thing to build if a later stage needs a *file* rather than a *frame* |

So `test_cohort.py` declares its own, **built on `hand_frame()` rather than replacing it**:

```python
_COHORT_RECORDS = [
    # a control at a treating centre with a DOCUMENTED reason — without this one, no centre has
    # both arms
    dict(case_id="COHORT-1", center=CENTER_CODES[0], ivt=0, ivt_contraindicated=0,
         contraindication_reason="Clinician decision", core_ml=8.0, tmax6_ml=50.0,
         penumbra_ml=42.0, onset_to_ivt_min=None, onset_to_groin_min=330.0),
    # the only ineligible record anywhere in the file: restriction 2 removes exactly this one
    dict(case_id="COHORT-2", center=CENTER_CODES[0], ivt=0, ivt_contraindicated=1,
         contraindication_reason="Anticoagulation", core_ml=12.0, tmax6_ml=70.0,
         penumbra_ml=58.0, onset_to_ivt_min=None, onset_to_groin_min=280.0),
    # a control with NO documented reason that SURVIVES both restrictions. Under DECISION 1a this
    # is the record that proves an undocumented patient is eligible rather than a third class, and
    # it is the one §7.3's last row counts. It is deliberately at a DIFFERENT centre from COHORT-1,
    # so the flow table's `of which` row reads 1 centre against 2 in the row above it
    dict(case_id="COHORT-3", center=CENTER_CODES[1], ivt=0, ivt_contraindicated=0,
         contraindication_reason=None, core_ml=3.0, tmax6_ml=40.0,
         penumbra_ml=37.0, onset_to_ivt_min=None, onset_to_groin_min=420.0),
]


def cohort_frame(**overrides: object) -> pd.DataFrame:
    """`hand_frame()` plus three records, so that a cohort with both arms exists. §12.0.

    Built as one DataFrame from records rather than by `pd.concat`, which emits a FutureWarning on
    pandas 2.3.3 when an operand carries an all-NA column — `onset_to_ivt_min` is None on all three
    records below — and whose dtype resolution for that case is documented as changing. Verified
    warning-free under `simplefilter("error")` (§18).
    """
    base = hand_frame()
    template = base[base["case_id"] == "HAND-4"].iloc[0].to_dict()
    records = [*base.to_dict("records"), *({**template, **r} for r in _COHORT_RECORDS)]
    df = pd.DataFrame(records)
    for column, dtype in ANALYSIS_DTYPES.items():
        df[column] = df[column].astype(dtype)
    for column, value in overrides.items():
        df[column] = value
    return df[sorted(config.ANALYSIS_NAMES)]
```

`hand_frame`, `CENTER_CODES`, `ANALYSIS_DTYPES`, `corrupt`, `hand_source` and `run` are **imported
from `test_data.py`, never re-declared**, for Stage 3 §12's reason; `config` and `pandas as pd` are
imported directly, as every other test module does — the body above uses all three. `HAND-4` is the
template because
it is the one record whose volumes are internally consistent and which carries no correction —
`HAND-1` carries the copy-paste penumbra that Stage 2 rewrites, and `HAND-2` the 999 groin time that
raises an observation.

Its verified classification and flow (§18):

```
  HAND-1   HUG     bridging  flag 0  no reason   → eligible     removed by nothing
  HAND-2   CHUV    bridging  flag 0  no reason   → eligible     removed by nothing
  HAND-3   Lugano  control   flag 0  reason      → eligible     removed by restriction 1
  HAND-4   USZ     control   flag 0  reason      → eligible     removed by restriction 1
  HAND-5   HUG     bridging  flag 0  no reason   → eligible     removed by nothing
  HAND-6   Lugano  control   flag 0  no reason   → eligible     removed by restriction 1
  COHORT-1 HUG     control   flag 0  reason      → eligible     RETAINED
  COHORT-2 HUG     control   flag 1  reason      → ineligible   removed by restriction 2
  COHORT-3 CHUV    control   flag 0  no reason   → eligible     RETAINED, undocumented

  9 records → restriction 1 keeps HUG and CHUV, removing 3 → restriction 2 removes 1 → 5
  cohort: HUG (2 bridging, 1 control), CHUV (1 bridging, 1 control)
  8 eligible, 1 ineligible;  in the cohort, 5 eligible and 0 ineligible
  §7.3's `of which EVT alone, no reason on file` row, in the cohort: 1 record at 1 centre
      — COHORT-3. Over ALL arms it would read 4 at 2 centres, which is why the row is
        arm-restricted; over all arms and before restriction 1 it would read 5 at 3.
  core_ml median: 10.0 over the 9 records, 8.0 over the cohort — so §12.7 has two distinct numbers
```

Four properties of that design are deliberate and must survive an edit to it: restriction 1 removes
records **at more than one centre**; restriction 2 removes exactly one record, at a centre that keeps
a control (so P2 does not fire); the flow table's `of which` row reads **1 record at 1 centre**
against **2 centres** in the row above it, and reads a different number again if the arm restriction
is dropped — which is how §12.8 tells that row is computed, arm-restricted, and not copied; and the
two medians differ, which is the only way §12.7 can tell where `derive_cohort` was called.

**Note what stopped being testable here, and where it went.** Under DECISION 1 this frame carried two
`indeterminate` records and the acceptance tests could assert a three-way partition. There is no
third class now, so the property that replaces it is narrower and more direct: `COHORT-3` has no
documented reason and is `eligible` — §12.2.

**And the two existing frames are kept, as positive tests of P2** (§12.4). Their inability to reach
Stage 5 is the guard working, not an obstacle to routing around.

### 12.0.1 Bare frame or normalised frame — the distinction every test below turns on

`hand_frame()` and `cohort_frame()` return the frame **as the reader would deliver it**, carrying the
raw centre codes `("1", "Lausanne", "Lugano", "USZ")` from `CENTER_RECODE`'s keys. `test_data.run()`
is what maps them to `CENTER_ORDER`'s labels. Both forms are used below and they exercise *different*
checks, so every acceptance test must say which it means:

```
  cohort_frame()                       raw codes    → C4 fires: no centre is in CENTER_ORDER
  run(cohort_frame(), hand_source(…))  labels       → the frame the restrictions are written for
```

Unless a section says otherwise, **every frame below is the normalised one** — `run(...)` followed by
`derive` and `classify`. The two deliberate exceptions are §12.3's C4 case, which uses the bare frame
precisely because it is one, and §12.11's, which stops after `classify`. A test that means to reach P2
and passes a bare frame gets a C4 raise instead and proves nothing: the failure looks like a passing
`pytest.raises(SchemaError)` and the message is never read.

### 12.1 The two restrictions

On `cohort_frame()` through `test_data.run()`, `derive`, `classify` and `build`: the cohort is
exactly `{HAND-1, HAND-2, HAND-5, COHORT-1, COHORT-3}`, asserted **by identifier, never by row
position or count alone**. `treating_centres` returns `("HUG", "CHUV")` — in `CENTER_ORDER` order,
not the frame's. A count-only assertion would pass on a cohort that removed `COHORT-1` and kept
`COHORT-2`, which is the inversion §12.6 exists to prevent and which would leave the count at 5.

### 12.2 An undocumented patient is eligible, and the reason column is not read [DECISION 1a]

`COHORT-3` has `ivt_contraindicated = 0` and no `contraindication_reason`, is classified `eligible`,
and is **in the cohort**. `COHORT-1` differs from it only in carrying a reason and is also in the
cohort. Both halves, because either alone passes against a broken restriction — the first against one
that retains everything, the second against one that still reads the reason column and drops the
blank.

The pointed version of the second half, and the one that would have failed under DECISION 1:
`cohort_frame(contraindication_reason=None)` — the reason column blanked on **every** record — gives
the **same cohort**, identifier for identifier. Under DECISION 1 it would have moved three records
into `indeterminate`; under DECISION 1a the column contributes nothing to any class and the cohort
cannot move. It is one line and it is the whole of the amendment.

### 12.3 C1-C4 fire, naming their cases

- **C1** — a frame with no `eligibility` column raises `SchemaError` naming the column, **and the
  audit is untouched**: `len(audit.entries)` is unchanged after the raise. That second half is what
  C1 is for; without it the test passes against a `build` that raises from `retained`'s `KeyError`
  three entries later.
- **C2** — a frame already carrying `core_above_median` raises, and again the audit is untouched.
  Reached by calling `build` twice on the same input, which is also the shape a bootstrap replicate
  would take (§9).
- **C3** — an exposure of `2`, and a **missing** exposure, each raise. Both must be corrupted *after*
  Stage 2, because A5 and A4 raise on them first; the construction note below is the same one Stage 4
  §12 records.
- **C4** — a centre of `"Bern"` raises, and so does a **missing** centre: `isin` returns `False` for
  `<NA>`, so `~isin` catches it, and the test asserts that cell of §3.1's second fact rather than
  assuming it. Both must be injected **after `classify`**, and the reason is stronger than C3's: two
  earlier checks fire on the same cell first, in two different modules. Stage 2's A3 rejects any
  centre outside `CENTER_RECODE.values()`, and Stage 4's E4 rejects any centre outside `CENTER_ORDER`
  — so a frame corrupted before either one never reaches C4, and a test that corrupts early passes
  while asserting nothing about this stage.
- A frame tripping several at once reports all of them in one message; assert the identifiers of two
  different checks in a single `SchemaError`. The **bare** `cohort_frame()` is the cheapest such
  frame: it has no `eligibility` column and no centre in `CENTER_ORDER`, so C1 and C4 fire together.

```
  C1  eligibility column absent      unreachable through classify(); drop the column after it
  C2  core_above_median present      reachable: call build twice
  C3  ivt = 2                        A5 raises first  — corrupt AFTER Stage 2
  C3  ivt missing                    A4 raises first  — corrupt AFTER Stage 2
  C4  centre = "Bern"                A3 and E4 raise first — corrupt AFTER classify
  C4  centre missing                 A3 and E4 raise first — corrupt AFTER classify
  C1 + C4 together                   the BARE cohort_frame(): raw codes, never classified
```

### 12.4 P1-P4 fire, and the two existing frames are what fire P2

- **P2 on `tests/fixture_schema.xlsx`** — `build` raises `SchemaError` naming `HUG` and the control
  arm. Data-gate-free: the fixture is on every checkout.
- **P2 on `run(hand_frame(), …)` through `derive` and `classify`** — raises naming **both** `HUG` and
  `CHUV`, which is what proves the message enumerates rather than reporting the first. Verified: the
  clause returns `[("HUG", "EVT alone"), ("CHUV", "EVT alone")]`. It must be the **normalised** frame
  (§12.0.1); the bare one raises C4 long before P2 and the test would pass on the wrong exception.
- **P2 on `cohort_frame(...)` with `COHORT-1` flagged contraindicated** — the frame in which HUG's
  only remaining control is removed by restriction 2. This is the one that reproduces the failure
  §5.2 is written for, on a frame that is otherwise a valid cohort, and its message must name the
  arm as well as the centre.
- **P1** — asserted on every frame here, and driven to fail by a locally reimplemented two-line
  restriction that keeps `ineligible` (never by monkeypatching `cohort.py`, per Stage 3 §12.3).
- **P3** — a frame with no treated patient at any centre: restriction 1 keeps nothing, and the raise
  says so rather than returning an empty frame.
- **P4** — driven by a locally reimplemented `build` whose restriction-1 mask is `df["center"].isin(keep)`
  on a frame carrying a missing centre with C4 removed: the record is selected by neither mask,
  every table reconciles, and P4 is the only thing that notices. This is §3.1's first fact and it is
  the reason P4 sums against the input frame. The missing centre is injected **after `classify`**, for
  §12.3's reason — A3 and E4 both raise on it earlier.

### 12.5 The order does not change the cohort, and restriction 1 is invariant to restriction 2

Two assertions, both on `cohort_frame()` and, data-gated, on the workbook:

1. Applying the restrictions in the reverse order gives the **same identifiers** — not the same
   count. On v7 both orders give the same 93 records while removing `22 then 11` and `19 then 14`
   respectively (§18), so a count-only assertion is satisfied by both and says nothing.
2. `treating_centres(df) == treating_centres(df[eligibility.retained(df)])`. This is the property
   §4.5 derives from DECISION 1, and it is asserted separately because it is what *explains* the
   first: the two would come apart on exactly the frame Stage 4's E3 forbids, and a test that
   checked only the sets would not say which property had broken.

### 12.6 The class labels and the reason column — two AST scans

Over `cohort.py`, following Stage 1 9.4's and Stage 4 §12.6's technique, each with a companion test
proving the scan fires against a synthetic snippet that does contain the construct:

```
  ast.Constant   any str equal to a member of ELIGIBILITY_ORDER          → none
  ast.Subscript  of C.ELIGIBILITY_ORDER by an integer constant           → none
```

**The `Compare` scan of the earlier draft is dropped, deliberately** (§4.3): with two classes
`!= C.INELIGIBLE` and `== C.ELIGIBLE` select the same patients, so it is redundant rather than
dangerous, and a scan whose stated reason has evaporated teaches the next reader something false.

**What replaces it is a scan the earlier draft did not need**, because DECISION 1a changes what
`contraindication_reason` is allowed to do here. The column may reach a *table* or a *`detail`
string* and must never reach a *mask that decides a class*. So the third scan asserts:

```
  every reference to `contraindication_reason` in cohort.py occurs inside one of
      _flow_table          §7.3's last row: the arm-restricted no-reason count
      _flow_detail         §7.3's {u}
      _eligibility_detail  §7.2's {r} and {c}
```

walking the AST function by function rather than over the module as a whole. **Three functions, not
one**, and the list is exactly the three that render a count of undocumented patients — §7.2's detail
interpolates both `{r}` and `{c}` from that column, which an allowlist of `_flow_table` alone would
fail on the first honest implementation. The companion test fires the scan by pasting a reference
into `build`, which is the function the restrictions live in (§3) and therefore the one place the
column must never appear.

The rule is *how* the column is read, not *where* it may be mentioned: any of the three may count it,
none of them may branch a class on it, and a fourth function that starts reading it has to be added
here deliberately rather than by drift.

This is the one place Stage 5 can catch a partial revert of DECISION 1a — an implementer who
half-remembers "a blank is never read as no contraindication" and reintroduces `.notna()` into the
restriction. That edit would silently drop 43 patients and every table would still reconcile.

### 12.7 `derive_cohort` is called on the cohort, and this is where it is proved

**Two cells of the `core_above_median` audit entry, not one.** §4.5 is the argument; this is the
assertion.

```
                                     records cell     median (mL) cell
  cohort_frame()   cohort   (5)          "5"               "8"
                   after 1  (6)          "6"               "10"     ← median alone: caught
                   as given (9)          "9"               "10"     ← median alone: caught

  v7               cohort   (93)        "93"               "5"
                   after 1  (104)      "104"               "5"      ← median alone: MISSED
                   as given (126)      "126"               "6"      ← median alone: caught
```

- **The `records` cell equals `len(cohort)`** — `93` on v7, `5` on `cohort_frame()`. This is the
  assertion that holds on every frame, including the one where the median does not move.
- **The `median (mL)` cell is the cohort's**, asserted against the other two frames' medians as well,
  so the test states what it is discriminating rather than merely matching a number.
- **Both cells are strings** — `data._fmt` renders through `%.6g`, so the median cell is `"5"` and
  `"8"`, never `"5.0"`. Compare as `float(cell)` or against the rendered form; a test written against
  `"5.0"` fails on a correct implementation, which is the worst kind of red.

And: `build` is the only call site of `derive.derive_cohort` in `extended_bridging/**/*.py` outside
`derive.py` and its own test file — an AST or grep scan, because `derive.py`'s "Stage 5 is
`derive_cohort`'s only caller" is a sentence until something checks it.

### 12.8 The three tables

- `restrict_centres` renders `len(CENTER_ORDER) + 1` rows plus the header on **every** frame,
  including one where a declared centre contributes nothing, so an absent centre reads as zeros
  rather than vanishing. Never `pd.crosstab`. Its **centre rows'** `status` cells partition into
  `retained` / `dropped` and agree with `treating_centres`; the `all` row's cell is the summary
  `"{n} removed"` and is excluded from that assertion, because it is a total and not a status.
- `restrict_eligibility` renders one row per **retained** centre plus `all`; its `removed` column
  equals its `ineligible` column; its `all` row equals the entry's `n`; and its `all / eligible` cell
  equals `len(cohort)` — 93 on v7 — which is the identity §7.2 notes.
- `cohort_flow` renders exactly four rows plus the header on every frame. **Every `centres` cell is
  `frame["center"].nunique()` at that step** — the centres that contribute at least one record —
  never `len(CENTER_ORDER)` and never the retained tuple's length. The three agree on v7 by
  coincidence and disagree on the fixture, so the rule is declared rather than inferred from the
  table. Each row's arm cells sum to
  its `records` cell and its class cells sum to it too; the `after restriction 2` row equals the
  returned frame; the `of which` row counts **control-arm** records with no `contraindication_reason`
  and its `centres` cell is the number of centres contributing one — **2 on v7 against 3 in the row
  above**, and **1 against 2** on `cohort_frame()`, asserted as different numbers because that
  difference is §7.3's finding. The arm restriction is asserted directly as well: the same count over
  all arms is **82 on v7 and 4 on `cohort_frame()`**, and the test names those numbers so that
  dropping the restriction fails rather than merely inflating a cell nobody checks.
- On `cohort_frame()` all three are asserted cell for cell against §12.0's block.

### 12.9 The audit inventory

Exactly four new entries in `build`'s own frame, plus `derive_cohort`'s two, in this order:

```
  ("cohort",      "restrict_centres")          n = removed,  names them
  ("cohort",      "restrict_eligibility")      n = removed,  names them
  ("cohort",      "cohort_flow")               n = len(cohort), names nobody
  ("derivation",  "core_above_median")         derive_cohort's
  ("derivation",  "constant_covariates")       derive_cohort's
  ("missingness", "absence_by_cohort_column")
```

**Position is asserted against an index captured immediately before the `build` call, never against a
remembered total** — Stage 4 §12.9's rule, and the same trap `TODOS.md` records for `test_derive.py`.
`_record_removal`'s length check is fired directly, with a frame whose `case_id` has been duplicated.
The entries render under `## Cohort construction`, and `data.KINDS` is asserted to be the declared
**eight** in the declared order, so a future implementer cannot quietly add a ninth. Reproduction is
byte-identical within a process and across interpreters launched with `PYTHONHASHSEED=0` and `=1`.

The two-seed driver is written out rather than sketched, and note what it imports:

```python
script = ("import sys, pathlib; sys.path.insert(0, '.')\n"
          "import derive, eligibility, cohort, test_data, test_cohort\n"
          "df, audit = test_data.run(test_cohort.cohort_frame(),\n"
          "                          test_data.hand_source(pathlib.Path('/tmp')))\n"
          "df = eligibility.classify(derive.derive(df, audit), audit)\n"
          "cohort.build(df, audit)\n"
          "sys.stdout.write(audit.to_markdown())\n")
```

**It imports a test module in a subprocess, and that is the price of §12.0.** Stage 4's equivalent
drove `data.load(data.FIXTURE)`, which needs no test import — but the fixture cannot reach this stage
(§12.0), so the frame has to come from Python. Verified importable: the modules are flat in
`extended_bridging/` and `test_derive.py` already imports from `test_data.py` (§18).

Three things about the script that are not decoration. It calls **`test_data.run` and
`test_data.hand_source`, the helpers that already exist** — there is no `data.run_frame` and none is
to be added, since a second frame-driving entry point in `data.py` would exist only for this test.
`hand_source`'s path **need not exist**: nothing after the read touches the filesystem, and the log
header prints `path.name` and the label rather than opening anything. And the `build` call stands on
its own line with the render after it, so a `SchemaError` from `build` surfaces as a traceback rather
than being swallowed by a boolean expression into an empty string that then compares equal across
both seeds.

### 12.10 The cohort's denominators

`absence_by_cohort_column` exists, its `n` is `len(cohort)` and not `len(source)`, its rows are
`[*sorted(ANALYSIS_NAMES), *sorted(DERIVED_NAMES)]` in that order, and it includes
`core_above_median` — which is the one row that could only be produced after `derive_cohort`, so it
also pins the call order a second way. Its per-centre columns come from `CENTER_ORDER`, so the
dropped centre renders as an all-zero column rather than vanishing. (Strike this section with §7.4 if
that decision is struck.)

### 12.11 Stage 5 needs Stage 4 and nothing else new

`build` runs on a frame that has been through `classify` but **not** through `derive` — C1 is about
`eligibility`, and neither restriction, neither removal entry and no cell of the three tables reads a
derived column. It gets all the way to `derive_cohort` and raises there.

**What it raises is `KeyError('onset_type')`, not `SchemaError`, and the test must assert exactly
that** (verified, §18). `_core_above_median` succeeds, because `core_ml` is a Stage 2 analysis column
— it is `constant_covariates` that fails, ranging over `PS_COVARIATES_FULL`, whose `onset_type` term
is Stage 3's. A test written against a `SchemaError` and a helpful message would never go green, and
the temptation would be to "fix" `build` by adding a precondition that Stage 5 does not want: this
stage genuinely does not need `derive`, and the boundary is documented rather than enforced.

Assert also that the raise happens **after** three `cohort` entries and one `derivation` entry have
been recorded. That is not a defect — it is the ordinary consequence of a postcondition failing late
— but it is the difference between this and C1/C2, which are placed early precisely so the log is
never written. A caller who wants to restrict a Stage 2 frame has a sequencing bug; a caller who
wants a cohort calls `derive` first.

### 12.12 No bare `assert`, no raw header, not exempt

An AST scan of `cohort.py` for `ast.Assert`, as Stage 2 §12.10, with the companion test proving the
scan fires. Neither `cohort.py` nor `test_cohort.py` appears in `EXEMPT_FROM_RAW_NAME_SCAN` — asserted
directly, so the exemption cannot be granted quietly.

### 12.13 Structural facts from the workbook

Data-gated, against v7 and the §18 record, on `test_cohort.py`'s **own** module-scoped `workbook`
fixture (`load → derive → classify`, `DATA_GATED`) — `test_data.py`'s and `test_eligibility.py`'s do
not cross files:

- `treating_centres` returns `("HUG", "CHUV", "Lugano")`;
- the flow is `126 → 104 → 93`, with `39` bridging and `54` control in the cohort;
- the classification is `107` eligible and `19` ineligible over 126 records, and `93` / `0` in the
  cohort — **and the retained set is identical, patient for patient, to the one DECISION 1's
  three-class rule produced**, which is the assertion that pins the amendment as a relabelling of
  the cohort rather than a change to it (§18, §20);
- `43` retained **control-arm** records carry no `contraindication_reason`, at `2` centres, and are
  `43/54` of the control arm; over all arms the same count is `82` at `3` centres;
- the three tables reproduce §7.1, §7.2 and §7.3 cell for cell, including the all-zero
  `USZ / bridging` cell;
- the frozen median is `5.0` mL, the input frame's is `6.0`, and the frame between the two
  restrictions is **also `5.0`** — the coincidence §4.5 is written about, asserted here so that a
  future workbook in which it stops being true is noticed rather than assumed;
- the `core_above_median` entry's `records` cell is `93`, which is the assertion that holds when the
  median does not move;
- `constant_covariates` over the cohort is **empty** — so no [§6] covariate is lost to the
  restriction, which is the premise Stage 6's design-matrix builder rests on and which is not true of
  every frame in this file (it is 10 covariates on `cohort_frame()`);
- and the two order-invariance properties of §12.5.

### 12.14 Stage 5 removes rows, adds one column, and changes nothing else

On `cohort_frame()` through Stage 4, and again data-gated on the workbook:

```python
before = df.copy()
out = cohort.build(df, audit)

assert list(out.columns) == [*df.columns, *C.COHORT_DEPENDENT_SUBGROUPS]  # appended, and last
assert out.index.is_unique and out.index.isin(df.index).all()            # preserved, not reset
assert out.drop(columns=list(C.COHORT_DEPENDENT_SUBGROUPS)).equals(df.loc[out.index])
assert df.equals(before)                                                 # caller's frame untouched
```

The third line is the one that fails if any value is edited or any dtype "tidied" on the way past —
it compares the returned frame against the *same rows* of the input. The fourth is the one that fails
if `df.copy()` is moved below a write, and nothing else in §12 would notice. The second is the one
that fails if a future edit adds `reset_index(drop=True)`, which would otherwise change nothing
observable (§0.2).

### Coverage map

```
cohort.py                                          test_cohort.py
├── treating_centres(df)
│   ├── computed, over CENTER_ORDER                ├── 12.1  order is CENTER_ORDER's
│   ├── > 0 bridging patients                      ├── 12.13 ("HUG","CHUV","Lugano") on v7
│   └── invariant to restriction 2                 └── 12.5  the DECISION 1 consequence
│
├── _assert_cohort_inputs(df)
│   ├── C1 no eligibility column   → SchemaError   ├── 12.3  audit untouched by the raise
│   ├── C2 subgroup already present→ SchemaError   ├── 12.3  build called twice
│   ├── C3 ivt missing or ∉ {0,1}  → SchemaError   ├── 12.3  corrupt AFTER Stage 2
│   ├── C4 centre ∉ CENTER_ORDER   → SchemaError   ├── 12.3  incl. a MISSING centre
│   └── several at once → one message              └── 12.3
│
├── the two restrictions, inline in build (§3)
│   ├── restriction 1 by computed predicate        ├── 12.1  by identifier, never by count
│   ├── restriction 2 by eligibility.retained      ├── 12.2  COHORT-1 and COHORT-3 both retained
│   ├── undocumented -> eligible  [DECISION 1a]    ├── 12.2  reason blanked -> SAME cohort
│   ├── order-invariant membership                 ├── 12.5  identifiers, not counts
│   └── no label literal, no ORDER[i]              └── 12.6  two AST scans + both fire
│
├── _record_removal(...)
│   └── names exactly n cases, else SchemaError    └── 12.9  fired on a duplicated case_id
│
├── _flow_table / the three tables
│   ├── every declared centre rendered             ├── 12.8  incl. the all-zero bridging cell
│   ├── retained centres only, in table 2          ├── 12.8  all/eligible == len(cohort)
│   ├── four flow rows on every frame              ├── 12.8
│   ├── the no-reason line, ARM-restricted         ├── 12.8  centres ≠ the row above; and the
│   │                                              │         all-arm count (82 / 4) named
│   ├── `centres` cell is nunique at that step     ├── 12.8  declared, not inferred
│   ├── reason column reaches only the three       ├── 12.6  AST scan, per function
│   │   rendering helpers                          │
│   └── rows computed, never subtracted            └── 12.8  cell for cell on cohort_frame()
│
├── derive.derive_cohort(after_2, audit)
│   ├── called AFTER both restrictions             ├── 12.7  records cell == len(cohort): 93 / 5.
│   │                                              │         The median (5.0 vs 6.0) misses the
│   │                                              │         mid-position on v7 — §4.5
│   └── build is its only caller                   └── 12.7  repository scan
│
├── absence_by_column(out, ...)                    └── 12.10 n = 93; includes core_above_median
│
├── _assert_cohort(out, source, removed)
│   ├── P1 no ineligible remains                   ├── 12.4  + a local reimplementation
│   ├── P2 every retained centre both-armed        ├── 12.4  the fixture, hand_frame, and the
│   │                                              │         flagged-COHORT-1 frame
│   ├── P3 cohort non-empty, ≥ 1 centre            ├── 12.4
│   └── P4 removals + cohort == source             └── 12.4  the silent-unselect failure
│
├── build(df, audit)                                        written out in §8
│   ├── asserts BEFORE it copies                   ├── 12.3  audit untouched
│   ├── rows removed, one column added             ├── 12.14 §0.2's promise, asserted
│   ├── index preserved, never reset               ├── 12.14
│   ├── six entries in the declared order          ├── 12.9  position vs captured index
│   ├── renders under ## Cohort construction       ├── 12.9  KINDS asserted as EIGHT
│   ├── byte-identical, same process               ├── 12.9
│   ├── byte-identical, PYTHONHASHSEED 0 vs 1      ├── 12.9  driver imports test_cohort
│   └── needs classify, not derive                 └── 12.11
│
└── no bare assert / no raw header / not exempt    └── 12.12 + the scan fires

data.py amendments                                 test_data.py / test_eligibility.py
├── KINDS gains "cohort" in position               ├── 12.9 + the two renamed literal pins
├── _HEADINGS gains "Cohort construction"          ├── 12.9
└── _MUST_NAME_CASES unchanged                     └── the existing kind-parametrised tests
                                                       cover the new kind unedited

Data-gated (skipif not DATA_XLSX.exists()):        └── 12.13  126→104→93, 39/54, 43 at 2 centres,
  own module-scoped `workbook` fixture                        the three tables, 5.0 vs 6.0, and
  (load → derive → classify)                                  an empty constant-covariate tuple

Every branch above has a test. No branch is untested.
```

### Failure modes

| Failure | Test | Error handling | Visible or silent |
|---|---|---|---|
| A partial revert of DECISION 1a — `.notna()` reintroduced into restriction 2, dropping 43 patients | 12.2, 12.6 | the reason column is asserted to reach only `_flow_table`; and blanking it must not move the cohort | visible |
| A class label written as a literal, or the rule keyed on a display order | 12.6 | two AST scans, both with companions | visible |
| A centre's whole control arm removed by restriction 2, leaving a non-positive stratum in the propensity model | 12.4 | **P2**, naming the centre and the arm | visible |
| `derive_cohort` called before **both** restrictions — the [§13] subgroup split on 6.0 mL instead of 5.0 | 12.7 | `build` is its only caller and calls it last; median and `records` cell both asserted | visible |
| `derive_cohort` called **between** the restrictions — the subgroup split on 104 patients' median, which on v7 is 5.0 and therefore identical to the cohort's | 12.7 | the median cannot see this on v7; the entry's `records` cell reads 104 instead of 93 and is what catches it [§4.5] | visible **only** via the `records` cell |
| `build` called twice, giving a bootstrap replicate its own cut-point | 12.3 | **C2**, before any row is removed; `derive_cohort`'s own guard behind it | visible |
| A record with a missing centre or exposure silently unselected by a mask and attributed to restriction 1 | 12.3, 12.4 | **C4** and **C3** at run time, naming the cases; **P4** independently, summing against the input frame | visible |
| A restriction reporting a count it cannot name | 12.9 | `_record_removal`'s length check | visible |
| An ineligible patient surviving into [§7] — invariant 2 | 12.4 | **P1** | visible |
| An empty cohort returned rather than raised | 12.4 | **P3** | visible |
| The index reset, breaking traceability to the [§14] frame | 12.14 | asserted preserved and unique | visible |
| `build` mutating or reshaping the caller's frame | 12.14 | the four assertions of §12.14; `df.copy()` before the first write | visible |
| A centre × arm combination absent from the data vanishing from a table | 12.8 | every declared cell rendered; `pd.crosstab` is not used | visible |
| The flow table reconciling with itself because rows were subtracted | 12.8 | every cell computed from its own frame | visible |
| A new audit kind added without a spec | 12.9 | `data.KINDS` asserted as the declared eight | visible |
| The log made non-reproducible across processes | 12.9 | two-seed render comparison | visible |
| `test_derive.py`'s tail slice silently spanning the new entries | §10 | the repair lands in this commit | visible |
| A check written as `assert`, run under `-O` | 12.12 | AST scan at test time | visible |
| A control at CHUV or Lugano whose `ivt_contraindicated = 0` was an unfilled default rather than an assessment — 43 patients | none | **not detectable in the data** [DECISION 1a]. §7.3's last row keeps the group countable; §13 carries the assumption forward and [§3] must state it | acknowledged |

One row is undetectable by construction and is DECISION 1a's acknowledged assumption; every other
failure is visible at the point it occurs.

## 13. Known gaps carried forward

- **DECISION 1a moves an assumption rather than removing one, and the new one is less visible.**
  Under DECISION 1 the 43 undocumented controls carried a label that said so, in every table, and
  [§3] required the assumption stated as a limitation. Under DECISION 1a they carry the same label as
  a patient who was assessed and documented, and the assumption becomes: **`ivt_contraindicated = 0`
  means the same thing at CHUV and Lugano, which never collected a reason, as it does at HUG and USZ,
  which collected one for every control.** That is a claim about how the workbook was filled in, and
  nothing in the data can test it — the flag is 0/1 and never missing on all 126 records (§18), which
  is consistent both with it having been assessed for everyone and with 0 being an unfilled default.
  §7.3's last row is the only thing in the pipeline that keeps the group countable, and it is the
  reason that row was kept rather than deleted with the class it used to report. **This is now the
  single largest untested assumption in the analysis** and [§3] must say so in the same words it
  currently uses for the retention assumption.
- **The [§13] sensitivity analysis that would bracket it is more valuable after DECISION 1a, not
  less.** Restricting to controls with a documented reason gives 11 patients at HUG and collapses the
  design to one centre (measured — it is the option costed on 2026-08-13). So the honest bracket is
  probably not a re-restriction but a quantitative-bias or E-value-style argument. Deferred by
  decision (roadmap Stage 11); recorded here because the amendment changes what the deferred work
  should be, and §16 records why lifting the pilots' `ELIGIBILITY_PRIMARY_RULE` switch is the wrong
  way to un-defer it.
- **Restriction 1 is one-directional and P2 raises rather than resolving.** If a workbook arrives in
  which a retained centre has no control arm, the pipeline stops and the resolution is a [§3]
  amendment, not a code change (§5.2, §14).
- **The cohort is 93 records and the [§6] covariate list is 9 terms** — `PS_COVARIATES`, verified —
  which is ~11 design-matrix parameters once `center` expands to two dummies over the cohort's three
  centres and `onset_type` to two over its three levels, and therefore the ~3.5-treated-per-parameter
  budget Stage 0 confirmed against the 39 treated. Stage 5 does not check it: the budget is
  a property of the covariate list and the treated arm, both of which are Stage 1's and Stage 6's. It
  is recorded here because this is the stage that fixes the denominator of that ratio, and a
  restriction that removed more patients than expected would move it silently.
- **No `tests/fixture_cohort.xlsx`.** §12.0 declines it and Stages 6-13 may want it — a *file* whose
  read produces a both-armed cohort would let the whole pipeline run end to end with no `data/`.
  `cohort_frame()` gives them a *frame*, which is enough for every test this spec writes. Build the
  file when a stage needs a `data.Source` rather than a DataFrame, and give it its own spec section
  when it arrives. **Filed in `TODOS.md` at review**, with the two costs this spec did not have to
  price: it needs a **generator script** (the schema fixture has none, which is why §12.0 declined
  editing that one), and it is a Stage 1/Stage 2 amendment as well as a new file, because
  `test_data.py`'s `test_sources_are_exactly_the_workbook_and_the_fixture` pins `data.SOURCES`
  against a literal of exactly two. Whoever builds it must give it a **control at a treating centre**
  — the property `COHORT-1` carries — or it hits P2 exactly as the schema fixture does.
- **`_record_removal`'s length check guards a Stage 2 property.** It can only fire on a duplicated or
  missing `case_id`, which A2 forbids. It stays because this is the stage whose output is
  unreconstructable without it, but a reader should know it is a belt over Stage 2's braces and not an
  independent guarantee.

## 14. What Stage 5 deliberately does not decide

Recorded so a later stage does not look here for an answer that was never placed here.

- **What to do about a retained centre with no control arm.** P2 raises. Both resolutions — extend
  [§3] restriction 1 to the control direction, or keep the centre and accept a non-positive stratum —
  are [§3] amendments and belong to the PI. §5.2 states the trade; the code states it at run time.
- **Whether the undocumented group gets a sensitivity analysis, and of what kind.** [§13], deferred
  by decision at roadmap Stage 11 — and §13's second bullet records that DECISION 1a changes what the
  deferred work should be, since the obvious re-restriction collapses the design to one centre. Stage
  5 supplies the counts such an analysis would be sized against, in §7.3's last row.
- **Which population [§14a] and [§14b] draw.** Stages 12 and 13, from the *unrestricted* frame. Stage
  5 neither produces nor blesses them, and `build` is not a step on their path (§9).
- **How the [§3] retention assumption is worded in the manuscript.** Stage 14 [§16].
- **Whether the per-estimate denominators of [§11] are reported from §7.4's table or rebuilt in the
  reporting layer.** Stage 14. Stage 5 provides one table over the cohort; it does not claim to be
  the last word on denominators.

## 15. NOT in scope for Stage 5

| Considered | Why deferred |
|---|---|
| A `NEVER_IVT_CENTERS` or `BOTH_ARM_CENTERS` constant in `config.py` | §4.2. The predicate is a property of the data; a constant makes a workbook change require someone to remember. `pilots/pilot_config.py:78` is the prior art and §16 declines it |
| A parameter letting a caller pass a different centre list or a different eligibility rule | `pilots/analysis.py:164,171` has both. One cohort, one rule: a `rule="exclude"` switch is the `== eligible` reading under another name, and Stage 4 §5.1 abolished the choice |
| Dropping centres with no *control* arm | [§3] names one direction. P2 raises instead (§5.2, §14) |
| `reset_index(drop=True)` | §0.2. Traceability to the frame Stages 12 and 13 analyse |
| `remove_unused_categories` | §3.1. `center` is a `string` by Stage 2 §6's decision, taken for exactly this reason |
| Dropping constant covariates found in the cohort | Stage 3 `constant_covariates` detects and logs; Stage 6's design matrix drops at the point of use, which is also the only place that can see a bootstrap replicate having emptied a factor level. Two stages cannot both own the deletion |
| Any weighting, fitting or balance diagnostic | Stages 6 and 7. [§3] is explicit that both restrictions are "by design, **before** any weighting" |
| A complete-case restriction on the [§6] covariates | [§11] is complete-case *per estimate*, with each estimate carrying its own denominator. Restricting the cohort to complete cases here would give every estimate one denominator and silently change the population |
| A CLI entry point | Stage 14 is the single entry point [§16]. `cohort.py` is a library module |
| Caching the cohort | 126 rows, two masks. A cache is a second source of truth |

## 16. What already exists, and what to lift

`pilots/analysis.py`'s `target_trial_cohort` is the only prior art and it implements the same two
restrictions. The roadmap's standing warning applies — "This should not be a ground source […] Some
implementations may be wrong" — and here it is load-bearing: four of the eight rows below are *do not
lift*.

| From `pilots` | Status |
|---|---|
| `target_trial_cohort`'s docstring — two restrictions, both by design, with the structural-non-positivity argument | **lift the argument**, which is [§3]'s own and is already quoted in §4.1 |
| its `steps` dict — `start`, `excluded_never_ivt_centre`, `excluded_ineligible`, `analysed`, `n_bridging` | **lift the shape, not the mechanism.** §7.3's flow table is the same decomposition, in the audit log where it is reproducible and reviewable, rather than a return value the caller may drop |
| `C.BOTH_ARM_CENTERS` / `C.NEVER_IVT_CENTERS = ["USZ"]` (`pilot_config.py:73-78`) | **do not lift.** A hardcoded list makes the [§3] restriction a decision taken once against one workbook. §4.2 computes it. The pilots' own comment — "USZ administered IVT to 0 of 22 patients" — is a *finding*, and a finding written as a constant is a finding that cannot be re-derived |
| `centers = centers or C.BOTH_ARM_CENTERS` — the caller may override | **do not lift.** One cohort. An override is how two stages come to analyse different populations while both calling the same function |
| `rule = rule or C.ELIGIBILITY_PRIMARY_RULE`, with `"exclude"` and `"include"` branches over the indeterminate group | **do not lift**, and after DECISION 1a there is no longer a class for it to branch on. The pilots' `_eligibility` also classifies from the **free text**, against two hand-maintained reason lists with an `assert not unknown` — the approach DECISION 1 declined for the typo and annotation problems Stage 4 §4.1 records, and which DECISION 1a moves further from still |
| `out.reset_index(drop=True)` (`analysis.py:184`) | **do not lift.** §0.2 |
| `out["center"].cat.remove_unused_categories()` (`analysis.py:185`) | **do not lift** — and know why: it is needed only because the pilots made `center` a `Categorical` (`pilot_data.py:99`). Stage 2 §6 chose `string` for exactly this reason, so there is no dead level to remove |
| its silence on what was removed — `steps` carries counts, never identifiers | **do not lift.** §6.3: a cohort the log cannot reconstruct is a cohort nobody can check. `_record_removal` names every removed patient |
| its absence of any both-arms postcondition | **do not lift the absence.** P2 (§5.2). The pilots restrict on the treated direction and never check the control one |

`stage0_data_inventory.py` establishes *which* centres contributed no bridging patient (roadmap Stage
0, item 1) and its cohort-flow numbers are what §18 reconciles against. It computes; it does not
declare. That posture is lifted.

## 17. Implementation tasks

Ordered. Each independently verifiable. T1 is the `data.py` amendment; T2-T5 are Stage 5's own; T6 is
the amendment sweep.

**Strictly sequential — there is no split worth making.** Every task after T1 funnels through
`cohort.py` and `test_cohort.py`, and T1 has to land first because `Audit.record("cohort", …)` raises
`ValueError` until `KINDS` carries the kind, so nothing below it can be run even once. Do not attempt
to parallelise this across worktrees: the coordination costs more than the six tasks do, and T2-T4
each verify against acceptance tests that T5 writes. The one genuinely independent piece is T6's
`data.py` prose corrections, which touch no code path — fold them in rather than branching for them.

- [ ] **T1 (P1)** — `data.py`: `KINDS` and `_HEADINGS` gain `cohort` in §6.2's position; the
      `_MUST_NAME_CASES` comment gains its sentence. `test_data.py` and `test_eligibility.py`: the
      three literal pins and the two test renames of §10. Verify: `uv run pytest -v` green with
      **no other test edited**, which is where §10's claim that the two kind-parametrised tests need
      no change is checked rather than believed.
- [ ] **T2 (P1)** — `cohort.py`: `treating_centres` and `_assert_cohort_inputs` with C1-C4 and the
      collected raise. Verify: acceptance 12.3, and 12.1's `treating_centres` half.
- [ ] **T3 (P1)** — the two restrictions, `_record_removal`, and `_assert_cohort` with P1-P4. Verify:
      acceptance 12.1, 12.2, 12.4, 12.5, 12.6.
- [ ] **T4 (P1)** — the three tables of §7.1-7.3, `build` as written in §8 including the
      `derive_cohort` call and §7.4's absence entry. Verify: acceptance 12.7, 12.8, 12.9, 12.10,
      12.11, 12.14.
- [ ] **T5 (P2)** — `test_cohort.py` remainder: `cohort_frame()` (§12.0), its own module-scoped
      `workbook` fixture, the four AST scans and their companions, the two-seed driver, and the
      data-gated §12.13 facts. Verify: `uv run pytest -v` green both with and without `data/`.
- [ ] **T6 (P1)** — the §10 sweep's remainder: `test_derive.py`'s tail slice, `data.py`'s three prose
      corrections (A9's message, the `_missingness` comment block, the `_MUST_NAME_CASES` count), and
      striking the `TODOS.md` item that defers the slice. Verify: `test_derive.py` green, and the
      repaired assertion **seen to fail** and then to pass — which is the check that the repair
      addresses the trap rather than merely moving it.

      **How to run that demonstration, because the obvious way does not work.**
      `_full_pipeline()` loads `data.FIXTURE`, and the fixture cannot reach the end of `build`: it
      collapses to HUG with one bridging patient and no control arm, so **P2 raises** (§12.0). A
      `build` call inserted there throws rather than failing an assertion — and inserted *above* the
      existing `derive_cohort(df, audit)` line it raises C2 instead, since `build` already called it.
      Either way the demonstration proves nothing. So: build a **throwaway local pipeline** on
      `test_cohort.cohort_frame()` — `run → derive → classify → build` — confirm the old
      `[-len(_STAGE_3_INVENTORY):]` slice now returns the three `cohort` entries plus
      `derive_cohort`'s two and `absence_by_cohort_column`, and therefore fails; confirm the repaired
      index-captured assertion passes on it; then discard the throwaway. **`_full_pipeline()` itself
      is not changed** — six tests in `test_derive.py` assert against its current shape, and widening
      this commit into that file buys nothing the throwaway does not.

      **The roadmap is not in this sweep.** Its `**Spec:**` line and its rewritten **Accept when**
      land with this document, as Stages 3 and 4 did (§19).

### Diagrams that belong in the code, not only here

Two, and no more — a diagram nobody maintains is worse than none, because it is believed. Keeping
them true is part of any change that touches them, in the same commit.

- **`cohort.py`'s module docstring** — §0's flow, and specifically the two arrows: the restricted path
  into Stages 6-11, and the *unrestricted* frame going past this stage into Stages 12 and 13. The
  second is why this module removes rows from a copy and not from the pipeline.
- **Above `build`** — §4.5's order block: the two orders, the counts each attributes, and the one
  sentence that the retained set does not depend on the order but the attribution does. That is the
  line an implementer "simplifying the two masks into one" deletes.

One further comment, not a diagram: **above the `derive_cohort` call, the three medians and what they
mean** —

```
#   126 records, as classified   6.0 mL
#   104, after restriction 1     5.0 mL   <-- equal to the cohort's on v7, so the median
#    93, the cohort             5.0 mL       CANNOT tell this position from the right one.
#                                            The audit entry's `records` cell can: 93. [§4.5]
```

It is the whole of why the call is on this line and not an earlier one, it is invisible from the call
itself, and the middle row is the one an implementer deletes as redundant. §21 R1 is why it is not.

### Definition of done

Stage 5 is complete when all of the following hold, and not before.

1. `uv run pytest -v` is green **with** `data/` present — every test, including the data-gated ones.
2. `uv run pytest -v` is green **with `data/` temporarily renamed** — the data-gated tests skip and
   nothing else fails or errors at collection.
3. `test_config.py`, `test_data.py`, `test_derive.py` and `test_eligibility.py` are still green, and
   neither `cohort.py` nor `test_cohort.py` has been added to `EXEMPT_FROM_RAW_NAME_SCAN`.
4. The pipeline's audit log has been produced through all five stages, against the **workbook**, and
   reproduces byte-identically across hash seeds:

   ```bash
   cd extended_bridging
   S='import sys, data, derive, eligibility, cohort
   df, audit = data.load(data.WORKBOOK)
   df = derive.derive(df, audit)
   df = eligibility.classify(df, audit)
   df = cohort.build(df, audit)
   sys.stdout.write(audit.to_markdown())'
   PYTHONHASHSEED=0 uv run python -c "$S" > /tmp/a.md
   PYTHONHASHSEED=1 uv run python -c "$S" > /tmp/b.md
   diff /tmp/a.md /tmp/b.md
   ```

   `diff` reports nothing, and the three rendered tables match §7.1, §7.2 and §7.3 cell for cell. The
   log names patients and must never be committed.
5. Acceptance 12.7 has been **seen to fail in both directions**, and the two are different tests of
   different cells (§4.5):
   - with the `derive_cohort` call moved to the **top of `build`**, above restriction 1: the median
     becomes 6.0 mL and the `records` cell 126;
   - with it moved **between the two restrictions**: on v7 the median stays 5.0 and only the
     `records` cell moves, to 104. Watch this one specifically. It is the failure the median alone
     does not see, and if the `records` assertion is ever dropped it is the failure that ships.

   Then seen to pass when the call is restored. A green test that has never been watched failing is
   not evidence.
6. Acceptance 12.6's third scan — the `contraindication_reason` allowlist — has been **seen to fail**
   against a reference pasted into `build`, and seen to pass with the three rendering helpers
   (`_flow_table`, `_flow_detail`, `_eligibility_detail`) reading the column as §7 specifies. And the
   partial-revert it exists for has been watched: `.notna()` reintroduced into restriction 2 takes
   the cohort from **93 to 50** on the workbook, with every table still reconciling.

   *(There is no `Compare` scan. §4.3 declines it deliberately: under two classes `!= INELIGIBLE` and
   `== ELIGIBLE` select the same 93 records — verified — so the construct is redundant rather than
   dangerous, and a scan whose stated reason has evaporated teaches the next reader something false.)*
7. Acceptance 12.4's P2 half has been **seen to fail** with P2 removed, on the flagged-`COHORT-1`
   frame: `build` returns a cohort in which HUG has a bridging arm and no control arm, every table
   reconciles, and nothing raises.
8. Acceptance 12.4's P4 half has been **seen to fail** with P4 removed, on the missing-centre frame
   with C4 also removed: one record leaves the analysis, the flow table reconciles, and its identifier
   appears among restriction 1's removals.
9. Acceptance 12.3's C1 and C2 halves have each been **seen to fail** with the check removed — C1
   raising three entries later from `retained`'s `KeyError`, C2 raising from `derive_cohort` after the
   log has already been written. Both leave an audit describing a cohort no caller received, which is
   the failure the two checks are placed early to prevent.
10. Acceptance 12.14 has been **seen to fail** with `df.copy()` moved below the first write, and
    separately with `reset_index(drop=True)` appended. Nothing else in §12 notices either, which is the
    point.
11. Acceptance 12.9's `_record_removal` check has been **seen to fire** on a duplicated `case_id`, and
    12.12's `assert` scan against a pasted `assert True`. Both pass trivially if they match nothing.
12. `test_derive.py`'s repaired inventory assertion has been **seen to fail** with the tail slice
    restored, on the throwaway `cohort_frame()` pipeline T6 describes — **not** by inserting `build`
    into `_full_pipeline()`, which raises P2 on the fixture and demonstrates nothing.
13. The audit log has been read end to end by a human, and the flow numbers agree with §18 and with
    `../out/stage0_data_inventory.md`.
14. `grep -rn "derive_cohort" extended_bridging/*.py` returns call sites in `cohort.py`, `derive.py`
    and the two test files only — and `derive.py`'s "Stage 5 is `derive_cohort`'s only caller" is
    therefore true of the repository rather than of the roadmap.

## 18. Verification record

Checked on 2026-08-13, before this spec was finalised, by a read-only probe reporting aggregate counts
only — no outcome by arm, in the posture Stages 0 and 4 used. Recorded so a reader can tell which
numbers were verified rather than carried, and so the checks are re-runnable after a workbook update.
Acceptance 12.13 re-checks them in code.

The probe drives `load → derive → classify` against `data.WORKBOOK` and applies §4's masks; it writes
nothing.

| Claim | Where used | Verified |
|---|---|---|
| Bridging patients per centre: HUG 30, CHUV 7, Lugano 2, **USZ 0** — so `treating_centres` is `("HUG", "CHUV", "Lugano")` | §4.2, §7.1, §12.13 | yes, run |
| Records per centre: HUG 52, CHUV 21, Lugano 31, USZ 22 | §7.1 | yes, run |
| The flow is **126 → 104 → 93**, removing 22 then 11 | §4.5, §7.3, §12.13 | yes, run |
| The cohort is **39 bridging and 54 control** at 3 centres | §7.3, §11, §12.13 | yes, run |
| Under DECISION 1a the classification is **107 eligible / 19 ineligible**, and the 43 formerly `indeterminate` records all move to `eligible` and nothing else moves | §4.1, §12.13, §20 | yes, run |
| **The retained set is unchanged by the amendment**: the primary cohort is the same 93 records under both rules, compared identifier for identifier — 39 bridging, 54 control, HUG 30/11, CHUV 7/14, Lugano 2/29 | Status, §12.13, §20 | yes, run, compared by identifier |
| The revealed-fact rule becomes **redundant** under DECISION 1a: no treated patient carries the flag, so `flag == 0` already classifies all 39 as eligible. Stage 4's E3 is what keeps it redundant | §20 | yes, run |
| **43** cohort **control-arm** records carry no `contraindication_reason`, at **2** centres (CHUV and Lugano) — 43/54 = 80% of the control arm. Over **all arms** the count is **82** at **3** centres, because no treated patient carries a reason either | §7.3, §12.8, §12.13, §13 | yes, run |
| A reason is recorded on **44** of 126 records — exactly the HUG and USZ controls — and after DECISION 1a that column is read by no classifier | §4.1, §7.3, §12.6 | yes, run |
| Restricting the cohort to controls with a **documented** reason gives **11** patients, all at HUG, leaving CHUV and Lugano single-armed and the design single-centre at n=41 | §13 | yes, run |
| §7.1's five-row table, cell for cell, including the all-zero `USZ / bridging` cell | §7.1, §12.13 | yes, run |
| §7.2's four-row table, cell for cell — HUG 41 eligible / 11 ineligible of 52, CHUV 21/0 of 21, Lugano 31/0 of 31, all **93/11 of 104** | §7.2, §12.13 | yes, re-run 2026-08-13 at review. *(The row this replaces read "HUG 41/0/11, CHUV 7/14/0, Lugano 2/29/0, all 50/43/11" — DECISION 1's three-class table, carried over from the earlier draft and contradicting §7.2's own rendering. §21 R5.)* |
| §7.3's four-row flow table, cell for cell | §7.3, §12.13 | yes, run |
| The reverse order gives `126 → 107 → 93`, removing 19 then 14, and **the same 93 identifiers** | §4.5, §12.5 | yes, run, compared by identifier |
| `treating_centres` returns the same tuple before and after restriction 2 | §4.5, §12.5 | yes, run |
| Per retained centre in the cohort: HUG (30, 11), CHUV (7, 14), Lugano (2, 29) — **every centre both-armed**, so P2 is unreachable on v7 | §5.2, §12.13 | yes, run |
| The core-volume median is **6.0 mL** over the 126-record frame and **5.0 mL** over the cohort | §0.1, §4.5, §12.7, §12.13 | yes, run |
| `derive_cohort` takes the cohort from 32 to **33** columns, 93 rows | §0, §9, §11 | yes, run |
| `constant_covariates` over the cohort is **empty** | §12.13 | yes, run |
| The pipeline holds **12** entries on the workbook and **11** on the fixture **before Stage 5** — `load` 7 / 6, `derive` 4, `classify` 1. (The 14 / 13 this replaces counted `derive_cohort`'s two, which only `build` triggers, so it was a count *after* Stage 5 labelled *before* it. §21 R5) | §9, §12.9 | yes, re-run 2026-08-13 at review |
| The core-volume median over the frame **between** the two restrictions (104 records) is **5.0 mL** — identical to the cohort's, so a `derive_cohort` misplaced there is invisible to the median on v7. On `cohort_frame()` the same three medians are 10.0 / 10.0 / 8.0 | §4.5, §12.7, §12.13, §21 R1 | yes, run at review |
| On a frame where restriction 1 keeps nothing, the empty frame reaches P3 without an earlier exception: `derive_cohort` returns 0 rows and 33 columns with `median (mL)` rendering `missing`, and `absence_by_column` guards its own division | §5.1 P3, §12.4 | yes, run at review |
| `build` on a classified-but-not-derived frame raises **`KeyError('onset_type')`** from `constant_covariates`, not `SchemaError`, after three `cohort` entries and one `derivation` entry are recorded | §12.11 | yes, run at review |
| Every `center` value on v7 is in `CENTER_ORDER`, and `ivt` is 0/1 and never missing — C4's and C3's premises, and why every branch of §4.4 is unreachable | §4.4, §12.13 | yes, run |
| No treated patient is `ineligible` on v7 — DECISION 1's consequence that §4.5 rests on | §4.5 | yes, run |

**Claims about the test frames**, all run against the committed repository, so they are re-checkable
with no `data/`:

| Claim | Where used | Verified |
|---|---|---|
| `tests/fixture_schema.xlsx` restricts to **1 record at HUG with no control arm** — restriction 1 keeps HUG only, and P2 fires | §12.0, §12.4 | yes, run |
| `hand_frame()` restricts to **3 records, HUG (2, 0) and CHUV (1, 0)** — both centres control-less, so P2 fires naming two centres | §12.0, §12.4 | yes, run |
| `cohort_frame()` as specified in §12.0 gives **9 → 6 → 5**, cohort `{HAND-1, HAND-2, HAND-5, COHORT-1, COHORT-3}`, HUG (2, 1) and CHUV (1, 1) — 8 eligible / 1 ineligible over the 9, and 5 / 0 in the cohort. The cohort is the same set under DECISION 1 and DECISION 1a | §12.0, §12.1, §12.8 | yes, run, both rules |
| Its `of which EVT alone, no reason on file` cell is **1 record at 1 centre** in the cohort, against **2 centres** in the row above; over all arms the same count is **4**, and before restriction 1 it is **5 at 3 centres** | §12.0, §12.8 | yes, run |
| Its core-volume median is **10.0** over the 9 records, **10.0** over the 6 after restriction 1, and **8.0** over the cohort — so unlike v7, this frame *does* distinguish all three call positions by median alone, and it is the frame §12.7's median half rests on | §12.0, §12.7, §4.5 | yes, run at review |
| `cohort_frame(contraindication_reason=None)` — the reason column blanked on **every** record — yields the **same cohort, identifier for identifier**: `{HAND-1, HAND-2, HAND-5, COHORT-1, COHORT-3}`. This is §12.2's pointed half, run rather than reasoned | §12.2, §12.6 | yes, run at review |
| It passes Stage 2's `_pipeline_after_read` unchanged — A1-A9 all green, `HAND-4` as template keeps the volumes internally consistent | §12.0 | yes, run |
| Its dtypes are identical to `hand_frame()`'s, and construction is **warning-free** under `warnings.simplefilter("error")`; the `pd.concat` form is **not**, emitting a `FutureWarning` about all-NA columns on pandas 2.3.3 | §12.0 | yes, run, both forms |
| `derive_cohort` on its cohort gives 33 columns and 5 rows; `constant_covariates` there is 10 covariates, unlike the workbook's empty tuple | §12.13 | yes, run |
| A subprocess launched from `extended_bridging/` can `import test_data` with `sys.path.insert(0, '.')` | §12.9 | yes, run |
| `tests/fixture_schema.xlsx` has **no generator script** in the repository | §12.0 | yes, searched |

**Claims checked by reading the committed repository:**

| Claim | Where used | Verified |
|---|---|---|
| `data.KINDS` is the seven of §6.2 in that order; `_HEADINGS` keys equal `KINDS`; `_MUST_NAME_CASES` is `{correction, observation}` | §6 | yes, read |
| `test_data.py:364` pins `KINDS` and `list(_HEADINGS.values())` against literals, in a test whose **name** contains the word `seven` | §10, T1 | yes, read |
| `test_eligibility.py:630`'s `test_data_kinds_is_unchanged_by_this_stage` pins the same tuple **and `_MUST_NAME_CASES`**, and its comment predicts this stage inserting a kind. (It was at :609 before the DECISION 1a sweep moved it) | §6.1, §10 | yes, re-read 2026-08-13 at review |
| `test_data.py:328,339`'s two kind-parametrised tests range over `data.KINDS` and need no edit; the `else` branch asserts an entry with `n = 126` and no identifiers is legal | §6.3, §10 | yes, read |
| `test_data.py:719`'s heading-adjacency test takes its successor from `KINDS` rather than by name, so the insertion moves nothing there | §6.1 | yes, read |
| `test_derive.py:714` reads `[-len(_STAGE_3_INVENTORY):]`, and `_full_pipeline()` is `load → derive → derive_cohort` with no `classify` and no `build` | §10, T6 | yes, read |
| `TODOS.md`'s second item defers that repair to "the commit that first builds the real pipeline order" | §10, §19 | yes, read |
| `derive.py`'s docstrings already carry the post-Stage-4 counts (31 → Stage 4 → Stage 5, `derive_cohort` "33 columns"), so Stage 5 amends none of them | §10 | yes, read |
| `derive_cohort` raises `SchemaError` if `COHORT_DEPENDENT_SUBGROUPS` are already present, and `constant_covariates` defaults to `PS_COVARIATES_FULL` and raises `KeyError` on a pre-`derive` frame | §4.4 C2, §12.11 | yes, read |
| `absence_by_column` is public, takes `(df, audit, columns, step, detail)`, and builds its per-centre columns from `CENTER_ORDER` | §7.4, §12.10 | yes, read |
| `Audit.record`'s signature is `(kind, step, n, detail, case_ids=(), table=None)` and sorts identifiers on the way in | §6.3, §7 | yes, read |
| `ANALYSIS_NAMES` is a `frozenset` and `DERIVED_NAMES` a tuple, so §7.4's list is `sorted` on the first | §7.4 | yes, read |
| `pilots/pilot_config.py:73-78` declares `BOTH_ARM_CENTERS` and `NEVER_IVT_CENTERS` as literals; `pilots/analysis.py:164,171,184,185` take a centre override, an eligibility-rule switch, `reset_index(drop=True)` and `remove_unused_categories` | §4.2, §15, §16 | yes, read |
| `pilots/analysis.py`'s `steps` dict records counts and never identifiers | §16 | yes, read |
| SAP [§3] names only the bridging direction of the centre restriction | §4.2, §5.2 | yes, read |
| `data.py`'s A9 message still says a blank "moves a patient from indeterminate to eligible [DECISION 1, §3]", and the `_missingness` comment block still calls the reason column's absence "the thing that separates eligible from indeterminate" — both describe a withdrawn class, and neither was touched by the DECISION 1a sweep | §10, §21 R7 | yes, read at review |
| `eligibility.classify` reads `contraindication_reason` to count undocumented controls for its own `detail`, so Stage 5 is **not** the only reader of that column after DECISION 1a — the rule is how it is read, not where | §7.3, §12.6 | yes, read at review |
| The whole suite is green at **490 tests** with `data/` present, on the tree carrying the landed DECISION 1a sweep | Status, §20 | yes, run at review |

## 19. What this spec changed elsewhere

Five entries, all landing with this document rather than with the implementation, because a
specification that contradicts the roadmap is worse than no specification.

| # | Change | Files | The question it settles |
|---|---|---|---|
| 1 | Roadmap Stage 5 gains its `**Spec:**` line | `implementation_roadmap.md` | Stage 5 was the next stage without one |
| 2 | Roadmap Stage 5's **Accept when** gains the both-arms criterion as a **runtime raise**, not a test-time check, with the note that it can only fire because of restriction 2 | `implementation_roadmap.md`, §5.2 here | Whether "every centre in the resulting cohort contains both arms" is something a test confirms about v7 or something the pipeline enforces about any workbook. It is the latter, and on v7 it is unreachable — so a test-time reading would have shipped an unenforced criterion |
| 3 | Roadmap Stage 5's **Accept when** gains the ordering requirement: `derive_cohort` is called by `build`, after both restrictions, and `build` is its only caller — enforced at review on the **row count** rather than only on the median | `implementation_roadmap.md`, §4.5 here | Whether the [§13] median is frozen on the cohort. Stage 3 already required it and left the enforcement to a caller nobody had written; this is that caller. 6.0 mL over 126 records — but **5.0 over both** the 104-record intermediate frame and the 93-record cohort, which is why the criterion names the `records` cell |
| 4 | The `cohort` audit kind is added to `data.py` — the question Stage 4 §9.1 deferred to this stage — and `_MUST_NAME_CASES` is deliberately not extended, with the naming rule enforced by `_record_removal` instead | §6 here; `data.py`, and the three literal pins in §10 | Whether a stage that removes rows gets its own heading, and whether a kind-keyed rule can express "the entries that remove patients must name them" when the kind also carries one that removes none. It cannot |

| 5 | Roadmap Stage 5's "that amendment has **not yet landed**" paragraph is **struck**, and `TODOS.md` gains the `tests/fixture_cohort.xlsx` item | `implementation_roadmap.md`, `TODOS.md` | Both landed at review on 2026-08-13, ahead of the implementation. The first was false within hours of being written (§21 R4); the second is a gap Stages 6-13 will hit and which lived only in §13 (§13, §21) |

A sixth change is DECISION 1a's, and it is upstream rather than here: §20.

**What this spec has now had, and what it still has not.** It has had an engineering review —
2026-08-13, §21, eight findings, all folded in, every numeric claim re-run against the workbook and
against `cohort_frame()`. Stage 4's review (its §20) found ten defects in a document of this shape,
all by running code rather than by reading; this one found eight the same way, which is the expected
yield and not a reason to think the document is now clean.

**It has still had no cross-model second opinion, and that is now a decision rather than an
accident.** The Codex pass was attempted at review on 2026-08-13 and failed on authentication (HTTP
401). Rather than block on re-authenticating, the question "is a second reader worth it here" was put
and answered, and the answer is recorded because the next stage will face it too.

**The evidence, from this repository rather than from principle.** The outside voice ran on
2026-08-08 (errored) and twice on 2026-08-10. It did **not** run for Stage 4 — the review log carries
no cross-model entry for 2026-08-12. Stage 4's review was single-model, and its §20 opens: *"Ten
defects, every one of them reproduced by running code … not by inspection."* This review has the same
profile: seven of its eight findings came from executing something, one from reading the document
against itself. And what Stage 4's review still missed was caught by **implementation**, not by
review — `TODOS.md`'s third item records two false statements in the Stage 4 spec surfaced *"while
working Definition of done 9 — the one item of the thirteen that could not be satisfied."* Both were
claims about pandas semantics, which is to say claims only execution can settle.

So the marginal value of a second *reader* is low for everything measurable in this document, because
every numeric claim in §7, §12 and §18 has now been re-run. **The generic pass is deferred, and the
Definition of done's seven seen-to-fail items are what stands in for it** — that mechanism is the one
that caught Stage 4's residual defects, and it is why those items are written as instructions to
break the code and watch, rather than as assertions to trust.

**What a second reader is still worth, and when.** Three questions, none of which execution can
settle, and all of which are better asked once `cohort.py` exists than of the spec alone — a reviewer
with the code in front of them is worth more than one with only the document:

- **§6.2's position argument.** It is a judgement about a rendered document, it is wrong about two
  entries by its own admission, and it is the one decision here that is cheap now and expensive after
  a log has been circulated. The review left it untouched.
- **§12.0's `cohort_frame()`**, which every acceptance test below it depends on. Its numbers are
  verified twice over (§18) but its *design* — which record carries which property — is what makes
  §12.1 through §12.8 able to distinguish a correct implementation from a plausible one.
- **Restriction 1 being one-directional** (§4.2, §5.2, §14). [§3] names only the bridging direction
  and P2 raises on the mirror case rather than resolving it. That is a methodological call and it is
  the PI's; it is listed here so a second reader knows it is *deliberately* open rather than
  overlooked.

One further item the review checked but did not settle: **§8's `build`**, specifically whether
`_record_removal` and `audit.record` sitting side by side reads as an inconsistency rather than as
§6.3's distinction. The `source` binding P4 rests on was checked and is right; the readability
question is a matter of taste and was left alone.

## 20. DECISION 1a — what it required upstream, and where each item landed

**The decision, PI, 2026-08-13.** Eligibility is classified from `ivt_contraindicated` alone. There
are **two** classes: flag = 1 → `ineligible`, flag = 0 → `eligible`. A patient whose contraindication
reason was never documented is **eligible**. `Contraindications_to_IVT` is not read at all — neither
its text (which DECISION 1 already forbade) nor its presence (which DECISION 1 used as one bit).

It amends DECISION 1 and it amends [§3]. It is recorded here because Stage 5 is the first document
written against it, and because §18 verified its consequences; it must be transcribed into
`../out/stage0_data_inventory.md` alongside DECISION 1, which is where the pipeline's decisions live.

### 20.1 What it does not change

**The cohort.** Verified identifier for identifier (§18): the retained set is the same 107 records,
the primary cohort the same 93 patients, the same 39 bridging and 54 control, the same per-centre
split. Every flow number in §7 and §18, the [§13] median of 5.0 mL, `constant_covariates` being
empty, and P2 being unreachable on v7 are all untouched.

This is because DECISION 1 already *retained* the indeterminate group. The amendment relabels 43
patients; it does not move them. Stage 4 §5.1's decision of 2026-08-10 — retained means not
ineligible — turns out to have chosen the same population DECISION 1a arrives at by a shorter route,
and that is why the amendment is cheap.

### 20.2 What it required, in dependency order — all nine landed 2026-08-13

**This is a record of work done, not a list of work to do.** Every row below is in the tree and
`uv run pytest -v` is green at 490 tests; verified at review (§18, §21). It is kept rather than
deleted for the reason DECISION 1 itself is kept in `../out/stage0_data_inventory.md`: a dated
account of what a protocol amendment cost across nine files is the thing a reviewer asks for in six
months, and it does not survive in a diff. Stage 5 was blocked on rows 3 and 4 in particular —
`cohort.py` reads `ELIGIBILITY_RETAINED` and `_flow_table` renders `ELIGIBILITY_ORDER` — and is not
blocked now.

| # | Where | Change |
|---|---|---|
| 1 | `statistical_analysis_plan.md` [§3] | The passage "Where the reason a patient did not receive IVT was never recorded, eligibility is *indeterminate*; those patients are retained. A blank is never read as 'no contraindication'…" is **replaced**. The new rule is the flag; the new limitation is §13's first bullet — that the flag's meaning at the two centres which never collected a reason is an assumption the data cannot test, and it governs 43 of 54 controls. The limitation does not disappear with the class; it changes shape and must be restated, or [§3] will read as though the problem were solved |
| 2 | `../out/stage0_data_inventory.md` | DECISION 1a recorded beside DECISION 1, dated, with §18's verification that the cohort is unchanged. DECISION 1 is **not** deleted — it is a dated record of a decision taken, and Stage 4 §13's third bullet already establishes that posture for superseded citations |
| 3 | `config.py` | `INDETERMINATE` removed; `ELIGIBILITY_ORDER` becomes `(ELIGIBLE, INELIGIBLE)`; `ELIGIBILITY_RETAINED` becomes `(ELIGIBLE,)`. The block comment of Stage 4 §5.2 is rewritten — its "three classes are declared individually" argument survives for two, but its reason for the third does not |
| 4 | `eligibility.py` | `_classify` collapses to one mask: `eligible` everywhere, `ineligible` where the flag is 1. **The revealed-fact mask is deleted** — verified redundant, since no treated patient carries the flag (§18) — and E3 is what keeps it redundant, so E3 becomes load-bearing where it was previously belt-and-braces. `contraindication_reason` is no longer read: the module stops reading it, and Stage 4 §4.5's "`.notna()` is the whole of the reason column" section goes with it |
| 5 | `eligibility.py` | **Stage 2's A9 stops being load-bearing for Stage 4.** Stage 4 §4.5 records A9 — the whitespace-only reason check — as "the only Stage 2 assertion that is", precisely because a blank moved a patient between classes. It no longer can. A9 stays as a data-quality check; the dependency note must be struck or it will be believed |
| 6 | `specs/stage4_eligibility_classification.md` | The larger edit. §4.1's four cases become two; §4.4's ordering argument and its 39-record claim go with the revealed-fact mask; §4.5 goes entirely; §5.1's retained-predicate decision becomes a two-class identity and its 43-record argument evaporates; §7.1's crosstab loses a column; §12.1, §12.2, §12.7 and §12.13 are rewritten; §18's rows are re-run. It is a dated, reviewed document, so the amendment wants a **§21 amendment section** rather than a silent rewrite — the same treatment [§8]'s rare-outcome amendment got |
| 7 | `test_eligibility.py`, `test_config.py`, `test_data.py` | Every assertion naming `indeterminate`, the 64/43/19 split, the 39-record revealed-fact test, and `ELIGIBILITY_ORDER`'s length. The 39-record test is the one to watch: it asserts a mask that is being deleted |
| 8 | `implementation_roadmap.md` | Stage 0's DECISION 1 and Stage 4's **Build** and **Accept when**. Stage 5's entry was rewritten with this document. **Its "that amendment has not yet landed" paragraph was struck at review** — it outlived the thing it described by a few hours and was the loudest false statement in the roadmap (§21 R4) |
| 9 | this document | Written against DECISION 1a throughout, and re-checked against the landed code at review |

### 20.3 The one thing to get right in that sweep

**The 43 patients must not become invisible.** Under DECISION 1 they carried a label, so every table
in the pipeline showed them and [§3] required a limitation about them. After the sweep they are
indistinguishable from documented-eligible controls in every column of the frame, and the only
things that still count them are §7.3's last row and whatever [§3] is rewritten to say.

Items 1 and 6 are therefore not bookkeeping. A sweep that deletes the class and forgets the
limitation converts a stated assumption into an unstated one, which is the failure the whole
three-class apparatus was built to prevent — and it would do it while every test passes and every
denominator reconciles.

## 21. Engineering review, 2026-08-13

Run before a line of `cohort.py` was written, on the tree carrying the landed DECISION 1a sweep.
**Method: every numeric claim in §7, §12 and §18 was re-run** — against the workbook through
`load → derive → classify`, against `tests/fixture_schema.xlsx`, against `hand_frame()` normalised
and bare, and against `cohort_frame()` built exactly as §12.0 specifies. Eight findings, all folded
into the sections above. Seven of the eight came from running something; one (R2) from reading the
document against itself.

Single-model. The cross-model pass failed on authentication and did not run; **deferring it rather
than blocking on it is a decision with its own reasoning, and §19 carries both** — the evidence being
that Stage 4 shipped after a single-model review of the same shape, and that what *that* review
missed was caught by implementation rather than by a second reader.

| # | Finding | Severity | Resolution |
|---|---|---|---|
| R1 | **The frozen-median guard missed the likeliest misplacement.** §12.7 claimed the median assertion "fails whichever direction the call moves". Measured: 6.0 mL over 126 records, **5.0 over 104**, 5.0 over 93 — so a `derive_cohort` called *between* the restrictions is invisible to it on v7, and that is the position a reader tidying the function would choose | P1 | §12.7 now asserts the entry's **`records` cell** (93, and 5 on `cohort_frame()`) as well as the median. A row count cannot coincide the way a median can. §4.5 carries the three-frame table; DoD-5 watches both directions fail |
| R2 | **`_restrict_centres` and `_restrict_eligibility` were named in §0, §3, §12.6 and the coverage map, and defined nowhere.** §8's `build` — presented as canonical — inlines both as plain masks. An implementer following §8 produced a module missing two of the test plan's scan targets | P1 | Both names struck. §3 says explicitly that neither exists and why: §8's twelve lines are the readable form, and extracting them would thread `after_1` / `after_2` out through return values for no gain. §12.6's companion now pastes into `build` |
| R3 | **The `contraindication_reason` scan allowlisted one function; §7 needs three.** §7.2's detail interpolates `{r}` and `{c}` and §7.3's `{u}`, all read from that column, so the scan would have failed against the first honest implementation — and the natural reaction is to weaken the one guard that catches a partial revert of DECISION 1a | P1 | Allowlist is `_flow_table`, `_flow_detail`, `_eligibility_detail`. §12.6 states the rule as *how* the column is read, not *where*: any of the three may count it, none may branch a class on it |
| R4 | **The document's loudest claim was false.** The Status block, §20.2 and the roadmap all said the DECISION 1a sweep had not landed and Stage 5 could not be implemented. All nine items were in the tree; 490 tests green | P1 | Status rewritten; §20 recast as a dated record of work done; the roadmap's paragraph struck |
| R5 | **Two rows of §18 were carried over from the DECISION 1 draft.** The §7.2 verification row recorded the *three-class* table (all 50/43/11), contradicting §7.2's own rendering; the entry-count row said 14 / 13 "before Stage 5" when the measured counts are 12 / 11, having silently included `derive_cohort`'s two | P2 | Both rows re-run and replaced, each carrying a note of what it replaced. DoD-6's `Compare` scan — deleted by §4.3, and citing a 50-record cohort that no longer arises — struck with it |
| R6 | **The tail-slice repair could not be demonstrated.** T6 and DoD-12 said to insert a `build` call into `test_derive.py`'s `_full_pipeline()`. That helper loads the fixture, and the fixture collapses to HUG with no control arm, so **P2 raises** — an exception, never a failing assertion. Inserted above the existing `derive_cohort` line it raises C2 instead | P1 | T6 rewritten around a throwaway `cohort_frame()` pipeline. `_full_pipeline()` is deliberately not touched: six tests assert against its current shape |
| R7 | **`data.py` still explained eligibility in the withdrawn language.** A9's **runtime error message** — shown to the data owner — says a blank "moves a patient from indeterminate to eligible [DECISION 1, §3]"; the `_missingness` comment block says the reason column separates "eligible from indeterminate". Neither was touched by the sweep | P2 | Three `data.py` prose rows added to §10 and to T6. The per-centre-columns argument survives; only the class does not |
| R8 | **Nine smaller defects, each of which stops or misleads an implementer.** §12.9's driver called a `data.run_frame` that does not exist; §5.1 left `_assert_cohort` unwritten while `build` and `_assert_cohort_inputs` were written out; §12.0/§12.4 needed the normalised `hand_frame()` and §12.3 the bare one, with nothing saying which; C4's cases had no construction note though A3 and E4 both raise earlier; §12.11's raise is `KeyError('onset_type')`, not `SchemaError`; the flow table's `centres` cell rule was unstated; §12.8's `status` partition fails on the `all` row; `cohort_frame()`'s import list omitted `config` and `pd`; §9 referred to a three-level column that no longer exists | P2 | All nine corrected in place. §12.0.1 is new and is the bare-versus-normalised rule |

**Two things the review checked and deliberately did not change.** §6.2's position for the `cohort`
kind is a judgement about a rendered document and stands. §7.4's cohort-wise missingness entry was
put to the PI as the one strikeable item and **kept**: the stage that creates the population is the
honest place to state what is missing for it, and the two existing tables are over a population no
estimate uses.

**One thing the review confirmed rather than found.** The performance section is a non-question at
this size — two boolean masks over ≤126 rows, ~31 columns × 5 reductions for the missingness table,
one `df.copy()`. §15's "a cache is a second source of truth" is right and §2's cost statement stands.
