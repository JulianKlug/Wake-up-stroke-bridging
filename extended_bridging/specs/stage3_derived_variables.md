# Stage 3 spec — derived variables

Implements roadmap Stage 3 [§5, §6, §13]. Section references in brackets are to
`statistical_analysis_plan.md`. Numbers and decisions referenced as DECISION *n* are established in
Stage 0 and recorded in `../out/stage0_data_inventory.md`. Constants named in `SMALL_CAPS` are Stage
1's and live in `config.py`; `stage1_config_and_data_contract.md` is their specification. The frame
this stage receives is specified by `stage2_load_and_clean.md` §11.

**Status.** Written 2026-08-10 against the frame as verified in §18, and **revised the same day after
an engineering review** that ran the pandas semantics and read the existing tests rather than taking
either on trust; §20 is the ledger of what that changed and why. It is the **sole source for the
Stage 3 implementation**: everything the implementer needs is here, and anything not here is not to
be invented.

Four claims in the first draft were false against this repository — one of them load-bearing for an
acceptance criterion that could not have failed — and one required test edit was missing from §10.
They are corrected in place, not appended, so the body reads as one document; §20 exists so a reader
comparing against the roadmap or an earlier draft can see which sentences moved and why.

**One [§13] amendment lands with this spec.** The target-mismatch subgroup is withdrawn (§6.1). The
amendment itself is in the SAP, dated and reasoned there; §6.1 records what it removes from the code
and what it deliberately leaves alone.

**Nothing about Stage 3 is settled anywhere but here.** If a decision was made about this stage and
is not in this file, it is not a decision.

**Goal.** Every variable the analysis needs that the workbook does not contain, built from the
registries rather than from literals, with missingness reimposed everywhere a comparison would
otherwise manufacture a zero.

**Not in scope.** Classifying eligibility (Stage 4), restricting the cohort (Stage 5), factor coding
and design matrices (Stage 6). Stage 3 **adds columns**. It drops no row, drops no column, and edits
no value that Stage 2 delivered.

---

## 0. Where Stage 3 sits

```
  data.load()  →  (df, audit)          25 analysis columns, 126 rows   [Stage 2 §11]
                       │
                       ▼
  ┌────────────────────────────────────────────────────────────────────────┐
  │  STAGE 3 — derive.py                     reads config.py, data.Audit   │
  │                                                                        │
  │   derive(df, audit)              row-wise: each record's own values    │
  │     ├─ _assert_onset_flags(df)   both flags positive → SchemaError     │
  │     ├─ onset_type                three levels, <NA> if a flag is <NA>  │
  │     ├─ unknown_onset             [§13], from onset_type                │
  │     ├─ _dichotomies(df)          four, through OUTCOMES and OPS        │
  │     └─ absence_by_column(…)      one row per derived column            │
  └────────────────────────────────────────────────────────────────────────┘
                       │
                       ▼
              31 columns  →  Stage 4 eligibility  →  Stage 5 cohort restrictions
                                                              │
                                                              ▼
  ┌────────────────────────────────────────────────────────────────────────┐
  │  STAGE 3 — derive.py, second entry point                               │
  │                                                                        │
  │   derive_cohort(df, audit)       cohort-wise: needs the whole frame    │
  │     ├─ core_above_median         [§13]; median OF THIS FRAME, frozen   │
  │     └─ constant_covariates(…)    detected and logged, never dropped    │
  └────────────────────────────────────────────────────────────────────────┘
                       │
                       ▼
              33 columns  →  Stage 6 propensity  →  Stage 10 bootstrap
                                                    resamples this column,
                                                    never recomputes it
```

Both entry points append to the `Audit` that `load()` created. Neither writes a file.

### 0.1 Why two entry points, and not one

Roadmap Stage 3 lists five things to build, and they are not all the same kind of thing. Three are
properties of a **record** — the onset factor, the four dichotomies, the unknown-onset flag — and can
be computed the moment the frame is read. Two are properties of a **cohort**: a median, and whether a
covariate varies. A cohort does not exist until Stage 5 has applied the [§3] restrictions, and Stage
1 §6 already refuses to freeze the median in the configuration for exactly this reason — it "is a
property of the cohort, computed by Stage 3 after the [§3] restrictions, and freezing it as a number
here would silently decouple it from the cohort it describes."

The numbers say the same thing. §18 measures the core-volume median at **6.0 mL over all 126
records** and **5.0 mL over the 93-record primary cohort**. One function computing both halves at
Stage 3's natural position would have used the first, and the [§13] subgroup would have been a split
of the wrong population — with nothing to notice, because both numbers are plausible and neither
appears anywhere else.

So the stage owns two functions and the pipeline calls them at two different points. Stage 5's spec
will call `derive_cohort`; Stage 3 owns the code and this specification owns the contract.

### 0.2 What Stage 3 does not touch

`derive` and `derive_cohort` both `df.copy()` before writing, as `_pipeline_after_read` does, so a
caller's frame is never mutated underneath it. No existing column is edited, no row is dropped, and
no value is corrected — Stage 2 owns corrections and its audit log is the record of them. If Stage 3
appears to need to fix a value, the fix belongs in Stage 2 §7 with its own audit entry.

The promise is about the two **entry points**, and it is worth being exact about that because one
private breaks it by design: `_dichotomies` writes into the frame it is passed (§5.1). Inside `derive`
that frame is already the copy, so the guarantee holds for every caller of the public surface. A test
reaching past the entry point into `_dichotomies` is outside the promise, and §5.1 states what it owes
instead.

## 1. Deliverables

| Path | Contents |
|---|---|
| `extended_bridging/derive.py` | `derive`, `derive_cohort`, `constant_covariates` |
| `extended_bridging/test_derive.py` | the acceptance tests in §12 |
| `extended_bridging/config.py` | **amended** — `ONSET_TYPE_FROM_FLAG`, `SUBGROUPS`, `COHORT_DEPENDENT_SUBGROUPS`, `DERIVED_DICHOTOMIES`; `DERIVED_NAMES` becomes a computed view; `TARGET_MISMATCH` and `MISMATCH_RATIO_UNDEFINED_COUNTS_AS_MET` are deleted (§6.1) |
| `extended_bridging/data.py` | **amended** — the log spans the pipeline (§9.1), `KINDS` gains `derivation` (§9.2), `_missingness` is generalised into `absence_by_column` (§8.2), and `_correct`'s target-mismatch sentence goes (§6.1) |
| `extended_bridging/test_config.py` | **amended** — the new registries; that no subgroup name reaches a covariate list; and 9.13's resolvable union and 9.8's immutability list, both of which name their constants explicitly (§12.13) |
| `extended_bridging/test_data.py` | **amended** — the renamed path in 12.8, the heading count in 12.16, one comment in `_HAND_RECORDS` (§10) — **and one repair**: the heading-adjacency assertion that a new `kind` breaks (§9.2) |
| `extended_bridging/statistical_analysis_plan.md` | **amended** — §13, the target-mismatch withdrawal. Landed with this spec, not with the implementation |
| `extended_bridging/implementation_roadmap.md` | **amended** — Stage 3's entry. Landed with this spec |
| `specs/stage1_config_and_data_contract.md`, `specs/stage2_load_and_clean.md` | **amended** — the sections the above make false (§10) |

`out/logs/audit_<label>.md` is an **output, not a deliverable**, as in Stage 2 §1: gitignored,
written only when a caller asks, and it names patients. Nothing under `specs/` may quote a case
identifier, and nothing in this file does.

## 2. Environment

Unchanged from Stage 2 §2.1. `uv`, Python 3.12, no new dependency — `derive.py` needs pandas and
nothing else. Commands run from `extended_bridging/`, flat module layout, `import config as C`.

```bash
cd extended_bridging && uv sync && uv run pytest -v
```

**What this costs to run.** Recorded for the same reason Stage 2 §2.3 recorded it. Every derivation
here is one vectorised pass over at most 126 values; `constant_covariates` is thirteen `nunique`
calls. The whole stage is faster than the `read_excel` that precedes it. There is nothing to cache,
and §15 says so as a standing instruction rather than as an observation.

## 3. Module shape

```python
# derive.py
from __future__ import annotations

from collections.abc import Sequence

import pandas as pd

import config as C
from data import Audit, absence_by_column
```

**No new types.** Stage 2 introduced three and §0.1 there defends each one; Stage 3 introduces none,
because everything it produces is a column or an audit entry and both already have types. A
`Subgroup` dataclass was considered and rejected: the two subgroups share no structure — one is a
relabelling of a factor, the other is a comparison against a statistic of the frame — so the type
would carry a name and a docstring and nothing else.

The public surface is three functions:

```python
def derive(df: pd.DataFrame, audit: Audit) -> pd.DataFrame
def derive_cohort(df: pd.DataFrame, audit: Audit) -> pd.DataFrame
def constant_covariates(
    df: pd.DataFrame, covariates: Sequence[str] = C.PS_COVARIATES_FULL) -> tuple[str, ...]
```

and the privates are `_assert_onset_flags`, `_onset_type`, `_unknown_onset`, `_dichotomies`,
`_core_above_median`. `constant_covariates` is public because Stage 6 calls it too, on each bootstrap
replicate's design frame, where the answer is different and is used rather than logged.

`derive.py` is **not** exempt from the Stage 1 §7 raw-name scan and must never become exempt. It
names no raw header; every column it reads it reads by analysis name.

### 3.1 The one pandas fact this stage turns on

Every derivation below is a comparison against a column that may be missing, so the whole stage rests
on what pandas does with `<NA>` in a comparison and in a mask. It is declared once, here, because it
is stated wrongly in the obvious way — "a comparison against a missing value returns false", which is
what roadmap Stage 3 says and what is true of numpy and of most other frameworks. It is **not** true
of the nullable dtypes this pipeline reads into. Verified on pandas 2.3.3, the pinned version:

```
                                              NA behaviour        mask still needed?
  Int64      s <= 2                           → <NA>   propagates      no, but see below
  string     s != "witnessed"                 → <NA>   propagates      no, but see below
  boolean    b & b                            → <NA>   propagates      no, but see below
  float64    s > 5.0        (NaN > 5.0)        → False  ABSORBS         YES
  int64      s > 5.0        (no NA possible)   → False  n/a             n/a
  object     s <= 2                            → False  ABSORBS         YES

  Series.mask(cond) where cond carries <NA>:  <NA> is treated as True — the row is masked
                                              and takes the replacement value.
```

Two consequences run through §4, §5 and §6, and they point in opposite directions:

- **Which columns are nullable is a Stage 1 decision, not a Stage 3 one.** `READ_DTYPES` is built
  from `Column.dtype`, so `mrs_90d`, `wake_up` and `unwitnessed` arrive as `Int64` and propagate,
  while `core_ml` and `tmax6_ml` declare **no** dtype (`config.py`'s contract entries for them pass
  `dtype=None`) and arrive as `float64` from the workbook — or as `int64` from a source whose values
  happen to all be integers, which is what `tests/fixture_schema.xlsx` does. So the dtype of the same
  analysis column differs by source, and a mask that is redundant on one frame is load-bearing on
  another.
- **Therefore every mask stays, and none of them is justified by "the comparison would return
  false".** They are justified by making the result independent of a dtype no derivation controls.
  Each site below says which of the two it is today, because a reader who is told a redundant call is
  load-bearing stops believing the ones that are.

## 4. `onset_type` [§5, §6]

### 4.1 The assertion comes first

Roadmap Stage 3: "Assert the two are never both positive." Before any level is assigned, because a
record positive on both belongs to none of the three levels and every way of assigning it one is a
clinical statement nobody has made.

```python
def _assert_onset_flags(df: pd.DataFrame) -> None:
    both = ((df["wake_up"] == 1) & (df["unwitnessed"] == 1)).fillna(False)
    if int(both.sum()):
        raise C.SchemaError(
            f"wake_up and unwitnessed are both 1 on {int(both.sum())} record(s): "
            f"{', '.join(sorted(df.loc[both, 'case_id']))}. [§5]'s onset factor has three levels "
            "and a record positive on both belongs to none of them. Correct the workbook or amend "
            "[§5]; do not choose a level here.")
```

**Why it raises rather than resolving.** A resolution — prefer wake-up, prefer unwitnessed, call it
witnessed — would apply silently to every future workbook, and `onset_type` is a [§6] covariate that
enters the propensity model, the outcome regression and the balance table. Stage 2 took the same
posture toward the 999 groin time (§7.3 there) and DECISION 2 toward the contradictory death flags: a
derivation rule where one is available, a loud stop where none is.

**`.fillna(False)` here is belt-and-braces, and the honest reason is not the obvious one.** Per §3.1
`(df[flag] == 1)` on an `Int64` column yields `<NA>` where the flag is missing, so `both` carries
`<NA>` rather than `False` — but `.sum()` on a `BooleanDtype` series skips `<NA>`, so
`int(both.sum())` is already the count of records positive on *both* and the `.fillna(False)` changes
no count today. It stays for one reason: it makes the mask total, so the expression means the same
thing if `wake_up` ever loses its `Int64` declaration and starts absorbing (§3.1). A record with one
missing flag is not evidence that both are positive, and it is handled in §4.3 rather than here.

What must **not** be written here is the plausible-sounding version of that argument — that a
`BooleanDtype` mask carrying `<NA>` "is not a mask". It is one, and §3.1 records what it actually
does: `<NA>` masks, taking the replacement. That misreading matters in §4.2 rather than here, and
§4.3 records which way each failure falls.

v7 has **0** such records (§18), so this branch is unreachable today. It is written for the workbook
that has not arrived yet, which is the same reason Stage 2 keeps A1.

### 4.2 The build, with the levels declared rather than written

`derive.py` may not write `"witnessed"`, `"unwitnessed"` or `"wake_up"`. `FACTOR_LEVELS` already
declares the level set and `REFERENCE_LEVELS` already declares which one is the baseline; a literal
in `derive.py` would give each level a second declaration, and the two would be free to drift on the
one edit — a renamed level — that nothing else would catch. One new constant closes the gap:

```python
# config.py — which flag produces which non-baseline level. The baseline is
# REFERENCE_LEVELS["onset_type"], which is already declared and is Stage 6's reference category, so
# it is deliberately not repeated here. test_config.py asserts that this tuple's levels plus that
# baseline are exactly FACTOR_LEVELS["onset_type"], so a level cannot be added in one place only.
ONSET_TYPE_FROM_FLAG: Final[tuple[tuple[str, str], ...]] = (
    ("unwitnessed", "unwitnessed"),
    ("wake_up", "wake_up"),
)
```

```python
def _onset_type(df: pd.DataFrame) -> pd.Series:
    onset = pd.Series(C.REFERENCE_LEVELS["onset_type"], index=df.index, dtype="string")
    unknown = pd.Series(False, index=df.index)
    for flag, level in C.ONSET_TYPE_FROM_FLAG:
        onset = onset.mask((df[flag] == 1).fillna(False), level)
        unknown = unknown | df[flag].isna()
    return onset.mask(unknown)
```

Three properties, each load-bearing:

- **The baseline is the declared reference level.** `witnessed` is the clinical baseline [§13] and
  Stage 6's reference category; it appears once, in `REFERENCE_LEVELS`.
- **The loop's order is not a precedence rule.** Both flags positive is impossible by §4.1. If the
  assertion were ever deleted, the last entry would win silently — which is precisely why the
  assertion is not optional, and why this sentence is here rather than a comment reading "wake-up
  wins".
- **`unknown` accumulates over the same tuple.** A third onset flag added to
  `ONSET_TYPE_FROM_FLAG` contributes both its level and its missingness in one edit. A separately
  written `df["wake_up"].isna() | df["unwitnessed"].isna()` would contribute only the level.
- **The two masks are independent guards and each one alone is wrong differently.** Verified, on a
  record whose `wake_up` is `<NA>`:

  ```
    both masks               → <NA>            correct
    fillna, no unknown mask  → "witnessed"     the baseline absorption §4.3 is about
    unknown mask, no fillna  → <NA>            correct by accident: the row is masked to the LEVEL
                                               by the <NA> cond (§3.1), then masked again to <NA>
    neither                  → "wake_up"       a fabricated *non-baseline* level
  ```

  So deleting `.fillna(False)` alone is invisible, and deleting it together with the `unknown` mask —
  the plausible "tidy the masks" edit — produces a fabricated `wake_up`, which is worse than the
  fabricated `witnessed` §4.3 warns about and reads more plausibly in a table. §12.3 pins all four
  cells of that square rather than only the first.

### 4.3 A missing flag yields `<NA>`, never the baseline

This is §5.2's trap in a different column, and it is the more dangerous of the two, because
`onset_type` is a [§6] **covariate** rather than an outcome. A record silently labelled `witnessed`
enters the propensity model, the outcome regression and the balance table carrying a fabricated
value. `<NA>` drops it by complete-case [§11], which is what [§11] prescribes and what the reported
denominator will then say.

Nothing upstream forbids a missing flag: Stage 2's A5 asserts that `BINARY_COLUMNS` take values in
`{0, 1}` **or missing**, and `wake_up` and `unwitnessed` are in that list. Stage 0 records both as
non-missing on all 126 records, which means **no test that runs against the real workbook can reach
this branch** — an implementation that omitted the mask would be green on v7 and wrong on v8. §12.3
reaches it on a hand-built frame.

The pilot gets this wrong, and it is worth naming because it is the obvious implementation:
`pilots/data.py:141` builds the factor with `np.select([...], [...], "witnessed")`, whose default
absorbs a missing flag into the baseline. §16 records it as **do not lift**.

### 4.4 `onset_type` is a `string`, not a `Categorical`

Stage 2 §6's argument for `center`, unchanged: Stage 5 restricts the cohort, and a categorical that
keeps a dead level makes every subsequent `groupby(observed=False)` resurrect it as an all-missing
row — in the balance table, the within-centre overlap table [§9] and every subgroup table. Stage 1
§5.4 assigns categorical construction to Stage 6, which builds it from `FACTOR_LEVELS` at the point
of use where the level set is chosen deliberately. `pilots/data.py:141` wraps the factor in
`pd.Categorical` at derivation time; that half is not lifted either.

## 5. The outcome dichotomies [§5]

### 5.1 One loop, driven by the registry

Stage 1 §3.2 already pins this and it is reproduced rather than restated, because a second wording is
a second specification:

```python
def _dichotomies(df: pd.DataFrame) -> pd.DataFrame:
    for key, o in C.OUTCOMES.items():
        if o.source is None:
            continue                                    # read directly; nothing to derive
        src = df[o.source]
        df[key] = C.OPS[o.op](src, o.threshold).astype("Int64").mask(src.isna())
    return df
```

Four columns result — `mrs_0_2_90d`, `mrs_0_1_90d`, `death_90d`, `mrs_5_6_90d` — all from `mrs_90d`.
`tici_2b_3`, `sich` and `ph2` are read directly from the workbook and are already columns of the
frame; the `continue` is what leaves them alone, and `Outcome.source is None` is the only test of
which is which.

**`_dichotomies` writes into the frame it is given, and that is deliberate.** It is the one function
in this module that does. `derive` has already called `df.copy()` (§0.2), so the frame it hands over
is `derive`'s own and an in-place write is free; three privates each copying a 126-row frame to
preserve an immutability nobody needs would be ceremony. The consequence is a contract on the caller
rather than on the function, so it is written down instead of inferred: **a test calling
`_dichotomies` directly passes a frame it owns and must not assert anything about that frame
afterwards.** §12.5, §12.6 and §12.7 all call it directly — they have to, to run under a patched
registry — and each builds its own frame for the purpose.

**No threshold, operator or outcome name is written in `derive.py`.** Stage 1 §3.2 gives the reason
in full: a threshold typed twice can have the manuscript print `mRS 0-1` while the analysis computes
`mRS 0-2`, with no test failing. §12.5 tests the property directly, by running the loop under a
registry whose `mrs_0_2_90d` threshold is 3 and asserting the output moves.

The iteration order is `OUTCOMES`' insertion order, which is [§5]'s table order. Any order would
compute the same four columns; a *declared* order is what makes §9.3's table byte-identical.

### 5.2 `.mask(src.isna())` is the whole of it — but not for the reason it looks like

Roadmap invariant 6 and roadmap Stage 3's acceptance criterion both reduce to this one call. The
roadmap states the reason as "a comparison against a missing value returns false in most frameworks
and would silently manufacture zeros". That is the right *failure* and the wrong *mechanism* for this
frame, and the difference decides whether the acceptance test can ever fail.

`mrs_90d` is declared `Int64` by the contract, so per §3.1 `mrs_90d <= 2` yields `<NA>` on the two
records with no 90-day mRS, `.astype("Int64")` carries that `<NA>` through, and **`.mask(src.isna())`
is a no-op on every frame this repository builds today** — the workbook's, the fixture's, and
`hand_frame()`'s, which casts to the declared dtypes (§12). The call is not decoration and is not
removed; what it buys is that invariant 6 holds for a dtype no derivation controls:

```
  source dtype   without the mask        with the mask
  Int64          <NA>   (correct)        <NA>
  float64        0      (FABRICATED)     <NA>
  object         0      (FABRICATED)     <NA>
```

The fabricated case is not hypothetical. It is what any frame built without `READ_DTYPES` delivers —
a `read_excel` whose `dtype=` argument was dropped, a frame assembled in a notebook, a test fixture
built from a dict and not cast. In that frame the two records with no 90-day mRS become **non-events
in every one of the four**: not-good-outcome, not-poor-outcome, and alive. Two things break at once
and [§11] names both — the denominator moves from 124 to 126, and a patient whose outcome is unknown
is counted as having had a specific one.

It is stated this plainly because the failure is invisible: the column is full, the counts are
plausible, and every downstream table reconciles.

**What this costs the acceptance criteria, and why it is not a licence to skip them.** A test that
only exercises an `Int64` source passes with the mask deleted, so it cannot be *seen* to fail — and a
guarantee whose test has never been watched failing is not evidence (Definition of done 5). §12.7
therefore runs the parametrised missingness assertion over **both** dtypes: the declared `Int64`
source, and the same frame with `mrs_90d` cast to `float64`. The second is the one that fails when
the mask goes.

### 5.3 The shipped dichotomies are never read, not even to compare

DECISION 2. The contract has already dropped all four shipped columns and Stage 0 compared them once
— `mRS56at90days` disagrees with the ordinal source on 3 records and `Deathat90days` on 2, both
recorded there by identifier as standing queries. Rebuilding the comparison here would require naming
a raw header, which Stage 1 §7 forbids and its 9.4 scan would flag; it would also reopen a question
DECISION 2 closed. If the data owner resolves the queries, the answer arrives as a new workbook and a
new `DATA_SHA256`, not as a reconciliation in `derive.py`.

### 5.4 `Int64`, not `boolean`

Every consumer treats these as 0/1 numerics — weighted means, risk differences, the Firth outcome
model, the augmented estimator's `Y`. `Int64` keeps `<NA>` and arithmetic; `boolean` keeps `<NA>` and
loses the arithmetic. The derived columns are deliberately **not** added to `BINARY_COLUMNS`: that
tuple is Stage 2's A5 domain check over the frame as read, and adding names to it would have Stage 2
assert columns that do not exist when it runs. §12.5 asserts the domain here instead.

## 6. Subgroups [§13]

### 6.1 Two, not three — the target-mismatch withdrawal

[§13] listed three subgroups. The amendment of 2026-08-10 withdraws target mismatch; the SAP carries
the decision and its reasoning, and this section records what it means for the code.

The short form of the argument: the cohort is CTP-selected by construction (§2), so target mismatch
is close to the criterion that *admitted* these patients rather than a contrast within them. Stage 0
records core volume as exactly 0 in 50 of 125 records, so most patients clear the core and ratio
criteria and the flag turns almost entirely on one volume threshold. A near-constant subgroup does
not generate hypotheses; it produces a table row that reads as a finding. The count is marginal and
predates any outcome examined by arm, as does the decision.

**What the withdrawal removes.**

| Removed | Because |
|---|---|
| `TARGET_MISMATCH` | its only consumer was the subgroup |
| `MISMATCH_RATIO_UNDEFINED_COUNTS_AS_MET` | declared solely to fix the `core == 0` convention inside that criterion |
| the `crossing` computation and its sentence in `data.py`'s `_correct` | it counts records whose recomputed penumbra crosses `TARGET_MISMATCH["min_penumbra_ml"]`, and reports that "[§13] target-mismatch membership changes for them" — a consequence that no longer exists |
| Stage 1 spec §6's second half, Stage 2 spec §7.2's consequence paragraph and §18's criteria row | they state that consequence |

**What it does not touch.** `penumbra_ml` keeps its [§6] exclusion from every model, keeps its
`BALANCE_ONLY` role, and keeps Stage 2 §7.2's recomputation — that correction rests on [§6] defining
penumbra as `tmax6_ml - core_ml`, which is independent of any subgroup. The per-case table in the
`penumbra_recomputed` entry already carries stored, recomputed and difference, so removing the
sentence costs the log nothing and costs no test (§10).

**A question this retires rather than answers.** The criterion needed three comparisons and no source
in the repository pinned whether each was strict or inclusive: DEFUSE-3 uses core **<** 70 mL, ratio
**≥** 1.8, mismatch volume **≥** 15 mL, while Stage 2 §7.2's prose wrote all three strictly. The
constants' names (`min_ratio`, `min_penumbra_ml`) imply inclusive. Nothing decided between them, and
the difference is real — it moves any patient sitting exactly on a boundary. Recorded here so that a
future amendment restoring the subgroup answers the question *before* writing the code, rather than
rediscovering it.

### 6.2 `unknown_onset`

```python
def _unknown_onset(onset: pd.Series) -> pd.Series:
    return (onset != C.REFERENCE_LEVELS["onset_type"]).astype("Int64").mask(onset.isna())
```

[§13] asks for "unknown versus witnessed onset". Witnessed is the baseline level, so unknown is its
complement over the three levels.

**Built from `onset_type`, never from the two flags again.** Two derivations of one concept drift,
and this pair would drift silently: both are correct on v7, and the difference appears only where a
flag is missing — where the flag route yields 0 and the `onset_type` route yields `<NA>`.

**`<NA>` propagates, and is never 1.** "We do not know the onset type" is not "the onset was
unknown". The distinction would be easy to lose in the count: §18 records 93 of 126 records as
unwitnessed or wake-up already, so one more would not look wrong.

### 6.3 `core_above_median`

```python
def _core_above_median(df: pd.DataFrame) -> tuple[pd.Series, float]:
    core = df["core_ml"]
    median = float(core.median())                       # skipna, the pandas default
    return (core > median).astype("Int64").mask(core.isna()), median
```

Three rules, each with a number behind it from §18.

- **The median is of the frame this is given**, which is why the function lives in `derive_cohort`.
  §18 measures 6.0 mL over all 126 records and 5.0 mL over the 93-record primary cohort. Computing it
  at Stage 3's first entry point would split the cohort on the wrong number, and both numbers are
  plausible.
- **Ties sit below.** The comparison is `>`, so a patient whose core volume equals the median is not
  above it. Not a formality: §18 records **3** records at exactly 5.0 mL in the primary cohort, and
  39 with a core of exactly 0. `>` and `>=` differ on those 3 patients, and nothing else in the
  repository would catch which was used.
- **A missing core yields `<NA>`, and this is the one mask in Stage 3 that is load-bearing today.**
  `core_ml` declares no dtype in the contract, so per §3.1 it arrives as `float64` from the workbook
  and `NaN > 5.0` is `False`, not `<NA>` — the comparison **absorbs**. Delete `.mask(core.isna())` and
  the one cohort record with no core volume is silently assigned to the below-median group, in the
  arm where it happens to sit, with no `<NA>` anywhere to notice. Contrast §5.2, where the same call
  on an `Int64` source is redundant today: this is the asymmetry §3.1 exists to keep straight, and it
  is why §12.10 asserts the missing-core cell directly rather than by analogy with the dichotomies.
  Complete-case [§11] then drops the record from the subgroup analysis, which is correct.

The audit entry prints the median, the frame's row count and all four counts — above, at, below,
missing — so the split's balance is legible without re-deriving it. §18's split is 44 above, 48
at-or-below, 1 missing.

**An all-missing `core_ml` degrades quietly, and the log is what catches it.** `Series.median()` on a
column with no non-missing values returns `NaN`, `core > NaN` is `False` everywhere, and the mask then
sets every row to `<NA>` — so the subgroup is empty rather than wrong, which is the right failure. It
is not silent: `data._fmt` renders the median as the literal `missing` in the audit entry, so the
`median (mL)` row reads `missing` and the `above`/`below` counts read 0. No guard is added for it —
Stage 5's cohort cannot be all-missing on a [§6] covariate without A6 and the balance table having
said so first — but §13 records it so it is not later diagnosed as a derivation bug.

### 6.4 The median is computed once and frozen

Decision of 2026-08-10, and the reason Stage 10 must be told about it here rather than discovering it
there.

Stage 10 refits everything in 2000 replicates, and the natural shape of a replicate is "resample,
then re-run the pipeline". Under that shape a `derive_cohort` that simply recomputed would give every
replicate its own cut-point and therefore its own subgroup: membership would move between replicates,
the estimand would differ between them, and **no test would fail** — the estimates would remain
plausible and the interval would simply be of a slightly different quantity than the point estimate.

So the guard is structural, not documentary:

```python
already = [c for c in C.COHORT_DEPENDENT_SUBGROUPS if c in df.columns]
if already:
    raise C.SchemaError(
        f"{', '.join(already)} already present. The [§13] median subgroup is computed once, on the "
        "analysis cohort, and frozen [Stage 3 §6.4]: a bootstrap replicate resamples this column, "
        "it does not recompute it. Calling derive_cohort twice would give the replicate its own "
        "cut-point and its own subgroup, silently.")
```

with the guard ranging over a declared constant rather than a literal:

```python
# config.py
COHORT_DEPENDENT_SUBGROUPS: Final[tuple[str, ...]] = ("core_above_median",)
```

**What Stage 10 must do instead:** resample the frame that already carries the column. §11 states it
as part of the handover, and Stage 10's own spec will carry the test.

### 6.5 The subgroup registry

```python
# config.py — the [§13] subgroups, after the 2026-08-10 amendment withdrew target mismatch.
# Neither is a covariate: test_config.py asserts no subgroup name appears in any covariate list, so
# a subgroup cannot become an adjustment variable by being convenient.
SUBGROUPS: Final[dict[str, str]] = {
    "unknown_onset":     "unwitnessed or wake-up onset versus witnessed [§13]",
    "core_above_median": "core volume above the cohort median [§13]",
}
```

Neither is post-time-zero: both are baseline attributes, so neither joins `POST_TIME_ZERO`. Both are
`Int64` 0/1 with `<NA>` where an input is missing, and Stage 11 will report each subgroup's cell
sizes beside its estimate — a subgroup whose smaller cell cannot support an interaction test is
reported as not estimable rather than estimated, which is the guard §6.1's argument implies but which
belongs to Stage 11.

## 7. Zero-variance covariates

### 7.1 Detected and logged here, dropped in Stage 6

```python
def constant_covariates(
        df: pd.DataFrame, covariates: Sequence[str] = C.PS_COVARIATES_FULL) -> tuple[str, ...]:
    """Names in `covariates` with at most one distinct non-missing value in `df`, in that order."""
    return tuple(c for c in covariates if df[c].nunique(dropna=True) <= 1)
```

`nunique(dropna=True) <= 1` catches three cases in one expression: a constant column, an all-missing
column (0 distinct), and a factor left with a single level after a restriction. It works unchanged on
`onset_type` and `center`, which are strings, and on the numeric covariates.

It ranges over the caller's sequence, so the order of the result is the caller's and the log is
deterministic — never over a set, per Stage 2 §9.3.

**Its default argument carries a precondition, and the docstring must say so.** `PS_COVARIATES_FULL`
contains `onset_type`, which the contract cannot supply, so `constant_covariates(df)` raises `KeyError`
on any frame `derive` has not run on. That is the correct failure — a caller asking which covariates
are constant before the covariates exist has a sequencing bug, and a silent skip of the missing name
would answer a question about a different covariate set. `derive_cohort` is downstream of `derive` by
construction, and Stage 6 passes its own design-frame column list, so neither reaches it; the sentence
is for the third caller.

**It returns names and changes nothing.** Not the frame, and not `covariates`. §12.8 asserts
`C.OUTCOME_COVARIATES is C.PS_COVARIATES` still holds after a call, because invariant 3 rests on
those two resolving to one object and a function that filtered a list in place would break it in a
way that only shows up as a covariate quietly missing from one of the two nuisance models.

**Why logged and not raised.** A constant covariate is legitimate. [§6] anticipates exactly one —
`first_image_mri`, "dropped automatically if constant" — and the contract has already dropped it by
kind, so it cannot reach this function at all. What is not legitimate is a *silent* drop, and the
division of labour follows from where each is knowable: Stage 3 reports what is constant in the
analysis cohort, Stage 6's design-matrix builder drops constant columns at the point of use, which is
also the only place that can see a bootstrap replicate having emptied a factor level.

### 7.2 Nothing mutates a covariate list

Roadmap Stage 3 says "drop any covariate with zero variance", and the roadmap's own Stage 6 says the
design-matrix builder drops constant columns. Two stages cannot both own the deletion. The reading
this spec takes — Stage 3 detects, Stage 6 drops — is the one that leaves invariant 3 intact:

> The propensity and outcome-regression covariate lists resolve to the same configuration entry […]
> Assert the override registry against that explicit list.

`OUTCOME_COVARIATES` **is** `PS_COVARIATES`, by object identity, and `STANDARDISATION_COVARIATES` and
`PS_COVARIATES_FULL` are computed from it. A Stage 3 that removed a name from the frame would leave
those lists naming a column that no longer exists; a Stage 3 that removed a name from the *lists*
would break the identity that invariant 3 is asserted on. Detecting and logging breaks neither.

## 8. Missingness of the derived columns [§11]

### 8.1 The requirement

Roadmap Stage 3: "Mark structurally non-applicable fields as distinct from missing ones [§11]." Stage
2 §10 built that classification — `structural`, `not recorded`, `missing`, `complete` — over the 25
analysis columns. The derived columns are not in it, because they did not exist when it ran. Stage 3
extends it to them.

### 8.2 One classification, called twice

Stage 2's `_missingness` is generalised rather than copied:

```python
# data.py, public because Stage 3 calls it
def absence_by_column(
        df: pd.DataFrame, audit: Audit, columns: Sequence[str], step: str, detail: str) -> None
```

Stage 2 keeps its structural loop over `STRUCTURALLY_NON_APPLICABLE` — one entry per declared column,
belonging to the stage that owns the declaration — and then calls
`absence_by_column(df, audit, sorted(C.ANALYSIS_NAMES), "absence_by_column", …)`. Stage 3 calls
`absence_by_column(df, audit, sorted(C.ROW_WISE_DERIVED), "absence_by_derived_column", …)`.

**Why generalise rather than write a second table.** The four-way `kind` classification is a rule
about columns, not about a stage. A second implementation would drift on the branch that matters —
`structural` versus `missing` — which is the exact misreading Stage 2 §10.2 exists to prevent. The
per-centre columns, the `CENTER_ORDER` header and the reconciliation that §12.11 there asserts all
come along unchanged.

### 8.3 No derived column is structural, and none is informative absence

Every absence a derived column carries is **inherited**: `onset_type` and `unknown_onset` from the
two onset flags, the four dichotomies from `mrs_90d`, `core_above_median` from `core_ml`. So every
row of Stage 3's table is `missing` or `complete`, and §12.9 asserts it.

That assertion is the point of the section. `structural` means "the value does not exist for this
patient", and it is the one label that tells a reader an absence is not data loss. A future derived
column given that label to make a table look better would be misrepresenting exactly what Stage 2
§10.2 built the distinction to represent.

`core_above_median` is absent from `derive`'s table, because the column does not exist at that point.
It is accounted for in `derive_cohort`'s own entry, which carries its missing count (§9.3). The two
row counts therefore differ by one **by design**, and §12.9 asserts the table's rows are exactly
`ROW_WISE_DERIVED` rather than `DERIVED_NAMES`.

### 8.4 The derived-name registries

All computed, none written twice:

```python
# config.py
DERIVED_DICHOTOMIES: Final[tuple[str, ...]] = tuple(
    key for key, o in OUTCOMES.items() if o.source is not None)
DERIVED_NAMES: Final[tuple[str, ...]] = ("onset_type", *DERIVED_DICHOTOMIES, *SUBGROUPS)
ROW_WISE_DERIVED: Final[tuple[str, ...]] = tuple(
    name for name in DERIVED_NAMES if name not in COHORT_DEPENDENT_SUBGROUPS)
```

`DERIVED_NAMES` today is the one-element literal `("onset_type",)`, written when only Stage 1
existed. `"onset_type"` remains the single literal — it is the only derived name that no other
registry declares — and §12.13 asserts `DERIVED_NAMES` is disjoint from `ANALYSIS_NAMES`, so a
derived column can never shadow one the contract delivers.

**Where each constant goes in `config.py`, because these are computed and the file executes top to
bottom.** `DERIVED_NAMES` currently sits in the covariates block, *above* `OUTCOMES`. Written there,
the version above is a `NameError` on import. The placement is part of the specification, exactly as
`POST_TIME_ZERO`'s already is ("Defined after `OUTCOMES`; the ordering is load-bearing"):

```
  config.py, in file order
  ├── covariates block
  │     PS_COVARIATES … PS_COVARIATES_FULL
  │     DERIVED_NAMES              ← the literal here is REPLACED BY A POINTER COMMENT
  │                                  ("computed in the derived-names block below, after the
  │                                   registries it reads"), so a reader who knows the old
  │                                   position is not left thinking the constant was deleted
  ├── OUTCOMES  [§5]
  │     DERIVED_DICHOTOMIES        ← immediately after OUTCOMES, which it reads
  ├── OUTCOME_MODEL_OVERRIDES, outcome_model_covariates, ELIGIBILITY_ORDER
  ├── POST_TIME_ZERO               ← unchanged; already declared after OUTCOMES
  ├── STRUCTURALLY_NON_APPLICABLE … BINARY_COLUMNS, SEX_LABELS
  ├── subgroups block  [§13]       ← TARGET_MISMATCH and
  │     SUBGROUPS                     MISMATCH_RATIO_UNDEFINED_COUNTS_AS_MET deleted from
  │     COHORT_DEPENDENT_SUBGROUPS     here (§6.1); MRS_THRESHOLDS stays
  └── derived-names block          ← NEW, after both OUTCOMES and SUBGROUPS
        DERIVED_NAMES, ROW_WISE_DERIVED
```

`ONSET_TYPE_FROM_FLAG` (§4.2) reads nothing, so it goes with the factor declarations beside
`FACTOR_LEVELS` and `REFERENCE_LEVELS`, which are what §12.13 asserts it against. No test is added for
the ordering: a mis-ordered module fails at import, on every test, with the offending name in the
traceback — the loudest failure available and one no assertion improves on.

## 9. The audit log becomes the pipeline's

### 9.1 The rename

Stage 2 §11 already states that "Stage 3 appends; it does not re-open" the same `Audit`. The artefact
is the only thing that still says otherwise:

| | Was | Is |
|---|---|---|
| path | `out/logs/stage2_audit_<label>.md` | `out/logs/audit_<label>.md` |
| title | `# Stage 2 — data audit` | `# Analysis audit` |
| footer | `Generated by extended_bridging/data.py.` | `Generated by the extended_bridging pipeline.` |

`_audit_path` keeps its shape and its label, so Stage 2 §9.5's guarantee holds unchanged: the
fixture's log can still never overwrite the workbook's. Everything else about the header is unchanged
— in particular its numbers are still read back from the entries, never recomputed.

The cost is three string literals in `data.py` and the clauses in acceptance 12.8 that name the path.

### 9.2 The new kind

```python
KINDS: Final[tuple[str, ...]] = (
    "provenance", "contract", "correction", "observation",
    "derivation",                                              # ← Stage 3
    "structural", "missingness")

_HEADINGS[...] = {..., "derivation": "Derivations", ...}
```

Inserted **before** `structural`, so the document still reads in pipeline order: what was read, what
the contract did, what was corrected, what was observed, what was derived, what is structurally
absent, what is missing. Appending it at the end would print the derivations after the missingness
table that describes them.

**This insertion breaks one existing Stage 2 test, and it is the only test in the repository that any
new `kind` breaks.** `test_data.py`'s `test_all_six_headings_appear_and_an_empty_section_says_none`
takes the rendered document between the `observation` heading and the `structural` heading and asserts
it strips to exactly `_none_`:

```python
    observations = rendered.split(f"## {data._HEADINGS['observation']}")[1]
    assert observations.split(f"## {data._HEADINGS['structural']}")[0].strip() == "_none_"
```

With `derivation` between them that slice is `_none_\n\n## Derivations\n\n_none_` and the assertion
fails. Every other heading test iterates `data.KINDS` and adapts by itself. The repair is not to
re-point the second split at `derivation` — that hardcodes the *new* neighbour and rots on the next
kind inserted, which is the same defect one position along. It is to take the neighbour **from
`KINDS`**: find `observation`'s index, split on the heading of the kind that follows it, and rename the
test off the word "six". §10 carries it as an amendment and T3 verifies it, so it is not discovered at
T8.

`_MUST_NAME_CASES` is unchanged. A derivation describes a column, not a patient. The one place Stage
3 names patients is §4.1's assertion, which raises — a frame with both onset flags positive never
reaches a log.

### 9.3 The entry inventory

Six entries, in this order. `step` strings are exact. This is the Stage 2 §9.6 discipline applied
here for the same reason: byte-identical reproduction holds for any log, including one whose contents
an implementer chose freely, so the contents are pinned.

| kind | step | recorded by | `n` counts | table | `case_ids` |
|---|---|---|---|---|---|
| `derivation` | `onset_type` | `derive` | records assigned a level | 4 rows + header | — |
| `derivation` | `unknown_onset` | `derive` | records flagged 1 | — | — |
| `derivation` | `dichotomies` | `derive` | dichotomies derived | 4 rows + header | — |
| `missingness` | `absence_by_derived_column` | `derive` | rows in the frame | one row per `ROW_WISE_DERIVED` + header | — |
| `derivation` | `core_above_median` | `derive_cohort` | records above the median | 5 rows + header | — |
| `derivation` | `constant_covariates` | `derive_cohort` | covariates found | `n` rows + header, only when `n > 0` | — |

`detail` templates, every number through `data._fmt`:

- `onset_type` — `three levels from wake_up and unwitnessed [§5]; the baseline is
  REFERENCE_LEVELS['onset_type'] = {ref!r}. {n_missing} record(s) carry <NA> because a source flag is
  missing — never the baseline [§4.3].` Table `level | n`, rows in `FACTOR_LEVELS["onset_type"]`
  order, then a final `<NA>` row — so four rows always, whatever the counts.
- `unknown_onset` — `[§13] subgroup: onset_type != {ref!r}. {n} of {n_rows} record(s); {n_missing}
  carry <NA>, inherited from onset_type and never read as 0.`
- `dichotomies` — `{k} dichotomies derived from their ordinal source [§5], each carrying exactly its
  source's missingness. No threshold is written outside OUTCOMES.` Table
  `outcome | rule | n | events | missing`, one row per `DERIVED_DICHOTOMIES` in registry order, with
  `rule` taken from `Outcome.rule` so the printed derivation and the applied one are the same fields.
- `absence_by_derived_column` — `one row per derived column; kind is missing or complete. No derived
  column is structurally non-applicable: every absence here is inherited from a source column [§8.3].`
- `core_above_median` — `[§13] subgroup: core_ml > the median of the frame this was computed on. The
  median is a property of the cohort [Stage 1 §6] and is frozen here — Stage 10 resamples this
  column, it does not recompute it [§6.4]. Ties sit below.` Table `statistic | value` with rows
  `median (mL)`, `records`, `above`, `at the median`, `below`, `missing`.
- `constant_covariates` — `{n} of {k} covariate(s) in PS_COVARIATES_FULL have at most one distinct
  non-missing value in this frame. Nothing is dropped here: Stage 6's design matrix drops constant
  columns at the point of use [§7.1].` Table `covariate | distinct values`, only when `n > 0`.

Five of the six are unconditional within their own entry point; only `constant_covariates`' table is
conditional, and the entry itself is recorded with `n = 0` when nothing is constant — the same shape
Stage 2 gives `penumbra_recomputed`. §12.12 asserts the inventory against this table, so an entry
cannot be dropped, renamed or reordered.

## 10. What Stage 3 amends in Stages 1 and 2

The full ledger, so that no amendment is discovered during implementation.

| File | Amendment | Why |
|---|---|---|
| `config.py` | add `ONSET_TYPE_FROM_FLAG` (§4.2), `SUBGROUPS` and `COHORT_DEPENDENT_SUBGROUPS` (§6.5, §6.4), `DERIVED_DICHOTOMIES` / `ROW_WISE_DERIVED` and `DERIVED_NAMES` as computed views, **each in the file position §8.4 pins** — `DERIVED_NAMES` moves out of the covariates block and leaves a pointer comment behind | the registries Stage 3 reads instead of literals; computed views must follow the registries they read or the module does not import |
| `config.py` | delete `TARGET_MISMATCH`, `MISMATCH_RATIO_UNDEFINED_COUNTS_AS_MET`, and the comment block that explains the `core == 0` convention | the [§13] amendment leaves them no consumer (§6.1) |
| `data.py` | `_audit_path`, the header title and footer (§9.1); `KINDS` and `_HEADINGS` (§9.2); `_missingness` split into its structural loop plus a public `absence_by_column` (§8.2) | one log across the pipeline, one classification |
| `data.py` | delete the `crossing` computation and its sentence from `_correct` | it reports a [§13] consequence that no longer exists (§6.1). The per-case table already carries stored, recomputed and difference, so the entry loses no information |
| `test_data.py` | 12.8's path assertions; 12.16's heading list gains `Derivations`; the `_HAND_RECORDS[0]` comment that explains the copy-paste record in terms of `TARGET_MISMATCH["min_penumbra_ml"]` | follows the three edits above |
| `test_data.py` — **a fix, not a follow-on** | `test_all_six_headings_appear_and_an_empty_section_says_none` asserts the text between the `observation` and `structural` headings is exactly `_none_`, so inserting `derivation` between them **fails it**. Take the neighbour from `KINDS` rather than naming `structural`, and rename the test off "six" (§9.2) | the one existing test any new `kind` breaks. Listed separately because it is a green test going red, not a consequence of an edit elsewhere, and T3 is where it is caught |
| `test_config.py` | `ONSET_TYPE_FROM_FLAG`'s levels plus the reference level equal `FACTOR_LEVELS["onset_type"]`; `DERIVED_NAMES` disjoint from `ANALYSIS_NAMES`; no `SUBGROUPS` key appears in any covariate list; `COHORT_DEPENDENT_SUBGROUPS ⊆ SUBGROUPS` | §12.13 |
| `test_config.py` | 9.13's `parametrize` union gains `SUBGROUPS`, `COHORT_DEPENDENT_SUBGROUPS`, `DERIVED_DICHOTOMIES` and `ROW_WISE_DERIVED`; 9.8's tuple list gains the four new tuples | both lists name their constants **explicitly**, so a new column-keyed constant is outside them until added — it would ship with no name-resolution and no immutability assertion (§12.13) |
| `specs/stage1_config_and_data_contract.md` §6, §13 | the target-mismatch half of §6 is replaced by a pointer to the amendment; §13's "lift `TARGET_MISMATCH` verbatim" row goes | they specify deleted constants |
| `specs/stage2_load_and_clean.md` §7.2, §9.6, §14, §18 | the target-mismatch consequence paragraph, the `penumbra_recomputed` detail template's second sentence, §14's "Whether `target_mismatch` is met" line, §18's criteria row | same |

**No test asserts the removed crossing sentence.** Verified by reading `test_data.py`: 12.6 asserts
the entry's `n`, `case_ids` and table, and 12.16 asserts the `(kind, step)` inventory. So that edit
costs no test, only the comment named above.

## 11. Data flow into Stages 4, 5, 6 and 10

```
  derive(df, audit) returns
   ├─ every column Stage 2 delivered, unchanged — no value edited, no row dropped
   ├─ onset_type        string, in FACTOR_LEVELS["onset_type"]; <NA> only where a flag is <NA>
   ├─ unknown_onset     Int64 0/1; <NA> exactly where onset_type is
   ├─ mrs_0_2_90d  mrs_0_1_90d  death_90d  mrs_5_6_90d
   │                    Int64 0/1; <NA> exactly where mrs_90d is        [invariant 6]
   └─ 31 columns: 25 analysis + 6 derived. Default RangeIndex, still not an identifier.

  derive_cohort(df, audit) returns
   ├─ core_above_median Int64 0/1; <NA> where core_ml is; ties below the median
   └─ 33 columns. Raises if called twice on the same frame.               [§6.4]

  audit
   └─ the same Audit object load() created, carrying Stage 2's seven entries and now
      Stage 3's six. Stage 4 appends; nobody re-opens. No file written unless asked.
```

What each later stage may rely on:

- **Stage 4** — `ivt_contraindicated` and `contraindication_reason` are untouched. Stage 3 reads
  neither, so DECISION 1's classifier meets the frame exactly as Stage 2 left it.
- **Stage 5** — must call `derive_cohort` **after** both [§3] restrictions and before any estimate.
  It is the only caller.
- **Stage 6** — `onset_type` is a plain `string`, so `FACTOR_LEVELS["onset_type"]` is a superset of
  the levels present and `pd.Categorical` construction at the point of use is well defined even if a
  level is empty.
- **Stage 10** — resamples `core_above_median`; never recomputes it. §6.4 gives the reason and the
  guard that makes the alternative raise.
- **Stage 11** — the two [§13] subgroups are `SUBGROUPS`' keys, with their cell sizes in the audit
  log, which is where the not-estimable gate reads them.

## 12. Acceptance criteria

All tests live in `test_derive.py`, with section banners matching these numbers, as `test_data.py`
does for Stage 2. Tests needing the private workbook reuse the `DATA_GATED` idiom; everything else
runs against hand-built frames. `test_derive.py` is not exempt from the raw-name scan and builds its
frames with analysis names.

The hand-built frame is `test_data.py`'s `hand_frame()`, imported rather than re-declared: two hand
frames drift, and this one is already contract-valid over six records, so any error a test sees comes
from the corruption that test applied. It needs **no extending** — every column Stage 3 reads is
already in it (`wake_up`, `unwitnessed` and `mrs_90d` in `_HAND_CONSTANT`, `core_ml` per record) — and
it already encodes the missing-volume and both-arm cases §6.3 and §7 need.

**What it needs is variation, and there is exactly one tool for that.** `hand_frame(**overrides)`
assigns `df[column] = value`, a scalar broadcast over all six rows, so it can make a column uniform
but cannot give two records different onset flags. The per-record tool is `test_data.corrupt(column,
value, where="HAND-N")`, which matches on `case_id` and never on row position. So §12.2's three-level
frame is built by `corrupt`-ing one record's `wake_up` to 1 and another's `unwitnessed` to 1, not by
passing overrides. Both helpers are imported from `test_data`; neither is reimplemented.

**The frame's derived expectations, pinned, because they fall out of the existing records and must not
be hand-tuned into a new frame.** `core_ml` over the six records is `(13, 20, None, 5, 0, 30)`:

```
  median (skipna, over 0, 5, 13, 20, 30)   13.0        ← equals HAND-1's own value
  above  (> 13.0)                          2           HAND-2 (20), HAND-6 (30)
  at     (== 13.0)                         1           HAND-1  → core_above_median = 0  [ties below]
  below  (< 13.0)                          2           HAND-5 (0), HAND-4 (5)
  missing                                  1           HAND-3  → core_above_median = <NA>
```

So §12.10's tie case and its missing case both already exist in the frame as it stands — the tie is
not a construction, it is what a median over an odd count of retained values does.

`_HAND_CONSTANT` sets the remaining inputs to one value each, and the consequences are worth writing
out because two of them look like they say something they do not:

```
  wake_up = 0, unwitnessed = 0 on every record
      → onset_type is uniformly the baseline, unknown_onset uniformly 0.
        §12.2 and §12.3 each corrupt a specific record's flag to move off it.

  mrs_90d = 2 on every record
      → the four dichotomies are NOT uniformly events. Verified:
            mrs_0_2_90d   2 <= 2   →  1  on all six
            mrs_0_1_90d   2 <= 1   →  0  on all six
            death_90d     2 == 6   →  0  on all six
            mrs_5_6_90d   2 >= 5   →  0  on all six
        One event column and three non-event columns, which is what makes the frame usable
        for §12.5 at all: a frame where every dichotomy read the same would not distinguish
        a registry-driven loop from a constant. §12.7 corrupts one record's mrs_90d to <NA>.

  ten of the thirteen PS_COVARIATES_FULL entries are single-valued
      → see §12.11. This is a property of the frame, not a defect in it.
```

1. **The onset assertion fires.** A frame with `wake_up == 1` and `unwitnessed == 1` on one record
   raises `SchemaError` naming that record's identifier. A frame where each flag is 1 on a
   *different* record does not raise. A frame where one flag is `<NA>` and the other is 1 does not
   raise — a missing flag is not evidence of both being positive (§4.1).
2. **The three levels partition the frame.** Every record carries exactly one of
   `FACTOR_LEVELS["onset_type"]`; no fourth value appears; the three counts sum to the row count
   minus the `<NA>` count. Asserted on a hand-built frame covering all three levels.
3. **A missing flag yields `<NA>`, never the baseline and never a level.** A record with
   `wake_up = <NA>` and `unwitnessed = 0` gets `<NA>` for `onset_type` — not `witnessed` — and `<NA>`
   for `unknown_onset`. This is the §4.3 branch that no data-gated test can reach. It also pins the
   two guards separately, because §4.2's square shows each one failing differently and the combined
   assertion sees neither: with the `unknown` mask removed the record must be `witnessed`, and with
   `.fillna(False)` **also** removed it must be `wake_up` — asserted against locally reimplemented
   two-line variants of the loop in the test, never by monkey-patching `derive.py`, so the shipped
   function keeps both guards and the test still documents what each one buys.
4. **The levels come from the configuration.** `derive.py` contains no string literal equal to any
   member of `FACTOR_LEVELS["onset_type"]` — an AST constant scan over the module, the same technique
   as Stage 1's 9.4 and Stage 2's 12.10, with a companion test proving the scan fires on a synthetic
   snippet that does contain one.
5. **The dichotomies.** Exactly one column per `OUTCOMES` entry with `source is not None`, and none
   for the three read directly. Each takes values in `{0, 1}` or `<NA>`. Parametrised over the
   registry: the column equals `OPS[o.op](src, o.threshold)` wherever the source is present.
6. **Registry-driven, not literal.** Run under a registry copy whose `mrs_0_2_90d` threshold is 3;
   the produced column follows the patched threshold. If it does not, a threshold is written in
   `derive.py`.
7. **Missingness is reimposed** — roadmap Stage 3's acceptance criterion and invariant 6.
   Parametrised over every derived dichotomy **× both source dtypes**: its `isna()` equals its
   source's `isna()`, exactly. And the trap in its own test: on a frame where `mrs_90d` is `<NA>` for
   one record, that record is `<NA>` in all four dichotomies and is **not** 0 in any of them.

   The second dtype is what makes this criterion testable at all. Per §5.2 the mask is a no-op while
   `mrs_90d` is `Int64`, so the `Int64` half of the parametrisation passes with `.mask(src.isna())`
   deleted. The `float64` half — the same frame with `mrs_90d` cast, `hand_frame().assign(...)` or an
   explicit `astype("float64")`, which is the shape any frame built without `READ_DTYPES` has — is the
   half that fails. Definition of done 5 is watched against that half, and the test's own comment says
   so, because a future reader who runs only the `Int64` case will otherwise conclude the mask is dead
   code and delete it.
8. **`unknown_onset` follows `onset_type`.** It equals `onset_type != REFERENCE_LEVELS["onset_type"]`
   with missingness carried; and it is derived from `onset_type` rather than from the flags —
   asserted by overwriting `onset_type` on a frame and confirming `unknown_onset` moves with it
   (§6.2).
9. **The derived missingness table.** Its rows are exactly `sorted(ROW_WISE_DERIVED)` — not
   `DERIVED_NAMES`, since `core_above_median` does not exist yet (§8.3). Every row's kind is `missing`
   or `complete`, never `structural` or `not recorded`. `n + n_absent` reconciles to the row count on
   every row, and the per-centre counts sum to the overall count on every row.
10. **`core_above_median`.** `derive` alone does not create it. `derive_cohort` on a frame that
    already has it raises `SchemaError` naming the column and §6.4 (the frozen-median guard). On the
    hand frame the whole split is asserted against §12's pinned table — median `13.0`, and
    `2` above / `1` at / `2` below / `1` missing — so ties-sit-below is checked on HAND-1, which equals
    the median, and the missing case on HAND-3, which has no core volume. Neither needs a purpose-built
    record. The median is that of the frame passed in — asserted by calling it on a subset and on the
    whole and confirming the two medians differ where the data make them differ. And the load-bearing
    mask (§6.3): with `.mask(core.isna())` removed, HAND-3 becomes 0 rather than `<NA>`, which — unlike
    12.7's `Int64` half — fails immediately, because `core_ml` is `float64` and absorbs.
11. **`constant_covariates`.** On the hand frame it returns **ten of the thirteen**
    `PS_COVARIATES_FULL` entries, not `()` — `_HAND_CONSTANT` gives each of them a single value. The
    expectation is pinned as a literal tuple rather than described, in `PS_COVARIATES_FULL` order,
    which is also what proves the caller's-order clause:

    ```python
    ("age", "sex", "prestroke_mrs", "nihss_baseline", "onset_type",
     "atrial_fib", "hypertension", "hyperlipidemia", "diabetes", "smoking")
    ```

    The three that are absent are the three that vary: `core_ml` (5 distinct), `tmax6_ml` (5) and
    `center` (4). Note what the frame gives for free — `onset_type` is in the tuple because every
    record is the baseline, so **"finds a factor reduced to one level" is covered by the baseline
    assertion itself**, on a `string` column, with no construction. The remaining discovery cases are
    asserted as *additions* to this tuple: a numeric column made constant through
    `hand_frame(core_ml=…)` joins it, and an all-missing column joins it too (`nunique(dropna=True)`
    is 0 there, per §7.1). And it mutates nothing: `C.OUTCOME_COVARIATES is C.PS_COVARIATES` still
    holds after a call, and the frame is unchanged (§7.1).

    **Do not "fix" the frame to make this return `()`.** It cannot be done through
    `hand_frame(**overrides)`, which assigns a scalar and therefore makes a column *more* constant,
    and doing it through ten `corrupt` calls would build a second hand frame by the back door — the
    thing §12 refuses. Asserting the ten is the stronger test anyway: it is the only assertion in
    §12 that would catch `constant_covariates` silently ranging over a set instead of the caller's
    sequence.
12. **The audit inventory.** The six §9.3 entries appear with the declared `kind` and `step` strings,
    in the declared order, `derive`'s three-plus-one before `derive_cohort`'s two. The rendered log
    carries a `## Derivations` heading. `constant_covariates` is present with `n = 0` and no table
    when nothing is constant. Reproduction stays byte-identical, within a process and across
    interpreters launched with `PYTHONHASHSEED=0` and `=1` — Stage 3 adds two dict-ordered loops
    (`OUTCOMES`, `SUBGROUPS`) and one sorted one, and none may become a set.

    The cross-process half needs a driver, and it is written out rather than sketched, for the reason
    Stage 2 §9.6 gives about pinning log contents: an implementer choosing the driver freely can
    choose one that renders no derivation at all and still see two identical outputs. It follows
    `test_data.py`'s existing two-seed test exactly, one `subprocess.run` per seed with
    `cwd=MODULE_DIR` and `env={**os.environ, "PYTHONHASHSEED": seed}`:

    ```python
    script = ("import sys, data, derive\n"
              "df, audit = data.load(data.FIXTURE)\n"
              "df = derive.derive(df, audit)\n"
              "derive.derive_cohort(df, audit)\n"
              "sys.stdout.write(audit.to_markdown())\n")
    ```

    `FIXTURE` and not `WORKBOOK`: the workbook is gitignored, so a data-gated reproduction test is no
    test on most checkouts. The fixture reaches every entry — verified against
    `tests/fixture_schema.xlsx` as committed: 2 records, `wake_up` 1 and 0 so `onset_type` renders two
    of its three level rows and its `<NA>` row, `mrs_90d` `<NA>` on one record so the dichotomies
    render a missing count, `core_ml` 0 and 12 giving a median of 6.0 with one record either side, and
    both onset flags present so §4.1 does not raise. Note the dtype asymmetry §3.1 names: the
    fixture's `core_ml` arrives as `int64` because both its values are integers, where the workbook's
    is `float64` — so this test exercises the `int64` path and 12.10 exercises `float64`. Neither
    replaces the other.
13. **Stage 1 and Stage 2 amendments.** In `test_config.py` and `test_data.py`, not here:
    `ONSET_TYPE_FROM_FLAG`'s levels plus `REFERENCE_LEVELS["onset_type"]` equal
    `FACTOR_LEVELS["onset_type"]`; `DERIVED_NAMES ∩ ANALYSIS_NAMES = ∅`;
    `COHORT_DEPENDENT_SUBGROUPS ⊆ SUBGROUPS`; no `SUBGROUPS` key appears in `PS_COVARIATES_FULL`,
    `BALANCE_ONLY` or any outcome model's covariates; the renamed audit path; the seven headings; and
    the heading-adjacency test repaired per §9.2.

    Two of Stage 1's existing sweeps have to be extended by hand, and this is the clause that says so
    because neither extends itself:

    - **9.13, `test_every_named_variable_can_be_produced`.** Its `parametrize` union names its six
      constants explicitly, so `SUBGROUPS`, `COHORT_DEPENDENT_SUBGROUPS`, `DERIVED_DICHOTOMIES` and
      `ROW_WISE_DERIVED` are outside it until added. A typo in any of them would otherwise surface as a
      `KeyError` inside `derive` or, worse, inside Stage 6.
    - **9.8, `test_declared_sequences_are_tuples`.** Same shape, an explicit twelve-name list. The four
      new tuples join it, so none can be softened to a list that a later stage appends to.

    **And one existing guarantee is weakened by this stage, which is why the disjointness assertion
    above is load-bearing rather than tidy.** 9.3's `test_every_ps_covariate_resolves` accepts a
    covariate that is `in ANALYSIS_NAMES or in DERIVED_NAMES`. Widening `DERIVED_NAMES` to carry the two
    subgroup keys means a subgroup wrongly added to `PS_COVARIATES` would now *pass* that test. The
    "no `SUBGROUPS` key in any covariate list" assertion is what restores the guarantee 9.3 used to give
    on its own; it is not a second opinion on the same property.
14. **No bare `assert` in `derive.py`.** An AST scan for `ast.Assert`, as Stage 2 §12.10, with the
    companion test that proves the scan fires.
15. **Structural facts from the workbook.** Data-gated, against v7 and the §18 record: onset levels
    33 witnessed / 36 unwitnessed / 57 wake-up over 126 records, with neither flag ever missing and
    no record positive on both; the four dichotomies reproduce Stage 0's full-frame counts — mRS 0–2:
    53 events, mRS 0–1: 31, death: 37, mRS 5–6: 42, each with 2 missing; `unknown_onset` is 1 on 93
    records; and over the 93-record primary cohort the core median is 5.0 mL with 44 above, 3 at, 45
    below and 1 missing, against 6.0 mL over all 126.

### Coverage map

```
derive.py                                          test_derive.py
├── _assert_onset_flags(df)
│   ├── both flags 1        → SchemaError          ├── 12.1  identifier named
│   ├── each on a different record → no raise      ├── 12.1
│   └── one flag <NA>       → no raise             └── 12.1
│
├── _onset_type(df)
│   ├── three levels, partition                    ├── 12.2  built with corrupt(), not overrides
│   ├── flag <NA>           → <NA>, not baseline   ├── 12.3  unreachable on v7
│   ├── unknown mask alone removed → "witnessed"   ├── 12.3  §4.2's square, cell 2
│   ├── both guards removed → "wake_up" (level!)   ├── 12.3  §4.2's square, cell 4
│   ├── baseline from REFERENCE_LEVELS             ├── 12.4  AST scan + fires
│   └── string, not Categorical                    └── 12.2
│
├── _unknown_onset(onset)
│   ├── complement of the baseline level           ├── 12.8
│   ├── <NA> carried, never 1                      ├── 12.8 / 12.3
│   └── derived from onset_type, not the flags     └── 12.8  overwrite test
│
├── _dichotomies(df)
│   ├── one per registry source, none otherwise    ├── 12.5
│   ├── values in {0, 1} or <NA>                   ├── 12.5
│   ├── threshold from OUTCOMES                    ├── 12.6  patched registry
│   ├── missingness, Int64 source                  ├── 12.7  passes without the mask — by design
│   ├── missingness, float64 source                ├── 12.7  THE one that fails; DoD 5 watches it
│   ├── a missing source is never a 0              ├── 12.7  the trap, its own test
│   └── writes in place; caller owns the frame     └── 12.5/12.6/12.7 each build their own
│
├── absence_by_column(…, ROW_WISE_DERIVED, …)
│   ├── rows are ROW_WISE_DERIVED, not DERIVED_    ├── 12.9
│   ├── kind is missing or complete only           ├── 12.9
│   └── reconciles per row and per centre          └── 12.9
│
├── _core_above_median(df)
│   ├── absent after derive alone                  ├── 12.10
│   ├── second call            → SchemaError       ├── 12.10  the frozen-median guard
│   ├── ties sit below                             ├── 12.10  HAND-1 == median 13.0; 3 in v7 [§18]
│   ├── missing core           → <NA>              ├── 12.10  HAND-3; float64 absorbs, so the
│   │                                              │          mask here IS load-bearing  [§6.3]
│   ├── all-missing core → median NaN → all <NA>   ├── none — degrades to empty, log says
│   │                                              │          `missing`  [§13]
│   └── median is of the frame given               └── 12.10 / 12.15
│
├── constant_covariates(df, covariates)
│   ├── 10 of 13 on the hand frame, pinned tuple   ├── 12.11  NOT () — _HAND_CONSTANT
│   ├── 0 of 13 on the fixture                     ├── 12.12  which is why n = 0 there
│   ├── one-level factor (string)                  ├── 12.11  free: onset_type is in the ten
│   ├── column made constant / all-missing         ├── 12.11  asserted as additions to the ten
│   ├── caller's order                             ├── 12.11  the pinned tuple IS the order test
│   ├── pre-derive frame → KeyError on onset_type  ├── none — docstring precondition  [§7.1]
│   └── mutates neither frame nor list             └── 12.11  identity assertion
│
└── the audit
    ├── six entries, declared kinds and steps      ├── 12.12
    ├── `## Derivations` heading                   ├── 12.12
    ├── constant_covariates with n = 0             ├── 12.12
    ├── byte-identical, same process               ├── 12.12
    ├── byte-identical, PYTHONHASHSEED 0 vs 1      ├── 12.12
    └── no bare assert in derive.py                └── 12.14  + the scan fires

config.py / data.py amendments                     test_config.py / test_data.py
├── ONSET_TYPE_FROM_FLAG resolves to FACTOR_LEVELS ├── 12.13
├── DERIVED_NAMES ∩ ANALYSIS_NAMES = ∅             ├── 12.13
├── COHORT_DEPENDENT_SUBGROUPS ⊆ SUBGROUPS         ├── 12.13
├── no subgroup in any covariate list              ├── 12.13  restores what widening
│                                                  │          DERIVED_NAMES weakens in 9.3
├── 4 new constants in 9.13's resolvable union     ├── 12.13  the union is explicit
├── 4 new tuples in 9.8's immutability list        ├── 12.13  the list is explicit
├── computed views follow their registries         ├── none — NameError at import  [§8.4]
├── the renamed audit path                         ├── 12.13 (12.8 there)
├── heading adjacency, neighbour from KINDS        ├── 12.13 — a RED test made green  [§9.2]
└── seven headings incl. Derivations               └── 12.13 (12.16 there)

Data-gated (skipif not DATA_XLSX.exists()):        └── 12.15  33/36/57, 53/31/37/42,
                                                            93 unknown, median 5.0
                                                            (44/3/45/1) vs 6.0
Every branch above has a test. No branch is untested.
```

### Failure modes

| Failure | Test | Error handling | Visible or silent |
|---|---|---|---|
| A record positive on both onset flags | 12.1 | `SchemaError` naming the cases | visible |
| A missing onset flag read as `witnessed` | 12.3 | the `unknown` mask; a fabricated covariate value becomes `<NA>` | visible |
| A missing onset flag read as `wake_up` — a fabricated *level*, not the baseline | 12.3 | `.fillna(False)`; §4.2's square, the failure the wrong rationale hid | visible |
| A missing mRS read as a non-event, `float64` source | 12.7 | `.mask(src.isna())`; invariant 6 | visible |
| A missing mRS read as a non-event, `Int64` source | 12.7 | unreachable — the dtype propagates (§3.1, §5.2). The mask is the guarantee that this stays true if the dtype changes | n/a today |
| A missing core volume swept into the below-median group | 12.10 | `.mask(core.isna())`, the one mask load-bearing today because `core_ml` is `float64` | visible |
| An all-missing `core_ml` producing an empty subgroup | none | degrades to all-`<NA>`; the log prints the median as `missing` and both counts as 0 | visible, untested — §13 |
| A threshold typed a second time in `derive.py` | 12.6 | patched-registry test diverges | visible |
| A level name typed in `derive.py` | 12.4 | AST constant scan | visible |
| `unknown_onset` rebuilt from the flags | 12.8 | overwrite test diverges | visible |
| The median taken over the wrong population | 12.10, 12.15 | unrepresentable through `derive`, which cannot compute it | visible |
| A bootstrap replicate recomputing the median | 12.10 | `SchemaError` from the §6.4 guard on the second call | visible |
| A covariate silently dropped for zero variance | 12.11 | nothing is dropped here; the entry names what is constant | visible |
| A covariate list mutated in place | 12.11 | identity assertion on `OUTCOME_COVARIATES is PS_COVARIATES` | visible |
| A derived column labelled `structural` | 12.9 | kind assertion over every row | visible |
| A derived column shadowing an analysis column | 12.13 | disjointness assertion | visible |
| A subgroup reaching a covariate list | 12.13 | membership assertion — 9.3's own resolution check no longer catches it (§12.13) | visible |
| A new registry constant outside 9.13's union or 9.8's list | 12.13 | both lists are explicit; the amendment adds all four names | visible |
| A computed view declared above the registry it reads | none | `NameError` at import, on every test, naming the constant | visible |
| `derivation` inserted between two headings an existing test reads as adjacent | 12.13 | the neighbour is taken from `KINDS`, not named (§9.2) | visible |
| A test asserting on a frame `_dichotomies` mutated under it | — | contract stated in §0.2 and §5.1; each test builds its own frame | convention |
| An audit entry dropped, renamed or reordered | 12.12 | inventory comparison against §9.3 | visible |
| The log made non-reproducible across processes | 12.12 | two-seed render comparison | visible |
| A check written as `assert`, run under `-O` | 12.14 | AST scan at test time | visible |
| A near-constant subgroup reported as a finding | none | Stage 11's not-estimable gate, not built yet | **deferred** — §13 |

One row is deferred rather than silent: the [§13] cell-size gate belongs to Stage 11, and §6.5 states
what Stage 11 must do. Everything else is visible at the point of failure, and two are listed as
unrepresentable — `derive` cannot compute a median, and a second `derive_cohort` cannot proceed.

## 13. Known gaps carried forward

- **A missing onset flag is permitted upstream and would drop a patient from every [§6] model.**
  Stage 2's A5 allows missing in every `BINARY_COLUMNS` column, so `onset_type` may legitimately be
  `<NA>`. v7 has none (§18). Not promoted to an assertion: a missing flag is ordinary missingness and
  [§11] is complete-case with the denominator reported. What §4.3 guarantees is that it is *counted*
  as missing rather than absorbed into the baseline.
- **The median split can move.** §18 records 3 records sitting exactly on the primary-cohort median
  and 39 with a core of exactly 0. A corrected volume in a future workbook can shift the median and
  reshuffle membership near it. That is a property of a median split, not a defect; it is recorded so
  it is not later diagnosed as one, and it is why the audit entry prints the median and the tie count.
- **`penumbra_ml` now feeds only the balance table.** After the [§13] withdrawal it is a
  `BALANCE_ONLY` variable and nothing else. Stage 2's recomputation of it stands on [§6]'s definition
  alone. If a future amendment restores the target-mismatch subgroup, §6.1's boundary question must
  be answered first, in writing.
- **The two records whose shipped death flag contradicts their mRS** remain a standing query
  [DECISION 2]. `death_90d` derives from the mRS, so it is unaffected unless the mRS itself changes —
  in which case the query returns as a new workbook and a new `DATA_SHA256`.
- **An all-missing `core_ml` would give an empty median subgroup rather than a wrong one.** §6.3
  records the mechanism: `median()` returns `NaN`, the comparison is `False` everywhere, the mask sets
  every row to `<NA>`. No guard and no test, deliberately — the audit entry prints the median as the
  literal `missing` and both counts as 0, so it is visible where it matters, and a cohort all-missing on
  a [§6] covariate cannot get past Stage 2's A6 and the balance table without being noticed first.
  Recorded so it is not later diagnosed as a derivation bug.
- **The dtype of `core_ml` and `tmax6_ml` is source-dependent** because the contract declares none for
  them (§3.1): `float64` from the workbook, `int64` from a source whose values are all integers, as
  `tests/fixture_schema.xlsx` is. Nothing in Stage 3 depends on which — every comparison is masked — but
  a future stage that *does* depend on it should declare the dtype in Stage 1 rather than cast at the
  point of use. `../../TODOS.md` carries that work, and **§15 defers it until after Stage 3 lands.**

  One caveat for whoever does it, because it makes a sentence in this spec stale rather than merely
  imprecise: the choice between `float64` and the nullable `Float64` is not a formality. Declaring
  `float64` keeps things as they are. Declaring **`Float64`** makes `core > median` *propagate* `<NA>`
  instead of absorbing it (§3.1), which turns `_core_above_median`'s `.mask(core.isna())` from the one
  load-bearing mask in Stage 3 into a redundant one — so §6.3's third bullet, §12.10's last clause,
  Definition of done 6's last sentence and the failure-modes table would all then be describing a
  guarantee that a dtype provides rather than a mask. The mask stays either way, per §3.1; what changes
  is which sentences are true about why. Update them in the same commit.
- **`sex` coding is undocumented**, carried from Stage 1 §10. `SEX_LABELS` stays `None`.

## 14. What Stage 3 deliberately does not decide

Recorded so a later stage does not look here for an answer that was never placed here.

- **Which patients are eligible.** Stage 4. Stage 3 reads neither `ivt_contraindicated` nor
  `contraindication_reason`.
- **Which centres and patients are dropped.** Stage 5, which is also `derive_cohort`'s only caller.
- **How factors are coded, and which constant columns are dropped from a design matrix.** Stage 6,
  from `FACTOR_LEVELS` and `REFERENCE_LEVELS`, at the point of use (§7.1).
- **Whether a subgroup is estimable.** Stage 11. Stage 3 supplies the flag and its cell sizes; the
  gate that refuses to report an interaction test on a degenerate split is Stage 11's (§6.5).
- **What happens to a bootstrap replicate.** Stage 10. Stage 3 states one requirement — resample the
  column, do not recompute it — and enforces it with the §6.4 guard.

## 15. NOT in scope for Stage 3

| Considered | Why deferred |
|---|---|
| Restoring the target-mismatch subgroup | Withdrawn by the [§13] amendment of 2026-08-10 (§6.1). Restoring it is an amendment, not an implementation decision. |
| Dropping a zero-variance covariate | Stage 6, at the point of use (§7). Detecting it is here; deleting it is not. |
| Any imputation | Never [§11]. Every derived column carries its source's missingness by construction. |
| Deriving `mrs_shift_90d`, `center_hug`, or any 24-hour NIHSS dichotomy | Not in the [§5] registry. The pilot has all of them; restoring any requires amending [§5] first. |
| Reading a shipped dichotomy to cross-check a derived one | DECISION 2, and the contract has dropped them (§5.3). |
| Making `onset_type` a `Categorical` | Stage 6, after the Stage 5 restrictions (§4.4). |
| Caching or memoising a derivation | 126 rows and one vectorised pass. A cache is a second source of truth (§2). |
| A CLI entry point | Stage 14 is the single entry point [§16]. `derive.py` is a library module. |
| Declaring an explicit dtype for `core_ml`, `tmax6_ml` and `penumbra_ml` | A Stage 1 amendment, captured in `../../TODOS.md` with its full context. Deferred **until after Stage 3 lands**, not because it is unimportant — §3.1 shows the inferred dtype is what decides whether a mask is load-bearing — but because §18 records the current dtypes as verified facts and §12.10 and §12.12 are written against them. Changing the contract mid-implementation would move the ground under the spec the implementation is checked against. See the caveat in §13 before doing it. |
| A guard for an all-missing `core_ml` | It degrades to an empty subgroup rather than a wrong one and the log says so (§6.3). Recorded in §13. |

## 16. What already exists, and what to lift

`pilots/data.py` is gitignored but present and is the closest prior art. `stage0_data_inventory.py`
is in the repository and its `derive` helper is already correct.

| From the pilot / Stage 0 | Status |
|---|---|
| `stage0_data_inventory.py:97`'s `derive(source, is_event)` — `is_event.astype("Int64").mask(source.isna())`, with the docstring explaining why | **lift the idiom.** It is what Stage 1 §3.2 pins and what §5.1 writes as a loop over the registry |
| `pilots/data.py:31`'s `_rebuild(df, col, truth, source, audit, rule)` | lift the shape — a derivation that logs — but drive it from `OUTCOMES` rather than from eight hand-written call sites, four of which build columns the [§5] registry does not contain |
| `pilots/data.py:140`'s `assert not (wake_up & unwitnessed).any()` | lift the **check**; not the statement form (Stage 2 §8.3), and not the message. §4.1 raises `SchemaError` and names the cases |
| `pilots/data.py:141`'s `np.select([...], [...], "witnessed")` | **do not lift.** The default absorbs a missing flag into the baseline, which is exactly §4.3 |
| `pilots/data.py:141`'s `pd.Categorical(...)` at derivation time | **do not lift.** §4.4 |
| `pilots/data.py:144`'s `unknown_onset = (onset_type != "witnessed").astype(int)` | lift the definition; not the `astype(int)`, which cannot hold `<NA>`, nor the literal level name (§6.2) |
| `pilots/data.py:150-154`'s `target_mismatch` | **do not lift.** Withdrawn (§6.1) |
| `mrs_shift_90d`, `center_hug`, the 24-hour NIHSS dichotomies | **do not lift.** Not in the [§5] registry (§15) |
| the pilot's habit of writing thresholds and level names at the call site | **do not lift.** §4.2, §5.1 |

## 17. Implementation tasks

Ordered. Each independently verifiable. T1–T3 touch Stage 1 and Stage 2 files; T4–T8 touch
`derive.py` and `test_derive.py`. The second lane depends on the first, so there is one lane in
practice and no parallelism worth taking.

- [x] **T1 (P1)** — `config.py`: add `ONSET_TYPE_FROM_FLAG`, `SUBGROUPS`,
      `COHORT_DEPENDENT_SUBGROUPS`, `DERIVED_DICHOTOMIES`, `ROW_WISE_DERIVED`; make `DERIVED_NAMES` a
      computed view. **Each in the file position §8.4 pins** — computed views after the registries they
      read, `DERIVED_NAMES` out of the covariates block with a pointer comment left behind. `test_config.py`:
      the §12.13 assertions, **including** the four names added to 9.13's `parametrize` union and the four
      tuples added to 9.8's list. Verify: acceptance 12.13 — and see the level-set assertion fail first,
      by temporarily renaming a level in `ONSET_TYPE_FROM_FLAG`, since nothing covered these constants
      before this task. If `import config` raises `NameError`, the placement is wrong, not the constant.
- [x] **T2 (P1)** — the [§13] withdrawal in code: delete `TARGET_MISMATCH` and
      `MISMATCH_RATIO_UNDEFINED_COUNTS_AS_MET` from `config.py`; delete the `crossing` computation and
      its sentence from `data.py`'s `_correct`; fix the `_HAND_RECORDS[0]` comment. Verify:
      `uv run pytest -v` green with no test edited — §10 records that none asserts the removed
      sentence, and this task is where that claim is checked rather than believed.
- [x] **T3 (P1)** — `data.py`: `KINDS` and `_HEADINGS` gain `derivation`; `_audit_path`, title and
      footer per §9.1; `_missingness` split into its structural loop plus the public
      `absence_by_column`. `test_data.py`: 12.8's path, 12.16's headings, **and the heading-adjacency
      repair of §9.2** — take `observation`'s successor from `KINDS` instead of naming `structural`, and
      rename the test off "six". Verify: acceptance 12.13's last clauses, and **Stage 2's own 12.11,
      12.13 and 12.16 all green** — 12.16 is where the adjacency test lives and it goes red the moment
      `KINDS` gains the kind, so run the whole file, not the two sections the generalisation touches.
      The generalisation must not change the table Stage 2 already produces.
- [x] **T4 (P1)** — `derive.py`: `_assert_onset_flags`, `_onset_type`, `_unknown_onset`, and their
      two audit entries. Verify: acceptance 12.1, 12.2, 12.3, 12.4, 12.8.
- [x] **T5 (P1)** — `_dichotomies` and its entry. Verify: acceptance 12.5, 12.6, 12.7.
- [x] **T6 (P1)** — `derive` wiring, including the `absence_by_derived_column` call. Verify:
      acceptance 12.9.
- [x] **T7 (P1)** — `_core_above_median`, `constant_covariates`, `derive_cohort` with the §6.4 guard,
      and their two entries. Verify: acceptance 12.10, 12.11.
- [x] **T8 (P2)** — `test_derive.py` remainder: the inventory and reproduction tests, the AST scans,
      and the data-gated §18 facts. Verify: `uv run pytest -v` green both with and without `data/`.

### Diagrams that belong in the code, not only here

Three, and no more — a diagram nobody maintains is worse than none, because it is believed. Keeping
them true is part of any change that touches them, in the same commit.

- **`derive.py`'s module docstring** — §0's two entry points and where each is called from, including
  the arrow from Stage 5 into `derive_cohort`. It is the first thing a reader sees and it is the one
  fact about this module that is not obvious from its functions.
- **Above `_onset_type`** — the three levels, the two flags, and the two ways a record leaves the
  baseline branch: a positive flag, or a missing one. The second is §4.3 and is the thing an
  implementer "tidying" the mask will delete.
- **Above `_core_above_median`** — the median, the tie rule, and the frozen arrow into Stage 10. The
  guard's error message says it too, but the guard is only read after it fires.

### Definition of done

Stage 3 is complete when all of the following hold, and not before.

1. `uv run pytest -v` is green **with** `data/` present — every test, including the data-gated ones.
2. `uv run pytest -v` is green **with `data/` temporarily renamed** — the data-gated tests skip and
   nothing else fails or errors at collection.
3. `test_config.py` and `test_data.py` are still green, and neither `derive.py` nor `test_derive.py`
   has been added to `EXEMPT_FROM_RAW_NAME_SCAN`.
4. The pipeline's audit log has been produced through both entry points, against the **workbook**, and
   reproduces byte-identically across hash seeds. The script is 12.12's, with `data.WORKBOOK` in place
   of `data.FIXTURE`:

   ```bash
   cd extended_bridging
   S='import sys, data, derive
   df, audit = data.load(data.WORKBOOK)
   df = derive.derive(df, audit)
   derive.derive_cohort(df, audit)
   sys.stdout.write(audit.to_markdown())'
   PYTHONHASHSEED=0 uv run python -c "$S" > /tmp/a.md
   PYTHONHASHSEED=1 uv run python -c "$S" > /tmp/b.md
   diff /tmp/a.md /tmp/b.md
   PYTHONHASHSEED=0 uv run python -c "$S
   audit.write()"
   git check-ignore -v out/logs/audit_v7_july26_with_abs_contra_indication.md
   ```

   `diff` reports nothing and `check-ignore` reports a match. Note that `derive_cohort` is called here
   on the **full 126-record frame**, not the cohort — this is a reproduction check, not an analysis, so
   its median is 6.0 and not the 5.0 of §18. The analysis path runs `derive_cohort` from Stage 5, after
   the [§3] restrictions, and 12.15 is what checks the cohort median. The log names patients and must
   never be committed.
5. Acceptance 12.7 has been **seen to fail** with `.mask(src.isna())` removed — specifically its
   `float64` half, per §5.2; the `Int64` half stays green because the dtype propagates, and watching only
   that half is the trap this item exists to close. Then seen to pass when the mask is restored. It is the
   roadmap's acceptance criterion and a green test that has never been watched failing is not evidence.
6. Acceptance 12.3 has been **seen to fail** with the `unknown` mask removed from `_onset_type` — and
   seen to fail *differently*, yielding a fabricated `wake_up` rather than a fabricated `witnessed`, with
   `.fillna(False)` also removed (§4.2's square). And 12.10's guard clause seen to fail when
   `derive_cohort` is called twice without it, and 12.10's missing-core cell seen to fail with
   `.mask(core.isna())` removed — that one fails on the first run, because `core_ml` is `float64`.
7. Acceptance 12.4 and 12.14 have been **seen to fail** against a pasted level literal and a pasted
   `assert True` respectively. Both scans pass trivially if they match nothing.
8. The audit log has been read end to end by a human, and the three [§13]-relevant counts in it — the
   onset levels, the unknown-onset count, the median split — agree with §18.
9. No `.py` file mentions `TARGET_MISMATCH` or `MISMATCH_RATIO_UNDEFINED_COUNTS_AS_MET` at all —
   `grep -rn` over `extended_bridging/*.py` returns nothing. Under `specs/`, every surviving mention is a
   record of the withdrawal or of an edit that performs it (§1, §6.1, §8.4, §10, §15, §16, §18, T2); none
   describes either constant as live, and none instructs an implementer to read one.

## 18. Verification record

Checked on 2026-08-10, before this spec was finalised, by a read-only probe reporting aggregate
counts only — no outcome by arm, in the posture Stage 0 used. Recorded so a reader can tell which
numbers were verified rather than carried, and so the checks are re-runnable after a workbook update.
Acceptance 12.15 re-checks them in code.

| Claim | Where used | Verified |
|---|---|---|
| `wake_up` and `unwitnessed` are never both 1, and neither is ever missing, across all 126 records | §4.1, §4.3, §13 | yes |
| Onset levels over 126 records: witnessed 33, unwitnessed 36, wake-up 57 — so the three partition the frame | §4.2, §12.2, §12.15 | yes |
| `unknown_onset` is 1 on 93 of 126 records | §6.2, §12.15 | yes |
| Primary cohort — zero-bridging centre dropped, `ivt_contraindicated == 1` dropped — is 93 records, matching Stage 0's flow table | §0.1, §6.3 | yes |
| `core_ml` median: **5.0 mL** over the 93-record primary cohort, **6.0 mL** over all 126 | §0.1, §6.3, §6.4 | yes |
| Primary-cohort split at 5.0 mL: 44 above, **3 exactly at the median**, 45 below, 1 missing | §6.3, §12.10, §13 | yes |
| `core_ml` is exactly 0 on 39 of the 93 cohort records (50 of 125 non-missing over the full frame, per Stage 0) | §6.1, §6.3, §13 | yes |
| Stage 0's dichotomy counts over the full frame: mRS 0–2 53 events, mRS 0–1 31, death 37, mRS 5–6 42, each with 2 missing | §5.2, §12.15 | carried from Stage 0 |
| No derived name collides with an analysis name | §8.4, §12.13 | yes |

**pandas semantics, checked by running them** on pandas 2.3.3 as pinned by `uv.lock`. §3.1 is the
declaration; these are the observations behind it. Re-check them on any pandas major bump, because two
of them decide whether a mask is redundant or load-bearing.

| Claim | Where used | Verified |
|---|---|---|
| `Int64` `<= 2` yields `<NA>` where the value is missing — it propagates, it does **not** return `False` | §3.1, §5.2, §12.7, DoD 5 | yes, run |
| `float64` `> x` on `NaN` yields `False` — it absorbs | §3.1, §6.3, §12.10 | yes, run |
| `string` `!=` on `<NA>` yields `<NA>` | §3.1, §6.2 | yes, run |
| `Series.mask(cond)` with `<NA>` in `cond` treats it as **True** and applies the replacement — so an unfilled `(flag == 1)` assigns the *level*, not the baseline | §3.1, §4.1, §4.2, §12.3 | yes, run |
| `BooleanDtype.sum()` skips `<NA>`, so §4.1's count is right with or without `.fillna(False)` | §4.1 | yes, run |
| `Series.median()` skips missing by default and returns `NaN` when every value is missing | §6.3, §13 | yes, run |
| `nunique(dropna=True)` is 0 on an all-missing column, so `<= 1` catches it | §7.1, §12.11 | yes, run |

**Claims checked against the committed repository**, so they are re-checkable with no `data/`:

| Claim | Where used | Verified |
|---|---|---|
| `test_data.py` asserts the `penumbra_recomputed` entry's `n`, `case_ids` and table, and 12.16 its `(kind, step)` — none asserts the crossing sentence | §6.1, §10, T2 | yes, read |
| `data.py`'s only readers of `TARGET_MISMATCH` are `_correct`'s `threshold` and `crossing`, both inside the one entry | §6.1, §10 | yes, read |
| `config.py`'s `DERIVED_NAMES` is the literal `("onset_type",)` | §8.4 | yes, read |
| **[REV]** `DERIVED_NAMES` is read in **three** places in `test_config.py` — 9.3's `test_every_ps_covariate_resolves`, 9.8's tuple list, and 9.13's `_RESOLVABLE`. An earlier draft of this spec recorded "no module reads it yet", which was wrong and hid the 9.3 weakening | §8.4, §12.13 | yes, read |
| `test_config.py` 9.13's `parametrize` union and 9.8's tuple list both name their constants explicitly, so a new constant is outside them until added by hand | §12.13, T1 | yes, read |
| `test_data.py`'s heading test asserts the rendered text between the `observation` and `structural` headings is exactly `_none_`, so a kind inserted between them fails it | §9.2, §10, T3 | yes, read |
| `hand_frame()` casts every column carrying a declared dtype, so its `mrs_90d` is `Int64` and its `core_ml` is `float64`; `hand_frame(**overrides)` assigns a scalar to a whole column, and `corrupt()` is the per-record tool | §5.2, §12 | yes, read |
| `hand_frame()`'s `core_ml` is `(13, 20, None, 5, 0, 30)`: median 13.0, with 2 above, 1 exactly at (HAND-1), 2 below and 1 missing (HAND-3) | §12, §12.10 | yes, computed |
| **[REV]** `hand_frame()` has `mrs_90d = 2` on every record, so the four dichotomies are `mrs_0_2_90d` 1 and the other three 0 — **not** uniformly events, which an earlier revision of §12 stated | §12, §12.5 | yes, run |
| **[REV]** `constant_covariates` on `hand_frame()` returns **10 of 13** names, not `()`: everything in `PS_COVARIATES_FULL` except `core_ml`, `tmax6_ml` and `center`, because `_HAND_CONSTANT` single-values the rest. `onset_type` is among the ten, so the one-level-factor case needs no construction | §12.11 | yes, run |
| `constant_covariates` on the fixture returns `()` — each of its 13 covariates has 2 distinct values over its 2 records — which is what makes 12.12's "present with `n = 0` and no table" true of that frame | §12.12, §9.3 | yes, run |
| `tests/fixture_schema.xlsx` passes `derive` and `derive_cohort`: 2 records, `wake_up` 1 and 0, `mrs_90d` `<NA>` on one, `core_ml` `int64` 0 and 12, median 6.0, neither onset flag missing | §12.12 | yes, run |
| `pilots/data.py:141` builds `onset_type` with `np.select(..., "witnessed")` and wraps it in `pd.Categorical`; `:144` uses `.astype(int)` | §4.3, §4.4, §16 | yes, read |
| `stage0_data_inventory.py:97`'s `derive` already carries the `.mask(source.isna())` idiom and its justification | §5.1, §16 | yes, read |
| `test_config.py` 9.4 walks `extended_bridging/**/*.py` minus three exemptions and pins that set, so `derive.py` is scanned from the moment it exists | §3, §12.13, DoD 3 | yes, read |

## 19. What this spec changed elsewhere

Two entries, both landing with this document rather than with the implementation, because a
specification that contradicts the SAP is worse than no specification.

| # | Change | Files | The question it settles |
|---|---|---|---|
| 1 | The [§13] target-mismatch subgroup is withdrawn | `statistical_analysis_plan.md` §13 (dated amendment), `implementation_roadmap.md` Stage 3, and §6.1 here | Whether a CTP-selected cohort can be subgrouped on close to its own selection criterion. It cannot usefully, and the boundary-operator question the criterion needed is retired rather than answered |
| 2 | Roadmap Stage 3 gains its `**Spec:**` line, the frozen-median rule, and the detect-don't-drop reading of zero variance | `implementation_roadmap.md` | Whether the median is recomputed per bootstrap replicate (no), and which of Stages 3 and 6 deletes a constant covariate (Stage 6) |

Everything else this spec asks for is an amendment to Stage 1 or Stage 2 in service of Stage 3, and
§10 is the ledger for those. None of them changes a statistical decision, a column's fate, or a
number.

## 20. Revision record — the engineering review of 2026-08-10

Split deliberately into **defects** and **decisions**. A defect is a sentence that was false about this
repository or about pandas; a decision is a choice the review put to the author and the author made.
The distinction matters because a later reader who disagrees with a decision may revisit it, whereas a
defect has no other side.

**Defects found and corrected.** Each was verified by running the code or reading the named file, not
by inspection.

| # | What the draft said | What is true | Corrected in |
|---|---|---|---|
| D1 | "`NA <= 2` evaluates to `False`, not to `NA`, so without the mask the two records with no 90-day mRS become non-events" | `mrs_90d` is `Int64`, which **propagates** `<NA>`. `.mask(src.isna())` is a no-op on every frame this repository builds, so Definition of done 5 — "seen to fail with the mask removed" — was **impossible to satisfy**, and an implementer who deleted the mask would have watched the test stay green | §3.1, §5.2, §12.7, DoD 5 |
| D2 | "a `BooleanDtype` mask carrying `<NA>` is not a mask" | It is one, and `<NA>` is treated as **True**: an unfilled `(flag == 1)` assigns the *level*. So the failure the missing `.fillna(False)` causes is a fabricated `wake_up`, not the fabricated `witnessed` §4.3 describes — a worse error, and one no assertion covered | §3.1, §4.1, §4.2, §12.3 |
| D3 | §10's ledger listed only "12.16's heading list gains `Derivations`" | `test_all_six_headings_appear_and_an_empty_section_says_none` asserts the rendered text between the `observation` and `structural` headings is exactly `_none_`. Inserting `derivation` between them **fails a currently green test**, and no ledger row or task verify-line named it | §9.2, §10, T3 |
| D4 | §18 recorded, as verified: "`DERIVED_NAMES` … no module reads it yet" | It is read three times in `test_config.py`. The consequence that was missed: widening it **weakens** 9.3's `test_every_ps_covariate_resolves`, so the disjointness assertion is what restores that guarantee rather than duplicating it. Separately, 9.13's union and 9.8's list name their constants explicitly, so all four new registries would have shipped outside both sweeps | §12.13, §18, §10, T1 |
| D5 | §8.4 gave the computed `DERIVED_NAMES` with no placement | It depends on `OUTCOMES` and `SUBGROUPS`, both declared *below* `DERIVED_NAMES`' current position. Written as specified, `config.py` raises `NameError` on import | §8.4, §10, T1 |
| D6 | §12: the frame is "`hand_frame()` extended with the columns Stage 3 reads" | Nothing needs extending — every column is already there. What was needed is per-record *variation*, and `hand_frame(**overrides)` cannot give it (it assigns a scalar to a whole column); `corrupt()` is the tool. The draft named the wrong helper for its own acceptance tests | §12 |
| D7 | 12.12 and DoD 4 gave the reproduction driver as `"…load, derive, derive_cohort…"` | An ellipsis in a sole-source spec is a decision made at the keyboard, and a freely chosen driver can render no derivation at all and still compare equal | §12.12, DoD 4 |
| D8 | 12.11: "Returns `()` on the hand frame" | It returns **10 of 13**. `_HAND_CONSTANT` single-values ten of `PS_COVARIATES_FULL`, so the test as written fails on its first run. The repair is better than the original: the pinned ten-name tuple is the only assertion in §12 that would catch `constant_covariates` ranging over a set instead of the caller's sequence, and it covers the one-level-factor case for free through `onset_type` | §12.11, §18 |
| D9 | An earlier pass of *this revision* wrote "every dichotomy is uniformly an event until a test varies one" | With `mrs_90d = 2`, only `mrs_0_2_90d` is 1; `mrs_0_1_90d`, `death_90d` and `mrs_5_6_90d` are all 0. Recorded because the mixed pattern is what makes the frame usable for 12.5 at all — a frame where all four read alike could not distinguish a registry-driven loop from a constant | §12, §18 |

**Decisions taken in review.** Options were weighed; these are the ones chosen.

| # | Question | Chosen |
|---|---|---|
| R1 | Where does the dtype fact live, given three sections state it? | Declared **once** in §3.1 with the verified table, and each site says which of the two cases it is. Three separate statements is how D1 and D2 became inconsistent |
| R2 | Keep the two redundant `.fillna(False)` calls, or delete them? | **Keep**, with the honest rationale (they make the guarantee independent of a dtype no derivation controls) and with §12.3 pinning all four cells of §4.2's square. Deleting them would have removed the guard that makes §4.2's loop safe if §4.1's assertion is ever dropped |
| R3 | Fix the heading test by re-pointing it at `derivation`? | **No** — take the neighbour from `KINDS`. Naming the new neighbour is the same defect one position along, and rots on the next kind |
| R4 | Should `_dichotomies` copy, to make §0.2's promise unconditional? | **No** — keep the in-place write and state the contract instead. Three privates copying a 126-row frame is ceremony; the caller obligation is written down in §0.2 and §5.1 rather than inferred |
| R5 | Is the all-missing-`core_ml` path worth a guard? | **No guard, recorded in §13.** It degrades to an empty subgroup rather than a wrong one, and the log prints the median as `missing` — visible where it matters |
| R6 | Should the Stage 1 dtype declaration be folded into Stage 3? | **No — deferred to `TODOS.md`, §15 carries the row.** Doing it now would invalidate §18's dtype rows and §12.10/§12.12 mid-implementation. §13 carries the `Float64`-versus-`float64` caveat for whoever picks it up |

**Alternatives that were considered and rejected.** Recorded because each is the obvious thing to
propose on reading the spec, and a reader who proposes one deserves the reason it was already declined
rather than a re-run of the argument.

| Instead of | Rejected because |
|---|---|
| Deleting `.mask(src.isna())` and pinning the source dtype with an assertion (the smallest code) | Invariant 6 would then rest on `READ_DTYPES` never changing *and* on every frame ever passed to `_dichotomies` being cast. Hand frames are trivially built uncast — a dict and a `DataFrame` call — and the failure is a fabricated non-event, which is the exact silent error §5.2 exists to prevent |
| Moving every derived view to the end of `config.py` as one block (§8.4) | A larger diff than the change needs, and it moves lines the Stage 1 spec documents by position. Pinning each constant next to the registry it reads keeps the `POST_TIME_ZERO` precedent |
| Re-pointing the heading test at `derivation`, the new neighbour (§9.2) | The same defect one position along: it rots on the next kind inserted. Taking the neighbour from `KINDS` is the fix that cannot recur |
| A purpose-built six-record frame in `test_derive.py` instead of composing `test_data`'s helpers (§12) | "Two hand frames drift" — §12's own argument. It would also duplicate a contract-valid record set that already encodes the missing-volume, both-arm, tie and copy-paste cases |
| Having `_dichotomies` return only its four new columns, so it touches nothing (§5.1) | The largest change to a snippet Stage 1 §3.2 already pins as the reproduced form. The in-place write with a stated caller contract costs one sentence and no code |
| Testing byte identity within one process only (§12.12) | Stage 3 adds two dict-ordered loops, which is exactly the shape that has broken this before; a same-process comparison cannot see iteration order at all. The cross-process test is cheap and this repository has already been bitten once |
| Building the hand frame so `constant_covariates` returns `()` (§12.11) | Unreachable through `hand_frame(**overrides)`, which makes columns *more* constant, and ten `corrupt` calls would be a second hand frame by the back door. The pinned ten is the stronger assertion |

**What the review did not change.** The architecture: two entry points, the frozen median, detect-don't-drop
for zero variance, the [§13] withdrawal, the registry-driven loops, and every number in §18's first table.
The review challenged the two-entry-point split and the audit-log rename and sustained both — §0.1's
6.0-versus-5.0 mL argument is the reason the split is right, and the rename is three string literals for a
log that now spans stages.

**What this review did not do, so it is not mistaken for having been done.** A cross-model challenge was
attempted and did not complete — the second reviewer authenticated, began, and failed partway with repeated
`401 Unauthorized`. So every finding above is single-reviewer, and the sections it did not reach are §6.1's
statistical argument for the [§13] withdrawal (which stands on the SAP amendment, not on this spec) and §9.3's
`detail` wording. Nothing in §20 depends on the outside voice; the two lists above are what one reviewer
verified by running the code and reading the tests. A later cross-model pass is worth having and would start
with §9.3 and §11.
