# Stage 7 spec — balance and overlap diagnostics

Implements roadmap Stage 7 [§9]. Section references in brackets are to `statistical_analysis_plan.md`.
Numbers and decisions referenced as DECISION *n* are established in Stage 0 and recorded in
`../out/stage0_data_inventory.md`. Constants named in `SMALL_CAPS` are Stage 1's and live in
`config.py`; `stage1_config_and_data_contract.md` is their specification. The frame this stage
receives is specified by `stage5_cohort_construction.md` §11 and the result object by
`stage6_propensity_and_weights.md` §11.

**Status.** Written 2026-08-18 against the landed Stages 1-6 and the workbook as verified in §18,
and **revised the same day after one independent review round, recorded in §20.** Every number below
was produced by running code; none is carried. It is the **sole source for the Stage 7
implementation**: everything the implementer needs is here, and anything not here is not to be
invented. Where the review changed a decision rather than a number, §20 names the decision and this
document carries only the outcome — there is no second version to reconcile.

**Goal.** Standardised mean differences before and after weighting exist for every [§6] confounder
and every negative control, against one unweighted pooled standard deviation that does not move
between the two columns; overlap is reported within each centre as well as pooled; a centre with
structural non-positivity is named as such rather than given an overlap plot; and every quantity
that is undefined is *reported* as undefined rather than as zero, because the difference between
those two is the difference between a diagnostic and a reassurance.

**Not in scope.** Any estimate of the treatment effect (Stages 8 and 9), the bootstrap (Stage 10),
the [§13] subgroups and sensitivity specifications (Stage 11), the [§14] populations (Stages 12 and
13), and every figure (Stage 14 [§16]). Stage 7 **fits nothing, refits nothing, and re-weights
nothing.** It adds no column to the cohort frame and edits no value; it returns a result object and
appends two audit entries.

**One thing this spec settles that no earlier stage could.** Balance is judged **per declared factor
level, not per design-matrix column** — and the two are different sets. `model.design` drops each
factor's reference dummy by name and then drops constant columns (Stage 6 §4.2), so a balance table
built from its columns silently omits `center_HUG`, `onset_type_witnessed` and `center_USZ`: three of
the seven indicators the [§6] factors expand to on this cohort, including the baseline every centre
coefficient is measured against. That is Stage 6 §7.5's finding in its second form, and this is the
stage where it decides a reported number rather than a table's completeness. §4.2 is the rule and
Stage 7 therefore does not call `model.design` at all.

**And one finding that arrives with the stage rather than being designed into it.** Stage 6 §13 left
the size of Firth's residual imbalance open and quoted the pilots' synthetic measurement, `0 < worst
|SMD| < 0.05`, as what to expect — "inside [§9]'s 0.10 but not zero". **On this cohort it is 0.211,
on `center = HUG`, and a second in-model row sits at 0.154** (§18). [§9]'s threshold is exceeded, by
covariates that are *in* the propensity model, and §6.4 measures that this is a property of the
prescribed estimator rather than a defect: with an unpenalised score every in-model row balances to
1.03e-15, and the [§7] Firth score does not. This is a finding for the PI under [§9] and [§13]. It
is **not** a licence to change the estimator, and §14 says so.

---

## 0. Where Stage 7 sits

```
  cohort.build(df, audit)  →  93 rows × 33 columns                   [Stage 5 §11]
  propensity.fit(cohort, audit)  →  Propensity(e, w, in_model, …)    [Stage 6 §11]
                       │
                       ▼
  ┌────────────────────────────────────────────────────────────────────────────┐
  │  STAGE 7 — balance.py     reads config, data.Audit, propensity              │
  │                           and NOT model — §4.2                              │
  │                                                                             │
  │   smd(x, a, w)                  [§9] one standardised mean difference       │
  │     ├─ weighted means, UNWEIGHTED pooled SD          (§5.1)                 │
  │     └─ four undefined branches, each its own         (§5.3)                 │
  │                                                                             │
  │   assess(df, ps, audit) -> Balance                                          │
  │     ├─ _assert_balance_inputs   B1,B2 THEN B3…B5 → SchemaError  (§4.5)      │
  │     ├─ the balance set, expanded to DECLARED LEVELS  (§4.2)                 │
  │     ├─ SMD before and after, over in_model's 92      (§4.4, §5)             │
  │     ├─ the within-centre overlap table + pooled row  (§6)                   │
  │     └─ two `model` audit entries                     (§7)                   │
  └────────────────────────────────────────────────────────────────────────────┘
                       │
                       ▼
   Balance(covariates, centres, pooled)  →  Stage 11 [§13], Stage 14 [§16]
   19 covariate rows over 14 declared names; `centres` is exactly CENTER_ORDER — 4 rows —
   and the `all (pooled)` row is its OWN field beside them, never a fifth centre  (§3, §6.1)
```

`assess` appends two entries of the existing `model` kind. It writes no file, as no stage before it
does.

### 0.1 Why a module of its own

**Stage 7 is the first stage that returns a verdict.** Stages 1-5 classify and restrict, Stage 6
fits; every one of them produces something a later stage consumes as input. Stage 7 produces
something only a *reader* consumes: [§9]'s threshold, applied to [§6]'s covariates, reported for
[§16]. Nothing downstream of it computes with its output, and that is why it is a module and not a
function on `propensity.py`.

Three consequences, all of them the point:

- **`balance.py` does not import `model`.** §4.2 is the reason: the design matrix is the wrong
  object to judge balance on, and a module that imports `design` will eventually call it. What it
  needs instead — indicators over `FACTOR_LEVELS` — is four lines, and they range over the
  declaration exactly as `propensity._design_table` does (Stage 6 §7.5).
- **`smd` is public, and it is Stage 12's as much as this stage's.** [§14a]'s support check compares
  IVT-treated patients with never-IVT-centre patients and asks for a baseline table; that is a
  standardised mean difference over a different population with unit weights. A second
  implementation there would be a second definition of the yardstick, and §9 hands it over
  explicitly — the same move Stage 6 §3 made with `ess`.
- **`assess` takes the `Propensity`, not the cohort alone.** It does not refit and cannot: there is
  no code path in this module that fits anything, and §12.10's scan asserts that `model.firth` and
  `propensity.fit` are named nowhere in it.

### 0.2 What Stage 7 does not touch

`assess` takes the cohort frame and the `Propensity` and gives back neither: **no column is added,
no value is edited, no row is removed, and no weight is recomputed.** Stage 6 §0.2 made the same
promise for the same reason one stage earlier, and here it is stronger, because the temptation is
different: the natural way to write a balance table is to attach `w` to the frame and group. A frame
carrying the point fit's weights, resampled into a Stage 10 replicate, is a replicate weighted by the
wrong score.

**No second effective sample size exists anywhere in this module.** Stage 6 §11's handover ends
"It does not refit, does not re-weight, and does not compute a second effective sample size", and
§6.3 is how that is honoured: `propensity.ess` is called per centre per arm, and a Kish sum appears
nowhere in `balance.py`.

## 1. Deliverables

| Path | Contents |
|---|---|
| `extended_bridging/balance.py` | the two public names `smd` and `assess`; the three types `Balance`, `CovariateBalance` (carrying `sd`) and `CentreOverlap`; and the privates §3 lists — including `_pooled_sd` and `_ratio`, which §5.5's one-yardstick rule requires, and `_worst` / `_over_threshold`, which the `Balance` methods delegate to |
| `extended_bridging/tests/test_balance.py` | the acceptance tests of §12, its **own** module-scoped `workbook` fixture returning `(df, ps, audit)` from **one** linear run against **one** `Audit` (§12.0.3), the hand-computed SMD vectors of §12.3 **written out there in full**, and §12.0.2's golden balance vector |
| `extended_bridging/config.py` | **amended** — two computed views, `NEGATIVE_CONTROLS` and `BALANCE_SET`, in the covariates block, and nothing else (§4.3, §10) |
| `extended_bridging/tests/test_config.py` | **amended** — four assertions pinning what those two views are computed from (§4.3, §10) |
| `extended_bridging/tests/reference/balance_psweight.R` | ~25 committed lines: read a CSV and a score, ask `PSweight::SumStat` for its balance table, write it back. No covariate list and no analysis logic (§16b) |
| `extended_bridging/tests/test_reference_r.py` | **amended** — the Stage 7 half of the R oracle behind the gate that already exists, plus `_r_user_library`'s candidate search and the assertion that the gate is measured OPEN rather than assumed open (§16b, T6) |
| `extended_bridging/implementation_roadmap.md` | **amended, and landing with this document** — Stage 7 gains its `**Spec:**` line, three **Accept when** items and one correction (§19) |

`out/logs/audit_<label>.md` is an **output, not a deliverable**, as in Stages 2-6: gitignored,
written only when a caller asks.

**Nothing under `specs/` may quote a case identifier, and nothing here does.** Every patient below is
a count, and every named record is a `HAND-N` or `COHORT-N` fixture record.

## 2. Environment

`uv`, Python 3.12, flat module layout, `import config as C`. Commands run from `extended_bridging/`.

```bash
cd extended_bridging && uv sync && uv run pytest -v
```

**No new dependency, and none is available.** `numpy` and `pandas` carry everything; the balance
tooling that would otherwise be worth importing — `cobalt`, `WeightIt` — is R-only, which Stage 6
§16b measured and §16 here does not re-litigate. The one thing R buys this stage is an oracle, and
§16b is that.

**What this costs to run.** Measured (§18): the 19-row pooled table takes **7.4 ms**, and the whole
stage — one pooled table plus one per retained centre — is about **30 ms**. That is three orders of
magnitude more than the fit it diagnoses (Stage 6 §2: tens of microseconds), and it does not matter
in the slightest, because **Stage 7 is not on the bootstrap path.** [§10] refits the propensity model
and the outcome regression in every replicate; it does not recompute balance, because balance is a
property of the point specification and [§13] asks for "ESS and worst residual |SMD|" once per
specification rather than once per replicate. This paragraph exists so that no later stage reaches
for a vectorised rewrite of a Python loop over nineteen rows.

## 3. Module shape

```python
# balance.py
from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass
from typing import Final

import numpy as np
import pandas as pd

import config as C
import propensity
from data import Audit, _fmt
```

`Final` is imported because §7.4 declares three module constants with it. `from __future__ import
annotations` makes a module-level annotation a string that is never evaluated, so omitting the import
would run — measured — and fail nothing but a type check. It is the one import in this block whose
absence is invisible at runtime, which is exactly why it is named here.

`_fmt` is imported for the reason `derive.py` and `propensity.py` give: it is the pipeline's *one*
float formatter, and a second one is a second way for two runs to disagree.

**Three types, and the third is what makes the first two safe to read.**

```python
@dataclass(frozen=True)
class CovariateBalance:
    """One row of the [§9] balance table: one declared level, or one linear covariate."""

    covariate: str              # "age", or "center = HUG" for a declared factor level
    role: str                   # §4.3 — computed, never declared
    n: int                      # records this covariate is present on, within in_model [§11]
    sd: float                   # the UNWEIGHTED pooled SD both columns are divided by (§5.1)
    unweighted: float           # SMD before weighting; nan where undefined (§5.3)
    weighted: float             # SMD after weighting;  nan where undefined (§5.3)
```

**`sd` is a field and not a local, and it earns its place three times.** [§9]'s denominator is the one
quantity in this stage that two rows can disagree about, and until it was a field nothing could read
it: §6.1's within-centre column needs it as a *yardstick from another population* (§5.5), §16b's
oracle comparison reconstructs `|ours| × our_pooled_unweighted_SD` and had to recompute it, and
§12.3's "the yardstick does not move when the weights do" was an algebraic reconstruction from two
returned ratios rather than an assertion on the number itself. One field turns all three into reads.

```python
@dataclass(frozen=True)
class CentreOverlap:
    """One row of the [§9] within-centre overlap table, or the pooled row."""

    centre: str
    in_cohort: int              # records at this centre in the [§3] cohort
    weighted: int               # of those, records inside in_model — NOT the same number (§6.1)
    n: dict[int, int]           # keyed by arm code, as TREATMENT_LABELS is
    e_range: dict[int, tuple[float, float]]
    ess: dict[int, float]       # propensity.ess per arm; nan where an arm is empty (§6.3)
    max_weight: float
    weight_share: float
    worst: float                # worst within-centre |SMD|; nan where none is defined
    status: str                 # §6.2's THREE, verbatim — _REPORTED | _STRUCTURAL | _NO_CONTRAST
```

**The `status` values are §6.2's three and nothing else**, written here as the constant names rather
than as prose so that this comment cannot drift from §7.4's literals. An earlier draft wrote
`"reported"` and `"structural non-positivity"`, and neither is what §7.4 declares — the strings are
`"overlap reported"` and `"structural non-positivity [§3]"`. Stage 14 reads this column by value
before it draws anything (§6.2, §9), so a docstring naming a string the code never produces is a
docstring that tells the next reader to match on the wrong one. **The pooled row is not a fourth
status**: it carries `_REPORTED` like any other both-arm row, and what makes it the pooled row is
that it lives in `Balance.pooled` rather than in `Balance.centres`.

The two verdicts are **module-level functions over a row sequence, and the dataclass delegates**:

```python
def _worst(rows: Sequence[CovariateBalance]) -> tuple[CovariateBalance | None, tuple[str, ...]]:
    """The largest |weighted| SMD among the DEFINED rows, AND the undefined ones' names.

    Both, from one call, and that is the whole design of this pair [§3.1]. `max` over a column
    containing `nan` skips it in pandas and poisons it in numpy, so a caller asking only for the
    worst gets either a number that ignores the rows nobody could judge or a `nan` that hides the
    rows anybody could. Returning the pair makes it impossible to report the first without being
    handed the second.
    """
    defined = [row for row in rows if np.isfinite(row.weighted)]
    undefined = tuple(row.covariate for row in rows if not np.isfinite(row.weighted))
    return (max(defined, key=lambda row: abs(row.weighted)) if defined else None), undefined


def _over_threshold(rows: Sequence[CovariateBalance]) -> tuple[str, ...]:
    """Covariates whose |weighted| SMD reaches SMD_THRESHOLD. Undefined rows are NOT here.

    The undefined rows are `_worst`'s second element, so the two collections partition the table
    with the defined-and-balanced ones and no row is in neither by accident.
    """
    return tuple(row.covariate for row in rows
                 if np.isfinite(row.weighted) and abs(row.weighted) >= C.SMD_THRESHOLD)
```

```python
@dataclass(frozen=True)
class Balance:
    """What Stage 7 returns. The frame and the Propensity come back untouched [§0.2]."""

    covariates: tuple[CovariateBalance, ...]
    centres: tuple[CentreOverlap, ...]   # exactly CENTER_ORDER, in CENTER_ORDER's order
    pooled: CentreOverlap                # the `all (pooled)` row — its own field, never a centre

    def worst(self) -> tuple[CovariateBalance | None, tuple[str, ...]]:
        """§3's pair, over this table's rows. The public name; `_worst` is the definition."""
        return _worst(self.covariates)

    def unbalanced(self) -> tuple[str, ...]:
        """The rows reaching SMD_THRESHOLD, over this table's rows. Undefined rows are not here."""
        return _over_threshold(self.covariates)
```

**Why the split, because it is not decoration.** §7's `_smd_detail` needs both verdicts to write its
`detail` string, and it has only the rows — the centres do not exist yet when the first audit entry is
recorded (§8). An earlier draft reached them by constructing a throwaway `Balance(tuple(rows), ())`,
which broke the moment `Balance` gained a third field, and — worse while it worked — it wrote the
threshold predicate a **second** time, inline, beside `unbalanced()`'s. Two spellings of "which rows
reach 0.10" in one module is one edit away from a `detail` string that names four rows above a table
that flags five. One definition, two callers, no temporary object.

**`pooled` is a field and not the last element of `centres`, and that is a decision.** An earlier
draft returned `centres + (pooled,)` as one tuple, which cost three things at once: `len(centres)`
stopped being `len(CENTER_ORDER)`, so no consumer could iterate the declared centres without knowing
to drop a row; the pooled row needed a status of its own — a fourth literal, `"pooled"`, that §6.2's
three-status table never declared and Stage 14 would have had to match on; and a caller counting
reported centres off `centres` would have counted the whole cohort as one of them. Measured on the
fixture, the rendered table's last row read `pooled` in the `status` column with nothing in §6.2
explaining it.

The split costs nothing, because **the rendered table is still one table**: §8 passes
`centres + (pooled,)` to `_overlap_table`, so [§9]'s "within each centre as well as pooled" is one
grid computed by one function (§6.1), and only the *return value* distinguishes the two kinds of row.
That is the right place for the distinction: the log has a `centre` column reading `all (pooled)`,
and a Python caller has a field name.

Public surface, and it is two names:

```python
# balance.py
def smd(x: pd.Series | np.ndarray, a: pd.Series | np.ndarray,
        w: pd.Series | np.ndarray) -> float: ...
def assess(df: pd.DataFrame, ps: propensity.Propensity, audit: Audit) -> Balance: ...
```

**Every `python` fence in this document is valid Python**, including the signature listings above —
hence the `: ...` bodies. Stage 6 §3 states the rule and §18c is the check that earned it; three
fences failed it there before the rule existed.

Privates are `_assert_balance_inputs`, `_levels`, `_role`, `_table` (the covariate rows), `_centre`,
the two verdicts `_worst` and `_over_threshold` that `Balance`'s methods delegate to, and §7's
`_status`, `_verdict`, `_range`, `_smd_detail`, `_smd_table`, `_overlap_detail` and `_overlap_table`.
There **is** a `_pooled_sd`, and §5.5 is why — an earlier draft of this
document forbade one, on the argument that a denominator computed anywhere else can be computed over a
different population than the means it divides. That argument is right about the risk and wrong about
the remedy: §6.1's within-centre column has to divide a centre's means by the **whole** weighted set's
SD, so the separation exists whether or not a function name admits it. Hiding it inside `smd` did not
prevent it; it only made it the one thing in this stage nobody could assert. So the denominator is one
function, applied to a foreign population in exactly one place, and that place says so (§5.5, §6.1).

**This module is not exempt from the Stage 1 §7 raw-name scan and must never become exempt.** It
names no raw header.

### 3.1 The numerical facts this stage turns on

Declared once here, verified by running them on pandas 2.3.3 / numpy 1.26.4 as pinned by `uv.lock`
(§18). Every one of them is a way this stage returns a plausible wrong number.

```
  np.average(x, weights=w)  with w summing to zero
      →  ZeroDivisionError("Weights sum to zero, can't be normalized").  RAISES.
         So an arm whose weights all vanish must be checked BEFORE the mean, not after.

  np.array([3.0]).var(ddof=1)   →  nan, with a RuntimeWarning.  An arm of one has no
  np.array([]).var(ddof=1)      →  nan, with a RuntimeWarning.  variance and no pooled SD.

  pd.Series([nan, 0.2]).abs().max()  →  0.2     pandas SKIPS the missing value
  np.abs(np.array([nan, 0.2])).max() →  nan     numpy POISONS the result
      Two libraries, two opposite wrong answers, from the same expression. §3's `worst()`
      returns the undefined names beside the maximum so that neither can be reported alone.

  np.average(np.full(3, 14.0), weights=w1)  →  14.0
  np.average(np.full(2, 14.0), weights=w0)  →  14.000000000000002
      The SAME constant, two weight vectors, differing by 1.776e-15. A zero-variance branch
      that compares the two WEIGHTED means for equality calls a covariate identical in both
      arms "undefined". Measured on the committed fixture, where it did exactly that. §5.3a.

  np.average(np.full(n, 1.0), weights=w)  →  exactly 1.0, always
  np.average(np.full(n, 0.0), weights=w)  →  exactly 0.0, always
      Indicators are immune to the fact above, which is why it hid: 17 of the fixture's 19
      rows are indicators or genuinely varying, and only the one integer-valued constant bit.
```

**The first fact and the second compose into this stage's F4.** An arm with no positive weight and an
arm with fewer than two records are different failures — the first raises inside `np.average`, the
second returns `nan` from `var` — and both mean "no standardised mean difference exists here". §5.3
gives each its own branch and its own return, because a `ZeroDivisionError` escaping `assess` would
stop a diagnostic mid-table, and [§9] wants the table.

## 4. The balance set [§6, §9]

### 4.1 The rule

```
  Balance is judged against the FULL [§6] confounder set, not only the covariates a given
  specification put in its propensity model [§9].

  The set is PS_COVARIATES + BALANCE_ONLY — 14 declared names — expanded to one row per
  DECLARED LEVEL of each factor and one row per linear covariate: 19 rows (§18).

  Every row carries a `role` computed from the declarations, never declared (§4.3).
```

The set does not change when the specification does. That is [§9]'s sentence and it is the whole
purpose of the `role` column: under [§13]'s full-covariate propensity model the four vascular risk
factors move from `negative control` to `propensity model` and **not one row moves in or out of the
table**, so the two specifications are read against the same fourteen names.

### 4.2 Declared levels, not design columns

```python
def _levels(df: pd.DataFrame, name: str) -> list[tuple[str, pd.Series]]:
    """One (label, indicator) pair per DECLARED level of a factor; one pair for a linear covariate.

    Ranges over FACTOR_LEVELS, never over the frame and never over a design matrix. A level nobody
    in the cohort has is a row reading 0, which is a fact about the population; a level that
    vanishes is a fact nobody sees [Stage 6 §7.5].

    The indicator carries the factor's own missingness — `mask(isna())` — rather than encoding an
    absent factor value as a zero in every level, which is Stage 6 §4.5 D4's failure in a table
    instead of in a design.
    """
    if name not in C.CATEGORICAL:
        return [(name, df[name].astype("Float64").astype(float))]
    return [(f"{name} = {level}",
             (df[name] == level).astype(float).mask(df[name].isna()))
            for level in C.FACTOR_LEVELS[name]]
```

**Why this is not `model.design`, stated at length because the alternative is one import away.**
`design` is built for a fit and does three things a balance table must not inherit:

- **It drops each factor's reference dummy by name** (Stage 6 §4.2). On the [§6] factors that is
  `center_HUG` and `onset_type_witnessed`. `center = HUG` is the **worst-imbalanced covariate in
  this cohort** — measured, |SMD| 1.343 before weighting and 0.211 after (§18) — so a table built on
  `design`'s columns would omit the row that carries the finding, and would report the worst residual
  imbalance among the covariates it did keep as though it were the worst overall.
- **It drops constant columns** (Stage 6 §4.3). On the cohort that is `center_USZ`, which is
  all-zero because [§3] restriction 1 removed the centre. A balance table with no `center = USZ` row
  cannot state that the restriction worked; with one, it reads 0 before and 0 after, which is what a
  declared level nobody has *should* read.
- **It raises D1-D4** on inputs a balance table has no reason to reject, and D3/D4 exist to protect a
  *fit* from an all-zero dummy row. Here the same value is a legitimate table cell.

So: three of the seven indicator rows the two [§6] factors expand to on this cohort are absent from
`design`'s output, and one of them is the finding. `balance.py` therefore does not import `model`, and
§12.10's scan asserts it — not because calling `design` would fail loudly, but because it would
succeed.

**The widths, measured, because §12.2 and DoD-6 assert one of them** (§18):

```
    declared levels over BALANCE_SET            19 rows   =  14 from PS_COVARIATES + 5
                                                             (7 linear + 3 onset_type + 4 center)
    design(sub, BALANCE_SET)[0].columns         16 cols   — omits EXACTLY the three indicators:
                                                             center_HUG, onset_type_witnessed,
                                                             center_USZ
    design(sub, PS_COVARIATES)[0].columns       11 cols   — omits those three AND all five
                                                             BALANCE_ONLY names: eight of nineteen
```

**§12.2's companion builds the `BALANCE_SET` form, and 16 is the number it asserts.** It is the
apples-to-apples comparison: the same fourteen declared names, differing from §4.2 only in how the
factors are encoded, so the entire gap is the three indicator rows this section is about. The
`PS_COVARIATES` form is what `pilots/analysis.py:249-265` actually does (§16) and it is worse by five
further rows, but a companion built on it would confound two omissions and let a reader attribute the
missing `center = HUG` to the covariate list rather than to the encoding.

An earlier draft of this document wrote **12** for both, in §4.2, §12.1, §12.2, §18 and DoD-6. Twelve
is `design_matrix`'s `n` — `X.shape[1] + 1`, the parameter count including the intercept
(`propensity.py:477`) — and it is the design matrix's *width plus one*, not any table's row count.
The 19-row total was right and its decomposition was written as 12 + 5, which is 17.

### 4.3 The `role` column is computed, and `BALANCE_ONLY` is not the negative-control set

`BALANCE_ONLY` holds **five** names; [§6] names **four** vascular risk factors as balance negative
controls. The fifth is `penumbra_ml`, which [§6] excludes because it is a deterministic function of
`core_ml` and `tmax6_ml` — including all three makes the design exactly singular — and which [§13]
gives a balance-table role of its own. It is in the table and it is **not** a negative control, and a
document or a table that equates the two sets is asserting that a collinear covariate is a covariate
weighting was never expected to fix.

The split is derived from declarations that already exist, never written out as a sixth literal:

```python
# config.py, in the covariates block — §10

# [§6] names four vascular risk factors as balance negative controls and adds them back in the
# [§13] full-covariate sensitivity propensity model. Those are the same four, so the set is
# COMPUTED from that identity rather than declared beside it: a fifth risk factor added to
# PS_COVARIATES_FULL becomes a negative control here in the same edit, and cannot fail to.
NEGATIVE_CONTROLS: Final[tuple[str, ...]] = tuple(
    c for c in PS_COVARIATES_FULL if c not in PS_COVARIATES)

# The [§9] balance set: the full [§6] confounder set plus everything BALANCE_ONLY carries, which
# is the four negative controls and penumbra_ml. NOT the propensity model's covariate list —
# [§9] judges balance against the full set regardless of what a specification fitted.
BALANCE_SET: Final[tuple[str, ...]] = PS_COVARIATES + BALANCE_ONLY
```

```python
def _role(name: str) -> str:
    """What this covariate is TO the specification being diagnosed. §4.3.

    Reads PS_COVARIATES for the same reason `propensity.fit` does [Stage 6 §6.4]: the estimand is
    indexed by the propensity model, so "in the model" is a property of the one prespecified
    specification and not an argument a caller may vary.
    """
    if name in C.PS_COVARIATES:
        return "propensity model"
    if name in C.NEGATIVE_CONTROLS:
        return "negative control"
    return "excluded [§6]"
```

Measured (§18): `NEGATIVE_CONTROLS` is `('hypertension', 'hyperlipidemia', 'diabetes', 'smoking')`
and `BALANCE_ONLY` minus it is `('penumbra_ml',)`, so the three roles partition the fourteen names
with nothing left over. §10 adds the assertions that keep it that way.

### 4.4 The population is `in_model`, in both columns

**Both the unweighted and the weighted column range over Stage 6's `in_model` mask — 92 of 93
records on v7 — and the pooled standard deviation is computed over the same 92.** Decided with the
PI, 2026-08-18.

[§9] asks for one yardstick: "common unweighted pooled SD, so the yardstick does not move". It says
that about the two *columns*, and the same argument settles the *rows*: a "before" column over 93
records and an "after" column over 92 are not two views of one contrast, they are two contrasts, and
the standard deviation dividing both would be a property of a population only one of them describes.
The weighted column has no choice — `w` is `nan` off `in_model` and Stage 6 §9 forbids filling it —
so the unweighted column follows it.

**What the decision costs, measured rather than asserted** (§18): the largest difference between an
unweighted SMD computed over the 93 and over the 92 is **0.0284**, on `onset_type = witnessed`, and
no row changes side of the 0.10 threshold. It is small, it is not zero, and it is recorded so that
nobody has to re-derive whether it mattered.

Two consequences the entry's `detail` must carry:

- **The [§3] cohort is 93 and this table describes 92.** That is Stage 6 §4.4's statement — the ATO
  population is the complete-case set — and Stage 7 restates the *pointer*, not the number: the
  excluded record is named in `covariate_completeness`, one section up in the same log.
- **Each row carries its own `n`**, because a `BALANCE_ONLY` covariate is not required to be
  complete: nothing masked on it, since it never entered the fit. [§11] wants the denominator per
  estimate, and here that is per *row*. On v7 every row reads 92 (§18) — `penumbra_ml` happens to be
  absent on exactly the record `core_ml` and `tmax6_ml` are absent on — so the column is a constant
  on this workbook and is rendered anyway, for the workbook that follows it.

**The per-row column is a constant on every frame this stage is tested against, and §12.2a is
therefore a constructed test rather than an observation.** Every workbook row reads 92 and every
fixture row reads 5, so nothing available exercises the mask that gives the column its meaning — and
a load-bearing column no frame varies is a column whose implementation could be global masking and
still pass. Measured (§18): blanking `penumbra_ml` on one `in_model` record of the fixture cohort
gives that row `n = 4` against every other row's 5, with a weighted SMD of 1.037070 computed over
those 4 records and a pooled SD over the same 4 — a different denominator from every other row in the
same table, which is what [§11] per-estimate means and what §12.2a asserts. `in_model` itself does not
move, because a `BALANCE_ONLY` covariate never entered the fit; that non-movement is the second half
of the test and the contrast that makes the first half mean something.

### 4.5 The preconditions

Five checks in **two phases**, each phase collected and raised together, following Stage 2's
`_assert_schema`, Stage 5's `_assert_cohort_inputs` and Stage 6's `_assert_design_inputs` and
`_assert_fit_inputs` — and departing from all four in exactly one way, which §4.5a is about.

```python
def _assert_balance_inputs(df: pd.DataFrame, ps: propensity.Propensity) -> None:
    # PHASE 1 — can this frame and this mask be READ? B1, B2, B3.
    # B4, B5 and everything in §8 read through all three, and pandas does not survive a mask that
    # is misaligned or three-valued, nor a column that is absent: it raises before any message here
    # is assembled (§4.5a, measured). So these three raise on their own rather than joining a
    # collection that never gets returned.
    bad: list[str] = []

    misaligned = [name for name, s in (("e", ps.e), ("w", ps.w), ("in_model", ps.in_model))
                  if not s.index.equals(df.index)]
    if misaligned:
        bad.append(
            f"B1  the Propensity is not aligned to this frame: {', '.join(misaligned)} "
            f"carr{'ies' if len(misaligned) == 1 else 'y'} a different index. "
            f"{len(ps.in_model)} mask row(s) against {len(df)} frame row(s), sharing "
            f"{len(ps.in_model.index.intersection(df.index))} index label(s). e, w and in_model are "
            "Series on the COHORT's index [Stage 6 §0.2], and a positional read of a misaligned pair "
            "produces a balance table for a population that does not exist, with every cell finite.")

    if ps.in_model.dtype != bool or ps.in_model.isna().any():
        bad.append(
            f"B2  in_model is {ps.in_model.dtype} and carries "
            f"{int(ps.in_model.isna().sum())} missing value(s). It is boolean and TOTAL by "
            "construction [Stage 6 §4.4]; a three-valued mask resolves by branch order at the "
            "first `if`, and here that decides who is in the denominator.")

    absent = [c for c in (*C.BALANCE_SET, C.TREATMENT) if c not in df.columns]
    if absent:
        bad.append(
            f"B3  {', '.join(absent)}: not a column of the frame. [§9] judges balance against the "
            "full [§6] confounder set, so a name this frame cannot supply is a balance table that "
            f"silently judges less than it claims — and {C.TREATMENT} is checked with them because "
            "B5 and every row of §5 read the arm from it.")

    if bad:
        raise C.SchemaError(
            "\n  " + "\n  ".join(bad)
            + f"\n  {len(bad)} balance assertion(s) failed over {len(df)} records. The frame or the "
            "mask cannot be read, so B4 and B5 were not run: they read through both (§4.5a).")

    # PHASE 2 — the frame and the mask are readable. These two describe the DATA, and they are
    # collected as Stages 2, 5 and 6 collect theirs.
    bad = []

    denied = [c for c in C.BALANCE_SET if c in C.POST_TIME_ZERO]
    if denied:
        bad.append(
            f"B4  {', '.join(denied)}: post-time-zero [invariant 4, §12]. Balance is a property of "
            "confounders measured at or before time zero; a post-exposure variable judged against "
            "[§9]'s threshold reads as a confounder that weighting failed to fix, when it is a "
            "variable no weighting should touch. onset_to_groin_min is reported by arm [§12].")

    present = [code for code in C.TREATMENT_LABELS
               if int((df.loc[ps.in_model, C.TREATMENT] == code).sum())]
    if len(present) < 2:
        bad.append(
            "B5  " + ", ".join(
                f"{C.TREATMENT_LABELS[c]}: "
                f"{int((df.loc[ps.in_model, C.TREATMENT] == c).sum())}"
                for c in C.TREATMENT_LABELS)
            + ". A standardised mean difference needs both arms, so with one every row of this "
            "table is undefined and the table is not a diagnostic. Stage 6's F2 established both "
            "arms were present before the fit; an arm missing HERE was lost by in_model.")

    if bad:
        raise C.SchemaError(
            "\n  " + "\n  ".join(bad)
            + f"\n  {len(bad)} balance assertion(s) failed over {len(df)} records, "
            f"{int(ps.in_model.sum())} of them weighted.")
```

- **B4 is invariant 4's cousin and not invariant 4 itself.** Stage 6 §4.5's D2 is the assertion the
  invariant asks for, on the one function every model's covariates pass through. Stage 7 builds no
  model, so nothing here can violate the invariant — but a post-time-zero variable *in a balance
  table* is a different error with the same cause, and this is the only place it can be caught.
  Measured (§18): `BALANCE_SET` and `POST_TIME_ZERO` are disjoint today, so B4 is unreachable on the
  declared set. Recorded rather than treated as a reason to omit it, exactly as Stage 5's C1-C4 and
  Stage 6's D1-D4 are.
- **B1 checks the index and not the length, on all THREE Series**, and **its message names which of
  them is misaligned and reports the shared label count.** Two Series of equal length on different
  indexes align to nothing under `.loc` and to the wrong rows under `.to_numpy()`, and the second is
  what this module does. Measured: on the fixture's shifted-index case the two lengths are both 5, so a
  message quoting only the lengths reads as though nothing were wrong — in the one case B1 exists for.
  The intersection is 0 there, and that is the number that says so. **`e` is checked because §8's
  `_centre` reads `ps.e.loc[weighted]`**: an earlier draft checked `in_model` and `w` only while its
  message claimed all three, and a `Propensity` whose `e` alone was misaligned reached `_centre` and
  raised pandas' `IndexingError` — **after `balance_smd` was already recorded**, leaving the log
  describing half a diagnostic that DoD-9 names for the `ess` case. Measured.
- **B3 covers `TREATMENT` as well as `BALANCE_SET`, and it is a phase-1 check.** `center` is in
  `BALANCE_SET` because [§6] makes it a confounder, so §6.1's per-centre loop is already protected;
  `ivt` is in neither list and is read by B5, by §8's `a`, and by every arm mask in §5 and §6. Naming
  it in B3 is necessary and was not sufficient: B5 reads the column two checks later, so B3 has to
  raise before B5 runs (§4.5a).
- **B5 is F2 one stage on**, and its message says which stage owns which guarantee, because the two
  are different: F2 is about the cohort, B5 is about the cohort *after* complete-casing.

### 4.5a Why the checks are in two phases, and why an earlier draft's single collection could not raise

Every assertion helper in Stages 2-6 collects its failures and raises once, and that pattern is right
wherever the checks are independent. **Here they are not: B4, B5 and everything in §8 read through the
frame and the mask that B1, B2 and B3 are about**, and pandas does not return a wrong answer for a
broken mask or an absent column — it raises, from inside the collection, before the `SchemaError` is
ever assembled.

Measured on pandas 2.3.3, against the frame and `Propensity` this stage receives, running earlier forms
of the function. All three land on B5's `df.loc[ps.in_model, C.TREATMENT]`:

```
  Propensity reindexed to a shifted index   →  B1 appended, then
      df.loc[ps.in_model, C.TREATMENT]      →  AssertionError('')      BARE, no message at all.

  in_model cast to object with one <NA>      →  B2 appended, then
      df.loc[ps.in_model, C.TREATMENT]      →  ValueError("Cannot mask with non-boolean array
                                                containing NA / NaN values")

  the TREATMENT column dropped               →  B3 appended, then
      df.loc[ps.in_model, C.TREATMENT]      →  KeyError('ivt')
```

Each message — B1's paragraph about a finite table for a population that does not exist, B2's about a
three-valued mask deciding the denominator, B3's naming the column — was assembled into `bad` and then
thrown away by an exception raised a check or two later. **The first case is the worst: an
`AssertionError` whose message is the empty string.** A caller sees no name, no explanation, and nothing
identifying which of the five assertions found anything, on the one failure mode B1 exists to describe.

**The third case is why the phase boundary sits where it does, and it arrived late.** A first repair put
B1 and B2 in phase 1 and left B3 with B4 and B5 — and added `TREATMENT` to B3's list, which named the
column without stopping B5 from reading it. So the same defect survived the fix, one check further
along, and was found only by running `assess` on a frame with the column dropped. The boundary is
therefore not "mask checks" against "the rest": it is **can this be read** against **is the data
judgeable**, and a column-presence check belongs to the first for exactly the reason a mask check does.

So the phases are not stylistic. §12.11 asserts both halves: that B1, B2 and B3 raise a `SchemaError`
naming themselves, and — the companions — that with the phase-1 raise removed each input produces the
bare pandas exception above. The companions are what stop the phases being tidied back into one
collection by a reader who sees five appends and one raise everywhere else in this repository.

**Two consequences for the messages.** Phase 1's raise says that B4 and B5 were not run, because a
`SchemaError` naming two of five assertions on a frame that fails three would otherwise read as a
clean bill of health for the other three. And Phase 2's raise keeps the `in_model` count in its
summary line, which Phase 1's cannot report honestly — `int(ps.in_model.sum())` on an object-dtype
mask carrying `<NA>` is not the number of weighted records.

## 5. The standardised mean difference [§9]

### 5.1 The rule

```
                       mean_w(x | treated)  −  mean_w(x | control)
    SMD(x, a, w)  =  ───────────────────────────────────────────────
                        √( var(x | treated) + var(x | control) ) / 2 )

    The NUMERATOR is weighted; the DENOMINATOR is not, and it is the same denominator in
    both columns [§9]. Sample variances, ddof = 1.
```

The unweighted column is this expression with `w ≡ 1`; the weighted column is this expression with
the [§7] overlap weights. **One function, called twice, with only the weights differing** — because
the alternative is two functions that can disagree about the denominator, which is precisely what
[§9]'s "so the yardstick does not move" forbids.

### 5.2 The estimator, written out

Written out in full, as Stage 5 §5.1 and Stage 6 §5.2 are: it is short, and every reported balance
number passes through it.

Three functions, and the split is §5.5's: `_pooled_sd` is the yardstick, `_ratio` applies a **given**
yardstick to a numerator, and `smd` is [§9]'s definition — the two over the same records.

```python
def _pooled_sd(x: np.ndarray, a: np.ndarray) -> float:
    """[§9]'s UNWEIGHTED pooled SD over the records given. nan where there is none. §5.1.

    The one definition of the yardstick. Sample variances, ddof = 1 — which §16b's oracle asserts
    rather than assumes, because a ddof = 0 convention differs by sqrt(n/(n-1)) and at n = 92 that
    is 0.5%, visible at six significant figures and invisible to the eye.

    An arm of fewer than two records has no sample variance, so there is no pooled SD: that is
    §5.3's branch 1 and it is a fact about the DENOMINATOR's population, which is why it lives here
    rather than in `_ratio`.
    """
    t, c = a == 1.0, a == 0.0
    if int(t.sum()) < 2 or int(c.sum()) < 2:
        return np.nan
    sd = float(np.sqrt((x[t].var(ddof=1) + x[c].var(ddof=1)) / 2.0))
    return sd if np.isfinite(sd) else np.nan
```

```python
def _ratio(x: np.ndarray, a: np.ndarray, w: np.ndarray, sd: float) -> float:
    """Weighted mean difference over a GIVEN yardstick. §5.1, and §6.1's foreign-SD case.

    Returns nan — never 0.0 — wherever the quantity is undefined, and §5.3 is the taxonomy.
    `pilots/analysis.py:194-208` returns 0.0 in one of those cases, and that case is a covariate
    which PERFECTLY SEPARATES the arms reported as perfectly balanced [§5.3a].

    **The two-record floor is a floor and no longer a variance requirement.** With `sd` supplied
    from another population, a weighted mean over one record is arithmetically fine — and a [§9]
    diagnostic standardising a single observation is not a diagnostic. So the guard stays, for a
    reason that changed: §5.3 branch 1 was about `var(ddof=1)`, this is about what a number means.
    Measured (§18): it is what keeps every within-centre row of the fixture cohort — one control per
    centre — reported as `missing` rather than as a number nobody should read.
    """
    t, c = a == 1.0, a == 0.0

    if int(t.sum()) < 2 or int(c.sum()) < 2:
        return np.nan                                  # §5.3, first — the FLOOR
    if not w[t].sum() or not w[c].sum():
        return np.nan                                  # §5.3, second — BEFORE np.average [§3.1]
    if not np.isfinite(sd):
        return np.nan                                  # §5.3, third
    if sd == 0.0:
        # §5.3a — the ARM CONSTANTS, never the weighted means. Both arms are constant whenever
        # sd is 0, so x[t][0] and x[c][0] ARE the arm means, exactly, with no summation in them.
        return 0.0 if x[t][0] == x[c][0] else np.nan   # §5.3, fourth

    return float((np.average(x[t], weights=w[t]) - np.average(x[c], weights=w[c])) / sd)
```

```python
def smd(x: pd.Series | np.ndarray, a: pd.Series | np.ndarray,
        w: pd.Series | np.ndarray) -> float:
    """The [§9] standardised mean difference: weighted means over an UNWEIGHTED pooled SD.

    Numerator and denominator over the SAME records — [§9]'s definition, and the only form Stage 12
    needs [§9, §14a]. Its signature is unchanged by §5.5's split and must stay unchanged: a `sd=`
    keyword here would make [§9]'s prescribed denominator look like a caller's option, which is
    Stage 6 §6.4's argument and §15's.

    Positional, not by keyword, and aligned by position rather than by index: the caller has
    already masked all three to the same records, and B1 is what makes that safe.
    """
    x = np.asarray(x, dtype=float)
    a = np.asarray(a, dtype=float)
    w = np.asarray(w, dtype=float)
    return _ratio(x, a, w, _pooled_sd(x, a))
```

Hand-verified (§18): `x = [1, 2, 3, 4]`, `a = [0, 0, 1, 1]`, `w ≡ 1` gives `sd = √((0.5 + 0.5)/2) =
0.70711` and `(3.5 − 1.5)/0.70711 = 2.8284271247461900`, against `2√2 = 2.8284271247461903` — the
roadmap's hand-computed acceptance criterion, agreeing to 4.4e-16.

### 5.3 The four ways it is undefined, and why each returns missing

The roadmap asks for one of these — "returns missing rather than zero for a variable observed in
only one arm" — and writing it revealed four, of which that is the third.

**Each branch belongs to one of the two populations §5.5 separates**, and the table says which: branch
1 fires on the **numerator's** records as a floor and, inside `_pooled_sd`, on the **denominator's** as
a variance requirement; branches 2 and 4 are the numerator's; branch 3 is the denominator's arriving as
a `nan`.

| # | Condition | Why no SMD exists | Returns |
|---|---|---|---|
| 1 | an arm has fewer than 2 records | in `_pooled_sd`, `var(ddof=1)` of one observation is `nan` [§3.1] so there is no yardstick; in `_ratio`, a floor — a standardised single observation is not a diagnostic (§5.2) | `nan` |
| 2 | an arm's weights sum to zero | the weighted mean is `0/0`; `np.average` **raises** rather than returning `nan` [§3.1], so the check must precede it | `nan` |
| 3 | a variable observed in only one arm | every value in the other arm is absent, the caller's mask leaves that arm empty, and branch 1 fires | `nan` |
| 4 | the pooled SD is exactly zero **and the arms differ** | the covariate is constant within each arm at two different values: it separates the arms perfectly, and the standardised difference is infinite rather than zero | `nan` |

**Two further routes into `nan` exist and are not branches, and they are recorded because `smd` is
public and Stage 12 calls it directly (§9).** `smd` does not mask; its docstring says the caller has
already masked all three vectors, and `_table` does exactly that per row (§8). A caller that does not:

| Input | What happens | Returns |
|---|---|---|
| a `nan` anywhere in `x` | that arm's `var(ddof=1)` is `nan`, so the pooled SD is not finite and **branch 3's own test fires** | `nan` |
| a `nan` anywhere in `w` | `w[t].sum()` is `nan`, `not nan` is `False`, so branch 2 does not fire and `np.average` propagates | `nan` |

Both measured (§18). Neither is a new branch and neither needs one: the first is already branch 3's
finiteness test doing what it was written for, and the second lands on `nan` by propagation rather than
by decision. What matters is that **neither returns a number** — a caller who forgets the mask gets the
undefined answer and not a plausible one, which is the property [§9] needs and the reason the guard
order in §5.2 is what it is. Stage 12 owes its own mask [§11]; this table is what it may assume if it
gets one wrong.

Branch 4's complement — the pooled SD is zero and the arms hold the **same** constant — is the one
case that legitimately returns `0.0`: the covariate is constant over the whole population, the two
arms are identical on it, and there is no difference to standardise. On v7 that is `center = USZ`,
which reads 0 before and 0 after (§18), and on the fixture cohort it is 14 of the 19 rows.

**`pilots/analysis.py:194-208` returns `0.0` for branch 4**, because its only guard is `sd > 0` and
its fall-through is `0.0`. A covariate on which every treated patient scores 1 and every control 0 —
complete separation on that one covariate, the strongest imbalance expressible — is reported as
`0.0`, passes `abs(smd) < SMD_THRESHOLD`, and appears in the balance table as the best-balanced row
in it. Measured (§18): `smd([0, 0, 1, 1], a=[0, 0, 1, 1], w=1)` is `nan` here and `0.0` there. This
is the second reading of the roadmap's own acceptance criterion and §19 lands it there.

### 5.3a The zero-variance branch compares the arm constants, and an earlier form did not

The natural way to write branch 4 is to compare the two **weighted means** for equality: if they
agree the covariate is constant everywhere, and if they do not it separates the arms. That form was
written, run, and is wrong.

Measured (§18), on the committed fixture cohort: `nihss_baseline` is **14.0 on every record in both
arms** — identical, no imbalance, the answer is 0.0 — and the two weighted means come back as `14.0`
and `14.000000000000002`. `np.average` divides a weighted sum by a weight sum, and the rounding of
that division depends on the weights, which differ by arm. The difference is 1.776e-15, the equality
fails, and a covariate with no imbalance at all is reported as **undefined**, which under §7's
rendering reads `missing` in the log and lands the row in `worst()`'s undefined list.

The fix is not a tolerance. It is that **whenever the pooled SD is exactly zero, both arm variances
are zero** — they are non-negative and they sum to zero — **so each arm is constant, and its constant
is `x[t][0]`**: the raw datum, with no summation, no weights and no rounding in it. Comparing those
two is exact by construction.

A tolerance would have been the wrong repair twice over: it would need a scale (the covariate's, not
1.0), and it would make branch 4 fire on covariates that are *nearly* constant, which is branch 4's
opposite. Measured after the correction (§18): every row of both frames is defined, the workbook's
19 rows are unchanged to every digit, and the constructed separation case still returns `nan`.

**Indicators are immune, which is why this hid.** `np.average` of a vector of exact `1.0`s returns
exactly `1.0` and of `0.0`s exactly `0.0` (§18), so 17 of the fixture's 19 rows could never have
shown it; the one integer-valued constant covariate did. A test suite built only from factor levels
would have passed.

### 5.4 The threshold is a verdict, and never a raise

`SMD_THRESHOLD` is 0.10 [§9] and Stage 1 declares it. Stage 7 **reports against it and does not
raise on it.** Residual imbalance is a finding for [§16] and, above the threshold, a question for the
PI under [§13]; it is not a malformed input and there is nothing for a caller to catch. `assess`
raises only on B1-B5, which are all statements about the *frame* rather than about the *data*.

This is the one place a reader might expect Stage 6's posture and not get it, so it is stated
plainly: Stage 6 raises because a fit that has gone wrong returns a number nobody can tell from a
right one, and Stage 7 does not because a balance table that has gone *badly* is a balance table
doing its job. §6.4 is what that looks like on this cohort.

### 5.5 One yardstick for the whole stage, including within a centre

[§9] says "common unweighted pooled SD, so the yardstick does not move". §4.4 applied that to the two
*columns* and to the *rows*. **It applies to §6.1's within-centre `worst |SMD|` column too, and an
earlier draft of this document did not apply it there** — `_centre` called the same `smd` over the
centre's records, so that column re-derived its denominator from each centre's own variance.

What that cost, measured (§18):

```
    centre    worst within-centre |SMD|      local SD      pooled in_model SD
    ─────────────────────────────────────────────────────────────────────────
    HUG       diabetes                        +0.5385            +0.4956
    CHUV      onset_type = unwitnessed        +2.0608            +1.4100
    Lugano    hyperlipidemia                  −1.5275            −1.1087
```

Up to 46% apart. Under the local SD no two rows of that column shared a yardstick, none shared one with
`SMD after` in the entry two sections above it, and CHUV's 2.06 against Lugano's 1.53 was **partly an
artefact of CHUV's covariate having the smaller within-centre variance** rather than the larger
imbalance. A table whose rows are not comparable with each other, in a document that argues six times
that one SD must serve two columns, is the same defect at ninety degrees.

**So the denominator is the pooled unweighted SD over `in_model`, everywhere in this stage.** The
numerator is the centre's; the yardstick is the whole weighted set's. `CovariateBalance.sd` carries it
(§3), §8 harvests it from the pooled table and hands it to every centre row, and §12.6 asserts that a
centre row's `sd` equals the pooled row's for the same covariate.

Two things this does **not** change, both measured (§18): the covariate identified as worst at each
centre is the same under either denominator — `diabetes`, `onset_type = unwitnessed`, `hyperlipidemia` —
and so is the order of the three centres. What changes is that the numbers can now be read against each
other and against `balance_smd`. **And §13's warning stands unaltered**: the arm sizes are still in the
same row, an |SMD| over an arm of two still carries almost no information, and a reader must still not
rank the centres by this column. One yardstick removes one of the two reasons not to; arm size is the
other and it is the bigger one.

## 6. Overlap [§9]

### 6.1 The within-centre table, and the pooled row inside it

[§9] asks for overlap "within each centre as well as pooled", and roadmap Stage 7 lists the columns:
numbers per arm, propensity range per arm, ESS, maximum weight, weight share, worst within-centre
|SMD|. **The pooled report is the `all (pooled)` row of that same table, and not a second table** — one
grid in the log, built by one function, from one code path.

**In the return value it is a field of its own, `Balance.pooled`, and not a fifth element of
`centres`** (§3). The two statements are compatible and the distinction is deliberate: the log's reader
wants one table with a labelled row, and a Python caller wants `centres` to mean the declared centres
and nothing else. §8 renders `centres + (pooled,)`; §3 records what the earlier one-tuple form cost.

Two reasons, and the second is the load-bearing one:

- Stage 6's `overlap_weights` entry already carries the pooled arm table — `n`, `Σw`, ESS, max `w`,
  weight share, per arm plus pooled (Stage 6 §7.2). A third rendering of those cells would be the
  duplication Stage 6 §7.4 declined for `absence_by_column`, in the same log, two sections apart.
- A pooled row **computed by the same code as the centre rows** cannot disagree with them about what
  a column means. §12.8 then asserts that its ESS and max weight equal Stage 6's to the digit, which
  makes the two entries reconcile without either recomputing the other.

What Stage 7's row adds over Stage 6's is the **propensity range per arm**, which is the quantity
overlap is actually about and which Stage 6 does not carry, and the **worst within-centre |SMD|**.

**`in_cohort` and `weighted` are two columns and not one.** A centre's cohort count and its count
inside `in_model` differ wherever a covariate-incomplete record sits, and on v7 they do: Lugano is
**31 in the cohort and 30 weighted** (§18). One column would have to choose, and either choice makes
one of the two entries in the log wrong.

### 6.2 Structural non-positivity is a status, not an absence

[§9]: "Centres with structural non-positivity are reported as such, not given an overlap plot."

The table therefore has **one row per declared centre**, ranging over `CENTER_ORDER` and never over
the cohort's observed centres, and a centre with no cohort record gets a row whose `status` says so
and whose overlap statistics are all `missing`. On v7 that is USZ, which [§3] restriction 1 removed
because it contributed no bridging patient (§18).

Three statuses, and each names a different fact:

| status | when | what it tells Stage 14 |
|---|---|---|
| `overlap reported` | the centre has both arms inside `in_model` | plot it |
| `structural non-positivity [§3]` | the centre has no record in the cohort | **do not plot it**; [§9] forbids it, and `restrict_centres` names which restriction removed it |
| `no weighted contrast` | the centre is in the cohort and `in_model` left one of its arms **empty** | do not plot it; this is not [§3]'s structural case but complete-casing, and the two must not read alike. An arm of **one** is not this case (§6.3): it keeps `overlap reported` and reports a `missing` worst |

**Three, and there is no fourth.** The `all (pooled)` row goes through the same `_status` call and
reads `overlap reported`, because it has both arms by construction — B5 raises otherwise. An earlier
draft gave it a status literal of its own, `"pooled"`, which put a fourth value in the column Stage 14
matches on while §6.2 still described three; the pooled row is now distinguished by being
`Balance.pooled` and by its `centre` cell reading `all (pooled)`, neither of which is a status (§3).

**Stage 7 does not distinguish "removed by restriction 1" from "contributed no patient at all", and
does not need to.** Both are declared centres absent from the cohort, both are structurally
non-positive for this analysis, and Stage 5's `restrict_centres` entry — three sections up the same
log — already names which patients each removed. Inventing a second answer here would need the
unrestricted frame as an argument, and a diagnostic that takes the pre-restriction frame in order to
label a row is a diagnostic that can be given the wrong one.

**Stage 7 draws nothing.** Every figure is Stage 14's [§16]; what this stage owes is the column that
tells Stage 14 which centre may be drawn, and §9 states the obligation on Stage 14's side too.

### 6.3 `ess` is called, never reimplemented — and never allowed to raise

`propensity.ess` is public precisely so this stage does not write a second Kish sum (Stage 6 §3), and
`assess` calls it per centre per arm.

But **Stage 7 checks the arm before it calls**, and that is a decision rather than defensive
padding. `ess` raises `model.FitError` on an empty arm, with a message that reads "the in_model mask
lost one" (Stage 6 §6.2) — correct for the [§7] fit, where an empty arm means something went wrong,
and wrong here, where an empty arm within one centre is a **finding to report**. Stage 5's P2
guarantees both arms at every retained centre; nothing guarantees both arms at every retained centre
*after complete-casing*, and on a workbook where one centre's only control is missing a CTP volume
the diagnostic would abort mid-table rather than saying so.

So: `assess` renders `missing` for that arm's ESS and range, sets the row's status to `no weighted
contrast`, and continues. Stage 6's raise stays the backstop for a caller that passes an empty arm
anyway, and §12.7 asserts both halves — that Stage 7 does not raise, and that `ess` still would.

**An arm of one is not an empty arm, and the two land in different columns.** Measured (§18): on the
fixture cohort every centre has exactly **one** control, so `ess` is perfectly computable there —
Kish of a single weight is 1.0 — and so is the propensity range, while **every within-centre SMD is
undefined**, because §5.3's first branch needs two records per arm to have a variance. Those rows
therefore read `overlap reported` with a `missing` worst |SMD|: the centre's overlap *is* reported
and its balance is not judgeable, which are two different facts in two different columns.

So `no weighted contrast` — the status — needs an arm of **zero**, and **no frame available to this
stage produces one**: Stage 5's P2 guarantees both arms at every retained centre, and on both the
workbook and the fixture `in_model` empties none of them. It is reached only by a constructed frame,
which §12.7 constructs by blanking a covariate on a centre's only control. Recorded here rather than
left to be discovered, in the posture Stage 6 §12.4a uses for a branch no natural input reaches: the
guard is what makes the status possible, the status is what makes it reportable, and neither has a
witness on this workbook.

### 6.4 What the residual imbalance is, measured — and whose it is

This is the section Stage 6 §13 asked for by name.

**On the workbook, after weighting, 5 of the 19 rows reach [§9]'s threshold** (§18):

```
    row                       role                before    after
    ──────────────────────────────────────────────────────────────
    center = HUG              propensity model    +1.343    +0.211
    center = Lugano           propensity model    −1.224    −0.154
    hypertension              negative control    +0.051    −0.139
    hyperlipidemia            negative control    −0.390    −0.301
    diabetes                  negative control    +0.219    +0.380
```

Eleven of the nineteen were above the threshold before weighting and five are after. **Two of the
five are in the propensity model**, and that is the finding: overlap weights have an exact-balance
property, and [§7]'s estimator does not deliver it.

**The MLE companion measures exactly whose the residual is.** With the propensity score replaced by
an unpenalised `statsmodels.Logit` fit on the same design (§18):

```
    every in-model row            |SMD| ≤ 1.03e-15        exactly balanced
    Σw treated vs control         12.9716 vs 12.9716      equal to 1.8e-15
    the four negative controls    worst 0.391             NOT balanced, and never were
```

Two things follow, and they are different things:

- **The in-model residual is Firth's, in full.** Overlap weights solve the balance equation for every
  covariate in the score *because that is the score equation*; Firth's modified score is
  `X'(y − p + h(0.5 − p)) = 0`, a different equation, so the property holds only approximately.
  **And this is why the MLE balances `center = HUG`, which is not a column of the design at all.** The
  exact-balance property applies to the intercept and to every design column, and `center_HUG` is
  `1 − center_CHUV − center_Lugano − center_USZ` with `center_USZ ≡ 0` on this cohort (§4.2) — so its
  balance is *implied* by the columns that are in the score. The 1.03e-15 on a row the design does not
  contain is therefore evidence about the estimator and not a coincidence, which is what makes it worth
  reporting; the same algebra is why `onset_type = witnessed` balances too.
  Stage 6 §6.3 established this on `Σw` — a difference of 0.484, which §18d there traced to
  `Σ h(e − 0.5)` — and this section is that same discrepancy expressed per covariate.
  **Stage 6 §13 predicted `0 < worst < 0.05` from the pilots' synthetic measurement. On this cohort
  it is 0.211.** The prediction was about synthetic data and is not withdrawn as a measurement; it is
  simply not what this cohort does, and the correction belongs here because this is the stage that
  can measure it.
- **The negative controls' residual is nobody's, and that is what they are for.** They are imbalanced
  under Firth (worst 0.380) and under the MLE (worst 0.391) alike, because exact balance is a
  property of the covariates *in* the score and of no others. [§6] retains them "where residual
  imbalance shows what weighting does not fix", and 0.38 on `diabetes` is that sentence with a number
  in it. `diabetes` is also **worse after weighting than before** (0.219 → 0.380), which is not a
  malfunction: the weights tilt the population toward equipoise, and a covariate outside the score is
  free to become more imbalanced in the tilted population than it was in the original one.

**Why `center` and not something else.** The two rows left above the threshold are the two that
started furthest from it: `center = HUG` at 1.343 and `center = Lugano` at 1.224, against at most
0.428 for everything else. Of the eight in-model rows whose unweighted |SMD| exceeds 0.10, weighting
removes between **66% and 94%** of it (§18) — the centre rows are not treated worse than the others
in proportion; they simply had further to go. `center` is where the treatment assignment lives in
this dataset, which is why [§3] restricts on it and why [§6] adjusts for it.

**What Stage 7 does about it: reports it.** [§7] prescribes one estimator, refitted identically in
every replicate, and forbids substituting another. The MLE companion above is an *oracle*, run in a
test, and swapping it into the pipeline would exchange a documented residual imbalance for the
mixture-of-two-estimators bootstrap [§7] was written against — Stage 6 §12.10's pair exists to catch
exactly that substitution. §13 files what is owed: the PI reads this under [§9], and [§13]'s deferred
sensitivity specifications are where a different propensity model would be a *declared* one.

## 7. The two `model` audit entries

### 7.1 No new kind, and that is Stage 6's declaration honoured

`data.py:144-149` says of `model`: "the last one this pipeline needs — Stages 7-13 all render under
`model` or under an existing kind", and Stage 6 §10 says the same in the ledger that landed it.
Stage 7 is the first stage to test that claim and it holds.

A balance table is not obviously a "fitted model", so the argument has to be made rather than
assumed. It is this: `model`'s entries are *claims about a fit* — the population it ran on, the
matrix it used, the coefficients it produced, the weights it implies — and a standardised mean
difference is a claim about the same fit, namely how well its weights did what they were for. It
renders under `## Fitted models`, immediately after `overlap_weights`, and the log then reads: this
is what was fitted, these are the weights, this is what the weights achieved. A separate heading
would put the verdict in a different section from the thing it judges.

So **`data.py` is not amended by this stage** — no `KINDS` member, no `_HEADINGS` entry, and none of
the four literal pins in `test_data.py`, `test_eligibility.py`, `test_cohort.py` and
`test_propensity.py` moves. §12.9 asserts the count is still nine.

### 7.2 The two entries

| step | kind | `n` | contents |
|---|---|---|---|
| `balance_smd` | model | records **judged** (92 on v7) | one row per declared level of every factor plus every linear covariate in `BALANCE_SET`: `covariate`, `role`, `n`, `SMD before`, `SMD after`, `\|SMD\| < 0.10`. 19 rows on v7 — **not** the `n` |
| `overlap_by_centre` | model | centres with a weighted **both-arm** contrast (3 on v7) | one row per declared centre plus `all (pooled)`: `n`, `n` per arm, `e` range per arm, ESS per arm, `max w`, `weight share`, `worst \|SMD\|`, `status` |

Four rules about the tables, all inherited:

- **Every declared thing is rendered whether or not the data fills it.** Every level of every factor,
  including the reference and including a level nobody has; every centre in `CENTER_ORDER`,
  including one the [§3] restriction removed. `pd.crosstab` is used nowhere in this repository for
  that reason.
- **Every cell is computed from the object it describes, never by subtracting another row** (Stage 5
  §7.3). The pooled row is computed over `in_model`, not by combining the centre rows.
- **`_fmt` is the only float formatter and `missing` the only rendering of an absent value**, so an
  undefined SMD renders as `missing` and never as `nan`, `0`, or a numpy repr. §12.9's byte-identity
  criterion rests on it. The `e` range is one cell built as `f"{_fmt(lo)}–{_fmt(hi)}"`, so both
  numbers still pass through `_fmt`.
- **The verdict column is a string, never a bool.** `_fmt(True)` is `float(True)` is `"1"`, so a
  boolean verdict would render as `1` beside an SMD and read as a number. The three values are
  `yes`, `no` and `undefined`, and the third is why the column cannot be two-valued: a row whose SMD
  is `missing` has not passed the threshold and has not failed it.

`n` deserves its own sentence in each case, because in both entries it is **not** the row count —
the same distinction Stage 6 §7.5 drew for `design_matrix`. `balance_smd`'s `n` is the records the
table is computed over, so it can be read against `propensity_fit`'s `n` immediately above it;
`overlap_by_centre`'s `n` is the number of centres that actually contributed a contrast, so a centre
reported as structurally non-positive is visible as a row *and* absent from the count.

### 7.3 What the tables deliberately do not carry

**No arm means.** Stage 6 §7.6's `overlap_weights` second table already carries the unweighted and
ATO-weighted mean per arm for every [§6] covariate, expanded per declared level — 14 of these 19
rows — and rendering them again two sections down would be the fourth-copy problem §7.4 there
declined. The five rows it does not cover are the negative controls and `penumbra_ml`, and their arm
means belong in [§16]'s baseline table, which is Stage 14's: what [§9] asks about a negative control
is its residual imbalance, not its prevalence. §15 records the decision.

**No variance ratios.** `pilots/analysis.py:249-265` reports one per covariate. [§9] specifies the
SMD and a threshold and says nothing about second moments, and a column with no threshold beside it
is a column every reader interprets privately. §15.

### 7.4 The two tables and the three statuses, written out

Written out rather than described, because Stage 6 §20's round 5 found that its own two prose-only
tables could not be built from the arguments the document gave their helpers (finding 51), and a
table specified in a column list is a table nobody has run.

```python
_REPORTED: Final[str] = "overlap reported"
_STRUCTURAL: Final[str] = "structural non-positivity [§3]"
_NO_CONTRAST: Final[str] = "no weighted contrast"


def _status(at: pd.Series, arms: dict[int, np.ndarray]) -> str:
    """§6.2's three, and the order is the specification. Called for the pooled row too.

    `at` is COHORT membership and `arms` is what survived `in_model`, so the two branches are
    different facts: a declared centre with no cohort record was removed by [§3] restriction 1,
    and a centre with records but an empty arm was emptied by complete-casing. Reading the second
    as the first would report a data-handling consequence as a design one.

    There is no pooled branch and no fourth constant. The pooled row's `at` is True everywhere and
    B5 has already established both arms, so it reaches `_REPORTED` — which is true of it. What
    marks it as pooled is `Balance.pooled` and its `centre` cell, never this column (§3, §6.2).
    """
    if not int(at.sum()):
        return _STRUCTURAL
    if any(not int(mask.sum()) for mask in arms.values()):
        return _NO_CONTRAST
    return _REPORTED
```

```python
def _verdict(value: float) -> str:
    """`yes`, `no` or `undefined` — a STRING, because `_fmt(True)` is `float(True)` is "1" [§7.2].

    Three-valued and not two: a row whose SMD is undefined has neither passed the threshold nor
    failed it, and rendering it `no` would put a covariate nobody could judge among the findings.
    """
    if not np.isfinite(value):
        return "undefined"
    return "yes" if abs(value) < C.SMD_THRESHOLD else "no"


def _smd_table(rows: Sequence[CovariateBalance]) -> tuple[tuple[str, ...], ...]:
    """One row per declared level, in BALANCE_SET order. §7.2.

    The threshold is interpolated into the header from `C.SMD_THRESHOLD` rather than written as
    `0.10`, so the column cannot claim a threshold the verdict was not computed against. It renders
    as `0.1`, because `_fmt` is six significant figures — the header says what the code compares.

    **No cell in this table contains a `|`.** `data._md_table` does no escaping (data.py:242-245):
    it joins cells with " | " and computes the separator's width from `len(rows[0])`, so one pipe
    inside a header cell gives a header row with more markdown cells than its own separator and
    body. Measured: the column was `|SMD| < 0.1` and rendered 8 cells against a 6-cell separator,
    in every log this stage writes. §12.9 asserts the absence of the character.
    """
    header = ("covariate", "role", "n", "SMD before", "SMD after",
              f"abs SMD < {_fmt(C.SMD_THRESHOLD)}")
    return (header, *(
        (row.covariate, row.role, str(row.n), _fmt(row.unweighted), _fmt(row.weighted),
         _verdict(row.weighted))
        for row in rows))
```

```python
def _range(bounds: tuple[float, float]) -> str:
    """The propensity range as ONE cell, both numbers through `_fmt`. §7.2.

    One cell rather than two columns per arm, which would take this table to sixteen. `missing`
    rather than `missing–missing` when the arm is empty, so an absent range reads as the single
    fact it is.

    **An arm of ONE renders as `0.257067–0.257067`, and that is deliberate.** Measured on the fixture
    cohort, where every centre has exactly one control. Collapsing it to a single number would make
    the cell's width carry information the `n` column already carries, and would make an arm of one
    indistinguishable from an arm of many whose propensities happen to coincide — which is a real
    state and a different one. A degenerate range that LOOKS degenerate is the honest rendering; the
    reader learns the arm size from the column two cells to the left [§9].
    """
    low, high = bounds
    return "missing" if not np.isfinite(low) else f"{_fmt(low)}–{_fmt(high)}"


def _overlap_table(rows: Sequence[CentreOverlap]) -> tuple[tuple[str, ...], ...]:
    """One row per declared centre, then `all (pooled)`. §6.1, §7.2.

    `worst abs SMD` renders the absolute value the column names, while `CentreOverlap.worst` keeps
    the sign for a caller that wants to know which arm it favours. It is standardised by the POOLED
    in_model SD, not by the centre's own — §5.5, and the column name says nothing about that because
    a column name cannot; `_overlap_detail` does.

    **No cell contains a `|`**, for `_smd_table`'s reason: the column was `worst |SMD|` and rendered
    a 15-cell header against a 13-cell separator and body. §12.9 asserts it.
    """
    labels = C.TREATMENT_LABELS
    header = ("centre", "in cohort", "weighted",
              *(f"n ({label})" for label in labels.values()),
              *(f"e ({label})" for label in labels.values()),
              *(f"ESS ({label})" for label in labels.values()),
              "max w", "weight share", "worst abs SMD", "status")
    return (header, *(
        (row.centre, str(row.in_cohort), str(row.weighted),
         *(str(row.n[code]) for code in labels),
         *(_range(row.e_range[code]) for code in labels),
         *(_fmt(row.ess[code]) for code in labels),
         _fmt(row.max_weight), _fmt(row.weight_share), _fmt(abs(row.worst)), row.status)
        for row in rows))
```

```python
def _smd_detail(df: pd.DataFrame, ps: propensity.Propensity,
                rows: Sequence[CovariateBalance]) -> str:
    worst, undefined = _worst(rows)              # §3 — never a throwaway Balance
    over = _over_threshold(rows)                 # §3 — the SAME predicate `unbalanced()` uses
    return (
        f"[§9] standardised mean differences over {len(rows)} declared level(s) of "
        f"{len(C.BALANCE_SET)} [§6] covariate(s), before and after the [§7] overlap weights, on "
        f"{int(ps.in_model.sum())} of {len(df)} cohort record(s) — the weighted set, in BOTH "
        "columns, so the pooled standard deviation is one yardstick over one population "
        "[Stage 7 §4.4]. The denominator is the UNWEIGHTED pooled SD in both, which is what makes "
        "the two columns comparable [§9]. Balance is judged against the FULL [§6] confounder set "
        f"regardless of what this specification fitted, and the {len(C.NEGATIVE_CONTROLS)} vascular "
        "risk factors are "
        "negative controls: residual imbalance on them shows what weighting does not fix [§6]. "
        f"{len(over)} row(s) reach |SMD| = {_fmt(C.SMD_THRESHOLD)} after weighting"
        + (": " + ", ".join(over) if over else "")
        + (f"; worst {_fmt(abs(worst.weighted))} on {worst.covariate}" if worst else "")
        # The trailing full stop is load-bearing and an earlier draft of this fence omitted it:
        # without it the clause runs into the next sentence — "…, penumbra_ml A row's `n` is its own
        # denominator" — which is visible only by reading a rendered `detail` as prose, on a frame
        # that HAS an undefined row. Neither review round had one (§18, §21.5).
        + (f". {len(undefined)} row(s) are undefined and are reported as such rather than as zero: "
           + ", ".join(undefined) + "." if undefined else ". No row is undefined.")
        + " A row's `n` is its own denominator [§11]; the record(s) outside the fit are named in "
        "`covariate_completeness` above. This entry reports balance and does not judge the "
        "estimator: [§7] prescribes one propensity model and there is no second one to try.")


def _overlap_detail(centres: Sequence[CentreOverlap]) -> str:
    absent = [c.centre for c in centres if c.status == _STRUCTURAL]
    thin = [c.centre for c in centres if c.status == _NO_CONTRAST]
    return (
        f"[§9] overlap within each centre as well as pooled, over {len(C.CENTER_ORDER)} declared "
        "centre(s). A pooled distribution can look acceptable while treatment is nearly determined "
        "by centre, and only the within-centre view distinguishes patient-level equipoise from "
        f"that. {len(absent)} centre(s) contribute no cohort record and are reported as "
        "structurally non-positive rather than given an overlap plot [§9]"
        + (": " + ", ".join(absent) if absent else "")
        + (f". {len(thin)} centre(s) lost an arm to complete-casing and carry no weighted contrast: "
           + ", ".join(thin) if thin else "")
        + ". `restrict_centres` above names which patients each [§3] restriction removed. The "
        "pooled row is computed by the same code as the centre rows and reconciles with "
        "`overlap_weights`; no second effective sample size is computed anywhere in this stage "
        "[Stage 6 §11]. Every `worst abs SMD` in this table is standardised by the pooled unweighted "
        "standard deviation over the WEIGHTED SET, not by the centre's own, so the column is "
        "comparable across these rows and against `balance_smd` above [Stage 7 §5.5]. It is still "
        "not a ranking: an abs SMD over an arm of two carries almost no information, which is why "
        "the arm sizes are in the same row [Stage 7 §13].")
```

## 8. `assess`, written out

Canonical, as Stage 5 §8 and Stage 6 §8 are.

```python
def assess(df: pd.DataFrame, ps: propensity.Propensity, audit: Audit) -> Balance:
    """The [§9] balance and overlap diagnostics over the [§3] cohort and its [§7] weights.

    Takes no covariate list and no threshold: the balance set is BALANCE_SET and the threshold is
    SMD_THRESHOLD, both prespecified [§9]. Adds no column to `df`, edits nothing, refits nothing
    and computes no second effective sample size [Stage 6 §11].

    Raises SchemaError on B1-B5 and on nothing else. Imbalance is a finding, not a failure [§5.4].
    """
    _assert_balance_inputs(df, ps)                             # B1, B2 then B3-B5 — §4.5

    sub = df.loc[ps.in_model]
    a = sub[C.TREATMENT].to_numpy(dtype=float)
    w = ps.w.loc[ps.in_model].to_numpy(dtype=float)

    covariates = _table(sub, a, w, C.BALANCE_SET)
    audit.record("model", "balance_smd", int(ps.in_model.sum()),
                 _smd_detail(df, ps, covariates), table=_smd_table(covariates))

    # §5.5 — the ONE yardstick, harvested from the table computed over in_model and handed to every
    # centre row. A centre's numerator is its own; its denominator is the whole weighted set's.
    sds = {row.covariate: row.sd for row in covariates}
    centres = tuple(_centre(df, ps, centre, sds) for centre in C.CENTER_ORDER)
    pooled = _centre(df, ps, None, sds)                        # the `all (pooled)` row — §6.1
    audit.record("model", "overlap_by_centre",
                 sum(1 for c in centres if c.status == _REPORTED),
                 _overlap_detail(centres), table=_overlap_table(centres + (pooled,)))

    return Balance(covariates=covariates, centres=centres, pooled=pooled)
```

```python
def _table(sub: pd.DataFrame, a: np.ndarray, w: np.ndarray, names: Sequence[str],
           sds: dict[str, float] | None = None) -> tuple[CovariateBalance, ...]:
    """One CovariateBalance per declared level, over the records `sub` already restricts to.

    The per-row mask is the covariate's own: [§11] is complete-case PER ESTIMATE, and a row of this
    table is an estimate. `n` is what survives it, which is why the column exists (§4.4), and
    §12.2a is what tells this form from one mask over the frame.

    `sds` is §5.5's yardstick. None means "compute each row's own", which is right for the table
    over `in_model` — that IS the population the yardstick is defined on. A dict means "use these",
    which is every within-centre call: the numerator is the centre's and the denominator is not.
    It is a positional-or-keyword parameter with a default rather than two functions, because the
    two calls differ in one argument and nothing else; and it is NOT on `smd`, whose signature
    stays [§9]'s (§5.2).
    """
    rows: list[CovariateBalance] = []
    for name in names:
        for label, indicator in _levels(sub, name):
            values = indicator.to_numpy(dtype=float)
            present = np.isfinite(values)
            x, arm, weight = values[present], a[present], w[present]
            sd = _pooled_sd(x, arm) if sds is None else sds.get(label, np.nan)
            rows.append(CovariateBalance(
                covariate=label, role=_role(name), n=int(present.sum()), sd=sd,
                unweighted=_ratio(x, arm, np.ones_like(weight), sd),
                weighted=_ratio(x, arm, weight, sd)))
    return tuple(rows)
```

**`sds.get(label, np.nan)` and not `sds[label]`**, for one reason worth stating: a label the pooled
table did not produce has no yardstick, and the honest answer for its row is `missing` rather than a
`KeyError` inside a diagnostic. It cannot happen today — every within-centre call passes a subset of
`BALANCE_SET`, so every label is present — and §12.6 asserts that, so the `get` is the branch that
stays inert while the assertion is what would notice.

```python
def _centre(df: pd.DataFrame, ps: propensity.Propensity, centre: str | None,
            sds: dict[str, float]) -> CentreOverlap:
    """One row of §6.1's table. `centre=None` builds the `all (pooled)` row from the same code.

    `ess` is called only where the arm is non-empty (§6.3): an arm complete-casing emptied is a
    finding this table reports, not an exception that stops it being written.

    `sds` is §5.5's one yardstick and is REQUIRED, with no default: a centre row computed against
    its own variances is the defect §5.5 records, and a defaulted `None` here is the one line that
    would reintroduce it silently.
    """
    at = pd.Series(True, index=df.index) if centre is None else df["center"] == centre
    weighted = at & ps.in_model
    sub = df.loc[weighted]
    e, w = ps.e.loc[weighted], ps.w.loc[weighted]
    arms = {code: (sub[C.TREATMENT] == code).to_numpy() for code in C.TREATMENT_LABELS}

    worst, defined = np.nan, [
        row.weighted for row in _table(
            sub, sub[C.TREATMENT].to_numpy(dtype=float), w.to_numpy(dtype=float),
            # `center` is dropped INSIDE a centre: its indicators are constant there by
            # construction, so those rows can only ever read 0.0 and a row that cannot vary is
            # not evidence. Measured: dropping them changes no centre's worst (§18).
            [c for c in C.BALANCE_SET if c != "center"] if centre is not None else C.BALANCE_SET,
            sds)                                    # §5.5 — the pooled yardstick, not the centre's
        if np.isfinite(row.weighted)]
    if defined:
        worst = float(max(defined, key=abs))

    return CentreOverlap(
        centre="all (pooled)" if centre is None else centre,
        in_cohort=int(at.sum()), weighted=int(weighted.sum()),
        n={code: int(mask.sum()) for code, mask in arms.items()},
        e_range={code: ((float(e[mask].min()), float(e[mask].max())) if mask.any()
                        else (np.nan, np.nan)) for code, mask in arms.items()},
        ess={code: (propensity.ess(w[mask]) if mask.any() else np.nan)
             for code, mask in arms.items()},
        max_weight=float(w.max()) if len(w) else np.nan,
        weight_share=float(w.sum() / ps.w[ps.in_model].sum()) if len(w) else np.nan,
        worst=worst,
        status=_status(at, arms))
```

Five things about the order and the shape, each of which is a failure if moved:

- **`_assert_balance_inputs` comes first**, before the mask and before any arithmetic — Stage 5 §4.4's
  rule. A misaligned `Propensity` otherwise produces a complete table for a population that never
  existed, and §4.5's phase 1 is what makes that a message rather than a bare `AssertionError`.
- **`sub`, `a` and `w` are bound once**, from `in_model`, and every row of every table is computed
  from those three. Two separate `.loc` reads would be two chances to align the weights to a
  different subset than the means.
- **`_table` is called by both entries, with the same signature and the same arguments**, so the
  pooled worst |SMD| and the balance table's worst are the same number computed by the same call.
  It is computed **twice** — once here, once inside `_centre(df, ps, None)` — and that is the price of
  the one-code-path property below, priced in §2 at 7.4 ms. Not "computed once": an earlier draft said
  so, and a reader who believed it would look for the caching that would have made it true. §12.8
  asserts the two agree, which is the property that matters and the one a shared computation would
  have made unfalsifiable.
- **`_centre(df, ps, None, sds)` builds the pooled row**, so the pooled and per-centre cells cannot
  mean different things (§6.1). Two differences are deliberate and both are at the point of use: the
  `center` indicators are judged pooled and dropped within a centre, and `sds` is threaded so that a
  centre's numerator is its own while its yardstick is the weighted set's (§5.5). It goes into
  `Balance.pooled` rather than into `centres`, so `centres` is exactly `CENTER_ORDER` (§3).
- **`sds` is harvested from `covariates`, after the first entry is recorded and before any centre
  row.** That ordering is the whole of §5.5's implementation: the yardstick exists as a value only
  because `CovariateBalance` carries `sd`, and it is read from the table computed over `in_model` —
  never recomputed, so the pooled row and the centre rows cannot disagree about it. §12.6 asserts a
  centre row's `sd` equals the pooled row's for the same covariate.
- **`max(defined, key=abs)` rather than a pandas `max`.** §3.1: `Series.abs().max()` skips missing
  values and `np.max` poisons on them. Filtering to the finite rows explicitly, and counting what was
  filtered, is what `Balance.worst()` exposes and what neither library expression can express.

## 9. Data flow into Stages 8-14

```
  assess(cohort, ps, audit) returns Balance
   ├─ covariates  19 CovariateBalance over 14 declared names        [§4.1, §4.2]
   ├─ centres     4 CentreOverlap — exactly CENTER_ORDER            [§6.1]
   ├─ pooled      1 CentreOverlap — the `all (pooled)` row          [§3, §6.1]
   ├─ worst()     (the largest defined |SMD|, the undefined names)  [§3]
   └─ unbalanced() the rows reaching SMD_THRESHOLD                  [§5.4]

  the cohort frame and the Propensity come back unchanged                    [§0.2]

  audit — the same object load() created:
    load 7 / derive 4 / classify 1 / build 6 / fit 4 / assess 2  =  24       [§18]
```

What each later stage may rely on, and what each owes:

- **Stages 8, 9 and 10** — read nothing from this stage. Balance does not enter an estimate, and no
  estimator is conditioned on it. In particular **Stage 10 does not recompute balance per replicate**
  (§2): the bootstrap resamples the population the point specification defines, and a per-replicate
  balance table would be a distribution of diagnostics with no null to test it against.
- **Stage 11 [§13]** — the deferred sensitivity specifications each report "ESS and worst residual
  |SMD|". The ESS comes from `Propensity.ess` and the worst from `Balance.worst()`, **which returns
  the undefined names beside it** (§3): a sensitivity row quoting a worst |SMD| computed over a
  covariate set some of which it could not judge is a row that reads better than the evidence. When
  those specifications arrive they call `assess` once each, with the `Propensity` from their own
  named specification, and the `role` column moves without the table changing shape (§4.1).
- **Stage 12 [§14a]** — calls **`smd` directly**, with unit weights, for the support check's baseline
  table comparing IVT-treated patients with never-IVT-centre patients. It does **not** call `assess`:
  that function reads a `Propensity`, and [§14] fits no propensity model anywhere. What it inherits
  is the yardstick and the four undefined branches; what it owes is its own population and its own
  denominator [§11], because a comparison across the restriction-1 boundary is not this stage's 92.
- **Stage 14 [§16]** — reports the balance table and the overlap table, and **reads the `status`
  column before drawing anything**: [§9] forbids an overlap plot for a structurally non-positive
  centre, and §6.2's column is how that reaches the reporting layer rather than being remembered. It
  matches on §6.2's three strings and needs no fourth, because the pooled row is `Balance.pooled` and
  not a member of `centres` (§3). A panel per declared centre is `for row in bal.centres`, with no
  row to drop and no literal to filter on.
  It also owes the [§16] baseline table, which is where the arm means §7.3 keeps out of this stage's
  tables belong.

**One rule every one of them owes, and it is Stage 6 §9's, unchanged.** `e` and `w` carry `nan` off
`in_model`; range over `in_model`, never over `notna()`, and never fill. Stage 7 is the first
consumer of that rule and it honours it in one place — `sub`, `a` and `w` are bound from
`ps.in_model` in §8 and nothing else reads `ps.w` — which is why the rule appears here as a fact
about §8's first three lines rather than as a warning.

## 10. What Stage 7 amends in Stages 1-6

The full ledger, so that no amendment is discovered during implementation. **Two files, and neither
is a shipped module of an earlier stage.**

| File | Amendment | Why |
|---|---|---|
| `config.py` | two computed views in the covariates block: `NEGATIVE_CONTROLS` and `BALANCE_SET` | §4.3. Stage 1's rule is that a covariate list lives in `config.py` and that none may be written out as a second literal. Both are computed from `PS_COVARIATES`, `PS_COVARIATES_FULL` and `BALANCE_ONLY`, which already exist; they live in `config.py` rather than in `balance.py` because **Stage 12 needs `NEGATIVE_CONTROLS` too** and would otherwise import a module named for a diagnostic in order to build a baseline table |
| `test_config.py` | four assertions: `NEGATIVE_CONTROLS` is exactly `PS_COVARIATES_FULL` minus `PS_COVARIATES`; it is disjoint from `PS_COVARIATES`; `BALANCE_SET` is `PS_COVARIATES + BALANCE_ONLY` with no duplicate; and `BALANCE_SET` is disjoint from `POST_TIME_ZERO` | §4.3, §4.5 B4. The last one is what makes B4 unreachable on the declared set, and asserting it is how that stays true |
| `implementation_roadmap.md` | Stage 7 gains its `**Spec:**` line, three **Accept when** items and one correction | §19; lands with this document |

**Five things that look like they need amending and do not.**

- **`data.py`.** §7.1: no new kind, no new heading, and therefore none of the four literal pins in
  the existing test modules moves. This is the first stage since Stage 2 to add audit entries without
  touching `data.py`, and it is what Stage 6 §10 predicted.
- **`propensity.py`.** `ess` is already public and already documented as Stage 7's (Stage 6 §6.2), and
  `Propensity` already carries everything this stage reads. Nothing is added to it.
- **`model.py`.** Not imported (§4.2), so not amended.
- **`TODOS.md`'s shared per-centre row-builder.** That item names its own trigger: "Do it when Stage
  7 or Stage 14 first wants the table under a non-`missingness` heading." **Stage 7 does not pull
  it.** The item is about `absence_by_column`'s four-way `structural`/`not recorded`/`missing`/
  `complete` classification and its per-centre *absence counts*; §6.1's table has one row per centre
  rather than one column, and its cells are propensity ranges, effective sample sizes and weight
  shares — not absence counts under a different heading. There is no loop here that
  `absence_by_column` could have supplied. The item stays filed, with its trigger unchanged, and this
  paragraph is the record that it was checked rather than forgotten.
- **`config.SMD_THRESHOLD`.** Declared at Stage 1 for this stage, unread until now, and read rather
  than written.

## 11. Handover to Stage 8

```
  cohort = cohort.build(eligibility.classify(derive.derive(*data.load()), ...), ...)
  ps     = propensity.fit(cohort, audit)
  bal    = balance.assess(cohort, ps, audit)

    19 covariate rows over 14 declared names       BALANCE_SET       [§4.1]
    92 records judged, every row's n = 92          in_model          [§4.4]
    worst |SMD| before  1.343 on center = HUG                        [§18]
    worst |SMD| after   0.380 on diabetes, a negative control        [§6.4]
    worst |SMD| after, in-model  0.211 on center = HUG   ABOVE 0.10  [§6.4]
    5 of 19 rows at or above SMD_THRESHOLD after weighting           [§18]
    0 rows undefined on this workbook                                [§5.3]
    3 centres reported, 1 structurally non-positive                  [§6.2]
    centres is 4 rows; pooled is a field beside it                   [§3, §6.1]
```

Stage 8 receives the cohort and the `Propensity`, exactly as Stage 6 handed them over. **It does not
read the `Balance`**: no estimate in this SAP is conditioned on a balance diagnostic, and §14 records
that the decision of what an unbalanced covariate implies is the PI's and not the estimator's.

## 12. Acceptance criteria

Tests live in `test_balance.py`, with section banners matching these numbers, as `test_cohort.py`,
`test_model.py` and `test_propensity.py` do. Tests needing the private workbook reuse the
`DATA_GATED` idiom and this file's **own** module-scoped `workbook` fixture; everything else runs on
a plain checkout with no `data/`.

### 12.0 The frames and fixtures this stage is tested on

Three kinds of input, and the split follows Stage 6 §12.0's:

**1. Hand-built vectors, for `smd`.** Plain `numpy` arrays with hand-computable answers. They carry
the roadmap's acceptance criterion and all four of §5.3's branches, and they need no frame:

```python
def hand_smd():
    """§12.3 — the roadmap's hand-computed value. sd = sqrt((0.5 + 0.5)/2), (3.5 - 1.5)/sd."""
    return (np.array([1.0, 2.0, 3.0, 4.0]), np.array([0.0, 0.0, 1.0, 1.0]), np.ones(4))


def separating():
    """§12.4a — constant within each arm, different between: the branch pilots returns 0.0 for."""
    return (np.array([0.0, 0.0, 1.0, 1.0]), np.array([0.0, 0.0, 1.0, 1.0]), np.ones(4))


def constant_everywhere():
    """§12.4a — the SAME constant in both arms, and NOT an indicator, so §3.1's 1-ulp fact bites.

    14.0 rather than 1.0 deliberately: np.average of a vector of exact 1.0s returns exactly 1.0
    under any weights, so an indicator cannot exercise §5.3a and 17 of 19 fixture rows are
    indicators.

    **The control arm's weights are (0.2, 0.8) and that is not free either.** What makes the two
    weighted means differ is whether an arm's Σ(w·x)/Σw division ROUNDS, which is a property of the
    weights and not of the arm size. Measured (§18): Σ[0.1, 0.9] and Σ[0.3, 0.7] both divide 14.0
    exactly, while Σ[0.2, 0.8] gives 14.000000000000002. An earlier draft of this fixture used
    (0.1, 0.9) against (0.3, 0.7) — two vectors that are both exact — so the two means agreed to
    every bit, §12.4a's rejected-form companion returned 0.0 rather than nan, and the fixture passed
    while exercising nothing. A constant covariate under exact weights cannot reach this branch.
    """
    return (np.full(4, 14.0), np.array([0.0, 0.0, 1.0, 1.0]),
            np.array([0.2, 0.8, 0.3, 0.7]))


def one_armed():
    """§12.4 — a covariate observed in only one arm, after the caller's own mask. Branch 1 and 3."""
    return (np.array([1.0, 2.0, 3.0]), np.array([0.0, 1.0, 1.0]), np.ones(3))
```

**2. `test_cohort.cohort_frame()`'s cohort**, taken through `derive → classify → build → fit`. Five
records, three bridging and two control, at two centres. Measured (§18) and it exercises more of this
stage than the workbook does:

- **14 of its 19 rows read exactly 0.0**, because the fixture's records are built from one template,
  so most covariates are constant across the whole cohort — §5.3's legitimate-zero branch, fourteen
  times.
- **`nihss_baseline` is the row §5.3a is about**: 14.0 in both arms, and the pre-correction form
  reported it undefined.
- **Every centre has exactly one control**, so **every within-centre SMD is undefined** while the
  per-centre ESS and propensity range are not — §6.3's arm-of-one case, which reads `missing` in one
  column of a row whose status is `overlap reported`. `no weighted contrast` is **not** reached here
  and is not reached on the workbook either (§6.3); §12.7 constructs it.
- Two of its four declared centres have no cohort record at all, so §6.2's structural row is
  exercised twice.

`test_balance.py` imports these from the modules that declare them and never re-declares them
(Stage 3 §12's rule). The modules, named because there are two of them and the natural guess is one:

```python
from test_cohort import cohort_frame            # tests/test_cohort.py:116
from test_data import hand_source, run          # tests/test_data.py:98, :107
```

`hand_frame()` and `tests/fixture_schema.xlsx` cannot reach this stage — Stage 5's P2 raises first,
measured (§18) — exactly as at Stage 6.

**3. The workbook, data-gated.** §12.12 only, through §12.0.3's fixture.

### 12.0.1 The bare-versus-normalised rule still applies

Stage 5 §12.0.1's distinction carries forward unchanged: `cohort_frame()` returns raw centre codes
and `run(...)` maps them to `CENTER_ORDER`'s labels. Unless a section says otherwise, every frame
below is the normalised one. A test that means to reach `assess` and passes a bare frame never gets
there — Stage 5's C4 raises first — and the failure looks like a passing `pytest.raises` whose
message is never read.

### 12.0.2 The golden balance vector

Pinned in `test_balance.py` as literals, from the synthetic frame so that no patient-derived number
enters git (Stage 6 §7.3). Measured (§18), on `cohort_frame()` normalised and taken through
`derive → classify → build → fit → assess`:

```
  row                       role                unweighted    weighted
  ─────────────────────────────────────────────────────────────────────
  core_ml                   propensity model      0.723747    0.635113
  tmax6_ml                  propensity model      2.248300    1.839089
  center = HUG              propensity model      0.258199    0.139181
  center = CHUV             propensity model     -0.258199   -0.139181
  penumbra_ml               excluded [§6]         1.630961    1.318924
  the other 14 rows                               0.000000    0.000000

  19 rows, 0 undefined, worst |weighted| 1.839089 on tmax6_ml
```

Asserted to 1e-6, for §12.0.2 of Stage 6's reason: this is a regression pin against an edit to our
own arithmetic, and pinning it at machine precision would make it fail on a numpy patch release
rather than on a mistake. **The assertion is by covariate name, one row at a time** — `center = HUG`
and `center = CHUV` differ only in sign, so an assertion against a sorted list of values would pass
while confusing them.

### 12.0.3 The workbook fixture returns three things from ONE audit

Module-scoped and this file's own, as Stage 5's and Stage 6's are — never shared across modules. Its
shape is specified rather than left to the pattern, because **§12.8 reads `overlap_weights` out of the
audit and compares it with Stage 7's own entries, so both must be in the same `Audit`**:

```python
@pytest.fixture(scope="module")
def workbook():
    """(df, ps, audit) from ONE linear run against ONE Audit. §12.0.3."""
    df, audit = data.load(data.WORKBOOK)
    df = derive.derive(df, audit)
    df = eligibility.classify(df, audit)
    df = cohort.build(df, audit)
    return df, propensity.fit(df, audit), audit
```

`test_propensity.py:593-596`'s `workbook_ps` deliberately does the opposite — it refits against a
**fresh** `Audit`, because that module's tests are about `fit`'s own four entries and a shared audit
would make their positions depend on what ran before. Copying that idiom here would put
`overlap_weights` in one `Audit` and `balance_smd` in another, and §12.8's reconciliation would have
nothing to read. Named here because the two fixtures look alike and differ in the one way that matters.

The audit this fixture returns carries **22** entries before `assess` and **24** after (§9, §18), which
is also what §12.9's captured-index assertions are counted against.

### 12.1 The balance set and the role column

- `BALANCE_SET` is the 14 declared names in `PS_COVARIATES + BALANCE_ONLY` order, and the table has
  **19** rows: **14** from `PS_COVARIATES`' expansion — 7 linear plus 3 `onset_type` levels plus 4
  `center` levels — and **5** from `BALANCE_ONLY`. Asserted against the declarations, never against a
  literal 19, and 14 + 5 is the decomposition an earlier draft wrote as 12 + 5 (§4.2).
- `_role` partitions the set: 14 rows `propensity model`, 4 `negative control`, 1 `excluded [§6]`,
  and the assertion is that the three counts sum to the row count with no row unlabelled.
- **`penumbra_ml` is asserted NOT to be a negative control**, by name, with the test naming [§6]'s
  reason. This is the roadmap correction of §19 in test form.
- Adding a fifth name to `PS_COVARIATES_FULL` is shown to move it into `NEGATIVE_CONTROLS`
  automatically, by patching the config value and rebuilding — so §4.3's "cannot fail to" is measured
  rather than described.

### 12.2 Declared levels, including the reference and the absent one

- Every declared level of every `CATEGORICAL` covariate has a row, parametrised over `CATEGORICAL`
  so a third factor is covered by the existing test.
- **`center = HUG` and `onset_type = witnessed` are present**, and the test names them as the
  reference levels `model.design` drops (§4.2).
- **`center = USZ` is present and reads 0.0 in both columns** on the workbook — a declared level
  nobody has, rendered rather than vanishing, via §5.3's legitimate-zero branch.
- A local reimplementation over `model.design(sub, C.BALANCE_SET)[0].columns` is shown to produce
  **16 rows rather than 19**, missing exactly `center = HUG`, `onset_type = witnessed` and
  `center = USZ` — and the test asserts that the row it misses carries the workbook's worst residual
  imbalance. This is the companion that makes §4.2 a measurement rather than an argument, and
  `BALANCE_SET` is the list it is built over for §4.2's reason: the same fourteen declared names, so
  the entire difference is the encoding. Measured (§18).
- A second assertion pins the worse form the pilots actually use (§16): over
  `model.design(sub, C.PS_COVARIATES)[0].columns` the table is **11 rows**, missing those three *and*
  all five `BALANCE_ONLY` names — eight of nineteen. Both numbers are asserted, because the first
  isolates the encoding and the second measures the shortcut that exists in code today.
- **`_levels`' `mask(df[name].isna())` is unreachable from `assess`, and the test says so rather than
  pretending otherwise.** Measured: `[c for c in C.CATEGORICAL if c not in C.PS_COVARIATES]` is empty —
  both declared factors are propensity covariates — so any record with an absent factor value is
  complete-cased out of `in_model` and never reaches `_levels`. Running `assess` against a `_levels`
  with the mask **deleted** gives byte-identical tables on the unmodified fixture, on a fixture with
  one `center` set to `pd.NA`, and on one with `onset_type` set to `pd.NA`. So the record's absence
  from every row's `n` comes from `in_model`, not from the mask, and an earlier draft's bullet — "that
  record contributes to no level's `n` rather than contributing a zero to all four" — was true of the
  output and false about the cause. The test asserts the **unreachability**: that every `CATEGORICAL`
  name is in `PS_COVARIATES`, which is what makes it so, and which fails if a third factor is declared
  outside the propensity model. The mask stays — it is correct, it is what Stage 12 will need over a
  population `complete_cases` did not build, and a guard that costs one method call is not worth
  removing to make a test honest.
- Note for the implementer, all measured on pandas 2.3.3 so they are not re-derived:
  `(string_series == level)` is nullable `boolean`, `.astype(float)` carries `<NA>` to `nan`, and
  `df.loc[<boolean with NA>]` treats `<NA>` as False. So §4.2's `_levels` and §8's `_centre` need no
  extra guard for an absent factor even in the population where one could appear.

### 12.2a The per-row denominator [§11]

The column §4.4 argues for, tested on a constructed frame because no available frame varies it — every
workbook row reads 92 and every fixture row reads 5 (§18).

- On `cohort_frame()`'s cohort with **one `in_model` record's `penumbra_ml` blanked**: that row's `n`
  is **4** against every other row's **5**, and its weighted SMD is **1.037070**. Measured (§18).
- **`propensity.fit`'s `in_model` is asserted UNCHANGED** by the same edit — a `BALANCE_ONLY`
  covariate never entered the fit, so it cannot move who was weighted. This is half the test: without
  it, an implementation that complete-cased on `BALANCE_SET` would pass the first assertion by
  dropping the record from every row.
- The row's SMD is asserted **equal to `smd` called directly** on the masked `(x, a, w)` triple, so the
  row is shown to be computed over its own records and its own pooled SD rather than over the table's.
- The contrast: blanking a `PS_COVARIATE` instead moves `in_model` **and** every row's `n` together.
  That is what makes the first three assertions mean something rather than merely pass.

Without this section the `n` column is a constant on every frame in the suite, and an implementation
that masked globally instead of per row would be green.

### 12.3 `smd` reproduces hand-computed values [roadmap]

- `hand_smd()` gives `2.8284271247461900`, asserted against `2 * np.sqrt(2)` to 1e-12 — the
  roadmap's acceptance criterion, computed rather than pinned.
- Two further hand values with **non-unit weights**, so the numerator's weighting is exercised
  independently of the denominator, and one asserting that the denominator does **not** change when
  the weights do — **read off `_pooled_sd` and off `CovariateBalance.sd` directly** (§3, §5.5), not
  reconstructed from two returned ratios as an earlier draft had to. That is [§9]'s "the yardstick does
  not move" as an assertion on the number itself rather than on an algebraic consequence of it.
- `_pooled_sd` is asserted to use **ddof = 1**, against a hand-computed variance, because §16b's oracle
  comparison turns on it: a ddof = 0 convention differs by `sqrt(n/(n-1))`, which at n = 92 is 0.5% and
  is visible at six significant figures while being invisible to the eye.
- `smd(x, a, w)` is asserted to equal `_ratio(x, a, w, _pooled_sd(x, a))` exactly, so §5.5's split
  cannot drift: the public function is the composition and nothing else.
- The unweighted column is asserted to equal `smd(x, a, np.ones_like(w))` exactly, not approximately:
  it is the same function and the same call.

### 12.4 The four undefined branches, each distinguished from the others

- Each of §5.3's four conditions returns `nan`, parametrised over the four with the branch named in
  the test id.
- **`one_armed()` returns `nan` and not `0.0`** — the roadmap's second acceptance criterion — and a
  companion runs `pilots/analysis.py`'s form on the same input to show it also returns `nan` there,
  so the test is honest about which of the two criteria the pilots already met.
- **An arm whose weights all sum to zero is shown to raise `ZeroDivisionError` from `np.average`
  when the guard is removed**, and to return `nan` with it. §3.1's fact, measured rather than
  asserted in a comment.
- An arm of one record returns `nan`, and the same data with one record added to that arm returns a
  number — so the branch is shown to be about the arm's size and not about the data's content.

### 12.4a The zero-variance split, both halves

- **`separating()` returns `nan`**, and the companion asserts `pilots/analysis.py`'s form returns
  **`0.0`** on the identical input. The test's comment states what that means: a covariate on which
  every treated patient scores 1 and every control 0 reported as the best-balanced row in the table.
- **`constant_everywhere()` returns `0.0`**, and a companion implementing §5.3a's rejected form —
  comparing the two weighted means — is shown to return `nan` on it. Both halves are required: the
  first alone passes for the wrong estimator, the second alone passes for one that never returns
  0.0 at all.
- The 1-ulp fact itself is asserted directly: `np.average(np.full(3, 14.0), weights=w1)` and
  `np.average(np.full(2, 14.0), weights=w0)` differ, and the same expression on a vector of `1.0`s
  does not. This is the test that explains the other two.
- On `cohort_frame()`'s cohort, `nihss_baseline` reads **0.0** and is asserted to be defined — the
  real-frame witness for §5.3a, running with no `data/`.

### 12.5 The threshold is a verdict and `assess` does not raise on it

- With a frame contrived so that every covariate is grossly imbalanced, `assess` returns normally and
  `Balance.unbalanced()` names the rows. No exception type is raised anywhere on that path.
- `unbalanced()` excludes undefined rows, and a companion asserts that a row whose weighted SMD is
  `nan` appears in `worst()`'s second element instead — the two collections partition the table with
  the defined-and-balanced rows.
- **`balance_smd`'s `detail` names exactly the rows `unbalanced()` returns**, asserted by comparing the
  two — one threshold predicate with two callers (§3). An earlier draft wrote the comparison twice, and
  the failure it invites is a `detail` string naming four rows above a table flagging five.
- `Balance.worst()` returns **both** the worst row and the undefined names from one call, and the
  test asserts the pair rather than the maximum: on a table with one undefined row, the worst is the
  largest **defined** |SMD| and the undefined name is returned beside it. A companion computes
  `pd.Series(...).abs().max()` and `np.abs(...).max()` over the same column and shows they give
  0.2-style skipping and `nan`-style poisoning respectively (§3.1) — the two wrong answers `worst()`
  exists to make unavailable.
- The threshold is read from `C.SMD_THRESHOLD` and never from a literal `0.10`, asserted by patching
  the constant and watching the verdict column move.

### 12.6 The within-centre table

- `Balance.centres` is asserted to be exactly `len(C.CENTER_ORDER)` rows in `CENTER_ORDER`'s order,
  and `Balance.pooled` to be a `CentreOverlap` **not** among them (§3). The rendered table is asserted
  separately to be those rows plus `all (pooled)`, in that order — the return value and the table have
  different shapes on purpose and both are pinned.
- `in_cohort` and `weighted` differ for at least one centre on the workbook (Lugano, 31 against 30),
  and the test asserts they are read from `at` and `at & in_model` rather than from each other.
- The `e` range per arm is asserted to lie inside the pooled range, and the pooled range to equal
  Stage 6's `propensity_fit` detail — one number in one place (Stage 6 §20.2's rule).
- **The `weight_share` cells that are DEFINED sum to 1.0 to 1e-12 over `Balance.centres`, and
  `Balance.pooled`'s share is exactly 1.0.** Both qualifications are load-bearing and an earlier draft
  had neither: a structurally non-positive centre's share is `nan`, not `0.0` — `len(w)` is 0, so §8
  returns `np.nan` and `_fmt` renders `missing` — and a `no weighted contrast` centre carries a
  positive share while not being a *reported* centre, so "over the reported centres" is the wrong
  index set. Measured on the fixture: the defined shares sum to exactly 1.0 with two structural rows
  reading `missing`. The pooled row is exact because §8 divides a sum by itself.
- The worst within-centre abs SMD is computed over `BALANCE_SET` **minus `center`**, and a companion
  shows that including the centre indicators changes no centre's worst on the workbook (§18) — so
  the exclusion is measured to be inert on this data and is kept for its reason rather than its
  effect.
- **§5.5's one yardstick, three assertions.** A centre row's `sd` for a covariate equals the pooled
  row's for the same covariate, exactly and by name. A companion recomputing the centre's worst against
  the centre's **own** variances is shown to differ — measured on the workbook: HUG +0.4956 against
  +0.5385, CHUV +1.4100 against +2.0608, Lugano −1.1087 against −1.5275 — so the choice is pinned by a
  test rather than by a comment. And the covariate *identified* as worst is asserted to be the same
  under either denominator (`diabetes`, `onset_type = unwitnessed`, `hyperlipidemia`), which is what
  makes the change one of scale and not of finding.
- Every label in a within-centre `_table` call is asserted to be a key of `sds`, so `_table`'s
  `sds.get(label, np.nan)` fallback is measured inert rather than assumed so (§8).

An earlier draft carried a companion here asserting that the defined shares sum to **less** than 1.0
when a record's `center` is `pd.NA`. **It is deleted, because the premise is impossible**: `center` is a
`PS_COVARIATES` member, so complete-casing removes such a record before Stage 7 sees it. Measured:
`in_model` drops from 5 to 4 and the defined shares still sum to exactly 1.0. There is no frame on which
a centre-less record is weighted, and a test asserting otherwise fails for the right reason on the wrong
premise, which is the worst kind of red. If the state is ever wanted it needs a hand-built `Propensity`
whose `in_model` includes such a record, in §12.7's declared "constructed, no frame reaches it" posture.

### 12.7 Structural non-positivity, and `ess` is never allowed to raise

- On the workbook, USZ has a row, its status is the structural one, and every overlap cell in it is
  `missing`. The test asserts the *status string*, because that is what Stage 14 reads (§6.2).
- **The set of status values over the whole table is asserted to be a subset of §6.2's three**, by
  comparing against the three module constants rather than against literals — so a fourth value cannot
  appear without failing here. `Balance.pooled.status` is asserted to be `_REPORTED`, which is what
  makes the three-value claim true of the rendered table and not only of `centres` (§6.2, §7.4).
- On `cohort_frame()`'s cohort, **two** centres carry the structural status and the other two carry
  `overlap reported` with a **`missing` worst |SMD|** — an arm of one, whose ESS and propensity range
  exist and whose balance cannot be judged (§6.3). The test asserts both cells of those rows, because
  asserting only the status would pass on an implementation that reported a worst of 0.0.
- **`no weighted contrast` needs a constructed frame and gets one**: blanking a [§6] covariate on a
  centre's only control empties that arm inside `in_model`, and the row is asserted to carry the
  status, a `missing` ESS for that arm, and a `missing` worst — with the centre still present in the
  table. Neither the workbook nor the fixture reaches it (§6.3), which is why it is constructed here
  rather than found.
- **`assess` does not raise on a centre whose arm `in_model` emptied**, and a companion asserts that
  `propensity.ess` on that same arm **does** raise `model.FitError` with the empty-arm message. Both
  halves: the guard is shown to be doing the work, and Stage 6's raise is shown to be intact behind
  it.
- `ess` is asserted to be called rather than reimplemented, by an AST scan asserting `balance.py`
  contains no division of a squared sum — with the companion proving the scan fires on a pasted Kish
  sum.

### 12.8 The pooled row reconciles with Stage 6

Every assertion here reads `Balance.pooled` and, on the other side, the `overlap_weights` entry pulled
out of the **same** `Audit` (§12.0.3) rather than recomputed.

**Every comparison against `overlap_weights` goes through `_fmt`, and that is not a formality.** An
`AuditEntry`'s `table` holds **strings** — Stage 6 built those cells with `_fmt` (`propensity.py:333-337`)
— so `Balance.pooled.ess[1] == float(cell)` compares a full-precision float against six significant
figures and fails on the seventh digit. The assertion is `_fmt(Balance.pooled.ess[code]) == cell`. §12.8's
first bullet says "exactly", and what is exact is the *rendering*: two entries of the same log agreeing to
the precision the log prints is the property a reader can check, and it is the strongest one available
without recomputing what Stage 6 already published.

- `Balance.pooled.ess` equals `Propensity.ess` for both arms, exactly, and equals the cells of Stage
  6's `overlap_weights` table read out of the audit entry **through `_fmt`**.
- `Balance.pooled.max_weight` equals `overlap_weights`' pooled `max w`, through `_fmt`.
- `Balance.pooled.worst` equals `Balance.worst()[0].weighted` in absolute value — one quantity,
  computed by the same `_table` call with the same arguments in two places (§8), which is the property
  this assertion exists to hold rather than to assume.
- `Balance.pooled.n` per arm equals `overlap_weights`' per-arm `n`, so the two entries cannot
  disagree about who was weighted.

### 12.9 The two audit entries, and the log

- The two `model` entries appear in `assess`'s declared order, at positions asserted against an index
  **captured before the call** — never a tail slice, which is the repair Stage 5 §10 landed and which
  Stage 6 §12.11 carried.
- **`data.KINDS` is asserted to be the declared nine, unchanged**, and `_HEADINGS` to render nothing
  new. This is §7.1's claim in test form, and it is the assertion that fails if someone adds a
  `balance` kind.
- Both entries render under `## Fitted models`, after `overlap_weights`.
- `balance_smd`'s `n` is the records judged and **not** the row count, and `overlap_by_centre`'s `n`
  is the count of reported centres — 3 on v7, with 4 rows plus a pooled row in the table. The test
  asserts the two are different numbers, because that is the property a reader needs.
- The verdict column contains only `yes`, `no` and `undefined`, and never `True`, `1` or `nan`.
- Every undefined SMD renders as the literal `missing`, and no cell contains a numpy repr: every
  float goes through `data._fmt`.
- **No cell of either table contains a `|`, and each table's rendered header, separator and body have
  the same markdown cell count.** `data._md_table` does no escaping and sizes the separator from
  `len(rows[0])` (data.py:242-245), so one pipe inside a header cell gives a header row with more cells
  than its own body — measured, before the columns were renamed: `balance_smd` rendered 8 against 6 and
  `overlap_by_centre` 15 against 13, in every log this stage writes. The companion pastes a `|` into a
  header and shows the counts diverge, so the assertion is known to fire. **Byte-identity does not catch
  this** — a table that is identically broken under both hash seeds passes §12.9's other criterion —
  which is why this is its own assertion and why it is the one §12.9 item that came from rendering the
  log rather than from reading the spec.
- The log renders byte-identically twice in one process, and under `PYTHONHASHSEED=0` versus `1` in a
  subprocess — the two-seed driver of Stage 5 §12.9, extended one stage.

### 12.10 Stage 7 adds nothing, refits nothing, and imports no fitter

- The frame passed in is returned unchanged, cell for cell, with the same columns, dtypes and index —
  the four assertions of Stage 5 §12.14 and Stage 6 §12.13, one stage on.
- The `Propensity` is unchanged: `e`, `w` and `in_model` are asserted equal to the objects passed in,
  and `assess` is asserted to add no attribute to the frame.
- An AST scan asserts `balance.py` imports neither `model` nor `statsmodels`, `sklearn` or `scipy`,
  and names neither `model.firth` nor `propensity.fit` anywhere. **The `model` half is the one that
  matters and is the reason the scan exists** (§4.2): importing it would not fail, it would silently
  produce a 12-row table. Companion test included for each.
- `assess` is asserted to take exactly `(df, ps, audit)` by `inspect` — no covariate list, no
  threshold, no centre list — so a defaulted keyword fails the test rather than passing silently.
- `balance.py` contains no bare `assert`, names no raw header, and is not in
  `EXEMPT_FROM_RAW_NAME_SCAN`, each with its scan-fires companion.

### 12.11 The preconditions

- B1 fires on a `Propensity` whose Series carry a different index, **raising `SchemaError` with B1's
  message**, and it is **parametrised over the three Series** — `e` alone, `w` alone, `in_model` alone —
  with the message asserted to name which one. The `e`-alone case is the one an earlier draft missed:
  measured, it reached `_centre`'s `ps.e.loc[weighted]` and raised pandas' `IndexingError` **after
  `balance_smd` had been recorded**, so the test also asserts that a failing `assess` appends **no**
  entry at all. Two companions: one shows that with the phase-1 raise removed the same input dies at
  B5's `df.loc[ps.in_model, C.TREATMENT]` with a bare `AssertionError('')` and B1's message discarded
  (§4.5a, measured); the other removes B1 entirely and shows the positional read producing a
  **complete, finite table** for the wrong records — the failure mode B1 exists for.
- B2 fires on a mask cast to `object` with a missing value in it, raising `SchemaError`; its companion
  shows that without the phase-1 raise the same input dies with pandas' `ValueError("Cannot mask with
  non-boolean array containing NA / NaN values")` instead (§4.5a).
- B3 fires on a frame missing a `BALANCE_ONLY` column, parametrised over `BALANCE_SET`, **and on a
  frame missing `C.TREATMENT`**; its companion shows that with B3 in phase 2 rather than phase 1 the
  `TREATMENT` case dies at B5 with `KeyError('ivt')` and B3's message discarded. That companion is the
  one that would have caught the incomplete first repair (§4.5a).
- **A frame failing B1 and B4 together reports only B1**, and the message says B4 and B5 were not run.
  Asserted, because the alternative reading — two of five assertions named on a frame that fails three
  — is a clean bill of health for the ones nobody checked (§4.5a).
- **A failing `assess` records nothing.** Asserted for every one of B1-B5 by capturing
  `len(audit.entries)` before the call and after the raise: a precondition that fires after an entry is
  written leaves a log describing half a diagnostic, which is DoD-9's failure arriving through the
  preconditions instead of through `ess`.
- B4 fires on a patched `BALANCE_SET` containing a `POST_TIME_ZERO` member, parametrised over the
  denylist — and the test records that it is unreachable on the declared set (§18), as Stage 6
  §12.6's D2 test does.
- B5 fires on a frame whose `in_model` leaves one arm, and its message is asserted to name both arm
  labels and their counts.
- Several at once produce one `SchemaError` naming each, as C1-C4 and D1-D4 do — **within a phase**.
  B3, B4 and B5 together give one error naming three; B1 and B2 together give one naming two; and the
  cross-phase case is the bullet above, which is a different assertion and not a weaker one.

### 12.12 Structural facts from the workbook [data-gated]

- 19 rows, every row's `n` is 92, and 0 rows are undefined.
- The worst unweighted |SMD| is `center = HUG` at 1.343, and the worst weighted is `diabetes` at
  0.380 — asserted to two decimals so a change in the estimator moves them.
- **Exactly 5 rows reach `SMD_THRESHOLD` after weighting, and exactly 2 of those are in the
  propensity model** — asserted as counts and by name. This is §6.4's finding pinned, so that a
  future workbook or a future specification cannot quietly change it.
- `center = USZ` reads 0.0 in both columns and its centre row carries the structural status.
- The overlap table's three reported centres carry the ESS, max weight and share of §18, each to
  four decimals, and their `weight_share` cells sum to 1.0 to 1e-12 with USZ's reading `missing`
  (§12.6).
- `statsmodels` is a declared project dependency (`pyproject.toml`, "retained for unpenalised
  cross-checks in tests"), so the MLE companion below needs no `importorskip` — unlike Stage 6's
  `firthlogist` oracle, which lives in the `reference` group. Checked rather than assumed.
- **The MLE companion**, data-gated: with `e` replaced by a `statsmodels.Logit` score on the same
  design, every in-model row's weighted |SMD| is below 1e-12 while the negative controls stay above
  0.1 — and under the [§7] Firth score the in-model worst is above 0.2. The pair is the assertion:
  the first half alone would pass on a silently reinstated MLE, which is the substitution Stage 6
  §12.10 exists to catch, and the second half is what makes this a test about Firth rather than about
  arithmetic.

### Coverage map

```
balance.py                                         test_balance.py
├── _pooled_sd / _ratio / smd                      written out in §5.2
│   ├── hand-computed value                        ├── 12.3  2*sqrt(2) to 1e-12  [roadmap]
│   ├── weighted numerator, UNWEIGHTED denominator ├── 12.3  sd read off the field, ddof = 1
│   ├── smd == _ratio(..., _pooled_sd(...))        ├── 12.3  the split cannot drift
│   ├── branch 1  arm of < 2 records               ├── 12.4
│   ├── branch 2  arm with zero total weight       ├── 12.4  + ZeroDivisionError without it
│   ├── branch 3  observed in one arm only         ├── 12.4  [roadmap]
│   ├── branch 4  sd == 0, arms differ             ├── 12.4a + pilots returns 0.0
│   ├── sd == 0, arms EQUAL → 0.0                  ├── 12.4a + the 1-ulp form returns nan
│   └── nan in x, nan in w → nan (not branches)    └── 12.4  §5.3's two recorded routes
│
├── _levels(df, name)
│   ├── one row per DECLARED level                 ├── 12.2  reference and absent levels present
│   ├── never model.design                         ├── 12.2  BALANCE_SET → 16 rows, PS_COV → 11
│   └── carries the factor's missingness           └── 12.2  + the three pandas facts checked
│
├── _role(name) / NEGATIVE_CONTROLS / BALANCE_SET
│   ├── three roles partition 14 names             ├── 12.1  14 + 4 + 1 = 19 rows
│   ├── penumbra_ml is NOT a negative control      ├── 12.1  [roadmap correction, §19]
│   └── computed, not declared                     └── 12.1  patch PS_COVARIATES_FULL
│
├── _table(sub, a, w, names, sds)
│   ├── 19 rows in BALANCE_SET order               ├── 12.1
│   ├── n is PER ROW, not the table's              ├── 12.2a CONSTRUCTED — no frame varies it
│   └── sds=None computes, a dict is imposed       └── 12.6  §5.5's one yardstick
│
├── _assert_balance_inputs   B1,B2,B3 | B4,B5      ├── 12.11 all five + several at once
│   ├── the phases are what let B1/B2/B3 raise     │         + AssertionError('') / ValueError /
│   │                                              │           KeyError companions, B1+B4 → B1
│   └── a failing assess records NOTHING           │           only, and e-alone misalignment
│
├── _centre(df, ps, centre, sds)
│   ├── centres is CENTER_ORDER; pooled is a field ├── 12.6  return shape vs table shape
│   ├── in_cohort != weighted                      ├── 12.6  Lugano 31 vs 30
│   ├── three statuses, and NO fourth              ├── 12.7  two on the fixture; the third is
│   │                                               │         CONSTRUCTED — no frame reaches it §6.3
│   ├── ess called, never reimplemented            ├── 12.7  AST scan + companion
│   ├── ess never allowed to raise                 ├── 12.7  + Stage 6's raise shown intact
│   ├── weight_share: DEFINED cells sum to 1       ├── 12.6  + nan for structural rows
│   ├── the POOLED yardstick, not the centre's     ├── 12.6  + the local-SD companion differs
│   └── center dropped WITHIN a centre             └── 12.6  measured inert on v7
│
├── _worst / _over_threshold  (Balance delegates)
│   ├── the pair, never the maximum alone          ├── 12.5  + pandas skips, numpy poisons
│   ├── undefined rows are in neither verdict      ├── 12.5
│   └── ONE threshold predicate, two callers       └── 12.5  _smd_detail's count == unbalanced()
│
└── assess(df, ps, audit)                          written out in §8
    ├── two `model` entries, in order              ├── 12.9  position vs captured index
    ├── KINDS is still the declared NINE           ├── 12.9  §7.1's claim
    ├── n is not the row count, in both            ├── 12.9
    ├── byte-identical, two hash seeds             ├── 12.9
    ├── NO cell holds a `|`; header == body cells  ├── 12.9  byte-identity does NOT catch it
    ├── pooled reconciles with Stage 6             ├── 12.8  ESS, max w, n, worst — one Audit
    ├── frame and Propensity unchanged             ├── 12.10
    ├── imports no fitter, and NOT model           ├── 12.10 AST scan + companion
    └── signature is exactly (df, ps, audit)       └── 12.10 by inspect

config.py amendments                               test_config.py
├── NEGATIVE_CONTROLS computed                     ├── four assertions (§10)
└── BALANCE_SET computed, disjoint POST_TIME_ZERO  └── the last is what makes B4 unreachable

R oracle                                           test_reference_r.py
└── PSweight's balance table                       └── §16b — numerators agree, denominators
                                                            differ BY DECLARATION

Data-gated (skipif not DATA_XLSX.exists()):        └── 12.12  19 rows, n = 92, 0 undefined,
  own module-scoped `workbook` fixture, §12.0.3               worst 1.343 / 0.380, 5 above 0.10
  (load → derive → classify → build → fit), ONE Audit        of which 2 in-model, the MLE pair

  smd 9 routes, _levels 3, _role 3, _table 3, _centre 10, worst() 3, assess 10, config 2  =  43
  43 with a test that can fail.  0 unreachable-by-construction, 3 unreachable ON EVERY FRAME
  THIS STAGE HAS — B4, disjoint by declaration (§4.5), tested by patching the set;
  `_NO_CONTRAST`, which needs an arm of zero and gets one only from a frame §12.7 builds
  (§6.3); and the PER-ROW `n`, which needs a BALANCE_ONLY covariate incomplete inside
  in_model and gets one only from the frame §12.2a builds (§4.4). All three are tested on
  constructed input rather than asserted unreachable, because all three are reachable on a
  workbook that has not arrived — and the third is reachable on a ONE-LINE edit to the frame,
  which is what makes leaving it untested the least defensible of the three.

  If the map and this count disagree, the map is authoritative and this block is stale.
```

### Failure modes

| Failure | Test | Error handling | Visible or silent |
|---|---|---|---|
| A covariate that perfectly separates the arms reported as perfectly balanced | 12.4a | §5.3 branch 4 returns `nan`; the verdict column reads `undefined` and `worst()` names it | visible **only** because of branch 4 — `0.0` and passing otherwise |
| A covariate identical in both arms reported as undefined, because two weighted means of one constant differ at 1 ulp | 12.4a | §5.3a compares the **arm constants**; measured on the committed fixture, where the rejected form did exactly this | visible in the table, and it was found by running it |
| The reference level of a factor — the workbook's worst-imbalanced covariate — absent from the balance table | 12.2 | §4.2 ranges over `FACTOR_LEVELS`; `balance.py` does not import `model` and a scan asserts it | **silent** without §4.2: the table is complete-looking and 16 rows long over `BALANCE_SET`, 11 over `PS_COVARIATES` |
| A declared level nobody has vanishing from the table, so the [§3] restriction cannot be read off it | 12.2 | the same declaration; `center = USZ` reads 0.0 in both columns | silent without it |
| The worst |SMD| reported over a table whose undefined rows were skipped | 12.5 | `worst()` returns the undefined names in the same call; `max(..., key=abs)` over the finite rows | **silent** in pandas, which skips `nan` in `max` |
| The worst |SMD| poisoned to `nan` by one undefined row | 12.5 | the same | visible, and the opposite error from the row above |
| A diagnostic aborting because complete-casing emptied one centre's arm | 12.7 | §6.3 checks the arm before calling `ess`; the row reports `no weighted contrast` | visible — it raised `FitError` before the guard |
| A structurally non-positive centre given an overlap plot [§9] | 12.7 | the `status` column, which Stage 14 reads before drawing (§6.2, §9) | visible in review |
| A balance table judged against the propensity model's covariates rather than the full [§6] set | 12.1 | `BALANCE_SET` is declared and `assess` takes no covariate list | visible |
| `penumbra_ml` reported as a negative control | 12.1 | `NEGATIVE_CONTROLS` is computed from `PS_COVARIATES_FULL`, which does not contain it | visible |
| An unweighted column over 93 and a weighted column over 92, divided by one SD belonging to neither | 12.12 | §4.4: both range over `in_model`, and every row carries its own `n` | visible |
| A misaligned `Propensity` producing a complete, finite table for the wrong records | 12.11 | **B1**, on the index and not the length | visible **only** because of B1 |
| B1 and B2 detected and then discarded, the caller seeing `AssertionError('')` from pandas instead of either message | 12.11 | §4.5's phase 1 raises before B3-B5 read through the mask; the companion asserts the bare `AssertionError` without it | visible — it raised either way, with **no message at all** |
| A `BALANCE_ONLY` covariate incomplete inside `in_model`, its row silently judged over the wrong denominator — or every row dropped to the narrowest one | 12.2a | §4.4's per-row mask, and §8's `_table` masking per row rather than per table | **silent** on every frame this stage has: the column reads a constant, so a global mask passes |
| Both audit tables rendered with a header row holding more markdown cells than its body, because a header cell contained a `\|` | 12.9, DoD-13 | the two columns are named `abs SMD < 0.1` and `worst abs SMD`, and §12.9 asserts no cell holds the character | **silent** past every other check: byte-identity across hash seeds passes on a table that is identically broken twice |
| The within-centre worst standardised by each centre's own variances, so no two rows of the column share a yardstick and none shares one with `balance_smd` | 12.6, DoD-14 | §5.5: `sds` is harvested from the pooled table and imposed on every centre row; `_centre` takes it with no default | **silent** — every cell is finite and plausible, and up to 46% out |
| A precondition firing after `balance_smd` is recorded, leaving a log describing half a diagnostic | 12.11 | B1 covers all three Series in phase 1, and §12.11 asserts a failing `assess` appends nothing | visible — it raised — but the log is already written |
| A fourth `status` value reaching Stage 14, which matches on §6.2's three | 12.7 | §6.2 has three and `_status` returns only those; the pooled row is `Balance.pooled`, not a status | visible in review only — Stage 14 would draw or skip on an unrecognised string |
| Residual imbalance above [§9]'s threshold read as a bug and repaired by substituting an MLE score | 12.12 | §6.4 measures it and §14 forbids it; Stage 6 §12.10's pair fails on the substitution | **structural** — no assertion in this stage can stop a reader, only inform one |
| A second effective sample size, or a second Kish sum | 12.7 | `propensity.ess` is called; an AST scan forbids the arithmetic | visible |
| Balance recomputed inside the bootstrap, producing a distribution of diagnostics | — | §2 and §9 state that Stage 7 is not on the bootstrap path | **structural, not asserted here** |

Every failure above is visible at the point it occurs, and **four** are visible only because of a
rule this stage adds: the separating covariate (§5.3 branch 4), the missing reference level (§4.2),
the skipped-or-poisoned worst (§3), and the misaligned `Propensity` (B1). The first two both return a
plausible number otherwise, and the first two are also the two the roadmap's single **Accept when**
turns out to have been about.

## 13. Known gaps carried forward

- **[§9]'s threshold is exceeded, by covariates in the propensity model, and this stage reports it
  rather than resolving it.** §6.4 is the measurement: `center = HUG` at 0.211 and `center = Lugano`
  at −0.154, against a threshold of 0.10, with the MLE companion establishing that the residual is
  the prescribed estimator's. **What is owed is a PI decision under [§9] and [§13]**, and the options
  are all [§13] amendments rather than code: report the residual with the estimate; add a [§13]
  sensitivity specification whose propensity model differs; or state in [§16] that the ATO population
  is balanced to 0.21 on centre and let the reader weigh it. Stage 7 does not choose, and §14 says
  why.
- **Three of the four negative controls are imbalanced after weighting**, worst 0.380 on `diabetes`,
  and one of them is worse after than before. [§6] anticipated this in words — "residual imbalance
  shows what weighting does not fix" — and the numbers are now attached to the sentence. They belong
  in [§16] beside the estimate, and the [§13] full-covariate sensitivity specification is where they
  stop being negative controls; when it runs, its balance table has the same 19 rows with four of
  them re-roled, which is what §4.1 is for.
- **The within-centre imbalances are large and this stage cannot say how much of that is noise.**
  Measured (§18), under §5.5's pooled yardstick: worst within-centre abs SMD is **0.50** at the largest
  centre, **1.41** at the second and **1.11** at the third, the last on 2 bridging patients against 28
  controls. A standardised mean difference over an arm of two is a number with almost no information in
  it, and this table reports it without a confidence interval because [§9] asks for the diagnostic and
  [§10]'s bootstrap is about estimates rather than about diagnostics. **A reader must not rank the
  centres by these numbers.** The `n` per arm is in the same row precisely so that they cannot be read
  alone. §5.5 removed one of the two reasons not to rank them — the yardstick no longer moves between
  rows — and this one, the arm sizes, is the larger and is unaffected.
- **No row is undefined on this workbook, so §5.3's four branches are exercised only by fixtures.**
  The zero-variance branch fires 14 times on the committed fixture cohort and the arm-size branch
  fires on every within-centre row there, so both have a real-frame witness — but neither has one on
  the workbook, and a future workbook with a covariate missing in one arm would reach branch 3 for
  the first time in production. §12.4 is what makes that a reported `missing` rather than a surprise.
- **The unweighted column's population was decided and its cost measured, not bounded.** §4.4: over
  the 92 rather than the 93, differing by at most 0.028 on v7. Nothing in this stage bounds that
  difference on the workbook that follows, and a workbook where twenty patients are covariate-
  incomplete would make the two columns describe visibly different populations while the table still
  reconciled. Stage 6 §13 files the same gap for the ATO population itself, and the two are one gap
  seen from two stages.
- **The R oracle is installed on this machine and currently skips, and the mechanism is not the one an
  earlier draft named.** Measured 2026-08-18 (§18): `logistf` 1.26.1 and `PSweight` 2.1.2 are present
  in `~/R/library` and load; `test_reference_r.py`'s gate reports them **missing**. The draft
  attributed that to `Rscript --vanilla` not reading `.Renviron` — **there is no `~/.Renviron` on this
  machine at all**, so nothing was being suppressed. What is actually true is simpler and worse:
  `R_LIBS_USER` is unset, R therefore falls back to its own compiled-in default
  (`~/R/x86_64-pc-linux-gnu-library/4.5`), **that directory does not exist**, and the packages sit in a
  directory R never looks at. So the draft's prescribed repair — "fall back to R's own default user
  library" — was a no-op: R already defaults there, with or without `--vanilla`, and there is nothing
  there. This is exactly the silently-skipping oracle Stage 6 §16b.2 was written to prevent, arriving
  by a route that section did not anticipate, and it was found only by running `Rscript` both ways and
  looking for the file. §16b records the repair that works.

## 14. What Stage 7 deliberately does not decide

- **What to do about the residual imbalance it measures.** §6.4 and §13. Reporting it, bracketing it
  with a [§13] specification, or accepting it are all PI decisions, and all three are [§13]
  amendments. The one option this stage removes is the silent one: after this document, the number is
  in the log and in the spec, and nobody can adopt the estimate without meeting it.
- **Whether an unbalanced covariate should change any estimate.** No estimator in this SAP is
  conditioned on a balance diagnostic, and Stage 8 does not read the `Balance` (§11). A pipeline in
  which a threshold silently selects an estimator is the mixture [§7] forbids, one level up.
- **Which centres get an overlap plot, beyond the ones [§9] forbids.** Stage 7 supplies the status
  and the numbers; how many panels a figure carries is Stage 14's [§16].
- **Whether the within-centre table should carry an interval.** §13. It would be a bootstrap of a
  diagnostic, and [§10]'s bootstrap is specified for estimates.
- **What the [§13] full-covariate specification's balance table will say.** Deferred by decision
  (roadmap Stage 11). §4.1 records only that it has the same 19 rows with four re-roled, so the two
  are readable against each other.
- **Whether `PSweight`'s SMD convention or [§9]'s is the better one.** §16b measures that they
  differ, exactly and recoverably. [§9] prescribes the unweighted pooled SD and gives the reason; a
  package choosing otherwise is a convention difference to record, not a disagreement to resolve.

## 15. NOT in scope for Stage 7

| Considered | Why deferred |
|---|---|
| Variance ratios per covariate (`pilots/analysis.py:249-265`) | §7.3. [§9] specifies the SMD and a threshold and says nothing about second moments; a column with no threshold beside it is interpreted privately by every reader |
| Arm means in the balance table | §7.3. Stage 6 §7.6 already carries them for 14 of the 19 rows, and [§16]'s baseline table is Stage 14's |
| Any figure — propensity histograms, love plots, per-centre overlap panels | Stage 14 is the single entry point for every figure [§16]. What this stage owes is §6.2's `status` column, which is what tells Stage 14 what it may draw |
| Bootstrapping the balance table | §13, §14. [§10]'s bootstrap is specified for estimates, and a per-replicate diagnostic has no null to be read against |
| Recomputing balance inside Stage 10's replicate loop | §2, §9. Balance is a property of the point specification |
| A `covariates=` or `threshold=` parameter on `assess` | §4.1, §5.4. [§9] fixes both, and a keyword makes a prespecified choice look like an option — Stage 6 §6.4's argument, one stage on |
| Substituting an unpenalised MLE score to obtain exact balance | §6.4, §14. [§7] prescribes one estimator and forbids a second; Stage 6 §12.10's pair fails on the substitution, which is what it is for |
| Trimming, truncating or re-normalising the weights to improve balance | Stage 6 §6.1 and [§15]: overlap weights are bounded by construction and this SAP has no trimmed row anywhere |
| A new `balance` audit kind and heading | §7.1. `data.py:144-149` declares `model` the last amendment this pipeline needs, and Stage 7 is the first stage to test it |
| Extracting the shared per-centre row-builder from `data.absence_by_column` | §10. Its `TODOS.md` entry names Stage 7 as a possible trigger and Stage 7 is not one: this table's cells are propensity ranges and effective sample sizes, not absence counts under a different heading |
| Removing the duplicate pooled `_table` pass | §8, §2. `assess` computes the 19 rows and `_centre(df, ps, None, sds)` computes them again. Sharing the result would break the one-code-path property §6.1 rests on and make §12.8's reconciliation unfalsifiable — a table that agrees with itself by construction has stopped being evidence [Stage 5 §7.3]. Priced at 7.4 ms |
| Escaping `\|` inside `data._md_table` | §7.4, §21.1. It is the durable fix and it is a **Stage 2** change: `_md_table` is what Stage 2's byte-identity criterion rests on, and amending it pulls four earlier test modules' pins into this stage's scope. Stage 7 renames its own two columns instead, and §12.9 asserts the character's absence so the next table cannot reintroduce it. Filed in `TODOS.md` |
| Bounding §4.4's 0.0284 cost on a future workbook | §13. The cost of the both-columns-over-`in_model` decision is measured on v7 and nothing here bounds it on the workbook that follows; a frame where twenty patients are covariate-incomplete would make the two columns describe visibly different populations. Bounding it needs the workbook that has not arrived |
| Promoting §18's per-centre table and §6.4's 66-94% to assertions | §12.12 asserts the ESS, max weight, share and the two worst |SMD| figures; the per-arm `e` ranges and the reduction percentages are verified by the §18 probe and by §21's re-derivation but are not in the suite. They are one loop over rows §12.12 already holds, and they are deferred rather than forgotten: §21.5 names them |
| A `pytest.importorskip("statsmodels")` on §12.12's MLE companion | §12.12. Checked rather than assumed: `statsmodels` is a declared project dependency, not a `reference`-group one, so the companion runs on any `uv sync`. Only `firthlogist` needs the guard, and that is Stage 6's |
| Reviewing §7.2's two `detail` strings as prose now | §21.5. Two review rounds have rendered the tables and read their cells; neither has checked that a `detail` sentence's counts match the table beside it or that its conditional clauses read as English at 0, 1 and many. Filed in `TODOS.md`, to run alongside T5 against real output rather than against this document |
| Distinguishing "removed by restriction 1" from "contributed no patient" in the centre status | §6.2. It would require the unrestricted frame as an argument, and Stage 5's `restrict_centres` entry already names which patients each restriction removed |

## 16. What already exists, and what to lift

`pilots/analysis.py` has both functions this stage needs — `smd` (194-208), `overlap_by_centre`
(210-247) and `balance_table` (249-265) — and the roadmap's standing warning applies. Here it splits
more finely than at Stage 6: **the formulae are right and are lifted; the undefined cases and the
column set are not.**

| From `pilots` | Status |
|---|---|
| `smd`'s numerator and its unweighted-pooled-SD denominator, and the docstring's reason ("the same yardstick is applied before and after weighting") | **lift both**, formula and comment. Verified (§18) to agree with this stage's `smd` on every row of both frames |
| `smd`'s `if not np.isfinite(sd): return np.nan` | **lift.** §5.3 branch 1 and 3, and the pilots got this one right — it is the roadmap's own acceptance criterion |
| `smd`'s `return … if sd > 0 else 0.0` | **do not lift.** §5.3a. It is branch 4, and it reports a perfectly separating covariate as perfectly balanced |
| `smd`'s `if w[t].sum() else np.nan` guard | **lift the guard, keep it before the mean.** §3.1: `np.average` raises rather than returning `nan`, so the order is load-bearing |
| `balance_table`'s covariate set, built by calling `design(df, covars)` | **do not lift.** §4.2. It is the design matrix, so it silently omits every reference level and every constant column — three of nineteen rows here, including the worst-imbalanced one |
| `balance_table`'s `variance_ratio` column | **do not lift.** §7.3, §15 |
| `balance_table`'s `balanced = abs(smd_weighted) < SMD_THRESHOLD` | **lift the comparison, not the dtype.** A boolean renders as `1` through `_fmt` (§7.2); the verdict is a three-valued string because an undefined row has neither passed nor failed |
| `overlap_by_centre`'s column list — n per arm, ps range per arm, ESS per arm, max weight, weight share, worst within-centre SMD | **lift in full.** It is roadmap Stage 7's list, and the pilots' is where the roadmap's came from |
| `overlap_by_centre`'s `for centre in df["center"].astype(str).unique()` | **do not lift.** It ranges over the data, so a centre absent from the cohort has no row and [§9]'s "reported as such" cannot happen. §6.2 ranges over `CENTER_ORDER` |
| `overlap_by_centre`'s `structural = n_t == 0 or n_c == 0`, and reporting it rather than plotting | **lift the concept.** §6.2 splits it into two statuses, because a centre [§3] removed and a centre complete-casing emptied are different facts |
| `overlap_by_centre`'s `judged = [c for c in covars if c != "center" and df.loc[k, c].nunique() > 1]` | **lift the `!= "center"`, decline the `nunique > 1`.** §8: dropping the centre indicators inside a centre is right and measured inert; dropping every covariate that is constant *there* silently narrows the judged set per centre, so the worst is a maximum over a different set at every centre |
| `overlap_by_centre`'s `round(…, 1)` and `round(…, 3)` on the reported cells | **do not lift.** `_fmt` is the pipeline's one formatter at six significant figures, and a second rounding is a second way for two runs to disagree (Stage 2 §11) |
| `overlap_by_centre`'s `ess(wk[dk == 1])` on a possibly-empty arm | **do not lift the call site.** §6.3: this repository's `ess` raises where the pilots' returns 0.0, and the guard belongs on this side of it |

`pilots/test_estimators.py`'s `test_firth_ato_balance_is_near_exact_but_not_exact` asserts
`0 < worst |SMD| < 0.05` on synthetic data. It is the origin of Stage 6 §6.3 and of §6.4 here, and
**its upper bound does not hold on this cohort** (0.211, §18). It is not lifted as a test; it is
recorded as the measurement that made the finding surprising.

## 16b. Validation against a reference implementation

Stage 6 §16b established the standing rule this stage inherits: *use a maintained library wherever
one exists; where none does, the estimator carries an oracle rather than a comment, and the oracle is
named in the spec.* For a standardised mean difference there is no Python library worth importing —
it is three lines — and there **is** an R one already installed, so this stage gets an oracle for
free.

**`PSweight` reports a per-covariate balance table on any supplied propensity score**, which makes it
an independent check on §5's arithmetic and on §4's covariate set at once. Stage 6's T6 installed it
(`PSweight` 2.1.2, `logistf` 1.26.1) and `test_reference_r.py` already carries the gate, the
environment fix and the CSV boundary. Stage 7 adds one script and reuses all of it.

**The convention differs by a DEFAULT, not by a limitation, and that changes what the oracle is worth.**
Measured 2026-08-18 (§18), on synthetic data through `SumStat(ps.estimate=…)`:

```
    PSweight's SMD  =  |mean_w(t) − mean_w(c)|  /  √( (SD_w,t² + SD_w,c²) / 2 )
    [§9]'s SMD      =   mean_w(t) − mean_w(c)   /  √( (SD_t²   + SD_c²  ) / 2 )

    same pooled form, WEIGHTED standard deviations against [§9]'s UNWEIGHTED ones,
    and PSweight reports the ABSOLUTE value where §5 keeps the sign.
```

**But `summary.SumStat` takes a `weighted.var` argument, and `weighted.var = FALSE` is [§9]'s
denominator.** Read out of the installed package (`~/R/library/PSweight/R/PSweight.rdb`, whose lazyload
database carries the literal strings `"ASD weighted var"` and `"ASD unweighted var"` alongside the
argument lists of `summary.SumStat(object, weighted.var, metric = "ASD")` and `plot.SumStat`). So an
earlier draft of this section asserted a convention *difference* where the package offers a switch, and
under-used the oracle as a result: **T6 runs the summary both ways** — `weighted.var = FALSE` for a
direct comparison against §5's weighted column, which is a strictly stronger assertion than
reconstructing a numerator, and `weighted.var = TRUE` for the recorded convention difference below.

**One thing T6 must assert rather than assume: PSweight's ddof under `weighted.var = FALSE`.** §5.2 uses
`ddof = 1`. A `ddof = 0` convention differs by `sqrt(n/(n−1))`, which at n ≈ 90 is about 0.5% — visible
at six significant figures, invisible to the eye, and enough to make an "agrees exactly" claim false in
the fourth digit. It is unresolved here because R would not run in the session that found this, and
§12.3 asserts our own ddof so that only one side of the comparison is ever in question.

Three consequences of the convention difference itself, all of them assertions rather than caveats:

- **The unweighted column agrees exactly**, because with `w ≡ 1` the weighted SD is the unweighted
  SD. Measured: PSweight's 0.420264 reproduced from its own printed means and SDs by §5.1's formula.
  So the oracle validates §5's arithmetic outright on half the table.
- **The weighted column's numerators agree**, and the test asserts that rather than the ratio:
  `|ours| × our_pooled_unweighted_SD == theirs × their_pooled_weighted_SD`. Measured on a design with
  a deliberately imperfect score, so the numerator is non-zero: ours 0.168377 × 1.023623 and theirs
  0.169416 × 1.017345 both give 0.172354. The denominators differ **by declaration** and the
  difference is reconstructible, which is Stage 6 §16b.2's "compare `e` and `w` first and the ESS
  separately" applied to a third quantity.
- **[§9] prescribes the unweighted denominator and gives the reason** — "so the yardstick does not
  move" — and PSweight recomputes the SD under each weighting scheme, so its before-and-after columns
  are divided by two different numbers. That is a legitimate convention and it is not this SAP's.
  §14 records that the oracle settles which convention each side uses and not which is better.

**One thing the oracle confirms that no Python check can.** Under `PSweight`'s own unpenalised `glm`
score, every in-model covariate's SMD is 0 to printed precision — the exact-balance property, from
the authors of the overlap-weight method, on their own implementation. §6.4's MLE companion measures
the same thing with `statsmodels`, and the two agreeing means §6.4's attribution of the residual to
Firth rests on two independent implementations rather than on one.

**And the gate must be made loud, because it is currently silent.** Measured 2026-08-18 (§18):

```
    ~/.Renviron                                       DOES NOT EXIST
    Rscript --vanilla, R_LIBS_USER unset  →  R_LIBS_USER resolves to
        ~/R/x86_64-pc-linux-gnu-library/4.5           DOES NOT EXIST  (absent from .libPaths())
    ~/R/library                                       81 packages, incl. PSweight + logistf
    Rscript --vanilla with R_LIBS_USER=~/R/library →  both load
    Rscript --no-init-file, R_LIBS_USER unset      →  same as --vanilla: still missing
```

So `--vanilla` is **not** the cause and there is nothing for it to suppress; the packages are simply in
a directory nothing points R at. **A fallback to "R's own default user library" is therefore a no-op**
— that default is exactly what R already uses, and it does not exist. An earlier draft of this section
prescribed that fallback and would have left the oracle skipping while claiming to have fixed it.

T6 fixes it by **searching** rather than inheriting:

```python
# test_reference_r.py — §16b
_R_USER_LIB_CANDIDATES: Final[tuple[Path, ...]] = (Path.home() / "R" / "library",)


def _r_user_library() -> str | None:
    """A user library that EXISTS, when R_LIBS_USER is unset. §16b.

    Measured 2026-08-18: R's own default (~/R/x86_64-pc-linux-gnu-library/<ver>) does not exist on
    this machine and the oracle's packages are in ~/R/library, so inheriting R's default finds
    nothing. There is no ~/.Renviron either — `--vanilla` suppresses nothing here and stays, because
    §16b.1 needs a run no user profile can change.
    """
    return next((str(p) for p in _R_USER_LIB_CANDIDATES if p.is_dir()), None)
```

`_r_environment` sets `R_LIBS_USER` from it when the variable is unset, and leaves an explicitly set
one alone — a caller pointing at a different library still wins. Two assertions, and the second is the
one that matters: the resolved path is the directory the packages are in, **and** `_r_has` returns True
for both packages under the environment `_r_environment` actually builds, so the gate is asserted open
rather than assumed open. When neither is true the skip reason names every candidate searched, so a
reader on a different machine is told where to install rather than that something is missing. That is a
change to Stage 6's file and §19 records it.

## 17. Implementation tasks

Ordered. **T1, T2, T3 and T4 are independently verifiable against their own privates; every acceptance
section that reads a `Balance` or an `Audit` belongs to T5**, and an earlier draft assigned §12.5 to T3
and §12.6-12.8 to T4 where none of them can run — §12.5 needs `assess` and a constructible `Balance`
(whose `pooled` field is a T4 deliverable), and §12.6-12.8 need the audit `assess` writes into. The
lanes below are unchanged; what moved is which task closes which acceptance section. T1 is the amendment
that unblocks the rest, T2-T5 are Stage
7's own, T6 is the R oracle, T7 is the log and the definition of done.

**T1 must land first** because `_role` and `_assert_balance_inputs` both read `C.NEGATIVE_CONTROLS`
and `C.BALANCE_SET`, so nothing below it imports. T2 and T3 are genuinely independent — `smd` knows
nothing about a frame and `_levels` knows nothing about weights — and are the one place in this stage
where two worktrees would not cost more than they save. **T6 depends on nothing but T2**: it reads
`smd`'s output through a CSV and changes no shipped code.

- [x] **T1 (P1)** — `config.py`: `NEGATIVE_CONTROLS` and `BALANCE_SET` in the covariates block, as
      §4.3 writes them. `test_config.py`: the four assertions of §10. Verify: `uv run pytest -v`
      green with **no other test edited**, which is where §10's claim that `data.py` needs no
      amendment is checked rather than believed — and where the four new assertions are confirmed to
      pass on the config as it stands, since it would be a landed-code defect if they did not.
- [x] **T2 (P1)** — `smd`, with all four branches of §5.3 and §5.3a's arm-constant comparison.
      Verify: acceptance 12.3, 12.4, 12.4a. **Two things to get exactly right**, each of which was a
      defect in an earlier draft of this document: the zero-total-weight guard comes **before**
      `np.average`, or the branch raises instead of returning (§3.1); and the `sd == 0` branch
      compares `x[t][0]` to `x[c][0]` and **not** the two weighted means, or a covariate identical in
      both arms is reported undefined (§5.3a, measured on the committed fixture).
- [x] **T3 (P1)** — `_levels`, `_role`, `_table`, `CovariateBalance`, `_worst` / `_over_threshold`, and
      the `Balance` methods that delegate to them. Verify: acceptance 12.1, 12.2, and **direct
      assertions on `_worst` / `_over_threshold` over a hand-built list of `CovariateBalance`**. Watch
      `_levels`' `mask(isna())`: it is correct and it is unreachable from `assess` while every
      `CATEGORICAL` name is a PS covariate (§12.2, measured) — keep it, and assert the unreachability
      rather than a behaviour no input produces. And watch `_table`'s mask: it is **per row** —
      `present = np.isfinite(values)` inside the loop — never one mask over the frame, which only
      §12.2a can tell apart because the `n` column is a constant on every unmodified frame (§4.4).
- [x] **T4 (P1)** — `_centre`, `CentreOverlap`, `_status`, the pooled row through the same function, and
      §5.5's `sds` threading. Verify: **direct assertions on `_centre` and `_status`** — one row per
      declared centre, the three statuses, the arm check before every `propensity.ess` call (§6.3), and a
      centre row's `sd` equal to the pooled row's. The arm check precedes every `ess` call; the pooled row
      is `_centre(df, ps, None, sds)` rather than a second code path; and it lands in `Balance.pooled`
      rather than in `centres`, with `_status` called for it like any other row — there is no `_POOLED`
      constant and §6.2 has exactly three values (§3, §6.2).
- [x] **T5 (P1)** — `assess` as written in §8, `_assert_balance_inputs` with B1-B5 **in §4.5's two
      phases**, and §7.2's two entries with their detail and table helpers. Verify: acceptance
      **12.2a, 12.5, 12.6, 12.7, 12.8**, 12.9, 12.10, 12.11 — everything that reads a `Balance` or the
      `Audit`. **The phase-1 raise is load-bearing, not stylistic**: B4, B5 and §8 read through the frame
      and the mask, so collecting all five into one raise means B1's, B2's and B3's messages are
      discarded by a bare `AssertionError('')`, a `ValueError` or a `KeyError` that pandas throws a check
      or two later (§4.5a, all three measured). Every other assertion helper in this repository collects
      and raises once, and a reader who tidies this one back into that shape reintroduces the defect —
      which is why §12.11 carries a companion for each of the three.
- [x] **T6 (P2)** — the R oracle: `tests/reference/balance_psweight.R`, the Stage 7 half of
      `test_reference_r.py`, and `_r_user_library`'s **candidate search** (§16b) — not a fallback to
      R's own default user library, which is measured to be a no-op: R already defaults there and the
      directory does not exist. Verify: green
      **with** the packages reachable — the comparison runs, the unweighted column agrees exactly and
      the weighted numerators agree — and green **without** them, where it skips and `pytest -rs`
      names every candidate searched. **`summary.SumStat` is called BOTH ways** (§16b): `weighted.var =
      FALSE` for a direct comparison against §5's weighted column, which is the strong assertion, and
      `weighted.var = TRUE` for the recorded convention difference. **Assert PSweight's ddof rather than
      assuming it** — a `ddof = 0` denominator differs from §5.2's by 0.5% at n ≈ 90, which is visible at
      six significant figures. The assertion that closes this task is that `_r_has` returns True
      under the environment `_r_environment` builds **with `R_LIBS_USER` unset in the parent shell** —
      that is the exact state the oracle was silently skipping in. **Not a completion blocker**, unlike
      Stage 6's T6: that gate existed
      because no independent check on the Firth kernel was available, and this stage's arithmetic is
      three lines with a hand-computed acceptance criterion. It is P2 for that reason and for no
      other.
- [x] **T7 (P1)** — the audit log produced through all seven stages against the workbook, read end to
      end by a human, and the definition of done below.

### What can be built in parallel, and what cannot

| Task | Modules touched | Depends on |
|---|---|---|
| T1 | `config.py`, `test_config.py` | — |
| T2 | `balance.py` (`smd`), `test_balance.py` | T1 |
| T3 | `balance.py` (the covariate rows), `test_balance.py` | T1 |
| T4 | `balance.py` (the centre rows, and §5.5's `sds` threading), `test_balance.py` | T2, T3 |
| T5 | `balance.py` (`assess`, the entries), `test_balance.py` | T3, T4 |
| T6 | `test_reference_r.py`, `tests/reference/balance_psweight.R` | T2 |
| T7 | none — it runs the pipeline and reads the log | T1-T5 |

```
  Lane A:  T1  →  T2  ┬→  T4  →  T5  →  T7
  Lane B:       T3   ─┘
  Lane C:       T6 ────────────────────┘   (joins only at the suite)
```

**T2 and T3 both touch `balance.py` and are still worth separating.** They share a file and no state:
`smd` takes three arrays and `_levels` takes a frame and a name, with no call between them until T4.
Merge T2 first — T3's `_table` calls `smd` twice per row and a conflict there means §5.1's
one-function rule was implemented twice.

### Diagrams and comments that belong in the code, not only here

Three, and no more — a diagram nobody maintains is worse than none, because it is believed.

- **`balance.py`'s module docstring** — §4.2's sentence: balance is judged per **declared level**, so
  this module does not import `model`, because `design` drops every reference level and every
  constant column and one of those carries this cohort's worst residual imbalance. That is the
  sentence that stops the next reader tidying the indicator loop into a `design` call.
- **Above `smd`'s zero-variance branch** — §5.3a's two lines: both arms are constant whenever `sd`
  is 0, so the comparison is of the raw arm constants; comparing the weighted means fails at 1 ulp
  and was measured doing so on the committed fixture.
- **Above `_centre`'s `ess` calls** — §6.3: the arm is checked first because `propensity.ess` raises
  on an empty arm by design (Stage 6 §6.2), and here an empty arm is a finding to report rather than
  an exception to propagate.

### Definition of done

Stage 7 is complete when all of the following hold, and not before.

1. `uv run pytest -v` is green with `data/` present.
2. Green **with `data/` temporarily renamed** — the data-gated tests skip and nothing else fails or
   errors at collection.
3. `test_config.py`, `test_data.py`, `test_derive.py`, `test_eligibility.py`, `test_cohort.py`,
   `test_model.py` and `test_propensity.py` are still green, **unedited except for
   `test_config.py`'s four new assertions**, and `balance.py` and `test_balance.py` are not in
   `EXEMPT_FROM_RAW_NAME_SCAN`.
4. The pipeline's audit log has been produced through all seven stages against the **workbook** and
   reproduces byte-identically across hash seeds:

   ```bash
   cd extended_bridging
   S='import sys, data, derive, eligibility, cohort, propensity, balance
   df, audit = data.load(data.WORKBOOK)
   df = derive.derive(df, audit)
   df = eligibility.classify(df, audit)
   df = cohort.build(df, audit)
   ps = propensity.fit(df, audit)
   balance.assess(df, ps, audit)
   sys.stdout.write(audit.to_markdown())'
   PYTHONHASHSEED=0 uv run python -c "$S" > /tmp/a.md
   PYTHONHASHSEED=1 uv run python -c "$S" > /tmp/b.md
   diff /tmp/a.md /tmp/b.md
   ```

   `diff` reports nothing, and the two new entries render under `## Fitted models` after
   `overlap_weights`. The log names patients and must never be committed.
5. **The zero-variance branch has been seen to matter, in both directions.** With §5.3's branch 4
   replaced by the pilots' `return … if sd > 0 else 0.0`, a covariate that perfectly separates the
   arms is shown to report `0.0` and to pass the threshold as the best-balanced row in the table.
   Then, with the branch present but comparing the two **weighted means**, `nihss_baseline` on
   `cohort_frame()`'s cohort is shown to report `missing` — a covariate identical in both arms,
   undefined. Watch this one specifically: it is the only defect in this stage that was found by
   running the code rather than by writing it.
6. **The design-matrix shortcut has been seen to fail.** With `_levels` replaced by a loop over
   `model.design(sub, C.BALANCE_SET)[0].columns`, the balance table is shown to have **16 rows rather
   than 19** and to omit `center = HUG`, which carries the workbook's worst residual imbalance — so the
   shortcut produces a shorter, complete-looking table whose reported worst is wrong. Nothing raises.
   Over `C.PS_COVARIATES` instead it is **11 rows**, omitting eight of nineteen. Both counts measured
   (§4.2, §18); neither is 12, which is `design_matrix`'s parameter count.
7. **`worst()` has been seen to be necessary.** On a table with one undefined row,
   `pd.Series(...).abs().max()` is shown to return a number that ignores it and
   `np.abs(...).max()` to return `nan`, and `worst()` to return the number **and** the name.
8. **B1 has been seen to fail, in both of its two ways.** With the assertion removed and a
   `Propensity` reindexed to a shifted index, `assess` is shown to return a complete table of finite
   numbers describing the wrong records, raising nothing. **And with the assertion present but §4.5's
   phase-1 raise removed**, the same input is shown to die at B5 with a bare `AssertionError('')` —
   B1's message assembled and discarded. The second is what the two phases exist for and is the one a
   later tidy-up would undo (§4.5a).
8a. **The per-row denominator has been seen to be per row.** With one `in_model` record's
    `penumbra_ml` blanked, that row's `n` is below every other row's and its SMD is the one `smd`
    returns for the masked triple; with `_table`'s mask hoisted out of the loop, every row's `n` drops
    together and the table still reconciles with itself (§4.4, §12.2a). The column is a constant on the
    workbook and on the fixture, so this is the only place either form can be told from the other.
9. **`ess`'s empty-arm raise has been seen behind the guard.** With §6.3's arm check removed and one
   centre's only control excluded from `in_model`, `assess` is shown to raise `model.FitError`
   mid-table, leaving `balance_smd` recorded and `overlap_by_centre` absent — a log describing half a
   diagnostic. With the guard, the row reads `no weighted contrast`.
10. **§6.4's finding has been reproduced and read.** The in-model worst weighted |SMD| is above 0.20
    on `center = HUG`, the MLE companion puts the same quantity below 1e-12, and both numbers are in
    the acceptance suite (§12.12). A reader who has not seen both has not seen the finding, and it is
    the one thing in this stage a manuscript depends on.
11. The audit log has been read end to end by a human, and its numbers agree with §18: the
    `balance_smd` entry's `n` equals `propensity_fit`'s `n`, and the pooled row of
    `overlap_by_centre` equals `overlap_weights`' pooled cells.
12. Run from `extended_bridging/`, `grep -nE "^\s*(import|from)\s+.*\bmodel\b" balance.py` returns
    nothing **and exits 1**, and `test_balance.py`'s AST scan for the same is green with its companion
    firing on a pasted import. Three things about that command, each of which was wrong in an earlier
    draft: **the path is relative to `extended_bridging/`**, because DoD-4 has already `cd`-ed there and
    `extended_bridging/balance.py` would then not exist — a grep that exits 2 on a missing file
    "returns nothing" and reads as a pass; **the pattern is `.*\bmodel\b` and not `model\b`**, so
    `import propensity, model` matches; and **it is import-shaped rather than `"model\."`**, because §17
    requires a module docstring explaining why `design` is not called and a text scan for `model.`
    collides with the one comment this stage is obliged to write. §12.10's AST scan is the authoritative
    check — it reads imports and attribute accesses rather than characters — and this grep is the cheap
    human-runnable version of it.
13. **Both audit tables render as well-formed markdown.** For each of the two entries, the rendered
    header, separator and body rows have the same number of markdown cells, and no cell of either table
    contains a `|` (§12.9). Measured before the columns were renamed: `balance_smd` rendered an 8-cell
    header against a 6-cell separator and `overlap_by_centre` 15 against 13, in every log this stage
    writes, and **byte-identity across hash seeds passed the whole time** — a table identically broken
    under both seeds is still identical. Look at the log, not only at its hash.
14. **The within-centre column is on the pooled yardstick.** A centre row's `sd` equals the pooled row's
    for the same covariate, and the companion computing the centre's worst against the centre's own
    variances is shown to differ (§5.5, §12.6): measured +0.4956 against +0.5385 at HUG, +1.4100 against
    +2.0608 at CHUV, −1.1087 against −1.5275 at Lugano.

## 18. Verification record

Checked on 2026-08-18, before this spec was finalised, by a read-only probe reporting aggregate
quantities only — no outcome by arm, in the posture Stages 0, 4, 5 and 6 used. The probe drives
`load → derive → classify → build → fit` against `data.WORKBOOK` and applies §4, §5 and §6; it writes
nothing. Acceptance 12.12 re-checks the workbook rows in code.

**This table is normative and every number elsewhere in this document points at it.** That rule exists
because of a measured failure of this document's own structure: the design-matrix width was written into
§4.2, §12.1, §12.2, this table and DoD-6 — five places — and when §20.2 established it was wrong, it had
to be corrected in five. The same shape of error produced the 12 + 5 = 19 decomposition, the "67%" lower
bound and the reversed per-arm ESS pair, all three found by later rounds rather than by the four earlier
readings that walked past them. **Twenty-eight hundred lines of prose around a hundred and fifty lines of
code means a load-bearing number has many homes and no address**, and the only cheap defence is that one
of the homes is declared authoritative. So: a number in prose cites the row here that produced it, this
table cites the code that produced it, and where the two disagree **this table is right and the prose is
stale** — the same precedence rule §12's coverage map states about itself. An implementer who finds a
conflict should fix the prose and not re-derive the number.

| Claim | Where used | Verified |
|---|---|---|
| `NEGATIVE_CONTROLS` is `('hypertension', 'hyperlipidemia', 'diabetes', 'smoking')` and `BALANCE_ONLY` minus it is `('penumbra_ml',)`, so the three roles partition the 14 names | §4.3, §12.1, §19 | yes, run |
| `BALANCE_SET` expands to **19 rows** over 14 declared names: **14** from `PS_COVARIATES` (7 linear + 3 `onset_type` levels + 4 `center` levels, the two factor names being replaced by their levels) and **5** from `BALANCE_ONLY`. An earlier draft of this row wrote 12 + 5, which is 17 | §4.1, §7.2, §12.1 | yes, run |
| **`design(sub, C.BALANCE_SET)[0]` is 16 columns and `design(sub, C.PS_COVARIATES)[0]` is 11**, both with `dropped == ('center_USZ',)`, on a frame retaining every declared level. So the design shortcut omits 3 of the 19 rows over `BALANCE_SET` and 8 over `PS_COVARIATES` — never 12, which is `design_matrix`'s parameter count `X.shape[1] + 1` | §4.2, §12.2, DoD-6 | yes, run |
| `BALANCE_SET` and `POST_TIME_ZERO` are **disjoint**, so B4 is unreachable on the declared set | §4.5, §10, §12.11 | yes, run |
| Every row's `n` is **92** on the workbook — `penumbra_ml` is absent on exactly the record the two CTP volumes are absent on — so §4.4's per-row denominator is a constant on this workbook and not in general | §4.4, §12.12 | yes, run |
| Worst |SMD| **before** weighting: **1.34303** on `center = HUG`. Worst **after**: **0.379832** on `diabetes`. Worst after among in-model covariates: **0.211121** on `center = HUG` | §4.2, §6.4, §11, §12.12, §13 | yes, run |
| **11 of 19** rows are at or above `SMD_THRESHOLD` before weighting and **5 of 19** after: `center = HUG` +0.211121, `center = Lugano` −0.153704, `hypertension` −0.139037, `hyperlipidemia` −0.300935, `diabetes` +0.379832 | §6.4, §12.12, §13 | yes, run |
| The four negative controls, before → after: hypertension +0.051189 → −0.139037; hyperlipidemia −0.390111 → −0.300935; diabetes +0.218912 → **+0.379832, worse after than before**; smoking −0.164367 → +0.045324 | §6.4, §13 | yes, run |
| `penumbra_ml` +0.078856 → +0.011317, and `center = USZ` **0 → 0** — a declared level nobody has, rendered rather than vanishing | §4.2, §5.3, §12.2 | yes, run |
| **The MLE companion.** With `e` from `statsmodels.Logit` on the same design (converged, max\|coef\| 4.546): every in-model row's weighted \|SMD\| ≤ **1.03e-15** and `Σw` is 12.9716 in both arms, equal to 1.8e-15; the negative controls stay imbalanced, worst **0.391258**. Under the [§7] Firth score the same quantities are 0.211121 and 0.379832 | §6.4, §12.12, §13, §16b | yes, run |
| Of the eight in-model rows above 0.10 before weighting, weighting removes between **66% and 94%** of the imbalance — minimum 66.5% on `center = CHUV` (0.2026 → 0.0679), maximum 93.6% on `onset_type = wake_up`; the two left above the threshold are the two that started above 1.2. An earlier draft rounded the lower bound UP to 67%, which is a claim the data does not support | §6.4, §12.12 | yes, run |
| Per centre, on the workbook, every per-arm pair written `{0: EVT alone, 1: bridging}` so no reader has to infer the order: HUG n=41, n `{0: 11, 1: 30}`, `e` `{0: [0.521465, 0.831072], 1: [0.496788, 0.931070]}`, ESS `{0: 10.6761, 1: 25.3441}`, max w 0.831072, share 0.566719; CHUV n=21, n `{0: 14, 1: 7}`, `e` `{0: [0.179797, 0.630395], 1: [0.393084, 0.627335]}`, ESS `{0: 12.2544, 1: 6.8295}`, max w 0.630395, share 0.275094; Lugano n=30, n `{0: 28, 1: 2}`, `e` `{0: [0.026972, 0.237312], 1: [0.035045, 0.224802]}`, ESS `{0: 19.8179, 1: 1.9765}`, max w 0.964955, share 0.158186; USZ **0 cohort records** | §6.1, §6.2, §12.6, §12.12, §13 | yes, run |
| **The worst within-centre abs SMD, under §5.5's POOLED yardstick**: HUG **+0.4956** on `diabetes`, CHUV **+1.4100** on `onset_type = unwitnessed`, Lugano **−1.1087** on `hyperlipidemia`, 15 of 15 rows defined at each. Under each centre's OWN variances the same three covariates are identified but read +0.5385, +2.0608 and −1.5275 — up to 46% apart, and the order of the three centres is the same either way | §5.5, §12.6, §13, DoD-14 | yes, run |
| **Lugano is 31 in the cohort and 30 weighted**, so §6.1's two columns differ on this workbook | §6.1, §12.6 | yes, run |
| The pooled row reproduces Stage 6 exactly: per-arm ESS `{0: 29.1990, 1: 30.3394}` — **written as a dict because an earlier draft wrote this pair bridging-first while §18b wrote it control-first, both unlabelled, and §12.12 pins them to four decimals against this table** — pooled ESS **59.4588**, max w **0.964955**, and the pooled `e` range [0.0269724, 0.9310700] is Stage 6 §11's. `TREATMENT_LABELS` is `{0: 'EVT alone', 1: 'bridging'}`, so every rendered per-arm column and every `CentreOverlap` dict is control-first | §6.1, §12.8 | yes, run |
| Dropping the `center` indicators inside a centre **changes no centre's worst** — 15 of 15 rows defined at each of the three against 19 of 19 with them, and the worst is `diabetes`, `onset_type = unwitnessed`, `hyperlipidemia` either way. Re-measured under §5.5's pooled yardstick and still inert | §8, §12.6 | yes, run |
| The unweighted SMD computed over the **93** rather than over `in_model`'s 92 differs by at most **0.0284**, on `onset_type = witnessed`, and no row changes side of the threshold | §4.4 | yes, run |
| `np.average` with weights summing to zero raises `ZeroDivisionError("Weights sum to zero, can't be normalized")`; it does not return `nan` | §3.1, §5.3, §12.4 | yes, run |
| `np.array([3.0]).var(ddof=1)` and `np.array([]).var(ddof=1)` are both `nan`, with a `RuntimeWarning` | §3.1, §5.3 | yes, run |
| `pd.Series([nan, 0.2]).abs().max()` is **0.2** and `np.abs(np.array([nan, 0.2])).max()` is **nan** — the same expression, two opposite wrong answers | §3.1, §3, §12.5 | yes, run |
| **`np.average` of the constant 14.0 gives 14.0 under one arm's weights and 14.000000000000002 under the other**, a difference of 1.776e-15; the same expression on a vector of `1.0`s or `0.0`s is exact | §3.1, §5.3a, §12.4a | yes, run |
| **What decides whether that difference appears is the WEIGHTS and not the arm size.** Σ[0.1, 0.9] and Σ[0.3, 0.7] both divide 14.0 exactly; Σ[0.2, 0.8] gives 14.000000000000002 and Σ[0.3, 0.7, 0.11] gives 13.999999999999998. Of 49 two-record weight pairs tried, 12 produce a difference. **§12.0's `constant_everywhere()` originally declared (0.1, 0.9) against (0.3, 0.7) — both exact — so it did not exercise §5.3a at all and §12.4a's rejected-form companion returned 0.0 on it.** Found at T2 by running the fixture; the control arm's weights are now (0.2, 0.8). An implementation report first attributed this to the arm sizes being 2 and 2 rather than 3 and 2, which is true of the case where the defect was found and false about the cause | §12.0, §12.4a | yes, run |
| An earlier form of §5.3's branch 4 compared the two weighted means and reported `nihss_baseline` — **14.0 on every record of both arms** on `cohort_frame()`'s cohort — as undefined. The arm-constant form returns 0.0, and the workbook's 19 rows are unchanged to every digit | §5.3a, §12.4a, DoD-5 | yes, run twice |
| `smd([1, 2, 3, 4], a=[0, 0, 1, 1], w=1)` is **2.8284271247461900** against `2√2 = 2.8284271247461903` — the roadmap's hand-computed criterion, to 4.4e-16 | §5.2, §12.3 | yes, run |
| `smd([0, 0, 1, 1], a=[0, 0, 1, 1], w=1)` is **nan** here and **0.0** under `pilots/analysis.py`'s form | §5.3, §12.4a, §16, §19 | yes, run |
| This stage's `smd` and the pilots' agree on **every row of both frames** — max difference 0.0 — so the divergence is confined to the constructed separation case and §12.4a must test it on a hand-built vector rather than on a frame | §5.3a, §12.4a, §16 | yes, run |
| `propensity.ess` on an empty arm and on an all-`nan` arm both give the **empty-arm** message, because `dropna()` empties the second; `ess` over the full 93-row `w` equals `ess` over the masked 92 | §6.3, §12.7 | yes, run |
| `cohort_frame()`'s cohort: 5 records, `in_model` 5, 19 rows, **0 undefined**, **14 rows reading exactly 0.0**, worst \|weighted\| **1.839089** on `tmax6_ml` and worst \|unweighted\| 2.2483. Every centre has one control, so **every within-centre SMD is undefined**, 0 of 15 defined at each of the two centres present, and the other two declared centres carry no cohort record | §12.0, §12.0.2, §12.7 | yes, run |
| §12.0.2's golden vector: `core_ml` 0.723747 → 0.635113, `tmax6_ml` 2.2483 → 1.839089, `center = HUG` 0.258199 → 0.139181, `center = CHUV` −0.258199 → −0.139181, `penumbra_ml` 1.630961 → 1.318924. **Every figure here is `_fmt`'s six significant figures, correctly ROUNDED** — `penumbra_ml`'s unweighted value is 1.6309612667875129 and an earlier draft wrote it 1.630960, truncated at the seventh digit, which is 1.27e-6 from the true value and therefore outside §12.0.2's own 1e-6 tolerance. Found at T3 by asserting the spec's literal | §12.0.2 | yes, run |
| `hand_frame()` still cannot reach this stage: Stage 5's P2 raises with "HUG contributes no EVT alone patient; CHUV contributes no EVT alone patient" | §12.0 | yes, run |
| Audit entries through the six landed stages: **22** (load 7 / derive 4 / classify 1 / build 6 / fit 4), confirming Stage 6 §9's ledger; Stage 7 adds 2 for **24** | §9, §12.9 | yes, run |
| The 19-row pooled table takes **7.4 ms** (100 repetitions in 0.738 s); the whole stage is about 30 ms | §2 | yes, run |
| pandas **2.3.3**, numpy **1.26.4** in the project environment — numpy at the version Stage 6 §18g's `reference` group resolved it to | §3.1, §18 | yes, read |
| **A misaligned `Propensity` makes B5's `df.loc[ps.in_model, C.TREATMENT]` raise `AssertionError('')` — bare, empty message — and an object-dtype mask carrying `<NA>` makes it raise `ValueError("Cannot mask with non-boolean array containing NA / NaN values")`. In both cases B1's and B2's assembled messages are discarded**, which is why §4.5 raises in two phases | §4.5, §4.5a, §12.11, DoD-8 | yes, run |
| Blanking `penumbra_ml` on one `in_model` record of the fixture cohort gives that row **`n = 4` against every other row's 5** and a weighted SMD of **1.037070**, with `in_model` unchanged at 5 — the per-row denominator, which no unmodified frame exercises | §4.4, §12.2a, DoD-8a | yes, run |
| `smd` returns `nan` for a `nan` anywhere in `x` (branch 3, via a non-finite pooled SD) and for a `nan` anywhere in `w` (`w[t].sum()` is `nan`, so branch 2 does not fire and `np.average` propagates). Neither returns a number | §5.3, §12.4 | yes, run |
| On pandas 2.3.3 `(string_series == level)` is nullable `boolean`; `.astype(float)` maps `<NA>` to `nan` rather than raising; and `df.loc[<boolean with NA>]` treats `<NA>` as **False** rather than raising. So §4.2's `_levels` and §8's `_centre` need no extra guard for a missing `center` | §4.2, §8, §12.2 | yes, run |
| `statsmodels` is a **declared project dependency** (`pyproject.toml`), not a `reference`-group one, so §12.12's MLE companion needs no `importorskip`; `firthlogist` is the group-gated one | §12.12 | yes, read |
| **B5's `df.loc[ps.in_model, C.TREATMENT]` raises `KeyError('ivt')` on a frame without the column**, discarding B3's message — so naming `TREATMENT` in B3 was necessary and not sufficient, and B3 had to move into phase 1 | §4.5, §4.5a, §12.11 | yes, run |
| **A `Propensity` whose `e` ALONE is misaligned reached `_centre` and raised pandas' `IndexingError` after `balance_smd` was recorded** — a log describing half a diagnostic, from an input B1 was meant to reject. B1 now checks all three Series and names which is misaligned; measured afterwards, every B1-B5 failure records **zero** entries | §4.5, §12.11 | yes, run |
| **A literal `|` inside a header cell breaks the rendered table**: `data._md_table` does no escaping and sizes the separator from `len(rows[0])`, so `balance_smd` rendered an **8-cell header against a 6-cell separator** and `overlap_by_centre` **15 against 13**. Byte-identity across hash seeds passed throughout. After renaming the two columns to `abs SMD < 0.1` and `worst abs SMD`, both render 6/6/6 and 13/13/13 | §7.4, §12.9, DoD-13 | yes, run |
| **`_levels`' `mask(isna())` is unreachable from `assess`**: every `C.CATEGORICAL` name is in `PS_COVARIATES`, so a factor-absent record is complete-cased out of `in_model` first. Deleting the mask leaves the tables byte-identical on the unmodified fixture and on fixtures with `center` or `onset_type` set to `pd.NA` | §4.2, §12.2 | yes, run |
| **A record whose `center` is `pd.NA` cannot be weighted**, so the declared centres' defined `weight_share` cells always sum to exactly 1.0: `in_model` drops 5 → 4 and the sum stays 1.0. An earlier §12.6 companion asserted the opposite and was deleted | §12.6 | yes, run |
| `summary.SumStat` and `plot.SumStat` take a **`weighted.var`** argument and the package carries the literal strings `"ASD weighted var"` and `"ASD unweighted var"` — read out of `~/R/library/PSweight/R/PSweight.rdb`'s lazyload database, R itself being unavailable in that session. So [§9]'s denominator is available from the oracle directly and §16b's "difference by declaration" is a default | §16b, §14, T6 | yes, read from the installed package |
| **`_smd_detail`'s undefined-rows clause ran into the sentence after it.** With any row undefined the entry read "…, penumbra_ml A row's `n` is its own denominator [§11]…" — the clause ended with `", ".join(undefined)` and no terminator. Reached only by the constructed frame §12.7 builds, because neither the workbook nor `cohort_frame()`'s cohort has an undefined row, so no rendering by §18b, §18c, §18d or either review round had ever executed the branch. Fence corrected at T5; §12.9a asserts every branch of both `detail` strings as a sentence | §7.4, §12.9a, §21.5a | yes, run |
| **The 18 `python` fences of this document all parse, and the 15 that form a module assemble and run**: `assess` on `cohort_frame()`'s cohort returns 19 covariate rows and 5 centre rows, appends exactly 2 entries taking the audit 20 → 22, leaves the frame equal cell for cell, and reproduces §12.0.2's golden vector to every digit printed. Re-run 2026-08-18 during the §20 review, against the landed Stages 1-6 | §18b, §20 | yes, run |

**Claims about the R oracle**, run 2026-08-18:

| Claim | Where used | Verified |
|---|---|---|
| `PSweight` **2.1.2** and `logistf` **1.26.1** are installed in `~/R/library` (81 packages) and both load when `R_LIBS_USER` names that directory | §16b | yes, run |
| **`test_reference_r.py`'s gate reports them absent, and the mechanism is NOT `.Renviron`.** There is no `~/.Renviron` on this machine; `R_LIBS_USER` is unset, so R falls back to its own compiled-in default `~/R/x86_64-pc-linux-gnu-library/4.5`, **which does not exist** and is absent from `.libPaths()`. Confirmed identical under `--vanilla` and under `--no-init-file`, so `--vanilla` suppresses nothing here. The oracle Stage 6 opened is skipping on the machine it was installed on | §13, §16b, T6 | yes, run |
| **A fallback to "R's own default user library" is therefore a no-op** — that default is what R already uses. The repair that works is a candidate search for a user library that EXISTS; `_r_has` returns True under `_r_environment` once `R_LIBS_USER` is set to `~/R/library`, with the parent shell's variable unset | §16b, §19, T6 | yes, run |
| `PSweight::SumStat`'s balance table is `Mean 0, Mean 1, Weighted SD 0, Weighted SD 1, SMD`, and its SMD is `\|Δmean_w\| / √((SD_w,0² + SD_w,1²)/2)` — the **weighted** pooled SD, and **absolute** | §16b, §14 | yes, run, and reconstructed from its own printed cells |
| On unit weights the two conventions **coincide**: PSweight's 0.420264 is reproduced by §5.1's formula to printed precision | §16b | yes, run |
| On a deliberately imperfect score the **numerators agree**: ours 0.168377 × 1.023623 and PSweight's 0.169416 × 1.017345 both give 0.172354, so the entire difference is the declared denominator | §16b, §12 coverage map | yes, run |
| Under `PSweight`'s own unpenalised `glm` score, every in-model covariate's overlap-weighted SMD is **0** to printed precision — the exact-balance property, from the authors of the method | §6.4, §16b | yes, run |
| `SumStat` requires `xname` alongside `ps.estimate`; it does not infer the covariates from the frame | T6 | yes, run |

**Claims checked by reading the committed repository:**

| Claim | Where used | Verified |
|---|---|---|
| `data.KINDS` is the declared **nine** and `data.py:144-149`'s comment states that `model` is the last amendment this pipeline needs, "Stages 7-13 all render under `model` or under an existing kind" | §7.1, §10, §12.9 | yes, read |
| `propensity.ess` is public, documented as Stage 7's, and raises `model.FitError` on an empty arm with a message naming the `in_model` mask | §6.3, §12.7 | yes, read |
| `Propensity` carries `e`, `w`, `in_model`, `ess`, `fit` and `dropped`, all a balance table needs, and Stage 6 §11 states that Stage 7 "does not refit, does not re-weight, and does not compute a second effective sample size" | §0.2, §3, §11 | yes, read |
| `model.design` drops the reference dummy by name **before** the constant-column rule, so neither `X.columns` nor `dropped` contains it — Stage 6 §7.5's premise, and §4.2's | §4.2, §12.2 | yes, read |
| `config.SMD_THRESHOLD` is 0.10 and is read by nothing before this stage | §5.4, §10 | yes, read |
| **The Stage 1 §7 raw-name scan enumerates the tree with `os.walk` (`tests/test_config.py:409`), so `balance.py` and `test_balance.py` are picked up with NO edit to `test_config.py`** — which is what makes DoD-3's "unedited except for `test_config.py`'s four new assertions" true on that axis rather than hopeful. Run against the assembled module: zero raw headers | §3, §10, §12.10, DoD-3 | yes, run and read |
| `pilots/analysis.py:194-208` is `smd` with the `sd > 0 else 0.0` fall-through; `:210-247` is `overlap_by_centre` ranging over `df["center"].unique()` with `round()` on its cells; `:249-265` is `balance_table` built on `design(df, covars)` and carrying `variance_ratio` | §15, §16 | yes, read |
| `pilots/test_estimators.py` asserts `0 < worst |SMD| < 0.05` under a Firth score on synthetic data, which Stage 6 §13 carried forward as the expectation | §6.4, §13, §16 | yes, read |
| `TODOS.md`'s shared row-builder item names "when Stage 7 or Stage 14 first wants the table under a non-`missingness` heading" as its trigger, and describes a four-way absence classification with per-centre absence counts | §10, §15 | yes, read |
| SAP [§9] specifies the unweighted pooled SD, the 0.10 threshold, the full confounder set and the within-centre report; [§6] names four vascular risk factors as negative controls and excludes `penumbra_ml` for collinearity | throughout | yes, read |

### 18b. This document's own code, executed — 2026-08-18

A specification that ships code is a specification whose code has to have been executed, which is
Stage 6 §18c's rule and the failure its §18b's first row found. **Every `python` fence in this
document except the signature listing of §3 and the `config.py` block of §4.3 was extracted
verbatim, assembled into one module, and run** against the workbook and against `cohort_frame()`,
with §4.3's block supplying `NEGATIVE_CONTROLS` and `BALANCE_SET` from the same text T1 will paste
into `config.py`.

| Block | Executed result |
|---|---|
| **All 18 fences parse and 16 assemble into a working module** | `ast.parse` on each, then one `exec`. The two excluded are §3's `: ...` signature listing and §4.3's `config.py` block, which is applied separately because it belongs to another file |
| `assess` end to end on the workbook | Returns a `Balance` of **19 covariate rows and 5 centre rows**, appends exactly **two** entries taking the audit from 22 to 24, and leaves the frame **equal cell for cell** to a deep copy taken beforehand (§0.2, §12.10) |
| `Balance.worst()` | `('diabetes', 0.379832)` with an **empty** undefined tuple, from one call (§3) |
| `Balance.unbalanced()` | `('center = HUG', 'center = Lugano', 'hypertension', 'hyperlipidemia', 'diabetes')` — §6.4's five, by name |
| `_smd_table` rendered | 20 lines, one header and 19 rows, the verdict column reading `yes`/`no` only, the threshold interpolated as `\|SMD\| < 0.1` from `C.SMD_THRESHOLD` |
| `_overlap_table` rendered | 6 lines: four declared centres, the pooled row, and a header of **13** columns. USZ's row reads `missing` in every overlap cell with the structural status; the pooled row reads `93 / 92`, arms `53 / 39`, ESS `29.199 / 30.3394`, max w `0.964955`, share `1` |
| `_range` on an empty arm | `missing`, once — not `missing–missing` |
| The two `detail` strings | Both build and both name their own counts: "19 declared level(s) of 14 [§6] covariate(s)", "5 row(s) reach \|SMD\| = 0.1", "1 centre(s) contribute no cohort record" |
| `Audit.to_markdown()` through all seven stages | Renders, 28,323 characters on the fixture, with `## Fitted models` carrying all six `model` entries |
| **`_status` on the fixture cohort returned `overlap reported`, not `no weighted contrast`** | **A correction to this document, found by running it.** Each fixture centre has an arm of **one**, not zero, so `ess` and the propensity range are computable and only the SMD is undefined. §6.3, §12.0 and §12.7 said all three statuses were exercised there; two of them are, and the third is reached by no available frame. The code is right and three passages were wrong; they now say so and §12.7 constructs the case |

**What this does not establish.** It ran the assembled fences against real frames, not `balance.py`
as a module with its own imports, its `Final` constants and the AST scans §12.10 asks for — that is
T1-T5's work and §12 is what checks it. What it does establish is that no block in this document is
of the kind Stage 6 §18b's first row found: called and undefined, or written and unrunnable.

### 18c. And again after the §20 revision — 2026-08-18

Re-run against the revised text, because a revision that changes three code fences is a revision whose
code has to be executed again. **22 fences now, all parse; the 16 that form the module assemble and
run**, and the two additions this document does not own — §12.0.3's fixture and §16b's
`_r_user_library` — are excluded with §3's signature listing, §4.3's config block and §12.0's test
imports and vectors.

| Claim | Executed result |
|---|---|
| `assess` on `cohort_frame()`'s cohort, revised | **19 covariate rows, `centres` of 4, `pooled` a field**; two entries, audit 20 → 22; frame equal cell for cell; §12.0.2's golden vector unchanged to every digit printed |
| The rendered overlap table | 6 lines, 13 columns, and **three distinct `status` values across all five rows** — `overlap reported` twice, `structural non-positivity [§3]` twice, and `overlap reported` on the pooled row. No fourth literal (§6.2) |
| `weight_share` | the four declared centres' **defined** cells sum to exactly 1.0, the two structural rows render `missing`, and `pooled` is exactly 1.0 (§12.6) |
| **B1 on a shifted index** | `SchemaError` carrying B1's message — `5 mask row(s) and 5 frame row(s), sharing 0 index label(s)`. The earlier single-collection form gave `AssertionError('')` on the same input |
| **B2 on an object mask with `<NA>`** | `SchemaError` carrying B2's message. The earlier form gave pandas' `ValueError` |
| The per-row denominator | `penumbra_ml` blanked on one `in_model` record: that row `n = 4`, SMD 1.037070; every other row `n = 5`; `in_model` still 5 (§12.2a) |
| `_NO_CONTRAST`, constructed | CHUV reads `no weighted contrast` with `ess` missing for the emptied arm and a `missing` worst, the centre still present, and `propensity.ess` on that same arm still raising `FitError` behind the guard (§6.3, §12.7) |
| `smd`, seven routes | 2.82842712474619; `nan` for separating, one-armed, arm-of-one, zero-weight arm, `nan`-in-`x` and `nan`-in-`w`; `0.0` for the same constant in both arms |
| **`_smd_detail` built a throwaway `Balance(tuple(rows), ())`** | **A defect the revision itself exposed.** §20's `pooled` field made it a `TypeError` on the first call — but while it worked it also wrote the threshold predicate a second time inline, beside `unbalanced()`'s. §3 now defines `_worst` and `_over_threshold` once and the dataclass delegates; measured green afterwards |
| **B1's message quoted only the two lengths** | Both are 5 on the shifted-index fixture, so the message read as though nothing were wrong in the one case B1 exists for. It now reports the shared label count, which is 0 there |

### 18d. And a third time, after the §21 revision — 2026-08-18

§21 changed six code fences — `smd` split into `_pooled_sd` / `_ratio` / `smd`, `CovariateBalance` gained
`sd`, `_table` and `_centre` gained `sds`, `_assert_balance_inputs` moved B3 into phase 1, and both table
headers were renamed — so the module was assembled and run a third time. **24 fences now, all parse; the
18 that form the module assemble and run.** The fence map was re-derived by *content* rather than by index
(the `: ...` signature listing, the `config.py` block, the two `test_balance.py` blocks and the
`test_reference_r.py` block are identified by their first line), because §21's additions shifted every
index and a hardcoded skip set is the one part of this harness that silently rots.

| Claim | Executed result |
|---|---|
| `assess` on `cohort_frame()`'s cohort | 19 covariate rows, `centres` of 4, `pooled` a field with status `overlap reported`; audit 20 → 22; frame equal cell for cell; **§12.0.2's golden vector unchanged to every digit** — the yardstick change does not move the pooled table, because its yardstick was already its own population's |
| `CovariateBalance.sd` | present and populated on every row; 14 of the fixture's 19 read 0.0, which is §5.3's legitimate-zero population being constant |
| Both rendered tables | `balance_smd` **6 / 6 / 6** markdown cells across header, separator and body; `overlap_by_centre` **13 / 13 / 13**; no cell of either contains a `\|`. Before §21 they were 8 / 6 / 6 and 15 / 13 / 13 |
| B1, parametrised over the three Series | `SchemaError` naming `in_model`, naming `e` alone, and naming `w` alone. The `e`-alone case previously reached `_centre` and raised `IndexingError` |
| B2, B3 | `SchemaError` naming themselves; B3 fires on a dropped `TREATMENT`, where it previously produced `KeyError('ivt')` from B5 |
| **Every precondition failure records nothing** | `len(audit.entries)` is 0 after each of the four raises above, so no failing `assess` leaves a partial log |
| The within-centre column | every centre row's `sd` equals the pooled row's for the same covariate; the fixture's within-centre worsts stay `missing` because `_ratio`'s two-record floor survived the split (§5.2), which is what keeps §12.7's arm-of-one assertions true |

**What this still does not establish**, unchanged from §18b: it ran the assembled fences, not `balance.py`
as a module with its own imports and the AST scans §12.10 asks for, and it has not read a `detail` string
as prose — which §21.5 names as where the next round starts.

## 19. What this spec changed elsewhere

**Five** entries, all landing with this document rather than with the implementation, because a
specification that contradicts the roadmap is worse than no specification. Precedent: Stages 3, 4, 5
and 6 landed their roadmap edits with the document.

| # | Change | Files | The question it settles |
|---|---|---|---|
| 1 | Roadmap Stage 7 gains its `**Spec:** specs/stage7_balance_and_overlap.md` line | `implementation_roadmap.md` | Stage 7 was the next stage without one, and the first since Stage 1 |
| 2 | Roadmap Stage 7's **Accept when** gains the **declared-level rule**: every declared level of every [§6] factor gets a balance row, including each factor's reference level and a level nobody has, because the design matrix carries neither | `implementation_roadmap.md`, §4.2 here | Whether "judged against the full [§6] confounder set" means the covariate *names* or the covariate *columns*. It means the names: measured, the design matrix omits three of the nineteen rows on this cohort, and one of the three carries the worst residual imbalance in it. The roadmap's own criterion could be met by a 12-row table |
| 3 | Roadmap Stage 7's **Accept when** gains the **zero-variance rule**: a covariate constant within each arm and different between them returns missing, not zero | `implementation_roadmap.md`, §5.3 here | Whether "returns missing rather than zero for a variable observed in only one arm" is the whole of the criterion. It is not: measured, there are four ways an SMD is undefined, the pilots' implementation already handles the one the roadmap names, and the one it gets wrong reports a perfectly separating covariate as perfectly balanced. The roadmap named the branch that was already right |
| 4 | Roadmap Stage 7's **Accept when** gains the **reconciliation**: the pooled row of the within-centre table is asserted against Stage 6's `overlap_weights` entry rather than recomputed, and no second effective sample size exists | `implementation_roadmap.md`, §6.1, §6.3 here | Where "a within-centre overlap table alongside the pooled one" puts the pooled one. Stage 6 already publishes those cells, so a second table would be a third rendering; making the pooled report a row of the same table, computed by the same code, is what stops the two disagreeing |
| 5 | Roadmap Stage 7's negative-control bullet is **corrected**: it reads "the four vascular risk factors reported as negative controls", and `BALANCE_ONLY` holds **five** names. The fifth, `penumbra_ml`, is excluded from the propensity model for collinearity and is not a negative control | `implementation_roadmap.md`, §4.3 here | Whether `BALANCE_ONLY` is the negative-control set. It is not, and an implementer reading the roadmap beside `config.py` would have made it one — which would report a covariate excluded for exact collinearity as a covariate weighting was expected not to fix. The `role` column is computed from `PS_COVARIATES_FULL` instead, which is [§6]'s own identity |

**One amendment to a Stage 6 file that is not a roadmap edit** (§16b, T6): `_r_environment` in
`test_reference_r.py` gains a **candidate search** for a user library that exists, when `R_LIBS_USER`
is unset. Stage 6 §16b.2's whole argument is that a silently-skipping oracle is worse than no oracle,
and measured on this machine the oracle it installed is silently skipping — not for the reason that
section anticipated, and not for the reason the first draft of §16b here gave either. Passing a
variable through does nothing when the variable is unset; **inheriting R's own default does nothing
when that default does not exist**, which is the case here and is why the fix searches.

**What this spec has had, and what it has not.** Every numeric claim in §4, §5, §6, §11, §12 and §18
was produced by running code against the workbook, against `cohort_frame()`, against hand-built
vectors and against `PSweight` — and one of them, §5.3a, is a defect in this document's own first
draft that only running it exposed: the estimator as first written reported a covariate identical in
both arms as undefined, on the committed fixture, and no amount of reading would have shown it.
Stage 6 §19 predicted exactly that yield for a second pass and was right about Stages 4 and 5. This
document has now had **two independent review rounds**, recorded in §20 and §21, which found seventeen
and twelve items respectively — and the two sets are **disjoint**, which is the whole argument for the
second round having been worth running. Stage 6 took four reviews to reach
four disjoint defect sets of seventeen each, and there is no reason to think this one is finished at
two — only that each round has been smaller than the last, which is what the extra care in §5 and §6
was for. **Round 1 told round 2 to start at §7.2 and §12.9, and that is where round 2's second P1 came
from** — a header cell holding a `|`, found by reading a cell back out of `to_markdown()`, which nothing
in §18b's two assemblies had done. §21.5 says where a third should start.

## 20. Review round 1 — 2026-08-18

An independent engineering review, run before any of Stage 7 was implemented. **Seventeen items: three
that changed behaviour, one asserted number, one environment fact, nine errata, and two the revision
found in itself.** Method: the 18 `python` fences of the pre-review text were extracted, parsed,
assembled into a module and run against `cohort_frame()`'s cohort
through the landed `derive → classify → build → fit`, and then against constructed frames aimed at the
branches §12 describes; the landed modules' claims were checked by reading `config.py`,
`propensity.py`, `model.py` and `data.py`; and the R environment was probed directly. Every item below
was **measured, not argued** — the two exceptions are marked.

**What the assembly confirmed, so it does not have to be re-derived.** 19 covariate rows, 5 rendered
overlap rows, 0 undefined, 14 rows exactly 0.0, §12.0.2's golden vector to every digit printed, the
frame returned equal cell for cell, exactly two audit entries appended, `smd` reproducing 2√2 and all
four §5.3 branches, `_NO_CONTRAST` reachable exactly as §12.7 says, and Stage 6's `FitError` intact
behind §6.3's guard. **The core of this stage was right and the review did not change it.**

### 20.1 The three that changed behaviour

| # | Finding | Where it landed |
|---|---|---|
| 1 | **`_assert_balance_inputs` collected B1 and B2 and then crashed before it could raise them**, because B3-B5 read through the mask B1 and B2 are about. Measured: a misaligned `Propensity` gave `AssertionError('')` — bare, no message — and an object-dtype mask with `<NA>` gave pandas' `ValueError`. B1's paragraph about a finite table for a population that does not exist was assembled into `bad` and discarded two checks later, and §12.11's first two tests could not have passed as written | §4.5 now raises in **two phases**; §4.5a records the measurement and why the phases must not be tidied back into one; §12.11 gains both companions and the cross-phase assertion; DoD-8 gains the second half |
| 2 | **`Balance.centres` carried the pooled row**, which cost three things: `len(centres)` was not `len(CENTER_ORDER)`, the pooled row needed a fourth status literal `"pooled"` that §6.2's three-status table never declared, and a caller counting reported centres off `centres` would have counted the whole cohort as one. Measured in the rendered table, whose last row read `pooled` in the `status` column | `Balance` gains a **`pooled` field**; `_POOLED` is deleted and `_status` is called for every row; §6.2 states there is no fourth; §3, §6.1, §9, §11, §12.6, §12.7 and §12.8 follow |
| 3 | **The per-row `n` — §4.4's [§11] denominator — had no test on any frame.** It reads 92 on every workbook row and 5 on every fixture row, so an implementation masking once over the frame instead of once per row would have been green. Measured reachable by a one-line edit: `penumbra_ml` blanked on one `in_model` record gives that row `n = 4` and SMD 1.037070 with `in_model` unchanged | new **§12.2a**, four assertions including the `PS_COVARIATE` contrast; §4.4 records the constructed-test reason; T3 and DoD-8a name the mask that has to be inside the loop |

### 20.2 The one that was an asserted number

**The design-matrix shortcut is 16 rows over `BALANCE_SET` and 11 over `PS_COVARIATES`, never 12.**
Measured. Twelve is `design_matrix`'s `n` — `X.shape[1] + 1`, the parameter count including the
intercept (`propensity.py:477`) — and it had propagated into §4.2, §12.1, §12.2, §18 and DoD-6, where
§12.2 *asserts* it and DoD-6 asks a human to reproduce it. The 19-row decomposition had the same shape
of error: written 12 + 5, which is 17, against the true 14 + 5. **§4.2's argument was unaffected** —
three of the seven indicators the two factors expand to are absent from `design`'s output, and one of
them carries this cohort's worst residual imbalance — so what changed is every count and none of the
reasoning. §12.2's companion is now built over `BALANCE_SET`, which isolates the encoding difference
from the covariate-list difference, and pins the `PS_COVARIATES` form beside it because that is what
`pilots/analysis.py:249-265` actually does.

### 20.3 The environment finding

**§16b's diagnosis of the skipping R oracle was wrong and its prescribed repair was a no-op.** The
draft blamed `Rscript --vanilla` for not reading `.Renviron`; measured, **there is no `~/.Renviron` on
this machine**, `R_LIBS_USER` is unset, R falls back to its own compiled-in default
`~/R/x86_64-pc-linux-gnu-library/4.5`, and **that directory does not exist** — the packages are in
`~/R/library`, which nothing points R at. Confirmed identical with and without `--vanilla`. So "fall
back to R's own default user library" changes nothing: R already defaults there. T6 now **searches**
candidate libraries and asserts that `_r_has` returns True under the environment `_r_environment`
builds with the parent's variable unset, which is the state the oracle was skipping in.

### 20.4 The nine errata

Each is a measured-wrong fact or an internal contradiction; none changed a decision.

| # | Where | What was wrong |
|---|---|---|
| E1 | §3 | `CentreOverlap.status`'s docstring named `"reported"` and `"structural non-positivity"`; §7.4 declares `"overlap reported"` and `"structural non-positivity [§3]"`. Stage 14 matches on these by value |
| E2 | §8 | "the same number **computed once**" — `_table` runs twice, once in `assess` and once inside `_centre(df, ps, None)`. The duplication is right and priced in §2; the claim was not |
| E3 | §3 | the import block omitted `from typing import Final` while §7.4 declares three `Final[str]` constants. Measured harmless at runtime under PEP 563, which is why it survived |
| E4 | DoD-12 | `grep "model\."` collides with the module docstring §17 *requires*. Replaced with an import-shaped grep; §12.10's AST scan was always the real check |
| E5 | §4.5 B3 | B3 checked `BALANCE_SET` only, and `C.TREATMENT` is in neither list while B5 and every arm mask read it — a bare `KeyError` from inside a collected assertion |
| E6 | §5.3 | two further routes into `nan` exist and were undocumented: a `nan` in `x` (reaches branch 3) and a `nan` in `w` (`not nan` is False, so `np.average` propagates). Both measured; both matter because `smd` is Stage 12's too |
| E7 | §12.6 | "`weight_share` over the reported centres sums to 1" is wrong twice: a structural centre's share is `nan`, not 0.0, and a `no weighted contrast` centre carries a share while not being reported |
| E8 | §1, §12.8 | the `workbook` fixture's shape was unspecified, and §12.8 needs `overlap_weights` and Stage 7's entries in **one** `Audit` — the opposite of `test_propensity.py`'s deliberately-fresh one. Now §12.0.3 |
| E9 | §12.0 | "from the modules that declare them" left the two modules unnamed: `hand_source` and `run` are `test_data.py`, `cohort_frame` is `test_cohort.py` |

### 20.4a Two the revision found in itself

Applying the round is a change to shipped code, so the fences were assembled and run again (§18c), and
that second run found two more:

- **`_smd_detail` built a throwaway `Balance(tuple(rows), ())`** to reach `worst()`, because the first
  audit entry is recorded before the centres exist (§8). The `pooled` field turned that into a
  `TypeError` on the first call — and while it worked, it also spelled the threshold predicate a
  **second** time, inline, beside `unbalanced()`'s. §3 now defines `_worst` and `_over_threshold` at
  module level and the dataclass delegates: one definition, two callers, no temporary object. The DRY
  problem predates the round and would have survived it.
- **B1's message quoted the two lengths and nothing else**, and on the shifted-index case both are 5 —
  so the message read as though nothing were wrong, in precisely the case B1 exists for. It now reports
  the count of shared index labels, which is 0 there.

### 20.5 What the round did not find, and what it could not

**Not found:** any error in §5's estimator, in §5.3a's arm-constant argument, in §6.3's guard, in
§6.4's attribution of the residual imbalance to Firth, in §4.4's population decision, or in §7.1's
claim that `model` needs no tenth kind. Each was checked — the first five by running, the last by
reading `data.py:153-158`.

**Could not:** the round did not run the R oracle's comparison (it established only that the gate can
be opened), did not read a cell back out of `Audit.to_markdown()`, and did not touch the workbook —
every measurement above is from `cohort_frame()`'s cohort, from constructed frames, from hand-built
vectors or from reading the landed modules. So §12.12's workbook numbers and §18's per-centre table
carry the same standing they had before this round: verified once, by the probe §18 describes, and
re-checked in code only when T7 runs. **Two claims in this document are still unverified by anything
but that probe** — the 0.0284 cost of §4.4's decision and the 67-94% imbalance reduction of §6.4 —
and both are recorded in §18 as `yes, run` rather than asserted in a test, which is the honest status
and not an oversight.

**Outside voice: none.** A cross-model pass was attempted and failed on an expired credential before
producing any finding. So this document has had one reviewer, not two, and §19's closing paragraph
says what a second should read first.

## 21. Review round 2 — 2026-08-18

A second independent review, by a reader with no context from round 1 beyond §20 as an exclusion list,
run before any of Stage 7 was implemented. **Twelve items, two of them P1, and the set is disjoint from
§20's seventeen.** Method: the fences were extracted and the module assembled *independently* — the fence
map was re-derived rather than inherited, and it agreed — then `assess` was run on the fixture cohort and
on the workbook, `Audit.to_markdown()` was rendered and **its cells read back out**, all five
preconditions were driven individually and in combination, byte-identity was checked across two hash
seeds in subprocesses, and `PSweight`'s source was read out of the installed package's lazyload database
because R would not start in that session.

**Round 1 predicted where this round would pay.** §19 pointed a second reviewer at §7.2 and §12.9,
because two fence assemblies had called `to_markdown()` and neither had read a cell out of it. That is
exactly where the second P1 came from. The prediction is recorded because it held, and because the same
reasoning names the next gap in §21.5.

### 21.1 The two P1s

| # | Finding | Where it landed |
|---|---|---|
| 1 | **§20's E5 was an incomplete repair, and the incompleteness was the same defect one check further along.** E5 added `TREATMENT` to B3's absent-column list; B5 still read `df.loc[ps.in_model, C.TREATMENT]` two lines later, so a frame without the column gave `KeyError('ivt')` and B3's message was discarded. §4.5a's own argument, inside phase 2, in the case E5 existed to fix | the phase boundary moved: phase 1 is now **B1, B2, B3** — *can this be read* — and phase 2 is **B4, B5** — *is the data judgeable*. §4.5a records why a column check belongs with the mask checks, and §12.11 gains the `KeyError` companion that would have caught it |
| 2 | **Both audit tables rendered as malformed markdown.** `f"\|SMD\| < …"` and `"worst \|SMD\|"` put a literal `\|` in a header cell and `data._md_table` does no escaping, so `balance_smd` rendered an 8-cell header against a 6-cell separator and `overlap_by_centre` 15 against 13 — in every log this stage writes. **Byte-identity across hash seeds passed throughout**: a table identically broken twice is still identical | the columns are `abs SMD < 0.1` and `worst abs SMD` (§7.4); §12.9 asserts no cell holds the character with a companion that pastes one in; DoD-13 makes a human check the rendered grid. `data.py` is untouched, so §10's no-amendment claim stands |

### 21.2 The one that was a decision, not a defect

**The within-centre `worst abs SMD` was standardised by each centre's own pooled SD**, so no two rows of
that column shared a yardstick, none shared one with `SMD after` in the entry above it, and CHUV's 2.06
against Lugano's 1.53 was partly an artefact of the smaller within-centre variance rather than the larger
imbalance. Measured: up to 46% apart. In a document that argues six times that one SD must serve two
columns, this was the same principle violated at ninety degrees — in the one table nobody had asked the
question about.

Resolved toward one yardstick everywhere (PI-facing decision, 2026-08-18): §5.5 is the rule, `_pooled_sd`
and `_ratio` are the split that makes it expressible, `CovariateBalance.sd` carries it, §8 harvests it
from the pooled table, and §12.6 and DoD-14 pin it. **What it did not change**: the covariate identified
as worst at each centre, and the order of the three centres. What it changed: the numbers are now
readable against each other. §13's warning about arm sizes is unaltered and remains the larger of the two
reasons not to rank the centres.

### 21.3 Two the round found in the spec's own tests

Both are tests §20 added or amended, and both were vacuous:

- **§12.6's `center = pd.NA` companion asserted an impossible premise.** `center` is a `PS_COVARIATES`
  member, so complete-casing removes such a record before Stage 7 sees it; measured, `in_model` drops
  from 5 to 4 and the defined shares still sum to exactly 1.0. Deleted, with the note that reaching the
  state needs a hand-built `Propensity` in §12.7's declared posture.
- **`_levels`' `mask(isna())` is unreachable from `assess`**, because every `CATEGORICAL` name is a PS
  covariate — so §12.2's bullet and T3's warning both described a failure no input can produce. Measured:
  deleting the mask leaves the tables byte-identical on three different fixtures. The mask stays and the
  test now asserts the **unreachability**, which fails if a factor is ever declared outside the
  propensity model.

### 21.4 The rest

| # | Finding | Where it landed |
|---|---|---|
| 3 | B1 checked `in_model` and `w` while its message claimed all three; a `Propensity` whose `e` alone was misaligned reached `_centre` and raised `IndexingError` **after `balance_smd` was recorded** | B1 covers all three and names which; §12.11 asserts a failing `assess` records nothing |
| 4 | T3's and T4's verification sets could not run at T3 and T4 — §12.5 needs `assess` and a `Balance` whose `pooled` field is a T4 deliverable; §12.6-12.8 need the audit | §17 reassigns every acceptance section that reads a `Balance` or an `Audit` to T5, and gives T3/T4 direct assertions on their own privates |
| 5 | `PSweight`'s `summary.SumStat` takes a **`weighted.var`** argument, so [§9]'s denominator is available from the oracle directly and §16b asserted a difference where there is a default | §16b: T6 runs it both ways, the `FALSE` pass being a direct check on §5's weighted column, and asserts PSweight's ddof rather than assuming it |
| 6 | §6.4's "between 67% and 94%" rounded its own lower bound **up**; the measured minimum is 66.5% on `center = CHUV` | 66%, with the full eight-row table in §18 and the bound now computed in §12.12 rather than quoted |
| 7 | `_smd_detail` hard-coded "the four vascular risk factors" against §4.3's computed-set promise and the very patch §12.1 applies | `len(C.NEGATIVE_CONTROLS)`, interpolated like every other count in that string |
| 8 | §18 wrote per-arm pairs bridging-first in one row and control-first in another, both unlabelled, while §12.12 pins them against §18 | every per-arm pair in §18 is now written `{0: …, 1: …}`, and the row states that `TREATMENT_LABELS` makes the rendered order control-first |
| 9 | DoD-12's grep used a path invalid after DoD-4's `cd`, so a missing-file exit 2 read as a pass, and its pattern missed `import propensity, model` | corrected on all three counts, with the reason for each recorded |

### 21.5 What round 2 could not do, and where a third should start

**It could not run R.** `Rscript` segfaulted inside its own startup in that session — `system("uname -a")`
returned status 139, so `utils` and `stats` failed to load — which is the PATH collision Stage 6 §16b.2
already documents, arriving again. So **every numeric claim in §16b remains unverified by both rounds**:
0.420264, the 0.168377 × 1.023623 = 0.169416 × 1.017345 = 0.172354 reconstruction, and PSweight's
exact-balance result under its own `glm` score. §16b has now been read twice and executed once, by the
§18 probe. **T6 is the only thing that closes it, and it is the one place a third round would add least**
— the work is running the oracle, not reading it again.

**A third round should start at §7.2's two `detail` strings.** Round 2 rendered both tables and read their
cells; neither round has read a `detail` string as a *reader* — checked that its counts match the table
beside it, that its prose is true of the data, and that its several conditional clauses read as English on
a workbook where the counts are 0, 1 and many. §12.9 asserts the tables and the entry `n`s; nothing
asserts that a sentence saying "5 row(s) reach abs SMD = 0.1" appears above a table with five `no` cells.
That is the same gap one level up, and it is the last part of this stage nothing has executed as prose.

### 21.5a — CLOSED at T5, and it found what §21.5 predicted it would

**The prediction held twice.** Round 1 pointed round 2 at §7.2 and §12.9, and that is where round 2's
second P1 came from; §21.5 then pointed the next reader at the two `detail` strings, and reading them
found a defect neither round could have seen.

**What it is.** `_smd_detail`'s undefined-rows clause ended with `", ".join(undefined)` and no
terminator, so whenever a row was undefined the log ran two sentences together: *"…, penumbra_ml A
row's `n` is its own denominator [§11]…"*. §7's fence now carries the full stop and says why.

**Why three assemblies and two review rounds walked past it.** §18b, §18c and §18d each called
`to_markdown()`, and round 2 read cells back out of the *tables* — but **the branch never rendered**,
because neither the workbook nor `cohort_frame()`'s cohort has an undefined row (§13). It required
the constructed frame §12.7 already builds for `_NO_CONTRAST`: complete-casing a centre's only
control leaves one control in the whole weighted set, so §5.3's first branch fires on all 19 rows and
the clause finally appears. **A conditional clause that no available frame reaches is a clause no
rendering has executed**, which is the `detail`-string form of §12's own point about `_NO_CONTRAST`
and the per-row `n`.

**Closed.** Both strings have now been read as prose against real rendered output under every branch
they have — 0 and many over-threshold rows, 0 and all-19 undefined rows, 0 and 1 thin centres — and
§12.9a asserts each branch's SENTENCE rather than only its count, so the terminator cannot be lost
again. **A fourth round should not start here.** What is left of §21.5 is its first half, and only
T6 closes that: run the oracle.

### 21.6 Applying a review is itself a step that can be incomplete

Recorded because it happened twice in one day, and the second time was worse than the first.

**§20's E5 named `TREATMENT` in B3 and left B5 reading it** — the repair reproduced the defect one check
along, and §21.1 found it. Then, applying §21, a sweep of the document against everything the round had
reported found **four more items that had been agreed and not written down**: §12.8's comparison having to
go through `_fmt`, because an `AuditEntry.table` holds strings and a float will not equal one; `_range`
rendering an arm of one as `0.257067–0.257067`, which is a decision and was sitting in a reviewer's
"probably fine, noted so it is a decision rather than an accident"; the raw-name scan's `os.walk`
enumeration, which is what makes DoD-3 true rather than hopeful; and the normative-numbers rule now at the
head of §18.

None of the four was a disagreement. Each was a thing everyone involved considered settled and nobody had
put in the file — which, for a document declaring itself **the sole source for the implementation**, is
indistinguishable from never having found it. An implementer reads §5 and §8, not a review transcript.

So the rule this document adds for its own maintenance, and it is the cheap one: **after applying a review
round, re-read the round against the file and check each item is in the file rather than in the argument
for it.** §20 and §21 exist to make that possible — they are the checklists, not the record — and §18's
precedence rule exists because a number agreed in a review still has to be written into all of its homes.
Both of the P1s in §21 would have been caught by an implementer running the code; neither would have been
caught by re-reading the prose, which is why §18b, §18c and §18d each end by naming what they did not run.
