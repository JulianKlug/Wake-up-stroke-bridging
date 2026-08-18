# Stage 6 spec — propensity score and weights

Implements roadmap Stage 6 [§7]. Section references in brackets are to `statistical_analysis_plan.md`.
Numbers and decisions referenced as DECISION *n* are established in Stage 0 and recorded in
`../out/stage0_data_inventory.md`. Constants named in `SMALL_CAPS` are Stage 1's and live in
`config.py`; `stage1_config_and_data_contract.md` is their specification. The frame this stage
receives is specified by `stage5_cohort_construction.md` §11.

**Status.** Written 2026-08-13 against the landed Stages 1-5 and the workbook as verified in §18;
**reviewed 2026-08-14 in three rounds** — an engineering review, the corrections it produced, and an
independent second reviewer with no prior context — a second independent reviewer after those corrections landed, and an audit of what that reviewer left —
for **52 findings in total**, all acted on. Two are retractions of earlier rounds' conclusions, both
found by re-running them.
**§20 is the ledger of what moved and why**, §18b-§18d the measurements. It is the
**sole source for the Stage 6 implementation**: everything the implementer needs is here, and anything
not here is not to be invented. Where this document previously said one thing and now says another, §20
says which and the earlier form is quoted rather than deleted — because a reader who remembers the old
sentence needs to know it was wrong, not to wonder whether they misread it.

**Goal.** A propensity score exists for every cohort patient the [§6] covariates are complete on; the
weights are the [§7] overlap weights and nothing else; the fit either converges or raises, with no
second estimator anywhere on the path; and the effective sample size, the weighted population and
[§7]'s statement that both are conditional on this specification are in the audit log rather than in
a reporting layer that may forget them.

**Not in scope.** Standardised mean differences and the within-centre overlap table (Stage 7 [§9]),
any outcome model or estimate (Stages 8 and 9), the bootstrap (Stage 10), and the [§14] populations
(Stages 12 and 13). Stage 6 **adds no column to the cohort frame and edits no value.** It returns a
result object; the frame it was given comes back unchanged, and §12.13 asserts it.

**This is the first stage that fits anything**, and that changes the character of the failure modes.
Stages 1-5 fail by mislabelling a patient or dropping a row; those failures are visible in a count.
Stage 6 fails by returning a number, and §3.1 records the two ways it returns a *plausible* one:
`nan` propagating silently out of a masked-integer column, and an unlisted factor level entering the
design as the reference category.

**One thing this spec settles that no earlier stage could.** Roadmap invariant 4 — "no post-time-zero
variable appears in any model's design matrix — assert against an explicit denylist, not by
inspection" — becomes enforceable here, because `model.design` is the one function every model's
covariates pass through. §4.5's D2 is that assertion, and Stages 8, 9 and 12 inherit it by calling
`design` rather than by remembering (§9).

**And one question that was put and answered before a line was written: why is the estimator not a
library?** It was investigated rather than assumed, and §16b is the record. In short: there is no
maintained Firth implementation for Python, and the one that exists cannot be a runtime dependency of
this project. What replaces the dependency is four oracles, not a comment.

---

## 0. Where Stage 6 sits

```
  cohort.build(df, audit)  →  93 rows × 33 columns                  [Stage 5 §11]
                       │
                       ▼
  ┌────────────────────────────────────────────────────────────────────────────┐
  │  STAGE 6a — model.py            reads config.py only                        │
  │                                                                             │
  │   design(df, covariates)        [§6] reference-coded design matrix          │
  │     ├─ _assert_design_inputs    D1…D4 → SchemaError                         │
  │     ├─ declared FACTOR_LEVELS → Categorical → reference-coded  (§4.2)       │
  │     ├─ constant columns dropped, AFTER dummying                (§4.3)       │
  │     └─ returns (X, dropped)                                                 │
  │                                                                             │
  │   complete_cases(df, covariates)   the [§11] mask, RETURNED not applied     │
  │                                                                             │
  │   firth(X, y)                   [§7] penalised logistic         (§5)        │
  │     ├─ _assert_fittable         F3 rank / F4 X / F6 y → FitError (§5.4)     │
  │     ├─ l(b) + ½log|I(b)|,  score X'(y − p + h(0.5 − p))                     │
  │     ├─ step-halving on the PENALISED likelihood                            │
  │     ├─ converge on Δpll  OR a flat modified score              (§5.3)      │
  │     └─ raise FitError — never a second estimator     [invariant 5]          │
  └────────────────────────────────────────────────────────────────────────────┘
        │                                        ▲                 ▲
        │                                        │                 │
        │                          Stage 9 [§8] m_a(X)   Stage 12 [§14a]
        │                          enter model.py DIRECTLY and never
        │                          through propensity.py            (§0.1, §9)
        ▼
  ┌────────────────────────────────────────────────────────────────────────────┐
  │  STAGE 6b — propensity.py       reads config.py, data.Audit, model          │
  │                                                                             │
  │   fit(df, audit) -> Propensity                                              │
  │     ├─ complete_cases(df, PS_COVARIATES)      92 of 93 on v7   (§4.4)      │
  │     ├─ design(df[in_model], PS_COVARIATES)    11 columns       (§4)        │
  │     ├─ firth(X, y)                            7 iterations     (§5)        │
  │     ├─ _assert_probabilities(e)               F5 → FitError    (§5.5)      │
  │     ├─ weights, per-arm ESS                                    (§6)        │
  │     └─ four `model` audit entries                              (§7)        │
  └────────────────────────────────────────────────────────────────────────────┘
                       │
                       ▼
   Propensity(e, w, in_model, …)  →  Stages 7-11
   e and w are Series over ALL 93 cohort rows, nan where in_model is False (§8)
```

`fit` appends four entries of the new `model` kind to the `Audit` that `load()` created. It writes no
file, as no stage before it does.

### 0.1 Why two modules

**`design` and `firth` are not the propensity model's; they are the pipeline's numerics.** Stage 9's
outcome regression `m_a(X)` is a Firth logistic fit on a design matrix over
`outcome_model_covariates(outcome)` — the same two functions, a different covariate list and a
different response. Stage 12's proportional-odds model is a design matrix over
`STANDARDISATION_COVARIATES`. If both lived in `propensity.py`, Stage 9 would import a module named
for the exposure in order to fit an **outcome** model, and the first reader to tidy that import would
be right to.

So `model.py` holds what is outcome-agnostic and `propensity.py` holds the [§7] specification. Three
consequences, all of them the point:

- **`model.py` imports `config` and nothing else** — not `data`, so it records no audit entry, and
  not `pandas` for anything but a frame in and a frame out. It is the first module in the pipeline
  with no `Audit` parameter, because it is the first that is a *function* rather than a *step*.
  Its caller logs; it computes.
- **`propensity.py` is the only module that knows what the exposure is.** It reads `TREATMENT` and
  `PS_COVARIATES`; `model.py` reads neither, and §12.12's scan asserts it.
- **The repository's unit is one module, one spec, one test file**, which Stages 1-5 all hold. Stage 6
  is two modules and two test files under one spec, and that is a deliberate departure: the roadmap
  stage is one stage, the [§7] decisions are one set, and splitting the *document* would put the
  Firth argument in one file and the estimand it serves in another.

### 0.2 What Stage 6 does not touch

`fit` takes the cohort frame and gives it back untouched: **no column is added, no value is edited,
no row is removed.** This is the first stage of which all three are true — Stage 3 added columns,
Stage 4 added one, Stage 5 removed rows.

It matters because of what a column would cost. A `propensity` column on the frame would be a second
place the score lives, and Stage 10 refits in every one of `N_BOOT` replicates: a frame carrying a
column from the point fit, resampled into a replicate, is a replicate silently weighted by the wrong
score. `Propensity` is a return value, so a caller that wants the point score and the replicate score
at once has to hold two objects and name them.

**The index is preserved.** `e`, `w` and `in_model` are Series on the cohort's index, which is the
unrestricted frame's index (Stage 5 §0.2). A record is traceable from a weight back to the [§14]
population without a join on `case_id`.

## 1. Deliverables

| Path | Contents |
|---|---|
| `extended_bridging/model.py` | `design`, `complete_cases`, `firth`, `Fit` (eight fields — §3.2), `FitError` |
| `extended_bridging/test_model.py` | the acceptance tests of §12.1-12.7 — **including §12.4a's safeguards and its one unreachability assertion** — and §12.12; the synthetic fixtures of §12.0 **written out there in full**, both the three answer-known ones and the hard-fit family; and **three** of §16b's four oracles — `firthlogist`, the `scipy` objective check and the large-n MLE. The fourth, §12.0.2's golden vector, is `test_propensity.py`'s, because it pins `fit`'s output rather than `firth`'s |
| `extended_bridging/propensity.py` | `fit`, `Propensity`, `ess` |
| `extended_bridging/test_propensity.py` | the acceptance tests of §12.8-12.14, its **own** module-scoped `workbook` fixture, and §16b's second oracle — the golden vector of §12.0.2 |
| `extended_bridging/test_reference_r.py` | the gated R oracle of §16b.2 — the **only** Python in the repository that knows R exists, imported by nothing |
| `extended_bridging/tests/reference/firth_logistf.R`, `ato_psweight.R` | ~20 committed lines each: read a CSV, fit, write coefficients. No analysis logic and no covariate list (§16b.2) |
| `extended_bridging/config.py` | **amended** — the seven declared `FIRTH_*` numerics of §5.6, and nothing else (§10) |
| `extended_bridging/test_config.py` | **amended** — one assertion that `REFERENCE_LEVELS[f] == FACTOR_LEVELS[f][0]`, which is the coupling §4.2's by-name drop rests on and which nothing asserted (§4.2, §10) |
| `extended_bridging/data.py` | **amended** — `KINDS` and `_HEADINGS` gain `model`; `_MUST_NAME_CASES` deliberately unchanged, for Stage 5 §6.3's reason (§7.1, §10) |
| `extended_bridging/pyproject.toml` | **amended** — `numpy` becomes a direct import of a shipped module; `scipy` becomes a **test-only** import (§2), and §12.7's scan forbids it in `model.py` and `propensity.py`; a new **dev-only** `reference` group holds `firthlogist` behind its scikit-learn ceiling (§2, §16b) |
| `extended_bridging/test_data.py`, `test_eligibility.py`, `test_cohort.py` | **amended** — the three literal pins of `KINDS` and `_HEADINGS.values()`, and the two test names carrying the word *eight* (§10) |
| `extended_bridging/implementation_roadmap.md` | **amended, and landing with this document** — Stage 6 gains its `**Spec:**` line and three **Accept when** items (§19) |

`out/logs/audit_<label>.md` is an **output, not a deliverable**, as in Stages 2-5: gitignored, written
only when a caller asks. Stage 6 makes that distinction carry more weight than before — see §7.3 on
where the fitted coefficients live.

**Nothing under `specs/` may quote a case identifier, and nothing here does.** Every patient below is
a `HAND-N` or `COHORT-N` fixture record, or a count.

## 2. Environment

`uv`, Python 3.12, flat module layout, `import config as C`. Commands run from `extended_bridging/`.

```bash
cd extended_bridging && uv sync && uv run pytest -v
```

**Three dependency facts, and the third is a decision.**

- **`numpy` becomes a direct import of a shipped module** for the first time. It is already declared
  (`numpy>=1.26,<2.1`) and already installed as a pandas dependency; the pin's upper bound exists
  because NEP 50 scalar promotion could perturb bootstrap percentiles, and that reason now applies to
  this stage's coefficients as well as to Stage 10's intervals. Verified at 2.0.2 (§18).
- **`scipy` becomes a test-only import** — §16b's third oracle. Already declared (`scipy>=1.13`).
- **`firthlogist` is added to a dev-only `reference` group and is never imported by a shipped
  module.** §16b is the argument. The group is not part of `uv sync`'s default set; the one test that
  uses it is guarded by `pytest.importorskip`, so a plain checkout is green without it. This is the
  one place in the repository where a dependency exists to disagree with the code rather than to run
  it.

**`statsmodels` keeps the scope its own declaration gives it** — "retained for unpenalised
cross-checks in tests, not for any reported estimate". §12.4 and §12.5 are those cross-checks.

**What this costs to run.** One 92 × 12 fit, converging in 7 iterations (§18): seven 12 × 12
inversions and seven hat-diagonal `einsum`s over 92 rows, plus **one SVD of the 92 × 12 design** for
F3's `matrix_rank`, once per fit before the loop. The whole stage is some tens of microseconds and
there is nothing to cache. §15 says so as a standing instruction, for Stage 5 §15's reason: a cache is
a second source of truth, and here it would be a second propensity score.

**And what it costs Stage 10, since that is the only place the number could matter.** [§10] refits both
nuisance models in each of `N_BOOT = 2000` replicates, so the totals are roughly 28,000 12 × 12
inversions, 28,000 `einsum`s over ~92 rows, and 4,000 rank SVDs of a 92 × 12 — all of it on matrices
small enough that the per-call Python overhead dominates the arithmetic. **Plus §9's local-maximum
check**, which multiplies the fit count by however many probes Stage 10 draws: measured, 2000 replicates
× (1 primary + 62 probe fits) of the 92 × 12 design took **290 seconds** single-threaded. That is the
real cost of the only unguarded failure mode this stage hands on (§5.2c), and it is minutes rather than
hours — which is why §9 asks for the check rather than for a sampling of it. Nothing here is a performance
consideration at any point in this pipeline, and that is worth stating once so no later stage reaches
for a cache, a `lru_cache` on `design`, or a vectorised-over-replicates rewrite of `firth` on
performance grounds. If Stage 10 is ever slow, the cause will be the 2000-fold Python loop and not
this stage's numerics.

## 3. Module shape

```python
# model.py
from __future__ import annotations

from dataclasses import dataclass
from collections.abc import Sequence

import numpy as np
import pandas as pd

import config as C
```

```python
# propensity.py
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd

import config as C
import model
from data import Audit
```

**Two new types, and both are justified where Stages 3-5 declined one.** Stage 3 §3's argument was
that everything those stages produced was a frame, a column or an audit entry, and all three already
had types. That stops being true here: a fit is a coefficient vector *and* a fitted probability vector
*and* an iteration count, and a propensity result is four aligned objects plus two scalars. The
alternative is a tuple of six, positionally unpacked at three call sites.

```python
@dataclass(frozen=True)
class Fit:
    """One Firth fit. `p` is aligned to the rows of the X it was given, positionally."""

    beta: np.ndarray            # (k + 1,), intercept first
    p: np.ndarray               # (n,), fitted probabilities
    iterations: int
    converged_on: str           # "likelihood" | "score"  — §5.3
    columns: tuple[str, ...]    # the design's column names, in beta's order after the intercept
    first_step_norm: float      # ‖step‖ at iteration 1, BEFORE any rescale   ┐
    rescales: int               # steps shortened by the trust region [§5.2a] ├─ §3.2
    halvings: int               # total step-halvings across all iterations   ┘
```

```python
@dataclass(frozen=True)
class Propensity:
    """The [§7] propensity score and overlap weights, over the frame `fit` was given."""

    e: pd.Series                # float64, the cohort's index, np.nan where not in_model (§8)
    w: pd.Series                # float64, same index and same missingness
    in_model: pd.Series         # boolean, TOTAL — never missing; this is where the intent lives
    ess: dict[int, float]       # keyed by arm code, as TREATMENT_LABELS is
    fit: Fit
    dropped: tuple[str, ...]    # design columns dropped as constant, in design order
```

Public surface, and it is small on purpose:

```python
# model.py
def design(df: pd.DataFrame, covariates: Sequence[str]) -> tuple[pd.DataFrame, tuple[str, ...]]: ...
def complete_cases(df: pd.DataFrame, covariates: Sequence[str]) -> pd.Series: ...
def firth(X: pd.DataFrame, y: np.ndarray) -> Fit: ...

# propensity.py
def fit(df: pd.DataFrame, audit: Audit) -> Propensity: ...
def ess(w: pd.Series | np.ndarray) -> float: ...
```

**Every `python` fence in this document is valid Python**, including the signature listings above — hence
the `: ...` bodies. It is a small rule with a specific purpose: this
document is the sole source, an implementer copies out of it, and a fence that does not parse is a fence
whose contents nobody has checked. §18c is the check; three fences failed it before this rule existed.
Where a fence is an *excerpt* rather than a definition, it is still parseable and its elisions are
comments (§7.1) rather than `...` in a position that breaks the grammar.

Privates are `model._assert_design_inputs`, `model._assert_fittable`, `model._penalised_loglik`,
`model._probabilities`, and `propensity._assert_fit_inputs`, `propensity._assert_probabilities`,
`_record_exclusion`, and §7's detail/table helpers — one `_*_detail` per entry, plus `_completeness_table`,
`_design_table` (§7.5), `_fit_table`, `_weights_table` and `_weighted_population_table` (§7.6), which is
five tables across four entries because `overlap_weights` carries two. **There is no `_step`** — the Newton
step is computed inline in §5.2's loop, because factoring it out would put the trust region on one
side of a call boundary and the step-halving that consumes it on the other.

**`ess` is public and lives in `propensity.py`, not in `model.py`.** It is three lines and Stage 7's
within-centre overlap table needs it per centre per arm [§9]; a second implementation there would be
a second definition of the effective sample size. It is not in `model.py` because a Kish sum has
nothing to do with fitting anything, and `model.py`'s no-exposure-no-weights rule (§0.1) is asserted
by §12.12.

**Neither module is exempt from the Stage 1 §7 raw-name scan and neither may become exempt.** Neither
names a raw header.

### 3.2 Why `Fit` carries three numbers about the journey, not just the answer

`first_step_norm`, `rescales` and `halvings` describe how the fit was reached rather than what it is,
which is a different kind of thing from `beta` and `p` and needs its own justification.

**The precedent is `converged_on`**, which is already a field nothing but a test and a future Stage 10
report reads (§5.3). These three extend it for the same reason: §5.2's trust region and step-halving are
the safeguards that fire when a bootstrap replicate is hard, [§10] refits `N_BOOT` times, and a
safeguard whose activation nobody can count is a safeguard nobody can evaluate. §13 asks Stage 10 to
aggregate exactly this across replicates.

**And without them §12.4a cannot be written.** The four safeguards it tests are properties of the
*iteration*, not of the result: whether the first step exceeded the trust radius, whether the accepted
step was shortened, by how much. Without them the alternative is to test the safeguards by side effect
(does it converge at all?), which passes a broken trust region that still happens to converge.

**They are counters, not state.** The loop computes all three anyway; the fields cost three assignments
and make the frozen dataclass a complete record of one fit. Nothing outside a test and Stage 10 reads
them, and no downstream estimate depends on them — so §0.2's argument about a second place the score
lives does not apply: these are not the score.

### 3.1 The numerical facts this stage turns on

Declared once here, verified by running them on pandas 2.3.3 / numpy 2.0.2 as pinned by `uv.lock`
(§18). The first two are how this stage returns a plausible wrong number, and they are why §4.4's
mask and §4.5's D3 exist.

```
  Series([1, pd.NA], dtype="Int64").to_numpy(dtype=float)
      →  array([1., nan])      SILENTLY. No raise, no warning.

  firth(X, y) with one nan anywhere in X
      →  beta is all-nan, p is all-nan, and the loop TERMINATES NORMALLY.
         Two RuntimeWarnings are emitted (logaddexp, slogdet) and warnings are not errors.

  firth(X, y) with one nan anywhere in y
      →  every candidate step is nan, so no step is ever accepted and the fit raises
         FitError("step-halving exhausted …") at iteration 1. LOUD — but the message
         names the estimator and [§7]'s no-second-estimator rule, so it reads as a sparse
         replicate rather than as a missing response. F6 (§5.4) is why the message is right.

  pd.Series(pd.NA, index=idx, dtype="float64")
      →  TypeError: float() argument must be a string or a real number, not 'NAType'.
         RAISES. np.nan is the sentinel a float64 Series takes; §8 uses it.

  pd.Categorical(s, categories=FACTOR_LEVELS[c])  with a value not in categories
      →  that element becomes NaN, and get_dummies then gives it an ALL-ZERO row —
         which is exactly the encoding of the REFERENCE level. A record at an
         unknown centre is silently modelled as though it were HUG.

  np.linalg.slogdet(singular)     →  (0.0, -inf).   No raise.
  np.linalg.inv(exactly singular) →  LinAlgError.   Raises.
  1 / (1 + np.exp(-500.0))        →  exactly 1.0.   The upper clip is reachable in float.
  1 / (1 + np.exp( 500.0))        →  7.12e-218.     The lower clip is not.
  w.sum() ** 2 / np.sum(w ** 2)   with all-zero w  →  nan, from 0/0.
  np.sum(np.asarray([]) ** 2)     →  0.0.           An EMPTY arm reaches the same branch
                                                    as an all-zero one (§6.2).
```

**The two `nan`-into-`firth` facts compose into the failure this stage is built against.** A missing
covariate is not a `KeyError` and not a `nan` guard trip: it becomes a `nan` in the design, and the fit
runs to completion and returns `nan` coefficients, `nan` probabilities and a plausible iteration count.
Every weighted estimate downstream is then `nan`, and the first place it surfaces is a confidence
interval of `[nan, nan]` in Stage 10 with nothing pointing back here. That is why §4.4's complete-case
mask is **returned and asserted** rather than applied by `dropna` inside the fit, and why F4 checks the
design for `nan` immediately before `firth` sees it even though `complete_cases` has already run. A
`nan` in the *response* does raise — but with the wrong diagnosis, which is F6's job to fix.

**The `Categorical` fact is quieter and worse.** It does not produce a `nan`; it produces a *number*,
computed as though the patient were at the reference centre. D3 is the only thing between that and an
estimate.

**`slogdet`'s `-inf` is listed because `_penalised_loglik` branches on it, and it is not reachable from
inside `firth`.** With `FIRTH_WEIGHT_FLOOR > 0` the working weight is floored, so `X'WX ≥ FIRTH_WEIGHT_FLOOR
· X'X`, which is nonsingular whenever `X` is full rank — and F3 has already established that it is.
Measured (§18): a design whose columns are scaled small enough for the determinant to underflow loses
rank under `matrix_rank` first, so F3 fires before `_penalised_loglik` is ever called. The branch is
kept, and §5.2's docstring records what it is for now that it is known not to fire.

**The `exp(-500)` fact is why F5 is an assertion and not a clip.** `pilots/analysis.py:82` clips `e` to
`[1e-8, 1 - 1e-8]`; §5.5 declines that. A fitted probability of exactly 1.0 means the linear predictor
overflowed, which under a penalised likelihood means the fit is degenerate — and clipping converts a
degenerate fit into a finite weight of 0.0 that every downstream sum accepts.

## 4. The design matrix [§6]

### 4.1 The rule

```
  Reference-coded dummies for the declared factors, linear terms for everything else,
  constant columns dropped [roadmap Stage 6].

  Factors are CATEGORICAL = ("onset_type", "center").
  Their level sets are FACTOR_LEVELS and their baselines are REFERENCE_LEVELS.
  No interaction, no spline, no transform — [§6] says "linear terms only".
```

`design` returns `(X, dropped)`: the matrix, and the column names it removed as constant, in design
order. It does **not** return the complete-case mask — that is `complete_cases`, a separate function,
for §4.4's reason.

### 4.2 Declared levels, and never `drop_first`

```python
def design(df: pd.DataFrame, covariates: Sequence[str]) -> tuple[pd.DataFrame, tuple[str, ...]]:
    """The [§6] design matrix over `covariates`, reference-coded on REFERENCE_LEVELS.

    No intercept column: `firth` prepends its own, so a design matrix is never ambiguous about
    whether it has one. Constant columns are dropped and named in the return value rather than
    logged here, because model.py takes no Audit [Stage 6 §0.1] and its caller is the step.
    """
    _assert_design_inputs(df, covariates)
    factors = [c for c in covariates if c in C.CATEGORICAL]
    sub = df[list(covariates)].copy()
    for c in factors:
        sub[c] = pd.Categorical(sub[c], categories=C.FACTOR_LEVELS[c])
    X = pd.get_dummies(sub, columns=factors, dtype=float)
    X = X.drop(columns=[f"{c}_{C.REFERENCE_LEVELS[c]}" for c in factors]).astype(float)
    dropped = tuple(c for c in X.columns if X[c].nunique(dropna=False) <= 1)
    return X.drop(columns=list(dropped)), dropped
```

Four properties, each load-bearing, and each one a thing `pilots/analysis.py:25-41` does differently.

- **The level set is declared, not observed.** `pd.Categorical(sub[c], categories=FACTOR_LEVELS[c])`
  makes the column set a function of `config.py` and not of the rows. Measured (§18): on the cohort
  this produces twelve columns of which `center_USZ` is all-zero, because restriction 1 removed USZ.
  Without the conversion, `get_dummies` on the bare `string` column emits
  `center_CHUV, center_HUG, center_Lugano` and **no `center_USZ` at all** — so the all-zero column
  never exists, the constant-column rule has nothing to do, and the width is a property of the
  sample.
- **The reference column is dropped by name, never by position — and the reason is not the one it looks
  like.** It would be natural to argue that `drop_first=True` picks whichever level sorts first among
  those *present*, so a replicate with no HUG patient would silently rebaseline on CHUV. **Measured
  (§18): that is false here, and it is false because of the bullet above.** Once the column is a
  `Categorical` with declared categories, `get_dummies` emits *every* declared level — including an
  all-zero `center_HUG` when no HUG patient exists — and `drop_first=True` drops the first **declared**
  category, which is `center_HUG`. So on a declared `Categorical` the two forms agree, in every
  replicate, and `f"center_{REFERENCE_LEVELS['center']}"` **does not raise `KeyError`**: the column is
  always there. The rebaselining hazard is real but it belongs to the *bare* column — measured, on a
  bare `string` centre with no HUG record, `drop_first=True` drops `center_CHUV` and leaves only
  `center_Lugano` — which makes it the same defect as the missing declaration, not a second one.

  What by-name actually buys is independence from a coupling nothing else states: `drop_first` is
  correct only while `REFERENCE_LEVELS[c] == FACTOR_LEVELS[c][0]`. That holds today for both factors —
  `HUG` leads `CENTER_ORDER`, `witnessed` leads its levels — and `test_config.py:592` asserts only that
  the reference is *in* the level set, not that it is first. So `drop_first` would be silently
  equivalent now and silently wrong the first time someone reorders `CENTER_ORDER` or names a
  non-leading reference. §10 adds that assertion to `test_config.py`, and with it in place the by-name
  form is provably the same model as `drop_first` rather than accidentally the same one.

  **What a replicate missing the reference level does instead: F3.** The all-zero `center_HUG` is
  dropped as the reference, the surviving centre dummies then sum to 1 on every row and are collinear
  with the intercept, and `_assert_fittable` raises F3 with the rank and the column names (§5.4).
  Measured on a frame with no HUG patient: rank *k*−1 of *k*, one dependency. Loud, and Stage 10 counts it
  as a failed replicate exactly as it would have counted a `KeyError`. `pilots/analysis.py:40` uses
  `drop_first=True` on a column it never declared, which is the combination that is wrong; §16 declines
  the parameter and §12.2 measures what it does on both a declared and a bare column.
- **`nunique(dropna=False)`, not `nunique()`.** A column that is `nan` on every row has
  `nunique() == 0` and `nunique(dropna=False) == 1`; both are ≤ 1 so both drop it, but a column that
  is 1.0 on one row and `nan` on the rest has `nunique() == 1` and would be dropped as constant while
  carrying a `nan` into nothing. The `dropna=False` form is the one whose predicate matches its name.
  On this pipeline the point is moot, because `complete_cases` runs first; it is written this way so
  that it stays moot if a caller ever forgets.
- **The drop happens after dummying**, so both an all-zero dummy (an absent level) and an all-one
  dummy (a single-level factor in a subset) are removed. A factor with one level present contributes
  no column at all, which is correct: it is collinear with the intercept.

**`remove_unused_categories` is not called and must not be.** `pilots/analysis.py:39` needs it because
`pilot_data.py:99` made `center` a `Categorical`; Stage 2 §6 chose `string` precisely so that no dead
level survives a restriction, and this function creates its `Categorical` fresh from the declared
levels on every call. Calling it here would delete the all-zero `center_USZ` column *before* the
constant-column rule could name it, and `dropped` — which §7.2's table publishes — would be silently
empty.

### 4.3 What the constant-column rule is actually for

It is not tidiness. Three distinct things reach it, and only the first is visible on v7:

1. **A declared level absent from the cohort.** `center_USZ`, always, by [§3] restriction 1. Measured:
   one column, every time.
2. **A level emptied by a bootstrap replicate.** Stage 10 resamples 93 patients with replacement
   stratified by centre, so a replicate cannot lose a whole centre — but Stage 11's subgroup analyses
   and [§13]'s single-centre sensitivity row can, and a replicate *within* a subgroup can lose an
   `onset_type` level. The dummy is then all-zero and the design is exactly singular.
3. **A covariate constant in the cohort by nature.** `constant_covariates` (Stage 3) detects and logs
   these; measured empty on the cohort for both `PS_COVARIATES` and `PS_COVARIATES_FULL` (§18). Stage
   3 §7.1 assigned the *deletion* to this stage, at the point of use, and this is that point. Two
   stages cannot both own it.

**Nothing here mutates a covariate list.** `PS_COVARIATES` is read; the dropping happens to the
matrix, and `dropped` reports it. Stage 3 §7.2's rule is unchanged.

### 4.4 The complete-case mask is returned, never applied

```python
def complete_cases(df: pd.DataFrame, covariates: Sequence[str]) -> pd.Series:
    """Boolean, total, on `df`'s index: True where every covariate in `covariates` is present.

    [§11] is complete-case *per estimate*, with the denominator reported for each. So this is
    returned for a caller to report and to log, never applied here: `design(df.dropna(...))` would
    make the exclusion a property of a call nobody can see, and the excluded patient would leave
    the analysis with every table still reconciling [Stage 5 §5.1 P4's failure, one stage on].
    """
    return df[list(covariates)].notna().all(axis=1)
```

**Measured on v7: 92 of 93.** One control at Lugano is missing `core_ml` and `tmax6_ml` (and
`penumbra_ml`, and therefore `core_above_median`); the other seven [§6] covariates are complete on
every cohort record, and so are all four of `PS_COVARIATES_FULL`'s additions (§18).

Four consequences, and they are the whole of decision 2 of this stage:

- **The record stays in the cohort and loses its weight.** `e` and `w` are `nan` for it, `in_model`
  is `False`, and every weighted sum in Stages 7-11 ranges over `in_model`. It is not deleted from the
  frame, because Stage 5's cohort is the [§3] population and a covariate that was never recorded is
  not a [§3] restriction.
- **Absent, not 0.0.** A weight of zero is a patient who was weighed and found irrelevant; an absent
  weight is a patient who was never weighed. The two are the same in every sum and different in every
  denominator, and [§11] requires the denominator. This is Stage 3 §5.2's rule — reimpose the
  missingness explicitly — applied to a number rather than to a dichotomy. **The sentinel is `np.nan`
  and the intent lives in `in_model`**, not in the sentinel: `float64` has exactly one missing value and
  `pd.NA` is not it (§3.1, §8), so there is no dtype in which the *deliberateness* of this absence can be
  expressed by the value alone. That is why `in_model` is returned rather than left to `notna()`.
- **`in_model` is total and never missing**, for the reason Stage 4 §6.2 gives about the eligibility
  column: a three-valued mask is a mask that resolves by branch order at the first `if`.
- **The ATO population is the 92, and that has to be said out loud.** [§7] describes the target
  population as `h(X) = e(X){1 − e(X)}`, which is defined only where `e` is. So the estimand's
  population is the complete-case set, the difference from the [§3] cohort is one patient on this
  workbook, and §7.1's audit entry names them. §13 carries the general form forward: the difference is
  one patient *on v7*, and nothing in this stage bounds it on the workbook that follows.

**Why the mask is not folded into `design`.** Because `design` is called by Stages 9 and 12 over
different covariate lists, and the complete-case set is a property of the list. An outcome model over
`outcome_model_covariates("tici_2b_3")` — treatment + `center` + `atrial_fib`, all complete —
excludes nobody, and a `design` that silently dropped rows would give Stage 9 a design matrix and a
response of different lengths, aligned by luck.

### 4.5 The preconditions [invariant 4]

Before anything is coded. Four checks, collected and raised together, following Stage 2's
`_assert_schema`, Stage 4's `_assert_classifier_inputs` and Stage 5's `_assert_cohort_inputs`:

```python
def _assert_design_inputs(df: pd.DataFrame, covariates: Sequence[str]) -> None:
    bad: list[str] = []

    absent = [c for c in covariates if c not in df.columns]
    if absent:
        bad.append(
            f"D1  {', '.join(absent)}: not a column of the frame. A design matrix over a covariate "
            "the frame does not carry is a model over a different specification than the one "
            "requested. Stage 3's derived columns exist only after derive(); the [§13] subgroups "
            "only after derive_cohort().")

    denied = [c for c in covariates if c in C.POST_TIME_ZERO]
    if denied:
        bad.append(
            f"D2  {', '.join(denied)}: post-time-zero [invariant 4, §12]. This is the assertion "
            "against the explicit denylist that invariant 4 requires, and this function is where it "
            "can be made: every model's covariates pass through here. onset_to_groin_min is "
            "reported by arm and never adjusted for [§12]; an outcome in a design matrix is the "
            "outcome conditioning on itself.")

    for c in (c for c in covariates if c in C.CATEGORICAL and c in df.columns):
        absent = df[c].isna()
        if absent.any():
            bad.append(
                f"D4  {c}: {int(absent.sum())} record(s) carry no value. A missing factor value is NOT "
                "caught downstream and this is the only place it can be: the Categorical conversion "
                "maps it to NaN, get_dummies emits an ALL-ZERO row, and that row is byte-identical to "
                "a record genuinely at the reference level — measured [Stage 6 §4.5]. There is no nan "
                "left for F4 to find and D3 excludes missing values from its own check by design. Use "
                "model.complete_cases() [§4.4]; every caller owns its own mask [§9].")

        offside = ~df[c].isin(C.FACTOR_LEVELS[c]) & df[c].notna()
        if not offside.any():
            continue
        # `case_id` is used to NAME the offenders and is not a precondition of building a design.
        # Indexing it unconditionally made `design` raise a bare KeyError on any frame without it —
        # including the single-centre and factor-only frames §12.1 describes, and anything Stage 9
        # or Stage 12 assembles [§20 round 3, finding 5]. The count is the assertion; the names are
        # a courtesy, and a courtesy may not impose a schema.
        named = (", ".join(sorted(df.loc[offside, "case_id"])) if "case_id" in df.columns
                 else f"index {sorted(df.index[offside])} — the frame carries no case_id")
        off = df.index[offside]
        if len(off):
            bad.append(
                f"D3  {c}: {len(off)} record(s) carry a level outside FACTOR_LEVELS[{c!r}]: "
                f"{named}. The Categorical conversion turns such a value into NaN "
                f"and get_dummies then encodes it as an all-zero row — which is the encoding of the "
                f"reference level, {C.REFERENCE_LEVELS[c]!r}. The record would be modelled as though "
                "it were at the reference, with no missing value anywhere and nothing to notice "
                "[Stage 6 §3.1].")

    if bad:
        raise C.SchemaError(
            "\n  " + "\n  ".join(bad)
            + f"\n  {len(bad)} design assertion(s) failed over {len(df)} records and "
            f"{len(covariates)} covariate(s).")
```

- **D2 is the one that earns this function.** Invariant 4 has been declared since Stage 1 and
  assertable nowhere: `POST_TIME_ZERO` existed, and no code path took every model's covariates. This
  one does, and Stages 8, 9 and 12 inherit the check by calling `design` — which is why §9 states that
  they must not build a matrix themselves, and §12.6 scans for it.
- **D3 excludes missing values from its own check** (`& df[c].notna()`) because they are D4's, not
  because they are harmless. Without the `notna` clause D3 would fire on every incomplete record and the
  two failures — a level outside the declared set, and no level at all — would report as one, with one
  message describing the wrong mechanism.
- **D4 exists because nothing else catches a missing factor value.** Measured (§18f): `pd.Categorical`
  maps `pd.NA` to `NaN`, `get_dummies` emits an all-zero row, and that row is **byte-identical** to a
  record genuinely at the reference level — both `[0. 0. 0.]` on the cohort's `center` dummies. So:
  - **F4 cannot see it.** There is no non-finite value left in the design; the missingness was consumed
    by the encoding.
  - **D3 cannot see it**, by its own `& notna()`.
  - **F5 cannot see it.** The fitted probability is finite and interior; it is simply the probability of
    a patient the model believes is at HUG.

  This is §3.1's "quieter and worse" failure in its third form, and the only guard that ever stood in its
  way is `complete_cases` — which §9 makes the caller's, and which `propensity.fit` runs but Stages 8, 9
  and 12 must each remember. An earlier version of this bullet asserted that F4 was the backstop; it is
  not, and D4 is what makes the sentence true rather than what makes it unnecessary. D4 never fires on
  the [§7] path, because `fit` masks first.
- **D1 fires before D3 can `KeyError` on a covariate.** The generator in D3's loop re-filters on
  `c in df.columns` for that reason; the checks are collected, so D1 must not prevent D3 from running.
- **D3 does not require a `case_id` column, and an earlier draft did.** It indexed
  `df.loc[..., "case_id"]` *before* its own `if len(off)` guard, so `design` raised a bare
  `KeyError: 'case_id'` on any frame carrying a declared factor and no identifier — measured, on a clean
  two-row frame. That is an undeclared schema requirement on the one function every model's covariates
  pass through (§9), and `design`'s signature and docstring promise nothing of the kind. §12.1's own
  single-centre and `PS_COVARIATES_FULL` frames would have hit it. The identifiers are now read only
  when there is something to name, and the message falls back to the index when the frame has no
  `case_id` — which is the right trade: a check that *names* patients is better than one that does not,
  but not at the price of refusing to run.
- **Every branch is unreachable on v7**, and that is recorded rather than treated as a reason to skip
  one (§18): every cohort `center` is in `CENTER_ORDER`, every `onset_type` is in its declared levels,
  and `PS_COVARIATES` is disjoint from `POST_TIME_ZERO` — which `test_config.py` already asserts.
  Like Stage 3's `_assert_onset_flags` and Stage 5's C1-C4, they are written for the workbook that has
  not arrived yet, and for the covariate list a later stage passes in.

## 5. The Firth fit [§7]

### 5.1 What is being maximised, and why it is not a library call

[§7] prescribes one estimator: Firth-penalised logistic regression, refitted identically in every
bootstrap replicate. The penalised log-likelihood is

```
    l*(b)  =  l(b)  +  ½ log|I(b)|,        I(b) = X' W X,  W = diag(p(1 − p))
```

whose stationary point solves the **modified score** equation

```
    X'(y − p + h(0.5 − p))  =  0,          h = diag(H),  H = W^½ X (X' W X)⁻¹ X' W^½
```

The `½ log|I(b)|` term is the Jeffreys prior, and `h(0.5 − p)` is what pulls the fit back from a
boundary an unpenalised likelihood would run to. That is the entire reason [§7] names it: **an
unpenalised MLE fails to converge in a minority of sparse replicates**, and Stage 10 refits 2000
times.

§16b is why this is implemented here rather than imported.

### 5.2 The estimator, written out

Written out in full, as Stage 5 §5.1 and §8 are: this is the function every reported number passes
through, and a reader reconstructing it from prose would be reconstructing the one thing that must not
be approximate.

```python
def firth(X: pd.DataFrame, y: np.ndarray) -> Fit:
    """Firth-penalised logistic regression [§7]. Raises FitError; never returns a fallback.

    `X` carries no intercept — one is prepended here, so a design matrix is never ambiguous about
    whether it has one, and `beta[0]` is always the intercept.
    """
    Xc = np.column_stack([np.ones(len(y)), X.to_numpy(dtype=float)])
    y = np.asarray(y, dtype=float)
    _assert_fittable(Xc, y, tuple(X.columns))            # §5.4
    beta = np.zeros(Xc.shape[1])
    ll_old = _penalised_loglik(Xc, y, beta)
    first_step_norm, rescales, halvings = 0.0, 0, 0      # §3.2; §12.4a asserts on all three

    for iteration in range(1, C.FIRTH_MAX_ITER + 1):
        eta = np.clip(Xc @ beta, -C.FIRTH_ETA_CLIP, C.FIRTH_ETA_CLIP)
        p = 1.0 / (1.0 + np.exp(-eta))
        w = np.clip(p * (1.0 - p), C.FIRTH_WEIGHT_FLOOR, None)
        Xw = Xc * np.sqrt(w)[:, None]
        info_inv = np.linalg.inv(Xw.T @ Xw)              # LinAlgError is NOT caught — §5.4
        h = np.clip(np.einsum("ij,jk,ik->i", Xw, info_inv, Xw), 0.0, 1.0)
        score = Xc.T @ (y - p + h * (0.5 - p))
        step = info_inv @ score

        norm = float(np.linalg.norm(step))
        if iteration == 1:
            first_step_norm = norm
        # RELATIVE to the current iterate, never absolute — §5.2a. A trust region, not a rejection:
        # the Newton direction is kept and only its length is bounded.
        size = norm / max(1.0, float(np.linalg.norm(beta)))
        if size > C.FIRTH_MAX_STEP:
            step = step * (C.FIRTH_MAX_STEP / size)
            rescales += 1

        for _ in range(C.FIRTH_MAX_HALVINGS):
            ll_new = _penalised_loglik(Xc, y, beta + step)
            if ll_new >= ll_old:
                break
            step = step / 2.0
            halvings += 1
        else:
            raise FitError(
                f"Firth: step-halving exhausted {C.FIRTH_MAX_HALVINGS} halvings at iteration "
                f"{iteration} without increasing the penalised likelihood. [§7] prescribes one "
                "estimator; there is no second one to try.")

        beta = beta + step
        moved = abs(ll_new - ll_old)                     # BEFORE ll_old is rebound — §5.2b
        if moved < C.FIRTH_TOL:
            return Fit(beta, _probabilities(Xc, beta), iteration, "likelihood", tuple(X.columns),
                       first_step_norm, rescales, halvings)
        if float(np.max(np.abs(score))) < C.FIRTH_SCORE_TOL:
            return Fit(beta, _probabilities(Xc, beta), iteration, "score", tuple(X.columns),
                       first_step_norm, rescales, halvings)
        ll_old = ll_new

    raise FitError(
        f"Firth: no convergence in {C.FIRTH_MAX_ITER} iterations. The penalised likelihood moved "
        f"{moved:.3g} on the last step against a tolerance of {C.FIRTH_TOL:g}, and the largest "
        f"modified score component was {float(np.max(np.abs(score))):.3g} against "
        f"{C.FIRTH_SCORE_TOL:g}. {rescales} step(s) were shortened by the trust region and {halvings} "
        f"halving(s) were taken [§5.2a]. [§7] prescribes one estimator, and [§10] drops and counts "
        "the replicate rather than substituting another.")
```

### 5.2a The trust region is relative to the iterate, and that is a correction

An earlier version of this document bounded the step **absolutely**: `if ‖step‖ > FIRTH_MAX_STEP`. That is
wrong, and the argument is not a matter of taste.

**A Firth estimate is equivariant under rescaling a covariate; an absolute step bound is not.** Replace a
covariate `x` by `1000·x` and the maximum-penalised-likelihood `β` for it is the old one divided by 1000,
with identical fitted probabilities, identical weights and an identical ESS. The estimate does not know
what units the data was recorded in. But `FIRTH_MAX_STEP × FIRTH_MAX_ITER = 5.0 × 200 = 1000` caps how far the
coefficient vector can travel in total, in absolute terms — so a covariate whose natural coefficient is
large *because the covariate itself is small* can need more travel than the budget allows, and the fit is
reported as a failure when a perfectly good finite optimum exists.

Measured (§18d), on the four-record design §12.0 1b calls `travels_far()`:

```
    optimum is finite:  β ≈ (0, 2007.9), penalised log-likelihood −8.6292176
    absolute radius:    FAILS at FIRTH_MAX_ITER = 200        ← and 1000 is exactly the travel budget
    the same data with the covariate × 1000:  converges in 9 iterations
    the same data with the covariate × 0.001: FAILS
```

**Why this is not a curiosity.** [§10] drops and counts failed replicates. If failure depends on covariate
units, *which* replicates are dropped is partly an artifact of how the data was recorded — selection on
the replicate, which is precisely the hazard §5.3 invokes to justify excluding the step norm from the
convergence test. Stage 11's subgroups and [§13]'s sensitivity specifications are where sparse designs
first arise, and nothing in the [§6] covariate list guarantees a benign scale forever.

**The fix, and what it costs.** `size = ‖step‖ / max(1, ‖β‖)`. Bounded relative to the current iterate, so
the budget grows with the coefficient vector it is bounding. Measured (§18d): **the workbook fit does not
move at all** — 7 iterations, `max|Δβ| = 0`, `e` identical to every digit §18 records, ESS
30.3394 / 29.1990 — and **§12.0.2's golden vector is unchanged to all six pinned decimals**, at 11
iterations and one rescale exactly as before. `separated()`, `well_behaved()` and `needs_halving()` are
untouched. What changes is that `travels_far()` now converges, at iteration 10, and the scale test
converges at every multiplier from 1e-3 to 1e6.

**A more principled alternative was measured and declined.** Bounding the **Newton decrement**,
`√(stepᵀ I step)`, is *exactly* invariant under any linear reparameterisation of the design rather than
approximately so, and it is the textbook trust region for a Newton method. Measured (§18d): it also fixes
both problems, in fewer iterations — but it perturbs the workbook fit at 6.0e-06 in the coefficients and
moves §12.0.2's golden vector in the fifth decimal, and it binds four times on `well_behaved()`, a
well-conditioned 4000-record design where no safeguard should fire. It buys exactness in a regime this
pipeline never enters, at the cost of re-pinning every number this document has measured. The relative
form is kept, the decrement is recorded here, and if a future workbook makes the distinction bite this
section is where to start.

**This is a [§13] amendment to the estimator, not a repair**, on §5.6's own terms: the trust region is
part of what is refitted in every replicate. It is recorded as an amendment with its measurement, made
before any estimate exists — which is the cheapest this could ever be done, and the reason to do it now
rather than to file it.

### 5.2b The non-convergence message reported zero, always

The loop's last statement is `ll_old = ll_new`. Falling out of the iteration range therefore reached the
final raise with the two equal, so `abs(ll_new - ll_old)` was **identically zero on every input** — and
the message read "the penalised likelihood moved 0 against a tolerance of 1e-08" inside a *non*-convergence
error, stating a movement below the tolerance it claims was not met. §12.4a asserts on that number, so the
assertion could not fail either.

`moved` is now bound before `ll_old` is rebound. The message also carries the rescale and halving counts,
because after §5.2a those are what distinguish "this replicate is genuinely separated" from "the trust
region ran out of budget" — which is the question Stage 10's operator is actually asking.

### 5.2c The penalised likelihood is not concave, and this loop can converge to a local maximum

**The one failure mode in this stage that is not guarded, stated plainly because it is not.** `l(b)` is
concave; `½log|I(b)|` is not, and their sum is not. So the iteration above — Newton from `beta = 0`,
monotone step-halving — finds *a* stationary point of the penalised likelihood, and [§7] prescribes
*the maximiser*. On some designs those differ.

**Measured (§18e), and both measurements matter.**

```
  the penalty is non-concave:  d²/ds²  of ½log|I| along a slope ray
      slope 0.5 → −0.66      slope 2 → +0.20      slope 5 → +0.019     slope 20 → +0.0047

  synthetic completely-separated designs, standardised covariates, n∈[20,60), k∈[1,4)
      4 of 600 land below what other starts of this same loop reach
      worst gap 0.29 nats,  max|Δe| 0.40

  2000 centre-stratified bootstrap replicates of the ACTUAL 92-record cohort
      1 of 2000,  gap 1.0e-02 nats        (a second, independent 2000-replicate sweep: also 1)
      the other 1999 sit at the optimum to ~1e-10; median gap 4.2e-11, 99th pct 1.2e-09
```

**Every guard passes on these fits.** `e` is strictly interior so F5 is silent; the design is full rank
so F3 is silent; `converged_on` is `"likelihood"`, `rescales` and `halvings` are zero. The three counters
§3.2 added so Stage 10 could audit its fits all read clean. There is no symptom.

**The v7 point fit is not affected.** A 29-start search puts it at the global optimum (gap 1.3e-10,
`max|Δβ| = 1.7e-07`), which is what a design at condition number 997 converging in 7 iterations should
do. Every number this document reports is sound. The exposure is entirely [§10]'s refits and Stage 11's
subgroups.

**Four detectors were designed and measured, and none works inside `firth`.** This is recorded so nobody
rebuilds them:

| Candidate check | Detected, on the clean 2000-replicate sweep |
|---|---|
| 5 structured starts — the OLS direction on `2y−1`, ±, ×1/×10/×100 | 0 of 1 |
| 48 coordinate-wise starts — `±c·e_j`, `c ∈ {1,10}` | 0 of 1 |
| 12 random starts seeded once from `C.SEED` | 0 of 1 |
| 12-start multistart used as a *corrector* rather than a check | recovered 1 of 4 synthetic cases; 12 of 12 starts converged to the same wrong point on the other three |

The basins are narrow and their locations depend on the design, so only probes that **vary per fit** find
them — and `firth` fits one design, holds no seed, and cannot draw varying probes without either
breaking determinism or inventing a seed from the data's own bytes, which is a worse defect than the one
it would catch. A check measured to fire 0 times out of 1 is not a check; it is a green light over an
unguarded path, which is the failure class every other section of this document exists to remove.

**So `firth` is unchanged and Stage 10 owns the detection** (§9, §13). That is where the replicate loop
and `C.SEED` already live, where the exposure entirely is, and where varying probes are natural rather
than contrived. Nothing about the estimator changes; §12.0.2's golden vector, §18's every figure and the
Definition of done stand exactly as measured.

**What this costs to leave unguarded, stated so the trade is visible.** Roughly one replicate in two
thousand carries a propensity score that is not the [§7] estimator, silently, inside a percentile
interval built from two thousand. The prevalence rests on two observed events, so a Poisson interval on
the rate spans about 0 to 6 per 2000. The effect on a percentile bound is far below its own Monte Carlo
error at `N_BOOT = 2000` — which is the reason this is a documented gap with an owner rather than a
blocker, and the reason it is written out at this length rather than in a footnote.

```python
def _penalised_loglik(Xc: np.ndarray, y: np.ndarray, beta: np.ndarray) -> float:
    """l(b) + ½log|I(b)|, or -inf where the information matrix is singular.

    -inf rather than a raise, so that step-halving can *reject* a step and try a shorter one; a raise
    there would turn a recoverable step into a failed fit. slogdet returns -inf for a singular matrix
    rather than raising [§3.1], so the branch is about a non-finite logdet and not an exception.

    THE BRANCH DOES NOT FIRE, and that is measured rather than hoped [§3.1, §18]. FIRTH_WEIGHT_FLOOR
    floors w, so X'WX >= FIRTH_WEIGHT_FLOOR * X'X, which is nonsingular whenever X is full rank — and F3
    established that before the loop began. A design scaled small enough to underflow the determinant
    loses rank under matrix_rank first, so F3 raises before this function is ever called. It is kept
    because it is the contract between the floor and the halving loop: remove the floor, or weaken F3
    to a warning, and this branch becomes the thing that stops a crash. Do not delete it as dead code
    without deleting the reason it is dead.
    """
    eta = np.clip(Xc @ beta, -C.FIRTH_ETA_CLIP, C.FIRTH_ETA_CLIP)
    ll = float(np.sum(y * eta - np.logaddexp(0.0, eta)))
    w = np.clip(1.0 / (1.0 + np.exp(-eta)) * (1.0 - 1.0 / (1.0 + np.exp(-eta))),
                C.FIRTH_WEIGHT_FLOOR, None)
    _, logdet = np.linalg.slogdet((Xc * w[:, None]).T @ Xc)
    return ll + 0.5 * logdet if np.isfinite(logdet) else -np.inf
```

```python
def _probabilities(Xc: np.ndarray, beta: np.ndarray) -> np.ndarray:
    """The fitted probabilities, on the SAME clipped linear predictor the loop used.

    Written out rather than inlined because it is called on both convergence routes and F5 asserts
    strict interiority on its output (§5.5): if this applied a different clip — or none — from the one
    inside the loop, F5 would be testing a quantity the weights are not computed from, and a fit whose
    linear predictor hit FIRTH_ETA_CLIP could pass F5 while producing a weight of exactly 0.0.
    """
    return 1.0 / (1.0 + np.exp(-np.clip(Xc @ beta, -C.FIRTH_ETA_CLIP, C.FIRTH_ETA_CLIP)))
```

Six numerical details, all of them lifted from `pilots/analysis.py:346-393` and all of them
load-bearing:

- **`np.logaddexp(0.0, eta)`, never `np.log(1 + np.exp(eta))`.** The second overflows at `eta ≈ 710`
  and loses precision long before; the clip at ±`FIRTH_ETA_CLIP` keeps `exp` in range, and `logaddexp`
  is exact at the ends.
- **`beta` starts at zero**, so the first `p` is 0.5 everywhere and the first information matrix is
  `X'X/4` — the best-conditioned it will ever be, and the iteration where a rank deficiency shows.
- **The hat diagonal is computed on the square-root-weighted design**, via `einsum`, not as
  `np.diag(H)`: the full `H` is n × n and only its diagonal is used.
- **The step is rescaled, not rejected, when its norm exceeds `FIRTH_MAX_STEP`.** A trust region; the
  direction is a Newton direction and is kept, only its length is bounded, so a first iteration far
  from the optimum cannot leap into a region where the likelihood is `-inf` and burn all thirty
  halvings.
- **Step-halving is on the penalised likelihood, not on the unpenalised one.** Halving on `l(b)`
  would drive the fit toward the boundary the penalty exists to hold it back from, and it would do so
  while converging.
- **`ll_old` is updated only when the loop continues**, so the acceptance test at the next iteration
  compares against the last *accepted* value. The pilots' asymmetry here is deliberate and is kept.

### 5.3 The convergence criterion, and why the step norm is not in it

Two criteria, either one sufficient, and the `Fit` records which fired:

```
    |Δ l*(b)| < FIRTH_TOL              →  converged_on = "likelihood"
    max |modified score| < FIRTH_SCORE_TOL  →  converged_on = "score"
```

**Excluding the step norm is what `logistf` is believed to do — carried from the pilots' comment and
NOT verified against the package (§16b.2, §13, DoD-16) — and the reason stands on its own merits.**
Under near-separation the penalised surface is genuinely flat: the fit is finite and correct, and the
Newton step keeps drifting along a ridge without changing the objective. Requiring `‖step‖ → 0` as well
would report that fit as a failure — and [§10] drops failed replicates, so the ones dropped would be
exactly the sparse ones, which is selection on the replicate. That argument does not depend on
`logistf`, which is why the unverified attribution is a sentence to strike rather than a decision to
revisit if the oracle disagrees.

**Two facts about the score route, both measured, and neither of them was known when the two criteria
were written down.**

- **The score is evaluated at the iterate the step was computed *from*, not at the returned `beta`.**
  In §5.2 `score` is bound at the top of the iteration from the pre-step `beta`, and both convergence
  tests run after `beta = beta + step`. So a `Fit` reporting `converged_on = "score"` is reporting a
  flat score one step behind its own coefficients. Harmless — when the score is below `FIRTH_SCORE_TOL` the
  step is negligible — but it is stated here because §12 asserts `converged_on` as a property of the
  returned fit, and because it is one of the two things the R oracle should compare (DoD-16).
- **The score route fires exactly where it was designed to: on an ill-conditioned information matrix.**
  It is reached when `‖I⁻¹‖` is large, which near-separation and small covariate scales both produce.
  Measured (§18d) under §5.2a's relative trust radius:

  ```
    design                        ‖I⁻¹‖      converged_on
    ────────────────────────────────────────────────────────
    separated(), x = 1…10          1.1e+01   likelihood
    big_first_step(), x ~ 5e-2     1.5e+03   likelihood
    WORKBOOK, 92 × 11              small     likelihood   (7 iterations)
    separated() with x × 1e-3      3.3e+05   SCORE        (11 iterations)
    travels_far(), x ~ 1e-3        3.3e+06   SCORE        (10 iterations)
    big_first_step() with x × 1e-3 1.4e+09   SCORE         (8 iterations)
  ```

  The mechanism, and it is worth writing down because an earlier draft of this bullet got it backwards.
  Near the optimum `Δl* ≈ ½·sᵀI⁻¹s`. That draft reasoned that `max|s| < FIRTH_SCORE_TOL = 1e-6` therefore
  implies `Δl*` of order `1e-12`, three orders below `FIRTH_TOL = 1e-8`, so the likelihood test — which
  §5.2 runs first — always fires first and the score route is dead. **The step it skipped is that
  `Δl* ~ ½·s²·‖I⁻¹‖`, and it assumed `‖I⁻¹‖ ~ 1.`** At `‖I⁻¹‖ = 1e6`, `|s| = 1e-6` gives `Δl* ~ 5e-7`,
  comfortably *above* `FIRTH_TOL` — so the likelihood is still moving when the score has gone flat, and the
  second criterion fires. A nearly-singular information matrix and a flat penalised surface are the same
  fact stated twice, which is precisely the condition [§7] chose Firth for.

**How the earlier claim survived measurement, which is the more useful lesson.** It was checked against
fourteen designs and held on all fourteen — because every one of them was well-scaled, and because the
one badly-scaled design in the set never converged at all under the **absolute** trust radius §5.2a
replaced. Fixing the trust region is what made the score route observable. A property confirmed on
fourteen fixtures that share a hidden characteristic is confirmed about the characteristic, not about
the property.

**Nothing about the criteria themselves is changed.** A tolerance is part of the estimator (§5.6) and so
is the order the two tests run in; [§10] refits in every one of `N_BOOT` replicates, so reordering them
would change the sampling distribution and is a [§13] amendment, not a repair. What changed is that the
document no longer claims a route is dead when it is not.

**`converged_on` is therefore a live diagnostic, not a tripwire.** It costs one string on a frozen
dataclass, and §13's instruction to Stage 10 — surface the distribution of the two routes across
replicates, because a shift toward `"score"` is the signature of near-separation becoming common — is a
real instruction with a real expected signal rather than a watch on a constant. Measured on v7:
`"likelihood"`, at iteration 7 (§18). A replicate that comes back `"score"` is telling you its design is
close to singular, and §12.4a's fixtures are how a reader learns to recognise that before Stage 10 does.

### 5.4 A singular information matrix is a rank failure, not a numerical hiccup

```python
def _assert_fittable(Xc: np.ndarray, y: np.ndarray, columns: tuple[str, ...]) -> None:
    """F3, F4 and F6 — the three things that make `firth` return, or blame, the wrong answer [§3.1].

    Ordered design → response → rank: F4 and F6 are O(n) and their failures are the common ones, and
    F3's matrix_rank is an SVD, so it goes last. Unlike D1-D4 these raise at the first failure rather
    than collecting, because each one makes the next check meaningless — a rank computed over a design
    containing nan is not a rank.
    """
    if not np.all(np.isfinite(Xc)):
        bad = [columns[j - 1] if j else "intercept"
               for j in np.unique(np.argwhere(~np.isfinite(Xc))[:, 1])]
        raise FitError(
            f"F4  the design carries a non-finite value in: {', '.join(bad)}. An Int64 column with "
            "pd.NA converts to nan SILENTLY [§3.1], the fit then runs to completion and returns "
            "all-nan coefficients with a plausible iteration count, and the first thing anyone sees "
            "is a [nan, nan] interval in Stage 10. Use model.complete_cases() [§4.4].")

    if not np.all(np.isfinite(y)):
        raise FitError(
            f"F6  the response carries {int((~np.isfinite(y)).sum())} non-finite value(s) of "
            f"{len(y)}. Every candidate step is then nan, no step is ever accepted, and the fit "
            "raises step-halving-exhausted at iteration 1 — a message about the estimator for a "
            "failure of the data [§3.1]. The response is complete-case per estimate [§11] and the "
            "mask is the caller's: this module does not know what the response is [Stage 6 §4.4, §9].")

    off = np.unique(y[~np.isin(y, (0.0, 1.0))])
    if off.size:
        raise FitError(
            f"F6  the response takes {off.size} value(s) outside {{0, 1}}: "
            f"{', '.join(f'{v:g}' for v in off[:10])}. This is a logistic fit; a response on any "
            "other scale produces confident coefficients for a model of something else, with no "
            "missing value and nothing to notice. An ordinal response belongs to a "
            "proportional-odds fit and does not come through here [Stage 6 §5.4, §9].")

    rank = int(np.linalg.matrix_rank(Xc))
    if rank < Xc.shape[1]:
        raise FitError(
            f"F3  the design has rank {rank} of {Xc.shape[1]} columns: "
            f"{', '.join(('intercept',) + columns)}. The information matrix is singular at every "
            "beta, so no penalised optimum is identified. The constant-column rule [§4.3] removes an "
            "absent factor level; exact collinearity between two covariates is a [§6] "
            "specification error — penumbra_ml with core_ml and tmax6_ml is the one [§6] names.")
```

**`np.linalg.inv` inside the loop is deliberately not wrapped.** `pilots/analysis.py:366-369` catches
`LinAlgError` and falls back to `np.linalg.pinv`, and that is a fallback of exactly the kind
[§7] forbids — not to a different estimator, but to a different *estimand*: the pseudo-inverse
silently picks the minimum-norm solution among infinitely many, so the fit returns coefficients for a
design that identifies none. Measured (§18): `inv` **does** raise on an exactly singular information
matrix, so the pinv branch is what hides it.

F3's rank check is what makes the raise informative rather than a bare `LinAlgError` from inside a
loop. It uses `matrix_rank`, which is tolerance-based, and so it also catches the *near*-collinearity
that `inv` would happily invert at a condition number of 1e17 — measured at 8.8e17 for an exactly
duplicated column (§18). The cohort's own design is rank 12 of 12 at a condition number of 997.

**F6 is the response's F4, and it exists because of where §0.1 put the module boundary.** `model.py` is
outcome-agnostic by design, so it is the one place a response arrives with nobody having checked it:
`propensity.fit`'s F1 checks `TREATMENT` before calling `firth`, but §9 states that Stages 8, 9 and 12
enter `model.py` **directly**, and Stage 9's responses are outcomes with genuine missingness —
`mrs_90d`, `sich` and `ph2` all carry it. Two failures, and only the second is silent:

- **A non-finite response raises today, with the wrong message.** Measured (§18): every candidate step
  is `nan`, so no step is ever accepted and the loop exhausts its halvings at iteration 1, reporting
  "[§7] prescribes one estimator; there is no second one to try." A reader follows that to Stage 10 and
  sparse replicates. F6 costs one line and names the actual cause.
- **A response outside {0, 1} does not raise at all.** Nothing anywhere in `model.py` or
  `propensity.py` establishes that `y` is a dichotomy, so `firth(X, mrs_90d)` fits, converges and
  returns coefficients — a logistic model of a seven-level ordinal, with no missing value and no
  warning. This is `mrs_90d`'s own hazard rather than a hypothetical: it is the [§8] primary outcome,
  it is the one variable in the pipeline most likely to be handed to a fitter by mistake, and [§8]'s
  proportional-odds model is a different function that does not exist yet.

`y` was already in `_assert_fittable`'s signature and unused, which is the shape this check was always
supposed to have.

**F6's messages name neither the exposure nor any outcome, and that is a constraint rather than a
stylistic choice.** Writing it revealed the conflict: the natural message for the first branch is "use
`propensity.fit`'s F1, or your outcome's own mask", and for the second "`mrs_90d` has seven levels and
would fit silently" — and both put the name of a specific response inside the module §0.1 defines as
outcome-agnostic, which §12.12's scan exists to prevent. A module that names `TREATMENT` in a string is
one grep away from being a module that knows what the exposure is, and the next reader to add a check
will follow the precedent. So the concrete examples live **here**, in the spec prose above, where they
are more useful anyway, and the runtime messages say "this module does not know what the response is".
§12.12 records the rule this produced.

### 5.5 Raise on failure — no fallback estimator anywhere [invariant 5]

```python
class FitError(RuntimeError):
    """The prespecified fit did not converge on this sample. [§7]: there is no second estimator."""
```

One exception type, from `model.py`, and it is what Stage 10 catches to drop and count a replicate.

**The evidence is in this repository, not in principle.** `pilots/analysis.py`'s docstrings record
that SAP v1.0's scikit-learn fallback fired in **16.6% of bootstrap replicates** — so one interval in
six was a mixture of two estimators' sampling distributions rather than one estimator's. [§7] is
written against that measurement: "switching estimators on failure makes the bootstrap a mixture of
two estimators rather than the sampling distribution of one."

Four things follow, and §12.7 asserts all four:

- **No `except` clause in `model.py` or `propensity.py` returns a value.** `_penalised_loglik`'s
  `-inf` branch is not an exception handler; it is a value the step-halving loop is designed to
  reject (§5.2).
- **`fit_ps_mle`, `fit_ps_ridge` and `_sklearn_ps` are not ported.** They cannot be — scikit-learn is
  not a dependency of this project, by Stage 1's decision — and §16 records that this is a happy
  coincidence rather than the reason.
- **The sensitivity specifications of [§13] do not get a different estimator.** They are deferred
  (roadmap Stage 11), and when they arrive they are different *covariate lists* under the same
  estimator, which is what makes them different target populations rather than competing estimates
  (§6.4).
- **A test asserts that no shipped module names a second estimator**, by AST scan over the imports of
  `model.py` and `propensity.py`: `sklearn`, `statsmodels` and `scipy` may not appear in either
  (§12.7, first bullet — **not** §12.12, which is the module-boundary scan).

### 5.6 The seven numerics move into `config.py`

Every tolerance above is a prespecified analytical choice, and Stage 1's rule is that those live in
`config.py`. Declared as a block, with the reason each has the value it has:

```python
# --- the Firth fit [§7, Stage 6 §5] -----------------------------------------------------------
#
# Prespecified, and they are part of the estimator: [§10] refits this model in every one of N_BOOT
# replicates, so a tolerance is a property of the sampling distribution and not a runtime knob. The
# ORDER the two convergence tests run in is prespecified for the same reason [Stage 6 §5.3].
#
# These are the PILOTS' values, carried from pilots/analysis.py. An earlier draft of this block said
# they are `logistf`'s defaults; that is NOT verified and the claim is withdrawn until DoD-16 records
# it, because `logistf.control()` must be passed these values explicitly for §16b's comparison to be
# like-for-like anyway — so what makes the comparison meaningful is the explicit pass, not a shared
# default [Stage 6 §16b.2].

FIRTH_MAX_ITER: Final[int] = 200          # measured: 7 on the cohort, 11 on the hand frame
FIRTH_TOL: Final[float] = 1e-8            # |Δ penalised log-likelihood| — the route that always fires
FIRTH_SCORE_TOL: Final[float] = 1e-6      # max |modified score|; fires on ill-conditioned I — §5.3
FIRTH_MAX_HALVINGS: Final[int] = 30       # exhausting them raises; it is not a convergence route
FIRTH_MAX_STEP: Final[float] = 5.0        # trust radius, RELATIVE to ‖beta‖ — see Stage 6 §5.2a
FIRTH_ETA_CLIP: Final[float] = 500.0      # keeps exp() in range; reachable at the upper end [§3.1]
FIRTH_WEIGHT_FLOOR: Final[float] = 1e-10  # floors p(1-p) so the information matrix stays invertible
```

**Which of the seven are reached, measured (§18, §18d, §18f) rather than assumed.** `FIRTH_TOL`,
`FIRTH_MAX_STEP` and `FIRTH_ETA_CLIP` bind on some tested design; **`FIRTH_SCORE_TOL` binds on the
small-scale fixtures** (§5.3); `FIRTH_MAX_ITER` and `FIRTH_MAX_HALVINGS` bind only when lowered by a test,
since after §5.2a no fixture exhausts either by data alone.

**`FIRTH_WEIGHT_FLOOR` is the one whose scope must be stated carefully.** Its clip is never active on any
fixture in `test_model.py` or on the cohort design — `min p(1−p)` is 0.0124 across the fixtures and 0.0262
on the workbook, eight orders above the floor. But it is **not** unreachable in general: measured (§18f),
on 800 synthetic completely-separated designs the floor binds on **68 (8.5%)**, and the smallest `p(1−p)`
observed is **exactly 0.0**. So §12.4a's assertion is scoped to the fixtures it names and is a statement
about them, not about the estimator — and `_penalised_loglik`'s `-inf` branch is unreachable *behind F3
and a floor that is holding*, which on those designs it is not.

Two consequences worth stating, because they compound with §5.2c. Where `p(1−p)` reaches exactly 0, `p` is
exactly 0 or 1, so **F5 fires** on the [§7] path and the failure is loud — but `model.py`'s direct callers
(§9) have no F5. And where the floor binds without saturating, the information matrix is built from
floored weights while the score is the gradient of the **unfloored** objective, so the iteration's fixed
point maximises neither; that is a plausible mechanism behind §5.2c's local-maximum cases, two of whose
four had `min p(1−p)` below the floor. Neither is acted on here: both are properties of designs this
pipeline does not fit, and §9 hands the general case to Stage 10 with the measurement.

**`FIRTH_MAX_STEP` is the one whose *meaning* changed after this document was first written**, and the
change is §5.2a's: it bounds `‖step‖ / max(1, ‖β‖)` rather than `‖step‖`. The number is unchanged and the
workbook fit is unchanged to every digit; what changed is that the bound no longer depends on the units
the covariates were recorded in. A reader who remembers the absolute form should read §5.2a rather than
assume a typo.

**One departure from the pilots, and it is deliberate.** `pilots/analysis.py` floors the working
weight at `1e-12` inside its likelihood and `1e-10` in its update — two different values for one
quantity, with no stated reason. One constant, used in both places. Verified (§18) that the coefficient
agreement with `firthlogist` is unaffected to 2.7e-15, so the pilots' asymmetry was inert; it is
removed because a reader cannot tell that from the code.

**`PROB_CLIP` is deliberately absent.** `pilots/analysis.py:82` clips `e` into
`[1e-8, 1 − 1e-8]`; §5.5's F5 assertion replaces it. §15 records the decision.

## 6. Weights, ESS and the weighted population [§7]

### 6.1 The overlap weights

```python
w = np.where(a == 1, 1.0 - e, e)          # [§7]
```

Two lines of argument, both from [§7], because this is the estimand and not an implementation detail:

- **`w = 1 − e` treated, `w = e` control**, which is the ATO tilt `h(X) = e(X){1 − e(X)}` split
  between the arms. Bounded in `[0, 1]` by construction, and largest at `e = 0.5` — patients near
  equipoise.
- **No trimming, no truncation, no stabilisation.** There is nothing to trim: the weights are bounded
  above by 1 and the pathology trimming exists for — an inverse-probability weight exploding as
  `e → 0` — cannot arise. `pilots/analysis.py:127-131` has a `scheme="iptw"` branch reading
  `IPTW_TRIM`; this SAP has no trimmed-IPTW row anywhere, `config.py` declares no such constant, and
  §16 declines the branch. Measured on v7: the largest weight is **0.965** (§18).

**`w` is a Series on the cohort's index with `nan` where `in_model` is False** (§4.4), so a caller
who forgets the mask gets `nan` in a sum rather than a plausible total. `np.where` above operates on
the complete-case subset; the assembly back onto the full index is in §8.

### 6.2 Kish effective sample size, per arm

```python
def ess(w: pd.Series | np.ndarray) -> float:
    """Kish effective sample size, (Σw)² / Σw². Raises on an empty arm and on a zero-weight arm.

    Reported per arm [§7]. `pilots/analysis.py:140` returns 0.0 when no weight is positive; that is
    declined. An arm whose weights all vanish is an arm in which the fit assigned every patient a
    propensity of exactly 0 or 1 — a structural non-positivity the [§3] restrictions were supposed
    to have removed — and 0.0 is a number every downstream ratio will accept [Stage 6 §3.1].

    `FitError` is qualified as `model.FitError`: it is declared in model.py (§5.5) and propensity.py
    imports the module, not the name (§3). Stage 7 calls this per centre per arm [§9], where the empty
    branch is the one that fires — a centre with no patient in one arm is exactly what
    `treating_centres` excludes, and this raise is what stops Stage 7 plotting it as a zero.
    """
    v = np.asarray(pd.Series(w).dropna(), dtype=float)
    if not v.size:
        raise model.FitError(
            "ESS: this arm has no weighted patient at all, so the effective sample size is 0/0 over "
            "an empty sum. This is NOT the all-zero-weight case below and must not report as it: an "
            "empty arm means the in_model mask lost one, since F2 established both arms were present "
            "before the fit [Stage 6 §8]. Measured: np.sum([] ** 2) is 0.0, so the two reach the "
            "same branch and only a separate check tells them apart [§3.1].")
    total = float(np.sum(v ** 2))
    if not total:
        raise model.FitError(
            "ESS: every weight in this arm is zero, so the effective sample size is 0/0. The fit "
            "put every patient in the arm at a propensity of exactly 0 or 1 — structural "
            "non-positivity, which [§3]'s restrictions remove by design [Stage 5 §5.2 P2].")
    return float(np.sum(v) ** 2 / total)
```

Measured on v7 (§18): **30.34 of 39** treated, **29.20 of 53** control — the control denominator being
53 rather than 54 because of §4.4's one excluded patient, which is exactly why [§11] requires the
denominator beside the estimate.

**Hand-computable, and §12.9 computes it by hand:** `ess([1]*10) = 10.0` exactly, and
`ess([1, 1, 1, 3]) = 6² / 12 = 3.0` exactly. A weighted sample of four with one weight of three is
worth three observations.

### 6.3 The weights do not balance exactly, and that is a property of Firth

Overlap weights have an exact-balance property: with the propensity score from an **unpenalised**
logistic MLE, the weighted means of every covariate in the score are identical across arms, because
that is the score equation. Firth's modified score is a different equation, so the property holds
only approximately.

Measured on v7 (§18): `Σw = 13.63` in the treated arm against `14.11` in the control arm. **Not
equal**, and §12.10 asserts the inequality rather than tolerating it — because a test that asserted
`Σw₁ == Σw₀` to a loose tolerance would pass under Firth today and pass under a silently reinstated
MLE tomorrow, which is the substitution [§7] forbids. The companion assertion is that the same
quantities *are* equal, to 1e-10, when the score comes from `statsmodels.Logit` — so the
near-exactness is measured against its own explanation.

This is Stage 7's business to report [§9], not Stage 6's to fix.

### 6.4 The estimand is indexed by the propensity model, and the log says so

[§7], and it is not a caveat:

> The ATO target population is `h(X) = e(X){1 − e(X)}`, so changing the propensity specification
> changes the *population*, not merely the precision of the estimate.

Two things follow that Stage 6 owns rather than defers:

- **The ESS and the weighted-population description are conditional on this specification**, and all
  three — the ESS, the description, and the conditionality statement itself — are recorded in §7.2's
  `overlap_weights` entry: the ESS in its first table, the description in §7.6's second table, and the
  statement in the `detail`. Roadmap Stage 6's "record in the output" is discharged there. **An earlier
  draft claimed the conditionality sentence alone discharged all of it**; §7.6 is the correction, and
  the reason it matters is that [§7] asks for a description of who the weighted population comprises and
  a caveat about it is not one.
- **`fit` takes no covariate list.** It reads `PS_COVARIATES`, and a caller cannot request a different
  propensity model by passing an argument. `pilots/analysis.py:953` takes `covars or C.PS_COVARIATES`;
  §16 declines it, for Stage 5 §15's reason one stage on: the [§13] propensity-specification rows are
  *different target populations*, and a parameter makes them look like options. When they are
  un-deferred they arrive as a named specification with a [§13] amendment, not as a keyword.

`design`, by contrast, **is** parameterised, because Stage 9's `m_a(X)` and Stage 12's standardisation
model genuinely need other lists — and neither of those is the [§7] estimand.

## 7. The `model` audit kind

### 7.1 `model` is added to `data.py`

Stage 5 §6 added `cohort`; this is the second such amendment and the argument is the same shape. A
fitted model is not a derived column and not a missing value: it is a *claim about the data* with
convergence properties, and it renders under no existing heading.

```python
KINDS: Final[tuple[str, ...]] = (
    "provenance", "contract", "correction", "observation",
    "derivation",
    "cohort",
    "model",
    "structural", "missingness")

_HEADINGS = {
    # … "provenance" through "derivation" unchanged …
    "cohort": "Cohort construction",
    "model":  "Fitted models",      # NEW, in this position
    # … "structural", "missingness" unchanged …
}
```

**Position: after `cohort`, before `structural`.** The document then reads: what was read, what the
contract did, what was corrected, what was observed, what was derived, **who is in the analysis**,
**what was fitted to them**, what is structurally absent, what is missing. The fit follows the
population it was fitted to, which is the only order in which the `n` of §7.1's entry can be checked
against the `n` of Stage 5's `cohort_flow` by eye.

**`_MUST_NAME_CASES` is deliberately unchanged**, and this is the second kind for which that is a
decision rather than an omission (Stage 5 §6.3). `model` carries four entries, of which exactly one
removes patients — §7.1's complete-case entry — and a rule keyed by *kind* would force the other three
either to misrepresent their `n` or to print the whole cohort's identifiers under a coefficient table.
So the rule lives in `propensity.py`, in `_record_exclusion`, where it can be stated exactly and more
strongly: the entry must name as many patients as it says it excluded. The comment above
`_MUST_NAME_CASES` gains one sentence and its count of kinds that legitimately name none goes from
six to seven.

### 7.2 The four entries

| step | kind | `n` | contents |
|---|---|---|---|
| `covariate_completeness` | model | records **excluded** | names every one; table of one row per [§6] covariate: `covariate`, `n_absent`, `excluded_by` — and **no per-centre columns**, for §7.4's reason |
| `design_matrix` | model | **parameters**, including the intercept | table of one row per **level of every declared factor plus every linear covariate**: `column`, `term` (linear / dummy), `status` (`kept` / `dropped: constant` / **`reference (baseline)`**), and for a dummy its factor and level. §7.5 is why the reference rows are there |
| `propensity_fit` | model | records **fitted** | `detail` carries the iteration count, the convergence route (§5.3) and the range of `e`; table of one row per parameter: `parameter`, `coefficient` |
| `overlap_weights` | model | records **weighted** | **two** tables. First, one row per arm plus an `all` row: `arm`, `n`, `sum w`, `ESS`, `max w`, `weight share`. Second, the **weighted-population table** of §7.6: one row per [§6] covariate with its unweighted and ATO-weighted mean in each arm. `detail` carries [§6.4]'s conditionality statement |

Four rules about the tables, three of them inherited:

- **Every declared thing is rendered whether or not the data fills it.** `design_matrix` has a row for
  `center_USZ` reading `dropped: constant`, because a column that vanishes from a table is a column
  nobody knows was considered. `overlap_weights` has a row per `TREATMENT_LABELS` entry, as Stage 4
  §7.1 and Stage 5 §7.1 do, and `pd.crosstab` is used nowhere in this repository for that reason.
- **Every cell is computed from the object it describes, never by subtracting another row.** Stage 5
  §7.3's rule: a table that reconciles with itself by construction has stopped being evidence.
- **`_fmt` is the only float formatter** and `missing` is the only rendering of an absent value, so
  §12.11's byte-identity criterion holds. Coefficients and ESS go through it at six significant
  figures.
- **The `all` row of `overlap_weights` carries `ESS` computed over both arms pooled, and it is
  labelled** — not the sum of the two arm ESSes, which is not a Kish sum of anything.

### 7.3 Where the coefficients live, and why not here

`propensity_fit`'s table carries the twelve fitted coefficients, and **the audit log is the only place
in the repository they appear.** `out/` is gitignored in full, for the reason its own comment gives:
every number in it is patient-derived. A fitted parameter vector is patient-derived, and a
committed-file literal of it would put that in git for the first time.

This is a line worth drawing explicitly, because it is not where Stage 5's line was. Stage 5's spec
quotes cohort counts and the 5.0 mL median, and this one quotes the ESS and the range of `e` — those
are one-number summaries of a population, of the same kind [§16] will publish in a manuscript table.
A twelve-parameter model of the treatment assignment is not a summary; it is the model. So:

- **the coefficients go in the log**, which is reproducible, reviewable and gitignored;
- **the regression oracle is pinned on the synthetic frame**, not on the workbook — §12.0.2's golden
  vector is five fitted probabilities of five `HAND`/`COHORT` fixture records, which contain no patient
  and can therefore live in a committed test and run with no `data/`;
- **the workbook fit is asserted on properties** — iteration count, convergence route, interiority,
  ESS to two decimals — which §12.14 does, data-gated.

### 7.4 What `covariate_completeness` does not render, and why

This entry does **not** carry the per-centre absence columns, because those numbers are already in the
log for this same frame, and `data.absence_by_column` exists to stop them being computed twice. Its own docstring (`data.py:677`) is the argument —

> "Public because Stage 3 calls it too … A second implementation over there would drift on the branch
> that matters — `structural` versus `missing` — which is the exact misreading §10.2 exists to prevent.
> So there is one classifier, called twice with different column lists."

— and it is now called three times, the third by `cohort.py:479` over
`[*sorted(ANALYSIS_NAMES), *sorted(DERIVED_NAMES)]` on the built cohort. Every [§6] covariate is in that
list, so **the per-centre absence counts this entry would render are already rendered**, under
`## Missingness and denominators`, over the same 93 records. Writing them again would be a fourth
hand-rolled copy of `(absent & (df["center"] == centre)).sum()` in the one repository whose acceptance
criterion is a byte-identical log.

So the table carries only what this entry uniquely owns:

| column | contents |
|---|---|
| `covariate` | one row per name in `PS_COVARIATES`, in that order — every declared covariate, including the complete ones |
| `n_absent` | absent on the cohort, which is the same number `absence_by_cohort_column` reports |
| `excluded_by` | whether this covariate is a reason a named patient is outside `in_model` — the fact no other entry carries |

`excluded_by` is the column that earns the entry. `absence_by_cohort_column` answers "what is missing";
only this entry answers "what did the missingness cost **this fit**", which is [§11]'s per-estimate
denominator and the thing §4.4's decision turns on. The entry's `detail` names
`absence_by_cohort_column` as where the per-centre breakdown lives, so a reader is sent one section down
rather than given a second copy, and §12.11 asserts the two agree on `n_absent` for every [§6]
covariate — so the tables cannot drift without a red suite, which is what the duplication was for.

**The alternative was extracting a shared row-builder from `absence_by_column`** — splitting it into a
public table function and a thin `missingness`-recording wrapper. That is the more thorough fix and it
was declined for blast radius: `absence_by_column` is landed, tested, and called from three modules, and
Stage 6 has no other reason to reopen `data.py` beyond §7.1's one-line `KINDS` amendment. It is filed in
TODOS.md rather than done here.

### 7.5 The reference level is a row, because the log is where the coefficients live

`design_matrix`'s table has a row for every **level of every declared factor**, not for every column of
the matrix — and the difference is two rows that an earlier draft omitted.

Reference dummies are dropped by name inside `design` (§4.2), **before** the constant-column rule runs.
So they appear in neither `X.columns` nor `dropped`, and §12.1 asserts exactly that: *"`center_HUG` and
`onset_type_witnessed` are absent from the returned matrix, and neither appears in `dropped`."* A table
built from `kept + dropped` therefore names every column of the model and **none of its baselines**.

That is a worse omission here than it would be anywhere else, because of §7.3: the audit log is the only
place in the repository the fitted coefficients ever appear. `center_CHUV = 0.41` means *relative to HUG*;
without HUG on the page, the number is uninterpretable, and Definition of done 14 asks a human to read
this log end to end. Stage 14 reports from these entries [§16], so a manuscript table of centre effects
would inherit the gap — which is §12.1's "twelve right numbers under eleven wrong labels" in its other
form: right numbers, right labels, missing referent.

It also restores §7.2's own first rule. *"Every declared thing is rendered whether or not the data fills
it"* is why `center_USZ` gets a `dropped: constant` row despite no patient being at USZ. The reference
level is the most declared thing in the design — it is named in `REFERENCE_LEVELS`, asserted by
`test_config.py`, and load-bearing for every other coefficient — and it was the one declared thing the
table left out.

```
  column                  term    factor      level        status
  ─────────────────────────────────────────────────────────────────────────
  age                     linear  —           —            kept
  …
  onset_type_witnessed    dummy   onset_type  witnessed    reference (baseline)
  onset_type_unwitnessed  dummy   onset_type  unwitnessed  kept
  onset_type_wake_up      dummy   onset_type  wake_up      kept
  center_HUG              dummy   center      HUG          reference (baseline)
  center_CHUV             dummy   center      CHUV         kept
  center_Lugano           dummy   center      Lugano       kept
  center_USZ              dummy   center      USZ          dropped: constant
```

`n` is unchanged: it remains the **parameter** count including the intercept (12 on v7), which is what the
[§13] degrees-of-freedom budget is computed against. A reference row is a level that contributes no
parameter, and the table now says so rather than staying silent about it.

**The helper takes the covariate list, and that is not incidental.** `_design_table(X, dropped, covariates)`
— three arguments, where §8's first draft passed two. From `X.columns` and `dropped` alone the reference
rows can only be *inferred*, by looking for `{factor}_{level}` name prefixes among the columns that
survived and guessing which declared level is missing from the set. That happens to work on v7 and it is
the wrong shape: a linear covariate that happened to be named `center_something` would be misread as a
dummy, and the reference level is by construction in neither list, so `REFERENCE_LEVELS` has to be
consulted regardless. With `covariates` in hand the helper ranges over
`[c for c in covariates if c in C.CATEGORICAL]`, reads `FACTOR_LEVELS[c]` and `REFERENCE_LEVELS[c]`, and
emits every level with a `status` — which is a lookup rather than an inference, and is the same
declared-not-observed rule §4.2 applies to the matrix itself.

```python
def _design_table(X: pd.DataFrame, dropped: tuple[str, ...],
                  covariates: Sequence[str]) -> tuple[tuple[str, ...], ...]:
    """One row per level of every declared factor, plus one per linear covariate. §7.5.

    Ranges over the DECLARATION, never over the matrix: a level absent from both X and `dropped` is a
    reference level, and a table built from `X.columns + dropped` would omit exactly those.
    """
    header = ("column", "term", "factor", "level", "status")
    rows: list[tuple[str, ...]] = []
    for c in covariates:
        if c not in C.CATEGORICAL:
            rows.append((c, "linear", "—", "—",
                         "dropped: constant" if c in dropped else "kept"))
            continue
        for level in C.FACTOR_LEVELS[c]:
            name = f"{c}_{level}"
            status = ("reference (baseline)" if level == C.REFERENCE_LEVELS[c]
                      else "dropped: constant" if name in dropped
                      else "kept")
            rows.append((name, "dummy", c, level, status))
    return (header, *rows)
```

Note the row order: linear covariates and factors appear in **`covariates`** order, which is *not* the
design matrix's column order (§12.1 — `get_dummies` appends the dummies). That is deliberate. The matrix's
order exists to map a coefficient to a name and is asserted in §12.1; this table's order exists to be read
by a human against `PS_COVARIATES`, and §7.2's rule is that every declared thing is rendered, in the order
it was declared. `propensity_fit`'s coefficient table keeps `Fit.columns`' order, so the two tables are
ordered by their two different purposes and §12.11 asserts both.

### 7.6 What "describe the weighted population" requires, and the table that discharges it

[§7] asks for two things and this stage owed only one of them:

> Report the effective sample size per arm and **describe the weighted population**. Because the ATO
> population is defined statistically, **state explicitly who it comprises.**

The ESS is discharged by §7.2's first table. The **description** needs its own, and a sentence saying the
population *depends on the propensity specification* does not supply it: that is a caveat **about** the
population, not a description **of** it.

The gap matters because of what overlap weighting does. The analysed population is not "the 92": it is a
tilt of them, up-weighting patients near equipoise and down-weighting those whose treatment was nearly
determined. A reader asking "who is this estimate about?" cannot answer it from an ESS.

So `overlap_weights` carries a second table — one row per [§6] covariate, four columns:

```
  covariate        unweighted (bridging / EVT alone)     ATO-weighted (bridging / EVT alone)
  ────────────────────────────────────────────────────────────────────────────────────────────
  age                        …            …                       …            …
  nihss_baseline             …            …                       …            …
  core_ml                    …            …                       …            …
  …
```

Every cell is computed from `e`, `w`, `in_model` and the frame, all of which `fit` holds at the point the
entry is written. Continuous covariates carry the weighted mean; the factors carry the weighted
proportion in each level, one row per level, so `center` and `onset_type` are described rather than
summarised.

```python
def _weighted_population_table(df: pd.DataFrame, in_model: pd.Series, w: pd.Series,
                               covariates: Sequence[str]) -> tuple[tuple[str, ...], ...]:
    """One row per [§6] covariate — per LEVEL for a factor — unweighted and ATO-weighted, by arm. §7.6.

    Ranges over the declaration for the same reason _design_table does: a level nobody in the weighted
    population has is a fact about the population, and a row that vanishes is a fact nobody sees.
    """
    header = ("covariate", *(f"{lab} ({kind})"
                             for kind in ("unweighted", "ATO-weighted")
                             for lab in C.TREATMENT_LABELS.values()))
    sub, ww = df.loc[in_model], w.loc[in_model].to_numpy(dtype=float)
    arms = {code: (sub[C.TREATMENT] == code).to_numpy() for code in C.TREATMENT_LABELS}

    def cells(v: np.ndarray) -> tuple[str, ...]:
        return (*(_fmt(v[a].mean()) for a in arms.values()),
                *(_fmt(np.average(v[a], weights=ww[a])) for a in arms.values()))

    rows: list[tuple[str, ...]] = []
    for c in covariates:
        if c in C.CATEGORICAL:
            for level in C.FACTOR_LEVELS[c]:
                rows.append((f"{c} = {level}",
                             *cells((sub[c] == level).to_numpy(dtype=float))))
        else:
            rows.append((c, *cells(sub[c].to_numpy(dtype=float))))
    return (header, *rows)
```

A continuous covariate contributes its mean, a factor level its proportion, and both go through `_fmt`
so §12.11's byte-identity criterion holds. `np.average(..., weights=ww[a])` is the weighted mean **within
an arm**, which is what "who does this arm comprise after weighting" asks; it is not a between-arm
contrast, and §14 keeps Stage 6 out of judging one.

**Three boundaries, and they are what keep this from becoming Stage 7's job.**

- **These are not standardised mean differences.** No pooled standard deviation, no threshold, no verdict.
  [§9]'s |SMD| < 0.10 balance assessment is Stage 7's and §14 says Stage 6 does not judge the balance it
  reports. This table describes a population; Stage 7 evaluates an adjustment.
- **The unweighted column is the reconciliation.** §12.11 asserts it against the cohort's own means, so
  the table cannot drift from the frame it describes — Stage 5 §7.3's rule that every cell is computed
  from the object it describes, never by subtracting another row.
- **Nothing here is a negative control.** `BALANCE_ONLY`'s covariates are deliberately outside the
  propensity model and belong to Stage 7's balance table [§9]; this table ranges over `PS_COVARIATES`,
  which is what the weighted population is defined by.

## 8. `fit`, written out

Canonical, as Stage 5 §8 is. A reader traces the complete-case mask, the design, the fit, the weights
and all four audit entries through one screen.

```python
def fit(df: pd.DataFrame, audit: Audit) -> Propensity:
    """The [§7] propensity score and overlap weights over the [§3] cohort.

    Takes no covariate list: the propensity specification is PS_COVARIATES and the estimand is
    indexed by it [§7, §6.4]. Returns a Propensity; adds no column to `df` and edits nothing.
    """
    _assert_fit_inputs(df)                                     # F1, F2

    in_model = model.complete_cases(df, C.PS_COVARIATES)
    _record_exclusion(df, in_model, audit)                     # §7.2, entry 1

    X, dropped = model.design(df.loc[in_model], C.PS_COVARIATES)
    audit.record("model", "design_matrix", X.shape[1] + 1,
                 _design_detail(X, dropped),
                 table=_design_table(X, dropped, C.PS_COVARIATES))

    a = df.loc[in_model, C.TREATMENT].to_numpy(dtype=float)
    fitted = model.firth(X, a)                                 # raises FitError; no fallback [§5.5]
    _assert_probabilities(fitted, df.loc[in_model, "case_id"])  # F5

    e = pd.Series(np.nan, index=df.index, dtype="float64")   # np.nan, NOT pd.NA — §3.1
    e.loc[in_model] = fitted.p
    w = pd.Series(np.nan, index=df.index, dtype="float64")
    w.loc[in_model] = np.where(a == 1.0, 1.0 - fitted.p, fitted.p)

    audit.record("model", "propensity_fit", int(in_model.sum()),
                 _fit_detail(fitted, e), table=_fit_table(fitted))

    per_arm = {code: ess(w[in_model & (df[C.TREATMENT] == code)]) for code in C.TREATMENT_LABELS}
    audit.record("model", "overlap_weights", int(in_model.sum()), _weights_detail(),
                 table=_weights_table(df, in_model, w, per_arm)
                 + _weighted_population_table(df, in_model, w, C.PS_COVARIATES)[1:])

    return Propensity(e=e, w=w, in_model=in_model, ess=per_arm, fit=fitted, dropped=dropped)
```

Six things about the order, each of which is a failure if moved:

- **`_assert_fit_inputs` comes first**, before the mask and before the design — Stage 5 §4.4's rule.
  F1 is the exposure check and F2 is both-arms-present; a frame failing either produces a fit that
  raises later with a message about rank.
- **`_record_exclusion` comes before the design**, so the log names the excluded patients even if
  `design` then raises on D2 or D3. Stage 5's C1/C2 argument in the other direction: there the
  assertion had to precede the log; here the exclusion is a *fact about the frame*, established before
  any modelling choice, and a log that names it is useful precisely when what follows fails.
- **`e` and `w` are built as all-`nan` `float64` Series and filled**, never by `reindex` on a shorter
  Series, and the sentinel is `np.nan` and **not** `pd.NA`. Measured (§3.1, §18):
  `pd.Series(pd.NA, index=idx, dtype="float64")` raises `TypeError` — a `float64` Series has one
  missing value and `pd.NA` is not it. The masked `Float64` dtype would take `pd.NA`, and is declined:
  every downstream `.to_numpy(dtype=float)` would convert it back to `nan` silently, which is §3.1's
  first hazard reintroduced on the one column most likely to be summed. So the absence carried in the
  data is `nan`, the **deliberateness** of the absence is carried by `in_model` — boolean, total, never
  missing — and every downstream sum ranges over `in_model` rather than over `notna()`. Wherever this
  document says `<NA>` of `e` or `w`, it means `nan` with `in_model` False beside it.
  `reindex` is still declined: it would produce the same `nan` by accident rather than by construction,
  and on a caller's shorter Series it would also silently accept a mismatched index.
- **`a` is bound once, from the complete-case subset**, and used for both the fit and the weights. Two
  separate `df.loc[...][TREATMENT]` reads would be two chances to align the weights to a different
  subset than the fit.
- **`per_arm` ranges over `TREATMENT_LABELS`, not over the arms observed.** An arm absent from the
  cohort would be a `FitError` from `ess`, which is right: F2 has already established both arms are
  present, so an empty one here means the mask lost one.
- **`_assert_probabilities` sits between the fit and the weights**, because a degenerate `e` becomes
  a weight of exactly 0.0 one line later and stops being visible (§5.5, §3.1).

```python
def _assert_fit_inputs(df: pd.DataFrame) -> None:
    bad: list[str] = []

    off_arm = df.loc[~df[C.TREATMENT].isin((0, 1)), "case_id"]
    if len(off_arm):
        bad.append(
            f"F1  {C.TREATMENT}: {len(off_arm)} record(s) whose exposure is missing or outside "
            f"{{0, 1}}: {', '.join(sorted(off_arm))}. The response of the propensity model is the "
            "exposure; a missing one is a row with no y, and Int64's pd.NA becomes nan silently "
            "[§3.1]. Stage 4's E5 and Stage 5's C3 check the same property upstream.")

    present = [code for code in C.TREATMENT_LABELS if int((df[C.TREATMENT] == code).sum())]
    if len(present) < 2:
        bad.append(
            "F2  " + ", ".join(f"{C.TREATMENT_LABELS[c]}: "
                               f"{int((df[C.TREATMENT] == c).sum())}" for c in C.TREATMENT_LABELS)
            + ". A propensity model needs both arms: with one, the fit is a model of a constant and "
            "every overlap weight is 0 or 1. Stage 5's P2 guarantees both arms at every retained "
            "centre, so this can only be a caller subsetting after build().")

    if bad:
        raise C.SchemaError(
            "\n  " + "\n  ".join(bad)
            + f"\n  {len(bad)} propensity assertion(s) failed against {len(df)} records.")
```

```python
def _assert_probabilities(fitted: model.Fit, case_ids: pd.Series) -> None:
    """F5 — the fitted probabilities are finite and strictly inside (0, 1).

    Not a clip. `pilots/analysis.py:82` clips into [1e-8, 1-1e-8]; a probability of exactly 1.0
    means the linear predictor hit FIRTH_ETA_CLIP [§3.1], which under a penalised likelihood means the
    fit is degenerate, and a clip turns that into a weight of 0.0 that every sum accepts [§5.5].
    """
    bad = ~np.isfinite(fitted.p) | (fitted.p <= 0.0) | (fitted.p >= 1.0)
    if bad.any():
        raise model.FitError(
            f"F5  {int(bad.sum())} fitted probability/ies are not strictly in (0, 1): "
            f"{', '.join(sorted(case_ids.to_numpy()[bad]))}. Firth's penalty exists to keep the fit "
            "off the boundary, so a boundary value is a degenerate fit and not a confident one. "
            "[§10] drops and counts the replicate; it does not clip it into range.")
```

```python
def _record_exclusion(df: pd.DataFrame, in_model: pd.Series, audit: Audit) -> None:
    """The complete-case entry, which must name as many patients as it says it excluded.

    This is the rule `_MUST_NAME_CASES` cannot express for this kind [§7.1], stated here where it can
    be stated exactly — and it is Stage 5's `_record_removal` (cohort.py:200-206) line for line, both
    mechanisms included, because both are what make the comparison able to fail at all:

      * the identifiers go through a `set`, so a duplicated case_id collapses to one name;
      * a missing case_id is dropped by `notna()` rather than stringified.

    Without either, `len(case_ids)` equals `n` by construction and the check is dead code. That is not
    hypothetical — `Audit.record` normalises with `tuple(sorted(str(c) for c in case_ids))`
    (data.py:255), which sorts but does not deduplicate, and `str(pd.NA)` is the four characters
    `<NA>`. An earlier draft of this function collected `sorted(df.loc[~in_model, "case_id"])` into a
    plain list and could never raise; §12.11's acceptance test could not have passed.

    The raise comes BEFORE `audit.record`, again as Stage 5 does it. Recording first and reading the
    entry back leaves an unreconcilable entry in the Audit when the check fires, and `Audit` has no
    removal path — the log would carry the very inconsistency the raise exists to prevent.
    """
    excluded = df.loc[~in_model & df["case_id"].notna(), "case_id"]
    case_ids = tuple(sorted(set(excluded)))
    n_excluded = int((~in_model).sum())
    if len(case_ids) != n_excluded:
        raise C.SchemaError(
            f"covariate_completeness excludes {n_excluded} record(s) and can name {len(case_ids)}. "
            "An estimate whose denominator the log cannot reconstruct is an estimate nobody can "
            "check [§11]. A duplicated or missing case_id is forbidden by Stage 2's A2, so this is a "
            "belt over those braces — written here because this is the first stage whose denominator "
            "differs from the cohort's [Stage 6 §4.4].")
    audit.record("model", "covariate_completeness", n_excluded,
                 _completeness_detail(df, in_model), case_ids=case_ids,
                 table=_completeness_table(df, in_model))
```

## 9. Data flow into Stages 7-14

```
  fit(cohort, audit) returns Propensity
   ├─ e          float64, 93 rows, np.nan where in_model is False      [§4.4, §8]
   ├─ w          float64, same index, same missingness                 [§6.1]
   ├─ in_model   boolean, TOTAL, 92 True on v7                         [§4.4]
   ├─ ess        {0: 29.20, 1: 30.34} on v7, keyed as TREATMENT_LABELS [§6.2]
   ├─ fit        Fit(beta, p, iterations=7, converged_on="likelihood")  [§5]
   └─ dropped    ("center_USZ",) on v7                                 [§4.3]

  the cohort frame comes back unchanged — no column added, no value edited  [§0.2]

  audit — the same object load() created:
    load() 7 / derive() 4 / classify() 1 / build() 6 / fit() 4          [§12.11]
```

What each later stage may rely on, and what each owes:

- **Stage 7 [§9]** — reads `e`, `w` and `in_model`, and computes standardised mean differences against
  the **full** [§6] confounder set plus `BALANCE_ONLY`'s negative controls, regardless of what went
  into the score. Its within-centre table calls `propensity.ess` per centre per arm rather than
  reimplementing a Kish sum, and it reports the centres `treating_centres` does *not* return as
  structurally non-positive rather than plotting them. It is Stage 7, not Stage 6, that owns §6.3's
  near-balance finding.
- **Stage 8 [§8]** — weights the proportional-odds likelihood by `w` over `in_model`. It does **not**
  call `model.firth` (a proportional-odds fit is not a logistic one) but it **does** call
  `model.design` if it needs a matrix, so invariant 4 keeps holding.
- **Stage 9 [§8]** — calls `model.design(df, C.outcome_model_covariates(name))` and `model.firth` for
  `m_a(X)`, and **never** `propensity.fit`. The covariate list comes from Stage 1's registry so that
  invariant 3 is a fact about one configuration entry; the TICI override reaches the reporting layer
  because it comes back out of that call. Its own complete-case set is per outcome and is *not* this
  stage's `in_model` — an outcome with missing values excludes more patients, and [§11] wants both
  denominators.
- **Stage 10 [§10]** — refits **both** nuisance models in every replicate, by calling `fit` and
  Stage 9's fitter on the resampled frame. It must not reuse `e` or `w` from the point fit; §0.2 is
  why they are on a return value rather than on the frame. `FitError` is what it catches to drop and
  count, and it may catch nothing else — a `SchemaError` from F1/F2 or D1-D4 is a bug in the resampler,
  not a sparse replicate.

  **And it owns the local-maximum check, which `firth` cannot do (§5.2c).** The penalised likelihood is
  not concave, so `firth`'s fit from `beta = 0` is a stationary point and not certainly the maximiser
  [§7] prescribes. Measured: 1 replicate in 2000 on this cohort, with every guard clean. Stage 10 is the
  only place this can be caught, because catching it needs probes that **vary per replicate** and Stage
  10 is the stage that has both the loop and `C.SEED`. What it owes:
  - refit each replicate from **k additional random starts**, drawn from a generator seeded off `C.SEED`
    and advanced per replicate, and compare the penalised likelihood at each against the primary fit;
  - **count** the replicates where a start beats the primary by more than a declared tolerance, and
    report that count beside the intervals — as it already reports the `FitError` drop count;
  - decide, as a [§13] matter and not a coding one, whether such a replicate is dropped, corrected to
    the better fit, or kept and counted. Stage 6 does not decide it (§14), because the right answer
    depends on the rate Stage 10 measures and this stage can only bound it at ~1 in 2000.

  Two warnings from §5.2c's four failed detectors. Random probes must be **redrawn per replicate** — a
  fixed set seeded once caught 0 of 1. And the probe generator must be **separate from the resampling
  generator**: an early measurement here consumed draws from one stream for both, which silently changed
  which replicates were drawn and made two sweeps incomparable.
- **Stages 12 and 13 [§14]** — call `model.design` over `STANDARDISATION_COVARIATES` and fit **no
  propensity model at all**, which [§14] states twice. Nothing in `propensity.py` is on their path,
  which is what §0.1's split buys them.
- **Stage 14 [§16]** — reports the ESS, the weighted-population description and [§6.4]'s statement
  from §7.3's entry, and the `dropped` tuple beside the specification.

**One rule every one of them owes, stated once here because no assertion in this stage can enforce it.**
`e` and `w` carry `nan` off `in_model`, which is the safe direction — a stage that forgets the mask gets
`nan` rather than a plausible wrong number, and `nan` propagates through every sum until someone looks.
That safety holds **only while nothing fills it.** A `fillna(0)` anywhere downstream converts "never
weighed" into "weighed and found irrelevant" (§4.4), silently, in the direction that changes an estimate
rather than breaking it — and it is the natural thing to reach for when a weighted mean comes back `nan`.
So: **range over `in_model`, never over `notna()`, and never fill.** Stage 6 returns the mask precisely
so that no later stage has to infer it, and this sentence is the whole of what §19's third open question
could be answered with.

**And one for `model.py`'s direct callers specifically.** Stages 8, 9 and 12 call `design` and `firth`
without passing through `propensity.py`, so they get D1-D4, F3, F4 and F6 and **nothing else**. In
particular they do not get F1's response check or F2's both-arms check — those are the [§7] exposure's,
and an outcome model's response has different requirements. What F6 guarantees them is that a response
reaching `firth` is finite and dichotomous; what it cannot guarantee is that it is the *right*
dichotomy, or that its complete-case set matches the design's. Each of those stages owes its own mask and
its own denominator [§11].

## 10. What Stage 6 amends in Stages 1-5

The full ledger, so that no amendment is discovered during implementation. Two files are new pairs;
of the **seven** amended, four are the literal-pin sweep that a `KINDS` member *forces* — the same cost
Stage 5 §10 priced and for the same reason: a pin that disagrees with `KINDS` is a red suite, and
pinning against literals is what those tests are for. The seventh, `test_config.py`, is not forced by
`KINDS`: it is one assertion this stage's §4.2 discovered it depends on.

| File | Amendment | Why |
|---|---|---|
| `config.py` | the seven-constant Firth block of §5.6, named **`FIRTH_*`** | Stage 1's rule: a prespecified analytical choice lives in `config.py`. [§10] refits in every replicate, so a tolerance is part of the estimator. The prefix is `FIRTH_` and not `PS_` because `model.py` reads them and `model.py` does not know what the exposure is (§0.1, §12.12) — a one-file change now, a four-stage change once Stage 9, 10 and 12 have landed against the old names |
| `data.py` | `KINDS` gains `"model"` between `"cohort"` and `"structural"`; `_HEADINGS` gains `"model": "Fitted models"` in the same position | §7.1. The second `data.py` amendment since Stage 2, and the last one this pipeline needs — Stages 7-13 all render under `model` or under an existing kind |
| `data.py` | `_MUST_NAME_CASES` — **no change to the frozenset**; the comment above it gains one sentence and its "the other **six** kinds legitimately name none" becomes **seven** | §7.1. Stage 5 §6.3 established the posture; a count that is off by one in the same comment is how a reader learns not to trust it |
| `test_data.py` | `test_the_headings_are_the_declared_eight_in_pipeline_order` — both literals gain the kind and heading, and the test is **renamed** (`_eight_` → `_nine_`) | it pins `KINDS` and `list(_HEADINGS.values())` against literals |
| `test_data.py` | the comment inside that test gains §7.1's position argument | it currently explains why `derivation` and `cohort` sit where they do |
| `test_eligibility.py` | the `KINDS` literal in the renamed kinds-pin test | Stage 5 renamed it; the literal moves again |
| `test_cohort.py` | §12.9's `KINDS`-as-**eight** assertion becomes **nine**, and the test is **renamed**: `test_data_kinds_is_the_declared_eight_in_order` → `_nine_` (`test_cohort.py:780`) | Stage 5's own coverage map records it as a literal pin. This is the second of the two test names §1 counts as carrying the word *eight* |
| `test_config.py` | its existing `@pytest.mark.parametrize("factor", list(REFERENCE_LEVELS))` test gains one line: `assert C.REFERENCE_LEVELS[f] == C.FACTOR_LEVELS[f][0]` | §4.2. The by-name reference drop is equivalent to `drop_first=True` **only** while the reference is the leading declared level, and `test_config.py:592` asserts only that it is *in* the level set. Without this line the equivalence is a coincidence that a reorder of `CENTER_ORDER` breaks silently; with it, `design`'s choice is provably the model it claims to be |
| `pyproject.toml` | `numpy` and `scipy` move from "installed anyway" to declared direct use; the `reference` dev group is added | §2, §16b |
| `implementation_roadmap.md` | Stage 6 gains its `**Spec:**` line and three **Accept when** items | §19; lands with this document |

**Four things that look like they need amending and do not.**

- **`cohort.py` and `eligibility.py`.** Stage 6 reads their output and calls neither. `treating_centres`
  is Stage 7's and Stage 12's caller, not this stage's — Stage 6 has no per-centre logic at all.
- **`derive.py`.** `constant_covariates` **detects and logs** at Stage 5; `design` **drops** at Stage 6.
  Stage 3 §7.1 assigned the two to different stages in those words, so neither changes. The two agree
  on v7 — both find nothing — and their disagreement is the interesting case: a covariate constant in a
  bootstrap replicate is dropped by `design` and was never seen by `constant_covariates`, which runs
  once.
- **`config.py`'s `PS_COVARIATES`, `FACTOR_LEVELS`, `REFERENCE_LEVELS`, `CATEGORICAL`.** All four
  exist, all four were declared for this stage — `config.py:299-316`'s comment is addressed to it by
  name — and all four are read rather than written.
- **`test_data.py`'s two kind-parametrised tests.** They range over `data.KINDS`, so `model` is covered
  the moment it is declared. Their `else` branch accepts an entry with `n > 0` and no identifiers,
  which is three of §7.2's four entries. Read, not assumed (§18) — and T1's verification step is where
  the claim is checked rather than believed.

## 11. Handover to Stage 7

```
  cohort = cohort.build(eligibility.classify(derive.derive(*data.load()), ...), ...)
  ps     = propensity.fit(cohort, audit)

    92 of 93 weighted                       in_model               [§4.4]
    e in (0.02697, 0.93107), strictly interior  F5                 [§5.5]
    11 design columns + intercept            center_USZ dropped    [§4.3]
    3.25 treated patients per parameter      39 / 12               [§13]
    ESS 30.34 treated, 29.20 control         Kish, per arm         [§6.2]
    max weight 0.965, no trimming                                  [§6.1]
    converged in 7 iterations on the likelihood                    [§5.3]
```

Stage 7 receives the cohort and this object. It does not refit, does not re-weight, and does not
compute a second effective sample size.

## 12. Acceptance criteria

Tests live in `test_model.py` (§12.1-12.7, §12.12) and `test_propensity.py` (§12.8-12.14), with
section banners matching these numbers, as `test_derive.py`, `test_eligibility.py` and
`test_cohort.py` do. Tests needing the private workbook reuse the `DATA_GATED` idiom; everything else
runs on a plain checkout with no `data/`.

### 12.0 The frames and fixtures this stage is tested on

Three kinds of input, and the split is not arbitrary: the numerics are testable *only* on synthetic
data, and the integration is testable *only* on a real cohort.

**1. Synthetic numeric fixtures, in `test_model.py`.** Plain `numpy` arrays with known answers. These
carry the roadmap's three acceptance criteria and every property of §5, and they need no frame at all:

**Written out in full, not by signature.** A fixture given only a name and a signature makes every number
§18 pins from it unreproducible by a reader of this document, and contradicts the standing rule at the
head of this spec that nothing absent from it is to be invented. All three are specified here.

```python
def well_behaved(n=4000, seed=C.SEED):
    """§12.4 — Firth ≈ MLE at large n. Three standard-normal covariates, no near-separation."""
    rng = np.random.default_rng(seed)
    X = rng.normal(size=(n, 3))
    eta = 0.3 + X @ np.array([0.5, -0.4, 0.2])
    y = (rng.random(n) < 1.0 / (1.0 + np.exp(-eta))).astype(float)
    return pd.DataFrame(X, columns=["x1", "x2", "x3"]), y


def separated():
    """§12.5 — complete separation. The pilots' benchmark, unchanged."""
    return pd.DataFrame({"x": np.arange(1.0, 11.0)}), np.array([0] * 5 + [1] * 5, dtype=float)


def collinear(n=200, seed=C.SEED):
    """§12.6 — F3. well_behaved's design with x3 replaced by an exact copy of x1."""
    X, y = well_behaved(n, seed)
    X = X.copy()
    X["x3"] = X["x1"]
    return X, y
```

Both seed from `C.SEED`, never from a literal, so a fixture cannot drift from the run summary's
recorded seed. Measured (§18): `well_behaved()` gives `max |β_firth − β_mle| = 0.000787` in 4
iterations; `collinear()` is rank 3 of 4 at a condition number of 5.03e15, and `np.linalg.inv` on its
information matrix raises `LinAlgError`. `n=200` for `collinear` rather than 4000 because F3 raises
before the first iteration and a 4000-row design buys nothing but runtime.

**1b. The hard-fit family, in `test_model.py`.** Four of §5.2's numerical safeguards are unreachable on
every input above — measured, on the workbook design and on `separated()`, at zero trust-region
rescales, zero step-halvings and zero exhausted loops. These fixtures reach them, and each one is
recorded with the branch it exists to drive and the measurement that shows it does:

```python
def big_first_step():
    """§12.4a — the TRUST REGION. A tiny-scale separated design, so the first Newton step is huge.

    Measured under §5.2a's relative radius: ‖step₁‖ = 48.3 against FIRTH_MAX_STEP = 5.0, 2 rescales,
    converging in 9 iterations. The scale is what does it: the coefficient that separates these four
    points is ~1/0.05, so the unbounded Newton step at beta = 0 overshoots by an order of magnitude.

    This is also the fixture §12.4a multiplies by 1e-3 … 1e6 to assert scale invariance.
    """
    return pd.DataFrame({"x": [-0.05, -0.02, 0.02, 0.05]}), np.array([0, 0, 1, 1], dtype=float)


def needs_halving():
    """§12.4a — STEP-HALVING. Two orthogonal separating columns over four records.

    Measured: 1 halving, converging in 25 iterations, 0 trust-region rescales — so this fixture
    isolates halving from the trust region rather than firing both at once. Under
    FIRTH_MAX_HALVINGS = 1 it raises at iteration 2, which is the other half of §12.4a.
    """
    return (pd.DataFrame({"x1": [-1.0, -1.0, 1.0, 1.0], "x2": [-1.0, 1.0, -1.0, 1.0]}),
            np.array([0, 0, 1, 1], dtype=float))


def travels_far():
    """§12.4a — the fixture §5.2a is ABOUT. Scale pushed a further 20x down.

    Its Firth optimum is finite and perfectly identified: beta ~ (0, 2007.9), penalised
    log-likelihood -8.6292176, confirmed against scipy.optimize on the same objective.

    Under an ABSOLUTE trust radius it fails at FIRTH_MAX_ITER, because FIRTH_MAX_STEP * FIRTH_MAX_ITER =
    5.0 * 200 = 1000 caps total coefficient travel below ||beta|| ~ 2008 — a fit reported as a
    failure because of the covariate's units, not its data [§5.2a]. Under §5.2a's RELATIVE radius it
    converges at iteration 10, on the score route.

    It was called `travels_far` in an earlier draft, on the false account that "the linear
    predictor never stops growing". It does stop; the algorithm could not reach it. The rename is
    deliberate: a fixture named for a behaviour it does not have is a fixture that teaches the wrong
    lesson to whoever reads the test next.
    """
    return pd.DataFrame({"x": [-1e-3, -5e-4, 5e-4, 1e-3]}), np.array([0, 0, 1, 1], dtype=float)
```

**Two raises need a monkeypatched constant rather than a fixture, and that is the honest arrangement.**
`FIRTH_MAX_HALVINGS` exhaustion: `needs_halving()` under `monkeypatch.setattr(C, "FIRTH_MAX_HALVINGS", 1)`
raises at iteration 2. `FIRTH_MAX_ITER` exhaustion: `travels_far()` under
`monkeypatch.setattr(C, "FIRTH_MAX_ITER", 5)` raises with §5.2b's corrected message. In both cases the
constant *is* the subject of the test, so patching it is the test rather than a shortcut around one —
and after §5.2a no fixture reaches `FIRTH_MAX_ITER` by data alone, which is the point of §5.2a.

**One branch no fixture reaches, and it is recorded rather than faked** (§5.6, §3.1):
`_penalised_loglik`'s `-inf` return, which `FIRTH_WEIGHT_FLOOR` and F3 together make unreachable. §12.4a
asserts that unreachability — that the floor is never active — which is a real assertion that fails if
the floor changes, rather than a test of a branch that cannot run. **`converged_on == "score"` is NOT
in that category**: it fires on the small-scale fixtures, and §12.4a parametrises over both routes
(§5.3).

**2. `test_cohort.cohort_frame()`'s cohort — and it works, which was not obvious.** Stage 5's
five-record cohort has 9 [§6] covariates over 5 patients, which looks unfittable. Measured (§18): the
**constant-column rule strips nine of the twelve candidate columns**, leaving
`core_ml`, `tmax6_ml`, `center_CHUV` and an intercept — four parameters over five records, full rank,
and the fit converges in 11 iterations with `e ∈ (0.257, 0.929)`.

So `test_propensity.py` runs the whole pipeline end to end with no `data/`, which no stage before it
could do past Stage 5's P2. It imports `cohort_frame`, `hand_source` and `run` from the test modules
that declare them, never re-declaring them (Stage 3 §12's rule). **And the frame is the best available
test of §4.3**, because nine dropped columns is a bigger signal than the workbook's one.

**3. The workbook, data-gated.** §12.14 only.

**`TODOS.md`'s two open items that name this stage, resolved.** Both were checked rather than assumed
(§18b), and neither blocks:

- **"Commit a `tests/fixture_cohort.xlsx` that reads into a both-armed cohort"** names "Stages 6-13" as
  what will need a `data.Source` rather than a DataFrame. **Stage 6 does not.** Every test above runs
  either on a frame (`cohort_frame()`) or on synthetic arrays, and `test_propensity.py` takes the whole
  pipeline end to end with no `data/` through frame 2 above. The TODO's own instruction — build it when a
  stage first needs a `Source`, not speculatively — is therefore satisfied by leaving it alone. Stage 10
  is the likelier trigger: a resampling test wants a reproducible file rather than a frame built in
  Python.
- **"Declare an explicit dtype for the CTP volume columns"** names "Stage 6's design matrix" as a reason
  to do it. **It is not one.** Measured (§18b): `.astype(float)` maps `Int64`+`pd.NA`, `Float64`+`pd.NA`
  and `float64`+`nan` to `nan` identically, so `design` behaves the same under any of the three and F4
  catches all three. The TODO stands on its Stage 3 reasons; Stage 6 neither blocks it nor needs it.

### 12.0.1 The bare-versus-normalised rule still applies

Stage 5 §12.0.1's distinction carries forward unchanged and every test below must say which frame it
means. `cohort_frame()` returns raw centre codes; `run(...)` maps them to `CENTER_ORDER`'s labels.
A test meaning to reach `design` and passing a bare frame gets **D3**, not the check it intended — and
the failure looks like a passing `pytest.raises(SchemaError)` whose message is never read. Unless a
section says otherwise, every frame below is the normalised one.

### 12.0.2 The golden vector

Pinned in `test_propensity.py` as literals, from the synthetic frame so that no patient-derived number
enters git (§7.3). Measured (§18), on `cohort_frame()` normalised and taken through
`derive → classify → build → fit`:

```
  case_id     centre  arm   core_ml  tmax6_ml         e         w
  HAND-1      HUG      1      13.0     109.0   0.929038  0.070962
  HAND-2      CHUV     1      20.0      60.0   0.711530  0.288470
  HAND-5      HUG      1       0.0      88.0   0.702237  0.297763
  COHORT-1    HUG      0       8.0      50.0   0.257067  0.257067
  COHORT-3    CHUV     0       3.0      40.0   0.288470  0.288470

  design: core_ml, tmax6_ml, center_CHUV      9 of 12 candidates dropped as constant
  iterations 11, converged_on "likelihood"
  ESS 2.4413 of 3 treated, 1.9934 of 2 control
```

**One trap in that table, and it must be written into the test as a comment.** `HAND-2`'s weight and
`COHORT-3`'s weight and propensity are all `0.288470`, because the fixture's two CHUV records happen to
have linear predictors of equal magnitude and opposite sign. A test asserting weights by value alone
would pass while confusing the two records, so **the assertion is by identifier**, one record at a
time, and never against a sorted list of values.

The vector is asserted to 1e-6, which is looser than the 2.7e-15 agreement §16b measures against
`firthlogist`, and deliberately so: this is a regression pin against an edit to *our* algorithm, and
pinning it at machine precision would make it fail on a numpy patch release rather than on a mistake.

### 12.1 `design` reference-codes on the declared levels

- `center_USZ` is a column of the design over `cohort_frame()`'s cohort **before** the constant drop,
  and is in `dropped` after it — although no record is at USZ.
- The reference columns `center_HUG` and `onset_type_witnessed` are **absent** from the returned
  matrix, and neither appears in `dropped`: they were removed as references, not as constants, and the
  table of §7.2 distinguishes the two.
- **The column order is asserted against a literal, and the literal is not what it looks like.**
  `pd.get_dummies(sub, columns=factors)` **drops the original columns and appends the dummies at the
  end** — it does not expand a factor in place. So the order is: the non-factor covariates, in
  `covariates`' order; then each factor's dummies, in `factors` order (which is `covariates`' order
  restricted to `CATEGORICAL`), each expanded in `FACTOR_LEVELS` order with the reference removed.
  Measured (§18) over `PS_COVARIATES`, whose fifth entry is `onset_type`:

  ```
    age, sex, prestroke_mrs, nihss_baseline, core_ml, tmax6_ml, atrial_fib,
    onset_type_unwitnessed, onset_type_wake_up, center_CHUV, center_Lugano, center_USZ
  ```

  `onset_type` sits at index 4 of `PS_COVARIATES` and its dummies come after `atrial_fib`. An earlier
  draft of this bullet described the order as "`covariates`' order with each factor expanded in
  `FACTOR_LEVELS` order", which would have produced a literal that fails. This matters beyond the test:
  `beta`'s order **is** this order, `Fit.columns` is what maps a coefficient to a name, and §7.2's
  `propensity_fit` table is the only place the coefficients are ever written down (§7.3) — so a wrong
  order here is twelve correctly-computed numbers under eleven wrong labels.
- A frame in which every record is at one centre yields **no** centre column at all: two levels absent
  and the third constant. The design is still full rank, and this is the replicate case of §4.3.
- Called with `PS_COVARIATES_FULL`, the four extra names appear **before** the dummies, per the order
  rule above, because they are non-factor names. On the **workbook** all four survive (§12.14, data-gated).
  On `cohort_frame()`'s cohort all four are **constant and therefore dropped**, so the kept set is
  identical to `PS_COVARIATES`' — measured, and the assertion is written that way rather than as "four
  more columns appear", which is true only where §12.0 does not let this test run.
- Called with a covariate list containing **no factor at all**, `design` returns the linear columns
  unchanged and `dropped` reflects only genuine constants. Measured: `pd.get_dummies(sub, columns=[])`
  is a no-op and the `.drop(columns=[])` of references is empty, so no branch special-cases it. Asserted
  because Stage 9's `outcome_model_covariates` registry may yet hold an all-linear list, and a reader
  should not have to reason about whether the reference-drop comprehension is safe when `factors` is
  empty.

### 12.2 `drop_first` and `remove_unused_categories` are not used, and it is measured what they would do

- An AST scan over `model.py` asserts neither `drop_first` nor `remove_unused_categories` appears, with
  a companion test proving the scan fires on a pasted call.
- **`drop_first`'s behaviour is measured on both a declared and a bare column, because the two differ
  and only one of them is the hazard** (§4.2). On a frame with no HUG record:
  - a local reimplementation using `drop_first=True` on the **declared `Categorical`** is shown to drop
    `center_HUG` — the same column `design` drops by name — so the two forms agree and this half of the
    test asserts *equivalence*, not divergence;
  - the same reimplementation on the **bare `string`** column is shown to drop `center_CHUV`, leaving
    `center_Lugano` alone and every remaining coefficient meaning something else. This is the pilots'
    combination (`pilots/analysis.py:40` on a column `pilot_data.py` never declared);
  - `design` itself is shown **not** to raise `KeyError` on that frame — the reference column exists and
    is all-zero — and to raise **F3** from `firth` instead, at rank *k*−1 of *k*, because the surviving centre
    dummies are collinear with the intercept. A test asserting `KeyError` here was in an earlier draft
    and cannot pass.
  - the equivalence is asserted to depend on `REFERENCE_LEVELS[f] == FACTOR_LEVELS[f][0]`, by
    parametrising over both factors and pointing at `test_config.py`'s new assertion (§10) — so the
    reader of this test learns *why* by-name is kept when `drop_first` would do.
- A local reimplementation without the `Categorical` conversion is shown to emit no `center_USZ` column
  at all, so `dropped` is empty and the width is a property of the sample (§18).

### 12.3 `complete_cases` is returned, and a missing covariate cannot reach the fit

- On `cohort_frame()`'s cohort the mask is all-True; blanking `core_ml` on one record makes it False
  for exactly that record and leaves every other value in the frame unchanged.
- `design` does **not** drop the incomplete row: the returned matrix has as many rows as the frame,
  and the incomplete cell is `nan`. This is the assertion that `complete_cases` is a separate function
  rather than a `dropna` inside `design`.
- `firth` on that matrix raises `FitError` naming **F4** and naming the column — and the same matrix,
  with `_assert_fittable` removed, is shown to return all-`nan` coefficients and a plausible iteration
  count without raising (§3.1). The necessity of F4 is measured, not asserted in a comment.
- An `Int64` covariate with `pd.NA` reaches F4 by the same route as a `float64` `nan`, which is the
  §3.1 fact the test names. A `Float64` covariate with `pd.NA` reaches it too — measured, all three
  dtypes convert to `nan` identically under `.astype(float)`, so `design` is dtype-agnostic and the
  open CTP-volume dtype question in `TODOS.md` does not block this stage.
- **A `nan` in `y` raises `FitError` naming F6**, and the same response with F6 removed is shown to
  raise the *step-halving-exhausted* message at iteration 1 instead — so what F6 buys is measured to be
  the diagnosis rather than the raise. This is the response-side half of §3.1 and it is the branch
  Stage 9's outcome models will reach first (§5.4).

### 12.4 Firth matches the unpenalised MLE on well-behaved data [roadmap]

On `well_behaved(n=4000)`, `max |beta_firth − beta_mle| < 0.01`, against `statsmodels.Logit`.
Measured: **0.000787**, Firth converging in 4 iterations (§18), against the fixture as §12.0 now
specifies it. Firth's penalty is `O(1)` and the likelihood is `O(n)`, so the two agree as `n` grows; a
test at n = 100 would not distinguish a correct penalty from a mis-scaled one. The threshold stays at
`0.01` rather than tightening to the measured value: it is a statement about the *penalty's scaling*,
and pinning it at the measurement would make it fail on a numpy patch release rather than on a mistake
(§12.0.2's argument, one section on).

`converged_on` is asserted `"likelihood"` **for this fixture**; §5.3 records when the score route fires
instead, and §12.4a asserts the route per fixture rather than globally.

### 12.4a The four safeguards that only fire on hard data [§5.2, §5.6]

None of these branches is reached by the cohort, by `cohort_frame()`'s cohort, by `separated()` or by
`separated()` or `well_behaved()` — measured at zero rescales and zero halvings on the workbook and on
`separated()`. (**`cohort_frame()`'s cohort is not in that set**: measured 1 rescale and 2 halvings, which
§18d records for the golden vector. An earlier draft listed it among the four.) Each
one is what runs when a bootstrap replicate is sparse, which is the condition [§7] cites to justify this
estimator and the condition Stage 10 meets `N_BOOT` times. Untested, their first execution would be
inside that loop, where a defect reads as a dropped-replicate count rather than as a bug.

All four are asserted through `Fit`'s `first_step_norm`, `rescales` and `halvings` (§3.2). An earlier
draft asserted them on quantities `Fit` did not carry — `‖step₁‖`, the accepted step's norm, the halving
count — none of which is recoverable from a `(beta, p, iterations, converged_on, columns)` result, so
the four acceptance tests DoD-10b gates on could not have been written at all.

- **The trust region.** On `big_first_step()`, `fit.first_step_norm` is asserted `> C.FIRTH_MAX_STEP` and
  `fit.rescales >= 1`; the fit converges and `e` is strictly interior. Measured under §5.2a's relative
  radius: **9 iterations, 2 rescales**, `‖step₁‖ = 48.3`. **No companion asserts what removing the
  rescale does.** Measured: removing it makes this fixture converge in **8** iterations, not more, and
  moves the answer by 4.8e-06 — so the natural companion ("the trust region costs iterations and buys
  the same answer") is false in both halves. What the trust region buys is stated in §5.2a and measured
  there on the scale test, which is the assertion that matters. The `-inf` branch is **not**
  demonstrated here: an earlier draft claimed the unrescaled first step reaches `-inf` and burns
  halvings, and it does not — measured, the unrescaled candidate is finite and *uphill*
  (`ll = −4.88` against `ll_old = −6.04`) and is accepted immediately with zero halvings. §5.2's
  docstring claim about that branch stands on §3.1's argument, not on a companion test.
- **Step-halving.** On `needs_halving()`, `fit.halvings >= 1` and `fit.rescales == 0`, so halving is
  isolated from the trust region. Measured: 25 iterations, 1 halving, 0 rescales.
- **`FIRTH_MAX_ITER` exhausted**, under `monkeypatch.setattr(C, "FIRTH_MAX_ITER", 5)` on `travels_far()` —
  after §5.2a no fixture reaches it by data alone, which is §5.2a's point (§12.0 1b, §5.6). The message
  is asserted to carry
  a **non-zero** likelihood movement — §5.2b's correction, since the earlier form reported `0` on every
  input and the assertion could not fail — plus the largest score component and the rescale and halving
  counts. Under §5.2a's relative radius `travels_far()` **converges** (at iteration 10), so this
  assertion runs on a fixture chosen for it rather than on that one; §12.0 1b names which.
- **`FIRTH_MAX_HALVINGS` exhausted.** `needs_halving()` with `monkeypatch.setattr(C, "FIRTH_MAX_HALVINGS", 1)`
  raises at iteration 2, and the message is asserted to say there is no second estimator to try
  [invariant 5]. The constant is the subject of the test, so patching it is the test rather than a
  shortcut around one.
- **The trust region is scale-invariant**, which is §5.2a's whole point and is asserted directly: the
  same fixture with its covariate multiplied by 1e-3, 1, 1e3 and 1e6 **converges in every case** — which
  is the assertion. Under the absolute radius the 1e-3 case fails at `FIRTH_MAX_ITER` and the 1e-1 case takes
  84 iterations (§18d), so this is what catches a reversion to `if ‖step‖ > FIRTH_MAX_STEP`.

  The fitted probabilities agree across the four to **1e-3**, not to machine precision — measured
  `max|Δp| = 4.1e-04` between the 1e-3 and 1e6 scalings. That is not a defect and the tolerance is not
  slack: the optimum for the small-scale design sits at `‖β‖ ≈ 2008` on a surface flat enough that the
  1e-3 case converges on the **score** route and the others on the likelihood (§5.3), so the four stop at
  genuinely different points on the same ridge. Asserting 1e-8 here — as an earlier draft did — would be
  asserting that a flat surface has a sharp optimum.

**Then the convergence route, which is a live assertion and not the unreachability claim an earlier draft
made:**

- **`converged_on` is asserted per fixture, both values**, and the split is the assertion: `"likelihood"`
  on `separated()`, `well_behaved()`, `big_first_step()`, `needs_halving()` and the data-gated cohort
  design; **`"score"`** on `travels_far()` and on the ×1e-3 variants of `separated()` and
  `big_first_step()`. Parametrised, so both routes are exercised and neither can quietly stop firing.
  The test's comment carries §5.3's mechanism — `Δl* ≈ ½·s²·‖I⁻¹‖`, so an ill-conditioned design keeps
  the likelihood moving while the score goes flat — because a reader who sees `"score"` in a Stage 10
  replicate needs to know it means *near-singular design*, not *bug*.
- The measured `‖I⁻¹‖` for each fixture is asserted on the right side of roughly `1e4`, which is what
  makes the split above a property of the design rather than a coincidence of iteration counts.

**And one assertion of unreachability — the honest form of coverage for a branch no input reaches:**

- `FIRTH_WEIGHT_FLOOR` is asserted **never active**: over every fixture in `test_model.py` plus the
  data-gated cohort design, `min p(1−p)` is asserted `> C.FIRTH_WEIGHT_FLOOR`. Measured **0.0124** across
  `well_behaved`, `separated`, `big_first_step` and `needs_halving` — eight orders above the floor, so
  the assertion has real headroom and is not a near-miss dressed as a check. Fails if the floor is raised
  or a fixture starts saturating, which is exactly when someone needs to know. Assert against
  `C.FIRTH_WEIGHT_FLOOR` rather than against `1e-10`: the point is the relationship, not the number.
  This is also what keeps `_penalised_loglik`'s `-inf` branch unreachable (§3.1), so the two travel
  together and one assertion covers both.

### 12.5 Firth stays finite under complete separation, where the MLE diverges [roadmap]

On `separated()`:

- `firth` returns finite coefficients with slope in `(0.5, 2.0)` — measured **0.9706**, and
  `firthlogist`'s to four decimals (§18) — and `e` strictly interior. Also measured: intercept
  **−5.3385** in **9 iterations** on the likelihood route. An earlier draft added "which is `logistf`'s
  value for this benchmark"; that is withdrawn with §5.6's parallel claim, for the same reason —
  `logistf` has not been run, `firthlogist` is a port of it so agreement is transitive, and DoD-16 is
  where the attribution is either earned or struck. The threshold stays the pilots' `(0.5, 2.0)`:
  what the test is about is *finite versus divergent*, not four decimals.
- **`statsmodels.Logit` on the same data is asserted to fail**, and this half is what makes the test
  about separation rather than about agreement: Newton reports `converged=False` with a slope of
  **71.5**, and BFGS returns **27.2** and climbing (§18). Both are the divergence, and asserting only
  that Firth is finite would pass on data that was never separated.

### 12.6 A rank-deficient design raises, and `pinv` is not reachable

- `collinear()` raises `FitError` naming **F3**, its rank and its columns. Measured: rank 3 of 4 at a
  condition number of 5.03e15.
- **F6 raises on a response outside {0, 1}**, parametrised over the response a caller is most likely to
  hand this function by mistake: `mrs_90d`'s 0-6 scale. A companion with F6 removed is shown to *fit,
  converge and return coefficients* on that response, with no missing value and no warning — a logistic
  model of a seven-level ordinal, which is the one silent failure `model.py`'s outcome-agnostic boundary
  makes possible (§5.4). The test names [§8]'s proportional-odds model as where an ordinal belongs.
- An AST scan asserts `pinv` appears nowhere in `model.py`, with a companion proving the scan fires.
- With `_assert_fittable` removed, `np.linalg.inv` is shown to raise `LinAlgError` from inside the loop
  on the same input (§18) — so F3's job is the message, not the catch, and a `pinv` fallback would
  replace both.
- `design`'s D2 raises on `POST_TIME_ZERO` members: `onset_to_groin_min`, `onset_to_ivt_min`, and each
  of the eight outcomes, parametrised over `C.POST_TIME_ZERO` so a new outcome is covered by existing.
  **Invariant 4's test lives here.**
- D1 and D3 raise with their cases named, D3 on a frame carrying an unlisted centre — and a companion
  shows that without D3 the record is encoded as an all-zero dummy row, indistinguishable from the
  reference level, and the fit returns a number (§3.1).
- **D4 raises on a frame whose factor carries `pd.NA`**, and its companion is the sharper of the two:
  without D4, the missing record's dummy row is asserted **element-for-element equal** to a real
  reference-level record's, F4 is shown not to fire (no non-finite value anywhere in the design), F5 is
  shown not to fire (the fitted probability is finite and interior), and the fit returns a propensity
  for a patient the model has placed at the reference centre. Parametrised over `CATEGORICAL`, so a
  third factor is covered by the existing test.
- D3 and D4 are asserted to be **distinguishable**: a frame carrying one unlisted level and one missing
  value produces both messages in one `SchemaError`, naming different records. They share a loop and
  they are different failures.

### 12.7 No fallback estimator exists on any path [invariant 5]

- An AST scan over `model.py` and `propensity.py` asserts that no `import` names `sklearn`,
  `statsmodels` or `scipy`, and that no `except` clause in either module returns a value. Both scans
  have companions.
- `FitError` is asserted to be the only exception type either module raises for a numerical failure,
  and `SchemaError` the only one for a precondition — so Stage 10 can catch the first and not the
  second (§9).
- The repository is scanned for `firthlogist`: it appears in `test_model.py` and `pyproject.toml` and
  nowhere else.

### 12.8 `fit` returns aligned Series with the missingness reimposed

- `e`, `w` and `in_model` all carry the cohort's index, and the index is asserted **equal**, not merely
  the same length.
- With one covariate blanked, `e` and `w` are absent for exactly that record and `in_model` is False
  for it. **The other four probabilities are NOT asserted unchanged.** An earlier draft asserted they
  match §12.0.2's golden vector to 1e-12 "which is what proves the excluded record left the design
  rather than entering it as a zero" — that is false, and measured: blanking `HAND-1`'s `core_ml` moves
  the four survivors by up to **0.0478**. A logistic fit over four records is a different fit from one
  over five, so the change is the correct behaviour and the assertion would have failed. Worse, the
  natural repair — loosening the tolerance until it passes — tests nothing at all.

  What *is* asserted, and what actually distinguishes "left the design" from "entered as a zero":
  - the design handed to `firth` has **four rows, not five** — `X.shape[0] == int(in_model.sum())`;
  - the blanked record's `case_id` is in `covariate_completeness`'s named exclusions, and its `e` and
    `w` are `nan`;
  - refitting the four-record frame **directly**, with no blanking, gives probabilities identical to
    the blanked five-record run to 1e-12. This is the real form of the property: the excluded record
    left the design entirely, rather than contributing a zero row to it. A zero row would change the
    fit, and this assertion is what would catch it.
- `in_model` is asserted `bool` dtype with no missing value — total, per §4.4.
- The frame passed in is asserted unchanged, cell for cell, and `fit` is asserted to add no column
  (§0.2, §12.13).

### 12.9 Weights and ESS

- `w == 1 − e` on the treated rows and `w == e` on the control rows, exactly, and asserted through
  `TREATMENT_LABELS`' codes rather than literal 0/1.
- Every weight is in `[0, 1]`; nothing is trimmed; `IPTW_TRIM` does not exist in `config.py` and a test
  asserts the name appears nowhere in `extended_bridging/`.
- `ess([1] * 10) == 10.0` and `ess([1, 1, 1, 3]) == 3.0`, both exactly — the hand computations of
  §6.2.
- `ess` on an all-zero arm raises `FitError`, and a companion shows the pilots' `0.0` return would
  produce a finite ratio downstream.
- **`ess` on an EMPTY arm raises with the empty-arm message, not the all-zero one**, and the two messages
  are asserted distinct. They share a branch by arithmetic — measured, `np.sum(np.asarray([]) ** 2)` is
  `0.0` — and they are different diagnoses: all-zero means the fit put every patient on the boundary,
  empty means the mask lost an arm that F2 had already established was present. Stage 7 calls `ess` per
  centre per arm [§9], where the empty case is the one that actually occurs.
- `ess` drops `<NA>` before summing, so passing the full 93-row `w` and the 92-row masked `w` give the
  same number — asserted, because a Kish sum over a `nan` gives `nan` and it would be tempting to
  "fix" that with `fillna(0)`.

### 12.10 The weights do not balance exactly under Firth, and do under an MLE score

Two assertions, and the pair is the point (§6.3) — **but they run on different frames, and an earlier
draft ran both on `cohort_frame()`, where the second cannot hold.**

- **Firth, on `cohort_frame()`'s cohort:** `Σw` differs between arms — asserted **non-equal** beyond 1e-6.
  Runs with no `data/`.
- **The MLE companion is data-gated onto the workbook design**, because `cohort_frame()`'s cohort is
  5 records over 3 columns and **completely separated**, so there is no maximum-likelihood estimate to
  compare against. Measured: `statsmodels.Logit` Newton reports `converged=False` with `max|coef| = 130`
  and `Σw₁ − Σw₀ = −4.4e-11`, which "passes" a 1e-10 assertion only because both sums have underflowed;
  BFGS converges to a different point with `max|coef| = 168` and `Σw₁ − Σw₀ = 3.9e-07`, which fails it.
  The assertion's verdict was a function of the optimiser. On the workbook design the MLE converges
  (`max|coef| 4.55`) and the property holds properly: measured `Σw₁ − Σw₀ = 1.8e-15`, and every design
  column's weighted mean equal across arms to 1e-8.

So the near-balance is measured against its own explanation — on a frame where the explanation exists.
A silently reinstated MLE fails the first assertion, which needs no `data/`, so the discrimination §6.3
depends on survives on a plain checkout even though its companion does not.

### 12.11 The four audit entries, and the log

- The four `model` entries appear in `fit`'s declared order, at positions asserted against an index
  **captured before the call** — never a tail slice, which is the repair Stage 5 §10 landed in
  `test_derive.py` and which this stage must not reintroduce.
- `data.KINDS` is asserted to be the declared **nine**, and `_HEADINGS` to render `## Fitted models`
  between `## Cohort construction` and `## Structural non-applicability`.
- `covariate_completeness` names as many identifiers as its `n`, and `_record_exclusion`'s check is
  shown to fire on a frame with a duplicated `case_id` — modelled on `test_cohort.py:759`, which is the
  working precedent. The test also asserts that **no entry was recorded** when it fires, which is what
  the raise-before-record ordering buys (§8). A companion reverts `_record_exclusion` to the plain
  `sorted(...)` list of an earlier draft and shows the check cannot fail at all: `Audit.record`
  normalises with `tuple(sorted(str(c) for c in case_ids))` (`data.py:255`), which sorts without
  deduplicating, so `len(case_ids) == n` by construction. The necessity of the `set` and the `notna` is
  measured, not asserted in a comment.
- **`covariate_completeness`'s `n_absent` agrees with `absence_by_cohort_column`'s, for every [§6]
  covariate**, read out of the two entries' tables rather than recomputed. This is what replaces the
  per-centre columns §7.4 removed: the numbers stay in one place and a drift between the two tables is a
  red suite. The `model` entry's own contribution, `excluded_by`, is asserted to name the covariates that
  actually cost a patient their weight — `core_ml` and `tmax6_ml` on v7 — and nothing else.
- **`design_matrix`'s table has one row per LEVEL of every declared factor plus every linear covariate**
  — 14 rows on v7, not 12 — and the reference rows are asserted present and carrying
  `status = "reference (baseline)"`, distinct from `dropped: constant` (§7.5). `center_USZ` is asserted
  present as `dropped: constant`, and `center_HUG` and `onset_type_witnessed` as references. Its `n` is
  the parameter count *including* the intercept — 12 — which is the number the [§13] degrees-of-freedom
  budget is computed against and is **not** the row count. The two being different is the point: a level
  that contributes no parameter still gets a row.
- **`overlap_weights`' second table is asserted** (§7.6): one row per [§6] covariate, with the
  **unweighted** column reconciled against the cohort's own means and proportions computed directly from
  the frame — so the description of the weighted population cannot drift from the population it
  describes. Factor covariates are asserted to contribute one row per level. The ATO-weighted column is
  asserted to differ from the unweighted one on at least one covariate, because a weighted population
  identical to the unweighted one would mean the weights did nothing and the table would be decoration.
- `overlap_weights`' table renders a row for every `TREATMENT_LABELS` entry plus `all`, and its
  `detail` contains [§6.4]'s conditionality sentence verbatim.
- The log renders byte-identically twice in one process, and under `PYTHONHASHSEED=0` versus `1` in a
  subprocess — the two-seed driver of Stage 5 §12.9, extended one stage.
- No entry's `detail` or table contains a numpy repr: every float goes through `data._fmt`.

### 12.12 The module boundary holds

- An AST scan asserts `model.py` names neither `C.TREATMENT` nor `C.PS_COVARIATES` nor any weight
  concept: it is outcome-agnostic (§0.1), and Stage 9 must be able to import it without importing an
  exposure. Companion test included.
- **And that no `C.` attribute it reads begins with `PS_`.** The seven tolerances were named `PS_*` in an
  earlier draft — for the propensity score, inside the module defined by knowing nothing about it — and
  they are `FIRTH_*` now (§5.6). The attribute scan above would have passed the old names, because they
  are neither `TREATMENT` nor `PS_COVARIATES`; this bullet is what makes the boundary a rule rather than
  a habit, and it is what stops Stage 9 reading a propensity-flavoured constant for an outcome fit.
  `propensity.py` is exempt: `PS_COVARIATES` is exactly what it is *for* (§0.1).
- **And a second scan, over the string literals inside `model.py`'s `raise` statements**, asserting that
  no runtime message names the exposure or any specific outcome — the words `ivt`, `TREATMENT`,
  `propensity`, and every key of `C.OUTCOMES`, parametrised over `C.OUTCOMES` so a ninth outcome is
  covered by the existing test. Companion test proving the scan fires on a pasted message.

  **Two things about the scope, and both are the point.** It is not belt-and-braces over the first scan:
  an attribute scan passes a module whose *error text* is full of exposure names, and F6's first draft was
  exactly that — "use `propensity.fit`'s F1, or your outcome's own mask", and "`mrs_90d` has seven levels
  and would fit silently" (§5.4). A module that says that in a raise is a module the next reader treats as
  knowing about the exposure, and the concrete examples belong in this spec, where they are more useful.

  And it is scoped to `raise` arguments rather than to all strings **because the module docstring must
  name `propensity.py`** — §17 requires the sentence "Stages 8, 9 and 12 come into this module directly
  and never through `propensity.py`", which is the one sentence that stops someone moving `design` there
  for tidiness. A whole-file text scan would forbid the comment this stage most wants to keep. So: the
  docstring may name the module it must not depend on; a raise may not. An AST walk over `ast.Raise`
  nodes collecting `ast.Constant` strings in the exception's arguments is the implementation, and the
  companion pastes a forbidden word into a raise rather than into a docstring, so the test proves the
  scope as well as the scan.
- `model.py` is asserted to take no `Audit` parameter anywhere, and to import neither `data` nor
  `propensity`.
- `propensity.fit` is asserted to take exactly `(df, audit)` — no covariate parameter, no `scheme`, no
  centre list (§6.4, §16). A companion asserts the signature by `inspect`, so adding a defaulted
  keyword fails the test rather than passing silently.
- Neither module contains a bare `assert` (they vanish under `-O`), neither names a raw header, and
  neither is in `EXEMPT_FROM_RAW_NAME_SCAN` — with the scan-fires companion for each.

### 12.13 Stage 6 adds nothing to the frame and needs Stage 5

- `fit` on a frame that has not been through `build` — classified but unrestricted, **9 rows on
  `cohort_frame()`** (126 on the workbook, which is §12.14's business and not this section's; §12.0
  item 3 says the workbook appears in §12.14 only, and an earlier draft of this bullet quoted 126 here,
  contradicting it and putting a data-gated number in a test DoD-2 requires green without `data/`) — is
  **not** asserted to raise: it is a legitimate call, and Stage 12 makes something like it. What is
  asserted is that it raises **F2** on a frame with one arm, and that on the unrestricted frame the
  ineligible patients are *in* the fit, which is the property Stage 12 needs and roadmap invariant 2
  forbids for [§7]. The test names invariant 2 and records that the guarantee is Stage 5's, not this
  stage's.
- `fit` on a frame with no `core_above_median` works, because Stage 6 never reads a subgroup.
- The returned frame is the frame given: same columns, same values, same dtypes, same index — the four
  assertions of Stage 5 §12.14, one stage on.

### 12.14 Structural facts from the workbook [data-gated]

- The cohort is 93 records and `in_model` is **92**; the excluded record is at Lugano, in the control
  arm, and is missing `core_ml` and `tmax6_ml`.
- The design is **11 columns**, `dropped == ("center_USZ",)`, and the parameter count including the
  intercept is **12** — so the [§13] budget is 39 / 12 = 3.25 treated patients per parameter, asserted
  as a number so a covariate added without an amendment moves it.
- The fit converges in **7 iterations** on `"likelihood"`, and `e` is strictly interior with
  `min > 0.02` and `max < 0.94`.
- ESS is **30.34** treated of 39 and **29.20** control of 53, each to two decimals.
- The largest weight is below 1.0 and above 0.9 — measured 0.965 — and no weight is 0.
- `constant_covariates(cohort, PS_COVARIATES)` is empty while `dropped` is not, which is §10's
  detect-versus-drop distinction asserted rather than described.

### Coverage map

```
model.py                                           test_model.py
├── design(df, covariates)
│   ├── declared FACTOR_LEVELS → Categorical        ├── 12.1  center_USZ exists, then drops
│   ├── reference dropped BY NAME                   ├── 12.1, 12.2  drop_first EQUIVALENT on a
│   │      (absent level → F3, NOT KeyError)        │         declared Categorical, divergent on a
│   │                                               │         bare column; coupling pinned in
│   │                                               │         test_config.py (§10)
│   ├── column ORDER: linear then dummies           ├── 12.1  literal, get_dummies APPENDS
│   ├── constant columns dropped after dummying     ├── 12.1  9 of 12 on the hand cohort
│   ├── covariate list with no factor at all        ├── 12.1  get_dummies(columns=[]) is a no-op
│   ├── no remove_unused_categories, no drop_first  ├── 12.2  two AST scans + both fire
│   └── returns (X, dropped)                        └── 12.1  dropped excludes the references
│
├── _assert_design_inputs(df, covariates)
│   ├── D1 covariate absent      → SchemaError      ├── 12.6
│   ├── D2 POST_TIME_ZERO        → SchemaError      ├── 12.6  parametrised over the denylist
│   │                                               │         — INVARIANT 4 lives here
│   ├── D3 unlisted factor level → SchemaError      ├── 12.6  + what it would encode without D3
│   ├── D4 MISSING factor value  → SchemaError      ├── 12.6  + the all-zero row is identical to
│   │                                               │         a real reference record — §4.5
│   └── several at once → one message               └── 12.6
│
├── complete_cases(df, covariates)
│   ├── returned, never applied                     ├── 12.3  design keeps the nan row
│   └── total, boolean                              └── 12.8
│
├── firth(X, y)                                     written out in §5.2
│   ├── ≈ MLE at large n                            ├── 12.4  0.000787 at n = 4000, 4 iterations
│   ├── finite under separation, MLE diverges       ├── 12.5  0.9706 vs 71.5 / 27.2
│   ├── p finite and in (0, 1)                      ├── 12.5, 12.8
│   ├── step-halving on the PENALISED likelihood    ├── 12.4  local reimpl. on l(b) diverges
│   ├── converges on the LIKELIHOOD                 ├── 12.4, 12.4a  asserted on every fixture
│   ├── raises, never falls back                    ├── 12.7  three scans + companions
│   ├── _assert_fittable  F3 / F4 / F6              ├── 12.3, 12.6  all three seen without the guard
│   ├── trust region, RELATIVE to ‖beta‖            ├── 12.4a  big_first_step(), 2 rescales — §5.2a
│   ├── trust region is SCALE-INVARIANT             ├── 12.4a  same fixture x 1e-3 … 1e6
│   ├── step-halving actually halves                ├── 12.4a  needs_halving(), 1 halving
│   ├── FIRTH_MAX_ITER exhausted → FitError            ├── 12.4a  travels_far() + monkeypatch
│   ├── the message reports a NON-ZERO movement     ├── 12.4a  §5.2b — it was always 0 before
│   ├── FIRTH_MAX_HALVINGS exhausted → FitError        ├── 12.4a  needs_halving() + monkeypatch to 1
│   ├── FIRTH_WEIGHT_FLOOR clip                        ├── 12.4a  asserted NEVER ACTIVE — §5.6
│   ├── converged_on, BOTH routes                   ├── 12.4a  parametrised: "likelihood" on the
│   │                                               │         well-scaled fixtures, "score" on the
│   │                                               │         ill-conditioned ones — §5.3
│   └── _penalised_loglik → -inf                    └── 12.4a  UNREACHABLE behind F3 + the floor;
│                                                             NOT demonstrable by removing the
│                                                             rescale — §3.1, §5.2, §20 rd 3
│
└── no bare assert / no raw header / not exempt     └── 12.12 + the scans fire

propensity.py                                      test_propensity.py
├── _assert_fit_inputs(df)   F1 arm / F2 both arms  ├── 12.13
├── fit(df, audit)                                  written out in §8
│   ├── e, w aligned, np.nan off in_model           ├── 12.8  index equality, not length
│   │      (pd.NA RAISES on float64 — §3.1)         │
│   ├── w = 1−e treated, e control                  ├── 12.9  through TREATMENT_LABELS
│   ├── no trimming; IPTW_TRIM nowhere              ├── 12.9  repository scan
│   ├── F5 probabilities strictly interior          ├── 12.8  + the clip declined
│   ├── frame returned unchanged, no column added   ├── 12.13
│   ├── four `model` entries in order               ├── 12.11 position vs captured index
│   ├── byte-identical, two hash seeds              ├── 12.11
│   ├── takes NO covariate list                     ├── 12.12 signature asserted by inspect
│   └── the golden vector                           └── 12.0.2 by identifier, not by value
├── ess(w)
│   ├── hand-computable                             ├── 12.9  10.0 and 3.0, exactly
│   ├── raises on an all-zero arm                   ├── 12.9  + the pilots' 0.0 declined
│   ├── raises on an EMPTY arm, DIFFERENT message   ├── 12.9  same branch by arithmetic
│   └── drops <NA> before summing                   └── 12.9
└── _record_exclusion(...)  names as many as its n  └── 12.11 fires on a duplicated case_id, and
                                                             records NOTHING when it fires; the
                                                             plain-list draft shown unfirable

data.py / config.py amendments                     the existing suites
├── KINDS gains "model" in position                ├── 12.11 + three renamed literal pins
├── _HEADINGS gains "Fitted models"                ├── 12.11
├── _MUST_NAME_CASES unchanged                     ├── the kind-parametrised tests, unedited
├── the seven Firth constants                      ├── test_config.py's declaration tests
└── REFERENCE_LEVELS[f] == FACTOR_LEVELS[f][0]     └── test_config.py, §4.2's coupling (§10)

§16b oracles                                       test_model.py
├── firthlogist, dev-only, importorskip            ├── 1e-12 on the cohort design and on
│                                                   │       separated()
├── scipy.optimize on the SAME objective           ├── the Newton optimum IS the optimum
├── statsmodels.Logit at large n                   ├── 12.4
└── the golden vector, no dependency               └── 12.0.2

Data-gated (skipif not DATA_XLSX.exists()):        └── 12.14  92 of 93, 11 columns + intercept,
  own module-scoped `workbook` fixture                        7 iterations, ESS 30.34 / 29.20,
  (load → derive → classify → build)                          max w 0.965, 3.25 per parameter

Every branch above is either tested or listed as UNREACHABLE with the constant that makes it so.
ONE is unreachable — `_penalised_loglik`'s -inf, behind F3 and FIRTH_WEIGHT_FLOOR (§3.1) — and it is
covered by an assertion that it STAYS unreachable, which fails if the floor moves. Two earlier
drafts of this block were wrong in opposite directions: the first read "No branch is untested"
while five branches had no input reaching them; the second listed `converged_on == "score"` as
unreachable, which it is not (§5.3) — it fires on any design with an ill-conditioned information
matrix, and §12.4a now parametrises over both routes.

  model.py       design 13, complete_cases 4, firth 19           =  36  ->  35 tested, 1 unreachable
  propensity.py  fit 11, ess 6, _record_exclusion 1              =  18  ->  18 tested, 0 unreachable
  amendments     KINDS, _HEADINGS, _MUST_NAME_CASES, 7 constants, the coupling  =   5  ->   5 pinned
  ------------------------------------------------------------------------------------------------
  59 paths total.  58 with a test that can fail.  1 with an assertion that it cannot run.
  0 paths with no coverage of either kind. That is the claim; the map above is how to check it.

  If the map and this block disagree, the map is authoritative and this block is stale. The count is
  stated at all so that "every branch" can be checked against a number rather than believed.
```

### Failure modes

| Failure | Test | Error handling | Visible or silent |
|---|---|---|---|
| A missing covariate reaching the fit — all-`nan` coefficients, plausible iteration count, `[nan, nan]` intervals in Stage 10 | 12.3 | **F4** before the first iteration, naming the column; `complete_cases` upstream of it | visible **only** because of F4 — silent without it (§3.1) |
| A record at an unlisted factor level modelled as the reference | 12.6 | **D3**, naming the cases and the level it would have become | visible **only** because of D3 |
| A record with **no** factor value modelled as the reference — the all-zero dummy row is byte-identical to a real reference record | 12.6 | **D4**. F4 cannot see it (no non-finite value survives the encoding), D3 excludes it by design, F5 cannot see it (the probability is finite and interior) | visible **only** because of D4 |
| A post-time-zero variable entering a design matrix — invariant 4 | 12.6 | **D2**, over the explicit `POST_TIME_ZERO` denylist, in the one function every model's covariates pass through | visible |
| A bootstrap replicate or subgroup losing the **reference** centre entirely | 12.2, 12.6 | the all-zero reference column is dropped, the surviving centre dummies are collinear with the intercept, and **F3** raises at rank *k*−1 of *k*, one dependency. **Not** `KeyError` — the declared `Categorical` guarantees the column exists (§4.2, measured) | visible |
| A design built on a **bare** factor column, where `drop_first` really does rebaseline on the first present level | 12.2 | the `Categorical` conversion is unconditional in `design`; the bare form is exercised only as a local reimplementation in the test, and measured | visible in review |
| A missing **outcome** reaching an outcome model's fit in Stage 9 — every step `nan`, halvings exhausted, and a message about [§7]'s single estimator | 12.3 | **F6** on the response, naming the count and pointing at the caller's own complete-case mask | raises either way; **F6 is what makes the diagnosis right** |
| An **ordinal** response handed to `firth` — `mrs_90d`'s 0-6 fitted as a logistic model of something else | 12.6 | **F6**'s `{0, 1}` check, naming the offending values and pointing at [§8]'s proportional-odds model | visible **only** because of F6 — silent without it, with no missing value anywhere |
| An empty arm reported as an all-zero arm, sending a reader after a boundary fit that never happened | 12.9 | `ess` branches on `v.size` before the sum; the two messages are distinct and asserted so | visible |
| A `"score"` convergence read as a bug rather than as the signal it is | 12.4a | `converged_on` asserted **per fixture, both routes** — `"likelihood"` on the well-scaled ones, `"score"` on the small-scale ones — with §5.3's `Δl* ≈ ½·s²·‖I⁻¹‖` mechanism in the test's comment | visible |
| An absent factor level making the design exactly singular | 12.1, 12.6 | the constant-column rule removes it; **F3** if anything else is rank-deficient | visible |
| A `pinv` fallback returning the minimum-norm solution for a design that identifies none | 12.6, 12.7 | `inv` is not wrapped; **F3** raises first with the rank and the columns | visible |
| A second estimator on the path — the pilots' 16.6% fallback rate | 12.7 | `FitError` only; import and `except` scans over both modules | visible |
| A degenerate fit clipped into range, becoming a weight of exactly 0.0 | 12.8 | **F5** asserts strict interiority instead of clipping | visible |
| An arm whose weights all vanish, reported as ESS 0.0 and divided by downstream | 12.9 | `ess` raises `FitError` | visible |
| A propensity model fitted on a different covariate list by a caller passing an argument | 12.12 | `fit` has no such parameter; the signature is asserted | visible |
| The point fit's `e` reused inside a bootstrap replicate | — | `Propensity` is a return value, not a column (§0.2); Stage 10 owns the loop | **structural, not asserted here** |
| An estimate over 92 patients reported against a denominator of 93 | 12.8, 12.11 | `in_model`; `covariate_completeness` names the excluded; [§11] per-estimate denominators | visible |
| A fitted coefficient vector committed to git | 12.0.2 | the golden vector is on the synthetic frame; the workbook fit is asserted on properties (§7.3) | visible in review |
| Firth's near-balance mistaken for exact balance, or an MLE silently reinstated | 12.10 | `Σw` asserted **unequal** under Firth and equal under an MLE score | visible |
| A tolerance changed as a runtime knob rather than as an amendment | 12.12 | the seven constants are in `config.py`; nothing takes them as an argument | visible |

Every failure above is visible at the point it occurs, and **three** are visible only because of an
assertion this stage adds: a `nan` in the design (F4), a record at an unlisted factor level (D3), and an
ordinal response (F6). All three return a plausible number otherwise, and the third was found by the
engineering review rather than by writing the spec (§20). That is the difference between this stage and
the five before it: they failed loudly by default.

**One further asymmetry worth naming.** F4, D3 and F6 all guard `model.py`, which is the module with no
`Audit` and no exposure — the outcome-agnostic one Stages 8, 9 and 12 enter directly (§0.1). Every guard
this stage adds for its own sake turns out to be a guard those stages inherit, and F6 exists *only*
because of them: the propensity model's response is checked by F1 in `propensity.py` and would never
have needed it.

## 13. Known gaps carried forward

- **The ATO population is 92, not 93, and nothing bounds that gap on the next workbook.** [§11] is
  complete-case per estimate and the [§3] cohort is 93, so the two differ by whoever is missing a [§6]
  covariate. On v7 that is one control at Lugano. A workbook in which `core_ml` is missing for twenty
  patients would silently shrink the target population by twenty, with `covariate_completeness`
  reporting it and nothing raising. Stage 6 does not set a threshold, because a threshold is a [§11]
  amendment and the PI's. **§7.2's entry is what makes it noticeable; it is not what makes it safe.**
- **3.25 treated patients per parameter, against the ~3.5 Stage 0 confirmed.** Measured: 39 treated
  over 12 parameters, because `center` expands to two dummies over the cohort's three centres and
  `onset_type` to two over its three levels. The [§6] budget is not violated — it was always
  approximate — but the number is now exact and it is the denominator every [§13] sensitivity
  specification would move. The full-covariate model of [§13] would take it to 16 parameters and 2.44
  treated per parameter, which is the reason that row is a *sensitivity* analysis and not the primary.
- **The near-balance of §6.3 is unquantified here.** `Σw` differing by 0.48 between arms says the Firth
  score does not solve the exact-balance equation; it says nothing about how much residual imbalance
  that leaves on any covariate. That is Stage 7 [§9] and its |SMD| threshold, and the two must be read
  together: `pilots/test_estimators.py` records the worst |SMD| under a Firth score as
  `0 < worst < 0.05` on synthetic data, which is inside [§9]'s 0.10 but is not zero.
- **No replicate-level evidence yet.** Every claim in §18 is about the point fit. [§7]'s stated reason
  for Firth — that an unpenalised MLE fails in a minority of sparse replicates — is measured in this
  repository only for the *pilot's* cohort and covariate list (16.6%). Stage 10 will measure it for
  this one, and if the rate is zero, [§7]'s choice is still correct and its justification will want a
  sentence.
- **`converged_on` is recorded and nothing reads it — and it is a live signal, not a constant.** Stage 10
  should surface the distribution of the two routes across replicates: a shift toward `"score"` is the
  signature of near-separation becoming common, which is the condition under which the ATO estimand
  itself is fragile. §5.3 measures **when** it fires — on a design whose information matrix is
  ill-conditioned, `‖I⁻¹‖` of order `1e4` and above, because `Δl* ≈ ½·s²·‖I⁻¹‖` keeps the likelihood
  moving while the score goes flat. On v7 every route is `"likelihood"`; a replicate that comes back
  `"score"` has a near-singular design.

  **Two earlier drafts of this bullet were wrong and the second was worse than the first.** The first said
  nothing reads the field. The second said its second value was structurally **unreachable**, from an
  argument that assumed `‖I⁻¹‖ ~ 1` and from fourteen fixtures that were all well-scaled — and it used
  that to file "should the criteria be reordered" as the open question. There is no such question: the
  criteria work as specified. What Stage 10 owes is the distribution, and what this stage owes is the
  fixtures that let a reader recognise a `"score"` result when they see one (§12.4a).
- **The estimator can converge to a local maximum, and nothing in this stage detects it.** §5.2c is the
  full record: `½log|I(b)|` is not concave, so `firth` returns a stationary point where [§7] prescribes
  the maximiser. Measured at **1 replicate in 2000** on this cohort and **4 of 600** synthetic separated
  designs, with every guard clean — F5, F3, `converged_on`, `rescales` and `halvings` all silent. The v7
  point fit is at the global optimum, so no reported number is affected; the exposure is [§10]'s refits
  and Stage 11's subgroups.

  **This is the only unguarded failure mode in Stage 6, and it is unguarded by decision rather than by
  omission** (2026-08-14). Four detectors were designed and measured to fire 0 of 1 (§5.2c's table); the
  basins are narrow and design-dependent, so only probes that vary per fit find them, and `firth` holds
  no seed. §9 assigns the check to Stage 10, which has the loop and `C.SEED`. Two things are owed back
  here: the **rate** Stage 10 measures — the prevalence above rests on two observed events, so a Poisson
  interval spans roughly 0 to 6 per 2000 — and, if that rate is materially higher than 1 in 2000, a
  [§13] decision on whether such replicates are dropped, corrected or counted. At the measured rate the
  effect on a percentile bound is far below its own Monte Carlo error, which is why this is filed rather
  than blocking.
- **Four safeguards in the estimator have no natural input and are reached only by contrived fixtures.**
  §12.4a's `big_first_step`, `needs_halving` and `travels_far` exist because the trust region,
  step-halving and the iteration cap fire on none of the cohort, `cohort_frame()`'s cohort,
  `separated()` or `well_behaved()` — measured at zero rescales and zero halvings on all four. They are
  now tested, but tested on synthetic designs chosen to reach them, which is weaker evidence than a real
  sparse replicate. **Stage 10 is where that evidence arrives**, and it should record how often each
  fires across `N_BOOT` replicates. If the answer is "never", the trust region and the halving loop are
  dead weight on this dataset and that is worth knowing; if it is "often", §5.6's constants are load-
  bearing on the reported intervals and deserve a sensitivity row of their own.
- **The R cross-check against `logistf` and `PSweight` has not been run, and it is the strongest
  evidence available.** §16b.1 is the argument: R 4.5.0 is installed on this machine, the four
  relevant packages are not, and every oracle this stage does have is either transitive
  (`firthlogist` is a port of `logistf`), self-referential (the golden vector, the `scipy` check on our
  own objective) or partial (the MLE agreement, which bites only at large n). A direct `logistf`
  comparison is the one independent check on the estimator, and `PSweight` — written by the authors of
  the overlap-weight method — is the only available check on the **whole [§7] pipeline** rather than on
  its first step. **It is a blocker, by decision of 2026-08-14, and Definition of done 16 is that
  blocker.** An earlier draft of this bullet and of §16b.2 called it "an open item rather than a
  blocker" while DoD-16 required it before Stage 6 could be called complete; the document argued both
  sides and the contradiction is resolved in favour of the gate. **Stage 6 is not complete until it has
  run**, and whoever runs it adds two rows to §18 with the package versions attached.

  **What the gate does NOT cover, corrected 2026-08-14.** An earlier version of this bullet said §5.3's
  convergence-rule claim is "unverifiable from this side" and folded it into the R gate. Two things were
  wrong with that. The route is **not** unreachable (§5.3), so our own criteria can be exercised here
  directly; and even where experiment could not settle the attribution, `logistf` is open source and its
  convergence test can be **read** from a CRAN tarball on any machine. Definition of done 15a is that read, and it is not gated
  on R. What is left for the run is `PSweight`'s ESS convention, which is a property of what the package
  reports rather than of what its source says, and a check that `logistf`'s implemented behaviour matches
  the defaults 15a read off. The distinction was the outside voice's (§20 round 3, finding 16) and it is
  worth stating plainly: "we cannot test this" and "we have not looked" are different sentences.
- **`firthlogist` is a dependency that exists to disagree with the code, and it will rot.** It was last
  released in 2022 and it is pinned behind a scikit-learn ceiling (§16b). When it stops installing, the
  oracle to keep is §12.0.2's golden vector and §16b's scipy check, both of which are dependency-free;
  the `firthlogist` test should then be deleted rather than repaired, and this bullet is the
  instruction to do that rather than to pin harder.

## 14. What Stage 6 deliberately does not decide

- **Whether the balance achieved is adequate.** Stage 7 [§9], against |SMD| < 0.10 over the full
  confounder set. Stage 6 reports the ESS and the weights; it does not judge them.
- **What to do about a replicate whose fit fails.** [§10] and Stage 10: dropped and counted, never
  substituted. Stage 6 raises `FitError` and says nothing about what a caller does with it.
- **Whether the excluded incomplete patient should be reported in the manuscript's flow diagram.**
  Stage 14 [§16]. Stage 5's `cohort_flow` ends at 93; this stage's `covariate_completeness` says 92.
  Whether those are one diagram or two is a reporting decision.
- **Which propensity specifications [§13] eventually runs, and in what order.** Deferred by decision
  (roadmap Stage 11). §6.4 records only that they arrive as named specifications with an amendment, not
  as arguments.
- **Whether Stage 8 uses R.** Stage 6's answer — oracle, not dependency (§16b.1) — is Stage 6's, and
  it was reached on a fact specific to this stage: the Firth kernel is already validated to 2.7e-15,
  so a reference implementation has nothing left to buy here. **Stage 8 is not in that position.** The
  weighted proportional-odds model of [§8] is the primary outcome analysis and has **no** Python
  implementation — `statsmodels.OrderedModel` takes no weights — where R's `ordinal::clm(weights = )`
  does it natively. The question should therefore be asked again there, on its own evidence, and a
  reader must not treat this section as having settled it. What Stage 6 does hand forward is the
  boundary design (§16b.2) and the observation that Stage 8 arrives with a strong oracle of its own:
  [§8]'s own acceptance criterion, that the weighted fit equals an unweighted fit when every weight is
  1, is checkable against `statsmodels.OrderedModel` directly.
- **Whether our convergence rule is `logistf`'s.** Both criteria work as specified and both fire on
  real designs (§5.3), so there is no reordering question — an intermediate draft of this bullet raised
  one on the false premise that the score route was dead. What remains open is the *attribution*: whether
  `logistf` tests the same two things, in the same order, against a raw or an information-scaled score.
  Definition of done 15a settles it by reading the source; Stage 6 changes neither tolerance nor order
  either way, because both are prespecified and part of the estimator (§5.6, §15).
- **Whether the `sex` coding is 0 = female or 0 = male.** `SEX_LABELS` is `None` pending a data query
  (Stage 1 §10). The coefficient is unaffected; its *label* is not, and Stage 14 must print `sex = 1`
  rather than invent one. Named here because this is the first stage that produces a coefficient for it.

## 15. NOT in scope for Stage 6

| Considered | Why deferred |
|---|---|
| Clipping `e` into `[1e-8, 1 − 1e-8]` (`pilots/analysis.py:82`) | §5.5. F5 asserts strict interiority instead. A clip converts a degenerate fit into a finite weight every downstream sum accepts |
| A `pinv` fallback when the information matrix is singular (`pilots/analysis.py:366`) | §5.4. It returns coefficients for a design that identifies none |
| `scheme="iptw"` / `"unadjusted"` weights and `IPTW_TRIM` (`pilots/analysis.py:127-133`) | §6.1. There is no trimmed-IPTW row in this SAP, overlap weights are bounded by construction, and [§15] lists unrestricted and trimmed IPTW among the approaches deliberately not used |
| A `covariates=` parameter on `fit` (`pilots/analysis.py:953`) | §6.4. The estimand is indexed by the propensity model, so a keyword makes two target populations look like two options |
| A `propensity` or `weight` column on the cohort frame | §0.2. Stage 10 refits per replicate; a column is a second place the score lives |
| Standardised mean differences, the overlap table, per-centre diagnostics | Stage 7 [§9]. `ess` is public so that stage does not reimplement a Kish sum |
| Cross-fitting, splines, interactions, or any flexible nuisance model | [§8] and [§15]: at this treated-arm size they add variance rather than robustness, and cross-fitting folds would be too small to serve their purpose |
| Caching the fit | 92 × 12, seven iterations. A cache is a second source of truth, and here it would be a second propensity score |
| A CLI entry point | Stage 14 is the single entry point [§16]. Both modules are libraries |
| `firthlogist` (or any library) as a runtime dependency | §16b |
| **Reordering the two convergence tests, or tightening `FIRTH_TOL`, so the score route becomes reachable** | §5.3. The order and the tolerances are part of the estimator, [§10] refits in every replicate, and changing either changes the sampling distribution. It is a [§13] amendment and DoD-16's `logistf` comparison is the evidence it wants |
| **Extracting a shared table row-builder from `data.absence_by_column`** | §7.4. It is the more thorough fix for the duplication and it reopens a landed module three stages call. Filed in `TODOS.md`; §7.4's narrowed table removes the duplication without it |
| **A `Float64` (masked) dtype for `e` and `w`**, so the absence is `pd.NA` rather than `nan` | §4.4, §8. It would make the sentinel say what §4.4 means, and it would reintroduce §3.1's first hazard on every downstream `.to_numpy(dtype=float)` — on the one column most likely to be summed. `in_model` carries the intent instead |

## 16. What already exists, and what to lift

`pilots/analysis.py` is the only prior art and it is substantial: `_firth_coefs` (346-393),
`firth_logit` (48-70), `fit_ps` (73-82), `make_weights` (123-135), `ess` (138-140) and `design`
(25-41). The roadmap's standing warning applies — "This should not be a ground source […] Some
implementations may be wrong" — and here it splits cleanly: **the numerics are almost entirely right
and are lifted; every parameterisation and every fallback is declined.**

| From `pilots` | Status |
|---|---|
| `_firth_coefs`'s penalised likelihood, modified score, hat diagonal, trust region and step-halving | **lift, nearly verbatim.** Verified correct to 2.7e-15 against `firthlogist` on the real cohort design (§16b, §18). This is the one place in the repository where the pilots' code is the specification |
| its convergence criterion, and the comment explaining why the step norm is excluded | **lift both**, comment included. §5.3. The comment is the only record of a decision that would otherwise look like an oversight |
| its two different working-weight floors (`1e-12` and `1e-10`) | **lift the mechanism, not the asymmetry.** One constant, `FIRTH_WEIGHT_FLOOR`. Verified inert (§18); it is unified because a reader cannot tell that from the code |
| `PSFitError` and `fit_ps`'s no-fallback docstring, including the measured 16.6% | **lift the type and the argument.** §5.5. The 16.6% is this repository's own evidence for [§7] and belongs in the spec, the docstring and the roadmap |
| `np.linalg.pinv` on `LinAlgError` (357-358) | **do not lift.** §5.4, §15 |
| `np.clip(e, 1e-8, 1 - 1e-8)` (82) | **do not lift.** §5.5, §15 |
| `fit_ps_mle`, `fit_ps_ridge`, `_sklearn_ps` (85-120) | **do not lift** — and note *why* it is not merely that scikit-learn is absent from this project. `fit_ps_mle` falls back to `_sklearn_ps` on any exception, which is the mixture [§7] forbids; the absence of the dependency is a convenience, not the argument |
| `ess` (138-140), and its `0.0` return for an all-zero arm | **lift the formula, decline the zero.** §6.2 |
| `make_weights`' `"ato"` branch (125-126) | **lift.** It is `np.where(d == 1, 1 - e, e)` and there is nothing else it could be |
| its `scheme=` parameter and `"iptw"` / `"unadjusted"` branches | **do not lift.** §6.1, §15 |
| `design`'s constant-column drop and its docstring's reason | **lift both.** The reason — a subset or a replicate can empty a level and make the design singular — is `config.py:301-305`'s reason, written independently, and the two agreeing is worth noting |
| `design`'s factor detection by **dtype** (`CategoricalDtype or object`) | **do not lift.** `CATEGORICAL` is declared. Detection by dtype means `center` becoming `object` through an unrelated edit changes the model silently, and `onset_type` being `string` (Stage 3 §4.4) is not caught by either branch |
| `design`'s `remove_unused_categories()` + `drop_first=True` | **do not lift, and know why:** the first is needed only because the pilots made `center` a `Categorical` (`pilot_data.py:99`) and the second silently rebaselines. §4.2 |
| `prepare`'s complete-casing with `reset_index(drop=True)` (951-963) | **lift the complete-case rule, decline the mechanism.** The reset destroys the traceability Stage 5 §0.2 preserves, and applying the mask inside the function makes the exclusion invisible (§4.4) |
| `sensitivity`'s `dropna(subset=[c for c in covars if c != "center"])` (915) | **do not lift.** A complete-case rule that excludes one covariate from its own subset is a different population with no statement of why |

`pilots/test_estimators.py` is lifted more directly than any test in this repository so far: its
separation benchmark, its large-n MLE comparison and its two hand-computed ESS values become §12.4,
§12.5 and §12.9. Its `test_firth_ato_balance_is_near_exact_but_not_exact` is the origin of §6.3, and
§12.10 strengthens it by pairing the inequality with the MLE equality.

## 16b. Validation against reference implementations

**The question, put directly before a line was written: rely on well-maintained modules for the
statistical models rather than implementing them?** It is the right question — a standardised
implementation cannot carry an implementation mistake — and it was investigated rather than answered
from principle. The answer differs by component, and this section is the record so that the question
is not re-litigated at every stage that fits something.

**Where a maintained implementation exists, it is used.** `pandas.get_dummies` over a declared
`pd.Categorical` *is* the standard design-matrix construction — pandas is the maintained library, and
the pilots' defect was not hand-rolling but failing to declare the levels. `numpy.linalg` carries
every matrix operation. `statsmodels.Logit` is the unpenalised cross-check. `scipy.optimize` is the
independent optimiser check below. Nothing in this stage reimplements a linear-algebra primitive, an
optimiser it could import, or a dummy encoder.

**For the Firth fit itself, there is no maintained option.** Measured 2026-08-13:

| Candidate | Finding |
|---|---|
| `statsmodels` 0.14.6 | **no Firth at all.** `grep -ril firth` over the installed package returns nothing |
| `firthlogist` 0.5.0 | **last released 2022-08-07.** Requires **scikit-learn**, which Stage 1 excluded from this project by decision; raises `AttributeError: 'super' object has no attribute '__sklearn_tags__'` against **scikit-learn ≥ 1.6**; and `uv run --with firthlogist` **silently downgraded `tabulate` from 0.9 to 0.8.10**, violating this project's own pin. Numerically correct where it runs, and unusable as a runtime dependency |
| `pyfirth` 0.0.6 (2024-02-17) | a **two-line stub**: `name = "pyfirth"` and no code |
| R `logistf` — the reference implementation | **R 4.5.0 is installed and working on this machine**, and `logistf` is not (nor are `brglm2`, `PSweight`, `WeightIt`, `cobalt`, `MatchIt` or `survey`). R is declined as a *runtime* dependency for the reasons in §16b.1 — not because it is unavailable — and is **recommended as the authoritative oracle** there instead |
| `statsmodels.OrderedModel` (Stages 8 and 12) | `__init__(endog, exog, offset, distr)` — **no weights parameter**, which is what `pyproject.toml`'s own comment already records, and why [§8]'s weighted proportional-odds fit is implemented directly too |

**The causal-inference frameworks were checked as a second route, and they fail on the estimand
rather than on the estimator.** The question is a fair one — a framework might supply the weighting,
the ESS and the augmented estimator while taking a propensity learner of our own — so it was put to
the ecosystem rather than dismissed. Measured 2026-08-13:

| Library | Last release | Finding |
|---|---|---|
| `causallib` 0.10.0 (IBM) | 2025-04-06 | The one genuinely promising API: `IPW(learner, clip_min, clip_max, use_stabilized)` takes any learner with `predict_proba`. But it implements **IPTW only** — stabilised, clipped, truncated. **No overlap weights**; the `overlap` hits in the tree are `contrib/bicause_tree/overlap_utils.py`, a diagnostic. No Firth. And it declares `torch` plus **two conflicting `faiss-cpu` pins** (`~=1.7.0` and `~=1.8.0`) |
| `dowhy` 0.14 | 2025-11-08 | Nine propensity weighting schemes, every one an IPS variant — `ips_weight`, `ips_stabilized_weight`, `ips_normalized_weight` and their `c`/`t` forms. **No overlap weights, no Firth.** Pulls `econml`, `cvxpy`, `cython`, `causal-learn`: 19 dependencies |
| `zepid` 0.9.1 | 2022-10-23 | **Does not import on Python 3.12** — `ModuleNotFoundError: pkg_resources`. No Firth, no overlap weights |
| `econml` 0.17.0, `causalml` 0.17.0, `DoubleML` 0.11.4 | 2026 | Actively maintained, and all three are heterogeneous-effect / double-ML frameworks requiring **scikit-learn ≥ 1.6**. [§15] rules out cross-fitted machine-learning nuisance models **by name**, and [§8] gives the reason: at this treated-arm size flexible nuisance models add variance rather than robustness |
| `CausalInference` 0.1.3 | 2019-12-17 | What `wake_up_bridging/` used. Unpenalised logit with internal stepwise covariate selection, and unmaintained for seven years |
| `PSweight` / `WeightIt` — the reference implementations of overlap weighting | — | **R only.** `psweight`, `PSweight`, `pyweightit` and `overlap-weights` all 404 on PyPI |

**The pattern is one fact, not six.** Python's causal-inference ecosystem grew around inverse-probability
weighting and around ML estimation of heterogeneous effects. Overlap weighting — Li, Morgan and
Zaslavsky's `h(X) = e(X){1 − e(X)}`, which is this SAP's estimand — is served in R and has no Python
implementation at all. So the library that would have to supply it does not exist, and the ones that
exist supply the estimator [§15] lists among the approaches deliberately not used.

**And even a library that did offer it would buy the wrong three lines.** The ATO weight is
`np.where(a == 1, 1 - e, e)` and the Kish sum is one expression; what is difficult here is estimating
`e` under near-separation with one prespecified estimator and no fallback, and a framework taking a
`learner` argument requires that estimator to be written before it can be passed in. The framework
buys the trivial part and cannot supply the hard part — while adding `torch`, or scikit-learn, or
nineteen dependencies, to the reproducibility path of a 2000-replicate bootstrap.

Corroborating, from this repository: `pilots/pyproject.toml` already carries `causallib>=0.9` and
`causalinference==0.1.3` in an **optional `crosscheck` extra** — so the pilot reached the same
conclusion by a different route, and used them to disagree with its own numbers rather than to produce
them. That is the posture the four oracles below formalise.

**One sharpening of the headline agreement, from §20 round 3.** At the returned `beta` on the workbook
design, `max|modified score|` is **6.1e-3** — not near zero. Convergence here is an absolute likelihood
tolerance on a genuinely flat surface (§5.3), so "the pilots' kernel and `firthlogist` agree to 2.7e-15"
means **two implementations of the same algorithm stop in the same place**, which is a weaker statement
than "both find the optimum". The `scipy` oracle's number — coefficients agreeing to 2.7e-05 against an
independent optimiser on the same objective — is the honest one, and §16b's fourth oracle is listed with
that caveat attached rather than as a footnote. It does not weaken the case for an in-repo estimator; it
sharpens what the oracles are evidence *of*.

**And the algorithm is not where the risk lives.** The pilots' kernel and `firthlogist` were run on
the real 92 × 12 cohort design: **7 iterations each**, coefficients agreeing to **2.7e-15**, fitted
probabilities to **4.4e-16**, and ESS identical at 30.3394 / 29.1990. Both implement `logistf`'s
algorithm, which is what [§7] and the roadmap prescribe. The exposure is therefore not "did we invent
a wrong estimator" — that question is answered, twice, by measurement — but **undetected drift in the
one hand-written kernel**, and the remedy for drift is an oracle, not a dependency.

**One qualification the review added, and it narrows this paragraph rather than contradicting it.** The
2.7e-15 agreement is agreement on the *optimum*, on one design, reached by both implementations along
the easy route. It says nothing about the convergence *rule*: §5.3's second criterion is now measured to
fire only on ill-conditioned designs, and the workbook is not one — so `firthlogist` agreeing with us
here does not establish that either of us tests convergence the way `logistf` does. The arithmetic is validated; the stopping rule is
not, and no Python oracle can validate it. That is the gap DoD-16 exists to close.

### 16b.1 R, which does have all of this — and is the oracle rather than the dependency

**R is where the methodologically authoritative implementations live, and this must be said plainly
rather than worked around.** The Python situation above is not a statement about the method; it is a
statement about one language's ecosystem.

| R package | What it is |
|---|---|
| `logistf` (Heinze & Ploner) | **The reference Firth implementation.** §5's *penalty and modified score* are `logistf`'s; `firthlogist` is a port of it, which is why the 2.7e-15 agreement below is transitive rather than independent. §5's **convergence rule** is a separate claim, carried from the pilots' comment and unverified — see §16b.2's three questions and DoD-16 |
| `brglm2` (Kosmidis) | Bias reduction in the wider sense — mean- and median-bias-reduced GLMs, actively developed, and arguably the better modern reference for the same penalty |
| `PSweight` (Zhou, Matsouaka & Li) | **The reference implementation of this SAP's exact estimand.** Overlap weights, balance, augmented estimators and bootstrap inference, written by the authors of the overlap-weight method. Nothing in Python touches it |
| `WeightIt` + `cobalt` (Greifer) | `estimand = "ATO"` with a large family of propensity models, and the balance tooling Stage 7 [§9] would otherwise write by hand |

**Verified 2026-08-13: R 4.5.0 is installed on this machine and none of those four packages is.**
(An earlier draft of this section recorded R as segfaulting on startup. That was wrong: R itself runs;
what fails is `system()` from inside the sandboxed shell the probe ran in, which prevents `utils` from
loading and so prevented the comparison below from being executed. The correction matters, because the
conclusion it supported — "R is unusable here" — is false, and the real conclusion is a trade rather
than an absence.)

**Decided 2026-08-13: R is an oracle and not a dependency — its role is to establish that this
implementation is correct, not to produce any number that is reported.** The question was put
directly, with the case for adoption stated at its strongest (`PSweight` is the reference
implementation of this SAP's own estimand, written by the method's authors), and the decision is
recorded here rather than left as a recommendation, so that a later reader finds a choice with a date
on it rather than an argument.

**Five reasons, in descending weight:**

1. **It doubles the reproducibility surface.** This project hashes the workbook, pins pandas and numpy
   with reasoned upper bounds, and asserts a byte-identical audit log. A second language runtime with
   an independent package manager, an unpinned CRAN, and its own version drift is the largest
   reproducibility liability available — and [§16]'s methods section would have to state provenance for
   two toolchains.
2. **The bridge is the fragile part, not the statistics.** `rpy2` requires R built with
   `--enable-R-shlib` and is historically brittle across R × Python × numpy combinations; the
   alternative is a subprocess and a serialisation boundary crossed on **every one of `N_BOOT`
   replicates**, refitting both nuisance models each time [§10].
3. **It would not be one call.** [§8]'s weighted proportional-odds fit has no Python implementation
   either, so the boundary would be crossed for Stage 8 and Stage 12 as well. That is not "use a
   library for the propensity model"; it is a second pipeline in a second language.
4. **Nothing is installed**, so adoption starts with a fresh install of four packages and their
   compilation toolchains — and then needs `renv` or a pinned CRAN snapshot to be reproducible at all.
5. **`PSweight` would take the estimation, not just the estimator.** Its value is that it owns the
   whole ATO pipeline — which means adopting it is adopting its conventions for weights, variance and
   balance, and [§7], [§8] and [§10] are already specified down to the p-value's definition.

**As an oracle, R is better than anything in §16b's list and is recommended as a fifth.** The four
oracles below are all transitive or partial: `firthlogist` is a port of `logistf`, `scipy` checks our
own objective, the golden vector checks us against ourselves, and the MLE check only bites at large n.
A direct `logistf` comparison is the one independent check on the estimator, and a `PSweight`
comparison is the only available check on the **whole [§7] pipeline** — score, weights and ESS
together — by the people who defined the estimand.

### 16b.2 How R fits a mainly-Python pipeline: a gated subprocess oracle

**The shape that keeps both properties** — R's authority and a Python-only runtime — is the one
`firthlogist` already has, one level further out: a test that shells out to `Rscript`, skips cleanly
when R or its packages are absent, and is never on any estimation path.

```
  extended_bridging/
    tests/reference/firth_logistf.R      committed, ~20 lines, reviewable in a diff
    tests/reference/ato_psweight.R       committed
    test_reference_r.py                  the only Python that knows R exists
```

```python
# test_reference_r.py — the whole gate
RSCRIPT = shutil.which("Rscript")
R_PACKAGES = ("logistf", "PSweight")

reference_r = pytest.mark.skipif(
    RSCRIPT is None or not _r_has(R_PACKAGES),
    reason=f"R oracle: needs Rscript and {', '.join(R_PACKAGES)}; see spec §16b.2")
```

Six properties, each of which is a decision rather than an implementation detail:

- **The boundary is a CSV written to `tmp_path` and a CSV read back**, never `rpy2`. No embedded
  interpreter, no shared-library build requirement, no R × Python × numpy compatibility matrix, and
  nothing to install into the project environment. The cost is one process spawn per test, in a test
  that runs once — not 2000 times inside a bootstrap, which is exactly the distinction §16b.1 turns on.
- **Floats cross at full precision.** `to_csv(float_format="%.17g")` and `options(digits = 17)` on the
  way back, so a disagreement is the estimator's and never the serialisation's. A comparison that
  round-trips at 6 significant figures cannot distinguish 2.7e-15 agreement from 1e-7 agreement, and
  the first is the claim being made.
- **`Rscript --vanilla`**, so no user or site profile can change the result, and the R script sets its
  own seed-free deterministic path. This is the `PYTHONHASHSEED` discipline of Stage 2 §11 applied
  across the boundary.
- **The R scripts are committed and short.** An oracle nobody can read is not evidence. Twenty lines
  each: read CSV, fit, write coefficients. They carry no analysis logic and no covariate list — the
  design matrix arrives already built by `model.design`, so the R side cannot disagree about the
  *specification* while appearing to disagree about the *estimator*.
- **The R package versions are captured and asserted non-empty**, then written into §18 by the person
  who runs it. An oracle that agreed at an unrecorded version is a claim with no date on it.
- **`test_reference_r.py` is imported by nothing** and appears in no shipped module. §12.7's scan is
  extended to assert that `subprocess`, `Rscript` and the two R package names appear in this file and
  nowhere else under `extended_bridging/`.

**And the skip must be loud, because a silently skipped oracle is the failure this repository already
guards against elsewhere.** Every AST scan in Stages 1-5 has a companion test proving it fires,
precisely so that a check matching nothing cannot pass as green. The same hazard applies here in a
worse form: the R oracle skips on a machine where R is absent, which is *most* machines, so the
default state is "not run". Three mitigations, all required:

1. The skip reason names what is missing and points at this section, so `pytest -rs` reads as an
   instruction rather than as noise.
2. §18 carries a row recording the versions it last agreed at — so "when did this last actually run"
   is answerable from the document.
3. **Definition of done items 15a and 16 split the evidence**: 15a reads `logistf`'s convergence rule from source and needs no R; 16 requires the comparison to have been seen to run and to pass at least once, on a
   machine with R, before Stage 6 is called complete. A gate that has never opened is not a gate.

**Two mismatches to expect, neither of which is a bug in either implementation**, and the oracle is
worth more for finding them than for confirming agreement:

- **`PSweight`'s effective sample size may not be Kish's.** [§7] asks for an ESS and §6.2 computes
  `(Σw)²/Σw²`; a package owning the whole pipeline may report a different one, or report it after its
  own normalisation of the weights. Compare `e` and `w` **first** and the ESS separately, so a
  convention difference cannot be read as a weighting error.
- **`logistf`'s control defaults are not §5.6's constants.** They must be passed explicitly through
  `logistf.control(...)` to make the comparison like-for-like, and §5.3's claim that our two
  convergence criteria are "what `logistf` does" is **carried from `pilots/analysis.py`'s comment and
  has not been checked against the package** — `logistf` may also test the coefficient step, which
  §5.3 deliberately excludes. **This is the single most valuable thing the R oracle would settle**, and
  it is a claim about our own convergence rule rather than about our arithmetic. Recorded in §13.

  Both our criteria fire (§5.3), so this side can exercise them — what it cannot establish is whether
  `logistf` uses the same pair, in the same order, against a raw or an information-scaled score, or
  whether it tests the coefficient step as a third.
  Three things to read off `logistf.control` and its source while the gate is open, all of them cheap
  once R is installed: **(a)** the default `lconv`, `gconv` and `xconv` and whether `xconv` is a
  coefficient-step test that §5.3 excludes; **(b)** the order the checks run in; **(c)** whether `gconv`
  is compared against the raw modified score or against a scaled version of it — a score normalised by
  the information matrix would be reachable where ours is not, and that single difference would explain
  the whole finding. If any of the three differs, §5.3's argument survives untouched — excluding the step
  norm is still right, for the near-separation reason — but the sentence attributing the rule to
  `logistf` is struck rather than softened, and §5.6's constants get a comment saying whose they are.

**Filed as a blocker, by decision of 2026-08-14.** It does not gate T2-T5 — those land and go green
without R — but it gates *completion*, which is what Definition of done 16 says and what §13 now says
too. An earlier draft of this paragraph called it "an open item rather than a blocker" while DoD-16
treated it as a gate; the document argued both sides and the gate wins. The reason it wins is §5.3: the
score-route attribution is unverified: our two criteria both fire (§5.3), but whether `logistf` uses the
same pair in the same order is a fact about its source, not about our runs. The four oracles below can
establish that our arithmetic is right; only `logistf` can establish that our *convergence rule* is the
one [§7] means, and Definition of done 15a is the cheap half of that.

**Four oracles, none of which may touch a reported estimate.**

1. **`firthlogist`, dev-only.** A `[dependency-groups] reference` entry pinned as
   `firthlogist>=0.5,<0.6` with `scikit-learn<1.6`, imported only inside `test_model.py` behind
   `pytest.importorskip("firthlogist")`. Asserted to agree with `model.firth` to **1e-12** on the
   cohort design (data-gated) and on `separated()` (always). The primary path never imports it, and
   §12.7's scan asserts that.
2. **A golden vector with no dependency at all.** §12.0.2 pins five fitted probabilities from the
   synthetic frame as literals, so an edit to the kernel is caught on a plain checkout with no
   `data/` and no reference package. This is the oracle that survives `firthlogist` rotting (§13).
3. **`scipy.optimize.minimize` on the *same* penalised objective.** Asserts that the Newton loop reaches
   the objective's optimum **on the designs this pipeline fits** — which validates the loop independently
   of any Firth implementation, and is the check no external package can provide, because an external
   package would be validating its own algorithm rather than our objective. **Scoped deliberately:** an
   earlier wording said the loop's optimum *is* the objective's optimum, without qualification, and that
   is false on separated designs where the penalised likelihood is multimodal (§5.2c). The oracle is
   asserted on the cohort design and on §12.0's fixtures, where it holds; §5.2c documents where it does
   not, and §9 assigns the general case to Stage 10. Measured (§18): BFGS reaches a
   penalised log-likelihood of **−19.0637927357** against the Newton loop's **−19.0637927357** on the
   cohort design, with coefficients agreeing to 2.7e-05 — the looser coefficient agreement being the
   flat surface §5.3 describes, which is why the assertion is on the **objective value** and not on
   the coefficients.
4. **`statsmodels.Logit` at large n**, §12.4 — the roadmap's own criterion.

**The standing rule this section establishes**, for Stages 8, 9 and 12 to inherit: use a maintained
library wherever one exists; where none does, the estimator carries an oracle rather than a comment,
and the oracle is named in the spec.

## 17. Implementation tasks

Ordered. Each independently verifiable. T1 is the amendment that unblocks everything, T2-T5 are Stage
6's own, T6 is the R oracle, T7 is the log and the definition of done.

**T1 must land first** because `Audit.record("model", …)` raises `ValueError` until `KINDS` carries the
kind, so nothing below it can be run even once. T2 and T3 are genuinely independent of each other —
`design` and `firth` share no state — and are the one place in this stage where two worktrees would
not cost more than they save. T4 depends on both. **T6 depends on nothing but T3** — it reads
`model.firth`'s output through a CSV and changes no shipped code — so it can be done in parallel with
T4 and T5. It cannot be *deferred*: it gates Definition of done 16, and §16b.2 and §13 both now call it
a blocker rather than an open item.

- [x] **T1 (P1)** — `data.py`: `KINDS` and `_HEADINGS` gain `model` in §7.1's position; the
      `_MUST_NAME_CASES` comment gains its sentence and its count goes six → seven. `config.py`: the
      seven-constant block of §5.6. `test_data.py`, `test_eligibility.py`, `test_cohort.py`: the three
      literal pins and the two renames (§10). `test_config.py`: the one
      `REFERENCE_LEVELS[f] == FACTOR_LEVELS[f][0]` assertion, added to the existing per-factor
      parametrised test (§4.2, §10). `pyproject.toml`: the `reference` group and the two
      direct-dependency comments. Verify: `uv run pytest -v` green with **no other test edited**,
      which is where §10's claim that the kind-parametrised tests need no change is checked rather
      than believed — and where the new `test_config.py` assertion is confirmed to pass on the config
      as it stands, since it would be a landed-code defect if it did not.
- [x] **T2 (P1)** — `model.design`, `model.complete_cases`, `_assert_design_inputs` with D1-D4.
      Verify: acceptance 12.1, 12.2, 12.6's D-half (all four), and 12.3's design-keeps-the-nan-row half. Watch
      §12.1's column-order literal specifically: `get_dummies` appends the dummies rather than
      expanding each factor in place, and the literal in §12.1 is the measured order, not the
      covariate order.
- [x] **T3 (P1)** — `model.firth`, `Fit` **with its eight fields** (§3.2's three counters included),
      `FitError`, `_penalised_loglik`, `_probabilities`, and `_assert_fittable` with **F3, F4 and F6**.
      Verify: acceptance 12.4, **12.4a**, 12.5, 12.6's F-half, 12.7. Three things to get exactly right,
      each of which was a defect in an earlier draft of this document: the trust radius is **relative**,
      `‖step‖ / max(1, ‖β‖)` and not `‖step‖` (§5.2a); `moved` is bound **before** `ll_old = ll_new`
      (§5.2b); and `_probabilities` applies the same `FIRTH_ETA_CLIP` the loop does, or F5 tests a quantity
      the weights are not built from (§5.2).
- [x] **T4 (P1)** — `propensity.fit` as written in §8, `Propensity`, `ess`, `_assert_fit_inputs`,
      `_assert_probabilities`, `_record_exclusion`, and §7.2's four entries with their table helpers —
      including the **reference-level rows** of §7.5 and the **weighted-population table** of §7.6, both
      of which are new since the review and neither of which any earlier draft specified.
      Verify: acceptance 12.8, 12.9, 12.10, 12.11, 12.13.
- [x] **T5 (P2)** — `test_model.py` and `test_propensity.py` remainder: the synthetic fixtures of
      §12.0 **written out as specified there**, the hard-fit family of §12.0's 1b and its §12.4a
      assertions, the golden vector of §12.0.2, `test_propensity.py`'s own module-scoped `workbook`
      fixture, every AST scan and its companion, the two-seed driver, §16b's four oracles, and the
      data-gated facts of §12.14. Verify: `uv run pytest -v` green both with and without `data/`, and
      green with the `reference` group **not** installed.
- [x] **T6 (P1, was P2)** — the R oracle of §16b.2: the two committed `.R` scripts,
      `test_reference_r.py` with its `Rscript`-and-packages gate, and §12.7's scan extended to pin R's
      blast radius to that one file. Verify: green **with** R and the packages present (the comparison
      runs), green **without** them (it skips, and `pytest -rs` names what is missing). **Independent
      of T2-T5** in its code — it reads `model.firth`'s output through a CSV and changes nothing — so it
      parallelises freely; but it is **P1 because DoD-16 gates completion on it**, and it now also
      carries §16b.2's three `logistf.control` questions, which are what settle §5.3. Installing
      `logistf` and `PSweight` is part of this task, not a precondition someone else supplies.
- [x] **T7 (P1)** — the audit log produced through all six stages against the workbook, read end to
      end by a human, and the definition of done below.

### What can be built in parallel, and what cannot

| Task | Modules touched | Depends on |
|---|---|---|
| T1 | `data.py`, `config.py`, `pyproject.toml`, four existing test modules | — |
| T2 | `model.py` (`design`, `complete_cases`), `test_model.py` | T1 |
| T3 | `model.py` (`firth` and friends), `test_model.py` | T1 |
| T4 | `propensity.py`, `test_propensity.py` | T2, T3 |
| T5 | `test_model.py`, `test_propensity.py` | T2, T3, T4 |
| T6 | `test_reference_r.py`, `tests/reference/*.R` | T3 |
| T7 | none — it runs the pipeline and reads the log | T1-T5 |

```
  Lane A:  T1  →  T2  ┐
  Lane B:          T3 ┼→  T4  →  T5  →  T7
  Lane C:          T6 ┘   (T6 joins only at the suite, not at T4)
```

**T2 and T3 both touch `model.py`, and they are still the one place two worktrees pay for themselves.**
They share a file but no state: `design` and `firth` have no call between them, no shared constant beyond
`config`, and disjoint test sections. The merge is two function bodies appended to one module. Everything
else is sequential — T4 needs both, and T5 needs T4's audit entries to exist before it can assert their
positions.

**T6 is the only lane that never merges into another.** It adds one test module and two `.R` files, imports
`model` and nothing else, and changes no shipped code — so it can run start-to-finish beside T4 and T5. It
is still P1, because DoD-16 gates completion on it (§16b.2, §13); "parallelisable" and "optional" are
different properties and this stage is the place that distinction is easiest to get wrong.

**No lane conflicts.** The only shared file between concurrent lanes is `model.py` between T2 and T3, which
is flagged above rather than avoided. If the two are run in separate worktrees, merge T2 first: T3's
`_assert_fittable` references `X.columns` from a design T2 defines the shape of, so a conflict there is a
signal that §12.1's column-order rule was implemented twice.

### Diagrams and comments that belong in the code, not only here

Three, and no more — a diagram nobody maintains is worse than none, because it is believed. Keeping
them true is part of any change that touches them, in the same commit.

- **`model.py`'s module docstring** — §0's two entry arrows: Stages 8, 9 and 12 come into this module
  directly and never through `propensity.py`. That is the sentence that stops someone moving `design`
  into `propensity.py` for tidiness.
- **`propensity.py`'s module docstring** — that `e` and `w` are Series over the *whole* cohort carrying
  `nan` off `in_model`, that the deliberateness of that absence lives in `in_model` and not in the
  sentinel, and that a weight of zero and an absent weight are different facts.
- **Above `firth`'s convergence test** — §5.3's two routes, the near-separation reason the step norm is
  not a third, **and the two measured facts**: that the score is evaluated at the pre-step iterate, and
  when the score route fires — on an ill-conditioned design, per §5.3's `Δl* ≈ ½·s²·‖I⁻¹‖`. Lift the pilots' comment for
  the first part; the rest is this stage's own and is the only thing standing between a future reader
  and the conclusion that the second criterion is decoration.

One further comment, not a diagram: **above `_assert_fittable`, the three silent failures** —

```
#   Int64 + pd.NA  --to_numpy(float)-->  nan   SILENTLY, no warning
#   one nan in X   --firth-->  all-nan beta, all-nan p, and the loop RETURNS NORMALLY
#                              with a plausible iteration count.        [Stage 6 §3.1]
#   one nan in y   --firth-->  raises, but with the step-halving message: a sentence
#                              about [§7]'s single estimator for a missing outcome.
#   y = mrs_90d    --firth-->  FITS. Converges. Returns coefficients for a logistic
#                              model of a seven-level ordinal, with nothing missing
#                              and nothing to notice.                   [Stage 6 §5.4]
#   So this function is not defensive. F4 is the only thing between a missing covariate
#   and a [nan, nan] confidence interval in Stage 10; F6 is the only thing between an
#   ordinal response and a table of confident nonsense. Stages 8, 9 and 12 call firth
#   directly [§0.1, §9] — there is no second checkpoint downstream of this one.
```

### Definition of done

Stage 6 is complete when all of the following hold, and not before.

1. `uv run pytest -v` is green **with** `data/` present and the `reference` group installed.
2. Green **with `data/` temporarily renamed** — the data-gated tests skip and nothing else fails or
   errors at collection.
3. Green **with the `reference` group not installed** — §16b's first oracle skips, the other three run.
4. `test_config.py`, `test_data.py`, `test_derive.py`, `test_eligibility.py` and `test_cohort.py` are
   still green, and neither `model.py`, `propensity.py`, `test_model.py` nor `test_propensity.py` has
   been added to `EXEMPT_FROM_RAW_NAME_SCAN`.
5. The pipeline's audit log has been produced through all six stages, against the **workbook**, and
   reproduces byte-identically across hash seeds:

   ```bash
   cd extended_bridging
   S='import sys, data, derive, eligibility, cohort, propensity
   df, audit = data.load(data.WORKBOOK)
   df = derive.derive(df, audit)
   df = eligibility.classify(df, audit)
   df = cohort.build(df, audit)
   propensity.fit(df, audit)
   sys.stdout.write(audit.to_markdown())'
   PYTHONHASHSEED=0 uv run python -c "$S" > /tmp/a.md
   PYTHONHASHSEED=1 uv run python -c "$S" > /tmp/b.md
   diff /tmp/a.md /tmp/b.md
   ```

   `diff` reports nothing, and the four `model` entries render under `## Fitted models` between
   `## Cohort construction` and `## Structural non-applicability`. The log names patients and must
   never be committed.
6. **F4 has been seen to fail, and then seen to be necessary.** With `_assert_fittable` removed and one
   `core_ml` value blanked, `fit` returns a `Propensity` whose `e` is `nan` throughout, whose
   `iterations` is a plausible small integer, and which raises nothing. Watch this one specifically: it
   is the failure mode this stage exists to make loud, and every acceptance test downstream of it would
   pass on `nan`.
7. **D3 has been seen to fail the same way.** With D3 removed and one record's centre set to a value
   outside `FACTOR_LEVELS`, the fit returns a finite propensity for that record, computed as though it
   were at HUG, with no missing value anywhere in the design.
8. **F3 has been seen to fail in both directions.** With `_assert_fittable` removed and an exactly
   duplicated column, `np.linalg.inv` raises `LinAlgError` from inside the loop; with a `pinv` fallback
   added, the fit returns finite coefficients for a design that identifies none. The second is the one
   to watch — it is what the pilots ship.
9. **F5 has been seen to matter.** With the assertion replaced by `np.clip(e, 1e-8, 1 - 1e-8)` on a
   frame contrived to overflow the linear predictor, `fit` returns weights of exactly 0.0 and an ESS
   that is finite and wrong.
10. **The reference-by-name drop has been seen to fail loudly on a frame with no HUG record — as F3, not
    as `KeyError`.** Watch all three halves, because the first is the one an earlier draft of this
    document got wrong: (a) `design` does **not** raise, because the declared `Categorical` emits an
    all-zero `center_HUG` that is dropped as the reference; (b) `firth` then raises **F3** at rank
    ***k*−1** of *k*, because the surviving centre dummies sum to 1 on every row and are collinear with
    the intercept — that is **one** dependency, so k−1. An earlier draft said k−2, generalising from a
    toy frame that happened to carry a second accidental collinearity (`onset_type_wake_up` identical to
    `center_Lugano`); (c) `drop_first=True`
    on the same **declared** frame drops the same column and is therefore *equivalent*, while on a
    **bare** `string` centre it drops the first present level and silently rebaselines. If (a) raises a
    `KeyError`, the `Categorical` conversion is missing and §4.2's first bullet has been dropped.
10a. **F6 has been seen to matter, and its ordinal half is the one to watch.** With F6 removed:
    passing a response containing one `nan` raises the step-halving message at iteration 1 — loud, wrong
    diagnosis — and passing `mrs_90d`'s 0-6 response **raises nothing at all**, converging and returning
    coefficients for a logistic model of a seven-level ordinal. The second is the only silent failure
    this stage adds a guard for that no earlier stage would have caught, and it is the one Stage 9 will
    meet first (§5.4).
10b. **§12.4a's safeguards have each been seen to fire, and the two unreachable branches seen not to.**
    `big_first_step()` rescales (2, measured), `needs_halving()` halves (1, with 0 rescales), and both
    raises fire under their monkeypatched constants — with `FIRTH_MAX_ITER`'s message carrying a
    **non-zero** likelihood movement, which is §5.2b's correction and the thing to watch: the earlier
    form printed `0` on every input. Then: over every fixture in `test_model.py` plus the workbook
    design, `min p(1−p)` is above `FIRTH_WEIGHT_FLOOR` — so the floor's clip is confirmed unreached rather
    than assumed so (§5.6). `converged_on` is asserted **per fixture against both routes** (§5.3), not
    globally as `"likelihood"`. If the floor assertion fails, a constant moved and §5.6 needs re-reading,
    not the test.
10c. **The trust region has been seen to be scale-invariant, and the absolute form seen not to be.**
    `big_first_step()` with its covariate multiplied by 1e-3, 1, 1e3 and 1e6 converges in every case
    under §5.2a's relative radius. The fitted probabilities agree to **1e-3**, not 1e-8 — the surface is
    flat and the four stop at different points on the same ridge (§12.4a). With the radius reverted
    to `if ‖step‖ > FIRTH_MAX_STEP`, the 1e-3 case fails at `FIRTH_MAX_ITER` and the 1e-1 case takes 84
    iterations. **Watch this one alongside DoD-6**: it is the second failure in this stage where a
    correct answer exists and the code reports a failure instead, and [§10] would have dropped those
    replicates and counted them.
10d. **The golden vector and the workbook fit have been confirmed unmoved by §5.2a.** §12.0.2's five
    probabilities and weights are identical to all six pinned decimals, at 11 iterations and one
    rescale; the workbook fit is 7 iterations with `max|Δβ| = 0` against the absolute radius. If either
    moves, the relative radius was implemented as something other than
    `‖step‖ / max(1, ‖β‖)` — most likely `‖step‖ / ‖β‖`, which is unbounded at `β = 0`.
11. **§16b's first oracle has been watched agreeing**, on the workbook design, at 1e-12 — and the
    `firthlogist` version and scikit-learn version it agreed at recorded in §18 rather than assumed
    from `pyproject.toml`.
12. **§12.10's pair has been seen to discriminate**: with `e` swapped to a `statsmodels.Logit` score,
    the Firth inequality assertion fails and the MLE equality assertion passes. A test that cannot fail
    when the estimator is substituted is not evidence about the estimator.
13. **§12.11's entry-position assertions use an index captured before the call**, verified by grep: no
    tail slice of the form `[-n:]` appears in either new test file. This is Stage 5 §10's repair, and
    this stage inserts four entries where it would bite.
14. The audit log has been read end to end by a human, and its numbers agree with §18 and with
    `../out/stage0_data_inventory.md`: the fit's `n` equals `cohort_flow`'s `records`, minus
    `covariate_completeness`'s `n`.
15. `grep -rn "sklearn\|firthlogist" extended_bridging/*.py` returns hits in `test_model.py` only, and
    `grep -rn "Rscript\|logistf\|PSweight" extended_bridging/*.py` in `test_reference_r.py` only.
15a. **`logistf`'s convergence rule has been READ, from source, and recorded in §18.** This does not need
    R installed: `logistf`'s CRAN tarball is a download and its convergence test is a few lines of R.
    Record four things in §18 with the package version: the default `lconv`, `gconv` and `xconv`;
    whether `xconv` is a coefficient-step test §5.3 deliberately excludes; the **order** the checks run
    in; and whether `gconv` compares the raw modified score or an information-scaled one — because a
    scaled score would be reachable where §5.3 measures ours is not, and that single difference would
    explain the whole finding. Then either §5.3's attribution is earned and restored, or it is struck.

    **This item exists because the review's outside voice was right and §13 was wrong.** §13 called the
    convergence rule "unverifiable from this side"; that is true of *experiment* and false of *reading*,
    and the distinction had been collapsed. Splitting it means §5.3 can be settled this week, on a
    laptop, while item 16 keeps gating what genuinely needs execution.
16. **The R oracle has been seen to run and to pass at least once**, on a machine with `Rscript`,
    `logistf` and `PSweight` installed — and the versions of all three are recorded in §18 on the day
    it ran. A gate that has never opened is not a gate, and this one's default state is *skipped*
    (§16b.2). **This item is a blocker and both §13 and §16b.2 now say so** — an earlier draft of those
    two sections called the R oracle "an open item rather than a blocker" while this item required it,
    and the contradiction was resolved in favour of the gate on 2026-08-14. Installing the two packages
    is T6's work, not a precondition someone else supplies. **What it gates was narrowed on 2026-08-14
    to the two things that genuinely need execution**, the convergence rule having moved to item 15a
    where reading settles it:
    - **whether `PSweight`'s reported ESS is Kish's.** A convention, not a document — no source read
      substitutes for running it, because the question is what the package *reports* after its own
      normalisation. Compare `e` and `w` first; a convention difference in the third quantity must not
      be reported as a weighting error.
    - **whether `logistf` and we agree on the coefficients** of the workbook design to the tolerance
      §16b.2's CSV boundary supports — and, while it is open, whether `logistf`'s *implemented*
      behaviour matches the defaults item 15a read off its source. Documented defaults and shipped code
      diverge often enough that the read and the run answer slightly different questions, which is why
      both items exist rather than one.

## 18. Verification record

Checked on 2026-08-13, before this spec was finalised, by a read-only probe reporting aggregate
quantities only — no outcome by arm, in the posture Stages 0, 4 and 5 used. Recorded so a reader can
tell which numbers were verified rather than carried, and so the checks are re-runnable after a
workbook update. Acceptance 12.14 re-checks them in code.

The probe drives `load → derive → classify → build` against `data.WORKBOOK` and applies §4 and §5; it
writes nothing.

| Claim | Where used | Verified |
|---|---|---|
| The cohort is complete on 9 [§6] covariates for **92 of 93** records; the one exception is a **control at Lugano** missing `core_ml` and `tmax6_ml` (and `penumbra_ml`, `core_above_median`) | §4.4, §11, §12.14, §13 | yes, run |
| All four of `PS_COVARIATES_FULL`'s additions are complete on every cohort record, so the [§13] full-covariate model excludes the same one patient | §13 | yes, run |
| The design over `PS_COVARIATES` is **12 columns before the constant drop**, `center_USZ` is all-zero and dropped, leaving **11 + intercept = 12 parameters**; rank **12 of 12**, condition number **997** | §4.2, §4.3, §5.4, §11, §12.14 | yes, run |
| 39 treated over 12 parameters is **3.25 treated per parameter** | §11, §12.14, §13 | yes, computed |
| Levels present in the cohort: `center` HUG 41 / Lugano 31 / CHUV 21, USZ absent; `onset_type` wake_up 40 / unwitnessed 31 / witnessed 22, none missing | §4.2, §4.5 D3, §13 | yes, run |
| `constant_covariates` over the cohort is **empty** for both `PS_COVARIATES` and `PS_COVARIATES_FULL`, while `design`'s `dropped` is not — the detect-versus-drop distinction | §4.3, §10, §12.14 | yes, run |
| The Firth fit converges in **7 iterations**, on the **likelihood** route; `e ∈ [0.0269724, 0.9310700]`, all finite and strictly interior | §5.3, §11, §12.14 | yes, run; **max corrected 2026-08-14** from 0.93110, which was wrong in the fifth decimal |
| Overlap weights: **max 0.965**, none zero; `Σw` = **13.626** treated against **14.110** control — not equal | §6.1, §6.3, §12.10, §12.14 | yes, run |
| Kish ESS **30.3394 of 39** treated, **29.1990 of 53** control | §6.2, §11, §12.14 | yes, run |
| The unpenalised MLE **converges on the point sample** (max\|coef\| 4.55) and differs from Firth by up to **0.79**, so separation is a bootstrap-replicate risk rather than a point-fit one | §5.1, §13 | yes, run |
| `firthlogist` 0.5.0 (under scikit-learn 1.5.2, numpy 1.26.4) agrees with the lifted kernel on the cohort design: **7 iterations each**, coefficients to **2.665e-15**, probabilities to **4.441e-16**, ESS identical at 30.3394 / 29.1990 | §16, §16b, §12.14 | yes, run |
| `firthlogist` under **scikit-learn ≥ 1.6** raises `AttributeError: 'super' object has no attribute '__sklearn_tags__'`; `uv run --with firthlogist` resolved **tabulate 0.8.10**, below this project's `tabulate>=0.9` | §2, §16b | yes, run |
| `statsmodels` 0.14.6 contains **no Firth**: `grep -ril firth` over the installed package returns nothing | §16b | yes, searched |
| `pyfirth` 0.0.6 is a two-line stub exporting only `name` | §16b | yes, installed and read |
| **R 4.5.0 is installed and working**; `logistf`, `brglm2`, `PSweight`, `WeightIt`, `cobalt`, `MatchIt` and `survey` are **all absent** from its three library paths. The comparison could not be run from the probe's sandboxed shell, where `system()` fails and `utils` therefore does not load — a property of the probe, not of the machine. *(This row replaces one reading "R and Rscript segfault on startup", which attributed the shell's failure to R and supported the false conclusion that R is unusable here. §16b.1.)* | §16b.1, §13 | yes, re-run 2026-08-13 |
| `statsmodels.miscmodels.ordinal_model.OrderedModel.__init__` takes `(endog, exog, offset, distr)` — **no weights**; `statsmodels.GLM` does take `freq_weights` / `var_weights` | §16b | yes, read |
| `scipy.optimize.minimize(method="BFGS")` on the **same** penalised objective reaches **−19.0637927357** against the Newton loop's **−19.0637927357**, coefficients agreeing to **2.710e-05** | §16b | yes, run |
| On `well_behaved()` **as §12.0 now specifies it** — n = 4000, three standard-normal covariates, β = (0.5, −0.4, 0.2), intercept 0.3, `default_rng(C.SEED)` — `max abs` difference between the Firth and `statsmodels.Logit` coefficients is **0.000787**, and Firth converges in **4 iterations** | §12.4 | yes, re-run 2026-08-14. *(This row replaces one reading 0.00116 in 5 iterations, measured on a generator the spec never wrote down — see §20.)* |
| On `x = 1…10`, `y = [0]*5 + [1]*5`: Firth gives intercept **−5.3385**, slope **0.9706** in 9 iterations, `e` interior; `statsmodels.Logit` Newton reports **`converged=False`** with slope **71.5**, and BFGS returns **27.2** | §12.5 | yes, run, both methods |
| `Series([1, pd.NA], dtype="Int64").to_numpy(dtype=float)` yields `[1., nan]` with **no raise and no warning** | §3.1, §4.4, §12.3 | yes, run |
| `firth` with one `nan` in `X` returns **all-nan** coefficients and probabilities and **terminates normally**, emitting two `RuntimeWarning`s | §3.1, §12.3, DoD-6 | yes, run |
| `pd.Categorical(s, categories=FACTOR_LEVELS["center"])` maps an unlisted value to `NaN`, and `get_dummies` then gives it an **all-zero row** — identical to the reference level's encoding | §3.1, §4.5 D3, §12.6 | yes, run |
| `get_dummies` on the **bare `string`** columns emits `center_CHUV, center_HUG, center_Lugano` and **no `center_USZ`**, so without the `Categorical` conversion no constant column exists to drop | §4.2, §12.2 | yes, run |
| `np.linalg.slogdet(singular)` returns `(0.0, -inf)` without raising; `np.linalg.inv` on an exactly singular information matrix **raises `LinAlgError`**, at condition number 8.839e17 | §3.1, §5.2, §5.4, §12.6 | yes, run |
| `1/(1+exp(-500))` is **exactly 1.0**; `1/(1+exp(500))` is 7.125e-218, so the clip is reachable at the upper end only | §3.1, §5.5 | yes, run |
| `(Σw)²/Σw²` with all-zero weights is `nan` from 0/0 | §3.1, §6.2 | yes, run |
| `ess([1]*10) = 10.0` and `ess([1,1,1,3]) = 3.0`, both exactly | §6.2, §12.9 | yes, run |
| pandas **2.3.3**, numpy **2.0.2**, statsmodels **0.14.6** in the project environment | §2, §3.1 | yes, read |

**Claims about the test frames**, run against the committed repository, so they are re-checkable with
no `data/`:

| Claim | Where used | Verified |
|---|---|---|
| `test_cohort.cohort_frame()`'s cohort is **5 records, complete on all 9 covariates**; the design's twelve candidate columns reduce to **`core_ml`, `tmax6_ml`, `center_CHUV`** — nine dropped as constant — giving 4 parameters over 5 records at **full rank** | §12.0, §12.1 | yes, run |
| That fit converges in **11 iterations**, `e ∈ (0.2571, 0.9290)` strictly interior, ESS **2.4413 of 3** treated and **1.9934 of 2** control | §12.0.2 | yes, run |
| The five fitted probabilities and weights of §12.0.2, by identifier, with the cohort's index preserved as `(0, 1, 4, 6, 8)` | §12.0.2 | yes, run |
| `HAND-2`'s weight and `COHORT-3`'s weight and propensity are all **0.288470**, because those two CHUV records have linear predictors of equal magnitude and opposite sign — so an assertion by value could confuse them | §12.0.2 | yes, run |
| `hand_frame()` and `tests/fixture_schema.xlsx` still cannot reach this stage: Stage 5's P2 raises before `fit` is called | §12.0 | yes, carried from Stage 5 §18, re-read |

**Claims checked by reading the committed repository:**

| Claim | Where used | Verified |
|---|---|---|
| `data.KINDS` is the eight of Stage 5 §6.2 in that order; `_HEADINGS` keys equal `KINDS`; `_MUST_NAME_CASES` is `{correction, observation}` and its comment says "the other **six** kinds legitimately name none" | §7.1, §10 | yes, read |
| `Audit.record`'s signature is `(kind, step, n, detail, case_ids=(), table=None)`, and `AuditEntry.__post_init__` raises `ValueError` on an undeclared kind | §7.2, §8, T1 | yes, read |
| `config.py:299-316` declares `CATEGORICAL`, `FACTOR_LEVELS` and `REFERENCE_LEVELS`, and its comment addresses Stage 6 by name, giving the deterministic-width argument this spec's §4.2 independently reaches | §4.2, §10 | yes, read |
| `PS_COVARIATES` is 9 names; `OUTCOME_COVARIATES is PS_COVARIATES` by object identity; `STANDARDISATION_COVARIATES` is `PS_COVARIATES` minus `center`; `POST_TIME_ZERO` is the two timing columns plus all eight outcomes | §4.5 D2, §6.4, §9 | yes, read |
| `test_config.py` already asserts `PS_COVARIATES` is disjoint from `POST_TIME_ZERO`, which is why §4.5's D2 is unreachable on the primary path | §4.5 | yes, read |
| `derive.constant_covariates(df, covariates=PS_COVARIATES_FULL)` detects and logs and mutates nothing; Stage 3 §7.1 assigns the deletion to Stage 6 "at the point of use" | §4.3, §10 | yes, read |
| `pyproject.toml` declares `scipy>=1.13` and `statsmodels>=0.14.2` with the comment scoping statsmodels to "unpenalised cross-checks in tests, not for any reported estimate", and no `scikit-learn` | §2, §16b | yes, read |
| `pilots/analysis.py:346-393` is `_firth_coefs` as quoted in §16; `:82` clips `e`; `:366-369` falls back to `pinv`; `:127-133` reads `IPTW_TRIM`; `:140` returns `0.0` for an all-zero arm; `:25-41` is `design` with `remove_unused_categories` and `drop_first=True`; `:951-963` complete-cases with `reset_index(drop=True)` | §15, §16 | yes, read |
| `pilots/analysis.py`'s `fit_ps` and `bootstrap` docstrings record the SAP v1.0 scikit-learn fallback firing in **16.6%** of replicates | §5.5, §16 | yes, read |
| `pilots/test_estimators.py` holds the separation benchmark (`0.5 < beta < 2.0`, comment "~0.97"), the n=4000 MLE comparison at `< 0.01`, two hand-computed ESS values, and `test_firth_ato_balance_is_near_exact_but_not_exact` asserting `0 < worst |SMD| < 0.05` | §6.3, §12.4, §12.5, §12.9, §13, §16 | yes, read |
| `test_cohort.py` pins `data.KINDS` as the declared eight | §10, T1 | yes, read |
| SAP [§7] states the estimand is indexed by the propensity model, [§8] that no library supports the required observation weights, [§11] complete-case per estimate with the denominator reported, and [§15] that unrestricted and trimmed IPTW are deliberately not used | §4.4, §6.1, §6.4, §16b | yes, read |

### 18b. The engineering review's measurements — 2026-08-14

Run by the engineering review of §20, by executing this spec's own §4.2, §5.2 and §8 code against the
landed pipeline and the workbook. Same read-only posture, aggregates only. Every row here either
corrected a statement in this document or added one it did not have.

| Claim | Where used | Verified |
|---|---|---|
| `pd.Series(pd.NA, index=idx, dtype="float64")` raises `TypeError: float() argument must be a string or a real number, not 'NAType'`. `np.nan` is the sentinel a `float64` Series takes | §3.1, §8 | yes, run — §8's canonical `fit` did not execute as written |
| `pd.get_dummies(sub, columns=factors)` **appends** the dummy columns after all non-factor columns; it does not expand a factor in place. Over `PS_COVARIATES`: `age, sex, prestroke_mrs, nihss_baseline, core_ml, tmax6_ml, atrial_fib, onset_type_unwitnessed, onset_type_wake_up, center_CHUV, center_Lugano, center_USZ` | §12.1, §7.2 | yes, run — `onset_type` is `PS_COVARIATES`' fifth entry and its dummies follow `atrial_fib` |
| On a frame with **no HUG record**, the declared `Categorical` still emits an all-zero `center_HUG`, so dropping the reference by name **does not raise `KeyError`**; the surviving centre dummies are collinear with the intercept and **F3** fires at rank ***k*−1 of *k*** — one dependency | §4.2, §12.2, §12.6, DoD-10, failure modes | yes, run — the `KeyError` claim was false in four places and in the roadmap |
| `drop_first=True` on a **declared `Categorical`** drops the first **declared** category (`center_HUG`) — identical to the by-name drop. On the **bare `string`** column it drops `center_CHUV`, leaving `center_Lugano` alone | §4.2, §12.2 | yes, run — the rebaselining hazard belongs to the missing declaration, not to `drop_first` |
| `test_config.py:592` asserts `REFERENCE_LEVELS[f] in FACTOR_LEVELS[f]` and **not** that it is `FACTOR_LEVELS[f][0]`. Both factors satisfy the stronger property today | §4.2, §10 | yes, read — the coupling the by-name drop rests on was unasserted |
| `firth` with one `nan` in **`y`** raises `FitError("step-halving exhausted …")` at iteration 1 — loud, with a message naming [§7]'s single estimator rather than the missing response | §3.1, §5.4 F6, §12.3 | yes, run |
| `firth` with a **non-binary** `y` fits, converges and returns coefficients. Nothing in `model.py` or `propensity.py` establishes that `y` is a dichotomy | §5.4 F6, §12.6, failure modes | yes, run — the one silent failure the review added a guard for |
| `.astype(float)` converts `Int64`+`pd.NA`, `Float64`+`pd.NA` and `float64`+`nan` to `nan` identically, so `design` is dtype-agnostic and `TODOS.md`'s open CTP-volume dtype question does not block this stage | §4.2, §12.3, TODOS.md | yes, run |
| A covariate list with **no factor** works unchanged: `pd.get_dummies(sub, columns=[])` is a no-op and the reference-drop comprehension is empty | §12.1 | yes, run |
| **Trust region:** on `big_first_step()`, `‖step₁‖ = 48.3` against `FIRTH_MAX_STEP = 5.0`; 8 rescales, converging in 12 iterations. **Zero** rescales on the workbook design and on `separated()` | §12.0 1b, §12.4a, §13 | yes, run |
| **Step-halving:** on `needs_halving()`, 1 halving, 0 rescales, converging in 25 iterations. **Zero** halvings on the workbook design and on `separated()` | §12.0 1b, §12.4a | yes, run |
| **`FIRTH_MAX_ITER` exhausted:** `travels_far()` reaches 200 iterations and 200 rescales without meeting either criterion — no monkeypatch required | §12.0 1b, §12.4a | yes, run |
| **`FIRTH_MAX_HALVINGS` exhausted:** `needs_halving()` under `FIRTH_MAX_HALVINGS = 1` raises at iteration 2 | §12.0 1b, §12.4a | yes, run |
| **`converged_on == "score"` fires on none of fourteen designs** — the workbook, `separated()`, `well_behaved()`, and eleven further near-separated and small-scale designs. §5.2 tests `\|Δl*\| < FIRTH_TOL` first, and near the optimum `Δl* ≈ ½·sᵀI⁻¹s`, so `max\|s\| < 1e-6` implies `Δl*` of order 1e-12 | §5.3, §5.6, §12.4a, §13, DoD-16 | yes, run |
| **`_penalised_loglik`'s `-inf` return is unreachable behind F3 and `FIRTH_WEIGHT_FLOOR`.** With the floor at 1e-10, `X'WX ≥ 1e-10·X'X` is nonsingular for full-rank `X`. Probed: at column scale 1e-3 the logdet is finite (−88.43); at 1e-30 and below `matrix_rank` gives 1 of 3 and F3 fires first | §3.1, §5.2, §5.6, §12.4a | yes, run |
| `FIRTH_WEIGHT_FLOOR`'s clip is never active: `min p(1−p)` stays above 1e-10 on every design tested, including the workbook's | §5.6, §12.4a | yes, run |
| `np.sum(np.asarray([]) ** 2)` is `0.0`, so an **empty** arm and an **all-zero** arm reach the same branch of `ess` and only a separate `v.size` check distinguishes them | §3.1, §6.2, §12.9 | yes, run |
| `collinear()` as §12.0 now specifies it is rank **3 of 4** at condition number **5.03e15**, and `np.linalg.inv` on its information matrix raises `LinAlgError` | §12.0, §12.6 | yes, run |
| `Audit.record` normalises identifiers with `tuple(sorted(str(c) for c in case_ids))` (`data.py:255`) — it **sorts without deduplicating**, and `str(pd.NA)` is `<NA>`. So a `_record_exclusion` collecting a plain `sorted()` list can never fail its own reconciliation check. `cohort.py:200-206` collects through `set()` on `notna()` and raises **before** recording, which is why Stage 5's version works | §8, §12.11 | yes, read and run |
| `cohort.py:479` already calls `absence_by_column(out, audit, [*sorted(ANALYSIS_NAMES), *sorted(DERIVED_NAMES)], …)` on the built cohort, so every [§6] covariate's per-centre absence count is already in the log | §7.4, §12.11 | yes, read |
| Audit entry counts through the five landed stages: **load 7, derive 4, classify 1, build 6**, total 18 — confirming §9's ledger exactly | §9, §12.11 | yes, run |
| `data.KINDS` is the declared eight; `test_data.py:328` and `:339` are the two kind-parametrised tests that range over it, so `model` is covered the moment it is declared — confirming §10's claim rather than assuming it | §10 | yes, read |
| `separated()` reproduces this document's numbers exactly: intercept **−5.3385**, slope **0.9706**, 9 iterations; `statsmodels.Logit` Newton `converged=False` at slope **71.5**, BFGS **27.2**. It is the one fixture the spec had written out in full, and it is the one whose numbers reproduced | §12.0, §12.5, §20 | yes, run |
| `cohort_frame()`'s golden vector reproduces to every digit §12.0.2 pins, including the `0.288470` coincidence between `HAND-2`'s weight and `COHORT-3`'s weight and propensity, and the index `(0, 1, 4, 6, 8)` | §12.0.2 | yes, run |
| Workbook: 93 records × 33 columns, `in_model` **92**, excluded record is a **control at Lugano** missing `core_ml` and `tmax6_ml`, design **11 kept of 12 candidates** with `dropped == ("center_USZ",)`, 7 iterations on the likelihood, `Σw` **13.6264** treated / **14.1103** control, ESS **30.3394 / 29.1990**, max `w` **0.964955**, rank 12 of 12, `cond(Xc)` **996.7**, `constant_covariates` empty | §4.3, §4.4, §6, §11, §12.14 | yes, re-run 2026-08-14 — every figure in §11's handover block confirmed |

### 18c. The revised code, executed — 2026-08-14

The review's own corrections were then run, because a specification that ships code is a specification
whose code has to have been executed — which is the failure §18b's first row records. Every block §20
**added or rewrote** was extracted verbatim from this document and run. Nothing below is a claim about the
implementation to come; it is a claim that the code written *here* does what the text beside it says.

| Block | Executed result |
|---|---|
| `_assert_fittable` (§5.4), all four branches and their **order** | Clean binary `y` passes. A `nan` in `X` → F4 naming the column. A `nan` in `y` → **F6**, "the response carries 1 non-finite value(s) of 6" — and it reaches F6 rather than F3, confirming the design → response → rank order. An ordinal `y` of `[0,1,2,3,4,6]` → **F6**, "takes 4 value(s) outside {0, 1}: 2, 3, 4, 6". A single stray `2` → F6 naming it. A duplicated column → **F3**, "rank 2 of 3 columns: intercept, x, x_copy" |
| F6's finiteness branch running **before** its domain branch | A `nan` response reports the non-finite message, not the outside-{0,1} one — `np.isin(nan, (0., 1.))` is False, so without the ordering a missing outcome would be reported as a domain error |
| `_probabilities` (§5.2) at a saturating linear predictor | `beta = (0, 1e6)` gives an unclipped `eta` of 5e6, clipped to `FIRTH_ETA_CLIP`, and `p` comes back **exactly 1.0** — so F5 fires, which is the property §5.5 rests on and the reason the clip must match the loop's |
| §8's `np.nan` sentinel | `pd.Series(np.nan, index=idx, dtype="float64")` builds, `dtype` is `float64`, and a partial `.loc` fill leaves `nan` in the unfilled positions |
| `_record_exclusion`'s check (§8), on the four cases that matter | One clean exclusion → names `('D',)`, no raise. A **duplicated** `case_id` among the excluded → raises, "excludes 2 record(s) and can name 1". A **missing** `case_id` among the excluded → raises, "excludes 1 record(s) and can name 0". Nobody excluded → names `()`, no raise. So the check fires on both mechanisms and does not misfire on an empty exclusion set — which the earlier plain-`sorted()` form could not do at all |
| §12.0's `well_behaved()` as specified | `max abs` difference from `statsmodels.Logit` = **0.000787** in **4 iterations** — reproducing §18b's row from the written-out generator, which is the whole point of writing it out |
| §12.0's `separated()` | intercept **−5.3385**, slope **0.9706**, **9 iterations**, likelihood route |
| §12.0's `collinear()` as specified | F3, "rank 3 of 4 columns: intercept, x1, x2, x3" |
| §12.0 1b's `big_first_step()` | **12 iterations**, likelihood, **8 rescales**, 0 halvings — matching §12.4a's claim, and isolating the trust region from halving |
| §12.0 1b's `needs_halving()` | **25 iterations**, likelihood, **1 halving**, **0 rescales** — the converse isolation |
| §12.0 1b's `travels_far()` | raises at 200 iterations, on data alone, no monkeypatch |
| `needs_halving()` under `FIRTH_MAX_HALVINGS = 1` | raises "halvings exhausted at iteration 2" |
| §12.4a's two unreachability assertions | `min p(1−p) = 0.0124` over all four fixtures, against a floor of 1e-10; `converged_on` takes exactly one value across all four, `"likelihood"` |
| **Every `python` fence in this document, parsed** | **24 of 24 parse** with `ast.parse`. Three did not before §3's rule was added: §3's signature listing (no `: ...` bodies), §6.1's one-line weight excerpt (leading indentation), and §7.1's `_HEADINGS` elision (`{..., "k": "v", ...}` mixes a set element with a dict item and is not valid Python). All three were pre-existing illustrative excerpts rather than definitions, and all three are now parseable |
| **§12.12's raise-string scan, run against this document's own `model.py` blocks** | **Passes.** 36 raise-argument strings across the nine `model.py` fences, scanned for `ivt`, `TREATMENT`, `propensity` and all eight `C.OUTCOMES` keys: zero hits. F6's reworded messages contain none of the eleven words, and the scan is implementable exactly as §12.12 describes — `ast.walk` over `ast.Raise` nodes collecting `ast.Constant` strings |

**What this does not establish.** It ran the blocks in isolation against synthetic frames and the
workbook, not `model.py` and `propensity.py` as modules with their imports, their audit calls, their
`config` constants and their AST scans — that is T2-T5's work and §12's acceptance criteria are what check
it. What it does establish is that no block in this document is of the kind §18b's first row found: called
and undefined, or written and unrunnable.

### 18d. The independent review's measurements — 2026-08-14

A second reviewer with no prior context read this document and ran its code (§20 round 3). Its findings
were then re-run independently before any of them was acted on; every row below was confirmed twice.

| Claim | Where used | Verified |
|---|---|---|
| The non-convergence message reported a likelihood movement of **exactly 0.0 on every input**, because the loop's last statement rebinds `ll_old = ll_new` before the raise reads them. Measured on the fixture that reaches it: "moved 0 against a tolerance of 1e-08" — a movement *below* the tolerance, inside a non-convergence error | §5.2, §5.2b, §12.4a | yes, run twice |
| `travels_far()` (then named `never_converges`) has a **finite, well-identified optimum**: `β ≈ (0, 2007.9)`, penalised log-likelihood **−8.6292176**, confirmed against `scipy.optimize` on the same objective. It failed only because `FIRTH_MAX_STEP × FIRTH_MAX_ITER = 1000` capped total travel below `‖β̂‖ ≈ 2008`. With `FIRTH_MAX_ITER = 1000` it converges at 403; with `FIRTH_MAX_STEP = 500` at **6** | §5.2a, §12.0 1b, §13 | yes, run twice |
| **The absolute trust radius makes convergence depend on covariate units.** Same four-record design, covariate multiplied by k: `k = 1e-3` fails at `FIRTH_MAX_ITER`; `1e-1` takes **84** iterations; `1` takes 12; `≥ 10` takes 8. A Firth estimate is equivariant under rescaling a covariate; this was not | §5.2a, §12.4a, §13 | yes, run twice |
| **The relative radius `‖step‖ / max(1, ‖β‖)` fixes it and moves nothing that is pinned.** Workbook: 7 iterations, `max|Δβ| = 0` against the absolute form, `e ∈ [0.0269724, 0.9310700]`, ESS 30.3394 / 29.1990 — identical. §12.0.2's golden vector: identical to all six pinned decimals, 11 iterations, 1 rescale. `separated()` 9 it, `well_behaved()` 4 it, `needs_halving()` 25 it / 1 halving — all unchanged. `travels_far()` now converges at 10; the scale test converges at every k from 1e-3 to 1e6 | §5.2a, §12.4a, DoD-10c, DoD-10d | yes, run |
| **The Newton-decrement radius `√(stepᵀ I step)` also fixes it and was declined.** Exactly reparameterisation-invariant, and faster — `travels_far()` at 7, the scale test at 4-8 — but it perturbs the workbook coefficients at **6.04e-06**, moves §12.0.2's golden vector in the **fifth decimal** (`0.929035` against `0.929038`), and binds 4 times on `well_behaved()`, a well-conditioned 4000-record design where no safeguard should fire | §5.2a | yes, run |
| **`design` raised `KeyError: 'case_id'` on any frame carrying a declared factor and no identifier**, because D3 indexed `df.loc[..., "case_id"]` before its own `if len(off)` guard. Measured on a clean two-row frame. A frame with no `CATEGORICAL` covariate was unaffected, so the failure was invisible to any test that happened to pass a cohort | §4.5 D3, §9, §12.1 | yes, run twice |
| **§12.8's "the other four probabilities are unchanged to 1e-12" is false.** Blanking `HAND-1`'s `core_ml` moves the four survivors by up to **0.047763** — correct behaviour, since a 4-record fit is a different fit from a 5-record one | §12.8 | yes, run twice |
| **`cohort_frame()`'s cohort is completely separated, so §12.10's MLE companion cannot hold on it.** `statsmodels.Logit` Newton: `converged=False`, `max|coef| = 129.7`, `Σw₁ − Σw₀ = −4.4e-11` — "passing" 1e-10 only because both sums underflowed. BFGS: converged to a different point, `max|coef| = 167.7`, `Σw₁ − Σw₀ = 3.9e-07`, failing it. On the **workbook** design the MLE converges (`max|coef| 4.55`) and `Σw₁ − Σw₀ = 1.8e-15` | §12.10, §6.3 | yes, run twice |
| **The reference dummies appear in neither `X.columns` nor `dropped`**, so the `design_matrix` table named no baseline. Measured on `cohort_frame()`: kept `core_ml, tmax6_ml, center_CHUV`; dropped nine; `center_HUG` and `onset_type_witnessed` in neither | §7.2, §7.5, §12.11 | yes, run twice |
| **§12.4a's `-inf` companion does not reach `-inf`.** With the rescale removed, `big_first_step()`'s first candidate step is finite and **uphill** — `ll = −4.88047` against `ll_old = −6.04068` — and is accepted immediately with zero halvings | §12.4a, §5.2 | yes, run |
| **F3 fires at rank *k*−1, not *k*−2**, when the reference level is absent: the surviving dummies sum to 1 on every row, which is one dependency. The earlier "3 of 5" came from a toy frame where `onset_type_wake_up` was accidentally identical to `center_Lugano` | DoD-10, failure modes | yes, run |
| Four of eight `pilots/analysis.py` line references were wrong: `drop_first=True` is at **:40** (not :38), `remove_unused_categories` at **:39** (not :34), the `pinv` fallback at **:366-369** (not :355-358 / :357), `covars or C.PS_COVARIATES` at **:953** (not :964, which is blank), and the `e` clip at **:82** (§3.1 said :81, four other citations said :82). Correct as cited: `:25-41`, `:346-393`, `:140`, `:127-133`, `:951-963` | §15, §16, §18 | yes, read |
| SAP §7 line 145 reads *"Report the effective sample size per arm and **describe the weighted population**. Because the ATO population is defined statistically, **state explicitly who it comprises**."* §7.2 specified the ESS and a conditionality sentence and no description | §6.4, §7.6 | yes, read |
| ~~Over 400 randomised designs the estimator produced **no false convergence**, worst gap 4.8e-9~~ — **WITHDRAWN.** True of the designs drawn, which were not completely separated. On separated designs the loop converges to a local maximum on 4 of 600; see §5.2c and §18e | §5.2c, §16b | withdrawn 2026-08-14 |
| **§6.3's explanation is right, not just its measurement**: `Σw₁ − Σw₀ = −0.483902` against `Σ h(e−0.5) = −0.483891`, the 1.15e-5 remainder being the unsolved intercept score | §6.3 | yes, run |
| **At the returned `beta` the workbook fit's `max|modified score|` is 6.1e-3**, not near zero — convergence is an absolute likelihood tolerance on a flat surface. So §16b's "agrees with `firthlogist` to 2.7e-15" is two implementations **stopping in the same place**, not agreement on the optimum; the `scipy` row's 2.7e-05 on the coefficients is the honest number | §16b, §5.3 | yes, run |

### 18e. Round 4's measurements — 2026-08-14

A second independent reviewer, again with no prior context, read the revised document and ran its code.
Its findings were re-run before any was acted on. §5.2c's non-concavity finding is the substantive one;
the rest are the document disagreeing with itself.

| Claim | Where used | Verified |
|---|---|---|
| **`½log|I(b)|` is not concave.** Numeric `d²/ds²` along a slope ray on a 1-D separated design: **−0.66** at slope 0.5, then **+0.20**, **+0.019**, **+0.0047** at slopes 2, 5, 20. So the penalised likelihood is not unimodal and Newton from `beta = 0` finds a stationary point, not certainly the maximiser | §5.2c, §16b | yes, run |
| **4 of 600** synthetic completely-separated designs (standardised covariates, n∈[20,60), k∈[1,4)) converge below what other starts of the same loop reach. Worst gap **0.29 nats**, `max|Δe|` **0.40**. On all four: `e` strictly interior (max 0.99999997), `converged_on = "likelihood"`, **0 rescales, 0 halvings** — F5, F3 and all three §3.2 counters clean | §5.2c | yes, run twice |
| **1 of 2000** centre-stratified bootstrap replicates of the actual 92-record cohort, gap **1.0e-02** nats; a second independent 2000-replicate sweep also found 1. Median gap **4.2e-11**, 99th percentile **1.2e-09**, 0 replicates dropped for rank, arm or convergence | §5.2c, §9, §13 | yes, run twice |
| **The v7 point fit is at the global optimum** — 29-start search, gap 1.3e-10, `max|Δβ| = 1.7e-07`. Every reported figure stands | §5.2c, §11, §18 | yes, run |
| **Four candidate detectors all fail**: 5 structured starts 0 of 1; 48 coordinate-wise `±c·e_j` 0 of 1; 12 starts seeded once from `C.SEED` 0 of 1; 12-start multistart as a corrector recovered 1 of 4 synthetic cases, with 12 of 12 starts agreeing on the wrong point in the other three | §5.2c, §9 | yes, run |
| **A 12-start multistart would move pinned numbers**: `cohort_frame()`'s own fit sits 3.6e-08 below the best found, shifting §12.0.2's golden vector in the 4th decimal (`HAND-5` 0.702237 → 0.702268). This is why §5.2c detects rather than corrects | §5.2c, §12.0.2 | yes, run |
| `propensity.fit` on `test_cohort.built()` fires **1 rescale and 2 halvings** — so `cohort_frame()`'s cohort is not among the fixtures that reach neither, which an earlier §12.4a sentence claimed | §12.4a, §13 | yes, run |
| Removing the trust-region rescale makes `big_first_step()` converge in **8** iterations rather than more, and moves the answer by **4.8e-06** — so the natural companion test is false in both halves | §12.4a | yes, run |
| The reference dummies appear in neither `X.columns` nor `dropped`, so a `kept + dropped` table names no baseline — confirming §7.5's premise against the revised `design` | §7.5, §12.11 | yes, run |
| On `cohort_frame()`'s cohort, `PS_COVARIATES_FULL`'s four extra covariates are **constant and all dropped**, so the design is identical to `PS_COVARIATES`' — the "four more columns appear" assertion holds only on the workbook | §12.1 | yes, run |
| The Firth kernel is right: the numeric gradient of `_penalised_loglik` matches the modified score `X'(y − p + h(0.5 − p))` to six digits, so the penalty, the score and the hat diagonal are all correct | §5.1, §5.2 | yes, run |
| Every workbook figure, §12.0.2's golden vector, §5.3's `‖I⁻¹‖` route table, §12.4a's fixture measurements, §12.8's respecified assertions, §12.10's data-gated MLE companion and §12.13's unrestricted-frame fit all reproduce exactly | throughout | yes, run |

### 18f. Round 5's measurements — 2026-08-16

The audit that followed round 4, checking the four of its findings that had not been acted on.

| Claim | Where used | Verified |
|---|---|---|
| **A missing factor value is encoded identically to a real reference-level record.** With one record's `center` set to `pd.NA` and `complete_cases` bypassed: the design carries **no nan at all**, and that record's centre dummies are `[0. 0. 0.]` — element-for-element the same as a real HUG record's. F4 has nothing to test, D3 excludes it by `& notna()`, F5 sees a finite interior probability | §4.5 D4, §12.6, failure modes | yes, run |
| **`FIRTH_WEIGHT_FLOOR` is reachable.** On 800 synthetic completely-separated designs the floor binds on **68 (8.5%)**, and the smallest `p(1−p)` observed is **exactly 0.0**. It is never active on any `test_model.py` fixture (0.0124) or on the cohort design (0.0262) — so §12.4a's assertion is a statement about those inputs, not about the estimator | §5.6, §12.4a, §5.2c | yes, run |
| `_design_table(X, dropped)` cannot build §7.5's reference rows from its two arguments: the reference level is by construction in neither `X.columns` nor `dropped`, so it must be read from `REFERENCE_LEVELS`, and inferring factor membership from `{factor}_{level}` name prefixes misreads any linear covariate that happens to share a prefix | §7.5, §8 | yes, run |
| The seven tolerances were named `PS_*` and read by `model.py`, which §0.1 defines as knowing nothing about the exposure; §12.12's attribute scan forbids `C.TREATMENT` and `C.PS_COVARIATES` and would have passed all seven | §5.6, §10, §12.12 | yes, read |

### 18g. The implementation's measurements — 2026-08-16

T1-T7 built against this document. Every figure §11, §12.14 and §18 pin **reproduced exactly** and is
not re-listed; what follows is only what this document did **not** predict, plus the two Definition of
done items that are **not satisfied**. Four rows correct a claim here, and three of the four are of the
same kind: a measurement §18 records that the code as specified does not produce.

| Claim | Where used | Verified |
|---|---|---|
| **The whole of §11's handover block reproduces**: cohort 93 x 33, `in_model` **92**, the excluded record a control at Lugano missing both CTP volumes, design **11 kept of 12** with `dropped == ("center_USZ",)`, rank 12 of 12, **7 iterations** on the likelihood, `e ∈ [0.0269724, 0.9310700]`, `Σw` **13.6264** / **14.1103**, ESS **30.3394 / 29.1990**, max `w` **0.964955**, `constant_covariates` empty, 3.25 treated per parameter | §11, §12.14 | yes, run |
| **§12.0's and §12.0 1b's every fixture measurement reproduces**: `well_behaved` 0.000787 in 4 iterations, `separated` −5.3385 / 0.9706 in 9, `collinear` rank 3 of 4, `big_first_step` 9 iterations with **2 rescales** and `‖step₁‖ = 48.3`, `needs_halving` 25 iterations with 1 halving and 0 rescales, `travels_far` converging at 10 on the **score** route, `min p(1−p)` = 0.0124, and §5.3's `‖I⁻¹‖` route table to its stated orders | §12.0, §12.4, §12.4a, §12.5, §5.3 | yes, run |
| **§12.0.2's golden vector reproduces to all six pinned decimals** (worst deviation 4.06e-07), at 11 iterations and 1 rescale; DoD-10d's workbook comparison against the absolute radius gives **`max\|Δβ\| = 0`** at 7 iterations, and DoD-10c's scale test reproduces §18d exactly — the absolute radius **fails** at 1e-3 and takes **84** iterations at 1e-1 where the relative one converges at every scale | §12.0.2, DoD-10c, DoD-10d | yes, run |
| **§9's audit ledger reproduces exactly**: load 7 / derive 4 / classify 1 / build 6 / fit 4 = 22 entries, the four `model` entries in `fit`'s declared order, rendering under `## Fitted models` between `## Cohort construction` and `## Structural non-applicability`, byte-identical across `PYTHONHASHSEED` 0 and 1 | §9, §12.11, DoD-5 | yes, run |
| **DoD-12's pair discriminates**: `\|Σw₁ − Σw₀\|` is **0.483902** under Firth and **1.78e-15** under a `statsmodels` MLE score, so substituting the estimator fails the inequality assertion | §6.3, §12.10, DoD-12 | yes, run |
| **§8's `overlap_weights` record does not execute as written.** It concatenates a **6-column** weights table with `_weighted_population_table(...)[1:]`, whose rows are **5** columns, so `data._md_table` raises `IndexError` on `rows[0]`'s width; and the `[1:]` discards the population header, which is what names its four arm-by-weighting columns — so the second block would render under `arm \| n \| sum w \| ESS \| max w \| weight share`. §7.2 asks one entry to carry two tables and `AuditEntry` holds one. **Fixed** by keeping the population header as a labelled body row and padding the short rows with `_design_table`'s not-applicable dash; nothing is dropped, nothing is unlabelled, and the rendering stays deterministic. §18c parsed this fence but states it did not run the blocks as modules with their audit calls, which is exactly where it surfaces | §7.2, §7.6, §8 | yes, run — **§8 corrected** |
| **§3.1 and §18's "a nan in `X` terminates normally" is wrong, and the correct behaviour is worse for Stage 10 rather than milder.** Measured: `slogdet` of a matrix containing nan returns `logabsdet = nan`, so `_penalised_loglik`'s own `isfinite` branch returns **−inf**; `-inf >= -inf` accepts every candidate step with **zero halvings**; `moved` is `abs(-inf − −inf)` = nan and the modified score is nan, so **neither convergence route can ever fire** and the loop raises *no convergence in 200 iterations*. It does not return all-nan coefficients and it does not terminate normally. The harm F4 prevents is therefore not "a `[nan, nan]` interval" but a **misdiagnosis**: a missing covariate reported as non-convergence, which [§10] catches to **drop and count the replicate** — the data defect absorbed into a sparse-replicate failure count. Confirmed identical under `pilots/analysis.py`'s kernel, which shares the branch | §3.1, §12.3, §18, DoD-6 | yes, run — **§3.1 corrected** |
| **§18's `firthlogist` agreement of 2.665e-15 is not reproducible, from either kernel.** Measured on the cohort design: **3.659e-06**. This implementation and `pilots/analysis.py`'s `_firth_coefs` agree to **`max\|Δβ\| = 0`**, so it is not a property of this implementation; the pilots' kernel disagrees with `firthlogist` by the same 3.659e-06. The cause is structural and is the row below. The two **do** agree on the penalised log-likelihood to **2e-10** — −19.0637927357 against −19.0637927355, reproducing §16b's recorded value — and `scipy.optimize` reaches −19.0637927356 from an independent direction. So §12.4's oracle asserts the **objective**, which is the form §16b itself calls honest for the `scipy` check, with a 1e-4 coefficient bound beside it | §16b, §18, §12.4, DoD-11 | yes, run — **§18 corrected** |
| **Definition of done 15a is DISCHARGED, and §5.3's attribution is STRUCK.** `logistf` 1.26.1's CRAN tarball was downloaded and read; no working R is needed for this. **(a) Defaults** (`R/logistf.control.R:41`): `maxit=25, maxhs=0, maxstep=5, lconv=1e-5, gconv=1e-5, xconv=1e-5` — none of them §5.6's values, and `maxhs=0` means no step-halving at all by default. **`xconv` IS a coefficient-step test**, the criterion §5.3 deliberately excludes. **(b) Order** (`src/logistf.c:339-341`): one conjunction, `maxabs(delta) <= xconv && maxabs(Ustar) < gconv && loglik_change < lconv`, OR'd with `iter >= maxit`. **(c)** `gconv` compares the **raw** modified score — `Ustar` is built by `XtY(x, w, Ustar, ...)` at `:332` and is not information-scaled. **So `logistf` requires all three CONJUNCTIVELY, including the step test, where §5.3 requires either of two DISJUNCTIVELY and excludes the third.** The rule is not `logistf`'s, in structure or in tolerance. Corroborated independently by `firthlogist`, a port, whose only test is `norm(coef_new - coef) < tol` (`firthlogist.py:307`) — `xconv` alone — with an absolute infinity-norm step bound and `lstsq` in place of `inv`. **§5.3's argument survives untouched and never depended on the attribution**, exactly as this document anticipated: excluding the step norm is right because the penalised surface is genuinely flat under near-separation, and [§10] drops failed replicates so the ones dropped would be the sparse ones. The sentence is struck in `config.py`'s §5.6 block, in `model.py`'s convergence comment, and pinned by a test | §5.3, §5.6, §13, §16b.2, DoD-15a | yes, read — **§5.3 corrected** |
| **Definition of done 16 is SATISFIED: the R gate has been opened and passed.** `logistf` **1.26.1** and `PSweight` **2.1.2** installed on R **4.5.0**; all nine tests in `test_reference_r.py` pass. This is the first independent — not transitive — check on the estimator this stage has had | §16b.1, §16b.2, §13, DoD-16 | yes, run |
| **`logistf` reaches the same penalised optimum as this implementation to 9.96e-10** on the workbook design: −19.0637927357 against −19.0637927367, with `max\|Δβ\| = 2.275e-05`. **This is the strongest evidence available for the estimator and it is worth stating precisely what it is evidence of.** `firthlogist` is a *port* of `logistf`, so its agreement was transitive; `logistf` is the reference implementation itself, and it stops by a **different rule** (xconv AND gconv AND lconv, see the DoD-15a row) yet lands on the same optimum. The coefficient spread of 2.3e-05 is the flat ridge §5.3 describes, and it sits between `firthlogist`'s 3.66e-06 and `scipy`'s 2.42e-05 — three independent optimisers strung along the same ridge, agreeing on the objective to 1e-9 or better. **The arithmetic — penalty, modified score and hat diagonal — is now validated against the reference implementation rather than against a port of it** | §16b, §16b.2, §5.3, DoD-16 | yes, run |
| **`PSweight`'s ATO weights reproduce the identity exactly, and its propensity score is a different score — as expected.** Its `e` satisfies `w = 1 − e` treated / `w = e` control to 1e-12, so the ESTIMAND agrees; `spearman(e_ours, e_PSweight) = 0.9997` and `max\|Δe\| = 0.0499`, because `PSweight` fits an unpenalised `glm` where [§7] prescribes Firth. Kish's formula on **its** scores gives 28.20 / 24.81 against our 30.34 / 29.19 — a difference driven by the SCORES, not by the ESS formula, which is why §16b.2's instruction to compare `e` and `w` first and the ESS separately is the right order. **And §16b.2's last open question is now CLOSED: the ESS `PSweight` reports IS Kish's.** Its `SumStat` object carries an `ess` field, and its overlap column equals `(Σw)²/Σw²` on its own weights to every digit it prints — measured on synthetic data, so the check is about the convention and not about a patient number. There is therefore **no convention difference to caveat**, and the whole 28.20 / 24.81 against 30.34 / 29.19 gap is the propensity model: [§7]'s Firth score against `PSweight`'s unpenalised `glm`. Firth shrinks the linear predictor toward zero, which pulls `e` toward 0.5 and *raises* the overlap weights, so a slightly larger effective sample is the expected direction rather than a surprise. **§16b.2's stated worry — "or report it after its own normalisation of the weights" — could not have produced a difference in any case, and that is worth recording as a small correction to the section: Kish is scale-invariant, since `(Σcw)²/Σ(cw)² = (Σw)²/Σw²`, so only a different FORMULA could ever have differed.** The committed script now emits the package's own `ess` and the test asserts it against `propensity.ess` on the same weights, so the convention cannot drift unnoticed | §6.2, §16b.2, DoD-16 | yes, run — **§16b.2 closed** |
| **`PSweight` independently confirms §6.3's explanation**, which no Python oracle could: on its unpenalised `glm` score the arm weight sums are **exactly equal** (47.3998 against 47.3998 on the synthetic frame), because that is the score equation — while [§7]'s Firth score gives 13.6264 against 14.1103. §12.10 asserts that pair against `statsmodels`; it now also holds against the reference implementation of the estimand, by the authors of the method | §6.3, §12.10 | yes, run |
| **The obstacle to DoD-16 was environmental, and two intermediate diagnoses were wrong.** R's startup runs `system("uname -a", intern = TRUE)`; PATH carried a second toolchain prefix whose `sh`/`uname` are linked against a foreign libc, so the call segfaulted (status 139) inside R's forked child, `utils` never loaded, and `install.packages` did not exist. **It is not the sandbox** — reproduced with the sandbox disabled and under `--vanilla --no-init-file`, which is what an intermediate draft of this row claimed §18 had got wrong — **and it is not a broken `r-base-core`**, which is what that same draft claimed instead. Prepending `/usr/bin:/bin` to PATH fixes it entirely. Two further prerequisites, both about libraries: a glibc-matched `cmake` for `nloptr` → `lme4` → `PSweight`, and `R_LIBS_SITE` excluding `/usr/local/lib/R/site-library`, which holds **75 packages built under R 3.6.3** that abort dependent builds with "installed before R 4.0.0". `_r_environment` in `test_reference_r.py` encodes the PATH fix and passes the library variables through, so **the gate cannot silently skip for this reason again** — which is the hazard §16b.2 names, and which did occur: the probe found both packages while every script died in R's startup, until `_run` was given the same environment as `_r_has` | §16b.2, §18, DoD-16 | yes, run — **§18 corrected** |
| **The `reference` group cannot be declared without moving the primary environment's numpy.** uv resolves dependency groups against the project's own dependencies into one lock, `firthlogist` 0.5.0 declares `numpy>=1.22.4,<2.0.0`, and declaring the group therefore resolved the **base** lock from numpy 2.0.2 to **1.26.4** and scipy 1.18.0 to 1.17.1. Constraining numpy to `>=2.0,<2.1` makes the resolution unsatisfiable outright — the oracle genuinely cannot run on numpy 2. Both versions are inside the declared `numpy>=1.26,<2.1`, so the pin is not violated, and **every figure above was re-measured at 1.26.4 and reproduces**, so the downgrade is inert here rather than assumed so. It also needs `override-dependencies = ["tabulate>=0.9"]`, because `firthlogist` caps `tabulate<0.9` against this project's `>=0.9` — safe for a measured reason, since it imports `tabulate` only for its `summary()` printer and the oracle calls `.fit()` | §2, §16b, §13 | yes, run |
| **DoD-9's weights are 1e-8 rather than "exactly 0.0", and the point stands.** With F5 replaced by `np.clip(e, 1e-8, 1 - 1e-8)` on a design whose linear predictor saturates `FIRTH_ETA_CLIP`, the weights become 1.0e-08 and the treated arm's ESS is a finite **2.0** — which reads as two whole observations for an arm the fit placed entirely on the boundary. The clip converts a degenerate fit into numbers every downstream sum accepts, which is §5.5's argument; the literal `0.0` is what an unclipped `1 - 1.0` gives | §5.5, DoD-9 | yes, run |
| pandas **2.3.3**, numpy **1.26.4**, scipy **1.17.1**, statsmodels **0.14.6**, `firthlogist` **0.5.0** under scikit-learn **1.5.2**, tabulate **0.10.0**, Python 3.12 | §2, §18, DoD-11 | yes, read |

**Definition of done, item by item.** 1 ✅, 2 ✅ (61 skips, nothing errors at collection), 3 ✅ (3 skips),
4 ✅, 5 ✅, 6 ✅ *(as corrected above)*, 7 ✅, 8 ✅, 9 ✅, 10 ✅, 10a ✅, 10b ✅, 10c ✅, 10d ✅, 11 ✅
*(at the objective; see the correction)*, 12 ✅, 13 ✅, 14 ✅, 15 ✅, 15a ✅, **16 ✅ — the gate was
opened on 2026-08-16 against `logistf` 1.26.1 and `PSweight` 2.1.2, and passed.**

**All sixteen items hold, and nothing in §16b.2 remains open, so Stage 6 is complete on this
document's own terms.** The ESS-convention question DoD-16 named was the last of them and is answered:
`PSweight` reports Kish.

## 19. What this spec changed elsewhere

**Five** entries, all landing with this document rather than with the implementation, because a
specification that contradicts the roadmap is worse than no specification. Precedent: Stages 3, 4 and
5 landed their roadmap edits with the document. The fifth was added by the engineering review of §20 and
is a **correction** rather than an addition: the roadmap carried the same false `KeyError` claim this
document did.

| # | Change | Files | The question it settles |
|---|---|---|---|
| 1 | Roadmap Stage 6 gains its `**Spec:**` line | `implementation_roadmap.md` | Stage 6 was the next stage without one |
| 2 | Roadmap Stage 6's **Accept when** gains the **complete-case rule**: the propensity model is fitted on the [§6]-complete subset, `e` and `w` are returned aligned to the whole cohort with the missingness reimposed, and the excluded patients are named in the log — one control on v7 | `implementation_roadmap.md`, §4.4 here | Whether a patient missing a covariate leaves the analysis silently. [§11] said complete-case per estimate and no stage had yet had an estimate; this is the first, and the roadmap did not say what the fit does with the gap |
| 3 | Roadmap Stage 6's **Accept when** gains the two **silent-failure** assertions — a non-finite design and an unlisted factor level — as runtime raises rather than test-time checks | `implementation_roadmap.md`, §3.1, §4.5, §5.4 here | Whether "the returned probabilities are finite and in (0, 1)" is enough. It is not: measured, a `nan` in the design produces all-`nan` output through a normally terminating loop, and an unlisted level produces a *finite* probability computed at the reference category. Neither is caught by a check on the output |
| 4 | Roadmap Stage 6's **Accept when** gains the `model` audit kind and the requirement that the ESS, the weighted-population description and [§7]'s conditionality statement live in the reproducible log | `implementation_roadmap.md`, §7 here | The roadmap said "record in the output" without saying where. `out/tables` is Stage 14's and does not exist yet, so "the output" would have been a return value the reporting layer might not print |
| 5 | Roadmap Stage 6's design-matrix bullet is **corrected**: "a replicate containing no HUG patient would otherwise rebaseline silently" is false, because the declared level set makes `drop_first` pick the reference too. The bullet now gives the real reason (independence from the reference-is-first coupling, now asserted in `test_config.py`) and names **F3** as what fires when the reference level is absent | `implementation_roadmap.md`, §4.2 here | Whether the roadmap's own justification for a landed decision is true. It was not — measured — and the decision it justifies is still right, which is the combination most likely to survive unexamined. Added by the engineering review of §20, finding 2 |

**What this spec has had, and what it has not.** Every numeric claim in §4, §5, §6, §11, §12 and §18
was produced by running code against the workbook, against `cohort_frame()`, against synthetic data and
against `firthlogist` — and the §16b investigation was run rather than reasoned, which is the only
reason its conclusion is trustworthy: three of the four candidate libraries fail for reasons no amount
of reading would have surfaced.

**It has now had both, and §20 is the record.** The engineering review is §20's rounds 1 and 2; the
independent second reader this section asked for is §20.1, obtained on 2026-08-14 after a first attempt
failed on authentication. A **second** independent reviewer ran after those corrections landed. The two found **seventeen defects
each**, none overlapping either the first seventeen or each other — so this paragraph's prediction was
right and its implied magnitude was badly low. Anyone tempted to treat the document as settled should
read that sentence twice: **four passes, all of which ran code, produced four disjoint defect sets of the
same size**, and the fourth still found one defect in the estimator itself.

**The Definition of done's items 6-12 remain instructions to break the code and watch**, and the review
did not weaken them — it added 10a and 10b, and it corrected 10, which told the implementer to watch for
a `KeyError` that cannot occur.

## 20. What the reviews changed — 2026-08-14 to 2026-08-16

Reviewed by executing this document's own §4.2, §5.2 and §8 code against the landed Stages 1-5, the
workbook and `cohort_frame()`, in the manner Stage 4's and Stage 5's reviews used. §18b is the
measurement record and §18c records the revised code being run in turn. **Seventeen findings, all of them
acted on**, and the yield of the exercise was in fact what §19 predicted rather than zero: three
statements in this document were false against pandas 2.3.3, two specified acceptance tests could not
have passed, and one silent failure had no guard.

**Findings 1-12 came from reading and running the original; 13-17 came from writing the fixes** — which
is worth separating, because it is the same pattern §19 noted about Stages 4 and 5: what a review misses
gets caught by implementation. Five things only became visible once the corrections existed as code; one
of them (13) is a conflict the review itself introduced, and one (17) was found by running the check that
13 produced. That last chain is the argument for §18c existing at all: a check written into a spec is
worth more once it has been run against the spec.

The three items §19 asked a reviewer to weigh, answered first:

- **§0.1's two-module split: keep it.** The seam is in the right place — `model.py` imports `config`
  alone, takes no `Audit`, and §12.12 asserts it names neither the exposure nor the covariate list, so
  Stage 9 can fit an outcome model without importing a module named for the treatment. The review's own
  F6 finding is evidence *for* the split rather than against it: the boundary is exactly what made a
  missing response check visible as a gap, and every guard added at that boundary is one Stages 8, 9 and
  12 inherit.
- **§7.3's line on where a coefficient vector may live: keep it.** It costs nothing that is not already
  paid — `out/` is gitignored in full — and the golden vector of §12.0.2 gives the regression pin the
  coefficients would otherwise have provided, on synthetic records, with no patient in git.
- **§4.4's `<NA>` weight: keep the decision, fix the word.** The sentinel is `np.nan`, not `pd.NA`
  (§18b), and `float64` is right — a masked `Float64` would reintroduce §3.1's first hazard on the one
  column most likely to be summed. On "can every downstream stage be trusted to range over `in_model`":
  the answer is that it cannot be *trusted*, and §9 now states per stage what each one owes. `nan` is the
  safe direction only while nothing fills it, and nothing in this stage can enforce that.

| # | Finding | Fix | Where |
|---|---|---|---|
| 1 | `pd.Series(pd.NA, index=df.index, dtype="float64")` raises `TypeError`. §8's canonical `fit` — the block a reader is told to trace through one screen — did not run | `np.nan`, with §3.1 carrying the measurement and §8's bullet rewritten around `(nan, in_model)` rather than around a sentinel the column cannot hold | §3.1, §8, §18b |
| 2 | The by-name reference drop **cannot** raise `KeyError`: the declared `Categorical` always emits the reference column. And `drop_first=True` on a declared `Categorical` drops the *first declared* category, so it is **equivalent** to by-name, not divergent. The rebaselining hazard belongs to the bare column | §4.2's bullet rewritten around the real reason — independence from `REFERENCE_LEVELS[c] == FACTOR_LEVELS[c][0]`, which nothing asserted — plus that assertion added to `test_config.py`, and F3 named as the actual failure in all five places that said `KeyError` | §4.2, §10, §12.2, DoD-10, failure modes, roadmap |
| 3 | `firth` validated `X` and never `y`. `_assert_fittable` took `y` and never read it. A `nan` response raises with a message about [§7]'s single estimator; an **ordinal** response fits silently and returns coefficients | **F6** added: `y` finite and in `{0, 1}`, ordered after F4 and before F3. §5.4 gains the argument, §12.3 and §12.6 the tests, and the `_assert_fittable` comment a third silent-failure line | §5.4, §12.3, §12.6, §17, DoD-10a |
| 4 | `_probabilities` was called twice in §5.2 and defined nowhere; §3's private list named a `_step` that does not exist and omitted `_probabilities` and `_assert_fittable` | `_probabilities` written out, with the reason its clip must match the loop's; the private list corrected and the absence of `_step` explained | §3, §5.2 |
| 5 | §6.2's `ess` raised a bare `FitError`, which is not in `propensity.py`'s namespace — §3's import block is `import model` | `model.FitError`, and the docstring says why it is qualified | §6.2 |
| 6 | `_record_exclusion`'s reconciliation check was dead code: a plain `sorted()` list makes `len(case_ids) == n` by construction, since `Audit.record` sorts without deduplicating. §12.11's acceptance test could not have passed | Mirrored `cohort.py:200-206` exactly — `set()`, `notna()`, and the raise **before** `audit.record` — with §12.11 gaining the companion that shows the earlier form unfirable | §8, §12.11 |
| 7 | §7.2's completeness table re-rendered `data.absence_by_column`'s per-centre columns a fourth time, over a frame `cohort.py:479` already covers | Narrowed to `covariate`, `n_absent`, `excluded_by`; §7.4 written to record the decision and the declined alternative; §12.11 asserts the two tables agree | §7.2, §7.4, §12.11 |
| 8 | Five of §5.2's six "load-bearing" numerical details were reached by **no** specified input — measured at zero trust-region rescales, zero halvings and zero exhausted loops on the workbook and on `separated()` — while §12's coverage map read "No branch is untested" | Three reachable branches got fixtures (`big_first_step`, `needs_halving`, `travels_far`) and a fourth a monkeypatched constant, all in §12.4a. **Two are unreachable by construction** and are recorded as such with assertions that they stay so | §5.6, §12.0 1b, §12.4a, coverage map, §13 |
| 9 | `converged_on == "score"` did not fire on any of fourteen designs, and the score is evaluated at the pre-step iterate | The pre-step observation stands and is recorded in §5.3. **The unreachability conclusion drawn from it was wrong and is retracted in round 3** — see finding 34. It is left in this table rather than deleted, because a review that records only its correct findings is not a record | §5.3 |
| 10 | `well_behaved()` and `collinear()` were signatures with no data-generating process, while §18 pinned a measured `0.00116` from one of them. `separated()` was fully written out — and it is the one whose numbers reproduced exactly | Both written out in §12.0 and re-measured: `0.000787` in 4 iterations, and `collinear()` at rank 3 of 4, condition number 5.03e15. §18's row replaced with the date and the reason | §12.0, §12.4, §18 |
| 11 | §12.1's column-order literal was specified as "`covariates`' order with each factor expanded in place". `get_dummies` **appends** the dummies, so the literal would have failed — and `Fit.columns` is what maps a coefficient to a name in the only table the coefficients ever appear in | The measured order written out, with the consequence stated: a wrong order is twelve right numbers under eleven wrong labels | §12.1, §7.2, §18b |
| 12 | The document called the R oracle "an open item rather than a blocker" in §13 and §16b.2 while Definition of done 16 made it a completion gate, and §18 records that neither R package is installed. It argued both sides | **Resolved in favour of the gate**, by decision of 2026-08-14. The softening language is struck from both sections, T6 moves to P1, and installing `logistf` and `PSweight` becomes T6's work. Finding 9 is the reason the gate is worth its cost | §13, §16b.2, §17 T6, DoD-16 |

| 13 | **A conflict the review introduced.** F6's first draft read "use `propensity.fit`'s F1, or your outcome's own mask" and "`mrs_90d` has seven levels and would fit silently" — putting the exposure's name and a specific outcome's name inside the module §0.1 defines as outcome-agnostic, which §12.12's scan exists to prevent. An attribute scan would have passed it | F6's messages reworded to name neither; the concrete examples moved into §5.4's prose, where they are more useful. §12.12 gains a **second scan over the string literals inside `raise` statements**, scoped to raises rather than to the whole file **because §17 requires the module docstring to name `propensity.py`** — the one sentence that stops someone moving `design` there for tidiness | §5.4, §12.12, §17 |
| 14 | §2's cost model counted the seven inversions and seven `einsum`s but not F3's `matrix_rank` SVD, and said nothing about what the stage costs when Stage 10 calls it 2000 times — which is the only context where the number could matter and the context in which someone would reach for a cache | §2 gains the SVD and the Stage 10 arithmetic (~28,000 inversions, ~4,000 rank SVDs), and the standing instruction that nothing here is ever a performance consideration — so no later stage adds an `lru_cache` on `design` or vectorises `firth` over replicates on performance grounds | §2, §15 |
| 15 | §12.5 attributed the 0.9706 separation slope to `logistf` — the same unverified attribution §5.6's constants carried, in a second place the review's first pass missed | Withdrawn on the same terms and for the same reason: `logistf` has not been run, `firthlogist` is a port of it so agreement is transitive, and DoD-16 either earns the attribution or strikes it. The measured intercept and iteration count added in its place | §12.5, §5.6, DoD-16 |
| 16 | The revised blocks had not themselves been executed — which is precisely the defect finding 1 was | §18c: every block §20 added or rewrote, extracted verbatim and run. All four `_assert_fittable` branches in order, F6's finiteness-before-domain ordering, `_probabilities` reaching exactly 1.0 so F5 fires, the `np.nan` sentinel, `_record_exclusion` firing on **both** a duplicated and a missing `case_id` and not misfiring on an empty exclusion set, and all seven §12.0 fixtures reproducing the numbers §12.4a and §18b claim for them | §18c |
| 17 | **Three `python` fences in this document did not parse** — §3's signature listing, §6.1's indented one-line excerpt, and §7.1's `{..., "k": "v", ...}` elision, which mixes a set element with a dict item and is not valid Python. Found by running finding 13's scan, which had to parse every fence to search it. A sole-source document an implementer copies out of should not contain a fence nobody can run | All three made parseable, and §3 states the rule: **every `python` fence in this document is valid Python**, elisions included, which is why the signature listing carries `: ...` bodies and §7.1's elision is comments rather than `...`. §18c asserts 24 of 24 | §3, §6.1, §7.1, §18c |

**Smaller corrections, verified and applied without a separate finding:** §1 said "six declared
numerics" where §5.6, §10 and the coverage map say seven; §0's diagram said `design` returns
`(X, kept, dropped)`; §18's `e` maximum was `0.93110` against a measured `0.9310700` and §11's handover
block rounded it to four places; §10's `test_cohort.py` row named the literal pin but not the rename §1
counts; §1's `test_model.py` row did not mention §12.4a or the hard-fit family; `ess` on an **empty** arm
reached the all-zero-arm message, which is a different diagnosis; `design` over a covariate list with no
factor at all was untested (it works — `get_dummies(columns=[])` is a no-op); the coverage map had no
count, so "every branch" could not be checked against a number; and every use of `<NA>` for `e` or `w`
was replaced with `nan`, since that is what a `float64` column carries.

### 20.1 Round 3 — the independent review, 2026-08-14

§19 asked for a cross-model second opinion and round 2 could not deliver one. It was then obtained from a
reviewer with **no prior context**, which read this document and ran its code. It produced **seventeen
further findings, none overlapping the first seventeen**, and verifying the fixes surfaced an eighteenth
that retracts one of round 2's conclusions, and every one it reported as measured was
re-run independently before being acted on. All held. §18d is that record.

That result is worth stating plainly rather than burying: **the first review, which ran code and found
seventeen defects, missed seventeen more of the same kind** — including two that reach the implementation
as bugs and three acceptance tests that could not have passed. §19's prediction about the expected yield
of a second reader was correct, and its estimate of the size was low.

| # | Finding | Fix | Where |
|---|---|---|---|
| 18 | The non-convergence message reports a likelihood movement of **exactly 0 on every input** — the loop rebinds `ll_old = ll_new` before the raise reads them. §12.4a asserted on that number, so the test could not fail either | `moved` bound before the rebind; the message also carries the rescale and halving counts, which after finding 19 are what separate "genuinely separated" from "ran out of travel budget" | §5.2, §5.2b, §12.4a |
| 19 | **The trust region was absolute, so convergence depended on covariate units.** A Firth estimate is equivariant under rescaling a covariate; `FIRTH_MAX_STEP × FIRTH_MAX_ITER = 1000` capped total travel absolutely. Measured: the same design fails at scale 1e-3 and converges in 9 iterations at scale 1e3. [§10] drops failed replicates, so which ones are dropped was partly an artifact of how the data was recorded | Radius made **relative** to the iterate. Measured to move nothing pinned — workbook `max|Δβ| = 0`, golden vector identical to six decimals — while fixing both the failure and the scale dependence. Recorded as a [§13] amendment with its measurement; the Newton-decrement alternative measured and declined | §5.2a, §5.6, §12.0 1b, §12.4a, DoD-10c/10d |
| 20 | `design` raised **`KeyError: 'case_id'`** on any frame carrying a declared factor and no identifier — D3 indexed the column before its own guard. An undeclared schema requirement on the one function invariant 4 is enforced through, and §12.1's own frames would have hit it | Identifiers read only when there is something to name, with an index-based fallback message. §4.5 states that a check which names patients may not impose a schema to do it | §4.5, §12.1 |
| 21 | §12.8's "the other four probabilities are **unchanged** to 1e-12" is false — a 4-record fit differs from a 5-record one; measured 0.0478. The natural repair, loosening the tolerance, tests nothing | Respecified onto what actually distinguishes "left the design" from "entered as a zero": the design has `int(in_model.sum())` rows, and a direct 4-record refit matches the blanked 5-record run to 1e-12 | §12.8 |
| 22 | §12.10's MLE companion — the half that gives §6.3 "its own explanation" — was specified on `cohort_frame()`, which is **completely separated**, so no MLE exists. Its verdict was a function of the optimiser: Newton "passed" by underflow, BFGS failed | The Firth inequality stays on `cohort_frame()` (no `data/` needed); the MLE equality is data-gated onto the workbook design, where the MLE converges and `Σw₁ − Σw₀ = 1.8e-15` | §12.10, §6.3 |
| 23 | The `design_matrix` table **never named the reference level**, because references are dropped before both `kept` and `dropped` — and §7.3 makes the log the only place the coefficients live. The one permanent record of the model was a dummy table with no baseline | The table now has one row per **level of every declared factor**, with `status = "reference (baseline)"`. §7.5 is the argument, and it restores §7.2's own "every declared thing is rendered" rule | §7.2, §7.5, §12.11 |
| 24 | SAP §7's *"describe the weighted population… state explicitly who it comprises"* was claimed discharged by a `detail` sentence saying the population is model-dependent — a caveat about the population, not a description of it, with an acceptance test that passed on one sentence | §7.6's weighted-population table: one row per [§6] covariate, unweighted and ATO-weighted means by arm. Explicitly not SMDs, which stay Stage 7's [§9] | §6.4, §7.2, §7.6, §12.11 |
| 25 | §12.4a asserted on `‖step₁‖`, the accepted step's norm and the halving count — **none recoverable from `Fit`** — while DoD-10b made them a completion gate, and the standing rule forbids the implementer inventing the fields | `Fit` gains `first_step_norm`, `rescales`, `halvings`; §3.2 justifies them on `converged_on`'s precedent and on §13's ask that Stage 10 aggregate them | §3, §3.2, §12.4a |
| 26 | §12.4a's `-inf` companion — described as "the **only** route by which the branch can be observed" — does not reach `-inf`; the unrescaled step is finite and uphill and is accepted with zero halvings | The claim struck. The branch's unreachability rests on §3.1's floor-and-rank argument, which is measured, rather than on a companion test that does not do what it says | §12.4a, §5.2 |
| 27 | The coverage map's total, added in round 2 so "every branch" could be checked against a number, was **wrong** — 53/51 against a map summing to 55/53 | Recomputed (59/57/2) and the block now says the map is authoritative if the two disagree, so the count cannot silently drift again | §12 coverage map |
| 28 | §12.13 specified a **126-row** frame — the workbook — in a section §12.0 restricts to §12.14, inside tests DoD-2 requires green with `data/` renamed | 9 rows, `cohort_frame()`'s classified frame; the workbook figure moved to §12.14's business | §12.13, §12.0 |
| 29 | Four of eight `pilots/analysis.py` line references were wrong, in the table whose purpose is separating verified from carried — and §16 instructs "lift, nearly verbatim" from lines that are not the ones cited | All corrected and re-read: `:40`, `:39`, `:366-369`, `:953`, `:82` | §3.1, §4.2, §5.4, §6.4, §12.2, §15, §16, §18 |
| 30 | §5.6's **heading** still said "the six numerics" above a block declaring seven — the same off-by-one round 2 fixed in §1 and left here | Heading corrected. This is the third recurrence of one fact stated in five places, which is what finding 33 is about | §5.6 |
| 31 | DoD-10 said F3 fires at rank *k*−2; the mechanism it describes is **one** dependency, so *k*−1. The earlier figure came from a toy frame with an accidental second collinearity | Corrected, with the toy frame's artifact named so the number is not "fixed" back | DoD-10, failure modes |
| 32 | DoD-16 gated completion on R for three questions that are **properties of published source** — `logistf`'s tolerances, its check order, and whether its score is information-scaled — all readable from a CRAN tarball. §13 called them "unverifiable from this side", conflating *cannot test* with *have not looked* | Split: **DoD-15a** requires reading `logistf`'s source and recording it, needs no R, and settles §5.3. **DoD-16** stays a blocker for what needs execution — `PSweight`'s ESS convention, and whether implemented behaviour matches documented defaults | §13, §16b.2, DoD-15a, DoD-16 |
| 34 | **Round 2's finding 9 was itself wrong, and this one is mine rather than the reviewer's.** I concluded from fourteen designs that `converged_on == "score"` was structurally unreachable, reasoning that `max\|s\| < 1e-6` forces `Δl*` three orders below `FIRTH_TOL`. That reasoning silently assumed `‖I⁻¹‖ ~ 1`. It is `Δl* ≈ ½·s²·‖I⁻¹‖`, and on an ill-conditioned design `‖I⁻¹‖` reaches 1e6-1e9, putting `Δl*` **above** `FIRTH_TOL` while the score is flat. The fourteen fixtures all shared a hidden characteristic — they were well-scaled — and the one that was not could not converge at all under the absolute trust radius, so fixing finding 19 is what made the route observable | §5.3 rewritten with the mechanism, the `‖I⁻¹‖` table, and how the wrong claim survived measurement. §12.4a now **parametrises over both routes**, asserting `"score"` on `travels_far()` and the ×1e-3 variants and `"likelihood"` on the well-scaled ones. §13's instruction to Stage 10 becomes a real signal; §14's "should the criteria be reordered" question is struck, having rested entirely on the false premise. The unreachability count in §12 drops from two branches to one | §5.3, §5.6, §12.4a, §13, §14, §16b, §16b.2, §12 coverage map |
| 33 | **The document's redundancy is what generates its contradictions.** 2878 lines for ~250 lines of code; every decision restated three to five times; findings 27, 30 and 31 are all one fact updated in some places and not others, as was round 2's `KeyError` in five | §20.2's single-source rule for the handful of facts that appear five or more times. Not a restructure: the depth does real work everywhere else, and a large edit is how two of round 2's findings were introduced | §20.2 |

**One thing the round-3 reviewer confirmed that is worth more than most of the findings.** Over 400
randomised designs it found **no false convergence** — worst gap against a `scipy` optimum 4.8e-9 — and it
verified §6.3's *explanation*, not merely its measurement. The estimator is right. What was wrong was
almost everything the document said *around* it.

**And finding 34 is the one to read twice, because nobody reported it.** It surfaced while *verifying the
fix* for finding 19: making the trust region scale-aware let `travels_far()` converge, and it converged on
the `"score"` route that round 2 had declared structurally unreachable and written into four sections, a
coverage map, a Definition-of-done item and a §14 open question. Two reviews had passed over it. The
mechanism that caught it was neither review but the habit of re-running a claim after changing the code
underneath it — which is the same habit that produced §18c, and the argument for keeping it.

It also carries a general lesson worth more than the specific correction: **a property confirmed on
fourteen fixtures that share a hidden characteristic is confirmed about the characteristic.** Every
fixture in §12.0 was well-scaled, so every fixture agreed, and the agreement read as structure. §12.4a's
×1e-3 variants exist now partly to keep that particular blind spot open to view.

### 20.1a Round 4 — the second independent review, 2026-08-14

A second reviewer with no prior context, run after the round-3 corrections landed. **17 further findings**,
again none overlapping. One is about the statistics; the rest are this document disagreeing with itself.

**The pattern is now the finding.** Round 4's own summary put it precisely: *six of its thirteen confirmed
findings are one fact fixed in one place and left in another* — the exact failure §20.2's rule was written
to close, on facts §20.2's table names as single-sourced. **The rule was written in round 3 and not
applied.** That is why §20.2 now carries a retraction policy as well as a single-source rule, and why
this round's corrections **delete** the superseded text rather than annotating it.

| # | Finding | Fix | Where |
|---|---|---|---|
| 35 | **The penalised likelihood is not concave and the loop can converge to a local maximum, silently.** Measured: 1 replicate in 2000 on this cohort, 4 of 600 synthetic separated designs, every guard clean. Contradicts §18d's "no false convergence" and §16b oracle 3's "the loop's optimum **is** the objective's optimum" | **§5.2c**, written out in full with both prevalences and all four failed detector designs. `firth` and `Fit` unchanged, so no pinned number moves; §9 assigns detection to Stage 10, which has the replicate loop and `C.SEED`; §13 files the rate and the [§13] decision it implies. §18d's row **withdrawn**, oracle 3 rescoped | §5.2c, §9, §13, §16b, §18d, §18e |
| 36 | Round 3's finding-34 retraction was applied to §5.3 and §12.4a and **left stale in 8 places**, two of them acceptance criteria (§12.0 1b, §12.4, DoD-10b) that contradict §12.4a on fixtures both name — the suite could not go green — plus §0's diagram, the failure-modes table, the coverage **count** block, §16b.2 and the roadmap | All eight corrected. §12's count block now defers to the map explicitly | §0, §12.0 1b, §12.4, §12.4a, §12, §16b.2, DoD-10b, roadmap |
| 37 | Round 3's finding-31 `k−1` correction **left as `k−2` or "3 of 5" in four places**, one of them §12.2's assertion | All four corrected to *k*−1 | §4.2, §12.2, failure modes, §18b |
| 38 | **DoD-10c requires probabilities to agree to 1e-8 across the scale variants; §12.4a says 1e-3** and calls 1e-8 "an earlier draft". DoD-10c *is* that draft — a completion gate that cannot be met | DoD-10c corrected to the measured 1e-3, with the flat-surface reason | DoD-10c |
| 39 | §12.4a's opening says the safeguards are reached by none of four fixtures "measured at zero rescales, zero halvings"; `propensity.fit` on `cohort_frame()` fires **1 rescale and 2 halvings**, contradicting §18d two sections away | Corrected; `cohort_frame()` removed from the set and its counts stated | §12.4a, §13 |
| 40 | §12.4a's trust-region companion asserts that removing the rescale takes **more** iterations and does not change the answer. Measured: **8** iterations rather than 9, and the answer moves 4.8e-06 — false in both halves | Companion struck; what the trust region buys is asserted in §5.2a's scale test instead | §12.4a |
| 41 | §12.4a says `FIRTH_MAX_ITER` exhaustion is "reached by data with no monkeypatch" while §12.0 1b, §5.6 and the coverage map all say no fixture reaches it after §5.2a | Corrected to the monkeypatched form | §12.4a |
| 42 | **§12.11 was never updated for §7.5 or §7.6.** Both of round 3's substantive fixes were specified in §7 and unasserted in §12 — and §12.11 still specified the 12-row candidate table §7.5 replaced and the one-sentence discharge §7.6 rejected. §7.6 even claims "§12.11 asserts it" | Both assertions written into §12.11: the reference rows and the 14-row count, and the weighted-population table with its unweighted column reconciled against the frame | §12.11 |
| 43 | §12.1's `PS_COVARIATES_FULL` bullet is false on the only frame §12.0 lets it run on — all four extra covariates are constant on `cohort_frame()` and are dropped | Rewritten to state both frames, with the four-column claim scoped to the data-gated workbook | §12.1 |
| 44 | §1 says `scipy` becomes "a direct import of a shipped module"; §2 says test-only; §12.7's scan forbids it in both shipped modules. An implementer following §1 fails §12.7 | §1 corrected to match §2 and §12.7 | §1 |
| 45 | The roadmap landing with this document still carried `<NA>` for `e`/`w`, and the score-route claim. §19's "five entries" ledger did not cover either | Both corrected; §19's ledger extended | roadmap, §19 |
| 46 | §20.2's closing paragraph was a **truncated sentence** with no subject, asserting that the review "is held only in a review artefact" — the opposite of its intent, in the paragraph telling the implementer §20 is optional reading | Restored | §20.2 |
| 47 | §5.5 cites §12.12 for the no-second-estimator import scan; that scan is §12.7's first bullet. §12.12 is the module-boundary scan | Corrected | §5.5 |

**One thing round 4 confirmed that is worth the whole exercise.** It checked the estimator itself rather
than the prose: the numeric gradient of `_penalised_loglik` matches the modified score to six digits. The
penalty, the score and the hat diagonal are right. Across four reviews, **no defect has ever been found in
the Firth kernel** — every one has been in what the document says around it, or in a guard that could not
fire.

### 20.1b Round 5 — the audit of what round 4 left, 2026-08-16

Round 4 produced seventeen findings and thirteen were acted on immediately. This round is the audit of the
remaining four, plus the improvements they implied. **Five entries**, and two of them add a runtime guard.

| # | Finding | Fix | Where |
|---|---|---|---|
| 48 | **A missing factor value is a third silent path to the reference level, and §4.5 claimed F4 caught it.** Measured: `pd.NA` → `NaN` → an all-zero dummy row, **byte-identical** to a real reference-level record. F4 sees no non-finite value, D3 excludes the case by its own `& notna()`, F5 sees a finite interior probability. Nothing caught it | **D4** added to `_assert_design_inputs` — a declared factor may not be missing — collected with D1-D3, never firing on the [§7] path because `fit` masks first, and firing for the direct callers §9 names. §4.5's false backstop sentence replaced by D4's argument; §12.6 gains the test and the element-for-element companion; the failure-modes table gains the row | §4.5, §12.6, §17 T2, failure modes |
| 49 | **`FIRTH_WEIGHT_FLOOR` is reachable** — 68 of 800 separated designs, smallest `p(1−p)` exactly **0.0** — so §5.6's "one constant is genuinely unreached" was true only of the tested fixtures, and §12.4a leaned on it for the `-inf` branch | §5.6 rewritten to scope the claim to the fixtures and record the general behaviour, including the two compounding consequences: `F5` fires where `p(1−p)` reaches 0 on the [§7] path but not for `model.py`'s direct callers, and a binding floor makes the information matrix and the score inconsistent — a plausible mechanism behind §5.2c, two of whose four cases had the floor binding | §5.6, §12.4a, §5.2c |
| 50 | **The seven tolerances were named `PS_*` and read by the outcome-agnostic module.** §12.12's attribute scan forbids `C.TREATMENT` and `C.PS_COVARIATES` and would have passed all seven — the letter of the boundary, not its spirit | Renamed **`FIRTH_*`** throughout: 102 occurrences across the spec and the roadmap. §12.12 gains a scan forbidding any `C.PS_`-prefixed attribute in `model.py` (`propensity.py` exempt — `PS_COVARIATES` is what it is for). A one-file change today; a four-stage change once Stages 9, 10 and 12 land against the old names, which is §0.1's own cheap-now argument applied to the names that cross the boundary | §5.6, §10, §12.12, roadmap |
| 51 | **§7.5's and §7.6's tables were specified in prose alone** — the only two in the document without code — and `_design_table(X, dropped)` cannot build the reference rows from its arguments: the reference level is by construction in neither list, and inferring factor membership from name prefixes misreads any linear covariate sharing one | Both written out. `_design_table` gains `covariates` and ranges over the **declaration** — `FACTOR_LEVELS`, `REFERENCE_LEVELS` — rather than inferring from the matrix, which is §4.2's declared-not-observed rule applied to the table. `_weighted_population_table` written out beside it, and §8's `overlap_weights` call updated to carry both | §7.5, §7.6, §8, §3 |
| 52 | §1 placed all four §16b oracles in `test_model.py`; the golden vector is §12.0.2's and lives in `test_propensity.py`, which §1's own next row says | §1 corrected: three oracles in `test_model.py`, the golden vector in `test_propensity.py`, with the reason — it pins `fit`'s output, not `firth`'s | §1 |

**What this round did not find.** No new defect in the estimator, and no new contradiction beyond the four
carried from round 4. That is the first round of five where the document did not disagree with itself in a
way nobody had already named — which is weak evidence of convergence and should not be read as more.

### 20.2 The single-source rule

Findings 27, 30, 31 and round 2's finding 2 are one failure repeated: a fact stated in five places,
updated in four. The remedy is not a restructure — the depth of this document does real work, and a
large rewrite is how findings 13 and 26 were introduced in the first place. It is a rule about the
handful of facts that are genuinely multiply-stated.

**Each of these has ONE authoritative section. Everywhere else cites it and states no number:**

| Fact | Lives in | Cited from |
|---|---|---|
| The R oracle's status as a blocker, and what it does and does not gate | **DoD-15a and DoD-16** | §13, §16b.2, §17 T6 |
| When the score route fires, and the `‖I⁻¹‖` mechanism | **§5.3** | §5.6, §12.4a, §13, §16b, §16b.2 |
| The trust region's form and the scale argument | **§5.2a** | §5.2, §5.6, §12.0 1b, §12.4a, §13, DoD-10c |
| The seven Firth constants and their count | **§5.6** | §1, §10, coverage map |
| The design's width — 12 candidates, 11 kept, 12 parameters | **§18** | §0, §4.2, §4.3, §11, §12.14 |
| Coverage totals | **§12's coverage map** | the count block beneath it, which defers to it explicitly |
| The workbook's fitted figures | **§18** | §9, §11, §12.14 |

**The rule for anyone editing this document:** if you are about to change a number that appears in more
than two sections, change it in the authoritative one and check that the others cite rather than restate
it. If a section restates it, that is the defect — fix the restatement, not just the number.

### 20.2a The retraction policy, and why it changed

Rounds 1-3 recorded corrections by **keeping the wrong statement visible** — "an earlier draft said X" —
so a reader who remembered X would learn it had changed. Round 4 established what that costs: across a
document this size, with each fact stated in several sections, corrections **accumulate beside** the
errors instead of replacing them. Findings 36, 37, 38, 41, 44 and 45 are all one correction applied in
one place and not the others, and two of them left acceptance criteria that contradict each other on
fixtures both name — a suite that cannot go green.

**From 2026-08-14 the policy is: correct in the body, record in §20.**

- The body states **only what is currently true.** A superseded sentence is deleted, not annotated.
- §20's ledger keeps the full before-and-after, so nothing is lost and a reader who remembers an earlier
  form finds it in exactly one place.
- A correction is not finished until the document has been searched for the retracted claim and every
  occurrence removed. §20.2's table names which facts are multiply-stated; those are where to look first.

Two exceptions, both deliberate. §18d and §18b keep **withdrawn rows marked as withdrawn** rather than
deleted, because §18 is a verification *record* and a measurement that was made and later disbelieved is
part of the record. And §5.2a, §5.2b, §5.3 and §7.5 each keep one sentence naming what the earlier
behaviour was, because in those four the correction **reverses** the document's apparent meaning and an
implementer working from memory would otherwise write the old form — a `KeyError` test, an absolute trust
radius, a `converged_on == "likelihood"` assertion, a 12-row design table. Everywhere else, the earlier
form is gone from the body and lives in §20.
**Where the intent of these reviews lives, for a reader who wants it in one place.** Nothing decided
during them is held only in a review artefact: §20 is the ledger, §18b-§18e are the measurements, and
every decision is written at the point of use — F6 in §5.4, the reference-drop argument in §4.2, the
completeness table in §7.4, the fixtures in §12.0, the safeguards in §12.4a, the unreachability in §5.3
and §5.6, the module-boundary string rule in §12.12, and the downstream obligations in §9. §17 and the
Definition of done were updated to match, so an implementer working the task list top to bottom reaches
every one of them without reading §20 at all. §20 is for the reader who wants to know **why** the document
says what it says; the rest of the document is sufficient to build from.

**What the review did not settle, and what it confirmed.** It confirmed §9's audit-entry ledger exactly
(load 7 / derive 4 / classify 1 / build 6), §10's claim that `test_data.py`'s two kind-parametrised tests
need no edit, §12.0.2's golden vector to every pinned digit including the `0.288470` coincidence, every
figure in §11's handover block, and the whole of §16b's library investigation. It did **not** get a
cross-model second opinion — the Codex pass failed on authentication and was not re-run — so §19's
request for a second reader is still open, and it is the one item of this review's scope that was not
delivered.
