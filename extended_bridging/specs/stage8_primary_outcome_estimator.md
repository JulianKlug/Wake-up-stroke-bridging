# Stage 8 spec — primary outcome estimator

Implements roadmap Stage 8 [§8]. Section references in brackets are to `statistical_analysis_plan.md`.
Numbers and decisions referenced as DECISION *n* are established in Stage 0 and recorded in
`../out/stage0_data_inventory.md`. Constants named in `SMALL_CAPS` are Stage 1's and live in
`config.py`; `stage1_config_and_data_contract.md` is their specification. The frame this stage
receives is specified by `stage5_cohort_construction.md` §11 and the result object it reads by
`stage6_propensity_and_weights.md` §11.

**Status.** Written 2026-08-21 against the landed Stages 1-7, and **revised the same day after its own
code was extracted and executed** (§20b), which moved five claims and one prespecified constant — §21
names each. Every number below was produced by running code; none is carried. It has since had **one
independent review rounds**, 2026-08-21/24, which together found **thirty-six** items: one logic defect in a
rendered table that no frame in the suite witnesses (§9.2, §14.0.5), one class of documentary defect
that the sole-source claim cannot survive (§14.0.4), one acceptance section that re-ran a 6100-fit
calibration on every commit (§14.7), an acceptance assertion that is false on correct code
and aimed at the workbook (§14.10, §14.12), two safeguards on the bootstrap path with no test that can
fail (§14.4), two defects the review introduced into its own fixtures and caught by executing them
(§14.0.4, §20c), and twenty-nine smaller — of which thirty-four are applied and one is recorded as
declined. §21.1 is the ledger and §20c is what was executed. It is the
**sole source for the Stage 8 implementation**: everything the implementer needs is here, and anything
not here is not to be invented.

**One number is deliberately absent, and its absence is the first decision this document takes.**
`β`, `exp(β)` and the `RD_k` **on the workbook** are not in this document, are not in any file under
`specs/`, and are not pinned in `test_outcome.py`. They are produced at implementation time into the
gitignored log. §4.3 is the argument; every design decision below was therefore made without any mRS
distribution split by arm having been examined, which is the posture Stages 0, 4, 5, 6 and 7 recorded
and the posture this stage is the last one able to hold.

**Goal.** A weighted proportional-odds model with treatment as the sole predictor returns `β` and
`exp(β)` oriented so that `> 1` favours bridging; weighted empirical cumulative risk differences
`RD_k`, `k = 0…5`, come from the weighted distributions rather than from six threshold models, so the
cumulative probabilities are ordered by construction rather than by luck; the fit raises rather than
returning whenever what it would return is not an odds ratio; and the [§11] denominator of the
estimate is named, counted and logged, because it is a third population and not the cohort's.

**Not in scope.** Every interval and every p-value (Stage 10 [§10]), the secondary and safety binary
estimators and the augmentation (Stage 9 [§8]), multiplicity, subgroups and the E-value (Stage 11
[§13]), the [§14] populations (Stages 12 and 13), and every figure and every table file (Stage 14
[§16]). Stage 8 **adds no column to the cohort frame, edits no value, removes no row, re-weights
nothing and writes no file**; it returns a result object and appends three audit entries. It does
**not** read the `Balance` (Stage 7 §11).

**One thing this spec settles that no earlier stage could.** Every stage so far could pin its own
output in its own acceptance suite, because no stage's output was the answer. Stage 7 pinned 1.343 and
0.380 into `test_balance.py` and was right to: a balance diagnostic is evidence *about* the analysis.
Stage 8's output is the analysis. §4.3 settles that the workbook's estimate is measured **at
implementation and not at specification**, that the spec's regression pins come from a synthetic
frame instead (§14.0.2), and that the data-gated acceptance tests therefore assert **properties**
rather than values. The cost is real and is stated in §15: after this document there is no pinned
regression number on the primary effect anywhere in git.

**And one finding that arrives with the stage rather than being designed into it.** The roadmap's
**Accept when** reads as though non-convergence were how an unpenalised ordinal fit fails on sparse
data. **It is not. Measured (§20): a perfectly separated 20-against-20 frame converges in 17
iterations on the score criterion, with every safeguard counter reading zero, every fitted
probability finite, and nothing anywhere out of range — and returns `β = 36.4`,
`exp(β) = 6.5e15`.** Separation does not present as failure; it presents as success. Four hundred
sparse replicates of the cohort's shape produced **zero** convergence failures, because there are
none to produce. §6 is what stands between that and a bootstrap percentile interval computed over a
`β` vector containing 36.4.

---

## 0. Where Stage 8 sits

```
  cohort.build(df, audit)         →  93 rows × 33 columns              [Stage 5 §11]
  propensity.fit(cohort, audit)   →  Propensity(e, w, in_model, …)     [Stage 6 §11]
  balance.assess(...)             →  Balance — NOT read here           [Stage 7 §11]
                       │
                       ▼
  ┌─────────────────────────────────────────────────────────────────────────────┐
  │  STAGE 8a — model.py, amended            reads config.py only               │
  │                                                                             │
  │   polr(X, y, w) -> PolrFit        [§8] weighted proportional odds   (§5)    │
  │     ├─ _assert_polr_fittable      O1…O6 → FitError                 (§5.5)  │
  │     ├─ _weighted_categories       the POSITIVE-WEIGHT level set     (§5.3)  │
  │     │                             ONE definition, TWO callers      (§5.3a)  │
  │     ├─ start: β = 0, α_k = logit of the weighted cumulative share   (§5.2)  │
  │     ├─ Newton on the WEIGHTED log-likelihood, analytic score + H            │
  │     ├─ relative trust region, halving on the log-likelihood                 │
  │     └─ raise FitError — never a second estimator     [invariant 5]          │
  └─────────────────────────────────────────────────────────────────────────────┘
                       │                                    ▲
                       │                        Stage 12 [§14a] enters HERE
                       ▼
  ┌─────────────────────────────────────────────────────────────────────────────┐
  │  STAGE 8b — outcome.py, new     reads config, data.Audit, model, propensity │
  │                                                                             │
  │   weighted_proportion(x, a, w)   the one weighted per-arm mean      (§8.1)  │
  │   cumulative_rd(y, a, w)         RD_k over MRS_THRESHOLDS           (§8)    │
  │                                                                             │
  │   primary(df, ps, audit) -> Primary                                         │
  │     ├─ _assert_primary_inputs    G1,G2,G3 THEN G4,G5 → SchemaError  (§4.4)  │
  │     ├─ in_estimate = in_model & outcome present     the [§11] third (§4.1)  │
  │     ├─ design(sub, (TREATMENT,))  → G6: the exposure SURVIVED       (§4.5)  │
  │     ├─ model.polr(X, y, w)                                          (§5)    │
  │     ├─ _assert_reportable        G7: |β| bounded → FitError         (§6)    │
  │     ├─ cumulative_rd from the weighted distributions                (§8)    │
  │     └─ three `model` audit entries                                  (§9)    │
  └─────────────────────────────────────────────────────────────────────────────┘
                       │
                       ▼
   Primary(beta, odds_ratio, alpha, cut_levels, rd, cumulative, in_estimate, fit)
                       │
                       ▼
   Stage 10 [§10] refits this in every replicate      Stage 14 [§16] reports it
```

`primary` appends three entries of the existing `model` kind. It writes no file, as no stage before it
does.

### 0.1 Why the fitter is `model.py`'s and the estimator is its own module

Stage 6 §0.1 split the numerics from the specification and reserved this boundary in writing:
*"Stages 8, 9 and 12 come into this module directly and never through `propensity.py`."*
`model.py:7-14` is that sentence in the code. Stage 8 is the first stage to arrive, and it arrives
exactly there.

**`polr` is general in its covariates from the first line, and that is not speculative generality.**
[§14a] prescribes `logit P(Y ≤ k | A, X) = α_k + βA + γᵀX` over `STANDARDISATION_COVARIATES`, which is
the same estimator with a wider design and unit weights. Writing a treatment-only fitter now and
generalising it at Stage 12 would mean two proportional-odds implementations in one repository, and
the [§14a] guard — `exp(β)` there is a **conditional** odds ratio and must never be labelled the
standardised marginal effect — is a statement about the *same* `β` this function returns. One
implementation, two callers, and Stage 12 inherits O1-O6 by calling rather than by remembering.

**`model.py` may take weights, and the sentence that seemed to forbid it does not.**
`propensity.ess`'s docstring (`propensity.py:107-109`) says a Kish sum *"has nothing to do with
fitting anything, so it would breach `model.py`'s no-exposure-no-weights rule"*. Read against
`model.py`'s own docstring, that clause is about `ess` and not about a fitter: `firth(X, y)` already
takes the **exposure** as an opaque `y` and is outcome-agnostic because it does not know what `y`
*is*. A non-negative observation-weight vector is opaque in exactly the same way. The rule being
honoured is "this module names neither the treatment nor any covariate list", and `polr(X, y, w)`
names none of the three. §12 records the one-sentence docstring amendment that says so in the file, so
the next reader does not have to re-derive it.

**`primary` is `outcome.py`'s and takes no outcome name**, for the reason `propensity.fit` takes no
covariate list (Stage 6 §6.4): [§8] names one primary quantity with one test, `PRIMARY_OUTCOME` is
computed from the registry (§12), and a parameter would make a prespecified choice look like an
option.

**And `outcome.py` is [§8]'s module, not Stage 8's.** Stage 9's weighted risk difference is
`weighted_proportion` differenced, its marginal odds ratio is the same two proportions, and its
augmentation tilts by `h = e(1 − e)` — the same weights. A `primary.py` beside a later `secondary.py`
would put the definition of "the weighted mean of `x` in arm `a`" in two files, and the two would be
one edit from disagreeing about a denominator. So §8.1's function is public and Stage 9 extends this
module.

### 0.2 What Stage 8 does not touch

`primary` takes the cohort frame and the `Propensity` and gives back neither changed: **no column is
added, no value is edited, no row is removed, and no weight is recomputed.** Stage 6 §0.2 and Stage 7
§0.2 made the same promise; here the temptation is a third one, and it is the strongest yet — the
natural way to write this is to attach `w` and `mrs_90d` to one frame and group. [§10] refits in every
one of `N_BOOT` replicates, and a frame carrying the point fit's weights, resampled into a replicate,
is a replicate weighted by the wrong score.

**Nothing here reads `Balance`.** Stage 7 §11 states it and §16 records why: no estimate in this SAP
is conditioned on a balance diagnostic, and a pipeline in which a threshold silently selects an
estimator is the mixture [§7] forbids one level up.

**No standard error, no interval, no p-value leaves this stage.** §5.6 is the argument and it is not
a division of labour but a correctness claim: the naive weighted-likelihood covariance treats the
overlap weights as frequencies, which they are not.

## 1. Deliverables

| Path | Contents |
|---|---|
| `extended_bridging/model.py` | **amended** — `PolrFit`, `polr`, and the five privates §3 lists; the module docstring gains §0.1's sentence on observation weights. No existing function changes |
| `extended_bridging/outcome.py` | the four public names `weighted_proportion`, `cumulative_rd`, `primary` and the type `Primary`; the privates §3 lists, including `_assert_primary_inputs`, `_assert_exposure_survived` and `_assert_reportable` |
| `extended_bridging/config.py` | **amended** — `MRS_LEVELS` beside `MRS_THRESHOLDS`, `PRIMARY_OUTCOME` after the outcome registry, and the `POLR_*` block after the `FIRTH_*` one (§12) |
| `extended_bridging/tests/test_model.py` | **amended** — the `polr` sections of §14, its own hand-built ordinal vectors, and **the six generators of §14.0.4**, which `test_outcome.py` imports by name exactly as it imports `cohort_frame` from `test_cohort` |
| `extended_bridging/tests/test_outcome.py` | the acceptance tests of §14, its **own** module-scoped `workbook` fixture returning `(df, ps, audit)` from **one** linear run against **one** `Audit` (§14.0.3), the fixture-derived `ordinal_cohort()` of §14.0.1 and `truncated_cohort()` of §14.0.5, and §14.0.2's golden vector |
| `extended_bridging/tests/test_config.py` | **amended** — six assertions pinning what the three new config views are computed from (§12) |
| `extended_bridging/tests/reference/polr_clm.R` | ~25 committed lines: read a CSV, ask `ordinal::clm` for a weighted proportional-odds fit, write the coefficients back. No covariate list and no analysis logic (§18b) |
| `extended_bridging/tests/test_reference_r.py` | **amended** — the Stage 8 half of the R oracle behind the gate that already exists, and the recorded fact that `Rscript` does not start on this machine (§18b) |
| `extended_bridging/implementation_roadmap.md` | **amended, and landing with this document** — Stage 8 gains its `**Spec:**` line, four **Accept when** items and one correction (§21) |

`out/logs/audit_<label>.md` is an **output, not a deliverable**, as in Stages 2-7: gitignored, written
only when a caller asks. **On this stage it is also the only place the primary estimate exists** (§4.3).

**Nothing under `specs/` may quote a case identifier, and nothing here does.** Every patient below is
a count, and every named record is a `HAND-N` or `COHORT-N` fixture record. **And nothing under
`specs/` may quote the primary effect measured on the workbook** — new at this stage, argued in §4.3,
and the reason §14.12 asserts properties where Stage 7 §12.12 asserted numbers.

## 2. Environment

`uv`, Python 3.12, flat module layout, `import config as C`. Commands run from `extended_bridging/`.

```bash
cd extended_bridging && uv sync && uv run pytest -v
```

**No new dependency, and the one that looked like a candidate is measured unusable.** `statsmodels` is
already a declared project dependency and it *has* a proportional-odds model —
`statsmodels.miscmodels.ordinal_model.OrderedModel`, version 0.14.6 in this environment. It has **no
observation-weight support of any kind**: its signature is `__init__(self, endog, exog, offset=None,
distr='probit', **kwds)` and there is no `weights`, no `freq_weights` and no `var_weights` (§20). The
roadmap anticipated this — *"most implementations have no observation-weight support"* — and it is
measured here rather than assumed, which is the same posture Stage 6 §16b took toward the absence of a
maintained Firth implementation.

`pilots/analysis.py:272-284` works around it by subclassing `OrderedModel` and overriding `loglike` to
weight `loglikeobs`, then fitting with `bfgs` and numerical derivatives. §18 records why that is not
lifted; the short version is that it inherits a `fit` whose failures the pilots catch and convert to
`nan`, and [§10] cannot count a `nan`.

So the estimator is implemented here, and it carries **five oracles instead of a comment** (§18b).

**What this costs to run, and unlike Stage 7 it matters.** Measured (§20): one 92-record
seven-category treatment-only weighted fit is **1.157 ms**, and `cumulative_rd` over the six declared
thresholds is **0.067 ms**. Stage 6's Firth fit is tens of microseconds (Stage 6 §2) and Stage 7 is
not on the bootstrap path at all (Stage 7 §2). **Stage 8 is**: [§10] refits it in every one of
`N_BOOT` replicates, so the fit alone is about **2.3 s** of the bootstrap, which is where the run's
time will be. That is not a reason to vectorise anything — it is a reason for §11 to state the cost
where Stage 10 will read it, rather than have Stage 10 discover it.

## 3. Module shape

```python
# model.py — the amendment. The existing imports already cover it.
from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass

import numpy as np
import pandas as pd

import config as C
```

```python
# outcome.py
from __future__ import annotations

from dataclasses import dataclass
from typing import Final

import numpy as np
import pandas as pd

import config as C
import model
import propensity
from data import Audit, _fmt
```

`_fmt` is imported for the reason `derive.py`, `propensity.py` and `balance.py` give: it is the
pipeline's *one* float formatter, and a second one is a second way for two runs to disagree. `Final`
is imported because this module declares **five** constants with it — §9.1's three step names and
§8.1's two arm codes.

**Two types, and the second is shorter than a reader expects for a reason.**

```python
@dataclass(frozen=True)
class PolrFit:
    """One weighted proportional-odds fit, in the `logit P(Y <= k) = alpha_k + x'beta` form (§5.1).

    `alpha` is ascending and `cut_index[j]` is the position in `categories` that `alpha[j]` cuts
    above, so a caller can say which declared level each cutpoint belongs to WITHOUT re-deriving the
    collapse (§5.3). The last four fields describe how the fit was reached rather than what it is,
    which is `Fit`'s precedent and earns its place for `Fit`'s reason (Stage 6 §3.2): [§10] refits
    N_BOOT times and a safeguard whose activation nobody can count is a safeguard nobody can
    evaluate.

    THERE IS NO STANDARD ERROR FIELD, and §5.6 is why. It is not an omission to be filled in by
    whoever needs an interval.
    """

    beta: np.ndarray            # (m,), one per design column — NO intercept, and no alpha
    alpha: np.ndarray           # (K,), ascending; K = len(categories) - 1
    categories: tuple[int, ...] # the response values FITTED, ascending — §5.3
    columns: tuple[str, ...]    # the design's column names, in beta's order
    iterations: int
    converged_on: str           # "likelihood" | "score"  — §5.2
    first_step_norm: float      # ‖step‖ at iteration 1, BEFORE any rescale   ┐
    rescales: int               # steps shortened by the trust region         ├─ §5.2
    halvings: int               # total step-halvings across all iterations   ┘
```

**The last three are `Fit`'s three, field for field** (`model.py:99-101`), and all three are asserted
rather than only the two that `_fit_detail` prints. Stage 6 §12.4a asserts on all three because a
counter nobody reads is a counter nobody maintains; §14.4 does the same here. `first_step_norm` is the
one that is *not* in the log, deliberately — a step norm is a number only a numerics reader can
interpret, and `_fit_detail` is written for the PI — so the assertion is the only thing keeping it
alive, which is why §14.4 names it explicitly instead of saying "the counters".

**`categories` is a field and not a length, and it earns its place three times.** §5.3 collapses the
response to the levels carrying positive weight, so the number of cutpoints is a property of the
sample rather than of `MRS_LEVELS`; §9.2's table renders one row per **declared** level and marks the
absent ones, which it cannot do from a count; and §11 hands Stage 10 a `β` whose comparability across
replicates depends on which categories each replicate had. A `n_categories: int` would carry the shape
and lose the identity, and the identity is what the log has to print.

```python
@dataclass(frozen=True)
class Primary:
    """The [§8] primary estimate. Point estimates only: no interval, no p-value, no standard error."""

    beta: float                              # the treatment coefficient — §5.1's parametrisation
    odds_ratio: float                        # exp(beta), oriented so > 1 favours bridging — §7
    alpha: tuple[float, ...]                 # the fitted cutpoints, ascending
    cut_levels: tuple[int, ...]              # the declared mRS level each alpha cuts above — §5.3
    rd: dict[int, float]                     # RD_k, keyed by MRS_THRESHOLDS — §8
    cumulative: dict[int, dict[int, float]]  # P(Y <= k | arm), keyed threshold then arm code — §8
    in_estimate: pd.Series                   # boolean, TOTAL — the [§11] population — §4.1
    fit: model.PolrFit
```

**`in_estimate` is boolean and total, exactly as `Propensity.in_model` is**, and for Stage 6 §9's
reason: the deliberateness of an absence lives in a mask, never in a value. Stage 8 adds no `nan`
column of its own and so has nothing to fill; what it has is a third population, and a caller that
wants the denominator reads the mask rather than counting `notna()` on something.

**There is no `favours_bridging()` method and no verdict of any kind.** Stage 7's `Balance` carries
`worst()` and `unbalanced()` because [§9] specifies a threshold and a verdict is what a threshold is
for. [§8] specifies no threshold on `β`: the null is tested at Stage 10 from the bootstrap
distribution, and a method here returning `odds_ratio > 1.0` would be a significance test with no
interval behind it, one attribute access away from a reporting layer.

Public surface, and it is four names:

```python
# model.py
def polr(X: pd.DataFrame, y: np.ndarray, w: np.ndarray | None = None) -> PolrFit: ...


# outcome.py
def weighted_proportion(x: np.ndarray, a: np.ndarray, w: np.ndarray) -> dict[int, float]: ...
def cumulative_rd(y: np.ndarray, a: np.ndarray,
                  w: np.ndarray) -> tuple[dict[int, float], dict[int, dict[int, float]]]: ...
def primary(df: pd.DataFrame, ps: propensity.Propensity, audit: Audit) -> Primary: ...
```

**Every `python` fence in this document is valid Python**, including the signature listings above —
hence the `: ...` bodies. Stage 6 §3 states the rule, Stage 7 §18d records that the fence map must be
derived by *content* and not by index, and §20b is this document's check.

`polr`'s `w` defaults to `None` and that is the one defaulted argument in this stage. It is not a
prespecified choice made to look like an option: it is the **unweighted** fit, which is Stage 12's
call and §18b's second oracle, and `None` rather than `np.ones(n)` because a caller who passes nothing
has said something different from a caller who passes ones — §14.3's companion turns on exactly that
distinction being invisible in the result and visible in the call.

Privates in `model.py` are `_weighted_categories`, `_ord_pieces`, `_ord_loglik`, `_ord_score_hess` and
`_assert_polr_fittable` — **five, and the first is new to this list**: §5.3a records why the
positive-weight rule is one function with two callers rather than one expression written twice.
Privates in `outcome.py` are `_assert_primary_inputs`, `_assert_exposure_survived`,
`_assert_reportable`, `_orientation`, `_record_outcome_completeness`, and §9's `_fit_detail`,
`_fit_table`, `_rd_detail`, `_rd_table` and `_completeness_table` — **ten, and `_orientation` was
missing from an earlier draft of this list** (§7.1 declares it and §14.5 asserts it raises), which
matters because §14.11's `inspect` assertion is written against this list and would have pinned a
surface the module does not have.

**Neither module is exempt from the Stage 1 §7 raw-name scan and neither may become exempt.** Neither
names a raw header. `test_config.py:410-412`'s scan enumerates the tree with `os.walk`, so both files are
picked up with no edit to the scan itself (Stage 7 §18).

### 3.1 The numerical facts this stage turns on

Declared once here, verified by running them on pandas 2.3.3 / numpy 1.26.4 / statsmodels 0.14.6 as
pinned by `uv.lock` (§20). Every one of them is a way this stage returns a plausible wrong number.

```
  A PERFECTLY SEPARATED FIT CONVERGES AND RETURNS.
      20 v 20, two categories, no overlap:
        beta 36.4058,  exp(beta) 6.469e+15,  17 iterations,  converged_on "score"
        halvings 0,  rescales 0,  every fitted probability finite,  nothing out of range
      Nothing is out of range. The likelihood FLATTENS as beta grows, the score goes below
      POLR_SCORE_TOL, and the loop returns normally.                              §6

  np.nan > 0.0  →  False
      So a `keep the categories whose total weight is positive` test SILENTLY DELETES a
      category whose weight vector contains one nan. Measured: a single nan weight took a
      four-category fit to three, changed beta from -1.775 to -1.205, and raised nothing. §5.5

  np.unique on a response containing nan returns nan as a value, and `w[y == nan]` is empty,
  so its total weight is 0.0 and the record is dropped by the same rule. Measured: 80
  records fitted as 79, with nothing missing anywhere in the output.               §5.5

  model.design(one_armed_frame, (TREATMENT,))  →  (0 columns, dropped ('ivt',))
      The constant-column rule is right for a covariate (Stage 6 §4.3) and CATASTROPHIC for
      the sole predictor: the estimand vanishes, `design` does not raise, and `polr` on the
      width-0 design FITS an intercept-only ordinal model and returns a `beta` of length 0.
      Measured; nothing raises on either side.                                    §4.5

  np.array([]).var(ddof=1) → nan; a weighted cumulative share of exactly 0 or 1 → logit ±inf
      Both are unreachable AS START VALUES because §5.3 collapses to the categories carrying
      positive weight, which makes every cumulative share strictly interior. The collapse is
      what makes §5.2's start values finite, and that coupling is the reason for both.  §5.2

  np.logaddexp(0.0, eta), never np.log(1 + np.exp(eta)): the second overflows at eta ~ 710.
      Measured: the largest |linear predictor| in the separated case above is 18.2, so
      POLR_ETA_CLIP at 500 is never approached — recorded, not relied on.
```

**The first fact is this stage's whole character.** Stage 6 chose Firth because *"an unpenalised MLE
fails to converge in a minority of sparse replicates"* [§7]. That sentence is true of the *propensity*
model, where the response is the exposure and a separated design genuinely stalls. It is **not** true
of the ordinal outcome fit, and the difference is not a subtlety: there the failure is loud and
[§10]'s drop-and-count handles it, here the failure is silent and drop-and-count never fires. §6 is
the guard, and it is the only reason Stage 8 has one.

## 4. The population [§11]

### 4.1 Three populations, and this is the third

```
    93   the [§3] primary cohort                                    [Stage 5 §11]
    92   in_model — covariate-complete, the ATO population          [Stage 6 §4.4]
    92   in_estimate — of those, the ones whose mrs_90d is present  [§11], §4.1
```

Stage 6 §4.4 stated the second rather than leaving it implicit — *"the ATO target population is this
complete-case set rather than the [§3] cohort"* — and the same sentence is owed here one level down.
[§11] is complete-case **per estimate** with the denominator reported for each, and the primary
estimate's denominator is neither the cohort's nor the propensity model's: it is the intersection.

**On v7 the third equals the second, and that is measured rather than assumed** (§20): `mrs_90d` is
present on all 93 cohort records, so `in_estimate` is `in_model` and the sum of weights is unchanged
to every digit. The two `N/A` cells `config.py:250-256` records in the primary outcome column are
outside the cohort. So on this workbook the outcome mask removes nobody — which is exactly the
condition under which an implementation that never built the mask would be green, and §14.2 is
therefore a constructed test for the reason Stage 7 §12.2a is.

**The mask is `in_model & present`, in that order and never `present` alone.** A record covariate-
incomplete but outcome-complete has no weight — `w` is `nan` there and Stage 6 §9 forbids filling it —
so it cannot enter a weighted mean, and a mask built from the outcome alone would put it in the
denominator of `Σw` while contributing `nan` to the numerator. Range over `in_estimate`, never over
`notna()`, and never fill: Stage 6 §9's rule, third consumer.

### 4.2 The denominator is reported per estimate, and there are three estimates

`β` is fitted on `in_estimate`. `RD_k` is computed on `in_estimate`. The two share a population **by
construction and not by coincidence**: §10 binds the mask once and both read it, because [§8]
presents `RD_k` as *"the absolute-scale presentation of the same contrast"*, and two presentations of
one contrast computed over two populations are two contrasts. §14.10 asserts they are the same object,
not merely the same number.

The per-arm weight totals are the `RD_k` denominators and they are in the log (§9.3), because a
cumulative probability whose denominator nobody can see is a proportion of nothing in particular.

### 4.3 Why the workbook's estimate is not in this document

Every stage of this pipeline records its decisions as taken before the outcome was examined by arm.
Stage 0: *"None requires looking at outcomes by arm."* DECISION 3, [§8]'s amendment, [§13]'s
amendment and Stage 4's DECISION 1a all carry the same clause. It is the strongest property this
analysis has, because it is the one a reader cannot verify from the code and can only take from the
record — and a record that says "fixed before any outcome was examined by arm" in eight places and
then quotes the treatment effect in the ninth has spent it.

**Stage 8 is where that ends, and the choice is only about where.** The estimate has to be computed
eventually. The question this section settles is whether it is computed **while the estimator is
being specified** — with every remaining decision in §5, §6, §7 and §8 then made by someone who has
seen it — or **after**, at implementation, against a specification already committed.

Three decisions in this document would have been contaminated by the first order:

- **§5.3's collapse.** Which categories to fit is a decision with a visible effect on `β`. Measured on
  the workbook (§20), all seven mRS levels are occupied, so the collapse does nothing there — but that
  is a fact about the marginal distribution, which is outcome-blind, and it was established *as* a
  marginal count.
- **§6's bound.** A bound on `|β|` chosen after seeing `β` is a bound chosen to admit it. §6 chooses
  it from the gap between 500 synthetic healthy fits and the separated constructions, and that gap is
  measured on data with no patient in it.
- **§8's thresholds.** `MRS_THRESHOLDS` is Stage 1's and prespecified, and the temptation at this
  stage is to report the `k` where the shift is largest. §8 reports all six because [§8] says so, and
  nothing in this document knows which one is largest.

So: the probe that produced §20's numbers **does not compute `β`, `exp(β)`, `RD_k`, or any mRS
distribution split by arm**, and §20 records that constraint as a property of the probe rather than as
an intention. What replaces the pinned workbook numbers:

- **§14.0.2's golden vector** is the regression pin, from a fixture-derived cohort, so an edit to this
  stage's arithmetic fails a test that contains no patient data at all.
- **§14.12's workbook assertions are properties**: the fit converged, `α` is ascending, `|β|` is inside
  the bound, `exp(β)` is finite and positive, `RD_k` is monotone in `k` within each arm, the
  orientation agrees between `β` and the sign of the `RD_k`, and the three audit entries reconcile
  with Stage 6's. Every one of them can fail; none of them quotes the answer.
- **The estimate lands in the gitignored log**, beside the coefficients Stage 6 already puts there
  (`propensity.py:56-57`), which is the file the PI reads and the file git never sees.

**What this costs, stated plainly and not minimised** (§15): there is no pinned regression number on
the primary effect anywhere in the repository. An edit that changes `β` in the fourth decimal without
changing any property above passes the whole suite. Stage 7 would have caught the analogous edit,
because it pinned 1.343. The mitigation is §14.0.2 — the same arithmetic, on a frame whose numbers may
be committed — and it is a real mitigation and not a complete one, because the golden vector's cohort
has five records and four cutpoints where the workbook's has ninety-two and six.

### 4.4 The preconditions, in two phases

Five checks in **two phases**, each phase collected and raised together, following Stage 7 §4.5 and
§4.5a — which established that the pattern every other stage uses cannot work when the later checks
read through the earlier ones' subject.

```python
def _assert_primary_inputs(df: pd.DataFrame, ps: propensity.Propensity) -> None:
    """G1-G5. Phase 1 is `can this be read`; phase 2 is `is the data estimable`. §4.4a.

    Stage 7 §4.5a is the precedent and the argument is unchanged: G4 and G5 and everything in §10
    read through the frame, the columns and the mask that G1, G2 and G3 are about, and pandas does
    not return a wrong answer for a broken mask — it raises, from inside the collection, before the
    SchemaError is ever assembled.
    """
    # PHASE 1 — can this frame, these columns and this mask be READ?
    bad: list[str] = []

    misaligned = [name for name, s in (("e", ps.e), ("w", ps.w), ("in_model", ps.in_model))
                  if not s.index.equals(df.index)]
    if misaligned:
        bad.append(
            f"G1  the Propensity is not aligned to this frame: {', '.join(misaligned)} "
            f"carr{'ies' if len(misaligned) == 1 else 'y'} a different index. "
            f"{len(ps.in_model)} mask row(s) against {len(df)} frame row(s), sharing "
            f"{len(ps.in_model.index.intersection(df.index))} index label(s). e, w and in_model are "
            "Series on the COHORT's index [Stage 6 §0.2]. §10 reads them by LABEL throughout, so the "
            "failure is not a positional read: it is that `w` arrives ordered by the Propensity's "
            "index while `y` and the arm vector come from the frame's, so each record is weighted by "
            "another record's score. Measured on the fixture cohort — a reversed index gives a common "
            "odds ratio of 58.01 against the correct 5.76, with nothing raising. A subset index is "
            "worse still: the missing labels resolve to False and leave the denominator in silence. "
            "This is Stage 7's B1 one stage on, and here the wrong number is the estimate itself.")

    if ps.in_model.dtype != bool or ps.in_model.isna().any():
        bad.append(
            f"G2  in_model is {ps.in_model.dtype} and carries "
            f"{int(ps.in_model.isna().sum())} missing value(s). It is boolean and TOTAL by "
            "construction [Stage 6 §4.4]. A missing entry does NOT raise downstream and does not "
            "propagate: `object & bool` COERCES it to False, so `in_estimate` comes back a clean "
            "boolean Series of the right length with that record silently outside the [§11] "
            "denominator. Measured for None, pd.NA and np.nan alike — 4 records estimated where 5 "
            "are in_model, nothing raised [§14.14].")

    absent = [c for c in (C.PRIMARY_OUTCOME, C.TREATMENT, "case_id") if c not in df.columns]
    if absent:
        bad.append(
            f"G3  {', '.join(absent)}: not a column of the frame. {C.PRIMARY_OUTCOME} is the [§5] "
            f"primary outcome, {C.TREATMENT} is the exposure, and case_id is what entry 1 names its "
            "removed records by; G4, §10's design and every arm mask read the first two, and "
            "`_record_outcome_completeness` dereferences the third. Stage 7's B3 had to move into "
            "phase 1 for exactly this reason [Stage 7 §4.5a, §21.1] and this check is written there "
            "from the start. case_id is in this list because without it `primary` raises a bare "
            "KeyError('case_id') from inside entry 1 — measured — which is neither of the two types "
            "§10's docstring declares, and `model.py:164-170` is the landed precedent for exactly "
            "this defect one module over.")

    if bad:
        raise C.SchemaError(
            "\n  " + "\n  ".join(bad)
            + f"\n  {len(bad)} primary-estimate assertion(s) failed over {len(df)} records. The "
            "frame, a column or the mask cannot be read, so G4 and G5 were not run: they read "
            "through all three (§4.4a).")

    # PHASE 2 — readable. These two describe the DATA, over the population §4.1 defines.
    bad = []
    in_estimate = ps.in_model & df[C.PRIMARY_OUTCOME].notna()

    present = [code for code in C.TREATMENT_LABELS
               if int((df.loc[in_estimate, C.TREATMENT] == code).sum())]
    if len(present) < 2:
        bad.append(
            "G4  " + ", ".join(
                f"{C.TREATMENT_LABELS[c]}: "
                f"{int((df.loc[in_estimate, C.TREATMENT] == c).sum())}"
                for c in C.TREATMENT_LABELS)
            + ". A treatment contrast needs both arms. With one, `design` drops the sole predictor "
            "as a constant column and `polr` fits an intercept-only ordinal model and RETURNS a beta "
            "of length zero — measured, and nothing raises on either side (§3.1, §4.5). Stage 6's F2 "
            "and Stage 7's B5 check the same property over their own populations; an arm missing "
            "HERE was lost by in_model or by the outcome mask.")

    nan_w = int((in_estimate & ps.w.isna()).sum())
    if nan_w:
        bad.append(
            f"G5  {nan_w} record(s) are in the estimation population and carry no weight. `w` is "
            "nan off in_model by construction [Stage 6 §9], so this is unreachable while "
            "in_estimate is built from in_model — and it is checked because a nan weight does not "
            "propagate here, it DELETES: `nan > 0.0` is False, so §5.3's positive-weight test drops "
            "the whole outcome category that record sits in. Measured: one nan weight took a "
            "four-category fit to three and moved beta by 0.57, with nothing raising (§3.1). Never "
            "fill it [Stage 6 §9]; the frame is wrong, not the value.")

    if bad:
        raise C.SchemaError(
            "\n  " + "\n  ".join(bad)
            + f"\n  {len(bad)} primary-estimate assertion(s) failed over {len(df)} records, "
            f"{int(in_estimate.sum())} of them estimable.")
```

- **G1 checks the index and not the length, on all three Series**, and names which is misaligned.
  Stage 7 §21.4 found that an `e`-only misalignment reached `_centre` and raised pandas'
  `IndexingError` after an entry had been recorded; the lesson is inherited rather than re-learned.
- **G5 is unreachable while §10 builds `in_estimate` from `in_model`, and it is not defensive.** It is
  the only check that names the deletion mechanism, and the mechanism is the one thing about this
  stage that a reader will not guess: a `nan` weight in this estimator does not produce a `nan`
  estimate, it produces a *different, finite, plausible* estimate over a coarser outcome scale.
  Recorded as unreachable-on-the-declared-path exactly as Stage 6's D1-D4 and Stage 7's B4 are.
- **G4's message names the consequence and not the rule.** "A propensity model needs both arms"
  (Stage 6's F2) is a statement about identifiability. Here the consequence is concrete and measured,
  and the message carries it because the failure it describes is silent on both sides of the call.

### 4.4a Why the phases, in this stage's own terms

Stage 7 §4.5a measured three pandas failures that discard an assembled message. All three apply here
unchanged, because G4 reads `df.loc[in_estimate, C.TREATMENT]` and `in_estimate` is built from
`ps.in_model` and `df[C.PRIMARY_OUTCOME]` — a misaligned mask, a three-valued mask and an absent
column each raise from inside the collection. This stage adds a fourth, and it is the reason G3 names
the **outcome** column and not only the exposure:

```
  the PRIMARY_OUTCOME column dropped   →  G3 appended, then
      ps.in_model & df[C.PRIMARY_OUTCOME].notna()   →  KeyError('mrs_90d')
```

Phase 2's first line builds the mask, so an absent outcome column raises there — before G4 and before
any message is returned. Stage 7's E5/§21.1 pair is the precedent: naming a column in an
absent-column check is necessary and is **not** sufficient unless that check raises before the column
is read. Both columns are therefore in G3, and G3 is in phase 1.

### 4.5 The exposure must survive the design, and `design` will not say so

`primary` builds its design through `model.design(sub, (C.TREATMENT,))` rather than by assembling a
column itself, and the reason is `model.py:47-50`: `design` is the one function every model's
covariates pass through, D2 is invariant 4's assertion against the post-time-zero denylist, and
*"Stages 8, 9 and 12 inherit it by calling `design` rather than by remembering"*. Stage 8 is the
first to inherit it, and it does so by calling.

**What it must not inherit is the constant-column rule.** Measured (§20):

```
    design(cohort,      (TREATMENT,))  →  1 column ('ivt',),  dropped ()
    design(one_armed,   (TREATMENT,))  →  0 columns,          dropped ('ivt',)
```

The second is correct behaviour for a covariate — Stage 6 §4.3 drops constant columns so a resampled
cohort that empties a factor level does not make the design singular — and it is a catastrophe for the
**sole predictor**: the design comes back with no treatment column, `polr` fits an intercept-only
ordinal model over the cutpoints alone, converges, and returns a `beta` of length zero. Neither
function raises. Whether the caller then gets an `IndexError` on `beta[0]` or a silent `1.0` from
`exp(beta.sum())` depends on how the reporting layer happens to be written, and both are worse than a
raise.

G4 catches the one-armed frame that causes it. `_assert_exposure_survived` catches everything else —
and there is an "everything else", because `dropped` is computed from the matrix and not from the arm
counts, so any future change to `design`'s rule reaches here first:

```python
def _assert_exposure_survived(X: pd.DataFrame, dropped: tuple[str, ...]) -> None:
    """The design handed to `polr` carries the exposure, and exactly it. §4.5.

    Not folded into G4. G4 is about the DATA — both arms are present — and this is about the
    MATRIX, and the two can come apart: `design`'s constant-column rule is keyed on
    `nunique(dropna=False)` over the built column, so a future rule change, a covariate list
    that grew, or an exposure column of one distinct value for a reason G4 does not model all
    arrive here. A design of width 0 is a model with no treatment term whose fit CONVERGES
    (§3.1, measured), so this is the last point at which the estimand is still checkable.
    """
    if tuple(X.columns) != (C.TREATMENT,) or dropped:
        raise model.FitError(
            f"G6  the [§8] design is {tuple(X.columns)} with {dropped or 'nothing'} dropped, and "
            f"[§8] prescribes treatment as the SOLE predictor: exactly ({C.TREATMENT!r},) and "
            "nothing dropped. A width-0 design fits an intercept-only ordinal model, converges, and "
            "returns a beta of length zero — the estimand vanishes and neither `design` nor `polr` "
            "raises (§4.5, measured). `design` drops constant columns by declaration "
            "[Stage 6 §4.3], which is right for a covariate and wrong for the exposure.")
```

**It raises `model.FitError` and not `SchemaError`**, and that placement is deliberate: [§10] catches
`FitError` to drop and count a replicate and *"may catch nothing else: a SchemaError here is a bug in
the resampler, not a sparse replicate"* (`propensity.py:468-469`). A replicate whose stratified
resample happened to draw one arm at every centre is a sparse replicate, not a bug — so it must be
droppable, which means it must raise the type Stage 10 catches.

## 5. The weighted proportional-odds model [§8]

### 5.1 The rule, and the sign

```
    logit P(Y <= k | A, X)  =  alpha_k  +  beta*A  +  gamma'X ,        k = 0 … K-1

    alpha_1 < alpha_2 < … < alpha_K        the cutpoints, ascending
    beta                                   the treatment coefficient
    exp(beta)                              the common odds ratio, > 1 FAVOURS BRIDGING

    The log-likelihood is WEIGHTED per observation:

        l(theta)  =  SUM_i  w_i * log( P(Y = y_i | A_i, X_i) )

    with P(Y = j) = P(Y <= j) - P(Y <= j-1),  P(Y <= -1) = 0,  P(Y <= K) = 1.
```

**The `+ beta*A` is the whole of §7 and half of this stage's failure surface.** [§14a] writes the model
in exactly this form — `logit P(Y ≤ k | A, X) = α_k + βA + γᵀX` — so the parametrisation is the SAP's
own and not a choice made here. Every reference implementation writes `α_k − x'β` instead. §7 is what
that costs and how it is pinned.

**The numerator is weighted and nothing else is.** There is no weighted variance, no weighted
information adjustment, and no design-effect correction anywhere in this estimator, because [§10] owns
inference and §5.6 is why that is a correctness statement rather than a scope statement.

### 5.2 The estimator, written out

Written out in full, as Stage 5 §5.1, Stage 6 §5.2 and Stage 7 §5.2 are: it is the only place in this
pipeline where the primary quantity is computed.

```python
def _weighted_categories(y: np.ndarray, w: np.ndarray) -> tuple[int, ...]:
    """The response levels carrying POSITIVE total weight, ascending. §5.3.

    One function and two callers — `_assert_polr_fittable`'s O5 counts what this returns and `polr`
    collapses to it — because the expression is not the obvious one and §5.5 turns on a subtlety
    inside it: `np.nan > 0.0` is False, so a single nan weight makes its level's total weight nan and
    DELETES the level rather than propagating. Written twice, the two copies are one edit from
    disagreeing about which levels a fit has, and the disagreement would be silent in both directions
    — O5 would count a different number of categories from the number `polr` then fits, and both
    numbers would be finite and plausible (§5.5, measured).

    `int(v)` is safe because O4 has already established that every response value is integral, and O4
    precedes O5 for that reason as well as for §5.5's. This function does NOT re-run O3: it is called
    from inside O5, so a nan weight has already raised by the time `polr` calls it a second time.
    """
    return tuple(int(v) for v in np.unique(y) if w[y == v].sum() > 0.0)
```

```python
def _ord_pieces(Xb: np.ndarray, alpha: np.ndarray, upper: np.ndarray, lower: np.ndarray,
                has_u: np.ndarray, has_l: np.ndarray) -> tuple[np.ndarray, ...]:
    """The two cumulative probabilities bracketing each observation's category, and two derivatives.

    `upper`/`lower` are cutpoint indices and `has_u`/`has_l` say whether each exists: an observation
    in the LOWEST category has no lower cutpoint (P(Y <= -1) = 0) and one in the HIGHEST has no upper
    one (P(Y <= K) = 1). Those two boundary conventions are what make the first and last categories
    contribute a single-sided derivative rather than a special case.

    The clip on the linear predictor is the one `firth` uses and for the same reason: it keeps `exp`
    in range. Measured (§20): the largest |linear predictor| this stage has produced is 18.2, in the
    separated construction of §6, so POLR_ETA_CLIP is never approached on any real input — recorded
    rather than relied on.
    """
    zu = np.where(has_u, alpha[np.clip(upper, 0, len(alpha) - 1)] + Xb, 0.0)
    zl = np.where(has_l, alpha[np.clip(lower, 0, len(alpha) - 1)] + Xb, 0.0)
    gu = np.where(has_u, 1.0 / (1.0 + np.exp(-np.clip(zu, -C.POLR_ETA_CLIP, C.POLR_ETA_CLIP))), 1.0)
    gl = np.where(has_l, 1.0 / (1.0 + np.exp(-np.clip(zl, -C.POLR_ETA_CLIP, C.POLR_ETA_CLIP))), 0.0)
    du = np.where(has_u, gu * (1.0 - gu), 0.0)              # d gamma_u / d z
    dl = np.where(has_l, gl * (1.0 - gl), 0.0)
    ddu = np.where(has_u, du * (1.0 - 2.0 * gu), 0.0)       # d2 gamma_u / d z2
    ddl = np.where(has_l, dl * (1.0 - 2.0 * gl), 0.0)
    return gu, gl, du, dl, ddu, ddl
```

```python
def _ord_loglik(Xn: np.ndarray, y_idx: np.ndarray, w: np.ndarray, alpha: np.ndarray,
                beta: np.ndarray, K: int) -> float:
    """The weighted log-likelihood, or -inf where any category probability is non-positive.

    -inf rather than a raise, so that step-halving can REJECT a step and try a shorter one; a raise
    there would turn a recoverable step into a failed fit. This is `_penalised_loglik`'s posture
    (Stage 6 §5.2) and it is the branch §5.4 is about: a crossed pair of cutpoints makes some
    observation's category probability non-positive, so the halving rejects the crossing rather than
    an assertion catching it.

    THE BRANCH IS MEASURED NOT TO FIRE on any input this stage has been given (§5.4, §20): 3000
    attempts at reaching a crossing with pure Newton from a crowded start never lost the ordering.
    It is kept because it is the CONTRACT between §5.3's collapse and the halving loop — remove the
    collapse, or start from an unordered alpha, and this branch becomes the thing that stops a
    crash. Do not delete it as dead code without deleting the reason it is dead.
    """
    upper, lower = y_idx, y_idx - 1
    has_u, has_l = y_idx <= K - 1, y_idx >= 1
    gu, gl = _ord_pieces(Xn @ beta, alpha, upper, lower, has_u, has_l)[:2]
    p = gu - gl
    if np.any(p <= 0.0) or not np.all(np.isfinite(p)):
        return -np.inf
    return float(np.sum(w * np.log(p)))
```

```python
def _ord_score_hess(Xn: np.ndarray, y_idx: np.ndarray, w: np.ndarray, alpha: np.ndarray,
                    beta: np.ndarray, K: int) -> tuple[np.ndarray, np.ndarray]:
    """The analytic weighted score and Hessian over (alpha, beta), in that order.

    Analytic and not differenced, for the reason `pilots/analysis.py:272-284` is not lifted (§18):
    that class inherits `OrderedModel`'s numerical derivatives, and [§10] refits N_BOOT times, so a
    finite-difference gradient is 2*(K+m) extra likelihood evaluations per Newton step forever.
    Verified against central differences to 2.8e-08 on both, and the Hessian is exactly symmetric
    and negative definite (§20).

    The derivation, once, because a wrong sign here is a plausible number: with u = P(Y <= j),
    v = P(Y <= j-1) and p = u - v, the per-observation log-likelihood is log p and

        d/dzu  = u'/p = du/p          d2/dzu2      = ddu/p - (du/p)^2
        d/dzl  = -v'/p = -dl/p        d2/dzl2      = -ddl/p - (dl/p)^2
                                      d2/dzu dzl   = du*dl/p^2

    and zu = alpha_ju + x'beta, zl = alpha_jl + x'beta, so d z_k / d alpha_m = 1{k = m} and
    d z_k / d beta = x. The beta block therefore collects Huu + 2*Hul + Hll, which is the sum of
    both single-sided terms plus twice the cross term, and the alpha-beta blocks collect one
    single-sided term plus the cross term each.
    """
    n, m = Xn.shape
    upper, lower = y_idx, y_idx - 1
    has_u, has_l = y_idx <= K - 1, y_idx >= 1
    gu, gl, du, dl, ddu, ddl = _ord_pieces(Xn @ beta, alpha, upper, lower, has_u, has_l)
    p = gu - gl

    A, B = du / p, -dl / p
    Huu = ddu / p - (du / p) ** 2
    Hll = -ddl / p - (dl / p) ** 2
    Hul = du * dl / p ** 2

    g = np.zeros(K + m)
    H = np.zeros((K + m, K + m))
    np.add.at(g, upper[has_u], (w * A)[has_u])
    np.add.at(g, lower[has_l], (w * B)[has_l])
    g[K:] = Xn.T @ (w * (A + B))

    np.add.at(H, (upper[has_u], upper[has_u]), (w * Huu)[has_u])
    np.add.at(H, (lower[has_l], lower[has_l]), (w * Hll)[has_l])
    both = has_u & has_l
    np.add.at(H, (upper[both], lower[both]), (w * Hul)[both])
    np.add.at(H, (lower[both], upper[both]), (w * Hul)[both])

    cu = np.where(both, Huu + Hul, np.where(has_u, Huu, 0.0))
    cl = np.where(both, Hll + Hul, np.where(has_l, Hll, 0.0))
    for r in range(m):
        column = np.zeros(K)
        np.add.at(column, upper[has_u], (w * cu * Xn[:, r])[has_u])
        np.add.at(column, lower[has_l], (w * cl * Xn[:, r])[has_l])
        H[:K, K + r] += column
        H[K + r, :K] += column
    H[K:, K:] = (Xn * (w * (Huu + 2.0 * Hul + Hll))[:, None]).T @ Xn
    return g, H
```

```python
def polr(X: pd.DataFrame, y: np.ndarray, w: np.ndarray | None = None) -> PolrFit:
    """Weighted proportional-odds regression [§8]. Raises FitError; never returns a fallback.

    `X` carries NO intercept and must not: the K cutpoints are the intercepts, and adding a column of
    ones makes the design exactly singular against their sum. That is the opposite of `firth`, which
    prepends its own — so the two functions differ in the one place a reader will assume they agree,
    and O6's rank check is what turns the mistake into a raise rather than into a pseudo-inverse.

    `w` is a non-negative observation-weight vector or None. None is the UNWEIGHTED fit and is not a
    synonym for ones (§3): §14.3's companion turns on the two being indistinguishable in the result.

    Four numerical details, each load-bearing:

    * **The response is collapsed to the categories carrying positive weight** (§5.3), before
      anything else numeric happens, and `categories` reports what was fitted.
    * **The start values are beta = 0 and alpha_k = logit of the weighted cumulative share.** At
      beta = 0 that is the exact maximiser of the intercept-only model, so iteration 1 starts at the
      solution of a nested model rather than at an arbitrary point — measured 4 to 5 iterations on
      every frame in §20. And they are FINITE BY CONSTRUCTION because of the collapse: every kept
      category carries positive weight, so every cumulative share is strictly inside (0, 1). The
      collapse and the start values are one decision seen twice.
    * **The step is rescaled, not rejected**, when it exceeds the trust radius, and the radius is
      RELATIVE to the iterate — `firth`'s correction (Stage 6 §5.2a), inherited for its reason: a
      proportional-odds estimate is equivariant under rescaling a covariate and an absolute bound is
      not, so with one, WHICH replicates [§10] drops would depend on how the data was recorded.
    * **`np.linalg.solve`, not `inv`, and it is deliberately NOT wrapped.** Stage 6 §5.2's argument
      applies unchanged: a `pinv` fallback silently picks the minimum-norm solution among infinitely
      many, so the fit returns coefficients for a design that identifies none. O6 establishes full
      rank before the loop begins.

    Unlike `firth` there is no penalty and no `½log|I|` term, so **the objective is concave** — the
    weighted sum of concave per-observation log-likelihoods, verified negative definite at the optimum
    and away from it (§20). Stage 6 §5.2c's one unguarded failure mode, "this loop finds *a*
    stationary point where [§7] prescribes *the* maximiser", therefore does not recur here: for this
    estimator they are the same point. That is the one respect in which Stage 8's numerics are safer
    than Stage 6's, and §6 is the respect in which they are far more dangerous.
    """
    Xn = np.asarray(X.to_numpy(dtype=float))
    y = np.asarray(y, dtype=float)
    w = np.ones(len(y)) if w is None else np.asarray(w, dtype=float)
    _assert_polr_fittable(Xn, y, w, tuple(X.columns))                      # O1-O6, §5.5

    categories = _weighted_categories(y, w)                                # §5.3, §5.3a
    keep = np.isin(y, np.asarray(categories, dtype=float))
    Xn, y, w = Xn[keep], y[keep], w[keep]
    K, m = len(categories) - 1, Xn.shape[1]
    y_idx = np.searchsorted(np.asarray(categories, dtype=float), y)

    share = np.array([w[y_idx <= k].sum() / w.sum() for k in range(K)])    # strictly interior
    par = np.concatenate([np.log(share / (1.0 - share)), np.zeros(m)])
    ll_old = _ord_loglik(Xn, y_idx, w, par[:K], par[K:], K)
    first_step_norm, rescales, halvings = 0.0, 0, 0

    for iteration in range(1, C.POLR_MAX_ITER + 1):
        g, H = _ord_score_hess(Xn, y_idx, w, par[:K], par[K:], K)
        step = -np.linalg.solve(H, g)                       # LinAlgError is NOT caught — §5.5
        norm = float(np.linalg.norm(step))
        if iteration == 1:
            first_step_norm = norm
        size = norm / max(1.0, float(np.linalg.norm(par)))
        if size > C.POLR_MAX_STEP:
            step = step * (C.POLR_MAX_STEP / size)
            rescales += 1

        for _ in range(C.POLR_MAX_HALVINGS):
            ll_new = _ord_loglik(Xn, y_idx, w, (par + step)[:K], (par + step)[K:], K)
            if ll_new >= ll_old:
                break
            step = step / 2.0
            halvings += 1
        else:
            raise model_fit_error(
                f"polr: step-halving exhausted {C.POLR_MAX_HALVINGS} halvings at iteration "
                f"{iteration} without increasing the weighted log-likelihood. [§8] prescribes one "
                "estimator; there is no second one to try [invariant 5].")

        par = par + step
        moved = abs(ll_new - ll_old)
        if moved < C.POLR_TOL:
            return PolrFit(par[K:], par[:K], categories, tuple(X.columns), iteration,
                           "likelihood", first_step_norm, rescales, halvings)
        if float(np.max(np.abs(g))) < C.POLR_SCORE_TOL:
            return PolrFit(par[K:], par[:K], categories, tuple(X.columns), iteration,
                           "score", first_step_norm, rescales, halvings)
        ll_old = ll_new

    raise model_fit_error(
        f"polr: no convergence in {C.POLR_MAX_ITER} iterations. The weighted log-likelihood moved "
        f"{moved:.3g} against a tolerance of {C.POLR_TOL:g}, and the largest score component was "
        f"{float(np.max(np.abs(g))):.3g} against {C.POLR_SCORE_TOL:g}. {rescales} step(s) were "
        f"shortened by the trust region and {halvings} halving(s) taken. [§8] prescribes one "
        "estimator, and [§10] drops and counts the replicate rather than substituting another.")
```

`model_fit_error` above is `FitError` — the name is spelled that way in this document's fences **only**
so that the assembled module of §20b can bind it without redeclaring `model.py`'s class; in
`model.py` the two raises read `raise FitError(...)`, the class that is already there
(`model.py:73`). No second exception type is declared by this stage, in either module.

**Two convergence routes, either one sufficient, and `PolrFit` records which fired** — `firth`'s rule
(Stage 6 §5.3), and the step norm is deliberately not a third, for Stage 6's reason. The order the two
run in is part of the estimator, because [§10] refits in every replicate: reordering them changes the
sampling distribution, which is a [§13] amendment and not a repair.

**And the score route is the one that fires on a separated fit.** Measured (§20): the perfectly
separated construction returns `converged_on == "score"` at iteration 17. Stage 6 §5.3 established
that a fit coming back "score" is telling you its design is close to singular; here it is telling you
something stronger and worse, and §6 is what reads it.

**`g` is bound at the top of the iteration and both tests run after the step, so a `PolrFit` reporting
"score" reports a flat score one iterate behind its own coefficients.** `firth` has exactly this
property and `model.py:492-495` states it rather than leaving it to be found — *"harmless — when the
score is below the tolerance the step is negligible — but stated because §12 asserts `converged_on` as
a property of the return"*. It is inherited here with the sentence, not without it: §14.4 asserts
`converged_on` on the returned object, and an implementer who recomputed the score after the step to
"fix" this would change which iteration the loop exits on, which changes the sampling distribution
[§10] takes percentiles of, which is a [§13] amendment and not a tidy-up.

**Both tolerances are ABSOLUTE on quantities that scale with `Σw`, and this stage is the first with two
callers at different weight scales.** `POLR_TOL` bounds `|Δℓ|` and `POLR_SCORE_TOL` bounds
`max |score|`; scaling every weight by `c` scales both by `c`, so neither is scale-free — which is the
same objection §5.2 raises against an absolute *trust radius* and resolves there by making the radius
relative. It is not resolved here, and that is a decision rather than an oversight. Three reasons, in
order of weight:

- **Changing them is a [§13] amendment.** `FIRTH_TOL` and `FIRTH_SCORE_TOL` are absolute on the same
  grounds and are landed; a `POLR_*` pair normalised by `Σw` would be the only tolerance in the
  pipeline with a different shape, and [§10] refits in every replicate so the shape decides which
  replicates converge on which route.
- **Both declared callers are measured, and both are comfortable.** §20 records `Σw` at
  **27.736623** on the workbook's estimation population under overlap weights, against `Σw = n` under
  Stage 12's unit weights over a larger population — so the two callers differ by roughly an order of
  magnitude and no more. `POLR_TOL = 1e-8` at `Σw = 27.74` is a tolerance of `3.6e-10` per unit of
  weight; at Stage 12's scale it is tighter still. The risk an absolute tolerance carries is being too
  *loose* where the objective is small, and one part in `2.8e9` of the total weight is not loose.
- **The loop cannot be starved by it.** Measured 4 to 5 iterations against `POLR_MAX_ITER = 200`
  (§20), so a tolerance that is effectively tighter at one scale than another costs iterations and not
  correctness. §14.4 asserts `iterations` on every frame in the suite, which is what would show a
  scale at which this stopped being true.

What is **not** known is `Σw` inside a [§10] replicate, because that needs the replicate loop. §11
hands it over.

### 5.3 The response is collapsed to the categories carrying positive weight

```python
categories = _weighted_categories(y, w)   # §5.2's line; the rule itself is §5.3a's function
```

**Not `MRS_LEVELS`, and not `np.unique(y)` either.** Three candidate rules, and only the third works:

- **The declared level set.** `MRS_LEVELS` is seven values, so K would be 6 always. A category nobody
  is in has `P(Y = j) = 0` for every observation, its two bracketing cutpoints are unidentified, and
  the fit runs one of them to a boundary. The number of parameters would be a property of `config.py`
  and the identifiability a property of the sample, which is the worst pairing available.
- **The observed level set**, `np.unique(y)`. Correct whenever every observed record has positive
  weight, and this estimator is general: a caller may pass a weight of zero. A category all of whose
  records carry zero weight contributes nothing to the likelihood, so its cutpoints are unidentified
  exactly as in the first case — and worse, §5.4's crossing argument fails, because a crossing between
  two cutpoints only makes a category probability non-positive **for observations in that category**,
  and if they all have zero weight the `-inf` never fires.
- **The positively-weighted level set**, which is what is written. It is the observed set on any input
  where `w > 0`, so on the [§8] path it changes nothing — overlap weights are strictly interior by
  Stage 6's F5, so `w = 1 − e` or `e` is strictly positive on every in-model record. And it is the set
  that makes both the start values (§5.2) and the crossing argument (§5.4) true rather than nearly
  true.

Measured (§20): a frame with a middle category emptied fits with **3** cutpoints against the full
frame's 4, `α` ascending in both, start values finite in both; and a frame with that category present
but carrying zero total weight gives a fit **identical to 1e-10** to the fit on the frame with those
rows removed.

**So `K` is a property of the sample, and `PolrFit.categories` is how a caller learns which.** Three
consequences, and the third is the one Stage 10 owes:

- On the workbook all seven mRS levels are occupied, so `K` is 6 and the collapse does nothing (§20,
  marginal count only). §14.1 asserts that with the *reason*, so the branch is known to be inert here
  rather than assumed to be.
- §9.2's audit table renders one row per **declared** level and marks the absent ones, because a
  category that vanishes from the table is a category nobody knows was considered — Stage 6 §7.5's
  rule, third application.
- **[§10] takes percentiles of `β` across replicates whose `K` may differ.** A replicate that misses
  mRS 0 entirely fits five cutpoints, and its `β` is a common odds ratio over a coarser scale. Under
  proportional odds those are the same parameter, which is exactly the assumption [§8] makes and does
  not test; §15 files it and §11 hands it to Stage 10 rather than leaving it to be discovered there.
  **And §6.3's band was measured only at `K = 6`**, which is the second thing §11 owes Stage 10 and is
  filed in §15.

### 5.3a One rule, one function, two callers

The positive-weight level set is computed in **two** places — O5 counts it (§5.5) and `polr` collapses
to it (§5.2) — and it is therefore **one function**, `_weighted_categories`, and not one expression
written twice. An earlier draft of this document wrote it out at both sites, with different casts
(`tuple(v for v in …)` at O5 and `tuple(int(v) for v in …)` in `polr`), which is the shape the rest of
this pipeline is careful to avoid: Stage 6 §3 made `ess` public rather than let two stages define the
effective sample size, Stage 7 §0.1 did the same for `smd`, and §8.1 does it here for the weighted
per-arm proportion.

**The duplication mattered more than the usual DRY argument**, and that is why it is a numbered
subsection rather than a note. §5.5's whole reason for ordering O4 before O5 is a subtlety *inside this
expression* — `np.nan > 0.0` is False, so a nan weight deletes a level rather than propagating. Two
copies of an expression whose correctness depends on a fact about `nan` are two chances to fix one and
not the other, and the failure would be silent in both directions: O5 would count a different number
of categories from the number `polr` then fits, and both numbers would be finite, plausible, and
reported. §14.8's O3-before-O5 pair tests the fact; this section is what stops there being two places
for the fact to be true of.

### 5.4 Why a crossed cutpoint is unreachable rather than caught

The cutpoints must satisfy `α_1 < … < α_K` or the model assigns a negative probability to some
category. Nothing in `polr` asserts it, and nothing needs to:

```
    the start values are ORDERED          (cumulative shares are non-decreasing, §5.2)
    a crossed pair gives p <= 0           for the observations between them
    p <= 0  →  _ord_loglik returns -inf   (§5.2)
    -inf < ll_old  →  the step is HALVED  and a shorter one is tried
```

so no accepted iterate is unordered, and the returned `α` is ascending because every one before it was.
This is a stronger guarantee than an assertion would give: an assertion detects a crossing after the
step has been taken, and this makes the step unavailable.

**Measured, and the measurement is a negative result that is recorded rather than smoothed** (§20):
3000 attempts to reach a crossing with **pure Newton from a deliberately crowded start** — the
cutpoints initialised to a span of 2e-3 — never lost the ordering on any of them. So the `-inf` branch
does not fire on any input this stage has been given, and §14.4 asserts the ordering as a
post-condition of the return rather than asserting that the branch fires. That is Stage 6
`_penalised_loglik`'s posture exactly: *"THE BRANCH DOES NOT FIRE, and that is measured rather than
hoped… it is kept because it is the contract"*. An earlier draft of this section claimed the crossing
was reachable by removing the halving; it is not, and the claim is struck rather than weakened.

### 5.5 O1-O6, and the order of two of them is load-bearing

```python
def _assert_polr_fittable(Xn: np.ndarray, y: np.ndarray, w: np.ndarray,
                          columns: tuple[str, ...]) -> None:
    """O1-O6 — the six things that make `polr` return, or blame, the wrong answer.

    Ordered cheap-to-expensive and raising at the FIRST failure rather than collecting, as
    `_assert_fittable` does (Stage 6 §5.4) and for its reason: each one makes the next meaningless.
    A rank computed over a design containing nan is not a rank.

    **O4 MUST precede O5, and that ordering is the finding of §3.1.** O5 counts the categories
    carrying positive weight. `np.nan > 0.0` is False, so a single nan weight makes its category's
    total weight nan, the test drops the category, and O5 counts one fewer — measured: a
    four-category fit became a three-category fit, beta moved from -1.775 to -1.205, and nothing
    raised. With O4 first, the nan is a raise; with O5 first, the nan is a different model.

    The messages name neither the exposure nor any outcome, which is a constraint and not a
    stylistic choice: this module is outcome-agnostic (Stage 6 §0.1) and §14.13 scans its raise
    strings for both. The concrete examples live in this spec instead.
    """
    if not np.all(np.isfinite(Xn)):
        bad = [columns[j] for j in np.unique(np.argwhere(~np.isfinite(Xn))[:, 1])]
        raise model_fit_error(
            f"O1  the design carries a non-finite value in: {', '.join(bad)}. An Int64 column with "
            "pd.NA converts to nan SILENTLY, and every candidate step is then nan, so no step is "
            "accepted and the fit blames step-halving for a failure of the data. Use "
            "model.complete_cases().")

    if not np.all(np.isfinite(y)):
        raise model_fit_error(
            f"O2  the response carries {int((~np.isfinite(y)).sum())} non-finite value(s) of "
            f"{len(y)}. It does NOT propagate here — np.unique returns nan as a level, its total "
            "weight is 0, and the record is DROPPED by the positive-weight rule with nothing "
            "missing anywhere in the result. Measured: 80 records fitted as 79. The response is "
            "complete-case per estimate and the mask is the caller's.")

    if not np.all(np.isfinite(w)) or np.any(w < 0.0):
        raise model_fit_error(
            f"O3  the weights carry {int((~np.isfinite(w)).sum())} non-finite and "
            f"{int(np.sum(w < 0.0))} negative value(s). A non-finite weight DELETES an outcome "
            "category rather than propagating (see O5's ordering); a negative one makes the "
            "objective non-concave, so the Hessian may be indefinite and a Newton step may ascend "
            "away from the maximiser while every counter reads clean.")

    off = y[np.not_equal(np.mod(y, 1.0), 0.0)]
    if off.size:
        raise model_fit_error(
            f"O4  the response takes {off.size} non-integer value(s): "
            f"{', '.join(f'{v:g}' for v in np.unique(off)[:10])}. This is an ORDINAL fit and the "
            "response is a category code: one cutpoint is estimated per distinct value, so a "
            "continuous response silently makes the parameter count the number of distinct values "
            "and every category a singleton. It fits. A continuous response belongs to no model in "
            "this pipeline.")

    categories = _weighted_categories(y, w)                                # §5.3a, one definition
    if len(categories) < 2:
        raise model_fit_error(
            f"O5  {len(categories)} response category/ies carry positive weight: "
            f"{tuple(float(v) for v in categories)}. A proportional-odds fit needs at least two: "
            "with one, there is no cutpoint to estimate and the likelihood is constant in every "
            "parameter. This is the branch the committed fixture cohort reaches — its outcome is "
            "constant across all five records — so it is the first thing an acceptance test sees "
            "(§14.0.1).")

    rank = int(np.linalg.matrix_rank(Xn)) if Xn.shape[1] else 0
    if Xn.shape[1] and rank < Xn.shape[1]:
        raise model_fit_error(
            f"O6  the design has rank {rank} of {Xn.shape[1]} columns: {', '.join(columns)}. The "
            "cutpoints ARE the intercepts, so a column of ones is exactly collinear with their sum "
            "and a design carrying its own intercept fails here — which is the opposite of the "
            "penalised logistic fit in this module, and the one place a reader will assume the two "
            "agree.")
```

**O6 does not fire on a width-0 design**, and that is why §4.5's `_assert_exposure_survived` exists on
the caller's side: a matrix with no columns has rank 0 of 0 and is not rank-deficient. `polr` on it
fits the cutpoints alone and returns a `beta` of length zero (§3.1, measured). A module that does not
know what the exposure is cannot know that the missing column was the estimand, so the check belongs
where the specification does — the same division Stage 6 made when it put F5 in `propensity.py` and
not in `model.py`.

**And O6 runs on the design BEFORE the collapse, so the collapsed design's rank is not checked.** This
is stated rather than repaired, and the reason it is safe on the declared path is narrow enough to be
worth writing down:

```
    _assert_polr_fittable(Xn, y, w, …)      O6 sees the FULL design, n rows
    categories = _weighted_categories(y, w)
    keep       = np.isin(y, categories)     drops rows in ZERO-WEIGHT levels only
    Xn         = Xn[keep]                   ← O6 never sees this matrix
```

`keep` is all-True whenever every weight is positive, which is every input either declared caller
gives: [§8]'s overlap weights are strictly interior by Stage 6's F5, and Stage 12 passes unit weights.
So on the declared path the two matrices are the same matrix and the ordering is invisible. It stops
being invisible for **a caller that passes a weight of exactly zero** — which §5.3 explicitly allows,
and which is the whole reason the rule is the positively-weighted set rather than the observed one. If
every row of some design column sits in a zero-weight level, that column becomes constant or the
matrix becomes rank-deficient after the collapse, and `polr` fits it: `np.linalg.solve` raises
`LinAlgError` on an exactly singular Hessian, which is **not** `FitError`, so [§10] would not catch it
(§4.5) and the caller would get a bare `numpy` exception.

**Not repaired, for two reasons, and the second is the one that decides it.** Re-running the rank check
after the collapse would put O6 in two places, which is exactly what §5.3a just took the category rule
out of; and no declared caller can reach it, so the check would be dead code guarding a caller that
does not exist while the *real* residual — a zero-weight caller getting `LinAlgError` instead of
`FitError` — would remain. §15 files it as a limitation of `w = 0` and §14.4's zero-weight-category
test states in its comment that the two fits it compares do **not** share an O6 path, so nobody reads
that assertion as covering this.

### 5.6 No standard error leaves this stage, and that is a correctness claim

`PolrFit` has no `se`, no `cov` and no `bse`; `Primary` has no interval and no p-value. Three reasons,
and the first alone is sufficient:

- **The naive weighted-likelihood covariance is wrong here.** `−H⁻¹` at the optimum is the inverse
  observed information *of a likelihood in which `w_i` counts observations*. The overlap weights are
  not frequencies: they are a tilting function of an **estimated** propensity score, so the sampling
  variability of `e` enters the estimate and `−H⁻¹` omits it entirely. [§7] says so directly — *"the
  bootstrap refits the propensity model in every replicate, so the uncertainty from estimating it
  enters the interval"* — and an interval that omits it is too narrow by an amount nothing here
  bounds.
- **[§10] prescribes the interval that is used**, percentile bootstrap over 2000 replicates, and
  prescribes the primary p-value as a function of the bootstrap distribution of `β`. A second interval
  existing anywhere is a second answer to a question with one prespecified answer.
- **A field is an invitation.** `Fit` carries no standard error either, but nobody would reach for one:
  a propensity model's coefficients are nuisance parameters and Stage 6 §7.3 keeps them in the log
  alone. `β` is the reported effect. A `se` attribute on the object holding it would be in a table
  within two stages, and §14.11's `inspect` assertion on the field list is what stops it.

**What `PolrFit` carries instead** is `iterations`, `converged_on` and the three safeguard counters —
properties of the *iteration*, which is what [§10] needs to report a replicate's fate, and which
`Fit` carries for the same reason (Stage 6 §3.2).

## 6. Separation does not present as failure, and §6 is the only reason this stage has a guard

### 6.1 What was measured

A perfectly separated frame — 20 control at mRS 5, 20 treated at mRS 0, no overlap:

```
    beta         36.4058              exp(beta)  6.469e+15
    alpha        [-18.2029]           iterations 17        converged_on  "score"
    rescales     0                    halvings   0
    min fitted category probability   1.000000    max  1.000000
    max |linear predictor|            18.2029     (POLR_ETA_CLIP is 500)
    max |score|                       2.487e-07   (POLR_SCORE_TOL is 1e-06)
    cond(-H)                          6.854       eigenvalues  [-6.5e-07, -9.0e-08]
```

**Nothing in that block is out of range.** It is not a rank failure: `cond(-H)` is 6.854 — and that
number is a 2×2 matrix's, this construction having two outcome categories, so §6.3 is where it is
compared against anything. It is not a halving exhaustion: zero halvings. It is not a
non-convergence: it converged, on the score criterion, in 17 iterations. Every fitted probability is
finite and in `[0, 1]`. The Hessian's eigenvalues are small in absolute terms — the information has
collapsed — but `cond` cannot see that, because conditioning is a ratio and both eigenvalues collapse
together.

The likelihood approaches its supremum as `β → ∞` and flattens; the score falls below tolerance; the
loop returns. `exp(β)` is `6.469e+15`.

**And it scales with `n` rather than being an artefact of a tiny frame** (§20):

```
    n per arm    10        20        40        80
    beta         34.41     36.41     38.41     38.41
    exp(beta)    8.76e+14  6.47e+15  4.78e+16  4.78e+16
```

**Near separation, which is what a replicate actually produces**, lands in the same place with a
softer landing: one crossover patient out of twenty gives `β = 21.15`, `exp(β) = 1.53e+09`; two gives
`20.42`; three gives `19.97`. Every one of those converged with clean counters.

**So the drop-and-count path never fires.** Four hundred bootstrap-shaped replicates of a 92-record
39-against-53 seven-category frame produced **zero** convergence failures (§20). The roadmap's
**Accept when** and [§10]'s *"replicates whose prespecified fit fails are dropped and counted"* are
both about a failure mode this estimator does not have.

### 6.2 What a percentile interval does with it

[§10] takes `p = 2·min{Pr(β̂* ≤ 0), Pr(β̂* ≥ 0)}` and a percentile interval from the `β̂*`. A
replicate that draws a near-separated sample contributes `β̂* ≈ 21` rather than failing. Two
consequences, and neither is a rounding matter:

- The **percentile interval's upper limit** is a quantile of a distribution with a tail at 1e+09 on the
  odds-ratio scale. At 2000 replicates a handful of such draws move the 97.5th percentile arbitrarily
  far.
- The **p-value is unaffected in the direction that matters and affected in the direction that does
  not**: sign counts are robust to magnitude, so `p` barely moves, while the interval it is supposed to
  agree with "by construction" ([§10]) moves a lot. An interval and a p-value that disagree is the one
  outcome [§10] wrote its p-value definition to prevent.

Dropping the replicate is what [§10] prescribes for a fit that fails. **So the fit must fail.**

### 6.3 The guard, and why it is a bound on `|β|`

```python
def _assert_reportable(fit: model.PolrFit) -> None:
    """G7 — the fitted treatment coefficient is inside the declared bound. §6.

    Stage 6's F5 is the precedent, structurally and in placement: F5 asserts the fitted
    probabilities are strictly interior and lives in `propensity.py` rather than in `model.py`,
    because "a boundary value is a degenerate fit and not a confident one" is a statement about the
    specification and not about the arithmetic. G7 is the same statement for this estimator, and it
    is on the CALLER's side for the same reason.

    It is a bound on beta and NOT on alpha. A cutpoint's magnitude legitimately grows when a
    category is rare — the separated construction's alpha is -18.2 and there is nothing wrong with a
    rare-category cutpoint of that size on a frame that is not separated — while beta is the
    reported effect and exp(beta) is the number a manuscript prints.
    """
    worst = float(np.max(np.abs(fit.beta)))
    if worst >= C.POLR_MAX_ABS_BETA:
        raise model.FitError(
            f"G7  |coefficient| reached {worst:.4g} against a bound of {C.POLR_MAX_ABS_BETA:g}, so "
            f"the common odds ratio is {float(np.exp(worst)):.4g}. This is separation, and it does "
            "NOT present as non-convergence: measured, a perfectly separated frame converges in 17 "
            "iterations on the score criterion with every safeguard counter at zero and every fitted "
            "probability finite [§6.1]. Over 4800 fits at twelve true effect sizes, no "
            "non-degenerate fit exceeded 11.04 at any cutpoint count and no degenerate one came below "
            "18.98; nothing at all lies between them, and every bound from 12 to 20 drops the same "
            "replicates [§6.3]. The "
            "likelihood flattens as the coefficient grows and the loop returns normally. [§10] "
            "drops and counts a replicate whose fit fails; without this raise there is no failure "
            "to count, and the percentile interval is a quantile of a distribution with a tail at "
            f"exp({worst:.1f}) [§6.2]. There is no second estimator [§8, invariant 5].")
```

**Why a magnitude bound and not something more principled.** Four detectors were designed against
§6.1's measurements. **One of them was written up wrongly in an earlier draft of this section, and
running it properly changed both the verdict and the bound.**

| Detector | On a degenerate fit | On a non-degenerate one | Verdict |
|---|---|---|---|
| non-convergence | converges, 17 iterations, score route | converges | **does not fire, ever** |
| `cond(-H)` above a bound | 7.9e7 to 3.8e8 at equal cutpoint count | 29.1 to 162.8 | **discriminates — and is still rejected, below** |
| fitted probability at a boundary | min = max = 1.000000 | max of maxes 0.4726 | fires on *perfect* separation, **misses near separation** (`min p` 0.05 at `β = 21`) |
| `max abs(beta)` | 18.98 to 38.41 | 0.0006 to 11.04 | **a sparse region between them, wide at every cutpoint count** |

**The `cond(-H)` row is a correction.** An earlier draft wrote "6.85 against 3 to 40, does not
discriminate", and both halves were wrong. The 6.854 is real but it is the condition number of a
**2×2** matrix — §6.1's separated construction has two outcome categories and therefore one cutpoint —
while the healthy fits it was compared against are **7×7**. That is not a comparison. Measured properly,
at equal cutpoint count (§20): a 92-record seven-category frame driven to near separation gives
`cond(-H)` of **7.9e7**, exact separation gives **3.8e8**, and 500 healthy seven-category fits give
**29.1 to 162.8**. So it discriminates by six orders of magnitude, and the earlier claim that it did not
was an artefact of comparing matrices of different sizes.

**It is still not the guard, for two reasons that survive the correction.** The first is the same size
dependence that produced the error: `cond(-H)` is a property of a matrix whose dimension is
`len(alpha) + len(beta)`, and §5.3 makes that a property of the sample, so a fixed bound would mean
different things in replicates with different numbers of occupied categories — and on §6.1's own
two-category construction it reads 6.854, *below* every healthy seven-category fit, so a
"cond above a bound" rule fires on nothing there. The second is that a condition-number threshold
cannot be justified against anything a reader can evaluate, where a bound on the reported coefficient
can: `exp(β)` is the number the manuscript prints.

The fitted-probability detector is the one that looks most principled and it fails where it matters:
at one crossover patient the probabilities are unremarkable (`min p = 0.05`) and `exp(β)` is still
`1.5e+09`. Near separation is the case a bootstrap draws, so a detector that only catches *perfect*
separation catches the case that does not arise and misses the case that does.

**The bound's value is prespecified, chosen from the measured band and not from the answer** (§4.3) —
and the band is narrower than an earlier draft of this section claimed. That draft measured 500 fits at
a **single** true effect size of 0.5, got `|β| ∈ [0.0006, 2.1316]`, and concluded the gap was
`(2.13, 19.97)`. A true effect of 0.5 is an assumption about the answer, which §4.3 forbids, so the
measurement was repeated across **twelve** true effect sizes from 0 to 6, 400 samples each — 4800 fits
(§20):

```
    |beta| in [ 0,  2) : 2723        |beta| in [ 8,  9) :    12
    |beta| in [ 2,  4) : 1024        |beta| in [ 9, 10) :     0
    |beta| in [ 4,  6) :  475        |beta| in [10, 18) :     0
    |beta| in [ 6,  7) :   234       |beta| in [18, 20) :   279
    |beta| in [ 7,  8) :   34        |beta| in [20, 25) :    17

    max below 10:  8.7873            min above 10:  18.8055
    values strictly inside (8.79, 18.81):  NONE   <- at SIX cutpoints only; see below
```

**The distribution is bimodal, and an earlier draft called the region between the modes a
"genuinely empty band". IT IS NOT EMPTY, and the correction is item 24.** The 4800 fits above are all
**seven-category** frames — six cutpoints — and §5.3 makes the cutpoint count a property of the sample,
so a [§10] replicate can have fewer. Re-measured 2026-08-24 across every cutpoint count a replicate can
produce, 1800 fits each:

```
    occupied levels     7        6        5        4        3
    cutpoints           6        5        4        3        2
    max |beta| below   8.82    11.04     8.53    11.03     7.63
    min |beta| above  20.67    34.99    19.50    35.33    18.98
    fits in (10, 18)      0        8        0        7         0
```

**At five and three cutpoints, legitimate fits reach 11.04 and land inside the region 8.79-18.81 that
was called empty.** The right description is a **sparse** region whose lower edge moves with the
cutpoint count and whose upper edge is stable near 19-21. Pooled over all four counts, 4800 fits:
**nothing lies between 11.037 and 20.487**, and across individual counts the safe interval is
**(11.04, 18.98)**.

**The bound stays at 14.0, and the reason is that the choice is behaviourally irrelevant across a wide
range.** Pooled over 3 to 6 cutpoints:

```
    bound        8.0     10.0    12.0    14.0    16.0    18.0    20.0    25.0
    dropped     28.6%   26.4%   26.1%   26.1%   26.1%   26.1%   26.1%   10.2%
    largest kept 7.99    9.99   11.04   11.04   11.04   11.04   11.04   21.27
    smallest dropped 8.00 10.01 20.49   20.49   20.49   20.49   20.49   34.97
```

**Every bound from 12 to 20 drops exactly the same replicates**, so 14.0 is a point on a plateau and
not a tuned value. What the re-measurement *does* settle is the other direction: **10.0 is now
positively excluded rather than merely superseded** — it rejects the legitimate fits at 10.01 to 11.04
that only appear at odd cutpoint counts, which the single-category-count calibration could not see.
8.0 bites well into the legitimate distribution; 25.0 admits a degenerate fit at `exp(21.3) = 1.7e9`.

**And the bound is nowhere near the range a stroke trial operates in.** Over 300 samples at each true
effect (§20):

```
    true OR          1.00    1.65    2.72    4.48    7.39
    max |beta|       1.44    2.49    3.28    4.26    5.71
    fraction >= 14      0%      0%      0%      0%      0%
```

The degenerate mode sits at 20.9-21.5 whatever the true effect. So the guard separates two
well-separated populations with a factor of two of clearance on each side, and §12's `test_config.py`
assertion is written against the **corrected** endpoints — `11.04 < POLR_MAX_ABS_BETA < 18.98` — so a
future edit back to 10.0 fails there.

**This also strengthens §6.3's rejection of `cond(-H)`.** That detector was rejected for depending on
the matrix dimension, which is `len(alpha) + len(beta)` and therefore a property of the sample. The
K-dependence just measured is the same defect showing up in the magnitude bound's own calibration —
which is survivable there, because the bound is flat and the plateau is wide, and would not be for a
condition-number threshold whose scale moves with dimension.

**And a legitimate fit can reach `|β| = 8.79`, which is an odds ratio of 6570.** That is not a
degenerate fit and this guard does not reject it. Whether a common odds ratio of that size is
*credible* in this cohort is a question for the PI and for [§9]'s balance diagnostics, not for a
runtime assertion: §16 records that this stage does not decide it, and the bound is set to catch
degeneracy rather than to adjudicate plausibility.

It lives in `config.py` and not in a default argument, for the reason `FIRTH_TOL` does
(`config.py:281-284`): [§10] refits this model in every one of `N_BOOT` replicates, so a bound is a
property of the sampling distribution and not a runtime knob. Changing it changes which replicates are
dropped, which is a [§13] amendment.

### 6.4 What this stage does not do about it

**It does not penalise.** [§8] prescribes an unpenalised weighted proportional-odds model and names no
penalty, where [§7] names Firth explicitly for the propensity model. A Firth-penalised ordinal fit
would keep sparse replicates finite and reduce the drop count, and it would also change the reported
quantity from a common odds ratio to a penalised one and change the null distribution [§10]'s p-value
is read against. That is a PI decision and an [§8] amendment, not a code choice, and §17 records it as
out of scope rather than as unconsidered.

**It does not trim, clip or winsorise `β`.** A clipped coefficient is a reported effect that is a
function of the bound, and [§10] would take percentiles of a distribution with an atom at 10.

**It does not report the drop rate.** Counting dropped replicates is [§10]'s, and §11 hands over the
measurement that Stage 10 will need to size it: on synthetic sparse replicates of this cohort's shape
the guard is what produces the drops, so Stage 10's failure counter will read a number that is zero
today only because nothing has fired it.

## 7. The orientation

### 7.1 The rule

```
    exp(beta) > 1  FAVOURS BRIDGING.

    Under §5.1's `logit P(Y <= k) = alpha_k + beta*A`, a positive beta raises P(Y <= k) at every
    k in the treated arm: mass moves to LOW mRS, which is a better outcome. So `> 1 favours
    bridging` requires `higher_is_better` to be False, and that is asserted rather than assumed.
```

`OUTCOMES[PRIMARY_OUTCOME].higher_is_better` is `False` for `mrs_90d` (`config.py:410-411`). The
orientation sentence in §9.2's audit entry is emitted **from that field**, and `primary` asserts it:

```python
def _orientation() -> str:
    """The [§8] orientation sentence, computed from the outcome registry. §7.1.

    Not a literal. If an ordinal outcome for which higher is better were ever registered as primary,
    a hardcoded sentence would print `> 1 favours bridging` beside a coefficient meaning the
    opposite, and nothing anywhere would disagree with it. The assertion is the point; the returned
    string is a convenience.
    """
    outcome = C.OUTCOMES[C.PRIMARY_OUTCOME]
    if outcome.kind != "ordinal" or outcome.higher_is_better:
        raise C.SchemaError(
            f"the [§5] primary outcome {C.PRIMARY_OUTCOME!r} is kind={outcome.kind!r} with "
            f"higher_is_better={outcome.higher_is_better}. §5.1's parametrisation gives "
            "`exp(beta) > 1 favours bridging` only for an ordinal outcome on which LOWER is better; "
            "for one on which higher is better the same coefficient favours the comparator, and the "
            "orientation statement [§8] requires would be exactly inverted with every number in the "
            "table still finite and still plausible.")
    return (f"exp(beta) is the common odds ratio and is oriented so that > 1 favours "
            f"{C.TREATMENT_LABELS[_TREATED]}: under logit P(Y <= k) = alpha_k + beta*A a positive beta "
            f"raises P(Y <= k) in the treated arm at every k, and lower {outcome.label} is better.")
```

### 7.2 Every reference implementation has the other sign, and one of them is in this repository

```
    this stage      logit P(Y <= k)  =  alpha_k  +  x'beta         [§14a]'s own form
    statsmodels     logit P(Y <= k)  =  theta_k  -  x'beta         OrderedModel
    R  polr / clm   logit P(Y <= k)  =  zeta_k   -  x'beta
```

**Measured (§20), and it is the cleanest measurement in this document.** On a 500-record two-covariate
synthetic frame, unweighted:

```
    ours.beta        [ 0.97283069,  0.58382310]
    statsmodels      [-0.97284722, -0.58383309]      ours.beta + sm.beta  ->  1.65e-05
    ours.alpha       [-1.54340634, -0.52258248,  0.37686850,  1.76824648]
    statsmodels      [-1.54343790, -0.52259457,  0.37684297,  1.76826674]   diff  3.16e-05
```

**The coefficients are negated and the cutpoints are not.** That asymmetry is what makes the mistake
survive review: a reader comparing the two fits sees four cutpoints agreeing to five decimals and
concludes the implementations agree.

**And `pilots/analysis.py:475` writes `float(np.exp(-res.params.iloc[0]))`**, with the comment *"A
positive coefficient means a shift towards higher (worse) mRS; invert so the common OR is oriented
like every other outcome"*. That is **correct for statsmodels' parametrisation and wrong for this
one.** The negation is not a convention the pilots chose; it is the conversion between two
parametrisations, and lifting the line without lifting the parametrisation inverts the primary result
of the study. §18 records it as the single most consequential do-not-lift in this stage, and §14.5
pins it with a companion that applies the pilots' negation to our `β` and asserts the direction
reverses.

### 7.3 How it is tested, since the criterion cannot fail on its own

The roadmap's criterion is *"the orientation test confirms `exp(β) > 1` when treatment shifts mRS
downward"*. On a construction where treatment shifts mRS downward, an implementation with the wrong
sign returns `exp(β) < 1` and the criterion catches it — so unlike §5.4's monotonicity and §8's
ordering, this one is genuinely falsifiable as written. What it needs is the construction, and a
companion so that the test is known to fire.

Measured (§20), on a 200-record construction with the treated arm's mass at mRS 0-3 and the control's
at 2-6 (mean category 1.11 against 3.82):

```
    beta        4.437745        exp(beta)  84.584014     ->  > 1, correct
    the alpha_k - beta*A form would report exp(beta)  0.011823
```

§14.5 asserts both numbers. The second is what makes the first a measurement: a test that only checks
`exp(β) > 1` on data where the effect is large passes for an implementation that returns
`exp(|β|)`.

## 8. The cumulative risk differences [§8]

### 8.1 The one weighted per-arm proportion

```python
def weighted_proportion(x: np.ndarray, a: np.ndarray, w: np.ndarray) -> dict[int, float]:
    """The weighted mean of `x` within each declared arm. Keyed by arm code, as TREATMENT_LABELS is.

    Public, and it is Stage 9's as much as this stage's — the [§8] weighted risk difference is this
    function differenced and the weighted marginal odds ratio is these two proportions — which is
    why it lives here rather than inside `cumulative_rd`. A second implementation in a later stage
    would be a second definition of the denominator, and [§11] requires the denominator. The same
    move Stage 6 §3 made with `ess` and Stage 7 §0.1 made with `smd`.

    `nan` where an arm carries no positive weight, never 0.0: `pilots/analysis.py:287-289` returns
    nan here too and it is the one guard in that file to lift unchanged. The check comes BEFORE the
    mean, because `np.average` with weights summing to zero RAISES ZeroDivisionError rather than
    returning nan [Stage 7 §3.1] — the order is load-bearing, and it is load-bearing in the same
    place one stage on.
    """
    out: dict[int, float] = {}
    for code in C.TREATMENT_LABELS:
        arm = a == float(code)
        total = float(w[arm].sum())
        out[code] = float(np.average(x[arm], weights=w[arm])) if total > 0.0 else np.nan
    return out
```

```python
def cumulative_rd(y: np.ndarray, a: np.ndarray,
                  w: np.ndarray) -> tuple[dict[int, float], dict[int, dict[int, float]]]:
    """RD_k and the two weighted cumulative distributions, over MRS_THRESHOLDS. [§8]

    Takes no threshold list: MRS_THRESHOLDS is prespecified at Stage 1 for this stage, and a
    keyword would make a prespecified choice look like an option [Stage 7 §15, Stage 6 §6.4].

    Returns BOTH the differences and the two distributions they came from, from one call, for
    Stage 7 §3's reason: `RD_k` alone cannot be checked against anything, and a reader given only
    the difference cannot see that P(Y <= k) reached 1 three thresholds early. §9.3's table prints
    all three columns because of this signature, not despite it.

    Each threshold's indicator carries the response's own missingness — a comparison against a
    missing value is False in numpy, which would silently count an absent outcome as a non-event
    [roadmap Stage 3] — and there is no missingness to carry on this path, because §4.1's mask
    removed it. The mask is the guard; this is the statement that it is.

    `share[TREATED] - share[COMPARATOR]` and not `share[1] - share[0]`: the arm codes are named once at
    module scope from TREATMENT_LABELS and read here, because this is the one line in the stage whose
    SIGN is a function of which literal goes first, and §7 is four hundred words about that sign. Two
    bare integers in a subtraction are the one place a reader checking the orientation has nothing to
    check against.
    """
    rd: dict[int, float] = {}
    cumulative: dict[int, dict[int, float]] = {}
    for k in C.MRS_THRESHOLDS:
        share = weighted_proportion((y <= float(k)).astype(float), a, w)
        cumulative[k] = share
        rd[k] = share[_TREATED] - share[_COMPARATOR]
    return rd, cumulative
```

**`_TREATED` and `_COMPARATOR` are two more module constants**, beside §9.1's three, and they are
derived rather than written:

```python
# outcome.py — the arm codes, named once. §8.1
#
# [§8]'s RD_k is P(. | bridging) - P(. | EVT alone) and the ORDER OF THE SUBTRACTION IS THE SIGN, which
# §7 is four hundred words about. `share[1] - share[0]` is correct and is two bare integers in the one
# expression where a reader checking the orientation most needs a name to check against. Derived from
# TREATMENT_LABELS rather than written as 1 and 0, for PRIMARY_OUTCOME's reason (§12): the registry
# already says which code is the treated arm, and a literal here is a second place for it to be said.
_COMPARATOR: Final[int] = min(C.TREATMENT_LABELS)
_TREATED: Final[int] = max(C.TREATMENT_LABELS)
```

`min` and `max` and not `0` and `1`: `config.py:316` documents the coding as `1 = IVT before EVT
(bridging), 0 = EVT alone`, so the treated arm is the higher code, and that is the property being read.
§14.6 asserts `_TREATED` is the key whose label is `TREATMENT_LABELS[1]` and that the two constants are
distinct, which is what fails if a third arm code is ever declared — at which point a two-arm
subtraction is the wrong shape and should fail rather than silently contrast the extremes.

### 8.2 From the weighted distributions, and never from six threshold models

[§8]'s sentence is *"taken from the weighted empirical distributions — ordered by construction, so no
monotonicity problem arises"*, and the reason it gives is that *"augmenting each threshold separately
would replace one primary number with six and can produce crossing cumulative probabilities"*.

**Ordering is a property of the construction, so the roadmap's monotonicity criterion cannot fail.**
`P_w(Y ≤ k | arm)` is a cumulative sum of non-negative weights over a denominator that does not depend
on `k`, so it is non-decreasing in `k` arithmetically. Measured (§20): 2000 samples of a 40-record
seven-category frame gave **0** non-monotone results by this route. An acceptance test asserting
monotonicity therefore asserts something that is true of the code, of any correct alternative
implementation, **and of a wrong one that reached the same numbers a different way** — which is why
§21's roadmap amendment 3 replaces the criterion rather than restating it.

**What must be asserted instead is that the route taken is this one.** §14.6 does it two ways: the
`RD_k` are asserted equal to a hand-computed weighted proportion difference on a frame small enough to
verify by hand, and a companion computes the same six numbers through six weighted threshold models
and asserts they **differ**. Different, not crossing — and that is the honest form, because:

**The crossing the SAP warns about is real and is rare, and the measurement says so** (§20). Over 2000
samples of a 40-record seven-category frame, the threshold-model route produced a strict decrease in 8
of the 1990 that fitted — 0.40% — and the largest decrease measured was below 5e-07, which is float
noise rather than a crossing anyone would see. **So the argument for the empirical route rests on
[§8]'s reasoning and on construction, not on a measured failure rate**, and an earlier draft of this
section that claimed 40 of 40 samples crossed was reporting a companion with a complementary-
probability bug: it computed `P(ind ≤ 0)`, which is `1 − P(Y ≤ k)`, so every sequence was decreasing
by construction. The claim is struck and the corrected rate is recorded.

### 8.3 What the thresholds are, and where they stop

`MRS_THRESHOLDS` is `(0, 1, 2, 3, 4, 5)` — six thresholds over seven declared levels — and §12 adds
`MRS_LEVELS` beside it with the assertion `MRS_LEVELS[:-1] == MRS_THRESHOLDS`, so the two cannot
drift.

**`RD_6` is structurally zero and is therefore not reported.** `P(Y ≤ 6) = 1` in both arms by
definition, so `RD_6 = 0` on every input. Measured (§20) on a construction with deaths in both arms:
`RD_5 = +0.049027` and `RD_6 = +0.000000` exactly. `RD_5` is the survival contrast —
`P(alive | bridging) − P(alive | EVT alone)` — and is **not** structurally zero; a reader who assumes
the last reported threshold is the trivial one has the wrong one.

**Ties in `RD_k` are expected and are not a defect.** Where a declared level is unoccupied, two
adjacent thresholds have the same cumulative probability and the same `RD`. §14.0.2's golden vector
has two such pairs, from a cohort missing mRS 2 and 5, and the assertion is non-decreasing rather than
increasing. An implementation asserting strict monotonicity would fail on a frame that is entirely
correct.

## 9. The three `model` audit entries

### 9.1 No new kind, and that is `data.py`'s own declaration honoured

`data.py:143-152` declares `model` and states what it is for: *"`model` is Stage 6's and describes
**what was fitted** to the population `cohort` describes — a claim about the data, with convergence
properties, which renders under no existing heading"*, and `data.py:151-152` adds that *"Stages 7-13
all render under `model` or under an existing kind"*. Stage 7 §7.1 tested that claim for the first
time and left `data.py` untouched; Stage 8 tests it a second time and reaches the same answer.

A weighted proportional-odds fit is a fitted model, its coefficient is the estimate, and the
cumulative risk differences are that estimate on an absolute scale. All three belong under
`## Fitted models`, after `overlap_weights` and after Stage 7's two diagnostic entries, which is the
order the document reads in: what was fitted to the population, whether it balanced, and what it
estimated.

**So `KINDS` stays nine, `_HEADINGS` gains nothing, and none of the literal pins in the six existing
test modules moves.** §14.9 asserts it, which is what fails if someone adds an `estimate` kind.

```python
# outcome.py — three of the five module constants §3 imports `Final` for; §8.1 declares the other two
_STEP_COMPLETENESS: Final[str] = "outcome_completeness"
_STEP_FIT: Final[str] = "primary_fit"
_STEP_RD: Final[str] = "cumulative_rd"
```

Named as constants rather than written at the three `audit.record` call sites because §14.9 asserts
the entries appear in this order at positions captured before the call, and a test comparing against a
string literal it also writes is a test of nothing.

### 9.2 The three entries

| # | step | `n` | Contents |
|---|---|---|---|
| 1 | `outcome_completeness` | records the outcome mask removed from `in_model` | the [§11] denominator chain: one row per population step with its count, and the removed records named |
| 2 | `primary_fit` | records estimated | one row per **declared** mRS level marking its cutpoint, its absence, or that it is the top of the **fitted** scale, then the treatment coefficient and `exp(β)`; the orientation sentence and the no-standard-error sentence in `detail` |
| 3 | `cumulative_rd` | `len(MRS_THRESHOLDS)` | one row per declared threshold: the two weighted cumulative probabilities and `RD_k` |

Four rules about the tables, all inherited:

- **Every declared thing is rendered whether or not the data fills it.** Entry 2 has a row for every
  `MRS_LEVELS` entry, including a level nobody is in, because a category that vanishes from the table
  is a category nobody knows was considered (Stage 6 §7.5, Stage 7 §4.2). Entry 3 has a row for every
  `MRS_THRESHOLDS` entry.
- **Every cell is computed from the object it describes, never by subtracting another row.** Entry 3's
  `RD_k` column is `Primary.rd[k]` and not the difference of the two printed cells — a table that
  reconciles with itself by construction has stopped being evidence (Stage 5 §7.3).
- **`_fmt` is the only float formatter** and `missing` the only rendering of an absent value.
- **No cell of any table contains a `|`.** `data._md_table` does no escaping and sizes its separator
  from `len(rows[0])` (`data.py:242-245`), so one pipe in a header cell gives a header row with more
  markdown cells than its body — and **byte-identity across hash seeds does not catch it**, because a
  table identically broken under both seeds is still identical. Stage 7 §21.1 found this by reading a
  cell back out of `to_markdown()`; this stage writes its headers accordingly from the start, and
  §14.9 asserts the character's absence. The column that would have carried one is entry 3's
  cumulative probability, whose natural header is `P(Y <= k | bridging)` — no pipe, and the
  mathematical bar is deliberately not used.

```python
def _completeness_table(df: pd.DataFrame, ps: propensity.Propensity,
                        in_estimate: pd.Series) -> tuple[tuple[str, ...], ...]:
    """The [§11] denominator chain, one row per population step. §4.1.

    Four rows, and THEY ARE NOT A CUMULATIVE CHAIN — row 3 is a count over the whole cohort, as its
    own third cell says, so adjacent differences are not restriction costs. On the workbook the four
    read 93, 92, 93, 92 and the differences are −1, +1, −1: a "restriction" that adds a record. An
    earlier draft called them cumulative and invited exactly the subtraction that does not work
    (§21.1 item 15). `cohort_flow` (Stage 5 §7.3) IS a nested chain; this is a chain of two with a
    parallel denominator printed beside it, because [§11] asks for the outcome's own completeness as
    well as the estimate's.

    Row 3 is kept and not made nested, because the quantity a reader wants when the last row is short
    is "how many outcomes are present at all", and deriving it by subtraction from a nested row is
    the reconciles-with-itself table Stage 5 §7.3 forbids. The last row is what the estimate below it
    is computed on, and it is the number [§11] asks to be reported for every estimate.
    """
    present = df[C.PRIMARY_OUTCOME].notna()
    header = ("population", "n", "what it is")
    return (
        header,
        ("[§3] primary cohort", str(len(df)), "both restrictions applied [Stage 5 §11]"),
        ("in_model", str(int(ps.in_model.sum())),
         "covariate-complete; the ATO population [Stage 6 §4.4]"),
        (f"{C.PRIMARY_OUTCOME} present (whole cohort)", str(int(present.sum())),
         "NOT a step in the chain — the outcome's own completeness, before the mask below"),
        ("in_estimate", str(int(in_estimate.sum())),
         "in_model AND outcome present — this estimate's [§11] denominator"),
    )
```

```python
def _record_outcome_completeness(df: pd.DataFrame, ps: propensity.Propensity,
                                 in_estimate: pd.Series, audit: Audit) -> None:
    """Entry 1, and it must name as many patients as it says it removed.

    `propensity._record_exclusion` (propensity.py:400-432) line for line, both mechanisms included,
    because both are what make the comparison able to fail at all: the identifiers go through a
    `set`, so a duplicated case_id collapses to one name, and a missing case_id is dropped by
    `notna()` rather than stringified. Without either, `len(case_ids)` equals `n` by construction and
    the check is dead code.

    The raise comes BEFORE `audit.record`, again as Stages 5 and 6 do it: recording first and reading
    the entry back leaves an unreconcilable entry in an `Audit` that has no removal path.

    On v7 this entry's `n` is 0 and it names nobody (§4.1, measured), and it is recorded anyway. An
    entry that appears only when it has something to say is an entry whose absence a reader has to
    interpret, and [§11] asks for the denominator of every estimate rather than for the denominators
    that happen to differ.
    """
    removed = ps.in_model & ~in_estimate
    excluded = df.loc[removed & df["case_id"].notna(), "case_id"]
    case_ids = tuple(sorted(set(excluded)))
    n_removed = int(removed.sum())
    if len(case_ids) != n_removed:
        raise C.SchemaError(
            f"{_STEP_COMPLETENESS} removes {n_removed} record(s) from the ATO population and can "
            f"name {len(case_ids)}. An estimate whose denominator the log cannot reconstruct is an "
            "estimate nobody can check [§11].")
    audit.record("model", _STEP_COMPLETENESS, n_removed,
                 f"[§11] complete-case per estimate. {int(in_estimate.sum())} of "
                 f"{int(ps.in_model.sum())} weighted record(s) carry the [§5] primary outcome and "
                 f"are estimated; {n_removed} do not, and every one of them is named above. Such a "
                 "record keeps its weight and loses its estimate — it is in the ATO population "
                 "[Stage 6 §4.4] and not in this estimate's denominator, which is why this is a "
                 "THIRD population and not a restatement of that one [§4.1]. The primary estimate "
                 "and the cumulative risk differences below are computed on the same mask, bound "
                 "once, because [§8] presents them as one contrast on two scales.",
                 case_ids=case_ids, table=_completeness_table(df, ps, in_estimate))
```

```python
def _fit_table(fit: model.PolrFit) -> tuple[tuple[str, ...], ...]:
    """One row per DECLARED mRS level, then the coefficient and the odds ratio. §9.2.

    Ranges over MRS_LEVELS and never over `fit.categories`, for `propensity._design_table`'s reason
    (propensity.py:242-261): a level absent from the fitted set is a fact about the population, and a
    row that vanishes is a fact nobody sees. The `status` column is where the collapse of §5.3
    becomes visible in the log rather than inferable from a parameter count.

    **The no-cutpoint row is keyed on `fit.categories[-1]` and NOT on `MRS_LEVELS[-1]`, and that is a
    correction.** There are two different reasons a declared level has no cutpoint — it is absent from
    the population, or it is the TOP of the fitted scale, P(Y <= it) = 1 — and only the second is a
    property of the fit. An earlier draft tested `level == MRS_LEVELS[-1]` for the second, which is
    right only while the top fitted category IS the top declared one. It is, on every frame this stage
    has: the workbook occupies all seven levels and `ordinal_cohort()`'s top category is 6, which is
    `MRS_LEVELS[-1]` by coincidence. It is NOT, in any [§10] replicate that draws no patient at the
    worst mRS — and mRS 6 is death, so that replicate is ordinary rather than exotic. On
    `truncated_cohort()` (§14.0.5), whose categories are (0, 1, 2, 3, 4), the draft rendered — measured,
    by running both branches:

        alpha at mRS <= 4 | 4 | missing | fitted     ← no cutpoint above the top FITTED level, so
                                                       `fitted` is empty there and the estimate
                                                       column reads `missing` under status `fitted`
        alpha at mRS <= 6 | 6 | missing | highest declared level — no cutpoint
                                                     ← true, and it HIDES that nobody is at mRS 6

    Both rows are wrong and neither raises. §14.9's assertion — absent rows render `missing` with the
    absence status — **passes on that table**, because the broken row is not an absent one. §14.0.5 is
    the frame that witnesses it, and §14.9 asserts both the top fitted level's status by name and the
    invariant that generalises it: a row reads `missing` if and only if its status says there is no
    cutpoint. Verified on all three frames under the branch below and violated on `truncated_cohort()`
    under the draft's.
    """
    fitted = {level: value for level, value in zip(fit.categories[:-1], fit.alpha)}
    header = ("parameter", "mRS level", "estimate", "status")
    rows: list[tuple[str, ...]] = []
    for level in C.MRS_LEVELS:
        if level not in fit.categories:
            status = ("absent from this population — no cutpoint [§5.3]"
                      + (", and it is the highest declared level"
                         if level == C.MRS_LEVELS[-1] else ""))
        elif level == fit.categories[-1]:
            status = "highest FITTED level — no cutpoint, P(Y <= it) = 1"
        else:
            status = "fitted"
        rows.append((f"alpha at mRS <= {level}", str(level),
                     _fmt(fitted.get(level)), status))
    for name, value in zip(fit.columns, fit.beta):
        rows.append((name, "—", _fmt(value), "treatment coefficient [§8]"))
        rows.append((f"exp({name})", "—", _fmt(float(np.exp(value))),
                     "common odds ratio — see the orientation above"))
    return (header, *rows)
```

```python
def _fit_detail(fit: model.PolrFit, n: int) -> str:
    absent = tuple(level for level in C.MRS_LEVELS if level not in fit.categories)
    return (
        f"[§8] weighted proportional-odds model with treatment as the sole predictor, over {n} "
        f"record(s). Converged in {fit.iterations} iteration(s) on the {fit.converged_on} "
        f"criterion; {fit.rescales} step(s) were shortened by the trust region and {fit.halvings} "
        f"halving(s) taken. {len(fit.categories)} outcome category/ies carried positive weight, so "
        f"{len(fit.alpha)} cutpoint(s) were estimated"
        + (f"; {len(absent)} declared level(s) are absent from this population and are marked as "
           f"such below: {', '.join(str(level) for level in absent)}. " if absent else ". ")
        + _orientation()
        + " NO STANDARD ERROR AND NO INTERVAL IS PRODUCED HERE, and that is a property of the "
        "estimator rather than a division of labour: the observation weights are a tilting function "
        "of an ESTIMATED propensity score and not frequencies, so the inverse observed information "
        "omits the variability of estimating it [§7]. [§10]'s percentile bootstrap over N_BOOT "
        "replicates is the prespecified interval and the source of the primary p-value, which is a "
        "test of this coefficient and NOT a risk-difference-scale test [§10]. The estimate is "
        "marginal in the overlap population [§8] and is conditional on the [§7] propensity "
        "specification, which defines that population [§7].")
```

```python
def _rd_table(rd: dict[int, float],
              cumulative: dict[int, dict[int, float]]) -> tuple[tuple[str, ...], ...]:
    """One row per DECLARED threshold: the two weighted cumulative probabilities and RD_k. §8

    No cell and no header contains a `|` — the columns are named `P(Y <= k | …)` in mathematics and
    `P(Y <= k), <arm>` here, because `data._md_table` does no escaping and a pipe in a header cell
    gives a header row with more markdown cells than its body [Stage 7 §21.1, §14.9].

    The RD column is read from `rd` and not computed from the two printed cells: a table that
    reconciles with itself by construction has stopped being evidence [Stage 5 §7.3].
    """
    header = ("threshold", *(f"P(Y <= k), {label}" for label in C.TREATMENT_LABELS.values()),
              "RD_k")
    rows = tuple(
        (f"mRS <= {k}", *(_fmt(cumulative[k][code]) for code in C.TREATMENT_LABELS),
         _fmt(rd[k]))
        for k in C.MRS_THRESHOLDS)
    return (header, *rows)
```

```python
def _rd_detail(totals: dict[int, float], n: int) -> str:
    return (
        f"[§8] weighted empirical cumulative risk differences, RD_k = P(mRS <= k | "
        f"{C.TREATMENT_LABELS[_TREATED]}) - P(mRS <= k | {C.TREATMENT_LABELS[_COMPARATOR]}), over "
        f"{len(C.MRS_THRESHOLDS)} declared threshold(s) and {n} record(s) — the same records and the "
        "same weights as the coefficient above [§4.2]. Taken from the weighted DISTRIBUTIONS and "
        "not from threshold-specific models, so the cumulative probabilities are non-decreasing in k "
        "by construction rather than by luck [§8]. Ties between adjacent thresholds are expected "
        "wherever a declared level is unoccupied and are not a defect [§8.3]. Per-arm weight totals, "
        "which are these proportions' denominators [§11]: "
        + ", ".join(f"{C.TREATMENT_LABELS[code]} {_fmt(total)}" for code, total in totals.items())
        + f". mRS <= {C.MRS_LEVELS[-1]} is not reported because P(Y <= {C.MRS_LEVELS[-1]}) is 1 in "
        f"both arms by definition, so its risk difference is structurally zero; mRS <= "
        f"{C.MRS_THRESHOLDS[-1]} is the survival contrast and is not [§8.3]. THESE CARRY CONFIDENCE "
        "INTERVALS BUT NO p-VALUES [§8]: they are the absolute-scale presentation of one contrast, "
        "and six threshold-wise tests of it would only invite selection of the most favourable cut. "
        "The primary outcome has exactly one test [§8, §10].")
```

### 9.3 What the entries deliberately do not carry

- **No interval, no p-value, no standard error**, in any of the three. §5.6.
- **No verdict.** [§8] specifies no threshold on `β` and Stage 8 declares none, so there is no
  `significant` column and nothing analogous to Stage 7's `abs SMD < 0.1`. The one bound this stage
  has is §6's, and a fit that reaches it raises rather than being reported with a flag.
- **No comparison against the [§14] estimates.** [§16] requires that wherever a [§14] estimate appears
  beside the primary one, the different-populations statement appears with it. That is Stage 14's
  guardrail over a table this stage does not build, and putting a hedge about Stage 12 into Stage 8's
  log would put it in the one place Stage 12's output never appears.
- **No arm-level mRS distribution beyond the cumulative probabilities.** The seven-category
  distribution per arm is [§16]'s baseline-table material and Stage 14's; the cumulative form is what
  [§8] asks for and what `RD_k` is computed from.

## 10. `primary`, written out

```python
def primary(df: pd.DataFrame, ps: propensity.Propensity, audit: Audit) -> Primary:
    """The [§8] primary outcome estimate over the [§3] cohort and the [§7] weights.

    Takes no outcome name: [§8] rests on there being ONE primary quantity with one test, and
    PRIMARY_OUTCOME is computed from the [§5] registry (§12). Takes no covariate list: treatment is
    the sole predictor by specification. Takes no threshold list: MRS_THRESHOLDS is Stage 1's. A
    keyword for any of the three would make a prespecified choice look like an option
    [Stage 6 §6.4, Stage 7 §15].

    Returns a Primary; adds no column to `df` and edits nothing (§0.2). Appends three `model`
    entries.

    Raises SchemaError on G1-G5 and on D1-D4, and model.FitError on G6, G7, O1-O6 or a fit that does
    not converge. [§10] catches the second to drop and count a replicate, and may catch nothing else:
    a SchemaError here is a bug in the resampler, not a sparse replicate (propensity.py:468-469).

    Six things about the order below, each of which is a failure if moved:

      * `_assert_primary_inputs` comes FIRST, in two phases, before the mask is built for anything
        but its own phase-2 checks (§4.4a).
      * `in_estimate` is bound ONCE and every subsequent line reads it — the design, the response,
        the weights, the arm vector and the cumulative distributions. Two separate reads would be two
        chances to compute the coefficient and the risk differences over different populations, which
        [§8] presents as one contrast (§4.2).
      * `_record_outcome_completeness` comes BEFORE the design, so the log names the [§11]
        denominator even if `design` then raises on D2 — the exclusion is a fact about the frame,
        established before any modelling choice, and a log that names it is useful precisely when
        what follows fails (propensity.py:440-443).
      * `_assert_exposure_survived` sits BETWEEN the design and the fit, because a width-0 design
        FITS and the estimand's absence stops being visible one line later (§4.5).
      * `_assert_reportable` sits BETWEEN the fit and everything that reads it, because a separated
        fit's `exp(beta)` is a finite float that every downstream table will accept (§6).
      * `cumulative_rd` runs AFTER the fit rather than before, so that a frame which cannot support
        the primary quantity does not produce the absolute-scale presentation of a contrast that has
        no estimate. [§8] makes RD_k a presentation of `beta`, not an alternative to it.
    """
    _assert_primary_inputs(df, ps)                                     # G1-G3 | G4-G5

    in_estimate = ps.in_model & df[C.PRIMARY_OUTCOME].notna()          # §4.1
    _record_outcome_completeness(df, ps, in_estimate, audit)           # entry 1

    sub = df.loc[in_estimate]
    X, dropped = model.design(sub, (C.TREATMENT,))                     # D1-D4 inherited
    _assert_exposure_survived(X, dropped)                              # G6, §4.5

    y = sub[C.PRIMARY_OUTCOME].to_numpy(dtype=float)
    a = sub[C.TREATMENT].to_numpy(dtype=float)
    w = ps.w.loc[in_estimate].to_numpy(dtype=float)

    fit = model.polr(X, y, w)                                          # raises FitError; §5
    _assert_reportable(fit)                                            # G7, §6

    n = int(in_estimate.sum())
    audit.record("model", _STEP_FIT, n, _fit_detail(fit, n), table=_fit_table(fit))

    rd, cumulative = cumulative_rd(y, a, w)                            # §8
    totals = {code: float(w[a == float(code)].sum()) for code in C.TREATMENT_LABELS}
    audit.record("model", _STEP_RD, len(C.MRS_THRESHOLDS), _rd_detail(totals, n),
                 table=_rd_table(rd, cumulative))

    beta = float(fit.beta[0])
    return Primary(
        beta=beta,
        odds_ratio=float(np.exp(beta)),
        alpha=tuple(float(v) for v in fit.alpha),
        cut_levels=tuple(int(level) for level in fit.categories[:-1]),
        rd=rd,
        cumulative=cumulative,
        in_estimate=in_estimate,
        fit=fit)
```

**`fit.beta[0]` and not `fit.beta.sum()`, and G6 is what makes the index safe.** A width-0 design gives
`beta` of length zero, so `beta[0]` would be an `IndexError` and `beta.sum()` would be `0.0` — an odds
ratio of exactly 1.0, which reads as "no effect" and is the estimand having vanished (§4.5, measured).
G6 raises before either can happen, and §14.11 asserts the ordering by removing G6 and showing which
of the two a caller gets.

**`cut_levels` is `fit.categories[:-1]` and that is the identity `alpha` needs.** `alpha[j]` is the
cutpoint between `categories[j]` and `categories[j+1]`, so the level it cuts *at or below* is
`categories[j]`. On a collapsed fit those are not `MRS_THRESHOLDS`: a cohort missing mRS 2 has a
cutpoint at `mRS <= 1` whose next category is 3, and reporting it against `MRS_THRESHOLDS[1]` would
label a cutpoint with a threshold it is not.

## 11. Data flow into Stages 9-14

```
  primary(cohort, ps, audit) returns Primary
   ├─ beta          the [§8] treatment coefficient, alpha_k + beta*A form      [§5.1]
   ├─ odds_ratio    exp(beta), > 1 favours bridging                            [§7]
   ├─ alpha         the fitted cutpoints, ascending                            [§5.2]
   ├─ cut_levels    the declared mRS level each cutpoint cuts at or below      [§5.3]
   ├─ rd            RD_k over MRS_THRESHOLDS                                   [§8]
   ├─ cumulative    P(Y <= k | arm), the distributions rd came from            [§8]
   ├─ in_estimate   boolean, TOTAL — the [§11] denominator                     [§4.1]
   └─ fit           PolrFit: categories, iterations, converged_on, 3 counters  [§3]

  the cohort frame and the Propensity come back unchanged                      [§0.2]

  audit — the same object load() created:
    load 7 / derive 4 / classify 1 / build 6 / fit 4 / assess 2 / primary 3  =  27
```

What each later stage may rely on, and what each owes:

- **Stage 9 [§8]** — reads `weighted_proportion` and nothing else from this stage. Its weighted risk
  difference is that function differenced, its marginal odds ratio is the same two proportions, and its
  augmentation tilts by `h = e(1 − e)`, which is `Propensity`'s and not this stage's. It does **not**
  read `Primary`: the binary estimators are their own estimates on their own [§11] denominators, and
  the primary ordinal quantity is not an input to any of them.
- **Stage 10 [§10]** — refits **this whole function** in every one of `N_BOOT` replicates, and four
  things are owed to it explicitly:
  - **It is deterministic and writes no file.** No seed is held here, no clock is read, and the only
    mutable object touched is the `Audit` a caller passes in. **A replicate loop must pass its own
    `Audit`, and "or accept 3 entries per replicate" understates it**: at `N_BOOT` = 2000 that is 6000
    entries, each carrying a rendered table of 3 to 10 rows and a `detail` string of several hundred
    characters, and one of the three carries a `case_ids` tuple. Reusing the point fit's `Audit` would
    make the log unreadable and the process's memory a function of `N_BOOT`. A throwaway `Audit` per
    replicate is the cheap answer; a `record`-suppressing mode is not offered, because a stage that can
    be asked not to log is a stage whose log is optional.
  - **`Σw` inside a replicate is not measured here and the tolerances are absolute in it.** §5.2 records
    `Σw = 27.736623` on the point estimate's population and argues that `POLR_TOL` and
    `POLR_SCORE_TOL` are comfortable at that scale and at Stage 12's. A stratified resample of 92
    records has a similar `Σw` by construction, so nothing is expected to move — but "expected" is not
    "measured", and Stage 10 is the first place it can be. If a replicate's `Σw` came in an order of
    magnitude low, the same absolute tolerance would be an order of magnitude looser relative to the
    objective, and the visible symptom would be `iterations` dropping rather than anything failing.
  - **`FitError` is the droppable failure and `SchemaError` is not**, exactly as Stage 6 handed it
    over. G6 and G7 raise `FitError` for that reason (§4.5, §6.3).
  - **The failure it counts will come from G7 and not from non-convergence.** Measured (§20): 400
    sparse replicates of this cohort's shape produced zero convergence failures, because separation
    converges (§6.1). A Stage 10 failure counter reading zero is therefore **not** evidence that no
    replicate was degenerate unless G7 is in the path — and Stage 10 should report the G7 count
    separately from the convergence count, because the two mean different things about the data.
  - **`β` is comparable across replicates only under proportional odds**, because `len(fit.alpha)` is a
    property of the replicate (§5.3). A replicate that misses a declared mRS level fits a coefficient on
    a coarser scale, and [§10] takes percentiles across them. §15 files this; it is [§8]'s own
    assumption and not a defect introduced here, but Stage 10 is where it becomes a distribution.
  - **§6.3's band was measured at `K = 6` ONLY, and G7's bound is what decides which replicates are
    dropped.** Every one of the 4800 fits behind the original band was a **seven-category** frame, and
    re-measuring across cutpoint counts moved its lower edge from 8.79 to **11.04** (§6.3, §21.1
    item 24).
    §5.3 makes `K` a property of the replicate, and a replicate that misses mRS 0 or mRS 6 fits five
    cutpoints — on which nothing establishes that the band is still empty, or still that band. This is
    a **nearer** version of §15's "the bound is calibrated on a one-column design and Stage 12's is
    wider": Stage 12 is three stages away and this is the next one. **Stage 10 should report the
    distribution of `len(fit.alpha)` across replicates alongside its G7 count**, because if that
    distribution is degenerate at 6 the gap is theoretical, and if it is not, the drop rate is partly a
    function of a bound calibrated on frames unlike the ones being dropped.
  - **The cost is 1.157 ms per fit** (§2), so the fit alone is about 2.3 s over `N_BOOT`. Stage 8 is
    the bootstrap's cost centre and Stage 6 is not.
- **Stage 11 [§13]** — takes `β` for the E-value. [§13] converts a common odds ratio with `RR ≈ √OR`
  and reports the E-value for the point estimate and for the confidence limit nearest the null, **with
  the rule that the limit's E-value is 1.0 whenever the interval spans the null**. The limit is Stage
  10's, not this stage's; what Stage 11 owes is the statement that the approximation was used, and what
  it inherits from here is that `odds_ratio` is a *common* odds ratio in the overlap population and not
  a risk ratio and not a conditional one.
- **Stage 12 [§14a]** — calls **`model.polr` directly**, with `design(eligible, STANDARDISATION_COVARIATES)`
  and unit weights. It does **not** call `primary`: that function reads a `Propensity` and [§14] fits no
  propensity model anywhere. What it inherits is the estimator, O1-O6, §5.3's collapse and §5.4's
  ordering guarantee. What it owes is its own population, its own [§11] denominator, **its own
  separation guard** — §6's G7 is `outcome.py`'s and Stage 12's design is wider, so a bound calibrated
  on a one-column design is not automatically its — and [§14a]'s guard that `exp(β)` there is a
  **conditional** odds ratio which *"may appear as a model parameter but never as the standardised
  marginal effect"*. That last one is a statement about the same `β` this function returns, which is why
  one implementation serving both stages is what keeps the two labels from being applied to two
  differently-computed numbers.
- **Stage 14 [§16]** — reports `odds_ratio` with its interval from Stage 10, the six `RD_k` **with
  intervals and without p-values** [§8], and the [§11] denominator from `in_estimate`. Three guardrails
  it owes and this stage cannot enforce: the primary p-value is labelled a test of the
  proportional-odds treatment coefficient and never a risk-difference-scale test [§10]; any [§14]
  estimate appearing beside this one carries the different-populations statement [§14a, §16]; and the
  ATO population is described rather than named, which is Stage 6's `overlap_weights` entry and not
  this stage's.

**One rule every one of them owes, and it is Stage 6 §9's, unchanged.** `e` and `w` carry `nan` off
`in_model`; range over the mask, never over `notna()`, and never fill. Stage 8 is the third consumer
and it honours it in one place — `w` is read once, as `ps.w.loc[in_estimate]` in §10 — which is why the
rule appears here as a fact about one line rather than as a warning.

## 12. What Stage 8 amends in Stages 1-7

The full ledger, so that no amendment is discovered during implementation. **Two shipped modules and
two test modules, and one of the shipped changes is a docstring.**

| File | Amendment | Why |
|---|---|---|
| `config.py` | `MRS_LEVELS`, immediately above `MRS_THRESHOLDS` | §8.3. Computed from `PLAUSIBLE_RANGES["mrs_90d"]`, which already declares the range, so the level set is not a second literal. `MRS_THRESHOLDS` is already there and is already commented `# cumulative RD_k [§8]`; it sits in the `[§13]` subgroups block, which is where it was declared and is not where it belongs — it is **not moved**, because moving a constant six test modules import is churn with a failure mode, and the assertion below makes the misfiling harmless |
| `config.py` | `PRIMARY_OUTCOME`, immediately after `DERIVED_DICHOTOMIES` | §0.1. Computed as the single `family == "primary"` key of `OUTCOMES`, so `outcome.py` cannot write `"mrs_90d"` as a literal. `OUTCOMES`' own comment already states that *"exactly one of them is primary: [§8] rests on there being one primary quantity with one test"* — this makes that sentence readable by code |
| `config.py` | the `POLR_*` block, after the `FIRTH_*` one | §5.2, §6.3. Seven constants. They are part of the estimator for `FIRTH_*`'s reason (`config.py:281-284`): [§10] refits in every replicate, so a tolerance and a bound are properties of the sampling distribution and not runtime knobs. The prefix is `POLR_` and not `OUTCOME_` because `model.py` reads them and `model.py` does not know what the outcome is [Stage 6 §0.1] |
| `test_config.py` | **six** assertions | `MRS_LEVELS` is exactly the inclusive `PLAUSIBLE_RANGES["mrs_90d"]` range **and that range's upper bound is not `None`** — the dict's declared type is `tuple[int, int \| None]` and `None + 1` is a `TypeError`, so the computed form in §12 rests on a fact the type does not guarantee; `MRS_LEVELS[:-1] == MRS_THRESHOLDS`, which is what stops the two drifting; `PRIMARY_OUTCOME` is the unique primary-family key and `OUTCOMES` holds exactly one; that outcome's `kind` is `"ordinal"` and its `higher_is_better` is `False`, which is §7.1's precondition asserted where the registry lives; and `POLR_MAX_ABS_BETA` is strictly inside §6.3's measured sparse region, **`11.04 < bound < 18.98`**, asserted against both endpoints as literals with §6.3 cited — so lowering the bound to the 10.0 an earlier draft used **fails here**, which it would not have under the superseded `8.79 < bound < 18.81` (§21.1 item 24). **That last one carries more weight after §14.7's split** (§21.1 item 3): it is now the only assertion in the repository guarding the 4800-fit calibration, and its docstring says so |
| `model.py` | `PolrFit`, `polr`, `_weighted_categories`, `_ord_pieces`, `_ord_loglik`, `_ord_score_hess`, `_assert_polr_fittable`; and **one sentence** in the module docstring | §0.1. The docstring currently says this module *"names neither the treatment nor any covariate list"* and reserves itself for Stages 8, 9 and 12. The added sentence says that an observation-weight vector is opaque in the same way `y` is, so a weighted fitter does not breach the rule — because `propensity.ess`'s docstring calls it a "no-exposure-no-weights rule" and the next reader to grep that phrase will otherwise conclude `polr` is in the wrong file. **No existing function changes** |
| `test_model.py` | the `polr` sections of §14, **and the six generators of §14.0.4** | §14.3, §14.4, §14.7, §14.8. The generators live here because four of the six are consumed here, and `test_outcome.py` imports them by name exactly as it imports `cohort_frame` from `test_cohort` — no new module, so Stage 1's layout rule and `test_config.py:409`'s scan are both untouched. Its existing sections are untouched |
| `test_reference_r.py` | the Stage 8 half of the R oracle, behind the gate Stage 7 built | §18b |
| `implementation_roadmap.md` | Stage 8 gains its `**Spec:**` line, four **Accept when** items and one correction | §21; lands with this document |

**Six things that look like they need amending and do not.**

- **`data.py`.** §9.1: no new kind, no new heading, and therefore none of the literal pins in the six
  existing test modules moves. `KINDS` stays nine and §14.9 asserts it. This is the second stage in a
  row to add audit entries without touching `data.py`, which is what `data.py:151-152` predicted.
- **`propensity.py`.** `Propensity` already carries `e`, `w` and `in_model`, which is everything this
  stage reads. Nothing is added to it. Its `ess` docstring's "no-exposure-no-weights" phrase is
  **not** edited — §0.1 reads it as being about `ess`, and the clarifying sentence goes in `model.py`
  where the rule lives, so that one file states the rule and the other is not made to restate it.
- **`balance.py`.** Not imported (§0.2), not amended, and §14.11's scan asserts `outcome.py` does not
  import it — not because it would fail, but because it would succeed and Stage 7 §11 states that
  Stage 8 does not read the `Balance`.
- **`model.Fit` and `model.firth`.** Untouched. `polr` is a second fitter beside `firth`, not a
  generalisation of it, and the two share only `FitError`. §14.7 asserts `firth`'s existing behaviour
  is unchanged by running the Stage 6 suite unedited.
- **`MRS_THRESHOLDS`.** Declared at Stage 1 for this stage, unread until now, and read rather than
  written — exactly as `SMD_THRESHOLD` was at Stage 7. Its declaration is not moved and its value is
  not touched.
- **`config.RARE_MINORITY_THRESHOLD`.** Declared for [§8]'s rare-outcome rule and it is **Stage 9's**,
  not this stage's: the primary outcome is ordinal and has no minority cell. Checked rather than
  assumed, and recorded here so that nobody wires it into the ordinal path looking for a use.

```python
# config.py — the three additions, as they are to be pasted. §12
#
# The [§5] primary outcome, computed from the registry rather than named a second time. OUTCOMES'
# own comment states that exactly one entry is primary because [§8] rests on there being one primary
# quantity with one test; this is that sentence made readable by code, so that outcome.py cannot
# write the key as a literal and cannot drift from the registry. Declared immediately after
# DERIVED_DICHOTOMIES, which is the first point at which OUTCOMES exists.
PRIMARY_OUTCOME: Final[str] = next(
    key for key, o in OUTCOMES.items() if o.family == "primary")

# The declared mRS level set, computed from the plausible range that already declares it, so the
# level set and the range cannot disagree. MRS_THRESHOLDS below is the K = J - 1 cutpoints and
# test_config.py asserts MRS_LEVELS[:-1] == MRS_THRESHOLDS: RD_k stops at 5 because P(Y <= 6) is 1 in
# both arms by definition [Stage 8 §8.3], and a seventh threshold would be a structural zero.
MRS_LEVELS: Final[tuple[int, ...]] = tuple(
    range(PLAUSIBLE_RANGES["mrs_90d"][0], PLAUSIBLE_RANGES["mrs_90d"][1] + 1))

# The weighted proportional-odds fit [§8, Stage 8 §5]. Prespecified for FIRTH_*'s reason: [§10]
# refits this model in every one of N_BOOT replicates, so a tolerance is a property of the sampling
# distribution and not a runtime knob, and the ORDER the two convergence tests run in is prespecified
# for the same reason [Stage 8 §5.2].
#
# POLR_MAX_ABS_BETA is not a tolerance and is the one constant here that changes an ANSWER: a fit
# reaching it raises, and [§10] drops and counts the replicate. Separation in this estimator does not
# present as non-convergence — measured, it converges in 17 iterations on the score criterion with
# every safeguard counter at zero and returns exp(beta) = 6.5e15 — so this bound is the only thing
# that turns a degenerate fit into a countable failure [Stage 8 §6].
#
# Its value is chosen from a MEASURED EMPTY BAND and not from any estimate. Over 4800 fits of a
# 92-record seven-category frame at twelve true effect sizes from 0 to 6, the largest non-degenerate
# |beta| was 8.7873 and the smallest degenerate one 18.8055, with NOTHING in between: a fit either
# identifies a moderate effect or separates. 14.0 sits near the middle of that band — 59% above the
# first and 25% below the second — so a workbook that pushes either mode needs no re-tuning.
# exp(14) = 1.2e6, five orders of magnitude above any reportable stroke odds ratio, so a rejected fit
# is one whose point estimate no manuscript could print. An earlier draft used 10.0, calibrated on 500
# fits at a single true effect of 0.5, which put it 14% above the largest fit that measurement had
# seen and at the band's lower edge [Stage 8 §6.3, §20].
POLR_MAX_ITER: Final[int] = 200          # measured: 4 to 5 on every frame in Stage 8 §20
POLR_TOL: Final[float] = 1e-8            # |Δ weighted log-likelihood|
POLR_SCORE_TOL: Final[float] = 1e-6      # max |score|; the route a SEPARATED fit returns on
POLR_MAX_HALVINGS: Final[int] = 30       # exhausting them raises; not a convergence route
POLR_MAX_STEP: Final[float] = 5.0        # trust radius, RELATIVE to ‖par‖ [Stage 6 §5.2a]
POLR_ETA_CLIP: Final[float] = 500.0      # keeps exp() in range; measured never approached
POLR_MAX_ABS_BETA: Final[float] = 14.0   # the separation guard [Stage 8 §6.3]
```

## 13. Handover to Stage 9

```
  cohort = cohort.build(eligibility.classify(derive.derive(*data.load()), ...), ...)
  ps     = propensity.fit(cohort, audit)
  bal    = balance.assess(cohort, ps, audit)        # not read by Stage 8
  est    = outcome.primary(cohort, ps, audit)

    93 / 92 / 92    cohort / in_model / in_estimate          [§4.1]
    7 of 7          declared mRS levels occupied, so 6 cutpoints        [§5.3]
    6               RD_k reported; RD_6 is structurally zero            [§8.3]
    3               audit entries, taking the ledger 24 -> 27           [§11]
    0               records the outcome mask removes on v7              [§4.1]
    1.157 ms        one fit; Stage 8 IS on the bootstrap path           [§2]

    beta, exp(beta) and every RD_k on the workbook are DELIBERATELY ABSENT
    from this ledger and from this document. They are in the gitignored log.   [§4.3]
```

Stage 9 receives the cohort and the `Propensity`, exactly as Stages 7 and 8 did. **It does not read
the `Primary`**, and it does not read the `Balance`: each estimate in [§8] stands on its own [§11]
denominator, and the ordinal quantity is an input to nothing.

What Stage 9 inherits and should not rebuild: `weighted_proportion` (§8.1), which is its risk
difference and its marginal odds ratio; `model.design` with `outcome_model_covariates(outcome)` for
`m_a(X)`, which is `config.py:441-455` and already written; `model.firth` for that fit, which is
Stage 6's; and `RARE_MINORITY_THRESHOLD`, which is Stage 1's and which §12 records as Stage 9's and
not this stage's.

## 14. Acceptance criteria

Tests live in `test_outcome.py`, with section banners matching these numbers, as `test_cohort.py`,
`test_model.py`, `test_propensity.py` and `test_balance.py` do. The `polr` sections — §14.3, §14.4,
§14.7, §14.8 — live in `test_model.py`, because that is where the function is. Tests needing the
private workbook reuse the `DATA_GATED` idiom and §14.0.3's fixture; everything else runs on a plain
checkout with no `data/`.

### 14.0 The frames and fixtures this stage is tested on

**And the first of them is a problem no earlier stage had.** Four kinds of input, and the split
follows Stage 6 §12.0's and Stage 7 §12.0's with one addition that is forced rather than chosen.

**1. Hand-built ordinal vectors, for `polr` and `cumulative_rd`.** Plain `numpy` arrays with
hand-computable answers, in `test_model.py` and `test_outcome.py` respectively. They carry the
roadmap's criteria and every branch of §5.5, and they need no frame.

**2. `ordinal_cohort()` — the committed fixture with one column varied**, because the committed
fixture **cannot reach this stage at all**:

```
    test_cohort.cohort_frame()  →  build  →  5 records, 3 bridging / 2 control
    mrs_90d on all five:  [2, 2, 2, 2, 2]
    polr  →  FitError: "O5  1 response category/ies carry positive weight: (2.0,)"
```

Measured (§20). `cohort_frame()` builds its three cohort records from the `HAND-4` template
(`tests/test_cohort.py:124-127`), so every analysis column not overridden is constant across the
frame — and `mrs_90d` is one of them. Stage 6 fitted a propensity model on it because the *exposure*
varies there by construction; Stage 7 computed nineteen balance rows on it, fourteen of them exactly
0.0, and said so. **Stage 8 is the first stage whose response is a column the fixture holds constant**,
and the constant is 2.

`hand_frame()` and `tests/fixture_schema.xlsx` cannot reach this stage either — Stage 5's P2 raises
first, measured at Stages 6 and 7 and unchanged.

So this stage declares one frame, and it declares it as an **override on the committed fixture rather
than as a new frame**, for Stage 3 §12's reason: a second hand-built cohort would be a second place
the pipeline's fixture data lives, and the two would diverge on the next Stage 1 amendment.

```python
# test_outcome.py — §14.0
def ordinal_cohort() -> pd.DataFrame:
    """`cohort_frame()` with mrs_90d spread across the declared levels. §14.0.

    The committed fixture holds mrs_90d constant at 2 on all nine records, so its cohort reaches
    O5 and not a fit (§14.0.1, measured). This is the ONE column Stage 8 has to vary, and it is
    varied by `cohort_frame`'s own `**overrides` rather than by a new frame, so Stage 5's
    restrictions, Stage 3's derivations and Stage 2's schema all still run over exactly the records
    they run over everywhere else.

    The spread is a `pd.Series` over the frame's index and covers 0 through 6 in order, so which of
    the nine records survive both restrictions decides which levels the cohort carries — measured
    (§20): the five survivors take mrs_90d 0, 1, 4, 6 and 3, so the cohort is missing levels 2 and 5
    and §5.3's collapse fires on it. That is the point: the collapse branch has a real-frame witness
    that runs with no `data/`, and the workbook has all seven levels and cannot provide one.

    The derived dichotomies come back recomputed from the ordinal source rather than carried, which
    §14.0.1 asserts, because that is Stage 3's contract and this override is the first thing in the
    repository to exercise it on a varying mrs_90d.

    The `[: len(base)]` slice is a ONE-DIRECTIONAL guard and stays one, which is a decision rather than
    an oversight: it is total if `cohort_frame()` ever shrinks and raises a bare pandas length mismatch
    if it GROWS, and a tenth fixture record is the likelier of the two edits. Cycling MRS_LEVELS over
    the index instead would be total in both directions — and would give a DIFFERENT spread in the
    surviving positions, so it would move §14.0.2's golden vector. The pin wins; §14.0.1's assertion
    that the survivors are the five case_ids test_cohort.py already pins is what names the cause when
    the mismatch fires.
    """
    base = cohort_frame()
    spread = pd.Series([0, 1, 2, 3, 4, 5, 6, 2, 3][: len(base)],
                       index=base.index, dtype="Int64")
    return cohort_frame(mrs_90d=spread)
```

**The alternative was measured before it was declined.** Cycling `MRS_LEVELS` over nine records gives
`[0, 1, 2, 3, 4, 5, 6, 0, 1]` against this literal's `[0, 1, 2, 3, 4, 5, 6, 2, 3]`, and the survivors
sit at index positions 0, 1, 4, 6 and 8 — so the cycled form's cohort would take `mrs_90d`
`[0, 1, 4, 6, 1]` and fit `categories` **`(0, 1, 4, 6)`**, four levels and three cutpoints, where this
one fits `(0, 1, 3, 4, 6)` with four cutpoints. Every number in §14.0.2 would move, and this fixture's
whole purpose is to be the frame those numbers are pinned on. A totality repair that changes the pin is
not a repair. The hazard is therefore recorded rather than removed, and §21.1 lists it as the one review
item raised and **not** applied.

`test_outcome.py` imports the fixture machinery from the modules that declare it and never
re-declares it (Stage 3 §12's rule). The modules, named because there are two:

```python
from test_cohort import built, cohort_frame, set_cell   # tests/test_cohort.py:116, :147, :152
from test_data import hand_source, run                  # tests/test_data.py:98, :107
```

**3. A synthetic separated construction**, for §6. Hand-built, no frame, no patient data: 20 control
records at one mRS level and 20 treated at another, which is the input §6.1 measured.

**4. The workbook, data-gated.** §14.12 only, through §14.0.3's fixture, and **asserting properties
rather than values** (§4.3).

### 14.0.1 What `ordinal_cohort()` exercises, measured

Measured (§20), on `ordinal_cohort()` normalised and taken through `derive → classify → build → fit`:

```
    5 records, 3 bridging / 2 control, in_model 5 of 5, in_estimate 5 of 5
    case_ids       HAND-1, HAND-2, HAND-5, COHORT-1, COHORT-3
    mrs_90d        0, 1, 4, 6, 3
    categories     (0, 1, 3, 4, 6)   ->  5 of 7 declared levels, 4 cutpoints
    absent         2 and 5           ->  §5.3's collapse fires, and §9.2 renders both rows
```

and the derived dichotomies come back recomputed from the ordinal source, which is what
`mrs_90d = [0, 1, 4, 6, 3]` requires and what a carried column would not give:

```
    mrs_0_2_90d    1, 1, 0, 0, 0        mrs_0_1_90d    1, 1, 0, 0, 0
    death_90d      0, 0, 0, 1, 0        mrs_5_6_90d    0, 0, 0, 1, 0
```

- **The bare-versus-normalised rule still applies.** Stage 5 §12.0.1's distinction carries forward:
  `cohort_frame()` returns raw centre codes and `run(...)` maps them. A test that means to reach
  `primary` and passes a bare frame never gets there — Stage 5's C4 raises first — and the failure
  looks like a passing `pytest.raises` whose message is never read.
- **`built(ordinal_cohort())` is asserted to raise nothing**, and the five-record cohort is asserted
  to be the same five `case_id`s `test_cohort.py` already pins. If the override changed which records
  survive, every number in §14.0.2 would be describing a different population.

### 14.0.2 The golden vector

Pinned in `test_outcome.py` as literals, from `ordinal_cohort()` so that **no patient-derived number
enters git** (Stage 6 §7.3) — and, at this stage, so that no *estimate* enters it either (§4.3).
Measured (§20):

```
  the fit
    categories        (0, 1, 3, 4, 6)          4 cutpoints over 5 declared-and-occupied levels
    beta              1.7516026138
    exp(beta)         5.7638324754
    alpha            [-4.0128444951, -1.8393869358, -0.7580715177, 0.4483367817]
    cut_levels        (0, 1, 3, 4)
    iterations        4       converged_on  "likelihood"     rescales 0    halvings 0

  the cumulative risk differences
    k      P(Y <= k), EVT alone     P(Y <= k), bridging          RD_k
    0            0.0000000000            0.1079772317      +0.1079772317
    1            0.0000000000            0.5469186443      +0.5469186443
    2            0.0000000000            0.5469186443      +0.5469186443
    3            0.5287822719            0.5469186443      +0.0181363724
    4            0.5287822719            1.0000000000      +0.4712177281
    5            0.5287822719            1.0000000000      +0.4712177281
```

Asserted to 1e-6, for Stage 6 §12.0.2's and Stage 7 §12.0.2's reason: this is a regression pin against
an edit to our own arithmetic, and pinning it at machine precision would make it fail on a numpy patch
release rather than on a mistake. **The assertion is by name and by threshold, one at a time** — `α`
against `cut_levels` and `RD_k` against `k` — because the α vector is monotone and a set comparison
would pass while permuting it.

**Three things in that block are load-bearing beyond being a pin**, and §14.6 and §14.4 assert each
separately so that a failure names its own cause:

- **`k = 1` and `k = 2` are identical, and so are `k = 4` and `k = 5`.** Levels 2 and 5 are unoccupied,
  so the cumulative probability cannot move across them. The assertion is **non-decreasing**, never
  strictly increasing (§8.3).
- **`cut_levels` is `(0, 1, 3, 4)` and not `(0, 1, 2, 3)`.** The third cutpoint sits above mRS 3
  because level 2 is absent; labelling it `MRS_THRESHOLDS[2]` would name a cutpoint with a threshold
  it is not (§10).
- **`P(Y <= 3)` in the control arm is 0.5287822719 while `P(Y <= 2)` is 0.0**, which is one control
  record's whole weight arriving at once on a two-record arm. A five-record cohort makes every
  cumulative probability a small sum of weights, which is what makes this vector reproducible by hand
  and what makes it a weak witness for anything about the workbook (§4.3, §15).

### 14.0.3 The workbook fixture returns three things from ONE audit

Module-scoped and this file's own, as Stages 5, 6 and 7's are — never shared across modules. Its shape
is specified rather than left to the pattern, because §14.10 reads `overlap_weights` out of the audit
and compares it with Stage 8's own entries, so both must be in the same `Audit`:

```python
@pytest.fixture(scope="module")
def workbook():
    """(df, ps, audit) from ONE linear run against ONE Audit. §14.0.3."""
    df, audit = data.load(data.WORKBOOK)
    df = derive.derive(df, audit)
    df = eligibility.classify(df, audit)
    df = cohort.build(df, audit)
    return df, propensity.fit(df, audit), audit
```

`test_propensity.py:593-596`'s `workbook_ps` deliberately does the opposite, and Stage 7 §12.0.3
records why. The audit this fixture returns carries **22** entries before `primary` and **25** after —
`balance.assess` is not called by it, so the ledger reaching 27 is the full-pipeline count of §11 and
not this fixture's. §14.9's captured-index assertions are counted against 22, never against a tail
slice.

### 14.0.4 The synthetic constructions, written out

**A document that pins a number produced by a construction it does not specify is not a sole source,
and this section is the repair.** Six constructions carry numbers that §7.2, §14.3, §14.5, §14.7,
§14.6 and §18b assert or quote, and an earlier draft described all six in prose only — "a 200-record
construction with the treated arm's mass at mRS 0-3 and the control's at 2-6", "a 120-record
five-category frame with weights drawn from `{1, 2, 3, 4}`". An implementer cannot reach
`exp(β) = 84.584014` from the first sentence or `1.332e-15` from the second. They are therefore
written out here as **deterministic generators**, in `test_model.py` beside the sections that consume
them, and `test_outcome.py` imports them by name exactly as it imports `cohort_frame` from
`test_cohort` (§14.0).

**Every one seeds from `C.SEED` and none holds a literal seed**, for `FIRTH_TOL`'s reason one level
down: a seed is what makes a measured number reproducible, so a second seed in a test module is a
second answer to a question `config.py` already answers. Stage 10 reads the same constant.

```python
# test_model.py — §14.0.4. Six generators, no assertions, no fixtures, no frames.

def hand_ordinal() -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """12 records, 6 against 6, all seven declared mRS levels. (y, a, w) — §14.3, §14.6.

    Hand-computable and hand-checkable: the weights are small integers summing to 21 in each arm, so
    every weighted cumulative probability in §14.6 is a sum of at most six of them over a denominator
    of 21. This is the frame §14.6 checks RD_k against a hand computation on, and three earlier drafts
    of this document called it "§14.0.2's twelve-record frame" — which it is not: §14.0.2's frame is
    `ordinal_cohort()`'s FIVE-record cohort. §21.1 item 2.

    **THE WEIGHTS ARE INTEGERS, AND THAT IS A REQUIREMENT AND NOT A CONVENIENCE.** §14.3 and §18b's
    oracle 1 run the integer-weight replication check on this frame, and "the unweighted fit on the
    frame with each row repeated `w` times" is undefined for a fractional weight — `np.repeat` raises
    `TypeError: Cannot cast array data from dtype('float64') to dtype('int64')`. An earlier draft of
    this generator used two-decimal weights and could not carry the oracle §14.3 assigns it (§21.1
    item 12).

    **And scaling fractional weights up to integers does NOT work**, which is measured rather than
    assumed: the polr MLE is invariant to a common weight scale in exact arithmetic, but `POLR_TOL`
    is ABSOLUTE on a log-likelihood that scales with the weight sum, so multiplying every weight by
    20 tightens the stopping criterion twentyfold and the two fits stop at different points —
    `|Δbeta| = 9.458e-10` measured, which FAILS §14.3's 1e-10 tolerance (§20c). Integers chosen from
    the start agree at **1.11e-16**. That is the first measured consequence of §5.2's absolute
    tolerances, and §14.3 now states the tolerance per weight scale because of it.
    """
    y = np.array([0., 1., 2., 3., 4., 5., 1., 2., 3., 4., 5., 6.])
    a = np.array([1.] * 6 + [0.] * 6)
    w = np.array([1., 2., 3., 4., 5., 6., 6., 5., 4., 3., 2., 1.])
    return y, a, w


def replication_frame(n: int = 120, levels: int = 5) -> tuple[pd.DataFrame, np.ndarray, np.ndarray]:
    """(X, y, w) with INTEGER weights in {1, 2, 3, 4}, for the replication oracle. §14.3, §18b.

    Integer by construction and not by rounding: the oracle is that a weighted fit equals the
    unweighted fit on the row-replicated frame, and a non-integer weight makes "replicated" undefined.
    Two covariates so the oracle covers a beta of length > 1, which the [§8] path never has and
    Stage 12 always does.

    **THE LOGISTIC NOISE TERM IS LOAD-BEARING AND WAS MISSING FROM THE FIRST DRAFT OF THIS
    GENERATOR.** This is the proportional-odds data-generating process: a latent variable
    `x'b + Logistic(0, 1)`, cut at its own quantiles. Without the noise, `y` is a MONOTONE function of
    `x'b` — every category is an interval of the linear predictor — and the frame is **perfectly
    separated in the ordinal sense**, so the likelihood has no interior maximum and `beta` runs to the
    trust region's limit. Measured on the noiseless draft (§20c): `beta = [-3281.7, -1652.8]` at n=120
    and `[-48828, -24382]` at n=500, both converging on the score criterion with clean counters, which
    is §6.1's finding arriving in a place nobody was looking for it. **Every oracle built on this
    frame is void on a separated one**: two optimisers of a likelihood whose maximum is at infinity
    need not agree at all, so §14.3's 1e-10 replication tolerance and §14.5's statsmodels comparison
    would both have been asserting agreement between two arbitrary stopping points. With the noise,
    `max|beta|` is 0.93 here and 0.76 on `reference_frame()` — comfortably interior (§20c).
    """
    rng = np.random.default_rng(C.SEED)
    X = pd.DataFrame({"x1": rng.normal(size=n), "x2": rng.normal(size=n)})
    latent = (0.8 * X["x1"].to_numpy() + 0.4 * X["x2"].to_numpy()
              + rng.logistic(size=n))                          # <- the noise. §20c.
    cuts = np.quantile(latent, np.linspace(0.0, 1.0, levels + 1)[1:-1])
    y = np.searchsorted(cuts, latent).astype(float)
    w = rng.integers(1, 5, size=n).astype(float)
    return X, y, w


def reference_frame(n: int = 500, levels: int = 5) -> tuple[pd.DataFrame, np.ndarray]:
    """(X, y), two covariates, unweighted — the frame the statsmodels sign oracle runs on. §7.2, §14.5.

    Unweighted deliberately: `OrderedModel` has no weight support (§2), so the only comparison
    available is the unweighted one, and this generator is what makes `ours.beta + sm.beta` a
    reproducible quantity rather than a remembered one.

    It delegates to `replication_frame` rather than repeating the construction, so it inherits the
    noise term that docstring is about — and it inherited the BUG for the same reason, at n = 500
    where separation is sharper. §14.5's sign oracle is the assertion that would have been meaningless
    (§20c).
    """
    X, y, _ = replication_frame(n=n, levels=levels)
    return X, y


def orientation_frame(per_arm: int = 100) -> tuple[pd.DataFrame, np.ndarray, np.ndarray]:
    """(X, y, w) where treatment shifts mRS DOWNWARD by construction. §7.3, §14.5.

    The construction whose truth is known, which is what the roadmap's orientation criterion needs
    and what prose cannot supply: the treated arm's mass sits at mRS 0-3 and the control's at 2-6, so
    the direction of the true effect is a property of these two literal tuples and not of a fitted
    number. Unit weights — the orientation of `beta` is a question about the parametrisation (§5.1)
    and not about the weighting, and §14.3's oracles own the weighting.
    """
    rng = np.random.default_rng(C.SEED)
    treated = rng.choice(np.array([0., 1., 2., 3.]), size=per_arm, p=[.40, .30, .20, .10])
    control = rng.choice(np.array([2., 3., 4., 5., 6.]), size=per_arm, p=[.20, .25, .25, .20, .10])
    y = np.concatenate([treated, control])
    a = np.concatenate([np.ones(per_arm), np.zeros(per_arm)])
    return pd.DataFrame({C.TREATMENT: a}), y, np.ones(len(y))


def separated_frame(per_arm: int = 20, crossovers: int = 0) -> tuple[pd.DataFrame, np.ndarray,
                                                                    np.ndarray]:
    """PERFECT separation at crossovers = 0, NEAR separation above it. §6.1, §14.7.

    Every treated record at mRS 0 and every control at mRS 5, with `crossovers` treated records moved
    to the control's level. This is the construction §6.1 measured — it holds no seed at all, which is
    why §6.1's numbers are exact rather than distributional: at crossovers = 0 and per_arm = 20 it
    returns beta 36.4058 in 17 iterations on the score criterion, every counter at zero, and it does
    so on every machine (§20).
    """
    y = np.concatenate([np.zeros(per_arm), np.full(per_arm, 5.0)])
    y[:crossovers] = 5.0
    a = np.concatenate([np.ones(per_arm), np.zeros(per_arm)])
    return pd.DataFrame({C.TREATMENT: a}), y, np.ones(len(y))


def band_samples(effects: tuple[float, ...], per_effect: int, per_arm: int = 46,
                 levels: int = 7) -> list[PolrFit]:
    """The FITS from `len(effects) * per_effect` draws of a cohort-shaped frame. §6.3, §14.7.

    It returns the fits and not a `|beta|` vector, which is a correction: §14.7 asserts things about
    `cond(-H)` and about the safeguard counters, and both are properties of the fit. An earlier draft
    returned magnitudes only and left those assertions with no subject (§21.1 item 21). The band
    assertion takes `max(abs(f.beta))` from each; nothing is lost by returning more.

    A treatment-only weighted proportional-odds sampler over a range of TRUE effect sizes, returning
    the |beta| vector and nothing else. The range is the point: §6.3's band was first measured at a
    single true effect of 0.5, which is an assumption about the answer and is what §4.3 forbids, so
    the generator takes the effects as an argument and both callers pass a range.

    TWO callers with different budgets, which is §14.7's split. The in-suite probe passes
    per_effect = 20 for 240 fits and about 0.3 s; §6.3's calibration passed per_effect = 400 for 4800
    fits and about 6 s, which is a calibration and not a test (§21.1 item 3).
    """
    rng = np.random.default_rng(C.SEED)
    out: list[PolrFit] = []
    for effect in effects:
        for _ in range(per_effect):
            a = np.concatenate([np.ones(per_arm), np.zeros(per_arm)])
            latent = effect * a + rng.normal(size=2 * per_arm)
            cuts = np.quantile(latent, np.linspace(0.0, 1.0, levels + 1)[1:-1])
            y = np.searchsorted(cuts, latent).astype(float)
            w = rng.uniform(0.05, 0.95, size=2 * per_arm)
            try:
                out.append(model.polr(pd.DataFrame({C.TREATMENT: a}), y, w))
            except model.FitError:
                continue                       # O5 on a degenerate draw; NOT a band member
    return out
```

**`BAND_EFFECTS` is a module constant in `test_model.py` and not a literal in two calls**, for the same
reason the generators seed from `C.SEED`:

```python
BAND_EFFECTS: Final[tuple[float, ...]] = (0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 6.0)
```

**And `band_samples` returns the FITS, not a `|beta|` vector.** An earlier draft returned
`np.asarray(out)` of magnitudes only, which left §14.7's `cond(-H)` assertions with no subject —
the Hessian is a property of the fit and the fits had been discarded (§21.1 item 21). It returns
`list[PolrFit]`; §14.7's band assertion takes `max(abs(f.beta))` from each, and its `cond(-H)`
assertion recomputes the Hessian from the fit it names.

### 14.0.4a The six that O1-O6, §14.4 and §14.6 need

**§14.0.4 repaired six constructions and there were more**, which an earlier DoD-18 asserted there were
not (§21.1 item 22). These are the rest: every one is hand-built, seed-free, and small, so their
numbers are **values** rather than tolerances.

```python
# test_model.py — §14.0.4a. Hand-built, no seed, so §14.8's and §14.4's numbers are exact.

def nan_weight_frame() -> tuple[pd.DataFrame, np.ndarray, np.ndarray]:
    """A four-category weighted frame, for O3's ordering. §14.8.

    Twenty records over four categories with strictly positive weights. The test sets `w[0]` to nan
    and asserts O3; the companion removes O3 and asserts the category DELETION (§5.5), which is the
    mechanism no other check in this stage names.
    """
    y = np.array([0., 1., 2., 3.] * 5)
    a = np.array([1., 0.] * 10)
    return pd.DataFrame({C.TREATMENT: a}), y, np.linspace(0.2, 0.9, 20)


def nan_response_frame(n: int = 80) -> tuple[pd.DataFrame, np.ndarray, np.ndarray]:
    """Eighty records over four categories, for O2 — and for O2's interaction with O4. §14.8."""
    rng = np.random.default_rng(C.SEED)
    a = np.array([1., 0.] * (n // 2))
    latent = 0.9 * a + rng.logistic(size=n)
    cuts = np.quantile(latent, [0.25, 0.50, 0.75])
    return (pd.DataFrame({C.TREATMENT: a}),
            np.searchsorted(cuts, latent).astype(float), np.ones(n))


def noninteger_response_frame() -> tuple[pd.DataFrame, np.ndarray, np.ndarray]:
    """THREE distinct half-integer values over sixty records, for O4. §14.8.

    Sized so that the with-O4-removed companion reaches a FIT. Six distinct values over 24 records
    makes the Hessian singular and `np.linalg.solve` raises `LinAlgError` — which is not `FitError`,
    so the companion would assert the wrong thing (§20c, measured).
    """
    y = np.array([0., 0.5, 1.] * 20)
    a = np.array([1., 0.] * 30)
    return pd.DataFrame({C.TREATMENT: a}), y, np.ones(60)


def emptied_category_frame(empty: int = 2) -> tuple[tuple, tuple]:
    """(full, reduced) — a five-category frame and the same frame with one category removed. §14.4.

    Returns BOTH, because §14.4's assertion is that the collapse gives the same fit as the removal.
    The zero-weight variant is the same frame with `w` zeroed on that category rather than the rows
    dropped, which is the one input separating §5.3's rule from the observed-set rule.
    """
    y = np.array([v for v in (0., 1., 2., 3., 4.) for _ in range(6)])
    a = np.array([1., 0.] * 15)
    w = np.linspace(0.3, 0.8, 30)
    keep = y != float(empty)
    return ((pd.DataFrame({C.TREATMENT: a}), y, w),
            (pd.DataFrame({C.TREATMENT: a[keep]}), y[keep], w[keep]))


def rare_category_frame() -> tuple[pd.DataFrame, np.ndarray, np.ndarray]:
    """Sixty-one records, ONE of them in the top category. §14.7.

    The frame on which `_assert_reportable` is shown to bound `beta` and not `alpha`: `max|alpha|` is
    4.1599 and `|beta|` 0.126998, so with the bound patched to 3.0 an alpha-bound would reject it and
    the real guard passes. No frame with `|alpha| >= 14` is constructible — alpha for a rare category
    grows like log n and 14 needs order 1e7 records (§20c).
    """
    y = np.array([0.] * 30 + [1.] * 30 + [2.])
    a = np.array(([1., 0.] * 15) * 2 + [1.])
    return pd.DataFrame({C.TREATMENT: a}), y, np.ones(61)


def deaths_both_arms_frame() -> tuple[pd.DataFrame, np.ndarray, np.ndarray, np.ndarray]:
    """Fourteen records, all seven levels in both arms, a death in each. §14.6, §8.3.

    Returns the arm vector as well, because §14.6 calls `cumulative_rd` on it directly. It is the
    frame on which RD_6 is shown structurally zero while RD_5 is not.
    """
    y = np.array([0., 1., 2., 3., 4., 5., 6.] * 2)
    a = np.array([1.] * 7 + [0.] * 7)
    w = np.array([3., 3., 3., 2., 2., 1., 1., 1., 2., 2., 3., 3., 3., 4.])
    return pd.DataFrame({C.TREATMENT: a}), y, a, w
```

Measured on them (§20c), and these are the numbers §14.8, §14.4 and §14.6 assert:

```
  nan_weight_frame       clean: 4 categories, beta +1.678900
                         w[0] = nan  ->  O3 raises
                         with O3 REMOVED: 4 categories -> 3, beta +1.678900 -> +0.194854,
                                          iterations 3, converged_on "likelihood", nothing raised
  nan_response_frame     with O2 AND O4 removed: 79 of 80 fitted, categories (0,1,2,3) UNCHANGED,
                                          beta -1.162237 -> -1.199291
                         with O2 alone removed: O4 raises  (§21.1 item 18)
  noninteger_response    O4 raises; with O4 removed, 3 categories and 2 cutpoints — one per
                                          distinct value
  emptied_category_frame full 4 cutpoints (0,1,2,3,4); emptied 3 cutpoints (0,1,3,4);
                         alpha ascending in both
                         zero-weight variant == rows-removed variant: categories identical,
                                          |dbeta| 0.0, |dalpha| 0.0
  rare_category_frame    alpha [0.0307, 4.1599], max|alpha| 4.1599, beta -0.126998, inside the bound
  deaths_both_arms_frame RD_5 +0.1555555556, RD_6 +0.0000000000 exactly,
                         P(Y <= 6) = 1.000000 in both arms, one death per arm
```

**`separated_cohort()` is `test_outcome.py`'s and not `test_model.py`'s**, because §14.7's "`primary`
on the same frame raises G7" needs a cohort frame and a `Propensity`, which `separated_frame()`'s three
arrays are not:

```python
# test_outcome.py — §14.0.4a
def separated_cohort() -> tuple[pd.DataFrame, propensity.Propensity, Audit]:
    """`ordinal_cohort()`'s cohort with the outcome made to separate the arms perfectly. §14.7.

    The Propensity is the REAL one, fitted before the outcome is overwritten, so the weights are a
    genuine [§7] score and only the response is constructed — which is what makes `primary` raise G7
    for the reason §6 gives rather than because the frame is unreachable.
    """
    df, audit = built(ordinal_cohort())
    ps = propensity.fit(df, audit)
    df = df.copy()
    df.loc[df[C.TREATMENT] == 1, C.PRIMARY_OUTCOME] = C.MRS_LEVELS[0]
    df.loc[df[C.TREATMENT] == 0, C.PRIMARY_OUTCOME] = C.MRS_LEVELS[-1]
    return df, ps, audit
```

**What this section does and does not settle.** It settles that every construction is reproducible and
that §14's assertions have a subject. It does **not** re-derive §20's measured numbers against these
generators: §20's values came from the probe harnesses that produced them, whose sampling differs in
detail from the code above, so **the numbers in §20 are the record of what was seen and not a
prediction of what these generators return.** Every §14 assertion that depends on a generator is
therefore written as a **tolerance or a property** rather than as a value — `assert abs(Δβ) < 1e-10`,
not `assert β == 0.0` — which is the form §14.3's replication oracle already had and which §21.1 item 1
extends to the other five. The two constructions that are seed-free, `hand_ordinal()` and
`separated_frame()`, are the exceptions: their numbers are exact and §14.6 and §14.7 pin them.

### 14.0.5 `truncated_cohort()` — the frame that witnesses §9.2's status branch

One more override on the committed fixture, and it exists for one branch:

```python
def truncated_cohort() -> pd.DataFrame:
    """`ordinal_cohort()` with the worst mRS removed, so the top FITTED level is not the top DECLARED
    one. §9.2, §14.9.

    `ordinal_cohort()`'s survivors take mrs_90d 0, 1, 4, 6 and 3, so its top fitted category is 6 —
    which is `MRS_LEVELS[-1]`, by coincidence and not by design. The workbook occupies all seven
    levels, so its top fitted category is 6 as well. **Every frame in this suite therefore has
    `fit.categories[-1] == MRS_LEVELS[-1]`, and §9.2's status branch is untested without this one.**

    Replacing the 6 with a 2 gives categories (0, 1, 2, 3, 4): five occupied levels, four cutpoints,
    mRS 5 and 6 both absent, and the top fitted level is 4. That is the frame on which the draft
    §9.2 corrects rendered `alpha at mRS <= 4` as `missing` under status `fitted`, and `alpha at
    mRS <= 6` as the highest declared level with no mention that nobody is in it.

    A [§10] replicate that draws no death produces exactly this shape, which is why the branch is
    worth a fixture rather than a note.
    """
    base = cohort_frame()
    spread = pd.Series([0, 1, 2, 3, 4, 5, 2, 2, 3][: len(base)],
                       index=base.index, dtype="Int64")
    return cohort_frame(mrs_90d=spread)
```

The value the **seventh** element takes is what does the work, and which element to change is derived
from a measurement rather than guessed: §14.0.1 records that `ordinal_cohort()`'s survivors carry
`mrs_90d` `[0, 1, 4, 6, 3]` against the spread `[0, 1, 2, 3, 4, 5, 6, 2, 3]`, so the surviving index
positions are 0, 1, 4, 6 and 8 — the seventh element is one of them and the sixth is not. Which
records survive does not depend on `mrs_90d` at all: the [§3] restrictions are on centre, arm and
eligibility, so the override moves the level set and never the population, which is §14.0's own
argument for overriding rather than building a second frame.

**§14.9 asserts `fit.categories == (0, 1, 2, 3, 4)` on this frame before it asserts anything about the
table.** Without that, a change to which records survive would silently turn this fixture back into one
whose top fitted level is the top declared one, and the branch would go back to being untested with
every test still green — which is the failure §14.2 exists to prevent one section earlier.

### 14.1 The population and the collapse

- `in_estimate` is asserted to be `in_model & mrs_90d.notna()`, boolean, total, and on the frame's
  index — and asserted to equal neither `in_model` nor `notna()` on a frame where they differ, which
  §14.2 constructs.
- **On the workbook all seven declared mRS levels are occupied and `len(fit.alpha)` is 6**, asserted
  **with the reason** — so the collapse branch is recorded as inert on this workbook rather than
  assumed to be. Data-gated.
- **On `ordinal_cohort()` the collapse fires**: `fit.categories` is `(0, 1, 3, 4, 6)`, `len(fit.alpha)`
  is 4, and `MRS_LEVELS` minus `categories` is `(2, 5)`. This is the real-frame witness the workbook
  cannot provide, and it runs with no `data/`.
- A frame whose response is constant is asserted to raise O5, parametrised over the constant, with the
  message asserted to name the category count — the branch `cohort_frame()`'s own cohort reaches
  (§14.0).

### 14.2 The [§11] denominator is per estimate, and no available frame varies it

The column §4.1 argues for, tested on a constructed frame because **nothing available varies it**:
every workbook record and every `ordinal_cohort()` record has its outcome present (§20).

- On `ordinal_cohort()` with **one `in_model` record's `mrs_90d` blanked** via `set_cell`:
  `in_estimate` is 4 against `in_model`'s 5, entry 1's `n` is 1, and the entry **names that record**.
- **`propensity.fit`'s `in_model` is asserted UNCHANGED** by the same edit. This is half the test:
  without it, an implementation that complete-cased the propensity model on the outcome too would pass
  the first assertion by removing the record from both.
- The fit is asserted to be computed over the four, by comparing `Primary.fit` against `model.polr`
  called directly on the masked triple — so the estimate is shown to be over its own records and its
  own weights rather than over the table's.
- **The contrast**: blanking a `PS_COVARIATE` instead moves `in_model` **and** `in_estimate` together.
  That is what makes the first three assertions mean something rather than merely pass.
- And the second half of §4.1's rule: a frame with a record whose outcome is present but whose
  covariates are not is asserted to have that record in **neither** mask, so `in_estimate` is the
  intersection and not the outcome mask.

Without this section the mask is a copy of `in_model` on every frame in the suite, and an
implementation that never built it would be green.

### 14.3 `polr` and the weights [roadmap, amended]

In `test_model.py`.

**Every assertion in this section is a tolerance or a property and none is a value**, which was already
true of the first bullet and is now true of all of them (§14.0.4, §21.1 item 1). The measured numbers
are kept beside them as the record of what was seen on the harnesses of §20, and an implementer who
gets a different one from the generators of §14.0.4 has found something rather than broken something.

- **Integer weights equal row replication.** A fit with integer weights `w` equals an unweighted fit on
  the frame with each row repeated `w` times, in `beta` and in every `alpha`. Asserted to **1e-10** on
  `replication_frame()` and on `hand_ordinal()`. **The tolerance is stated per weight scale, and that
  is measured rather than stylistic**: `POLR_TOL` is absolute on a log-likelihood that scales with
  `Σw` (§5.2), so how closely two fits of the same likelihood agree depends on the weight sum. At
  `Σw = 307` (`replication_frame()`) they agree at **0.0**; at `Σw = 42` (`hand_ordinal()`) at
  **1.11e-16**; and at `Σw = 138` — `hand_ordinal()`'s weights scaled by 20, which is mathematically
  the same MLE — at **9.458e-10**, which would FAIL a 1e-10 assertion (§20c, §21.1 item 12). So 1e-10
  holds on both declared fixtures and is **not** a scale-free guarantee, and a future fixture at a very
  different `Σw` must re-measure rather than inherit it. Measured (§20) at **0.000e+00** on `beta` and
  **1.332e-15** on `alpha` over a 120-record five-category frame, and at **1.11e-16** on `beta` over
  `hand_ordinal()` — which is `hand_ordinal()` and **not** §14.0.2's, whose cohort has
  five records; three earlier drafts of this document called it §14.0.2's and §21.1 item 2 records the
  correction. **This is the roadmap's weighting criterion in the form that can fail**, and §21's
  amendment 2 is why.
- **The companion that makes it necessary.** `polr(X, y, np.ones(n))` is asserted **exactly equal** to
  `polr(X, y, None)` — measured, bit for bit — and the test's comment states the consequence: an
  implementation that never reads `w` at all satisfies the roadmap's "equals an unweighted fit when all
  weights are 1" criterion. The pair is the test. A third assertion pins that the weighted fit
  **differs** from the unweighted one on the same data, **by more than 1e-3 in `beta`** on
  `hand_ordinal()`, so "reads `w`" is asserted positively rather than only negatively. Asserted as a
  floor and not as a value: measured **1.8268** on `hand_ordinal()` — where the weighted and unweighted
  fits have **opposite signs**, −0.8696 against +0.9573 — and 0.0566 on `replication_frame()`. The
  difference is a property of each frame's weights, and a magnitude assertion here would pin the
  fixture rather than the behaviour. An earlier draft quoted 0.2781, which is the §20 harness's frame
  and not either of these (§21.1 item 12).
- The analytic score and Hessian are asserted against **central differences** of `_ord_loglik` and of
  `_ord_score_hess`, to 1e-6 — measured 2.784e-08 and 2.812e-08. The Hessian is asserted **symmetric to 1e-12 and negative definite**, which is §5.2's concavity claim as
  a measurement rather than a citation. **A tolerance and not `array_equal`, and that is measured:**
  `H` is exactly symmetric at `m = 1` — the [§8] path — and asymmetric by **1.776e-15** at `m = 2`,
  because the `beta` block is formed as `Xᵀ(wH)X` and that matmul is not bitwise symmetric for more
  than one column (§20c). An `exactly symmetric` assertion passes on every frame this stage has and
  fails the first time Stage 12 calls `polr` with a real covariate list. **Checked at three parameter values and on two frames**, not at the optimum alone:
  §21's own review note is that a wrong sign in the `alpha`-`beta` cross block which happens to vanish
  at one point would survive everything in §14, and one extra loop is what closes it. `hand_ordinal()`
  and `replication_frame()`, at `par = 0`, at the returned optimum, and at the optimum perturbed by
  0.5 in every coordinate.
- `polr` is asserted to equal a `scipy.optimize` maximiser of the same weighted objective to 1e-6 —
  measured 1.062e-07 with identical log-likelihoods to ten decimals (§18b). `scipy` is test-only by
  policy and this is its second use in the repository, Stage 6 §16b being the first.

### 14.4 The cutpoints, the collapse and the start values

In `test_model.py`.

- `fit.alpha` is asserted **strictly ascending** on every frame in the suite — the post-condition §5.4
  argues is guaranteed by construction, asserted on the return rather than asserted to be enforced.
- **The start values are asserted finite on a frame with an emptied middle category**, and the test
  states the coupling: they are finite *because* §5.3 collapses to positively-weighted categories, so
  every cumulative share is strictly interior. A companion computing the same start values over
  `MRS_LEVELS` instead is shown to produce `-inf` for the absent level's cutpoint.
- **A category present but carrying zero total weight is dropped**, and the resulting fit is asserted
  **equal to 1e-10** to the fit on the frame with those rows removed — measured. This is the
  difference between §5.3's rule and the observed-set rule, on the one input that separates them. **The
  test's comment states that the two fits do NOT share an O6 path**: O6 runs before the collapse, so
  the left-hand fit's rank was checked on the full design and the right-hand fit's on the reduced one
  (§5.5). Nobody should read this assertion as covering the collapsed design's rank, and §15 files
  what is left uncovered.
- **`_weighted_categories` is asserted to be the one definition**, by calling it directly and by an
  AST scan asserting that `np.unique` appears in `model.py` inside that function and nowhere else in
  the ordinal path. §5.3a is the argument; without the scan the two copies come back on the next edit
  and nothing fails.
- **The `-inf` branch of `_ord_loglik` is asserted to return `-inf`** on a hand-built unordered `alpha`,
  and the test records that **no input reaches it through `polr`**: 3000 attempts at reaching a
  crossing with pure Newton from a crowded start never lost the ordering (§20). Asserted directly on
  the private, in the posture Stage 7 §12.11 asserts B4 — unreachable on the declared path, tested by
  calling it.
- `converged_on` is asserted to be one of the two declared strings on every frame, and `iterations` in
  `[1, POLR_MAX_ITER]` — **closed at the top**, because `range(1, POLR_MAX_ITER + 1)` can return 200
  and a half-open assertion would fail on the last iteration that still converges.
- **`rescales` and `halvings` are asserted as VALUES on frames that drive them non-zero, not as
  integers on frames where they are zero.** This is the one place §14.4 departed from its own
  precedent, and the departure was measured (§21.1 item 14): with the assertions as drafted, **both
  safeguards can be deleted from `polr` and every frame in the suite returns bit-identically** —
  §5.2's trust region and its step-halving loop are on the bootstrap path, refit `N_BOOT` times, and
  had no test that could fail. Stage 6 §12.4a does not have this gap; `test_model.py:93`'s
  `big_first_step()` and `:105`'s `needs_halving()` are the fixtures it uses, and this stage needs its
  own two. They are one line each, because a proportional-odds estimate is equivariant under
  rescaling a covariate (§5.2) and the trust radius is relative, so **scaling the exposure is enough
  to drive both counters**. Measured on `orientation_frame()` (§20c):

  ```
      exposure scaled by      1.0        0.1        0.01
      rescales                  0          1           2
      halvings                  0          0           1
      iterations                6          6           8
      first_step_norm      4.8078    31.7083    315.0171
      beta * scale       4.0663825235  4.0663825234  4.0663825216
  ```

  So `polr_big_first_step()` is `orientation_frame()` with the exposure at 0.1 and
  `polr_needs_halving()` at 0.01, and the assertions are `rescales == 1`, `(rescales, halvings) ==
  (2, 1)`, `first_step_norm > C.POLR_MAX_STEP` on both, and `(0, 0)` on every well-behaved and every
  separated frame. **The same pair is the scale-invariance test**: `beta * scale` is asserted equal
  across the three to 1e-8 — not tighter, because `POLR_TOL` is absolute and the three fits stop at
  slightly different points (§5.2), which is itself the measurement.
- **The halving-exhausted raise is asserted under a monkeypatched `POLR_MAX_HALVINGS = 1`**, which is
  Stage 6 §12.4a's fourth assertion and was absent here. Without it the `for…else` branch of §5.2 is
  unreachable in the suite.
- **`first_step_norm` is asserted strictly positive and finite, and greater than 1e-8 on every frame
  that takes more than one iteration.** Strictly-positive alone is too weak: measured, a frame whose
  arms have identical weighted distributions (`y = [0,1,2,0,1,2]`, `a = [1,1,1,0,0,0]`) returns
  `iterations = 1`, `converged_on = "likelihood"`, `beta = 0.0` and `first_step_norm = 1.86e-16` — a
  value that passes the assertion and carries no information. `first_step_norm` is named rather than
  folded into "the counters" because it is the one field of the four that `_fit_detail` does not print
  (§3), so this assertion is the only thing keeping it alive.
- **§14.7's band probe is asserted to contain fits with a non-zero counter**, which is a free witness:
  measured, **10 of its 240 fits** carry one. An earlier draft of this section said "zero of both
  counters on every frame in §20, which is what makes a non-zero counter in a future run a signal
  rather than noise" — that sentence was contradicted by the suite's own probe, and it is struck.

### 14.5 The orientation, and the sign that every reference implementation has [roadmap]

- **The roadmap's criterion, on a construction whose truth is known.** On `orientation_frame()`
  (§14.0.4) — the treated arm's mass at mRS 0-3, the control's at 2-6 — `exp(beta)` is asserted `> 1`,
  and the construction's own truth is asserted first: the treated arm's mean mRS is asserted strictly
  below the control's, so the test states what "downward" means on its own data before it asks the fit
  about it. Measured **84.584014** on a mean category of 1.11 against 3.82 (§20).
- **The companion, which is what makes it a measurement.** The same fit read under the `alpha_k - beta*A`
  convention is asserted to be **strictly below 1** and to equal `1 / exp(beta)` to 1e-12 — the two
  readings are reciprocals by algebra, and asserting the reciprocal identity is what shows the
  companion is reading the *same fit* rather than a second one. Measured **0.011823** (§20), and
  asserted as a side-of-1 property rather than as that value, because the value is a property of
  `orientation_frame()`'s sampling and the direction is not (§14.0.4, §21.1 item 1). And a second
  companion applies `pilots/analysis.py:475`'s `np.exp(-res.params.iloc[0])` to **our** `beta` and
  asserts the direction reverses — the single most consequential do-not-lift in this stage (§7.2, §18).
- **The unweighted fit against `statsmodels.OrderedModel(distr="logit")`, with the asymmetry pinned by
  name.** On `reference_frame()`: `ours.beta + sm.beta` is asserted to be 0 to 1e-4 and
  `ours.alpha - sm.thresh` to be 0 to 1e-4 — measured 1.65e-05 and 3.16e-05 (§20). Both halves: **the
  coefficients are negated and the cutpoints are not**, and a test asserting only that the two fits
  "agree" would pass on a comparison of the cutpoints alone, which is exactly how the mistake survives
  a review. **A third assertion is what makes the first two able to fail**: `ours.beta - sm.beta` is
  asserted to be **larger than 1e-2**, so a frame on which both coefficients happened to be near zero
  cannot pass the sum test by having nothing to negate. The tolerance is `statsmodels`' `bfgs`, not
  ours, and the test says so — §14.3's replication oracle is the one that holds at machine precision.
- **`statsmodels` needs no `importorskip`**: it is a declared project dependency and not a `reference`-
  group one, checked rather than assumed (§2), as Stage 7 §12.12 checked the same thing.
- `_orientation()` is asserted to **raise** on a patched registry whose primary outcome has
  `higher_is_better=True`, and on one whose primary outcome is `kind="binary"` — so §7.1's assertion is
  shown to fire rather than assumed to.

### 14.6 The cumulative risk differences [roadmap, amended]

- **`RD_k` against a hand-computed weighted proportion difference**, on **`hand_ordinal()`**
  (§14.0.4), threshold by threshold. Six assertions, each naming its `k`, each against a literal
  computed by hand from that generator's two-decimal weights — which is what makes this the one
  assertion in the section that is a *value* rather than a tolerance, and it can be because
  `hand_ordinal()` holds no seed. Three earlier drafts said "§14.0.2's twelve-record frame"; §14.0.2's
  frame is `ordinal_cohort()`'s five-record cohort and its `RD_k` are already pinned there. §21.1
  item 2.
- **Non-decreasing in `k` within each arm**, on every frame in the suite, with the comment stating that
  this is **satisfied by construction** and therefore cannot detect the forbidden route (§8.2, §21
  amendment 3).
- **The companion that can fail**: the same six numbers computed through six weighted threshold models
  and asserted to **differ** from `cumulative_rd`'s. Different, not crossing — and the test's comment
  carries the measurement that makes "crossing" the wrong assertion: over 2000 samples of a 40-record
  seven-category frame the threshold route produced a strict decrease in **8 of 1990** fitted samples,
  0.40%, with the largest decrease below 5e-07, while the empirical route was non-monotone in **0 of
  2000** (§20).
- **`RD_6` is asserted structurally zero and asserted not to be reported**: `MRS_THRESHOLDS` stops at 5,
  `Primary.rd` has exactly six keys, and a direct computation at `k = 6` gives `+0.000000` while `k = 5`
  gives a non-zero value on a construction with deaths in both arms — measured `+0.049027` (§20). The
  pair is the assertion: without the second, a reader concludes the last reported threshold is the
  trivial one.
- **Ties are asserted where a level is unoccupied.** On `ordinal_cohort()`, `RD_1 == RD_2` and
  `RD_4 == RD_5` exactly, and the test names levels 2 and 5 as the reason. An implementation asserting
  strict monotonicity fails here and is wrong to.
- `weighted_proportion` returns `nan` and **not** `0.0` for an arm with no positive weight, and a
  companion shows `np.average` raising `ZeroDivisionError` when the guard is removed — Stage 7 §3.1's
  fact, one stage on, on the same expression.
- `cumulative_rd` is asserted to take exactly `(y, a, w)` by `inspect` — no threshold list — so a
  defaulted keyword fails the test rather than passing silently.
- **`_TREATED` and `_COMPARATOR` are asserted against the registry, not against 1 and 0**: `_TREATED` is
  asserted to be the code whose `TREATMENT_LABELS` value is the bridging label, `_COMPARATOR` the
  comparator's, and the two asserted distinct — so a third declared arm code fails here rather than
  silently making `RD_k` a contrast of the extremes. §8.1 is the argument; the sign of `RD_k` is the one
  thing in this module a reader cannot check by reading (§7).

### 14.7 Separation converges, and G7 is what turns it into a failure

The section that carries this stage's finding. In `test_model.py` for the fit's behaviour and in
`test_outcome.py` for the guard.

**This section's budget is stated because an earlier draft did not have one.** Three of its assertions
ranged over a **calibration** rather than a sample: the 4800-fit sweep of §6.3, 500 healthy `cond(-H)`
fits, and 400 sparse replicates run twice — about **6100 fits**, or **7 seconds** at §2's measured
1.157 ms, before generation. That roughly triples a suite that currently returns in about four, and the
cost is paid on every commit for the remaining six stages. Note what it is *not*: §14.4's 3000 crossing
attempts and §14.6's 2000 threshold samples were already comment-only, carried as measurements rather
than re-run, which is the posture this section did not have.

**The split is: calibrations live in §20 and in `test_config.py`'s literal, probes live in the suite**
(§21.1 item 3). Everything below runs in **under one second**, about 340 fits, and every number that
took a calibration to establish is guarded by the one assertion that actually catches its regression —
§12's `11.04 < POLR_MAX_ABS_BETA < 18.98` in `test_config.py`, asserted against both endpoints as
literals with §6.3 cited. A calibration re-run on every commit is not a stronger test than that
comparison; it is the same guarantee at a thousand times the cost, and it is the kind of cost that gets
paid by someone eventually deleting the test.

- **A perfectly separated frame is asserted to CONVERGE**, not to raise, on `separated_frame()`:
  `polr` returns, `converged_on` is `"score"`, `iterations` is 17, `rescales` and `halvings` are both 0,
  and `|beta|` exceeds 19. Every one of those is asserted as a **value**, and it can be because
  `separated_frame()` holds no seed (§14.0.4) — the finding is that nothing looks wrong, and a
  tolerance would be a weaker way of saying so. **40 fits.**
- **And the three detectors that do not work are asserted not to work**, each with its measured number:
  the fit converges; `cond(-H)` **is asserted to discriminate at equal cutpoint count and to be
  rejected anyway** — asserted as **orders of magnitude and not as ranges**, `cond(-H)` on a
  near-separated seven-category frame above `1e6` against a healthy seven-category one below `1e4`,
  with `separated_frame()`'s two-category `cond(-H)` asserted to sit *below* the healthy
  seven-category value, which is §6.3's size-dependence argument as a measurement. Measured 29.1 to
  162.8 healthy, 7.9e7 near-separated and 6.854 on the two-category construction (§20); asserted as
  thresholds because the healthy range is a property of the sampler. **20 healthy fits, not 500.**
  And the minimum fitted category probability is asserted to be 1.0 on perfect separation and is then
  asserted **not** to fall below 0.01 at one crossover, where `exp(beta)` is nonetheless above 1e8 —
  measured `min p = 0.05` and `1.53e+09` (§20). That last pair is the reason §6.3's guard is a
  magnitude bound, and it is asserted rather than argued.
- **The empty band is asserted as a band**, not as two endpoints, and **at the probe's budget rather
  than the calibration's**: `band_samples(BAND_EFFECTS, per_effect=20)` gives 240 fits in about 0.3 s,
  and the test asserts that the count landing strictly inside `(POLR_MAX_ABS_BETA - 4, POLR_MAX_ABS_BETA
  + 4)` is **0**, that some `|beta|` falls below the bound and some above it, and that
  `POLR_MAX_ABS_BETA` is strictly between the observed max-below and min-above. That is the band as a
  band, at 5% of the fits. **§6.3's 4800-fit calibration is not re-run here**: its endpoints 8.7873 and
  18.8055 are recorded in §20 and pinned as literals in `test_config.py` (§12), which is what fails if
  the bound is edited back to the 10.0 an earlier draft used — the failure mode the 4800 fits were
  guarding against, caught by one comparison instead of six seconds of fitting.
- **`primary` on the same frame is asserted to raise `model.FitError` naming G7**, with the message
  asserted to contain the odds ratio it rejected — so a reader of a Stage 10 failure log can tell a
  separated replicate from a non-converged one.
- **It raises `FitError` and not `SchemaError`**, asserted by type, because [§10] catches only the
  first to drop and count a replicate (§4.5). A companion asserts the same for G6.
- **The bound is read from `C.POLR_MAX_ABS_BETA` and never from a literal**, asserted by patching the
  constant low and watching a healthy fit start raising.
- **`_assert_reportable` bounds `beta` and not `alpha`**, asserted by **patching the bound low** rather
  than by finding a frame whose `alpha` exceeds 14.0 — because no such frame is constructible. `alpha`
  for a rare category grows like `log n`: measured (§20c), one rare record in 61 gives `max|alpha|`
  **4.16**, in 601 gives **6.40**, in 6001 gives **8.70**, in 60001 gives **11.00**, so reaching 14
  needs order 1e7 records. An earlier draft asserted that the test "constructs one", and it cannot
  (§21.1 item 20). The constructible form: on `rare_category_frame()` — `max|alpha|` **4.1599**,
  `|beta|` **0.126998** — the guard is asserted to **pass** with `POLR_MAX_ABS_BETA` patched to
  **3.0**, where a bound applied to `alpha` would raise, and to raise only when the bound is patched
  below `|beta|`. That asserts the field the guard reads, which is what the bullet is about.
- **Sparse replicates of the cohort's shape produce zero convergence failures**, asserted as a count —
  the measurement (§20) that makes §11's warning to Stage 10 concrete rather than cautionary. **Forty
  replicates in the suite, not four hundred**: the assertion is that the count of non-convergence
  failures is exactly 0 while the count of G7 raises with the guard in the path is greater than 0, and
  40 draws at 12 true effect sizes is enough for the second to be non-empty — measured over 400 in §20,
  which is where the number four hundred belongs. The point of this assertion is that the two counters
  mean different things, and that does not need four hundred samples to state.

### 14.8 O1-O6, and the two whose order is the point

In `test_model.py`. Each of the six raises `model.FitError` naming itself, parametrised over the six.

- **O3 before O5, asserted as an ordering.** A frame with **one `nan` weight** is asserted to raise O3;
  and a companion with O3 removed is asserted to **return a fit with one category fewer** — measured:
  four categories became three and `beta` moved from −1.775143 to −1.205433, with `iterations` 4 and
  `converged_on` `"likelihood"`. Both halves are required: the first alone passes for an
  implementation that raises somewhere, and the second is what shows what it is raising *instead of*.
- **O2's mechanism is asserted, not just its raise — and the companion must remove O2 AND O4.** With
  **both** removed, a `nan` response is asserted to produce a fit over **79 of 80** records with
  `fit.categories` unchanged and `beta` moving from −1.162237 to −1.199291 (§20c) — the record is
  dropped, not reported. **Removing O2 alone does not reach the fit**: `np.mod(nan, 1.0)` is `nan` and
  `nan != 0.0` is `True`, so **O4 catches a non-finite response too**, and the companion as drafted
  would have asserted a silent drop while watching O4 raise (§21.1 item 18). That interaction is worth
  stating in its own right: O2 is not the sole detector of a `nan` response, it is the one with the
  right message and the earlier position, and O4 is a backstop nobody designed as one. `np.nan > 0.0`
  is asserted `False` directly, which is the fact that explains both this and O3.
- **O4 fires on a non-integer response**, and a companion with it removed is shown to estimate one
  cutpoint per distinct value — on `noninteger_response_frame()`, which is **3 distinct values over 60
  records** and not the wider frame an earlier draft implied. Measured (§20c): with O4 removed that
  frame fits 3 categories and 2 cutpoints, one per distinct value; a frame with 6 distinct values over
  24 records instead raises `LinAlgError` from `np.linalg.solve`, which is not `FitError` and would
  have made the companion assert the wrong thing (§21.1 item 19). The generator is sized so the
  companion reaches a fit.
- **O6 fires on a design carrying its own intercept**, with the test's comment naming the asymmetry
  with `firth`, which prepends one. And O6 is asserted **not** to fire on a width-0 design — rank 0 of
  0 is not rank-deficient — which is what makes §4.5's `_assert_exposure_survived` necessary rather
  than redundant.
- Several at once are **not** collected into one error, unlike D1-D4 and B1-B5: the test asserts that
  a frame failing O1 and O6 reports **O1 only**, and the comment gives §5.5's reason.

### 14.9 The three audit entries, and the log

- The three `model` entries appear in `primary`'s declared order, at positions asserted against an
  index **captured before the call** — never a tail slice, which is the repair Stage 5 §10 landed and
  Stages 6 and 7 carried.
- **`data.KINDS` is asserted to be the declared nine, unchanged**, and `_HEADINGS` to render nothing
  new. This is §9.1's claim in test form, and it is the assertion that fails if someone adds an
  `estimate` kind.
- All three render under `## Fitted models`, after `overlap_weights`.
- **Each entry's `n` is asserted to be a different quantity**: entry 1's is records removed by the
  outcome mask, entry 2's is records estimated, entry 3's is `len(MRS_THRESHOLDS)`. The test asserts
  the three are not the same number on a frame where they differ, because that is the property a
  reader needs.
- **Entry 1 is asserted to be recorded even when its `n` is 0** — the v7 case — and asserted to name as
  many records as it says it removed on the §14.2 frame where it removes one.
- **Entry 2's table has one row per `MRS_LEVELS` entry**, including the absent ones, and the absent
  rows are asserted to render `missing` in the estimate column with the absence status. On
  `ordinal_cohort()` two such rows are asserted by level.
- **The status of the top FITTED level is asserted by name, on `truncated_cohort()`** (§14.0.5) — the
  one frame in the suite where the top fitted level is not the top declared one. `fit.categories` is
  asserted `(0, 1, 2, 3, 4)` first; then `mRS <= 4` is asserted to carry the highest-fitted status and
  **not** `"fitted"`, and `mRS <= 6` to carry the absence status **and** to say it is the highest
  declared level. This is the branch §9.2 corrects, and it is worth its own bullet because the previous
  bullet **passes on the broken table**: the row that rendered `missing` under status `fitted` is not
  an absent row, so an assertion about absent rows never looks at it.
- **The estimate column and the status column are asserted to agree, as an invariant**: a row reads
  `missing` in the estimate column **if and only if** its status is one of the three no-cutpoint
  strings. That one assertion catches every miswiring of §9.2's branch, including ones nobody
  anticipated, and it holds on all four frames.
- The status column is asserted to contain only the **six** declared strings — the three no-cutpoint
  variants, `"fitted"`, and the two coefficient-row statuses — and never `True`, `1` or `nan`. **Six,
  not four**: an earlier draft of this bullet said four while `_fit_table` emitted five, and §9.2's
  correction makes it six (§21.1 item 4). The count is asserted against a literal tuple of the six, so
  a seventh status added without a test fails here.
- **No cell of any of the three tables contains a `|`, and each table's rendered header, separator and
  body have the same markdown cell count.** The companion pastes a `|` into a header and shows the
  counts diverge, so the assertion is known to fire. Byte-identity does not catch this (Stage 7 §21.1),
  which is why it is its own assertion.
- **Both `detail` strings are asserted as SENTENCES and not only as counts**, under every branch they
  have: entry 2's absent-levels clause with 0 and with 2 absent levels, and entry 3's per-arm totals
  clause. Stage 7 §21.5a found an unterminated clause in exactly this position that three fence
  assemblies and two review rounds walked past, because the branch never rendered — so this stage
  asserts the branch's prose from the start, on the frame that reaches it.
- The log renders byte-identically twice in one process, and under `PYTHONHASHSEED=0` versus `1` in a
  subprocess — the two-seed driver of Stage 5 §12.9, extended one stage.

### 14.10 The estimate and the risk differences are one contrast

- `Primary.in_estimate` is asserted to be the **same object's** mask as the one `cumulative_rd` was
  called over, by asserting `Primary.fit` equals `model.polr` on the masked triple **and** `Primary.rd`
  equals `cumulative_rd` on the same three arrays. §4.2's "bound once" as an assertion rather than a
  code comment.
- **The orientation is asserted to agree between the two scales — through the weighted mean, NOT
  threshold by threshold.** `sign(beta)` is asserted to equal `−sign(mean_w(Y | treated) −
  mean_w(Y | comparator))`, over the same mask and the same weights: under §5.1's parametrisation a
  positive `beta` moves mass to LOW mRS in the treated arm, so the treated arm's weighted mean mRS
  must be the lower one. A sign error in either scale alone breaks this, and neither alone would be
  caught by §14.5 or §14.6. Measured on both fixture frames (§20c): `ordinal_cohort()` beta +1.7516
  with means 2.2513 against 4.4137, `truncated_cohort()` beta +0.3448 with means 2.2513 against
  2.5288 — the relation holds on both.

  **The per-threshold form of this assertion is FALSE and an earlier draft asserted it** (§21.1
  item 13). It read: "on a frame where `beta > 0`, the `RD_k` are non-negative at every `k` where the
  two arms' cumulative probabilities differ." That is not a property of a proportional-odds fit. `beta`
  is a single summary of a shift the data need not exhibit at every cut, which is exactly the
  assumption §15 files as untested — so where proportional odds fails, `beta` and an individual `RD_k`
  legitimately disagree in sign. Measured: on `truncated_cohort()`, `beta = +0.3448` while
  `RD_3 = −0.4531`, because both control records sit at mRS ≤ 3 so `P(Y ≤ 3 | comparator)` is exactly
  1.0. It fails on `hand_ordinal()` too — `beta = −0.8696` with `RD_0` and `RD_5` both `+0.0476`. **The
  suite's own fixtures are counter-examples**, so the assertion as drafted would have failed at T6 on
  correct code, and the two natural repairs — restrict it to `ordinal_cohort()`, or weaken it to a
  majority-sign check — both discard the property it exists to test.
- **The strict per-`k` form is kept, on `orientation_frame()` only**, where the construction's truth is
  known and proportional odds holds by design: there every `RD_k` is asserted to carry `beta`'s sign.
  That is the frame on which the strong form means something, and §14.5 already fits it.
- **Entry 3's per-arm weight totals are asserted to equal Stage 6's `overlap_weights` per-arm `sum w`**
  on the workbook, read out of the **same** `Audit` (§14.0.3) and compared **through `_fmt`** — Stage 7
  §12.8's rule, because an `AuditEntry.table` holds strings and a full-precision float will not equal
  one. Data-gated. On v7 the two agree because `in_estimate` is `in_model`; on a frame where they differ
  the test asserts they differ **and** that entry 1 explains the gap.

### 14.11 Stage 8 adds nothing, refits nothing, and reads no diagnostic

- The frame passed in is returned unchanged, cell for cell, with the same columns, dtypes and index —
  the four assertions of Stage 5 §12.14, Stage 6 §12.13 and Stage 7 §12.10, one stage on.
- The `Propensity` is unchanged: `e`, `w` and `in_model` are asserted equal to the objects passed in,
  and `primary` is asserted to add no attribute to the frame.
- An AST scan asserts `outcome.py` imports neither `balance` nor `statsmodels`, `sklearn` or `scipy`,
  and names `propensity.fit` nowhere. **The `balance` half is the one that matters**: importing it
  would not fail, and Stage 7 §11 states that Stage 8 does not read the `Balance`. Companion test
  included for each.
- **`primary` is asserted to take exactly `(df, ps, audit)` by `inspect`** — no outcome name, no
  covariate list, no threshold list, no bound — so a defaulted keyword fails the test rather than
  passing silently.
- **`Primary` and `PolrFit` are asserted to carry exactly their declared fields**, by `inspect`, and
  asserted to carry **no field whose name contains `se`, `cov`, `ci`, `pval` or `interval`**. §5.6 is
  the argument and this is the assertion that makes it survive a later edit.
- **`Primary` is asserted to have no method returning a verdict** — no `significant`, no
  `favours_bridging` — with §3's reason in the comment.
- With G6 removed, a one-armed frame is asserted to produce **either** an `IndexError` on `beta[0]`
  **or** an odds ratio of exactly 1.0 depending on how the result is read, and the test asserts both
  readings on the same input to show that neither is a failure the caller can see (§4.5, §10).
- `outcome.py` and `model.py` contain no bare `assert`, name no raw header, and are not in
  `EXEMPT_FROM_RAW_NAME_SCAN`, each with its scan-fires companion.
- **`model.py`'s existing behaviour is unchanged**: `test_model.py`'s Stage 6 sections and
  `test_propensity.py` are asserted green **unedited**, which is DoD-3 and is where §12's
  "no existing function changes" is checked rather than believed.

### 14.12 Properties from the workbook [data-gated]

**This section asserts properties and not values, and §4.3 is why.** Stage 7 §12.12 pinned 1.343 and
0.380; the analogous numbers here are the study's result. Every assertion below can fail.

- 93 / 92 / 92 for the three populations, and entry 1's `n` is 0. These are **counts**, not estimates,
  and they are pinned.
- All seven declared mRS levels are occupied, so `len(fit.alpha)` is 6 and `fit.categories` is
  `MRS_LEVELS`. Pinned, with §14.1's reason.
- **The fit converged**, `converged_on` is one of the two declared strings, and `iterations` is below
  `POLR_MAX_ITER`.
- **`alpha` is strictly ascending** and has exactly 6 elements.
- **`|beta|` is strictly below `POLR_MAX_ABS_BETA`**, so G7 does not fire on v7 — asserted, because a
  guard that fires on the real data would mean the primary analysis has no estimate, and that is a
  finding for the PI and not a passing test.
- **`odds_ratio` is finite and strictly positive**, and equals `exp(beta)` to 1e-12.
- **`RD_k` is non-decreasing in `k` within each arm**, all six thresholds, both arms.
- **Each arm's `P(Y <= 5)` is strictly below 1.0 or exactly 1.0 according to whether that arm has a
  death**, which ties §8.3's statement to the workbook without quoting a difference.
- **The orientation agrees between the scales, through the weighted mean**: `sign(beta)` equals
  `−sign(mean_w(mRS | bridging) − mean_w(mRS | EVT alone))` over `in_estimate`. This is the strongest
  property available without quoting a magnitude, and it is the one that catches a sign error
  introduced downstream of §14.5's synthetic construction. **It reports no magnitude and no
  distribution split by arm**: a sign is not an estimate, and §4.3 permits it for the reason §20's
  probe may count arms.

  **NOT the per-threshold form**, and this matters more here than in §14.10 because this section runs
  on the workbook. A draft asserted that every `RD_k` carries `beta`'s sign wherever the arms differ;
  that is false wherever proportional odds fails (§14.10, measured on both fixture frames). Asserted
  on v7 it could go **red for a legitimate reason** — and its failure message would say "sign error"
  when the content is a failed modelling assumption, on the one stage where §4.3 makes that
  impossible to distinguish in advance. A data-gated test that can fail for a reason the analysis
  already declares untested is worse than no test: it converts a known limitation into an alarm at
  the exact moment the PI is reading the estimate for the first time (§21.1 item 13).
- **Entry 3's per-arm totals reconcile with Stage 6's `overlap_weights` through `_fmt`** (§14.10).
- **The three entries render, and their `detail` strings' interpolated counts match their tables** —
  entry 2's "6 cutpoint(s)" against a table with six `fitted` rows, entry 3's threshold count against
  six rows.
- **No assertion in this section quotes `beta`, `exp(beta)` or any `RD_k`**, and a meta-assertion in
  `test_outcome.py` scans this section's source for a float literal outside an allowlist. That scan is
  what keeps §4.3 true as the file is edited, and it has a companion that fires on a pasted magnitude.
  **Its boundary and its allowlist are specified here rather than left to the implementer**, because a
  scan whose scope is "this section" is either unwritable or vacuous:

  - **Boundary.** The section is the source between its own banner comment, `# --- 14.12`, and the
    **next `# ---` banner at column 0** — next in file order and not next in section number, because
    `test_balance.py:261,364` shows the banners are not in numeric order and a scan keyed on `14.13`
    would silently run to the end of the file. That banner convention already exists in
    `test_cohort.py`, `test_model.py`, `test_propensity.py` and `test_balance.py`, so the delimiter is
    a fact about the repository rather than something this test invents. **The scan asserts both
    delimiters were found, and that the range is non-empty, before it asserts anything about what is
    in it** — a scan that silently matched nothing is the green-test failure
    `test_the_raw_name_scan_actually_fires` (Stage 1 §7) exists to prevent, three modules earlier.
  - **Allowlist**, and it is short by design. An `ast.Constant` float is permitted only when it is
    (a) the value of a `rel=` or `abs=` keyword, or the second positional argument of `pytest.approx`;
    or (b) `0.0` or `1.0`, which are the two boundary values §14.12 asserts `odds_ratio` and the
    cumulative probabilities against, and which carry no magnitude. Everything else fails. Integers are
    unrestricted: 93 / 92 / 92, the 6 cutpoints and the 7 levels are counts and §14.12 pins them
    deliberately.
  - **It is an AST walk, so a float in a `#` comment is invisible to it and a float in a docstring is
    not.** That asymmetry is deliberate and is stated because the natural place to paste `exp(β)` while
    reading the log is a docstring, and catching it there is the point. A comment is not covered, which
    is a real hole and a small one: §15's first bullet and DoD-17's read-by-eye are what stand behind
    it, and a comment carrying the estimate does not make an assertion depend on it.

### 14.13 The module boundary

- An AST scan asserts `model.py` names neither the treatment nor any outcome: no `C.` attribute
  beginning `PS_`, `OUTCOME_` or `PRIMARY_`, and no string literal containing `ivt` or `mrs`. Stage 6
  §12.12's scan, extended to the constants this stage adds — `POLR_*` passes it and an
  `OUTCOME_MAX_ABS_BETA` would not, which is §12's reason for the prefix.
- `_assert_polr_fittable`'s six raise strings are asserted to name neither an exposure nor an outcome,
  parametrised over the six, with the companion asserting the scan fires on a pasted mention.
- `outcome.py` is asserted to be the only module naming `PRIMARY_OUTCOME`, and `model.py` asserted not
  to.

### 14.14 The preconditions

- G1 fires on a `Propensity` whose Series carry a different index, **raising `SchemaError` with G1's
  message**, and it is **parametrised over the three Series** — `e` alone, `w` alone, `in_model`
  alone — with the message asserted to name which one. Stage 7 §21.4's `e`-alone finding is inherited
  rather than re-learned. **It is also parametrised over the four kinds of misalignment** — a superset
  index, a subset, a permutation and a dtype change — all four measured to fail `index.equals` (§20c),
  so G1 is complete for §10 rather than assumed to be.
- **And a companion shows what G1 PREVENTS, which no earlier draft did.** Every other guard in this
  stage has a with-it-removed witness producing the wrong number — O2, O3, O4, G6, G7 — and G1 had
  none, while its own message claims it is the difference between an estimate and "a finite common
  odds ratio for a population that does not exist". The witness is a **permutation**, not `e`:
  measured (§20c), with G1 removed a `Propensity` carrying the frame's labels in reversed order gives
  `odds_ratio` **58.01** against the correct **5.76** on `ordinal_cohort()`, with nothing raising,
  because `w` arrives ordered by the `Propensity`'s index while `y` and the arm vector come from the
  frame's. A **subset** index is the second witness and is quieter: the missing labels resolve to
  `False` and the record leaves the [§11] denominator in silence. **A superset is benign** — measured,
  it returns the correct estimate — so G1 is stricter than §10 strictly needs, which is the right side
  to err on and is asserted as such. `e` alone is the one case that CANNOT produce a wrong estimate,
  because §10 never reads `e`; the test says so where it parametrises over the three (§21.1 item 17).
- G2 fires on a mask cast to `object` with a missing value in it, **parametrised over `None`, `pd.NA`
  and `np.nan`**; and its companion asserts what actually happens without the phase-1 raise, which is
  **nothing**. Measured (§20c): `object & bool` coerces the missing entry to `False` for all three, so
  `in_estimate` comes back a clean boolean Series of the frame's length with **4 records estimated
  against `in_model`'s 5** and no exception anywhere. An earlier draft said the input "dies with
  pandas' `ValueError`"; it does not, and the companion asserts the silent exclusion instead (§21.1
  item 16). **That makes G2 stronger, not weaker** — a guard whose absence is silent is a guard worth
  having, and one whose absence raises anyway is nearly redundant.
- **G3 fires on a frame missing `PRIMARY_OUTCOME`** and on one missing `TREATMENT`, and its companion
  shows that with G3 in phase 2 the outcome case dies at phase 2's first line with
  `KeyError('mrs_90d')` and G3's message discarded (§4.4a). That companion is the one that would have
  caught Stage 7's incomplete first repair, applied here before the fact.
- G4 fires on a frame whose `in_estimate` leaves one arm, with its message asserted to name both arm
  labels and their counts.
- G5 fires on a **hand-built `Propensity`** whose `in_model` includes a record with a `nan` weight —
  unreachable from `propensity.fit`, so constructed, in the posture Stage 7 §12.7 declares.
- G6 fires on a one-armed frame's design and on a hand-built width-0 design, raising **`model.FitError`
  and not `SchemaError`**, asserted by type (§4.5).
- **A failing `primary` records the [§11] denominator and nothing else — never a partial estimate.**
  Two halves, and the bullet is titled for both because an earlier draft titled it "records nothing"
  and then explained in its own body that it does not (§21.1 item 5). Asserted for every one of G1-G7
  by capturing `len(audit.entries)` before the call and reading it after the raise:
  - **G1-G5 append nothing at all.** They are preconditions and they run before any `audit.record`.
    Delta is exactly 0, asserted for each of the five.
  - **G6 and G7 append exactly entry 1, and never entries 2 or 3.** Delta is exactly 1, the one entry
    is asserted to be `outcome_completeness` by step name, and `audit.entry("model", _STEP_FIT)` is
    asserted `None`. Both fire after entry 1 would naturally be written, so §10 orders
    `_record_outcome_completeness` before the design deliberately: the [§11] denominator is a fact
    about the frame established before any modelling choice, and a log that names it is useful
    precisely when what follows fails (`propensity.py:440-443`). This is the one place this stage
    departs from Stage 7's "records nothing", and the departure is the point rather than a leak.
- **A frame failing G1 and G4 together reports only G1**, and the message says G4 and G5 were not run.
- Several at once produce one `SchemaError` naming each, **within a phase**: G1 and G2 together give
  one naming two; G4 and G5 together give one naming two; and the cross-phase case is the bullet above.

### Coverage map

```
model.py, amended                                  test_model.py, amended
├── _weighted_categories(y, w)                      the ONE positive-weight rule (§5.3a)
│   ├── the rule, called directly                   ├── 14.4
│   └── np.unique appears ONLY inside it            └── 14.4  AST scan; O5 and polr share it
│
├── _ord_pieces / _ord_loglik / _ord_score_hess     written out in §5.2
│   ├── analytic score vs central differences       ├── 14.3  2.784e-08, 3 points x 2 frames
│   ├── analytic Hessian vs central differences     ├── 14.3  2.812e-08, symmetric, neg definite
│   ├── -inf on a non-positive category probability ├── 14.4  direct; UNREACHABLE via polr (§5.4)
│   └── the clip is never approached                └── 14.4  max |eta| 18.2 vs POLR_ETA_CLIP 500
│
├── _assert_polr_fittable  O1…O6
│   ├── O1 design non-finite                        ├── 14.8
│   ├── O2 response non-finite → DROPS a record     ├── 14.8  + 80 fitted as 79 without it
│   ├── O3 weight non-finite → DELETES a category   ├── 14.8  + 4 categories became 3, beta moved
│   ├── O3 BEFORE O5, asserted as an ordering       ├── 14.8  `nan > 0.0` is False
│   ├── O4 non-integer response                     ├── 14.8  one cutpoint per distinct value
│   ├── O5 fewer than two weighted categories       ├── 14.1  the committed fixture reaches it
│   ├── O6 rank, and the intercept asymmetry        ├── 14.8  + does NOT fire on width 0
│   └── first failure raises, never collected       └── 14.8  O1+O6 reports O1 only
│
└── polr(X, y, w)                                   written out in §5.2
    ├── integer weights == row replication          ├── 14.3  0.000e+00  [roadmap, amended]
    ├── all-ones == None, bit for bit               ├── 14.3  the companion that makes it needed
    ├── weighted != unweighted                      ├── 14.3  1.8268, OPPOSITE SIGNS
    ├── vs scipy on the same objective              ├── 14.3  1.062e-07
    ├── vs statsmodels unweighted, SIGNS            ├── 14.5  beta NEGATED, alpha NOT
    ├── alpha strictly ascending on return          ├── 14.4  every frame
    ├── start values finite BECAUSE of the collapse ├── 14.4  + -inf over MRS_LEVELS
    ├── zero-weight category == rows removed        ├── 14.4  1e-10; NOT the same O6 path (§5.5)
    ├── all four iteration diagnostics recorded     ├── 14.4  incl. first_step_norm (§3)
    ├── SEPARATION CONVERGES, 17 iters, score       ├── 14.7  counters 0, nothing out of range
    ├── the three detectors that do not work        ├── 14.7  as thresholds, 20 fits not 500
    └── the empty band IS a band                    └── 14.7  240-fit probe; 4800 is §20's

outcome.py                                         test_outcome.py
├── weighted_proportion(x, a, w)
│   ├── nan for a zero-weight arm, not 0.0          ├── 14.6  + ZeroDivisionError without it
│   └── keyed by TREATMENT_LABELS, both arms        └── 14.6
│
├── cumulative_rd(y, a, w)
│   ├── RD_k vs hand-computed, six times            ├── 14.6
│   ├── non-decreasing — BY CONSTRUCTION            ├── 14.6  0 of 2000; cannot detect the route
│   ├── the threshold-model route DIFFERS           ├── 14.6  the companion that can fail
│   ├── ties where a level is unoccupied            ├── 14.6  RD_1 == RD_2, RD_4 == RD_5
│   ├── RD_6 structurally 0 and not reported        ├── 14.6  + RD_5 non-zero with deaths
│   └── takes no threshold list                     └── 14.6  by inspect
│
├── _assert_primary_inputs  G1…G5                   ├── 14.14 all five + several at once
│   ├── the phases are what let G1/G2/G3 raise      │         + the KeyError companion on the
│   │                                                │           OUTCOME column (§4.4a)
│   └── a failing primary records NOTHING            │
│
├── _assert_exposure_survived  G6                   ├── 14.14 width-0 design; FitError not Schema
├── _assert_reportable  G7                          ├── 14.7  the separation guard
├── _orientation()                                  ├── 14.5  raises on a patched registry
│
├── the three audit entries                         ├── 14.9  order, KINDS is NINE, three n's,
│                                                    │         no `|`, detail as SENTENCES,
│                                                    │         SIX statuses, estimate-iff-no-cutpoint,
│                                                    │         top FITTED level (§14.0.5)
└── primary(df, ps, audit)                          written out in §10
    ├── in_estimate = in_model & present            ├── 14.1
    ├── the [§11] denominator is PER ESTIMATE       ├── 14.2  CONSTRUCTED — no frame varies it
    ├── one contrast, two scales, one mask          ├── 14.10
    ├── frame and Propensity unchanged              ├── 14.11
    ├── no se / cov / ci / pval field anywhere      ├── 14.11 by inspect  (§5.6)
    ├── signature is exactly (df, ps, audit)        ├── 14.11 by inspect
    ├── imports no diagnostic, and NOT balance      ├── 14.11 AST scan + companion
    └── workbook: PROPERTIES, never values          └── 14.12 + the float-literal scan (§4.3)

config.py amendments                               test_config.py
├── PRIMARY_OUTCOME computed, unique                ├── six assertions (§12)
├── MRS_LEVELS computed; [:-1] == MRS_THRESHOLDS    ├── the second is what stops the drift
└── POLR_MAX_ABS_BETA inside the sparse region      └── 11.04 < 14.0 < 18.98, over 2-6 cutpoints

R oracle                                           test_reference_r.py
└── ordinal::clm with weights                       └── §18b — SPECIFIED, and skipping: Rscript
                                                            segfaults at startup on this machine

  polr 16 routes, _weighted_categories 2, _ord_* 4, O1-O6 6, weighted_proportion 2,
  cumulative_rd 6, G1-G7 9, entries 12, primary 8, config 6  =  71
  71 with a test that can fail.  0 unreachable-by-construction, 5 unreachable ON EVERY FRAME
  THIS STAGE HAS — G5, unreachable while in_estimate is built from in_model (§4.4), tested by
  a hand-built Propensity; the PER-ESTIMATE denominator, which needs an outcome-incomplete
  in_model record and gets one only from the frame §14.2 builds; `_ord_loglik`'s -inf branch,
  measured unreachable through polr in 3000 attempts (§5.4), tested on the private directly;
  G6, which needs a design `design` narrowed and gets one only from a one-armed frame; and
  §9.2's top-FITTED-level status, which needs a fit whose highest occupied category is not
  MRS_LEVELS[-1] and gets one only from §14.0.5's truncated_cohort().
  All five are tested on constructed input rather than asserted unreachable, because all five
  are reachable on a workbook that has not arrived — and the second and the fifth are reachable
  on a ONE-LINE edit to the frame, which makes leaving them untested the least defensible. The
  fifth is also reachable in an ORDINARY [§10] REPLICATE, one that happens to draw no death,
  which is what took it from a note to a fixture (§21.1 item 6).

  If the map and this count disagree, the map is authoritative and this block is stale.
```

### Failure modes

| Failure | Test | Error handling | Visible or silent |
|---|---|---|---|
| A separated fit reported as an odds ratio of 6.5e15 | 14.7 | **G7** (§6.3). Nothing else fires: it converges in 17 iterations on the score criterion, every counter is 0 and every fitted probability is finite | **silent** past every other check, and the only stage-8 failure that would reach a manuscript |
| A bootstrap percentile interval whose upper limit is a quantile of a distribution with a tail at exp(21) | 14.7 | G7 raises `FitError`, which is the type [§10] drops and counts | **silent** — the p-value barely moves and the interval moves a lot, so the two disagree where [§10] says they agree |
| Near separation missed by a guard that only catches perfect separation | 14.7 | §6.3's bound is on the coefficient, not on the fitted probability: at one crossover `min p` is 0.05 and `exp(beta)` is 1.5e+09 | **silent** — and it is the case a bootstrap actually draws |
| The primary result reported with its sign inverted, by lifting `exp(-beta)` from the pilots | 14.5 | §7.1's parametrisation is [§14a]'s own; §7.2 records that the coefficients are negated and the cutpoints are not | **silent** — four cutpoints agree to five decimals, so the two fits look identical |
| One `nan` weight silently deleting an outcome category and moving `beta` | 14.8 | **O3 before O5** (§5.5). `nan > 0.0` is False, so the positive-weight rule drops the category | **silent** — measured: 4 categories became 3, `beta` −1.775 → −1.205, `iterations` 4, nothing raised |
| One `nan` response silently dropping a record | 14.8 | O2. `np.unique` returns `nan` as a level and its total weight is 0 | **silent** — 80 records fitted as 79 |
| A one-armed frame producing an intercept-only fit with no treatment coefficient | 14.14, 14.11 | **G4** on the data and **G6** on the matrix. `design`'s constant-column rule drops the sole predictor and `polr` fits the cutpoints alone | **silent on both sides** — neither `design` nor `polr` raises, and the caller gets `IndexError` or `1.0` |
| The coefficient and the risk differences computed over different populations | 14.10 | §4.2: the mask is bound once in §10 and both read it | **silent** — both numbers are finite and each is correct for its own population |
| An outcome-incomplete weighted record in the `Σw` denominator but not the numerator | 14.2 | §4.1's mask is `in_model & present`, never `present` alone | **silent** on v7, where no record is one |
| A model-based standard error reaching a table | 14.11 | §5.6: no field exists, and `inspect` asserts the field list | visible in review only — and one attribute access away without the assertion |
| `RD_k` from six threshold models rather than the weighted distributions | 14.6 | §8.2. The monotonicity criterion cannot catch it — asserting the numbers against a hand computation, and against the threshold route, can | **silent** — measured, the threshold route crosses in 0.40% of samples and by less than 5e-07, so a monotonicity check passes on it |
| A cutpoint labelled with a threshold it is not, on a collapsed fit | 14.0.2 | `cut_levels` is `fit.categories[:-1]`, not `MRS_THRESHOLDS` (§10) | **silent** — the numbers are right and the labels are wrong |
| An absent mRS level vanishing from the log, so nobody knows it was considered | 14.9 | §9.2's table ranges over `MRS_LEVELS` and marks the absent ones | silent without it |
| `RD_5` read as the structural zero, so the survival contrast is discarded | 14.6 | §8.3: `RD_6` is the structural zero and is not reported | visible in review only |
| A table rendered with a header row holding more markdown cells than its body, from a `\|` in a header | 14.9 | Headers are written without the character and §14.9 asserts its absence | **silent** past byte-identity across hash seeds (Stage 7 §21.1) |
| A `detail` clause running into the sentence after it, on a branch no available frame renders | 14.9 | §14.9 asserts each branch's SENTENCE, on the frame that reaches it | **silent** — Stage 7 §21.5a is the precedent, found after three assemblies and two review rounds |
| A misaligned `Propensity` producing a finite common odds ratio for the wrong records | 14.14 | **G1**, on the index and not the length | visible **only** because of G1 |
| G1/G2/G3 detected and then discarded, the caller seeing a bare pandas exception | 14.14 | §4.4's phase 1 raises before G4 builds the mask | visible — it raised either way, with the wrong message |
| A precondition firing after entry 1 is recorded, leaving a log describing half an estimate | 14.14 | All five preconditions run before any `audit.record`; §14.14 asserts a delta of exactly 0 on G1-G5 and exactly 1 — entry 1 alone, by step name — on G6 and G7 | visible — it raised — but the log is already written |
| A log carrying the coefficient of a fit that was then rejected | 14.14 | §10 orders `_assert_reportable` between the fit and `audit.record`, so a G7 fit is never written; the G6/G7 delta assertion is what pins it | **silent** — the entry would render under `## Fitted models` beside a `missing` estimate and read as a partial run |
| Entry 2 rendering `missing` under status `fitted`, on a fit whose top category is not the top declared level | 14.9 | §9.2's branch keys on `fit.categories[-1]`, never on `MRS_LEVELS[-1]`; §14.9's estimate-column-iff-no-cutpoint invariant and §14.0.5's frame | **silent past every other assertion** — no frame in the suite witnessed it, and §14.9's absent-row assertion passes on the broken table |
| Two copies of the positive-weight category rule disagreeing, so O5 counts a different `K` from the one `polr` fits | 14.4 | §5.3a: `_weighted_categories` is one function with two callers, and an AST scan asserts `np.unique` appears in the ordinal path only inside it | **silent** — both counts finite, both plausible, and the nan subtlety §5.5 turns on lives in the expression |
| An edit moving `beta` in the fourth decimal, caught by nothing | 14.0.2 | §14.0.2's golden vector, on a 5-record cohort with 4 cutpoints | **structural, and NOT fully mitigated** — §4.3's cost, stated in §15 |
| The proportional-odds assumption failing, so `beta` is a summary of a non-constant effect | — | Nothing here. [§8] specifies no test of it and §15 files it | **structural, not asserted anywhere** |
| **A test asserting that `beta` and every `RD_k` share a sign** — which is the assumption above, asserted as though it were arithmetic | 14.10, 14.12 | §14.10 asserts the agreement through the WEIGHTED MEAN instead, which is a monotone functional of the same distributions. Measured: `truncated_cohort()` has `beta` +0.3448 and `RD_3` −0.4531 | **it fails on CORRECT code**, and on the workbook it would fail with a message saying "sign error" for a reason §15 already declares untested |
| A `nan` response reaching the fit with O2 removed, believed to be O2's own mechanism | 14.8 | O4 catches it too — `np.mod(nan, 1.0)` is `nan` — so the companion removes **both**. O2 is the better message and the earlier position, not the sole detector | visible in review only, and it made a specified companion unwritable |
| A safeguard on the bootstrap path with no test that can fail | 14.4 | §14.4 asserts `rescales` and `halvings` as **values** on fixtures that drive them, not as integers on fixtures where they are zero. Measured: with the drafted assertions, deleting both safeguards left every frame bit-identical | **silent** — the suite stays green while `polr` loses its trust region and its halving loop |
| A `Propensity` whose labels are permuted rather than absent, producing a plausible estimate for mis-weighted records | 14.14 | **G1**, on `index.equals`. The companion measures what it prevents: `odds_ratio` 58.01 against 5.76 with G1 removed | **silent** — every number finite, and §10 reads by label so nothing misbehaves visibly |
| A three-valued mask silently removing a record from the [§11] denominator | 14.14 | **G2**. `object & bool` coerces the missing entry to `False`, so without G2 the record leaves the denominator with **nothing raising** — measured for `None`, `pd.NA` and `np.nan` alike | **silent**, and an earlier draft expected a `ValueError` that does not occur |
| `beta` compared across bootstrap replicates with different numbers of cutpoints | — | §5.3, §11. Stage 10's, and it is handed over rather than left to be found | **structural** — under proportional odds they are the same parameter, which is the assumption above |

Every failure above is visible at the point it occurs, and **five** are visible only because of a rule
this stage adds: the separated fit (G7), the deleted category (O3's ordering), the vanished exposure
(G6), the two-population contrast (§4.2's single mask), and the misaligned `Propensity` (G1). The
first is the one the roadmap's **Accept when** did not anticipate, and the one §21's amendment 1
exists for.

## 15. Known gaps carried forward

- ~~**There is no pinned regression number on the primary effect anywhere in git, and that is §4.3's
  cost.**~~ **CLOSED 2026-08-27, and not by the route this bullet predicted.** The gap was real: an
  edit that moves `β` in the fourth decimal without changing convergence, ordering, orientation or
  monotonicity passed the whole suite, and §14.0.2's golden vector is a partial mitigation only —
  five records and four cutpoints against the workbook's ninety-two and six, so it exercises the same
  arithmetic at a different scale and cannot witness a defect that appears only at the workbook's.

  **What this bullet got wrong was the FORM of the fix.** It said the numbers *"become legitimate
  regression pins"* once the analysis is locked, which read "pin" as "write the number down" — and on
  that reading the gap needed a PI decision, an [§15] amendment and the dismantling of §14.12a, the
  meta-assertion built to keep the estimate out of git. None of that follows. The pin is a **content
  hash**, which is `config.DATA_SHA256`'s own pattern one stage on: `test_outcome.py`'s
  locked-estimate section renders the primary through `data._fmt` in a declared order and asserts a
  SHA-256 of it. It fires on any drift at the precision the audit log publishes, it localises to the
  stage, **no estimate enters version control**, and §14.12a needs no change because there is no float
  literal to scan — asserted, so that stays true.

  **So the pin did not need the lock and could have existed since this stage landed.** DECISION 9 and
  [§15]'s amendment of 2026-08-27 stand on their own — [§16] was citing [§15] for a rule [§15] did not
  contain — but they are not what this bullet was waiting for. The baseline the digest is checked
  against lives beside the log in `../out/stage0_data_inventory.md`, which is where the numbers stay.
- **The proportional-odds assumption is not tested, and this stage adds no test of it.** [§8] specifies
  a common odds ratio and gives its own reason for not augmenting each threshold separately; it does not
  ask whether one `β` summarises the shift at all six cutpoints. A Brant test or a partial-proportional-
  odds fit would answer it and would also be a second estimator on the primary outcome, which [§8] and
  invariant 5 forbid. The absolute-scale presentation is the honest mitigation and [§8] says so — the six
  `RD_k` show *where* the shift sits, and a reader comparing them against a single `exp(β)` is doing
  informally what a formal test would do. Recorded so that the absence is a decision.
- **`β`'s comparability across bootstrap replicates rests on that same assumption**, and Stage 10 is
  where it becomes a distribution rather than a parameter. `len(fit.alpha)` is a property of the
  replicate (§5.3), so a replicate missing a declared mRS level contributes a coefficient on a coarser
  scale. Under proportional odds those are the same parameter. §11 hands this to Stage 10 explicitly
  rather than leaving it to be discovered there, and nothing in this stage bounds how often it happens
  on the workbook's shape — that needs the replicate loop.
- **§6.3's bound is calibrated on a one-column design and Stage 12's is wider.** 4800 synthetic fits and
  a handful of separated constructions, all with treatment as the sole predictor, gave the empty band
  `(11.04, 18.98)`. A [§14a] standardisation model carries eight covariates, some continuous, and a
  legitimate coefficient on `age` in years is small while one on a rare factor level may not be. §11
  states that Stage 12 owes its own guard; what is *not* known is whether `POLR_MAX_ABS_BETA` is the
  right number there, and this stage does not find out.
- ~~**And it is calibrated at `K = 6` only**~~ — **MEASURED 2026-08-24, and the answer is that the band
  narrows but the bound holds** (§6.3, §21.1 item 24). Across every cutpoint count a replicate can
  produce, the legitimate maximum reaches **11.04** rather than 8.79 and the degenerate minimum falls to
  **18.98** rather than 18.81, so the safe interval is `(11.04, 18.98)` and 14.0 sits inside it with
  27% clearance below and 26% above. Every bound from 12 to 20 drops the same replicates. **What this
  retires**: the claim that the region between the modes is empty, and the value 10.0, which is now
  positively excluded. **What it leaves open**: nothing about `K` — but the sweep is still a sampler of
  this cohort's *shape* and not this cohort, so the rate at which [§10] actually meets each cutpoint
  count is still Stage 10's to report, which is why §11 asks for the distribution of `len(fit.alpha)`
  and `TODOS.md` keeps its item.
- **`polr`'s O6 rank check runs before §5.3's collapse, so the collapsed design's rank is unchecked.**
  Unreachable from either declared caller — overlap weights are strictly interior and Stage 12 passes
  unit weights, so the collapse drops no row — and reachable for **a caller that passes a weight of
  exactly zero**, which §5.3 explicitly permits. Such a caller can collapse away every row of a design
  column and get `np.linalg.solve`'s `LinAlgError`, which is not `FitError`, so [§10] would not catch it
  (§5.5). Not repaired, because repairing it puts the rank check in two places for a caller that does
  not exist; recorded so that §14.4's zero-weight assertion is not read as covering it.
- **The separation guard is measured to fire on constructions and has never fired on data.** Every
  number in §6 comes from synthetic frames. Whether v7 is anywhere near the bound is asserted in §14.12
  as a property — `|β| < POLR_MAX_ABS_BETA` — and that assertion is the first thing that will tell
  anyone, which is the right place for it and is not the same as knowing.
- **The threshold-model crossing [§8] warns about was not reproducible at a material magnitude.** 8 of
  1990 samples produced a strict decrease and the largest was below 5e-07 (§20). So [§8]'s reason for
  the empirical route is sound and its *urgency* is not measured here; the construction that would show
  a visible crossing is presumably sparser or more separated than the one tried, and finding it was
  stopped rather than completed. §14.6 asserts the difference between the two routes instead, which is
  falsifiable and does not depend on the rate.
- ~~**The R oracle is specified and cannot run on this machine.**~~ **STRUCK 2026-08-24 — it runs**
  (§18b, §21.1 item 23). `ordinal::clm` validates the weighted fit at 1e-12 on the coefficient and
  1e-9 on the cutpoints through `test_reference_r.py`'s `R_ENV`. What follows is the measurement that
  was mistaken for a blocker, kept because the failure it describes is real from a bare shell and the
  next person to run `Rscript` by hand will see it: measured 2026-08-21, **`Rscript --vanilla`
  segfaults during its own startup** — `system("uname -a")` returns status 139,
  `utils` and `stats` fail to load, and `requireNamespace` returns `FALSE` for every package including
  ones that are present. This is the PATH collision Stage 6 §16b.2 documented and Stage 7 §21.5 hit
  again, arriving a third time, and it is now the *only* thing blocking an independent oracle at two
  consecutive stages. `MASS` is absent from that library, so `polr` is not an alternative to `clm`
  there. §18b specifies the script and the gate; T6 is P2 for that reason.

## 16. What Stage 8 deliberately does not decide

- **Whether the ordinal fit should be penalised.** §6.4. [§8] names no penalty and [§7] names Firth for
  the propensity model, so the asymmetry is the SAP's. Penalising would reduce the replicates [§10]
  drops and would change the reported quantity to a penalised common odds ratio and the null
  distribution the primary p-value is read against. PI decision, [§8] amendment.
- **What `POLR_MAX_ABS_BETA` should be on a wider design.** §15. This stage declares the constant and
  measures the gap on its own design; Stage 12 owes its own calibration.
- **What to do with a replicate the guard rejects.** [§10] says drop and count. Whether a drop *rate*
  above some level invalidates the interval is [§10]'s question and the PI's, and this stage's
  contribution is to make the rate exist at all.
- **Whether the primary estimate should be reported at all if [§9]'s balance threshold is exceeded.**
  Stage 7 §13 and §14 record that the residual imbalance it measured is a PI decision under [§9] and
  [§13], and that *"it is not a licence to substitute an estimator"*. Stage 8 does not read the
  `Balance` and does not condition on it (§0.2). A pipeline in which a diagnostic silently selects an
  estimator is the mixture [§7] forbids.
- **Which of the six `RD_k` is the headline.** [§8] reports all six with intervals and no p-values, for
  the stated reason that threshold-wise tests invite selection of the most favourable cut. Nothing in
  this document knows which `k` carries the largest difference, and §4.3 is why that is deliberate.
- **How the estimate is presented beside the [§14] ones.** [§16]'s guardrail, Stage 14's to enforce.

## 17. NOT in scope for Stage 8

| Considered | Why deferred |
|---|---|
| A Firth or ridge penalty on the ordinal likelihood | §6.4, §16. [§8] names none; it is an amendment and not a repair, and it changes the reported quantity |
| Any standard error, covariance matrix or Wald interval | §5.6. The naive weighted-likelihood covariance treats the overlap weights as frequencies, and [§10] prescribes the interval that is used |
| A Brant test or any test of proportional odds | §15. It would be a second estimator on the primary outcome, and [§8]'s absolute-scale presentation is the mitigation it prescribes instead |
| A partial-proportional-odds or non-proportional-odds fit | §15, and [§15]'s list of approaches deliberately not used — the primary quantity is one number with one test |
| Augmenting the ordinal estimate | [§8] states there is no simple augmented form of a proportional-odds fit and that augmenting each threshold separately would replace one primary number with six. [§13] lists the augmented cumulative differences as a deferred sensitivity analysis, which Stage 11 owns |
| Threshold-specific outcome models for `RD_k` | §8.2. [§8] prescribes the weighted empirical distributions and gives the reason; the threshold route is [§13]'s deferred sensitivity row |
| Bootstrapping anything | Stage 10 [§10]. This stage returns point estimates and must be re-runnable inside a replicate loop, which §11 is the contract for |
| **Handing the fit to `scipy.optimize` and deleting `_ord_score_hess`** | Considered and declined, and recorded because it is the one alternative that would shrink this stage's highest-risk code. The Hessian's forty lines of chain rule are §21's own review-risk item 1, no standard error is produced from them (§5.6), and a gradient-only quasi-Newton needs only the score. Three reasons against, in order: `scipy` is **test-only by policy** (§14.3, Stage 6 §16b) and this would make it a runtime dependency of a shipped module; DoD-4's byte-identical reproduction across hash seeds would then rest on a third-party optimiser's iteration path, which is exactly the argument `data._md_table` gives for not using `tabulate`; and §5.2's two prespecified convergence routes with `converged_on` recording which fired are [§10]'s to read, and a `scipy` result object does not have them. The cost of declining is 40 lines and one review risk; §14.3's three-point derivative check is the mitigation and §18b's oracle 3 keeps `scipy` as the check it is good at |
| **Normalising `POLR_TOL` and `POLR_SCORE_TOL` by `Σw`** | §5.2. Both are absolute on quantities that scale linearly in `Σw`, which is the same objection §5.2 upholds against an absolute *trust radius* — and it is declined here because `FIRTH_TOL` and `FIRTH_SCORE_TOL` are absolute on identical grounds and are landed, so a normalised `POLR_*` pair would be the only tolerance in the pipeline with a different shape. Both declared callers are measured comfortable (`Σw = 27.736623` and `Σw = n`); a replicate's `Σw` is not, and §11 hands that to Stage 10. Changing the shape is a [§13] amendment, not a repair |
| A `covariates=`, `outcome=`, `thresholds=` or `bound=` parameter on `primary` | §10, §14.11. [§8] fixes all four, and a keyword makes a prespecified choice look like an option — Stage 6 §6.4's argument, two stages on |
| Reporting the drop rate or a failure counter | §6.4, §11. [§10] owns the replicate loop and therefore the counter; this stage owns the raise |
| A new `estimate` audit kind and heading | §9.1. `data.py:151-152` declares `model` the last amendment this pipeline needs and that Stages 7-13 render under it; Stage 7 tested that once and this stage tests it twice |
| Reading `Balance` to gate the estimate | §0.2, §16. Stage 7 §11 states Stage 8 does not, and §14.11's scan asserts it |
| Vectorising or caching the fit across replicates | §2. 1.157 ms per fit and about 2.3 s over `N_BOOT`; caching would make a replicate's estimate a function of what ran before it |
| Moving `MRS_THRESHOLDS` out of the `[§13]` subgroups block | §12. It is misfiled and the misfiling is harmless once `MRS_LEVELS[:-1] == MRS_THRESHOLDS` is asserted; moving a constant six test modules import is churn with a failure mode |
| Reporting the seven-category weighted distribution per arm | §9.3. [§8] asks for the cumulative form and `RD_k`; the per-category distribution is [§16]'s baseline material and Stage 14's |
| Pinning the workbook's estimate in this document or in the suite | §4.3, and §15 records the trigger under which it becomes correct to do |

## 18. What already exists, and what to lift

`pilots/analysis.py` has all three pieces this stage needs, and the split is sharper than at Stages 6
and 7: **the cumulative risk differences are right and are lifted whole; the ordinal fit is wrong in
three separate ways and none of it is lifted.**

| From `pilots` | Status |
|---|---|
| `cumulative_mrs_rd` (`:442-457`) — the weighted empirical route, and the docstring's reason ("threshold-specific outcome models can produce internally inconsistent cumulative predictions") | **lift both**, formula and comment. It is [§8]'s own route and [§8]'s own argument, and §8.2 verified the numbers agree |
| `_wmean`'s `if w.sum() else np.nan` (`:287-289`) | **lift the guard, and keep it before the mean.** Stage 7 §3.1: `np.average` with weights summing to zero RAISES rather than returning nan, so the order is load-bearing — the same lift Stage 7 made of the same idiom |
| `weighted_rd` (`:291-293`) — the difference of two weighted means | **lift as `weighted_proportion` differenced** (§8.1), made public and keyed by arm, because Stage 9 needs the two proportions and not only their difference |
| `WeightedOrderedModel` (`:272-284`) — subclass `OrderedModel`, override `loglike` to weight `loglikeobs` | **do not lift.** §2, §5.2. It inherits `OrderedModel`'s **numerical** derivatives, so every Newton step costs `2*(K+m)` extra likelihood evaluations and [§10] pays it 2000 times; and the class it inherits from has no weight support at all, which is why the override exists — a subclass whose parent would silently ignore the argument is one refactor from being replaced by the parent |
| `weighted_common_or`'s `except Exception: return np.nan` (`:471-472`) | **do not lift, and it is the most consequential line in this comparison.** A fit that fails returns a missing value rather than raising, so [§10] cannot drop and count the replicate — it silently receives a `nan` draw, and *dropping the nan draws* is selection on the outcome, which is exactly what [§10]'s p-value section forbids for a different quantity. [§8] prescribes one estimator and invariant 5 forbids a fallback; `nan` is a fallback |
| `weighted_common_or`'s `if y.nunique() < 3: return np.nan` (`:463`) | **do not lift.** Three is the wrong number — a proportional-odds fit needs **two** categories, not three — so a frame with exactly two silently returns `nan` where a perfectly well-defined fit exists. §5.5's O5 is the corrected form and it raises |
| `weighted_common_or`'s `float(np.exp(-res.params.iloc[0]))` (`:475`) | **do not lift, and do not lift the comment either.** §7.2. The negation is **correct** for `statsmodels`' `θ_k − x'β` and **wrong** for §5.1's `α_k + βA`, which is [§14a]'s own form. Lifting the line without lifting the parametrisation inverts the primary result of the study, with every number finite and the four cutpoints agreeing to five decimals |
| `weighted_common_or`'s `pd.Categorical(y, categories=sorted(y.unique()), ordered=True)` | **lift the concept, not the call.** Collapsing to the observed categories is right in spirit and §5.3 corrects it to the **positively-weighted** set, which is the set that makes the start values finite and the crossing argument true |
| `weighted_or`'s Haldane-Anscombe correction (`:296-312`) | **not this stage's.** It keeps a *marginal* odds ratio finite when a weighted proportion reaches 0 or 1, which is Stage 9's [§8] problem. `β` is a fitted model parameter and is defined whenever the model converges, which [§10] states explicitly |
| `sandwich_or` (`:478`) — an HC0 sandwich SE as a cross-check | **do not lift.** §5.6. It is a second interval on the primary quantity and it is the specific thing a reporting layer would find and print |
| `bootstrap_cumulative_mrs` (`:632-655`) | **not this stage's.** Stage 10 |

`pilots/test_estimators.py` has no test of the ordinal fit's sign, of its weighting, or of separation.
That is recorded rather than treated as a gap in the pilots: they are a prototype, and §21's amendment
1 exists because the roadmap inherited the prototype's assumption that a sparse ordinal fit fails
loudly.

## 18b. Validation against reference implementations

Stage 6 §16b established the standing rule this stage inherits: *use a maintained library wherever one
exists; where none does, the estimator carries an oracle rather than a comment, and the oracle is named
in the spec.* §2 measured that no library exists for a **weighted** proportional-odds fit in Python. So
this stage carries **five** oracles, and unlike Stage 6 the strongest of them is dependency-free.

**1. Integer weights equal row replication.** The strongest and the cheapest. A weighted likelihood
with integer weights is by definition the unweighted likelihood of the frame with each row repeated,
so the two fits must agree to machine precision — and no library, no subprocess and no tolerance
negotiation is involved. Measured (§20): **0.000e+00** on `β` and **1.332e-15** on `α` over a
120-record five-category frame with weights drawn from `{1, 2, 3, 4}` — which is
`replication_frame()` of §14.0.4 — and **0.000e+00** on `β` over the twelve-record hand frame,
`hand_ordinal()`, which is **not** §14.0.2's five-record cohort as three earlier drafts said (§21.1
item 2). **This is the oracle that validates the weighting**, which is the one
thing about this estimator that no unweighted comparison can check, and it is why §21's amendment 2
promotes it over the roadmap's all-weights-1 criterion.

**2. `statsmodels.OrderedModel`, unweighted.** A declared project dependency (§2), so no
`importorskip`. It validates the **unweighted** optimum and, more usefully, the sign convention: the
coefficients come back negated and the cutpoints do not (§7.2, measured 1.65e-05 and 3.16e-05). The
tolerance is `statsmodels`' `bfgs` and not ours, which is why oracle 1 carries the precision claim and
this one carries the convention claim.

**3. `scipy.optimize` on the same weighted objective.** Test-only by policy, Stage 6 §16b being the
precedent. It checks that the Newton loop reaches the maximiser of the function it says it is
maximising, which neither of the first two can: oracle 1 compares two of our own fits and oracle 2
compares a different objective. Measured (§20): Nelder-Mead from a deliberately poor start agrees to
**1.062e-07** on all six parameters with log-likelihoods identical to ten decimals. It is a check on
*ourselves against ourselves* in the sense Stage 6 §16b names — and its value is that it is the only
one that would catch a sign error in `_ord_score_hess` that happened to be self-consistent.

**4. The analytic derivatives against central differences.** Not usually called an oracle, and it is
one here: the score and Hessian are forty lines of chain rule (§5.2) and a wrong sign in the
`alpha`-`beta` cross block is the failure that would produce a slightly-wrong maximiser with a
plausible iteration count. Measured to **2.784e-08** and **2.812e-08**, with the Hessian exactly
symmetric and negative definite (§20).

**5. R's `ordinal::clm(..., weights = )` — SPECIFIED, AND IT RUNS.** `clm` takes observation weights
natively, which makes it the only **independent** implementation of the quantity this stage computes:
oracles 1, 3 and 4 check us against ourselves and oracle 2 checks an unweighted special case. It is
therefore the one that matters most, and an earlier draft of this section recorded it as unrunnable on
this machine. **That was measured wrongly and the correction is item 23.**

Measured 2026-08-24, through `test_reference_r.py`'s own `R_ENV` — Stage 7's `_r_environment()`, which
puts `/usr/bin:/bin` FIRST on `PATH`:

```
    _r_has(("ordinal",))  ->  True          _r_has(("MASS",))  ->  True
    ordinal::clm(ordered(y) ~ a, weights = w, link = "logit")  on hand_ordinal()

      ours.beta            -0.869564412763
      clm coefficient      +0.869564412758
      |ours + clm|          4.773e-12       the coefficients ARE NEGATED  [§7.2]
      |ours - clm|          1.739           so the sum test has something to negate

      |ours.alpha - clm thresholds|  max 7.177e-10   the cutpoints are NOT negated
```

**So the WEIGHTED fit is validated against an independent implementation at 1e-12 on the coefficient
and 1e-9 on the cutpoints**, and §7.2's asymmetry — coefficients negated, cutpoints not — is confirmed
by a second reference implementation rather than only by `statsmodels`. This is the strongest oracle
this stage has and it was recorded as unavailable.

**Why the earlier measurement said otherwise.** It ran `Rscript --vanilla` from a bare shell, where R
does segfault during startup: `system("uname -a")` returns status 139, `utils` and `stats` fail to load,
and `requireNamespace` then returns `FALSE` for packages that are present. All of that is real and
reproducible. What the draft missed is that **Stage 7 already repaired it** — `_r_environment()`
(`test_reference_r.py:89-130`) documents this exact PATH collision and fixes it — and the draft
explicitly wrote that Stage 7's repair "does not help", which is the opposite of what is measured.
Fourteen R oracle tests from Stages 6 and 7 pass in the suite today; the evidence was one `pytest -q`
away.

T6 therefore does three things and none of them is "get R working":

- Commits `tests/reference/polr_clm.R`, about 25 lines: read a CSV of `(y, a, w)`, call
  `clm(ordered(y) ~ a, weights = w, link = "logit")`, write the coefficient and the thresholds back at
  `%.17g`. No covariate list and no analysis logic, exactly as `firth_logistf.R`, `ato_psweight.R` and
  `balance_psweight.R` are shaped.
- Adds the Stage 8 half of `test_reference_r.py` behind the gate that already exists, comparing `clm`'s
  coefficient against **`-ours`** and its thresholds against **ours** — the convention §7.2 measured,
  asserted rather than assumed, because `clm` uses `zeta_k − x'β` like `polr` and `OrderedModel`.
- **Makes the skip name the segfault.** The gate currently reports "missing", which on this machine is
  false. The skip reason must distinguish *package absent* from *R does not start*, because Stage 6
  §16b.2's whole argument is that a silently-skipping oracle is worse than no oracle, and a skip
  reason that misattributes the cause sends the next reader to install a package that is already
  installed. Stage 7 §21.5 predicted a third round would find nothing new here and would only need to
  run the oracle; what a third round found instead is that it still cannot.

**Not a completion blocker**, for Stage 7 T6's reason and one more: four oracles do run, one of them at
machine precision, and the quantity they check is a maximum-likelihood fit rather than a penalised one
whose kernel nothing else implements. It is P2, filed with the segfault as the blocker rather than as a
mystery.

## 19. Implementation tasks

Ordered. **T1 is the amendment that unblocks everything, T2-T3 are `model.py`'s and are independently
verifiable against their own privates, T4-T6 are `outcome.py`'s, T7 is the R oracle and T8 is the log
and the definition of done.** Every acceptance section that reads a `Primary` or an `Audit` belongs to
T6 — Stage 7 §17 records what it cost to assign those to earlier tasks that could not run them.

- [ ] **T1 (P1)** — `config.py`: `PRIMARY_OUTCOME`, `MRS_LEVELS` and the `POLR_*` block, as §12 writes
      them. `test_config.py`: the six assertions of §12. Verify: `uv run pytest -v` green with **no
      other test edited**, which is where §12's claim that `data.py` needs no amendment is checked
      rather than believed. `PRIMARY_OUTCOME` uses `next(...)` over a generator and there is exactly
      one primary entry, which the second assertion pins — if a second were ever added, `next` would
      silently take the first and only the assertion would notice.
- [ ] **T2 (P1)** — `model.py`: `PolrFit`, `_weighted_categories`, `_ord_pieces`, `_ord_loglik`,
      `_ord_score_hess`. `test_model.py`: **§14.0.4's six generators, §14.0.4a's six, and
      `BAND_EFFECTS`, written before anything that asserts on them** — every later task's numbers have these as their subject,
      and an earlier draft of this document specified five of the six in prose only, which is the
      defect §21.1 item 1 records. Verify: **direct assertions on the privates** — acceptance 14.3's
      derivative checks, **at three parameter values on two frames and not at the optimum alone**, and
      14.4's `-inf` branch. **Two things to get exactly right**, each of which is a wrong sign away from
      a plausible number: the boundary conventions, where an observation in the lowest category has no
      lower cutpoint and one in the highest has no upper one, so `gu` defaults to 1.0 and `gl` to 0.0
      and **not** to the clipped expit of a clipped index; and the `beta` block of the Hessian, which
      collects `Huu + 2*Hul + Hll` — both single-sided terms **plus twice** the cross term. Assert
      against central differences before writing anything that calls these.
- [ ] **T3 (P1)** — `model.py`: `_assert_polr_fittable` with O1-O6 in §5.5's order, and `polr`. Verify:
      acceptance 14.3, 14.4, 14.7's fit half, 14.8. **O3 must precede O5** or a `nan` weight deletes an
      outcome category instead of raising (§5.5, measured). **The collapse must be to the
      positively-weighted set** and not to `np.unique(y)`, or §5.4's ordering guarantee and §5.2's
      finite start values both become nearly-true — and it must go through `_weighted_categories`, from
      which O5 also reads, because the rule is one function and not one expression written twice
      (§5.3a). And `polr` must **not** prepend an intercept — the opposite of `firth` in the same
      module, which is the one place a reader will assume they agree. **§14.7's assertions are probes
      and not calibrations**: 240 fits and not 4800, 20 healthy `cond(-H)` fits and not 500, 40 sparse
      replicates and not 400. If this task's tests take more than a second, they are re-deriving §20
      (§21.1 item 3).
- [ ] **T4 (P1)** — `outcome.py`: `weighted_proportion`, `cumulative_rd`. Verify: acceptance 14.6.
      Independent of T2 and T3 — neither function fits anything — and the one place in this stage where
      two worktrees would not cost more than they save. **The zero-weight guard comes before the mean**
      (§8.1), and `cumulative_rd` takes no threshold list.
- [ ] **T5 (P1)** — `outcome.py`: `Primary`, `_assert_primary_inputs` with G1-G5 **in §4.4's two
      phases**, `_assert_exposure_survived`, `_assert_reportable`, `_orientation`. Verify: **direct
      assertions on each** — acceptance 14.5's `_orientation` half, 14.7's guard half, 14.14. **The
      phase-1 raise is load-bearing, not stylistic**: phase 2's first line builds the mask from
      `ps.in_model` and `df[C.PRIMARY_OUTCOME]`, so an absent outcome column raises `KeyError` there
      and G3's message is discarded (§4.4a). Stage 7 §21.1 found exactly this defect *inside its own
      first repair*; this stage writes the boundary correctly from the start and §14.14 carries the
      companion that proves it.
- [ ] **T6 (P1)** — `outcome.py`: §9's three entries with their detail and table helpers, and `primary`
      as written in §10. `test_outcome.py`: `truncated_cohort()` (§14.0.5) alongside `ordinal_cohort()`.
      `test_outcome.py` also gains §14.0.4a's `separated_cohort()`, which is what gives §14.7's
      "`primary` on the same frame raises G7" a subject — `separated_frame()`'s three arrays are not a
      cohort and a `Propensity`. Verify: acceptance **14.1, 14.2, 14.9, 14.10, 14.11, 14.12** —
      everything that reads a `Primary` or the `Audit`. Watch four things: the mask is bound **once** (§4.2); entry 1 is recorded
      **before** the design so a G6 or G7 failure still leaves the [§11] denominator in the log
      (§14.14); every `detail` branch must be **read as prose against rendered output**, not only
      asserted as a count — Stage 7 §21.5a is the precedent and it cost a review round to find; and
      **`_fit_table`'s no-cutpoint branch keys on `fit.categories[-1]` and never on `MRS_LEVELS[-1]`**
      (§9.2). That last one is the defect this stage's review found in its own draft: the two are the
      same on every frame the suite had, they come apart in an ordinary [§10] replicate that draws no
      death, and the wrong version renders `missing` under status `fitted` with nothing raising
      (§21.1 item 6). `truncated_cohort()` exists for it.
- [ ] **T7 (P1)** — the R oracle: `tests/reference/polr_clm.R` and the Stage 8 half of
      `test_reference_r.py`, behind the gate Stage 7 built. **P1 and not P2, because it runs** (§18b,
      §21.1 item 23): it is the only INDEPENDENT check on the weighted fit, and every other oracle
      compares us against ourselves or against an unweighted special case. Verify: `clm`'s coefficient
      matches **`-ours`** to 1e-9 and its thresholds match **ours** to 1e-8 — measured 4.773e-12 and
      7.177e-10 on `hand_ordinal()` — and the gate is asserted **measured open** rather than assumed,
      as `test_reference_r.py:388` already does for Stage 6's. Keep the skip path working for a machine
      without `ordinal`, and make its reason distinguish a missing package from an R that does not
      start, because the bare-shell segfault of §15 is real and the next reader will meet it.
- [ ] **T8 (P1)** — the audit log produced through all eight stages against the workbook, read end to
      end by a human, and the definition of done below. **This is the task at which the primary
      estimate is computed for the first time** (§4.3), and it is deliberately last.

### What can be built in parallel, and what cannot

| Task | Modules touched | Depends on |
|---|---|---|
| T1 | `config.py`, `test_config.py` | — |
| T2 | `model.py` (the ordinal privates), `test_model.py` (§14.0.4's generators) | T1 |
| T3 | `model.py` (`polr`, O1-O6), `test_model.py` | T2 |
| T4 | `outcome.py` (the two weighted functions), `test_outcome.py` | T1 |
| T5 | `outcome.py` (the types and the five guards), `test_outcome.py` | T1 |
| T6 | `outcome.py` (`primary`, the entries), `test_outcome.py` | T3, T4, T5 |
| T7 | `test_reference_r.py`, `tests/reference/polr_clm.R` | T3 |
| T8 | none — it runs the pipeline and reads the log | T1-T6 |

```
  Lane A:  T1  →  T2  →  T3  ┬→  T6  →  T8
  Lane B:      T4  ──────────┤
  Lane C:      T5  ──────────┘
  Lane D:              T7 ───────────────┘   (joins only at the suite)
```

**T4 and T5 both touch `outcome.py` and are genuinely independent**: `weighted_proportion` takes three
arrays and knows nothing about a frame, and the guards take a frame and a `Propensity` and compute
nothing. Merge T4 first — T6's `primary` calls both, and a conflict there means the mask was bound
twice.

**The one dependency the table does not show is `test_model.py`.** T2 writes §14.0.4's generators and
T4, T5 and T6 all import them, so **T2's generator block must merge before Lane B or Lane C writes a
test that uses one**. It is a dozen lines and no logic; write it first, even before T2's privates.

### Diagrams and comments that belong in the code, not only here

Five, and no more — a diagram nobody maintains is worse than none, because it is believed.

- **`outcome.py`'s module docstring** — §6's sentence: separation in this estimator does not present as
  non-convergence, it converges on the score criterion with clean counters and returns an odds ratio
  of 1e15, so G7 is the only thing that turns a degenerate fit into a failure [§10] can count. That is
  the sentence that stops the next reader deleting the bound as a magic number.
- **`model.py`'s `polr` docstring** — §5.1's parametrisation, spelled out, with the note that every
  reference implementation uses `α_k − x'β` and that the coefficients differ in sign while the
  cutpoints do not. That is the sentence that stops the next reader "fixing" the sign to match a
  textbook.
- **Above `polr`'s collapse line** — §5.3's two reasons in two lines: positively-weighted and not
  merely observed, because that is what makes the start values finite and the crossing unreachable.
- **Above `_assert_polr_fittable`'s O3** — §5.5's ordering: O3 precedes O5 because `nan > 0.0` is
  False, so without it a `nan` weight deletes an outcome category instead of raising, measured.
- **Above `_fit_table`'s status branch** — §9.2's correction in two lines: the no-cutpoint test is
  `level == fit.categories[-1]` and never `level == MRS_LEVELS[-1]`, because the two are the same on
  every frame this suite has and come apart in any replicate that draws no death. That is the sentence
  that stops the next reader simplifying it back, and it is the only one of the five that guards against
  a defect this document actually contained.

### Definition of done

Stage 8 is complete when all of the following hold, and not before.

1. `uv run pytest -v` is green with `data/` present.
2. Green **with `data/` temporarily renamed** — the data-gated tests skip and nothing else fails or
   errors at collection.
3. `test_config.py`, `test_data.py`, `test_derive.py`, `test_eligibility.py`, `test_cohort.py`,
   `test_model.py`, `test_propensity.py` and `test_balance.py` are still green, **unedited except for
   `test_config.py`'s five new assertions and `test_model.py`'s new sections**, and neither
   `outcome.py` nor `test_outcome.py` is in `EXEMPT_FROM_RAW_NAME_SCAN`.
4. The pipeline's audit log has been produced through all eight stages against the **workbook** and
   reproduces byte-identically across hash seeds:

   ```bash
   cd extended_bridging
   S='import sys, data, derive, eligibility, cohort, propensity, balance, outcome
   df, audit = data.load(data.WORKBOOK)
   df = derive.derive(df, audit)
   df = eligibility.classify(df, audit)
   df = cohort.build(df, audit)
   ps = propensity.fit(df, audit)
   balance.assess(df, ps, audit)
   outcome.primary(df, ps, audit)
   sys.stdout.write(audit.to_markdown())'
   PYTHONHASHSEED=0 uv run python -c "$S" > /tmp/a.md
   PYTHONHASHSEED=1 uv run python -c "$S" > /tmp/b.md
   diff /tmp/a.md /tmp/b.md
   ```

   `diff` reports nothing, the three new entries render under `## Fitted models` after
   `overlap_by_centre`, and the ledger reads 27. **The log carries the primary estimate and names
   patients, and must never be committed.**
5. **Separation has been seen to converge.** With `_assert_reportable` removed and a perfectly
   separated frame, `primary` is shown to **return** an `odds_ratio` above 1e15 with
   `converged_on == "score"`, `iterations == 17`, `rescales == 0`, `halvings == 0` and every fitted
   probability finite. Then, with the guard present, the same input is shown to raise `FitError` naming
   G7 and the rejected odds ratio. **Watch this one specifically: it is the only finding in this stage
   that was found by running the code rather than by writing it, and it is the one the roadmap's
   Accept when got wrong.**
6. **The sign has been seen to invert.** With §5.1's `+ beta*A` replaced by `- beta*A`, the orientation
   construction of §14.5 is shown to report `exp(beta) = 0.0118` where it reported 84.58, with `alpha`
   unchanged to five decimals and nothing raising. Then `pilots/analysis.py:475`'s `np.exp(-...)` is
   applied to our `beta` and shown to do the same thing. A reader who has not seen the cutpoints agree
   while the coefficient flips has not seen why §7.2 is long.
7. **The weighting has been seen to be read.** Integer weights are shown to equal row replication to
   0.000e+00 on `beta`; then `w` is ignored entirely inside `polr` and the all-weights-1 test is shown
   to **still pass** while the replication test fails. That pair is §21's amendment 2 in executable
   form.
8. **The category deletion has been seen.** With O3 removed, one `nan` weight is shown to take a
   four-category fit to three and to move `beta` from −1.775143 to −1.205433, with `iterations == 4`,
   `converged_on == "likelihood"` and nothing raising. `np.nan > 0.0` is confirmed `False` at the
   prompt.
9. **The vanished exposure has been seen.** With G6 removed and a one-armed frame,
   `model.design(sub, (TREATMENT,))` is shown to return 0 columns with `dropped == ('ivt',)`, `polr` on
   it is shown to **converge**, and `primary` is shown to give either an `IndexError` on `beta[0]` or an
   odds ratio of exactly 1.0 depending on how the result is read. Neither is a failure the caller can
   see.
10. **The [§11] denominator has been seen to be per estimate.** With one `in_model` record's `mrs_90d`
    blanked on `ordinal_cohort()`, `in_estimate` is shown to be 4 against `in_model`'s unchanged 5,
    entry 1's `n` to be 1 and the record to be named; and blanking a `PS_COVARIATE` instead is shown to
    move both together. On the workbook the column is a constant, so this is the only place either form
    can be told from the other.
11. **The collapse has been seen on a real frame.** `ordinal_cohort()` is shown to fit 4 cutpoints over
    `(0, 1, 3, 4, 6)` with levels 2 and 5 marked absent in entry 2's table, and the workbook is shown
    to fit 6 over all seven. §14.0.2's golden vector reproduces to every digit printed.
12. **The top FITTED level has been seen not to be the top DECLARED one.** On `truncated_cohort()`
    (§14.0.5), `fit.categories` is `(0, 1, 2, 3, 4)`, and entry 2's table is **read by eye**: `mRS <= 4`
    reads `missing` under the highest-**fitted** status, `mRS <= 5` and `mRS <= 6` read `missing` under
    the absence status with 6 also saying it is the highest declared level, and no row reads `missing`
    under status `fitted`. Then §9.2's branch is reverted to `level == MRS_LEVELS[-1]` and the same
    table is shown to render `alpha at mRS <= 4 | 4 | missing | fitted` with nothing raising and the
    whole suite green except §14.9. **Watch this one: it is the defect the review round found in this
    document's own draft, and it is invisible on the workbook and on `ordinal_cohort()` alike.**
13. **The acceptance suite has been timed.** `uv run pytest --durations=15` is read, and no Stage 8
    test is above **1 s**. §14.7's assertions are probes at 240 / 20 / 40 fits; if one of them is
    re-deriving §6.3's 4800-fit calibration or §20's 500 healthy `cond(-H)` fits, it shows here as a
    multi-second test and the fix is to move the calibration to §20 rather than to raise the threshold
    (§21.1 item 3).
14. **Both `detail` strings have been read as prose against rendered output, under every branch they
    have** — entry 2's absent-levels clause at 0 and at 2, entry 3's per-arm totals clause — and each
    reads as English with its counts matching the table beside it. Stage 7 §21.5a found an unterminated
    clause in exactly this position that three fence assemblies and two review rounds walked past;
    §14.9 asserts the sentences, and this item is the human half.
15. **Both tables render as well-formed markdown.** For each of the three entries the rendered header,
    separator and body rows have the same number of markdown cells, and no cell contains a `|`. Look at
    the log, not only at its hash: a table identically broken under both hash seeds is still identical
    (Stage 7 §21.1).
16. Run from `extended_bridging/`, `grep -nE "^\s*(import|from)\s+.*\bbalance\b" outcome.py` returns
    nothing **and exits 1**, and `test_outcome.py`'s AST scan for the same is green with its companion
    firing on a pasted import. The path is relative to `extended_bridging/` because DoD-4 has already
    `cd`-ed there, and a grep that exits 2 on a missing file "returns nothing" and reads as a pass
    (Stage 7 §21.4).
17. **No float literal outside §14.12's declared allowlist appears in §14.12's test section**, checked
    by the meta-assertion §14.12 specifies — banner-delimited, both delimiters asserted found, the range
    asserted non-empty — and confirmed by eye. Its companion is shown to fire on a pasted magnitude
    **and** on one pasted into a docstring. This is what keeps §4.3 true as the file is edited.
18. **Every §14 assertion has a subject that exists in the repository.** Each of §14.0.4's six
    generators and §14.0.4a's seven is called from at least one test, and no §14 assertion names a
    construction that is only described in prose. **Checked by enumeration and not by impression**:
    walk §14 bullet by bullet, and for each one that names a frame, name the generator it comes from.
    An earlier draft of this item asserted the property was already true after §14.0.4 repaired six
    constructions; **it was not — eight more assertions still named prose-only frames**, which is
    §21.1 item 22 and is why §14.0.4a exists. The lesson is the item's own: a completeness claim
    written in the same pass as the repair is a claim about the repair's intent, not its extent.
19. The audit log has been read end to end by a human, and its numbers agree with §20: entry 2's `n`
    equals `propensity_fit`'s `n` on v7, entry 3's per-arm totals equal `overlap_weights`' per-arm
    `sum w`, and the ledger reads 27.
20. **The estimate has been produced, read, and handed to the PI** — and §15's `TODOS.md` item is filed
    at the same time, so that the pins §4.3 declined become available the moment the analysis is locked.

## 20. Verification record

Checked on 2026-08-21, before this spec was finalised, by a read-only probe and four synthetic
harnesses. The probe drives `load → derive → classify → build → fit` against `data.WORKBOOK`, reports
aggregate quantities only, and writes nothing. **The probe does not compute `β`, `exp(β)`, `RD_k`, or
any mRS distribution split by arm**, which is §4.3's constraint expressed as a property of the
instrument rather than as an intention; every fit measured below is on synthetic or fixture-derived
data. Acceptance 14.12 re-checks the workbook's properties in code.

**This table is normative and every number elsewhere in this document points at it.** Stage 7 §18
earned that rule the hard way — a wrong design-matrix width propagated into five places — and the
argument is the same here and stronger, because this document has more measured claims than that one.
A number in prose cites the row here that produced it, this table cites the code that produced it, and
where the two disagree **this table is right and the prose is stale**. An implementer who finds a
conflict should fix the prose and not re-derive the number.

| Claim | Where used | Verified |
|---|---|---|
| **`statsmodels.miscmodels.ordinal_model.OrderedModel` has no weight support**: its signature is `__init__(self, endog, exog, offset=None, distr='probit', **kwds)` — no `weights`, `freq_weights` or `var_weights`. statsmodels **0.14.6** | §2, §18, §18b | yes, run |
| The three populations on the workbook: cohort **93**, `in_model` **92**, `mrs_90d` present on **93**, `in_estimate` **92**. The outcome mask removes **0** records, and `Σw` is unchanged to every digit. Arms inside `in_estimate`: 53 EVT alone, 39 bridging | §4.1, §4.2, §13, §14.12 | yes, run |
| **`Σw` over the estimation population is 27.736623**, identical over `in_model` to every digit (the mask removes nobody), with individual weights in **[0.0269724, 0.964955]** — strictly interior, which is Stage 6's F5 confirmed one stage on. This is the scale at which `POLR_TOL = 1e-8` and `POLR_SCORE_TOL = 1e-6` are absolute tolerances, and it is measured because §5.2 argues from it. An aggregate weight sum with no arm split and no outcome in it, so §4.3 is untouched | §5.2, §11, §17 | yes, run |
| All **7** declared mRS levels are occupied in the workbook's estimation population, so `K` is **6** and §5.3's collapse is inert there. Marginal count only — no arm split | §5.3, §14.1, §14.12 | yes, run |
| `MRS_LEVELS` is `(0,…,6)` from `PLAUSIBLE_RANGES["mrs_90d"]` and `MRS_LEVELS[:-1] == MRS_THRESHOLDS` is **True** | §8.3, §12 | yes, run |
| `design(cohort, (TREATMENT,))` gives **1 column** `('ivt',)` with nothing dropped; on a one-armed frame it gives **0 columns** with `dropped == ('ivt',)`, and `polr` on that width-0 design **converges and returns a `beta` of length 0**. Neither function raises | §3.1, §4.5, §10, §14.11, DoD-9 | yes, run |
| **The analytic score and Hessian agree with central differences to 2.784e-08 and 2.812e-08**; the Hessian is symmetric to **0.000e+00** and its eigenvalues are all negative (max −7.799) — §5.2's concavity claim, measured | §5.2, §14.3, §18b | yes, run |
| **Integer weights equal row replication**: max \|Δbeta\| **0.000e+00** and max \|Δalpha\| **1.332e-15** on a 120-record five-category frame (n replicated 302 against Σw 302); **0.000e+00** on beta on §14.0.2's frame (n replicated 23) | §14.3, §18b, §21 | yes, run |
| **`polr(X, y, np.ones(n))` equals `polr(X, y, None)` bit for bit**, so an implementation that never reads `w` satisfies the roadmap's all-weights-1 criterion; and the weighted fit differs from the unweighted by **1.8268** in beta on `hand_ordinal()`, with opposite signs | §14.3, §21 | yes, run |
| **`scipy.optimize` (Nelder-Mead) on the same weighted objective agrees to 1.062e-07** on all six parameters, with log-likelihoods identical at −207.4848223118 to ten decimals | §14.3, §18b | yes, run |
| **`statsmodels` unweighted: the coefficients are NEGATED and the cutpoints are NOT.** `ours.beta + sm.beta` max **1.65e-05**; `ours.alpha − sm.thresh` max **3.16e-05**. ours.beta `[0.97283069, 0.58382310]` against sm `[-0.97284722, -0.58383309]`; ours.alpha `[-1.54340634, -0.52258248, 0.37686850, 1.76824648]` against sm `[-1.54343790, -0.52259457, 0.37684297, 1.76826674]`. The tolerance is statsmodels' bfgs, not ours | §7.2, §14.5, §18b | yes, run |
| **A PERFECTLY SEPARATED FIT CONVERGES AND RETURNS**: 20 v 20, two categories — beta **36.4058**, exp(beta) **6.469e+15**, alpha **[−18.2029]**, **17** iterations, `converged_on` **"score"**, rescales **0**, halvings **0**, min and max fitted category probability both **1.000000**, max \|linear predictor\| **18.2029** against POLR_ETA_CLIP 500, max \|score\| **2.487e-07** against POLR_SCORE_TOL 1e-06, **cond(−H) 6.854** on a 2×2 (this construction has two outcome categories), eigenvalues [−6.5e-07, −9.0e-08] | §3.1, §6.1, §14.7, DoD-5 | yes, run |
| It scales with n rather than being an artefact of a small frame: n per arm 10 / 20 / 40 / 80 gives beta **34.41 / 36.41 / 38.41 / 38.41** and exp(beta) 8.76e+14 / 6.47e+15 / 4.78e+16 / 4.78e+16, in 16 / 17 / 18 / 18 iterations | §6.1 | yes, run |
| **NEAR separation lands in the same place with clean counters**: 1 / 2 / 3 crossover patients out of 20 give beta **21.1512 / 20.4154 / 19.9711**, exp(beta) 1.534e+09 / 7.35e+08 / 4.713e+08, and min fitted probability **0.05 / 0.10 / 0.15** — so a detector keyed on the fitted probability catches perfect separation and MISSES the case a bootstrap draws | §6.1, §6.3, §14.7 | yes, run |
| **400 sparse replicates of a 92-record 39/53 seven-category frame produced ZERO convergence failures**, because separation converges rather than failing | §3.1, §6.1, §11, §14.7 | yes, run |
| **500 well-behaved 92-record seven-category fits at a true beta of 0.5** give \|beta\| min **0.0006**, median **0.5162**, max **2.1316**, and fitted category probabilities min-of-mins 0.0062 / max-of-maxes 0.4726. **This measurement is superseded as a basis for the bound**: a single true effect size is an assumption about the answer, which §4.3 forbids | §6.3 | yes, run |
| **THE EMPTY BAND. 4800 fits** of a 92-record seven-category frame over **twelve** true effect sizes from 0 to 6, 400 samples each: **max \|beta\| below 10 is 8.7873, min above 10 is 18.8055, and NOTHING lies in (8.79, 18.81)**. Histogram: [0,2) 2723, [2,4) 1024, [4,6) 475, [6,7) 234, [7,8) 34, [8,9) 12, [9,18) **0**, [18,20) 279, [20,25) 17. The distribution is bimodal — a fit either identifies a moderate effect or separates. `POLR_MAX_ABS_BETA = 14.0` sits 59% above the first endpoint and 25% below the second; exp(14) = 1.2e6. **An earlier draft set it to 10.0 from the single-effect-size measurement above, which put it 14% above that measurement's maximum and at the band's lower edge** | §6.3, §12, §14.7, §21 | yes, run |
| **`cond(-H)` DOES discriminate at equal cutpoint count, and the earlier claim that it did not was a size artefact.** On seven-category 92-record frames: healthy **29.1 to 162.8** (500 fits), near-separated (\|beta\| ≈ 20) **7.9e7**, exactly separated **3.8e8**. §6.1's two-category construction reads **6.854**, which is a 2×2's and sits *below* the entire healthy seven-category range — so the earlier "6.85 against 3 to 40, does not discriminate" compared matrices of different sizes and both of its numbers were wrong. It is rejected as the guard for size dependence and for being unjustifiable against a reportable quantity, not for failing to discriminate | §6.3, §14.7, §21 | yes, run |
| **One `nan` weight DELETES an outcome category**: a four-category fit became three, `categories` (0,1,2,3) → (0,1,2), beta **−1.775143 → −1.205433**, iterations 4, `converged_on` "likelihood", nothing raised. `np.nan > 0.0` is **False** | §3.1, §5.5, §14.8, DoD-8 | yes, run |
| **One `nan` response DROPS a record**: 80 records fitted as **79** with `categories` unchanged, nothing missing anywhere in the output | §3.1, §5.5, §14.8 | yes, run |
| **An absent middle category collapses the fit**: a five-category frame with category 2 emptied fits **3** cutpoints against the full frame's 4, alpha ascending and start values finite in both. And a category present but carrying **zero total weight** gives a fit identical to **1e-10** to the fit on the frame with those rows removed | §5.3, §14.4 | yes, run |
| **A cutpoint crossing was NOT reachable**: 3000 attempts with pure Newton from a start spanning 2e-3 never lost the ordering, so `_ord_loglik`'s −inf branch does not fire on any input this stage has. An earlier draft claimed it was reachable by removing the halving; the claim is struck | §5.2, §5.4, §14.4 | yes, run |
| **The orientation**: on a 200-record construction with treated-arm mean category 1.11 against 3.82, beta **4.437745** and exp(beta) **84.584014**, so > 1; the `alpha_k − beta*A` reading of the same fit gives **0.011823** | §7.3, §14.5, DoD-6 | yes, run |
| **`RD_6` is structurally zero and `RD_5` is not**: on a construction with deaths in both arms (18 and 12), RD_5 = **+0.049027** and RD_6 = **+0.000000** exactly, with P(Y≤6) = 1.000000 in both arms | §8.3, §14.6 | yes, run |
| **The threshold-model route crosses rarely and marginally, and an earlier draft's 40-of-40 was a bug.** Over 2000 samples of a 40-record seven-category frame, 10 had a failing threshold fit and **8 of the remaining 1990 (0.40%)** showed a strict decrease > 1e-9, with the largest decrease **below 5e-07** — float noise, not a visible crossing. The empirical route was non-monotone in **0 of 2000**. The earlier 40-of-40 measurement computed `P(ind ≤ 0)`, which is `1 − P(Y ≤ k)`, so every sequence decreased by construction | §8.2, §14.6, §15, §21 | yes, run twice |
| **The committed fixture cohort cannot reach this stage**: `cohort_frame()`'s five-record cohort has `mrs_90d` **[2, 2, 2, 2, 2]** — one category — and `polr` raises O5. `cohort_frame()` builds its cohort records from the `HAND-4` template, so every non-overridden column is constant | §14.0, §14.0.1, §5.5 | yes, run |
| **`ordinal_cohort()` fixes it and exercises the collapse**: the five survivors are HAND-1, HAND-2, HAND-5, COHORT-1, COHORT-3 with `mrs_90d` **[0, 1, 4, 6, 3]**, `in_model` 5 of 5, `categories` **(0, 1, 3, 4, 6)** — levels 2 and 5 absent — and **4** cutpoints. The derived dichotomies come back recomputed: `mrs_0_2_90d` [1,1,0,0,0], `mrs_0_1_90d` [1,1,0,0,0], `death_90d` [0,0,0,1,0], `mrs_5_6_90d` [0,0,0,1,0] | §14.0, §14.0.1, §14.0.2 | yes, run |
| **§14.0.2's golden vector** on `ordinal_cohort()`: beta **1.7516026138**, exp(beta) **5.7638324754**, alpha **[−4.0128444951, −1.8393869358, −0.7580715177, 0.4483367817]**, `cut_levels` (0,1,3,4), **4** iterations, `converged_on` "likelihood", rescales 0, halvings 0, alpha ascending, \|beta\| inside the bound. RD_k: 0 → +0.1079772317, 1 → +0.5469186443, 2 → +0.5469186443, 3 → +0.0181363724, 4 → +0.4712177281, 5 → +0.4712177281, with P(Y≤k) control [0, 0, 0, 0.5287822719, 0.5287822719, 0.5287822719] and bridging [0.1079772317, 0.5469186443, 0.5469186443, 0.5469186443, 1, 1] | §14.0.2, §14.6 | yes, run |
| A second hand-built twelve-record frame (6 v 6, all seven levels) gives beta **3.8970544631**, exp **49.2571464698**, 6 cutpoints, 5 iterations; unweighted beta **3.9927437960**; integer-weight replication beta **4.08833138004692** both ways. **This is the §20 harness's frame and is NOT `hand_ordinal()`** — an earlier draft identified the two and they are different frames with different weights (§21.1 item 12). It is also not §14.0.2's, whose cohort has five records (§21.1 item 2). `hand_ordinal()`'s own measured numbers are the row below | §21.1 | yes, run |
| **`hand_ordinal()` of §14.0.4, measured**: weights `[1,2,3,4,5,6,6,5,4,3,2,1]`, Σw **21 in each arm**. Weighted beta **−0.8695644128**, exp **0.4191340789**, **6** cutpoints, **4** iterations on "likelihood", rescales 0, halvings 0; alpha `[−3.3804054388, −0.9361654608, 0.0247795185, 0.8447848942, 1.8057298736, 4.2499698516]`. Unweighted beta **0.9572618608**, so weighted and unweighted differ by **1.8268** and have **opposite signs** — the weights ascend across the treated arm's worsening outcomes and descend across the control's. Integer-weight replication agrees at **1.11e-16** on beta and **8.88e-16** on alpha. RD_k: 0 → **+0.0476190476**, 1 → **−0.1428571429**, 2 → **−0.2380952381**, 3 → **−0.2380952381**, 4 → **−0.1428571429**, 5 → **+0.0476190476** | §14.3, §14.6, §18b, §20c | yes, run |
| **The numbers in this table came from the probe harnesses, not from §14.0.4's generators.** The generators were written afterwards, to make each construction reproducible, and their sampling differs in detail from the harnesses' — so every §14 assertion that ranges over a seeded generator is a **tolerance or a property** and not one of these values. The two seed-free constructions, `hand_ordinal()` and `separated_frame()`, are the exceptions and §14.6 and §14.7 pin them as values | §14.0.4, §14.3, §14.5, §14.7, §21.1 | stated |
| **Timing**: one 92-record seven-category treatment-only weighted fit **1.157 ms**; `cumulative_rd` over six thresholds **0.067 ms**; the fit alone over `N_BOOT` = 2000 is about **2.31 s**. Stage 6's Firth fit is tens of microseconds | §2, §11 | yes, run |
| pandas **2.3.3**, numpy **1.26.4**, statsmodels **0.14.6**, Python **3.12.12** in the project environment | §3.1, throughout | yes, read |
| Audit entries through the seven landed stages: **22** after `fit` (load 7 / derive 4 / classify 1 / build 6 / fit 4), confirming Stage 6 §9's and Stage 7 §9's ledgers; Stage 7 adds 2 and Stage 8 adds 3, for **27** | §11, §13, §14.0.3, §14.9 | yes, run |

**Claims about the R oracle**, run 2026-08-21:

| Claim | Where used | Verified |
|---|---|---|
| `ordinal` **2023.12-4.1** is present in `~/R/library` and takes `weights` on `clm`. **`MASS` is present too** — `_r_has(("MASS",))` is True — so an earlier claim that it was absent is struck (§21.1 item 23) | §18b, §15 | yes, run |
| **`Rscript` segfaults during its own startup on this machine**: `Segmentation fault (core dumped)`, `system("uname -a")` status **139**, `utils` and `stats` not found, and `requireNamespace` returns `FALSE` for every package **including ones that are present**. `R.version.string` prints (R 4.5.0), so the binary runs. Identical with `R_LIBS_USER` set to the directory the packages are in | §15, §18b, T7 | yes, run |
| **A `FALSE` from a BARE-SHELL probe is not evidence of absence here** — but `test_reference_r.py`'s `R_ENV` is not a bare shell, and through it every probe returns True. An earlier draft concluded that Stage 7's repair "does not help"; measured 2026-08-24, **it helps completely** and the fourteen R oracle tests of Stages 6 and 7 pass in the suite today (§21.1 item 23) | §15, §18b, §21.1 | yes, run |
| **`ordinal::clm` REPRODUCES OUR WEIGHTED FIT.** On `hand_ordinal()` through `R_ENV`: clm's coefficient **+0.869564412758** against ours **−0.869564412763**, so `\|ours + clm\|` is **4.773e-12** and `\|ours − clm\|` is **1.739**; clm's thresholds match `ours.alpha` to **7.177e-10**, unnegated. The coefficients are negated and the cutpoints are not — §7.2's asymmetry, confirmed by a SECOND reference implementation and this time on the WEIGHTED fit, which no other oracle can check | §7.2, §14.5, §18b, T7 | yes, run 2026-08-24 |

**Claims checked by reading the committed repository:**

| Claim | Where used | Verified |
|---|---|---|
| `data.KINDS` is the declared **nine** and `data.py:143-152` states that `model` describes what was fitted and that "Stages 7-13 all render under `model` or under an existing kind" | §9.1, §12, §14.9 | yes, read |
| `model.py:7-14` reserves this module for Stages 8, 9 and 12 entering directly, and `model.py:47-50` states that they inherit invariant 4's denylist assertion by **calling** `design` | §0.1, §4.5 | yes, read |
| `propensity.ess`'s docstring (`propensity.py:107-109`) is where the phrase "no-exposure-no-weights rule" appears, and it appears in an argument about a Kish sum rather than about a fitter | §0.1, §12 | yes, read |
| `propensity.py:400-432` is `_record_exclusion` with the `set` and the `notna()` that make its names-its-cases check able to fail; `propensity.py:440-443` states why the exclusion entry is recorded before the design | §9.2 | yes, read |
| `propensity.py:468-469` states that [§10] catches `FitError` to drop and count a replicate and "may catch nothing else: a SchemaError here is a bug in the resampler" | §4.5, §11, §14.7 | yes, read |
| `config.py:409-426` holds eight outcomes with exactly one `family == "primary"`, and its comment states that [§8] rests on there being one; `config.py:410-411` gives `mrs_90d` `kind="ordinal"` and `higher_is_better=False` | §7.1, §12 | yes, read |
| `config.py:281-304` gives the `FIRTH_*` block's argument that a tolerance is part of the estimator because [§10] refits in every replicate, and that the prefix is chosen so `model.py` reads no exposure-flavoured constant | §12, §6.3 | yes, read |
| `config.py:553` declares `MRS_THRESHOLDS` as `(0,…,5)` commented "cumulative RD_k [§8]", inside the `[§13]` subgroups block; `config.py:523-533` declares `PLAUSIBLE_RANGES["mrs_90d"]` as `(0, 6)` | §8.3, §12 | yes, read |
| `config.py:250-256` records that the workbook holds the literal `'N/A'` in `mRSscoreat90days` in **2 cells**, and `NA_VALUES` is load-bearing on the primary outcome | §4.1 | yes, read |
| `tests/test_cohort.py:116-133` is `cohort_frame(**overrides)`, building its cohort records from the `HAND-4` template; `:147` is `built` and `:152` is `set_cell` | §14.0, §14.2 | yes, read |
| `pilots/analysis.py:272-284` is `WeightedOrderedModel` overriding `loglike`; `:287-289` is `_wmean` with the zero-weight guard; `:291-293` is `weighted_rd`; `:442-457` is `cumulative_mrs_rd` on the weighted empirical route; `:460-475` is `weighted_common_or` with `nunique() < 3`, `except Exception: return np.nan`, and `np.exp(-res.params.iloc[0])` | §18, §7.2 | yes, read |
| Stage 7 §11 states that Stage 8 "does not read the `Balance`", and Stage 7 §12.0.2's rule that a golden vector comes from the synthetic frame so no patient-derived number enters git | §0.2, §14.0.2 | yes, read |
| SAP [§8] prescribes the weighted proportional-odds model with treatment as the sole predictor, `RD_k` from the weighted empirical distributions with intervals and **no p-values**, and no penalty; [§14a] writes the model as `α_k + βA + γᵀX`; [§10] defines the primary p-value from the bootstrap distribution of `β` | throughout | yes, read |

### 20b. This document's own code, executed — 2026-08-21

A specification that ships code is a specification whose code has to have been executed, which is
Stage 6 §18c's rule and Stage 7 §18b's practice. Every `python` fence in this document except the
signature listings of §3, the `config.py` block of §12 and the two `test_outcome.py` blocks of §14 was
extracted verbatim, assembled into one module against the landed Stages 1-7, and run — with §12's
block supplying `PRIMARY_OUTCOME`, `MRS_LEVELS` and the `POLR_*` constants from the same text T1 will
paste into `config.py`. **The fence map is derived by content and not by index** (Stage 7 §18d's rule).

| Block | Executed result |
|---|---|
| **All 29 fences parse and the 21 that form a module assemble and run** | `ast.parse` on each, then one `exec`. The 8 excluded are §3's two import blocks and its one public-surface signature listing, §5.3's one-line excerpt of §5.2, §9.1's three module constants, §12's `config.py` block, and §14's two `test_outcome.py` blocks — each identified by its first line, not by its index |
| §5.3's one-line excerpt | **A correction to this document, found by running it.** It was written indented, as it appears inside `polr`, so it did not parse standalone and broke §3's "every `python` fence in this document is valid Python" rule at the first fence that quoted another. Dedented, with a `# §5.2's line` marker so it is still readable as a quotation |
| `polr` on `ordinal_cohort()`'s cohort | 4 cutpoints over `(0, 1, 3, 4, 6)`, beta 1.7516026138, `converged_on` "likelihood" at iteration 4, rescales 0, halvings 0, alpha ascending — §14.0.2's golden vector to every digit printed |
| `primary` end to end on the same cohort | Returns a `Primary` — beta 1.7516026138, `odds_ratio` 5.7638324754, `cut_levels` (0, 1, 3, 4), all six `RD_k` matching §14.0.2; appends exactly **three** `model` entries in the declared order `outcome_completeness`, `primary_fit`, `cumulative_rd`; `in_estimate` is `bool`, total, and 5 of 5; and leaves the frame **equal cell for cell** to a deep copy taken beforehand (§0.2, §14.11) |
| `_assert_polr_fittable`, six routes | Each raises naming itself; O3 fires on the `nan` weight that would otherwise delete a category |
| `_assert_exposure_survived` | Raises `FitError` naming G6 on the width-0 design a one-armed frame produces |
| `_assert_reportable` | Raises `FitError` naming G7 on the separated fit, and passes on every healthy one |
| `_orientation()` | Returns the sentence, and raises on a patched registry with `higher_is_better=True` |
| The three rendered tables | `outcome_completeness` **3 / 3 / 3** markdown cells across header, separator and body; `primary_fit` **4 / 4 / 4** over 10 rows; `cumulative_rd` **4 / 4 / 4** over 7 rows. No cell of any of the three holds a `\|`. Entry 2 renders one row per `MRS_LEVELS` entry, with levels **2 and 5** reading `missing` under `absent from this population — no cutpoint [§5.3]` and level 6 under `highest declared level — no cutpoint` |
| The three `detail` strings | All three build and all three terminate. Entry 2's absent-levels clause renders as *"…so 4 cutpoint(s) were estimated; 2 declared level(s) are absent from this population and are marked as such below: 2, 5."* — read as prose, with the count matching the two `absent` rows in the table beside it. This is the branch Stage 7 §21.5a lost a full stop in, exercised here because `ordinal_cohort()` reaches it |
| `Audit.to_markdown()` through all eight stages | Renders **33,658** characters on `ordinal_cohort()`, with `## Fitted models` carrying all **nine** `model` entries — 4 from `fit`, 2 from `assess`, 3 from `primary`. Its ledger is **25 = 20 after `fit` + 2 + 3**, where the workbook's is **27 = 22 + 2 + 3** (§11): the two totals differ because the fixture's earlier stages record fewer entries, not because this stage's three change. The coincidence that 25 is also §14.0.3's workbook-fixture count — 22 + 3, that fixture not calling `assess` — is a coincidence and is named here so nobody reconciles the two |
| `data.KINDS` | still **nine** — §9.1's claim, and `data.py` is untouched |

**What this does not establish.** It ran the assembled fences against fixture-derived and synthetic
frames, not `outcome.py` and `model.py` as modules with their own imports, their `Final` constants and
the AST scans §14.11 and §14.13 ask for — that is T1-T6's work and §14 is what checks it. **And it has
not run against the workbook**, deliberately: §4.3 forbids this document's own harness from computing
the estimate, so DoD-4 and T8 are the first execution that does.

### 20c. The review round's own fences — 2026-08-21/22

§21.1 changed four fences and added five, so §20b's assembly was stale on nine — and **§20b was
therefore re-run in full rather than stood in for.** It found one defect, in code the round itself had
added: §14.0.4's two covariate frames were **perfectly separated**, and every oracle standing on them
was void. That is item 11, and it is the reason this section exists as a re-run and not as a parse.

| Check | Result |
|---|---|
| **All 34 `python` fences parse** | `ast.parse` on each, run after the round's edits. **0 failures.** The count moves **29 → 34** on five additions: `_weighted_categories` (§5.3a), §8.1's arm-code constants, §14.0.4's generator block, `BAND_EFFECTS`, and §14.0.5's `truncated_cohort()`. The module-assemblable count moves **21 → 23** — the first two additions belong in the module, the last three are `test_model.py` and `test_outcome.py` blocks and join §20b's exclusion list, which moves 8 → 11 |
| **§9.2's corrected status branch, both versions, on all three frames** | Run standalone. The branch below renders correctly on the workbook's `(0,…,6)`, on `ordinal_cohort()`'s `(0,1,3,4,6)` and on `truncated_cohort()`'s `(0,1,2,3,4)`, and the **estimate-column-iff-no-cutpoint invariant holds on all three**. The draft's branch **violates it on the third**, rendering `alpha at mRS <= 4 \| 4 \| missing \| fitted` — which is §21.1 item 6 measured rather than argued |
| **Exactly six distinct status strings across the three frames plus the two coefficient rows** | Enumerated. §14.9's count is asserted against that six, and **no status contains a `\|`** — checked, including `P(Y <= it) = 1` |
| **`truncated_cohort()`'s spread gives `(0, 1, 2, 3, 4)`** | Survivors at index positions 0, 1, 4, 6, 8 take `[0, 1, 4, 2, 3]`. And the same arithmetic reproduces `ordinal_cohort()`'s `(0, 1, 3, 4, 6)` from §14.0.1's recorded `[0, 1, 4, 6, 3]`, which is how the position inference was checked rather than assumed |
| **The cycling alternative to §14.0's spread** | Gives `[0, 1, 4, 6, 1]` and `categories (0, 1, 4, 6)`, moving §14.0.2's golden vector. Declined; §21.1 records it |
| **`Σw` on the workbook's estimation population** | **27.736623**, identical over `in_model`, weights in `[0.0269724, 0.964955]`. Read by a probe that reports three aggregates and no arm split and no outcome, so §4.3 is untouched. This is §5.2's new argument's one measured input |
| **93 / 92 / 92 and seven occupied levels, re-confirmed** | The same probe, independently of §20's. Both agree |
| **THE ASSEMBLY COMPOSES.** 24 of the 34 fences executed as two modules against the landed Stages 1-7 — 7 into `model.py`, 16 into a fresh `outcome` module, plus §12's `config.py` block supplying `PRIMARY_OUTCOME`, `MRS_LEVELS` and the `POLR_*` constants; 10 excluded | `_weighted_categories` resolves for **both** callers, `model_fit_error` binds to `FitError` as §5.2 says it must, and `_TREATED = 1` / `_COMPARATOR = 0` derive from `TREATMENT_LABELS` as §8.1 says. §12's two in-file assertions hold: `MRS_LEVELS[:-1] == MRS_THRESHOLDS`, and `8.79 < 14.0 < 18.81` |
| **§5.2's derivatives, at THREE parameter values on TWO frames** — §21's review risk 1, closed by measurement | Worst `\|Δscore\|` **4.260e-08**, worst `\|ΔHessian\|` **3.790e-08** against central differences, worst asymmetry **1.776e-15**, and the Hessian **negative definite at all six points** with eigenvalues comfortably away from zero (max −0.11 on `hand_ordinal()`, max −21.4 on `replication_frame()`). Points: the §5.2 start values, the returned optimum, and the optimum + 0.5. Frames: `hand_ordinal()` (K=6, m=1) and `replication_frame()` (K=4, **m=2** — the two-covariate case the [§8] path never has and Stage 12 always does). §14.3's tolerance is 1e-6, so both PASS with two orders of margin. **The `alpha`-`beta` cross block is correct**, which the round had also re-derived by hand; the measurement is what makes it a finding rather than a second opinion. **And a near-zero maximum eigenvalue on the FIRST run is what exposed item 11** — a collapsed information matrix is §6.1's signature, and it had no business appearing on a well-behaved oracle frame |
| **§14.0.2's golden vector reproduces to EVERY DIGIT** | `primary()` end to end on `ordinal_cohort()`: beta **1.7516026138**, `odds_ratio` **5.7638324754**, alpha **[−4.0128444951, −1.8393869358, −0.7580715177, 0.4483367817]**, `cut_levels` (0,1,3,4), `categories` (0,1,3,4,6), 4 iterations on "likelihood", rescales 0, halvings 0, `first_step_norm` 2.51298, and all six `RD_k` — `+0.1079772317, +0.5469186443, +0.5469186443, +0.0181363724, +0.4712177281, +0.4712177281`. The five survivors are HAND-1, HAND-2, HAND-5, COHORT-1, COHORT-3 with `mrs_90d` [0,1,4,6,3] |
| **`truncated_cohort()` fits `(0, 1, 2, 3, 4)` and the corrected branch renders right** | Measured through the real `_fit_table` and `data._md_table`, not standalone: `mRS <= 4` reads `missing` under **"highest FITTED level — no cutpoint"**, `mRS <= 5` under the absence status, and `mRS <= 6` under **"absent … and it is the highest declared level"**. The estimate-iff-no-cutpoint invariant holds over all 7 alpha rows on **both** fixture frames. beta 0.3447607100 there, 3 iterations |
| **All three tables render uniform, on both frames** | `outcome_completeness` 3 columns over 5 rows, `primary_fit` 4 over 10, `cumulative_rd` 4 over 7; header, separator and body agree in every case and **no cell of any table on either frame contains a `\|`** |
| **All seven guards fire, with the right type, naming themselves, appending nothing** | G1 (an `e` reindexed alone), G2 (an `object` mask carrying `None`), G3 (`mrs_90d` dropped), G4 (one arm), G5 (a hand-built `nan` weight inside `in_estimate`) → `SchemaError`; G6 (a width-0 design) and G7 (`separated_frame()`) → **`model.FitError`**, which is the type [§10] drops on. Every one appended **0** entries to a fresh `Audit` |
| **`separated_frame()` reproduces §6.1 and §20 EXACTLY** | crossovers 0: beta **36.4058**, `exp(beta)` **6.469e+15**, **17** iterations on **"score"**, rescales **0**, halvings **0**, min fitted probability **1.0000**. crossovers 1/2/3: **21.1512 / 20.4154 / 19.9711** with min p **0.05 / 0.10 / 0.15**. Every digit §6.1 records, from a generator written independently of the harness that first produced them — which is what makes §14.7's value assertions safe to keep as values |
| **`orientation_frame()` is downward by construction and the reciprocal identity holds** | Treated mean mRS **1.0400** against control **3.7700**; beta 4.066383, `exp(beta)` **58.35 > 1**; the `alpha_k − beta*A` reading **0.017139 < 1** and equal to `1/exp(beta)` to 1e-12; `\|beta\|` inside the bound. §20's harness measured 1.11/3.82 and 84.58 — different sampling, same properties, which is exactly why §14.5 asserts the properties |
| **§14.7's band probe: 240 fits in 0.69 s** | 153 below the bound, 87 at or above, max below **8.5534**, min above **20.7879**, **0** in `(10.0, 18.0)`, both modes populated, bound strictly inside the observed gap. Bimodal as §6.3 describes. **`per_effect=400` — §6.3's calibration — measures at 13.9 s**, worse than the 7 s the round estimated, so the split of §21.1 item 3 is confirmed by the thing it was arguing about |
| **The three `detail` strings build, terminate and read as prose** | 542, 1267 and 1059 characters, all ending in a full stop. Entry 2's absent-levels clause renders *"…so 4 cutpoint(s) were estimated; 2 declared level(s) are absent from this population and are marked as such below: 2, 5."* — the Stage 7 §21.5a position, with the count matching the two `absent` rows beside it. Entry 3's per-arm totals clause renders *"…denominators [§11]: EVT alone 0.545537, bridging 0.657195."* |
| **Three entries, declared order, `KINDS` still nine, frame untouched** | `outcome_completeness`, `primary_fit`, `cumulative_rd` at positions captured before the call on both frames; `data.KINDS` still the declared nine; the frame **equal cell for cell** to a deep copy with identical columns; `in_estimate` boolean, total, on the frame's index |

**What 20c does not establish.** It ran against fixture-derived and synthetic frames only. It did not
run `outcome.py` and `model.py` as files on disk with their own imports, their `Final` constants and the
AST scans §14.11 and §14.13 ask for — the `outcome` module here is assembled in memory, so a stray
import or a missing `__all__` would not show. **And it has not run against the workbook**, deliberately:
§4.3 forbids this document's harness from computing the estimate, so DoD-4 and T8 remain the first
execution that does. T1-T6 close the rest; DoD-18 is the item that says so.

## 21. What this spec changed elsewhere

**Nine** entries, all landing with this document rather than with the implementation, because a
specification that contradicts the roadmap is worse than no specification. Precedent: Stages 3, 4, 5,
6 and 7 landed their roadmap edits with the document.

The count reached eight and not six because applying this section to the roadmap added two the table
below did not originally carry — items 7 and 8 — which is Stage 7 §21.6's rule turned on this document
itself: *"after applying a review round, re-read the round against the file and check each item is in
the file rather than in the argument for it."* Both were decisions this spec had already taken in §0.1
and §4.3 and neither was in the ledger, which for a document declaring itself the sole source is
indistinguishable from never having found them. **The ninth arrived from the review round of §21.1**,
and it is the only one of the review's thirty-six items that changes a file outside this stage: it
carries §11's handover into Stage 10's own roadmap entry, where Stage 10's implementer will read it.

| # | Change | Files | The question it settles |
|---|---|---|---|
| 1 | Roadmap Stage 8's **Accept when** gains the **separation rule**: a fit whose treatment coefficient reaches a declared bound raises rather than returning, because an unpenalised proportional-odds fit on a separated sample **converges** rather than failing | `implementation_roadmap.md`, §6 here | Whether [§10]'s "replicates whose prespecified fit fails are dropped and counted" covers a degenerate ordinal fit. **It does not**, and the roadmap's silence reads as though it did. Measured: a perfectly separated frame converges in 17 iterations on the score criterion with every safeguard counter at zero, every fitted probability finite and **nothing anywhere out of range**, and returns `exp(β) = 6.5e15`; 400 sparse replicates of the cohort's shape produced **zero** failures. Without a bound there is no failure for [§10] to count, and the percentile interval is a quantile of a distribution with a tail at exp(21). **The conditioning is deliberately not cited** — an earlier draft of this row said "the information matrix conditioned at 6.85", which §6.3 strikes as a 2×2-against-7×7 artefact, and the roadmap as landed carried the excision as a dangling clause: §21.1 item 7 |
| 2 | Roadmap Stage 8's **Accept when** replaces "the weighted fit equals an unweighted fit when all weights are 1" as the primary weighting criterion with **integer weights equal row replication**, keeping the first as its companion | `implementation_roadmap.md`, §14.3 here | Whether the criterion as written can detect an implementation that ignores the weights. **It cannot**: measured, `polr(X, y, ones)` equals `polr(X, y, None)` bit for bit, so a function that never reads `w` passes it. The replication form agrees at **0.000e+00** and fails for such an implementation, and it is the only oracle available that checks the weighting at machine precision |
| 3 | Roadmap Stage 8's **Accept when** replaces "the cumulative probabilities are monotone in `k` within each arm" with the assertion that the **weighted empirical route was taken** — the `RD_k` against a hand computation, and against the six-threshold-model route, shown to differ | `implementation_roadmap.md`, §8.2 here | Whether monotonicity can detect the route [§8] forbids. **It cannot, in two directions.** The empirical route is monotone by construction — 0 of 2000 samples non-monotone — so the criterion cannot fail on a correct implementation; and the threshold route crosses in only **8 of 1990** samples and by less than **5e-07**, so it does not reliably fail on an incorrect one either. An earlier draft of this section claimed 40 of 40 samples crossed and was reporting a companion that computed `1 − P(Y ≤ k)`; the claim is struck and the corrected rate recorded |
| 4 | Roadmap Stage 8's **Accept when** gains the **orientation companion**: the test needs a construction whose truth is known **and** a sign-flipped reimplementation shown to report the opposite direction | `implementation_roadmap.md`, §7 here | Whether "the orientation test confirms `exp(β) > 1`" is enough. It is the one criterion of the four that can fail as written, and it still needs its companion: every reference implementation parametrises `α_k − x'β`, measured — the coefficients negated and **the cutpoints not** — and `pilots/analysis.py:475` carries the negation that is correct for `statsmodels` and inverts the primary result under [§14a]'s own parametrisation. A test that only checks `exp(β) > 1` on a large effect passes for an implementation returning `exp(\|β\|)` |
| 5 | Roadmap Stage 8 gains the **[§11] denominator** item: the estimate's population is a third population after the cohort's 93 and the ATO's 92, and it is named, counted and logged | `implementation_roadmap.md`, §4.1 here | Where [§11]'s "complete-case, with the denominator reported for every estimate" lands at this stage. Stage 6 §4.4 stated the second population rather than leaving it implicit; the third was in neither the roadmap nor any spec. On v7 it **equals** the second, which is exactly the condition under which an implementation that never built the mask is green |
| 6 | Roadmap Stage 8's "weighting the log-likelihood contribution per observation is sufficient" is **corrected**: sufficient for the point estimate, **not** for a standard error, because the weights are a tilting function of an estimated propensity score and not frequencies | `implementation_roadmap.md`, §5.6 here | Whether "sufficient" licenses a model-based interval. It does not, and the sentence as written is the one a reader would take it from: `−H⁻¹` at the optimum is the inverse observed information of a likelihood in which `w_i` counts observations, and it omits the variability of estimating `e` entirely — which [§7] states enters the interval. [§10]'s percentile bootstrap is the prespecified interval and there is no second one |
| 7 | Roadmap Stage 8's **Build** gains the **two-module split**: the fitter is `model.py`'s and is general in its covariates from the first line, and the [§8] specification is `outcome.py`'s, which Stage 9 extends | `implementation_roadmap.md`, §0.1 here | Where the estimator lives, which the roadmap left unstated where it stated it for Stage 6. [§14a] prescribes the same proportional-odds model with a wider design and unit weights, so a treatment-only fitter written now and generalised at Stage 12 means two implementations of one estimator — and [§14a]'s guard that `exp(β)` there is a *conditional* odds ratio is a statement about the same `β` this stage returns. `model.py:7-14` already reserved the boundary; this records that Stage 8 is where it is first used |
| 8 | Roadmap Stage 8 gains **"On the estimate not being in the specification"**: the spec does not measure or quote the workbook's `β`, `exp(β)` or `RD_k`, and its data-gated tests assert properties rather than values | `implementation_roadmap.md`, §4.3 here | Whether the house rule "every number was produced by running code" applies to the answer itself. Every earlier stage recorded its decisions as taken before any outcome was examined by arm, and this is the stage where that ends — so the question is whether it ends while the estimator is being specified or after. Three decisions in this document would have been contaminated by the first order: §5.3's collapse, §6's bound and §8's thresholds. The roadmap needed to carry it because a future reader comparing this spec against Stages 1-7 would otherwise read the missing numbers as an omission |
| 9 | **Roadmap Stage 10 gains the three things Stage 8 hands it** — two failure counters reported separately, the distribution of `len(fit.alpha)` across replicates, and `FitError`-not-`SchemaError` — with an **Accept when** item requiring the separation count to be surfaced and exercised | `implementation_roadmap.md`, §11 and §15 here | Whether §11's handover survives being in a spec that Stage 10's implementer may not read. It is the ninth entry because §21.1 item 10 found a gap §11 could state and the roadmap could not: Stage 8's bound is calibrated at seven occupied categories and §5.3 makes the category count a property of the replicate, so **[§10] is where that bites and Stage 12 is only where §15 had filed it.** Reporting the cutpoint-count distribution answers it from replicates already drawn rather than from a second sweep, and it costs Stage 10 one column |

**And one roadmap addition that is not an amendment**: Stage 8 gains its
`**Spec:** specs/stage8_primary_outcome_estimator.md` line, as Stages 1-7 have.

**What this spec has had, and what it has not.** Every numeric claim in §3.1, §5, §6, §7, §8, §14 and
§20 was produced by running code — against fixture-derived cohorts, against hand-built vectors, against
some 6000 synthetic fits, and against `statsmodels` and `scipy`. **Five of them are corrections to this
document's own earlier drafts, and every one was found by running it rather than by reading it:**

1. **The roadmap's implied convergence failure** — §21's amendment 1, and the largest finding in this
   stage. An unpenalised ordinal fit on a separated sample converges; [§10] had nothing to count.
2. **The `cond(-H)` detector "does not discriminate"** (§6.3). It does, by six orders of magnitude. The
   earlier claim compared a 2×2 matrix against 7×7 ones and got both of its numbers wrong. The detector
   is still rejected, and now for reasons that survive being measured.
3. **`POLR_MAX_ABS_BETA = 10.0`** (§6.3, §12). Calibrated on 500 fits at a **single** true effect size
   of 0.5 — itself an assumption about the answer, which §4.3 forbids — it sat 14% above that
   measurement's maximum and at the lower edge of the real band. Re-measured over twelve effect sizes
   and 4800 fits, the region is `(8.79, 18.81)` at six cutpoints and `(11.04, 18.98)` across all of
   them (item 24), and the bound is **14.0**.
4. **The 40-of-40 threshold-model crossing rate** (§8.2). A complementary-probability bug in the
   companion: it computed `1 − P(Y ≤ k)`, so every sequence decreased by construction. The corrected
   rate is 0.40% with the largest decrease below float noise, which changes what §21's amendment 3 can
   claim.
5. **That a cutpoint crossing was reachable by removing the halving** (§5.4). 3000 attempts could not
   produce one. Struck rather than weakened.

**Two of those five moved a number that another section asserted**, which is why §20 is declared
normative and why the second and third are recorded as corrections rather than quietly fixed: a bound
chosen against one side of a band is the kind of error that reads as rigour.

Stage 6 §19 predicted that yield for a second pass and Stage 7 §19 recorded that two independent rounds
found seventeen and twelve **disjoint** items. **This document has now had one**, and §21.1 is its
ledger. What it did *not* have before that round is stated rather than implied: §20b ran its own fences
and §6.3's re-measurement was this document checking its own arithmetic, neither of which is a review.

The round was pointed at the three places the previous version of this paragraph named as where a first
one should start, and **two of the three were productive**:

- **§5.2's Hessian assembly** — the derivation was re-derived independently and every block is correct:
  `Huu`, `Hll`, `Hul`, the `β` block's `Huu + 2·Hul + Hll`, and the `α`-`β` cross block's per-side
  `Huu + Hul` / `Hll + Hul` under the `both` mask. Nothing moved. What moved is the **test**: §14.3 now
  checks the derivatives at three parameter values on two frames rather than at the optimum on one,
  because the original objection stands even though the arithmetic did — a wrong sign that vanishes at
  one point would have survived (§21.1 item 8).
- **§14.12's property assertions** — no assertion was found insufficient, and the section's *mechanism*
  was: the float-literal scan had no boundary, no allowlist and no non-empty check, which is §21.1
  item 9.
- **§6.3's band, on a design that is not one column** — unchanged, and the round found a **nearer**
  version of it that the paragraph had not: the band is calibrated at `K = 6` only, and §5.3 makes `K`
  a property of the replicate, so [§10] is where it bites first and Stage 12 second (§21.1 item 10).

**Three places a second round should start**, and none of them is where the first one did:

- **§5.2's `polr` loop against §14.4's counter assertions.** The loop's `moved` and `g` are both bound
  before the step, `halvings` is incremented on the 30th rejection that is never re-evaluated, and the
  `for…else` raise is reached without a final likelihood evaluation. None of that is wrong; whether
  §14.4's four diagnostic assertions actually pin the values the loop produces, or merely their types,
  is not established.
- **§9's four table helpers against `data._md_table`'s padding.** §14.9 asserts cell counts and the
  absence of `|`. It does not assert that `_fmt` is the only formatter reaching those cells, and
  `_fit_table` calls `_fmt` on `fitted.get(level)` — a `None` — while `_rd_table` calls it on floats,
  so the two go through different branches of `_fmt` and only one of them is exercised on every frame.
- **§10's `primary` against a `Propensity` whose index is a superset of the frame's.** G1 asserts
  `index.equals`, which is the right check; what is not established is that every later `.loc` in §10
  would fail loudly rather than silently if G1 were relaxed, and §14.14 tests G1 firing rather than
  what it protects.

### 21.1 The independent review rounds — 2026-08-21/24

**Thirty-six items: thirty-five applied to this document, one raised and declined.** Twenty-three are
numbered below because they changed a decision, a specification or a test; eleven were prose,
cross-reference and rendering corrections; one added code; one was declined. **Items 12-22 come from a
second, independent round** (§21.2) run against items 1-11's own output. Eight of the eleven came from
that round and three more were found while verifying its findings — and **four of the eleven are
defects in text items 1-11 had themselves added** (12, 19-20's sizing, 21, 22), which is the pattern
worth naming: a repair pass writes new code, and new code is unreviewed code. The last three groups are listed after the
table. **Six of the twenty-two would have reached the implementation.** Item 6 is a logic error in a shipped
table; items 1 and 2 together mean three acceptance sections named a fixture that does not exist, with
the numbers of one that does; item 11 was introduced by this review and caught only by running it —
two oracle frames that were perfectly separated, on which the oracles asserted nothing; **item 13 is an
acceptance assertion that is false and that §14.12 aims at the workbook**, where it could fail for a
reason the analysis already declares untested; item 14 leaves two safeguards on the bootstrap path with
no test that can fail; and item 12 leaves an oracle unable to run at all.

| # | Item | Where | What it changed |
|---|---|---|---|
| 1 | **Five of the six synthetic constructions were specified in prose only, and §14 asserted their numbers.** "A 200-record construction with the treated arm's mass at mRS 0-3 and the control's at 2-6" does not yield `exp(β) = 84.584014`, and an implementer cannot reach it. A document that declares itself the sole source and pins a number produced by a construction it does not specify has a hole exactly the shape of the number | **§14.0.4** new; §7.2, §14.3, §14.5, §14.7, §18b, §20, T2 | Six deterministic generators written out, all seeding from `C.SEED`, in `test_model.py` beside the sections that consume them. Every §14 assertion resting on a seeded generator becomes a **tolerance or a property**; §20's values become the record of what was seen rather than a prediction. The two seed-free constructions keep their values |
| 2 | **§14.3, §14.6, §18b and §20 all called the twelve-record hand frame "§14.0.2's".** §14.0.2's frame is `ordinal_cohort()`'s cohort and has **five** records with four cutpoints. §14.6's six hand-computed `RD_k` assertions would have been written against the wrong frame, and §14.0.2 already pins its own `RD_k`, so both would have looked plausible | §14.3, §14.6, §18b, §20 | The frame is named `hand_ordinal()` and written out in §14.0.4 with its twelve weights as two-decimal literals, so §14.6's assertion is hand-checkable and is the one value assertion in that section |
| 3 | **§14.7 re-ran §6.3's calibration on every commit.** Its assertions ranged over 4800 fits, 500 healthy `cond(-H)` fits and 400 sparse replicates run twice — about 6100 fits, 7 s at §2's measured rate, roughly tripling a suite that returns in four | §14.7, §12, T3, DoD-13 | Calibrations stay in §20; the suite gets probes at 240 / 20 / 40 fits, under a second. `cond(-H)` is asserted as orders of magnitude rather than as the sampler's ranges. §12's `test_config.py` literal `8.79 < POLR_MAX_ABS_BETA < 18.81` is now stated to be the only guard on the calibration, and DoD-13 times the suite |
| 4 | **§9.2's status column had five distinct strings and §14.9 asserted "the four declared strings".** The two coefficient rows carry statuses too | §14.9 | Six after item 6's fix, asserted against a literal tuple of the six |
| 5 | **§14.14's bullet was titled "A failing `primary` records nothing" and its own body said it does not** — entry 1 is present after a G6 or G7 failure, deliberately | §14.14, Failure modes | Retitled and split into two halves with explicit deltas: exactly 0 on G1-G5, exactly 1 on G6 and G7 with the step name asserted and entries 2 and 3 asserted absent |
| 6 | **`_fit_table`'s no-cutpoint branch tested `level == MRS_LEVELS[-1]` where it must test `level == fit.categories[-1]`.** The two are equal on every frame this stage has — the workbook occupies all seven levels and `ordinal_cohort()`'s top category is 6 by coincidence — and come apart in any [§10] replicate that draws no death. Measured on §14.0.5's frame by running both branches: the draft renders `alpha at mRS <= 4 \| 4 \| missing \| fitted` and marks mRS 6 as the highest declared level while hiding that nobody is in it. Nothing raises, and §14.9's absent-row assertion **passes on the broken table** because the broken row is not an absent one | **§9.2**, §14.0.5 new, §14.9, T6, DoD-12, Failure modes, §19's code comments | Branch rewritten to key on the fitted set, with a third no-cutpoint status distinguishing an absent declared maximum. `truncated_cohort()` added as the frame that witnesses it. §14.9 gains an estimate-column-iff-no-cutpoint invariant, which catches this class of error rather than this instance of it |
| 7 | **§21 item 1 cited "the information matrix conditioned at 6.85" as evidence a separated fit looks healthy — the exact claim §6.3 strikes** as a 2×2-against-7×7 artefact. And the roadmap as landed shows the excision: `implementation_roadmap.md` reads "every fitted probability finite and the information matrix / and nothing out of range", a dangling clause. Stage 7 §21.5a's failure mode, third occurrence | §21 item 1, `implementation_roadmap.md` | The conditioning is no longer cited in either place, and §21 item 1 says why it is not, so the next reader does not restore it from §6.1's block |
| 8 | **The Hessian's `α`-`β` cross block was checked against central differences at one parameter value on one frame** — which the previous §21 named as review risk 1 and did not close | §14.3, T2 | Three parameter values on two frames: `par = 0`, the returned optimum, and the optimum perturbed by 0.5 in every coordinate. The derivation itself was re-derived in review and is correct |
| 9 | **§14.12's float-literal meta-assertion had no boundary, no allowlist and no non-empty check.** "Scans this section's source" is unwritable as stated, and §14.12 itself needs `1e-12`, `0.0` and `1.0` | §14.12, DoD-17 | Boundary is banner-to-next-banner in **file order** (`test_balance.py`'s banners are not in numeric order); both delimiters and a non-empty range are asserted found; allowlist is `pytest.approx` / `rel=` / `abs=` values plus `0.0` and `1.0`; the docstring-versus-comment asymmetry is stated |
| 10 | **§6.3's band was measured at `K = 6` only, and §15 filed only the wider-design version of that gap.** §5.3 makes `K` a property of the sample and [§10] draws 2000 samples, so the nearer question is whether the band is empty at `K = 5` — and G7 with a fixed bound is what decides which replicates are dropped | §5.3, §11, §15 | Filed in §15 as its own bullet, and §11 asks Stage 10 to report the distribution of `len(fit.alpha)` beside its G7 count, so the question is answered with the replicates it already draws rather than by a second sweep |
| 11 | **INTRODUCED BY THIS ROUND AND CAUGHT BY RUNNING IT.** §14.0.4's `replication_frame()` and `reference_frame()` — added by item 1 to make the oracles reproducible — built `y` as `searchsorted(quantiles(eta), eta)` with **no noise term**. That makes `y` a monotone function of `x'b`, every category an interval of the linear predictor, and the frame **perfectly separated in the ordinal sense**: the likelihood has no interior maximum. Measured (§20c): `beta = [-3281.7, -1652.8]` at n=120 and `[-48828, -24382]` at n=500, both converging on the score criterion with clean counters and both far past `POLR_MAX_ABS_BETA`. **Every oracle standing on those two frames was void** — two optimisers of a likelihood whose maximum is at infinity need not agree anywhere, so §14.3's 1e-10 replication tolerance and §14.5's `ours.beta + sm.beta ~ 0` would have been asserting agreement between two arbitrary stopping points | §14.0.4, §20c | A `rng.logistic(size=n)` term, which is the actual proportional-odds data-generating process: latent = `x'b + Logistic(0,1)`, cut at its own quantiles. With it, `max\|beta\|` is **0.93** and **0.76**, and re-measured (§20c): replication agrees at **0.000e+00** on beta and **8.327e-16** on alpha; `statsmodels` gives `\|ours.beta + sm.beta\|` **1.425e-05** and `\|ours.alpha - sm.thresh\|` **8.947e-06**, with `\|ours.beta - sm.beta\|` **1.52** so the sum test has something to negate. **§6.1's finding turning up uninvited in a test fixture is the argument for §20b's rule in one sentence:** a specification that ships code has to have run it, and this fence would have parsed, read plausibly, and produced a green-looking oracle that asserted nothing |
| 12 | **`hand_ordinal()`'s weights were not integers, so it could not carry the oracle §14.3 assigns it — and §20's row misattributed the frame.** §14.3 and §18b run integer-weight row replication "on `replication_frame()` and on `hand_ordinal()`"; `np.repeat` with a float weight raises `TypeError`. Worse, an earlier draft identified §20's twelve-record hand frame WITH `hand_ordinal()` — they are different frames, and §20's `beta 3.897` is not `hand_ordinal()`'s `−0.8696`, which is not even the same sign | §14.0.4, §14.3, §20, §20c | Integer weights `[1,2,3,4,5,6,6,5,4,3,2,1]`, Σw 21 per arm, chosen from the start. **Scaling the fractional weights up does NOT work and that is the finding inside the finding**: the MLE is invariant to a common weight scale in exact arithmetic, but `POLR_TOL` is absolute on a likelihood that scales with `Σw`, so ×20 gives `\|Δbeta\| = 9.458e-10` and FAILS §14.3's 1e-10. Integers from the start agree at **1.11e-16**. §20 gains `hand_ordinal()`'s own measured row; §14.3 now states its tolerance **per weight scale** — the first measured consequence of §5.2's absolute tolerances |
| 13 | **§14.10's and §14.12's cross-scale orientation assertion is FALSE, and §14.12 asserts it on the workbook.** "On a frame where `beta > 0`, the `RD_k` are non-negative at every `k` where the arms differ" is not a property of a proportional-odds fit — `beta` summarises a shift the data need not exhibit at every cut, which is the assumption §15 files as untested. Measured: `truncated_cohort()` has `beta` **+0.3448** and `RD_3` **−0.4531**; `hand_ordinal()` has `beta` **−0.8696** with `RD_0` and `RD_5` both **+0.0476`**. The suite's own fixtures are counter-examples | §14.10, §14.12, Failure modes | Asserted through the **weighted mean** instead: `sign(beta) == −sign(mean_w(Y \| treated) − mean_w(Y \| comparator))`, which holds on both fixtures (§20c) and is what §5.1's parametrisation actually implies. The strict per-`k` form is kept on `orientation_frame()`, where proportional odds holds by construction. **The workbook half is the severe one**: asserted on v7 it could go red for a legitimate reason, with a message saying "sign error", on the one stage where §4.3 makes that indistinguishable in advance |
| 14 | **§14.4 pinned only the TYPES of `rescales` and `halvings`, so both of §5.2's safeguards had no test that could fail.** Measured: deleting the trust region **and** the step-halving loop from `polr` leaves **every frame in the suite bit-identical**. Both counters are 0 everywhere the suite looks, and §14.4 cited Stage 6 §12.4a as precedent — where the same two safeguards ARE driven non-zero by dedicated fixtures | §14.4, §20c | Two fixtures, one line each, exploiting §5.2's own equivariance argument: `orientation_frame()` with the exposure scaled by 0.1 gives `rescales = 1`, by 0.01 gives `(2, 1)`. Asserted as values, with the halving-exhausted raise asserted under a patched `POLR_MAX_HALVINGS = 1`, and `beta * scale` asserted invariant to 1e-8 so the pair doubles as the scale-invariance test. §14.4's "zero of both counters on every frame" is struck — **10 of §14.7's own 240 band-probe fits carry a non-zero counter** |
| 15 | **`_completeness_table`'s docstring calls its four rows cumulative; they are not.** Row 3 counts the outcome over the WHOLE cohort, so adjacent differences are not restriction costs: on the workbook the four read 93, 92, **93**, 92 and the differences are −1, **+1**, −1 — a "restriction" that adds a record | §9.2 | The claim is struck, the row is retitled `mrs_90d present (whole cohort)` with its cell saying it is not a step in the chain, and the docstring says why it is kept parallel rather than made nested: deriving "how many outcomes are present at all" by subtraction is the reconciles-with-itself table Stage 5 §7.3 forbids |
| 16 | **§14.14's G2 companion is false.** It claimed that without the phase-1 raise the input "dies with pandas' `ValueError`". Measured for `None`, `pd.NA` and `np.nan` alike: `object & bool` **coerces the missing entry to `False`**, `in_estimate` comes back clean and boolean, and the record leaves the [§11] denominator with nothing raising | §4.4's G2 message, §14.14 | The companion asserts the silent exclusion — 4 estimated against `in_model`'s 5, no exception — which makes **G2 stronger, not weaker**. G2's message also named the wrong mechanism, "a three-valued mask resolves by branch order at the first `if`"; there is no `if` on the mask in §10, it resolves inside `&`, and the message now says so |
| 17 | **§14.14 tested that G1 fires and never what G1 prevents, and its executed witness was the one Series §10 never reads.** Every other guard here has a with-it-removed witness producing the wrong number; G1 had none, while its message claims it is the difference between an estimate and "a finite common odds ratio for a population that does not exist". `ps.e` is read nowhere in §10 | §4.4's G1 message, §14.14, §20c | A **permutation** companion: with G1 removed, a `Propensity` carrying the frame's labels reversed gives `odds_ratio` **58.01** against **5.76**, nothing raising, because `w` arrives in the `Propensity`'s order while `y` and the arm vector come from the frame's. A **subset** is the quieter second witness. A **superset is benign** — measured — so G1 is stricter than §10 needs, which is stated. G1's message no longer claims a "positional read"; §10 reads by label throughout |
| 18 | **§14.8's O2 companion cannot be written as specified: removing O2 does not reach a fit, because O4 catches a `nan` response too.** `np.mod(nan, 1.0)` is `nan` and `nan != 0.0`, so O4 fires on any non-finite response — measured, and the same holds for `inf` | §14.8, §14.0.4a | The companion removes **O2 and O4**, and then measures what the draft claimed: 79 of 80 fitted, `categories` unchanged, `beta` −1.162237 → −1.199291. The interaction is stated in its own right: O2 is the better message and the earlier position, not the sole detector, and O4 is a backstop nobody designed as one |
| 19 | **§14.8's O4 companion is unwritable at the size implied.** With O4 removed, six distinct half-integer values over 24 records makes the Hessian singular and `np.linalg.solve` raises `LinAlgError` — not `FitError`, and not the "one cutpoint per distinct value" fit the companion asserts | §14.8, §14.0.4a | `noninteger_response_frame()` is **3 distinct values over 60 records**, sized so the companion reaches a fit: measured, 3 categories and 2 cutpoints, one per distinct value |
| 20 | **§14.7's "a bound applied to `alpha` would reject a legitimately rare category, and the test constructs one" — no such frame is constructible.** `alpha` for a rare category grows like `log n`: measured, 1 rare record in 61 gives `max\|alpha\|` **4.16**, in 601 **6.40**, in 6001 **8.70**, in 60001 **11.00**. Reaching 14.0 needs order **1e7** records | §14.7, §14.0.4a | The assertion is made by **patching the bound**, not by finding the frame: on `rare_category_frame()` (`max\|alpha\|` 4.1599, `\|beta\|` 0.126998) the guard passes with `POLR_MAX_ABS_BETA` patched to **3.0**, where a bound on `alpha` would raise, and raises only when the bound goes below `\|beta\|`. That asserts the field the guard reads, which is the bullet's actual subject |
| 21 | **`band_samples` returned a `\|beta\|` vector, so §14.7's `cond(-H)` assertions had no subject** — the Hessian is a property of the fit and the fits were discarded | §14.0.4, §14.7 | It returns `list[PolrFit]`. The band assertion takes `max(abs(f.beta))` from each; nothing is lost, and §14.7's `cond(-H)` and counter assertions gain the objects they are about. Measured unchanged afterwards: 240 fits in 0.69 s, 153 below and 87 above, `max below 8.5534`, `min above 20.7879`, 0 in the gap, and **10 fits carrying a non-zero safeguard counter** |
| 22 | **§14.0.4 repaired six constructions and eight more assertions still named prose-only frames — and DoD-18 asserted the opposite.** The completeness claim was written in the same pass as the repair, so it described the repair's intent rather than its extent | §14.0.4a, DoD-18 | §14.0.4a adds six more hand-built generators plus `separated_cohort()` in `test_outcome.py`, all seed-free so their numbers are values; §20c records each one's measured output. DoD-18 now requires the check to be done **by enumeration** — walk §14 bullet by bullet and name the generator each frame comes from |
| 23 | **§15 and §18b recorded the R oracle as unrunnable, and it runs.** The draft measured `Rscript --vanilla` from a bare shell, where R genuinely segfaults during startup and `requireNamespace` returns `FALSE` for packages that are present — and concluded that Stage 7's `_r_environment()` repair "does not help". **It helps completely.** Through `R_ENV`, `_r_has(("ordinal",))` and `_r_has(("MASS",))` are both True, and the fourteen R oracle tests of Stages 6 and 7 pass in the suite today — the evidence was one `pytest -q` away | §18b, §15, §20, §19's T7 | `ordinal::clm(ordered(y) ~ a, weights = w)` on `hand_ordinal()` gives `\|ours + clm\|` **4.773e-12** and thresholds matching to **7.177e-10**. So the **weighted** fit is validated against an independent implementation — the thing §18b says only this oracle can do — and §7.2's negated-coefficients / unnegated-cutpoints asymmetry is confirmed by a second reference. **T7 moves P2 → P1**, §15's bullet is struck, and §20's `MASS is absent` row is corrected. The lesson is the review's own: a blocker measured once and filed is a blocker nobody re-measures |

**And the eleven smaller ones**, each applied: the positive-weight category rule was written out twice
with different casts and is now `_weighted_categories` with an AST scan (§5.3a — promoted to a numbered
subsection because §5.5's O4-before-O5 argument is a fact about that expression); `_orientation` was
missing from §3's private list that §14.11 asserts against; `first_step_norm` was carried, unlogged and
unasserted, and §14.4 now names all four diagnostics as Stage 6 §12.4a does; §5.2 did not record that
`g` is bound before the step, which `model.py:492-495` states for `firth`; §5.2 did not record that
`POLR_TOL` and `POLR_SCORE_TOL` are absolute on quantities scaling in `Σw`, now argued with `Σw`
measured at 27.736623 and §17 carrying the declined alternative; §5.5 did not record that O6 runs before
the collapse, now stated with the `w = 0` residual filed in §15; §11's "or accept 3 entries per
replicate" understated 6000 entries with tables; §17 was missing the scipy-gradient-only alternative,
which is the one option that would shrink this stage's highest-risk code; §12's `test_config.py` row
gained a sixth assertion, that `PLAUSIBLE_RANGES["mrs_90d"]`'s upper bound is not `None` (the declared
type permits it and `None + 1` is a `TypeError` in §12's computed `MRS_LEVELS`); §21's table rendered as
two tables, items 7 and 8 headerless, from a stray blank line; and §19's parallelism block did not show
that T4, T5 and T6 all depend on T2's generator block in `test_model.py`.

**One item was raised and NOT applied**, and it is recorded because a review whose declined items are
invisible is a review that looks unanimous. `ordinal_cohort()`'s `[: len(base)]` is a one-directional
guard: total if the committed fixture shrinks, a bare pandas length mismatch if it grows. Cycling
`MRS_LEVELS` over the index is total in both directions, and — measured — gives the surviving records
`mrs_90d` `[0, 1, 4, 6, 1]` instead of `[0, 1, 4, 6, 3]`, so `categories` becomes `(0, 1, 4, 6)` and
every number in §14.0.2's golden vector moves. A totality repair that moves the pin the fixture exists
to carry is not a repair. §14.0 records the hazard instead.

**And one item was raised and answered by measurement rather than by change.** `cumulative_rd`'s
`share[1] - share[0]` was flagged as two bare literals in a document that computes `PRIMARY_OUTCOME`,
`MRS_LEVELS` and the orientation sentence from registries. It is now `share[_TREATED] -
share[_COMPARATOR]` with both derived from `TREATMENT_LABELS` (§8.1) — the only one of the twenty-one
that added code rather than a test or a sentence, and it earns that because the order of that
subtraction *is* the sign §7 spends four hundred words on.

**The coverage-map count moves 62 → 71 and the unreachable-on-every-frame count 4 → 5.** The fifth is
item 6's branch, and it is the only one of the five that is reachable in an ordinary [§10] replicate
rather than only on a workbook that has not arrived — which is what took it from a note to a fixture.

### 21.2 The second round, and what it got wrong

Items 12-22 come from an **independent reviewer with no context on the first round**, given five scoped
targets and the standing list of settled decisions (§4.3, §5.6, §6.4, §17's declined rows, §14.7's
split, §14.0's declined change) with instructions not to re-litigate them. It was told where to look and
**not** what was suspected, so its conclusions are its own. It could read the repository and execute
code, and was held to the same §4.3 constraint: it did not compute the workbook's estimate.

**It was scoped this way because a second unscoped pass over 4,000 lines would have re-derived the first
round.** Stage 6 §19 and Stage 7 §19 record that two independent rounds found seventeen and twelve
**disjoint** items; this one found eight, and every one of them was in a place the first round had
flagged as unresolved or had itself just written. That is the argument for scoping a second reader
rather than repeating the first.

**What it got wrong, recorded because a review whose errors are invisible reads as an oracle.** Two
things:

- **It reported `test_balance.py:261,364` as a bad citation**, claiming line 364 is blank and the `12.1`
  banner is at 366. Checked: **364 is the banner** and the citation is correct. The finding was
  rejected and the citation left alone.
- **It proposed fixing item 12 by scaling `hand_ordinal()`'s fractional weights by 20**, reasoning that
  the polr MLE is invariant to a common weight scale. True in exact arithmetic and **false in this
  implementation**: `POLR_TOL` is absolute, so the scaled fit stops elsewhere and replication comes back
  at `9.458e-10`, failing §14.3's own 1e-10. Measured, and the fix used instead is integers chosen from
  the start. That correction became the more interesting half of item 12 — the first measured
  consequence of §5.2's absolute tolerances, which §17 had declined to change on the argument that the
  consequence was theoretical.

**Its negative results are recorded too, because they close questions.** It could not construct a
`Propensity` that G1 misses — superset, subset, permutation and dtype all caught, and a superset is
benign — so G1 is complete for §10 and stricter than it needs to be. And over some 6,000 fits it found
no input on which `polr` returns a non-ascending `alpha`, stops at a non-stationary point, or produces a
materially wrong number. It reproduced §20c's derivative measurements independently to three digits.

**Three places a third round should start**, and none of them is where either of the first two did:
`_ord_score_hess`'s `np.add.at` accumulation on a design with repeated response levels across arms;
`_record_outcome_completeness`'s `set`-and-`notna` pair on a frame with duplicate `case_id` values,
which Stage 2's A2 forbids but this stage does not re-check; and whether §14.12's surviving property
assertions can catch a regression that §14.0.2's five-record golden vector cannot — which is §15's first
bullet, still open, and now the largest unresolved question in the document.
