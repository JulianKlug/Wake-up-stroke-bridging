# Stage 10 spec — bootstrap engine

Implements roadmap Stage 10 [§10]. Section references in brackets are to `statistical_analysis_plan.md`.
Numbers and decisions referenced as DECISION *n* are established in Stage 0 and recorded in
`../out/stage0_data_inventory.md`. Constants named in `SMALL_CAPS` are Stage 1's and live in
`config.py`; `stage1_config_and_data_contract.md` is their specification. The frame this stage
resamples is specified by `stage5_cohort_construction.md` §11, the result object it refits by
`stage6_propensity_and_weights.md` §11, and the two estimators it refits by
`stage8_primary_outcome_estimator.md` §11 and `stage9_secondary_binary_estimators.md` §14 — the
latter of which is written as a handover *to this document* and is answered item by item in §14.

**Status.** Written 2026-08-25 against the landed Stages 1-9. **Twelve probes were run before any
section was drafted, and four of them changed what this document says** — §5.4 (the resampler blocker
Stage 9 handed over is one of *two*, and the second is not a resampler problem at all), §7.4 (the
`len(fit.alpha)` distribution Stage 8 asked for is **not** degenerate, and the answer arrives with a
second symptom Stage 8 predicted and could not measure), §8.2 (`numpy`'s default percentile definition
breaks [§10]'s own claim that the p-value agrees with the interval) and §10.2 (Stage 9's `max|beta|`
tail is real, reproducible, and attributed to an outcome the specified estimator fits no model for).
§21 names each. **Two of its decisions were put to the PI and confirmed on 2026-08-25 — DECISION 6
(S8 becomes a `FitError`) and DECISION 7 (per-outcome replicate sets) — and the same round answered a
standing open decision: the drop rate a penalised ordinal fit would reduce is measured at zero.**
Every number below was produced by running code; none is carried — where a Stage 8 or
Stage 9 number is quoted for comparison it is labelled as theirs and re-measured here. It is the **sole
source for the Stage 10 implementation**: everything the implementer needs is here, and anything not
here is not to be invented.

**Goal.** One interval and, where [§10] prescribes one, one p-value for every point estimate Stages 8
and 9 produce: the primary common odds ratio, the six cumulative `RD_k`, and the seven binary
outcomes' risk differences, marginal odds ratios and model-assisted augmented risk differences.
Patient-level nonparametric bootstrap, stratified by centre, `N_BOOT` replicates, recorded seed,
percentile limits. Plus the three diagnostics Stages 8 and 9 asked this stage to report because they
are answerable only from replicates.

**Not in scope.** The Benjamini-Hochberg correction within the two families, the subgroup estimates and
the E-value (Stage 11 [§13]); the standardisation analyses, which resample the same way and refit
something else (Stages 12 and 13 [§14a, §14b]); the reporting layer and every label it owes
(Stage 14 [§16]); and any second interval — BCa, bootstrap-t, an analytic standard error — which [§10]
does not prescribe and Stage 8 §5.6 establishes cannot be computed from the fitter's own information
matrix (§18).

**One thing this spec settles that no earlier stage could.** Every stage from 6 onward was written
against one frame — the workbook's cohort — and prespecified its behaviour on the population of frames
[§10] would feed it. This is the first stage that produces that population, so it is the first that can
check the prespecification against it, and **two of the four preconditions that were written as belts
over Stage 2's braces turn out to be assertions about the workbook rather than about the data**:
`propensity._record_exclusion`'s identifier count (Stage 9 §12.3 found this and handed it over) and
`outcome._assert_estimable`'s S8, which nobody had found because no frame before this one had a
constant outcome on it. §5 settles both, and the second is not the kind of problem the first is.

**And one finding that arrives with the stage rather than being designed into it.** Stage 8 §11 asked
Stage 10 to report the distribution of `len(fit.alpha)` across replicates, and said what turned on it:
*"if that distribution is degenerate at 6 the gap is theoretical, and if it is not, the drop rate is
partly a function of a bound calibrated on frames unlike the ones being dropped."* Measured: **it is
not degenerate — 4.6% of replicates fit five cutpoints rather than six** (§7.4). The drop rate happens
to be zero, so the specific consequence Stage 8 feared does not arise; but the same measurement carries
a second symptom Stage 8 predicted from a different direction and could not test, and that one *is*
present. Stage 8 §11 wrote that a replicate's `Σw` "is not measured here", that nothing was expected to
move, and that if it did *"the visible symptom would be `iterations` dropping rather than anything
failing"*. Measured: `Σw` runs from **4.90** to 37.48 against the point estimate's 27.74, and
`polr`'s iteration count runs from **1** against Stage 8's measured 4 to 5 on every frame it had. The
symptom is present, in the direction predicted, and §7.5 says what it does and does not imply.

---

## 0. Where Stage 10 sits

```
  cohort.build(...)          -> DataFrame, 93 x 33, one row per patient      [Stage 5 §11]
  propensity.fit(cohort, a)  -> Propensity: e, w, in_model (92)              [Stage 6 §11]
  outcome.primary(...)       -> Primary:   beta, alpha, rd, cumulative       [Stage 8 §11]
  outcome.secondary(...)     -> Secondary: seven BinaryEstimate              [Stage 9 §11]
  balance.assess(...)        -> Balance   — NOT read here                    [Stage 7 §11]

  +------------------------------------------------------------------------------------+
  |  bootstrap.run(cohort, ps, est, sec, audit)                                        |
  |                                                                                    |
  |    _assert_run_inputs(cohort, ps, est, sec)      R1-R7, collected         (§4.4)   |
  |    paths = {k: e.augmented_path for k, e in sec.estimates.items()}        (§6.3)   |
  |                                    READ off the point estimate, never re-derived   |
  |                                                                                    |
  |    for b in range(N_BOOT):              rng = default_rng(SEED), ONE stream (§4.3) |
  |        draw = resample(cohort, rng, BOOT_STRATUM)   renames each row      (§5.3)   |
  |        _replicate(draw, paths, source)          the body, below          (§6.1)    |
  |          ├─ propensity.fit(draw, throwaway Audit)                        (§6.2)   |
  |          ├─ outcome.primary(...)      -> beta, six RD_k, len(alpha)      (§6.1)   |
  |          └─ per outcome, in BINARY_OUTCOMES order:                       (§6.4)   |
  |                 rd, odds_ratio, augmented, max_abs_beta                           |
  |                 each independently droppable                            (§7.3)   |
  |                                                                                    |
  |    per estimand:  draws -> percentile_ci(...)                            (§8)     |
  |                          -> bootstrap_p(...)  where [§10] prescribes one (§9)     |
  |    _record_replicates(...)   one entry: counters + three distributions   (§11)    |
  +------------------------------------------------------------------------------------+

                       -> Bootstrap: draws, intervals, diagnostics             (§3.1)

  Stage 11 [§13] reads the p-values and corrects within families              (§12.1)
  Stages 12, 13 [§14a, §14b] reuse resample / replicates / percentile_ci      (§12.2)
  Stage 14 [§16] reports every interval with its surviving-draw count         (§12.3)
```

### 0.1 Why this is a new module and not an extension of `outcome.py`

Stage 9 §0.1 extended `outcome.py` rather than adding a module, and the reason it gave was that the
weighted per-arm proportion should have one home. The opposite conclusion holds here, for a reason that
is about the import graph and not about taste.

Stage 10 **refits Stages 6, 8 and 9**. It calls `propensity.fit`, `outcome.primary` and the Stage 9
estimators, so it must import all three. `outcome.py` already imports `propensity` and `model`; a
bootstrap living inside `outcome.py` would be a function in a module that imports two of the three
things it drives and is itself driven by nothing — which is fine — but it would also put the replicate
loop in the same module as the estimators it resamples, and that is the arrangement in which a future
edit can reach from the loop into an estimator's private. `_augmentable` and `_augmented_path` are
`outcome.py` privates, and §6.3's whole rule is that Stage 10 **reads the decided path and never
re-derives it**. A module boundary is what makes "never re-derives it" checkable by scan (§15.11)
rather than by review.

The second reason is Stage 12. [§14a] and [§14b] each prescribe *"patient-level bootstrap stratified by
centre"* over machinery that fits no propensity model at all, so they need this stage's resampler and
its percentile machinery and none of its estimators. A `bootstrap.py` whose three general functions
take the replicate body as a callable serves them; a bootstrap inside `outcome.py` would be imported by
Stage 12 for its resampler and would drag the [§7] estimators in behind it. **This is Stage 8 §0.1's
own move**: `model.polr` was made general in its covariates from the first line because [§14a]
prescribes the same estimator with a wider design, so that Stage 12 reuses it rather than a second
proportional-odds implementation existing. Here the reused thing is the resampling and the percentile
definition, and the argument is identical.

### 0.2 What Stage 10 does not touch

- **`Balance`.** Not read. Stage 7's diagnostics are not refit and no interval is put on an SMD.
  [§9] is a statement about the realised sample and DECISION 4 requires the exceedance reported beside
  the primary estimate — which is Stage 14's job with Stage 7's number, not a bootstrap quantity.
  §15.11 asserts `balance` is not imported, for Stage 7 §14.11's reason: not because importing it would
  fail, but because it would succeed.
- **The estimators.** `model.polr`, `model.firth`, `model.design`, `model.predict`,
  `outcome.weighted_proportion`, `outcome.cumulative_rd`, `outcome.weighted_rd`,
  `outcome.marginal_odds_ratio`, `outcome.augmented_rd`, `outcome.outcome_model` and
  `propensity.fit` are **called and not edited**. §15.12 asserts the Stages 1-9 suite passes unedited.
- **`propensity.py`.** Not amended, and that is the point of §5.3 rather than an accident of it. The
  guard Stage 9 §12.3 identified keeps its exact meaning and its full strength on the point estimate;
  what changes is the resampler, which is this stage's own code (§13).
- **The point estimate.** `run` reads `Primary` and `Secondary`; it does not recompute them, and every
  interval is centred on nothing — a percentile interval is a quantile of the replicate distribution
  and is not built around the point estimate. §15.8 asserts the limits are order statistics of the
  draws and are not a function of the point value.

---

## 1. Deliverables

| Path | Contents |
|---|---|
| `extended_bridging/bootstrap.py` | **new.** `Draws`, `Interval`, `Diagnostics`, `Bootstrap`; `resample`, `replicates`, `percentile_ci`, `bootstrap_p`, `run`; and **the privates §3.2 lists and counts** — the count is stated there and nowhere else |
| `extended_bridging/config.py` | **amended.** `CI_LEVEL`, `BOOT_STRATUM`, `PERCENTILE_METHOD`, `FAILURE_BUCKETS`, and `ci_min_draws()`; each with its comment block, which is shipping text and not commentary (§13) |
| `extended_bridging/outcome.py` | **amended, and these are the two Stage 9 changes this stage's requirement licenses.** S8 becomes `FitError` (§5.4); `secondary` gains `collect` so a per-outcome failure is captured rather than raised past the other six (§7.3). No other function changes |
| `extended_bridging/tests/fixtures_stage10.py` | **new.** §15.0's constructions as code — **six functions and five constants**, enumerated in §13 and counted there and nowhere else, for §3.2's reason |
| `extended_bridging/tests/test_bootstrap.py` | **new.** The §15 sections, banner-commented `# --- 15.x` |
| `extended_bridging/tests/test_outcome.py` | **amended.** S8's reclassification and `collect`'s branches (§15.5, §15.7) |
| `extended_bridging/tests/test_config.py` | **amended.** Six assertions (§13) |
| `extended_bridging/implementation_roadmap.md` | **amended.** Stage 10 gains its `**Spec:**` line, five corrections and two additions (§22) |
| `TODOS.md` | **amended.** Two items close, three are answered and rewritten, two are new (§16) |

`out/logs/audit_<label>.md` is an **output, not a deliverable**, and for this stage it is the *only*
place any interval exists (§4.5). It is gitignored.

**Nothing under `specs/` may quote a case identifier, and nothing here does.** Stage 8 §4.3 set that
rule and §4.5 states why it binds harder at this stage than at any before it.

---

## 2. Environment

Unchanged from Stage 9: `uv`, Python 3.12.12, pandas 2.3.3, numpy 1.26.4, statsmodels 0.14.6. `scipy`
is test-only by policy. **No dependency is added, and in particular no parallelism library**: §18 gives
the reason and §16 files the trigger under which that is revisited.

```bash
cd extended_bridging && uv sync && uv run pytest -v
```

**Cost, measured (§21).** Medians over 60 replicates, per replicate, on the cohort's shape:

```
  propensity.fit                                     13.4 ms   [Stage 9 §12.1, re-measured]
  outcome.primary                                     4.7 ms   [Stage 9 §12.1, re-measured]
  the seven binaries, whole `secondary`              69.9 ms   -> 139.9 s over N_BOOT
  the seven binaries, per outcome (§7.3's route)      61.5 ms   -> 123.0 s over N_BOOT
```

**The per-outcome route is CHEAPER than calling `secondary` whole, and that was not the expected
result.** §7.3 chooses per-outcome replicate sets on an inference argument — a failure in one outcome
must not thin another's draws — and the cost of that choice was expected to be the price paid for it.
Measured, it is negative: the per-outcome route skips the three audit tables `secondary` renders on
every call, which at `N_BOOT` are 6000 entries nobody reads (§6.2). So §7.3's decision costs **−16.9 s**
and needs no cost argument at all.

**`model.design` is the one worthwhile optimisation and it is worth 27.2 s.** The four full-list
outcomes share an identical design matrix, so building it once per replicate rather than four times
takes 18.12 ms to 4.50 ms. Stage 9 §2 estimated this at "roughly 22 s"; measured over `N_BOOT` it is
**27.2 s**. §6.4 states the condition under which the four designs are equal, which Stage 9 asserted
without one, and §15.9 asserts it rather than assuming it.

**Total, serial, single-threaded: about 145 s** for the whole engine with §6.4's optimisation in place.
That is small enough that §18's refusal to parallelise costs nothing worth having.

---

## 3. Module shape

### 3.1 What the stage returns

```python
@dataclass(frozen=True)
class Draws:
    """One estimand's replicate draws, and the counters that explain its denominator.

    **`draws` and `n_attempted` are BOTH fields and neither is derivable from the other**, which is
    the whole of [§10]'s "dropped and counted": `len(draws)` is what the percentiles were taken over
    and `n_attempted` is what was asked for, and a reader given only the first cannot tell a
    complete interval from a thinned one. Stage 14 prints both beside every interval (§12.3).

    `failures` is keyed by BUCKET and not by message, and §7.2 is why there is more than one key:
    Stage 8 §11 requires the separation count reported separately from the convergence count,
    because it measured that the second never fires on this estimator and a single counter reading
    zero therefore says nothing. The buckets are C.FAILURE_BUCKETS' values and the mapping from a
    raised FitError to one of them is §7.2's, asserted by scan in §15.6 rather than trusted.

    **The sum of `failures.values()` plus `len(draws)` equals `n_attempted`, and that is an
    assertion and not a comment** (§15.4). A replicate is either a draw or exactly one counted
    failure; an implementation that drops silently anywhere makes the three numbers stop reconciling,
    which is the one property of this dataclass that no wrong implementation satisfies by accident.
    """

    quantity: str                  # "beta", "rd_2", "sich.rd", ... — the estimand's key
    draws: np.ndarray              # (n,) float64, finite, in replicate order — §7.1
    n_attempted: int               # replicates in which this estimand was ATTEMPTED — §7.3
    failures: dict[str, int]       # bucket -> count; sums with len(draws) to n_attempted — §7.2


@dataclass(frozen=True)
class Interval:
    """A percentile interval, its level, and the p-value where [§10] prescribes one.

    **`p` is `None` where [§10] prescribes no test, and that is a field rather than an omission.**
    [§8] gives the six cumulative RD_k intervals and NO p-values, so `p is None` on all six, and
    Stage 14 printing a p-value for one of them is then a `None` reaching a formatter rather than a
    number nobody questioned (§9.4).

    **`method` travels with the interval.** §8.2 measured that numpy's default percentile definition
    disagrees with [§10]'s p-value at exactly the boundary that decides significance, so which
    definition produced these two numbers is part of what they are. It is pinned in `config.py`, so
    this field records rather than chooses.
    """

    lo: float
    hi: float
    level: float                   # C.CI_LEVEL
    method: str                    # C.PERCENTILE_METHOD — §8.2
    n_draws: int                   # the surviving draws these limits are order statistics of
    p: float | None                # §9; None where [§8] prescribes intervals only


@dataclass(frozen=True)
class Diagnostics:
    """The three things Stages 8 and 9 asked this stage to report, and one it asked for itself.

    None of these is an estimate and none gets an interval. They exist because a counter reading
    zero is not evidence of anything unless something beside it says what the replicates looked
    like — which is Stage 8 §11's argument for `len(fit.alpha)` and Stage 9 §9.6's for `max|beta|`,
    made twice for the same reason.

    `sum_w` is this stage's own addition and Stage 8 §11 is why: it wrote that a replicate's Sw was
    not measured there, that nothing was expected to move, and that the visible symptom if it did
    would be `polr`'s iteration count dropping. Both are measured in §7.5 and both moved, so the
    pair is reported rather than the reassurance.
    """

    n_alpha: dict[int, int]               # cutpoint count -> replicates — §7.4
    polr_iterations: dict[int, int]       # iteration count -> replicates — §7.5
    sum_w: np.ndarray                     # (n,) Sw over in_model per replicate — §7.5
    max_abs_beta: dict[str, np.ndarray]   # outcome -> max|beta| of m_a(X); AUGMENTED ONLY — §10.2
    n_in_model: dict[int, int]            # in_model size -> replicates — §7.5


@dataclass(frozen=True)
class Bootstrap:
    """Every interval [§10] prescribes, the draws behind them, and the diagnostics beside them.

    **The seed and the replicate count are fields because [§16] requires them reported**, and they
    are recorded here rather than looked up from `config.py` at print time: a run summary that reads
    its seed from the configuration reports the configuration and not the run.

    There is no `favours_bridging()` and no `significant` field, for Stage 8 §3's reason carried one
    stage on. [§10] closes with "estimation, not testing, is the reportable output"; a boolean here
    would be the dichotomisation that sentence declines, computed once and then quoted forever.
    """

    seed: int
    n_boot: int
    draws: dict[str, Draws]           # keyed as §3.3's estimand keys
    intervals: dict[str, Interval]    # the same keys
    diagnostics: Diagnostics
```

### 3.2 Public surface, and it is five names

```python
def resample(df: pd.DataFrame, rng: np.random.Generator, stratum: str) -> pd.DataFrame: ...
def replicates(df: pd.DataFrame, body: Callable[[pd.DataFrame], object],
               n: int, seed: int, stratum: str) -> tuple[object, ...]: ...
def percentile_ci(draws: np.ndarray, level: float) -> tuple[float, float]: ...
def bootstrap_p(draws: np.ndarray) -> float: ...
def run(df: pd.DataFrame, ps: propensity.Propensity, est: outcome.Primary,
        sec: outcome.Secondary, audit: Audit) -> Bootstrap: ...
```

The first four are **general and know nothing about [§7] or [§8]**: `resample` takes a frame and a
stratum column, `replicates` takes a callable, `percentile_ci` and `bootstrap_p` take an array. That is
§0.1's argument as a signature — Stages 12 and 13 call these four and supply their own body (§12.2).
`run` is the [§10] instantiation and is the only one of the five that names an estimator.

Privates in `bootstrap.py` are `_assert_run_inputs`, `_replicate`, `_bucket`, `_collect`,
`_estimand_keys`, `_shared_design`, `_counters_table`, `_diagnostics_table`, `_replicates_detail`,
`_record_replicates` — **ten**, against the fourteen Stage 9 added. **This count is stated once, here,
and §1 and §13 cite this sentence rather than repeating the number**, because Stage 9 §22.3 item 5
records three drafts carrying three different counts of its own.

### 3.3 The numerical facts this stage turns on

```
  * np.percentile's DEFAULT method is "linear", and it interpolates between order statistics.
    [§10] claims the p-value "is the smallest level at which the percentile interval for beta
    excludes the null" and "agrees by construction with the reported interval". THAT CLAIM IS FALSE
    UNDER THE DEFAULT. Measured over 30000 constructed draw sets at B = 2000: "linear" disagrees in
    4.707%, "higher" in 4.707%, "nearest" in 4.707%, "midpoint" in 0.117%, and "inverted_cdf" and
    "lower" in 0.000%. Every disagreement is at a smaller-tail count of exactly 50, where p is
    exactly 0.0500 and the interpolated limit lands on the far side of zero (§8.2).

  * SO THE METHOD IS PINNED IN config.py AND IS NOT A DEFAULT. `inverted_cdf` is the empirical-CDF
    quantile — the smallest order statistic whose cumulative weight reaches the level — which is the
    same object the p-value's Pr(beta* <= 0) is computed from, and that identity is why the two agree.

  * A REPLICATE'S ROWS ARE DUPLICATED BY CONSTRUCTION and two landed preconditions read that as a
    defect. propensity._record_exclusion (propensity.py:400-432) compares DISTINCT excluded case_ids
    against excluded ROWS; outcome._assert_estimable's S8 (outcome.py:924-930) raises when an
    outcome is constant on its population. The first is about names and is the resampler's to fix;
    THE SECOND IS ABOUT DATA AND IS NOT (§5).

  * model.FitError CARRIES NO CODE. It is `class FitError(RuntimeError)` with no attributes
    (model.py:89), and all SIXTEEN raise sites — fourteen in model.py, two in outcome.py — distinguish
    themselves by the FIRST TOKEN OF THE MESSAGE: G6, G7, O1-O6, F3, F4, F6, "polr:", "Firth:".
    Stage 8 §11 requires the separation count reported SEPARATELY from the convergence count, so
    Stage 10 must classify, and classification is by prefix (§7.2).

  * EVERY QUANTITY DERIVED FROM CI_LEVEL BY FLOAT ARITHMETIC MUST BE SNAPPED, BECAUSE `1.0 - 0.95`
    IS `0.050000000000000044`. Two consequences, both measured, both of which broke a fence in this
    document before they were fixed (§21b): `100.0 * ((1.0 - 0.95) / 2.0)` is `2.500000000000002`
    and not `2.5`, which at the §8.2 boundary moves the order-statistic index up by one and FLIPS
    THE LIMIT'S SIGN -- measured, `-0.107095` at the literal 2.5 against `+0.101248` at the computed
    one, on the same draws; and `int(np.ceil(2.0 / (1.0 - 0.90)))` is 21 where 20 is intended.
    Neither is a rounding cosmetic: the first destroys the agreement §8.2 exists to guarantee, at
    exactly the tail count the pin was chosen for.

  * AND A PREFIX MAP IS ONLY SAFE IF SOMETHING FAILS WHEN A MESSAGE CHANGES. §15.6 scans both
    modules for FitError raise sites and asserts every leading token is in C.FAILURE_BUCKETS, so a
    reworded message is a test failure and not a silently misfiled counter.

  * S7 — an arm carrying zero total weight — DID NOT FIRE IN N_BOOT REPLICATES. Stage 9 §4.4 argues
    S7 means the mask and the weights disagree about the same rows, which is a bug rather than a
    sparse replicate, and §5.4 therefore reclassifies S8 ALONE. That argument is untouched by this
    measurement: 0 in 2000 is consistent with "it means a bug" and is not evidence for it. What the
    measurement does establish is that S7 is not a second blocker, which is the question that had to
    be asked once S8 turned out to be one.

  * A DUPLICATED INDEX PASSES THROUGH THE WHOLE PIPELINE UNCHANGED, AND AN EARLIER DRAFT OF THIS
    BULLET SAID IT DOES NOT. `.iloc` with a repeated integer array returns repeated rows and a
    repeated index, and the draft asserted that `ps.w.loc[mask]` in Stage 9's estimators would then
    return more rows than the mask selected. Measured on a replicate with 59 distinct labels over 93
    rows: `propensity.fit`, `outcome.primary` and `outcome.secondary` all run, `len(ps.w.loc[mask])`
    is 93 and not more, and every count is right — because those reads are BOOLEAN and boolean `.loc`
    is positional in effect, not label-based. §5.3 resets the index anyway and the reason it gives is
    not this one.

```

---

## 4. The resample [§10]

### 4.1 The rule, and the design [§10] declines

[§10], in full on this point: *"Nonparametric bootstrap, 2000 replicates, stratified by centre,
refitting the propensity model in every replicate. Percentile 95% confidence intervals; seed recorded.
Replicates whose prespecified fit fails are dropped and counted, never substituted."* And:
*"Resampling is at patient level within centre. A centre-level cluster bootstrap is not used: a handful
of centres cannot support cluster-bootstrap consistency. **Inference is conditional on the participating
centres and their observed treatment practices.**"*

Three things follow that this stage does not get to reconsider. The resampling unit is the patient. The
stratification is by centre and the centre totals are therefore **fixed across replicates** — which is
what makes the inference conditional on these centres, and is why a replicate cannot lose a centre.
And a failed replicate is dropped, never replaced by a different estimator: §7 is that sentence as
code, and roadmap invariant 5 — *"no estimate is ever produced by a fallback estimator"* — is what
§15.10 asserts.

### 4.2 The strata, and the one that decides §5

```
  stratum   n     bridged / EVT alone     note
  HUG      41     30 / 11                 [§13 amendment: treatment is nearly determined by centre]
  Lugano   31      2 / 29                 carries the ONE covariate-incomplete record  -> §5.1
  CHUV     21                             
  ---------------------------------------------------------------------------------------
  total    93                             USZ contributes none: never-IVT, removed at Stage 5
```

Lugano's size is the number §5.1 turns on: the analytic probability that a *given* row is drawn at
least twice from a stratum of 31 is **0.2642**, and that is what the blocker's measured 26.0% rate is.
The HUG and Lugano treatment splits are quoted from [§13]'s amendment rather than re-derived; they are
here because they are why no replicate may lose a centre.

### 4.3 The seed, one stream, and what determinism does and does not mean

`C.SEED` = 20260807 and `C.N_BOOT` = 2000 are already declared (`config.py:269-270`) and are read, not
written. **One `numpy.random.Generator` is created once and consumed by every replicate in order**,
rather than one seeded per replicate: a per-replicate seed derived from `b` is reproducible too, but it
makes the draw a function of an index that a later edit to the loop can renumber, and the failure is
silent. Measured: the same seed returns a bit-identical frame and a different seed does not (§21).

`np.random.default_rng`, not `RandomState`. The legacy class is not deprecated but its stream is a
compatibility guarantee rather than a design, and nothing here needs the old stream.

**What determinism does not mean.** `resample` ranges over the stratum's **sorted unique values**, so
the draw is a function of the seed *and of the stratum label alphabet*. Renaming a centre — which
`C.CENTER_RECODE` does — would consume the stream in a different order and produce a different, equally
valid set of replicates. That is stated rather than guarded against, because the alternative is to
range over `C.CENTER_ORDER`, and `resample` is general (§12.2): Stage 12 stratifies the same way over
a frame this stage never sees, and a general function cannot name this study's centres. §16 files it.

### 4.4 The preconditions, in two phases

Stage 9 §4.4a's boundary, and its wording is the rule rather than a paraphrase: *"the boundary is
can-this-be-READ vs is-the-DATA-judgeable — NOT mask-checks vs the-rest."* A column-presence check is a
precondition of every read that follows it.

```
  phase 1 — can this be read at all; collected, all failures reported at once
  R1  BOOT_STRATUM is a column of `df`
  R2  `case_id` is a column of `df`                                   -- §5.3 rewrites it
  R3  ps.in_model, ps.e, ps.w are indexed like `df`: SAME LABELS, SAME ORDER   [Stage 9 §12]
  R4  ps.in_model is boolean and total
  R5  sec.estimates' keys are exactly C.BINARY_OUTCOMES               -- so §6.3's map is complete

  phase 2 — is the data resamplable; collected
  R6  no stratum is empty, and the strata partition `df` with no missing label
  R7  every sec.estimates[k].augmented_path is one of the three declared strings
  R8  C.N_BOOT >= the §8.3 floor, so a percentile limit is an order statistic
```

**R1-R5 are `SchemaError` and R6-R8 are `SchemaError` too, and this stage adds no `FitError` of its
own.** That asymmetry is deliberate and it is the opposite of the one Stage 9 needed: every condition
above is a property of the *call*, not of a replicate, so all eight are the caller's bug. The only
`FitError` this stage ever sees is one raised by an estimator inside a replicate, and §7.1 is what it
does with it.

A stratum of size 1 satisfies R6 and contributes no variability — it resamples to itself in every
replicate. That is not an error and is not guarded: it is what stratifying on a near-determining
variable means, and [§10] chose it knowing so. `USZ` is not such a stratum; it is absent, because
Stage 5 removed the never-IVT centre from the [§7] cohort entirely.

### 4.5 Why no interval is in this document

Stage 8 §4.3 established the rule and gave the reason: every stage before it recorded its decisions as
taken before any outcome was examined by arm, Stage 8 is where that necessarily ends, and it ends after
the specification is committed rather than during it. The regression pin is a synthetic golden vector,
the estimate lands in the gitignored log, and the cost — no pinned regression number on the primary
effect anywhere in git — is stated rather than minimised.

**At this stage the rule binds harder, because this is where the number a manuscript prints comes into
existence.** A `β` in a committed document is a point estimate somebody could have chosen a
specification to obtain; an *interval* is that plus the claim of significance, and [§16] requires — as
DECISION 5's amendment puts it — that *"no number on the primary effect be pinned before the analysis
is locked"*. So:

- **No `β`, no `exp(β)`, no `RD_k`, no `rd`, no `odds_ratio`, no `tau`, no interval limit and no
  p-value from the workbook appears anywhere below.** Not in prose, not in a fence, not in §21.
- **Rates, counts and distributional shape do appear**, and they are what every measured claim here is
  made of: how often a guard fired, how many cutpoints a replicate fitted, how the nuisance
  coefficient's magnitude is distributed, how long the loop took. Stage 8 §4.3's line is *"a count may
  appear here; a weighted quantity may not"*, and `Σw` and `max|beta|` are admitted under it on the
  same footing Stage 8 §5.2 and Stage 9 §9.6 admitted theirs: neither is an effect estimate.
- **Every regression pin is measured on a synthetic fixture** (§15.0), and the acceptance criteria
  assert *properties* — the counters reconcile, the limits are order statistics, the p agrees with the
  interval, coverage is near nominal on data whose truth is known.

---

## 5. A resampler the pipeline does not raise on

### 5.1 The blocker Stage 9 handed over, reproduced

Stage 9 §12.3 measured it and could not fix it. Reproduced here as the first probe, because a handover
that is not re-measured is a handover taken on trust:

```
  400 stratified replicates, the [§10] draw with NO renaming, through propensity.fit:

    SchemaError                                              104   26.0%
      "covariate_completeness excludes 2 record(s) and can name 1."   75
      "covariate_completeness excludes 3 record(s) and can name 1."   22
      "covariate_completeness excludes 4 record(s) and can name 1."    7
    ok                                                       296   74.0%

  analytic P(a given row drawn >= 2x from a stratum of 31)         0.2642
```

104 of 400 is 26.0%, which is Stage 9's number to the decimal, and the analytic 0.2642 identifies the
cause exactly: the cohort has exactly **one** covariate-incomplete record, it is in Lugano, and
`_record_exclusion` (`propensity.py:400-432`) compares the number of **distinct** excluded `case_id`s
against the number of excluded **rows**. Draw that row twice and two excluded rows collapse to one
name.

The guard's own docstring states the premise it rests on: *"A duplicated or missing case_id is
forbidden by Stage 2's A2, so this is a belt over those braces."* A2 forbids duplicates **in the
workbook**. A patient-level bootstrap replicate contains duplicates by construction — that is what
sampling with replacement *is* — so the premise does not hold for the population of frames [§10] feeds
this function, and the guard is measuring the resampler rather than the data.

And Stage 8 §11 prespecified that Stage 10 may not catch it: *"`FitError` is the droppable failure and
`SchemaError` is not — a `SchemaError` from `primary` or `propensity.fit` is a bug in the resampler,
not a sparse replicate, and catching it would drop replicates for a reason that is not about the data."*

### 5.2 The two candidates, and why the choice is not a toss-up

`TODOS.md` named two and said choosing between them *"is a decision about what a replicate's log
means"*. Both work; they differ in what they cost.

| | Candidate 1 — the resampler renames each drawn row | Candidate 2 — the guard learns about resamples |
|---|---|---|
| Where the change lands | `bootstrap.resample`, this stage's own new code | `propensity._record_exclusion`, a Stage 6 amendment |
| The guard on the point estimate | **unchanged, full strength** | weakened: it can no longer compare a count against a count |
| What a replicate's `case_ids` mean | distinct draw names, `patient#2` | patient names, with the duplicate silently tolerated |
| Blast radius | one function nobody else calls | a function on the path of every `propensity.fit` call ever made |

**Candidate 1, and the deciding argument is that the guard is doing no work in a replicate anyway.**
Its purpose, stated in its own message, is that *"an estimate whose denominator the log cannot
reconstruct is an estimate nobody can check [§11]"* — a statement about a log somebody reads. A
replicate's `Audit` is a throwaway that is never written (§6.2), so in a replicate the guard protects a
document that does not exist. Weakening it there would weaken it on the point estimate too, where it
protects a document that does. Renaming leaves it exactly as strong where it matters and satisfies it
truthfully where it does not: the k-th draw of a patient **is** a distinct row of the replicate, and
naming it distinctly is a true statement about the frame rather than a concession.

**PI decision, 2026-08-25**, recorded in §17 with its reversible alternative, because neither [§10] nor
any amendment specifies how a replicate's identifiers are formed.

Measured, after the change: **0 of 400** raises, and **0 of `N_BOOT` = 2000** in the full loop. The
replicate carries 93 rows and 93 distinct `case_id`s, and `is_unique` is `True` (§21).

### 5.3 `resample`, written out

```python
def resample(df: pd.DataFrame, rng: np.random.Generator, stratum: str) -> pd.DataFrame:
    """One [§10] replicate: patient-level draw with replacement WITHIN each stratum.

    General in its frame and its stratum from the first line, because [§14a] and [§14b] prescribe
    the same resampling over a population this stage never sees (§12.2). It names no covariate, no
    outcome and no centre, and `stratum` is a parameter rather than `"center"` for that reason.

    Three things about it, each of which is a failure if changed:

      * THE STRATUM TOTALS ARE FIXED. `size=len(g)` per stratum, never a draw over the whole frame,
        because [§10]'s inference is "conditional on the participating centres" and a replicate that
        lost a centre would not be (§4.1).
      * THE ORDER IS `sorted(unique)`, NOT frame order and not `groupby`'s default. The draw is a
        function of the seed AND of the order the strata are visited in; frame order would make it a
        function of how the workbook happened to be sorted, which is the property Stage 2's shuffle
        test exists to deny (§4.3).
      * EACH DRAWN ROW GETS A DISTINCT case_id, and §5.2 is the whole argument. The k-th appearance
        of a patient becomes `case_id#k` counting from 1, so the FIRST appearance is `id#1` and not
        the bare id -- uniform, because a mixed scheme makes "was this row drawn once" a question
        about string formatting.

    The index is reset before the rename. **This is NOT because a duplicated index breaks anything**
    -- measured, `propensity.fit`, `outcome.primary` and `outcome.secondary` all run correctly on a
    replicate with 59 distinct labels over 93 rows, because every read of `e` and `w` is boolean
    (§3.3). It is reset because this function is general and Stages 12 and 13's bodies are unwritten,
    and a frame with a duplicated index is a frame on which a future label-based `.loc` is wrong in a
    way that returns a number.

    Deterministic given `rng`'s state. Adds no column, and the frame it returns has the same columns
    and dtypes as the frame it was given -- asserted in §15.3, because a resampler that quietly
    changes a dtype changes `model.design` (TODOS' Stage 1 dtype item).
    """
    parts = [
        g.iloc[rng.integers(0, len(g), size=len(g))]
        for _, g in sorted(df.groupby(stratum, sort=True, observed=True), key=lambda kv: str(kv[0]))
    ]
    out = pd.concat(parts, axis=0).reset_index(drop=True)
    occurrence = out.groupby("case_id", sort=False, observed=True).cumcount() + 1
    out["case_id"] = [f"{cid}#{k}" for cid, k in zip(out["case_id"], occurrence)]
    return out
```

### 5.4 The second blocker, and it is not a resampler problem

**The renaming fix does not make the pipeline resamplable, because there is a second `SchemaError` on
the path and it is about the data rather than about the names.** Nobody had found it, because no frame
before this one had a constant outcome on it.

`outcome._assert_estimable`'s S8 (`outcome.py:924-930`) raises when a binary outcome takes one value on
its [§11] population:

> `S8  {key}: constant at {y[0]!r} on its [§11] population of {y.size}. The risk difference is 0 and
> both weighted proportions are degenerate in the SAME direction, so §6.3's correction returns a finite
> odds ratio near 1 for an outcome with no contrast in it — a plausible number rather than a legible
> failure [Stage 9 §6.2, §16 item 3].`

Measured over 1998 live replicates:

```
  outcome        S8 raises      rate
  sich                  16      0.8%
  tici_2b_3              3      0.2%
  every other            0      0.0%
  ------------------------------------
  replicates affected   19      1.0%
```

`sich` has 5 events in 92, so a stratified resample drawing none of them is neither rare nor a bug. And
by the rule immediately above, **Stage 10 may not catch it** — so as landed, the three options are the
same three Stage 9 enumerated for the first blocker, and all three are still prohibited.

**S8 becomes `FitError` — DECISION 6 (PI, 2026-08-25) — and the argument for it is one Stage 9 wrote
itself.** `_assert_estimable`'s
docstring already names this as open: *"S8 is `SchemaError` and §16 item 3 is the open question about
whether it should be FitError"* (`outcome.py:897`). `TODOS.md` states the criterion and the trigger —
*"The first Stage 10 replicate that hits it. If the rate is non-negligible, `FitError` is almost
certainly right and the change is one line plus Stage 9 §4.4 and §15.1."* The rate is 0.8% on one
outcome and 0.2% on another; the trigger has fired and the criterion is met.

It is also the right classification on the merits, and Stage 9's own reason for the other choice is
what shows it. Stage 9 chose `SchemaError` because a constant outcome yields *"a plausible number
rather than a legible failure"* — the objection is to **returning** a number, and `FitError` does not
return one. It drops the replicate and counts it, which is [§10]'s mechanism for exactly this, and
which makes the count visible in the log where `SchemaError` makes it a crash. The distinction
`SchemaError` is for is *"a bug in the resampler, not a sparse replicate"*, and a safety outcome with
five events in ninety-two having none in a resample is the definition of a sparse replicate.

**S6 and S7 stay `SchemaError` and only S8 moves.** S6 is an empty population and S7 is an arm with no
record or no weight; Stage 9 §4.4 argues both mean the mask and the weights disagree about the same
rows, which is a bug. Neither fired in 2000 replicates, which is not evidence for that argument but does
establish that neither is a third blocker (§3.3).

### 5.5 What stays unchanged in Stage 6, and that is the point

`propensity.py` is **not amended**. `_record_exclusion` keeps both of its mechanisms — the `set` that
collapses a duplicate and the `notna()` that drops a missing id — and keeps comparing a name count
against a row count, which on the workbook is the check its docstring says it is. Stage 9 §13 declined
to touch it on the grounds that *"that is a Stage 6 amendment on the strength of a Stage 10 requirement,
and it does not belong in this commit"*; this is that commit, and the requirement turns out not to need
the amendment.

One thing found while checking that claim and **not** fixed here: `propensity.py:413` cites
`data.py:255` for `Audit.record`'s identifier normalisation, and the normalisation is at `data.py:271`.
`data.py` is unchanged since the Stage 9 spec commit, so the citation drifted between Stage 6 and now.
It is a stale line reference in shipped code, it is in the one module this stage's central claim says is
untouched, and touching it would cost that claim for a comment. §16 files it with its trigger.

---

## 6. The replicate body

### 6.1 What is refit, and in what order

[§10]: *"refitting the propensity model in every replicate"*, and the roadmap's Stage 10 entry expands
it: *"In every replicate refit the propensity model *and* the outcome regression, recompute weights,
recompute every estimate."* So the whole of Stages 6, 8 and 9 runs inside the loop; nothing is carried
in from the point estimate except §6.3's paths.

```
  _replicate(draw, paths, source):
      audit = Audit(source)                    throwaway, never written             (§6.2)
      ps    = propensity.fit(draw, audit)      refits e, w, in_model                [Stage 6]
      primary:                                 beta, six RD_k, len(alpha),
                                               polr iterations, Sw                  (§7.4, §7.5)
      per outcome in C.BINARY_OUTCOMES order:   rd, odds_ratio, augmented,
                                               max|beta| of m_a(X)                  (§7.3)
```

The order is the point estimate's order and it is not an optimisation target. `propensity.fit` comes
first because everything downstream reads `e`, `w` and `in_model`; `primary` before the binaries
because that is the order the point estimate ran in and a log that interleaves them differently is a
log that cannot be diffed against it.

**No estimate is computed from another replicate's fit.** `e` and `w` are the *replicate's*, refit from
the replicate's rows, and the one rule every consumer owes — Stage 6 §9's, restated by Stage 8 §11 and
Stage 9 §12 — is honoured unchanged: `e` and `w` carry `nan` off `in_model`, so range over the mask,
never over `notna()`, and never fill.

### 6.2 The throwaway `Audit`, and the `Source` it carries

Stage 8 §11 handed this over and gave the arithmetic: *"A replicate loop must pass its own `Audit`, and
'or accept 3 entries per replicate' understates it — at `N_BOOT` = 2000 that is 6000 entries, each
carrying a rendered table of 3 to 10 rows and a `detail` string of several hundred characters, and one
of the three carries a `case_ids` tuple. Reusing the point fit's `Audit` would make the log unreadable
and the process's memory a function of `N_BOOT`."* With Stage 9's three entries it is nine per
replicate, or 18 000.

`Audit.__init__` takes a `Source` (`data.py:259-264`), so a throwaway needs one, and **it is the point
estimate's own `audit.source` rather than a synthetic one.** A second module-level `Source` is pinned
against by `data.SOURCES` and a test, and inventing one here would either break that pin or lie about
provenance: the replicate *is* drawn from the workbook, and the label in a log nobody writes should
still be true.

A `record`-suppressing mode is not offered, for the reason Stage 8 gave: *"a stage that can be asked
not to log is a stage whose log is optional."* The throwaway is discarded when the replicate returns.

### 6.3 The frozen augmentation paths

Stage 9 §9.5 decided this and §11.1 built the parameter for it; Stage 10's job is to use it and not to
re-derive it. The rule: **the augmentation path is decided once, on the point estimate, and passed into
every replicate.**

```python
paths = {key: est.augmented_path for key, est in sec.estimates.items()}   # complete, never partial
```

Read off `Secondary`. Never computed from a minority cell, never from `_augmentable`, never from
`C.RARE_MINORITY_THRESHOLD`. §15.11 asserts by scan that `bootstrap.py` names none of those three.

Why it matters, measured on this run over 1998 live replicates:

```
  outcome        point minority   frozen path     replicates whose OWN cell
                                                  would decide differently
  ph2                         9   unaugmented     812      40.6%
  sich                        5   unaugmented      52       2.6%
  every other                     unchanged         0       0.0%
```

`ph2` sits one below `RARE_MINORITY_THRESHOLD` = 10 (Stage 9 §9.4 looked at that and left it alone), so
a resample crosses it four times in ten. A rule re-evaluated per replicate would make `ph2`'s
percentile interval a quantile over a near-even mixture of two different estimators — [§10]'s *"never
substituted with a different estimator"* violated by a route that raises nothing.

Stage 9 measured 46.1% for `ph2` and 3.8% for `sich`; this run measures 40.6% and 2.6%. Both are the
same phenomenon at the same order and neither supersedes the other: Stage 9's were taken over 400
replicates with its scratchpad renaming and a different stream, these over 1998 with the landed
`resample`. §21 records both and the difference is sampling, not disagreement.

### 6.4 One design per replicate, and the condition Stage 9 asserted without

Stage 9 §2 and §14 item 5: *"The four full-list outcomes share an **identical** design matrix, so
building it once per replicate is worth roughly 22 s."* Measured over `N_BOOT`, the saving is **27.2 s**
of about 145 s, which makes it the one optimisation worth taking.

**But "share an identical design matrix" is a claim with a condition, and Stage 9 stated it without
one.** `model.design` is called on `df.loc[in_estimate]`, and `in_estimate` is `ps.in_model &
df[key].notna()` — *per outcome*. Four outcomes share a design only if they are present on the same
rows. The condition holds, and it holds for a reason that is asserted elsewhere rather than by luck:

```
  outcome        source     C.outcome_model_covariates(...)
  mrs_0_2_90d    mrs_90d    the full [§6] list          ┐
  mrs_0_1_90d    mrs_90d    the full [§6] list          ├─ all four derived from mrs_90d
  death_90d      mrs_90d    the full [§6] list          │
  mrs_5_6_90d    mrs_90d    the full [§6] list          ┘
  tici_2b_3      None       center + atrial_fib   -- the [§8] override, a different list
  sich, ph2      None       unaugmented           -- no design is built at all
```

All four full-list outcomes have `source = "mrs_90d"` in the registry, and roadmap invariant 6 is
*"derived dichotomies carry exactly the missingness of their ordinal source"* — asserted in the test
suite since Stage 3. So their `notna()` masks are identical **by an invariant**, not by coincidence, and
the optimisation rests on something that fails loudly if it stops being true.

Measured anyway, because an invariant cited is not an invariant checked at this call site: over 300
replicates the four masks differ in **0**, and the four design matrices compare equal with
`DataFrame.equals`. §15.9 asserts both — the equality *and* a companion in which one outcome's
missingness is perturbed and the shared design is shown to be wrong, because an optimisation whose
precondition is never violated in a test is an optimisation nobody has tested.

`_shared_design` builds it once for the four and lets `tici_2b_3` build its own. It is a cache keyed by
nothing: one design, one replicate, discarded.

---

## 7. The failure taxonomy, and what a dropped replicate is

### 7.1 `FitError` is droppable and `SchemaError` is not

Handed over twice, by Stage 6 to Stage 8 and by Stage 8 §11 to here, and it is the rule the whole of §5
was about honouring:

- **`model.FitError` from inside a replicate is caught, the replicate's affected estimand is dropped,
  and the drop is counted in a named bucket.** This is [§10]'s *"replicates whose prespecified fit fails
  are dropped and counted"*.
- **`C.SchemaError` is never caught.** It propagates out of `run` and terminates the analysis. A
  `SchemaError` means a frame that cannot be read, which for a frame this stage constructed means this
  stage constructed it wrongly.
- **`assert not issubclass(model.FitError, C.SchemaError)`** is itself a test (§15.4), because the whole
  taxonomy is one `except` clause away from collapsing.

§5.4's reclassification of S8 is what makes this rule *satisfiable* rather than merely stated: before
it, 1.0% of replicates raised a `SchemaError` that no correct implementation could catch and no correct
resampler could prevent.

### 7.2 The buckets, and why prefix classification needs a scan test

Stage 8 §11 requires more than a count: *"Stage 10 should report the G7 count separately from the
convergence count, because the two mean different things about the data."* To report them separately,
Stage 10 must classify the `FitError` it caught — and **`FitError` carries nothing to classify by**. It
is `class FitError(RuntimeError)` with no attributes (`model.py:89`), and all **sixteen** raise sites —
fourteen in `model.py`, two in `outcome.py` — identify themselves by the first token of the message:

```
  bucket              codes                     raised by
  separation          G7                        outcome._assert_reportable        (outcome.py:346)
  constant_outcome    S8                        outcome._assert_estimable         (§5.4)
  nonconvergence      "polr:", "Firth:"         model.polr, model.firth
  degenerate_design   G6, O1-O6, F3, F4, F6     outcome._assert_exposure_survived (outcome.py:320)
                                                and model's fittability guards
```

**G6 is in that table because the scan found it and an earlier draft of this section did not.** It is
`_assert_exposure_survived` — the [§8] design lost its treatment column to `design`'s constant-column
drop — and it is a genuine sparse-replicate mode rather than a caller bug: a draw in which one arm
vanishes produces it. Stage 8 §11 names *"G6 and G7"* together as the two `FitError` routes on the
caller's side, and this document's first draft carried only G7 — which is exactly the omission §15.6's
scan exists to catch, caught by running it. It is bucketed as `degenerate_design` because it is
literally a statement about the design matrix, and because Stage 8 §11 asks for G7 against everything
else rather than for a bucket per code.

So `_bucket` reads the leading token and maps it, and `C.FAILURE_BUCKETS` is that map declared in one
place. **A prefix map is only safe if something fails when a message changes**, so §15.6 scans
`model.py` and `outcome.py` for `raise FitError(`/`raise model.FitError(` sites, extracts the first
token of each message, and asserts every one is a key of `C.FAILURE_BUCKETS`. A reworded message is
then a test failure rather than a counter that silently reads zero.

**The alternative is better and is deferred rather than adopted.** A `code` field on `FitError`, set at
each raise site, would make this structural instead of textual. It is not taken here because it is a
fourteen-site amendment to `model.py` — the module Stages 6, 8, 9 and 12 all depend on — on the strength
of a *diagnostic*, and because it would not remove the need for §15.6's scan so much as move what it
asserts. §16 files it with the trigger: the first time a bucket is wrong.

`"polr:"` covers two different failures — step-halving exhausted (`model.py:945`) and no convergence in
`POLR_MAX_ITER` (`model.py:975`) — and they land in one bucket. That is accepted rather than overlooked:
Stage 8 §11 asks for G7 against everything else, Stage 8 §5.2 records that exhausting the halvings *"is
not a convergence route"*, and both are measured below at zero.

### 7.3 Per-outcome replicate sets

[§10] says *"replicates whose prespecified fit fails are dropped and counted"* and does not say whether
a failure drops that replicate for every estimand or only for the one that failed. It is silent, so this
is a decision, and Stage 9 §14 constrains it in exactly one place: *"never take percentiles of `tau` and
of `rd` from different replicate sets, since the two are reported as a comparison."*

**Each estimand group keeps its own replicate set** — DECISION 7 (PI, 2026-08-25). The primary's `β`
and six `RD_k` are one group —
they come from one `polr` fit and a comparison across them is [§16]'s constant-shift statement, so they
must share draws. Each binary outcome's `rd`, `odds_ratio` and `augmented` are one group, which is Stage
9's rule satisfied exactly. Different outcomes do not share a set.

Measured over 1998 live replicates, this is what the choice is worth:

```
  calling `secondary` whole      19 replicates lose ALL SEVEN outcomes to one S8 raise
                                  1 replicate  loses ALL SEVEN to a Firth FitError on death_90d
                                 -> 20 replicates, 1.0%, x 7 outcomes

  per outcome                    sich  16 dropped   tici_2b_3  3 dropped
                                 death_90d 1 dropped   every other 0
                                 -> each outcome's set is thinned only by its own failures
```

**And it is cheaper, which was not expected.** §2: 61.5 ms per replicate against 69.9 ms for the whole-
`secondary` route, because the per-outcome route skips the three audit tables `secondary` renders on
every call. So the inference argument and the cost argument point the same way and there is no trade-off
to adjudicate.

**The mechanism is a Stage 9 amendment and not a re-implementation of Stage 9's loop.** Stage 9 §14 is
explicit that Stage 10's second option — driving the estimators individually — risks the thing §6.3
forbids: *"`_augmentable`'s rule, which it must **read** rather than re-implement, since
re-implementing it is how item 2 gets violated by accident."* Re-driving the seven-outcome loop in
`bootstrap.py` would put a copy of Stage 9's ordering, masking and path logic in a second module. So
instead `outcome.secondary` gains one parameter:

```python
def secondary(df, ps, audit, paths=None, collect=False) -> Secondary:
    """... `collect` turns a per-outcome FitError into a recorded failure instead of a raise.

    Default False, which is the point-estimate contract unchanged: on the workbook a FitError is
    fatal and must be, because there is no replicate to drop. Stage 10 passes True, and then a
    single outcome's unfittable m_a(X) leaves the other six estimated -- which is [§10]'s "dropped
    and counted" at the granularity [§10] left open and Stage 10 §7.3 chose.

    `collect` NEVER catches SchemaError, at either setting. S6 and S7 stay fatal (§5.4).
    """
```

One loop, one home, and the paths are still read rather than derived. §15.7 asserts both settings: that
`collect=False` still raises, and that `collect=True` returns six estimates and one recorded failure on
a frame where exactly one outcome cannot be fitted.

### 7.4 `len(fit.alpha)` is not degenerate, and Stage 8 said what that means

Stage 8 §11 asked for this distribution and stated the stakes: *"if that distribution is degenerate at 6
the gap is theoretical, and if it is not, the drop rate is partly a function of a bound calibrated on
frames unlike the ones being dropped."* Measured over 1998 replicates:

```
  cutpoints   replicates    rate
      6           1907     95.4%     the point estimate's shape
      5             91      4.6%     a declared mRS level carried no positive weight
```

**It is not degenerate.** Stage 8 §5.3 collapses the response to the categories carrying positive
weight, so a replicate missing a declared mRS level fits a coefficient on a coarser scale, and [§10]
takes percentiles across them. Under proportional odds they are the same parameter — which is [§8]'s own
untested assumption and now, per DECISION 5's [§16] amendment, an assumption the report must state.
**4.6% of the primary's draws are on a five-cutpoint scale**, and that is the number that statement is
about.

The specific consequence Stage 8 feared does **not** arise, and saying so is as much the finding as the
distribution is. Its worry was that G7's bound was calibrated at seven occupied categories only — *"if
that distribution is not degenerate at six cutpoints, the drop rate is partly a function of a bound
measured on frames unlike the ones being dropped"* — but the drop rate is **zero** (§7.5), so the bound
drops nothing at any cutpoint count and there is no rate for it to be a function of. What survives is
the comparability question, not the calibration one. `TODOS.md`'s "Re-measure the separation band at
fewer than seven occupied mRS categories" is answered to that extent and no further (§16).

### 7.5 The two counters both read zero — and `Σw` is where the replicates show

```
  primary, over 1998 live replicates
    G7 (separation)                    0     0.0%
    nonconvergence                     0     0.0%
    degenerate_design                  0     0.0%
    SchemaError                        0     0.0%
  propensity.fit
    FitError                           2     0.1%
    SchemaError                        0     0.0%     <- §5.2, was 26.0%
```

**Both of Stage 8's counters read zero, and Stage 8's warning about that is why the pair is reported
anyway.** Its words: *"A Stage 10 failure counter reading zero is therefore not evidence that no
replicate was degenerate unless G7 is in the path."* Here G7 *is* in the path, it is exercised by a
deliberate construction (§15.5), and it reads zero on the data. That is a stronger statement than either
counter alone, and it is the statement only the pair can make.

Stage 8 also predicted that its own expectation about `Σw` was untested and named the symptom that would
show if it were wrong: *"A stratified resample of 92 records has a similar `Σw` by construction, so
nothing is expected to move — but 'expected' is not 'measured' ... If a replicate's `Σw` came in an order
of magnitude low, the same absolute tolerance would be an order of magnitude looser relative to the
objective, and the visible symptom would be `iterations` dropping rather than anything failing."*
Measured:

```
  Sw over in_model, per replicate     min 4.9045   p1 10.9815   median 23.5118
                                      p99 33.0853  max 37.4814
                                      point estimate 27.736623            [Stage 8 §5.2]

  polr iterations                     min 1   median 3   max 4
                                      Stage 8 measured 4 to 5 on every frame it had  [Stage 8 §2]

  n in_model per replicate            88: 1   89: 26   90: 120   91: 374   92: 741   93: 736
```

**`Σw` moved and the predicted symptom is present.** The minimum is 4.90 against the point estimate's
27.74 — a factor of 5.7, not an order of magnitude, but well outside "similar by construction" — and the
iteration count runs down to 1 where Stage 8 never saw below 4. `POLR_TOL` and `POLR_SCORE_TOL` are
absolute (`config.py:376-377`), so at `Σw` = 4.90 the same tolerance is 5.7 times looser relative to the
objective, and a fit that stops after one iteration is a fit that met an absolute criterion on a small
objective rather than one that converged hard.

**What this does and does not establish.** It does not establish that any replicate's `β` is wrong: the
score criterion is on the gradient, every fit reported `converged_on` and nothing failed, and a
one-iteration Newton step from a good start value is a normal thing on a smooth small problem. It
establishes that Stage 8's "nothing is expected to move" was wrong as stated, that the diagnostic pair it
asked for is the one that shows it, and that a relative-tolerance form of `POLR_TOL` is now a question
with evidence behind it rather than a style preference. §16 files it with the trigger; it is **not**
changed here, because changing a convergence tolerance is changing the estimator in every replicate and
that is a [§8] decision, not a [§10] one.

`n in_model` reaching **93** in 736 replicates is the covariate-incomplete row not being drawn at all,
which is the same coin as §5.1's 26.0% and is why the two blockers were always going to be found at the
same time.

---

## 8. The percentile interval [§10]

### 8.1 The rule

[§10]: *"Percentile 95% confidence intervals; seed recorded."* So for each estimand, over that
estimand's surviving draws (§7.3), the limits are the `(1 − level)/2` and `1 − (1 − level)/2` quantiles
— the 2.5th and 97.5th percentiles at `C.CI_LEVEL` = 0.95.

Nothing is centred, reflected, bias-corrected or accelerated. A percentile interval is a pair of order
statistics of the replicate distribution and is **not** a function of the point estimate — §15.8 asserts
that, because "the interval is centred on the estimate" is what a reader assumes and an implementation
that quietly did it would produce plausible numbers. §18 records why no second interval exists.

### 8.2 The percentile definition, and `numpy`'s default breaks [§10]'s own claim

[§10] makes a claim about the relationship between its p-value and its interval:

> *"so it is the smallest level at which the percentile interval for `β` excludes the null, and it
> agrees by construction with the reported interval for the common odds ratio."*

**"By construction" is true of one percentile definition and false of `numpy`'s default.** `np.percentile`
defaults to `method="linear"`, which *interpolates between* order statistics; the p-value is built from
`Pr(β̂* ≤ 0)`, which counts them. At the boundary the two disagree.

Measured: 30 000 constructed draw sets at `B` = 2000, smaller-tail count drawn uniformly on 40..60, no
ties at zero —

```
  method            disagreements with p < 0.05      rate      of which tail count == 50
  linear                                  1412      4.707%                          1412
  higher                                  1412      4.707%                          1412
  nearest                                 1412      4.707%                          1412
  midpoint                                  35      0.117%                            35
  inverted_cdf                               0      0.000%                             0
  lower                                      0      0.000%                             0
```

**Every disagreement is at a smaller-tail count of exactly 50**, and that is not a coincidence to be
smoothed over — it is the decision boundary itself. At `B` = 2000 the p-value is a multiple of
`2/B` = 0.001, so `p < 0.05` iff the smaller tail holds at most 49 draws. At exactly 50, `p = 0.0500`,
which is **not** less than 0.05, so the interval must include the null — and `linear` interpolates the
2.5th percentile between the 50th draw (non-positive) and the 51st (positive), landing on whichever side
the arithmetic falls. The constructed sweep at exact tail counts shows it directly:

```
  tail count   p        p < 0.05    linear   inverted_cdf   lower   higher   nearest   midpoint
      49       0.0490     yes        excl        excl        excl     excl     excl      excl
      50       0.0500      no        EXCL        incl        incl     EXCL     EXCL      incl
      51       0.0510      no        incl        incl        incl     incl     incl      incl
```

**And the method is not enough on its own: the quantile ARGUMENT has to be exact.**
`100.0 * ((1.0 − level) / 2.0)` is `2.500000000000002` at level 0.95, because `1.0 − 0.95` is
`0.050000000000000044`. `inverted_cdf` takes the order statistic at index `ceil(q/100 · n)`, so at
`n` = 2000 that 2e-15 excess moves the index from 50 to 51 — from the largest negative draw to the
smallest positive one. Measured on one tail-50 draw set: **`−0.107095` at the literal 2.5 against
`+0.101248` at the computed quantile**, opposite signs and therefore opposite verdicts, at exactly the
tail count the pin was chosen for. §11's `percentile_ci` snaps both quantiles with `round(..., 9)` and
says so in its docstring; §21b records that executing this document's own fence is what found it.

**So `C.PERCENTILE_METHOD = "inverted_cdf"` is prespecified, and it is not a default being written
down.** `inverted_cdf` is the empirical-CDF quantile — the smallest order statistic whose cumulative
proportion reaches the level — which is the *same object* `Pr(β̂* ≤ 0)` is computed from. That identity
is why the two agree, and it is why the pin is `inverted_cdf` rather than `lower`, which agrees here by
arithmetic coincidence at these levels rather than by being the same estimator of the same quantile.

**On the workbook, all three agree.** Over the cohort's own 2000 draws, `linear`, `inverted_cdf` and
`lower` give the same verdict for `β` and for all seven risk differences (§21). So the workbook does not
witness this and the pin is justified by construction rather than by the data — which is the honest
version, and is why §15.8 tests it on a constructed fixture with the tail count set to 50 rather than on
the cohort.

**Ties at exactly zero make the test conservative and do not break the agreement.** If any draw is
exactly 0.0 then `Pr(≤0) + Pr(≥0) > 1` and `2·min` counts those draws in the smaller tail, inflating p.
Measured across five tie constructions (2 to 100 exact zeros), the interval included the null in every
one and the verdicts agreed under every method. The condition is reachable in principle — a safety
outcome with no weighted events in either arm gives a weighted risk difference of exactly 0.0 — but is
**unreachable on this cohort**, because §5.4's S8 raises on a constant outcome before an estimate is
produced: measured, `rd == 0.0` in **0** of 1998 replicates for all seven outcomes. §16 files the
interaction, because reclassifying S8 to `FitError` keeps it unreachable and a future decision to
estimate a constant outcome instead of dropping it would make it reachable.

### 8.3 The floor below which no interval is emitted, and it is derived

At `level` and `n` draws, `inverted_cdf` puts the lower limit at order statistic `ceil((1−level)/2 · n)`,
one-based. For that index to exceed 1 — for the limit to be an *interior* order statistic rather than the
sample minimum — `n` must exceed `2/(1−level)`, which at 0.95 is **40**.

**At exactly 40 the lower limit is the sample minimum, and the upper limit is NOT the maximum** — an
earlier draft said "the two limits are the sample extremes" and that is wrong on the upper half.
Measured on `np.arange(40.0)`: `inverted_cdf` returns `[0.0, 38.0]`, because `ceil(0.025 · 40) = 1`
gives index 0 while `ceil(0.975 · 40) = 39` gives index 38, the second-largest. The floor is a
statement about the **lower** limit, which is the one that stops carrying information first; below 40
it is the minimum and `np.percentile` returns it while still calling it a percentile.

`C.ci_min_draws` is that expression and **the function is given once, in §13, rather than twice.** An
earlier draft carried it here as well; the two copies were identical, which is the condition under which
they stop being identical later, and this document's own single-source rule (§3.2) applies to its
fences and not only to its counts.

**And it is not `int(np.ceil(2.0 / (1.0 - level)))`, which is what that draft said and is wrong.**
`1.0 - 0.90` is `0.09999999999999998`, so the quotient is `20.000000000000004` and the ceiling is **21**
against an intended 20. At `CI_LEVEL` = 0.95 the error happens to fall the other way — `1.0 - 0.95` is
`0.050000000000000044`, the quotient is `39.99999999999996`, and the ceiling is the intended 40 — so the
bug is invisible at the one level this study uses and appears at the first sensitivity level anybody
tries. §13's form snaps to the nearest integer when the quotient is within `1e-9` of one and takes the
ceiling otherwise; §21b records that executing the fence is what found this.

An estimand whose surviving draws fall below the floor gets **no interval** — `run` records the estimand
with its `Draws` and its counters and no `Interval`, rather than emitting a pair of extremes labelled as
a percentile. `C.N_BOOT` = 2000 against a floor of 40 means this needs 98% of replicates to fail, and
the worst measured drop rate is `sich`'s 0.8%, so **the branch is unreachable on v7** — the same status
Stage 9 §6.4 recorded for its correction, and stated for the same reason: an unreachable branch that is
specified is a branch a future cohort does not discover the hard way. §15.8 reaches it with a fixture.

**The surviving-draw count is reported beside every interval regardless of whether the floor was
approached** (§3.1's `Interval.n_draws`), because [§10]'s "dropped and counted" is only informative
where the count is printed.

### 8.4 Coverage, measured

The roadmap's first acceptance clause: *"On synthetic data with a known effect the interval covers the
truth at roughly the nominal rate."* Two designs, because they establish different things, and both are
specified as code in §15.0.

```
  design U — unconfounded: A depends on Z and X; Y depends on A ALONE
    truth                       0.700000    closed form: the marginal odds ratio IS the coefficient
    n = 400, M_outer = 300, B_inner = 300, method = "inverted_cdf"
    coverage of 95%             0.9600      MC se 0.0113  ->  0.9378 .. 0.9822
    usable intervals            300 / 300   inner fits dropped: 1
    mean interval width         0.7622

  design C — confounded: Y depends on A, Z and X
    truth                       0.659272    n = 200,000 reference fit of the SAME estimator
    n = 400, M_outer = 300, B_inner = 300, method = "inverted_cdf"
    coverage of 95%             0.9400      MC se 0.0137  ->  0.9131 .. 0.9669
    usable intervals            300 / 300   inner fits dropped: 0
    mean interval width         0.7366
```

Nominal 0.95 is inside both Monte Carlo intervals. **The standard error is quoted because a coverage
estimate from 300 outer simulations is itself a binomial proportion**, and 0.9400 quoted bare reads as
undercoverage when it is 0.8 standard errors from nominal.

**Design C's truth is 0.659272 and not 0.700000, and that gap is the design working.** A logistic model
is not collapsible, so the marginal ATO common odds ratio is attenuated relative to the conditional
coefficient the data were generated from. A coverage test that used 0.7 as design C's truth would report
about 88% coverage and read as a bootstrap defect; the truth has to be the estimand, obtained the way
Stage 9 §8.5 obtained its ATO — by running the estimator once at very large `n`.

**What this does not establish.** `B_inner` is 300 and not `C.N_BOOT` = 2000, and `n` is 400 and not 93:
300 × 300 × 2 designs is already 180 000 fits and 540 s. So coverage is measured for the *procedure* at a
size where the Monte Carlo error is tolerable, not for this cohort at this `B`. It does not establish
coverage at n = 93, where Stage 9 §8.5 measured that the per-replicate spread swamps effects of the size
at issue; it does not exercise `resample`'s stratification, because these designs have no strata; and it
does not exercise the [§7] Firth propensity model on a near-separated design, which is the condition
[§13]'s amendment says this cohort is actually in. §16 files the last of those as the one that matters.

---

## 9. p-values [§10]

### 9.1 The primary, and it is one test of one parameter

[§10] specifies it exactly, and the specification includes what it may be called:

```
  p = 2 · min{ Pr(beta* <= 0), Pr(beta* >= 0) },   floored at 1 / (B + 1)
```

- **`β`, not the odds ratio, not a risk difference.** The null is `β = 0`.
- **`B` is the number of SURVIVING draws, not `C.N_BOOT`.** The floor is the smallest p-value a
  bootstrap of that size can express, and expressing it as `1/(N_BOOT + 1)` on a thinned set would claim
  a resolution the draws do not have. On the primary the two are equal on this run (nothing was
  dropped), which is exactly the condition under which a wrong implementation is green — so §15.8
  asserts it on a fixture where they differ.
- **The floor is a maximum, not a clamp on both ends**: `max(2·min(...), 1/(B+1))`. A p of 0 is not
  reportable from 2000 replicates.
- **It is labelled a test of the proportional-odds treatment coefficient and never a
  risk-difference-scale test.** [§10] says it *"must not be described as one"*. The label is Stage 14's
  to print and this stage's to make printable: `Interval.p` on the primary carries the estimand key
  `"beta"`, and there is no field anywhere that pairs a p-value with an `RD_k`.

### 9.2 Secondary and safety, on the risk-difference scale

[§10]'s reasoning is quoted rather than restated, because it is an argument about selection and a
paraphrase loses it:

> *"`β` is a fitted model parameter and is defined whenever the model converges, so no scale problem
> arises. That is not true of the *marginal* odds ratio for a binary outcome, which is a ratio of
> weighted proportions and becomes undefined when one reaches 0 or 1; taking p from those draws would
> condition on the replicates where the ratio happens to exist — selection on the outcome, which badly
> biases rare-event p-values. **For the secondary binary and safety outcomes, therefore, p is computed
> on the risk-difference scale**, which is defined in every replicate. Testing `RD = 0` is the same null
> as `OR = 1`."*

So each binary outcome gets **one p-value, from its `rd` draws**, by §9.1's formula. Specifically:

- **No p-value on `odds_ratio`**, ever, and the reason is [§10]'s above. This matters more than it
  reads: measured, the continuity correction fires in **16.1%** of `sich`'s replicates, **13.4%** of
  `tici_2b_3`'s and **2.5%** of `ph2`'s (§10.2). Those are precisely the draws [§10] is refusing to
  condition on, and they are a sixth of one outcome's.
- **No p-value on `augmented`.** [§10] names the risk-difference scale and the model-assisted estimator
  is on that scale, so a p could be computed — and is not, because [§10] prescribes one test per
  outcome and two would be a multiplicity Stage 11's Benjamini-Hochberg is not told about. The
  augmented estimate gets an **interval**, and the comparison [§10.3] draws is between two intervals on
  one scale. PI decision, §17.
- **The correction's own firing rate travels with the interval**, because Stage 9 §6.3 requires the
  correction *"flagged and printed beside the estimate ... wherever the estimate appears"* and a
  bootstrap interval is somewhere it appears. `Draws` carries no `or_corrected` count field; the rate
  goes in the audit entry (§10.1), and §12.3 records it as Stage 14's to print.

### 9.3 No p-value on the cumulative `RD_k`

[§8], via the roadmap's Stage 10 entry: *"Cumulative `RD_k`: intervals only, **no p-values**."* All six
carry `Interval.p is None`. §15.8 asserts it as a positive check on all six keys rather than as an
absence, because "we did not compute one" and "we computed one and it is `None`" look identical from the
outside and only the first is what [§8] asks for.

The reason is [§16]'s, as amended by DECISION 5: the six threshold differences exist to let a reader see
whether the constant-shift assumption behind one odds ratio is doing work. Six p-values beside them would
turn a diagnostic into six tests nobody prespecified and nobody corrects for.

### 9.4 The agreement is asserted, not assumed

[§10] claims the p-value and the interval agree by construction. §8.2 measured that the claim is
definition-dependent, so the agreement is a test rather than a remark:

- On the cohort's own draws, for `β` and all seven `rd` groups, agreement holds under `inverted_cdf`
  (and, on these draws, under `linear` and `lower` too) — 8 of 8 (§21).
- On a constructed fixture with the smaller-tail count set to exactly **50**, agreement holds under
  `inverted_cdf` and **fails under `linear`** — which is the test that can fail, and the one §15.8
  requires to have been seen failing before the pin is trusted.

---

## 10. The diagnostics, and the questions they close

### 10.1 One audit entry, and no new kind

`data.py`'s `KINDS` stays nine, as it has for three stages. The counters and the distributions go in a
single `model` entry, `bootstrap_replicates`, taking the ledger **30 → 31**:

```
  load 7 / derive 4 / classify 1 / build 6 / fit 4 / assess 2
        / primary 3 / secondary 3 / bootstrap 1   =  31
```

One entry and not three, unlike Stages 8 and 9. The reason is that all of it is one fact — what the
2000 replicates looked like — and Stage 9 §10.1's precedent for splitting was that its three entries
described three different populations. Here there is one population of replicates, described by a
counters table and a distributions table, and `AuditEntry` holds one table: so they are concatenated
into one grid with the second block's header kept as a labelled row and short rows padded, which is
exactly what `propensity._pad` already does for Stage 6's two-table entry (`propensity.py:380-398`).
That function is `propensity.py`'s private and is **not** imported — §0.2 says `propensity.py` is
untouched, and reaching into it for a helper would be a module-boundary violation to save six lines.
`bootstrap.py` builds its own grid.

**`case_ids` is empty on this entry.** Nine `model` entries name cases and this one removes no patient;
naming the 2000 × 93 drawn rows would be a log the size of the data. The excluded-record accounting is
`propensity.fit`'s, per replicate, in a throwaway nobody reads — which is §6.2's point.

### 10.2 `max|beta|`, and Stage 9 attributes its tail to an outcome that fits no model

Stage 9 §9.6 prescribed **no bound** on the nuisance coefficient and asked Stage 10 for a distribution
instead. The roadmap carries the request and Stage 9's reason for it:

> *"Stage 9 prescribes no bound, so nothing turns a separated nuisance fit into a countable failure: the
> `FitError` rate is 0.3% for `sich` and 0% elsewhere while 15.1% of `sich`'s fits exceed `|β| = 14`.
> The counters will read near-zero on outcomes whose nuisance models are degenerate a sixth of the
> time."*

**`sich` is unaugmented. The specified estimator fits it no `m_a(X)` at all**, so it has no nuisance
coefficient, no `max|beta|`, and no fits to be degenerate a sixth of the time. Stage 9 §9.3 establishes
this itself — five augmented (four full, one reduced) and two unaugmented, `sich` and `ph2` — and §9.6's
measurement is of a model the analysis does not fit.

The measurement **reproduces**, which is what makes this a misattribution rather than an error. Fitting
`m_a(X)` for the two unaugmented outcomes anyway, over 400 replicates:

```
  outcome    n     min      median      max        > 8       > 14
  sich      400   0.6848    6.6661    31.5907     41.8%     16.2%    Stage 9: 42.7% / 15.1%
  ph2       400   1.2082    3.7834    64.3557     12.2%      2.5%
```

41.8% against Stage 9's 42.7% and 16.2% against its 15.1%: the same distribution. So Stage 9 measured
something real and correctly; the roadmap sentence built on it describes a counterfactual model.

**The diagnostic is still needed, still for Stage 9's reason, and the tail is on different outcomes.**
Over 1998 replicates, for the five outcomes that actually have an `m_a(X)`:

```
  outcome           n      min     median      p95        max      > 8      > 14
  mrs_0_2_90d    1998   0.8016    3.8252    7.8317    20.8773     4.3%     0.4%
  mrs_0_1_90d    1998   0.6916    2.4818    6.2976    33.8315     1.9%     0.5%
  tici_2b_3      1995   1.1391    2.5965    5.2081     9.6145     0.2%     0.0%
  death_90d      1997   0.8935    4.4245    8.9904    59.5360     8.4%     0.7%
  mrs_5_6_90d    1998   1.4198    4.9896    9.6896    52.6980    11.6%     1.0%

  FitError from m_a(X):  death_90d 1, every other 0.  So 0.05% and not Stage 9's 0.3%.
```

Stage 9's structural argument is untouched and is why no bound is added here: `m_a(X)` is a nuisance
whose *predictions* enter `tau`, predictions are probabilities, so every term is bounded and there is no
`exp(β)`-style tail in the estimate however large the coefficient gets. What is lost is precision, and
the tail is real — `death_90d` reaches 59.54 and `mrs_5_6_90d` 52.70, both safety outcomes, with 8.4%
and 11.6% above 8. So the sentence Stage 9 wanted written is true with two names substituted, and
`TODOS.md`'s *"Establish whether `sich`'s `max|beta|` tail leaves its interval usable"* is asking about a
fit that does not exist. §16 rewrites it against `death_90d` and `mrs_5_6_90d`.

`max_abs_beta` is keyed by outcome and holds entries **only for augmented outcomes** — `sich` and `ph2`
are absent from the dict rather than present with an empty array, because an empty array reads as "we
looked and found none" where absence reads as "there is nothing here to look at", and the second is
true.

### 10.3 The rendering rule the roadmap already stated

The roadmap's Stage 10 entry: *"The name is ASCII deliberately: it becomes a column in a rendered
markdown table, `data._md_table` does no escaping, and a literal pipe in a cell breaks the table in a
way byte-identity checks do not catch."*

`data._md_table` (`data.py:242-245`) computes its column widths from `rows[0]` and builds its separator
from `len(widths)`, so a cell containing `|` yields a body row with more cells than the separator has
dashes — which renders as a broken table while remaining a byte-identical string. So:

- The column is **`max_abs_beta`**, not `max|beta|`, in every rendered cell and header.
- No cell in either block may contain `|`. §15.13 asserts it over the whole grid, as Stage 9 §15.14 does
  for its three tables.
- Every row of the concatenated grid has the same cell count, including the padded short rows and the
  kept second header (§10.1). §15.13 asserts that too, because §10.1's concatenation is the exact
  operation Stage 6 §7.2 got wrong first time.

Prose in this document writes `max|beta|` freely; the constraint is on rendered output.

---

## 11. `run`, written out

```python
def run(df: pd.DataFrame, ps: propensity.Propensity, est: outcome.Primary,
        sec: outcome.Secondary, audit: Audit) -> Bootstrap:
    """The [§10] bootstrap: N_BOOT stratified replicates, percentile limits, prespecified p-values.

    Takes the point estimate's `Primary` and `Secondary` and does not recompute them. It reads
    exactly two things from them -- `sec.estimates[k].augmented_path` for §6.3's frozen map, and the
    estimand keys -- and no estimate, because a percentile limit is an order statistic of the draws
    and is not a function of the point value (§8.1, §15.8).

    Takes no `n`, no `seed` and no `level`: C.N_BOOT, C.SEED and C.CI_LEVEL are the prespecified
    ones and [§10] refits this in every replicate, so they are not runtime knobs. This is Stage 6
    §6.4's choice, made for its reason -- `propensity.fit` takes no covariate list because the
    specification is not an option -- and `replicates` below IS parameterised, because that one is
    general and this one is [§10].

    Five things about the order below, each of which is a failure if moved:

      * `_assert_run_inputs` runs FIRST, so a caller passing a `Secondary` over a different frame
        reports R3/R5 rather than producing 2000 replicates of a misalignment (§4.4).
      * `paths` is built ONCE, before the loop, from `sec`. Building it inside would make it a
        function of the replicate, which is §6.3's whole prohibition.
      * ONE Generator, created once, consumed in order. Not one per replicate (§4.3).
      * The per-replicate `Audit` is constructed INSIDE the loop and discarded at the end of it, so
        the process's memory is not a function of N_BOOT (§6.2).
      * `_record_replicates` runs AFTER the loop and before the intervals are built, so a log exists
        naming the counters even if a percentile call then raises on an estimand nobody expected to
        be empty.

    NO SchemaError IS CAUGHT ANYWHERE IN THIS FUNCTION, and that is §7.1 as code. The only `except`
    is on `model.FitError`, inside `_replicate`, per estimand group.

    Deterministic given C.SEED, and writes no file. Appends ONE `model` entry (§10.1). The only
    mutable object touched is the `Audit` passed in -- the per-replicate ones are its own.
    """
    _assert_run_inputs(df, ps, est, sec)                         # R1-R8
    paths = {key: e.augmented_path for key, e in sec.estimates.items()}
    rng = np.random.default_rng(C.SEED)

    collected = [
        _replicate(resample(df, rng, C.BOOT_STRATUM), paths, audit.source)
        for _ in range(C.N_BOOT)
    ]

    draws = _collect(collected, _estimand_keys(sec))             # §7.3, per group
    diagnostics = _diagnostics(collected)                       # §7.4, §7.5, §10.2
    _record_replicates(draws, diagnostics, audit)               # §10.1

    intervals = {}
    for key, d in draws.items():
        if len(d.draws) < C.ci_min_draws():                     # §8.3
            continue
        lo, hi = percentile_ci(d.draws, C.CI_LEVEL)
        p = bootstrap_p(d.draws) if _tested(key) else None       # §9.1, §9.3
        intervals[key] = Interval(lo, hi, C.CI_LEVEL, C.PERCENTILE_METHOD, len(d.draws), p)

    return Bootstrap(C.SEED, C.N_BOOT, draws, intervals, diagnostics)
```

And the two arithmetic functions, which are the whole of §8 and §9 and are separate from `run` so that
§15.8 can drive them with arrays it wrote:

```python
def percentile_ci(draws: np.ndarray, level: float = C.CI_LEVEL) -> tuple[float, float]:
    """The [§10] percentile limits, at the ONE definition §8.2 pins.

    `method=C.PERCENTILE_METHOD` is passed explicitly and is never left to default: numpy's default
    is "linear", which interpolates between order statistics and disagrees with `bootstrap_p` below
    at exactly the tail count that decides significance -- measured, 4.707% of constructed draw sets
    (§8.2). The two functions are one decision and this argument is where it is recorded.

    Raises on fewer than `ci_min_draws(level)` draws rather than returning an extreme under a
    percentile's name (§8.3). `run` checks the count before calling, so the raise is a caller bug.

    **THE QUANTILES ARE SNAPPED AND `100.0 * ((1.0 - level) / 2.0)` IS WRONG.** That expression
    returns 2.500000000000002 at level 0.95, not 2.5, because `1.0 - 0.95` is
    0.050000000000000044. The excess is 2e-15 and it is not cosmetic: `inverted_cdf` takes the order
    statistic at index ceil(q/100 * n), so at n = 2000 the index moves from 50 to 51 and the limit
    moves from the largest negative draw to the smallest positive one. Measured on one draw set with
    a tail count of 50: `-0.107095` at the literal 2.5 against `+0.101248` at the computed quantile
    -- opposite signs, so opposite verdicts, at exactly the tail count §8.2 chose this method for.
    An earlier form of this fence had it, and §21b is where running the fence found it.
    """
    if len(draws) < ci_min_draws(level):
        raise C.SchemaError(
            f"percentile_ci: {len(draws)} draw(s) against a floor of {ci_min_draws(level)} at level "
            f"{level:g}. Below the floor the lower limit is the sample minimum and np.percentile "
            "returns it while still calling it a percentile [Stage 10 §8.3].")
    q_lo = round(50.0 * (1.0 - level), 9)          # 2.5 EXACTLY at level 0.95 -- see the docstring
    q_hi = round(100.0 - q_lo, 9)                  # 97.5 exactly
    lo, hi = np.percentile(draws, [q_lo, q_hi], method=C.PERCENTILE_METHOD)
    return float(lo), float(hi)


def bootstrap_p(draws: np.ndarray) -> float:
    """[§10]'s two-sided bootstrap p, floored at 1/(B+1) where B is the SURVIVING draw count.

    `B` is `len(draws)` and never C.N_BOOT: the floor is the smallest p a bootstrap of this size can
    express, and using N_BOOT on a thinned set claims a resolution the draws do not have. On this
    cohort nothing was dropped from the primary, so the two are equal -- which is exactly the
    condition under which an implementation reading N_BOOT is green, and why §15.8 asserts it on a
    fixture where they differ.

    Draws exactly equal to 0.0 are counted in BOTH tails, so `Pr(<=0) + Pr(>=0) > 1` and the test is
    conservative. That is deliberate and measured (§8.2): the alternative -- splitting ties -- makes
    p a function of a tie-breaking rule nobody prespecified.
    """
    p_le = float(np.mean(draws <= 0.0))
    p_ge = float(np.mean(draws >= 0.0))
    return max(2.0 * min(p_le, p_ge), 1.0 / (len(draws) + 1))
```

---

## 12. Data flow into Stages 11-14

```
  run(cohort, ps, est, sec, audit) returns Bootstrap
   ├─ seed, n_boot          [§16] requires both reported                       (§3.1)
   ├─ draws                 estimand key -> Draws
   │    ├─ draws            the surviving replicate values                     (§7.3)
   │    ├─ n_attempted      what was asked for                                 (§7.3)
   │    └─ failures         bucket -> count; sums to n_attempted               (§7.2)
   ├─ intervals             estimand key -> Interval(lo, hi, level, method,
   │                                                 n_draws, p)               (§8, §9)
   └─ diagnostics           n_alpha, polr_iterations, sum_w, max_abs_beta,
                            n_in_model                                         (§7.4, §7.5, §10.2)

  the cohort frame, the Propensity, the Primary and the Secondary come back unchanged  (§0.2)

  audit: load 7 / derive 4 / classify 1 / build 6 / fit 4 / assess 2
         / primary 3 / secondary 3 / bootstrap 1  =  31                        (§10.1)

  every interval limit and every p-value on the workbook is DELIBERATELY ABSENT
  from this document. They are in the gitignored log.                          (§4.5)
```

### 12.1 Stage 11 [§13]

- **It reads `intervals[key].p` and does not recompute one.** [§13]'s Benjamini-Hochberg runs *within*
  the secondary and safety families, and Stage 9 §3.1 already gives it `Secondary.by_family()` so it
  does not partition the outcomes itself. What it owes: raw and adjusted p both reported, and the
  primary uncorrected.
- **Seven p-values, not fourteen.** §9.2: one per binary outcome, on the risk-difference scale. There is
  no p on any `odds_ratio` and none on any `augmented`, so the family sizes are three and four and
  Stage 11 must not discover a different count by counting intervals.
- **The E-value reads the primary interval's limit nearest the null**, and the roadmap's rule applies to
  a number this stage produces: *"the limit's E-value is 1.0 whenever the interval already spans the
  null"*. `Interval.lo` and `.hi` are what that rule is evaluated on, and `Interval.p` agrees with them
  by §9.4 — so Stage 11 can test null-spanning either way and get the same answer, which is the property
  §8.2 exists to guarantee.
- **The full-covariate sensitivity arm [§13, DECISION 4] needs its own bootstrap and this stage does not
  provide it.** Stage 11's arm refits [§7] over `PS_COVARIATES_FULL`, so its interval is a second `run`
  over a second `Propensity` — which needs the covariate seam Stage 11 §11 already identifies as its one
  design decision. Nothing here blocks it and nothing here builds it.

### 12.2 Stages 12 and 13 [§14a, §14b]

They reuse **`resample`, `replicates` and `percentile_ci`** and none of `run`. [§14] fits no propensity
model anywhere, so there is no `e`, no `w`, no `h` and no augmentation; the replicate body is the
standardisation, and `replicates` takes it as a callable (§3.2).

Two things they inherit and one they must decide for themselves. Inherited: the stratification is
patient-level within centre, over their own population; and the percentile definition is
`C.PERCENTILE_METHOD`, because §8.2's argument is about arithmetic and not about this cohort. Theirs to
decide: **their strata are not these strata.** [§14a] runs over all eligible patients at all centres,
which includes USZ — a centre with no treated patients — so a stratum can exist there that contributes
no contrast, and whether that is resampled or held fixed is [§14a]'s question and not this document's.
§17 records that this stage does not decide it.

### 12.3 Stage 14 [§16]

Five things it owes that this stage cannot enforce, each because a prespecified rule says so:

- **Every interval prints its `n_draws` and its failure buckets.** [§10]'s "dropped and counted" is a
  reporting requirement, not just a bookkeeping one, and `Draws.failures` is keyed so that the
  separation count appears beside the convergence count per Stage 8 §11.
- **The primary p-value is labelled a test of the proportional-odds treatment coefficient**, never a
  risk-difference-scale test [§10]. And the six `RD_k` carry intervals and no p-values [§8]; `p is None`
  on those six is what makes printing one a `None` reaching a formatter (§9.3).
- **The constant-shift statement is computed, not written** [§16 amendment, DECISION 5] — and §7.4 gives
  it a second clause this stage measured: **4.6% of the primary's draws fit five cutpoints rather than
  six**, so "the same parameter across replicates" is [§8]'s proportional-odds assumption doing work in
  the interval as well as in the point estimate.
- **Stage 9 §6.3's correction is named wherever the estimate appears**, and a bootstrap interval is
  somewhere it appears: the firing rates are 16.1% (`sich`), 13.4% (`tici_2b_3`) and 2.5% (`ph2`) of
  replicates, against **0.0%** on the point estimate. So an interval on a marginal odds ratio for those
  three is an interval a sixth of whose draws are corrected while the point estimate is not — which is
  Stage 9 §6.4's *"materially a choice about two intervals and barely one about any point estimate"*
  arriving where it lands.
- **"Model-assisted", never "doubly robust"**, and the reduced specification printed beside the TICI
  estimate — both inherited from Stage 9 §12.2 unchanged.

---

## 13. What Stage 10 amends in Stages 1-9

The full ledger, so that no amendment is discovered during implementation. **Two shipped modules and
three test modules, and both shipped changes are in `outcome.py`.**

| File | Amendment | Why |
|---|---|---|
| `config.py` | `CI_LEVEL`, `BOOT_STRATUM`, `PERCENTILE_METHOD`, `FAILURE_BUCKETS`, `ci_min_draws()` — in the inference block beside the existing `SEED` and `N_BOOT` (`config.py:269-270`) | §8, §9, §7.2. Prespecified for `FIRTH_*`'s reason: [§10] refits in every replicate, so a percentile definition that changes an answer is a property of the sampling distribution and not a runtime knob. `ci_min_draws` is a function because it is **derived from** `CI_LEVEL` and a second constant could disagree with it (§8.3) |
| `outcome.py` | **S8 raises `FitError` and not `SchemaError`** (`outcome.py:924-930`) | §5.4. One line plus the docstring paragraph that already anticipates it (`outcome.py:897`). S6 and S7 unchanged. Without it 1.0% of replicates raise something Stage 8 §11 prohibits catching |
| `outcome.py` | **`secondary` gains `collect: bool = False`** | §7.3. The per-outcome replicate sets the drop-granularity decision requires, expressed where Stage 9's seven-outcome loop already lives rather than re-implemented in `bootstrap.py` — which Stage 9 §14 names as how §6.3 gets violated by accident. Default `False` is the point-estimate contract unchanged. **No other function changes** |
| `tests/test_outcome.py` | S8's reclassification, and `collect` at both settings | §15.7. Stage 9's S8 assertions change error class; its other sections are untouched |
| `tests/test_config.py` | **six** assertions | `0.0 < CI_LEVEL < 1.0`; `BOOT_STRATUM` is a column the [§6] balance set knows about; `PERCENTILE_METHOD` is one numpy accepts, asserted by calling `np.percentile` with it rather than against a string list that can drift from numpy; `ci_min_draws(CI_LEVEL) == 40` and `N_BOOT >= ci_min_draws(CI_LEVEL)`, which is R8 as a static check; and `FAILURE_BUCKETS`' values are exactly the four §7.2 names, because a fifth bucket appearing by typo would silently collect the counts Stage 8 §11 asked to be separated |
| `tests/fixtures_stage10.py` | **new file, and this row is the one place its contents are enumerated.** Six functions — `two_centre_frame`, `separable_ordinal_frame`, `constant_outcome_frame`, `unfittable_nuisance_frame`, `boundary_draws`, `known_effect_population` — and five constants: `COVERAGE_SEED`, `M_OUTER`, `B_INNER`, `BETA_TRUE`, `ALPHA_TRUE` | §15.0. Every constructed fixture in this document, as code, for Stage 9 §22.3 item 1's reason: three of its fixtures were pinned to ten decimals and described rather than given, which made its acceptance criteria unperformable |
| `implementation_roadmap.md` | Stage 10 gains its `**Spec:**` line, five corrections and two additions | §22. Lands with this document |
| `TODOS.md` | two items close, three are answered and rewritten, two are new | §16 |

**Five things that look like they need amending and do not.**

- **`propensity.py`.** §5.5, and it is the central claim of §5 rather than an incidental one. Not
  amended, not imported for its privates (§10.1), and its stale `data.py:255` citation deliberately left
  (§16).
- **`data.py`.** §10.1: no new kind, no new heading, so no literal pin in any existing test module moves
  except the ledger count. `KINDS` stays nine and §15.13 asserts it. Fourth stage running.
- **`model.py`.** Untouched. `FitError` gains no `code` field, and §7.2 argues the deferral rather than
  assuming it; `polr`, `firth`, `design` and `predict` are called and not edited. **`POLR_MAX_ABS_BETA`
  is read and not rewritten** — §7.4 establishes that the cutpoint-count question Stage 8 raised about it
  does not bite, and §7.5 establishes that a different constant, `POLR_TOL`, is the one with new evidence
  against it; neither is changed by a [§10] document.
- **`balance.py`.** Not read and not imported (§0.2). No interval is put on an SMD.
- **`RARE_MINORITY_THRESHOLD` and `OUTCOME_MODEL_OVERRIDES`.** Read by Stage 9 and **not read here at
  all**, which is stronger than unchanged: §6.3's whole rule is that the path is read off the point
  estimate, so `bootstrap.py` naming either constant would be the violation. §15.11 asserts by scan that
  it names neither.

```python
# config.py — the five additions, as they are to be pasted. Stage 10 §13
#
# --- [§10] inference: the bootstrap's definitions ----------------------------------------------
#
# The [§10] stratification variable. A NAME and not the column itself, because `bootstrap.resample`
# is general in its stratum (Stage 10 §12.2): [§14a] and [§14b] resample the same way over a
# population that includes a centre with no treated patients, and a function that hard-coded
# "center" would still be right there while being right for the wrong reason.
BOOT_STRATUM: Final[str] = "center"

# The [§10] confidence level. 0.95 is [§10]'s own "percentile 95% confidence intervals"; it is here
# rather than as a literal so that `ci_min_draws` below can be derived from it and cannot disagree.
CI_LEVEL: Final[float] = 0.95

# The percentile DEFINITION, and this constant changes an ANSWER rather than a tolerance
# [Stage 10 §8.2]. numpy's default is "linear", which interpolates between order statistics, and
# [§10] claims its p-value "is the smallest level at which the percentile interval for beta excludes
# the null" and "agrees by construction with the reported interval". THAT CLAIM IS FALSE UNDER THE
# DEFAULT: measured over 30000 constructed draw sets at B = 2000, "linear" disagrees with the
# p-value in 4.707% of them, and EVERY disagreement is at a smaller-tail count of exactly 50, where
# p is exactly 0.0500 and an interpolated limit lands on whichever side of zero the arithmetic falls.
#
# "inverted_cdf" is the empirical-CDF quantile -- the smallest order statistic whose cumulative
# proportion reaches the level -- which is the SAME object Pr(beta* <= 0) is computed from. That
# identity is why the two agree, and it is why this is not "lower", which agrees at these levels by
# arithmetic coincidence rather than by estimating the same quantile.
#
# Measured: on this workbook's own 2000 draws all three methods give the same verdict for beta and
# for all seven risk differences, so the cohort does not witness this and the pin is justified by
# construction. Prespecified for FIRTH_*'s reason: [§10] refits in every one of N_BOOT replicates.
PERCENTILE_METHOD: Final[str] = "inverted_cdf"

# The FitError buckets Stage 8 §11 requires reported separately, keyed by the leading token of the
# raised message [Stage 10 §7.2]. `model.FitError` carries no code -- it is a bare RuntimeError
# subclass (model.py:89) -- and all SIXTEEN raise sites -- fourteen in model.py, two in outcome.py --
# themselves by that token, so classification is textual and test_bootstrap.py §15.6 SCANS both
# modules and asserts every token found is a key here. A reworded message is then a test failure
# rather than a counter that silently reads zero.
#
# "polr:" covers two failures -- step-halving exhausted and no convergence in POLR_MAX_ITER -- and
# they share a bucket. Stage 8 §11 asks for G7 against everything else, and Stage 8 §5.2 records
# that exhausting the halvings is not a convergence route; both measured 0 in N_BOOT replicates.
FAILURE_BUCKETS: Final[dict[str, str]] = {
    "G7": "separation",
    "S8": "constant_outcome",
    "G6": "degenerate_design",
    "polr:": "nonconvergence",
    "Firth:": "nonconvergence",
    "O1": "degenerate_design", "O2": "degenerate_design", "O3": "degenerate_design",
    "O4": "degenerate_design", "O5": "degenerate_design", "O6": "degenerate_design",
    "F3": "degenerate_design", "F4": "degenerate_design", "F6": "degenerate_design",
}


def ci_min_draws(level: float = CI_LEVEL) -> int:
    """The fewest draws at which a `level` percentile limit is an order statistic at all.

    DERIVED, not chosen. Under PERCENTILE_METHOD the lower limit is order statistic
    ceil((1-level)/2 * n), one-based; for that index to exceed 1 -- for the limit to be interior
    rather than the sample minimum -- n must exceed 2/(1-level), which is 40 at CI_LEVEL = 0.95.
    At exactly the floor the two limits ARE the sample extremes, which is the weakest interval the
    definition can produce; below it np.percentile returns an extreme while still calling it a
    percentile [Stage 10 §8.3].

    A function rather than a constant so it cannot disagree with CI_LEVEL. N_BOOT = 2000 against a
    floor of 40 means an estimand needs 98% of its replicates to fail before it loses its interval:
    measured, the worst per-outcome drop rate on this cohort is 0.8%, so the branch is unreachable
    on v7 and is specified anyway.

    THE SNAP IS NOT DEFENSIVE PROGRAMMING AND `int(np.ceil(2.0 / (1.0 - level)))` IS WRONG.
    `1.0 - 0.90` is 0.09999999999999998, so that expression returns 21 where 20 is intended. At
    CI_LEVEL = 0.95 it happens to return the intended 40 -- `1.0 - 0.95` errs the other way and the
    quotient is 39.99999999999996 -- so the bug is invisible at the only level this study uses and
    appears at the first [§13] sensitivity level anybody tries. A level whose quotient is genuinely
    non-integral still takes the ceiling: at 0.93 the quotient is 28.571 and the answer is 29.
    """
    exact = 2.0 / (1.0 - level)
    nearest = round(exact)
    return int(nearest if abs(exact - nearest) < 1e-9 else np.ceil(exact))
```

---

## 14. Handover to Stage 11

```
  cohort = cohort.build(eligibility.classify(derive.derive(*data.load()), ...), ...)
  ps     = propensity.fit(cohort, audit)
  est    = outcome.primary(cohort, ps, audit)
  sec    = outcome.secondary(cohort, ps, audit)
  boot   = bootstrap.run(cohort, ps, est, sec, audit)

    2000 / 20260807  N_BOOT / SEED, both reported [§16]                      [§4.3]
    0                SchemaError from propensity.fit, was 26.0%              [§5.2]
    19 / 1.0%        replicates S8 raised on before reclassification         [§5.4]
    0 / 0            G7 / nonconvergence on the primary, and the PAIR is
                     the finding rather than either number                   [§7.5]
    1907 / 91        replicates fitting 6 / 5 cutpoints -- NOT degenerate    [§7.4]
    4.90 .. 37.48    Sw per replicate against the point estimate's 27.74,
                     with polr iterations down to 1 from Stage 8's 4-to-5    [§7.5]
    40.6% / 2.6%     ph2 / sich replicates whose own cell would decide the
                     augmentation path differently -- hence frozen           [§6.3]
    59.54 / 52.70    max|beta| of m_a(X) on death_90d / mrs_5_6_90d, the two
                     outcomes the tail is actually on                        [§10.2]
    31               audit entries, taking the ledger 30 -> 31               [§10.1]
    ~145 s           serial wall-clock with §6.4's shared design             [§2]

    every interval limit and every p-value on the workbook is DELIBERATELY
    ABSENT from this ledger and from this document. They are in the
    gitignored log.                                                          [§4.5]
```

**Four things Stage 11 must build that this stage cannot.**

1. **The Benjamini-Hochberg correction over seven p-values, not fourteen.** §9.2: one p per binary
   outcome, on the risk-difference scale, families of three and four from
   `Secondary.by_family()`. There is no p on any `odds_ratio` and none on any `augmented`, so counting
   `intervals` to find the family sizes gives the wrong answer.
2. **The E-value, from `Interval.lo`/`.hi` on the primary.** The roadmap's rule — the limit's E-value is
   1.0 whenever the interval spans the null — is evaluated on this stage's limits, and §9.4 guarantees
   that testing null-spanning via `p` and via the limits agree.
3. **A second `run` for the [§13, DECISION 4] full-covariate arm**, over a `Propensity` refit on
   `PS_COVARIATES_FULL`. This stage provides the engine and not the seam; Stage 11's roadmap entry
   already names that seam as its one design decision.
4. **The subgroup interaction tests**, which are [§13]'s and are hypothesis-generating. Nothing here
   resamples within a subgroup, and a subgroup interval taken by re-slicing these draws would be
   conditional on the full-cohort resample rather than stratified within the subgroup.

**What Stage 11 must not do**, each because a prespecified rule says so: describe the primary p as a
risk-difference-scale test [§10]; attach a p-value to any of the six `RD_k` [§8]; correct the primary for
multiplicity [§13]; take the E-value from a null-crossing limit without the 1.0 rule (roadmap); or
recompute any p-value or limit, since `Interval` carries the definition that produced it and a second
computation under a different one is §8.2's disagreement re-introduced downstream.

---

## 15. Acceptance criteria

Tests live in `test_bootstrap.py`, with section banners matching these numbers, as `test_cohort.py`,
`test_model.py`, `test_propensity.py`, `test_balance.py` and Stages 8's and 9's sections do. The S8 and
`collect` sections — §15.7 — live in `test_outcome.py`, because that is where the changed function is.
Tests needing the workbook are gated as Stage 7 built the gate and are tagged `[data-gated]`.

### 15.0 The frames and fixtures this stage is tested on

Every one is given as code, with its seeds and constants written out, and every pin in this document is
measured against these and only these. Six functions and five constants (§13).

```python
COVERAGE_SEED: Final[int] = 20260807
M_OUTER: Final[int] = 300
B_INNER: Final[int] = 300
BETA_TRUE: Final[float] = 0.7
ALPHA_TRUE: Final[tuple[float, ...]] = (-2.0, -1.0, 0.0, 1.0, 2.0)   # 6 categories, 5 cutpoints


def two_centre_frame() -> pd.DataFrame:
    """A cohort-shaped frame with two strata and EXACTLY ONE covariate-incomplete row.

    The witness for §5.1 and §5.2, and the shape is what makes it one: the blocker needs a stratum
    small enough that a duplicate draw is likely and exactly one excluded row, so that two excluded
    rows collapse to one name. Stratum sizes 8 and 6 give P(a given row drawn >= 2x) of 0.37 and
    0.42, so the naive resampler raises within a handful of replicates rather than needing hundreds.

    Carries C.PS_COVARIATES, C.TREATMENT, `case_id`, `center` and the primary outcome. No value is
    taken from the workbook.
    """


def separable_ordinal_frame(n: int = 20) -> tuple[pd.DataFrame, np.ndarray, np.ndarray]:
    """Treatment determines mRS entirely: `a = 1` -> Y = 0, `a = 0` -> Y = 6, unit weights.

    Drives G7, and Stage 8 §6.1's measurement is reproduced on it (§21): converges on the SCORE
    criterion in 17 iterations with rescales and halvings both 0, every fitted quantity finite,
    beta = 36.4058, exp(beta) = 6.469e15, alpha = -18.2029. That is the whole point of the fixture:
    separation does not present as failure, so G7 is the only thing that turns it into one, and a
    "failures" counter without G7 in it reads zero here (§15.5).
    """


def constant_outcome_frame(key: str) -> pd.DataFrame:
    """A frame on which `key` takes one value across its whole [§11] population.

    Drives S8, at both its classifications: §15.7 asserts `FitError` after §5.4's change and that
    `C.SchemaError` is NOT raised, which is the assertion that fails on the landed code and is
    therefore the one that has to exist.
    """


def unfittable_nuisance_frame() -> pd.DataFrame:
    """A frame on which EXACTLY ONE augmented outcome's `m_a(X)` cannot be fitted and six can.

    Drives `collect` (§7.3, §15.7). "Exactly one" is the specification: the test asserts six
    estimates and one recorded failure, and a frame on which two fail would pass a wrong
    implementation that gives up after the first.
    """


def boundary_draws(tail: int, n: int = C.N_BOOT, rng=None) -> np.ndarray:
    """`tail` draws strictly below zero and `n - tail` strictly above. NO draw equal to zero.

    The §8.2 fixture, and `tail` is a parameter because the finding is entirely about ONE value of
    it: at n = 2000 the p-value is a multiple of 2/n, so p < 0.05 iff tail <= 49, and tail = 50 is
    where p is exactly 0.0500 and `method="linear"` interpolates the limit onto the wrong side. The
    test sweeps 49, 50, 51 and asserts the verdict at each -- 50 alone distinguishes the methods, and
    a test written at 45 or 55 passes under every one of them (§15.8).
    """


def known_effect_population(n: int, rng: np.random.Generator, confounded: bool) -> pd.DataFrame:
    """§8.4's coverage generator. `A` depends on Z and X; `Y` depends on A alone, or on all three.

    `confounded=False` is design U: the marginal odds ratio IS BETA_TRUE in closed form, so coverage
    measures the machinery. `confounded=True` is design C, where the conditional coefficient is NOT
    the marginal ATO odds ratio -- logistic models are not collapsible -- and the truth must be
    obtained by running the estimator once at n = 200,000, as Stage 9 §8.5 obtained its ATO.
    Measured: 0.659272 against BETA_TRUE = 0.7, and a coverage test using 0.7 there reports about
    88% and reads as a bootstrap defect (§8.4).
    """
```

### 15.1 The preconditions, one frame per branch

```
  branch  phase  the frame                                       asserts
  ------  -----  ---------------------------------------------   ------------------------------
  R1        1    drop the BOOT_STRATUM column                    SchemaError naming R1
  R2        1    drop the `case_id` column                       SchemaError naming R2
  R3        1    a Propensity whose e/w index is PERMUTED but
                 equal as a set                                  SchemaError naming R3
  R4        1    ps.in_model with dtype int8                     SchemaError naming R4
  R5        1    a Secondary missing one BINARY_OUTCOMES key      SchemaError naming R5
  R6        2    a frame with a NaN stratum label                SchemaError naming R6
  R7        2    a Secondary with augmented_path = "partial"      SchemaError naming R7
  R8        2    level such that ci_min_draws exceeds N_BOOT      SchemaError naming R8
```

**R3's frame is permuted-but-equal and not merely different**, for Stage 9 §22.3 item 19's reason: Stage
9 §12 measured that an order mismatch pairs each record's weight with another record's fitted value and
returns a different number with no raise, and `index.equals` is order-sensitive where a set comparison is
not. This stage constructs a `Propensity` per replicate, so it is the caller Stage 9's S3a was written
against.

**Phase 1 failures are collected**, so a frame missing two columns reports both. §15.1 asserts a
two-failure frame produces one `SchemaError` naming both codes, because an implementation raising on the
first satisfies every single-branch test above.

### 15.2 `resample`

- **The stratum totals are exact.** For every replicate, `draw[stratum].value_counts()` equals the input
  frame's, on a frame with unequal strata. **And the companion is what makes it a test**: a resampler
  drawing `len(df)` rows from the whole frame ignoring strata passes a total-row-count check and fails
  this one.
- **Every drawn row's `case_id` is unique**, and the count equals `len(df)`. Measured on the cohort: 93
  rows, 93 distinct, `is_unique` True.
- **The first appearance is `id#1` and not the bare id.** Asserted positively, because the natural
  implementation — rename only the duplicates — is the one that leaves "was this row drawn once" a
  question about string formatting (§5.3).
- **`propensity.fit` does not raise on 400 consecutive replicates** of `two_centre_frame()`
  `[data-gated: no]`, and **the companion is the whole of §5**: the same 400 draws through a resampler
  without the rename raise `C.SchemaError` on a materially non-zero fraction, asserted as `> 0` rather
  than pinned, because the fixture's rate is a property of its stratum sizes and pinning it makes the
  test brittle for nothing.
- **Determinism**: the same seed gives a frame equal under `DataFrame.equals`; a different seed does not.
- **The frame's columns and dtypes are unchanged** from the input, asserted column by column. A
  resampler that silently promotes an `Int64` to `float64` changes `model.design` and hence every
  estimate, which is `TODOS.md`'s Stage 1 dtype item arriving by another route.

### 15.3 The replicate body, and what it may not read

- **`_replicate` refits `propensity.fit` and does not reuse the point estimate's `Propensity`.**
  Asserted by giving it a frame whose `e` would differ and checking that it does — a body that
  accidentally closed over the outer `ps` returns the point estimate's weights and every interval
  collapses toward zero width, which is a failure that produces plausible numbers.
- **The per-replicate `Audit` is discarded**: after `run`, the `Audit` passed in holds exactly one new
  entry (§10.1), not `9 × N_BOOT + 1`.
- **The per-replicate `Audit` carries the caller's `Source`** and not a synthetic one, asserted against
  `audit.source` identity (§6.2).

### 15.4 The failure taxonomy reconciles

- **`len(draws) + sum(failures.values()) == n_attempted`**, for every estimand key. This is the one
  property no wrong implementation satisfies by accident: a silent drop anywhere makes the three numbers
  stop reconciling. Asserted over a synthetic replicate sequence with a known number of injected
  failures, so the expected counts are known rather than read back from the thing under test.
- **`assert not issubclass(model.FitError, C.SchemaError)`**, and the converse. The whole taxonomy is one
  `except` clause away from collapsing.
- **A `SchemaError` raised inside a replicate propagates out of `run`.** Asserted with a body that raises
  one, because §7.1's rule is a negative — "is never caught" — and the only way to test a negative is to
  raise the thing and require it to arrive.

### 15.5 G7 is exercised, and the pair is asserted `[roadmap]`

The roadmap's own clause: *"A test drives a replicate loop over a deliberately separable frame and
asserts the G7 count is non-zero while the convergence count stays zero — the pair, because Stage 8
measured that the second never fires on this estimator and a single 'failures' counter therefore reads
zero on data that is degenerate throughout."*

- Over `separable_ordinal_frame()`, `failures["separation"] > 0` **and**
  `failures["nonconvergence"] == 0`.
- And the fixture's own witness first: `model.polr` on that frame **converges** — 17 iterations on the
  score criterion, every safeguard counter zero, every fitted value finite — so the assertion is about a
  fit that succeeded and was rejected, not about a fit that failed. Measured (§21), reproducing Stage 8
  §6.1 exactly.
- **And the counters read 0 and 0 on the cohort** `[data-gated]`, which is the measurement that makes the
  pair worth reporting rather than the one that makes it worth asserting.

### 15.6 The bucket map is scanned, not trusted

- Scan `model.py` and `outcome.py` for `raise FitError(` and `raise model.FitError(`, extract the first
  whitespace-delimited token of each message, and assert every one is a key of `C.FAILURE_BUCKETS`.
  Measured: **sixteen** sites -- fourteen in `model.py`, two in `outcome.py` -- with tokens `G6`,
  `G7`, `O1`-`O6`, `F3`, `F4`, `F6`, `polr:`, `Firth:`, plus `S8` once §5.4 lands. **G6 was missing
  from an earlier draft of `FAILURE_BUCKETS` and this scan is what found it** (§7.2).
- **The scan is by content and not by line number**, as Stage 7 §18d requires of every scan in this
  repository: a raise site that moves must not fail this test and a raise site whose *message* changes
  must.
- `_bucket` on an unrecognised token **raises rather than defaulting to a catch-all bucket**. A default
  would make §15.6 cosmetic: the scan would pass, the map would be incomplete, and the counter Stage 8
  §11 asked to be separate would be silently merged.

### 15.7 S8 and `collect` `[in test_outcome.py]`

- **S8 raises `model.FitError` and does not raise `C.SchemaError`** on `constant_outcome_frame(key)`,
  for a key in each family. The second half is the assertion that fails on the landed code.
- **S6 and S7 still raise `C.SchemaError`**, asserted on their own frames, because §5.4 changes one
  precondition and a change that moved all three would pass a test written only for S8.
- **`collect=False` raises** on `unfittable_nuisance_frame()` — the point-estimate contract unchanged.
- **`collect=True` returns six estimates and one recorded failure** on the same frame. Both halves:
  that the six are present *and* that the seventh is recorded rather than absent, because an
  implementation dropping the key entirely gives six estimates too.
- **`collect` never catches `C.SchemaError` at either setting**, asserted with a frame that trips S6.
- **The frozen `paths` still governs under `collect=True`**: an outcome whose path is `"unaugmented"`
  fits no `m_a(X)` and so cannot contribute a nuisance `FitError` at all.

### 15.8 The interval, the p-value, and the definition `[roadmap, amended]`

- **The limits are order statistics of the draws and are not a function of the point estimate.** Asserted
  by computing an interval from a draws array, then again with the point estimate moved far away, and
  requiring both limits identical. "The interval is centred on the estimate" is what a reader assumes.
- **The percentile definition is pinned and the pin is exercised at the boundary.** On
  `boundary_draws(50)`: `percentile_ci` under `C.PERCENTILE_METHOD` gives an interval **containing** zero
  while `bootstrap_p` gives exactly 0.0500, so the two agree; **and the companion is the test** — the
  same draws under `method="linear"` give an interval **excluding** zero, which is the disagreement §8.2
  measured, asserted rather than described. `boundary_draws(49)` and `(51)` agree under both, which is
  why a test written at either would pass on the wrong pin.
- **`bootstrap_p`'s floor uses the surviving count.** Asserted on a draws array of length 500 inside a
  `Draws` whose `n_attempted` is 2000: the floor must be `1/501` and not `1/2001`. On the cohort the two
  are equal, which is exactly the condition under which an implementation reading `C.N_BOOT` is green.
- **Ties at exactly zero are counted in both tails.** On draws holding `t` exact zeros, `Pr(≤0) + Pr(≥0)`
  exceeds 1 and p is conservative; asserted for `t` in a few values, with the interval containing zero
  in each.
- **The floor refuses rather than returning a degenerate limit.** `percentile_ci` on 39 draws raises;
  on exactly 40 it returns the **sample minimum** as its lower limit and the second-largest draw as its
  upper, and **both halves are asserted positively** — on `np.arange(40.0)` the answer is `(0.0, 38.0)`
  and not `(0.0, 39.0)`. A test asserting the maximum passes on a wrong percentile definition and fails
  on the pinned one (§8.3).
- **The six `RD_k` carry `p is None`**, asserted on all six keys positively (§9.3). And no key anywhere
  in `intervals` pairs a p-value with an `odds_ratio` or an `augmented` estimand (§9.2).

### 15.9 The shared design, and the companion that makes it a test

- **The four full-list outcomes' designs compare equal** under `DataFrame.equals` on a replicate, and
  their `in_estimate` masks are identical. Measured over 300 replicates of the cohort `[data-gated]`: 0
  differ.
- **And the companion is what makes it a test**: on a frame where one of the four has an extra missing
  value — violating roadmap invariant 6 rather than relying on it — the shared design is **shown to be
  wrong**, and `_shared_design` is required to detect the mask mismatch and fall back to per-outcome
  designs rather than silently reusing one. An optimisation whose precondition is never violated in a
  test is an optimisation nobody has tested.
- **`tici_2b_3` builds its own design**, asserted, because its covariate list is the [§8] override and
  sharing the full-list design would silently estimate the wrong nuisance model.

### 15.10 No fallback estimator, ever `[roadmap invariant 5]`

- **A replicate whose augmented fit fails contributes nothing to that outcome's `augmented` draws and is
  not substituted with its unaugmented `rd`.** Asserted by injecting a nuisance failure and checking the
  `augmented` draw count fell by exactly one while the `rd` count did not — which also asserts §7.3's
  granularity in the one place it is observable.
- **A replicate whose `polr` hits G7 contributes nothing to `beta` and nothing to the six `RD_k`.** The
  primary is one group (§7.3), so all seven counts fall together, and that is asserted as a group rather
  than key by key.
- **No code path in `bootstrap.py` calls an estimator that the point estimate did not use.** Scanned:
  the module names `propensity.fit`, `outcome.primary` and Stage 9's public estimators, and no other
  fitter.

### 15.11 The module boundary

- **`bootstrap.py` does not import `balance`** (§0.2).
- **`bootstrap.py` names none of `_augmentable`, `_augmented_path`, `RARE_MINORITY_THRESHOLD`,
  `OUTCOME_MODEL_OVERRIDES`.** By AST scan, and this is §6.3's rule as a test: the path is read off the
  point estimate and the only way to guarantee it is not re-derived is to assert the vocabulary of
  re-deriving it is absent.
- **`bootstrap.py` reaches into no other module's privates**, and `propensity._pad` is the named case
  (§10.1) — it does exactly what this stage's grid needs and is not imported.
- **`resample`, `replicates`, `percentile_ci` and `bootstrap_p` name no covariate, no outcome and no
  centre.** By scan. That is §0.1's generality claim made checkable, and it is what Stages 12 and 13
  depend on.

### 15.12 Stage 10 adds nothing, refits nothing, and breaks nothing earlier

- The Stages 1-9 suite passes **unedited** except for §15.7's S8 sections, and the diff to
  `test_outcome.py` touches only those.
- `run` returns the cohort frame, the `Propensity`, the `Primary` and the `Secondary` unmutated —
  asserted by comparing before and after, because none of them is frozen against in-place mutation of
  its `pd.Series` fields.
- The audit ledger is **31**, and `data.KINDS` is still nine (§10.1).

### 15.13 The audit entry and its rendering

- **One entry, `("model", "bootstrap_replicates")`**, with `case_ids` empty (§10.1).
- **Every row of the concatenated grid has the same cell count**, including padded short rows and the
  kept second header. This is the operation Stage 6 §7.2 got wrong first time and the reason `_pad`
  exists there.
- **No cell in the grid contains `|`**, asserted over every cell, and the `max_abs_beta` column is
  spelled in ASCII (§10.3).
- **The log is byte-identical across two hash seeds.** Spelled out as a heredoc, as Stages 7, 8 and 9
  spell theirs out:
  ```bash
  cd extended_bridging
  S='import data, derive, eligibility, cohort, propensity, outcome, bootstrap
  df, a = data.load()
  df = cohort.build(eligibility.classify(derive.derive(df, a), a), a)
  ps = propensity.fit(df, a)
  bootstrap.run(df, ps, outcome.primary(df, ps, a), outcome.secondary(df, ps, a), a)
  print(a.to_markdown())'
  PYTHONHASHSEED=0 uv run python -c "$S" > /tmp/s10-a.md
  PYTHONHASHSEED=1 uv run python -c "$S" > /tmp/s10-b.md
  diff /tmp/s10-a.md /tmp/s10-b.md && echo IDENTICAL
  ```

### 15.14 Coverage `[roadmap]`

- **Design U covers at roughly nominal**, with the Monte Carlo standard error asserted rather than the
  point coverage: the test requires nominal 0.95 inside `coverage ± 1.96 · se`, which at
  `M_OUTER` = 300 is a band of about ±0.022. Measured 0.9600, se 0.0113.
- **Design C covers at roughly nominal against the estimand and not against `BETA_TRUE`.** Measured
  0.9400, se 0.0137, truth 0.659272. **And the companion is the point**: the same run scored against
  `BETA_TRUE` = 0.7 gives about 88% coverage, and the test asserts *that* too — so a future edit that
  "fixes" design C's truth to the generating coefficient fails loudly instead of reading as a bootstrap
  defect (§8.4).
- Both are `M_OUTER` × `B_INNER` = 90 000 fits and take about 270 s each, so both are marked slow and
  gated behind the same mechanism the R oracles use. §16 files the reduction from `C.N_BOOT`.

### Coverage map

```
  resample                    15.2, 15.3
  replicates                  15.2 (determinism), 15.4 (the injected-failure sequence)
  percentile_ci               15.8 (definition, boundary, floor, order-statistic property)
  bootstrap_p                 15.8 (formula, floor-on-surviving-count, ties)
  run                         15.3, 15.4, 15.10, 15.12, 15.13; 15.14 for the whole loop
  _assert_run_inputs          15.1, one frame per branch R1-R8, plus the two-failure frame
  _replicate                  15.3, 15.5, 15.10
  _bucket                     15.6, including the unrecognised-token raise
  _collect                    15.4 (reconciliation), 15.10 (granularity)
  _estimand_keys              15.8 (the RD_k keys carry no p), 15.12 (the ledger)
  _shared_design              15.9, both the equality and the violated-precondition companion
  _counters_table             15.13
  _diagnostics_table          15.13
  _replicates_detail          15.13
  _record_replicates          15.13, 15.12

  config.ci_min_draws         15.8 (39 raises, 40 returns the extremes), test_config.py
  outcome._assert_estimable   15.7 (S8 reclassified; S6, S7 unchanged)
  outcome.secondary           15.7 (collect at both settings)

  THE MAP IS THE COUNT, and Stage 9's rule applies unchanged: there is no total, because a total not
  derived from the rows above cannot be checked against them. Every route listed has a test that can
  fail; a route not listed is not covered.

  UNREACHABLE ON v7, SPECIFIED ANYWAY, AND REACHED BY A FIXTURE INSTEAD:
    the ci_min_draws floor          needs 98% of replicates to fail; worst measured 0.8%   §8.3
    an rd of exactly 0.0            S8 raises first; measured 0 of 1998                    §8.2
    the percentile disagreement     all three methods agree on the cohort's own draws      §8.2
    G7 firing                       0 of 1998                                              §7.5
    nonconvergence                  0 of 1998                                              §7.5
```

### Failure modes

| Failure | Test | Error handling | Visible or silent |
|---|---|---|---|
| `np.percentile` left at its default `method` | 15.8 | none | **Silent, and it changes a significance verdict in 4.7% of boundary cases** (§8.2). This is the failure the pin exists to make catchable |
| `bootstrap_p` floors at `1/(N_BOOT+1)` instead of `1/(B+1)` | 15.8 | none | **Silent on this cohort**, where nothing was dropped from the primary and the two are equal |
| The augmentation path re-derived per replicate | 15.11 | none | **Silent, and it mixes two estimators in one interval for `ph2` in 40.6% of replicates** (§6.3) |
| A `SchemaError` caught inside the loop | 15.4 | raises | Visible — but only because §5.4 removed the 1.0% of replicates that made catching it tempting |
| The point estimate's `Propensity` closed over instead of refit | 15.3 | none | **Silent, and every interval collapses toward zero width** |
| A single "failures" counter instead of the pair | 15.5, 15.6 | none | **Silent: it reads zero on a frame that is separated throughout** (Stage 8 §6.1) |
| `_bucket` defaults an unrecognised token to a catch-all | 15.6 | raises | Visible, and the scan is what keeps it so |
| The shared design reused when the four masks differ | 15.9 | detects and falls back | Would be **silent**: four outcomes estimated on one wrong design |
| Design C's coverage scored against `BETA_TRUE` | 15.14 | none | **Silent, and it reads as a bootstrap defect at about 88%** (§8.4) |
| An `Interval` emitted below the draw floor | 15.8 | raises | Visible |
| A pipe in a rendered cell | 15.13 | none | **Silent: byte-identical and broken** (§10.3) |

---

## 16. Known gaps carried forward

1. **`POLR_TOL` and `POLR_SCORE_TOL` are absolute, and `Σw` moves by a factor of 5.7 across
   replicates.** §7.5. Stage 8 predicted the symptom and this stage measured it: iterations run down to
   1 where Stage 8 never saw below 4. Not changed here, because a convergence tolerance is part of the
   estimator and [§10] does not get to alter [§8]'s. **Trigger:** a replicate that fails to converge, or
   a Stage 12 design where `Σw` is larger still. → `TODOS.md`.
2. **`death_90d` and `mrs_5_6_90d` carry a `max|beta|` tail with no bound and no countable failure.**
   §10.2. 8.4% and 11.6% of their nuisance fits exceed 8, reaching 59.54 and 52.70, while the `FitError`
   rate is 0.05%. Stage 9 §9.6's structural argument — predictions are bounded, so `tau` has no
   `exp(β)`-style tail — still holds and is why nothing is added. **Trigger:** the replicate-level
   correlation between an outcome's `tau` draws and its `max|beta|`, computable from draws already
   taken. → `TODOS.md`, rewritten from `sich` (§22).
3. **`FitError` carries no code and classification is textual.** §7.2. §15.6's scan makes a drifted
   message fail loudly, which is a mitigation and not a fix. **Trigger:** the first bucket that is
   wrong, or the first stage that needs to branch on a failure kind rather than count it. → `TODOS.md`.
4. **Coverage is measured at n = 400 and B = 300, not at n = 93 and `N_BOOT`.** §8.4. Neither design has
   strata, so `resample`'s stratification is unexercised by the coverage test, and neither is
   near-separated, which is the condition [§13]'s amendment says this cohort is actually in. **Trigger:**
   a coverage design carrying centre-like near-determination — which is the one that would say something
   about this study rather than about the procedure. → `TODOS.md`.
5. **4.6% of the primary's draws are on a five-cutpoint scale.** §7.4. Under proportional odds they are
   the same parameter; [§8] does not test that and [§15] declines to. Reported rather than resolved, and
   [§16]'s constant-shift statement is where a reader meets it. **Trigger:** none — this is a property of
   the estimator, filed so it is not rediscovered.
6. **The draw depends on the stratum label alphabet.** §4.3. Renaming a centre gives different, equally
   valid replicates. **Trigger:** a `CENTER_RECODE` edit, at which point the run summary's seed no longer
   identifies the draw by itself. → `TODOS.md`.
7. **`propensity.py:413` cites `data.py:255` for a normalisation now at `data.py:271`.** §5.5. Left
   deliberately: it is in the one module §5 claims is untouched. **Trigger:** the next commit that edits
   `propensity.py` for any reason. → `TODOS.md`.
8. **No parallelism.** §18. At ~145 s serial the gap is not worth the determinism risk. **Trigger:**
   Stage 12's bootstrap, which resamples a larger population, or the [§13] sensitivity suite, which is
   five more `run` calls.

**Two items close and three are rewritten** (§22). `TODOS.md` gains items 1, 3, 4, 6, 7 and 8 above,
loses "Make the pipeline resamplable" and "Decide whether a constant outcome on a replicate is
`SchemaError` or `FitError`", and rewrites "Establish whether `sich`'s `max|beta|` tail leaves its
interval usable" against the two outcomes that have one and "Re-measure the separation band at fewer
than seven occupied mRS categories" against §7.4's measurement.

---

## 17. What Stage 10 deliberately does not decide

- **How a replicate's identifiers are formed.** §5.2 renames the k-th draw `case_id#k`. Neither [§10]
  nor any amendment specifies it. **PI decision, 2026-08-25. PI-reversible**: the alternative is to
  amend `propensity._record_exclusion`, and the consequence of reversing is that the guard weakens on
  the point estimate too, where §5.2 argues it is doing real work.
- **Whether a failure drops a replicate for one estimand or for all.** §7.3 keeps per-outcome sets.
  [§10] is silent and Stage 9 §14 constrains only that `tau` and `rd` share a set. **DECISION 7 (PI,
  2026-08-25). PI-reversible**: reversing costs 20 replicates in 1998 — 1.0% — across all seven binary
  outcomes rather than across the one that failed, and costs 8.4 ms per replicate more (§2). Reversing
  it would also remove the need for `secondary`'s `collect` parameter entirely, leaving S8 as the only
  change to `outcome.py`.
- **Whether the augmented estimate gets its own p-value.** §9.2 gives it an interval and no p. [§10]
  names the risk-difference scale and one test per outcome; two would be a multiplicity [§13] is not
  told about. **PI decision, PI-reversible**, and the pilot's implementation did give it one (§19).
- **Whether S8 should exist at all rather than be reclassified.** §5.4 makes it a `FitError`. Two
  alternatives were put and both declined: estimating the constant outcome, because Stage 9 §6.2
  measured that the correction then returns a finite odds ratio near 1 for an outcome with no contrast
  in it; and rejecting-and-redrawing such replicates, which conditions the resample on the outcome and
  biases every estimate rather than only the affected one. **DECISION 6 (PI, 2026-08-25).
  PI-reversible**, and §8.2 records that reversing it makes exact-zero draws reachable and the tie
  behaviour live.
- **How Stages 12 and 13 stratify.** §12.2. [§14a]'s population includes a centre with no treated
  patients, and whether such a stratum is resampled is that stage's question.
- **Nothing about the balance diagnostics.** No interval on an SMD, and DECISION 4's exceedance is
  reported from Stage 7's number by Stage 14 (§0.2).

---

## 18. NOT in scope for Stage 10

| Considered | Why deferred |
|---|---|
| BCa and bootstrap-t intervals | [§10] prescribes percentile intervals and one interval. BCa needs an acceleration constant from a jackknife over 93 patients and a bias correction from the same draws; bootstrap-t needs a standard error per replicate, which Stage 8 §5.6 establishes cannot be had from the fitter's information matrix — *"no standard error and no interval leaves this stage"*. Both would be second estimators of the same quantity, and [§10] prescribes one |
| An analytic or sandwich standard error | Stage 8 §5.6, and it is a correctness claim rather than a preference: the weights are a tilting function of an *estimated* propensity score, so the inverse observed information omits the variability of estimating `e`, which [§7] states enters the interval |
| A centre-level cluster bootstrap | [§10] declines it explicitly: *"a handful of centres cannot support cluster-bootstrap consistency"*, and inference is conditional on these centres |
| Parallelism across replicates | §16 item 8. ~145 s serial. Every parallel implementation makes the stream a function of the scheduling unless each replicate is separately seeded, which §4.3 declines for its own reason. Deferred, not refused |
| Bootstrapping the balance diagnostics | §0.2. [§9] is a statement about the realised sample |
| An interval on the E-value | Stage 11's, and [§13] prescribes the E-value at the point estimate and at the limit nearest the null — both functions of numbers that already have intervals |
| Caching the draws to disk | The whole `Bootstrap` is about 60 000 floats. Persisting it would create a second source for numbers the log already carries, and `out/` is where the log goes |

---

## 19. What already exists, and what to lift

| From `pilots/` | Status |
|---|---|
| `analysis.bootstrap` (`pilots/analysis.py:535-630`) — the overall shape: stratified draw, refit everything, percentile limits, count failures | **Structure lifted, implementation not.** It is the right shape and it is the direct ancestor of [§10]'s text; six specific things in it are what §5 through §9 exist to do differently |
| `except Exception: failures["other"] += 1` (`analysis.py:595-596`) | **Not lifted, and it is the sharpest reason this stage has a taxonomy.** It counts every exception as a dropped replicate — a `KeyError` from a typo, a `SchemaError` from a malformed frame, an unfittable model — so the failure count is the union of "sparse data" and "bug". Stage 8 §11's rule is the correction, and §7.1 is it as code |
| Two failure buckets, `ps_fit` and `other` (`analysis.py:560-561`) | **Not lifted.** This is exactly the single-counter failure Stage 8 §11 warns about: no G7 bucket exists, so a separated ordinal fit that converges is not counted anywhere. §7.2 |
| `if m.sum() < 10 or len(np.unique(d[m])) < 2: continue` (`analysis.py:573-574`) | **Not lifted, twice over.** It re-evaluates a rare-outcome rule **per replicate**, which Stage 9 §9.5 and §6.3 forbid — `ph2` would cross in 40.6% of replicates. And it `continue`s, so the outcome is dropped from that replicate **without being counted**, which is [§10]'s "dropped and counted" with the second half missing |
| `p = 2 * min(frac, 1 - frac)` where `frac = mean(draws <= null)` (`analysis.py:608-609`) | **Not lifted.** [§10] specifies `2·min{Pr(β̂* ≤ 0), Pr(β̂* ≥ 0)}`, and `1 − Pr(≤0)` is `Pr(>0)`, not `Pr(≥0)`. Identical for draws with no ties at zero and different when there are: the pilot's form is anti-conservative where [§10]'s is conservative (§8.2, §9.1) |
| `np.percentile(draws, [2.5, 97.5])` (`analysis.py:615`) | **Not lifted.** Default `method="linear"` — the §8.2 defect, in the ancestor implementation |
| `if len(draws) < 50: return np.nan` (`analysis.py:605-606`) | **Not lifted.** The right instinct with an undeclared constant: 50 is a magic number, and §8.3's floor of 40 is derived from `CI_LEVEL`. And it returns `nan` where §8.3 refuses, so a `nan` limit flows onward as a number |
| `n_boot=max(200, n_boot // 4)` for subgroups and sensitivity (`run_all.py:140, 188, 193`) | **Not lifted.** A reduced replicate count with no stated reason and no record in the output that it was reduced. §8.4 reduces `B_INNER` for coverage and says so, in the document and in the fixture |
| `point = one(df)` recomputed inside `bootstrap` (`analysis.py:587`) | **Not lifted.** `run` takes `Primary` and `Secondary` in (§11), so the point estimate the intervals accompany is the one Stages 8 and 9 produced and audited, not a second computation that could differ |
| `_strata_indices` | **Not lifted**, but it is the same idea as `resample`'s per-stratum draw. The difference is the rename (§5.3) and the sorted stratum order (§4.3) |

The roadmap's standing note applies: *"Some code has already been built inside pilots. This should not be
a ground source, but can be used for inspiration. Some implementations may be wrong."* Nine rows above
say which.

### 19b. Validation against a reference implementation

**No R oracle is added, and the reason is that there is nothing here for one to check.** Stages 6, 7, 8
and 9 each added one (`ato_psweight.R`, `balance_psweight.R`, `polr_clm.R`, `aug_psweight.R`) because
each implemented an *estimator* with a published counterpart. This stage implements a resampling loop
around estimators that are already oracle-checked, plus two arithmetic functions.

R's `boot::boot.ci(type = "perc")` would be checkable, and it is declined for a specific reason rather
than for cost: `boot.ci`'s percentile type uses `quantile(..., type = 7)`, which is the **linear
interpolation** §8.2 measures as disagreeing with [§10]'s p-value. An oracle agreeing with it would
assert the defect. The empirical-CDF quantile is R's `quantile(..., type = 1)`, and the comparison worth
making is against that.

**What is checked instead, and it is stronger for this stage.** `percentile_ci` under
`inverted_cdf` is asserted against the order statistic computed by hand — `sorted(draws)[ceil(q·n) − 1]`
— on fixtures where the index is unambiguous (§15.8). A closed-form oracle beats a cross-language one
where a closed form exists.

**What 19b does not establish.** That the *loop* is right. Nothing external validates that 2000
replicates of this pipeline have the sampling distribution they should; §8.4's coverage measurement is
the only evidence for that, and §8.4 states its own limits.

---

## 20. Implementation tasks

- [ ] **T0 (P1)** `config.py`: §13's five additions with their comment blocks; `test_config.py`'s six
      assertions. **Gates everything** — no module compiles against `CI_LEVEL` or `FAILURE_BUCKETS`
      until it lands.
- [ ] **T1 (P1)** `outcome.py`: S8 → `FitError` (§5.4), plus `test_outcome.py`'s reclassified sections.
      One line and a docstring paragraph. **Gates T5**: without it 1.0% of replicates raise something
      the loop may not catch, so the loop cannot be run end to end.
- [ ] **T2 (P1)** `tests/fixtures_stage10.py`: §15.0's six functions and five constants. **Gates T6,
      T7, T8, T9** — every pin in this document is measured against these, so nothing that asserts a
      number can be written before it.
- [ ] **T3 (P1)** `bootstrap.py`: `resample`, and §15.2's sections. The blocker's fix is here and it is
      the first thing that has to work.
- [ ] **T4 (P2)** `bootstrap.py`: `percentile_ci`, `bootstrap_p`, `ci_min_draws`'s caller; §15.8. Pure
      arithmetic over arrays, so it needs only T0 and T2.
- [ ] **T5 (P2)** `outcome.py`: `secondary`'s `collect` parameter; §15.7. Depends on T1.
- [ ] **T6 (P2)** `bootstrap.py`: `_bucket`, `_collect`, the `Draws`/`Interval` dataclasses; §15.4,
      §15.6. The reconciliation property is the one to write first.
- [ ] **T7 (P2)** `bootstrap.py`: `_replicate`, `_shared_design`, `replicates`; §15.3, §15.9, §15.10.
      Depends on T3 and T5.
- [ ] **T8 (P3)** `bootstrap.py`: `run`, `_assert_run_inputs`, the diagnostics; §15.1, §15.5, §15.11,
      §15.12.
- [ ] **T9 (P3)** `bootstrap.py`: the audit entry and its grid; §15.13, including the byte-identity
      heredoc.
- [ ] **T10 (P3)** §15.14's coverage tests, marked slow.
- [ ] **T11 (P3)** `implementation_roadmap.md` and `TODOS.md` per §22 and §16. **Lands with the spec
      commit, not with the implementation** — this row is here so the ledger is complete, and it is the
      one task already done when the implementation starts.

`Lanes:` `T0 → T1` first and alone. Then `{T2 → T4}`, `{T3}`, `{T5}` in parallel. Then `T6`, then
`{T7 → T8 → T9}`, then `T10`.

### Definition of done

Complete when all of the following hold, and not before.

1. `uv run pytest -v` is green, and the Stages 1-9 sections are **unedited** except `test_outcome.py`'s
   S8 and `collect` sections (§15.12).
2. **`propensity.py` has a zero-line diff.** §5.5 is the central claim of this stage and this is it as a
   check.
3. **The full loop runs to completion on the workbook with zero `SchemaError`.** Measured before this
   document: 26.0% from `_record_exclusion` and 1.0% from S8. Both must read zero.
4. **The percentile pin has been seen to fail.** `boundary_draws(50)` under `method="linear"` gives an
   interval excluding zero against a p of exactly 0.0500, and the test asserting agreement fails when
   the pin is changed to `"linear"` (§15.8).
5. **The G7 pair has been seen to fail** as a single counter: replacing the two buckets with one
   "failures" count makes §15.5 pass on `separable_ordinal_frame()`, which is the failure Stage 8 §6.1
   describes, and §15.5 as written must fail on that mutation.
6. **The frozen-path scan has been seen to fail**: adding a read of `C.RARE_MINORITY_THRESHOLD` to
   `bootstrap.py` fails §15.11.
7. **The reconciliation property has been seen to fail**: dropping a replicate without incrementing a
   bucket fails §15.4.
8. **The shared-design companion has been seen to fail**: reusing one design when the four masks differ
   fails §15.9.
9. The audit log is **byte-identical across two hash seeds**, by §15.13's heredoc.
10. The audit ledger reads **31** and `data.KINDS` is nine.
11. `_bucket` raises on an unrecognised token, and §15.6's scan finds **sixteen** raise sites and no
    token outside `C.FAILURE_BUCKETS`.
12. **Every count stated in this document has been re-derived before the commit** — the five public
    names, the ten privates, the six fixture functions, the five constants, the sixteen raise sites,
    the thirty-one audit entries.
13. **§21 has no row reading `stated` that a run could have produced.**
14. `TODOS.md`'s two closed items are closed with their measured rates, and the three rewritten items
    name the outcomes and numbers §16 gives.
15. **No interval limit, p-value, `β` or `RD_k` from the workbook appears anywhere under `specs/`**, by
    grep over this document at commit time (§4.5).

---

## 21. Verification record

**Normative.** A number in prose cites the row here that produced it. Where the two disagree **this
table is right and the prose is stale**, and an implementer who finds a conflict should fix the prose
and not re-derive the number.

| Claim | Where used | Verified |
|---|---|---|
| cohort 93 rows; strata HUG 41, Lugano 31, CHUV 21; USZ absent | §4.2 | yes, run |
| exactly one covariate-incomplete record, and it is in Lugano | §5.1 | yes, run |
| analytic P(a given row drawn ≥2× from 31) = 0.2642 | §4.2, §5.1 | yes, run |
| naive resampler: 104 of 400 raise `SchemaError`, 26.0%; message split 75 / 22 / 7 | §5.1 | yes, run |
| renaming resampler: 0 of 400, and 0 of `N_BOOT` = 2000 | §5.2, §7.5 | yes, run |
| a replicate holds 93 rows, 93 distinct `case_id`, `is_unique` True | §5.2, §15.2 | yes, run |
| same seed → identical frame; different seed → different | §4.3, §15.2 | yes, run |
| a duplicated index passes `propensity.fit`, `primary` and `secondary` unchanged (59 distinct labels over 93 rows) | §3.3, §5.3 | yes, run — **this refuted an earlier draft's claim** |
| S8 raises on 16 `sich` and 3 `tici_2b_3` replicates of 1998; 19 replicates, 1.0% | §5.4, §7.3 | yes, run |
| S6 and S7 fire 0 times in 2000 replicates | §3.3, §5.4 | yes, run |
| `secondary` whole loses all seven on 19 S8 + 1 Firth `FitError` = 20 of 1998 | §7.3 | yes, run |
| `outcome.py:897` already names S8's classification as open | §5.4 | yes, read (`outcome.py:897`) |
| G7 = 0, nonconvergence = 0, degenerate_design = 0, `SchemaError` = 0 on the primary over 1998 | §7.5 | yes, run |
| `propensity.fit` `FitError` 2 of 2000, 0.1% | §7.5 | yes, run |
| `len(fit.alpha)`: 1907 at six cutpoints, 91 at five — 95.4% / 4.6% | §7.4, §12.3 | yes, run |
| `polr` iterations min 1, median 3, max 4 | §7.5 | yes, run |
| `Σw` per replicate: min 4.9045, p1 10.9815, median 23.5118, p99 33.0853, max 37.4814 | §7.5 | yes, run |
| the point estimate's `Σw` is 27.736623 | §7.5 | yes, read (Stage 8 §5.2) + run |
| `n in_model` across replicates: 88:1, 89:26, 90:120, 91:374, 92:741, 93:736 | §7.5 | yes, run |
| separated 20v20: 17 iterations on the score criterion, rescales 0, halvings 0, β = 36.4058, exp(β) = 6.469e15, α = −18.2029, all finite | §15.0, §15.5 | yes, run — reproduces Stage 8 §6.1 |
| G7 fires on it, message prefix `'G7  '` | §7.2, §15.5 | yes, run |
| `model.FitError` has no attributes beyond `RuntimeError`; **sixteen** raise sites (14 `model.py`, 2 `outcome.py`); tokens G6, G7, O1-O6, F3, F4, F6, `polr:`, `Firth:` | §3.3, §7.2, §13, §15.6 | yes, read (`model.py:89`) + run — **the count and G6 both corrected by the scan** |
| percentile sweep, 30 000 sets at B = 2000: linear 4.707%, higher 4.707%, nearest 4.707%, midpoint 0.117%, inverted_cdf 0.000%, lower 0.000% | §3.3, §8.2 | yes, run |
| every disagreement is at a smaller-tail count of exactly 50 | §8.2 | yes, run |
| tail counts 49 / 50 / 51 give p = 0.0490 / 0.0500 / 0.0510, and only 50 splits the methods | §8.2, §15.0 | yes, run |
| ties at exactly 0.0: `Pr(≤0)+Pr(≥0) > 1`, interval contains 0 in all five constructions | §8.2, §9.1 | yes, run |
| on the cohort's own 2000 draws, linear / inverted_cdf / lower agree for β and all seven `rd` — 8 of 8 | §8.2, §9.4 | yes, run |
| `rd == 0.0` in 0 of 1998 for all seven outcomes | §8.2 | yes, run |
| `nan` odds ratios: 0 of 1998, all seven | §8.2 | yes, run |
| `ci_min_draws(0.95)` = 40 | §8.3, §13 | yes, run |
| coverage design U: truth 0.700000, coverage 0.9600, se 0.0113, 300/300 usable, 1 inner fit dropped, mean width 0.7622 | §8.4, §15.14 | yes, run |
| coverage design C: truth 0.659272, coverage 0.9400, se 0.0137, 300/300 usable, 0 dropped, mean width 0.7366 | §8.4, §15.14 | yes, run |
| design C scored against `BETA_TRUE` = 0.7 gives about 88% | §8.4, §15.14 | yes, run |
| `or_corrected` fires: `sich` 319/1982 = 16.1%, `tici_2b_3` 267/1995 = 13.4%, `ph2` 49/1998 = 2.5%, others 0 | §9.2, §12.3 | yes, run |
| the correction fires on 0.0% of the point estimate's seven outcomes | §12.3 | yes, run (Stage 9 §6.4 agrees) |
| frozen paths: 4 full, 1 reduced (`tici_2b_3`), 2 unaugmented (`sich`, `ph2`) | §6.3, §10.2 | yes, run |
| path crossing: `ph2` 812/1998 = 40.6%, `sich` 52/1982 = 2.6%, others 0.0% | §6.3 | yes, run |
| Stage 9 measured 46.1% and 3.8% for the same two | §6.3 | yes, read (Stage 9 §9.5) |
| `max|beta|` on the five augmented outcomes: maxima 20.8773, 33.8315, 9.6145, 59.5360, 52.6980; >8 rates 4.3 / 1.9 / 0.2 / 8.4 / 11.6% | §10.2 | yes, run |
| `m_a(X)` `FitError`: `death_90d` 1, others 0 — 0.05% | §10.2 | yes, run |
| `sich` and `ph2` are unaugmented, so the estimator fits them no `m_a(X)` | §10.2 | yes, run + read (Stage 9 §9.3) |
| fitting one anyway: `sich` 41.8% > 8 and 16.2% > 14 over 400; `ph2` max 64.3557 | §10.2 | yes, run — reproduces Stage 9's 42.7% / 15.1% |
| `propensity.fit` 13.4 ms, `primary` 4.7 ms, `secondary` whole 69.9 ms, per outcome 61.5 ms | §2, §7.3 | yes, run |
| the four full-list designs: 18.12 ms built four times, 4.50 ms once — 27.2 s over `N_BOOT` | §2, §6.4 | yes, run |
| all four full-list outcomes have `source = "mrs_90d"` | §6.4 | yes, read (`config.py:479-497`) |
| the four masks differ in 0 of 300 replicates; the four designs compare equal | §6.4, §15.9 | yes, run |
| `tici_2b_3`'s list is `("center", "atrial_fib")`; the full list is nine covariates | §6.4 | yes, run |
| whole loop, both routes, 2000 replicates: 236.8 s, 118.4 ms/rep | §21b | yes, run |
| `propensity.py:413` cites `data.py:255`; the normalisation is at `data.py:271`; `data.py` unchanged since `5fe1750` | §5.5, §16 | yes, read + run (`git log`) |
| `data._md_table` sizes widths from `rows[0]` and does no escaping | §10.3 | yes, read (`data.py:242-245`) |
| the pilot catches bare `Exception`, uses two buckets, `continue`s on a per-replicate rare-outcome rule, uses `1 − Pr(≤0)`, defaults `np.percentile`, and floors at a magic 50 | §19 | yes, read (`pilots/analysis.py:535-630`) |
| `C.SEED` = 20260807, `C.N_BOOT` = 2000 | §4.3 | yes, read (`config.py:269-270`) |

### 21b. This document's own code, executed — 2026-08-25

Every `python` fence in this document was extracted **by content, not by index** — Stage 7 §18d's rule —
and processed in four passes.

```
  pass 1 — extract and parse                                    9 / 9 fences parse
    fence 0  Draws, Interval, Diagnostics, Bootstrap                              §3.1
    fence 1  resample, replicates, percentile_ci, bootstrap_p, run  (stubs)       §3.2
    fence 2  resample                                                             §5.3
    fence 3  paths                                                                §6.3
    fence 4  secondary                                             (signature)    §7.3
    fence 5  run                                                                  §11
    fence 6  percentile_ci, bootstrap_p                                           §11
    fence 7  ci_min_draws, BOOT_STRATUM, CI_LEVEL, PERCENTILE_METHOD,
             FAILURE_BUCKETS                                                      §13
    fence 8  the six fixture functions and five constants                         §15.0

  pass 2 — execute                                               8 / 9 execute
    fence 3 is `paths = {k: e.augmented_path for k, e in sec.estimates.items()}`, a one-line
    excerpt that needs a `Secondary`. Not self-contained by design, and the only one.

  pass 3 — the identities this document asserts                 25 / 25 pass
    ci_min_draws at 0.95, 0.90 and 0.93, plus the naive form's 21 at 0.90 as its companion
    the tail 49 / 50 / 51 sweep: p exact, interval verdict, and p-vs-interval agreement
    tail 50 under "linear" EXCLUDES 0 and under the pin INCLUDES it -- the §8.2 companion
    the floor: raises at 39; at 40 the lower limit is the minimum and the upper is NOT the max
    bootstrap_p's floor on the surviving count; ties at 0.0 counted in both tails
    the limits are order statistics: lo == sorted[ceil(0.025n) - 1]
    the naive quantile 100*((1-level)/2) != 2.5 while the snapped one == 2.5
    resample: stratum totals, distinct case_ids, row count, id#1, determinism, dtypes, index

  pass 4 — the counts, re-derived                                6 / 6 reconcile
    public names 5   privates 10   fixture functions 6   fixture constants 5
    FitError raise sites 16, with no token outside FAILURE_BUCKETS
    audit ledger 7+4+1+6+4+2+3+3+1 = 31
```

**Three defects in this document were found by running it, and all three were in its own code.**
`ci_min_draws` returned 21 instead of 20 at level 0.90 (§8.3); `percentile_ci` computed its quantile as
`100.0 * ((1 − level) / 2.0)` = 2.500000000000002, which flips the limit's sign at exactly the tail
count §8.2 pins the method for (§8.2, §11); and `FAILURE_BUCKETS` was missing `G6`, with the raise-site
count stated as fourteen against a measured sixteen (§7.2). The third is the one §15.6's scan is
specified to catch, caught by the scan, before the scan existed as a test.



**What this does not establish.** That the specification is right. A fence that parses and whose asserted
identities hold is a fence that is internally consistent; §15's criteria are what make it checkable
against an implementation, and the implementation does not exist yet.

---

## 22. What this spec changed elsewhere

| # | Change | Files | The question it settles |
|---|---|---|---|
| 1 | The Stage 10 "Build" bullet gains the renaming rule and the measured post-fix rate | `implementation_roadmap.md` | How the resampler stops `propensity.fit` raising, which the entry named as a blocker and did not resolve |
| 2 | **The four-item Stage 9 handover block gains a fifth**: S8 is a second `SchemaError` on the path, at 1.0% of replicates, and it is not a resampler problem | `implementation_roadmap.md` | The entry says the first blocker is *"the first thing Stage 10 will hit"*. It is, and it is not the only one |
| 3 | **The `max|beta|` bullet is rewritten against `death_90d` and `mrs_5_6_90d`.** As it read, it attributed the tail to `sich` — *"the counters will read near-zero on outcomes whose nuisance models are degenerate a sixth of the time"* — and `sich` is unaugmented, so the estimator fits it no nuisance model | `implementation_roadmap.md` | Which outcomes the diagnostic is actually about. The measurement reproduces; the attribution does not |
| 4 | The `len(fit.alpha)` bullet gains its answer: not degenerate, 4.6% at five cutpoints, and the drop-rate consequence Stage 8 feared does not arise because the drop rate is zero | `implementation_roadmap.md` | What the requested distribution turned out to be, and which half of Stage 8's worry survives |
| 5 | The "Accept when" clause gains the per-outcome-replicate-set rule and the pinned percentile definition | `implementation_roadmap.md` | Two things [§10] left open that an implementer would otherwise decide silently |
| 6 | **An addition, not an amendment**: the entry gains a note that `Σw` moves by a factor of 5.7 across replicates and `polr`'s iteration count runs down to 1 — Stage 8 §11's predicted symptom, measured | `implementation_roadmap.md` | Stage 8 asked for this and the roadmap had nowhere to record the answer |
| 7 | **And one roadmap addition that is not an amendment**: Stage 10 gains its `**Spec:**` line, as Stages 1-9 have | `implementation_roadmap.md` | It was the only stage of 1-10 without one |
| 8 | Two items closed, three rewritten, six added | `TODOS.md` | §16 |

**No SAP amendment.** Every substantive rule above is already [§10]'s; what changed is the roadmap's
transcription of it and two gaps [§10] leaves by silence. The percentile *definition* (§8.2) is the one
place this document comes close, because [§10] makes a claim — that its p-value agrees with its interval
by construction — that is true only under one definition. That is not a defect in [§10]: the claim is
correct, and this document pins the definition that makes it so. Recording it as a `config.py` constant
with the measurement beside it is the [§8]-amendment-free way to keep it correct. The drop-granularity
rule (§7.3), the identifier scheme (§5.2), the augmented estimate's lack of a p-value (§9.2) and S8's
reclassification (§5.4) are recorded as **PI decisions in §17** rather than as protocol changes, because
[§10] and [§8] delegate all four by silence rather than specifying them wrongly — which is Stage 9 §22's
own test for the difference. **Two of the four were put to the PI and confirmed on 2026-08-25** and are
now **DECISION 6** (S8's reclassification) and **DECISION 7** (per-outcome replicate sets) in
`../out/stage0_data_inventory.md`. The other two — the identifier scheme and the augmented estimate's
lack of a p-value — stand as this document's decisions and remain PI-reversible.

### 22.1 The probe round — 2026-08-25

Twelve probes, run before any section was drafted. **Four changed what this document says and two more
changed a number.**

| # | Probe | Finding |
|---|---|---|
| P1 | Reproduce Stage 9 §12.3's blocker | 104/400 = 26.0%, matching to the decimal. Handover confirmed rather than trusted |
| P2 | The renaming fix | 0/400, then 0/2000. §5.2 |
| P3 | The full loop at `N_BOOT` | **Found the second blocker**: S8 raises on 1.0% of replicates and Stage 8 §11 forbids catching it. Neither Stage 8 nor Stage 9 had a frame that could show this. §5.4 |
| P4 | G7 versus non-convergence | Stage 8 §6.1 reproduced exactly. **And `FitError` carries no code**, so Stage 8 §11's "report them separately" needs a prefix map and a scan. §7.2 |
| P5 | `nan` odds ratios and exact-zero risk differences | **The hypothesis was wrong.** 0 `nan` and 0 exact zeros in 1998, because S8 raises first. The planned "`nan` is neither a failure nor a draw" section does not exist; §8.2 records the interaction instead |
| P6 | The percentile definition | **The first attempt measured nothing** — it drew from `normal(≈0, 1)`, where half the draws are ≤ 0 and p ≈ 1, so the boundary was never approached and "0 disagreements" was vacuous. Re-run across the boundary: `linear` disagrees in 4.707%, all at tail count 50. §8.2 |
| P7 | The frozen paths | `ph2` 40.6%, `sich` 2.6%, against Stage 9's 46.1% and 3.8%. Same phenomenon, different stream; §6.3 records both |
| P8 | The S8 rate | 0.8% `sich`, 0.2% `tici_2b_3` — the measurement `TODOS.md` named as the trigger for reclassifying it |
| P9 | The requested distributions | `len(fit.alpha)` **not degenerate** (§7.4); and `Σw` down to 4.90 with iterations down to 1 — Stage 8's predicted symptom, present (§7.5) |
| P10 | Cost and determinism | The per-outcome route is **cheaper** than `secondary` whole, so §7.3's decision needs no cost argument (§2) |
| P11 | Coverage | Both designs bracket nominal. Design C's truth is 0.659272 and not 0.7 — non-collapsibility — and a test using 0.7 reads as an 88% defect (§8.4) |
| C1 | Does `sich` have a nuisance model? | **No.** Stage 9's tail measurement reproduces (41.8% vs 42.7% above 8) but is of a model the estimator does not fit. §10.2, §22 item 3 |
| F1 | Execute this document's own fences (§21b) | **Three defects, all in this document's code.** `ci_min_draws` off by one at level 0.90; `percentile_ci`'s quantile inexact, flipping the limit's sign at the §8.2 boundary; `FAILURE_BUCKETS` missing `G6` with the raise-site count stated as fourteen against sixteen |

**What the probe round did not find, and what it could not.** It did not exercise `resample` on any
frame but this cohort's, so nothing here establishes the stratification is right for [§14a]'s population
(§12.2). It could not measure coverage at `n` = 93 and `B` = `N_BOOT` — 90 000 fits per design is
already 270 s — so §8.4's numbers are about the procedure and not about this study. And it could not test
the one condition that most characterises this cohort: none of the coverage designs is near-separated,
which is what [§13]'s amendment says the [§7] propensity model is actually fitting. §16 item 4 files it
as the gap that matters.

### 22.2 What did not change, which is the more important half

Every handover from Stages 8 and 9 was re-measured rather than carried, and **all but two reproduced
exactly**:

- Stage 8 §6.1's separated-frame witness: 17 iterations, score criterion, counters at zero,
  `exp(β) = 6.469e15`. Identical.
- Stage 9 §12.3's 26.0%: identical, 104 of 400.
- Stage 9 §6.4's correction rates: 16.1 / 13.4 / 2.5% here against 17.0 / 13.0 / 2.3% there.
- Stage 9 §9.3's five-augmented, two-unaugmented split: identical.
- Stage 8 §11's `FitError`/`SchemaError` rule: unchanged, and §5 exists to make it satisfiable.
- Stage 6 §9's `nan`-off-`in_model` rule: unchanged, fourth consumer.

The two that did not: Stage 9's `max|beta|` attribution (§22 item 3, the measurement reproduces and the
attribution does not) and Stage 8's expectation that `Σw` would not move (§7.5, it moved).

---

## GSTACK REVIEW REPORT

| Review | Trigger | Why | Runs | Status | Findings |
|--------|---------|-----|------|--------|----------|
| CEO Review | `/plan-ceo-review` | Scope & strategy | 0 | — | — |
| Codex Review | `/codex review` | Independent 2nd opinion | 0 | — | — |
| Eng Review | `/plan-eng-review` | Architecture & tests (required) | 0 | **not run** | — |
| Design Review | `/plan-design-review` | UI/UX gaps | 0 | — | — |
| DX Review | `/plan-devex-review` | Developer experience gaps | 0 | — | — |
| Probe round | twelve probes before drafting (§22.1) | Whether the handovers hold and the rules are implementable | 1 | issues_found | 12 probes; 4 changed the document, 2 changed a number |
| Fence execution | §21b, after drafting | Executability of the document's own code | 1 | issues_found | 3 defects, all in this document's fences; all fixed and recorded |

**Findings by origin.** Ten of the twelve probe findings in §22.1 came from running the landed code
against frames it had never seen. Five were the author's own errors caught before or just after they
reached a section: P5's wrong hypothesis, P6's vacuous first measurement, and §21b's three fence
defects. **The three fence defects are the ones worth noting**, because two of them — an off-by-one in a
derived constant and an inexact quantile — would have shipped as code and are exactly the class of thing
a document that only *described* its fixtures would have carried into the implementation.

**Verification performed rather than asserted.** §21 carries 48 rows and every one reads `yes, run` or
`yes, read`; §20's Definition of done item 13 forbids a row reading `stated` that a run could have
produced. Two claims in an earlier draft were refuted by their own verification and are recorded as such
rather than deleted: the duplicated-index claim (§3.3) and P5's `nan` hypothesis (§22.1).

**What is still not established.** Coverage at this cohort's size and at `N_BOOT`; coverage under
near-separation, which is this cohort's actual condition; that the loop is right in any sense beyond
coverage; and everything §15 asserts, since the implementation does not exist yet.

**VERDICT: SPEC COMPLETE, ENG REVIEW NOT RUN.** The document is implementable from itself — every
fixture is code, every pin is measured against a fixture in the repository, every number cites a §21 row
— and it has not been read by anyone but its author.

NO UNRESOLVED DECISIONS
