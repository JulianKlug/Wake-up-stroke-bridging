# Stage 9 spec — secondary binary estimators

Implements roadmap Stage 9 [§8]. Section references in brackets are to `statistical_analysis_plan.md`.
Numbers and decisions referenced as DECISION *n* are established in Stage 0 and recorded in
`../out/stage0_data_inventory.md`. Constants named in `SMALL_CAPS` are Stage 1's and live in
`config.py`; `stage1_config_and_data_contract.md` is their specification. The frame this stage
receives is specified by `stage5_cohort_construction.md` §11, the result object it reads by
`stage6_propensity_and_weights.md` §11, and the function it extends by
`stage8_primary_outcome_estimator.md` §13.

**Status.** Written 2026-08-24 against the landed Stages 1-8, revised 2026-08-25 by review round 3
(§22.3), and **amended 2026-08-25 by the implementation itself (§22.4)** — the stage is built, the suite
is green, and five defects in this document's margins were found by building it. Every pin below
reproduced unchanged; §22.4 says what did not. Every number below was produced by running code; none is carried. Ten probes were run before
any section was drafted, and **three of them changed what this document says** — §8.4 (the roadmap's
stated rationale for its own acceptance test does not hold), §9.5 (the guard has to be frozen before
the bootstrap, and the rate that makes it necessary is 46.1%, not a rounding error) and §12.3
(`propensity.fit` raises on 26.0% of stratified replicates, which is a Stage 10 blocker discovered
here). §21 names each. It is the **sole source for the Stage 9 implementation**: everything the
implementer needs is here, and anything not here is not to be invented.

**What round 3 changed, because it is the round that made the previous sentence true.** Rounds 1 and 2
checked citations and claims. Round 3 checked *executability*, and found that the document could not be
implemented from itself: **three of its six fixtures carried pins to ten decimals and no construction**,
so `tau = 0.2252774207` was an instruction to reproduce a number from a frame that appears nowhere.
§15.0 now carries every fixture as code, and **every pin in this document was re-measured against those
constructions** — so the pins changed, and §21 records the new ones. Round 3 also found that
`model.predict` as specified did not reproduce `Fit.p` (§7.3), that the odds-ratio correction can move
an estimate *across* the null rather than toward it (§6.4), and that `secondary`'s signature could not
express the frozen-path contract §9.5 requires of it (§11). §22.3 is the full list; twelve findings
changed behaviour or a number.

**Goal.** Seven binary outcomes — three secondary and four safety — each with a weighted risk
difference, a weighted marginal odds ratio, and, where the minority cell supports a nuisance model, a
**model-assisted** augmented risk difference tilted by `h = e(1 − e)`. Each estimate on its own [§11]
denominator. Point estimates only.

**Not in scope.** Every interval and every p-value (Stage 10 [§10]), the Benjamini–Hochberg correction
within the two families and the subgroup estimates (Stage 11 [§13]), the standardisation analysis
(Stage 12 [§14a]), the reporting layer (Stage 14 [§16]), and the primary ordinal quantity — which
Stage 8 owns and which this stage **does not read** (§0.2).

**One thing this spec settles that no earlier stage could.** Stage 8's `[§11]` denominator was a third
population that on v7 **equalled** the second — 92 of 92 — which Stage 8 §4.1 recorded as "exactly the
condition under which an implementation that never built the mask is green". At Stage 9 the mask bites:
`tici_2b_3` is present on **91** of the 92, so seven denominators are not one denominator, and the
workbook itself witnesses it. The assertion Stage 8 could only write against a synthetic frame has a
live case here (§4.2, §15.2).

**And one finding that arrives with the stage rather than being designed into it.** Stage 8 needed
`POLR_MAX_ABS_BETA` because separation in the ordinal fitter converges silently, and it could calibrate
that bound because the fits were **bimodal with an empty band** between 8.79 and 18.81. Both halves were
re-asked of `model.firth` here. The first reproduces exactly — `separated_frame()` (§15.0.4a) converges in
**6 iterations** with every safeguard counter at zero and returns `exp(β_treatment) = 441.005397`, while
the unpenalised MLE on the same frame has no maximum at all — allowed to reach it, the Newton route
raises `LinAlgError: Singular matrix`; stopped at statsmodels' default cap it returns a plausible
`max_abs_beta` of 67.87 instead (§9.6). The
second does **not**: `sich`'s `max_abs_beta` across replicates runs continuously from **0.8558 to
80.9825** with 42.7% above 8 and no gap anywhere. So Stage 8's remedy does not transfer, and this stage
prescribes **no bound** — for a structural reason rather than a measured one, which §9.6 states and §16
files.

---

## 0. Where Stage 9 sits

```
  cohort.build(...)          -> DataFrame, 93 x 33, one row per patient      [Stage 5 §11]
  propensity.fit(cohort, a)  -> Propensity: e, w, in_model (92)              [Stage 6 §11]
  balance.assess(...)        -> Balance   — NOT read here                    [Stage 7 §11]
  outcome.primary(...)       -> Primary   — NOT read here                    [Stage 8 §11]

  +------------------------------------------------------------------------------------+
  |  outcome.secondary(cohort, ps, audit)                                              |
  |                                                                                    |
  |    _assert_secondary_inputs(df, ps)          S1-S5, collected            (§4.4)    |
  |                                                                                    |
  |    for each key in C.BINARY_OUTCOMES:            seven, registry order   (§4.1)    |
  |        in_estimate = ps.in_model & df[key].notna()                       (§4.2)    |
  |        minority    = min(events, non-events)                             (§9.1)    |
  |                                                                                    |
  |    _record_estimable(...)      one entry, seven rows -- BEFORE any        (§10.2)   |
  |                                estimate, so a raise below still logs             |
  |                                                                                    |
  |    for each key in C.BINARY_OUTCOMES:                                              |
  |        _assert_estimable(key, y, a, w)       S6-S8, per outcome          (§4.4)    |
  |        share  = weighted_proportion(y, a, w)      STAGE 8'S, reused      (§5.1)    |
  |        rd     = weighted_rd(y, a, w)              share differenced      (§5.2)    |
  |        or_    = marginal_odds_ratio(y, a, w)      + the [§8] guard       (§6)      |
  |                                                                                    |
  |        if _augmentable(key, minority):                                   (§9.2)    |
  |            fit    = outcome_model(df, in_estimate, key)   model.firth    (§7)      |
  |            m1, m0 = _counterfactuals(fit, X)              model.predict  (§7.3)    |
  |            aug    = augmented_rd(y, a, w, h, m1, m0)                     (§8)      |
  |                                                                                    |
  |    _record_estimates(...)      one entry, seven rows                     (§10.3)   |
  |    _record_models(...)         one entry, one row per augmented outcome  (§10.4)   |
  +------------------------------------------------------------------------------------+

                       -> Secondary: estimates, a dict of seven BinaryEstimate  (§3)

  Stage 10 [§10] refits ALL of it, N_BOOT times, and cannot yet                (§14)
  Stage 11 [§13] groups by family for Benjamini-Hochberg                       (§12)
  Stage 14 [§16] reports each estimate with its reduced specification beside it (§12)
```

### 0.1 Why this extends `outcome.py` and adds no module

Stage 8 §0.1 already decided it, and the decision is honoured rather than revisited: *"Stage 9's binary
estimators extend that module so the weighted per-arm proportion has one home."* [§11] requires the
denominator of every estimate, and `weighted_proportion` (`outcome.py:146-164`) **is** the denominator.
A second implementation in a second module would be a second definition of it. The same move Stage 6 §3
made with `ess` and Stage 7 §0.1 made with `smd`.

The split that *is* made is the one Stages 6 and 8 both made: the **numerics** go in `model.py` and the
**[§8] specification** stays in `outcome.py`. Stage 9 adds exactly one function to `model.py` —
`predict` (§7.3) — because `model.firth` returns fitted probabilities at the *observed* treatment and
the augmentation needs both counterfactuals. `predict` names neither the treatment nor any covariate
list, so `model.py`'s docstring rule survives it; the treatment-column-setting stays in `outcome.py`,
where the exposure is already named.

### 0.2 What Stage 9 does not touch

- **`Primary`.** Not read. Stage 8 §11 states it: *"the binary estimators are their own estimates on
  their own [§11] denominators, and the primary ordinal quantity is not an input to any of them."*
  §15.11 asserts it by scan, for Stage 7 §14.11's reason — not because reading it would fail, but
  because it would succeed.
- **`Balance`.** Not read, not imported, and asserted not to be (§15.11). Stage 7 §11 says Stage 8 does
  not read it; the same holds one stage on, and DECISION 4's requirement that the residual imbalance be
  named beside the primary estimate is Stage 14's obligation and not a data dependency here.
- **`model.firth` and `model.design`.** Called, never edited. `firth` is Stage 6's and is used here
  exactly as `propensity.fit` uses it. **`model.polr` is neither called nor edited** — it is Stage 8's
  ordinal fitter and no binary estimate goes near it, which is recorded so that nobody wires it in
  looking for a use. §15.12 asserts the Stage 6 and Stage 8 suites pass unedited.
- **`data.KINDS`.** Nine, unchanged. Stage 9's three entries render under the existing `model` kind, as
  Stages 7 and 8's did. This is the third stage in a row to add audit entries without touching
  `data.py`, which `data.py:151-152` predicted.
- **The propensity model.** Not refit, not reweighted, not reread beyond `e`, `w` and `in_model`.

---

## 1. Deliverables

| Path | Contents |
|---|---|
| `extended_bridging/outcome.py` | **amended.** `Secondary`, `BinaryEstimate`, `secondary`, `weighted_rd`, `marginal_odds_ratio`, `augmented_rd`, `outcome_model`, and **the privates §3.2 lists and counts** — the count is stated there and nowhere else. No existing function changes — asserted, not assumed (§15.12) |
| `extended_bridging/model.py` | **amended.** `predict`, and one sentence in the module docstring (§7.3). No existing function changes |
| `extended_bridging/tests/fixtures_stage9.py` | **new.** §15.0's constructions as code — **seven functions and four constants**, enumerated in §13 and counted there and nowhere else, for §3.2's reason. Every pin in this document is measured against these. The count is SEVEN and an earlier draft of this row said "four" over a list of five (§22.4) |
| `extended_bridging/config.py` | **amended.** `BINARY_OUTCOMES`, and `OR_CONTINUITY` (§6.3, §13) |
| `extended_bridging/tests/test_outcome.py` | **amended.** The §15 sections, banner-commented to match these numbers — **and `secondary_cohort`, the cohort-level fixture of §15.0.7**, which lives here rather than in `fixtures_stage9.py` because it is the one construction this document did not specify (§15.0.7, §22.4 item 5) |
| `extended_bridging/tests/test_model.py` | **amended.** `predict`'s sections and the separated-frame witness (§15.7) |
| `extended_bridging/tests/test_config.py` | **amended.** Five assertions (§13) |
| `extended_bridging/tests/reference/aug_psweight.R` | The R oracle for the augmented estimator, behind the gate Stage 7 built (§19b) |
| `extended_bridging/tests/test_reference_r.py` | **amended.** The Stage 9 half of the oracle |
| `extended_bridging/implementation_roadmap.md` | **amended.** Stage 9 gains its `**Spec:**` line, four corrections and one addition (§22) |
| `TODOS.md` | **amended.** Four items with their triggers (§16) |

`out/logs/audit_<label>.md` is an **output, not a deliverable**: it is written only when a caller asks,
and it is gitignored. It is the only place the estimates exist (§4.3).

**Nothing under `specs/` may quote a case identifier, and nothing here does.**

---

## 2. Environment

Unchanged from Stage 8: `uv`, Python 3.12.12, pandas 2.3.3, numpy 1.26.4, statsmodels 0.14.6. `scipy`
is test-only by policy. No dependency is added — `model.firth` and `model.design` are already here, and
`predict` is an inverse logit.

```bash
cd extended_bridging && uv sync && uv run pytest -v
```

**Cost, measured (§21) — and the first draft of this section got it wrong twice.** Medians of 20-30
runs, per augmented outcome, split into the two calls `outcome_model` makes:

```
  outcome         model.design   model.firth
  mrs_0_2_90d          3.81 ms       1.23 ms
  mrs_0_1_90d          3.78 ms       1.73 ms
  tici_2b_3            2.34 ms       0.89 ms
  death_90d            3.68 ms       1.40 ms
  mrs_5_6_90d          3.78 ms       1.40 ms
  ------------------------------------------
  the five augmented  17.40 ms       6.65 ms      -> 24.05 ms, and design is 72% of it
```

**`model.design` is the cost, not the Firth fit.** An earlier draft of this section timed
`model.firth(X, y)` alone, got 1.08-3.09 ms per outcome, and reported ~8 ms for `secondary` — which
omitted the design construction that `outcome_model` performs on every call and which turns out to
dominate it three to one (§21).

In context, per replicate:

```
  propensity.fit   13.39 ms          secondary, the 5 augmented   24.05 ms
  primary           4.73 ms          -------------------------------------
                                     fitting, per replicate      ~42.2 ms
  measured wall-clock in a replicate loop:  71 ms
    the difference is resampling, .loc slicing and frame construction, not fitting
```

**Stage 9 is the bootstrap's cost centre** — 24.05 ms of the 42.2 ms — and Stage 10 should size its
expectations from here and not from Stage 8 §2's *"1.157 ms; one fit"*, which is one `polr` fit and not
one `primary` call. Over `N_BOOT` = 2000 the measured loop is about **142 seconds**.

**And there is an optimisation Stage 10 should take, stated here because it is a fact about these
designs.** The four full-list outcomes are all fitted on the same 92-row population with the same
[§6] covariate list, so `model.design` returns the **identical** matrix four times — about 11 ms of the
17.40. Building it once per replicate and reusing it for those four is worth roughly 22 s over
`N_BOOT`. `tici_2b_3` cannot share it: its population is 91 rows and its covariate list is the override
(§4.2, §7.2). This is a Stage 10 change, not a Stage 9 one — `secondary` computing a design per outcome
is what makes each estimate independent of the others, and §15.12 asserts that independence.

---

## 3. Module shape

`outcome.py`'s imports do not change. `model.py` gains none: `predict` uses `numpy` and `pandas`, both
already imported.

Every `python` fence in this document is valid Python. Stage 6 §3 states the rule, Stage 7 §18d records
that the fence map must be derived by content and not by index, and §21b records this document's own
extraction.

### 3.1 What the stage returns

```python
@dataclass(frozen=True)
class BinaryEstimate:
    """One [§8] binary outcome's three estimates, on its own [§11] population.

    **Four fields move as a block and `reduced` is NOT one of them.** `augmented`, `covariates`,
    `fit` and `dropped` are `None`/`None`/`None`/`()` exactly when `augmented_path` is
    `"unaugmented"`, and all four are populated otherwise: an outcome is either augmented or it is
    not, and §9.2 is the only thing that decides which. A reader must not have to infer the path from
    which fields are absent, so `augmented_path` states it as one of three declared strings and
    §15.9 asserts the block.

    `reduced` is **orthogonal to the path**, and an earlier draft of this docstring listed it in the
    block, which is false on four of the seven outcomes: `reduced` is `True` exactly when the key is
    in `OUTCOME_MODEL_OVERRIDES`, so on v7 it is `True` for `tici_2b_3` alone and `False` for the four
    `full`-path outcomes *whose other three fields are populated*. A block-move test written from the
    earlier sentence fails on `mrs_0_2_90d`, `mrs_0_1_90d`, `death_90d` and `mrs_5_6_90d`. §15.9
    asserts the block and `reduced` separately, and §22.3 item 4 records the correction.

    **There is no verdict and no `favours_bridging()`**, for Stage 8 §3's reason: [§8] specifies no
    threshold on a risk difference, the null is tested at Stage 10 from the bootstrap distribution,
    and a method here returning `rd > 0.0` would be a significance test with no interval behind it.

    `or_corrected` means A WEIGHTED PROPORTION REACHED 0 OR 1, which is very nearly but not exactly
    "a cell was empty". The guard tests the float, so a cell carrying a weight small enough that the
    proportion rounds to exactly 1.0 sets the flag while the cell is not empty -- measured, a
    non-event weight of 1e-18 gives p1 == 1.0 and `or_corrected` True with an uncorrected odds ratio
    of +inf. One-directional: an empty cell always gives exactly 0.0 or 1.0, so there are no false
    negatives, and the correction is right in both cases because +inf is what it exists to avoid.
    §10.3's column therefore reads "corrected" and not "empty cell", and §15.5 asserts the
    distinction rather than leaving a reader to infer a cell count from a flag.

    `or_corrected` is a FLAG AND NOT A MODE. [§8] gives no degenerate-cell rule and the roadmap says
    only "kept finite"; §6.3 prescribes the correction and requires that it be printed wherever the
    estimate is. The two uncorrected proportions travel with it in `proportion` so that a reader can
    see which cell was empty — the correction is never the only record that it fired.
    """

    outcome: str                             # the C.OUTCOMES key
    family: str                              # "secondary" | "safety" — Stage 11 groups on this
    minority: int                            # min(events, non-events) on in_estimate — §9.1
    rd: float                                # the weighted risk difference — §5.2
    odds_ratio: float                        # the weighted marginal odds ratio — §6
    or_corrected: bool                       # §6.3's correction fired
    proportion: dict[int, float]             # the two weighted proportions, UNCORRECTED — §5.1
    augmented_path: str                      # "full" | "reduced" | "unaugmented" — §9.2
    augmented: float | None                  # tau — §8; None when unaugmented
    covariates: tuple[str, ...] | None       # m_a(X)'s declared list, treatment excluded — §7.1
    reduced: bool                            # the list is an OUTCOME_MODEL_OVERRIDES entry — §7.2
    dropped: tuple[str, ...]                 # design columns dropped as constant — §7.1
    in_estimate: pd.Series                   # boolean, TOTAL — the [§11] population — §4.2
    fit: model.Fit | None                    # the m_a(X) fit; None when unaugmented


@dataclass(frozen=True)
class Secondary:
    """The seven [§8] binary estimates, keyed by outcome in C.BINARY_OUTCOMES order.

    A wrapper rather than a bare dict for one reason: `by_family` is Stage 11's grouping and [§13]
    applies Benjamini-Hochberg WITHIN each family. A consumer that groups by reading
    `C.OUTCOMES[k].family` itself is a second place the family partition is computed, and [§13]'s
    correction is wrong if the two disagree. This is not a verdict method — it partitions, it does
    not judge.
    """

    estimates: dict[str, BinaryEstimate]

    def by_family(self) -> dict[str, tuple[BinaryEstimate, ...]]:
        """The estimates grouped by [§5] family, each tuple in C.BINARY_OUTCOMES order."""
        out: dict[str, list[BinaryEstimate]] = {}
        for key in C.BINARY_OUTCOMES:
            out.setdefault(self.estimates[key].family, []).append(self.estimates[key])
        return {fam: tuple(v) for fam, v in out.items()}
```

### 3.2 Public surface, and it is seven names

```python
def secondary(df: pd.DataFrame, ps: propensity.Propensity, audit: Audit,
              paths: dict[str, str] | None = None) -> Secondary: ...


def weighted_rd(y: np.ndarray, a: np.ndarray, w: np.ndarray) -> float: ...


def marginal_odds_ratio(y: np.ndarray, a: np.ndarray,
                        w: np.ndarray) -> tuple[float, bool, dict[int, float]]: ...


def outcome_model(df: pd.DataFrame, in_estimate: pd.Series,
                  outcome: str) -> tuple[model.Fit, pd.DataFrame, tuple[str, ...]]: ...


def augmented_rd(y: np.ndarray, a: np.ndarray, w: np.ndarray,
                 h: np.ndarray, m1: np.ndarray, m0: np.ndarray) -> float: ...
```

plus `Secondary` and `BinaryEstimate`. **Each of the three estimators is public and separately
callable**, which is not the shape Stage 8 chose — `primary` exposed one entry point and kept its
pieces private. The reason is [§8]'s own: the three are *separately reported* ("augmented risk
difference, with the unaugmented weighted risk difference and the weighted marginal odds ratio
alongside"), and §15.8's test that a constant outcome model reduces `augmented_rd` to `weighted_rd`
has to call both. A test that can only reach an estimator through a seven-outcome loop over a workbook
is a test of the loop.

`predict` is `model.py`'s and is its sixth public name (§7.3).

Privates in `outcome.py` are `_assert_readable`, `_assert_secondary_inputs`, `_assert_estimable`,
`_tilt`, `_augmentable`, `_augmented_path`, `_counterfactuals`, `_estimable_table`, `_estimates_table`,
`_models_table`, `_estimable_detail`, `_record_estimable`, `_record_estimates`, `_record_models` —
**fourteen**, against the ten Stage 8 added. **This count is stated once, here, and §1 and §13 cite this
sentence rather than repeating the number**, because three earlier drafts of this document carried three
different counts (§22.3 item 5).

Two of the fourteen are round 3's: `_assert_readable` is §4.4's first assertion phase, split out because
a presence check is a precondition of every read that follows it (§4.4a); and `_tilt` is the [§7]
tilting function `h = e(1 − e)`, which had no home and was written inline in `secondary` (§8.2a).
`_or_detail` is **deleted** — it appeared in no section, had no call site and was in the coverage map by
arithmetic only; `_estimable_detail` is the detail-string function §10.2 actually requires.

### 3.3 The numerical facts this stage turns on

```
  * np.average with weights summing to zero RAISES ZeroDivisionError, it does not return nan.
    weighted_proportion already checks the total BEFORE the mean and the order is load-bearing
    [Stage 7 §3.1, Stage 8 §8.1]. Stage 9 calls that function and inherits the check; it does not
    repeat it.

  * A comparison against a missing value is False in numpy, so `y == 1` on a column carrying nan
    counts an absent outcome as a non-event. §4.2's mask is what removes it; there is nothing to
    carry ON THIS PATH and this line is the statement that the mask is the guard.

  * model.Fit.p is the fitted probability at the OBSERVED treatment, because the X handed to `firth`
    carries the observed treatment column. It is NOT m1 and it is NOT m0, and an implementation that
    reads it for either is wrong on every row whose treatment differs from the counterfactual being
    asked for -- which is every row, for one of the two. §7.3 is why `predict` exists.

  * model.design drops constant columns AND RETURNS THEIR NAMES. On the workbook it drops
    `center_USZ` from every one of the seven designs, because USZ contributes no records to the ATO
    population: the strata are HUG 41, Lugano 31, CHUV 21, summing to the cohort's 93. So a
    counterfactual frame built by copying the design and overwriting the treatment column is the
    right shape, and one built from the covariate list is not.

  * Firth's penalty does not turn separation into a failure. Measured on separated_frame() (§15.0.4a):
    converges in 6 iterations, rescales 0, halvings 0, exp(beta_treatment) = 441.005397, fitted
    probabilities strictly interior. The UNPENALISED MLE on the same frame does not converge at all --
    the Newton route raises LinAlgError: Singular matrix AT maxiter >= 200, and BFGS reaches
    max_abs_beta 33.0927 at 60 iterations while REPORTING SUCCESS. So the penalty does not shrink a
    large coefficient; it turns a fit with no maximum into one that converges quietly, and hides the
    separation completely [Stage 8 §6.1, one fitter over].

  * AND THE ITERATION CAP IS PART OF THAT MEASUREMENT, which an earlier draft omitted. statsmodels'
    Logit defaults to maxiter 35, and at that cap the Newton route does NOT raise: it returns
    max_abs_beta 67.87 with converged False. At 100 it returns 67.45, at 200 it raises. Three caps,
    three different answers, none of them a maximum -- which is a SHARPER statement of "no maximum
    exists" than the raise alone, because a raise reads as a numerical accident and a coefficient that
    moves with the cap cannot. BFGS at 60 is worse still: converged True at 33.0927, a REPORTED
    SUCCESS at a value that is a function of the cap. §15.7 asserts all three [§22.4 item 3].

  * 1/(1+exp(-eta)) is the better-conditioned inverse logit and it is NOT unconditionally safe.
    At the POSITIVE tail it is safe where exp(eta)/(1+exp(eta)) is not: measured, at eta = 500 the
    two agree at 1.0; at eta = 800 and 1000 the first is 1.0 and the second is nan. At the NEGATIVE
    tail it OVERFLOWS: at eta = -800 it emits RuntimeWarning: overflow encountered in exp and
    returns exactly 0.0, where the landed fitter returns 7.124576406741285e-218. An earlier draft
    of this line said the form is "safe at BOTH tails"; it measured one tail. §22.3 item 9.

  * SO THE CLIP IS NOT OPTIONAL, AND IT IS NOT STAGE 8'S. model._probabilities (model.py:320) is
    1.0/(1.0+np.exp(-np.clip(Xc @ beta, -C.FIRTH_ETA_CLIP, C.FIRTH_ETA_CLIP))) -- the SAME form,
    clipped, at FIRTH_ETA_CLIP = 500.0 (config.py:310). Its docstring gives the reason and the
    reason is load-bearing (model.py:317-318): "if this applied a different clip -- or none -- from
    the one inside the loop, F5 would be testing a quantity the weights are not computed from, and
    a fit whose linear predictor hit FIRTH_ETA_CLIP could pass F5 while producing a weight of
    exactly 0.0." An earlier draft attributed the clip to Stage 8's POLR_ETA_CLIP and concluded
    `predict` needed none. Both halves were wrong. §7.3 clips at FIRTH_ETA_CLIP and §15.6 asserts
    the consequence: predict(fit, X) == fit.p BIT FOR BIT.

  * AND THE INTERCEPT MUST BE PREPENDED AS A COLUMN, NOT ADDED AS A SCALAR. `firth` computes
    eta = Xc @ beta over an intercept-prepended Xc. `beta[0] + X @ beta[1:]` is the same
    mathematics and a different floating-point summation order, so it does not reproduce Fit.p:
    measured, max|difference| 1.11e-16 with 10 of 60 rows differing. §15.6's oracle is an EQUALITY
    and an implementation that adds the intercept as a scalar fails it. §7.3's fence is the form
    that passes.
```

**`config.RARE_MINORITY_THRESHOLD` is read, never compared against a literal.** It is 10
(`config.py:277`), declared at Stage 1 for this stage and recorded by Stage 8 §12 as *"Stage 9's, not
this stage's: the primary outcome is ordinal and has no minority cell"*. §9 is the first read.
```

---

## 4. The populations [§11]

### 4.1 Seven outcomes, computed from the registry and never listed

`C.OUTCOMES` holds eight entries and exactly one is primary. The other seven are `kind == "binary"`,
across the two [§5] families:

```
  secondary   mrs_0_2_90d    mRS 0-2 at 90 days          higher_is_better = True
              mrs_0_1_90d    mRS 0-1 at 90 days          higher_is_better = True
              tici_2b_3      TICI 2b-3                   higher_is_better = True

  safety      sich           Symptomatic ICH             higher_is_better = False
              ph2            Parenchymal haematoma t2    higher_is_better = False
              death_90d      Death at 90 days            higher_is_better = False
              mrs_5_6_90d    mRS 5-6 at 90 days          higher_is_better = False
```

The set is **computed**, as `PRIMARY_OUTCOME` is (`config.py:471-481`): `BINARY_OUTCOMES` is a §13
addition, so `outcome.py` cannot write seven keys as literals and cannot drift from the registry. A
ninth outcome added to `OUTCOMES` is estimated by this stage in the same edit and cannot fail to be.

**`higher_is_better` is carried into no arithmetic here.** Every risk difference is
`P(event | bridging) − P(event | EVT alone)`, in that order, on every outcome — the sign convention is
the *contrast's*, not the outcome's, and orienting the sign per outcome would make a table of seven
numbers in which the meaning of "positive" varies by row. Stage 14 [§16] is where a direction is
attached to a number, and `higher_is_better` is what it reads. §15.3 asserts the sign convention holds
uniformly and §16 records the alternative not taken.

### 4.2 Seven denominators, and on this workbook one of them differs

Each estimate's [§11] population is that outcome's own: `ps.in_model & df[key].notna()`. Measured:

```
  outcome         family      n [§11]   lost   events   non-events   minority
  mrs_0_2_90d     secondary        92      0       42           50         42
  mrs_0_1_90d     secondary        92      0       24           68         24
  tici_2b_3       secondary        91      1       85            6          6
  sich            safety           92      0        5           87          5
  ph2             safety           92      0        9           83          9
  death_90d       safety           92      0       25           67         25
  mrs_5_6_90d     safety           92      0       28           64         28

  cohort 93  ->  in_model 92 [Stage 6 §4.4]  ->  seven populations, six of 92 and one of 91
```

**This is the assertion Stage 8 could not make against the workbook.** Stage 8 §4.1 built the outcome
mask, found the resulting population equal to `in_model`'s 92, and recorded that as *"exactly the
condition under which an implementation that never built the mask is green"*. Here `tici_2b_3` is
present on 91, so an implementation that ranges over `in_model` instead of over the mask computes
TICI's weighted proportions over a row whose outcome is absent — where `y == 1` is False on a missing
value (§3.3) and the record is silently counted as a non-event, moving the minority cell from 6 to 7
and the denominator from 91 to 92. §15.2 is the test, and on v7 it fails on a wrong implementation
rather than passing vacuously.

`in_estimate` is boolean and **total**, indexed like the cohort frame, exactly as `Propensity.in_model`
and `Primary.in_estimate` are, and for Stage 6 §9's reason: the deliberateness of an absence lives in a
mask, never in a value.

### 4.3 Why no estimate is in this document

Stage 8 §4.3 spent forty-five lines arguing that `β`, `exp(β)` and every `RD_k` stay out of `specs/`,
and reached the position that outcome-blindness ends *after* the estimator is specified rather than
during. Stage 9 inherits that position and **narrows one part of it**, because the situation is not
identical and inheriting the argument unexamined would misstate the repository:

- **What is already public.** `../out/stage0_data_inventory.md:202` is committed and quotes *crude*
  arm contrasts for two of these seven — `mrs_0_2_90d` at +0.069 favouring bridging and `death_90d` at
  +0.086 favouring EVT alone — and names Stage 9 as what will estimate them properly. So for two
  outcomes the direction is on the record before this document is written, and pretending otherwise
  would be a fiction.
- **What is nonetheless absent here.** No *weighted* proportion, no risk difference, no odds ratio, no
  `tau`. The counts in §4.2 are marginal — they do not decompose by arm — and every table below that
  needs a value uses a synthetic frame. The regression pin is the §15.0.2 golden vector, computed on a
  synthetic construction, exactly as Stage 8's was.
- **What that costs.** The same cost Stage 8 named: no pinned regression number on any reported
  estimate anywhere in git. `TODOS.md` carries the trigger (§16).

The rule, stated once: **a count may appear here; a weighted quantity may not.**

### 4.4 The preconditions, in three phases

Three collections, not two. Stage 7 §4.5a's reason gives the frame-level/per-outcome split; **round 3
adds the split *within* the frame-level phase**, and §4.4a is why.

`_assert_readable`, first, collecting every failure into one `SchemaError` and **raising before any
value is read**:

```
  S1  every C.BINARY_OUTCOMES key is a column of `df`
  S5a C.TREATMENT is a column of `df`                           [Stage 4 §4.2]
  S3a ps.in_model, ps.e and ps.w are indexed like df -- SAME LABELS IN THE SAME ORDER
  S3b ps.in_model is boolean and total                          -- see below; this is phase 1's
```

`_assert_secondary_inputs`, second, over what `_assert_readable` passed, collecting into one
`SchemaError`:

```
  S2  every C.BINARY_OUTCOMES column is 0/1-or-missing   -- C.BINARY_COLUMNS' rule, re-asserted on
                                                            the four DERIVED_DICHOTOMIES too, which
                                                            BINARY_COLUMNS does not cover
  S4  ps.e and ps.w are nan off in_model and finite on it       [Stage 6 §9]
  S5b C.TREATMENT is 0/1 and never missing                      [Stage 4 §4.2]
```

**S3b is in phase 1 and an earlier draft of this section put it in phase 2, which reproduced §4.4a's own
defect one check further on.** S4 reads `ps.e[ps.in_model]`. If `in_model` is an integer Series, pandas
treats it as a sequence of **labels** rather than as a mask, so that read raises
`KeyError: "None of [Index([-2, -2, ...])] are in the [index]"` — measured, by running it — before the
collected `SchemaError` naming S3b is ever built. A boolean-dtype mask is a precondition of every masked
read exactly as a column's existence is a precondition of every column read, so it belongs on the
readability side of §4.4a's boundary. This was found by calling `secondary` end to end (§21b) and not by
reading, which is the third time in this document's history that the same rule has had to be re-applied
one check later.

`_assert_estimable(key, ...)`, per outcome, inside the loop:

```
  S6  the [§11] population is non-empty
  S7  both arms are represented in it, and each carries positive total weight
  S8  the outcome is not constant on it
```

**S7 and S8 raise `SchemaError` and not `FitError`, and the distinction is Stage 10's.** A replicate in
which an arm carries no weight, or in which an outcome is constant, is a *sparse replicate* and not a
resampler bug — which argues for `FitError`. The argument the other way wins: `weighted_proportion`
already returns `nan` rather than raising when an arm has no positive weight (`outcome.py:155-157`), so
S7 firing means the mask and the weights disagree about the same rows, which is a bug. A constant
outcome (S8) makes the risk difference 0 and the odds ratio undefined at both ends at once; it is
reachable under resampling for `sich` and it is the one case §16 files as possibly wanting `FitError`
instead. **This is the sharpest open question the stage leaves** (§16, §14).

Every one of S1-S8 is unreachable on the workbook. That is recorded rather than treated as a reason to
skip one — Stage 6 §4.5's position, and Stage 10 is the population of frames they are written for.
**Every one of them also has a test, one frame per branch, in §15.1a.** An earlier draft pointed T7 at
§15.1, which specified a test for S1 alone, so seven guards written for Stage 10 were exercised by
nothing while the coverage map claimed 71 of 71 routes were covered (§22.3 item 6).

### 4.4a Why the frame-level phase splits, and it is a defect one stage back repeating

A collection is only a collection if every check in it can run. S2 reads `df[key]`; S1 is the check that
`df[key]` exists. In one collecting phase, a frame missing an outcome column has S1 record its failure
and then S2 raise a **bare `KeyError`** before the collected `SchemaError` is ever built — so the caller
sees a pandas traceback instead of the message this document wrote, and §15.1's promise that a missing
column *"raises S1 rather than being skipped"* is false.

**This is Stage 7's defect, one stage on.** Its learnings entry states the rule: *"the boundary is
can-this-be-READ vs is-the-DATA-judgeable — NOT mask-checks vs the-rest. Stage 7 split B1/B2 (mask) into
phase 1 and left B3 (column presence) with B4/B5; B5 then read the column B3 had just reported missing
and produced `KeyError(ivt)`, reintroducing the identical defect one check later. A column-presence
check is a precondition of every read that follows it, exactly as a mask check is."*

**And Stage 6 solved it the other way, which is recorded because the two are not equally safe.**
`model._assert_design_inputs` keeps one phase and re-filters inside it (`model.py:152-154`): *"D1 fires
before D3 can KeyError on a covariate, so the generator re-filters on `c in df.columns` — the checks are
collected, so D1 must not prevent D3 and D4 from running."* That is correct and it is correct *by
vigilance*: every future check added to the helper must remember to re-filter. The three-phase split is
correct *by structure* — a check placed in phase 2 cannot see a column phase 1 rejected, because phase 1
has already raised. Stage 7 had Stage 6's precedent available and shipped the defect anyway, which is
the argument for structure over vigilance.

`_assert_readable` therefore raises rather than returning, and §15.1a's companions assert the failure
mode it prevents: with the phase-1 raise removed, each phase-1 branch produces a bare pandas exception
instead of a `SchemaError`.

### 4.4b The three helpers, written out

```python
def _assert_readable(df: pd.DataFrame, ps: propensity.Propensity) -> None:
    """S1, S5a, S3a, S3b: what must hold BEFORE ANY VALUE IS READ. Raises; never returns a verdict.

    Phase 1 of three (§4.4a). Every check here is a precondition of a READ that a later phase
    performs, so this function raising is what makes the later phases' collection honest.

    `.index.equals(df.index)` and NOT a length check or a set comparison, which is S3a's whole
    content: `.loc[boolean_series]` returns rows in the SERIES' order, while `m1` and `m0` come from
    `model.design(df.loc[mask])` in the FRAME's order. A permuted-but-equal index therefore pairs
    each record's weight with a different record's fitted value -- no raise, no nan, a different
    number. `propensity.fit` builds e and w from df so the orders agree there, and nothing in S4
    constrains a Propensity a caller assembles by hand [Stage 9 §12].

    S3b IS HERE AND NOT IN PHASE 2, because S4's `ps.e[ps.in_model]` is a masked read: an integer
    in_model is interpreted as LABELS and raises KeyError before any collected SchemaError exists
    (§4.4). Measured, not reasoned -- it came out of calling `secondary` end to end.
    """
    failures: list[str] = []
    if ps.in_model.dtype != bool:
        failures.append(
            f"S3b ps.in_model has dtype {ps.in_model.dtype} and not bool. Every read below masks "
            "with it, and pandas reads a non-boolean Series as LABELS -- so this raises KeyError "
            "rather than reporting a mask problem [Stage 9 §4.4a].")
    if len(ps.in_model) != len(df):
        failures.append(
            f"S3b ps.in_model has {len(ps.in_model)} entries against the frame's {len(df)}. The "
            "mask is TOTAL [Stage 6 §9]: absence lives in a False, never in a missing row.")
    absent = [k for k in C.BINARY_OUTCOMES if k not in df.columns]
    if absent:
        failures.append(
            f"S1  {', '.join(absent)}: not a column of the frame. Every C.BINARY_OUTCOMES key is "
            "estimated by this stage, so a registry entry with no column is a frame that cannot "
            "answer what [§8] asks of it.")
    if C.TREATMENT not in df.columns:
        failures.append(
            f"S5a {C.TREATMENT}: not a column of the frame. Every estimate here is a contrast "
            "between the two declared arms [Stage 4 §4.2].")
    for name, series in (("in_model", ps.in_model), ("e", ps.e), ("w", ps.w)):
        if not series.index.equals(df.index):
            failures.append(
                f"S3a ps.{name} is not indexed like the frame. `.loc` returns the SERIES' order, so "
                "a permuted index pairs each record's weight with another record's fitted value and "
                "returns a different number without raising [Stage 9 §4.4b].")
    if failures:
        raise C.SchemaError("\n".join(["outcome.secondary cannot read this frame:", *failures]))


def _assert_secondary_inputs(df: pd.DataFrame, ps: propensity.Propensity) -> None:
    """S2, S4, S5b: whether the DATA is judgeable. Collected into one SchemaError.

    Phase 2 of three. Runs only on what `_assert_readable` passed, which is what lets it read
    df[key], df[C.TREATMENT], ps.e and ps.w without re-filtering -- Stage 6's
    `_assert_design_inputs` re-filters instead (model.py:152-154) and §4.4a is why this one does not.

    S2 is re-asserted on the four DERIVED_DICHOTOMIES, which C.BINARY_COLUMNS does not cover: those
    columns are Stage 3's and their 0/1-ness is a property of the derivation rather than of the
    workbook's contract.
    """
    failures: list[str] = []
    for key in C.BINARY_OUTCOMES:
        offside = df.loc[df[key].notna(), key]
        bad = sorted(set(offside[~offside.isin((0.0, 1.0))].tolist()))
        if bad:
            failures.append(
                f"S2  {key}: carries {bad} outside {{0, 1}} and missing. `y == 1` is False on any "
                "other value, so a third level would be counted as a non-event silently "
                "[Stage 9 §3.3].")
    for name, series in (("e", ps.e), ("w", ps.w)):
        inside, outside = series[ps.in_model], series[~ps.in_model]
        if not np.isfinite(inside.to_numpy(dtype=float)).all():
            failures.append(
                f"S4  ps.{name} is not finite on in_model. Every weighted quantity divides by a sum "
                "of these [Stage 6 §9].")
        if not outside.isna().all():
            failures.append(
                f"S4  ps.{name} carries {int(outside.notna().sum())} value(s) OFF in_model. Stage 6 "
                "§9's rule is that a deliberate absence lives as nan and is never filled; a finite "
                "value off the mask is a value some consumer will range over.")
    arm = df[C.TREATMENT]
    if arm.isna().any():
        failures.append(
            f"S5b {C.TREATMENT}: {int(arm.isna().sum())} record(s) carry no arm. Treatment is "
            "assigned for every cohort member by Stage 4 [§4.2]; a missing arm is a bug upstream.")
    bad_arm = sorted(set(arm.dropna()[~arm.dropna().isin([float(c) for c in C.TREATMENT_LABELS])]))
    if bad_arm:
        failures.append(
            f"S5b {C.TREATMENT}: carries {bad_arm} outside the declared arm codes "
            f"{sorted(C.TREATMENT_LABELS)}.")
    if failures:
        raise C.SchemaError("\n".join(["outcome.secondary cannot estimate on this frame:", *failures]))


def _assert_estimable(key: str, y: np.ndarray, a: np.ndarray, w: np.ndarray) -> None:
    """S6, S7, S8, on ONE outcome's [§11] population. Names the outcome in every message.

    Phase 3 of three, and the only one that runs inside the loop -- so its messages carry `key`,
    because "the outcome is constant" against seven outcomes costs a bisect (§15.1a).

    S7 RAISES rather than tolerating the nan `weighted_proportion` would return, and §4.4 is the
    argument: that function returns nan when an arm has no positive weight, so S7 firing means the
    MASK and the WEIGHTS disagree about the same rows, which is a bug and not a sparse replicate.

    S8 is `SchemaError` and §16 item 3 is the open question about whether it should be FitError. The
    reason it is SchemaError is NOT that a constant outcome yields nan -- measured, it yields a
    finite 1.3529411765 (§6.2) -- but that it yields a PLAUSIBLE number for an outcome nobody
    observed an event in, which is worse than a nan a reader would notice.
    """
    if y.size == 0:
        raise C.SchemaError(
            f"S6  {key}: its [§11] population is empty. `in_model & notna()` selected no record, so "
            "there is no denominator and every estimate below would divide by zero.")
    for code in C.TREATMENT_LABELS:
        arm = a == float(code)
        if not arm.any():
            raise C.SchemaError(
                f"S7  {key}: arm {code} ({C.TREATMENT_LABELS[code]}) has no record on its [§11] "
                "population. Every estimate here is a contrast and a contrast needs both arms.")
        if not float(w[arm].sum()) > 0.0:
            raise C.SchemaError(
                f"S7  {key}: arm {code} ({C.TREATMENT_LABELS[code]}) carries zero total weight over "
                f"{int(arm.sum())} record(s). weighted_proportion returns nan rather than raising "
                "here (outcome.py:155-157), so reaching this means the mask and the weights "
                "disagree about the same rows -- a bug, not a sparse replicate [Stage 9 §4.4].")
    if len(set(y.tolist())) < 2:
        raise C.SchemaError(
            f"S8  {key}: constant at {y[0]!r} on its [§11] population of {y.size}. The risk "
            "difference is 0 and both weighted proportions are degenerate in the SAME direction, so "
            "§6.3's correction returns a finite odds ratio near 1 for an outcome with no contrast in "
            "it -- a plausible number rather than a legible failure [Stage 9 §6.2, §16 item 3].")
```

**S7 and S8 raise on the first failure rather than collecting**, unlike the two frame-level phases, and
the asymmetry is deliberate: a frame-level phase reports on a frame the caller can fix in one edit,
where a per-outcome failure means this outcome is unestimable and the remaining six are a separate
question. Collecting across the loop would report seven problems for one bad mask.

---

## 5. The weighted risk difference [§8]

### 5.1 The one weighted per-arm proportion, and it is Stage 8's

`weighted_proportion(x, a, w)` (`outcome.py:146-164`) is called and not reimplemented. Its docstring
already says why: *"Public, and it is Stage 9's as much as this stage's — the [§8] weighted risk
difference is this function differenced and the weighted marginal odds ratio is these two
proportions."* Stage 9 is the consumer that sentence was written for.

Two properties are inherited rather than restated: `nan` where an arm carries no positive weight and
never `0.0`; and the zero-total check *before* the mean, because `np.average` raises
`ZeroDivisionError` when the weights sum to zero rather than returning `nan` (§3.3).

### 5.2 The estimator

```python
def weighted_rd(y: np.ndarray, a: np.ndarray, w: np.ndarray) -> float:
    """The [§8] weighted risk difference: P_w(Y = 1 | bridging) - P_w(Y = 1 | EVT alone).

    `share[_TREATED] - share[_COMPARATOR]` and not `share[1] - share[0]`, for the reason
    `cumulative_rd` gives at outcome.py:190-194: the arm codes are named once at module scope from
    TREATMENT_LABELS, and this is one of the two lines in the module whose SIGN is a function of
    which literal goes first. Two bare integers in a subtraction are the one place a reader checking
    the orientation has nothing to check against.

    `nan` propagates from `weighted_proportion` rather than being caught. An arm with no positive
    weight makes the difference undefined and that is what nan means; S7 is what decides whether
    such a frame reaches here at all (§4.4).
    """
    share = weighted_proportion(y, a, w)
    return share[_TREATED] - share[_COMPARATOR]
```

Two lines, and it is public for §3.2's reason: §15.8's identity test calls it directly against
`augmented_rd`, and a test that can only reach it through a seven-outcome loop over a workbook is a
test of the loop.

---

## 6. The weighted marginal odds ratio, and the cell that makes it infinite [§8]

### 6.1 The rule

    OR_w  =  [p1 / (1 - p1)]  /  [p0 / (1 - p0)]

with `p1`, `p0` the same two weighted proportions §5.1 returns. **Marginal, not conditional**: [§8]'s
first line is *"All estimates are marginal in the overlap population. No covariate adjustment of the
reported effect"*, and this odds ratio is a contrast of two marginal proportions and never a
coefficient from a model. §15.4 asserts it differs from the treatment coefficient of a weighted
logistic fit on the same rows, because the two are the quantities most easily confused and the wrong
one is one `firth` call away in this very module.

### 6.2 The four ways `p1/(1-p1)` fails, and only one of them is [§8]'s

```
  p in (0, 1)   the odds are finite and positive                  -- the ordinary case
  p == 0        odds 0.  UNCORRECTED OR is 0 if p1, +inf if p0    -- §6.3's branch takes it first
  p == 1        odds +inf. UNCORRECTED OR is +inf if p1, 0 if p0  -- §6.3's branch takes it first
  p is nan      an arm carried no positive weight; OR is nan      -- inherited, §5.1, NOT corrected
```

**The middle two rows describe the uncorrected arithmetic and NOT what the function returns**, and the
distinction is round 3's. `marginal_odds_ratio` tests the two proportions *before* it divides, so the
degenerate branch is taken and no `0`, no `+inf` and no `nan`-from-`0/0` is ever emitted from it. An
earlier draft of this table read as though those were return values.

The nan row is **not** corrected and must not be. A missing proportion is not a boundary proportion:
there is no cell to add half of anything to, and a correction there would manufacture an odds ratio for
an arm with no data. The correction applies to the two boundary rows only, and §15.5 asserts nan
survives it.

**`0 / 0` and `inf / inf` are unreachable, and the earlier claim that they are reachable was the stated
reason for S8.** An earlier draft said: *"`0 / 0` and `inf / inf` are both reachable when p1 and p0 are
degenerate in the same direction — a constant outcome, S8's case — and both give `nan` in numpy,
silently. This is why S8 exists as a precondition rather than as a branch here."* Measured on a constant
outcome, the branch fires before any division and the function returns a **finite number**:

```
  all-zero outcome on the population   p1 = 0.0, p0 = 0.0   ->  OR 1.3529411765   corrected True
  all-one  outcome on the population   p1 = 1.0, p0 = 1.0   ->  OR 0.7391304348   corrected True
```

So S8's real justification is the opposite of the one recorded, and it is **stronger**: without S8, a
binary outcome with no events at all — or no non-events at all — on its [§11] population reports a
**plausible finite odds ratio near 1** rather than a `nan` a reader would notice. `nan` is a legible
failure; `1.3529411765` is not. §16 item 3 is re-argued on this premise, because the version it was
argued on is false and the correction cuts the same way: softening S8 to `FitError` makes the estimator
return 1.353 for an outcome nobody observed an event in. §22.3 item 8.

### 6.3 The correction — Haldane–Anscombe, only when it fires, named wherever the estimate is

Neither [§8] nor any amendment specifies a degenerate-cell rule. The roadmap said only that the odds
ratio is *"kept finite when a weighted proportion reaches 0 or 1"* — naming the requirement and not the
method. **That text is superseded by §22 item 7 and is quoted here as it read before this document, so
it carries no line reference: the line now states the rule below.** **PI decision, 2026-08-24, recorded here because this document is where
it was needed.** The rule:

```python
def marginal_odds_ratio(y: np.ndarray, a: np.ndarray,
                        w: np.ndarray) -> tuple[float, bool, dict[int, float]]:
    """The [§8] weighted marginal odds ratio, its correction flag, and the raw proportions.

    Returns THREE things from one call, for Stage 8 §8.1's reason: an odds ratio alone cannot be
    checked against anything, and a reader given only the corrected value cannot see which cell was
    empty. §10.3's table prints all three because of this signature, not despite it.

    The correction is Haldane-Anscombe: OR_CONTINUITY added to all four weighted pseudo-counts --
    the per-arm weight-sums of events and non-events. It fires ONLY when a proportion is degenerate,
    so an interior estimate is bit-for-bit what §6.1 specifies and is not perturbed to buy
    smoothness it does not need.

    **IT DOES NOT NECESSARILY SHRINK TOWARD THE NULL, AND AN EARLIER FORM OF THIS DOCSTRING SAID IT
    DOES.** The four pseudo-counts are WEIGHT-SUMS, and the two arms' weight totals are not equal
    across replicates, so the same additive 0.5 is a different relative nudge in each arm. When the
    arm holding the empty cell is the HEAVIER one, the corrected odds ratio can land on the far side
    of 1 from the uncorrected limit -- measured, an empty treated event cell with Sw1 = 17.007
    against Sw0 = 22.768 and p0 = 0.00647 returns 1.0201 where the uncorrected value is 0.0. §6.4
    carries the measurements and §17 files what it implies for a decision that is PI-reversible.
    What IS guaranteed is only that the result is finite and positive.
    """
    share = weighted_proportion(y, a, w)
    p1, p0 = share[_TREATED], share[_COMPARATOR]
    if np.isnan(p1) or np.isnan(p0):
        return np.nan, False, share
    if min(p1, p0) > 0.0 and max(p1, p0) < 1.0:
        return (p1 / (1.0 - p1)) / (p0 / (1.0 - p0)), False, share
    c = C.OR_CONTINUITY
    odds = {}
    for code in C.TREATMENT_LABELS:
        arm = a == float(code)
        ev = float(w[arm][y[arm] == 1.0].sum()) + c
        non = float(w[arm][y[arm] == 0.0].sum()) + c
        odds[code] = ev / non
    return odds[_TREATED] / odds[_COMPARATOR], True, share
```

`OR_CONTINUITY` is `0.5` and is prespecified in `config.py` for `FIRTH_*`'s reason (`config.py:281-284`):
[§10] refits this in every replicate, so a correction that changes an answer is a property of the
sampling distribution and not a runtime knob (§13).

**The pseudo-counts are weight-sums and not row counts**, which is what makes the correction's scale a
choice rather than a convention. §6.4 is that argument.

### 6.4 What the correction costs, stated rather than minimised

Three measured facts, and the third is the one that matters most:

- **`Σw = 27.736623` over the ATO population** [Stage 8 §5.2, re-measured here], and it splits
  **13.626360 treated / 14.110263 control** — near-even, which is not something to assume, since `w` is
  `1 − e` in one arm and `e` in the other. The row counts are **not** near-even: **39 treated, 53
  control**. The conventional Haldane–Anscombe `0.5` is calibrated against *row counts*, so its
  relative size here is `0.5/13.63 = 0.0367` and `0.5/14.11 = 0.0354` against the `0.5/39 = 0.0128` and
  `0.5/53 = 0.0094` it was designed for — **2.9× more aggressive in the treated arm and 3.8× in the
  control arm**, measured, rather than the single "three times" a symmetric reading would give. The
  correction is not scaled to `Σw` to compensate, because a correction whose magnitude is a function of
  the weights is a correction whose magnitude is a function of the propensity model, and [§10] refits
  that in every replicate.
- **It is unreachable on the point estimate.** Measured: on v7, for every one of the seven outcomes,
  neither arm holds an empty event cell or an empty non-event cell. The branch **never fires** on the
  workbook. This is Stage 8 §4.1's condition again — an implementation that never wrote the branch is
  green on the data — so §15.5 tests it on synthetic frames and the coverage map records it as
  unreachable-on-every-frame-this-stage-has-except-the-constructed-ones.
- **And it is reachable, often, in the bootstrap.** Measured over 400 stratified replicates: the branch
  fires for `sich` in **17.0%**, for `tici_2b_3` in **13.0%** and for `ph2` in **2.3%**. So the choice
  of correction rule is materially a choice about `sich`'s and `tici_2b_3`'s **intervals**, and barely
  a choice about any point estimate. That is the honest description of what was decided, and it is not
  what the roadmap's one clause suggests is at stake.

- **And the asymmetry does not merely widen the estimate, it can invert it.** This is round 3's, and it
  is the consequence the three facts above imply without stating. Because the pseudo-counts are
  weight-sums and the arms' totals differ, the additive `0.5` is a *different relative* nudge per arm;
  when the arm carrying the empty cell is the heavier one, the correction carries the odds ratio across
  1 rather than toward it. Measured, at realistic ATO arm totals:

```
    empty cell   Sw treated   Sw control   p other arm   uncorrected OR   corrected OR   crosses 1
    treated ev      17.007       22.768       0.00647              0.0         1.0201      YES
    treated ev      12.404       20.945       0.01283              0.0         1.0674      YES
    treated ev      13.626       14.110       0.02000              0.0         0.6484      no
```

  The third row is the workbook's own near-even split (§6.4 first bullet) and does **not** cross; the
  crossing needs the empty arm to be the heavier one, which resampling produces. Over ATO-shaped
  92-row replicates with an empty treated event cell, the corrected odds ratio exceeded 1 in **0.9%**
  of them. Against §6.4's measured firing rate of 17.0% for `sich`, that is roughly **3 replicates in
  an `N_BOOT` = 2000 bootstrap in which a safety outcome with no bridging events at all is reported as
  favouring EVT alone**. It is a small number and it is the wrong direction on the wrong family, so it
  is stated here and filed in §17 rather than buried. §15.5 asserts it rather than asserting the
  opposite, which is what an earlier draft did.

**The alternative, and why it is recorded rather than adopted.** Returning the odds ratio as `nan` when
a cell is empty is what the rest of this repository's conventions point at: `weighted_proportion`
returns `nan` and never `0.0`, and Stage 6 §9's rule is that a deliberate absence lives in a mask and is
never filled. It was not adopted because the roadmap's instruction is to keep the value finite, and
because of the third fact above — under `nan`, 17% of `sich` replicates would return no odds ratio and
Stage 10 would have a count of undefined replicates in place of an interval. It is filed in §17 as the
PI-reversible option with that consequence attached.

**Why not squeeze the proportions instead.** Clipping the offending proportion into `[ε, 1−ε]` leaves
the odds ratio at the **largest finite value the bound permits**, so an empty haemorrhage cell in one
arm reports as a dramatic odds ratio where Haldane–Anscombe reports a modest one. On a safety outcome
that is the wrong direction to err in. It is also not a named procedure, where Haldane–Anscombe is one
a reader recognises. Neither escapes the tuning constant — under both, the reported value is roughly
proportional to one over the constant — and §17 says so rather than claiming the chosen rule is
principled where the other is arbitrary.

---

## 7. The outcome regression `m_a(X)` [§8]

### 7.1 The rule

[§8]: *"treatment main effect plus the §6 covariate set, linear terms, fitted by Firth logistic. No
interactions, no splines, no cross-fitting — at this treated-arm size flexible nuisance models add
variance rather than robustness, and cross-fitting folds would be too small to serve their purpose."*

```python
def outcome_model(df: pd.DataFrame, in_estimate: pd.Series,
                  outcome: str) -> tuple[model.Fit, pd.DataFrame, tuple[str, ...]]:
    """m_a(X) for `outcome`: a Firth logistic fit of treatment plus its declared covariate list.

    Returns the fit, the design it was fitted on, and the columns `design` dropped as constant. The
    design comes back because §7.3's counterfactuals are built by COPYING it and overwriting one
    column -- never by rebuilding from the covariate list, which would not know what was dropped.

    The covariate list comes from `C.outcome_model_covariates(outcome)` and never from
    C.OUTCOME_COVARIATES directly (config.py:489-503). That accessor defaults to the shared [§6]
    entry and returns an override only for an outcome that declares one, which is what makes the two
    nuisance models unable to drift apart: [§8] says "the cross-reference, rather than a second list,
    keeps the two sets from drifting apart in later revisions".

    Treatment is inserted FIRST and by name from C.TREATMENT. `design` never returns it -- the
    covariate list does not contain it, and [§8] calls it a main effect added to that list -- so
    inserting it here is the whole of "treatment main effect". `firth` prepends the intercept, so a
    design matrix is never ambiguous about whether it has one [Stage 6 §5.2].
    """
    sub = df.loc[in_estimate]
    X, dropped = model.design(sub, C.outcome_model_covariates(outcome))
    X = X.copy()
    X.insert(0, C.TREATMENT, sub[C.TREATMENT].to_numpy(dtype=float))
    return model.firth(X, sub[outcome].to_numpy(dtype=float)), X, dropped
```

**`model.firth` is Stage 6's and raises `FitError`, never a fallback** (`model.py:393`). Its docstring
records that `np.linalg.inv` is deliberately not wrapped, because the pilots' `pinv` fallback *"is a
fallback not to a different estimator but to a different estimand"*. Stage 9 inherits that and adds
nothing: an `m_a(X)` that cannot be fitted is a `FitError` and [§10] drops and counts the replicate.

### 7.2 The one declared override, and the roadmap sentence it contradicts

`C.OUTCOME_MODEL_OVERRIDES` (`config.py:485`) holds exactly one entry —
`"tici_2b_3": ("center", "atrial_fib")` — and `outcome_model_covariates` defaults to the shared list
for every other key. The roadmap's design requirement is met by construction: *"an outcome can only
diverge by being named — never by an edit to the default quietly applying to one outcome."*

**The roadmap and the SAP disagree about what happens to `tici_2b_3`, and the SAP is right.**

- **The roadmap, as it read before §22 item 1**: *"outcomes whose minority cell —
  `min(events, non-events)` — is below 10 are **not augmented**"*, with TICI named as the case the rule
  is written for. Quoted without a line reference because this document replaced it; the line now
  carries the corrected rule.
- `statistical_analysis_plan.md:229-233`, the amendment of 2026-08-07, item 2: *"Where the minority cell
  cannot support the full covariate set but can support a smaller one, the outcome is augmented with a
  **declared reduced `m_a(X)`** rather than dropped from augmentation. For TICI 2b–3 that model is
  treatment + `center` + `atrial_fib`."*

TICI's minority cell is 6 (§4.2). Under the roadmap's sentence it is dropped from augmentation; under
the amendment it is augmented with a five-parameter model. **The SAP is the protocol and the roadmap
sentence is stale** — and the roadmap contradicted *itself* twice over before it contradicted the SAP,
which is what settles it:

- Its own adjacent bullet requires that *"the override must reach the reporting layer and be printed
  beside the TICI estimate"* — an instruction about an augmented TICI estimate that the guard bullet had
  just removed.
- And **Stage 0's DECISION 3, in the same file, already states the rule correctly**: *"the [§8]
  rare-outcome threshold applies to `min(events, non-events)`, and `TICI 2b–3` is augmented with a
  reduced `m_a(X)` — treatment + `center` + `atrial_fib`, 5 parameters against 6 non-events — rather
  than dropped from augmentation. ... `sICH` and `PH2` remain unaugmented."*

So the correction is not this document choosing between two documents. It is one stale sentence in the
Stage 9 section against the SAP, the Stage 0 section and the Stage 9 section's own neighbouring bullet.
§9.2 is the rule as this stage implements it and §22 item 1 is the roadmap amendment.

**And the reduced model's parameter count is contingent on a fact nobody declared.** The amendment says
*"five parameters, chosen on procedural grounds"*. Measured, that is right — and it is right by
accident. `center` has **four** declared levels (`FACTOR_LEVELS["center"] = ('HUG', 'CHUV', 'Lugano',
'USZ')`, reference `HUG`), which is three dummies; `design` then drops `center_USZ` as constant on every
one of the seven designs, because USZ contributes **no records** to the ATO population — the strata are
HUG 41, Lugano 31, CHUV 21. So the fitted design is treatment + 2 centre dummies + `atrial_fib` = 4
columns, plus the intercept `firth` prepends = **5 parameters, measured**, against 6 non-events.

If USZ ever contributes a single record, the same declared model becomes **6 parameters against 6
non-events** — saturated — and the amendment's stated arithmetic silently stops holding. §15.9 asserts
the count against the measured 5 *and* asserts what it becomes on a frame where USZ is present, so the
contingency is visible in the suite rather than in a comment. §16 files it and `TODOS.md` carries the
trigger.

### 7.3 `model.predict`, and why `Fit.p` is not it

`model.Fit.p` is documented as *"fitted probabilities"* aligned positionally to the rows of the `X` it
was given (`model.py:104`). The `X` carries the **observed** treatment column, so `p` is `m_A(X)` — the
counterfactual matching each row's actual arm. The augmentation needs `m_1` and `m_0` **on every row**,
which `p` is for one arm and is not for the other. An implementation reading `p` for either is wrong on
every row assigned to the other arm.

```python
def predict(fit: Fit, X: pd.DataFrame) -> np.ndarray:
    """The fitted probabilities of `fit` on `X`. `X` carries no intercept; one is prepended.

    Named neither for treatment nor for any covariate: this is `Fit` evaluated somewhere, and the
    somewhere is the caller's business [Stage 6 §0.1]. `outcome.py` is what knows that one of the
    columns is an exposure and that setting it to a constant makes the result a counterfactual.

    The column check is BY NAME and is not a length check. `design` drops constant columns and
    returns their names, so a frame built from a covariate list has a different width from the one a
    fit saw -- and a positional dot product against a mismatched design returns a number rather than
    an error. Comparing `fit.columns` to `X.columns` as a SEQUENCE catches a reordering too, which a
    set comparison would not: `firth`'s beta is positional after the intercept.

    **THIS IS `_probabilities` EVALUATED SOMEWHERE ELSE, AND IT IS THE SAME EXPRESSION ON PURPOSE.**
    Two details, both of which an earlier draft of this function got wrong, and both of which §15.6
    now asserts as an EQUALITY against `fit.p` rather than as a tolerance:

      * The intercept is PREPENDED AS A COLUMN, never added as a scalar. `firth` computes
        `eta = Xc @ beta` over an intercept-prepended `Xc`; `beta[0] + X @ beta[1:]` is the same
        mathematics in a different summation order and does not reproduce `Fit.p` -- measured,
        max|difference| 1.11e-16 on 10 of 60 rows (§3.3).
      * The linear predictor IS CLIPPED, at FIRTH_ETA_CLIP, because `_probabilities` clips at
        FIRTH_ETA_CLIP (model.py:320) and model.py:317-318 gives the reason: a reader of `Fit` that
        applies "a different clip -- or none -- from the one inside the loop" is reading a quantity
        the fit was not computed from. 1/(1+exp(-eta)) is better conditioned than
        exp(eta)/(1+exp(eta)) at the POSITIVE tail only; at eta = -800 it overflows and returns
        exactly 0.0 where the clipped form returns 7.124576406741285e-218 (§3.3). The clip for this
        form is FIRTH_ETA_CLIP and has nothing to do with Stage 8's POLR_ETA_CLIP.
    """
    if tuple(X.columns) != fit.columns:
        raise C.SchemaError(
            f"predict was given columns {tuple(X.columns)} for a fit on {fit.columns}. "
            "The coefficient vector is positional after the intercept, so a mismatched or "
            "reordered design returns a number instead of an error [Stage 9 §7.3].")
    Xc = np.column_stack([np.ones(len(X)), X.to_numpy(dtype=float)])
    eta = np.clip(Xc @ fit.beta, -C.FIRTH_ETA_CLIP, C.FIRTH_ETA_CLIP)
    return 1.0 / (1.0 + np.exp(-eta))
```

`model.py` already imports `config as C` and `numpy as np`, so `predict` adds no import — §3's claim
survives the correction. Measured on both golden-frame fits (§15.0.2): `predict(fit, X)` is
`np.array_equal` to `fit.p`, exactly, on the very design the fit was made from.

The counterfactual frames are built in `outcome.py`, by copying the design `outcome_model` returned:

```python
def _counterfactuals(fit: model.Fit, X: pd.DataFrame) -> tuple[np.ndarray, np.ndarray]:
    """m_1 and m_0 on every row of `X`: the design with the exposure set to each declared arm.

    A COPY with one column overwritten, never a rebuild from the covariate list, because `design`
    dropped constant columns and the list does not know which (§7.1). Every other column keeps its
    observed value: [§8]'s m_a(X) is the conditional mean given the covariates AT their observed
    values, and only the exposure is intervened on.
    """
    out = []
    for code in (_TREATED, _COMPARATOR):
        frame = X.copy()
        frame[C.TREATMENT] = float(code)
        out.append(model.predict(fit, frame))
    return out[0], out[1]
```

One sentence is added to `model.py`'s module docstring, for the reason Stage 8 §12 added one: the
docstring says the module *"names neither the treatment nor any covariate list"*, and a reader who finds
`predict` being used to build counterfactuals will want to know whether the rule was bent. It was not —
`predict` cannot tell an exposure from a covariate — and the sentence says so where the rule lives.

### 7.4 The seven fits, measured

Every one of the seven converges on the workbook. None raises. None triggers a safeguard:

```
  outcome         path          rows  cols  iters   max|beta|  rescales  halvings   dropped
  mrs_0_2_90d     full            92    12      6      3.0933         0         0   center_USZ
  mrs_0_1_90d     full            92    12      9      1.0268         0         0   center_USZ
  tici_2b_3       reduced         91     4      7      2.2840         0         0   center_USZ
  sich            unaugmented     92    12     17      6.6522         0         0   center_USZ
  ph2             unaugmented     92    12     11      2.0272         0         0   center_USZ
  death_90d       full            92    12      7      3.8288         0         0   center_USZ
  mrs_5_6_90d     full            92    12      7      4.2117         0         0   center_USZ
```

`sich` and `ph2` are fitted **by the probe and not by this stage** — §9.2 does not augment them, so
`outcome_model` is never called for either. They are in the table because the question §9.6 answers is
what the fitter does on the sparsest designs, and the answer is that the sparsest one takes 17
iterations to reach `max|β| = 6.6522` and reports nothing unusual while doing it.

---

## 8. The augmentation [§8]

### 8.1 The estimator

```
  tau  =  Σ₁ w(Y − m₁)/Σ₁ w  −  Σ₀ w(Y − m₀)/Σ₀ w  +  Σ h(m₁ − m₀)/Σ h        with h = e(1 − e)
```

```python
def augmented_rd(y: np.ndarray, a: np.ndarray, w: np.ndarray,
                 h: np.ndarray, m1: np.ndarray, m0: np.ndarray) -> float:
    """The [§8] MODEL-ASSISTED augmented risk difference. Never "doubly robust" (§8.3).

    `h` is passed in and is not computed here, because a tilting function is not this function's to
    choose: `augmented_rd` is arithmetic over six arrays and §8.2a's `_tilt` is the one place
    h = e(1-e) is written. An earlier draft of this docstring said h "is `Propensity`'s", which is
    false -- `Propensity` carries e, w and in_model and no h -- and the falsehood mattered, because
    it let the expression live inline at one call site with nothing asserting it (§8.2a).

    The three terms are summed in the order [§8] writes them. Two arm terms, each normalised by its
    OWN arm's weight total, then one correction term normalised over BOTH arms -- h is not
    arm-specific, so its denominator is the whole [§11] population and not either arm's. An
    implementation normalising the correction by an arm total is caught by §15.8 item 1, not by item
    4: measured on golden_frame(), the arm-normalised form returns tau - rd = +0.0812746102 on the
    two-constant test, against item 1's tolerance of 1e-15 -- thirteen orders of magnitude past it,
    and a decisive failure. An earlier draft of this docstring cited "§15.8's third test" for this, which
    is the tilt test, and §15.8 item 4 turns out to test nothing item 1 does not (§22.3 item 13).

    THE THREE DENOMINATORS ARE CHECKED BEFORE THEY DIVIDE, for `weighted_proportion`'s reason
    (§3.3): this function is public and separately callable (§3.2), so S7 is not between it and a
    caller, and a raw divide returns nan with a RuntimeWarning rather than saying which total was
    empty. An earlier draft divided raw and §5.1's claim that Stage 9 "inherits the check" was true
    of `weighted_rd` and false of this function (§22.3 item 12).
    """
    t, c = a == float(_TREATED), a == float(_COMPARATOR)
    totals = {"treated": float(np.sum(w[t])), "control": float(np.sum(w[c])),
              "tilt": float(np.sum(h))}
    empty = sorted(k for k, v in totals.items() if not v > 0.0)
    if empty:
        raise C.SchemaError(
            f"augmented_rd was given zero total weight for: {', '.join(empty)}. "
            f"Totals are {totals}. np.sum over an empty selection is 0.0 and the division would "
            "return nan with a warning rather than naming which term was undefined "
            "[Stage 9 §8.1]. S7 is the guard on the `secondary` path; this is the guard on the "
            "public one.")
    return (float(np.sum(w[t] * (y[t] - m1[t])) / totals["treated"])
            - float(np.sum(w[c] * (y[c] - m0[c])) / totals["control"])
            + float(np.sum(h * (m1 - m0)) / totals["tilt"]))
```

### 8.2 Why `h` and not `w`

[§8]: *"The augmentation uses the same tilting function `h = e(1 − e)` as the weights, so it targets the
ATO rather than the ATE."* `w` is arm-specific — `1 − e` treated, `e` control [§7] — and `h` is not;
`h` is the ATO's own tilting function and the target population `[§7]` defines. Measured on the
workbook: `Σw = 27.736623` over the ATO population and `Σh = 14.797901`. The two differ by roughly a
factor of two and are **not** proportional, which is what makes §8.4's third test able to separate them
at all.

### 8.2a `_tilt`, and why the ATO's tilting function needs a home

`h = e(1 − e)` is the [§7] tilting function that *defines the estimand*. §8.3 is four paragraphs on
what happens when it is wrong, §8.4 exists because no constant-model test can detect it being wrong,
§15.8 item 3 is described as *"the only test in the suite that pins `h`"*, and the Failure modes table
calls tilting by `w` *"the failure the spec exists to make catchable"*.

**And an earlier draft wrote it as a bare sub-expression at one call site**, `e * (1.0 - e)`, inside
`secondary`'s loop, with `augmented_rd`'s docstring asserting it belonged to a class that does not
carry it. Nothing named it and nothing tested the call site: every test that distinguishes `h` from `w`
calls `augmented_rd` directly with arrays the test constructed, so an implementer typing `w` on that one
line shipped the wrong estimand with the whole suite green. §22.3 item 2.

```python
def _tilt(e: np.ndarray) -> np.ndarray:
    """The [§7] ATO tilting function, h = e(1 - e). One home, one citation, one place to be wrong.

    §0.1's rule applied to the tilting function: `weighted_proportion` is here rather than inside
    `cumulative_rd` because "a second implementation in a later stage would be a second definition
    of the denominator", and Stage 6 §3 and Stage 7 §0.1 made the same move for `ess` and `smd`.
    `h` is the one quantity in [§8] that was left as a sub-expression, and it is the quantity that
    decides WHICH weighted average treatment effect is being estimated -- so it is the one where a
    second definition costs most (§8.3).

    NOT a method on Propensity and not a field of it: `propensity.py` is not amended by this stage
    (§13), and h is read by exactly one consumer. If a second consumer appears, moving it to
    Propensity is the right change and §16 records the trigger.

    Takes `e` and not a Propensity, so it is callable on the raw array §11 already sliced and on a
    synthetic one a test constructs.
    """
    return e * (1.0 - e)
```

Measured on the ATO population, `_tilt(ps.e[in_model]).sum()` is `Σh = 14.797901` against
`Σw = 27.736623` (§8.2) — the two differ by roughly a factor of two and are not proportional, which is
what makes §8.4's third test able to separate them at all. §15.8 item 6 asserts that `secondary`'s
returned `augmented` equals `augmented_rd` recomputed with `_tilt(e)` and **differs** from the same call
made with `w`, which is the assertion that reaches the call site rather than the function.

### 8.3 Not doubly robust, and the word is never used

[§8]'s own paragraph (`statistical_analysis_plan.md:186-200`) is the argument, and this document cites
rather than restates it: the ATO estimand *is itself indexed by the true propensity score*, because the
target population is defined by `h(X) = e(X){1 − e(X)}` and corresponds to no observable subset of the
data. If `e(X)` is misspecified the tilting function and therefore the **estimand** are wrong even when
`m_a(X)` is exactly right. The estimator stays consistent for a weighted average treatment effect — for
the wrong weighting.

**"model-assisted", never "doubly robust", in code, docstrings, audit output and this document**
[roadmap]. §15.13 is a scan of the shipped modules and of the rendered audit text for the forbidden
phrase, because the failure mode is a word in a docstring and no behavioural test can catch it.

The reduced TICI model weakens the bias-reduction argument further, exactly as the amendment says it
does, and the guard is the comparison [§8] already requires: the unaugmented weighted risk difference
is reported alongside, and material disagreement between the two is **evidence about the outcome
model**, not confirmation of either. §10.3's table prints them in adjacent columns for that reason.

### 8.4 The roadmap's first acceptance test is right, and its stated reason for it is wrong

The roadmap required, **as it read before §22 items 2 and 3**: *"With the outcome model set to a
constant, the augmented estimate equals the unaugmented weighted risk difference exactly. If it does
not, the augmentation is not built on the weights' tilting function and is targeting a different
estimand."* No line reference, for §6.3's reason — this document replaced the sentence.

The identity holds. **The inference does not.** For constants `m₁ ≡ c₁`, `m₀ ≡ c₀` and **any**
normalised weight `u`:

```
  Σu(c₁ − c₀)/Σu  =  (c₁ − c₀) · Σu/Σu  =  c₁ − c₀
```

The correction term is independent of `u`. So the identity `tau = p₁ − p₀` holds whether the
augmentation tilts by `h`, by `w`, or by unit weights, and **passing the test establishes nothing about
the tilting function.** Measured, over a 92-record ATO-shaped frame:

```
  c1 = 0.62, c0 = 0.23         tau − rd
    h = e(1−e)                 −2.498e-16
    w                          −1.943e-16
    unit                       −2.776e-17          all three pass, at machine precision
```

Three consequences, and all three are in §15.8 rather than one:

1. **The test is kept, with the correct rationale.** What it checks is that the correction term is
   *normalised* and *cancels* — a real property, and the one an implementation summing `Σ(m₁ − m₀)`
   without a denominator, or normalising an arm term by the wrong total, gets wrong.
2. **The two constants must differ.** With `c₁ = c₀` the test is passed by an implementation that drops
   the correction term **entirely**: measured, `tau − rd = −1.665e-16` for the dropped form, identical
   to the correct one. With `c₁ = 0.62, c₀ = 0.23` the dropped form comes back at `rd − 0.39` — exactly
   `c₀ − c₁`, and a decisive failure. The roadmap's phrasing does not say the constants must differ and
   the obvious reading of "set to a constant" is that they do not.
3. **A third test is needed and it is the one that pins `h`.** With a *non-constant* `m_a(X)`, the
   estimate under `h` and the estimate under `w` differ, and the assertion is against both — the `h`
   value asserted, the `w` value shown and shown to differ. This is Stage 8 §8.2's move: there, the
   empirical `RD_k` was asserted *and* the six-threshold-model route was computed and shown to differ,
   because an assertion against one number cannot distinguish two routes that agree on it.

**The third test's power is a property of the fixture and has to be chosen by measurement, not by
reasoning.** The gap between two tilts is `Σu·d/Σu − Σv·d/Σv` with `d = m₁ − m₀`, which is zero exactly
when `d` is constant — so the power comes from the *variation* in `d`, i.e. from treatment-effect
heterogeneity on the probability scale. Reasoning stops there; measurement does not agree with the
obvious next step:

```
  fixture                                  sd(d)   range(d)     |tau_h − tau_w|
  plain logistic, beta_x = 0.5            0.0205     0.0854          2.201e-04
  plain logistic, beta_x = 2.0            0.0680     0.1955          1.037e-03
  wide x (scale 3), beta_x = 2.0          0.0647     0.1970          1.738e-03
  interaction 1.5                         0.1915     0.5435          5.105e-03
  interaction 3.0, wide x                 0.1807     0.5472          9.416e-03
  interaction 3.0, wide x, steep ps       0.1854     0.5487          2.459e-04   <-- large sd, no gap
```

The last row has the largest `sd(d)` in the table and the second-smallest gap. **The gap is not monotone
in the heterogeneity**: it depends on how the *difference between the two normalised tilts* lines up
against `d`, and a steeper propensity model can align them. So a fixture chosen on the argument "make
the effect heterogeneous and the tilts must separate" can land at `2.5e-04` and produce a test that
passes on the wrong tilting function within any tolerance loose enough to survive a numpy version bump.

**That table is the round-1 exploration and is retained for its argument, not for its numbers.** It was
measured on six ad-hoc constructions none of which this document specified — which is the defect §15.0
now fixes. The fixture §15.8 actually uses is `tilt_frame` (§15.0.4), which is the
`interaction 3.0, wide x` shape written out as code with a pinned seed, and **re-measured against that
construction the gap is `1.030664e-02`**, not the `9.416e-03` an earlier draft carried. The
non-monotonicity argument is unaffected: it is algebra about `Σu·d/Σu − Σv·d/Σv`, not a property of any
one frame.

And the sanity rows, since the algebra above predicts them exactly — re-measured on `tilt_frame`, with
`d` constant at 0.00, 0.05 and 0.40 the gap is `0.000e+00`, `0.000e+00` and `0.000e+00`. (The earlier
draft's `5.551e-17` on the third row was that frame's rounding, not a property of the identity; on this
construction all three cancel exactly.)

### 8.5 The second acceptance test, and the size it has to be run at

The roadmap, **as it read before §22 item 5**: *"A heterogeneous-effect scenario with a correct outcome
model and a misspecified propensity model does not recover the true ATO. This must be asserted, not
merely allowed."* No line reference; §22 item 5 rewrote it.

§8.3's argument predicts something sharper than non-recovery, and the sharper version is what §15.10
asserts: the estimator should recover the ATO indexed by the **misspecified** score, accurately, while
missing the one indexed by the true score. Three quantities, on a population where the outcome model is
the **oracle** — the true conditional means, not a fit, which is the strongest available form of
"correct" and removes nuisance-estimation error from the result:

```
  n = 200,000, seed 20260825, propensity model omits x2, effect modified by x2   [§15.0.5]

    ATO(e_true), the [§8] estimand                     +0.320705
    tau_hat, augmented, misspecified e                 +0.290456     err vs estimand  −0.030249
    ATO indexed by the MISSPECIFIED score              +0.295884     err vs estimand  −0.024821
    tau_hat − ATO(misspecified score)                                                 −0.005428

    tau_hat, augmented, correctly specified e          +0.316556     err vs estimand  −0.004149
```

The estimator misses the estimand by −0.0302, a 9.4% relative error, and sits within 0.0054 of the
population the wrong model defines — against −0.0041 when `e` is correct. That is [§8]'s paragraph as a
number.

**And it is worth naming why the two errors are so close in size.** With an *oracle* `m_a`, the two arm
terms of §8.1 have mean zero and the correction term is `Σh(m₁−m₀)/Σh` evaluated at the misspecified
`h` — which *is* `ATO(e_mis)` exactly. So the estimator is not approximately recovering the wrong
population, it is recovering it by construction, and the residual −0.0054 is sampling error in the arm
terms. That is a sharper statement than "does not recover the true ATO" and it is what §15.10 asserts
positively.

**The bias tracks the variation in `tau(X)`, and the companion has to be constant on the
risk-difference scale.** The roadmap's reason for insisting the assertion be positive is that *"a test
built on a homogeneous effect will appear to show that it is [doubly robust], because every weighted
average treatment effect coincides in that case"*. That is true when `tau(X)` is constant — and a
logistic outcome model with no interaction term is **not** that, because the link's curvature makes a
constant log-odds effect a non-constant risk difference. Measured, R = 200 replicates, misspecified `e`
and oracle `m_a` throughout:

```
  construction                          sd(tau(X))   n         mean err        se     |t|
  heterogeneous, mod = 1.5                  0.2435   20,000    −0.025338  0.000453    55.9
  no-interaction logit                      0.0523   20,000    −0.009901  0.000465    21.3
  constant risk difference, DELTA 0.12      0.0000   20,000    +0.000463  0.000536     0.9
```

The middle row is the one a reasonable draft would have written as the companion, and it is
**detectably biased** — `|t| = 21.3`. Only the constructed constant-risk-difference arm is
indistinguishable from zero. §15.10's companion is built as `m₀` from a logistic model and
`m₁ = m₀ + DELTA` with `DELTA = 0.12`, and `m₀` bounded into `(0.25, 0.70)` so that `m₁` stays a
probability: without the bound the difference is truncated at the top and the construction is no longer
constant where it matters most. `DELTA`'s value is `0.12` and is written in `§15.0.5`'s fence — an
earlier draft named the constant and never gave it (§22.3 item 1).

**And the test cannot be run at the cohort's own size.** At n = 92 the per-replicate spread swamps the
quantity being asserted: measured over 400 replicates, mean error −0.025665 with **sd 0.087669**, range
−0.291374 to +0.228369, and **40.5% of replicates have the wrong sign**. A test written on a 92-row
fixture measures noise. The `n` is therefore chosen from a measured separation, and this is the only
tuning constant §15.10 has. Measured on `dr_population` (§15.0.5), twelve seeds per cell:

```
  n          heterogeneous: min |err|   constant: max |err|   separation   cost / 24 draws
  200,000                    0.022950              0.004929   4.66x            0.31 s
```

**§15.10 uses n = 200,000** with the seed pinned at `20260825` and the twelve seeds written out, where
the heterogeneous arm's errors span **[0.022950, 0.030249]** and the companion's span
**[0.000176, 0.004929]**. The band between **0.004929 and 0.022950 is empty across all twelve seeds**,
and the two asserted bounds — companion below `0.010`, heterogeneous above `0.015` — both sit inside it,
verified. This is Stage 8 §6.3's empty-band calibration, applied to a tolerance instead of to a
coefficient bound. Cost is **13 ms per draw**, so the test is two draws and about 0.03 s.

**One claim an earlier draft made about smaller `n` is withdrawn, because it was arithmetically false.**
It carried a four-row table and concluded *"At n = 20,000 the two arms overlap — 0.0164 against 0.0134
— so the obvious 'large enough' choice is not large enough."* Those two numbers are
`heterogeneous: min |err|` and `constant: max |err|`, and `0.0134 < 0.0164`: the arms did **not** overlap,
the band was empty, and the table's own `ratio` column said `1.22x`. The defensible claim was that a
1.22× separation is too narrow to survive a numpy version bump, which is a different argument. Round 2
re-derived §8.5 and reported finding no error in it; this is the error it did not find. §22.3 item 3.
The row above is re-measured on the specified construction at the `n` the test actually uses, and no
smaller `n` is now claimed about.

---

## 9. The rare-minority guard, and why it is decided once

### 9.1 The minority cell, and the population it is counted on

`min(events, non-events)`, and **not** the event count. `config.RARE_MINORITY_THRESHOLD`'s own comment
says why (`config.py:273-276`): *"A constant named `RARE_EVENT_THRESHOLD` would invite the exact
misreading the amendment corrects: TICI 2b-3 has 114 events and 6 non-events, so it passes an
event-count rule while failing that rule's rationale."*

**The population it is counted on is the outcome's own [§11] population, and that has to be said
because the two available answers differ.** `config.py`'s comment counts 114 events for TICI; §4.2
measures 85 on `in_model & notna()`. Both are correct about different frames — the comment is about the
workbook, §4.2 is about the ATO population — and the non-event count is 6 in both, which is exactly the
coincidence that would let a wrong choice go unnoticed. The rule is applied on the [§11] population
because that is **the population the nuisance model is fitted on**, and the constraint the rule encodes
is that a model cannot carry more parameters than its own fitting sample supports. §15.9 asserts the
count against 6 *and* asserts that the population used is the masked one, on the frame where the two
differ.

### 9.2 The rule

```python
def _augmentable(outcome: str, minority: int) -> bool:
    """Whether `outcome` is augmented [§8, as amended 2026-08-07].

    Two ways in, and the second is not a special case of the first:

      * the minority cell supports the shared [§6] covariate list, or
      * the outcome DECLARES a reduced m_a(X) in C.OUTCOME_MODEL_OVERRIDES.

    The amendment's item 2 is explicit that a declared reduction is augmented "rather than dropped
    from augmentation", so an override is a route IN and not a modifier of an outcome that qualified
    anyway. The roadmap said such outcomes "are not augmented"; that sentence was stale and §22
    item 1 amends it (§7.2).

    Reading the registry rather than a second list of names is what makes this unable to drift: an
    outcome cannot become reduced except by being named in config.py, and cannot be named there
    without becoming augmentable here.
    """
    return (minority >= C.RARE_MINORITY_THRESHOLD
            or outcome in C.OUTCOME_MODEL_OVERRIDES)
```

The declared path is recorded as a string on the estimate, not inferred from which fields are
populated:

```python
def _augmented_path(outcome: str, minority: int) -> str:
    """"full" | "reduced" | "unaugmented" -- the three declared paths, named once."""
    if not _augmentable(outcome, minority):
        return "unaugmented"
    return "reduced" if outcome in C.OUTCOME_MODEL_OVERRIDES else "full"
```

### 9.3 The seven, measured, and the SAP's own prediction confirmed

```
  outcome         minority   threshold   override   path          m_a(X) covariates
  mrs_0_2_90d           42     >= 10        no      full          the [§6] nine
  mrs_0_1_90d           24     >= 10        no      full          the [§6] nine
  tici_2b_3              6     BELOW        YES     reduced       center, atrial_fib
  sich                   5     BELOW        no      unaugmented   --
  ph2                    9     BELOW        no      unaugmented   --
  death_90d             25     >= 10        no      full          the [§6] nine
  mrs_5_6_90d           28     >= 10        no      full          the [§6] nine

  five augmented (four full, one reduced), two unaugmented
```

The amendment's closing sentence is *"No other outcome is reduced; symptomatic ICH and parenchymal
haematoma type 2 have a minority cell below the threshold and remain unaugmented"*. Measured: `sich` 5
and `ph2` 9, both below 10, both unaugmented, and no other outcome overridden. **The SAP's prediction
about this workbook is confirmed rather than assumed**, which is worth one line because it is the only
prespecified numerical claim about the outcome distribution that the protocol makes.

Two of the four **safety** outcomes are augmented with the full [§6] list. [§8] says the rare-outcome
rule *"will normally apply to the safety family, which §10 already treats as descriptive"* — "normally"
carries it, and `death_90d` at 25 and `mrs_5_6_90d` at 28 are not close to the threshold. This is
recorded because a reader who takes "the safety family is unaugmented" from that sentence will find two
augmented safety estimates in the output.

### 9.4 `ph2` sits one below the threshold, and that is not a curiosity

`ph2`'s minority cell is **9** against a threshold of **10**. One additional parenchymal haematoma in
the ATO population flips it from unaugmented to fully augmented — a different estimator, not a
different number from the same one. Three things follow and only the third is actionable here:

- The threshold is prespecified at Stage 1 and is **not** revisited. It was fixed before any outcome
  was examined by arm, the amendment says so explicitly, and moving it now because a measured count
  landed at 9 is exactly the kind of choice the prespecification exists to prevent. §17 records that it
  was looked at and left alone.
- The margin is disclosed. §10.2's table prints the minority cell beside the path for every outcome, so
  a reader can see that one of the seven is one event from a different estimator without reading this
  document.
- **It makes §9.5 necessary rather than tidy**, which is the next section.

### 9.5 Eligibility is decided once, on the point estimate, and carried into every replicate

A minority cell is a property of the **sample**. [§10] resamples. So a replicate can carry `ph2` over
the threshold, and if `_augmentable` were re-evaluated per replicate, some replicates would contribute
an augmented estimate and others an unaugmented one to the same percentile interval. That is a mixture
of two estimators, which is what [§10] forbids in terms: *"Replicates whose prespecified fit fails are
dropped and counted; **never substituted with a different estimator**"*
(`implementation_roadmap.md:496`).

**How often it would happen, measured** over 400 stratified replicates with the propensity model refit
in each:

```
  outcome         point minority   point path      replicates crossing the threshold
  mrs_0_2_90d               42     full                                       0.0%
  mrs_0_1_90d               24     full                                       0.0%
  tici_2b_3                  6     reduced                                    0.0%
  sich                       5     unaugmented                                3.8%
  ph2                        9     unaugmented                              46.1%
  death_90d                 25     full                                       0.0%
  mrs_5_6_90d               28     full                                       0.0%
```

**`ph2` crosses in 46.1% of replicates.** A per-replicate rule would make its interval a quantile over a
near-even mixture of two estimators — the worst available case of the thing [§10] prohibits, and not a
tail event to be waved at. This is the measurement that turns §9.5 from a tidiness argument into a
requirement.

Two details of the rule, both load-bearing:

- **`tici_2b_3` crosses in 0.0%, by construction.** Its override makes `_augmentable` return `True`
  whatever the minority cell does, so the one outcome whose cell is smallest is also the one immune to
  resampling variation in it. That is a property of the registry design worth noticing rather than a
  coincidence: naming an outcome fixes its estimator.
- **`sich` at 3.8% is small and is not zero.** A rule that only looked at `ph2`'s number would conclude
  the problem is one outcome's; it is two.

So `secondary` computes the path on the frame it is given, and **Stage 10 passes the point estimate's
paths in** rather than letting each replicate decide. The mechanism is Stage 10's to build and §14
states the contract; what this stage owes is that `augmented_path` is on the returned object, per
outcome, so there is something to pass.

**What happens to a replicate whose declared model then cannot be fitted.** It raises `FitError` from
`model.firth` and [§10] drops and counts it. It does **not** fall back to the unaugmented form — that
would be the substitution the whole section exists to prevent, arrived at by a different route.
Measured, the rate is small: `model.firth` raised on **0.3%** of replicates for `sich` and **0%** for
every other outcome, over the same 400.

### 9.6 There is no bound on the nuisance coefficient, and Stage 8's reason for having one does not transfer

Stage 8 §6 established that separation in `polr` converges silently and therefore needs
`POLR_MAX_ABS_BETA` to become a countable failure. Both halves of that finding were re-asked of
`model.firth`.

**The first half reproduces exactly.** A perfectly separated 40-record frame handed to `model.firth`
**converges in 6 iterations**, rescales 0, halvings 0, and returns `exp(β_treatment) = 441.005397` with
`max|β| = 6.089057` and fitted probabilities strictly interior in `[4.545e-02, 0.954546]`. The
unpenalised MLE on the same frame does not merely grow — **it has no maximum to reach, and every route
says so differently**. The Newton route raises `LinAlgError: Singular matrix` **once it is allowed to
reach the optimum, at `maxiter >= 200`**; at statsmodels' own default of 35 it stops without raising, at
`max|β| = 67.87` with `converged = False`, and at 100 it returns 67.45. A BFGS route capped at 60
reaches `max|β| = 33.0927` and **reports success**. The cap is therefore part of the measurement and an
earlier draft of this paragraph omitted it (§22.4 item 3) — and stating it makes the point stronger
rather than weaker: a raise reads as a numerical accident, where a coefficient that changes with the
iteration cap is the definition of an optimum that is not there. Firth's penalty turns a fit that has no maximum into one that
converges in six iterations with **every safeguard counter at zero**, and **hides the separation just as
completely** — nothing in the returned `Fit` is out of range.

The frame is `separated_frame` (§15.0.4a) and is written out as code. An earlier draft pinned this to
`exp(β) = 765.443` and an unpenalised `max|β| = 73.12` against a frame that appeared nowhere in the
document, and DoD-9 gated completion on observing the first of those two numbers (§22.3 item 1).

**The second half does not.** Stage 8 could choose 14.0 because the fits were bimodal with nothing
between `|β| = 8.79` and `18.81` (later re-measured to 11.04, Stage 8 §21.1 item 24). Measured here over
400 replicates:

```
  outcome           n   raised      min   median      p95      p99       max    >8   >14
  mrs_0_2_90d     399        0   0.9052   3.8312   7.5160   9.1537   14.2152    15     1
  mrs_0_1_90d     399        0   0.9776   2.6035   5.9665  10.3211   37.4910     9     4
  tici_2b_3       399        0   1.2696   2.7202   5.2411   6.1702    8.4047     1     0
  sich            398        1   0.8558   6.9554  19.7819  36.4966   80.9825   170    60
  ph2             399        0   1.1267   3.9003  11.3722  19.7882   31.6525    47    11
  death_90d       399        0   1.1743   4.6312   9.0089  13.3151   52.7679    35     3
  mrs_5_6_90d     399        0   1.4498   5.1203   9.8685  13.3317   15.5013    53     4
```

`sich` runs **continuously from 0.8558 to 80.9825**, with 42.7% of replicates above 8 and 15.1% above 14
and **no gap anywhere**. There is no empty band, so there is no non-arbitrary bound: any threshold drops
legitimate fits, and one placed where Stage 8's is would drop 15% of `sich`'s replicates for a reason
that is not about the data.

**Why no bound is needed anyway, and the reason is structural rather than empirical.** At Stage 8, `β`
*is* the estimand, so a degenerate `β` is a degenerate answer — `exp(21)` in a percentile interval. Here
`m_a(X)` is a **nuisance** and only its *predictions* enter `tau`, and predictions are probabilities,
bounded in `[0, 1]` whatever `β` does. Every term of `augmented_rd` is therefore bounded, and there is
no `exp(β)`-style tail to protect against. Measured on synthetic frames of 92 records at three event
rates, the augmented estimate's error stayed inside `[−0.344, +0.298]` across 1778 fits.

**What is lost instead, and what it costs.** A separated `m_a` is an over-fitted one: its fitted values
reproduce `Y`, so the correction term removes signal rather than residual confounding. That is a
variance-and-bias cost, not an infinity. On the synthetic frames — where the truth is a **constant**
risk difference, so the outcome model has nothing legitimate to exploit — augmentation raised RMSE
against the unaugmented estimator in all three arms:

```
  events   augmented RMSE   unaugmented RMSE   ratio
      ~5         0.056678           0.047799   1.19x
      ~9         0.071635           0.064810   1.11x
     ~25         0.104487           0.098285   1.06x
```

This is [§8]'s own caveat as a number — *"at this treated-arm size flexible nuisance models add variance
rather than robustness"* — and it is a **lower bound on the cost, not an estimate of it**, because the
construction is the case where the augmentation can only hurt.

**And the honest limitation of that measurement.** The synthetic construction did **not** reproduce the
condition it was built to study: its `max|β|` reached only 9.458, 11.475 and 3.466 against the
workbook's 80.98, because independent standard-normal covariates carry none of the centre structure that
drives the workbook's near-separation — treatment is nearly determined by centre [Stage 7 §6.4] and USZ
is absent entirely. So the RMSE table above is measured on frames unlike the ones the tail comes from,
and it does **not** establish that `sich`'s heavy `max|β|` tail leaves its interval usable. §16 files
that, `TODOS.md` carries the trigger, and §14 hands Stage 10 the diagnostic that would answer it: **the
distribution of `max|β|` across replicates, per outcome, reported beside the failure counters** — the
same move Stage 8 §11 made with the distribution of `len(fit.alpha)`, and for the same reason.

---

## 10. The three `model` audit entries

### 10.1 No new kind, and that is `data.py`'s own declaration honoured

`data.KINDS` stays nine. Stage 7 added two entries under the existing `model` kind, Stage 8 added three,
Stage 9 adds three, and none of the three stages touched `data.py` — which `data.py:151-152` predicted.
The ledger:

```
  load 7 / derive 4 / classify 1 / build 6 / fit 4 / assess 2 / primary 3 / secondary 3  =  30
```

Measured: 22 entries after `propensity.fit`, which reconciles with Stage 8 §13's 24-after-`assess` and
27-after-`primary`. §15.14 asserts 30 and asserts `len(KINDS) == 9`.

### 10.2 Entry 1 — `binary_estimable`: the seven populations and the path each takes

One entry, seven rows, recorded **before** any estimate is computed, so that a failure in the estimator
leaves a log that still names the populations. This is `propensity._record_exclusion`'s ordering rule
(`propensity.py:441-443`) applied one stage on: *"the exclusion is a fact about the frame, established
before any modelling choice, and a log that names it is useful precisely when what follows fails."*

```
  | outcome      | family    | n [§11] | lost | minority | path        | m_a(X)             |
  |--------------|-----------|---------|------|----------|-------------|--------------------|
  | mrs_0_2_90d  | secondary |      92 |    0 |       42 | full        | [§6] set, 9        |
  | ...          |           |         |      |          |             |                    |
  | tici_2b_3    | secondary |      91 |    1 |        6 | reduced     | center, atrial_fib |
  | sich         | safety    |      92 |    0 |        5 | unaugmented | --                 |
```

`n` on the entry is the **cohort's** 93, and the per-outcome denominators are in the table — because
`Audit.record`'s `n` is one number and there are seven here. The `detail` string states the three
populations in the Stage 8 §4.1 form: the cohort's 93, the ATO's 92, and then per outcome.

**So `_record_estimable` takes `ps`, and an earlier draft's signature could not build this table.** It
was `_record_estimable(df, populations, minorities, audit)`, and both the `lost` column and the detail
string need `ps.in_model`: `lost` is `|in_model| − |in_model & notna|`, which is **not**
`df[key].isna().sum()` — that counts over the cohort, not over the ATO population, and the two differ
whenever a covariate-incomplete record is also outcome-missing. Nor is `in_model` recoverable from
`populations`, which holds only the intersections, and it is never a column on the frame
(`propensity.py:84`: *"A return value and never a column on the frame"*). The signature is
`_record_estimable(df, ps, populations, minorities, audit)`. §22.3 item 10.

**No `case_ids`.** `_MUST_NAME_CASES` does not cover this kind and this entry removes no patient from
anything — the seven masks are subsets of `in_model`, whose exclusion `propensity.py` already named. An
entry naming the same excluded patient a second time is not a second fact. §15.14 asserts the tuple is
empty, and asserts that the one record `tici_2b_3` loses is **not** named — for Stage 8 §9.3's reason,
that what an entry deliberately omits needs asserting as much as what it carries.

### 10.3 Entry 2 — `binary_estimates`: the three estimates per outcome, side by side

```
  | outcome      | n  | p1_w | p0_w | RD_w | OR_w | corrected | tau (model-assisted) | m_a(X) |
  |--------------|----|------|------|------|------|-----------|----------------------|--------|
```

**The last column is `m_a(X)` and not `reduced`, which is what an earlier draft called it.** `reduced` is
the name of a *boolean field* on `BinaryEstimate` (§3.1), and this column renders a *tuple of covariate
names*. A column and a field sharing a name while carrying different types is the kind of collision that
survives review and then costs an afternoon; `m_a(X)` is also what §10.2 already calls the same content,
so the three tables now agree. §22.3 item 4.

Seven things about this table, each of which is a reporting requirement rather than a layout choice:

- **The unaugmented risk difference and `tau` are adjacent columns.** [§8]: *"the augmented estimate is
  therefore reported as model-assisted, always alongside the unaugmented weighted risk difference;
  material disagreement between the two is evidence about the outcome model."* Adjacency is what makes
  the comparison the protocol requires available without arithmetic.
- **The two raw weighted proportions are printed**, uncorrected, always. §6.3 returns them for this.
- **`corrected` is printed even when it is false for all seven**, which on v7 it is (§6.4). A column that
  disappears when nothing triggers it is a column a reader cannot tell was checked.
- **`m_a(X)` names the covariate list**, not a boolean, wherever it is non-empty. The amendment:
  *"Every reduced specification is reported beside its estimate."* A `True` is not the specification.
- **No cell in any of the three tables may contain a `|`.** `data._md_table` (`data.py:242-245`) does no
  escaping and sizes its separator from `len(rows[0])`, so a literal pipe inside any cell yields a header
  with more markdown cells than its separator and the table renders crooked. Stage 7 shipped this twice
  (`|SMD| < 0.1` and `worst |SMD|`, rendering 8/6 and 15/13) and §10.4 was about to ship it a third time.
  §15.14 asserts the structural property per entry, on the rendered text — **not** via the hash, because
  DoD-4's byte-identity check passes on a table that is identically broken twice.
- **`tau` is `--` and never `0.0` for an unaugmented outcome.** Stage 6 §9's rule: the deliberateness of
  an absence lives in a mask, never in a value. `_fmt` (`data.py`) is what renders `None`.
- **The header says "model-assisted".** Not "augmented (doubly robust)", not "AIPW". §15.13 scans the
  rendered text.

### 10.4 Entry 3 — `binary_outcome_models`: one row per augmented outcome

Five rows on the workbook, not seven. The two unaugmented outcomes have no model and a row asserting
that would be a row about nothing.

```
  | outcome     | covariates         | k | rows | iters | max_abs_beta | rescales | halvings | dropped    |
  |-------------|--------------------|---|------|-------|--------------|----------|----------|------------|
  | mrs_0_2_90d | [§6] set           | 12|   92 |     6 |       3.0933 |        0 |        0 | center_USZ |
  | tici_2b_3   | center, atrial_fib | 4 |   91 |     7 |       2.2840 |        0 |        0 | center_USZ |
```

**The column is `max_abs_beta` and not `max|beta|`.** Two literal pipes in a header cell break
`data._md_table` (`data.py:242-245`): it does no escaping and builds its separator from
`len(rows[0])`, so `| ... | max|beta| | ... |` presents **eleven** markdown cells against a separator of
**nine** and the whole table renders misaligned. This is Stage 7's shipped defect — its learnings entry
records `|SMD| < 0.1` and `worst |SMD|` rendering 8/6 and 15/13 — and an earlier draft of this section
reintroduced it. It is invisible to DoD-4, because a table that is identically broken under both hash
seeds is still byte-identical. §15.14 now asserts the structural property directly, and the ASCII name is
what makes the assertion pass rather than what makes it unnecessary. §22.3 item 7.

`max_abs_beta` is in the table for §9.6's reason: no bound rejects a fit, so the only thing standing
between a separated nuisance model and an unremarked estimate is a number in the log that a reader can
see. The
three safeguard counters are there for Stage 6 §3.2's reason — *"a safeguard whose activation nobody can
count is a safeguard nobody can evaluate"* — and `dropped` is there because `center_USZ` disappearing
from every design is a fact about the cohort that should not have to be rediscovered.

### 10.5 The seven audit privates, written out

Three tables, one detail string, three recorders. `_fmt` (`data.py:222-231`) is the one float formatter
and renders `None` as `missing`; §10.3's *"`tau` renders as `--`"* is that rule with the table supplying
the dash for a field that is *structurally* absent rather than missing — an unaugmented outcome has no
`tau`, which is different from a `tau` nobody could compute.

```python
def _estimable_table(df: pd.DataFrame, ps: propensity.Propensity,
                     populations: dict[str, pd.Series],
                     minorities: dict[str, int]) -> tuple[tuple[str, ...], ...]:
    """§10.2's seven rows. Takes `ps` because `lost` is not derivable without in_model (§10.2).

    `lost` is |in_model| - |in_model & notna|, and it is NOT df[key].isna().sum(): that counts over
    the cohort rather than over the ATO population, and the two differ whenever a
    covariate-incomplete record is also outcome-missing. On v7 they happen to agree for six of the
    seven, which is exactly the coincidence that would let the wrong expression pass unnoticed.

    Ranges over C.BINARY_OUTCOMES so the row order is the registry's, which is what makes the
    rendered log byte-identical across hash seeds (§15.14).
    """
    in_model = int(ps.in_model.sum())
    header = ("outcome", "family", "n [§11]", "lost", "minority", "path", "m_a(X)")
    rows: list[tuple[str, ...]] = [header]
    for key in C.BINARY_OUTCOMES:
        n = int(populations[key].sum())
        path = _augmented_path(key, minorities[key])
        covariates = ("--" if path == "unaugmented"
                      else ", ".join(C.outcome_model_covariates(key)))
        rows.append((key, C.OUTCOMES[key].family, str(n), str(in_model - n),
                     str(minorities[key]), path, covariates))
    return tuple(rows)


def _estimable_detail(df: pd.DataFrame, ps: propensity.Propensity,
                      populations: dict[str, pd.Series]) -> str:
    """The three populations in the Stage 8 §4.1 form: the cohort's, the ATO's, then per outcome.

    Stated as a RANGE over the seven rather than seven numbers, because the table already carries
    them per outcome and a detail string that repeats a table is a second place to disagree with it.
    """
    sizes = sorted({int(populations[k].sum()) for k in C.BINARY_OUTCOMES})
    span = str(sizes[0]) if len(sizes) == 1 else f"{sizes[0]}-{sizes[-1]}"
    return (f"{len(df)} cohort record(s); {int(ps.in_model.sum())} in the [§7] overlap population; "
            f"{span} on the seven [§8] outcome populations, one per outcome in the table. Each "
            f"estimate uses its OWN denominator [§11] and no outcome is estimated on another's.")


def _estimates_table(estimates: dict[str, BinaryEstimate]) -> tuple[tuple[str, ...], ...]:
    """§10.3's seven rows: both proportions, RD, OR, the correction flag, tau, the m_a(X) list.

    `corrected` is rendered for every row even when it is False for all seven, which on v7 it is: a
    column that disappears when nothing triggers it is a column a reader cannot tell was checked
    (§10.3). `tau` renders `--` and never 0.0 for an unaugmented outcome, because the deliberateness
    of an absence lives in a mask and never in a value [Stage 6 §9].

    NO CELL MAY CONTAIN A PIPE and none does: `data._md_table` does no escaping (§10.3, §15.14).
    """
    header = ("outcome", "n", "p1_w", "p0_w", "RD_w", "OR_w", "corrected",
              "tau (model-assisted)", "m_a(X)")
    rows: list[tuple[str, ...]] = [header]
    for key in C.BINARY_OUTCOMES:
        est = estimates[key]
        rows.append((
            key, str(int(est.in_estimate.sum())),
            _fmt(est.proportion[_TREATED]), _fmt(est.proportion[_COMPARATOR]),
            _fmt(est.rd), _fmt(est.odds_ratio), str(est.or_corrected),
            "--" if est.augmented is None else _fmt(est.augmented),
            "--" if est.covariates is None else ", ".join(est.covariates)))
    return tuple(rows)


def _models_table(estimates: dict[str, BinaryEstimate]) -> tuple[tuple[str, ...], ...]:
    """§10.4's one row per AUGMENTED outcome -- five on v7, not seven.

    The two unaugmented outcomes have no model and a row asserting that would be a row about
    nothing. `max_abs_beta` is spelled in ASCII: `max|beta|` would put two pipes in a header cell and
    break the rendered table in a way byte-identity does not catch (§10.4, §15.14).

    The three safeguard counters are here for Stage 6 §3.2's reason -- "a safeguard whose activation
    nobody can count is a safeguard nobody can evaluate" -- and `dropped` because `center_USZ`
    vanishing from every design is a fact about the cohort nobody should rediscover.
    """
    header = ("outcome", "covariates", "k", "rows", "iters", "max_abs_beta",
              "rescales", "halvings", "dropped")
    rows: list[tuple[str, ...]] = [header]
    for key in C.BINARY_OUTCOMES:
        est = estimates[key]
        if est.fit is None:
            continue
        rows.append((
            key, ", ".join(est.covariates or ()), str(len(est.fit.columns)),
            str(int(est.in_estimate.sum())), str(est.fit.iterations),
            _fmt(float(np.max(np.abs(est.fit.beta)))),
            str(est.fit.rescales), str(est.fit.halvings),
            ", ".join(est.dropped) if est.dropped else "--"))
    return tuple(rows)


def _record_estimable(df: pd.DataFrame, ps: propensity.Propensity,
                      populations: dict[str, pd.Series], minorities: dict[str, int],
                      audit: Audit) -> None:
    """§10.2's entry. `n` is the COHORT's, because Audit.record takes one number and there are seven.

    No `case_ids`: `_MUST_NAME_CASES` does not cover the `model` kind (data.py:196, 215) and this
    entry removes no patient from anything -- the seven masks are subsets of in_model, whose
    exclusion `propensity._record_exclusion` already named. An entry naming the same excluded patient
    a second time is not a second fact (§10.2).
    """
    audit.record("model", "binary_estimable", len(df),
                 _estimable_detail(df, ps, populations),
                 table=_estimable_table(df, ps, populations, minorities))


def _record_estimates(estimates: dict[str, BinaryEstimate], audit: Audit) -> None:
    """§10.3's entry. `n` is the number of [§8] binary estimates, which is seven and not a population.

    The detail names what the table cannot: that the unaugmented risk difference and tau sit in
    ADJACENT columns because [§8] requires the comparison, and that "model-assisted" is the licensed
    description (§8.3). §15.13 scans this string for the forbidden phrase.
    """
    augmented = sum(1 for e in estimates.values() if e.augmented is not None)
    corrected = sum(1 for e in estimates.values() if e.or_corrected)
    audit.record(
        "model", "binary_estimates", len(estimates),
        f"{len(estimates)} [§8] binary outcome(s), each on its own [§11] denominator. "
        f"{augmented} carry a model-assisted augmented risk difference beside the unaugmented one, "
        "in adjacent columns, because [§8] requires the comparison and calls material disagreement "
        "between them evidence about the outcome model rather than confirmation of either. "
        f"§6.3's continuity correction fired on {corrected} of them. No interval and no p-value: "
        "[§10] owns those.",
        table=_estimates_table(estimates))


def _record_models(estimates: dict[str, BinaryEstimate], audit: Audit) -> None:
    """§10.4's entry. `n` is the number of nuisance models FITTED -- five on v7, not seven.

    No coefficient vector, for §10.4's reason: `m_a(X)` is a nuisance and [§8]'s first line forbids
    reporting adjusted effects, so printing its treatment coefficient invites the misreading that a
    conditional odds ratio is the estimate. No standard error and no interval anywhere.
    """
    fitted = [e for e in estimates.values() if e.fit is not None]
    audit.record(
        "model", "binary_outcome_models", len(fitted),
        f"{len(fitted)} Firth nuisance model(s) for the augmented outcomes; "
        f"{len(estimates) - len(fitted)} outcome(s) are unaugmented under §9.2 and have none. "
        "max_abs_beta is reported because no bound rejects a fit [§9.6], so it is the only thing "
        "between a separated nuisance model and an unremarked estimate. No coefficient vector, no "
        "standard error, no interval: only the predictions of these models enter any estimate.",
        table=_models_table(estimates))
```

`_fmt` is imported from `data` alongside `Audit`, which is the import `outcome.py` already has — §3's
claim that its imports do not change survives these seven functions.

**What the entries deliberately do not carry.** No coefficient vector — `m_a(X)` is a nuisance and
[§8]'s first line forbids reporting adjusted effects, so printing its treatment coefficient invites
exactly the misreading that a conditional odds ratio is the estimate. No standard error and no interval,
anywhere: [§10]'s percentile bootstrap is the prespecified one and there is no second (Stage 8 §5.6,
identically). No p-value.

---

## 11. `secondary`, written out

```python
def secondary(df: pd.DataFrame, ps: propensity.Propensity, audit: Audit,
              paths: dict[str, str] | None = None) -> Secondary:
    """The [§8] binary estimates: seven outcomes, each on its own [§11] denominator.

    `paths` is Stage 10's frozen point-estimate decision (§9.5, §11.1). None means "decide from this
    frame", which is correct on the point estimate and WRONG in a replicate: `ph2` crosses the
    threshold in 46.1% of replicates and [§10] forbids mixing estimators across one interval. An
    earlier draft omitted the parameter entirely, which made §12.1's and §14's "the paths must be
    passed IN" unsatisfiable through this function (§22.3 item 11).

    Six things about the order below, each of which is a failure if moved:

      * `_assert_readable` runs FIRST, raising before any value is read, so a frame missing an
        outcome column reports S1 rather than a bare pandas KeyError from S2 (§4.4a).
      * `_assert_secondary_inputs` runs SECOND and over the whole frame, before any outcome is
        touched, so a frame that fails S2-S5 reports that rather than a message about one outcome.
      * `_record_estimable` runs BEFORE the estimator loop, so a raise inside an estimate still
        leaves a log naming the seven populations (§10.2).
      * `in_estimate` is built per outcome and `e`, `w` are sliced to IT, never to `in_model`. This is
        the one place either is read, which is how Stage 6 §9's "range over the mask, never over
        notna()" is honoured -- as a fact about two lines rather than as a warning (§12).
      * The minority cell is counted on `in_estimate` and not on the cohort (§9.1).
      * `_augmentable` is consulted BEFORE `outcome_model` is called, so an unaugmented outcome never
        fits a nuisance model. Reversing them would make `sich`'s 17-iteration fit (§7.4) happen on
        every call and be discarded, and would put a FitError on a path [§8] says has no model.

    NO FitError IS CAUGHT HERE, and that is the whole of §9.5's no-substitution rule as code. A
    `model.firth` failure on a declared nuisance model propagates, [§10] drops and counts the
    replicate, and the estimate is ABSENT rather than silently downgraded to the unaugmented form.
    A try/except around `outcome_model` is the most natural thing an implementer writes when a fit
    fails on 0.3% of `sich` replicates and it is the thing that turns this stage into a two-estimator
    mixture. §15.7a asserts it and shows the caught version passing every other test.

    Deterministic, and writes no file. No seed is held, no clock is read, and the only mutable object
    touched is the `Audit` passed in. Stage 10 must pass a throwaway one per replicate, for Stage 8
    §11's reason and more so: three entries per replicate at N_BOOT = 2000 is 6000 entries, and one
    of them carries a seven-row rendered table.
    """
    _assert_readable(df, ps)
    _assert_secondary_inputs(df, ps)
    if paths is not None and set(paths) != set(C.BINARY_OUTCOMES):
        raise C.SchemaError(
            f"secondary was given frozen paths for {sorted(paths)} against "
            f"{sorted(C.BINARY_OUTCOMES)}. A partial map would let some outcomes carry the "
            "point estimate's path and others re-decide from the replicate, which is the "
            "two-estimator mixture §9.5 forbids arrived at one outcome at a time.")

    populations: dict[str, pd.Series] = {}
    minorities: dict[str, int] = {}
    for key in C.BINARY_OUTCOMES:
        populations[key] = ps.in_model & df[key].notna()
        y = df.loc[populations[key], key].to_numpy(dtype=float)
        minorities[key] = min(int((y == 1.0).sum()), int((y == 0.0).sum()))
    _record_estimable(df, ps, populations, minorities, audit)

    estimates: dict[str, BinaryEstimate] = {}
    for key in C.BINARY_OUTCOMES:
        in_estimate = populations[key]
        sub = df.loc[in_estimate]
        y = sub[key].to_numpy(dtype=float)
        a = sub[C.TREATMENT].to_numpy(dtype=float)
        w = ps.w.loc[in_estimate].to_numpy(dtype=float)
        e = ps.e.loc[in_estimate].to_numpy(dtype=float)
        _assert_estimable(key, y, a, w)

        odds_ratio, corrected, share = marginal_odds_ratio(y, a, w)
        path = paths[key] if paths is not None else _augmented_path(key, minorities[key])
        tau: float | None = None
        fit: model.Fit | None = None
        covariates: tuple[str, ...] | None = None
        dropped: tuple[str, ...] = ()
        if path != "unaugmented":
            covariates = C.outcome_model_covariates(key)
            fit, X, dropped = outcome_model(df, in_estimate, key)
            m1, m0 = _counterfactuals(fit, X)
            tau = augmented_rd(y, a, w, _tilt(e), m1, m0)

        estimates[key] = BinaryEstimate(
            outcome=key, family=C.OUTCOMES[key].family, minority=minorities[key],
            rd=weighted_rd(y, a, w), odds_ratio=odds_ratio, or_corrected=corrected,
            proportion=dict(share), augmented_path=path, augmented=tau, covariates=covariates,
            reduced=key in C.OUTCOME_MODEL_OVERRIDES, dropped=dropped,
            in_estimate=in_estimate, fit=fit)

    _record_estimates(estimates, audit)
    _record_models(estimates, audit)
    return Secondary(estimates=estimates)
```

`proportion=dict(share)` and not `proportion=share`: `BinaryEstimate` is `frozen=True`, but a `dict`
field is mutable regardless and `share` is the very object `weighted_proportion` returned, so
`estimate.proportion[1] = 99.0` succeeds and changes a frozen estimate — and, before the copy, changed
the object `marginal_odds_ratio` had also returned. The copy is one word and §15.14 asserts the estimate
is unchanged after a caller mutates what it was handed. §22.3 item 14.

`weighted_rd(y, a, w)` recomputes the two proportions `marginal_odds_ratio` already returned in `share`.
That is deliberate and it is not an oversight: it is two `np.average` calls on at most 92 rows against
`model.design`'s 17.40 ms (§2), and routing the reported risk difference through the **public** estimator
is what makes §15.3's hand computation a test of the number this function returns rather than of a
sibling expression. §22.3 item 15 records that it was considered.

### 11.1 The frozen-path parameter, and why it belongs here rather than in Stage 10

§9.5 closes with *"The mechanism is Stage 10's to build and §14 states the contract; what this stage owes
is that `augmented_path` is on the returned object, so there is something to pass."* Round 3 disagrees
with the first half of that sentence, and the reason is that there was nothing to pass it **to**.

The contract §12.1 and §14 both state — paths decided once, carried into every replicate, never
recomputed — is a contract about a *call*. A function whose signature cannot express it leaves Stage 10
two options: call `secondary` and recompute the paths, which is the 46.1% mixture; or reimplement the
seven-outcome loop, which §14 itself warns is *"how item 2 gets violated by accident"*. §14 then made the
contradiction explicit by saying a replicate loop *"can call `secondary` whole"*.

Three properties of the parameter, each chosen against the alternative:

- **Optional, defaulting to `None`.** The point-estimate call is the one that *decides* the paths, so it
  must be able to run without them. A required parameter would make the first call impossible.
- **Complete or absent, never partial.** The guard rejects a map that does not cover
  `C.BINARY_OUTCOMES` exactly. A partial map is the mixture arrived at one outcome at a time, and it is
  the failure mode a permissive `paths.get(key, _augmented_path(...))` would introduce while looking
  like a convenience.
- **A `dict[str, str]` of the same three declared strings `_augmented_path` returns**, so Stage 10
  builds it by reading `augmented_path` off the point estimate's seven `BinaryEstimate`s and passing it
  straight back. It never re-derives the rule, which is what §14 requires of it.

`_augmentable` and `_augmented_path` stay public-facing in the sense §14 means — Stage 10 **reads** them
rather than reimplementing them — and with `paths` in the signature it does not need to call them at all.

---

## 12. Data flow into Stages 10-14

```
  secondary(cohort, ps, audit) returns Secondary
   └─ estimates      dict of seven BinaryEstimate, C.BINARY_OUTCOMES order      [§4.1]
       ├─ outcome         the C.OUTCOMES key                                    [§4.1]
       ├─ family          "secondary" | "safety" -- Stage 11 groups on this     [§4.1]
       ├─ rd              the weighted risk difference                          [§5.2]
       ├─ odds_ratio      the weighted marginal odds ratio                      [§6]
       ├─ or_corrected    §6.3's correction fired -- FALSE on all seven on v7   [§6.4]
       ├─ proportion      the two weighted proportions, uncorrected             [§5.1]
       ├─ augmented       tau, or None on the two unaugmented outcomes          [§8]
       ├─ augmented_path  "full" | "reduced" | "unaugmented"                    [§9.2]
       ├─ covariates      m_a(X)'s declared list, or None                       [§7.1]
       ├─ reduced         the list is an override -- Stage 14 must print it     [§7.2]
       ├─ dropped         design columns dropped as constant                    [§7.1]
       ├─ minority        min(events, non-events) on in_estimate                [§9.1]
       ├─ in_estimate     boolean, TOTAL -- one of the seven [§11] denominators [§4.2]
       └─ fit             model.Fit, or None

  the cohort frame and the Propensity come back unchanged                       [§0.2]

  audit: load 7 / derive 4 / classify 1 / build 6 / fit 4 / assess 2
         / primary 3 / secondary 3  =  30
```

**One rule every consumer owes, and it is Stage 6 §9's, unchanged.** `e` and `w` carry `nan` off
`in_model`; range over the mask, never over `notna()`, and never fill. Stage 9 is the fourth consumer
and honours it in **two** lines — `ps.w.loc[in_estimate]` and `ps.e.loc[in_estimate]` in §11 — where
Stage 8 needed one. Both slice the *outcome's* mask, which is a subset of `in_model`, so the rule is
satisfied a fortiori and §15.11 asserts no other read of either exists in the module.

**And the a-fortiori claim is about MEMBERSHIP, not ORDER, which S3a is what covers.** `ps.w.loc[mask]`
returns rows in **`ps.w`'s** order while `m1` and `m0` come from `model.design(df.loc[mask])` in **the
frame's** order. If the two indexes hold the same labels in different sequences, every arm term pairs
each record's weight with another record's fitted value: no raise, no `nan`, a different number.
Measured on a 92-row construction, reversing `ps.e` and `ps.w` moved `rd` from `+0.224831` to `+0.219968`
and `tau` from `+0.227928` to `+0.235780`. `propensity.fit` builds both from `df` so the orders agree on
anything it produces — but nothing in S4 constrained a `Propensity` a caller assembles, and Stage 10
assembles one per replicate. §4.4b's S3a is `index.equals`, which is order-sensitive, and §15.1a's S3a
frame is a permuted-but-equal index (§22.3 item 19).

### 12.1 Stage 10 [§10]

Refits this whole function in every replicate. Five things are owed to it explicitly, and the first is a
blocker.

- **`propensity.fit` raises `SchemaError` on 26.0% of stratified replicates, and Stage 10 has
  prespecified that it must not catch it.** This is §12.3, below. It is not Stage 9's to fix and it is
  Stage 9's to hand over, because this stage is the first to resample the pipeline and find out.
- **The augmentation paths must be passed IN, not recomputed.** §9.5: `ph2` crosses the threshold in
  46.1% of replicates and `sich` in 3.8%. A replicate that decides its own path contributes a different
  estimator to the same interval.
- **Two failure sources, and they are not the ones Stage 8 named.** Stage 8 handed over `FitError` from
  `polr` (separation) versus non-convergence, and measured that the second never fires. Here the
  `FitError` comes from `model.firth` on `m_a(X)` and its measured rate is **0.3% for `sich`, 0% for
  every other outcome**. So Stage 9 contributes almost no droppable replicates — but see the next item,
  because that is not the same as contributing no degenerate ones.
- **`max|β|` per outcome per replicate should be reported as a distribution**, beside the counters. §9.6:
  there is no bound, so nothing turns a separated nuisance fit into a countable failure, and `sich`'s
  `max|β|` reaches 80.98 with 15.1% of replicates above 14. The counters will read near-zero while a
  sixth of one outcome's nuisance models are degenerate. This is Stage 8 §11's `len(fit.alpha)` move,
  and it is the diagnostic that answers §16's open question about whether `sich`'s interval is usable.
- **Stage 9 is the cost centre, and `model.design` is 72% of Stage 9.** Measured: 24.05 ms per
  replicate against `propensity.fit`'s 13.39 ms and `primary`'s 4.73 ms, with the measured loop at
  71 ms wall-clock and about **142 s** over `N_BOOT` = 2000. The four full-list outcomes share an
  identical design matrix, so building it once per replicate is worth roughly 22 s (§2).

### 12.2 Stages 11, 12 and 14

- **Stage 11 [§13]** — applies Benjamini–Hochberg *within* the secondary and safety families. It reads
  `Secondary.by_family()` and does not partition the outcomes itself (§3.1): [§13]'s correction is wrong
  if two places compute the family partition and disagree. What Stage 11 owes is the raw and adjusted p
  both reported, and that the primary stays uncorrected. The E-value is the *primary* estimate's and
  reads Stage 8's `β`, not anything here.
- **Stage 12 [§14a]** — reads nothing from this stage. [§14] fits no propensity model, so there is no
  `h`, no `w` and no augmentation. Recorded so that nobody wires `augmented_rd` into the
  standardisation path looking for a use.
- **Stage 14 [§16]** — reports each estimate with its interval from Stage 10, and owes four guardrails
  this stage cannot enforce: every reduced specification printed beside its estimate [§8 amendment];
  §6.3's correction named wherever it fired; "model-assisted" and never "doubly robust"; and the
  direction of each risk difference attached from `higher_is_better` rather than from the sign
  convention, which is uniform here (§4.1).

### 12.3 `propensity.fit` cannot be bootstrapped as landed, and this is where it was found

`propensity._record_exclusion` (`propensity.py:400-432`) compares the number of **distinct** excluded
`case_id`s against the number of excluded rows and raises `SchemaError` when they differ. Its docstring
states the premise: *"a duplicated or missing case_id is forbidden by Stage 2's A2, so this is a belt
over those braces."*

A2 forbids duplicates **in the workbook**. A patient-level bootstrap replicate contains duplicates by
construction — that is what sampling with replacement is — so the premise does not hold for the
population of frames [§10] feeds this function.

The workbook has exactly **one** covariate-incomplete record, in Lugano (stratum size 31). Whenever it
is drawn twice or more, two excluded rows collapse to one name and the guard fires. Measured over 400
stratified replicates: **104 raised, 26.0%**. The analytic probability that a given row is drawn at least
twice from a stratum of 31 is **26.4%**, and the match identifies the cause exactly.

And Stage 8 §11 handed over to Stage 10 that *"`FitError` is the droppable failure and `SchemaError` is
NOT — a `SchemaError` from `primary` or `propensity.fit` is a bug in the resampler, not a sparse
replicate, and catching it would drop replicates for a reason that is not about the data."* So as landed,
Stage 10 has three options and all three are prohibited: crash on a quarter of its replicates; catch a
`SchemaError` it has prespecified it must not catch; or drop 26% of replicates for a reason that is not
about the data.

**Not fixed here, and the reasons are stated rather than assumed.** The fix is a change to Stage 6's
audit contract — either the guard learns that a replicate's rows are distinct draws, or the resampler
gives each drawn row a distinct name — and both are decisions about what a replicate's log means, which
is [§10]'s subject and not [§8]'s. Making it here would put a Stage 6 amendment inside a Stage 9 commit
on the strength of a Stage 10 requirement that has not been specified yet. What Stage 9 does instead:
records the measurement, hands it over in §14, files it in §16, and puts it in `TODOS.md` with the
trigger. **The probes behind §9.5, §9.6 and §12.1 all had to work around it** by renaming each drawn row
— which is one of the two candidate fixes, adopted provisionally in a scratchpad and not in the
repository. §21 records that those three measurements are therefore conditional on a workaround, which
is a limitation of them and is stated where they are used.

---

## 13. What Stage 9 amends in Stages 1-8

The full ledger, so that no amendment is discovered during implementation. **Three shipped modules and
three test modules, and one of the shipped changes is a docstring sentence.**

| File | Amendment | Why |
|---|---|---|
| `config.py` | `BINARY_OUTCOMES`, immediately after `PRIMARY_OUTCOME` | §4.1. Computed from `OUTCOMES` as the `kind == "binary"` keys, so `outcome.py` cannot write seven keys as literals. Declared beside `PRIMARY_OUTCOME` because the two partition the registry and a reader checking that should not have to look in two places |
| `config.py` | `OR_CONTINUITY`, in the inference-and-thresholds block beside `RARE_MINORITY_THRESHOLD` | §6.3. Prespecified for `FIRTH_*`'s reason (`config.py:281-284`): [§10] refits in every replicate, so a correction that changes an answer is a property of the sampling distribution and not a runtime knob. Placed beside `RARE_MINORITY_THRESHOLD` because both are [§8]'s degenerate-data rules and neither is a tolerance |
| `test_config.py` | **five** assertions | `BINARY_OUTCOMES` has seven entries, every one `kind == "binary"`; it is disjoint from `PRIMARY_OUTCOME` and their union is `OUTCOMES`, which is what stops an outcome being estimated by neither stage or both; every key's `family` is in `("secondary", "safety")`, because §3.1's `by_family` and [§13]'s Benjamini–Hochberg both assume two families and a third would silently get its own correction group; `OUTCOME_MODEL_OVERRIDES`' keys are a subset of `BINARY_OUTCOMES`, so an override cannot be declared for the ordinal outcome or for a typo; and `0.0 < OR_CONTINUITY <= 0.5`, asserted against both endpoints with §6.4 cited |
| `model.py` | `predict`; and **one sentence** in the module docstring | §7.3. The docstring says the module *"names neither the treatment nor any covariate list"*; the sentence records that `predict` cannot tell an exposure from a covariate, so a reader who finds it building counterfactuals does not conclude the rule was bent. **No existing function changes** |
| `outcome.py` | `Secondary`, `BinaryEstimate`, `secondary`, `weighted_rd`, `marginal_odds_ratio`, `outcome_model`, `augmented_rd`, and **the privates §3.2 lists and counts** | §3, §5-§11. The count lives in §3.2 alone; three earlier drafts carried three different numbers here, in §1 and in §3.2 (§22.3 item 5). **No existing function changes** — `weighted_proportion`, `cumulative_rd` and `primary` are read, not edited, and §15.12 asserts the Stage 8 suite passes unedited |
| `tests/fixtures_stage9.py` | **new file, and this row is the one place its contents are enumerated.** Seven functions — `golden_frame`, `golden_arrays`, `tilt_frame`, `separated_frame`, `dr_population`, `ato`, `renamed_replicate` — and four constants: `TILT_SEED`, `DR_SEED`, `DR_SEEDS_12`, `DELTA` | §15.0. Every constructed fixture in this document, as code. Round 3's largest change: three of them were described and never given, while their pins were asserted to ten decimals (§22.3 item 1). **`golden_arrays` is the implementation's and is not in §15.0's fences**: it returns `(y, a, w, e)` from the golden frame with `w` built as [§7]'s overlap weight — `1 − e` treated, `e` control — because that construction is a two-line expression every §15 test needs and a test that writes it the other way round measures a different estimand while every assertion in it still reads plausibly. One construction, one place to be wrong (§22.4) |
| `test_outcome.py` | the §15 sections, banner-commented `# --- 15.x` | Its existing sections are untouched |
| `test_model.py` | `predict`'s sections and the separated-frame witness | §15.6, §15.7. The witness lives here because the function does |
| `test_reference_r.py`, `tests/reference/aug_psweight.R` | the Stage 9 oracle, behind the gate Stage 7 built | §19b |
| `implementation_roadmap.md` | Stage 9 gains its `**Spec:**` line, four corrections and one addition | §22. Lands with this document |
| `TODOS.md` | four items with their triggers | §16 |

**Five things that look like they need amending and do not.**

- **`data.py`.** §10.1: no new kind, no new heading, so none of the literal pins in the existing test
  modules moves except the ledger count. `KINDS` stays nine and §15.14 asserts it. Third stage running.
- **`propensity.py`.** `Propensity` already carries `e`, `w` and `in_model`, which is everything this
  stage reads. Nothing is added. **And `_record_exclusion` is deliberately not touched** despite §12.3 —
  that is a Stage 6 amendment on the strength of a Stage 10 requirement, and it does not belong in this
  commit (§12.3's closing paragraph).
- **`RARE_MINORITY_THRESHOLD`.** Declared at Stage 1 for this stage, recorded by Stage 8 §12 as Stage
  9's, and **read rather than written** — exactly as `SMD_THRESHOLD` was at Stage 7 and `MRS_THRESHOLDS`
  at Stage 8. Its value is not touched, and §9.4 records that `ph2` landing one below it was looked at
  and left alone.
- **`OUTCOME_MODEL_OVERRIDES` and `outcome_model_covariates`.** Written at Stage 8 §12 for this stage,
  unread until now, and correct as written — the accessor already excludes treatment, already defaults,
  and already raises on an unregistered key. §15.9 asserts the default and the override both, so the
  first read of a function written one stage early is a tested read.
- **`model.Fit` and `model.firth`.** Untouched. `predict` is a reader of `Fit`, not a change to it, and
  `firth`'s six numerical details are inherited unaltered.

```python
# config.py — the two additions, as they are to be pasted. §13
#
# The [§5] binary outcomes, computed from the registry rather than named a second time. PRIMARY_OUTCOME
# above takes the one primary entry; this takes the seven binary ones, and test_config.py asserts the
# two partition OUTCOMES exactly — so an outcome cannot be estimated by neither Stage 8 nor Stage 9,
# and cannot be estimated by both. Declared immediately after PRIMARY_OUTCOME for that reason: the two
# are one decision about the registry and a reader checking the partition should find them together.
BINARY_OUTCOMES: Final[tuple[str, ...]] = tuple(
    key for key, o in OUTCOMES.items() if o.kind == "binary")

# The [§8] weighted marginal odds ratio's continuity correction [Stage 9 §6.3]. Haldane-Anscombe:
# added to all four weighted pseudo-counts, and ONLY when a weighted proportion reaches 0 or 1, so an
# interior estimate is bit-for-bit uncorrected. PI decision, 2026-08-24 — neither [§8] nor any
# amendment specifies a rule and the roadmap says only "kept finite".
#
# Prespecified for FIRTH_*'s reason (config.py:281-284): [§10] refits this in every one of N_BOOT
# replicates, so a correction that changes an ANSWER is a property of the sampling distribution and
# not a runtime knob.
#
# The scale is stated rather than implied, and PER ARM rather than averaged. These are WEIGHT-SUMS,
# not row counts: Sw = 27.736623 over the ATO population, splitting 13.626360 treated / 14.110263
# control against row counts of 39 treated / 53 control. The conventional 0.5 is calibrated against
# ROWS, so against these weight-sums it is 2.9x more aggressive in the treated arm and 3.8x in the
# control arm — measured per arm, because averaging the two row counts to "about 46" and reporting
# "about three times" hides that the two arms are corrected by different relative amounts, which is
# exactly the asymmetry that lets the correction cross the null [Stage 9 §6.4].
#
# It is NOT scaled to Sw to compensate, because a correction whose magnitude is a function of the
# weights is a correction whose magnitude is a function of the propensity model, and [§10] refits
# that in every replicate.
#
# CONSEQUENCE, MEASURED, AND IT IS NOT A ROUNDING EFFECT: because the arms are corrected unequally,
# an empty cell in the HEAVIER arm can carry the odds ratio ACROSS 1 rather than toward it — e.g.
# Sw1 = 17.007 against Sw0 = 22.768 with p0 = 0.00647 returns 1.0201 where the uncorrected value is
# 0.0. About 3 replicates in N_BOOT = 2000 report a safety outcome with no bridging events as
# favouring EVT alone. PI-reversible; [Stage 9 §6.4, §17].
#
# Measured: the branch is unreachable on v7 — no arm holds an empty cell for any of the seven
# outcomes — and fires in 17.0% of stratified replicates for sich, 13.0% for tici_2b_3 and 2.3% for
# ph2. So this constant is materially a choice about two INTERVALS and barely a choice about any point
# estimate [Stage 9 §6.4].
OR_CONTINUITY: Final[float] = 0.5
```

---

## 14. Handover to Stage 10

```
  cohort = cohort.build(eligibility.classify(derive.derive(*data.load()), ...), ...)
  ps     = propensity.fit(cohort, audit)
  bal    = balance.assess(cohort, ps, audit)        # not read by Stage 8 or Stage 9
  est    = outcome.primary(cohort, ps, audit)       # not read by Stage 9
  sec    = outcome.secondary(cohort, ps, audit)

    93 / 92          cohort / in_model                                        [§4.2]
    92 x6, 91 x1     the seven [§11] denominators; tici_2b_3 is the 91        [§4.2]
    5 / 2            augmented (4 full + 1 reduced) / unaugmented             [§9.3]
    3                audit entries, taking the ledger 27 -> 30                [§10.1]
    0.0%             replicates in which the OR correction fires on v7        [§6.4]
    42 ms / 71 ms    fitting / wall-clock per replicate; STAGE 9 IS THE
                     COST CENTRE, and model.design is 72% of Stage 9         [§2]

    every weighted proportion, RD, OR and tau on the workbook is DELIBERATELY
    ABSENT from this ledger and from this document. They are in the
    gitignored log.                                                          [§4.3]
```

**Five things Stage 10 must build that this stage cannot, the first four in the order they will bite.**

1. **A resampler `propensity.fit` does not raise on.** §12.3: 26.0% of stratified replicates raise
   `SchemaError`, and [§10] has prespecified that a `SchemaError` is a resampler bug and must not be
   caught. Either the guard learns that a replicate's rows are distinct draws, or each drawn row gets a
   distinct name. **This is a blocker and it is the first thing Stage 10 will hit.**
2. **The augmentation paths, passed in — through the parameter §11 now has.** §9.5: computed once on the
   point estimate, carried into every replicate, never recomputed. `ph2` crosses in 46.1% of replicates.
   Stage 10 reads `augmented_path` off the point estimate's seven `BinaryEstimate`s and passes the map
   back as `secondary(..., paths=frozen)`. It does not re-derive the rule (§11.1).
3. **Two counters and one distribution.** `FitError` from `model.firth` per outcome (measured: 0.3% for
   `sich`, 0% elsewhere) — and the **distribution of `max|β|` per outcome**, because §9.6 prescribes no
   bound and the counters therefore read near-zero while 15.1% of `sich`'s nuisance fits are above 14.
4. **A throwaway `Audit` per replicate.** Three entries per replicate at `N_BOOT` = 2000 is 6000, one
   carrying a seven-row rendered table.
5. **One design per replicate, not four.** The four full-list outcomes share an identical design matrix
   and `model.design` is 72% of Stage 9's cost; reusing it is worth about 22 s over `N_BOOT` (§2).

What Stage 10 inherits and should not rebuild: every estimator here is public and separately callable
(§3.2), so a replicate loop can call **`secondary` whole, passing `paths`**, or the three estimators
individually; `predict` for any re-derivation of `m_a`; and `_augmentable`'s rule, which it must **read**
rather than re-implement, since re-implementing it is how item 2 gets violated by accident. An earlier
draft said the replicate loop *"can call `secondary` whole"* without qualification, at a point when doing
so necessarily recomputed the paths — the two sentences contradicted each other and §11.1 is the
resolution.

**What Stage 10 must not do**, each because a prespecified rule says so: substitute an unaugmented
estimate when an augmented fit fails [§10]; catch `SchemaError` [Stage 8 §11]; correct an odds ratio
that is `nan` rather than degenerate (§6.2); or take percentiles of `tau` and of `rd` from different
replicate sets, since the two are reported as a comparison (§10.3) and a comparison across two different
resamples is not one.

---

## 15. Acceptance criteria

Tests live in `test_outcome.py`, with section banners matching these numbers, as `test_cohort.py`,
`test_model.py`, `test_propensity.py`, `test_balance.py` and Stage 8's sections do. The `predict`
sections — §15.6 and §15.7 — live in `test_model.py`, because that is where the function is. Tests
needing the workbook are gated as Stage 7 built the gate and are tagged `[data-gated]`.

### 15.0 The frames and fixtures this stage is tested on

**Every constructed fixture is code in this document, and that is round 3's largest change.** An earlier
draft *described* fixtures 1, 4 and 5 and pinned quantities measured on them to ten decimal places — so
`tau = 0.2252774207` was an instruction to reproduce a number from a frame that appeared nowhere, while
§21's own rule forbade re-deriving it and §0's status block claimed anything not in the document must not
be invented. The frames live in `tests/fixtures_stage9.py`, every fence below is their specification, and
**every pin in this document was re-measured against them** — so the pins are not the earlier draft's.
§22.3 items 1 and 2.

**Eight, not seven, and the eighth is the one round 3 said it could not write.** `secondary` takes a
cohort-level frame, and every fixture round 3 specified is outcome-level — arrays with a frame around
them. §15.0.7 is that gap closed, and it is the fixture T7 and T9 could not have been written without.
§22.4 item 5.

**1. The golden frame (§15.0.2)** — 24 records, every value a literal, no RNG anywhere. Three centres,
one continuous covariate, one binary covariate, a written-out propensity score. It carries no fit, so
the vector is a property of this document rather than of a numpy version's random stream.

**2. Hand-built weighted frames** — four to eight rows with weights and outcomes chosen so that the
weighted proportions are exact decimals, for checking `weighted_rd` and `marginal_odds_ratio` against
arithmetic done by hand rather than against another implementation. Built inline in the test that uses
them, because each is three lines and naming them centrally would hide which arithmetic each checks.

**3. Degenerate frames (§15.0.3)** — an arm with no events, an arm with no non-events, a constant
outcome, an arm with zero total weight, a proportion that is `nan`, and — round 3's addition — **an empty
cell in the heavier arm**, which is the frame on which the correction crosses the null (§6.4). One frame
per branch, built inline.

**4. The tilt-separation frame (§15.0.4)** — 92 records with a deliberately heterogeneous `m_a(X)`, and
**4a. the separated frame (§15.0.4a)** — 40 records, perfectly separated on treatment.

**5. The double-robustness population (§15.0.5)** — 200,000 records, three arms of the construction,
seeded, with `DELTA` and the twelve seeds written out.

**6. The workbook, data-gated** — properties only, never values (§4.3).

**7. The renamed replicate (§15.0.6)** — the `case_id`-renaming helper §12.3's measurements were taken
through, shipped so that §21's four caveated rows are reproducible rather than resting on a scratchpad.

**8. The cohort-level frame (§15.0.7)** — 40 records carrying the nine [§6] covariates and all seven
outcome columns, which is the only shape `secondary` itself can be called on. **Round 3 recorded this as
the fixture it could not supply** — *"building a cohort-level one means either the workbook or a fixture
this document still does not specify. So the six ordering constraints, the three audit entries and S1-S8
remain specified and unexecuted"* (§22.3, closing). It lives in `test_outcome.py` rather than in
`fixtures_stage9.py` for the reason §15.0.7 gives, and §22.4 item 5 records what building it found.

#### 15.0.2 The golden frame and the golden vector

```python
def golden_frame() -> pd.DataFrame:
    """24 records, every value a literal, no RNG. The Stage 9 regression pin [§15.0.2].

    Three centres and NOT four: USZ is absent, so `model.design` drops `center_USZ` exactly as it
    does on the workbook (§3.3), and the fixture exercises the constant-column rule rather than
    stepping around it. `age` is the one continuous covariate and `atrial_fib` the one binary one,
    which is the minimum that makes the two covariate lists of §15.0.2's pair differ by a column.

    `e` is WRITTEN OUT rather than fitted. The frame therefore carries no propensity model, and the
    golden vector is a property of this document rather than of a numpy version's random stream
    [Stage 7 §12.0.2]. `w` is [§7]'s overlap weight, `1 - e` treated and `e` control, computed by
    the caller from this column.

    `age` is deliberately near-separating on `y` -- 80/76 on the events, 53/57 on the non-events --
    because the pair's whole purpose is to show a six-parameter model over-fitting a minority cell
    of 11 (§15.0.2's closing paragraphs). A frame on which the two fits agreed would pin the
    arithmetic and demonstrate nothing.
    """
    return pd.DataFrame({
        "case_id":    [f"G{i:02d}" for i in range(1, 25)],
        "center":     ["HUG"] * 9 + ["CHUV"] * 8 + ["Lugano"] * 7,
        "age":        [80.0, 53.0, 57.0, 76.0, 53.0, 80.0, 76.0, 53.0, 57.0,
                       80.0, 53.0, 80.0, 76.0, 53.0, 57.0, 76.0, 80.0,
                       53.0, 76.0, 53.0, 57.0, 53.0, 80.0, 53.0],
        "atrial_fib": [0.0, 1.0, 0.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0,
                       1.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 0.0,
                       0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 1.0],
        C.TREATMENT:  [1.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0,
                       1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 0.0, 1.0,
                       0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0],
        "y":          [1.0, 0.0, 0.0, 1.0, 0.0, 1.0, 1.0, 0.0, 0.0,
                       1.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0,
                       0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0],
        "e":          [0.62, 0.55, 0.41, 0.70, 0.33, 0.58, 0.46, 0.66, 0.29,
                       0.74, 0.37, 0.61, 0.52, 0.44, 0.26, 0.49, 0.68,
                       0.31, 0.72, 0.39, 0.64, 0.22, 0.57, 0.35],
    })
```

Two pins, one per augmentation path, on that frame. Its minority cell is 11, so it is augmentable under
§9.2, and its two covariate lists differ by one continuous column:

```
  the frame          n 24, events 11, non-events 13, minority 11
                     rows 12 treated / 12 control
                     Sw  8.8900000000        Sh  5.4257000000
                     p1_w  0.6191536748      p0_w  0.3340909091
                     RD_w  0.2850627657      OR_w  3.2404025938   or_corrected False

  covariates = ("center", "atrial_fib")            -- the REDUCED shape, TICI's list
                     design 4 cols (ivt, atrial_fib, center_CHUV, center_Lugano)
                     dropped ("center_USZ",)
                     iterations 7, converged_on "likelihood", max_abs_beta 1.5877919570
                     rescales 0, halvings 0, first_step_norm 2.4935249391
                     tau           0.2233300473
                     tau under w   0.2240676700      gap 7.376226e-04
                     predict(fit, X) == fit.p        EXACTLY

  covariates = ("center", "atrial_fib", "age")     -- one continuous column added
                     design 5 cols (ivt, atrial_fib, age, center_CHUV, center_Lugano)
                     dropped ("center_USZ",)
                     iterations 10, converged_on "likelihood", max_abs_beta 11.8282581834
                     rescales 1, halvings 0, first_step_norm 10.7967589095
                     tau           0.0402879733
                     tau under w   0.0402999125      gap 1.193917e-05
                     predict(fit, X) == fit.p        EXACTLY
```

**The two pins are a pair and the second one is the point.** Adding one continuous covariate to a
24-record frame moves `tau` from `+0.2233` — within 0.062 of the unaugmented `+0.2851` — to **`+0.0403`,
a seven-fold collapse toward zero**. Six parameters against a minority cell of 11 over-fits: the fitted
values approach `Y`, the correction term stops carrying signal, `max_abs_beta` goes from **1.588 to
11.828**, the iteration count goes from 7 to 10, and **the trust region fires once**.

This is [§8]'s *"material disagreement between the two is evidence about the outcome model, not
confirmation of either"* as an executable demonstration, and it is why §10.3 prints `RD_w` and `tau` in
adjacent columns. It is also a live argument for the rare-outcome rule the amendment tightened: the
failure mode that rule guards against is visible at a minority cell of 11, which is *above* the
threshold of 10.

The pair also gives the suite one fit with `rescales == 0` and one with `rescales == 1`, so the trust
region is exercised in both states by the regression pin rather than only by §15.7's constructions.

**What changed from the earlier draft, and what did not.** The earlier draft's pins were `RD_w 0.2266`,
`tau +0.2253` and `tau −0.2050` — a *sign flip* on the second fit — measured on a frame it did not
contain. Those numbers are gone, because no frame in this document produces them and a pin nobody can
reproduce is not a pin. The **property** the pair demonstrates survives and is if anything cleaner: the
earlier frame showed over-fitting by inverting the estimate, this one shows it by collapsing the estimate
while the trust region fires and `max_abs_beta` moves by an order of magnitude. What is *not* claimed is
a sign flip — §21 records the measured collapse, and an implementer who reproduces `+0.0403` has
reproduced the fixture.

#### 15.0.3 The degenerate frames

Six frames, one per branch of §6.2 plus round 3's addition. Built inline in `test_outcome.py` because
each is a handful of rows and naming them centrally would hide which branch each reaches — but the
**arm weight totals are the point** on **three** of them, so they are specified here rather than left
to whoever writes the test. **Round 3 said two; the third was found by building frame 7 and watching
it return `1.0` (§22.4 item 2).**

```
  #  frame                                          Sw treated / control   reaches
  1  no events in the treated arm                        even, ~4 / ~4     p1 == 0
  2  no non-events in the treated arm                    even, ~4 / ~4     p1 == 1
  3  no events in the control arm                        even, ~4 / ~4     p0 == 0
  4  no non-events in the control arm                    even, ~4 / ~4     p0 == 1
  5  zero total weight in one arm                        0 / ~4            p is nan, NOT corrected
  6  no events in the treated arm, treated arm HEAVIER   17.007 / 22.768   the NULL CROSSING (§6.4)
  7  outcome constant on the population                  8 / 11            finite OR, not nan (§6.2)
```

**Frame 6 is round 3's and its weights are not free.** The crossing needs the arm holding the empty cell
to carry more total weight than the other, and the effect is a function of *how much* more: at the
workbook's near-even split (13.626 / 14.110) the corrected odds ratio is `0.6484` and does not cross,
while at 17.007 / 22.768 with `p0 = 0.00647` it is `1.0201` and does. So a frame built with round
weights would test the branch and miss the property. The two rows §6.4 tabulates are the two this frame
reproduces.

**And frame 7's weights are not free either, which round 3 got wrong and only building it found.** That
draft listed it as "even, ~4 / ~4" while §6.2 pins the values it produces to ten decimals — and the two
cannot both hold. On a constant outcome **every event pseudo-count is `OR_CONTINUITY` alone**, so the
degenerate branch reduces to

```
  all-zero outcome:  OR = (0.5 / (Sw0 + 0.5)) / (0.5 / (Sw1 + 0.5))  =  (Sw0 + 0.5) / (Sw1 + 0.5)
  all-one  outcome:  OR = ((Sw1 + 0.5) / 0.5) / ((Sw0 + 0.5) / 0.5)  =  (Sw1 + 0.5) / (Sw0 + 0.5)
```

— a pure function of the two arm totals and of nothing else in the frame. **Even totals therefore give
exactly `1.0`**, which is a degenerate special case that hides the very asymmetry §6.4 exists to expose,
and no frame with "even, ~4 / ~4" can produce §6.2's numbers. `Sw1 = 8` and `Sw0 = 11` give
`11.5 / 8.5 = 1.3529411765` and `8.5 / 11.5 = 0.7391304348`, which are §6.2's two measured values —
and **the fact that they are exact reciprocals is itself the statement** that the two directions of
degeneracy are one piece of arithmetic read twice. §15.5 asserts the reciprocal identity for that reason
rather than asserting the two numbers independently.

Frames 1-5 need only that the arms be non-degenerate in weight, so their totals are unconstrained and
stated as "even" to make clear that nothing rests on them.

#### 15.0.4 The tilt-separation frame, and why it was measured rather than reasoned

§8.4 established that the `h`-versus-`w` gap is driven by the variation in `d = m₁ − m₀` and that it is
**not monotone** in that variation. So the frame cannot be argued into existence; it has to be measured.

```python
TILT_SEED: Final[int] = 20260824

def tilt_frame(seed: int = TILT_SEED, n: int = 92,
               x_scale: float = 3.0, interaction: float = 3.0):
    """The §15.8 item 3 fixture: 92 records on which the h- and w-tilts measurably differ.

    Returns (x, e, a, m1, m0, y). `m1` and `m0` are the ORACLE conditional means and not a fit,
    because the test pins the tilting function and a nuisance-estimation error between the two calls
    would be indistinguishable from a tilt difference.

    `x_scale = 3.0` and `interaction = 3.0` are the `wide x, interaction 3.0` shape of §8.4's
    six-candidate table -- the row with the largest measured gap that was not the non-monotone
    outlier. The seed is pinned because the gap is a property of the draw as well as of the
    construction: §8.4's table shows one construction with the LARGEST sd(d) producing the
    SECOND-SMALLEST gap, so "make it heterogeneous" is not a specification.

    `e` is clipped into [0.05, 0.95] so that h = e(1-e) is bounded away from zero and Sh does not
    collapse onto a handful of records -- which would make the gap a property of those records.
    """
    g = np.random.default_rng(seed)
    x = g.normal(0.0, x_scale, n)
    e = np.clip(1.0 / (1.0 + np.exp(-(0.8 * x))), 0.05, 0.95)
    a = (g.random(n) < e).astype(float)
    m1 = 1.0 / (1.0 + np.exp(-(-0.3 + 0.5 * x + interaction * 0.5)))
    m0 = 1.0 / (1.0 + np.exp(-(-0.3 + 0.5 * x - interaction * 0.5)))
    y = np.where(a == 1.0, (g.random(n) < m1).astype(float),
                 (g.random(n) < m0).astype(float))
    return x, e, a, m1, m0, y
```

Measured on it:

```
  seed 20260824, n 92, x ~ N(0, 3^2), interaction 3.0
    sd(d)  0.1698        range(d)  0.5737
    Sw    25.913563      Sh       13.083855
    tau under h   0.5112006369
    tau under w   0.5008939927
    |tau_h - tau_w|  =  1.030664e-02

  sanity, d constant at 0.00 / 0.05 / 0.40:  gap 0.000e+00 / 0.000e+00 / 0.000e+00
```

The test asserts the `h` value and **computes and prints the `w` value**, asserting they differ by more
than `1e-03`. Two properties of that bound: it is **10× below the measured gap**, and the three
constant-`d` sanity rows cancel **exactly**, at `0.000e+00`, so there is no floating-point floor to clear.
A tolerance loose enough to survive a numpy bump but tight enough to fail on the wrong tilt has to sit in
that window, and the window is measured rather than assumed.

#### 15.0.4a The separated frame

```python
def separated_frame(n: int = 40) -> pd.DataFrame:
    """40 records, PERFECTLY separated on treatment: y == C.TREATMENT on every row. §9.6, §15.7.

    `center` alternates independently of treatment and `atrial_fib` runs on a period of four, so the
    design is FULL RANK. A first draft made `center` track the arms -- HUG control, CHUV treated --
    and `firth` raised F3 "the design has rank 3 of 4 columns" rather than converging, which is a
    collinearity failure and not the separation this fixture exists to exhibit.

    Only two of the four declared centres appear, so `design` drops `center_Lugano` AND `center_USZ`
    and the frame exercises the constant-column rule twice.
    """
    return pd.DataFrame({
        "case_id":    [f"S{i:02d}" for i in range(n)],
        "center":     ["HUG", "CHUV"] * (n // 2),
        "atrial_fib": [0.0, 0.0, 1.0, 1.0] * (n // 4),
        C.TREATMENT:  [0.0] * (n // 2) + [1.0] * (n // 2),
        "y":          [0.0] * (n // 2) + [1.0] * (n // 2),
    })
```

Measured, with `covariates = ("center", "atrial_fib")`:

```
  design 3 cols (ivt, atrial_fib, center_CHUV), dropped ("center_Lugano", "center_USZ")
  iterations 6, converged_on "likelihood", rescales 0, halvings 0
  beta_treatment 6.0890571141      exp(beta_treatment) 441.005397
  max_abs_beta   6.089057
  fitted p strictly interior, in [4.545e-02, 0.954546]

  the UNPENALISED MLE on the same frame, and the CAP IS PART OF THE RESULT:
    Newton, maxiter 200            raises LinAlgError: Singular matrix
    Newton, maxiter 35 (the        max_abs_beta 67.8705, converged False -- NO RAISE
      statsmodels default)
    Newton, maxiter 100            max_abs_beta 67.4529, converged False
    BFGS capped at 60 iterations   max_abs_beta 33.0927, converged TRUE -- a reported success
```

#### 15.0.5 The double-robustness population

```python
DR_SEED: Final[int] = 20260825
DR_SEEDS_12: Final[tuple[int, ...]] = tuple(20260825 + 7 * i for i in range(12))
DELTA: Final[float] = 0.12          # the constant RISK DIFFERENCE of the companion arm

def dr_population(n: int, seed: int, arm: str):
    """§15.10's population. `arm` is "heterogeneous" | "constant" | "no_interaction".

    Returns (e_true, e_mis, a, m1, m0, y) with m1, m0 the ORACLE conditional means -- the strongest
    available form of "a correct outcome model", which removes nuisance-estimation error from a
    result about the PROPENSITY model being wrong (§8.5).

    `e_mis` OMITS x2 and `e_true` does not, which is the misspecification. The three arms differ in
    how tau(X) = m1 - m0 varies, and that is the whole design:

      * "heterogeneous"  -- x2 enters m1 only, so tau(X) varies strongly. sd(tau(X)) 0.2435.
      * "constant"       -- m1 = m0 + DELTA exactly, so tau(X) is CONSTANT ON THE RISK-DIFFERENCE
                            SCALE. m0 is bounded into (0.25, 0.70) so m1 stays a probability;
                            without the bound the difference truncates at the top and the arm is no
                            longer constant where it matters. sd(tau(X)) 0.0000.
      * "no_interaction" -- a constant LOG-ODDS effect and therefore a NON-constant risk difference.
                            This is the companion a reasonable draft would have written and it is
                            DETECTABLY BIASED at |t| = 21.3 (§8.5). It exists to be asserted biased.

    The three arms share x1, x2 and the treatment draw for a given seed, because they differ in the
    outcome model and not in the population.
    """
    g = np.random.default_rng(seed)
    x1, x2 = g.normal(0, 1, n), g.normal(0, 1, n)
    e_true = 1.0 / (1.0 + np.exp(-(0.5 * x1 + 1.5 * x2)))
    e_mis = 1.0 / (1.0 + np.exp(-(0.5 * x1)))
    a = (g.random(n) < e_true).astype(float)
    if arm == "heterogeneous":
        m1 = 1.0 / (1.0 + np.exp(-(-0.5 + 0.6 * x1 + 1.6 + 1.5 * x2)))
        m0 = 1.0 / (1.0 + np.exp(-(-0.5 + 0.6 * x1)))
    elif arm == "constant":
        m0 = np.clip(1.0 / (1.0 + np.exp(-(-0.5 + 0.6 * x1 + 0.5 * x2))), 0.25, 0.70)
        m1 = m0 + DELTA
    else:
        m1 = 1.0 / (1.0 + np.exp(-(-0.5 + 0.6 * x1 + 0.5 * x2 + 1.6)))
        m0 = 1.0 / (1.0 + np.exp(-(-0.5 + 0.6 * x1 + 0.5 * x2)))
    y = np.where(a == 1.0, (g.random(n) < m1).astype(float),
                 (g.random(n) < m0).astype(float))
    return e_true, e_mis, a, m1, m0, y


def ato(e: np.ndarray, m1: np.ndarray, m0: np.ndarray) -> float:
    """The [§7] ATO indexed by whichever score it is given -- which is §8.3's whole point."""
    h = e * (1.0 - e)
    return float(np.sum(h * (m1 - m0)) / np.sum(h))
```

§8.5 carries the measurements. Cost is **13 ms per draw** at `n = 200,000`, so §15.10 is two draws and
about 0.03 s.

#### 15.0.6 The renamed replicate, shipped so §21 is reproducible

§12.3's blocker means the replicate-based measurements behind §9.5, §9.6, §6.4 and §12.1 could not be
taken through the landed pipeline. They were taken with each drawn row given a distinct `case_id`, which
is one of the two candidate fixes. **That helper is now in the repository** rather than in a scratchpad,
so §21's four caveated rows can be re-run by anyone:

```python
def renamed_replicate(df: pd.DataFrame, seed: int,
                      stratum: str = "center") -> pd.DataFrame:
    """One stratified patient-level replicate, with every drawn row given a distinct case_id.

    NOT A FIX and not proposed as one. `propensity._record_exclusion` raises SchemaError on 26.0% of
    stratified replicates because a resample contains duplicate case_ids by construction (§12.3),
    and choosing between the two candidate repairs is a decision about what a replicate's log means
    -- which is [§10]'s subject. This helper adopts the renaming candidate FOR MEASUREMENT ONLY, and
    it lives in tests/ rather than in a shipped module for exactly that reason.

    It is here because four rows of §21 rest on it. A measurement whose apparatus is not in the
    repository is not a measurement anyone can check, and an earlier draft left it in a scratchpad.

    If Stage 10 adopts the OTHER candidate -- teaching the guard that a replicate's rows are distinct
    draws -- these numbers should be unchanged, because renaming touches only the audit entry and not
    e, w or any mask. That is an argument and §21 labels it as one.
    """
    g = np.random.default_rng(seed)
    out = []
    for _, rows in df.groupby(stratum, sort=True):
        take = rows.iloc[g.integers(0, len(rows), len(rows))].copy()
        out.append(take)
    rep = pd.concat(out, ignore_index=True)
    rep["case_id"] = [f"{cid}#{i}" for i, cid in enumerate(rep["case_id"])]
    return rep
```

`groupby(..., sort=True)` and a stratum-ordered loop, so a replicate is a function of the seed alone and
two runs of the same experiment draw the same frames — the `separate-rng-streams` lesson from Stage 6,
which found that sharing one generator between the resampler and a per-replicate probe made two runs of
one measurement incomparable.

#### 15.0.7 The cohort-level frame, which is the one round 3 could not write

`secondary` takes a **cohort-level** frame: the nine [§6] covariates, `C.TREATMENT`, `case_id` and all
seven outcome columns. Every fixture above is **outcome-level** — `golden_frame()` carries one response
called `y` and two covariates; `tilt_frame()` returns six arrays. So none of them can be passed to the
function this stage exists to specify, which is why §22.3's closing paragraph records the six ordering
constraints, the three audit entries and S1-S8 as *"specified and unexecuted"*.

It lives in **`test_outcome.py` and not in `fixtures_stage9.py`**, which is a departure from every other
fixture here and is deliberate: `fixtures_stage9.py` is the file this document specifies fence by fence,
and this frame is the one construction the document did **not** specify. Putting it there would blur
which fixtures are the spec's and which are the implementation's. §15.1a's ten broken frames all derive
from it, and they live beside it for the same reason.

Four properties, each chosen against an alternative:

- **Deterministic, no RNG.** Every column is a literal or a cycle over one, so the frame is a property of
  the test file rather than of a numpy version's random stream [Stage 7 §12.0.2].
- **Three centres, USZ absent.** `model.design` therefore drops `center_USZ` from every design exactly as
  it does on the workbook (§3.3), so the fixture exercises the constant-column rule rather than stepping
  around it. §15.9's USZ companion is the frame that puts one back.
- **The outcome mask bites.** One record's `tici_2b_3` is missing, so its [§11] population is 39 against
  the ATO's 40 — Stage 8 §4.1's condition, under which an implementation that never built the mask is
  green (§4.2, §15.2).
- **The seven minority cells straddle the threshold**, at 20, 12, 4, 4, 8, 12 and 20, so all three of
  §9.2's paths are reached without the workbook.

**AND THE COVARIATE PATTERNS HAD TO BE CHOSEN BY MEASURING THE DESIGN RANK, WHICH IS THIS FIXTURE'S ONE
NON-OBVIOUS CONSTRAINT.** A rank-deficient design does not produce a bad estimate — `model.firth` raises
F3 and no estimate exists at all, so every augmented outcome fails at once with a message about
collinearity. Three dependencies had to be broken and **all three looked innocent**:

```
  sex = i % 2                         IS the treatment, which is also i % 2
  prestroke_mrs = i % 3               is a linear combination of the two centre dummies,
                                        because `center` also cycles on three
  atrial_fib = (i // 2) % 2           with an onset_type of period four, makes
                                        ivt - atrial_fib - onset_unwitnessed + onset_wake_up == 0
```

Hence the periods 5, 7 and 3. Measured, the resulting design is rank 13 of 13. **The response patterns
carry the same constraint one level on**: `ivt` alternates on `i % 2`, so an outcome written as
`i % 2 == 0` is perfectly separated on the arm — every event in one arm, `p1 = 0`, `p0 = 1` — and every
estimate on it takes the degenerate branch. Two of the seven wanted a 20/20 split and the obvious way to
write one is exactly that mistake, so the periods are 4, 5, 10 and 20 and every outcome is asserted to
reach both arms in both cells. §22.4 item 5.

### 15.1 The seven outcomes and the registry

- `secondary` estimates exactly `C.BINARY_OUTCOMES`, in that order, and the returned dict's keys equal
  it — so an outcome added to `OUTCOMES` appears here without an edit, and one removed disappears.
- **A companion shows what the computed set PREVENTS**: a frame carrying an eighth binary outcome column
  that is *not* in `OUTCOMES` is estimated for none of it, and an `OUTCOMES` entry whose column is
  missing from the frame raises S1 rather than being skipped.
- `by_family()` partitions the seven into exactly two families, and the union of its values is the
  seven. `[roadmap]`

### 15.1a The preconditions, one frame per branch

Eight branches, eight frames, each derived from `golden_frame()` by breaking exactly one thing. Every one
asserts a `SchemaError` **whose message names its own condition** — not merely that something raised —
because a collected assertion that reports the wrong branch is worse than one that reports nothing.

An earlier draft pointed T7 at §15.1, which specified a test for S1 alone, so seven guards written for
Stage 10 were exercised by nothing (§22.3 item 6). They are unreachable on the workbook (§4.4), which is
precisely why the suite is the only thing that can reach them.

```
  branch  phase  the frame                                     asserts
  ------  -----  --------------------------------------------  ------------------------------
  S1        1    drop the `mrs_0_1_90d` column                 SchemaError naming S1 + the key
  S5a       1    drop the C.TREATMENT column                   SchemaError naming S5a
  S3a       1    ps.e and ps.w REVERSED -- the same labels     SchemaError naming S3a
                 in the other order, which is §12's case
  S3b       1    ps.in_model as int64 0/1 instead of bool      SchemaError naming S3b
  S2        2    set one outcome cell to 2.0                   SchemaError naming S2 + the key
  S4        2    ps.e finite on a row that is OFF in_model     SchemaError naming S4
  S5b       2    C.TREATMENT = nan on one row                  SchemaError naming S5b
  S6        3    ps.in_model all False, AND e/w all nan        SchemaError naming S6
  S7        3    w = 0.0 on every treated row                  SchemaError naming S7 + the arm
  S8        3    an outcome constant on its [§11] population   SchemaError naming S8 + the key
```

**Three frames have to be built more carefully than "break one thing", and all three were found by
running them rather than by reading.** Every one of them is the same trap: an earlier phase's check
firing on damage the frame did incidentally, so the branch under test is never reached.

- **`S6` must blank `e` and `w` as well as `in_model`.** With `in_model` all-False and `e` left
  populated, every row sits off the mask carrying a finite value, so S4 fires in phase 2 — measured, the
  raise reads `S4  ps.e carries 39 value(s) OFF in_model`.
- **`S7` must zero `w` only on rows that are IN `in_model`.** A blanket `w = 0.0` on every treated row
  also overwrites the off-mask `nan`, so again S4 fires first — measured, `S4  ps.w carries 1 value(s)
  OFF in_model`. The correct mutation is `w.mask((treated) & in_model, 0.0)`.
- **`S3a` must be a permutation, not a fresh `RangeIndex`.** Relabelling would be caught for the wrong
  reason; a reversal keeps the labels and changes only the order, which is the failure §12 measures
  moving `rd` from +0.224831 to +0.219968.

**The general rule, since it will bite whoever adds a ninth precondition:** a frame that tests phase-`k`
must be *pristine* for phases 1 through `k−1`. §15.1a asserts not merely that a `SchemaError` is raised
but that **the message names the intended branch**, which is the only assertion that catches this — a
test asserting `pytest.raises(SchemaError)` passes on all three of the misconstructed frames above.

Four properties of the section, and the last two are the ones that make it a test rather than a list:

- **The two frame-level phases collect.** A frame breaking S2, S3b and S5b at once produces **one**
  `SchemaError` naming all three, asserted by counting the named conditions in the message. A phase that
  raises on the first failure reports one problem per run and makes fixing a bad frame an iterative
  guessing game (Stage 6 §4.5).
- **The per-outcome phase names the outcome.** S6, S7 and S8 fire inside the loop, so their messages
  carry the key. Asserted, because a message reading "the outcome is constant" against seven outcomes is
  a message that costs a bisect.
- **`_assert_readable` raises BEFORE the second phase, and the companions prove it matters.** For each of
  S1, S5a and S3a: with the phase-1 raise removed, the frame produces a **bare pandas exception** —
  `KeyError` for S1 and S5a, and a misaligned-index failure for S3a — instead of a `SchemaError`. This is
  the companion form the Stage 7 learning prescribes (*"Test each phase-1 check with a companion that
  removes the phase-1 raise and asserts the bare pandas exception"*), and it is what stops the split from
  being refactored away by someone who sees two helpers where one would do (§4.4a).
- **S8's frame is the one §16 item 3 turns on**, so it asserts both halves: that `SchemaError` is raised,
  **and** that with S8 removed `marginal_odds_ratio` returns the finite `1.3529411765` of §6.2 rather
  than `nan`. That is the measurement the open question should be decided on, and an earlier draft argued
  it from a `nan` the function cannot produce.

### 15.2 The seven [§11] denominators, and the one the workbook witnesses `[data-gated]`

- Six of the seven populations have 92 members and `tici_2b_3` has **91**; the masks are boolean, total,
  and indexed like the cohort frame.
- **The companion is the assertion Stage 8 could not make.** An implementation ranging over `in_model`
  instead of over the outcome mask is constructed and shown to differ: TICI's minority cell moves from
  **6 to 7** and its denominator from 91 to 92, because `y == 1` is `False` on a missing value and the
  absent record is counted as a non-event (§3.3). Stage 8 §14.2 recorded that *"no available frame
  varies it"*; one does now.
- Every mask is a subset of `in_model`, asserted, so §12's a-fortiori claim about Stage 6 §9's rule is
  tested and not merely argued.

### 15.3 The risk difference, and the sign

- `weighted_rd` against a hand computation on fixture 2, to machine precision.
- **The sign convention is uniform across all seven outcomes**: `P(event | bridging) − P(event | EVT
  alone)` on every one, with `higher_is_better` entering no arithmetic. A companion flips one outcome's
  sign per `higher_is_better` and shows that the resulting table has two meanings of "positive" in it
  (§4.1).
- A frame with the arm labels reversed gives the negated risk difference and not the same one — the
  Stage 8 §14.5 orientation move, on the other scale.
- `nan` where an arm carries no positive weight, and **not** `0.0`.

### 15.4 The marginal odds ratio is marginal

- `marginal_odds_ratio` against a hand computation on fixture 2.
- **It differs from the treatment coefficient of a weighted logistic fit on the same rows**, computed and
  shown, because those are the two quantities most easily confused and `model.firth` is one call away in
  this module (§6.1). This is the §8.2-shaped companion: assert the route taken, and show the route not
  taken giving a different number.
- The three returned values are consistent: recomputing the odds ratio from the returned proportions
  reproduces it exactly whenever `or_corrected` is `False`.

### 15.5 The degenerate cell, on all four of §6.2's rows

- `p1 == 0`, `p1 == 1`, `p0 == 0`, `p0 == 1`: the correction fires, the result is finite, and
  `or_corrected` is `True`.
- **`nan` is NOT corrected**: a frame where an arm carries no positive weight returns `nan` for the odds
  ratio with `or_corrected` `False`. A companion shows a wrong implementation — one branching on
  `not (0 < p < 1)`, which `nan` satisfies — manufacturing a finite odds ratio for an arm with no data.
- An interior estimate is **bit-for-bit uncorrected**: the corrected and uncorrected paths are compared
  on a frame where both are computable and shown to agree at `0.000e+00`, so the "only when it fires"
  clause is asserted rather than assumed.
- **The correction's direction is NOT toward the null, and this bullet asserts the opposite of what an
  earlier draft asserted.** The earlier bullet read *"the corrected odds ratio is strictly between 1 and
  the uncorrected limit — i.e. it shrinks toward the null"*, which is false and would have been a test
  that fails once someone built the frame for it. What is asserted: the result is **finite and positive**,
  and `or_corrected` is `True`. Then two companions, both on §15.0.3's frames:
  - **near-even arm weights** (`Sw` 13.626 / 14.110, the workbook's own split): an empty treated event
    cell gives `0.6484`, on the same side of 1 as the uncorrected `0.0`;
  - **the empty cell in the HEAVIER arm** (`Sw` 17.007 / 22.768, `p0 = 0.00647`): the corrected value is
    `1.0201`, **across** 1 from the uncorrected `0.0`. Asserted `> 1.0`, so the crossing is a fact the
    suite records rather than a surprise Stage 10 discovers in a percentile interval (§6.4).
  `[roadmap, amended]`
- **A constant outcome returns a finite number and not `nan`**, which is §6.2's corrected claim: with S8
  bypassed, an all-zero outcome on the population gives `1.3529411765` and an all-one outcome gives
  `0.7391304348`, both with `or_corrected` `True`. The degenerate branch is taken before any division, so
  `0/0` is unreachable. This is the frame §16 item 3 should be decided on.
- **`[data-gated]`: on v7 the branch fires for none of the seven.** Asserted as a property of the
  workbook, and recorded in the coverage map as the reason the branch's coverage comes from fixture 3
  alone.

### 15.6 `predict`

- Against a hand-computed inverse logit on a 3-column design.
- **`predict(fit, X)` on the very `X` the fit was made from is `np.array_equal` to `fit.p`** — an
  EQUALITY and not a tolerance, on both golden-frame fits (§15.0.2). This is the only available oracle
  that `predict` and `firth` agree about the link, and it is the assertion that pins both of §7.3's
  details at once. **Two companions, because it is the equality that makes them fail:**
  - an implementation adding the intercept as a scalar (`beta[0] + X @ beta[1:]`) instead of prepending
    it as a column returns `allclose` but **not** `array_equal` — measured, max difference `1.11e-16` on
    10 of 60 rows. A test written with `allclose` passes on it, which is why this one is not.
  - an implementation omitting the `FIRTH_ETA_CLIP` clip agrees on the golden frame and **diverges from
    `fit.p` above `|eta| = 500`**, shown on a constructed high-`eta` design.
- **`1/(1+exp(-eta))` is asserted safe at the positive tail and UNSAFE at the negative one**, which is
  §3.3's corrected fact: at `eta = 1000` it is `1.0` where `exp(eta)/(1+exp(eta))` is `nan`; at
  `eta = -800` the *unclipped* form emits `RuntimeWarning: overflow encountered in exp` and returns
  exactly `0.0` where `model._probabilities` returns `7.124576406741285e-218`. Both computed alongside,
  so the reason `predict` clips is in the suite and not only in a fence. An earlier draft asserted the
  form was safe at both tails and tested one.
- A design whose columns are **reordered** raises `SchemaError`, and so does one with a column dropped.
  **The companion is the point**: with the check removed, the reordered design returns a full array of
  finite probabilities — a number instead of an error (§7.3).

### 15.7 `m_a(X)`, the counterfactuals, and separation

- `_counterfactuals` returns two arrays differing in exactly the treatment column's contribution:
  `logit(m1) − logit(m0)` is constant across rows and equals the treatment coefficient. This is the
  no-interaction property [§8] prescribes, asserted rather than assumed.
- Every non-treatment column keeps its observed value: a companion that rebuilds the frame from the
  covariate list instead of copying the design produces a **different width** on a frame where a column
  was dropped, and is shown to raise rather than to return a number (§7.3).
- **Separation converges rather than failing.** `separated_frame()` (§15.0.4a) handed to `model.firth`
  returns rather than raising, with `iterations == 6`, `rescales == 0`, `halvings == 0`,
  `exp(β_treatment) == 441.005397` and fitted probabilities strictly interior. This is asserted
  **positively**, exactly as Stage 8 §14.7 asserts its analogue: a test written as "a separated frame
  raises" would pass on a broken fitter and fails on the correct one.
- The unpenalised MLE on the same frame **fails outright**, computed alongside so the penalty's effect
  is a measured contrast rather than a claim — and **all three caps are asserted, not just the raise**:
  `LinAlgError: Singular matrix` from Newton at `maxiter = 200`; `max_abs_beta 67.87` with
  `converged = False` at statsmodels' default of 35, which does **not** raise; and `max_abs_beta
  33.0927` with `converged = True` from BFGS at 60 — a reported success at a value that is a function
  of the cap. The contrast is stronger than the earlier draft's `73.12` and stronger than the raise
  alone: the unpenalised fit does not merely grow, it has no maximum, and the evidence for that is
  three caps giving three answers rather than one exception.
- **No bound is asserted, and §16 records that as the deliberate absence it is** (§9.6).

### 15.7a The `FitError` contract: it propagates and it is never substituted

§9.5's no-substitution rule is the thing this stage argues hardest for and an earlier draft tested it
nowhere. Stage 8 tested the same contract for its own fitter one stage earlier
(`tests/test_outcome.py:673`), which is the shape ported here.

- On a frame where an augmented outcome's `m_a(X)` cannot be fitted, **`secondary` raises
  `model.FitError`** — it does not return a `Secondary` in which that outcome is quietly unaugmented.
  Asserted with `pytest.raises`, and asserted that no `BinaryEstimate` was produced for it.
- **`assert not issubclass(model.FitError, config.SchemaError)`**, which is Stage 8 §11's handover as an
  assertion: `FitError` is the droppable failure and `SchemaError` is not, and Stage 10's whole
  drop-and-count rule is wrong if the two are ever related by inheritance.
- **The companion is what makes it a test.** With `outcome_model`'s call wrapped in `try/except
  model.FitError` and the path downgraded to `"unaugmented"`, `secondary` returns seven estimates, one
  of them silently on a different estimator — and **every other test in §15 still passes**. Asserted, so
  the suite records that the most natural thing an implementer writes when a fit fails on 0.3% of `sich`
  replicates is invisible to everything else here.
- `[data-gated]`: on v7 no outcome raises, which is why the frame is constructed.

### 15.8 The augmentation: the roadmap's two acceptance tests become three, asserted in seven items `[roadmap, amended]`

The heading counts three *acceptance tests* (§22 items 2-4 and §8.4's three consequences) and seven
*assertions*, because they are different things and an earlier heading conflated them — it read "the
three tests the roadmap's two become" over a list of five, and the coverage map said "5 assertions, 4
companions" (§22.3 item 17). Test 1 is items 1-2, test 2 is items 3 and 6, and the third acceptance test
is §15.10. Items 4, 5 and 7 are this document's own.

1. **Two different constants.** `m₁ ≡ 0.62`, `m₀ ≡ 0.23`: `augmented_rd` equals `weighted_rd` to
   `1.943e-16` on the golden frame. **And the companion is what makes it a test**: an implementation
   dropping the correction term returns `rd − 0.390000`, exactly `c₀ − c₁`.
2. **One constant, and it is kept as a companion and asserted to be insufficient.** `m₁ = m₀ = 0.45`:
   the identity holds at `8.327e-17` — and the dropped-term implementation **also** passes, at
   `+0.000000`. Asserted as a *pair*, so the suite records that the roadmap's phrasing admits a form the
   test cannot fail.
3. **The tilt.** On `tilt_frame()` (§15.0.4), `augmented_rd` with `h` is asserted; the same call with `w`
   is computed and shown to differ by `1.030664e-02`, against a bound of `1e-03`. Neither the
   constant-model tests nor any tolerance on them can distinguish these two (§8.4), so this is the test
   that pins `h` **as an argument**.
4. The correction term's denominator is the **whole** [§11] population and not an arm's: a companion
   normalising it by `Σ₁w` is computed and shown to differ. **Recorded as redundant and kept anyway**:
   the arm-normalised implementation already fails item 1 decisively, so this asserts nothing item 1
   does not. **The SIZE of that failure is a property of the frame and the `+0.482` an earlier draft
   carried is not `golden_frame()`'s** — that number was measured on a 92-record ATO-shaped
   construction, and on the frame item 1 is actually specified against the deviation is
   `+0.0812746102` (§22.4 item 4). Both are decisive; the assertion is written against **item 1's own
   `1e-15` tolerance** rather than against either measured gap, because that is the quantity which
   decides whether item 1 catches it. It stays because it localises the failure — item 1 says "the
   augmentation is wrong", item 4 says "the correction term's denominator is wrong" — and §8.1's
   docstring no longer claims item 3 is what catches it (§22.3 item 13).
5. `h` is passed in and never recomputed inside `augmented_rd`: asserted by giving it an `h` that is not
   `e(1−e)` and checking the result changes, which is what makes it the caller's value.
6. **AND THE CALL SITE IS PINNED, WHICH IS WHAT ITEMS 3 AND 5 DO NOT DO.** Items 1-5 all call
   `augmented_rd` directly with arrays the test constructed, so they prove the *function* distinguishes
   `h` from `w` and say nothing about what `secondary` passes it. On the golden frame:
   `secondary(...).estimates[k].augmented` **equals** `augmented_rd(y, a, w, _tilt(e), m1, m0)` exactly,
   and **differs** from `augmented_rd(y, a, w, w, m1, m0)`. Without this, an implementation typing `w` on
   §11's one tilt line ships the wrong estimand with items 1-5 green — which is what the Failure modes
   table calls *"the failure the spec exists to make catchable"* (§8.2a, §22.3 item 2).
7. `_tilt(e)` is `e(1 − e)` elementwise, against a hand computation, and `_tilt` is what §11 calls —
   asserted by scan, so the expression cannot migrate back inline.

### 15.9 The rare-minority guard `[roadmap, amended]`

- `_augmentable` is `True` for a minority cell at the threshold and `False` one below it, asserted
  against `C.RARE_MINORITY_THRESHOLD` and never against `10`.
- **An override is a route IN.** `tici_2b_3` with a minority cell of 6 is augmentable; the same cell on
  an outcome with no override is not. This is the roadmap-versus-SAP contradiction of §7.2 asserted in
  the direction the amendment settles.
- `_augmented_path` returns the three declared strings, and **the block is `augmented`, `covariates`,
  `fit` and `dropped`** — `None`/`None`/`None`/`()` exactly when the path is `"unaugmented"`, all four
  populated otherwise, asserted across all seven outcomes (§3.1).
- **`reduced` is asserted SEPARATELY and is not in the block.** It is `True` exactly when the key is in
  `OUTCOME_MODEL_OVERRIDES` — so on v7, `True` for `tici_2b_3` and `False` for the four `full`-path
  outcomes *whose other four fields are populated*. A test written from the earlier draft's docstring,
  which listed `reduced` in the block, fails on `mrs_0_2_90d`, `mrs_0_1_90d`, `death_90d` and
  `mrs_5_6_90d`; that this is what would have happened is asserted by a companion (§22.3 item 4).
- **`augmented_path == "reduced"` and `reduced == True` agree on all seven**, which is the invariant that
  makes the two fields non-contradictory even though only one is in the block.
- **A frozen path overrides a crossed threshold** (§11.1): on a frame where `ph2`'s minority cell has
  crossed to 12, `secondary(..., paths={... "ph2": "unaugmented" ...})` leaves it unaugmented, while the
  same frame with `paths=None` augments it. This is §9.5's 46.1% as an assertion rather than a
  measurement. A partial `paths` map raises `SchemaError`.
- **The minority cell is counted on the masked population.** On a frame where the mask and `in_model`
  differ, the count from the mask is asserted and the count from `in_model` is computed and shown to
  differ (§9.1) — the same construction §15.2 uses, reused so the two facts are checked on one frame.
- `outcome_model_covariates` returns the shared [§6] list for six outcomes and the two-name override for
  `tici_2b_3`, and raises `KeyError` for an unregistered outcome. First read of a function written at
  Stage 8 (§13).
- **TICI's reduced design is 5 parameters, and the count's contingency is asserted with it.** The fitted
  design is 4 columns plus the intercept, matching the amendment's *"five parameters"*. A companion
  builds a frame in which USZ contributes one record and shows the same declared model becoming **6
  parameters against 6 non-events** — saturated — so §7.2's fragility fails a test rather than living in
  a paragraph. `[data-gated]` for the 5; the companion is synthetic.
- **`[data-gated]`:** the seven paths on v7 are four `full`, one `reduced`, two `unaugmented`, and
  `sich` and `ph2` are the two — which is the SAP amendment's own closing prediction (§9.3).

### 15.10 Not doubly robust `[roadmap, amended]`

On `dr_population` (§15.0.5), at `n = 200,000`, `DR_SEED = 20260825`:

- **Heterogeneous arm**: the augmented estimate with an oracle `m_a(X)` and a misspecified `e(X)` misses
  `ATO(e_true)` by more than `0.015`. Measured across `DR_SEEDS_12`: `[0.022950, 0.030249]`.
- **And it recovers what §8.3 predicts it recovers**: the same estimate is within `0.010` of the ATO
  indexed by the *misspecified* score — measured `0.005428` at `DR_SEED`. This is the positive half, and
  it is what distinguishes "not doubly robust" from "wrong". §8.5 explains why it is close by
  construction and not by luck: with an oracle `m_a` the correction term *is* `ATO(e_mis)` exactly.
- **Constant-risk-difference companion**: the same misspecified `e(X)` recovers the estimand to within
  `0.010`. Measured across `DR_SEEDS_12`: `[0.000176, 0.004929]`.
- The two bounds sit in a **band that is empty across all twelve seeds** — nothing between `0.004929`
  and `0.022950`, and both `0.010` and `0.015` verified inside it — which is Stage 8 §6.3's empty-band
  calibration applied to a tolerance (§8.5).
- **The no-interaction logistic companion is asserted to be biased**, at `|t| = 21.3` over 200
  replicates at `n = 20,000`, so the suite records that the obvious "homogeneous" construction is not
  homogeneous on the risk-difference scale and would not have shown what the roadmap says it shows
  (§8.5). The constant-risk-difference arm at the same size is `|t| = 0.9`.
- **A note in the test, not an assertion**: the same measurement at `n = 92` gives a mean error of
  `−0.025665` with `sd 0.087669` and the wrong sign in `40.5%` of replicates. The test does not run
  there, and the docstring says why.
- **No claim is made about intermediate `n`.** An earlier draft carried a four-row sizing table and drew
  a conclusion from its `n = 20,000` row that its own numbers contradicted (§8.5, §22.3 item 3).

### 15.11 The module boundary

- `outcome.py` does not import `balance`, and `secondary` does not read `Primary` — by scan of the
  module source, for Stage 7 §14.11's reason: not because it would fail, but because it would succeed.
- `ps.e` and `ps.w` are read in exactly **two** places, both `.loc[in_estimate]`, and no read of either
  ranges over `notna()` or over `in_model` (§12). By scan.
- No `nan`-**filling** anywhere: `fillna`, `.fillna(0`, `interpolate`, `ffill` and `bfill` appear
  nowhere in the added code.

  **`dropna` is NOT on that list, and an earlier draft of this bullet put it there — against this
  document's own §4.4b fence.** S5b is
  `sorted(set(arm.dropna()[~arm.dropna().isin(...)]))`, and the `dropna` is not filling anything: S5b
  asks the treatment column **two** questions — is any value missing, and is any value outside the
  declared arm codes — and the second cannot be asked of a `nan`, because `isin` against one is
  `False` and the record would be reported as an out-of-range arm rather than as an absent one. The
  missing count is already reported by the branch immediately above it. This is Stage 6's D3/D4 split
  exactly (`model.py`): *"D3 excludes missing values from its own check because they are D4's, not
  because they are harmless."*

  Since DoD-20 makes the fences authoritative and the fence is right, the **bullet** was the defect.
  §15.11 therefore asserts the filling verbs are absent **and** that the one `dropna` is located in
  `_assert_secondary_inputs` — a location rather than a waiver, so a second one appearing anywhere
  else fails. §22.4 item 1.

### 15.12 Stage 9 adds nothing, refits nothing, and breaks nothing earlier

- The Stage 1-8 suites pass **unedited** except for the ledger count (§10.1) and the five `test_config`
  additions (§13). Asserted by running them.
- `weighted_proportion`, `cumulative_rd` and `primary` are byte-identical to Stage 8's, and `primary`
  returns the same `Primary` before and after `secondary` runs on the same frame — so the two are
  independent and order does not matter.
- `model.firth` and `model.design` are unchanged; the Stage 6 suite is the assertion.

### 15.13 "model-assisted", never "doubly robust" `[roadmap]`

A scan for `doubly robust`, `doubly-robust`, `AIPW` and `DR estimator`, case-insensitive, asserting none
appears; and for `model-assisted`, asserting it does. This is a test of a word, which no behavioural test
can be — the failure mode is a docstring, and [§8]'s reason for the prohibition is that the word licenses
a conclusion the estimator does not support (§8.3).

**The scope is every Stage 9 deliverable that ships text, which is wider than an earlier draft's.** That
draft scanned `outcome.py`, `model.py`, `config.py` and the rendered audit text, and omitted
`tests/fixtures_stage9.py`, the three amended test modules and `tests/reference/aug_psweight.R` — all
Stage 9 deliverables, and §19b's own prose discusses the R script in AIPW terms, so the omission was not
hypothetical (§22.3 item 18):

```
  scanned   outcome.py, model.py, config.py
            tests/fixtures_stage9.py, tests/test_outcome.py, tests/test_model.py, tests/test_config.py
            tests/reference/aug_psweight.R
            the RENDERED audit text of all three §10 entries
  asserted  none of the four forbidden phrases appears; `model-assisted` does
```

**And `specs/` is scanned with one declared exemption.** §8.3 bans the phrase *"in code, docstrings,
audit output **and this document**"*, and this document uses `AIPW` and `double-robustness` descriptively
about a dozen times — §15.0.5's population, T10's title, §19's pilot survey, §19b's argument about why no
package's routine is an oracle. Those are the *name of the thing the estimator is not*, which is the one
use the ban cannot forbid without making itself unstatable. The exemption is therefore: **a line that
names the prohibition, quotes it, or explains what the estimator is not, is exempt; a line that describes
this estimator is not.** Asserted mechanically by exempting lines containing `never`, `not`, `NOT`,
`forbidden` or a `§8.3` citation, and §22.3 item 18 records that the exemption is a judgement rather than
a rule a scanner discovers.

### 15.14 The three audit entries, and the log

- Three entries, kind `model`, steps `binary_estimable`, `binary_estimates`, `binary_outcome_models`, in
  that order, at positions captured before the call.
- `len(data.KINDS) == 9` and the ledger is **30** (§10.1).
- `binary_outcome_models` has **five** rows on v7 and not seven `[data-gated]`.
- `case_ids` is empty on all three, and the record `tici_2b_3` loses is **not** named (§10.2).
- `tau` renders as `--` and never `0.0` for an unaugmented outcome; `corrected` renders even when false
  for all seven; a reduced specification renders as its covariate names and not as `True` (§10.3).
- **The log is byte-identical across `PYTHONHASHSEED=0` and `1`**, which is where the seven-outcome
  dict's iteration order is actually asserted.
- **STRUCTURAL, PER ENTRY, ON THE RENDERED TEXT: header cells == separator cells == every body row's
  cells, and no cell contains a `|`.** `data._md_table` (`data.py:242-245`) does no escaping and sizes
  its separator from `len(rows[0])`, so a pipe inside any cell silently adds markdown columns to that row
  alone. Asserted on the rendered string and **not** via the hash, because a table that is identically
  broken under both seeds passes the byte-identity check above — which is exactly how Stage 7 shipped
  `|SMD| < 0.1` and `worst |SMD|` rendering 8/6 and 15/13. The assertion covers all three Stage 9
  entries, so §10.4's `max_abs_beta` is checked rather than trusted (§22.3 item 7).
- **`proportion` cannot be mutated through the returned estimate.** `BinaryEstimate` is `frozen=True`,
  which stops attribute rebinding and does nothing about a `dict` field: asserted that
  `estimate.proportion[1] = 99.0` leaves `estimate.rd` and the rendered table unchanged, because §11
  stores `dict(share)` rather than the object `weighted_proportion` returned. A companion storing the
  object directly is shown to let a caller change a frozen estimate — and, before the copy, to change
  the dict `marginal_odds_ratio` returned to a different caller (§22.3 item 14).

### Coverage map

```
  outcome.py
    weighted_rd ................................... 15.3        hand computation + sign + nan
    marginal_odds_ratio
      interior branch ............................. 15.4, 15.5  bit-for-bit uncorrected
      degenerate branch (4 routes) ................ 15.5        fixture 3 ONLY — see below
      the null-CROSSING route ..................... 15.5        empty cell in the HEAVIER arm
      constant-outcome route (finite, not nan) .... 15.5, 15.1a §6.2's corrected claim
      nan branch .................................. 15.5        + the companion that corrupts it
    outcome_model ................................. 15.7, 15.9  both covariate lists
      FitError propagates, never substituted ...... 15.7a       + the try/except companion
    _counterfactuals .............................. 15.7        + the rebuild companion
    _tilt ......................................... 15.8 item 7 hand computation + the scan
    augmented_rd .................................. 15.8        7 assertions, 5 companions
      the zero-total guard (3 routes) ............. 15.8        treated / control / tilt
    _augmentable / _augmented_path ................ 15.9        threshold, override, block-move
    _assert_readable          S1, S5a, S3a ........ 15.1a       one frame per branch + 3 companions
    _assert_secondary_inputs  S2, S3b, S4, S5b .... 15.1a       one frame per branch + collection
    _assert_estimable         S6-S8 ............... 15.1a       one frame per branch, outcome named
    the three _record_* / _*_table ................ 15.14       + the structural pipe/cell assertion
    _estimable_detail ............................. 15.14
    secondary (ordering, 6 constraints) ........... 15.1, 15.12
      the frozen-paths parameter (3 routes) ....... 15.9        None / complete / partial-raises
  model.py
    predict ....................................... 15.6        5 routes + 2 oracle companions

  THE MAP IS THE COUNT. An earlier draft closed with "= 71 / 71 routes with a test that can fail",
  which was arithmetic over a list that included `_or_detail` — a function no section defined and
  nothing called — and that pointed S1-S8 at §15.1, which specified S1 alone. A total not derived
  from the rows above cannot be checked against them, so there is no total. Every route listed has
  a test that can fail; a route not listed is not covered.

  5 unreachable ON THE WORKBOOK, covered by constructed frames alone: the four degenerate-cell
    routes (§6.4 — the branch fires for none of the seven on v7) and the nan branch (S7 makes it
    unreachable through `secondary`; it is reachable through `marginal_odds_ratio` directly).
  1 unreachable on the workbook, covered by a constructed frame: the FitError route (§15.7a).
  8 unreachable on the workbook BY DESIGN: S1-S8, which exist for [§10]'s frames (§4.4) and which
    §15.1a is the only thing that reaches.
  0 unreachable-by-construction.
```

### Failure modes

| Failure | Test | Error handling | Visible or silent |
|---|---|---|---|
| Ranges over `in_model` instead of the outcome mask | 15.2 | none — wrong number | **Silent on six of seven**, and on `tici_2b_3` it moves the denominator 91→92 and the minority cell 6→7. Nothing raises |
| Reads `Fit.p` as `m₁` or `m₀` | 15.7 | none | **Silent.** Wrong on every row assigned to the other arm; the result is a plausible number |
| Augmentation tilted by `w` instead of `h` **inside `augmented_rd`** | 15.8 item 3 | none | **Silent, and both roadmap tests pass** (§8.4). This is the failure the spec exists to make catchable |
| Augmentation tilted by `w` **at the call site** — `secondary` passing `w` where §11 passes `_tilt(e)` | 15.8 item 6 | none | **Silent, and items 1-5 all pass** — they call the estimator directly with their own arrays. Untested until round 3 (§8.2a) |
| `FitError` caught and the outcome downgraded to unaugmented | 15.7a | none | **Silent, and every other test in §15 passes.** Mixes two estimators in one interval — §9.5's prohibition reached by the most natural `try/except` an implementer writes |
| A `\|` in any audit table cell | 15.14 | none — the table renders crooked | **Silent, and byte-identity across hash seeds does NOT catch it**: identically broken twice is still identical. Stage 7 shipped it twice |
| `predict` adds the intercept as a scalar, or omits the clip | 15.6 | none | **Silent on ordinary designs** — `allclose` to `fit.p` but not equal, and divergent above \|eta\| = 500. Only an EQUALITY oracle fails on it |
| A caller mutates `estimate.proportion` | 15.14 | none | Silent; `frozen=True` does not freeze a `dict` field, and before §11's copy the mutation also reached `marginal_odds_ratio`'s other caller |
| A missing outcome column reaching S2 | 15.1a | bare `KeyError`, not S1's `SchemaError` | Visible but **wrong** — a pandas traceback in place of the message this document wrote. Stage 7's shipped defect (§4.4a) |
| A partial `paths` map | 15.9 | `SchemaError` | Visible. Without the guard it would be **silent**, and the two-estimator mixture would arrive one outcome at a time (§11.1) |
| Correction term normalised by an arm total | 15.8 item 4 | none | Silent; passes the constant-model tests |
| Correction term dropped entirely | 15.8 items 1-2 | none | Silent **iff** the two constants are equal |
| `predict` on a reordered design | 15.6 | `SchemaError` | Would be **silent** without the by-name check: returns finite probabilities |
| Eligibility recomputed per replicate | Stage 10's | none here | **Silent**, and mixes two estimators in 46.1% of `ph2`'s replicates (§9.5) |
| `nan` proportion "corrected" to a finite OR | 15.5 | none | Silent; manufactures an odds ratio for an arm with no data |
| Separated `m_a(X)` | 15.7 | **none, by decision** | **Silent by design** (§9.6). Mitigated only by `max|β|` in §10.4's table and by Stage 10 reporting its distribution |
| `propensity.fit` on a resample | §12.3 | `SchemaError` | Visible, loudly, on 26.0% of replicates — and it is the wrong exception for Stage 10 to see |
| An outcome in `OUTCOMES` but not in the frame | 15.1 | `SchemaError` (S1) | Visible |
| A third `family` value added to the registry | 15.1, §13 | `test_config` assertion | Would be **silent** in `by_family` and would get its own Benjamini–Hochberg group |
| `OR_CONTINUITY` changed | §13 | `test_config` bound | Visible: asserted against both endpoints |

---

## 16. Known gaps carried forward

1. **`propensity.fit` raises on 26.0% of stratified replicates** (§12.3). The blocker. Not fixed here
   because the fix is a Stage 6 amendment justified by a Stage 10 requirement that is not yet
   specified. `TODOS.md`; trigger: the first line of Stage 10.
2. **Whether `sich`'s `max|β|` tail leaves its interval usable** (§9.6). No bound is prescribed, and the
   measurement that would settle whether the 15.1%-above-14 tail contaminates the interval requires a
   synthetic construction reproducing the cohort's centre structure — which the one used here does not
   (its `max|β|` reached 11.475 against the workbook's 80.98). `TODOS.md`; trigger: Stage 10 reporting
   the `max|β|` distribution, which is when the question becomes answerable from replicates already
   drawn.
3. **S8 may want to be `FitError` rather than `SchemaError`** (§4.4). A constant outcome on a replicate is
   arguably a sparse replicate and not a bug, and `sich` is where it would happen. **The reason recorded
   in an earlier draft was false and the true reason is stronger.** That draft argued S8 must be
   `SchemaError` because *"§6.2's `0/0` gives `nan` silently"* — but the degenerate branch fires before
   any division, so `0/0` is unreachable and a constant outcome returns a **finite** odds ratio:
   measured, `1.3529411765` all-zero and `0.7391304348` all-one (§6.2). So softening S8 does not trade a
   raise for a `nan` a reader would notice; it trades a raise for **a plausible number near 1 on an
   outcome nobody observed an event in**. That makes the case for `SchemaError` stronger than the case
   originally written for it, and it is still not decisive, because a sparse replicate is a real thing
   and dropping it is a real option. **This is the sharpest open question the stage leaves**, and §15.1a
   asserts both halves so the decision is made on measurements rather than on the false premise.
   `TODOS.md`; trigger: the first Stage 10 replicate that hits it.
4. **TICI's reduced model is five parameters only while USZ contributes no records** (§7.2). A single
   USZ record makes it six against six non-events. §15.9's companion tests the contingency; nothing
   prevents it. `TODOS.md`; trigger: any workbook revision that adds a USZ record.
5. **No pinned regression number on any reported estimate exists in git** (§4.3). Stage 8's cost, and
   Stage 9's identically. The golden vector is synthetic. `TODOS.md`, on Stage 8's existing item.
6. **The augmentation's variance cost is measured only where it can only hurt** (§9.6). The RMSE table is
   a lower bound on the cost, from a constant-risk-difference truth. What it does not measure is the
   case [§8] hopes for — an informative outcome model — and that case cannot be constructed without
   assuming the thing being estimated.
7. **`by_family` assumes two families** (§13). A third registry family would get its own
   Benjamini–Hochberg group silently. `test_config` asserts two; that is the only guard. **No `TODOS.md`
   item**: the guard is adequate and the trigger is a registry edit that `test_config` fails on.
8. **The SAP's reduction is conditional and `_augmentable`'s override is not** (§9.2). The amendment reads
   *"**Where** the minority cell **cannot support** the full covariate set but can support a smaller
   one..."*; `_augmentable` returns `True` for a key in `OUTCOME_MODEL_OVERRIDES` **whatever its minority
   cell is**, so an override fixes the estimator in both directions. §9.5 presents that as a design virtue
   — *"naming an outcome fixes its estimator"* — and it is one, for the 46.1% reason; but it also means
   that if a workbook revision took `tici_2b_3`'s minority cell to 40, it would still be fitted with the
   two-name reduced model rather than the [§6] set. Moot on v7 at a minority cell of 6, and §7.2 spent
   forty lines on a roadmap/SAP conflict about TICI while this one went unnamed. **No `TODOS.md` item**:
   it is contingent on a workbook revision that item 4 already carries a trigger for.
9. **`_tilt` is `outcome.py`'s and arguably `Propensity`'s** (§8.2a). One consumer, so one home here. If a
   second consumer appears — Stage 10 re-deriving `h` for any reason — moving it onto `Propensity` is the
   right change. **No `TODOS.md` item**; trigger: a second reader of `h`.

**Which of these reach `TODOS.md`, stated once so DoD-13 has something to check.** Items 1, 2, 3 and 4
are new `TODOS.md` entries. Item 5 is an amendment to Stage 8's **existing** pinned-regression item, not
a new one. Items 6, 7, 8 and 9 are carried here only, each with the reason above. **Four new items and
one amendment** — an earlier DoD-13 required *"§16's seven items"* against a §16 that then had seven
entries of which four were destined for the file (§22.3 item 5).

## 17. What Stage 9 deliberately does not decide

- **Whether the degenerate odds ratio should be `nan` instead of corrected.** §6.4 is the argument both
  ways. The PI-reversible alternative is: return `nan`, report the risk difference alone, name the empty
  cell. Consequence if reversed: Stage 10 reports a count of undefined replicates in place of an
  interval for `sich` (17.0% of replicates) and `tici_2b_3` (13.0%). **PI decision, [§8] amendment.**

  **And round 3 adds a fact that bears on this decision and was not available when it was made.** The
  correction does not only widen the estimate — because the pseudo-counts are weight-sums and the arms'
  totals differ, **it can move the odds ratio across the null**. Measured: an empty treated event cell
  with `Sw` 17.007 against 22.768 returns `1.0201` where the uncorrected value is `0.0`, and about
  **0.9% of empty-cell replicates cross**, i.e. roughly **3 replicates in an `N_BOOT` = 2000 bootstrap in
  which a safety outcome with no bridging events at all is reported as favouring EVT alone** (§6.4). The
  `nan` alternative cannot do this. So the trade is no longer "finite but wide" against "undefined": it
  is "finite, wide, and occasionally sign-inverted on the safety family" against "undefined". That is a
  different decision and it is the PI's.
- **`OR_CONTINUITY`'s value.** `0.5` is the conventional constant against row counts; §6.4 records that
  it is **2.9× more aggressive in the treated arm and 3.8× in the control arm** against these
  weight-sums — per arm, because the asymmetry between those two numbers is what permits the crossing
  above. Scaling it to `Σw` was rejected for a stated reason and not overlooked, and that rejection is
  what leaves the asymmetry in place. **PI decision.**
- **`RARE_MINORITY_THRESHOLD`'s value**, with `ph2` measured one below it (§9.4). Prespecified before any
  outcome was examined by arm; looked at and left alone.
- **Whether `higher_is_better` should orient the risk differences.** §4.1 keeps the contrast's sign
  uniform and leaves direction to Stage 14. The alternative makes a seven-row table in which "positive"
  varies by row.
- **Whether the augmented or the unaugmented estimate is the headline for an augmented outcome.** [§8]
  says both are reported and says which is *"model-assisted"*; it does not rank them. Stage 14's.

## 18. NOT in scope for Stage 9

| Considered | Why deferred |
|---|---|
| Any interval, standard error or p-value | [§10]. Stage 8 §5.6's position, unchanged: the weights are a tilting function of an *estimated* score, so an inverse-information standard error omits a term [§7] says enters the interval |
| Benjamini–Hochberg within families | Stage 11 [§13]. `by_family` is the handover, not the correction |
| Cross-fitting or sample-splitting for `m_a(X)` | [§8] forbids it by name: *"cross-fitting folds would be too small to serve their purpose"* |
| Interactions or splines in `m_a(X)` | [§8]: *"linear terms... no interactions, no splines"* |
| A bound on the nuisance coefficient | §9.6. No empty band exists to calibrate one from, and the predictions are bounded regardless |
| Augmenting the ordinal outcome | [§8]: *"There is no simple augmented form of a proportional-odds fit"* |
| An E-value for any binary estimate | [§6] specifies it for the primary only |
| Subgroup estimates of these outcomes | Stage 11 [§13], labelled hypothesis-generating |

## 19. What already exists, and what to lift

| From pilots | Status |
|---|---|
| `pilots/analysis.py`'s binary contrasts | **Not lifted.** They are unweighted crude differences; [§8] requires weighted marginal estimates in the overlap population |
| Any augmented / AIPW implementation | **None exists in the pilots.** The augmentation is new code, which is why §15.8 carries five assertions and four companions rather than a comparison |
| `model.firth` | Lifted whole, unmodified, from Stage 6 |
| `model.design` | Lifted whole, unmodified |
| `weighted_proportion` | Lifted whole from Stage 8, which wrote it for this stage |

### 19b. Validation against a reference implementation

`tests/reference/aug_psweight.R`, behind the gate Stage 7 built: read a CSV, ask the library, write it
back. Two rows are checked and the second is the one that matters.

- **The weighted marginal proportions and risk difference** against `PSweight`'s overlap-weight output on
  the same frame. This checks the weighting, and it is the same oracle Stage 6 §18b used one function
  over.
- **The augmented estimate is checked against a hand-assembled R computation and NOT against a package's
  AIPW routine**, because every maintained one targets the ATE or the ATT and therefore tilts by
  `1/e(1−e)`-style weights or by `e/(1−e)`, not by `h = e(1−e)`. An oracle computing a different
  estimand is not an oracle. The R script assembles §8.1's formula term by term, which makes it a
  cross-language check on the *arithmetic* and explicitly not on the *estimand* — and §8.4's third test
  is what checks the estimand, in Python, because only this repository knows which tilt [§8] means.

**What 19b does not establish**, stated because the distinction is easy to lose: agreement with a
hand-assembled formula in another language shows the two transcriptions agree. It does not show the
formula is [§8]'s. §8.4 exists for that.

## 20. Implementation tasks

- [x] **T0 (P1)** `tests/fixtures_stage9.py`: §15.0's constructions and four constants — seven
      functions, enumerated in §13. **Gates
      T3, T4, T5, T7, T9, T10 and T11** — every pin in this document is measured against these, so
      nothing that asserts a number can be written before it.
- [x] **T1 (P1)** `config.py`: `BINARY_OUTCOMES`, `OR_CONTINUITY`; `test_config.py`'s five assertions.
      No other file compiles against these until it lands.
- [x] **T2 (P1)** `model.py`: `predict` (§7.3's fence, prepended intercept and `FIRTH_ETA_CLIP`), the
      docstring sentence; §15.6's tests including the reordered-design companion and the two oracle
      companions.
- [x] **T3 (P1)** `outcome.py`: `weighted_rd`, `marginal_odds_ratio`, the two dataclasses; §15.3-§15.5.
      Depends on T1 (`marginal_odds_ratio` reads `C.OR_CONTINUITY`; `Secondary.by_family` reads
      `C.BINARY_OUTCOMES`) and on T0 (§15.5's frames).
- [x] **T4 (P1)** `outcome.py`: `outcome_model`, `_counterfactuals`; §15.7. Depends on T2 and T0.
- [x] **T5 (P1)** `outcome.py`: `_tilt`, `augmented_rd` with its zero-total guard; §15.8's seven
      assertions and five companions. Depends on T0 (`tilt_frame`, `golden_frame`) and, for items 3, 6
      and 7, on T4.
- [x] **T6 (P1)** `outcome.py`: `_augmentable`, `_augmented_path`; §15.9. Depends on T1.
- [x] **T7 (P1)** `outcome.py`: `_assert_readable` and `_assert_secondary_inputs` (§4.4's two frame-level
      phases) and `_assert_estimable`; **§15.1a**, one frame per branch, with the three phase-1
      companions. Depends on T1 (S1 and S2 range over `C.BINARY_OUTCOMES`) and T0.
- [x] **T8 (P1)** `outcome.py`: the three audit entries, their tables and `_estimable_detail`; §15.14
      including the structural pipe/cell assertion. Depends on T3 (the tables take
      `dict[str, BinaryEstimate]`).
- [x] **T9 (P1)** `outcome.py`: `secondary` with its `paths` parameter, wiring T3-T8 in §11's order;
      §15.2, §15.7a, §15.9's frozen-path routes, §15.12.
- [x] **T10 (P2)** §15.10's double-robustness tests. Depends on **T0 and T5** — it asserts properties of
      `augmented_rd` on `dr_population`. An earlier draft called it "independent of everything above".
- [x] **T11 (P2)** §19b's R oracle. Depends on **T5** for the Python side of the comparison, and on T0
      for the frame it writes out. Stage 7's `_r_environment` already searches candidate library
      directories when `R_LIBS_USER` is unset, so the gate needs no work here.
- [x] **T12 (P1)** the roadmap amendments (§22) and `TODOS.md` (§16's four items plus the amendment to
      Stage 8's existing one). Lands in the same commit.
- [x] **T13 (P2)** §15.13's word scan, over the shipped modules **and** the amended test modules and
      `tests/reference/aug_psweight.R` (§15.13).

**What can be built in parallel, and what cannot. The earlier draft's graph was wrong in four places and
this one is measured.** `T0` gates everything that asserts a number. `T1` gates **T3, T6, T7 and T9** —
not just T6 and T9: verified by execution, `marginal_odds_ratio` raises `AttributeError: module 'config'
has no attribute 'OR_CONTINUITY'` without T1, and `Secondary.by_family` raises the same for
`BINARY_OUTCOMES`. `T2` gates `T4`, which gates `T5`'s non-constant-model tests — though T5's
*constant*-model tests need neither, since they pass arrays directly. `T3` gates **T8**, whose tables take
`BinaryEstimate`. `T5` gates **T10 and T11**, which the earlier draft called independent of the module
entirely.

Lanes: `T0 → T1` first and alone. Then `{T2 → T4}`, `{T3 → T8}`, `{T6}`, `{T7}` in parallel. Then `T5`.
Then `{T10}`, `{T11}`, `{T13}` in parallel with `T9`. `T12` any time.

The one serialisation that matters: **T9 last**, because `secondary`'s six ordering constraints (§11) can
only be tested once the pieces exist, and writing it first produces a function whose order is asserted by
nothing. **And §22.2's closing advice is withdrawn**: it recommended *"write §15.1 and §15.14 early —
before the estimators"*, which cannot be done, because both call `secondary` and `secondary` is T9. What
it was reaching for is served by T0 instead — the fixtures are the thing nothing had exercised.

### Definition of done

Complete when all of the following hold, and not before.

1. `uv run pytest -v` green, with `data/` present.
2. `uv run pytest -v` green with `data/` absent — every data-gated test skipped, none failing.
3. The Stage 1-8 suites are **unedited** except the ledger count and `test_config`'s five additions, and
   the raw-name scan in `test_config.py:450` still passes.
4. The audit log is **byte-identical** across two hash seeds. Spelled out as a heredoc, as Stages 7 and
   8 spell theirs out, because an earlier draft ran `import runner` and **no `runner.py` exists** —
   `find . -name "runner*"` outside `.venv` returns nothing (§22.3 item 16):
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
   outcome.secondary(df, ps, audit)
   sys.stdout.write(audit.to_markdown())'
   PYTHONHASHSEED=0 uv run python -c "$S" > /tmp/s9-a.md
   PYTHONHASHSEED=1 uv run python -c "$S" > /tmp/s9-b.md
   diff /tmp/s9-a.md /tmp/s9-b.md && echo IDENTICAL
   ```

   **This is Stage 8 §15's driver with one line added**, deliberately rather than rewritten: `data.load`
   returns `(DataFrame, Audit)` and there is no `data.audit()`, the pipeline threads one `Audit` through
   every step, and the log goes to **stdout** redirected by the shell rather than through
   `Audit.write`, which derives its own path from the source label and would put both runs in the same
   file. An earlier draft of this item got all three wrong while fixing a fourth defect — it had called
   a nonexistent `runner` module (§22.3 items 16 and 25).
5. **The tilt test has been seen to fail** when `h` is replaced by `w` in `augmented_rd`, **and the
   call-site test has been seen to fail** when `_tilt(e)` is replaced by `w` in §11 — two assertions,
   because the first passes on the second bug (§15.8 items 3 and 6).
6. **The constant-model test has been seen to pass** on an implementation with the correction term
   deleted, when the two constants are made equal — so the reason §15.8 uses different ones is observed
   and not taken on trust.
7. **The `tici_2b_3` denominator test has been seen to fail** on an implementation ranging over
   `in_model`, giving minority 7 and denominator 92.
8. **`predict`'s reorder check has been seen to return a plausible array** with the check removed; and
   **`predict(fit, X) == fit.p` has been seen to fail** on the scalar-intercept form, at `allclose` but
   not `array_equal`.
9. **A separated frame has been seen to converge** through `model.firth` on `separated_frame()` and
   return `exp(β_treatment) = 441.005397` in 6 iterations with all three safeguard counters at zero,
   while the unpenalised MLE on the same frame raises `LinAlgError`.
10. `grep -ri "doubly.robust\|AIPW" outcome.py model.py config.py tests/ ` returns nothing outside the
    passages that quote the prohibition, and the rendered audit text contains `model-assisted` (§15.13).
11. The five augmented outcomes on v7 are the four full plus `tici_2b_3`, and `sich` and `ph2` are
    unaugmented — the SAP amendment's own prediction, observed.
12. `max_abs_beta` appears in the rendered `binary_outcome_models` table for every augmented outcome,
    **and the structural assertion of §15.14 passes on all three entries** — equal cell counts, no pipe
    in any cell. Item 12 without the second clause passes on a crooked table (§22.3 item 7).
13. **§16 items 1-4 are four new `TODOS.md` entries, each with its trigger; §16 item 5 is an amendment to
    Stage 8's existing pinned-regression item; §16 items 6-9 are carried in §16 only.** Verified by
    counting entries added, not by counting §16's paragraphs.
14. The roadmap's Stage 9 section carries its `**Spec:**` line and §22's corrections.
15. This document's §21 has no row reading `stated` that a run could have produced.
16. **Every fixture this document pins a number against is in `tests/fixtures_stage9.py`, and every pin
    reproduces from it.** Checked by re-running §21b's pass 3. This is the criterion an earlier draft
    could not have met: three of its fixtures existed only in a scratchpad (§22.3 item 1).
17. **Every count stated in this document has been re-derived before the commit** — the privates of §3.2,
    the audit ledger of §10.1, the `KINDS` length, and the `TODOS.md` entries of item 13. Three earlier
    drafts disagreed with themselves about the first of those (§22.3 item 5).
18. **`§15.1a` passes, all eight branches**, and each of S1, S5a and S3a has been seen to produce a bare
    pandas exception with its phase-1 raise removed (§4.4a).
19. **`§15.7a` passes**, and the `try/except model.FitError` companion has been seen to make `secondary`
    return seven estimates with every other test in §15 still green.
20. **All fourteen privates of §3.2 are implemented from the fences in §4.4b, §9.2, §10.5 and §8.2a** —
    none is written from prose. Ten of them had no fence in an earlier draft (§22.3 item 20).
21. **Each of §15.1a's ten branches raises a `SchemaError` whose message names ITS OWN condition**, not
    merely a `SchemaError`. Verified by asserting on the message: three of the ten frames were
    mis-specified in a way that only a message assertion catches (§15.1a, §22.3 item 24).
22. **`secondary` has been called end to end** on a cohort-level frame before any data-gated test runs,
    and the three audit entries' rendered tables pass §15.14's structural check. This is DoD-16's
    counterpart for the plumbing: DoD-16 says the pins reproduce, this says the function runs.

## 21. Verification record

**Normative.** A number in prose cites the row here that produced it. Where the two disagree **this
table is right and the prose is stale**, and an implementer who finds a conflict should fix the prose
and not re-derive the number.

| Claim | Where used | Verified |
|---|---|---|
| `Σw = 27.736623`, `Σh = 14.797901` over the ATO population | §6.4, §8.2, §13 | yes, run |
| cohort 93, `in_model` 92, one covariate-incomplete row, in Lugano | §4.2, §12.3 | yes, run |
| strata HUG 41 / Lugano 31 / CHUV 21; USZ absent; `center_USZ` dropped from all seven designs | §3.3, §7.2, §7.4 | yes, run |
| the seven [§11] denominators: 92 ×6, `tici_2b_3` 91 | §4.2, §14 | yes, run |
| the seven minority cells: 42, 24, 6, 5, 9, 25, 28 | §4.2, §9.3 | yes, run |
| paths: four full, one reduced, two unaugmented; `sich` and `ph2` are the two | §9.3 | yes, run |
| TICI's reduced design is 4 columns + intercept = 5 parameters | §7.2, §9.4 | yes, run |
| `center` has 4 declared levels, reference `HUG` | §7.2 | yes, read (`config.py`) |
| the seven Firth fits converge; iterations 6-17; `max|β|` 1.0268-6.6522; 0 rescales, 0 halvings | §7.4 | yes, run |
| `model.design` 17.40 ms and `model.firth` 6.65 ms for the five augmented; `propensity.fit` 13.39 ms; `primary` 4.73 ms; 71 ms wall-clock per replicate | §2, §12.1 | yes, run |
| the four full-list designs are identical, worth ~11 ms of the 17.40 | §2, §12.1 | yes, read + run |
| `Σw` splits 13.626360 treated / 14.110263 control; rows split 39 / 53 | §6.4 | yes, run |
| the correction is 2.9× (treated) and 3.8× (control) more aggressive than against rows | §6.4 | yes, run |
| audit 22 entries after `propensity.fit`, reconciling with Stage 8's 24 / 27 | §10.1 | yes, run |
| constant-model identity holds for `h`, `w` and unit tilts at 2.5e-16, 1.9e-16, 2.8e-17 | §8.4 | yes, run |
| dropped correction term: `−0.390000 = c₀ − c₁` with different constants; `+0.000000` with equal ones | §8.4, §15.8 | yes, run — re-measured on `golden_frame` |
| the six-candidate tilt-gap table; the non-monotone row at `2.459e-04` | §8.4 | yes, run — **round-1 exploration, retained for its ARGUMENT only; the fixture it chose is superseded by `tilt_frame`** |
| `tilt_frame` at seed 20260824: sd(d) 0.1698, Σw 25.913563, Σh 13.083855, `|tau_h − tau_w| = 1.030664e-02` | §8.4, §15.0.4, §15.8 | yes, run |
| constant-`d` sanity on `tilt_frame`: gap `0.000e+00` at d = 0.00, 0.05 and 0.40 | §8.4, §15.0.4 | yes, run |
| n=200,000 seed 20260825: estimand +0.320705, `tau_hat` +0.290456, ATO(mis) +0.295884, `tau_hat` − ATO(mis) −0.005428; correct-`e` err −0.004149 | §8.5, §15.10 | yes, run |
| bias tracks `sd(tau(X))`: −0.025338 / −0.009901 / +0.000463 at `|t|` 55.9 / 21.3 / 0.9 | §8.5, §15.10 | yes, run |
| n=92: mean −0.025665, sd 0.087669, 40.5% wrong sign | §8.5, §15.10 | yes, run |
| twelve seeds at n=200,000: heterogeneous [0.022950, 0.030249], constant [0.000176, 0.004929]; **empty band, both bounds inside**; 13 ms/draw | §8.5, §15.10 | yes, run |
| `ph2` crosses the threshold in 46.1% of replicates; `sich` 3.8%; others 0.0% | §9.5, §12.1 | yes, run — **with §12.3's workaround** |
| OR correction unreachable on v7; fires in 17.0% / 13.0% / 2.3% of replicates | §6.4, §13, §17 | yes, run — **with §12.3's workaround** |
| `model.firth` raises on 0.3% of `sich` replicates, 0% elsewhere | §9.5, §12.1 | yes, run — **with §12.3's workaround** |
| `sich` `max|β|`: 0.8558 to 80.9825, continuous, 42.7% > 8, 15.1% > 14; the seven-row table | §9.6 | yes, run — **with §12.3's workaround** |
| `separated_frame()`: `firth` converges in **6** iterations, rescales 0, halvings 0, `exp(β_treatment) = 441.005397`, `max_abs_beta` 6.089057, fitted p interior in [4.545e-02, 0.954546] | §3.3, §9.6, §15.7, §15.0.4a | yes, run |
| the unpenalised MLE on the same frame **raises `LinAlgError: Singular matrix`** (Newton, `maxiter` 200); **at the default 35 it does NOT raise** — `max_abs_beta` 67.8705, `converged` False; 67.4529 at 100; BFGS at 60 gives 33.0927 with `converged` **True** | §0, §3.3, §9.6, §15.7, §15.0.4a | yes, run — re-measured at implementation, §22.4 item 3 |
| `golden_frame()`: Σw 8.89, Σh 5.4257, p1_w 0.6191536748, p0_w 0.3340909091, RD_w 0.2850627657, OR_w 3.2404025938, `or_corrected` False | §15.0.2 | yes, run |
| `golden_frame()` 2-covariate fit: 4 cols, dropped `center_USZ`, 7 iterations, `max_abs_beta` 1.5877919570, rescales 0, tau **0.2233300473**, tau under w 0.2240676700 | §15.0.2 | yes, run |
| `golden_frame()` 3-covariate fit: 5 cols, 10 iterations, `max_abs_beta` **11.8282581834**, **rescales 1**, tau **0.0402879733** — a seven-fold collapse, not a sign flip | §15.0.2 | yes, run |
| `predict(fit, X) == fit.p` **exactly** (`np.array_equal`) on both `golden_frame` fits, with the prepended intercept and `FIRTH_ETA_CLIP` | §7.3, §15.6, §21b | yes, run |
| the scalar-intercept form is `allclose` but NOT `array_equal` to `fit.p`: max difference 1.11e-16 on 10 of 60 rows | §3.3, §7.3, §15.6 | yes, run |
| `1/(1+exp(-η))` unclipped OVERFLOWS at η = −800, returning 0.0 where `_probabilities` returns 7.124576406741285e-218 | §3.3, §7.3, §15.6 | yes, run |
| the correction can cross the null: Sw1 17.007 / Sw0 22.768, p0 0.00647 → **1.0201**; 12.404 / 20.945, p0 0.01283 → **1.0674**; near-even 13.626 / 14.110, p0 0.02 → 0.6484 (does not cross) | §6.3, §6.4, §15.5, §17 | yes, run |
| a constant outcome returns a **finite** OR and not `nan`: all-zero 1.3529411765, all-one 0.7391304348, both `corrected` True — so `0/0` is unreachable | §6.2, §15.1a, §16 | yes, run |
| without `predict`'s column check a reordered design returns finite probabilities `[0.6278, 0.3722, …]` rather than raising | §7.3, §15.6, §21b | yes, run |
| `augmented_rd`'s zero-total guard raises `SchemaError` on an all-zero weight vector | §8.1, §15.8 | yes, run |
| T1 gates T3 and T7 as well as T6 and T9: without it `marginal_odds_ratio` and `by_family` raise `AttributeError` on `OR_CONTINUITY` / `BINARY_OUTCOMES` | §20 | yes, run |
| `runner.py` does not exist anywhere outside `.venv`, so the earlier DoD-4 was unrunnable | DoD-4 | yes, run (`find`) |
| `data._md_table` does no escaping and sizes its separator from `len(rows[0])`, so `max\|beta\|` presents 11 header cells against a 9-cell separator | §10.3, §10.4, §15.14 | yes, read (`data.py:242-245`) + run |
| `model._probabilities` applies `FIRTH_ETA_CLIP = 500.0` to the same `1/(1+exp(-η))` form, for the reason at `model.py:317-318` | §3.3, §7.3 | yes, read (`model.py:312-320`, `config.py:310`) |
| Stage 6 solved the assertion-phasing problem by re-filtering inside one collected phase | §4.4a | yes, read (`model.py:152-154`) |
| `Propensity` carries `e`, `w`, `in_model` and **no `h`**, so `augmented_rd`'s earlier "it is `Propensity`'s" was false | §8.2a | yes, read (`propensity.py`) |
| augmented error bounded in [−0.344, +0.298] over 1778 synthetic fits | §9.6 | yes, run |
| augmentation RMSE ratios 1.19x / 1.11x / 1.06x on a constant-RD truth | §9.6 | yes, run |
| the synthetic construction reached `max|β|` 9.458 / 11.475 / 3.466 — it does **not** reproduce the tail | §9.6, §16 | yes, run |
| `propensity.fit` raises `SchemaError` on 104/400 = 26.0%; analytic 26.4% for a stratum of 31 | §12.3, §14 | yes, run |
| the golden vector, both covariate paths, to 10 decimals, **against the frame §15.0.2 now contains** | §15.0.2, §21b | yes, run — 30 / 30 |
| `1/(1+exp(-η))` is 1.0 at η = 800 and 1000 where `exp/(1+exp)` is `nan`; both 1.0 at η = 500 | §3.3, §15.6 | yes, run |
| `np.average` raises `ZeroDivisionError` on zero weights | §3.3, §5.1 | yes, run |
| `np.exp(-1000)` underflows to 0.0 silently | §3.3 | yes, run |
| `OUTCOME_MODEL_OVERRIDES` holds one entry; `outcome_model_covariates` excludes treatment and defaults | §7.1, §7.2 | yes, read (`config.py:485-505`) |
| `RARE_MINORITY_THRESHOLD` is 10 and its comment names the minority-cell reading | §9.1 | yes, read (`config.py:273-277`) |
| `Fit.p` is at the observed treatment | §3.3, §7.3 | yes, read (`model.py:104`) |
| `_record_exclusion`'s premise is Stage 2's A2 | §12.3 | yes, read (`propensity.py:400-432`) |
| [§10]'s "never substituted with a different estimator" | §9.5 | yes, read (`implementation_roadmap.md:496`) |
| the roadmap/SAP contradiction on TICI, and the roadmap's self-contradiction | §7.2, §9.2, §22 | yes, read — the roadmap **at commit `ae1416c`**, since §22 item 1 replaced the text; its Stage 0 DECISION 3; `statistical_analysis_plan.md:229-233` |
| Stage 0 quotes crude arm contrasts for two outcomes | §4.3 | yes, read (`out/stage0_data_inventory.md:202`) |
| a constant outcome's OR reduces to `(Sw0 + c)/(Sw1 + c)`, so EVEN arm totals give exactly 1.0 and cannot produce §6.2's pins; `Sw1 = 8`, `Sw0 = 11` give 1.3529411765 and its exact reciprocal 0.7391304348 | §6.2, §15.0.3, §15.5 | yes, run — §22.4 item 2 |
| the arm-normalised correction term deviates from `rd` by **+0.0812746102** on `golden_frame()`, against item 1's `1e-15` tolerance | §8.1, §15.8 item 4 | yes, run — §22.4 item 4 |
| `secondary_cohort()`'s design is rank **13 of 13**; three innocent-looking covariate patterns were exactly collinear (`sex == ivt`; `prestroke_mrs` in the span of the centre dummies; `ivt − atrial_fib − onset_unwitnessed + onset_wake_up == 0`) | §15.0.7 | yes, run — §22.4 item 5 |
| `secondary` on the workbook reproduces §7.4 and §9.3 exactly: iterations 6 / 9 / 7 / 7 / 7 and `max_abs_beta` 3.0933 / 1.0268 / 2.2840 / 3.8288 / 4.2117 across the five fitted models, 0 rescales and 0 halvings, `center_USZ` dropped from every design | §7.4, §9.3, §14 | yes, run — the first end-to-end run on v7, §22.4 |
| the audit ledger reaches **30** after the full pipeline and the log is byte-identical across `PYTHONHASHSEED` 0 and 1 | §10.1, DoD-4 | yes, run — §22.4 |
| `PSweight` agrees with `weighted_proportion` and `weighted_rd` on `golden_frame()` to **1e-12**, and the hand-assembled R augmentation reproduces `tau = 0.2233300473` | §19b | yes, run — the oracle, §22.4 |
| the `try/except model.FitError` downgrade leaves **1220 of 1224** non-data-gated tests green; only §15.7a's four catch it | §9.5, §15.7a, DoD-19 | yes, run — mutation, §22.4 |
| tilting by `w` **inside** `augmented_rd` fails §15.8 items 3 and 6; tilting by `w` **at the call site** fails item 6 and **passes item 3** | §8.2a, §15.8, DoD-5 | yes, run — mutation, §22.4 |

**Four rows carry a caveat and it is not decorative.** The replicate-based measurements — §9.5's
crossing rates, §6.4's correction rates, §9.6's `max_abs_beta` distributions and §12.1's `FitError`
rates — could not be taken through the landed pipeline, because `propensity.fit` raises on a quarter of
replicates (§12.3). They were taken with each drawn row renamed, which is one of the two candidate fixes.

**Round 3 changed what that caveat costs, without changing the caveat.** The renaming helper is now
`renamed_replicate` in `tests/fixtures_stage9.py` (§15.0.6) rather than a scratchpad, so the four rows
are **reproducible**: anyone can re-run them from the repository. What remains conditional is the
*inference*, not the measurement. **If Stage 10 adopts the other fix** — teaching the guard that a
replicate's rows are distinct draws — the numbers should be unchanged, because renaming affects only the
audit entry and not `e`, `w` or any mask. That is an argument, not a measurement, and it is the one place
this table's `yes, run` rests on one. An earlier draft made the same argument while the apparatus behind
it existed nowhere, which made the four rows unreproducible as well as conditional.

### 21b. This document’s own code, executed — 2026-08-25, round 3

Every `python` fence extracted **by content** and not by index [Stage 7 §18d], in **four** passes — the
fourth is `secondary` called end to end, which round 2 recorded as impossible and which became possible
once §4.4b and §10.5 supplied the ten privates it needs. The fence count changed from 12 to 20 because
§15.0 carries the fixtures and §4.4b and §10.5 carry the privates, and the map is derived by walking each
fence's AST for the names it defines — a hardcoded index set would have rotted on exactly this edit, which
is why Stage 7 §18d requires content.

```
  pass 1  ast.parse each fence alone                              20 / 20 parse
  pass 2  the §13 config additions executed FIRST (they are a
          paste-ready fence and every other fence reads them),
          then each definition fence into a namespace carrying
          the landed modules. The §3.2 signature-stub fence is
          classified out by content, not skipped by index.        18 / 18 execute, 0 failures

          injected  C.BINARY_OUTCOMES = ('mrs_0_2_90d', 'mrs_0_1_90d', 'tici_2b_3',
                                         'sich', 'ph2', 'death_90d', 'mrs_5_6_90d')
          injected  C.OR_CONTINUITY   = 0.5
          defined   BinaryEstimate, Secondary, _augmentable, _augmented_path,
                    _counterfactuals, _tilt, augmented_rd, marginal_odds_ratio,
                    outcome_model, predict, secondary, weighted_rd,
                    golden_frame, tilt_frame, separated_frame, dr_population, ato,
                    renamed_replicate, and DELTA / DR_SEED / DR_SEEDS_12 / TILT_SEED

  pass 3  the assembled definitions CALLED against §15.0.2's
          golden frame -- WHICH IS NOW IN THE DOCUMENT --
          and checked against every pin it states               30 / 30 checks pass
```

Pass 3 in full, since a pass count is not a result:

```
  Sw                                8.89                     pin 8.8900000000
  Sh                                5.425700000000001        pin 5.4257000000
  p1_w                              0.619153674832962        pin 0.6191536748
  p0_w                              0.3340909090909091       pin 0.3340909091
  weighted_rd                       0.2850627657420529       pin 0.2850627657
  marginal_odds_ratio               3.2404025937860506       pin 3.2404025938
    or_corrected                    False                    pin False
  augmented_rd, 2 covariates        0.2233300473286825       pin 0.2233300473
    max_abs_beta                    1.5877919570244747       pin 1.5877919570
    rescales / iterations / cols    0 / 7 / 4                pin 0 / 7 / 4
    predict(fit, X) == fit.p        True                     pin True
  augmented_rd, 3 covariates        0.04028797332778189      pin 0.0402879733
    max_abs_beta                    11.82825818339737        pin 11.8282581834
    rescales / iterations / cols    1 / 10 / 5               pin 1 / 10 / 5
    predict(fit, X) == fit.p        True                     pin True
  identity, c1=0.62 c0=0.23         5.551115123125783e-17
  identity, c1=0.45 c0=0.45         1.6653345369377348e-16
  augmented_rd zero-total guard     raised SchemaError
  predict reorder guard             raised SchemaError
    and WITHOUT the check            [0.6278, 0.3722, 0.6278, 0.3722, 0.8147,
                                      0.6278, 0.3722, 0.1853] -- finite, not an error
  _augmented_path x5                (42,full) (9,unaug) (tici 6,reduced)
                                    (sich 5,unaug) (ph2 10,full)
  by_family                         secondary 3, safety 4
```

The `__init__` of `BinaryEstimate` accepted every field §11 constructs it with, in the order §11 writes
them, which is the check that §3.1's field list and §11's constructor call have not drifted — a
disagreement no prose review catches and no reader notices.

**Two things pass 3 establishes that the round-2 pass could not.** `predict(fit, X) == fit.p` is asserted
as an **equality** and holds on both fits, which is the check that caught §7.3's two defects; and the
golden-frame pins are checked against **a frame this document contains**, so a reader can reproduce them
rather than take them on trust. The round-2 pass reported `20 / 20` against a frame that existed only in
a probe, which is why it could report success on numbers nobody else could obtain.

**And a fourth pass, which is the one round 2 said could not be done.**

```
  pass 4  `secondary` CALLED END TO END on a 40-record cohort-level
          frame with all seven outcome columns, the nine [§6]
          covariates and a hand-assembled Propensity              all four checks pass

    7 estimates, keys == C.BINARY_OUTCOMES in order                        yes
    by_family: secondary 3, safety 4                                        yes
    3 audit entries, in order: binary_estimable, binary_estimates,
      binary_outcome_models                                                 yes
    one outcome-missing record -> a 38-member population against
      the ATO's 39, so the §4.2 mask BITES on this frame too                yes

  pass 4b the three rendered tables, structurally (§15.14)

    binary_estimable        8 rows   header/sep/body cells 7/7/{7}   equal, no pipe
    binary_estimates        8 rows   header/sep/body cells 9/9/{9}   equal, no pipe
    binary_outcome_models   8 rows   header/sep/body cells 9/9/{9}   equal, no pipe

  pass 4c all ten precondition branches, each naming ITS OWN condition

    S1 S5a S3a S3b S2 S4 S5b S6 S7 S8       10 / 10 name the intended branch

  pass 4d the frozen-path parameter (§11.1)

    a complete map is accepted and the paths are carried unchanged           yes
    a partial map raises SchemaError                                         yes
```

**Pass 4 found two defects in §4.4b that no amount of reading found**, and both are recorded where they
belong rather than only here: S3b was specified in phase 2 and produces a bare `KeyError` from phase 2's
own masked read (§4.4), and three of §15.1a's ten frames were mis-specified in a way that made an earlier
phase fire instead of the branch under test (§15.1a). The second is why §15.1a asserts the *message*
rather than the exception type.

**What this does not establish.** Two things:

- **The workbook is not among the frames.** Pass 4's cohort is synthetic and 40 records; the seven
  denominators, the five augmentation paths and the audit ledger's 30 are asserted against v7 only in
  §15's data-gated tests. So `secondary` is now known to *run and to reject correctly*; it is not yet
  known to produce §14's numbers.
- **A fence that executes is not a fence that is right.** `augmented_rd` reproducing a pin computed by
  the same assembled formula shows the transcription is self-consistent, exactly as §19b's R oracle does
  one language over. §8.4's third test and §15.8 item 6 are what check the formula is [§8]'s.
- **The `config.py` additions were executed against `OUTCOMES` and injected onto the live module**, so
  pass 3's `_augmented_path` results depend on an injection rather than on a landed constant. T1 is what
  makes them real.

## 22. What this spec changed elsewhere

| # | Change | Files | The question it settles |
|---|---|---|---|
| 1 | Stage 9's **"minority cell below 10 → not augmented"** corrected to the amendment's rule: below the threshold **and no declared override** → not augmented | `implementation_roadmap.md` | Whether `tici_2b_3` is augmented. One stale sentence against three agreeing sources — the SAP amendment, the roadmap's own Stage 0 DECISION 3, and the adjacent bullet requiring TICI's override to be printed beside an estimate the guard had removed (§7.2) |
| 2 | Acceptance test 1's rationale corrected: it checks that the correction term is normalised and cancels, **not** that the tilting function is `h` | `implementation_roadmap.md` | Whether passing it establishes the estimand. It does not: the identity holds for every normalised tilt (§8.4) |
| 3 | Acceptance test 1 gains **"and the two constants must differ"** | `implementation_roadmap.md` | Whether the test can fail. With equal constants a dropped correction term passes it (§8.4) |
| 4 | A **third** acceptance test added: a non-constant `m_a(X)`, `h` asserted, `w` computed and shown to differ | `implementation_roadmap.md` | What pins the tilting function, since test 1 cannot (§8.4) |
| 5 | Acceptance test 2 gains the size it must be run at and the **constant-risk-difference** form of its companion | `implementation_roadmap.md` | Whether the test measures bias or noise. At n = 92 it measures noise (41.0% wrong sign); a no-interaction logistic companion is itself biased at `|t| = 10.8` (§8.5) |
| 6 | Stage 9 gains **"the augmentation paths are decided once, on the point estimate, and passed into the bootstrap"** | `implementation_roadmap.md` | Whether a replicate may choose its own estimator. `ph2` crosses the threshold in 46.1% of replicates (§9.5) |
| 7 | Stage 9 gains the odds-ratio rule: Haldane–Anscombe, only when degenerate, named beside the estimate | `implementation_roadmap.md` | What "kept finite" means (§6.3) |
| 8 | Stage 10 gains **`propensity.fit` raises `SchemaError` on 26.0% of stratified replicates** as a named blocker, and the `max|β|` distribution as a required diagnostic | `implementation_roadmap.md` | Whether Stage 10 can resample the landed pipeline. It cannot (§12.3) |
| 9 | **And one roadmap addition that is not an amendment**: Stage 9 gains its `**Spec:**` line, as Stages 1-8 have | `implementation_roadmap.md` | — |

**No SAP amendment.** Every substantive rule above is already the SAP's; what changed is the roadmap's
transcription of it, in two places, plus three test-design corrections that are this document's own and
belong to no protocol. §6.3's odds-ratio rule is the one genuinely new decision, and it is recorded as a
PI decision in §17 rather than as an amendment, because [§8] delegates it by silence rather than
specifying it wrongly.

### 22.1 Review round 1 — citations and quotations, 2026-08-24

Mechanical rather than by reading, because Stage 8's own rounds found stale line numbers by reading and
missed others. Every `file.py:N` and `file.md:N-M` in the document was extracted by regex and the cited
lines printed for comparison against what the document claims is there. **Seven were wrong.**

| # | Citation | Was | Is |
|---|---|---|---|
| 1 | `config.py:270-274` → `273-276` | started on `N_BOOT` and `SMD_THRESHOLD`, two lines above the comment quoted | the minority-cell comment block |
| 2 | `config.py:270-275` → `273-277` | same, in §21's row | the block plus the constant |
| 3 | `config.py:275` → `277` | the third line of a comment | `RARE_MINORITY_THRESHOLD: Final[int] = 10` |
| 4 | `config.py:493-505` → `489-503` | started mid-docstring | `def outcome_model_covariates` through its `return` |
| 5 | `model.py:113` → `104` | the design-matrix section banner | `p: np.ndarray  # (n,), fitted probabilities` |
| 6 | `test_config.py:409` → `450` | a blank line — **inherited from Stage 8 §14 and stale there too** | `def test_the_raw_name_scan_actually_fires()` |
| 7 | `implementation_roadmap.md:429-430` → `496` | shifted by this document's own amendment | [§10]'s "never substituted with a different estimator" |

**And one finding that is not a wrong number but a wrong kind of citation.** Four passages quoted
roadmap text *by line* that §22 replaces — the "kept finite" clause, the "not augmented" sentence and
both original acceptance tests. After the amendment those lines say something else, so the citation sent
a reader to text that contradicted the quotation. **A line reference to text the same commit deletes is
always wrong**, however carefully the number is maintained. All four now quote the passage as it read
before the amendment, name the §22 item that replaced it, and carry no line number; §21's row for the
contradiction cites the roadmap **at commit `ae1416c`** instead.

**What round 1 could not do.** It checked that each cited location contains what the document says it
contains. It did not check that the document's *reading* of that content is right — item 6 had been
stale in Stage 8 for a release without anyone noticing, which is what a mechanical check catches and a
careful reader does not.

### 22.2 Review round 2 — the claims themselves, 2026-08-24

Adversarial, against the SAP, the roadmap, the landed code and the document's own arithmetic. **Nine
findings; two changed a conclusion and one changed a fence.**

**The two that changed a conclusion, and both are in §2.**

1. **The cost figures timed the wrong thing, and the correction reverses which call is expensive.** The
   first draft timed `model.firth(X, y)` alone — 1.08 to 3.09 ms per outcome — and reported that
   `secondary` costs *"about 8 ms of nuisance fitting"*. But `outcome_model` calls `model.design` first,
   and re-measured as medians of 20-30 runs, **`design` is 17.40 ms against `firth`'s 6.65 ms for the
   five augmented outcomes — 72% of the cost.** The estimator is not the expensive part of the
   estimator.
2. **And the comparison it was drawn against was not like-for-like.** *"against `primary`'s measured
   1.157 ms [Stage 8 §2]"* compares five design-plus-fit cycles to **one `polr` fit**, not to one
   `primary` call. Measured, `primary` costs 4.73 ms and `propensity.fit` 13.39 ms, putting Stage 9's
   24.05 ms at 57% of a replicate's 42.2 ms of fitting — so the conclusion "Stage 9 is the cost centre"
   survives, but the factor it rested on was wrong by roughly an order of magnitude, in the flattering
   direction. The corrected §2 also yields something the wrong version could not see: **the four
   full-list outcomes share an identical design matrix**, so ~11 ms of the 17.40 is the same computation
   performed four times, and Stage 10 can have 22 s back over `N_BOOT`.

**The one that changed a fence.**

3. **§0's diagram contradicted §11 and §10.2 about when `_record_estimable` runs.** The diagram had it
   after the estimator loop; §11's code and §10.2's argument both require it **before**, so that a raise
   inside an estimate still leaves a log naming the seven populations. Three readers of the same document
   would have implemented two different orders. §0 now shows the two-loop shape §11 actually specifies.

**Four counting and completeness errors.**

4. **"Public surface, and it is six names"** listed five functions plus two dataclasses. Seven.
5. **The private list said "eleven"** and omitted `_augmented_path`, which §9.2 defines and §11 calls.
   Twelve.
6. **§0.2 said `model.polr` is "called, never edited".** It is not called — no binary estimate goes near
   the ordinal fitter. Now stated as untouched *and* uncalled, which is the Stage 8 §12 move for
   `RARE_MINORITY_THRESHOLD`: record the thing nobody should wire in.
7. **§12's data-flow tree omitted `outcome`, `family`, `reduced` and `dropped`** — two of which are read
   downstream, `family` by Stage 11's Benjamini–Hochberg grouping and `reduced` by Stage 14's obligation
   to print every reduced specification. A handover section that omits a field a later stage is required
   to read is the failure mode the section exists to prevent.

**Two claims that were true but stated more strongly than the measurement supported.**

8. **§6.4's "the row count per arm is around 46"** averaged two arms that are **39 and 53**, and the
   "roughly three times more aggressive" that followed is **2.9× in the treated arm and 3.8× in the
   control arm**. Measured, and now stated per arm. The related assumption — that `Σw` splits evenly —
   was checked rather than assumed and does hold: 13.626360 against 14.110263, which is not obvious,
   since `w` is `1 − e` in one arm and `e` in the other.
9. **§21's four caveated rows.** Round 2 confirmed the caveat is correctly scoped: the replicate
   measurements rest on a `case_id`-renaming workaround (§12.3), and the argument that the other
   candidate fix would give the same numbers — renaming touches only the audit entry, not `e`, `w` or any
   mask — is an argument and is labelled as one. It is the one place §21's `yes, run` rests on reasoning,
   and that is now said in §21 rather than only here.

**What round 2 did not find, and what it could not.**

- It found **no error in §8.4 or §8.5**, the two sections that correct the roadmap. Their algebra was
  re-derived independently and the six-candidate and four-`n` tables re-read against the probe output.
- It could **not** check `secondary` end to end, for §21b's reason: the function is specified and
  unexecuted, so its five ordering constraints, its three audit entries and S1-S8 are prose. §15's suite
  is the first place they run, and **the highest-value thing a third round could do is write §15.1 and
  §15.14 early** — before the estimators — because those are the two sections nothing has exercised.
- It could not evaluate the §6.3 decision, only its statement. Whether Haldane–Anscombe is the right
  rule is a PI question; §17 is where it sits, and §6.4's third measurement — that the choice moves two
  *intervals* and no point estimate — is the fact that decision should be made on. **Round 3 found a
  fourth measurement that also bears on it: the correction can invert the sign of a safety signal
  (§6.4).**

### 22.3 Review round 3 — executability, 2026-08-25

Rounds 1 and 2 checked citations and claims. **Round 3 assembled the fences and ran them**, which is the
only method that finds the class of defect below — and the method matters: measured on Stage 7, a
read-only review pass produced 0 findings where an executing pass produced 12 including 2 P1s, and 5 of 5
behaviour defects came from execution rather than reading. **Twenty-five findings; seventeen changed
behaviour or a number, and three changed what this document is for.**

Items 17-25 came from a **completeness audit run after the first sixteen were folded**, asking what the
revised document still claimed about itself that was not true. That audit is why item 20 exists, and item
20 is the one that most directly contradicted §0: ten of the fourteen privates were named, counted,
assigned to tasks and pointed at by the coverage map with **no body anywhere in the document**. A review
that had stopped at "all findings resolved" would have shipped a spec asserting it contained everything
an implementer needs while omitting ten functions.

**The three that changed what the document is for.**

| # | Finding |
|---|---|
| 1 | **Three of six fixtures had pins to ten decimals and no construction.** The golden frame was *"24 records, every value a literal"* and the 24 rows appeared nowhere; the tilt frame was named by a label from §8.4's table with no seed, no `n` and no generative model; the double-robustness population said *"seed pinned"* without writing the seed and needed a `DELTA` whose value never appeared. A fourth — the separated 40-record frame — was pinned in five places and in **DoD-9**. So T3, T5, T10 and T11 were unstartable and DoD 5, 6 and 9 unperformable, while §0 claimed the document was the sole source and §21 forbade re-deriving its numbers. §15.0 now carries all five constructions as code, **every pin was re-measured against them**, and §21b's pass 3 checks 30 of 30 against a frame a reader can run. The pins are therefore *not* the earlier draft's: `RD_w` is 0.2850627657 rather than 0.2265702811, and the second golden fit collapses tau to +0.0403 rather than inverting it to −0.2050. |
| 2a | **Ten of the fourteen privates had no body** (item 20 below, listed here because it belongs with item 1 in kind). Named, counted, assigned to tasks T7 and T8, pointed at by the coverage map — and specified nowhere. §4.4b and §10.5 are the bodies. |
| 2 | **The tilting function `h = e(1 − e)` had no home and no call-site test.** It was `e * (1.0 - e)`, inline, at §11's one tilt line, with `augmented_rd`'s docstring claiming it *"is `Propensity`'s"* — which is false, `Propensity` carries `e`, `w` and `in_model`. Every test that distinguishes `h` from `w` called `augmented_rd` directly with its own arrays, so an implementer typing `w` on that line shipped the wrong estimand with the whole suite green — the exact failure the Failure modes table calls *"the failure the spec exists to make catchable"*. §8.2a's `_tilt` is the home and §15.8 item 6 is the assertion that reaches the call site. |

**Ten that changed behaviour or a number.**

| # | Finding |
|---|---|
| 3 | **§8.5 concluded the opposite of its own table.** *"At n = 20,000 the two arms overlap — 0.0164 against 0.0134"*: those are `heterogeneous: min |err|` and `constant: max |err|`, and `0.0134 < 0.0164`, so the arms did **not** overlap and the table's own `ratio` column said `1.22x`. The defensible claim — 1.22× is too narrow to survive a numpy bump — is a different one. Round 2 re-derived §8.5 and reported no error in it. |
| 4 | **§3.1 said `reduced` moves as a block with `augmented`/`covariates`/`fit`.** It does not: `reduced` is `True` only for override outcomes, so it is `False` on the four `full`-path outcomes whose other fields are populated. §15.9 already said *"the four dependent fields"* and then listed three. The real block is `augmented`, `covariates`, `fit`, `dropped`. Also renamed §10.3's `reduced` **column**, which rendered a tuple under the name of a boolean field. |
| 5 | **Four counts disagreed with themselves.** The privates read *eight* in §1, *twelve* in §3.2 and *eleven* in §13 — round 2 fixed §3.2 alone and recorded it as fixed. §14's list was numbered 1, 2, 3, **5, 4**. DoD-13 required *"§16's seven items"* in `TODOS.md` against a deliverable of four. And §13's **paste-ready `config.py` comment** still carried the arithmetic round 2 had corrected everywhere else — *"13.9 per arm against a row count near 46"*, *"about three times more aggressive"* — so the wrong numbers were the ones shipping into the analysis code. The count now lives in §3.2 alone and DoD-17 requires every count re-derived. |
| 6 | **Seven of eight preconditions had no test.** T7 said *"the preconditions S1-S8; §15.1"*; §15.1 specified S1. The coverage map claimed *"71 / 71 routes with a test that can fail"* and annotated S1-S5 and S6-S8 with *"one frame per branch"* against frames the document did not specify. All eight are unreachable on the workbook by design, so nothing but the suite could ever reach them. §15.1a is one frame per branch, and the 71/71 total is gone — the map is now the count. |
| 7 | **§10.4's audit table header contained `max\|beta\|`.** `data._md_table` does no escaping and sizes its separator from `len(rows[0])`, so that header presents 11 markdown cells against a 9-cell separator. Stage 7 shipped this exact defect twice. It is invisible to DoD-4, because a table identically broken under both hash seeds is still byte-identical. Renamed to `max_abs_beta`, and §15.14 now asserts cell counts and pipe-freedom per entry on the rendered text. |
| 8 | **§6.2's stated reason for S8 was false.** *"`0 / 0` and `inf / inf` are both reachable ... and both give `nan` in numpy, silently. This is why S8 exists."* The degenerate branch fires **before** any division, so a constant outcome returns a finite `1.3529411765`. §16 item 3 — *"the sharpest open question the stage leaves"* — was argued on that premise, and the correction makes the case for `SchemaError` **stronger**: softening S8 yields a plausible number near 1 for an outcome with no events, not a legible `nan`. |
| 9 | **§7.3's `predict` did not reproduce `Fit.p`, and its reason for not clipping was wrong twice.** It added the intercept as a scalar where `firth` multiplies an intercept-prepended matrix — same mathematics, different summation order, `allclose` but not equal on 10 of 60 rows — so §15.6's oracle bullet was false. And it declined to clip on the grounds that *"the safe form has no overflow to clip — Stage 8's `POLR_ETA_CLIP` exists for the other form's problem"*: `model._probabilities` applies **`FIRTH_ETA_CLIP`** to exactly that form, for a load-bearing reason at `model.py:317-318`, and the form overflows at `η = −800`. §3.3 measured one tail and generalised. |
| 10 | **`_record_estimable`'s signature could not build §10.2's table.** It took `(df, populations, minorities, audit)`; the `lost` column and the detail string both need `ps.in_model`, which is not derivable from `populations` (those are intersections) and is never a column on the frame. `df[key].isna().sum()` is not `lost` — it counts over the cohort, not the ATO population. |
| 11 | **`secondary`'s signature could not express the contract §9.5 argues hardest for.** §12.1 and §14 both required the augmentation paths *"passed IN, not recomputed"*; the function took `(df, ps, audit)`. §14 then said a replicate loop *"can call `secondary` whole"*, which necessarily recomputed them — a direct contradiction, on the rule `ph2`'s 46.1% crossing rate exists to justify. §11's `paths` parameter and §11.1 are the resolution. |
| 12 | **`augmented_rd` divided raw**, while §3.3 called the check-before-the-mean *"load-bearing"* and §5.1 claimed Stage 9 *"inherits"* it. True of `weighted_rd`, false of `augmented_rd`, which §3.2 makes public and separately callable — so S7 is not between it and a caller. Three denominators, three guarded. |

**Four that changed only prose, and are recorded because §21's rule is that a number cites the row that produced it.**

| # | Finding |
|---|---|
| 13 | `augmented_rd`'s docstring cited *"§15.8's third test"* for the arm-normalised failure. That is the tilt test; the arm-normalisation is item **4** — and item 4 turns out to assert nothing item 1 does not, since the arm-normalised form fails item 1 decisively. Kept for localisation, and now described as such. (The `+0.482` this item quoted is a 92-row measurement and not `golden_frame()`'s `+0.0812746102` — §22.4 item 4.) |
| 14 | `proportion` is a mutable `dict` on a `frozen=True` dataclass **and aliased the object `weighted_proportion` returned**, so `estimate.proportion[1] = 99.0` changed a frozen estimate and the dict another caller held. §11 stores `dict(share)`. |
| 15 | §11 computes `weighted_proportion` twice per outcome — once inside `marginal_odds_ratio`, once inside `weighted_rd`. Deliberate: two `np.average` calls on ≤92 rows against `model.design`'s 17.40 ms, and routing the reported number through the public estimator is what makes §15.3 a test of it. Recorded as considered rather than left to be rediscovered. |
| 16 | DoD-4 ran `import runner`. **No `runner.py` exists** outside `.venv`; Stages 7 and 8 both spell the pipeline out in a heredoc. It also wrote `a.md`/`b.md` into `extended_bridging/`. |
| 17 | §15.8's heading read *"the three tests the roadmap's two become"* over a list of five items, while the coverage map said *"5 assertions, 4 companions"* — three different counts of two different things. Acceptance **tests** and **assertions** are now counted separately and the heading says which items realise which test. |
| 18 | §15.13's scan covered three shipped modules and the audit text, and omitted `tests/fixtures_stage9.py`, the three amended test modules and `tests/reference/aug_psweight.R` — all Stage 9 deliverables, and §19b's prose discusses the R script in AIPW terms. §8.3 also bans the phrase in *"this document"*, which the document then uses descriptively about a dozen times; the scan now has a **declared exemption** for lines that name or quote the prohibition, recorded as a judgement rather than a rule. |
| 19 | **S4 pinned `e` and `w`'s values and not their index ORDER.** `.loc[mask]` returns the series' order while `m1`/`m0` come from the frame's, so a permuted-but-equal index pairs each record's weight with another's fitted value — measured, `rd` moves from +0.224831 to +0.219968 and `tau` from +0.227928 to +0.235780, with nothing raising. §12's a-fortiori claim covered membership only. S3a is now `index.equals`, which is order-sensitive. |
| 20 | **Ten of the fourteen privates had no body.** `_assert_readable`, `_assert_secondary_inputs`, `_assert_estimable`, the three `_*_table`s, `_estimable_detail` and the three `_record_*`s were named, counted, assigned to tasks and pointed at by the coverage map, with no fence anywhere — against §0's *"everything the implementer needs is here"*. §4.4b and §10.5 are the bodies. This is the finding that most directly contradicted the document's own claim about itself, and it was raised by the outside voice rather than by the review sections. |
| 21 | `or_corrected` was documented as recording that a **cell was empty**; the guard tests the **float**, so a cell carrying weight 1e-18 sets the flag with the cell non-empty. One-directional and harmless — an empty cell always gives exactly 0.0 or 1.0 — but §10.3's column is now described as what it is. |
| 22 | §15.0.3's degenerate frames were listed as *"one frame each"* with no weights, and **two of the six only work at specific arm weight totals**: the null-crossing frame needs 17.007 / 22.768 to cross, where the workbook's near-even 13.626 / 14.110 gives 0.6484 and does not. A frame built with round weights would exercise the branch and miss the property. |
| 23 | **S3b was specified in phase 2 and belongs in phase 1** — found by pass 4, not by reading. S4's `ps.e[ps.in_model]` is a masked read, and pandas treats a non-boolean Series as *labels*, so an int64 `in_model` raises `KeyError: "None of [Index([-2, -2, ...])] are in the [index]"` before phase 2's collected `SchemaError` exists. This is §4.4a's own defect one check further on — the third time this document has had to re-apply the same rule — and it is the strongest available argument that the readability/judgeability boundary is structural rather than stylistic. |
| 25 | **DoD-4's replacement driver was itself wrong in three ways**, written while fixing the `import runner` defect: it called a nonexistent `data.audit()`, unpacked `data.load()` as `(df, src)` when it returns `(DataFrame, Audit)`, and routed the log through `Audit.write` — which derives its path from the source label, so both hash-seed runs would have written to one file. Stage 8 §15's driver was the correct form all along and DoD-4 is now that driver plus one line. Found by checking every name it calls against the modules, which is the check that should have accompanied the original fix. |
| 24 | **Three of §15.1a's ten frames were mis-specified**, each so that an earlier phase fired on incidental damage and the branch under test was never reached: `S6` left `e` populated while blanking `in_model` (S4 fires), `S7` zeroed `w` on off-mask rows too (S4 fires), `S3a` relabelled rather than permuted. Caught only because §15.1a asserts the **message names the branch**; `pytest.raises(SchemaError)` passes on all three. §15.1a now states the general rule: a frame testing phase *k* must be pristine for phases 1..*k*−1. |

**Two the round found and deliberately did not act on.**

- **`_augmentable`'s override is unconditional where the SAP's reduction is conditional** (*"**Where** the
  minority cell cannot support the full covariate set..."*). So a workbook revision taking `tici_2b_3`'s
  minority cell to 40 would still fit the two-name reduced model. Moot on v7, and §9.5's argument for
  unconditionality — naming an outcome fixes its estimator across replicates — is the stronger
  consideration. Filed as §16 item 8.
- **The §12.3 blocker itself.** Still not fixed here, for §12.3's stated reasons, which round 3 did not
  disturb. What changed is that the measurement apparatus is now in the repository (§15.0.6), so §21's
  four caveated rows are reproducible rather than resting on a scratchpad.

**What round 3 could not do.** It could not call `secondary` end to end: §15.0's fixtures are
outcome-level, and building a cohort-level one means either the workbook or a fixture this document still
does not specify. So the six ordering constraints, the three audit entries and S1-S8 remain **specified
and unexecuted** — §21b says so, and T7 and T9 are where they first run. Round 2 named the same gap; the
difference is that its recommended remedy (*"write §15.1 and §15.14 early"*) is impossible under §20's
task order, and §20 now says so and points at T0 instead.

**§22.4 is that gap closed.** T7 and T9 were written, `secondary` ran on the workbook, and the five
findings below are what running it produced.

### 22.4 Implementation — 2026-08-25

Rounds 1, 2 and 3 checked citations, claims and executability of the *fences*. **This round built the
stage**, which is the only method that finds the class of defect below — every one of the five is a place
where the document contradicted itself and no amount of reading the fences in isolation would have shown
it, because in each case **the fence was right and the prose about the fence was wrong**.

**What did NOT change, which is the more important half.** Every pin reproduced: §15.0.2's thirty golden
checks, `tilt_frame`'s `1.030664e-02` gap with the constant-`d` rows cancelling at exactly `0.000e+00`,
`separated_frame`'s six iterations at `exp(β) = 441.005397`, §8.5's twelve-seed bands `[0.022950,
0.030249]` and `[0.000176, 0.004929]` with the band between them empty. On the workbook, §4.2's seven
denominators, §4.2's seven minority cells, §7.4's five iteration counts and `max_abs_beta` values, §9.3's
five paths and §10.1's ledger of 30 all came back as written. **Round 3's claim that the document is
implementable from itself holds**; what follows are five defects in its margins, not in its arithmetic.

| # | Finding |
|---|---|
| 1 | **§15.11 forbade `dropna` in the added code, against §4.4b's own fence.** That fence's S5b is `arm.dropna()[~arm.dropna().isin(...)]`, and the `dropna` is not filling anything — it stops a missing arm being reported as an *out-of-range* arm, which is Stage 6's D3/D4 split one module over. DoD-20 makes the fences authoritative, so the bullet was the defect. §15.11 now forbids the FILLING verbs and asserts the one `dropna` is **located** in `_assert_secondary_inputs`, which is a stronger assertion than the ban it replaces: a second one anywhere fails. |
| 2 | **§15.0.3 said frame 7's arm totals are "even, ~4 / ~4" while §6.2 pins the numbers it produces to ten decimals, and the two cannot both hold.** On a constant outcome every event pseudo-count is `OR_CONTINUITY` alone, so the whole odds ratio is `(Sw0 + c)/(Sw1 + c)` — a pure function of the two arm totals. **Even totals give exactly `1.0`**, the one value that hides the asymmetry §6.4 exists to expose. `Sw1 = 8`, `Sw0 = 11` give `1.3529411765` and `0.7391304348`, and **those two are exact reciprocals**, which §15.5 now asserts as an identity rather than as two independent numbers. Round 3's §22.3 item 22 said two of the six frames have non-free weights; it is three. |
| 3 | **"The Newton route raises `LinAlgError`" is true only at `maxiter >= 200`.** At statsmodels' own default of 35 it does **not** raise: it returns `max_abs_beta 67.8705` with `converged = False`; at 100, `67.4529`; and BFGS at 60 returns `33.0927` with `converged = **True**` — a reported success at a value that is a function of the cap. DoD-9 gated completion on observing the raise, so the omission was load-bearing. **Stating the cap makes §9.6's point stronger, not weaker**: a raise reads as a numerical accident, where a coefficient that moves with the iteration cap is what "there is no maximum" actually looks like. Amended at six sites. |
| 4 | **§15.8 item 4's `+0.482` is not `golden_frame()`'s.** It was measured on a 92-record ATO-shaped construction, while item 1 — the test item 4 is said to be redundant against — is specified on the golden frame, where the arm-normalised form deviates by `+0.0812746102`. Both are decisive against item 1's `1e-15`, and §15.8 now asserts against **that tolerance** rather than against either measured gap, since the tolerance is what decides whether item 1 catches it. |
| 5 | **The cohort-level fixture round 3 recorded as unwritable is now §15.0.7, and building it found a constraint the document does not state anywhere.** `secondary` takes a cohort-level frame and every specified fixture is outcome-level, so T7 and T9 had no input. Writing one is not routine: a rank-deficient design does not give a bad estimate, it makes `model.firth` raise F3 and every augmented outcome fail at once — and **three innocent-looking covariate patterns were exactly collinear** (`sex = i % 2` *is* the treatment; `prestroke_mrs = i % 3` lies in the span of the centre dummies because `center` also cycles on three; `atrial_fib = (i // 2) % 2` with a period-four `onset_type` gives `ivt − atrial_fib − onset_unwitnessed + onset_wake_up == 0`). The **response** patterns carry the same trap one level on: `ivt` alternates on `i % 2`, so the obvious way to write a 20/20 outcome split is perfectly separated on the arm. §15.0.7 states both. |

**And two things this round verified that no earlier round could**, both by mutation rather than by
assertion, because the claims are about what a *wrong* implementation does:

- **§15.8's Failure-modes table is exactly right about the two tilt bugs.** Tilting by `w` inside
  `augmented_rd` fails items 3 and 6. Tilting by `w` **at the call site** fails item 6 and **passes item
  3** — which is what §8.2a predicted and is the reason item 6 exists.
- **§9.5's `try/except` downgrade is invisible to everything else.** With it in place `secondary` returns
  seven estimates and **1220 of 1224** non-data-gated tests still pass; only §15.7a's four catch it. That
  is the claim §15.7a makes about itself, observed rather than argued.

---

## GSTACK REVIEW REPORT

| Review | Trigger | Why | Runs | Status | Findings |
|--------|---------|-----|------|--------|----------|
| CEO Review | `/plan-ceo-review` | Scope & strategy | 0 | — | — |
| Codex Review | `/codex review` | Independent 2nd opinion | 0 | — | — |
| Eng Review | `/plan-eng-review` | Architecture & tests (required) | 1 | issues_found | 24 findings, all folded; 0 critical gaps |
| Design Review | `/plan-design-review` | UI/UX gaps | 0 | — | — |
| DX Review | `/plan-devex-review` | Developer experience gaps | 0 | — | — |
| Implementation | building the stage (§22.4) | Executability of the DOCUMENT, not of its fences | 1 | issues_found | 5 findings, all folded; every pin reproduced |

**Findings by origin (§22.3, §22.4):** 8 from the four review sections, 4 from the outside voice after
independent verification, 1 documentation-sync consequence, **8 more from a completeness audit run after
the first sixteen were folded**, and **5 from building the stage (§22.4) — every one of them a place
where the FENCE was right and the PROSE about the fence was wrong**, which is the class of defect a
review of fences cannot reach — including item 20, that ten of the fourteen privates had no body
anywhere in a document claiming to contain everything an implementer needs, and items 23-24, two defects
in the newly written material that only calling `secondary` end to end could find.

**OUTSIDE VOICE:** Codex was `ready` at preflight but failed at runtime with `401 Unauthorized` on token
refresh (`codex login` restores it). A Claude subagent ran the fallback with explicit extract-and-execute
instructions and returned 19 findings, 12 of them behaviour defects observed by running the fences. Three
were re-verified here before being acted on: `predict ≠ fit.p` (confirmed, and worse than reported —
`model.py:320` clips the identical form at `FIRTH_ETA_CLIP`), the Haldane–Anscombe null-crossing
(reproduced at 1.0201 and 1.0674), and the absent `runner.py` (confirmed).

**CROSS-MODEL TENSION:** one, on the twice-computed `weighted_proportion`. The performance review declined
to raise it; the outside voice raised it and attached a defect the review had missed — `proportion` was a
mutable `dict` on a `frozen=True` dataclass that aliased the returned object. The aliasing was fixed
(§11's `dict(share)`); the recomputation was kept and §22.3 item 15 records why.

**Prior learnings applied — three, and two changed the outcome:**
`pipe-in-audit-header-breaks-md-table` (10/10) caught §10.4's `max|beta|` before it shipped a third time.
`phase-assertion-helpers-by-readability-not-by-kind` (10/10) caught S1/S2 sharing a phase — **and then
caught S3b doing the same thing in the fix**, which is the third application of one rule in this
document's history. `outside-voice-must-be-told-to-execute` (9/10) is why the fallback was told to run the
fences. A fourth, `r-oracle-needs-R_LIBS_USER-search-not-fallback`, was checked and needs no action.

**Verification performed rather than asserted (§21b, four passes):** 20/20 fences parse, 18/18 definition
fences execute against the landed modules, **30/30 pins reproduce** from the fixtures now written into
§15.0, and **`secondary` runs end to end** — 7 estimates in registry order, 3 audit entries whose rendered
tables pass the structural pipe/cell check, **10/10 precondition branches naming their own condition**, and
the frozen-path parameter accepting a complete map and rejecting a partial one. Twenty-one cited line
references were checked against the code. Every number this revision introduces was measured.

**VERDICT:** ENG CLEARED — the spec is implementable from itself, and **as of 2026-08-25 it has been
implemented from itself**. All fourteen privates, all five fixtures and all three audit entries were
specified as code and were built from those fences without re-derivation. `uv run pytest` is green at
1330 tests with `data/` present and at 1224 passed / 106 skipped with it absent; the audit log is
byte-identical across two hash seeds; the R oracle agrees to 1e-12; and all 22 Definition-of-Done items
hold, including the five that require a wrong implementation to be *seen* to fail, which were run as
source mutations against a copied tree rather than asserted.

**What round 3 could not establish and the implementation did.** Pass 4's cohort was synthetic and 40
records, so `secondary` was known to run and to reject correctly but not to produce §14's numbers on v7.
It does: the seven denominators, the seven minority cells, the five iteration counts and `max_abs_beta`
values of §7.4, the five augmentation paths of §9.3 and the ledger of 30 all reproduce.

**What is still not established**, stated because §22.4 closed one gap and not this one: no reported
estimate on the workbook is pinned anywhere in git, by §4.3's design. The golden vector is synthetic.
`TODOS.md` carries the trigger (§16 item 5).

NO UNRESOLVED DECISIONS
