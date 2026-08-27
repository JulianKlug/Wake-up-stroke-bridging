# Stage 12 spec — all-centre standardisation

Implements roadmap Stage 12 [§14a]. Section references in brackets are to
`statistical_analysis_plan.md`. Numbers and decisions referenced as DECISION *n* are established in
Stage 0 and recorded in `../out/stage0_data_inventory.md`. Constants named in `SMALL_CAPS` are
Stage 1's and live in `config.py`; `stage1_config_and_data_contract.md` is their specification. The
classified frame this stage draws its population from is `stage4_eligibility_classification.md` §11;
the restriction it deliberately does **not** apply is `stage5_cohort_construction.md` §3; the ordinal
fitter it calls is `stage8_primary_outcome_estimator.md` §5; the standardised mean difference it
reuses is `stage7_balance_and_overlap.md` §5; and the resampler, the percentile definition and the
failure taxonomy it inherits are `stage10_bootstrap_engine.md` §4, §8 and §7 — whose §12.2 is written
as a handover *to this document*, names two things Stage 12 inherits and one it must decide for
itself, and is answered in §13.1.

**Precondition, and it is the one structural fact about where this document sits.** It is written
against the **landed Stages 1–10** and against the **specified-but-unlanded Stage 11**. Stage 12
depends on Stage 11 in exactly three places, each named where it occurs: `bootstrap.bucket` and
`bootstrap.collect` are public (Stage 11 §5.4, §13); `outcome.estimation_population` exists (Stage 11
§8.5) and is **not** used here, for the reason §4.1 gives; and the interval loop this stage makes a
third caller of is the one Stage 11 §16 left duplicated. Nothing else. Stage 12 fits no propensity
model, so `Propensity.spec`, the `Specification` registry and `balance._role`'s re-roling — Stage 11's
whole seam — are invisible to it.

**Status.** Written 2026-08-27 against the landed Stages 1–10. **Twenty-two probes were run before
any section was drafted, and six of them changed what this document says.**

- **§7.3 — two of the four quantities [§14a] lists as outputs are not new numbers.** The standardised
  mRS 0–2 risk difference *is* `RD_2` and the standardised mortality difference *is* `−RD_5`, to
  8.3e-17. This document's first draft specified four computations; there are two, and the other two
  are namings.
- **§7.4 — but their *intervals* are not reflections, and the reason is Stage 10 §8.2's pin.** Under
  `PERCENTILE_METHOD = "inverted_cdf"` the limits of `−X` and the reflected limits of `X` disagree by
  up to **1.7e-3** on the risk-difference scale — measured over 2000 draws — while under numpy's
  default `"linear"` they agree to 5.6e-16. The pin that makes [§10]'s p-value agree with its interval
  is the same pin that breaks reflection symmetry, and the mortality key therefore carries its own
  draws.
- **§9 — separation does not behave here as it does at Stage 8, and the guard moves.** The reported
  quantities are averaged probabilities, bounded in [0, 1] *by construction*: measured over 2000
  replicates every standardised probability is in [0, 1] and every risk difference in [−1, 1], while
  Stage 8's reported `exp(β)` on a separated fit is 6.5e15. `POLR_MAX_ABS_BETA` therefore guards a
  quantity this stage reports as a *model parameter* and nothing else, and it is specified at key
  granularity.
- **§12.2 — the hierarchical arm's quadrature must be ADAPTIVE, and this is measured rather than
  taken from the literature.** Non-adaptive Gauss–Hermite is off by **0.48 log-likelihood units at 31
  nodes** at σ = 3 and is **non-monotone in the node count**; adaptive quadrature is off by 8.9e-4 at
  **three** nodes. And 24.0% of replicates fit σ̂ above 1, where the non-adaptive rule is not
  converged — so this is not a refinement, it is what makes the arm computable at all.
- **§12.4, §12.8 — the fitter's shape is decided by a cost measurement, not by taste.** The same fit
  costs **39.0 s** under Newton with a central-difference Hessian and **0.76 s** under BFGS with an
  analytic gradient — 2000 replicates being 21 hours against **17.5 min**. The analytic gradient is
  therefore part of the specification and not an optimisation.
- **§12.5 — σ̂ reaches the boundary in 22.0% of replicates**, and the conditional-mode Newton
  overflows on `1/σ²` when it does. A prespecified floor is required, and it is a constant that
  changes an answer.

Two probes closed open `TODOS.md` items against their own stated triggers and **both answers are
negative**, which is recorded rather than quietly dropped: `POLR_MAX_ABS_BETA`'s calibration does not
bite here (§9.2), and the absolute convergence tolerance's change of scale does not either (§13.4).
One probe **falsified the stated reason** of a third item while confirming its trigger (§18, `TODOS`
row 3).

Every number below was produced by running code; none is carried — where a Stage 8, 10 or 11 number
is quoted for comparison it is labelled as theirs. It is the **sole source for the Stage 12
implementation**: everything the implementer needs is here, and anything not here is not to be
invented.

**Goal.** [§14a] whole. One pooled proportional-odds model over all eligible patients at all four
centres with `X` = the [§6] covariates minus `center`; standardisation by g-computation to the two
mRS distributions and the quantities [§14a] takes from them; the support check [§14a] requires and
the two sensitivity analyses it prescribes; and percentile intervals from a patient-level bootstrap
stratified by centre. Labelled an **ATE in the all-centre eligible population**, never the [§7] ATO.

**Not in scope.** [§14b]'s feasible-policy contrast (Stage 13), which reuses this stage's machinery
and is handed over in §19; the reporting layer and every label it owes (Stage 14 [§16]); any [§13]
sensitivity row (Stage 11); and any recomputation of a [§7]/[§8]/[§10] quantity — this stage refits
nothing Stages 6–10 fitted and reads none of their results.

**One thing this spec settles that no earlier stage could — and it settles five promissory notes at
once.** Every one of them was written *to this stage*, by name, and each is discharged in a numbered
section rather than acknowledged:

| Note | Where it was written | Discharged |
|---|---|---|
| `STANDARDISATION_COVARIATES` exists, is `PS_COVARIATES` minus `center`, and is labelled `[§14a]` | `config.py:562` | §5.1 |
| `cohort.treating_centres` is public *only* because Stage 12 needs its **complement** | `cohort.py:96-123`, Stage 5 §4.2 | §10.2 |
| `balance.smd` is public and is "Stage 12's as much as this stage's" | `balance.py:252`, Stage 7 §5.5 | §10.3 |
| Stage 12 takes the **unrestricted classified frame** and `build` is not on its path | Stage 5 §9, §10 | §4.1 |
| Whether the never-IVT stratum is resampled "is [§14a]'s question and not this document's" | Stage 10 §12.2, §17 | §13.1 |

**And one finding that arrives with the stage rather than being designed into it.** [§14a] lists four
things to report from the standardised distributions — the cumulative `RD_k`, the mRS 0–2 risk
difference and the mortality difference — and two of the last three are the first one under other
names (§7.3). That is not a defect in [§14a]: they are the quantities a clinical reader wants named.
It *is* a defect in any implementation that computes them twice, because two computations of one
number are two things that can disagree after an edit. They are therefore **derived from the
distributions and asserted equal**, in the manner [§16]'s constant-shift statement is computed rather
than written — and their **intervals are not** derived from each other, which §7.4 measures and is
the more surprising half.

---

## 0. Where Stage 12 sits

```
  data.load()                     ->  (df, audit)      25 cols, 126 rows   [Stage 2 §11]
  derive.derive(df, audit)        ->  df               31 cols             [Stage 3 §11]
  eligibility.classify(df, audit) ->  df               32 cols, 126 rows   [Stage 4 §11]
        │
        ├──────────────────────────────────────────────>  cohort.build(...)  93 rows
        │                                                 Stages 6-11, the [§7] ATO
        ▼   THE UNRESTRICTED FRAME.  restriction 1 is NOT applied  [Stage 5 §9]
  +--------------------------------------------------------------------------------------+
  |  standardise.py                                                                       |
  |                                                                                       |
  |  population(df, audit)               -> DataFrame  126 -> 107 -> 104         (§4.1)  |
  |     eligibility.retained             [§3 restriction 2] and NOT restriction 1        |
  |     model.complete_cases(., STANDARDISATION_COVARIATES)  &  outcome observed  (§4.1) |
  |                                                                                       |
  |  all_centre(df, audit)               -> Standardisation                      (§16)   |
  |     X, dropped = model.design(pop, (TREATMENT,) + STANDARDISATION_COVARIATES) (§5.3) |
  |     T1  the exposure column must survive                                     (§5.3)  |
  |     fit = model.polr(X, y)                UNWEIGHTED, one fit, no propensity  (§5.2) |
  |     regime  ->  X with the treatment column overwritten, NEVER a new design   (§7.1) |
  |     model.ordinal_probabilities(fit, X_a)  the ONE new name in model.py       (§6)   |
  |     expand onto MRS_LEVELS; unoccupied level -> structural zero               (§6.3) |
  |     average over the population; RD_k, mrs_0_2 = RD_2, mortality = -RD_5      (§7.3) |
  |                                                                                       |
  |  support(pop, audit)                 -> Support                              (§10)   |
  |     the treated-arm box over the non-categorical covariates                  (§10.1) |
  |     the complement of cohort.treating_centres, never a literal               (§10.2) |
  |     balance.smd with UNIT weights, treated vs never-IVT-centre               (§10.3) |
  |                                                                                       |
  |  hierarchical(pop, audit)            -> Hierarchical           [§14a] sens 2 (§12)   |
  |     model.polr_ri(X, y, groups)      A NEW ESTIMATOR, and the only one here   (§12)  |
  |     adaptive Gauss-Hermite, POLR_RI_NODES = 11, analytic gradient, BFGS      (§12.4) |
  |     the sigma -> 0 boundary is reached in 22.0% of replicates                (§12.5) |
  |                                                                                       |
  |  inference(pop, audit)               -> bootstrap.Bootstrap                  (§13)   |
  |     bootstrap.replicates(pop, body, N_BOOT, SEED, BOOT_STRATUM)                       |
  |     ALL FOUR STRATA, USZ included -- Stage 10 §12.2 left this here           (§13.1) |
  |     three failure groups, one per arm                                        (§13.2) |
  |     bootstrap.intervals(draws, lambda key: False)   NO p on ANY key          (§13.3) |
  +--------------------------------------------------------------------------------------+

  model.py gains TWO public names:   ordinal_probabilities  (§6),  polr_ri  (§12)
  bootstrap.py gains ONE:            intervals              (§13.3, the third caller)
  eleven `C.SchemaError` / `FitError` identifiers, and they are the T series   (§14)

                  ->  Standardisation, Support, Hierarchical, Bootstrap          (§3.1)

  the classified frame comes back unchanged                                      (§0.2)

  audit: the Stage 2-4 path's 12 entries, then TEN of this stage's              (§15)

  every standardised probability, every risk difference, every interval limit and
  every centre intercept on the workbook is DELIBERATELY ABSENT from this
  document. They are in the gitignored log.                                      (§4.3)
```

### 0.1 Why this is one module and not two, and the argument against it is real

`standardise.py`, new. It implements [§14a] whole — the model, the standardisation, the support check
and both sensitivity analyses — and it is written so that [§14b] adds a *regime* rather than a second
g-computation (§19).

**The rejected alternative is a split between the machinery and the analysis**, and it has a genuine
argument: `ordinal_probabilities` and the averaging step name no outcome, no covariate and no
population, and Stage 13 needs both unchanged. A `gcomp.py` holding the transport machinery and a
`standardise.py` holding [§14a]'s choices would put the reusable half behind a boundary that Stage 13
would then be unable to break by accident.

It is declined, and the reason is that **the boundary already exists in the right place and is not a
module.** The genuinely general half is *category prediction from a fitted proportional-odds model*,
and that belongs in `model.py` beside `predict` (§6.1) — where Stage 6 §11's rule already puts every
estimator-shaped function, and where Stage 9's `m_a(X)` found its counterfactual helper. What is left
after that extraction is the duplication-and-average loop, which is nine lines and names the
treatment column; a module boundary around nine lines that Stage 13 will call with one argument
different is a boundary drawn where the variation is, not around it. **The cost is named rather than
denied**: `standardise.py` will be imported by Stage 13 for its regime seam, so a Stage 13 change can
reach [§14a]'s code. §19 fixes what may and may not move.

The second half of the argument decided it. [§14] is one section of the plan, its two subsections
share one estimator and one set of assumptions, and the assumption paragraph — *"the conditional
effect learned in IVT-using centres transports to never-IVT centres … not testable from these
data"* — governs both. A split between "the machinery" and "[§14a]" would put the assumption in
neither file.

### 0.2 What Stage 12 does not touch

**It fits no propensity model, and this is the first stage of which that is true since Stage 5.**
[§14] says so in its own first paragraph — *"No propensity model is fitted anywhere in §14"* — and it
is not a simplification but the reason the section exists: a propensity model cannot produce a
contrast at a centre with no treated patients, so the estimand had to change. There is no `e`, no
`w`, no `h`, no ESS and no overlap weight anywhere in this stage. `propensity.py` is not imported.

**It reads no Stage 6–11 result.** It takes the classified frame and nothing else. In particular it
does not read `Primary`, `Secondary`, `Bootstrap`, `Balance` or any `Propensity`: §15.12's equivalent
here is that the Stage 1–11 pipeline run with and without this stage produces byte-identical numbers
and an identical audit prefix (§20.11).

**It does not call `cohort.build`, and Stage 5 §9 already fixed that.** Restriction 1 removes the
never-IVT centre, which is precisely the population [§14a] exists to include. `build` offers no seam
to skip it and must not acquire one — Stage 5 §15 declines a `centres=` parameter by name, and §4.1
here declines it again from the other side.

**It adds no [§7] estimand and reports no ATO.** Every output of this stage carries the
different-populations statement, which Stage 14 [§16] is required to print and §17.2 hands over as a
computed clause rather than a sentence.

---

## 1. Deliverables

| Path | Contents |
|---|---|
| `extended_bridging/standardise.py` | **new.** [§14a] whole: five public names (§3.2), four returned records (§3.1), seventy bootstrap estimand keys (§3.3) |
| `extended_bridging/model.py` | **amended, and it gains a new estimator.** `ordinal_probabilities` — category prediction from a `PolrFit`, the counterpart of `predict` (§6); and `polr_ri` — the random-centre-intercept proportional-odds fit, adaptive Gauss–Hermite, analytic gradient, BFGS (§12). `RIFit`, its frozen record. `polr`, `firth`, `design`, `predict` and `complete_cases` are untouched, and §20.11 asserts that by number |
| `extended_bridging/bootstrap.py` | **amended.** `intervals(draws, tested)` becomes the ninth public name — the loop Stage 11 §16 left duplicated across two callers, this stage being the third (§13.3). Both Stage 11 callers move onto it. `resample`, `replicates`, `percentile_ci` and `bootstrap_p` are unchanged |
| `extended_bridging/config.py` | **amended.** `POLR_RI_NODES`, `POLR_RI_SIGMA_FLOOR`, `POLR_RI_MAX_ITER`, `SUPPORT_COVARIATES` and three `FAILURE_BUCKETS` tokens, in the blocks §18 names (§18, fence) |
| `extended_bridging/tests/test_standardise.py` | **new.** §20 |
| `extended_bridging/tests/fixtures_stage12.py` | **new.** §20.0, and that row is the one place its contents are enumerated |
| `extended_bridging/tests/reference/polr_ri_clmm.R` | **new.** The `ordinal::clmm` oracle for `polr_ri`, gated as `polr_clm.R` is (§12.7) |
| `extended_bridging/tests/test_model.py`, `test_bootstrap.py`, `test_config.py`, `test_reference_r.py` | **amended.** The two new estimators, the ninth public name and the surface count that moves with it, the four new constants, the second R gate (§20) |
| `extended_bridging/implementation_roadmap.md` | **amended.** Stage 12 gains its `**Spec:**` line; five corrections and three additions (§27) |
| `TODOS.md` | **amended.** Two items close, one has its stated reason corrected while its trigger fires, **six** are new (§21) |

`out/logs/audit_<label>.md` is an **output, not a deliverable**. For this stage it is the only place
any standardised probability, any risk difference, any interval limit and any centre intercept
exists. It is gitignored.

**Nothing under `specs/` may quote a case identifier, and nothing here does.** Stage 8 §4.3 set that
rule, Stage 10 §4.5 tightened it to interval limits and p-values and Stage 11 §1 to E-values and
adjusted p-values. This document tightens it once more, to **standardised probabilities and centre
random intercepts**, and §4.3 is the reason.

---

## 2. Environment

Unchanged from Stage 11: `uv`, Python 3.12.12, pandas 2.3.3, numpy 1.26.4, statsmodels 0.14.6.
**No dependency is added, and one was specifically checked for and not needed.** `polr_ri` requires
Gauss–Hermite nodes and a quasi-Newton loop; `numpy.polynomial.hermite.hermgauss` supplies the first
and §12.4 writes the second, so `scipy` stays test-only by policy. §12.4 records what a `scipy`
dependency would have bought and why it is not worth the reproducibility surface.

```bash
cd extended_bridging && uv sync && uv run pytest -v
```

**Cost, measured (§26).** This stage adds one `bootstrap.replicates` call to the pipeline, whose body
runs three arms over each drawn frame:

```
   the pooled arm and the support-restricted arm, together      44 - 50 s / 2000   (§13.4)
   the hierarchical arm                                        1050   s / 2000   (§12.8)
                                                              ─────────────────
   Stage 12's bootstrap, one call, three arms                  ~18.4 min
```

**And the hierarchical arm is 95% of it**, which is the number that makes §12.4's optimiser choice a
specification and not an optimisation: the same arm under the first prototype this document wrote is
**21 hours** (§12.8).

---

## 3. Module shape

### 3.1 What the stage returns

Four frozen dataclasses, none carrying a verdict, following Stage 8 §3's rule that a record reports
what it computed and a caller decides what it means.

```python
@dataclass(frozen=True)
class Standardisation:
    """One [§14a] standardisation: one fit, one averaging population, two regimes."""

    population: str                      # "all eligible" | "treated support" — §11
    n_average: int                       # the averaging population's size — §7.2, [§11]
    n_fit: int                           # the FIT's population; equal unless §11 changes it
    distribution: dict[int, dict[int, float]]   # arm code -> mRS level -> P, over MRS_LEVELS
    cumulative: dict[int, dict[int, float]]     # arm code -> threshold -> P(Y <= k)
    rd: dict[int, float]                 # RD_k, keyed by MRS_THRESHOLDS — §7.3
    mrs_0_2: float                       # == rd[2]; [§14a]'s name for it — §7.3
    mortality: float                     # == -rd[5]; [§14a]'s name for it — §7.3
    conditional_log_odds: float          # beta. NOT called an odds ratio, and §8 is why
    conditional_odds_ratio: float        # exp(beta), and it is CONDITIONAL — §8
    measure: str                         # the guard string, always _CONDITIONAL — §8
    fit: model.PolrFit
    dropped: tuple[str, ...]             # design columns dropped as constant — §5.3
```

**There is no field called `odds_ratio`**, and that is a decision rather than a naming preference:
`Primary.odds_ratio` (`outcome.py:155`) is the [§8] marginal common odds ratio, and a Stage 14
formatter reaching for the same attribute on two records would print a conditional quantity under a
marginal label. §8 is the argument; the absent name is the enforcement.

```python
@dataclass(frozen=True)
class Support:
    """[§14a]'s support check: the treated-arm box, who is outside it, and the baseline table."""

    box: dict[str, tuple[float, float]]  # covariate -> (treated min, treated max) — §10.1
    inside: pd.Series                    # boolean, TOTAL on the population's index — §10.1
    outside_by_centre: dict[str, int]    # centre -> records outside the box
    never_ivt: tuple[str, ...]           # the complement of treating_centres — §10.2
    n_never_ivt: int
    n_never_ivt_outside: int             # [§14a]'s "the proportion … falling outside it"
    baseline: tuple[balance.CovariateBalance, ...]   # treated vs never-IVT centre — §10.3
```

```python
@dataclass(frozen=True)
class Hierarchical:
    """[§14a]'s sensitivity 2: one common treatment effect, a random centre intercept."""

    standardisation: Standardisation     # the same record, from the b_hat-conditioned prediction
    sigma: float                         # sigma_hat, the between-centre SD — §12.1
    at_floor: bool                       # sigma_hat reached POLR_RI_SIGMA_FLOOR — §12.5
    intercepts: dict[str, float]         # centre -> b_hat, the conditional mode — §12.6
    posterior_sd: dict[str, float]       # centre -> the curvature at the mode — §12.6
    nodes: int                           # POLR_RI_NODES, recorded because it changes an answer
    fit: model.RIFit
```

```python
@dataclass(frozen=True)
class RIFit:            # model.py, beside PolrFit
    """A random-intercept proportional-odds fit. No standard error, for PolrFit's reason."""

    beta: np.ndarray            # (m,), one per design column — NO intercept
    alpha: np.ndarray           # (K,), ascending
    sigma: float                # the between-group SD, > 0
    groups: tuple[str, ...]     # the group labels, in the order b_hat is keyed
    b: np.ndarray               # (G,), the conditional modes — §12.6
    b_sd: np.ndarray            # (G,), the curvature at each mode
    categories: tuple[int, ...]
    columns: tuple[str, ...]
    nodes: int                  # the quadrature node count this fit used — §12.3
    iterations: int
    converged_on: str           # "likelihood" | "score"
    first_step_norm: float
    rescales: int
    halvings: int
    at_floor: bool              # sigma reached POLR_RI_SIGMA_FLOOR — §12.5
```

`RIFit` carries `nodes` and `PolrFit` does not, because the node count **changes the answer** (§12.3)
while `POLR_TOL` changes only how precisely the same answer is found. A record that omitted it would
let two fits with different numerical content compare equal on every field.

### 3.2 The public names

`standardise.py` exposes five: `population`, `all_centre`, `support`, `hierarchical`, `inference`.
Everything else is private. `model.py` goes from six public names to eight
(`ordinal_probabilities`, `polr_ri`); `bootstrap.py` from eight to nine (`intervals`).

### 3.3 The estimand keys — seventy, and the count is derived

Twenty-three per arm, three arms, plus one:

```
   beta                                            the CONDITIONAL log-odds — §8      1
   rd_0 … rd_5                                     len(C.MRS_THRESHOLDS)              6
   mrs_0_2, mortality                              [§14a]'s two namings — §7.3        2
   dist1_0 … dist1_6, dist0_0 … dist0_6            2 * len(C.MRS_LEVELS)             14
                                                                                    ───
                                                              per arm                23

   arms:  ""  (pooled)   "support."   "hier."                                 23 * 3 = 69
   plus:  hier.sigma                                the between-centre SD — §12.5     1
                                                                                    ───
                                                                                     70
```

**The counts are computed from `C.MRS_LEVELS` and `C.MRS_THRESHOLDS` and never written as 6, 14 or
23.** Stage 8 §12 made the level set derived from the plausible range so the two cannot disagree;
this stage's key set is the next link in that chain, and `test_standardise.py` §20.3 asserts the total
against the config rather than against 70.

**No key carries a p-value, and the rule is passed to `bootstrap.intervals` as `lambda key: False`.**
[§14]'s Inference paragraph prescribes *"percentile intervals"* and nothing else; [§10]'s p is defined
for the treatment coefficient of the **[§8] weighted** proportional-odds model, which is a different
estimator on a different population, and a p-value here would be a test [§14] does not ask for
against a null [§14] does not state. Stage 11 §5.5 passes the same rule for the same shape of reason,
and §13.3 is why the two now pass it to one function.

---

## 4. The population [§14a, §11] — and it is not `cohort.build`'s

### 4.1 The ledger, measured

```
   the classified frame                                                126        [Stage 4 §11]
   ├─ [§3] restriction 2, eligibility        eligibility.retained      107   -19
   │        restriction 1 is NOT applied.  USZ stays.                            [§14a]
   ├─ [§11] covariate-complete on STANDARDISATION_COVARIATES           106    -1
   └─ [§11] the primary outcome observed                               104    -2
```

Three records leave for [§11] reasons and **the spec names what each is, because the two causes fall
on different centres and that is a fact about this workbook worth a reader's attention**:

- The **one** covariate-incomplete record is incomplete on `core_ml` **and** `tmax6_ml`, both of which
  are in `STANDARDISATION_COVARIATES` and in `PS_COVARIATES`. It is the same record Stage 11 §5.2
  identified as the [§7] cohort's single covariate-incomplete row — measured, and it is why Stage 11's
  full-covariate arm loses nobody the primary did not.
- The **two** outcome-missing records are **both at USZ**. The centre that contributes no contrast is
  also the one whose outcome ascertainment is incomplete, so [§14a]'s target population loses 2 of the
  14 patients its whole reason for existing is to include. That is reported per centre in the
  population audit entry (§15) and it is a denominator [§11] requires stated.

The resulting per-centre × arm table:

| centre | EVT alone | bridging | n |
|---|---|---|---|
| HUG | 11 | 30 | 41 |
| CHUV | 14 | 7 | 21 |
| Lugano | 28 | 2 | 30 |
| USZ | 12 | 0 | 12 |
| **all** | **65** | **39** | **104** |

**`treating_centres(frame)` returns `('HUG', 'CHUV', 'Lugano')` and its complement is `('USZ',)`,
which equals `EXPECTED_NEVER_IVT` — measured, and asserted rather than assumed** (§10.2).

**The mask is spelled here and not taken from `outcome.estimation_population`.** Stage 11 §8.5
extracts the [§11] mask `ps.in_model & df[outcome].notna()` as a public one-liner, and it is the wrong
one: its first conjunct is a `Propensity`'s covariate-completeness over `PS_COVARIATES`, and this
stage has no `Propensity` and a different covariate list. Calling it would require constructing a
`Propensity` this stage is forbidden to fit. The two masks are the same *shape* and different
*content*, and a shared helper parameterised over both would be one function whose two callers agree
on nothing but the `&`.

**`derive.derive_cohort` is not called, and must not be.** It computes `core_above_median` on the
frame it is given, and Stage 5 §4.5 established that the median is the *cohort's* — 5.0 mL over 93
records against 6.0 mL over 126. `core_above_median` is a [§13] subgroup variable, Stage 11's and not
this stage's; [§14a] names no subgroup. Stage 3's guard raises if a caller tries, and this stage never
does.

### 4.2 The audit step names, and Stage 11 §4.2's collision applies here in a second form

`Audit.entry` is first-match (`data.py:273-275`). Stage 11 §4.2 measured that two `propensity.fit`
calls over one cohort record byte-identical step names, so every programmatic read returns the first
while the rendered log looks complete — which is why its seam became a named record.

**Stage 12 meets the same hazard from the other direction**, and it is worse because it is silent in
both directions: in a Stage 14 driver, `cohort.build` and `standardise.population` both remove rows
from the same classified frame under `kind="cohort"`, and both are row-removal ledgers with case
identifiers. A step name reused between them would make `audit.entry("cohort", …)` return whichever
ran first, and the two populations differ by 11 records and one whole centre.

Every Stage 12 step name is therefore **distinct by construction and asserted distinct** (§15,
§20.9). The step names `cohort.build` owns and this stage may not use are `restrict_centres`,
`restrict_eligibility`, `cohort_flow`, `absence_by_cohort_column`, `core_above_median` and
`constant_covariates` — measured against the landed module, and `test_standardise.py` §20.9 asserts
the intersection of the two step-name sets is empty rather than asserting a list.

### 4.3 What may not be quoted, and this stage tightens the rule

Stage 8 §4.3 forbade case identifiers in `specs/`; Stage 10 §4.5 added interval limits and p-values;
Stage 11 §1 added E-values and adjusted p-values. **Stage 12 adds standardised probabilities and
centre random intercepts.**

A standardised probability is worse than an interval limit for the reason an E-value is: it reads as
a directly clinical number — *"this many patients in 100 would be independent under bridging"* — and
a reader who meets it in a version-controlled design document will quote it. A centre random
intercept is worse still, because it is a **centre-level outcome contrast** on four named hospitals,
and this analysis is not powered to make one and does not claim to. Both belong in the gitignored log
and in the manuscript, and neither belongs here.

What this document *does* quote is structural: counts, rates, iteration counts, convergence routes,
timings, numerical error magnitudes, and — in §12.2 and §12.3 alone — the fitted between-centre SD,
which appears there as **the quantity whose numerical stability is being demonstrated** and nowhere as
a result. Its interval, and any reading of between-centre heterogeneity, are the log's.

---

## 5. The pooled proportional-odds model [§14a]

### 5.1 The covariate set, and centre is absent for a reason that is not parsimony

    logit P(Y <= k | A, X)  =  alpha_k  +  beta*A  +  gamma'X ,    k = 0 … 5

with `X = C.STANDARDISATION_COVARIATES`, which `config.py:562` already declares as `PS_COVARIATES`
minus `center` and labels `[§14a]`. The list is not written out here and must not be written out
anywhere: Stage 1's rule is that a covariate list lives in `config.py` and that none may appear as a
second literal, and roadmap invariant 3 exists because two hand-maintained lists drift.

**Centre is omitted because the identifying assumption is exactly that potential outcomes do not
depend on centre given `X`.** [§14a] states it and the consequence is worth spelling: including
`center` would make the model contradict the assumption the analysis runs on, *and* it would make the
standardisation uncomputable in the direction that matters — a `center_USZ` coefficient is estimated
from a stratum with no treated patient, so the IVT counterfactual for a USZ patient would be an
extrapolation along a coefficient no USZ data can identify. Omitting centre is what makes the
transport explicit rather than hidden inside a dummy.

**[§15] permits this only inside [§14].** *"Omitting centre and modelling across centres are permitted
inside §14 only, where extrapolating to centres that never used IVT is the declared purpose and is
labelled as such."* `test_config.py` asserts `STANDARDISATION_COVARIATES` is `PS_COVARIATES` minus
exactly `{"center"}` (§20.2), so the permission cannot be widened by editing a tuple.

**Contraindication status is not a covariate**, [§14a] having already restricted to eligible patients;
the `eligibility` column is the restriction and never a regressor. `penumbra_ml`, `hir`,
`first_image_mri` and the four vascular risk factors keep their [§6] exclusions unchanged — [§14a]
says *"exclusions otherwise as §6"*, and `BALANCE_SET` is what the support check's baseline table
ranges over (§10.3) precisely so that the excluded ones are still *seen*.

### 5.2 One fit, unweighted, and it is `model.polr` unchanged

    fit = model.polr(X, y)          # w is not passed, and None is not a synonym for ones

`polr`'s docstring is explicit that `w=None` is the unweighted fit and *"is not a synonym for ones: a
caller who passes nothing has said something different from a caller who passes ones"*. [§14a]
weights nobody — there is no propensity model to weight by — so this stage passes nothing, and
§20.4 asserts the call site rather than the arithmetic.

**The parametrisation is [§14a]'s own**, which is the happy accident `polr`'s docstring already
records: `alpha_k + x'beta`, not the `theta_k − x'beta` of `statsmodels`, R's `polr` and `ordinal::clm`.
Stage 8 adopted it because [§14a] writes it that way; this is the stage [§14a] was talking about, so
no sign is converted anywhere in this module. §12.7's oracle is where the negation is handled, once,
and `polr_clm.R`'s banner already explains why the asymmetry is asserted by name.

**Measured on the estimation population**: 10 design columns, nothing dropped, rank 10; 6 cutpoints;
16 parameters on 104 records; converged in **5 iterations on the likelihood criterion**, 0 rescales, 0
halvings. For comparison Stage 8's fit is **1** column and 6 cutpoints on 92 weighted records — which
is the difference §9.2 and §13.4 both turn on.

### 5.3 The design, and T1: the exposure must survive

    X, dropped = model.design(pop, (C.TREATMENT,) + C.STANDARDISATION_COVARIATES)

`design` drops constant columns after dummying and returns their names; it does not raise on them.
Stage 8 §5 wraps the same call in G6, which asserts the exposure column survived, because a design
that lost the treatment column fits happily and returns a result in which `beta` does not exist.
**T1 is that assertion for this design**, and it differs from G6 in one way that matters: G6 also
asserts nothing else was dropped, because Stage 8's design has one column. Here nine covariates may
legitimately lose a level in a replicate — `onset_type_wake_up` going constant is a real event on a
resample — and dropping one is a change in the *model*, not a failure. So:

- **T1 raises `model.FitError` when `C.TREATMENT` is not in `X.columns`.** `FitError` and not
  `SchemaError`, because [§10] must drop and count such a replicate rather than crash: this is Stage 9
  §12.3's S8 reclassification applied at the point it is created rather than after 1.0% of replicates
  hit it.
- **A dropped covariate is recorded and not raised on.** `Standardisation.dropped` carries the names
  and the fit audit entry prints them, following `outcome.py:1342`'s rule that *"`center_USZ`
  vanishing from every design is a fact about the cohort nobody should"* have to rediscover.

**Measured over 2000 replicates: `dropped` is empty in 100.00% of them.** T1 is therefore specified
and never fires on this workbook, which is recorded rather than treated as a reason to skip it —
Stage 5's `_assert_cohort_inputs` and Stage 4's E1–E5 are written to the same rule, for the workbook
that has not arrived yet.

---

## 6. Category prediction from a `PolrFit` — the one genuinely new capability

### 6.1 Nothing in this pipeline can do it, and the gap is exact

`model.predict(fit, X)` (`model.py:541`) evaluates a **Firth binary** `Fit` and returns one
probability per row. `model._ord_pieces` (`model.py:653`) computes the two bracketing cumulative
probabilities for **one** category per row — the row's own observed category — and is private and used
only inside the log-likelihood, the score and the Hessian. `tests/test_model.py:1469` reaches into it
to build a per-row probability of the observed category, for separation detection.

None of those is `P(Y = j | x)` for **every** `j`, which is what g-computation averages. It is
therefore new, and it goes in `model.py` beside `predict` for the reason §0.1 gives.

```python
def ordinal_probabilities(fit: PolrFit, X: pd.DataFrame) -> np.ndarray:
    """P(Y = c | x) for every c in `fit.categories`: (len(X), len(fit.categories)).

    The counterpart of `predict` for the ordinal fitter, and it inherits three of its four rules:
    the column check is BY NAME against `fit.columns` and a reorder raises `C.SchemaError`; the
    linear predictor is clipped at POLR_ETA_CLIP; and the function names neither treatment nor
    covariate, so a caller supplies the counterfactual meaning by overwriting a column (§7.1).

    The fourth rule is where it differs and the difference is load-bearing. `predict` prepends its
    own intercept because `firth` fits one; this must not, because THE CUTPOINTS ARE THE INTERCEPTS
    and O6 has already established that `X` carries none.

        P(Y <= k | x) = expit(alpha_k + x'beta),   with P(Y <= -1) = 0 and P(Y <= K) = 1 EXACTLY
        P(Y = c_j | x) = P(Y <= j) - P(Y <= j-1)

    The two structural bounds are written as literal 0.0 and 1.0 columns and never as an expit of a
    large number: `expit(POLR_ETA_CLIP)` is 1.0 in float64 but that is an accident of the clip, and
    a row's probabilities must sum to exactly 1.0 by CONSTRUCTION for §20.5's assertion to be about
    the arithmetic rather than about the clip.
    """
```

**Measured: the rows sum to 1.0 with a worst absolute error of 6.7e-16 across every arm of every one
of 2000 replicates**, and every value is in [0, 1]. That is the differencing of a monotone sequence
bracketed by exact 0 and exact 1, so it is a property of the construction and not of this cohort —
which is why §20.5 asserts it on a constructed fit as well as on the workbook.

### 6.2 The returned level set is the FITTED one, and 4.15% of replicates have a smaller one

`polr` collapses the response to the categories carrying positive weight before it fits anything, and
`PolrFit.categories` reports what was fitted. So `ordinal_probabilities` returns
`len(fit.categories)` columns — **not** `len(C.MRS_LEVELS)` — and on a frame missing an mRS level it
returns six columns where the declared level set has seven.

This is not hypothetical here and it is the single most consequential measurement in this document:

| | estimation population | across 2000 replicates |
|---|---|---|
| mRS levels occupied | all 7 | 6 in **4.15%** (83 of 2000) |
| cutpoints | 6 | 5 in 4.15% |
| which level is lost | — | **mRS 5, in all 83** |

The cause is in the point estimate's own counts: mRS 5 carries **3 patients of 104** — 2 at HUG, 1 at
Lugano — and a resample stratified by centre must miss all three, which it does at the measured rate.
Stage 10 §7.4 measured the analogous collapse at **4.6%** on the 92-record one-column [§8] design, so
the two are close; but Stage 10's consequence was a comparability question about `beta` across
replicates, and **this stage's consequence is arithmetic**: a six-column distribution cannot be
averaged with a seven-column one, and a naive implementation would either raise deep inside numpy or
— worse — broadcast and silently misalign every level above the missing one.

### 6.3 The re-expansion rule, and it belongs to the caller

`ordinal_probabilities` returns the fitted level set because that is what the fit is about; **the
caller re-expresses it on `C.MRS_LEVELS`, and an unoccupied level becomes a structural zero.**

    P(Y = j) = 0   for every declared j not in fit.categories

This is prescribed rather than convenient, and the argument is the one Stage 7 §5.3 makes about
`nan`: a level no record occupies is a level for which the fit provides no cutpoint, and the honest
value for it is the one the model implies — zero mass. The alternative of interpolating a cutpoint is
inventing a parameter; the alternative of returning a shorter vector propagates the misalignment into
the average.

**Three consequences, each specified:**

1. **The distribution still sums to 1** — adding zeros changes no sum — so §20.5's assertion holds on
   collapsed replicates too. Measured across all 2000: worst 6.7e-16.
2. **`RD_k` is still defined for every `k` in `MRS_THRESHOLDS`**, because a structural zero leaves the
   cumulative sequence monotone. Measured: monotone in every arm of every surviving replicate.
3. **`RD_4` and `RD_5` coincide in a replicate that lost level 5**, exactly and by construction. This
   is a real narrowing of the sampling distribution of `RD_5`, it is not a bug, and **it is counted
   and reported**: the collapse count is a Stage 12 diagnostic and Stage 14 prints it beside the
   `RD_5` interval, in the same way Stage 10 §7.4's cutpoint distribution reaches a reader through
   [§16]'s constant-shift statement.

---

## 7. The standardisation [§14a]

### 7.1 The duplication is on the DESIGN, and never on the frame

[§14a] says *"duplicate every patient with `A = 1` and `A = 0`"*. The obvious implementation — copy
the frame, set the treatment column, re-run `model.design` — **is wrong, and it fails silently**:
`design` drops constant columns, so a frame in which every patient has `A = 1` yields a design with
**no treatment column at all**, `polr` fits the remaining nine, and the "counterfactual" prediction is
made under a model with no exposure in it. Nothing raises. T1 would catch it only if the
counterfactual design were the one fitted, and it is not.

So the rule is:

    X, dropped = model.design(pop, cols)      # ONCE, on the OBSERVED frame
    fit        = model.polr(X, y)             # ONCE
    for a in (1, 0):
        Xa = X.copy(); Xa[C.TREATMENT] = a    # overwrite the COLUMN of the fitted design
        p_a = model.ordinal_probabilities(fit, Xa)

**This is the pattern `outcome._counterfactuals` (`outcome.py:1087-1104`) already uses** for the Stage
9 binary nuisance models, and Stage 9 chose it for the same reason: `model.predict`'s by-name column
check passes because the columns are the fitted ones, and the overwrite cannot change which columns
exist. `ordinal_probabilities` inherits that check (§6.1), so the two counterfactual designs are
verified identical in shape to the fitted one at every call.

The arm codes are `max(C.TREATMENT_LABELS)` and `min(C.TREATMENT_LABELS)`, taken from the registry
exactly as `outcome.py:118-119` takes them, and never written as 1 and 0.

### 7.2 The averaging, and its denominator is [§11]'s

    P_hat(Y = j | A = a) = (1/N) * sum_i P(Y = j | A = a, X_i)

over the `N = n_average` records of the averaging population, unweighted. `Standardisation.n_average`
and `n_fit` are both carried because §11's sensitivity makes them differ, and [§11] requires the
denominator reported for every estimate.

**Patients at never-IVT centres contribute their covariates and their share of the target
population**, which is [§14a]'s own sentence and the whole point: their IVT counterfactual is
predicted from relationships learned where IVT was observed. Nothing in the code distinguishes them —
`center` is not in the model — and that is what §10's support check exists to make visible rather
than to fix.

### 7.3 The outputs, and two of the four are not new numbers

[§14a] asks for the two standardised distributions and, from them, the cumulative `RD_k`, the
standardised mRS 0–2 risk difference, and the standardised mortality difference. Written out:

```
   RD_k       =  P(Y <= k | A=1) - P(Y <= k | A=0)          k in MRS_THRESHOLDS   (6)
   mrs_0_2    =  P(Y <= 2 | A=1) - P(Y <= 2 | A=0)                       == RD_2
   mortality  =  P(Y  = 6 | A=1) - P(Y  = 6 | A=0)                       == -RD_5
```

The second identity is one line of algebra — `P(Y=6) = 1 − P(Y<=5)` in both arms, and the two 1s
cancel — and the first is not even that. **Measured on the estimation population: `|mrs_0_2 − RD_2|`
is 0.0 exactly and `|mortality + RD_5|` is 8.3e-17.**

They are **reported under [§14a]'s names, computed from the distributions, and asserted equal to their
`RD` counterparts** (§20.6). They are not computed a second way. The rule is [§16]'s: a statement
derived from the result object cannot disagree with the table after an edit, and two independent
computations of one number are two things that can.

The orientation is Stage 8 §7's, unchanged: `arm[1] − arm[0]`, so a positive `RD_k` favours bridging
on the functional scale and a positive `mortality` is **worse**. That sign flip between the two
reported quantities is exactly why `mortality` is named and not left as `−RD_5` for a reader to
negate.

### 7.4 The mortality interval is NOT the reflection of `RD_5`'s, and Stage 10 §8.2's pin is why

Given §7.3, the natural implementation is to compute one set of draws for `RD_5` and obtain the
mortality interval by swapping and negating its limits. **That is wrong, and this document measured
how wrong.**

`C.PERCENTILE_METHOD` is pinned to `"inverted_cdf"` — the empirical-CDF quantile, the smallest order
statistic whose cumulative proportion reaches the level. Stage 10 §8.2 pinned it because it is the
same object `Pr(beta* <= 0)` is computed from, which is what makes [§10]'s p-value agree with its
interval *by construction*. It is **not** a symmetric quantile definition: the lower limit of `−X` is
not the negated upper limit of `X`, because the order statistic reached from below at level `q` is not
the mirror of the one reached from below at `1 − q`.

Measured over the 2000 surviving draws:

| | `|lo(−X) + hi(X)|` | `|hi(−X) + lo(X)|` |
|---|---|---|
| `PERCENTILE_METHOD = "inverted_cdf"` (pinned) | **9.8e-5** | **1.7e-3** |
| numpy's default `"linear"` | 1.1e-16 | 5.6e-16 |

**1.7e-3 on a risk-difference scale is 0.17 percentage points**, which is a difference a manuscript
prints. So the two intervals are genuinely different objects and the specification is:

- **`mortality` carries its own draws and its own `percentile_ci` call.** The draw vector is the
  negation of `rd_5`'s draw-for-draw — measured agreement 8.9e-16 — and the *interval* is computed
  from it independently.
- **`test_standardise.py` §20.7 asserts the two intervals are NOT related by reflection**, at the
  boundary this measurement found, so a later "simplification" that derives one from the other fails a
  test rather than shifting a limit by 0.17 points.

This is the second time Stage 10 §8.2's pin has had a downstream consequence its own document could
not have anticipated — Stage 11 §7.2 was the first, where the pin is what licenses reporting an
odds-ratio interval at all — and both are recorded in §21 as evidence that the pin is load-bearing
rather than a tolerance.

---

## 8. `exp(beta)` is conditional, and the guard is structural rather than a label

The roadmap's Guard for this stage is one sentence: *"`exp(β)` here is a conditional odds ratio. It
may be emitted as a model parameter and must never be labelled as the standardised marginal effect."*
[§14a] says the same at more length and adds the reason — a conditional common odds ratio is not
generally equal to a marginal one in the standardised population, and the non-collapsibility gap is
not small at this covariate count.

A comment saying so would be a label, and labels are what [§16] amendment DECISION 5 exists because of.
The enforcement is therefore three structural things:

1. **There is no field called `odds_ratio` on any Stage 12 record** (§3.1). `Primary.odds_ratio` is
   the [§8] marginal quantity; a formatter written against one record and pointed at the other would
   get an `AttributeError` rather than a number.
2. **`Standardisation.measure` carries the guard string as data**, not as a docstring, so Stage 14
   prints it from the record. It is a `Final` constant with one value, `_CONDITIONAL`, and §20.8
   asserts every constructed `Standardisation` carries it — a second value can only arrive with a
   [§14] amendment behind it.
3. **The estimand key is `beta` and not `odds_ratio`**, and the interval is on `beta`. `exp` of a
   percentile limit *is* the percentile limit of `exp` under `PERCENTILE_METHOD` — Stage 11 §7.2
   established that and it holds here — so a Stage 14 reader wanting the ratio scale exponentiates the
   limits and gets the right interval. The key stays on the log scale so that the record never holds a
   number a careless read turns into a reported effect.

**And the reported [§14a] effect is the set of risk differences, not `beta`.** That is [§14a]'s
choice, made for exactly this reason: the standardised quantities are marginal by construction, and
they are on the absolute scale where non-collapsibility does not arise.

---

## 9. Separation, and this stage's guard is at key granularity

### 9.1 The asymmetry Stage 8 does not have

Stage 8 §6 established that a separated proportional-odds fit **converges and returns**: measured
`beta` 36.4, `exp(beta)` 6.5e15, every safeguard counter at zero. `POLR_MAX_ABS_BETA = 14.0` exists
because that fit's *reported quantity* is `exp(beta)`, and 6.5e15 is a number no manuscript could
print. Stage 8 §6.3 calibrated 14.0 from a measured empty band on a **one-column** design.

Here the reported quantities are **averaged probabilities**. A separated fit drives `beta` to
infinity, but `expit` of anything is in [0, 1], so every standardised probability and every risk
difference stays in range no matter how degenerate the fit is. **Measured over 2000 replicates:**

| | min | max | finite | in range |
|---|---|---|---|---|
| `rd_0` | −0.2081 | +0.1724 | all | all in [−1, 1] |
| `rd_2` | −0.3235 | +0.2743 | all | all in [−1, 1] |
| `rd_5` | −0.2460 | +0.1857 | all | all in [−1, 1] |
| `dist1_0` | 0.0157 | 0.3397 | all | all in [0, 1] |
| `dist1_6` | 0.0937 | 0.4873 | all | all in [0, 1] |

So the bound is not protecting the reported estimand here; it is protecting `beta`, which this stage
reports as a **model parameter** (§8) and nothing else.

**The guard is therefore at key granularity**: a fit reaching `POLR_MAX_ABS_BETA` on the treatment
coefficient drops the `beta` key and **keeps** the standardised keys, which are unaffected and
correct. Stage 10 §7 already supports per-key failure — its `_replicate` costs one binary outcome's
three keys without touching the others — so this needs no new mechanism, only the declaration of a
third granularity.

**The departure from Stage 8 §6 is declared and it is not a weakening.** Stage 8 drops the whole
replicate because the whole replicate's reported content is `beta` and the six `RD_k` derived from the
*weighted empirical distributions*, which a separated fit does not touch either — Stage 8's coupling
is a simplification it could afford at one estimand. Here the coupling would discard 22 valid keys to
suppress one invalid one, which is a bias in the surviving-draw set, not a safeguard.

### 9.2 And on this workbook it never fires, which closes a `TODOS` item negatively

`TODOS.md` files that `POLR_MAX_ABS_BETA` *"is calibrated on a one-column design while [§14a]'s
standardisation model carries eight covariates … that one bites at Stage 12."* Measured over 2000
replicates of the ten-column design:

```
   |beta_treatment|   min 0.0000   median 0.3219   max 2.0016
   max|gamma|         min 0.2638   median 0.9798   max 3.2119
   POLR_MAX_ABS_BETA  14.0
   replicates reaching it, on either                                    0
```

**The bound is nowhere near.** The largest coefficient anywhere in 2000 fits is 3.21 against a bound
of 14.0, and Stage 8 §6.3's measured empty band ran from 8.79 to 18.81 — so this design's fits sit
entirely below the band's *lower* edge and the calibration question does not arise. The item closes
with a negative answer and a re-stated trigger (§21).

**`gamma` is not bounded and that is deliberate**, following Stage 11 §9's G9 rule: a nuisance
coefficient may grow legitimately when a covariate is nearly collinear in a resample, and bounding it
would drop replicates for a reason the reported quantity does not care about. The measurement above
is a diagnostic, recorded in the replicate audit entry, not a guard.

---

## 10. The support check [§14a]

[§14a]'s framing is the one worth keeping in the code's comments: *"Standardisation avoids infinite
weights, not extrapolation."* There is no weight to blow up here, so nothing in the arithmetic
complains when a prediction is made far outside the region where any comparison was observed; the
support check is the only thing that makes it visible.

### 10.1 The box, and the definition is measured to be inert between two readings

[§14a] asks for *"the treated-arm distribution of each continuous covariate"* and *"the proportion of
never-IVT-centre patients outside it"*. Two decisions:

**The support is the per-covariate treated-arm `[min, max]` box.** Not a percentile trim, not a convex
hull. The box is the weakest possible reading of "outside the treated support" — a patient outside it
is outside on a *single* covariate's observed range, which no smoothing assumption can be argued to
cover — and it is the only definition under which "outside" needs no tuning parameter. A convex hull
in seven dimensions on 39 treated patients is a set that almost every point is outside of, which
would make the diagnostic report 100% and mean nothing.

**"Continuous" is resolved to "not in `C.CATEGORICAL`", and the choice is measured to be inert.** The
non-categorical covariates in `STANDARDISATION_COVARIATES` are `age`, `sex`, `prestroke_mrs`,
`nihss_baseline`, `core_ml`, `tmax6_ml`, `atrial_fib` — of which `sex` and `atrial_fib` are binary and
`prestroke_mrs` is an integer 0–5 modelled linearly. A narrow reading would box only `age`,
`nihss_baseline`, `core_ml` and `tmax6_ml`. Measured, the two readings give the **same set**:

| covariate | treated min | treated max | records outside | of which at USZ |
|---|---|---|---|---|
| `age` | 52 | 97 | 4 | 0 |
| `sex` | 0 | 1 | 0 | 0 |
| `prestroke_mrs` | 0 | 3 | 0 | 0 |
| `nihss_baseline` | 3 | 25 | 3 | 1 |
| `core_ml` | 0 | 62 | 4 | 2 |
| `tmax6_ml` | 13 | 401 | 0 | 0 |
| `atrial_fib` | 0 | 1 | 0 | 0 |

`sex`, `prestroke_mrs` and `atrial_fib` exclude nobody — the treated arm covers both levels of each
binary and reaches `prestroke_mrs = 3` — so the wider reading is chosen **because it cannot be wrong
on a workbook where the narrow one is right**, and the inertness is asserted rather than assumed
(§20.10). `SUPPORT_COVARIATES` is a computed view in `config.py`, not a fourth literal list (§18).

The box is over **declared factor levels for the categorical covariate too**, in the weak sense that
`onset_type`'s three declared levels are all present in the treated arm — measured — so no patient is
outside on it. That is recorded in the audit entry rather than folded into `inside`, because a
categorical level unseen in the treated arm is a *different* kind of unsupported prediction — its
coefficient does not exist — and it would be caught by `design` dropping the column and reported
through `Standardisation.dropped`.

### 10.2 Who falls outside, and the never-IVT set is computed and never named

```
   inside the box on all seven covariates          93 of 104
   outside                                         11
     by centre    HUG 4, CHUV 3, USZ 3, Lugano 1
     by arm       11 EVT alone, 0 bridging
   never-IVT-centre patients                       12
     of which outside the box                       3   (25.0%)   <- [§14a]'s requested proportion
```

**Everyone outside the box is a control patient**, which is arithmetic rather than a finding — the box
is the treated arm's own range, so no treated patient can be outside it — and it is stated because a
reader meeting "11 outside" will otherwise wonder about the split.

**The never-IVT set is `tuple(c for c in C.CENTER_ORDER if c not in cohort.treating_centres(df))`,
computed at the one place that needs it.** `cohort.py:112-119` wrote that expression into its own
docstring as the thing Stage 12 should use, and Stage 5 §4.2's argument is why: a centre's treatment
availability is a property of the data, and a workbook in which USZ starts administering IVT must
change this diagnostic rather than require someone to remember a list. `EXPECTED_NEVER_IVT` is
asserted against the computed set and is never the operative rule — measured, the two agree at
`('USZ',)`.

### 10.3 The baseline table, and its centre rows are the grouping variable

[§14a] asks for *"a baseline table comparing IVT-treated with never-IVT-centre patients."* Two
disjoint groups: the 39 treated patients anywhere, and the 12 eligible patients at USZ. The 53
control patients at treating centres are in **neither** — they are not what [§14a] asked to compare —
and their absence is printed in the table's caption, with all three counts, so the denominator is
never inferred.

**It calls `balance.smd` directly, with unit weights**, which is the note Stage 7 §5.5 wrote to this
stage: *"[§14a]'s support check is a standardised mean difference over a different population with
unit weights, and a second implementation there would be a second definition of the yardstick."* The
rows range over `C.BALANCE_SET` — the [§6] confounders plus the five balance-only covariates —
because [§9]'s rule that balance is judged against the **full** confounder set is about what a reader
must be shown, and it does not stop applying because the comparison is not a treatment contrast.

`smd`'s five documented `nan` routes (Stage 7 §5.3) apply unchanged and the caller owns the mask. One
of them **fires here by construction and is reported rather than suppressed**: `center = USZ` is
constant at 1 in one group and 0 in the other, so the pooled SD is exactly zero with the groups
differing, which is Stage 7's branch 4 and returns `nan`. Measured: `center = USZ` is `nan`.

**The four `center` rows are the grouping variable and are labelled as such in the table**, not read
as imbalance. This is the one place the table would mislead: the never-IVT group *is* USZ, so the
centre rows report a definition rather than a difference. They stay in the table — dropping a declared
row is what Stage 7 §4.1's role column exists to avoid, and a reader must see that the comparison is
between disjoint centre sets — and they carry a role of `"grouping"` rather than a [§6] role.

The remaining rows are the ones [§14a] wants read, and **they are not small**: on `BALANCE_SET`
excluding the centre rows, the largest |SMD| between the two groups is on `prestroke_mrs` and
`hypertension`, both above 0.8, against `SMD_THRESHOLD = 0.10`. That is the transport assumption's
size stated as a number, and §17.2 hands it to Stage 14 as something that must appear beside the
[§14a] estimate rather than in a supplement — [§9]'s 2026-08-24 amendment applied to this stage's own
diagnostic.

---

## 11. Sensitivity 1 — restrict the standardisation population to the treated support

[§14a]: *"Sensitivity: restrict the standardisation population to patients inside the treated
support."*

**The averaging population is restricted; the fit is not.** `all_centre` takes `over=` and passes it
to the averaging step only, so `n_fit` stays 104 and `n_average` becomes 93.

**The alternative was considered and is declined with its argument recorded**, in the manner of Stage
5 §15 and Stage 11 §17. Refitting on the in-support subset is a defensible reading — the 11 excluded
patients do influence `gamma`, and a reader could ask for an estimate untouched by them. It is
declined because it changes **two** things at once: the target population *and* the fitted model, so a
difference between the arms could not be attributed to either. [§14a]'s sentence names the
standardisation population, which is the averaging set; the one-thing-at-a-time reading is also the
literal one. §22 records that Stage 12 does not decide whether a refit-on-support arm should exist —
that is a [§14] amendment.

**The box is recomputed inside every replicate and is not held at the point estimate's.** The support
is a function of the treated arm's observed range, so it is a statistic; holding it fixed would
condition the bootstrap on one realisation of it and understate the interval. Measured across 2000
replicates the in-support averaging population ranges **69 to 102 with a median of 91**, against 93 at
the point estimate — a spread wide enough that the choice is not cosmetic.

This arm's keys carry the `support.` prefix and there are 23 of them (§3.3). `Standardisation.population`
carries `"treated support"` so the record says which it is without a caller tracking the prefix.

---

## 12. Sensitivity 2 — the random centre intercept, and it is a new estimator

[§14a]: *"Add a random centre intercept `b_c ~ N(0, σ²_C)` with **one common treatment effect**; no
treatment-by-centre interaction and no random slope, neither being identifiable in any useful sense
with four centres."*

This is the only new estimator in Stage 12, and **nothing in the environment can fit it.**
`model.polr` is fixed-effects; `statsmodels` has `MixedLM` for the linear case and Bayesian mixed GLMs
for binomial and Poisson, and **no ordinal mixed model at any API**; `scipy` is test-only by policy.
So it is written, and this section is its specification.

### 12.1 The model

    logit P(Y <= k | x, b_c)  =  alpha_k + x'beta + b_c ,      b_c ~ N(0, sigma^2)

    log L_c  =  log INTEGRAL  phi(b; 0, sigma^2) * PROD_i P(Y_i = y_i | x_i, b)  db

    log L    =  SUM_c log L_c

One common `beta` across centres — no `A × centre` term and no random slope, both forbidden by
[§14a] and by [§15]. `x` is the same design `polr` fits (§5.3), including the treatment column; the
grouping variable is `center`, which is **not** in the design and could not be: a fixed `center`
effect and a random one are the same parameter twice.

**The parameter vector is `(alpha, beta, log sigma)`** — `log sigma` and not `sigma`, so the positivity
constraint is structural. §12.5 is what that costs at the boundary.

### 12.2 The quadrature must be ADAPTIVE, and this is measured rather than assumed

The integral is one-dimensional per centre and Gauss–Hermite is the standard tool. **Plain
(non-adaptive) Gauss–Hermite places its nodes at `b = sigma*sqrt(2)*t_q`, which is right when the
integrand is close to the prior and increasingly wrong as the data pull the posterior away from it.**
Adaptive Gauss–Hermite finds each centre's conditional mode `b_hat_c` and curvature `tau_c` and places
the nodes at `b_hat_c + sqrt(2)*tau_c*t_q`.

Measured on the estimation population, absolute error in the negative log-likelihood against a
121-node adaptive reference:

| σ | nodes | non-adaptive error | adaptive error |
|---|---|---|---|
| 0.5 | 3 | 2.8e-2 | 9.8e-5 |
| 0.5 | 7 | 6.7e-3 | 1.5e-8 |
| 0.5 | 31 | 4.4e-7 | 0 |
| 1.0 | 7 | 9.6e-2 | 2.0e-7 |
| 1.0 | 31 | 1.9e-3 | 0 |
| 2.0 | 9 | 4.6e-1 | 1.0e-7 |
| 2.0 | 31 | 4.5e-2 | 0 |
| 3.0 | 7 | 1.1e0 | 2.4e-6 |
| 3.0 | 31 | **4.8e-1** | 0 |

**Three things this settles, none of which is a preference:**

1. **Non-adaptive quadrature does not converge usefully here.** At σ = 3 it is off by 0.48 log-
   likelihood units with 31 nodes — a deviance error of nearly 1, on an objective whose whole
   comparison against the pooled model is a deviance difference below 2.
2. **It is non-monotone in the node count.** At σ = 2 the error at 9 nodes (4.6e-1) is *worse* than at
   7 (3.5e-1) and worse than at 5 (4.0e-1). A node count chosen by "increase it until the answer stops
   moving" would have stopped at the wrong place, and this document's first prototype did.
3. **Adaptive quadrature's accuracy is essentially independent of σ.** Seven adaptive nodes beat
   thirty-one non-adaptive ones by five orders of magnitude at *every* σ tested.

**And point 3 is not a refinement, it is what makes the arm computable**, because §12.5 measures that
**24.0% of replicates fit σ̂ above 1** — squarely in the region where the non-adaptive rule is not
converged. A fixed non-adaptive node count would give a quarter of the bootstrap draws a different
numerical answer from the other three quarters, and the difference would be invisible.

### 12.3 The node count is a prespecified constant, and it changes an ANSWER

    POLR_RI_NODES = 11

`POLR_RI_NODES` belongs in `config.py` under `PERCENTILE_METHOD`'s argument, not `POLR_TOL`'s: it is
not a tolerance on how precisely a fixed answer is found, it is **part of the definition of the
objective**, and [§10] refits this model in every one of `N_BOOT` replicates. Measured, refitting at
each node count:

| nodes | σ̂, non-adaptive | σ̂, adaptive |
|---|---|---|
| 3 | — | 0.54236775 |
| 5 | 0.72036516 | 0.54242568 |
| 7 | 0.80328763 | 0.54242682 |
| 9 | 0.83240300 | 0.54242685 |
| 11 | 0.56136873 | 0.54242685 |
| 15 | 0.54029101 | 0.54242685 |
| 21 | 0.54193807 | 0.54242685 |
| 31 | 0.54242685 | 0.54242685 |
| 61 | — | 0.54242685 |

*(σ̂ appears in this table as the quantity whose numerical stability is being demonstrated. Its
reported value, its interval and any reading of between-centre heterogeneity are the log's and the
manuscript's — §4.3.)*

**Non-adaptive quadrature at 5, 7 and 9 nodes gives an answer 33% to 53% too large.** Adaptive
quadrature is stable to eight significant figures from **9** nodes. `POLR_RI_NODES = 11` is chosen one
step above that — the measured cost of the extra two nodes is 0.5 s per 2000-replicate arm — so a
workbook whose posterior is slightly less Gaussian than this one needs no re-tuning. `test_config.py`
asserts it is odd and at least 9 (§20.2); odd because a symmetric rule with an odd node count places a
node at the mode, which is where the mass is.

### 12.4 The optimiser, and its shape is a cost measurement

`polr`'s loop is Newton on the analytic Hessian, and it is excellent there: 5 iterations, 0 halvings
(§5.2). It is not available here — the marginal log-likelihood's Hessian in `(alpha, beta, log sigma)`
has no compact closed form once the quadrature weights depend on the parameters — so the choice is
between a numerical Hessian and a quasi-Newton update. Measured, on the same fit:

| prototype | quadrature | gradient | curvature | iterations | one fit | at `N_BOOT` |
|---|---|---|---|---|---|---|
| C2 | non-adaptive, 31 | numerical | central-difference Hessian | 6 | **39.0 s** | **21.7 h** |
| D2 | non-adaptive, 31 | analytic | BFGS | 25 | 0.82 s | 27.3 min |
| E2 | adaptive, 11 | numerical | BFGS | 25 | 3.91 s | 130 min |
| **F2** | **adaptive, 11** | **analytic** | **BFGS** | **25** | **0.76 s** | **17.5 min** |

**The specification is F2**, and the two orthogonal findings are:

- **The analytic gradient is worth 5×**, and it is what makes the arm affordable at all. A numerical
  gradient is `2 * n_par = 34` objective evaluations; the analytic one is one. It is therefore part of
  the specification and §20.12 asserts it against central differences rather than leaving it to a
  future optimisation.
- **Newton with a numerical Hessian is worth −51×.** It converges in 6 iterations rather than 25, and
  each costs `4 * n_par^2 = 1156` objective evaluations. The iteration count is the tempting number
  and the wrong one.

**The loop is `polr`'s discipline, transplanted**, so a reader moving between the two fitters is not
made to learn a second convention: BFGS inverse-Hessian update; a **relative** trust region at
`POLR_MAX_STEP` for Stage 6 §5.2a's reason (a proportional-odds estimate is equivariant under
rescaling a covariate and an absolute bound is not, so with one, *which* replicates [§10] drops would
depend on how the data was recorded); step-halving to `POLR_MAX_HALVINGS` **per iteration**, raising
`FitError` on exhaustion; and the same two convergence routes in the same order, with `converged_on`
recorded. The iteration cap is its own constant, `POLR_RI_MAX_ITER`, because 25 iterations against
`polr`'s 5 means the two estimators are not on the same scale and one shared cap would be either loose
for `polr` or tight here.

**Two departures from `polr`, both declared:**

1. **`POLR_MAX_HALVINGS` is per iteration and the total is much larger.** Measured on the point
   estimate: 25 iterations, **85 halvings total**, 2 rescales, no single iteration exhausting. Across
   200 replicates the total ranges 75 to 95. BFGS proposes long steps early and the line search
   shortens them; that is the algorithm working, not a warning, and the total is recorded as a
   diagnostic so a future change that makes it grow is visible.
2. **The quadrature nodes are recomputed at the current parameters and then held FIXED for that
   evaluation's gradient**, so the gradient is exact for the frozen-node quadrature sum and neglects
   the derivative of the node positions. Measured against central differences on the true objective:

   | σ | relative error of the frozen-node gradient |
   |---|---|
   | 0.20 | 5.0e-8 |
   | 0.54 | 9.3e-9 |
   | 1.50 | 9.5e-9 |
   | 3.00 | 6.0e-8 |

   Six to eight orders of magnitude below `POLR_SCORE_TOL = 1e-6`, at every σ the bootstrap reaches.
   The approximation is named rather than hidden, and §20.12 pins the measurement so that a change to
   the mode-finding step which made it worse fails a test.

### 12.5 The σ → 0 boundary is reached in 22.0% of replicates, and it needs a floor

[§14a] anticipates the boundary: *"if outcomes truly do not differ by centre given `X` then σ²_C = 0
and it collapses to the pooled model."* It is not a hypothetical. Measured over 200 replicates:

```
   sigma_hat    min 0.000021   median 0.65923   max 1.78734
   below 0.01, the variance collapsing to the pooled model     44 / 200   22.0%
   above 1.0                                                   48 / 200   24.0%
   converged                                                  200 / 200   all on "likelihood"
   iterations   min 23   median 27   max 47      (POLR_RI_MAX_ITER 200)
```

**Under a `log sigma` parametrisation the boundary is at −∞**, and the conditional-mode Newton
(§12.2) evaluates `1/sigma^2`, which **overflows** as the search approaches it — observed directly in
the probe, as an overflow warning in the curvature and a `-inf` in a quadrature weight. The fit still
returned, because the affected nodes are dominated in the log-sum-exp, but "still returned" is not a
specification.

    POLR_RI_SIGMA_FLOOR = 1e-4

- **`log sigma` is clipped at `log(POLR_RI_SIGMA_FLOOR)` at every evaluation**, so `1/sigma^2` is at
  most 1e8 and every quadrature weight is finite.
- **`RIFit.at_floor` records whether the returned `sigma` is at it**, and the fit is **not** a failure
  when it is: σ̂ = 0 is a legitimate answer that [§14a] names, and it means the arm has collapsed to the
  pooled model. Dropping those replicates would select the bootstrap on the value of a variance
  parameter, which is selection on exactly the quantity the arm exists to examine.
- **The `at_floor` count is reported beside the arm's interval**, as a rate. A percentile interval on
  `hier.sigma` whose lower limit sits at a floor is a limit that means "the boundary", not "a value",
  and Stage 14 must print the rate for the limit to be readable (§17.2).
- **1e-4 is a floor on the SD, so 1e-8 on the variance**, five orders of magnitude below the smallest
  σ̂ that is not at the boundary in the measured 200. It separates "collapsed" from "small" without
  being reachable by a genuinely small non-zero variance.

**The floor changes an answer at the boundary and is therefore a `config.py` constant**, prespecified
for `PERCENTILE_METHOD`'s reason, not a local literal.

### 12.6 What is standardised under it — the empirical-Bayes intercepts

A random-intercept model offers two standardisations and they answer different questions:

1. **Marginal** — integrate the prediction over `b ~ N(0, sigma^2)`, giving the distribution for a
   patient at a *randomly drawn* centre.
2. **Conditional on the centre's own intercept** — predict each patient with `b_hat_c`, the
   conditional mode for the centre they are actually at.

**[§14a] prescribes the second, in a sentence that is easy to read past**: *"Patients at a never-IVT
centre inform that centre's intercept through their direct-EVT outcomes while the IVT effect is
borrowed."* A patient can only inform their own centre's intercept if the prediction uses it. Under
the marginal reading, USZ's 12 patients would inform `sigma` and nothing else, and the sentence would
be false.

So `Hierarchical.standardisation` is built from predictions at `alpha_k + x'beta + b_hat_c`, with
`b_hat_c` the conditional mode the adaptive quadrature already computes at every evaluation — it is
free, and `RIFit.b` carries it. The marginal alternative is recorded in §22 as a [§14] question rather
than implemented, because it is a *different estimand* and [§14a] named one.

**And the arm's honest limitation is a measurement, not a caveat.** The conditional-mode precisions,
at the fitted model:

| centre | n | treated | posterior SD of `b_hat` |
|---|---|---|---|
| HUG | 41 | 30 | 0.2517 |
| CHUV | 21 | 7 | 0.3156 |
| Lugano | 30 | 2 | 0.2823 |
| **USZ** | **12** | **0** | **0.3810** |

**USZ's intercept is the least precisely estimated of the four and it is the one the whole arm rests
on**, because it is the centre whose IVT counterfactual is entirely borrowed. Its posterior SD is 51%
larger than HUG's, from 12 patients none of whom received the treatment. That is [§14a]'s *"a centre
variance from four clusters is fragile"* made specific, it is reported in the arm's audit entry, and
Stage 14 prints it beside the arm (§17.2). The intercept **values** are centre-level outcome contrasts
and are the log's (§4.3).

### 12.7 The oracle: `ordinal::clmm`

`polr_ri` is a new maximiser for a model with a published reference implementation, so it gets an
oracle on Stage 8 §18b's pattern. `tests/reference/polr_clm.R` already uses the `ordinal` package;
`clmm` is in the same package, so **no new R dependency is added** and the existing
`ORDINAL_PACKAGES = ("ordinal",)` gate covers it.

```r
ordinal::clmm(ordered(y) ~ . - center + (1 | center), data = frame, link = "logit", nAGQ = 11)
```

Four things the oracle must assert, and the first is the one that makes it worth having:

1. **The coefficients come back NEGATED and the thresholds do not**, exactly as `polr_clm.R`'s banner
   documents for `clm`, because `clmm` uses `zeta_k − x'beta`. **`sigma` does not come back negated
   either** — it is a scale — so the oracle asserts a three-way asymmetry and never that "the fits
   agree", which would pass on the thresholds alone.
2. **`nAGQ` is passed and set to `POLR_RI_NODES`.** `clmm`'s default is `nAGQ = 1`, the Laplace
   approximation, which §12.2's table shows is a different objective; comparing against it would
   measure the approximation and not the maximiser.
3. **`sigma` is compared on the SD scale**, `clmm` reporting the standard deviation of the random
   effect, and the tolerance is separate from the coefficients' because a variance parameter at this
   cluster count is the least well determined thing in the fit.
4. **The gate is loud.** `pytest -rs` reads as an instruction, as Stage 6 §16b.2 established, and this
   machine's R currently **segfaults on startup** for the `uname`/PATH reason `test_reference_r.py`'s
   banner documents at length — so `_r_environment` is what this oracle must inherit, and it must not
   grow its own subprocess call.

### 12.8 The cost, and it is 95% of Stage 12's bootstrap

```
   one fit, point estimate                                        0.76 s
   200 replicates, end to end                                   105.0 s     525 ms each
   projected at N_BOOT = 2000                                    17.5 min
   the same arm under the first prototype (C2)                   21.7 hours
```

The projection is from a measured 200-replicate run of the full body — `model.design`, `model.polr`
for the start values, then `polr_ri` — and not from the point fit's 0.76 s, because the start values
matter: `polr_ri` starts from the pooled fit's `alpha` and `beta` with `sigma = 0.5`, which is a
nested model's exact maximiser in every coordinate but one, and it is why 25 iterations suffice.
`polr`'s docstring makes the same argument for its own start values.

**`TODOS.md` records that Stage 12's bootstrap is now the sole trigger for parallelising the replicate
loop**, with a seed *sequence* — `SeedSequence.spawn(N_BOOT)` — as the right design. 18.4 minutes on a
130–145 s pipeline is a real change in its character, and §21 keeps the item open with the trigger now
fired rather than pending.

---

## 13. Inference [§10, §14]

*"Patient-level bootstrap stratified by centre, as §10: resample, refit the ordinal model, predict
every resampled patient under both regimes, average, recompute the risk differences, percentile
intervals."*

### 13.1 All four strata are resampled, and Stage 10 §12.2 left this decision here

Stage 10 §12.2, in as many words: *"[§14a] runs over all eligible patients at all centres, which
includes USZ — a centre with no treated patients — so a stratum can exist there that contributes no
contrast, and whether that is resampled or held fixed is [§14a]'s question and not this document's."*

**It is resampled. All four strata, USZ included, with `bootstrap.resample` unchanged.**

The argument is about what the estimand is. [§14a]'s target is `E[Y^IVT − Y^NoIVT | eligible]`
averaged over **the full eligible cohort at all centres**, and USZ's 12 patients are members of that
population: they contribute their covariates to the average (§7.2) and their observed outcomes to the
fit. The uncertainty in an average over a sample includes the uncertainty in *which* sample was
drawn. Holding USZ fixed would produce an interval conditional on those 12 covariate vectors, which is
a different and narrower inferential statement than the one [§14a] makes — and it would be inconsistent
with the other three strata, which are resampled for exactly the reason USZ would not be.

**The counter-argument is recorded because it is not empty.** USZ contributes no contrast, so
resampling it adds variance to the *target population* without adding information about the *effect*;
a reader might reasonably want the effect's uncertainty separated from the population's. That is a
real distinction and it is not [§14a]'s: [§14a] defines its estimand as an average over a population
and asks for a bootstrap of it. §22 records that Stage 12 does not decide whether a
population-conditional variant should exist.

**Nothing in `bootstrap.resample` changes.** It is general in frame and stratum by Stage 10 §4's
design — *"a function that hard-coded `center` would still be right there while being right for the
wrong reason"* — and it groups with `sort=True, observed=True`, so a fourth stratum appears simply
because the frame has one. `C.BOOT_STRATUM` is `"center"` and this stage passes it unchanged.

**One property is inherited and re-measured**: the same seed, the same frame and the same stratum give
the byte-identical sequence of drawn frames regardless of what the body does or raises (Stage 10
§15.2, re-confirmed across three bodies at Stage 11 §14). Stage 12's body is the fourth and §20.13
asserts it again here, because this stage's body is the first that runs **three arms** off one draw.

### 13.2 The replicate body, and it has three failure groups

One `bootstrap.replicates` call. One drawn frame feeds all three arms, so the three are comparable
draw for draw and the resample stream is paid for once.

```
   body(draw):
     X, dropped = model.design(draw, cols)
     T1         -> FitError; costs EVERY key                    the design has no exposure
     fit        = model.polr(X, y)
     FitError   -> costs every key EXCEPT hier.*                the pooled fit is both arms' fit
     G-bound    -> costs `beta` and `support.beta` ONLY         §9.1, key granularity
     pooled     = standardise(...)                              23 keys
     support    = standardise(..., over=in_support(draw))       23 keys, box RECOMPUTED  §11
     ri         = model.polr_ri(X, y, draw[BOOT_STRATUM])
     FitError   -> costs hier.* ONLY                            24 keys
     hier       = standardise(..., b_hat)                       §12.6
```

**The granularity is Stage 10 §7's, extended by one level and not reinvented.** Stage 10's
`_replicate` already has three: a propensity failure costs every key, a primary failure costs `beta`
and the six `rd_k` together, a per-outcome failure costs that outcome's three. Stage 12's are the same
idea over its own dependency graph, and the rule is the one Stage 10 states — a failure costs exactly
the keys that could not be computed and no others — with §9.1's key-granularity case as the one
addition.

**`bootstrap.bucket` classifies and `bootstrap.collect` reconciles**, both public at Stage 11 §5.4.
Nothing here re-implements either, which is what Stage 9 §14 names as how a taxonomy gets violated by
accident. Three tokens are added to `C.FAILURE_BUCKETS` (§14).

**`Replicate.n_alpha` and `.polr_iterations` carry the POOLED fit's**, which is the model [§14a] names,
and `polr_ri`'s diagnostics are **not** forced into them. `TODOS.md` files that `Replicate` cannot
describe a multi-fit body and names [§14a] as the candidate trigger; §18 records that the trigger has
fired and that the fix is still declined, because widening two scalars to dicts touches every Stage 10
assertion over them to gain what `hier.sigma` already gives as a first-class key (§3.3).

### 13.3 The interval loop becomes `bootstrap.intervals`, and this stage is the third caller

`TODOS.md`: *"the same five-line loop over `Draws` — `percentile_ci`, the `ci_min_draws` floor, the
p-rule — exists in two places … **Trigger.** A third caller — Stage 12's standardisation bootstrap is
the candidate."*

**It is the third caller, so the item closes here rather than acquiring a fourth copy.**

```python
def intervals(draws: dict[str, Draws], tested: Callable[[str], bool]) -> dict[str, Interval]:
    """The [§10] percentile interval for every key with enough surviving draws, and a p for those
    `tested` returns True for. The p RULE is the parameter, because it is the only thing the three
    callers disagree about: `run` uses `_tested`, Stage 11's arm `lambda key: False`, its subgroups
    `_tested_subgroup`, and Stage 12 `lambda key: False` (§3.3)."""
```

Stage 11's two callers move onto it in the same change. `bootstrap.py`'s public surface goes from
eight names to nine, and `run`'s own loop is replaced by a call — which is a change to the function
whose twenty-six intervals Stage 10's Definition of done asserts, so **§20.11 asserts those
twenty-six are byte-identical before and after**, which is the guard Stage 11 §16 said the refactor
would need.

`C.ci_min_draws(C.CI_LEVEL)` is 40 and `N_BOOT` is 2000. Measured, every Stage 12 key has **2000**
draws, so the floor is nowhere near and no key is dropped for thinness.

### 13.4 The drop rate is zero, and that closes a second `TODOS` item negatively

```
   replicates attempted                                    2000
   survived                                                2000       100.00%
   failure buckets                                         none
   T1 (design lost the exposure)                              0
   separation (POLR_MAX_ABS_BETA)                             0        §9.2
   nonconvergence, pooled                                     0
   nonconvergence, polr_ri                                    0        200/200 measured, §12.5
   design columns dropped                                     0        100.00% empty
```

**Nothing fails.** That is worth stating plainly against the run of this pipeline's history: Stage 10
measured a zero drop rate on the [§8] estimator, Stage 11 measured **12.5%** on the `unknown_onset`
subgroup and called it the largest the pipeline had produced. [§14a]'s design is wider (ten columns
against one) but its population is larger (104 against 92) and — the reason that matters — **it is
unweighted**, so no replicate can concentrate its effective sample on a handful of records the way an
overlap-weighted fit can.

**And this closes the `TODOS` item about the convergence tolerance's scale, negatively.** The item
records that `POLR_TOL = 1e-8` is an **absolute** tolerance on a log-likelihood whose scale is `Σw`,
and files *"a Stage 12 design where `Σw` is materially larger — [§14a] uses unit weights over a bigger
population, so its objective is on a different scale again"* as its trigger. Measured: `Σw = n = 104`
here against roughly 23 for the [§7] overlap-weighted cohort, so the objective **is** on a different
scale — and the fit converges in **4 to 6 iterations in every one of 2000 replicates, all on the
likelihood criterion**, with `first_step_norm` between 1.01 and 22.8. The trigger fired and the
consequence did not occur. §21 keeps the item open, narrowed, because it remains a real property of an
absolute tolerance and this is evidence rather than proof.

---

## 14. The failure taxonomy — the T series, and there are eleven

`C`, `D`, `E`, `F`, `G`, `H`, `O`, `P`, `R`, `S` and `B` are taken by Stages 2–11. Stage 12's
identifiers are **T**, and the split between `C.SchemaError` and `model.FitError` follows the rule
Stage 9 §12.3 established and Stage 10 §7 depends on: **a `FitError` is droppable-and-countable inside
a replicate; a `SchemaError` is a contract break and is never caught.**

| id | raised by | class | meaning |
|---|---|---|---|
| **T1** | `_assert_exposure_survived` | `FitError` | the [§14a] design lost the treatment column (§5.3) |
| **T2** | `polr_ri` | `FitError` | step-halving exhausted at an iteration (§12.4) |
| **T3** | `polr_ri` | `FitError` | no convergence in `POLR_RI_MAX_ITER` (§12.4) |
| **T4** | `polr_ri` | `FitError` | the grouping variable has fewer than two groups with records |
| **T5** | `ordinal_probabilities` | `SchemaError` | the columns are not the fit's, or are reordered (§6.1) |
| **T6** | `population` | `SchemaError` | the frame carries no `eligibility` column — `classify` was not run |
| **T7** | `population` | `SchemaError` | the [§14a] population is empty, or one arm of it is |
| **T8** | `all_centre` | `SchemaError` | `over=` is not boolean, or is not indexed like the population |
| **T9** | `support` | `SchemaError` | `treating_centres`' complement is not `EXPECTED_NEVER_IVT` (§10.2) |
| **T10** | `_assert_distribution` | `SchemaError` | a standardised distribution does not sum to 1, or is not monotone (§20.5) |
| **T11** | `inference` | `SchemaError` | a key's draws disagree with its point estimate's identity (§7.3) |

**T2, T3 and T4 lead with `polr_ri:` and that token goes in `C.FAILURE_BUCKETS`**, mapping to
`nonconvergence`, `nonconvergence` and `degenerate_design`. `test_bootstrap.py`'s scan asserts every
`FitError` raise-site token in `model.py` and `outcome.py` is a `FAILURE_BUCKETS` key; **the scan
scope grows to include `standardise.py`** (§18), so T1's token is covered by the same mechanism and a
reworded message is a test failure rather than a counter that silently reads zero.

**T9 is a `SchemaError` and not a finding.** A workbook in which USZ acquires a treated patient, or in
which a fifth centre appears, changes what [§14a]'s support check is *about* — the never-IVT set is
the thing the check is defined against — and the pipeline must stop rather than report a diagnostic
whose subject moved. `EXPECTED_NEVER_IVT`'s own comment (`config.py:509-511`) describes exactly this
role: an assertion target, never an operative rule.

**T10 and T11 are structural self-checks and both are unreachable on this workbook** — measured worst
`|sum − 1|` of 6.7e-16 across every arm of every replicate, and the two identities agreeing to 8.3e-17
— which is recorded rather than treated as a reason to omit them. They are the [§14a] Accept-when
conditions expressed as code, and §20.5 is the argument for asserting a property that cannot currently
fail.

---

## 15. The audit entries, and there are ten

All under existing `data.KINDS`; **no kind is added**, which is the first stage since Stage 5 of which
that is worth saying, and the reason is that `cohort` and `model` between them already describe what
this stage does. Every step name is distinct from `cohort.build`'s for §4.2's reason.

| # | kind | step | contents |
|---|---|---|---|
| 1 | `cohort` | `standardisation_population` | 126 → 107 → 104, per centre × arm × eligibility; case identifiers for every removal, `_record_removal`'s exact-naming discipline (§4.1) |
| 2 | `missingness` | `absence_by_standardisation_column` | `data.absence_by_column` over the [§14a] population |
| 3 | `model` | `standardisation_fit` | the design's columns, what was dropped, the fitted category set, `len(alpha)`, iterations, `converged_on`, the coefficient magnitudes against `POLR_MAX_ABS_BETA` (§5.2, §9.2) |
| 4 | `model` | `standardised_distributions` | the two distributions over `MRS_LEVELS`, the cumulative table, `RD_k`, and the two identities of §7.3 **as computed differences**, so the log shows they were checked |
| 5 | `model` | `treated_support` | the box per covariate, records outside by centre and arm, the never-IVT proportion [§14a] asks for (§10.1, §10.2) |
| 6 | `model` | `support_baseline` | the treated-vs-never-IVT-centre table over `BALANCE_SET`, three group counts in the caption, centre rows roled `grouping` (§10.3) |
| 7 | `model` | `standardisation_support_restricted` | sensitivity 1: `n_fit`, `n_average`, and the same output block (§11) |
| 8 | `model` | `hierarchical_fit` | sensitivity 2: `POLR_RI_NODES`, σ̂, `at_floor`, iterations, halvings, `converged_on`, and the per-centre posterior SDs (§12) |
| 9 | `model` | `standardised_hierarchical` | sensitivity 2's standardisation, and the statement that it conditions on `b_hat_c` (§12.6) |
| 10 | `model` | `standardisation_replicates` | the diagnostic grid: surviving draws, buckets, `len(alpha)` distribution, iteration distribution, `at_floor` rate, `n_average` spread, coefficient-magnitude spread (§13.4) |

**Entry 10's grid labels its denominators**, which is a `TODOS` item filed against Stage 10's
equivalent: it renders distributions over three different populations — replicates attempted,
replicates surviving, and replicates reaching `polr_ri` — and each block names its own. Stage 12 does
not fix Stage 10's grid; it declines to add a second unlabelled one.

The ledger: the Stage 2–4 path records **12** entries before this stage (`load` 7, `derive` 4,
`classify` 1 — measured, and the same 12 Stage 10 §12 counts), and Stage 12 adds **10**. In a Stage 14
driver running Stages 1–12, the total is Stage 11's **46** plus these 10 = **56**, the 12 shared
entries being recorded once. §20.9 asserts the count and the step-name disjointness together.

---

## 16. The entry points, written out

```python
def population(df: pd.DataFrame, audit: Audit) -> pd.DataFrame:
    """The [§14a] estimation population: eligible, covariate-complete, outcome-observed. 126 -> 104.

    Restriction 1 is NOT applied and `cohort.build` is not called: [§14a]'s population is the one
    restriction 1 removes a centre from [Stage 5 §9]. T6, T7.
    """


def all_centre(pop: pd.DataFrame, audit: Audit,
               over: pd.Series | None = None) -> Standardisation:
    """One [§14a] standardisation. `over` restricts the AVERAGING population, never the fit (§11).

    T1, T8, T10. Records entries 3 and 4 (or 7 when `over` is given).
    """


def support(pop: pd.DataFrame, audit: Audit) -> Support:
    """[§14a]'s support check: the treated-arm box, the never-IVT proportion, the baseline table.

    T9. Records entries 5 and 6. Calls `cohort.treating_centres` for its COMPLEMENT and
    `balance.smd` with unit weights; implements neither.
    """


def hierarchical(pop: pd.DataFrame, audit: Audit) -> Hierarchical:
    """[§14a]'s sensitivity 2: one common treatment effect over a random centre intercept.

    Records entries 8 and 9. The one caller of `model.polr_ri`.
    """


def inference(pop: pd.DataFrame, audit: Audit) -> bootstrap.Bootstrap:
    """The [§14] bootstrap: one `replicates` call, three arms per draw, seventy keys, no p-value.

    ALL FOUR STRATA (§13.1). Records entry 10. Calls `bootstrap.replicates`, `bootstrap.bucket`,
    `bootstrap.collect` and `bootstrap.intervals`; implements none of them.
    """
```

Five functions, and **`inference` recomputes the three arms inside its body rather than taking the
point-estimate records**, for Stage 10 §11's reason: a replicate must run the whole procedure, and a
body that took a fitted object would be resampling around a fixed fit.

---

## 17. Data flow into Stages 13 and 14

```
  population(df, audit)   -> DataFrame, 104 rows                       (§4.1)
  all_centre(pop, audit)  -> Standardisation   n_average 104           (§7)
  all_centre(pop, audit, over=sup.inside)
                          -> Standardisation   n_average  93           (§11)
  support(pop, audit)     -> Support           12 never-IVT, 3 outside (§10)
  hierarchical(pop, audit)-> Hierarchical      sigma, b_hat, at_floor  (§12)
  inference(pop, audit)   -> Bootstrap         70 keys, 2000 draws each, NO p
                                                                       (§13)
  the classified frame comes back unchanged                            (§0.2)

  every standardised probability, every risk difference, every limit and every
  centre intercept is DELIBERATELY ABSENT from this document           (§4.3)
```

### 17.1 Stage 13 [§14b]

Handed over in §19. In one line: **it reuses `population`, `all_centre`'s internals and `inference`'s
body shape, and it changes exactly one thing — the regime.**

### 17.2 Stage 14 [§16]

Six things Stage 12 owes the reporting layer that this stage cannot enforce, each because a
prespecified rule says so:

- **The different-populations statement, and it must be computed.** [§16]: *"Wherever a §14 estimate
  appears beside the primary, state that they target different populations and are not a like-for-like
  comparison, and label the §14 estimates as model-based transported results."* The statement is
  derivable from the two records — `Standardisation.n_average` is 104 and `Primary.in_estimate.sum()`
  is 92, at four centres against three — so Stage 14 prints the numbers, not an adjective. This is
  [§16] amendment DECISION 5's rule applied to a second statement.
- **`exp(beta)` is conditional wherever it appears**, and `Standardisation.measure` carries the string
  (§8). It must never be placed in a column headed with `Primary.odds_ratio`.
- **The support check appears beside the estimate, not in a supplement.** [§9]'s 2026-08-24 amendment
  established that an exceedance in a supplement lets a reader adopt the estimate without meeting it;
  the same argument applies with more force here, because the [§14a] transport assumption is *"much
  stronger than observing that unadjusted mRS distributions look alike across centres"* and it is
  **not testable from these data**. The numbers that must travel with the estimate: 12 never-IVT
  patients, 3 outside the treated box, and the largest non-grouping |SMD| in the baseline table
  (§10.3), which is above 0.8 against a threshold of 0.10.
- **The cutpoint-collapse rate is printed beside the `RD_5` interval.** 4.15% of draws lost mRS 5, in
  which `RD_5` equals `RD_4` exactly (§6.3). A reader comparing the two intervals is entitled to know
  that 83 draws of 2000 made them identical by construction.
- **The hierarchical arm prints its `at_floor` rate and USZ's posterior SD.** 22.0% of replicates
  collapse to the pooled model, and the interval on `hier.sigma` has a limit that means "the boundary"
  in those draws (§12.5); USZ's intercept is the least precise of the four and the arm rests on it
  (§12.6).
- **No [§14a] quantity carries a p-value, and none is in a [§13] family.** §3.3 is why, and Stage 11's
  Benjamini–Hochberg must not acquire a fifth family by counting Stage 12's keys — which is the
  mistake Stage 10 §12.1 warned about in its own direction ("seven p-values, not fourteen").

---

## 18. What Stage 12 amends in Stages 1–11

The full ledger, so that no amendment is discovered during implementation. **Three shipped modules and
four test modules; the two new estimators are both in `model.py` and nothing in `outcome.py`,
`propensity.py`, `balance.py`, `cohort.py`, `eligibility.py`, `derive.py` or `data.py` changes at
all.**

| File | Amendment | Why |
|---|---|---|
| `model.py` | `ordinal_probabilities` (§6) and `polr_ri` + `RIFit` (§12). Public surface 6 → 8 | §0.1: an estimator-shaped function lives in `model.py`. `polr`, `firth`, `design`, `predict`, `complete_cases` and `Fit`/`PolrFit` are untouched, and §20.11 asserts it by number |
| `bootstrap.py` | `intervals(draws, tested)` becomes the ninth public name; `run` and Stage 11's two callers move onto it | §13.3. The `TODOS` trigger is a third caller and this is it. `run`'s twenty-six intervals are asserted byte-identical across the change (§20.11) |
| `config.py` | `POLR_RI_NODES`, `POLR_RI_SIGMA_FLOOR`, `POLR_RI_MAX_ITER` in the model-tolerance block; `SUPPORT_COVARIATES` as a computed view in the covariates block; three `FAILURE_BUCKETS` tokens | §12.3, §12.5, §12.4, §10.1, §14. The fence below |
| `tests/test_model.py` | `ordinal_probabilities` against a hand-built fit and against `_ord_pieces`; `polr_ri` against its analytic gradient, its node-count stability and its floor | §20.5, §20.12 |
| `tests/test_bootstrap.py` | `intervals`' three p-rules; the surface count 8 → 9; the raise-site scan scope grows to `standardise.py` | §13.3, §14 |
| `tests/test_config.py` | **five** assertions: `POLR_RI_NODES` odd and ≥ 9; `0 < POLR_RI_SIGMA_FLOOR < 1e-2`; `STANDARDISATION_COVARIATES == tuple(c for c in PS_COVARIATES if c != "center")` — [§15]'s permission as a static check (§5.1); `SUPPORT_COVARIATES` is `STANDARDISATION_COVARIATES` minus `CATEGORICAL`; the three new `FAILURE_BUCKETS` tokens | §20.2 |
| `tests/test_reference_r.py` | the `clmm` oracle, under the existing `reference_r_ordinal` gate and `_r_environment` | §12.7 |
| `implementation_roadmap.md` | Stage 12 gains its `**Spec:**` line; five corrections and three additions | §27 |
| `TODOS.md` | two close, one has its stated reason corrected, six are new | §21 |

```python
# config.py — the additions, as they are to be pasted. Stage 12 §18

# --- [§14a] the random-intercept fit's definition ------------------------------------------------
#
# POLR_RI_NODES is NOT a tolerance and it CHANGES AN ANSWER, so it sits under PERCENTILE_METHOD's
# argument and not under POLR_TOL's: it is part of the definition of the objective, and [§10] refits
# this model in every one of N_BOOT replicates. Measured (Stage 12 §12.3), refitting at each count:
# NON-adaptive quadrature returns sigma_hat 0.720 / 0.803 / 0.832 at 5 / 7 / 9 nodes against
# 0.542426850 at 31 -- 33% to 53% too large -- while ADAPTIVE quadrature is stable to eight
# significant figures from NINE. That gap is why Stage 12 §12.2 specifies adaptive and why the
# constant lives here: a node count is a modelling decision at low counts and an arithmetic detail
# only above them.
#
# 11 is one step above the measured plateau. Odd, because a symmetric rule with an odd node count
# places a node at the conditional mode, which is where the mass is; test_config.py asserts both.
# The measured cost of the two nodes above the plateau is 0.5 s per 2000-replicate arm.
POLR_RI_NODES: Final[int] = 11

# The sigma -> 0 boundary is REACHED -- 22.0% of replicates, measured (Stage 12 §12.5) -- and under a
# log-sigma parametrisation it is at -inf, where the conditional-mode Newton's 1/sigma^2 overflows.
# This is a floor on the SD, so 1e-8 on the variance: five orders of magnitude below the smallest
# non-boundary sigma_hat in 200 replicates, so it separates "collapsed to the pooled model" from
# "small" without being reachable by a genuinely small non-zero variance. A fit at the floor is NOT a
# failure -- [§14a] names sigma^2_C = 0 as a legitimate answer -- and RIFit.at_floor records it.
POLR_RI_SIGMA_FLOOR: Final[float] = 1e-4

# Its own cap and not POLR_MAX_ITER's, because the two estimators are not on the same scale: `polr`
# converges in 4 to 6 Newton iterations and `polr_ri` in 23 to 47 BFGS ones (measured, Stage 12
# §12.5), so one shared cap would be either loose for the first or tight for the second.
POLR_RI_MAX_ITER: Final[int] = 200


# --- [§14a] the support check's covariates -------------------------------------------------------
#
# A computed view, never a fourth literal list, for OUTCOME_COVARIATES' reason and invariant 3's.
# [§14a] says "each continuous covariate"; this resolves it to "not a declared factor", which on v7
# is the WIDER of the two readings and is MEASURED to be inert -- `sex`, `prestroke_mrs` and
# `atrial_fib` exclude nobody, because the treated arm covers both levels of each binary and reaches
# prestroke_mrs = 3 (Stage 12 §10.1). The wider reading is chosen because it cannot be wrong on a
# workbook where the narrow one is right.
SUPPORT_COVARIATES: Final[tuple[str, ...]] = tuple(
    c for c in STANDARDISATION_COVARIATES if c not in CATEGORICAL)          # [§14a]
```

Three `FAILURE_BUCKETS` tokens are added beside the existing sixteen: `"T1"` →
`"degenerate_design"`, `"polr_ri:"` → `"nonconvergence"`, `"T4"` → `"degenerate_design"`. `"T2"` and
`"T3"` both lead with `polr_ri:` and share a bucket, which is `"polr:"`'s arrangement and Stage 10
§7.2's argument for it, unchanged.

---

## 19. Handover to Stage 13 [§14b]

```
  pop  = standardise.population(classified, audit)                    104 rows, all 4 centres
  std  = standardise.all_centre(pop, audit)                           [§14a]
  sup  = standardise.support(pop, audit)
  hier = standardise.hierarchical(pop, audit)
  boot = standardise.inference(pop, audit)                            70 keys, 2000 draws

    126 -> 107 -> 104   the ledger, and the 2 outcome-missing are BOTH at USZ     [§4.1]
    12 / 3              never-IVT patients, and how many are outside the box      [§10.2]
    4.15%               replicates losing mRS 5, in which RD_5 == RD_4 exactly    [§6.3]
    0 / 2000            failures. Nothing in this stage drops a replicate         [§13.4]
    22.0% / 24.0%       polr_ri replicates at the sigma floor / above sigma = 1   [§12.5]
    ~18.4 min           the bootstrap, of which the hierarchical arm is 95%       [§2]
    10                  audit entries, taking a Stages 1-12 ledger to 56          [§15]
```

**Four things Stage 13 must know, and three of them are prohibitions.**

1. **[§14b]'s population is NOT this one, and `population()` must not be reused unchanged.** [§14a] is
   restricted to eligible patients; [§14b] is *"the only §14 analysis whose population includes
   contraindicated patients"*. Stage 13 takes the same unrestricted classified frame and applies
   **neither** restriction — so it needs `population`'s covariate-completeness and outcome-observed
   half without its eligibility half. The clean seam is a keyword on `population`, and Stage 12 does
   not add one speculatively (§22).
2. **The regime is the one thing that changes, and it is a function of `eligibility`, not a
   constant.** [§14b]'s active regime assigns bridging to eligible patients and direct EVT to
   contraindicated ones. In §7.1's terms that is a *vector* written into the treatment column instead
   of a scalar — the same overwrite, the same fitted design, the same `ordinal_probabilities` call.
   Nothing in §6 or §7.1 names the value being written, deliberately.
3. **DECISION 1a already closed [§14b]'s open question and Stage 13 must not look for a third value.**
   The `eligibility` column is **binary** — `ELIGIBILITY_ORDER` has two entries — so "assign treatment
   as a function of eligibility" is a binary decision over a binary column. Stage 5 §9 records this;
   it is repeated because the [§14b] text predates the amendment.
4. **The Accept-when is a prohibition and it is checkable.** *"No contraindicated patient is ever
   assigned a predicted IVT outcome."* Under §7.1's design-overwrite that is an assertion on the
   treatment column of the counterfactual design — `Xa[TREATMENT]` must be 0 on every row where
   `eligibility` is `INELIGIBLE` — and it is one line, checked before the prediction rather than after
   the average.

**And one thing Stage 13 inherits that Stage 12 measured for it**: the fit will be over a larger and
more heterogeneous population, so the cutpoint-collapse rate, the drop rate and `polr_ri`'s boundary
rate are all Stage 13's to re-measure. Nothing here transfers except the shapes.

---

## 20. Acceptance criteria

### 20.0 `tests/fixtures_stage12.py` — this row is the one place its contents are enumerated

Five functions and four constants, following Stage 10 §15.0's rule that a constructed fixture is
given as code and not described, because Stage 9 §22.3 found three fixtures pinned to ten decimals and
described rather than given, which made its acceptance criteria unperformable.

- `four_centre_frame()` — an eligible population over four declared centres, one of which has no
  treated patient, all seven mRS levels occupied.
- `collapsed_level_frame()` — the same with mRS 5 absent, so `polr` fits five cutpoints and §6.3's
  re-expansion is exercised on a frame rather than on a mock.
- `no_exposure_frame()` — every patient treated, so `design` drops the treatment column and T1 fires.
- `two_group_ri_frame()` — a frame with a known non-zero between-group intercept, for `polr_ri`.
- `flat_ri_frame()` — a frame generated with `sigma = 0`, so `polr_ri` reaches the floor and
  `at_floor` is exercised.
- `RI_SEED`, `RI_SIGMA_TRUE`, `RI_BETA_TRUE`, `COLLAPSE_LEVEL`.

### 20.1 The population is [§14a]'s and not the cohort's

- `population(classified, audit)` returns **104** rows over **four** centres, and `cohort.build` on the
  same frame returns **93** over **three**. Both asserted in one test, so the difference is the
  assertion.
- Every USZ record in the result has `ivt == 0` — no treated patient exists at a never-IVT centre by
  definition, and a workbook in which one does must fail T9 before this line is reached.
- The three [§11] losses are asserted **by cause**: one covariate-incomplete on `core_ml` and
  `tmax6_ml`, two outcome-missing, **both at USZ** (§4.1).
- `derive.derive_cohort` is not called on this path, asserted by scanning `standardise.py` for the
  name — the Stage 1 §7 scan's pattern — because the failure it prevents is a subgroup median computed
  on the wrong population and is silent.

### 20.2 The config additions

Five assertions, listed in §18. The load-bearing one is
`STANDARDISATION_COVARIATES == tuple(c for c in PS_COVARIATES if c != "center")`, which makes [§15]'s
*"permitted inside §14 only"* a static check rather than a sentence: a covariate added to the [§14a]
model without being a [§6] confounder fails at import.

### 20.3 The estimand keys are derived and not counted

`len(keys) == 3 * (1 + len(C.MRS_THRESHOLDS) + 2 + 2 * len(C.MRS_LEVELS)) + 1`, asserted against the
config and never against 70. And `all(interval.p is None for interval in boot.intervals.values())` —
§3.3's no-p rule, asserted over every key rather than spot-checked.

### 20.4 The fit is unweighted, and the call site is the assertion

`model.polr` is called with two positional arguments and no `w`, asserted by inspecting the call
rather than by comparing numbers — `w=None` and `w=ones` are numerically identical and semantically
different, and `polr`'s docstring says the acceptance suite must turn on the difference being visible
in the call.

### 20.5 The distribution is a distribution, structurally

- On `four_centre_frame()` and on the workbook: every row of `ordinal_probabilities` sums to 1.0 to
  1e-12, and every value is in [0, 1].
- **Constructed, not measured**: a hand-built `PolrFit` with known `alpha` and `beta` is compared
  against `expit` differences computed by hand, so the assertion is about the arithmetic and not about
  the cohort.
- The structural bounds are **exact**: the first column equals `expit(alpha[0] + eta)` exactly and the
  last equals `1 - expit(alpha[-1] + eta)` exactly, which fails if a future edit replaces the literal
  0.0/1.0 columns with an `expit` of a clipped large number (§6.1).
- On `collapsed_level_frame()`: `ordinal_probabilities` returns **six** columns, the re-expansion gives
  **seven** with a structural zero at `COLLAPSE_LEVEL`, the sum is still 1.0, the cumulative sequence
  is still monotone, and `RD_4 == RD_5` exactly (§6.3).
- T10 fires on a deliberately corrupted distribution, so the check is known to be live.

### 20.6 The two identities

`abs(std.mrs_0_2 - std.rd[2]) == 0.0` and `abs(std.mortality + std.rd[5]) < 1e-15`, on the workbook and
on both fixtures. And a test asserts `mortality` is computed from the **distributions** and not as
`-rd[5]`, by monkeypatching `rd` and requiring `mortality` to be unaffected — because an
implementation that derived one from the other would pass the identity trivially.

### 20.7 The mortality interval is NOT the reflection of `RD_5`'s

The measurement of §7.4, as a test: over `boundary_draws()`-style constructed draws under
`C.PERCENTILE_METHOD`, `percentile_ci(-x)` and the reflected `percentile_ci(x)` **differ**, and the
test asserts the difference is non-zero rather than asserting a magnitude. A companion asserts they
**agree** under `method="linear"`, so the test documents that the disagreement is the pin's and not an
arithmetic error — which is the shape Stage 10 §15.8 used for the same pin.

### 20.8 The conditional-measure guard

`Standardisation.measure is _CONDITIONAL` for every record the module constructs; `hasattr(std,
"odds_ratio")` is **False**; and a scan asserts the string `"odds_ratio"` does not appear in
`standardise.py` outside a comment. §8's three enforcements, each as an assertion.

### 20.9 The audit entries and their step names

Ten entries in the declared order; the step-name set is **disjoint** from `cohort.build`'s, asserted as
a set intersection rather than against a list (§4.2); entry 1 names exactly as many case identifiers as
its `n` claims, `cohort._record_removal`'s discipline; and the Stages 1–12 ledger totals 56.

### 20.10 The support check

- The box's two readings agree on this workbook — the set of records outside `SUPPORT_COVARIATES`'s box
  equals the set outside the box over `("age", "nihss_baseline", "core_ml", "tmax6_ml")` — asserted, so
  §10.1's "inert" is a check and not a claim.
- 93 inside, 11 outside, all 11 in the comparator arm; 12 never-IVT patients of whom 3 are outside.
- T9 fires on a frame in which the never-IVT set is not `EXPECTED_NEVER_IVT`.
- `center = USZ`'s baseline SMD is `nan`, and the test names Stage 7 §5.3's branch 4 as the reason, so
  a future change to `smd`'s nan routes surfaces here.
- The four centre rows carry role `"grouping"` and the other fifteen do not.

### 20.11 Stage 12 changes nothing upstream, asserted by number

- The Stage 1–11 pipeline is run with and without this stage; every landed number and the first 46
  audit entries are byte-identical. Stage 11 §15.12's test, extended by one stage.
- **`bootstrap.run`'s twenty-six intervals are byte-identical across the `intervals` extraction**
  (§13.3), which is the specific guard Stage 11 §16 said the refactor would need before it could be
  done.
- `model.polr`, `model.firth`, `model.design`, `model.predict` and `model.complete_cases` are
  untouched, asserted by their behaviour on the existing Stage 6–10 fixtures rather than by diff.

### 20.12 `polr_ri`

- **The analytic gradient against central differences**, at σ ∈ {0.2, 0.54, 1.5, 3.0}, relative error
  below 1e-6 at every one — §12.4's frozen-node approximation pinned at the boundary where it would
  start to matter.
- **Node-count stability**: σ̂ at 9, 11, 15 and 31 adaptive nodes agrees to 1e-7; and σ̂ at 5, 7 and 9
  **non-adaptive** nodes does **not**, which is the measurement §12.2 rests on, asserted so that a
  change to non-adaptive quadrature fails rather than shifts an answer.
- **The floor**: on `flat_ri_frame()`, `at_floor` is True and `sigma == POLR_RI_SIGMA_FLOOR`, and the
  fit does **not** raise.
- **Recovery**: on `two_group_ri_frame()` with a known `RI_SIGMA_TRUE`, σ̂ is within a stated tolerance —
  the only test here that checks the estimator estimates the right thing rather than that it runs.
- **The collapse**: with `RI_SIGMA_TRUE = 0`, `polr_ri`'s `alpha` and `beta` agree with `model.polr`'s
  on the same frame to 1e-6, which is [§14a]'s own *"it collapses to the pooled model"* as an
  assertion.
- **`ordinal::clmm`**, gated (§12.7): coefficients negated, thresholds not, `sigma` not, `nAGQ`
  passed.
- T2, T3 and T4 each fire on a constructed frame, and each token is a `C.FAILURE_BUCKETS` key.

### 20.13 The replicate loop

- The drawn-frame sequence is byte-identical to the one a body returning `None` produces, at the same
  seed, frame and stratum — Stage 10 §15.2's property, now over a three-arm body (§13.1).
- **All four strata appear in every drawn frame** and each stratum's size is exact, USZ's 12 included.
- The three failure groups each cost exactly their own keys, driven by three constructed failures, and
  `Draws.n_attempted` reconciles for every one of the seventy keys.
- A `POLR_MAX_ABS_BETA` breach costs `beta` and `support.beta` and **nothing else** — §9.1's key
  granularity, which is the one place this stage departs from Stage 8 and therefore the one that most
  needs a test.

---

## 21. Known gaps carried forward

1. **The [§14a] transport assumption is not testable and nothing here tests it.** [§14a] says so and
   adds that differences between centre-specific estimates are *not* a test of it. Stage 12 reports
   the support check and the baseline table; it does not offer a diagnostic that could be mistaken for
   evidence. **Status: by design, not a gap to close.**
2. **`POLR_MAX_ABS_BETA`'s calibration — CLOSED, negatively (§9.2).** The `TODOS` item's trigger was a
   Stage 12 design; measured, the largest coefficient in 2000 fits of the ten-column design is 3.21
   against a bound of 14.0 and below Stage 8's measured band's *lower* edge. **Trigger, rewritten:** a
   workbook or a Stage 13 population on which any Stage 12 coefficient exceeds 8.79 — Stage 8 §6.3's
   largest measured non-degenerate `|beta|` — at which point the band's calibration on a one-column
   design becomes load-bearing here.
3. **The absolute convergence tolerance — CLOSED, negatively (§13.4).** `Σw` is 104 here against
   roughly 23 for the [§7] weighted cohort, which is the change of scale the item names, and every one
   of 2000 fits converged in 4 to 6 iterations on the likelihood criterion. **Status: open, narrowed** —
   it remains a real property of an absolute tolerance, and this is evidence rather than proof.
4. **`Replicate` cannot describe a multi-fit body — the trigger has FIRED and the fix is still
   declined, and the item's stated reason was wrong.** It reads *"[§14a]'s standardisation is the
   candidate, since it fits per arm."* Standardisation does **not** fit per arm — it fits once and
   predicts twice, which is §7.1's whole point — so the stated mechanism never occurs. What does occur
   is that the hierarchical arm adds a **second estimator** to the body. The fix is declined because
   widening `n_alpha` and `polr_iterations` to dicts touches every Stage 10 assertion over them to gain
   what `hier.sigma` already provides as a first-class key. **Trigger, rewritten:** a stage wanting the
   *cutpoint distribution* of a second fit, which `hier.sigma` does not carry.
5. **Parallelising the replicate loop — the trigger has fired.** `TODOS` names Stage 12's bootstrap as
   the sole remaining trigger and `SeedSequence.spawn(N_BOOT)` as the right design. Measured: 18.4
   minutes, of which the hierarchical arm is 95%. **Status: open, and it is now the largest single
   cost in the pipeline.**
6. **NEW — the marginal-versus-conditional standardisation under the random-intercept model.** §12.6
   implements [§14a]'s reading; the marginal one is a different estimand and would need a [§14]
   amendment. **Trigger:** a reviewer asking what the arm says for a patient at an unobserved centre.
7. **NEW — a percentile interval on a parameter at a boundary.** 22.0% of `hier.sigma`'s draws are at
   `POLR_RI_SIGMA_FLOOR`, so its lower limit is the floor and means "the boundary". §17.2 requires the
   rate printed beside it; a bounded parameter's interval is a known-hard problem and this pipeline
   does not solve it. **Trigger:** any reading of that interval as a range for the between-centre SD.
8. **NEW — `RD_5` and `RD_4` coincide in 4.15% of draws** (§6.3), which narrows `RD_5`'s sampling
   distribution by an amount nobody has quantified. **Trigger:** a manuscript reporting `RD_5` and
   `RD_4` as separate findings.
9. **NEW — the two outcome-missing records are both at USZ** (§4.1). The centre whose inclusion is
   [§14a]'s reason for existing loses 2 of 14 to outcome ascertainment, so the transported population
   is 12 and not 14. **Trigger:** a workbook in which the never-IVT centre's outcome completeness is
   materially worse.
10. **NEW — the frozen-node gradient's neglected term** (§12.4). Measured at 1e-8 relative, six orders
    below `POLR_SCORE_TOL`. **Trigger:** a change to the conditional-mode search, or a population where
    a centre's posterior is far from Gaussian.
11. **NEW — Stage 12's spec is written against an unlanded Stage 11.** The three dependencies are named
    in the header. **Trigger:** Stage 11 landing with `bucket`/`collect` private, or `intervals` not
    extracted.

---

## 22. What Stage 12 deliberately does not decide

- **Whether a refit-on-support arm should exist** beside §11's averaging-only restriction. It is a
  [§14] amendment and the PI's; §11 records the argument on both sides.
- **Whether a population-conditional bootstrap variant should exist**, holding the never-IVT stratum
  fixed (§13.1). It is a different inferential statement, not a different implementation.
- **Whether the hierarchical arm's standardisation should be marginal** rather than conditional on
  `b_hat_c` (§12.6, §21 item 6).
- **What [§14b]'s population seam looks like.** Stage 13 needs `population` without its eligibility
  half (§19 item 1); Stage 12 does not add a keyword speculatively, for Stage 5 §15's reason — a
  parameter shaped like a future requirement is a parameter that will be used for something else
  first.
- **Whether the pipeline is parallelised** (§21 item 5). Stage 12 measures the cost that makes it a
  question.
- **How [§14a]'s estimate is placed beside the primary's in the manuscript.** Stage 14 [§16], and
  §17.2 lists what must travel with it.

---

## 23. NOT in scope for Stage 12

| Considered | Why declined |
|---|---|
| A `centres=` or `restrict=` keyword on `cohort.build` | Stage 5 §15 declined it once from the cohort side; §0.2 declines it again from this side. `build` applies [§3]'s two restrictions and [§14a] is not a variant of them |
| A shared [§11] mask helper across `outcome` and `standardise` | §4.1. The two conjuncts differ in content and agree only in shape; one function with two callers agreeing on nothing but the `&` |
| A `gcomp.py` / `standardise.py` split | §0.1. The general half is `ordinal_probabilities` and it goes in `model.py`; the rest is nine lines |
| `scipy.optimize` for `polr_ri` | §2, §12.4. `scipy` is test-only by policy and the hand-written BFGS is measured at 0.76 s; a runtime dependency for a 40-line loop doubles the reproducibility surface for nothing |
| `statsmodels.OrderedModel` for the pooled fit | Stage 8 §18b measured that it has no weight support and uses the other sign convention. It is an oracle at Stage 8 and not an estimator, and that does not change here |
| A treatment × centre interaction, or a random slope | [§14a] and [§15] both forbid them at four centres, in as many words |
| A test of the proportional-odds assumption in the [§14a] model | [§15] declines one for [§8]'s model and the argument transfers: it would be a second estimator on the primary outcome. [§16]'s constant-shift statement covers it |
| A p-value on any [§14a] quantity | §3.3. [§14] prescribes percentile intervals and states no null |
| Adding [§14a]'s keys to a [§13] family | §17.2. [§13] scopes Benjamini–Hochberg to the secondary and safety outcome families |

---

## 24. What already exists, and what to lift

| Need | Exists | Lift |
|---|---|---|
| The [§14a] covariate list | `config.py:562` | Call it. Do not write it out |
| The population predicate | `eligibility.retained` (`eligibility.py:290`) | Call it. Stage 5's restriction 2 and this stage's are the same predicate |
| Covariate completeness | `model.complete_cases` (`model.py:254`) | Call it with `STANDARDISATION_COVARIATES` |
| The design | `model.design` (`model.py:204`) | Call it **once**, on the observed frame (§7.1) |
| The ordinal fit | `model.polr` (`model.py:844`) | Call it with no `w` |
| The counterfactual pattern | `outcome._counterfactuals` (`outcome.py:1087`) | Copy the *shape* — overwrite the column of the fitted design — not the function; it is Firth-specific |
| The never-IVT set | `cohort.treating_centres` (`cohort.py:96`) | Take its **complement**. `cohort.py:112-119` writes the expression |
| The standardised mean difference | `balance.smd` (`balance.py:252`) | Call it with unit weights. Stage 7 §5.5 made it public for this |
| The balance row list and labels | `balance._levels`, `C.BALANCE_SET` | `_levels`' rule — range over declared levels, never over observed — is what §10.3 needs; the function is private, so §10.3 builds its rows the same way rather than importing it |
| The resampler | `bootstrap.resample` (`bootstrap.py:290`) | Call it via `replicates`. General in stratum by design |
| The replicate driver | `bootstrap.replicates` (`bootstrap.py:347`) | Call it. One-argument body |
| Failure classification and reconciliation | `bootstrap.bucket`, `bootstrap.collect` (public at Stage 11 §5.4) | Call them. Do not re-implement |
| The percentile interval | `bootstrap.percentile_ci` (`bootstrap.py:378`) | Via `intervals` (§13.3) |
| The audit ledger | `data.Audit`, `data.absence_by_column` | Call them. Distinct step names (§4.2) |
| Gauss–Hermite nodes | `numpy.polynomial.hermite.hermgauss` | Call it. No `scipy` |
| The R ordinal gate | `test_reference_r.py`'s `reference_r_ordinal`, `_r_environment` | Reuse both. `clmm` is in the package already gated |

**Nothing worth importing exists for the two new functions.** There is no maintained Python
implementation of a random-intercept proportional-odds model — `statsmodels` has none at any API, and
Stage 8 §18b already measured that its ordinal fitter has no weight support and the other sign
convention — which is why §12 is a specification and not a call site.

---

## 25. Implementation tasks

In order, each landing green before the next begins.

1. **`config.py`** — the four constants and the three `FAILURE_BUCKETS` tokens (§18 fence).
   `test_config.py`'s five assertions.
2. **`model.ordinal_probabilities`** (§6) with its hand-built and `_ord_pieces` tests (§20.5). It is
   independent of everything else here and it is what the rest rests on.
3. **`bootstrap.intervals`** (§13.3), with Stage 11's two callers and `run` moved onto it, and
   §20.11's byte-identity guard on `run`'s twenty-six intervals. Done before Stage 12 has any keys, so
   the refactor is verified against the landed stage rather than against new code.
4. **`standardise.population`** (§4) and its audit entries. §20.1, §20.9.
5. **`standardise.all_centre`** (§5, §7, §8) — the design-overwrite, the re-expansion, the two
   identities, the conditional-measure guard. §20.4–20.6, §20.8.
6. **`standardise.support`** (§10). §20.10.
7. **`model.polr_ri`** (§12) — the largest single piece, and in this order: the adaptive quadrature
   objective; the analytic gradient with §20.12's central-difference test **before** the optimiser
   exists; the BFGS loop; the floor; `RIFit`. Then the `clmm` oracle.
8. **`standardise.hierarchical`** (§12.6). §20.12.
9. **`standardise.inference`** (§13) — the three-arm body, three failure groups, seventy keys.
   §20.3, §20.13.
10. **The roadmap and `TODOS.md` amendments** (§21, §27), landing with the code.

**Step 3 before step 5 is the one ordering that is not obvious and is deliberate**: extracting the
interval loop while Stage 12 has no callers means its byte-identity guard runs against Stage 10's
landed twenty-six intervals and nothing else, so a discrepancy is unambiguously the refactor's.

---

## 26. Verification record

Every number in this document, and where it came from. Twenty-two probes, run 2026-08-27 against the
landed Stages 1–10 on the v7 workbook (`DATA_SHA256` unchanged), Python 3.12.12 / numpy 1.26.4 /
pandas 2.3.3.

| Probe | What it measured | Result | Section |
|---|---|---|---|
| A1 | the population ledger | 126 → 107 → 104; 1 covariate-incomplete on `core_ml`+`tmax6_ml`; 2 outcome-missing, both USZ | §4.1 |
| A1b | `treating_centres`' complement against `EXPECTED_NEVER_IVT` | `('USZ',)`, agree | §10.2 |
| A2 | the mRS level set | all 7 occupied; mRS 5 carries 3 (HUG 2, Lugano 1) | §6.2 |
| A3 | the design | 10 columns, 0 dropped, rank 10, 16 parameters on 104 | §5.2 |
| A4 | the pooled fit | 5 iterations, likelihood, 0 rescales, 0 halvings, max abs coefficient 0.73 | §5.2 |
| A5 | the treated support | box per covariate; 93 inside / 11 outside, all comparator; 3 of 12 USZ outside | §10.1, §10.2 |
| A5b | the two "continuous" readings | identical outside-sets; `sex`, `prestroke_mrs`, `atrial_fib` exclude nobody | §10.1 |
| A6 | the baseline table | 39 / 12 / 53; `center = USZ` is `nan`; largest non-grouping SMD above 0.8 | §10.3 |
| A7 | the audit step names | 12 entries before this stage; `cohort.build`'s six names enumerated | §4.2, §15 |
| B0 | the structural Accept-when, point estimate | sums to 1 exactly, monotone, all in [0, 1] | §20.5 |
| B0b | the two identities | `|mrs_0_2 − RD_2|` 0.0; `|mortality + RD_5|` 8.3e-17 | §7.3 |
| B1 | 2000 replicates, three-arm body | 2000/2000 survived, no failures, 43.9 and 50.1 s for two arms over two runs | §13.4, §2 |
| B2 | the cutpoint collapse | 83 / 2000 = **4.15%**, always mRS 5 | §6.2 |
| B3 | convergence across replicates | iterations 4–6, all likelihood, `first_step_norm` 1.01–22.8 | §13.4 |
| B4 | separation on ten columns | max `|beta|` 2.00, max `|gamma|` 3.21, 0 reaching 14.0; every reported quantity in range | §9.1, §9.2 |
| B5 | dropped design columns | 0 in 100.00% of replicates | §5.3 |
| B6 | the support arm's averaging population | 69–102, median 91, point 93 | §11 |
| B7 | the reflection check | `inverted_cdf` 9.8e-5 / **1.7e-3**; `linear` 1.1e-16 / 5.6e-16 | §7.4 |
| B9 | the Accept-when across replicates | worst `|sum − 1|` 6.7e-16; monotone in every one | §20.5 |
| C1–C4 | the first `polr_ri` prototype | non-adaptive, numerical Hessian: 39.0 s per fit, **21.7 h** at `N_BOOT` | §12.4, §12.8 |
| D1–D4 | analytic gradient, non-adaptive | gradient to 7.8e-9; 0.82 s per fit; σ̂ 0.72/0.80/0.83 at 5/7/9 nodes | §12.3, §12.4 |
| E1–E4 | adaptive against non-adaptive | non-adaptive error 4.8e-1 at 31 nodes at σ=3, non-monotone; adaptive 8.9e-4 at **3** | §12.2 |
| F1–F4 | the specified shape | frozen-node gradient 1e-8 relative; 0.76 s per fit; 200/200 converge; **22.0%** at the floor, **24.0%** above σ=1; **17.5 min** at `N_BOOT` | §12.4, §12.5, §12.8 |

**Six probes changed what this document says**, each named in the Status block: B0b (§7.3), B7 (§7.4),
B4 (§9), E1 (§12.2), C4/D2/F2 (§12.4), F4 (§12.5).

**Two probes closed `TODOS` items negatively** — B4 (§9.2) and B3 (§13.4) — and one falsified the
stated reason of a third while confirming its trigger (§21 item 4).

**The probe suite was re-run end to end after this document was drafted and every number above
reproduced exactly**, the seed being `C.SEED` and the workbook unchanged — the sole exception being
wall-clock timings, which are machine state and are quoted as ranges or as both observations. The
Stages 1–10 suite was run at the same time: **1588 passed, 2 skipped**, the two skips being the R
oracles that this machine's segfaulting R cannot open (§12.7). This document changes no shipped code
and the run is recorded as the guard against an accidental edit rather than as evidence about Stage 12.

**What was NOT measured, and is therefore not claimed.** The `ordinal::clmm` oracle has not been run:
R segfaults on this machine at startup for the `uname`/PATH reason `test_reference_r.py`'s banner
documents, and opening that gate is an implementation task (§25 step 7) and not a spec task. So
`polr_ri` is verified here against its own analytic gradient, against node-count stability, against the
σ = 0 collapse to `model.polr`, and against nothing external. §12.7 specifies the oracle; §21 item 11
is not it. **This is the one estimator in this pipeline whose specification ships without an
independent implementation having agreed with it**, and that is stated plainly rather than left to be
discovered at step 7.

---

## 27. What this spec changed elsewhere

`implementation_roadmap.md`, Stage 12 — five corrections and three additions, landing with this
document:

1. **Its `**Spec:**` line**, as Stages 8–11 carry.
2. **Correction.** *"Outputs: … cumulative `RD_k`; standardised mRS 0–2 risk difference; standardised
   mortality difference"* reads as four computations; the last two are `RD_2` and `−RD_5` (§7.3). The
   roadmap keeps all four names and gains the identity.
3. **Correction.** *"Inference: patient-level bootstrap stratified by centre"* did not say whether the
   never-IVT stratum is resampled, which Stage 10 §12.2 had left open. It is (§13.1).
4. **Correction.** *"Sensitivity: random centre intercept"* is one bullet for what is a **new
   estimator** with its own quadrature rule, its own optimiser, its own boundary behaviour and 95% of
   the stage's compute (§12). The roadmap gains three sentences and the node-count constant.
5. **Correction.** *"Support check: treated-arm distribution of each continuous covariate"* leaves
   "continuous" undefined; it resolves to "not a declared factor" and the two readings are measured
   equivalent here (§10.1).
6. **Correction.** The Accept-when's *"the standardised category probabilities sum to 1 … the
   cumulative probabilities are monotone"* are true **by construction** of the differencing (§6.1), so
   they are self-checks against a coding error rather than properties of a fit. Worth stating so that
   a reader does not read them as diagnostics that could inform about the model.
7. **Addition.** The cutpoint collapse: 4.15% of replicates lose mRS 5, in which `RD_5` equals `RD_4`
   exactly (§6.3). Nothing in the roadmap anticipated it and it shapes the module.
8. **Addition.** Stage 12 fits **no** propensity model **and reads no Stage 6–11 result** — the
   roadmap says the first and not the second, and the second is what makes the stage runnable
   independently.

`TODOS.md` — two items close (§21 items 2, 3 — both **negatively**, which is recorded rather than
allowing a closed item to read as a solved problem), one has its stated reason corrected while its
trigger fires (item 4), one has its trigger fire and stays open as the largest cost in the pipeline
(item 5), and six are new (items 6–11).
