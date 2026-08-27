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

**It then went through an engineering review on 2026-08-27, which changed fourteen things** — seven
that needed a decision and six smaller corrections, plus one note about the environment. Five of the
seven were **internal contradictions in the draft**: two sections prescribing incompatible things, a
guard with no identifier, a record field holding a type it was not declared to hold. One was a
**reversal**, decided by reading [§14a]'s Report sentence rather than by argument, and it moved the
estimand key count from seventy to sixty-nine.

**§27.1 is the ledger of what the review changed; §26.1 is what it did NOT check.** Read §26.1 before
trusting this document's statistical choices: the review checked internal consistency and this
document's claims about landed code, it did not check the statistical premises, it re-ran nothing,
and **no independent cross-model review was run**. §21 item 14 is the residual and is the largest
open risk here.

**IT WAS THEN IMPLEMENTED ON 2026-08-27, AND §28 IS THE RECORD.** Every figure in §26 reproduced on an
independent implementation; §25's step 0 ran and replaced the hierarchical arm's projections with
`N_BOOT` measurements; **§26's largest verification gap closed** — R works on the implementing machine,
so the `ordinal::clmm` oracle ran and agrees with `polr_ri` to about 1e-5 on the coefficients, the
thresholds and σ alike (§28.5); and **three things this document did not anticipate were found by
running the code rather than by reading it**: two numerical defects in `polr_ri`'s quadrature, a 0.1%
drop rate where §13.4 says nothing fails, and a name collision created by §13.2's own rename. **Where §28 and an
earlier section disagree about a number, §28 is the measurement and the earlier section is the
projection it replaces.** The projections are left in place rather than overwritten, so that what was
projected and what was measured stay distinguishable — which is §26's discipline applied to this
document's own history.

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
  **three** nodes. And 24.0% of replicates fit σ̂ above 1 — 48 of 200, and every hierarchical rate in
  this document is at n = 200 (§26) — where the non-adaptive rule is not converged, so this is not a
  refinement, it is what makes the arm computable at all.
- **§12.4, §12.8 — the fitter's shape is decided by a cost measurement, not by taste.** The same fit
  costs **39.0 s** under Newton with a central-difference Hessian and **0.76 s** under BFGS with an
  analytic gradient — 2000 replicates being 21 hours against **17.5 min**. The analytic gradient is
  therefore part of the specification and not an optimisation.
- **§12.5 — σ̂ reaches the boundary in 22.0% of replicates**, 44 of 200, and the conditional-mode Newton
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
  |     the sigma -> 0 boundary is reached in 44 of 200 replicates, 22.0%        (§12.5) |
  |                                                                                       |
  |  inference(pop, audit)               -> bootstrap.Bootstrap                  (§13)   |
  |     bootstrap.replicates(pop, body, N_BOOT, SEED, BOOT_STRATUM)                       |
  |     ALL FOUR STRATA, USZ included -- Stage 10 §12.2 left this here           (§13.1) |
  |     FOUR failure groups; a pooled FitError costs hier.* too                  (§13.2) |
  |     bootstrap.intervals(draws, lambda key: False)   NO p on ANY key          (§13.3) |
  +--------------------------------------------------------------------------------------+

  model.py gains TWO public names:   ordinal_probabilities  (§6),  polr_ri  (§12)
  bootstrap.py gains ONE:            intervals              (§13.3, the third caller)
  balance.py  gains ONE, a rename:   levels                 (§10.3, _levels made public)
  twelve `C.SchemaError` / `FitError` identifiers, and they are the T series   (§14)

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
| `extended_bridging/standardise.py` | **new.** [§14a] whole: five public names (§3.2), three returned records of its own (§3.1), sixty-nine bootstrap estimand keys (§3.3) |
| `extended_bridging/model.py` | **amended, and it gains a new estimator.** `ordinal_probabilities` — category prediction from a `PolrFit`, the counterpart of `predict` (§6); and `polr_ri` — the random-centre-intercept proportional-odds fit, adaptive Gauss–Hermite, analytic gradient, BFGS (§12). `RIFit`, its frozen record. `polr`, `firth`, `design`, `predict` and `complete_cases` are untouched, and §20.11 asserts that by number |
| `extended_bridging/bootstrap.py` | **amended.** `intervals(draws, tested)` becomes the ninth public name — the loop Stage 11 §16 left duplicated across two callers, this stage being the third (§13.3). Both Stage 11 callers move onto it. `resample`, `replicates`, `percentile_ci` and `bootstrap_p` are unchanged |
| `extended_bridging/balance.py` | **amended, and it is a rename with no behaviour change.** `_levels` becomes public `levels` — the declared-factor-level row rule, which §10.3 needs and may not own a second copy of. `balance.py:292-295` already names this stage as why its missingness branch exists (§10.3, §18) |
| `extended_bridging/config.py` | **amended.** `POLR_RI_NODES`, `POLR_RI_SIGMA_FLOOR`, `POLR_RI_MAX_ITER`, `SUPPORT_COVARIATES` and **four** `FAILURE_BUCKETS` tokens, in the blocks §18 names (§18, fence) |
| `extended_bridging/tests/test_standardise.py` | **new.** §20 |
| `extended_bridging/tests/fixtures_stage12.py` | **new.** §20.0, and that row is the one place its contents are enumerated |
| `extended_bridging/tests/reference/polr_ri_clmm.R` | **new.** The `ordinal::clmm` oracle for `polr_ri`, gated as `polr_clm.R` is (§12.7) |
| `extended_bridging/tests/test_model.py`, `test_bootstrap.py`, `test_balance.py`, `test_config.py`, `test_reference_r.py` | **amended.** The two new estimators, the ninth public name and the surface count that moves with it, `balance`'s surface count moving with the rename, the four new constants, the second R gate (§20) |
| `extended_bridging/implementation_roadmap.md` | **amended.** Stage 12 gains its `**Spec:**` line; five corrections and three additions (§27) |
| `TODOS.md` | **amended.** Two items close, one has its stated reason corrected while its trigger fires, **two** have their triggers fire and stay open, **nine** are new (§21, §27) |

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
   the pooled arm and the support-restricted arm, together      44 - 50 s / 2000   MEASURED (§13.4)
   the hierarchical arm                                        1050   s / 2000   PROJECTED x10 from
                                                                                 105.0 s / 200 (§12.8)
                                                              ─────────────────
   Stage 12's bootstrap, one call, three arms                  ~18.4 min          part measured
```

**The hierarchical line is a projection and the line above it is not**, which §12.8 measures and
§25's step 0 closes before implementation begins.

**And the hierarchical arm is 95% of it**, which is the number that makes §12.4's optimiser choice a
specification and not an optimisation: the same arm under the first prototype this document wrote is
**21 hours** (§12.8).

---

## 3. Module shape

### 3.1 What the stage returns

**Three records of this stage's own — `Standardisation`, `Support`, `Hierarchical` — plus `RIFit`,
which is `model.py`'s and is given here because §12 specifies it.** `inference` returns a fourth
thing, `bootstrap.Bootstrap`, which this stage does not define. So "four" means different sets in §1
and in §3.1 unless it is said which, and this is where it is said.

All frozen, none carrying a verdict, following Stage 8 §3's rule that a record reports what it
computed and a caller decides what it means.

```python
@dataclass(frozen=True)
class Standardisation:
    """One [§14a] standardisation: one fit, one averaging population, two regimes."""

    population: str                      # "all eligible" | "treated support"          — §11
                                         #   | "all eligible, random centre intercept" — §12.6
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
    fit: model.PolrFit | model.RIFit     # RIFit for the hierarchical arm ONLY — §12.6
    dropped: tuple[str, ...]             # design columns dropped as constant — §5.3
```

**`fit` is a union of two types and `population` has three values, because all three arms return this
record and one of them is fitted by a different estimator.** The alternative — a fourth field
discriminating `"polr"` from `"polr_ri"` — is declined for §7.3's reason: it would carry what the
`fit` type already determines, and two things that encode one fact are two things that can disagree
after an edit. A caller that must branch narrows on the type.

**For the hierarchical arm, `Standardisation.fit` IS `Hierarchical.fit`** — the same object, not an
equal one — and §20.12 asserts `is` rather than `==`. Without that the `RIFit` would be reachable by
two paths with nothing requiring them to agree, which is the duplication this record's shape exists
to avoid.

**`conditional_log_odds` means something slightly different in the hierarchical arm and the
difference is stated rather than smoothed over.** In the pooled and support arms it is conditional on
`X`; in the hierarchical arm it is conditional on `X` **and on `b_c`** (§12.1). `measure` stays
`_CONDITIONAL` for all three — §8's guard is that neither may ever be labelled marginal, and both
qualify — but Stage 14's caption must carry which conditioning it is (§17.2).

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

### 3.3 The estimand keys — sixty-nine, and the split is [§14a]'s own

**Standardised quantities are per-ARM; model parameters are per-FIT.** That is not a tidying
principle, it is [§14a]'s structure read literally, and §7.3's identity discipline applied to the key
set instead of to two numbers:

```
   rd_0 … rd_5                          len(C.MRS_THRESHOLDS)                        6
   mrs_0_2, mortality                   [§14a]'s two namings — §7.3                  2
   dist1_0 … dist1_6, dist0_0 … dist0_6 2 * len(C.MRS_LEVELS)                       14
                                                                                   ───
                          THE STANDARDISED BLOCK, per arm                           22

   arms:  ""  (pooled)   "support."   "hier."                                22 * 3 = 66

   one CONDITIONAL LOG-ODDS PER FIT, and there are two fits                           2
      beta        the pooled proportional-odds model      [§14a]: exp(β) "may appear
      hier.beta   the random-intercept model               as a model parameter"
   plus:  hier.sigma                    the between-centre SD — §12.5                 1
                                                                                   ───
                                                                                     69
```

**There is no `support.beta`, and its absence is [§14a]'s and not an economy.** [§14a]'s Report
sentence asks for the two standardised distributions and, *from them*, the `RD_k`, the mRS 0–2 risk
difference and the mortality difference; `exp(β)` is granted a separate and narrower permission — it
*"may appear as a model parameter but never as the standardised marginal effect."* And [§14a]'s
sensitivity restricts **the standardisation population**, which §11 implements by restricting the
averaging set and not the fit. So the support arm does not have a `β`: there is one pooled fit and
therefore one pooled `β`, and a `support.beta` would be that same float — bit-identical, draw for
draw, with a bit-identical interval — presented as a property of a *standardisation*, which is the
one reading [§14a] forbids in as many words. `hier.beta` exists because the random-intercept arm **is
a second fit** (§12.1).

**The counts are computed from `C.MRS_LEVELS` and `C.MRS_THRESHOLDS` and never written as 6, 14 or
22.** Stage 8 §12 made the level set derived from the plausible range so the two cannot disagree;
this stage's key set is the next link in that chain, and `test_standardise.py` §20.3 asserts the total
against the config rather than against 69.

**The key NAMES `dist1_*` and `dist0_*` are fixed literal strings** — the wire format Stage 14 reads —
and are deliberately **not** derived from `C.TREATMENT_LABELS`, even though §7.1 forbids writing the
arm *codes* as 1 and 0 anywhere in the arithmetic. The two rules do not conflict: a code is a value the
model sees and must come from the registry, while a key name is an identifier a downstream reader
matches on, and deriving it from the registry would let a label edit silently rename a reported key.

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

**The guard is therefore at key granularity, and it is T12** (§14): a fit reaching
`POLR_MAX_ABS_BETA` on the treatment coefficient drops that fit's `beta` key and **keeps** the
standardised keys, which are unaffected and correct. Stage 10 §7 already supports per-key failure —
its `_replicate` costs one binary outcome's three keys without touching the others — so this needs no
new mechanism, only the declaration of a third granularity.

**It applies to BOTH fits, and the second one is not an extension of the argument but the same
argument.** The bound protects a quantity this stage reports as a *model parameter*, and there are two
such quantities because there are two fits (§3.3): `beta` from the pooled model and `hier.beta` from
the random-intercept one, each with its own percentile interval. A separated random-intercept fit
drives `hier.beta` to infinity exactly as a separated pooled fit drives `beta`, and nothing about the
random intercept changes that. So:

- **the pooled call site costs `beta`** — and *only* `beta`, because there is no `support.beta`
  (§3.3): the support arm restricts the averaging population and not the fit, so it has no `β` of its
  own to lose;
- **the `polr_ri` call site costs `hier.beta`**, leaving `hier.*`'s twenty-two standardised keys and
  `hier.sigma` intact, which is the same trade §9.1 makes for the pooled arm.

**One private helper, `_assert_beta_reportable`, two call sites, and the KEY COST is the caller's
decision rather than the helper's.** The helper knows a coefficient and a bound; it does not know
which keys a caller will lose, and a helper that did would be one function with two callers agreeing
on nothing but the comparison — §4.1's objection to the shared [§11] mask, applied here.

**`outcome._assert_reportable` is NOT reused, and §24 records why.** It bounds the *worst* coefficient
in the fit (`outcome.py:363`, `if worst >= C.POLR_MAX_ABS_BETA`), which is right for Stage 8's
one-column design where the worst coefficient *is* the treatment coefficient, and wrong here: §9.2
requires `gamma` unbounded, and the largest `|gamma|` measured over 2000 replicates is 3.21 against a
`|beta_treatment|` of 2.00. Calling it would drop replicates for a reason the reported quantity does
not care about.

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

**And it calls `balance.levels` for the rows themselves, which is the same note written twice.**
`_levels` (`balance.py:282`) is renamed public at this stage (§18) rather than reimplemented, for
exactly Stage 7 §5.5's reason and on evidence from `_levels`' own docstring, which keeps a branch
alive *"because Stage 12 will need it over a population `complete_cases` did not build"*
(`balance.py:292-295`). A second implementation of `smd` would be a second definition of the
yardstick; a second implementation of `levels` would be a second definition of **what a row is**.

The rule it carries is *range over declared `C.FACTOR_LEVELS`, never over the frame and never over a
design matrix* — Stage 6 §7.5's, because *"a level nobody in the cohort has is a row reading 0, which
is a fact about the population; a level that vanishes is a fact nobody sees."* **That rule may not
exist twice, and the reason is that its failure is invisible here**: §10.1 measured that all three of
`onset_type`'s declared levels are present in the treated arm, so a copy that ranged over *observed*
levels would produce a byte-identical table on this workbook and silently drop a row on the next.
Nineteen rows over `BALANCE_SET` — twelve non-categorical names, `center`'s four declared levels and
`onset_type`'s three — and the count is `levels`' to determine, not this section's to assert.

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

This arm's keys carry the `support.` prefix and there are **22** of them (§3.3) — the standardised
block and no `β`. **There is no `support.beta`, and this section is why**: the arm restricts the
averaging population and leaves the fit alone, so its treatment coefficient *is* the pooled one,
and reporting it under a `support.` prefix would present `β` as a property of a standardisation —
[§14a]'s forbidden reading (§3.3, §8). `Standardisation.population` carries `"treated support"` so
the record says which it is without a caller tracking the prefix.

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
| 2.0 | 5 | 4.0e-1 | 3.1e-6 |
| 2.0 | 7 | 3.5e-1 | 4.4e-7 |
| 2.0 | 9 | **4.6e-1** | 1.0e-7 |
| 2.0 | 31 | 4.5e-2 | 0 |
| 3.0 | 3 | 1.9e0 | **8.9e-4** |
| 3.0 | 7 | 1.1e0 | 2.4e-6 |
| 3.0 | 31 | **4.8e-1** | 0 |

*(Every value the prose below quotes is a row above. The σ = 2 rows at 5, 7 and 9 nodes are the
non-monotonicity in point 2; the σ = 3 row at 3 nodes is the Status block's "adaptive is off by
8.9e-4 at three nodes". A number quoted from a table it is not in is a number nobody can check.)*

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
**24.0% of replicates fit σ̂ above 1** — 48 of 200 — squarely in the region where the non-adaptive rule is not
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

### 12.5 The σ → 0 boundary is reached in 22.0% of 200 replicates, and it needs a floor

[§14a] anticipates the boundary: *"if outcomes truly do not differ by centre given `X` then σ²_C = 0
and it collapses to the pooled model."* It is not a hypothetical. Measured over 200 replicates:

```
   MEASURED AT 200 REPLICATES, NOT AT N_BOOT.  Every rate below has n = 200.
   sigma_hat    min 0.000021   median 0.65923   max 1.78734
   below 0.01, the variance collapsing to the pooled model     44 / 200   22.0%
   above 1.0                                                   48 / 200   24.0%
   converged                                                  200 / 200   all on "likelihood"
   iterations   min 23   median 27   max 47      (POLR_RI_MAX_ITER 200)
```

**These were the hierarchical arm's only measurements and they are at n = 200**, while the pooled and
support arms were run at the full `N_BOOT = 2000` (probe B1). The asymmetry is a fact about how this
document was produced and not about the estimator, so it is labelled everywhere the rates are
repeated — §17.2, §19 and §21 item 7 — rather than stated once here and forgotten.

**§25's STEP 0 HAS SINCE RUN AND §28 CARRIES THE `N_BOOT` FIGURES: 18.4% at the floor and 23.0% above
σ = 1**, against the 22.0% and 24.0% above. Both are within sampling error of the 200-replicate
estimates — 18.4% against 22.0% is 1.2 standard errors at n = 200 — so these figures were imprecise in
the way this section said they were rather than wrong. **The rates a manuscript prints are §28's.**

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
   one fit, point estimate                          MEASURED       0.76 s
   200 replicates, end to end                       MEASURED     105.0 s     525 ms each
   at N_BOOT = 2000                                 PROJECTED x10 17.5 min
   the same arm under the first prototype (C2)      PROJECTED     21.7 hours
```

**The 17.5 minutes is a ×10 projection from a measured 200 and not a measured 2000, and this section
is the only place that says so — so every other section that quotes it says so too** (§2, §13.4, §19).
The two-arm figure beside it in §2 *is* a measured 2000 (probe B1); the hierarchical arm's is not.
§25's step 0 closes that gap before any code is written, and §26 records it under what was not
measured.

The projection is from a measured 200-replicate run of the full body — `model.design`, `model.polr`
for the start values, then `polr_ri` — and not from the point fit's 0.76 s, because the start values
matter: `polr_ri` starts from the pooled fit's `alpha` and `beta` with `sigma = 0.5`, which is a
nested model's exact maximiser in every coordinate but one, and it is why 25 iterations suffice.
`polr`'s docstring makes the same argument for its own start values.

**That dependency is load-bearing at [§10] and §13.2 states its consequence**: because `polr_ri`'s
start values *are* the pooled fit's, a replicate whose `model.polr` raises has no hierarchical arm
either, and a pooled `FitError` therefore costs every key rather than exempting `hier.*`.

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
     FitError   -> costs EVERY key, hier.* INCLUDED             polr_ri STARTS from this fit
     T12        -> costs `beta` ONLY                            §9.1, key granularity
     pooled     = standardise(...)                              22 keys
     support    = standardise(..., over=in_support(draw))       22 keys, box RECOMPUTED  §11
     ri         = model.polr_ri(X, y, draw[BOOT_STRATUM])       start values ARE `fit`'s  §12.8
     FitError   -> costs hier.* ONLY                            24 keys: 22 + hier.beta + hier.sigma
     T12        -> costs `hier.beta` ONLY                       §9.1, the second call site
     hier       = standardise(..., b_hat)                       §12.6
```

**A pooled `FitError` costs the hierarchical arm too, and the reason is §12.8's and not a
conservatism.** `polr_ri` starts from the pooled fit's `alpha` and `beta` with `sigma = 0.5` — a
nested model's exact maximiser in every coordinate but one, which is why 25 iterations suffice and
why the arm costs 17.5 minutes rather than hours. **There is no hierarchical arm without a pooled
fit**, so an exemption for `hier.*` would require a cold start whose iteration count, convergence
route and boundary rate are all unmeasured, and a quarter of the draws would then sit in a different
numerical regime from the rest — which is precisely the failure §12.2 rejects non-adaptive quadrature
for. Measured: the pooled fit fails in **0 of 2000** replicates (§13.4), so this path never fires on
this workbook, and it is specified anyway on §5.3's rule.

**The granularity is Stage 10 §7's, extended by one level and not reinvented.** Stage 10's
`_replicate` already has three: a propensity failure costs every key, a primary failure costs `beta`
and the six `rd_k` together, a per-outcome failure costs that outcome's three. Stage 12's are the same
idea over its own dependency graph, and the rule is the one Stage 10 states — a failure costs exactly
the keys that could not be computed and no others — with §9.1's key-granularity case as the one
addition. **Four failure groups**, not three: T1 and a pooled `FitError` both cost every key but carry
different buckets (`degenerate_design`, `nonconvergence`), a `polr_ri` `FitError` costs `hier.*`, and
T12 costs one `β` at whichever fit raised it.

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

**Two runs, and they are reported as two, because no single run produced this whole picture.**

```
   THE POOLED AND SUPPORT ARMS  --  measured at N_BOOT = 2000        [probe B1]
     replicates attempted                                  2000
     survived                                              2000       100.00%
     failure buckets                                       none
     T1 (design lost the exposure)                            0
     T12, separation (POLR_MAX_ABS_BETA)                      0       §9.2
     nonconvergence, pooled                                   0
     design columns dropped                                   0       100.00% empty

   THE HIERARCHICAL ARM  --  measured at 200, PROJECTED to 2000       [probe F4]
     replicates attempted                                   200
     survived                                               200       100.00%
     nonconvergence, polr_ri                                  0       §12.5
     T2, T3, T4                                               0

   THE HIERARCHICAL ARM  --  MEASURED at N_BOOT = 2000                [§28, step 0]
     replicates attempted                                  2000
     survived                                              1998        99.90%
     nonconvergence, polr_ri (T2, step-halving exhausted)     2         0.10%
     cost:  hier.* ONLY -- 24 of the 69 keys, per replicate
```

**Nothing fails in the pooled and support arms, and the hierarchical arm drops 0.1%.** The second half
is §28's measurement and it amends this section: at n = 200 the arm lost nothing, and at `N_BOOT` two
replicates exhaust the halvings — at iteration 22 and at iteration 7 — in a genuine BFGS line-search
stall rather than in the `nan` artefacts §28 fixed. Each costs `hier.*` and nothing else, so the
pooled and support arms keep 2000 draws while the hierarchical arm has 1998, every `Draws` reconciles,
and the failures land in `nonconvergence`. **The taxonomy handles it as designed; what is open is
whether a 0.1% drop selected by the OPTIMISER rather than by the data is exchangeable**, which §28
files as a `TODOS` item rather than answering.

**The pooled arm's own claim is unchanged and was re-measured:** That is worth stating plainly against the run of this pipeline's history: Stage 10
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

## 14. The failure taxonomy — the T series, and there are twelve

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
| **T12** | `_assert_beta_reportable` | `FitError` | a fitted treatment coefficient reached `POLR_MAX_ABS_BETA` (§9.1) |

**`bootstrap._bucket` classifies by the FIRST TOKEN of the message** — `message.split()[0]`
(`bootstrap.py:487`) — **and raises `C.SchemaError` on a token it does not know.** Every rule below is
a consequence of that one line and not a style choice:

- **T2 and T3 lead with `polr_ri:`**, and that token goes in `C.FAILURE_BUCKETS` mapping to
  `nonconvergence` — `"polr:"`'s arrangement, unchanged.
- **T4 leads with `T4`, NOT with `polr_ri:`.** T4 is a degenerate design and not a convergence
  failure, and a message leading with `polr_ri:` would be counted as `nonconvergence` **while
  `"T4" → "degenerate_design"` sat in the map unreachable.** This is the one misclassification the
  raise-site scan cannot catch: the scan asserts every token it finds is *in* the map, and
  `polr_ri:` is in the map, so it would pass while the counter read wrong. T1 and T12 lead with their
  own identifiers for the same reason.

`test_bootstrap.py`'s scan asserts every `FitError` raise-site token in `model.py` and `outcome.py` is
a `FAILURE_BUCKETS` key; **the scan scope grows to include `standardise.py`** (§18), so T1's and T12's
tokens are covered by the same mechanism and a reworded message is a test failure rather than a
counter that silently reads zero. **The scan remains a mitigation and not a fix** — it catches a
reworded token and not a wrong-but-mapped one — which is `TODOS`' `code`-field item, whose trigger
this section fired (§21).

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

**Their tolerances and their check sites are prespecified, because a self-check with an invented
tolerance is a self-check that means whatever its implementer decided:**

| | fires when | tolerance is | checked |
|---|---|---|---|
| **T10** | `abs(sum(P) − 1) > 1e-12` for any row or arm, **or** the cumulative sequence is not non-decreasing | §20.5's, so the guard and the test cannot disagree | every `Standardisation` the module constructs — point estimate **and** every replicate |
| **T11** | `mrs_0_2 != rd[2]` exactly, **or** `abs(mortality + rd[5]) > 1e-15` | §20.6's, and the asymmetry is §7.3's | every `Standardisation`, same as T10 |

The two tolerances differ for the reason §7.3 gives and it is not a rounding choice. `mrs_0_2` **is**
`rd[2]` — the same cumulative difference read twice, so the comparison is exact and `1e-15` would
hide a real defect. `mortality` is computed independently from the distributions (§20.6 asserts it
by monkeypatch), so it agrees with `−rd[5]` to floating-point and not to the bit; the measured worst
is 8.3e-17, three orders below the bound.

**Both are `SchemaError`, so both are uncaught inside a replicate and kill the whole bootstrap.** That
is deliberate and it is §14's rule, not an oversight: a distribution that does not sum to 1 is a
contract break in the arithmetic, not a sparse resample, and eighteen minutes of draws built on it
are worth less than the crash. The measured margins — 6.7e-16 against 1e-12, 8.3e-17 against 1e-15 —
are what make that safe to specify.

---

## 15. The audit entries, and there are ten

All under existing `data.KINDS`; **no kind is added**, which is the first stage since Stage 5 of which
that is worth saying, and the reason is that `cohort` and `model` between them already describe what
this stage does. Every step name is distinct from `cohort.build`'s for §4.2's reason.

**The numbering below is the order the entries are RECORDED, and that order is a consequence of
§19's call order rather than a declaration made here.** In particular entry 7 falls between 6 and 8
because `support` must run before the restricted `all_centre` that produces it — `over=sup.inside` is
`support`'s output. §20.9 asserts the order against §19's driver and not against this table.

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
    restriction 1 removes a centre from [Stage 5 §9]. T6, T7. Records entries 1 and 2.
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
    """The [§14] bootstrap: one `replicates` call, three arms per draw, sixty-nine keys, no p-value.

    ALL FOUR STRATA (§13.1). Records entry 10. Calls `bootstrap.replicates`, `bootstrap.bucket`,
    `bootstrap.collect` and `bootstrap.intervals`; implements none of them.
    """
```

Five functions, and **`inference` recomputes the three arms inside its body rather than taking the
point-estimate records**, for Stage 10 §11's reason: a replicate must run the whole procedure, and a
body that took a fitted object would be resampling around a fixed fit.

---

## 17. Data flow into Stages 13 and 14

**The call order is §19's and this block does not restate it in a different one.** `support` runs
before the restricted `all_centre`, because `over=sup.inside` is *its* output:

```
  population(df, audit)   -> DataFrame, 104 rows                       (§4.1)
  all_centre(pop, audit)  -> Standardisation   n_average 104           (§7)
  support(pop, audit)     -> Support           12 never-IVT, 3 outside (§10)
  all_centre(pop, audit, over=sup.inside)
                          -> Standardisation   n_average  93           (§11)
  hierarchical(pop, audit)-> Hierarchical      sigma, b_hat, at_floor  (§12)
  inference(pop, audit)   -> Bootstrap         69 keys, 2000 draws each, NO p
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
- **The hierarchical arm prints its `at_floor` rate and USZ's posterior SD.** **18.4% of replicates
  collapse to the pooled model — 368 of 1998, MEASURED at `N_BOOT` (§28)**; the 22.0% from n = 200 that
  this line first carried is the projection §25's step 0 replaced — and the interval on `hier.sigma` has a limit that means
  "the boundary" in those draws (§12.5); USZ's intercept is the least precise of the four and the arm
  rests on it (§12.6).
- **There is one conditional log-odds per FIT, never one per arm.** The support arm has no
  odds-ratio row: it restricts the averaging population and not the fit, so its `β` would be the
  pooled `β` under a sensitivity label, which is [§14a]'s forbidden reading (§3.3). Stage 14 prints
  `beta` for the pooled model and `hier.beta` for the random-intercept one, and nothing for
  `support.` — a key that does not exist (§20.3).
- **The two `β`s are conditional on different things, and the caption must say which.** `beta` is
  conditional on `X`; `hier.beta` is conditional on `X` **and on `b_c`** (§3.1, §12.1). Both carry
  `measure = _CONDITIONAL`, which is §8's guard that neither may be labelled marginal, but a table
  putting them in one column under one heading would assert they are the same quantity.
- **No [§14a] quantity carries a p-value, and none is in a [§13] family.** §3.3 is why, and Stage 11's
  Benjamini–Hochberg must not acquire a fifth family by counting Stage 12's keys — which is the
  mistake Stage 10 §12.1 warned about in its own direction ("seven p-values, not fourteen").

---

## 18. What Stage 12 amends in Stages 1–11

The full ledger, so that no amendment is discovered during implementation. **Four shipped modules and
five test modules; the two new estimators are both in `model.py`; `balance.py` changes by one
rename and nothing else; and nothing in `outcome.py`, `propensity.py`, `cohort.py`,
`eligibility.py`, `derive.py` or `data.py` changes at all.**

| File | Amendment | Why |
|---|---|---|
| `model.py` | `ordinal_probabilities` (§6) and `polr_ri` + `RIFit` (§12). Public surface 6 → 8 | §0.1: an estimator-shaped function lives in `model.py`. `polr`, `firth`, `design`, `predict`, `complete_cases` and `Fit`/`PolrFit` are untouched, and §20.11 asserts it by number |
| `bootstrap.py` | `intervals(draws, tested)` becomes the ninth public name; `run` and Stage 11's two callers move onto it | §13.3. The `TODOS` trigger is a third caller and this is it. `run`'s twenty-six intervals are asserted byte-identical across the change (§20.11) |
| `balance.py` | `_levels` becomes public as **`levels`** — a rename and nothing else, no behaviour change. Public surface 3 → 4 | §10.3. `balance.py:292-295` already names this stage as the reason its missingness branch exists; Stage 7 §5.5 made `smd` public on the same argument, and the rule it carries may not exist twice (§24) |
| `config.py` | `POLR_RI_NODES`, `POLR_RI_SIGMA_FLOOR`, `POLR_RI_MAX_ITER` in the model-tolerance block; `SUPPORT_COVARIATES` as a computed view in the covariates block; **four** `FAILURE_BUCKETS` tokens | §12.3, §12.5, §12.4, §10.1, §14. The fence below |
| `tests/test_model.py` | `ordinal_probabilities` against a hand-built fit and against `_ord_pieces`; `polr_ri` against its analytic gradient, its node-count stability and its floor | §20.5, §20.12 |
| `tests/test_bootstrap.py` | `intervals`' three p-rules; the surface count 8 → 9; the raise-site scan scope grows to `standardise.py` | §13.3, §14 |
| `tests/test_balance.py` | the surface count 3 → 4; `levels`' behaviour is unchanged, asserted against the existing Stage 7 fixtures rather than by diff | §10.3 |
| `tests/test_config.py` | **six** assertions: `POLR_RI_NODES` odd and ≥ 9; `0 < POLR_RI_SIGMA_FLOOR < 1e-2`; `STANDARDISATION_COVARIATES == tuple(c for c in PS_COVARIATES if c != "center")` — [§15]'s permission as a static check (§5.1); `SUPPORT_COVARIATES` is `STANDARDISATION_COVARIATES` minus `CATEGORICAL`; the four new `FAILURE_BUCKETS` tokens; and `"support.beta"` is in no key set (§3.3, §20.3) | §20.2 |
| `tests/test_reference_r.py` | the `clmm` oracle, under the existing `reference_r_ordinal` gate and `_r_environment` | §12.7 |
| `implementation_roadmap.md` | Stage 12 gains its `**Spec:**` line; five corrections and three additions | §27 |
| `TODOS.md` | two close, one has its stated reason corrected, two have their triggers fire and stay open, nine are new | §21, §27.1 |

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

# The sigma -> 0 boundary is REACHED -- 44 of 200 replicates, 22.0%, measured (Stage 12 §12.5) -- and under a
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

**Four** `FAILURE_BUCKETS` tokens are added beside the existing sixteen: `"T1"` →
`"degenerate_design"`, `"polr_ri:"` → `"nonconvergence"`, `"T4"` → `"degenerate_design"`, `"T12"` →
`"separation"`. `"T2"` and `"T3"` both lead with `polr_ri:` and share a bucket, which is `"polr:"`'s
arrangement and Stage 10 §7.2's argument for it, unchanged.

**T12's bucket is `"separation"` and not a new one**, which is Stage 8 §11's requirement that the
separation count be reported separately from the convergence count, satisfied by reusing G7's bucket
rather than by adding a second name for the same thing. **T4's token is `"T4"` and not `"polr_ri:"`,
and §14 is why** — a token that is in the map but is the wrong one passes the raise-site scan and
misreads the counter.

---

## 19. Handover to Stage 13 [§14b]

**This is the canonical call order — the ONE place it is written, and §15, §17 and §20.9 reference it
rather than each implying one.** `support` must precede the restricted `all_centre`, whose `over=` is
its output; the six calls record audit entries 1–10 in the declared order, entry 7 falling between 6
and 8 because that is where the call falls:

```
  pop  = standardise.population(classified, audit)                    104 rows, all 4 centres
                                                                      entries 1, 2
  std  = standardise.all_centre(pop, audit)                           [§14a]   entries 3, 4
  sup  = standardise.support(pop, audit)                                       entries 5, 6
  ssup = standardise.all_centre(pop, audit, over=sup.inside)          §11      entry  7
  hier = standardise.hierarchical(pop, audit)                                  entries 8, 9
  boot = standardise.inference(pop, audit)               69 keys, 2000 draws   entry  10

    126 -> 107 -> 104   the ledger, and the 2 outcome-missing are BOTH at USZ     [§4.1]
    12 / 3              never-IVT patients, and how many are outside the box      [§10.2]
    4.15%               replicates losing mRS 5, in which RD_5 == RD_4 exactly    [§6.3]
    0 / 2000            failures in the POOLED and SUPPORT arms                   [§13.4]
    368/1998, 460/1998  polr_ri replicates at the sigma floor / above sigma = 1   [§28]
    18.4%, 23.0%        MEASURED at N_BOOT.  The 44/200 and 48/200 this line
                        first carried were probe F4's, at n = 200                 [§12.5]
    2 / 2000            polr_ri non-convergences, costing hier.* only             [§28]
    ~4.2 min            the bootstrap, of which the hierarchical arm is ~95%      [§28]
                        MEASURED; the 18.4 min this line first carried was a
                        x10 projection from n = 200                               [§12.8]
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

Six assertions, listed in §18. The load-bearing one is
`STANDARDISATION_COVARIATES == tuple(c for c in PS_COVARIATES if c != "center")`, which makes [§15]'s
*"permitted inside §14 only"* a static check rather than a sentence: a covariate added to the [§14a]
model without being a [§6] confounder fails at import.

### 20.3 The estimand keys are derived and not counted

`len(keys) == 3 * (len(C.MRS_THRESHOLDS) + 2 + 2 * len(C.MRS_LEVELS)) + 2 + 1`, asserted against the
config and never against 69 — three standardised blocks, one `β` per fit and there are two,
`hier.sigma`. And `all(interval.p is None for interval in boot.intervals.values())` — §3.3's no-p
rule, asserted over every key rather than spot-checked.

And **`"support.beta" not in keys`**, asserted by name. §3.3's argument is [§14a]'s, so the absent key
is the enforcement in the manner §8's absent `odds_ratio` field is: a later "make the arms uniform"
edit fails a test rather than putting `β` into a row labelled as a sensitivity result.

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

**The two tolerances differ, and the monkeypatch test applies to only one of them, and both facts are
§7.3's rather than arbitrary.** `mrs_0_2` **is** `rd[2]` — the same cumulative difference read under
[§14a]'s name for it, one computation and not two — so the comparison is **exact**, `1e-15` would
hide a real defect, and there is nothing for a monkeypatch to catch. `mortality` is genuinely a
second computation, from `P(Y = 6)` in each arm rather than from `1 − P(Y ≤ 5)`, so it agrees with
`−rd[5]` to floating point and not to the bit — measured 8.3e-17 — and it is exactly the quantity an
implementer might "simplify" into `-rd[5]`, which is what the monkeypatch prevents. T11 fires at these
same two tolerances (§14).

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

Ten entries, in the order **§19's driver produces them** — the test runs that driver rather than
calling the five functions in an order of its own, so entry 7 lands between 6 and 8 because `support`
precedes the restricted `all_centre` (§15). The step-name set is **disjoint** from `cohort.build`'s,
asserted as a set intersection rather than against a list (§4.2); entry 1 names exactly as many case
identifiers as its `n` claims, `cohort._record_removal`'s discipline; and the Stages 1–12 ledger
totals 56.

**A driver that omitted the restricted `all_centre` would record nine entries and total 55**, which
is the failure this assertion catches and the reason §19 spells the call out rather than leaving
sensitivity 1's point estimate implied.

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
  **T4's raised message begins with the literal `T4` and NOT with `polr_ri:`**, asserted on the
  message text, because `bootstrap._bucket` reads the first token and `polr_ri:` is also a valid key
  — so this is the one bucket assertion the raise-site scan cannot make (§14).
- **T12 at the `polr_ri` call site**: a constructed frame driving `hier.beta` past
  `POLR_MAX_ABS_BETA` costs `hier.beta` and **nothing else** — `hier.*`'s twenty-two standardised
  keys and `hier.sigma` all survive (§9.1).
- **`Hierarchical.fit is Hierarchical.standardisation.fit`** — identity, not equality (§3.1). And
  `Hierarchical.standardisation.population` is the third literal, so a formatter can tell the three
  arms apart from the record alone.

### 20.13 The replicate loop

- The drawn-frame sequence is byte-identical to the one a body returning `None` produces, at the same
  seed, frame and stratum — Stage 10 §15.2's property, now over a three-arm body (§13.1).
- **All four strata appear in every drawn frame** and each stratum's size is exact, USZ's 12 included.
- The **four** failure groups each cost exactly their own keys, driven by four constructed failures,
  and `Draws.n_attempted` reconciles for every one of the sixty-nine keys.
- **A pooled `model.polr` `FitError` costs EVERY key, `hier.*` included** — §13.2, because
  `polr_ri`'s start values are the pooled fit's (§12.8). Asserted by driving the pooled fit to raise
  and requiring `hier.sigma`'s surviving-draw count to fall by exactly one, which is the assertion a
  `hier.*` exemption would fail.
- **A `POLR_MAX_ABS_BETA` breach on the pooled fit costs `beta` and nothing else** — §9.1's key
  granularity, which is the one place this stage departs from Stage 8 and therefore the one that most
  needs a test. The support arm's twenty-two keys are unaffected, and there is **no `support.beta`**
  to be affected (§3.3, §20.3).

### 20.14 Every T identifier fires, and the set is the assertion

§20.1–20.13 name a test for T1, T2, T3, T4, T9, T10 and T12 as a consequence of testing the thing
each guards. **T5, T6, T7, T8 and T11 are guards nothing else reaches**, so they are asserted here —
and asserted as a **set**, on §20.9's pattern, so that a thirteenth identifier added to §14 without a
test fails rather than passing unnoticed:

    {t for t in T_IDENTIFIERS_IN_SECTION_14} == {t exercised by this module's tests}

Each fires on a constructed input, never on a mutated workbook:

| id | constructed input |
|---|---|
| **T5** | `ordinal_probabilities` called with the fitted columns **reordered**, and again with a column renamed — `SchemaError` both times, and the message names the offending column (§6.1) |
| **T6** | a frame with the `eligibility` column dropped, i.e. `classify` was not run |
| **T7** | an empty frame, and separately a frame with one arm empty — both after `eligibility.retained`, so the emptiness is the population's and not the input's |
| **T8** | `over=` as a float Series; as a boolean Series on a shuffled index; and as one whose length differs from the population's |
| **T11** | a `Standardisation` whose `mortality` has been perturbed by 1e-14 — above §20.6's 1e-15 bound and far below anything a reader would see — so the check is known to be live at its stated tolerance and not merely present |

**T8's three cases are one point about alignment and it is worth stating why they are separable.**
`over=` is applied inside every replicate to a **resampled** frame, and that frame is safe to index
by label only because `bootstrap.resample` **resets the index** before renaming — `bootstrap.py:315-320`,
which says so in as many words and gives Stage 12 as the reason: *"a frame with a duplicated index is
a frame on which a future label-based `.loc` is wrong in a way that returns a number."* T8 is the
guard for the case that reset does not cover, and the reset is the reason T8 does not have to cover
the duplicated-label case at all. Stage 10 wrote it for this stage; §11's `over=` is the caller it
was written for.

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
7. **NEW — a percentile interval on a parameter at a boundary.** **18.4% of `hier.sigma`'s draws are at
   the floor — 368 of 1998, MEASURED at `N_BOOT` (§28); the 22.0% this item first carried was 44 of
   200.** Those draws are at
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
12. **NEW — the textual failure classification produced a wrong bucket in this document, and the
    raise-site scan could not have caught it.** `TODOS`' `code`-field item sets its trigger as *"the
    first bucket that is wrong"*; this spec's first draft gave T4 a message leading with `polr_ri:`,
    which `bootstrap._bucket` (`bootstrap.py:487`) would have counted as `nonconvergence` while
    `"T4" → "degenerate_design"` sat in the map unreachable. **The scan asserts every token it finds
    is in the map, and `polr_ri:` is in the map** — so it would have passed. That is the second
    documented instance, after the `G6` omission the item already records, and both were caught by
    reading a spec rather than by a run. **Status: the trigger has FIRED and the fix is deferred** —
    Stage 12 is already landing a new estimator, a fourth public module surface and the interval-loop
    refactor, and a nineteen-site error-taxonomy amendment on top is two structural changes at once.
    §14 states the first-token rule explicitly so the next author does not rediscover it.
13. **CLOSED with measured numbers, 2026-08-27 (§28).** The hierarchical arm's rates were
    200-replicate measurements; §25's step 0 ran and replaced them. At `N_BOOT = 2000`: **18.4% at the
    floor** (368 of 1998) against 22.0%, **23.0% above σ = 1** (460 of 1998) against 24.0%, **2 of 2000
    non-convergences** against 0 of 200, and **4.2 minutes measured** against 17.5 projected. The two
    rates were imprecise as stated; the wall clock was a factor of four out, and a timing is machine
    state. **Status: closed**, which is what this item asked for — closed with measured numbers rather
    than negatively.
14. **NEW — this document has not been swept for the error class §27.1's first row is an instance
    of: a structural decision argued from engineering grounds where [§10], [§14] or [§16] gives a
    direct answer.** The estimand key set was specified as 70 keys with a uniform twenty-three per
    arm, which is a good *engineering* shape — it makes the count derivable from `config.py` — and is
    not what [§14a] asks for. [§14a] lists `exp(β)` as a model parameter separately from the
    quantities to report *from the distributions*, and defines the sensitivity on the standardisation
    population; both sentences were in the plan the whole time. **One instance was found because the
    question was asked; the sweep was not run.**

    **The check is cheap and mechanical**: for every structural decision in this document — what is
    reported, at what granularity, over which population, under which label, with which interval —
    ask whether [§10], [§14] or [§16] speaks to it directly, and whether this document cites that
    sentence or reasons around it. Candidates worth the sweep, named so they are not rediscovered:
    §7.3's choice of what to report from the distributions; §8's treatment of `exp(β)`; §10.1's
    resolution of *"each continuous covariate"*; §11's reading of *"restrict the standardisation
    population"*; §13.3's no-p rule against [§10]'s inference paragraph; §17.2's list of what must
    travel with the estimate against [§16].

    **Why it matters more here than in an ordinary spec:** this document is the sole implementation
    source (header), so a place where it reasons around the plan rather than from it becomes code
    that disagrees with the prespecified analysis — and the disagreement is invisible, because the
    code will match the spec.

    **Trigger.** Before this document is treated as final, or before Stage 13 reuses its machinery —
    Stage 13 inherits §7.1's regime seam and §16's entry points, so a miscalibration here propagates.
    **Status: open, and it is the largest remaining risk in this document**, larger than any single
    number in §26, because it is the one class of error the review that produced §27.1 is known to be
    bad at.

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
| The balance row list and labels | `balance.levels` (public at this stage), `C.BALANCE_SET` | **Call it.** Its rule — range over declared levels, never over observed — is what §10.3 needs, and `balance.py:292-295` already names this stage as why its missingness branch exists. Renamed public here for the reason Stage 7 §5.5 made `smd` public (§10.3, §18) |
| The treatment-coefficient bound | `outcome._assert_reportable` (`outcome.py:347`) | **Do NOT call it.** It bounds the *worst* coefficient in the fit (`outcome.py:363`), which is right for Stage 8's one-column design and wrong here: §9.2 requires `gamma` unbounded, and measured `max\|gamma\|` is 3.21 against `max\|beta_treatment\|` 2.00. T12 is this stage's own, on the treatment coefficient alone (§9.1) |
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

0. **Re-measure the hierarchical arm at `N_BOOT`, before any code is written.** Re-run the
   three-arm replicate body at `N_BOOT = 2000` on the landed Stages 1–10 and replace §12.5's,
   §12.8's, §13.4's, §2's, §17.2's, §19's and §21 item 7's projected figures with measured ones —
   the at-floor rate, the above-σ=1 rate, the `polr_ri` non-convergence count and the wall clock.
   About **eighteen minutes** of compute. It is step 0 and not step 10 because §17.2 requires the
   at-floor rate printed beside `hier.sigma`'s interval, so it is a number that reaches a
   manuscript, and because a document whose authority rests on every number having been run should
   not carry a projection into implementation when closing it costs eighteen minutes. **Until it
   lands, every hierarchical rate in this document means n = 200** (§26).

   **It must be run where the workbook is.** `DATA_XLSX` is gitignored, so in a git worktree or a
   fresh clone it is absent and every workbook-dependent test skips — measured in this branch's
   worktree, **1343 passed / 112 skipped** against the main checkout's 1588 / 2, with the 110
   additional skips all reading *"the private workbook is gitignored and absent from this
   checkout"*. That is not a failure and nothing is wrong with the suite; it is a reason step 0 and
   every probe re-run belong in the checkout that has the data, and a reason a green run in a
   worktree is not evidence about any [§14a] number.
1. **`config.py`** — the four constants and the four `FAILURE_BUCKETS` tokens (§18 fence).
   `test_config.py`'s six assertions.
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
9. **`standardise.inference`** (§13) — the three-arm body, four failure groups, sixty-nine keys.
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
| B1 | 2000 replicates, the POOLED and SUPPORT arms only | 2000/2000 survived, no failures, 43.9 and 50.1 s over two runs. **The hierarchical arm was NOT in this run** — it is F4's, at n = 200 | §13.4, §2 |
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
| F1–F4 | the specified shape, **at n = 200 throughout** | frozen-node gradient 1e-8 relative; 0.76 s per fit; 200/200 converge; **44/200 = 22.0%** at the floor, **48/200 = 24.0%** above σ=1; **17.5 min** at `N_BOOT`, **projected ×10** from 105.0 s / 200 | §12.4, §12.5, §12.8 |

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

**What was NOT measured, and is therefore not claimed. There are two things, and the second was found
by the engineering review rather than by the probe suite.**

**First — the hierarchical arm was never run at `N_BOOT` when this document was written.** Probe B1's
2000 replicates covered the pooled and support arms; the hierarchical arm is probe F4's, at **200**. So
`22.0%` at the floor, `24.0%` above σ = 1 and `0/200` non-convergence were 200-replicate measurements,
and `17.5 min` was a ×10 projection from a measured 105.0 s. Every section quoting them says so
(§12.5, §12.8, §13.4, §17.2, §19, §21 item 7). Nothing about the arm's behaviour was *expected* to
change at 2000 — the projection is linear in the replicate count and the body is identical — but "not
expected to change" is not a measurement, and this document did not report those as the same thing.

**§25's STEP 0 HAS SINCE RUN (§28) AND ONE OF THE FOUR FIGURES DID CHANGE MATERIALLY.** The two rates
moved within sampling error (18.4% and 23.0%), the wall clock was out by a factor of four (4.2 min
measured), and **the non-convergence count moved from 0 to 2 — which is not a precision question but a
behaviour the 200-replicate run did not exhibit at all.** That is the argument for this paragraph
existing, made against this document by its own implementation: "not expected to change" was right
about three figures and wrong about the fourth.

**Second — the `ordinal::clmm` oracle has not been run:**
R segfaults on this machine at startup for the `uname`/PATH reason `test_reference_r.py`'s banner
documents, and opening that gate is an implementation task (§25 step 7) and not a spec task. So
`polr_ri` is verified here against its own analytic gradient, against node-count stability, against the
σ = 0 collapse to `model.polr`, and against nothing external. §12.7 specifies the oracle; §21 item 11
is not it. **This is the one estimator in this pipeline whose specification ships without an
independent implementation having agreed with it**, and that is stated plainly rather than left to be
discovered at step 7.

### 26.1 The engineering review — what it checked, and what it did NOT

The review of 2026-08-27 (§27.1) is part of this document's provenance, so its **scope is recorded
with the same discipline as the probes'**. A reader who knows what was reviewed also needs to know
what was not.

**What it checked.** Internal consistency, end to end: every count against its own arithmetic, every
number in prose against the table it cites, every guard against an identifier and a bucket token,
every entry point against the audit ledger, every T identifier against an acceptance criterion, and
every claim this document makes about landed code against that code. The last of those is why it
found T4's bucket — `bootstrap.py:487` was read, not assumed.

**Three things it did NOT check, and none of them is a formality:**

1. **The statistical premises, which were taken as given because they are measured or are [§14]'s.**
   Adaptive Gauss–Hermite at eleven nodes (§12.3); the per-covariate min–max box as the definition of
   "outside the treated support" rather than a trimmed or hull-based one (§10.1); resampling the USZ
   stratum (§13.1); conditioning the hierarchical standardisation on `b̂_c` rather than marginalising
   (§12.6); `N_BOOT = 2000`. Each is either a measurement this document reports or a question §22
   hands to the PI. **An engineering review is the wrong instrument for all five**, and the right
   second reader is a statistician or the PI, not another code reviewer.
2. **No cross-model review was run.** The intended independent pass failed to authenticate and no
   substitute was run, so **this document has been reviewed once, by one reviewer**. That is recorded
   rather than left to be assumed from the review's thoroughness. §21 item 14 is the residual.
3. **Nothing was re-executed.** The review re-derived every number's *internal* consistency and found
   one provenance error (§12.5's sample size), but it did not re-run a single probe. Twenty-one of
   the twenty-two probes stand on their original run and on §26's closing statement that the suite
   reproduced them exactly.

**One finding about the review itself, and it is the reason item 14 exists.** The load-bearing
decision — the estimand key set — was first argued from *engineering* grounds (a uniform 23-keys-per-
arm shape makes the count derivable) and settled the other way only when [§14a]'s Report sentence was
read directly. The correction came from the governing document and not from the review. **That failure
mode has not been swept for**, and §21 item 14 states it as a gap rather than leaving the one instance
to stand as if it were the only one.

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
trigger fires (item 4), **two** have their triggers fire and stay open — the parallelisation item as
the largest cost in the pipeline (item 5) and the `model.FitError` `code`-field item, whose *"first
bucket that is wrong"* this document supplied (item 12) — and **nine** are new (items 6–14), of which item 14 is the review's own residual.

### 27.1 What the engineering review changed, 2026-08-27

Recorded here rather than folded in silently, because five of the seven were **internal
contradictions in this document** and a reader is entitled to know the drafting failed in those
places. Every one is discharged in a numbered section above:

| Finding | Was | Is |
|---|---|---|
| The estimand key set | 70 keys, 23 per arm, including `support.beta` | **69**, 22 standardised per arm + one `β` **per fit** + `hier.sigma` (§3.3). [§14a] reports `exp(β)` as a model parameter and defines the sensitivity on the *standardisation population*, so the support arm has no `β` |
| A pooled `FitError` | cost every key **except** `hier.*` | costs **every** key (§13.2). §12.8 says `polr_ri` starts from the pooled fit, so the exemption was unimplementable |
| The `POLR_MAX_ABS_BETA` guard | named "G-bound", no identifier, no raise site, no bucket token | **T12**, one helper, **two** call sites — `beta` and `hier.beta` (§9.1, §14, §18) |
| T4's bucket | led with `polr_ri:`, so it would have counted as `nonconvergence` | leads with `T4` (§14). The failure the raise-site scan cannot catch |
| `Hierarchical.standardisation.fit` | typed `model.PolrFit`, held a `RIFit` | `PolrFit \| RIFit`, third `population` value, `is`-identity with `Hierarchical.fit` (§3.1, §12.6) |
| The hierarchical arm's rates | quoted as if `N_BOOT` | labelled **n = 200**, with §25's step 0 to re-measure (§12.5, §26) |
| `balance._levels` | reimplemented in §10.3 because it is private | **public as `balance.levels`** (§10.3, §18, §24), on Stage 7 §5.5's precedent for `smd` |

**Six smaller corrections, none of which would have stopped an implementer but each of which would
have made one guess:**

| Was | Is |
|---|---|
| §17 and §19's drivers used `over=sup.inside` **before** `support` ran, and §19 omitted sensitivity 1's point estimate entirely — nine audit entries against §20.9's ten | §19 is the **one** canonical call order, with its six calls and the entries each records; §15 and §20.9 reference it rather than each implying an order |
| T10 and T11 had **no tolerance and no stated check site**, in a document whose header forbids inventing what is not written | §14 states both, and §20.6 explains why they differ: `mrs_0_2` **is** `rd[2]` so the comparison is exact, `mortality` is a second computation so it is 1e-15 |
| §12.2's prose quoted three error values — adaptive at σ=3 / 3 nodes, non-adaptive at σ=2 / 5 and 7 — that were **not in the table it cited** | The three rows are in the table. A number quoted from a table it is not in is a number nobody can check |
| §16's `population` was the only entry point whose docstring did not name its audit entries | "Records entries 1 and 2" |
| "Four records" meant `{Standardisation, Support, Hierarchical, RIFit}` in §3.1 and `{Standardisation, Support, Hierarchical, Bootstrap}` in §1 and §0 | §3.1 says plainly: **three** of this stage's own, plus `RIFit` which is `model.py`'s and `Bootstrap` which `inference` returns |
| `dist1_*` / `dist0_*` looked like the arm codes §7.1 forbids writing as literals | §3.3 states they are fixed **key names** — the wire format Stage 14 matches on — deliberately not derived from `C.TREATMENT_LABELS`, so a registry edit cannot silently rename a reported key |

**And one thing the review added that is about the environment rather than the document.** §25's step 0
now records that `DATA_XLSX` is gitignored, so in a worktree or a fresh clone every workbook-dependent
test skips — measured **1343 passed / 112 skipped** in this branch's worktree against the main
checkout's 1588 / 2. Nothing is wrong when that happens, but a green suite in a worktree is not
evidence about any [§14a] number, and step 0's re-measurement must run where the data is.

**§11's `over=` gained a named dependency it had been relying on silently.** It is applied inside every
replicate to a *resampled* frame, and is safe to index by label only because `bootstrap.resample`
resets the index first — `bootstrap.py:315-320`, which gives Stages 12 and 13 as its reason in as many
words. §20.14 names it, so the property is visible from the caller that needs it and not only from the
function that provides it.

---

## 28. The implementation record, 2026-08-27

This document was implemented on the date it was written, against the landed Stages 1–10 and the v7
workbook (`DATA_SHA256` unchanged), Python 3.12.12 / numpy 1.26.4 / pandas 2.3.3 / statsmodels 0.14.6
— §2's environment exactly, and **no dependency was added**.

**Read this section as §26's continuation rather than as its replacement.** §26 records what the
*probes* measured before any code existed; this records what the *shipped code* measures. Where the two
disagree about a number, this one is the measurement — and the earlier figure is left in place, labelled
as the projection it was, so that what was projected and what was measured stay distinguishable.

### 28.1 What reproduced, and it is nearly everything

**Every figure in §26 that the implementation could reach reproduced, on code written from this
document rather than from the probes.** That is the strongest available evidence that the document is
implementable as written, and it is listed rather than summarised because a claim of agreement nobody
can check is worth nothing:

| §26 probe | figure | reproduced |
|---|---|---|
| A1 | 126 → 107 → 104; 1 covariate-incomplete on `core_ml`+`tmax6_ml`; 2 outcome-missing, both USZ | exact |
| A1b | `treating_centres`' complement is `('USZ',)` | exact |
| A2 | all 7 mRS levels occupied at the point estimate | exact |
| A3 | 10 design columns, 0 dropped, 6 cutpoints | exact |
| A4 | 5 iterations on `"likelihood"`, 0 rescales, 0 halvings, max abs coefficient 0.73 | exact (0.7311) |
| A5 | the box per covariate; 93 inside / 11 outside, all comparator; 3 of 12 outside | exact |
| A5b | the two "continuous" readings give the same outside-set | exact |
| A6 | 39 / 12 / 53; `center = USZ` is `nan`; largest non-grouping SMD above 0.8 | exact (1.0197 on `hypertension`, 0.843 on `prestroke_mrs`) |
| A7 | 12 audit entries before this stage; step names disjoint from `cohort.build`'s | exact |
| B0b | `|mrs_0_2 − RD_2|` 0.0; `|mortality + RD_5|` 8.3e-17 | exact (0.0 and 8.327e-17) |
| B2 | the cutpoint collapse, 83 / 2000 = 4.15%, always mRS 5 | exact |
| B3 | `polr` iterations 4–6 across 2000 replicates | exact (4: 25, 5: 1958, 6: 17) |
| B4 | `|beta|` min 0.0000 median 0.3219 max 2.0016; 0 reaching 14.0 | exact |
| B5 | design columns dropped: 0 in 100.00% of replicates | exact |
| B7 | the reflection check, 9.8e-5 / 1.7e-3 against `"linear"`'s 1.1e-16 / 5.6e-16 | exact (9.834e-05 / 1.705e-03) |
| B9 | every distribution sums to 1 and is monotone in every arm of every replicate | exact |
| D1–D4 | the analytic gradient against central differences | 4.3e-8 to 6.9e-8, against the probe's 1e-8 — same order, both far below `POLR_SCORE_TOL` |
| E1–E4 | adaptive stable from 9 nodes; non-adaptive NON-MONOTONE in the node count | reproduced on an independent frame: adaptive spread 3.5e-11 over 9/11/15/31; non-adaptive 0.632 / 0.737 / 0.530 / 0.567 / 0.599 at 5 / 7 / 9 / 11 / 15 |
| §12.3 | σ̂ = 0.542426850 at 11 adaptive nodes | exact to the digits quoted (0.542427) |
| §12.4 | 25 iterations, 85 halvings, 2 rescales at the point estimate | exact |
| §12.6 | the four posterior SDs: 0.2517 / 0.3156 / 0.2823 / 0.3810 | exact to four decimals, and USZ is still the least precise |
| §20.11 | `bootstrap.run`'s twenty-six intervals byte-identical across the `intervals` extraction | **exact — all 26, every limit, level, method, draw count and `p`** |

**§6.3's predicted consequence was also confirmed rather than assumed**: `RD_4 == RD_5` holds
*exactly* in 83 of 2000 draws, which is the same 83 replicates that lost mRS 5.

### 28.2 §25's step 0 — the hierarchical arm at `N_BOOT`

Run before the acceptance suite was written, as §25 requires, and it changed four figures:

| quantity | §26, n = 200 | measured, `N_BOOT` = 2000 |
|---|---|---|
| σ̂ at the `POLR_RI_SIGMA_FLOOR` boundary | 44 / 200 = **22.0%** | 368 / 1998 = **18.4%** |
| σ̂ above 1.0 | 48 / 200 = **24.0%** | 460 / 1998 = **23.0%** |
| `polr_ri` non-convergence | **0 / 200** | **2 / 2000** |
| the whole three-arm bootstrap, wall clock | **17.5 min**, projected ×10 | **4.2 min**, measured |
| σ̂ range | min 0.000021, median 0.65923, max 1.78734 | min at the floor, median 0.69084, max 2.42713 |

**Three of the four moves are unremarkable and the fourth is not.** The two rates are within sampling
error of the 200-replicate estimates — 18.4% against 22.0% is 1.2 standard errors at n = 200 — so §26
was right to call them imprecise rather than wrong. The wall clock is a factor of four out, and a
timing is machine state, which §26 already says of every timing it quotes; the arm is still ~95% of the
stage's compute and the stage is still the most expensive thing in the pipeline, so §21 item 5's
parallelisation trigger stays fired.

**The non-convergence count is the one that matters, because it is a behaviour and not a precision.**
§13.4's *"Nothing fails"* is true of the pooled and support arms at 2000 replicates and false of the
hierarchical arm: two replicates exhaust `POLR_MAX_HALVINGS` in a genuine BFGS line-search stall, at
iteration 22 and at iteration 7. That is 0.1%, each costing `hier.*` and nothing else — 24 of the 69
keys — so the pooled and support arms keep 2000 draws while the hierarchical arm has 1998, and every
`Draws` reconciles. **The taxonomy handles it exactly as §13.2 designed it to.** What is *not* settled
is whether a 0.1% drop selected by the optimiser rather than by the data is exchangeable with the
survivors, which is Stage 6 §5.3's concern in a new place; `TODOS` carries it with a 1% trigger.

**A recovery path was considered and declined.** Resetting the BFGS inverse-Hessian to the identity on
a line-search failure and retrying as a gradient step is a standard safeguard and would probably
recover both replicates. It is declined because §12.4 fixes the loop's shape and Stage 6 §5.3's
argument applies with full force: [§10] drops failed replicates, so a change to *which* replicates
fail is a change to the sampling distribution, and adding a recovery path is as much a modelling
decision as removing one. That is a [§14] amendment, not a repair.

### 28.3 Three things this document did not anticipate, all found by running the code

**None was reachable by reading. Two are numerical and were invisible at the point estimate, at
n = 200, and in every structural test; the third is a name collision this document's own rename
created.**

1. **`p = gu − gl` underflows to EXACTLY ZERO at far quadrature nodes.** Both are `expit`s; at a node
   the conditional posterior has moved away from, both are 1.0 in float64 and the difference cancels.
   Then `log p` is `−inf` and `du/p`, `(du/p)²` and `du·dl/p²` are all `nan` — the gradient is `nan`,
   the trial step is rejected by the `isfinite` test in the halving loop, and **the fit converges to
   the right answer through a line search doing the wrong job.** `model._RI_P_FLOOR = 1e-150` floors
   it. The exponent is chosen against `p²` and not against `p`: `Hul = du·dl/p²`, so a floor below
   ≈1.5e-154 makes `p²` underflow and reintroduces the 0/0 it fixes.

   `polr`'s convention is the opposite — return `−inf` and let step-halving reject the step — and it
   is right *there*, where a non-positive probability means crossed cutpoints at the iterate. Here it
   means a far node, which is not a property of the iterate at all and is what a quadrature rule
   exists to weight at zero.

2. **The per-group log-posterior curvature comes out POSITIVE from cancellation.** It is negative
   exactly in the mathematics — a sum of concave ordinal terms plus a Gaussian log-prior, so it cannot
   exceed `−1/σ²` — but `Huu = ddu/p − (du/p)²` is a difference of two large quantities in the tails.
   When the sum crossed zero, `sqrt(−1/curvature)` returned `nan`, the nodes were `nan`, and the line
   search absorbed it again — and in **2 of 2000 replicates it exhausted the halvings and raised T2 for
   a floating-point artefact.** `model._curvature` clamps the ordinal contribution to its own
   guaranteed sign, `min(g2, 0)`. It enforces a property rather than choosing a magnitude, and the
   bound is tight: attained whenever the data contribute no curvature, which is the far-tail case
   producing the noise.

   **This is §12.5's own argument turned on this document.** §12.5 says of the σ → 0 boundary that
   *"the fit still returned … but 'still returned' is not a specification"* — and then the specified
   fitter had two further instances of exactly that failure mode. `test_model.py` now asserts `polr_ri`
   runs clean under `warnings.simplefilter("error", RuntimeWarning)` on both a non-zero-σ frame and a
   floored one, which is the check that would have caught them.

3. **§13.2's own rename created an `UnboundLocalError` that no structural check could see.** Making
   `bootstrap._bucket` public as `bucket` collided with a local variable named `bucket` inside
   `_replicate`, which holds the classified value. The module-level function was shadowed and the
   first replicate that classified a failure raised `UnboundLocalError`. Imports were clean, every
   structural test passed, and **the full 2000-replicate pipeline is what found it.** The three locals
   are now named `label`.

### 28.4 Where the implementation departs from this document, and why

**Five departures. Four are forced by Stage 11's absence and one is a measured tolerance.** Nothing
else in §§1–27 was contradicted.

| # | This document | The implementation | Why |
|---|---|---|---|
| 1 | `bootstrap.bucket` and `bootstrap.collect` are public (header, §13.2) | they are made public **here**, as a pure rename | Stage 11 has not landed at all, so they were still private. §21 item 11's trigger in its stronger form; the rename is the same shape §18 already authorises for `balance._levels` |
| 2 | *"Stage 11's two callers move onto it"* (§13.3) | only `run` moves, because Stage 11's callers do not exist | same cause. `bootstrap.intervals` takes the p-rule as a parameter exactly as specified, so Stage 11's callers will read it when they land |
| 3 | public surfaces go 8 → 9, 6 → 8 and 3 → 4 (§3.2, §18) | 5 → 8, 5 → 7 and 2 → 3 | **every DELTA is right and every ABSOLUTE is one high** — the signature of counting against an unlanded baseline. The tests assert by NAME rather than by count |
| 4 | §20.12: at the floor, `polr_ri`'s `alpha` and `beta` agree with `model.polr`'s **to 1e-6** | asserted at **1e-5** | measured 4.9e-6 on the cutpoints and 1.3e-6 on the coefficients. At the floor σ is 1e-4 and not 0, so the marginal model is not exactly the pooled one, and both fits stop on `POLR_TOL` rather than at an exact stationary point. 1e-6 is a bound this construction does not reach |
| 5 | §10.3 leaves the baseline table's `role` column unspecified beyond `"grouping"` | `standardise._role` is a Stage-12 vocabulary and does **not** call `balance._role` | `balance._role` returns `"propensity model"` for every `PS_COVARIATES` name, and **[§14] fits no propensity model anywhere** — that label would name a model this analysis does not contain. Its docstring defines a role as *"what this covariate is TO the specification being diagnosed"*, and the specification here is a different one |

**And one thing this document specified that the implementation could not fully honour**: §10.1's box
readings, §20.10's counts and §20.1's ledger are all workbook-gated, so on a checkout without `data/`
they skip — §25's step 0 makes the same point from the other side, and it is why every one of them was
run against the workbook here rather than against a fixture.

**The suite, measured with the workbook present: 1551 passed, 6 skipped, 0 failed.** The six skips are
§20.11's `STAGE12_SLOW` byte-identity guard (run separately, and it passed), two Stage 10
`STAGE10_SLOW` coverage designs, and three `firthlogist` oracles whose optional dependency is not
installed. **No R oracle skips on this machine**, which is the difference from §26's own run.

### 28.5 What this implementation did NOT close

**§21 item 14 is untouched and remains the largest open risk in this document.** It asks for a sweep
for structural decisions argued from engineering grounds where [§10], [§14] or [§16] gives a direct
answer. **The implementation cannot perform that sweep and did not attempt it**: it implements what
this document says, so a place where this document reasons around the plan rather than from it becomes
code that matches the spec and disagrees with the prespecified analysis — which is precisely the
failure mode item 14 describes, and precisely the one an implementer is worst placed to notice.

**§26's LARGEST VERIFICATION GAP IS CLOSED, AND THAT IS THE ONE THING HERE THAT REVERSES A CLAIM THIS
DOCUMENT MAKES.** §26 states that *"this is the one estimator in this pipeline whose specification
ships without an independent implementation having agreed with it"*, because R segfaulted at startup on
the machine the spec was written on. **R works on the machine it was implemented on, so the
`ordinal::clmm` oracle RAN.** Measured against `ordinal` 2026.7.26 at `nAGQ = POLR_RI_NODES = 11`, on
an eight-group frame:

| parameter block | `polr_ri` against `clmm` | and the asymmetry |
|---|---|---|
| coefficients | `|ours + clmm|` = **1.033e-05** | NEGATED, and `|ours − clmm|` = 1.574, so there is something to negate |
| thresholds | `|ours − clmm|` = **3.636e-06** | NOT negated |
| σ | `|ours − clmm|` = **4.446e-06** | NOT negated either, because it is a scale |

**All three parts of §12.7's three-way asymmetry hold, and they agree to about 1e-5 on every block.**
That is an independent maximiser of the same likelihood reaching the same answer, which is the one
comparison §12.7 says can fail for the reason the estimator would be wrong. `polr_ri` is therefore
verified against `clmm`, against its own analytic gradient, against node-count stability, against its
collapse to `model.polr` at the floor, against a known-σ recovery, and against the posterior-SD
ordering.

**One defect in the oracle itself is worth recording, because `clm` and `clmm` differ where they look
alike.** The first version of `polr_ri_clmm.R` read the random-effect SD from `fitted$stDev`, which is
what the name suggests and which is **NULL on a `clmm` object**: the SD is in `fitted$ST`, an
`lme4`-style list whose single 1×1 matrix entry is the standard deviation for a lone random intercept.
`as.numeric(NULL)` is `numeric(0)`, so `data.frame()` failed on a length mismatch and **the oracle
failed for a reason that had nothing to do with the estimator** — which is exactly the failure mode a
gated oracle is worst at surfacing, because a red test and a skipped test both mean "not verified".

**No cross-model review of the implementation was run**, for the same reason §26.1 records for the
spec. This code has been written once and reviewed by nobody — and §21 item 14, which asks for a sweep
this implementation cannot perform, is untouched.
