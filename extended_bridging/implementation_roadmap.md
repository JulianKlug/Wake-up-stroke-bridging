# Implementation roadmap

Build order for `statistical_analysis_plan.md`. Each stage lists what to build, the SAP section it
implements, and the acceptance criteria that must pass before the next stage starts. Section
references in brackets are to the SAP.

Stages 1–6 are prerequisites for everything. Stages 7–11 produce the primary results. Stages 12–13 are
independent of the propensity machinery and can be built in parallel with 7–11.

Some code has already build inside pilots. This should not be a ground source, but can be used for inspiration. Some implementations may be wrong. 
---

## Stage 0 — Facts to establish before coding

These instantiate the plan; the code cannot be written against guesses.

- [x] Which centres contributed zero bridging patients [§3]. Determines the primary cohort.
- [x] Whether a contraindication-reason field exists, and the full set of distinct values it takes
      [§3, §11]. Determines whether three-category eligibility is possible at all.
- [x] Per-centre completeness of that field. Any centre where it was never collected produces the
      *indeterminate* group.
- [x] Event counts for each binary outcome, to know which fall under the <10-event rule [§8].
- [x] Treated-arm size, to confirm the ~3.5-per-parameter degrees-of-freedom budget in [§6].

Record all six in a short data-inventory note. None requires looking at outcomes by arm.

**Done.** `stage0_data_inventory.py` establishes all five from the workbook and writes
`../out/stage0_data_inventory.md` (gitignored: every number in it is patient-derived). Re-running it
reproduces the note byte-identically.

**Decisions taken (PI, 2026-08-06), before any outcome was examined by arm.** Stage 1 encodes these;
the note carries the consequences in full.

1. **`IVT_contraindicated_binary` is the [§4] classifier. The free text of `Contraindications_to_IVT`
   is not read.** The reason column still contributes one bit — whether a reason was recorded —
   because the flag is 0 both for patients documented as having no absolute contraindication and for
   patients at a centre that never collected the field, and [§3] forbids reading that blank as "no
   contraindication". So: treated → eligible; flag = 1 → ineligible; flag = 0 with a reason recorded
   → eligible; flag = 0 with none recorded → indeterminate.

   Two knock-on effects. Because no string is matched, the [§4] requirement to assert on unrecognised
   values has nothing to range over; what Stage 4 must assert instead is that the flag is 0/1 and
   never missing, and that **no treated patient carries a 1** — the case the revealed-fact rule would
   otherwise silently absorb. And the free-text normalisation problem disappears: the `Clnician` typo
   and the parenthetical annotations can no longer break anything.

2. **The 90-day mRS is ground truth for vital status; death is `mRS == 6`.** `Deathat90days` and
   `mRS56at90days` are dropped in the Stage 1 mapping rather than left available. Death in the primary
   cohort is 25 events, not the shipped 26. The two records whose shipped death flag contradicts an
   mRS of 1 and 4 are read as alive — a derivation rule, not a resolution; they remain a query for the
   data owner and are listed by identifier in the note.

3. **The [§8] rare-outcome threshold applies to `min(events, non-events)`, and `TICI 2b–3` is
   augmented with a reduced `m_a(X)`** — treatment + `center` + `atrial_fib`, 5 parameters against 6
   non-events — rather than dropped from augmentation. Recorded as a dated amendment in [§8], since
   it is a declared exception to both nuisance models sharing one covariate list. Covariates chosen
   on procedural grounds, never from the data. The propensity model is untouched: its budget comes
   from the treated arm, not from any outcome. `sICH` and `PH2` remain unaugmented.

**Nothing open.** Stages 1 onward can be built against the note.

---

## Stage 1 — Configuration and data contract

**Spec:** `specs/stage1_config_and_data_contract.md`.

**Build.** A single configuration module holding: file paths; raw→analysis column mapping; the [§6]
covariate list; the outcome registry (name, family, type per [§5]); the eligibility reason
classification [§11]; bootstrap replicate count and random seed [§10]; the |SMD| threshold [§9]; the
rare-outcome event threshold [§8]. Nothing outside this module may reference a raw column name.
Setup an UV env. 

**Accept when.** Every raw column is either mapped, explicitly dropped with a reason, or fails an
assertion. No silent passthrough.

## Stage 2 — Load and clean [§11]

**Build.** Reader that treats `'N/A'` and empty strings as missing. All corrections content-driven —
matched on values, never on row index — and each appended to an audit log carrying the affected case
identifiers. Schema assertions that fail loudly on drift: row count, unique identifiers, treatment
non-missing, plausible ranges for age, NIHSS and mRS, non-negative volumes and times.

**Accept when.** The audit log reproduces byte-identically on a re-run, and every correction in it
names the cases it touched.

## Stage 3 — Derived variables [§5, §6, §13]

**Build.**
- `onset_type`, a three-level factor from the wake-up and unwitnessed flags. Assert the two are never
  both positive.
- All outcome dichotomies derived from their ordinal source [§5]. **Reimpose missingness explicitly**
  — a comparison against a missing value returns false in most frameworks and would silently
  manufacture zeros.
- Subgroup variables [§13]: unknown vs witnessed onset, target mismatch, core volume above/below
  median.
- Drop any covariate with zero variance.
- Mark structurally non-applicable fields as distinct from missing ones [§11].

**Accept when.** A test asserts that every derived dichotomy has exactly the missingness of its
ordinal source, and that the three onset levels partition the cohort.

## Stage 4 — Eligibility classification [§3, §11]

**Build.** A function mapping each patient to `ineligible` / `eligible` / `indeterminate`, under
DECISION 1: treated → eligible by revealed fact; `ivt_contraindicated` = 1 → ineligible; flag = 0
with a reason recorded → eligible; flag = 0 with none recorded → indeterminate. The free text of
`Contraindications_to_IVT` contributes one bit — whether a reason was recorded — and is never read.

**Assert the flag, not the text** [DECISION 1]. No string is matched, so there is no set of
recognised reasons to check a new value against, and the earlier "every observed reason string is
classified" requirement has nothing to range over. Assert instead that `ivt_contraindicated` is 0/1
and never missing, and that **no treated patient carries a 1** — the case the revealed-fact rule
would otherwise silently absorb. This also retires the free-text normalisation problem: the Stage 0
note records a `Clnician` typo and parenthetical annotations that would have made an exact-string
classifier raise on four of five spellings of one reason.

**Accept when.** A cross-tabulation of eligibility by centre and arm is produced, and a test confirms
no blank is ever read as "no contraindication".

## Stage 5 — Cohort construction [§2, §3]

**Build.** Two restrictions applied in order, each recording what it removed:
1. Drop centres with zero treated patients.
2. Drop ineligible patients.

Patients of indeterminate eligibility are retained [§3]. Emit a cohort-flow table that reports their
count separately, so the size of the retained-but-undocumented group is visible.

**Accept when.** Every centre in the resulting cohort contains both arms; no ineligible patient
remains; the cohort-flow table shows the indeterminate count as its own line.

## Stage 6 — Propensity score and weights [§7]

**Build.**
- A design-matrix builder: reference-coded dummies for factors, constant columns dropped (a resampled
  or subset cohort can leave a factor level empty and make the design singular).
- **Firth-penalised logistic regression.** Maximise `l(b) + ½ log|I(b)|`; modified score
  `X'(y − p + h(0.5 − p))` with `h` the hat-matrix diagonal. Step-halving on the penalised likelihood;
  converge on the penalised likelihood change or a flat modified score. **Raise on failure — never
  fall back to a different estimator** [§7].
- Overlap weights `w = 1 − e` treated, `w = e` control.
- Kish effective sample size per arm.

**Accept when.** Tests confirm: at large n on well-behaved data the coefficients match an unpenalised
maximum-likelihood fit; under complete separation the fit stays finite where the unpenalised one
diverges; the returned probabilities are finite and in (0, 1).

**Record in the output**, per [§7]: ESS per arm, a description of the weighted population, and the
statement that both are conditional on this propensity specification because the ATO target population
is `h(X) = e(X){1 − e(X)}`.

## Stage 7 — Balance and overlap diagnostics [§9]

**Build.**
- Standardised mean differences before and after weighting, using the **unweighted** pooled SD as a
  common denominator so the yardstick does not move.
- Judged against the full [§6] confounder set regardless of which covariates a given specification put
  in its own propensity model.
- The four vascular risk factors reported as negative controls [§6].
- A **within-centre** overlap table alongside the pooled one: per centre, numbers per arm, propensity
  range per arm, ESS, maximum weight, weight share, worst within-centre |SMD|. Centres with structural
  non-positivity are reported as such rather than given an overlap plot.

**Accept when.** The SMD function reproduces a hand-computed value, and returns missing rather than
zero for a variable observed in only one arm.

## Stage 8 — Primary outcome estimator [§8]

**Build.**
- A **weighted proportional-odds model with treatment as the sole predictor**. Most implementations
  have no observation-weight support; weighting the log-likelihood contribution per observation is
  sufficient. Return `β` and `exp(β)`, oriented so > 1 favours bridging.
- **Weighted empirical cumulative risk differences** `RD_k`, `k = 0…5`, taken from the weighted
  distributions rather than from threshold-specific models, so the cumulative probabilities stay
  ordered by construction.

**Accept when.** The weighted fit equals an unweighted fit when all weights are 1; the cumulative
probabilities are monotone in `k` within each arm; the orientation test confirms `exp(β) > 1` when
treatment shifts mRS downward.

## Stage 9 — Secondary binary estimators [§8]

**Build.**
- Weighted risk difference: difference of weighted means.
- Weighted marginal odds ratio, kept finite when a weighted proportion reaches 0 or 1.
- An outcome regression `m_a(X)`: treatment main effect plus **the [§6] covariate set**, linear terms,
  Firth logistic. Take the covariate list from the same configuration entry the propensity model uses,
  so the two cannot drift apart.
- **One declared per-outcome override** of that list, per the [§8] amendment: `tici_2b_3` →
  treatment + `center` + `atrial_fib`. Build it as a registry keyed by outcome that *defaults* to the
  shared [§6] entry, so an outcome can only diverge by being named — never by an edit to the default
  quietly applying to one outcome. The override must reach the reporting layer and be printed beside
  the TICI estimate.
- The **model-assisted augmented** estimator:

      tau = Σ₁ w(Y − m₁)/Σ₁ w  −  Σ₀ w(Y − m₀)/Σ₀ w  +  Σ h(m₁ − m₀)/Σ h

  with `h = e(1 − e)`, the same tilting function as the weights.
- A guard: **outcomes whose minority cell — `min(events, non-events)` — is below 10 are not
  augmented** [§8, as amended]. Keying this on the event count alone would let TICI 2b–3 through with
  6 non-events.

**Accept when.** Two tests pass:
1. With the outcome model set to a constant, the augmented estimate equals the unaugmented weighted
   risk difference exactly. If it does not, the augmentation is not built on the weights' tilting
   function and is targeting a different estimand.
2. A heterogeneous-effect scenario with a **correct** outcome model and a **misspecified** propensity
   model does *not* recover the true ATO. This must be asserted, not merely allowed: the estimator is
   not doubly robust [§8], and a test built on a homogeneous effect will appear to show that it is,
   because every weighted average treatment effect coincides in that case.

**Naming.** "model-assisted", never "doubly robust", in code, docstrings and output.

## Stage 10 — Bootstrap engine [§10]

**Build.** Patient-level resampling stratified by centre, 2000 replicates, recorded seed. In every
replicate refit the propensity model *and* the outcome regression, recompute weights, recompute every
estimate. Percentile 95% intervals. **Replicates whose prespecified fit fails are dropped and counted;
never substituted with a different estimator.** Surface the failure count and rate in the run log.

**p-values.**
- Primary outcome: from the bootstrap distribution of `β`,
  `p = 2·min{Pr(β̂* ≤ 0), Pr(β̂* ≥ 0)}`, floored at `1/(B+1)`. Label it a test of the
  proportional-odds treatment coefficient — **not** a risk-difference-scale test.
- Secondary binary and safety outcomes: on the risk-difference scale.
- Cumulative `RD_k`: intervals only, **no p-values** [§8].

**Accept when.** On synthetic data with a known effect the interval covers the truth at roughly the
nominal rate; the failure counter is exercised by a deliberately degenerate input; and no code path
can return an estimate from a fallback estimator.

## Stage 11 — Multiplicity, subgroups, E-value [§6, §13]

**Build.**
- Benjamini–Hochberg within the secondary and safety families; primary uncorrected; raw and adjusted
  p both reported [§13].
- Subgroup estimates with a treatment × subgroup interaction test, labelled hypothesis-generating
  [§13].
- **E-value for the primary estimate** [§6]. The primary effect measure is a common odds ratio, so
  convert with the usual approximation for a common outcome (`RR ≈ √OR`) and state that the
  approximation was used. Report the E-value for the point estimate and for the confidence limit
  nearest the null, with the rule that **the limit's E-value is 1.0 whenever the interval already
  spans the null** — feeding a null-crossing limit into the formula measures how far past the null the
  interval reaches and returns a large number that reads as robustness while meaning the opposite.

**Deferred.** The sensitivity-analysis suite in [§13] is out of scope for now by decision. Leave the
hook but do not build the specifications.

## Stage 12 — All-centre standardisation [§14a]

Independent of Stages 6–10; **no propensity model is fitted anywhere in this stage.**

**Build.**
- One proportional-odds model over all eligible patients at all centres:
  `logit P(Y ≤ k | A, X) = α_k + βA + γᵀX`, with `X` = the [§6] covariates **minus centre**.
- Standardisation: duplicate every patient with `A = 1` and `A = 0`, predict the full category
  distribution `P(Y = j | A = a, X)` for `j = 0…6`, average over all N eligible patients.
- Outputs: the two standardised mRS distributions; cumulative `RD_k`; standardised mRS 0–2 risk
  difference; standardised mortality difference.
- **Support check:** treated-arm distribution of each continuous covariate, the proportion of
  patients from never-IVT centres falling outside it, and a baseline table comparing IVT-treated
  patients with never-IVT-centre patients.
- Sensitivity: restrict the standardisation population to patients inside the treated support.
- Sensitivity: random centre intercept, one common treatment effect, no interaction and no random
  slope.
- Inference: patient-level bootstrap stratified by centre — resample, refit, predict everyone under
  both regimes, average, recompute, percentile intervals.

**Accept when.** The standardised category probabilities sum to 1 within each regime; the cumulative
probabilities are monotone; the estimate is labelled an **ATE in the all-centre eligible population**
and carries the statement that it is not the [§7] ATO and the two are not a like-for-like comparison.

**Guard.** `exp(β)` here is a *conditional* odds ratio. It may be emitted as a model parameter and must
never be labelled as the standardised marginal effect.

## Stage 13 — Feasible-policy analysis [§14b]

**Build.** The same standardisation machinery over all patients at all centres, but the active regime
assigns treatment **as a function of eligibility** — bridge if eligible, direct thrombectomy if
contraindicated — against a comparator of direct thrombectomy for everyone.

**Accept when.** No contraindicated patient is ever assigned a predicted IVT outcome, and the output is
labelled an operational policy contrast, distinct from both [§14a] and [§7].

## Stage 14 — Outputs and guardrails [§16]

**Build.** A single entry point writing every table, figure and log, plus a run summary recording the
seed, replicate count, and bootstrap failure counts.

**Guardrails to enforce in the reporting layer**, so they cannot be lost in editing:
- Any output placing a [§14] estimate beside the primary carries the different-populations statement.
- The augmented estimator is labelled model-assisted everywhere.
- Safety outcomes are labelled descriptive.
- Every estimate carries its own denominator [§11].
- The manuscript checklist from [§16] is emitted as a file, not left implicit.

---

## Invariants worth asserting in the test suite

Independent of any single stage, these should fail loudly if a future edit breaks them:

1. No never-IVT centre appears in the primary cohort.
2. No ineligible patient appears in any [§7] or [§14a] population.
3. The propensity and outcome-regression covariate lists resolve to the same configuration entry, for
   every outcome except those holding a declared override — currently `tici_2b_3` alone [§8
   amendment]. Assert the override registry against that explicit list, so a second outcome cannot
   acquire a reduced model without the assertion failing.
4. No post-time-zero variable appears in any model's design matrix — assert against an explicit
   denylist, not by inspection.
5. No estimate is ever produced by a fallback estimator.
6. Derived dichotomies carry exactly the missingness of their ordinal source.
7. The augmented estimator collapses to the unaugmented one under a constant outcome model.
