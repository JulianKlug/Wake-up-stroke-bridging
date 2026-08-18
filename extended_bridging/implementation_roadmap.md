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
      [§3, §11]. Determined that three-category eligibility was possible; DECISION 1a has since
      withdrawn it, and the field is now read by nothing.
- [x] Per-centre completeness of that field. Under DECISION 1 any centre where it was never collected
      produced the *indeterminate* group; under DECISION 1a it identifies the 43 retained controls
      whose eligibility rests on the flag alone, which Stage 5's cohort-flow table counts.
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

   **Amended by DECISION 1a (PI, 2026-08-13): the reason column contributes no bit either, and
   eligibility is two-valued.** Flag = 1 → ineligible; flag = 0 → eligible, including every patient
   whose reason was never documented. `Contraindications_to_IVT` is read by nothing. [§3]'s amendment
   of the same date is the protocol record and `../out/stage0_data_inventory.md` carries DECISION 1a
   in full. It **does not change the population** — the 43 `indeterminate` records were already
   retained, so the primary cohort is the same 93 patients — and it makes the no-treated-patient
   assertion load-bearing rather than belt-and-braces, because the revealed-fact rule it guarded is
   deleted. The substantive calls in the two bullets above stand: they follow from taking the flag as
   the classifier, which DECISION 1a does more strictly than DECISION 1 did.

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

**Spec:** `specs/stage2_load_and_clean.md`.

**Build.** Reader that treats `'N/A'` and empty strings as missing. All corrections content-driven —
matched on values, never on row index — and each appended to an audit log carrying the affected case
identifiers. Schema assertions that fail loudly on drift: row count, unique identifiers, treatment
non-missing, plausible ranges for age, NIHSS and mRS, non-negative volumes and times.

**Accept when.** The audit log reproduces byte-identically on a re-run, and every correction in it
names the cases it touched.

## Stage 3 — Derived variables [§5, §6, §13]

**Spec:** `specs/stage3_derived_variables.md`.

**Build.**
- `onset_type`, a three-level factor from the wake-up and unwitnessed flags. Assert the two are never
  both positive.
- All outcome dichotomies derived from their ordinal source [§5]. **Reimpose missingness explicitly**
  — a comparison against a missing value returns false in most frameworks and would silently
  manufacture zeros.
- Subgroup variables [§13]: unknown vs witnessed onset, and core volume above/below median. The
  target-mismatch subgroup is withdrawn by the [§13] amendment of 2026-08-10; the cohort is
  CTP-selected, so it is close to the criterion that admitted these patients.
- The median is a property of the cohort, so it is computed after the Stage 5 restrictions and then
  **frozen**: subgroup membership is a fixed patient attribute the Stage 10 bootstrap resamples along
  with the patient, never recomputed inside a replicate.
- **Detect and log** any covariate with zero variance. Nothing is dropped here and no covariate list
  is mutated — Stage 6's design-matrix builder drops constant columns at the point of use, which is
  also where a bootstrap replicate can empty a factor level.
- Mark structurally non-applicable fields as distinct from missing ones [§11], extending the Stage 2
  classification to the derived columns.

**Accept when.** A test asserts that every derived dichotomy has exactly the missingness of its
ordinal source, and that the three onset levels partition the cohort.

## Stage 4 — Eligibility classification [§3, §11]

**Spec:** `specs/stage4_eligibility_classification.md`.

**Build.** A function mapping each patient to `eligible` / `ineligible`, under **DECISION 1a**:
`ivt_contraindicated` = 1 → ineligible; flag = 0 → eligible, including every patient whose
contraindication reason was never documented. `Contraindications_to_IVT` is read by nothing — neither
its text nor its presence.

*Built and landed under DECISION 1, which made eligibility three-valued; the module and its tests were
amended on 2026-08-13. The spec's §21 records what moved.*

**Assert the flag, not the text** [DECISION 1a]. No string is matched and no presence is tested, so
there is no set of recognised reasons to check a new value against. Assert instead that
`ivt_contraindicated` is 0/1 and never missing, and that **no treated patient carries a 1**. That last
assertion is **load-bearing** under DECISION 1a: it is what makes the deleted revealed-fact rule
redundant, so without it a treated-and-flagged record is classified ineligible and silently leaves the
cohort. The free-text normalisation problem is retired twice over — the Stage 0 note records a
`Clnician` typo and parenthetical annotations that would have made an exact-string classifier raise on
four of five spellings of one reason.

**Retained means not ineligible, and under two classes that is the same as `== eligible`.** The PI's
predicate decision of 2026-08-10 chose the same population DECISION 1a reaches directly. The code
still offers exactly one predicate, `eligibility.retained()`, over one declared tuple — not because
the two readings can diverge today, but because that tuple is the one seam a third class would come
back through, and a [§13] arm over the undocumented controls would reopen it.

**Accept when.** A cross-tabulation of eligibility by centre and arm is produced, every declared
centre × arm cell rendered whether or not the data fills it; a test confirms that **blanking the
reason column on every record leaves the classification unmoved**, which is what catches a partial
revert of DECISION 1a; and the retained predicate is asserted to follow a patched registry rather than
a comparison, so it stays registry-driven rather than a coincidence of the current two classes.

## Stage 5 — Cohort construction [§2, §3]

**Spec:** `specs/stage5_cohort_construction.md`.

**Build.** Two restrictions applied in order, each recording what it removed **by naming the cases**:
1. Drop centres with zero treated patients. The predicate is **computed from the frame**, never
   declared as a list of centre names — a centre's treatment availability is a property of the data.
2. Drop ineligible patients, through `eligibility.retained` and never through a comparison of this
   stage's own.

**Written against DECISION 1a (PI, 2026-08-13): eligibility is two classes, from
`ivt_contraindicated` alone, and a patient whose reason was never documented is `eligible`.** That
amendment **landed the same day**, across nine files; `specs/stage5_cohort_construction.md` §20 is
the record of where each item went. It does not change the cohort: the same 93 patients either way,
verified identifier for identifier.

Emit a cohort-flow table with a line reporting the retained **control-arm** patients who carry no
documented contraindication reason — the 43 whose eligibility rests on the flag alone — so the group
DECISION 1a makes indistinguishable stays countable.

Then call `derive_cohort` — **here and nowhere else** — so the [§13] median is frozen on the cohort.

**Accept when.** No ineligible patient remains; the cohort-flow table carries the no-reason-on-file
count as its own line, arm-restricted; and:

- **Every centre in the resulting cohort contains both arms, enforced as a runtime raise rather than
  confirmed by a test.** It can only fail through restriction 2 removing a centre's whole control arm
  — the mirror of the [§3] restriction-1 violation, in the direction [§3] does not name. It is
  unreachable on v7, so a test-time reading would ship an unenforced criterion. Dropping such a centre
  would extend [§3] and keeping it puts a structurally non-positive stratum in the propensity model:
  raise, and amend [§3] with the PI.
- **`derive_cohort` is called by the cohort builder, after both restrictions, and the builder is its
  only caller.** A call placed early produces a full, plausible `core_above_median` describing a
  different subgroup, with nothing failing [Stage 3 §6.4]. Enforced on the **row count**, not only on
  the median: the entry's `records` cell must equal the cohort's size. The median is 6.0 mL over the
  full frame but 5.0 over *both* the post-restriction-1 frame and the cohort, so a call misplaced
  between the two restrictions is invisible to a median assertion on this workbook
  [Stage 5 §4.5, §21 R1].

## Stage 6 — Propensity score and weights [§7]

**Spec:** `specs/stage6_propensity_and_weights.md`.

**Two modules, not one.** `model.py` holds what is outcome-agnostic — the design-matrix builder and
the Firth fitter — and `propensity.py` holds the [§7] specification. Stage 9's `m_a(X)` and Stage 12's
standardisation model enter `model.py` directly, so neither imports a module named for the exposure in
order to fit an outcome model [Stage 6 §0.1].

**Build.**
- A design-matrix builder: reference-coded dummies for factors, constant columns dropped (a resampled
  or subset cohort can leave a factor level empty and make the design singular). The level set is
  **declared** (`FACTOR_LEVELS`) and the reference column is dropped **by name**, never `drop_first`.
  The reason is not the obvious one, and an earlier version of this bullet gave the obvious one: on a
  *declared* `Categorical`, `drop_first=True` drops the first **declared** level, which is the
  reference, so the two forms are equivalent and neither rebaselines. Measured. By-name is kept because
  it does not depend on `REFERENCE_LEVELS[c] == FACTOR_LEVELS[c][0]` — a coupling nothing asserted until
  Stage 6 added it to `test_config.py`. A replicate that loses the reference level entirely fails as a
  **rank deficiency**, not as a `KeyError` [Stage 6 §4.2, §20 finding 2].
- **Firth-penalised logistic regression.** Maximise `l(b) + ½ log|I(b)|`; modified score
  `X'(y − p + h(0.5 − p))` with `h` the hat-matrix diagonal. Step-halving on the penalised likelihood;
  converge on the penalised likelihood change or a flat modified score — both criteria fire, the second on
  designs whose information matrix is ill-conditioned. **The rule is deliberately not the reference
  implementation's, which is now read rather than guessed:** `logistf` requires the likelihood change,
  the score *and* the coefficient step conjunctively, where this one takes either of the first two and
  excludes the third — because under near-separation the surface is flat, so a step-norm criterion would
  fail a finite correct fit and [§10] would drop exactly the sparse replicates [Stage 6 §5.3, §18g]. **Raise on failure — never fall
  back to a different estimator** [§7].
- Overlap weights `w = 1 − e` treated, `w = e` control.
- Kish effective sample size per arm.

**Accept when.** Tests confirm: at large n on well-behaved data the coefficients match an unpenalised
maximum-likelihood fit; under complete separation the fit stays finite where the unpenalised one
diverges; the returned probabilities are finite and in (0, 1). And:

- **Complete-case on the [§6] covariates, with the missingness reimposed rather than resolved** [§11].
  The fit runs on the covariate-complete subset — 92 of the 93 cohort patients on v7, the exception a
  control at Lugano missing both CTP volumes — and `e` and `w` come back as Series over the **whole**
  cohort, `nan` off an explicit `in_model` mask, with the excluded patients named in the log. A
  weight of zero is a patient who was weighed; an absent weight is a patient who was not, and only the
  second gives [§11] its denominator. The ATO population is therefore the complete-case set, which is
  stated rather than left implicit [Stage 6 §4.4].
- **Three silent failures are runtime raises, not test-time checks**, because none is caught by a
  check on the output. A `pd.NA` in an `Int64` covariate becomes `nan` in the design **silently**, and
  the fit then runs to completion returning all-`nan` coefficients with a plausible iteration count;
  a factor value outside `FACTOR_LEVELS` becomes an all-zero dummy row, which is the encoding of the
  *reference* level, so the record is modelled as though it were at the reference and the fit returns
  a finite number; and a **response** outside `{0, 1}` — an ordinal `mrs_90d` handed to the logistic
  fitter — fits, converges and returns coefficients for a model of something else, with nothing missing
  anywhere. The third guards the module boundary Stages 8, 9 and 12 enter directly, and was added by the
  engineering review. All three measured [Stage 6 §3.1, §4.5, §5.4, §20 finding 3].
- **A new `model` audit kind carries the [§7] record**, so "record in the output" has a place: ESS per
  arm, the weighted-population description, and the statement that both are conditional on this
  propensity specification because the ATO target population is `h(X) = e(X){1 − e(X)}` — in the
  reproducible log, beside the fitted coefficients, rather than in a return value a reporting layer may
  drop [Stage 6 §7].

**On the estimator being in-repo.** There is no maintained Firth implementation for Python:
`statsmodels` has none, `firthlogist` was last released in 2022 and requires scikit-learn (excluded at
Stage 1), `pyfirth` is a two-line stub, and R's `logistf` cannot be a dependency of a 2000-replicate
bootstrap. Measured rather than assumed. So the estimator is implemented here and carries **five
oracles** instead: `firthlogist` as a dev-only reference, a dependency-free golden vector, a
`scipy.optimize` check on the same penalised objective, the unpenalised MLE at large n
[Stage 6 §16b], and a gated subprocess comparison against R's `logistf` and `PSweight`. The fifth is
the only *independent* one — the first is a port of `logistf`, the second and third check us against
ourselves, and the fourth bites only at large n — and it agrees with this implementation on the
penalised optimum to **1e-9** [Stage 6 §16b.2, §18g]. **Raise on failure — never fall back to a different estimator** [§7] is unchanged and
is now asserted by import scan.

## Stage 7 — Balance and overlap diagnostics [§9]

**Spec:** `specs/stage7_balance_and_overlap.md`.

**Build.**
- Standardised mean differences before and after weighting, using the **unweighted** pooled SD as a
  common denominator so the yardstick does not move.
- Judged against the full [§6] confounder set regardless of which covariates a given specification put
  in its own propensity model.
- The **four** vascular risk factors reported as negative controls [§6]. `BALANCE_ONLY` holds
  **five** names and the fifth is not one of them: `penumbra_ml` is excluded from the propensity
  model because it is core and Tmax>6 s restated, and a covariate excluded for exact collinearity is
  not a covariate weighting was expected to fix. The role column is computed from
  `PS_COVARIATES_FULL`, which is [§6]'s own identity, rather than from `BALANCE_ONLY`
  [Stage 7 §4.3].
- A **within-centre** overlap table alongside the pooled one: per centre, numbers per arm, propensity
  range per arm, ESS, maximum weight, weight share, worst within-centre |SMD|. Centres with structural
  non-positivity are reported as such rather than given an overlap plot.

**Accept when.** The SMD function reproduces a hand-computed value, and returns missing rather than
zero for a variable observed in only one arm. And:

- **One balance row per declared factor *level*, not per design-matrix column** — including each
  factor's reference level and a declared level nobody has. `model.design` drops the reference dummy
  by name and then drops constant columns, so a table built from its columns is three rows shorter on
  this cohort and omits `center_HUG`, `onset_type_witnessed` and `center_USZ`. The first of those
  carries the **worst residual imbalance in the cohort**, so the shortcut produces a complete-looking
  table whose reported worst is wrong, with nothing raising. Stage 7 therefore builds indicators from
  `FACTOR_LEVELS` and does not import `model` [Stage 7 §4.2, §12.2].
- **A covariate constant within each arm and different between them returns missing, not zero.** That
  is a covariate which perfectly separates the arms, and `pilots/analysis.py:194-208` reports it as
  `0.0` — passing the threshold as the best-balanced row in the table — because its only guard is
  `sd > 0`. This is the second reading of the criterion above, and it is the branch the pilots get
  wrong; the branch the criterion names, they already get right. The comparison is of the two **arm
  constants** and never of the two weighted means: measured, `np.average` of one constant under two
  weight vectors differs at 1 ulp, and the weighted-mean form reports a covariate identical in both
  arms as undefined [Stage 7 §5.3, §5.3a].
- **The pooled overlap report is a row of the within-centre table, reconciled against Stage 6 rather
  than recomputed.** Stage 6's `overlap_weights` entry already carries the pooled per-arm `n`, ESS,
  max weight and weight share, so a second table would be a third rendering of them in the same log.
  Building the pooled row with the same function as the centre rows is what stops the two disagreeing
  about what a column means, and no second effective sample size exists anywhere: `propensity.ess` is
  called per centre per arm and a Kish sum appears nowhere in the module [Stage 7 §6.1, §6.3].

**What this stage found, and it is not a code defect.** [§9]'s |SMD| < 0.10 threshold is **exceeded
after weighting, by covariates that are in the propensity model**: measured 0.211 on `center = HUG`
and 0.154 on `center = Lugano`, where Stage 6 §13 carried the pilots' synthetic expectation of
`0 < worst < 0.05`. With an unpenalised score every in-model row balances to 1e-15, so the residual
is the [§7] estimator's own — Firth's modified score solves a different equation from the one the
exact-balance property rests on. It is reported and handed to the PI under [§9] and [§13]; it is
**not** a licence to substitute an estimator [§7], and Stage 6 §12.10's pair fails if one is
silently reinstated [Stage 7 §6.4, §13].

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
