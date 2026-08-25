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

**Spec:** `specs/stage8_primary_outcome_estimator.md`.

**Build.**
- A **weighted proportional-odds model with treatment as the sole predictor**. Most implementations
  have no observation-weight support — measured: `statsmodels`' `OrderedModel` has none at all, and
  there is no maintained weighted alternative for Python — and weighting the log-likelihood
  contribution per observation is sufficient **for the point estimate and not for a standard error**,
  because the weights are a tilting function of an *estimated* propensity score rather than
  frequencies. The inverse observed information therefore omits the variability of estimating `e`,
  which [§7] states enters the interval, so no standard error and no interval leaves this stage:
  [§10]'s percentile bootstrap is the prespecified one and there is no second
  [Stage 8 §5.6, §21 item 6].
- **Weighted empirical cumulative risk differences** `RD_k`, `k = 0…5`, taken from the weighted
  distributions rather than from threshold-specific models, so the cumulative probabilities stay
  ordered by construction.
- Two modules, as at Stage 6: the fitter is `model.py`'s and is **general in its covariates from the
  first line**, because [§14a] prescribes the same estimator with a wider design and unit weights, so
  Stage 12 reuses it rather than a second proportional-odds implementation existing. The [§8]
  specification is `outcome.py`'s, and Stage 9's binary estimators extend that module so the weighted
  per-arm proportion has one home [Stage 8 §0.1].

**Accept when.** The orientation test confirms `exp(β) > 1` when treatment shifts mRS downward — **with
a sign-flipped companion shown to report the opposite direction**, because every reference
implementation parametrises `α_k − x'β` where [§14a] writes `α_k + βA`: measured, the coefficients come
back negated and the cutpoints do not, so two fits that disagree about the primary result look
identical in four of five printed quantities, and `pilots/analysis.py:475` carries the negation that is
correct for `statsmodels` and inverts this study's result [Stage 8 §7]. And:

- **A fit whose treatment coefficient reaches a declared bound raises rather than returning**, because
  an unpenalised proportional-odds fit on a separated sample **converges**. This is not what the
  criterion below originally assumed and it is the finding the stage arrived with: measured, a
  perfectly separated 20-against-20 frame converges in **17 iterations on the score criterion** with
  every safeguard counter at zero, every fitted probability finite and nothing anywhere out of range —
  and returns `exp(β) = 6.5e15`. Four hundred sparse replicates of the cohort's
  shape produced **zero** convergence failures, because there are none to produce. So [§10]'s "replicates
  whose prespecified fit fails are dropped and counted" has nothing to count unless the bound raises,
  and without it a percentile interval is a quantile of a distribution with a tail at `exp(21)`
  [Stage 8 §6]. The bound is prespecified in `config.py` for the reason the Firth tolerances are, and
  its value is chosen from a measured **empty band** rather than from any estimate: over 4800 fits at
  twelve true effect sizes, no non-degenerate fit exceeded `|β| = 8.79`, no degenerate one came below
  `18.81`, and nothing at all landed between the two.
- **Integer weights equal an unweighted fit on the row-replicated frame**, and *that* is the weighting
  criterion. "The weighted fit equals an unweighted fit when all weights are 1" is kept as its
  companion and is not sufficient alone: measured, `polr(X, y, ones)` equals `polr(X, y, None)` bit for
  bit, so an implementation that never reads the weights satisfies it. The replication form agrees at
  **0.000e+00** and is the only available oracle that checks the weighting at machine precision
  [Stage 8 §14.3, §18b].
- **The cumulative probabilities are monotone in `k` within each arm — by construction, so this is not
  the assertion that protects the route.** `P_w(Y ≤ k)` is a cumulative sum of non-negative weights over
  a `k`-independent denominator, so it cannot fail on a correct implementation: measured non-monotone in
  **0 of 2000** samples. And it does not reliably fail on an incorrect one either — the
  threshold-specific route crosses in only **8 of 1990** samples and by less than 5e-07. What must be
  asserted instead is that the route taken *is* the empirical one: the `RD_k` against a hand computation,
  and against the six-threshold-model route, shown to **differ** [Stage 8 §8.2].
- **The estimate's [§11] denominator is a third population and is named, counted and logged** — the
  cohort's 93, then the ATO's covariate-complete 92 [Stage 6 §4.4], then those of *those* whose mRS is
  present. On v7 the third equals the second, which is exactly the condition under which an
  implementation that never built the mask is green [Stage 8 §4.1, §14.2].

**On the estimate not being in the specification.** Stage 8's spec deliberately does **not** measure or
quote `β`, `exp(β)` or any `RD_k` on the workbook, and its data-gated acceptance tests assert
properties — converged, ascending cutpoints, bounded coefficient, monotone `RD_k`, orientation
agreeing across the two scales — rather than values. Every earlier stage recorded its decisions as
taken before any outcome was examined by arm, and Stage 8 is where that necessarily ends; the choice is
whether it ends while the estimator is being specified or after, against a specification already
committed. It ends after. The regression pin is a synthetic golden vector instead, the estimate lands
in the gitignored log, and the cost — no pinned regression number on the primary effect anywhere in
git — is stated rather than minimised, with the trigger under which the pins become correct filed in
`TODOS.md` [Stage 8 §4.3, §15].

## Stage 9 — Secondary binary estimators [§8]

**Spec:** `specs/stage9_secondary_binary_estimators.md`.

**Build.**
- Weighted risk difference: difference of weighted means.
- Weighted marginal odds ratio, **and "kept finite" is a rule and not an instruction**: when a weighted
  proportion reaches 0 or 1, `OR_CONTINUITY` = 0.5 is added to all four weighted pseudo-counts —
  Haldane–Anscombe — and the correction is **flagged and printed beside the estimate**, with the two
  uncorrected proportions, wherever the estimate appears. It fires only on a degenerate proportion, so
  an interior estimate is bit-for-bit uncorrected, and it is never applied to a `nan` proportion: an
  arm carrying no positive weight has no cell to correct. Measured: the branch is **unreachable on v7**
  and fires in 17.0% of stratified replicates for `sich` and 13.0% for `tici_2b_3`, so this is
  materially a choice about two *intervals* and barely one about any point estimate. The alternative —
  report the odds ratio as undefined — is the [§8] convention this repository follows elsewhere and is
  recorded as PI-reversible [Stage 9 §6.3, §6.4, §17].
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
- A guard: **outcomes whose minority cell — `min(events, non-events)` — is below 10 *and which declare
  no override* are not augmented** [§8, as amended]. Keying this on the event count alone would let
  TICI 2b–3 through with 6 non-events. **An earlier form of this bullet said "below 10 are not
  augmented" and that was wrong**: the [§8] amendment's item 2 says such an outcome is augmented with a
  declared reduced `m_a(X)` *"rather than dropped from augmentation"*, and **the roadmap contradicted
  itself twice over**: the override bullet above requires TICI's reduced list to be printed beside its
  estimate, which is an instruction about an augmented estimate this bullet had removed — and Stage 0's
  DECISION 3 already states the rule correctly, that TICI *"is augmented with a reduced `m_a(X)` ...
  rather than dropped from augmentation"*. So an override is a route **in**, and the two paths are
  disjoint: TICI is reduced-and-augmented, `sich` and `ph2` are unaugmented [Stage 9 §7.2, §9.2].
- The guard is applied on **that outcome's own [§11] population**, not the workbook — it is the
  population the nuisance model is fitted on. The two differ: `tici_2b_3`'s events are 114 on the
  workbook and 85 on the ATO population, while its non-event count is 6 on both, which is exactly the
  coincidence that would let a wrong choice go unnoticed [Stage 9 §9.1].
- **The augmentation path is decided ONCE, on the point estimate, and passed into the bootstrap.** A
  minority cell is a property of the sample, so a replicate can cross the threshold — measured, `ph2`
  crosses in **46.1%** of stratified replicates and `sich` in 3.8%. A rule re-evaluated per replicate
  would make `ph2`'s percentile interval a quantile over a near-even mixture of two estimators, which
  is [§10]'s "never substituted with a different estimator" violated by a different route. A replicate
  whose declared model then fails to fit is a `FitError`, dropped and counted, and never a silent
  downgrade to the unaugmented form [Stage 9 §9.5].
- **And the mechanism for passing them in is Stage 9's, not Stage 10's**: `outcome.secondary` takes an
  optional `paths` map, complete or absent and never partial, so the contract above is expressible as a
  call rather than as an instruction. Stage 10 reads `augmented_path` off the point estimate's seven
  estimates and hands the map back; it never re-derives the rule. A function whose signature could not
  express this left Stage 10 only two options, both prohibited: recompute the paths, or reimplement the
  seven-outcome loop [Stage 9 §11.1].
- **No bound on the nuisance coefficient, and Stage 8's reason for having one does not transfer.**
  Separation in `model.firth` converges silently exactly as it does in `polr` — measured, a separated
  40-record frame converges in 6 iterations with every safeguard counter at zero and returns
  `exp(β_treatment) = 441.005397`, while the *unpenalised* MLE on the same frame has no maximum at all
  and raises `LinAlgError: Singular matrix` — but Stage 8 could calibrate a bound from an *empty band*
  and there is none here: `sich`'s `max_abs_beta` runs continuously from 0.8558 to 80.9825 across
  replicates, 42.7% above 8, with no gap. A bound is not needed anyway, because `m_a(X)` is a nuisance
  whose *predictions* enter `tau` and predictions are probabilities: every term is bounded, so there is
  no `exp(β)`-style tail. What is lost instead is precision, and the replacement for a guard is a
  diagnostic — Stage 10 reports the `max_abs_beta` distribution [Stage 9 §9.6].
- **The augmentation tilt has one home and the call site is tested.** `h = e(1 − e)` is `outcome._tilt`
  and not a sub-expression: every test that distinguishes `h` from `w` calls the estimator directly with
  its own arrays, so without an assertion on what `secondary` *passes*, tilting by `w` at the call site
  ships the wrong estimand with the whole suite green [Stage 9 §8.2a, §15.8 item 6].
- **`model.predict` is the fitter's own expression evaluated elsewhere**: the intercept is prepended as a
  column and the linear predictor is clipped at `FIRTH_ETA_CLIP`, so `predict(fit, X) == fit.p` **exactly**
  on the design the fit was made from. A scalar-intercept form is `allclose` but not equal, and an
  unclipped one diverges from the fitter above `|η| = 500` [Stage 9 §7.3, §15.6].
- **The odds-ratio correction can cross the null, not merely widen.** The four pseudo-counts are
  weight-sums and the arms' totals differ, so an empty cell in the *heavier* arm can carry the estimate
  past 1 — measured, `Sw` 17.007 against 22.768 with `p0 = 0.00647` returns 1.0201 where the uncorrected
  value is 0.0, in about 0.9% of empty-cell replicates. Roughly 3 replicates in `N_BOOT` = 2000 report a
  safety outcome with no bridging events as favouring EVT alone. PI-reversible [Stage 9 §6.4, §17].

**Accept when.** Three tests pass, and the first one's original rationale was wrong:
1. With the outcome model set to **two different constants**, the augmented estimate equals the
   unaugmented weighted risk difference exactly. **This checks that the correction term is normalised
   and cancels — it does *not* check the tilting function**, and an earlier form of this item claimed it
   did. For constants `m₁ ≡ c₁`, `m₀ ≡ c₀` and any normalised weight `u`, `Σu(c₁−c₀)/Σu = c₁−c₀`
   independently of `u`, so the identity holds whether the augmentation tilts by `h`, by `w` or by unit
   weights — measured, all three at ~2e-16. The constants must **differ**: with `c₁ = c₀` a dropped
   correction term passes at +0.000000, where with `c₁ = 0.62, c₀ = 0.23` it returns `rd − 0.39`
   [Stage 9 §8.4].
2. **The tilting function is pinned by two tests and neither is test 1.** With a *non-constant* `m_a(X)`,
   the estimate under `h` is asserted and the estimate under `w` is computed and shown to differ —
   measured `1.030664e-02` against a `1e-03` bound, with the constant-`d` sanity cases cancelling
   exactly at `0.000e+00`. **And the call site is asserted separately**, because test 2 exercises the
   estimator rather than what `secondary` hands it. The fixture must be **written out with its seed**:
   the gap is driven by the variation in `m₁ − m₀` but is *not monotone* in it — one of six measured
   candidates had the largest spread and the second-smallest gap, at 2.5e-04 — so "make the effect
   heterogeneous" is not a specification [Stage 9 §8.4, §15.0.4, §15.8 items 3 and 6].
3. A heterogeneous-effect scenario with a **correct** outcome model and a **misspecified** propensity
   model does *not* recover the true ATO — and *does* recover the ATO indexed by the misspecified score,
   which is the positive half that distinguishes "not doubly robust" from "wrong". This must be
   asserted, not merely allowed: the estimator is not doubly robust [§8], and a test built on a
   homogeneous effect will appear to show that it is, because every weighted average treatment effect
   coincides in that case. Two corrections to how it is run: **"homogeneous" means constant on the
   risk-difference scale** — a logistic outcome model with no interaction is not that, and is itself
   detectably biased at `|t| = 21.3` — and **the test cannot be run at the cohort's size**, where the
   per-replicate spread is sd 0.0877 against a bias of 0.026 and 40.5% of replicates have the wrong
   sign. Run at n = 200,000 with the seed written out, where twelve seeds put the heterogeneous arm's
   error in [0.02295, 0.030249] and the constant-effect companion's in [0.000176, 0.004929], with the
   band between them empty and both asserted bounds (0.010, 0.015) inside it [Stage 9 §8.5, §15.10].

**And every fixture these three tests run on is specified as code**, with its seeds and constants written
out. This is not a style note: the spec's first two drafts pinned quantities on three fixtures to ten
decimal places and described rather than gave them, which made tests 1-3 unstartable and the acceptance
criteria unperformable [Stage 9 §15.0, §22.3 item 1].

**Naming.** "model-assisted", never "doubly robust", in code, docstrings and output.

**What Stage 9 found that is Stage 10's, and it is a blocker.** `propensity.fit` raises `SchemaError` on
**26.0%** of patient-level stratified replicates: `_record_exclusion` (`propensity.py:400-432`) compares
distinct excluded `case_id`s against excluded rows, and its premise is that duplicates are forbidden by
Stage 2's A2 — which forbids them *in the workbook*, not in a resample. The cohort has exactly one
covariate-incomplete record, in Lugano's stratum of 31, and the analytic probability of drawing a given
row twice from 31 is 26.4%, which matches. Since Stage 8 §11 hands over that a `SchemaError` is a
resampler bug and **must not be caught**, Stage 10's three options as landed are all prohibited: crash,
catch what it prespecified it would not, or drop a quarter of its replicates for a reason that is not
about the data. Either the guard learns that a replicate's rows are distinct draws, or the resampler
names them distinctly. This is the first thing Stage 10 will hit [Stage 9 §12.3, §14].

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

**Three things Stage 8 hands over rather than leaving to be found here** [Stage 8 §11].
- **Two failure counters, reported separately.** A dropped replicate is either a non-convergence or a
  separation guard (`G7`), and the two mean different things about the data. Measured at Stage 8: 400
  sparse replicates of this cohort's shape produced **zero** convergence failures, because an
  unpenalised ordinal fit on a separated sample converges rather than failing. **A failure counter
  reading zero is therefore not evidence that no replicate was degenerate** unless the G7 count is
  reported beside it.
- **The distribution of `len(fit.alpha)` across replicates**, alongside those counters. The number of
  fitted cutpoints is a property of the replicate — Stage 8 collapses the response to the levels
  carrying positive weight — so a replicate missing a declared mRS level contributes a `β` on a coarser
  scale, and percentiles are taken across them. Under proportional odds they are the same parameter,
  which is [§8]'s own untested assumption. **And Stage 8's separation bound was calibrated at seven
  occupied categories only**, so if that distribution is not degenerate at six cutpoints, the drop rate
  is partly a function of a bound measured on frames unlike the ones being dropped. Reporting the
  distribution is what makes that answerable from the replicates already drawn
  [Stage 8 §5.3, §6.3, §15].
- **`FitError` is the droppable failure and `SchemaError` is not.** A `SchemaError` from `primary` or
  `propensity.fit` is a bug in the resampler, not a sparse replicate, and catching it would drop
  replicates for a reason that is not about the data (`propensity.py:468-469`).

**And four things Stage 9 hands over** [Stage 9 §14].
- **`propensity.fit` raises `SchemaError` on 26.0% of stratified replicates as landed**, and by the rule
  immediately above, Stage 10 may not catch it. Fix the resampler or the guard before anything else;
  Stage 9's own replicate measurements had to work around it and say so.
- **The augmentation paths are passed in, not recomputed** — `ph2` crosses the rare-minority threshold
  in 46.1% of replicates.
- **The distribution of `max_abs_beta` from `m_a(X)`, per outcome, reported beside the failure counters.**
  The name is ASCII deliberately: it becomes a column in a rendered markdown table, `data._md_table` does
  no escaping, and a literal pipe in a cell breaks the table in a way byte-identity checks do not catch
  [Stage 9 §10.4, §15.14].
  Stage 9 prescribes no bound, so nothing turns a separated nuisance fit into a countable failure: the
  `FitError` rate is 0.3% for `sich` and 0% elsewhere while 15.1% of `sich`'s fits exceed `|β| = 14`. The
  counters will read near-zero on outcomes whose nuisance models are degenerate a sixth of the time.
  Same move as `len(fit.alpha)` above, same reason.
- **Stage 9 is the cost centre, not Stage 8 — and `model.design` is 72% of Stage 9.** Measured:
  24.05 ms per replicate for the five augmented outcomes against `propensity.fit`'s 13.39 ms and
  `primary`'s 4.73 ms; 71 ms wall-clock, about 142 s over `N_BOOT`. The four full-list outcomes share
  an **identical** design matrix, so building it once per replicate is worth roughly 22 s.

**Accept when.** On synthetic data with a known effect the interval covers the truth at roughly the
nominal rate; the failure counter is exercised by a deliberately degenerate input; and no code path
can return an estimate from a fallback estimator. And:

- **The separation count and the convergence count are both surfaced, and the separation one is
  exercised.** A test drives a replicate loop over a deliberately separable frame and asserts the G7
  count is non-zero while the convergence count stays zero — the pair, because Stage 8 measured that the
  second never fires on this estimator and a single "failures" counter therefore reads zero on data that
  is degenerate throughout [Stage 8 §6.1, §14.7].

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

- **The full-covariate propensity sensitivity analysis** [§13, DECISION 4]. Refit the [§7] propensity
  model over `config.PS_COVARIATES_FULL` — the [§6] set plus the four vascular risk factors — re-derive
  the overlap weights, and re-run the [§8] primary estimator on them. Report the common odds ratio and
  the six `RD_k` beside the primary ones, **with ESS and worst residual |SMD|** as [§13] requires, which
  is the axis that prompted the arm.
  - **It reuses every landed stage and adds no estimator.** `propensity.fit` takes no covariate list —
    Stage 6 §6.4 made that a prespecified choice rather than an option — so this arm needs a seam for
    the covariate set, and that seam is the one design decision here. `outcome.primary` then runs
    unchanged on the resulting `Propensity`, because it reads `e`, `w` and `in_model` and nothing else.
  - **It spends the negative controls, by construction and not by accident.**
    `config.NEGATIVE_CONTROLS` is computed as `PS_COVARIATES_FULL` minus `PS_COVARIATES`, so the four
    covariates this arm adjusts for are exactly the four that stop being controls in it. Say so where
    the arm's balance table is rendered; `penumbra_ml` is the only balance-table covariate left outside
    every model.
  - **Its own balance table has the same 19 rows with four of them re-roled**, which is what Stage 7
    §4.1's role column exists for.

**Accept when.** The arm's `in_model` population, ESS and worst residual |SMD| are reported beside the
primary's, and the two estimates appear together with the statement that they differ in the propensity
specification and in nothing else. **Balance is still judged against the full [§6] confounder set**
[§9], so the arm's table is comparable with the primary's row for row.

**Still deferred.** The other **five** [§13] sensitivity rows — propensity model without centre;
unadjusted; largest centre alone; pre-stroke mRS ≤ 2; model-assisted augmented cumulative mRS risk
differences. Leave the hook but do not build the specifications. **The without-centre row is deferred
deliberately rather than by omission** (DECISION 4): dropping centre removes the near-separation at its
source and would balance the clinical covariates nearly exactly, while abandoning adjustment for the
variable that most nearly determines treatment — 30 of 41 bridged at HUG against 2 of 31 at Lugano — so
it is very likely a *more* confounded estimand and agreement with the primary must not be read as
reassurance.

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
- **Any output carrying `exp(β)` carries the constant-shift statement** [§16 amendment, DECISION 5].
  The primary measure is a *common* odds ratio; the assertion is that the assumption is untested, and
  — computed from `Primary.rd`, not asserted by hand — that the six `RD_k` do not all share a sign
  where that is true of the run. Deriving it from the result object is what stops the statement and
  the table disagreeing after an edit.
- **Any table of raw event rates by arm is labelled descriptive and not an unadjusted estimate**
  [§16 amendment, DECISION 5]. Treatment is nearly determined by centre in this cohort, so a pooled
  crude comparison can differ from the within-centre one in **direction**; a reader taking the crude
  column as an unadjusted estimate can get the sign backwards. The table still appears — STROBE and
  RECORD require it — and the guardrail is on its label.

**Accept when.** Both statements are emitted from the result objects rather than written into a
template, and a test constructs a run in which the `RD_k` **do** all share a sign and asserts the
one-directional clause is absent — so the sentence is known to track the data rather than being
printed unconditionally.

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
