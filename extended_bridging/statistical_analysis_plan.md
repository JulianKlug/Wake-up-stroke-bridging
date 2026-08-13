# Statistical analysis plan

**IVT before thrombectomy versus thrombectomy alone in late-window, CT-perfusion-selected stroke**

## 1. Objective

To estimate the effect of intravenous thrombolysis before thrombectomy (bridging) versus thrombectomy
alone on functional outcome, among late-window patients eligible for either strategy. No randomised
evidence exists here: all bridging RCTs and the IRIS meta-analysis randomised within 4.5 h.

## 2. Target trial

| Element | Specification |
|---|---|
| Eligibility | LVO, late window, CTP-selected, eligible for thrombectomy, no absolute contraindication to IVT, at a centre using both strategies |
| Time zero | The IVT decision point, after initial imaging |
| Strategies | (a) IVT then thrombectomy; (b) thrombectomy alone |
| Assignment | Not randomised; emulated by conditioning on §6 covariates measured at or before time zero |
| Outcome | 90-day mRS |
| Estimand | ATO — average treatment effect in the overlap population |
| Follow-up | Time zero to 90 days |

Eligibility, assignment and follow-up are aligned at time zero. No covariate measured after time zero
enters any model.

**Known limitation.** Entry requires receipt of thrombectomy, which is downstream of the IVT decision.
Patients whose occlusion resolves after IVT so completely that thrombectomy is cancelled never enter
the cohort. This biases against bridging and must be stated in the manuscript.

## 3. Population and cohort restrictions

Two restrictions define the primary cohort, both applied **by design, before any weighting**:

1. **Exclude centres that contributed no bridging patients** (expected: USZ). There,
   `P(IVT = 1 | centre, X) = 0` structurally. No weighting recovers a contrast that was never
   available, and a penalised score would report shrinkage as though it were treatment availability.
2. **Exclude patients with an absolute contraindication to IVT.** They were never candidates for
   bridging; retaining them makes "no IVT" a marker of contraindications and their prognosis.

Eligibility is classified using only information available at time zero, applied symmetrically across
arms, and fixed before outcomes are examined. **The classifier is the recorded absolute-contraindication
flag, and it takes two values.** A patient carrying it is ineligible and is removed by restriction 2;
a patient not carrying it is eligible. The free-text reason field is not used — neither its content nor
whether it was completed. No patient who received IVT carries the flag; this is asserted rather than
assumed, since a patient recorded as both having received IVT and having an absolute contraindication
to it is a contradiction to resolve with the data owner, not a class to assign.

Restriction 1 is **one-directional by design**: it names centres with no bridging patient, which is the
violation this dataset presents. A retained centre left with no thrombectomy-alone patient by
restriction 2 would be the mirror violation, `P(IVT = 0 | centre, X) = 0`. It does not arise here; should
it arise, it is a protocol amendment and not a choice made in analysis code.

Cohort flow is reported, including the number of retained thrombectomy-alone patients for whom no
contraindication reason was recorded.

**Amendment, 2026-08-13 — eligibility is classified from the recorded contraindication flag alone, and
the *indeterminate* category is withdrawn.** The original text made eligibility three-valued: where the
reason a patient did not receive IVT was never recorded, eligibility was *indeterminate*, those patients
were retained, and a blank was never read as "no contraindication". The reason field does not support
that role. It was completed for every thrombectomy-alone patient at two of the four centres and for none
at the other two, so a category keyed on its presence classifies **documentation practice at the centre**
rather than eligibility, and yields an indicator that is a near-perfect function of `center` — itself a
§6 covariate. The contraindication flag is recorded for all 126 patients and is the field in which the
question was actually asked.

**The population does not change.** The indeterminate group was already retained, so the primary cohort
is the same 93 patients — 39 bridging, 54 thrombectomy alone — at the same centres. Forty-three patients
change label; none changes arm, centre or cohort. Established before any outcome was examined by arm, as
was this decision.

**The assumption is not removed by this amendment. It is relocated, and in its new form it is less
visible, so it is stated here in full.** Eligibility now rests on `flag = 0` meaning the same thing at the
two centres that never recorded a reason as at the two that recorded one for every thrombectomy-alone
patient. The flag is 0 or 1 on all 126 records and never missing, which is equally consistent with its
having been assessed for every patient and with 0 being an uncompleted default; nothing in the data
distinguishes the two. This governs **43 of the 54 thrombectomy-alone patients in the primary cohort**.
If any of those 43 in fact had an absolute contraindication, they were never candidates for bridging and
their prognosis is worse on exactly the grounds restriction 2 gives, so retaining them depresses the
comparator arm and **biases in favour of bridging** — the opposite direction from §2's known limitation,
which biases against it. Both are stated in the manuscript.

**It is stated rather than bracketed, and the reason is quantitative.** The obvious sensitivity analysis
— restrict the comparator arm to patients with a documented reason — leaves 11 of them, all at one
centre, and turns the primary comparison into a single-centre contrast of 30 against 11. That is a
different study, not a sensitivity analysis of this one. The cohort-flow table therefore reports the
group's size by centre so that a reader can size the assumption directly. No analysis brackets it; it is
stated as a limitation, exactly as the category it replaces was.

One deferred decision closes with this amendment. §14b's active regime assigns treatment as a function of
eligibility, which under the original text was a binary decision over a three-valued column with no
assignment defined for the third value. It is now binary over a binary column.

## 4. Exposure

IVT before thrombectomy versus thrombectomy alone, as recorded at time zero.

## 5. Outcomes

| Family | Outcome | Type |
|---|---|---|
| **Primary** | 90-day mRS | ordinal |
| Secondary | mRS 0–2 · mRS 0–1 · TICI 2b–3 | binary |
| Safety | symptomatic ICH · parenchymal haematoma type 2 · death at 90 days · mRS 5–6 | binary |

All dichotomies are derived from the ordinal source, not from shipped indicator columns.

TICI 2b–3 is partly intermediate between IVT and functional outcome. It is a secondary outcome and no
analysis conditions on it. Outcome ascertainment is not blinded; state as a limitation.

## 6. Confounders

`age` · `sex` · `prestroke_mrs` · `nihss_baseline` · `onset_type` · `core_ml` · `tmax6_ml` ·
`atrial_fib` · `center`. Linear terms only, no interactions.

Excluded, with the reason for each:

- `onset_to_groin_min` — post-exposure (§12). It is the only timing variable available and it cannot
  be used for adjustment.
- `penumbra_ml` — deterministic function of core and Tmax>6 s volumes; including all three makes the
  design exactly singular.
- `hir` — excluded a priori; not used in any model or sensitivity analysis.
- `first_image_mri` — a centre-pathway attribute rather than a patient characteristic, largely
  collinear with `center`. Dropped automatically if constant.
- `hypertension`, `hyperlipidemia`, `diabetes`, `smoking` — weak confounders that do not meaningfully
  predict the IVT decision, and the degrees-of-freedom budget is ~3.5 treated patients per parameter
  already. Retained as balance negative controls, where residual imbalance shows what weighting does
  not fix, and added back in a full-covariate sensitivity propensity model (§13).
- Anything measured after time zero.

Confounders absent from the dataset — ASPECTS, occlusion site, collateral status, transfer status,
calendar period, and **time from onset to presentation**, since the only timing variable is
post-exposure — are named in the manuscript, and an E-value is reported for the primary estimate.

## 7. Estimand, propensity score and weighting

**Estimand: ATO**, among IVT-eligible patients at both-arm centres.

**Propensity score:** Firth-penalised logistic regression on the §6 covariates. One prespecified
estimator, refitted identically in every bootstrap replicate. An unpenalised MLE fails to converge in
a minority of sparse replicates, and switching estimators on failure makes the bootstrap a mixture of
two estimators rather than the sampling distribution of one.

**Weights:** `w = 1 − e` treated, `w = e` control. Bounded, and targeted at patients near equipoise.

Report the effective sample size per arm and describe the weighted population. Because the ATO
population is defined statistically, state explicitly who it comprises.

**The estimand is indexed by the propensity model.** The ATO target population is `h(X) = e(X){1 − e(X)}`,
so changing the propensity specification changes the *population*, not merely the precision of the
estimate. Two consequences: the reported ESS and the description of the weighted population are
conditional on the working model; and the propensity-specification rows in §13 must be read as
different target populations rather than as competing estimates of one quantity. The bootstrap refits
the propensity model in every replicate, so the uncertainty from estimating it enters the interval.

## 8. Outcome analysis

All estimates are marginal in the overlap population. No covariate adjustment of the reported effect:
where an outcome regression appears it enters only as a correction term.

**Primary (ordinal mRS) — weighting only, unaugmented.**

- Weighted proportional-odds model with treatment as the sole predictor; common odds ratio
  `exp(β)`, oriented so > 1 favours bridging. This is the primary quantity, and the primary p-value
  is the test of `β = 0` defined in §10.
- Cumulative risk differences `RD_k = P(mRS ≤ k | bridging) − P(mRS ≤ k | EVT alone)`, `k = 0…5`,
  taken from the weighted empirical distributions — ordered by construction, so no monotonicity
  problem arises — reported as the absolute-scale presentation of the same contrast. **With
  confidence intervals but no p-values:** they describe where on the scale the shift sits, and six
  threshold-wise tests of one contrast would only invite selection of the most favourable cut. The
  primary outcome has exactly one test, `β = 0`.

There is no simple augmented form of a proportional-odds fit, and the alternative of augmenting each
threshold separately would replace one primary number with six and can produce crossing cumulative
probabilities. The primary analysis therefore rests on the propensity model alone.

**Secondary binary outcomes — model-assisted augmented overlap weighting.**

The augmentation uses the same tilting function `h = e(1 − e)` as the weights, so it targets the ATO
rather than the ATE. With `h` throughout:

    tau = Σ₁ w(Y − m₁)/Σ₁ w  −  Σ₀ w(Y − m₀)/Σ₀ w  +  Σ h(m₁ − m₀)/Σ h

Setting the outcome model to a constant reduces this exactly to the unaugmented weighted risk
difference, which is how the two are kept comparable.

**This estimator is not doubly robust, and is not described as such.** The usual "consistent if either
nuisance model is correct" property does not transfer to the ATO, because the ATO *estimand is itself
indexed by the true propensity score*: the target population is defined by `h(X) = e(X){1 − e(X)}` and
corresponds to no observable subset of the data. The ATT and ATC differ in exactly this respect —
their target populations are `{A = 1}` and `{A = 0}`, which do not move when the propensity model
changes. Here, if `e(X)` is misspecified, the tilting function and therefore the estimand are wrong
even when `m_a(X)` is exactly right; the estimator remains consistent for a weighted average treatment
effect, but for the wrong weighting. The formula above additionally treats the propensity score as
known and omits terms the efficient influence function requires once it is estimated.

What the augmentation does buy is reduced residual bias and improved precision when the working
outcome model is informative. Consistent estimation of the ATO still requires an adequate propensity
model. The augmented estimate is therefore reported as **model-assisted**, always alongside the
unaugmented weighted risk difference; material disagreement between the two is evidence about the
outcome model, not confirmation of either.

**Outcome regression `m_a(X)`:** treatment main effect plus **the §6 covariate set**, linear terms,
fitted by Firth logistic. No interactions, no splines, no cross-fitting — at this treated-arm size
flexible nuisance models add variance rather than robustness, and cross-fitting folds would be too
small to serve their purpose.

The two nuisance models carry the same covariates deliberately. The augmentation can only remove
residual bias that `m_a(X)` actually captures, so an outcome model omitting a variable §6 declares a
confounder — one already asserted to predict the outcome — disables the correction term exactly where
it is needed. This is a bias-reduction argument, not a double-robustness argument; see above. The
cross-reference, rather than a second list, keeps the two sets from drifting apart in later revisions.

Reported per binary outcome: augmented risk difference, with the unaugmented weighted risk difference
and the weighted marginal odds ratio alongside.

**Rare outcomes.** Outcomes with fewer than 10 events are **not** augmented: a nuisance model with more
parameters than events cannot be fitted meaningfully. For these, report the unaugmented weighted risk
difference and the weighted marginal odds ratio only. This will normally apply to the safety family,
which §10 already treats as descriptive.

**Amendment, 2026-08-07 — the rare-outcome rule reads the minority cell, and one declared
per-outcome reduction of `m_a(X)`.** The stated reason for the rule is that a nuisance model cannot
carry more parameters than the data support, and that constraint binds on whichever cell is small,
not on the event cell specifically. TICI 2b–3 has many events and few non-events, so it passes the
rule as written while failing its rationale. Two changes follow, both fixed before any outcome was
examined by arm:

1. The threshold is applied to `min(events, non-events)`.
2. Where the minority cell cannot support the full covariate set but can support a smaller one, the
   outcome is augmented with a **declared reduced `m_a(X)`** rather than dropped from augmentation.
   For TICI 2b–3 that model is treatment + `center` + `atrial_fib` — five parameters, chosen on
   procedural grounds (recanalisation tracks device, technique and operator volume, carried by
   centre, and clot composition, proxied by atrial fibrillation) and never from the data.

This is a deliberate exception to the rule above that both nuisance models carry the same covariate
list, and it weakens the bias-reduction argument for TICI exactly as that paragraph describes: the
correction can only remove residual bias the outcome model captures. It is not a double-robustness
claim, which was never available here. The guard is the comparison this section already requires —
the unaugmented weighted risk difference is reported alongside, and material disagreement between
the two is evidence about the reduced outcome model. **Every reduced specification is reported beside
its estimate.** No other outcome is reduced; symptomatic ICH and parenchymal haematoma type 2 have a
minority cell below the threshold and remain unaugmented.

## 9. Balance and overlap

Standardised mean differences before and after weighting, common unweighted pooled SD, threshold
|SMD| < 0.10. Balance is judged against the **full** confounder set, not only the covariates a given
specification put in its propensity model.

Overlap is reported **within each centre as well as pooled**: a pooled distribution can look
acceptable while treatment is nearly determined by centre. Centres with structural non-positivity are
reported as such, not given an overlap plot.

## 10. Inference

Nonparametric bootstrap, 2000 replicates, stratified by centre, refitting the propensity model in
every replicate. Percentile 95% confidence intervals; seed recorded. Replicates whose prespecified fit
fails are dropped and counted, never substituted.

Resampling is at patient level within centre. A centre-level cluster bootstrap is not used: a handful
of centres cannot support cluster-bootstrap consistency. **Inference is conditional on the
participating centres and their observed treatment practices.**

**p-value for the primary outcome.** From the **treatment coefficient `β` of the weighted
proportional-odds model** — the parameter whose exponential is the reported common odds ratio. The
null is `β = 0`. Let `β̂*₁ … β̂*_B` be the coefficient from the refitted model in each replicate; the
two-sided p is

    p = 2 · min{ Pr(β̂* ≤ 0), Pr(β̂* ≥ 0) },   floored at 1/(B + 1)

so it is the smallest level at which the percentile interval for `β` excludes the null, and it agrees
by construction with the reported interval for the common odds ratio. This is a test of the
proportional-odds treatment coefficient, **not** a risk-difference-scale test, and must not be
described as one. Replicates in which the model does not converge are dropped and counted, as above.

`β` is a fitted model parameter and is defined whenever the model converges, so no scale problem
arises. That is not true of the *marginal* odds ratio for a binary outcome, which is a ratio of
weighted proportions and becomes undefined when one reaches 0 or 1; taking p from those draws would
condition on the replicates where the ratio happens to exist — selection on the outcome, which badly
biases rare-event p-values. **For the secondary binary and safety outcomes, therefore, p is computed
on the risk-difference scale**, which is defined in every replicate. Testing `RD = 0` is the same null
as `OR = 1`.

**Estimation, not testing, is the reportable output.** Safety outcomes resting on few events are
descriptive only.

## 11. Data handling

Complete-case, with the denominator reported for every estimate. No imputation. Structural
non-applicability is distinguished from missingness and is never imputed. Dichotomies rebuilt from
their ordinal source with missingness reimposed explicitly. Every correction is content-driven, never
indexed by row, and logged with the affected case identifiers.

## 12. Post-exposure variables

Onset-to-groin time is excluded from every model and reported by arm instead. Bridging can itself
delay puncture, so conditioning on it would remove part of the total effect. Any change in the
treatment estimate when it is added to an outcome model is **descriptive** and must not be read as a
mediated proportion.

## 13. Multiplicity, subgroups and sensitivity analyses

**Multiplicity.** Primary outcome uncorrected. Benjamini–Hochberg FDR within the secondary and safety
families; raw and adjusted p-values both reported.

**Subgroups** (hypothesis-generating, with an interaction test): unknown versus witnessed onset;
core volume above/below median.

**Amendment, 2026-08-10 — the target-mismatch subgroup is withdrawn.** It was listed as a third
subgroup and is removed before any subgroup estimate was produced. The cohort is CTP-selected by
construction (§2), so target mismatch is close to the criterion that admitted these patients rather
than a contrast within them: the subgroup is expected to be near-constant, and an interaction test on
a split of a few patients against the rest is not hypothesis-generating but hypothesis-shaped noise,
reported in a subgroup table where it reads as a finding. The supporting count — core volume is
exactly 0 in 50 of 125 records, so most patients clear the core and ratio criteria and the flag turns
almost entirely on one volume threshold — is marginal and was established in Stage 0 before any
outcome was examined by arm, as was this decision. The two remaining subgroups are unaffected:
onset is unwitnessed or on waking in 93 of the 126 records against 33 witnessed, and the median
split is balanced by construction. `penumbra_ml` keeps its §6 exclusion and its balance-table role,
both of which stand independently of this subgroup.

**Sensitivity analyses** on the primary outcome (Leave out for now), each reporting ESS and worst residual |SMD|:
full-covariate propensity model (adding the four vascular risk factors); propensity model without
centre; unadjusted; largest centre alone; pre-stroke mRS ≤ 2; model-assisted augmented cumulative mRS
risk differences (with a monotonicity check).

## 14. Secondary analyses

**Neither re-runs the primary specification.** A propensity model cannot produce a contrast for a
centre with no treated patients, and the estimands differ from §7. Both use pooled ordinal outcome
standardisation (parametric g-computation). No propensity model is fitted anywhere in §14.

### 14a. All centres, eligible patients

**Estimand:** `E[Y^IVT − Y^NoIVT | eligible]`, averaged over the full eligible cohort at all centres.
This is an **ATE in the all-centre eligible population, not the ATO of §7.** The two must not be read
as one effect estimated with and without extra centres — they target different populations by
construction. State this wherever both appear.

**Model.** One proportional-odds model over all eligible patients from all centres:

    logit P(Y ≤ k | A, X) = α_k + βA + γᵀX,   k = 0…5

with `X` = the §6 covariates **minus `center`**. Centre is omitted because the identifying assumption
is exactly that potential outcomes do not depend on centre given `X`; including it would contradict
the assumption the analysis runs on. Exclusions otherwise as §6; contraindication status is not a
covariate, the cohort being already restricted to eligible patients.

**Standardisation.** Duplicate every patient with `A = 1` and `A = 0`, predict `P(Y = j | A = a, X)`
for `j = 0…6`, average over all N eligible patients. Patients at never-IVT centres contribute their
covariates, their observed direct-EVT outcomes and their share of the target population; their IVT
counterfactual is predicted from relationships learned where IVT was observed.

**Report** the two standardised mRS distributions and, from them: cumulative risk differences `RD_k`,
`k = 0…5`; the standardised mRS 0–2 risk difference; the standardised mortality difference
`P(Y = 6 | IVT) − P(Y = 6 | no IVT)`. `exp(β)` is a **conditional** common odds ratio, not generally
equal to a marginal odds ratio in the standardised population; it may appear as a model parameter but
never as the standardised marginal effect.

**Assumptions.** `Y^a ⟂ C | X` and `Y^a ⟂ A | X` in the pooled eligible population — the centre-level
determinants of IVT use are not independently associated with outcome given `X`. Much stronger than
observing that unadjusted mRS distributions look alike across centres, and it additionally assumes the
conditional effect learned in IVT-using centres transports to never-IVT centres. **Not testable from
these data.** Differences between centre-specific estimates are not a test of it: they are equally
consistent with true heterogeneity, residual case-mix differences, sparse-data variation and model
instability.

**Support check.** Standardisation avoids infinite weights, not extrapolation. Report the treated-arm
distribution of each continuous covariate, the proportion of never-IVT-centre patients outside it, and
a baseline table comparing IVT-treated with never-IVT-centre patients. Predictions far outside the
treated support are driven by the assumed linear form, not by observed comparisons. Sensitivity:
restrict the standardisation population to patients inside the treated support.

**Sensitivity — hierarchical.** Add a random centre intercept `b_c ~ N(0, σ²_C)` with **one common
treatment effect**; no treatment-by-centre interaction and no random slope, neither being identifiable
in any useful sense with four centres. Patients at a never-IVT centre inform that centre's intercept
through their direct-EVT outcomes while the IVT effect is borrowed. This slightly relaxes the
assumption above: if outcomes truly do not differ by centre given `X` then `σ²_C = 0` and it collapses
to the pooled model. A centre variance from four clusters is fragile, and a Bayesian fit would be
steadier but prior-driven — so the pooled model is the implementation and this is a check, not an
upgrade.

### 14b. All patients including contraindicated — a feasible-policy estimand

Contraindicated patients are **not** added to the 14a standardisation: predicting an IVT outcome for
someone for whom IVT was impossible standardises over people who were never in the target population.

The question is instead posed as a policy contrast: **bridge every eligible patient and give direct
EVT to contraindicated patients, versus direct EVT for everyone.** Same standardisation machinery, all
centres, all patients, but treatment under the active regime is a function of eligibility rather than
a constant.

This answers a different question from 14a and from §7, and is reported with that stated. It is the
only §14 analysis whose population includes contraindicated patients, and its interpretation is
operational — what adopting a bridging policy would deliver in this cohort — not the biological effect
of IVT.

### Inference for §14

Patient-level bootstrap stratified by centre, as §10: resample, refit the ordinal model, predict every
resampled patient under both regimes, average, recompute the risk differences, percentile intervals.

## 15. Approaches deliberately not used

Unrestricted or trimmed IPTW; adjustment for post-exposure variables; centre as an instrumental
variable; matching across centres; cross-fitted machine-learning nuisance models; omitting centre from
the propensity model as the primary analysis; centre-specific treatment effects or random treatment
slopes, which four centres cannot identify in any useful sense (§14a).

Omitting centre and modelling across centres are permitted **inside §14 only**, where extrapolating to
centres that never used IVT is the declared purpose and is labelled as such.

## 16. Reporting

STROBE with the RECORD extension. State explicitly: the target trial protocol (§2); cohort flow with
both design restrictions (§3); the estimand and who the ATO population comprises (§7); realised
balance and within-centre overlap (§9); that inference is conditional on these centres (§10); the
absent confounders with an E-value (§6); and the entry-conditioning selection effect (§2).

Wherever a §14 estimate appears beside the primary, state that they target different populations and
are not a like-for-like comparison, and label the §14 estimates as model-based transported results.
