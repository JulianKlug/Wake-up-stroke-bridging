# Stage 13 spec — the feasible-policy contrast

Implements roadmap Stage 13 [§14b]. Section references in brackets are to
`statistical_analysis_plan.md`. Numbers and decisions referenced as DECISION *n* are established in
Stage 0 and recorded in `../out/stage0_data_inventory.md`. Constants named in `SMALL_CAPS` are
Stage 1's and live in `config.py`. The classified frame this stage draws its population from is
`stage4_eligibility_classification.md` §11; the two [§3] restrictions it applies **neither** of are
`stage5_cohort_construction.md` §3; the ordinal fitter is `stage8_primary_outcome_estimator.md` §5;
the category prediction, the re-expansion rule, the two namings and the failure taxonomy it inherits
are `stage12_all_centre_standardisation.md` §6, §6.3, §7.3 and §14 — whose §19 is written as a
handover *to this document*, names four things Stage 13 must know, and is answered in §0.3.

**Precondition.** Written against the **landed Stages 1–12**. Stage 13 depends on Stage 12 in exactly
four places, each named where it occurs: `standardise.gcompute` — the regime seam this document
makes public, returning the new arithmetic-only `GComputation` record (§6.1); `standardise.record_removal`
— a rename and one message edit (§4.2); `standardise.diagnostics` — a rename (§10.2); and
`model.ordinal_probabilities` (§6). **It reaches into no private name of any module** (§0.1, §16.15). It
fits no propensity model and reads no Stage 6–12 result, for the same reason Stage 12 did not.

**Status.** Written 2026-08-27 against the landed Stages 1–12; **revised 2026-08-29 after the
engineering review** — nine findings, all folded, listed in §22.1. **Nine probe groups were run before any
section was drafted (§22), and three of them changed what this document says:**

- **§7.3 — the policy contrast is the eligible-population contrast diluted by the eligible share, to
  machine precision, and this is an identity rather than a finding.** A contraindicated patient is
  assigned EVT alone under *both* regimes, so their per-patient contrast is exactly 0 — measured
  `0.0`, not `1e-17` — and `RD_k^policy = (n_eligible / N) · RD_k^eligible` holds to **8.9e-17** at the
  point estimate and to 5.3e-16 across every one of 2000 replicates. Any implementation that computes
  the two sides independently must assert them equal (U7); any report that prints `RD_k^policy`
  without the factor prints a number whose size is mostly a proportion.
- **§9 — the contraindication indicator separates in 1 of 2000 replicates, its coefficient reaches
  −23.6 past `POLR_MAX_ABS_BETA`, the fit converges in 21 iterations — and it changes no reported
  quantity.** Ten of the nineteen contraindicated patients are mRS 6, so a resample in which *every*
  drawn contraindicated patient is mRS 6 exists and occurred once. In that draw the eligible
  contrast under the policy fit equals an eligible-only fit's **exactly** (gap 0.0): a separated
  indicator absorbs the patients it indexes and leaves `α_k`, `β` and `γ` to the others. So U8 binds
  the treatment coefficient alone, as Stage 12's T12 does, and a bound on the worst coefficient would
  drop this replicate for a reason no reported quantity cares about.
- **§5.1 — the fit is not [§14a]'s, and the difference is measured rather than assumed.** Pooling
  the nineteen contraindicated patients moves the shared coefficients: the eligible-subpopulation
  contrast under this fit differs from Stage 12's by 2.1e-3 at the point estimate and by a median
  1.1e-2 (max 9.7e-2) across replicates. Same sign at every threshold. That gap is what DECISION 10
  buys and costs, and §18 hands the question of printing it to Stage 14.

**One decision was taken for this stage, and it is the PI's.** [§14b] says *"same standardisation
machinery, all centres, all patients"* and does not say what the model is fitted on. **DECISION 10
(PI, 2026-08-27):** one proportional-odds model over the whole [§14b] population, with
`X = STANDARDISATION_COVARIATES` plus a contraindication indicator (§5.1). Taken before any [§14b]
estimate was examined, recorded in `../out/stage0_data_inventory.md`, and written into [§14b] as a
dated amendment in the same commit as this document. Of the three readings put to the PI, one was
concordant with the plan, one strained and one discordant; §5.1 records all three. **A fourth question
was put to the PI on 2026-08-29 — whether the contraindicated patients' distributions are the model's
prediction or their observed outcome — and the PI confirmed the model's prediction, which is what
[§14b] as written prescribes (§18). DECISION 10 stands unamended; nothing is held.**

Every number below was produced by running code; none is carried — where a Stage 12 number is quoted
for comparison it is labelled as theirs. It is the **sole source for the Stage 13 implementation**:
everything the implementer needs is here, and anything not here is not to be invented.

**Goal.** [§14b] whole. One pooled proportional-odds model over **all** classified patients at all
four centres — eligible and contraindicated — with the [§14a] covariates plus a contraindication
indicator; g-computation to two regimes, the **active** one assigning bridging to every eligible
patient and EVT alone to every contraindicated one, the **comparator** assigning EVT alone to
everyone; the two standardised mRS distributions and the quantities [§14a]'s Report sentence takes
from them, under [§14b]'s label; and percentile intervals from a patient-level bootstrap stratified
by centre. **Labelled an operational policy contrast** — what adopting a bridging policy would have
delivered in this cohort — never [§14a]'s ATE, never [§7]'s ATO, and never the biological effect of
IVT.

**Not in scope.** Any [§14a] arm — the support check, the treated-support restriction, the random
centre intercept — for reasons §19 gives by name; the reporting layer and every label it owes
(Stage 14 [§16]); and any recomputation of a Stage 6–12 quantity.

---

## 0. Where Stage 13 sits

```
  data.load()                     ->  (df, audit)      25 cols, 126 rows   [Stage 2 §11]
  derive.derive(df, audit)        ->  df               31 cols             [Stage 3 §11]
  eligibility.classify(df, audit) ->  df               32 cols, 126 rows   [Stage 4 §11]
        │
        ├──────────────────────────────>  cohort.build(...)              93 rows   Stages 6-11, [§7]
        ├──────────────────────────────>  standardise.population(...)   104 rows   Stage 12,   [§14a]
        ▼   THE UNRESTRICTED FRAME.  NEITHER [§3] restriction is applied            [§14b]
  +--------------------------------------------------------------------------------------+
  |  policy.py                                                                            |
  |                                                                                       |
  |  population(df, audit)            -> DataFrame   126 -> 123                  (§4.1)  |
  |     [§11] covariate-complete on STANDARDISATION_COVARIATES  &  outcome observed        |
  |     eligibility.retained is READ and never applied as a filter               (§4.1)  |
  |     + the CONTRAINDICATED indicator, derived from `eligibility` alone        (§5.2)  |
  |                                                                                       |
  |  contrast(pop, audit)             -> Policy                                  (§13)   |
  |     X, dropped = model.design(pop, (TREATMENT,) + POLICY_COVARIATES)         (§5.3)  |
  |     U1  the exposure column must survive                                     (§5.3)  |
  |     fit = model.polr(X, y)              UNWEIGHTED, one fit, no propensity   (§5.1)  |
  |     active     = X with TREATMENT := 1[eligible]   a VECTOR, not a scalar   (§6.2)  |
  |     comparator = X with TREATMENT := 0                                       (§6.2)  |
  |     U2  no contraindicated row carries 1 in the active design, BEFORE predict (§6.3) |
  |     standardise.gcompute(...)      Stage 12's g-computation, regimes passed in (§6.1) |
  |     rd_k, mrs_0_2, mortality;  eligible_rd_k;  share_eligible;  U7 identity  (§7)    |
  |                                                                                       |
  |  inference(pop, audit)            -> bootstrap.Bootstrap                     (§10)   |
  |     bootstrap.replicates(pop, body, N_BOOT, SEED, BOOT_STRATUM)  ALL FOUR STRATA       |
  |     ONE arm per draw; 30 estimand keys; NO p on ANY key                      (§10.3) |
  +--------------------------------------------------------------------------------------+

  standardise.py gains THREE public names and ONE record, all lifted from Stage 12 privates:
      gcompute        (§6.1)   the g-computation, regimes as per-row vectors -> GComputation
      GComputation    (§6.1)   the arithmetic-only record gcompute returns; _standardise wraps it
      record_removal  (§4.2)   the exact-naming removal ledger, message made caller-neutral
      diagnostics     (§10.2)  the generic Replicate -> Diagnostics tally
  config.py gains a column name, a computed covariate view and two FAILURE_BUCKETS tokens  (§15)
  eight `C.SchemaError` / `FitError` identifiers, and they are the U series             (§11)
  U1 and U8 are policy.py's OWN literal raise sites — the raise-site scan reads tokens    (§11)

                  ->  Policy, Bootstrap                                                (§3.1)

  audit: the Stage 2-4 path's 12 entries, then FIVE of this stage's                    (§12)

  every standardised probability, every risk difference, every interval limit, the
  treatment coefficient and the contraindication coefficient on the workbook are
  DELIBERATELY ABSENT from this document. They are in the gitignored log.            (§4.3)
```

### 0.1 Why this is a second module and not a fourth function in `standardise.py`

`policy.py`, new. Stage 12 §0.1 declined a machinery/analysis split and named the cost: *"`standardise.py`
will be imported by Stage 13 for its regime seam, so a Stage 13 change can reach [§14a]'s code. §19
fixes what may and may not move."* This document takes the other side of that boundary and keeps it
narrow.

**Why not a `policy()` function inside `standardise.py`.** [§14] is one section with two subsections,
but the two make different statements over different populations under different labels, and [§14b]
is *"the only §14 analysis whose population includes contraindicated patients"*. A module whose
docstring opens with *"the estimand is an ATE in the all-centre eligible population"* (`standardise.py:11`)
cannot also house a function whose population is not eligible-only without one of the two sentences
becoming false. Stage 12 §4.2 also measured that `Audit.entry` is first-match and that two row-removal
ledgers under `kind="cohort"` from one frame collide silently; a third population inside the same
module would be a third such ledger behind one name.

**What crosses the boundary, and it is three names and one record.** `standardise.gcompute` — the
arithmetic half of Stage 12's `_standardise`, with its two regimes passed in as per-row vectors rather
than fixed inside, returning the new `GComputation` record (§6.1); `standardise.record_removal` (§4.2);
and `standardise.diagnostics` (§10.2). `_standardise` **stays private** and becomes a wrapper: it calls
`gcompute` with scalar arms and adds `beta`, `measure` and the population label to make a
`Standardisation`, so `all_centre` and `hierarchical` do not change a line, and Stage 12's sixty-nine
intervals and ten audit entries are asserted byte-identical across the change by two digests captured
before it (§16.11). Nothing else in `standardise.py` is touched, and `test_standardise.py`'s
private-cross-module-call assertion (`private == {"balance._pooled_sd"}`) is extended to `policy.py`:
the set of private names it reaches into is **empty** (§16.15). That is why U1 and U8 are `policy.py`'s
own raise sites and not calls to Stage 12's T1 and T12 guards (§11).

### 0.2 What Stage 13 does not touch

**It fits no propensity model** — [§14]'s own first sentence — and imports neither `propensity.py`
nor `cohort.py`. Stage 12 needed `cohort.treating_centres` for the never-IVT complement; this stage
has no support check (§19) and no reason to know which centres treat.

**It reads no Stage 6–12 result.** It takes the classified frame and nothing else. §16.11's guard is
that the Stage 1–12 pipeline run with and without this stage produces byte-identical numbers and an
identical 56-entry audit prefix.

**It calls neither `cohort.build` nor `standardise.population`.** The first applies both [§3]
restrictions; the second applies restriction 2. [§14b]'s population is the one *neither* removes
anyone from, and Stage 12 §22 declined to give `population` an eligibility keyword speculatively.
§4.1 is why this document does not add one either.

**It does not read `ivt_contraindicated`.** The contraindication indicator is derived from the
`eligibility` column and from nothing else (§5.2). Stage 4 is the one classifier; a second reading of
the flag would be a second classifier that agrees with the first today.

### 0.3 Stage 12 §19's four items, answered

| §19 item | Answer |
|---|---|
| 1. `population()` must not be reused unchanged; the clean seam is a keyword Stage 12 did not add | **No keyword is added.** `policy.population` applies the [§11] masks over the unrestricted frame and reads `eligibility` without filtering on it (§4.1). A keyword whose two values select two populations under one step name is Stage 12 §4.2's collision with a switch on it |
| 2. The regime is a vector written into the treatment column of the fitted design | **Yes**, and it is the whole change: `gcompute` takes `active` and `comparator` as per-row arrays; `all_centre` passes broadcast scalars (§6.1, §6.2) |
| 3. DECISION 1a made `eligibility` binary; do not look for a third value | The indicator is `eligibility == INELIGIBLE`; U4 asserts the column takes no value outside `ELIGIBILITY_ORDER` before it is read (§5.2) |
| 4. The Accept-when is one line on the active design's treatment column, before prediction | **U2**, a `SchemaError` (§6.3). Measured: 0 violations at the point estimate and across 2000 replicates, which is what a contract should measure |

**And the thing §19 said Stage 13 must re-measure — the cutpoint-collapse rate, the drop rate — is
re-measured**: **0.65%** and **0** (§6.4, §10.4). Both are lower than Stage 12's 4.15% and 0.1%,
and §6.4 says why the first is.

---

## 1. Deliverables

| Path | Contents |
|---|---|
| `extended_bridging/policy.py` | **new.** [§14b] whole: three public names (§3.2), one record (§3.1), thirty estimand keys (§3.3) |
| `extended_bridging/standardise.py` | **amended: one extraction, one new record, three renames, one message edit.** The arithmetic of `_standardise` moves into public `gcompute`, which takes `active`/`comparator` per-row arguments and `keys=` and returns the new frozen `GComputation`; `_standardise` stays private as the `Standardisation`-building wrapper (§6.1). `_regime_design` accepts a vector. `_record_removal` → public `record_removal`, its message made caller-neutral (§4.2). `_diagnostics` → public `diagnostics` (§10.2). Public surface 5 → 8. `all_centre`, `hierarchical`, `support`, `population`, `inference` do not change, and §16.11 asserts it by two digests |
| `extended_bridging/config.py` | **amended.** `CONTRAINDICATED`, `POLICY_COVARIATES`, two `FAILURE_BUCKETS` tokens (§15, fence) |
| `extended_bridging/tests/test_policy.py` | **new.** §16 |
| `extended_bridging/tests/fixtures_stage13.py` | **new.** §16.0, the one place its contents are enumerated |
| `extended_bridging/tests/test_standardise.py`, `test_config.py`, `test_bootstrap.py` | **amended.** The surface list 5 → 8 and the privates list; the two pinned digests captured before the refactor (§16.11); the `GComputation` and vector-regime tests (§16.15); the config assertions (§16.2); the raise-site scan scope grows to `policy.py` and its per-module pin gains `policy.py: ["U1", "U8"]`, 29 sites (§11, §16.17) |
| `extended_bridging/statistical_analysis_plan.md` | **landed with this document in `37051fc`.** [§14b]'s dated amendment of 2026-08-27 stating the model, the two regimes, the reported quantities and the dilution identity (DECISION 10). Nothing further to write |
| `extended_bridging/implementation_roadmap.md` | **landed with this document in `37051fc`.** Stage 13's `**Spec:**` line, the model, the identity, the U series and the measured rates (§23). Nothing further to write |
| `TODOS.md` | **four items landed in `37051fc`** (§17 items 2, 3, 5, 7). **Two existing items amended by the review of 2026-08-29** (§17 items 9, 10): the row-builder item gains Stage 13's three renderers as instances; the `code`-field item gains the U1/U8 guard duplication as a recorded cost |
| `../out/stage0_data_inventory.md` | **amended, gitignored.** DECISION 10 (§5.1). The one document deliverable that cannot be verified from the repository |

`out/logs/audit_<label>.md` is an **output, not a deliverable**, and for this stage it is the only
place any standardised probability, any risk difference, any interval limit, the treatment
coefficient or the contraindication coefficient exists.

**Nothing under `specs/` may quote a case identifier, and nothing here does.** Stage 12 §4.3 extended
the rule to standardised probabilities and centre intercepts; this document extends it to **the
contraindication coefficient**, which is a prognosis contrast on nineteen named patients' worth of
records and is not a result of this analysis (§4.3).

---

## 2. Environment

Unchanged from Stage 12: `uv`, Python 3.12.12, pandas 2.3.3, numpy 1.26.4, statsmodels 0.14.6.
**No dependency is added and no new estimator is written.** Every fit is `model.polr`; every
prediction is `model.ordinal_probabilities`.

```bash
cd extended_bridging && uv sync && uv run pytest -v
```

**Cost, measured (§22, P6).** One `bootstrap.replicates` call, one arm per draw, on 123 records with
an eleven-column design: **27.6 s for 2000 replicates**. That is about a ninth of Stage 12's bootstrap
(4.2 min measured at §28 there), because there is no `polr_ri` here, and it does not move the
parallelisation `TODOS` item's trigger in either direction.

---

## 3. Module shape

### 3.1 What the stage returns

**One record of this stage's own — `Policy`.** `inference` returns `bootstrap.Bootstrap`, which this
stage does not define. Frozen, carrying no verdict, on Stage 8 §3's rule.

```python
@dataclass(frozen=True)
class Policy:
    """The [§14b] feasible-policy contrast: one fit, one averaging population, two REGIMES."""

    n_average: int                       # the averaging population's size — §7.1, [§11]
    n_fit: int                           # equal to n_average; carried for Stage 14's denominators
    n_eligible: int                      # patients the active regime bridges — §7.3
    share_eligible: float                # n_eligible / n_average, the dilution factor — §7.3
    distribution: dict[str, dict[int, float]]   # "active" | "comparator" -> mRS level -> P
    cumulative: dict[str, dict[int, float]]     # regime -> threshold -> P(Y <= k)
    rd: dict[int, float]                 # RD_k, keyed by MRS_THRESHOLDS — §7.2
    mrs_0_2: float                       # == rd[2]; [§14a]'s name for it, under [§14b]'s label
    mortality: float                     # == -rd[5]; likewise
    eligible_rd: dict[int, float]        # the UNDILUTED eligible-subpopulation contrast — §7.3
    conditional_log_odds: float          # beta, conditional on X AND on contraindication — §8
    conditional_odds_ratio: float        # exp(beta). CONDITIONAL. — §8
    contraindication_log_odds: float     # the indicator's coefficient; a NUISANCE, log only — §9
    measure: str                         # the guard string, always _CONDITIONAL_14B — §8
    label: str                           # the guard string, always _OPERATIONAL — §8
    fit: model.PolrFit
    dropped: tuple[str, ...]             # design columns dropped as constant — §5.3
```

**`distribution` and `cumulative` are keyed by regime NAME and not by arm code**, where Stage 12's are
keyed by `C.TREATMENT_LABELS`' codes. That is not inconsistency, it is the difference between the two
stages: the active regime is **not an arm**. Under it 104 patients receive bridging and 19 receive
EVT alone, and a dict keyed `1` would say every patient was bridged. `"active"` and `"comparator"`
are fixed literals (`_ACTIVE`, `_COMPARATOR`) and are the estimand-key stems (§3.3).

**Two guard strings, not one.** Stage 12 §8 carries `measure` so that `exp(β)` is never labelled
marginal. This stage has a second thing that must never be mislabelled — the *estimand* — and it
carries `label` for it (§8). Both are `Final` module strings and a second value can only arrive with
a [§14] amendment.

**There is no field called `odds_ratio`**, for Stage 12 §8's reason, inherited unchanged. And there
is no field called `rd_eligible_share_adjusted` or any other pre-multiplied quantity: `rd` **is** the
policy contrast, `eligible_rd` is what it dilutes, `share_eligible` is the factor, and U7 asserts the
three agree. Three fields that encode one identity are kept because two of them are what [§14b] and
[§14a] each report and the third is what makes the relationship readable (§7.3).

### 3.2 The public names

`policy.py` exposes three: `population`, `contrast`, `inference`. **Its privates are enumerated here
and `test_policy.py` asserts the list**, on `test_standardise.py:991`'s pattern: `_regimes` (§5.2),
`_assert_exposure_survived` (U1, §5.3), `_assert_beta_reportable` (U8, §9.2), `_assert_dilution` (U7,
§7.3), `_estimand_keys` (§3.3), `_values` (the thirty draws of one `Policy`), `_replicate` (§10.2),
`_fit_table`, `_contrast_table`, `_replicates_table` (§12), **and three the implementation added
(2026-08-29)**: `_assert_regime` (U2, called once per `(design, active)` pair, §6.3),
`_population_table` (entry 1's per-centre ledger, §12) and `_estimate` (the fit-to-`Policy` path
`contrast` and the body share, so the two cannot drift; U8 stays the caller's, §10.2). Nothing else.

`standardise.py` goes from five public names to **eight** (`gcompute`, `record_removal`, `diagnostics`)
and gains one public record, `GComputation` (§6.1); its private list loses `_diagnostics` and keeps
`_standardise`. `model.py`, `bootstrap.py`, `balance.py`, `cohort.py`, `eligibility.py` are untouched.

### 3.3 The estimand keys — thirty, and the split is [§14b]'s own read through [§14a]'s Report sentence

```
   rd_0 … rd_5                          len(C.MRS_THRESHOLDS)                        6
   mrs_0_2, mortality                   the two namings — §7.2                       2
   active_0 … active_6                  len(C.MRS_LEVELS)                            7
   comparator_0 … comparator_6          len(C.MRS_LEVELS)                            7
                                                                                   ───
                          [§14a]'s Report sentence, under [§14b]'s label            22

   eligible_rd_0 … eligible_rd_5        the undiluted contrast — §7.3                 6
   share_eligible                       the dilution factor — §7.3                    1
   beta                                 one conditional log-odds, one fit — §8        1
                                                                                   ───
                                                                                     30
```

**The first twenty-two are [§14b]'s *"same standardisation machinery"* read literally**: the machinery
reports two standardised distributions and, from them, `RD_k`, the mRS 0–2 difference and the mortality
difference ([§14a] Report). `beta` is the *"may appear as a model parameter"* permission, one per fit,
and there is one fit.

**The other seven are the decomposition, and the argument for them is [§14b]'s and not an
engineering one** — this is the check Stage 12 §21 item 14 asks for, run on this decision before it
was made. [§14b] says the contrast's *"interpretation is operational — what adopting a bridging policy
would deliver in this cohort"*. §7.3 measures that `rd_k` is `share_eligible · eligible_rd_k` to
machine precision. A reader given `rd_k` alone therefore reads a number that is 84.6% proportion and
the rest effect, with no way to tell which moved; the same reader given the three can. And
`share_eligible` is not a constant of the design: it is a **statistic** of the resampled cohort —
measured 0.748 to 0.935 across replicates — so its interval is part of what "in this cohort" means
under [§14]'s bootstrap. **If the PI reads [§14b]'s sentence more narrowly, the seven keys move to the
log and the count becomes twenty-three**; §18 records that as the one PI-reversible choice in the key
set.

**There is no `contraindication` key.** The indicator's coefficient is a nuisance parameter (§9); it
is printed in the fit entry with its role and gets no interval, because an interval invites a reading
as a prognostic effect of contraindication that this analysis is not designed to estimate (§4.3).

**The counts are computed from `C.MRS_LEVELS` and `C.MRS_THRESHOLDS` and never written as 6, 7 or 22**;
§16.3 asserts `len(keys) == 2*len(MRS_THRESHOLDS) + 2 + 2*len(MRS_LEVELS) + 2` against the config.

**The key NAMES `active_*` and `comparator_*` are fixed literal strings** — the wire format Stage 14
matches on — and are deliberately **not** Stage 12's `dist1_*`/`dist0_*`, because a Stage 14 formatter
that matched `dist1_` across both stages would print the policy regime's distribution under an arm
heading (§3.1).

**No key carries a p-value.** [§14]'s Inference paragraph prescribes percentile intervals and nothing
else, which is Stage 12 §3.3's argument unchanged, and `bootstrap.intervals` is passed
`lambda key: False`.

---

## 4. The population [§14b, §11] — and it is nobody else's

### 4.1 The ledger, measured (P1)

```
   the classified frame                                                126        [Stage 4 §11]
   ├─ [§3] restriction 1 NOT applied.  USZ stays.                                 [§14b]
   ├─ [§3] restriction 2 NOT applied.  the 19 contraindicated stay.               [§14b]
   ├─ [§11] covariate-complete on STANDARDISATION_COVARIATES           125    -1
   └─ [§11] the primary outcome observed                               123    -2
```

**The three [§11] losses are the same three records Stage 12 §4.1 names** — one covariate-incomplete
on `core_ml` and `tmax6_ml`, two outcome-missing at USZ — and **all three are eligible**. Every one of
the nineteen contraindicated patients is complete on the covariates and has an observed outcome. So
the [§14b] population is Stage 12's 104 plus exactly the 19 the [§3] flag removes, and §16.1 asserts
that by set difference rather than by count.

| centre | eligible, EVT alone | eligible, bridging | contraindicated, EVT alone | n |
|---|---|---|---|---|
| HUG | 11 | 30 | 11 | 52 |
| CHUV | 14 | 7 | 0 | 21 |
| Lugano | 28 | 2 | 0 | 30 |
| USZ | 12 | 0 | 8 | 20 |
| **all** | **65** | **39** | **19** | **123** |

**Two facts about that table govern the rest of this document.** First, **no contraindicated patient
is bridged** — the third column has no bridging partner, because Stage 4's load-bearing assertion
(*"no treated patient carries a 1"*) already holds on the classified frame, and U3 re-checks it at
this stage's boundary rather than trusting it (§11). Second, **the contraindicated patients are at two
centres only**, HUG and USZ — a stratified resample therefore draws its contraindicated count from
two strata of 52 and 20, which is why the count runs 8 to 31 across replicates (§10.4) and why the
indicator is never constant on a draw.

**The mask is spelled here and is `complete & observed` with no eligibility conjunct.** Stage 12 §4.1
declined a shared [§11] helper because two callers agreed on nothing but the `&`; this stage is a
third caller whose mask differs from *both* in content — no `in_model`, no `retained` — and the
argument holds a third time.

**`eligibility.retained` is READ and never applied.** It supplies the regime (§6.2) and the indicator
(§5.2). A frame on which it filtered anyone would be [§14a]'s.

### 4.2 The removal ledger, and Stage 12's copy becomes public

Stage 12 §4.1 wrote its own `_record_removal` because `cohort._record_removal` is private and
`cohort.py` was not to change. **That function becomes public as `standardise.record_removal`**, a
rename **and one message edit**: its `SchemaError` text reads *"A [§14a] population that cannot name
what it removed…"* (`standardise.py:190-194`) and would describe a [§14b] removal as [§14a]'s. The
message becomes `"{step} removed N record(s) and names M. A population that cannot name what it removed
makes the analysis unreconstructable from the log."`, and the docstring's *"it is a second copy on
purpose"* paragraph — now false — is replaced by one sentence saying it is the shared public rule with
two callers. So this stage calls it rather than writing a third copy — the move Stage 12 §10.3 made for
`balance.levels`, on the same argument: the rule it carries (an entry must name exactly as many
patients as it says it removed) may not exist twice. §16.15 asserts the message names no section. `TODOS` gains an item: there are now two
implementations, `cohort`'s private and `standardise`'s public, and the third caller is the trigger
for one of them to move to `data.py` (§17).

**Every Stage 13 step name is distinct by construction and asserted disjoint** from `cohort.build`'s
six (Stage 12 §4.2 enumerates them) and from `standardise`'s ten (Stage 12 §15). Three row-removal
ledgers under `kind="cohort"` from one classified frame — 93, 104 and 123 rows — in a Stage 14 driver
is exactly the first-match hazard Stage 12 §4.2 describes, with a third party (§12).

### 4.3 What may not be quoted, and this stage adds one thing

Stage 8 §4.3 forbade case identifiers; Stage 10 §4.5 interval limits and p-values; Stage 11 §1 E-values
and adjusted p-values; Stage 12 §4.3 standardised probabilities and centre intercepts. **Stage 13 adds
the contraindication coefficient.**

It is a nuisance parameter and would be read as a finding: *"having an absolute contraindication to
IVT shifts the 90-day mRS odds by this much"*. Nineteen patients at two centres cannot support that
sentence, this analysis does not claim it, and [§14b] does not ask for it. What this document quotes
is its **sign** (negative on the `P(Y ≤ k)` scale, the direction [§3] predicts), its rank (the largest
of the eleven coefficients in magnitude at the point estimate, by more than a factor of two over the
next), and its behaviour across replicates (§9). Its value is the log's.

---

## 5. The pooled proportional-odds model [§14b, DECISION 10]

### 5.1 The fit, and the decision behind it

    logit P(Y <= k | A, X, D)  =  alpha_k  +  beta*A  +  gamma'X  +  delta*D ,    k = 0 … 5

with `X = C.STANDARDISATION_COVARIATES` — the [§6] covariates minus `center`, exactly [§14a]'s — and
`D = 1[contraindicated]`. Fitted **once, unweighted, on all 123 records**, with `model.polr` unchanged.

**[§14b] left the model open and three readings were concordant with the plan or nearly so.** They
were put to the PI on 2026-08-27, with the following ledger, and the first was taken:

| Reading | Concordance with the plan | Decision |
|---|---|---|
| **All patients, `X` + a contraindication indicator** | Honours *"all patients"*. Adjusts the bias [§3] restriction 2 names — *"retaining them makes 'no IVT' a marker of contraindications and their prognosis"*. Respects [§14a]'s exclusion on its own terms, which is **conditional**: *"contraindication status is not a covariate, **the cohort being already restricted** to eligible patients"* — a clause whose premise [§14b] lifts | **DECISION 10** |
| All patients, `X` unchanged | The most literal *"same machinery"*. **Discordant with [§3]**: it builds, inside the one analysis that lifts restriction 2, the confounding restriction 2 was written to remove — 19 EVT-alone patients with worse prognosis enter `α_k` and `γ` and bias `β` | declined |
| Reuse the [§14a] fit, average over all patients | Concordant with *"same machinery"* and with the IVT prohibition. Strained on *"all patients"*: the population is 123 and the model saw 104, and the contraindicated patients' EVT level is extrapolated from people unlike them. Makes [§14b] a re-weighting of [§14a] | declined |

**What DECISION 10 changes and what it does not, measured (P5).** The nineteen extra records move the
*shared* coefficients — `|γ|` differs from Stage 12's fit by up to 0.29 on the workbook and by a
median 0.36 (max 1.21) on a common draw across replicates — and therefore move the eligible patients'
counterfactuals: the eligible-subpopulation contrast under this fit differs from Stage 12's `all_centre`
by **2.1e-3** at the point estimate and by a median 1.1e-2 (max 9.7e-2) across 2000 replicates, with
the **same sign at every threshold** on the workbook. So [§14b]'s eligible contrast is not [§14a]'s
and is not reported as though it were; §18 hands Stage 14 the question of whether to print the gap.
What DECISION 10 does *not* change: `β` is identified from eligible patients only, because `A` is
constant among the contraindicated — which is why the indicator can separate without touching it (§9).

**Centre is omitted for [§14a]'s reason, and [§15]'s permission covers [§14] whole**: *"Omitting
centre and modelling across centres are permitted inside §14 only."* `test_config.py` holds
`POLICY_COVARIATES` to `STANDARDISATION_COVARIATES + (CONTRAINDICATED,)` statically (§16.2), so the
permission cannot widen by editing a tuple, and `CONTRAINDICATED ∉ PS_COVARIATES` is asserted beside
it — a contraindication indicator in a [§7] propensity model would be [§3] restriction 2 undone.

**Linear terms only, no interaction, which is [§6]'s rule and is a stated assumption here.** `D` shifts
the cutpoints and nothing else: the model assumes a contraindicated patient's covariates act on their
outcome as an eligible patient's do. That is the same class of assumption [§14a] makes about centre —
untestable here for the same reason, nineteen records at two centres — and §17 carries it.

### 5.2 The indicator is derived from `eligibility`, and from nothing else — in ONE private helper

```python
def _regimes(frame: pd.DataFrame) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """(D, active, comparator) for `frame`, from `eligibility.retained` and nothing else. U4.

    THE ONE PLACE the contraindication predicate is spelled. `population` stores D as the
    CONTRAINDICATED column; `contrast` and the replicate body call this on THEIR frame and assert the
    stored column agrees — so a frame whose column and eligibility disagree raises rather than fitting
    one thing and standardising another. `D + active/treated == 1` on every row by construction.
    """
    if C.ELIGIBILITY not in frame.columns or not frame[C.ELIGIBILITY].isin(C.ELIGIBILITY_ORDER).all():
        raise C.SchemaError("U4  ...")
    retained = eligibility.retained(frame).to_numpy()
    treated, control = max(C.TREATMENT_LABELS), min(C.TREATMENT_LABELS)
    D = (~retained).astype(float)
    active = retained.astype(float) * treated + (~retained).astype(float) * control
    comparator = np.full(len(frame), float(control))
    return D, active, comparator
```

`population` writes `pop[C.CONTRAINDICATED] = D`; `contrast` and the body call `_regimes` on the frame
they are given and raise `SchemaError` (under U4) if `frame[C.CONTRAINDICATED]` differs from the `D`
they compute. The review of 2026-08-29 found the first draft derived this predicate three times —
`== INELIGIBLE` in `population`, again in the body, and `eligibility.retained` for the regime — which
is the "two things that can disagree after an edit" §7.3 argues against, applied to a definition.

**Not from `ivt_contraindicated`.** Stage 4 is the one classifier — DECISION 1a made the flag *its*
input — and a stage that re-read the flag would be a second classifier that agrees with the first on
this workbook and could stop agreeing on the next without anything raising. `eligibility` is a string
column in `ELIGIBILITY_ORDER`, and **U4 asserts every value is in that tuple before the comparison**,
because `== INELIGIBLE` on a column carrying a third value would silently code it as *not
contraindicated* (§11).

**The column is added to the population frame by `policy.population` from `_regimes` and is a `float` 0/1**, which is
how `model.design` receives `atrial_fib` and `sex`: it is not in `C.CATEGORICAL`, so `design`
neither dummies it nor drops a reference level, and it enters as one column. Measured (P2): the
design is **11 columns, rank 11, nothing dropped**. `_assert_design_inputs` accepts it — measured
rather than assumed, because that guard checks factor values and missingness and a float indicator
with no `nan` passes both.

### 5.3 The design, and U1

    X, dropped = model.design(pop, (C.TREATMENT,) + C.POLICY_COVARIATES)

**U1 is Stage 12's T1 over this design, and it is `policy.py`'s OWN raise site**: `FitError` with the
literal leading token `U1` when `C.TREATMENT` is not in `X.columns`, a dropped covariate recorded in
`Policy.dropped` and printed, never raised on. It is a three-line private, `_assert_exposure_survived`,
and **not** a call to `standardise._assert_exposure_survived` with a token argument — because
`test_bootstrap.py`'s raise-site scan (`raise_sites()`, `test_bootstrap.py:788-793`) reads the first
token of every `raise FitError(` **statically** and asserts it is a string literal; a parameterised
token is an `ast.FormattedValue`, the scan's own assertion fails, and the whole suite goes red. The
scan is the repository's one defence against the wrong-but-mapped bucket (T4's story, Stage 12 §14),
and a shared guard would weaken it to save three lines. `TODOS`' `code`-field item records the cost
(§17 item 9). **The indicator falls under the covariate rule and not the exposure rule**,
deliberately: a draw with no contraindicated patient drops `D` as constant, the model on that draw *is*
[§14a]'s, the active regime bridges everyone, `share_eligible` is 1 and `eligible_rd_k == rd_k` — a
legitimate draw of a cohort with no contraindicated patient, and the record says so through `dropped`
and `share_eligible`. **§16.16 asserts exactly that on `no_contraindicated_frame()`** rather than
leaving it as prose. It is also **unreachable on the workbook**: measured over 2000 stratified
replicates the contraindicated count is never below **8** (§10.4), because 11 of HUG's 52 and 8 of
USZ's 20 are contraindicated and both strata would have to miss all of theirs at once.

**Measured on the estimation population (P2)**: converged in **5 iterations on the likelihood
criterion**, 0 rescales, 0 halvings, first step norm 8.2, 6 cutpoints over all seven levels,
17 parameters on 123 records. Stage 12's fit is 5 iterations on 104 records with 16 parameters.

---

## 6. The regimes [§14b]

### 6.1 The seam: `standardise.gcompute`, and the two regimes are per-row vectors

Stage 12's `_standardise(X, fit, over, label, dropped, groups)` loops `for arm in (treated, control)`,
overwrites the treatment column with a scalar, and builds a `Standardisation` — a record whose
`population` field takes one of three [§14a] labels (`test_standardise.py:1345-1347`), whose `measure`
is the `_CONDITIONAL` string *"conditional odds ratio [§14a] — conditional on X"*, and whose
`distribution` is typed `dict[int, ...]` by arm code. **None of those three is true of a [§14b]
g-computation**, and the review of 2026-08-29 declined to have Stage 13 receive a record carrying a
label about a fit it was not made for, even transiently: §8's rule is that a guard string is never
false. So the seam is cut one level lower.

**The arithmetic becomes public as `gcompute`, returning a new arithmetic-only record:**

```python
@dataclass(frozen=True)
class GComputation:
    """What a g-computation COMPUTES, and nothing about what it means. Stage 12 §7 / Stage 13 §6.1.

    Two averaged distributions keyed by whatever `gcompute` was told to call them, the cumulative
    tables, RD_k and the two namings. NO `beta`, NO `measure`, NO population label: those belong to
    the record that interprets this one — `Standardisation` for [§14a], `Policy` for [§14b] — and a
    record that carried both would be a record one of its two callers has to lie in.
    """

    n_average: int                              # over.sum()
    n_fit: int                                  # len(X)
    distribution: dict[object, dict[int, float]]   # key -> mRS level -> P, over MRS_LEVELS
    cumulative: dict[object, dict[int, float]]     # key -> threshold -> P(Y <= k)
    rd: dict[int, float]                        # cumulative[keys[0]] - cumulative[keys[1]]
    mrs_0_2: float                              # == rd[2], one computation
    mortality: float                            # from the distributions; asserted == -rd[5]


def gcompute(X: pd.DataFrame, fit: model.PolrFit | model.RIFit, over: np.ndarray,
             active: np.ndarray, comparator: np.ndarray, *, label: str,
             groups: pd.Series | None = None,
             keys: tuple[object, object] = (1, 0)) -> GComputation:
    """Overwrite the treatment column of the FITTED design with `active`, then with `comparator`,
    predict, re-expand onto MRS_LEVELS, average over `over`, run T10 and T11 under `label`.
    `active` and `comparator` are per-row arrays of length len(X) — a length mismatch raises
    SchemaError. `_standardise` passes `np.full(len(X), treated)` and `np.full(len(X), control)`;
    Stage 13 passes `1[eligible]` and zeros. `keys` names the two distributions: Stage 12's arm codes,
    Stage 13's regime literals. `keys[0]` is the minuend of every RD_k."""


def _standardise(X, fit, over, label, dropped, groups=None) -> Standardisation:
    """UNCHANGED SIGNATURE, UNCHANGED CALLERS. Calls `gcompute` with the two arm codes broadcast to
    vectors and wraps the result with `beta`, `_CONDITIONAL` and the [§14a] label."""
```

Everything after the overwrite — `_arm_probabilities`, `_expanded`, the averaging, T10, T11, the two
namings — is **unchanged and shared**, which is what Stage 12 §0.1 said the reusable half was.
`all_centre`, `hierarchical`, `support`, `population` and `inference` **do not change a line**: they
call `_standardise` as before. §16.11 asserts Stage 12's sixty-nine intervals and ten audit entries are
byte-identical across the change against digests captured before it; a scalar broadcast to a vector
is the same float in every row. §16.15 asserts `gcompute` with broadcast scalars equals `_standardise`'s
arithmetic field for field, and that a mismatched regime length raises.

**`_regime_design(X, arm)` accepts the vector**: `Xa[C.TREATMENT] = np.asarray(arm, dtype=float)` with
a length check. A scalar still broadcasts. Nothing in `ordinal_probabilities` sees the difference —
its by-name column check passes because the columns are the fitted ones, which is Stage 12 §7.1's
whole argument for overwriting rather than re-designing, and it holds a fortiori for a vector.

### 6.2 The two regimes, written out

    active     :  A_i = 1  if eligibility_i == ELIGIBLE      bridge every eligible patient
                  A_i = 0  if eligibility_i == INELIGIBLE    direct EVT for the contraindicated
    comparator :  A_i = 0  for every i                        direct EVT for everyone

The active vector is `_regimes(pop)`'s second element (§5.2), built from `eligibility.retained` — **the
same predicate that defines [§14a]'s population**, so the two subsections of [§14] agree about who is
eligible by construction and not by coincidence, and the indicator `D` and the regime come from one
expression rather than two. The arm codes come from `C.TREATMENT_LABELS` as `outcome.py:118-119` takes
them; `1[eligible]` is written as `retained.astype(float) * treated + (~retained) * control`, never as
a literal.

**Under the active regime 104 patients are bridged and 19 are not (P3).** Under Stage 12's arm 1 all
104 are bridged. The two designs differ in exactly the nineteen rows the [§14b] population added.

### 6.3 U2 — the roadmap's Accept-when, as a contract checked before prediction

    stored = frame[C.CONTRAINDICATED].to_numpy() == 1          # the STORED column, not _regimes' D
    if (active[stored] != control).any(): raise C.SchemaError("U2 ...")

*"No contraindicated patient is ever assigned a predicted IVT outcome."* Under §6.1's design-overwrite
that is one line on the regime vector, and it runs **before `ordinal_probabilities` is called**, not
after the average — an assertion on the average could not tell a violated regime from a small effect.

**U2 runs on EVERY `(design, active)` pair `contrast` and the body hand to `gcompute`, and it compares
two sources.** The outside voice of 2026-08-29 made two points that shaped this. First, a check of
`_regimes`' `active` against `_regimes`' own `D` is `x[~x] != 0` — an expression checked against itself.
So U2 reads the **stored** `CONTRAINDICATED` column (written by `population`, carried through
`resample`) and the freshly built vector: two sources that can disagree after an edit. Second, **the
`eligible_rd` call must not hand the full design to `gcompute` with every row set to `treated`**: that
would have `ordinal_probabilities` compute an IVT outcome for the nineteen contraindicated rows and
`over` discard it afterwards — [§14b]'s first paragraph happening inside the module whose contract is
that it never happens. So the second call receives **`X[eligible]`** (§7.3 item 3), a design on which
no contraindicated row exists, and U2 over it is trivially satisfied because there is nothing to
violate. No prediction of an IVT outcome for a contraindicated patient exists anywhere in `policy.py`,
even transiently, and §16.6's spy asserts the second call's design has `n_eligible` rows. `SchemaError` and not `FitError`: the regime is a *definition*, so a violation is a
bug in the code that built the vector and never a property of a resample, and Stage 12 §14's rule is
that a contract break kills the bootstrap rather than being counted.

Measured: **0 violations** at the point estimate and in every one of 2000 replicates, which is the
only number a contract should produce and is recorded so that U2 is not later read as unexercised —
§16.14 fires it on a constructed vector.

### 6.4 The re-expansion rule is inherited, and the collapse rate is a sixth of Stage 12's

`polr` fits the occupied levels, `ordinal_probabilities` returns them, and the caller re-expands onto
`C.MRS_LEVELS` with a structural zero at any unoccupied level — Stage 12 §6.2–6.3, unchanged and
called through `gcompute`. **Measured over 2000 replicates (P6): 13 lose a level, 0.65%, always
mRS 5; in those 13, `RD_4 == RD_5` exactly.** Stage 12 measured 4.15% on 104 records.

The reason is in the population: mRS 5 carries **5 patients of 123** here against 3 of 104 there,
because **two of the nineteen contraindicated patients are mRS 5**. So [§14b]'s inclusion of the
contraindicated makes the primary outcome's thinnest level less thin, and the collapse — which
narrows `RD_5`'s sampling distribution by construction (Stage 12 §6.3) — happens a sixth as often.
Counted, and Stage 14 prints it beside the `RD_5` interval as it does Stage 12's (§14).

---

## 7. The standardisation and its outputs

### 7.1 The averaging, and its denominator is the whole population

    P_hat(Y = j | regime)  =  (1/N) * SUM_i  P(Y = j | A_i(regime), X_i, D_i)

over all `N = 123` records, unweighted, under each regime in turn. `n_average == n_fit` here — there
is no `over=` restriction in [§14b] (§19) — and both are carried so that Stage 14 prints [§11]'s
denominator from the record and not from an assumption that the two are equal.

### 7.2 The two namings, inherited

`rd_k = P(Y ≤ k | active) − P(Y ≤ k | comparator)`; `mrs_0_2 = rd[2]`, one computation; `mortality =
P(Y = 6 | active) − P(Y = 6 | comparator)`, a second computation asserted equal to `−rd[5]`. Stage 12
§7.3's identity discipline and tolerances, unchanged, run inside `gcompute` as T10 and T11 on every
record it builds. Measured (P3): `|mrs_0_2 − rd[2]|` is 0.0 and `|mortality + rd[5]|` is 1.1e-16.

The orientation is Stage 8 §7's: a positive `rd_k` favours the bridging policy on the functional
scale and a positive `mortality` is **worse**.

### 7.3 The dilution identity — the structural finding of this stage

For a contraindicated patient `A_i` is `control` under both regimes, so

    P(Y = j | active, X_i, D_i) − P(Y = j | comparator, X_i, D_i)  =  0      exactly, per row

— measured `0.0` on every one of the nineteen rows, not `1e-17`: the two predictions are the same
function of the same floats. Summing over the population,

    rd_k  =  (n_eligible / N) · eligible_rd_k ,     eligible_rd_k = mean_{eligible} [P(Y≤k|1,X_i,0) − P(Y≤k|0,X_i,0)]

**Measured: the identity holds to 8.9e-17 at the point estimate and to 5.3e-16 across every replicate
(P3, P6).** `share_eligible` is 104/123 = 0.8455 on the workbook and runs **0.748 to 0.935** across
stratified replicates.

**Three consequences, each specified.**

1. **The policy contrast is an eligible-population effect scaled by a proportion, and the report says
   so** — `eligible_rd_k` and `share_eligible` travel with `rd_k` as estimand keys (§3.3) and Stage 14
   prints the factor beside every `rd_k` (§14). Without them a reader compares [§14b]'s `rd_k` with
   [§14a]'s and sees a smaller effect where there is a smaller *share*.
2. **`eligible_rd_k` is not [§14a]'s `RD_k`**, and §5.1 measures by how much. It is the eligible
   contrast under the [§14b] fit. The two are printed under different labels and never in one column.
3. **U7 asserts the identity on every `Policy` the module constructs, at 1e-12** — T10's tolerance
   and T10's argument. The first draft said 1e-14, "two orders above the measured worst" (5.3e-16);
   the outside voice of 2026-08-29 pointed out that a naive `N·ε` rounding bound on 123 rows summed
   in two different orders is already 1.4e-14, so a tolerance below it is a contract that can fail on
   arithmetic that is correct, on a workbook that has not arrived. 1e-12 is still eight orders below
   anything printed and cannot fire on rounding below ~10⁴ rows. The two sides are computed
   independently, **as two `gcompute` calls on the same fit**: the left from `gcompute(X, fit,
   over=all, active=1[eligible], comparator=zeros, keys=(_ACTIVE, _COMPARATOR))`; the right from
   `gcompute(X[eligible], fit, over=all-True, active=full(treated), comparator=full(control))` — the
   [§14a]-shaped averaging over the eligible rows under the [§14b] fit, on a design that **contains no
   contraindicated row** (§6.3). Only the second call's `rd` is read; its `n_fit` is `n_eligible` and
   is not carried. `eligible_rd` is that `rd`; `share_eligible` is `n_eligible / n_average`; and the
   three are therefore three things that can disagree after an edit. `_assert_dilution` raises
   `SchemaError`, for T10's reason: an identity that fails is arithmetic broken, not a sparse draw.
   §16.6 spies on `gcompute` to assert both calls happen with those arguments, so `eligible_rd`
   cannot quietly become `rd / share_eligible`.

**The absolute levels are NOT an identity, and the fit choice moves them.** `active_j` and
`comparator_j` include the contraindicated patients' predicted EVT-alone distribution, which DECISION
10 models with their own indicator and the declined third reading would have extrapolated from
eligible patients. Two fits that agree on every `rd_k` — because the contraindicated cancel — can
disagree on every `comparator_j`. That is why the distributions are reported and why the fit was a
decision.

---

## 8. Two labels, both structural

**`exp(β)` is a conditional odds ratio, conditional on `X` and on contraindication status**, and
Stage 12 §8's three enforcements are inherited: no `odds_ratio` field, `measure` carried on the record
as a `Final` string, and a scan asserting the string `"odds_ratio"` is absent from `policy.py`. The
string is `_CONDITIONAL_14B`, distinct from Stage 12's `_CONDITIONAL`, because the conditioning set
differs and a Stage 14 caption must say which — the argument Stage 12 §3.1 makes for `hier.beta`.

**The estimand is an operational policy contrast**, and that is the second label, `_OPERATIONAL`:
`"operational policy contrast [§14b] — what adopting a bridging policy would have delivered in this
cohort; not [§14a]'s ATE in the eligible population, not [§7]'s ATO, and not the biological effect
of IVT"`. It is data on the record, Stage 14 prints it from the record, and §16.8 asserts every
`Policy` carries it. The roadmap's Accept-when — *"labelled an operational policy contrast, distinct
from both [§14a] and [§7]"* — is thereby a field and not an adjective.

**`Policy` has no field called `rd_ate` or `ate`**, and `eligible_rd` is deliberately not called that
either: it is the eligible-subpopulation contrast under a model fitted on a wider population, and
[§14a]'s ATE is Stage 12's record.

---

## 9. Separation, and the guard is on the treatment coefficient alone

### 9.1 The contraindication indicator separates, once in 2000, and it is benign

Ten of the nineteen contraindicated patients are mRS 6; the other nine spread over mRS 2–5. A
stratified resample in which **every drawn contraindicated patient is mRS 6** exists — measured, it
occurred in **1 of 2000** replicates (P6b), with 16 contraindicated rows all at mRS 6. In that draw
`D` perfectly predicts `Y = 6` among the rows it indexes, `δ` runs to **−23.6** on the `P(Y ≤ k)`
scale — past `POLR_MAX_ABS_BETA = 14.0` — and `polr` **converges in 21 iterations on the likelihood
criterion**, against 4–6 on every other draw. Stage 8 §6 measured the same thing on the [§8]
estimator: a separated proportional-odds fit converges and returns.

**The magnitude −23.6 is `POLR_TOL`'s and not the data's.** Under separation the likelihood gain
from pushing `δ` further is roughly `n · e^δ`; the loop stops when that falls below
`POLR_TOL = 1e-8`, which happens near |δ| ≈ 23 whatever the sixteen rows look like. The data facts
in that draw are the **sign**, the **21 iterations** and the **16 rows at one level**; the number is
where the criterion stopped, and every place this document or the log prints it says so (§10.2's grid
row). Stage 14 must not print it as a coefficient.

**And it changes nothing this stage reports.** In that draw the eligible-subpopulation contrast under
the policy fit equals the contrast from an eligible-only fit on the same rows — **measured as 0.0 by
P6b**. The outside voice of 2026-08-29 doubts the bit-exactness: two Newton runs on different designs
stopping on a 1e-8 likelihood tolerance would be expected to agree to ~1e-6–1e-9, and a 0.0 can be a
probe comparing one object to itself. The workbook is not in the reviewing checkout, so **task 4
re-measures this and records the value in §22**; if it is not 0.0, the mechanism below — not the
number — is what the claim rests on, and §16.7 asserts the mechanism at 1e-6. The mechanism is the
one §5.1 relies on: a separated `δ` absorbs the contraindicated rows
entirely — their fitted `P(Y = 6)` is 1 — so `α_k`, `β` and `γ` are determined by the eligible rows
alone, which is the [§14a]-shaped fit. The contraindicated rows' contrast is zero regardless (§7.3),
so every `rd_k` is the eligible contrast times the share, and every `active_j`/`comparator_j` carries a
contraindicated block that is a point mass at mRS 6 — which is what the drawn data say. `|β|` in that
draw is 0.35, unremarkable.

### 9.2 So U8 is T12: the treatment coefficient, key granularity, and `gamma` and `delta` unbounded

`standardise._assert_beta_reportable`'s **argument** transfers whole (Stage 12 §9.1); its **code** is
re-stated as `policy._assert_beta_reportable`, a private with the literal leading token `U8`, for the
raise-site-scan reason §5.3 gives for U1. The reported quantities are averaged probabilities bounded
in [0, 1] by construction — measured, every `rd_k` in [−1, 1] across all draws — and the one quantity
a separated fit makes unprintable is `exp(β)`. **U8 binds `fit.beta[TREATMENT]` at
`POLR_MAX_ABS_BETA`, costs `beta` alone, and does not look at `δ`.** A bound
on the worst coefficient — `outcome._assert_reportable`'s rule, which Stage 12 §24 already declined —
would drop the §9.1 replicate for a nuisance parameter that touched no reported number.

Measured (P6): `|β|` across 2000 replicates runs from below 5e-5 at its minimum to **2.16**, median
0.32; U8 fires **0** times.
`δ` is negative in 1964 of 2000 draws and positive in 36 — a small-sample indicator on 8–31 records
crosses zero — with quantiles −3.0 (1%) to +0.14 (99%) and the one −23.6. The fit entry prints the
coefficient table with `δ`'s role stated as *"contraindication indicator — a NUISANCE, unbounded [§9],
not a result [§4.3]"*.

**`POLR_MAX_ABS_BETA`'s calibration item in `TODOS` has its trigger checked again and it does not
fire**: Stage 12 §21 item 2 rewrote the trigger as *"a Stage 13 population on which any Stage 12
coefficient exceeds 8.79"*, meaning a coefficient the bound is applied to. `|β|` here maxes at 2.16.
`δ` exceeds 14 once and is not bound, so the item stays closed.

---

## 10. Inference [§10, §14]

### 10.1 All four strata are resampled, and the contraindicated are members of the population

Stage 12 §13.1's argument transfers with one extra step. [§14b]'s estimand is a contrast averaged over
*"all patients"*, so the contraindicated patients are members of the target population in the same
sense USZ's eligible patients are: their covariates enter the average, their outcomes enter the fit,
and their *share* is a property of the cohort that an interval over "this cohort" must carry.
`share_eligible` running 0.748 to 0.935 across draws is that uncertainty made visible, and it is why
`share_eligible` is a key with an interval rather than a constant (§3.3).

`bootstrap.resample` is unchanged and `C.BOOT_STRATUM` is passed as is; the stratum is `center`, not
eligibility, so a draw's contraindicated count varies (§10.4) exactly as its per-arm counts do in
Stage 10.

### 10.2 The replicate body — one arm, three failure groups

```
   body(draw):
     D, active, comparator = _regimes(draw)               §5.2, U4; D == draw[CONTRAINDICATED] or raise
     X, dropped = model.design(draw, cols)
     U1         -> FitError; costs EVERY key                the design has no exposure
     fit        = model.polr(X, y)
     FitError   -> costs EVERY key                          nonconvergence
     U8         -> costs `beta` ONLY                        §9.2, key granularity
     U2 on the active vector, before predict               §6.3
     gcompute(over=all,      active, comparator)           22 keys
     gcompute(over=eligible, treated, control)             eligible_rd_k;  share_eligible;  U7
```

**Three failure groups**, Stage 12's four minus the `polr_ri` group: U1 and a pooled `FitError` cost
every key under different buckets (`degenerate_design`, `nonconvergence`); U8 costs `beta`. Nothing
else can fail inside a replicate — U2, U4, U6 and U7 are `SchemaError` and are not caught, on Stage 12
§14's rule.

`bootstrap.bucket` classifies and `bootstrap.collect` reconciles; `Replicate.n_alpha`,
`.polr_iterations`, `.sum_w = n` and `.n_in_model = n` carry the one fit's, and `max_abs_beta` is
keyed `{"policy": max|coef|}` — which **includes `δ`**, so the diagnostic grid shows the §9.1 tail
(23.6 against a bound of 14.0 that does not apply to it) with the row label saying so **and saying
that a separated magnitude reflects `POLR_TOL`, not the data [§9.1]**.

**The aggregator is `standardise.diagnostics`, Stage 12's `_diagnostics` made public — and it is the
only one that can be.** `bootstrap.diagnostics` (`bootstrap.py:832-834`) keeps `max_abs_beta` entries
only for keys in `C.BINARY_OUTCOMES`, so a `"policy"` key is **silently dropped** — the §9.1 tail would
never reach the grid and nothing would raise — and its `or_corrected` comes back keyed by seven outcomes
this stage does not have. Stage 12 §13.2 declined it for the same reason and wrote a generic tally;
that tally has no Stage 12 content in it and is the one Stage 13 calls. The review of 2026-08-29 found
the first draft named no aggregator at all.

### 10.3 `bootstrap.intervals` with `lambda key: False`

Thirty keys, thirty intervals, no p. `C.ci_min_draws(C.CI_LEVEL)` is 40; every key has 2000 draws
(P6), so no interval is withheld for thinness.

**The `mortality` interval is not the reflection of `rd_5`'s** — Stage 12 §7.4 re-measured here (P7):
under `inverted_cdf` the limits differ by **4.3e-4 and 6.3e-5**; under `"linear"` by 2.6e-16.
`mortality` carries its own draws.

### 10.4 The drop rate is zero, measured

```
   2000 replicates, stratified by centre, C.SEED                        [P6, P6b]
     attempted                                          2000
     survived, every key                                2000     100.00%
     U1 (design lost the exposure)                         0
     nonconvergence, polr                                  0
     U8 (POLR_MAX_ABS_BETA on beta)                        0
     design columns dropped                                0     100.00% empty
     level lost (always mRS 5)                            13       0.65%     §6.4
     polr iterations                    4: 15   5: 1949   6: 35   21: 1     §9.1
     route                              likelihood, all 2000
     first_step_norm                    1.68 to 21.6
     contraindicated per draw           8 to 31, median 19
     share_eligible                     0.748 to 0.935
     wall clock                         27.6 s
```

**The one 21-iteration fit is §9.1's separated draw and it converged.** Stage 12 §13.4's reasoning
for a zero drop rate — an unweighted fit on a larger population cannot concentrate its effective
sample — holds here on 123 records with `Σw = n = 123`, and the absolute-tolerance `TODOS` item's
trigger fires a second time with the same negative result: 4–6 iterations on the likelihood criterion
in 1999 of 2000 draws.

---

## 11. The failure taxonomy — the U series, and there are eight

`C`–`S`, `B` and `T` are taken by Stages 2–12. Stage 13's identifiers are **U**. Same rule: a
`FitError` is droppable-and-countable inside a replicate; a `SchemaError` is a contract break and is
never caught.

| id | raised by | class | meaning |
|---|---|---|---|
| **U1** | `policy._assert_exposure_survived`, called by `contrast` and the body | `FitError` | the [§14b] design lost the treatment column (§5.3). Token `U1` → `degenerate_design`. **A literal raise site in `policy.py`** |
| **U2** | `contrast`, the body | `SchemaError` | a contraindicated row carries the treated code in the **active** design — the roadmap's Accept-when (§6.3) |
| **U3** | `population` | `SchemaError` | a contraindicated patient is treated — Stage 4's assertion re-checked at this boundary (§4.1) |
| **U4** | `_regimes`, called by `population`, `contrast` and the body | `SchemaError` | the frame carries no `eligibility` column, a value outside `ELIGIBILITY_ORDER`, **or a stored `CONTRAINDICATED` column that disagrees with `eligibility`** (§5.2) |
| **U5** | `population` | `SchemaError` | the population is empty, or has an empty arm among the eligible (T7's argument). **A population with no contraindicated patient is NOT U5** — see below |
| **U6** | `gcompute` (T10, inherited) | `SchemaError` | a standardised distribution does not sum to 1 or is not monotone |
| **U7** | `policy._assert_dilution`, called by `contrast` and the body | `SchemaError` | `rd_k ≠ share_eligible · eligible_rd_k` at 1e-12 for any `k` (§7.3), or T11's two namings |
| **U8** | `policy._assert_beta_reportable`, called by `contrast` and the body | `FitError` | `|β_treatment| ≥ POLR_MAX_ABS_BETA` (§9.2). Token `U8` → `separation`. **A literal raise site in `policy.py`** |

**A workbook with no contraindicated patient is a data fact, recorded, and not a contract break
(decided 2026-08-29, on the outside voice's finding).** The first draft made it U5's third clause. The
argument against: a `SchemaError` is a bug in the code or the contract, and "this cohort has nobody to
withhold IVT from" is neither — it is the same event §5.3 calls legitimate inside a replicate, and a
Stage 14 driver killed by it would lose Stages 1–12's output for a fact about Stage 13's subject. So
`population` **returns** the frame, entry 1 carries the sentence *"[§14b] has no contraindicated
patient in this cohort: the active regime bridges everyone, `share_eligible` is 1, and this policy
contrast coincides with a [§14a]-shaped contrast under a fit with no indicator; read it with that
stated"*, `contrast` proceeds — `D` dropped as constant, `share_eligible == 1.0`, `rd == eligible_rd`
(§16.16) — and **`Policy.label` still says operational policy contrast**, because that is what it is:
a policy over a cohort in which the policy happens to reach everyone. What §8 forbids is a [§14b]
number printed under a [§14a] label or the reverse; the entry-1 sentence and `share_eligible == 1` on
the record are how Stage 14 tells the reader, and §14 owes that sentence wherever `share_eligible`
prints as 1.

**T10 and T11 are called by their Stage 12 names inside `gcompute` and are not renamed U6/U7a.** They
are the same raise sites; this table lists U6 so that the U series is complete as a reader's
checklist, and §16.14's set assertion ranges over U1–U5, U7, U8 — the sites this module owns.

**`bootstrap.bucket` reads the first token** (`bootstrap.py:519`): U1's and U8's messages lead with
their own identifiers **as string literals at raise sites inside `policy.py`**, because the scan that
keeps `FAILURE_BUCKETS` honest reads them statically (§5.3). Three amendments to `test_bootstrap.py`
follow (§16.17): `FITERROR_MODULES` gains `"policy.py"`; the per-module pin gains
`"policy.py": 2` with tokens `["U1", "U8"]` and the total becomes **29**; and
`test_every_declared_bucket_token_is_ACTUALLY_RAISED_somewhere` is what makes the two new
`FAILURE_BUCKETS` entries reachable rather than decorative. `standardise.py`'s pin stays `["T1", "T12"]`.
The `code`-field `TODOS` item's rewritten trigger — *"the first stage that needs to branch on a failure
kind rather than count it"* — is examined and does **not** fire: this stage counts. **What the stage
paid instead is recorded on that item** (§17 item 9): two three-line guards exist twice because a
`code` attribute would have let them be one.

---

## 12. The audit entries, and there are five

All under existing `data.KINDS`; no kind is added. Every step name is disjoint from `cohort.build`'s
six and `standardise`'s ten, asserted as two empty intersections (§16.9).

| # | kind | step | contents |
|---|---|---|---|
| 1 | `cohort` | `policy_population` | 126 → 125 → 123, per centre × eligibility × arm (the §4.1 table); case identifiers for every removal through `standardise.record_removal`; the statement that **neither** [§3] restriction is applied and why; the contraindicated count per centre; **when that count is 0, §11's no-contraindicated sentence** |
| 2 | `missingness` | `absence_by_policy_column` | `data.absence_by_column` over `POLICY_COVARIATES + (PRIMARY_OUTCOME,)` on the population |
| 3 | `model` | `policy_fit` | the design's columns and what was dropped, the fitted level set, `len(alpha)`, iterations, `converged_on`; the coefficient table with `β` roled *"EXPOSURE — the guarded one [U8]"*, `γ` *"unbounded"*, and `δ` *"contraindication indicator — a NUISANCE, unbounded [§9], not a result [§4.3]"*; `|β|` against `POLR_MAX_ABS_BETA` |
| 4 | `model` | `policy_contrast` | the two regime distributions over `MRS_LEVELS`, the cumulative table, `rd_k`; `eligible_rd_k` and `share_eligible` beside them; **the dilution identity and the two namings as computed differences**; `n_eligible` of `n_average`; the `_OPERATIONAL` label printed from the record |
| 5 | `model` | `policy_replicates` | §10.4's grid with each block labelling its denominator; the collapse count; `max_abs_coef` with the note that it includes `δ` and that the bound does not |

Ledger: 12 before this stage, **5** added — 17 alone; in a Stage 14 driver running Stages 1–13,
Stage 12's 56 plus these 5 = **61**.

**The three tables are `policy.py`'s own renderers — `_fit_table`, `_contrast_table`,
`_replicates_table` — and not Stage 12's.** Stage 12's `_fit_table` hard-codes the roles *"EXPOSURE
— the guarded one [T12]"* and *"gamma, unbounded"* and has no row for `δ`; its `_distribution_table`
is keyed by arm code with `[§14a]` in its row labels and `bridging` / `EVT alone` in its headers, where
this stage's columns are `active` / `comparator`; its `_replicates_table` indexes `draws["hier.sigma"]`
and would `KeyError` on thirty keys. Renderers carry their labels, and a renderer whose labels are
wrong for a caller is not shared with it. The shapes recur — a coefficient table with a role column, a
grid with a labelled denominator per block — and `TODOS`' row-builder item records this stage as the
next instance (§17 item 10). **Two rendering rules, both asserted by §16.17: no cell of any table
contains a `|`** (`data._md_table` does no escaping — Stage 7 shipped `|SMD|` and got an 8/6 header;
so the grid's column is `max_abs_coef`, never `max|coef|`), **and the string `[§14a]` appears in no
Stage 13 table.**

---

## 13. The entry points, written out

```python
def population(df: pd.DataFrame, audit: Audit) -> pd.DataFrame:
    """The [§14b] population: EVERY classified patient, covariate-complete, outcome-observed,
    with the CONTRAINDICATED indicator added. 126 -> 123. Neither [§3] restriction is applied.

    U3, U4, U5. Records entries 1 and 2.
    """


def contrast(pop: pd.DataFrame, audit: Audit) -> Policy:
    """The [§14b] feasible-policy contrast: one unweighted fit on all patients, two regimes.

    U1, U2, U7 (and T10, T11 through `standardise.gcompute`). Records entries 3 and 4.
    """


def inference(pop: pd.DataFrame, audit: Audit) -> bootstrap.Bootstrap:
    """The [§14] bootstrap: one `replicates` call, one arm per draw, thirty keys, no p-value.

    ALL FOUR STRATA. Records entry 5. Calls `bootstrap.replicates`, `bucket`, `collect`,
    `intervals`; implements none of them.
    """
```

The canonical call order — the one place it is written:

```
  pop  = policy.population(classified, audit)        123 rows, all 4 centres     entries 1, 2
  pol  = policy.contrast(pop, audit)                 [§14b]                      entries 3, 4
  boot = policy.inference(pop, audit)                30 keys, 2000 draws, NO p   entry  5
```

---

## 14. Data flow into Stage 14 [§16]

What this stage owes the reporting layer and cannot enforce:

- **The label is printed from `Policy.label`, never written.** [§14b]: *"reported with that stated"*.
- **Wherever the policy contrast appears beside [§14a] or [§7], the three populations are printed as
  numbers** — `Policy.n_average` (123, four centres, 19 contraindicated), `Standardisation.n_average`
  (104, four centres) and `Primary.in_estimate.sum()` (92, three centres) — which is Stage 12 §17.2's
  computed clause with a third party.
- **`share_eligible` is printed beside every `rd_k`**, with `eligible_rd_k`, so the dilution is
  visible (§7.3). A table of `rd_k` alone is a table of a proportion. **Where `share_eligible` prints
  as 1, the entry-1 no-contraindicated sentence prints with it** (§11).
- **[§14a]'s `RD_k` (Stage 12) is THE [§14] eligible-population estimate. `eligible_rd_k` is a
  DIAGNOSTIC — the factor that decomposes [§14b]'s `rd_k` — and Stage 14 prints it in that role
  only** (decided 2026-08-29, on the outside voice's finding that two "eligible-population effects"
  would otherwise sit on adjacent pages). It is never placed in a column with Stage 12's `RD_k`, never
  given a heading that reads as an estimate in its own right, and its interval is printed as the
  interval of a factor. The gap between the two — 2.1e-3 at the point estimate, up to 9.7e-2 in a
  draw (§5.1) — is printed once, as a diagnostic beside the decomposition, with the sentence that the
  two are fitted on different populations; Stage 14 decides the layout and nothing else.
- **A Stage 14 driver does not need to treat Stage 13 as skippable**: a cohort without
  contraindicated patients is reported, not refused (§11).
- **`exp(β)` is conditional on `X` and on contraindication status**, and `Policy.measure` says so; it
  is never in a column headed by `Primary.odds_ratio` or by Stage 12's `beta`.
- **`δ` is printed in the fit table as a nuisance and nowhere else.** No interval, no sentence about
  contraindicated prognosis (§4.3).
- **The collapse count (13 of 2000) is printed beside the `rd_5` interval**, as Stage 12's is.
- **No [§14b] key carries a p-value and none joins a [§13] family.**

---

## 15. What Stage 13 amends in Stages 1–12

**Two shipped modules and three test modules; no estimator changes; `cohort.py`, `eligibility.py`,
`model.py`, `bootstrap.py`, `balance.py`, `outcome.py`, `propensity.py`, `derive.py`, `data.py` and
`sensitivity.py` are untouched.**

| File | Amendment | Why |
|---|---|---|
| `standardise.py` | the arithmetic of `_standardise` → public `gcompute` returning the new `GComputation`, regimes as per-row arguments, `keys=` for the distribution dict, a length check; `_standardise` stays private as the wrapper; `_regime_design` accepts a vector; `_record_removal` → `record_removal` with a caller-neutral message; `_diagnostics` → `diagnostics`. Public surface 5 → 8 | §6.1, §4.2, §10.2. `all_centre`/`hierarchical`/`support`/`population`/`inference` unchanged; §16.11 asserts byte-identity of the sixty-nine intervals and all ten audit entries against two pinned digests |
| `config.py` | `CONTRAINDICATED`, `POLICY_COVARIATES`, `"U1"`, `"U8"` in `FAILURE_BUCKETS` | the fence below |
| `tests/test_standardise.py` | surface list 5 → 8 and the privates list (`_diagnostics` out); `STAGE12_INTERVALS_SHA256` and `STAGE12_ENTRIES_SHA256` captured **before** the refactor and the two `STAGE12_SLOW` tests that check them; `GComputation`/vector-regime/length-check/message tests; the private-call set extended to `policy.py` and asserted empty there | §0.1, §16.11, §16.15 |
| `tests/test_config.py` | `POLICY_COVARIATES == STANDARDISATION_COVARIATES + (CONTRAINDICATED,)`; `CONTRAINDICATED not in PS_COVARIATES`; `CONTRAINDICATED not in CATEGORICAL`; the two tokens | §5.1, §16.2 |
| `tests/test_bootstrap.py` | `FITERROR_MODULES += ("policy.py",)`; per-module pin gains `"policy.py": 2`, total 29, tokens `["U1", "U8"]` | §11, §16.17 |
| `statistical_analysis_plan.md` | [§14b] amendment, 2026-08-27 — **landed in `37051fc`** | DECISION 10, §5.1 |
| `implementation_roadmap.md` | Stage 13 block — **landed in `37051fc`** | §23 |
| `TODOS.md` | four new items — **landed in `37051fc`**; two existing items amended by the review of 2026-08-29 | §17 items 9, 10 |
| `../out/stage0_data_inventory.md` | DECISION 10 | §5.1 |

```python
# config.py — the additions, as they are to be pasted. Stage 13 §15

# --- [§14b] the contraindication indicator and the policy model's covariates ---------------------
#
# The NAME of the 0/1 column `policy.population` derives from `eligibility` — never from the raw
# `ivt_contraindicated` flag, because Stage 4 is the one classifier and a stage re-reading the flag
# would be a second one that agrees with the first on this workbook (Stage 13 §5.2). A float 0/1
# like `atrial_fib`, NOT in CATEGORICAL, so `model.design` takes it as one column.
CONTRAINDICATED: Final[str] = "contraindicated"                                # [§14b]

# A computed view, never a fifth literal list. DECISION 10: [§14b]'s model is [§14a]'s covariates plus
# the indicator, and test_config.py holds it to exactly that — [§15]'s "inside §14 only" permission
# for omitting centre covers both subsections, and adding the indicator to PS_COVARIATES would be [§3]
# restriction 2 undone inside the propensity model, which the same test forbids.
POLICY_COVARIATES: Final[tuple[str, ...]] = STANDARDISATION_COVARIATES + (CONTRAINDICATED,)   # [§14b]
```

**Two `FAILURE_BUCKETS` tokens**: `"U1"` → `"degenerate_design"`, `"U8"` → `"separation"`, beside Stage
12's four, with the same first-token comment. Stage 13 adds **no** token leading with `polr:`. Both
tokens are **string literals at raise sites in `policy.py`** (§5.3, §9.2, §11); no token anywhere in
the repository is a parameter, because the raise-site scan cannot read one.

---

## 16. Acceptance criteria

### 16.0 `tests/fixtures_stage13.py` — the one place its contents are enumerated

Given as code:

- `mixed_frame()` — a classified population over four declared centres, one with no treated patient,
  contraindicated patients at two centres and none treated, all seven mRS levels occupied.
- `all_eligible_frame()` — `mixed_frame()` with every patient eligible: `population` returns it with
  entry 1's no-contraindicated sentence and `contrast` proceeds with `share_eligible == 1` (§11,
  §16.10).
- `treated_contraindicated_frame()` — one contraindicated patient with the treated code, for U3.
- `separated_indicator_frame()` — every contraindicated patient at mRS 6, so `δ` separates while `β`
  does not; §9.1 on a fixture.
- `no_contraindicated_frame()` — `mixed_frame()` with its contraindicated rows removed, **for the
  replicate body and `contrast` only, never `population`** (which raises U5 on it): the §5.3 draw in
  which `D` is dropped as constant, `share_eligible == 1` and `rd == eligible_rd` (§16.16).
- `POLICY_SEED`, `N_CONTRAINDICATED`.

Five functions and two constants. Every frame declares all four centres so that
`bootstrap.resample` under `C.BOOT_STRATUM` draws from four strata, as the workbook does.

### 16.1 The population is [§14b]'s and nobody else's

`population(classified)` returns **123** rows over **four** centres with **19** contraindicated;
`standardise.population` returns 104 and `cohort.build` 93 on the same frame — all three in one test,
**`DATA_GATED`** (the committed fixtures cannot reach `cohort.build`'s both-arms postcondition; only
the workbook can).
**The 123 is the 104 plus exactly the 19 records `eligibility.retained` is False for**, asserted as
set equality on `case_id`. Every contraindicated record has `ivt == 0`. The three [§11] losses are
Stage 12's three, asserted by identifier equality with `standardise.population`'s removals. Neither
`cohort`, `propensity` nor `derive` is imported, by scan.

### 16.2 The config additions

The three relations and two tokens of §15. The load-bearing one is `CONTRAINDICATED not in
PS_COVARIATES`.

### 16.3 The estimand keys are derived and not counted

`len(keys) == 2*len(C.MRS_THRESHOLDS) + 2 + 2*len(C.MRS_LEVELS) + 2`; every interval's `p is None`;
`"contraindication" not in keys` by name; no key starts with `dist`.

### 16.4 The fit is unweighted, on all rows, and the call site is the assertion

`model.polr` called with two positional arguments and no `w`; `len(X) == len(pop)` — no row is
excluded from the fit.

### 16.5 The regimes

On `mixed_frame()`: the active vector has `1[eligible].sum()` ones; the comparator is all zeros; both
have `len(X)` rows. **U2 fires** when a contraindicated row of the active design is set to the treated
code by the test, before any prediction runs — asserted by the message and by `ordinal_probabilities`
not having been called (monkeypatched counter).

### 16.6 The dilution identity and the two namings

`abs(rd[k] - share_eligible*eligible_rd[k]) < 1e-12` for every `k`, on the workbook and on
`mixed_frame()`; `mrs_0_2 == rd[2]` exactly; `abs(mortality + rd[5]) < 1e-15`. **And the per-row
contrast of every contraindicated row is exactly `0.0`**, asserted on the arrays before averaging.
**U7 fires** when `_assert_dilution` is handed an `eligible_rd[2]` perturbed by 1e-13. **And the two
sides are shown to be two computations, not one and a division**: a spy wrapping the real
`standardise.gcompute` (monkeypatched onto the module `policy` imports it from) records every call
`contrast` makes on `mixed_frame()` and asserts there are **exactly two** — one with `over` all-True,
`active` equal to `_regimes`' vector and `comparator` all-`control`, keyed `(_ACTIVE, _COMPARATOR)`;
one whose **design has exactly `n_eligible` rows and no row with `CONTRAINDICATED == 1`**, `over`
all-True, and both regimes scalar broadcasts of `treated` and `control` — and that `Policy.rd` is the
first call's `rd` and `Policy.eligible_rd` the second's, by identity of values. A `rd / share_eligible`
shortcut makes the second call disappear and the test fail; a full-design second call fails the row
count. And **U2 fires** when the stored `CONTRAINDICATED` column is set to 1 on an eligible row after
`population` ran and `_regimes` is bypassed — the two-source check of §6.3.

### 16.7 The absolute levels move with the fit and the contrast does not

On `mixed_frame()`, the contraindicated rows' `comparator_j` block under the [§14b] fit differs from
the prediction a [§14a]-shaped fit (eligible rows only, no indicator) makes for those rows; the
`rd_k` under the two fits differ too (the shared coefficients move) — **but on
`separated_indicator_frame()` the eligible contrast under the two fits agrees to 1e-6**, which is
§9.1's mechanism as an assertion. 1e-6 and not tighter because the two fits are two Newton runs on
different designs stopping on `POLR_TOL = 1e-8`; a bit-level assertion here would test the optimiser's
last step, not the mechanism.

### 16.8 The two labels

`Policy.measure is _CONDITIONAL_14B`, `Policy.label is _OPERATIONAL`, `hasattr(pol, "odds_ratio")`
is False, and `"odds_ratio"` does not appear in `policy.py` outside a comment.

### 16.9 The audit entries and their step names

Five entries in §13's order; the step-name set is disjoint from `cohort.build`'s and from
`standardise`'s, both as empty intersections; entry 1 names exactly as many identifiers as its `n`;
the Stages 1–13 ledger totals 61 in a driver that runs Stages 12 and 13 — **`STAGE12_SLOW`-gated**,
because Stage 12's entry 10 is its 4.2-minute bootstrap. The 17-entry Stage 13-alone count is the
ungated form of the same assertion.

### 16.10 The indicator is derived from `eligibility`

Blanking `ivt_contraindicated` on every record after `classify` has run leaves `population`'s output
unchanged; **U4 fires** on a frame with a third eligibility value and on one with no `eligibility`
column; **U3 fires** on `treated_contraindicated_frame()`; **U5 fires** on an empty frame and on a
frame whose eligible patients are all in one arm. **U5 does NOT fire on `all_eligible_frame()`**:
`population` returns every row, entry 1's detail contains the no-contraindicated sentence and reports
0 contraindicated per centre, and `contrast` on it returns `share_eligible == 1.0`, `dropped ==
(C.CONTRAINDICATED,)` and `rd == eligible_rd` exactly, with `label is _OPERATIONAL`.

### 16.11 Stage 13 changes nothing upstream, asserted by number — and the mechanism is two digests

The Stage 1–12 pipeline run with and without this stage: every landed number, the first 56 audit
entries, and **Stage 12's sixty-nine intervals — every limit, level, method and draw count — and its
ten audit entries, byte-identical across the `gcompute` extraction and the two renames.**

**A test cannot compare against code that no longer exists, so the comparison is against two SHA-256
digests captured on this branch BEFORE task 2 touches `standardise.py`** — `STAGE12_INTERVALS_SHA256`
over the sixty-nine `Interval`s canonicalised exactly as `test_standardise.py:1046-1053` canonicalises
Stage 10's twenty-six (`{key: [repr(lo), repr(hi), level, method, n_draws, repr(p)]}`,
`json.dumps(sort_keys=True, separators=(",", ":"))`), and `STAGE12_ENTRIES_SHA256` over the rendered
markdown of entries 1–10 of a `standardise`-only run, concatenated in ledger order. Both pinned in
`test_standardise.py` beside `STAGE10_INTERVALS_SHA256`, on its precedent and for its reason: **a hash
quotes no interval limit**, so §4.3 and Stage 10 §4.5 hold. Two digests and not one because entry 4's
table is where a wrong `keys=` default or a swapped minuend would show first, and the existing 56-entry
prefix test checks count and order, not content. Both tests `STAGE12_SLOW`-gated (the run is ~4.2 min)
and `DATA_GATED`, as Stage 12 §20.11 is. **The review of 2026-08-29 found the first draft named the
assertion and no mechanism**; §21 task 2a is the capture step, and it is ordered before task 2 because
a digest captured after the refactor proves only that the refactor agrees with itself.

### 16.12 Separation

On `separated_indicator_frame()`: `polr` converges; `|δ| > POLR_MAX_ABS_BETA`; **U8 does not fire**;
every `rd_k` is finite and in [−1, 1]; and a companion applying `outcome._assert_reportable`'s
worst-coefficient rule to the same fit **does** raise — so the test documents that the guard's scope
is a choice and names the rule it declined. **U8 fires** on a constructed fit whose treatment
coefficient is set past the bound, costing `beta` and nothing else across a replicate.

### 16.13 The replicate loop

Drawn-frame sequence byte-identical to a `None`-body's at the same seed, frame and stratum; all four
strata in every draw at exact size; the contraindicated count varies across draws and is never zero
on the workbook (asserted `min > 0`, data-gated); the three failure groups each cost exactly their
keys, driven by constructed failures; `n_attempted` reconciles for all thirty keys.

### 16.14 Every U identifier fires, and the set is the assertion

`{U1, U2, U3, U4, U5, U7, U8}` — the sites `policy.py` owns — equals the set exercised, on §16.9's
pattern. Each fires on a constructed input; U1 on a frame where every patient is treated (Stage 12's
`no_exposure_frame()` pattern), U2 as §16.5, U3–U5 as §16.10, U7 as §16.6, U8 as §16.12. The set is
collected by scanning `policy.py`'s `raise` sites for `C.SchemaError` and `model.FitError` messages
whose first token matches `U\d`, so an identifier added to the code and not to a test fails here.

### 16.15 The `standardise.py` refactor (added 2026-08-29)

- **`GComputation` is arithmetic only**: the record has exactly the seven fields §6.1 lists, by
  `dataclasses.fields`; `hasattr(g, "measure")`, `hasattr(g, "population")` and
  `hasattr(g, "conditional_log_odds")` are all False.
- **Scalar broadcast equals the wrapper**: on `fixtures_stage12.pooled_frame()`'s fit,
  `gcompute(X, fit, over, np.full(n, treated), np.full(n, control), label=..., keys=(1, 0))` equals
  `_standardise(X, fit, over, label, dropped)` on every one of the seven shared fields, by `==` on
  the dicts and floats — not `isclose`; the same floats in the same order.
- **A regime vector of the wrong length raises `SchemaError`** naming both lengths; a scalar passed
  where a vector is expected also raises (the vector form is the only form).
- **`keys` reaches the record**: `gcompute(..., keys=("a", "b")).distribution` has exactly the keys
  `{"a", "b"}` and `rd` is `cumulative["a"] − cumulative["b"]` at every threshold.
- **`record_removal`'s message names no section**: on a frame that removes 2 and names 1 the
  `SchemaError` text contains the step name and neither `"§14a"` nor `"§14b"`.
- **The surface**: `test_the_public_surface_is_FIVE_NAMES…` becomes `…EIGHT_NAMES…` with the
  five plus `gcompute`, `record_removal`, `diagnostics`, each renamed in place so the source order
  is the sections' — `["record_removal", "population", "gcompute", "all_centre", "support",
  "hierarchical", "diagnostics", "inference"]` — and the privates list minus `_diagnostics`;
  the module's dataclasses are `[Standardisation, Support, Hierarchical, GComputation]`.
- **The private-call set**: the assertion at `test_standardise.py:962` is run over `policy.py`'s
  source and the set is `{}`.
- **The two digests**, §16.11.

### 16.16 `_regimes` and the no-contraindicated draw (added 2026-08-29)

- On `mixed_frame()` and on every resample of it: `D + active / treated == 1` on every row,
  `comparator` is all-`control`, `D` equals `population`'s stored `CONTRAINDICATED` column exactly.
- **U4 fires** when the stored column is flipped on one row after `population` ran and `contrast` is
  called — the disagreement clause — and the message names the count of disagreeing rows.
- On `no_contraindicated_frame()`, run through the replicate body (not `population`): `dropped ==
  (C.CONTRAINDICATED,)`, `share_eligible == 1.0` exactly, `rd == eligible_rd` at every threshold
  exactly, `active` is all-`treated`, and no `U` identifier fires. This is §5.3's *"legitimate draw"*
  as an assertion, and the same outcome `population` + `contrast` produce on `all_eligible_frame()`
  (§16.10) — the two paths agree that a cohort with nobody to withhold IVT from is reported, not
  refused.

### 16.17 The tally, the renderers and the scan (added 2026-08-29)

- **`standardise.diagnostics` keeps every `max_abs_beta` key**: over three constructed `Replicate`s
  with `max_abs_beta={"policy": v}` it returns `Diagnostics.max_abs_beta == {"policy": array}` of
  length 3; the same three through `bootstrap.diagnostics` return `{}` — the test documents why the
  generic tally and not the bootstrap one is the one called.
- **`_fit_table`**: one row per design column plus one per dropped; the `TREATMENT` row's role
  contains `"U8"`; the `CONTRAINDICATED` row's role contains `"NUISANCE"` and `"[§4.3]"`; every other
  covariate's contains `"unbounded"`.
- **`_contrast_table`**: columns headed by `_ACTIVE` and `_COMPARATOR`, never by an arm label; rows
  for `eligible_rd_k`, `share_eligible` and the dilution residual `|rd_k − share·eligible_rd_k|` as
  computed differences.
- **`_replicates_table`** renders on a `Bootstrap` whose keys are the thirty of §3.3 (no `hier.*`),
  every block row labelling its denominator, with a `max_abs_coef` row noting the bound does not
  apply to `δ`.
- **Rendering rules over all three tables and all five entries**: no cell contains `|`; header,
  separator and every body row have equal cell counts in the rendered log; the string `[§14a]`
  appears in no Stage 13 table cell and in no detail line outside `_OPERATIONAL` and
  `_NO_CONTRAINDICATED` — the two sentences this document prescribes verbatim and requires printed,
  each of which names [§14a] to say what the contrast is *not*.
- **The scan**: `FITERROR_MODULES` includes `"policy.py"`; `raise_sites()` returns 29; the
  per-module pin is Stage 12's with `"policy.py": 2` added; `sorted(t for m, _, t in sites if m ==
  "policy.py") == ["U1", "U8"]`; `standardise.py`'s stays `["T1", "T12"]`.

---

## 17. Known gaps carried forward

1. **The [§14a] transport assumption is inherited whole and nothing here tests it.** By design.
2. **NEW — the indicator is a main effect.** DECISION 10 lets contraindication shift the cutpoints
   and nothing else; a contraindicated patient's `age` or `nihss_baseline` is assumed to act as an
   eligible patient's does. Nineteen records at two centres cannot test that and this document does
   not try. **Trigger:** a reviewer asking whether the contraindicated patients' EVT outcome model is
   the same model. Recorded as a [§14] question.
3. **NEW — the absolute levels depend on how the contraindicated are modelled; the contrast does
   not** (§7.3). A hybrid plug-in — use each contraindicated patient's *observed* outcome under both
   regimes, since `Y^0 = Y` for them — would give identical `rd_k` and different `active_j`,
   `comparator_j`. **Put to the PI 2026-08-29 and CLOSED the same day: the PI confirmed the model
   prediction (§18).** The outside voice called the plug-in the natural estimator for those fourteen
   keys; the SAP's *"predict every resampled patient under both regimes"* and its *"model-based
   transported results"* label say otherwise, and the PI kept the plan. No longer a gap.
4. **The eligible contrast under the [§14b] fit differs from [§14a]'s by up to 9.7e-2 in a draw.**
   DECIDED 2026-08-29 (§14): [§14a]'s `RD_k` is the estimate, `eligible_rd_k` is a diagnostic factor,
   the gap is printed once as a diagnostic. No longer open.
5. **NEW — the separated-indicator draw (§9.1).** One replicate in 2000 has `δ = −23.6`; it is benign
   for every reported quantity and is measured to be so, but a workbook where the contraindicated are
   more often all-mRS 6 in a stratum would make it frequent and the fit slower (21 against 5
   iterations). **Trigger:** the rate exceeding 1%, or a workbook where a whole centre's
   contraindicated patients share one outcome.
6. **`RD_5` and `RD_4` coincide in 0.65% of draws** — Stage 12 §21 item 8's item gains a second
   instance at a lower rate, for the reason §6.4 gives.
7. **The three removal-ledger implementations** (§4.2): `cohort._record_removal` private,
   `standardise.record_removal` public, and this stage as the third caller. **Trigger:** fired — the
   next change to either moves the rule into `data.py`.
8. **Every rate here is the workbook's.** Stage 12 §19's closing sentence applies to Stage 14 in turn:
   nothing transfers except the shapes.
9. **NEW (review, 2026-08-29) — two three-line guards exist twice.** `policy._assert_exposure_survived`
   and `policy._assert_beta_reportable` restate Stage 12's T1 and T12 with the tokens `U1` and `U8`,
   because `test_bootstrap.py`'s raise-site scan can only read a literal first token (§5.3). A
   `FitError` with a `code` attribute would let each be one public guard taking the code. Recorded on
   `TODOS`' `code`-field item as a cost that item's fix removes. **Trigger:** that item's own.
10. **NEW (review, 2026-08-29) — three renderer shapes exist twice.** `policy._fit_table`,
    `_contrast_table` and `_replicates_table` are Stage 12's three with this stage's labels (§12).
    Recorded on `TODOS`' row-builder item as the next instance. **Trigger:** the next stage adding a
    coefficient table or a replicate grid, which is Stage 14 if it renders any.

---

## 18. What Stage 13 deliberately does not decide

- **Whether the seven decomposition keys are reported or logged** (§3.3). The argument for reporting
  is [§14b]'s *"operational"* sentence; the PI may read it more narrowly. PI-reversible: 30 → 23 keys,
  no code path.
- ~~Whether Stage 14 prints the gap between `eligible_rd_k` and Stage 12's `RD_k`~~ — **decided
  2026-08-29**, §14: printed once, as a diagnostic; `RD_k` is the estimate.
- ~~The contraindicated patients' distributions: model prediction or observed outcome?~~ — **decided
  2026-08-29, PI: model prediction, DECISION 10 confirmed unamended.** Put to the PI because the
  outside voice called the observed-outcome plug-in the natural estimator for the fourteen
  `active_*`/`comparator_*` keys (`Y^0 = Y` is observed for the nineteen). The SAP as written
  prescribes the prediction three times over — *"same standardisation machinery"*, *"predict every
  resampled patient under both regimes"*, *"label the §14 estimates as model-based transported
  results"* — and the PI confirmed it. The plug-in would change fourteen keys and no `rd_k`; it is not
  built and needs a [§14b] amendment to exist (§19). Recorded beneath DECISION 10 in
  `../out/stage0_data_inventory.md` as its dated confirmation (task 6).
- **Whether the contraindication indicator should interact with anything** (§17 item 2). [§6] says
  linear terms only; the PI may amend [§14b].
- **Whether the pipeline is parallelised.** Not moved by this stage (27.6 s).

---

## 19. NOT in scope for Stage 13

| Considered | Why declined |
|---|---|
| A support check for [§14b] | [§14a]'s support check exists because never-IVT-centre patients receive a *predicted IVT outcome* whose covariates may lie outside the treated range. No contraindicated patient receives one (U2), and the eligible patients' support is Stage 12 §10's, already reported. A second box over the same 39 treated patients would be the same table under a second heading |
| The treated-support-restricted arm | Prescribed by [§14a] as a sensitivity *on the standardisation population*; [§14b] prescribes no sensitivity. Adding one is a [§14] amendment. Stage 12 §21 item 14 is the standing warning against adding it from symmetry |
| The random centre intercept arm | Same: unprescribed for [§14b], 95% of Stage 12's compute, and [§14a] itself calls it *"a check, not an upgrade"* |
| An `eligible_only=` keyword on `standardise.population` | §0.3 item 1. Two populations under one step name with a switch |
| Reading `ivt_contraindicated` | §5.2. One classifier |
| A p-value on any key; adding keys to a [§13] family | Stage 12 §3.3, §17.2, unchanged |
| A hybrid plug-in for the contraindicated | §17 item 3. **Put to the PI 2026-08-29; the PI confirmed the model prediction (§18).** Not built; would need a [§14b] amendment |
| Bounding `δ` | §9.2. A nuisance coefficient that touches no reported number |

---

## 20. What already exists, and what to lift

| Need | Exists | Lift |
|---|---|---|
| The [§14a] covariates | `C.STANDARDISATION_COVARIATES` | Extend by one name in `config.py`; never write the list out |
| The eligibility predicate | `eligibility.retained` (`eligibility.py:290`) | Call it for the regime and the indicator; **never filter on it** |
| Covariate completeness | `model.complete_cases` | Call it with `STANDARDISATION_COVARIATES` — the indicator is complete by construction |
| The design, the fit, the prediction | `model.design`, `model.polr`, `model.ordinal_probabilities` | Call them; nothing in `model.py` changes |
| The g-computation, re-expansion, T10, T11 | `standardise._standardise`, `_expanded`, `_arm_probabilities` | **Extract the arithmetic into public `gcompute` → `GComputation`, pass the regimes in; `_standardise` stays as the wrapper** (§6.1). Do not copy |
| The removal ledger | `standardise._record_removal` | **Rename → `record_removal`, genericise the message** (§4.2) |
| The exposure-survived guard | `standardise._assert_exposure_survived` | **Do NOT call it.** Restate as `policy._assert_exposure_survived` with literal token `U1` (§5.3) — the raise-site scan reads literals |
| The treatment-coefficient bound | `standardise._assert_beta_reportable` | **Do NOT call it.** Restate as `policy._assert_beta_reportable` with literal token `U8` (§9.2), same reason |
| The `Replicate` → `Diagnostics` tally | `standardise._diagnostics` | **Rename → `diagnostics`** and call it (§10.2). `bootstrap.diagnostics` drops non-outcome `max_abs_beta` keys and must not be called |
| The three audit tables | `standardise._fit_table`, `_distribution_table`, `_replicates_table` | **Do NOT call them** (§12). Write `policy._fit_table`, `_contrast_table`, `_replicates_table` with this stage's roles, keys and labels |
| The resampler, buckets, reconciliation, intervals | `bootstrap.replicates`, `bucket`, `collect`, `intervals` | Call them |
| The audit ledger | `data.Audit`, `data.absence_by_column` | Call them. Distinct step names |
| The worst-coefficient rule | `outcome._assert_reportable` | **Do NOT call it** (§9.2); §16.12 uses it as the companion that must raise |
| The digest pattern | `test_standardise.py:1006-1053`, `STAGE10_INTERVALS_SHA256` | Copy the canonicalisation exactly for the two Stage 12 digests (§16.11) |

---

## 21. Implementation tasks

In order, each landing green before the next begins.

1. **`config.py`** — `CONTRAINDICATED`, `POLICY_COVARIATES`, two tokens (§15). `test_config.py`'s
   assertions. **Note:** adding `"U1"` and `"U8"` to `FAILURE_BUCKETS` makes
   `test_every_declared_bucket_token_is_ACTUALLY_RAISED_somewhere` red until task 4 lands the raise
   sites; land tasks 1–5 on one branch and accept that test red in between, or land the two tokens
   with task 4. The spec prefers the second: **move the two tokens to task 4.**
2. **2a — capture the two digests FIRST.** With `standardise.py` untouched, run the Stage 1–12
   pipeline at `N_BOOT`, compute `STAGE12_INTERVALS_SHA256` and `STAGE12_ENTRIES_SHA256` exactly as
   §16.11 canonicalises them, and commit the two constants and their two gated tests to
   `test_standardise.py`. They pass trivially at this commit; that is the point. **The gate is a
   run, not a test**: `STAGE12_SLOW=1 uv run pytest -k SHA256 -rs` must report the two as passed and
   zero as skipped — a skipped gated test is green and proves nothing — here and again after 2b.
   **2b — the `standardise.py` refactor** — `GComputation`; `gcompute` with regime arguments, `keys=`
   and the length check; `_standardise` reduced to the wrapper; the vector `_regime_design`;
   `record_removal` with the neutral message; `diagnostics`. Surface and privates lists in
   `test_standardise.py`; §16.15. **Run the two digest tests here, before `policy.py` exists**, so a
   discrepancy is unambiguously the refactor's — Stage 12 §25's step-3 argument.
3. **`policy.population`** and `_regimes` (§4, §5.2) and entries 1–2. §16.1, §16.9, §16.10, §16.16.
4. **`policy.contrast`** (§5–§9), the three guards, the two renderers for entries 3–4, and the two
   `FAILURE_BUCKETS` tokens with `test_bootstrap.py`'s scan amendments. §16.3–16.8, §16.12, §16.17.
   Not gated: the PI confirmed DECISION 10 on 2026-08-29 (§18). Task 4 also re-measures §9.1's
   separated-draw gap and records it in §22.
5. **`policy.inference`** (§10), `_replicate`, `_replicates_table`, entry 5. §16.13, §16.14, §16.17.
6. **DECISION 10 in `../out/stage0_data_inventory.md`** (§5.1), if not already written, **and its
   dated confirmation of 2026-08-29** (§18). The SAP
   amendment, the roadmap block and the four `TODOS` items landed with this document in `37051fc`;
   the two `TODOS` amendments of §17 items 9–10 landed with this revision. Nothing else to write.

---

## 22. Verification record

Nine probe groups, run 2026-08-27 against the landed Stages 1–12 on the v7 workbook (`DATA_SHA256`
unchanged), Python 3.12.12 / numpy 1.26.4 / pandas 2.3.3, in the checkout that holds `data/`.

| Probe | What it measured | Result | Section |
|---|---|---|---|
| P1 | the population ledger | 126 → 125 → 123; the 3 losses all eligible and identical to Stage 12's; 19 contraindicated at HUG (11) and USZ (8), none treated; mRS 5 carries 5 of 123, 2 of them contraindicated; 10 of 19 contraindicated are mRS 6 | §4.1, §6.4, §9.1 |
| P2 | the design and fit | 11 columns, rank 11, 0 dropped; 5 iterations, likelihood, 0 rescales/halvings, first step 8.2; `δ` negative and the largest coefficient by >2× | §5.2, §5.3, §4.3 |
| P3 | the regimes and the identity | active bridges 104 of 123; U2 holds; per-row contraindicated contrast exactly 0.0; dilution identity to 8.9e-17; sums to 1 exactly; monotone; `|mortality + rd5|` 1.1e-16 | §6.2, §6.3, §7.3, §7.2 |
| P5 | against Stage 12's fit | eligible contrast differs by 2.1e-3 max, same sign at every threshold; `|γ|` differs by up to 0.29 | §5.1 |
| P6 | 2000 stratified replicates, one arm | 2000/2000 survive; 0 U1 / 0 nonconvergence / 0 U8; 0 dropped columns; 13 lose mRS 5 (0.65%) with `rd_4 == rd_5` in exactly those 13; iterations 4–6 in 1999, 21 in 1; contraindicated 8–31 per draw; `share_eligible` 0.748–0.935; identity to 5.3e-16; 30 intervals, all `p is None`; 27.6 s | §6.4, §10.4, §3.3, §2 |
| P6b | the indicator's tail | `δ` beyond 8.79 in 1 draw (−23.6, 21 iterations, all 16 contraindicated at mRS 6); positive in 36; in that draw the eligible contrast equals an eligible-only fit's to 0.0; across draws the two fits' eligible contrasts differ by median 1.1e-2, max 9.7e-2; `|β|` 0.35 there, 2.16 max overall | §9, §5.1 |
| P7 | the reflection check | `inverted_cdf` 4.3e-4 / 6.3e-5; `linear` 2.6e-16 | §10.3 |
| P8 | audit step names | 12 entries before this stage; `cohort`'s and `standardise`'s step sets enumerated; the five new names disjoint from both | §12 |

**Three probes changed what this document says** — P3 (§7.3), P6b (§9), P5 (§5.1) — and each is named
in the Status block.

**Re-measured at task 4 (2026-08-29), through the landed `policy.py`** — 2000 stratified replicates,
32.2 s: 2000/2000 survive, 0 U1 / 0 nonconvergence / 0 U8; 13 lose mRS 5 with `rd_4 == rd_5` in
exactly those 13; iterations 4: 15, 5: 1949, 6: 35, 21: 1; `share_eligible` 0.748–0.935; one draw
with `|δ| > POLR_MAX_ABS_BETA` (−23.61, 21 iterations, 16 contraindicated rows, `|β|` 0.35). **In that
draw the eligible contrast under the policy fit and under an eligible-only fit differ by 1.1e-16 at the
worst threshold — not the 0.0 the probe reported.** The outside voice was right that the original
0.0 was not bit-exactness of two Newton runs; the mechanism §9.1 states is what the claim rests on,
and §16.7 asserts it at 1e-6 on the fixture.

**What was NOT measured, and is therefore not claimed.** No `Policy` record and no `policy.py` exist:
every number above comes from the landed `standardise` privates and `model` functions driven by probe
scripts, so the entry points, the U series and the audit entries are specified and not yet run. The
probes were re-run end to end after drafting and every number reproduced exactly, timings excepted.

### 22.1 The engineering review — what it checked, and what it did NOT

Reviewed once, 2026-08-27, by the author's own consistency pass: every count against its arithmetic
(30 keys, 8 identifiers, 5 entries, 17/61 ledger), every quoted number against the probe output,
every claim about landed code against the code (`design`'s factor handling at `model.py:245-256`,
`bucket`'s first-token rule at `bootstrap.py:519`, `retained` at `eligibility.py:290`). **No
independent cross-model review was run, and the §21-item-14 sweep was run on this document's own
decisions with the following result**: the key set (§3.3), the population (§4.1), the fit (§5.1), the
regime (§6.2), the arms declined (§19) and what travels to Stage 14 (§14) were each checked against
[§3], [§14], [§15] and [§16] for a sentence that speaks directly; §5.1 and §19 are where the plan
speaks and this document cites it, and §3.3's seven decomposition keys are the one place this
document reasons *from* the plan's sentence rather than *reading* it — which is why they are
PI-reversible (§18) rather than settled.

### 22.2 The engineering review of 2026-08-29 — nine findings, all folded

Run as `/plan-eng-review` against the landed Stages 1–12 with every claim checked against the code
it names. Scope accepted as written (twelve files is the stage pattern, not sprawl). Findings, each
with where it landed:

| # | Finding | Verified at | Folded in |
|---|---|---|---|
| 1 | A parameterised token on `_assert_beta_reportable` crashes the raise-site scan (`assert isinstance(first, ast.Constant)`), and §0.1's empty-private-set promise contradicted §20's "call it" | `test_bootstrap.py:788-793` | U1 and U8 are `policy.py`'s literal raise sites; §5.3, §9.2, §11, §15, §20 |
| 2 | `gcompute` returning `Standardisation` hands Stage 13 a record whose `measure` says [§14a] and whose `population` has three declared values | `standardise.py:86-128`, `test_standardise.py:1345` | `GComputation`; `_standardise` stays as wrapper; §6.1 |
| 3 | No diagnostics aggregator or table renderer named; `bootstrap.diagnostics` silently drops a `"policy"` `max_abs_beta` key; Stage 12's renderers would `KeyError` or mislabel | `bootstrap.py:832-834`, `standardise.py:640-684, 1327-1372` | `standardise.diagnostics` public; three own renderers; §10.2, §12, §20 |
| 4 | SAP, roadmap and four `TODOS` items listed as pending had landed in `37051fc` | `git show --stat 37051fc` | §1, §15, §21 task 6, §23 |
| 5 | The contraindication predicate derived three ways (`population`, the body, the regime) | §5.2, §6.2, §10.2 of the first draft | `_regimes`, U4's disagreement clause; §5.2, §16.16 |
| 6 | `record_removal`'s message says [§14a] | `standardise.py:190-194` | Caller-neutral message; §4.2, §16.15 |
| 7 | §16.11 asserted byte-identity with no mechanism to compare against pre-refactor code | `test_standardise.py:1006-1053` precedent | Two digests captured before the refactor; §16.11, §21 task 2a |
| 8 | §16.6's "monkeypatching `rd`" has no patch point | — | Spy on `gcompute`, two calls; §16.6 |
| 9 | Nine paths created by findings 1–8 had no acceptance criterion | the coverage map | §16.15–16.17, `no_contraindicated_frame()` |

Performance: no finding — the second `gcompute` call is two predictions on a 104×11 design per
replicate and was inside the measured 27.6 s. Two `TODOS` amendments (§17 items 9, 10). One prior
learning applied: the `|`-in-cell rule (Stage 7's `|SMD|` header) is why §12 renames `max|coef|` to
`max_abs_coef` and §16.17 asserts cell counts.

**The outside voice** — a second reader with no context of this review, given the full document and
the code — returned fourteen findings. Four repeated the review's own (1, 6, 8, 11). The other ten were
each put to the author and decided:

| OV | Finding | Decision | Folded in |
|---|---|---|---|
| 2 | The `eligible_rd` call predicts an IVT outcome for contraindicated rows and discards it; U2 checks an expression against itself | **Accepted**: second call on `X[eligible]`; U2 on every design, against the stored column | §6.3, §7.3, §16.6 |
| 3 | "gap exactly 0.0" and §16.7's 1e-10 not credible for two Newton runs | **Accepted**: §16.7 at 1e-6 with the argument; the 0.0 kept as measured and re-measured at task 4 | §9.1, §16.7, §21 |
| 4 | δ = −23.6 is `POLR_TOL`'s magnitude, not the data's | **Accepted** | §9.1, §10.2 |
| 5 | §6.4 heading "a third" vs body "a sixth" | **Accepted** | §6.4 |
| 7 | U7 at 1e-14 is below a naive `N·ε` bound | **Accepted**: 1e-12, T10's | §7.3, §11, §16.6 |
| 9 | The hybrid plug-in is the natural estimator for the fourteen distribution keys; ask the PI now | **Accepted, against the review's recommendation**: put to the PI 2026-08-29; **the PI confirmed the model prediction the same day**, DECISION 10 unamended | §18, §17 item 3, §19, §21 |
| 10 | Two "eligible-population effects" on adjacent pages | **Accepted, against the review's recommendation**: [§14a]'s `RD_k` is the estimate; `eligible_rd_k` a diagnostic | §14, §17 item 4, §18 |
| 12 | The byte-identity gate is a skippable test | **Accepted**: run non-skipped, stated as a gate | §21 task 2a |
| 13 | U5's third clause makes a data fact a `SchemaError` | **Accepted, against the review's recommendation**: downgraded to an audit statement; `contrast` proceeds with `share_eligible == 1` | §11, §12, §14, §16.0, §16.10, §16.16 |
| 14 | Status wording; `|β|` 0.0000 display | **Accepted** | Status, §9.2 |

The outside voice found nothing further in `bootstrap.replicates/collect/intervals`, `model.design`'s
handling of the float indicator, or `resample`'s index reset.

---

## 23. What this spec changed elsewhere

**Landed with this document in `37051fc` (2026-08-27):** `implementation_roadmap.md`, Stage 13 — the
`**Spec:**` line; the model under DECISION 10; the dilution identity; the U series; the measured
collapse, drop and cost figures; the arms declined by name. `statistical_analysis_plan.md` [§14b] —
the amendment of 2026-08-27. `TODOS.md` — four new items (§17 items 2, 3, 5, 7), one trigger fired
(the removal-ledger duplication), three examined and not fired (`POLR_MAX_ABS_BETA`'s calibration, the
`code`-field item, the `_role` vocabulary — this stage has no balance table).
`../out/stage0_data_inventory.md` — DECISION 10 (gitignored; task 6).

**Landed with the revision of 2026-08-29:** `TODOS.md` — the row-builder item gains Stage 13's three
renderers as its next instance (§17 item 10); the `code`-field item gains the U1/U8 guard duplication
as a recorded cost of the deferral (§17 item 9). The roadmap's Stage 13 block needs no change: it names
`standardise.gcompute` and the U series, both of which survive the revision.

---

## GSTACK REVIEW REPORT

| Review | Trigger | Why | Runs | Status | Findings |
|--------|---------|-----|------|--------|----------|
| CEO Review | `/plan-ceo-review` | Scope & strategy | 0 | — | — |
| Codex Review | `/codex review` | Independent 2nd opinion | 1 | issues_found (claude subagent, fresh context) | 14 findings: 4 duplicates of the eng review, 10 decided (9 folded, 1 open with the PI) |
| Eng Review | `/plan-eng-review` | Architecture & tests (required) | 1 | clean | 19 issues (9 review + 10 outside voice), 0 critical gaps, scope accepted as-is |
| Design Review | `/plan-design-review` | UI/UX gaps | 0 | — | — |
| DX Review | `/plan-devex-review` | Developer experience gaps | 0 | — | — |

- **CROSS-MODEL:** the outside voice (Claude subagent, cleared context; Codex was not used per the author) agreed with the review on the raise-site scan, `record_removal`'s message, the §16.6 mechanism and the test-count pins, and added ten findings the review missed — the strongest being that the `eligible_rd` call predicted an IVT outcome for contraindicated rows before discarding it (§6.3). Three outside-voice findings were accepted against the review's own recommendation (§22.2, OV 9, 10, 13).
- **VERDICT:** ENG CLEARED — all findings folded; the one question raised to the PI (§18) was answered the same day (model prediction, DECISION 10 confirmed). Ready to implement, tasks 1–6 in order.

NO UNRESOLVED DECISIONS
