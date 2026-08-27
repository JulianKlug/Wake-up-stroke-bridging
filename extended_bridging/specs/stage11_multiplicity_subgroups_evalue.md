# Stage 11 spec — multiplicity, subgroups, E-value

Implements roadmap Stage 11 [§6, §13]. Section references in brackets are to
`statistical_analysis_plan.md`. Numbers and decisions referenced as DECISION *n* are established in
Stage 0 and recorded in `../out/stage0_data_inventory.md`; **DECISION 4 is the one this stage exists
to discharge**. Constants named in `SMALL_CAPS` are Stage 1's and live in `config.py`;
`stage1_config_and_data_contract.md` is their specification. The propensity fit this stage refits
under a second specification is `stage6_propensity_and_weights.md` §11, the balance table it re-roles
`stage7_balance_and_overlap.md` §11, the primary estimator it re-runs
`stage8_primary_outcome_estimator.md` §11, the family partition it corrects within
`stage9_secondary_binary_estimators.md` §14, and the intervals and p-values it reads
`stage10_bootstrap_engine.md` §14 — the last of which is written as a handover *to this document*,
names four things Stage 11 must build and one it must not, and is answered item by item in §14.

**Status.** Written 2026-08-27 against the landed Stages 1-10. **Thirteen probes were run before any
section was drafted, and four of them changed what this document says** — §5.2 (the roadmap's
Accept-when clause was drafted here as needing a correction, and it does not: the two specifications
share an `in_model` population exactly, measured, and the reason is a property of *this workbook* and
not of the design, which is why the clause is derived rather than asserted); §8.3 (the way the
treatment × subgroup interaction column degenerates on this cohort is **not** the constant-column
route this document first specified — it is rank deficiency between two columns that both survive,
which `model.py`'s O6 already catches and already buckets); §8.4 (the separation bound on the
subgroup fit was drafted as a specified-and-unreachable guard and it is **neither** — it fires in
12.0% of the replicates that reach it, and the 12.5% drop rate that follows is the largest this
pipeline has produced); and §4.2 (the audit step names of a second `propensity.fit` over the same
cohort collide **exactly**, which turned the seam's payload from a covariate tuple into a named
record before any of §4 was drafted). §21 names each.

**And an engineering review, also 2026-08-27, which changed four structural things and closed fourteen
gaps** (§22.2). The four: `Subgroups` carries no `bootstrap.Diagnostics` and `Arm` does, because
`Replicate`'s cutpoint fields are scalars and the subgroup body runs two `polr` fits (§3.1, §8.5); a
`model.FitError` from a **point-estimate** subgroup fit is `C.SchemaError` H8, because the
droppable-and-counted argument is about the replicate loop and the point estimate is not in one
(§8.5); the [§11] estimation mask becomes the public `outcome.estimation_population` rather than being
spelled a second time here (§8.5); and this stage's own `C.SchemaError` identifiers are declared as
the **H series**, ten of them, because §15.1 required matching on identifiers that did not exist
(§9.1). The fourteen closed gaps are in §22.2 and the largest was that `bootstrap.replicates` takes a
one-argument body while both bodies here take three — so neither run's call site could be written
from this document as it stood.

**And a literature check, the same day, on the one claim in this document that was about the
literature rather than about this code** (§19c). §7.5 said no published bounding factor exists for a
common odds ratio from a proportional-odds model; it holds, and the check sharpened two things — the
`sqrt(OR)` conversion is cited to **VanderWeele 2017, Epidemiology 28(6):e58** wherever it is stated,
including inside `EValue.approximation`, and its boundary is a documented **15% prevalence** rather
than this document's earlier "a baseline near 0.5", which turns the caveat into a check against
`Primary.cumulative` (§12.2 item 4a). **A cross-model review pass was declined and the reason is a
measurement**: this document is 208 KB, the pass truncates to 30 KB, and 30 KB ends inside §3.1 —
before §4, §5, §6, §7 and §8, which is every decision.

Every number below was produced by running code; none is carried — where a Stage 6, 7, 8 or 10 number
is quoted for comparison it is labelled as theirs and re-measured here. It is the **sole source for
the Stage 11 implementation**: everything the implementer needs is here, and anything not here is not
to be invented.

**Goal.** The four things [§6] and [§13] prescribe on top of the estimates Stages 8, 9 and 10 have
already produced, and no fifth. Benjamini–Hochberg within the secondary and safety families, primary
uncorrected, raw and adjusted p both reported. The E-value for the primary estimate and for its
confidence limit nearest the null. Subgroup estimates on the primary outcome with a treatment ×
subgroup interaction test, labelled hypothesis-generating. And the full-covariate propensity
sensitivity arm [§13, DECISION 4] — one refit of [§7] over `PS_COVARIATES_FULL`, its own weights, its
own balance table, its own bootstrap, reported beside the primary with ESS and worst residual |SMD|.

**Not in scope.** The other **five** [§13] sensitivity rows, which stay deferred and one of which is
deferred *deliberately* (§18, roadmap Stage 11); the standardisation analyses (Stages 12 and 13
[§14a, §14b]); the reporting layer and every label it owes (Stage 14 [§16]); and any recomputation of
an interval, a limit or a p-value that Stage 10 already produced, which Stage 10 §14 forbids by name
and §6.1 restates as a precondition.

**One thing this spec settles that no earlier stage could.** Every stage from 6 onward was written
against **one** propensity specification, and said so: `propensity.fit` takes no covariate list,
`balance._role` reads `PS_COVARIATES` directly, and Stage 6 §6.4 argued that a keyword would make
[§13]'s propensity rows "look like options" when they are *different target populations*. That
argument was correct and it left a promissory note — *"when they are un-deferred they arrive as a
named specification with a [§13] amendment, not as a keyword"* (`propensity.py:462-463`). This is the
stage at which one of them is un-deferred, so it is the first stage that can pay the note, and §4 is
the payment. **The seam is one frozen record, not one tuple**, and §4.2 is why the difference is
load-bearing rather than stylistic.

**And one finding that arrives with the stage rather than being designed into it.** The roadmap's
Accept-when says the primary estimate and the arm *"differ in the propensity specification and in
nothing else"*, and this document's first draft treated that as an overstatement to be corrected:
`in_model` is `complete_cases(df, covariates)`, the four vascular risk factors carry their own
missingness, so the arm's population can only be a subset and its ESS can only fall. Measured: the
two `in_model` masks are **identical at 92**, because all four risk factors are complete on all
**93** cohort rows and the single covariate-incomplete record is incomplete on `core_ml` and
`tmax6_ml` — both of which are in *both* specifications (§5.2). The roadmap is right. But it is right
about **this workbook** and not about the design, so the clause is specified as a number **derived
from the two `Propensity` objects and printed**, never as a sentence asserted — which is the same rule
[§16]'s constant-shift statement is under.

---

## 0. Where Stage 11 sits

```
  cohort.build(...)          -> DataFrame, 93 x 33, one row per patient      [Stage 5 §11]
  propensity.fit(cohort, a)  -> Propensity: e, w, in_model (92), spec        [Stage 6 §11, §4]
  balance.assess(...)        -> Balance: 19 rows, roles off ps.spec          [Stage 7 §11, §4.4]
  outcome.primary(...)       -> Primary:   beta, alpha, rd, cumulative       [Stage 8 §11]
  outcome.secondary(...)     -> Secondary: seven BinaryEstimate, by_family() [Stage 9 §11]
  bootstrap.run(...)         -> Bootstrap: 26 draws, intervals, diagnostics  [Stage 10 §11]

  +--------------------------------------------------------------------------------------+
  |  sensitivity.py                                                                      |
  |                                                                                      |
  |  multiplicity(sec, boot, audit)          -> Multiplicity                    (§6)     |
  |    families = sec.by_family()      3 and 4, NEVER counted off `intervals`   (§6.1)   |
  |    p        = boot.intervals[f"{k}.rd"].p     READ, never recomputed        (§6.1)   |
  |    benjamini_hochberg(p) per family       step-up + descending running min  (§6.2)   |
  |    H3 H4 H5   a None p / a p below THAT KEY'S floor / no `beta` interval    (§6.1)   |
  |                                                                                      |
  |  e_value_primary(est, bal, boot, audit)  -> EValue                          (§7)     |
  |    RR = sqrt(est.odds_ratio)                 the approximation is a FIELD   (§7.1)   |
  |    spans = boot.intervals["beta"].lo <= 0 <= .hi   three routes, one answer (§7.2)   |
  |    H6 H7      a non-finite input / the estimate outside its own interval    (§7.3)   |
  |                                                                                      |
  |  subgroups(cohort, ps, est, audit)       -> Subgroups                       (§8)     |
  |    per S in SUBGROUPS: polr on (A, S, A x S), ps.w, NO [§6] covariate       (§8.2)   |
  |    G8 before `polr`, G9 after it;  H8 when the POINT ESTIMATE raises        (§8.5)   |
  |    replicates(...)  6 keys, 1 run, TWO fits per replicate  -> NO Diagnostics (§8.5)  |
  |    _intervals(draws, _tested_subgroup)      p on `.gamma` ONLY               (§8.5)   |
  |                                                                                      |
  |  full_covariate(cohort, ps, audit)       -> Arm             [§13, DECISION 4] (§5)   |
  |    propensity.fit_full(cohort, a)   -> Propensity(spec=PROPENSITY_FULL)      (§4)    |
  |    balance.assess(...)              -> the SAME 19 rows, four re-roled       (§4.4)  |
  |    outcome.primary(...)             -> UNCHANGED; it reads e, w, in_model    (§5.1)  |
  |    replicates(...)   7 keys, its own body, NOT `run`, ONE fit -> Diagnostics (§5.4)  |
  |    _intervals(draws, lambda k: False)       no p on ANY arm key              (§5.5)  |
  +--------------------------------------------------------------------------------------+

  bootstrap.py's surface goes 5 -> 8:  bucket, collect, diagnostics             (§5.4)
  outcome.py gains ONE public name:    estimation_population                    (§8.5)
  ten `C.SchemaError` identifiers, and they are the H series                     (§9.1)

                  -> Multiplicity, EValue, Subgroups, Arm                       (§3.1)

  the cohort frame, the Propensity, the Primary, the Secondary and the Bootstrap
  come back unchanged                                                            (§0.2)

  audit: load 7 / derive 4 / classify 1 / build 6 / fit 4 / assess 2 / primary 3
         / secondary 3 / bootstrap 1 / sensitivity 15  =  46                     (§10)

  every interval limit, every p-value and every E-value on the workbook is
  DELIBERATELY ABSENT from this document. They are in the gitignored log.        (§1)
```

### 0.1 Why this is one module and not two, and the argument against it is real

`sensitivity.py`, new. It implements [§13] whole — multiplicity, subgroups, sensitivity analyses, the
three clauses of that section in that order — plus [§6]'s E-value, which is prescribed one section
away and reported beside the primary estimate the [§13] arm is compared against.

**The rejected alternative had a genuine import-graph argument and it is recorded rather than
dismissed.** `multiplicity` and `e_value_primary` fit nothing. They are pure functions of numbers
Stages 9 and 10 already produced: a `Secondary`, a `Bootstrap`, a `Primary`. `subgroups` and
`full_covariate` refit — they import `propensity`, `balance` and `model` and they resample. A split
into a `multiplicity.py` that imports no estimator and a `sensitivity.py` that does would be exactly
the move Stage 10 §0.1 made when it refused to put the bootstrap in `outcome.py`, and it would buy
one real thing: two of the four deliverables would live in a module testable with no frame at all.

It is declined for one module per stage-shaped concern, and **the cost is named rather than denied**:
`benjamini_hochberg` and `e_value` are pure arithmetic living in a module that imports `propensity`,
so a test of the step-up drags the whole estimator stack in behind it. That is bearable because both
are public (§3.2) and both are tested as functions over float arrays with no frame constructed
(§15.4, §15.6) — the import is a cost at collection time and not a coupling in the tests. §16 files
the trigger under which the split is revisited: **a second consumer of `benjamini_hochberg`**, which
would mean a family partition existing outside [§13].

The second half of the argument is the one that decided it. [§13] is one section of the plan and its
three clauses share one subject — *what is reported alongside the primary estimate, and under what
qualification*. The multiplicity correction, the hypothesis-generating label and the sensitivity arm
are three answers to that one question, and Stage 14 prints them on one page. A module boundary drawn
between "fits something" and "does not" would split [§13] down an axis [§13] does not have.

### 0.2 What Stage 11 does not touch

**It refits nothing that Stage 10 already refit, and it recomputes no interval.** Stage 10 §14's
prohibition is explicit — *"recompute any p-value or limit, since `Interval` carries the definition
that produced it and a second computation under a different one is §8.2's disagreement re-introduced
downstream"* — and §6.1 turns it into a precondition rather than a convention: `multiplicity` reads
`Interval.p` and raises if one is `None`, because a `None` there is a Stage 10 contract break and not
a sparse family.

**It adds no estimator.** The roadmap says so of the arm — *"it reuses every landed stage and adds no
estimator"* — and it is true of three of the four deliverables literally: `full_covariate` calls
`propensity.fit_full`, `balance.assess` and `outcome.primary`, all landed. The exception is
`subgroups`, and the exception is exact and bounded: **one `model.polr` call over a three-column
design**, which is the landed ordinal fitter with one column added. Nothing here writes a maximiser,
a percentile, a resampler or a p-value definition; `bootstrap.resample`, `bootstrap.replicates`,
`bootstrap.percentile_ci` and `bootstrap.bootstrap_p` are called, and none is reimplemented.

**It does not touch the primary estimate or the [§7] specification.** [§13]'s amendment of 2026-08-24
says this in as many words and DECISION 4 records that it was established *before any sensitivity
estimate existed and before the [§10] interval existed*. §15.12 asserts it by running the whole
Stage 1–10 pipeline with and without this stage and requiring every landed number and the first 31
audit entries to be byte-identical.

**It does not decide the remaining five [§13] rows and does not leave a parameter shaped like them.**
§4's seam is a **registry of two declared specifications**, not a covariate-list argument: a third row
arrives as a third `Specification` with a [§13] amendment behind it, and §15.10 asserts that the
registry has exactly the two.

---

## 1. Deliverables

| Path | Contents |
|---|---|
| `extended_bridging/sensitivity.py` | **new.** [§13] whole and [§6]'s E-value: six public names (§3.2), four returned records (§3.1), thirteen bootstrap estimand keys (§3.3) |
| `extended_bridging/config.py` | **amended.** The `Specification` record, the two declared instances and their registry; one `FAILURE_BUCKETS` token; the three role literals (§13) |
| `extended_bridging/propensity.py` | **amended.** `_fit` private, `fit` and `fit_full` public and named, `Propensity.spec`; eight privates take the record; four audit step names go through `spec.step` (§4.3) |
| `extended_bridging/balance.py` | **amended.** `_role` and `_table` take the specification off the `Propensity`; `assess`'s signature is unchanged; `Balance.spec`; the roles become `Final` constants (§4.4) |
| `extended_bridging/bootstrap.py` | **amended.** R9; `_bucket`, `_collect` and `_diagnostics` become public with no private aliases kept; the public surface goes from five names to eight (§4.5, §5.4) |
| `extended_bridging/model.py` | **amended.** One clause on the design-matrix banner, because "no interaction" is about the [§6] confounder set and reads as a global prohibition (§8.3) |
| `extended_bridging/outcome.py` | **amended in four places and none of them is the estimator.** `primary`'s three audit step names go through `ps.spec.step` so a second call over the same cohort does not collide (§4.4, §10); and the [§11] estimation mask at `:634` is extracted as the public one-line `estimation_population(df, ps)`, which `primary` then calls and `_subgroup_replicate` reads on a drawn frame (§8.5). The estimator is untouched and §15.12 asserts that by number rather than by diff |
| `extended_bridging/tests/test_sensitivity.py` | **new.** §15 |
| `extended_bridging/tests/fixtures_stage11.py` | **new.** §15.0 |
| `extended_bridging/tests/test_propensity.py`, `test_balance.py`, `test_bootstrap.py`, `test_config.py`, `test_outcome.py` | **amended.** The seam's wrong-wiring tests, the ledger and surface counts that move, and `estimation_population`'s equivalence to the mask `primary` binds (§15) |
| `extended_bridging/implementation_roadmap.md` | **amended.** Stage 11 gains its `**Spec:**` line; three corrections and two additions (§22) |
| `TODOS.md` | **amended.** Two items close, one is rewritten because this stage falsifies its trigger, **eight** are new (§16) |

`out/logs/audit_<label>.md` is an **output, not a deliverable**. For this stage it is the only place
any adjusted p, any E-value, any subgroup odds ratio and any arm limit exists. It is gitignored.

**Nothing under `specs/` may quote a case identifier, and nothing here does.** Stage 8 §4.3 set that
rule and Stage 10 §4.5 tightened it to interval limits and p-values. This document tightens it once
more, to **E-values and adjusted p-values**, for the reason that makes those two worse than a limit:
an E-value is a single number a reader will quote as the robustness of the finding, and an adjusted p
is the number a reader will read as significance. Both belong to the log and to the manuscript, and
neither belongs to a document in version control that describes how they are computed.

---

## 2. Environment

Unchanged from Stage 10: `uv`, Python 3.12.12, pandas 2.3.3, numpy 1.26.4, statsmodels 0.14.6.
`scipy` is test-only by policy and **`statsmodels` becomes test-only in practice at this stage** —
§6.4 is the one place in this pipeline where a shipped module could plausibly have imported it, and
does not. **No dependency is added.**

```bash
cd extended_bridging && uv sync && uv run pytest -v
```

**Cost, measured (§21).** This stage adds two bootstrap runs to the one Stage 10 already pays for.
Neither is a `run`: both are `bootstrap.replicates` over a body that refits far less than [§10]'s
does, and §5.4 and §8.5 are the arguments for that shape rather than the cost being its motivation.

The numbers are in §21 and the conclusion they support is in §16: **the trigger `TODOS.md` records
for parallelising the bootstrap — *"the [§13] sensitivity suite, which is five more `run` calls and
turns 145 s into something over ten minutes"* — is falsified twice over.** There is one arm and not
five, and it is not a `run`.
---

## 3. Module shape

### 3.1 What the stage returns

Four frozen dataclasses, one per [§13]/[§6] deliverable, none of them carrying a verdict. Stage 8 §3's
rule — *"there is no `favours_bridging()` and no verdict of any kind"* — binds harder here than
anywhere before it, because three of the four records are the ones a reader most wants a verdict from:
an adjusted p, an E-value, and an interaction test.

```python
@dataclass(frozen=True)
class FamilyCorrection:
    """One [§13] family's Benjamini-Hochberg result, and the counters that explain its denominator.

    **`m_declared` and `m_used` are BOTH fields and neither is derivable from the other**, which is
    the whole of §6.5: `m_declared` is the family's size in the [§5] registry and `m_used` is the
    number of p-values that entered the procedure. A reader given only the second cannot tell a
    complete family from a thinned one, and BH's scaling factor IS `m_used`, so the two numbers
    together are what makes the adjusted p checkable by hand.

    `absent` names the outcomes that had no p, with the reason, because "we tested this outcome and
    it had no interval" and "this outcome is not in this family" are the two things §6.5's disclosure
    has to tell apart. It is empty on this workbook and the field exists anyway (§6.5).

    **The reason string's shape is specified and is not free text**, because §6.5 requires the
    rendered row to carry the surviving-draw count and the floor it fell below, and a formatter
    cannot recover either from prose:

        f"no interval: {len(boot.draws[key].draws)} surviving draw(s) of "
        f"{boot.draws[key].n_attempted} against a floor of {C.ci_min_draws()}"

    Both numbers come off the `Draws` that Stage 10 kept when it withheld the `Interval`
    (`bootstrap.py:1099-1103`), so the disclosure is a read and never a recomputation.

    There is no `significant` field and no threshold anywhere in this record. [§13] prescribes a
    correction, not a decision rule, and [§10] closes with "estimation, not testing, is the reportable
    output".
    """

    family: str                         # "secondary" | "safety" — Stage 9's `by_family()` key
    m_declared: int                     # |family| in the [§5] registry — 3 and 4
    m_used: int                         # p-values that entered the step-up — §6.5
    raw: dict[str, float]               # outcome -> the p Stage 10 produced, READ not recomputed
    adjusted: dict[str, float]          # outcome -> the BH-adjusted p; same keys as `raw`
    absent: tuple[tuple[str, str], ...]  # (outcome, why) for members with no p — §6.5
    descriptive_only: bool              # [§10]'s label; True for "safety" — §6.6


@dataclass(frozen=True)
class Multiplicity:
    """The two families, corrected separately, and the primary's exemption stated as a field.

    **`primary_uncorrected` is a field holding the primary's raw p and not an omission.** [§13]'s
    first sentence is "Primary outcome uncorrected", and an implementation that satisfies it by
    leaving the primary out of the record satisfies it in a way Stage 14 cannot print. The rule is
    that the primary appears beside the corrected families carrying the statement that it is
    uncorrected — so the exemption is visible where the correction is.
    """

    families: dict[str, FamilyCorrection]   # keyed as `Secondary.by_family()` is
    primary_uncorrected: float              # `boot.intervals["beta"].p` — §6.1


@dataclass(frozen=True)
class EValue:
    """The [§6] E-value for the primary estimate and for the limit nearest the null.

    **`approximation` is a FIELD carrying the string, not a literal at the print site.** The roadmap
    requires that the `RR ~ sqrt(OR)` conversion be *stated* wherever the E-value is reported, and
    `Interval.method`'s docstring gives the design rule this follows: which definition produced a
    number is part of what the number is. A field makes printing the number without the statement a
    thing Stage 14 has to do deliberately rather than by forgetting.

    **`e_limit` is `None` only when there is no interval at all** — an estimand below
    `C.ci_min_draws()`, which is unreachable on v7 (§7.4). It is NOT `None` when the interval spans
    the null: that case is `1.0`, exactly, and §7.3 is why the difference matters.

    **`n_draws` is the field that makes the `None` readable.** §7.4 requires that when no interval
    exists the draw count be recorded rather than the E-value silently omitted, and a record whose
    only statement is `e_limit is None` cannot tell "the primary lost too many replicates" from "the
    caller passed no interval". It is `boot.intervals["beta"].n_draws` when an interval exists and
    `len(boot.draws["beta"].draws)` when one does not, so it is populated on both branches.

    `spans_null` is a field and not a method for the same reason `Interval.method` is a field: three
    different routes decide it and §7.2 asserts they agree, so which answer this record carries must
    be inspectable rather than recomputed by a caller that might take a fourth route.
    """

    measure: str            # "common odds ratio [§8]" — what OR is, said once
    odds_ratio: float       # est.odds_ratio, READ off the point estimate — §7.1
    risk_ratio: float       # sqrt(odds_ratio) — the approximation, applied
    e_point: float          # the point estimate's E-value
    spans_null: bool        # decided on intervals["beta"] against 0 — §7.2
    limit: float | None     # the limit nearest the null, on the BETA scale; None when spanning
    e_limit: float | None   # 1.0 exactly when spans_null; None only when no interval exists
    n_draws: int            # the draws behind the limit — populated on BOTH branches
    approximation: str      # "RR ~ sqrt(OR) for a common outcome [VanderWeele 2017,
                            #  Epidemiology 28(6):e58; roadmap Stage 11]" — §7.5 cites the source
                            #  and §19c is why the citation is in the string and not only in a spec


@dataclass(frozen=True)
class SubgroupEstimate:
    """One [§13] subgroup: the two level-specific common odds ratios and the interaction.

    **`hypothesis_generating` is `True` and it is a field rather than a docstring.** [§13] labels
    these estimates and Stage 14 has to print the label; a boolean that is structurally always True
    looks redundant until one asks where the label would otherwise live, which is a sentence in a
    reporting module that nothing asserts.

    **`n_interaction_tests` is carried on every subgroup and is the count over the WHOLE stage**, not
    this subgroup's own. §8.8 declines to put these p-values into [§13]'s Benjamini-Hochberg and owes
    the reader the denominator instead; a per-subgroup field would print `1` beside each of two tests.
    **It is `len(C.SUBGROUPS)` and never the literal 2**, which is §3.3's rule for the six estimand
    keys applied to the one number that would otherwise contradict them: a third subgroup restored to
    the registry gains its keys and its interval automatically and must not go on printing that two
    interaction tests were performed.

    `or_level` and `or_ratio` are `exp` of coefficients from ONE fit (§8.2), so they cannot disagree
    about the model. There is no field holding the primary's odds ratio for comparison, and §8.7 is
    why: the two level ORs need not bracket it, a reader will expect them to, and a field inviting
    the comparison is the wrong place to answer that.
    """

    subgroup: str                       # a C.SUBGROUPS key
    n: int                              # records in the fit — the [§11] denominator
    n_by_level_arm: dict[str, int]      # "S=0,A=1" -> count; the four cells §8.3 turns on
    or_level: dict[int, float]          # level -> exp(beta) and exp(beta + gamma) — §8.2
    or_ratio: float                     # exp(gamma) — the interaction, on the OR scale
    gamma: float                        # the interaction coefficient; null 0 — §8.5
    cumulative: dict[int, dict[int, dict[int, float]]]   # level -> threshold -> arm -> P(Y<=k) — §8.6
    hypothesis_generating: bool         # True — [§13]
    n_interaction_tests: int            # 2 — the multiplicity [§13] does not address — §8.8


@dataclass(frozen=True)
class Subgroups:
    """The [§13] subgroups, their draws and their intervals. Keyed as C.SUBGROUPS is.

    `draws` and `intervals` are Stage 10's own dataclasses, unmodified and not re-declared: an
    interval produced here is the same kind of object as an interval produced there, carrying the
    same `method`, the same `n_draws` and the same floor rule. A second `Interval` type would be a
    second percentile definition one edit away (§8.5).

    **`intervals` may hold fewer keys than `draws`**, which is Stage 10 §3.3's rule inherited rather
    than restated: an estimand below `C.ci_min_draws()` keeps its `Draws` and gets no `Interval`.

    **There is NO `diagnostics` field here and `Arm` has one, which is a difference in the bodies
    and not an oversight.** `bootstrap.Diagnostics` tallies `n_alpha` and `polr_iterations` as
    scalars per replicate (`bootstrap.py:259-260`) because [§10]'s body and the arm's each run
    exactly ONE `polr` fit. The subgroup body runs TWO (§8.5), so a single scalar cannot say which
    fit it came from and a tally mixing them is a distribution of nothing. Everything
    `subgroup_replicates` prints therefore comes off `Draws` — `n_attempted`, `len(draws)` and the
    bucket counts, which is where §8.4's 12.5% already lives — and `_subgroup_replicate` sets
    `n_alpha` and `polr_iterations` to `None` deliberately (§8.5).
    """

    estimates: dict[str, SubgroupEstimate]
    seed: int
    n_boot: int
    draws: dict[str, bootstrap.Draws]        # six keys — §3.3; the counters live HERE, not in a Diagnostics
    intervals: dict[str, bootstrap.Interval]  # the same keys; p on `.gamma` only — §8.5


@dataclass(frozen=True)
class Arm:
    """The [§13, DECISION 4] full-covariate propensity sensitivity analysis.

    **It carries its own `Propensity`, its own `Balance` and its own `Primary`, and it carries the
    PRIMARY's alongside them.** Every number in this record exists to be read beside another number,
    and the roadmap's Accept-when is a statement about a pair: the arm's `in_model`, ESS and worst
    residual |SMD| "are reported beside the primary's". A record holding only the arm's half makes
    that pairing Stage 14's job, done from two objects it has to be trusted to match up.

    **`__post_init__` raises unless `ps.spec is C.PROPENSITY_FULL` and `reference.spec is
    C.PROPENSITY_PRIMARY`**, by identity and not by equality (§4.5). This is the last of the four
    places the arm's numbers can be labelled as the primary's, and it is the one a caller reaches by
    building the record by hand.

    `differs_only_in_specification` is COMPUTED from the two `in_model` masks and never asserted.
    On this workbook it is True and §5.2 measures why; on a workbook where any of the four vascular
    risk factors is missing on any row it would be False, and the clause a reader sees must be the
    one the data supports.

    **`diagnostics` is `bootstrap.Diagnostics`, unmodified and not re-declared**, and it is here
    because the arm's body runs exactly ONE `polr` fit per replicate — so `n_alpha` and
    `polr_iterations` mean what they mean in [§10]'s own run, and §5.4's cutpoint counts are a
    counter in the log rather than a number this document measured once. `Subgroups` has no such
    field and §3.1's `Subgroups` docstring is why. `max_abs_beta` and `or_corrected` come back empty:
    the arm produces no binary estimate and no `m_a(X)` (§5.4), which is `Diagnostics`' own
    "AUGMENTED ONLY" rule satisfied by there being nothing augmented.
    """

    ps: propensity.Propensity           # spec is C.PROPENSITY_FULL — asserted, §4.5
    balance: balance.Balance            # the same 19 rows, four re-roled — §4.4
    estimate: outcome.Primary           # `outcome.primary`, UNCHANGED — §5.1
    reference: propensity.Propensity    # the [§7] fit, for the pair the Accept-when is about
    seed: int
    n_boot: int
    draws: dict[str, bootstrap.Draws]         # seven keys — §3.3
    intervals: dict[str, bootstrap.Interval]  # the same keys; `p is None` on all seven — §5.5
    diagnostics: bootstrap.Diagnostics         # ONE polr fit per replicate, so the tallies mean it
    differs_only_in_specification: bool        # COMPUTED from the two masks — §5.2
```

### 3.2 Public surface, and it is six names

`benjamini_hochberg`, `e_value`, `multiplicity`, `e_value_primary`, `subgroups`, `full_covariate`.
**This count is stated once, here.** §15.11 asserts it by module scan and §1 cites this sentence.

The first two are pure functions over floats and they are public for the reason `propensity.ess`,
`balance.smd`, `outcome.weighted_proportion` and `bootstrap.percentile_ci` are: each is a
**definition** rather than a step, and a definition with one caller is still a definition somebody
will want to check against an oracle without constructing a frame. `benjamini_hochberg` takes an
array and returns an array; `e_value` takes a risk ratio and returns a float. Neither names an
outcome, a family, a covariate or a centre, and §15.11 scans for that.

The six privates that matter: `_arm_replicate`, `_subgroup_fit`, `_subgroup_replicate`,
`_assert_arm_inputs`, `_tested_subgroup` and `_intervals`. There is no `_bucket`, no `_collect` and
no `_diagnostics` here, because §5.4 makes Stage 10's three public (§13) rather than writing second
ones.

**`_intervals` exists so the draws-to-intervals loop is written once and not twice.** Both runs turn
a `dict[str, Draws]` into a `dict[str, Interval]` under the same three rules — the
`C.ci_min_draws()` floor, `bootstrap.percentile_ci` at `C.CI_LEVEL`, and `bootstrap_p` where and
only where the stage prescribes a test — and the two runs differ in exactly one of them:

```python
def _intervals(draws: dict[str, bootstrap.Draws],
               tested: Callable[[str], bool]) -> dict[str, bootstrap.Interval]:
    """§8.3's rule of Stage 10, applied twice from one place. The floor withholds the Interval and
    keeps the Draws; `tested` is the ONLY thing that differs between the arm and the subgroups.

    A second copy of this loop is a second place the floor could be spelled, a second place
    `C.CI_LEVEL` could be defaulted, and the [§13] arm and the [§13] subgroups would then be one
    edit away from disagreeing about what a percentile interval is. `bootstrap.run` keeps its own
    copy (`bootstrap.py:1105-1111`) and is not refactored onto this one: Stage 10's DoD asserts its
    numbers and this stage refits nothing Stage 10 refit (§0.2).
    """
    out: dict[str, bootstrap.Interval] = {}
    for key, d in draws.items():
        if len(d.draws) < C.ci_min_draws():                  # Stage 10 §8.3 — Draws kept, no Interval
            continue
        lo, hi = bootstrap.percentile_ci(d.draws, C.CI_LEVEL)
        p = bootstrap.bootstrap_p(d.draws) if tested(key) else None
        out[key] = bootstrap.Interval(lo, hi, C.CI_LEVEL, C.PERCENTILE_METHOD, len(d.draws), p)
    return out
```

The arm passes `lambda key: False` — §5.5's rule as an argument rather than as a convention — and the
subgroups pass `_tested_subgroup`, which is `key.endswith(".gamma")` and nothing else (§8.5).

**`_tested` is NOT reused for either**, and §5.5 is why: `bootstrap._tested("beta")` is `True`
(`bootstrap.py:456`), so an implementation that reached for it would emit a p on the arm's `beta`
that no prespecified section asks for.

### 3.3 The estimand keys, and there are thirteen

Two disjoint sets over two independent `bootstrap.replicates` calls. They are never merged into one
dict and §5.4 says why that prohibition needs stating.

```
  the arm (§5.4) — seven, and the key NAMES ARE THE PRIMARY'S:
    beta                    the [§13] common odds ratio's coefficient
    rd_0 rd_1 rd_2 rd_3 rd_4 rd_5      the six cumulative RD_k, off est_full.rd

  the subgroups (§8.5) — six, two subgroups x three quantities:
    unknown_onset.beta                 log OR at S = 0   (witnessed onset)
    unknown_onset.gamma                the interaction; THE ONLY TESTED KEY
    unknown_onset.beta_plus_gamma      log OR at S = 1   (unwitnessed or wake-up)
    core_above_median.beta             log OR at S = 0   (core at or below the median)
    core_above_median.gamma            the interaction; THE ONLY TESTED KEY
    core_above_median.beta_plus_gamma  log OR at S = 1   (core above the median)

  tested: 2 of 13.  `.gamma` and nothing else (§5.5, §8.5).
```

**The arm's keys are deliberately the primary's, and that is a decision with a cost.** The two tables
are compared key for key — `beta` against `beta`, `rd_2` against `rd_2` — which is what the roadmap's
"reported beside the primary ones" means when the reader is a person and not a formatter. The cost is
that `{**boot.draws, **arm.draws}` is a silent overwrite: seven keys collide, nothing in `Draws`
notices, and the result is a dict of 26 keys in which seven describe a different specification.
**The two dicts are never merged and §15.9 asserts the prohibition by constructing the merge and
showing what it costs**, because a rule stated in prose about a dict operation is a rule nobody runs.

The six subgroup keys are computed from `C.SUBGROUPS` and never listed, which is Stage 9 §4.1's move
two stages on: a third subgroup restored to the registry gains its three keys, its interval and its
p in the same edit, and cannot fail to.

**The arm's six `rd_k` are computed from `est_full.rd` and NOT from `C.MRS_THRESHOLDS`**, which is
`bootstrap._estimand_keys`' own argument taken rather than re-derived: *"the six `RD_k` keys live on
the POINT ESTIMATE and not on the configuration: `Primary.rd` is `dict[int, float]` and its keys are
what the primary actually fitted"* (`bootstrap.py:431-441`). The arm therefore keys off its own
`Primary` — the one `outcome.primary` returned under `C.PROPENSITY_FULL` — so if a future cohort
fitted a different set of thresholds under one specification than the other, the mismatch surfaces
as a key that is absent rather than as a `Draws` reconciliation that cannot be explained.

### 3.4 The numerical facts this stage turns on

Every one measured in §21; the ones carried from an earlier stage are labelled as that stage's and
re-measured here.

```
  92 / 92        in_model under [§7] and under [§13]'s full covariate set — IDENTICAL,
                 and the roadmap's Accept-when is therefore TRUE on v7           [§5.2]
  93 / 0         cohort rows / rows the four vascular risk factors are missing on [§5.2]
  1 record       covariate-incomplete, and it is incomplete on core_ml and tmax6_ml,
                 which are in BOTH specifications                                 [§5.2]
  12 -> 16       propensity design width, [§7] -> [§13]                           [§5.1]
  5 -> 2         balance rows reaching SMD_THRESHOLD, [§7] -> [§13]                [§5.3]
  0.380 -> 0.052 diabetes |SMD| after weighting, [§7] -> [§13]                     [§5.3]
  0.211 -> 0.275 center = HUG |SMD| after weighting, [§7] -> [§13] — it gets WORSE [§5.3]
  1986 / 2000    the arm's surviving replicates; 13 degenerate_design, 1 nonconvergence,
                 against the primary's 1998 measured by Stage 10                   [§5.4]
  69.3 s         the arm's whole bootstrap, against Stage 10's 130-145 s for `run`  [§2, §5.4]
  9 / 2000       subgroup replicates dropped by O6 rank deficiency, ALL of them the
                 replicates that drew no treated patient at witnessed onset        [§8.3]
  0 / 2000       subgroup replicates dropped by the constant-interaction route the
                 first draft of this document specified as the reachable one       [§8.3]
  5 / 41         treated patients at witnessed onset, on the point estimate, and this
                 one number is why §8.3 and §8.4 both exist                        [§8.1]
  3 and 4        the two [§13] family sizes, off `Secondary.by_family()`            [§6.1]
  7 of 26        Stage 10 keys carrying a p — so counting `intervals` gives 26      [§6.1]
  1000           distinct p-values `bootstrap_p` can return at B = 1998, which is
                 why §6.2's tie rule is ordinary rather than theoretical            [§6.2]
  40             `ci_min_draws()`, unchanged and never re-derived                   [§6.5]

  every adjusted p, every E-value, every subgroup odds ratio and every arm limit on
  the workbook is DELIBERATELY ABSENT from this document.                          [§1]
```

---

## 4. The covariate seam, and it is the one design decision the roadmap left here

The roadmap says so directly: *"this arm needs a seam for the covariate set, and that seam is the one
design decision here."* Everything else in the arm is a call to a landed function.

### 4.1 The three candidates, and the one that cannot be taken is the obvious one

**A `covariates=` keyword on `propensity.fit`.** This is what the pilot does —
`pilots/analysis.py:953` takes `covars or C.PS_COVARIATES` — and it is what Stage 6 §6.4 declined in
advance, with `propensity.fit`'s own docstring carrying the refusal: *"the [§13]
propensity-specification rows are different target populations and a parameter makes them look like
options."*

Two things make it not merely inadvisable but unavailable. First, `test_propensity.py:519-524`
asserts by `inspect.signature` that `fit` takes exactly `["df", "audit"]` **and that no parameter has
a default**, with the comment *"so adding a defaulted keyword fails the test rather than passing
silently"* — the keyword requires deleting the guard whose first customer this arm is. Second, and
worse, its failure mode is silent in the direction that matters: under `covars or C.PS_COVARIATES` a
call site that forgets the keyword produces the **primary** specification, and every number
downstream reads as the arm's.

**A second named function in `sensitivity.py`, duplicating `fit`'s body.** Declined against the
roadmap's own sentence — *"it reuses every landed stage and adds no estimator"* — and against what
duplication costs here specifically: F1 and F2 (`propensity.py:139-163`), F5
(`propensity.py:166-182`), the nan-sentinel construction (`propensity.py:485-488`),
`_record_exclusion`'s identifier guard (`propensity.py:400-432`) and the per-arm ESS loop would each
exist twice. It would also put a **second `model.firth` call site** in the repository, which is what
makes roadmap invariant 5 — no estimate is ever produced by a fallback estimator — verifiable by the
scan that currently proves it (`test_model.py:816-833`) rather than by reading two functions.

**A named record on the `Propensity`, and two named entry points.** Taken. §4.2.

### 4.2 `Specification`, and why the payload is a record and not a tuple

The obvious form of the winner is a `covariates: tuple[str, ...]` field on `Propensity`. It is wrong,
and the way it is wrong is measurable: the seam has four jobs and a tuple does one.

| job | a bare tuple | a named record |
|---|---|---|
| `balance._role` needs the adjusted-for set | yes | yes |
| the audit step names must not collide (§4.4) | **no** — `assess` would have to compute a *name* from a *value*, and any third specification falls back to the primary's name silently | yes, the name is declared |
| Stage 14 needs a printable label bound to these numbers | no | yes |
| a third [§13] row cannot appear without a [§13] amendment | no | yes — the registry is pinned (§15.10) |

The second row is decisive and §4.4 measures it. It is also the house rule read literally: *a
prespecified specification is a NAMED specification*, and a tuple of column names is not a name, it is
a value. The precedent is `config.OUTCOMES` and `config.outcome_model_covariates`
(`config.py:637-655`), whose own argument is that *"an outcome can only diverge by being named."*

```python
@dataclass(frozen=True)
class Specification:
    """One prespecified propensity specification. There are two and there is a registry (§15.10).

    `suffix` is what keeps two fits over the same cohort apart in ONE audit log, and it is a
    declared string rather than a name computed from `covariates` for the reason `Source` is a
    declared object in `data.py`: a label derived from a value is a label that changes when the
    value does, silently, and this one is what `Audit.entry` matches on.

    **The primary's suffix is the EMPTY STRING, and that is load-bearing rather than tidy.** Every
    step name Stages 6 and 7 already record is `spec.step(base)` with `suffix = ""`, so it is
    byte-identical and the four landed ledger assertions do not move on a run without the arm.

    `covariates` is the REQUESTED list and never the surviving one. §4.3.
    """

    label: str                    # what a table prints beside the numbers
    covariates: tuple[str, ...]   # the REQUESTED list — §4.3
    suffix: str                   # the audit step suffix; "" for [§7] — §4.4
    sap: str                      # the section that prescribes it

    def step(self, base: str) -> str:
        return base + self.suffix
```

**`covariates` is the requested list and never the post-`design` survivors**, for three reasons in
this order:

1. **`in_model` is defined by the requested list.** `model.complete_cases(df, spec.covariates)` runs
   *before* `model.design`, so the requested list is what decided the [§11] denominator. A field
   holding the survivors would describe a different set from the one the denominator came from, and
   the object would be internally inconsistent.
2. **`dropped` already exists and is a different granularity.** It is design-*column*-level —
   measured, `center_USZ` on both specifications — and `_design_table` (`propensity.py:242-275`)
   already renders it per declared level. Collapsing the two is a category error.
3. **A covariate dropped as constant is still adjusted for.** A survivors-based field would re-role
   such a covariate from `propensity model` to `excluded [§6]` in the balance table — a false claim,
   and a violation of Stage 7 §4.1's *"not one row moves"* that nothing would assert.

### 4.3 What the seam costs in `propensity.py`, and the claim it spends

```python
def _fit(df: pd.DataFrame, spec: C.Specification, audit: Audit) -> Propensity: ...   # private

def fit(df: pd.DataFrame, audit: Audit) -> Propensity:        # [§7] — the estimand
    return _fit(df, C.PROPENSITY_PRIMARY, audit)

def fit_full(df: pd.DataFrame, audit: Audit) -> Propensity:   # [§13, DECISION 4]
    return _fit(df, C.PROPENSITY_FULL, audit)
```

Both public entry points keep `(df, audit)` and no defaults, so `test_propensity.py:519-524` survives
**verbatim** and is parametrised over the two names rather than rewritten. That test is the door, and
the seam is built so as not to take it off its hinges.

`Propensity` gains `spec: C.Specification`, **appended and with no default**. A default of
`C.PROPENSITY_PRIMARY` would keep the direct construction sites in `test_balance.py` and
`test_outcome.py` green, and that is exactly the hazard: an arm-derived `Propensity` that forgot the
field would role-label as the primary's. Without a default those sites raise `TypeError` — loud, and
each a one-line fix. `dataclasses.replace` (the pattern at `test_bootstrap.py:192-194`) carries the
field correctly and stays the safe way to bend a `Propensity`.

Eight privates take the record because each reads `C.PS_COVARIATES` or asserts "[§6]" today:
`_completeness_detail` (`:200`), `_completeness_table` (`:214`, hardcode at `:225`), `_design_detail`
(`:229`), `_design_table` (`:242`), `_fit_detail` (`:278`), `_weights_detail` (`:302`),
`_weighted_population_table` (`:341`), `_record_exclusion` (`:400`). Three of their `detail` strings
are false in the arm and must name the specification rather than say "the [§6] covariates" — and one
of them is false in a way that has no referent at all: `_weights_detail` says *"conditional on this
specification"*, and with two fits in one log "this" points at nothing.

**The claim this spends, and it is Stage 10's.** Stage 10 §5's central result is that `propensity.py`
has a **zero-line diff** — the resampler was fixed on its own side precisely so that Stage 6's guard
kept its meaning — and Stage 10's Definition of done item 2 checks it. That claim is about *Stage
10's commit* and stays historically true; this document must say so, because a reader who runs
Stage 10's DoD item 2 as a live check after this stage lands will watch it fail and conclude
something is broken.

**And one thing the seam pays for.** `TODOS.md` carries *"Fix the stale `data.py:255` citation in
`propensity.py`"* with the trigger *"the next commit that edits `propensity.py` for any reason"*.
The citation is at `propensity.py:413` and the correct line is **`data.py:271`** — verified. The
trigger is met and the item closes (§16, §22).

### 4.4 `balance._role` reads the specification off the `Propensity`, and `assess` keeps its signature

```python
def _role(name: str, spec: C.Specification) -> str: ...
def _table(sub, a, w, names, spec, sds=None) -> tuple[CovariateBalance, ...]: ...
def assess(df, ps: propensity.Propensity, audit: Audit) -> Balance: ...   # UNCHANGED
```

`assess` reads `ps.spec` and passes it down. **This is not a workaround and Stage 7 promised it in
these words** (`stage7_balance_and_overlap.md:1458-1462`): *"When those specifications arrive they
call `assess` once each, with the `Propensity` from their own named specification, and the `role`
column moves without the table changing shape (§4.1)."*

The branch order is unchanged and it is what produces the arm's table: the four vascular risk factors
hit `name in spec.covariates` first and read `propensity model`; `penumbra_ml` still reads
`excluded [§6]`; the `negative control` branch is dead in the arm and live in the primary. Measured
(§21): **19 rows in both, identical labels in identical order, roles 14 / 4 / 1 in the primary and
18 / 0 / 1 in the arm, and exactly four rows differ** — `hypertension`, `hyperlipidemia`, `diabetes`,
`smoking`. That is the roadmap's *"the same 19 rows with four of them re-roled"* as a measurement.

Two further `balance.py` changes, each with its own reason:

- **The three role strings become module-level `Final` constants**, as `_REPORTED`, `_STRUCTURAL` and
  `_NO_CONTRAST` already are (`balance.py:367-369`). Stage 14 now holds two `Balance` objects whose
  only structural difference is the role column, so it matches on `row.role` **by value**;
  `CentreOverlap`'s docstring (`balance.py:93-105`) makes exactly this argument for `status` — *"Stage
  14 reads this column BY VALUE before it draws anything"* — and it transfers unchanged.
- **`Balance` gains `spec`**, for `Interval.method`'s stated reason. Two balance tables will sit side
  by side and the role column is the only thing telling them apart today.

**Measured, and this is the probe that changed §4.2 (§21).** A second `propensity.fit` over the same
cohort records the same four `model` step names — `covariate_completeness`, `design_matrix`,
`propensity_fit`, `overlap_weights` — and `assess` the same two, `balance_smd` and
`overlap_by_centre`. `Audit.entry` is **first-match** (`data.py:273-275`), so with reused names every
`audit.entry("model", ...)` read in `test_propensity.py`, `test_balance.py` and Stage 14 returns the
**primary's** entry while the rendered log looks complete — and every existing test stays green,
because on a run without the arm they *are* the same entry. `spec.step` is what defuses it, and
`test_config.py` asserts that exactly one declared specification has an empty suffix and that all
suffixes are distinct, so a third specification cannot collide by being added.

### 4.5 R9, and the wrong wiring nothing today would notice

`bootstrap._replicate` calls `propensity.fit(draw, audit)` unconditionally (`bootstrap.py:578`). So:

```python
boot_full = bootstrap.run(cohort, ps_full, est_full, sec_full, audit)     # <- silently wrong
```

returns twenty-six intervals whose **point estimates are the arm's** and whose **2000 replicates are
the primary's specification**. Every number is finite. `Draws.__post_init__` reconciles. R3 passes —
the index is aligned. R5 passes — the `Secondary` is complete. Nothing raises, and the arm's interval
is a lie about the arm.

**R9, new, in `_assert_run_inputs` phase 1:** `C.SchemaError` unless `ps.spec is
C.PROPENSITY_PRIMARY`. Identity and not equality, because the registry holds exactly two instances
and identity is the comparison a typo cannot satisfy. §15.2 asserts it with a mutation companion —
the call succeeds with R9 removed — because a guard whose absence is never demonstrated is a guard
nobody can price.

`Arm.__post_init__` is the same check from the other side, on the record rather than on the call
(§3.1), and it is the one a caller reaches by building the record by hand.

---

## 5. The full-covariate arm [§13, DECISION 4]

### 5.1 What it refits, and what it does not

Three landed calls and one loop:

```
  ps_full = propensity.fit_full(cohort, audit)        [Stage 6 §11, through §4's seam]
  bal_full = balance.assess(cohort, ps_full, audit)   [Stage 7 §11, roles off ps.spec]
  est_full = outcome.primary(cohort, ps_full, audit)  [Stage 8 §11, UNCHANGED]
  draws, intervals = replicates(...)                  [Stage 10 §3.2, own body — §5.4]
```

**`outcome.primary` is not amended and the roadmap's reason is exactly right.** It reads `e`, `w` and
`in_model` off the `Propensity` and nothing else, and its own design matrix is `model.design(sub,
(C.TREATMENT,))` — treatment alone (`outcome.py:639`). So the outcome side of the arm is not a
different model; it is the same model on different weights. §15.12 asserts `outcome.py` has a
zero-line diff at this stage, which is the roadmap's claim tested rather than quoted.

**`C.OUTCOME_COVARIATES` does not move, and that is why the arm is one change and not two.** It is
bound to the `PS_COVARIATES` object at import (`config.py:561`), so a reader may reasonably expect
the [§13] specification to widen the outcome regressions too. It does not: the arm changes the
propensity specification and the estimand's weighting, and [§13] asks for nothing else. Stage 9's
`m_a(X)` models are not refit at all, because the arm produces no binary estimate (§5.4).

Measured (§21): the propensity design goes from **12 columns to 16**; `dropped` is `center_USZ` in
both, so the four added covariates all survive as columns.

### 5.2 The population, and the roadmap's Accept-when is true on v7 — measured, not reasoned

The roadmap requires the two estimates to appear together *"with the statement that they differ in
the propensity specification and in nothing else."* **This document's first draft treated that as an
overstatement and specified a correction to it.** The argument was straightforward: `in_model` is
`model.complete_cases(df, spec.covariates)`, the four vascular risk factors carry their own
missingness, so `in_model_full` can only be a subset of `in_model` and the arm's ESS and every row's
`n` can only fall.

Measured, on the cohort:

```
  93        cohort rows
  92 / 92   in_model under [§7] and under [§13] — IDENTICAL masks, not merely equal counts
  0         rows lost to the four vascular risk factors
  0 / 0 / 0 / 0   missing values on hypertension / hyperlipidemia / diabetes / smoking
  1         covariate-incomplete record, incomplete on core_ml and tmax6_ml — BOTH of which
            are in BOTH specifications, so it is excluded from both for the same reason
```

The roadmap is right. **It is right about this workbook and not about the design**, and the
difference is the whole of the specification here: on a workbook where any one of the four risk
factors is missing on any row that clears the [§6] mask, the two populations diverge and the clause
becomes false while nothing in the code changes.

So `differs_only_in_specification` is a **field computed from the two masks** (§3.1) and the printed
sentence is selected from it, with **both `in_model` counts printed either way**. It is never a
sentence asserted in a formatter. This is the rule [§16]'s constant-shift statement is already
under — Stage 10 §12.3, *"computed, not written"* — applied to the one clause the roadmap states as
prose.

The ESS moves and the direction is not the one the completeness argument predicted, because the
populations are identical and only the weights differ. Measured (§21): control **29.199 -> 29.308**,
treated **30.339 -> 28.472**. The treated arm's effective size falls by about two patients — the
wider model separates treated from control a little more sharply, so the overlap weights concentrate
— and that is a fact about the *weights*, which is what a sensitivity analysis on the propensity
specification is for.

### 5.3 It spends the negative controls, and this is what that bought

`config.NEGATIVE_CONTROLS` is computed as `PS_COVARIATES_FULL` minus `PS_COVARIATES`
(`config.py:571-572`), so the four covariates this arm adjusts for are exactly the four that stop
being controls in it. The roadmap requires that said where the arm's balance table is rendered;
§10's entry says it, and this is the number behind the sentence.

Measured (§21), |SMD| after weighting, [§7] then [§13]:

```
  hypertension     -0.139  ->  +0.026
  hyperlipidemia   -0.301  ->  -0.082
  diabetes         +0.380  ->  +0.052        the primary's WORST row
  smoking          +0.045  ->  -0.051

  center = HUG     +0.211  ->  +0.275        the arm's worst row — it gets WORSE
  center = Lugano  -0.154  ->  -0.199        also worse

  rows reaching SMD_THRESHOLD:  5  ->  2,  and both survivors are `center`
  worst |SMD|:  0.380 (diabetes)  ->  0.275 (center = HUG)
  penumbra_ml (in neither model):  +0.011  ->  +0.025
```

The first four numbers reproduce [§13]'s amendment of 2026-08-24 exactly — it quotes `diabetes`
0.380, `hyperlipidemia` −0.301, `hypertension` −0.139, `center = HUG` 0.211 and `center = Lugano`
−0.154 — which is a cross-check on this stage's arithmetic rather than a coincidence, and §21
records it as one.

**Two things follow, and the second is the one a reader will not expect.**

The arm does what DECISION 4 promoted it to do: the three risk factors that exceeded [§9]'s threshold
fall below it, the count of exceeding rows goes from five to two, and the worst residual |SMD| falls.
That is the axis that prompted the arm, measured on the axis.

And **centre gets worse on both rows**, which is [§13]'s own mechanism running in the direction
[§13] describes. Its amendment says overlap weights balance an in-model covariate exactly only under
maximum likelihood, that Firth solves the score equations *plus a penalty*, and that *"the deviation
is largest exactly where the design is closest to separated"*. A 16-column design on 92 records is
closer to separated than a 12-column one, so the penalty does more work and the deviation on the
near-determining covariate grows. This is not a defect in the arm and it is not grounds to prefer
either specification: it is the trade-off [§13] already made, visible from the other side. **What it
forbids is reading the arm's lower worst-|SMD| as uniformly better balance**, and §10's entry prints
both directions rather than the summary statistic alone.

`penumbra_ml` is the only balance-table covariate outside every model, in both specifications, and
after this arm it is the only negative control the analysis has left. [§13] accepts that: *"the
controls have already served their purpose."*
### 5.4 Its bootstrap is `replicates` and not `run`, and Stage 10 §14 said `run`

Stage 10 §12.1 and §14 item 3 both say the arm needs *"a second `run`"*. **It cannot be, and this
document records the departure rather than quietly taking a different route.** `run`'s `_replicate`
names `propensity.fit` at `bootstrap.py:578` — R9 (§4.5) now makes that a raise rather than a silent
mislabel — and `run` also drives `outcome.secondary` at `bootstrap.py:614`.

Two reasons, in this order, and the first is not about cost:

1. **[§13] gives this arm the primary only.** A `run` would produce seven binary outcomes' risk
   differences, marginal odds ratios and augmented estimates under the [§13] specification — nineteen
   more intervals that no prespecified section asks for. A table of unreported estimates is a
   multiplicity waiting to be read, and [§13]'s own Benjamini-Hochberg is scoped to seven p-values
   that are not these.
2. **Cost, measured.** Stage 10 measures its seven binaries at 123.0 s of a 130-145 s run. The arm's
   whole bootstrap, over `replicates` with a primary-only body, is **69.3 s** (§21).

So the arm owns a private body over the public `bootstrap.replicates`:

```python
def _arm_replicate(draw: pd.DataFrame, keys: tuple[str, ...], source: object) -> bootstrap.Replicate:
    """One [§13, DECISION 4] replicate: refit [§7] over PS_COVARIATES_FULL, then [§8], and nothing else.

    It is `bootstrap._replicate` with two things removed and one changed, and each is a decision:
    `outcome.secondary` is not called (§5.4), `_shared_design`'s invariant check is not run because
    no binary outcome is estimated here, and `propensity.fit` becomes `propensity.fit_full`.

    The `Audit` is constructed here and discarded, which is Stage 10 §6.2 unchanged: nine entries per
    replicate at N_BOOT would be memory that is a function of the replicate count, and no replicate's
    log is ever written.

    `C.SchemaError` is never caught. The only `except` is on `model.FitError`, which is Stage 10 §7.1
    as code and for its reason: a frame this stage constructed that cannot be read means this stage
    constructed it wrongly.

    The seven keys fail as ONE group. They come from one `polr` fit and a comparison across them is
    [§16]'s constant-shift statement, so they must share draws — Stage 10 §7.3, inherited rather than
    re-argued.
    """
    audit = Audit(source)
    try:
        ps = propensity.fit_full(draw, audit)               # the seam, INSIDE the loop
    except model.FitError as failure:
        return bootstrap.Replicate(
            values={}, failures={key: bootstrap.bucket(str(failure)) for key in keys},
            n_alpha=None, polr_iterations=None, sum_w=None, n_in_model=None)

    values: dict[str, float] = {}
    failures: dict[str, str] = {}
    n_alpha = polr_iterations = None
    try:
        est = outcome.primary(draw, ps, audit)              # UNCHANGED — §5.1
        values["beta"] = est.beta
        values.update({f"rd_{k}": v for k, v in est.rd.items() if f"rd_{k}" in keys})
        n_alpha, polr_iterations = len(est.fit.alpha), est.fit.iterations
    except model.FitError as failure:
        failures = {key: bootstrap.bucket(str(failure)) for key in keys}

    return bootstrap.Replicate(
        values=values, failures=failures, n_alpha=n_alpha, polr_iterations=polr_iterations,
        sum_w=float(ps.w[ps.in_model].sum()), n_in_model=int(ps.in_model.sum()))
```

**And the call site, because `replicates` takes a ONE-ARGUMENT body and `_arm_replicate` takes
three.** `bootstrap.replicates(df, body, n, seed, stratum)` declares
`body: Callable[[pd.DataFrame], object]` (`bootstrap.py:347-348`) and calls it as
`body(resample(df, rng, stratum))`, so the keys and the audit source are closed over and not passed.
This is `run`'s own wiring (`bootstrap.py:1096-1100`) and it is written out because a spec that
described the body but not its binding would leave the one line where the two runs' shapes are
visible to an implementer's guess:

```python
keys = ("beta", *(f"rd_{k}" for k in est_full.rd))            # §3.3 — off the ARM's Primary
collected = bootstrap.replicates(
    df,
    lambda draw: _arm_replicate(draw, keys, audit.source),    # keys and source CLOSED OVER
    C.N_BOOT, C.SEED, C.BOOT_STRATUM,                         # no second seed — §5.4
)
draws = bootstrap.collect(collected, keys)                    # the reconciling loop, not a second one
diagnostics = bootstrap.diagnostics(collected)                # ONE polr fit per replicate — §3.1
intervals = _intervals(draws, lambda key: False)              # no p on any arm key — §5.5
```

`keys` is built **before** `replicates` and not inside the body, which is Stage 10's stated reason
rather than a style: §7.1 requires a `propensity.fit` failure counted against *every* key, and a
body that named its own keys could not count a failure that happened before it named them.

**`bootstrap._bucket`, `bootstrap._collect` and `bootstrap._diagnostics` become public** — `bucket`,
`collect` and `diagnostics` — rather than being written again here. A second classifier would be a
second definition of the separation-versus-convergence split that Stage 8 §11 asked to be kept
apart; a second `Draws` loop would be a second implementation of the reconciliation property
`Draws.__post_init__` exists to enforce; and a second diagnostics tally would be a second answer to
"how many replicates fitted five cutpoints", which is the question [§16]'s constant-shift statement
turns on. This is the move `propensity.ess`, `balance.smd` and `outcome.weighted_proportion` each
made. Stage 10 §3.2's *"five names ... this count is stated once, here"* becomes **eight**, and
Stage 10 §1 and §13 cite that sentence (§22).

**No private alias is kept for any of the three.** `_bucket`, `_collect` and `_diagnostics` cease to
exist rather than becoming one-line forwarders, so `test_bootstrap.py`'s references move in the same
commit and a reader cannot find two spellings of one function. The rename is mechanical and the
suite is the proof it was complete: T5 lands red-to-green in one edit (§20).

`C.SEED`, `C.N_BOOT` and `C.BOOT_STRATUM` are reused unchanged and there is no second seed. So the
arm sees the **same sequence of resampled frames** as Stage 10's run: `replicates` creates one
`default_rng(seed)` and calls `resample` outside `body` with no `try`, so the stream is a function of
the seed alone. Measured (§21): two independent 200-replicate drives produce identical `case_id`
fingerprints, and both match a hand-driven `resample` loop over `default_rng(C.SEED)`.

**The pairing is asserted and never used.** Nothing computes a contrast between an arm draw and a
primary draw, so no false precision arises from the two runs not being independent; and if a later
analysis wanted such a contrast, the pairing is what would make it legitimate. §15.9 asserts it, so
it is a property rather than something a future reader discovers.

Measured (§21), and it is the number a reader of the arm's interval needs:

```
  1986 / 2000    surviving replicates on all seven keys
  13             degenerate_design
   1             nonconvergence
   0             separation — G7 never fires in the arm
  0.70%          drop rate, against Stage 10's measured 0.10% for the primary
  1898 / 91      replicates fitting 6 / 5 cutpoints — the same 4.6% Stage 10 measured
  88 .. 93       in_model per replicate; 93 is the replicates that drew the one
                 covariate-incomplete record zero times
  69.3 s         serial wall-clock
```

**The arm's drop rate is seven times the primary's, and that is the finding rather than the
footnote.** Both bodies refit the same two stages on the same 2000 frames; the only difference is
four covariates in the propensity design. Sixteen columns on 92 records with treatment nearly
determined by centre is closer to a degenerate design than twelve, so `propensity.fit`'s own guards —
F3 through F6 — fire more often. 0.70% is far above `ci_min_draws()`'s reach (98% of replicates would
have to fail before an interval is withheld), so every one of the seven intervals is emitted; and it
is the same fact §5.3 measures from the balance side, which is that the wider design is the more
strained one. §16 files it.

### 5.5 No p-value on any arm key

`bootstrap._tested("beta")` returns `True` (`bootstrap.py:456`), so an implementation that reuses
`_tested` to decide which arm keys get a p produces one for `beta`. It must not.

[§13] prescribes sensitivity analyses that *report* — "each reporting ESS and worst residual |SMD|" —
and prescribes no test of any of them. [§13]'s multiplicity clause is scoped to the secondary and
safety families; a p on the arm's `beta` would be an eighth test that [§13]'s correction is not told
about, and the primary's own p is [§10]'s and is uncorrected by prescription. `Arm.intervals[k].p is
None` for all seven keys, asserted (§15.3).

The six `rd_k` carry no p in the arm for the reason they carry none in the primary: [§8] gives the
cumulative risk differences intervals only, and Stage 10 §9.3 is the argument.

---

## 6. Benjamini-Hochberg within the two families [§13]

### 6.1 Seven p-values, and where they come from

```python
families = sec.by_family()                      # {"secondary": (...3), "safety": (...4)}
raw = {e.outcome: boot.intervals[f"{e.outcome}.rd"].p for e in families[fam]}
```

Three rules, each of which is a way to get the family sizes wrong:

- **The partition comes from `Secondary.by_family()` and never from `C.OUTCOMES[k].family`.**
  `outcome.py:757-761` states the reason as [§13]'s: a consumer that groups by reading the registry
  itself is a second place the partition is computed, and the correction is wrong if the two
  disagree. Measured (§21): `secondary` = `mrs_0_2_90d`, `mrs_0_1_90d`, `tici_2b_3`; `safety` =
  `sich`, `ph2`, `death_90d`, `mrs_5_6_90d`. **Three and four.**
- **The p-values come from `Interval.p` and are never recomputed.** Stage 10 §14 forbids it by name.
  `Interval` carries `method`, and a second computation under a different percentile definition is
  Stage 10 §8.2's disagreement re-introduced one stage downstream.
- **The family size is never obtained by counting `intervals`.** Stage 10 emits 26 keys of which 8
  carry a p — `beta` and the seven `<outcome>.rd` — so `len(boot.intervals)` is 26 and the count of
  p-carrying keys is 8, and neither is 3 or 4. Stage 10 §12.1 says *"Stage 11 must not discover a
  different count by counting intervals"*; §15.4 asserts the count comes from `by_family()` by
  constructing a `Bootstrap` with extra keys and requiring the answer not to move.

Three preconditions, each `C.SchemaError` and none repairable:

- **H3 — every family member's `.rd` key is present and its `p` is not `None`.** `bootstrap._tested`
  returns `True` for every `.rd`, so a `None` here is a Stage 10 contract break and not a sparse
  family. The absent-interval case is different and is §6.5's.
- **H4 — every raw p is in `[1/(n_draws + 1), 1.0]`, and `n_draws` is THAT KEY'S.** `bootstrap_p`
  floors at `1.0/(len(draws)+1)` over the draws it was given (`bootstrap.py:411`), and the draws it
  was given are that outcome's surviving ones — which is `Interval.n_draws` and is **not** `C.N_BOOT`.
  The distinction is not pedantic: the worst-off outcome carries about 1982 draws, so its floor is
  `1/1983` while `1/(C.N_BOOT + 1)` is `1/2001`, and a check written against `C.N_BOOT` is
  **strictly weaker than the floor it is checking** — it would accept a p below the smallest value
  `bootstrap_p` can return for that key, which is the one thing this precondition exists to catch. A
  p of exactly 0, or above 1, means something upstream produced it under a different definition.
- **H5 — `intervals["beta"]` is present and its `p` is not `None`.** `Multiplicity.primary_uncorrected`
  reads it (§3.1), so without this check an absent key is a bare `KeyError` from a reported field's
  initialiser and a `None` is an unlabelled `None` in the one number [§13]'s first sentence
  prescribes. The primary is **uncorrected, not unreported**: [§13] exempts it from the correction
  and requires it printed beside the corrected families, so a missing primary p is a failure of this
  deliverable and not an absence to disclose. Unreachable on v7 — the primary lost two replicates of
  2000 — and it is the same class of Stage 10 contract break H3 covers, from the one key that is not
  in a family.

### 6.2 The step-up, written out

```python
def benjamini_hochberg(p: np.ndarray) -> np.ndarray:
    """[§13]'s Benjamini-Hochberg adjusted p-values, over ONE family. Returns them in INPUT order.

    Public, and it is a definition rather than a step (§3.2): it names no outcome and no family, and
    it is checked against an oracle over arrays with no frame constructed (§15.4).

    Four properties, each of which is a test:

      * `kind="stable"` on the argsort, so tied p-values keep [§5] registry order. Ties are ORDINARY
        here and not theoretical: `bootstrap_p` at B = 1998 can return 1000 distinct values (§21), so
        a family of four is drawn from a small grid.
      * The DESCENDING RUNNING MIN is the monotonicity enforcement and there is no other. `m*p/i` is
        not monotone in `i` -- measured on sorted (0.01, 0.03, 0.04) the unenforced values are
        (0.03, 0.045, 0.04), so without it the middle outcome reports a LARGER adjusted p than the
        one below it and two of three secondary outcomes swap rank.
      * The UNSORT is `out[order] = q` and never `q`. An implementation that returns `q` returns the
        family sorted, which on a frame keyed by outcome silently reattaches every p to the wrong
        outcome.
      * The cap at 1.0 is INERT and is written anyway. q_(m) = min(1, m*p_(m)/m) = p_(m) for any
        valid p, and every other q_(i) is at most that, so the cap can bind only on an input §6.1
        already rejects. Measured over 10 000 random families of three and four drawn from the
        achievable grid: the pre-cap running min exceeds 1.0 in ZERO of them (§21). It stays for the
        reason `ci_min_draws`'s snap and `_shared_design`'s check stay: a specified-and-unreachable
        branch is one a future cohort does not discover the hard way.
    """
    p = np.asarray(p, dtype=float)
    m = p.size
    order = np.argsort(p, kind="stable")
    q = np.minimum(1.0, p[order] * m / np.arange(1, m + 1))
    q = np.minimum.accumulate(q[::-1])[::-1]
    out = np.empty_like(q)
    out[order] = q
    return out
```

### 6.3 One vector catches both defects, and it is the only one that does

Registry-order `p = (0.04, 0.01, 0.03)` at `m = 3`. Sorted it is `(0.01, 0.03, 0.04)`; unenforced
`(0.03, 0.045, 0.04)`; after the running min `(0.03, 0.04, 0.04)`; unsorted back
**`(0.04, 0.03, 0.04)`**.

- An implementation missing the running min returns `0.045` in the second position.
- An implementation missing the unsort returns `(0.03, 0.04, 0.04)`.
- An implementation missing both returns `(0.03, 0.045, 0.04)`.

Three distinct wrong answers, one input, and the assertion is to exact float equality on the products
rather than to a tolerance. The companion case `p = (0.01, 0.02, 0.03)` -> `(0.03, 0.03, 0.03)` is
also asserted and is **labelled in the test as the case a wrong implementation also gets right**, so
that a suite consisting only of it cannot look like coverage.

### 6.4 Implement it; `statsmodels` is the oracle and may not be imported

`pyproject.toml:26-31` declares `statsmodels` *"Retained for unpenalised cross-checks in tests, not
for any reported estimate"* and `scipy` *"TEST-ONLY, and deliberately so"*. **An adjusted p-value is
a reported estimate** — it goes in a manuscript table — so importing `multipletests` into a shipped
module to produce one is the thing that policy exists to forbid, for an algorithm with a four-line
closed form. It would also put `multipletests`' own sort/unsort and its `is_sorted` handling between
the bootstrap and the paper.

Its correct role is the one `pyproject.toml` already assigns it, and it is Stage 6 §16b's pattern
verbatim: **the oracle**. Measured (§21): the implementation above agrees with
`multipletests(..., method="fdr_bh")[1]` on all six hand-built cases, on the workbook's own two
families, and on **10 000 random families of size three and four drawn from the achievable p grid**,
at `rtol=1e-12, atol=0`. The same 10 000 confirm permutation invariance (10 000/10 000) and validity
— `adjusted >= raw` elementwise and `adjusted <= 1` (10 000/10 000).

The grid is `{1/(B+1)} u {2k/B}` and not `U(0, 1)`, because ties are the reachable condition and a
uniform draw never produces one. Stage 6 §12.7's AST scan is extended rather than duplicated:
**no shipped module imports `statsmodels` or `scipy`** (§15.11).

### 6.5 The family denominator when a p is absent

`Secondary.estimates` is complete at the point-estimate contract — `failures` is empty unless
`collect=True`, which is Stage 10's replicate call and not this one — so `by_family()` cannot come
back short *here*. **One route exists**: `run` emitted no `Interval` for `<outcome>.rd` because that
outcome's surviving draws fell below `C.ci_min_draws()` (`bootstrap.py:1099-1103`, Stage 10 §8.3).

**Decision, and [§13] does not fix it: `m` is the number of p-values that entered the procedure**, so
`m_used` and not `m_declared` is BH's scaling factor. Benjamini-Hochberg controls the false discovery
rate over tests *performed*, and a hypothesis with no p was not tested. The alternative — holding `m`
at the declared family size — is also valid and is conservative, and it is rejected rather than
overlooked: it inflates every adjusted p in the family to pay for a test nobody ran. §17 records it
PI-reversible; reversing it changes at most seven numbers and no code path.

The disclosure is typed and is not a footnote. `FamilyCorrection` carries `m_declared`, `m_used` and
`absent`, and §10's table prints the excluded outcome **as a row**, with its surviving-draw count and
the floor it fell below in the raw-p and adjusted-p cells — never as a blank, and never absent from
the table. §3.1 pins the reason string's shape so the formatter reads those two numbers rather than
parsing prose for them.

**If `m_used` is zero, `raw` and `adjusted` are both empty, `absent` holds every declared member, and
that is the whole of "says so".** There is no sentinel and no sentence field: a reader given
`m_declared = 4`, `m_used = 0` and four `absent` rows has been told everything, and the statement
Stage 14 prints is selected from `m_used == 0` (§12.2) rather than stored. `benjamini_hochberg` is
**not called at all** on that branch — it is never handed an empty array, so `m = p.size = 0` and the
`m*p/i` division by an empty range are both unreachable rather than defended against.

**Unreachable on v7.** Stage 10 measured the worst per-outcome drop rate at 0.8%, so the worst-off
outcome carries about 1982 draws against a floor of 40; an outcome needs 98% of its replicates to
fail to lose its interval. Specified anyway, and §15.5 reaches it with a constructed `Bootstrap`
rather than with data.

### 6.6 The safety family: [§10] and [§13] disagree, and this is what is done about it

[§10] closes with *"estimation, not testing, is the reportable output"* and the plan describes safety
outcomes resting on few events as descriptive. [§13] prescribes Benjamini-Hochberg **within the
safety family**. Computing an FDR correction over four p-values from a family the plan elsewhere
calls descriptive is prescribed by one section and deprecated by another, and neither section
mentions the other.

**Decision: compute it, because [§13] prescribes it, and carry [§10]'s label in the same row.**
`FamilyCorrection.descriptive_only` is a field, `True` for `safety`, so the adjusted p exists and
cannot be printed without the sentence that says what it is for. This is not inferable from either
section; §17 records it and `../out/questions_for_the_pi.md` carries it as a decision taken on the
documents rather than on the numbers, in the shape DECISION 4 was taken.

---

## 7. The E-value for the primary estimate [§6]

### 7.1 The formula, and the conversion

VanderWeele and Ding's bounding factor, on the risk-ratio scale, away from the null:

```python
def e_value(risk_ratio: float) -> float:
    """The [§6] E-value for a risk ratio. Public, and it is a definition (§3.2).

    ROOT FIRST, THEN INVERT, and the order is pinned rather than left to the implementer:
    `sqrt(1/OR)` and `1/sqrt(OR)` are mathematically equal and NOT bit-identical for every input.
    Measured on OR in {0.25, 0.3, 0.7, 0.9, 0.99} the two agree exactly (§21), so the pin is
    currently inert -- which is why it is a pin and not an assertion of exactness: §15.6 asserts the
    two forms agree to rtol=1e-15 and does not claim they agree to the bit.

    Takes a RISK RATIO and not an odds ratio, so the sqrt(OR) conversion happens at ONE call site
    (§7.5) rather than inside a function whose name does not mention it.
    """
    rr = float(risk_ratio)
    if rr < 1.0:
        rr = 1.0 / rr
    return rr + math.sqrt(rr * (rr - 1.0))
```

The primary effect measure is a **common odds ratio** from `model.polr`, so the roadmap prescribes
`RR ~ sqrt(OR)`:

```
  odds_ratio = est.odds_ratio          READ off the point estimate, never from the draws
  risk_ratio = sqrt(odds_ratio)
  e_point    = e_value(risk_ratio)
```

`est.odds_ratio` is `exp(est.beta)` (`outcome.py:148`) and the point estimate is not a function of
the draws — Stage 10 §8.1's converse, and §7.3's sign guard is the one place that matters.

### 7.2 Null-spanning is decided on `beta`, and the three routes agree for three different reasons

Stage 10 emits **no primary `odds_ratio` key**: `_estimand_keys` gives `beta` and the six `rd_k`, and
the `<outcome>.odds_ratio` keys are the binaries' (`bootstrap.py:431-454`). So the interval this
stage has is on the log-odds scale and the reported odds-ratio limits are `exp(lo), exp(hi)`.

**Spanning is evaluated on `boot.intervals["beta"]` against 0**, which is what Stage 10 §14 item 2
names. Three routes reach the same verdict and each is a different guarantee:

1. `lo <= 0 <= hi` **iff** `exp(lo) <= 1 <= exp(hi)`. `exp` is strictly increasing and
   `exp(0.0) == 1.0` exactly in IEEE-754 — measured (§21) — so this is an identity and not a
   tolerance.
2. **`exp` of the percentile limit IS the percentile limit of `exp` of the draws.**
   `C.PERCENTILE_METHOD = "inverted_cdf"` selects an *order statistic* — an actual data value at
   index `ceil(q*n)` (Stage 10 §8.2, §8.3) — and a strictly monotone transform commutes with
   order-statistic selection exactly. Measured (§21) at n in {41, 97, 400, 1998, 2000}: equal in all
   five under `inverted_cdf`, and equal in **none** of the five under `linear`, which interpolates.
   **So Stage 10 §8.2's pin is what licenses reporting an odds-ratio interval as `exp` of a `beta`
   interval at all** — a second consequence of that pin which §8.2 does not draw.
3. `boot.intervals["beta"].p >= 1.0 - C.CI_LEVEL` **iff** the interval includes 0, which is
   Stage 10 §9.4's asserted agreement. **The comparison is `>=` and not `>`, and the operator is
   pinned here rather than left to the implementer**, because the boundary case is the one this
   route exists to decide: `bootstrap_p` counts a draw of exactly 0.0 in both tails
   (`bootstrap.py:420-422`), so a p landing exactly on `1 - C.CI_LEVEL` is the arithmetic image of a
   limit landing exactly on 0 — which the tie rule below calls spanning. With `>` the two rules
   disagree at exactly that point and route 3 contradicts route 1.

§15.6 asserts all three give the same verdict on Stage 10's own tail-count-50 fixture rather than on
the workbook. **That fixture is the one where the three routes can be told apart, not one where
route 3 fails**: Stage 10 §9.4 measured agreement holding there under `inverted_cdf` and **failing
under `linear`**, so it is the frame on which a wrong percentile method — the only thing that breaks
routes 2 and 3 — is visible. Asserting agreement on the workbook's own draws would pass under
`linear` too (Stage 10 §9.4, 8 of 8) and would therefore assert nothing.

**Ties at the null count as spanning.** `lo == 0.0`, `hi == 0.0` or both are spanning cases. This is
consistent with `bootstrap_p`'s deliberate conservatism (`bootstrap.py:420-422`: draws exactly equal
to 0.0 are counted in both tails), and inconsistency here would break route 3.

### 7.3 The limit's E-value is 1.0, and the point estimate at the null is 1.0 without a branch

**Hard rule: `e_limit` is exactly `1.0` whenever the interval spans the null.** Not approximately,
not `nan`, not the formula applied to the limit.

The argument, rather than the assertion. The E-value answers *what is the minimum strength of
association, on the risk-ratio scale, that an unmeasured confounder would need with both treatment
and outcome to explain this result away*. If the interval already includes the null, the answer is
**none** — no confounding at all is needed for these data to be compatible with no effect — and the
minimum strength is `RR = 1`, which the formula returns as `1 + sqrt(1*0) = 1.0` exactly. Feeding a
null-crossing limit in instead answers a *different* question, "how far past the null does the
interval reach", and the mechanism that makes the answer misleading is that the formula is minimised
at 1 and increases in **both** directions once inverted: a limit on the far side of the null scores
as strong evidence. The roadmap says it returns *"a large number that reads as robustness while
meaning the opposite"*; this is why.

**The limit nearest the null is selected by sign and never by `min(|lo|, |hi|)`:**

```
  if lo > 0.0:   limit = lo
  elif hi < 0.0: limit = hi
  else:          spans the null -> e_limit = 1.0, limit = None
```

A `min(|.|)` implementation picks the wrong limit **only** in the spanning case, where the 1.0 rule
then masks it — so the defect is invisible until the first analysis whose interval excludes the null,
which is the analysis nobody wants to find it in. §15.6 asserts the selection on a deliberately
asymmetric non-spanning fixture where taking `hi` instead of `lo` gives a visibly different number.

**A point estimate exactly at the null gives exactly 1.0, with no branch, and a branch must not be
written.** `beta == 0.0` -> `odds_ratio == exp(0.0) == 1.0` -> `rr == 1.0` -> `rr*(rr-1) == 0.0` ->
`E == 1.0`, all exact. The pilot has an `if rd == 0: return 1.0` branch (`pilots/analysis.py:761`):
the right answer reached by the wrong route, which is worse than either, because it makes the reader
believe the formula needs help there.

**One guard the pilot does not have, and it raises.** Stage 10 §8.1 states plainly that a percentile
interval *"is not a function of the point estimate"*, so `beta < 0` with `lo > 0` is arithmetically
reachable — the point estimate outside its own interval. That is a finding to look at, not a number
to print: **`C.SchemaError`**, naming the estimate and both limits. The pilot returns `1.0`
(`pilots/analysis.py:761` tests `np.sign(rd_limit) != np.sign(rd)`), which reports it as
"the interval spans the null" — a different, benign condition. §15.6 asserts the raise and §21
records that it is unreachable on the workbook's own draws.

### 7.4 Every numerical edge case is unreachable, and G7 is why

- **Overflow is impossible, and the `sqrt(OR)` conversion is exactly what makes it so.** The naive
  form `rr + sqrt(rr*(rr-1))` overflows when `rr^2 > DBL_MAX`. Here `rr = sqrt(OR)`, so
  `rr^2 = OR <= DBL_MAX` always. No `sqrt(rr)*sqrt(rr-1)` rewrite is needed and none is written; the
  reachability argument goes in the docstring instead of a defensive form.
- **`OR` is bounded far below that, by a guard two stages back.** `outcome._assert_reportable`
  (`outcome.py:346`) raises `FitError` at `|beta| >= C.POLR_MAX_ABS_BETA = 14.0`, so
  `OR` is in `(8.3e-7, 1.2e6)`, `RR` in `(9.1e-4, 1097)` and `E <= 2192.77` — measured (§21).
  **G7 is what makes every numerical edge case in this deliverable unreachable**, which is worth
  stating because a reader will otherwise assume the E-value is unbounded.
- **`OR` extremely close to 1**: `E - 1 ~ sqrt(sqrt(OR) - 1)`, so `dE/dOR -> inf` as `OR -> 1`. E is
  **not smooth at the null**. Measured (§21): at `OR = 1 + 1e-12` — that is, `e_value(sqrt(OR))` — `E - 1 = 7.07e-07`, so a change of `1e-12`
  in the odds ratio moves E by `7e-07` — six orders of magnitude of amplification. Consequence for
  reporting: two decimals, and the statement that an E-value of 1.05 and one of 1.15 differ by about
  a hundredfold in the odds ratio behind them. §15.6 pins it as a property so it is not rediscovered
  as a surprise.
- **A non-finite limit or point estimate**: `C.SchemaError`, never `nan`. Unreachable by
  construction — `bootstrap.collect` raises on a non-finite draw and a percentile of finite draws is
  finite — so a non-finite input means a caller passed something Stage 10 cannot produce. The pilot
  returns `np.nan` (`pilots/analysis.py:740, 759`), which prints as `nan` in a manuscript table.
- **No `Interval` at all** (below `ci_min_draws()`): `e_point` is reported, `e_limit` and `limit` are
  `None`, `spans_null` is `False`, and nothing is substituted. **The draw count is recorded in
  `EValue.n_draws`** (§3.1) — read off `boot.draws["beta"]`, which Stage 10 keeps when it withholds
  the interval — so the `None` says *"too few surviving replicates, and here is how few"* rather than
  standing alone. Unreachable on v7 — the primary lost two replicates of 2000.

**There is no interval on the E-value itself.** Both E-values are deterministic functions of numbers
that already carry intervals, and Stage 10 §18 already declined a second interval on a transform of
an estimate for this reason.

### 7.5 What `sqrt(OR)` costs on an ordinal outcome, argued rather than glossed

The roadmap requires the approximation be *stated*. Stating it is not the same as pricing it, and
this is the section that prices it.

1. **The bound is derived for a risk ratio of a binary event, and the negative claim is now
   evidenced rather than asserted.** The bounding factor `RR_AU*RR_UY / (RR_AU + RR_UY - 1)` is
   derived for a single dichotomous outcome with a definable baseline risk. **There is no published
   bounding factor for a common odds ratio from a proportional-odds model**, and §19c is the
   literature check behind that sentence: the reference implementation — CRAN `EValue` 4.1.4, dated
   2026-05-07, authored by Mathur, Smith, Ding and VanderWeele — ships `evalues.RR`, `evalues.OR`,
   `evalues.HR`, `evalues.RD`, `evalues.MD`, `evalues.OLS` and `evalues.IC`, **and nothing for an
   ordinal or proportional-odds outcome.** The number this stage reports is a heuristic on an
   approximated scale, not the bound the derivation licenses, and `EValue.measure` carries
   `"common odds ratio [§8]"` so the object itself says which quantity was converted.
2. **`sqrt(OR)` is a binary-outcome correction with a DOCUMENTED validity range, applied to a
   quantity that has no single baseline risk.** The conversion is VanderWeele's, *"On a square-root
   transformation of the odds ratio for a common outcome"*, Epidemiology 2017;28(6):e58 — the
   citation `evalues.OR` itself carries — and the reference implementation encodes its boundary as an
   argument: `rare = 1` for an outcome under **15%** at end of follow-up, `rare = 0` for one over
   15%, with the square-root conversion applied only in the second case. So the approximation is not
   licensed by the outcome being "ordinal" or "common" in the abstract; it is licensed by a
   **prevalence**, and `exp(beta)` here is a cumulative-odds shift *assumed constant across all six
   thresholds*, each of which has its own baseline risk and therefore its own OR-to-RR factor.
   **This is checkable from a table the stage already computes**: `Primary.cumulative` holds
   `P(Y <= k)` per threshold per arm (`outcome.py:152`), so the control-arm value at each of
   `C.MRS_THRESHOLDS` can be read against the 15% boundary and the thresholds falling outside it
   named. §12.2 item 4 makes that naming a Stage 14 rule rather than a remark here, and no number
   from it appears in this document (§1).
3. **The direction of the error is statable and it is the conservative one.** `E` is strictly
   increasing in `RR` away from the null and `sqrt(OR) < OR` for `OR > 1`, so the prescribed
   conversion yields a **smaller** E-value than treating the odds ratio as a risk ratio. One sentence
   covers this. **No second E-value at `RR = OR` is reported as a bracket**: that is a second number
   for one estimand, and monotonicity already tells a reader which way it moves.
4. **The pilot's route is unavailable here, not merely different.** `pilots/analysis.py:729` converts
   a *risk difference* to a risk ratio at the *observed control risk*: `rr = (p_control + rd)/p_control`.
   That is a legitimate use of the same tool for a genuinely binary outcome. For the ordinal common
   odds ratio **neither ingredient exists** — there is no control risk of "mRS", and the contrast is
   not a difference of two proportions. So this is not the pilot's E-value computed differently; the
   pilot's cannot be computed at all.
5. **The constant-shift caveat is inside the E-value too.** Stage 10 §7.4 measured 4.6% of the
   primary's replicates fitting five cutpoints rather than six, and [§16] requires the constant-shift
   statement *computed, not written*. The E-value is a function of one number that assumes one shift
   across six thresholds; the assumption does not become safer for being transformed twice.

### 7.6 It is never printed alone

**The E-value bounds *unmeasured* confounding, and this analysis has measured confounding it has not
fixed.** [§9]'s threshold is exceeded after weighting on `center = HUG` and `center = Lugano` in the
primary specification, and §5.3 measures the arm making both worse. An E-value printed beside a
residual |SMD| of 0.275 invites the reading *"residual confounding would have to be strong"* when a
covariate that is **in the model** is visibly not balanced.

So the rule is structural rather than editorial: **the E-value is reported with `Balance.worst()`
beside it**, from the same specification, and `worst()` returns the undefined covariate names with it
(`balance.py:127-138`) so a worst |SMD| computed over rows some of which could not be judged cannot
be quoted alone. [§9]'s amendment of 2026-08-24 already requires the exceedance named at the point of
reading the estimate; the E-value is such a point, and §15.6 asserts the pairing by scanning the
rendered entry.

DECISION 4 draws the same line from the other direction and it is quoted here because it is the
sentence that keeps the two apart: *"It does not pre-empt the E-value [§6, §11], which addresses
unmeasured confounding. What the negative controls show is measured and visible; the two are
different questions and both are reported."*
---

## 8. The subgroups [§13]

### 8.1 The four candidates, and the estimand argument that decides it

[§13] names two subgroups — *"unknown versus witnessed onset; core volume above/below median"* — and
asks for *"subgroup estimates ... with an interaction test"*, hypothesis-generating. It names no
model, no scale, no estimator and no reported quantity. **All four are decided here and recorded as
decisions rather than inferred** (§17).

The candidates, ranked:

**(a) One weighted proportional-odds fit per subgroup, over the whole overlap-weighted population,
with a treatment × subgroup interaction. Chosen.** §8.2.

**(a′) Refit `outcome.primary` within each level using the same frozen weights and bootstrap
`β₁ − β₀`. Runner-up, and its advantage is real.** It is *saturated in S*, so each level's odds ratio
rests only on within-level proportional odds, where (a) additionally assumes the subgroup main effect
is itself a single proportional-odds shift — implausible for core volume, which moves the whole mRS
distribution rather than shifting it. Four things decide against it: the tested quantity under (a)
is **one coefficient with one null**, which is exactly what `bootstrap_p` tests and what
Stage 10 §9.1's *"one test of one parameter"* language is about; the reported level odds ratios and the test
come from one fit and therefore cannot disagree; two fits compound two independent failure surfaces
where one fit has one; and the pooled fit estimates its cutpoints on all 92 records, where two fits
at n ≈ 46 would each face Stage 10 §7.4's collapse more often. The tempting argument — that
collapsing categories changes `β` — is **not** used, because proportional odds is invariant to
collapsing adjacent categories and Stage 10 §16 item 5 says so.

**(b) Refit the propensity model within each level and bootstrap the difference. Rejected on an
estimand argument, not a power argument.** [§8] states that the ATO estimand *"is itself indexed by
the true propensity score"* — the target population is defined by `h(X) = e(X){1 − e(X)}` and
corresponds to no observable subset of the data — and contrasts it with the ATT and ATC, whose target
populations *"do not move when the propensity model changes"*. A within-level `e` defines a different
tilting function, so the two level effects are effects **in two different weighted populations** and
their difference is not an interaction. That is decisive on its own. Two further costs make it
unattractive even to someone who does not accept it: `model.design`'s constant-column rule would drop
the Lugano dummy in any level or replicate containing no Lugano patient, so the two levels would be
adjusted for different covariate sets and `Propensity.dropped` would record it only after the
estimate existed; and the degrees-of-freedom budget, which [§6] already records as tight at ~3.5
treated per parameter over 41 treated and 11 columns, falls to about 1.8 in a level — and **to 0.45
in the witnessed-onset level, which carries 5 treated patients** (§21).

**(c) Re-slice Stage 10's existing draws. Wrong, and Stage 10 §14 item 4 already says so.** Three
reasons and the first is fatal on its own. A Stage 10 replicate survives as **26 floats**: the drawn
frame, the weights and the fit are constructed inside `_replicate` and discarded (`bootstrap.py:578`,
Stage 10 §6.2), so there is no subgroup quantity in `Bootstrap` to slice and `beta`'s draws are draws
of the **full-cohort** common odds ratio. Second, even if each replicate had stored subgroup values,
`resample` fixes the **centre** totals and not the subgroup margin (`bootstrap.py:290`), so a
replicate's level size is random and a percentile interval over subgroup values would carry the
sampling variability of *level membership count* inside the effect's interval — and whether that
inflates or deflates depends on the level-size distribution, so it is not even conservative. Third,
`beta`'s surviving set is defined by the full-cohort primary fit converging, which is the wrong
denominator for a failure surface that is measured here at 12.5% (§8.4).

**The one number that governs this whole section**, measured on the point estimate (§21):

```
  unknown_onset            S = 0 (witnessed)      S = 1 (unwitnessed or wake-up)
    control  A = 0                16                        37
    treated  A = 1                 5                        34
      by centre, treated at S = 0:  HUG 4, CHUV 0, Lugano 1, USZ 0

  core_above_median        S = 0 (at or below)    S = 1 (above)
    control  A = 0                24                        29
    treated  A = 1                24                        15
```

**Five treated patients at witnessed onset, four of them at one centre.** [§13]'s amendment of
2026-08-10 justified keeping this subgroup on the cohort split — *"onset is unwitnessed or on waking
in 93 of the 126 records against 33 witnessed"* — which is the right number for the question it was
answering and is not the number that governs estimability. §8.4 is what that 5 turns into.

### 8.2 The model, and the one departure from [§8]

For each `S` in `C.SUBGROUPS`, over `in_estimate & df[S].notna()`, weighted by the **one prespecified**
`Propensity`'s `w`:

```
  logit P(Y <= k | A, S)  =  alpha_k  +  beta*A  +  delta*S  +  gamma*(A x S)

    exp(beta)          the common odds ratio at S = 0
    exp(beta + gamma)  the common odds ratio at S = 1
    exp(gamma)         the ratio of the two — the interaction, on the odds-ratio scale
    gamma              the tested coefficient; null gamma = 0
```

**No [§6] covariate enters, and this is not an omission.** [§8]: all estimates are marginal in the
overlap population; no covariate adjustment of the reported effect. Adding the confounders here would
make `exp(beta)` a **conditional** odds ratio — precisely what the roadmap forbids Stage 12 from
mislabelling: *"`exp(β)` here is a conditional odds ratio ... must never be labelled as the
standardised marginal effect."* Confounding control is the weights' job and it has already been done.

**`S` in the linear predictor is the one declared departure from [§8]'s weighting-only rule.** It has
to be there — without the main effect, `gamma` is not an interaction — and it is minimal: one column,
the subgroup variable itself, no confounder. `config.py:772-774` already asserts that no subgroup name
appears in any covariate list, so `S` cannot arrive as an adjustment variable by convenience, and
that assertion stays true because this design is built here rather than by widening a covariate list.

**Level coding is arithmetic and not `C.REFERENCE_LEVELS`.** Both subgroup columns are `Int64` 0/1
with `<NA>` (`derive.py:200`, `derive.py:260`); level 1 is the named condition. This must be written
down because the pilot takes its reference level from `pd.get_dummies(series.astype(str),
drop_first=True)` (`pilots/analysis.py:827`), which makes the **sign of `gamma` a function of string
sort order**.

**The frozen median is resampled and never recomputed.** `core_above_median` is cohort-dependent:
`derive.derive_cohort` computes it once after the [§3] restrictions and raises if asked twice
(`derive.py:388-405`), for the reason Stage 3 §6.4 gives — recomputing per replicate would give every
replicate its own cut-point and its own subgroup. `resample` copies whole rows, so the column travels
with the frame. §15.11 asserts by AST scan that `sensitivity.py` names neither `derive_cohort` nor
`_core_above_median`, which is the instrument Stage 10 §15.11 uses for `_augmented_path`.

**The subgroup mask is a no-op on v7 and is written anyway.** `unknown_onset` derives from
`onset_type` and `core_above_median` from `core_ml`, both [§6] covariates, so
`model.complete_cases(df, PS_COVARIATES)` already excludes any row missing either. Measured (§21):
`in_estimate & df[S].notna()` equals `in_estimate` exactly, at 92, for both subgroups. Without the
mask an `Int64` `<NA>` converts to `nan`, reaches `polr` and raises O1 (`model.py:792`) — **droppable
and counted**, so the failure would be absorbed as a sparse replicate rather than reported as a fact
about the frame. §15.7 asserts the equality positively.

### 8.3 The interaction column, and the two ways it can vanish

Built here, on a local copy, and passed through `model.design` as an ordinary covariate. The whole
fit, written out, because the two assertions and their **placement** are the specification:

```python
def _subgroup_fit(df: pd.DataFrame, ps: propensity.Propensity, in_estimate: pd.Series,
                  s: str) -> tuple[float, float, model.PolrFit, pd.Series]:
    """One subgroup's [§13] interaction fit. Returns (beta, gamma, fit, mask).

    Three columns and no [§6] covariate (§8.2). The product is built HERE, on a local copy, and
    passed to `model.design` as an ordinary covariate, so D1-D4 and the constant-column rule apply to
    it uniformly and `dropped` is reported. Building it after `design` would exempt it from the one
    rule that catches G8.

    Two assertions, and each is a failure if MOVED rather than merely if removed:

      * G8 sits BETWEEN `design` and `polr`, because a two-column design FITS -- `gamma` simply
        ceases to exist while every number returned is finite (§8.3 route one). This is
        `outcome._assert_exposure_survived`'s placement and its reason, on a different column.
      * G9 sits BETWEEN `polr` and everything that reads it, because a separated ordinal fit
        CONVERGES and its `exp(beta)` is a finite float every downstream table will accept
        (`model.py`'s `polr` docstring measures it at 36.4). This is `outcome._assert_reportable`'s
        placement and its reason.

    Both raise `model.FitError` and not `C.SchemaError`: a stratified resample that drew no treated
    patient at a level is a sparse replicate, which [§10] drops and counts. Neither is ever
    substituted with a different estimator [roadmap invariant 5]. **That reasoning is about the
    replicate loop, and this function is also called once on the POINT ESTIMATE, where there is
    nothing to drop into — so the caller converts a `FitError` from the point-estimate call into
    `C.SchemaError` H8 (§8.5). The granularity decision lives here; the point-estimate decision does
    not.**

    **G8 asserts the EXACT column tuple and not "nothing was lost", which is §15.7's assertion
    written once instead of twice.** `design` ends `X = ... .astype(float)` after a `get_dummies`
    over the declared factors only (`model.py`, `design`), and none of these three columns is in
    `C.CATEGORICAL` — so the returned columns are the passed tuple, in the passed order, and the
    equality is the strongest available check rather than a hopeful one. A weaker predicate over
    set membership would pass a future `design` that reordered columns while `beta` and `gamma` are
    read back BY NAME below, so the two would silently stop describing the same fit.
    """
    mask = in_estimate & df[s].notna()
    sub = df.loc[mask].copy()                      # local; the caller's frame is not touched
    ix = f"{C.TREATMENT}_x_{s}"
    sub[ix] = sub[C.TREATMENT].astype(float) * sub[s].astype(float)

    declared = (C.TREATMENT, s, ix)
    X, dropped = model.design(sub, declared)
    if tuple(X.columns) != declared or dropped:                  # G8 -- the EXACT tuple, §15.7's own
        raise model.FitError(
            f"G8  the treatment x subgroup design came back as {tuple(X.columns)} against a declared "
            f"{declared}, dropping {', '.join(dropped) or 'nothing'}. `design` removes a constant "
            "column SILENTLY, after which `polr` fits the remainder, converges, and returns a result "
            "in which the interaction does not exist -- so this is checked BEFORE the fit and not "
            "after it [Stage 11 §8.3].")

    fit = model.polr(X, sub[C.PRIMARY_OUTCOME].to_numpy(dtype=float),
                     ps.w.loc[mask].to_numpy(dtype=float))     # raises O1-O6; O6 is route two
    coefficients = dict(zip(X.columns, (float(v) for v in fit.beta)))
    beta, gamma = coefficients[C.TREATMENT], coefficients[ix]

    reported = max(abs(beta), abs(gamma), abs(beta + gamma))
    if reported >= C.POLR_MAX_ABS_BETA:
        raise model.FitError(
            f"G9  a reported quantity reached {reported:.3f} against a bound of "
            f"{C.POLR_MAX_ABS_BETA}. The fit CONVERGED; a separated proportional-odds fit does. The "
            "bound is on max(|beta|, |gamma|, |beta+gamma|) and NOT on the subgroup main effect, "
            "which is a nuisance that grows legitimately when a level's outcome distribution is "
            "concentrated [Stage 11 §8.4].")
    return beta, gamma, fit, mask
```

Through `design` rather than assembled by hand, so that D1-D4 (`model.py:129-201`) and the
constant-column rule apply to the product uniformly, and `dropped` is reported. Building the product
after `design` would exempt it from the one rule that catches the first failure below.

**`model.py:120-126`'s design-matrix banner says "No interaction, no spline, no transform — [§6] says
'linear terms only'".** That sentence is about the [§6] confounder set and it reads as a global
prohibition; it gains one clause naming Stage 11 and stating that the product column is built by the
caller and passed in as an ordinary covariate, so `design` itself stays innocent of interactions
(§13). [§6]'s "linear terms only, no interactions" and [§13]'s interaction test look contradictory and
are not: the first is about confounders and the second about the exposure × subgroup contrast.

**Route one: `A × S` is constant and `design` drops it.** It is all-zero exactly when no treated
patient has `S = 1`. Then `design` removes it silently, `polr` fits a two-column model, converges,
and returns — `gamma` does not exist while every number in the output is finite and plausible. This
is `outcome._assert_exposure_survived`'s failure mode (`outcome.py:320`) on a different column and it
gets the same treatment: an assertion **between `design` and `polr`**, requiring
`tuple(X.columns) == (C.TREATMENT, s, ix)` and `dropped == ()`, raising **`model.FitError` and not
`C.SchemaError`** — because a stratified resample that drew no treated patient in a level is a sparse
replicate and must be droppable, which is `outcome.py:320`'s own stated reason.

**Route two, and it is the one this cohort actually takes.** When the `(A = 1, S = 0)` cell is empty,
`A × S` is **not** constant — it equals `A` — so `design` keeps both columns and the design is
**rank-deficient**, not degenerate. `model._assert_polr_fittable`'s **O6** catches it
(`model.py:837-843`), raises `FitError`, and `C.FAILURE_BUCKETS["O6"]` already buckets it as
`degenerate_design`. Nothing new is needed for the route that fires.

**Measured, and this is the finding that changed §8.3 (§21):**

```
                          route one (constant)   route two (O6 rank)   propensity failure
  unknown_onset                    0                     9                    2
  core_above_median                0                     0                    2

  the 9 are EXACTLY the 9 replicates that drew zero treated patients at witnessed onset
  the 2 are EXACTLY Stage 10's measured two propensity.fit FitErrors, in both subgroups
```

The first draft of this document specified route one as the reachable one and route two not at all.
Route one fires **zero** times in 2000 replicates and is specified anyway — with its own
`C.FAILURE_BUCKETS` token, `G8`, bucketed `degenerate_design` — for the reason `ci_min_draws`'s snap
and Stage 10 §5's preconditions are kept: it is unreachable on **this** cohort and the assertion is
what would notice on the next one. The two-propensity-failure agreement is a cross-check on the
pairing of §5.4 rather than a coincidence, and §21 records it as one.

### 8.4 The bound is on the reported quantities, and it is not inert

`model.polr`'s own docstring records the fact this section exists for: **a separated proportional-odds
fit converges and returns** — measured there at `beta` 36.4, `exp(beta)` 6.5e15, 17 iterations, every
safeguard clean. `outcome._assert_reportable` (`outcome.py:346`) is the guard on the primary, and it
bounds `max|fit.beta|` over **all** coefficients at `C.POLR_MAX_ABS_BETA = 14.0`.

**Decision: bound `max(|beta|, |gamma|, |beta + gamma|)` — the three reported quantities — and not
`|delta|`.** This is G7's own stated logic, that the bound is on the reported effect and not on a
nuisance which can legitimately be large: `delta` is the subgroup main effect and grows when a
level's outcome distribution is concentrated, and dropping a replicate on it would lose the
interaction estimate for a reason unrelated to the interaction. It raises `model.FitError` with a new
token, `G9`, bucketed `separation` — the bucket Stage 8 §11 asked to be kept separate from
`nonconvergence` precisely so that a counter reading zero says something.

**Measured (§21), and the draft's assumption that the guard would be inert on v7 is wrong:**

```
                                   unknown_onset      core_above_median
  replicates fitted                     1989                1998
  max |beta|                          18.368               5.326
  max |gamma|                         19.296               4.166
  max |beta + gamma|                   2.383               3.519
  max |delta|                          5.257               5.520
  |reported| median / p90 / p99      2.317 / 16.187 / 18.176    0.964 / 2.024 / 3.096
  replicates with |reported| >= 14       239                   0
  replicates with |delta|    >= 14         0                   0
  polr iterations in those 239        10 .. 15   (against 3-4 typical)
  surviving after the bound             1750                1998
  drop rate against N_BOOT             12.5%                0.1%
```

Four things follow, and the fourth is the one that has to be said out loud.

- **The witnessed-onset level is separated or near-separated in about one replicate in eight.**
  `|beta|` and `|gamma|` blow up together while `|beta + gamma|` stays under 2.4 — the signature of
  the `S = 0` level being unidentified while the `S = 1` level is fine, which is what five treated
  patients buys. Without the bound, the interval on `exp(beta)` would be a quantile over a mixture
  that includes odds ratios of order 1e8.
- **The bound choice is inert on v7 and is still the right rule.** `|delta|` never reaches 14, so
  bounding everything and bounding the reported quantities give **the same 1750 survivors** — the
  count of replicates where only `delta` exceeds is zero. The decision is recorded as principled
  rather than as load-bearing, which is the honest description of a rule that costs nothing here.
- **A cell-count rule would not substitute for it.** The 239 over-bound replicates carry between 1
  and 9 treated patients at `S = 0`, and only 86 of them have 2 or fewer; a threshold on the cell
  would drop replicates that fit and keep replicates that did not. The bound is on the fit.
- **The atomic group means a separated `S = 0` level costs the `S = 1` estimate too.** Stage 10 §7.3
  is inherited unchanged: the three keys come from one fit and are dropped together. That is right —
  a likelihood flat in one direction is not a fit whose other coefficients are trustworthy — and it
  is expensive here, because `beta + gamma` was well determined in every one of those 239.

**And the drop is SELECTION ON THE ESTIMATE, so the surviving interval is narrower and not merely
noisier. This is the sentence the drop rate needs and it is not the same sentence as "250 replicates
were lost."** G9's condition is `max(|beta|, |gamma|, |beta + gamma|) >= C.POLR_MAX_ABS_BETA` — a
condition on the quantity being estimated, not on the frame. So the 1750 survivors are a sample
selected on the estimate; `percentile_ci` over them takes quantiles of a **truncated** sampling
distribution; and truncation at `|.| >= 14` removes mass from both tails of `beta` and `gamma` and
none from the middle, so the limits are **systematically tighter** than the untruncated ones. The
direction is knowable even though the magnitude is not, and a reader told only the count will read
the interval as unbiased-but-thinned, which it is not.

Three things follow and none of them is that the bound is wrong.

- **The bound is the least-bad option and not a good one.** Without it the interval on `exp(beta)` is
  a quantile over a mixture containing `exp(18.368)`, which is above 9e7 (§21). That is not a
  conservative interval, it is an uninterpretable one. Stage 8 §6.2 makes the same argument for G7.
- **The truncation bites hardest on the level that was not the problem.** `max|beta + gamma|` over all
  1989 fitted replicates is **2.383**, so the `S = 1` level's odds ratio was estimable in every one
  of the 239 dropped replicates. Stage 10 §7.3's atomic-group rule drops it anyway, correctly — a
  likelihood flat in one direction is not a fit whose other coefficients are trustworthy — and this
  is the precise cost of that correctness.
- **Nothing here is repaired in code**, and §17 does not list a reversible decision for it, because
  the three options are not implementation choices: report with the mechanism stated, report the
  `S = 1` level only, or report no interval at all. The second would require relaxing the atomic-group
  rule for one subgroup, which changes how failures are counted rather than what is printed.
  **The implementation does the first**, and `../out/questions_for_the_pi.md` Question 4 item 1a puts
  all three, because choosing among them is a reporting decision on a hypothesis-generating quantity
  and belongs to the PI (§16 item 1).

**A 12.5% drop rate is above anything this pipeline has produced and it is not silently accepted.**
Stage 10's worst was 0.8%. `../out/questions_for_the_pi.md` carries a standing question deferred with
*"ask after Stage 10"* — **whether a drop rate above some level invalidates the interval** — and
records that it was premature because no drop rate existed. **This stage produces the first one that
makes the question live**, and §16 and §17 put it with the number rather than answering it here. What
is not deferred: the surviving count is printed beside every subgroup interval, which is
Stage 10 §12.3's rule, and the witnessed-onset level's treated count is printed beside its odds ratio.

### 8.5 The replicate body, one run for both subgroups

```python
def _subgroup_replicate(draw: pd.DataFrame, keys: tuple[str, ...],
                        source: object) -> bootstrap.Replicate:
    """One replicate: refit [§7] over the WHOLE population, then one interaction fit per subgroup.

    **The propensity model is refit in every replicate and over the whole `in_model` population of
    the draw, never within a level.** [§10] requires the refit and §8.1(b) is why the population is
    not the level's. This is the single line that separates the chosen design from the rejected one.

    Failure granularity mirrors Stage 10 §7.1 and §7.3 exactly: a `propensity.fit` failure fails the
    whole replicate and is counted against all six keys; a subgroup's own fit failure costs that
    subgroup's THREE keys atomically and costs the other subgroup nothing. Measured: the propensity
    route fires twice, in both subgroups, and they are Stage 10's same two replicates (§8.3).

    One propensity fit serves both subgroups. It is 13.4 ms of Stage 10's measured 79.6 ms per
    replicate, so a second run for the second subgroup would double the expensive part for nothing.

    **`n_alpha` and `polr_iterations` are `None` and that is the decision, not a shortcut.** They are
    scalars on `Replicate` (`bootstrap.py:148-149`) because [§10]'s body and the arm's each run one
    `polr` fit; this body runs two, so any single value would silently describe one subgroup while
    being tallied as the replicate's. `Subgroups` therefore carries no `bootstrap.Diagnostics` (§3.1)
    and every counter it reports comes off `Draws`, which is where §8.4's 12.5% already lives.
    `sum_w` and `n_in_model` ARE populated, because there is exactly one propensity fit and they are
    unambiguous; `subgroup_replicates` prints the `n_in_model` range over the collected replicates
    (§10) and that is their one consumer.

    **`propensity.fit` and never `fit_full`.** The subgroups are estimated under the ONE prespecified
    specification (§8.2); the [§13] arm is a different deliverable and combining them would produce a
    subgroup estimate under a sensitivity specification that no section asks for.
    """
    audit = Audit(source)
    try:
        ps = propensity.fit(draw, audit)                   # the [§7] specification — NOT fit_full
    except model.FitError as failure:
        return bootstrap.Replicate(                        # the whole replicate — §7.1, inherited
            values={}, failures={key: bootstrap.bucket(str(failure)) for key in keys},
            n_alpha=None, polr_iterations=None, sum_w=None, n_in_model=None)

    in_estimate = outcome.estimation_population(draw, ps)  # ONE definition of the [§11] mask — §8.5
    values: dict[str, float] = {}
    failures: dict[str, str] = {}
    for s in C.SUBGROUPS:                                  # declared order; no set iteration
        group = (f"{s}.beta", f"{s}.gamma", f"{s}.beta_plus_gamma")   # THREE, atomically — §7.3
        try:
            beta, gamma, _, _ = _subgroup_fit(draw, ps, in_estimate, s)
        except model.FitError as failure:                  # G8, G9, O1-O6 — this subgroup only
            failures.update({key: bootstrap.bucket(str(failure)) for key in group})
            continue
        values.update({group[0]: beta, group[1]: gamma, group[2]: beta + gamma})

    return bootstrap.Replicate(
        values=values, failures=failures,
        n_alpha=None, polr_iterations=None,                # TWO fits, one scalar field — see above
        sum_w=float(ps.w[ps.in_model].sum()), n_in_model=int(ps.in_model.sum()))
```

**`outcome.estimation_population` is extracted for this call and it is the one addition to
`outcome.py` that is not an audit step name.** `outcome.primary` computes the [§11] population as
`ps.in_model & df[C.PRIMARY_OUTCOME].notna()` and binds it once (`outcome.py:634`); a replicate needs
the same mask on the drawn frame and cannot have the `Primary` that carries it without paying for a
primary fit it does not read. The mask becomes a public one-line definition, `primary` calls it, and
this body calls it — which is `outcome.weighted_proportion`'s own argument applied to the quantity
that argument names: *"a second implementation in a later stage would be a second definition of the
denominator, and [§11] requires the denominator."* The alternative — spelling the mask again here
with a comment pointing at `outcome.py` — is the shape §6.1 rejects for the family partition, two
computations of one prespecified thing with nothing asserting they agree. **No number moves**, and
§15.12 asserts that by running the pipeline twice rather than by reading the diff.

**The point-estimate call is not in a loop, and a `FitError` there is `C.SchemaError` H8.**
`_subgroup_fit` runs once per subgroup on the cohort before any replicate (§8.2), and the
droppable-and-counted argument that makes G8, G9 and O6 `model.FitError` is an argument about the
replicate loop — there is nothing for the point estimate to be dropped into. So `subgroups` wraps the
two point-estimate calls and re-raises as `C.SchemaError` H8, naming the subgroup and the leading
token. A subgroup whose point estimate cannot be fitted is a fact about the cohort and about [§13]'s
choice of subgroups, and the two wrong answers are both silent: omitting that subgroup from
`estimates` is `pilots/analysis.py:811-813`'s uncounted `continue` (§19), and letting the `FitError`
out contradicts §11's contract and hands the caller a token whose bucket map is about replicates.
Unreachable on v7 — the point-estimate fits reach `|γ|` 1.841 and 0.562 against a bound of 14.0
(§21) — and reached by `separable_subgroup_frame()` instead (§15.8).

Six keys, one `bootstrap.replicates` call, `C.SEED` / `C.N_BOOT` / `C.BOOT_STRATUM` unchanged, and
the same wiring §5.4 writes out for the arm:

```python
keys = tuple(f"{s}.{q}" for s in C.SUBGROUPS                  # §3.3 — six, off the registry
             for q in ("beta", "gamma", "beta_plus_gamma"))
collected = bootstrap.replicates(
    df, lambda draw: _subgroup_replicate(draw, keys, audit.source),
    C.N_BOOT, C.SEED, C.BOOT_STRATUM)
draws = bootstrap.collect(collected, keys)
intervals = _intervals(draws, _tested_subgroup)               # p on `.gamma` ONLY — §8.5
```

No `bootstrap.diagnostics` call here, for §3.1's reason. The subgroup run therefore draws the **same
sequence of frames** as Stage 10's run and as the arm's, and §15.9 asserts it. Measured (§21):
**68.3 s**.

**A p-value on `.gamma` and on nothing else.** The rule is a *function* and not a frozen set —
`_tested_subgroup(key) -> key.endswith(".gamma")` — which is `bootstrap._tested`'s own move
(`bootstrap.py:456`) for its stated reason: a third subgroup restored to `C.SUBGROUPS` gains its p in
the same edit that gains its estimate. The two level odds ratios get intervals and `p is None`,
because [§13] asks for **an** interaction test and two level-specific p-values would be two tests of
one contrast — which is [§8]'s own objection to six threshold-wise tests, one stage on.

`exp` is applied to the interval **limits** and the coefficients are never exponentiated before the
percentile: §7.2's order-statistic equivariance is what makes those two the same numbers, and doing
it on the limits keeps one array of draws per estimand rather than two.

No new interval kind, no new percentile method, no new floor. `C.ci_min_draws()` governs and an
estimand below it keeps its `Draws` and gets no `Interval` — Stage 10 §8.3's rule, inherited.
Measured: all six keys are far above the floor, the worst at 1750.

### 8.6 Not per level: the six cumulative `RD_k`, and what is reported instead

**The six cumulative risk differences are not reported per level, with intervals or without.** [§13]
does not ask for them; twelve more intervals from a hypothesis-generating analysis on levels of ~21
and ~71 records is twelve more numbers a reader will read as findings; and [§8]'s own argument
against threshold-wise reporting — that it invites selection of the most favourable cut — applies
with more force per level, not less.

**What is reported instead is the pair of weighted cumulative distributions per level, point values
only, no interval and no p**, from the landed `outcome.weighted_proportion` and
`outcome.cumulative_rd` (`outcome.py:159`, `outcome.py:180`). It adds no estimator, and it is the one
diagnostic that adjudicates (a)'s extra assumption over (a′): it is [§16]'s constant-shift statement
evaluated per level, so a reader can see whether a single proportional-odds shift for `delta` does
violence to the data. `SubgroupEstimate.cumulative` is that table (§3.1).

### 8.7 Non-collapsibility, and the assertion nobody may write

A weighted proportional-odds model is **not collapsible**. The primary `exp(beta)` is therefore *not*
a weighted average of `exp(beta)` and `exp(beta + gamma)` from a subgroup fit, and the two level odds
ratios **need not bracket it**.

Stage 10 §8.4 already measured this from the other side: its design C has a true marginal ATO odds
ratio of 0.659 against a conditional generating coefficient of 0.7, and *"a coverage test that used
0.7 as design C's truth would report about 88% coverage and read as a bootstrap defect."*

Two consequences, and the second is a prohibition:

- The caveat is a printed field on the record and travels with the numbers.
- **No code anywhere asserts the bracketing.** A reviewer's instinct is to add
  `assert min(or_level) <= primary_or <= max(or_level)` as a sanity check, and it would fail on
  correct output. §15.7 asserts the prohibition by scanning the module and the test file for the
  comparison.

### 8.8 `gamma`'s p does not enter Benjamini-Hochberg

Argued from [§13]'s own structure and recorded as a decision (§17):

1. [§13] puts Benjamini-Hochberg in its **Multiplicity** clause and scopes it *"within the secondary
   and safety families"* — those are `C.OUTCOMES` families. The subgroups are a separate clause and
   they are on the **primary** outcome, whose first prescription is *"Primary outcome uncorrected."*
   Correcting an interaction test on the primary outcome inside an outcome-family correction
   contradicts both halves.
2. Stage 10 §12.1 pins the family sizes at three and four and warns that Stage 11 *"must not discover
   a different count"*. Adding one or two interaction p-values changes those denominators for a
   reason that has nothing to do with the outcome families.
3. [§13] labels the subgroups hypothesis-generating. An FDR correction controls a discovery rate
   among hypotheses one intends to act on; applying it here dresses the analysis up as confirmatory.

**The disclosure that replaces it is mandatory, and it is the honest half.** The raw `gamma` p is
reported with the hypothesis-generating label, with **the number of interaction tests performed — two
— printed beside it** (`SubgroupEstimate.n_interaction_tests`, §3.1), and with the statement that
those two are themselves an uncorrected multiplicity that [§13] does not address. That last sentence
is where [§13] is silent and where this stage must not be.

---

## 9. The failure taxonomy this stage adds, and it is two tokens

Stage 10 §7.2 classifies a `model.FitError` by the leading token of its message, against
`C.FAILURE_BUCKETS`, and `bootstrap.bucket` **raises** on an unrecognised token rather than
defaulting — which is what makes the token map a specification and not a lookup. Stage 10 shipped one
defect of exactly this kind: `propensity.py` raised three tokens its §7.2 scan scope did not cover,
and F5 was *"the only `FitError` the workbook itself produces"* (`config.py:327-337`).

So this stage adds its two tokens to `C.FAILURE_BUCKETS` **and extends the scan scope in the same
edit**:

```
  G8  ->  degenerate_design    the interaction column was dropped as constant (§8.3, route one)
  G9  ->  separation           a reported quantity reached POLR_MAX_ABS_BETA (§8.4)
```

`G8` is `degenerate_design` because the design lost a column, which is what that bucket already means
for O1-O6 and F3-F6. `G9` is `separation` because it is G7's guard on a different set of
coefficients, and putting it in `nonconvergence` would merge the two counters Stage 8 §11 required
kept apart — *"it measured that the second never fires on this estimator and a single counter reading
zero therefore says nothing."*

**Everything else this stage can raise is already bucketed and nothing is re-bucketed.** O6 covers the
rank-deficiency route that actually fires (§8.3); the propensity failures come through
`propensity.fit`'s own tokens; the arm's fourteen drops are `degenerate_design` and `nonconvergence`
and are Stage 8's and Stage 6's guards unchanged.

`test_bootstrap.py` §15.6 scans `model.py` and `outcome.py` for leading tokens and asserts each is a
`C.FAILURE_BUCKETS` key. **Its scope grows to `sensitivity.py` in the same commit that adds G8 and
G9** (§15.8). The alternative — putting the two raises in `outcome.py` so the existing scan covers
them — is declined, because Stage 8 §12 makes `outcome.py` the only shipped module that may name the
[§5] primary outcome and a subgroup estimator in it would be that rule spent for a scan's
convenience.

### 9.1 The ten `C.SchemaError` identifiers, and they are the H series

`FitError` tokens are classified; `SchemaError` identifiers are matched. §15.1 requires *"each
`C.SchemaError` message matched by its identifier"*, which is unsatisfiable unless the identifiers
exist and are declared somewhere a test can quote. Every prior stage declares its own block — D for
`design`, F for `propensity.fit`, G for `outcome`, O for `polr`, R for `bootstrap.run` — and this
stage's `FitError` tokens continue G at G8 and G9, so its `SchemaError` identifiers take the next
free letter, **H**, and a reader can tell the droppable from the fatal by the letter alone.

```
  H1   _assert_arm_inputs   `ps.spec is not C.PROPENSITY_PRIMARY` — the REFERENCE half of the
                            pair (§4.5). The arm's own half is checked by Arm.__post_init__ (H10)
  H2   _assert_arm_inputs   `ps.in_model.index` is not `df.index` — the alignment class R3 checks
                            in `run`, on the one function here that takes a frame and a Propensity
  H3   multiplicity         a family member's `<outcome>.rd` key is absent, or its `p` is None (§6.1)
  H4   multiplicity         a raw p is outside [1/(n_draws + 1), 1.0], per KEY (§6.1)
  H5   multiplicity         `intervals["beta"]` is absent, or its `p` is None (§6.1)
  H6   e_value_primary      a non-finite odds ratio or interval limit (§7.4)
  H7   e_value_primary      the point estimate is outside its own interval (§7.3)
  H8   subgroups            a POINT-ESTIMATE `_subgroup_fit` raised G8, G9 or O1-O6 (§8.5)
  H9   subgroups            `est.in_estimate` is not a subset of `ps.in_model`, or the indices
                            disagree — the caller error that would make the fits and the weights
                            describe different rows
  H10  Arm.__post_init__    `ps.spec is not C.PROPENSITY_FULL` or `reference.spec is not
                            C.PROPENSITY_PRIMARY` — by identity, never equality (§3.1, §4.5)
```

**No H identifier is in `C.FAILURE_BUCKETS` and none ever reaches `bootstrap.bucket`.** `bucket`
raises on an unrecognised leading token rather than defaulting (§9), so an H token arriving there
would present as a taxonomy gap; it cannot, because `C.SchemaError` is never caught anywhere in this
module (§5.4) and the only `except` in either replicate body is on `model.FitError`. §15.8's scan
asserts the separation from both sides: every leading token in a `sensitivity.py` `FitError` is a
`C.FAILURE_BUCKETS` key, and no H identifier is.

---

## 10. The audit entries, and there are fifteen

No new `data.KINDS` entry. All fifteen are `model` — claims about what was fitted to a population —
which is the kind `data.py:143-148` declares for exactly that, and `KINDS` has been nine for **four**
stages: Stage 10 §10.1 recorded three and this is the fourth. (§21 is normative on this and an
earlier draft of this section said five.)

```
  full_covariate  ->  10
    covariate_completeness_full_covariate    §4.3 — propensity.fit_full's four, suffixed
    design_matrix_full_covariate                    so `Audit.entry` cannot return the
    propensity_fit_full_covariate                   primary's (§4.4)
    overlap_weights_full_covariate
    balance_smd_full_covariate               §4.4 — assess's two; the role column moves
    overlap_by_centre_full_covariate                and the table does not
    outcome_completeness_full_covariate      §5.1 — outcome.primary's three, through
    primary_fit_full_covariate                      ps.spec.step; the ESTIMATOR is unchanged
    cumulative_rd_full_covariate
    sensitivity_arm                          §5.4 — the counters and the pair the
                                                    Accept-when is about

  subgroups       ->   3
    subgroup_fit_unknown_onset               §8.2 — one per C.SUBGROUPS key, and the cells
    subgroup_fit_core_above_median                  §8.1 turns on are its table
    subgroup_replicates                      §8.5 — counters, buckets, the 12.5%

  multiplicity    ->   1
    multiplicity                             §6 — both families, raw and adjusted, m_declared
                                                  and m_used, and the primary's exemption

  e_value         ->   1
    e_value_primary                          §7 — the approximation with its citation, both
                                                  E-values, §7.6's worst residual |SMD|, and
                                                  the per-threshold baseline risk against the
                                                  15% band (§7.5, §12.2 item 4a)
```

**31 -> 46.** The base is Stage 10's measured 31 and the fifteen are counted from the `audit.record`
calls this document prescribes; §15.13 asserts the total and each step name, because a ledger count
derived by arithmetic rather than by running is exactly the kind Stage 10 §10.1 required be asserted.

Six things about them, each a decision:

- **`sensitivity_arm`'s `n` is the arm's `in_model` and its table carries BOTH specifications' rows
  side by side.** The Accept-when is about a pair (§5.2) and an entry showing one half makes the
  comparison the reader's arithmetic. Its counters come from `Arm.draws` and `Arm.diagnostics`
  (§3.1) — the surviving count and buckets per key from the first, the cutpoint and iteration
  tallies and the `in_model` range from the second — so every number in the entry is a field of the
  returned record and none is recomputed at the print site.
- **`subgroup_replicates` reports off `Draws` and there is no `Diagnostics` behind it** (§3.1): per
  key, `n_attempted`, the surviving count and the bucket counts, which is what makes §8.4's 12.5%
  and the 1750 readable in the log. It also prints the `n_in_model` range over the collected
  replicates, tallied in `subgroups` from the one propensity fit each replicate made (§8.5).
- **`case_ids` is empty on all fifteen.** Nine `model` entries name cases and these remove no
  patient; the arm's own excluded record is named by `_record_exclusion` under the suffixed step,
  which is `propensity.py:400-432`'s rule applying to the arm unchanged.
- **`subgroup_fit_*` prints the four cells and the treated count per level**, so the 5 that governs
  §8.4 is visible to a reader of the log without this document.
- **`multiplicity` prints the absent-outcome row rather than omitting it** (§6.5), and
  `e_value_primary` prints the worst residual |SMD| in the same table (§7.6). Both are the rules of
  their sections rendered rather than restated.
- **`e_value_primary`'s table also carries one row per `C.MRS_THRESHOLDS` entry: the control-arm
  `P(Y <= k)` and whether it clears the 15% prevalence boundary the `sqrt(OR)` conversion is licensed
  for** (§7.5, §19c). It is read off `est.cumulative`, which the function already receives, so no
  signature moves and nothing is recomputed. **This is the only place those numbers exist**, because
  §1 keeps them out of this document and §12.2 item 4a requires Stage 14 to name the thresholds
  outside the band — a requirement Stage 14 cannot meet from a log that does not print them. A
  reader of the log can therefore see at which cut-points the scale the E-value was converted on is
  trustworthy, without this spec and without recomputing a cumulative distribution.

---

## 11. The four entry points, written out

```python
def multiplicity(sec: outcome.Secondary, boot: bootstrap.Bootstrap, audit: Audit) -> Multiplicity:
    """[§13]'s Benjamini-Hochberg within the secondary and safety families. Primary uncorrected.

    Takes no family list and no method name: the partition is `sec.by_family()` and the procedure is
    [§13]'s, so a keyword for either would make a prespecified choice look like an option — Stage 6
    §6.4's rule, and this stage's §4.1 is what it costs to break it.

    Reads `boot.intervals[f"{k}.rd"].p` and recomputes nothing (§6.1). Raises `C.SchemaError` H3 on
    a `None` p, H4 on a p outside `[1/(n_draws + 1), 1]` where `n_draws` is THAT KEY'S, and H5 on an
    absent or untested `intervals["beta"]` — and on nothing else: an absent `<outcome>.rd` INTERVAL
    is §6.5's disclosure and not an error, which is the one distinction this function's preconditions
    exist to draw.

    Appends one `model` entry (§10).
    """


def e_value_primary(est: outcome.Primary, bal: balance.Balance,
                    boot: bootstrap.Bootstrap, audit: Audit) -> EValue:
    """[§6]'s E-value for the primary estimate and for the limit nearest the null.

    **`bal` is a parameter and it is the whole of §7.6.** The E-value bounds UNMEASURED confounding
    and this analysis has measured confounding it has not fixed; the worst residual |SMD| is printed
    in the same table, from the same specification, and a signature that did not require it would
    make the pairing a convention a reporting layer could drop.

    Takes the `Primary` for `odds_ratio` and the `Bootstrap` for `intervals["beta"]`, and nothing
    from the draws: the point estimate is not a function of them (Stage 10 §8.1) and neither is the
    limit a function of the point estimate — §7.3's sign guard is where that stops being an academic
    remark.

    Raises `C.SchemaError` H6 on a non-finite input and H7 on an estimate outside its own interval
    (§7.3). Returns `1.0` for `e_limit` when the interval spans the null; never `nan`, ever. Fills
    `n_draws` on both branches, so an absent interval reports how few draws there were rather than
    only that there were too few (§3.1, §7.4).

    Appends one `model` entry (§10).
    """


def subgroups(df: pd.DataFrame, ps: propensity.Propensity, est: outcome.Primary,
              audit: Audit) -> Subgroups:
    """[§13]'s two subgroup estimates on the primary outcome, with one interaction test each.

    Takes no subgroup list: `C.SUBGROUPS` is the [§13] registry after the amendment of 2026-08-10
    withdrew the third, and a third restored there gains its estimate, its three keys, its interval
    and its p in one edit (§3.3).

    Takes `ps` and refits it per replicate over the WHOLE population — never within a level (§8.1).
    Takes `est` for `in_estimate`, which is the [§11] population the point-estimate fits run on, so
    the subgroup denominators are the primary's and not a fourth mask; each replicate recovers the
    same mask on its own drawn frame through `outcome.estimation_population`, which is the one
    definition both paths read (§8.5).

    Raises `C.SchemaError` H9 on the caller error `bootstrap.run`'s R-series covers, and H8 when a
    POINT-ESTIMATE `_subgroup_fit` raises — that call is outside the replicate loop and has nothing
    to be dropped into (§8.5). Inside the loop `model.FitError` never leaves it: G8, G9, O1-O6 and
    the propensity tokens are all droppable and are counted (§9).

    Returns no `bootstrap.Diagnostics` and that is §3.1's decision, not an omission: two `polr` fits
    per replicate cannot fill one scalar `n_alpha`.

    Appends three `model` entries (§10). Hypothesis-generating, and the label is a field (§3.1).
    """


def full_covariate(df: pd.DataFrame, ps: propensity.Propensity, audit: Audit) -> Arm:
    """The [§13, DECISION 4] full-covariate propensity sensitivity analysis.

    Takes no covariate list, for the reason `propensity.fit` takes none: the specification is
    `C.PROPENSITY_FULL` and it is named, not passed (§4). `ps` is the PRIMARY's fit and is taken as
    the reference half of the pair the roadmap's Accept-when is about (§3.1, §5.2) — it is read for
    its `in_model`, its `ess` and its `spec`, and it is never refitted here.

    Four things about the order below, each of which is a failure if moved:

      * `_assert_arm_inputs` first — H1 and H2 — so a caller passing the ARM's Propensity as `ps`
        reports it rather than producing a record comparing the arm with itself (§4.5).
      * `balance.assess` BEFORE the bootstrap, so the log carries the 19 re-roled rows even if the
        loop then raises — `propensity.py:440-443`'s rule, two stages on.
      * `differs_only_in_specification` is computed from the two `in_model` masks AFTER both exist
        and is never a literal (§5.2).
      * The loop is `bootstrap.replicates`' and not this function's, for Stage 10 §11's reason: a
        second loop would put the seeding and the ordering in two places. Its wiring — the closed-over
        keys, `collect`, `diagnostics` and `_intervals` — is written out in §5.4 rather than left to
        the implementer, because `replicates` takes a one-argument body and this body takes three.

    Appends ten `model` entries (§10).
    """
```
---

## 12. Data flow into Stages 12-14

### 12.1 Stages 12 and 13 [§14a, §14b]

**They read nothing from this stage and they inherit one thing from it.** [§14] fits no propensity
model anywhere, so `Specification`, `fit_full` and the whole of §4 are inapplicable there: the
standardisation covariate set is `C.STANDARDISATION_COVARIATES` (`config.py:562-563`), which is a
covariate *list* passed to `model.polr`, and that is not a seam — `polr` was made general in its
covariates at Stage 8 §0.1 precisely so [§14a] could widen the design.

What they inherit is **`bucket`, `collect` and `diagnostics` becoming public** (§5.4). Stage 10 §12.2
already gives them `resample`, `replicates` and `percentile_ci`; the three additions mean a Stage 12
replicate body classifies its failures against the same map, builds the same reconciling `Draws` and
tallies the same cutpoint distribution, rather than writing a third taxonomy. **`_intervals` is
`sensitivity.py`-private and is NOT offered to them** (§3.2): it is two `bootstrap` calls under one
floor rule, and a fourth caller wanting it is the trigger to move it into `bootstrap.py` beside the
loop it duplicates — §16 files that. Stage 10 §12.2's open question — whether a stratum at a centre with no
treated patients is resampled or held fixed — is theirs and this stage does not touch it.

### 12.2 Stage 14 [§16]

Stage 14 receives four records and owes eight guardrails this stage cannot enforce, because each is a
property of what is *printed*:

1. **The primary is printed uncorrected and beside the corrected families**, with the statement that
   it is uncorrected. `Multiplicity.primary_uncorrected` is a field so that this is a rendering of
   data rather than a sentence (§3.1).
2. **Raw and adjusted p are both printed**, [§13]'s explicit requirement, with `m_declared` and
   `m_used` beside them so the scaling factor is checkable by hand (§6.5).
3. **The safety family's adjusted p never appears without [§10]'s descriptive-only label**
   (§6.6). `FamilyCorrection.descriptive_only` is the field that carries it.
4. **The E-value never appears without the `RR ~ sqrt(OR)` statement and without the worst residual
   |SMD| beside it** (§7.5, §7.6). `EValue.approximation` and `Balance.worst()` are the two fields
   that make it structural, and the statement carries VanderWeele 2017's citation because the
   conversion is his and not this plan's (§19c).
4a. **And it names the thresholds whose baseline risk falls outside the conversion's documented
   range.** The square-root conversion is licensed for an outcome above **15%** prevalence (§7.5,
   §19c); `Primary.cumulative` holds the control-arm `P(Y <= k)` at each of `C.MRS_THRESHOLDS`, so
   Stage 14 reads it and names the thresholds outside the band. One E-value is reported, on the
   common odds ratio, and this is the sentence that says at which cut-points the scale it was
   converted on is trustworthy. Without it a reader takes the approximation as uniformly valid across
   an ordinal scale whose extreme thresholds it is not valid for.
5. **Every subgroup number carries the hypothesis-generating label, the count of interaction tests,
   and its surviving-draw count** (§8.8, §8.4). The 12.5% drop rate on `unknown_onset` is not a
   footnote: it is the denominator of that subgroup's interval.
6. **The arm and the primary are printed together**, with both `in_model` counts, both ESS, both
   worst residual |SMD|, and the `differs_only_in_specification` clause **selected from the field**
   and never written as prose (§5.2). The roadmap's Accept-when is this row.
7. **The arm's balance table names the spent negative controls where it is rendered** (§5.3), and it
   prints centre getting worse alongside the risk factors getting better — the summary statistic
   alone reads as uniformly better balance and it is not (§5.3).

**Two things Stage 14 must not do**, each because a prespecified rule says so: correct the primary
for multiplicity [§13]; and print the two subgroup level odds ratios as bracketing the primary's,
which they need not do (§8.7).

---

## 13. What Stage 11 amends in Stages 1-10

Seven shipped modules. The config additions are paste-ready and their comments are **shipping text
and not commentary**.

```python
# --- the propensity specifications [§7, §13] --------------------------------------------------
#
# TWO, and there is a registry, because [§13] lists SIX propensity-and-population sensitivity rows
# and DECISION 4 (PI, 2026-08-24) promoted exactly ONE of them from deferred to prespecified. A
# third arrives as a third Specification WITH a [§13] amendment behind it; test_config.py asserts
# the registry has these two and that PROPENSITY_FULL.covariates is PROPENSITY_PRIMARY's plus
# NEGATIVE_CONTROLS, so "it spends the negative controls, by construction and not by accident"
# [roadmap Stage 11] is a static check rather than a sentence.
#
# `suffix` is what keeps two fits over ONE cohort apart in ONE audit log. Audit.entry is first-match
# (data.py:273-275), so two fits recording `propensity_fit` would make every programmatic read
# return the primary's while the rendered log looked complete — measured, and it is why the payload
# is a record and not a tuple [Stage 11 §4.2, §4.4]. The PRIMARY'S SUFFIX IS THE EMPTY STRING, so
# every step name Stages 6 and 7 already record is byte-identical and the four landed ledger
# assertions do not move. test_config.py asserts exactly one declared specification has an empty
# suffix and that all suffixes are distinct.
#
# `covariates` is the REQUESTED list and never the post-design survivors: in_model is
# complete_cases(df, covariates), so the requested list is what decided the [§11] denominator, and a
# covariate dropped as constant is still adjusted for [Stage 11 §4.2].
@dataclass(frozen=True)
class Specification:
    label: str
    covariates: tuple[str, ...]
    suffix: str
    sap: str

    def step(self, base: str) -> str:
        return base + self.suffix


PROPENSITY_PRIMARY: Final[Specification] = Specification(
    "[§7] prespecified", PS_COVARIATES, "", "[§7]")
PROPENSITY_FULL: Final[Specification] = Specification(
    "[§13] full-covariate sensitivity", PS_COVARIATES_FULL, "_full_covariate",
    "[§13], DECISION 4")
PROPENSITY_SPECIFICATIONS: Final[tuple[Specification, ...]] = (
    PROPENSITY_PRIMARY, PROPENSITY_FULL)
```

```python
# Two tokens added for Stage 11, and the scan scope grows with them [Stage 11 §9].
#
# G8 is the treatment x subgroup interaction column dropped as CONSTANT — no treated patient at the
# named level — which model.design removes silently, after which polr fits a two-column model,
# converges, and returns a result in which gamma does not exist. Measured 0 of 2000 on this cohort
# and specified anyway, because the route that DOES fire is a different one: an empty (A=1, S=0)
# cell makes A and A x S identical rather than constant, which is rank deficiency and which O6
# already catches and already buckets. Both are degenerate_design: the design lost a column, one
# visibly and one to collinearity.
#
# G9 is a REPORTED quantity reaching POLR_MAX_ABS_BETA — max(|beta|, |gamma|, |beta+gamma|), and
# deliberately NOT |delta|, which is a nuisance that grows legitimately when a level's outcome
# distribution is concentrated. It is `separation` and not `nonconvergence` because it is G7's guard
# on a different set of coefficients, and Stage 8 §11 requires the two counters kept apart. Measured
# 239 of 1989 fitted replicates on `unknown_onset` and 0 on `core_above_median`: this bound is NOT
# inert, and five treated patients at witnessed onset is why [Stage 11 §8.4].
FAILURE_BUCKETS: Final[dict[str, str]] = {
    # the fifteen landed entries are unchanged and are not repeated here
    "G8": "degenerate_design",
    "G9": "separation",
}
```

| Module | What changes | Why |
|---|---|---|
| `propensity.py` | `_fit(df, spec, audit)` private; `fit` and `fit_full` public and named; `Propensity.spec`, appended, **no default**; eight privates take the record; four `audit.record` calls go through `spec.step` | §4.3. `fit`'s docstring at `:459-463` carries a promissory note this design pays and must be rewritten in the past tense pointing at `fit_full`. The stale `data.py:255` citation at `:413` becomes `data.py:271` — its `TODOS.md` trigger is met (§16) |
| `balance.py` | `_role(name, spec)`, `_table(..., spec, ...)`, `Balance.spec`; the three role strings become `Final`; `assess`'s signature **unchanged** | §4.4, and Stage 7 §11 promised it in these words. `_smd_detail:539-560` claims four negative controls, which is false in the arm; `:562-564` says *"[§7] prescribes one propensity model and there is no second one to try"*, which stops being true at this stage |
| `outcome.py` | Three `audit.record` step names go through `ps.spec.step`; the `:634` mask becomes the public `estimation_population(df, ps)` and `primary` calls it. **The estimator is untouched** | §10, §8.5. `primary` is called twice over one cohort and its three steps collide exactly as `propensity.fit`'s do (§4.4). The extraction gives the [§11] denominator ONE definition, which is `weighted_proportion`'s own argument. §15.12 asserts the estimator is unchanged **by number**, not by diff |
| `bootstrap.py` | R9 in `_assert_run_inputs`; `_bucket -> bucket`, `_collect -> collect`, `_diagnostics -> diagnostics`, no aliases kept; §3.2's *"five names"* becomes **eight** and its **twelve** privates become **nine** | §4.5, §5.4. §0.2's *"`propensity.py` is not amended — §5's central claim"* is a claim about Stage 10 and becomes false read as present tense; it is rescoped, and Stage 10's DoD item 2 is marked historical (§4.3). Stage 10 states these two counts in its §1, its §3.2, its acceptance criteria and its §21, and every one of them moves — which is §16 item 6's trigger firing on the first change to this surface |
| `model.py` | One clause on the design-matrix banner at `:120-126` | §8.3. *"No interaction, no spline, no transform"* is about the [§6] confounder set and reads as a global prohibition |
| `config.py` | The two fences above | §4.2, §9 |
| `data.py` | **Nothing.** `KINDS` stays nine, for a fourth stage | §10. Fifteen `model` entries, no new kind |

**Spec documents that gain amendment notes** (§22): `stage6_...md` §6.4 (the bullet is discharged),
§7.2 (step names are now suffixed), §11's `fit` fence, §14's *"which propensity specifications [§13]
eventually runs"* — answered for one; `stage7_...md` §4.1 and §4.3 (future tense becomes present),
§11 (fulfilled), §14's *"what the [§13] full-covariate specification's balance table will say.
Deferred"* — answered; `stage9_...md` §3.2 (`outcome.py`'s public count gains
`estimation_population`); `stage10_...md` §3.2 (five names becomes **eight**, twelve privates become
nine), §5 and DoD item 2 (rescoped), §12.1 and §14 item 3 (both say *"a second `run`"*, and §5.4 is
why it is not).

---

## 14. Handover to Stage 12

```
  cohort = cohort.build(eligibility.classify(derive.derive(*data.load()), ...), ...)
  ps     = propensity.fit(cohort, audit)                 spec = PROPENSITY_PRIMARY
  bal    = balance.assess(cohort, ps, audit)
  est    = outcome.primary(cohort, ps, audit)
  sec    = outcome.secondary(cohort, ps, audit)
  boot   = bootstrap.run(cohort, ps, est, sec, audit)    R9: ps.spec must be the primary's
  mult   = sensitivity.multiplicity(sec, boot, audit)            H3-H5 on the p's it reads
  ev     = sensitivity.e_value_primary(est, bal, boot, audit)
  subs   = sensitivity.subgroups(cohort, ps, est, audit)
  arm    = sensitivity.full_covariate(cohort, ps, audit)

    92 / 92        in_model, [§7] and [§13] — the Accept-when is TRUE on v7        [§5.2]
    5 -> 2         balance rows reaching SMD_THRESHOLD, and both survivors
                   are `center`, which gets WORSE                                   [§5.3]
    1986 / 2000    the arm's surviving replicates, 69.3 s                           [§5.4]
    1750 / 2000    `unknown_onset`'s, after G9 — a 12.5% drop rate, and the
                   first this pipeline has produced that makes the standing
                   PI question about drop rates live                                [§8.4]
    1998 / 2000    `core_above_median`'s, 68.3 s for both subgroups together        [§8.5]
    9 / 0          O6 rank-deficiency drops / constant-interaction drops, and
                   the second is the route the first draft specified                [§8.3]
    3 and 4        the [§13] family sizes, off `Secondary.by_family()`              [§6.1]
    15             audit entries, taking the ledger 31 -> 46                        [§10]

    every adjusted p, every E-value, every subgroup odds ratio and every arm
    limit is DELIBERATELY ABSENT from this ledger and from this document.
    They are in the gitignored log.                                                 [§1]
```

**Three things Stage 12 must know, and none of them is a number.**

1. **`Propensity` now carries a specification and `bootstrap.run` checks it.** [§14] fits no
   propensity model, so Stage 12 constructs no `Propensity` and R9 never fires for it — but Stage 12
   *does* build replicate bodies, and `bucket`, `collect` and `diagnostics` are now the public
   functions to build them with (§12.1). `_intervals` is not offered and §12.1 says why.
2. **`replicates` has now been driven by three different bodies** — [§10]'s, the arm's and the
   subgroups' — and the property Stage 10 §15.2 asserted about it holds across all three: same seed,
   same frame, same stratum gives the byte-identical sequence of drawn frames regardless of what the
   body does or raises. Measured (§21). Stage 12's standardisation body inherits that.
3. **A 12.5% drop rate exists now.** Stage 12's [§14a] population is larger and includes a centre
   with no treated patients, so its own degenerate cases will differ; what transfers is that this
   pipeline's answer to a high drop rate is to **count it, print it and ask**, never to substitute an
   estimator (roadmap invariant 5).

---

## 15. Acceptance criteria

`tests/test_sensitivity.py`, with fixtures in `tests/fixtures_stage11.py`. Banner comments are the
`### 15.x` headings below, lower-cased, verbatim. Sections whose subject lives in another module go
in that module's test file and are tagged.

### 15.0 The frames and fixtures this stage is tested on

Five fixtures, and each exists because a branch of this document is unreachable on the workbook.

- **`separable_subgroup_frame()`** — a cohort-shaped frame in which the `(A = 1, S = 0)` cell is
  **empty**, so O6 fires (§8.3 route two); and a companion in which every patient has `S = 1`, so
  `A × S` is constant and G8 fires (§8.3 route one, measured 0 of 2000 on the workbook).
- **`near_separated_subgroup_frame()`** — one treated patient at `S = 0`, chosen so the fit
  **converges** and returns `|gamma|` above `POLR_MAX_ABS_BETA`. *This is the fixture the whole of
  §8.4 rests on*: the specification is that a fit which succeeds is rejected, and a frame on which
  the fit failed instead would pass a wrong implementation.
- **`family_fixtures()`** — `Secondary` and `Bootstrap` pairs: one complete; one with an outcome
  whose `.rd` interval was not emitted (§6.5); one whose `intervals` carries extra keys, so that an
  implementation counting intervals gets a different `m` and fails (§6.1); one with a `None` p, which
  must raise H3; one carrying two `Interval`s with **different `n_draws`** and a raw p at exactly the
  smaller one's floor, which passes while the same p at the larger `n_draws` raises H4 (§15.4); and
  one with **no `"beta"` key at all**, which must raise H5. The last two are hand-built `Interval`s:
  neither is a thing `bootstrap_p` or `run` can produce, and the tests say so.
- **`interval_fixtures()`** — `Interval`s on `beta`: spanning with `lo < 0 < hi`; `lo == 0.0`;
  `hi == 0.0`; both zero; a wide one; a non-spanning asymmetric one where taking `hi` instead of `lo`
  gives a visibly different E-value; and one whose limits exclude the point estimate (§7.3's raise).
- **Stage 10's `tail_count_50` fixture, reused and not rebuilt** — the one fixture where route 3 of
  §7.2 can disagree with routes 1 and 2.

### 15.1 The preconditions, one frame per branch

**Ten H identifiers and R9, and the count comes from §9.1 rather than from this list.** `H1`, `H2`
(`_assert_arm_inputs`); `H3`, `H4`, `H5` (`multiplicity`); `H6`, `H7` (`e_value_primary`); `H8`, `H9`
(`subgroups`); `H10` (`Arm.__post_init__`); and R9 in `test_bootstrap.py` (§15.2). Each is asserted on
a frame that reaches it and nothing else, and each `C.SchemaError` message is matched by its
identifier — `pytest.raises(C.SchemaError, match="H4")` and not on prose, so a reworded message does
not break the suite and a renumbered guard does.

**Two of them need a hand-built input and are the ones a suite is most likely to skip.** H4 needs an
`Interval` whose `p` is below its own `1/(n_draws + 1)` floor, which no `bootstrap_p` call can
produce — so the fixture constructs the `Interval` directly, and the test is labelled as asserting a
check against a value the pipeline cannot generate. H5 needs a `Bootstrap` whose `intervals` has no
`"beta"` key at all, which `run` cannot emit above the floor; `family_fixtures()` carries it.

### 15.2 R9, and the mutation companion `[in test_bootstrap.py]`

`bootstrap.run(df, ps_full, ...)` raises `C.SchemaError` matching `R9`. **And the companion**: with
R9 removed the call **succeeds** and returns 26 intervals — asserted by monkeypatching the check out,
because a guard whose absence is never demonstrated is a guard nobody can price (§4.5).

### 15.3 The arm reuses every landed stage and adds no estimator

`sensitivity.py` calls `propensity.fit` in exactly one body and `fit_full` in exactly one (AST scan,
the pattern at `test_balance.py:1425`) — the subgroups' and the arm's respectively, which is the one
scan that catches the two bodies being wired to each other's specification; the arm's
`Propensity.spec is C.PROPENSITY_FULL`; `Arm.__post_init__` raises H10 on either half of the pair
being wrong; all seven `Arm.intervals[k].p is None` (§5.5) **and the companion that `_intervals` is
what enforces it** — the same helper with `_tested_subgroup` produces two p's on the same draws, so
the `lambda key: False` is asserted to be load-bearing rather than decorative; `Arm.diagnostics`
reconciles, with `n_alpha` summing to the replicates that reached `outcome.primary` and
`max_abs_beta` and `or_corrected` both empty (§3.1); and the behavioural falsifier — the arm's
`draws["beta"]` differs element-wise from Stage 10's, **with the companion that the two runs' k-th
drawn frames are identical**, so the difference is attributable to the specification and to nothing
else.

### 15.4 Benjamini-Hochberg, and the vector that catches both defects

§6.3's `(0.04, 0.01, 0.03) -> (0.04, 0.03, 0.04)` to exact float equality; the companion case
labelled as one a wrong implementation also passes; `m = 1` bit-identical; §6.2's cap case;
permutation invariance over 1000 permutations; tie equality bit-identical; validity and the
statsmodels oracle over 10 000 random families from the achievable grid at `rtol=1e-12, atol=0`; and
the family sizes taken from `by_family()` on a `Bootstrap` carrying extra keys.

**And H4's denominator, which is a test of the check and not of the procedure.** Two `Interval`s on
one family with different `n_draws` — 1998 and 1982 — and a raw p at exactly `1/1983`: it passes,
and the same p with `n_draws` 1998 raises H4. A check written against `C.N_BOOT` passes both and
therefore passes a p below the smallest value `bootstrap_p` could have returned for that key (§6.1).

### 15.5 The absent p, and the family denominator `[unreachable on v7]`

On `family_fixtures()`' thinned pair: `m_used < m_declared`, the adjusted p computed with `m_used`,
the absent outcome present in `absent` with its reason, and the rendered table carrying it **as a
row** with its draw count and the floor.

### 15.6 The E-value

`E(1.0) == 1.0` and `E(sqrt(4)) == 2 + sqrt(2)` bit-exact; `E(4) == E(0.25)` — measured exactly equal
and asserted to `rtol=1e-15` because exactness is not claimed; strict monotonicity on a grid to
`RR = 1e6`; the near-null sensitivity, `e_value(sqrt(1 + 1e-12)) - 1 > 1e-7`; the 1.0 rule on all five spanning
fixtures to exact equality; nearest-limit selection on the asymmetric fixture; the three routes
agreeing on `tail_count_50` **and route 3 asserted at the boundary** — a p exactly equal to
`1.0 - C.CI_LEVEL` is spanning, which is the case `>` gets wrong and `>=` gets right (§7.2);
`exp(percentile_ci(x)) == percentile_ci(exp(x))` bit-exact under `C.PERCENTILE_METHOD` **and not**
under `"linear"`; `C.SchemaError` H6 on non-finite input and H7 on the estimate-outside-its-interval
frame; `n_draws` populated on both branches, including the no-interval one where it is read off
`boot.draws["beta"]` (§7.4); and a scan asserting `approximation` is read from the field.

### 15.7 The subgroup model, and what it may not contain

`tuple(X.columns) == (C.TREATMENT, S, ix)` — three columns, no [§6] covariate — asserted by
inspecting the design the fit ran on, and it is **the same predicate G8 raises on** (§8.3) rather
than a second spelling of it; `outcome.estimation_population` returns the mask `outcome.primary`
binds at `:634`, asserted by comparing it against `est.in_estimate` on the cohort so the extraction
is proved equivalent and not merely intended to be (§8.5); the subgroup mask equals `in_estimate` on
the workbook `[data-gated]`; `exp(beta)` and `exp(beta + gamma)` come from one fit; the level coding is
arithmetic, asserted by constructing a frame whose subgroup labels sort the other way and requiring
`gamma`'s sign not to move; **and the prohibition** — a module and test-file scan for a bracketing
assertion between a subgroup odds ratio and the primary's (§8.7).

### 15.8 The two degeneracies, and the bucket map is scanned not trusted

G8 on the constant-interaction fixture, G9 on the near-separated one, O6 on the empty-cell one; each
raises `model.FitError` and not `C.SchemaError`; each is dropped and counted rather than
substituted **[roadmap invariant 5]**; each classifies into the bucket §9 declares. **And the
atomicity**: on a frame where one subgroup fails and the other does not, that subgroup's three keys
carry a bucket and the other subgroup's three carry draws, in one `Replicate`, which is §8.5's
granularity rule as an assertion rather than a docstring.

**The point-estimate half, which is the same fixtures through a different door.** The same three
frames driven through `subgroups` rather than through a replicate raise `C.SchemaError` H8 naming the
subgroup and the token, and the message is asserted to carry the G8/G9/O6 token so a reader of the
failure knows which degeneracy it was (§8.5).

And the scan, from both sides: every leading token raised in `sensitivity.py` is a
`C.FAILURE_BUCKETS` key — which is `test_bootstrap.py` §15.6's scan with its scope grown (§9) — and
no `H` identifier is one, because H is `C.SchemaError`'s block and `bucket` raises on anything it
does not know (§9.1).

### 15.9 The three runs draw the same frames, and the two dicts are never merged

Same seed, same frame, same stratum: Stage 10's, the arm's and the subgroups' drawn `case_id`
sequences are identical over 200 replicates, and identical to a hand-driven `resample` loop.
**And the prohibition**: the test constructs `{**boot.draws, **arm.draws}`, asserts it has 26 keys
rather than 33, and asserts that seven of them now describe the [§13] specification — so the cost of
the merge is a measured fact in the suite rather than a warning in a docstring (§3.3).

### 15.10 The registry has exactly two, and a third cannot appear un-spec'd `[in test_config.py]`

`PROPENSITY_SPECIFICATIONS == (PROPENSITY_PRIMARY, PROPENSITY_FULL)`;
`PROPENSITY_FULL.covariates == PROPENSITY_PRIMARY.covariates + NEGATIVE_CONTROLS`; exactly one
declared specification has an empty suffix; all suffixes distinct.

### 15.11 The module boundary, and the public surface is six names

Six public names in `sensitivity.py` and **eight in `bootstrap.py`**, both asserted by scan (§3.2,
§5.4). `sensitivity.py` names no `outcome.py` private, does not name `derive_cohort` or
`_core_above_median`, and imports neither `statsmodels` nor `scipy` — the last extending Stage 6
§12.7's scan rather than duplicating it. `benjamini_hochberg` and `e_value` name no outcome, family,
covariate or centre. And the rename is asserted complete: `bootstrap` has no attribute `_bucket`,
`_collect` or `_diagnostics`, so a forwarding alias left behind fails the test rather than living on
(§5.4).

### 15.12 Stage 11 adds nothing, refits nothing, and breaks nothing earlier `[data-gated]`

The whole Stage 1-10 pipeline runs with and without this stage, and every landed number and the first
**31** audit entries are byte-identical. This is how `outcome.py`'s estimator is asserted unchanged
(§13) — by number rather than by diff, because four things in it do change: three audit step names
and the extraction of `estimation_population` (§8.5). The extraction is the one that could move a
number and does not, and this test is the reason that claim is checkable: `primary` calls the
extracted function instead of computing the mask inline, so if the extraction were not equivalent
every downstream number would move at once and none of them does.

### 15.13 The fifteen audit entries, their step names, and the ledger

Fifteen `model` entries, the exact step names of §10, the ledger at **46**, no new `data.KINDS`
entry, `case_ids` empty on all fifteen, and the six rendering rules of §10 — including that
`e_value_primary`'s table carries a worst residual |SMD| (§7.6), that it carries **one row per
`C.MRS_THRESHOLDS` entry with the control-arm `P(Y <= k)` and its 15%-band verdict** (§7.5, §12.2
item 4a) read off `est.cumulative` and equal to it row for row, and that `multiplicity`'s carries the
absent-outcome row (§6.5). The step-name collision test builds one `Audit` through the primary and
the arm and asserts the six suffixed names are disjoint from the six unsuffixed, and that
`audit.entry("model", "propensity_fit").n` is the **primary's** `in_model`.

### 15.14 `[data-gated]` — the workbook, end to end

All four entry points run on the workbook with **zero** `C.SchemaError`; the arm's and subgroups'
surviving counts are above `ci_min_draws()`; the four re-roled rows are exactly
`C.NEGATIVE_CONTROLS`; the two `in_model` masks are equal and `differs_only_in_specification` is
`True`. No number this test produces appears in this document (§1).

### Coverage map

```
  benjamini_hochberg        15.4
  e_value                   15.6
  multiplicity              15.1, 15.4, 15.5, 15.13, 15.14
  e_value_primary           15.1, 15.6, 15.13, 15.14
  subgroups                 15.1, 15.7, 15.8, 15.9, 15.13, 15.14
  full_covariate            15.1, 15.3, 15.9, 15.13, 15.14
  _arm_replicate            15.3, 15.8, 15.9
  _subgroup_fit             15.7, 15.8
  _subgroup_replicate       15.8, 15.9
  _assert_arm_inputs        15.1
  _tested_subgroup          15.3, 15.8
  _intervals                15.3, 15.8            both callers, and the p rule is the difference
  propensity._fit           15.10, 15.12, 15.13   [in test_propensity.py]
  propensity.fit_full       15.3, 15.13, 15.14    [in test_propensity.py]
  balance._role             15.3, 15.14           [in test_balance.py]
  outcome.estimation_population  15.7, 15.12      [in test_outcome.py]
  bootstrap R9              15.2                  [in test_bootstrap.py]
  bootstrap.bucket/collect/diagnostics  15.8, 15.11   [in test_bootstrap.py]

  THE MAP IS THE COUNT. There is no total, because a total not derived from the rows
  above cannot be checked against them.

  UNREACHABLE ON v7, SPECIFIED ANYWAY, AND REACHED BY A FIXTURE INSTEAD:
    G8, the constant-interaction route            0 of 2000        15.8
    an absent family p (below ci_min_draws)       0 of 7 outcomes  15.5
    m_used == 0                                   unreachable      15.5
    H4, a raw p below its own key's floor         unconstructible  15.1, 15.4
    H5, no `beta` interval at all                unreachable      15.1
    H7, a point estimate outside its own interval 0 on the draws   15.6
    H6, a non-finite limit reaching e_value       unreachable      15.6
    H8, a POINT-ESTIMATE subgroup fit failing     0 of 2 subgroups 15.8
    e_limit is None (no interval at all)          unreachable      15.6
```

### Failure modes

| Failure | Test | Error handling | Visible or silent |
|---|---|---|---|
| `run` given the arm's `Propensity` | 15.2 | `SchemaError` R9 | **was silent — 26 intervals, no raise** |
| Two fits, colliding audit step names | 15.13 | `spec.step` | **was silent — first-match returns the primary's** |
| `_role` left reading `PS_COVARIATES` | 15.3 | AST scan | **was silent — the arm's adjusted covariates read as failed controls** |
| The two `draws` dicts merged | 15.9 | prohibition + test | **silent — `Draws` cannot notice** |
| `A × S` dropped as constant | 15.8 | `FitError` G8 | **was silent — a two-column fit converges** |
| A separated subgroup level | 15.8 | `FitError` G9 | **was silent — `polr` converges and returns** |
| BH without the running min | 15.4 | oracle | silent — plausible numbers, wrong ranking |
| BH without the unsort | 15.4 | oracle | silent — every p on the wrong outcome |
| The E-value from a null-crossing limit | 15.6 | the 1.0 rule | **silent — a large number reading as robustness** |
| `min(|lo|,|hi|)` for the nearest limit | 15.6 | sign selection | **silent until an interval excludes the null** |
| A `Propensity` built without `spec` | 15.10 | `TypeError` | visible |
| A p on an arm key | 15.3 | `_intervals(_, lambda: False)` | silent |
| A third specification, un-spec'd | 15.10 | registry pin | visible |
| No `beta` interval, so the primary's p is unreported | 15.1 | `SchemaError` H5 | **was a bare `KeyError` from a reported field** |
| The raw-p floor checked against `C.N_BOOT` | 15.1, 15.4 | `SchemaError` H4, per key | **silent — the check is weaker than the floor it checks** |
| A point-estimate subgroup fit failing | 15.8 | `SchemaError` H8 | **was a `FitError` escaping a function whose contract says it cannot** |
| `_diagnostics` on the subgroup run | 15.11 | no `Diagnostics` field | **would have printed "0 replicates reached it" for 1750 that did** |
| A forwarding alias left on `_bucket` | 15.11 | attribute scan | silent — two spellings of one function |

---

## 16. Known gaps carried forward

1. **A 12.5% drop rate exists, it is selection on the estimate, and nothing prescribes what to do
   about it.** `unknown_onset` loses 250 of 2000 replicates — 2 propensity, 9 O6, 239 G9 (§8.4).
   G9's condition is on the estimated coefficients, so the survivors are selected on the estimate and
   the percentile limits are quantiles of a truncated distribution: **narrower**, not merely thinned.
   `../out/questions_for_the_pi.md` carried the standing question — *"whether a drop rate above some
   level invalidates the interval"* — deferred as premature because no drop rate existed. It exists,
   and Question 4 item 1a now puts the mechanism and the three reporting options with it.
   **Trigger: now.** → `TODOS.md`, `questions_for_the_pi.md`.
2. **The witnessed-onset level rests on five treated patients and the subgroup was retained on a
   different number.** [§13]'s amendment of 2026-08-10 justified keeping it on 33 witnessed of 126
   records; the number that governs estimability is 5 of 41 treated, in the ATO population, and
   nobody had measured it. Not a reason to withdraw the subgroup — that decision was taken before any
   subgroup estimate existed and must not be revisited on one — but the PI should see the number.
   **Trigger: with item 1.** → `questions_for_the_pi.md`.
3. **The interaction test's power is not measured.** "Hypothesis-generating" is a word; the number
   that would make it mean something is the power of the `gamma` test at this cohort's shape, in
   Stage 10 §8.4's style with a Monte Carlo standard error. It would also close Stage 10 §16 item
   4's trigger — a coverage design carrying centre-like near-determination — for one parameter.
   **Cost:** of the order of Stage 10 §8.4's two designs. **Trigger:** before the subgroup table is
   read as evidence of no interaction. → `TODOS.md`.
4. **The arm's drop rate is seven times the primary's and only the balance side of that was
   predicted.** [§13]'s amendment predicted the Firth deviation growing with the design; §5.3
   measures it on balance and §5.4 measures it on estimability, and no section of the plan
   anticipated the second. **Trigger:** a third specification. → `TODOS.md`.
5. **`Primary` carries no specification, so two `Primary` objects are structurally identical.**
   `Arm` binds each to its `Propensity` and `__post_init__` checks the pair, but the binding lives
   outside `outcome.py` by the roadmap's choice. **Trigger: the second sensitivity arm** — with three
   `Primary` objects circulating the field should move onto `Primary`. → `TODOS.md`.
6. **`bucket`, `collect` and `diagnostics` are public now and Stage 12 will be their third caller.**
   Stage 10 states its public and private counts in four places and this stage moves all four
   (§13); a fifth statement added later can drift from them.
   **Trigger:** the next change to `bootstrap.py`'s surface. → `TODOS.md`.
7. **`_intervals` duplicates `bootstrap.run`'s draws-to-intervals loop and is not shared with it.**
   One floor rule, one percentile call and one p rule now exist in two modules: `bootstrap.py:1105-1111`
   and `sensitivity._intervals` (§3.2). They are deliberately not refactored together — Stage 10's
   numbers are asserted and §0.2 forbids touching what it refit — so the divergence risk is real and
   filed rather than denied. **Trigger:** a third caller, or any change to `C.ci_min_draws` or
   `C.PERCENTILE_METHOD`, at which point the loop moves into `bootstrap.py` as a ninth public name
   and both callers read it. → `TODOS.md`.
8. **`Replicate` cannot describe a body with more than one `polr` fit.** `n_alpha` and
   `polr_iterations` are scalars, so the subgroup body sets them to `None` and `Subgroups` carries no
   `Diagnostics` (§3.1, §8.5) — which loses the per-subgroup cutpoint distribution that [§16]'s
   constant-shift statement would want per level. Not a defect of this stage: it is a shape
   `bootstrap.py` fixed at Stage 10 for a one-fit body. **Trigger:** the first stage that wants
   diagnostics from a multi-fit body — [§14a] is a candidate — at which point the two fields become
   `dict[str, int | None]` keyed by fit. → `TODOS.md`.

**Two items close.** The stale `data.py:255` citation at `propensity.py:413` (§4.3), whose trigger
was the next commit editing that file. And `TODOS.md`'s parallelise-the-bootstrap trigger — *"the
[§13] sensitivity suite, which is five more `run` calls and turns 145 s into something over ten
minutes"* — is **falsified twice**: there is one arm and not five, and it is `replicates` and not
`run`. Measured, this stage adds **137.6 s** to a 130-145 s pipeline. The item stays open with a
rewritten trigger, because Stage 12's larger population is still ahead of it.

---

## 17. What Stage 11 deliberately does not decide

Each is **PI-reversible** and each names what reversing costs.

| Decision | §  | Reversing costs |
|---|---|---|
| The subgroup model is one pooled weighted proportional-odds fit with an interaction | 8.1, 8.2 | The whole deliverable; (a′) is the named alternative and its assumption advantage is stated |
| `S` enters the linear predictor and no [§6] covariate does | 8.2 | Not reversible without abandoning the interaction |
| `gamma`'s p is not corrected; the count of tests is printed instead | 8.8 | Two numbers, no code path |
| BH is computed on the safety family, carrying [§10]'s descriptive-only label | 6.6 | Four numbers |
| `m` is the number of tests performed, not the declared family size | 6.5 | Seven numbers; unreachable on v7 |
| The bound covers the reported quantities and not `\|delta\|` | 8.4 | Nothing on v7 — measured identical — and it is recorded as principled rather than load-bearing |
| The subgroup stratum stays `centre`, not `centre × subgroup` | 8.1 | `resample` is already general in its stratum. Declined because fixing a random patient attribute's margin conditions `gamma`'s null, and six strata over 92 records puts cells at three |
| `C.SEED` is reused, so all three runs draw identical frames; the pairing is asserted and never used | 5.4, 8.5 | Nothing; it is a property, not an input |
| A point estimate outside its own interval raises rather than reporting 1.0 | 7.3 | One branch; unreachable on the workbook's draws |
| No E-value beyond the primary; no `RR = OR` bracket; no p on the level odds ratios | 7.5, 8.5 | One number each |
| The six `RD_k` are not reported per subgroup level | 8.6 | Twelve intervals |

**And one thing this stage records rather than decides.** The two questions
`../out/questions_for_the_pi.md` deferred with *"ask after Stage 10"* are now both answerable.
Whether the ordinal fit should be penalised: Stage 10 measured the drop rate a penalised fit would
reduce at **zero on the primary**, which was the answer that made the question moot — and this stage
produces a drop rate of **12.5% on a subgroup level**, which makes the same question live again for
a quantity [§13] labels hypothesis-generating and [§8] does not cover. It is put, not answered.

---

## 18. NOT in scope for Stage 11

| Considered | Why deferred |
|---|---|
| The other five [§13] sensitivity rows | Deferred at DECISION 4, which promoted exactly one. The **without-centre** row is deferred *deliberately*: dropping centre removes the near-separation at its source and would balance the clinical covariates nearly exactly while abandoning adjustment for the variable that most nearly determines treatment — 30 of 41 bridged at HUG against 2 of 31 at Lugano — so it is very likely a *more* confounded estimand and agreement with the primary must not be read as reassurance |
| A propensity model without centre, as a diagnostic rather than an estimate | Same row, same reason. A "diagnostic" that produces an odds ratio is an estimate |
| An interval on the E-value | Both E-values are deterministic functions of numbers that already carry intervals (§7.4) |
| Bootstrapping the adjusted p | BH is a deterministic function of seven p-values that already exist. A bootstrap of it would resample the *outcomes*, which is a different procedure |
| Subgroup estimates of the seven binary outcomes | Sixteen interaction tests on 41 treated patients, and a multiplicity [§13] is not told about. [§13]'s sensitivity clause is explicitly *"on the primary outcome"* and Stage 9 §12.2 hands this stage only the family partition and the primary's `beta` |
| A third subgroup | Withdrawn by [§13]'s amendment of 2026-08-10, before any subgroup estimate was produced, and its two constants were deleted rather than commented out |
| Parallelising the two new runs | §16. 137.6 s does not justify the determinism risk, and Stage 10 §4.3's argument against per-replicate seeding is unchanged |

---

## 19. What already exists, and what to lift

| From `pilots/` | Status |
|---|---|
| `e_value`'s core, `rr + sqrt(rr*(rr-1))` with `if rr < 1: rr = 1/rr` (`analysis.py:745-748`) | **Lifted.** It is VanderWeele-Ding's bound and the away-from-null inversion, correct as written. It is the only thing in this file lifted as code |
| `e_value_ci`'s docstring argument for the 1.0 rule (`analysis.py:750-758`) | **Lifted as an argument, verbatim in substance.** It is the roadmap's own reasoning and it is right |
| `e_value`'s RD → RR at the observed control risk (`analysis.py:742-744`) | **Not lifted — wrong estimand, and unavailable rather than declined.** There is no control risk of an ordinal outcome and no event (§7.5) |
| `e_value`'s clip `min(max(p_control + rd, 1e-6), 1 - 1e-6)` (`analysis.py:742`) | **Not lifted.** It repairs an out-of-range implied risk into a sentinel and returns a finite, reportable E-value from it — a fallback producing a number, which invariant 5 forbids |
| `e_value` / `e_value_ci` returning `np.nan` (`analysis.py:740, 759`) | **Not lifted.** `nan` prints as `nan` in a manuscript table. §7.4 raises instead |
| `e_value_ci`'s `np.sign(rd_limit) != np.sign(rd)` (`analysis.py:761`) | **Not lifted, and it is the sharpest single reason §7.3 has a raise.** It conflates *the interval spans the null* with *the estimate is outside its own interval*, and reports the second as the first |
| `add_fdr`: one `multipletests` over `family.isin(["secondary","safety"])` (`analysis.py:789-792`) | **Not lifted, and this is the important one.** One call over the union is BH at m = 7, not [§13]'s "within the secondary and safety families" at m = 3 and m = 4. A substantive departure from the plan, and the exact error §6.1 exists to prevent |
| `add_fdr`'s `res["outcome"].map(lambda o: outcomes[o]["family"])` (`analysis.py:787`) | **Not lifted.** A second computation of the family partition; `outcome.py:757-761` says why that makes the correction wrong |
| `subgroup_analysis`'s comment *"Interaction on the logit scale, weighted by the full-cohort ATO weights"* (`analysis.py:821`) | **Lifted — the right instinct, and the one thing this function gets right.** It is §8.2's weight rule |
| `subgroup_analysis`'s docstring, *"39 treated patients split several ways"* (`analysis.py:801-803`) | **Lifted as framing.** §8.1 makes it a number: five, at witnessed onset |
| `if sub[C.TREATMENT].nunique() < 2 or len(sub) < 20: continue` (`analysis.py:811-813`) | **Not lifted, twice.** An undeclared inclusion rule, and a `continue` that omits a subgroup level **without counting it** — [§10]'s "dropped and counted" with the second half missing |
| `r = bootstrap(sub, covars, spec, ...)` — refitting the propensity model within the level (`analysis.py:814`) | **Not lifted.** This is §8.1's option (b) and the estimand argument is [§8]'s own |
| `sm.GLM(..., freq_weights=w).fit(cov_type="HC0")` and its Wald test (`analysis.py:830-837`) | **Not lifted, four ways.** A binary model where the primary is ordinal; a sandwich SE treating `e` as known, which Stage 8 §5.6 and Stage 10 §18 both refuse; unpenalised where [§7] prescribes Firth; and `except Exception: p_int = np.nan`, a fallback writing `nan` into a reported cell |
| `pd.get_dummies(series.astype(str), drop_first=True)` (`analysis.py:827`) | **Not lifted.** The reference level comes from string sort order, so the sign of `gamma` depends on how a label is spelled (§8.2) |
| `n_boot=max(200, n_boot // 4)` for subgroups and sensitivity (`run_all.py:140, 188, 193`) | **Not lifted.** A reduced replicate count with no stated reason and no record in the output. `C.N_BOOT` is prespecified and this stage uses it |

Net: **one function body and one comment.** Everything touching an estimand, a family partition, a
percentile, a reference level or a failure path is a rewrite — the same verdict Stages 6 and 10
reached about the same file.

### 19b. Validation against a reference implementation

**One oracle is added and it is `statsmodels`, in the role `pyproject.toml` already assigns it**:
`multipletests(..., method="fdr_bh")` checks `benjamini_hochberg` (§6.4). It is the strongest kind of
oracle available here — a widely used independent implementation of a published algorithm, over
inputs that can be enumerated — and it is why §6.4 refuses to *import* it.

**No R oracle is added, and the reason differs per deliverable.** `benjamini_hochberg`: R's
`p.adjust(method = "BH")` would agree with statsmodels and add nothing a second Python-side check
does not already give. `e_value`: the `EValue` R package computes the same closed form, and the part
of §7 that could be wrong is not the formula — it is the null-spanning rule, the limit selection and
the raise, none of which any package implements because they are this plan's. `subgroups`: the fit is
`model.polr`, which Stage 8 already validated against `MASS::polr` and `ordinal::clm`
(`tests/reference/polr_clm.R`); a three-column design is not a new estimator. `full_covariate`: it is
`propensity.fit` and `outcome.primary`, both already oracle-checked against
`tests/reference/ato_psweight.R` and `polr_clm.R`, over a wider covariate list.

**What 19b does not establish.** That the subgroup *estimand* is the right one. Nothing external
validates §8.1's choice; the argument is [§8]'s own text and §17 records it as reversible.

### 19c. The literature check on the E-value's standing — 2026-08-27

§7.5 opened with a claim about the *literature* rather than about this code — *"there is no published
bounding factor for a common odds ratio from a proportional-odds model"* — and neither the probe round
nor the first engineering review checked it. It is load-bearing: it is why the E-value is labelled a
heuristic rather than a bound, and if it were false §7 would change. **Checked, and it holds.** Four
findings: two sharpened §7.5 rather than merely confirming it, one corroborated a §7.3 rule that had
been argued from the formula's shape alone, and one is the confirmation itself.

| What was checked | Found |
|---|---|
| Does the reference implementation cover an ordinal or proportional-odds outcome? | **No.** CRAN `EValue` 4.1.4 (2026-05-07; Mathur, Smith, Ding, VanderWeele) ships `evalues.RR`, `evalues.OR`, `evalues.HR`, `evalues.RD`, `evalues.MD`, `evalues.OLS`, `evalues.IC` — read from its own contents. No ordinal function, no proportional-odds function. That is the strongest available form of a negative claim: the canonical implementation, authored by the method's authors, does not do it |
| Is there a published extension? | **None found.** A methods-forum question asking exactly this — how to compute an E-value for a common odds ratio from an ordinal logistic regression — states in its opening post that no formula appears to exist, and the searches surfaced no derivation for a proportional-odds bounding factor |
| What licenses `sqrt(OR)`, and where is its boundary? | **VanderWeele 2017, Epidemiology 28(6):e58**, *"On a square-root transformation of the odds ratio for a common outcome"* — the citation `evalues.OR` carries. Its boundary is encoded in the implementation as a **prevalence**: `rare = 1` under 15% at end of follow-up, `rare = 0` over 15%, and the square-root conversion applies only in the second case |
| Which limit does the reference implementation take? | **The one nearer the null**, in as many words: `evalues.OR` *"returns a data frame containing point estimates, the lower confidence limit, and the upper confidence limit on the risk ratio scale ... as well as E-values for the point estimate and the confidence interval limit closer to the null."* This is **independent corroboration of §7.3's sign-based selection**, which was argued from the formula's shape rather than from a source. It does not corroborate the 1.0 rule or the raise (§7.3), and §19b already says those are this plan's |

**Two things this changed.** The approximation is now cited to VanderWeele 2017 wherever it is stated,
in `EValue.approximation` itself and not only in a spec section (§3.1, §7.5) — the conversion is his
and attributing it to this plan's roadmap overstated whose choice it is. And §7.5's *"implicitly
assumes a baseline near 0.5"* becomes the documented 15% boundary, which turns a rhetorical caveat
into a **check against a table the stage already computes**: `Primary.cumulative` carries the
control-arm `P(Y <= k)` per threshold, so §12.2 item 4a makes Stage 14 name the thresholds outside the
band instead of asserting that some are.

**What was read from a primary source, and what was not — because a spec that cites a boundary owes
the provenance of the number.** The 15% threshold, the function list and the "limit closer to the
null" sentence were all read from the `EValue` 4.1.4 manual itself. **Tighter bands are reported in
secondary summaries of the method — that the conversion works "fairly well" for outcome probabilities
in 0.1-0.9 and "very well" in 0.2-0.8 — and those were NOT verified against VanderWeele 2017 or the
technical-considerations paper.** They are therefore not used anywhere in this document, and §7.5 and
§12.2 item 4a pin the check to 15%, which is the number the shipped implementation encodes. If a
manuscript wants the tighter band it needs the primary source read first.

**Why no cross-model review pass accompanies this, measured rather than asserted.** This document is
**208 KB** and the review path truncates to **30 KB** — 14% — ending inside §3.1's `Arm` dataclass,
before §4, §5, §6, §7 and §8. Sections 6 through 9 alone are **60 KB**, so even the decision block
does not fit in one pass. A scoped pass over §6 and §7 would fit and is still declined: those are the
two deliverables that already have oracles — `multipletests` for the step-up and a closed form for
the E-value (§6.4, §19b) — so a second model has nothing to check that a running oracle does not.
The parts with no oracle are §8.1's estimand and §8.4's truncation, and both are put to the PI
(§16 items 1 and 2, `../out/questions_for_the_pi.md` Question 4) rather than to a model.

**What 19c does not establish.** That `sqrt(OR)` is *appropriate* for a cumulative-odds shift. It
establishes that the conversion has a named source and a documented prevalence boundary, and that
nobody has published a bound for this estimand — which is exactly the situation §7.5 describes and
§17 records. No literature can settle it, because the quantity is not one the literature covers.

---

## 20. Implementation tasks

- [ ] **T0 (P1)** `config.py`: §13's two fences with their comment blocks; `test_config.py`'s four
      assertions (§15.10). **Gates everything** — no module compiles against `Specification` or the
      two new buckets until it lands.
- [ ] **T1 (P1)** `propensity.py`: `_fit`, `fit`, `fit_full`, `Propensity.spec`, the eight privates,
      the four step names, the docstring at `:459-463`, and the `data.py:271` citation.
      `test_propensity.py` parametrised over the two entry points. **Gates T2, T5, T7.**
      **The no-default field breaks 17 test sites and the count is given so it is not discovered on
      the third one** (§21): `Propensity` is constructed directly at **5** sites in
      `test_balance.py` and **12** in `test_outcome.py`, and every one raises `TypeError` the moment
      `spec` is appended without a default — which is §4.3's intent, not an accident, because the
      alternative is an arm-derived `Propensity` that role-labels as the primary's. Each is a
      one-line fix adding `spec=C.PROPENSITY_PRIMARY`. **The 11 `dataclasses.replace` sites need no
      edit at all** — `replace` carries the field, which is why `test_bootstrap.py:192-194`'s pattern
      is the safe way to bend a `Propensity` (§4.3).
- [ ] **T2 (P1)** `balance.py`: `_role`, `_table`, `Balance.spec`, the three `Final` role constants,
      `_smd_detail`'s two false sentences; `test_balance.py`'s scan for `PS_COVARIATES` (§15.3).
      Depends on T1.
- [ ] **T3 (P1)** `tests/fixtures_stage11.py`: §15.0's five fixtures. **Gates T6, T7, T8** — every
      pin in this document is measured against these, so nothing asserting a number precedes it.
      `near_separated_subgroup_frame` first: §8.4 is unfalsifiable without it.
- [ ] **T4 (P2)** `sensitivity.py`: `benjamini_hochberg`, `e_value`, and §15.4 and §15.6. Pure
      arithmetic over arrays; needs only T0 and T3.
- [ ] **T5 (P2)** `bootstrap.py`: R9; `_bucket -> bucket`, `_collect -> collect`,
      `_diagnostics -> diagnostics` with **no aliases kept**; the two surface counts in all four
      places Stage 10 states them; §15.2, §15.8's scan, §15.11's attribute scan. Depends on T1.
      **The rename is one mechanical edit and the suite is the proof it was complete** — a
      half-swept call site is an `AttributeError` at import, not a wrong number. **The sweep is 15
      lines and the count is given so "complete" is checkable rather than felt** (§21): three `def`
      lines in `bootstrap.py`, plus **12** call sites — **5** in `bootstrap.py` and **7** in
      `test_bootstrap.py`, of which `_bucket` carries 8, `_collect` 3 and `_diagnostics` 1.
      `_diagnostics_table` is a DIFFERENT name and is not renamed; a sweep matching `_diagnostics`
      without anchoring the open paren will rewrite it wrongly.
- [ ] **T6 (P2)** `sensitivity.py`: `multiplicity` and `e_value_primary` with their audit entries,
      **H3-H7**, and `EValue.n_draws`; §15.1, §15.4's H4 case, §15.5, §15.6, §15.13. Depends on T4.
- [ ] **T7 (P3)** `sensitivity.py`: `_subgroup_fit`, `_tested_subgroup`, `_subgroup_replicate`,
      `subgroups`, **H8 and H9**; §15.7, §15.8. **The G8 assertion between `design` and `polr` comes
      before the fit, G9 comes after it, and the point-estimate `FitError` becomes H8** — §8.3's and
      §8.4's placements are each a failure if moved, and §8.5's H8 is the one the replicate loop
      cannot absorb. Depends on T9 for `estimation_population`.
- [ ] **T8 (P3)** `sensitivity.py`: `_intervals`, `_arm_replicate`, `_assert_arm_inputs`,
      `full_covariate`, **H1, H2 and H10**, and `Arm.diagnostics`; §15.3, §15.9. Depends on T1, T2,
      T5. `_intervals` lands here and T7 consumes it, so it is written once and reviewed once.
- [ ] **T9 (P2)** `outcome.py`: the three step names through `ps.spec.step`, and the extraction of
      `estimation_population` with `primary` calling it; §15.7's equivalence assertion, §15.12,
      §15.13. **Promoted from P3 and moved ahead of T7**, which reads the extracted function.
- [ ] **T10 (P3)** `model.py`: the one banner clause (§8.3). Documentation only, and it lands with
      T7 so the sentence and the thing it permits arrive together.
- [ ] **T11 (P3)** §15.14's `[data-gated]` end-to-end test, marked slow.
- [ ] **T12 (P3)** `implementation_roadmap.md`, `TODOS.md` and `questions_for_the_pi.md` per §22 and
      §16. **Lands with the spec commit, not with the implementation.**

`Lanes:` `T0` first and alone. Then `{T1 → T2}`, `{T3 → T4}` and `{T9}` in parallel — T9 touches only
`outcome.py` and is on nobody else's critical path but T7's. Then `{T5}`, `{T6}`. Then `{T8}`, then
`{T7 → T10}`, then `T11`. **T8 before T7** because `_intervals` is T8's and T7 reads it; the earlier
ordering had T7 first and would have had one of them writing the loop twice.

### Definition of done

1. `uv run pytest` green, and the four landed ledger assertions unchanged on a run without this
   stage (§15.12).
2. **R9's absence has been seen to produce 26 intervals** (§15.2). A guard whose absence is never
   demonstrated is not priced.
3. **The step-name collision has been seen** — one `Audit`, both specifications, and
   `audit.entry("model", "propensity_fit")` returning the primary's (§15.13).
4. **`benjamini_hochberg` without the running min has been seen to fail** on §6.3's vector, and
   without the unsort has been seen to fail on the same one.
5. **The E-value taken from a null-crossing limit has been seen to produce a large number**, so the
   1.0 rule's value is a measured fact in the suite.
6. **G9 has been seen to reject a fit that converged** (§15.0's `near_separated_subgroup_frame`).
7. The workbook runs end to end with zero `C.SchemaError`; the arm's and both subgroups' surviving
   counts are recorded in the log and exceed `ci_min_draws()`.
8. `implementation_roadmap.md` carries Stage 11's `**Spec:**` line.
9. No shipped module imports `statsmodels` or `scipy` (§15.11).
10. No case identifier, point estimate, interval limit, adjusted p or E-value appears anywhere under
    `specs/`.
11. **The raw-p floor written against `C.N_BOOT` has been seen to accept a p that H4 rejects**
    (§15.4). The check is worthless unless the weaker version is demonstrated weaker.
12. **A point-estimate subgroup fit failing has been seen to raise `C.SchemaError` H8 and not
    `model.FitError`** (§15.8), on the same fixture the replicate loop drops and counts — so the two
    doors into `_subgroup_fit` are seen to behave differently on one frame.
13. **`bootstrap` has no attribute `_bucket`, `_collect` or `_diagnostics`** (§15.11), and
    `outcome.estimation_population(cohort, ps)` equals `outcome.primary(...).in_estimate` (§15.7).
    Both are one-line assertions and both are the kind an incomplete refactor passes without.
14. **`Arm.diagnostics` reconciles and `Subgroups` has no such field**, asserted rather than
    incidental (§3.1, §15.3): the arm's `n_alpha` sums to the replicates that reached
    `outcome.primary`, and the subgroup run's counters come off `Draws`.
---

## 21. Verification record

**This table is normative.** Where it and the prose disagree, **this table is right and the prose is
stale**. Every row was produced by running code against the landed Stages 1-10 on 2026-08-27, at
`C.SEED` and `C.N_BOOT` unless a row says otherwise. Rows that refuted an earlier draft of this
document are kept and marked.

| Claim | Where used | Verified |
|---|---|---|
| The cohort is 93 rows; `in_model` under [§7] is 92 | §3.4, §5.2 | yes, run |
| `in_model` under [§13]'s full covariate set is **92**, and the two masks are **identical** | §5.2, §15.14 | yes, run — **this refuted an earlier draft's claim** that the arm's population must be a strict subset |
| `hypertension`, `hyperlipidemia`, `diabetes`, `smoking` are missing on **0** of 93 rows | §5.2 | yes, run |
| The one covariate-incomplete record is incomplete on `core_ml` and `tmax6_ml`, both of which are in **both** specifications | §5.2 | yes, run |
| Propensity design width 12 → 16; `dropped` is `center_USZ` in both | §3.4, §5.1 | yes, run |
| ESS per arm: control 29.199 → 29.308, treated 30.339 → 28.472 | §5.2 | yes, run |
| Balance table: 19 rows in both, identical labels in identical order | §4.4, §15.3 | yes, run |
| Roles 14 / 4 / 1 → 18 / 0 / 1, and exactly four rows differ | §4.4 | yes, run |
| The four re-roled rows are exactly `C.NEGATIVE_CONTROLS` | §4.4, §15.14 | yes, run |
| \|SMD\| after weighting: `hypertension` −0.139 → +0.026, `hyperlipidemia` −0.301 → −0.082, `diabetes` +0.380 → +0.052, `smoking` +0.045 → −0.051 | §5.3 | yes, run |
| \|SMD\| after weighting: `center = HUG` +0.211 → +0.275, `center = Lugano` −0.154 → −0.199 — **both worse** | §5.3 | yes, run |
| The four primary-specification figures reproduce [§13]'s amendment of 2026-08-24 exactly | §5.3 | yes, run — cross-check on this stage's arithmetic |
| Rows reaching `SMD_THRESHOLD`: 5 → 2, and both survivors are `center` | §3.4, §5.3 | yes, run |
| Worst \|SMD\|: 0.380 (`diabetes`) → 0.275 (`center = HUG`) | §5.3 | yes, run |
| `penumbra_ml` is `excluded [§6]` in both and is the only such row | §4.4, §5.3 | yes, run |
| A second `propensity.fit` over the same cohort records the **same four** `model` step names, and `assess` the **same two** | §4.4 | yes, run — **this changed §4.2 from a tuple to a record** |
| `Audit.entry` is first-match | §4.4 | yes, read (`data.py:273-275`) |
| The arm: 1986 of 2000 surviving on all seven keys; 13 `degenerate_design`, 1 `nonconvergence`, 0 `separation` | §3.4, §5.4 | yes, run |
| The arm's drop rate 0.70% against Stage 10's measured 0.10% for the primary | §5.4 | yes, run |
| The arm's cutpoint counts 1898 / 91 — the same 4.6% Stage 10 measured at 1907 / 91 | §5.4 | yes, run |
| The arm's `in_model` per replicate runs 88 to 93 | §5.4 | yes, run |
| The arm's whole bootstrap: **69.3 s** | §2, §3.4, §5.4 | yes, run |
| Every one of the arm's seven keys is far above `ci_min_draws()` = 40 | §5.4 | yes, run |
| Family sizes are 3 and 4, off `Secondary.by_family()`; members as §6.1 lists them | §3.4, §6.1 | yes, run |
| Stage 10 emits 26 keys, 8 of which carry a p | §3.4, §6.1 | yes, read (`bootstrap.py:431-467`) |
| `bootstrap_p` at B = 1998 can return **1000** distinct values | §3.4, §6.2 | yes, run |
| BH on `(0.04, 0.01, 0.03)` gives `(0.04, 0.03, 0.04)`; unenforced it is `(0.03, 0.045, 0.04)` sorted | §6.2, §6.3 | yes, run |
| BH agrees with `multipletests(method="fdr_bh")` on 6 hand cases and **10 000 / 10 000** random families at `rtol=1e-12, atol=0` | §6.4, §15.4 | yes, run |
| BH permutation invariance **10 000 / 10 000**; validity **10 000 / 10 000** | §6.4 | yes, run |
| The cap at 1.0 binds in **0** of 10 000 families after the running min | §6.2 | yes, run |
| `exp(0.0) == 1.0` exactly | §7.2 | yes, run |
| `exp(percentile_ci(x))` equals `percentile_ci(exp(x))` in 5 of 5 sizes under `inverted_cdf`, and in **0 of 5** under `linear` | §7.2, §15.6 | yes, run |
| `E(1.0) == 1.0` exactly; `E(RR=2) == 2 + sqrt(2)` | §7.3, §15.6 | yes, run |
| `E(OR=4) == E(OR=0.25)` exactly — relative difference 0.0 | §15.6 | yes, run |
| `1/sqrt(v)` equals `sqrt(1/v)` on 5 of 5 tested values, so §7.1's pin is currently **inert** | §7.1 | yes, run |
| `e_value(sqrt(1 + 1e-12)) − 1 = 7.07e-07` | §7.4 | yes, run |
| `POLR_MAX_ABS_BETA` = 14.0 bounds `OR` to 1.2e6 and `E` to **2192.77** | §7.4 | yes, run |
| `rr² = OR ≤ DBL_MAX` always, so the plain form cannot overflow | §7.4 | yes, read |
| Point-estimate subgroup cells: `unknown_onset` 16 / 5 / 37 / 34; `core_above_median` 24 / 24 / 29 / 15 | §8.1 | yes, run |
| The 5 treated at witnessed onset are HUG 4, CHUV 0, Lugano 1, USZ 0 | §8.1 | yes, run |
| `in_estimate & df[S].notna()` equals `in_estimate` at 92 for both subgroups | §8.2, §15.7 | yes, run |
| Point-estimate fits: `unknown_onset` \|β\| 1.719, \|δ\| 0.519, \|γ\| 1.841, \|β+γ\| 0.123, 4 iterations; `core_above_median` 0.366 / 0.496 / 0.562 / 0.196, 3 iterations; `dropped` empty in both | §8.2, §8.3 | yes, run |
| Subgroup drops: `unknown_onset` 2 propensity + 9 O6 + **0** constant-interaction; `core_above_median` 2 propensity + 0 + 0 | §8.3 | yes, run — **this refuted an earlier draft's claim** that the constant-column route is the reachable one |
| The 9 O6 replicates are **exactly** the 9 that drew zero treated patients at witnessed onset | §8.3 | yes, run |
| The 2 propensity failures are the same in both subgroups and are Stage 10's measured two | §8.3, §8.5 | yes, run |
| O6 is a rank check raising `FitError`, and `C.FAILURE_BUCKETS["O6"]` is `degenerate_design` | §8.3, §9 | yes, read (`model.py:837-843`, `config.py:344-347`) |
| `unknown_onset`: max \|β\| 18.368, \|γ\| 19.296, \|β+γ\| 2.383, \|δ\| 5.257 | §8.4 | yes, run |
| `unknown_onset` \|reported\| median 2.317, p90 16.187, p99 18.176 | §8.4 | yes, run |
| **239** of 1989 fitted replicates reach `POLR_MAX_ABS_BETA` on a reported quantity; **0** on \|δ\| | §8.4 | yes, run — **the guard is not inert, which an earlier draft assumed it would be** |
| Bounding everything and bounding the reported quantities give the **same 1750** survivors | §8.4 | yes, run |
| Those 239 carry 1 to 9 treated at `S = 0`, of which only 86 carry 2 or fewer | §8.4 | yes, run |
| `polr` runs 10 to 15 iterations in those 239, against 3-4 typical | §8.4 | yes, run |
| `unknown_onset` surviving 1750 of 2000 — a **12.5%** drop rate; `core_above_median` 1998 | §3.4, §8.4 | yes, run |
| `core_above_median`: max \|reported\| 5.326, max \|δ\| 5.520, **0** over bound | §8.4 | yes, run |
| `core_above_median`'s cutpoint counts are **1907 / 91** — identical to Stage 10's primary, over the same 1998 replicates | §8.5, §15.9 | yes, run — **and it is a probe-round measurement and NOT a shipped counter**: `Subgroups` carries no `Diagnostics` (§3.1), so this number is not recoverable from the log and this row is where it exists |
| `bootstrap.replicates` declares `body: Callable[[pd.DataFrame], object]` and calls it with one argument, and `run` binds its four-argument body with a lambda | §5.4, §8.5 | yes, read (`bootstrap.py:347-348`, `:1096-1100`) — **this stage's engineering review added the call-site fences it forced** |
| `bootstrap.Diagnostics` carries `n_alpha` and `polr_iterations` as **scalars per replicate**, so a two-fit body cannot fill them | §3.1, §8.5 | yes, read (`bootstrap.py:259-260`, `:148-149`) |
| `model.design` returns `X` as `.astype(float)` with the passed covariates as its columns in the passed order, none of the subgroup design's three being in `C.CATEGORICAL` | §8.3, §15.7 | yes, read (`model.py`, `design`) — which is what makes G8's exact-tuple predicate the same check as §15.7's |
| `outcome.primary` binds the [§11] mask as `ps.in_model & df[C.PRIMARY_OUTCOME].notna()` at one place | §8.5 | yes, read (`outcome.py:634`) |
| `bootstrap._diagnostics` takes `collected` alone, so §5.4's `bootstrap.diagnostics(collected)` is the landed signature | §5.4, T5 | yes, read (`bootstrap.py:718`) |
| `Propensity` is constructed directly at **17** test sites — 5 in `test_balance.py`, 12 in `test_outcome.py` — each of which raises `TypeError` once `spec` lands without a default; the **11** `dataclasses.replace` sites do not | §4.3, T1 | yes, counted |
| The three privates going public are 15 lines: 3 `def`s plus **12** call sites, **5** in `bootstrap.py` and **7** in `test_bootstrap.py` (`_bucket` 8, `_collect` 3, `_diagnostics` 1); `_diagnostics_table` is not among them | §5.4, T5 | yes, counted |
| `bootstrap_p`'s floor is `1.0/(len(draws)+1)` over the draws it was given, so the floor is per key and not `C.N_BOOT`'s | §6.1 | yes, read (`bootstrap.py:411`) |
| Stage 10 §9.4 measures p-versus-interval agreement holding on the tail-count-50 fixture under `inverted_cdf` and **failing under `linear`** | §7.2, §15.6 | yes, read (`stage10_bootstrap_engine.md` §9.4) — **this corrected §7.2's claim that the fixture is where route 3 can fail** |
| `data.KINDS` has been nine for **four** stages, Stage 10 §10.1 having recorded three | §10 | yes, read — **this corrected §10's "five stages"** |
| Both subgroups, one run: **68.3 s** | §2, §8.5 | yes, run |
| All six subgroup keys exceed `ci_min_draws()`; the worst is 1750 | §8.5 | yes, run |
| Two independent 200-replicate drives give identical `case_id` fingerprints, matching a hand-driven `resample` loop; 200 distinct frames | §5.4, §8.5, §15.9 | yes, run |
| `replicates` calls `resample` outside `body` with no `try` | §5.4 | yes, read (`bootstrap.py:372-373`) |
| `bootstrap._tested("beta")` is `True` | §5.5 | yes, read (`bootstrap.py:456-467`) |
| A `Bootstrap` for the arm's `Propensity` passes R3 and R5 and raises nothing today | §4.5 | yes, read (`bootstrap.py:578`, `_assert_run_inputs`) |
| `outcome.primary`'s design is `(C.TREATMENT,)` alone | §5.1 | yes, read (`outcome.py:639`) |
| `C.OUTCOME_COVARIATES` is bound to the `PS_COVARIATES` object at import | §5.1 | yes, read (`config.py:561`) |
| `C.NEGATIVE_CONTROLS` is computed as `PS_COVARIATES_FULL` minus `PS_COVARIATES` | §5.3, §13 | yes, read (`config.py:571-572`) |
| `test_propensity.py:519-524` asserts `fit`'s signature and that no parameter has a default | §4.1, §4.3 | yes, read |
| `data.KINDS` is nine and has been for four stages | §10 | yes, run |
| The point-estimate pipeline records 30 entries before the bootstrap; Stage 10 measured 31 with it | §10 | yes, run |
| The stale citation at `propensity.py:413` says `data.py:255`; `Audit.record`'s normalisation is at `data.py:271` | §4.3, §16 | yes, read |
| This stage adds **137.6 s** to the pipeline, against `TODOS.md`'s predicted "over ten minutes" for five `run` calls | §2, §16 | yes, run |

### 21b. This document's own code, executed — 2026-08-27

**All fenced Python blocks parse.** As of the probe round there were fourteen, five carrying
executable bodies, and all five were run against the landed modules; the run is reproducible from
this document alone. The engineering review added three more and they were **not** run — the note
below the table is where that is recorded.

| Fence | § | What running it established |
|---|---|---|
| `Specification` | 13 | Constructs; the registry has two; `PROPENSITY_FULL.covariates == PROPENSITY_PRIMARY.covariates + C.NEGATIVE_CONTROLS`; exactly one suffix is empty and both are distinct; **`PROPENSITY_PRIMARY.step("propensity_fit") == "propensity_fit"`**, so no landed step name moves |
| `benjamini_hochberg` | 6.2 | §6.3's vector exactly; the four hand cases; and **10 000 / 10 000** against `multipletests`, permutation-invariant and valid on the same 10 000 |
| `e_value` | 7.1 | `E(1) == 1.0` and `E(RR=2) == 2 + sqrt(2)` bit-exact; symmetric under inversion; monotone over 400 points to `RR = 1e6`; `E` at the `POLR_MAX_ABS_BETA` bound is 2192.77 |
| `_subgroup_fit` | 8.3 | Runs on the workbook's cohort for **both** subgroups: three columns, `dropped` empty, mask equal to `in_estimate` at 92, and the reported quantities inside the bound. **And on a constructed one-level frame it raises `G8` before the fit** — so the assertion's placement is demonstrated and not merely described |
| `_arm_replicate` | 5.4 | 60 replicates over `bootstrap.replicates`; seven keys; `Draws.__post_init__` reconciles on all seven |

Two things this round established that the prose did not have before it:

- **`G8` is a token `bootstrap.bucket` currently raises on.** The degenerate-frame run raised `G8`,
  and classifying it needed `C.FAILURE_BUCKETS` extended by §13's fence first. That is §9's argument
  turned into a demonstration: the token and the map entry must land in one edit.
- **§7.4's near-null figure is on the odds-ratio scale and `e_value` takes a risk ratio.** The first
  draft of §7.4 wrote `E(1 + 1e-12)`, which reads as a risk ratio and is a different number. It is
  now written as `e_value(sqrt(1 + 1e-12))`. Caught by executing the fence.

The remaining nine fences are dataclass and signature declarations; they parse and are checked by
construction, since six of them hold Stage 10's and Stage 6's types.

**Three fences were added by the engineering review round and NOT executed, and that is stated rather
than left to be assumed** (§22.2): `_intervals` (§3.2), `_subgroup_replicate`'s body (§8.5) and the
two `bootstrap.replicates` call sites (§5.4, §8.5). They were written against signatures read from
the landed modules — `replicates`' one-argument body, `Replicate`'s six required fields, `collect`'s
and `diagnostics`' parameters, `Interval`'s constructor order — each of which §21 now records as a
`yes, read` row. **They are the one part of this document whose arithmetic has not been run**, and T5
through T8 are where they first execute. The three fences that were executed before this round —
`_subgroup_fit`, `_arm_replicate` and `benjamini_hochberg` — are unchanged in body by the review
except for `_subgroup_fit`'s G8 predicate, which is two lines and is asserted by §15.7 from both
sides.

### 21c. What was measured and deliberately not stated

The primary's interval limits, its p-value, its E-values, the seven raw and adjusted p-values, the
two subgroups' odds ratios and interaction p-values, and the arm's seven limits **were all produced**
by the probe round. None appears in this document, in any form, including as a range or a sign
(§1). They are in the gitignored log and they are the manuscript's.

---

## 22. What this spec changed elsewhere

| # | Change | Files | The question it settles |
|---|---|---|---|
| 1 | Stage 11 gains its `**Spec:**` line | `implementation_roadmap.md` | Where the specification is |
| 2 | The Accept-when's "and in nothing else" is kept and gains "and in what that specification's completeness costs, both counts printed" | `implementation_roadmap.md` | It is true on v7 and is not true by construction (§5.2) |
| 3 | The subgroup bullet gains the primary-outcome scope and the interaction-test count | `implementation_roadmap.md` | [§13] names no outcome; §18 records why it is the primary alone |
| 4 | Stage 10's "a second `run`" becomes "a second `replicates` body" in §12.1 and §14 item 3 | `specs/stage10_bootstrap_engine.md` | §5.4 — `run` drives `outcome.secondary` and names `propensity.fit` |
| 5 | Stage 10 §3.2's "five names" becomes **eight** and its "twelve privates" becomes nine, in each of the four places it states them; its Definition of done item 2 — *"`propensity.py` has a zero-line diff"* — is marked historical, since running it as a live check would now fail and is not supposed to pass | `specs/stage10_bootstrap_engine.md` | §4.3, §5.4 |
| 6 | Stage 6 §6.4's bullet is marked discharged and points at `fit_full` | `specs/stage6_propensity_and_weights.md` | The promissory note is paid (§4) |
| 7 | Stage 7 §11 is marked fulfilled and §14's deferred question answered | `specs/stage7_balance_and_overlap.md` | The role column moved and the table did not (§4.4) |
| 8 | Two items close, one trigger is rewritten, eight are new | `TODOS.md` | §16 |
| 9 | The 12.5% drop rate and the five treated patients are put to the PI | `../out/questions_for_the_pi.md` | §16 items 1 and 2; the standing "ask after Stage 10" question is now live |
| 10 | Stage 9 §3.2's public-surface count for `outcome.py` gains `estimation_population` | `specs/stage9_secondary_binary_estimators.md` | §8.5 — the [§11] mask becomes a named definition with two callers, which is the rule `weighted_proportion` is already under |
| 11 | Stage 10 §3.2's *"the private list becomes twelve"* becomes nine, in the same four places item 5 touches | `specs/stage10_bootstrap_engine.md` | §5.4 — three privates became public and the count is stated more than once |
| 12 | Question 4 gains item **1a**: the drop is selection on the estimate, so the interval is narrower and not merely thinned; three reporting options are put and the specified one is named | `../out/questions_for_the_pi.md` | §8.4 — the mechanism, which "250 replicates were lost" does not convey |
| 13 | Question 4's third drop route is corrected from 241 to **239**, the two propensity failures having been folded into it | `../out/questions_for_the_pi.md` | §21 is normative and says 2 + 9 + 239 |

### 22.1 The probe round — 2026-08-27

Thirteen probes, four of which changed the document. §21 marks each with
**"this refuted an earlier draft's claim"**. In order of how much they changed:

1. **The subgroup degeneracy route.** The draft specified the constant-interaction column as the
   reachable failure and did not mention rank deficiency. Measured: constant-interaction 0 of 2000,
   O6 rank deficiency 9 of 2000, and the 9 are exactly the empty-treated-cell replicates. The route
   that fires was already caught and already bucketed by a Stage 8 guard; the route the draft
   specified is unreachable here and is kept as a specified-and-unreachable branch (§8.3).
2. **The subgroup bound is not inert.** The draft assumed the G7-analogue would be a
   specified-and-unreachable guard, as the `|δ|` distribution suggested it might be. Measured: 239 of
   1989 fitted replicates on `unknown_onset` reach `POLR_MAX_ABS_BETA` on a reported quantity, and
   the resulting 12.5% drop rate is the largest this pipeline has produced (§8.4, §16 item 1).
3. **The arm's population.** The draft specified a correction to the roadmap's Accept-when. The
   correction is not needed on v7 and the clause becomes a computed field rather than either a
   correction or an assertion (§5.2).

4. **The audit step names collide exactly.** A second `propensity.fit` over the same cohort records
   the same four `model` step names and `assess` the same two. `Audit.entry` is first-match, so every
   programmatic read would have returned the primary's entry while the rendered log looked complete —
   and every existing test would have stayed green, because on a run without the arm they *are* the
   same entry. This is what turned §4.2's payload from a covariate tuple into a named record, before
   any of §4 was drafted.

The nine that confirmed rather than changed are in §21 without a marker. Two are worth naming: the
three runs draw identical frames, which was expected from `replicates`' construction and is now
asserted rather than assumed; and the four primary-specification |SMD| figures reproduce [§13]'s
amendment of 2026-08-24 exactly, which is a cross-check on this stage's arithmetic against a number
written five days before it.

### 22.2 The engineering review round — 2026-08-27

Eighteen findings, all folded, over three passes: the review itself (1-16), then a sole-source audit
that re-checked every folded finding against the file rather than against the review's own claims
(17-18). **Four changed a structure and are named in the Status note above**; the rest were gaps and
internal inconsistencies a sole-source document cannot carry. None changed a measured number, and §21
gained eight `yes, read` rows for the code these passes read.

**Findings 17 and 18 are the ones worth noting as a class**, because both were things said in review
conversation and not written down — a corroboration and a required log row. A document that is the
sole source cannot hold a finding in the reviewer's head, and the only way to catch that is to audit
the file against the claims rather than trusting the claims. Both were found that way.

| # | Finding | What it changed |
|---|---|---|
| 1 | `bootstrap._diagnostics` was left private while both bodies populated `Replicate`'s diagnostic fields, and neither `Arm` nor `Subgroups` had anywhere to put a `Diagnostics` — so §5.4's and §8.5's cutpoint counts had no home | `diagnostics` becomes public (surface 5 → **8**); `Arm` gains the field; `Subgroups` deliberately does not, and §3.1 says why (§3.1, §5.4, §8.5) |
| 2 | `_subgroup_replicate` was a docstring with no body, and `Replicate.n_alpha` is a scalar against this body's **two** `polr` fits | The body is written out; the two cutpoint fields are `None` by decision; `sum_w` and `n_in_model` keep their one consumer (§8.5) |
| 3 | `bootstrap.replicates` takes a **one-argument** body and both bodies here take three; the binding was never written | Both call sites are written out — closed-over keys, `collect`, `diagnostics`, `_intervals` (§5.4, §8.5) |
| 4 | The draws-to-intervals loop was described twice in prose and would have been written twice, differing in one line | `_intervals(draws, tested)`, one helper, with the p rule as its argument (§3.2) |
| 5 | A `FitError` from a **point-estimate** subgroup fit had no specified behaviour, and §11 promised one could not escape | `C.SchemaError` H8, with the argument for why the loop's rule does not transfer (§8.5) |
| 6 | The replicate needed the [§11] mask on a drawn frame, and `outcome.primary` computes it privately at `:634` | `outcome.estimation_population` extracted; §13's `outcome.py` row and §15.12's framing both move (§8.5) |
| 7 | §15.1 required each `C.SchemaError` matched by its identifier and no identifiers existed | The **H series**, ten of them, declared with the note that none is a `FAILURE_BUCKETS` key (§9.1) |
| 8 | `multiplicity` read `intervals["beta"].p` with no precondition on it | **H5**, and the reason the primary is uncorrected rather than unreported (§6.1) |
| 9 | The raw-p range check said `1/(B+1)` without saying what `B` was; against `C.N_BOOT` it is **weaker than the floor it checks** | **H4** pinned to that key's `Interval.n_draws`, with the two-`n_draws` fixture that demonstrates the weaker version passing (§6.1, §15.4) |
| 10 | §7.4 required the draw count recorded when no interval exists and `EValue` had no field for it | `EValue.n_draws`, populated on both branches (§3.1, §7.4) |
| 11 | §7.2 called the tail-count-50 fixture *"the one fixture where route 3 can fail"*; Stage 10 §9.4 measures agreement **holding** there under `inverted_cdf` and failing under `linear`. Route 3's operator was also unpinned | The claim is corrected to what §9.4 measured, and `>=` is pinned with the tie argument (§7.2, §15.6) |
| 12 | `_subgroup_fit`'s G8 predicate (`lost or dropped`) was not §15.7's assertion (the exact column tuple) | One predicate, asserted in both places, with the reason column order is load-bearing (§8.3) |
| 13 | `n_interaction_tests` was the literal 2 while §3.3 computes the keys from the registry | `len(C.SUBGROUPS)` (§3.1) |
| 14 | `FamilyCorrection.absent`'s reason string had no shape, though §6.5 requires the row to carry the draw count and the floor; and `m_used == 0` "says so" nowhere | The string is specified; the zero case is defined as empty dicts plus a full `absent`, with `benjamini_hochberg` not called at all (§3.1, §6.5) |
| 15 | Three doc-internal contradictions: the probe count (11 against 13), the findings count (three against four), and `data.KINDS`' stage count (five against §21's four) | All three resolved to §21's numbers; §21 is normative and said so already |
| 16 | §13's `FAILURE_BUCKETS` fence used `...,` as its elision, which is a set element followed by dict entries and **does not parse** — against §21b's claim that every fence does | The elision is a comment; every fence in this document now compiles, checked by extracting all seventeen and calling `compile` on each |
| 17 | §19c's literature check found that `evalues.OR` takes *"the confidence interval limit closer to the null"* in as many words — external corroboration of §7.3's sign-based selection, which had been argued from the formula's shape alone — and it was reported in review conversation but **not written into the document** | §19c gains the row. A corroboration nobody can find is not one |
| 18 | §12.2 item 4a required Stage 14 to name the thresholds outside the 15% band, and §1 keeps those numbers out of this document — so the **log had to print them and §10's `e_value_primary` entry did not** | The entry gains one row per `C.MRS_THRESHOLDS` entry, read off `est.cumulative`; §10's decision count goes to six and §15.13 asserts the rows equal the source table |

Two of these were found by reading code rather than by reading the document — findings 3 and 12 —
and one of them, 3, is the reason this round happened before implementation rather than during it:
it is not a wrong sentence, it is a missing one, and it is missing in the place where the two runs'
shapes are decided.

### 22.3 What did not change, which is the more important half

Nothing in Stages 1-10's **estimates** moves. `outcome.py`'s three amended lines are audit step
names; §15.12 asserts the numbers by running the pipeline twice. `propensity.fit` and
`balance.assess` keep their signatures. `data.KINDS` stays nine. The [§7] specification, the primary
estimate and [§13]'s deferral of the other five rows are untouched — which is what DECISION 4
established, *before any sensitivity estimate existed and before the [§10] interval existed*, and is
the reason this stage may report a sensitivity arm at all.

---

## GSTACK REVIEW REPORT

| Review | Trigger | Why | Runs | Status | Findings |
|---|---|---|---|---|---|
| CEO | new stage spec | Scope: four deliverables against [§6] and [§13] | 1 | done | Subgroup scope narrowed to the primary outcome; recorded §18 |
| Codex | — | — | 0 | not run | No cross-model second opinion at this stage; `TODOS.md` carries the standing item |
| Eng | new module + seven amended | The seam touches five landed modules, and this document is the sole implementation source | 2 | done | Round 1: R9; the record-not-tuple payload; `replicates` not `run`. Round 2 (2026-08-27): 18 findings over three passes, all folded — 4 structural, the rest sole-source gaps and internal inconsistencies; findings 17-18 came from auditing the file against the review's own claims (§22.2) |
| Design | — | — | 0 | n/a | No user-facing surface |
| DX | — | — | 0 | n/a | No developer-facing tooling |
| Outside Voice | §7.5's claim was about the literature, not the code | A cross-model pass truncates to 30 KB of this document's 208 KB and ends inside §3.1 — it reaches no decision | 1 | done | **Cross-model: declined, measured.** **Literature check: run** (§19c) — §7.5's negative claim holds; the `sqrt(OR)` conversion gains VanderWeele 2017's citation and its 15% prevalence boundary, which became §12.2 item 4a |
| Probe round | before drafting | House process since Stage 9 | 13 | done | Four findings changed the document (§22.1) |
| Fence execution | before landing | §21b | 1 | done | Five executable fences run at the probe round; the review's three added fences are **not** run and §21b says so |

**VERDICT: ready to implement.** The one design decision the roadmap left here is made and its
failure modes are guarded (§4); the two statistical procedures have closed-form oracles (§6.4, §19b);
the subgroup estimand is chosen on [§8]'s own argument and recorded as PI-reversible (§8.1, §17); and
the one number that should worry a reader — a 12.5% drop rate on a subgroup level resting on five
treated patients — is measured, printed beside its interval, and put to the PI rather than absorbed
(§8.4, §16). The second engineering round closed the gaps that mattered for a sole-source document:
both replicate bodies and both call sites are now written out, this stage's `C.SchemaError`
identifiers exist (§9.1), and the two places a scalar diagnostic field could not describe a two-fit
body are decided rather than left to the implementer (§3.1, §8.5).

**One thing an implementer should know before T5:** §21b's three new fences are the only part of this
document that has not been executed. They were written against signatures read from the landed
modules and §21 records each read, but T5 through T8 are where they first run.

NO UNRESOLVED DECISIONS
