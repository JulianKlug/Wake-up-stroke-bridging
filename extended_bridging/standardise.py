"""Stage 12 — [§14a] all-centre standardisation.

One pooled proportional-odds model over all eligible patients at all four centres; standardisation
by g-computation to the two mRS distributions and the quantities [§14a] takes from them; the support
check [§14a] requires and the two sensitivity analyses it prescribes; and percentile intervals from a
patient-level bootstrap stratified by centre.

Section references in brackets are to ``statistical_analysis_plan.md``. The specification for this
module is ``specs/stage12_all_centre_standardisation.md``; nothing here is invented outside it.

**THE ESTIMAND IS AN ATE IN THE ALL-CENTRE ELIGIBLE POPULATION AND IS NEVER THE [§7] ATO.** Every
output of this module carries the different-populations statement, which Stage 14 [§16] is required
to print. The two are not comparable: [§7]'s is an overlap-weighted contrast over the 93 patients at
three treating centres, and this one is an average over the 104 eligible patients at four.

**No propensity model is fitted anywhere here, and that is [§14]'s own first sentence rather than a
simplification** — *"No propensity model is fitted anywhere in §14"*. A propensity model cannot
produce a contrast at a centre with no treated patient, which is why the section exists at all. There
is no ``e``, no ``w``, no ``h``, no ESS and no overlap weight in this module, and ``propensity.py`` is
not imported.

**It does not call ``cohort.build``, and Stage 5 §9 already fixed that.** [§3] restriction 1 removes
the never-IVT centre, which is precisely the population [§14a] exists to include. ``build`` offers no
seam to skip it and must not acquire one, so this module applies restriction 2 itself and takes the
UNRESTRICTED classified frame [§4.1].

The identifying assumption, once, because it governs every function below and a module split would
have put it in neither half [§0.1]: *the conditional effect learned in IVT-using centres transports
to never-IVT centres given X, and that is not testable from these data*. ``support`` is what makes
the extrapolation visible; nothing here tests the assumption, and [§14a] says differences between
centre-specific estimates would not be a test of it either.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Final

import numpy as np
import pandas as pd

import balance
import bootstrap
import cohort
import config as C
import eligibility
import model
from data import Audit, absence_by_column

# --- the guard string, and it is DATA rather than a docstring [§8] --------------------------------
#
# The roadmap's Guard for this stage is one sentence: "`exp(beta)` here is a conditional odds ratio.
# It may be emitted as a model parameter and must never be labelled as the standardised marginal
# effect." A comment saying so would be a label, and labels are what [§16] amendment DECISION 5
# exists because of. So it travels on the record, Stage 14 prints it FROM the record, and §20.8
# asserts every constructed `Standardisation` carries it. A second value can only arrive with a [§14]
# amendment behind it.
_CONDITIONAL: Final[str] = (
    "conditional odds ratio [§14a] — conditional on X, and never the standardised marginal effect")

# The role the four `center` rows of §10.3's baseline table carry. They report a DEFINITION and not a
# difference: the never-IVT group IS USZ, so an |SMD| of 1 on `center = USZ` is the grouping variable
# and a reader who takes it for imbalance has misread the table. They stay in it — dropping a
# declared row is what Stage 7 §4.1's role column exists to avoid — and they carry this instead of a
# [§6] role.
_GROUPING: Final[str] = "grouping"

# The three values `Standardisation.population` takes, one per arm, so a Stage 14 formatter can tell
# the arms apart from the record alone (§20.12) without tracking an estimand-key prefix.
_ALL_ELIGIBLE: Final[str] = "all eligible"
_TREATED_SUPPORT: Final[str] = "treated support"
_RANDOM_INTERCEPT: Final[str] = "all eligible, random centre intercept"

# T10's and T11's tolerances, and they are §20.5's and §20.6's own so that the guard and the test
# cannot disagree. The asymmetry between the two halves of T11 is §7.3's and is not a rounding
# choice: `mrs_0_2` IS `rd[2]` — the same cumulative difference read under [§14a]'s name for it, one
# computation and not two — so its comparison is EXACT and 1e-15 would hide a real defect, while
# `mortality` is genuinely computed a second way, from `P(Y = 6)` in each arm rather than from
# `1 - P(Y <= 5)`, so it agrees with `-rd[5]` to floating point and not to the bit (measured 8.3e-17).
_DISTRIBUTION_TOL: Final[float] = 1e-12
_MORTALITY_TOL: Final[float] = 1e-15


# --- what the stage returns [§3.1] ----------------------------------------------------------------

@dataclass(frozen=True)
class Standardisation:
    """One [§14a] standardisation: one fit, one averaging population, two regimes.

    Frozen and carrying no verdict, following Stage 8 §3's rule that a record reports what it
    computed and a caller decides what it means.

    **`fit` is a union of two types and `population` has three values, because all three arms return
    this record and one of them is fitted by a different estimator.** The alternative — a fourth
    field discriminating `"polr"` from `"polr_ri"` — is declined for §7.3's reason: it would carry
    what the `fit` type already determines, and two things that encode one fact are two things that
    can disagree after an edit. A caller that must branch narrows on the type.

    **`conditional_log_odds` means something slightly different in the hierarchical arm and the
    difference is stated rather than smoothed over.** In the pooled and support arms it is conditional
    on `X`; in the hierarchical arm it is conditional on `X` **and on `b_c`** (§12.1). `measure` stays
    `_CONDITIONAL` for all three — §8's guard is that neither may ever be labelled marginal, and both
    qualify — but Stage 14's caption must carry which conditioning it is.

    **THERE IS NO FIELD CALLED `odds_ratio`, and that is a decision rather than a naming
    preference.** `Primary.odds_ratio` is the [§8] MARGINAL common odds ratio; a Stage 14 formatter
    reaching for the same attribute on two records would print a conditional quantity under a
    marginal label. §8 is the argument and the absent name is the enforcement — a formatter pointed
    at the wrong record gets an `AttributeError` rather than a number.

    **And the reported [§14a] effect is the set of risk differences, not `beta`.** That is [§14a]'s
    own choice, made because the standardised quantities are marginal by construction and are on the
    absolute scale, where non-collapsibility does not arise.
    """

    population: str                              # _ALL_ELIGIBLE | _TREATED_SUPPORT | _RANDOM_INTERCEPT
    n_average: int                               # the averaging population's size — §7.2, [§11]
    n_fit: int                                   # the FIT's population; equal unless §11 changes it
    distribution: dict[int, dict[int, float]]    # arm code -> mRS level -> P, over MRS_LEVELS
    cumulative: dict[int, dict[int, float]]      # arm code -> threshold -> P(Y <= k)
    rd: dict[int, float]                         # RD_k, keyed by MRS_THRESHOLDS — §7.3
    mrs_0_2: float                               # == rd[2]; [§14a]'s name for it — §7.3
    mortality: float                             # == -rd[5]; [§14a]'s name for it — §7.3
    conditional_log_odds: float                  # beta. NOT called an odds ratio, and §8 is why
    conditional_odds_ratio: float                # exp(beta), and it is CONDITIONAL — §8
    measure: str                                 # the guard string, always _CONDITIONAL — §8
    fit: model.PolrFit | model.RIFit             # RIFit for the hierarchical arm ONLY — §12.6
    dropped: tuple[str, ...]                     # design columns dropped as constant — §5.3


@dataclass(frozen=True)
class Support:
    """[§14a]'s support check: the treated-arm box, who is outside it, and the baseline table.

    [§14a]'s framing is worth keeping beside the code: *"Standardisation avoids infinite weights, not
    extrapolation."* There is no weight to blow up here, so nothing in the arithmetic complains when
    a prediction is made far outside the region where any comparison was observed. This record is the
    only thing that makes it visible.
    """

    box: dict[str, tuple[float, float]]          # covariate -> (treated min, treated max) — §10.1
    inside: pd.Series                            # boolean, TOTAL on the population's index — §10.1
    outside_by_centre: dict[str, int]            # centre -> records outside the box
    never_ivt: tuple[str, ...]                   # the complement of treating_centres — §10.2
    n_never_ivt: int
    n_never_ivt_outside: int                     # [§14a]'s "the proportion … falling outside it"
    baseline: tuple[balance.CovariateBalance, ...]   # treated vs never-IVT centre — §10.3


@dataclass(frozen=True)
class Hierarchical:
    """[§14a]'s sensitivity 2: one common treatment effect, a random centre intercept.

    **`fit` IS `standardisation.fit` — the same object, not an equal one** — and §20.12 asserts `is`
    rather than `==`. Without that the `RIFit` would be reachable by two paths with nothing requiring
    them to agree, which is the duplication this record's shape exists to avoid.

    **`intercepts` and `posterior_sd` are the arm's honest limitation as a measurement rather than a
    caveat** (§12.6). USZ's intercept is the least precisely estimated of the four and it is the one
    the whole arm rests on, because it is the centre whose IVT counterfactual is entirely borrowed —
    12 patients, none of whom received the treatment. That is [§14a]'s *"a centre variance from four
    clusters is fragile"* made specific. The intercept VALUES are centre-level outcome contrasts and
    belong in the gitignored log, never in a version-controlled document (§4.3).
    """

    standardisation: Standardisation             # from the b_hat-conditioned prediction — §12.6
    sigma: float                                 # sigma_hat, the between-centre SD — §12.1
    at_floor: bool                               # sigma_hat reached POLR_RI_SIGMA_FLOOR — §12.5
    intercepts: dict[str, float]                 # centre -> b_hat, the conditional mode — §12.6
    posterior_sd: dict[str, float]               # centre -> the curvature at the mode — §12.6
    nodes: int                                   # POLR_RI_NODES, recorded because it changes an answer
    fit: model.RIFit


@dataclass(frozen=True)
class GComputation:
    """What a g-computation COMPUTES, and nothing about what it means. §7 / Stage 13 §6.1.

    Two averaged distributions keyed by whatever `gcompute` was told to call them, the cumulative
    tables, RD_k and the two namings. NO `beta`, NO `measure`, NO population label: those belong to
    the record that interprets this one — `Standardisation` for [§14a], `Policy` for [§14b] — and a
    record that carried both would be a record one of its two callers has to lie in.
    """

    n_average: int                                  # over.sum()
    n_fit: int                                      # len(X)
    distribution: dict[object, dict[int, float]]    # key -> mRS level -> P, over MRS_LEVELS
    cumulative: dict[object, dict[int, float]]      # key -> threshold -> P(Y <= k)
    rd: dict[int, float]                            # cumulative[keys[0]] - cumulative[keys[1]]
    mrs_0_2: float                                  # == rd[2], one computation
    mortality: float                                # from the distributions; asserted == -rd[5]


# --- the population [§14a, §11, §4] ---------------------------------------------------------------

def record_removal(audit: Audit, step: str, removed: pd.DataFrame, detail: str,
                   table: tuple[tuple[str, ...], ...]) -> None:
    """The one way a §14 population's row removals reach the log. Names every record removed, or raises.

    `cohort._record_removal`'s discipline, made public here as the shared rule with two callers —
    [§14a]'s `population` and Stage 13's `policy.population` [Stage 13 §4.2]. What is shared is the
    ASSERTION — an entry must name exactly as many patients as it says it removed — and it is
    stronger than `data.py`'s kind-keyed rule, which asks only that *some* case be named.

    A population that reported 22 removals and named 3 would satisfy `_MUST_NAME_CASES` and would be
    a population nobody could reconstruct from the log, which matters more here than at Stage 5: this
    population is not `cohort.build`'s and no other entry in the log describes it.
    """
    named = removed.loc[removed["case_id"].notna(), "case_id"]
    case_ids = tuple(sorted(set(named)))
    if len(case_ids) != len(removed):
        raise C.SchemaError(
            f"{step} removed {len(removed)} record(s) and names {len(case_ids)}. A population that "
            "cannot name what it removed makes the analysis unreconstructable from the log.")
    audit.record("cohort", step, len(removed), detail, case_ids=case_ids, table=table)


def _population_table(source: pd.DataFrame, retained: pd.DataFrame,
                      pop: pd.DataFrame) -> tuple[tuple[str, ...], ...]:
    """The per-centre x arm x eligibility ledger, over CENTER_ORDER and never over the frame.

    Ranges over the DECLARED centres for `data.absence_by_column`'s reason: from the data, a run in
    which one centre contributes no rows renders a narrower table that still reconciles, and the
    missing centre is the information. **That is not hypothetical for this stage** — the never-IVT
    centre is the one whose row a reader checks first.
    """
    header = ("centre", "classified", "eligible", "in [§14a] population",
              "EVT alone", "bridging")
    rows = []
    for centre in C.CENTER_ORDER:
        at = pop[pop["center"] == centre]
        rows.append((
            centre,
            str(int((source["center"] == centre).sum())),
            str(int((retained["center"] == centre).sum())),
            str(len(at)),
            str(int((at[C.TREATMENT] == min(C.TREATMENT_LABELS)).sum())),
            str(int((at[C.TREATMENT] == max(C.TREATMENT_LABELS)).sum())),
        ))
    rows.append((
        "all", str(len(source)), str(len(retained)), str(len(pop)),
        str(int((pop[C.TREATMENT] == min(C.TREATMENT_LABELS)).sum())),
        str(int((pop[C.TREATMENT] == max(C.TREATMENT_LABELS)).sum()))))
    return (header, *rows)


def population(df: pd.DataFrame, audit: Audit) -> pd.DataFrame:
    """The [§14a] estimation population: eligible, covariate-complete, outcome-observed. 126 -> 104.

    **Restriction 1 is NOT applied and `cohort.build` is not called**: [§14a]'s population is the one
    restriction 1 removes a centre from [Stage 5 §9]. Restriction 2 — eligibility — IS applied, and
    through `eligibility.retained`, which is the one place the difference between `== eligible` and
    `!= ineligible` is decided.

    **The [§11] mask is spelled here and NOT taken from `outcome.estimation_population`.** That
    helper's first conjunct is a `Propensity`'s covariate-completeness over `PS_COVARIATES`, and this
    stage has no `Propensity` and a different covariate list — calling it would require constructing
    one this stage is forbidden to fit. The two masks are the same SHAPE and different CONTENT, and a
    shared helper parameterised over both would be one function whose two callers agree on nothing
    but the `&` (§4.1).

    **`derive.derive_cohort` is not called, and must not be.** It computes `core_above_median` on the
    frame it is given, and Stage 5 §4.5 established that the median is the COHORT's — 5.0 mL over 93
    records against 6.0 mL over 126. That is a [§13] subgroup variable, Stage 11's and not this
    stage's, and [§14a] names no subgroup. `derive.py` is not imported and §20.1 asserts the name
    appears nowhere in this module.

    T6, T7. Records entries 1 and 2 (§15).
    """
    if C.ELIGIBILITY not in df.columns:
        raise C.SchemaError(
            f"T6  the frame carries no {C.ELIGIBILITY!r} column, so `eligibility.classify` has not "
            "run on it. [§14a]'s population is the eligible one and restriction 2 is applied here "
            "rather than by `cohort.build`; a caller asking for it before anyone has been classified "
            "has a sequencing bug, and a predicate that answered anyway would answer about a "
            "different population.")

    retained = df[eligibility.retained(df)]
    complete = model.complete_cases(retained, C.STANDARDISATION_COVARIATES)
    observed = retained[C.PRIMARY_OUTCOME].notna()
    pop = retained[complete & observed]

    if len(pop) == 0 or any(int((pop[C.TREATMENT] == a).sum()) == 0
                            for a in C.TREATMENT_LABELS):
        counts = {C.TREATMENT_LABELS[a]: int((pop[C.TREATMENT] == a).sum())
                  for a in sorted(C.TREATMENT_LABELS)}
        raise C.SchemaError(
            f"T7  the [§14a] population holds {len(pop)} record(s), by arm {counts}. The "
            "standardisation duplicates every patient under both regimes and averages the two "
            "distributions, so an empty arm makes the treatment coefficient unidentified while "
            "every average remains computable — the arithmetic returns a risk difference for a "
            "contrast the data cannot support. The emptiness is the POPULATION's, after "
            "`eligibility.retained` and the [§11] masks, and not the input frame's.")

    # ENTRY 1. The ledger, and the two [§11] causes are named separately because they fall on
    # different centres and that is a fact about this workbook worth a reader's attention: the
    # outcome-missing records are BOTH at the never-IVT centre, so [§14a]'s target population loses 2
    # of the 14 patients its whole reason for existing is to include (§4.1).
    removed = retained[~(complete & observed)]
    ineligible = df[~eligibility.retained(df)]
    detail = "\n".join([
        f"[§14a]'s population: [§3] restriction 2 applied, restriction 1 NOT — the never-IVT "
        f"centre stays, because it is the population [§14a] exists to include [Stage 5 §9].",
        f"{len(df)} classified -> {len(retained)} eligible -> {len(pop)} in the [§14a] population.",
        f"restriction 2, ineligible [§3]: {len(ineligible)} removed.",
        f"[§11] covariate-incomplete on STANDARDISATION_COVARIATES: "
        f"{int((~complete).sum())} removed.",
        f"[§11] primary outcome not observed: "
        f"{int((complete & ~observed).sum())} removed, by centre "
        f"{ {c: int(((complete & ~observed) & (retained['center'] == c)).sum()) for c in C.CENTER_ORDER} }.",
        "This population is NOT `cohort.build`'s and its step name is disjoint from that module's, "
        "because both remove rows from the same classified frame under kind='cohort' and "
        "`Audit.entry` is first-match [Stage 11 §4.2, Stage 12 §4.2].",
    ])
    record_removal(audit, "standardisation_population",
                    pd.concat([ineligible, removed]), detail,
                    _population_table(df, retained, pop))

    # ENTRY 2.
    absence_by_column(
        pop, audit, list(C.STANDARDISATION_COVARIATES) + [C.PRIMARY_OUTCOME],
        "absence_by_standardisation_column",
        "Absence over the [§14a] population, which is complete on every covariate by construction: "
        "this table exists so that the complete-case denominator [§11] requires is READ rather than "
        "inferred from the ledger above.")
    return pop


# --- the guards [§14] -----------------------------------------------------------------------------

def _assert_exposure_survived(X: pd.DataFrame, dropped: tuple[str, ...]) -> None:
    """T1 — the [§14a] design still carries the treatment column. `FitError`, never `SchemaError`.

    `model.design` drops constant columns and returns their names; it does not raise on them. A
    design that lost the treatment column **fits happily and returns a result in which `beta` does
    not exist**, and the "counterfactual" prediction is then made under a model with no exposure in
    it. Nothing else notices.

    **This differs from Stage 8's G6 in one way that matters**: G6 also asserts nothing else was
    dropped, because Stage 8's design has one column. Here nine covariates may legitimately lose a
    level in a replicate — `onset_type_wake_up` going constant is a real event on a resample — and
    dropping one is a change in the MODEL, not a failure. So a dropped covariate is recorded in
    `Standardisation.dropped` and printed in the fit entry, and only the exposure raises.

    `FitError` and not `SchemaError`, because [§10] must drop and count such a replicate rather than
    crash: this is Stage 9 §12.3's S8 reclassification applied at the point it is created rather than
    after 1.0% of replicates hit it. Measured over 2000 replicates of the workbook, `dropped` is
    empty in 100.00% of them, so T1 never fires here — recorded rather than treated as a reason to
    skip it, on Stage 4's E1-E5 rule and for the workbook that has not arrived yet.
    """
    if C.TREATMENT not in X.columns:
        raise model.FitError(
            f"T1  the [§14a] design lost the treatment column; {len(dropped)} column(s) were dropped "
            f"as constant: {', '.join(dropped) or 'none'}. A design with no exposure FITS, and "
            "returns a standardisation in which both regimes are the same regime — every risk "
            "difference is zero and nothing raises. [§10] drops and counts the replicate.")


def _assert_beta_reportable(value: float, key: str) -> None:
    """T12 — a fitted treatment coefficient below `POLR_MAX_ABS_BETA`. `FitError`.

    **The bound protects a quantity this stage reports as a MODEL PARAMETER and nothing else**, and
    that is the whole difference from Stage 8 (§9.1). Stage 8's reported quantity is `exp(beta)`, and
    a separated fit returns 6.5e15 — a number no manuscript could print. Here the reported quantities
    are averaged probabilities, bounded in [0, 1] BY CONSTRUCTION: measured over 2000 replicates,
    every standardised probability is in [0, 1] and every risk difference in [-1, 1] no matter how
    degenerate the fit is, because `expit` of anything is in [0, 1].

    So the guard is at KEY granularity: it costs the caller's `beta` key and keeps the standardised
    keys, which are unaffected and correct. Coupling them would discard 22 valid keys to suppress one
    invalid one, which is a bias in the surviving-draw set rather than a safeguard.

    **THE KEY COST IS THE CALLER'S DECISION AND NOT THIS HELPER'S.** It knows a coefficient and a
    bound; it does not know which keys a caller will lose. A helper that did would be one function
    with two call sites agreeing on nothing but the comparison — §4.1's objection to the shared [§11]
    mask, applied here. The two call sites are the pooled fit (costing `beta`) and `polr_ri` (costing
    `hier.beta`), and each names its own key.

    **`outcome._assert_reportable` is NOT reused**, and §24 records why: it bounds the WORST
    coefficient in the fit, which is right for Stage 8's one-column design where the worst
    coefficient IS the treatment coefficient, and wrong here. §9.2 requires `gamma` unbounded — a
    nuisance coefficient may grow legitimately when a covariate is nearly collinear in a resample —
    and the largest |gamma| measured over 2000 replicates is 3.21 against a |beta_treatment| of 2.00.
    """
    if abs(value) >= C.POLR_MAX_ABS_BETA:
        raise model.FitError(
            f"T12  the fitted treatment coefficient for {key!r} is {value:.6g}, reaching "
            f"POLR_MAX_ABS_BETA = {C.POLR_MAX_ABS_BETA:g}. A separated proportional-odds fit "
            "CONVERGES AND RETURNS [Stage 8 §6], so this bound is the only thing that turns a "
            "degenerate fit into a countable failure. It costs this key alone: the standardised "
            "quantities are averaged probabilities and stay in range whatever beta does (§9.1).")


def _assert_distribution(distribution: dict[int, dict[int, float]],
                         cumulative: dict[int, dict[int, float]], label: str) -> None:
    """T10 — each arm's standardised distribution sums to 1 and its cumulative sequence is monotone.

    [§14a]'s Accept-when conditions expressed as code. **Unreachable on this workbook** — measured
    worst `|sum - 1|` of 6.7e-16 across every arm of every one of 2000 replicates, against a bound of
    1e-12 — which is recorded rather than treated as a reason to omit it (§20.5 is the argument for
    asserting a property that cannot currently fail).

    `SchemaError`, so it is NOT caught inside a replicate and kills the whole bootstrap. That is
    deliberate: a distribution that does not sum to 1 is a contract break in the arithmetic, not a
    sparse resample, and eighteen minutes of draws built on it are worth less than the crash. The
    measured margin — 6.7e-16 against 1e-12 — is what makes that safe to specify.
    """
    for arm, levels in distribution.items():
        total = sum(levels.values())
        if not np.isfinite(total) or abs(total - 1.0) > _DISTRIBUTION_TOL:
            raise C.SchemaError(
                f"T10  the standardised distribution for arm {arm} of {label!r} sums to {total!r}, "
                f"off 1 by more than {_DISTRIBUTION_TOL:g}. `ordinal_probabilities` brackets its "
                "cumulative sequence with EXACT 0.0 and EXACT 1.0 and differences it, and the "
                "re-expansion adds structural zeros — so a row sums to 1 by construction and an "
                "average of such rows does too. This is arithmetic, not a property of the cohort.")
        sequence = [cumulative[arm][k] for k in C.MRS_THRESHOLDS]
        if any(b < a - _DISTRIBUTION_TOL for a, b in zip(sequence, sequence[1:])):
            raise C.SchemaError(
                f"T10  the cumulative sequence for arm {arm} of {label!r} is not non-decreasing: "
                f"{sequence}. `alpha` is ascending, so `P(Y <= k)` is monotone per row and an "
                "average of monotone sequences is monotone — including on a replicate that lost a "
                "level, because a structural zero leaves the sequence monotone (§6.3).")


def _assert_identities(rd: dict[int, float], mrs_0_2: float, mortality: float,
                       label: str) -> None:
    """T11 — [§14a]'s two namings agree with the `RD_k` they are named for. §7.3, §20.6.

    Two of the four quantities [§14a] lists as outputs are not new numbers: the standardised mRS 0-2
    risk difference IS `RD_2` and the standardised mortality difference IS `-RD_5`. That is not a
    defect in [§14a] — they are the quantities a clinical reader wants named — but it IS a defect in
    any implementation that computes them twice, because two computations of one number are two
    things that can disagree after an edit.

    **The two tolerances differ and the asymmetry is §7.3's rather than arbitrary.** `mrs_0_2` is
    `rd[2]`, the same cumulative difference read under [§14a]'s name for it, so the comparison is
    EXACT. `mortality` is genuinely a second computation — from `P(Y = 6)` in each arm rather than
    from `1 - P(Y <= 5)` — so it agrees to floating point and not to the bit; measured 8.3e-17
    against a bound of 1e-15.
    """
    if mrs_0_2 != rd[2]:
        raise C.SchemaError(
            f"T11  {label!r}: mrs_0_2 is {mrs_0_2!r} against rd[2] = {rd[2]!r}. [§14a]'s mRS 0-2 "
            "risk difference IS the cumulative RD at threshold 2, read under its clinical name, so "
            "the two are ONE computation and the comparison is exact. A difference here means a "
            "second computation was introduced (§7.3).")
    if not np.isfinite(mortality) or abs(mortality + rd[5]) > _MORTALITY_TOL:
        raise C.SchemaError(
            f"T11  {label!r}: mortality is {mortality!r} against -rd[5] = {-rd[5]!r}, differing by "
            f"more than {_MORTALITY_TOL:g}. `P(Y = 6) = 1 - P(Y <= 5)` in both arms and the two 1s "
            "cancel, so the identity is one line of algebra. `mortality` is computed from the "
            "DISTRIBUTIONS and this is what holds the two computations together (§20.6).")


# --- the standardisation [§7] ---------------------------------------------------------------------

def _regime_design(X: pd.DataFrame, arm: np.ndarray | float) -> pd.DataFrame:
    """`X` with the treatment COLUMN overwritten by `arm` — a per-row vector, or a scalar that
    broadcasts. Never a new design, and §7.1 is why.

    [§14a] says *"duplicate every patient with `A = 1` and `A = 0`"*. The obvious implementation —
    copy the frame, set the treatment column, re-run `model.design` — **is wrong and it fails
    silently**: `design` drops constant columns, so a frame in which every patient has `A = 1` yields
    a design with NO TREATMENT COLUMN AT ALL, `polr` fits the remaining nine, and the
    "counterfactual" prediction is made under a model with no exposure in it. Nothing raises, and T1
    would not catch it either, because the counterfactual design is not the one fitted.

    **This is the pattern `outcome._counterfactuals` already uses** for the Stage 9 binary nuisance
    models, and Stage 9 chose it for the same reason: the by-name column check passes because the
    columns ARE the fitted ones, and an overwrite cannot change which columns exist.

    The arm codes come from `C.TREATMENT_LABELS` exactly as `outcome.py:118-119` takes them, and are
    never written as 1 and 0 in the arithmetic.
    """
    regime = np.asarray(arm, dtype=float)
    if regime.ndim == 1 and len(regime) != len(X):
        raise C.SchemaError(
            f"a regime vector of length {len(regime)} against a design of {len(X)} row(s). The "
            "regime is written INTO the fitted design row by row [Stage 13 §6.1]; a length mismatch "
            "would broadcast or misalign silently.")
    Xa = X.copy()
    Xa[C.TREATMENT] = regime
    return Xa


def _arm_probabilities(fit: model.PolrFit | model.RIFit, Xa: pd.DataFrame,
                       groups: pd.Series | None) -> np.ndarray:
    """P(Y = c | x) per row over `fit.categories`, under whichever of the two fitters produced `fit`.

    **A RANDOM INTERCEPT IS A PER-CENTRE SHIFT OF THE CUTPOINTS, and that is the whole of the
    hierarchical branch** — no offset parameter, no second prediction function::

        expit(alpha_k + x'beta + b_c)  ==  expit((alpha_k + b_c) + x'beta)

    So each centre's rows are predicted through `model.ordinal_probabilities` UNCHANGED, against a
    `PolrFit` whose `alpha` is this centre's shifted cutpoints. [§14a] prescribes the conditional
    standardisation rather than the marginal one, in a sentence that is easy to read past — *"Patients
    at a never-IVT centre inform that centre's intercept through their direct-EVT outcomes while the
    IVT effect is borrowed"* — and a patient can only inform their own centre's intercept if the
    prediction uses it (§12.6). Under the marginal reading USZ's 12 patients would inform `sigma` and
    nothing else, and that sentence would be false.

    The marginal alternative is a DIFFERENT ESTIMAND and is recorded in §22 as a [§14] question
    rather than implemented.
    """
    if isinstance(fit, model.PolrFit):
        return model.ordinal_probabilities(fit, Xa)

    if groups is None:
        raise C.SchemaError(
            "a random-intercept standardisation needs the grouping vector, because the prediction is "
            "conditional on each centre's own b_hat [§14a, §12.6]. Passing None would silently "
            "standardise at b = 0, which is neither of [§14a]'s two readings.")
    labels = np.asarray([str(g) for g in np.asarray(groups)])
    out = np.full((len(Xa), len(fit.categories)), np.nan)
    for index, label in enumerate(fit.groups):
        rows = labels == label
        if not rows.any():
            continue
        shifted = model.PolrFit(
            beta=fit.beta, alpha=fit.alpha + fit.b[index], categories=fit.categories,
            columns=fit.columns, iterations=fit.iterations, converged_on=fit.converged_on,
            first_step_norm=fit.first_step_norm, rescales=fit.rescales, halvings=fit.halvings)
        out[rows] = model.ordinal_probabilities(shifted, Xa[rows])
    if not np.all(np.isfinite(out)):
        unseen = sorted(set(labels.tolist()) - set(fit.groups))
        raise C.SchemaError(
            f"the standardisation frame carries {len(unseen)} group(s) the fit has no intercept for: "
            f"{unseen}. A patient at a centre the fit never saw has no `b_hat`, so their conditional "
            "prediction does not exist — and averaging over a nan would make every reported "
            "probability nan rather than raising here.")
    return out


def _expanded(probabilities: np.ndarray, categories: tuple[int, ...]) -> np.ndarray:
    """`probabilities` re-expressed on `C.MRS_LEVELS`; an unoccupied level becomes a STRUCTURAL ZERO.

    `polr` collapses the response to the categories carrying positive weight, so
    `ordinal_probabilities` returns `len(fit.categories)` columns — six where the declared level set
    has seven, in 4.15% of this stage's replicates, always mRS 5 (§6.2). **A six-column distribution
    cannot be averaged with a seven-column one**, and a naive implementation would either raise deep
    inside numpy or — worse — broadcast and silently misalign every level above the missing one.

    The rule is prescribed rather than convenient and the argument is Stage 7 §5.3's about `nan`: a
    level no record occupies is a level for which the fit provides no cutpoint, and the honest value
    is the one the model implies — zero mass. Interpolating a cutpoint would be inventing a parameter;
    returning a shorter vector would propagate the misalignment into the average.

    Three consequences, each specified (§6.3): the distribution still sums to 1, because adding zeros
    changes no sum; `RD_k` is still defined at every threshold, because a structural zero leaves the
    cumulative sequence monotone; and **`RD_4` and `RD_5` coincide exactly in a replicate that lost
    level 5**, which is a real narrowing of `RD_5`'s sampling distribution, is not a bug, and is
    counted and reported beside the interval.
    """
    position = {level: j for j, level in enumerate(categories)}
    out = np.zeros((probabilities.shape[0], len(C.MRS_LEVELS)))
    for column, level in enumerate(C.MRS_LEVELS):
        if level in position:
            out[:, column] = probabilities[:, position[level]]
    return out


def gcompute(X: pd.DataFrame, fit: model.PolrFit | model.RIFit, over: np.ndarray,
             active: np.ndarray, comparator: np.ndarray, *, label: str,
             groups: pd.Series | None = None,
             keys: tuple[object, object] = (max(C.TREATMENT_LABELS), min(C.TREATMENT_LABELS)),
             ) -> GComputation:
    """The g-computation: overwrite the treatment column of the FITTED design with `active`, then
    with `comparator`, predict, re-expand onto MRS_LEVELS, average over `over`. §7, Stage 13 §6.1.

        P_hat(Y = j | regime)  =  (1/N) * SUM_i  P(Y = j | A_i(regime), X_i)

    over the `N = over.sum()` records of the averaging population, UNWEIGHTED. `over` restricts the
    averaging set and never the fit, which is what makes §11's sensitivity one change rather than two.

    **Patients at never-IVT centres contribute their covariates and their share of the target
    population**, which is [§14a]'s own sentence and the whole point: their IVT counterfactual is
    predicted from relationships learned where IVT was observed. Nothing in the code distinguishes
    them — `center` is not in the model — and that is what `support` exists to make visible rather
    than to fix.

    The orientation is Stage 8 §7's, unchanged: `arm[1] - arm[0]`, so a positive `RD_k` favours
    bridging on the functional scale and **a positive `mortality` is WORSE**. That sign flip between
    the two reported quantities is exactly why `mortality` is named rather than left as `-RD_5` for a
    reader to negate.

    T10 and T11 run on every record this function builds — the point estimate AND every replicate —
    because a self-check that runs only on the point estimate is a self-check the bootstrap does not
    have.

    **`active` and `comparator` are per-row arrays of length `len(X)`; a length mismatch or a scalar
    raises `SchemaError`.** `_standardise` passes the two arm codes broadcast to vectors; Stage 13
    passes `1[eligible]` and zeros, which is the whole of what makes a REGIME differ from an ARM.
    `keys` names the two distributions — Stage 12's arm codes by default, Stage 13's regime literals —
    and `keys[0]` is the minuend of every RD_k. The record carries no `beta`, `measure` or population
    label: those are the interpreting record's (`GComputation`).
    """
    for name, regime in (("active", active), ("comparator", comparator)):
        shape = np.shape(regime)
        if len(shape) != 1 or shape[0] != len(X):
            raise C.SchemaError(
                f"`{name}` has shape {shape} against a design of {len(X)} row(s). A regime is a "
                "per-row vector of length len(X) and the vector form is the only form: a scalar "
                "here would hide which rows a caller meant to assign [Stage 13 §6.1].")

    distribution: dict[object, dict[int, float]] = {}
    cumulative: dict[object, dict[int, float]] = {}
    for key, regime in zip(keys, (active, comparator)):
        per_row = _expanded(
            _arm_probabilities(fit, _regime_design(X, regime), groups), fit.categories)
        averaged = per_row[over].mean(axis=0)                     # UNWEIGHTED — [§14a] weights nobody
        distribution[key] = {level: float(p) for level, p in zip(C.MRS_LEVELS, averaged)}
        running = np.cumsum(averaged)
        cumulative[key] = {k: float(running[j]) for j, k in enumerate(C.MRS_THRESHOLDS)}

    first, second = keys
    rd = {k: cumulative[first][k] - cumulative[second][k] for k in C.MRS_THRESHOLDS}
    _assert_distribution(distribution, cumulative, label)         # T10

    # [§14a]'s TWO NAMINGS. `mrs_0_2` IS `rd[2]` — one computation read under its clinical name, not a
    # second one — and `mortality` is computed from the DISTRIBUTIONS, which is what makes it a
    # genuine second computation and what §20.6's monkeypatch test protects. Deriving it as `-rd[5]`
    # would pass the identity trivially and would be the "simplification" that hides a real defect.
    top = C.MRS_LEVELS[-1]
    mrs_0_2 = rd[2]
    mortality = distribution[first][top] - distribution[second][top]
    _assert_identities(rd, mrs_0_2, mortality, label)             # T11

    return GComputation(
        n_average=int(over.sum()), n_fit=len(X), distribution=distribution,
        cumulative=cumulative, rd=rd, mrs_0_2=mrs_0_2, mortality=mortality)


def _standardise(X: pd.DataFrame, fit: model.PolrFit | model.RIFit, over: np.ndarray,
                 label: str, dropped: tuple[str, ...],
                 groups: pd.Series | None = None) -> Standardisation:
    """One [§14a] standardisation: `gcompute` with the two ARM codes broadcast to vectors, wrapped
    with `beta`, `_CONDITIONAL` and the [§14a] population label. Signature and callers unchanged
    across the Stage 13 extraction, which §16.11 of that spec asserts by two digests."""
    treated, control = max(C.TREATMENT_LABELS), min(C.TREATMENT_LABELS)
    g = gcompute(X, fit, over, np.full(len(X), float(treated)), np.full(len(X), float(control)),
                 label=label, groups=groups, keys=(treated, control))

    beta = float(fit.beta[fit.columns.index(C.TREATMENT)])
    return Standardisation(
        population=label, n_average=g.n_average, n_fit=g.n_fit,
        distribution=g.distribution, cumulative=g.cumulative, rd=g.rd,
        mrs_0_2=g.mrs_0_2, mortality=g.mortality,
        conditional_log_odds=beta, conditional_odds_ratio=float(np.exp(beta)),
        measure=_CONDITIONAL, fit=fit, dropped=dropped)


def _over_mask(over: pd.Series | None, pop: pd.DataFrame) -> np.ndarray:
    """T8 — `over` is boolean and indexed like the population, or it raises. §20.14.

    Three cases, one point about alignment. `over` is applied inside every replicate to a RESAMPLED
    frame, and that frame is safe to index by label only because `bootstrap.resample` RESETS THE
    INDEX before renaming (`bootstrap.py:315-320`, which gives Stage 12 as the reason: *"a frame with
    a duplicated index is a frame on which a future label-based `.loc` is wrong in a way that returns
    a number"*). T8 is the guard for the case that reset does not cover, and the reset is why T8 does
    not have to cover the duplicated-label case at all.
    """
    if over is None:
        return np.ones(len(pop), dtype=bool)
    # BOTH BOOLEAN DTYPES ARE ACCEPTED AND THE MISSING-VALUE CHECK IS THE REAL GUARD. A comparison
    # against a nullable column returns pandas' `boolean` and not numpy's `bool` — `_inside` produces
    # exactly that, because the [§6] covariates are Float64 — so a dtype test against `bool` alone
    # rejects this stage's own mask. What actually has to be caught is a NON-boolean mask, which
    # indexes by POSITION rather than by membership, and a boolean mask carrying `pd.NA`, which
    # `astype(bool)` silently turns into True.
    if not pd.api.types.is_bool_dtype(over.dtype):
        raise C.SchemaError(
            f"T8  `over=` has dtype {over.dtype!r} and must be boolean. A float or Int64 mask indexes "
            "by POSITION where a boolean one indexes by membership, so a non-boolean mask silently "
            "averages over a different set of patients and reports the size it averaged over as if "
            "it were the population [§11].")
    if bool(over.isna().any()):
        raise C.SchemaError(
            f"T8  `over=` carries {int(over.isna().sum())} missing value(s). `astype(bool)` maps "
            "`pd.NA` to True, so a mask with a gap in it silently INCLUDES the patients it could not "
            "decide about — the opposite of what a support restriction means [§11].")
    if not over.index.equals(pop.index):
        raise C.SchemaError(
            f"T8  `over=` carries {len(over)} label(s) against the population's {len(pop)}, and its "
            "index does not match — either in content or in ORDER. `over` restricts the averaging "
            "population, so a misaligned mask averages the right number of the wrong patients. The "
            "index is compared as a SEQUENCE and not as a set, because a shuffled index reorders the "
            "mask against the design's rows.")
    return over.to_numpy(dtype=bool)


# --- the [§14a] model and its standardisation [§5, §7, §16] ---------------------------------------

_DESIGN_COLUMNS: Final[tuple[str, ...]] = (C.TREATMENT,) + C.STANDARDISATION_COVARIATES


def _fit_table(fit: model.PolrFit, dropped: tuple[str, ...]) -> tuple[tuple[str, ...], ...]:
    """The fit entry's coefficient table: every design column, its coefficient, and the guard margin.

    The magnitudes are reported against `POLR_MAX_ABS_BETA` because §9.2's diagnostic is that the
    bound is NOWHERE NEAR on this design — largest coefficient 3.21 against 14.0 over 2000 replicates,
    against Stage 8 §6.3's measured empty band running from 8.79 up — and a diagnostic nobody can read
    is not one. `gamma` is deliberately unbounded (§9.2): a nuisance coefficient may grow legitimately
    when a covariate is nearly collinear in a resample, and bounding it would drop replicates for a
    reason the reported quantity does not care about.
    """
    header = ("design column", "coefficient", "|coef|", "role")
    rows = [(name, f"{value:+.6f}", f"{abs(value):.6f}",
             "EXPOSURE — the guarded one [T12]" if name == C.TREATMENT else "gamma, unbounded [§9.2]")
            for name, value in zip(fit.columns, fit.beta)]
    for name in dropped:
        rows.append((name, "—", "—", "DROPPED as constant [§5.3]"))
    return (header, *rows)


def _distribution_table(std: Standardisation) -> tuple[tuple[str, ...], ...]:
    """The two standardised distributions, the cumulative table and `RD_k`, level by level.

    One table and not three, because a reader compares the arms level by level and the cumulative
    column is what `RD_k` is a difference OF. The last two rows carry [§14a]'s two namings **as
    computed differences**, so the log shows they were checked rather than asserting it in prose
    (§15, entry 4).
    """
    treated, control = max(C.TREATMENT_LABELS), min(C.TREATMENT_LABELS)
    header = ("mRS", f"P | {C.TREATMENT_LABELS[treated]}", f"P | {C.TREATMENT_LABELS[control]}",
              "P(Y<=k) | bridging", "P(Y<=k) | EVT alone", "RD_k")
    rows = []
    for level in C.MRS_LEVELS:
        cumulative = (
            (f"{std.cumulative[treated][level]:.6f}", f"{std.cumulative[control][level]:.6f}",
             f"{std.rd[level]:+.6f}") if level in std.rd else ("1.000000", "1.000000", "—"))
        rows.append((str(level), f"{std.distribution[treated][level]:.6f}",
                     f"{std.distribution[control][level]:.6f}", *cumulative))
    rows.append(("mRS 0-2 [§14a]", "", "", "", "", f"{std.mrs_0_2:+.6f}"))
    rows.append(("mortality [§14a]", "", "", "", "", f"{std.mortality:+.6f}"))
    rows.append(("|mrs_0_2 - RD_2|", "", "", "", "",
                 f"{abs(std.mrs_0_2 - std.rd[2]):.3e}"))
    rows.append(("|mortality + RD_5|", "", "", "", "",
                 f"{abs(std.mortality + std.rd[5]):.3e}"))
    return (header, *rows)


def all_centre(pop: pd.DataFrame, audit: Audit,
               over: pd.Series | None = None) -> Standardisation:
    """One [§14a] standardisation. `over` restricts the AVERAGING population, never the fit (§11).

        logit P(Y <= k | A, X)  =  alpha_k  +  beta*A  +  gamma'X ,    k = 0 … 5

    with `X = C.STANDARDISATION_COVARIATES`, which `config.py` declares as `PS_COVARIATES` minus
    `center`. The list is not written out here and must not be written out anywhere: Stage 1's rule is
    that a covariate list lives in `config.py` and that none may appear as a second literal.

    **`center` is omitted because the identifying assumption is exactly that potential outcomes do not
    depend on centre given `X`**, and the consequence is worth spelling. Including `center` would make
    the model contradict the assumption the analysis runs on, AND it would make the standardisation
    uncomputable in the direction that matters: a `center_USZ` coefficient is estimated from a stratum
    with no treated patient, so the IVT counterfactual for a USZ patient would be an extrapolation
    along a coefficient no USZ data can identify. Omitting centre is what makes the transport EXPLICIT
    rather than hidden inside a dummy. [§15] permits this *"inside §14 only"* and `test_config.py`
    holds `STANDARDISATION_COVARIATES` to that permission statically.

    **The fit is UNWEIGHTED and `w` is not passed.** `polr`'s docstring is explicit that `w=None` is
    the unweighted fit and *"is not a synonym for ones: a caller who passes nothing has said something
    different from a caller who passes ones"*. [§14a] weights nobody — there is no propensity model to
    weight by — so this passes nothing, and §20.4 asserts the CALL SITE rather than the arithmetic.

    **Contraindication status is not a covariate**, [§14a] having already restricted to eligible
    patients: the `eligibility` column is the restriction and never a regressor.

    T1, T8, T10, T11. Records entries 3 and 4, or entry 7 when `over` is given (§15).
    """
    X, dropped = model.design(pop, _DESIGN_COLUMNS)
    _assert_exposure_survived(X, dropped)                        # T1
    mask = _over_mask(over, pop)                                 # T8
    y = pop[C.PRIMARY_OUTCOME].to_numpy(dtype=float)
    fit = model.polr(X, y)                                       # UNWEIGHTED — no `w`, §5.2
    label = _ALL_ELIGIBLE if over is None else _TREATED_SUPPORT
    std = _standardise(X, fit, mask, label, dropped)              # T10, T11

    if over is None:
        # ENTRY 3.
        audit.record(
            "model", "standardisation_fit", len(pop), "\n".join([
                f"[§14a]'s pooled proportional-odds model: {len(fit.columns)} design column(s), "
                f"{len(dropped)} dropped as constant, {len(fit.alpha)} cutpoint(s) over "
                f"{len(fit.categories)} fitted mRS level(s) {fit.categories}.",
                f"Converged in {fit.iterations} iteration(s) on the {fit.converged_on!r} criterion; "
                f"{fit.rescales} trust-region rescale(s), {fit.halvings} halving(s); first step norm "
                f"{fit.first_step_norm:.6g}.",
                "UNWEIGHTED: [§14a] weights nobody and no propensity model is fitted anywhere in "
                "[§14]. `center` is NOT a covariate — [§15] permits omitting it inside [§14] only, "
                "where extrapolating to centres that never used IVT is the declared purpose.",
                f"The fitted level set is the SAMPLE's: {len(C.MRS_LEVELS) - len(fit.categories)} "
                "declared level(s) unoccupied, each a structural zero in the standardisation [§6.3].",
                f"Largest |coefficient| {float(np.max(np.abs(fit.beta))):.4f} against "
                f"POLR_MAX_ABS_BETA = {C.POLR_MAX_ABS_BETA:g}; largest |beta_treatment| "
                f"{abs(float(fit.beta[fit.columns.index(C.TREATMENT)])):.4f}.",
                "**T12 IS NOT CHECKED HERE AND THAT IS §9.1's, NOT AN OMISSION.** The guard has TWO "
                "call sites and both are in the replicate body: a separated fit converges and "
                "returns [Stage 8 §6], so in the bootstrap the bound is the only thing that turns "
                "one into a countable failure, and at the POINT ESTIMATE this table is what makes it "
                "visible instead. `gamma` is unbounded by design [§9.2] — a nuisance coefficient may "
                "grow legitimately when a covariate is nearly collinear in a resample.",
            ]), table=_fit_table(fit, dropped))
        # ENTRY 4.
        audit.record(
            "model", "standardised_distributions", std.n_average, "\n".join([
                f"[§14a]'s g-computation: every patient duplicated under both regimes on the FITTED "
                f"design, averaged over {std.n_average} record(s) of {std.n_fit} fitted, unweighted.",
                f"Patients at never-IVT centres contribute their covariates and their share of the "
                f"target population — their IVT counterfactual is predicted from relationships "
                f"learned where IVT was observed. `support` is what makes that visible.",
                f"conditional log-odds {std.conditional_log_odds:+.6f}, exp() "
                f"{std.conditional_odds_ratio:.6f} — {std.measure}.",
                "The last two rows of the table are [§14a]'s two namings as COMPUTED DIFFERENCES "
                "against their RD_k, so the log shows the identities were checked [§7.3, T11].",
            ]), table=_distribution_table(std))
    else:
        # ENTRY 7. Sensitivity 1, and it falls between 6 and 8 because `support` produces its `over`.
        audit.record(
            "model", "standardisation_support_restricted", std.n_average, "\n".join([
                f"[§14a] sensitivity 1: the standardisation population restricted to the treated "
                f"support. n_fit {std.n_fit}, n_average {std.n_average}.",
                "THE AVERAGING POPULATION IS RESTRICTED AND THE FIT IS NOT. Refitting on the "
                "in-support subset is a defensible reading and is DECLINED: it changes the target "
                "population AND the fitted model at once, so a difference between the arms could not "
                "be attributed to either. [§14a]'s sentence names the standardisation population, "
                "which is the averaging set, and the one-thing-at-a-time reading is also the literal "
                "one [§11]. Whether a refit-on-support arm should exist is a [§14] amendment [§22].",
                f"conditional log-odds {std.conditional_log_odds:+.6f} — IDENTICAL to the pooled "
                "arm's by construction, which is why there is no `support.beta` estimand key [§3.3].",
            ]), table=_distribution_table(std))
    return std


# --- the support check [§14a, §10] ----------------------------------------------------------------

def _box(pop: pd.DataFrame) -> dict[str, tuple[float, float]]:
    """The per-covariate treated-arm `[min, max]` box over `C.SUPPORT_COVARIATES`. §10.1.

    **Not a percentile trim and not a convex hull.** The box is the WEAKEST possible reading of
    "outside the treated support" — a patient outside it is outside on a SINGLE covariate's observed
    range, which no smoothing assumption can be argued to cover — and it is the only definition under
    which "outside" needs no tuning parameter. A convex hull in seven dimensions on 39 treated
    patients is a set almost every point is outside of, which would make the diagnostic report 100%
    and mean nothing.

    "Continuous" is resolved to "not a declared factor" by `C.SUPPORT_COVARIATES`, which is the WIDER
    of the two readings and is MEASURED to be inert: `sex`, `prestroke_mrs` and `atrial_fib` exclude
    nobody, because the treated arm covers both levels of each binary and reaches `prestroke_mrs = 3`.
    The wider reading is chosen because it cannot be wrong on a workbook where the narrow one is
    right, and §20.10 asserts the inertness rather than claiming it.
    """
    treated = pop[pop[C.TREATMENT] == max(C.TREATMENT_LABELS)]
    return {c: (float(treated[c].min()), float(treated[c].max()))
            for c in C.SUPPORT_COVARIATES}


def _inside(pop: pd.DataFrame, box: dict[str, tuple[float, float]]) -> pd.Series:
    """Boolean, TOTAL on `pop`'s index: True where every boxed covariate is within the treated range.

    **Everyone outside is a control patient, and that is arithmetic rather than a finding** — the box
    is the treated arm's own range, so no treated patient can be outside it. It is stated because a
    reader meeting "11 outside" will otherwise wonder about the split.
    """
    inside = pd.Series(True, index=pop.index)
    for covariate, (low, high) in box.items():
        inside &= (pop[covariate] >= low) & (pop[covariate] <= high)
    # `bool` and not the nullable `boolean` the comparisons produce: the boxed covariates are complete
    # on this population by construction — `population` complete-cases on them — so there is no NA to
    # lose, and returning the numpy dtype keeps `Support.inside` a mask every caller can index with.
    # `_over_mask` is what checks the assumption rather than trusting it.
    return inside.astype(bool)


def _never_ivt(pop: pd.DataFrame) -> tuple[str, ...]:
    """The complement of `cohort.treating_centres`, computed at the one place that needs it. T9.

    `cohort.py` wrote this expression into its own docstring as the thing Stage 12 should use, and
    Stage 5 §4.2's argument is why: a centre's treatment availability is a PROPERTY OF THE DATA, and a
    workbook in which USZ starts administering IVT must change this diagnostic rather than require
    someone to remember a list. `EXPECTED_NEVER_IVT` is an assertion target and is never the operative
    rule — its own comment in `config.py` says exactly that.

    **T9 is a `SchemaError` and not a finding.** A workbook in which the never-IVT centre acquires a
    treated patient, or in which a fifth centre appears, changes what [§14a]'s support check is ABOUT
    — the never-IVT set is the thing the check is defined against — and the pipeline must stop rather
    than report a diagnostic whose subject moved.

    **This is why `support` is NOT called inside a replicate.** A resample of Lugano's 30 records
    misses both its bridging patients with probability 0.13, which would put Lugano in this set and
    kill the bootstrap on a `SchemaError`. The replicate body recomputes the BOX (§11) and never the
    never-IVT set, because the box is a statistic and this is a contract.
    """
    keep = cohort.treating_centres(pop)
    never = tuple(c for c in C.CENTER_ORDER if c not in keep)
    if never != C.EXPECTED_NEVER_IVT:
        raise C.SchemaError(
            f"T9  the centres contributing no bridging patient are {never} against "
            f"EXPECTED_NEVER_IVT = {C.EXPECTED_NEVER_IVT}. [§14a]'s support check is DEFINED against "
            "the never-IVT set — it reports the proportion of those patients outside the treated "
            "support and compares them with the treated in a baseline table — so a change here "
            "changes what the diagnostic is about. The pipeline stops rather than reporting a "
            "diagnostic whose subject moved.")
    return never


def _role(name: str) -> str:
    """What this covariate is TO [§14a]'s specification. §10.3.

    **`balance._role` is NOT reused, and the reason is not access but CORRECTNESS.** It returns
    `"propensity model"` for every `PS_COVARIATES` name, which is right for the specification it was
    written to diagnose and wrong here: **this stage fits no propensity model anywhere** ([§14]'s own
    first sentence), so a row labelled "propensity model" in a [§14a] table would name a model that
    does not exist in this analysis. Its docstring says what a role is — *"what this covariate is TO
    the specification being diagnosed"* — and the specification being diagnosed here is a different
    one. So this is a different rule, not a second copy of the same rule.

    `center` is the GROUPING variable and carries `_GROUPING`: it is in `BALANCE_SET` and therefore in
    the table, and it is not in `STANDARDISATION_COVARIATES` at all.
    """
    if name == C.BOOT_STRATUM:
        return _GROUPING
    if name in C.STANDARDISATION_COVARIATES:
        return "standardisation model [§14a]"
    if name in C.NEGATIVE_CONTROLS:
        return "negative control"
    return "excluded [§6]"


def _baseline(pop: pd.DataFrame, never: tuple[str, ...]) -> tuple[balance.CovariateBalance, ...]:
    """[§14a]'s baseline table: the treated patients anywhere against the never-IVT centre's. §10.3.

    Two DISJOINT groups — every treated patient, and every eligible patient at a never-IVT centre —
    and **the control patients at treating centres are in NEITHER**. They are not what [§14a] asked to
    compare, and their absence is printed in the caption with all three counts so the denominator is
    never inferred.

    **It calls `balance.smd` directly, with UNIT weights**, which is the note Stage 7 §5.5 wrote to
    this stage: *"[§14a]'s support check is a standardised mean difference over a different population
    with unit weights, and a second implementation there would be a second definition of the
    yardstick."* And it calls `balance.levels` for the rows themselves, which is the same note written
    twice — a second implementation of that would be a second definition of WHAT A ROW IS, and its
    failure would be invisible here, because all three of `onset_type`'s declared levels are present
    in the treated arm and a copy ranging over OBSERVED levels would produce a byte-identical table on
    v7 (§10.1, §10.3).

    `balance._pooled_sd` supplies the `sd` COLUMN and nothing else. That column is the denominator
    `smd` divided by, so computing it here a second way would be the second definition of the
    yardstick §10.3 forbids; reading it from the module that owns it is what keeps there being one.

    The rows range over `C.BALANCE_SET` — the [§6] confounders plus the five balance-only covariates —
    because [§9]'s rule that balance is judged against the FULL confounder set is about what a reader
    must be shown, and it does not stop applying because the comparison is not a treatment contrast.

    **One of `smd`'s five documented `nan` routes fires here BY CONSTRUCTION and is reported rather
    than suppressed**: `center = USZ` is constant at 1 in one group and 0 in the other, so the pooled
    SD is exactly zero with the groups differing — Stage 7 §5.3's branch 4 — and the row is `nan`.
    """
    treated = pop[C.TREATMENT] == max(C.TREATMENT_LABELS)
    at_never = pop["center"].isin(never)
    sub = pop[treated | at_never]
    arm = treated[treated | at_never].to_numpy(dtype=float)
    unit = np.ones(len(sub))

    rows: list[balance.CovariateBalance] = []
    for name in C.BALANCE_SET:
        for label, indicator in balance.levels(sub, name):
            values = indicator.to_numpy(dtype=float)
            present = np.isfinite(values)
            x, group, weight = values[present], arm[present], unit[present]
            rows.append(balance.CovariateBalance(
                covariate=label, role=_role(name), n=int(present.sum()),
                sd=balance._pooled_sd(x, group),
                unweighted=balance.smd(x, group, weight),
                weighted=balance.smd(x, group, weight)))
    return tuple(rows)


def _box_table(pop: pd.DataFrame, box: dict[str, tuple[float, float]],
               inside: pd.Series, never: tuple[str, ...]) -> tuple[tuple[str, ...], ...]:
    """One row per boxed covariate: the treated range, who is outside on it, and where they are."""
    header = ("covariate", "treated min", "treated max", "outside", "of which never-IVT centre")
    rows = []
    for covariate, (low, high) in box.items():
        outside = ~((pop[covariate] >= low) & (pop[covariate] <= high))
        rows.append((covariate, f"{low:g}", f"{high:g}", str(int(outside.sum())),
                     str(int((outside & pop["center"].isin(never)).sum()))))
    rows.append(("ALL (the box)", "—", "—", str(int((~inside).sum())),
                 str(int((~inside & pop["center"].isin(never)).sum()))))
    return (header, *rows)


def _baseline_table(rows: tuple[balance.CovariateBalance, ...]) -> tuple[tuple[str, ...], ...]:
    """The baseline table, with `SMD_THRESHOLD` as the verdict column and `nan` printed as `missing`.

    The verdict is against `C.SMD_THRESHOLD` for the non-grouping rows only: the four centre rows
    report a DEFINITION and not a difference, so flagging them would be flagging the grouping variable
    for being the grouping variable.
    """
    header = ("covariate", "role", "n", "pooled SD", "SMD", f"|SMD| > {C.SMD_THRESHOLD:g}")
    out = []
    for row in rows:
        if row.role == _GROUPING:
            verdict = "the GROUPING variable — a definition, not imbalance [§10.3]"
        elif not np.isfinite(row.unweighted):
            verdict = "undefined [Stage 7 §5.3]"
        else:
            verdict = "yes" if abs(row.unweighted) > C.SMD_THRESHOLD else "no"
        out.append((
            row.covariate, row.role, str(row.n),
            "missing" if not np.isfinite(row.sd) else f"{row.sd:.6g}",
            "missing" if not np.isfinite(row.unweighted) else f"{row.unweighted:+.6f}",
            verdict))
    return (header, *out)


def support(pop: pd.DataFrame, audit: Audit) -> Support:
    """[§14a]'s support check: the treated-arm box, the never-IVT proportion, the baseline table.

    [§14a]'s framing is the one worth keeping in mind here: *"Standardisation avoids infinite weights,
    not extrapolation."* There is no weight to blow up, so nothing in the arithmetic complains when a
    prediction is made far outside the region where any comparison was observed; this is the only
    thing that makes it visible.

    T9. Records entries 5 and 6. Calls `cohort.treating_centres` for its COMPLEMENT and `balance.smd`
    and `balance.levels` for the table; implements none of them.
    """
    never = _never_ivt(pop)                                      # T9
    box = _box(pop)
    inside = _inside(pop, box)
    at_never = pop["center"].isin(never)
    baseline = _baseline(pop, never)

    outside_by_centre = {
        centre: int((~inside & (pop["center"] == centre)).sum()) for centre in C.CENTER_ORDER}
    n_never_ivt = int(at_never.sum())
    n_never_outside = int((~inside & at_never).sum())

    # ENTRY 5.
    treated_outside = int((~inside & (pop[C.TREATMENT] == max(C.TREATMENT_LABELS))).sum())
    audit.record(
        "model", "treated_support", len(pop), "\n".join([
            f"[§14a]'s support check over C.SUPPORT_COVARIATES: the per-covariate treated-arm "
            f"[min, max] box. {int(inside.sum())} inside, {int((~inside).sum())} outside.",
            f"Outside by centre: {outside_by_centre}. {treated_outside} of them are treated, which "
            "is arithmetic and not a finding — the box is the treated arm's own range, so no treated "
            "patient can be outside it.",
            f"[§14a]'S REQUESTED PROPORTION: {n_never_outside} of {n_never_ivt} patients at the "
            f"never-IVT centre(s) {never} are outside the treated support"
            f"{f' — {100.0 * n_never_outside / n_never_ivt:.1f}%' if n_never_ivt else ''}.",
            "The box is the WEAKEST reading of 'outside the treated support' and needs no tuning "
            "parameter. A convex hull in seven dimensions on this many treated patients is a set "
            "almost every point is outside of, which would report 100% and mean nothing [§10.1].",
            f"All {len(C.FACTOR_LEVELS['onset_type'])} declared `onset_type` level(s) present in the "
            "treated arm: a categorical level unseen there is a DIFFERENT kind of unsupported "
            "prediction — its coefficient does not exist — and `model.design` dropping the column is "
            "what reports it, through `Standardisation.dropped` [§10.1].",
        ]), table=_box_table(pop, box, inside, never))

    # ENTRY 6.
    n_treated = int((pop[C.TREATMENT] == max(C.TREATMENT_LABELS)).sum())
    others = len(pop) - n_treated - n_never_ivt
    worst = max((abs(r.unweighted) for r in baseline
                 if r.role != _GROUPING and np.isfinite(r.unweighted)), default=float("nan"))
    audit.record(
        "model", "support_baseline", n_treated + n_never_ivt, "\n".join([
            f"[§14a]'s baseline table: {n_treated} IVT-treated patient(s) at any centre against "
            f"{n_never_ivt} eligible patient(s) at the never-IVT centre(s) {never}.",
            f"THE OTHER {others} CONTROL PATIENT(S) AT TREATING CENTRES ARE IN NEITHER GROUP. They "
            "are not what [§14a] asked to compare, and this line is here so the denominator is never "
            "inferred from the table.",
            f"Rows range over C.BALANCE_SET ({len(baseline)} declared level(s)) and the SMD is "
            "`balance.smd` with UNIT weights — [§9] judges balance against the full confounder set "
            "regardless of what a specification fitted, and there is no weighting here to judge.",
            f"The four `center` rows carry role {_GROUPING!r}: the never-IVT group IS that centre, so "
            "those rows report a definition rather than a difference. They stay in the table because "
            "a reader must see that the comparison is between disjoint centre sets [§10.3].",
            f"`center = {never[0] if never else '?'}` is `nan` by construction — constant in each "
            "group with the groups differing, so the pooled SD is exactly zero: Stage 7 §5.3's "
            "branch 4, reported rather than suppressed.",
            f"Largest non-grouping |SMD|: {worst:.4f} against SMD_THRESHOLD = {C.SMD_THRESHOLD:g}. "
            "THAT IS THE TRANSPORT ASSUMPTION'S SIZE STATED AS A NUMBER, and [§9]'s amendment of "
            "2026-08-24 requires it beside the [§14a] estimate rather than in a supplement.",
        ]), table=_baseline_table(baseline))

    return Support(
        box=box, inside=inside, outside_by_centre=outside_by_centre, never_ivt=never,
        n_never_ivt=n_never_ivt, n_never_ivt_outside=n_never_outside, baseline=baseline)


# --- sensitivity 2, the random centre intercept [§14a, §12] ---------------------------------------

def _intercept_table(pop: pd.DataFrame, fit: model.RIFit) -> tuple[tuple[str, ...], ...]:
    """One row per centre: its size, its treated count, and the precision of its own intercept.

    **The posterior SDs are here and the intercept VALUES are not, and §4.3 is why.** A centre random
    intercept is a centre-level outcome contrast on four named hospitals, and this analysis is not
    powered to make one and does not claim to — so the values belong in the gitignored log and in the
    manuscript, and the *precisions* are what a reader needs to judge the arm.

    This is where [§14a]'s *"a centre variance from four clusters is fragile"* becomes specific: the
    never-IVT centre's intercept is the least precisely estimated of the four and it is the one the
    whole arm rests on, because it is the centre whose IVT counterfactual is entirely borrowed.
    """
    header = ("centre", "n", "treated", "posterior SD of b_hat")
    rows = []
    for index, centre in enumerate(fit.groups):
        at = pop[pop[C.BOOT_STRATUM] == centre]
        rows.append((centre, str(len(at)),
                     str(int((at[C.TREATMENT] == max(C.TREATMENT_LABELS)).sum())),
                     f"{fit.b_sd[index]:.4f}"))
    return (header, *rows)


def hierarchical(pop: pd.DataFrame, audit: Audit) -> Hierarchical:
    """[§14a]'s sensitivity 2: one common treatment effect over a random centre intercept.

        logit P(Y <= k | x, b_c)  =  alpha_k + x'beta + b_c ,      b_c ~ N(0, sigma^2)

    [§14a] prescribes *"one common treatment effect; no treatment-by-centre interaction and no random
    slope, neither being identifiable in any useful sense with four centres"*, and both prohibitions
    are structural rather than checked: neither term is ever built.

    **The standardisation is CONDITIONAL on each centre's own `b_hat_c`, not marginal over
    `N(0, sigma^2)`** (§12.6). A random-intercept model offers both and they answer different
    questions; [§14a] prescribes the second in a sentence that is easy to read past — *"Patients at a
    never-IVT centre inform that centre's intercept through their direct-EVT outcomes while the IVT
    effect is borrowed"* — and under the marginal reading those patients would inform `sigma` and
    nothing else, which would make that sentence false. The marginal alternative is a different
    estimand and is a [§14] question (§22).

    **A fit at `POLR_RI_SIGMA_FLOOR` is an ANSWER and not a failure**, and this function does not
    branch on it: [§14a] names sigma^2_C = 0 as legitimate — *"it collapses to the pooled model"* —
    and `at_floor` travels on the record so Stage 14 can print the rate beside the interval, which is
    what makes a lower limit sitting at a floor readable as "the boundary" rather than as a value.

    Records entries 8 and 9. The one caller of `model.polr_ri`.
    """
    X, dropped = model.design(pop, _DESIGN_COLUMNS)
    _assert_exposure_survived(X, dropped)                        # T1
    y = pop[C.PRIMARY_OUTCOME].to_numpy(dtype=float)
    groups = pop[C.BOOT_STRATUM]
    fit = model.polr_ri(X, y, groups)
    std = _standardise(X, fit, np.ones(len(X), dtype=bool), _RANDOM_INTERCEPT, dropped,
                       groups=groups)                            # T10, T11

    # ENTRY 8.
    audit.record(
        "model", "hierarchical_fit", len(pop), "\n".join([
            f"[§14a] sensitivity 2: a random centre intercept over {len(fit.groups)} centre(s), one "
            f"common treatment effect, no interaction and no random slope — [§14a] and [§15] forbid "
            f"both at this cluster count and neither term is ever built.",
            f"sigma_hat {fit.sigma:.6g}; at the POLR_RI_SIGMA_FLOOR = "
            f"{C.POLR_RI_SIGMA_FLOOR:g} boundary: {fit.at_floor}.",
            f"A FIT AT THE FLOOR IS NOT A FAILURE: [§14a] names sigma^2_C = 0 as a legitimate answer "
            "and it means the arm has collapsed to the pooled model. Dropping such replicates would "
            "select the bootstrap on the value of the parameter this arm exists to examine [§12.5].",
            f"ADAPTIVE Gauss-Hermite quadrature at POLR_RI_NODES = {fit.nodes} node(s). Non-adaptive "
            "quadrature is off by 0.48 log-likelihood units at 31 nodes at sigma = 3 and is "
            "NON-MONOTONE in the node count; the node count is part of the definition of the "
            "objective and is prespecified for that reason [§12.2, §12.3].",
            f"Converged in {fit.iterations} BFGS iteration(s) on the {fit.converged_on!r} criterion; "
            f"{fit.rescales} rescale(s), {fit.halvings} halving(s) in total across the fit — BFGS "
            "proposes long steps early and the line search shortens them, which is the algorithm "
            "working [§12.4].",
            "The intercept VALUES are centre-level outcome contrasts on named hospitals and are "
            "deliberately absent from any version-controlled document; their PRECISIONS are the "
            "table below, because they are what a reader needs to judge the arm [§4.3, §12.6].",
        ]), table=_intercept_table(pop, fit))

    # ENTRY 9.
    audit.record(
        "model", "standardised_hierarchical", std.n_average, "\n".join([
            f"[§14a] sensitivity 2's standardisation, over {std.n_average} record(s).",
            "IT CONDITIONS ON b_hat_c — each patient is predicted with the conditional mode of the "
            "centre they are actually at, which is [§14a]'s prescribed reading and not the marginal "
            "one. A random intercept is a per-centre shift of the cutpoints, so this is "
            "`ordinal_probabilities` unchanged against a per-centre alpha [§12.6].",
            f"conditional log-odds {std.conditional_log_odds:+.6f} — conditional on X AND on b_c "
            f"here, which the pooled arm's is not: {std.measure}.",
        ]), table=_distribution_table(std))

    return Hierarchical(
        standardisation=std, sigma=fit.sigma, at_floor=fit.at_floor,
        intercepts={c: float(fit.b[i]) for i, c in enumerate(fit.groups)},
        posterior_sd={c: float(fit.b_sd[i]) for i, c in enumerate(fit.groups)},
        nodes=fit.nodes, fit=fit)


# --- inference [§10, §14, §13] ---------------------------------------------------------------------

_ARMS: Final[tuple[str, ...]] = ("", "support.", "hier.")


def _estimand_keys() -> tuple[str, ...]:
    """§3.3's sixty-nine keys. **Standardised quantities are per-ARM; model parameters are per-FIT.**

    That is not a tidying principle, it is [§14a]'s structure read literally::

        rd_0 … rd_5                            len(C.MRS_THRESHOLDS)                      6
        mrs_0_2, mortality                     [§14a]'s two namings                       2
        dist1_0 … dist1_6, dist0_0 … dist0_6   2 * len(C.MRS_LEVELS)                     14
                                                                                        ───
                               THE STANDARDISED BLOCK, per arm                          22
        arms: "" (pooled), "support.", "hier."                                    22 * 3 = 66
        beta, hier.beta        one conditional log-odds PER FIT, and there are two        2
        hier.sigma             the between-centre SD                                     1
                                                                                        ───
                                                                                        69

    **THERE IS NO `support.beta`, AND ITS ABSENCE IS [§14a]'s RATHER THAN AN ECONOMY.** [§14a]'s
    Report sentence asks for the two standardised distributions and, FROM THEM, the `RD_k`, the mRS
    0-2 risk difference and the mortality difference; `exp(beta)` is granted a separate and narrower
    permission — it *"may appear as a model parameter but never as the standardised marginal
    effect."* And [§14a]'s sensitivity restricts THE STANDARDISATION POPULATION, which §11 implements
    by restricting the averaging set and not the fit. So the support arm has no `beta` of its own: a
    `support.beta` would be the pooled `beta`, bit-identical draw for draw with a bit-identical
    interval, presented as a property of a STANDARDISATION — the one reading [§14a] forbids in as many
    words. `hier.beta` exists because the random-intercept arm IS a second fit.

    **The counts are computed from `C.MRS_LEVELS` and `C.MRS_THRESHOLDS` and never written as 6, 14 or
    22.** Stage 8 §12 made the level set derived from the plausible range so the two cannot disagree;
    this key set is the next link in that chain, and §20.3 asserts the total against the config rather
    than against 69.

    **The key NAMES `dist1_*` and `dist0_*` are fixed literal strings** — the wire format Stage 14
    reads — and are deliberately NOT derived from `C.TREATMENT_LABELS`, even though §7.1 forbids
    writing the arm CODES as 1 and 0 anywhere in the arithmetic. The two rules do not conflict: a code
    is a value the model sees and must come from the registry, while a key name is an identifier a
    downstream reader matches on, and deriving it from the registry would let a label edit silently
    rename a reported key.
    """
    block = (
        *(f"rd_{k}" for k in C.MRS_THRESHOLDS),
        "mrs_0_2", "mortality",
        *(f"dist1_{level}" for level in C.MRS_LEVELS),
        *(f"dist0_{level}" for level in C.MRS_LEVELS),
    )
    keys = [prefix + key for prefix in _ARMS for key in block]
    keys.extend(("beta", "hier.beta", "hier.sigma"))
    return tuple(keys)


def _standardised_values(std: Standardisation, prefix: str) -> dict[str, float]:
    """One arm's twenty-two standardised draws, keyed as `_estimand_keys` names them."""
    treated, control = max(C.TREATMENT_LABELS), min(C.TREATMENT_LABELS)
    values = {f"{prefix}rd_{k}": v for k, v in std.rd.items()}
    values[f"{prefix}mrs_0_2"] = std.mrs_0_2
    values[f"{prefix}mortality"] = std.mortality
    for level in C.MRS_LEVELS:
        values[f"{prefix}dist1_{level}"] = std.distribution[treated][level]
        values[f"{prefix}dist0_{level}"] = std.distribution[control][level]
    return values


def _replicate(draw: pd.DataFrame, keys: tuple[str, ...]) -> bootstrap.Replicate:
    """One replicate: three arms off ONE drawn frame, four failure groups. §13.2.

    One drawn frame feeds all three arms, so the three are comparable draw for draw and the resample
    stream is paid for once.

    **FOUR FAILURE GROUPS, and the granularity is Stage 10 §7's extended by one level rather than
    reinvented:**

      * T1 or a pooled `polr` `FitError` costs EVERY key, `hier.*` INCLUDED. The two carry different
        buckets — `degenerate_design` and `nonconvergence` — which is why they are two groups and not
        one.
      * A `polr_ri` `FitError` costs `hier.*` only: twenty-two standardised keys plus `hier.beta` and
        `hier.sigma`.
      * T12 costs ONE `beta` at whichever fit raised it, and nothing else (§9.1).

    **A POOLED FAILURE COSTS THE HIERARCHICAL ARM, AND THE REASON IS §12.8's RATHER THAN A
    CONSERVATISM.** `polr_ri` starts from the pooled fit's `alpha` and `beta` — a nested model's exact
    maximiser in every coordinate but one, which is why 25 iterations suffice — so **there is no
    hierarchical arm without a pooled fit.** An exemption would require a cold start whose iteration
    count, convergence route and boundary rate are all unmeasured, and a quarter of the draws would
    then sit in a different numerical regime from the rest, which is precisely the failure §12.2
    rejects non-adaptive quadrature for. Measured: the pooled fit fails in 0 of 2000 replicates, so
    this path never fires on this workbook, and it is specified anyway on §5.3's rule.

    `bootstrap.bucket` classifies and `bootstrap.collect` reconciles; nothing here re-implements
    either, which is what Stage 9 §14 names as how a taxonomy gets violated by accident.

    **`support` IS NOT CALLED HERE AND THE BOX IS RECOMPUTED INSTEAD** (§11, `_never_ivt`). The
    support is a function of the treated arm's observed range, so it is a STATISTIC; holding it at the
    point estimate's would condition the bootstrap on one realisation of it and understate the
    interval. Measured across 2000 replicates the in-support averaging population ranges 69 to 102
    against 93 at the point estimate — a spread wide enough that the choice is not cosmetic.
    """
    hierarchical_keys = tuple(k for k in keys if k.startswith("hier."))
    values: dict[str, float] = {}
    failures: dict[str, str] = {}
    y = draw[C.PRIMARY_OUTCOME].to_numpy(dtype=float)

    try:
        X, dropped = model.design(draw, _DESIGN_COLUMNS)
        _assert_exposure_survived(X, dropped)                    # T1 -> degenerate_design
        fit = model.polr(X, y)                                   # -> nonconvergence
    except model.FitError as failure:
        label = bootstrap.bucket(str(failure))
        return bootstrap.Replicate(
            values={}, failures={key: label for key in keys},
            n_alpha=None, polr_iterations=None, sum_w=None, n_in_model=None)

    treatment = fit.columns.index(C.TREATMENT)
    try:
        _assert_beta_reportable(float(fit.beta[treatment]), "beta")   # T12 -> separation
        values["beta"] = float(fit.beta[treatment])
    except model.FitError as failure:
        failures["beta"] = bootstrap.bucket(str(failure))

    values.update(_standardised_values(
        _standardise(X, fit, np.ones(len(X), dtype=bool), _ALL_ELIGIBLE, dropped), ""))
    values.update(_standardised_values(
        _standardise(X, fit, _inside(draw, _box(draw)).to_numpy(dtype=bool),
                     _TREATED_SUPPORT, dropped), "support."))

    # `n_alpha` and `polr_iterations` carry the POOLED fit's, which is the model [§14a] names, and
    # `polr_ri`'s are deliberately NOT forced into them (§13.2). `sum_w` is `Sw = n` because [§14a]
    # weights nobody — it is the number that closed the absolute-tolerance `TODOS` item negatively
    # (§13.4). `max_abs_beta` is keyed by FIT, which is Stage 10's own meaning of that field over a
    # different index: "the largest |beta| of a fitted model, keyed by what was fitted".
    diagnostics = dict(n_alpha=len(fit.alpha), polr_iterations=fit.iterations,
                       sum_w=float(len(draw)), n_in_model=len(draw),
                       max_abs_beta={"pooled": float(np.max(np.abs(fit.beta)))})
    try:
        ri = model.polr_ri(X, y, draw[C.BOOT_STRATUM])
    except model.FitError as failure:
        label = bootstrap.bucket(str(failure))
        failures.update({key: label for key in hierarchical_keys})
        return bootstrap.Replicate(values=values, failures=failures, **diagnostics)
    diagnostics["max_abs_beta"] = dict(
        diagnostics["max_abs_beta"], hier=float(np.max(np.abs(ri.beta))))

    try:
        _assert_beta_reportable(float(ri.beta[treatment]), "hier.beta")   # T12, the second site
        values["hier.beta"] = float(ri.beta[treatment])
    except model.FitError as failure:
        failures["hier.beta"] = bootstrap.bucket(str(failure))
    values["hier.sigma"] = ri.sigma
    values.update(_standardised_values(
        _standardise(X, ri, np.ones(len(X), dtype=bool), _RANDOM_INTERCEPT, dropped,
                     groups=draw[C.BOOT_STRATUM]), "hier."))
    return bootstrap.Replicate(values=values, failures=failures, **diagnostics)


def diagnostics(collected: tuple[object, ...]) -> bootstrap.Diagnostics:
    """`bootstrap.Diagnostics` over a §14 stage's replicates. Public: Stage 13 calls it (its §10.2).

    `bootstrap._diagnostics` is not reused and is not merely private: it aggregates `max_abs_beta` and
    `or_corrected`, which are Stage 9's augmented-estimator and continuity-correction diagnostics, and
    **this stage has neither** — no augmented outcome and no odds ratio to correct. Both come back
    empty here, which is a fact about [§14a] rather than a gap.

    `sum_w` is the replicate's ROW COUNT, because [§14a] weights nobody: `Sw = n`. It is carried
    anyway because Stage 8 §11 asked for it and because it is the number that closed the absolute-
    tolerance `TODOS` item negatively — `Sw` is 104 here against roughly 23 for the [§7]
    overlap-weighted cohort, which is the change of scale that item names as its trigger (§13.4).
    """
    live = [r for r in collected if r.n_alpha is not None]
    magnitudes: dict[str, list[float]] = {}
    for replicate in live:
        for key, value in replicate.max_abs_beta.items():
            magnitudes.setdefault(key, []).append(value)
    return bootstrap.Diagnostics(
        n_alpha=_tally(r.n_alpha for r in live),
        polr_iterations=_tally(r.polr_iterations for r in live),
        sum_w=np.asarray([r.sum_w for r in live], dtype=float),
        max_abs_beta={key: np.asarray(values, dtype=float)
                      for key, values in magnitudes.items()},
        n_in_model=_tally(r.n_in_model for r in live),
        or_corrected={})


def _tally(values) -> dict[int, int]:
    """value -> replicate count, ascending. `Diagnostics`' three distributions have this shape."""
    counts: dict[int, int] = {}
    for value in values:
        counts[int(value)] = counts.get(int(value), 0) + 1
    return dict(sorted(counts.items()))


def _replicates_table(draws: dict[str, bootstrap.Draws],
                      diagnostics: bootstrap.Diagnostics,
                      collected: tuple[object, ...]) -> tuple[tuple[str, ...], ...]:
    """Entry 10's diagnostic grid, and **every block labels its own denominator.**

    That is a `TODOS` item filed against Stage 10's equivalent: this grid renders distributions over
    three different populations — replicates attempted, replicates surviving the pooled fit, and
    replicates reaching `polr_ri` — and a grid whose blocks share one unlabelled denominator invites
    a reader to divide by the wrong one. Stage 12 does not fix Stage 10's grid; it declines to add a
    second unlabelled one.
    """
    attempted = max((d.n_attempted for d in draws.values()), default=0)
    live = len(diagnostics.sum_w)
    hier = len(draws["hier.sigma"].draws)
    at_floor = sum(1 for r in collected if r.values.get("hier.sigma") == C.POLR_RI_SIGMA_FLOOR)
    buckets = tuple(sorted(set(C.FAILURE_BUCKETS.values())))

    rows: list[tuple[str, ...]] = [("block / denominator", "quantity", "value")]
    rows.append((f"attempted = {attempted}", "estimand keys", str(len(draws))))
    for name in buckets:
        total = sum(d.failures.get(name, 0) for d in draws.values())
        rows.append((f"attempted = {attempted}", f"key-failures bucketed {name!r}", str(total)))
    rows.append((f"reached the pooled fit = {live}", "cutpoint counts len(alpha)",
                 str(diagnostics.n_alpha)))
    rows.append((f"reached the pooled fit = {live}", "polr iterations",
                 str(diagnostics.polr_iterations)))
    rows.append((f"reached the pooled fit = {live}", "Sw = n (unit weights)",
                 f"{diagnostics.sum_w.min():g} to {diagnostics.sum_w.max():g}"
                 if live else "—"))
    rows.append((f"reached polr_ri = {hier}", "at the sigma floor",
                 f"{at_floor} ({100.0 * at_floor / hier:.1f}%)" if hier else "—"))
    rows.append((f"reached polr_ri = {hier}", "sigma_hat range",
                 f"{draws['hier.sigma'].draws.min():.6g} to "
                 f"{draws['hier.sigma'].draws.max():.6g}" if hier else "—"))
    for key, spread in sorted(diagnostics.max_abs_beta.items()):
        rows.append((f"reached the {key} fit = {len(spread)}", f"max |coefficient|, {key}",
                     f"{spread.min():.4f} to {spread.max():.4f} against POLR_MAX_ABS_BETA = "
                     f"{C.POLR_MAX_ABS_BETA:g}"))
    # THE SUPPORT ARM's `n_average` SPREAD IS NOT HERE, and its absence is `Replicate`'s rather than
    # this table's: that record carries the POOLED fit's diagnostics by §13.2, and the averaging size
    # of a second arm is not one of its fields. §21 item 4 declines widening it, so the spread is
    # measured by a probe (69 to 102 over 2000 replicates, median 91) and stated in §11 rather than
    # rendered here. A row reading "not carried" would be a row nobody can act on.
    return tuple(rows)


def inference(pop: pd.DataFrame, audit: Audit) -> bootstrap.Bootstrap:
    """The [§14] bootstrap: one `replicates` call, three arms per draw, sixty-nine keys, NO p-value.

    *"Patient-level bootstrap stratified by centre, as §10: resample, refit the ordinal model, predict
    every resampled patient under both regimes, average, recompute the risk differences, percentile
    intervals."*

    **ALL FOUR STRATA ARE RESAMPLED, the never-IVT centre included, and Stage 10 §12.2 left this
    decision here** in as many words: *"whether that is resampled or held fixed is [§14a]'s question
    and not this document's."* The argument is about what the estimand is. [§14a]'s target is an
    average over the FULL eligible cohort at all centres, and those patients are members of that
    population — they contribute their covariates to the average and their observed outcomes to the
    fit. The uncertainty in an average over a sample includes the uncertainty in WHICH sample was
    drawn. Holding the stratum fixed would produce an interval conditional on those covariate vectors,
    which is a different and narrower inferential statement, and it would be inconsistent with the
    other three strata, which are resampled for exactly the reason this one would not be.

    The counter-argument is recorded because it is not empty: that centre contributes no contrast, so
    resampling it adds variance to the TARGET POPULATION without adding information about the EFFECT.
    That is a real distinction and it is not [§14a]'s, which defines its estimand as an average over a
    population and asks for a bootstrap of it. §22 records that Stage 12 does not decide whether a
    population-conditional variant should exist.

    **NO KEY CARRIES A p-VALUE, and the rule is passed to `bootstrap.intervals` as
    `lambda key: False`.** [§14]'s Inference paragraph prescribes *"percentile intervals"* and nothing
    else; [§10]'s p is defined for the treatment coefficient of the [§8] WEIGHTED proportional-odds
    model, which is a different estimator on a different population, so a p-value here would be a test
    [§14] does not ask for against a null [§14] does not state.

    **Nothing in `bootstrap.resample` changes.** It is general in frame and stratum by Stage 10 §4's
    design, and it groups with `sort=True, observed=True`, so a fourth stratum appears simply because
    the frame has one.

    Recomputes the three arms inside the replicate body rather than taking the point-estimate records,
    for Stage 10 §11's reason: a replicate must run the whole procedure, and a body that took a fitted
    object would be resampling around a fixed fit.

    Records entry 10. Calls `bootstrap.replicates`, `bootstrap.bucket`, `bootstrap.collect` and
    `bootstrap.intervals`; implements none of them.
    """
    keys = _estimand_keys()
    collected = bootstrap.replicates(
        pop, lambda draw: _replicate(draw, keys), C.N_BOOT, C.SEED, C.BOOT_STRATUM)
    draws = bootstrap.collect(collected, keys)
    diag = diagnostics(collected)
    limits = bootstrap.intervals(draws, lambda key: False)       # §3.3 — no p on ANY key

    # ENTRY 10.
    audit.record(
        "model", "standardisation_replicates", C.N_BOOT, "\n".join([
            f"[§14] inference: {C.N_BOOT} patient-level replicate(s) stratified by "
            f"{C.BOOT_STRATUM!r}, ALL FOUR STRATA including the never-IVT centre — Stage 10 §12.2 "
            f"left that decision to [§14a] and §13.1 is where it is made.",
            f"{len(keys)} estimand key(s); {len(limits)} interval(s) at level {C.CI_LEVEL:g} under "
            f"the {C.PERCENTILE_METHOD!r} percentile definition. NO KEY CARRIES A p-VALUE: [§14] "
            "prescribes percentile intervals and states no null.",
            f"Seed {C.SEED}; the floor below which an estimand keeps its draws and gets no interval "
            f"is {C.ci_min_draws(C.CI_LEVEL)}.",
            "Three arms off ONE drawn frame, so the three are comparable draw for draw. Four failure "
            "groups: T1 and a pooled non-convergence each cost every key (a pooled failure costs "
            "`hier.*` too, because `polr_ri`'s start values ARE the pooled fit's [§12.8, §13.2]); a "
            "`polr_ri` failure costs `hier.*`; T12 costs one `beta` at whichever fit raised it.",
            "The treated-support box is RECOMPUTED in every replicate and is not held at the point "
            "estimate's: the support is a statistic, and holding it fixed would condition the "
            "bootstrap on one realisation of it and understate the interval [§11].",
            "`mortality` carries its OWN draws and its own percentile call rather than the reflected "
            f"limits of `rd_5`. Under {C.PERCENTILE_METHOD!r} the limits of -X and the reflected "
            "limits of X disagree by up to 1.7e-3 on the risk-difference scale — 0.17 percentage "
            "points, which a manuscript prints — while under numpy's default they agree to 5.6e-16. "
            "The pin that makes [§10]'s p-value agree with its interval is the same pin that breaks "
            "reflection symmetry [§7.4].",
        ]), table=_replicates_table(draws, diag, collected))

    return bootstrap.Bootstrap(C.SEED, C.N_BOOT, draws, limits, diag)
