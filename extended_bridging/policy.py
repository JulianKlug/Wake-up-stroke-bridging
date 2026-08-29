"""Stage 13 — [§14b] the feasible-policy contrast.

One pooled proportional-odds model over ALL classified patients at all four centres — eligible and
contraindicated — with [§14a]'s covariates plus a contraindication indicator (DECISION 10);
g-computation to two REGIMES, the active one bridging every eligible patient and giving EVT alone to
every contraindicated one, the comparator giving EVT alone to everyone; and percentile intervals from
a patient-level bootstrap stratified by centre.

Section references in brackets are to ``statistical_analysis_plan.md``. The specification for this
module is ``specs/stage13_feasible_policy.md``; nothing here is invented outside it.

**THE ESTIMAND IS AN OPERATIONAL POLICY CONTRAST — what adopting a bridging policy would have
delivered in this cohort — and is never [§14a]'s ATE, never [§7]'s ATO, and never the biological
effect of IVT.** That sentence travels on every `Policy` as `label` (§8) and Stage 14 prints it from
the record.

**The population is nobody else's.** NEITHER [§3] restriction is applied: the never-IVT centre stays
and the contraindicated patients stay. `cohort.build` applies both and `standardise.population`
applies the second, so neither is called (§0.2). `eligibility.retained` is READ — it supplies the
regime and the indicator — and never applied as a filter.

**The structural fact of this stage (§7.3).** A contraindicated patient receives EVT alone under BOTH
regimes, so their per-patient contrast is exactly zero and every `rd_k` is the eligible-subpopulation
contrast diluted by the eligible share::

    rd_k  =  (n_eligible / N) * eligible_rd_k

U7 asserts it on every record this module builds, from two independent `gcompute` calls.

No propensity model is fitted anywhere here — [§14]'s own first sentence — and `propensity.py`,
`cohort.py` and `derive.py` are not imported.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Final

import numpy as np
import pandas as pd

import bootstrap
import config as C
import eligibility
import model
import standardise
from data import Audit, absence_by_column, coefficient_rows, replicate_rows

# --- the two guard strings, and they are DATA rather than docstrings [§8] -------------------------
#
# Stage 12 carries `measure` so that `exp(beta)` is never labelled marginal. This stage has a second
# thing that must never be mislabelled — the ESTIMAND — and carries `label` for it. Both are `Final`
# and a second value can only arrive with a [§14] amendment behind it.
_CONDITIONAL_14B: Final[str] = (
    "conditional odds ratio [§14b] — conditional on X AND on contraindication status, and never "
    "the standardised marginal effect")
_OPERATIONAL: Final[str] = (
    "operational policy contrast [§14b] — what adopting a bridging policy would have delivered in "
    "this cohort; not [§14a]'s ATE in the eligible population, not [§7]'s ATO, and not the "
    "biological effect of IVT")

# The two REGIME names. Not arm codes: under the active regime 104 patients receive bridging and 19
# receive EVT alone, and a dict keyed `1` would say every patient was bridged (§3.1). They are also
# the estimand-key stems (§3.3) — the wire format Stage 14 matches on, deliberately distinct from
# Stage 12's `dist1_`/`dist0_`.
_ACTIVE: Final[str] = "active"
_COMPARATOR: Final[str] = "comparator"

# U7's tolerance, and it is T10's rather than "two orders above the measured worst" (5.3e-16): a
# naive `N * eps` rounding bound on 123 rows summed in two different orders is already 1.4e-14, so a
# tighter contract could fail on arithmetic that is correct (§7.3). 1e-12 is eight orders below
# anything printed.
_DILUTION_TOL: Final[float] = 1e-12

_DESIGN_COLUMNS: Final[tuple[str, ...]] = (C.TREATMENT,) + C.POLICY_COVARIATES

# The sentence entry 1 carries when the cohort has nobody to withhold IVT from (§11). A data fact,
# recorded, and NOT a contract break: `contrast` proceeds with `share_eligible == 1`.
_NO_CONTRAINDICATED: Final[str] = (
    "[§14b] has no contraindicated patient in this cohort: the active regime bridges everyone, "
    "`share_eligible` is 1, and this policy contrast coincides with a [§14a]-shaped contrast under a "
    "fit with no indicator; read it with that stated.")

# The label `gcompute`'s T10/T11 messages name this stage's records by.
# Public alias for the reporting layer: T18 prints this sentence when `share_eligible` is 1 [Stage 14].
NO_CONTRAINDICATED: Final[str] = _NO_CONTRAINDICATED

_GCOMPUTE_LABEL: Final[str] = "[§14b] policy"


# --- what the stage returns [§3.1] ----------------------------------------------------------------

@dataclass(frozen=True)
class Policy:
    """The [§14b] feasible-policy contrast: one fit, one averaging population, two REGIMES.

    Frozen, carrying no verdict, on Stage 8 §3's rule.

    `distribution` and `cumulative` are keyed by regime NAME (`_ACTIVE`, `_COMPARATOR`) and not by arm
    code, because the active regime is not an arm (§3.1).

    **There is no field called `odds_ratio`**, for Stage 12 §8's reason. And there is no
    pre-multiplied field: `rd` IS the policy contrast, `eligible_rd` is what it dilutes,
    `share_eligible` is the factor, and U7 asserts the three agree (§7.3). `eligible_rd` is
    deliberately not called an ATE: it is the eligible-subpopulation contrast under a model fitted on
    a wider population, and [§14a]'s ATE is Stage 12's record (§8).

    `contraindication_log_odds` is a NUISANCE parameter: printed in the fit table with its role, given
    no interval, never quoted (§4.3, §9). It is `nan` when the indicator was dropped as constant.
    """

    n_average: int                       # the averaging population's size — §7.1, [§11]
    n_fit: int                           # equal to n_average; carried for Stage 14's denominators
    n_eligible: int                      # patients the active regime bridges — §7.3
    share_eligible: float                # n_eligible / n_average, the dilution factor — §7.3
    distribution: dict[str, dict[int, float]]   # _ACTIVE | _COMPARATOR -> mRS level -> P
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


# --- the indicator and the two regimes, from ONE predicate [§5.2, §6.2] ---------------------------

def _regimes(frame: pd.DataFrame) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """(D, active, comparator) for `frame`, from `eligibility.retained` and nothing else. U4.

    THE ONE PLACE the contraindication predicate is spelled. `population` stores D as the
    CONTRAINDICATED column; `contrast` and the replicate body call this on THEIR frame, and a stored
    column that disagrees with `eligibility` raises here rather than fitting one thing and
    standardising another. `D + active / treated == 1` on every row by construction.

    Not from `ivt_contraindicated`: Stage 4 is the one classifier, and a stage that re-read the flag
    would be a second one that agrees with the first on this workbook and could stop agreeing on the
    next without anything raising (§5.2).
    """
    if C.ELIGIBILITY not in frame.columns or not frame[C.ELIGIBILITY].isin(C.ELIGIBILITY_ORDER).all():
        raise C.SchemaError(
            f"U4  the frame carries no {C.ELIGIBILITY!r} column, or a value outside "
            f"ELIGIBILITY_ORDER = {C.ELIGIBILITY_ORDER}. The contraindication indicator is derived "
            "from that column and from nothing else; a comparison against a column carrying a third "
            "value would silently code it as NOT contraindicated [Stage 13 §5.2].")

    retained = eligibility.retained(frame).to_numpy(dtype=bool)
    treated, control = max(C.TREATMENT_LABELS), min(C.TREATMENT_LABELS)
    D = (~retained).astype(float)
    active = retained.astype(float) * treated + (~retained).astype(float) * control
    comparator = np.full(len(frame), float(control))

    if C.CONTRAINDICATED in frame.columns:
        disagree = int((frame[C.CONTRAINDICATED].to_numpy(dtype=float) != D).sum())
        if disagree:
            raise C.SchemaError(
                f"U4  the stored {C.CONTRAINDICATED!r} column disagrees with `eligibility` on "
                f"{disagree} row(s). The column is what U2 checks the regime against and the "
                "predicate is what builds the regime; two sources that disagree would fit one thing "
                "and standardise another [Stage 13 §5.2, §6.3].")
    return D, active, comparator


# --- the guards [§11] -----------------------------------------------------------------------------

def _assert_exposure_survived(X: pd.DataFrame, dropped: tuple[str, ...]) -> None:
    """U1 — the [§14b] design still carries the treatment column. `FitError`, never `SchemaError`.

    Stage 12's T1 restated with this stage's token, and NOT a call to it: `test_bootstrap.py`'s
    raise-site scan reads the first token of every `raise FitError(` statically and a parameterised
    token would crash it (§5.3). The indicator falls under the covariate rule: a draw with no
    contraindicated patient drops `D` as constant and is a legitimate draw, recorded through
    `Policy.dropped`.
    """
    if C.TREATMENT not in X.columns:
        raise model.FitError(
            f"U1  the [§14b] design lost the treatment column; {len(dropped)} column(s) were dropped "
            f"as constant: {', '.join(dropped) or 'none'}. A design with no exposure FITS, and "
            "returns a contrast in which both regimes are the same regime. [§10] drops and counts "
            "the replicate.")


def _assert_beta_reportable(value: float) -> None:
    """U8 — the fitted TREATMENT coefficient is below `POLR_MAX_ABS_BETA`. `FitError`, costs `beta`.

    Stage 12's T12 restated with this stage's token, for U1's reason. It binds `beta` alone and does
    not look at the contraindication coefficient (§9.2): a separated indicator absorbs the rows it
    indexes and leaves `alpha`, `beta` and `gamma` to the others, so it changes no reported quantity,
    and a bound on the worst coefficient would drop such a replicate for a reason no reported
    quantity cares about.
    """
    if abs(value) >= C.POLR_MAX_ABS_BETA:
        raise model.FitError(
            f"U8  the fitted treatment coefficient is {value:.6g}, reaching POLR_MAX_ABS_BETA = "
            f"{C.POLR_MAX_ABS_BETA:g}. A separated proportional-odds fit CONVERGES AND RETURNS "
            "[Stage 8 §6]; this bound turns it into a countable failure costing `beta` alone — the "
            "standardised quantities are averaged probabilities and stay in range whatever beta "
            "does [Stage 13 §9.2].")


def _assert_dilution(rd: dict[int, float], eligible_rd: dict[int, float],
                     share_eligible: float) -> None:
    """U7 — `rd_k == share_eligible * eligible_rd_k` at every threshold, at `_DILUTION_TOL`. §7.3.

    The two sides are computed independently, as two `gcompute` calls on the same fit, so this is a
    check between two computations and not of one against itself. `SchemaError`, for T10's reason: an
    identity that fails is arithmetic broken, not a sparse draw.
    """
    for k in C.MRS_THRESHOLDS:
        residual = abs(rd[k] - share_eligible * eligible_rd[k])
        if not np.isfinite(residual) or residual > _DILUTION_TOL:
            raise C.SchemaError(
                f"U7  rd_{k} = {rd[k]!r} against share_eligible * eligible_rd_{k} = "
                f"{share_eligible * eligible_rd[k]!r}, off by {residual:.3e} (> {_DILUTION_TOL:g}). "
                "A contraindicated patient receives EVT alone under BOTH regimes, so their per-row "
                "contrast is exactly zero and the policy contrast is the eligible contrast times the "
                "eligible share — an identity, and a failed identity is broken arithmetic "
                "[Stage 13 §7.3].")


def _assert_regime(frame: pd.DataFrame, active: np.ndarray) -> None:
    """U2 — no row the STORED column marks contraindicated carries the treated code in `active`.

    The roadmap's Accept-when — *"no contraindicated patient is ever assigned a predicted IVT
    outcome"* — as a contract on the regime vector, checked BEFORE `ordinal_probabilities` runs and
    against the stored column rather than `_regimes`' own `D`, so it compares two sources (§6.3).
    """
    control = min(C.TREATMENT_LABELS)
    stored = frame[C.CONTRAINDICATED].to_numpy(dtype=float) == 1.0
    violations = int((active[stored] != control).sum())
    if violations:
        raise C.SchemaError(
            f"U2  {violations} contraindicated row(s) carry the treated code in the ACTIVE regime. "
            "No contraindicated patient is ever assigned a predicted IVT outcome [§14b]; the regime "
            "is a definition, so this is a bug in the code that built the vector and never a "
            "property of a resample [Stage 13 §6.3].")


# --- the population [§14b, §11, §4] ---------------------------------------------------------------

def _population_table(source: pd.DataFrame, pop: pd.DataFrame) -> tuple[tuple[str, ...], ...]:
    """The per-centre x eligibility x arm ledger (§4.1's table), over CENTER_ORDER and never the frame."""
    treated, control = max(C.TREATMENT_LABELS), min(C.TREATMENT_LABELS)
    header = ("centre", "classified", "in [§14b] population",
              f"eligible, {C.TREATMENT_LABELS[control]}", f"eligible, {C.TREATMENT_LABELS[treated]}",
              f"contraindicated, {C.TREATMENT_LABELS[control]}")

    def row(name: str, source_at: pd.DataFrame, at: pd.DataFrame) -> tuple[str, ...]:
        contraindicated = at[C.CONTRAINDICATED] == 1.0
        return (
            name, str(len(source_at)), str(len(at)),
            str(int((~contraindicated & (at[C.TREATMENT] == control)).sum())),
            str(int((~contraindicated & (at[C.TREATMENT] == treated)).sum())),
            str(int((contraindicated & (at[C.TREATMENT] == control)).sum())))

    rows = [row(c, source[source["center"] == c], pop[pop["center"] == c]) for c in C.CENTER_ORDER]
    rows.append(row("all", source, pop))
    return (header, *rows)


def population(df: pd.DataFrame, audit: Audit) -> pd.DataFrame:
    """The [§14b] population: EVERY classified patient, covariate-complete, outcome-observed,
    with the CONTRAINDICATED indicator added. 126 -> 123. Neither [§3] restriction is applied.

    The [§11] mask is `complete & observed` with no eligibility conjunct — a third caller whose mask
    differs from both `cohort.build`'s and `standardise.population`'s in content, so Stage 12 §4.1's
    argument against a shared helper holds a third time.

    U3, U4, U5. Records entries 1 and 2.
    """
    D_all, _, _ = _regimes(df)                                   # U4
    complete = model.complete_cases(df, C.STANDARDISATION_COVARIATES)
    observed = df[C.PRIMARY_OUTCOME].notna()
    mask = complete & observed
    pop = df[mask].copy()
    pop[C.CONTRAINDICATED] = D_all[mask.to_numpy(dtype=bool)]

    treated, control = max(C.TREATMENT_LABELS), min(C.TREATMENT_LABELS)
    contraindicated = pop[C.CONTRAINDICATED] == 1.0
    n_treated_contraindicated = int((contraindicated & (pop[C.TREATMENT] == treated)).sum())
    if n_treated_contraindicated:
        raise C.SchemaError(
            f"U3  {n_treated_contraindicated} contraindicated patient(s) carry the treated code. "
            "Stage 4's load-bearing assertion — no treated patient carries an absolute "
            "contraindication — is re-checked at this boundary rather than trusted, because the "
            "active regime assigns EVT alone to every contraindicated patient and a treated one "
            "would make that regime contradict the data [Stage 13 §4.1].")

    eligible = pop[~contraindicated]
    arms = {C.TREATMENT_LABELS[a]: int((eligible[C.TREATMENT] == a).sum())
            for a in sorted(C.TREATMENT_LABELS)}
    if len(pop) == 0 or any(count == 0 for count in arms.values()):
        raise C.SchemaError(
            f"U5  the [§14b] population holds {len(pop)} record(s), the eligible by arm {arms}. "
            "The treatment coefficient is identified from eligible patients alone, and an empty "
            "eligible arm leaves it unidentified while every average remains computable "
            "[Stage 13 §11]. A population with NO contraindicated patient is not this error.")

    # ENTRY 1.
    per_centre = {c: int((contraindicated & (pop["center"] == c)).sum()) for c in C.CENTER_ORDER}
    lines = [
        "[§14b]'s population: NEITHER [§3] restriction is applied. The never-IVT centre stays "
        "(restriction 1) and the contraindicated patients stay (restriction 2), because [§14b] is "
        "the one §14 analysis whose population includes them — this frame is nobody else's.",
        f"{len(df)} classified -> {int(complete.sum())} covariate-complete on "
        f"STANDARDISATION_COVARIATES -> {len(pop)} with the primary outcome observed.",
        f"[§11] covariate-incomplete: {int((~complete).sum())} removed; primary outcome not "
        f"observed: {int((complete & ~observed).sum())} removed.",
        f"Contraindicated per centre: {per_centre}; {int(contraindicated.sum())} in all, none "
        "treated (U3). The indicator is derived from `eligibility` and never from the raw flag.",
        "`eligibility.retained` is READ — it supplies the regime and the indicator — and never "
        "applied as a filter. Step names are disjoint from `cohort.build`'s and `standardise`'s "
        "because all three remove rows from one classified frame under kind='cohort' and "
        "`Audit.entry` is first-match.",
    ]
    if int(contraindicated.sum()) == 0:
        lines.append(_NO_CONTRAINDICATED)
    standardise.record_removal(audit, "policy_population", df[~mask], "\n".join(lines),
                               _population_table(df, pop))

    # ENTRY 2.
    absence_by_column(
        pop, audit, list(C.POLICY_COVARIATES) + [C.PRIMARY_OUTCOME], "absence_by_policy_column",
        "Absence over the [§14b] population, complete on every covariate by construction and on "
        "the indicator by definition: the complete-case denominator [§11] requires is READ rather "
        "than inferred from the ledger above.")
    return pop


# --- the fit and the two regimes [§5, §6, §7] -----------------------------------------------------

def _estimate(pop: pd.DataFrame) -> Policy:
    """One `Policy` off one frame: the fit, U1, U2, the two `gcompute` calls, U7. No audit, no U8.

    Shared by `contrast` and the replicate body so the two cannot drift. U8 is the CALLER's, because
    the body must cost `beta` alone on it while every other field survives (§9.2, §10.2).
    """
    D, active, comparator = _regimes(pop)                        # U4
    X, dropped = model.design(pop, _DESIGN_COLUMNS)
    _assert_exposure_survived(X, dropped)                        # U1
    y = pop[C.PRIMARY_OUTCOME].to_numpy(dtype=float)
    fit = model.polr(X, y)                                       # UNWEIGHTED — no `w`, §5.1

    # THE POLICY CONTRAST: every patient, active against comparator, U2 on the stored column first.
    everyone = np.ones(len(X), dtype=bool)
    _assert_regime(pop, active)                                  # U2
    g = standardise.gcompute(X, fit, everyone, active, comparator, label=_GCOMPUTE_LABEL,
                             keys=(_ACTIVE, _COMPARATOR))        # T10, T11

    # THE UNDILUTED ELIGIBLE CONTRAST, on a design that CONTAINS NO CONTRAINDICATED ROW (§6.3): the
    # [§14a]-shaped averaging under the [§14b] fit. No IVT outcome is ever predicted for a
    # contraindicated patient, even transiently.
    eligible = D == 0.0
    treated, control = max(C.TREATMENT_LABELS), min(C.TREATMENT_LABELS)
    Xe, pop_e = X[eligible], pop[eligible]
    n_eligible = int(eligible.sum())
    _assert_regime(pop_e, np.full(n_eligible, float(treated)))  # U2, trivially — nothing to violate
    ge = standardise.gcompute(Xe, fit, np.ones(n_eligible, dtype=bool),
                              np.full(n_eligible, float(treated)),
                              np.full(n_eligible, float(control)), label=_GCOMPUTE_LABEL)
    share_eligible = n_eligible / g.n_average
    _assert_dilution(g.rd, ge.rd, share_eligible)               # U7

    beta = float(fit.beta[fit.columns.index(C.TREATMENT)])
    delta = (float(fit.beta[fit.columns.index(C.CONTRAINDICATED)])
             if C.CONTRAINDICATED in fit.columns else float("nan"))
    return Policy(
        n_average=g.n_average, n_fit=g.n_fit, n_eligible=n_eligible, share_eligible=share_eligible,
        distribution=g.distribution, cumulative=g.cumulative, rd=g.rd,
        mrs_0_2=g.mrs_0_2, mortality=g.mortality, eligible_rd=ge.rd,
        conditional_log_odds=beta, conditional_odds_ratio=float(np.exp(beta)),
        contraindication_log_odds=delta, measure=_CONDITIONAL_14B, label=_OPERATIONAL,
        fit=fit, dropped=dropped)


def _fit_table(fit: model.PolrFit, dropped: tuple[str, ...]) -> tuple[tuple[str, ...], ...]:
    """Entry 3's coefficient table: every design column, its coefficient, its ROLE. No `|` in a cell."""
    def role(name: str) -> str:
        if name == C.TREATMENT:
            return "EXPOSURE — the guarded one [U8]"
        if name == C.CONTRAINDICATED:
            return "contraindication indicator — a NUISANCE, unbounded [§9], not a result [§4.3]"
        return "gamma, unbounded [§9.2]"

    return coefficient_rows(fit.columns, fit.beta, dropped, role)


def _contrast_table(pol: Policy) -> tuple[tuple[str, ...], ...]:
    """Entry 4's table: the two regime distributions, the cumulative table, `rd_k`, and the
    decomposition beside it — `eligible_rd_k`, `share_eligible`, the dilution residual — as
    COMPUTED differences. Columns are headed by REGIME name, never by an arm label."""
    header = ("mRS", f"P(Y=j) {_ACTIVE}", f"P(Y=j) {_COMPARATOR}", f"P(Y<=k) {_ACTIVE}",
              f"P(Y<=k) {_COMPARATOR}", "rd_k", "eligible_rd_k", "abs(rd_k - share x eligible_rd_k)")
    rows = []
    for level in C.MRS_LEVELS:
        tail = (
            (f"{pol.cumulative[_ACTIVE][level]:.6f}", f"{pol.cumulative[_COMPARATOR][level]:.6f}",
             f"{pol.rd[level]:+.6f}", f"{pol.eligible_rd[level]:+.6f}",
             f"{abs(pol.rd[level] - pol.share_eligible * pol.eligible_rd[level]):.3e}")
            if level in pol.rd else ("1.000000", "1.000000", "—", "—", "—"))
        rows.append((str(level), f"{pol.distribution[_ACTIVE][level]:.6f}",
                     f"{pol.distribution[_COMPARATOR][level]:.6f}", *tail))
    blank = ("", "", "", "")
    rows.append(("share_eligible", *blank, f"{pol.share_eligible:.6f}",
                 f"{pol.n_eligible} of {pol.n_average}", ""))
    rows.append(("mRS 0-2", *blank, f"{pol.mrs_0_2:+.6f}", "", ""))
    rows.append(("mortality", *blank, f"{pol.mortality:+.6f}", "", ""))
    rows.append(("abs(mrs_0_2 - rd_2)", *blank, f"{abs(pol.mrs_0_2 - pol.rd[2]):.3e}", "", ""))
    rows.append(("abs(mortality + rd_5)", *blank, f"{abs(pol.mortality + pol.rd[5]):.3e}", "", ""))
    return (header, *rows)


def contrast(pop: pd.DataFrame, audit: Audit) -> Policy:
    """The [§14b] feasible-policy contrast: one unweighted fit on all patients, two regimes.

        logit P(Y <= k | A, X, D)  =  alpha_k + beta*A + gamma'X + delta*D ,    k = 0 … 5

    with `X = C.STANDARDISATION_COVARIATES` and `D = 1[contraindicated]` (DECISION 10). `center` is
    omitted for [§14a]'s reason and [§15]'s permission covers [§14] whole. Linear terms only, no
    interaction — [§6]'s rule, and a stated assumption (§5.1).

    U1, U2, U7, U8 (and T10, T11 through `standardise.gcompute`). Records entries 3 and 4.
    """
    pol = _estimate(pop)
    _assert_beta_reportable(pol.conditional_log_odds)            # U8
    fit = pol.fit

    # ENTRY 3.
    audit.record(
        "model", "policy_fit", len(pop), "\n".join([
            f"[§14b]'s pooled proportional-odds model over ALL {len(pop)} patient(s), eligible and "
            f"contraindicated: {len(fit.columns)} design column(s), {len(pol.dropped)} dropped as "
            f"constant ({', '.join(pol.dropped) or 'none'}), {len(fit.alpha)} cutpoint(s) over "
            f"{len(fit.categories)} fitted mRS level(s) {fit.categories}.",
            f"Converged in {fit.iterations} iteration(s) on the {fit.converged_on!r} criterion; "
            f"{fit.rescales} trust-region rescale(s), {fit.halvings} halving(s); first step norm "
            f"{fit.first_step_norm:.6g}.",
            "UNWEIGHTED, one fit, no propensity model. DECISION 10 (PI, 2026-08-27): "
            "STANDARDISATION_COVARIATES plus a contraindication indicator as a MAIN EFFECT — it "
            "shifts the cutpoints and nothing else. `center` is not a covariate ([§15], inside §14 "
            "only).",
            f"abs(beta_treatment) {abs(pol.conditional_log_odds):.4f} against POLR_MAX_ABS_BETA = "
            f"{C.POLR_MAX_ABS_BETA:g} (U8, the ONLY bounded coefficient); largest abs(coefficient) "
            f"{float(np.max(np.abs(fit.beta))):.4f}, which may be the indicator's — a NUISANCE, "
            "unbounded [§9], not a result [§4.3], and given no interval.",
            "The treatment coefficient is identified from eligible patients alone: treatment is "
            "constant among the contraindicated, so a separated indicator cannot move it [§9.1].",
        ]), table=_fit_table(fit, pol.dropped))

    # ENTRY 4.
    audit.record(
        "model", "policy_contrast", pol.n_average, "\n".join([
            f"[§14b]'s g-computation: every patient predicted under the ACTIVE regime (bridging for "
            f"the {pol.n_eligible} eligible, EVT alone for the {pol.n_average - pol.n_eligible} "
            f"contraindicated) and under the COMPARATOR (EVT alone for everyone), averaged over "
            f"{pol.n_average} of {pol.n_fit} fitted, unweighted.",
            f"{pol.label}.",
            f"rd_k = share_eligible x eligible_rd_k, share_eligible = {pol.n_eligible}/"
            f"{pol.n_average} = {pol.share_eligible:.6f}: a contraindicated patient's per-row "
            "contrast is exactly zero, so the policy contrast is the eligible contrast diluted by "
            "the eligible share (U7, asserted from two independent computations). Read rd_k with "
            "the factor beside it.",
            "eligible_rd_k is the eligible-subpopulation contrast UNDER THIS FIT — a diagnostic "
            "factor, not Stage 12's estimate: the two are fitted on different populations and are "
            "never placed in one column [§14, §7.3].",
            f"conditional log-odds {pol.conditional_log_odds:+.6f}, exp() "
            f"{pol.conditional_odds_ratio:.6f} — {pol.measure}.",
            "The last rows of the table are the two namings and the dilution residual as COMPUTED "
            "DIFFERENCES, so the log shows the identities were checked [T11, U7].",
        ]), table=_contrast_table(pol))
    return pol


# --- inference [§10, §14] -------------------------------------------------------------------------

def _estimand_keys() -> tuple[str, ...]:
    """§3.3's thirty keys, counted from the config and never written as 6, 7 or 22::

        rd_0 … rd_5, mrs_0_2, mortality, active_0 … active_6, comparator_0 … comparator_6    22
        eligible_rd_0 … eligible_rd_5, share_eligible, beta                                    8

    No `contraindication` key: the indicator's coefficient is a nuisance and gets no interval (§9).
    """
    return (
        *(f"rd_{k}" for k in C.MRS_THRESHOLDS),
        "mrs_0_2", "mortality",
        *(f"{_ACTIVE}_{level}" for level in C.MRS_LEVELS),
        *(f"{_COMPARATOR}_{level}" for level in C.MRS_LEVELS),
        *(f"eligible_rd_{k}" for k in C.MRS_THRESHOLDS),
        "share_eligible", "beta",
    )


def _values(pol: Policy) -> dict[str, float]:
    """The thirty draws of one `Policy`, keyed as `_estimand_keys` names them."""
    values = {f"rd_{k}": v for k, v in pol.rd.items()}
    values["mrs_0_2"] = pol.mrs_0_2
    values["mortality"] = pol.mortality
    for level in C.MRS_LEVELS:
        values[f"{_ACTIVE}_{level}"] = pol.distribution[_ACTIVE][level]
        values[f"{_COMPARATOR}_{level}"] = pol.distribution[_COMPARATOR][level]
    values.update({f"eligible_rd_{k}": v for k, v in pol.eligible_rd.items()})
    values["share_eligible"] = pol.share_eligible
    values["beta"] = pol.conditional_log_odds
    return values


def _replicate(draw: pd.DataFrame, keys: tuple[str, ...]) -> bootstrap.Replicate:
    """One replicate: one fit off one drawn frame, three failure groups (§10.2).

    U1 and a `polr` `FitError` cost EVERY key under `degenerate_design` and `nonconvergence`; U8
    costs `beta` alone. U2, U4, U6 and U7 are `SchemaError` and are not caught: a contract break
    kills the bootstrap rather than being counted (Stage 12 §14).

    `max_abs_beta` is keyed `{"policy": max abs(coef)}` and INCLUDES the indicator, so the grid shows
    §9.1's tail beside a bound that does not apply to it.
    """
    try:
        pol = _estimate(draw)                                    # U1 -> degenerate_design; polr -> nonconvergence
    except model.FitError as failure:
        label = bootstrap.bucket(str(failure))
        return bootstrap.Replicate(
            values={}, failures={key: label for key in keys},
            n_alpha=None, polr_iterations=None, sum_w=None, n_in_model=None)

    values = _values(pol)
    failures: dict[str, str] = {}
    try:
        _assert_beta_reportable(pol.conditional_log_odds)        # U8 -> separation, `beta` only
    except model.FitError as failure:
        del values["beta"]
        failures["beta"] = bootstrap.bucket(str(failure))

    fit = pol.fit
    return bootstrap.Replicate(
        values=values, failures=failures, n_alpha=len(fit.alpha), polr_iterations=fit.iterations,
        sum_w=float(len(draw)), n_in_model=len(draw),
        max_abs_beta={"policy": float(np.max(np.abs(fit.beta)))})


def _replicates_table(draws: dict[str, bootstrap.Draws], diagnostics: bootstrap.Diagnostics,
                      collected: tuple[object, ...]) -> tuple[tuple[str, ...], ...]:
    """Entry 5's grid. Every block labels its own denominator; no cell contains a `|`."""
    attempted = max((d.n_attempted for d in draws.values()), default=0)
    live = len(diagnostics.sum_w)
    collapsed = sum(count for n_alpha, count in diagnostics.n_alpha.items()
                    if n_alpha < len(C.MRS_THRESHOLDS))
    share = draws["share_eligible"].draws
    share_one = sum(1 for r in collected
                    if r.n_alpha is not None and r.values.get("share_eligible") == 1.0)

    rows: list[tuple[str, str, str]] = []
    rows.append((f"attempted = {attempted}", "estimand keys", str(len(draws))))
    for name in sorted(set(C.FAILURE_BUCKETS.values())):
        total = sum(d.failures.get(name, 0) for d in draws.values())
        rows.append((f"attempted = {attempted}", f"key-failures bucketed {name!r}", str(total)))
    rows.append((f"reached the fit = {live}", "cutpoint counts len(alpha)", str(diagnostics.n_alpha)))
    rows.append((f"reached the fit = {live}", "replicates that lost a level (rd_4 == rd_5 there)",
                 f"{collapsed} ({100.0 * collapsed / live:.2f}%)" if live else "—"))
    rows.append((f"reached the fit = {live}", "polr iterations", str(diagnostics.polr_iterations)))
    rows.append((f"reached the fit = {live}", "Sw = n (unit weights)",
                 f"{diagnostics.sum_w.min():g} to {diagnostics.sum_w.max():g}" if live else "—"))
    rows.append((f"reached the fit = {live}", "share_eligible",
                 f"{share.min():.4f} to {share.max():.4f}" if len(share) else "—"))
    rows.append((f"reached the fit = {live}", "draws with NO contraindicated patient (share 1)",
                 str(share_one)))
    for key, spread in sorted(diagnostics.max_abs_beta.items()):
        rows.append((f"reached the {key} fit = {len(spread)}", "max_abs_coef, INCLUDING the indicator",
                     f"{spread.min():.4f} to {spread.max():.4f}; POLR_MAX_ABS_BETA = "
                     f"{C.POLR_MAX_ABS_BETA:g} binds beta ONLY, not the indicator — a separated "
                     "magnitude reflects POLR_TOL, not the data [§9.1]"))
    return replicate_rows([(denominator, [(quantity, value)]) for denominator, quantity, value in rows])


def inference(pop: pd.DataFrame, audit: Audit) -> bootstrap.Bootstrap:
    """The [§14] bootstrap: one `replicates` call, one arm per draw, thirty keys, no p-value.

    ALL FOUR STRATA are resampled and the contraindicated patients are members of the population:
    their covariates enter the average, their outcomes the fit, and their SHARE is a statistic of the
    resampled cohort — which is why `share_eligible` is a key with an interval (§10.1).

    Records entry 5. Calls `bootstrap.replicates`, `bucket`, `collect`, `intervals` and
    `standardise.diagnostics`; implements none of them. `bootstrap.diagnostics` is NOT called: it
    keeps `max_abs_beta` entries for `C.BINARY_OUTCOMES` only and would drop the `"policy"` key
    silently (§10.2).
    """
    keys = _estimand_keys()
    collected = bootstrap.replicates(
        pop, lambda draw: _replicate(draw, keys), C.N_BOOT, C.SEED, C.BOOT_STRATUM)
    draws = bootstrap.collect(collected, keys)
    diag = standardise.diagnostics(collected)
    limits = bootstrap.intervals(draws, lambda key: False)       # §3.3 — no p on ANY key

    # ENTRY 5.
    audit.record(
        "model", "policy_replicates", C.N_BOOT, "\n".join([
            f"[§14] inference for [§14b]: {C.N_BOOT} patient-level replicate(s) stratified by "
            f"{C.BOOT_STRATUM!r}, ALL FOUR STRATA. The contraindicated count varies by draw because "
            "the stratum is centre and not eligibility.",
            f"{len(keys)} estimand key(s); {len(limits)} interval(s) at level {C.CI_LEVEL:g} under "
            f"the {C.PERCENTILE_METHOD!r} percentile definition. NO KEY CARRIES A p-VALUE.",
            f"Seed {C.SEED}; the floor below which an estimand keeps its draws and gets no interval "
            f"is {C.ci_min_draws(C.CI_LEVEL)}.",
            "One arm off one drawn frame. Three failure groups: U1 and a non-convergence each cost "
            "every key; U8 costs `beta` alone. U2, U4, U6, U7 are contract breaks and kill the run.",
            "`share_eligible` and `eligible_rd_k` carry intervals because the eligible share is a "
            "statistic of the resampled cohort, not a constant of the design [§7.3, §10.1]. "
            "`mortality` carries its own draws and is not the reflection of `rd_5`'s interval.",
            "max_abs_coef in the grid INCLUDES the contraindication indicator, which the bound does "
            "not apply to; a separated magnitude there reflects POLR_TOL and not the data [§9.1].",
        ]), table=_replicates_table(draws, diag, collected))

    return bootstrap.Bootstrap(C.SEED, C.N_BOOT, draws, limits, diag)
