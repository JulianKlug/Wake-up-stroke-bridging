"""Stage 14 — [§16] outputs and guardrails.

One entry point runs Stages 1–13 in the canonical order and writes every table, figure and log under
``config.OUT``. Section references in brackets are to ``statistical_analysis_plan.md``. The
specification for this module is ``specs/stage14_outputs_and_guardrails.md``; nothing here is
invented outside it.

THIS MODULE COMPUTES NO ESTIMATE, INTERVAL OR P-VALUE. Every labelling rule [§16] imposes is derived
from a result object — a clause that depends on the data is a function of the record, never a
sentence in a template — so an edit cannot drop it without a test noticing.

    run(source) ─► Run (every landed result object, held not recomputed)
    write(run)  ─► out/tables/Txx_*.md + .csv, out/figures/Fxx_*.svg, out/logs/audit_*.md, Manifest
"""

from __future__ import annotations

import csv
import hashlib
import io
from collections.abc import Callable, Sequence
from dataclasses import dataclass, fields
from pathlib import Path
from typing import Final

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402  — after the backend is fixed
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import balance  # noqa: E402
import bootstrap  # noqa: E402
import cohort  # noqa: E402
import config as C  # noqa: E402
import data  # noqa: E402
import derive  # noqa: E402
import eligibility  # noqa: E402
import outcome  # noqa: E402
import policy  # noqa: E402
import propensity  # noqa: E402
import sensitivity  # noqa: E402
import standardise  # noqa: E402
from data import _fmt, _md_table  # noqa: E402

Rows = tuple[tuple[str, ...], ...]
Table = tuple[Rows, tuple[str, ...]]          # (header-first rows, clauses printed beneath)

_STEP_WRITTEN: Final[str] = "report_written"
_DASH: Final[str] = "—"
_ONSET_TO_GROIN: Final[str] = "onset_to_groin_min"
_CONTROL: Final[int] = min(C.TREATMENT_LABELS)
_TREATED: Final[int] = max(C.TREATMENT_LABELS)
_ARMS: Final[tuple[int, ...]] = (_CONTROL, _TREATED)


# --- the [§16] statements, as constants ------------------------------------------------------------

CONSTANT_SHIFT: Final[str] = (
    "The common odds ratio rests on a constant-shift (proportional-odds) assumption that is NOT "
    "tested [§16a, DECISION 5]; the six cumulative risk differences RD_k are the prescribed "
    "mitigation and are read beside it.")
DIRECTION: Final[str] = (
    "The six threshold differences RD_k do NOT all point in the same direction: a single odds ratio "
    "summarising them is not a null result [§16a].")
DESCRIPTIVE_CRUDE: Final[str] = (
    "DESCRIPTIVE. Raw pooled event rates by arm are not unadjusted estimates: treatment is nearly "
    "determined by centre, so a pooled crude comparison can differ from the within-centre one in "
    "direction [§16b, DECISION 5].")
THINNING: Final[str] = (
    "Replicates were removed by a bound on the estimated coefficient: this is selection on the "
    "ESTIMATE, not on the frame, and the percentile limits are quantiles of a truncated sampling "
    "distribution — systematically narrower than the untruncated one, not merely noisier [§16c].")
E_VALUE_HEURISTIC: Final[str] = (
    "The E-value is a heuristic on an approximated scale, not the bound the derivation licenses "
    "[§16d]; the direction of the approximation's error is conservative (sqrt(OR) < OR away from "
    "the null).")
MODEL_ASSISTED: Final[str] = "model-assisted"
DESCRIPTIVE_SAFETY: Final[str] = "descriptive"
HYPOTHESIS_GENERATING: Final[str] = "hypothesis-generating [§13]"
SINGLE_SHIFT: Final[str] = (
    "The interaction is a single proportional-odds shift delta, assumed constant across thresholds "
    "[§13]; the two levels' weighted cumulative distributions are the diagnostic.")
NON_COLLAPSIBLE: Final[str] = (
    "A weighted proportional-odds model is not collapsible: the two level-specific odds ratios need "
    "not bracket the primary and must not be read as though they do [§13].")
TRANSPORTED: Final[str] = (
    "Model-based transported results [§14, §16]: an estimate in a population other than [§7]'s ATO, "
    "under an untestable transport assumption.")
POPULATIONS: Final[str] = (
    "The [§7] primary, the [§14a] all-centre and the [§14b] policy estimates target DIFFERENT "
    "populations and are not a like-for-like comparison [§16].")
ATO_DESCRIPTION: Final[str] = (
    "The weighted columns describe the [§7] overlap (ATO) population — patients whose treatment was "
    "genuinely in equipoise at the participating centres — not the cohort.")
UNCORRECTED_PRIMARY: Final[str] = "The primary p-value is uncorrected [§13]."


# --- what the stage returns ----------------------------------------------------------------------

@dataclass(frozen=True)
class Run:
    """Every landed result object, held not recomputed. V4 if a stage was skipped."""

    source: data.Source
    audit: data.Audit
    df: pd.DataFrame                       # the classified frame, Stages 12 and 13's input
    cohort: pd.DataFrame                   # the [§3] cohort, Stages 6–11's input
    ps: propensity.Propensity
    balance: balance.Balance
    primary: outcome.Primary
    secondary: outcome.Secondary
    boot: bootstrap.Bootstrap
    multiplicity: sensitivity.Multiplicity
    e_value: sensitivity.EValue
    subgroups: sensitivity.Subgroups
    arm: sensitivity.Arm
    std_population: pd.DataFrame
    std: standardise.Standardisation
    support: standardise.Support
    std_support: standardise.Standardisation
    hier: standardise.Hierarchical
    std_boot: bootstrap.Bootstrap
    pol_population: pd.DataFrame
    pol: policy.Policy
    pol_boot: bootstrap.Bootstrap

    def __post_init__(self) -> None:
        missing = [f.name for f in fields(self) if getattr(self, f.name) is None]
        if missing:
            raise C.SchemaError(
                f"V4  Run is missing {missing}: a stage was skipped. A partial run would print T22 "
                "without its comparators [Stage 14 §3.1].")


@dataclass(frozen=True)
class Manifest:
    paths: tuple[Path, ...]        # relative to `out`, sorted
    sha256: dict[str, str]         # relative path -> digest of the written bytes


# --- the canonical call order, written once ------------------------------------------------------

def _stages_1_10(source: data.Source):
    df, audit = data.load(source)
    df = derive.derive(df, audit)
    df = eligibility.classify(df, audit)
    classified = df
    df = cohort.build(df, audit)
    ps = propensity.fit(df, audit)
    bal = balance.assess(df, ps, audit)
    est = outcome.primary(df, ps, audit)
    sec = outcome.secondary(df, ps, audit)
    boot = bootstrap.run(df, ps, est, sec, audit)
    return classified, df, audit, ps, bal, est, sec, boot


def _stage_11(df, ps, bal, est, sec, boot, audit):
    mult = sensitivity.multiplicity(sec, boot, audit)
    ev = sensitivity.e_value_primary(est, bal, boot, audit)
    sub = sensitivity.subgroups(df, ps, est, audit)
    arm = sensitivity.full_covariate(df, ps, audit)
    return mult, ev, sub, arm


def _stage_12(classified: pd.DataFrame, audit: data.Audit, with_inference: bool = True):
    """Stage 12 §19's order: `support` precedes the restricted `all_centre`, whose `over=` it is."""
    pop = standardise.population(classified, audit)
    std = standardise.all_centre(pop, audit)
    sup = standardise.support(pop, audit)
    ssup = standardise.all_centre(pop, audit, over=sup.inside)
    hier = standardise.hierarchical(pop, audit)
    boot = standardise.inference(pop, audit) if with_inference else None
    return pop, std, sup, ssup, hier, boot


def _stage_13(classified: pd.DataFrame, audit: data.Audit, with_inference: bool = True):
    pop = policy.population(classified, audit)
    pol = policy.contrast(pop, audit)
    boot = policy.inference(pop, audit) if with_inference else None
    return pop, pol, boot


def run(source: data.Source = data.WORKBOOK) -> Run:
    """Stages 1–13 over `source`, in the canonical order. Writes nothing."""
    classified, df, audit, ps, bal, est, sec, boot = _stages_1_10(source)
    mult, ev, sub, arm = _stage_11(df, ps, bal, est, sec, boot, audit)
    std_pop, std, sup, ssup, hier, std_boot = _stage_12(classified, audit)
    pol_pop, pol, pol_boot = _stage_13(classified, audit)
    return Run(source, audit, classified, df, ps, bal, est, sec, boot, mult, ev, sub, arm,
               std_pop, std, sup, ssup, hier, std_boot, pol_pop, pol, pol_boot)


# --- the data-derived clauses --------------------------------------------------------------------

def _direction_clause(rd: dict[int, float]) -> str:
    """Present iff the six RD_k do not share a sign — computed, never asserted [§16a]."""
    return DIRECTION if len({np.sign(v) for v in rd.values()}) > 1 else ""


def _thinning_clause(draws: dict[str, bootstrap.Draws]) -> str:
    """Present iff a separation bound removed a replicate; names the surviving counts [§16c]."""
    thinned = {key: d for key, d in draws.items() if d.failures.get("separation", 0) > 0}
    if not thinned:
        return ""
    counts = "; ".join(f"{key}: {len(d.draws)} of {d.n_attempted} survive" for key, d in thinned.items())
    return f"{THINNING} Surviving counts — {counts}."


def _exceedance_clause(bal: balance.Balance) -> str:
    """Names every covariate whose residual |SMD| reaches the threshold, with its magnitude [§9]."""
    over = bal.unbalanced()
    if not over:
        return f"No covariate exceeds |SMD| {_fmt(C.SMD_THRESHOLD)} after weighting [§9]."
    by_name = {row.covariate: row.weighted for row in bal.covariates}
    named = ", ".join(f"{name} ({by_name[name]:+.3f})" for name in over)
    return f"Residual imbalance after weighting, |SMD| >= {_fmt(C.SMD_THRESHOLD)} [§9]: {named}."


def _e_value_range(cumulative: dict[int, dict[int, float]]) -> str:
    """Names the thresholds whose control-arm baseline risk lies outside the conversion's licence."""
    outside = [k for k in C.MRS_THRESHOLDS
               if min(cumulative[k][_CONTROL], 1.0 - cumulative[k][_CONTROL]) < C.E_VALUE_PREVALENCE_FLOOR]
    if not outside:
        return (f"Every threshold's control-arm baseline risk lies within the conversion's documented "
                f"licence (prevalence above {_fmt(C.E_VALUE_PREVALENCE_FLOOR)}) [§16d].")
    return (f"Thresholds whose control-arm baseline risk lies outside the conversion's documented "
            f"licence (prevalence within {_fmt(C.E_VALUE_PREVALENCE_FLOOR)} of 0 or 1) [§16d]: "
            + ", ".join(f"RD_{k}" for k in outside) + ".")


def _populations_clause(run: Run) -> str:
    """The different-populations statement, printed as numbers [Stage 12 §17.2, Stage 13 §14]."""
    n_primary = int(run.primary.in_estimate.sum())
    centres = len(cohort.treating_centres(run.cohort))
    return (f"{POPULATIONS} Primary [§7]: {n_primary} patients at {centres} centres. "
            f"[§14a]: {run.std.n_average} patients at {run.std_population['center'].nunique()} centres. "
            f"[§14b]: {run.pol.n_average} patients at {run.pol_population['center'].nunique()} centres, "
            f"of whom {run.pol.n_average - run.pol.n_eligible} contraindicated.")


# --- cell helpers --------------------------------------------------------------------------------

def _ci(intervals: dict[str, bootstrap.Interval], key: str,
        transform: Callable[[float], float] = lambda x: x) -> tuple[str, str, str, str]:
    """(lo, hi, p, n_draws) as cells; dashes where Stage 10 withheld the interval."""
    i = intervals.get(key)
    if i is None:
        return (_DASH, _DASH, _DASH, "0")
    return (_fmt(transform(i.lo)), _fmt(transform(i.hi)), _DASH if i.p is None else _fmt(i.p),
            str(i.n_draws))


def _count(mask: pd.Series) -> str:
    return str(int(mask.sum()))


def _arm_mask(df: pd.DataFrame, arm: int) -> pd.Series:
    return df[C.TREATMENT] == arm


def _weighted_mean(x: pd.Series, w: pd.Series) -> float:
    ok = x.notna() & w.notna()
    if not ok.any() or float(w[ok].sum()) == 0.0:
        return float("nan")
    return float((x[ok] * w[ok]).sum() / w[ok].sum())


def _pmf(cumulative: dict[int, dict[int, float]], arm: int) -> list[float]:
    """P(Y=j) for every MRS level, from the cumulative table; the last level is the complement."""
    out, previous = [], 0.0
    for k in C.MRS_THRESHOLDS:
        out.append(cumulative[k][arm] - previous)
        previous = cumulative[k][arm]
    out.append(1.0 - previous)
    return out


def _interval_note(boot: bootstrap.Bootstrap) -> str:
    return (f"Intervals: {_fmt(C.CI_LEVEL)} percentile ({C.PERCENTILE_METHOD}), {boot.n_boot} "
            f"replicates, seed {boot.seed}; the surviving draw count is printed beside each [§10].")


# --- tables --------------------------------------------------------------------------------------

def _t01(run: Run) -> Table:
    flow = run.audit.entry("cohort", "cohort_flow").table
    width = len(flow[0])
    cohort_df = run.cohort

    def denominator(label: str, mask: pd.Series) -> tuple[str, ...]:
        row = (f"[§11] {label}", str(int(cohort_df.loc[mask, "center"].nunique())), _count(mask))
        return row + (_DASH,) * (width - len(row))

    rows = list(flow)
    rows.append(denominator("covariate-complete (in_model)", run.ps.in_model))
    rows.append(denominator("primary outcome present (in_estimate)", run.primary.in_estimate))
    for key in C.BINARY_OUTCOMES:
        est = run.secondary.estimates.get(key)
        if est is not None:
            rows.append(denominator(f"{key} in_estimate", est.in_estimate))
    return tuple(rows), ("Every estimate carries its own denominator [§11]; the three nested "
                         "populations are cohort, covariate-complete, outcome-present.",)


def _t02(run: Run) -> Table:
    df = run.cohort[run.ps.in_model]
    w = run.ps.w[run.ps.in_model]
    header = ("covariate", "n", *(f"unweighted {C.TREATMENT_LABELS[a]}" for a in _ARMS),
              *(f"weighted {C.TREATMENT_LABELS[a]}" for a in _ARMS))
    rows = []
    for name in C.BALANCE_SET:
        for label, x in balance.levels(df, name):
            x = x.astype(float)
            rows.append((label, _count(x.notna()),
                         *(_fmt(x[_arm_mask(df, a)].mean()) for a in _ARMS),
                         *(_fmt(_weighted_mean(x[_arm_mask(df, a)], w[_arm_mask(df, a)])) for a in _ARMS)))
    rows.append(("effective sample size", _count(run.ps.in_model), *(_DASH,) * 2,
                 *(_fmt(run.ps.ess[a]) for a in _ARMS)))
    return (header, *rows), (ATO_DESCRIPTION,)


def _t03(run: Run) -> Table:
    df = run.cohort

    def cell(frame: pd.DataFrame, key: str, arm: int) -> str:
        y = frame.loc[_arm_mask(frame, arm), key].dropna()
        return f"{int(y.sum())} of {len(y)} ({_fmt(y.mean() if len(y) else float('nan'))})" \
            if len(y) else "0 of 0"

    header = ("outcome", "rule", *(f"pooled {C.TREATMENT_LABELS[a]}" for a in _ARMS),
              *(f"{c} {C.TREATMENT_LABELS[a]}" for c in C.CENTER_ORDER for a in _ARMS))
    rows = []
    for key in C.BINARY_OUTCOMES:
        o = C.OUTCOMES[key]
        rows.append((f"{DESCRIPTIVE_SAFETY}: {o.label}" if o.family == "safety" else o.label,
                     o.rule or "read directly", *(cell(df, key, a) for a in _ARMS),
                     *(cell(df[df["center"] == c], key, a) for c in C.CENTER_ORDER for a in _ARMS)))
    return (header, *rows), (DESCRIPTIVE_CRUDE,)


def _t04(run: Run) -> Table:
    df = run.cohort
    header = ("arm", "n", "median", "q1", "q3", "min", "max")
    rows = []
    for a in _ARMS:
        x = df.loc[_arm_mask(df, a), _ONSET_TO_GROIN].dropna().astype(float)
        q = x.quantile([0.25, 0.5, 0.75]) if len(x) else pd.Series([np.nan] * 3, index=[0.25, 0.5, 0.75])
        rows.append((C.TREATMENT_LABELS[a], str(len(x)), _fmt(q[0.5]), _fmt(q[0.25]), _fmt(q[0.75]),
                     _fmt(x.min() if len(x) else np.nan), _fmt(x.max() if len(x) else np.nan)))
    return (header, *rows), ("Onset-to-groin time is post-exposure and is reported by arm, never "
                             "modelled [§12]; any change in the estimate when it is added is "
                             "descriptive.",)


def _smd_rows(rows: Sequence[balance.CovariateBalance]) -> Rows:
    header = ("covariate", "role", "n", "pooled sd", "SMD unweighted", "SMD weighted", "exceeds")
    return (header, *((r.covariate, r.role, str(r.n), _fmt(r.sd), _fmt(r.unweighted), _fmt(r.weighted),
                       "yes" if np.isfinite(r.weighted) and abs(r.weighted) >= C.SMD_THRESHOLD else "no")
                      for r in rows))


def _t05(run: Run) -> Table:
    return _smd_rows(run.balance.covariates), (
        f"Balance is judged against the full [§6] confounder set at |SMD| < {_fmt(C.SMD_THRESHOLD)}; "
        f"roles are to {run.balance.spec.label}.", _exceedance_clause(run.balance))


def _overlap_rows(centres: Sequence[balance.CentreOverlap]) -> Rows:
    header = ("centre", "in cohort", "weighted", *(f"n {C.TREATMENT_LABELS[a]}" for a in _ARMS),
              *(f"e range {C.TREATMENT_LABELS[a]}" for a in _ARMS),
              *(f"ESS {C.TREATMENT_LABELS[a]}" for a in _ARMS), "max weight", "weight share",
              "worst SMD", "status")
    rows = []
    for c in centres:
        rows.append((c.centre, str(c.in_cohort), str(c.weighted), *(str(c.n.get(a, 0)) for a in _ARMS),
                     *(f"{_fmt(c.e_range[a][0])} to {_fmt(c.e_range[a][1])}" if a in c.e_range else _DASH
                       for a in _ARMS),
                     *(_fmt(c.ess.get(a, np.nan)) for a in _ARMS), _fmt(c.max_weight),
                     _fmt(c.weight_share), _fmt(c.worst), c.status))
    return (header, *rows)


def _t06(run: Run) -> Table:
    structural = [c.centre for c in run.balance.centres if c.status != run.balance.pooled.status]
    note = ("Centres with structural non-positivity or no weighted contrast are reported as such and "
            f"given no overlap plot [§9]: {', '.join(structural) if structural else 'none'}.")
    return _overlap_rows((*run.balance.centres, run.balance.pooled)), (note,)


def _t07(run: Run) -> Table:
    n = _count(run.primary.in_estimate)
    header = ("quantity", "estimate", "ci lo", "ci hi", "p", "draws", "denominator")
    rows = [("common odds ratio exp(beta)", _fmt(run.primary.odds_ratio),
             *_ci(run.boot.intervals, "beta", np.exp), n)]
    rows.extend((f"RD_{k}  P(mRS<=k) difference", _fmt(run.primary.rd[k]),
                 *_ci(run.boot.intervals, f"rd_{k}"), n) for k in C.MRS_THRESHOLDS)
    clauses = (CONSTANT_SHIFT, _direction_clause(run.primary.rd), _exceedance_clause(run.balance),
               _interval_note(run.boot),
               "The p-value tests the proportional-odds treatment coefficient, not a risk-difference "
               "scale quantity [§10]; RD_k carry no p-value [§8].")
    return (header, *rows), tuple(c for c in clauses if c)


def _distribution_rows(cumulative: dict[int, dict[int, float]], arms: Sequence[tuple[str, int]]) -> Rows:
    header = ("mRS", *(f"P(Y=j) {label}" for label, _ in arms), *(f"P(Y<=k) {label}" for label, _ in arms))
    pmf = {arm: _pmf(cumulative, arm) for _, arm in arms}
    rows = []
    for j in C.MRS_LEVELS:
        rows.append((str(j), *(_fmt(pmf[arm][j]) for _, arm in arms),
                     *(_fmt(cumulative[j][arm]) if j in cumulative else "1" for _, arm in arms)))
    return (header, *rows)


def _t08(run: Run) -> Table:
    arms = tuple((C.TREATMENT_LABELS[a], a) for a in _ARMS)
    return _distribution_rows(run.primary.cumulative, arms), (
        "Weighted per-category mRS distribution in the [§7] population; the six RD_k are differences "
        "of the cumulative columns.",)


def _t09(run: Run) -> Table:
    ev = run.e_value
    header = ("quantity", "value")
    rows = (("measure", ev.measure), ("odds ratio", _fmt(ev.odds_ratio)),
            ("risk ratio (approximated)", _fmt(ev.risk_ratio)), ("E-value, point", _fmt(ev.e_point)),
            ("interval spans the null", "yes" if ev.spans_null else "no"),
            ("limit nearest the null (beta scale)", _fmt(ev.limit)),
            ("E-value, limit", _fmt(ev.e_limit)), ("draws behind the limit", str(ev.n_draws)))
    return (header, *rows), (ev.approximation, E_VALUE_HEURISTIC,
                             _e_value_range(run.primary.cumulative), _exceedance_clause(run.balance))


def _binary_rows(run: Run, family: str) -> Table:
    corr = run.multiplicity.families[family]
    header = ("outcome", "label", "denominator", "RD", "ci lo", "ci hi", "p raw", "p adjusted",
              "odds ratio", "ci lo", "ci hi", f"augmented ({MODEL_ASSISTED})", "ci lo", "ci hi",
              "augmentation path", "m_a(X) covariates")
    rows = []
    for est in run.secondary.by_family()[family]:
        key = est.outcome
        label = DESCRIPTIVE_SAFETY if family == "safety" else "secondary"
        lo, hi, p, _ = _ci(run.boot.intervals, f"{key}.rd")
        or_lo, or_hi, _, _ = _ci(run.boot.intervals, f"{key}.odds_ratio")
        odds = _fmt(est.odds_ratio) + (" (corrected)" if est.or_corrected else "")
        aug = (_fmt(est.augmented), *_ci(run.boot.intervals, f"{key}.augmented")[:2]) \
            if est.augmented_path != "unaugmented" else (_DASH, _DASH, _DASH)
        covariates = ", ".join(est.covariates) + (" (reduced [§8])" if est.reduced else "") \
            if est.covariates else _DASH
        rows.append((C.OUTCOMES[key].label, label, _count(est.in_estimate), _fmt(est.rd), lo, hi, p,
                     _fmt(corr.adjusted.get(key)), odds, or_lo, or_hi, *aug, est.augmented_path,
                     covariates))
    for key, why in corr.absent:
        rows.append((C.OUTCOMES[key].label, "no p-value", *(_DASH,) * 13, why))
    for key, message in run.secondary.failures.items():
        rows.append((C.OUTCOMES[key].label, "not estimable", *(_DASH,) * 13, message))
    clauses = [f"Benjamini-Hochberg within the {family} family: m declared {corr.m_declared}, "
               f"m used {corr.m_used} [§13]. {UNCORRECTED_PRIMARY}", _interval_note(run.boot)]
    if family == "safety":
        clauses.insert(0, f"{DESCRIPTIVE_SAFETY.upper()} ONLY: safety outcomes resting on few events "
                          "are estimation, not testing, and every row carries the label [§10, §13].")
    else:
        clauses.insert(0, f"The augmented estimate is {MODEL_ASSISTED}, not doubly robust, and is "
                          "always read beside the unaugmented weighted RD [§8].")
    return (header, *rows), tuple(clauses)


def _t10(run: Run) -> Table:
    return _binary_rows(run, "secondary")


def _t11(run: Run) -> Table:
    return _binary_rows(run, "safety")


def _t12(run: Run) -> Table:
    arm = run.arm
    worst, undefined = arm.balance.worst()
    header = ("row", "specification", "in_model", *(f"ESS {C.TREATMENT_LABELS[a]}" for a in _ARMS),
              "worst residual SMD", "exp(beta)", "ci lo", "ci hi", "draws", "status")
    primary_worst, _ = run.balance.worst()
    rows = [
        ("primary [§7]", run.ps.spec.label, _count(run.ps.in_model), *(_fmt(run.ps.ess[a]) for a in _ARMS),
         _fmt(primary_worst.weighted if primary_worst else np.nan), _fmt(run.primary.odds_ratio),
         *_ci(run.boot.intervals, "beta", np.exp)[:2], _ci(run.boot.intervals, "beta")[3], "reported"),
        ("full-covariate arm [§13]", arm.ps.spec.label, _count(arm.ps.in_model),
         *(_fmt(arm.ps.ess[a]) for a in _ARMS), _fmt(worst.weighted if worst else np.nan),
         _fmt(arm.estimate.odds_ratio), *_ci(arm.intervals, "beta", np.exp)[:2],
         _ci(arm.intervals, "beta")[3], "reported; no p-value [§13]"),
        ("other [§13] sensitivity rows", _DASH, *(_DASH,) * 8, "deferred (DECISION 4)"),
    ]
    clauses = (f"The arm differs from the primary only in specification: "
               f"{'yes' if arm.differs_only_in_specification else 'no'}. Undefined SMD rows: "
               f"{', '.join(undefined) if undefined else 'none'}.",
               "Agreement with the primary must not be read as reassurance about unmeasured "
               "confounding [§13].")
    return (header, *rows), clauses


def _t13(run: Run) -> Table:
    sub = run.subgroups
    header = ("subgroup", "n", "cells S,A", "OR at S=0", "ci lo", "ci hi", "OR at S=1", "ci lo",
              "ci hi", "interaction OR exp(gamma)", "ci lo", "ci hi", "p", "draws", "label")
    rows = []
    for key, est in sub.estimates.items():
        cells = "; ".join(f"{k}={v}" for k, v in est.n_by_level_arm.items())
        g_lo, g_hi, g_p, g_n = _ci(sub.intervals, f"{key}.gamma", np.exp)
        g_p = _ci(sub.intervals, f"{key}.gamma")[2]
        rows.append((C.SUBGROUPS[key], str(est.n), cells, _fmt(est.or_level[0]),
                     *_ci(sub.intervals, f"{key}.beta", np.exp)[:2], _fmt(est.or_level[1]),
                     *_ci(sub.intervals, f"{key}.beta_plus_gamma", np.exp)[:2], _fmt(est.or_ratio),
                     g_lo, g_hi, g_p, g_n,
                     HYPOTHESIS_GENERATING if est.hypothesis_generating else "prespecified"))
    n_tests = next(iter(sub.estimates.values())).n_interaction_tests if sub.estimates else 0
    clauses = (HYPOTHESIS_GENERATING, SINGLE_SHIFT, NON_COLLAPSIBLE,
               f"{n_tests} interaction test(s) performed; their p-values are an uncorrected "
               "multiplicity the plan does not address [§13].",
               _thinning_clause(sub.draws),
               f"Intervals: {_fmt(C.CI_LEVEL)} percentile ({C.PERCENTILE_METHOD}), {sub.n_boot} "
               f"replicates, seed {sub.seed}.")
    return (header, *rows), tuple(c for c in clauses if c)


def _t14(run: Run) -> Table:
    header = ("subgroup", "level", "threshold k", *(f"P(Y<=k) {C.TREATMENT_LABELS[a]}" for a in _ARMS))
    rows = []
    for key, est in run.subgroups.estimates.items():
        for level, by_k in est.cumulative.items():
            for k in C.MRS_THRESHOLDS:
                rows.append((C.SUBGROUPS[key], str(level), str(k), *(_fmt(by_k[k][a]) for a in _ARMS)))
    return (header, *rows), (SINGLE_SHIFT,)


def _std_block(label: str, std: standardise.Standardisation, intervals: dict[str, bootstrap.Interval],
               prefix: str) -> list[tuple[str, ...]]:
    rows = [(label, f"n_average", str(std.n_average), _DASH, _DASH, _DASH),
            (label, "n_fit", str(std.n_fit), _DASH, _DASH, _DASH)]
    for arm in _ARMS:
        for j in C.MRS_LEVELS:
            rows.append((label, f"P(Y={j}) {C.TREATMENT_LABELS[arm]}", _fmt(std.distribution[arm][j]),
                         *_ci(intervals, f"{prefix}dist{arm}_{j}")[:2], _DASH))
    for k in C.MRS_THRESHOLDS:
        rows.append((label, f"RD_{k}", _fmt(std.rd[k]), *_ci(intervals, f"{prefix}rd_{k}")[:2],
                     _ci(intervals, f"{prefix}rd_{k}")[3]))
    rows.append((label, "mRS 0-2 difference", _fmt(std.mrs_0_2), *_ci(intervals, f"{prefix}mrs_0_2")[:2],
                 _ci(intervals, f"{prefix}mrs_0_2")[3]))
    rows.append((label, "mortality difference", _fmt(std.mortality),
                 *_ci(intervals, f"{prefix}mortality")[:2], _ci(intervals, f"{prefix}mortality")[3]))
    if prefix != "support.":
        beta_key = f"{prefix}beta"
        rows.append((label, f"conditional odds ratio — {std.measure}", _fmt(std.conditional_odds_ratio),
                     *_ci(intervals, beta_key, np.exp)[:2], _ci(intervals, beta_key)[3]))
    return rows


def _collapse_count(draws: dict[str, bootstrap.Draws], prefix: str) -> str:
    a, b = draws.get(f"{prefix}rd_4"), draws.get(f"{prefix}rd_5")
    if a is None or b is None or len(a.draws) != len(b.draws):
        return _DASH
    return f"{int(np.sum(a.draws == b.draws))} of {len(a.draws)}"


def _floor_rate(draws: bootstrap.Draws | None) -> str:
    """Replicates whose `hier.sigma` came back AT `POLR_RI_SIGMA_FLOOR`, over the surviving draws.

    An EXACT equality and that is safe rather than fragile: `model._ri_fit` snaps `sigma` to the
    constant at the boundary precisely so this comparison means "collapsed" (model.py:1471-1477), and
    `standardise._replicates_table` counts it the same way. Computed from the retained draws for
    `_collapse_count`'s reason — the rate is a property of the interval's own order statistics, so a
    second field carrying it could disagree with them.
    """
    if draws is None:
        return _DASH
    return f"{int(np.sum(draws.draws == C.POLR_RI_SIGMA_FLOOR))} of {len(draws.draws)}"


def _t15(run: Run) -> Table:
    header = ("population", "quantity", "estimate", "ci lo", "ci hi", "draws")
    blocks = ((run.std, ""), (run.std_support, "support."), (run.hier.standardisation, "hier."))
    rows = []
    for std, prefix in blocks:
        rows.extend(_std_block(std.population, std, run.std_boot.intervals, prefix))
    rows.append((run.hier.standardisation.population, "sigma (centre intercept sd)", _fmt(run.hier.sigma),
                 *_ci(run.std_boot.intervals, "hier.sigma")[:2], _ci(run.std_boot.intervals, "hier.sigma")[3]))
    for centre in C.CENTER_ORDER:
        rows.append((run.hier.standardisation.population, f"posterior sd of intercept, {centre}",
                     _fmt(run.hier.posterior_sd.get(centre)), _DASH, _DASH, _DASH))
    clauses = (TRANSPORTED,
               "exp(beta) is a CONDITIONAL odds ratio and a model parameter, never the standardised "
               "marginal effect [§14a]; `beta` is conditional on X, `hier.beta` on X and the centre "
               "intercept. The treated-support arm restricts the averaging population, not the fit, "
               "and has no odds ratio of its own.",
               f"Draws in which RD_5 equals RD_4 exactly (a lost mRS level): all eligible "
               f"{_collapse_count(run.std_boot.draws, '')}; support {_collapse_count(run.std_boot.draws, 'support.')}; "
               f"random intercept {_collapse_count(run.std_boot.draws, 'hier.')}.",
               f"sigma is bounded below at POLR_RI_SIGMA_FLOOR = {C.POLR_RI_SIGMA_FLOOR:g}, and "
               f"[§14a] names sigma^2_C = 0 a legitimate answer, so a limit sitting AT the floor "
               f"reports a variance that collapsed and is not a value for the between-centre SD — "
               f"an interval so limited is not a range for centre heterogeneity [Stage 12 §12.5]. "
               f"Point estimate at the floor: {'yes' if run.hier.at_floor else 'no'}; "
               f"replicates at the floor: {_floor_rate(run.std_boot.draws.get('hier.sigma'))}.",
               _interval_note(run.std_boot))
    return (header, *rows), clauses


def _t16(run: Run) -> Table:
    sup = run.support
    header = ("quantity", "value")
    rows = [(f"treated-arm box, {name}", f"{_fmt(lo)} to {_fmt(hi)}") for name, (lo, hi) in sup.box.items()]
    rows.extend((f"never-IVT patients outside the box, {c}", str(sup.outside_by_centre.get(c, 0)))
                for c in C.CENTER_ORDER)
    rows.append(("never-IVT patients", str(sup.n_never_ivt)))
    rows.append(("never-IVT patients outside the box", str(sup.n_never_ivt_outside)))
    return (header, *rows), ("The [§14a] transport assumption is not testable from these data; the "
                             "support check sizes the extrapolation and does not test it [§14a].",)


def _t17(run: Run) -> Table:
    return _smd_rows(run.support.baseline), (
        "IVT-treated patients against never-IVT-centre patients [§14a]; `center` rows are the "
        "grouping variable, not imbalance.",)


def _policy_role(name: str) -> str:
    if name == C.TREATMENT:
        return "exposure"
    if name == C.CONTRAINDICATED:
        return "contraindication indicator — a nuisance parameter, not a result [§14b]"
    return "covariate"


def _t18(run: Run) -> Table:
    pol, iv = run.pol, run.pol_boot.intervals
    header = ("row", "estimate", "ci lo", "ci hi", "draws", "share_eligible", "ci lo", "ci hi",
              "factor: eligible contrast under this model", "ci lo", "ci hi")
    blank = (_DASH,) * 6
    rows = []
    for k in C.MRS_THRESHOLDS:
        rows.append((f"rd_{k}", _fmt(pol.rd[k]), *_ci(iv, f"rd_{k}")[:2], _ci(iv, f"rd_{k}")[3],
                     _fmt(pol.share_eligible), *_ci(iv, "share_eligible")[:2],
                     _fmt(pol.eligible_rd[k]), *_ci(iv, f"eligible_rd_{k}")[:2]))
    rows.append(("mRS 0-2 difference", _fmt(pol.mrs_0_2), *_ci(iv, "mrs_0_2")[:2], _ci(iv, "mrs_0_2")[3],
                 *blank))
    rows.append(("mortality difference", _fmt(pol.mortality), *_ci(iv, "mortality")[:2],
                 _ci(iv, "mortality")[3], *blank))
    rows.append(("n_average", str(pol.n_average), *(_DASH,) * 3, *blank))
    rows.append(("n_eligible", str(pol.n_eligible), *(_DASH,) * 3, *blank))
    # The fit, coefficient by coefficient with its role: the indicator is a nuisance and has no interval.
    for name, value, _, role in data.coefficient_rows(pol.fit.columns, pol.fit.beta, pol.dropped,
                                                      _policy_role)[1:]:
        rows.append((f"coefficient {name} — {role}", value, *(_DASH,) * 3, *blank))
    rows.append((f"conditional odds ratio — {pol.measure}", _fmt(pol.conditional_odds_ratio),
                 *_ci(iv, "beta", np.exp)[:2], _ci(iv, "beta")[3], *blank))
    clauses = [pol.label, TRANSPORTED,
               "rd_k = share_eligible x eligible_rd_k: `eligible_rd_k` is the FACTOR that decomposes "
               "the policy contrast and is a diagnostic; [§14a]'s RD_k is THE eligible-population "
               "estimate and is fitted on a different population. No p-value [§14b].",
               f"Draws in which rd_5 equals rd_4 exactly: {_collapse_count(run.pol_boot.draws, '')}.",
               _interval_note(run.pol_boot)]
    if pol.share_eligible == 1.0:
        clauses.append(policy.NO_CONTRAINDICATED)
    return (header, *rows), tuple(clauses)


def _t19(run: Run) -> Table:
    buckets = tuple(sorted(set(C.FAILURE_BUCKETS.values())))
    header = ("block", "denominator", "estimand", "attempted", "draws", *buckets)
    blocks = (("primary and secondary [§10]", run.boot.draws), ("subgroups [§13]", run.subgroups.draws),
              ("sensitivity arm [§13]", run.arm.draws), ("[§14a]", run.std_boot.draws),
              ("[§14b]", run.pol_boot.draws))
    rows = []
    for label, draws in blocks:
        for key, d in draws.items():
            rows.append((label, f"attempted = {d.n_attempted}", key, str(d.n_attempted), str(len(d.draws)),
                         *(str(d.failures.get(b, 0)) for b in buckets)))
    clauses = ("Every row names its denominator; `draws` is a column and is never derived by "
               "subtraction [§10].", _thinning_clause(run.subgroups.draws))
    return (header, *rows), tuple(c for c in clauses if c)


def _t20(run: Run) -> Table:
    buckets = tuple(sorted(set(C.FAILURE_BUCKETS.values())))
    all_draws = {**run.boot.draws, **{f"sub.{k}": v for k, v in run.subgroups.draws.items()},
                 **{f"arm.{k}": v for k, v in run.arm.draws.items()},
                 **{f"std.{k}": v for k, v in run.std_boot.draws.items()},
                 **{f"pol.{k}": v for k, v in run.pol_boot.draws.items()}}
    header = ("key", "value")
    rows = [("source", run.source.label), ("DATA_SHA256", C.DATA_SHA256 if run.source.sha256 else "not pinned"),
            ("N_RECORDS_EXPECTED", str(C.N_RECORDS_EXPECTED)), ("SEED", str(C.SEED)),
            ("n_boot primary/secondary", str(run.boot.n_boot)), ("n_boot subgroups", str(run.subgroups.n_boot)),
            ("n_boot sensitivity arm", str(run.arm.n_boot)), ("n_boot [§14a]", str(run.std_boot.n_boot)),
            ("n_boot [§14b]", str(run.pol_boot.n_boot)), ("CI_LEVEL", _fmt(C.CI_LEVEL)),
            ("PERCENTILE_METHOD", C.PERCENTILE_METHOD)]
    rows.extend((f"key-failures bucketed {b}", str(sum(d.failures.get(b, 0) for d in all_draws.values())))
                for b in buckets)
    rows.append(("audit entries before this report", str(sum(
        1 for e in run.audit.entries if not (e.kind == "provenance" and e.step == _STEP_WRITTEN)))))
    rows.append(("outputs", str(len(C.OUTPUT_IDS))))
    return (header, *rows), ()


def _t21(run: Run) -> Table:
    header = ("item", "outputs")
    rows = []
    for item, ids in C.CHECKLIST:
        if not ids:
            raise C.SchemaError(f"V3  checklist item {item!r} maps to no output [Stage 14 §9].")
        rows.append((item, ", ".join(f"{i} {C.OUTPUT_IDS[i]}" for i in ids)))
    return (header, *rows), ("STROBE with the RECORD extension [§16]; each item names the output that "
                             "discharges it.",)


def _t22(run: Run) -> Table:
    columns = (("primary [§7] ATO", run.primary.rd, run.boot.intervals, "",
                {2: run.primary.rd[2], 5: -run.primary.rd[5]}),
               ("[§14a] all eligible", run.std.rd, run.std_boot.intervals, "",
                {2: run.std.mrs_0_2, 5: run.std.mortality}),
               (f"[§14b] policy", run.pol.rd, run.pol_boot.intervals, "",
                {2: run.pol.mrs_0_2, 5: run.pol.mortality}))
    header = ("quantity", *(f"{label}" for label, *_ in columns), *(f"{label} ci" for label, *_ in columns))
    n_row = ("population", _count(run.primary.in_estimate), str(run.std.n_average), str(run.pol.n_average),
             *(_DASH,) * 3)

    def ci_cell(intervals, key: str) -> str:
        lo, hi, _, n = _ci(intervals, key)
        return f"{lo} to {hi} ({n} draws)"

    rows = [n_row]
    for k in C.MRS_THRESHOLDS:
        rows.append((f"RD_{k}", *(_fmt(rd[k]) for _, rd, *_ in columns),
                     *(ci_cell(iv, f"rd_{k}") for _, _, iv, *_ in columns)))
    rows.append(("mRS 0-2 difference", *(_fmt(named[2]) for *_, named in columns),
                 *(ci_cell(iv, "rd_2" if i == 0 else "mrs_0_2") for i, (_, _, iv, *_) in enumerate(columns))))
    rows.append(("mortality difference", *(_fmt(named[5]) for *_, named in columns),
                 *(ci_cell(iv, "rd_5" if i == 0 else "mortality") for i, (_, _, iv, *_) in enumerate(columns))))
    sup = run.support
    grouping = [r for r in sup.baseline if not r.covariate.startswith("center")]
    largest = max((abs(r.weighted) for r in grouping if np.isfinite(r.weighted)), default=np.nan)
    clauses = (_populations_clause(run), TRANSPORTED,
               f"[§14a] support check: {sup.n_never_ivt} never-IVT patients, {sup.n_never_ivt_outside} "
               f"outside the treated-arm box; largest non-grouping baseline |SMD| {_fmt(largest)} "
               f"against a threshold of {_fmt(C.SMD_THRESHOLD)} [§14a]. The transport assumption is "
               "not testable from these data.",
               "The primary's mortality difference is -RD_5 by definition; no odds ratio is shown for "
               "the [§14] columns, whose exp(beta) is conditional [§14a, §14b].")
    return (header, *rows), clauses


_TABLES: Final[dict[str, Callable[[Run], Table]]] = {
    "T01": _t01, "T02": _t02, "T03": _t03, "T04": _t04, "T05": _t05, "T06": _t06, "T07": _t07,
    "T08": _t08, "T09": _t09, "T10": _t10, "T11": _t11, "T12": _t12, "T13": _t13, "T14": _t14,
    "T15": _t15, "T16": _t16, "T17": _t17, "T18": _t18, "T19": _t19, "T20": _t20, "T21": _t21,
    "T22": _t22}


# --- figures -------------------------------------------------------------------------------------

def _bars(ax, rows: Sequence[tuple[str, Sequence[float]]], title: str) -> None:
    left = np.zeros(len(rows))
    for j in C.MRS_LEVELS:
        values = np.array([pmf[j] for _, pmf in rows])
        ax.barh(range(len(rows)), values, left=left, label=f"mRS {j}", color=plt.cm.viridis(j / 6))
        left += values
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels([label for label, _ in rows])
    ax.set_xlim(0, 1)
    ax.set_xlabel("P(mRS = j)")
    ax.set_title(title)


def _f01(run: Run):
    fig, ax = plt.subplots(figsize=(8, 2.5))
    _bars(ax, [(C.TREATMENT_LABELS[a], _pmf(run.primary.cumulative, a)) for a in _ARMS],
          "Weighted mRS distribution by arm, [§7] population")
    ax.legend(ncol=7, fontsize=7, loc="upper center", bbox_to_anchor=(0.5, -0.35))
    return fig, "Weighted per-category mRS distribution; the six RD_k are cumulative differences."


def _forest(ax, rows: Sequence[tuple[str, float, float, float, str]], xlabel: str, null: float) -> None:
    for i, (label, point, lo, hi, marker) in enumerate(rows):
        ax.plot([lo, hi], [i, i], color="black")
        ax.plot([point], [i], marker, color="black", markerfacecolor="black" if marker == "o" else "white")
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels([r[0] for r in rows])
    ax.axvline(null, color="grey", linestyle=":")
    ax.set_xlabel(xlabel)
    ax.invert_yaxis()


def _f02(run: Run):
    rows = []
    for family, marker in (("secondary", "o"), ("safety", "s")):
        for est in run.secondary.by_family()[family]:
            iv = run.boot.intervals.get(f"{est.outcome}.rd")
            if iv is None:
                continue
            suffix = f" ({DESCRIPTIVE_SAFETY})" if family == "safety" else ""
            rows.append((C.OUTCOMES[est.outcome].label + suffix, est.rd, iv.lo, iv.hi, marker))
            aug = run.boot.intervals.get(f"{est.outcome}.augmented")
            if aug is not None and est.augmented is not None:
                rows.append((f"  {MODEL_ASSISTED}", est.augmented, aug.lo, aug.hi, "^"))
    fig, ax = plt.subplots(figsize=(8, 0.4 * len(rows) + 1.5))
    _forest(ax, rows, "risk difference (bridging - EVT alone)", 0.0)
    ax.set_title(f"Binary outcomes: filled = secondary, hollow = {DESCRIPTIVE_SAFETY} safety, "
                 f"triangle = {MODEL_ASSISTED}", fontsize=9)
    return fig, (f"Secondary and safety risk differences with percentile intervals; safety rows are "
                 f"{DESCRIPTIVE_SAFETY}; augmented points are {MODEL_ASSISTED}.")


def _f03(run: Run):
    rows = []
    for key, est in run.subgroups.estimates.items():
        for level, k in ((0, "beta"), (1, "beta_plus_gamma")):
            iv = run.subgroups.intervals.get(f"{key}.{k}")
            if iv is not None:
                rows.append((f"{C.SUBGROUPS[key]} S={level}", est.or_level[level], np.exp(iv.lo), np.exp(iv.hi), "o"))
    fig, ax = plt.subplots(figsize=(8, 0.5 * len(rows) + 1.5))
    _forest(ax, rows, "common odds ratio (log scale)", 1.0)
    ax.set_xscale("log")
    ax.set_title(f"Subgroups — {HYPOTHESIS_GENERATING}", fontsize=9)
    return fig, f"Level-specific odds ratios, {HYPOTHESIS_GENERATING}. {NON_COLLAPSIBLE}"


def _f04(run: Run):
    reported = [c for c in run.balance.centres if c.status == run.balance.pooled.status]
    omitted = [f"{c.centre}: {c.status}" for c in run.balance.centres if c.status != run.balance.pooled.status]
    df = run.cohort
    fig, axes = plt.subplots(1, max(len(reported), 1), figsize=(3 * max(len(reported), 1), 3), squeeze=False)
    bins = np.linspace(0, 1, 21)
    for ax, c in zip(axes[0], reported):
        at = (df["center"] == c.centre) & run.ps.in_model
        for a in _ARMS:
            ax.hist(run.ps.e[at & _arm_mask(df, a)].dropna(), bins=bins, alpha=0.5, label=C.TREATMENT_LABELS[a])
        ax.set_title(c.centre)
        ax.set_xlabel("propensity")
    axes[0][0].legend(fontsize=7)
    caption = "Propensity overlap by centre. Omitted: " + ("; ".join(omitted) if omitted else "none") + " [§9]."
    fig.suptitle(caption, fontsize=8)
    return fig, caption


def _f05(run: Run):
    fig, ax = plt.subplots(figsize=(8, 4))
    rows = [(f"primary [§7] {C.TREATMENT_LABELS[a]}", _pmf(run.primary.cumulative, a)) for a in _ARMS]
    rows += [(f"[§14a] {C.TREATMENT_LABELS[a]}", [run.std.distribution[a][j] for j in C.MRS_LEVELS]) for a in _ARMS]
    rows += [(f"[§14b] {regime}", [run.pol.distribution[regime][j] for j in C.MRS_LEVELS])
             for regime in run.pol.distribution]
    _bars(ax, rows, "mRS distributions: primary beside the [§14] transported analyses")
    ax.legend(ncol=7, fontsize=7, loc="upper center", bbox_to_anchor=(0.5, -0.2))
    return fig, f"{_populations_clause(run)} {TRANSPORTED}"


_FIGURES: Final[dict[str, Callable[[Run], tuple[object, str]]]] = {
    "F01": _f01, "F02": _f02, "F03": _f03, "F04": _f04, "F05": _f05}


# --- writing -------------------------------------------------------------------------------------

def _render_md(rows: Rows, clauses: Sequence[str]) -> str:
    data._assert_no_pipe(rows)  # noqa: SLF001 — V1 is data.py's rule, checked at the one write site
    return "\n".join([_md_table(rows), "", *(f"{c}\n" for c in clauses)])


def _render_csv(rows: Rows, clauses: Sequence[str]) -> str:
    buffer = io.StringIO()
    writer = csv.writer(buffer, lineterminator="\n")
    writer.writerows(rows)
    for clause in clauses:
        buffer.write(f"# {clause}\n")
    return buffer.getvalue()


def _render_svg(fig, caption: str) -> bytes:
    plt.rcParams.update(C.SVG_RC)
    buffer = io.BytesIO()
    fig.savefig(buffer, format="svg", bbox_inches="tight",
                metadata={"Date": None, "Creator": None, "Title": caption})
    plt.close(fig)
    return buffer.getvalue()


def _digest(content: bytes) -> str:
    return hashlib.sha256(content).hexdigest()


def write(run: Run, out: Path = C.OUT) -> Manifest:
    """Every output in `C.OUTPUT_IDS`, once each, then the audit log. Returns the manifest."""
    written: dict[str, bytes] = {}

    def put(relative: str, content: bytes) -> None:
        if relative in written:
            raise C.SchemaError(f"V2  {relative} rendered twice [Stage 14 §9].")
        written[relative] = content

    for oid, slug in C.OUTPUT_IDS.items():
        if oid in _TABLES:
            rows, clauses = _TABLES[oid](run)
            put(f"tables/{oid}_{slug}.md", _render_md(rows, clauses).encode("utf-8"))
            put(f"tables/{oid}_{slug}.csv", _render_csv(rows, clauses).encode("utf-8"))
            continue
        fig, caption = _FIGURES[oid](run)
        put(f"figures/{oid}_{slug}.svg", _render_svg(fig, caption))

    for relative, content in written.items():
        path = out / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(content)

    digests = {relative: _digest(content) for relative, content in sorted(written.items())}
    # Idempotent: a second `write` of the same Run replaces its entry rather than appending one.
    run.audit.entries[:] = [e for e in run.audit.entries
                            if not (e.kind == "provenance" and e.step == _STEP_WRITTEN)]
    run.audit.record("provenance", _STEP_WRITTEN, len(digests),
                     f"[§16] Stage 14 wrote {len(digests)} output file(s): "
                     f"{len(_TABLES)} tables as Markdown and CSV, {len(_FIGURES)} figures as SVG. "
                     "Every labelling statement is derived from a result object; the checklist is "
                     f"{C.OUTPUT_IDS['T21']}.",
                     table=(("file", "sha256"), *((k, v) for k, v in digests.items())))
    log = run.audit.write(out / "logs" / f"audit_{run.source.label}.md")
    digests[str(log.relative_to(out))] = _digest(log.read_bytes())
    return Manifest(tuple(Path(k) for k in sorted(digests)), digests)


def main() -> None:
    manifest = write(run())
    for path in manifest.paths:
        print(f"{manifest.sha256[str(path)]}  {path}")


if __name__ == "__main__":
    main()
