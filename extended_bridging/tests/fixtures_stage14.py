"""Stage 14 §12.0 — a synthetic `report.Run` built from dataclasses, with NO fitting.

Two knobs decide which data-derived clauses fire: `direction` ("one" | "mixed") sets whether the
six RD_k share a sign, `separation` sets how many subgroup replicates a bound removed. Every other
field carries a shape the renderers read and a value nothing asserts.
"""

from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd

import balance
import bootstrap
import config as C
import data
import outcome
import policy
import propensity
import report
import sensitivity
import standardise

N: int = 40
SEED: int = 14


def _frame(rng: np.random.Generator) -> pd.DataFrame:
    df = pd.DataFrame({
        C.TREATMENT: rng.integers(0, 2, N),
        "center": rng.choice(C.CENTER_ORDER[:3], N),
        "age": rng.normal(70, 10, N), "sex": rng.integers(0, 2, N),
        "prestroke_mrs": rng.integers(0, 3, N), "nihss_baseline": rng.integers(2, 25, N),
        "onset_type": rng.choice(C.FACTOR_LEVELS["onset_type"], N),
        "core_ml": rng.gamma(2.0, 10.0, N), "tmax6_ml": rng.gamma(3.0, 30.0, N),
        "atrial_fib": rng.integers(0, 2, N), "hypertension": rng.integers(0, 2, N),
        "hyperlipidemia": rng.integers(0, 2, N), "diabetes": rng.integers(0, 2, N),
        "smoking": rng.integers(0, 2, N), "penumbra_ml": rng.gamma(3.0, 20.0, N),
        "onset_to_groin_min": rng.integers(200, 900, N),
        C.PRIMARY_OUTCOME: rng.integers(0, 7, N),
    }).astype({"center": "string", "onset_type": "string"})
    for key in C.BINARY_OUTCOMES:
        df[key] = rng.integers(0, 2, N)
    return df


def _draws(keys: list[str], rng: np.random.Generator, n: int = 60,
           separation: dict[str, int] | None = None) -> tuple[dict, dict]:
    separation = separation or {}
    draws, intervals = {}, {}
    for key in keys:
        lost = separation.get(key, 0)
        values = rng.normal(0.0, 0.1, n - lost)
        draws[key] = bootstrap.Draws(key, values, n, {"separation": lost} if lost else {})
        intervals[key] = bootstrap.Interval(float(values.min()), float(values.max()), C.CI_LEVEL,
                                            C.PERCENTILE_METHOD, n - lost, 0.5 if key.endswith("beta") or key.endswith(".rd") or key.endswith("gamma") else None)
    return draws, intervals


def _cumulative(rd: dict[int, float]) -> dict[int, dict[int, float]]:
    base = np.linspace(0.1, 0.9, len(C.MRS_THRESHOLDS))
    return {k: {0: float(base[i]), 1: float(min(0.99, base[i] + rd[k]))}
            for i, k in enumerate(C.MRS_THRESHOLDS)}


def _distribution(cumulative: dict[int, dict[int, float]], arms) -> dict:
    out = {}
    for arm in arms:
        previous, pmf = 0.0, {}
        for k in C.MRS_THRESHOLDS:
            pmf[k] = cumulative[k][arm] - previous
            previous = cumulative[k][arm]
        pmf[C.MRS_LEVELS[-1]] = 1.0 - previous
        out[arm] = pmf
    return out


def _balance(spec, exceed: bool) -> balance.Balance:
    rows = []
    for name in C.BALANCE_SET:
        for label, _ in balance.levels(pd.DataFrame({name: pd.Series([0], dtype="string" if name in C.CATEGORICAL else float)}), name):
            rows.append(balance.CovariateBalance(label, "propensity model", N, 1.0, 0.3,
                                                 0.25 if exceed and label == "diabetes" else 0.02))
    overlap = balance.CentreOverlap("HUG", 20, 20, {0: 10, 1: 10}, {0: (0.2, 0.8), 1: (0.3, 0.9)},
                                    {0: 8.0, 1: 8.0}, 0.9, 0.2, 0.05, balance._REPORTED)
    structural = balance.CentreOverlap("USZ", 0, 0, {0: 0, 1: 0}, {}, {}, float("nan"), 0.0,
                                       float("nan"), balance._STRUCTURAL)
    centres = tuple(structural if c == "USZ" else balance.CentreOverlap(
        c, 20, 20, {0: 10, 1: 10}, {0: (0.2, 0.8), 1: (0.3, 0.9)}, {0: 8.0, 1: 8.0}, 0.9, 0.2, 0.05,
        balance._REPORTED) for c in C.CENTER_ORDER)
    return balance.Balance(tuple(rows), centres, overlap, spec)


def _standardisation(label: str, rd: dict[int, float]) -> standardise.Standardisation:
    cumulative = _cumulative(rd)
    return standardise.Standardisation(label, N, N, _distribution(cumulative, (0, 1)), cumulative, rd,
                                       rd[2], -rd[5], 0.1, float(np.exp(0.1)), standardise._CONDITIONAL,
                                       None, ())


def synthetic_run(direction: str = "mixed", separation: int = 0) -> report.Run:
    rng = np.random.default_rng(SEED)
    df = _frame(rng)
    source = data.Source(Path("synthetic.xlsx"), None, None, "synthetic")
    audit = data.Audit(source)
    audit.record("cohort", "cohort_flow", N, "synthetic flow",
                 table=(("step", "centres", "records", "EVT alone", "bridging", "eligible", "ineligible"),
                        ("as classified", "3", str(N), "20", "20", str(N), "0")))

    in_model = pd.Series(True, index=df.index)
    e = pd.Series(rng.uniform(0.1, 0.9, N), index=df.index)
    ps = propensity.Propensity(e, 1.0 - e, in_model, {0: 15.0, 1: 14.0}, None, (), C.PROPENSITY_PRIMARY)
    ps_full = propensity.Propensity(e, 1.0 - e, in_model, {0: 15.0, 1: 14.0}, None, (), C.PROPENSITY_FULL)

    signs = [1, -1, 1, -1, 1, -1] if direction == "mixed" else [1] * 6
    rd = {k: 0.05 * s for k, s in zip(C.MRS_THRESHOLDS, signs)}
    cumulative = _cumulative(rd)
    in_estimate = in_model & df[C.PRIMARY_OUTCOME].notna()
    primary = outcome.Primary(-0.1, float(np.exp(-0.1)), (0.0,) * 6, C.MRS_THRESHOLDS, rd, cumulative,
                              in_estimate, None)

    estimates = {}
    for key in C.BINARY_OUTCOMES:
        o = C.OUTCOMES[key]
        augmented = o.family == "secondary"
        estimates[key] = outcome.BinaryEstimate(
            key, o.family, 5, 0.02, 1.1, key == "sich", {0: 0.3, 1: 0.32},
            ("reduced" if key in C.OUTCOME_MODEL_OVERRIDES else "full") if augmented else "unaugmented",
            0.03 if augmented else None, ("age",) if augmented else None,
            key in C.OUTCOME_MODEL_OVERRIDES, (), in_model, None)
    secondary = outcome.Secondary(estimates)

    keys = ["beta", *(f"rd_{k}" for k in C.MRS_THRESHOLDS)]
    for key in C.BINARY_OUTCOMES:
        keys += [f"{key}.rd", f"{key}.odds_ratio"]
        if estimates[key].augmented_path != "unaugmented":
            keys.append(f"{key}.augmented")
    draws, intervals = _draws(keys, rng)
    diagnostics = bootstrap.Diagnostics({6: 60}, {4: 60}, np.full(60, 40.0), {}, {40: 60},
                                        {k: 0 for k in C.BINARY_OUTCOMES})
    boot = bootstrap.Bootstrap(C.SEED, 60, draws, intervals, diagnostics)

    families = {}
    for family, members in secondary.by_family().items():
        raw = {m.outcome: 0.4 for m in members}
        families[family] = sensitivity.FamilyCorrection(family, len(members), len(members), raw,
                                                        dict(raw), (), family == "safety")
    multiplicity = sensitivity.Multiplicity(families, 0.5)
    e_value = sensitivity.EValue("common odds ratio [§8]", primary.odds_ratio, 0.95, 1.2, True, None,
                                 1.0, 60, sensitivity._APPROXIMATION)

    sub_keys = [f"{s}.{q}" for s in C.SUBGROUPS for q in ("beta", "gamma", "beta_plus_gamma")]
    first = next(iter(C.SUBGROUPS))
    sub_draws, sub_intervals = _draws(sub_keys, rng, separation={f"{first}.{q}": separation for q in (
        "beta", "gamma", "beta_plus_gamma")})
    sub_estimates = {s: sensitivity.SubgroupEstimate(
        s, N, {"S=0,A=0": 10, "S=0,A=1": 10, "S=1,A=0": 10, "S=1,A=1": 10}, {0: 0.9, 1: 1.1}, 1.2, 0.18,
        {0: cumulative, 1: cumulative}, True, len(C.SUBGROUPS)) for s in C.SUBGROUPS}
    subgroups = sensitivity.Subgroups(sub_estimates, C.SEED, 60, sub_draws, sub_intervals)

    arm_draws, arm_intervals = _draws(keys[:7], rng)
    arm_intervals = {k: bootstrap.Interval(i.lo, i.hi, i.level, i.method, i.n_draws, None)
                     for k, i in arm_intervals.items()}
    arm = sensitivity.Arm(ps_full, _balance(C.PROPENSITY_FULL, False), primary, ps, C.SEED, 60,
                          arm_draws, arm_intervals, diagnostics, True)

    std = _standardisation(standardise._ALL_ELIGIBLE, rd)
    std_support = _standardisation(standardise._TREATED_SUPPORT, rd)
    hier = standardise.Hierarchical(_standardisation(standardise._RANDOM_INTERCEPT, rd), 0.2, False,
                                    {c: 0.0 for c in C.CENTER_ORDER}, {c: 0.3 for c in C.CENTER_ORDER},
                                    C.POLR_RI_NODES, None)
    std_keys = []
    for prefix in ("", "support.", "hier."):
        std_keys += [f"{prefix}rd_{k}" for k in C.MRS_THRESHOLDS] + [f"{prefix}mrs_0_2", f"{prefix}mortality"]
        std_keys += [f"{prefix}dist{a}_{j}" for a in (1, 0) for j in C.MRS_LEVELS]
    std_keys += ["beta", "hier.beta", "hier.sigma"]
    std_draws, std_intervals = _draws(std_keys, rng)
    std_intervals = {k: bootstrap.Interval(i.lo, i.hi, i.level, i.method, i.n_draws, None)
                     for k, i in std_intervals.items()}
    std_boot = bootstrap.Bootstrap(C.SEED, 60, std_draws, std_intervals, diagnostics)
    support = standardise.Support({"age": (50.0, 90.0)}, in_model, {c: 0 for c in C.CENTER_ORDER},
                                  pd.Series(False, index=df.index), 0, 0, _balance(C.PROPENSITY_PRIMARY, False).covariates)

    regimes = (policy._ACTIVE, policy._COMPARATOR)
    pol_cumulative = {k: {policy._ACTIVE: v[1], policy._COMPARATOR: v[0]} for k, v in cumulative.items()}
    pol = policy.Policy(N, N, N - 5, (N - 5) / N, _distribution(pol_cumulative, regimes), pol_cumulative,
                        rd, rd[2], -rd[5], {k: v / ((N - 5) / N) for k, v in rd.items()}, 0.1,
                        float(np.exp(0.1)), 0.5, policy._CONDITIONAL_14B, policy._OPERATIONAL,
                        SimpleNamespace(columns=(C.TREATMENT, "age", C.CONTRAINDICATED), beta=np.array([0.1, 0.01, 0.5])), ())
    pol_keys = [f"rd_{k}" for k in C.MRS_THRESHOLDS] + [f"eligible_rd_{k}" for k in C.MRS_THRESHOLDS]
    pol_keys += ["mrs_0_2", "mortality", "share_eligible", "beta"]
    pol_keys += [f"dist_{r}_{j}" for r in regimes for j in C.MRS_LEVELS]
    pol_draws, pol_intervals = _draws(pol_keys, rng)
    pol_intervals = {k: bootstrap.Interval(i.lo, i.hi, i.level, i.method, i.n_draws, None)
                     for k, i in pol_intervals.items()}
    pol_boot = bootstrap.Bootstrap(C.SEED, 60, pol_draws, pol_intervals, diagnostics)

    return report.Run(source, audit, df, df, ps, _balance(C.PROPENSITY_PRIMARY, True), primary, secondary,
                      boot, multiplicity, e_value, subgroups, arm, df, std, support, std_support, hier,
                      std_boot, df, pol, pol_boot)
