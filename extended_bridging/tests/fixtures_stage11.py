"""The constructed frames and records Stage 11 is tested on — §15.0 of the Stage 11 spec.

**Every fixture here exists because a branch of that document is unreachable on the workbook**, which
is `fixtures_stage9.py`'s and `fixtures_stage10.py`'s reason carried one stage on. A specification
that describes a guard it cannot reach has described an intention; a fixture that reaches it is what
turns the description into an assertion.

Six functions and one constant (§15.0):

    subgroup_frame()                  a cohort-shaped frame carrying the [§13] subgroup columns
                                      and the four vascular risk factors                    §15.0
    separable_subgroup_frame(kind)    the (A=1, S=0) cell empty -> O6;  no treated at S=1 -> G8
    near_separated_subgroup_frame()   ONE treated at S = 0, and the fit CONVERGES            §8.4
    family_fixtures()                 six (Secondary, Bootstrap) pairs, one per §6 branch    §6.5
    interval_fixtures()               seven `Interval`s on `beta`, one per §7 branch         §7.3
    tail_count_50()                   Stage 10's own fixture, REUSED and not rebuilt         §7.2

    SUBGROUP  — the C.SUBGROUPS key the constructed subgroup frames are built around

Nothing here reads the workbook, and no patient-derived value appears in it.

**One departure from §15.0's own prose, measured and recorded rather than absorbed.** §15.0 describes
the G8 companion as *"a frame in which every patient has `S = 1`, so `A x S` is constant"*. Measured:
with `S = 1` everywhere the product equals `A` exactly and is NOT constant — what is constant is `S`
itself, so `design` drops `S`, keeps `A` and `A x S` as two identical columns, and G8 fires on the
column tuple rather than on the product. The frame that makes the PRODUCT all-zero is the opposite
one — no treated patient at `S = 1` — and it is what `separable_subgroup_frame("constant_product")`
builds. Both reach G8 and only the second reaches it by the route §8.3 describes as route one, so
both are built and the tests say which is which.
"""
from __future__ import annotations

from typing import Final

import numpy as np
import pandas as pd

import bootstrap
import config as C
import outcome
from fixtures_stage10 import COVERAGE_SEED, boundary_draws, two_centre_frame

# The subgroup the constructed frames are built around. `unknown_onset` and not
# `core_above_median`, because it is the one whose treated cell the workbook itself makes small —
# five of forty-one — and a fixture built around the other would exercise the same code on a shape
# the cohort does not have.
SUBGROUP: Final[str] = "unknown_onset"


# --- 15.0  the cohort-shaped frame the subgroup fits run on -------------------------------------------

def subgroup_frame(n_a: int = 32, n_b: int = 28, seed: int = COVERAGE_SEED) -> pd.DataFrame:
    """`fixtures_stage10.two_centre_frame` plus what Stage 11 needs and Stage 10 did not.

    Two additions and nothing else:

      * **the four vascular risk factors and `penumbra_ml`**, so `C.PROPENSITY_FULL` has covariates
        to fit and `balance.assess`'s B3 has the whole of `C.BALANCE_SET` to judge. They are
        drawn rather than constant, because a constant column is dropped by `model.design` and the
        [§13] specification would then be the [§7] one with four names in `dropped` — a frame on
        which the arm is indistinguishable from the primary would pass a wrong implementation.
      * **the two [§13] subgroup columns**, as `Int64` with `<NA>` available, which is the dtype
        `derive.py` produces. `unknown_onset` is derived from `onset_type` through the same identity
        `derive._unknown_onset` uses, and `core_above_median` from this frame's own `core_ml` median
        with ties BELOW, which is `derive._core_above_median`'s rule. Neither is imported: Stage 11
        §8.2 requires that `sensitivity.py` name neither `derive_cohort` nor `_core_above_median`,
        and a fixture that imported them would be asserting the pipeline against itself.

    Everything Stage 10's frame guarantees is inherited unchanged — two strata, exactly one
    covariate-incomplete row, the seven [§5] binary outcomes, and treatment drawn from a mild logistic
    so that a replicate's propensity fit converges.
    """
    df = two_centre_frame(n_a, n_b, seed)
    g = np.random.default_rng(seed + 11)
    for name in C.NEGATIVE_CONTROLS:
        df[name] = g.integers(0, 2, len(df)).astype(float)
    # `penumbra_ml` is in `C.BALANCE_ONLY` and therefore in `C.BALANCE_SET`, so `balance.assess`'s
    # B3 raises without it — and this stage's arm calls `assess` (Stage 11 §5.1). It is in NEITHER
    # propensity specification, which is what makes it the one row reading `excluded [§6]` in both
    # tables and the only negative control the analysis has left after the arm spends the other four.
    df["penumbra_ml"] = np.round(g.uniform(10.0, 200.0, len(df)), 1)
    df[SUBGROUP] = (df["onset_type"] != C.REFERENCE_LEVELS["onset_type"]).astype("Int64")
    median = df["core_ml"].median()
    df["core_above_median"] = (df["core_ml"] > median).astype("Int64")
    return df


# --- 15.0  the two degenerate subgroup frames ---------------------------------------------------------

def separable_subgroup_frame(kind: str = "empty_cell",
                             subgroup: str = SUBGROUP) -> pd.DataFrame:
    """A frame on which one of §8.3's two routes fires. `kind` selects which, and they are different.

    `"empty_cell"` — **no treated patient at `S = 0`**, which is the route this cohort actually
    takes. `A x S` is then equal to `A` rather than constant, so `model.design` keeps BOTH columns
    and the design is rank-deficient, not degenerate: `model._assert_polr_fittable`'s **O6** catches
    it, raises `FitError`, and `C.FAILURE_BUCKETS["O6"]` already buckets it `degenerate_design`.
    Nothing new was needed for the route that fires, which is why §8.3 exists as a correction to an
    earlier draft rather than as a design.

    `"constant_product"` — **no treated patient at `S = 1`**, so `A x S` is all-zero and
    `model.design` drops it SILENTLY, after which `polr` would fit a two-column model, converge, and
    return a result in which `gamma` does not exist. This is §8.3's route one, measured at 0 of 2000
    replicates on the workbook and specified anyway. **G8 is the only thing that turns it into a
    failure**, and the assertion's PLACEMENT is what makes it one: between `design` and `polr`.

    `"one_level"` — **every patient at `S = 1`**, which is §15.0's own description of the G8
    companion and does NOT make the product constant (see this module's docstring). It reaches G8 by
    the other half of the same predicate: `S` is constant and is dropped, so the column tuple is not
    the declared one. It is built so that the difference between the two routes is a measured fact in
    the suite rather than a sentence here.
    """
    df = subgroup_frame()
    treated = df[C.TREATMENT] == 1.0
    values = pd.Series(pd.array([pd.NA] * len(df), dtype="Int64"), index=df.index)
    if kind == "empty_cell":
        values[:] = 1
        values[~treated] = (np.arange(len(df))[~treated] % 2).astype(int)
    elif kind == "constant_product":
        values[:] = 0
        values[~treated] = (np.arange(len(df))[~treated] % 2).astype(int)
    elif kind == "one_level":
        values[:] = 1
    else:
        raise AssertionError(
            f"{kind!r} is not one of 'empty_cell', 'constant_product', 'one_level'. §8.3 has exactly "
            "two routes and §15.0 describes the second two ways; a fixture reached by a typo is a "
            "branch nobody tested.")
    df[subgroup] = values
    return df


def near_separated_subgroup_frame(subgroup: str = SUBGROUP) -> pd.DataFrame:
    """ONE treated patient at `S = 0`, chosen so the fit CONVERGES and a reported quantity exceeds
    `C.POLR_MAX_ABS_BETA`.

    **This is the fixture the whole of §8.4 rests on.** The specification is that a fit which
    SUCCEEDS is rejected: `model.polr`'s own docstring records a separated proportional-odds fit
    converging in 17 iterations with every safeguard clean, `beta` 36.4 and `exp(beta)` 6.5e15. A
    frame on which the fit *failed* instead would pass a wrong implementation that has no G9 at all,
    because O1-O6 would already have raised.

    The construction: the `S = 0` level is made separated in treatment — its single treated patient
    takes the best outcome and every control at that level takes the worst — while the `S = 1` level
    keeps a spread. `|beta|` and `|gamma|` then blow up together while `|beta + gamma|` stays small,
    which is the signature §8.4 measures on the workbook: the `S = 0` level unidentified while the
    `S = 1` level is fine.

    The bound is on `max(|beta|, |gamma|, |beta + gamma|)` and deliberately NOT on `|delta|`, so this
    frame must not reach it through the subgroup MAIN EFFECT — §15.8 asserts which quantity crossed.
    """
    df = subgroup_frame().copy()
    treated = (df[C.TREATMENT] == 1.0).to_numpy()
    level = np.ones(len(df), dtype=int)
    # exactly one treated patient at S = 0, and a block of controls to separate it from
    level[np.flatnonzero(treated)[0]] = 0
    level[np.flatnonzero(~treated)[:12]] = 0
    df[subgroup] = pd.array(level, dtype="Int64")

    at_zero = level == 0
    outcomes = df[C.PRIMARY_OUTCOME].to_numpy(dtype=float)
    outcomes[at_zero & treated] = float(C.MRS_LEVELS[0])       # the one treated patient: best
    outcomes[at_zero & ~treated] = float(C.MRS_LEVELS[-1])     # every control there: worst
    df[C.PRIMARY_OUTCOME] = outcomes
    for key in C.DERIVED_DICHOTOMIES:                          # kept derived, never remembered
        entry = C.OUTCOMES[key]
        df[key] = C.OPS[entry.op](df[entry.source], float(entry.threshold)).astype(float)
    return df


# --- 15.0  the family fixtures: six (Secondary, Bootstrap) pairs --------------------------------------

def _binary_estimate(key: str) -> outcome.BinaryEstimate:
    """A minimal `BinaryEstimate` — this stage reads `outcome` and `family` and nothing else.

    Built by hand rather than fitted, because §6's whole subject is the p-values Stage 10 attached to
    these outcomes and not the estimates themselves. Every other field is a declared placeholder and
    §15.4 reads none of them.
    """
    entry = C.OUTCOMES[key]
    return outcome.BinaryEstimate(
        outcome=key, family=entry.family, minority=10, rd=0.0, odds_ratio=1.0, or_corrected=False,
        proportion={code: 0.5 for code in C.TREATMENT_LABELS}, augmented_path="unaugmented",
        augmented=None, covariates=None, reduced=False, dropped=(),
        in_estimate=pd.Series([], dtype=bool), fit=None)


def _draws(key: str, n: int) -> bootstrap.Draws:
    return bootstrap.Draws(quantity=key, draws=np.linspace(-1.0, 1.0, n), n_attempted=n,
                           failures={})


def _interval(p: float | None, n_draws: int = C.N_BOOT - 2) -> bootstrap.Interval:
    return bootstrap.Interval(-0.5, 0.5, C.CI_LEVEL, C.PERCENTILE_METHOD, n_draws, p)


def family_fixtures() -> dict[str, tuple[outcome.Secondary, bootstrap.Bootstrap]]:
    """Six (Secondary, Bootstrap) pairs, one per branch of §6, keyed by what each is for.

        complete      every family member carries a p                                    §6.1
        thinned       one member's `.rd` INTERVAL was not emitted — the DRAWS are kept    §6.5
        extra_keys    `intervals` carries keys no family declares, so an implementation
                      counting intervals gets a different m and fails                     §6.1
        none_p        an `Interval` whose `p` is None — a Stage 10 contract break, H3      §6.1
        floor         two Intervals with DIFFERENT `n_draws` and a raw p at the LARGER
                      one's floor: it passes there, and the SAME p on the smaller key is
                      below ITS floor and raises H4                                        §6.1
        no_beta       no `"beta"` key at all, which must raise H5                          §6.1

    **The last three are hand-built `Interval`s and none is a thing `bootstrap_p` or `run` can
    produce**, which is the point: H4 checks a value the pipeline cannot generate, and a suite that
    only drove the pipeline would never reach it.
    """
    members = {key: _binary_estimate(key) for key in C.BINARY_OUTCOMES}
    sec = outcome.Secondary(estimates=members)
    keys = [f"{key}.rd" for key in C.BINARY_OUTCOMES]

    def bootstrapped(intervals: dict[str, bootstrap.Interval],
                     draws: dict[str, bootstrap.Draws] | None = None) -> bootstrap.Bootstrap:
        base = {key: _draws(key, C.N_BOOT - 2) for key in [*keys, "beta"]}
        return bootstrap.Bootstrap(C.SEED, C.N_BOOT, {**base, **(draws or {})}, intervals,
                                   bootstrap.diagnostics(()))

    complete = {key: _interval(0.02 + 0.01 * i) for i, key in enumerate(keys)}
    complete["beta"] = _interval(0.03)

    thinned = dict(complete)
    starved = keys[0]
    del thinned[starved]

    extra = dict(complete)
    extra.update({f"{key}.odds_ratio": _interval(None) for key in C.BINARY_OUTCOMES})

    none_p = dict(complete)
    none_p[keys[1]] = _interval(None)

    # **§15.4 states the two `n_draws` the wrong way round and this is the corrected construction.**
    # It gives "a raw p at exactly 1/1983 ... passes, and the same p with n_draws 1998 raises H4",
    # and the floor moves the OTHER way: 1/(n+1) SHRINKS as n grows, so 1/1983 is comfortably ABOVE
    # the 1998-draw key's floor of 1/1999 and could never raise there. The construction that does
    # what §15.4 is about — one key passing, the other raising, and a check written against
    # `C.N_BOOT` accepting BOTH — is a p at the LARGER key's floor: 1/1999 is exactly achievable at
    # 1998 draws, is BELOW 1/1983 at 1982 draws, and is above `1/(C.N_BOOT + 1)` = 1/2001 either way.
    floored = dict(complete)
    floored[keys[0]] = bootstrap.Interval(
        -0.5, 0.5, C.CI_LEVEL, C.PERCENTILE_METHOD, 1998, 1.0 / 1999)
    floored[keys[1]] = bootstrap.Interval(
        -0.5, 0.5, C.CI_LEVEL, C.PERCENTILE_METHOD, 1982, 1.0 / 1999)

    no_beta = {key: interval for key, interval in complete.items() if key != "beta"}

    return {
        "complete": (sec, bootstrapped(complete)),
        "thinned": (sec, bootstrapped(thinned, {starved: _draws(starved, C.ci_min_draws() - 1)})),
        "extra_keys": (sec, bootstrapped(extra)),
        "none_p": (sec, bootstrapped(none_p)),
        "floor": (sec, bootstrapped(floored)),
        "no_beta": (sec, bootstrapped(no_beta)),
    }


# --- 15.0  the interval fixtures: seven `Interval`s on `beta` -----------------------------------------

def interval_fixtures() -> dict[str, bootstrap.Interval]:
    """Seven `Interval`s on `beta`, one per branch of §7.2 and §7.3, keyed by what each is for.

        spanning        lo < 0 < hi                                     e_limit is 1.0
        lo_at_null      lo == 0.0 exactly                               a TIE, and it spans
        hi_at_null      hi == 0.0 exactly                               a TIE, and it spans
        both_at_null    lo == hi == 0.0                                 a TIE, and it spans
        wide            lo << 0 << hi                                   e_limit is 1.0
        asymmetric      0 < lo << hi, so taking `hi` gives a VISIBLY different E-value
        outside         limits that EXCLUDE the point estimate          H7's raise

    **`asymmetric` is the fixture that makes the nearest-limit selection falsifiable.** A
    `min(|lo|, |hi|)` implementation picks the wrong limit ONLY in the spanning case, where the 1.0
    rule then masks it — so the defect is invisible until the first analysis whose interval excludes
    the null, which is the analysis nobody wants to find it in.
    """
    def interval(lo: float, hi: float, p: float | None = 0.05) -> bootstrap.Interval:
        return bootstrap.Interval(lo, hi, C.CI_LEVEL, C.PERCENTILE_METHOD, C.N_BOOT - 2, p)

    return {
        "spanning": interval(-0.4, 0.6, 0.30),
        "lo_at_null": interval(0.0, 0.6, 1.0 - C.CI_LEVEL),
        "hi_at_null": interval(-0.4, 0.0, 1.0 - C.CI_LEVEL),
        "both_at_null": interval(0.0, 0.0, 1.0),
        "wide": interval(-4.0, 4.0, 0.90),
        "asymmetric": interval(0.10, 2.50, 0.01),
        "outside": interval(0.80, 2.50, 0.01),
    }


# --- 15.0  Stage 10's own fixture, REUSED and not rebuilt ---------------------------------------------

def tail_count_50() -> np.ndarray:
    """`fixtures_stage10.boundary_draws(50)`, and it is imported rather than rebuilt.

    **The one fixture where route 3 of §7.2 can be told apart from routes 1 and 2.** Stage 10 §9.4
    measured p-versus-interval agreement HOLDING there under `inverted_cdf` and FAILING under
    `linear`, so it is the frame on which a wrong percentile method — the only thing that breaks
    routes 2 and 3 — is visible. Asserting the three routes agree on the workbook's own draws would
    pass under `linear` too and would therefore assert nothing.
    """
    return boundary_draws(50)
