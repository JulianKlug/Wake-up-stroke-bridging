"""The constructed frames Stage 10 is tested on — §15.0 of `specs/stage10_bootstrap_engine.md`.

**Every fixture here is the specification, not a convenience**, which is `fixtures_stage9.py`'s reason
carried one stage on: Stage 9 §22.3 item 1 records three fixtures pinned to ten decimals and described
rather than given, which made its acceptance criteria unperformable. Every pin in Stage 10 §15 is
measured against the code below.

Six functions and five constants (§13):

    two_centre_frame()          two strata, EXACTLY ONE covariate-incomplete row      §5.1, §5.2
    separable_ordinal_frame()   treatment determines mRS entirely; drives G7          §15.5
    constant_outcome_frame()    one outcome constant on its [§11] population; S8      §15.7
    unfittable_nuisance_frame() exactly ONE augmented outcome's m_a(X) cannot fit     §15.7
    boundary_draws()            `tail` draws below zero and the rest above            §8.2
    known_effect_population()   §8.4's coverage generator, confounded or not          §15.14

    COVERAGE_SEED  M_OUTER  B_INNER  BETA_TRUE  ALPHA_TRUE

Nothing here reads the workbook, and no patient-derived value appears in it [§4.5].

**Two departures from §15.0's own prose, both measured, both recorded rather than absorbed.**

1. `two_centre_frame` is 60 records at stratum sizes 32 and 28, not 8 and 6, and it carries the seven
   [§5] binary outcomes as well as the primary. §15.0 gives 8 and 6 and says the two give
   `P(a given row drawn >= 2x)` of 0.37 and 0.42. Neither half survives execution:

     * The probability is **0.2637 at a stratum of 8 and 0.2632 at 6** — the same 0.264 the cohort's
       stratum of 31 gives (§4.2's 0.2642), because `1 - (1-1/m)^m - m(1/m)(1-1/m)^(m-1)` is flat in
       `m` above about 5. So the blocker fires in about one replicate in four at ANY stratum size and
       the fixture does not need small strata to reach it. Measured here: **0 of 400** replicates
       raise `SchemaError` through the renaming resampler and **98 of 400 (24.5%)** through one
       without the rename, which is §15.2's assertion and its companion.
     * 8 + 6 = 14 records **cannot support the [§7] propensity model at all**. `model.design` over
       `C.PS_COVARIATES` on this frame is 11 parameters including the intercept, and a bootstrap
       replicate of 14 rows holds about 9 distinct ones — so `model.firth` raises F3 on a rank
       deficiency in essentially every replicate, and §15.2's "`propensity.fit` does not raise on 400
       consecutive replicates" is unsatisfiable. 60 records leave about 38 distinct rows per
       replicate against those 11 parameters, and 400 of 400 fit.

   The frame carries the binary outcomes because `bootstrap.run` refits Stages 6, 8 AND 9 in every
   replicate (§6.1), so a frame carrying only the primary cannot reach `_replicate` at all — and a
   second cohort-level fixture beside this one would be a second place this stage's fixture data
   lives, which is Stage 8 §14.0's reason for declining one.

2. `known_effect_population`'s design C attenuates to **0.5603** and not §8.4's 0.659272, and the
   confounding is deliberately STRONGER than §8.4's. §8.4 gives the design as *"A depends on Z and X;
   Y depends on A alone, or on all three"* and does not give the coefficients, so the estimand is a
   property of a construction that document never wrote down. It is therefore **measured at
   n = 200 000 by §15.14 rather than pinned here**, which is what §8.4 itself prescribes — *"the
   truth has to be the estimand, obtained the way Stage 9 §8.5 obtained its ATO, by running the
   estimator once at very large n"*.

   The coefficient is 1.0 rather than the value that reproduces §8.4's 0.659272, and that is a
   measurement rather than a preference: at an attenuation of 0.035 — which is §8.4's own gap — the
   two scorings agreed on 40 of 40 intervals at `M_OUTER` = 40, so §15.14's companion would be a
   coin flip rather than a test. `known_effect_population` carries the sweep. The property the tests
   assert is the one §8.4 is about: the truth is materially below `BETA_TRUE`, because a logistic
   model is not collapsible, and scoring design C against `BETA_TRUE` covers materially worse.
   Measured at M_OUTER = B_INNER = 300: design C covers **0.9500** against its own estimand and
   **0.8600** against `BETA_TRUE`, which is §8.4's "about 88%" reproduced.
"""
from __future__ import annotations

from typing import Final

import numpy as np
import pandas as pd

import config as C

# --- 15.0  the five constants ---------------------------------------------------------------------

COVERAGE_SEED: Final[int] = 20260807
M_OUTER: Final[int] = 300
B_INNER: Final[int] = 300
BETA_TRUE: Final[float] = 0.7
ALPHA_TRUE: Final[tuple[float, ...]] = (-2.0, -1.0, 0.0, 1.0, 2.0)   # 6 categories, 5 cutpoints

# The mRS composition of `two_centre_frame`, level by level, and it is a SPECIFICATION rather than a
# spread that looked reasonable. The four derived dichotomies' minority cells are computed from it,
# and three of them have to clear C.RARE_MINORITY_THRESHOLD for §9.2's `full` path while `sich` and
# `ph2` must fall below it — so that the frame reproduces the workbook's own five-augmented,
# two-unaugmented split (Stage 9 §9.3) and therefore Stage 10 §3.3's twenty-six estimand keys.
#
#   <= 1  18 / 42     <= 2  27 / 33     == 6  12 / 48     >= 5  18 / 42
#
# A uniform spread over the seven levels gives `death_90d` a minority of about 9 — one below the
# threshold — which silently turns it unaugmented and takes the key count to 25. Measured.
_MRS_COMPOSITION: Final[tuple[int, ...]] = (9, 9, 9, 8, 7, 6, 12)     # levels 0..6, summing to 60


# --- 15.0  the two-stratum cohort ------------------------------------------------------------------

def two_centre_frame(n_a: int = 32, n_b: int = 28, seed: int = COVERAGE_SEED) -> pd.DataFrame:
    """A cohort-shaped frame with two strata and EXACTLY ONE covariate-incomplete row.

    The witness for §5.1 and §5.2. The shape that makes it one is the single incomplete row, not the
    stratum size: `propensity._record_exclusion` compares the number of DISTINCT excluded `case_id`s
    against the number of excluded ROWS, so two excluded rows collapsing to one name is what raises,
    and with exactly one incomplete row that happens precisely when it is drawn twice — probability
    0.264 in every stratum above about five members (see this module's docstring for why §15.0's
    0.37/0.42 does not reproduce).

    Carries `C.PS_COVARIATES`, `C.TREATMENT`, `case_id`, `center`, the [§5] primary outcome and the
    seven [§5] binary outcomes. No value is taken from the workbook.

    **The treatment is drawn from a MILD logistic and not from an alternating pattern**, and that is
    what makes 400 of 400 replicates fit. `eta ~ N(0, 0.6)` puts every propensity near 0.5, so no
    covariate subset separates the arms in a replicate; a first draft alternated `ivt` on the row
    index with block-structured covariates and produced F5 — a fitted probability on the boundary —
    in 2 replicates of 400, which is a degenerate fit rather than the blocker this fixture is for.

    `case_id` and `center` are `string` dtype, as `COLUMN_CONTRACT` declares them, because §15.2
    asserts `resample` returns the input frame's dtypes column by column and `case_id` is the one
    that fails without §5.3's cast.
    """
    n = n_a + n_b
    g = np.random.default_rng(seed)
    mrs = np.repeat(np.arange(len(_MRS_COMPOSITION), dtype=float),
                    _rescaled(_MRS_COMPOSITION, n))
    g.shuffle(mrs)
    eta = g.normal(0.0, 0.6, n)
    df = pd.DataFrame({
        "case_id": pd.Series([f"TC-{i:02d}" for i in range(n)], dtype="string"),
        "center": pd.Series(["HUG"] * n_a + ["CHUV"] * n_b, dtype="string"),
        "onset_type": pd.Series(
            g.choice(["witnessed", "unwitnessed", "wake_up"], size=n, p=[0.6, 0.2, 0.2]),
            dtype="string"),
        "age": g.integers(50, 90, n).astype(float),
        "sex": g.integers(0, 2, n).astype(float),
        "prestroke_mrs": g.integers(0, 3, n).astype(float),
        "nihss_baseline": g.integers(4, 24, n).astype(float),
        "core_ml": np.round(g.uniform(0.0, 60.0, n), 1),
        "tmax6_ml": np.round(g.uniform(30.0, 140.0, n), 1),
        "atrial_fib": g.integers(0, 2, n).astype(float),
        C.TREATMENT: (g.random(n) < 1.0 / (1.0 + np.exp(-eta))).astype(float),
        C.PRIMARY_OUTCOME: mrs,
    })
    # The four derived dichotomies are DERIVED, through the registry's own operator and threshold,
    # so their missingness is `mrs_90d`'s exactly — which is roadmap invariant 6 and is the condition
    # §6.4's shared design rests on. Written as a registry loop rather than four expressions for
    # `config.DERIVED_DICHOTOMIES`' own reason: an outcome cannot be derived by being remembered.
    for key in C.DERIVED_DICHOTOMIES:
        entry = C.OUTCOMES[key]
        df[key] = C.OPS[entry.op](df[entry.source], float(entry.threshold)).astype(float)
    # The three read-directly outcomes. `tici_2b_3`'s minority is 8 and it is augmented ANYWAY,
    # because it is the one entry of OUTCOME_MODEL_OVERRIDES and an override is a route IN to
    # augmentation rather than a modifier of an outcome that qualified (outcome.py `_augmentable`).
    # `sich` at 9 and `ph2` at 8 fall below RARE_MINORITY_THRESHOLD and declare no override, so both
    # are unaugmented — which is the workbook's own split. `sich`'s period is 7 and not 6, and the
    # difference is one event: at 6 it lands on exactly 10, which IS the threshold, and `_augmentable`
    # compares `>=` — so the outcome the split needs unaugmented comes back `full` and the key count
    # is 27 rather than §3.3's 26. Measured.
    df["tici_2b_3"] = 1.0 - _every_kth(n, 8)
    df["sich"] = _every_kth(n, 7)
    df["ph2"] = _every_kth(n, 8)
    # EXACTLY ONE covariate-incomplete row, and it is in the second stratum. `core_ml` and not a
    # factor: a missing factor level raises D4 out of `model.design` before `complete_cases` is ever
    # consulted, so it would break the frame rather than shrink the fitted set by one.
    df.loc[n - 3, "core_ml"] = np.nan
    return df


def _rescaled(composition: tuple[int, ...], n: int) -> np.ndarray:
    """`composition` scaled to sum to `n`, the remainder going to the largest levels.

    A private helper and not one of §15.0's six functions: it exists so `two_centre_frame` takes an
    `n_a`/`n_b` a test can move without the mRS composition silently becoming a different one.
    """
    total = sum(composition)
    counts = np.array([c * n // total for c in composition])
    for j in np.argsort(-np.asarray(composition))[:n - int(counts.sum())]:
        counts[j] += 1
    return counts


def _every_kth(n: int, k: int) -> np.ndarray:
    """A 0/1 vector with a 1 every `k`-th position — n // k events, deterministically.

    Not a draw: the minority cells of the three read-directly outcomes decide §9.2's path for each,
    so they are counted rather than sampled.
    """
    return (np.arange(n) % k == 0).astype(float)


# --- 15.0  the separated ordinal frame -------------------------------------------------------------

def separable_ordinal_frame(n: int = 20) -> tuple[pd.DataFrame, np.ndarray, np.ndarray]:
    """Treatment determines mRS entirely: `a = 1` -> Y = 0, `a = 0` -> Y = 6, unit weights.

    Returns `(X, y, w)` — the design `model.polr` takes, the response and the weights — because that
    is what drives the fitter directly, and §15.5's whole point is that the fit SUCCEEDS.

    Drives G7, and Stage 8 §6.1's measurement is reproduced on it (§21): converges on the SCORE
    criterion in 17 iterations with rescales and halvings both 0, every fitted quantity finite,
    beta = 36.4058, exp(beta) = 6.469e15, alpha = -18.2029. That is the whole point of the fixture:
    separation does not present as failure, so G7 is the only thing that turns it into one, and a
    "failures" counter without G7 in it reads zero here (§15.5).

    The upper level is `C.MRS_LEVELS[-1]` and Stage 8's own `separated_frame` uses 5. The two give
    IDENTICAL fits and the reason is worth stating rather than leaving as a coincidence: `polr`
    collapses the response to the categories carrying positive weight (Stage 8 §5.3), so a two-category
    response is fitted on its own ordering and the numeric labels never enter the likelihood. §15.5
    asserts the pinned numbers, so if that ever stops being true the fixture says so rather than the
    pin quietly moving.
    """
    a = np.concatenate([np.ones(n), np.zeros(n)])
    y = np.where(a == 1.0, 0.0, float(C.MRS_LEVELS[-1]))
    return pd.DataFrame({C.TREATMENT: a}), y, np.ones(len(y))


# --- 15.0  the two frames that break one estimate ---------------------------------------------------

def constant_outcome_frame(key: str) -> pd.DataFrame:
    """A frame on which `key` takes one value across its whole [§11] population.

    Drives S8, at both its classifications: §15.7 asserts `model.FitError` after §5.4's change and
    that `C.SchemaError` is NOT raised, which is the assertion that fails on the landed code and is
    therefore the one that has to exist.

    Constant at **0.0** and not at 1.0, because that is the shape a [§10] replicate reaches: `sich`
    has five events in ninety-two and a stratified resample drawing none of them gives an all-zero
    column, which is the 0.8% §5.4 measured.
    """
    df = two_centre_frame()
    if key not in C.BINARY_OUTCOMES:
        raise AssertionError(
            f"{key!r} is not a C.BINARY_OUTCOMES key, and S8 is a per-outcome precondition of the "
            f"[§8] binary estimators alone. Known: {list(C.BINARY_OUTCOMES)}.")
    df[key] = 0.0
    return df


def unfittable_nuisance_frame() -> pd.DataFrame:
    """A frame on which EXACTLY ONE augmented outcome's `m_a(X)` cannot be fitted and six can.

    Drives `collect` (§7.3, §15.7). "Exactly one" is the specification: the test asserts six
    estimates and one recorded failure, and a frame on which two fail would pass a wrong
    implementation that gives up after the first.

    **The outcome is `tici_2b_3` and the mechanism is a rank deficiency the constant-column rule
    cannot repair.** Its [§11] population is restricted to the rows where `atrial_fib == ivt`, so the
    treatment column `outcome_model` INSERTS is identical to a column `model.design` already built —
    and `design` drops CONSTANT columns, not collinear ones, so both survive and `model.firth` raises
    F3 on the rank. Three properties make this the right one outcome of the seven:

      * `tici_2b_3` is the one entry of `OUTCOME_MODEL_OVERRIDES`, so its `m_a(X)` is four parameters
        over its own population and the deficiency is exact rather than marginal;
      * it is augmented through the override rather than through its minority cell, so restricting
        its population cannot accidentally turn it unaugmented and make the failure disappear;
      * `mrs_90d`'s own missingness is untouched, so the four full-list outcomes still share one
        design (§6.4) and §15.9's optimisation is not perturbed by this frame.

    S6, S7 and S8 all pass on it — the population is non-empty, both arms are present and the
    response is not constant — which is what makes the failure a `FitError` reaching `collect`
    rather than a `SchemaError` that `collect` may not catch at either setting.
    """
    df = two_centre_frame()
    collinear = df["atrial_fib"].to_numpy(dtype=float) == df[C.TREATMENT].to_numpy(dtype=float)
    df.loc[~collinear, "tici_2b_3"] = np.nan
    return df


# --- 15.0  the two arrays ---------------------------------------------------------------------------

def boundary_draws(tail: int, n: int = C.N_BOOT, rng: np.random.Generator | None = None):
    """`tail` draws strictly below zero and `n - tail` strictly above. NO draw equal to zero.

    The §8.2 fixture, and `tail` is a parameter because the finding is entirely about ONE value of
    it: at n = 2000 the p-value is a multiple of 2/n, so p < 0.05 iff tail <= 49, and tail = 50 is
    where p is exactly 0.0500 and `method="linear"` interpolates the limit onto the wrong side. The
    test sweeps 49, 50, 51 and asserts the verdict at each -- 50 alone distinguishes the methods, and
    a test written at 45 or 55 passes under every one of them (§15.8).

    `rng` defaults to a generator seeded on `COVERAGE_SEED`, so the array is a property of this
    module. The magnitudes are drawn and the SIGNS are not: exactly `tail` negatives is the whole
    construction, so they are assigned rather than sampled, and `+ 0.01` keeps every draw strictly
    off zero — a draw at exactly 0.0 is counted in both tails (§9.1) and would move the p-value the
    fixture exists to pin.
    """
    if not 0 <= tail <= n:
        raise AssertionError(f"tail must be in 0..{n}; got {tail}")
    g = np.random.default_rng(COVERAGE_SEED) if rng is None else rng
    magnitude = g.uniform(0.01, 1.0, n) + 0.01
    sign = np.where(np.arange(n) < tail, -1.0, 1.0)
    return sign * magnitude


# --- 15.0  the coverage generator -------------------------------------------------------------------

def known_effect_population(n: int, rng: np.random.Generator, confounded: bool) -> pd.DataFrame:
    """§8.4's coverage generator. `A` depends on Z and X; `Y` depends on A alone, or on all three.

    `confounded=False` is design U: the marginal odds ratio IS BETA_TRUE in closed form, because `Y`
    is generated from `logit P(Y <= k) = ALPHA_TRUE[k] + BETA_TRUE * A` and nothing else, so every
    weighting of the population estimates the same coefficient and coverage measures the machinery.
    `confounded=True` is design C, where the conditional coefficient is NOT the marginal ATO odds
    ratio -- logistic models are not collapsible -- and the truth must be obtained by running the
    estimator once at very large `n`, as Stage 9 §8.5 obtained its ATO. Measured at n = 200 000:
    **0.5603** against BETA_TRUE = 0.7. That is a larger attenuation than §8.4's 0.659272 and the
    comment on `eta` below says why it is chosen rather than inherited; §15.14 recomputes the truth
    rather than pinning it, so the number here is a record and not an input.

    Returned as a **cohort-shaped frame** rather than as arrays, because §8.4's estimator is the [§8]
    one -- `propensity.fit` then `outcome.primary` -- and both read columns by name. `Z` and `X` are
    carried as `age` and `nihss_baseline`, which are the two continuous entries of `C.PS_COVARIATES`;
    the other seven covariates are held CONSTANT, so `model.design` drops every one of them as a
    constant column and the propensity model is exactly `intercept + Z + X`. That is the design being
    correctly specified, stated as a fact about the frame rather than hoped for.

    **These designs have no strata** and §8.4 says so: `center` is constant, so `resample`'s
    stratification is not exercised by the coverage measurement and §16 item 4 files that as the gap
    that matters.
    """
    z = rng.normal(0.0, 1.0, n)
    x = rng.normal(0.0, 1.0, n)
    a = (rng.random(n) < 1.0 / (1.0 + np.exp(-(0.5 * z + 0.5 * x)))).astype(float)
    # The confounder coefficient is 1.0 and it is CHOSEN, not conventional. §8.4's own design
    # attenuates 0.700000 to 0.659272 -- a gap of 0.041 against a mean interval width of 0.7366 --
    # and §15.14's companion asks that scoring against BETA_TRUE cover MATERIALLY worse. Measured at
    # M_OUTER = 40: at a coefficient of 0.5 the attenuation is 0.035 and the two scorings agree on
    # 40 of 40 intervals, so the companion is a coin flip rather than a test. Measured at n = 200 000:
    # 0.5 -> 0.6652, 1.0 -> 0.5603, 1.5 -> 0.4594, 2.0 -> 0.3813. 1.0 puts the gap at 0.140 against
    # a width of about 0.73, which is the smallest value at which the companion is unmistakable.
    eta = BETA_TRUE * a + (1.0 * (z + x) if confounded else 0.0)
    # logit P(Y <= k) = ALPHA_TRUE[k] + eta, inverted on a uniform draw. The cumulative
    # probabilities are non-decreasing in k by construction, so `searchsorted` over them per row is
    # the ordinal inverse-CDF and never produces a category outside 0..len(ALPHA_TRUE).
    cumulative = 1.0 / (1.0 + np.exp(-(np.asarray(ALPHA_TRUE)[None, :] + eta[:, None])))
    y = (rng.random(n)[:, None] > cumulative).sum(axis=1).astype(float)
    return pd.DataFrame({
        "case_id": pd.Series([f"KE-{i:05d}" for i in range(n)], dtype="string"),
        "center": pd.Series([C.CENTER_ORDER[0]] * n, dtype="string"),
        "onset_type": pd.Series([C.REFERENCE_LEVELS["onset_type"]] * n, dtype="string"),
        "age": z,                       # Z, the first confounder
        "nihss_baseline": x,            # X, the second
        "sex": np.zeros(n),
        "prestroke_mrs": np.zeros(n),
        "core_ml": np.zeros(n),
        "tmax6_ml": np.zeros(n),
        "atrial_fib": np.zeros(n),
        C.TREATMENT: a,
        C.PRIMARY_OUTCOME: y,
    })
