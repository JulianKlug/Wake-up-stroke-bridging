"""The constructed frames Stage 12 is tested on — §20.0 of `specs/stage12_all_centre_standardisation.md`.

**Every fixture here is the specification, not a convenience**, which is `fixtures_stage10.py`'s reason
carried two stages on: Stage 9 §22.3 item 1 records three fixtures pinned to ten decimals and
described rather than given, which made its acceptance criteria unperformable. Every pin in Stage 12
§20 is measured against the code below.

Five functions and four constants (§20.0):

    four_centre_frame()       four declared centres, ONE with no treated patient, all 7 mRS levels
    collapsed_level_frame()   the same with mRS 5 absent, so `polr` fits five cutpoints      §6.3
    no_exposure_frame()       every patient treated, so `design` drops the exposure; T1      §5.3
    two_group_ri_frame()      a known non-zero between-group intercept, for `polr_ri`        §12
    flat_ri_frame()           generated with sigma = 0, so `polr_ri` reaches the floor       §12.5

    RI_SEED  RI_SIGMA_TRUE  RI_BETA_TRUE  COLLAPSE_LEVEL

Nothing here reads the workbook, and no patient-derived value appears in it [§4.3].

**Every frame carries `C.ELIGIBILITY`, and that is what makes it a [§14a] frame rather than a [§3]
cohort one.** `standardise.population` applies restriction 2 and NOT restriction 1 (§4.1), so a
fixture without the column fails T6 before any [§14a] arithmetic runs — which is the assertion §20.1
makes, not an inconvenience to work around.

**`two_group_ri_frame` takes `n_groups` and the name is §20.0's rather than a bound.** Its default is
the two groups the name says, which is the smallest frame on which `sigma` is identified at all and
therefore the right one for the structural tests. **Recovering `sigma` is a different question from
fitting it**, and it is not answerable at two clusters: the recovery assertion in §20.12 drives the
same generator at twelve, which is where the estimate has something to be within a tolerance OF.
Measured, at `RI_SIGMA_TRUE = 0.8` and 60 records per group: sigma_hat 0.7405 at twelve groups
against 0.3862 at two — the second being an honest reading of two numbers, not a defect.
"""
from __future__ import annotations

from typing import Final

import numpy as np
import pandas as pd

import config as C

# --- §20.0  the four constants ---------------------------------------------------------------------

RI_SEED: Final[int] = 20260827
RI_SIGMA_TRUE: Final[float] = 0.8
RI_BETA_TRUE: Final[float] = 0.7

# The level `collapsed_level_frame` removes, and it is mRS 5 rather than any other for a measured
# reason: that is the level the WORKBOOK loses, in all 83 of the 2000 replicates that lose one
# (Stage 12 §6.2), because it carries 3 patients of 104. A fixture that removed a level the cohort
# never loses would exercise the re-expansion on a case the bootstrap does not produce.
COLLAPSE_LEVEL: Final[int] = 5

# The cutpoints every generated frame is drawn from, ascending, in [§14a]'s own parametrisation:
# `logit P(Y <= k) = alpha_k + x'beta + b_c`, so a LARGER linear predictor means a SMALLER mRS. Six
# cutpoints for the seven declared levels.
_ALPHA_TRUE: Final[tuple[float, ...]] = (-2.5, -1.5, -0.6, 0.2, 1.0, 2.0)

# `four_centre_frame`'s per-centre intercepts, and they are DELIBERATELY LARGE ENOUGH THAT
# `polr_ri` FITS A NON-ZERO SIGMA ON THIS FRAME. A fixture whose between-centre SD collapses to
# POLR_RI_SIGMA_FLOOR cannot exercise §12.6's conditional standardisation at all — the intercepts are
# then all zero and the hierarchical arm is arithmetically the pooled one, so a test asserting the two
# differ would be asserting nothing. Measured at these values: sigma_hat is well clear of the floor
# and `at_floor` is False. `flat_ri_frame` is the fixture for the opposite case.
_CENTRE_EFFECT: Final[tuple[float, ...]] = (0.95, -0.65, 0.55, -0.85)

# `four_centre_frame`'s per-centre sizes and treated counts. The SHAPE is the workbook's own and the
# numbers are not: one centre with no treated patient at all is what makes this a [§14a] fixture
# rather than a [§3] cohort one, because it is the configuration `cohort.build` exists to remove and
# [§14a] exists to include (§4.1, §10.2).
_CENTRES: Final[tuple[tuple[str, int, int], ...]] = (
    ("HUG", 40, 24),        # centre, records, treated
    ("CHUV", 24, 8),
    ("Lugano", 28, 4),
    ("USZ", 16, 0),         # THE NEVER-IVT CENTRE. Its complement is what §10.2 asserts against
)


def _ordinal_draw(g: np.random.Generator, eta: np.ndarray,
                  alpha: tuple[float, ...] = _ALPHA_TRUE) -> np.ndarray:
    """One ordinal response per row, drawn from `logit P(Y <= k) = alpha_k + eta_i`.

    The inverse-CDF draw, and the cumulative sequence is bracketed by exact 0 and exact 1 for
    `model.ordinal_probabilities`' reason: the category probabilities must sum to exactly 1 for the
    draw to be a draw from a distribution rather than from a distribution-shaped array.
    """
    cumulative = np.column_stack([
        np.zeros(len(eta)),
        1.0 / (1.0 + np.exp(-(np.asarray(alpha)[None, :] + eta[:, None]))),
        np.ones(len(eta)),
    ])
    return np.asarray(
        (g.random(len(eta))[:, None] > cumulative[:, 1:-1]).sum(axis=1), dtype=float)


def _covariates(g: np.random.Generator, n: int) -> dict[str, object]:
    """The [§6] covariate columns plus `C.BALANCE_ONLY`'s five, which §10.3's table ranges over.

    `BALANCE_ONLY` is carried because [§9]'s rule is that balance is judged against the FULL
    confounder set, and §10.3 applies it to [§14a]'s support-check table — so a frame without those
    five columns cannot produce the nineteen rows §20.10 counts.
    """
    return {
        "onset_type": pd.Series(
            g.choice(["witnessed", "unwitnessed", "wake_up"], size=n, p=[0.6, 0.2, 0.2]),
            dtype="string"),
        "age": g.integers(50, 95, n).astype(float),
        "sex": g.integers(0, 2, n).astype(float),
        "prestroke_mrs": g.integers(0, 4, n).astype(float),
        "nihss_baseline": g.integers(3, 26, n).astype(float),
        "core_ml": np.round(g.uniform(0.0, 62.0, n), 1),
        "tmax6_ml": np.round(g.uniform(13.0, 401.0, n), 1),
        "atrial_fib": g.integers(0, 2, n).astype(float),
        "hypertension": g.integers(0, 2, n).astype(float),
        "hyperlipidemia": g.integers(0, 2, n).astype(float),
        "diabetes": g.integers(0, 2, n).astype(float),
        "smoking": g.integers(0, 2, n).astype(float),
        "penumbra_ml": np.round(g.uniform(5.0, 200.0, n), 1),
    }


def four_centre_frame(seed: int = RI_SEED) -> pd.DataFrame:
    """An eligible [§14a] population over four declared centres, ONE of which never treats.

    108 records, 36 treated, all seven mRS levels occupied — asserted by the tests rather than
    assumed here, because the mRS draw is stochastic and a fixture that claimed a level set it did
    not have would make §6.2's collapsed-frame contrast meaningless.

    `case_id` and `center` are `string` dtype, as `COLUMN_CONTRACT` declares them, because
    `bootstrap.resample` asserts it returns the input frame's dtypes column by column and `case_id`
    is the one that fails without Stage 10 §5.3's cast.

    Every record is `C.ELIGIBLE`: [§14a]'s population is the eligible one and restriction 2 is the
    only [§3] restriction on this stage's path (§4.1). `no_exposure_frame` below is the one fixture
    that varies anything else.
    """
    g = np.random.default_rng(seed)
    centres = np.concatenate([[name] * n for name, n, _ in _CENTRES])
    treated = np.concatenate([
        np.concatenate([np.ones(k), np.zeros(n - k)]) for _, n, k in _CENTRES])
    n = len(centres)
    columns = _covariates(g, n)
    # A mild covariate effect and a per-centre intercept, so the frame is not separable and the
    # standardisation has something to transport. The centre effect is what makes USZ's counterfactual
    # a genuine extrapolation rather than an average of nothing.
    eta = (RI_BETA_TRUE * treated
           + 0.02 * (columns["age"] - 70.0)
           - 0.04 * (columns["nihss_baseline"] - 14.0)
           + np.repeat(_CENTRE_EFFECT, [n for _, n, _ in _CENTRES]))
    return pd.DataFrame({
        "case_id": pd.Series([f"S12-{i:03d}" for i in range(n)], dtype="string"),
        "center": pd.Series(centres, dtype="string"),
        **columns,
        C.TREATMENT: treated,
        C.PRIMARY_OUTCOME: _ordinal_draw(g, eta),
        C.ELIGIBILITY: pd.Series([C.ELIGIBLE] * n, dtype="string"),
    })


def collapsed_level_frame(seed: int = RI_SEED) -> pd.DataFrame:
    """`four_centre_frame` with `COLLAPSE_LEVEL` absent, so `polr` fits five cutpoints of six.

    The witness for §6.2 and §6.3: `ordinal_probabilities` returns SIX columns where
    `C.MRS_LEVELS` declares seven, the caller re-expands with a structural zero at the missing level,
    the distribution still sums to 1, the cumulative sequence is still monotone, and `RD_4 == RD_5`
    exactly.

    **The level is moved DOWN and not deleted**, so the record count, the per-centre sizes and the
    treated counts are `four_centre_frame`'s exactly. A fixture that dropped the rows would change
    two things at once and a difference between the two frames could not be attributed to the
    collapse.
    """
    df = four_centre_frame(seed).copy()
    at_level = df[C.PRIMARY_OUTCOME] == float(COLLAPSE_LEVEL)
    df.loc[at_level, C.PRIMARY_OUTCOME] = float(COLLAPSE_LEVEL - 1)
    return df


def no_exposure_frame(seed: int = RI_SEED) -> pd.DataFrame:
    """`four_centre_frame` with every patient treated, so `model.design` drops the exposure. T1.

    The failure this drives is SILENT without T1 (§5.3, §7.1): `design` drops the constant treatment
    column and returns its name, `polr` fits the remaining nine covariates happily, and the
    "counterfactual" prediction is then made under a model with no exposure in it. Nothing raises,
    and `beta` does not exist in a fit that reports one.
    """
    df = four_centre_frame(seed).copy()
    df[C.TREATMENT] = 1.0
    return df


def two_group_ri_frame(n_groups: int = 2, per_group: int = 60, sigma: float = RI_SIGMA_TRUE,
                       seed: int = RI_SEED) -> pd.DataFrame:
    """A frame with a known non-zero between-group intercept SD, for `polr_ri` [§14a sensitivity 2].

    Returns the design-shaped columns `polr_ri` needs and nothing else — `x`, `group`, `y` — because
    this fixture drives the ESTIMATOR and not the stage: `polr_ri` takes a design matrix, a response
    and a grouping vector, and knows nothing about centres, treatment or eligibility (§12.1).

    `sigma = 0` is a legitimate argument and is what `flat_ri_frame` passes. At `sigma = 0` every
    `b_c` is exactly zero, so the generating model IS the pooled model and [§14a]'s own *"it
    collapses to the pooled model"* becomes checkable rather than quotable.

    **The group intercepts are drawn ONCE from N(0, sigma^2) and are a property of the frame**, not
    re-drawn per row: `b_c` is a group-level parameter, and a per-row draw would generate
    overdispersion rather than a random intercept — which fits, and returns a sigma that estimates
    the wrong thing.
    """
    g = np.random.default_rng(seed)
    n = n_groups * per_group
    labels = np.repeat([f"G{i}" for i in range(n_groups)], per_group)
    b = g.normal(0.0, sigma, n_groups) if sigma > 0.0 else np.zeros(n_groups)
    x1 = g.normal(0.0, 1.0, n)
    x2 = (g.random(n) < 0.5).astype(float)
    eta = RI_BETA_TRUE * x1 + 0.4 * x2 + np.repeat(b, per_group)
    return pd.DataFrame({
        "x1": x1,
        "x2": x2,
        "group": pd.Series(labels, dtype="string"),
        "y": _ordinal_draw(g, eta),
    })


def flat_ri_frame(n_groups: int = 6, per_group: int = 60, seed: int = RI_SEED) -> pd.DataFrame:
    """`two_group_ri_frame` at `sigma = 0`, so `polr_ri` reaches `POLR_RI_SIGMA_FLOOR`.

    The witness for §12.5: `at_floor` is True, `sigma == C.POLR_RI_SIGMA_FLOOR`, and **the fit does
    not raise** — [§14a] names sigma^2_C = 0 as a legitimate answer, so a floored fit is an answer
    and not a failure. Six groups rather than two, because the point is that the boundary is reached
    when the data put it there and not when the group count leaves sigma unidentified.
    """
    return two_group_ri_frame(n_groups, per_group, sigma=0.0, seed=seed)
