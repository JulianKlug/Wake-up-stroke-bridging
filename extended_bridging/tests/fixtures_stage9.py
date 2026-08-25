"""The constructed frames Stage 9 is tested on — §15.0 of `specs/stage9_secondary_binary_estimators.md`.

**Every fixture here is the specification, not a convenience.** Stage 9's spec pinned quantities on
three of these frames to ten decimal places while the frames themselves appeared nowhere, so
`tau = 0.2252774207` was an instruction to reproduce a number from a construction that did not exist
(§22.3 item 1). This module is that defect repaired: every pin in §15 is measured against the code
below, and a pin nobody can reproduce is not a pin.

Five constructions and four constants:

    golden_frame()        24 records, EVERY VALUE A LITERAL, no RNG           §15.0.2
    tilt_frame()          92 records on which the h- and w-tilts differ       §15.0.4
    separated_frame()     40 records, perfectly separated on treatment        §15.0.4a
    dr_population()       the three arms of the double-robustness design      §15.0.5
    ato()                 the [§7] ATO indexed by whichever score it is given §15.0.5
    renamed_replicate()   one stratified replicate, every row renamed         §15.0.6

`separated_frame` COLLIDES BY NAME with `test_model.py`'s, and the collision is deliberate rather than
an oversight: the two exist for the same reason one fitter apart — Stage 8's is a separated ORDINAL
frame for `polr`, this one a separated BINARY frame for `firth`. Importers must alias. `test_model.py`
imports this one as `separated_binary_frame` for exactly that reason.

The degenerate frames of §15.0.3 are deliberately NOT here: each is a handful of rows and naming them
centrally would hide which branch of §6.2 each reaches, so they are built inline in `test_outcome.py`.
What §15.0.3 does specify — and what is not free — is the ARM WEIGHT TOTALS on two of them, and those
are stated in the tests that build them.

Nothing in this file reads the workbook, and no patient-derived value appears in it [§4.3].
"""
from __future__ import annotations

from typing import Final

import numpy as np
import pandas as pd

import config as C


# --- 15.0.2  the golden frame ---------------------------------------------------------------------

def golden_frame() -> pd.DataFrame:
    """24 records, every value a literal, no RNG. The Stage 9 regression pin [§15.0.2].

    Three centres and NOT four: USZ is absent, so `model.design` drops `center_USZ` exactly as it
    does on the workbook (§3.3), and the fixture exercises the constant-column rule rather than
    stepping around it. `age` is the one continuous covariate and `atrial_fib` the one binary one,
    which is the minimum that makes the two covariate lists of §15.0.2's pair differ by a column.

    `e` is WRITTEN OUT rather than fitted. The frame therefore carries no propensity model, and the
    golden vector is a property of this document rather than of a numpy version's random stream
    [Stage 7 §12.0.2]. `w` is [§7]'s overlap weight, `1 - e` treated and `e` control, computed by
    the caller from this column.

    `age` is deliberately near-separating on `y` -- 80/76 on the events, 53/57 on the non-events --
    because the pair's whole purpose is to show a six-parameter model over-fitting a minority cell
    of 11 (§15.0.2's closing paragraphs). A frame on which the two fits agreed would pin the
    arithmetic and demonstrate nothing.
    """
    return pd.DataFrame({
        "case_id":    [f"G{i:02d}" for i in range(1, 25)],
        "center":     ["HUG"] * 9 + ["CHUV"] * 8 + ["Lugano"] * 7,
        "age":        [80.0, 53.0, 57.0, 76.0, 53.0, 80.0, 76.0, 53.0, 57.0,
                       80.0, 53.0, 80.0, 76.0, 53.0, 57.0, 76.0, 80.0,
                       53.0, 76.0, 53.0, 57.0, 53.0, 80.0, 53.0],
        "atrial_fib": [0.0, 1.0, 0.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0,
                       1.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 0.0,
                       0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 1.0],
        C.TREATMENT:  [1.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0,
                       1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 0.0, 1.0,
                       0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0],
        "y":          [1.0, 0.0, 0.0, 1.0, 0.0, 1.0, 1.0, 0.0, 0.0,
                       1.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0,
                       0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0],
        "e":          [0.62, 0.55, 0.41, 0.70, 0.33, 0.58, 0.46, 0.66, 0.29,
                       0.74, 0.37, 0.61, 0.52, 0.44, 0.26, 0.49, 0.68,
                       0.31, 0.72, 0.39, 0.64, 0.22, 0.57, 0.35],
    })


def golden_arrays(frame: pd.DataFrame | None = None):
    """(y, a, w, e) from `golden_frame()`, with `w` the [§7] overlap weight built from `e`.

    Written here rather than in each test for one reason: `w` is `1 - e` in the TREATED arm and `e`
    in the comparator, and a test that builds it the other way round measures a different estimand
    while every assertion in it still reads plausibly. One construction, one place to be wrong.
    """
    df = golden_frame() if frame is None else frame
    a = df[C.TREATMENT].to_numpy(dtype=float)
    e = df["e"].to_numpy(dtype=float)
    w = np.where(a == 1.0, 1.0 - e, e)
    return df["y"].to_numpy(dtype=float), a, w, e


# --- 15.0.4  the tilt-separation frame ------------------------------------------------------------

TILT_SEED: Final[int] = 20260824


def tilt_frame(seed: int = TILT_SEED, n: int = 92,
               x_scale: float = 3.0, interaction: float = 3.0):
    """The §15.8 item 3 fixture: 92 records on which the h- and w-tilts measurably differ.

    Returns (x, e, a, m1, m0, y). `m1` and `m0` are the ORACLE conditional means and not a fit,
    because the test pins the tilting function and a nuisance-estimation error between the two calls
    would be indistinguishable from a tilt difference.

    `x_scale = 3.0` and `interaction = 3.0` are the `wide x, interaction 3.0` shape of §8.4's
    six-candidate table -- the row with the largest measured gap that was not the non-monotone
    outlier. The seed is pinned because the gap is a property of the draw as well as of the
    construction: §8.4's table shows one construction with the LARGEST sd(d) producing the
    SECOND-SMALLEST gap, so "make it heterogeneous" is not a specification.

    `e` is clipped into [0.05, 0.95] so that h = e(1-e) is bounded away from zero and Sh does not
    collapse onto a handful of records -- which would make the gap a property of those records.
    """
    g = np.random.default_rng(seed)
    x = g.normal(0.0, x_scale, n)
    e = np.clip(1.0 / (1.0 + np.exp(-(0.8 * x))), 0.05, 0.95)
    a = (g.random(n) < e).astype(float)
    m1 = 1.0 / (1.0 + np.exp(-(-0.3 + 0.5 * x + interaction * 0.5)))
    m0 = 1.0 / (1.0 + np.exp(-(-0.3 + 0.5 * x - interaction * 0.5)))
    y = np.where(a == 1.0, (g.random(n) < m1).astype(float),
                 (g.random(n) < m0).astype(float))
    return x, e, a, m1, m0, y


# --- 15.0.4a  the separated frame -----------------------------------------------------------------

def separated_frame(n: int = 40) -> pd.DataFrame:
    """40 records, PERFECTLY separated on treatment: y == C.TREATMENT on every row. §9.6, §15.7.

    `center` alternates independently of treatment and `atrial_fib` runs on a period of four, so the
    design is FULL RANK. A first draft made `center` track the arms -- HUG control, CHUV treated --
    and `firth` raised F3 "the design has rank 3 of 4 columns" rather than converging, which is a
    collinearity failure and not the separation this fixture exists to exhibit.

    Only two of the four declared centres appear, so `design` drops `center_Lugano` AND `center_USZ`
    and the frame exercises the constant-column rule twice.
    """
    return pd.DataFrame({
        "case_id":    [f"S{i:02d}" for i in range(n)],
        "center":     ["HUG", "CHUV"] * (n // 2),
        "atrial_fib": [0.0, 0.0, 1.0, 1.0] * (n // 4),
        C.TREATMENT:  [0.0] * (n // 2) + [1.0] * (n // 2),
        "y":          [0.0] * (n // 2) + [1.0] * (n // 2),
    })


# --- 15.0.5  the double-robustness population -----------------------------------------------------

DR_SEED: Final[int] = 20260825
DR_SEEDS_12: Final[tuple[int, ...]] = tuple(20260825 + 7 * i for i in range(12))
DELTA: Final[float] = 0.12          # the constant RISK DIFFERENCE of the companion arm


def dr_population(n: int, seed: int, arm: str):
    """§15.10's population. `arm` is "heterogeneous" | "constant" | "no_interaction".

    Returns (e_true, e_mis, a, m1, m0, y) with m1, m0 the ORACLE conditional means -- the strongest
    available form of "a correct outcome model", which removes nuisance-estimation error from a
    result about the PROPENSITY model being wrong (§8.5).

    `e_mis` OMITS x2 and `e_true` does not, which is the misspecification. The three arms differ in
    how tau(X) = m1 - m0 varies, and that is the whole design:

      * "heterogeneous"  -- x2 enters m1 only, so tau(X) varies strongly. sd(tau(X)) 0.2435.
      * "constant"       -- m1 = m0 + DELTA exactly, so tau(X) is CONSTANT ON THE RISK-DIFFERENCE
                            SCALE. m0 is bounded into (0.25, 0.70) so m1 stays a probability;
                            without the bound the difference truncates at the top and the arm is no
                            longer constant where it matters. sd(tau(X)) 0.0000.
      * "no_interaction" -- a constant LOG-ODDS effect and therefore a NON-constant risk difference.
                            This is the companion a reasonable draft would have written and it is
                            DETECTABLY BIASED at |t| = 21.3 (§8.5). It exists to be asserted biased.

    The three arms share x1, x2 and the treatment draw for a given seed, because they differ in the
    outcome model and not in the population.
    """
    g = np.random.default_rng(seed)
    x1, x2 = g.normal(0, 1, n), g.normal(0, 1, n)
    e_true = 1.0 / (1.0 + np.exp(-(0.5 * x1 + 1.5 * x2)))
    e_mis = 1.0 / (1.0 + np.exp(-(0.5 * x1)))
    a = (g.random(n) < e_true).astype(float)
    if arm == "heterogeneous":
        m1 = 1.0 / (1.0 + np.exp(-(-0.5 + 0.6 * x1 + 1.6 + 1.5 * x2)))
        m0 = 1.0 / (1.0 + np.exp(-(-0.5 + 0.6 * x1)))
    elif arm == "constant":
        m0 = np.clip(1.0 / (1.0 + np.exp(-(-0.5 + 0.6 * x1 + 0.5 * x2))), 0.25, 0.70)
        m1 = m0 + DELTA
    else:
        m1 = 1.0 / (1.0 + np.exp(-(-0.5 + 0.6 * x1 + 0.5 * x2 + 1.6)))
        m0 = 1.0 / (1.0 + np.exp(-(-0.5 + 0.6 * x1 + 0.5 * x2)))
    y = np.where(a == 1.0, (g.random(n) < m1).astype(float),
                 (g.random(n) < m0).astype(float))
    return e_true, e_mis, a, m1, m0, y


def ato(e: np.ndarray, m1: np.ndarray, m0: np.ndarray) -> float:
    """The [§7] ATO indexed by whichever score it is given -- which is §8.3's whole point."""
    h = e * (1.0 - e)
    return float(np.sum(h * (m1 - m0)) / np.sum(h))


# --- 15.0.6  the renamed replicate ----------------------------------------------------------------

def renamed_replicate(df: pd.DataFrame, seed: int,
                      stratum: str = "center") -> pd.DataFrame:
    """One stratified patient-level replicate, with every drawn row given a distinct case_id.

    NOT A FIX and not proposed as one. `propensity._record_exclusion` raises SchemaError on 26.0% of
    stratified replicates because a resample contains duplicate case_ids by construction (§12.3),
    and choosing between the two candidate repairs is a decision about what a replicate's log means
    -- which is [§10]'s subject. This helper adopts the renaming candidate FOR MEASUREMENT ONLY, and
    it lives in tests/ rather than in a shipped module for exactly that reason.

    It is here because four rows of §21 rest on it. A measurement whose apparatus is not in the
    repository is not a measurement anyone can check, and an earlier draft left it in a scratchpad.

    If Stage 10 adopts the OTHER candidate -- teaching the guard that a replicate's rows are distinct
    draws -- these numbers should be unchanged, because renaming touches only the audit entry and not
    e, w or any mask. That is an argument and §21 labels it as one.
    """
    g = np.random.default_rng(seed)
    out = []
    for _, rows in df.groupby(stratum, sort=True):
        take = rows.iloc[g.integers(0, len(rows), len(rows))].copy()
        out.append(take)
    rep = pd.concat(out, ignore_index=True)
    rep["case_id"] = [f"{cid}#{i}" for i, cid in enumerate(rep["case_id"])]
    return rep
