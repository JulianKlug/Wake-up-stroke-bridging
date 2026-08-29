"""The constructed frames Stage 13 is tested on — §16.0 of `specs/stage13_feasible_policy.md`.

**Every fixture here is the specification, not a convenience**, `fixtures_stage12.py`'s reason carried
one stage on. Every pin in Stage 13 §16 is measured against the code below.

Five functions and two constants (§16.0):

    mixed_frame()                     four declared centres, one with no treated patient,
                                      contraindicated patients at TWO centres, none treated,
                                      all seven mRS levels occupied                          §4.1
    all_eligible_frame()              `mixed_frame` with every patient eligible: `population`
                                      returns it with the no-contraindicated sentence         §11
    treated_contraindicated_frame()   one contraindicated patient with the treated code; U3   §4.1
    separated_indicator_frame()       every contraindicated patient at mRS 6, so `delta`
                                      separates while `beta` does not                         §9.1
    no_contraindicated_frame()        `mixed_frame` minus its contraindicated rows — for the
                                      replicate body and `contrast` ONLY, never `population`  §5.3

    POLICY_SEED  N_CONTRAINDICATED

Nothing here reads the workbook, and no patient-derived value appears in it [§4.3].

**Every frame declares all four centres** so that `bootstrap.resample` under `C.BOOT_STRATUM` draws
from four strata, as the workbook does. **Every frame carries `C.ELIGIBILITY` with BOTH declared
values** where the workbook does: the contraindicated rows are `C.INELIGIBLE`, which is the one thing
that makes this a [§14b] frame rather than a [§14a] one.
"""
from __future__ import annotations

from typing import Final

import numpy as np
import pandas as pd

import config as C
from fixtures_stage12 import _covariates, _ordinal_draw

POLICY_SEED: Final[int] = 20260829

# `mixed_frame`'s per-centre sizes, treated counts and contraindicated counts. The SHAPE is the
# workbook's own and the numbers are not: one centre with no treated patient, contraindicated
# patients at exactly two centres — the first and the last in CENTER_ORDER, as on the workbook — and
# none of them treated (§4.1).
_CENTRES: Final[tuple[tuple[str, int, int, int], ...]] = (
    ("HUG", 40, 20, 8),         # centre, records, treated, contraindicated
    ("CHUV", 24, 8, 0),
    ("Lugano", 28, 4, 0),
    ("USZ", 16, 0, 6),          # the never-IVT centre, and the second contraindicated stratum
)
N_CONTRAINDICATED: Final[int] = sum(k for _, _, _, k in _CENTRES)

# The contraindicated rows' shift on the `logit P(Y <= k)` scale. NEGATIVE, in the direction [§3]
# predicts — a contraindicated patient has a worse prognosis — and large enough that the indicator
# carries a coefficient the fit table can show, without separating the fit.
_CONTRAINDICATION_EFFECT: Final[float] = -1.2
_TREATMENT_EFFECT: Final[float] = 0.7
_CENTRE_EFFECT: Final[tuple[float, ...]] = (0.6, -0.4, 0.3, -0.5)


def mixed_frame(seed: int = POLICY_SEED) -> pd.DataFrame:
    """A classified population over four declared centres with contraindicated patients at two.

    108 records, 32 treated, 14 contraindicated — all EVT alone, which is Stage 4's load-bearing
    assertion and what U3 re-checks. All seven mRS levels occupied, asserted by the tests rather
    than assumed here, because the mRS draw is stochastic.

    Within each centre's block the treated rows come first and the contraindicated rows last, and
    the two never overlap: treated rows are eligible by construction, contraindicated rows are
    never treated.
    """
    g = np.random.default_rng(seed)
    centres = np.concatenate([[name] * n for name, n, _, _ in _CENTRES])
    treated = np.concatenate([
        np.concatenate([np.ones(t), np.zeros(n - t)]) for _, n, t, _ in _CENTRES])
    contraindicated = np.concatenate([
        np.concatenate([np.zeros(n - k), np.ones(k)]) for _, n, _, k in _CENTRES])
    n = len(centres)
    columns = _covariates(g, n)
    eta = (_TREATMENT_EFFECT * treated
           + _CONTRAINDICATION_EFFECT * contraindicated
           + 0.02 * (columns["age"] - 70.0)
           - 0.04 * (columns["nihss_baseline"] - 14.0)
           + np.repeat(_CENTRE_EFFECT, [n for _, n, _, _ in _CENTRES]))
    eligibility = np.where(contraindicated == 1.0, C.INELIGIBLE, C.ELIGIBLE)
    return pd.DataFrame({
        "case_id": pd.Series([f"S13-{i:03d}" for i in range(n)], dtype="string"),
        "center": pd.Series(centres, dtype="string"),
        **columns,
        C.TREATMENT: treated,
        C.PRIMARY_OUTCOME: _ordinal_draw(g, eta),
        C.ELIGIBILITY: pd.Series(eligibility, dtype="string"),
    })


def _contraindicated_rows(df: pd.DataFrame) -> pd.Series:
    return df[C.ELIGIBILITY] == C.INELIGIBLE


def all_eligible_frame(seed: int = POLICY_SEED) -> pd.DataFrame:
    """`mixed_frame` with every patient eligible. The workbook-without-contraindicated case (§11).

    The rows are kept and only the class changes, so the record count and the arm sizes are
    `mixed_frame`'s exactly: the difference between the two frames is the indicator and nothing else.
    """
    df = mixed_frame(seed).copy()
    df[C.ELIGIBILITY] = pd.Series([C.ELIGIBLE] * len(df), dtype="string")
    return df


def treated_contraindicated_frame(seed: int = POLICY_SEED) -> pd.DataFrame:
    """`mixed_frame` with ONE contraindicated patient carrying the treated code. U3's witness."""
    df = mixed_frame(seed).copy()
    first = df.index[_contraindicated_rows(df)][0]
    df.loc[first, C.TREATMENT] = float(max(C.TREATMENT_LABELS))
    return df


def separated_indicator_frame(seed: int = POLICY_SEED) -> pd.DataFrame:
    """`mixed_frame` with every contraindicated patient at the top mRS level. §9.1 on a fixture.

    The indicator then perfectly predicts `Y = 6` among the rows it indexes: `delta` runs past
    `POLR_MAX_ABS_BETA`, `polr` converges, and `beta` — identified from the eligible rows alone — does
    not move. U8 must not fire.
    """
    df = mixed_frame(seed).copy()
    df.loc[_contraindicated_rows(df), C.PRIMARY_OUTCOME] = float(C.MRS_LEVELS[-1])
    return df


def no_contraindicated_frame(seed: int = POLICY_SEED) -> pd.DataFrame:
    """`mixed_frame` minus its contraindicated rows, with the CONTRAINDICATED column already present.

    **For the replicate body and `contrast` only, never `population`.** It stands in for the §5.3
    draw in which a stratified resample happens to contain no contraindicated patient: `design` drops
    the indicator as constant, `share_eligible` is 1 and `rd == eligible_rd`. The column is written
    here because the body receives frames `population` has already stamped.
    """
    df = mixed_frame(seed)
    out = df[~_contraindicated_rows(df)].reset_index(drop=True).copy()
    out[C.CONTRAINDICATED] = 0.0
    return out
