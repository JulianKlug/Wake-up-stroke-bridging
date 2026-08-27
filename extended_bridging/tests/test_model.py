"""Acceptance tests for Stage 6a — §12.1-12.7 and §12.12 of `specs/stage6_propensity_and_weights.md`.

`model.py` is the outcome-agnostic half of Stage 6, so almost everything here runs on synthetic arrays
with known answers: the numerics are testable **only** on synthetic data, and the integration **only**
on a real cohort, which is `test_propensity.py`'s business.

Tests needing the private workbook are marked `skipif(not DATA_XLSX.exists())` and take this file's
**own** module-scoped `workbook` fixture. Everything else runs on a plain checkout with no `data/` and
with the `reference` dependency group not installed.

This file is **not** exempt from the Stage 1 §7 raw-name scan and must not become exempt.

**The three silent failures this stage exists to make loud**, each tested with a companion that shows
what happens *without* the guard — because all three return a plausible number otherwise, which is what
separates Stage 6 from the five before it (they failed loudly by default)::

    a nan in the design      F4   without it: all-nan beta, all-nan p, a NORMAL return  §12.3
    an unlisted factor level D3   without it: a finite propensity, computed as though   §12.6
                                  the record were at the reference centre
    a missing factor value   D4   without it: an all-zero dummy row BYTE-IDENTICAL to a §12.6
                                  real reference-level record — F4, D3 and F5 all blind
    an ordinal response      F6   without it: fits, converges, returns coefficients for §12.6
                                  a logistic model of a seven-level ordinal

**The fixtures split into three families** (§12.0), and the second exists because four of §5.2's
safeguards are reached by nothing natural — measured at zero rescales and zero halvings on the workbook
design and on `separated()`. Untested, their first execution would be inside Stage 10's `N_BOOT` loop,
where a defect reads as a dropped-replicate count rather than as a bug.
"""
from __future__ import annotations

import ast
import dataclasses
import os
import warnings
from pathlib import Path
from typing import Final

import numpy as np
import pandas as pd
import pytest

import config
import model

DATA_GATED = pytest.mark.skipif(
    not config.DATA_XLSX.exists(),
    reason="the private workbook is gitignored and absent from this checkout")

MODULE = Path(model.__file__).resolve()
MODULE_DIR = MODULE.parent
SOURCE = MODULE.read_text(encoding="utf-8")


def message_of(excinfo) -> str:
    return str(excinfo.value)


# --- 12.0  the frames and fixtures this stage is tested on ------------------------------------------
#
# Written out in full rather than by signature: a fixture given only a name makes every number §18 pins
# from it unreproducible, and the spec's standing rule is that nothing absent from it is to be
# invented. Both seed from `config.SEED`, never from a literal, so a fixture cannot drift from the run
# summary's recorded seed.


def well_behaved(n=4000, seed=config.SEED):
    """§12.4 — Firth ≈ MLE at large n. Three standard-normal covariates, no near-separation."""
    rng = np.random.default_rng(seed)
    X = rng.normal(size=(n, 3))
    eta = 0.3 + X @ np.array([0.5, -0.4, 0.2])
    y = (rng.random(n) < 1.0 / (1.0 + np.exp(-eta))).astype(float)
    return pd.DataFrame(X, columns=["x1", "x2", "x3"]), y


def separated():
    """§12.5 — complete separation. The pilots' benchmark, unchanged."""
    return pd.DataFrame({"x": np.arange(1.0, 11.0)}), np.array([0] * 5 + [1] * 5, dtype=float)


def collinear(n=200, seed=config.SEED):
    """§12.6 — F3. well_behaved's design with x3 replaced by an exact copy of x1.

    n=200 rather than 4000 because F3 raises before the first iteration and a 4000-row design buys
    nothing but runtime.
    """
    X, y = well_behaved(n, seed)
    X = X.copy()
    X["x3"] = X["x1"]
    return X, y


# --- 12.0 1b  the hard-fit family --------------------------------------------------------------------


def big_first_step():
    """§12.4a — the TRUST REGION. A tiny-scale separated design, so the first Newton step is huge.

    Measured under §5.2a's relative radius: ‖step₁‖ = 48.3 against FIRTH_MAX_STEP = 5.0, 2 rescales,
    converging in 9 iterations. The scale is what does it: the coefficient that separates these four
    points is ~1/0.05, so the unbounded Newton step at beta = 0 overshoots by an order of magnitude.

    This is also the fixture §12.4a multiplies by 1e-3 … 1e6 to assert scale invariance.
    """
    return pd.DataFrame({"x": [-0.05, -0.02, 0.02, 0.05]}), np.array([0, 0, 1, 1], dtype=float)


def needs_halving():
    """§12.4a — STEP-HALVING. Two orthogonal separating columns over four records.

    Measured: 1 halving, converging in 25 iterations, 0 trust-region rescales — so this fixture
    isolates halving from the trust region rather than firing both at once. Under
    FIRTH_MAX_HALVINGS = 1 it raises at iteration 2, which is the other half of §12.4a.
    """
    return (pd.DataFrame({"x1": [-1.0, -1.0, 1.0, 1.0], "x2": [-1.0, 1.0, -1.0, 1.0]}),
            np.array([0, 0, 1, 1], dtype=float))


def travels_far():
    """§12.4a — the fixture §5.2a is ABOUT. Scale pushed a further 20x down.

    Its Firth optimum is finite and perfectly identified: beta ~ (0, 2007.9), penalised
    log-likelihood -8.6292176, confirmed against scipy.optimize on the same objective.

    Under an ABSOLUTE trust radius it fails at FIRTH_MAX_ITER, because
    FIRTH_MAX_STEP * FIRTH_MAX_ITER = 5.0 * 200 = 1000 caps total coefficient travel below
    ‖beta‖ ~ 2008 — a fit reported as a failure because of the covariate's units, not its data
    [§5.2a]. Under §5.2a's RELATIVE radius it converges at iteration 10, on the score route.
    """
    return pd.DataFrame({"x": [-1e-3, -5e-4, 5e-4, 1e-3]}), np.array([0, 0, 1, 1], dtype=float)


HARD_FITS = {"big_first_step": big_first_step, "needs_halving": needs_halving,
             "travels_far": travels_far}


def scaled(fixture, k):
    """`fixture`'s design with every covariate multiplied by k. A Firth estimate is equivariant."""
    X, y = fixture()
    return X * k, y


# --- frames, for the design-matrix half ---------------------------------------------------------------
#
# `design` takes a frame rather than an array, so these are the smallest frames that carry a declared
# factor. They deliberately carry NO `case_id`: D3 must name records when it can and must not require
# an identifier to run at all, and a frame without one is exactly what Stage 9 or Stage 12 may assemble.

LINEAR = ("age", "nihss_baseline")
WITH_FACTORS = ("age", "onset_type", "center")


def frame(centres=("HUG", "CHUV"), onsets=("witnessed", "wake_up"), **overrides):
    """A minimal frame over the declared factors, one record per `centres` entry."""
    df = pd.DataFrame({
        "case_id": [f"F-{i}" for i in range(len(centres))],
        "age": [70.0 + i for i in range(len(centres))],
        "nihss_baseline": [10.0 + i for i in range(len(centres))],
        "core_ml": [5.0 + i for i in range(len(centres))],
        "center": pd.Series(list(centres), dtype="string"),
        "onset_type": pd.Series([onsets[i % len(onsets)] for i in range(len(centres))],
                                dtype="string"),
    })
    for column, value in overrides.items():
        df[column] = value
    return df


# --- 12.1  `design` reference-codes on the declared levels --------------------------------------------

def test_an_absent_declared_level_becomes_a_column_and_is_then_dropped():
    # The declared level set is what makes this possible at all: without the Categorical conversion
    # `get_dummies` emits only the levels present, so an absent level's all-zero column never exists,
    # the constant-column rule has nothing to do, and the width is a property of the sample (§12.2).
    df = frame(centres=("HUG", "CHUV", "Lugano"), onsets=("witnessed", "wake_up", "unwitnessed"))
    X, dropped = model.design(df, WITH_FACTORS)
    assert "center_USZ" in dropped
    assert "center_USZ" not in X.columns


def test_the_reference_columns_are_absent_and_are_not_reported_as_dropped():
    # They were removed as REFERENCES, not as constants, and §7.2's table distinguishes the two — which
    # is why `_design_table` ranges over the declaration rather than over `X.columns + dropped`.
    df = frame(centres=("HUG", "CHUV", "Lugano"), onsets=("witnessed", "wake_up", "unwitnessed"))
    X, dropped = model.design(df, WITH_FACTORS)
    for reference in ("center_HUG", "onset_type_witnessed"):
        assert reference not in X.columns
        assert reference not in dropped


@pytest.mark.parametrize("covariates,expected", [
    (config.PS_COVARIATES, [
        "age", "sex", "prestroke_mrs", "nihss_baseline", "core_ml", "tmax6_ml", "atrial_fib",
        "onset_type_unwitnessed", "onset_type_wake_up",
        "center_CHUV", "center_Lugano", "center_USZ"]),
])
def test_the_column_order_is_the_measured_one_and_not_the_covariate_order(covariates, expected):
    """`pd.get_dummies(sub, columns=factors)` DROPS the original columns and APPENDS the dummies.

    It does not expand a factor in place. So the order is: the non-factor covariates in `covariates`'
    order, then each factor's dummies in `factors` order, each expanded in FACTOR_LEVELS order with the
    reference removed. `onset_type` is PS_COVARIATES' FIFTH entry and its dummies come after
    `atrial_fib`, which is the seventh.

    This matters beyond the test: `beta`'s order IS this order, `Fit.columns` is what maps a
    coefficient to a name, and the audit log is the only place the coefficients are ever written down
    (§7.3) — so a wrong order here is twelve correctly-computed numbers under eleven wrong labels.
    """
    df = pd.DataFrame({
        "case_id": ["A", "B", "C", "D"],
        "age": [70.0, 71.0, 72.0, 73.0], "sex": [0.0, 1.0, 0.0, 1.0],
        "prestroke_mrs": [0.0, 1.0, 0.0, 1.0], "nihss_baseline": [10.0, 11.0, 12.0, 13.0],
        "core_ml": [5.0, 6.0, 7.0, 8.0], "tmax6_ml": [50.0, 60.0, 70.0, 80.0],
        "atrial_fib": [0.0, 1.0, 0.0, 1.0],
        "onset_type": pd.Series(["witnessed", "unwitnessed", "wake_up", "witnessed"],
                                dtype="string"),
        "center": pd.Series(["HUG", "CHUV", "Lugano", "HUG"], dtype="string"),
    })
    # Before the constant drop, which is what `expected` lists: build it and re-add what was dropped.
    X, dropped = model.design(df, covariates)
    assert list(X.columns) + list(dropped) == [c for c in expected if c in
                                               list(X.columns) + list(dropped)]
    assert list(X.columns) == [c for c in expected if c not in dropped]


def test_a_single_centre_frame_yields_no_centre_column_at_all():
    # Two levels absent and the third constant, so the factor contributes nothing — which is correct,
    # because a single-level factor is collinear with the intercept. This is §4.3's replicate case.
    df = frame(centres=("HUG", "HUG"), onsets=("witnessed", "wake_up"))
    X, dropped = model.design(df, WITH_FACTORS)
    assert not [c for c in X.columns if c.startswith("center_")]
    assert {c for c in dropped if c.startswith("center_")} == {
        "center_CHUV", "center_Lugano", "center_USZ"}


def test_a_covariate_list_with_no_factor_at_all_is_a_no_op():
    # Asserted because Stage 9's outcome_model_covariates registry may yet hold an all-linear list, and
    # a reader should not have to reason about whether the reference-drop comprehension is safe when
    # `factors` is empty. Measured: get_dummies(columns=[]) is a no-op and .drop(columns=[]) is empty.
    df = frame()
    X, dropped = model.design(df, LINEAR)
    assert list(X.columns) == list(LINEAR)
    assert dropped == ()


def test_a_genuinely_constant_linear_covariate_is_dropped():
    df = frame(centres=("HUG", "CHUV"), age=70.0)
    X, dropped = model.design(df, LINEAR)
    assert "age" in dropped
    assert "age" not in X.columns


# --- 12.2  `drop_first` and `remove_unused_categories` are not used ------------------------------------

def _names_in_source(source: str, names: set[str]) -> list[tuple[int, str]]:
    """Every keyword argument or attribute call in `source` whose name is in `names`."""
    found: list[tuple[int, str]] = []
    for node in ast.walk(ast.parse(source)):
        if isinstance(node, ast.keyword) and node.arg in names:
            found.append((node.lineno, node.arg))
        if isinstance(node, ast.Attribute) and node.attr in names:
            found.append((node.lineno, node.attr))
    return found


FORBIDDEN_PANDAS = {"drop_first", "remove_unused_categories"}


def test_neither_drop_first_nor_remove_unused_categories_appears():
    assert _names_in_source(SOURCE, FORBIDDEN_PANDAS) == []


def test_that_scan_actually_fires():
    # A scan that silently matches nothing would otherwise pass as a green test.
    pasted = "X = pd.get_dummies(sub, columns=f, drop_first=True)\ns = s.cat.remove_unused_categories()\n"
    assert {name for _, name in _names_in_source(pasted, FORBIDDEN_PANDAS)} == FORBIDDEN_PANDAS


@pytest.mark.parametrize("factor", list(config.REFERENCE_LEVELS))
def test_drop_first_is_EQUIVALENT_on_a_declared_categorical(factor):
    """On a DECLARED Categorical the two forms agree, and that is measured rather than assumed.

    It would be natural to argue that `drop_first=True` picks whichever level sorts first among those
    *present*, so a frame with no HUG record would silently rebaseline on CHUV. That is FALSE here, and
    it is false because the level set is declared: `get_dummies` emits *every* declared level —
    including an all-zero `center_HUG` when no HUG patient exists — and `drop_first=True` drops the
    first *declared* category, which is the reference.

    So what by-name actually buys is independence from a coupling nothing else stated:
    `REFERENCE_LEVELS[f] == FACTOR_LEVELS[f][0]`. `test_config.py` now asserts it, and with it in place
    the by-name form is *provably* the same model as `drop_first` rather than accidentally the same
    one. This test is parametrised over both factors to point at that assertion.
    """
    assert config.REFERENCE_LEVELS[factor] == config.FACTOR_LEVELS[factor][0], (
        "test_config.py's coupling assertion is what makes the two forms equivalent")
    df = frame(centres=("CHUV", "Lugano"), onsets=("unwitnessed", "wake_up"))
    sub = df[[factor]].copy()
    sub[factor] = pd.Categorical(sub[factor], categories=config.FACTOR_LEVELS[factor])
    by_position = pd.get_dummies(sub, columns=[factor], dtype=float, drop_first=True)
    by_name = pd.get_dummies(sub, columns=[factor], dtype=float).drop(
        columns=[f"{factor}_{config.REFERENCE_LEVELS[factor]}"])
    assert list(by_position.columns) == list(by_name.columns)


def test_drop_first_on_a_BARE_column_rebaselines_and_that_is_the_pilots_combination():
    # `pilots/analysis.py:40` uses drop_first=True on a column `pilot_data.py` never declared, and it
    # is the COMBINATION that is wrong rather than either half. Measured: on a bare string centre with
    # no HUG record, drop_first drops center_CHUV and leaves center_Lugano alone, so every remaining
    # coefficient means something else.
    df = frame(centres=("CHUV", "Lugano"))
    bare = pd.get_dummies(df[["center"]], columns=["center"], dtype=float, drop_first=True)
    assert list(bare.columns) == ["center_Lugano"]
    assert "center_HUG" not in bare.columns


def test_without_the_categorical_conversion_no_constant_column_exists_to_drop():
    df = frame(centres=("HUG", "CHUV", "Lugano"))
    bare = pd.get_dummies(df[["center"]], columns=["center"], dtype=float)
    assert "center_USZ" not in bare.columns          # the width is a property of the SAMPLE
    X, dropped = model.design(df, ("center",))
    assert "center_USZ" in dropped                    # and of config.py, once declared


def test_a_frame_missing_the_reference_level_does_not_KeyError_and_raises_F3_instead():
    """DoD-10, all three halves. The first is the one an earlier draft of the spec got wrong.

    (a) `design` does NOT raise, because the declared Categorical emits an all-zero `center_HUG` that
        is dropped as the reference — the column is always there;
    (b) `firth` then raises F3 at rank k−1 of k, because the surviving centre dummies sum to 1 on every
        row and are collinear with the intercept. That is ONE dependency, so k−1 and not k−2.
    """
    df = frame(centres=("CHUV", "Lugano", "CHUV", "Lugano"))
    X, dropped = model.design(df, ("center",))        # (a) does not raise
    assert list(X.columns) == ["center_CHUV", "center_Lugano"]
    with pytest.raises(model.FitError) as e:
        model.firth(X, np.array([0.0, 1.0, 0.0, 1.0]))
    assert "F3" in message_of(e)
    assert f"rank {X.shape[1]} of {X.shape[1] + 1}" in message_of(e)     # k−1 of k


# --- 12.3  `complete_cases` is returned, and a missing covariate cannot reach the fit -----------------

def test_the_mask_is_total_boolean_and_flips_for_exactly_the_blanked_record():
    df = frame(centres=("HUG", "CHUV", "Lugano"))
    assert model.complete_cases(df, LINEAR).all()
    df.loc[df["case_id"] == "F-1", "age"] = np.nan
    mask = model.complete_cases(df, LINEAR)
    assert mask.dtype == bool and not mask.isna().any()
    assert list(mask) == [True, False, True]


def test_design_does_not_drop_the_incomplete_row():
    # This is the assertion that `complete_cases` is a separate function rather than a `dropna` inside
    # `design`: a design that silently dropped rows would give Stage 9 a matrix and a response of
    # different lengths, aligned by luck.
    df = frame(centres=("HUG", "CHUV", "Lugano"))
    df.loc[df["case_id"] == "F-1", "age"] = np.nan
    X, _ = model.design(df, LINEAR)
    assert len(X) == len(df)
    assert X["age"].isna().sum() == 1


@pytest.mark.parametrize("dtype,blank", [("float64", np.nan), ("Int64", pd.NA), ("Float64", pd.NA)])
def test_a_nan_in_the_design_raises_F4_naming_the_column(dtype, blank):
    # All three dtypes convert to nan identically under .astype(float) — measured — so `design` is
    # dtype-agnostic and TODOS.md's open CTP-volume dtype question does not block this stage.
    df = frame(centres=("HUG", "CHUV", "Lugano"))
    df["age"] = df["age"].astype(dtype)
    df.loc[df["case_id"] == "F-1", "age"] = blank
    X, _ = model.design(df, LINEAR)
    with pytest.raises(model.FitError) as e:
        model.firth(X, np.array([0.0, 1.0, 0.0]))
    assert "F4" in message_of(e) and "age" in message_of(e)


def test_WITHOUT_F4_the_failure_is_MISDIAGNOSED_as_non_convergence(monkeypatch):
    """The necessity of F4 is MEASURED, not asserted in a comment. This is DoD-6 as a test.

    **And what it measures is not what §3.1 predicted, which is recorded here rather than smoothed
    over.** §3.1 and §18 state that a nan anywhere in `X` makes the loop "terminate normally" and
    return all-nan coefficients with a plausible iteration count. It does not, and the mechanism is
    `_penalised_loglik`'s own `-inf` branch — which is in the spec's canonical code:

        slogdet(a matrix containing nan)  ->  (sign 1.0, logabsdet NAN), measured
        so `np.isfinite(logdet)` is False and `_penalised_loglik` returns -INF, not nan
        `-inf >= -inf` is True, so every candidate step is ACCEPTED with zero halvings
        `moved` is `abs(-inf - -inf)` = nan, and `nan < FIRTH_TOL` is False
        the modified score is nan, and `nan < FIRTH_SCORE_TOL` is False
        so NEITHER convergence route can ever fire, and the loop runs out of iterations

    The harm is therefore different from the one §3.1 describes, and for Stage 10 it is arguably
    worse rather than milder: a missing covariate is reported as a `FitError` saying the estimator did
    not converge — which is exactly what [§10] catches to **drop and count the replicate**. The data
    defect is absorbed into a sparse-replicate failure count and never surfaces as a data defect at
    all. F4 is what turns that into a message naming the column.

    Verified on both kernels: this implementation and `pilots/analysis.py`'s `_firth_coefs` share the
    `-inf` branch and both behave this way.
    """
    monkeypatch.setattr(model, "_assert_fittable", lambda *a, **k: None)
    df = frame(centres=("HUG", "CHUV", "Lugano"))
    df.loc[df["case_id"] == "F-1", "age"] = np.nan
    X, _ = model.design(df, LINEAR)
    with np.errstate(all="ignore"):
        with pytest.raises(model.FitError) as e:
            model.firth(X, np.array([0.0, 1.0, 0.0]))
    text = message_of(e)
    assert "no convergence" in text                  # a message about the ESTIMATOR ...
    assert "age" not in text                         # ... naming nothing about the DATA
    assert "moved nan" in text                       # the tell, if anyone reads it
    # and the mechanism, asserted directly so this test explains itself
    Xc = np.column_stack([np.ones(3), X.to_numpy(dtype=float)])
    assert model._penalised_loglik(Xc, np.array([0.0, 1.0, 0.0]), np.zeros(3)) == -np.inf


def test_a_nan_in_the_response_raises_F6_and_not_the_step_halving_message():
    X, y = separated()
    y = y.copy()
    y[3] = np.nan
    with pytest.raises(model.FitError) as e:
        model.firth(X, y)
    assert "F6" in message_of(e) and "non-finite" in message_of(e)
    # F6's own text *explains* the step-halving message it replaces, so the discriminator is the
    # actual raise, not the words: the exhausted-halvings error names its count.
    assert "exhausted 30 halvings" not in message_of(e)


def test_WITHOUT_F6_a_nan_response_raises_the_WRONG_diagnosis(monkeypatch):
    # What F6 buys is measured to be the DIAGNOSIS rather than the raise. Without it every candidate
    # step is nan, no step is ever accepted, and the fit reports a sentence about [§7]'s single
    # estimator — which a reader follows to Stage 10 and sparse replicates, not to a missing outcome.
    monkeypatch.setattr(model, "_assert_fittable", lambda *a, **k: None)
    X, y = separated()
    y = y.copy()
    y[3] = np.nan
    with pytest.raises(model.FitError) as e:
        model.firth(X, y)
    assert "step-halving exhausted" in message_of(e)
    assert "F6" not in message_of(e)


# --- 12.4  Firth matches the unpenalised MLE on well-behaved data [roadmap] ----------------------------

def test_firth_agrees_with_the_unpenalised_mle_at_large_n():
    """Firth's penalty is O(1) and the likelihood is O(n), so the two agree as n grows.

    A test at n = 100 would not distinguish a correct penalty from a mis-scaled one. The threshold
    stays at 0.01 rather than tightening to the measured 0.000787: it is a statement about the
    PENALTY'S SCALING, and pinning it at the measurement would make it fail on a numpy patch release
    rather than on a mistake.
    """
    import statsmodels.api as sm
    X, y = well_behaved()
    fit = model.firth(X, y)
    mle = sm.Logit(y, sm.add_constant(X.to_numpy(dtype=float))).fit(disp=0)
    assert float(np.max(np.abs(fit.beta - mle.params))) < 0.01
    assert fit.converged_on == "likelihood"
    assert fit.iterations == 4


def test_step_halving_is_on_the_PENALISED_likelihood():
    # Halving on l(b) would drive the fit toward the boundary the penalty exists to hold it back from —
    # and it would do so WHILE CONVERGING, which is why this is asserted rather than reviewed. On
    # separated data the unpenalised likelihood is monotone in the slope, so its optimum diverges;
    # Firth's stays finite.
    X, y = separated()
    fit = model.firth(X, y)
    eta = np.clip(np.column_stack([np.ones(len(y)), X.to_numpy(float)]) @ fit.beta, -500, 500)
    unpenalised = float(np.sum(y * eta - np.logaddexp(0.0, eta)))
    bigger = fit.beta * 10.0
    eta2 = np.clip(np.column_stack([np.ones(len(y)), X.to_numpy(float)]) @ bigger, -500, 500)
    assert float(np.sum(y * eta2 - np.logaddexp(0.0, eta2))) > unpenalised     # l(b) still climbing
    assert model._penalised_loglik(                                            # l*(b) does not
        np.column_stack([np.ones(len(y)), X.to_numpy(float)]), y, bigger) < model._penalised_loglik(
        np.column_stack([np.ones(len(y)), X.to_numpy(float)]), y, fit.beta)


# --- 12.4a  the four safeguards that only fire on hard data [§5.2, §5.6] --------------------------------

def test_the_trust_region_fires_on_a_huge_first_step():
    fit = model.firth(*big_first_step())
    assert fit.first_step_norm > config.FIRTH_MAX_STEP
    assert fit.rescales >= 1
    assert fit.iterations == 9 and fit.rescales == 2
    assert np.all((fit.p > 0.0) & (fit.p < 1.0))


def test_step_halving_fires_in_ISOLATION_from_the_trust_region():
    fit = model.firth(*needs_halving())
    assert fit.halvings >= 1
    assert fit.rescales == 0            # isolated: this fixture drives halving and nothing else
    assert fit.iterations == 25 and fit.halvings == 1


def test_neither_safeguard_fires_on_the_well_behaved_or_separated_fixtures():
    for fixture in (well_behaved, separated):
        fit = model.firth(*fixture())
        assert (fit.rescales, fit.halvings) == (0, 0)


def test_FIRTH_MAX_ITER_exhausted_reports_a_NON_ZERO_movement(monkeypatch):
    """§5.2b's correction, and the thing to watch: the earlier form printed 0 on EVERY input.

    The loop's last statement is `ll_old = ll_new`, so falling out of the iteration range reached the
    final raise with the two equal — and the message read "moved 0 against a tolerance of 1e-08" inside
    a *non*-convergence error, stating a movement below the tolerance it claims was not met. This
    assertion could not have failed.

    `travels_far()` CONVERGES under §5.2a's relative radius, which is §5.2a's whole point, so the cap
    is reached by monkeypatching the constant — and the constant IS the subject of the test.
    """
    monkeypatch.setattr(config, "FIRTH_MAX_ITER", 5)
    with pytest.raises(model.FitError) as e:
        model.firth(*travels_far())
    text = message_of(e)
    assert "no convergence in 5 iterations" in text
    moved = float(text.split("moved ")[1].split(" ")[0])
    assert moved > 0.0
    assert "score component" in text and "shortened by the trust region" in text


def test_FIRTH_MAX_HALVINGS_exhausted_says_there_is_no_second_estimator(monkeypatch):
    monkeypatch.setattr(config, "FIRTH_MAX_HALVINGS", 1)
    with pytest.raises(model.FitError) as e:
        model.firth(*needs_halving())
    assert "halvings at iteration 2" in message_of(e)
    assert "there is no second one to try" in message_of(e)      # invariant 5


@pytest.mark.parametrize("k", [1e-3, 1.0, 1e3, 1e6])
def test_the_trust_region_is_SCALE_INVARIANT(k):
    """§5.2a's whole point, asserted directly: the same fixture converges at every scale.

    A Firth estimate is equivariant under rescaling a covariate — replace x by 1000x and beta is the
    old one divided by 1000, with identical fitted probabilities. An ABSOLUTE step bound is not: under
    it the 1e-3 case fails at FIRTH_MAX_ITER and the 1e-1 case takes 84 iterations. This is what
    catches a reversion to `if ‖step‖ > FIRTH_MAX_STEP`.
    """
    fit = model.firth(*scaled(big_first_step, k))
    assert fit.iterations <= config.FIRTH_MAX_ITER
    assert np.all((fit.p > 0.0) & (fit.p < 1.0))


def test_the_scaled_fits_agree_on_the_probabilities_to_1e_3_and_no_further():
    """Not a defect and the tolerance is not slack.

    The optimum for the small-scale design sits at ‖beta‖ ≈ 2008 on a surface flat enough that the
    1e-3 case converges on the SCORE route and the others on the likelihood, so the four stop at
    genuinely different points on the same ridge. Asserting 1e-8 here — as an earlier draft did — would
    be asserting that a flat surface has a sharp optimum.
    """
    ps = [model.firth(*scaled(big_first_step, k)).p for k in (1e-3, 1.0, 1e3, 1e6)]
    spread = max(float(np.max(np.abs(p - ps[0]))) for p in ps)
    assert spread < 1e-3
    assert spread > 1e-8


@pytest.mark.parametrize("name,fixture,route", [
    ("separated", separated, "likelihood"),
    ("well_behaved", well_behaved, "likelihood"),
    ("big_first_step", big_first_step, "likelihood"),
    ("needs_halving", needs_halving, "likelihood"),
    ("travels_far", travels_far, "score"),
    ("separated x1e-3", lambda: scaled(separated, 1e-3), "score"),
    ("big_first_step x1e-3", lambda: scaled(big_first_step, 1e-3), "score"),
])
def test_converged_on_is_asserted_per_fixture_and_BOTH_routes_fire(name, fixture, route):
    """The split IS the assertion, and neither route can quietly stop firing.

    The mechanism, because a reader who sees "score" in a Stage 10 replicate needs to know what it
    means: near the optimum Δl* ≈ ½·sᵀI⁻¹s ~ ½·s²·‖I⁻¹‖. An earlier draft assumed ‖I⁻¹‖ ~ 1 and
    concluded the score route was dead. At ‖I⁻¹‖ = 1e6, |s| = 1e-6 gives Δl* ~ 5e-7, comfortably ABOVE
    FIRTH_TOL — so the likelihood is still moving when the score has gone flat. A nearly-singular
    information matrix and a flat penalised surface are the same fact stated twice, which is precisely
    the condition [§7] chose Firth for. "score" means NEAR-SINGULAR DESIGN, not bug.
    """
    assert model.firth(*fixture()).converged_on == route


@pytest.mark.parametrize("name,fixture,ill_conditioned", [
    ("separated", separated, False),
    ("big_first_step", big_first_step, False),
    ("travels_far", travels_far, True),
    ("separated x1e-3", lambda: scaled(separated, 1e-3), True),
])
def test_the_route_follows_the_conditioning_of_the_information_matrix(name, fixture,
                                                                     ill_conditioned):
    """What makes the split above a property of the DESIGN rather than a coincidence of iterations."""
    X, y = fixture()
    fit = model.firth(X, y)
    Xc = np.column_stack([np.ones(len(y)), X.to_numpy(dtype=float)])
    w = np.clip(fit.p * (1.0 - fit.p), config.FIRTH_WEIGHT_FLOOR, None)
    Xw = Xc * np.sqrt(w)[:, None]
    norm = float(np.linalg.norm(np.linalg.inv(Xw.T @ Xw)))
    assert (norm > 1e4) is ill_conditioned
    assert (fit.converged_on == "score") is ill_conditioned


def test_FIRTH_WEIGHT_FLOOR_is_never_active_on_any_fixture_here():
    """The honest form of coverage for a branch no input reaches (§5.6, §12.4a).

    This is a REAL assertion that fails if the floor is raised or a fixture starts saturating — not a
    test of a branch that cannot run. It is also what keeps `_penalised_loglik`'s -inf branch
    unreachable, so the two travel together and this one assertion covers both.

    **Scoped deliberately.** The floor is NOT unreachable in general: measured, on 800 synthetic
    completely-separated designs it binds on 68 (8.5%) and the smallest p(1−p) observed is exactly 0.0.
    So this is a statement about these fixtures, not about the estimator.

    Asserted against `config.FIRTH_WEIGHT_FLOOR` rather than against 1e-10: the point is the
    relationship, not the number.
    """
    worst = min(float(np.min(fit.p * (1.0 - fit.p)))
                for fit in (model.firth(*f()) for f in
                            (well_behaved, separated, big_first_step, needs_halving)))
    assert worst > config.FIRTH_WEIGHT_FLOOR
    assert worst > 1e-3            # eight orders of headroom: not a near-miss dressed as a check


# --- 12.5  Firth stays finite under complete separation, where the MLE diverges [roadmap] --------------

def test_firth_is_finite_under_complete_separation():
    fit = model.firth(*separated())
    assert 0.5 < float(fit.beta[1]) < 2.0            # the pilots' threshold; measured 0.9706
    assert np.all(np.isfinite(fit.beta))
    assert np.all((fit.p > 0.0) & (fit.p < 1.0))
    assert fit.iterations == 9 and fit.converged_on == "likelihood"
    assert round(float(fit.beta[0]), 4) == -5.3385
    assert round(float(fit.beta[1]), 4) == 0.9706


def test_the_unpenalised_mle_DIVERGES_on_the_same_data():
    # This half is what makes the test about SEPARATION rather than about agreement: asserting only
    # that Firth is finite would pass on data that was never separated.
    import statsmodels.api as sm
    X, y = separated()
    Xc = sm.add_constant(X.to_numpy(dtype=float))
    newton = sm.Logit(y, Xc).fit(disp=0, warn_convergence=False)
    assert not newton.mle_retvals["converged"]
    assert abs(float(newton.params[1])) > 10.0                  # measured 71.5
    bfgs = sm.Logit(y, Xc).fit(method="bfgs", disp=0, warn_convergence=False)
    assert abs(float(bfgs.params[1])) > 10.0                    # measured 27.2, and climbing


# --- 12.6  a rank-deficient design raises, and `pinv` is not reachable ---------------------------------

def test_a_collinear_design_raises_F3_with_its_rank_and_its_columns():
    with pytest.raises(model.FitError) as e:
        model.firth(*collinear())
    assert "F3" in message_of(e)
    assert "rank 3 of 4" in message_of(e)
    assert "intercept, x1, x2, x3" in message_of(e)


def test_WITHOUT_the_guard_numpy_raises_LinAlgError_from_inside_the_loop(monkeypatch):
    # So F3's job is the MESSAGE, not the catch — and a pinv fallback would replace both. DoD-8's
    # first direction.
    monkeypatch.setattr(model, "_assert_fittable", lambda *a, **k: None)
    with pytest.raises(np.linalg.LinAlgError):
        model.firth(*collinear())


def test_a_pinv_fallback_would_return_coefficients_for_a_design_that_identifies_none():
    """DoD-8's second direction, and the one to watch — it is what the pilots ship.

    `pilots/analysis.py:366-369` catches LinAlgError and falls back to np.linalg.pinv. That is not a
    fallback to a different estimator but to a different ESTIMAND: the pseudo-inverse silently picks
    the minimum-norm solution among infinitely many.
    """
    X, y = collinear()
    Xc = np.column_stack([np.ones(len(y)), X.to_numpy(dtype=float)])
    p = np.full(len(y), 0.5)
    w = p * (1.0 - p)
    info = (Xc * w[:, None]).T @ Xc
    with pytest.raises(np.linalg.LinAlgError):
        np.linalg.inv(info)
    assert np.all(np.isfinite(np.linalg.pinv(info)))      # finite, and for an unidentified design


def test_pinv_appears_nowhere_in_the_module():
    assert _names_in_source(SOURCE, {"pinv"}) == []


def test_the_pinv_scan_fires():
    assert _names_in_source("inv = np.linalg.pinv(info)\n", {"pinv"}) == [(1, "pinv")]


@pytest.mark.parametrize("denied", sorted(config.POST_TIME_ZERO))
def test_D2_raises_on_every_post_time_zero_name(denied):
    """INVARIANT 4's test lives here.

    Parametrised over the denylist itself, so a new outcome is covered by existing: POST_TIME_ZERO is
    built from the outcome registry, and `design` is the one function every model's covariates pass
    through — which is what made this assertable at all, after being declared since Stage 1 and
    assertable nowhere.
    """
    df = frame()
    df[denied] = 1.0
    with pytest.raises(config.SchemaError) as e:
        model.design(df, ("age", denied))
    assert "D2" in message_of(e) and denied in message_of(e)


def test_D1_raises_naming_the_absent_covariate():
    with pytest.raises(config.SchemaError) as e:
        model.design(frame(), ("age", "core_above_median"))
    assert "D1" in message_of(e) and "core_above_median" in message_of(e)


def test_D3_raises_on_an_unlisted_factor_level_and_names_the_case():
    df = frame(centres=("HUG", "CHUV"))
    df.loc[df["case_id"] == "F-1", "center"] = "Bern"
    with pytest.raises(config.SchemaError) as e:
        model.design(df, ("age", "center"))
    assert "D3" in message_of(e) and "F-1" in message_of(e) and "HUG" in message_of(e)


def test_D3_does_not_require_a_case_id_column():
    # An earlier draft indexed df.loc[..., "case_id"] BEFORE its own guard, so `design` raised a bare
    # KeyError on any frame carrying a declared factor and no identifier — an undeclared schema
    # requirement on the one function every model's covariates pass through. The names are a courtesy,
    # and a courtesy may not impose a schema.
    df = frame(centres=("HUG", "CHUV")).drop(columns=["case_id"])
    df.loc[df.index[0], "center"] = "Bern"
    with pytest.raises(config.SchemaError) as e:
        model.design(df, ("age", "center"))
    assert "D3" in message_of(e) and "carries no case_id" in message_of(e)
    df_clean = frame(centres=("HUG", "CHUV")).drop(columns=["case_id"])
    model.design(df_clean, ("age", "center"))          # and a clean frame does not raise at all


def test_WITHOUT_D3_the_record_is_encoded_as_the_REFERENCE_and_the_fit_returns_a_number(monkeypatch):
    # DoD-7. The Categorical conversion turns an unlisted value into NaN and get_dummies then encodes
    # it as an all-zero row — which IS the encoding of the reference level. The record is modelled as
    # though it were at HUG, with no missing value anywhere and nothing to notice.
    monkeypatch.setattr(model, "_assert_design_inputs", lambda *a, **k: None)
    df = frame(centres=("HUG", "CHUV", "Lugano", "CHUV"))
    df.loc[df["case_id"] == "F-0", "center"] = "Bern"
    X, _ = model.design(df, ("center",))
    assert not X.isna().any().any()                    # no missing value survives the encoding
    assert (X.iloc[0] == 0.0).all()                    # byte-identical to a real HUG record
    fit = model.firth(X, np.array([0.0, 1.0, 1.0, 0.0]))
    assert np.isfinite(fit.p[0]) and 0.0 < fit.p[0] < 1.0


@pytest.mark.parametrize("factor", list(config.CATEGORICAL))
def test_D4_raises_on_a_missing_factor_value(factor):
    df = frame(centres=("HUG", "CHUV", "Lugano"))
    df.loc[df["case_id"] == "F-1", factor] = pd.NA
    with pytest.raises(config.SchemaError) as e:
        model.design(df, ("age", factor))
    assert "D4" in message_of(e) and factor in message_of(e)


def test_WITHOUT_D4_the_missing_row_is_BYTE_IDENTICAL_to_a_real_reference_record(monkeypatch):
    """The sharper of the two companions, and the reason D4 exists at all.

    F4 cannot see it — the missingness was consumed by the encoding, so no non-finite value is left.
    D3 cannot see it, by its own `& notna()`. F5 cannot see it — the fitted probability is finite and
    interior; it is simply the probability of a patient the model believes is at the reference centre.
    """
    monkeypatch.setattr(model, "_assert_design_inputs", lambda *a, **k: None)
    df = frame(centres=("HUG", "CHUV", "Lugano", "HUG"))
    df.loc[df["case_id"] == "F-1", "center"] = pd.NA
    X, _ = model.design(df, ("center",))
    missing_row, real_reference_row = X.iloc[1], X.iloc[0]
    assert list(missing_row) == list(real_reference_row)          # element for element
    assert not X.isna().any().any()                              # F4 has nothing to test
    fit = model.firth(X, np.array([0.0, 1.0, 1.0, 0.0]))
    assert 0.0 < fit.p[1] < 1.0                                  # F5 has nothing to test


def test_D3_and_D4_are_DISTINGUISHABLE_in_one_message():
    # They share a loop and they are different failures: a level outside the declared set, and no level
    # at all. Without D3's `& notna()` clause the two would report as one, with one message describing
    # the wrong mechanism.
    df = frame(centres=("HUG", "CHUV", "Lugano"))
    df.loc[df["case_id"] == "F-1", "center"] = "Bern"
    df.loc[df["case_id"] == "F-2", "center"] = pd.NA
    with pytest.raises(config.SchemaError) as e:
        model.design(df, ("age", "center"))
    text = message_of(e)
    assert "D3" in text and "D4" in text
    assert "F-1" in text                                          # named by D3
    assert "2 design assertion(s) failed" in text


def test_F6_raises_on_an_ordinal_response():
    # Parametrised over the response a caller is most likely to hand this function by mistake: a 0-6
    # functional scale. It is the [§8] primary outcome, and [§8]'s proportional-odds model is a
    # different function that does not exist yet.
    X = pd.DataFrame({"x": np.arange(1.0, 8.0)})
    ordinal = np.array([0, 1, 2, 3, 4, 5, 6], dtype=float)
    with pytest.raises(model.FitError) as e:
        model.firth(X, ordinal)
    assert "F6" in message_of(e) and "outside {0, 1}" in message_of(e)
    assert "proportional-odds" in message_of(e)


def test_WITHOUT_F6_an_ordinal_response_FITS_CONVERGES_and_returns_coefficients(monkeypatch):
    """DoD-10a's ordinal half — the only SILENT failure this stage adds a guard for.

    Nothing anywhere in `model.py` establishes that y is a dichotomy, so this fits and returns
    coefficients for a logistic model of a seven-level ordinal, with no missing value and no warning.
    It is the one Stage 9 will meet first.
    """
    monkeypatch.setattr(model, "_assert_fittable", lambda *a, **k: None)
    X = pd.DataFrame({"x": np.arange(1.0, 8.0)})
    fit = model.firth(X, np.array([0, 1, 2, 3, 4, 5, 6], dtype=float))
    assert np.all(np.isfinite(fit.beta))
    assert fit.converged_on in ("likelihood", "score")
    assert fit.iterations >= 1


# --- 12.7  no fallback estimator exists on any path [invariant 5] ---------------------------------------

SHIPPED = ("model.py", "propensity.py")
FORBIDDEN_IMPORTS = {"sklearn", "scikit-learn", "statsmodels", "scipy"}


def _imported_roots(source: str) -> set[str]:
    roots: set[str] = set()
    for node in ast.walk(ast.parse(source)):
        if isinstance(node, ast.Import):
            roots.update(alias.name.split(".")[0] for alias in node.names)
        if isinstance(node, ast.ImportFrom) and node.module:
            roots.add(node.module.split(".")[0])
    return roots


@pytest.mark.parametrize("name", SHIPPED)
def test_no_shipped_module_imports_a_second_estimator(name):
    assert _imported_roots((MODULE_DIR / name).read_text(encoding="utf-8")) & FORBIDDEN_IMPORTS == set()


def test_the_import_scan_fires():
    assert _imported_roots("import statsmodels.api as sm\nfrom sklearn import linear_model\n") & \
        FORBIDDEN_IMPORTS == {"statsmodels", "sklearn"}


def _returning_except_clauses(source: str) -> list[int]:
    return [handler.lineno
            for node in ast.walk(ast.parse(source)) if isinstance(node, ast.Try)
            for handler in node.handlers
            if any(isinstance(n, ast.Return) and n.value is not None
                   for n in ast.walk(handler))]


@pytest.mark.parametrize("name", SHIPPED)
def test_no_except_clause_returns_a_value(name):
    # `_penalised_loglik`'s -inf branch is NOT an exception handler; it is a value the step-halving
    # loop is designed to reject. An except that returns is the pilots' 16.6% fallback in miniature.
    assert _returning_except_clauses((MODULE_DIR / name).read_text(encoding="utf-8")) == []


def test_the_except_returns_scan_fires():
    pasted = "def f():\n    try:\n        return g()\n    except Exception:\n        return 0.0\n"
    assert _returning_except_clauses(pasted) == [4]


@pytest.mark.parametrize("name", SHIPPED)
def test_only_FitError_and_SchemaError_are_raised(name):
    """So Stage 10 can catch the first to drop a replicate and not the second (§9).

    A SchemaError from F1/F2 or D1-D4 is a bug in the resampler, not a sparse replicate.
    """
    source = (MODULE_DIR / name).read_text(encoding="utf-8")
    raised = set()
    for node in ast.walk(ast.parse(source)):
        if isinstance(node, ast.Raise) and isinstance(node.exc, ast.Call):
            target = node.exc.func
            raised.add(target.attr if isinstance(target, ast.Attribute) else
                       target.id if isinstance(target, ast.Name) else None)
    assert raised <= {"FitError", "SchemaError"}


def test_firthlogist_appears_only_in_this_file_and_in_pyproject():
    offenders = [str(path.relative_to(MODULE_DIR))
                 for path in _python_files()
                 if "firthlogist" in path.read_text(encoding="utf-8")
                 and path.name != "test_model.py"]
    assert offenders == []


def _python_files() -> list[Path]:
    found: list[Path] = []
    for dirpath, dirnames, filenames in os.walk(MODULE_DIR):
        dirnames[:] = sorted(d for d in dirnames if d not in {".venv", "__pycache__"})
        found.extend(Path(dirpath) / n for n in sorted(filenames) if n.endswith(".py"))
    return found


def test_no_shipped_module_imports_subprocess():
    # Several test modules do, and legitimately: every stage since Stage 2 has a two-seed
    # byte-identity driver that spawns an interpreter. What must stay true is that nothing on an
    # ESTIMATION path can shell out — which is what keeps §16b.1's oracle-not-dependency line drawn.
    offenders = [name for name in SHIPPED
                 if "subprocess" in _imported_roots((MODULE_DIR / name).read_text(encoding="utf-8"))]
    assert offenders == []


def test_IPTW_TRIM_appears_nowhere_in_the_repository():
    # [§15] lists unrestricted and trimmed IPTW among the approaches deliberately not used, overlap
    # weights are bounded by construction, and `config.py` declares no such constant.
    # `pilots/analysis.py:127-133` has a `scheme="iptw"` branch reading it; §16 declines the branch.
    #
    # This file is excluded from its own scan, for the reason `test_firthlogist_appears_only_...` is:
    # a test that names the forbidden string in order to search for it would otherwise be its own
    # only offender, and the repair would be to gut the scan.
    assert not hasattr(config, "IPTW_TRIM")
    offenders = [str(p.relative_to(MODULE_DIR)) for p in _python_files()
                 if p.name != "test_model.py" and "IPTW_TRIM" in p.read_text(encoding="utf-8")]
    assert offenders == []


# --- 12.12  the module boundary holds --------------------------------------------------------------------

def _config_attributes(source: str) -> set[str]:
    """Every `C.<NAME>` / `config.<NAME>` attribute read in `source`."""
    return {node.attr for node in ast.walk(ast.parse(source))
            if isinstance(node, ast.Attribute) and isinstance(node.value, ast.Name)
            and node.value.id in {"C", "config"}}


def test_model_names_neither_the_exposure_nor_the_covariate_list():
    # §0.1: Stage 9 must be able to import this module to fit an OUTCOME model without importing an
    # exposure. Companion below.
    assert {"TREATMENT", "PS_COVARIATES"} & _config_attributes(SOURCE) == set()


def test_no_config_attribute_read_by_model_begins_with_PS_():
    """The seven tolerances were named PS_* in an earlier draft — for the propensity score, inside the
    module defined by knowing nothing about it — and they are FIRTH_* now.

    The attribute scan above would have passed the old names, because they are neither TREATMENT nor
    PS_COVARIATES. This is what makes the boundary a rule rather than a habit, and what stops Stage 9
    reading a propensity-flavoured constant for an outcome fit. `propensity.py` is exempt:
    PS_COVARIATES is exactly what it is for.
    """
    assert [name for name in _config_attributes(SOURCE) if name.startswith("PS_")] == []


def test_the_attribute_scan_fires():
    assert _config_attributes("a = C.TREATMENT\nb = config.PS_COVARIATES\n") == {
        "TREATMENT", "PS_COVARIATES"}


def _raise_strings(source: str) -> list[str]:
    """Every string constant appearing inside a `raise` statement's arguments."""
    return [n.value
            for node in ast.walk(ast.parse(source)) if isinstance(node, ast.Raise)
            for n in ast.walk(node)
            if isinstance(n, ast.Constant) and isinstance(n.value, str)]


FORBIDDEN_IN_RAISES = ("ivt", "TREATMENT", "propensity", *config.OUTCOMES)


@pytest.mark.parametrize("word", FORBIDDEN_IN_RAISES)
def test_no_runtime_message_in_model_names_the_exposure_or_an_outcome(word):
    """Scoped to `raise` arguments rather than to all strings, and BOTH halves are the point.

    It is not belt-and-braces over the attribute scan: an attribute scan passes a module whose error
    *text* is full of exposure names, and F6's first draft was exactly that. A module that says that in
    a raise is one the next reader treats as knowing about the exposure.

    And it is scoped to raises because THE MODULE DOCSTRING MUST NAME `propensity.py` — the sentence
    "Stages 8, 9 and 12 come into this module directly and never through propensity.py" is the one
    thing that stops someone moving `design` there for tidiness. A whole-file text scan would forbid
    the comment this stage most wants to keep. So: the docstring may name the module it must not depend
    on; a raise may not.

    Parametrised over `C.OUTCOMES`, so a ninth outcome is covered by the existing test.
    """
    assert [s for s in _raise_strings(SOURCE) if word in s] == []


def test_the_raise_string_scan_fires_and_is_SCOPED_to_raises():
    pasted = '"""A docstring naming propensity.py."""\nraise ValueError("the ivt column is missing")\n'
    assert [s for s in _raise_strings(pasted) if "ivt" in s] == ["the ivt column is missing"]
    assert [s for s in _raise_strings(pasted) if "propensity" in s] == []   # the docstring is exempt


def test_model_takes_no_audit_parameter_and_imports_neither_data_nor_propensity():
    assert _imported_roots(SOURCE) & {"data", "propensity"} == set()
    assert "Audit" not in SOURCE.split('"""')[2]        # not in any signature after the docstring


@pytest.mark.parametrize("name", SHIPPED)
def test_neither_module_contains_a_bare_assert(name):
    # Python strips `assert` under -O, so under a flag nobody remembers setting, a module whose entire
    # purpose is to fail loudly would succeed silently.
    source = (MODULE_DIR / name).read_text(encoding="utf-8")
    assert [n.lineno for n in ast.walk(ast.parse(source)) if isinstance(n, ast.Assert)] == []


def test_the_bare_assert_scan_fires():
    assert [n.lineno for n in ast.walk(ast.parse("def f(x):\n    assert x\n"))
            if isinstance(n, ast.Assert)] == [2]


@pytest.mark.parametrize("name", SHIPPED)
def test_neither_module_is_exempt_from_the_raw_name_scan(name):
    assert name not in config.EXEMPT_FROM_RAW_NAME_SCAN


# --- 16b  the oracles ------------------------------------------------------------------------------------

def test_scipy_reaches_the_same_penalised_optimum_as_the_newton_loop():
    """§16b's third oracle — an INDEPENDENT optimiser on the SAME objective.

    This validates the loop independently of any Firth implementation, and it is the check no external
    package can provide: an external package would be validating its own algorithm rather than our
    objective.

    **Scoped deliberately.** An earlier wording said the loop's optimum *is* the objective's optimum
    without qualification, and that is false on separated designs where the penalised likelihood is
    multimodal (§5.2c). It is asserted on the fixtures where it holds. The assertion is on the
    OBJECTIVE VALUE and not on the coefficients, because the surface is flat (§5.3) — on the cohort
    design the coefficients agree only to 2.7e-05 while the objective agrees to ten decimals.
    """
    from scipy.optimize import minimize
    X, y = well_behaved(n=400)
    fit = model.firth(X, y)
    Xc = np.column_stack([np.ones(len(y)), X.to_numpy(dtype=float)])
    ours = model._penalised_loglik(Xc, y, fit.beta)
    theirs = minimize(lambda b: -model._penalised_loglik(Xc, y, b),
                      np.zeros(Xc.shape[1]), method="BFGS")
    assert ours >= -float(theirs.fun) - 1e-6


def _firthlogist_beta(X, y):
    """Fit `firthlogist` with OUR tolerances passed explicitly, and return one beta vector.

    Passed explicitly rather than relying on shared defaults, for §16b.2's reason: what makes a
    comparison like-for-like is the explicit pass, not a default two packages happen to share.
    """
    firthlogist = pytest.importorskip("firthlogist")
    fitted = firthlogist.FirthLogisticRegression(
        max_iter=config.FIRTH_MAX_ITER, max_halfstep=config.FIRTH_MAX_HALVINGS,
        tol=config.FIRTH_TOL)
    fitted.fit(X.to_numpy(dtype=float), y)
    return np.concatenate([[float(fitted.intercept_)], np.asarray(fitted.coef_, dtype=float)])


def _agrees_with_firthlogist(X, y):
    """§16b's first oracle: assert on the OBJECTIVE, and on the coefficients only to 1e-4.

    **The spec asks for 1e-12 on the coefficients and that is not achievable, for a reason worth
    more than the agreement would have been** — which is exactly what §16b.2 predicts an oracle is
    for. Read from `firthlogist/firthlogist.py:307`, its convergence test is

        if iter > 1 and np.linalg.norm(coef_new - coef) < tol:

    a **coefficient-step test and nothing else** — the `xconv` criterion §5.3 deliberately
    excludes — with no likelihood-change test and no score test. It also bounds the step by an
    ABSOLUTE infinity norm (`max|step| / max_stepsize`, the form §5.2a replaced) and solves with
    `lstsq` rather than `inv` (the pseudo-inverse §5.4 declines).

    So the two implementations stop by different rules on a surface §5.3 describes as genuinely flat,
    and they stop in different places: measured 3.66e-06 apart on the cohort design and 7.75e-05 on
    `separated()`, while agreeing on the penalised log-likelihood to 2e-10. **`pilots/analysis.py`'s
    kernel — the one §18 recorded as agreeing to 2.665e-15 — is reproduced here to `max|Δβ| = 0` and
    disagrees with `firthlogist` by the same 3.66e-06**, so that row is not reproducible from either
    kernel rather than being a property of this implementation.

    The assertion is therefore on the objective value, which is the form §16b itself calls the honest
    one for the `scipy` oracle and for the same reason. It still fails on a real kernel defect: a
    wrong penalty, score or hat diagonal moves the optimum, not merely where on the ridge one stops.
    """
    ours = model.firth(X, y)
    theirs = _firthlogist_beta(X, y)
    Xc = np.column_stack([np.ones(len(y)), X.to_numpy(dtype=float)])
    assert abs(model._penalised_loglik(Xc, y, ours.beta)
               - model._penalised_loglik(Xc, y, theirs)) < 1e-9
    assert float(np.max(np.abs(ours.beta - theirs))) < 1e-4


def test_firthlogist_agrees_with_our_kernel_on_the_separation_benchmark():
    """§16b's first oracle, dev-only and behind importorskip so a plain checkout is green without it.

    `firthlogist` is a PORT of R's reference implementation, so this agreement is transitive rather
    than independent —
    which is why DoD-16's direct comparison against the R original is a separate gate.
    """
    _agrees_with_firthlogist(*separated())


def test_firthlogists_convergence_rule_is_the_one_the_spec_deliberately_EXCLUDES():
    """Pinned as a test, because it is the finding rather than an incidental detail (DoD-15a).

    §5.3 carries a claim from `pilots/analysis.py`'s comment — "converge on the penalised likelihood
    or on a flat modified score, **as the R reference implementation does**" — which the spec marks
    unverified and
    withdraws pending a read of the source. `firthlogist` is a port of that implementation, and its
    rule is a
    coefficient-step test **only**. That is evidence *against* the attribution, not for it.

    §5.3's ARGUMENT is untouched either way and does not depend on the attribution: excluding the step
    is right because the penalised surface is genuinely flat under near-separation, so a step-norm
    criterion would report a finite correct fit as a failure — and [§10] drops failed replicates, so
    the ones dropped would be exactly the sparse ones. What this kills is the attribution, not the
    rule. If a future reader restores that attribution to §5.3 or to `config.py`, this test is what
    should stop them.
    """
    firthlogist = pytest.importorskip("firthlogist")
    source = Path(firthlogist.__file__).with_name("firthlogist.py").read_text(encoding="utf-8")
    assert "np.linalg.norm(coef_new - coef) < tol" in source      # xconv, and only xconv
    assert "loglike_new - loglike" not in source                  # no lconv
    assert "U_star" in source and "np.max(np.abs(U_star))" not in source   # score built, never tested


# --- 12.14 (model's half)  structural facts from the workbook [data-gated] --------------------------------

@pytest.fixture(scope="module")
def workbook():
    """This module's OWN fixture, as Stage 5's test file has its own. Never shared across modules."""
    import cohort
    import data
    import derive
    import eligibility
    df, audit = data.load(data.WORKBOOK)
    df = derive.derive(df, audit)
    df = eligibility.classify(df, audit)
    return cohort.build(df, audit), audit


@DATA_GATED
def test_the_cohort_design_is_eleven_columns_and_center_USZ_is_dropped(workbook):
    df, _ = workbook
    X, dropped = model.design(df, config.PS_COVARIATES)
    assert X.shape[1] == 11
    assert dropped == ("center_USZ",)
    assert X.shape[1] + 1 == 12                       # parameters, including the intercept


@DATA_GATED
def test_the_cohort_fit_converges_in_seven_iterations_on_the_likelihood(workbook):
    df, _ = workbook
    mask = model.complete_cases(df, config.PS_COVARIATES)
    X, _ = model.design(df.loc[mask], config.PS_COVARIATES)
    fit = model.firth(X, df.loc[mask, config.TREATMENT].to_numpy(dtype=float))
    assert fit.iterations == 7
    assert fit.converged_on == "likelihood"
    assert (fit.rescales, fit.halvings) == (0, 0)
    assert float(fit.p.min()) > 0.02 and float(fit.p.max()) < 0.94
    assert float(np.min(fit.p * (1 - fit.p))) > config.FIRTH_WEIGHT_FLOOR


@DATA_GATED
def test_firthlogist_agrees_with_our_kernel_on_the_cohort_design(workbook):
    df, _ = workbook
    mask = model.complete_cases(df, config.PS_COVARIATES)
    X, _ = model.design(df.loc[mask], config.PS_COVARIATES)
    _agrees_with_firthlogist(X, df.loc[mask, config.TREATMENT].to_numpy(dtype=float))


# =====================================================================================================
#  STAGE 8 — the weighted proportional-odds fit. `specs/stage8_primary_outcome_estimator.md` §14.
#
#  `polr`'s sections live here because that is where the function is: §14.3, §14.4, §14.7's fit half
#  and §14.8. The sections that read a `Primary` or an `Audit` are `test_outcome.py`'s.
# =====================================================================================================

# --- 14.0.4  the synthetic constructions, written out ------------------------------------------------
#
# Six deterministic generators, no assertions, no fixtures, no frames. They are written out rather than
# described because a document that pins a number produced by a construction it does not specify has a
# hole exactly the shape of the number: "a 200-record construction with the treated arm's mass at
# mRS 0-3" does not yield exp(beta) = 84.58, and an implementer cannot reach it from that sentence.
#
# EVERY ONE SEEDS FROM `config.SEED` AND NONE HOLDS A LITERAL SEED, for FIRTH_TOL's reason one level
# down: a seed is what makes a measured number reproducible, so a second seed in a test module is a
# second answer to a question config.py already answers. Stage 10 reads the same constant.
#
# `test_outcome.py` imports them by name, exactly as it imports `cohort_frame` from `test_cohort`.


def hand_ordinal():
    """12 records, 6 against 6, all seven declared mRS levels. (y, a, w) — §14.3, §14.6.

    Hand-computable and hand-checkable: the weights are small integers summing to 21 in each arm, so
    every weighted cumulative probability in §14.6 is a sum of at most six of them over a denominator
    of 21.

    **THE WEIGHTS ARE INTEGERS, AND THAT IS A REQUIREMENT AND NOT A CONVENIENCE.** §14.3's and §18b's
    first oracle run the integer-weight replication check on this frame, and "the unweighted fit on
    the frame with each row repeated `w` times" is undefined for a fractional weight — `np.repeat`
    raises `TypeError: Cannot cast array data from dtype('float64') to dtype('int64')`.

    **And scaling fractional weights up to integers does NOT work**, which is measured rather than
    assumed: the polr MLE is invariant to a common weight scale in exact arithmetic, but `POLR_TOL` is
    ABSOLUTE on a log-likelihood that scales with the weight sum, so multiplying every weight by 20
    tightens the stopping criterion twentyfold and the two fits stop at different points —
    |Δbeta| = 9.458e-10 measured, which FAILS §14.3's 1e-10 tolerance. Integers chosen from the start
    agree at 1.11e-16. That is the first measured consequence of §5.2's absolute tolerances.
    """
    y = np.array([0., 1., 2., 3., 4., 5., 1., 2., 3., 4., 5., 6.])
    a = np.array([1.] * 6 + [0.] * 6)
    w = np.array([1., 2., 3., 4., 5., 6., 6., 5., 4., 3., 2., 1.])
    return y, a, w


def replication_frame(n: int = 120, levels: int = 5):
    """(X, y, w) with INTEGER weights in {1, 2, 3, 4}, for the replication oracle. §14.3, §18b.

    Integer by construction and not by rounding: the oracle is that a weighted fit equals the
    unweighted fit on the row-replicated frame, and a non-integer weight makes "replicated" undefined.
    Two covariates so the oracle covers a beta of length > 1, which the [§8] path never has and
    Stage 12 always does.

    **THE LOGISTIC NOISE TERM IS LOAD-BEARING.** This is the proportional-odds data-generating
    process: a latent variable `x'b + Logistic(0, 1)`, cut at its own quantiles. Without the noise,
    `y` is a MONOTONE function of `x'b` — every category is an interval of the linear predictor — and
    the frame is **perfectly separated in the ordinal sense**, so the likelihood has no interior
    maximum and `beta` runs to the trust region's limit. Measured on a noiseless draft:
    `beta = [-3281.7, -1652.8]` at n=120 and `[-48828, -24382]` at n=500, both converging on the
    score criterion with clean counters. **Every oracle built on this frame is void on a separated
    one**: two optimisers of a likelihood whose maximum is at infinity need not agree at all. With the
    noise, `max|beta|` is 0.93 here and 0.76 on `reference_frame()` — comfortably interior.
    """
    rng = np.random.default_rng(config.SEED)
    X = pd.DataFrame({"x1": rng.normal(size=n), "x2": rng.normal(size=n)})
    latent = (0.8 * X["x1"].to_numpy() + 0.4 * X["x2"].to_numpy()
              + rng.logistic(size=n))                          # <- the noise
    cuts = np.quantile(latent, np.linspace(0.0, 1.0, levels + 1)[1:-1])
    y = np.searchsorted(cuts, latent).astype(float)
    w = rng.integers(1, 5, size=n).astype(float)
    return X, y, w


def reference_frame(n: int = 500, levels: int = 5):
    """(X, y), two covariates, unweighted — the frame the statsmodels sign oracle runs on. §7.2, §14.5.

    Unweighted deliberately: `OrderedModel` has no weight support of any kind, so the only comparison
    available is the unweighted one, and this generator is what makes `ours.beta + sm.beta` a
    reproducible quantity rather than a remembered one.

    It delegates to `replication_frame` rather than repeating the construction, so it inherits the
    noise term that docstring is about — and it inherited the BUG for the same reason, at n = 500
    where separation is sharper.
    """
    X, y, _ = replication_frame(n=n, levels=levels)
    return X, y


def orientation_frame(per_arm: int = 100):
    """(X, y, w) where treatment shifts mRS DOWNWARD by construction. §7.3, §14.5.

    The construction whose truth is known, which is what the roadmap's orientation criterion needs and
    what prose cannot supply: the treated arm's mass sits at mRS 0-3 and the control's at 2-6, so the
    direction of the true effect is a property of these two literal tuples and not of a fitted number.
    Unit weights — the orientation of `beta` is a question about the parametrisation and not about the
    weighting, and §14.3's oracles own the weighting.
    """
    rng = np.random.default_rng(config.SEED)
    treated = rng.choice(np.array([0., 1., 2., 3.]), size=per_arm, p=[.40, .30, .20, .10])
    control = rng.choice(np.array([2., 3., 4., 5., 6.]), size=per_arm, p=[.20, .25, .25, .20, .10])
    y = np.concatenate([treated, control])
    a = np.concatenate([np.ones(per_arm), np.zeros(per_arm)])
    return pd.DataFrame({config.TREATMENT: a}), y, np.ones(len(y))


def separated_frame(per_arm: int = 20, crossovers: int = 0):
    """PERFECT separation at crossovers = 0, NEAR separation above it. §6.1, §14.7.

    Every treated record at mRS 0 and every control at mRS 5, with `crossovers` treated records moved
    to the control's level. This is the construction §6.1 measured — it holds no seed at all, which is
    why §6.1's numbers are exact rather than distributional: at crossovers = 0 and per_arm = 20 it
    returns beta 36.4058 in 17 iterations on the score criterion, every counter at zero, and it does
    so on every machine.
    """
    y = np.concatenate([np.zeros(per_arm), np.full(per_arm, 5.0)])
    y[:crossovers] = 5.0
    a = np.concatenate([np.ones(per_arm), np.zeros(per_arm)])
    return pd.DataFrame({config.TREATMENT: a}), y, np.ones(len(y))


BAND_EFFECTS: Final[tuple[float, ...]] = (0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 6.0)


def band_samples(effects, per_effect: int, per_arm: int = 46, levels: int = 7):
    """The FITS from `len(effects) * per_effect` draws of a cohort-shaped frame. §6.3, §14.7.

    It returns the fits and not a `|beta|` vector, which is a correction: §14.7 asserts things about
    `cond(-H)` and about the safeguard counters, and both are properties of the fit. An earlier draft
    returned magnitudes only and left those assertions with no subject. The band assertion takes
    `max(abs(f.beta))` from each; nothing is lost by returning more.

    A treatment-only weighted proportional-odds sampler over a range of TRUE effect sizes. The range
    is the point: §6.3's band was first measured at a single true effect of 0.5, which is an
    assumption about the answer, so the generator takes the effects as an argument and both callers
    pass a range.

    TWO callers with different budgets, which is §14.7's split. The in-suite probe passes
    per_effect = 20 for 240 fits and under a second; §6.3's calibration passed per_effect = 400 for
    4800 fits and about 14 s, which is a calibration and not a test.
    """
    rng = np.random.default_rng(config.SEED)
    out: list[model.PolrFit] = []
    for effect in effects:
        for _ in range(per_effect):
            a = np.concatenate([np.ones(per_arm), np.zeros(per_arm)])
            latent = effect * a + rng.normal(size=2 * per_arm)
            cuts = np.quantile(latent, np.linspace(0.0, 1.0, levels + 1)[1:-1])
            y = np.searchsorted(cuts, latent).astype(float)
            w = rng.uniform(0.05, 0.95, size=2 * per_arm)
            try:
                out.append(model.polr(pd.DataFrame({config.TREATMENT: a}), y, w))
            except model.FitError:
                continue                       # O5 on a degenerate draw; NOT a band member
    return out


# --- 14.0.4a  the six that O1-O6, §14.4 and §14.6 need -------------------------------------------------
#
# Hand-built, seed-free and small, so their numbers are VALUES rather than tolerances. §14.0.4's
# generators seed from config.SEED and their sampling differs in detail from the harnesses that first
# produced the spec's §20 numbers, so every assertion resting on one of those is a tolerance or a
# property. These six and `separated_frame()` are the exceptions.


def nan_weight_frame():
    """A four-category weighted frame, for O3's ordering. §14.8.

    Twenty records over four categories with strictly positive weights. The test sets `w[0]` to nan
    and asserts O3; the companion removes O3 and asserts the category DELETION, which is the mechanism
    no other check in this stage names.
    """
    y = np.array([0., 1., 2., 3.] * 5)
    a = np.array([1., 0.] * 10)
    return pd.DataFrame({config.TREATMENT: a}), y, np.linspace(0.2, 0.9, 20)


def nan_response_frame(n: int = 80):
    """Eighty records over four categories, for O2 — and for O2's interaction with O4. §14.8."""
    rng = np.random.default_rng(config.SEED)
    a = np.array([1., 0.] * (n // 2))
    latent = 0.9 * a + rng.logistic(size=n)
    cuts = np.quantile(latent, [0.25, 0.50, 0.75])
    return (pd.DataFrame({config.TREATMENT: a}),
            np.searchsorted(cuts, latent).astype(float), np.ones(n))


def noninteger_response_frame():
    """THREE distinct half-integer values over sixty records, for O4. §14.8.

    Sized so that the with-O4-removed companion reaches a FIT. Six distinct values over 24 records
    makes the Hessian singular and `np.linalg.solve` raises `LinAlgError` — which is not `FitError`,
    so the companion would assert the wrong thing (measured).
    """
    y = np.array([0., 0.5, 1.] * 20)
    a = np.array([1., 0.] * 30)
    return pd.DataFrame({config.TREATMENT: a}), y, np.ones(60)


def emptied_category_frame(empty: int = 2):
    """(full, reduced) — a five-category frame and the same frame with one category removed. §14.4.

    Returns BOTH, because §14.4's assertion is that the collapse gives the same fit as the removal.
    The zero-weight variant is the same frame with `w` zeroed on that category rather than the rows
    dropped, which is the one input separating §5.3's rule from the observed-set rule.
    """
    y = np.array([v for v in (0., 1., 2., 3., 4.) for _ in range(6)])
    a = np.array([1., 0.] * 15)
    w = np.linspace(0.3, 0.8, 30)
    keep = y != float(empty)
    return ((pd.DataFrame({config.TREATMENT: a}), y, w),
            (pd.DataFrame({config.TREATMENT: a[keep]}), y[keep], w[keep]))


def rare_category_frame():
    """Sixty-one records, ONE of them in the top category. §14.7.

    The frame on which `_assert_reportable` is shown to bound `beta` and not `alpha`: `max|alpha|` is
    4.1599 and `|beta|` 0.126998, so with the bound patched to 3.0 an alpha-bound would reject it and
    the real guard passes. No frame with `|alpha| >= 14` is constructible — alpha for a rare category
    grows like log n and 14 needs order 1e7 records.
    """
    y = np.array([0.] * 30 + [1.] * 30 + [2.])
    a = np.array(([1., 0.] * 15) * 2 + [1.])
    return pd.DataFrame({config.TREATMENT: a}), y, np.ones(61)


def deaths_both_arms_frame():
    """Fourteen records, all seven levels in both arms, a death in each. §14.6, §8.3.

    Returns the arm vector as well, because §14.6 calls `cumulative_rd` on it directly. It is the
    frame on which RD_6 is shown structurally zero while RD_5 is not.
    """
    y = np.array([0., 1., 2., 3., 4., 5., 6.] * 2)
    a = np.array([1.] * 7 + [0.] * 7)
    w = np.array([3., 3., 3., 2., 2., 1., 1., 1., 2., 2., 3., 3., 3., 4.])
    return pd.DataFrame({config.TREATMENT: a}), y, a, w


# --- 14.4  the two fixtures that drive §5.2's safeguards non-zero --------------------------------------
#
# Measured: with `rescales` and `halvings` asserted only as integers on frames where they are zero,
# BOTH of §5.2's safeguards can be deleted from `polr` and every frame in the suite returns
# bit-identically. Stage 6 §12.4a does not have that gap — `big_first_step()` and `needs_halving()` are
# the fixtures it uses — and this stage needs its own two.
#
# They are one line each, because a proportional-odds estimate is equivariant under rescaling a
# covariate and the trust radius is RELATIVE, so scaling the exposure is enough to drive both counters:
#
#     exposure scaled by      1.0        0.1        0.01
#     rescales                  0          1           2
#     halvings                  0          0           1
#     first_step_norm      4.8078    31.7083    315.0171
#
# The same pair is the scale-invariance test: `beta * scale` is equal across the three to 1e-8 — not
# tighter, because POLR_TOL is absolute and the three fits stop at slightly different points, which is
# itself the measurement.


def scaled_exposure(scale: float):
    """`orientation_frame()` with the exposure column multiplied by `scale`. §14.4."""
    X, y, w = orientation_frame()
    return X * scale, y, w


def polr_big_first_step():
    """§14.4 — the TRUST REGION, in isolation: 1 rescale and 0 halvings."""
    return scaled_exposure(0.1)


def polr_needs_halving():
    """§14.4 — both safeguards: 2 rescales and 1 halving."""
    return scaled_exposure(0.01)


def cohort_shaped_frame(effect: float, per_arm: int = 46, levels: int = 7):
    """One (X, y, w) draw of `band_samples`'s frame — the FRAME, which `band_samples` discards. §14.7.

    `cond(-H)` is a property of the (frame, fit) PAIR: the Hessian is a function of the design, the
    response indices and the weights, and a `PolrFit` carries none of the three. `band_samples`
    returns fits so the band assertion and the counter assertion have their subject; this returns the
    first draw at one true effect so the conditioning assertion has its own. Same sampler, same seed,
    same shape — 46 against 46 over seven declared levels, which is the workbook cohort's.

    At `effect = 0.5` it is a healthy seven-category frame (`cond(-H)` 58.4); at `effect = 6.0` it is
    near-separated (`cond(-H)` 2.8e8, `|beta|` 21.3). The gap is six orders of magnitude, which is
    §6.3's measurement — and the reason that detector is still rejected is not that it fails to
    discriminate but that the matrix's DIMENSION is a property of the sample.
    """
    rng = np.random.default_rng(config.SEED)
    a = np.concatenate([np.ones(per_arm), np.zeros(per_arm)])
    latent = effect * a + rng.normal(size=2 * per_arm)
    cuts = np.quantile(latent, np.linspace(0.0, 1.0, levels + 1)[1:-1])
    y = np.searchsorted(cuts, latent).astype(float)
    w = rng.uniform(0.05, 0.95, size=2 * per_arm)
    return pd.DataFrame({config.TREATMENT: a}), y, w


def ord_pieces_of(X, y, w, fit):
    """(Xn, y_idx, w, K) for a fitted frame — `polr`'s collapse, redone so a test can call a private.

    Written once here rather than inlined at five call sites, and it is deliberately NOT a second
    copy of the collapse rule: it reads `fit.categories`, which is what `polr` already decided. A
    helper that recomputed the rule would be exactly the duplication §5.3a took out of `model.py`.
    """
    categories = np.asarray(fit.categories, dtype=float)
    keep = np.isin(y, categories)
    Xn, y_kept, w_kept = X.to_numpy(dtype=float)[keep], y[keep], w[keep]
    return Xn, np.searchsorted(categories, y_kept), w_kept, len(fit.categories) - 1


def category_probabilities(X, y, w, fit):
    """P(Y = y_i) under `fit`, per observation — §14.7's fitted-probability detector."""
    Xn, y_idx, w_kept, K = ord_pieces_of(X, y, w, fit)
    upper, lower = y_idx, y_idx - 1
    gu, gl = model._ord_pieces(Xn @ fit.beta, fit.alpha, upper, lower,
                               y_idx <= K - 1, y_idx >= 1)[:2]
    return gu - gl


def cond_of(X, y, w, fit):
    """`cond(-H)` at the returned optimum. §14.7."""
    Xn, y_idx, w_kept, K = ord_pieces_of(X, y, w, fit)
    _, H = model._ord_score_hess(Xn, y_idx, w_kept, fit.alpha, fit.beta, K)
    return float(np.linalg.cond(-H))


def row_replicated(X, y, w):
    """The frame with each row repeated `w` times. Integer weights only — §14.3, §18b's oracle 1."""
    counts = np.asarray(w, dtype=int)
    assert np.array_equal(counts.astype(float), np.asarray(w, dtype=float)), (
        "row replication is undefined for a fractional weight; the generator must supply integers")
    return (pd.DataFrame(np.repeat(X.to_numpy(dtype=float), counts, axis=0), columns=X.columns),
            np.repeat(np.asarray(y, dtype=float), counts))


# The three frames every "on every frame in the suite" assertion below ranges over: a weighted
# treatment-only fit, a two-covariate weighted fit, and an unweighted one. Named once so a fourth is
# one line and so no assertion says "every frame" while looking at one.
def _hand_ordinal_triple():
    y, a, w = hand_ordinal()
    return pd.DataFrame({config.TREATMENT: a}), y, w


SUITE_FRAMES = {
    "hand_ordinal": _hand_ordinal_triple,
    "replication_frame": replication_frame,
    "orientation_frame": orientation_frame,
}

# The real O1-O6, captured once at import, so a test that removes ONE check can delegate the other
# five rather than reimplementing them — a reimplemented guard is a second copy of the thing under
# test, and it is how a companion ends up asserting its own behaviour.
_real_assert = model._assert_polr_fittable


# --- 14.3  `polr` and the weights [roadmap, amended] --------------------------------------------------
#
# EVERY ASSERTION IN THIS SECTION IS A TOLERANCE OR A PROPERTY AND NONE IS A VALUE. The measured
# numbers are kept beside them as the record of what was seen, and an implementer who gets a different
# one from these generators has found something rather than broken something.

@pytest.mark.parametrize("name,frame", [
    ("replication_frame", replication_frame),          # Sigma w = 307, measured agreement 0.0
    ("hand_ordinal", _hand_ordinal_triple),            # Sigma w = 42,  measured agreement 1.11e-16
])
def test_INTEGER_WEIGHTS_EQUAL_ROW_REPLICATION(name, frame):
    """§18b's first oracle, and the strongest and cheapest this stage has.

    A weighted likelihood with integer weights IS by definition the unweighted likelihood of the
    frame with each row repeated, so the two fits must agree to machine precision — and no library,
    no subprocess and no tolerance negotiation is involved. **This is the roadmap's weighting
    criterion in the form that can fail**: the all-weights-1 form below cannot, because an
    implementation that never reads `w` at all satisfies it.

    **The tolerance is stated per weight scale, and that is measured rather than stylistic.**
    `POLR_TOL` is ABSOLUTE on a log-likelihood that scales with `Sigma w`, so how closely two fits of
    the same likelihood agree depends on the weight sum. At `Sigma w = 307` they agree at 0.0; at
    `Sigma w = 42` at 1.11e-16; and at `Sigma w = 138` — `hand_ordinal()`'s weights scaled by 20,
    which is mathematically the SAME MLE — at 9.458e-10, which would FAIL this 1e-10 assertion. So
    1e-10 holds on both declared fixtures and is **not** a scale-free guarantee: a future fixture at a
    very different `Sigma w` must re-measure rather than inherit it.
    """
    X, y, w = frame()
    weighted = model.polr(X, y, w)
    replicated = model.polr(*row_replicated(X, y, w))
    assert float(np.max(np.abs(weighted.beta - replicated.beta))) < 1e-10
    assert float(np.max(np.abs(weighted.alpha - replicated.alpha))) < 1e-10
    assert weighted.categories == replicated.categories


def test_all_weights_one_equals_None_BIT_FOR_BIT_which_is_why_it_cannot_be_the_criterion():
    """The companion that makes the oracle above necessary, and the pair IS the test.

    An implementation that never reads `w` at all satisfies the roadmap's "equals an unweighted fit
    when all weights are 1" criterion — trivially, because it fits the unweighted model either way.
    Measured here bit for bit, so the criterion is shown to be unable to fail rather than argued to
    be. `None` is not a synonym for ones in the CALL (§3) and is indistinguishable in the RESULT,
    which is exactly the distinction this asserts.
    """
    X, y, w = _hand_ordinal_triple()
    ones, absent = model.polr(X, y, np.ones(len(y))), model.polr(X, y, None)
    assert np.array_equal(ones.beta, absent.beta)
    assert np.array_equal(ones.alpha, absent.alpha)
    assert (ones.iterations, ones.converged_on) == (absent.iterations, absent.converged_on)


@pytest.mark.parametrize("name,frame", [
    ("hand_ordinal", _hand_ordinal_triple),            # measured 1.8268, OPPOSITE SIGNS
    ("replication_frame", replication_frame),          # measured 0.0566
])
def test_the_weighted_fit_DIFFERS_from_the_unweighted_one_on_the_same_data(name, frame):
    """"Reads `w`" asserted positively rather than only negatively.

    A FLOOR and not a value: the size of the difference is a property of each frame's weights, so a
    magnitude assertion here would pin the fixture rather than the behaviour. On `hand_ordinal()` the
    weights ascend across the treated arm's worsening outcomes and descend across the control's, so
    the weighted and unweighted fits have **opposite signs** — -0.8696 against +0.9573.
    """
    X, y, w = frame()
    assert float(np.max(np.abs(model.polr(X, y, w).beta - model.polr(X, y).beta))) > 1e-3


DERIVATIVE_FRAMES = {
    "hand_ordinal": _hand_ordinal_triple,      # K = 6, m = 1 — the [§8] path's shape
    "replication_frame": replication_frame,    # K = 4, m = 2 — the shape Stage 12 always has
}


def _parameter_points(X, y, w, fit):
    """The three points §14.3 checks the derivatives at, never the optimum alone.

    A wrong sign in the `alpha`-`beta` cross block that happens to VANISH at one parameter value would
    survive a check made only at the optimum — which is exactly where the score is zero and the
    Hessian is best behaved. `par = 0` in the `beta` block with the start cutpoints, the returned
    optimum, and the optimum perturbed by 0.5 in every coordinate.
    """
    Xn, y_idx, w_kept, K = ord_pieces_of(X, y, w, fit)
    share = np.array([w_kept[y_idx <= k].sum() / w_kept.sum() for k in range(K)])
    start = np.concatenate([np.log(share / (1.0 - share)), np.zeros(Xn.shape[1])])
    optimum = np.concatenate([fit.alpha, fit.beta])
    return (Xn, y_idx, w_kept, K), (start, optimum, optimum + 0.5)


@pytest.mark.parametrize("name,frame", list(DERIVATIVE_FRAMES.items()))
def test_the_analytic_score_agrees_with_CENTRAL_DIFFERENCES_at_three_points(name, frame):
    """§18b's fourth oracle. Forty lines of chain rule, and a wrong sign is a plausible number.

    Measured worst |delta score| 4.260e-08 across both frames and all three points, against a
    tolerance of 1e-6 — two orders of margin.
    """
    X, y, w = frame()
    (Xn, y_idx, w_kept, K), points = _parameter_points(X, y, w, model.polr(X, y, w))
    for par in points:
        g, _ = model._ord_score_hess(Xn, y_idx, w_kept, par[:K], par[K:], K)
        numeric = np.zeros_like(par)
        for j in range(len(par)):
            step = np.zeros_like(par)
            step[j] = 1e-5
            numeric[j] = (model._ord_loglik(Xn, y_idx, w_kept, (par + step)[:K], (par + step)[K:], K)
                          - model._ord_loglik(Xn, y_idx, w_kept, (par - step)[:K], (par - step)[K:],
                                              K)) / 2e-5
        assert float(np.max(np.abs(g - numeric))) < 1e-6


@pytest.mark.parametrize("name,frame", list(DERIVATIVE_FRAMES.items()))
def test_the_analytic_hessian_agrees_with_central_differences_and_is_NEGATIVE_DEFINITE(name, frame):
    """The concavity claim as a measurement rather than a citation, at three points on two frames.

    **Symmetric to a TOLERANCE and not by `array_equal`, and that is measured**: `H` is exactly
    symmetric at m = 1 — the [§8] path — and asymmetric by 1.776e-15 at m = 2, because the `beta`
    block is formed as `X'(wH)X` and that matmul is not bitwise symmetric for more than one column.
    An `exactly symmetric` assertion passes on every frame this stage has and fails the first time
    Stage 12 calls `polr` with a real covariate list.

    Measured worst |delta Hessian| 3.790e-08 and worst asymmetry 1.776e-15, with every eigenvalue
    negative at all six points.
    """
    X, y, w = frame()
    (Xn, y_idx, w_kept, K), points = _parameter_points(X, y, w, model.polr(X, y, w))
    for par in points:
        _, H = model._ord_score_hess(Xn, y_idx, w_kept, par[:K], par[K:], K)
        numeric = np.zeros_like(H)
        for j in range(len(par)):
            step = np.zeros_like(par)
            step[j] = 1e-5
            up, _ = model._ord_score_hess(Xn, y_idx, w_kept, (par + step)[:K], (par + step)[K:], K)
            down, _ = model._ord_score_hess(Xn, y_idx, w_kept, (par - step)[:K], (par - step)[K:], K)
            numeric[:, j] = (up - down) / 2e-5
        assert float(np.max(np.abs(H - numeric))) < 1e-6
        assert float(np.max(np.abs(H - H.T))) < 1e-12
        assert float(np.max(np.linalg.eigvalsh(H))) < 0.0


def test_polr_reaches_the_same_optimum_as_an_INDEPENDENT_optimiser():
    """§18b's third oracle, and the only one that would catch a sign error in `_ord_score_hess` that
    happened to be self-consistent.

    Oracle 1 compares two of our own fits and oracle 2 compares a different objective; this checks
    that the Newton loop reaches the maximiser of the function it says it is maximising. Nelder-Mead
    from a deliberately poor start, because a gradient method would use a gradient this is checking.
    `scipy` is TEST-ONLY by policy and this is its second use in the repository (Stage 6 §16b being
    the first). Measured 1.062e-07 on all six parameters with log-likelihoods identical to ten
    decimals.
    """
    from scipy.optimize import minimize
    X, y, w = replication_frame()
    fit = model.polr(X, y, w)
    Xn, y_idx, w_kept, K = ord_pieces_of(X, y, w, fit)
    ours = np.concatenate([fit.alpha, fit.beta])
    # Poor but FEASIBLE: an all-zero start ties every cutpoint, so every category probability is 0,
    # `_ord_loglik` returns -inf across the whole initial simplex and Nelder-Mead has no gradient of
    # information to move on — measured, it exhausts maxfev and returns the start. "Deliberately poor"
    # has to mean far from the optimum, not outside the objective's domain.
    start = np.concatenate([np.linspace(-1.0, 1.0, K), np.zeros(Xn.shape[1])])
    theirs = minimize(lambda p: -model._ord_loglik(Xn, y_idx, w_kept, p[:K], p[K:], K),
                      start, method="Nelder-Mead",
                      options={"xatol": 1e-10, "fatol": 1e-12, "maxfev": 100000, "maxiter": 100000})
    assert theirs.success or theirs.status == 0 or float(np.max(np.abs(ours - theirs.x))) < 1e-6
    assert float(np.max(np.abs(ours - theirs.x))) < 1e-6            # measured 7.417e-08
    assert model._ord_loglik(Xn, y_idx, w_kept, ours[:K], ours[K:], K) >= -float(theirs.fun) - 1e-9


# --- 14.4  the cutpoints, the collapse and the start values -------------------------------------------

@pytest.mark.parametrize("name,frame", list(SUITE_FRAMES.items()))
def test_alpha_is_STRICTLY_ASCENDING_on_every_frame_in_the_suite(name, frame):
    """The post-condition §5.4 argues is guaranteed by construction, asserted on the RETURN.

    Not asserted to be *enforced*: nothing in `polr` checks the ordering and nothing needs to. The
    start values are ordered, a crossed pair gives a non-positive category probability for the
    observations between the two cutpoints, `_ord_loglik` returns -inf, and the step is halved — so no
    accepted iterate is unordered and the returned alpha is ascending because every one before it was.
    That is stronger than an assertion would be: an assertion detects a crossing after the step has
    been taken, and this makes the step unavailable.
    """
    X, y, w = frame()
    fit = model.polr(X, y, w)
    assert np.all(np.diff(fit.alpha) > 0.0)
    assert len(fit.alpha) == len(fit.categories) - 1


def test_the_START_VALUES_are_finite_BECAUSE_of_the_collapse():
    """The coupling stated as an assertion: the collapse and the start values are ONE decision.

    Every kept category carries positive weight, so every weighted cumulative share is strictly
    inside (0, 1) and its logit is finite. The companion computes the same start values over
    `MRS_LEVELS` instead and gets a NON-FINITE cutpoint — which is the fit running a cutpoint to a
    boundary, with the parameter count a property of `config.py` and the identifiability a property
    of the sample.

    **Both signs, and which one you get is a fact about which declared level is unoccupied.** A
    declared level BELOW every occupied one has a cumulative share of exactly 0 and gives `-inf`; one
    ABOVE every occupied one has a share of exactly 1 and gives `+inf`. This frame produces both at
    once, so the companion does not rest on a single direction: `emptied_category_frame`'s response
    occupies 0-4, so declared levels 5 and 6 are above it, and zeroing category 0's weight puts a
    declared level below it.
    """
    (X, y, w), _ = emptied_category_frame(empty=2)
    fit = model.polr(X, y, w)
    Xn, y_idx, w_kept, K = ord_pieces_of(X, y, w, fit)
    share = np.array([w_kept[y_idx <= k].sum() / w_kept.sum() for k in range(K)])
    assert np.all((share > 0.0) & (share < 1.0))
    assert np.all(np.isfinite(np.log(share / (1.0 - share))))

    # the companion: the same arithmetic over the DECLARED level set, which 5 and 6 are absent from
    zeroed = np.where(y == 0.0, 0.0, w)
    declared = np.array([zeroed[y <= float(level)].sum() / zeroed.sum()
                         for level in config.MRS_LEVELS[:-1]])
    with np.errstate(divide="ignore"):
        over_declared = np.log(declared / (1.0 - declared))
    assert not np.all(np.isfinite(over_declared))
    assert float(np.min(over_declared)) == -np.inf          # level 0, below every occupied one
    assert float(np.max(over_declared)) == np.inf           # levels 5 and 6, above every one
    # and the collapse's own start values on the SAME weights stay finite, which is the contrast
    collapsed = model.polr(X, y, zeroed)
    assert collapsed.categories == (1, 2, 3, 4)
    assert np.all(np.isfinite(collapsed.alpha))


def test_a_ZERO_WEIGHT_category_is_dropped_and_the_fit_equals_the_rows_removed_one():
    """The one input separating §5.3's rule from the observed-set rule.

    `np.unique(y)` would keep a category all of whose records carry zero weight: it contributes
    nothing to the likelihood, so its two cutpoints are unidentified — and worse, §5.4's crossing
    argument fails there, because a crossing only drives a category probability non-positive FOR
    OBSERVATIONS IN THAT CATEGORY, and if they all have zero weight the -inf never fires.

    **The two fits do NOT share an O6 path**, and that is stated so nobody reads this assertion as
    covering the collapsed design's rank: O6 runs BEFORE the collapse, so the left-hand fit's rank was
    checked on the full design and the right-hand fit's on the reduced one. §15 files what that leaves
    uncovered — a caller passing a weight of exactly zero can collapse away every row of a design
    column and get `np.linalg.solve`'s `LinAlgError`, which is not `FitError`.
    """
    (X, y, w), (X_removed, y_removed, w_removed) = emptied_category_frame(empty=2)
    zeroed = np.where(y == 2.0, 0.0, w)
    collapsed = model.polr(X, y, zeroed)
    removed = model.polr(X_removed, y_removed, w_removed)
    assert 2 not in collapsed.categories
    assert collapsed.categories == removed.categories == (0, 1, 3, 4)
    assert float(np.max(np.abs(collapsed.beta - removed.beta))) < 1e-10
    assert float(np.max(np.abs(collapsed.alpha - removed.alpha))) < 1e-10


def test_the_collapse_and_the_full_frame_differ_in_their_cutpoint_count():
    # The other half: without it, the assertion above would pass on a frame where nothing collapsed.
    (X, y, w), _ = emptied_category_frame(empty=2)
    assert len(model.polr(X, y, w).alpha) == 4
    assert len(model.polr(X, y, np.where(y == 2.0, 0.0, w)).alpha) == 3


def test_weighted_categories_is_the_rule_and_is_callable_directly():
    y = np.array([0., 1., 2., 3.])
    assert model._weighted_categories(y, np.array([1., 1., 1., 1.])) == (0, 1, 2, 3)
    assert model._weighted_categories(y, np.array([1., 0., 1., 1.])) == (0, 2, 3)
    # and the nan fact the whole ordering of O3 and O5 turns on, asserted where the rule lives
    assert (np.nan > 0.0) is False
    assert model._weighted_categories(y, np.array([np.nan, 1., 1., 1.])) == (1, 2, 3)


def _unique_on_the_bare_response(source: str) -> list[tuple[str, int]]:
    """Every `np.unique(y)` call — the CATEGORY RULE's expression — with its enclosing function.

    Scoped to the argument and not to the name, and the scoping is the point. `np.unique` legitimately
    appears three more times in this module, all of them formatting a message: `_assert_fittable`'s F4
    and `_assert_polr_fittable`'s O1 call it on `np.argwhere(...)` output to name the offending
    columns, and O4 calls it on the off-domain values to list them. A scan keyed on the bare name
    would forbid those and the repair would be to gut the scan (Stage 6 §12.7's rule).

    What §5.3a is actually about is one expression: `np.unique(y)`, the positive-weight level set's
    own subject. Two copies of it are one edit from disagreeing about which levels a fit has, and the
    disagreement is silent in both directions — O5 counts one number of categories and `polr` fits
    another, both finite and both plausible.
    """
    found: list[tuple[str, int]] = []
    for parent in ast.walk(ast.parse(source)):
        if not isinstance(parent, (ast.FunctionDef, ast.AsyncFunctionDef)):
            continue
        for node in ast.walk(parent):
            if (isinstance(node, ast.Call) and isinstance(node.func, ast.Attribute)
                    and node.func.attr == "unique" and len(node.args) == 1
                    and isinstance(node.args[0], ast.Name) and node.args[0].id == "y"):
                found.append((parent.name, node.lineno))
    return found


def test_the_positive_weight_rule_is_ONE_function_and_not_one_expression_written_twice():
    # §5.3a. Without this the two copies come back on the next edit and nothing fails.
    assert [name for name, _ in _unique_on_the_bare_response(SOURCE)] == ["_weighted_categories"]
    assert "np.unique" not in ast.get_source_segment(SOURCE, _polr_ast()) 


def _polr_ast():
    return next(n for n in ast.walk(ast.parse(SOURCE))
                if isinstance(n, ast.FunctionDef) and n.name == "polr")


def test_that_duplication_scan_fires():
    pasted = ("def polr(X, y, w):\n"
              "    categories = tuple(int(v) for v in np.unique(y) if w[y == v].sum() > 0.0)\n")
    assert [name for name, _ in _unique_on_the_bare_response(pasted)] == ["polr"]


def test_the_MINUS_INF_branch_of_ord_loglik_returns_minus_inf_on_a_CROSSED_alpha():
    """Asserted directly on the private, in the posture Stage 7 §12.11 asserts B4: unreachable on the
    declared path, tested by calling it.

    **No input reaches it through `polr`**, and that is a measured negative result rather than an
    assumption: 3000 attempts at reaching a crossing with pure Newton from a deliberately crowded
    start — the cutpoints initialised to a span of 2e-3 — never lost the ordering. An earlier draft
    of the spec claimed the crossing was reachable by removing the halving; it is not, and the claim
    is struck rather than weakened.
    """
    X, y, w = _hand_ordinal_triple()
    fit = model.polr(X, y, w)
    Xn, y_idx, w_kept, K = ord_pieces_of(X, y, w, fit)
    assert np.isfinite(model._ord_loglik(Xn, y_idx, w_kept, fit.alpha, fit.beta, K))
    crossed = fit.alpha.copy()
    crossed[0], crossed[1] = crossed[1] + 1.0, crossed[0] - 1.0      # a crossed PAIR
    assert model._ord_loglik(Xn, y_idx, w_kept, crossed, fit.beta, K) == -np.inf


@pytest.mark.parametrize("name,frame", list(SUITE_FRAMES.items()))
def test_POLR_ETA_CLIP_is_never_approached_on_any_frame_here(name, frame):
    """The honest form of coverage for a branch no input reaches, and a REAL assertion.

    Measured: the largest |linear predictor| this estimator has produced is 18.2, in the perfectly
    separated construction, against a clip of 500. Asserted against the constant rather than against
    500, because the point is the relationship.
    """
    X, y, w = frame()
    fit = model.polr(X, y, w)
    Xn, y_idx, w_kept, K = ord_pieces_of(X, y, w, fit)
    eta = np.abs(Xn @ fit.beta) + float(np.max(np.abs(fit.alpha)))
    assert float(np.max(eta)) < config.POLR_ETA_CLIP / 10.0


@pytest.mark.parametrize("name,frame", list(SUITE_FRAMES.items()))
def test_converged_on_is_declared_and_iterations_is_CLOSED_at_the_top(name, frame):
    # Closed at the top because `range(1, POLR_MAX_ITER + 1)` can return POLR_MAX_ITER, and a
    # half-open assertion would fail on the last iteration that still converges.
    fit = model.polr(*frame())
    assert fit.converged_on in ("likelihood", "score")
    assert 1 <= fit.iterations <= config.POLR_MAX_ITER


def test_the_TRUST_REGION_fires_in_isolation_from_step_halving():
    # Asserted as VALUES on a frame that drives the counter, not as an integer on a frame where it is
    # zero: measured, with the type-only form BOTH safeguards can be deleted from `polr` and every
    # frame in the suite returns bit-identically.
    fit = model.polr(*polr_big_first_step())
    assert (fit.rescales, fit.halvings) == (1, 0)
    assert fit.first_step_norm > config.POLR_MAX_STEP


def test_STEP_HALVING_fires_together_with_the_trust_region_at_a_smaller_scale():
    fit = model.polr(*polr_needs_halving())
    assert (fit.rescales, fit.halvings) == (2, 1)
    assert fit.first_step_norm > config.POLR_MAX_STEP


@pytest.mark.parametrize("name,frame", [
    *SUITE_FRAMES.items(),
    ("separated_frame", separated_frame),
])
def test_neither_safeguard_fires_on_a_well_behaved_or_a_separated_frame(name, frame):
    fit = model.polr(*frame())
    assert (fit.rescales, fit.halvings) == (0, 0)


def test_the_two_safeguard_fixtures_are_also_the_SCALE_INVARIANCE_test():
    """A proportional-odds estimate is equivariant under rescaling a covariate, which is §5.2's own
    argument for a RELATIVE trust radius — so the same pair does double duty.

    1e-8 and not tighter, because `POLR_TOL` is absolute and the three fits therefore stop at slightly
    different points. That looseness is itself the measurement (§5.2), not slack.
    """
    betas = [float(model.polr(*scaled_exposure(scale)).beta[0]) * scale
             for scale in (1.0, 0.1, 0.01)]
    assert max(abs(b - betas[0]) for b in betas) < 1e-8


def test_POLR_MAX_HALVINGS_exhausted_says_there_is_no_second_estimator(monkeypatch):
    # Stage 6 §12.4a's fourth assertion, one stage on. Without it the `for…else` branch of §5.2 is
    # unreachable in the suite.
    monkeypatch.setattr(config, "POLR_MAX_HALVINGS", 1)
    with pytest.raises(model.FitError) as e:
        model.polr(*polr_needs_halving())
    assert "halvings at iteration" in message_of(e)
    assert "there is no second one to try" in message_of(e)          # invariant 5


def test_POLR_MAX_ITER_exhausted_reports_a_NON_ZERO_movement(monkeypatch):
    monkeypatch.setattr(config, "POLR_MAX_ITER", 2)
    with pytest.raises(model.FitError) as e:
        model.polr(*polr_needs_halving())
    text = message_of(e)
    assert "no convergence in 2 iterations" in text
    assert float(text.split("moved ")[1].split(" ")[0]) > 0.0
    assert "score component" in text and "shortened by the trust region" in text


@pytest.mark.parametrize("name,frame", list(SUITE_FRAMES.items()))
def test_first_step_norm_is_positive_finite_and_INFORMATIVE(name, frame):
    """`first_step_norm` is the one field of the four that `_fit_detail` does not print, so this
    assertion is the only thing keeping it alive — which is why it is named rather than folded into
    "the counters".

    Strictly-positive alone is too weak: measured, a frame whose arms have identical weighted
    distributions returns `iterations = 1`, `beta = 0.0` and `first_step_norm = 1.86e-16`, a value
    that passes a positivity assertion and carries no information. The floor is therefore conditional
    on the fit having taken more than one iteration.
    """
    fit = model.polr(*frame())
    assert np.isfinite(fit.first_step_norm) and fit.first_step_norm > 0.0
    if fit.iterations > 1:
        assert fit.first_step_norm > 1e-8


def test_that_first_step_norm_floor_would_NOT_hold_on_a_one_iteration_fit():
    # The measurement behind the conditional above, asserted so the condition is not mistaken for
    # timidity: two arms with identical weighted distributions converge at iteration 1 with beta 0.
    X = pd.DataFrame({config.TREATMENT: np.array([1., 1., 1., 0., 0., 0.])})
    fit = model.polr(X, np.array([0., 1., 2., 0., 1., 2.]), np.ones(6))
    assert fit.iterations == 1 and fit.converged_on == "likelihood"
    assert float(fit.beta[0]) == pytest.approx(0.0, abs=1e-12)
    assert 0.0 < fit.first_step_norm < 1e-8


# --- 14.7  separation converges, and G7 is what turns it into a failure --------------------------------
#
# THE SECTION THAT CARRIES THIS STAGE'S FINDING. `polr`'s half is here; the guard's half is
# `test_outcome.py`'s, because G7 is the caller's.
#
# ITS BUDGET IS STATED because an earlier draft did not have one: three of its assertions ranged over a
# CALIBRATION rather than a sample — a 4800-fit sweep, 500 healthy cond(-H) fits and 400 sparse
# replicates run twice, about 6100 fits and 14 seconds, on every commit for the remaining six stages.
# The split is: calibrations live in the spec's §20 and in `test_config.py`'s literal, probes live
# here. Everything below is about 300 fits and under a second, and the calibration's endpoints are
# guarded by `test_POLR_MAX_ABS_BETA_is_strictly_inside_the_measured_sparse_region`, which is the one
# assertion that actually catches their regression.

def test_A_PERFECTLY_SEPARATED_FIT_CONVERGES_AND_RETURNS():
    """The finding, asserted as VALUES — and it can be, because `separated_frame()` holds no seed.

    Nothing in this block is out of range. It is not a rank failure, not a halving exhaustion and not
    a non-convergence: it converged, on the score criterion, in 17 iterations, with every safeguard
    counter at zero and every fitted probability finite. The likelihood approaches its supremum as
    beta grows and FLATTENS; the score falls below tolerance; the loop returns. `exp(beta)` is 6.5e15.

    A tolerance would be a weaker way of saying "nothing looks wrong", which is the whole point.
    """
    X, y, w = separated_frame()
    fit = model.polr(X, y, w)
    assert fit.converged_on == "score"
    assert fit.iterations == 17
    assert (fit.rescales, fit.halvings) == (0, 0)
    assert float(np.max(np.abs(fit.beta))) > 19.0
    assert np.all(np.isfinite(fit.beta)) and np.all(np.isfinite(fit.alpha))
    p = category_probabilities(X, y, w, fit)
    assert np.all(np.isfinite(p)) and np.all((p > 0.0) & (p <= 1.0))


def test_separation_scales_with_n_rather_than_being_an_artefact_of_a_tiny_frame():
    betas = [float(np.max(np.abs(model.polr(*separated_frame(per_arm=n)).beta)))
             for n in (10, 20, 40, 80)]
    # Non-decreasing to a tolerance and NOT `== sorted(betas)`: measured 34.41 / 36.41 / 38.41 /
    # 38.41, and the last two agree only to 2.9e-14 — the likelihood is flat out there, so an exact
    # ordering assertion tests float noise on a ridge rather than the scaling.
    assert all(later >= earlier - 1e-9 for earlier, later in zip(betas, betas[1:]))
    assert min(betas) > 30.0
    assert [round(b, 2) for b in betas] == [34.41, 36.41, 38.41, 38.41]


def test_the_NON_CONVERGENCE_detector_does_not_fire_ever():
    # Detector 1 of the four §6.3 designed and measured. It does not discriminate because it does not
    # fire at all: both a degenerate fit and a healthy one converge.
    assert model.polr(*separated_frame()).iterations <= config.POLR_MAX_ITER
    assert model.polr(*cohort_shaped_frame(0.5)).iterations <= config.POLR_MAX_ITER


def test_the_FITTED_PROBABILITY_detector_catches_PERFECT_and_misses_NEAR_separation():
    """Detector 3, and the pair is the reason §6.3's guard is a magnitude bound.

    At perfect separation every fitted category probability is 1.0 and the detector fires. At ONE
    crossover patient out of twenty the probabilities are unremarkable — min 0.05 — while `exp(beta)`
    is still 1.5e+09. **Near separation is the case a bootstrap actually draws**, so a detector that
    only catches perfect separation catches the case that does not arise and misses the case that
    does.
    """
    X, y, w = separated_frame()
    perfect = category_probabilities(X, y, w, model.polr(X, y, w))
    assert float(np.min(perfect)) == pytest.approx(1.0, abs=1e-6)

    Xn, yn, wn = separated_frame(crossovers=1)
    near = model.polr(Xn, yn, wn)
    assert float(np.min(category_probabilities(Xn, yn, wn, near))) > 0.01
    assert float(np.exp(np.max(np.abs(near.beta)))) > 1e8


def test_the_COND_H_detector_DOES_discriminate_at_equal_cutpoint_count_and_is_rejected_anyway():
    """Detector 2 — and the whole of this test is a CORRECTION.

    An earlier draft wrote "6.85 against 3 to 40, does not discriminate", and both halves were wrong:
    the 6.854 is real but it is a **2x2** matrix's, `separated_frame()` having two outcome categories
    and therefore one cutpoint, while the healthy fits it was compared against are 7x7. That is not a
    comparison. Measured properly, at equal cutpoint count, it discriminates by six orders of
    magnitude.

    **Asserted as orders of magnitude and not as ranges**, because the healthy range is a property of
    the sampler: measured 29.1 to 162.8 healthy, 2.8e8 near-separated, 6.854 on the two-category
    construction.

    It is still not the guard, for two reasons that survive the correction. The first is the same size
    dependence that produced the error — `cond(-H)` is a property of a matrix whose dimension is
    `len(alpha) + len(beta)`, which §5.3 makes a property of the SAMPLE — and the third assertion
    below is that dependence as a measurement: the two-category separated construction reads BELOW
    every healthy seven-category fit, so a "cond above a bound" rule fires on nothing there. The
    second is that a condition-number threshold cannot be justified against anything a reader can
    evaluate, where a bound on the reported coefficient can: `exp(beta)` is the number the manuscript
    prints.
    """
    healthy_X, healthy_y, healthy_w = cohort_shaped_frame(0.5)
    healthy_fit = model.polr(healthy_X, healthy_y, healthy_w)
    healthy = cond_of(healthy_X, healthy_y, healthy_w, healthy_fit)

    sep_X, sep_y, sep_w = cohort_shaped_frame(6.0)
    sep_fit = model.polr(sep_X, sep_y, sep_w)
    near = cond_of(sep_X, sep_y, sep_w, sep_fit)

    assert len(healthy_fit.alpha) == len(sep_fit.alpha) == 6      # EQUAL cutpoint count, or nothing
    assert healthy < 1e4                                          # measured 58.4
    assert near > 1e6                                             # measured 2.8e8
    assert float(np.max(np.abs(sep_fit.beta))) > config.POLR_MAX_ABS_BETA

    # the size dependence, which is why it is rejected: a 2x2's condition number is not comparable
    two_category_X, two_category_y, two_category_w = separated_frame()
    two_category = cond_of(two_category_X, two_category_y, two_category_w,
                           model.polr(two_category_X, two_category_y, two_category_w))
    assert two_category < healthy


def test_THE_BAND_IS_A_BAND_at_the_probes_budget_and_not_the_calibrations():
    """§6.3's distribution asserted as a band rather than as two endpoints, at 5% of the fits.

    240 fits in about 0.7 s across twelve true effect sizes. §6.3's 4800-fit calibration is NOT re-run
    here: its endpoints are recorded in the spec's §20 and pinned as literals in `test_config.py`,
    which is what fails if the bound is edited back to the 10.0 an earlier draft used — the failure
    mode the 4800 fits were guarding against, caught by one comparison instead of fourteen seconds of
    fitting.

    Measured on this probe: 153 below the bound and 87 at or above it, max below 8.5534, min above
    20.7879, and nothing at all within +-4 of the bound.
    """
    fits = band_samples(BAND_EFFECTS, per_effect=20)
    magnitudes = np.array([float(np.max(np.abs(f.beta))) for f in fits])
    below = magnitudes[magnitudes < config.POLR_MAX_ABS_BETA]
    above = magnitudes[magnitudes >= config.POLR_MAX_ABS_BETA]
    assert below.size and above.size                       # both modes populated
    assert int(((magnitudes > config.POLR_MAX_ABS_BETA - 4.0)
                & (magnitudes < config.POLR_MAX_ABS_BETA + 4.0)).sum()) == 0
    assert float(below.max()) < config.POLR_MAX_ABS_BETA < float(above.min())


def test_the_band_probe_contains_fits_with_a_NON_ZERO_safeguard_counter():
    """A free witness, and it strikes a sentence an earlier draft carried.

    That draft said "zero of both counters on every frame in §20, which is what makes a non-zero
    counter in a future run a signal rather than noise". The suite's own probe contradicts it:
    measured, 10 of these 240 fits carry one.
    """
    fits = band_samples(BAND_EFFECTS, per_effect=20)
    assert sum(1 for f in fits if f.rescales or f.halvings) > 0


def test_SPARSE_REPLICATES_produce_ZERO_convergence_failures_and_that_is_the_point():
    """§11's warning to Stage 10, made concrete rather than cautionary.

    The two counters mean different things about the data, and that does not need four hundred
    samples to state: **48 replicates here, not 400** — four draws at each of the twelve declared true
    effect sizes; the measurement over 400 belongs in the spec's §20. What is asserted is that the
    non-convergence count is exactly 0 while the count of fits that a G7-shaped bound would reject is
    greater than 0 — measured 0 and 10. A Stage 10 failure counter reading zero is therefore NOT
    evidence that no replicate was degenerate unless the guard is in the path.

    **The effect range is what makes the second counter non-empty**, which is why it ranges over all
    twelve rather than over the low end: at true effects of 0 to 1.5 nothing separates, so a probe
    restricted to those would assert `degenerate > 0` and fail on correct code.
    """
    failures, degenerate = 0, 0
    for effect in BAND_EFFECTS:
        for draw in range(4):
            rng = np.random.default_rng(config.SEED + draw)
            a = np.concatenate([np.ones(39), np.zeros(53)])
            latent = effect * a + rng.normal(size=92)
            cuts = np.quantile(latent, np.linspace(0.0, 1.0, 8)[1:-1])
            y = np.searchsorted(cuts, latent).astype(float)
            try:
                fit = model.polr(pd.DataFrame({config.TREATMENT: a}), y,
                                 rng.uniform(0.05, 0.95, size=92))
            except model.FitError as raised:
                failures += "no convergence" in str(raised)
                continue
            degenerate += float(np.max(np.abs(fit.beta))) >= config.POLR_MAX_ABS_BETA
    assert failures == 0
    assert degenerate > 0


# --- 14.8  O1-O6, and the two whose order is the point -------------------------------------------------

def test_O1_fires_on_a_non_finite_design_and_names_the_column():
    X, y, w = _hand_ordinal_triple()
    X = X.copy()
    X.loc[0, config.TREATMENT] = np.nan
    with pytest.raises(model.FitError) as e:
        model.polr(X, y, w)
    assert "O1" in message_of(e) and config.TREATMENT in message_of(e)


def test_O2_fires_on_a_non_finite_response():
    X, y, w = nan_response_frame()
    y = y.copy()
    y[3] = np.nan
    with pytest.raises(model.FitError) as e:
        model.polr(X, y, w)
    assert "O2" in message_of(e) and "non-finite" in message_of(e)


def test_WITHOUT_O2_AND_O4_a_nan_response_silently_DROPS_the_record(monkeypatch):
    """The mechanism asserted, not just the raise — and the companion must remove O2 **and** O4.

    `np.mod(nan, 1.0)` is `nan` and `nan != 0.0` is True, so **O4 catches a non-finite response too**,
    and the companion as first drafted would have asserted a silent drop while watching O4 raise.
    That interaction is worth stating in its own right: O2 is not the sole detector of a nan response,
    it is the one with the right message and the earlier position, and O4 is a backstop nobody
    designed as one.

    With both removed: 79 of 80 records fitted, `categories` UNCHANGED, and beta moves — the record is
    dropped, not reported. `np.unique` returns nan as a level, `w[y == nan]` is empty so its total
    weight is 0.0, and the positive-weight rule drops it with nothing missing anywhere in the output.
    """
    assert np.isnan(np.mod(np.nan, 1.0)) and (np.mod(np.nan, 1.0) != 0.0)
    X, y, w = nan_response_frame()
    clean = model.polr(X, y, w)
    damaged = y.copy()
    damaged[3] = np.nan

    # O4 alone still fires, which is what makes the two-removal companion necessary
    monkeypatch.setattr(model, "_assert_polr_fittable",
                        lambda Xn, yv, wv, cols: None if not np.all(np.isfinite(yv)) else
                        _real_assert(Xn, yv, wv, cols))
    fit = model.polr(X, damaged, w)
    assert fit.categories == clean.categories
    assert float(np.max(np.abs(fit.beta - clean.beta))) > 1e-3
    kept = np.isin(damaged, np.asarray(fit.categories, dtype=float))
    assert int(kept.sum()) == len(y) - 1                        # 79 of 80, and nothing missing


def test_removing_O2_ALONE_does_not_reach_a_fit_because_O4_catches_it_too(monkeypatch):
    """The interaction stated in its own right, and it is why the companion above removes BOTH.

    O2 is the better message and the earlier position; it is not the sole detector. O4 is a backstop
    nobody designed as one.
    """
    def without_O2(Xn, y, w, columns):
        if not np.all(np.isfinite(Xn)):
            raise model.FitError("O1  the design carries a non-finite value in: —")
        if not np.all(np.isfinite(w)) or np.any(w < 0.0):
            raise model.FitError("O3  the weights carry a non-finite value")
        if y[np.not_equal(np.mod(y, 1.0), 0.0)].size:
            raise model.FitError("O4  the response takes non-integer value(s)")

    monkeypatch.setattr(model, "_assert_polr_fittable", without_O2)
    X, y, w = nan_response_frame()
    y = y.copy()
    y[3] = np.nan
    with pytest.raises(model.FitError) as e:
        model.polr(X, y, w)
    assert "O4" in message_of(e)


def test_O3_fires_on_a_nan_weight_BEFORE_O5_can_delete_the_category():
    X, y, w = nan_weight_frame()
    w = w.copy()
    w[0] = np.nan
    with pytest.raises(model.FitError) as e:
        model.polr(X, y, w)
    assert "O3" in message_of(e) and "DELETES an outcome" in message_of(e)


def test_WITHOUT_O3_one_nan_weight_DELETES_A_CATEGORY_and_moves_beta(monkeypatch):
    """The ordering asserted as an ordering, and BOTH halves are required.

    The first alone passes for an implementation that raises somewhere; this is what shows what it is
    raising *instead of*. `np.nan > 0.0` is False, so the category's total weight is nan, the
    positive-weight test drops it, and the fit comes back over a COARSER SCALE — finite, plausible,
    with a clean iteration count and nothing raised.

    Measured: four categories became three and beta moved from -1.6789 to +0.1949 on this frame.
    """
    X, y, w = nan_weight_frame()
    clean = model.polr(X, y, w)
    damaged = w.copy()
    damaged[0] = np.nan

    def without_O3(Xn, yv, wv, columns):
        return _real_assert(Xn, yv, np.where(np.isfinite(wv), wv, 1.0), columns)

    monkeypatch.setattr(model, "_assert_polr_fittable", without_O3)
    with np.errstate(invalid="ignore"):
        fit = model.polr(X, y, damaged)
    assert len(fit.categories) == len(clean.categories) - 1
    assert clean.categories == (0, 1, 2, 3) and fit.categories == (1, 2, 3)   # w[0] sits at y = 0
    assert round(float(clean.beta[0]), 6) == 1.678900
    assert round(float(fit.beta[0]), 6) == 0.194854
    assert fit.converged_on == "likelihood"        # nothing looks wrong


def test_O3_fires_on_a_NEGATIVE_weight_too():
    X, y, w = nan_weight_frame()
    w = w.copy()
    w[0] = -0.5
    with pytest.raises(model.FitError) as e:
        model.polr(X, y, w)
    assert "O3" in message_of(e) and "negative" in message_of(e)


def test_O4_fires_on_a_non_integer_response():
    with pytest.raises(model.FitError) as e:
        model.polr(*noninteger_response_frame())
    assert "O4" in message_of(e) and "non-integer" in message_of(e)


def test_WITHOUT_O4_one_cutpoint_is_estimated_PER_DISTINCT_VALUE(monkeypatch):
    """The generator is sized so the companion reaches a FIT: three distinct values over sixty
    records, not the wider frame an earlier draft implied.

    With O4 removed, six distinct half-integer values over 24 records makes the Hessian singular and
    `np.linalg.solve` raises `LinAlgError` — which is not `FitError`, so the companion would assert
    the wrong thing (measured).
    """
    monkeypatch.setattr(
        model, "_assert_polr_fittable",
        lambda Xn, y, w, columns: _real_assert(Xn, np.round(y * 2.0), w, columns))
    X, y, w = noninteger_response_frame()
    fit = model.polr(X, y, w)
    assert len(fit.categories) == len(np.unique(y)) == 3
    assert len(fit.alpha) == 2


def test_O5_fires_on_a_response_with_ONE_weighted_category():
    X = pd.DataFrame({config.TREATMENT: np.array([1., 0., 1., 0.])})
    with pytest.raises(model.FitError) as e:
        model.polr(X, np.array([2., 2., 2., 2.]), np.ones(4))
    assert "O5" in message_of(e) and "1 response category" in message_of(e)


def test_O6_fires_on_a_design_carrying_ITS_OWN_INTERCEPT():
    """The asymmetry with `firth`, which PREPENDS one — the one place a reader will assume the two
    functions agree. The cutpoints ARE the intercepts, so a design that carries one of its own is
    over-parametrised.

    **The witness is an intercept beside a COMPLETE set of arm dummies**, and that is a correction
    found by running the check rather than by writing it: `matrix_rank` sees the DESIGN and not the
    (alpha, beta) parameter vector, so it can only detect an intercept that is collinear with other
    design columns. `[1, treated, control]` is rank 2 of 3 and fires; the companion below measures
    what happens to `[1, treated]`, which does not.
    """
    y, a, w = hand_ordinal()
    X = pd.DataFrame({"intercept": np.ones(len(y)), "treated": a, "control": 1.0 - a})
    with pytest.raises(model.FitError) as e:
        model.polr(X, y, w)
    assert "O6" in message_of(e) and "rank 2 of 3" in message_of(e)
    assert "opposite of the penalised logistic fit" in message_of(e)


def test_what_O6_does_NOT_catch_and_what_it_costs_measured():
    """An intercept beside ONE other column is full rank, so O6 passes and `polr` fits.

    Recorded rather than repaired, because the measurement says what it costs and the answer is
    "the cutpoints, not the estimate". The design `[1, a]` makes the intercept collinear with a
    uniform SHIFT of all K cutpoints, which is a dependency in the (alpha, beta) parameter vector and
    not in the design matrix — so `matrix_rank(Xn)` cannot see it. Measured: `np.linalg.solve` does
    NOT raise on the resulting near-singular Hessian; the loop converges to an arbitrary point on the
    ridge, the cutpoints come back shifted by the intercept coefficient, and **the treatment
    coefficient is unchanged to every printed digit**.

    So the failure is confined to parameters this pipeline does not report as effects, and the [§8]
    path cannot reach it at all: `primary` builds its design through `model.design(sub, (TREATMENT,))`
    and `_assert_exposure_survived` asserts the columns are exactly `(TREATMENT,)`, so an extra
    intercept column raises G6 before `polr` is called. It is filed here so nobody reads O6's message
    as a stronger guarantee than O6 gives.
    """
    y, a, w = hand_ordinal()
    plain = model.polr(pd.DataFrame({config.TREATMENT: a}), y, w)
    with_intercept = model.polr(
        pd.DataFrame({"intercept": np.ones(len(y)), config.TREATMENT: a}), y, w)
    assert float(with_intercept.beta[1]) == pytest.approx(float(plain.beta[0]), abs=1e-8)
    shift = float(with_intercept.beta[0])
    assert abs(shift) > 1e-3                                   # the ridge was travelled
    assert np.allclose(with_intercept.alpha + shift, plain.alpha, atol=1e-8)


def test_O6_does_NOT_fire_on_a_width_zero_design():
    """Rank 0 of 0 is not rank-deficient — which is what makes `_assert_exposure_survived` necessary
    on the caller's side rather than redundant.

    A module that does not know what the exposure is cannot know that the missing column was the
    estimand, so the check belongs where the specification does. `polr` on the width-0 design fits the
    cutpoints alone, CONVERGES, and returns a beta of length zero.
    """
    y, a, w = hand_ordinal()
    X = pd.DataFrame(index=range(len(y)))
    fit = model.polr(X, y, w)                      # does not raise
    assert fit.beta.shape == (0,)
    assert len(fit.alpha) == 6 and fit.columns == ()


def test_the_six_checks_raise_at_the_FIRST_failure_and_are_never_collected():
    # Unlike D1-D4 and B1-B5: each one makes the next meaningless, and a rank computed over a design
    # containing nan is not a rank. A frame failing O1 and O6 reports O1 only.
    y, a, w = hand_ordinal()
    X = pd.DataFrame({"intercept": np.ones(len(y)), config.TREATMENT: a})
    X.loc[0, config.TREATMENT] = np.nan
    with pytest.raises(model.FitError) as e:
        model.polr(X, y, w)
    assert "O1" in message_of(e)
    assert "O6" not in message_of(e)


@pytest.mark.parametrize("name", ["O1", "O2", "O3", "O4", "O5", "O6"])
def test_no_polr_precondition_message_names_an_exposure_or_an_outcome(name):
    # §14.13, parametrised over the six. `model.py` is outcome-agnostic and its raise strings are
    # where that most easily stops being true.
    ordinal_raises = [s for s in _raise_strings(SOURCE) if s.strip().startswith(name)]
    assert ordinal_raises, f"{name} has no raise string to scan"
    for text in ordinal_raises:
        assert not [word for word in FORBIDDEN_IN_RAISES if word in text]


# --- 15.6  `predict` — Stage 9 §15.6, and it lives here because the function does -----------------------
#
# `model.Fit.p` is the fitted probability at the OBSERVED treatment, because the X handed to `firth`
# carries the observed treatment column. It is m_A(X) — the counterfactual matching each row's actual
# arm — so it is NOT m_1 and NOT m_0, and an implementation reading it for either is wrong on every row
# assigned to the other arm. `predict` is why that is not necessary [Stage 9 §7.3].
#
# THE ORACLE IS AN EQUALITY AND NOT A TOLERANCE, and that is the whole design of this section: both of
# §7.3's details — the intercept prepended as a COLUMN, and the FIRTH_ETA_CLIP clip — are invisible to
# `allclose` and visible to `array_equal`. A test written with a tolerance passes on both defects.

from fixtures_stage9 import golden_arrays, golden_frame                        # noqa: E402
from fixtures_stage9 import separated_frame as separated_binary_frame          # noqa: E402
from fixtures_stage12 import (RI_BETA_TRUE, RI_SIGMA_TRUE,                    # noqa: E402
                              collapsed_level_frame, flat_ri_frame, two_group_ri_frame)

# THE ALIAS IS NECESSARY AND IS NOT TIDINESS. This module already defines `separated_frame` — Stage 8's
# perfectly separated ORDINAL frame for `polr` — and Stage 9's is a separated BINARY frame for `firth`.
# Same word, same idea, two fitters, two frames. An unaliased import would silently rebind this file's
# own fixture and every §14.7 test below would start asserting Stage 9's numbers.

GOLDEN_COVARIATE_PAIRS: Final[tuple[tuple[str, ...], ...]] = (
    ("center", "atrial_fib"), ("center", "atrial_fib", "age"))


def golden_fit(covariates: tuple[str, ...]):
    """(fit, X) for `golden_frame()` under one of §15.0.2's two covariate lists.

    The design is built exactly as `outcome.outcome_model` builds it — `model.design`, then treatment
    inserted first by name — because the oracle below is that `predict` reproduces `fit.p` on THE VERY
    X THE FIT WAS MADE FROM, and a differently assembled X would be testing a different claim.
    """
    df = golden_frame()
    y, a, _, _ = golden_arrays(df)
    X, _ = model.design(df, covariates)
    X = X.copy()
    X.insert(0, config.TREATMENT, a)
    return model.firth(X, y), X


def test_predict_is_the_inverse_logit_on_a_hand_computed_three_column_design():
    """Against arithmetic done by hand, not against another implementation.

    Two rows, three columns and a beta of (1, 2, 3, 4) chosen so the linear predictors are integers
    and the expected probabilities are writable: eta = 1 + 2*x1 + 3*x2 + 4*x3.
    """
    fit = model.Fit(beta=np.array([1.0, 2.0, 3.0, 4.0]), p=np.array([np.nan, np.nan]),
                    iterations=1, converged_on="likelihood", columns=("x1", "x2", "x3"),
                    first_step_norm=0.0, rescales=0, halvings=0)
    X = pd.DataFrame({"x1": [0.0, 1.0], "x2": [0.0, -1.0], "x3": [0.0, 0.5]})
    eta = np.array([1.0, 1.0 + 2.0 - 3.0 + 2.0])
    assert model.predict(fit, X) == pytest.approx(1.0 / (1.0 + np.exp(-eta)))


@pytest.mark.parametrize("covariates", GOLDEN_COVARIATE_PAIRS)
def test_predict_reproduces_fit_p_BIT_FOR_BIT_on_the_design_the_fit_was_made_from(covariates):
    """`np.array_equal`, not `allclose`, and Stage 9 §15.6 is emphatic about the difference.

    This is the only available oracle that `predict` and `firth` agree about the link, and it pins
    BOTH of §7.3's details at once. The two companions below are what make the strength of the
    assertion load-bearing rather than fussy.
    """
    fit, X = golden_fit(covariates)
    assert np.array_equal(model.predict(fit, X), fit.p)


@pytest.mark.parametrize("covariates", GOLDEN_COVARIATE_PAIRS)
def test_WITHOUT_the_prepended_intercept_column_predict_is_ALLCLOSE_BUT_NOT_EQUAL(covariates):
    """The companion for §7.3's first detail, and it is the reason the oracle is an equality.

    `beta[0] + X @ beta[1:]` is the same mathematics as `Xc @ beta` over an intercept-prepended `Xc`
    and a different floating-point summation order. Measured on the workbook's shape: max|difference|
    1.11e-16 on 10 of 60 rows. A test written with `allclose` passes on this form; `array_equal`
    fails on it, which is what turns a stylistic choice into a tested one.
    """
    fit, X = golden_fit(covariates)
    scalar_intercept = 1.0 / (1.0 + np.exp(-np.clip(
        fit.beta[0] + X.to_numpy(dtype=float) @ fit.beta[1:],
        -config.FIRTH_ETA_CLIP, config.FIRTH_ETA_CLIP)))
    assert scalar_intercept == pytest.approx(fit.p)          # allclose: passes
    assert not np.array_equal(scalar_intercept, fit.p)       # array_equal: does not
    assert 0.0 < float(np.max(np.abs(scalar_intercept - fit.p))) < 1e-15


def test_WITHOUT_the_clip_predict_AGREES_on_the_golden_frame_and_DIVERGES_above_eta_500():
    """The companion for §7.3's second detail. Agreement on ordinary data is the point.

    `model._probabilities` clips at FIRTH_ETA_CLIP (model.py:320) and model.py:317-318 gives the
    reason: a reader of `Fit` applying "a different clip — or none — from the one inside the loop" is
    reading a quantity the fit was not computed from. On the golden frame no |eta| comes near 500, so
    an unclipped implementation is bit-for-bit identical and no test on that frame can see the defect.
    """
    fit, X = golden_fit(GOLDEN_COVARIATE_PAIRS[0])
    Xc = np.column_stack([np.ones(len(X)), X.to_numpy(dtype=float)])
    unclipped = 1.0 / (1.0 + np.exp(-(Xc @ fit.beta)))
    assert np.array_equal(unclipped, fit.p)                       # agrees where it cannot matter
    assert float(np.max(np.abs(Xc @ fit.beta))) < config.FIRTH_ETA_CLIP

    # and on a design that reaches the tail, it does not. One column, beta = (0, 1), so eta IS x.
    tail = model.Fit(beta=np.array([0.0, 1.0]), p=np.array([np.nan]), iterations=1,
                     converged_on="likelihood", columns=("x",), first_step_norm=0.0,
                     rescales=0, halvings=0)
    frame = pd.DataFrame({"x": [-800.0]})
    with np.errstate(over="ignore"):
        unclipped_tail = 1.0 / (1.0 + np.exp(-frame["x"].to_numpy(dtype=float)))
    assert unclipped_tail[0] == 0.0
    assert model.predict(tail, frame)[0] == 7.124576406741285e-218


def test_the_inverse_logit_form_is_safe_at_the_POSITIVE_tail_and_UNSAFE_at_the_NEGATIVE_one():
    """Stage 9 §3.3's corrected fact. An earlier draft asserted it was safe at both and tested one.

    Both computed alongside, so the reason `predict` clips is in the suite and not only in a fence.
    """
    with np.errstate(over="ignore"):
        for eta in (800.0, 1000.0):
            assert 1.0 / (1.0 + np.exp(-eta)) == 1.0
            assert np.isnan(np.exp(eta) / (1.0 + np.exp(eta)))
        assert 1.0 / (1.0 + np.exp(-500.0)) == np.exp(500.0) / (1.0 + np.exp(500.0)) == 1.0
        # the negative tail, where the "better-conditioned" form is the one that fails
        assert 1.0 / (1.0 + np.exp(800.0)) == 0.0
    assert model._probabilities(np.array([[1.0]]), np.array([-800.0]))[0] == 7.124576406741285e-218


def test_the_unclipped_negative_tail_actually_WARNS_and_the_clipped_one_does_not():
    """The overflow is a `RuntimeWarning`, which is the only thing that would ever have surfaced it —
    and it is emitted by a form that returns a plausible 0.0 rather than raising."""
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        1.0 / (1.0 + np.exp(np.float64(800.0)))
    assert any(issubclass(w.category, RuntimeWarning) for w in caught)

    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        model._probabilities(np.array([[1.0]]), np.array([-800.0]))
    assert not [w for w in caught if issubclass(w.category, RuntimeWarning)]


@pytest.mark.parametrize("mutate", ["reorder", "drop", "rename"])
def test_predict_rejects_a_design_that_is_not_the_fits_by_NAME_and_by_ORDER(mutate):
    """A set comparison would pass on `reorder`, and a width check would pass on `rename`.

    `firth`'s beta is positional after the intercept, so the column names have to agree as a
    SEQUENCE. `design` drops constant columns and returns their names, which is how a frame built
    from a covariate list acquires a different width from the one a fit saw.
    """
    fit, X = golden_fit(GOLDEN_COVARIATE_PAIRS[0])
    broken = {"reorder": X[list(X.columns)[::-1]],
              "drop": X.drop(columns=[X.columns[-1]]),
              "rename": X.rename(columns={X.columns[-1]: "elsewhere"})}[mutate]
    with pytest.raises(config.SchemaError) as e:
        model.predict(fit, broken)
    assert "predict was given columns" in message_of(e)


def test_WITHOUT_the_column_check_a_REORDERED_design_returns_a_NUMBER_instead_of_an_error():
    """The companion is the point (Stage 9 §7.3). A mismatched positional dot product does not fail
    — it succeeds, and returns a full array of finite, plausible probabilities."""
    fit, X = golden_fit(GOLDEN_COVARIATE_PAIRS[0])
    reordered = X[list(X.columns)[::-1]]
    Xc = np.column_stack([np.ones(len(reordered)), reordered.to_numpy(dtype=float)])
    silent = 1.0 / (1.0 + np.exp(-np.clip(Xc @ fit.beta,
                                          -config.FIRTH_ETA_CLIP, config.FIRTH_ETA_CLIP)))
    assert np.all(np.isfinite(silent)) and np.all((silent > 0.0) & (silent < 1.0))
    assert not np.array_equal(silent, fit.p)


# --- 15.7 (model's half)  separation converges, and NO bound turns that into a failure -----------------
#
# Stage 8 §6 established that separation in `polr` converges silently and therefore needs
# POLR_MAX_ABS_BETA to become a countable failure. Both halves of that finding were re-asked of
# `model.firth` for Stage 9. **The first reproduces exactly** and is asserted here. **The second does
# not**: `sich`'s max_abs_beta across replicates runs continuously from 0.8558 to 80.9825 with 42.7%
# above 8 and no gap anywhere, so there is no empty band to calibrate a bound from — and Stage 9 §9.6
# prescribes none, for a STRUCTURAL reason rather than a measured one. At Stage 8 `beta` IS the
# estimand, so a degenerate one is a degenerate answer; here `m_a(X)` is a nuisance and only its
# PREDICTIONS enter tau, and predictions are probabilities, bounded in [0, 1] whatever beta does.


def test_a_perfectly_separated_binary_frame_CONVERGES_through_firth():
    """Asserted POSITIVELY, exactly as §14.7 asserts its ordinal analogue, and for the same reason: a
    test written as "a separated frame raises" passes on a broken fitter and fails on the correct one.

    Firth's penalty does not shrink a large coefficient into a small one and it does not turn
    separation into a failure — it turns a fit with NO MAXIMUM into one that converges in six
    iterations with every safeguard counter at zero, and hides the separation completely. Nothing in
    the returned `Fit` is out of range.
    """
    frame = separated_binary_frame()
    assert list(frame["y"]) == list(frame[config.TREATMENT])          # perfectly separated
    X, dropped = model.design(frame, ("center", "atrial_fib"))
    X = X.copy()
    X.insert(0, config.TREATMENT, frame[config.TREATMENT].to_numpy(dtype=float))
    fit = model.firth(X, frame["y"].to_numpy(dtype=float))

    # only two of the four declared centres appear, so the constant-column rule fires TWICE
    assert dropped == ("center_Lugano", "center_USZ")
    assert fit.columns == (config.TREATMENT, "atrial_fib", "center_CHUV")
    assert (fit.iterations, fit.converged_on) == (6, "likelihood")
    assert (fit.rescales, fit.halvings) == (0, 0)
    beta_treatment = float(fit.beta[1 + fit.columns.index(config.TREATMENT)])
    assert beta_treatment == pytest.approx(6.0890571141, abs=1e-6)
    assert float(np.exp(beta_treatment)) == pytest.approx(441.005397, abs=1e-4)
    assert float(np.max(np.abs(fit.beta))) == pytest.approx(6.089057, abs=1e-6)
    assert float(fit.p.min()) > 0.0 and float(fit.p.max()) < 1.0        # strictly interior
    assert float(fit.p.min()) == pytest.approx(4.545e-02, abs=1e-5)
    assert float(fit.p.max()) == pytest.approx(0.954546, abs=1e-6)


def test_the_UNPENALISED_mle_on_the_same_frame_FAILS_OUTRIGHT():
    """Computed alongside, so the penalty's effect is a measured contrast rather than a claim.

    The unpenalised fit does not merely grow — it has no maximum, and all three routes below say so
    in a different way. That contrast is stronger than "the coefficient gets large", which is what a
    reader takes from Stage 8's ordinal analogue.

    **THE ITERATION CAP IS PART OF THE MEASUREMENT AND `maxiter` IS THEREFORE NOT A DETAIL.** Stage 9
    §9.6 records the Newton route as raising `LinAlgError: Singular matrix`; measured here, it does so
    only once it is allowed to reach the singularity. At statsmodels' own default of 35 it stops
    without raising, at max_abs_beta 67.87 and `converged=False`. Both facts are asserted, because the
    second is the one that shows what an unpenalised fit looks like from the outside: a number.
    """
    statsmodels = pytest.importorskip("statsmodels.api")
    frame = separated_binary_frame()
    X, _ = model.design(frame, ("center", "atrial_fib"))
    X = X.copy()
    X.insert(0, config.TREATMENT, frame[config.TREATMENT].to_numpy(dtype=float))
    Xc = np.column_stack([np.ones(len(X)), X.to_numpy(dtype=float)])
    y = frame["y"].to_numpy(dtype=float)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        # allowed to run to the optimum, the information matrix goes singular and it raises
        with pytest.raises(np.linalg.LinAlgError) as e:
            statsmodels.Logit(y, Xc).fit(method="newton", maxiter=200, disp=0)
        assert "Singular matrix" in str(e.value)

        # stopped early, it returns a plausible number instead — and a different one at each cap,
        # which IS the statement that there is no maximum to return
        capped = statsmodels.Logit(y, Xc).fit(method="newton", disp=0)
        bfgs = statsmodels.Logit(y, Xc).fit(method="bfgs", maxiter=60, disp=0)

    assert capped.mle_retvals["converged"] is False
    assert float(np.max(np.abs(capped.params))) > 60.0

    # BFGS at 60 iterations reports SUCCESS, at a coefficient that is a function of the cap and not
    # of the data. Firth reaches 6.089057 on the same frame, in six iterations, and is right to.
    assert float(np.max(np.abs(bfgs.params))) == pytest.approx(33.0927, abs=1e-3)
    assert bfgs.mle_retvals["converged"] is True
    assert float(np.max(np.abs(bfgs.params))) > 5 * 6.089057


def test_NO_BOUND_ON_THE_FIRTH_COEFFICIENT_EXISTS_and_that_is_a_decision():
    """Stage 9 §9.6, asserted as a deliberate absence so that adding one is a visible change.

    Stage 8 declares POLR_MAX_ABS_BETA and `outcome._assert_reportable` reads it. `firth` has no
    counterpart: it returned the separated fit above, with max_abs_beta 6.09, and nothing anywhere
    rejected it. The two reasons are §9.6's — there is no empty band to calibrate a bound from, since
    `sich`'s replicate max_abs_beta runs continuously from 0.8558 to 80.9825; and at Stage 9 the
    coefficient is a NUISANCE whose predictions are probabilities, so every term of the augmented
    estimate is bounded whatever beta does. `outcome.py`'s half of this is §15.7's.
    """
    assert not hasattr(config, "FIRTH_MAX_ABS_BETA")
    assert "MAX_ABS_BETA" not in SOURCE          # model.py reads no such bound on any path
    assert hasattr(config, "POLR_MAX_ABS_BETA")  # and Stage 8's exists, so this is a contrast


# --- Stage 12 §20.5  ordinal_probabilities --------------------------------------------------------
#
# `P(Y = j | x)` for every `j`, which is what g-computation averages and which nothing in this
# pipeline could do before: `predict` evaluates a Firth binary `Fit`, and `_ord_pieces` computes the
# two cumulative probabilities bracketing ONE category per row — the row's own observed one.


def _hand_polr_fit() -> model.PolrFit:
    """A `PolrFit` with KNOWN alpha and beta, so §20.5's assertions are about the arithmetic.

    Constructed and not fitted, which is Stage 12 §20.5's own requirement: a property asserted only
    against a fitted object is a property asserted against this cohort.
    """
    return model.PolrFit(
        beta=np.array([0.5, -0.25]), alpha=np.array([-1.0, 0.0, 1.5]),
        categories=(0, 1, 2, 3), columns=("a", "b"), iterations=3,
        converged_on="likelihood", first_step_norm=1.0, rescales=0, halvings=0)


def _hand_design() -> pd.DataFrame:
    return pd.DataFrame({"a": [1.0, 0.0, 2.0, -1.5], "b": [0.0, 1.0, -1.0, 0.5]})


def test_ordinal_probabilities_agrees_with_EXPIT_DIFFERENCES_COMPUTED_BY_HAND():
    """Stage 12 §20.5. Constructed, not measured, so the assertion is about the arithmetic."""
    fit, X = _hand_polr_fit(), _hand_design()
    got = model.ordinal_probabilities(fit, X)

    eta = X["a"].to_numpy() * fit.beta[0] + X["b"].to_numpy() * fit.beta[1]
    def expit(z):
        return 1.0 / (1.0 + np.exp(-z))
    want = np.column_stack([
        expit(fit.alpha[0] + eta),
        expit(fit.alpha[1] + eta) - expit(fit.alpha[0] + eta),
        expit(fit.alpha[2] + eta) - expit(fit.alpha[1] + eta),
        1.0 - expit(fit.alpha[2] + eta),
    ])
    assert got.shape == (len(X), len(fit.categories))
    assert np.allclose(got, want, atol=0.0, rtol=1e-15)


def test_the_rows_sum_to_EXACTLY_ONE_and_the_STRUCTURAL_BOUNDS_ARE_EXACT():
    """Stage 12 §20.5, §6.1. The two bounds are literal 0.0 and 1.0 columns, never an expit.

    `expit(POLR_ETA_CLIP)` is 1.0 in float64, but that is an accident of the clip. A row's
    probabilities must sum to exactly 1.0 BY CONSTRUCTION for the sum assertion to be about the
    arithmetic rather than about the clip — so this test fails if a future edit replaces either
    bracketing column with an expit of a large number.
    """
    fit, X = _hand_polr_fit(), _hand_design()
    got = model.ordinal_probabilities(fit, X)

    assert np.all(got.sum(axis=1) == 1.0)                    # EXACTLY, not to a tolerance
    assert np.all((got >= 0.0) & (got <= 1.0))

    eta = X.to_numpy() @ fit.beta
    first = 1.0 / (1.0 + np.exp(-(fit.alpha[0] + eta)))
    last = 1.0 - 1.0 / (1.0 + np.exp(-(fit.alpha[-1] + eta)))
    assert np.array_equal(got[:, 0], first)
    assert np.array_equal(got[:, -1], last)


def test_ordinal_probabilities_agrees_with_ord_pieces_on_the_OBSERVED_category():
    """Stage 12 §6.1. The new public name and the private one used inside the likelihood agree.

    `_ord_pieces` returns the two cumulative probabilities bracketing each observation's own
    category, and their difference IS that observation's `P(Y = y_i)`. So the column of the new
    function selected at each row's observed category must equal it — which is what makes the two
    one definition of the same probability rather than two.
    """
    fit, X = _hand_polr_fit(), _hand_design()
    got = model.ordinal_probabilities(fit, X)

    y_idx = np.array([0, 1, 3, 2])
    K = len(fit.alpha)
    upper, lower = y_idx, y_idx - 1
    gu, gl = model._ord_pieces(
        X.to_numpy(dtype=float) @ fit.beta, fit.alpha, upper, lower,
        y_idx <= K - 1, y_idx >= 1)[:2]
    assert np.allclose(got[np.arange(len(X)), y_idx], gu - gl, rtol=1e-14, atol=0.0)


def test_ordinal_probabilities_returns_the_FITTED_level_set_and_not_the_declared_one():
    """Stage 12 §6.2. Six columns where `C.MRS_LEVELS` declares seven, on a collapsed frame.

    `polr` collapses the response to the categories carrying positive weight before it fits anything.
    Re-expressing the result on a declared level set is the CALLER's, because a level for which the
    fit provides no cutpoint has no probability this function could invent (§6.3).
    """
    frame = collapsed_level_frame()
    X, _ = model.design(frame, (config.TREATMENT,) + config.STANDARDISATION_COVARIATES)
    fit = model.polr(X, frame[config.PRIMARY_OUTCOME].to_numpy(dtype=float))
    got = model.ordinal_probabilities(fit, X)

    assert len(fit.categories) == len(config.MRS_LEVELS) - 1
    assert got.shape[1] == len(fit.categories)
    assert got.shape[1] < len(config.MRS_LEVELS)
    assert np.allclose(got.sum(axis=1), 1.0, atol=1e-15)


def test_ordinal_probabilities_PREPENDS_NO_INTERCEPT_which_is_where_it_differs_from_predict():
    """Stage 12 §6.1's fourth rule. THE CUTPOINTS ARE THE INTERCEPTS.

    `predict` prepends its own because `firth` fits one; this must not, because O6 has already
    established that `X` carries none. The assertion is that the result is invariant to nothing being
    added — i.e. that a design widened by a column of ones would be a T5 rather than a silent shift.
    """
    fit, X = _hand_polr_fit(), _hand_design()
    with_intercept = X.copy()
    with_intercept.insert(0, "intercept", 1.0)
    with pytest.raises(config.SchemaError) as excinfo:
        model.ordinal_probabilities(fit, with_intercept)
    assert str(excinfo.value).startswith("T5")

    source = Path(model.__file__).resolve().read_text(encoding="utf-8")
    node = next(n for n in ast.parse(source).body
                if isinstance(n, ast.FunctionDef) and n.name == "ordinal_probabilities")
    body = ast.unparse(node)
    assert "np.column_stack([np.ones(len(X))" not in body
    assert "np.zeros(len(X))" in body and "np.ones(len(X))" in body   # the two EXACT bounds


# --- Stage 12 §20.12  polr_ri ---------------------------------------------------------------------
#
# The random-centre-intercept proportional-odds fit: adaptive Gauss-Hermite quadrature, an analytic
# gradient, BFGS with a relative trust region and per-iteration step-halving, and a floor on sigma.
#
# Two of its choices are MEASUREMENTS rather than conventions and both are asserted here, because
# Stage 12 §12.2 and §12.4 rest on them: adaptive quadrature converges where non-adaptive quadrature
# does not, and the analytic gradient is what makes the arm affordable at all.


def _ri_design(frame: pd.DataFrame):
    return frame[["x1", "x2"]].astype(float), frame["y"].to_numpy(dtype=float), frame["group"]


def test_the_ANALYTIC_GRADIENT_agrees_with_central_differences_at_every_sigma_the_bootstrap_reaches():
    """Stage 12 §20.12, §12.4 departure 2. The frozen-node approximation, pinned.

    The quadrature nodes are recomputed at the current parameters and then held FIXED for that
    evaluation's gradient, so the gradient is exact for the frozen-node sum and neglects the
    derivative of the node positions. **Measured relative error 4.3e-8 to 6.9e-8 across
    sigma in {0.2, 0.54, 1.5, 3.0}** — two orders below `POLR_SCORE_TOL = 1e-6` — so the
    approximation is named rather than hidden, and a change to the mode-finding step that made it
    worse fails HERE rather than shifting an estimate.

    A numerical gradient is `2 * n_par` objective evaluations against the analytic one's one, which
    Stage 12 §12.4 measured as a 5x difference in the whole arm's cost. It is part of the
    specification, not an optimisation.
    """
    frame = two_group_ri_frame(n_groups=6, per_group=60)
    X, y, groups = _ri_design(frame)
    Xn = X.to_numpy(dtype=float)
    categories = tuple(int(v) for v in np.unique(y))
    y_idx = np.searchsorted(np.asarray(categories, dtype=float), y)
    labels = np.asarray([str(v) for v in groups])
    unique = tuple(sorted(set(labels.tolist())))
    gidx = np.searchsorted(np.asarray(unique), labels)
    nodes, weights = np.polynomial.hermite.hermgauss(config.POLR_RI_NODES)
    start = model.polr(X, y)
    K, m = len(categories) - 1, Xn.shape[1]

    for sigma in (0.20, 0.54, 1.50, 3.00):
        par = np.concatenate([start.alpha, start.beta, [np.log(sigma)]])
        _, gradient, _, _ = model._ri_objective(
            par, Xn, y_idx, K, m, gidx, len(unique), nodes, np.log(weights), True)
        numerical = np.zeros_like(gradient)
        for j in range(len(par)):
            step = 1e-6 * max(1.0, abs(par[j]))
            up, down = par.copy(), par.copy()
            up[j] += step
            down[j] -= step
            numerical[j] = (
                model._ri_objective(up, Xn, y_idx, K, m, gidx, len(unique), nodes,
                                    np.log(weights), True)[0]
                - model._ri_objective(down, Xn, y_idx, K, m, gidx, len(unique), nodes,
                                      np.log(weights), True)[0]) / (2.0 * step)
        relative = float(np.max(np.abs(gradient - numerical)
                                / np.maximum(np.abs(numerical), 1.0)))
        assert relative < 1e-6, (sigma, relative)


def test_ADAPTIVE_quadrature_is_STABLE_in_the_node_count_and_NON_ADAPTIVE_IS_NOT():
    """Stage 12 §20.12, §12.2, §12.3. The measurement the whole arm rests on.

    **Non-adaptive Gauss-Hermite is NON-MONOTONE in the node count** — measured here, sigma_hat runs
    0.6317 / 0.7373 / 0.5297 / 0.5668 / 0.5986 at 5 / 7 / 9 / 11 / 15 nodes against an adaptive
    reference of 0.5710. A node count chosen by "increase it until the answer stops moving" would
    have stopped at the wrong place, and Stage 12 records that its first prototype did.

    Adaptive quadrature is stable to 3.5e-11 from NINE nodes, which is why `POLR_RI_NODES = 11` sits
    one step above the plateau: a workbook whose posterior is slightly less Gaussian needs no
    re-tuning.

    Asserted so that a change to non-adaptive quadrature FAILS rather than shifting an answer.
    """
    frame = two_group_ri_frame(n_groups=6, per_group=60)
    X, y, groups = _ri_design(frame)

    adaptive = {n: model._polr_ri(X, y, groups, n, adaptive=True).sigma
                for n in (9, 11, 15, 31)}
    reference = adaptive[31]
    assert max(abs(value - reference) for value in adaptive.values()) < 1e-7

    non_adaptive = {n: model._polr_ri(X, y, groups, n, adaptive=False).sigma
                    for n in (5, 7, 9)}
    assert min(abs(value - reference) for value in non_adaptive.values()) > 1e-3

    # And the non-monotonicity itself, which is the half a "just use more nodes" reading misses.
    ordered = [non_adaptive[n] for n in (5, 7, 9)]
    assert not (ordered == sorted(ordered) or ordered == sorted(ordered, reverse=True))


def test_polr_ri_records_the_NODE_COUNT_because_it_CHANGES_THE_ANSWER():
    """Stage 12 §3.1, §12.3. `RIFit` carries `nodes` and `PolrFit` has no counterpart.

    `POLR_TOL` changes only how precisely the same answer is found; the node count changes the answer.
    A record that omitted it would let two fits with different numerical content compare equal on
    every field.
    """
    frame = two_group_ri_frame(n_groups=6, per_group=60)
    X, y, groups = _ri_design(frame)
    fit = model.polr_ri(X, y, groups)
    assert fit.nodes == config.POLR_RI_NODES
    assert {f.name for f in dataclasses.fields(model.RIFit)} - {
        f.name for f in dataclasses.fields(model.PolrFit)} >= {"nodes", "sigma", "b", "b_sd",
                                                               "groups", "at_floor"}
    assert not hasattr(model.polr(X, y), "nodes")
    # ODD and at least 9, which `test_config.py` also asserts: a symmetric rule with an odd node
    # count places a node AT the conditional mode, which is where the mass is.
    assert config.POLR_RI_NODES % 2 == 1 and config.POLR_RI_NODES >= 9


def test_polr_ri_REACHES_THE_FLOOR_on_a_flat_frame_and_DOES_NOT_RAISE():
    """Stage 12 §20.12, §12.5. sigma^2_C = 0 is a legitimate answer that [§14a] names.

    Dropping such replicates would select the bootstrap on the value of the very parameter the arm
    exists to examine. `sigma` is the floor EXACTLY, because `exp(log(1e-4))` is
    0.00010000000000000009 and a fit reporting that would make `at_floor` True while
    `sigma == POLR_RI_SIGMA_FLOOR` was False — two fields disagreeing about one fact.
    """
    X, y, groups = _ri_design(flat_ri_frame())
    fit = model.polr_ri(X, y, groups)                 # does not raise
    assert fit.at_floor is True
    assert fit.sigma == config.POLR_RI_SIGMA_FLOOR
    assert np.all(np.isfinite(fit.b)) and np.all(np.isfinite(fit.b_sd))


def test_at_the_floor_polr_ri_COLLAPSES_TO_model_polr_which_is_14a_s_own_sentence():
    """Stage 12 §20.12. *"if outcomes truly do not differ by centre given X … it collapses to the
    pooled model"* — as an assertion.

    **The tolerance is 1e-5 and not the specification's 1e-6, and the difference is measured rather
    than assumed.** At the floor `sigma` is 1e-4 and not 0, so the marginal model is not exactly the
    pooled one; and both fits stop on the `POLR_TOL` likelihood criterion rather than at an exact
    stationary point. Measured on `flat_ri_frame`: 4.9e-6 on the cutpoints and 1.3e-6 on the
    coefficients. 1e-5 is the smallest round bound above what the two estimators can agree to, and
    pinning 1e-6 would be pinning a number this construction does not reach.
    """
    X, y, groups = _ri_design(flat_ri_frame())
    hierarchical = model.polr_ri(X, y, groups)
    pooled = model.polr(X, y)

    assert hierarchical.at_floor
    assert hierarchical.categories == pooled.categories
    assert hierarchical.columns == pooled.columns
    assert float(np.max(np.abs(hierarchical.alpha - pooled.alpha))) < 1e-5
    assert float(np.max(np.abs(hierarchical.beta - pooled.beta))) < 1e-5


def test_polr_ri_RECOVERS_a_known_sigma_when_there_are_enough_groups_to_estimate_one():
    """Stage 12 §20.12. The one test here that checks the estimator estimates the RIGHT THING.

    **Recovering `sigma` is a different question from fitting it, and it is not answerable at two
    clusters.** Measured at `RI_SIGMA_TRUE = 0.8` with 60 records per group: sigma_hat 0.681 at
    twelve groups and 1.168 at two. The tolerance below is stated against the twelve-group case for
    that reason, and the two-group case is asserted only to FIT — which is what [§14a]'s *"a centre
    variance from four clusters is fragile"* means as a number.
    """
    X, y, groups = _ri_design(two_group_ri_frame(n_groups=12, per_group=60))
    fit = model.polr_ri(X, y, groups)
    assert abs(fit.sigma - RI_SIGMA_TRUE) < 0.25
    assert not fit.at_floor
    assert abs(fit.beta[0] - RI_BETA_TRUE) < 0.25
    assert len(fit.groups) == len(fit.b) == len(fit.b_sd) == 12
    assert fit.groups == tuple(sorted(fit.groups))            # deterministic keying

    X2, y2, groups2 = _ri_design(two_group_ri_frame(n_groups=2, per_group=60))
    fitted_at_two = model.polr_ri(X2, y2, groups2)             # fits, and that is all that is claimed
    assert fitted_at_two.sigma > config.POLR_RI_SIGMA_FLOOR


def test_the_posterior_SD_is_LARGER_where_a_group_carries_FEWER_records():
    """Stage 12 §12.6. The arm's honest limitation as a property rather than a caveat.

    On the workbook the never-IVT centre's intercept is the least precisely estimated of the four and
    it is the one the whole arm rests on, because it is the centre whose IVT counterfactual is
    entirely borrowed. Here the same relationship is asserted on a construction: a group with a
    quarter of the records has a visibly wider posterior.
    """
    frame = pd.concat([
        two_group_ri_frame(n_groups=1, per_group=120, sigma=0.0, seed=1).assign(group="BIG"),
        two_group_ri_frame(n_groups=1, per_group=30, sigma=0.0, seed=2).assign(group="SMALL"),
    ], ignore_index=True)
    frame["group"] = frame["group"].astype("string")
    X, y, groups = _ri_design(frame)
    fit = model.polr_ri(X, y, groups)
    sd = dict(zip(fit.groups, fit.b_sd))
    assert sd["SMALL"] > sd["BIG"]


def test_T4_fires_on_FEWER_THAN_TWO_GROUPS_and_its_token_is_T4_and_NOT_polr_ri():
    """Stage 12 §20.12, §14. **The one bucket assertion the raise-site scan CANNOT make.**

    `bootstrap.bucket` classifies by the FIRST TOKEN of the message and raises on a token it does not
    know. `polr_ri:` IS a valid key — T2 and T3 use it — so a T4 message leading with `polr_ri:`
    would be counted as `nonconvergence` **while `"T4" -> "degenerate_design"` sat in the map
    unreachable**, and the scan would pass because the token it found was in the map. That is the
    second documented instance of a wrong-but-mapped bucket, and this is the assertion for it.

    With one group, `b_c` is exactly collinear with the cutpoints and `sigma` is not identified: the
    likelihood is flat in it, so the fit returns whatever the start value was — finite and plausible.
    """
    frame = two_group_ri_frame(n_groups=6, per_group=40)
    X, y, _ = _ri_design(frame)
    single = pd.Series(["ONE"] * len(y), dtype="string")
    with pytest.raises(model.FitError) as excinfo:
        model.polr_ri(X, y, single)
    message = str(excinfo.value)
    assert message.startswith("T4")
    assert not message.startswith("polr_ri:")
    assert message.split()[0] == "T4"
    assert config.FAILURE_BUCKETS[message.split()[0]] == "degenerate_design"


def test_T2_and_T3_lead_with_polr_ri_and_share_the_nonconvergence_bucket(monkeypatch):
    """Stage 12 §14, §20.12. `"polr:"`'s arrangement, unchanged.

    T2 is step-halving exhausted and T3 is no convergence in `POLR_RI_MAX_ITER`; both are convergence
    failures, so they share a bucket by sharing a leading token.
    """
    frame = two_group_ri_frame(n_groups=6, per_group=40)
    X, y, groups = _ri_design(frame)

    # T3: one iteration is not enough for any frame.
    monkeypatch.setattr(config, "POLR_RI_MAX_ITER", 1)
    with pytest.raises(model.FitError) as excinfo:
        model.polr_ri(X, y, groups)
    assert str(excinfo.value).startswith("polr_ri:")
    assert "no convergence" in str(excinfo.value)
    assert config.FAILURE_BUCKETS[str(excinfo.value).split()[0]] == "nonconvergence"

    # T2: every trial step is made worse than the iterate, so no halving can accept one.
    #
    # **`POLR_MAX_HALVINGS = 0` is NOT how to reach it, and what happens instead is worth recording**:
    # `polr_ri` calls `model.polr` for its start values (§12.8), so that constant makes the POOLED fit
    # raise `polr:` first and the message never reaches this fitter. That is the start-value
    # dependency §13.2 rests on, visible from the other side — there is no hierarchical arm without a
    # pooled fit — so the assertion below is made on the objective rather than on a shared constant.
    monkeypatch.setattr(config, "POLR_RI_MAX_ITER", 200)
    monkeypatch.setattr(config, "POLR_MAX_HALVINGS", 0)
    with pytest.raises(model.FitError) as excinfo:
        model.polr_ri(X, y, groups)
    assert str(excinfo.value).startswith("polr:")            # the POOLED fit, not this one

    monkeypatch.setattr(config, "POLR_MAX_HALVINGS", 30)
    honest = model._ri_objective
    seen = {"calls": 0}

    def never_uphill(*args, **kwargs):
        value, gradient, modes, tau = honest(*args, **kwargs)
        seen["calls"] += 1
        return (value if seen["calls"] == 1 else -np.inf), gradient, modes, tau

    monkeypatch.setattr(model, "_ri_objective", never_uphill)
    with pytest.raises(model.FitError) as excinfo:
        model.polr_ri(X, y, groups)
    assert str(excinfo.value).startswith("polr_ri:")
    assert "step-halving" in str(excinfo.value)
    assert config.FAILURE_BUCKETS[str(excinfo.value).split()[0]] == "nonconvergence"
    assert seen["calls"] > config.POLR_MAX_HALVINGS          # every halving was tried


def test_polr_ri_inherits_O1_to_O6_and_raises_on_a_MISALIGNED_grouping_vector():
    """Stage 12 §12. The design preconditions are `_assert_polr_fittable`'s, unweighted.

    `groups` is aligned by POSITION with `y`, as `smd`'s arguments are, so a length mismatch is a
    grouping nobody specified — and it would return a between-group variance for it.
    """
    frame = two_group_ri_frame(n_groups=4, per_group=40)
    X, y, groups = _ri_design(frame)

    broken = X.copy()
    broken.iloc[0, 0] = np.nan
    with pytest.raises(model.FitError) as excinfo:
        model.polr_ri(broken, y, groups)
    assert str(excinfo.value).startswith("O1")

    with pytest.raises(config.SchemaError) as excinfo:
        model.polr_ri(X, y, groups.iloc[:-5])
    assert str(excinfo.value).startswith("T4")
    assert "aligned by POSITION" in str(excinfo.value)


def test_polr_ri_has_NO_STANDARD_ERROR_FIELD_for_PolrFits_reason_and_one_more():
    """Stage 12 §3.1. A correctness claim rather than an omission.

    At four clusters a variance parameter's asymptotic standard error is the least trustworthy number
    in the output, and [§10]'s percentile bootstrap is the prespecified interval. `b_sd` is NOT a
    standard error for `sigma` — it is the curvature of each group's conditional posterior at its
    mode, which is what the adaptive quadrature needs to place its nodes.
    """
    fields = {f.name for f in dataclasses.fields(model.RIFit)}
    assert "se" not in fields and "std_err" not in fields and "sigma_se" not in fields
    assert "b_sd" in fields
    assert not any("se" == f or f.endswith("_se") for f in fields)


def test_the_curvature_is_held_to_the_SIGN_THE_MATHEMATICS_GUARANTEES():
    """Stage 12's own finding, recorded as a test because it cost two replicates in 2000.

    The per-group log-posterior is a sum of concave ordinal terms plus a Gaussian log-prior, so
    `d2/db2` is negative EXACTLY and the total cannot exceed `-1/sigma^2`. In float64 it sometimes
    does: `Huu = ddu/p - (du/p)**2` is a difference of two large quantities in the tails and
    cancellation can return a small positive value. When the sum crossed zero, `sqrt` returned `nan`,
    the quadrature nodes were `nan`, and the line search rejected the step — so the fit STILL
    RETURNED THE RIGHT ANSWER through a halving loop doing the wrong job, and in 2 of 2000 replicates
    exhausted the halvings and raised T2 for a floating-point artefact.

    The clamp is on the ordinal contribution's own guaranteed sign, so it enforces a property rather
    than choosing a magnitude, and the bound it produces is TIGHT: attained whenever the data
    contribute no curvature, which is the far-tail case that produces the noise.
    """
    positive = np.array([1e-13, 0.0, -3.0])
    clamped = model._curvature(positive, 1.0 / 0.5 ** 2)
    assert np.all(clamped <= -(1.0 / 0.5 ** 2))
    assert np.all(np.isfinite(np.sqrt(-1.0 / clamped)))
    # And it is a no-op wherever the computed curvature already has the right sign.
    honest = np.array([-1.0, -20.0])
    assert np.array_equal(model._curvature(honest, 4.0), honest - 4.0)


def test_the_far_node_probability_FLOOR_keeps_the_derivatives_finite():
    """Stage 12's second numerical finding. A far node is NEGLIGIBLE, not invalid.

    At a quadrature node the posterior has moved away from, `gu` and `gl` are both 1.0 in float64 and
    `p = gu - gl` is EXACTLY ZERO by cancellation — so `log p` is -inf and `du/p`, `(du/p)**2` and
    `du*dl/p**2` are all `nan`. `polr`'s convention is the opposite — return -inf and let step-halving
    reject the step — and it is right THERE, where a non-positive probability means crossed cutpoints
    at the iterate. Here it means a far node, which is not a property of the iterate at all.

    The exponent is chosen against `p**2` and not against `p`: `Hul = du*dl/p**2`, so a floor below
    about 1.5e-154 makes `p**2` underflow and reintroduces the 0/0 this fixes.
    """
    assert model._RI_P_FLOOR ** 2 > 0.0                      # the reason for the exponent
    assert model._RI_P_FLOOR ** 2 > np.finfo(float).tiny
    alpha = np.array([-2.0, 0.0, 2.0])
    far = np.array([[900.0], [-900.0]])                      # beyond any expit's resolution
    p, A, B, d2 = model._ord_b_derivatives(far, alpha, np.array([[0], [3]]), len(alpha))
    for array in (p, A, B, d2):
        assert np.all(np.isfinite(array))
    assert np.all(p > 0.0)


def test_polr_ri_runs_CLEAN_with_numpy_warnings_as_errors():
    """Stage 12. The two numerical fixes above, asserted together at the level that matters.

    A fit that produces the right answer while emitting overflow and invalid-value warnings is a fit
    whose line search is absorbing arithmetic it should never have seen. Stage 12 §12.5 makes exactly
    that point about the sigma -> 0 boundary — *"the fit still returned … but 'still returned' is not
    a specification"* — and this is that sentence as a test.
    """
    X, y, groups = _ri_design(two_group_ri_frame(n_groups=6, per_group=60))
    with warnings.catch_warnings():
        warnings.simplefilter("error", RuntimeWarning)
        fit = model.polr_ri(X, y, groups)
    assert np.isfinite(fit.sigma)

    X, y, groups = _ri_design(flat_ri_frame())
    with warnings.catch_warnings():
        warnings.simplefilter("error", RuntimeWarning)
        floored = model.polr_ri(X, y, groups)
    assert floored.at_floor


def test_model_py_IMPORTS_NO_OPTIMISER_and_scipy_stays_test_only():
    """Stage 12 §2, §23. `polr_ri` is a 40-line BFGS loop and not a `scipy.optimize` call.

    The AST scan Stage 6 §12.7 established forbids an optimiser in this module: a shipped module
    importing one is a shipped module one edit from a second estimator. `polr_ri` needs Gauss-Hermite
    nodes, and `numpy.polynomial.hermite.hermgauss` supplies them — so no dependency is added, and
    Stage 12 §12.4 records what a `scipy` dependency would have bought and why the reproducibility
    surface is not worth it.
    """
    source = Path(model.__file__).resolve().read_text(encoding="utf-8")
    tree = ast.parse(source)
    imported = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            imported.update(a.name for a in node.names)
        elif isinstance(node, ast.ImportFrom) and node.module:
            imported.add(node.module)
    assert not [name for name in imported if name.split(".")[0] == "scipy"]
    assert "hermgauss" in source
    assert "minimize" not in source
