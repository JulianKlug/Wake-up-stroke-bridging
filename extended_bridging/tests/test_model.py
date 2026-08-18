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
import os
from pathlib import Path

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
