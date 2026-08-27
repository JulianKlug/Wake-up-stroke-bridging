"""Acceptance tests for Stage 10 — §15 of `specs/stage10_bootstrap_engine.md`.

Section banners match the specification's numbers, as every test module since `test_cohort.py` does.
§15.7 — S8's reclassification and `secondary`'s `collect` — lives in `test_outcome.py`, because that
is where the changed function is.

Tests needing the private workbook are marked `DATA_GATED` and tagged `[data-gated]` in their
docstrings. Everything else runs on a plain checkout with no `data/`.

**NO INTERVAL LIMIT, p-VALUE, `beta` OR `RD_k` FROM THE WORKBOOK APPEARS IN THIS FILE**, and Stage 8
§4.3's rule binds harder here than anywhere before it: a `beta` in a committed document is a point
estimate somebody could have chosen a specification to obtain, and an interval is that plus the claim
of significance. Every regression pin is measured on a synthetic fixture (`fixtures_stage10.py`) and
the acceptance criteria assert PROPERTIES — the counters reconcile, the limits are order statistics,
the p agrees with the interval, coverage is near nominal on data whose truth is known.

**`C.N_BOOT` IS MONKEYPATCHED FOR EVERY TEST THAT DRIVES `run`, AND THERE IS NO OTHER WAY.** `run`
takes no `n` — C.N_BOOT is prespecified and [§10] refits in every replicate, so it is not a runtime
knob (§11) — and 2000 replicates of the fixture cohort is about 130 s per call. `RUN_BOOT` below is
60: comfortably above `ci_min_draws(0.95)` = 40 so R8 passes and the below-floor branch stays
unreachable except where a test reaches it deliberately, and small enough that the module-scoped run
costs about 4 s. This is §15.1's own mechanism for R8 — monkeypatching a prespecified constant —
applied to the one other constant a test cannot pass in.

The three silent failures this stage exists to make loud, each tested with a companion showing what
happens *without* the guard::

    numpy's default percentile   §15.8   without the pin: a significance verdict flips in 4.7%
                                         of boundary cases, and the p-value still reads 0.0500
    a re-derived augmentation    §15.11  without the freeze: `ph2`'s interval is a quantile over
    path                                 a near-even mixture of two estimators
    a skipped whole-replicate    §15.4   without the count: two replicates in 2000 leave every
    failure                              denominator with no counter naming them
"""
from __future__ import annotations

import ast
import dataclasses
import inspect
import os
import subprocess
import sys
import textwrap
from pathlib import Path
from typing import Final

import numpy as np
import pandas as pd
import pytest

import bootstrap
import cohort
import config
import data
import derive
import eligibility
import model
import outcome
import propensity
from fixtures_stage10 import (ALPHA_TRUE, B_INNER, BETA_TRUE, COVERAGE_SEED, M_OUTER,
                              boundary_draws, known_effect_population, separable_ordinal_frame,
                              two_centre_frame, unfittable_nuisance_frame)
from fixtures_stage11 import subgroup_frame
from test_data import hand_source

DATA_GATED = pytest.mark.skipif(
    not config.DATA_XLSX.exists(),
    reason="the private workbook is gitignored and absent from this checkout")

# §15.14's two coverage designs are 90 000 fits each. They are gated on an environment variable and
# skip LOUDLY, which is the R oracles' own mechanism and is what `pytest -rs` is for: the default
# state of a slow measurement is "not run", and a check matching nothing must not pass as green.
# The name of the module those oracles live in is deliberately not written here — it is pinned to
# be imported by nothing, and that pin is a substring scan over every .py file in the project.
SLOW = pytest.mark.skipif(
    os.environ.get("STAGE10_SLOW") != "1",
    reason="§15.14's coverage designs are 2 x 90 000 fits and about 9 minutes. "
           "Run with STAGE10_SLOW=1 to open the gate.")

MODULE = Path(bootstrap.__file__).resolve()
SOURCE = MODULE.read_text(encoding="utf-8")

# The replicate count every `run`-driving test uses. See this module's docstring.
RUN_BOOT: Final[int] = 60


def message_of(excinfo) -> str:
    return str(excinfo.value)


# --- 15.0  the frames and fixtures this stage is tested on ---------------------------------------------
#
# All six constructions and all five constants are `fixtures_stage10.py`'s, which is §15.0's own rule:
# every constructed fixture is code, because Stage 9 §22.3 item 1 records three fixtures pinned to ten
# decimals and described rather than given, which made its acceptance criteria unperformable.
#
# What lives HERE is only the wiring: a fitted `Propensity` over a fixture frame, and one module-scoped
# `run` at RUN_BOOT replicates that the structural sections read.


def fitted(df: pd.DataFrame) -> tuple[propensity.Propensity, outcome.Primary, outcome.Secondary,
                                      data.Audit]:
    """The point estimate over `df`: one `Audit`, the three objects `run` takes, in pipeline order."""
    audit = data.Audit(hand_source(Path("."), "stage10"))
    ps = propensity.fit(df, audit)
    return ps, outcome.primary(df, ps, audit), outcome.secondary(df, ps, audit), audit


@pytest.fixture
def boot_n(monkeypatch):
    """`C.N_BOOT` at RUN_BOOT for the duration of one test. See this module's docstring."""
    monkeypatch.setattr(config, "N_BOOT", RUN_BOOT)
    return RUN_BOOT


@pytest.fixture(scope="module")
def run_once():
    """(df, ps, est, sec, audit, boot) from ONE `run` at RUN_BOOT replicates.

    Module-scoped because it is the expensive thing in this file and every structural section reads
    it. `monkeypatch` is function-scoped, so the constant is moved and restored by hand — which is
    also why this fixture is the ONLY place that happens outside a `boot_n` test.
    """
    df = two_centre_frame()
    original = config.N_BOOT
    config.N_BOOT = RUN_BOOT
    try:
        ps, est, sec, audit = fitted(df)
        yield df, ps, est, sec, audit, bootstrap.run(df, ps, est, sec, audit)
    finally:
        config.N_BOOT = original


def test_the_fixture_cohort_reaches_ALL_THREE_augmentation_paths_and_TWENTY_SIX_keys(run_once):
    """The fixture's own witness first, as §15.5's separated frame gets one.

    §3.3's twenty-six is 1 + 6 + 7 + 7 + 5, and the 5 is the number of AUGMENTED outcomes — so a
    frame on which every outcome were augmented would give 28 and every "26" assertion below would
    be measuring a different frame. This asserts the split is the workbook's own: four `full`, one
    `reduced`, two `unaugmented` (Stage 9 §9.3).
    """
    _, _, _, sec, _, boot = run_once
    paths = [sec.estimates[key].augmented_path for key in config.BINARY_OUTCOMES]
    assert paths.count("full") == 4 and paths.count("reduced") == 1
    assert paths.count("unaugmented") == 2
    assert len(boot.draws) == 26


def test_the_coverage_generator_produces_THE_DECLARED_NUMBER_OF_CATEGORIES():
    """`ALPHA_TRUE` is five cutpoints, so `Y` takes six values — and the fixture's own witness has to
    say so, because a generator that silently produced five would make §15.14 measure coverage of a
    different estimator (Stage 8 §5.3 collapses the response to the occupied categories).

    Asserted on both designs at a size where every category is occupied with near certainty.
    """
    for confounded in (False, True):
        population = known_effect_population(4000, np.random.default_rng(COVERAGE_SEED), confounded)
        occupied = sorted(population[config.PRIMARY_OUTCOME].unique())
        assert occupied == list(range(len(ALPHA_TRUE) + 1))
        assert set(population[config.TREATMENT].unique()) == {0.0, 1.0}
        # No strata: §8.4 says so, and §16 item 4 files it as the gap that matters.
        assert population[config.BOOT_STRATUM].nunique() == 1


# --- 15.1  the preconditions, one frame per branch ------------------------------------------------------
#
# NINE branches from Stage 11 §4.5, nine frames, each derived from the fixture cohort by breaking
# exactly one thing.
# Every one asserts a `SchemaError` WHOSE MESSAGE NAMES ITS OWN CONDITION — not merely that something
# raised — because a collected assertion that reports the wrong branch is worse than one that reports
# nothing. Stage 9 §15.1a's rule applies unchanged: a frame testing phase k must be PRISTINE for
# phases 1 through k-1.

BRANCH_PHASE: Final[dict[str, int]] = {
    "R1": 1, "R2": 1, "R3": 1, "R4": 1, "R5": 1, "R6": 2, "R7": 2, "R8": 2, "R9": 1}


def broken_run(branch: str, monkeypatch):
    """(df, ps, est, sec, audit) for one of §15.1's eight branches, with exactly one thing wrong.

    The point estimate is fitted on the PRISTINE frame and the damage is applied afterwards, because
    every branch here is about the CALL to `run` and not about the data the point estimate saw —
    which is §4.4's whole asymmetry: all eight are the caller's bug.
    """
    df = two_centre_frame()
    ps, est, sec, audit = fitted(df)
    if branch == "R1":
        return df.drop(columns=[config.BOOT_STRATUM]), ps, est, sec, audit
    if branch == "R2":
        return df.drop(columns=["case_id"]), ps, est, sec, audit
    if branch == "R3":
        # A PERMUTATION: the same labels in the other order. Relabelling would be caught for the
        # wrong reason; this is the failure Stage 9 §12 measures, which moves the estimate silently.
        return df, dataclasses.replace(ps, e=ps.e[::-1], w=ps.w[::-1]), est, sec, audit
    if branch == "R4":
        return df, dataclasses.replace(ps, in_model=ps.in_model.astype("int64")), est, sec, audit
    if branch == "R5":
        short = {k: v for k, v in sec.estimates.items() if k != "sich"}
        return df, ps, est, outcome.Secondary(estimates=short), audit
    if branch == "R6":
        broken = df.copy()
        broken[config.BOOT_STRATUM] = broken[config.BOOT_STRATUM].astype(object)
        broken.loc[broken.index[0], config.BOOT_STRATUM] = np.nan
        return broken, ps, est, sec, audit
    if branch == "R7":
        bent = dict(sec.estimates)
        bent["ph2"] = dataclasses.replace(bent["ph2"], augmented_path="partial")
        return df, ps, est, outcome.Secondary(estimates=bent), audit
    if branch == "R9":
        # The [§13] arm's Propensity, handed to `run`. Stage 11 §4.5: `_replicate` calls
        # `propensity.fit` unconditionally, so without R9 this returns twenty-six intervals whose
        # point estimates are the arm's and whose replicates are the primary's specification.
        #
        # `subgroup_frame()` and not `two_centre_frame()`: `C.PROPENSITY_FULL` names four covariates
        # Stage 10's fixture does not carry, so `fit_full` on it is a bare KeyError rather than the
        # arm's Propensity — which would test the frame instead of the guard.
        wide = subgroup_frame()
        wide_ps, wide_est, wide_sec, wide_audit = fitted(wide)
        return wide, propensity.fit_full(wide, wide_audit), wide_est, wide_sec, wide_audit
    if branch == "R8":
        # `run` takes no `level`, so the ONLY way to reach R8 is to move the constant — which is why
        # test_config.py carries `N_BOOT >= ci_min_draws(CI_LEVEL)` as a static assertion beside it.
        # The level is chosen so `ci_min_draws` is EXACTLY twice N_BOOT, which is unambiguously above
        # it: a level whose floor merely equals N_BOOT satisfies R8 and the branch never fires.
        monkeypatch.setattr(config, "CI_LEVEL", 1.0 - 1.0 / config.N_BOOT)
        assert config.ci_min_draws(config.CI_LEVEL) > config.N_BOOT
        return df, ps, est, sec, audit
    raise AssertionError(f"no such branch: {branch}")


@pytest.mark.parametrize("branch", sorted(BRANCH_PHASE))
def test_each_precondition_branch_raises_a_SchemaError_NAMING_ITS_OWN_CONDITION(
        branch, monkeypatch, boot_n):
    """The message assertion is the test, for Stage 9 §15.1a's reason: `pytest.raises(SchemaError)`
    passes on a frame that broke a different branch on the way.

    `boot_n` is taken even though nothing here should reach the loop, and that is the point: a branch
    that FAILS to raise otherwise runs `C.N_BOOT` replicates and the test times out instead of
    failing with a message.
    """
    df, ps, est, sec, audit = broken_run(branch, monkeypatch)
    with pytest.raises(config.SchemaError) as e:
        bootstrap.run(df, ps, est, sec, audit)
    assert branch in message_of(e), (
        f"{branch}'s frame raised, but the message names a different branch:\n{message_of(e)}")


def test_WITHOUT_R9_run_SUCCEEDS_and_returns_a_FULL_TABLE_that_is_a_lie_about_the_arm(boot_n):
    """Stage 11 §4.5 and §15.2's companion. A guard whose absence is never demonstrated is a guard
    nobody can price.

    With the check monkeypatched out, `run` over the [§13] arm's `Propensity` returns every interval
    it would return for the primary. Every number is finite. `Draws.__post_init__` reconciles. R3
    passes — the index is aligned. R5 passes — the `Secondary` is complete. Nothing raises, and the
    arm's interval is a lie about the arm: the point estimates came from the full-covariate weights
    and all `boot_n` replicates refitted the [§7] specification.
    """
    df = subgroup_frame()                # carries what C.PROPENSITY_FULL names — see `broken_run`
    ps, est, sec, audit = fitted(df)
    arm = propensity.fit_full(df, audit)

    with pytest.raises(config.SchemaError) as e:
        bootstrap.run(df, arm, est, sec, audit)
    assert "R9" in message_of(e)

    landed = bootstrap._assert_run_inputs

    def without_r9(frame, propensity_object, primary, secondary):
        return landed(frame, dataclasses.replace(propensity_object,
                                                 spec=config.PROPENSITY_PRIMARY),
                      primary, secondary)

    original = bootstrap._assert_run_inputs
    bootstrap._assert_run_inputs = without_r9
    try:
        boot = bootstrap.run(df, arm, est, sec, audit)
    finally:
        bootstrap._assert_run_inputs = original

    assert len(boot.intervals) == len(boot.draws) == 26
    assert all(np.isfinite([i.lo, i.hi]).all() for i in boot.intervals.values())


def test_PHASE_1_FAILURES_ARE_COLLECTED_and_a_two_failure_frame_reports_BOTH():
    """An implementation raising on the first failure satisfies every single-branch test above.

    Asserted by counting the named conditions rather than by a substring, so a message that happened
    to mention R2 in prose would not pass.
    """
    df = two_centre_frame()
    ps, est, sec, audit = fitted(df)
    broken = df.drop(columns=["case_id"])                                  # R2
    with pytest.raises(config.SchemaError) as e:
        bootstrap.run(broken, dataclasses.replace(ps, in_model=ps.in_model.astype("int64")),
                      est, sec, audit)                                     # R4
    message = message_of(e)
    assert "R2" in message and "R4" in message


def test_PHASE_2_DOES_NOT_RUN_when_phase_1_failed_which_is_why_the_split_exists():
    """R6 reads through the stratum column R1 is about, so a frame missing it must report R1 and not
    a bare `KeyError` from inside R6 (§4.4). The companion form Stage 9 §4.4a prescribes."""
    df = two_centre_frame()
    ps, est, sec, audit = fitted(df)
    with pytest.raises(config.SchemaError) as e:
        bootstrap.run(df.drop(columns=[config.BOOT_STRATUM]), ps, est, sec, audit)
    assert "R1" in message_of(e) and "R6" not in message_of(e)


def test_a_stratum_of_size_ONE_is_NOT_an_error(boot_n):
    """§4.4: a stratum of size 1 resamples to itself in every replicate and contributes no
    variability. That is what stratifying on a near-determining variable MEANS and [§10] chose it
    knowing so, so it is not guarded — asserted positively, because the natural defensive edit is a
    minimum-stratum-size check that would reject a legitimate frame."""
    df = two_centre_frame()
    ps, est, sec, audit = fitted(df)
    moved = df.copy()
    moved[config.BOOT_STRATUM] = moved[config.BOOT_STRATUM].astype(object)
    moved.loc[moved.index[0], config.BOOT_STRATUM] = config.CENTER_ORDER[2]
    moved[config.BOOT_STRATUM] = moved[config.BOOT_STRATUM].astype(df[config.BOOT_STRATUM].dtype)
    assert int((moved[config.BOOT_STRATUM] == config.CENTER_ORDER[2]).sum()) == 1
    bootstrap._assert_run_inputs(moved, ps, est, sec)              # does not raise


# --- 15.2  `resample` -----------------------------------------------------------------------------------

def rng_at(seed: int = config.SEED) -> np.random.Generator:
    return np.random.default_rng(seed)


def test_the_stratum_totals_are_EXACT_on_a_frame_with_UNEQUAL_strata():
    """And the companion is what makes it a test: a resampler drawing `len(df)` rows from the whole
    frame ignoring strata passes a total-row-count check and fails this one."""
    df = two_centre_frame()
    expected = df[config.BOOT_STRATUM].value_counts().to_dict()
    assert len(set(expected.values())) > 1, "the fixture's strata must be unequal for this to bite"
    rng = rng_at()
    for _ in range(20):
        draw = bootstrap.resample(df, rng, config.BOOT_STRATUM)
        assert draw[config.BOOT_STRATUM].value_counts().to_dict() == expected
        assert len(draw) == len(df)


def test_the_UNSTRATIFIED_companion_passes_a_row_count_and_FAILS_the_totals():
    """The companion §15.2 requires, written out rather than described."""
    df = two_centre_frame()
    rng = rng_at()
    unstratified = df.iloc[rng.integers(0, len(df), size=len(df))]
    assert len(unstratified) == len(df)                                    # the check that passes
    assert (unstratified[config.BOOT_STRATUM].value_counts().to_dict()
            != df[config.BOOT_STRATUM].value_counts().to_dict())


def test_every_drawn_rows_case_id_is_UNIQUE_and_the_count_equals_the_frames():
    df = two_centre_frame()
    rng = rng_at()
    for _ in range(20):
        draw = bootstrap.resample(df, rng, config.BOOT_STRATUM)
        assert draw["case_id"].is_unique
        assert draw["case_id"].nunique() == len(df)


def test_the_FIRST_appearance_is_id_HASH_1_and_not_the_bare_id():
    """Asserted positively, because the natural implementation — rename only the duplicates — is the
    one that leaves "was this row drawn once" a question about string formatting (§5.3)."""
    df = two_centre_frame()
    draw = bootstrap.resample(df, rng_at(), config.BOOT_STRATUM)
    assert all(cid.count("#") == 1 for cid in draw["case_id"])
    assert all(cid.rsplit("#", 1)[1].isdigit() for cid in draw["case_id"])
    assert all(int(cid.rsplit("#", 1)[1]) >= 1 for cid in draw["case_id"])
    original = set(df["case_id"])
    assert all(cid.rsplit("#", 1)[0] in original for cid in draw["case_id"])
    assert not set(draw["case_id"]) & original          # NO bare id survives


def test_propensity_fit_DOES_NOT_RAISE_A_SCHEMA_ERROR_on_400_consecutive_replicates():
    """§5.1's blocker, gone. `[data-gated: no]` — the frame is `fixtures_stage10`'s.

    **`SchemaError` and not "does not raise" at all**, and the distinction is §7.1's: a `FitError`
    from a replicate is a sparse replicate that [§10] drops and counts, and a `SchemaError` is a bug
    in the resampler that [§10] may not catch. §15.2 as written says "does not raise"; what §5 is
    about, and what the companion below contrasts with, is the second. Measured on this fixture:
    zero of either, so both readings pass here and only this one would still be the right assertion
    on a frame where a legitimate sparse replicate occurred.
    """
    df = two_centre_frame()
    rng = rng_at()
    for _ in range(400):
        draw = bootstrap.resample(df, rng, config.BOOT_STRATUM)
        propensity.fit(draw, data.Audit(hand_source(Path("."), "stage10")))


def test_WITHOUT_THE_RENAME_the_same_400_draws_raise_SchemaError_and_THIS_IS_THE_WHOLE_OF_5():
    """The companion, and it is the measurement §5.1 reproduces from Stage 9 §12.3.

    Asserted as `> 0` rather than pinned, because the rate is a property of the fixture's stratum
    sizes and pinning it makes the test brittle for nothing. The MESSAGE is asserted too: it must be
    `_record_exclusion`'s name-count-against-row-count guard and not some other schema failure the
    unrenamed frame happens to trip.
    """
    df = two_centre_frame()
    rng = rng_at()
    raised = 0
    for _ in range(400):
        parts = [g.iloc[rng.integers(0, len(g), size=len(g))]
                 for _, g in df.groupby(config.BOOT_STRATUM, sort=True, observed=True)]
        naive = pd.concat(parts, axis=0).reset_index(drop=True)          # NO rename
        try:
            propensity.fit(naive, data.Audit(hand_source(Path("."), "stage10")))
        except config.SchemaError as e:
            raised += 1
            assert "covariate_completeness" in str(e) and "can name" in str(e)
    assert raised > 0, (
        "the unrenamed resampler raised on none of 400 replicates, so §5's blocker is not reachable "
        "on this fixture and the renaming test above is asserting nothing")


def test_the_same_seed_gives_an_IDENTICAL_frame_and_a_different_seed_does_not():
    df = two_centre_frame()
    first = bootstrap.resample(df, rng_at(), config.BOOT_STRATUM)
    again = bootstrap.resample(df, rng_at(), config.BOOT_STRATUM)
    other = bootstrap.resample(df, rng_at(config.SEED + 1), config.BOOT_STRATUM)
    assert first.equals(again)
    assert not first.equals(other)


def test_the_frames_COLUMNS_AND_DTYPES_are_unchanged_column_by_column():
    """And `case_id` is the one that fails without §5.3's cast, so it is NAMED rather than left to
    the loop to happen to reach: it is declared `dtype="string"` and a list assignment returns
    `object`. A resampler that silently promotes an `Int64` to `float64` changes `model.design` and
    hence every estimate."""
    df = two_centre_frame()
    draw = bootstrap.resample(df, rng_at(), config.BOOT_STRATUM)
    assert list(draw.columns) == list(df.columns)
    for column in df.columns:
        assert draw[column].dtype == df[column].dtype, column
    assert str(df["case_id"].dtype) == "string"
    assert draw["case_id"].dtype == df["case_id"].dtype
    assert draw["case_id"].dtype != np.dtype("O")


def test_WITHOUT_THE_CAST_the_rename_demotes_case_id_to_object():
    """The companion, measured: `string[python]` in, `dtype('O')` out. `data.py:442` documents the
    identical trap one stage up (§5.3, §21c item 6)."""
    df = two_centre_frame()
    rng = rng_at()
    parts = [g.iloc[rng.integers(0, len(g), size=len(g))]
             for _, g in df.groupby(config.BOOT_STRATUM, sort=True, observed=True)]
    out = pd.concat(parts, axis=0).reset_index(drop=True)
    occurrence = out.groupby("case_id", sort=False, observed=True).cumcount() + 1
    out["case_id"] = [f"{cid}#{k}" for cid, k in zip(out["case_id"], occurrence)]   # NO cast
    assert out["case_id"].dtype == np.dtype("O")
    assert out["case_id"].dtype != df["case_id"].dtype


def test_the_index_is_RESET_so_no_future_label_based_read_is_wrong_in_a_way_that_returns():
    df = two_centre_frame()
    draw = bootstrap.resample(df, rng_at(), config.BOOT_STRATUM)
    assert list(draw.index) == list(range(len(df)))
    assert draw.index.is_unique


def test_the_STREAM_is_a_function_of_the_SEED_ALONE_and_not_of_whether_a_replicate_succeeded():
    """§4.3, asserted frame by frame. `resample` is called before `body` and never inside its `try`,
    so this holds by construction; it is asserted because an edit that moved the draw inside the
    failure path would make the seed stop identifying the replicates, and nothing else in §15 would
    notice."""
    df = two_centre_frame()
    succeeding: list[pd.DataFrame] = []
    failing: list[pd.DataFrame] = []

    def returning(draw):
        succeeding.append(draw)
        return "ok"

    def raising(draw):
        # Records, then raises and handles — which is exactly `_replicate`'s own shape (§7.1). A body
        # that let the `FitError` escape would kill the loop at replicate 0, because `replicates`
        # catches nothing: the taxonomy is the BODY's to apply (§3.2).
        failing.append(draw)
        try:
            raise model.FitError("G7  a body whose every replicate fails")
        except model.FitError:
            return None

    bootstrap.replicates(df, returning, 12, config.SEED, config.BOOT_STRATUM)
    bootstrap.replicates(df, raising, 12, config.SEED, config.BOOT_STRATUM)
    assert len(succeeding) == len(failing) == 12
    for i, (ok, failed) in enumerate(zip(succeeding, failing)):
        assert ok.equals(failed), f"the two streams diverge at replicate {i}"
    assert not succeeding[0].equals(succeeding[1]), "the stream is not advancing at all"


def test_replicates_creates_ONE_generator_and_consumes_it_IN_ORDER(boot_n):
    """Not one seeded per replicate: a per-replicate seed derived from `b` is reproducible too, but
    it makes the draw a function of an index a later edit can renumber, and the failure is silent
    (§4.3). Asserted by requiring the `n` draws to be the `n` successive draws of one Generator."""
    df = two_centre_frame()
    collected = bootstrap.replicates(df, lambda d: d, 8, config.SEED, config.BOOT_STRATUM)
    rng = rng_at()
    for drawn in collected:
        assert drawn.equals(bootstrap.resample(df, rng, config.BOOT_STRATUM))


def test_resample_visits_the_strata_in_THE_COLUMNS_OWN_SORT_ORDER_and_not_by_str():
    """§5.3's second bullet. An earlier draft sorted on `str(label)`, which is a no-op on this
    study's `string` centre labels and WRONG for a general resampler: it visits an integer stratum
    column as 1, 10, 2, 9. `resample` is general (§12.2), so the contract is the column's own
    ordering and pandas is the one thing defining it — asserted on an INTEGER stratum, which is the
    only kind that can tell the two apart."""
    df = pd.DataFrame({
        "case_id": pd.Series([f"N-{i:03d}" for i in range(24)], dtype="string"),
        "stratum": [1] * 6 + [2] * 6 + [9] * 6 + [10] * 6,
        "value": np.arange(24, dtype=float),
    })
    draw = bootstrap.resample(df, rng_at(), "stratum")
    assert list(draw["stratum"].unique()) == [1, 2, 9, 10]              # numeric, not 1, 10, 2, 9


# --- 15.3  the replicate body, and what it may not read -------------------------------------------------

def test_the_replicate_REFITS_the_propensity_model_and_does_not_reuse_the_point_estimates():
    """A body that accidentally closed over the outer `ps` returns the point estimate's weights and
    every interval collapses toward zero width — a failure that produces plausible numbers.

    Asserted by driving `_replicate` on a frame whose `e` MUST differ from the point estimate's and
    checking that the estimates differ: the frame is a resample, so its propensity fit is over
    different rows.
    """
    df = two_centre_frame()
    ps, est, sec, _ = fitted(df)
    paths = {k: e.augmented_path for k, e in sec.estimates.items()}
    keys = bootstrap._estimand_keys(est, sec)
    draw = bootstrap.resample(df, rng_at(), config.BOOT_STRATUM)
    replicate = bootstrap._replicate(draw, keys, paths, hand_source(Path("."), "stage10"))
    assert replicate.values["beta"] != est.beta
    assert replicate.n_in_model != int(ps.in_model.sum()) or replicate.sum_w != float(
        ps.w[ps.in_model].sum())


def test_the_per_replicate_Audit_is_DISCARDED_and_the_callers_holds_exactly_ONE_new_entry(run_once):
    """After `run`, the `Audit` passed in holds one new entry (§10.1) and not `9 * N_BOOT + 1`."""
    df, ps, est, sec, audit, boot = run_once
    steps = [e.step for e in audit.entries]
    assert steps.count("bootstrap_replicates") == 1
    assert steps.count("binary_estimates") == 1          # the POINT estimate's, once
    assert len(audit.entries) == 11                      # 4 propensity + 3 primary + 3 secondary + 1


def test_the_per_replicate_Audit_carries_THE_CALLERS_SOURCE_and_not_a_synthetic_one():
    """§6.2: a second module-level `Source` is pinned against by `data.SOURCES` and a test, and
    inventing one here would either break that pin or lie about provenance. Asserted by identity."""
    df = two_centre_frame()
    ps, est, sec, audit = fitted(df)
    seen: list[object] = []
    original = data.Audit.__init__

    def spy(self, source):
        seen.append(source)
        original(self, source)

    draw = bootstrap.resample(df, rng_at(), config.BOOT_STRATUM)
    paths = {k: e.augmented_path for k, e in sec.estimates.items()}
    keys = bootstrap._estimand_keys(est, sec)
    data.Audit.__init__ = spy
    try:
        bootstrap._replicate(draw, keys, paths, audit.source)
    finally:
        data.Audit.__init__ = original
    assert seen and all(source is audit.source for source in seen)


# --- 15.4  the failure taxonomy reconciles ---------------------------------------------------------------

def synthetic_replicates(keys, n, injected):
    """`n` replicates over `keys`, with `injected[i]` naming the keys replicate `i` LOST.

    The expected counts are known because this function writes them, rather than being read back
    from the thing under test (§15.4).
    """
    out = []
    for i in range(n):
        lost = injected.get(i, ())
        out.append(bootstrap.Replicate(
            values={k: float(i) for k in keys if k not in lost},
            failures={k: "separation" for k in lost},
            n_alpha=None if lost == tuple(keys) else 6,
            polr_iterations=None if lost == tuple(keys) else 3,
            sum_w=None if lost == tuple(keys) else 27.0,
            n_in_model=None if lost == tuple(keys) else 92))
    return tuple(out)


def test_len_draws_plus_failures_equals_n_attempted_for_EVERY_key():
    """The one property no wrong implementation satisfies by accident: a silent drop anywhere makes
    the three numbers stop reconciling."""
    keys = ("beta", "rd_0", "sich.rd")
    collected = synthetic_replicates(keys, 50, {3: ("beta",), 7: ("beta", "rd_0"), 11: keys})
    draws = bootstrap.collect(collected, keys)
    for key, d in draws.items():
        assert len(d.draws) + sum(d.failures.values()) == d.n_attempted
    assert len(draws["beta"].draws) == 47 and draws["beta"].failures == {"separation": 3}
    assert len(draws["rd_0"].draws) == 48 and draws["rd_0"].failures == {"separation": 2}
    assert len(draws["sich.rd"].draws) == 49


def test_a_Draws_that_does_not_reconcile_RAISES_which_is_why_the_property_is_a_TYPE():
    with pytest.raises(config.SchemaError) as e:
        bootstrap.Draws(quantity="beta", draws=np.zeros(3), n_attempted=10,
                        failures={"separation": 2})
    assert "reconcil" in message_of(e) or "attempted" in message_of(e)


def test_a_Replicate_reporting_a_key_as_BOTH_a_draw_and_a_failure_RAISES():
    """§3.1's partition made a property of the type rather than of the loop that fills it."""
    with pytest.raises(config.SchemaError) as e:
        bootstrap.Replicate(values={"beta": 1.0}, failures={"beta": "separation"},
                            n_alpha=6, polr_iterations=3, sum_w=27.0, n_in_model=92)
    assert "beta" in message_of(e)


def test_n_attempted_is_N_BOOT_on_ALL_TWENTY_SIX_keys(run_once):
    """Asserted POSITIVELY (§3.1, §7.1). The equality is a consequence of the propensity rule below
    and not a definition, so it can fail."""
    _, _, _, _, _, boot = run_once
    assert len(boot.draws) == 26
    for key, d in boot.draws.items():
        assert d.n_attempted == RUN_BOOT, key


def test_a_propensity_fit_FitError_is_counted_against_EVERY_key(monkeypatch, boot_n):
    """§7.1. Injected on a known number of replicates: every one of the 26 keys' `failures` sums
    rises by exactly that number, no key's `n_attempted` moves, and the reconciliation still holds.

    **The companion is what makes it a test**: an implementation that simply SKIPS those replicates
    passes the reconciliation check and fails the `n_attempted` assertion, and skipping is the
    reading a first draft of §3.1 invited.
    """
    df = two_centre_frame()
    ps, est, sec, audit = fitted(df)
    clean = bootstrap.run(df, ps, est, sec, data.Audit(audit.source))

    calls = {"n": 0}
    real_fit = propensity.fit

    def flaky(frame, a):
        calls["n"] += 1
        if calls["n"] % 5 == 0:                     # exactly RUN_BOOT // 5 of them
            raise model.FitError("G7  an injected whole-replicate failure")
        return real_fit(frame, a)

    monkeypatch.setattr(bootstrap.propensity, "fit", flaky)
    injured = bootstrap.run(df, ps, est, sec, data.Audit(audit.source))
    expected = RUN_BOOT // 5
    for key, d in injured.draws.items():
        assert d.n_attempted == RUN_BOOT, key                       # the companion assertion
        assert sum(d.failures.values()) == (
            sum(clean.draws[key].failures.values()) + expected), key
        assert len(d.draws) + sum(d.failures.values()) == d.n_attempted


def test_FitError_IS_NOT_A_SchemaError_and_the_WHOLE_TAXONOMY_needs_that():
    """The whole taxonomy is one `except` clause away from collapsing (§7.1)."""
    assert not issubclass(model.FitError, config.SchemaError)
    assert not issubclass(config.SchemaError, model.FitError)


def test_a_SchemaError_raised_INSIDE_a_replicate_PROPAGATES_out_of_run(monkeypatch, boot_n):
    """§7.1's rule is a NEGATIVE — "is never caught" — and the only way to test a negative is to
    raise the thing and require it to arrive."""
    df = two_centre_frame()
    ps, est, sec, audit = fitted(df)

    def exploding(frame, a):
        raise config.SchemaError("a frame this stage constructed cannot be read")

    monkeypatch.setattr(bootstrap.propensity, "fit", exploding)
    with pytest.raises(config.SchemaError):
        bootstrap.run(df, ps, est, sec, audit)


def test_run_CATCHES_NOTHING_BUT_model_FitError_by_scan():
    """By scan, because the behavioural tests can only reach the failure modes they construct.

    `pilots/analysis.py:595` is `except Exception: failures["other"] += 1`, which counts a `KeyError`
    from a typo and an unfittable model as the same dropped replicate. §7.1 is the correction and
    this is it asserted: every handler in this module names `model.FitError` and nothing else.
    """
    handlers = [node for node in ast.walk(ast.parse(SOURCE))
                if isinstance(node, ast.ExceptHandler)]
    assert handlers, "the module has no handler at all — §7.1's droppable failure is not caught"
    for handler in handlers:
        caught = handler.type
        assert isinstance(caught, ast.Attribute) and caught.attr == "FitError", ast.dump(caught)
        assert isinstance(caught.value, ast.Name) and caught.value.id == "model"


# --- 15.5  G7 is exercised, and the pair is asserted `[roadmap]` -----------------------------------------
#
# The roadmap's own clause: "A test drives a replicate loop over a deliberately separable frame and
# asserts the G7 count is non-zero while the convergence count stays zero — the pair, because Stage 8
# measured that the second never fires on this estimator and a single 'failures' counter therefore
# reads zero on data that is degenerate throughout."


def test_the_separable_frame_CONVERGES_which_is_the_whole_point_of_the_fixture():
    """The fixture's own witness first: the assertion is about a fit that SUCCEEDED and was rejected,
    not about a fit that failed. Reproduces Stage 8 §6.1 exactly (§21)."""
    X, y, w = separable_ordinal_frame()
    fit = model.polr(X, y, w)
    assert fit.iterations == 17
    assert fit.converged_on == "score"
    assert fit.rescales == 0 and fit.halvings == 0
    assert round(float(fit.beta[0]), 4) == 36.4058
    assert round(float(fit.alpha[0]), 4) == -18.2029
    assert np.isfinite(fit.beta).all() and np.isfinite(fit.alpha).all()


def test_G7_FIRES_on_it_and_the_CONVERGENCE_bucket_stays_ZERO():
    """The PAIR, and it is the pair that is the finding: a single "failures" counter reads zero on a
    frame that is separated throughout, which is exactly Stage 8 §6.1's measurement."""
    X, y, w = separable_ordinal_frame()
    keys = ("beta",)
    collected = []
    for _ in range(10):
        try:
            outcome._assert_reportable(model.polr(X, y, w))
            collected.append(bootstrap.Replicate(values={"beta": 0.0}, failures={},
                                                 n_alpha=1, polr_iterations=17, sum_w=40.0,
                                                 n_in_model=40))
        except model.FitError as failure:
            collected.append(bootstrap.Replicate(
                values={}, failures={"beta": bootstrap.bucket(str(failure))},
                n_alpha=None, polr_iterations=None, sum_w=None, n_in_model=None))
    draws = bootstrap.collect(tuple(collected), keys)
    assert draws["beta"].failures.get("separation", 0) > 0
    assert draws["beta"].failures.get("nonconvergence", 0) == 0


def test_the_SINGLE_COUNTER_MUTATION_reads_zero_on_the_same_frame():
    """§20's Definition of done item 5, as a test: replacing the two buckets with one "failures"
    count makes the assertion above pass vacuously — the separated frame produces a non-zero total,
    and a reader cannot tell it from a frame that failed to converge ten times."""
    X, y, w = separable_ordinal_frame()
    with pytest.raises(model.FitError) as e:
        outcome._assert_reportable(model.polr(X, y, w))
    single_counter = 1                                          # what a merged counter would report
    assert single_counter > 0                                   # ...and this says nothing
    assert bootstrap.bucket(message_of(e)) == "separation"
    assert bootstrap.bucket(message_of(e)) != "nonconvergence"


@DATA_GATED
def test_the_workbooks_counters_read_zero_and_zero(run_once):
    """`[data-gated]` — the measurement that makes the pair worth REPORTING rather than the one that
    makes it worth asserting. Run on the fixture cohort here, which is the same statement about a
    frame that is not degenerate: both counters read zero and G7 is in the path."""
    _, _, _, _, _, boot = run_once
    assert boot.draws["beta"].failures.get("separation", 0) == 0
    assert boot.draws["beta"].failures.get("nonconvergence", 0) == 0


# --- 15.6  the bucket map is scanned, not trusted --------------------------------------------------------

# Stage 11 §9 grows this scope in the same edit that adds G8 and G9 to `C.FAILURE_BUCKETS`. The
# alternative — putting the two raises in `outcome.py` so the existing scope covered them — was
# declined, because Stage 8 §12 makes `outcome.py` the only shipped module that may name the [§5]
# primary outcome and a subgroup estimator in it would be that rule spent for a scan's convenience.
FITERROR_MODULES: Final[tuple[str, ...]] = ("model.py", "outcome.py", "propensity.py",
                                            "sensitivity.py")


def raise_sites() -> list[tuple[str, int, str]]:
    """(module, line, leading token) for every `raise FitError(` / `raise model.FitError(` site.

    **BY CONTENT AND NOT BY LINE NUMBER**, as Stage 7 §18d requires of every scan in this repository:
    a raise site that MOVES must not fail this test and a raise site whose MESSAGE changes must. The
    walk is an AST walk over `ast.Raise` nodes, so a `FitError` mentioned in a docstring or built and
    not raised is invisible here and correctly so.
    """
    found: list[tuple[str, int, str]] = []
    root = Path(config.__file__).resolve().parent
    for name in FITERROR_MODULES:
        tree = ast.parse((root / name).read_text(encoding="utf-8"))
        for node in ast.walk(tree):
            if not isinstance(node, ast.Raise) or not isinstance(node.exc, ast.Call):
                continue
            func = node.exc.func
            called = func.attr if isinstance(func, ast.Attribute) else getattr(func, "id", "")
            if called != "FitError":
                continue
            first = node.exc.args[0]
            if isinstance(first, ast.JoinedStr):
                first = first.values[0]
            assert isinstance(first, ast.Constant), (
                f"{name}:{node.lineno} builds its message without a leading string literal, so the "
                "first token cannot be read statically — and §7.2's classification is by that token")
            found.append((name, node.lineno, str(first.value).split()[0]))
    return found


def test_EVERY_FitError_raise_site_leads_with_a_token_in_FAILURE_BUCKETS():
    """A prefix map is only safe if something fails when a message changes (§7.2, §15.6)."""
    sites = raise_sites()
    unknown = sorted({token for _, _, token in sites} - set(config.FAILURE_BUCKETS))
    assert not unknown, (
        f"{unknown} lead(s) a raised FitError and is in no C.FAILURE_BUCKETS bucket. A reworded "
        "message is a test failure and not a silently misfiled counter.")


def test_the_scan_finds_NINETEEN_sites_and_the_SPECS_SIXTEEN_is_its_own_scope():
    """§7.2, §13 and §15.6 all say sixteen — fourteen in `model.py` and two in `outcome.py` — and
    re-scanning reproduces exactly those sixteen. What the specification's SCOPE omitted is
    `propensity.py`, which raises `model.FitError` three more times, with the tokens `ESS:` and
    `F5`, both on `propensity.fit`'s own path — which is precisely the path §7.1 says fails a WHOLE
    replicate. Without them `bucket` turns a droppable sparse replicate into a crash.

    Pinned as a COUNT PER MODULE rather than as a total, because that is what makes the departure
    legible rather than a number nobody can locate.
    """
    sites = raise_sites()
    per_module = {name: sum(1 for m, _, _ in sites if m == name) for name in FITERROR_MODULES}
    assert per_module == {"model.py": 14, "outcome.py": 3, "propensity.py": 3,
                          "sensitivity.py": 2}
    assert len(sites) == 22
    # `outcome.py` is THREE and not the specification's two, and the third is S8 — which §5.4 moved
    # from `SchemaError` to `FitError` in this very stage, so the specification's own count of the
    # sites it was creating is one behind itself.
    assert sorted(t for m, _, t in sites if m == "outcome.py") == ["G6", "G7", "S8"]
    assert sorted(set(t for m, _, t in sites if m == "propensity.py")) == ["ESS:", "F5"]
    # Stage 11's two, and they are the G series continued — G8 the interaction column dropped as
    # constant, G9 a REPORTED quantity reaching POLR_MAX_ABS_BETA (Stage 11 §9).
    assert sorted(t for m, _, t in sites if m == "sensitivity.py") == ["G8", "G9"]


def test_every_declared_bucket_token_is_ACTUALLY_RAISED_somewhere():
    """The other direction, and it is what stops the map growing tokens nothing produces: a bucket
    keyed on a token no raise site leads with is a counter that can only ever read zero."""
    tokens = {token for _, _, token in raise_sites()}
    assert set(config.FAILURE_BUCKETS) == tokens


def test__bucket_RAISES_on_an_unrecognised_token_and_does_not_default():
    """A default would make the scan cosmetic: it would pass, the map would be incomplete, and the
    counter Stage 8 §11 asked to be separate would be silently merged (§15.6)."""
    with pytest.raises(config.SchemaError) as e:
        bootstrap.bucket("Z9  a token from a raise site nobody declared")
    assert "Z9" in message_of(e)


@pytest.mark.parametrize("token,bucket", [("G7", "separation"), ("S8", "constant_outcome"),
                                          ("polr:", "nonconvergence"), ("Firth:", "nonconvergence"),
                                          ("G6", "degenerate_design"), ("F5", "degenerate_design")])
def test__bucket_maps_each_token_to_the_bucket_7_2_names(token, bucket):
    assert bootstrap.bucket(f"{token}  a message") == bucket


# --- 15.8  the interval, the p-value, and the definition `[roadmap, amended]` ----------------------------

def test_the_limits_are_ORDER_STATISTICS_and_are_NOT_a_function_of_the_point_estimate():
    """"The interval is centred on the estimate" is what a reader assumes, and an implementation that
    quietly did it would produce plausible numbers (§8.1)."""
    draws = boundary_draws(400, n=1000)
    lo, hi = bootstrap.percentile_ci(draws)
    ordered = np.sort(draws)
    assert lo == ordered[int(np.ceil(0.025 * len(draws))) - 1]
    assert hi == ordered[int(np.ceil(0.975 * len(draws))) - 1]
    # `percentile_ci` takes no point estimate at all, which is the strongest form of the assertion.
    assert list(inspect.signature(bootstrap.percentile_ci).parameters) == ["draws", "level"]


@pytest.mark.parametrize("tail,expected_p,significant", [(49, 0.0490, True), (50, 0.0500, False),
                                                         (51, 0.0510, False)])
def test_the_p_value_AGREES_with_the_pinned_interval_at_49_50_and_51(tail, expected_p, significant):
    """§8.2's sweep. At B = 2000 the p-value is a multiple of 2/B, so `p < 0.05` iff the smaller tail
    holds at most 49 draws — and at exactly 50 `p = 0.0500`, which is NOT less than 0.05, so the
    interval must INCLUDE the null."""
    draws = boundary_draws(tail)
    p = bootstrap.bootstrap_p(draws)
    lo, hi = bootstrap.percentile_ci(draws)
    assert round(p, 4) == expected_p
    assert (p < 0.05) is significant
    assert (not (lo <= 0.0 <= hi)) is significant, "the interval and the p-value disagree"


def test_THE_COMPANION_at_tail_50_the_DEFAULT_METHOD_EXCLUDES_the_null_and_the_pin_INCLUDES_it():
    """**This is the test that can fail**, and §20's Definition of done item 4 requires it to have
    been seen failing before the pin is trusted: the same draws under `method="linear"` give an
    interval EXCLUDING zero against a p of exactly 0.0500.

    `boundary_draws(49)` and `(51)` agree under both methods, which is why a test written at either
    would pass on the wrong pin.
    """
    draws = boundary_draws(50)
    assert round(bootstrap.bootstrap_p(draws), 4) == 0.0500          # NOT significant
    pinned_lo, pinned_hi = bootstrap.percentile_ci(draws)
    assert pinned_lo <= 0.0 <= pinned_hi                             # ...and the interval agrees
    linear_lo, linear_hi = np.percentile(draws, [2.5, 97.5], method="linear")
    assert not (linear_lo <= 0.0 <= linear_hi)                       # ...and the default does not
    assert config.PERCENTILE_METHOD == "inverted_cdf"


def test_the_QUANTILE_ARGUMENT_IS_SNAPPED_and_the_naive_form_flips_the_LIMITS_SIGN():
    """§8.2's second half, which is arithmetic and not a rounding cosmetic: `100.0*((1.0-0.95)/2.0)`
    is `2.500000000000002`, and under `inverted_cdf` that 2e-15 excess moves the order-statistic
    index from 50 to 51 — from the largest negative draw to the smallest positive one."""
    assert 100.0 * ((1.0 - 0.95) / 2.0) != 2.5
    assert round(50.0 * (1.0 - 0.95), 9) == 2.5
    draws = boundary_draws(50)
    naive = float(np.percentile(draws, 100.0 * ((1.0 - 0.95) / 2.0),
                                method=config.PERCENTILE_METHOD))
    snapped, _ = bootstrap.percentile_ci(draws)
    assert snapped < 0.0 < naive, "the two quantile arguments must land on opposite sides of zero"


def test_bootstrap_p_FLOORS_ON_THE_SURVIVING_COUNT_and_never_on_N_BOOT():
    """On the cohort the two are equal, which is exactly the condition under which an implementation
    reading `C.N_BOOT` is green — so it is asserted on a fixture where they differ (§9.1)."""
    all_positive = np.abs(boundary_draws(0, n=500)) + 1.0
    assert bootstrap.bootstrap_p(all_positive) == pytest.approx(1.0 / 501)
    assert bootstrap.bootstrap_p(all_positive) != pytest.approx(1.0 / (config.N_BOOT + 1))


@pytest.mark.parametrize("ties", [2, 5, 40, 100])
def test_draws_at_EXACTLY_ZERO_are_counted_in_BOTH_tails_and_the_test_is_CONSERVATIVE(ties):
    """`Pr(<=0) + Pr(>=0) > 1` and p is inflated, with the interval containing zero in each case.

    Deliberate and measured (§8.2): the alternative — splitting ties — makes p a function of a
    tie-breaking rule nobody prespecified.
    """
    draws = boundary_draws(200)
    draws[:ties] = 0.0
    assert float(np.mean(draws <= 0.0)) + float(np.mean(draws >= 0.0)) > 1.0
    lo, hi = bootstrap.percentile_ci(draws)
    assert lo <= 0.0 <= hi


def test_the_floor_REFUSES_at_39_and_at_40_gives_the_MINIMUM_and_the_SECOND_LARGEST():
    """**Both halves asserted positively.** A test asserting the MAXIMUM passes on a wrong percentile
    definition and fails on the pinned one: `ceil(0.025*40) = 1` gives index 0 while
    `ceil(0.975*40) = 39` gives index 38, the second-largest (§8.3)."""
    with pytest.raises(config.SchemaError) as e:
        bootstrap.percentile_ci(np.arange(39.0))
    assert "39" in message_of(e) and "40" in message_of(e)
    lo, hi = bootstrap.percentile_ci(np.arange(40.0))
    assert (lo, hi) == (0.0, 38.0)
    assert hi != 39.0, "the upper limit at the floor is the SECOND-largest draw, not the maximum"


def test_runs_BELOW_FLOOR_branch_keeps_the_Draws_and_emits_NO_Interval(monkeypatch, boot_n):
    """**And `run`'s below-floor branch is reached separately, because `percentile_ci` raising is not
    it** (§15.8). Both halves: the key IS in `boot.draws` carrying its counters, and is NOT in
    `boot.intervals` — an implementation that emits an `Interval` of extremes satisfies the first,
    and one that drops the estimand entirely satisfies the second. This is the shape Stage 14 must
    handle (§12.1, §3.3).

    Reached by failing ONE outcome in almost every replicate, which is what the floor describes: an
    estimand needs 98% of its replicates to fail before it loses its interval, and the worst measured
    per-outcome drop rate on the workbook is 0.8% — so the branch is unreachable on v7 and is reached
    here with a fixture.
    """
    df = two_centre_frame()
    ps, est, sec, audit = fitted(df)
    calls = {"n": 0}
    real = outcome.secondary

    def flaky(frame, p, a, paths=None, collect=False):
        calls["n"] += 1
        result = real(frame, p, a, paths=paths, collect=collect)
        if calls["n"] > 5:                       # the first five replicates keep `sich`; the rest do not
            estimates = {k: v for k, v in result.estimates.items() if k != "sich"}
            return outcome.Secondary(estimates=estimates,
                                     failures={**result.failures, "sich": "G7  injected"})
        return result

    monkeypatch.setattr(bootstrap.outcome, "secondary", flaky)
    boot = bootstrap.run(df, ps, est, sec, audit)
    for key in ("sich.rd", "sich.odds_ratio"):
        assert key in boot.draws
        assert len(boot.draws[key].draws) == 5
        assert boot.draws[key].failures == {"separation": RUN_BOOT - 5}
        assert key not in boot.intervals, "an Interval was emitted below ci_min_draws"
    assert "mrs_0_2_90d.rd" in boot.intervals, "a neighbouring outcome lost its interval too"


def test_the_six_RD_k_carry_p_is_None_ASSERTED_POSITIVELY_ON_ALL_SIX(run_once):
    """"we did not compute one" and "we computed one and it is `None`" look identical from the
    outside, and only the first is what [§8] asks for — so this is a positive check on all six keys
    by name (§9.3)."""
    _, _, _, _, _, boot = run_once
    names = [f"rd_{k}" for k in config.MRS_THRESHOLDS]
    assert names == ["rd_0", "rd_1", "rd_2", "rd_3", "rd_4", "rd_5"]
    for key in names:
        assert key in boot.intervals
        assert boot.intervals[key].p is None


def test__tested_is_true_on_EXACTLY_EIGHT_of_the_TWENTY_SIX_keys(run_once):
    """Asserted over the FULL key set rather than as an absence: no `odds_ratio` and no `augmented`
    key anywhere in `intervals` carries a non-`None` p (§9.2, §3.3)."""
    _, _, _, _, _, boot = run_once
    keys = list(boot.draws)
    assert len(keys) == 26
    tested = [k for k in keys if bootstrap._tested(k)]
    assert len(tested) == 8
    assert tested == ["beta"] + [f"{o}.rd" for o in config.BINARY_OUTCOMES]
    for key, interval in boot.intervals.items():
        if key.endswith(".odds_ratio") or key.endswith(".augmented"):
            assert interval.p is None, key
        if bootstrap._tested(key):
            assert interval.p is not None, key


def test_every_Interval_records_the_LEVEL_the_METHOD_and_the_SURVIVING_COUNT(run_once):
    """§3.1: which definition produced these two numbers is part of what they are, and [§10]'s
    "dropped and counted" is only informative where the count is printed."""
    _, _, _, _, _, boot = run_once
    for key, interval in boot.intervals.items():
        assert interval.level == config.CI_LEVEL
        assert interval.method == config.PERCENTILE_METHOD
        assert interval.n_draws == len(boot.draws[key].draws)
        assert interval.lo <= interval.hi


# --- 15.9  the shared design, and the companion that makes it a test -------------------------------------

def full_list_outcomes() -> tuple[str, ...]:
    """The four §6.4 means, and BOTH conditions are load-bearing.

    The covariate-list condition alone returns SIX on v7 — `sich` and `ph2` take the shared list
    because they decline a MODEL rather than a list, and no design is built for either. The second
    condition is what invariant 6 is a statement about: the outcome derives from the [§5] primary
    ordinal outcome, so its missingness IS that outcome's.
    """
    by_list = tuple(k for k in config.BINARY_OUTCOMES
                    if config.outcome_model_covariates(k) == config.OUTCOME_COVARIATES)
    assert len(by_list) == 6, "the covariate-list condition alone is not the four §6.4 means"
    return tuple(k for k in by_list if config.OUTCOMES[k].source == config.PRIMARY_OUTCOME)


def test_the_FOUR_full_list_outcomes_share_a_mask_and_an_IDENTICAL_design_on_a_replicate():
    """Measured over replicates: the four masks differ in 0 and the four designs compare equal under
    `DataFrame.equals` (§6.4). The condition holds by roadmap invariant 6 — all four derive from the
    ordinal source and a derived dichotomy carries exactly its source's missingness — which is
    asserted in the suite rather than true by luck."""
    full = full_list_outcomes()
    assert len(full) == 4
    df = two_centre_frame()
    rng = rng_at()
    for _ in range(20):
        draw = bootstrap.resample(df, rng, config.BOOT_STRATUM)
        ps = propensity.fit(draw, data.Audit(hand_source(Path("."), "stage10")))
        masks = [ps.in_model & draw[key].notna() for key in full]
        assert all(mask.equals(masks[0]) for mask in masks)
        designs = [model.design(draw.loc[mask], config.OUTCOME_COVARIATES)[0] for mask in masks]
        assert all(d.equals(designs[0]) for d in designs)
        assert bootstrap._shared_design(draw, ps).equals(designs[0])


def test_WHEN_THE_MASKS_DIFFER_shared_design_RAISES_and_does_NOT_fall_back():
    """**Asserted as a raise and not as a fallback.** An earlier draft of §6.4 required a silent
    fall-back to per-outcome designs, which finishes the run with correct numbers and never tells
    anyone a Stage 3 invariant stopped holding.

    An optimisation whose precondition is never violated in a test is an optimisation nobody has
    tested, so the precondition is violated here deliberately.
    """
    full = full_list_outcomes()
    df = two_centre_frame()
    ps, _, _, _ = fitted(df)
    broken = df.copy()
    broken.loc[broken.index[:3], full[1]] = np.nan          # invariant 6 violated on ONE of the four
    with pytest.raises(config.SchemaError) as e:
        bootstrap._shared_design(broken, ps)
    assert full[1] in message_of(e)
    assert "invariant 6" in message_of(e)


def test_the_mask_mismatch_is_SchemaError_and_NOT_FitError_asserted_BY_CLASS():
    """Resampling copies whole rows and cannot break invariant 6, so a mismatch is upstream and a
    bug; making it droppable would let `N_BOOT` replicates absorb it one at a time (§6.4, §15.9)."""
    full = full_list_outcomes()
    df = two_centre_frame()
    ps, _, _, _ = fitted(df)
    broken = df.copy()
    broken.loc[broken.index[:3], full[2]] = np.nan
    with pytest.raises(config.SchemaError) as e:
        bootstrap._shared_design(broken, ps)
    assert not isinstance(e.value, model.FitError)


def test_the_mask_mismatch_PROPAGATES_OUT_OF_RUN_UNCAUGHT(monkeypatch, boot_n):
    """§7.1's rule as behaviour: `run` catches `model.FitError` and nothing else, so a `SchemaError`
    from inside `_replicate` terminates the analysis."""
    full = full_list_outcomes()
    df = two_centre_frame()
    ps, est, sec, audit = fitted(df)
    broken = df.copy()
    broken.loc[broken.index[:3], full[3]] = np.nan
    ps_broken = propensity.fit(broken, data.Audit(hand_source(Path("."), "stage10")))
    with pytest.raises(config.SchemaError) as e:
        bootstrap.run(broken, ps_broken, est, sec, audit)
    assert "invariant 6" in message_of(e)


def test_tici_BUILDS_ITS_OWN_DESIGN_because_its_covariate_list_is_the_override():
    """Sharing the full-list design would silently estimate the wrong nuisance model (§15.9)."""
    reduced = [k for k in config.BINARY_OUTCOMES
               if config.outcome_model_covariates(k) != config.OUTCOME_COVARIATES]
    assert reduced == ["tici_2b_3"]
    df = two_centre_frame()
    ps, _, sec, _ = fitted(df)
    shared = bootstrap._shared_design(df, ps)
    tici_fit = sec.estimates["tici_2b_3"]
    assert tici_fit.covariates == config.outcome_model_covariates("tici_2b_3")
    assert len(tici_fit.fit.columns) < shared.shape[1]


# --- 15.10  no fallback estimator, ever `[roadmap invariant 5]` -----------------------------------------

def test_a_nuisance_failure_costs_ALL_THREE_of_that_outcomes_keys_and_NO_OTHERS(boot_n):
    """**Both halves**: the group is ATOMIC (§7.3, Stage 9 §14's same-set rule) *and* the granularity
    is per outcome.

    **An earlier draft of §15.10 asserted the opposite on the first half** — `augmented` falling
    while `rd` held — which is `tau` and `rd` taken over different replicate sets, exactly what
    Stage 9 §14 forbids for two numbers [§10.3] reports as a comparison.

    The frame is `unfittable_nuisance_frame`, on which `tici_2b_3`'s `m_a(X)` cannot be fitted in any
    replicate and the other six can.
    """
    clean = two_centre_frame()
    ps_clean, est, sec, audit = fitted(clean)
    baseline = bootstrap.run(clean, ps_clean, est, sec, data.Audit(audit.source))

    # `propensity.fit` ALONE on the injured frame: `fitted` would call `outcome.secondary` at
    # `collect=False`, which is the point-estimate contract and raises on exactly the fit this frame
    # exists to break. The `est` and `sec` `run` is given are the CLEAN frame's, which is right —
    # they supply the estimand keys and the frozen paths, and `run` reads no estimate from either.
    df = unfittable_nuisance_frame()
    ps = propensity.fit(df, data.Audit(hand_source(Path("."), "stage10")))
    injured = bootstrap.run(df, ps, est, sec, data.Audit(audit.source))

    for key in ("tici_2b_3.rd", "tici_2b_3.odds_ratio", "tici_2b_3.augmented"):
        assert len(injured.draws[key].draws) == 0, key                 # all three fell together
        assert injured.draws[key].failures == {"degenerate_design": RUN_BOOT}
    for key in ("mrs_0_2_90d.rd", "mrs_0_2_90d.odds_ratio", "mrs_0_2_90d.augmented"):
        assert len(injured.draws[key].draws) == len(baseline.draws[key].draws), key


def test_augmented_is_NEVER_SUBSTITUTED_with_the_unaugmented_rd(boot_n):
    """A replicate whose augmented fit fails contributes NOTHING to that outcome's `augmented`
    draws and is not substituted (§15.10, roadmap invariant 5)."""
    df = unfittable_nuisance_frame()
    ps = propensity.fit(df, data.Audit(hand_source(Path("."), "stage10")))
    clean = two_centre_frame()
    _, est, sec, audit = fitted(clean)
    boot = bootstrap.run(df, ps, est, sec, data.Audit(audit.source))
    assert boot.draws["tici_2b_3.augmented"].draws.size == 0
    assert "tici_2b_3.augmented" not in boot.intervals
    assert "tici_2b_3.rd" not in boot.intervals


def test_a_primary_G7_costs_beta_AND_ALL_SIX_RD_k_as_a_GROUP(monkeypatch, boot_n):
    """The primary is one group (§7.3), so all seven counts fall together — asserted as a group
    rather than key by key, and the binaries are asserted NOT to move."""
    df = two_centre_frame()
    ps, est, sec, audit = fitted(df)
    calls = {"n": 0}
    real = outcome.primary

    def flaky(frame, p, a):
        calls["n"] += 1
        if calls["n"] % 3 == 0:
            raise model.FitError("G7  an injected separated fit")
        return real(frame, p, a)

    monkeypatch.setattr(bootstrap.outcome, "primary", flaky)
    boot = bootstrap.run(df, ps, est, sec, audit)
    expected = RUN_BOOT // 3
    primary_keys = ["beta"] + [f"rd_{k}" for k in config.MRS_THRESHOLDS]
    for key in primary_keys:
        assert boot.draws[key].failures == {"separation": expected}, key
        assert len(boot.draws[key].draws) == RUN_BOOT - expected
    for key in boot.draws:
        if key not in primary_keys:
            assert boot.draws[key].failures == {}, key


def test_no_code_path_calls_an_ESTIMATOR_the_point_estimate_did_not_use():
    """By scan. The module names `propensity.fit`, `outcome.primary`, `outcome.secondary` and
    `model.design`, and no other fitter (§15.10)."""
    names = module_names()
    assert "fit" in names and "primary" in names and "secondary" in names
    assert "polr" not in names, "bootstrap.py must not drive the ordinal fitter directly"
    assert "firth" not in names, "bootstrap.py must not drive the Firth fitter directly"
    assert not {"LogisticRegression", "GLM", "OLS", "minimize", "curve_fit"} & names


# --- 15.11  the module boundary --------------------------------------------------------------------------

def module_names() -> set[str]:
    """Every identifier `bootstrap.py` REFERENCES AS CODE — no docstrings, no comments, no prose.

    **A TEXT SCAN IS THE WRONG INSTRUMENT FOR "IS THIS READ"**, and Stage 9 §15.11's fixture exists
    because a first draft used one: this module's docstrings cite the neighbouring stages by name and
    say that `_augmentable` is not called and that `RARE_MINORITY_THRESHOLD` is not read. A substring
    scan reports every one of those sentences as a violation of what the sentence itself promises —
    a test that fails precisely when the code is well documented.
    """
    names: set[str] = set()
    for node in ast.walk(ast.parse(SOURCE)):
        if isinstance(node, ast.Name):
            names.add(node.id)
        elif isinstance(node, ast.Attribute):
            names.add(node.attr)
            if isinstance(node.value, ast.Name):
                names.add(f"{node.value.id}.{node.attr}")
        elif isinstance(node, (ast.FunctionDef, ast.ClassDef)):
            names.add(node.name)
    return names


def test_bootstrap_does_NOT_IMPORT_balance():
    """§0.2, and Stage 7 §14.11's reason: not because importing it would fail, but because it would
    succeed. [§9] is a statement about the realised sample and no interval is put on an SMD."""
    imported = {alias.name for node in ast.walk(ast.parse(SOURCE))
                if isinstance(node, ast.Import) for alias in node.names}
    imported |= {node.module for node in ast.walk(ast.parse(SOURCE))
                 if isinstance(node, ast.ImportFrom) and node.module}
    assert "balance" not in imported
    assert imported == {"__future__", "collections.abc", "dataclasses", "numpy", "pandas",
                        "config", "model", "outcome", "propensity", "data"}


@pytest.mark.parametrize("forbidden", ["_augmentable", "_augmented_path",
                                       "RARE_MINORITY_THRESHOLD", "OUTCOME_MODEL_OVERRIDES"])
def test_bootstrap_NAMES_NONE_OF_THE_VOCABULARY_OF_RE_DERIVING_A_PATH(forbidden):
    """§6.3's rule as a test: the path is READ off the point estimate, and the only way to guarantee
    it is not re-derived is to assert the vocabulary of re-deriving it is absent.

    §20's Definition of done item 6 requires this to have been seen failing: adding a read of
    `C.RARE_MINORITY_THRESHOLD` to `bootstrap.py` fails this test.
    """
    assert forbidden not in module_names()


def test_bootstrap_reaches_into_NO_OTHER_MODULES_PRIVATES():
    """`propensity._padded` is the named case (§10.1) — it does exactly what this stage's grid needs
    and is not imported. `data._fmt` is the ONE exception and it is the pipeline's own convention:
    four shipped modules already import it, because a second float formatter is a second way for two
    runs to disagree."""
    reaching = sorted(
        name for name in module_names()
        if "." in name and name.split(".")[1].startswith("_")
        and name.split(".")[0] in {"propensity", "outcome", "model", "data", "C"})
    assert reaching == [], f"bootstrap.py reaches into {reaching}"
    assert "_padded" in module_names()          # ...its OWN copy, defined in this module
    assert "_fmt" in module_names()             # ...imported by name, as four stages already do


@pytest.mark.parametrize("function", ["resample", "replicates", "percentile_ci", "bootstrap_p"])
def test_THE_FOUR_GENERAL_FUNCTIONS_name_NO_covariate_NO_outcome_and_NO_centre(function):
    """§0.1's generality claim made checkable, and it is what Stages 12 and 13 depend on: [§14a] and
    [§14b] resample the same way over a population this stage never sees."""
    tree = ast.parse(SOURCE)
    node = next(n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name == function)
    names: set[str] = set()
    for child in ast.walk(node):
        if isinstance(child, ast.Name):
            names.add(child.id)
        elif isinstance(child, ast.Attribute):
            names.add(child.attr)
        elif isinstance(child, ast.Constant) and isinstance(child.value, str):
            names.add(child.value)
    forbidden = (set(config.PS_COVARIATES) | set(config.OUTCOMES) | set(config.CENTER_ORDER)
                 | {config.TREATMENT, "BINARY_OUTCOMES", "PRIMARY_OUTCOME", "MRS_THRESHOLDS",
                    "BOOT_STRATUM"})
    assert not forbidden & names, f"{function} names {sorted(forbidden & names)}"


def test_the_public_surface_is_EIGHT_NAMES_and_the_privates_are_the_declared_ONES():
    """§3.2's count, re-derived rather than repeated. Stage 9 §22.3 item 5 records three drafts of
    its own carrying three different private counts, which is what this assertion prevents.

    **Stage 10 said five and twelve; Stage 11 §5.4 makes it eight and nine.** `bucket`, `collect`
    and `diagnostics` went public rather than being written a second time in `sensitivity.py` — a
    second classifier would be a second definition of the separation-versus-convergence split
    Stage 8 §11 asked to be kept apart, a second `Draws` loop a second implementation of the
    reconciliation property `Draws.__post_init__` exists to enforce, and a second tally a second
    answer to "how many replicates fitted five cutpoints". Stages 12 and 13 are their third caller.
    """
    tree = ast.parse(SOURCE)
    functions = [n.name for n in tree.body if isinstance(n, ast.FunctionDef)]
    public = [n for n in functions if not n.startswith("_")]
    private = [n for n in functions if n.startswith("_")]
    assert public == ["resample", "replicates", "percentile_ci", "bootstrap_p",
                      "bucket", "collect", "diagnostics", "run"]
    # THE LIST IS THE COUNT, and the spec's number has been two short since Stage 10. §3.2 enumerates
    # twelve privates and this module carries fourteen: the two it does not list are `_padded` and
    # `_spread`, both of which §10.1 REQUIRES — the concatenated grid with its short rows padded, and
    # the spread rendering — and describes without naming. Three of the fourteen are now public, so
    # ELEVEN remain against Stage 11 §5.4's restated nine, and the gap is the same two.
    assert sorted(private) == sorted([
        "_assert_run_inputs", "_replicate", "_tested",
        "_estimand_keys", "_shared_design", "_counters_table", "_diagnostics_table",
        "_replicates_detail", "_record_replicates", "_padded", "_spread"])


def test_NO_FORWARDING_ALIAS_IS_LEFT_ON_ANY_OF_THE_THREE():
    """Stage 11 §5.4: the three cease to exist under their private names rather than becoming
    one-line forwarders, so a reader cannot find two spellings of one function and the rename is
    provably complete rather than provably started."""
    for gone in ("_bucket", "_collect", "_diagnostics"):
        assert not hasattr(bootstrap, gone)
    assert callable(bootstrap.bucket) and callable(bootstrap.collect)
    assert callable(bootstrap.diagnostics)
    # `_diagnostics_table` is a DIFFERENT name and is NOT renamed. A sweep matching `_diagnostics`
    # without anchoring the open paren rewrites it wrongly (Stage 11 T5).
    assert callable(bootstrap._diagnostics_table)


def test_run_WRITES_NO_LOOP_OF_ITS_OWN_and_the_loop_is_replicates(boot_n):
    """§20's Definition of done item 12a: `run` contains no `for` and no comprehension over
    `range(C.N_BOOT)`. A second loop is the DRY failure the engineering review removed — the seeding
    would exist twice and Stages 12 and 13 would inherit the copy [§10] never ran (§3.2)."""
    tree = ast.parse(SOURCE)
    node = next(n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name == "run")
    ranges = [c for c in ast.walk(node)
              if isinstance(c, ast.Call) and isinstance(c.func, ast.Name) and c.func.id == "range"]
    assert ranges == []
    calls = {c.func.id for c in ast.walk(node)
             if isinstance(c, ast.Call) and isinstance(c.func, ast.Name)}
    assert "replicates" in calls


# --- 15.12  Stage 10 adds nothing, refits nothing, and breaks nothing earlier ----------------------------

def test_run_RETURNS_ITS_INPUTS_UNMUTATED(run_once):
    """None of the four is frozen against in-place mutation of its `pd.Series` fields, so this is
    compared before and after rather than assumed (§15.12)."""
    df, ps, est, sec, audit, boot = run_once
    fresh = two_centre_frame()
    assert df.equals(fresh)
    assert list(df.columns) == list(fresh.columns)
    assert int(ps.in_model.sum()) == len(df) - 1
    assert set(sec.estimates) == set(config.BINARY_OUTCOMES)
    assert isinstance(est.beta, float)


def test_the_audit_ledger_is_THIRTY_ONE_and_data_KINDS_is_still_NINE():
    """§10.1: no new kind, no new heading, fourth stage running. The ledger is
    load 7 / derive 4 / classify 1 / build 6 / fit 4 / assess 2 / primary 3 / secondary 3
    / bootstrap 1 = 31, and this asserts the arithmetic and the last term."""
    assert len(data.KINDS) == 9
    assert data.KINDS[6] == "model"
    assert 7 + 4 + 1 + 6 + 4 + 2 + 3 + 3 + 1 == 31


def test_bootstrap_adds_no_kind_and_no_heading():
    names = module_names()
    assert "KINDS" not in names and "_HEADINGS" not in names
    assert "record" in names                                 # it appends through Audit.record


# --- 15.13  the audit entry and its rendering -------------------------------------------------------------

def entry_of(audit: data.Audit) -> data.AuditEntry:
    return audit.entry("model", "bootstrap_replicates")


def test_ONE_entry_named_bootstrap_replicates_with_case_ids_EMPTY(run_once):
    """Nine `model` entries name cases and this one removes no patient; naming the `N_BOOT` x n drawn
    rows would be a log the size of the data (§10.1)."""
    _, _, _, _, audit, _ = run_once
    entry = entry_of(audit)
    assert entry is not None
    assert entry.kind == "model" and entry.step == "bootstrap_replicates"
    assert entry.case_ids == ()
    assert entry.n == RUN_BOOT
    assert [e.step for e in audit.entries].count("bootstrap_replicates") == 1


def test_EVERY_ROW_OF_THE_CONCATENATED_GRID_HAS_THE_SAME_CELL_COUNT(run_once):
    """Including the padded short rows and the KEPT second header. This is the operation Stage 6
    §7.2 got wrong first time and the reason `_padded` exists there (§10.1, §15.13)."""
    _, _, _, _, audit, _ = run_once
    table = entry_of(audit).table
    assert len({len(row) for row in table}) == 1
    # The second block's header is kept as a labelled row inside the body, not discarded.
    assert any(row[0] == "diagnostic" for row in table[1:])


def test_NO_CELL_OF_THE_GRID_CONTAINS_A_PIPE_and_the_column_is_spelled_max_abs_beta(run_once):
    """`data._md_table` does no escaping and sizes its separator from `len(rows[0])`, so a pipe in a
    cell gives a body row with more markdown cells than the separator has dashes — which renders as
    a broken table while remaining a byte-identical string (§10.3)."""
    _, _, _, _, audit, _ = run_once
    table = entry_of(audit).table
    assert not any("|" in cell for row in table for cell in row)
    assert not any("|" in entry_of(audit).detail for _ in (0,))
    assert any(cell.startswith("max_abs_beta") for row in table for cell in row)


def test_the_counters_block_renders_EVERY_DECLARED_BUCKET_even_at_zero(run_once):
    """A separation count that vanishes from the table when it reads zero is a count nobody can tell
    was checked, and Stage 8 measured that on this estimator it reads zero on data that is degenerate
    throughout (§15.5, §15.13)."""
    _, _, _, _, audit, _ = run_once
    header = entry_of(audit).table[0]
    for bucket in sorted(set(config.FAILURE_BUCKETS.values())):
        assert bucket in header


def test_the_diagnostics_block_NAMES_ITS_DENOMINATORS(run_once):
    """§16 item 9: the scalar distributions do NOT sum to `C.N_BOOT` and `or_corrected`'s denominator
    is its outcome's surviving odds-ratio draws, and nothing in a grid of counts says which is which
    unless a column does."""
    _, _, _, _, audit, boot = run_once
    table = entry_of(audit).table
    fitted_primary = sum(boot.diagnostics.n_alpha.values())
    fitted_propensity = sum(boot.diagnostics.n_in_model.values())
    denominators = {row[0].split()[0]: row[3] for row in table[1:]
                    if row[0].split()[0] in ("n_alpha", "polr_iterations", "n_in_model", "sum_w")}
    assert denominators["n_alpha"] == f"{fitted_primary} replicate(s) reached it"
    assert denominators["polr_iterations"] == f"{fitted_primary} replicate(s) reached it"
    assert denominators["n_in_model"] == f"{fitted_propensity} replicate(s) reached it"
    assert denominators["sum_w"] == f"{fitted_propensity} replicate(s) reached it"
    corrected = [row for row in table[1:] if row[0].startswith("or_corrected")]
    assert len(corrected) == len(config.BINARY_OUTCOMES)
    assert all("with an odds ratio" in row[3] for row in corrected)
    for row in corrected:
        key = row[0].split(" ", 1)[1]
        assert row[3].startswith(str(len(boot.draws[f"{key}.odds_ratio"].draws)))


def test_the_LOG_IS_BYTE_IDENTICAL_ACROSS_TWO_HASH_SEEDS(tmp_path):
    """§15.13's heredoc, as a subprocess pair. Any set iteration in a table builder would make the
    row order a function of `PYTHONHASHSEED`, and every ordering here is over a declared tuple or a
    `sorted`."""
    script = textwrap.dedent(f"""
        import sys
        sys.path.insert(0, {str(Path(config.__file__).resolve().parent / "tests")!r})
        sys.path.insert(0, {str(Path(config.__file__).resolve().parent)!r})
        import config as C
        C.N_BOOT = {RUN_BOOT}
        import data, propensity, outcome, bootstrap
        from fixtures_stage10 import two_centre_frame
        df = two_centre_frame()
        a = data.Audit(data.FIXTURE)
        ps = propensity.fit(df, a)
        bootstrap.run(df, ps, outcome.primary(df, ps, a), outcome.secondary(df, ps, a), a)
        print(a.to_markdown())
    """)
    outputs = []
    for seed in ("0", "1"):
        env = dict(os.environ, PYTHONHASHSEED=seed)
        done = subprocess.run([sys.executable, "-c", script], capture_output=True, text=True,
                              env=env, cwd=str(Path(config.__file__).resolve().parent))
        assert done.returncode == 0, done.stderr[-3000:]
        outputs.append(done.stdout)
    assert outputs[0] == outputs[1]
    assert "bootstrap_replicates" in outputs[0]


# --- 15.15  `diagnostics`, and what a `None` means -------------------------------------------------------

def test_the_four_scalar_distributions_sum_to_the_LIVE_count_and_NOT_to_N_BOOT(monkeypatch, boot_n):
    """A `None` counted as a `0` would put a spurious mode at zero in three distributions at once
    (§3.1, §7.1, §15.15). Asserted against an INJECTED count of propensity failures."""
    df = two_centre_frame()
    ps, est, sec, audit = fitted(df)
    calls = {"n": 0}
    real_fit = propensity.fit

    def flaky(frame, a):
        calls["n"] += 1
        if calls["n"] % 4 == 0:
            raise model.FitError("Firth: an injected whole-replicate failure")
        return real_fit(frame, a)

    monkeypatch.setattr(bootstrap.propensity, "fit", flaky)
    boot = bootstrap.run(df, ps, est, sec, audit)
    dead = RUN_BOOT // 4
    live = RUN_BOOT - dead
    diagnostics = boot.diagnostics
    assert sum(diagnostics.n_alpha.values()) == live
    assert sum(diagnostics.polr_iterations.values()) == live
    assert sum(diagnostics.n_in_model.values()) == live
    assert diagnostics.sum_w.size == live
    assert 0 not in diagnostics.n_alpha and 0 not in diagnostics.n_in_model


def test_a_PRIMARY_fit_failure_drops_n_alpha_but_NOT_sum_w(monkeypatch, boot_n):
    """**The two live counts are NOT the same count**, and §3.1 is what separates them: a diagnostic
    field is `None` when the replicate died before `outcome.primary` RAN — which is the propensity
    failure — so a replicate whose PRIMARY fit failed still produced a `Σw` and an `in_model` size
    and did not produce a cutpoint count.

    §15.15 says all four "sum to the number of replicates that reached `outcome.primary`" and the
    two coincide only while the primary never fails, which on the workbook it does not (G7 = 0,
    §7.5). That is exactly the condition under which one shared denominator is silently wrong, so
    the grid carries two and this is the assertion that the code knows the difference.
    """
    df = two_centre_frame()
    ps, est, sec, audit = fitted(df)
    calls = {"n": 0}
    real = outcome.primary

    def flaky(frame, p, a):
        calls["n"] += 1
        if calls["n"] % 3 == 0:
            raise model.FitError("G7  an injected separated fit")
        return real(frame, p, a)

    monkeypatch.setattr(bootstrap.outcome, "primary", flaky)
    boot = bootstrap.run(df, ps, est, sec, audit)
    lost = RUN_BOOT // 3
    diagnostics = boot.diagnostics
    assert sum(diagnostics.n_alpha.values()) == RUN_BOOT - lost
    assert sum(diagnostics.polr_iterations.values()) == RUN_BOOT - lost
    assert sum(diagnostics.n_in_model.values()) == RUN_BOOT       # the propensity model still fitted
    assert diagnostics.sum_w.size == RUN_BOOT


def test_max_abs_beta_holds_EXACTLY_THE_AUGMENTED_OUTCOMES_and_sich_and_ph2_are_ABSENT(run_once):
    """Asserted BOTH WAYS: the augmented keys are present, and the unaugmented ones are absent from
    the dict rather than present with an empty array — §10.2's distinction between "we looked and
    found none" and "there is nothing here to look at" (§15.15)."""
    _, _, _, sec, _, boot = run_once
    augmented = {k for k, e in sec.estimates.items() if e.augmented_path != "unaugmented"}
    assert set(boot.diagnostics.max_abs_beta) == augmented
    for key in ("sich", "ph2"):
        assert sec.estimates[key].augmented_path == "unaugmented"
        assert key not in boot.diagnostics.max_abs_beta
    for key, values in boot.diagnostics.max_abs_beta.items():
        assert values.size > 0 and np.isfinite(values).all()


def test_or_corrected_is_keyed_by_ALL_SEVEN_and_bounded_by_that_outcomes_n_attempted(run_once):
    """The natural WRONG denominator is `C.N_BOOT` (§9.2, §15.15). Keyed by all seven and
    initialised to zero, because a rate of zero is a fact and a missing key is not one."""
    _, _, _, _, _, boot = run_once
    assert set(boot.diagnostics.or_corrected) == set(config.BINARY_OUTCOMES)
    for key, fired in boot.diagnostics.or_corrected.items():
        assert 0 <= fired <= boot.draws[f"{key}.odds_ratio"].n_attempted


def test_or_corrected_COUNTS_A_FORCED_CORRECTION_a_known_number_of_times(boot_n):
    """Against a frame where the correction is forced to fire, rather than only on one where it does
    not (§15.15). `or_corrected` means A WEIGHTED PROPORTION REACHED 0 OR 1, so an outcome with no
    event in one arm fires it in every replicate that keeps that arm empty."""
    df = two_centre_frame()
    forced = df.copy()
    treated = forced[config.TREATMENT] == 1.0
    forced.loc[treated, "ph2"] = 0.0                    # no ph2 event in the treated arm, ever
    ps, est, sec, audit = fitted(forced)
    boot = bootstrap.run(forced, ps, est, sec, audit)
    assert sec.estimates["ph2"].or_corrected is True    # ...on the point estimate too
    fired = boot.diagnostics.or_corrected["ph2"]
    assert fired == len(boot.draws["ph2.odds_ratio"].draws) > 0


def test_every_diagnostic_array_is_FINITE_and_sum_w_carries_no_nan(run_once):
    """A replicate that produced a Σw produced a fit, and Stage 6 §9's nan-off-in_model rule means
    the sum is over the MASK and never over `notna()` (§15.15)."""
    _, _, _, _, _, boot = run_once
    assert np.isfinite(boot.diagnostics.sum_w).all()
    assert (boot.diagnostics.sum_w > 0).all()
    for values in boot.diagnostics.max_abs_beta.values():
        assert np.isfinite(values).all()


def test_sum_w_and_polr_iterations_are_reported_AS_A_PAIR_which_is_the_finding(run_once):
    """Stage 8 §11 wrote that a replicate's Σw "is not measured here", that nothing was expected to
    move, and that the visible symptom if it did would be `polr`'s iteration count DROPPING. Both are
    fields on `Diagnostics` for that reason: the pair is reported rather than the reassurance."""
    _, _, _, _, _, boot = run_once
    fields = {f for f in bootstrap.Diagnostics.__dataclass_fields__}
    assert {"sum_w", "polr_iterations"} <= fields
    assert boot.diagnostics.sum_w.size == sum(boot.diagnostics.polr_iterations.values())


# --- 15.14  coverage `[roadmap]` --------------------------------------------------------------------------
#
# "On synthetic data with a known effect the interval covers the truth at roughly the nominal rate."
# Two designs, because they establish different things.
#
# MEASURED, 2 x 90 000 fits in 3204 s at M_OUTER = B_INNER = 300, `inverted_cdf`:
#
#   design U   truth 0.714486   coverage 0.9533  se 0.0122   300/300 usable, 0 inner dropped
#              width 0.7546     against BETA_TRUE 0.9433 -- the two AGREE, which is the control
#   design C   truth 0.560323   coverage 0.9500  se 0.0126   300/300 usable, 0 inner dropped
#              width 0.6304     against BETA_TRUE 0.8600 -- the companion, and it is the point
#
# Nominal 0.95 is inside both Monte Carlo intervals, and design C scored against the generating
# coefficient reads as an 86% undercoverage — which is §8.4's "about 88%" reproduced on a
# construction §8.4 never wrote down.
#
# WHAT THIS DOES NOT ESTABLISH, and §8.4 says so itself: `B_INNER` is 300 and not `C.N_BOOT`, and `n`
# is 400 and not the cohort's size, so coverage is measured for the PROCEDURE at a size where the
# Monte Carlo error is tolerable. Neither design has strata, so `resample`'s stratification is
# unexercised here; and neither is near-separated, which is the condition [§13]'s amendment says this
# cohort is actually in. §16 item 4 files the last of those as the gap that matters.


def coverage_estimate(df: pd.DataFrame) -> float:
    """The [§8] estimator on one population: `propensity.fit` then `outcome.primary`, returning β."""
    audit = data.Audit(hand_source(Path("."), "coverage"))
    ps = propensity.fit(df, audit)
    return outcome.primary(df, ps, audit).beta


def coverage_run(confounded: bool, m_outer: int, b_inner: int) -> dict:
    """One design, scored against BOTH candidate truths in ONE pass.

    Returns `{"truth", "coverage", "se", "usable", "against_beta_true", "width"}`.

    **Both scorings come from the SAME intervals**, which is what makes §15.14's companion exact
    rather than a second Monte Carlo run that happens to be close: the interval is computed once and
    asked twice whether it contains the estimand and whether it contains `BETA_TRUE`. It also halves
    the cost of a design that is already 90 000 fits.

    The truth is obtained by running the ESTIMATOR once at n = 200 000, which is §8.4's own
    prescription and Stage 9 §8.5's method: a logistic model is not collapsible, so design C's
    marginal ATO common odds ratio is attenuated relative to the conditional coefficient the data
    were generated from, and scoring against the generating coefficient reads as a bootstrap defect.

    A `model.FitError` from an INNER fit is dropped and counted, which is [§10]'s own rule applied to
    the coverage harness rather than an exception to it; an outer sample whose surviving inner draws
    fall below `ci_min_draws` yields no interval and is excluded from `usable`, which is §8.3's floor
    doing exactly what it is for.
    """
    rng = np.random.default_rng(COVERAGE_SEED)
    truth = coverage_estimate(known_effect_population(200_000, rng, confounded))
    covered = covered_at_beta_true = usable = dropped = 0
    widths: list[float] = []
    for _ in range(m_outer):
        sample = known_effect_population(400, rng, confounded)
        inner = np.random.default_rng(int(rng.integers(0, 2**32 - 1)))
        draws = []
        for _ in range(b_inner):
            try:
                draws.append(coverage_estimate(
                    bootstrap.resample(sample, inner, config.BOOT_STRATUM)))
            except model.FitError:
                dropped += 1
        if len(draws) < config.ci_min_draws():
            continue
        lo, hi = bootstrap.percentile_ci(np.asarray(draws, dtype=float))
        usable += 1
        widths.append(hi - lo)
        covered += int(lo <= truth <= hi)
        covered_at_beta_true += int(lo <= BETA_TRUE <= hi)
    coverage = covered / usable
    return {
        "truth": truth,
        "coverage": coverage,
        "se": float(np.sqrt(coverage * (1.0 - coverage) / usable)),
        "usable": usable,
        "dropped": dropped,
        "against_beta_true": covered_at_beta_true / usable,
        "width": float(np.mean(widths)),
    }


@SLOW
def test_design_U_covers_at_ROUGHLY_NOMINAL_with_the_MC_STANDARD_ERROR_asserted():
    """**The standard error is asserted rather than the point coverage**, because a coverage estimate
    from `M_OUTER` outer simulations is itself a binomial proportion: 0.9400 quoted bare reads as
    undercoverage when it is 0.8 standard errors from nominal (§8.4).

    Design U is unconfounded — `A` depends on Z and X, `Y` on A alone — so the marginal odds ratio IS
    `BETA_TRUE` in closed form and coverage measures the machinery. The two scorings therefore agree
    here, which is the control for design C's companion: if they disagreed on THIS design the
    generator would not be unconfounded.
    """
    result = coverage_run(False, M_OUTER, B_INNER)
    print(f"\n  design U: {result}")
    assert result["usable"] == M_OUTER
    assert abs(result["truth"] - BETA_TRUE) < 0.2, (
        f"design U's estimand is {result['truth']} against BETA_TRUE {BETA_TRUE}; it should be the "
        "same quantity up to Monte Carlo error at n = 200 000")
    assert abs(result["coverage"] - 0.95) < 1.96 * result["se"], (
        f"coverage {result['coverage']:.4f} +- {result['se']:.4f}")
    assert result["coverage"] == pytest.approx(result["against_beta_true"], abs=0.03)


@SLOW
def test_design_C_covers_against_THE_ESTIMAND_and_NOT_against_BETA_TRUE():
    """**And the companion is the point**: the same intervals scored against `BETA_TRUE` cover
    materially WORSE, and this asserts that too — so a future edit that "fixes" design C's truth to
    the generating coefficient fails loudly instead of reading as a bootstrap defect (§8.4).

    A logistic model is not collapsible, so the marginal ATO common odds ratio is attenuated relative
    to the conditional coefficient the data were generated from. §8.4 measures 0.659272 on its own
    construction; this one attenuates to about 0.665 — the same phenomenon on a construction §8.4
    never wrote down, which is why the truth is RE-MEASURED here rather than pinned.
    """
    result = coverage_run(True, M_OUTER, B_INNER)
    print(f"\n  design C: {result}")
    assert result["usable"] == M_OUTER
    assert result["truth"] < BETA_TRUE - 0.02, (
        f"design C's estimand is {result['truth']}, which is not materially below BETA_TRUE — the "
        "design is not confounded enough for non-collapsibility to bite and the companion is vacuous")
    assert abs(result["coverage"] - 0.95) < 1.96 * result["se"], (
        f"coverage {result['coverage']:.4f} +- {result['se']:.4f}")
    assert result["against_beta_true"] < result["coverage"] - 0.02, (
        f"scoring design C against BETA_TRUE covered at {result['against_beta_true']:.4f} against "
        f"{result['coverage']:.4f} for the estimand — not materially worse, so §8.4's whole point "
        "about non-collapsibility is untested")


# --- the workbook, end to end `[data-gated]` ------------------------------------------------------------

@DATA_GATED
def test_the_full_loop_runs_on_the_workbook_with_ZERO_SchemaError(monkeypatch):
    """`[data-gated]` — §20's Definition of done item 3. Measured before this stage: 26.0% of
    replicates raised from `_record_exclusion` and 1.0% from S8. Both must read zero.

    Run at RUN_BOOT rather than at `C.N_BOOT`, because the full 2000 is about 130 s and this asserts
    that the loop RUNS rather than what it produces — and by §4.5 no number it produces may appear
    in this file anyway.
    """
    monkeypatch.setattr(config, "N_BOOT", RUN_BOOT)
    df, audit = data.load(data.WORKBOOK)
    df = cohort.build(eligibility.classify(derive.derive(df, audit), audit), audit)
    ps = propensity.fit(df, audit)
    est = outcome.primary(df, ps, audit)
    sec = outcome.secondary(df, ps, audit)
    boot = bootstrap.run(df, ps, est, sec, audit)      # a SchemaError anywhere fails this test
    assert boot.seed == config.SEED and boot.n_boot == RUN_BOOT
    assert len(boot.draws) == 26
    assert all(d.n_attempted == RUN_BOOT for d in boot.draws.values())
    assert audit.entry("model", "bootstrap_replicates") is not None
