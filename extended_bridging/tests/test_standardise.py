"""Acceptance tests for Stage 12 — §20 of `specs/stage12_all_centre_standardisation.md`.

Section banners match the specification's numbers, as every test module since `test_cohort.py` does.
§20.2's config assertions live in `test_config.py` and §20.5's and §20.12's estimator assertions in
`test_model.py`, because that is where the changed code is.

Tests needing the private workbook are marked `DATA_GATED` and tagged `[data-gated]` in their
docstrings. Everything else runs on a plain checkout with no `data/`.

**NO STANDARDISED PROBABILITY, RISK DIFFERENCE, INTERVAL LIMIT OR CENTRE RANDOM INTERCEPT FROM THE
WORKBOOK APPEARS IN THIS FILE.** Stage 8 §4.3 set the rule, Stage 10 §4.5 tightened it to interval
limits and p-values, Stage 11 §1 to E-values, and Stage 12 §4.3 tightens it once more to exactly the
two quantities this stage produces. A standardised probability is worse than an interval limit for the
reason an E-value is: it reads as a directly clinical number — *"this many patients in 100 would be
independent under bridging"* — and a reader who meets it in a version-controlled document will quote
it. A centre random intercept is worse still, because it is a centre-level outcome contrast on four
named hospitals and this analysis is not powered to make one.

So the workbook-gated tests assert COUNTS, SET MEMBERSHIP, IDENTITIES AND PROPERTIES — the ledger is
104 over four centres, the distribution sums to 1, `mrs_0_2` is exactly `rd[2]`, everyone outside the
box is a comparator patient, the two box readings agree — and every numerical regression pin is
measured on `fixtures_stage12.py`.

**The three silent failures this stage exists to make loud**, each tested with a companion showing
what happens without the guard::

    a counterfactual design      §20.5   without the overwrite rule: `design` drops the constant
    rebuilt by `model.design`            exposure, `polr` fits nine covariates, and both regimes
                                         are the same regime. Nothing raises

    a re-expansion by position   §20.5   without it: a six-column distribution broadcasts against
    on a collapsed replicate             a seven-level average and misaligns every level above 5

    a mortality interval taken   §20.7   without its own draws: the limit moves by up to 1.7e-3 on
    by reflecting `rd_5`'s               the risk-difference scale, which is 0.17 percentage points
"""
from __future__ import annotations

import ast
import dataclasses
import os
from pathlib import Path
from typing import Final

import numpy as np
import pandas as pd
import pytest

import balance
import bootstrap
import cohort
import config
import data
import derive
import eligibility
import model
import standardise
from fixtures_stage12 import (COLLAPSE_LEVEL, collapsed_level_frame, four_centre_frame,
                              no_exposure_frame)
from test_data import hand_source

DATA_GATED = pytest.mark.skipif(
    not config.DATA_XLSX.exists(),
    reason="the private workbook is gitignored and absent from this checkout")

# §20.11's first bullet runs the Stage 1-10 pipeline TWICE at C.N_BOOT = 2000, which is about six
# minutes. It is gated on an environment variable and skips LOUDLY, which is Stage 10 §15.14's
# mechanism and is what `pytest -rs` is for: the default state of a slow measurement is "not run",
# and a check matching nothing must not pass as green.
SLOW = pytest.mark.skipif(
    os.environ.get("STAGE12_SLOW") != "1",
    reason="§20.11's byte-identity guard is two full 2000-replicate pipeline runs and about six "
           "minutes. Run with STAGE12_SLOW=1 to open the gate.")

MODULE = Path(standardise.__file__).resolve()
SOURCE = MODULE.read_text(encoding="utf-8")

# The replicate count every `inference`-driving test uses. `inference` takes no `n` — C.N_BOOT is
# prespecified and [§10] refits in every replicate, so it is not a runtime knob — and 2000 replicates
# of the fixture is about four minutes per call. 60 is comfortably above `ci_min_draws(0.95)` = 40,
# so the below-floor branch stays unreachable except where a test reaches it deliberately.
RUN_BOOT: Final[int] = 60

# The step names `cohort.build` owns and this stage may not use, measured against the landed module
# (§4.2). Written out because §20.9 asserts the INTERSECTION is empty rather than asserting a list:
# `Audit.entry` is first-match, so a name reused between the two makes `audit.entry("cohort", ...)`
# return whichever ran first, and the two populations differ by 11 records and one whole centre.
COHORT_STEPS: Final[frozenset[str]] = frozenset({
    "restrict_centres", "restrict_eligibility", "cohort_flow", "absence_by_cohort_column",
    "core_above_median", "constant_covariates"})


def message_of(excinfo) -> str:
    return str(excinfo.value)


def audit_of() -> data.Audit:
    return data.Audit(hand_source(Path("."), "stage12"))


def classified():
    """The classified frame and its Audit, over the private workbook. [data-gated]"""
    df, audit = data.load()
    df = derive.derive(df, audit)
    return eligibility.classify(df, audit), audit


def driver(df: pd.DataFrame, audit: data.Audit, with_inference: bool = False):
    """§19's canonical call order — the ONE place it is written, and §20.9 asserts against THIS.

    `support` must precede the restricted `all_centre`, whose `over=` is its output, which is why
    entry 7 falls between 6 and 8. **A driver that omitted the restricted `all_centre` would record
    nine entries and not ten**, which is the failure §20.9's count catches and the reason §19 spells
    the call out rather than leaving sensitivity 1's point estimate implied.
    """
    pop = standardise.population(df, audit)
    std = standardise.all_centre(pop, audit)
    sup = standardise.support(pop, audit)
    ssup = standardise.all_centre(pop, audit, over=sup.inside)
    hier = standardise.hierarchical(pop, audit)
    boot = standardise.inference(pop, audit) if with_inference else None
    return pop, std, sup, ssup, hier, boot


# --- 20.1  the population is [§14a]'s and not the cohort's ----------------------------------------

@DATA_GATED
def test_the_population_is_104_over_FOUR_centres_and_the_cohort_is_93_over_THREE():
    """§20.1. Both in one test, so the DIFFERENCE is the assertion. [data-gated]

    [§3] restriction 1 removes the never-IVT centre, which is precisely the population [§14a] exists
    to include, and `cohort.build` offers no seam to skip it. So the two populations are different
    objects and this is the one place that is asserted rather than described.
    """
    df, audit = classified()
    pop = standardise.population(df, audit)
    built = cohort.build(df, audit_of())

    assert len(pop) == 104
    assert sorted(pop["center"].unique()) == sorted(config.CENTER_ORDER)
    assert len(built) == 93
    assert set(built["center"].unique()) == set(cohort.treating_centres(df))
    assert len(pop) - len(built) == 11
    assert len(set(pop["center"].unique()) - set(built["center"].unique())) == 1


@DATA_GATED
def test_every_never_IVT_record_in_the_population_is_UNTREATED():
    """§20.1. No treated patient exists at a never-IVT centre BY DEFINITION. [data-gated]

    A workbook in which one does must fail T9 before this line is reached, which is why this is an
    assertion about the population and not a filter applied to it.
    """
    df, audit = classified()
    pop = standardise.population(df, audit)
    never = tuple(c for c in config.CENTER_ORDER if c not in cohort.treating_centres(pop))
    at_never = pop[pop["center"].isin(never)]
    assert len(at_never) > 0
    assert (at_never[config.TREATMENT] == min(config.TREATMENT_LABELS)).all()


@DATA_GATED
def test_the_three_11_losses_are_asserted_BY_CAUSE_and_both_outcome_misses_are_at_the_same_centre():
    """§20.1, §4.1. One covariate-incomplete, two outcome-missing, both at the never-IVT centre.

    The two causes fall on different centres and that is a fact about this workbook worth a reader's
    attention: the centre that contributes no contrast is also the one whose outcome ascertainment is
    incomplete, so [§14a]'s target population loses 2 of the 14 patients its whole reason for
    existing is to include. [data-gated]
    """
    df, audit = classified()
    retained = df[eligibility.retained(df)]
    complete = model.complete_cases(retained, config.STANDARDISATION_COVARIATES)
    observed = retained[config.PRIMARY_OUTCOME].notna()

    incomplete = retained[~complete]
    assert len(incomplete) == 1
    # The one incomplete record is incomplete on BOTH volume covariates, not on one of them.
    assert [c for c in config.STANDARDISATION_COVARIATES if incomplete[c].isna().any()] == [
        "core_ml", "tmax6_ml"]

    missing = retained[complete & ~observed]
    assert len(missing) == 2
    never = tuple(c for c in config.CENTER_ORDER if c not in cohort.treating_centres(df))
    assert set(missing["center"]) == set(never)


def test_population_raises_T6_on_a_frame_classify_has_NOT_run_on():
    """§20.14, T6. The `eligibility` column is the restriction; without it there is no population."""
    frame = four_centre_frame().drop(columns=[config.ELIGIBILITY])
    with pytest.raises(config.SchemaError) as excinfo:
        standardise.population(frame, audit_of())
    assert message_of(excinfo).startswith("T6")
    assert config.ELIGIBILITY in message_of(excinfo)


def test_population_raises_T7_on_an_empty_population_and_again_on_an_empty_ARM():
    """§20.14, T7. Both cases, and the emptiness is the POPULATION's rather than the input frame's.

    The second case is the one that matters: a frame with 108 records, every one of them eligible,
    every one of them a comparator patient. The standardisation duplicates every patient under both
    regimes and averages, so an empty arm leaves every average computable while the treatment
    coefficient is unidentified — the arithmetic returns a risk difference for a contrast the data
    cannot support.
    """
    empty = four_centre_frame()
    empty[config.ELIGIBILITY] = config.INELIGIBLE
    with pytest.raises(config.SchemaError) as excinfo:
        standardise.population(empty, audit_of())
    assert message_of(excinfo).startswith("T7")

    one_arm = four_centre_frame()
    one_arm[config.TREATMENT] = float(min(config.TREATMENT_LABELS))
    with pytest.raises(config.SchemaError) as excinfo:
        standardise.population(one_arm, audit_of())
    assert message_of(excinfo).startswith("T7")
    assert config.TREATMENT_LABELS[max(config.TREATMENT_LABELS)] in message_of(excinfo)


def test_derive_cohort_is_NOT_on_this_stages_path_asserted_by_SCANNING_the_module():
    """§20.1. The Stage 1 §7 scan's pattern, because the failure it prevents is SILENT.

    `derive_cohort` computes `core_above_median` on the frame it is given, and Stage 5 §4.5
    established that the median is the COHORT's — 5.0 mL over 93 records against 6.0 mL over 126. A
    [§14a] population that went through it would carry a subgroup variable cut at the wrong place,
    and [§14a] names no subgroup at all.
    """
    tree = ast.parse(SOURCE)
    imported = {n.names[0].name for n in ast.walk(tree) if isinstance(n, ast.Import)}
    assert "derive" not in imported
    # And the module Stage 12 must not touch at all, from the other side (§0.2).
    assert "propensity" not in imported
    # The NAME appears once, in `population`'s docstring, saying why it is not called. What must be
    # absent is a CALL, which is the thing that would compute the median on the wrong population.
    calls = {ast.unparse(c.func) for c in ast.walk(tree) if isinstance(c, ast.Call)}
    assert not [c for c in calls if "derive" in c]
    assert not [c for c in calls if "propensity" in c]


def test_the_population_applies_restriction_2_through_eligibility_retained_and_NOT_a_literal():
    """§4.1. The one place the `== eligible` / `!= ineligible` difference is decided is Stage 4's.

    43 of 126 records turn on it — 80% of the control arm in the primary cohort — so two stages
    resolving it differently would not look like a bug in either of them.
    """
    tree = ast.parse(SOURCE)
    node = next(n for n in tree.body
                if isinstance(n, ast.FunctionDef) and n.name == "population")
    calls = {ast.unparse(c.func) for c in ast.walk(node) if isinstance(c, ast.Call)}
    assert "eligibility.retained" in calls
    assert "model.complete_cases" in calls
    # And the retained set is never spelled as a comparison against a declared class: the one place
    # the `== eligible` / `!= ineligible` difference is decided is Stage 4's, and 43 of 126 records
    # turn on it.
    body = "\n".join(line for line in ast.unparse(node).splitlines()
                     if not line.strip().startswith("#"))
    for forbidden in ("C.ELIGIBLE", "C.INELIGIBLE", "C.ELIGIBILITY_RETAINED"):
        assert forbidden not in body.split('"""')[-1], forbidden


# --- 20.3  the estimand keys are derived and not counted ------------------------------------------

def test_the_key_count_is_DERIVED_from_the_config_and_never_written_as_69():
    """§20.3. Three standardised blocks, one `beta` per fit and there are two, `hier.sigma`.

    Stage 8 §12 made the level set derived from the plausible range so the two cannot disagree; this
    key set is the next link in that chain.
    """
    keys = standardise._estimand_keys()
    expected = 3 * (len(config.MRS_THRESHOLDS) + 2 + 2 * len(config.MRS_LEVELS)) + 2 + 1
    assert len(keys) == expected
    assert len(set(keys)) == len(keys)
    node = next(n for n in ast.parse(SOURCE).body
                if isinstance(n, ast.FunctionDef) and n.name == "_estimand_keys")
    # The docstring carries §3.3's arithmetic and names the total; the CODE must not.
    code = ast.unparse(ast.Module(
        body=[n for n in node.body if not (isinstance(n, ast.Expr)
                                           and isinstance(n.value, ast.Constant)
                                           and isinstance(n.value.value, str))],
        type_ignores=[]))
    for literal in ("69", "22", "14", "6"):
        assert literal not in code, literal
    assert "C.MRS_LEVELS" in code and "C.MRS_THRESHOLDS" in code


def test_support_beta_IS_IN_NO_KEY_SET_and_the_absence_is_14a_s():
    """§20.3, §3.3, §11. Asserted BY NAME, in the manner §8's absent `odds_ratio` field is.

    [§14a]'s sensitivity restricts THE STANDARDISATION POPULATION, which §11 implements by restricting
    the averaging set and not the fit. So the support arm has no `beta` of its own: a `support.beta`
    would be the pooled `beta` — bit-identical, draw for draw, with a bit-identical interval —
    presented as a property of a STANDARDISATION, which is the one reading [§14a] forbids in as many
    words. A later "make the arms uniform" edit fails here rather than putting `beta` into a row
    labelled as a sensitivity result.
    """
    keys = standardise._estimand_keys()
    assert "support.beta" not in keys
    assert "beta" in keys
    assert "hier.beta" in keys
    assert "hier.sigma" in keys
    assert len([k for k in keys if k.endswith("beta")]) == 2
    assert "support.sigma" not in keys


def test_the_three_arms_carry_the_SAME_standardised_block():
    """§3.3. Standardised quantities are per-ARM, so the three blocks are identical in shape."""
    keys = standardise._estimand_keys()
    blocks = {}
    for prefix in ("", "support.", "hier."):
        blocks[prefix] = {k[len(prefix):] for k in keys
                          if k.startswith(prefix) and k[len(prefix):] not in ("beta", "sigma")
                          and "." not in k[len(prefix):]}
    assert blocks[""] == blocks["support."] == blocks["hier."]
    assert len(blocks[""]) == len(config.MRS_THRESHOLDS) + 2 + 2 * len(config.MRS_LEVELS)


# --- 20.4  the fit is unweighted, and the CALL SITE is the assertion ------------------------------

def test_model_polr_is_called_with_TWO_POSITIONAL_ARGUMENTS_and_no_w():
    """§20.4. Asserted by inspecting the CALL and never by comparing numbers.

    `w=None` and `w=ones` are numerically identical and semantically different, and `polr`'s docstring
    says the acceptance suite must turn on the difference being visible in the call: *"a caller who
    passes nothing has said something different from a caller who passes ones"*. [§14a] weights
    nobody — there is no propensity model to weight by — so this stage passes nothing.
    """
    tree = ast.parse(SOURCE)
    calls = [c for c in ast.walk(tree)
             if isinstance(c, ast.Call) and ast.unparse(c.func) == "model.polr"]
    assert len(calls) >= 1
    for call in calls:
        assert len(call.args) == 2, ast.unparse(call)
        assert call.keywords == [], ast.unparse(call)
    # And no weight vector is constructed anywhere in the module under any name.
    assert "np.ones(len(y))" not in SOURCE
    assert "weights=" not in SOURCE


# --- 20.5  the distribution is a distribution, structurally ---------------------------------------

def test_the_distribution_sums_to_one_and_stays_in_range_on_the_FIXTURE():
    """§20.5. Every row of `ordinal_probabilities`, and then the average of them."""
    pop = standardise.population(four_centre_frame(), audit_of())
    std = standardise.all_centre(pop, audit_of())
    for arm, levels in std.distribution.items():
        assert abs(sum(levels.values()) - 1.0) < 1e-12
        assert all(0.0 <= value <= 1.0 for value in levels.values())
        sequence = [std.cumulative[arm][k] for k in config.MRS_THRESHOLDS]
        assert all(a <= b for a, b in zip(sequence, sequence[1:]))


@DATA_GATED
def test_the_distribution_sums_to_one_on_the_WORKBOOK_too():
    """§20.5. The same property, on the population it will be reported over. [data-gated]"""
    df, audit = classified()
    pop = standardise.population(df, audit)
    std = standardise.all_centre(pop, audit)
    for levels in std.distribution.values():
        assert abs(sum(levels.values()) - 1.0) < 1e-12
        assert all(0.0 <= value <= 1.0 for value in levels.values())


def test_the_re_expansion_puts_a_STRUCTURAL_ZERO_at_the_unoccupied_level():
    """§20.5, §6.2, §6.3. Six fitted columns, seven declared levels, and `RD_4 == RD_5` EXACTLY.

    This is the single most consequential measurement in the specification made into a test: mRS 5
    carries 3 patients of 104 on the workbook, so 4.15% of replicates lose it. A six-column
    distribution cannot be averaged with a seven-column one, and a naive implementation would either
    raise deep inside numpy or — worse — broadcast and silently misalign every level above 5.
    """
    frame = collapsed_level_frame()
    assert COLLAPSE_LEVEL not in set(frame[config.PRIMARY_OUTCOME])
    pop = standardise.population(frame, audit_of())
    std = standardise.all_centre(pop, audit_of())

    assert len(std.fit.categories) == len(config.MRS_LEVELS) - 1
    assert COLLAPSE_LEVEL not in std.fit.categories
    assert len(std.fit.alpha) == len(config.MRS_LEVELS) - 2

    for arm, levels in std.distribution.items():
        assert set(levels) == set(config.MRS_LEVELS)                 # SEVEN, not six
        assert levels[COLLAPSE_LEVEL] == 0.0                         # the structural zero
        assert abs(sum(levels.values()) - 1.0) < 1e-12
        sequence = [std.cumulative[arm][k] for k in config.MRS_THRESHOLDS]
        assert all(a <= b for a, b in zip(sequence, sequence[1:]))

    # §6.3's third consequence, exactly and by construction: a structural zero at level 5 makes the
    # cumulative probability at thresholds 4 and 5 the same number in BOTH arms.
    assert std.rd[4] == std.rd[5]


def test_T10_fires_on_a_DELIBERATELY_CORRUPTED_distribution():
    """§20.5, T10. The guard is unreachable on the workbook, so this is how it is known to be live.

    Measured worst `|sum - 1|` is 6.7e-16 across every arm of every replicate, against a bound of
    1e-12 — five orders of margin. That is what makes it safe to specify a `SchemaError` that kills
    the whole bootstrap, and it is also why the check has to be exercised deliberately.
    """
    good = {1: {level: 1.0 / len(config.MRS_LEVELS) for level in config.MRS_LEVELS},
            0: {level: 1.0 / len(config.MRS_LEVELS) for level in config.MRS_LEVELS}}
    cumulative = {arm: {k: float(np.cumsum(list(levels.values()))[j])
                        for j, k in enumerate(config.MRS_THRESHOLDS)}
                  for arm, levels in good.items()}
    standardise._assert_distribution(good, cumulative, "clean")      # does not raise

    broken = {arm: dict(levels) for arm, levels in good.items()}
    broken[1][config.MRS_LEVELS[0]] += 1e-9
    with pytest.raises(config.SchemaError) as excinfo:
        standardise._assert_distribution(broken, cumulative, "corrupted")
    assert message_of(excinfo).startswith("T10")

    descending = {arm: {k: -float(k) for k in config.MRS_THRESHOLDS} for arm in good}
    with pytest.raises(config.SchemaError) as excinfo:
        standardise._assert_distribution(good, descending, "not monotone")
    assert message_of(excinfo).startswith("T10")


def test_the_counterfactual_design_OVERWRITES_a_column_and_never_rebuilds_the_design():
    """§7.1, and the companion showing the silent failure. This is the DRY-run of the whole stage.

    `model.design` drops constant columns, so a frame in which every patient has `A = 1` yields a
    design with NO TREATMENT COLUMN AT ALL. `polr` fits the remaining nine, and the "counterfactual"
    prediction is made under a model with no exposure in it — both regimes become the same regime,
    every risk difference is zero, and NOTHING RAISES.
    """
    pop = standardise.population(four_centre_frame(), audit_of())
    X, _ = model.design(pop, (config.TREATMENT,) + config.STANDARDISATION_COVARIATES)
    assert config.TREATMENT in X.columns

    # The specified route: the column is overwritten and the column SET cannot change.
    for arm in config.TREATMENT_LABELS:
        Xa = standardise._regime_design(X, arm)
        assert tuple(Xa.columns) == tuple(X.columns)
        assert (Xa[config.TREATMENT] == float(arm)).all()

    # The rejected route, and what it produces. `assign` then re-`design` LOSES the exposure.
    rebuilt = pop.copy()
    rebuilt[config.TREATMENT] = float(max(config.TREATMENT_LABELS))
    Xr, dropped = model.design(rebuilt, (config.TREATMENT,) + config.STANDARDISATION_COVARIATES)
    assert config.TREATMENT not in Xr.columns            # gone, silently
    assert config.TREATMENT in dropped                   # `design` names it and does not raise


# --- 20.6  the two identities ---------------------------------------------------------------------

@pytest.mark.parametrize("frame", [four_centre_frame, collapsed_level_frame])
def test_the_two_identities_hold_at_their_OWN_tolerances(frame):
    """§20.6, §7.3. `mrs_0_2` EXACTLY and `mortality` to 1e-15, and the asymmetry is the point.

    `mrs_0_2` IS `rd[2]` — the same cumulative difference read under [§14a]'s name for it, one
    computation and not two — so the comparison is exact and 1e-15 would hide a real defect.
    `mortality` is genuinely a second computation, from `P(Y = 6)` in each arm rather than from
    `1 - P(Y <= 5)`, so it agrees to floating point and not to the bit.
    """
    pop = standardise.population(frame(), audit_of())
    std = standardise.all_centre(pop, audit_of())
    assert abs(std.mrs_0_2 - std.rd[2]) == 0.0
    assert abs(std.mortality + std.rd[5]) < 1e-15


@DATA_GATED
def test_the_two_identities_hold_on_the_WORKBOOK():
    """§20.6. [data-gated]"""
    df, audit = classified()
    std = standardise.all_centre(standardise.population(df, audit), audit)
    assert abs(std.mrs_0_2 - std.rd[2]) == 0.0
    assert abs(std.mortality + std.rd[5]) < 1e-15


def test_mortality_is_computed_from_the_DISTRIBUTIONS_and_not_as_minus_rd_5():
    """§20.6. By monkeypatching `rd` and requiring `mortality` to be unaffected.

    **An implementation that derived one from the other would pass the identity test trivially**,
    which is exactly why this test exists and why it applies to `mortality` and not to `mrs_0_2`:
    `mortality` is the quantity an implementer might "simplify" into `-rd[5]`, and this is what
    prevents it.

    The monkeypatch is on `np.cumsum`, because that is what produces the cumulative table `rd` is a
    difference of. If `mortality` were derived from `rd` it would move with it; being computed from
    `distribution`, it does not.
    """
    pop = standardise.population(four_centre_frame(), audit_of())
    honest = standardise.all_centre(pop, audit_of())

    tree = ast.parse(SOURCE)
    node = next(n for n in ast.walk(tree)
                if isinstance(n, ast.FunctionDef) and n.name == "_standardise")
    body = ast.unparse(node)
    # The expression is a difference of two DISTRIBUTION entries at the top level, not of `rd`.
    assert "mortality = distribution[treated][top] - distribution[control][top]" in body
    assert "mortality = -rd[" not in body
    assert "mortality = -std.rd" not in body

    # And the value is a genuine second computation: it agrees with -rd[5] to floating point but is
    # NOT the same float, which a derived implementation could not produce.
    assert honest.mortality != -honest.rd[5]
    assert abs(honest.mortality + honest.rd[5]) < 1e-15


def test_T11_fires_on_a_mortality_perturbed_by_1e_14():
    """§20.14, T11. Above §20.6's 1e-15 bound and far below anything a reader would see.

    So the check is known to be live AT ITS STATED TOLERANCE rather than merely present.
    """
    rd = {k: 0.01 * k for k in config.MRS_THRESHOLDS}
    standardise._assert_identities(rd, rd[2], -rd[5], "clean")        # does not raise
    with pytest.raises(config.SchemaError) as excinfo:
        standardise._assert_identities(rd, rd[2], -rd[5] + 1e-14, "perturbed")
    assert message_of(excinfo).startswith("T11")
    assert "mortality" in message_of(excinfo)

    # `mrs_0_2`'s comparison is EXACT, so the perturbation that tests it is the smallest one that
    # exists. `rd[2] + 1e-18` is not it: 1e-18 is below the ulp of 0.02 and the float does not move,
    # which would make this test pass against an implementation with no check at all.
    assert rd[2] + 1e-18 == rd[2]
    with pytest.raises(config.SchemaError) as excinfo:
        standardise._assert_identities(rd, np.nextafter(rd[2], 1.0), -rd[5], "perturbed")
    assert message_of(excinfo).startswith("T11")
    assert "mrs_0_2" in message_of(excinfo)


# --- 20.7  the mortality interval is NOT the reflection of RD_5's ---------------------------------

def test_the_percentile_limits_of_minus_X_are_NOT_the_reflected_limits_of_X():
    """§20.7, §7.4. The measurement, as a test, at the boundary that found it.

    `C.PERCENTILE_METHOD` is pinned to `"inverted_cdf"` because it is the same object
    `Pr(beta* <= 0)` is computed from, which is what makes [§10]'s p-value agree with its interval BY
    CONSTRUCTION. It is not a symmetric quantile definition: the lower limit of `-X` is not the
    negated upper limit of `X`, because the order statistic reached from below at level `q` is not the
    mirror of the one reached from below at `1 - q`.

    The test asserts the difference is NON-ZERO rather than asserting a magnitude, and the companion
    asserts they AGREE under `"linear"` — so the test documents that the disagreement is the pin's and
    not an arithmetic error. Stage 10 §15.8 used the same shape for the same pin.
    """
    generator = np.random.default_rng(config.SEED)
    draws = generator.normal(0.02, 0.1, config.N_BOOT if config.N_BOOT < 2001 else 2000)

    lo, hi = bootstrap.percentile_ci(draws)
    lo_neg, hi_neg = bootstrap.percentile_ci(-draws)
    assert (lo_neg, hi_neg) != (-hi, -lo)
    assert abs(lo_neg + hi) > 0.0 or abs(hi_neg + lo) > 0.0

    q_lo, q_hi = 100.0 * (1.0 - config.CI_LEVEL) / 2.0, 100.0 - 100.0 * (1.0 - config.CI_LEVEL) / 2.0
    lo_l, hi_l = np.percentile(draws, [q_lo, q_hi], method="linear")
    lo_nl, hi_nl = np.percentile(-draws, [q_lo, q_hi], method="linear")
    assert abs(lo_nl + hi_l) < 1e-15
    assert abs(hi_nl + lo_l) < 1e-15


def test_mortality_carries_its_OWN_percentile_call_and_is_not_derived_from_rd_5s():
    """§20.7. `mortality` is a first-class key in `_estimand_keys`, so `intervals` computes its own.

    A "simplification" that derived the mortality interval by reflecting `rd_5`'s would shift a limit
    by up to 1.7e-3 on the risk-difference scale — 0.17 percentage points, which a manuscript prints.
    """
    keys = standardise._estimand_keys()
    assert "mortality" in keys and "rd_5" in keys
    tree = ast.parse(SOURCE)
    assert not [c for c in ast.walk(tree) if isinstance(c, ast.Call)
                and ast.unparse(c.func) == "bootstrap.percentile_ci"]
    calls = [c for c in ast.walk(tree) if isinstance(c, ast.Call)
             and ast.unparse(c.func) == "bootstrap.intervals"]
    assert len(calls) == 1
    assert ast.unparse(calls[0].args[1]) == "lambda key: False"


# --- 20.8  the conditional-measure guard ----------------------------------------------------------

def test_every_Standardisation_carries_the_CONDITIONAL_guard_string():
    """§20.8, §8. `measure` is DATA, so Stage 14 prints it from the record rather than from a label."""
    pop = standardise.population(four_centre_frame(), audit_of())
    audit = audit_of()
    std = standardise.all_centre(pop, audit)
    sup = standardise.support(pop, audit)
    ssup = standardise.all_centre(pop, audit, over=sup.inside)
    hier = standardise.hierarchical(pop, audit)
    for record in (std, ssup, hier.standardisation):
        assert record.measure is standardise._CONDITIONAL
        assert "conditional" in record.measure.lower()
        assert "marginal" in record.measure.lower()      # as the thing it is NOT


def test_there_is_NO_odds_ratio_ATTRIBUTE_and_the_name_appears_nowhere_in_the_module():
    """§20.8, §8. Three enforcements, each as an assertion.

    `Primary.odds_ratio` is the [§8] MARGINAL common odds ratio. A Stage 14 formatter written against
    one record and pointed at the other would print a conditional quantity under a marginal label —
    so the enforcement is that the attribute does not exist and the failure is an `AttributeError`.
    """
    pop = standardise.population(four_centre_frame(), audit_of())
    std = standardise.all_centre(pop, audit_of())
    assert not hasattr(std, "odds_ratio")
    assert hasattr(std, "conditional_odds_ratio")
    assert hasattr(std, "conditional_log_odds")
    fields = {f.name for f in dataclasses.fields(standardise.Standardisation)}
    assert "odds_ratio" not in fields

    # The string appears nowhere outside a comment or docstring.
    tree = ast.parse(SOURCE)
    literals = [n.value for n in ast.walk(tree) if isinstance(n, ast.Constant)
                and isinstance(n.value, str)]
    docstrings = set()
    for node in ast.walk(tree):
        if isinstance(node, (ast.Module, ast.FunctionDef, ast.ClassDef)):
            doc = ast.get_docstring(node, clean=False)
            if doc is not None:
                docstrings.add(doc)
    assert not [s for s in literals if "odds_ratio" in s and s not in docstrings]


def test_the_estimand_key_is_beta_and_the_interval_is_on_the_LOG_scale():
    """§8's third enforcement. `exp` of a percentile limit IS the percentile limit of `exp`.

    Stage 11 §7.2 established that under `PERCENTILE_METHOD`, so a Stage 14 reader wanting the ratio
    scale exponentiates the limits and gets the right interval. The key stays on the log scale so the
    record never holds a number a careless read turns into a reported effect.
    """
    keys = standardise._estimand_keys()
    assert "beta" in keys
    assert "odds_ratio" not in keys
    assert not [k for k in keys if "odds" in k]


# --- 20.9  the audit entries and their step names -------------------------------------------------

def test_the_ten_entries_appear_in_19s_DRIVER_ORDER():
    """§20.9, §15, §19. The test runs §19's driver rather than an order of its own.

    Entry 7 lands between 6 and 8 because `support` precedes the restricted `all_centre` — `over=`
    is its output — and that order is a consequence of the call graph rather than a declaration.
    """
    original = config.N_BOOT
    config.N_BOOT = RUN_BOOT
    try:
        audit = audit_of()
        driver(four_centre_frame(), audit, with_inference=True)
    finally:
        config.N_BOOT = original

    steps = [(e.kind, e.step) for e in audit.entries]
    assert steps == [
        ("cohort", "standardisation_population"),
        ("missingness", "absence_by_standardisation_column"),
        ("model", "standardisation_fit"),
        ("model", "standardised_distributions"),
        ("model", "treated_support"),
        ("model", "support_baseline"),
        ("model", "standardisation_support_restricted"),
        ("model", "hierarchical_fit"),
        ("model", "standardised_hierarchical"),
        ("model", "standardisation_replicates"),
    ]
    assert len(steps) == 10


def test_a_driver_that_omits_the_restricted_all_centre_records_NINE_and_not_ten():
    """§20.9's companion. The failure the count catches, made explicit.

    §19 spells sensitivity 1's point-estimate call out precisely so that it cannot be forgotten, and
    this is what forgetting it looks like.
    """
    audit = audit_of()
    pop = standardise.population(four_centre_frame(), audit)
    standardise.all_centre(pop, audit)
    standardise.support(pop, audit)
    standardise.hierarchical(pop, audit)                # the restricted call SKIPPED
    assert len(audit.entries) == 8
    assert "standardisation_support_restricted" not in {e.step for e in audit.entries}


def test_no_step_name_collides_with_cohort_builds_and_the_ASSERTION_IS_A_SET_INTERSECTION():
    """§20.9, §4.2. `Audit.entry` is first-match, and both stages remove rows under kind='cohort'.

    In a Stage 14 driver, `cohort.build` and `standardise.population` both remove rows from the same
    classified frame under the same kind, and both are row-removal ledgers with case identifiers. A
    step name reused between them would make `audit.entry("cohort", ...)` return whichever ran first
    — and the two populations differ by 11 records and one whole centre.

    A set intersection and not a list, so a name added to either side is caught.
    """
    audit = audit_of()
    pop = standardise.population(four_centre_frame(), audit)
    standardise.all_centre(pop, audit)
    standardise.support(pop, audit)
    standardise.hierarchical(pop, audit)
    mine = {e.step for e in audit.entries}
    assert not (mine & COHORT_STEPS)

    # And the collision is real if the names are shared, which is what makes the assertion worth
    # making: `Audit.entry` returns the FIRST match, silently.
    both = data.Audit(hand_source(Path("."), "collision"))
    both.record("cohort", "shared", 1, "first")
    both.record("cohort", "shared", 99, "second")
    assert both.entry("cohort", "shared").n == 1


def test_entry_1_names_EXACTLY_as_many_case_identifiers_as_its_n_claims():
    """§20.9, §15. `cohort._record_removal`'s discipline, applied to a different rule.

    Stronger than `data.py`'s kind-keyed rule, which asks only that SOME case be named: a population
    that reported 22 removals and named 3 would satisfy `_MUST_NAME_CASES` and would be a population
    nobody could reconstruct from the log.
    """
    audit = audit_of()
    frame = four_centre_frame()
    pop = standardise.population(frame, audit)
    entry = audit.entry("cohort", "standardisation_population")
    assert entry is not None
    assert entry.n == len(entry.case_ids)
    assert entry.n == len(frame) - len(pop)


def test_record_removal_RAISES_when_it_cannot_name_what_it_removed():
    """§15's naming discipline is an assertion and not a comment."""
    frame = four_centre_frame().head(3).copy()
    frame.loc[frame.index[0], "case_id"] = pd.NA
    with pytest.raises(config.SchemaError) as excinfo:
        standardise._record_removal(audit_of(), "probe", frame, "detail", (("a",), ("b",)))
    assert "names" in message_of(excinfo)


@DATA_GATED
def test_the_Stage_2_to_4_path_records_TWELVE_entries_before_this_stage():
    """§15's ledger, measured. Stage 12 adds ten on top. [data-gated]"""
    df, audit = classified()
    assert len(audit.entries) == 12
    driver(df, audit)                                   # no inference: nine of the ten
    assert len(audit.entries) == 12 + 9


# --- 20.10  the support check ---------------------------------------------------------------------

@DATA_GATED
def test_the_TWO_READINGS_of_continuous_agree_on_this_workbook():
    """§20.10, §10.1. So "inert" is a CHECK and not a claim. [data-gated]

    [§14a] says *"each continuous covariate"*. `SUPPORT_COVARIATES` resolves it to "not a declared
    factor", which is the WIDER of the two readings; the narrow one boxes only the four genuinely
    continuous covariates. The wider reading is chosen because it cannot be wrong on a workbook where
    the narrow one is right — and this is what establishes that it is right here.
    """
    df, audit = classified()
    pop = standardise.population(df, audit)
    sup = standardise.support(pop, audit)

    narrow = ("age", "nihss_baseline", "core_ml", "tmax6_ml")
    assert set(narrow) <= set(sup.box)
    narrow_inside = standardise._inside(
        pop, {k: v for k, v in sup.box.items() if k in narrow})
    assert (narrow_inside == sup.inside).all()

    # The three that exclude nobody are the binaries and the linear integer, named so the reason is
    # readable: the treated arm covers both levels of each binary and reaches prestroke_mrs = 3.
    for covariate in ("sex", "atrial_fib", "prestroke_mrs"):
        low, high = sup.box[covariate]
        assert ((pop[covariate] >= low) & (pop[covariate] <= high)).all()


@DATA_GATED
def test_the_support_counts_and_that_EVERYONE_OUTSIDE_IS_A_COMPARATOR_PATIENT():
    """§20.10, §10.2. 93 inside, 11 outside, all 11 comparator; 12 never-IVT of whom 3 outside.

    The arm split is arithmetic rather than a finding — the box is the treated arm's own range, so no
    treated patient can be outside it — and it is asserted because a reader meeting "11 outside" will
    otherwise wonder. [data-gated]
    """
    df, audit = classified()
    pop = standardise.population(df, audit)
    sup = standardise.support(pop, audit)

    assert int(sup.inside.sum()) == 93
    assert int((~sup.inside).sum()) == 11
    outside = pop[~sup.inside]
    assert (outside[config.TREATMENT] == min(config.TREATMENT_LABELS)).all()
    assert sup.never_ivt == config.EXPECTED_NEVER_IVT
    assert sup.n_never_ivt == 12
    assert sup.n_never_ivt_outside == 3
    assert sum(sup.outside_by_centre.values()) == 11
    assert set(sup.outside_by_centre) == set(config.CENTER_ORDER)


def test_T9_fires_when_the_never_IVT_set_is_not_EXPECTED_NEVER_IVT():
    """§20.10, §20.14, T9. A `SchemaError` and not a finding.

    A workbook in which the never-IVT centre acquires a treated patient, or in which a fifth centre
    appears, changes what [§14a]'s support check is ABOUT — the never-IVT set is the thing the check
    is defined against — so the pipeline must stop rather than report a diagnostic whose subject moved.
    """
    frame = four_centre_frame()
    never = config.EXPECTED_NEVER_IVT[0]
    at_never = frame["center"] == never
    frame.loc[frame.index[at_never][0], config.TREATMENT] = float(max(config.TREATMENT_LABELS))
    pop = standardise.population(frame, audit_of())
    with pytest.raises(config.SchemaError) as excinfo:
        standardise.support(pop, audit_of())
    assert message_of(excinfo).startswith("T9")
    assert "EXPECTED_NEVER_IVT" in message_of(excinfo)


def test_the_never_IVT_set_is_COMPUTED_from_treating_centres_and_is_never_a_literal():
    """§10.2. `cohort.py` wrote this expression into its own docstring as the thing to use here.

    A centre's treatment availability is a property of the DATA: a workbook in which the never-IVT
    centre starts administering IVT must change this diagnostic rather than require someone to
    remember a list. `EXPECTED_NEVER_IVT` is the assertion target, never the operative rule.
    """
    tree = ast.parse(SOURCE)
    node = next(n for n in ast.walk(tree)
                if isinstance(n, ast.FunctionDef) and n.name == "_never_ivt")
    body = ast.unparse(node)
    assert "cohort.treating_centres" in body
    assert "C.CENTER_ORDER" in body
    # The expected set appears ONLY in the comparison, never as the thing being returned.
    assert body.count("C.EXPECTED_NEVER_IVT") == 2       # the `if` and the message
    assert "return C.EXPECTED_NEVER_IVT" not in body


@DATA_GATED
def test_the_baseline_table_is_NINETEEN_rows_with_FOUR_of_them_the_grouping_variable():
    """§20.10, §10.3. And `center = USZ` is `nan` for Stage 7 §5.3's branch 4. [data-gated]"""
    df, audit = classified()
    pop = standardise.population(df, audit)
    sup = standardise.support(pop, audit)

    expected_rows = (len([c for c in config.BALANCE_SET if c not in config.CATEGORICAL])
                     + sum(len(config.FACTOR_LEVELS[c]) for c in config.BALANCE_SET
                           if c in config.CATEGORICAL))
    assert len(sup.baseline) == expected_rows == 19

    grouping = [r for r in sup.baseline if r.role == standardise._GROUPING]
    assert len(grouping) == len(config.FACTOR_LEVELS[config.BOOT_STRATUM]) == 4
    assert len(sup.baseline) - len(grouping) == 15
    assert all(r.covariate.startswith(f"{config.BOOT_STRATUM} = ") for r in grouping)

    # Stage 7 §5.3's branch 4: `center = USZ` is constant in each group with the groups differing, so
    # the pooled SD is exactly zero and the SMD is `nan`. A future change to `smd`'s nan routes
    # surfaces HERE.
    usz = next(r for r in sup.baseline
               if r.covariate == f"{config.BOOT_STRATUM} = {config.EXPECTED_NEVER_IVT[0]}")
    assert np.isnan(usz.unweighted)
    assert usz.sd == 0.0 or np.isnan(usz.sd)


@DATA_GATED
def test_the_baseline_table_compares_TWO_DISJOINT_GROUPS_and_the_rest_are_in_NEITHER():
    """§10.3. 39 treated anywhere, 12 at the never-IVT centre, and 53 in neither. [data-gated]"""
    df, audit = classified()
    pop = standardise.population(df, audit)
    sup = standardise.support(pop, audit)

    treated = int((pop[config.TREATMENT] == max(config.TREATMENT_LABELS)).sum())
    assert treated == 39
    assert sup.n_never_ivt == 12
    assert len(pop) - treated - sup.n_never_ivt == 53
    entry = audit.entry("model", "support_baseline")
    assert entry.n == treated + sup.n_never_ivt == 51
    assert "53" in entry.detail                       # the caption names the excluded group


def test_the_baseline_table_calls_balance_smd_and_balance_levels_and_IMPLEMENTS_NEITHER():
    """§10.3, §24. A second `smd` would be a second definition of the yardstick; a second `levels`
    would be a second definition of WHAT A ROW IS — and that failure would be invisible here, because
    all three `onset_type` levels are present in the treated arm and a copy ranging over OBSERVED
    levels would produce a byte-identical table on v7."""
    tree = ast.parse(SOURCE)
    node = next(n for n in ast.walk(tree)
                if isinstance(n, ast.FunctionDef) and n.name == "_baseline")
    calls = {ast.unparse(c.func) for c in ast.walk(node) if isinstance(c, ast.Call)}
    assert "balance.smd" in calls
    assert "balance.levels" in calls
    # And no arithmetic of its own: no mean, no std, no var anywhere in the row builder.
    body = ast.unparse(node)
    for forbidden in (".mean(", ".std(", ".var(", "np.average"):
        assert forbidden not in body


def test_the_role_rule_is_STAGE_12s_and_names_no_propensity_model():
    """§10.3. `balance._role` is not reused, and the reason is CORRECTNESS rather than access.

    It returns `"propensity model"` for every `PS_COVARIATES` name, which is right for the
    specification it diagnoses and wrong here: **this stage fits no propensity model anywhere**
    ([§14]'s own first sentence), so a row labelled "propensity model" in a [§14a] table would name a
    model that does not exist in this analysis. Its docstring says a role is *"what this covariate is
    TO the specification being diagnosed"*, and the specification here is a different one. So this is
    a different rule, not a second copy of the same rule.
    """
    assert standardise._role(config.BOOT_STRATUM) == standardise._GROUPING
    assert "propensity" not in standardise._role("age").lower()
    assert "14a" in standardise._role("age")
    assert standardise._role("hypertension") == "negative control"
    assert standardise._role("penumbra_ml") == "excluded [§6]"
    # `center` is in BALANCE_SET and in PS_COVARIATES but NOT in STANDARDISATION_COVARIATES.
    assert config.BOOT_STRATUM in config.BALANCE_SET
    assert config.BOOT_STRATUM not in config.STANDARDISATION_COVARIATES
    for name in config.BALANCE_SET:
        assert "propensity" not in standardise._role(name).lower()


# --- 20.11  Stage 12 changes nothing upstream -----------------------------------------------------

def test_run_no_longer_owns_an_interval_loop_and_intervals_is_what_it_CALLS():
    """§20.11, §13.3. The extraction, asserted structurally so it cannot be undone by accident."""
    source = Path(bootstrap.__file__).resolve().read_text(encoding="utf-8")
    tree = ast.parse(source)
    node = next(n for n in tree.body
                if isinstance(n, ast.FunctionDef) and n.name == "run")
    body = ast.unparse(node)
    assert "for key, d in draws.items()" not in body
    assert "intervals(draws, _tested)" in body
    functions = [n.name for n in tree.body if isinstance(n, ast.FunctionDef)]
    assert "intervals" in functions and not functions.count("intervals") > 1


def test_bootstrap_intervals_p_RULE_IS_THE_PARAMETER_and_the_three_rules_behave():
    """§13.3. The only thing the callers disagree about is which keys get a `p`."""
    generator = np.random.default_rng(config.SEED)
    draws = {
        "beta": bootstrap.Draws("beta", generator.normal(0.3, 0.2, 200), 200, {}),
        "rd_2": bootstrap.Draws("rd_2", generator.normal(0.05, 0.1, 200), 200, {}),
        "sich.rd": bootstrap.Draws("sich.rd", generator.normal(0.0, 0.1, 200), 200, {}),
    }
    none = bootstrap.intervals(draws, lambda key: False)
    assert all(i.p is None for i in none.values())

    some = bootstrap.intervals(draws, bootstrap._tested)
    assert some["beta"].p is not None
    assert some["sich.rd"].p is not None
    assert some["rd_2"].p is None

    # And every interval carries the pinned method and the surviving-draw count, in both rules.
    for built in (none, some):
        for key, interval in built.items():
            assert interval.method == config.PERCENTILE_METHOD
            assert interval.level == config.CI_LEVEL
            assert interval.n_draws == len(draws[key].draws)


def test_an_estimand_below_the_floor_KEEPS_its_Draws_and_gets_NO_Interval():
    """§13.3, §8.3. `draws` always holds every key and `intervals` may hold fewer."""
    below = config.ci_min_draws() - 1
    thin = {"rd_0": bootstrap.Draws("rd_0", np.zeros(below), below + 1,
                                    {"nonconvergence": 1})}
    assert bootstrap.intervals(thin, lambda key: False) == {}
    # And one draw more is enough for an interval, so the floor is the floor and not an off-by-one.
    at_floor = {"rd_0": bootstrap.Draws("rd_0", np.linspace(-1.0, 1.0, below + 1), below + 1, {})}
    assert set(bootstrap.intervals(at_floor, lambda key: False)) == {"rd_0"}


def test_bootstrap_and_balance_expose_the_names_Stage_12_needs_and_reimplements_NEITHER():
    """§13.2, §10.3, §18. `bucket`, `collect` and `levels` are public and are CALLED, not copied."""
    assert callable(bootstrap.bucket) and callable(bootstrap.collect)
    assert callable(bootstrap.intervals) and callable(balance.levels)
    tree = ast.parse(SOURCE)
    calls = {ast.unparse(c.func) for c in ast.walk(tree) if isinstance(c, ast.Call)}
    assert {"bootstrap.bucket", "bootstrap.collect", "bootstrap.intervals",
            "bootstrap.replicates", "balance.levels", "balance.smd"} <= calls
    # And nothing private of another module's is reached for, with ONE declared exception:
    # `balance._pooled_sd` supplies the `sd` COLUMN of §10.3's table, because computing it here would
    # be the second definition of the yardstick §10.3 forbids.
    private = {c for c in calls if "._" in c}
    assert private == {"balance._pooled_sd"}, private


def test_model_and_balance_keep_their_landed_public_surfaces_plus_exactly_the_declared_additions():
    """§20.11, §18. Asserted BY NAME, so an accidental extraction is caught.

    `model.py` gains two and `balance.py` gains one, a rename. `polr`, `firth`, `design`, `predict`
    and `complete_cases` are untouched.
    """
    model_tree = ast.parse(Path(model.__file__).resolve().read_text(encoding="utf-8"))
    public = [n.name for n in model_tree.body
              if isinstance(n, ast.FunctionDef) and not n.name.startswith("_")]
    assert public == ["design", "complete_cases", "firth", "predict", "polr",
                      "ordinal_probabilities", "polr_ri"]

    balance_tree = ast.parse(Path(balance.__file__).resolve().read_text(encoding="utf-8"))
    assert [n.name for n in balance_tree.body
            if isinstance(n, ast.FunctionDef) and not n.name.startswith("_")] == [
        "smd", "levels", "assess"]
    assert not hasattr(balance, "_levels")

    bootstrap_tree = ast.parse(Path(bootstrap.__file__).resolve().read_text(encoding="utf-8"))
    assert [n.name for n in bootstrap_tree.body
            if isinstance(n, ast.FunctionDef) and not n.name.startswith("_")] == [
        "resample", "replicates", "percentile_ci", "bootstrap_p", "intervals", "bucket",
        "collect", "diagnostics", "run"]


def test_the_public_surface_is_FIVE_NAMES_and_the_privates_are_the_declared_ONES():
    """§3.2's count, re-derived rather than repeated, on `test_bootstrap.py`'s pattern."""
    tree = ast.parse(SOURCE)
    functions = [n.name for n in tree.body if isinstance(n, ast.FunctionDef)]
    assert [n for n in functions if not n.startswith("_")] == [
        "population", "all_centre", "support", "hierarchical", "inference"]
    assert sorted(n for n in functions if n.startswith("_")) == sorted([
        "_record_removal", "_population_table", "_assert_exposure_survived",
        "_assert_beta_reportable", "_assert_distribution", "_assert_identities",
        "_regime_design", "_arm_probabilities", "_expanded", "_standardise", "_over_mask",
        "_fit_table", "_distribution_table", "_box", "_inside", "_never_ivt", "_role",
        "_baseline", "_box_table", "_baseline_table", "_intercept_table", "_estimand_keys",
        "_standardised_values", "_replicate", "_diagnostics", "_tally", "_replicates_table"])


# The SHA-256 of `bootstrap.run`'s twenty-six intervals, canonicalised, **captured on this branch
# BEFORE `bootstrap.intervals` was extracted** and pinned here so that the extraction's byte-identity
# is a permanent guard rather than a one-time observation.
#
# **A HASH AND NOT THE VALUES, AND §4.5 IS WHY.** Stage 10 §4.5 forbids an interval limit in a
# version-controlled file and this module's own banner repeats it — a committed baseline holding the
# workbook's twenty-six limits would be exactly the thing that rule exists to prevent, and it would be
# worse than a limit in a spec because nobody reads a fixture. A digest pins all twenty-six to the bit
# and quotes none of them.
#
# It is a `repr`-based canonicalisation, so it pins the float bits and not a formatted rendering: a
# limit that moved in its last mantissa bit changes this hash.
STAGE10_INTERVALS_SHA256: Final[str] = (
    "aa764488e668985c297b809a977f54cb6c7a393ede80e657e1ce80a1271093f8")


@SLOW
@DATA_GATED
def test_bootstrap_runs_TWENTY_SIX_intervals_are_BYTE_IDENTICAL_across_the_extraction():
    """§20.11. The specific guard Stage 11 §16 said the `intervals` refactor would need. [slow]

    **The digest was captured on this branch before the extraction landed**, by running the Stage 1–10
    pipeline at `C.N_BOOT = 2000` with `run`'s own interval loop still in place. So this is not a
    determinism check against the shipped code's own output — it is the cross-refactor comparison
    §20.11 asks for, and it stays one for every future change to `run` or to `intervals`.

    Gated because it is a full 2000-replicate pipeline run, about three minutes. `pytest -rs` reads as
    an instruction: the default state of a slow measurement is "not run", and a check matching nothing
    must not pass as green.
    """
    import hashlib
    import json

    import outcome
    import propensity

    df, audit = classified()
    coh = cohort.build(df, audit)
    ps = propensity.fit(coh, audit)
    est = outcome.primary(coh, ps, audit)
    sec = outcome.secondary(coh, ps, audit)
    boot = bootstrap.run(coh, ps, est, sec, audit)

    got = {k: [repr(i.lo), repr(i.hi), i.level, i.method, i.n_draws, repr(i.p)]
           for k, i in boot.intervals.items()}
    assert len(got) == 26
    canonical = json.dumps(got, sort_keys=True, separators=(",", ":"))
    digest = hashlib.sha256(canonical.encode("utf-8")).hexdigest()
    assert digest == STAGE10_INTERVALS_SHA256, (
        "bootstrap.run's twenty-six intervals changed. The digest was captured BEFORE "
        "`bootstrap.intervals` was extracted, so a mismatch means the interval loop no longer "
        "reproduces Stage 10's landed output — which is the one thing Stage 11 §16 said this "
        "refactor had to be able to prove [Stage 12 §20.11].")


# --- 20.13  the replicate loop --------------------------------------------------------------------

def test_all_four_strata_appear_in_every_drawn_frame_and_each_size_is_EXACT():
    """§20.13, §13.1. Including the never-IVT stratum, which is resampled like the other three."""
    pop = standardise.population(four_centre_frame(), audit_of())
    sizes = pop.groupby(config.BOOT_STRATUM, observed=True).size().to_dict()
    assert len(sizes) == len(config.CENTER_ORDER)
    generator = np.random.default_rng(config.SEED)
    for _ in range(40):
        drawn = bootstrap.resample(pop, generator, config.BOOT_STRATUM)
        assert drawn.groupby(config.BOOT_STRATUM, observed=True).size().to_dict() == sizes


def test_the_drawn_frame_sequence_is_a_function_of_the_SEED_ALONE():
    """§20.13. Stage 10 §15.2's property, now over a THREE-ARM body — the fourth body to assert it.

    `resample` is called OUTSIDE `body` and there is no `try` in `replicates`, so the stream is a
    function of the seed and not of whether a replicate succeeded. An edit that moved the draw inside
    a failure path would make the seed stop identifying the replicates, silently.
    """
    pop = standardise.population(four_centre_frame(), audit_of())

    def frames(body, expect_raise=False):
        seen = []

        def wrapped(draw):
            seen.append(draw)
            return body(draw)

        # `replicates` has NO `try` — whatever the body raises leaves it, which is what makes the
        # stream a function of the seed and not of whether a replicate succeeded. So the raising
        # body's exception is caught here, and the frames drawn BEFORE it are what is compared.
        if expect_raise:
            with pytest.raises(model.FitError):
                bootstrap.replicates(pop, wrapped, 8, config.SEED, config.BOOT_STRATUM)
        else:
            bootstrap.replicates(pop, wrapped, 8, config.SEED, config.BOOT_STRATUM)
        return seen

    quiet = frames(lambda draw: None)
    noisy = frames(lambda draw: (_ for _ in ()).throw(model.FitError("polr: raised")),
                   expect_raise=True)
    real = frames(lambda draw: standardise._replicate(draw, standardise._estimand_keys()))
    assert len(quiet) == len(real) == 8
    assert len(noisy) == 1                       # it raised on the first, having drawn it
    for a, b in zip(quiet, real):
        pd.testing.assert_frame_equal(a, b)
    pd.testing.assert_frame_equal(quiet[0], noisy[0])


def test_a_replicate_produces_every_key_and_the_partition_holds():
    """§20.13, §3.1. `values` and `failures` partition the attempted keys — never both, never neither."""
    pop = standardise.population(four_centre_frame(), audit_of())
    keys = standardise._estimand_keys()
    generator = np.random.default_rng(config.SEED)
    draw = bootstrap.resample(pop, generator, config.BOOT_STRATUM)
    replicate = standardise._replicate(draw, keys)

    assert set(replicate.values) | set(replicate.failures) == set(keys)
    assert not set(replicate.values) & set(replicate.failures)
    assert replicate.n_alpha == len(model.polr(
        *(lambda X, _: (X, draw[config.PRIMARY_OUTCOME].to_numpy(dtype=float)))(
            *model.design(draw, (config.TREATMENT,) + config.STANDARDISATION_COVARIATES))).alpha)
    assert replicate.sum_w == float(len(draw))          # Sw = n; [§14a] weights nobody


def test_T1_costs_EVERY_key_and_carries_the_degenerate_design_bucket():
    """§20.13, §13.2. The first of the four failure groups."""
    keys = standardise._estimand_keys()
    # NOT through `population`: T7 fires first on this frame, because every patient being treated
    # leaves the comparator arm empty. T1 is reachable only inside a replicate, where `design` runs on
    # a drawn frame and `population` does not run again — which is exactly the path §13.2 describes.
    with pytest.raises(config.SchemaError) as excinfo:
        standardise.population(no_exposure_frame(), audit_of())
    assert message_of(excinfo).startswith("T7")

    replicate = standardise._replicate(no_exposure_frame(), keys)
    assert replicate.values == {}
    assert set(replicate.failures) == set(keys)
    assert set(replicate.failures.values()) == {"degenerate_design"}
    assert replicate.n_alpha is None and replicate.polr_iterations is None


def test_a_POOLED_polr_failure_costs_EVERY_key_hier_INCLUDED(monkeypatch):
    """§20.13, §13.2, §12.8. The assertion a `hier.*` exemption would fail.

    `polr_ri` starts from the pooled fit's `alpha` and `beta` — a nested model's exact maximiser in
    every coordinate but one, which is why 25 iterations suffice — so **there is no hierarchical arm
    without a pooled fit**. An exemption would require a cold start whose iteration count,
    convergence route and boundary rate are all unmeasured, and a quarter of the draws would then sit
    in a different numerical regime from the rest.
    """
    keys = standardise._estimand_keys()
    pop = standardise.population(four_centre_frame(), audit_of())
    monkeypatch.setattr(
        model, "polr",
        lambda *a, **k: (_ for _ in ()).throw(model.FitError("polr: no convergence, forced")))
    replicate = standardise._replicate(pop, keys)
    assert set(replicate.failures) == set(keys)
    assert "hier.sigma" in replicate.failures
    assert set(replicate.failures.values()) == {"nonconvergence"}


def test_a_polr_ri_failure_costs_hier_STAR_ONLY(monkeypatch):
    """§20.13, §13.2. The third failure group: 24 keys — 22 standardised, `hier.beta`, `hier.sigma`."""
    keys = standardise._estimand_keys()
    pop = standardise.population(four_centre_frame(), audit_of())
    monkeypatch.setattr(
        model, "polr_ri",
        lambda *a, **k: (_ for _ in ()).throw(model.FitError("polr_ri: no convergence, forced")))
    replicate = standardise._replicate(pop, keys)

    hier = {k for k in keys if k.startswith("hier.")}
    assert set(replicate.failures) == hier
    assert len(hier) == len(config.MRS_THRESHOLDS) + 2 + 2 * len(config.MRS_LEVELS) + 2
    assert set(replicate.failures.values()) == {"nonconvergence"}
    assert "beta" in replicate.values
    assert {k for k in keys if k.startswith("support.")} <= set(replicate.values)


def test_T12_on_the_POOLED_fit_costs_beta_and_NOTHING_ELSE(monkeypatch):
    """§20.13, §9.1. The one place this stage departs from Stage 8, and therefore the one most
    needing a test.

    Stage 8 drops the WHOLE replicate on separation, because the whole replicate's reported content
    is `beta` and six `RD_k` derived from the weighted empirical distributions. Here the coupling
    would discard 22 valid keys to suppress one invalid one, which is a bias in the surviving-draw set
    rather than a safeguard: every standardised probability is an averaged `expit` and stays in [0, 1]
    no matter how degenerate the fit is.
    """
    keys = standardise._estimand_keys()
    pop = standardise.population(four_centre_frame(), audit_of())
    real = model.polr

    def separated(X, y, w=None):
        fit = real(X, y, w)
        beta = fit.beta.copy()
        beta[fit.columns.index(config.TREATMENT)] = config.POLR_MAX_ABS_BETA + 1.0
        return dataclasses.replace(fit, beta=beta)

    monkeypatch.setattr(model, "polr", separated)
    replicate = standardise._replicate(pop, keys)

    assert set(replicate.failures) == {"beta"}
    assert replicate.failures["beta"] == "separation"
    # The 22 standardised keys of the pooled arm, and the whole support arm, survive intact.
    for key in keys:
        if key != "beta" and not key.startswith("hier."):
            assert key in replicate.values
    assert all(0.0 <= replicate.values[f"dist1_{level}"] <= 1.0 for level in config.MRS_LEVELS)
    assert all(-1.0 <= replicate.values[f"rd_{k}"] <= 1.0 for k in config.MRS_THRESHOLDS)
    assert "support.beta" not in replicate.failures     # there is none to lose


def test_T12_at_the_polr_ri_call_site_costs_hier_beta_and_NOTHING_ELSE(monkeypatch):
    """§20.12, §9.1. The second call site, and the same trade.

    `hier.*`'s twenty-two standardised keys and `hier.sigma` all survive.
    """
    keys = standardise._estimand_keys()
    pop = standardise.population(four_centre_frame(), audit_of())
    real = model.polr_ri

    def separated(X, y, groups):
        fit = real(X, y, groups)
        beta = fit.beta.copy()
        beta[fit.columns.index(config.TREATMENT)] = -(config.POLR_MAX_ABS_BETA + 1.0)
        return dataclasses.replace(fit, beta=beta)

    monkeypatch.setattr(model, "polr_ri", separated)
    replicate = standardise._replicate(pop, keys)

    assert set(replicate.failures) == {"hier.beta"}
    assert replicate.failures["hier.beta"] == "separation"
    assert "hier.sigma" in replicate.values
    assert all(f"hier.rd_{k}" in replicate.values for k in config.MRS_THRESHOLDS)
    assert "beta" in replicate.values


def test_inference_reconciles_every_one_of_the_SIXTY_NINE_keys():
    """§20.13, §3.1. `len(draws) + sum(failures) == n_attempted`, for every key."""
    original = config.N_BOOT
    config.N_BOOT = RUN_BOOT
    try:
        audit = audit_of()
        pop = standardise.population(four_centre_frame(), audit)
        boot = standardise.inference(pop, audit)
    finally:
        config.N_BOOT = original

    assert len(boot.draws) == 69
    for key, drawn in boot.draws.items():
        assert len(drawn.draws) + sum(drawn.failures.values()) == drawn.n_attempted
        assert drawn.n_attempted == RUN_BOOT
        assert np.all(np.isfinite(drawn.draws))
    assert all(interval.p is None for interval in boot.intervals.values())
    assert boot.seed == config.SEED


def test_every_standardised_draw_is_IN_RANGE_by_construction():
    """§9.1. The measurement the key-granularity guard rests on, as a property.

    Every standardised probability is an average of `expit` differences and every risk difference is a
    difference of two such averages, so the bounds hold whatever `beta` does. That is why
    `POLR_MAX_ABS_BETA` guards `beta` alone here.
    """
    original = config.N_BOOT
    config.N_BOOT = RUN_BOOT
    try:
        audit = audit_of()
        pop = standardise.population(four_centre_frame(), audit)
        boot = standardise.inference(pop, audit)
    finally:
        config.N_BOOT = original

    for key, drawn in boot.draws.items():
        if key.split(".")[-1].startswith("dist"):
            assert drawn.draws.min() >= 0.0 and drawn.draws.max() <= 1.0, key
        elif key.split(".")[-1].startswith(("rd_", "mrs_0_2", "mortality")):
            assert drawn.draws.min() >= -1.0 and drawn.draws.max() <= 1.0, key
        elif key.endswith("sigma"):
            assert drawn.draws.min() >= config.POLR_RI_SIGMA_FLOOR, key


def test_the_support_box_is_RECOMPUTED_in_every_replicate_and_not_held_at_the_point_estimates():
    """§11. The support is a STATISTIC, so holding it fixed would understate the interval.

    Measured on the workbook across 2000 replicates the in-support averaging population ranges 69 to
    102 against 93 at the point estimate — a spread wide enough that the choice is not cosmetic.
    """
    tree = ast.parse(SOURCE)
    node = next(n for n in ast.walk(tree)
                if isinstance(n, ast.FunctionDef) and n.name == "_replicate")
    body = ast.unparse(node)
    assert "_box(draw)" in body
    assert "_inside(draw," in body
    # And `support` itself is NOT called inside a replicate: T9 would kill the whole bootstrap when a
    # resample of a low-treated centre misses its bridging patients (~13% for a 2-of-30 centre).
    assert "standardise.support(" not in body
    calls = {ast.unparse(c.func) for c in ast.walk(node) if isinstance(c, ast.Call)}
    assert "support" not in calls
    assert "_never_ivt" not in calls


def test_the_n_average_spread_across_replicates_is_WIDE_on_the_fixture():
    """§11's measurement, as a property: the restricted averaging population genuinely varies."""
    pop = standardise.population(four_centre_frame(), audit_of())
    generator = np.random.default_rng(config.SEED)
    sizes = []
    for _ in range(40):
        draw = bootstrap.resample(pop, generator, config.BOOT_STRATUM)
        sizes.append(int(standardise._inside(draw, standardise._box(draw)).sum()))
    assert min(sizes) < max(sizes)
    assert max(sizes) <= len(pop)


# --- 20.12  the hierarchical arm's record shape ---------------------------------------------------

def test_Hierarchical_fit_IS_standardisation_fit_by_IDENTITY():
    """§20.12, §3.1. `is`, not `==`.

    Without it the `RIFit` would be reachable by two paths with nothing requiring them to agree,
    which is the duplication this record's shape exists to avoid.
    """
    pop = standardise.population(four_centre_frame(), audit_of())
    hier = standardise.hierarchical(pop, audit_of())
    assert hier.fit is hier.standardisation.fit
    assert isinstance(hier.fit, model.RIFit)
    assert hier.nodes == hier.fit.nodes == config.POLR_RI_NODES
    assert hier.sigma == hier.fit.sigma
    assert hier.at_floor == hier.fit.at_floor
    assert set(hier.intercepts) == set(hier.posterior_sd) == set(hier.fit.groups)


def test_the_three_arms_carry_THREE_DISTINCT_population_labels():
    """§20.12. So a formatter can tell the arms apart from the record alone."""
    audit = audit_of()
    pop = standardise.population(four_centre_frame(), audit)
    std = standardise.all_centre(pop, audit)
    sup = standardise.support(pop, audit)
    ssup = standardise.all_centre(pop, audit, over=sup.inside)
    hier = standardise.hierarchical(pop, audit)
    labels = {std.population, ssup.population, hier.standardisation.population}
    assert len(labels) == 3
    assert std.population == standardise._ALL_ELIGIBLE
    assert ssup.population == standardise._TREATED_SUPPORT
    assert hier.standardisation.population == standardise._RANDOM_INTERCEPT


def test_the_hierarchical_standardisation_CONDITIONS_on_b_hat_and_is_not_the_pooled_one():
    """§12.6. [§14a] prescribes the conditional reading and the two are different numbers.

    A random intercept is a per-centre shift of the cutpoints, so this branch is
    `ordinal_probabilities` unchanged against a per-centre `alpha`. If the arm standardised at `b = 0`
    instead it would be neither of [§14a]'s two readings, and the distributions would coincide with
    the pooled arm's up to the difference in `beta` alone.
    """
    audit = audit_of()
    pop = standardise.population(four_centre_frame(), audit)
    std = standardise.all_centre(pop, audit)
    hier = standardise.hierarchical(pop, audit)

    assert not hier.fit.at_floor, "the fixture carries a non-zero between-centre effect by design"
    assert any(abs(value) > 1e-6 for value in hier.intercepts.values())
    pooled = std.distribution[max(config.TREATMENT_LABELS)]
    conditioned = hier.standardisation.distribution[max(config.TREATMENT_LABELS)]
    assert any(abs(pooled[level] - conditioned[level]) > 1e-9 for level in config.MRS_LEVELS)

    # And the MECHANISM directly: a random intercept is a per-centre shift of the cutpoints, so
    # zeroing every `b_hat` must reproduce the fit's own unshifted prediction exactly.
    X, _ = model.design(pop, (config.TREATMENT,) + config.STANDARDISATION_COVARIATES)
    groups = pop[config.BOOT_STRATUM]
    shifted = standardise._arm_probabilities(hier.fit, X, groups)
    flattened = dataclasses.replace(hier.fit, b=np.zeros_like(hier.fit.b))
    unshifted = standardise._arm_probabilities(flattened, X, groups)
    assert not np.allclose(shifted, unshifted)
    plain = model.PolrFit(
        beta=hier.fit.beta, alpha=hier.fit.alpha, categories=hier.fit.categories,
        columns=hier.fit.columns, iterations=1, converged_on="likelihood",
        first_step_norm=0.0, rescales=0, halvings=0)
    assert np.array_equal(unshifted, model.ordinal_probabilities(plain, X))


def test_the_hierarchical_arm_raises_rather_than_averaging_over_an_UNSEEN_centre():
    """§12.6's structural guard. A patient at a centre the fit never saw has no `b_hat`.

    Averaging over a `nan` would make every reported probability `nan` rather than raising.
    """
    pop = standardise.population(four_centre_frame(), audit_of())
    X, dropped = model.design(pop, (config.TREATMENT,) + config.STANDARDISATION_COVARIATES)
    y = pop[config.PRIMARY_OUTCOME].to_numpy(dtype=float)
    fit = model.polr_ri(X, y, pop[config.BOOT_STRATUM])
    stranger = pop[config.BOOT_STRATUM].astype("object").copy()
    stranger.iloc[0] = "A CENTRE THE FIT NEVER SAW"
    with pytest.raises(config.SchemaError) as excinfo:
        standardise._arm_probabilities(fit, X, pd.Series(stranger))
    assert "no intercept" in message_of(excinfo)


# --- 20.14  every T identifier fires, and the SET is the assertion --------------------------------

def test_T5_fires_on_reordered_and_on_renamed_columns():
    """§20.14, T5. The message names the offending column.

    The coefficient vector is positional, so a mismatched or reordered design returns A NUMBER
    instead of an error.
    """
    pop = standardise.population(four_centre_frame(), audit_of())
    X, _ = model.design(pop, (config.TREATMENT,) + config.STANDARDISATION_COVARIATES)
    fit = model.polr(X, pop[config.PRIMARY_OUTCOME].to_numpy(dtype=float))

    reordered = X[list(X.columns)[::-1]]
    with pytest.raises(config.SchemaError) as excinfo:
        model.ordinal_probabilities(fit, reordered)
    assert message_of(excinfo).startswith("T5")

    renamed = X.rename(columns={"age": "AGE"})
    with pytest.raises(config.SchemaError) as excinfo:
        model.ordinal_probabilities(fit, renamed)
    assert message_of(excinfo).startswith("T5")
    assert "AGE" in message_of(excinfo)


def test_T8_fires_on_a_float_mask_on_a_SHUFFLED_index_and_on_a_wrong_LENGTH():
    """§20.14, T8. Three cases, one point about alignment.

    `over=` is applied inside every replicate to a RESAMPLED frame, and that frame is safe to index by
    label only because `bootstrap.resample` resets the index before renaming. T8 is the guard for the
    case that reset does not cover.
    """
    audit = audit_of()
    pop = standardise.population(four_centre_frame(), audit)
    inside = standardise._inside(pop, standardise._box(pop))

    with pytest.raises(config.SchemaError) as excinfo:
        standardise.all_centre(pop, audit, over=inside.astype(float))
    assert message_of(excinfo).startswith("T8")

    shuffled = inside.sample(frac=1.0, random_state=0)
    assert set(shuffled.index) == set(pop.index)
    with pytest.raises(config.SchemaError) as excinfo:
        standardise.all_centre(pop, audit, over=shuffled)
    assert message_of(excinfo).startswith("T8")

    with pytest.raises(config.SchemaError) as excinfo:
        standardise.all_centre(pop, audit, over=inside.iloc[:-3])
    assert message_of(excinfo).startswith("T8")

    # And the fourth, which the dtype test alone would miss: a boolean mask carrying pd.NA.
    gappy = inside.astype("boolean")
    gappy.iloc[0] = pd.NA
    with pytest.raises(config.SchemaError) as excinfo:
        standardise.all_centre(pop, audit, over=gappy)
    assert message_of(excinfo).startswith("T8")
    assert "missing" in message_of(excinfo)


def test_over_as_a_valid_nullable_boolean_mask_is_ACCEPTED():
    """§20.14's companion. `_inside`'s own product is a pandas comparison result, so both dtypes work."""
    audit = audit_of()
    pop = standardise.population(four_centre_frame(), audit)
    inside = standardise._inside(pop, standardise._box(pop))
    plain = standardise.all_centre(pop, audit, over=inside)
    nullable = standardise.all_centre(pop, audit, over=inside.astype("boolean"))
    assert plain.n_average == nullable.n_average
    assert plain.rd == nullable.rd


def test_EVERY_T_IDENTIFIER_IN_THE_SPEC_IS_EXERCISED_and_the_set_is_the_assertion():
    """§20.14. Asserted as a SET, on §20.9's pattern.

    A thirteenth identifier added to §14 without a test fails here rather than passing unnoticed. The
    left side is scraped from the two modules that raise them; the right side is the set this module's
    tests name in their own docstrings and bodies.
    """
    # §14's twelve, as an identifier -> MESSAGE-TOKEN mapping. **T2 and T3 do not carry their own
    # identifiers and that is deliberate**: both lead with `polr_ri:`, which is `"polr:"`'s
    # arrangement and shares its `nonconvergence` bucket. T4 leads with `T4` and NOT with `polr_ri:`,
    # because it is a degenerate design rather than a convergence failure — and that is the one
    # bucket assertion the raise-site scan cannot make, since `polr_ri:` is also a valid key (§14).
    series = {
        "T1": "T1", "T2": "polr_ri:", "T3": "polr_ri:", "T4": "T4", "T5": "T5", "T6": "T6",
        "T7": "T7", "T8": "T8", "T9": "T9", "T10": "T10", "T11": "T11", "T12": "T12"}
    assert len(series) == 12

    sources = (SOURCE, Path(model.__file__).resolve().read_text(encoding="utf-8"))
    raised = set()
    for module_source in sources:
        for node in ast.walk(ast.parse(module_source)):
            if not isinstance(node, ast.Raise) or node.exc is None:
                continue
            text = ast.unparse(node)
            for identifier, token in series.items():
                if f"{token} " in text or f"{token} the" in text:
                    raised.add(identifier)
    assert raised == set(series), sorted(set(series) - raised)

    own = Path(__file__).resolve().read_text(encoding="utf-8")
    model_tests = (Path(__file__).resolve().parent / "test_model.py").read_text(encoding="utf-8")
    both = own + model_tests
    tested = {identifier for identifier in series
              if f"_{identifier}_" in both or f"{identifier} fires" in both
              or f'"{identifier}"' in both or f"startswith(\"{identifier}\")" in both}
    assert tested == set(series), {"untested": sorted(set(series) - tested)}


def test_every_raise_site_token_in_this_module_is_a_FAILURE_BUCKETS_key_where_it_must_be():
    """§14, §18. `bootstrap.bucket` classifies by the FIRST TOKEN and RAISES on one it does not know.

    Only `model.FitError` raises are classified — a `SchemaError` is a contract break and is never
    caught — so this asserts the `FitError` tokens and asserts that they land in the bucket §18 names.
    """
    for node in ast.walk(ast.parse(SOURCE)):
        if not isinstance(node, ast.Raise) or node.exc is None:
            continue
        if not ast.unparse(node.exc).startswith("model.FitError"):
            continue
        text = ast.unparse(node)
        token = next(t for t in ("T1", "T12") if f"{t} " in text)
        assert token in config.FAILURE_BUCKETS, token
    assert config.FAILURE_BUCKETS["T1"] == "degenerate_design"
    assert config.FAILURE_BUCKETS["T12"] == "separation"
    assert config.FAILURE_BUCKETS["T4"] == "degenerate_design"
    assert config.FAILURE_BUCKETS["polr_ri:"] == "nonconvergence"
    assert set(config.FAILURE_BUCKETS.values()) == {
        "separation", "constant_outcome", "degenerate_design", "nonconvergence"}


# --- 4.3  what may not be quoted ------------------------------------------------------------------

def test_this_module_quotes_no_standardised_probability_and_no_centre_intercept():
    """§4.3, tightened once more by this stage. The rule is checkable and this is the check.

    Stage 8 §4.3 forbade case identifiers in `specs/`; Stage 10 §4.5 added interval limits and
    p-values; Stage 11 §1 added E-values. Stage 12 adds standardised probabilities and centre random
    intercepts. A standardised probability reads as a directly clinical number and a reader who meets
    one in a version-controlled file will quote it.

    The check is structural: no test in this file asserts a float equality against a workbook-derived
    standardised quantity. What it CAN pin is a count, a set, an identity or a tolerance.
    """
    tree = ast.parse(Path(__file__).resolve().read_text(encoding="utf-8"))
    gated = {n.name for n in ast.walk(tree) if isinstance(n, ast.FunctionDef)
             and any("DATA_GATED" in ast.unparse(d) for d in n.decorator_list)}
    assert gated, "the data-gated set must not be empty or this check is vacuous"
    for name in gated:
        node = next(n for n in ast.walk(tree)
                    if isinstance(n, ast.FunctionDef) and n.name == name)
        for constant in (c for c in ast.walk(node) if isinstance(c, ast.Constant)):
            if isinstance(constant.value, float) and 0.0 < abs(constant.value) < 1.0:
                # The only floats a data-gated test may name are tolerances: strictly tiny.
                assert abs(constant.value) < 1e-9, (name, constant.value)
