"""Acceptance tests for Stage 1 — §9 of `specs/stage1_config_and_data_contract.md`.

Tests that need the private workbook are marked `skipif(not DATA_XLSX.exists())`; every other test
runs on any checkout, against `tests/fixture_schema.xlsx`.

This file is exempt from the §7 raw-name scan, because it has to construct malformed headers to test
the contract check. That exemption is why the fixture generator lives here too:

    uv run python test_config.py --write-fixture

The generator builds the header row from `COLUMN_CONTRACT` keys rather than from its own list, so
the fixture cannot drift away from the contract by construction.
"""
from __future__ import annotations

import ast
import dataclasses
import hashlib
import os
import sys
from datetime import date
from pathlib import Path

import pandas as pd
import pytest

import config

DATA_GATED = pytest.mark.skipif(
    not config.DATA_XLSX.exists(),
    reason="the private workbook is gitignored and absent from this checkout")

MODULE_DIR = Path(config.__file__).resolve().parent


# --- the fixture ------------------------------------------------------------------------------
#
# Two invented rows, no patient data. What makes the fixture worth committing is that it reproduces
# the three structural hazards of the real workbook: three headers carry trailing whitespace; the
# final column is headerless, so a header-only read returns 41 columns while a full read returns 42;
# and one cell holds the literal three-character string 'N/A', which the workbook uses as a
# missingness sentinel and NA_VALUES turns into <NA>.

_HEADERLESS_PREFIX = "Unnamed: "

_FIXTURE_ROWS: list[dict[str, object]] = [
    {
        "case_id": "FIX-001", "center": 1, "age": 71, "sex": 0,
        "hypertension": 1, "hyperlipidemia": 0, "diabetes": 0, "smoking": 1, "atrial_fib": 0,
        "prestroke_mrs": 0, "wake_up": 1, "unwitnessed": 0, "nihss_baseline": 14,
        "ivt": 1, "onset_to_ivt_min": 180, "onset_to_groin_min": 250,
        "tmax6_ml": 88.0, "core_ml": 0.0, "penumbra_ml": 88.0,
        # The literal 'N/A' sentinel, in the column the workbook carries three of them in. It goes
        # here and not in the mRS column because `test_the_fixture_reads_under_the_declared_dtypes`
        # pins mRS at exactly one missing value, contributed by the second row.
        "mrs_90d": 2, "sich": 0, "ph2": 0, "tici_2b_3": "N/A",
        "contraindication_reason": None, "ivt_contraindicated": 0,
    },
    {
        # case_id is a bare integer here and a prefixed string above: the formats really do differ
        # by centre, which is why the contract types the column as `string`.
        "case_id": 7001, "center": "Lausanne", "age": 64, "sex": 1,
        "hypertension": 0, "hyperlipidemia": 1, "diabetes": 1, "smoking": 0, "atrial_fib": 1,
        "prestroke_mrs": 1, "wake_up": 0, "unwitnessed": 0, "nihss_baseline": 9,
        "ivt": 0, "onset_to_ivt_min": None, "onset_to_groin_min": 410,
        "tmax6_ml": 52.0, "core_ml": 12.0, "penumbra_ml": 40.0,
        # left missing on purpose, so the Int64-with-<NA> path is exercised by the fixture
        "mrs_90d": None, "sich": 0, "ph2": 0, "tici_2b_3": 0,
        "contraindication_reason": "Clinician decision", "ivt_contraindicated": 0,
    },
]

# The headerless final column carries a contact date at HUG only.
_FIXTURE_LAST_COLUMN = [date(2026, 1, 15), None]


def _write_fixture(path: Path = config.FIXTURE_XLSX) -> Path:
    """Rebuild `tests/fixture_schema.xlsx` from COLUMN_CONTRACT. Run with --write-fixture."""
    from openpyxl import Workbook

    headers = list(config.COLUMN_CONTRACT)
    by_analysis_name = {c.name: raw for raw, c in config.COLUMN_CONTRACT.items()
                        if c.name is not None}

    wb = Workbook()
    ws = wb.active
    ws.title = config.SHEET

    # The header row stops at 'Status ': the final column has no header at all, and the two empty
    # columns before it have neither header nor data. That is what makes a header-only read return
    # 41 columns here, exactly as it does against the real workbook.
    for j, raw in enumerate(headers, start=1):
        if not raw.startswith(_HEADERLESS_PREFIX):
            ws.cell(row=1, column=j, value=raw)

    for i, row in enumerate(_FIXTURE_ROWS, start=2):
        for analysis_name, value in row.items():
            if value is None:
                continue
            ws.cell(row=i, column=headers.index(by_analysis_name[analysis_name]) + 1, value=value)
        last = _FIXTURE_LAST_COLUMN[i - 2]
        if last is not None:
            ws.cell(row=i, column=len(headers), value=last)

    path.parent.mkdir(parents=True, exist_ok=True)
    wb.save(path)
    return path


def _read_fixture(**kwargs) -> pd.DataFrame:
    return pd.read_excel(config.FIXTURE_XLSX, sheet_name=config.SHEET, **kwargs)


@pytest.fixture(scope="module")
def fixture_columns() -> list[str]:
    return list(_read_fixture().columns)


# --- 9.1  contract exhaustiveness ---------------------------------------------------------------

def test_contract_accepts_a_full_read_of_the_fixture(fixture_columns):
    assert len(fixture_columns) == config.N_COLUMNS_EXPECTED
    assert fixture_columns == list(config.COLUMN_CONTRACT)
    assert config.assert_column_contract(fixture_columns) is None


def test_contract_rejects_an_added_column(fixture_columns):
    with pytest.raises(config.SchemaError) as e:
        config.assert_column_contract([*fixture_columns, "NewColumn"])
    assert "'NewColumn'" in str(e.value)


def test_contract_rejects_a_removed_column(fixture_columns):
    with pytest.raises(config.SchemaError) as e:
        config.assert_column_contract([c for c in fixture_columns if c != "Age"])
    assert "'Age'" in str(e.value)


def test_contract_rejects_a_trailing_whitespace_variant_with_a_hint(fixture_columns):
    stripped = [c.strip() if c == "PrestrokemRS " else c for c in fixture_columns]
    with pytest.raises(config.SchemaError) as e:
        config.assert_column_contract(stripped)
    message = str(e.value)
    assert "did you mean 'PrestrokemRS '" in message
    assert "note trailing space" in message


def test_contract_names_the_header_only_read_as_the_cause():
    header_only = list(_read_fixture(nrows=0).columns)
    assert len(header_only) == config.N_COLUMNS_EXPECTED - 1
    with pytest.raises(config.SchemaError) as e:
        config.assert_column_contract(header_only)
    message = str(e.value)
    assert "header-only read" in message
    assert "nrows=0" in message
    assert "Unnamed" not in message      # the cause, not a generic missing-column report


def test_contract_rejects_an_empty_iterable():
    with pytest.raises(config.SchemaError):
        config.assert_column_contract([])


def test_contract_accepts_a_generator(fixture_columns):
    assert config.assert_column_contract(c for c in fixture_columns) is None


@DATA_GATED
def test_fixture_headers_equal_the_real_workbook(fixture_columns):
    real = pd.read_excel(config.DATA_XLSX, sheet_name=config.SHEET)
    assert list(real.columns) == fixture_columns


# --- 9.2  no silent passthrough -----------------------------------------------------------------

_EMPTY_REASONS = {"", "-", "n/a", "tbd", "todo", "see above", "dropped"}


@pytest.mark.parametrize("raw", list(config.COLUMN_CONTRACT))
def test_every_column_carries_a_substantial_reason(raw):
    reason = config.COLUMN_CONTRACT[raw].reason
    assert len(reason) >= 20, f"{raw!r} has a {len(reason)}-character reason"
    assert reason.strip().lower() not in _EMPTY_REASONS


def test_the_contract_covers_the_whole_sheet():
    assert len(config.COLUMN_CONTRACT) == config.N_COLUMNS_EXPECTED == 42


def test_the_mapped_dropped_split():
    assert len(config.RENAME) == 25
    assert len(config.DROPPED) == 17
    assert len(config.RENAME) + len(config.DROPPED) == len(config.COLUMN_CONTRACT)


def test_analysis_names_are_unique():
    # Two raw headers mapping to one analysis name would make df.rename silently collapse them,
    # and one column's data would vanish with no error anywhere.
    assert len(config.RENAME) == len(config.ANALYSIS_NAMES) == 25


def test_no_analysis_name_is_also_a_raw_header():
    assert config.ANALYSIS_NAMES & set(config.COLUMN_CONTRACT) == set()


# --- 9.3  registry integrity ---------------------------------------------------------------------

def test_outcomes_registry_has_the_eight_section_5_entries():
    assert set(config.OUTCOMES) == {
        "mrs_90d", "mrs_0_2_90d", "mrs_0_1_90d", "tici_2b_3",
        "sich", "ph2", "death_90d", "mrs_5_6_90d"}


def test_exactly_one_outcome_is_primary():
    # [§8] rests on there being one primary quantity with one test. The pilot had two.
    assert sum(o.family == "primary" for o in config.OUTCOMES.values()) == 1


def test_the_three_families_partition_the_registry():
    families = {o.family for o in config.OUTCOMES.values()}
    assert families == {"primary", "secondary", "safety"}


@pytest.mark.parametrize("key", list(config.OUTCOMES))
def test_outcome_source_is_none_or_a_mapped_analysis_name(key):
    source = config.OUTCOMES[key].source
    assert source is None or source in config.ANALYSIS_NAMES


@pytest.mark.parametrize("key", list(config.OUTCOMES))
def test_op_and_threshold_are_none_exactly_when_source_is(key):
    o = config.OUTCOMES[key]
    assert (o.op is None) == (o.threshold is None) == (o.source is None)


@pytest.mark.parametrize("key", list(config.OUTCOMES))
def test_every_op_has_an_ops_entry(key):
    o = config.OUTCOMES[key]
    assert o.op is None or o.op in config.OPS


def test_outcome_rule_round_trips():
    assert config.OUTCOMES["mrs_0_2_90d"].rule == "mrs_90d <= 2"
    assert config.OUTCOMES["death_90d"].rule == "mrs_90d == 6"
    assert config.OUTCOMES["mrs_5_6_90d"].rule == "mrs_90d >= 5"
    assert config.OUTCOMES["mrs_90d"].rule is None


@pytest.mark.parametrize("covariate", config.PS_COVARIATES)
def test_every_ps_covariate_resolves(covariate):
    # Weaker than it was, and deliberately not repaired here. DERIVED_NAMES used to be the single
    # name `onset_type`; Stage 3 widened it to carry the two [§13] subgroup keys, so a subgroup
    # wrongly added to PS_COVARIATES would now pass this test. What restores the guarantee is
    # `test_no_subgroup_is_a_covariate` below — not a second opinion on the same property.
    assert covariate in config.ANALYSIS_NAMES or covariate in config.DERIVED_NAMES


# --- 9.3b  the Stage 3 registries [Stage 3 §12.13] --------------------------------------------------

def test_the_onset_flags_and_the_reference_level_are_exactly_the_declared_levels():
    # The one assertion that stops a level being renamed in ONSET_TYPE_FROM_FLAG alone. Stage 3
    # writes no level name at all, so if these two declarations drift, the factor silently gains a
    # level that FACTOR_LEVELS does not know about and Stage 6's categorical drops to <NA>.
    from_flags = tuple(level for _, level in config.ONSET_TYPE_FROM_FLAG)
    assert set(from_flags) | {config.REFERENCE_LEVELS["onset_type"]} == set(
        config.FACTOR_LEVELS["onset_type"])
    assert len(from_flags) == len(set(from_flags)) == len(config.FACTOR_LEVELS["onset_type"]) - 1


@pytest.mark.parametrize("flag, _level", config.ONSET_TYPE_FROM_FLAG)
def test_every_onset_flag_is_a_mapped_analysis_name(flag, _level):
    assert flag in config.ANALYSIS_NAMES
    assert flag in config.BINARY_COLUMNS      # so Stage 2's A5 covers its domain


def test_the_reference_level_is_not_produced_by_a_flag():
    # The baseline is what a record with no positive flag gets. A flag producing it too would make
    # the loop's last entry decide, silently.
    assert config.REFERENCE_LEVELS["onset_type"] not in {
        level for _, level in config.ONSET_TYPE_FROM_FLAG}


def test_derived_names_cannot_shadow_an_analysis_column():
    assert set(config.DERIVED_NAMES) & config.ANALYSIS_NAMES == set()
    assert len(config.DERIVED_NAMES) == len(set(config.DERIVED_NAMES))


def test_the_derived_dichotomies_are_exactly_the_outcomes_with_a_source():
    assert config.DERIVED_DICHOTOMIES == tuple(
        key for key, o in config.OUTCOMES.items() if o.source is not None)
    assert set(config.DERIVED_DICHOTOMIES) == {
        "mrs_0_2_90d", "mrs_0_1_90d", "death_90d", "mrs_5_6_90d"}


def test_the_cohort_dependent_subgroups_are_declared_subgroups():
    assert set(config.COHORT_DEPENDENT_SUBGROUPS) <= set(config.SUBGROUPS)


def test_row_wise_derived_is_derived_names_minus_the_cohort_dependent_ones():
    assert set(config.ROW_WISE_DERIVED) == set(config.DERIVED_NAMES) - set(
        config.COHORT_DEPENDENT_SUBGROUPS)
    assert config.ROW_WISE_DERIVED == tuple(
        n for n in config.DERIVED_NAMES if n in config.ROW_WISE_DERIVED)      # order preserved


def test_every_subgroup_is_a_derived_name():
    assert set(config.SUBGROUPS) <= set(config.DERIVED_NAMES)


@pytest.mark.parametrize("subgroup", sorted(config.SUBGROUPS))
def test_no_subgroup_is_a_covariate(subgroup):
    # A subgroup is a variable the analysis is *split* on, never one it adjusts for. Adjusting for
    # it inside its own subgroup analysis conditions on the split. This is also what restores the
    # guarantee 9.3 gave before DERIVED_NAMES widened to carry these two keys.
    everywhere = (
        set(config.PS_COVARIATES) | set(config.PS_COVARIATES_FULL)
        | set(config.STANDARDISATION_COVARIATES) | set(config.BALANCE_ONLY)
        | {c for override in config.OUTCOME_MODEL_OVERRIDES.values() for c in override})
    assert subgroup not in everywhere


@pytest.mark.parametrize("subgroup", sorted(config.SUBGROUPS))
def test_every_subgroup_carries_a_substantial_description(subgroup):
    assert len(config.SUBGROUPS[subgroup]) >= 20


def test_no_subgroup_is_post_time_zero():
    # Both are baseline attributes. A post-time-zero subgroup would condition on the future.
    assert set(config.SUBGROUPS) & config.POST_TIME_ZERO == set()


def test_the_subgroup_registry_is_exactly_the_two_that_survived_the_amendment():
    # Asserted against the literal, in the way this file pins OUTCOME_MODEL_OVERRIDES and data.py
    # pins SOURCES: the [§13] amendment of 2026-08-10 withdrew the third subgroup, and a third key
    # reappearing here is a statistical decision that must fail a test rather than pass quietly.
    assert set(config.SUBGROUPS) == {"unknown_onset", "core_above_median"}


# --- 9.3c  the Stage 4 eligibility registry [Stage 4 §12.12] ----------------------------------------
#
# Nothing covered these constants before Stage 4: ELIGIBILITY_ORDER existed but was read only by 9.8's
# immutability list, which asserts its type and not its contents.

def test_the_order_is_exactly_the_two_declared_classes():
    # Computed from the two, so this is the assertion that the computation was not rewritten as a
    # literal — and that a third class cannot appear in the display order alone. DECISION 1a
    # (2026-08-13) withdrew `indeterminate`; [§3]'s amendment of the same date is the protocol change.
    assert config.ELIGIBILITY_ORDER == (config.ELIGIBLE, config.INELIGIBLE)
    assert len(set(config.ELIGIBILITY_ORDER)) == len(config.ELIGIBILITY_ORDER) == 2


def test_the_retained_set_partitions_the_order_with_ineligible_alone_outside_it():
    # [§3]'s partition. Under DECISION 1a nothing sits in the difference between `!= ineligible` and
    # `== eligible` — the two coincide — so this now asserts the partition itself: a third class
    # added to either tuple without a decision behind it must fail here rather than move a cohort
    # quietly. It is the seam a [§13] arm over the undocumented controls would come back through.
    assert set(config.ELIGIBILITY_RETAINED) < set(config.ELIGIBILITY_ORDER)
    assert set(config.ELIGIBILITY_ORDER) - set(config.ELIGIBILITY_RETAINED) == {config.INELIGIBLE}


def test_the_eligibility_column_shadows_no_name_any_other_stage_produces():
    assert config.ELIGIBILITY not in config.ANALYSIS_NAMES
    assert config.ELIGIBILITY not in config.DERIVED_NAMES


def test_eligibility_is_not_a_covariate():
    # [§14a] states it plainly: "contraindication status is not a covariate, the cohort being already
    # restricted to eligible patients". A patient's eligibility decides who is in the population and
    # nothing else, so it must never adjust anything.
    #
    # This is also why the column stays out of DERIVED_NAMES [Stage 4 §8]: 9.3's resolution check
    # accepts any covariate that is in ANALYSIS_NAMES *or* DERIVED_NAMES, so widening the latter would
    # widen what may legally appear in a covariate list, and this assertion is what covers the gap.
    everywhere = (
        set(config.PS_COVARIATES) | set(config.PS_COVARIATES_FULL)
        | set(config.STANDARDISATION_COVARIATES) | set(config.BALANCE_ONLY)
        | {c for outcome in config.OUTCOMES for c in config.outcome_model_covariates(outcome)}
        | {c for override in config.OUTCOME_MODEL_OVERRIDES.values() for c in override})
    assert config.ELIGIBILITY not in everywhere


def test_eligibility_is_not_post_time_zero():
    # [§3] classifies "using only information available at time zero", so it is a baseline attribute.
    # Invariant 4's denylist is about variables that would be adjusted for; this one is barred from
    # adjustment by the assertion above instead, which names the property rather than the timing.
    assert config.ELIGIBILITY not in config.POST_TIME_ZERO


# --- 9.4  no raw names outside config.py ---------------------------------------------------------

def _raw_names_in_source(source: str) -> list[tuple[int, str]]:
    """Every string literal in `source` that is exactly a contract key, with its line number.

    A raw header can only enter code as a string literal, so that is the only thing worth checking.
    Comments and identifiers are not in the AST and are therefore invisible to this. f-strings are
    `ast.JoinedStr` with `ast.Constant` parts, so `f"{df['Age']}"` is still caught; a header
    assembled at runtime from concatenated fragments is not, and that residual hole is accepted.

    Exact and case-sensitive on purpose: `TICI_2b_3` -> `tici_2b_3` and `Center` -> `center` differ
    from their own analysis names only by case, so a case-insensitive scan would flag every
    legitimate downstream use of them and the repair would be to gut the scan.
    """
    contract = set(config.COLUMN_CONTRACT)
    return [(node.lineno, node.value)
            for node in ast.walk(ast.parse(source))
            if isinstance(node, ast.Constant) and isinstance(node.value, str)
            and node.value in contract]


def _python_files_under_scan() -> list[Path]:
    found: list[Path] = []
    for dirpath, dirnames, filenames in os.walk(MODULE_DIR):
        dirnames[:] = sorted(d for d in dirnames if d not in {".venv", "__pycache__"})
        for name in sorted(filenames):
            if name.endswith(".py") and name not in config.EXEMPT_FROM_RAW_NAME_SCAN:
                found.append(Path(dirpath) / name)
    return found


def test_every_test_module_lives_in_the_tests_directory():
    """The layout is a rule, not a convention, and this is what keeps it one.

    The shipped modules are flat at the project root — `import config as C` throughout — and every
    acceptance test lives in `tests/`, beside the fixtures it already shared a directory with
    (`fixture_schema.xlsx`, `reference/*.R`). A `test_*.py` reappearing at the root would still be
    collected, because `testpaths` is a default rather than a restriction, so it would work and the
    split would quietly stop being true.

    Two things this does NOT assert, deliberately. It does not forbid a nested directory under
    `tests/`, because a future stage may want one. And it does not require the shipped modules to be
    flat — that is Stage 1's decision and `pyproject.toml`'s `pythonpath` records it.
    """
    root = MODULE_DIR
    misplaced = sorted(p.name for p in root.glob("test_*.py"))
    assert misplaced == [], (
        "acceptance tests belong in tests/, not at the project root: " + ", ".join(misplaced))
    assert sorted(p.name for p in (root / "tests").rglob("test_*.py")), (
        "tests/ contains no test module — the collection root moved without the tests")


def test_no_raw_column_names_outside_config():
    offenders = [
        f"{path.relative_to(MODULE_DIR)}:{lineno}: {header!r}"
        for path in _python_files_under_scan()
        for lineno, header in _raw_names_in_source(path.read_text(encoding="utf-8"))]
    assert not offenders, (
        "raw workbook headers referenced outside config.py:\n  " + "\n  ".join(offenders))


def test_the_raw_name_scan_actually_fires():
    # A scan that silently matches nothing would otherwise pass as a green test.
    snippet = "prestroke = frame['PrestrokemRS ']\n"
    assert _raw_names_in_source(snippet) == [(1, "PrestrokemRS ")]


def test_the_raw_name_scan_is_case_sensitive():
    assert _raw_names_in_source("y = frame['tici_2b_3']\nz = frame['center']\n") == []


def test_the_exemption_list_is_exactly_the_three_declared_files():
    # The keys are BASENAMES, and the scan matches on basename too, so an exemption survives a file
    # moving between directories — which is what happened when the acceptance tests moved into
    # `tests/`. The existence check therefore searches the tree rather than joining onto MODULE_DIR:
    # the earlier `(MODULE_DIR / path).exists()` form asserted a location the exemption never claimed,
    # and it began failing for `test_config.py` the moment this file moved one level down.
    assert set(config.EXEMPT_FROM_RAW_NAME_SCAN) == {
        "config.py", "stage0_data_inventory.py", "test_config.py"}
    present = {path.name for path in MODULE_DIR.rglob("*.py")
               if ".venv" not in path.parts and "__pycache__" not in path.parts}
    for path, reason in config.EXEMPT_FROM_RAW_NAME_SCAN.items():
        assert path in present, f"{path} is exempted but does not exist"
        assert reason.strip(), f"{path} is exempted without a stated reason"


# --- 9.5  invariant 3 -----------------------------------------------------------------------------

@pytest.mark.parametrize("key", [k for k in config.OUTCOMES if k != "tici_2b_3"])
def test_outcome_model_covariates_defaults_to_the_shared_set(key):
    assert config.outcome_model_covariates(key) is config.PS_COVARIATES


def test_tici_declares_the_only_override():
    assert config.outcome_model_covariates("tici_2b_3") == ("center", "atrial_fib")
    # asserted against the literal, so a second override fails here rather than passing quietly
    assert set(config.OUTCOME_MODEL_OVERRIDES) == {"tici_2b_3"}


def test_outcome_model_covariates_rejects_an_unregistered_name():
    with pytest.raises(KeyError) as e:
        config.outcome_model_covariates("mrs_0_3_90d")
    assert "mrs_0_2_90d" in str(e.value)      # the message lists the known outcomes


@pytest.mark.parametrize("key", list(config.OUTCOMES))
def test_treatment_is_never_an_outcome_model_covariate(key):
    # Every estimator adds the treatment main effect itself; TICI's reduced model is
    # intercept + treatment + center (2 df) + atrial_fib = 5 parameters against 6 non-events.
    assert config.TREATMENT not in config.outcome_model_covariates(key)


# --- 9.6  invariant 4 -----------------------------------------------------------------------------

@pytest.mark.parametrize("name", [
    "PS_COVARIATES", "PS_COVARIATES_FULL", "STANDARDISATION_COVARIATES"])
def test_no_covariate_set_contains_a_post_time_zero_variable(name):
    assert set(getattr(config, name)) & config.POST_TIME_ZERO == set()


def test_the_denylist_covers_every_outcome():
    # Built from the registry, so a new outcome is denylisted by existing rather than by being
    # remembered.
    assert set(config.OUTCOMES) <= config.POST_TIME_ZERO
    assert {"onset_to_ivt_min", "onset_to_groin_min"} <= config.POST_TIME_ZERO


# --- 9.7  derived, not duplicated -----------------------------------------------------------------

def test_outcome_covariates_is_the_same_object_as_ps_covariates():
    assert config.OUTCOME_COVARIATES is config.PS_COVARIATES


def test_standardisation_covariates_drop_centre_and_nothing_else():
    assert "center" not in config.STANDARDISATION_COVARIATES
    assert config.STANDARDISATION_COVARIATES == tuple(
        c for c in config.PS_COVARIATES if c != "center")
    assert len(config.STANDARDISATION_COVARIATES) == len(config.PS_COVARIATES) - 1


def test_the_full_covariate_set_is_a_prefix_plus_exactly_four():
    n = len(config.PS_COVARIATES)
    assert config.PS_COVARIATES_FULL[:n] == config.PS_COVARIATES
    assert config.PS_COVARIATES_FULL[n:] == (
        "hypertension", "hyperlipidemia", "diabetes", "smoking")


def test_the_negative_controls_are_the_full_set_minus_the_propensity_model():
    # [Stage 7 §4.3] The [§6] negative controls and the [§13] full-covariate additions are the same
    # four names, so the set is computed from that identity rather than declared beside it. Asserted
    # against the two declarations it is computed from, never against a literal of four: a fifth risk
    # factor added to PS_COVARIATES_FULL must become a negative control in the same edit.
    assert config.NEGATIVE_CONTROLS == tuple(
        c for c in config.PS_COVARIATES_FULL if c not in config.PS_COVARIATES)
    assert set(config.NEGATIVE_CONTROLS) & set(config.PS_COVARIATES) == set()


def test_the_balance_set_is_the_confounders_plus_the_balance_only_names():
    # [Stage 7 §4.1, §4.5] The second assertion is what makes Stage 7's B4 unreachable on the declared
    # set: balance is a property of confounders measured at or before time zero, and a post-exposure
    # variable judged against [§9]'s threshold reads as a confounder that weighting failed to fix.
    assert config.BALANCE_SET == config.PS_COVARIATES + config.BALANCE_ONLY
    assert len(set(config.BALANCE_SET)) == len(config.BALANCE_SET)
    assert set(config.BALANCE_SET) & config.POST_TIME_ZERO == set()


# --- 9.8  immutability ------------------------------------------------------------------------------

@pytest.mark.parametrize("name", [
    "PS_COVARIATES", "BALANCE_ONLY", "CENTER_ORDER", "ELIGIBILITY_ORDER", "MRS_THRESHOLDS",
    "PS_COVARIATES_FULL", "STANDARDISATION_COVARIATES", "BINARY_COLUMNS", "DERIVED_NAMES",
    "CATEGORICAL", "EXPECTED_NEVER_IVT", "NA_VALUES",
    # Stage 3's four. This list is explicit rather than discovered, so a new declared sequence is
    # outside it until it is added by hand — which is the point, but it means adding it is part of
    # declaring one.
    "ONSET_TYPE_FROM_FLAG", "COHORT_DEPENDENT_SUBGROUPS", "DERIVED_DICHOTOMIES",
    "ROW_WISE_DERIVED",
    # Stage 4's one. ELIGIBILITY_ORDER is already above, from Stage 1.
    "ELIGIBILITY_RETAINED"])
def test_declared_sequences_are_tuples(name):
    # A later module doing `covs = config.PS_COVARIATES; covs.append(...)` would rewrite the
    # specification for every importer thereafter, including all 2000 refits in Stage 10, and no
    # import-time test could catch it. Tuples make that an AttributeError instead.
    assert isinstance(getattr(config, name), tuple)


def test_factor_levels_are_tuples():
    assert all(isinstance(v, tuple) for v in config.FACTOR_LEVELS.values())


def test_the_denylist_is_a_frozenset():
    assert isinstance(config.POST_TIME_ZERO, frozenset)
    assert isinstance(config.DROPPED, frozenset)
    assert isinstance(config.ANALYSIS_NAMES, frozenset)


def test_column_rejects_attribute_assignment():
    with pytest.raises(dataclasses.FrozenInstanceError):
        config.COLUMN_CONTRACT["Age"].name = "years"


def test_outcome_rejects_attribute_assignment():
    with pytest.raises(dataclasses.FrozenInstanceError):
        config.OUTCOMES["mrs_90d"].threshold = 3


# --- 9.9  read contract ------------------------------------------------------------------------------

def test_read_dtypes_name_only_mapped_contract_columns():
    assert set(config.READ_DTYPES) <= set(config.COLUMN_CONTRACT)
    assert set(config.READ_DTYPES) & config.DROPPED == set()


@pytest.mark.parametrize("raw", list(config.READ_DTYPES))
def test_every_declared_dtype_is_supported(raw):
    assert config.READ_DTYPES[raw] in {"string", "Int64", "Float64"}


def test_the_fixture_reads_under_the_declared_dtypes():
    df = _read_fixture(dtype=config.READ_DTYPES, na_values=list(config.NA_VALUES))
    assert df.shape == (len(_FIXTURE_ROWS), config.N_COLUMNS_EXPECTED)
    assert df["CaseID"].dtype == "string"
    assert str(df["Age"].dtype) == "Int64"
    # Int64 keeps missingness as <NA> rather than forcing the column to float, which is what makes
    # roadmap Stage 3's "every derived dichotomy has exactly the missingness of its source" clean.
    assert df["mRSscoreat90days"].isna().sum() == 1


def test_the_fixture_carries_a_literal_na_sentinel():
    # The workbook holds the three-character string 'N/A' in mRSscoreat90days, TICI_2b_3 and the six
    # 24-hour NIHSS columns. A default read cannot see them — 'N/A' is already in pandas' default NA
    # list — so the fixture has to carry one for any test of NA_VALUES to mean anything.
    raw = _read_fixture(dtype=object, keep_default_na=False)
    assert raw["TICI_2b_3"].iloc[0] == "N/A"


def test_the_literal_na_sentinel_reads_as_missing_under_the_contract():
    df = _read_fixture(dtype=config.READ_DTYPES, na_values=list(config.NA_VALUES))
    assert str(df["TICI_2b_3"].dtype) == "Int64"      # not coerced to object by the string
    assert df["TICI_2b_3"].isna().sum() == 1


# --- 9.10  factor declarations ------------------------------------------------------------------------

def test_the_factor_declarations_agree_with_each_other():
    assert set(config.FACTOR_LEVELS) == set(config.CATEGORICAL) == set(config.REFERENCE_LEVELS)


@pytest.mark.parametrize("factor", list(config.REFERENCE_LEVELS))
def test_every_reference_level_is_one_of_its_own_levels(factor):
    assert config.REFERENCE_LEVELS[factor] in config.FACTOR_LEVELS[factor]
    # And it is the LEADING one, which is the coupling Stage 6 §4.2's by-name reference drop rests on
    # and which nothing asserted until Stage 6 discovered it depended on it. On a declared Categorical,
    # `pd.get_dummies(..., drop_first=True)` drops the first *declared* level — so `design`'s by-name
    # drop and `drop_first` are the same model only while the reference IS that first level. Without
    # this line the equivalence is a coincidence that reordering CENTER_ORDER, or naming a non-leading
    # reference, breaks silently: the two forms would then fit different models under the same name.
    assert config.REFERENCE_LEVELS[factor] == config.FACTOR_LEVELS[factor][0]


def test_centre_levels_are_exactly_the_recode_targets():
    assert config.CENTER_ORDER == config.FACTOR_LEVELS["center"]
    assert set(config.CENTER_ORDER) == set(config.CENTER_RECODE.values())


def test_the_zero_bridging_centre_is_a_declared_centre():
    assert set(config.EXPECTED_NEVER_IVT) <= set(config.CENTER_ORDER)


# --- 9.11  constants sanity ---------------------------------------------------------------------------

def test_the_data_hash_is_well_formed():
    assert len(config.DATA_SHA256) == 64
    assert all(c in "0123456789abcdef" for c in config.DATA_SHA256)


def test_inference_constants_are_in_range():
    assert config.N_BOOT >= 1000
    assert 0 < config.SMD_THRESHOLD < 1
    assert config.RARE_MINORITY_THRESHOLD >= 1
    assert isinstance(config.SEED, int) and config.SEED > 0


def test_cumulative_thresholds_track_the_mrs_range():
    # So the [§8] cumulative thresholds cannot fall out of step with the mRS range.
    assert config.MRS_THRESHOLDS == tuple(range(config.PLAUSIBLE_RANGES["mrs_90d"][1]))


# --- Stage 8 §12  the three computed views the [§8] estimator reads ---------------------------------
#
# Six assertions, and each one pins what a Stage 8 constant is COMPUTED FROM rather than what it
# happens to equal. Written here rather than in `test_outcome.py` because the registry lives here.

def test_MRS_LEVELS_is_the_inclusive_plausible_range_for_the_primary_outcome():
    low, high = config.PLAUSIBLE_RANGES["mrs_90d"]
    assert config.MRS_LEVELS == tuple(range(low, high + 1))
    assert config.MRS_LEVELS == (0, 1, 2, 3, 4, 5, 6)


def test_the_mrs_plausible_range_has_an_upper_bound_at_all():
    """The declared type of PLAUSIBLE_RANGES is `tuple[int, int | None]` and `None + 1` is a
    TypeError, so §12's computed `MRS_LEVELS` rests on a fact the type does not guarantee.

    Three entries in that dict — core_ml, tmax6_ml, penumbra_ml — genuinely carry `None` above, so
    this is not hypothetical: it is the one column of the nine for which the computed form works, and
    nothing else says so.
    """
    assert config.PLAUSIBLE_RANGES["mrs_90d"][1] is not None


def test_MRS_LEVELS_and_MRS_THRESHOLDS_cannot_drift():
    """Stage 8 §8.3. RD_k stops at 5 because P(Y <= 6) is 1 in both arms by definition, so a seventh
    threshold would be a structural zero — and this identity is what stops the two constants being
    edited apart. MRS_THRESHOLDS is misfiled inside the [§13] subgroups block and is deliberately not
    moved; this assertion is what makes the misfiling harmless (Stage 8 §12, §17).
    """
    assert config.MRS_LEVELS[:-1] == config.MRS_THRESHOLDS


def test_PRIMARY_OUTCOME_is_the_unique_primary_family_key():
    # `next(...)` over a generator would take the FIRST of two silently, so the uniqueness half is
    # the assertion that notices. `test_exactly_one_outcome_is_primary` above pins the count; this
    # pins that PRIMARY_OUTCOME is that one, computed and never written as a literal.
    primary = [key for key, o in config.OUTCOMES.items() if o.family == "primary"]
    assert primary == [config.PRIMARY_OUTCOME]
    assert config.PRIMARY_OUTCOME == "mrs_90d"


def test_the_primary_outcome_is_ORDINAL_and_LOWER_IS_BETTER():
    """Stage 8 §7.1's precondition, asserted where the registry lives.

    Under `logit P(Y <= k) = alpha_k + beta*A` a positive beta raises P(Y <= k) in the treated arm at
    every k, so mass moves to LOW mRS. `exp(beta) > 1 favours bridging` therefore holds only for an
    ordinal outcome on which lower is better. For one on which higher is better the same coefficient
    favours the comparator, and the orientation statement [§8] requires would be exactly inverted
    with every number in the table still finite and still plausible. `outcome._orientation` raises on
    a registry that breaks this; this is the same fact asserted at the declaration.
    """
    outcome = config.OUTCOMES[config.PRIMARY_OUTCOME]
    assert outcome.kind == "ordinal"
    assert outcome.higher_is_better is False


def test_POLR_MAX_ABS_BETA_is_strictly_inside_the_measured_sparse_region():
    """Stage 8 §6.3, and after §14.7's split this is the ONLY assertion in the repository guarding
    that calibration.

    §14.7's band probe runs 240 fits and asserts the band *is* a band; the 4800-fit sweep that
    measured its endpoints, and the 1800-fits-per-cutpoint-count re-measurement that corrected them,
    are recorded in the spec's §20 and are deliberately not re-run on every commit. What replaces
    them is this comparison against both endpoints as literals.

    The endpoints are the CORRECTED ones — 11.037 and 18.98, measured across every cutpoint count a
    [§10] replicate can produce — and not the superseded (8.79, 18.81), which was measured at six
    cutpoints only. That matters in one direction: an edit back to the 10.0 an earlier draft used
    passes the superseded pair and FAILS this one, because legitimate fits reach 11.04 at odd cutpoint
    counts (§21.1 item 24).
    """
    assert 11.04 < config.POLR_MAX_ABS_BETA < 18.98
    assert config.POLR_MAX_ABS_BETA == 14.0


# --- Stage 9 §13  the two computed views the [§8] BINARY estimators read ----------------------------
#
# Five assertions, written here for the Stage 8 §12 block's reason: the registry lives here, and each
# one pins what a Stage 9 constant is COMPUTED FROM rather than what it happens to equal.

def test_BINARY_OUTCOMES_is_the_seven_binary_registry_keys_in_registry_order():
    binary = [key for key, o in config.OUTCOMES.items() if o.kind == "binary"]
    assert list(config.BINARY_OUTCOMES) == binary
    assert len(config.BINARY_OUTCOMES) == 7
    assert all(config.OUTCOMES[key].kind == "binary" for key in config.BINARY_OUTCOMES)


def test_BINARY_OUTCOMES_and_PRIMARY_OUTCOME_PARTITION_the_registry():
    """Stage 9 §4.1. The two together are every outcome, and they share none.

    This is the assertion that stops an outcome being estimated by NEITHER Stage 8 nor Stage 9, and
    the one that stops it being estimated by BOTH. A ninth `OUTCOMES` entry of some third kind would
    be silently unestimated without it — `secondary` ranges over BINARY_OUTCOMES and `primary` over
    PRIMARY_OUTCOME, so nothing else anywhere would notice.
    """
    assert config.PRIMARY_OUTCOME not in config.BINARY_OUTCOMES
    assert set(config.BINARY_OUTCOMES) | {config.PRIMARY_OUTCOME} == set(config.OUTCOMES)
    assert len(config.BINARY_OUTCOMES) + 1 == len(config.OUTCOMES)


def test_every_binary_outcome_declares_one_of_the_TWO_families():
    """Stage 9 §3.1 and [§13]. `Secondary.by_family` and Benjamini-Hochberg both assume two families.

    A third value in the registry would get its OWN correction group, silently: `by_family` builds its
    dict from whatever `family` says, so a typo or a new family partitions the seven into three groups
    and [§13]'s within-family correction is applied to a group of one. Nothing raises. Stage 9 §16
    item 7 records this as the only guard, which is why it is asserted rather than assumed.
    """
    assert {config.OUTCOMES[key].family for key in config.BINARY_OUTCOMES} == {"secondary", "safety"}


def test_OUTCOME_MODEL_OVERRIDES_can_only_name_a_BINARY_outcome():
    """Stage 9 §9.2. An override is a route INTO augmentation, so a key outside BINARY_OUTCOMES is an
    override nothing reads — for the ordinal outcome, which has no `m_a(X)`, or for a typo, which
    `outcome_model_covariates` would answer with the shared list and no complaint."""
    assert set(config.OUTCOME_MODEL_OVERRIDES) <= set(config.BINARY_OUTCOMES)


def test_OR_CONTINUITY_is_a_positive_pseudo_count_no_larger_than_the_conventional_half():
    """Stage 9 §6.3, §6.4. Asserted against BOTH endpoints, and neither is decoration.

    Zero or below is not a continuity correction — the odds ratio stays 0 or +inf, which is the thing
    the constant exists to prevent. Above 0.5 is past the conventional Haldane-Anscombe value, and
    §6.4 measured that 0.5 is ALREADY 2.9x more aggressive in the treated arm and 3.8x in the control
    arm than it is against row counts, because the four pseudo-counts here are WEIGHT-sums: Sw
    27.736623 over the ATO population against 92 rows. Raising it further widens an asymmetry that is
    measured to carry the odds ratio across the null in about 0.9% of empty-cell replicates.
    """
    assert 0.0 < config.OR_CONTINUITY <= 0.5
    assert config.OR_CONTINUITY == 0.5


def test_sex_labels_stay_none_until_the_query_is_answered():
    assert config.SEX_LABELS is None


@DATA_GATED
def test_the_workbook_hashes_to_the_declared_version():
    digest = hashlib.sha256(config.DATA_XLSX.read_bytes()).hexdigest()
    assert digest == config.DATA_SHA256, (
        f"expected {config.DATA_SHA256[:8]}…, actual {digest[:8]}…. The workbook changed. "
        "Re-run stage0_data_inventory.py, review ../out/stage0_data_inventory.md, and update "
        "DATA_SHA256 in the same commit.")


# --- 9.12  exclusions hold -------------------------------------------------------------------------

@pytest.mark.parametrize("excluded", ["penumbra_ml", "hir"])
def test_excluded_variables_appear_in_no_covariate_tuple(excluded):
    everywhere = (
        set(config.PS_COVARIATES) | set(config.PS_COVARIATES_FULL)
        | set(config.STANDARDISATION_COVARIATES)
        | {c for override in config.OUTCOME_MODEL_OVERRIDES.values() for c in override})
    assert excluded not in everywhere


def test_hir_has_no_analysis_name_at_all():
    assert "hir" not in config.ANALYSIS_NAMES


def test_balance_only_variables_are_not_adjusted_for():
    assert set(config.BALANCE_ONLY) & set(config.PS_COVARIATES) == set()


# --- 9.13  names resolve ------------------------------------------------------------------------------

# Stage 4's column is in the union although it is deliberately not a derived name [Stage 4 §8], so a
# future column-keyed constant naming it still resolves here rather than being read as a typo.
_RESOLVABLE = (config.ANALYSIS_NAMES | set(config.DERIVED_NAMES) | set(config.OUTCOMES)
               | {config.ELIGIBILITY})


@pytest.mark.parametrize("name", sorted(
    set(config.PS_COVARIATES_FULL) | set(config.BALANCE_ONLY) | set(config.BINARY_COLUMNS)
    | set(config.PLAUSIBLE_RANGES) | set(config.STRUCTURALLY_NON_APPLICABLE)
    | set(config.INFORMATIVE_ABSENCE) | set(config.POST_TIME_ZERO)
    # Stage 3's four. Like 9.8's list this union is explicit, so a typo in a new column-keyed
    # constant would otherwise surface as a KeyError inside derive() — or, worse, inside Stage 6.
    | set(config.SUBGROUPS) | set(config.COHORT_DEPENDENT_SUBGROUPS)
    | set(config.DERIVED_DICHOTOMIES) | set(config.ROW_WISE_DERIVED)))
def test_every_named_variable_can_be_produced(name):
    # A constant naming a variable no stage can produce is a typo that would otherwise surface as a
    # KeyError deep in Stage 6 or Stage 9, hours into a bootstrap.
    assert name in _RESOLVABLE


def test_the_two_kinds_of_absence_are_disjoint():
    # Stage 2's missingness table has one branch per constant. A column in both would be classified
    # by branch order rather than by declaration, and the one column at stake carries the
    # [DECISION 1] bit that separates eligible from indeterminate.
    assert set(config.STRUCTURALLY_NON_APPLICABLE) & set(config.INFORMATIVE_ABSENCE) == set()


def test_plausible_ranges_are_ordered():
    for name, (low, high) in config.PLAUSIBLE_RANGES.items():
        assert high is None or low < high, name


# --- 9.14  structural claims from Stage 0 -----------------------------------------------------------

@pytest.fixture(scope="module")
def workbook() -> pd.DataFrame:
    return pd.read_excel(config.DATA_XLSX, sheet_name=config.SHEET)


@DATA_GATED
def test_the_workbook_is_126_by_42(workbook):
    assert workbook.shape == (config.N_RECORDS_EXPECTED, config.N_COLUMNS_EXPECTED)
    assert config.assert_column_contract(workbook.columns) is None


@DATA_GATED
def test_a_header_only_read_of_the_workbook_loses_the_last_column():
    header_only = pd.read_excel(config.DATA_XLSX, sheet_name=config.SHEET, nrows=0)
    assert len(header_only.columns) == config.N_COLUMNS_EXPECTED - 1


@DATA_GATED
def test_case_identifiers_are_unique(workbook):
    assert not workbook["CaseID"].duplicated().any()


@DATA_GATED
def test_every_centre_code_is_recodable(workbook):
    assert set(workbook["Center"].astype(str)) <= set(config.CENTER_RECODE)


@DATA_GATED
def test_the_first_image_column_is_constant(workbook):
    assert workbook["FirstbrainimageMRI"].nunique(dropna=True) == 1


if __name__ == "__main__":
    if "--write-fixture" in sys.argv:
        print(f"wrote {_write_fixture()}")
    else:
        print(__doc__)
