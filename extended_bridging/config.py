"""Stage 1 — the configuration and data contract.

Every fact about the source workbook and every prespecified constant lives here, so that no later
module ever names a raw column, and so that a change to the workbook fails loudly here rather than
propagating silently into an estimate.

This stage declares; it does not compute. Two exceptions, both pure functions of their arguments and
neither touching the filesystem: :func:`assert_column_contract` and :func:`outcome_model_covariates`.
Both are called by later stages.

Section references in brackets are to ``statistical_analysis_plan.md``. DECISION *n* refers to the
decisions recorded in ``../out/stage0_data_inventory.md``. The specification for this module is
``specs/stage1_config_and_data_contract.md``; nothing here is invented outside it.

Standard library only, deliberately: importing the configuration is free, and the pure-logic
acceptance tests run without a scientific stack. Where a value exists to be handed to pandas
(``READ_DTYPES``, ``FACTOR_LEVELS``) it is a plain string or tuple that the consuming stage converts.
"""
from __future__ import annotations

import operator
from collections.abc import Iterable
from dataclasses import dataclass
from pathlib import Path
from typing import Final


# --- exceptions -------------------------------------------------------------------------------

class SchemaError(Exception):
    """The workbook's columns do not match COLUMN_CONTRACT."""


class DataVersionError(Exception):
    """The workbook's content hash does not match DATA_SHA256."""


# --- structured types -------------------------------------------------------------------------

@dataclass(frozen=True)
class Column:
    """One column of the source workbook, keyed in COLUMN_CONTRACT by its verbatim raw header."""

    name: str | None            # analysis name, or None if dropped
    reason: str                 # why it is dropped, or what it is. Never empty, ever.
    dtype: str | None = None    # pandas dtype passed to read_excel; None means infer


@dataclass(frozen=True)
class Outcome:
    """One entry of the [§5] outcome registry."""

    label: str                  # for tables and figures
    kind: str                   # "ordinal" | "binary"
    family: str                 # "primary" | "secondary" | "safety"
    source: str | None          # analysis name of the ordinal source, or None if read directly
    op: str | None              # "<=" | ">=" | "==" ; None iff source is None
    threshold: int | None       # None iff source is None
    higher_is_better: bool

    @property
    def rule(self) -> str | None:
        """Human-readable derivation, e.g. 'mrs_90d <= 2'. Derived, never stored.

        The threshold is stored once, as data. Stage 3 applies it through ``OPS`` and this prose is
        computed from the same fields, so the manuscript cannot print `mRS 0-1` while the analysis
        computes `mRS 0-2`.
        """
        if self.source is None:
            return None
        return f"{self.source} {self.op} {self.threshold}"


OPS = {"<=": operator.le, ">=": operator.ge, "==": operator.eq}


# --- the column contract  [Stage 1 acceptance criterion] --------------------------------------
#
# Keyed by the *verbatim* raw header, including its whitespace. Three headers carry trailing spaces
# ('PrestrokemRS ', 'mRS56at90days   ', 'Status ') and two carry interior spaces. The contract must
# not strip them and the loader must not strip headers before matching: stripping either side turns
# a schema change into a silent rename.
#
# The final column is headerless, so pandas invents 'Unnamed: 41' for it — and only once rows are
# read. See assert_column_contract for why that matters.
#
# Entries are in **workbook order**, not grouped by fate: the three invented 'Unnamed: N' keys
# encode their own position in the sheet, so any other ordering would be internally inconsistent.
# A dropped column is the ones written `Column(None, …)`; there are 17 of them and 25 mapped.

COLUMN_CONTRACT: Final[dict[str, Column]] = {
    "CaseID": Column(
        "case_id", "patient identifier; formats differ by centre (SSR-HUG-…, L-…, bare integers), "
        "so inference yields mixed int/str and it must be read as text", "string"),
    "Center": Column(
        "center", "recruiting centre; HUG is the integer 1 and the others are names, recoded via "
        "CENTER_RECODE [§6 covariate, and the [§3] restriction is applied on it]", "string"),
    "Age": Column(
        "age", "[§6] propensity and outcome model covariate; age in years at admission", "Int64"),
    "Sex": Column(
        "sex", "[§6] propensity and outcome model covariate; 0/1 coding is undocumented, see §10 "
        "and SEX_LABELS", "Int64"),
    "MedHistHypertension": Column(
        "hypertension", "[§6] vascular risk factor; deliberately excluded from the propensity "
        "model and reported in the balance table as a negative control", "Int64"),
    "MedHistHyperlipidemia": Column(
        "hyperlipidemia", "[§6] vascular risk factor; deliberately excluded from the propensity "
        "model and reported in the balance table as a negative control", "Int64"),
    "MedHistDiabetes": Column(
        "diabetes", "[§6] vascular risk factor; deliberately excluded from the propensity model "
        "and reported in the balance table as a negative control", "Int64"),
    "MedHistSmoking": Column(
        "smoking", "[§6] vascular risk factor; deliberately excluded from the propensity model "
        "and reported in the balance table as a negative control", "Int64"),
    "MedHistAtrialFibr": Column(
        "atrial_fib", "[§6] propensity and outcome model covariate; also the clot-composition "
        "proxy in the TICI outcome model [DECISION 3]", "Int64"),
    "PrestrokemRS ": Column(
        "prestroke_mrs", "[§6] propensity and outcome model covariate; pre-stroke functional "
        "baseline. Trailing space in the raw header is deliberate", "Int64"),
    "Wakeupstroke": Column(
        "wake_up", "source of the derived onset_type factor [§5]; never positive at the same time "
        "as unwitnessed, which Stage 3 asserts", "Int64"),
    "Unwitnessedstroke": Column(
        "unwitnessed", "source of the derived onset_type factor [§5]; never positive at the same "
        "time as wake_up, which Stage 3 asserts", "Int64"),
    "NIHSSonadmission": Column(
        "nihss_baseline", "[§6] propensity and outcome model covariate; stroke severity at "
        "admission, before time zero", "Int64"),
    "FirstbrainimageMRI": Column(
        None, "constant 0 across all 126 records; a centre-pathway attribute largely collinear "
        "with center, not a patient characteristic, and excluded by kind rather than by variance "
        "[§6]"),
    "IVTwithrtPA": Column(
        "ivt", "the treatment: 1 = IVT before EVT (bridging), 0 = EVT alone. Named by TREATMENT",
        "Int64"),
    "TimefromONSETtoIVTmin": Column(
        "onset_to_ivt_min", "descriptive only; structurally non-applicable in the control arm, "
        "which never received IVT [§11], and denylisted as post-time-zero"),
    "TimefromONSETtogroinmin": Column(
        "onset_to_groin_min", "post-exposure — bridging can itself delay groin puncture [§12] — so "
        "it is reported by arm and denylisted, never adjusted for"),
    "HypoperfusedtissuevolumeTmax6sml": Column(
        "tmax6_ml", "[§6] propensity and outcome model covariate; hypoperfused tissue volume, "
        "Tmax > 6 s, in millilitres"),
    "IschemiccorevolumeCBF30ml": Column(
        "core_ml", "[§6] propensity and outcome model covariate; ischaemic core volume, "
        "CBF < 30%, in millilitres"),
    "Penumbravolumeml": Column(
        "penumbra_ml", "excluded from every model [§6]; balance-table and [§13] subgroup use only, "
        "because it is core and Tmax>6 s restated"),
    "HIR": Column(
        None, "excluded a priori from every model and every sensitivity analysis [§6]; it has no "
        "analysis name at all so that no covariate list can name it"),
    "NIHSSscoreat24h": Column(
        None, "post-time-zero, and not in the [§5] outcome registry. Analysed in the pilot; "
        "restoring it requires amending [§5] first"),
    "NIHSSscore02at24 h": Column(
        None, "post-time-zero, and not in the [§5] outcome registry. Analysed in the pilot; "
        "restoring it requires amending [§5] first"),
    "8pointsNIHSSreductionat24 h": Column(
        None, "post-time-zero, and not in the [§5] outcome registry. Analysed in the pilot; "
        "restoring it requires amending [§5] first"),
    "NIHSS01at24h": Column(
        None, "post-time-zero, and not in the [§5] outcome registry. Analysed in the pilot; "
        "restoring it requires amending [§5] first"),
    "EarlyNeurologicalRecovery": Column(
        None, "post-time-zero, and not in the [§5] outcome registry. Analysed in the pilot; "
        "restoring it requires amending [§5] first"),
    "ChangeinNIHSSscorefrombaselineto24h": Column(
        None, "post-time-zero, and not in the [§5] outcome registry. Analysed in the pilot; "
        "restoring it requires amending [§5] first"),
    "mRSscoreat90days": Column(
        "mrs_90d", "the primary outcome, and the ordinal source of all four derived dichotomies. "
        "Ground truth for vital status [DECISION 2]", "Int64"),
    "mRSscore01at90days": Column(
        None, "shipped dichotomy; rebuilt from mrs_90d [§5]. Agrees with the ordinal source on "
        "every record, and is still not read"),
    "Deathat7days": Column(
        None, "not in the [§5] outcome registry. Analysed in the pilot; restoring it requires "
        "amending [§5] first"),
    "mRS02at90days": Column(
        None, "shipped dichotomy; rebuilt from mrs_90d [§5]. Agrees with the ordinal source on "
        "every record, and is still not read"),
    "Deathat90days": Column(
        None, "shipped dichotomy; rebuilt as mrs_90d == 6 [DECISION 2]. Disagrees with the ordinal "
        "source on 2 records, which is a standing data query, not something to resolve silently"),
    "mRS56at90days   ": Column(
        None, "shipped dichotomy; rebuilt as mrs_90d >= 5 [DECISION 2]. Disagrees with the ordinal "
        "source on 3 records. Trailing space in the raw header is deliberate"),
    "Symptomaticintracranialhaemorrhage": Column(
        "sich", "[§5] safety outcome, read directly; no ordinal source exists for it", "Int64"),
    "Parenchymalhaematomatype2": Column(
        "ph2", "[§5] safety outcome, read directly; no ordinal source exists for it", "Int64"),
    "TICI_2b_3": Column(
        "tici_2b_3", "[§5] secondary outcome, read directly; no ordinal source exists, and it "
        "carries the only outcome-model override [DECISION 3]", "Int64"),
    "Contraindications_to_IVT": Column(
        "contraindication_reason", "presence only — the free text is never read [DECISION 1]. Its "
        "one bit is whether a reason was recorded at all, which separates eligible from "
        "indeterminate [§3]", "string"),
    "IVT_contraindicated_binary": Column(
        "ivt_contraindicated", "the [§3] eligibility classifier [DECISION 1]; Stage 4 asserts it "
        "is 0/1, never missing, and never 1 for a treated patient", "Int64"),
    "Unnamed: 38": Column(
        None, "empty in every record; carried by the sheet's dimensions, not by any data"),
    "Unnamed: 39": Column(
        None, "empty in every record; carried by the sheet's dimensions, not by any data"),
    "Status ": Column(
        None, "free-text administrative field (DCD, CG, dated contact notes). Vital-status-adjacent "
        "and deliberately unread: the 90-day mRS is the sole ground truth for vital status "
        "[DECISION 2]. Trailing space in the raw header is deliberate"),
    "Unnamed: 41": Column(
        None, "headerless HUG-only contact date, 13 values; post-time-zero administrative data. "
        "pandas invents this name, and only on a full read"),
}


# --- derived views of the contract -------------------------------------------------------------
#
# Every consumer reads one of these, never COLUMN_CONTRACT directly. All four are computed from the
# contract; none may be written out as a second literal.

RENAME: Final[dict[str, str]] = {
    raw: c.name for raw, c in COLUMN_CONTRACT.items() if c.name is not None}
DROPPED: Final[frozenset[str]] = frozenset(
    raw for raw, c in COLUMN_CONTRACT.items() if c.name is None)
ANALYSIS_NAMES: Final[frozenset[str]] = frozenset(RENAME.values())
READ_DTYPES: Final[dict[str, str]] = {
    raw: c.dtype for raw, c in COLUMN_CONTRACT.items() if c.dtype is not None}


# --- paths and data identity --------------------------------------------------------------------
#
# Every path is anchored to this module file, never to the working directory: Stage 14 is a single
# entry point that a reader will plausibly invoke from the repository root.

ROOT: Final[Path] = Path(__file__).resolve().parents[1]          # the repository root
DATA_XLSX: Final[Path] = ROOT / "data" / (
    "Excel_bridging_EXTEND_paper_HUG_CHUV_LUGANO_USZ_def_v7_july26_"
    "with_abs_contra_indication.xlsx")
SHEET: Final[str] = "Feuil1"
OUT: Final[Path] = ROOT / "out"
TABLES, FIGURES, LOGS = OUT / "tables", OUT / "figures", OUT / "logs"

FIXTURE_XLSX: Final[Path] = Path(__file__).resolve().parent / "tests" / "fixture_schema.xlsx"

N_RECORDS_EXPECTED: Final[int] = 126
N_COLUMNS_EXPECTED: Final[int] = 42
# Load-bearing, not a guard. The workbook holds the literal three-character string 'N/A' in
# mRSscoreat90days (2 cells), TICI_2b_3 (3) and each of the six 24-hour NIHSS columns (1 each) — so
# it is doing real work on the primary outcome and on a secondary one. The Stage 0 note reports
# those cells as blanks because 'N/A' is already in pandas' default NA list, which is also why
# `keep_default_na` must stay True: setting it False would leave 'N/A' as a string in an Int64
# column and crash the read. The repair for that crash is to name the specific sentinel here, never
# to widen this tuple blindly.
NA_VALUES: Final[tuple[str, ...]] = ("N/A", "")

# The column contract catches schema changes only. It catches nothing if the data owner returns a
# corrected file with identical headers and two amended mRS values — and §10 records exactly two
# open queries, so that is expected. Stage 2 hashes DATA_XLSX and raises DataVersionError on
# mismatch. The point is not to forbid a new workbook: it is to make adopting one a deliberate
# commit that also re-runs Stage 0.
DATA_SHA256: Final[str] = "54934fbb2ae22647a9c0a2cbaff7ac425ed8948bcaeabe36d69aca00df657371"


# --- inference and thresholds ---------------------------------------------------------------

SEED: Final[int] = 20260807       # recorded in the run summary alongside DATA_SHA256
N_BOOT: Final[int] = 2000
SMD_THRESHOLD: Final[float] = 0.10     # [§9]

# Named for the minority cell, min(events, non-events), not for events [§8 amendment, DECISION 3].
# A constant named RARE_EVENT_THRESHOLD would invite the exact misreading the amendment corrects:
# TICI 2b-3 has 114 events and 6 non-events, so it passes an event-count rule while failing that
# rule's rationale.
RARE_MINORITY_THRESHOLD: Final[int] = 10


# --- treatment and centres -----------------------------------------------------------------

TREATMENT: Final[str] = "ivt"                       # 1 = IVT before EVT (bridging), 0 = EVT alone
TREATMENT_LABELS: Final[dict[int, str]] = {0: "EVT alone", 1: "bridging"}

CENTER_RECODE: Final[dict[str, str]] = {
    "1": "HUG", "Lausanne": "CHUV", "Lugano": "Lugano", "USZ": "USZ"}

# Display and factor order. Observed sizes are HUG 52, CHUV 21, Lugano 31, USZ 22, so this is *not*
# a size ordering: HUG leads because it is the largest and the declared reference level, and the
# remaining three keep the pilot's order so tables stay comparable.
CENTER_ORDER: Final[tuple[str, ...]] = ("HUG", "CHUV", "Lugano", "USZ")

# An assertion target, never an operative rule. Stage 5 derives the zero-bridging set from the data
# and asserts it equals this; hardcoding the exclusion instead would keep excluding USZ if a
# corrected workbook gave it treated patients.
EXPECTED_NEVER_IVT: Final[tuple[str, ...]] = ("USZ",)


# --- factor levels and reference categories ---------------------------------------------------
#
# Stage 6 builds each factor as a categorical with the declared level set, so a bootstrap replicate
# missing a level yields an all-zero column that the "constant columns dropped" rule then removes —
# deterministically, in every replicate. Without the declaration, pd.get_dummies(drop_first=True)
# would silently baseline on CHUV, and a replicate containing no Lugano patients would produce a
# design matrix of different width whose columns mean different things.

CATEGORICAL: Final[tuple[str, ...]] = ("onset_type", "center")

FACTOR_LEVELS: Final[dict[str, tuple[str, ...]]] = {
    "center": CENTER_ORDER,
    "onset_type": ("witnessed", "unwitnessed", "wake_up"),
}
REFERENCE_LEVELS: Final[dict[str, str]] = {
    "center": "HUG",             # largest centre, n = 52
    "onset_type": "witnessed",   # the clinical baseline
}

# Which onset flag produces which non-baseline level of `onset_type` [§5]. The baseline is
# REFERENCE_LEVELS["onset_type"], which is already declared and is Stage 6's reference category, so
# it is deliberately not repeated here. test_config.py asserts that this tuple's levels plus that
# baseline are exactly FACTOR_LEVELS["onset_type"], so a level cannot be added in one place only.
#
# Stage 3 loops over this tuple rather than writing the level names, and accumulates each flag's
# missingness over the same loop: a third onset flag then contributes both its level and its
# missingness in one edit. Its order is not a precedence rule — Stage 3 §4.1 asserts no record is
# positive on both flags, so no record can match twice.
ONSET_TYPE_FROM_FLAG: Final[tuple[tuple[str, str], ...]] = (
    ("unwitnessed", "unwitnessed"),
    ("wake_up", "wake_up"),
)


# --- covariates [§6] --------------------------------------------------------------------------

PS_COVARIATES: Final[tuple[str, ...]] = (
    "age", "sex", "prestroke_mrs", "nihss_baseline", "onset_type",
    "core_ml", "tmax6_ml", "atrial_fib", "center")

# Deliberately excluded from every propensity model; reported in the balance table as negative
# controls so residual imbalance on them is visible.
BALANCE_ONLY: Final[tuple[str, ...]] = (
    "hypertension", "hyperlipidemia", "diabetes", "smoking", "penumbra_ml")

# Bound to the same object rather than repeated, and the two below are computed from it. None may be
# written out as a second literal: invariant 3 exists because two hand-maintained lists drift.
OUTCOME_COVARIATES: Final[tuple[str, ...]] = PS_COVARIATES
STANDARDISATION_COVARIATES: Final[tuple[str, ...]] = tuple(
    c for c in PS_COVARIATES if c != "center")                                  # [§14a]
PS_COVARIATES_FULL: Final[tuple[str, ...]] = PS_COVARIATES + (
    "hypertension", "hyperlipidemia", "diabetes", "smoking")                    # [§13]

# DERIVED_NAMES used to be declared here, as the literal ("onset_type",). It is now computed, in the
# derived-names block at the foot of this file — after OUTCOMES and SUBGROUPS, the two registries it
# reads. It was not deleted; it moved, and it had to move, because this file executes top to bottom
# and both of those registries are declared below this line.


# --- outcome registry [§5] --------------------------------------------------------------------
#
# Exactly eight entries, and exactly one of them is primary: [§8] rests on there being one primary
# quantity with one test. The pilot had fourteen entries and two primaries; it is not lifted.

OUTCOMES: Final[dict[str, Outcome]] = {
    "mrs_90d": Outcome(
        "mRS at 90 days", "ordinal", "primary", None, None, None, False),
    "mrs_0_2_90d": Outcome(
        "mRS 0-2 at 90 days", "binary", "secondary", "mrs_90d", "<=", 2, True),
    "mrs_0_1_90d": Outcome(
        "mRS 0-1 at 90 days", "binary", "secondary", "mrs_90d", "<=", 1, True),
    "tici_2b_3": Outcome(
        "TICI 2b-3", "binary", "secondary", None, None, None, True),
    "sich": Outcome(
        "Symptomatic intracranial haemorrhage", "binary", "safety", None, None, None, False),
    "ph2": Outcome(
        "Parenchymal haematoma type 2", "binary", "safety", None, None, None, False),
    "death_90d": Outcome(
        "Death at 90 days", "binary", "safety", "mrs_90d", "==", 6, False),
    "mrs_5_6_90d": Outcome(
        "mRS 5-6 at 90 days", "binary", "safety", "mrs_90d", ">=", 5, False),
}

# The outcomes Stage 3 builds, as against the ones read directly from the workbook. `source is None`
# is the only test of which is which, here and in Stage 3's loop, so an outcome cannot be derived by
# being remembered. Declared immediately after OUTCOMES, which it reads.
DERIVED_DICHOTOMIES: Final[tuple[str, ...]] = tuple(
    key for key, o in OUTCOMES.items() if o.source is not None)


# --- outcome-model overrides [§8 amendment, DECISION 3] --------------------------------------

OUTCOME_MODEL_OVERRIDES: Final[dict[str, tuple[str, ...]]] = {
    "tici_2b_3": ("center", "atrial_fib")}


def outcome_model_covariates(outcome: str) -> tuple[str, ...]:
    """m_a(X) covariates for `outcome`: the shared §6 set unless it declares an override.

    Treatment is NOT in the returned tuple. Every estimator adds the treatment main effect itself,
    so TICI's reduced model is intercept + treatment + center (2 df) + atrial_fib = 5 parameters
    against 6 non-events, as [§8 amendment] states.

    Every consumer calls this; nobody reads OUTCOME_COVARIATES directly for an outcome model. A
    registry that defaults to the shared entry means an outcome can only diverge by being named.
    """
    if outcome not in OUTCOMES:
        raise KeyError(
            f"{outcome!r} is not in OUTCOMES. Known outcomes: {sorted(OUTCOMES)}. "
            "An outcome model may only be requested for a registered outcome.")
    return OUTCOME_MODEL_OVERRIDES.get(outcome, OUTCOME_COVARIATES)


# --- eligibility [§3, DECISION 1] ---------------------------------------------------------------
#
# There is deliberately no list of contraindication reason strings. Under DECISION 1 the classifier
# is `ivt_contraindicated` plus whether `contraindication_reason` is present, so no string is ever
# matched and there is nothing to enumerate. What Stage 4 asserts instead is that the flag is 0/1,
# never missing, and never 1 for a treated patient [Stage 4 §4.2].
#
# The three classes are declared individually because eligibility.py writes them by name. Indexing
# into ELIGIBILITY_ORDER would couple the classification rule to a display order, so reordering a
# table's columns would rewrite the rule; a literal in eligibility.py would give each class a second
# declaration. test_config.py asserts the order below holds exactly these three.
ELIGIBLE: Final[str] = "eligible"
INDETERMINATE: Final[str] = "indeterminate"
INELIGIBLE: Final[str] = "ineligible"

# Display and table order [Stage 4 §7], computed so the three above are declared once.
ELIGIBILITY_ORDER: Final[tuple[str, ...]] = (ELIGIBLE, INDETERMINATE, INELIGIBLE)

# The classes a patient is RETAINED on [§3]. Decided by the PI on 2026-08-10: [§3] and [§14a] use
# one rule, and only ineligible patients are ever dropped. This is 43 of 126 records — 80% of the
# primary cohort's control arm — so `== ELIGIBLE` and `!= INELIGIBLE` meaning different things in
# two stages would not look like a bug in either. eligibility.retained() is the only reader.
ELIGIBILITY_RETAINED: Final[tuple[str, ...]] = (ELIGIBLE, INDETERMINATE)

# The column Stage 4 produces. Named here for the reason TREATMENT is: every consumer refers to the
# column through this constant, and it is deliberately NOT in DERIVED_NAMES [Stage 4 §8].
ELIGIBILITY: Final[str] = "eligibility"


# --- post-time-zero denylist [invariant 4] ----------------------------------------------------
#
# Built from the outcome registry so a new outcome is denylisted by existing, not by being
# remembered. Defined after OUTCOMES; the ordering is load-bearing.
POST_TIME_ZERO: Final[frozenset[str]] = frozenset(
    {"onset_to_ivt_min", "onset_to_groin_min", *OUTCOMES})


# --- structural non-applicability and plausible ranges ---------------------------------------

STRUCTURALLY_NON_APPLICABLE: Final[dict[str, str]] = {
    "onset_to_ivt_min": "control arm — no IVT was given, so the time does not exist [§11]",
}

# Absence that carries information rather than data loss. `contraindication_reason` is missing for 82
# of 126 records, and that absence is exactly the bit [DECISION 1] reads: no reason was recorded,
# which [§3] forbids reading as "no contraindication" and Stage 4 maps to `indeterminate`. Listing
# those 82 beside genuine data loss misrepresents both. Stage 2 labels them; nothing imputes them.
#
# Disjoint from STRUCTURALLY_NON_APPLICABLE, which test_config.py asserts: Stage 2's missingness
# table has one branch per constant, so an overlap would resolve by branch order rather than by
# declaration — on the column that carries the [DECISION 1] bit.
INFORMATIVE_ABSENCE: Final[dict[str, str]] = {
    "contraindication_reason":
        "no reason was recorded — the one bit [DECISION 1] reads [§3]",
}

# (low, high), inclusive; None as an upper bound means unbounded above. For Stage 2's assertions.
PLAUSIBLE_RANGES: Final[dict[str, tuple[int, int | None]]] = {
    "age":             (18, 110),
    "nihss_baseline":  (0, 42),
    "prestroke_mrs":   (0, 5),       # 6 is death; impossible pre-stroke
    "mrs_90d":         (0, 6),
    "core_ml":         (0, None),
    "tmax6_ml":        (0, None),
    "penumbra_ml":     (0, None),
    "onset_to_ivt_min":   (0, None),
    "onset_to_groin_min": (0, None),
}

# Stage 2 asserts every value in these columns is in {0, 1} or missing. Two matter more than the
# rest: `ivt` is the exposure, and a 2 or a -1 there would flow into every weight, every estimate
# and every arm label without a single error; `ivt_contraindicated` is the [§3] classifier.
BINARY_COLUMNS: Final[tuple[str, ...]] = (
    "sex", "hypertension", "hyperlipidemia", "diabetes", "smoking", "atrial_fib",
    "wake_up", "unwitnessed", "ivt", "sich", "ph2", "tici_2b_3", "ivt_contraindicated")


# --- labels pending a data query ---------------------------------------------------------------
#
# §10 records that the `sex` coding is undocumented (0/1, 74 vs 52). Adjustment is unaffected, but
# the baseline table cannot label the levels. Stage 14 must print `sex = 0` / `sex = 1` while this
# is None, and must never invent labels.
SEX_LABELS: dict[int, str] | None = None


# --- subgroups [§13] ----------------------------------------------------------------------------

MRS_THRESHOLDS: Final[tuple[int, ...]] = (0, 1, 2, 3, 4, 5)      # cumulative RD_k [§8]

# The [§13] subgroups, after the amendment of 2026-08-10 withdrew the third one. That amendment is
# in statistical_analysis_plan.md §13 and its consequences for the code are in Stage 3 §6.1; the two
# constants that served only the withdrawn subgroup were declared here and are deleted rather than
# commented out. A withdrawn subgroup with a live constant is an invitation to restore it without
# the amendment, and Stage 3's definition of done greps this repository's Python for their names.
#
# Neither key is a covariate: test_config.py asserts no subgroup name appears in any covariate list,
# so a subgroup cannot become an adjustment variable by being convenient. Neither is post-time-zero
# either — both are baseline attributes — so neither joins POST_TIME_ZERO.
SUBGROUPS: Final[dict[str, str]] = {
    "unknown_onset":     "unwitnessed or wake-up onset versus witnessed [§13]",
    "core_above_median": "core volume above the cohort median [§13]",
}

# The subgroups whose definition needs the whole cohort rather than one record. Stage 3's second
# entry point, derive_cohort, computes these after the [§3] restrictions and raises if asked twice:
# the median is a property of the cohort, so freezing it as a number here would silently decouple it
# from the cohort it describes, and recomputing it per bootstrap replicate would give every replicate
# its own cut-point and its own subgroup [Stage 3 §6.4].
COHORT_DEPENDENT_SUBGROUPS: Final[tuple[str, ...]] = ("core_above_median",)


# --- derived names ------------------------------------------------------------------------------
#
# Produced by Stage 3, so they are analysis names that the contract cannot supply. Declared here, at
# the foot of the file, because both are computed from registries above: DERIVED_DICHOTOMIES from
# OUTCOMES, and the subgroup keys from SUBGROUPS. Written in the covariates block where
# DERIVED_NAMES used to sit, they would raise NameError on import.
#
# "onset_type" is the one literal, because it is the only derived name no other registry declares.
# test_config.py asserts DERIVED_NAMES is disjoint from ANALYSIS_NAMES, so a derived column can never
# shadow one the contract delivers.

DERIVED_NAMES: Final[tuple[str, ...]] = ("onset_type", *DERIVED_DICHOTOMIES, *SUBGROUPS)

# The derived columns that exist after Stage 3's first, row-wise entry point. `core_above_median` is
# not among them — it does not exist until derive_cohort has run — which is why Stage 3's derived
# missingness table ranges over this and not over DERIVED_NAMES [Stage 3 §8.3].
ROW_WISE_DERIVED: Final[tuple[str, ...]] = tuple(
    name for name in DERIVED_NAMES if name not in COHORT_DEPENDENT_SUBGROUPS)


# --- the no-raw-names rule ---------------------------------------------------------------------
#
# "Nothing outside this module may reference a raw column name" is tested, not left as a convention.
# read_excel returns the whole sheet, so 'Deathat90days', the trailing-space mRS 5-6 column, the
# administrative status field and the free text of the contraindication column are all sitting in
# the DataFrame with no analysis name. Nothing but this rule stops a later module from reading the
# shipped death flag — which would move the primary safety outcome from 25 events to 26 — or from
# reading the free text that DECISION 1 rules out.
#
# The scan itself lives in test_config.py; this module owns only the exemptions. The comparison is
# exact and case-sensitive, and must stay that way: TICI_2b_3 -> tici_2b_3 and Center -> center
# differ from their own analysis names only by case, so a case-insensitive scan would flag every
# legitimate downstream use of them.
EXEMPT_FROM_RAW_NAME_SCAN: Final[dict[str, str]] = {
    "config.py":                "owns the contract",
    "stage0_data_inventory.py": "predates the configuration; reads raw headers by design",
    "test_config.py":           "must construct malformed headers to test the contract check",
}


# --- the contract check -------------------------------------------------------------------------

_HEADER_ONLY_ABSENT = "Unnamed: 41"


def _whitespace_note(contract_key: str) -> str:
    """How `contract_key` differs from its own stripped form, for the strip-match hint."""
    if contract_key.rstrip() != contract_key and contract_key.lstrip() == contract_key:
        return "note trailing space"
    if contract_key.lstrip() != contract_key and contract_key.rstrip() == contract_key:
        return "note leading space"
    if contract_key.strip() != contract_key:
        return "note surrounding whitespace"
    return "they differ only in whitespace"


def assert_column_contract(raw_columns: Iterable[str]) -> None:
    """Raise SchemaError on any symmetric difference between the contract and the workbook.

    Must be called against a **full** read, never a header-only read: the workbook's last column is
    headerless, and pandas materialises it only once rows are read.
    """
    cols = list(raw_columns)      # materialise, so a generator can be passed and still be named
    contract = set(COLUMN_CONTRACT)

    if len(cols) == N_COLUMNS_EXPECTED - 1 and set(cols) == contract - {_HEADER_ONLY_ABSENT}:
        raise SchemaError(
            f"{len(cols)} columns — this is the signature of a header-only read (nrows=0). "
            "pandas materialises the headerless final column only once rows are read. "
            "Re-read with rows.")

    missing = contract - set(cols)
    extra = set(cols) - contract
    if not missing and not extra:
        return None

    lines = [
        f"the workbook's columns do not match COLUMN_CONTRACT "
        f"({len(cols)} given, {N_COLUMNS_EXPECTED} expected)"]
    if missing:
        lines.append("  in the contract, absent from the workbook: "
                     + ", ".join(repr(m) for m in sorted(missing)))
    if extra:
        lines.append("  in the workbook, absent from the contract: "
                     + ", ".join(repr(e) for e in sorted(extra)))
    for e in sorted(extra):
        for key in sorted(k for k in contract if k.strip() == e.strip()):
            lines.append(f"  {e!r} is not in the contract; did you mean {key!r} "
                         f"({_whitespace_note(key)})?")
    raise SchemaError("\n".join(lines))
