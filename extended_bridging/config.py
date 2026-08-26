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
``math`` is standard library and ``ci_min_draws`` is the one function that needs it; the Stage 10
spec's fence writes ``np.ceil`` there, and ``math.ceil`` is the same value without making the import
of this module depend on numpy [Stage 10 §13].
"""
from __future__ import annotations

import math
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


# --- [§10] inference: the bootstrap's definitions ----------------------------------------------
#
# The [§10] stratification variable. A NAME and not the column itself, because `bootstrap.resample`
# is general in its stratum (Stage 10 §12.2): [§14a] and [§14b] resample the same way over a
# population that includes a centre with no treated patients, and a function that hard-coded
# "center" would still be right there while being right for the wrong reason.
BOOT_STRATUM: Final[str] = "center"

# The [§10] confidence level. 0.95 is [§10]'s own "percentile 95% confidence intervals"; it is here
# rather than as a literal so that `ci_min_draws` below can be derived from it and cannot disagree.
CI_LEVEL: Final[float] = 0.95

# The percentile DEFINITION, and this constant changes an ANSWER rather than a tolerance
# [Stage 10 §8.2]. numpy's default is "linear", which interpolates between order statistics, and
# [§10] claims its p-value "is the smallest level at which the percentile interval for beta excludes
# the null" and "agrees by construction with the reported interval". THAT CLAIM IS FALSE UNDER THE
# DEFAULT: measured over 30000 constructed draw sets at B = 2000, "linear" disagrees with the
# p-value in 4.707% of them, and EVERY disagreement is at a smaller-tail count of exactly 50, where
# p is exactly 0.0500 and an interpolated limit lands on whichever side of zero the arithmetic falls.
#
# "inverted_cdf" is the empirical-CDF quantile -- the smallest order statistic whose cumulative
# proportion reaches the level -- which is the SAME object Pr(beta* <= 0) is computed from. That
# identity is why the two agree, and it is why this is not "lower", which agrees at these levels by
# arithmetic coincidence rather than by estimating the same quantile.
#
# Measured: on this workbook's own 2000 draws all three methods give the same verdict for beta and
# for all seven risk differences, so the cohort does not witness this and the pin is justified by
# construction. Prespecified for FIRTH_*'s reason: [§10] refits in every one of N_BOOT replicates.
PERCENTILE_METHOD: Final[str] = "inverted_cdf"

# The FitError buckets Stage 8 §11 requires reported separately, keyed by the leading token of the
# raised message [Stage 10 §7.2]. `model.FitError` carries no code -- it is a bare RuntimeError
# subclass (model.py:89) -- and every raise site IDENTIFIES itself by that token, so classification
# is textual and test_bootstrap.py §15.6 SCANS the modules and asserts every token found is a key
# here. A reworded message is then a test failure rather than a counter that silently reads zero.
#
# "polr:" covers two failures -- step-halving exhausted and no convergence in POLR_MAX_ITER -- and
# they share a bucket. Stage 8 §11 asks for G7 against everything else, and Stage 8 §5.2 records
# that exhausting the halvings is not a convergence route; both measured 0 in N_BOOT replicates.
#
# **THERE ARE NINETEEN RAISE SITES AND NOT SIXTEEN, AND THE SCAN IS WHAT FOUND THE OTHER THREE.**
# Stage 10 §7.2, §13 and §15.6 all say sixteen -- fourteen in model.py and two in outcome.py -- and
# re-scanning by AST reproduces those sixteen exactly. What that scan's SCOPE omitted is
# `propensity.py`, which raises `model.FitError` three more times: twice from `ess` with the token
# `ESS:` and once from `_assert_probabilities` with `F5`. Both are on `propensity.fit`'s own path,
# which is the path Stage 10 §7.1 says fails a WHOLE replicate, so both are reachable from inside a
# replicate and would have hit `_bucket`'s unrecognised-token raise -- turning a droppable sparse
# replicate into a crash.
#
# **AND F5 IS NOT HYPOTHETICAL: IT IS THE ONLY `FitError` THE WORKBOOK ITSELF PRODUCES.** Stage 10
# §7.5 measures `propensity.fit` raising `FitError` in 2 of N_BOOT replicates and does not say which
# token. Measured over the [§10] draw at C.SEED: both are **F5** -- a fitted probability on the
# boundary -- and neither is any of the sixteen tokens §7.2 scanned for. Without these two entries the
# prespecified run terminates on replicate ~700 of 2000 with an unrecognised-token SchemaError.
#
# Both are bucketed `degenerate_design`, which keeps this map's VALUES the four §7.2 names that
# test_config.py pins. F5 is a fitted probability on the boundary, which Stage 6 §5.5 calls "a
# degenerate fit and not a confident one"; `ESS:` is an arm with no weighted patient or with every
# weight at zero, which Stage 6 §3.1 calls a structural non-positivity. Neither is a separation
# guard and neither is a convergence route, so the two remaining buckets would both be wrong.
FAILURE_BUCKETS: Final[dict[str, str]] = {
    "G7": "separation",
    "S8": "constant_outcome",
    "G6": "degenerate_design",
    "polr:": "nonconvergence",
    "Firth:": "nonconvergence",
    "O1": "degenerate_design", "O2": "degenerate_design", "O3": "degenerate_design",
    "O4": "degenerate_design", "O5": "degenerate_design", "O6": "degenerate_design",
    "F3": "degenerate_design", "F4": "degenerate_design", "F6": "degenerate_design",
    "F5": "degenerate_design", "ESS:": "degenerate_design",
}


def ci_min_draws(level: float = CI_LEVEL) -> int:
    """The fewest draws at which a `level` percentile limit is an order statistic at all.

    DERIVED, not chosen. Under PERCENTILE_METHOD the lower limit is order statistic
    ceil((1-level)/2 * n), one-based; for that index to exceed 1 -- for the limit to be interior
    rather than the sample minimum -- n must exceed 2/(1-level), which is 40 at CI_LEVEL = 0.95.
    **AT EXACTLY THE FLOOR THE LOWER LIMIT IS THE SAMPLE MINIMUM AND THE UPPER IS THE SECOND-LARGEST
    DRAW, NOT THE MAXIMUM**, and the asymmetry is the point: ceil(0.025 * 40) = 1 gives index 0
    while ceil(0.975 * 40) = 39 gives index 38. Measured on np.arange(40.0): (0.0, 38.0). The floor
    is a statement about the LOWER limit, which is the one that stops carrying information first;
    below it np.percentile returns the minimum while still calling it a percentile
    [Stage 10 §8.3, §15.8].

    A function rather than a constant so it cannot disagree with CI_LEVEL. N_BOOT = 2000 against a
    floor of 40 means an estimand needs 98% of its replicates to fail before it loses its interval:
    measured, the worst per-outcome drop rate on this cohort is 0.8%, so the branch is unreachable
    on v7 and is specified anyway.

    THE SNAP IS NOT DEFENSIVE PROGRAMMING AND `int(ceil(2.0 / (1.0 - level)))` IS WRONG.
    `1.0 - 0.90` is 0.09999999999999998, so that expression returns 21 where 20 is intended. At
    CI_LEVEL = 0.95 it happens to return the intended 40 -- `1.0 - 0.95` errs the other way and the
    quotient is 39.99999999999996 -- so the bug is invisible at the only level this study uses and
    appears at the first [§13] sensitivity level anybody tries. A level whose quotient is genuinely
    non-integral still takes the ceiling: at 0.93 the quotient is 28.571 and the answer is 29.

    `math.ceil` and not `np.ceil`, which is what Stage 10 §13's fence writes. The two return the
    same integer here and this module imports only the standard library [Stage 1 §3], which is a
    property Stage 10 §13 does not list among what it amends.
    """
    exact = 2.0 / (1.0 - level)
    nearest = round(exact)
    return int(nearest if abs(exact - nearest) < 1e-9 else math.ceil(exact))

# Named for the minority cell, min(events, non-events), not for events [§8 amendment, DECISION 3].
# A constant named RARE_EVENT_THRESHOLD would invite the exact misreading the amendment corrects:
# TICI 2b-3 has 114 events and 6 non-events, so it passes an event-count rule while failing that
# rule's rationale.
RARE_MINORITY_THRESHOLD: Final[int] = 10

# The [§8] weighted marginal odds ratio's continuity correction [Stage 9 §6.3]. Haldane-Anscombe:
# added to all four weighted pseudo-counts, and ONLY when a weighted proportion reaches 0 or 1, so an
# interior estimate is bit-for-bit uncorrected. PI decision, 2026-08-24 — neither [§8] nor any
# amendment specifies a rule and the roadmap says only "kept finite".
#
# Prespecified for FIRTH_*'s reason (config.py:315-317, which this insertion moved down from the
# 281-284 the Stage 9 spec cites): [§10] refits this in every one of N_BOOT replicates, so a
# correction that changes an ANSWER is a property of the sampling distribution and not a runtime knob.
#
# The scale is stated rather than implied, and PER ARM rather than averaged. These are WEIGHT-SUMS,
# not row counts: Sw = 27.736623 over the ATO population, splitting 13.626360 treated / 14.110263
# control against row counts of 39 treated / 53 control. The conventional 0.5 is calibrated against
# ROWS, so against these weight-sums it is 2.9x more aggressive in the treated arm and 3.8x in the
# control arm — measured per arm, because averaging the two row counts to "about 46" and reporting
# "about three times" hides that the two arms are corrected by different relative amounts, which is
# exactly the asymmetry that lets the correction cross the null [Stage 9 §6.4].
#
# It is NOT scaled to Sw to compensate, because a correction whose magnitude is a function of the
# weights is a correction whose magnitude is a function of the propensity model, and [§10] refits
# that in every replicate.
#
# CONSEQUENCE, MEASURED, AND IT IS NOT A ROUNDING EFFECT: because the arms are corrected unequally,
# an empty cell in the HEAVIER arm can carry the odds ratio ACROSS 1 rather than toward it — e.g.
# Sw1 = 17.007 against Sw0 = 22.768 with p0 = 0.00647 returns 1.0201 where the uncorrected value is
# 0.0. About 3 replicates in N_BOOT = 2000 report a safety outcome with no bridging events as
# favouring EVT alone. PI-reversible; [Stage 9 §6.4, §17].
#
# Measured: the branch is unreachable on v7 — no arm holds an empty cell for any of the seven
# outcomes — and fires in 17.0% of stratified replicates for sich, 13.0% for tici_2b_3 and 2.3% for
# ph2. So this constant is materially a choice about two INTERVALS and barely a choice about any point
# estimate [Stage 9 §6.4].
OR_CONTINUITY: Final[float] = 0.5


# --- the Firth fit [§7, Stage 6 §5] -----------------------------------------------------------
#
# Prespecified, and they are part of the estimator: [§10] refits this model in every one of N_BOOT
# replicates, so a tolerance is a property of the sampling distribution and not a runtime knob. The
# ORDER the two convergence tests run in is prespecified for the same reason [Stage 6 §5.3].
#
# These are the PILOTS' values, carried from pilots/analysis.py. An earlier draft of this block called
# them the R reference implementation's defaults. **That claim is struck, not merely withdrawn**: the
# reference implementation's source has since been read (DoD-15a, recorded in Stage 6 §18g) and its
# defaults are none of these — and, more importantly, its stopping RULE is a different rule, testing
# the likelihood change, the score AND the coefficient step conjunctively where §5.3 tests the first
# two disjunctively and excludes the third deliberately.
#
# §5.3's argument for excluding the coefficient step is untouched by that and never depended on the
# attribution: under near-separation the penalised surface is genuinely flat, so a step-norm criterion
# would report a finite correct fit as a failure — and [§10] drops failed replicates, so the ones
# dropped would be exactly the sparse ones. That is selection on the replicate.
#
# The R package's name is deliberately not written in any shipped module [Stage 6 DoD-15]; it lives in
# the spec and in the one test module that knows R exists.
#
# The prefix is FIRTH_ and not PS_ because model.py reads them, and model.py does not know what the
# exposure is [Stage 6 §0.1, §12.12] — whose attribute scan forbids any `C.` name beginning `PS_`
# there. A PS_-prefixed tolerance would be a propensity-flavoured constant read by an outcome fit.

FIRTH_MAX_ITER: Final[int] = 200          # measured: 7 on the cohort, 11 on the hand frame
FIRTH_TOL: Final[float] = 1e-8            # |Δ penalised log-likelihood| — the route that always fires
FIRTH_SCORE_TOL: Final[float] = 1e-6      # max |modified score|; fires on ill-conditioned I — §5.3
FIRTH_MAX_HALVINGS: Final[int] = 30       # exhausting them raises; it is not a convergence route
FIRTH_MAX_STEP: Final[float] = 5.0        # trust radius, RELATIVE to ‖beta‖ — see Stage 6 §5.2a
FIRTH_ETA_CLIP: Final[float] = 500.0      # keeps exp() in range; reachable at the upper end [§3.1]
FIRTH_WEIGHT_FLOOR: Final[float] = 1e-10  # floors p(1-p) so the information matrix stays invertible


# --- the weighted proportional-odds fit [§8, Stage 8 §5] ---------------------------------------
#
# Prespecified for FIRTH_*'s reason: [§10] refits this model in every one of N_BOOT replicates, so a
# tolerance is a property of the sampling distribution and not a runtime knob, and the ORDER the two
# convergence tests run in is prespecified for the same reason [Stage 8 §5.2].
#
# POLR_MAX_ABS_BETA is not a tolerance and is the one constant here that changes an ANSWER: a fit
# reaching it raises, and [§10] drops and counts the replicate. Separation in this estimator does not
# present as non-convergence — measured, it converges in 17 iterations on the score criterion with
# every safeguard counter at zero and returns exp(beta) = 6.5e15 — so this bound is the only thing
# that turns a degenerate fit into a countable failure [Stage 8 §6].
#
# Its value is chosen from a MEASURED SPARSE REGION and not from any estimate. Over 4800 fits of a
# 92-record seven-category frame at twelve true effect sizes from 0 to 6, the largest non-degenerate
# |beta| was 8.7873 and the smallest degenerate one 18.8055. Re-measured across every cutpoint count a
# [§10] replicate can produce — 2 to 6 cutpoints, 1800 fits each — the legitimate maximum reaches
# 11.04 and the degenerate minimum falls to 18.98, so the safe interval is (11.04, 18.98) and the
# region between the modes is SPARSE rather than empty. 14.0 sits inside it with 27% clearance below
# and 26% above, and every bound from 12 to 20 drops exactly the same replicates — so 14.0 is a point
# on a plateau and not a tuned value. exp(14) = 1.2e6, five orders of magnitude above any reportable
# stroke odds ratio, so a rejected fit is one whose point estimate no manuscript could print. An
# earlier draft used 10.0, calibrated on 500 fits at a single true effect of 0.5; it is now positively
# excluded, because it rejects the legitimate fits at 10.01 to 11.04 that appear at odd cutpoint
# counts [Stage 8 §6.3, §20].
#
# The prefix is POLR_ and not OUTCOME_ because model.py reads them and model.py does not know what the
# outcome is [Stage 6 §0.1, Stage 8 §14.13].

POLR_MAX_ITER: Final[int] = 200          # measured: 4 to 5 on every frame in Stage 8 §20
POLR_TOL: Final[float] = 1e-8            # |Δ weighted log-likelihood|
POLR_SCORE_TOL: Final[float] = 1e-6      # max |score|; the route a SEPARATED fit returns on
POLR_MAX_HALVINGS: Final[int] = 30       # exhausting them raises; not a convergence route
POLR_MAX_STEP: Final[float] = 5.0        # trust radius, RELATIVE to ‖par‖ [Stage 6 §5.2a]
POLR_ETA_CLIP: Final[float] = 500.0      # keeps exp() in range; measured never approached
POLR_MAX_ABS_BETA: Final[float] = 14.0   # the separation guard [Stage 8 §6.3]


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

# [§6] names four vascular risk factors as balance negative controls and adds them back in the
# [§13] full-covariate sensitivity propensity model. Those are the same four, so the set is
# COMPUTED from that identity rather than declared beside it: a fifth risk factor added to
# PS_COVARIATES_FULL becomes a negative control here in the same edit, and cannot fail to.
NEGATIVE_CONTROLS: Final[tuple[str, ...]] = tuple(
    c for c in PS_COVARIATES_FULL if c not in PS_COVARIATES)

# The [§9] balance set: the full [§6] confounder set plus everything BALANCE_ONLY carries, which
# is the four negative controls and penumbra_ml. NOT the propensity model's covariate list —
# [§9] judges balance against the full set regardless of what a specification fitted.
BALANCE_SET: Final[tuple[str, ...]] = PS_COVARIATES + BALANCE_ONLY

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

# The [§5] primary outcome, computed from the registry rather than named a second time. OUTCOMES'
# own comment states that exactly one entry is primary because [§8] rests on there being one primary
# quantity with one test; this is that sentence made readable by code, so that outcome.py cannot
# write the key as a literal and cannot drift from the registry. Declared immediately after
# DERIVED_DICHOTOMIES, which is the first point at which OUTCOMES exists.
#
# `next` over a generator takes the FIRST primary entry, so it would be silent about a second one.
# test_config.py asserts there is exactly one, which is the assertion that notices [Stage 8 §12, T1].
PRIMARY_OUTCOME: Final[str] = next(
    key for key, o in OUTCOMES.items() if o.family == "primary")

# The [§5] binary outcomes, computed from the registry rather than named a second time. PRIMARY_OUTCOME
# above takes the one primary entry; this takes the seven binary ones, and test_config.py asserts the
# two partition OUTCOMES exactly — so an outcome cannot be estimated by neither Stage 8 nor Stage 9,
# and cannot be estimated by both. Declared immediately after PRIMARY_OUTCOME for that reason: the two
# are one decision about the registry and a reader checking the partition should find them together.
BINARY_OUTCOMES: Final[tuple[str, ...]] = tuple(
    key for key, o in OUTCOMES.items() if o.kind == "binary")


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


# --- eligibility [§3, DECISION 1a] --------------------------------------------------------------
#
# There is deliberately no list of contraindication reason strings. Under DECISION 1a (PI,
# 2026-08-13) the classifier is `ivt_contraindicated` ALONE: flag = 1 -> ineligible, flag = 0 ->
# eligible. No string is ever matched and — unlike under DECISION 1 — neither is the presence of one,
# so `contraindication_reason` feeds no class at all and there is nothing to enumerate.
#
# What Stage 4 asserts instead is that the flag is 0/1, never missing, and never 1 for a treated
# patient [Stage 4 §4.2]. That third assertion is **load-bearing** under DECISION 1a where it was
# belt-and-braces under DECISION 1: it is what makes the revealed-fact rule redundant, and with that
# rule deleted it is the only thing standing between a treated-and-flagged record and a contradiction
# absorbed without a trace.
#
# The two classes are declared individually because eligibility.py writes them by name. Indexing into
# ELIGIBILITY_ORDER would couple the classification rule to a display order, so reordering a table's
# columns would rewrite the rule; a literal in eligibility.py would give each class a second
# declaration. test_config.py asserts the order below holds exactly these two.
ELIGIBLE: Final[str] = "eligible"
INELIGIBLE: Final[str] = "ineligible"

# Display and table order [Stage 4 §7], computed so the two above are declared once.
ELIGIBILITY_ORDER: Final[tuple[str, ...]] = (ELIGIBLE, INELIGIBLE)

# The classes a patient is RETAINED on [§3]. Under DECISION 1a eligibility is two-valued, so
# `== ELIGIBLE` and `!= INELIGIBLE` coincide and the 43-record gap DECISION 1's third class opened is
# closed — the PI's decision of 2026-08-10 turns out to have chosen the same population by a longer
# route. The tuple stays for two reasons that survive the amendment: it is the one seam a third class
# would come back through (a [§13] sensitivity arm over the undocumented controls would reopen
# exactly that gap), and a predicate reading a registry rather than writing a comparison costs
# nothing. eligibility.retained() is the only reader.
ELIGIBILITY_RETAINED: Final[tuple[str, ...]] = (ELIGIBLE,)

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

# The declared mRS level set, computed from the plausible range that already declares it, so the
# level set and the range cannot disagree. MRS_THRESHOLDS below is the K = J - 1 cutpoints and
# test_config.py asserts MRS_LEVELS[:-1] == MRS_THRESHOLDS: RD_k stops at 5 because P(Y <= 6) is 1 in
# both arms by definition [Stage 8 §8.3], and a seventh threshold would be a structural zero.
#
# MRS_THRESHOLDS sits in this [§13] subgroups block, which is where it was declared and is not where
# it belongs. It is deliberately NOT moved: six test modules import it, and the assertion above makes
# the misfiling harmless [Stage 8 §12, §17].
MRS_LEVELS: Final[tuple[int, ...]] = tuple(
    range(PLAUSIBLE_RANGES["mrs_90d"][0], PLAUSIBLE_RANGES["mrs_90d"][1] + 1))

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
