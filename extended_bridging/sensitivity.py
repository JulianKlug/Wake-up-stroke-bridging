"""Stage 11 — [§13] multiplicity, subgroups, the E-value and the full-covariate sensitivity arm.

Four deliverables and no fifth. [§13] whole — its multiplicity clause, its subgroup clause and the
one propensity-specification sensitivity row DECISION 4 promoted from deferred — plus [§6]'s E-value,
which is prescribed one section away and is reported beside the primary estimate the [§13] arm is
compared against.

Section references in brackets are to ``statistical_analysis_plan.md``. The specification for this
module is ``specs/stage11_multiplicity_subgroups_evalue.md``; nothing here is invented outside it.

::

    propensity.fit(cohort, a)  -> Propensity: e, w, in_model, spec       [Stage 6 §11]
    balance.assess(...)        -> Balance:    19 rows, roles off ps.spec [Stage 7 §11]
    outcome.primary(...)       -> Primary:    beta, alpha, rd, cumulative [Stage 8 §11]
    outcome.secondary(...)     -> Secondary:  seven BinaryEstimate        [Stage 9 §11]
    bootstrap.run(...)         -> Bootstrap:  26 draws, intervals         [Stage 10 §11]
                       │
                       ▼
    ┌──────────────────────────────────────────────────────────────────────────────────┐
    │  multiplicity(sec, boot, audit)        -> Multiplicity                  (§6)     │
    │    families = sec.by_family()          3 and 4, NEVER off `intervals`   (§6.1)   │
    │    p        = boot.intervals[k].p      READ, never recomputed           (§6.1)   │
    │    benjamini_hochberg(p) per family    step-up + descending running min (§6.2)   │
    │                                                                                  │
    │  e_value_primary(est, bal, boot, a)    -> EValue                        (§7)     │
    │    RR = sqrt(est.odds_ratio)           the approximation is a FIELD     (§7.1)   │
    │    spans = lo <= 0 <= hi               three routes, one answer         (§7.2)   │
    │                                                                                  │
    │  subgroups(cohort, ps, est, audit)     -> Subgroups                     (§8)     │
    │    per S: polr on (A, S, A x S), ps.w, NO [§6] covariate                (§8.2)   │
    │    G8 before `polr`, G9 after it; H8 when the POINT ESTIMATE raises     (§8.5)   │
    │    replicates(...) 6 keys, TWO fits per replicate -> NO Diagnostics     (§8.5)   │
    │                                                                                  │
    │  full_covariate(cohort, ps, audit)     -> Arm       [§13, DECISION 4]   (§5)     │
    │    propensity.fit_full -> balance.assess -> outcome.primary (UNCHANGED) (§5.1)   │
    │    replicates(...) 7 keys, its own body, NOT `run` -> Diagnostics       (§5.4)   │
    │    _intervals(draws, lambda k: False)  no p on ANY arm key              (§5.5)   │
    └──────────────────────────────────────────────────────────────────────────────────┘
                       │
                       ▼
      Multiplicity, EValue, Subgroups, Arm  →  Stage 14

**Why this is one module and not two, and the argument against it is real** [Stage 11 §0.1].
``multiplicity`` and ``e_value_primary`` fit nothing: they are pure functions of numbers Stages 9 and
10 already produced. ``subgroups`` and ``full_covariate`` refit, and they import ``propensity``,
``balance`` and ``model``. A split would put two of the four deliverables in a module testable with
no frame at all, and **the cost of declining it is named rather than denied**: ``benjamini_hochberg``
and ``e_value`` are pure arithmetic living in a module that imports ``propensity``, so a test of the
step-up drags the whole estimator stack in behind it at collection time. It is declined because [§13]
is ONE section of the plan whose three clauses share one subject — what is reported alongside the
primary estimate, and under what qualification — and a boundary drawn between "fits something" and
"does not" would split [§13] down an axis [§13] does not have. The trigger to revisit it is a second
consumer of ``benjamini_hochberg``, which would mean a family partition existing outside [§13].

**What this module does not do** [Stage 11 §0.2]. It refits nothing Stage 10 refit and recomputes no
interval, no limit and no p-value: Stage 10 §14 forbids it by name, and ``multiplicity`` turns that
into a precondition rather than a convention. It adds no estimator — ``bootstrap.resample``,
``bootstrap.replicates``, ``bootstrap.percentile_ci`` and ``bootstrap.bootstrap_p`` are called and
none is reimplemented, and the one fit written here is ``model.polr`` over a three-column design,
which is the landed ordinal fitter with one column added. It does not touch the primary estimate or
the [§7] specification. And it does not decide the remaining five [§13] rows: the seam is a registry
of two declared specifications in ``config.py``, not a covariate-list argument.

**No number this module produces appears in any document under ``specs/``.** Every adjusted p, every
E-value, every subgroup odds ratio and every arm limit exists in the gitignored audit log and nowhere
else in version control [Stage 11 §1].

This module is **not** exempt from the Stage 1 §7 raw-name scan and must never become exempt.
"""
from __future__ import annotations

import math
from collections.abc import Callable
from dataclasses import dataclass

import numpy as np
import pandas as pd

import balance
import bootstrap
import config as C
import model
import outcome
import propensity
# `_fmt` is private to data.py and is imported anyway, for the reason derive.py, propensity.py,
# balance.py, outcome.py and bootstrap.py give: it is the pipeline's *one* float formatter, and a
# second one is a second way for two runs to disagree.
from data import Audit, _fmt


# The `RR ~ sqrt(OR)` statement, as a constant, so that `EValue.approximation` carries one string and
# no print site can carry a second [Stage 11 §3.1, §7.5, §19c]. The citation is VanderWeele's and not
# this plan's: the conversion is his, and attributing it to the roadmap overstates whose choice it is.
_APPROXIMATION: str = (
    "RR ~ sqrt(OR) for a common outcome [VanderWeele 2017, Epidemiology 28(6):e58; roadmap "
    "Stage 11]. Licensed by a PREVALENCE and not by the outcome being ordinal: the reference "
    "implementation applies the square root only above 15% at end of follow-up, and the primary "
    "estimand is a cumulative-odds shift ASSUMED CONSTANT across all six [§8] thresholds, each with "
    "its own baseline risk. There is no published bounding factor for a common odds ratio from a "
    "proportional-odds model, so this is a heuristic on an approximated scale and not the bound the "
    "derivation licenses.")

_MEASURE: str = "common odds ratio [§8]"

# The audit step names. Constants rather than literals at the `audit.record` call sites, for
# `outcome.py`'s reason (outcome.py, the three step bases): the acceptance suite asserts each entry
# appears under this name, and a test comparing against a string literal it also writes is a test of
# nothing.
_STEP_ARM: str = "sensitivity_arm"
_STEP_SUBGROUP_FIT: str = "subgroup_fit"
_STEP_SUBGROUP_REPLICATES: str = "subgroup_replicates"
_STEP_MULTIPLICITY: str = "multiplicity"
_STEP_EVALUE: str = "e_value_primary"

# The three quantities each subgroup contributes to §3.3's key set, in the order the replicate body
# fills them. Computed into keys from `C.SUBGROUPS` and never listed as six literals, which is
# Stage 9 §4.1's move two stages on: a third subgroup restored to the registry gains its three keys,
# its interval and its p in the same edit, and cannot fail to.
_SUBGROUP_QUANTITIES: tuple[str, ...] = ("beta", "gamma", "beta_plus_gamma")


# --- what the stage returns [§3.1] -----------------------------------------------------------------

@dataclass(frozen=True)
class FamilyCorrection:
    """One [§13] family's Benjamini-Hochberg result, and the counters that explain its denominator.

    **`m_declared` and `m_used` are BOTH fields and neither is derivable from the other**, which is
    the whole of §6.5: `m_declared` is the family's size in the [§5] registry and `m_used` is the
    number of p-values that entered the procedure. A reader given only the second cannot tell a
    complete family from a thinned one, and BH's scaling factor IS `m_used`, so the two numbers
    together are what makes the adjusted p checkable by hand.

    `absent` names the outcomes that had no p, with the reason, because "we tested this outcome and
    it had no interval" and "this outcome is not in this family" are the two things §6.5's disclosure
    has to tell apart. It is empty on this workbook and the field exists anyway.

    **The reason string's shape is specified and is not free text**, because §6.5 requires the
    rendered row to carry the surviving-draw count and the floor it fell below, and a formatter
    cannot recover either from prose. Both numbers come off the `Draws` Stage 10 kept when it
    withheld the `Interval`, so the disclosure is a read and never a recomputation.

    There is no `significant` field and no threshold anywhere in this record. [§13] prescribes a
    correction, not a decision rule, and [§10] closes with "estimation, not testing, is the
    reportable output".
    """

    family: str                          # "secondary" | "safety" — Stage 9's `by_family()` key
    m_declared: int                      # |family| in the [§5] registry — 3 and 4
    m_used: int                          # p-values that entered the step-up — §6.5
    raw: dict[str, float]                # outcome -> the p Stage 10 produced, READ not recomputed
    adjusted: dict[str, float]           # outcome -> the BH-adjusted p; same keys as `raw`
    absent: tuple[tuple[str, str], ...]  # (outcome, why) for members with no p — §6.5
    descriptive_only: bool               # [§10]'s label; True for "safety" — §6.6


@dataclass(frozen=True)
class Multiplicity:
    """The two families, corrected separately, and the primary's exemption stated as a field.

    **`primary_uncorrected` is a field holding the primary's raw p and not an omission.** [§13]'s
    first sentence is "Primary outcome uncorrected", and an implementation that satisfies it by
    leaving the primary out of the record satisfies it in a way Stage 14 cannot print. The rule is
    that the primary appears beside the corrected families carrying the statement that it is
    uncorrected — so the exemption is visible where the correction is.
    """

    families: dict[str, FamilyCorrection]   # keyed as `Secondary.by_family()` is
    primary_uncorrected: float              # `boot.intervals["beta"].p` — §6.1


@dataclass(frozen=True)
class EValue:
    """The [§6] E-value for the primary estimate and for the limit nearest the null.

    **`approximation` is a FIELD carrying the string, not a literal at the print site.** The roadmap
    requires that the `RR ~ sqrt(OR)` conversion be *stated* wherever the E-value is reported, and
    `bootstrap.Interval.method`'s docstring gives the design rule this follows: which definition
    produced a number is part of what the number is. A field makes printing the number without the
    statement a thing Stage 14 has to do deliberately rather than by forgetting.

    **`e_limit` is `None` only when there is no interval at all** — an estimand below
    `C.ci_min_draws()`, unreachable on v7. It is NOT `None` when the interval spans the null: that
    case is `1.0`, exactly, and §7.3 is why the difference matters.

    **`n_draws` is the field that makes the `None` readable.** A record whose only statement is
    `e_limit is None` cannot tell "the primary lost too many replicates" from "the caller passed no
    interval", so the draw count is populated on BOTH branches.

    `spans_null` is a field and not a method for the same reason `Interval.method` is a field: three
    different routes decide it and §7.2 asserts they agree, so which answer this record carries must
    be inspectable rather than recomputed by a caller that might take a fourth route.
    """

    measure: str            # "common odds ratio [§8]" — what OR is, said once
    odds_ratio: float       # est.odds_ratio, READ off the point estimate — §7.1
    risk_ratio: float       # sqrt(odds_ratio) — the approximation, applied
    e_point: float          # the point estimate's E-value
    spans_null: bool        # decided on intervals["beta"] against 0 — §7.2
    limit: float | None     # the limit nearest the null, on the BETA scale; None when spanning
    e_limit: float | None   # 1.0 exactly when spans_null; None only when no interval exists
    n_draws: int            # the draws behind the limit — populated on BOTH branches
    approximation: str      # _APPROXIMATION — §7.5, and §19c is why it carries the citation


@dataclass(frozen=True)
class SubgroupEstimate:
    """One [§13] subgroup: the two level-specific common odds ratios and the interaction.

    **`hypothesis_generating` is `True` and it is a field rather than a docstring.** [§13] labels
    these estimates and Stage 14 has to print the label; a boolean that is structurally always True
    looks redundant until one asks where the label would otherwise live, which is a sentence in a
    reporting module that nothing asserts.

    **`n_interaction_tests` is carried on every subgroup and is the count over the WHOLE stage**,
    not this subgroup's own. §8.8 declines to put these p-values into [§13]'s Benjamini-Hochberg and
    owes the reader the denominator instead; a per-subgroup field would print `1` beside each of two
    tests. **It is `len(C.SUBGROUPS)` and never the literal 2**: a third subgroup restored to the
    registry gains its keys and its interval automatically and must not go on printing that two
    interaction tests were performed.

    `or_level` and `or_ratio` are `exp` of coefficients from ONE fit (§8.2), so they cannot disagree
    about the model. **There is no field holding the primary's odds ratio for comparison**, and §8.7
    is why: a weighted proportional-odds model is not collapsible, so the two level odds ratios need
    not bracket the primary's, a reader will expect them to, and a field inviting the comparison is
    the wrong place to answer that. No code anywhere asserts the bracketing.
    """

    subgroup: str                       # a C.SUBGROUPS key
    n: int                              # records in the fit — the [§11] denominator
    n_by_level_arm: dict[str, int]      # "S=0,A=1" -> count; the four cells §8.3 turns on
    or_level: dict[int, float]          # level -> exp(beta) and exp(beta + gamma) — §8.2
    or_ratio: float                     # exp(gamma) — the interaction, on the OR scale
    gamma: float                        # the interaction coefficient; null 0 — §8.5
    cumulative: dict[int, dict[int, dict[int, float]]]   # level -> threshold -> arm -> P(Y<=k)
    hypothesis_generating: bool         # True — [§13]
    n_interaction_tests: int            # len(C.SUBGROUPS) — the multiplicity [§13] does not address


@dataclass(frozen=True)
class Subgroups:
    """The [§13] subgroups, their draws and their intervals. Keyed as C.SUBGROUPS is.

    `draws` and `intervals` are Stage 10's own dataclasses, unmodified and not re-declared: an
    interval produced here is the same kind of object as an interval produced there, carrying the
    same `method`, the same `n_draws` and the same floor rule. A second `Interval` type would be a
    second percentile definition one edit away.

    **`intervals` may hold fewer keys than `draws`**, which is Stage 10 §3.3's rule inherited rather
    than restated: an estimand below `C.ci_min_draws()` keeps its `Draws` and gets no `Interval`.

    **There is NO `diagnostics` field here and `Arm` has one, which is a difference in the bodies
    and not an oversight.** `bootstrap.Diagnostics` tallies `n_alpha` and `polr_iterations` as
    scalars per replicate because [§10]'s body and the arm's each run exactly ONE `polr` fit. The
    subgroup body runs TWO (§8.5), so a single scalar cannot say which fit it came from and a tally
    mixing them is a distribution of nothing. Everything `subgroup_replicates` prints therefore comes
    off `Draws` — `n_attempted`, `len(draws)` and the bucket counts — and `_subgroup_replicate` sets
    `n_alpha` and `polr_iterations` to `None` deliberately.
    """

    estimates: dict[str, SubgroupEstimate]
    seed: int
    n_boot: int
    draws: dict[str, bootstrap.Draws]         # six keys — §3.3; the counters live HERE
    intervals: dict[str, bootstrap.Interval]  # the same keys; p on `.gamma` only — §8.5


@dataclass(frozen=True)
class Arm:
    """The [§13, DECISION 4] full-covariate propensity sensitivity analysis.

    **It carries its own `Propensity`, its own `Balance` and its own `Primary`, and it carries the
    PRIMARY's alongside them.** Every number in this record exists to be read beside another number,
    and the roadmap's Accept-when is a statement about a pair: the arm's `in_model`, ESS and worst
    residual |SMD| "are reported beside the primary's". A record holding only the arm's half makes
    that pairing Stage 14's job, done from two objects it has to be trusted to match up.

    **`__post_init__` raises unless `ps.spec is C.PROPENSITY_FULL` and `reference.spec is
    C.PROPENSITY_PRIMARY`**, by identity and not by equality (§4.5). This is the last of the four
    places the arm's numbers can be labelled as the primary's, and it is the one a caller reaches by
    building the record by hand.

    `differs_only_in_specification` is COMPUTED from the two `in_model` masks and never asserted. On
    this workbook it is True; on a workbook where any of the four vascular risk factors is missing on
    any row it would be False, and the clause a reader sees must be the one the data supports.

    **`diagnostics` is `bootstrap.Diagnostics`, unmodified and not re-declared**, and it is here
    because the arm's body runs exactly ONE `polr` fit per replicate — so `n_alpha` and
    `polr_iterations` mean what they mean in [§10]'s own run. `Subgroups` has no such field and its
    docstring is why. `max_abs_beta` and `or_corrected` come back empty: the arm produces no binary
    estimate and no `m_a(X)`, which is `Diagnostics`' own "AUGMENTED ONLY" rule satisfied by there
    being nothing augmented.
    """

    ps: propensity.Propensity           # spec is C.PROPENSITY_FULL — asserted, §4.5
    balance: balance.Balance            # the same 19 rows, four re-roled — §4.4
    estimate: outcome.Primary           # `outcome.primary`, UNCHANGED — §5.1
    reference: propensity.Propensity    # the [§7] fit, for the pair the Accept-when is about
    seed: int
    n_boot: int
    draws: dict[str, bootstrap.Draws]         # seven keys — §3.3
    intervals: dict[str, bootstrap.Interval]  # the same keys; `p is None` on all seven — §5.5
    diagnostics: bootstrap.Diagnostics        # ONE polr fit per replicate, so the tallies mean it
    differs_only_in_specification: bool       # COMPUTED from the two masks — §5.2

    def __post_init__(self) -> None:
        wrong = []
        if self.ps.spec is not C.PROPENSITY_FULL:
            wrong.append(
                f"`ps` carries {self.ps.spec.label!r} and this record is the "
                f"{C.PROPENSITY_FULL.label!r} arm")
        if self.reference.spec is not C.PROPENSITY_PRIMARY:
            wrong.append(
                f"`reference` carries {self.reference.spec.label!r} and it is the half of the pair "
                f"that must be {C.PROPENSITY_PRIMARY.label!r}")
        if wrong:
            raise C.SchemaError(
                "H10  " + "; ".join(wrong) + ". Every number in this record exists to be read "
                "BESIDE another number, and the roadmap's Accept-when is a statement about a pair "
                "[Stage 11 §3.1, §5.2]. A record built by hand with the halves swapped, or with the "
                "arm's own `Propensity` in both slots, reports the arm compared with itself and "
                "every cell is finite. Compared by IDENTITY and not by equality: the registry holds "
                "exactly two instances, so a copied record cannot satisfy this.")


# --- Benjamini-Hochberg [§13] ------------------------------------------------------------------------

def benjamini_hochberg(p: np.ndarray) -> np.ndarray:
    """[§13]'s Benjamini-Hochberg adjusted p-values, over ONE family. Returns them in INPUT order.

    Public, and it is a DEFINITION rather than a step (§3.2): it names no outcome, no family, no
    covariate and no centre, and it is checked against an oracle over arrays with no frame
    constructed. `statsmodels.stats.multitest.multipletests(..., method="fdr_bh")` is that oracle and
    **may not be imported here**: `pyproject.toml` declares statsmodels *"retained for unpenalised
    cross-checks in tests, not for any reported estimate"*, an adjusted p-value IS a reported
    estimate, and this is a four-line closed form (§6.4).

    Four properties, each of which is a test:

      * `kind="stable"` on the argsort, so tied p-values keep [§5] registry order. Ties are ORDINARY
        here and not theoretical: `bootstrap_p` at B = 1998 can return 1000 distinct values, so a
        family of four is drawn from a small grid.
      * The DESCENDING RUNNING MIN is the monotonicity enforcement and there is no other. `m*p/i` is
        not monotone in `i` — on sorted (0.01, 0.03, 0.04) the unenforced values are
        (0.03, 0.045, 0.04), so without it the middle outcome reports a LARGER adjusted p than the
        one below it and two of three secondary outcomes swap rank.
      * The UNSORT is `out[order] = q` and never `q`. An implementation that returns `q` returns the
        family sorted, which on a table keyed by outcome silently reattaches every p to the wrong
        outcome.
      * The cap at 1.0 is INERT and is written anyway. q_(m) = min(1, m*p_(m)/m) = p_(m) for any
        valid p, and every other q_(i) is at most that, so the cap can bind only on an input §6.1
        already rejects. It stays for the reason `ci_min_draws`' snap and `_shared_design`'s check
        stay: a specified-and-unreachable branch is one a future cohort does not discover the hard
        way.

    It is never handed an empty array: §6.5 defines the `m_used == 0` case as empty dicts plus a full
    `absent` and does not call this function at all, so `m = 0` and the division by an empty range
    are unreachable rather than defended against.
    """
    p = np.asarray(p, dtype=float)
    m = p.size
    order = np.argsort(p, kind="stable")
    q = np.minimum(1.0, p[order] * m / np.arange(1, m + 1))
    q = np.minimum.accumulate(q[::-1])[::-1]
    out = np.empty_like(q)
    out[order] = q
    return out


# --- the E-value [§6] --------------------------------------------------------------------------------

def e_value(risk_ratio: float) -> float:
    """The [§6] E-value for a risk ratio. Public, and it is a definition (§3.2).

    VanderWeele and Ding's bounding factor, on the risk-ratio scale, away from the null. Lifted from
    `pilots/analysis.py:745-748`, which is the only thing in that file lifted as code: it is the
    bound and the away-from-null inversion, correct as written.

    ROOT FIRST, THEN INVERT, and the order is pinned rather than left to the implementer:
    `sqrt(1/OR)` and `1/sqrt(OR)` are mathematically equal and NOT bit-identical for every input. On
    the odds ratios this pipeline can produce the two agree exactly, so the pin is currently inert —
    which is why it is a pin and not a claim of exactness.

    **Takes a RISK RATIO and not an odds ratio**, so the `sqrt(OR)` conversion happens at ONE call
    site (§7.5) rather than inside a function whose name does not mention it.

    **A risk ratio of exactly 1 gives exactly 1.0, with no branch, and a branch must not be
    written**: `rr*(rr-1) == 0.0` and `1 + sqrt(0) == 1.0`, all exact. `pilots/analysis.py:761` has
    an `if rd == 0: return 1.0` branch — the right answer reached by the wrong route, which is worse
    than either, because it makes the reader believe the formula needs help there.

    **Overflow is impossible, and the `sqrt(OR)` conversion is exactly what makes it so.** The naive
    form overflows when `rr**2 > DBL_MAX`; here `rr = sqrt(OR)`, so `rr**2 = OR <= DBL_MAX` always.
    No rewrite is needed and none is written; the reachability argument goes here instead of into a
    defensive form. `outcome._assert_reportable` bounds `|beta|` below `C.POLR_MAX_ABS_BETA`, so on
    this pipeline the E-value is bounded far below that again — which is worth stating because a
    reader will otherwise assume it is unbounded.

    **E is NOT smooth at the null.** `E - 1 ~ sqrt(sqrt(OR) - 1)`, so `dE/dOR -> inf` as `OR -> 1`:
    a change of 1e-12 in the odds ratio moves E by about 7e-07. Consequence for reporting — two
    decimals, and the statement that an E-value of 1.05 and one of 1.15 differ by about a hundredfold
    in the odds ratio behind them (§7.4).
    """
    rr = float(risk_ratio)
    if rr < 1.0:
        rr = 1.0 / rr
    return rr + math.sqrt(rr * (rr - 1.0))


# --- the draws-to-intervals loop, written ONCE and driven twice [§3.2] ------------------------------

def _tested_subgroup(key: str) -> bool:
    """Does [§13] prescribe a p-value for this subgroup estimand? `.gamma` and nothing else (§8.5).

    A FUNCTION and not a frozen set, which is `bootstrap._tested`'s own move for its stated reason: a
    third subgroup restored to `C.SUBGROUPS` gains its p in the same edit that gains its estimate.

    The two level odds ratios get intervals and `p is None`, because [§13] asks for **an** interaction
    test and two level-specific p-values would be two tests of one contrast — which is [§8]'s own
    objection to six threshold-wise tests, one stage on.
    """
    return key.endswith(f".{_SUBGROUP_QUANTITIES[1]}")


def _intervals(draws: dict[str, bootstrap.Draws],
               tested: Callable[[str], bool]) -> dict[str, bootstrap.Interval]:
    """Stage 10 §8.3's rule, applied twice from ONE place. `tested` is the only thing that differs.

    Both of this stage's runs turn a `dict[str, Draws]` into a `dict[str, Interval]` under the same
    three rules — the `C.ci_min_draws()` floor withholds the `Interval` and keeps the `Draws`,
    `bootstrap.percentile_ci` at `C.CI_LEVEL`, and `bootstrap.bootstrap_p` where and only where the
    stage prescribes a test — and they differ in exactly the last of them.

    A second copy of this loop is a second place the floor could be spelled, a second place
    `C.CI_LEVEL` could be defaulted, and the [§13] arm and the [§13] subgroups would then be one edit
    away from disagreeing about what a percentile interval is.

    **`bootstrap.run` keeps its own copy and is NOT refactored onto this one**, and that is a decision
    rather than an oversight: Stage 10's Definition of done asserts its numbers and this stage refits
    nothing Stage 10 refit. The divergence risk is real, is filed in Stage 11 §16 item 7, and its
    trigger is a third caller or any change to `C.ci_min_draws` or `C.PERCENTILE_METHOD` — at which
    point the loop moves into `bootstrap.py` as a ninth public name and both callers read it.

    **`bootstrap._tested` is not reused for either caller**, and §5.5 is why: `_tested("beta")` is
    `True`, so an implementation that reached for it would emit a p on the arm's `beta` that no
    prespecified section asks for.
    """
    out: dict[str, bootstrap.Interval] = {}
    for key, d in draws.items():
        if len(d.draws) < C.ci_min_draws():          # Stage 10 §8.3 — Draws kept, no Interval
            continue
        lo, hi = bootstrap.percentile_ci(d.draws, C.CI_LEVEL)
        p = bootstrap.bootstrap_p(d.draws) if tested(key) else None
        out[key] = bootstrap.Interval(lo, hi, C.CI_LEVEL, C.PERCENTILE_METHOD, len(d.draws), p)
    return out


# ======================================================================================================
# §6 — Benjamini-Hochberg within the two families [§13]
# ======================================================================================================

def _absent_reason(key: str, boot: bootstrap.Bootstrap) -> str:
    """§6.5's disclosure string, whose SHAPE is specified and is not free text.

    A formatter cannot recover the surviving-draw count or the floor from prose, and §6.5 requires
    the rendered row to carry both. Both numbers come off the `Draws` Stage 10 kept when it withheld
    the `Interval`, so this is a read and never a recomputation.
    """
    d = boot.draws[key]
    return (f"no interval: {len(d.draws)} surviving draw(s) of {d.n_attempted} against a floor of "
            f"{C.ci_min_draws()}")


def _assert_multiplicity_inputs(families: dict[str, tuple[outcome.BinaryEstimate, ...]],
                                boot: bootstrap.Bootstrap) -> None:
    """H3, H4 and H5, collected. Every one is a Stage 10 contract break and none is repairable.

    **H3 and §6.5's disclosure are the one distinction these preconditions exist to draw.** An
    absent `<outcome>.rd` INTERVAL is a sparse family and is disclosed; an absent `<outcome>.rd`
    DRAWS key, or an interval whose `p` is `None`, is a Stage 10 contract break — `bootstrap._tested`
    returns `True` for every `.rd`, so there is no route by which Stage 10 emits one untested.

    **H4's denominator is THAT KEY'S `Interval.n_draws` and never `C.N_BOOT`.** `bootstrap_p` floors
    at `1/(len(draws) + 1)` over the draws it was given, and those are that outcome's surviving ones.
    The distinction is not pedantic: at 1982 surviving draws the floor is 1/1983 while
    `1/(C.N_BOOT + 1)` is 1/2001, so a check written against `C.N_BOOT` is **strictly weaker than
    the floor it is checking** — it accepts a p below the smallest value `bootstrap_p` can return for
    that key, which is the one thing this precondition exists to catch.

    **H5 is why the primary is uncorrected rather than unreported.** [§13]'s first sentence exempts
    the primary from the correction and requires it printed beside the corrected families, so a
    missing primary p is a failure of this deliverable and not an absence to disclose. Without this
    check an absent key is a bare `KeyError` from a reported field's initialiser and a `None` is an
    unlabelled `None` in the one number [§13]'s first sentence prescribes.
    """
    bad: list[str] = []
    for family, members in families.items():
        for estimate in members:
            key = f"{estimate.outcome}.rd"
            if key not in boot.draws:
                bad.append(
                    f"H3  {key}: not a key of the Bootstrap at all, so the {family} family cannot "
                    "be corrected over the tests [§13] prescribes. The partition comes from "
                    "`Secondary.by_family()` and the p-values from `Interval.p`, and a member with "
                    "neither draws nor an interval is a [§5] registry the bootstrap did not run "
                    "[Stage 11 §6.1].")
                continue
            interval = boot.intervals.get(key)
            if interval is None:
                continue                     # §6.5's disclosure, NOT an error
            if interval.p is None:
                bad.append(
                    f"H3  {key}: the Interval exists and its `p` is None. `bootstrap._tested` "
                    "returns True for every `<outcome>.rd`, so this is a Stage 10 contract break "
                    "and not a sparse family — the sparse case keeps its Draws and gets NO Interval "
                    "[Stage 11 §6.1, §6.5].")
                continue
            floor = 1.0 / (interval.n_draws + 1)
            if not (floor <= interval.p <= 1.0):
                bad.append(
                    f"H4  {key}: a raw p of {interval.p!r} is outside "
                    f"[{floor!r}, 1.0]. The floor is `1/(n_draws + 1)` over THIS KEY'S "
                    f"{interval.n_draws} surviving draw(s) and NOT over C.N_BOOT — a check written "
                    "against N_BOOT is strictly WEAKER than the floor it checks, because it accepts "
                    "a p below the smallest value `bootstrap_p` can return for this key. A p of "
                    "exactly 0, or above 1, means something upstream produced it under a different "
                    "definition [Stage 11 §6.1].")
    beta = boot.intervals.get("beta")
    if beta is None or beta.p is None:
        bad.append(
            "H5  intervals['beta'] is " + ("absent" if beta is None else "present with p None")
            + ". [§13]'s first sentence is \"Primary outcome uncorrected\", which means UNCORRECTED "
            "and not UNREPORTED: the primary appears beside the corrected families carrying the "
            "statement that it is uncorrected, so the exemption is visible where the correction is. "
            "Without this check an absent key is a bare KeyError from a reported field's initialiser "
            "[Stage 11 §6.1].")
    if bad:
        raise C.SchemaError(
            "\n  " + "\n  ".join(bad)
            + f"\n  {len(bad)} multiplicity precondition(s) failed over "
            f"{sum(len(m) for m in families.values())} family member(s). Every one is a Stage 10 "
            "contract break and none is repairable here: this stage READS `Interval.p` and "
            "recomputes nothing, because `Interval` carries the definition that produced it and a "
            "second computation under a different one is Stage 10 §8.2's disagreement re-introduced "
            "downstream [Stage 10 §14].")


def _multiplicity_table(families: dict[str, FamilyCorrection],
                        primary: float) -> tuple[tuple[str, ...], ...]:
    """One row per family member, in [§5] registry order, plus the primary's exempt row.

    **An absent outcome is printed AS A ROW and never as a blank**, with its surviving-draw count and
    the floor it fell below in the raw-p and adjusted-p cells (§6.5). A row that vanishes is a test a
    reader does not know was declared.
    """
    header = ("outcome", "family", "m declared", "m used", "raw p", "adjusted p", "label")
    rows: list[tuple[str, ...]] = [
        (C.PRIMARY_OUTCOME, "primary", "1", "1", _fmt(primary), "uncorrected [§13]",
         "the [§13] exemption, stated where the correction is")]
    for family, correction in families.items():
        absent = dict(correction.absent)
        for key in correction.raw:
            rows.append((key, family, str(correction.m_declared), str(correction.m_used),
                         _fmt(correction.raw[key]), _fmt(correction.adjusted[key]),
                         "descriptive only [§10]" if correction.descriptive_only else "—"))
        for key, why in absent.items():
            rows.append((key, family, str(correction.m_declared), str(correction.m_used),
                         why, why, "no test was performed, so it entered no denominator"))
    return (header, *rows)


def _multiplicity_detail(families: dict[str, FamilyCorrection], primary: float) -> str:
    used = ", ".join(f"{name} m_used={c.m_used} of m_declared={c.m_declared}"
                     for name, c in families.items())
    return (
        f"[§13] Benjamini-Hochberg WITHIN each family and never over their union: {used}. "
        "`pilots/analysis.py:789-792` makes one `multipletests` call over "
        "`family.isin(['secondary','safety'])`, which is BH at m = 7 rather than [§13]'s within-family "
        "correction — a substantive departure from the plan, and the exact error these counters exist "
        "to make visible [Stage 11 §6.1, §19]. The partition is `Secondary.by_family()`'s and is "
        "never recomputed from the [§5] registry, because a consumer that groups by reading the "
        "registry itself is a second place the partition is computed and the correction is wrong if "
        "the two disagree. Every raw p is READ off the [§10] `Interval` that produced it and none is "
        "recomputed: `Interval` carries its percentile definition, and a second computation under a "
        f"different one is Stage 10 §8.2's disagreement re-introduced downstream. The primary's p, "
        f"{_fmt(primary)}, is reported UNCORRECTED and is in no family's denominator — [§13]'s first "
        "sentence, and it is printed here rather than omitted, because an implementation that "
        "satisfies the exemption by leaving the primary out satisfies it in a way no report can "
        "print. m is the number of tests PERFORMED and not the declared family size, so a member "
        "whose interval fell below the draw floor is disclosed as a row rather than paid for by "
        "every other adjusted p in its family [Stage 11 §6.5]. The safety family's adjusted p-values "
        "exist because [§13] prescribes them and carry [§10]'s descriptive-only label because [§10] "
        "deprecates them: the two sections disagree, neither mentions the other, and both are "
        "reported rather than one being chosen [Stage 11 §6.6]. There is no threshold anywhere in "
        "this entry: [§13] prescribes a correction, not a decision rule.")


def multiplicity(sec: outcome.Secondary, boot: bootstrap.Bootstrap, audit: Audit) -> Multiplicity:
    """[§13]'s Benjamini-Hochberg within the secondary and safety families. Primary uncorrected.

    Takes no family list and no method name: the partition is `sec.by_family()` and the procedure is
    [§13]'s, so a keyword for either would make a prespecified choice look like an option — Stage 6
    §6.4's rule, and Stage 11 §4.1 is what breaking it costs.

    Reads `boot.intervals[f"{k}.rd"].p` and recomputes nothing (§6.1). Raises `C.SchemaError` H3 on
    an absent draws key or a `None` p, H4 on a p outside `[1/(n_draws + 1), 1]` where `n_draws` is
    THAT KEY'S, and H5 on an absent or untested `intervals["beta"]` — and on nothing else: an absent
    `<outcome>.rd` INTERVAL is §6.5's disclosure and not an error.

    **The family size is never obtained by counting `intervals`.** Stage 10 emits 26 keys of which 8
    carry a p, so `len(boot.intervals)` is 26 and the count of p-carrying keys is 8, and neither is
    3 or 4 (§6.1).

    Appends one `model` entry (§10).
    """
    families = sec.by_family()
    _assert_multiplicity_inputs(families, boot)                      # H3, H4, H5

    corrections: dict[str, FamilyCorrection] = {}
    for family, members in families.items():
        raw: dict[str, float] = {}
        absent: list[tuple[str, str]] = []
        for estimate in members:
            key = f"{estimate.outcome}.rd"
            interval = boot.intervals.get(key)
            if interval is None:                                     # §6.5 — disclosed, not dropped
                absent.append((estimate.outcome, _absent_reason(key, boot)))
                continue
            raw[estimate.outcome] = float(interval.p)
        # `m` is the number of tests PERFORMED (§6.5). BH is not called at all when none was, so the
        # empty array and the division by an empty range are unreachable rather than defended against.
        adjusted = (dict(zip(raw, benjamini_hochberg(np.asarray(list(raw.values()), dtype=float))))
                    if raw else {})
        corrections[family] = FamilyCorrection(
            family=family, m_declared=len(members), m_used=len(raw), raw=raw,
            adjusted={k: float(v) for k, v in adjusted.items()}, absent=tuple(absent),
            descriptive_only=family == "safety")

    primary = float(boot.intervals["beta"].p)
    audit.record("model", _STEP_MULTIPLICITY, sum(c.m_used for c in corrections.values()),
                 _multiplicity_detail(corrections, primary),
                 table=_multiplicity_table(corrections, primary))
    return Multiplicity(families=corrections, primary_uncorrected=primary)


# ======================================================================================================
# §7 — the E-value for the primary estimate [§6]
# ======================================================================================================

def _e_value_table(ev: EValue, bal: balance.Balance,
                   est: outcome.Primary) -> tuple[tuple[str, ...], ...]:
    """The two E-values, the worst residual |SMD| beside them, and the 15%-band verdict per threshold.

    **The worst residual |SMD| is in this table and that is structural rather than editorial**
    (§7.6). The E-value bounds UNMEASURED confounding, and this analysis has measured confounding it
    has not fixed: [§9]'s threshold is exceeded after weighting on centre. An E-value printed beside
    a residual |SMD| of that size invites the reading "residual confounding would have to be strong"
    when a covariate that is IN THE MODEL is visibly not balanced. `Balance.worst()` returns the
    undefined covariate names with the number, so a worst |SMD| computed over rows some of which
    could not be judged cannot be quoted alone.

    **And one row per `C.MRS_THRESHOLDS` entry** (§7.5, §12.2 item 4a): the control-arm `P(Y <= k)`
    and whether it clears the 15% prevalence boundary the `sqrt(OR)` conversion is licensed for. It
    is read off `est.cumulative`, which this function already receives, so nothing is recomputed —
    and **this is the only place those numbers exist**, because Stage 11 §1 keeps them out of the
    specification and Stage 14 cannot name the thresholds outside the band from a log that does not
    print them.
    """
    worst, undefined = bal.worst()
    header = ("quantity", "value", "note")
    rows: list[tuple[str, ...]] = [
        ("measure", ev.measure, "what the odds ratio IS, said once"),
        ("odds ratio", _fmt(ev.odds_ratio), "READ off the point estimate; never from the draws"),
        ("risk ratio", _fmt(ev.risk_ratio), ev.approximation),
        ("E-value (point)", _fmt(ev.e_point), "the minimum RR with BOTH treatment and outcome an "
                                              "unmeasured confounder would need"),
        ("interval spans the null", "yes" if ev.spans_null else "no",
         "decided on intervals['beta'] against 0; ties at the null count as spanning"),
        ("limit nearest the null", "not applicable" if ev.limit is None else _fmt(ev.limit),
         "selected BY SIGN and never by the smaller absolute limit, which picks the wrong one only in the "
         "spanning case where the 1.0 rule then masks it"),
        ("E-value (limit)", "not computed" if ev.e_limit is None else _fmt(ev.e_limit),
         "EXACTLY 1.0 when the interval spans the null: no confounding at all is needed for these "
         "data to be compatible with no effect. Never nan, and never the formula applied to a "
         "null-crossing limit, which returns a large number reading as robustness"),
        ("surviving draws", str(ev.n_draws), "the denominator behind the limit"),
        ("worst residual abs SMD",
         "undefined on every row" if worst is None else _fmt(abs(worst.weighted)),
         "REPORTED HERE and not elsewhere: the E-value bounds UNMEASURED confounding and this "
         + (f"analysis has measured confounding it has not fixed — {worst.covariate} under "
            f"{bal.spec.label} [§9]" if worst is not None else "table has no judgeable row")),
        ("rows no SMD could be computed on", str(len(undefined)),
         ", ".join(undefined) if undefined else "none"),
    ]
    for k in C.MRS_THRESHOLDS:
        control = est.cumulative[k][min(C.TREATMENT_LABELS)]
        rows.append((f"control P(Y <= {k})", _fmt(control),
                     "above the 15% prevalence the sqrt(OR) conversion is licensed for"
                     if control > 0.15 else
                     "BELOW 15%: the sqrt(OR) conversion is NOT licensed at this threshold"))
    return (header, *rows)


def _e_value_detail(ev: EValue, bal: balance.Balance) -> str:
    worst, undefined = bal.worst()
    return (
        f"[§6] E-value for the {ev.measure} and for the confidence limit nearest the null, under the "
        f"{bal.spec.label} specification. {ev.approximation} The direction of the error is statable "
        "and it is the conservative one: E is strictly increasing in RR away from the null and "
        "sqrt(OR) < OR for OR > 1, so the prescribed conversion yields a SMALLER E-value than "
        "treating the odds ratio as a risk ratio. No second E-value at RR = OR is reported as a "
        "bracket — that is a second number for one estimand, and monotonicity already tells a reader "
        "which way it moves [Stage 11 §7.5]. E is NOT smooth at the null: E - 1 ~ sqrt(sqrt(OR) - 1), "
        "so an E-value of 1.05 and one of 1.15 differ by about a hundredfold in the odds ratio behind "
        "them, and two decimals is the resolution this number has [Stage 11 §7.4]. There is no "
        "interval on the E-value itself: both are deterministic functions of numbers that already "
        "carry intervals [Stage 10 §18]. "
        + ("The interval SPANS the null, so the limit's E-value is exactly 1.0 — no confounding at "
           "all is needed for these data to be compatible with no effect, and feeding a "
           "null-crossing limit into the formula instead answers a different question and returns a "
           "large number that reads as robustness while meaning the opposite [Stage 11 §7.3]. "
           if ev.spans_null else
           "The interval excludes the null, so the limit nearest it is the one selected BY SIGN and "
           "its E-value is the formula's [Stage 11 §7.3]. ")
        + "**It is never printed alone.** The E-value bounds UNMEASURED confounding and this analysis "
        "has measured confounding it has not fixed, so the worst residual |SMD| from the SAME "
        "specification is in the table above"
        + (f" — {_fmt(abs(worst.weighted))} on {worst.covariate}" if worst is not None else "")
        + (f", over {len(undefined)} row(s) no SMD could be computed on: " + ", ".join(undefined)
           if undefined else ", and no row is undefined")
        + ". [§9]'s amendment of 2026-08-24 requires the exceedance named at the point of reading the "
        "estimate, and this is such a point [Stage 11 §7.6]. DECISION 4 draws the same line from the "
        "other direction: the negative controls show confounding that is measured and visible, the "
        "E-value addresses confounding that is not, and the two are different questions. The table "
        "also carries the control-arm cumulative probability at each of the six [§8] thresholds "
        "against the conversion's documented 15% boundary, because `exp(beta)` is a cumulative-odds "
        "shift assumed constant across all six and each threshold has its own baseline risk "
        "[Stage 11 §7.5, §12.2 item 4a].")


def e_value_primary(est: outcome.Primary, bal: balance.Balance,
                    boot: bootstrap.Bootstrap, audit: Audit) -> EValue:
    """[§6]'s E-value for the primary estimate and for the limit nearest the null.

    **`bal` is a parameter and it is the whole of §7.6.** The E-value bounds UNMEASURED confounding
    and this analysis has measured confounding it has not fixed; the worst residual |SMD| is printed
    in the same table, from the same specification, and a signature that did not require it would
    make the pairing a convention a reporting layer could drop.

    Takes the `Primary` for `odds_ratio` and the `Bootstrap` for `intervals["beta"]`, and nothing
    from the draws: the point estimate is not a function of them (Stage 10 §8.1) and neither is the
    limit a function of the point estimate — §7.3's sign guard is where that stops being an academic
    remark.

    **Null-spanning is decided on `beta` against 0**, and three routes reach the same verdict for
    three different reasons: `lo <= 0 <= hi` iff `exp(lo) <= 1 <= exp(hi)`, because `exp` is strictly
    increasing and `exp(0.0) == 1.0` exactly; `exp` of a percentile limit IS the percentile limit of
    `exp` of the draws, because `C.PERCENTILE_METHOD` selects an ORDER STATISTIC and a strictly
    monotone transform commutes with order-statistic selection exactly; and `p >= 1 - C.CI_LEVEL` iff
    the interval includes 0, which is Stage 10 §9.4's asserted agreement. **Ties at the null count as
    spanning**, consistently with `bootstrap_p` counting a draw of exactly 0.0 in both tails.

    Raises `C.SchemaError` H6 on a non-finite input and H7 on an estimate outside its own interval
    (§7.3). Returns `1.0` for `e_limit` when the interval spans the null; never `nan`, ever. Fills
    `n_draws` on both branches, so an absent interval reports HOW FEW draws there were rather than
    only that there were too few.

    Appends one `model` entry (§10).
    """
    odds_ratio = float(est.odds_ratio)
    interval = boot.intervals.get("beta")
    n_draws = interval.n_draws if interval is not None else len(boot.draws["beta"].draws)

    _assert_e_value_inputs(est, interval)                            # H6, H7

    if interval is None:
        spans_null, limit, e_limit = False, None, None
    elif interval.lo > 0.0:
        spans_null, limit = False, float(interval.lo)
        e_limit = e_value(math.sqrt(math.exp(limit)))
    elif interval.hi < 0.0:
        spans_null, limit = False, float(interval.hi)
        e_limit = e_value(math.sqrt(math.exp(limit)))
    else:
        # §7.3's hard rule. The interval already includes the null, so the minimum strength of
        # unmeasured confounding needed to explain the result away is NONE, and the formula returns
        # `1 + sqrt(1*0) = 1.0` for exactly that. Not approximately, not nan, and NEVER the formula
        # applied to a null-crossing limit — which is minimised at 1 and increases in BOTH directions
        # once inverted, so a limit on the far side of the null scores as strong evidence.
        spans_null, limit, e_limit = True, None, 1.0

    ev = EValue(measure=_MEASURE, odds_ratio=odds_ratio, risk_ratio=math.sqrt(odds_ratio),
                e_point=e_value(math.sqrt(odds_ratio)), spans_null=spans_null, limit=limit,
                e_limit=e_limit, n_draws=int(n_draws), approximation=_APPROXIMATION)
    audit.record("model", _STEP_EVALUE, int(est.in_estimate.sum()),
                 _e_value_detail(ev, bal), table=_e_value_table(ev, bal, est))
    return ev


def _assert_e_value_inputs(est: outcome.Primary,
                           interval: bootstrap.Interval | None) -> None:
    """H6 and H7, collected.

    **H6 — a non-finite odds ratio or interval limit.** `C.SchemaError`, never `nan`. Unreachable by
    construction — `bootstrap.collect` raises on a non-finite draw and a percentile of finite draws
    is finite, and `outcome._assert_reportable` bounds `|beta|` — so a non-finite input means a
    caller passed something Stage 10 cannot produce. `pilots/analysis.py:740, 759` returns `np.nan`
    instead, which prints as `nan` in a manuscript table.

    **H7 — the point estimate is outside its own interval.** Stage 10 §8.1 states plainly that a
    percentile interval *"is not a function of the point estimate"*, so `beta < 0` with `lo > 0` is
    arithmetically reachable. That is a finding to look at, not a number to print.
    `pilots/analysis.py:761` tests `np.sign(rd_limit) != np.sign(rd)` and returns 1.0, which reports
    it as "the interval spans the null" — a different, benign condition, and it is the sharpest
    single reason this guard raises.
    """
    bad: list[str] = []
    if not np.isfinite(est.odds_ratio) or est.odds_ratio <= 0.0:
        bad.append(
            f"H6  the point estimate's odds ratio is {est.odds_ratio!r}, which is not a positive "
            "finite number. The E-value is a function of it and `nan` prints as `nan` in a "
            "manuscript table, so this raises rather than returning one [Stage 11 §7.4].")
    if interval is not None:
        if not (np.isfinite(interval.lo) and np.isfinite(interval.hi)):
            bad.append(
                f"H6  intervals['beta'] carries limits ({interval.lo!r}, {interval.hi!r}). "
                "`bootstrap.collect` raises on a non-finite draw and a percentile of finite draws is "
                "finite, so this is a limit Stage 10 cannot produce [Stage 11 §7.4].")
        elif not (interval.lo <= est.beta <= interval.hi):
            bad.append(
                f"H7  the point estimate beta = {est.beta!r} is outside its own interval "
                f"[{interval.lo!r}, {interval.hi!r}]. A percentile interval is an order statistic of "
                "the replicate distribution and is NOT a function of the point value [Stage 10 §8.1], "
                "so this is arithmetically reachable and it is a finding to look at rather than a "
                "number to print. It is NOT the interval spanning the null, which is a different and "
                "benign condition the pilot reports this as [Stage 11 §7.3].")
    if bad:
        raise C.SchemaError(
            "\n  " + "\n  ".join(bad)
            + f"\n  {len(bad)} E-value precondition(s) failed. Neither is repairable and neither is "
            "clipped: `pilots/analysis.py:742` repairs an out-of-range implied risk into a sentinel "
            "and returns a finite, reportable E-value from it, which is a fallback producing a "
            "number and roadmap invariant 5 forbids it.")


# ======================================================================================================
# §8 — the subgroups [§13]
# ======================================================================================================
#
# ONE weighted proportional-odds fit per subgroup, over the whole overlap-weighted population, with a
# treatment × subgroup interaction [Stage 11 §8.1(a)]. [§13] names two subgroups and asks for
# "subgroup estimates ... with an interaction test", hypothesis-generating; it names no model, no
# scale, no estimator and no reported quantity, and all four are decided in §8.1 and recorded as
# PI-reversible in §17.
#
#     logit P(Y <= k | A, S)  =  alpha_k  +  beta*A  +  delta*S  +  gamma*(A x S)
#
#       exp(beta)          the common odds ratio at S = 0
#       exp(beta + gamma)  the common odds ratio at S = 1
#       exp(gamma)         the ratio of the two — the interaction, on the odds-ratio scale
#       gamma              the tested coefficient; null gamma = 0
#
# **The propensity model is NEVER refit within a level**, and that is an estimand argument rather than
# a power argument. [§8] states that the ATO estimand *"is itself indexed by the true propensity
# score"* — the target population is defined by h(X) = e(X){1 - e(X)} and corresponds to no observable
# subset of the data — and contrasts it with the ATT and ATC, whose target populations *"do not move
# when the propensity model changes"*. A within-level `e` defines a different tilting function, so two
# level effects would be effects in two DIFFERENT weighted populations and their difference is not an
# interaction [Stage 11 §8.1(b)].
#
# **No [§6] covariate enters, and this is not an omission.** [§8]: all estimates are marginal in the
# overlap population; no covariate adjustment of the reported effect. Adding the confounders here
# would make `exp(beta)` a CONDITIONAL odds ratio — precisely what the roadmap forbids Stage 12 from
# mislabelling. Confounding control is the weights' job and it has already been done. `S` in the
# linear predictor is the ONE declared departure from [§8]'s weighting-only rule: without the main
# effect `gamma` is not an interaction, and it is minimal — one column, the subgroup variable itself.
#
# **Level coding is arithmetic and not `C.REFERENCE_LEVELS`.** Both subgroup columns are `Int64` 0/1
# with `<NA>` and level 1 is the named condition. This is written down because the pilot takes its
# reference level from `pd.get_dummies(series.astype(str), drop_first=True)`
# (`pilots/analysis.py:827`), which makes the SIGN of `gamma` a function of string sort order.
#
# **The frozen median is resampled and never recomputed.** `core_above_median` is cohort-dependent;
# `derive.derive_cohort` computes it once after the [§3] restrictions and raises if asked twice,
# because recomputing per replicate would give every replicate its own cut-point and its own
# subgroup. `bootstrap.resample` copies whole rows, so the column travels with the frame. Nothing in
# this module names `derive_cohort` or `_core_above_median`, and §15.11 asserts that by AST scan.


def _subgroup_fit(df: pd.DataFrame, ps: propensity.Propensity, in_estimate: pd.Series,
                  s: str) -> tuple[float, float, model.PolrFit, pd.Series]:
    """One subgroup's [§13] interaction fit. Returns (beta, gamma, fit, mask).

    Three columns and no [§6] covariate (§8.2). The product is built HERE, on a local copy, and
    passed to `model.design` as an ORDINARY covariate, so D1-D4 and the constant-column rule apply to
    it uniformly and `dropped` is reported. Building it after `design` would exempt it from the one
    rule that catches G8.

    Two assertions, and each is a failure if MOVED rather than merely if removed:

      * **G8 sits BETWEEN `design` and `polr`**, because a two-column design FITS — `gamma` simply
        ceases to exist while every number returned is finite. This is
        `outcome._assert_exposure_survived`'s placement and its reason, on a different column.
      * **G9 sits BETWEEN `polr` and everything that reads it**, because a separated ordinal fit
        CONVERGES and its `exp(beta)` is a finite float every downstream table will accept —
        `model.polr`'s own docstring measures it at 36.4. This is `outcome._assert_reportable`'s
        placement and its reason.

    Both raise `model.FitError` and not `C.SchemaError`: a stratified resample that drew no treated
    patient at a level is a sparse replicate, which [§10] drops and counts, and neither is ever
    substituted with a different estimator [roadmap invariant 5]. **That reasoning is about the
    replicate loop, and this function is also called once on the POINT ESTIMATE, where there is
    nothing to drop into — so the caller converts a `FitError` from the point-estimate call into
    `C.SchemaError` H8 (§8.5). The granularity decision lives here; the point-estimate decision does
    not.**

    **G8 asserts the EXACT column tuple and not "nothing was lost".** `design` ends
    `X = ....astype(float)` after a `get_dummies` over the declared factors only, and none of these
    three columns is in `C.CATEGORICAL` — so the returned columns are the passed tuple, in the passed
    order, and the equality is the strongest available check rather than a hopeful one. A weaker
    predicate over set membership would pass a future `design` that reordered columns while `beta` and
    `gamma` are read back BY NAME below, so the two would silently stop describing the same fit.

    **The subgroup mask is a no-op on v7 and is written anyway.** Both subgroup columns derive from
    [§6] covariates, so `model.complete_cases` has already excluded any row missing either. Without
    the mask an `Int64` `<NA>` converts to `nan`, reaches `polr` and raises O1 — droppable and
    counted — so the failure would be absorbed as a sparse replicate rather than reported as a fact
    about the frame.

    **The bound is on `max(|beta|, |gamma|, |beta + gamma|)` and NOT on `|delta|`**, which is G7's
    own stated logic: the bound is on the REPORTED effect and not on a nuisance that can legitimately
    be large. `delta` is the subgroup main effect and grows when a level's outcome distribution is
    concentrated, and dropping a replicate on it would lose the interaction estimate for a reason
    unrelated to the interaction (§8.4).
    """
    mask = in_estimate & df[s].notna()
    sub = df.loc[mask].copy()                      # local; the caller's frame is not touched
    ix = f"{C.TREATMENT}_x_{s}"
    sub[ix] = sub[C.TREATMENT].astype(float) * sub[s].astype(float)

    declared = (C.TREATMENT, s, ix)
    X, dropped = model.design(sub, declared)
    if tuple(X.columns) != declared or dropped:                  # G8 — the EXACT tuple
        raise model.FitError(
            f"G8  the treatment x subgroup design came back as {tuple(X.columns)} against a declared "
            f"{declared}, dropping {', '.join(dropped) or 'nothing'}. `design` removes a constant "
            "column SILENTLY, after which `polr` fits the remainder, converges, and returns a result "
            "in which the interaction does not exist — so this is checked BEFORE the fit and not "
            "after it [Stage 11 §8.3].")

    fit = model.polr(X, sub[C.PRIMARY_OUTCOME].to_numpy(dtype=float),
                     ps.w.loc[mask].to_numpy(dtype=float))     # raises O1-O6; O6 is route two
    coefficients = dict(zip(X.columns, (float(v) for v in fit.beta)))
    beta, gamma = coefficients[C.TREATMENT], coefficients[ix]

    reported = max(abs(beta), abs(gamma), abs(beta + gamma))
    if reported >= C.POLR_MAX_ABS_BETA:                          # G9
        raise model.FitError(
            f"G9  a reported quantity reached {reported:.3f} against a bound of "
            f"{C.POLR_MAX_ABS_BETA}. The fit CONVERGED; a separated proportional-odds fit does. The "
            "bound is on max(|beta|, |gamma|, |beta+gamma|) and NOT on the subgroup main effect, "
            "which is a nuisance that grows legitimately when a level's outcome distribution is "
            "concentrated [Stage 11 §8.4].")
    return beta, gamma, fit, mask


def _subgroup_replicate(draw: pd.DataFrame, keys: tuple[str, ...],
                        source: object) -> bootstrap.Replicate:
    """One replicate: refit [§7] over the WHOLE population, then one interaction fit per subgroup.

    **The propensity model is refit in every replicate and over the whole `in_model` population of
    the draw, never within a level.** [§10] requires the refit and §8.1(b) is why the population is
    not the level's. This is the single line that separates the chosen design from the rejected one.

    **`propensity.fit` and never `fit_full`.** The subgroups are estimated under the ONE prespecified
    specification (§8.2); the [§13] arm is a different deliverable and combining them would produce a
    subgroup estimate under a sensitivity specification that no section asks for.

    Failure granularity mirrors Stage 10 §7.1 and §7.3 exactly: a `propensity.fit` failure fails the
    whole replicate and is counted against all six keys; a subgroup's own fit failure costs that
    subgroup's THREE keys atomically and costs the other subgroup nothing.

    One propensity fit serves both subgroups, because it is the expensive part and a second run for
    the second subgroup would double it for nothing.

    **`n_alpha` and `polr_iterations` are `None` and that is the decision, not a shortcut.** They are
    SCALARS on `Replicate` because [§10]'s body and the arm's each run one `polr` fit; this body runs
    two, so any single value would silently describe one subgroup while being tallied as the
    replicate's. `Subgroups` therefore carries no `bootstrap.Diagnostics` and every counter it reports
    comes off `Draws`. `sum_w` and `n_in_model` ARE populated, because there is exactly one propensity
    fit and they are unambiguous.

    The `Audit` is constructed here and discarded, which is Stage 10 §6.2 unchanged: nine entries per
    replicate at `N_BOOT` would be memory that is a function of the replicate count, and no
    replicate's log is ever written.

    `C.SchemaError` is never caught. The only `except` is on `model.FitError`, which is Stage 10 §7.1
    as code and for its reason: a frame this stage constructed that cannot be read means this stage
    constructed it wrongly.
    """
    audit = Audit(source)
    try:
        ps = propensity.fit(draw, audit)                   # the [§7] specification — NOT fit_full
    except model.FitError as failure:
        return bootstrap.Replicate(                        # the whole replicate — §7.1, inherited
            values={}, failures={key: bootstrap.bucket(str(failure)) for key in keys},
            n_alpha=None, polr_iterations=None, sum_w=None, n_in_model=None)

    in_estimate = outcome.estimation_population(draw, ps)  # ONE definition of the [§11] mask
    values: dict[str, float] = {}
    failures: dict[str, str] = {}
    for s in C.SUBGROUPS:                                  # declared order; no set iteration
        group = tuple(f"{s}.{q}" for q in _SUBGROUP_QUANTITIES)        # THREE, atomically — §7.3
        try:
            beta, gamma, _, _ = _subgroup_fit(draw, ps, in_estimate, s)
        except model.FitError as failure:                  # G8, G9, O1-O6 — this subgroup only
            failures.update({key: bootstrap.bucket(str(failure)) for key in group})
            continue
        values.update({group[0]: beta, group[1]: gamma, group[2]: beta + gamma})

    return bootstrap.Replicate(
        values=values, failures=failures,
        n_alpha=None, polr_iterations=None,                # TWO fits, one scalar field — see above
        sum_w=float(ps.w[ps.in_model].sum()), n_in_model=int(ps.in_model.sum()))


def _assert_subgroup_inputs(df: pd.DataFrame, ps: propensity.Propensity,
                            est: outcome.Primary) -> None:
    """H9 — the caller error `bootstrap.run`'s R-series covers, on the one function here that takes a
    frame, a `Propensity` and a `Primary` together.

    `est.in_estimate` is the [§11] population the point-estimate fits run on and `ps.w` is what
    weights them, so if the two describe different rows the fits and the weights describe different
    patients — and `.loc[boolean_series]` returns rows in the SERIES' order, so a permuted-but-equal
    index pairs each record's weight with another record's outcome and returns a different number with
    no raise (Stage 9 §12, R3's own argument).
    """
    bad: list[str] = []
    if not est.in_estimate.index.equals(df.index) or not ps.in_model.index.equals(df.index):
        bad.append(
            "H9  the Primary's `in_estimate` or the Propensity's `in_model` carries a different index "
            "from the frame. `.index.equals` and NOT a length check or a set comparison: a "
            "permuted-but-equal index pairs each record's weight with another record's outcome and "
            "returns a different number with no raise [Stage 9 §12, Stage 11 §9.1].")
    elif not bool((est.in_estimate <= ps.in_model).all()):
        bad.append(
            f"H9  `in_estimate` holds {int(est.in_estimate.sum())} record(s) and is not a subset of "
            f"`in_model`'s {int(ps.in_model.sum())}: "
            f"{int((est.in_estimate & ~ps.in_model).sum())} record(s) are in the estimate and carry "
            "no weight. [§11]'s population is `in_model & outcome present` [Stage 8 §4.1], so a "
            "record outside `in_model` reaching a weighted fit is a `nan` weight in a sum, and "
            "Stage 6 §9 forbids filling it.")
    if bad:
        raise C.SchemaError(
            "\n  " + "\n  ".join(bad)
            + f"\n  {len(bad)} subgroup precondition(s) failed against {len(df)} records. Both are "
            "properties of the CALL and not of a replicate, so both are the caller's bug — which is "
            "why they are `C.SchemaError` and the H series, and why no `bootstrap.bucket` will ever "
            "see them [Stage 11 §9.1].")


def _cells(df: pd.DataFrame, mask: pd.Series, s: str) -> dict[str, int]:
    """The four (level, arm) cell counts, keyed `"S=<level>,A=<code>"`, over the fitted records.

    **This is the one number that governs §8.3 and §8.4 both**, and it is rendered per subgroup so a
    reader of the log can see it without the specification: an interaction estimated over a cell of
    five treated patients is separated or near-separated in about one replicate in eight, and the
    bound that drops those replicates is selection on the estimate (§8.4).

    Ranges over the DECLARED levels and the DECLARED arms, so an empty cell is a row reading 0 rather
    than a row nobody sees — Stage 6 §7.5's rule.
    """
    sub = df.loc[mask]
    return {f"S={level},A={code}":
            int(((sub[s] == level) & (sub[C.TREATMENT] == code)).sum())
            for level in (0, 1) for code in C.TREATMENT_LABELS}


def _subgroup_estimate(df: pd.DataFrame, ps: propensity.Propensity, s: str,
                       beta: float, gamma: float, mask: pd.Series) -> SubgroupEstimate:
    """The record for one subgroup, from ONE fit's two coefficients (§8.2) plus §8.6's table.

    **The six cumulative `RD_k` are NOT reported per level, with intervals or without** (§8.6). [§13]
    does not ask for them; twelve more intervals from a hypothesis-generating analysis is twelve more
    numbers a reader will read as findings; and [§8]'s own argument against threshold-wise reporting
    — that it invites selection of the most favourable cut — applies with more force per level.

    What is reported instead is the pair of weighted cumulative distributions per level, POINT VALUES
    ONLY, no interval and no p, from the landed `outcome.cumulative_rd`. It adds no estimator, and it
    is the one diagnostic that adjudicates §8.1(a)'s extra assumption over (a′)'s: it is [§16]'s
    constant-shift statement evaluated per level, so a reader can see whether a single
    proportional-odds shift for `delta` does violence to the data.
    """
    sub = df.loc[mask]
    cumulative: dict[int, dict[int, dict[int, float]]] = {}
    for level in (0, 1):
        at = sub[s] == level
        _, distribution = outcome.cumulative_rd(
            sub.loc[at, C.PRIMARY_OUTCOME].to_numpy(dtype=float),
            sub.loc[at, C.TREATMENT].to_numpy(dtype=float),
            ps.w.loc[mask][at].to_numpy(dtype=float))
        cumulative[level] = distribution
    return SubgroupEstimate(
        subgroup=s, n=int(mask.sum()), n_by_level_arm=_cells(df, mask, s),
        or_level={0: float(np.exp(beta)), 1: float(np.exp(beta + gamma))},
        or_ratio=float(np.exp(gamma)), gamma=float(gamma), cumulative=cumulative,
        hypothesis_generating=True, n_interaction_tests=len(C.SUBGROUPS))


def _subgroup_fit_table(estimate: SubgroupEstimate) -> tuple[tuple[str, ...], ...]:
    """The four cells, the two level odds ratios, the interaction, and §8.6's cumulative table.

    **No row holds the primary's odds ratio for comparison** (§8.7). A weighted proportional-odds
    model is not collapsible, so the primary `exp(beta)` is NOT a weighted average of these two and
    they need not bracket it; a reviewer's instinct is to add
    `assert min(or_level) <= primary_or <= max(or_level)` as a sanity check, and it would fail on
    correct output. Nothing in this module or its tests writes that comparison.
    """
    header = ("quantity", "value", "note")
    rows: list[tuple[str, ...]] = [
        ("n fitted", str(estimate.n), "the [§11] denominator this subgroup's fit ran on"),
        ("interaction tests performed", str(estimate.n_interaction_tests),
         "over the WHOLE stage, and they are an UNCORRECTED multiplicity [§13] does not address"),
        ("label", "hypothesis-generating [§13]",
         "the label is a FIELD on the record, not a sentence in a reporting module"),
    ]
    for cell, count in estimate.n_by_level_arm.items():
        rows.append((f"n {cell}", str(count),
                     "an interaction rests on the SMALLEST of these four, and a level with few "
                     "treated patients is separated or near-separated in a large fraction of "
                     "replicates [Stage 11 §8.4]"))
    rows.extend([
        ("odds ratio at S = 0", _fmt(estimate.or_level[0]), "exp(beta), from ONE fit"),
        ("odds ratio at S = 1", _fmt(estimate.or_level[1]), "exp(beta + gamma), from the SAME fit"),
        ("interaction odds ratio", _fmt(estimate.or_ratio),
         "exp(gamma) — the ratio of the two; NOT a weighted average of them, because a weighted "
         "proportional-odds model is not collapsible and the two levels need not bracket the "
         "primary's odds ratio [Stage 11 §8.7]"),
        ("gamma", _fmt(estimate.gamma), "the tested coefficient; the null is gamma = 0"),
    ])
    for level, distribution in estimate.cumulative.items():
        for threshold, by_arm in distribution.items():
            rows.append((
                f"S = {level}: P(Y <= {threshold})",
                " / ".join(f"{label} {_fmt(by_arm[code])}"
                           for code, label in C.TREATMENT_LABELS.items()),
                "point values only, no interval and no p — [§16]'s constant-shift statement "
                "evaluated PER LEVEL [Stage 11 §8.6]"))
    return (header, *rows)


def _subgroup_fit_detail(estimate: SubgroupEstimate, s: str) -> str:
    treated_at_zero = estimate.n_by_level_arm[f"S=0,A={max(C.TREATMENT_LABELS)}"]
    treated_at_one = estimate.n_by_level_arm[f"S=1,A={max(C.TREATMENT_LABELS)}"]
    return (
        f"[§13] subgroup estimate on the [§5] primary outcome: {C.SUBGROUPS[s]}. ONE weighted "
        "proportional-odds fit over the whole overlap-weighted population, with a treatment x "
        f"subgroup interaction, on {estimate.n} record(s). Three columns — treatment, {s}, and their "
        "product — and NO [§6] covariate: [§8] makes every estimate marginal in the overlap "
        "population, and adding the confounders here would make exp(beta) a CONDITIONAL odds ratio. "
        "Confounding control is the weights' job and it has already been done. The subgroup variable "
        "itself is the one declared departure from [§8]'s weighting-only rule, and it has to be "
        "there: without the main effect, gamma is not an interaction [Stage 11 §8.2]. The propensity "
        "model is the ONE prespecified specification's and is never refit within a level — the ATO "
        "estimand is indexed by the propensity score, so a within-level score defines a different "
        "target population and the difference of two such effects is not an interaction [§8, Stage "
        f"11 §8.1]. HYPOTHESIS-GENERATING [§13], and {estimate.n_interaction_tests} interaction "
        "test(s) were performed across this stage — themselves an uncorrected multiplicity that "
        "[§13] does not address, which is where [§13] is silent and this analysis must not be. The "
        "interaction p is NOT in [§13]'s Benjamini-Hochberg: that clause is scoped to the secondary "
        "and safety OUTCOME families, these estimates are on the primary outcome whose first "
        "prescription is \"uncorrected\", and an FDR correction controls a discovery rate among "
        "hypotheses one intends to act on [Stage 11 §8.8]. "
        f"The level cells are {treated_at_zero} treated at S = 0 and {treated_at_one} at S = 1, and "
        "they are printed because an interaction rests on the smallest of the four: a level with few "
        "treated patients is separated or near-separated in a large fraction of replicates, and the "
        "bound that drops those replicates is a condition on the ESTIMATE rather than on the frame "
        "[Stage 11 §8.4]. The two level odds ratios come from one fit and cannot disagree about the "
        "model; they need NOT bracket the primary's, because a weighted proportional-odds model is "
        "not collapsible [Stage 11 §8.7]. The six cumulative RD_k are not reported per level "
        "[Stage 11 §8.6]; the two weighted cumulative distributions above are, as point values with "
        "no interval and no p.")


def _subgroup_replicates_table(draws: dict[str, bootstrap.Draws],
                               n_in_model: dict[int, int]) -> tuple[tuple[str, ...], ...]:
    """Per key: `n_attempted`, the surviving count and the bucket counts — all off `Draws`.

    **There is no `Diagnostics` behind this table and that is `Subgroups`' own decision.** A two-fit
    body cannot fill `Replicate`'s scalar `n_alpha`, so a tally of it would be a distribution of
    nothing; every counter here comes off `Draws`, which is where the drop rate lives. The
    `n_in_model` range is tallied from the ONE propensity fit each replicate made, which is the one
    consumer those two fields have.
    """
    buckets = sorted(set(C.FAILURE_BUCKETS.values()))
    header = ("estimand", "attempted", "surviving", "drop rate", "tested", *buckets)
    rows: list[tuple[str, ...]] = []
    for key, d in draws.items():
        dropped = d.n_attempted - len(d.draws)
        rows.append((
            key, str(d.n_attempted), str(len(d.draws)),
            _fmt(dropped / d.n_attempted) if d.n_attempted else "missing",
            "yes" if _tested_subgroup(key) else "no",
            *(str(d.failures.get(bucket, 0)) for bucket in buckets)))
    rows.append((
        "in_model per replicate",
        f"{min(n_in_model)} .. {max(n_in_model)}" if n_in_model else "missing",
        str(sum(n_in_model.values())), "—", "—", *("—" for _ in buckets)))
    return (header, *rows)


def _subgroup_replicates_detail(draws: dict[str, bootstrap.Draws],
                                intervals: dict[str, bootstrap.Interval]) -> str:
    worst = min((len(d.draws) for d in draws.values()), default=0)
    attempted = max((d.n_attempted for d in draws.values()), default=0)
    return (
        f"[§13] subgroup bootstrap: {C.N_BOOT} stratified replicates at seed {C.SEED}, stratified "
        f"within {C.BOOT_STRATUM}, driving {len(draws)} estimand(s) — three per subgroup, which fail "
        "as ONE group because they come from one fit [Stage 10 §7.3]. The propensity model is refit "
        "in every replicate and over the whole in_model population of the draw, never within a level "
        "[Stage 11 §8.1]. ONE fit serves both subgroups, so a propensity failure costs all "
        f"{len(draws)} keys and a subgroup's own failure costs that subgroup's three and the other "
        "subgroup nothing [Stage 11 §8.5]. "
        f"The worst-off estimand keeps {worst} of {attempted} replicate(s), against a floor of "
        f"{C.ci_min_draws()}: {len(intervals)} of {len(draws)} estimand(s) carry an interval, and one "
        "below the floor keeps its Draws and gets NO Interval [Stage 10 §8.3]. "
        "**A drop through G9 is SELECTION ON THE ESTIMATE and the surviving interval is therefore "
        "narrower and not merely noisier**, which is not the same statement as a count of lost "
        "replicates: G9's condition is on max(|beta|, |gamma|, |beta+gamma|), a condition on the "
        "quantity being estimated rather than on the frame, so the survivors are a sample selected "
        "on the estimate and the percentile limits are quantiles of a TRUNCATED sampling "
        "distribution. Truncation removes mass from both tails and none from the middle, so the "
        "limits are systematically tighter than the untruncated ones. The direction is knowable even "
        "though the magnitude is not, and a reader told only the count will read the interval as "
        "unbiased-but-thinned, which it is not [Stage 11 §8.4]. Nothing is repaired: the alternative "
        "to the bound is a quantile over a mixture containing odds ratios of order 1e8, which is not "
        "a conservative interval but an uninterpretable one, and the choice among reporting options "
        "belongs to the PI. "
        "A p-value is reported on the interaction coefficient and on NOTHING else: [§13] asks for AN "
        "interaction test, and two level-specific p-values would be two tests of one contrast "
        "[Stage 11 §8.5]. `exp` is applied to the interval LIMITS and the coefficients are never "
        "exponentiated before the percentile, which is exact under this percentile method because it "
        "selects an order statistic and a strictly monotone transform commutes with order-statistic "
        "selection [Stage 11 §7.2]. No new interval kind, no new percentile method and no new floor: "
        "every one is `bootstrap.py`'s.")


def subgroups(df: pd.DataFrame, ps: propensity.Propensity, est: outcome.Primary,
              audit: Audit) -> Subgroups:
    """[§13]'s two subgroup estimates on the primary outcome, with one interaction test each.

    Takes no subgroup list: `C.SUBGROUPS` is the [§13] registry after the amendment of 2026-08-10
    withdrew the third, and a third restored there gains its estimate, its three keys, its interval
    and its p in one edit (§3.3).

    Takes `ps` and refits it per replicate over the WHOLE population — never within a level (§8.1).
    Takes `est` for `in_estimate`, which is the [§11] population the point-estimate fits run on, so
    the subgroup denominators are the primary's and not a fourth mask; each replicate recovers the
    same mask on its own drawn frame through `outcome.estimation_population`, which is the one
    definition both paths read (§8.5).

    Raises `C.SchemaError` H9 on the caller error `bootstrap.run`'s R-series covers, and H8 when a
    POINT-ESTIMATE `_subgroup_fit` raises — that call is outside the replicate loop and has nothing
    to be dropped into (§8.5). Inside the loop `model.FitError` never leaves it: G8, G9, O1-O6 and
    the propensity tokens are all droppable and are counted (§9).

    Returns no `bootstrap.Diagnostics` and that is §3.1's decision, not an omission: two `polr` fits
    per replicate cannot fill one scalar `n_alpha`.

    Appends three `model` entries (§10). Hypothesis-generating, and the label is a field (§3.1).
    """
    _assert_subgroup_inputs(df, ps, est)                              # H9

    estimates: dict[str, SubgroupEstimate] = {}
    for s in C.SUBGROUPS:                                             # declared order
        try:
            beta, gamma, _, mask = _subgroup_fit(df, ps, est.in_estimate, s)
        except model.FitError as failure:
            # §8.5. The droppable-and-counted argument is about the REPLICATE LOOP and this call is
            # not in one — there is nothing for the point estimate to be dropped into. Omitting the
            # subgroup from `estimates` is `pilots/analysis.py:811-813`'s uncounted `continue`, and
            # letting the `FitError` out contradicts §11's contract and hands the caller a token
            # whose bucket map is about replicates.
            raise C.SchemaError(
                f"H8  the POINT-ESTIMATE interaction fit for {s!r} raised: "
                f"{str(failure).splitlines()[0][:300]}\n  A subgroup whose point estimate cannot be "
                "fitted is a fact about the cohort and about [§13]'s choice of subgroups, not a "
                "sparse replicate. Inside the replicate loop this same token is dropped and counted; "
                "here there is nothing to drop into, so it raises and names the subgroup and the "
                "leading token [Stage 11 §8.5]. It is not substituted with a different estimator "
                "[roadmap invariant 5] and the subgroup is not silently omitted."
            ) from failure
        estimates[s] = _subgroup_estimate(df, ps, s, beta, gamma, mask)
        audit.record("model", f"{_STEP_SUBGROUP_FIT}_{s}", estimates[s].n,
                     _subgroup_fit_detail(estimates[s], s),
                     table=_subgroup_fit_table(estimates[s]))

    keys = tuple(f"{s}.{q}" for s in C.SUBGROUPS for q in _SUBGROUP_QUANTITIES)   # §3.3 — six
    collected = bootstrap.replicates(                                 # ONE loop, and it is [§10]'s
        df,
        lambda draw: _subgroup_replicate(draw, keys, audit.source),   # keys and source CLOSED OVER
        C.N_BOOT, C.SEED, C.BOOT_STRATUM,                             # no second seed
    )
    draws = bootstrap.collect(collected, keys)                        # the reconciling loop
    intervals = _intervals(draws, _tested_subgroup)                   # p on `.gamma` ONLY — §8.5

    # Tallied here rather than through `bootstrap.diagnostics`, because `Subgroups` carries no
    # `Diagnostics` (§3.1) and this is the one field of `Replicate` a two-fit body fills unambiguously.
    n_in_model: dict[int, int] = {}
    for replicate in collected:
        if replicate.n_in_model is not None:
            n_in_model[replicate.n_in_model] = n_in_model.get(replicate.n_in_model, 0) + 1

    audit.record("model", _STEP_SUBGROUP_REPLICATES, len(draws),
                 _subgroup_replicates_detail(draws, intervals),
                 table=_subgroup_replicates_table(draws, n_in_model))
    return Subgroups(estimates=estimates, seed=C.SEED, n_boot=C.N_BOOT, draws=draws,
                     intervals=intervals)


# ======================================================================================================
# §5 — the full-covariate arm [§13, DECISION 4]
# ======================================================================================================

def _arm_replicate(draw: pd.DataFrame, keys: tuple[str, ...], source: object) -> bootstrap.Replicate:
    """One [§13, DECISION 4] replicate: refit [§7] over PS_COVARIATES_FULL, then [§8], and nothing else.

    It is `bootstrap._replicate` with two things removed and one changed, and each is a decision:
    `outcome.secondary` is not called (§5.4), `_shared_design`'s invariant check is not run because
    no binary outcome is estimated here, and `propensity.fit` becomes `propensity.fit_full`.

    The `Audit` is constructed here and discarded, which is Stage 10 §6.2 unchanged: nine entries per
    replicate at `N_BOOT` would be memory that is a function of the replicate count, and no
    replicate's log is ever written.

    `C.SchemaError` is never caught. The only `except` is on `model.FitError`, which is Stage 10 §7.1
    as code and for its reason: a frame this stage constructed that cannot be read means this stage
    constructed it wrongly.

    The seven keys fail as ONE group. They come from one `polr` fit and a comparison across them is
    [§16]'s constant-shift statement, so they must share draws — Stage 10 §7.3, inherited rather than
    re-argued.
    """
    audit = Audit(source)
    try:
        ps = propensity.fit_full(draw, audit)               # the seam, INSIDE the loop
    except model.FitError as failure:
        return bootstrap.Replicate(
            values={}, failures={key: bootstrap.bucket(str(failure)) for key in keys},
            n_alpha=None, polr_iterations=None, sum_w=None, n_in_model=None)

    values: dict[str, float] = {}
    failures: dict[str, str] = {}
    n_alpha = polr_iterations = None
    try:
        est = outcome.primary(draw, ps, audit)              # UNCHANGED — §5.1
        values["beta"] = est.beta
        values.update({f"rd_{k}": v for k, v in est.rd.items() if f"rd_{k}" in keys})
        n_alpha, polr_iterations = len(est.fit.alpha), est.fit.iterations
    except model.FitError as failure:
        failures = {key: bootstrap.bucket(str(failure)) for key in keys}

    return bootstrap.Replicate(
        values=values, failures=failures, n_alpha=n_alpha, polr_iterations=polr_iterations,
        sum_w=float(ps.w[ps.in_model].sum()), n_in_model=int(ps.in_model.sum()))


def _assert_arm_inputs(df: pd.DataFrame, ps: propensity.Propensity) -> None:
    """H1 and H2, collected. Both are properties of the CALL and both are the caller's bug.

    **H1 is the REFERENCE half of the pair** (§4.5). `full_covariate` takes the PRIMARY's fit as the
    half of the pair the roadmap's Accept-when is about, and a caller passing the arm's own
    `Propensity` here would produce a record comparing the arm with itself — every cell finite, every
    ESS reconciling, and the printed clause a tautology. `Arm.__post_init__`'s H10 is the same check
    from the other side, on the record rather than on the call, and it is the one a caller reaches by
    building the record by hand.

    **H2 is R3's alignment class**, on the one function here that takes a frame and a `Propensity`
    together: `.loc[boolean_series]` returns rows in the SERIES' order, so a permuted-but-equal index
    pairs each record's weight with another record's fitted value and returns a different number with
    no raise.
    """
    bad: list[str] = []
    if ps.spec is not C.PROPENSITY_PRIMARY:
        bad.append(
            f"H1  `ps` carries the {ps.spec.label} specification {ps.spec.sap}, and this function "
            f"takes the {C.PROPENSITY_PRIMARY.label} fit as the REFERENCE half of the pair the "
            "roadmap's Accept-when is about — the arm's own fit is made HERE, by `fit_full`. Passing "
            "the arm's `Propensity` produces a record comparing the arm with itself, in which every "
            "cell is finite and the printed clause is a tautology [Stage 11 §4.5, §5.2]. Compared by "
            "IDENTITY and not by equality: the registry holds exactly two instances.")
    misaligned = [name for name, s in (("e", ps.e), ("w", ps.w), ("in_model", ps.in_model))
                  if not s.index.equals(df.index)]
    if misaligned:
        bad.append(
            f"H2  the Propensity is not aligned to this frame: {', '.join(misaligned)} "
            f"carr{'ies' if len(misaligned) == 1 else 'y'} a different index. `.index.equals` and "
            "NOT a length check or a set comparison: `.loc[boolean_series]` returns rows in the "
            "SERIES' order, so a permuted-but-equal index pairs each record's weight with another "
            "record's fitted value and returns a different number with no raise [Stage 9 §12]. This "
            "is `bootstrap.run`'s R3, on the one function here that takes a frame and a Propensity.")
    if bad:
        raise C.SchemaError(
            "\n  " + "\n  ".join(bad)
            + f"\n  {len(bad)} sensitivity-arm precondition(s) failed against {len(df)} records.")


def _arm_table(arm_ps: propensity.Propensity, arm_balance: balance.Balance,
               reference: propensity.Propensity, reference_balance: balance.Balance | None,
               draws: dict[str, bootstrap.Draws],
               diagnostics: bootstrap.Diagnostics, differs: bool) -> tuple[tuple[str, ...], ...]:
    """BOTH specifications' rows, side by side. The Accept-when is about a PAIR (§5.2).

    An entry showing one half makes the comparison the reader's arithmetic. Every number here is a
    field of the returned record — the surviving counts and buckets from `Arm.draws`, the cutpoint and
    iteration tallies and the `in_model` range from `Arm.diagnostics` — and none is recomputed at the
    print site.
    """
    worst_arm, undefined_arm = arm_balance.worst()
    header = ("quantity", C.PROPENSITY_PRIMARY.label, C.PROPENSITY_FULL.label, "note")
    rows: list[tuple[str, ...]] = [
        ("propensity covariates", str(len(reference.spec.covariates)),
         str(len(arm_ps.spec.covariates)),
         "the REQUESTED list and never the post-design survivors: in_model is complete_cases over "
         "it, so it is what decided the [§11] denominator [Stage 11 §4.2]"),
        ("design columns", str(len(reference.fit.columns)), str(len(arm_ps.fit.columns)),
         "dropped as constant: "
         f"{', '.join(reference.dropped) or 'none'} / {', '.join(arm_ps.dropped) or 'none'}"),
        ("in_model", str(int(reference.in_model.sum())), str(int(arm_ps.in_model.sum())),
         "the [§11] denominator of each fit"),
        ("the two in_model masks are IDENTICAL", "yes" if differs else "no", "—",
         "COMPUTED from the two masks and never asserted: the roadmap's Accept-when says the two "
         "estimates differ in the propensity specification AND IN NOTHING ELSE, which is true on a "
         "workbook where the added covariates are complete and false on one where they are not "
         "[Stage 11 §5.2]"),
    ]
    for code, label in C.TREATMENT_LABELS.items():
        rows.append((f"ESS ({label})", _fmt(reference.ess[code]), _fmt(arm_ps.ess[code]),
                     "the [§7] effective sample size per arm; it moves because the WEIGHTS move, "
                     "which is what a sensitivity analysis on the propensity specification is for"))
    if reference_balance is not None:
        worst_reference, _ = reference_balance.worst()
        rows.append((
            "worst residual abs SMD",
            "undefined" if worst_reference is None else
            f"{_fmt(abs(worst_reference.weighted))} ({worst_reference.covariate})",
            "undefined" if worst_arm is None else
            f"{_fmt(abs(worst_arm.weighted))} ({worst_arm.covariate})",
            "NOT a verdict: the summary statistic alone reads as uniformly better balance and it is "
            "not — see the per-row block below [Stage 11 §5.3]"))
        rows.append((
            f"rows reaching abs SMD {_fmt(C.SMD_THRESHOLD)}",
            str(len(reference_balance.unbalanced())), str(len(arm_balance.unbalanced())),
            ", ".join(reference_balance.unbalanced()) + " / " + ", ".join(arm_balance.unbalanced())))
        before = {row.covariate: row.weighted for row in reference_balance.covariates}
        for row in arm_balance.covariates:
            if row.covariate in C.NEGATIVE_CONTROLS or row.covariate.startswith("center = "):
                rows.append((f"abs SMD after weighting: {row.covariate}",
                             _fmt(before[row.covariate]), _fmt(row.weighted),
                             "SPENT as a negative control by this arm [§6, DECISION 4]"
                             if row.covariate in C.NEGATIVE_CONTROLS else
                             "in BOTH models, and it moves — overlap weights balance an in-model "
                             "covariate exactly only under maximum likelihood, and Firth solves the "
                             "score equations PLUS a penalty whose deviation is largest exactly "
                             "where the design is closest to separated [§13]"))
    rows.append((
        "rows no SMD could be computed on", "—", str(len(undefined_arm)),
        ", ".join(undefined_arm) if undefined_arm else "none"))

    buckets = sorted(set(C.FAILURE_BUCKETS.values()))
    for key, d in draws.items():
        rows.append((
            f"replicates surviving: {key}", "—", f"{len(d.draws)} of {d.n_attempted}",
            ", ".join(f"{bucket} {d.failures.get(bucket, 0)}" for bucket in buckets)))
    rows.append((
        "cutpoints fitted", "—",
        ", ".join(f"{k}: {v}" for k, v in diagnostics.n_alpha.items()),
        "[§16]'s constant-shift statement, computed and not written: a replicate fitting five "
        "cutpoints rather than six fitted a collapsed response"))
    rows.append((
        "polr iterations", "—",
        ", ".join(f"{k}: {v}" for k, v in diagnostics.polr_iterations.items()),
        "reported AS A PAIR with sum_w, which is Stage 10 §7.5's finding rather than its "
        "reassurance"))
    rows.append((
        "in_model per replicate", "—",
        (f"{min(diagnostics.n_in_model)} .. {max(diagnostics.n_in_model)}"
         if diagnostics.n_in_model else "missing"),
        "the upper end is the replicates that drew the covariate-incomplete record zero times"))
    rows.append((
        "augmented estimates", "—", str(len(diagnostics.max_abs_beta)),
        "EMPTY by construction: [§13] gives this arm the primary only, so there is no m_a(X) and "
        "nothing augmented [Stage 11 §5.4]"))
    return (header, *rows)


def _arm_detail(arm_ps: propensity.Propensity, reference: propensity.Propensity,
                arm_balance: balance.Balance, draws: dict[str, bootstrap.Draws],
                intervals: dict[str, bootstrap.Interval], differs: bool) -> str:
    worst, _ = arm_balance.worst()
    surviving = min((len(d.draws) for d in draws.values()), default=0)
    attempted = max((d.n_attempted for d in draws.values()), default=0)
    return (
        f"[§13, DECISION 4] full-covariate propensity sensitivity analysis: the [§7] fit refitted "
        f"over {len(arm_ps.spec.covariates)} covariate(s) instead of "
        f"{len(reference.spec.covariates)}, its own overlap weights, its own [§9] balance table and "
        f"its own bootstrap, reported BESIDE the primary and never instead of it. "
        "**It spends the negative controls, and that is said where the balance table is rendered.** "
        "`C.NEGATIVE_CONTROLS` is computed as PS_COVARIATES_FULL minus PS_COVARIATES, so the "
        f"{len(C.NEGATIVE_CONTROLS)} covariate(s) this arm adjusts for are exactly the ones that stop "
        "being controls in it — after this arm the only covariate in no propensity model is "
        "`penumbra_ml`, and [§13] accepts that: the controls have already served their purpose. "
        "**A lower worst residual |SMD| here is not uniformly better balance**: the risk factors "
        "improve and centre gets WORSE, which is [§13]'s own mechanism running in the direction [§13] "
        "describes — overlap weights balance an in-model covariate exactly only under maximum "
        "likelihood, Firth solves the score equations plus a penalty, and the deviation is largest "
        "exactly where the design is closest to separated. A wider design on the same records is "
        "closer to separated, so the deviation on the near-determining covariate grows. That is a "
        "trade-off made visible from the other side and it is not grounds to prefer either "
        "specification [Stage 11 §5.3]. "
        + ("The two `in_model` masks are IDENTICAL on this workbook, so the estimates differ in the "
           "propensity specification and in nothing else — and that clause is COMPUTED from the two "
           "masks rather than written, because it is a property of THIS workbook and not of the "
           "design: on a workbook where any added covariate is missing on any row the two "
           "populations diverge and the clause becomes false while nothing in the code changes "
           "[Stage 11 §5.2]. " if differs else
           "The two `in_model` masks DIFFER, so the estimates differ in the propensity specification "
           "AND in what that specification's completeness costs — both counts are printed above and "
           "the roadmap's Accept-when clause does not hold on this workbook [Stage 11 §5.2]. ")
        + f"The worst residual |SMD| after weighting is "
        + ("undefined on every row" if worst is None else
           f"{_fmt(abs(worst.weighted))} on {worst.covariate}")
        + f". The bootstrap is `bootstrap.replicates` over a primary-only body and NOT a second "
        f"`bootstrap.run`, and that is a departure this entry records rather than takes quietly: "
        "`run`'s replicate body calls `propensity.fit` unconditionally — R9 now makes that a raise "
        "rather than a silent mislabel — and `run` also drives `outcome.secondary`, which would "
        "produce nineteen more intervals under the [§13] specification that no prespecified section "
        "asks for. A table of unreported estimates is a multiplicity waiting to be read, and [§13]'s "
        f"own Benjamini-Hochberg is scoped to seven p-values that are not these [Stage 11 §5.4]. "
        f"The worst-off estimand keeps {surviving} of {attempted} replicate(s) and "
        f"{len(intervals)} of {len(draws)} estimand(s) carry an interval, against a floor of "
        f"{C.ci_min_draws()}. **NO p-value is reported on any of the {len(draws)} keys**: [§13] "
        "prescribes sensitivity analyses that REPORT, each with its ESS and worst residual |SMD|, and "
        "prescribes no test of any of them; a p on this arm's `beta` would be an eighth test [§13]'s "
        "correction is not told about, and the primary's own p is [§10]'s and is uncorrected by "
        "prescription [Stage 11 §5.5]. The six cumulative RD_k carry no p here for the reason they "
        "carry none in the primary: [§8] gives them intervals only. The seed, the replicate count and "
        "the stratum are the prespecified ones and there is no second seed, so this run draws the "
        "SAME sequence of resampled frames as [§10]'s — a property that is asserted and never used, "
        "since nothing computes a contrast between an arm draw and a primary draw [Stage 11 §5.4].")


def full_covariate(df: pd.DataFrame, ps: propensity.Propensity, audit: Audit) -> Arm:
    """The [§13, DECISION 4] full-covariate propensity sensitivity analysis.

    Takes no covariate list, for the reason `propensity.fit` takes none: the specification is
    `C.PROPENSITY_FULL` and it is NAMED, not passed (§4). `ps` is the PRIMARY's fit and is taken as
    the reference half of the pair the roadmap's Accept-when is about — it is read for its
    `in_model`, its `ess` and its `spec`, and it is never refitted here.

    Four things about the order below, each of which is a failure if moved:

      * `_assert_arm_inputs` FIRST — H1 and H2 — so a caller passing the ARM's `Propensity` as `ps`
        reports it rather than producing a record comparing the arm with itself (§4.5).
      * `balance.assess` BEFORE the bootstrap, so the log carries the 19 re-roled rows even if the
        loop then raises — `propensity.py`'s own rule, two stages on.
      * `differs_only_in_specification` is computed from the two `in_model` masks AFTER both exist
        and is never a literal (§5.2).
      * The loop is `bootstrap.replicates`' and not this function's, for Stage 10 §11's reason: a
        second loop would put the seeding and the ordering in two places.

    `keys` is built BEFORE `replicates` and not inside the body, which is Stage 10's stated reason
    rather than a style: §7.1 requires a `propensity.fit` failure counted against EVERY key, and a
    body that named its own keys could not count a failure that happened before it named them. The
    six `rd_k` come off the ARM's own `Primary` and not from `C.MRS_THRESHOLDS`, which is
    `bootstrap._estimand_keys`' argument taken rather than re-derived: if a future cohort fitted a
    different set of thresholds under one specification than the other, the mismatch surfaces as a
    key that is absent rather than as a `Draws` reconciliation nobody can explain.

    Appends ten `model` entries (§10) — four from `fit_full`, two from `assess`, three from `primary`,
    all under `_full_covariate`-suffixed names, and one of its own.
    """
    _assert_arm_inputs(df, ps)                                   # H1, H2

    ps_full = propensity.fit_full(df, audit)                     # [Stage 6 §11], through §4's seam
    bal_full = balance.assess(df, ps_full, audit)                # [Stage 7 §11], roles off ps.spec
    est_full = outcome.primary(df, ps_full, audit)               # [Stage 8 §11], UNCHANGED

    keys = ("beta", *(f"rd_{k}" for k in est_full.rd))           # §3.3 — off the ARM's Primary
    collected = bootstrap.replicates(
        df,
        lambda draw: _arm_replicate(draw, keys, audit.source),   # keys and source CLOSED OVER
        C.N_BOOT, C.SEED, C.BOOT_STRATUM,                        # no second seed — §5.4
    )
    draws = bootstrap.collect(collected, keys)                   # the reconciling loop
    diagnostics = bootstrap.diagnostics(collected)               # ONE polr fit per replicate — §3.1
    intervals = _intervals(draws, lambda key: False)             # no p on any arm key — §5.5

    # COMPUTED from the two masks, after both exist, and never a literal (§5.2). The roadmap's
    # Accept-when is a statement about this workbook and not about the design.
    differs = bool(ps.in_model.equals(ps_full.in_model))

    # The primary's balance table, recomputed here over a THROWAWAY `Audit` so the pipeline's log
    # gains no duplicate entries. §10 requires this entry to carry BOTH specifications' rows side by
    # side, because the Accept-when is a statement about a PAIR and an entry showing one half makes
    # the comparison the reader's arithmetic — and `full_covariate`'s signature is `(df, ps, audit)`,
    # so the primary's `Balance` is not passed in. Stage 7 §8 already prices a second `_table` call
    # at 7.4 ms and takes it for the same kind of reason: a table that agrees with itself by
    # construction has stopped being evidence. The COUNTERS in this entry are not recomputed — every
    # one comes off `Arm.draws` and `Arm.diagnostics`.
    reference_balance = balance.assess(df, ps, Audit(audit.source))
    audit.record("model", _STEP_ARM, int(ps_full.in_model.sum()),
                 _arm_detail(ps_full, ps, bal_full, draws, intervals, differs),
                 table=_arm_table(ps_full, bal_full, ps, reference_balance, draws, diagnostics,
                                  differs))
    return Arm(ps=ps_full, balance=bal_full, estimate=est_full, reference=ps, seed=C.SEED,
               n_boot=C.N_BOOT, draws=draws, intervals=intervals, diagnostics=diagnostics,
               differs_only_in_specification=differs)
