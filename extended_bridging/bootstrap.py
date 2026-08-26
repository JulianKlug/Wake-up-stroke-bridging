"""Stage 10 — the [§10] bootstrap engine.

One interval and, where [§10] prescribes one, one p-value for every point estimate Stages 8 and 9
produce: the primary common odds ratio, the six cumulative ``RD_k``, and the seven binary outcomes'
risk differences, marginal odds ratios and model-assisted augmented risk differences. Patient-level
nonparametric bootstrap, stratified by centre, ``N_BOOT`` replicates, recorded seed, percentile
limits. Plus the three diagnostics Stages 8 and 9 asked this stage to report because they are
answerable only from replicates.

Section references in brackets are to ``statistical_analysis_plan.md``. The specification for this
module is ``specs/stage10_bootstrap_engine.md``; nothing here is invented outside it, and the four
places where it could not be followed as written are marked ``[SPEC]`` and argued at the point of
departure.

**Why this is a module and not an extension of ``outcome.py``** [Stage 10 §0.1]. Stage 10 REFITS
Stages 6, 8 and 9, so it imports all three; a bootstrap living inside ``outcome.py`` would put the
replicate loop in the same module as the estimators it resamples, which is the arrangement in which a
future edit can reach from the loop into an estimator's private. ``_augmentable`` and
``_augmented_path`` are ``outcome.py`` privates and §6.3's whole rule is that this stage **reads the
decided path and never re-derives it** — a module boundary is what makes that checkable by scan
(§15.11) rather than by review. And Stages 12 and 13 [§14a, §14b] need this stage's resampler and its
percentile machinery and none of its estimators, so the four general functions take the replicate
body as a callable.

::

    cohort.build(...)          -> DataFrame, one row per patient              [Stage 5 §11]
    propensity.fit(cohort, a)  -> Propensity: e, w, in_model                  [Stage 6 §11]
    outcome.primary(...)       -> Primary:   beta, alpha, rd, cumulative      [Stage 8 §11]
    outcome.secondary(...)     -> Secondary: seven BinaryEstimate             [Stage 9 §11]
    balance.assess(...)        -> Balance   — NOT read here                   [Stage 7 §11]
                       │
                       ▼
    ┌──────────────────────────────────────────────────────────────────────────────────┐
    │  run(cohort, ps, est, sec, audit)                                                │
    │    ├─ _assert_run_inputs(...)          R1-R8, collected            (§4.4)        │
    │    ├─ paths = {k: e.augmented_path}    READ off the point estimate  (§6.3)        │
    │    ├─ keys  = _estimand_keys(est, sec) twenty-six of them           (§3.3)        │
    │    ├─ replicates(...)                  ONE Generator, in order      (§4.3)        │
    │    │     resample(...) -> _replicate(...) -> Replicate                            │
    │    ├─ _collect / _diagnostics          per estimand GROUP           (§7.1, §7.3)  │
    │    ├─ _record_replicates(...)          one `model` entry            (§10.1)       │
    │    └─ percentile_ci / bootstrap_p      per estimand                 (§8, §9)      │
    └──────────────────────────────────────────────────────────────────────────────────┘
                       │
                       ▼
      Bootstrap(seed, n_boot, draws, intervals, diagnostics)  →  Stages 11-14

``run`` appends ONE entry of the ``model`` kind, taking the ledger 30 → 31. It writes no file, and the
gitignored audit log is the only place any interval on this workbook exists [Stage 10 §4.5].

**What this module does not touch** [Stage 10 §0.2]. ``balance`` is not imported and no interval is
put on an SMD. ``propensity.py`` is not amended — §5's central claim — and none of its privates is
reached into. The estimators are called and not edited. And the point estimate is READ and never
recomputed: a percentile limit is an order statistic of the replicate distribution and is not a
function of the point value.

This module is **not** exempt from the Stage 1 §7 raw-name scan and must never become exempt.

--------------------------------------------------------------------------------------------------
THE FOUR DEPARTURES FROM THE SPECIFICATION, each measured and each marked ``[SPEC]`` below.

1. **``_estimand_keys`` is called BEFORE the loop and ``_replicate`` takes ``keys``.** §11's fence
   builds the keys after ``replicates`` returns and passes ``_replicate(draw, paths, source)``, but
   §7.1 requires that a ``propensity.fit`` failure be counted against *every* key — which
   ``_replicate`` cannot do without knowing them. Both `est` and `sec` are in hand at the top of
   ``run``, so the call moves up and the lambda gains an argument. Nothing else about §11's order
   changes.

2. **§6.4's shared design is a per-replicate INVARIANT CHECK and its 27.2 s saving is not
   realised.** §6.4 wants one design built per replicate instead of four; §7.3 requires the seven
   binary outcomes be driven through ``outcome.secondary`` rather than individually — *"re-driving
   the seven-outcome loop in `bootstrap.py` would put a copy of Stage 9's ordering, masking and path
   logic in a second module"* — and ``secondary`` builds its own designs and takes no way to be
   given one. The two sections cannot both be honoured without a second parameter §13 does not
   license. §7.3's is the INFERENCE decision and wins; ``_shared_design`` keeps everything §15.9
   asserts of it — the mask equality, the single build, and the loud ``SchemaError`` when roadmap
   invariant 6 stops holding — and costs about 4.5 ms per replicate instead of saving 13.6.

3. **``C.FAILURE_BUCKETS`` carries nineteen tokens' worth of raise site and not sixteen.**
   §7.2's scan scope is ``model.py`` and ``outcome.py``; ``propensity.py`` raises ``model.FitError``
   three more times, with the tokens ``ESS:`` and ``F5``, both on ``propensity.fit``'s own path —
   which is precisely the path §7.1 says fails a whole replicate. `config.py` carries the argument.

4. **``Secondary`` gains a ``failures`` field.** §7.3 specifies ``collect`` and requires that
   ``collect=True`` return *"six estimates and one recorded failure"*, and §15.7 asserts the failure
   is RECORDED rather than absent — which needs somewhere to record it. `outcome.py` carries it.
"""
from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass, field

import numpy as np
import pandas as pd

import config as C
import model
import outcome
import propensity
# `_fmt` is private to data.py and is imported anyway, for the reason derive.py, propensity.py,
# balance.py and outcome.py give: it is the pipeline's *one* float formatter, and a second one is a
# second way for two runs to disagree. Every number below reaches the log through it.
from data import Audit, _fmt


# The three declared augmentation paths, as a literal, and §15.11 is why it is one. Stage 9 declares
# them inside `outcome._augmented_path`, and R7 has to check a `Secondary` against them — but §6.3's
# whole rule is that this stage READS the decided path and never re-derives it, and §15.11 asserts by
# AST scan that this module names neither `_augmented_path` nor `_augmentable`. So the three strings
# are written here, where R7 can compare against them, and nowhere else in this module.
_DECLARED_PATHS: tuple[str, ...] = ("full", "reduced", "unaugmented")

# The audit step name. A constant rather than a literal at the one `audit.record` call site, for
# `outcome.py`'s reason (outcome.py:121-123): the acceptance suite asserts the entry appears under
# this name, and a test comparing against a string literal it also writes is a test of nothing.
_STEP_REPLICATES: str = "bootstrap_replicates"


# --- what the stage returns [§3.1] ------------------------------------------------------------------

@dataclass(frozen=True)
class Replicate:
    """What ONE replicate produced: every estimand it reached, and every one it lost.

    **`values` and `failures` PARTITION the attempted keys** — a key is in exactly one of them,
    never both and never neither. That partition is the whole of §15.4's reconciliation property:
    `_collect` sums over replicates key by key, and if a replicate can lose a key silently the
    three numbers in `Draws` stop adding up. A dict-of-value plus a dict-of-bucket makes the
    partition a property of the type rather than of the loop that fills it — which is what
    `__post_init__` below turns from a sentence into a raise.

    It is declared here and not left to the implementation because `_collect` and `_diagnostics`
    both destructure it, and every other structure this stage carries is a documented frozen
    dataclass. This one carries every number the stage reports.

    The five diagnostic fields are `None` when the replicate died before `outcome.primary` ran —
    which is exactly the propensity-fit failure §7.1 describes — so a `None` is "we never got
    there" and is dropped from the distribution rather than counted as a zero.

    `max_abs_beta` holds entries only for outcomes with an `m_a(X)`, so `sich` and `ph2` are
    ABSENT rather than present with a nan (§10.2). `or_corrected` is keyed by every outcome that
    produced an odds ratio, because the correction is a property of the odds ratio (§9.2).
    """

    values:          dict[str, float]   # estimand key -> the draw            — §3.3's keys
    failures:        dict[str, str]     # estimand key -> C.FAILURE_BUCKETS value — §7.2
    n_alpha:         int | None         # len(fit.alpha)                      — §7.4
    polr_iterations: int | None         # the primary fit's iteration count   — §7.5
    sum_w:           float | None       # Σw over in_model                    — §7.5
    n_in_model:      int | None         # |in_model|                          — §7.5
    max_abs_beta:    dict[str, float] = field(default_factory=dict)   # AUGMENTED ONLY — §10.2
    or_corrected:    dict[str, bool] = field(default_factory=dict)    # §6.3's correction fired

    def __post_init__(self) -> None:
        overlap = sorted(set(self.values) & set(self.failures))
        if overlap:
            raise C.SchemaError(
                f"a replicate reports {len(overlap)} estimand(s) as BOTH a draw and a failure: "
                f"{', '.join(overlap)}. `values` and `failures` partition the attempted keys "
                "[Stage 10 §3.1]: a key counted twice makes `len(draws) + sum(failures)` exceed "
                "`n_attempted`, and §15.4's reconciliation is the one property no wrong "
                "implementation satisfies by accident.")


@dataclass(frozen=True)
class Draws:
    """One estimand's replicate draws, and the counters that explain its denominator.

    **`draws` and `n_attempted` are BOTH fields and neither is derivable from the other**, which is
    the whole of [§10]'s "dropped and counted": `len(draws)` is what the percentiles were taken over
    and `n_attempted` is what was asked for, and a reader given only the first cannot tell a
    complete interval from a thinned one. Stage 14 prints both beside every interval (§12.3).

    `failures` is keyed by BUCKET and not by message, and §7.2 is why there is more than one key:
    Stage 8 §11 requires the separation count reported separately from the convergence count,
    because it measured that the second never fires on this estimator and a single counter reading
    zero therefore says nothing. The buckets are C.FAILURE_BUCKETS' values and the mapping from a
    raised FitError to one of them is §7.2's, asserted by scan in §15.6 rather than trusted.

    **The sum of `failures.values()` plus `len(draws)` equals `n_attempted`, and that is an
    assertion and not a comment** (§15.4). A replicate is either a draw or exactly one counted
    failure; an implementation that drops silently anywhere makes the three numbers stop reconciling,
    which is the one property of this dataclass that no wrong implementation satisfies by accident.

    **`n_attempted` is `C.N_BOOT` for every key this stage produces**, and §7.1 is why: a
    `propensity.fit` failure fails the whole replicate and is counted against EVERY estimand rather
    than making those replicates disappear from every denominator uncounted. The field is kept
    rather than derived because `replicates` is general (§12.2) and Stage 12's bodies may attempt
    an estimand selectively; on [§10]'s body the equality `n_attempted == C.N_BOOT` is itself an
    assertion (§15.4) and not a definition.
    """

    quantity: str                  # "beta", "rd_2", "sich.rd", ... — §3.3 enumerates all 26
    draws: np.ndarray              # (n,) float64, finite, in replicate order — §7.1
    n_attempted: int               # replicates in which this estimand was ATTEMPTED — §7.1, §7.3
    failures: dict[str, int]       # bucket -> count; sums with len(draws) to n_attempted — §7.2

    def __post_init__(self) -> None:
        reconciled = len(self.draws) + sum(self.failures.values())
        if reconciled != self.n_attempted:
            raise C.SchemaError(
                f"{self.quantity}: {len(self.draws)} draw(s) plus {sum(self.failures.values())} "
                f"counted failure(s) is {reconciled} against {self.n_attempted} attempted. A "
                "replicate is either a draw or exactly one counted failure, and an implementation "
                "that drops silently anywhere makes these three numbers stop reconciling "
                "[Stage 10 §3.1, §15.4].")


@dataclass(frozen=True)
class Interval:
    """A percentile interval, its level, and the p-value where [§10] prescribes one.

    **`p` is `None` where [§10] prescribes no test, and that is a field rather than an omission.**
    [§8] gives the six cumulative RD_k intervals and NO p-values, so `p is None` on all six, and
    Stage 14 printing a p-value for one of them is then a `None` reaching a formatter rather than a
    number nobody questioned (§9.4).

    **`method` travels with the interval.** §8.2 measured that numpy's default percentile definition
    disagrees with [§10]'s p-value at exactly the boundary that decides significance, so which
    definition produced these two numbers is part of what they are. It is pinned in `config.py`, so
    this field records rather than chooses.
    """

    lo: float
    hi: float
    level: float                   # C.CI_LEVEL
    method: str                    # C.PERCENTILE_METHOD — §8.2
    n_draws: int                   # the surviving draws these limits are order statistics of
    p: float | None                # §9; None where [§8] prescribes intervals only


@dataclass(frozen=True)
class Diagnostics:
    """The three things Stages 8 and 9 asked this stage to report, and one it asked for itself.

    None of these is an estimate and none gets an interval. They exist because a counter reading
    zero is not evidence of anything unless something beside it says what the replicates looked
    like — which is Stage 8 §11's argument for `len(fit.alpha)` and Stage 9 §9.6's for `max|beta|`,
    made twice for the same reason.

    `sum_w` is this stage's own addition and Stage 8 §11 is why: it wrote that a replicate's Σw was
    not measured there, that nothing was expected to move, and that the visible symptom if it did
    would be `polr`'s iteration count dropping. Both are measured in §7.5 and both moved, so the
    pair is reported rather than the reassurance.

    **`or_corrected` is a FIELD and not only a rendered cell**, and §9.2 is the reason: Stage 9 §6.3
    requires the correction named *"wherever the estimate appears"*, §12.3 establishes that a
    bootstrap interval is somewhere it appears, and the only interface that satisfies that rule is a
    typed one. Every other number Stage 14 needs travels as a field; routing this one through
    `data._md_table` would make Stage 14 parse a markdown cell it printed itself. It is still
    rendered in the grid (§10.1) — the field is what Stage 14 reads.

    Every distribution here is over the replicates that REACHED it. `n_alpha`, `polr_iterations`,
    `sum_w` and `n_in_model` therefore sum to the live-replicate count and not to `C.N_BOOT`; the
    replicates §7.1 kills before `outcome.primary` contribute a `None` and are dropped (§15.15).
    """

    n_alpha: dict[int, int]               # cutpoint count -> replicates — §7.4
    polr_iterations: dict[int, int]       # iteration count -> replicates — §7.5
    sum_w: np.ndarray                     # (n,) Σw over in_model per replicate — §7.5
    max_abs_beta: dict[str, np.ndarray]   # outcome -> max|beta| of m_a(X); AUGMENTED ONLY — §10.2
    n_in_model: dict[int, int]            # in_model size -> replicates — §7.5
    or_corrected: dict[str, int]          # outcome -> replicates the correction fired in — §9.2
                                          #   denominator is that outcome's Draws.n_attempted


@dataclass(frozen=True)
class Bootstrap:
    """Every interval [§10] prescribes, the draws behind them, and the diagnostics beside them.

    **The seed and the replicate count are fields because [§16] requires them reported**, and they
    are recorded here rather than looked up from `config.py` at print time: a run summary that reads
    its seed from the configuration reports the configuration and not the run.

    There is no `favours_bridging()` and no `significant` field, for Stage 8 §3's reason carried one
    stage on. [§10] closes with "estimation, not testing, is the reportable output"; a boolean here
    would be the dichotomisation that sentence declines, computed once and then quoted forever.
    """

    seed: int
    n_boot: int
    draws: dict[str, Draws]           # keyed as §3.3's estimand keys
    intervals: dict[str, Interval]    # the same keys
    diagnostics: Diagnostics


# --- the resample [§10, §5.3] -----------------------------------------------------------------------

def resample(df: pd.DataFrame, rng: np.random.Generator, stratum: str) -> pd.DataFrame:
    """One [§10] replicate: patient-level draw with replacement WITHIN each stratum.

    General in its frame and its stratum from the first line, because [§14a] and [§14b] prescribe
    the same resampling over a population this stage never sees (§12.2). It names no covariate, no
    outcome and no centre, and `stratum` is a parameter rather than `"center"` for that reason.

    Three things about it, each of which is a failure if changed:

      * THE STRATUM TOTALS ARE FIXED. `size=len(g)` per stratum, never a draw over the whole frame,
        because [§10]'s inference is "conditional on the participating centres" and a replicate that
        lost a centre would not be (§4.1).
      * THE ORDER IS THE STRATUM COLUMN'S OWN SORT ORDER -- `groupby(sort=True)` -- and NOT frame
        order. The draw is a function of the seed AND of the order the strata are visited in; frame
        order would make it a function of how the workbook happened to be sorted, which is the
        property Stage 2's shuffle test exists to deny (§4.3). **There is no second sort on top of
        `groupby`'s**, and an earlier draft had one keyed on `str(label)`: that is a no-op on this
        study's `string` centre labels and WRONG for a general resampler, because it visits an
        integer stratum column as 1, 10, 2, 9 -- measured. `resample` is general (§12.2), so the
        contract is the column's own ordering and pandas is the one thing defining it.
      * EACH DRAWN ROW GETS A DISTINCT case_id, and §5.2 is the whole argument. The k-th appearance
        of a patient becomes `case_id#k` counting from 1, so the FIRST appearance is `id#1` and not
        the bare id -- uniform, because a mixed scheme makes "was this row drawn once" a question
        about string formatting.

    The index is reset before the rename. **This is NOT because a duplicated index breaks anything**
    -- measured, `propensity.fit`, `outcome.primary` and `outcome.secondary` all run correctly on a
    replicate with 59 distinct labels over 93 rows, because every read of `e` and `w` is boolean
    (§3.4). It is reset because this function is general and Stages 12 and 13's bodies are unwritten,
    and a frame with a duplicated index is a frame on which a future label-based `.loc` is wrong in a
    way that returns a number.

    Deterministic given `rng`'s state. Adds no column, and the frame it returns has the same columns
    and dtypes as the frame it was given -- asserted in §15.2, because a resampler that quietly
    changes a dtype changes `model.design` (TODOS' Stage 1 dtype item).

    **THE CAST BACK TO `case_id`'S OWN DTYPE IS NOT COSMETIC, AND WITHOUT IT §15.2 FAILS.**
    `case_id` is declared `dtype="string"` in `COLUMN_CONTRACT`, and assigning a plain Python list
    to a column returns `object` -- measured, `string[python]` in and `dtype('O')` out on the fence
    as an earlier draft wrote it. `data.py:442` already carries this warning verbatim for the same
    reason one stage up: *"The cast back to `string` is not cosmetic: `.map` returns object dtype"*.
    The cast reads the INPUT frame's dtype rather than naming `"string"`, because this function is
    general and Stage 12's identifier column is not this study's.
    """
    parts = [
        g.iloc[rng.integers(0, len(g), size=len(g))]
        for _, g in df.groupby(stratum, sort=True, observed=True)
    ]
    out = pd.concat(parts, axis=0).reset_index(drop=True)
    occurrence = out.groupby("case_id", sort=False, observed=True).cumcount() + 1
    out["case_id"] = pd.Series(
        [f"{cid}#{k}" for cid, k in zip(out["case_id"], occurrence)],
        index=out.index,
    ).astype(df["case_id"].dtype)
    return out


def replicates(df: pd.DataFrame, body: Callable[[pd.DataFrame], object],
               n: int, seed: int, stratum: str) -> tuple[object, ...]:
    """`n` replicates of `body` over stratified resamples of `df`. ONE Generator, consumed in order.

    General in its body, which is §0.1's argument as a signature: [§14a] and [§14b] resample the same
    way over machinery that fits no propensity model at all, so the replicate body is a callable and
    this function knows nothing about [§7] or [§8] (§12.2).

    **ONE `numpy.random.Generator` is created once and consumed by every replicate in order**, rather
    than one seeded per replicate: a per-replicate seed derived from `b` is reproducible too, but it
    makes the draw a function of an index that a later edit to the loop can renumber, and the failure
    is silent (§4.3). `np.random.default_rng` and not `RandomState`: the legacy class is not
    deprecated but its stream is a compatibility guarantee rather than a design.

    **`resample` is called OUTSIDE `body` and there is no `try` here at all**, which is what makes
    the stream a function of the seed alone and not of whether a replicate succeeded. §15.2 asserts
    it by driving this function with a body that returns and again with one that raises, and
    requiring the same sequence of drawn frames: an edit that moved the draw inside a failure path
    would make the seed stop identifying the replicates, and nothing else in §15 would notice.

    Whatever `body` raises leaves this function. §7.1's taxonomy is the BODY's to apply — `run`'s
    body catches `model.FitError` per estimand group and never catches `C.SchemaError` — because a
    general loop that decided which failures were droppable would be deciding it for Stages 12 and
    13 as well.
    """
    rng = np.random.default_rng(seed)
    return tuple(body(resample(df, rng, stratum)) for _ in range(n))


# --- the percentile interval and the p-value [§8, §9, §11] -------------------------------------------

def percentile_ci(draws: np.ndarray, level: float = C.CI_LEVEL) -> tuple[float, float]:
    """The [§10] percentile limits, at the ONE definition §8.2 pins.

    `method=C.PERCENTILE_METHOD` is passed explicitly and is never left to default: numpy's default
    is "linear", which interpolates between order statistics and disagrees with `bootstrap_p` below
    at exactly the tail count that decides significance -- measured, 4.707% of constructed draw sets
    (§8.2). The two functions are one decision and this argument is where it is recorded.

    Raises on fewer than `C.ci_min_draws(level)` draws rather than returning an extreme under a
    percentile's name (§8.3). `run` checks the count before calling, so the raise is a caller bug.
    The call is QUALIFIED: `ci_min_draws` lives in `config.py`, which this module imports as `C`,
    and an earlier draft of this fence called it bare -- a NameError in the shipped module.

    **THE QUANTILES ARE SNAPPED AND `100.0 * ((1.0 - level) / 2.0)` IS WRONG.** That expression
    returns 2.500000000000002 at level 0.95, not 2.5, because `1.0 - 0.95` is
    0.050000000000000044. The excess is 2e-15 and it is not cosmetic: `inverted_cdf` takes the order
    statistic at index ceil(q/100 * n), so at n = 2000 the index moves from 50 to 51 and the limit
    moves from the largest negative draw to the smallest positive one. Measured on one draw set with
    a tail count of 50: `-0.107095` at the literal 2.5 against `+0.101248` at the computed quantile
    -- opposite signs, so opposite verdicts, at exactly the tail count §8.2 chose this method for.
    An earlier form of this fence had it, and §21b is where running the fence found it.
    """
    if len(draws) < C.ci_min_draws(level):
        raise C.SchemaError(
            f"percentile_ci: {len(draws)} draw(s) against a floor of {C.ci_min_draws(level)} at level "
            f"{level:g}. Below the floor the lower limit is the sample minimum and np.percentile "
            "returns it while still calling it a percentile [Stage 10 §8.3].")
    q_lo = round(50.0 * (1.0 - level), 9)          # 2.5 EXACTLY at level 0.95 -- see the docstring
    q_hi = round(100.0 - q_lo, 9)                  # 97.5 exactly
    lo, hi = np.percentile(draws, [q_lo, q_hi], method=C.PERCENTILE_METHOD)
    return float(lo), float(hi)


def bootstrap_p(draws: np.ndarray) -> float:
    """[§10]'s two-sided bootstrap p, floored at 1/(B+1) where B is the SURVIVING draw count.

    `B` is `len(draws)` and never C.N_BOOT: the floor is the smallest p a bootstrap of this size can
    express, and using N_BOOT on a thinned set claims a resolution the draws do not have. On this
    cohort nothing was dropped from the primary, so the two are equal -- which is exactly the
    condition under which an implementation reading N_BOOT is green, and why §15.8 asserts it on a
    fixture where they differ.

    Draws exactly equal to 0.0 are counted in BOTH tails, so `Pr(<=0) + Pr(>=0) > 1` and the test is
    conservative. That is deliberate and measured (§8.2): the alternative -- splitting ties -- makes
    p a function of a tie-breaking rule nobody prespecified.
    """
    p_le = float(np.mean(draws <= 0.0))
    p_ge = float(np.mean(draws >= 0.0))
    return max(2.0 * min(p_le, p_ge), 1.0 / (len(draws) + 1))


# --- the estimand keys [§3.3] -------------------------------------------------------------------------

def _estimand_keys(est: outcome.Primary, sec: outcome.Secondary) -> tuple[str, ...]:
    """§3.3's twenty-six keys, in replicate-body order: the primary group, then outcome by outcome.

    It takes `est` because the six `RD_k` keys live on the POINT ESTIMATE and not on the
    configuration: `Primary.rd` is `dict[int, float]` and its keys are what the primary actually
    fitted. An earlier draft of the specification passed only `sec`, which left `est` an unused
    parameter of `run` and the six keys unsourced (§3.3, §21c item 5).

    The order is `_replicate`'s own — `beta`, the `rd_<k>`, then each outcome's `rd`, `odds_ratio`
    and `augmented` in `C.BINARY_OUTCOMES` order — rather than §3.3's table order, which groups by
    KIND. Grouped by outcome is what §7.3's atomicity is about, and it is what makes the rendered
    counters table read down the estimand groups a failure moves together.

    `<outcome>.augmented` exists ONLY where the frozen path is not `"unaugmented"`, so `sich` and
    `ph2` are ABSENT rather than present-and-None (§3.3). The path is READ off `Secondary` and never
    recomputed, which is §6.3.
    """
    keys = ["beta", *(f"rd_{k}" for k in est.rd)]
    for key in C.BINARY_OUTCOMES:
        keys.extend((f"{key}.rd", f"{key}.odds_ratio"))
        if sec.estimates[key].augmented_path != _DECLARED_PATHS[-1]:
            keys.append(f"{key}.augmented")
    return tuple(keys)


def _tested(key: str) -> bool:
    """Does [§10] prescribe a p-value for this estimand? TRUE for 8 of §3.3's 26 keys.

    `beta` by §9.1 and each `<outcome>.rd` by §9.2. NOT the six `rd_<k>` ([§8] gives them intervals
    only, §9.3), NOT any `odds_ratio` ([§10] refuses to condition on the replicates where the ratio
    exists, §9.2) and NOT any `augmented` (one test per outcome; a second is a multiplicity [§13] is
    not told about, §9.2).

    A function and not a frozen set, so that an eighth outcome added to the [§5] registry gains its
    p-value in the same edit that gains its estimate -- Stage 9 §4.1's move, one stage on.
    """
    return key == "beta" or key.endswith(".rd")


def _bucket(message: str) -> str:
    """The `C.FAILURE_BUCKETS` entry for a raised `FitError`, by the FIRST TOKEN of its message.

    `model.FitError` carries no code -- it is `class FitError(RuntimeError)` with no attributes
    (model.py:89) -- and Stage 8 §11 requires the separation count reported SEPARATELY from the
    convergence count, so Stage 10 must classify and classification is textual (§7.2).

    **It RAISES on an unrecognised token rather than defaulting to a catch-all bucket**, and §15.6 is
    why: a default would make that section's scan cosmetic. The scan would pass, the map would be
    incomplete, and the counter Stage 8 §11 asked to be separate would be silently merged into
    whatever the default was.

    A `code` field on `FitError` set at each raise site would make this structural instead of
    textual, and §7.2 defers it rather than adopting it: it is a nineteen-site amendment to the
    module Stages 6, 8, 9 and 12 all depend on, on the strength of a diagnostic, and it would move
    what §15.6 asserts rather than removing the need for it.
    """
    token = message.split()[0] if message.split() else ""
    if token not in C.FAILURE_BUCKETS:
        raise C.SchemaError(
            f"a FitError leading with {token!r} is in no C.FAILURE_BUCKETS bucket. Known tokens: "
            f"{sorted(C.FAILURE_BUCKETS)}. Classification is by the first token of the message "
            "because FitError carries no code [Stage 10 §7.2], so a reworded or new raise site "
            "must be added to that map -- and this raises rather than defaulting, because a "
            "catch-all bucket would silently merge the separation count Stage 8 §11 asked to be "
            f"reported separately. The message was: {message.splitlines()[0][:200]}")
    return C.FAILURE_BUCKETS[token]


# --- the replicate body [§6] ---------------------------------------------------------------------------

def _shared_design(draw: pd.DataFrame, ps: propensity.Propensity) -> pd.DataFrame:
    """The ONE design the full-list [§8] outcomes share, and the check that they may share it.

    §6.4: *"the four full-list outcomes share an identical design matrix"* -- a claim with a
    CONDITION Stage 9 asserted without. `model.design` is called on `df.loc[in_estimate]` and
    `in_estimate` is `ps.in_model & df[key].notna()` PER OUTCOME, so four outcomes share a design
    only if they are present on the same rows. They are, and for a reason asserted elsewhere rather
    than by luck: all four derive from the ordinal primary outcome and roadmap invariant 6 is that a
    derived dichotomy carries exactly the missingness of its ordinal source.

    **A mismatch RAISES `C.SchemaError` and there is no fall-back.** An earlier draft of §6.4 had it
    detect the mismatch and quietly build four designs, which is the opposite of that section's own
    "fails loudly if it stops being true": the run would finish, every number would be right, and
    nobody would learn that a Stage 3 guarantee had stopped holding. Resampling copies whole rows
    and cannot break invariant 6 by itself, so a mismatch inside a replicate means a derived
    dichotomy has stopped carrying its ordinal source's missingness -- upstream, and a bug. That is
    `SchemaError` by §7.1's taxonomy and explicitly NOT `FitError`: making it droppable would let
    `N_BOOT` replicates absorb a broken invariant one at a time (§15.9).

    **The four are picked by TWO conditions and an earlier draft used only the first.** Asking
    `C.outcome_model_covariates` which outcomes take the shared [§6] list returns SIX on v7, not four:
    two of them decline no covariate list, they decline a MODEL — they are unaugmented, read directly
    from the workbook, and carry their own missingness rather than an ordinal source's. The second
    condition is the one §6.4's argument actually rests on: the outcome must be a DERIVED DICHOTOMY,
    which is what roadmap invariant 6 is a statement about. Measured: six against four, and the two
    extra are exactly the two unaugmented ones.

    The membership test is `C.DERIVED_DICHOTOMIES`, which `config.py` computes from the registry as
    the outcomes whose `source is not None`. That constant and not the source's NAME, because Stage 8
    §12 makes `outcome.py` the only shipped module that may name the [§5] primary outcome and this
    module is not it. The sources are then asserted to be ONE, because invariant 6 is a statement
    about an outcome and ITS OWN source: two derived dichotomies of different ordinal columns would
    have different missingness and no shared design.

    The accessor is asked rather than `C.OUTCOME_MODEL_OVERRIDES` being read -- §15.11 asserts by
    scan that this module names neither that constant nor `_augmentable` nor `_augmented_path`.
    `tici_2b_3` is the one outcome the accessor answers differently for, and it builds its own design
    inside `outcome.secondary`, which is right: sharing the full-list design would silently estimate
    the wrong nuisance model.

    **[SPEC] This is a per-replicate INVARIANT CHECK and §6.4's 27.2 s saving is not realised.** The
    design is built once here, and `outcome.secondary` then builds its own four inside §7.3's loop,
    because §7.3 forbids re-driving that loop in this module and `secondary` takes no design. The
    two sections cannot both be honoured without a parameter §13 does not license; §7.3's is the
    inference decision and wins. What survives is everything §15.9 asserts -- the masks agree, one
    design is built, and the invariant fails loudly -- at about 4.5 ms per replicate.
    """
    full_list = tuple(key for key in C.BINARY_OUTCOMES
                      if key in C.DERIVED_DICHOTOMIES
                      and C.outcome_model_covariates(key) == C.OUTCOME_COVARIATES)
    sources = sorted({str(C.OUTCOMES[key].source) for key in full_list})
    if not full_list or len(sources) > 1:
        raise C.SchemaError(
            f"the [§6] full-list derived dichotomies are {list(full_list)} over {len(sources)} "
            f"ordinal source(s) {sources}. A shared design needs exactly one source: roadmap "
            "invariant 6 says a derived dichotomy carries exactly the missingness of ITS OWN "
            "source, so two dichotomies of different columns have different masks and share "
            "nothing, and none at all means there is no shared design to build "
            "[Stage 10 §6.4].")
    masks = {key: ps.in_model & draw[key].notna() for key in full_list}
    sizes = {key: int(mask.sum()) for key, mask in masks.items()}
    reference = full_list[0]
    differing = sorted(key for key in full_list if not masks[key].equals(masks[reference]))
    if differing:
        raise C.SchemaError(
            f"the [§6] full-list outcomes do not share an estimation population: {', '.join(differing)} "
            f"differ(s) from {reference}. Sizes: "
            + ", ".join(f"{key} {sizes[key]}" for key in full_list)
            + ". All of them derive from the same ordinal source, and roadmap invariant 6 is that a "
            "derived dichotomy carries exactly the missingness of that source -- so their notna() "
            "masks are identical BY AN INVARIANT and not by coincidence. Resampling copies whole "
            "rows and cannot break it, so this is upstream and it is a bug: SchemaError and not "
            "FitError, because making it droppable would let N_BOOT replicates absorb a broken "
            "Stage 3 guarantee one at a time [Stage 10 §6.4, §15.9].")
    return model.design(draw.loc[masks[reference]], C.OUTCOME_COVARIATES)[0]


def _replicate(draw: pd.DataFrame, keys: tuple[str, ...], paths: dict[str, str],
               source: object) -> Replicate:
    """One replicate: refit Stages 6, 8 and 9 on `draw` and return every estimand it reached.

    [§10] refits *"the propensity model in every replicate"* and the roadmap's Stage 10 entry expands
    it to the outcome regression, the weights and every estimate, so the whole of Stages 6, 8 and 9
    runs inside this function and nothing is carried in from the point estimate except `paths`.

    The order is the POINT ESTIMATE's order and it is not an optimisation target: `propensity.fit`
    first because everything downstream reads `e`, `w` and `in_model`, `primary` before the binaries
    because that is the order the point estimate ran in and a log that interleaves them differently
    is a log that cannot be diffed against it (§6.1).

    **Three failure granularities, and they are §7.1's and §7.3's:**

      * a `model.FitError` from `propensity.fit` fails the WHOLE replicate and is counted against
        EVERY key -- not against none. Nothing downstream runs without `e`, `w` and `in_model`, so
        there is no estimand the failure does not affect, and a reading in which those replicates
        are simply never *attempted* is precisely the silent drop §15.4 exists to make impossible.
        DECISION 7's own arithmetic assumes this rule: it records that the primary keeps 1998 draws,
        which is `C.N_BOOT` less exactly the two propensity failures §7.5 measured.
      * a `FitError` from `outcome.primary` costs `beta` and the six `RD_k` TOGETHER. They come from
        one `polr` fit and a comparison across them is [§16]'s constant-shift statement, so they
        must share draws (§7.3).
      * a `FitError` inside one binary outcome costs that outcome's `rd`, `odds_ratio` AND
        `augmented` together, and costs no other outcome anything. **The group is ATOMIC and that
        includes `augmented`**: dropping `augmented` alone while `rd` survives is `tau` and `rd`
        taken over different replicate sets, which is exactly what Stage 9 §14 prohibits for two
        numbers [§10.3] reports as a comparison (§7.3, §15.10).

    **`C.SchemaError` is never caught here**, which is §7.1 as code. The only `except` is on
    `model.FitError`. A `SchemaError` means a frame that cannot be read, and for a frame this stage
    constructed that means this stage constructed it wrongly.

    The `Audit` is constructed HERE and discarded when this function returns, so the process's memory
    is not a function of `N_BOOT` (§6.2): nine entries per replicate at 2000 replicates is 18 000,
    each carrying a rendered table and a detail string of several hundred characters. It carries the
    POINT ESTIMATE's own `Source` rather than a synthetic one -- a second module-level `Source` is
    pinned against by `data.SOURCES` and a test, and the replicate really is drawn from the workbook,
    so the label in a log nobody writes should still be true.

    `keys` is a parameter and §11's fence does not have it. **[SPEC]** §7.1 requires the propensity
    failure counted against every key and this function cannot name them; `run` has both `est` and
    `sec` in hand before the loop, so `_estimand_keys` moves above it.
    """
    audit = Audit(source)
    try:
        ps = propensity.fit(draw, audit)                       # [Stage 6] — refit, never reused
    except model.FitError as failure:
        bucket = _bucket(str(failure))
        return Replicate(values={}, failures={key: bucket for key in keys},
                         n_alpha=None, polr_iterations=None, sum_w=None, n_in_model=None)

    values: dict[str, float] = {}
    failures: dict[str, str] = {}
    n_alpha: int | None = None
    polr_iterations: int | None = None
    sum_w = float(ps.w[ps.in_model].sum())                     # over the MASK, never notna() [§9]
    n_in_model = int(ps.in_model.sum())

    primary_keys = tuple(k for k in keys if k == "beta" or k.startswith("rd_"))
    try:
        est = outcome.primary(draw, ps, audit)                 # [Stage 8] — one group (§7.3)
        values["beta"] = est.beta
        values.update({f"rd_{k}": v for k, v in est.rd.items() if f"rd_{k}" in primary_keys})
        n_alpha = len(est.fit.alpha)
        polr_iterations = est.fit.iterations
    except model.FitError as failure:
        bucket = _bucket(str(failure))
        failures.update({key: bucket for key in primary_keys})

    _shared_design(draw, ps)                                   # invariant 6, per replicate (§6.4)

    # `collect=True` is §7.3's mechanism and `paths` is §6.3's frozen map. Neither is optional: with
    # `collect=False` one outcome's unfittable m_a(X) would raise past the other six, and with
    # `paths=None` `ph2` would re-decide its own augmentation path in 40.6% of replicates and the
    # interval would be a quantile over a near-even mixture of two estimators.
    sec = outcome.secondary(draw, ps, audit, paths=paths, collect=True)     # [Stage 9]
    max_abs_beta: dict[str, float] = {}
    or_corrected: dict[str, bool] = {}
    for key in C.BINARY_OUTCOMES:
        group = tuple(k for k in keys if k.startswith(f"{key}."))
        if key in sec.failures:
            failures.update({k: _bucket(sec.failures[key]) for k in group})
            continue
        estimate = sec.estimates[key]
        values[f"{key}.rd"] = float(estimate.rd)
        values[f"{key}.odds_ratio"] = float(estimate.odds_ratio)
        if f"{key}.augmented" in group:
            values[f"{key}.augmented"] = float(estimate.augmented)
        or_corrected[key] = bool(estimate.or_corrected)
        if estimate.fit is not None:
            max_abs_beta[key] = float(np.max(np.abs(estimate.fit.beta)))

    return Replicate(values=values, failures=failures, n_alpha=n_alpha,
                     polr_iterations=polr_iterations, sum_w=sum_w, n_in_model=n_in_model,
                     max_abs_beta=max_abs_beta, or_corrected=or_corrected)


# --- collecting the replicates [§7] ---------------------------------------------------------------------

def _collect(collected: tuple[object, ...], keys: tuple[str, ...]) -> dict[str, Draws]:
    """One `Draws` per estimand key, over the replicates that ATTEMPTED it (§7.1, §7.3).

    `n_attempted` is counted and not assumed: it is the number of replicates in which the key
    appeared in `values` or in `failures`, which on [§10]'s body is `C.N_BOOT` for all twenty-six —
    an ASSERTION (§15.4) rather than a definition, because `replicates` is general and Stage 12's
    bodies may attempt an estimand selectively.

    Every draw is checked FINITE before it enters the array. A `nan` odds ratio is unreachable on
    this pipeline — `marginal_odds_ratio` returns `nan` only when a weighted proportion is `nan`, and
    S7 raises on that first, measured 0 of 1998 — and it is checked anyway because a `nan` in a draws
    array is a percentile that comes back `nan` and an interval nobody can read as a failure (§8.2).
    """
    out: dict[str, Draws] = {}
    for key in keys:
        drawn: list[float] = []
        counts: dict[str, int] = {}
        attempted = 0
        for replicate in collected:
            if key in replicate.values:
                value = float(replicate.values[key])
                if not np.isfinite(value):
                    raise C.SchemaError(
                        f"{key}: a replicate produced the non-finite draw {value!r}. A percentile "
                        "over an array holding one comes back nan, which is an interval a reader "
                        "cannot tell from a failure [Stage 10 §8.2]. This is unreachable on the "
                        "[§10] body -- S7 raises before a weighted proportion can be nan -- so it "
                        "means an estimator returned something its own preconditions forbid.")
                drawn.append(value)
                attempted += 1
            elif key in replicate.failures:
                bucket = replicate.failures[key]
                counts[bucket] = counts.get(bucket, 0) + 1
                attempted += 1
        out[key] = Draws(quantity=key, draws=np.asarray(drawn, dtype=float),
                         n_attempted=attempted, failures=counts)
    return out


def _diagnostics(collected: tuple[object, ...]) -> Diagnostics:
    """`Replicate`'s eight fields aggregated into `Diagnostics`' six (§15.15).

    **A `None` is DROPPED and never counted as a zero**, which is the one way each of these four
    aggregations can be wrong that nothing else here would catch: the replicates §7.1 kills before
    `outcome.primary` ran contribute `None` to `n_alpha`, `polr_iterations`, `sum_w` and
    `n_in_model`, and counting those as zeros would put a spurious mode at zero in three
    distributions at once. So the four sum to the LIVE replicate count and not to `C.N_BOOT`.

    `max_abs_beta` holds exactly the augmented outcomes: `sich` and `ph2` are ABSENT from the dict
    rather than present with an empty array, because an empty array reads as "we looked and found
    none" where absence reads as "there is nothing here to look at", and the second is true (§10.2).

    `or_corrected` is keyed by ALL SEVEN outcomes and initialised to zero, because a rate of zero is
    a fact and a missing key is not one — Stage 9 §6.3 requires the correction named wherever the
    estimate appears, and an interval is somewhere it appears (§9.2). It is a COUNT and not a rate,
    so Stage 14 forms whichever denominator [§16] asks for; `_diagnostics_table` renders it over the
    surviving odds-ratio draws, which is §9.2's own rationale — the correction can only fire in a
    replicate that produced an odds ratio — and is what reproduces §9.2's own quoted percentages.
    """
    def tally(attribute: str) -> dict[int, int]:
        counts: dict[int, int] = {}
        for replicate in collected:
            value = getattr(replicate, attribute)
            if value is not None:
                counts[int(value)] = counts.get(int(value), 0) + 1
        return dict(sorted(counts.items()))

    sum_w = np.asarray([r.sum_w for r in collected if r.sum_w is not None], dtype=float)
    max_abs_beta: dict[str, list[float]] = {}
    or_corrected: dict[str, int] = {key: 0 for key in C.BINARY_OUTCOMES}
    for replicate in collected:
        for key, value in replicate.max_abs_beta.items():
            max_abs_beta.setdefault(key, []).append(float(value))
        for key, fired in replicate.or_corrected.items():
            or_corrected[key] += int(bool(fired))
    return Diagnostics(
        n_alpha=tally("n_alpha"),
        polr_iterations=tally("polr_iterations"),
        sum_w=sum_w,
        max_abs_beta={key: np.asarray(max_abs_beta[key], dtype=float)
                      for key in C.BINARY_OUTCOMES if key in max_abs_beta},
        n_in_model=tally("n_in_model"),
        or_corrected=or_corrected)


# --- the one audit entry [§10.1] ------------------------------------------------------------------------
#
# `data.py`'s KINDS stays nine, as it has for four stages. The counters and the distributions go in a
# single `model` entry, `bootstrap_replicates`, taking the ledger 30 -> 31:
#
#   load 7 / derive 4 / classify 1 / build 6 / fit 4 / assess 2 / primary 3 / secondary 3
#       / bootstrap 1  =  31
#
# ONE entry and not three, unlike Stages 8 and 9. All of it is one fact -- what the replicates looked
# like -- and Stage 9 §10.1's precedent for splitting was that its three entries described three
# different populations. `AuditEntry` holds one table, so the two blocks are concatenated into one
# grid with the second block's header kept as a labelled row and short rows padded, which is exactly
# what `propensity._padded` already does for Stage 6's two-table entry. That function is
# `propensity.py`'s private and is NOT imported: §0.2 says `propensity.py` is untouched, and reaching
# into it for a helper would be a module-boundary violation to save six lines.
#
# NO CELL OF EITHER BLOCK CONTAINS A `|`, which is why §10.2's column is `max_abs_beta` and not
# `max|beta|`. `data._md_table` (data.py:242-245) does no escaping and sizes its separator from
# `len(rows[0])`, so a pipe in a cell gives a body row with more markdown cells than the separator has
# dashes -- and byte-identity across hash seeds does not catch it, because a table identically broken
# under both seeds is still identical (§10.3, §15.13).


def _padded(table: tuple[tuple[str, ...], ...], width: int) -> tuple[tuple[str, ...], ...]:
    """Widen every row to `width` with the not-applicable dash the other stages' tables use.

    A six-line copy of `propensity._padded` rather than an import of it, and §0.2 is why: this stage's
    central claim is that `propensity.py` is untouched, and reaching across a module boundary for a
    private helper is the kind of coupling that makes "untouched" stop being checkable.
    """
    return tuple(row + ("—",) * (width - len(row)) for row in table)


def _counters_table(draws: dict[str, Draws]) -> tuple[tuple[str, ...], ...]:
    """One row per estimand key: what was attempted, what survived, and why the rest did not.

    The bucket columns are `C.FAILURE_BUCKETS`' VALUES, sorted, so every declared bucket has a column
    whether or not anything landed in it -- which is Stage 6 §7.5's rule ("every declared thing is
    rendered whether or not the data fills it") and Stage 8 §11's requirement in one move: a
    separation count that vanishes from the table when it reads zero is a count nobody can tell was
    checked, and Stage 8 measured that on this estimator it reads zero on data that is degenerate
    throughout.

    `draws` is a column and is NOT derived by subtracting the buckets from `attempted`: a table that
    reconciles with itself by construction has stopped being evidence [Stage 5 §7.3]. `Draws`'
    `__post_init__` is where the reconciliation is checked.
    """
    buckets = tuple(sorted(set(C.FAILURE_BUCKETS.values())))
    header = ("estimand", "attempted", "draws", "tested", *buckets)
    rows = tuple(
        (key, str(d.n_attempted), str(len(d.draws)), "yes" if _tested(key) else "no",
         *(str(d.failures.get(bucket, 0)) for bucket in buckets))
        for key, d in draws.items())
    return (header, *rows)


def _spread(values: np.ndarray) -> str:
    """`min / median / max` in one cell, through `_fmt`. No pipe, by §10.3.

    A range and not a full quantile row, because the point of these distributions is Stage 8 §11's
    and Stage 9 §9.6's question -- did anything move, and how far -- and a reader who needs the whole
    distribution has `Bootstrap.diagnostics` as a typed field.
    """
    if not values.size:
        return "none"
    return (f"min {_fmt(float(values.min()))} / median {_fmt(float(np.median(values)))} "
            f"/ max {_fmt(float(values.max()))}")


def _diagnostics_table(diagnostics: Diagnostics,
                       draws: dict[str, Draws]) -> tuple[tuple[str, ...], ...]:
    """The four distributions, the two safeguard-adjacent counters, and their DENOMINATORS.

    **The denominator is a column and not a footnote**, and §16 item 9 is why it has to be: the four
    scalar distributions sum to the LIVE replicate count while `or_corrected`'s denominator is that
    outcome's SURVIVING odds-ratio draws, and nothing in a grid of counts says which is which.
    §15.15 asserts that the four sum to the live count; this column is what tells a reader, and it is
    the item §16 files as having already fired once with nobody able to read the grid.

    Row order is declared throughout -- `sorted` over the integer-keyed tallies and
    `C.BINARY_OUTCOMES` over the outcomes -- so the rendered log is byte-identical across hash seeds
    (§15.13). No set iteration anywhere.
    """
    # **THERE ARE TWO LIVE COUNTS AND NOT ONE**, and §3.1 is what makes them differ: the five
    # diagnostic fields are `None` "when the replicate died before `outcome.primary` ran -- which is
    # exactly the propensity-fit failure §7.1 describes". So `sum_w` and `n_in_model` survive a
    # PRIMARY-fit failure and `n_alpha` and `polr_iterations` do not, because those two are
    # properties of a fit that did not happen. §15.15 says all four "sum to the number of replicates
    # that reached `outcome.primary`" and the two counts coincide only while the primary never
    # fails -- which is the workbook's case, measured at G7 = 0 (§7.5), and is exactly the condition
    # under which one shared denominator is silently wrong.
    fitted_primary = int(sum(diagnostics.n_alpha.values()))
    fitted_propensity = int(sum(diagnostics.n_in_model.values()))
    header = ("diagnostic", "value", "replicates", "denominator")
    rows: list[tuple[str, ...]] = []
    for label, tally, live in (("n_alpha", diagnostics.n_alpha, fitted_primary),
                               ("polr_iterations", diagnostics.polr_iterations, fitted_primary),
                               ("n_in_model", diagnostics.n_in_model, fitted_propensity)):
        if not tally:
            rows.append((label, "none", "0", f"{live} replicate(s) reached it"))
        for value, count in tally.items():
            rows.append((label, str(value), str(count), f"{live} replicate(s) reached it"))
    rows.append(("sum_w", _spread(diagnostics.sum_w), str(int(diagnostics.sum_w.size)),
                 f"{fitted_propensity} replicate(s) reached it"))
    for key in C.BINARY_OUTCOMES:
        if key in diagnostics.max_abs_beta:
            values = diagnostics.max_abs_beta[key]
            rows.append((f"max_abs_beta {key}", _spread(values), str(int(values.size)),
                         "m_a(X) fits — AUGMENTED outcomes only"))
    for key in C.BINARY_OUTCOMES:
        # **THE DENOMINATOR IS THE SURVIVING DRAW COUNT AND NOT `n_attempted`**, and §9.2's own
        # rationale is why: "the correction can only fire in a replicate that produced an odds
        # ratio". §9.2 names `Draws.n_attempted` and then quotes rates -- 16.1%, 13.4%, 2.5% -- whose
        # denominators are 1982, 1995 and 1998, which are surviving counts and not `n_attempted`,
        # since §3.1 holds `n_attempted` at `C.N_BOOT` for every key. The two halves of §9.2 disagree;
        # the rationale is the operative one, and it is the one that reproduces its own numbers.
        # `Diagnostics.or_corrected` is the COUNT either way, so Stage 14 can form whichever rate
        # [§16] asks for -- this cell is what a reader of the log sees.
        surviving = len(draws[f"{key}.odds_ratio"].draws)
        rows.append((f"or_corrected {key}", "fired", str(diagnostics.or_corrected[key]),
                     f"{surviving} with an odds ratio"))
    return (header, *rows)


def _replicates_detail(draws: dict[str, Draws], diagnostics: Diagnostics) -> str:
    """The caption both blocks share: what was resampled, what was frozen, and what is NOT here.

    It names the seed and the replicate count because [§16] requires both reported, and it names the
    two things a reader of a counters table would otherwise have to infer: that a `propensity.fit`
    failure is counted against every estimand rather than vanishing (§7.1), and that each estimand
    group keeps its own replicate set (§7.3).
    """
    live = int(sum(diagnostics.n_alpha.values()))
    tested = sum(1 for key in draws if _tested(key))
    return (
        f"[§10] patient-level nonparametric bootstrap, {C.N_BOOT} replicate(s), stratified by "
        f"{C.BOOT_STRATUM}, seed {C.SEED}. The stratum totals are FIXED across replicates, so no "
        "replicate can lose a centre and the inference is conditional on the participating centres "
        "and their observed treatment practices [§10]. Every replicate refits the [§7] propensity "
        "model, the [§8] primary estimator and the seven [§8] binary estimators; nothing is carried "
        "in from the point estimate except the augmentation paths, which are DECIDED ONCE on the "
        "point estimate and frozen, because a rule re-evaluated per replicate would make one "
        f"interval a quantile over a mixture of two estimators [Stage 9 §9.5]. {len(draws)} "
        f"estimand(s), of which {tested} carry a p-value: [§10] prescribes one for the [§8] "
        "treatment coefficient and one per binary outcome on the RISK-DIFFERENCE scale, and none "
        "for a cumulative RD_k [§8], none for a marginal odds ratio -- which becomes undefined when "
        "a weighted proportion reaches 0 or 1, so taking p from those draws would condition on the "
        "replicates where the ratio happens to exist -- and none for an augmented estimate, since "
        f"[§10] prescribes one test per outcome. {live} replicate(s) reached the primary fit; a "
        "replicate whose propensity model fails is counted as a failure against EVERY estimand "
        "rather than disappearing from every denominator uncounted, and each estimand group keeps "
        "its own surviving set, so one outcome's failure never thins another's draws [§10]. "
        "Replicates whose prespecified fit fails are dropped and counted, NEVER substituted with a "
        "different estimator [§10, invariant 5]. The percentile limits and the p-values are not in "
        f"this table: they are {C.CI_LEVEL:g} percentile intervals under the "
        f"{C.PERCENTILE_METHOD} definition, reported by Stage 14 with the surviving-draw count "
        "beside each. Estimation, not testing, is the reportable output [§10].")


def _record_replicates(draws: dict[str, Draws], diagnostics: Diagnostics, audit: Audit) -> None:
    """The one `model` entry, and `case_ids` is EMPTY on it.

    Nine `model` entries name cases and this one removes no patient; naming the `N_BOOT` x n drawn
    rows would be a log the size of the data. The excluded-record accounting is `propensity.fit`'s,
    per replicate, in a throwaway nobody reads -- which is §6.2's point (§10.1).

    `n` is the replicate count, which is what the entry is about. The two blocks are concatenated
    with the second header kept as a labelled row and the short rows padded to the first block's
    width; `data._md_table` sizes its columns from `rows[0]` and raises `IndexError` on a narrower
    row, which is the defect Stage 6 §7.2 shipped first time (§10.1, §15.13).
    """
    counters = _counters_table(draws)
    distributions = _diagnostics_table(diagnostics, draws)
    audit.record("model", _STEP_REPLICATES, C.N_BOOT,
                 _replicates_detail(draws, diagnostics),
                 table=counters + _padded(distributions, len(counters[0])))


# --- the preconditions [§4.4] ---------------------------------------------------------------------------

def _assert_run_inputs(df: pd.DataFrame, ps: propensity.Propensity, est: outcome.Primary,
                       sec: outcome.Secondary) -> None:
    """R1-R8, in two phases, collected. Every one is `C.SchemaError` and this stage adds no FitError.

    Stage 9 §4.4a's boundary is the rule rather than a paraphrase: *"the boundary is
    can-this-be-READ vs is-the-DATA-judgeable -- NOT mask-checks vs the-rest."* A column-presence
    check is a precondition of every read that follows it, so R1 and R2 are phase 1 and R6, which
    reads the stratum column, is phase 2.

    **All eight are `SchemaError` and that asymmetry is the opposite of the one Stage 9 needed.**
    Every condition here is a property of the CALL and not of a replicate, so all eight are the
    caller's bug. The only `FitError` this stage ever sees is one raised by an estimator inside a
    replicate, and §7.1 is what `_replicate` does with it.

    A stratum of size 1 satisfies R6 and contributes no variability -- it resamples to itself in
    every replicate. That is not an error and is not guarded: it is what stratifying on a
    near-determining variable means, and [§10] chose it knowing so (§4.4).
    """
    # PHASE 1 -- can this frame, this Propensity and this Secondary be READ?
    bad: list[str] = []
    absent = [c for c in (C.BOOT_STRATUM,) if c not in df.columns]
    if absent:
        bad.append(
            f"R1  {C.BOOT_STRATUM}: not a column of the frame. [§10] resamples at patient level "
            "WITHIN it, so without it there is no stratification and the inference stops being "
            "conditional on the participating centres [§10, §4.1].")
    if "case_id" not in df.columns:
        bad.append(
            "R2  case_id: not a column of the frame. `resample` REWRITES it -- the k-th draw of a "
            "patient becomes `case_id#k`, which is what stops `propensity._record_exclusion` "
            "reading a replicate's duplicated rows as a defect [§5.2, §5.3].")
    misaligned = [name for name, s in (("e", ps.e), ("w", ps.w), ("in_model", ps.in_model))
                  if not s.index.equals(df.index)]
    if misaligned:
        bad.append(
            f"R3  the Propensity is not aligned to this frame: {', '.join(misaligned)} "
            f"carr{'ies' if len(misaligned) == 1 else 'y'} a different index. `.index.equals` and "
            "NOT a length check or a set comparison: `.loc[boolean_series]` returns rows in the "
            "SERIES' order, so a permuted-but-equal index pairs each record's weight with another "
            "record's fitted value and returns a different number with no raise [Stage 9 §12]. "
            "This stage constructs a Propensity per replicate, so it is the caller Stage 9's S3a "
            "was written against.")
    if ps.in_model.dtype != bool or ps.in_model.isna().any():
        bad.append(
            f"R4  in_model is {ps.in_model.dtype} and carries "
            f"{int(ps.in_model.isna().sum())} missing value(s). It is boolean and TOTAL by "
            "construction [Stage 6 §4.4], and pandas reads a non-boolean Series as LABELS -- so a "
            "masked read raises KeyError rather than reporting a mask problem [Stage 9 §4.4a].")
    if set(sec.estimates) != set(C.BINARY_OUTCOMES):
        bad.append(
            f"R5  the Secondary carries estimates for {sorted(sec.estimates)} against "
            f"{sorted(C.BINARY_OUTCOMES)}. §6.3's frozen path map is built from these keys and must "
            "be COMPLETE: a partial map would let some outcomes carry the point estimate's path and "
            "others re-decide from the replicate, which is the two-estimator mixture [§10] forbids "
            "arrived at one outcome at a time. A `Secondary` short of a key is also a point "
            "estimate that was computed with `collect=True`, which is a replicate's contract and "
            "not the workbook's [§7.3].")
    if bad:
        raise C.SchemaError(
            "\n  " + "\n  ".join(bad)
            + f"\n  {len(bad)} bootstrap assertion(s) failed against {len(df)} records. The frame, "
            "the Propensity or the Secondary cannot be read, so the phase 2 checks were not run: "
            "they read through the columns and the object this phase is about (§4.4).")

    # PHASE 2 -- readable. Is the DATA resamplable, and does the configuration support an interval?
    bad = []
    strata = df[C.BOOT_STRATUM]
    if strata.isna().any():
        bad.append(
            f"R6  {C.BOOT_STRATUM}: {int(strata.isna().sum())} record(s) carry no stratum label. "
            "`groupby` drops them silently, so the replicate would be SHORTER than the frame and "
            "every weighted denominator would be computed over a population nobody chose [§4.4].")
    empty = [str(label) for label, size in strata.value_counts(dropna=True).items() if size == 0]
    if empty:
        bad.append(
            f"R6  {C.BOOT_STRATUM}: {', '.join(empty)} carr{'ies' if len(empty) == 1 else 'y'} no "
            "record. A stratum of size 0 cannot be resampled to its own size [§4.1].")
    off_path = sorted(
        f"{key}={estimate.augmented_path!r}" for key, estimate in sec.estimates.items()
        if estimate.augmented_path not in _DECLARED_PATHS)
    if off_path:
        bad.append(
            f"R7  augmented_path is not one of {list(_DECLARED_PATHS)} for: {', '.join(off_path)}. "
            "The path is READ off the point estimate and frozen into every replicate (§6.3), so an "
            "undeclared one is a value every replicate would carry and no estimator would "
            "recognise.")
    if C.N_BOOT < C.ci_min_draws(C.CI_LEVEL):
        bad.append(
            f"R8  N_BOOT is {C.N_BOOT} against a floor of {C.ci_min_draws(C.CI_LEVEL)} at level "
            f"{C.CI_LEVEL:g}. Below the floor the lower percentile limit is the sample MINIMUM and "
            "np.percentile returns it while still calling it a percentile, so no replicate could "
            "produce an interval whose lower limit is an order statistic at all [§8.3]. "
            "test_config.py carries this as a static assertion; this is its runtime restatement.")
    if bad:
        raise C.SchemaError(
            "\n  " + "\n  ".join(bad)
            + f"\n  {len(bad)} bootstrap assertion(s) failed against {len(df)} records in "
            f"{strata.nunique(dropna=True)} stratum/strata.")


# --- the stage [§10, §11] -------------------------------------------------------------------------------

def run(df: pd.DataFrame, ps: propensity.Propensity, est: outcome.Primary,
        sec: outcome.Secondary, audit: Audit) -> Bootstrap:
    """The [§10] bootstrap: N_BOOT stratified replicates, percentile limits, prespecified p-values.

    Takes the point estimate's `Primary` and `Secondary` and does not recompute them. It reads
    exactly two things from them -- `sec.estimates[k].augmented_path` for §6.3's frozen map, and the
    estimand keys, whose `rd_<k>` half comes from `est.rd` -- and no estimate, because a percentile
    limit is an order statistic of the draws and is not a function of the point value (§8.1, §15.8).
    `est` is a parameter BECAUSE of that second read: an earlier draft passed only `sec` to
    `_estimand_keys`, which left `est` unused and the six `RD_k` keys sourced from nowhere (§3.3).

    Takes no `n`, no `seed` and no `level`: C.N_BOOT, C.SEED and C.CI_LEVEL are the prespecified
    ones and [§10] refits this in every replicate, so they are not runtime knobs. This is Stage 6
    §6.4's choice, made for its reason -- `propensity.fit` takes no covariate list because the
    specification is not an option -- and `replicates` above IS parameterised, because that one is
    general and this one is [§10].

    Six things about the order below, each of which is a failure if moved:

      * `_assert_run_inputs` runs FIRST, so a caller passing a `Secondary` over a different frame
        reports R3/R5 rather than producing 2000 replicates of a misalignment (§4.4).
      * `paths` is built ONCE, before the loop, from `sec`. Building it inside would make it a
        function of the replicate, which is §6.3's whole prohibition.
      * THE LOOP IS `replicates`' AND NOT THIS FUNCTION'S. One Generator, created once inside
        `replicates`, consumed in order; not one per replicate (§4.3). `run` writing its own loop
        would put the seeding and the ordering in two places, and Stages 12 and 13 would inherit
        the copy [§10] never ran (§3.2).
      * The per-replicate `Audit` is constructed INSIDE `_replicate` and discarded when it returns,
        so the process's memory is not a function of N_BOOT (§6.2).
      * `_record_replicates` runs AFTER the loop and before the intervals are built, so a log exists
        naming the counters even if a percentile call then raises on an estimand nobody expected to
        be empty.
      * An estimand below §8.3's floor keeps its `Draws` and gets NO `Interval`. `draws` therefore
        always holds every key and `intervals` may hold fewer (§3.3) -- asserted, both halves,
        because an implementation emitting a pair of extremes under a percentile's name passes any
        test that only checks the `Draws` side (§15.8).

    NO SchemaError IS CAUGHT ANYWHERE IN THIS FUNCTION, and that is §7.1 as code. The only `except`
    is on `model.FitError`, inside `_replicate`, per estimand group.

    Deterministic given C.SEED, and writes no file. Appends ONE `model` entry (§10.1). The only
    mutable object touched is the `Audit` passed in -- the per-replicate ones are its own.

    **[SPEC]** `_estimand_keys` is called BEFORE `replicates` and not after, because §7.1 requires a
    `propensity.fit` failure counted against every key and `_replicate` cannot name them otherwise.
    """
    _assert_run_inputs(df, ps, est, sec)                         # R1-R8
    paths = {key: e.augmented_path for key, e in sec.estimates.items()}
    keys = _estimand_keys(est, sec)                              # §3.3 -- twenty-six of them

    collected = replicates(                                      # §3.2 -- ONE loop, one home
        df,
        lambda draw: _replicate(draw, keys, paths, audit.source),
        C.N_BOOT, C.SEED, C.BOOT_STRATUM,
    )

    draws = _collect(collected, keys)                            # §7.1, §7.3, per group
    diagnostics = _diagnostics(collected)                        # §7.4, §7.5, §9.2, §10.2
    _record_replicates(draws, diagnostics, audit)                # §10.1

    intervals: dict[str, Interval] = {}
    for key, d in draws.items():
        if len(d.draws) < C.ci_min_draws():                      # §8.3 -- no Interval, Draws kept
            continue
        lo, hi = percentile_ci(d.draws, C.CI_LEVEL)
        p = bootstrap_p(d.draws) if _tested(key) else None       # §9.1, §9.3 -- 8 of 26
        intervals[key] = Interval(lo, hi, C.CI_LEVEL, C.PERCENTILE_METHOD, len(d.draws), p)

    return Bootstrap(C.SEED, C.N_BOOT, draws, intervals, diagnostics)
