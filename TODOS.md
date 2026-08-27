# TODOS

Deferred work, with enough context that someone picking it up in three months knows why it exists.
Each item names what blocks it and what it unblocks. Nothing here is a reminder to be tidy; if an
item has no consequence, it does not belong.

---

## Declare an explicit dtype for the CTP volume columns in Stage 1's contract

**Surfaced by:** the engineering review of `specs/stage3_derived_variables.md`, 2026-08-10 (§20 R5
and §13 there).

**What.** Give `core_ml`, `tmax6_ml` and `penumbra_ml` a declared `dtype` in
`config.COLUMN_CONTRACT`, so they enter `READ_DTYPES` like every other numeric analysis column, and
decide between `float64` and the nullable `Float64` as part of it.

**Why.** All three currently pass `dtype=None`, so pandas infers, and the inferred dtype depends on
the file rather than on the contract:

```
  core_ml from data/…v7….xlsx            float64   (one record has no core volume)
  core_ml from tests/fixture_schema.xlsx  int64     (both its values are integers)
```

That asymmetry is not cosmetic — it decides which Stage 3 masks do work. Per Stage 3 §3.1, `float64`
**absorbs** missing in a comparison (`NaN > 5.0` is `False`), so `_core_above_median`'s
`.mask(core.isna())` is the one mask in Stage 3 that is load-bearing today; the identical-looking call
in `_dichotomies` is redundant because `mrs_90d` is declared `Int64` and propagates. A reader who does
not know which column is which will draw the wrong conclusion about both. Declaring the dtype makes
Stage 3 §3.1's table a property of the contract rather than of whichever workbook was loaded.

**Pros.** One dtype per analysis column regardless of source. Removes the class of surprise where a
test frame and the workbook exercise different branches of the same line. Makes a future stage that
*does* depend on the dtype (Stage 6's design matrix, Stage 10's resampling) safe to write against the
contract.

**Cons.** It is a Stage 1 amendment, so it touches `READ_DTYPES`, Stage 2's A6 plausible-range checks,
and Stage 1's own §9.9 read-contract tests. `Float64` versus `float64` is a real decision and not a
detail: nullable `Float64` would make `core > median` propagate `<NA>` instead of absorbing, which
changes `_core_above_median`'s mask from load-bearing to redundant and makes Stage 3 §6.3's
justification stale. Whichever is chosen, Stage 3 §3.1, §6.3, §13 and §18's dtype rows must be updated
in the same commit.

**Depends on / blocked by.** Nothing technically. Deliberately deferred until **after Stage 3 lands**:
Stage 3 §18 records the current dtypes as verified facts and §12.10/§12.12 are written against them, so
changing the contract mid-implementation would invalidate the spec the implementation is being checked
against. Do it as its own commit, with the Stage 1, Stage 2 and Stage 3 spec edits included.

**Status: unblocked.** Stage 3 landed on 2026-08-10. The tests this will move are now concrete rather
than anticipated — `test_derive.py`'s `test_without_the_mask_the_missing_core_joins_the_below_median_group`
asserts `core_ml` is `float64` as its stated premise, and `test_the_fixture_reaches_every_entry` asserts
the fixture's is `int64`. Both are deliberate records of today's asymmetry, and both must be revisited by
whoever declares the dtype — the first one's whole point disappears if `Float64` is chosen.

---

## Correct Stage 4 §4.2 and Definition of done 9: a fill on E5 is a no-op, not a bug

**Surfaced by:** implementing `specs/stage4_eligibility_classification.md`, 2026-08-12, while working
Definition of done 9 — the one item of the thirteen that could not be satisfied.

**What.** Two statements in the Stage 4 spec are false against pandas 2.3.3, and they contradict the
spec's own §3.1 and §18:

- **§4.2** — "`.fillna(False)` on E5 would convert 'this patient's arm was never recorded' from
  *caught* into *ignored*", and E4/E5 are "the one place in this module where a fill would be a bug
  rather than belt-and-braces".
- **Definition of done 9** — "A `.fillna(False)` pasted onto E5's mask has been **seen to make a test
  fail** — 12.3's missing-exposure half."

Both rest on E5's mask being able to carry `<NA>`. It cannot. Measured:

```
  pd.Series([0, 1, 2, pd.NA], dtype="Int64")
    ~s.isin((0, 1))            → [False, False, True, True]   dtype boolean, NO <NA>
    ~s.isin((0, 1)).fillna(…)  → identical
```

So the fill is a **pure no-op**: pasting it onto E5 changes no behaviour and makes no test fail. This
was confirmed by running DoD 9's cycle — six of the seven break/restore cycles behaved as the spec
says; this one could not, because there is nothing to break.

**Why it happened.** §3.1's *third* fact — `isin` returns `False` for `<NA>` and never `<NA>` — was
added by the §20 engineering review of 2026-08-12 ("The one substantive addition to §3.1 is a third
fact"). §4.2's sentence and DoD 9 predate it and were not reconciled with it. §18 already records the
correct behaviour, so the document disagrees with itself rather than with the world.

**What the implementation did.** Nothing to the code: E4 and E5 carry no fill, exactly as specified.
The *reason* is now stated correctly in `eligibility.py`'s comment above them — a fill would assert
that the mask can carry `<NA>` when it cannot, which is what makes a reader stop believing the three
fills in `_classify` that are load-bearing; and it would become a real bug the moment either check is
rewritten as a comparison chain. `test_eligibility.py`'s
`test_the_isin_masks_carry_no_na_on_either_dtype_they_are_used_over` pins the property both checks
rest on, over both dtypes, so a pandas bump that made `isin` propagate fails loudly.

**Pros of doing it.** The spec is the sole source for Stage 4 and currently ships a false mechanism in
the section a reader consults before touching the fill discipline. Stage 5 and Stage 12 both add masks
over these columns and will read §4.2 for precedent.

**Cons.** It edits a spec that has already been reviewed and dated, and DoD 9 would be struck rather
than repaired — there is no cheap substitute cycle that demonstrates the same thing, because the code
is robust either way.

**Depends on / blocked by.** Nothing technical. It is the PI's document, so the amendment wants their
sign-off rather than a unilateral edit: **§4.2's mechanism sentence and DoD 9 both go, replaced by a
pointer to §3.1's third fact and to the test named above.** Nothing else in the spec moves — §3.1,
§18, §12.3 and the shipped code are all already correct.

---

## Commit a `tests/fixture_cohort.xlsx` that reads into a both-armed cohort

**Surfaced by:** the engineering review of `specs/stage5_cohort_construction.md`, 2026-08-13 (§12.0
and §13 there).

**What.** A second committed fixture workbook whose read produces a cohort with both arms at at least
one centre, plus the **generator script** that builds it — `tests/fixture_schema.xlsx` has none, which
is why extending it was declined. It gets its own `data.Source` declaration and its own spec section
when it arrives.

**Why.** Neither existing test input can reach the end of Stage 5, and this is measured rather than
anticipated:

```
                        records   after restriction 1   per retained centre (treated, control)
  fixture_schema.xlsx      2        keeps HUG only      HUG (1, 0)   → Stage 5's P2 RAISES
  hand_frame()             6        keeps HUG, CHUV     HUG (2, 0)   → P2 RAISES
                                                        CHUV (1, 0)
```

Both were built for stages that never removed a row, so both collapse to a cohort with no control arm
and Stage 5's both-arms postcondition stops them. That is the guard working, not a defect. Stage 5
works around it with `test_cohort.cohort_frame()` — a nine-record **frame** built in Python — which is
enough for every test Stage 5 writes. It is not enough for a test that needs a **file**: Stages 6-13
fit models, resample, and will want the whole pipeline runnable from a checkout with no `data/`
directory, which means from a `data.Source`.

**Pros.** Makes the full pipeline reproducible end to end with no patient data present — useful for
CI, for a fresh machine, and for anyone reproducing the analysis. A generator script makes the fixture
reviewable in a diff, which the current binary is not. Removes the need for each later stage to
rediscover why the schema fixture stops at Stage 5.

**Cons.** A new binary artifact plus a new `data.Source`, and `test_data.py` pins `data.SOURCES`
against a literal of exactly two (`test_sources_are_exactly_the_workbook_and_the_fixture`), so it is a
Stage 1/Stage 2 amendment as well as a new file. No test written before Stage 6 needs it.

**Depends on / blocked by.** Nothing technical. Build it when a stage first needs a `data.Source`
rather than a DataFrame — do **not** build it speculatively with Stage 5, whose acceptance tests are
all frame-driven. When it is built, give it a control at a treating centre (the property `cohort_frame`
adds as `COHORT-1`) or it will hit P2 exactly as the schema fixture does.

**Status update, 2026-08-14.** **Stage 6 is not the trigger**, checked rather than assumed during its
engineering review. Every Stage 6 acceptance test runs on a frame or on synthetic arrays, and
`test_propensity.py` takes the whole pipeline end to end with no `data/` through `cohort_frame()` — see
`specs/stage6_propensity_and_weights.md` §12.0. **Stage 10 is the likelier trigger**: a resampling test
wants a reproducible file rather than a frame built in Python. Re-evaluate there.

---

## Extract a shared table row-builder from `data.absence_by_column`

**Surfaced by:** the engineering review of `specs/stage6_propensity_and_weights.md`, 2026-08-14
(finding 7; the decision it produced is recorded in that spec's §7.4).

**What.** Split `data.absence_by_column` (`data.py:677`) into two pieces: a public row-builder returning
`(header, rows)` for a per-column, per-centre absence table, and a thin wrapper that records the
`missingness` audit entry from it. Any stage wanting that table under a *different* audit kind then calls
the builder instead of rewriting the loop.

**Why.** The function's own docstring is the argument for it existing at all:

> "Public because Stage 3 calls it too … A second implementation over there would drift on the branch
> that matters — `structural` versus `missing` — which is the exact misreading §10.2 exists to prevent.
> So there is one classifier, called twice with different column lists."

It is now called three times — `data.py`'s `_missingness`, `derive.py:381`, and `cohort.py:479` over the
built cohort — and every call renders under `missingness`, because recording that kind is baked into the
function. Stage 6's spec originally asked for a fourth rendering inside its new `model` kind, which would
have been a hand-rolled copy of the per-centre loop
(`(absent & (df["center"] == centre)).sum()` over `CENTER_ORDER`). Stage 6 solved that by **narrowing its
own table** to `covariate`, `n_absent`, `excluded_by` and pointing at `absence_by_cohort_column` for the
per-centre breakdown (§7.4). That is right for Stage 6 and does not generalise: Stage 7's within-centre
overlap table [§9] and Stage 14's flow diagram [§16] both want per-centre counts under headings that are
not `missingness`.

**Pros.** One implementation of the four-way `structural` / `not recorded` / `missing` / `complete`
classification and of the per-centre columns, for every stage rather than for every stage that records
`missingness`. Removes the reason a future stage has to re-derive Stage 6 §7.4's decision. The centre
columns keep coming from `CENTER_ORDER` rather than from the data, in one place, which is the property
that comment block exists to protect.

**Cons.** It reopens a landed, tested module that three modules import, in a repository whose acceptance
criterion is a byte-identical audit log — so `test_data.py`'s rendering tests are in the blast radius and
the split must be provably output-preserving. **No stage needs it today**: Stage 6 §7.4's narrowing means
Stage 6 does not, and Stages 7 and 14 do not exist yet.

**Depends on / blocked by.** Nothing technical. Do it when Stage 7 or Stage 14 first wants the table under
a non-`missingness` heading — that is the point where the duplication becomes real rather than
hypothetical, and where the refactor can be verified against a second real caller instead of an imagined
one. Whoever does it should read Stage 6 §7.4 first: it records what was declined and why, so the
refactor can be judged against the argument rather than against the code.

---

## Get the cross-model second opinion on the Stage 6 spec

**Surfaced by:** the engineering review of `specs/stage6_propensity_and_weights.md`, 2026-08-14 — as the
one piece of that review's own scope it failed to deliver. Recorded in that spec's §19 and §20.

**What.** Re-run an independent second-model review of the Stage 6 spec and fold anything it finds into
that document's §20 as a second round.

**Why.** §19 of the Stage 6 spec asked for two things that execution cannot settle: an engineering review
and a cross-model second opinion. The review of 2026-08-14 delivered the first — twelve findings, §20 —
and could not deliver the second. The `codex exec` pass failed with `ERROR: Failed to refresh token: 401
Unauthorized` **mid-run**, after it had already streamed most of the document, and no substitute was
dispatched. So half of what §19 asked for is genuinely outstanding, which §19 and §20 both now state
rather than glossing.

The expected yield is not zero and the precedent is in this repository: Stage 4's review found ten
defects in a document of this shape, Stage 5's found eight, this one found twelve, and in all three cases
some of what the review missed was caught during implementation. A second reader is cheapest **now**,
before `model.py`, `propensity.py` and two test files are written against the spec.

**Pros.** Two models disagreeing on a 2700-line specification is a stronger signal than one model
agreeing with itself. The failed run got far enough to show the pass itself is viable — the failure was
auth, not capability or context length.

**Cons.** It is a re-run of work already attempted, so the marginal yield over §20's twelve findings is
unknown and may be zero. Needs `codex login`, a manual step. §20's findings plus Definition of done items
6-12 — written as break-it-and-watch instructions precisely because Stages 4 and 5 showed reviews miss
things — already cover the classes of defect a second reader usually catches.

**Context for whoever picks it up.** Two practical notes from the failed attempt: the installed `codex`
rejects `--enable web_search_cached` (`error: unexpected argument '--enable' found`), so drop that flag;
and the prompt is worth reusing — it fed the spec's first 30KB, named this review's nine principal
findings, and instructed the second model **not** to repeat them but to find what they missed. Reviewing
the *revised* spec rather than the original is the point, so the finding list should be updated to §20's
twelve before re-running.

**Depends on / blocked by.** `codex login`, or any second model given the spec. Do it **before T2 of
Stage 6 starts** — after that, a finding costs a code change as well as a spec change.

## Second review round on Stage 7 §7.2 and §12.9 — the audit entries and their rendering

**What.** A focused review pass over `specs/stage7_balance_and_overlap.md` §7.2 (the two `model`
entries, their `n`, their four table rules) and §12.9 (their acceptance criteria), reading the
rendered log rather than the spec. Fold anything found into the spec, and into `balance.py` if T5 has
already landed.

**Why.** §18b and §18c ran every code fence in that document — twice, before and after the review
round of 2026-08-18 — and both runs called `Audit.to_markdown()` and confirmed it renders. **Neither
read a cell back out of it.** The two `detail` strings were checked only for building, and no assertion
compared a rendered cell against the object it describes. That is the one part of Stage 7 no execution
has actually exercised, and §19 now points a second reviewer at it by name.

The precedent is exact and it is the same section number. Stage 6 §7.2 shipped an `overlap_weights`
entry carrying two tables of different widths; `propensity._padded` (`propensity.py:379-397`) exists
because `data._md_table` computes its widths from `rows[0]` and raised `IndexError` on the narrower
rows, while the spec's `[1:]` silently discarded the second table's header. That defect was found by
running the audit call, not by reading the spec — and Stage 6 §18c had already parsed every fence
without catching it. Stage 7's `_overlap_table` is 13 columns wide and `_smd_table` is 6, and both go
through the same renderer.

**Pros.** Cheap now: the entries are 60 lines of helper and two `detail` strings. The failure class is
known, documented, and has already cost this repository one unplanned helper. Catching it before T7
means a human reads a correct log rather than debugging a renderer.

**Cons.** The entries are genuinely easier to review against real rendered output than against a spec,
so part of this may be better done during T5 than before it. Some of it duplicates Definition of done
items 4 and 11, which already require the log to be produced under two hash seeds and read end to end.

**Context for whoever picks it up.** Start from the two rendered tables in §18c, which are the real
output on `cohort_frame()`'s cohort: `_overlap_table` is 6 lines and 13 columns, `_smd_table` 20 lines
and 6. The specific things no test yet asserts: that `_fmt` is the only formatter reaching either grid
(§7.2), that the verdict column holds only `yes`/`no`/`undefined` and never `1` (§7.2's `_fmt(True)`
argument), that an undefined SMD renders as the literal `missing`, that `_range` gives one `missing`
rather than `missing–missing` on an empty arm, and that both `n` values are what §7.2 says they are
rather than row counts. §12.9 specifies all of these; none of them has been run.

**Depends on / blocked by.** Nothing to start. Best done **before T7**, and ideally alongside T5 so
the assertions are written against a real `Audit` rather than an assembled module.

---

## Three errata in the Stage 7 spec, all found by running it — all fixed, 2026-08-18

**Surfaced by:** implementing `specs/stage7_balance_and_overlap.md`, 2026-08-18. None changes a
decision. The first two are in §12's fixtures and changed no shipped code; the third is in a §7
code fence and changed one character of `balance.py`. All three are recorded in §18, which the
document declares normative, and this entry is the account of how each was found.

**Kept rather than deleted** because each names a way a specification of this shape fails: a
fixture that exercises nothing while passing, a number truncated below its own stated tolerance,
and a conditional clause no available frame reaches.

**1. §12.0's `constant_everywhere()` did not exercise §5.3a — FIXED IN THE SPEC, 2026-08-18.** The
declared fixture was four records, two per arm, weights `(0.1, 0.9)` and `(0.3, 0.7)`. Both weighted
means come back as exactly 14.0, so the rejected weighted-mean form returns `0.0` on it and §12.4a's
companion — "a companion implementing §5.3a's rejected form is shown to return `nan` on it" — could
not have passed as written. The fixture passed while exercising nothing.

**The cause is the weights and not the arm size**, and an implementation report said the opposite
first. What decides it is whether an arm's `Σ(w·x)/Σw` division rounds:

```
  Σ[0.1, 0.9] = 1.0     14.0 / 1.0        →  14.0                  exact
  Σ[0.3, 0.7] = 1.0     14.0 / 1.0        →  14.0                  exact
  Σ[0.2, 0.8] = 1.0     …                 →  14.000000000000002    rounds
  Σ[0.3, 0.7, 0.11]     …                 →  13.999999999999998    rounds
```

Of 49 two-record weight pairs tried, 12 produce a difference — so two records per arm bite or do not
bite depending on which weights they carry, and the 3-versus-2 arm sizes of `cohort_frame()`'s cohort
are where the defect happened to be found rather than why. **Resolved by one number**: the control
arm's weights are now `(0.2, 0.8)`, the fixture stays four records, §12.4a's bullet stands unedited,
and §18 gains the measurement row. `test_only_SOME_weight_vectors_round_and_that_is_what_decides_the_branch`
pins the mechanism so the fixture cannot drift back to two exact vectors.

**2. §12.0.2 and §18 wrote `penumbra_ml`'s unweighted SMD as 1.630960; it is 1.630961 — FIXED IN
THE SPEC, 2026-08-18.** Measured 1.6309612667875129, so the spec's figure was the same number
truncated rather than rounded at the seventh digit. Not merely cosmetic: it is 1.27e-6 from the true
value, and §12.0.2 asserts the vector **to 1e-6**, so a test pinning the spec's own literal at the
spec's own tolerance failed. Every other row of the golden vector reproduces to every digit printed.
Both homes now read 1.630961 and §18's row states the convention — `_fmt`'s six significant figures,
correctly rounded — because this is the one place in the document where that slipped.

**What is owed.** Both are one-line spec edits and neither is a disagreement, so they want the PI's
sign-off rather than a unilateral edit — the posture the Stage 4 §4.2 item above takes for the same
reason. Whoever makes them should also strike §12.4a's "a companion implementing §5.3a's rejected
form … is shown to return `nan` on it [`constant_everywhere()`]" and repoint it at the uneven-arm
vector, so the spec and the suite describe the same fixture.

**3. §7's `_smd_detail` fence produced a run-on sentence whenever a row was undefined — FIXED IN THE
SPEC AND IN THE CODE, 2026-08-18.** The undefined-rows clause ended with `", ".join(undefined)` and no
terminator, so the log read:

```
  … 19 row(s) are undefined and are reported as such rather than as zero: age, …, penumbra_ml A
  row's `n` is its own denominator [§11] …
```

Found by doing what §21.5 said a third round should do — reading a rendered `detail` as prose — on
the constructed frame §12.7 already builds. It needed that frame: neither the workbook nor
`cohort_frame()`'s cohort has an undefined row, so the branch had never rendered in §18b, §18c, §18d
or either review round. §7's fence and `balance.py` both carry the full stop now, and §18 records it.

**§21.5 is CLOSED on its second half**, recorded in the spec as a new §21.5a. Both `detail` strings
have been read as prose under every branch they have — 0 and many over-threshold rows, 0 and all-19
undefined rows, 0 and 1 thin centres — and `test_balance.py` §12.9a asserts each branch's sentence
rather than only its count. §21.5's FIRST half stands: §16b's numbers are now verified by T6 having
run, so nothing there is owed either. The separate `TODOS` item above asking for a focused round on
§7.2 and §12.9 is therefore closed by this work, not deferred.

---

## Pin the primary AND the secondary estimates once the analysis is locked

**Surfaced by:** writing `extended_bridging/specs/stage8_primary_outcome_estimator.md`, 2026-08-21
(§4.3 and §15 there). **Extended to Stage 9, 2026-08-25** — see the Stage 9 half at the end.

**What.** Add the workbook's `β`, `exp(β)` and six `RD_k` as literal regression pins in
`extended_bridging/tests/test_outcome.py` §14.12, and to that spec's §20 verification table, replacing
the property assertions that stand in for them today.

**Why.** Stage 8's spec deliberately does not measure or quote the primary effect, and §4.3 is the
argument: every earlier stage recorded its decisions as taken before any outcome was examined by arm,
and specifying the estimator *after* seeing the estimate would have put three decisions — §5.3's
category collapse, §6's coefficient bound and §8's thresholds — on the wrong side of that line. So
§14.12 asserts properties (converged, ascending cutpoints, `|β|` inside the bound, monotone `RD_k`,
orientation agreeing across the two scales) and nothing else.

**The cost, which this item exists to repay.** There is no pinned regression number on the primary
effect anywhere in git. An edit that moves `β` in the fourth decimal without breaking convergence,
ordering, orientation or monotonicity passes the whole suite. §14.0.2's golden vector is a partial
mitigation only: its cohort has five records and four cutpoints against the workbook's ninety-two and
six, so it exercises the same arithmetic at a different scale and cannot witness a defect that appears
only at the workbook's.

**Why it is not a code change.** The objection in §4.3 is to computing the estimate *while specifying
the estimator*, not to pinning it afterwards. Once the PI has seen the estimate and the analysis is
locked, the numbers are legitimate pins and the objection has expired.

**Depends on / blocked by.** ~~Stage 8's T8, which computes the estimate for the first time~~ — **T8 has
run, 2026-08-24**, and the estimate is in `out/logs/audit_v7_july26_with_abs_contra_indication.md`,
which is gitignored and is the only place it exists. So the first half of this dependency is
discharged and **the PI is now the whole of it.** Do it as its own commit, touching `test_outcome.py`
and the spec's §20, §14.12 and §15 together — and delete §15's first bullet in the same commit, since
it describes a gap that will no longer exist.

**What T8 leaves for whoever picks this up.** §14.12's twelve property assertions are green on v7 and
each of them can fail: the fit converged on the likelihood criterion in 3 iterations, all seven
declared mRS levels are occupied so `fit.categories` is `MRS_LEVELS` and there are 6 cutpoints, `α` is
strictly ascending, `|β|` is inside `POLR_MAX_ABS_BETA` so **G7 did not fire on the real data**,
`odds_ratio` equals `exp(β)` to 1e-12, every `RD_k` is monotone in `k` within each arm, and the
orientation agrees between the two scales through the weighted mean. What none of them does is pin a
magnitude, which is exactly the hole this item closes.

### The Stage 9 half, added 2026-08-25 (Stage 9 §4.3, §16 item 5)

**The same hole, seven more estimates.** Stage 9 inherits Stage 8 §4.3's position and narrows it: no
*weighted* quantity appears in `specs/` — no weighted proportion, no risk difference, no odds ratio, no
`tau` — so there is no pinned regression number on any of the **seven binary** estimates either. The
marginal counts in Stage 9 §4.2 do appear, because they do not decompose by arm.

**What to add when the analysis locks.** For each of `C.BINARY_OUTCOMES`: the two weighted proportions,
`RD_w`, `OR_w`, and `tau` on the five augmented outcomes, as literal pins in `test_outcome.py`'s §15
sections and in Stage 9's §21 table.

**The mitigation that exists today, and its limit.** Stage 9 §15.0.2's golden vector is a 24-record
constructed frame — round 3 wrote it into the spec as code, so it is reproducible — pinning `RD_w`
0.2850627657, `OR_w` 3.2404025938 and `tau` 0.2233300473 / 0.0402879733 across two covariate lists.
Like Stage 8's, it exercises the same arithmetic at a different scale: 24 records against 92, a minority
cell of 11 against six values from 5 to 42, and a written-out propensity score against a fitted one. It
cannot witness a defect that appears only at the workbook's scale.

**Why the two halves are one item.** Both are the same decision by the same person at the same moment —
the PI seeing the estimates and locking the analysis — and both are discharged by one commit over
`test_outcome.py`. Splitting them would put two entries in this file whose triggers are identical.

**Status: blocked on the PI**, and on nothing else. Both halves.

---

## Re-measure the separation band at fewer than seven occupied mRS categories

**Surfaced by:** the first independent review round on
`extended_bridging/specs/stage8_primary_outcome_estimator.md`, 2026-08-21 (§21.1 item 10; filed in
§15 and handed to Stage 10 in §11 and in the roadmap's Stage 10 entry).

**Amended 2026-08-24 — the synthetic half of this item is DONE and the answer is that the bound
holds.** Stage 8 §6.3 re-measured across every cutpoint count a replicate can produce, 1800 fits each:
the largest legitimate `|β|` reaches **11.04** rather than 8.79 and the smallest degenerate one falls
to **18.98** rather than 18.81, so the safe interval is `(11.04, 18.98)` and 14.0 sits inside it with
27% clearance below and 26% above. Every bound from 12 to 20 drops exactly the same replicates, so
14.0 is a point on a plateau. **What this retired:** the claim that the region between the modes is
*empty* — it is sparse, and legitimate fits do land in the old `(8.79, 18.81)` at odd cutpoint counts
— and the value **10.0**, which is now positively excluded and which `test_config.py`'s
`11.04 < POLR_MAX_ABS_BETA < 18.98` is written to reject. **What is left open is the rate**, below.

**What.** Establish **how often [§10] actually meets each cutpoint count**, which decides whether the
K-dependence above is theoretical on this cohort or is deciding which replicates get dropped. Then
either close this item as theoretical, or amend `POLR_MAX_ABS_BETA` with the measurement behind the
change.

**Why.** Stage 8 §5.3 collapses the response to the levels carrying positive weight, so the number of
cutpoints is a property of the *sample*, and [§10] draws 2000 samples. A replicate that happens to draw
no death fits five cutpoints — and the re-measurement above shows the legitimate maximum genuinely
moves with the cutpoint count, from 8.82 at six to 11.04 at five. `G7` is what turns a degenerate fit
into a countable failure, so the bound partly decides **which replicates the percentile interval is
computed from**, and a bound calibrated on frames unlike the ones being dropped is the kind of error
that reads as rigour.

**Why the rate cannot be measured at Stage 8.** The sweep above varies the occupied level set
*synthetically*; how often **this cohort's own resamples** do it needs the replicate loop, which is
Stage 10. **And T8 makes the question sharper rather than softer**: measured on v7, all seven declared
mRS levels are occupied and the point fit has six cutpoints, so the collapse is inert on the full
cohort and every replicate starts from a frame that has every level. Whether a 92-record stratified
resample loses one is exactly what nobody knows. Stage 8 §11 therefore asks Stage 10 to report the
**distribution of `len(fit.alpha)` across replicates** beside its two failure counters, so the question
is answered from replicates already drawn rather than by a second synthetic sweep. If that
distribution turns out degenerate at six, this item closes as theoretical and should be closed
explicitly rather than left open.

**Not the same as the Stage 12 gap.** Stage 8 §15 already files that the bound is calibrated on a
one-column design while [§14a]'s standardisation model carries eight covariates. That one bites at
Stage 12; this one bites at Stage 10, which is next.

**Depends on / blocked by.** Stage 10's replicate loop and its cutpoint-count report.

**ANSWERED IN PART 2026-08-25 — the cutpoint-count report exists (Stage 10 §7.4) and the distribution is
NOT degenerate: 1907 replicates of 1998 fit six cutpoints and 91 fit five, 4.6%.** Stage 8 §11's stated
consequence was that *"if that distribution is not degenerate at six cutpoints, the drop rate is partly a
function of a bound measured on frames unlike the ones being dropped"*. **That half is closed and the
answer is that it does not bite: the drop rate is zero.** G7 fired 0 times in 1998 replicates and so did
non-convergence, so the bound drops nothing at any cutpoint count and there is no rate for its
calibration to be a function of.

What is *not* closed is the comparability question, which is a different one and is [§8]'s rather than
this item's: 4.6% of the primary's draws are a `β` on a five-cutpoint scale, and percentiles are taken
across them. Under proportional odds they are the same parameter; [§8] prescribes no test of that and
[§15] declines one. It is filed as Stage 10 §16 item 5 and it reaches a reader through [§16]'s
constant-shift statement, not through a re-measurement of the band.

**ANSWERED IN FULL 2026-08-27 — the Stage 12 half is measured too, and it is also negative.** The
paragraph above said the calibration gap *"bites at Stage 12"*, where the standardisation model carries
nine covariates rather than one. Measured over 2000 replicates of the ten-column [§14a] design
(Stage 12 §9.2): the largest `|β_treatment|` is **2.0016** and the largest `|γ|` **3.2119**, against
`POLR_MAX_ABS_BETA = 14.0` and against Stage 8 §6.3's measured empty band running from 8.7873 to
18.8055. **Every fit sits entirely below the band's lower edge**, so the bound is not merely unfired,
it is not approached, and the question of whether a band calibrated on one column transfers to ten does
not arise on this workbook.

**And Stage 12 found the reason the question is differently shaped there**, which is worth recording
because it is not a property of the data: the [§14a] reported quantities are averaged probabilities,
bounded in [0, 1] *by construction*, so a separated fit produces no unprintable estimate the way
Stage 8's `exp(β)` produces 6.5e15. Stage 12 §9.1 therefore applies the bound at **key granularity** —
it drops `beta`, which is reported as a model parameter, and keeps the twenty-two standardised keys,
which are unaffected. The bound protects a different thing there, so its calibration matters less.

**Trigger, rewritten.** A cohort or a design on which any coefficient exceeds **8.7873**, Stage 8
§6.3's largest measured non-degenerate `|beta|` — at which point a fit is inside the band and the
band's calibration decides whether it is dropped. G7 firing at all is the older and looser form of the
same trigger.

**Status: open**, narrowed twice — both halves measured, both consequences absent.

## Make the pipeline resamplable: `propensity.fit` raises on 26.0% of stratified replicates

**What.** `propensity._record_exclusion` (`propensity.py:400-432`) compares the number of **distinct**
excluded `case_id`s against the number of excluded rows and raises `SchemaError` when they differ. Its
docstring states the premise: *"a duplicated or missing case_id is forbidden by Stage 2's A2, so this is
a belt over those braces."* A2 forbids duplicates **in the workbook**. A patient-level bootstrap
replicate contains duplicates by construction, so the premise does not hold for the population of frames
[§10] feeds this function.

**Measured** (Stage 9 §12.3): the cohort has exactly one covariate-incomplete record, in Lugano's stratum
of 31. Over 400 stratified replicates, **104 raised — 26.0%** — against an analytic 26.4% for the
probability of drawing a given row at least twice from 31. The match identifies the cause exactly.

**Why it blocks.** Stage 8 §11 hands over to Stage 10 that *"`FitError` is the droppable failure and
`SchemaError` is NOT — a `SchemaError` from `primary` or `propensity.fit` is a bug in the resampler, not
a sparse replicate, and catching it would drop replicates for a reason that is not about the data."* So
Stage 10's three options as landed are all prohibited: crash on a quarter of its replicates, catch a
`SchemaError` it has prespecified it must not catch, or drop 26% of replicates for a reason that is not
about the data.

**Two candidate fixes**, and choosing between them is a decision about what a replicate's log means:
teach the guard that a replicate's rows are distinct draws, or give each drawn row a distinct name. Stage
9's replicate measurements adopted the second **provisionally, in a scratchpad**, which is why four rows
of its §21 carry a caveat.

**Not fixed at Stage 9** because it is a Stage 6 amendment justified by a Stage 10 requirement that is
not yet specified, and putting it in the Stage 9 commit would land a change to Stage 6's audit contract
on the strength of a rule nobody has written down.

**Depends on / blocked by.** Nothing. This is actionable now and is the first thing Stage 10 hits.

**RESOLVED 2026-08-25 by the Stage 10 spec (§5.2, §5.3), and the second candidate was chosen against.**
Recorded as a PI decision in Stage 10 §17; it is not a numbered DECISION, because it is a choice about
this stage's own code rather than about the analysis.
The resampler gives each drawn row a distinct `case_id`; `propensity.py` is not amended and
`_record_exclusion` keeps its full strength on the point estimate. The deciding argument is that the
guard protects a log somebody reads, and a replicate's `Audit` is a throwaway that is never written — so
weakening it there would weaken it where it does real work, while renaming satisfies it truthfully: the
k-th draw of a patient *is* a distinct row. Measured: 26.0% of replicates raising before, **0 of 400 and
0 of 2000** after. Recorded as a PI decision in Stage 10 §17 with its reversible alternative.

**And it was not the only blocker.** Stage 10 §5.4 found a second `SchemaError` on the same path — S8 —
which the rename does not touch; see the item below.

**Status: closed.**

## Decide whether a constant outcome on a replicate is `SchemaError` or `FitError` (Stage 9 S8)

**What.** Stage 9 §4.4's precondition S8 raises `SchemaError` when a binary outcome is constant on its
[§11] population. Stage 9 §16 item 3 calls this *"the sharpest open question the stage leaves"*.

**The case both ways.** For `FitError`: a constant outcome on a resample is a sparse replicate, not a
bug, and [§10] already has the mechanism — drop and count. `sich` has 5 events in 92 and is where it
would happen. For `SchemaError`, which is what Stage 9 chose: a constant outcome makes the odds ratio
undefined at both ends simultaneously, and Stage 9 §6.2 measured that `0/0` and `inf/inf` both come back
as `nan` from numpy **silently**, so the frame produces no estimate by any route and calling it sparse
understates it.

**Trigger.** The first Stage 10 replicate that hits it. If the rate is non-negligible, `FitError` is
almost certainly right and the change is one line plus Stage 9 §4.4 and §15.1.

**RESOLVED 2026-08-25 — DECISION 6 (PI): the trigger fired and the rate is not negligible.**
Measured over 1998 replicates: `sich` 16 (0.8%), `tici_2b_3` 3 (0.2%), **19 replicates in total, 1.0%**
— and by Stage 8 §11 a `SchemaError` may not be caught, so as landed this was a second blocker of the
same shape as the one above. S8 becomes a `FitError`. Stage 9's own reason for the other choice is what
settles it: it objected to *returning* a plausible number, and `FitError` returns nothing — it drops the
replicate and counts it, which is [§10]'s mechanism for exactly this. S6 and S7 stay `SchemaError`;
neither fired in 2000 replicates.

**Status: closed.**

## Establish whether `death_90d`'s and `mrs_5_6_90d`'s `max|beta|` tail leaves their intervals usable (Stage 9 §9.6)

**What.** Stage 9 prescribes **no bound** on the `m_a(X)` coefficient, on two grounds: there is no empty
band to calibrate one from — `sich`'s `max|β|` runs continuously from 0.8558 to 80.9825 across
replicates, 42.7% above 8, 15.1% above 14 — and none is needed, because `m_a(X)` is a nuisance whose
*predictions* enter `tau` and predictions are bounded in [0, 1], so there is no `exp(β)`-style tail.

**What is unestablished.** A separated `m_a` is an over-fitted one: its fitted values reproduce `Y`, so
the correction term removes signal rather than residual confounding. Stage 9 measured the cost on
synthetic frames — augmentation raised RMSE by 1.06x to 1.19x against the unaugmented estimator on a
constant-risk-difference truth — but **that construction did not reproduce the condition it was built to
study**: its `max|β|` reached only 9.458, 11.475 and 3.466, because independent standard-normal
covariates carry none of the centre structure that drives the workbook's near-separation. So the RMSE
table is measured on frames unlike the ones the tail comes from.

**AMENDED 2026-08-25 — this item named the wrong outcome, and Stage 10 §10.2 is why.** `sich` is
**unaugmented**: the specified estimator fits it no `m_a(X)`, so it has no nuisance coefficient and no
tail. Stage 9's measurement reproduces almost exactly — 41.8% above 8 and 16.2% above 14, against its
42.7% and 15.1% — but it is a measurement of a model the analysis does not fit, obtained by fitting one
anyway. The real tail is on the two augmented **safety** outcomes: `death_90d` reaches 59.54 with 8.4%
above 8, and `mrs_5_6_90d` reaches 52.70 with 11.6%. The `m_a(X)` `FitError` rate is 0.05%, not 0.3%.

Stage 9's structural argument is untouched and is still why no bound is added: `m_a(X)`'s *predictions*
enter `tau`, predictions are bounded in [0, 1], so there is no `exp(β)`-style tail in the estimate
however large the coefficient gets. What is unestablished is unchanged too — a separated `m_a` is
over-fitted, so its correction term removes signal rather than residual confounding.

**Trigger.** The replicate-level correlation between an outcome's `tau` draws and its `max_abs_beta`,
per outcome, computable from draws Stage 10 already takes. If `death_90d`'s and `mrs_5_6_90d`'s `tau` is
uncorrelated with their `max_abs_beta`, this closes as theoretical — and should be closed explicitly.

**Status: open**, and now actionable: Stage 10 reports the distribution (Stage 10 §10.2).

## TICI's reduced `m_a(X)` is five parameters only while USZ contributes no records

**What.** The [§8] amendment of 2026-08-07 specifies TICI 2b–3's reduced model as *"treatment + `center`
+ `atrial_fib` — five parameters"*. Measured (Stage 9 §7.2), that is correct — and correct by accident.
`center` has **four** declared levels (`FACTOR_LEVELS["center"]`, reference `HUG`), which is three
dummies; `model.design` then drops `center_USZ` as constant on every design because USZ contributes no
records to the ATO population (strata: HUG 41, Lugano 31, CHUV 21). So the fitted design is 4 columns
plus the intercept `firth` prepends = 5.

**Why it matters.** If USZ ever contributes a single record, the same declared model becomes **6
parameters against 6 non-events** — saturated — and the amendment's stated arithmetic silently stops
holding. Stage 9 §15.9's companion asserts the count against both frames so the contingency fails a
test rather than living in a paragraph, but nothing prevents it.

**Trigger.** Any workbook revision that adds a USZ record, or any change to `FACTOR_LEVELS["center"]`.
The response is a PI decision about the reduced specification, not a code change.

**Status: open**, contingent.

## Make `POLR_TOL` and `POLR_SCORE_TOL` relative to the weight total (Stage 10 §7.5)

**Surfaced by:** the Stage 10 spec's probe round, 2026-08-25. Stage 8 §11 predicted this exactly and
could not measure it: *"A stratified resample of 92 records has a similar `Σw` by construction, so
nothing is expected to move — but 'expected' is not 'measured' ... If a replicate's `Σw` came in an order
of magnitude low, the same absolute tolerance would be an order of magnitude looser relative to the
objective, and the visible symptom would be `iterations` dropping rather than anything failing."*

**Measured** (Stage 10 §7.5): `Σw` over `in_model` runs from **4.9045** to 37.4814 across 1998
replicates, median 23.5118, against the point estimate's 27.736623 — a factor of 5.7 at the low end, not
an order of magnitude but well outside "similar by construction". And `polr`'s iteration count runs
**1** to 4, where Stage 8 measured 4 to 5 on every frame it had. The predicted symptom is present in the
predicted direction.

**What is and is not established.** Nothing shows any replicate's `β` is wrong: every fit reported
`converged_on`, the score criterion is on the gradient, and a one-iteration Newton step from a good start
value is normal on a smooth small problem. What is established is that "nothing is expected to move" was
wrong as stated, and that a relative-tolerance form is now a question with evidence rather than a style
preference.

**Why it was not done at Stage 10.** A convergence tolerance is part of [§8]'s estimator and is refit in
every replicate; changing it changes every estimate. That is a [§8] decision, not a [§10] one.

**Trigger.** A replicate that fails to converge, or a Stage 12 design where `Σw` is materially larger —
[§14a] uses unit weights over a bigger population, so its objective is on a different scale again.

**THE SECOND HALF OF THE TRIGGER FIRED 2026-08-27, AND THE CONSEQUENCE DID NOT OCCUR.** [§14a]'s
objective is over `Σw = n = 104` unit weights against roughly 23 for the [§7] overlap-weighted cohort —
a factor of about 4.5, which is the change of scale this item names. Measured over 2000 replicates
(Stage 12 §13.4): **every fit converged, in 4 to 6 iterations, all on the likelihood criterion**, with
`first_step_norm` between 1.009 and 22.76 and zero non-convergence. An absolute tolerance on a
log-likelihood four and a half times larger did not become either loose or tight enough to matter here.

That is evidence and not proof — the concern is structural, and a scale change large enough to bite
would be one where the *ratio* of `POLR_TOL` to the per-observation contribution moved by orders of
magnitude rather than by 4.5. **Trigger, narrowed:** a population an order of magnitude larger than the
cohort, or the first replicate anywhere in the pipeline that exits on `POLR_MAX_ITER`.

**Status: open**, narrowed, and it is still a PI/[§8] decision rather than a code change.

## Give `model.FitError` a `code` field instead of classifying by message prefix (Stage 10 §7.2)

**What.** Stage 8 §11 requires the separation count reported separately from the convergence count, so
Stage 10 must classify a caught `FitError`. `FitError` carries nothing to classify by: it is
`class FitError(RuntimeError)` with no attributes (`model.py:89`), and all **sixteen** raise sites —
fourteen in `model.py`, two in `outcome.py` — identify themselves by the first token of the message
(`G6`, `G7`, `O1`-`O6`, `F3`, `F4`, `F6`, `polr:`, `Firth:`, plus `S8`). So `C.FAILURE_BUCKETS` maps
prefixes and `test_bootstrap.py` §15.6 scans both modules to assert every token found is in the map.

**Why the scan is a mitigation and not a fix.** It makes a *reworded* message fail loudly, which is the
failure mode that matters. It does not make the classification structural, and a stage that needed to
branch on a failure kind rather than count it would still be parsing English.

**That this is a real risk and not a hypothetical:** an earlier draft of the Stage 10 spec stated fourteen
raise sites and omitted `G6` from the bucket map entirely. The scan found it — before the scan existed as
a test, while it was still a paragraph being checked.

**Why it was not done at Stage 10.** Sixteen raise sites across `model.py` — the module Stages 6, 8, 9
and 12 all depend on — amended on the strength of a diagnostic.

**Trigger.** The first bucket that is wrong, or the first stage that needs to branch on a failure kind.

**Status: open.**

## Measure bootstrap coverage under near-separation (Stage 10 §8.4, §16 item 4)

**What.** Stage 10 §8.4 measures percentile-interval coverage on two synthetic designs — 0.9600
(MC se 0.0113) unconfounded and 0.9400 (se 0.0137) confounded — both bracketing nominal. Neither design
has **strata**, so `resample`'s stratification is unexercised by the coverage test, and neither is
**near-separated**.

**Why that is the gap that matters.** [§13]'s amendment of 2026-08-24 establishes that treatment in this
cohort is nearly determined by centre — 30 of 41 bridged at HUG against 2 of 31 at Lugano — which is why
[§7] prescribes a Firth-penalised propensity model and why overlap weights do not balance exactly. So the
one condition that most characterises what the [§7] model is actually fitting is the one condition the
coverage measurement does not reproduce. A percentile bootstrap over a near-separated propensity model
is where coverage would degrade if it degrades anywhere.

**Cost.** Stage 10's two designs are 90 000 fits each and 270 s each. A stratified, near-separated design
is the same order.

**Trigger.** Before the primary interval is reported, if anyone wants a coverage claim about *this*
analysis rather than about the procedure. Stage 10 §8.4 states its own limits, so nothing currently
overclaims.

**Status: open.**

## The bootstrap draw depends on the stratum label alphabet (Stage 10 §4.3)

**What.** `bootstrap.resample` visits strata in the stratum column's own `groupby(sort=True)` order,
so the sequence in which the
generator is consumed — and therefore the entire set of replicates — is a function of the seed **and** of
the stratum labels. Renaming a centre through `C.CENTER_RECODE` produces different, equally valid
replicates.

**Why it is not guarded.** The alternative is to range over `C.CENTER_ORDER`, a declared tuple, which is
what `data.py`'s byte-identity rules prescribe for column loops. But `resample` is general in its stratum
because [§14a] and [§14b] resample the same way over a population this stage never sees, and a general
function cannot name this study's centres.

**Consequence.** The run summary records the seed, and [§16] requires it. After a `CENTER_RECODE` edit
the seed alone no longer identifies the draw.

**Trigger.** Any change to `CENTER_RECODE` or `CENTER_ORDER`. The cheap response is to record the sorted
stratum labels in the run summary beside the seed.

**Status: open**, low consequence, stated so it is not rediscovered as a reproducibility failure.

## Fix the stale `data.py:255` citation in `propensity.py` — CLOSES with Stage 11

**What.** `propensity.py:413`, inside `_record_exclusion`'s docstring, cites `data.py:255` for
`Audit.record`'s identifier normalisation. That normalisation is at `data.py:271`. `data.py` is unchanged
since the Stage 9 spec commit (`5fe1750`), so the reference drifted between Stage 6 and now.

**Why it was not fixed at Stage 10.** Stage 10 §5's central claim is that `propensity.py` has a
**zero-line diff** — the blocker was fixed in the resampler precisely so that Stage 6's guard keeps its
meaning and its strength untouched — and Stage 10's Definition of done item 2 checks it. Spending that
claim on a comment is a bad trade.

**Trigger.** The next commit that edits `propensity.py` for any reason.

**Fired, 2026-08-27.** Stage 11 §4.3 builds the covariate seam in `propensity.py` — `_fit`, `fit`,
`fit_full`, `Propensity.spec` and eight privates — so the trigger is met and the citation is corrected
in the same commit. The correct line is `data.py:271`, verified. Stage 10's zero-line-diff claim is
about *Stage 10's commit* and stays historically true; its Definition of done item 2 is marked
historical rather than left as a live check that now fails.

**Status: closed** with the Stage 11 implementation commit.

## Parallelise the bootstrap (Stage 10 §18)

**What.** Stage 10's engine is serial and single-threaded: about **145 s** for `N_BOOT` = 2000 with the
shared-design optimisation, of which the seven binary outcomes are 123 s and `model.design` was 27.2 s
before it was shared once per replicate.

**Why not now.** 145 s does not justify the determinism risk. Every parallel implementation makes the
stream a function of the scheduling unless each replicate is separately seeded, and Stage 10 §4.3
declines per-replicate seeding for its own reason: it makes the draw a function of an index a later edit
to the loop can renumber, and the failure is silent.

**Trigger, rewritten 2026-08-27 because Stage 11 falsified the old one.** It read *"the [§13]
sensitivity suite, which is five more `run` calls and turns 145 s into something over ten minutes."*
Both halves are wrong. DECISION 4 promoted **one** row and not five; and that one is not a `run` — a
`run` also drives `outcome.secondary`, so Stage 11 §5.4 uses `bootstrap.replicates` with a primary-only
body. Measured: the arm 69.3 s, the two subgroups together 68.3 s, **137.6 s added** to a 130-145 s
pipeline. The trigger is now **Stage 12's bootstrap alone**, which resamples a larger population; the
right design there is still a seed *sequence* — `SeedSequence.spawn(N_BOOT)` — reproducible under any
scheduling and a different decision from the one §4.3 declined.

**THE TRIGGER FIRED 2026-08-27 AND THIS IS NOW THE LARGEST SINGLE COST IN THE PIPELINE.** Measured
(Stage 12 §2, §12.8): Stage 12's bootstrap is **~18.4 minutes** — one `replicates` call whose body runs
three arms per draw — of which the **hierarchical arm is 95%**, at 17.5 minutes for 2000 fits of
`model.polr_ri`. Against a 130–145 s Stages 1–10 pipeline plus Stage 11's 137.6 s, that is roughly an
eightfold increase in total wall clock, and it is concentrated in one loop over one estimator.

**Two things make it the good case for the seed sequence rather than a reason to panic.** The
hierarchical body is pure — one `model.design`, one `model.polr` for start values, one `polr_ri` — so
it parallelises without sharing state; and Stage 12 §12.4 already measured that the 21.7-hour version
of the same arm was an optimiser choice, so the remaining 17.5 minutes is not slack to be tuned away.
`SeedSequence.spawn(N_BOOT)` remains the right design, for §4.3's reason.

**Status: open**, and it is now actionable rather than deferred.

## Label the denominator on each row of the bootstrap audit grid (Stage 10 §16 item 9)

**What.** `_diagnostics_table`'s block in the `bootstrap_replicates` audit entry renders four
distributions and one per-outcome counter, over **three different denominators**, none of them labelled.

**Why it is a trap.** Stage 10 §15.15 asserts that `n_alpha`, `polr_iterations`, `n_in_model` and
`sum_w` sum to the **live** replicate count and not to `N_BOOT`, because a replicate killed by a
`propensity.fit` failure never reaches `outcome.primary` and contributes a `None` that is dropped rather
than counted as a zero. `or_corrected`'s denominator is that outcome's `Draws.n_attempted`, which is
different again. A reader comparing the six-cutpoint count of 1907 against `N_BOOT` = 2000 concludes 93
replicates fitted something else; the answer is 91 at five cutpoints plus however many never reached
`polr`.

**Why it is invisible today.** The propensity `FitError` rate is 2 of 2000 — 0.1% — so the three
denominators currently differ by 2. That is exactly the condition under which the first run where they
differ materially gets misread.

**Cost.** One column of labels. Stage 10 §10.3's constraints apply to it: no `|` in any cell, and every
row of the concatenated grid keeps the same cell count (§15.13).

**Depends on / blocked by.** Stage 10 T9, the audit entry and its grid.

**Status: open**, cheap, filed so the grid is not shipped ambiguous.


## A 12.5% bootstrap drop rate exists, and nothing prescribes what to do about it (Stage 11 §8.4, §16 item 1)

**What.** Stage 11's `unknown_onset` subgroup loses **250 of 2000** replicates: 9 to O6 rank deficiency
when a replicate draws no treated patient at witnessed onset, 241 to G9 when the fit converges and
returns a reported coefficient at or above `POLR_MAX_ABS_BETA`. Measured `|gamma|` up to **19.296**,
against a bound of 14.0, with `polr` taking 10-15 iterations where 3-4 is typical.

**Why it matters now and did not before.** `../out/questions_for_the_pi.md` carries a question deferred
with *"ask after Stage 10"* — whether a drop rate above some level invalidates the interval — recorded
as premature because **no drop rate existed**. Stage 10 then measured the primary's at 0.10% and the
worst per-outcome at 0.8%, which answered it by making it moot. This is the first number that makes it
live, and it is on a quantity [§13] labels hypothesis-generating.

**The cause, measured.** Five treated patients at witnessed onset, of 41, four of them at one centre.
[§13]'s amendment of 2026-08-10 retained the subgroup on the cohort split — 33 witnessed of 126 records
— which is the right number for the question it answered and is not the number that governs
estimability.

**What is NOT proposed.** Withdrawing the subgroup. That decision was taken before any subgroup estimate
existed and must not be revisited on one. Nor relaxing the bound: a separated fit's `exp(beta)` is a
finite float every downstream table accepts, which is why the bound exists.

**Trigger.** Now. The surviving count is printed beside the interval and the treated count beside the
odds ratio either way, so nothing overclaims while the question is open.

**Status: open**, put to the PI.

## Measure the interaction test's power (Stage 11 §16 item 3)

**What.** "Hypothesis-generating" is a word. The number that would make it mean something is the power
of the `gamma` test at this cohort's shape — 92 records, 41 treated, treatment near-determined by
centre — in Stage 10 §8.4's style, with the Monte Carlo standard error quoted.

**Why.** It converts *"no interaction was detected"* into *"the test had X% power to detect a doubling
of the odds ratio between levels"*, which is what the label is for. It would also close Stage 10 §16
item 4's trigger — a coverage design carrying centre-like near-determination, *"the one that would say
something about this study rather than about the procedure"* — for one parameter.

**Cost.** Of the order of Stage 10 §8.4's two designs: 90 000 fits and about 270 s each.

**Trigger.** Before the subgroup table is read as evidence of no interaction.

**Status: open.**

## `Primary` carries no specification, so two of them are structurally identical (Stage 11 §16 item 5)

**What.** After Stage 11 the pipeline holds two `outcome.Primary` objects — the [§7] estimate and the
[§13] arm's — with no field distinguishing them. `sensitivity.Arm` binds each to its `Propensity` and
its `__post_init__` checks the pair, but the binding lives outside `outcome.py`, by the roadmap's own
choice that *"`outcome.primary` then runs unchanged"*.

**Why it is not fixed now.** Two objects, one of which is always inside an `Arm` record that asserts
the pairing. A `spec` field on `Primary` would be a third place the specification is recorded and a
third that can disagree.

**Trigger. The second sensitivity arm.** With three `Primary` objects circulating the field should move
onto `Primary`, and `Arm` should read it rather than assert it.

**Status: open**, deferred with a named trigger.

## The draws-to-intervals loop exists in two modules (Stage 11 §3.2, §16 item 7)

**What.** `bootstrap.run` turns a `dict[str, Draws]` into a `dict[str, Interval]` at
`bootstrap.py:1105-1111`, and `sensitivity._intervals` does the same for the [§13] arm and the [§13]
subgroups. One `C.ci_min_draws()` floor, one `percentile_ci` call at `C.CI_LEVEL` and one p rule now
exist in two places. The two are not the same function because the p rule differs — `run` uses
`_tested`, the arm passes `lambda key: False` and the subgroups pass `_tested_subgroup` — which is
why `_intervals` takes the rule as an argument and `run` does not.

**Why it is not fixed now.** Refactoring `run` onto the shared helper would edit the function whose
twenty-six intervals Stage 10's Definition of done asserts, and Stage 11 §0.2 refits nothing Stage 10
refit. The duplication is two callers of a five-line loop, and both are covered: §15.3 asserts the
arm emits no p and §15.8 asserts the subgroups emit exactly two.

**Trigger.** A third caller — Stage 12's standardisation bootstrap is the candidate — or any change to
`C.ci_min_draws` or `C.PERCENTILE_METHOD`. At that point the loop moves into `bootstrap.py` as a ninth
public name taking the p rule as a parameter, and both callers read it.

**CLOSES WITH STAGE 12 — the third caller is real and the trigger's own prescription is what Stage 12
§13.3 specifies.** `standardise.inference` passes `lambda key: False` over seventy keys, so the loop
would be in three modules rather than two. `bootstrap.intervals(draws, tested)` becomes the ninth
public name, `run` and Stage 11's two callers move onto it, and the p rule stays a parameter for the
reason this item already gives: it is the only thing the three callers disagree about.

**The guard this item asked for is specified too.** *"Refactoring `run` onto the shared helper would
edit the function whose twenty-six intervals Stage 10's Definition of done asserts"* — Stage 12 §20.11
asserts those twenty-six are **byte-identical** before and after, and Stage 12 §25 orders the
extraction as implementation step 3, *before* Stage 12 has any keys of its own, so a discrepancy is
unambiguously the refactor's rather than the new stage's.

**Status: closed** with the Stage 12 implementation commit.

## `Replicate` cannot describe a body with more than one `polr` fit (Stage 11 §3.1, §16 item 8)

**What.** `bootstrap.Replicate.n_alpha` and `.polr_iterations` are scalars (`bootstrap.py:148-149`),
which is right for a body running one ordinal fit — [§10]'s and the [§13] arm's. Stage 11's subgroup
body runs two, one per `C.SUBGROUPS` key, so it sets both to `None` and `sensitivity.Subgroups`
carries no `bootstrap.Diagnostics` at all. The cost is real and named: the per-subgroup cutpoint
distribution is not recoverable from the log, so [§16]'s constant-shift statement cannot be evaluated
per subgroup from a shipped counter — Stage 11 §21 carries the one measurement that exists.

**Why it is not fixed now.** Widening the two fields to `dict[str, int | None]` keyed by fit touches
`Replicate`, `Diagnostics`, `_diagnostics`, `_diagnostics_table` and every Stage 10 assertion over
them, to gain a diagnostic for a hypothesis-generating analysis that Stage 11 §8.8 already labels as
such. Stage 11 reports what it needs off `Draws`, which carries the surviving count and the buckets.

**Trigger.** The first stage that wants diagnostics from a multi-fit replicate body. [§14a]'s
standardisation is the candidate, since it fits per arm.

**THE TRIGGER FIRED 2026-08-27 — AND THE REASON WRITTEN ABOVE IS WRONG, WHICH IS CORRECTED HERE
RATHER THAN QUIETLY REPLACED.** *"[§14a]'s standardisation … fits per arm"* describes a mechanism that
does not occur. Standardisation fits **once** and predicts **twice**, by overwriting the treatment
column of the already-fitted design — Stage 12 §7.1, where re-fitting per arm is shown to be not merely
wasteful but silently wrong, because `model.design` drops the treatment column as constant on a frame
where every patient is treated. So the pooled and support-restricted arms share one `polr` fit and one
`n_alpha`, and this item's stated route to a multi-fit body never opens.

**What does open it is the hierarchical arm**, which adds a **second estimator** — `model.polr_ri` —
to the same body (Stage 12 §12). Two fits, two convergence stories, one scalar pair on `Replicate`.

**And the fix is still declined.** Widening `n_alpha` and `polr_iterations` to dicts touches
`Replicate`, `Diagnostics`, `_diagnostics`, `_diagnostics_table` and every Stage 10 assertion over
them, to gain what Stage 12 §3.3 already provides better: `hier.sigma` is a **first-class estimand
key**, so the random-intercept fit's own parameter has draws, a surviving count and an interval,
recoverable from `Draws` rather than from a diagnostic. `at_floor`'s rate is reported the same way.

**Trigger, rewritten.** A stage wanting the **cutpoint distribution** of a second fit — which
`hier.sigma` does not carry and no key can, since it is a property of the fit's shape rather than a
number the fit produces.

**Status: open**, trigger fired, fix declined with the reason restated.

## A percentile interval on a parameter sitting at its boundary (Stage 12 §12.5, §21 item 7)

**What.** `model.polr_ri`'s between-centre SD reaches `POLR_RI_SIGMA_FLOOR` in **22.0%** of replicates
— measured over 200 — because [§14a] names `σ²_C = 0` as a legitimate answer and four clusters
frequently give it. The `hier.sigma` key therefore has a percentile interval whose lower limit **is the
floor**, and a floor is not a value: it means "the variance collapsed", which is a qualitatively
different statement from "the variance is small".

**Why it is not fixed now.** There is no cheap correct answer. The literature's options — a
likelihood-ratio interval against the boundary mixture `0.5·χ²₀ + 0.5·χ²₁`, a parametric bootstrap, or
a Bayesian fit with a prior on σ — are respectively a second inference procedure, a second bootstrap
inside the first, and the prior-driven route [§14a] explicitly declines (*"a Bayesian fit would be
steadier but prior-driven"*). Stage 12's answer is this pipeline's standing one: **count it, print it,
and put it to the PI.** Stage 12 §17.2 requires the boundary rate printed beside the interval.

**Trigger.** Any reading of `hier.sigma`'s interval as a range for the between-centre SD, or a
manuscript reviewer asking for a test of `σ²_C = 0`.

**Status: open**, put to the PI, and it is a reporting decision before it is a code one.

## `RD_5` and `RD_4` coincide in 4.15% of the [§14a] draws (Stage 12 §6.3, §21 item 8)

**What.** mRS 5 carries three patients of 104, so a stratified resample loses the level entirely in
**4.15%** of replicates — measured, 83 of 2000, and it is always that level. In those replicates the
fit has five cutpoints, the standardised distribution carries a structural zero at mRS 5, and
`RD_5 == RD_4` **exactly**. The sampling distribution of `RD_5` is therefore narrowed by an amount
nobody has quantified, in a direction that makes it look more precise than it is.

**Why it is not fixed now.** The structural zero is the correct value — a level no record occupies is
a level for which the fit provides no cutpoint, and interpolating one invents a parameter (Stage 12
§6.3). The alternative of dropping collapsed replicates would select the bootstrap on the outcome
distribution, which is worse. Stage 12 counts and prints the rate.

**Trigger.** A manuscript reporting `RD_4` and `RD_5` as separate findings, or any comparison of their
interval widths.

**Status: open**, and it is a reporting caveat rather than a defect.

## `polr_ri` ships without an independent implementation having agreed with it (Stage 12 §12.7, §26)

**What.** `model.polr_ri` is a new maximiser and Stage 12 §12.7 specifies its `ordinal::clmm` oracle in
full — the three-way sign asymmetry, `nAGQ = POLR_RI_NODES`, the separate tolerance on σ. **The oracle
has not been run.** R segfaults at startup on this machine for the `uname`/PATH reason
`test_reference_r.py`'s banner documents at length, so the spec's verification of `polr_ri` rests on
its own analytic gradient, its node-count stability, and its collapse to `model.polr` at σ = 0 — all
internal.

**Why it matters more here than for the other four oracles.** `logistf`, `PSweight`, `clm` and
`OrderedModel` check estimators whose specifications were already fixed by [§7] and [§8]; this one
checks an estimator **this repository designed**, including its quadrature rule and its boundary
handling. It is the only estimator in the pipeline in that position.

**Trigger.** Stage 12's implementation step 7. Opening the gate is the deliverable, not a nice-to-have,
and `_r_environment` is the fix the environment needs.

**Status: open**, and it blocks calling the hierarchical arm verified.

## The [§14a] transported population is 12 and not 14 (Stage 12 §4.1, §21 item 9)

**What.** USZ contributes 14 eligible patients and **two of them have no 90-day mRS**, so the [§11]
estimation population carries 12. Both outcome-missing records in the whole [§14a] population are at
USZ — measured. The centre whose inclusion is [§14a]'s entire reason for existing loses one seventh of
its contribution to outcome ascertainment, and the transport rests on the 12 that remain, of whom
three are outside the treated support (Stage 12 §10.2).

**Why it is not fixed now.** [§11] is complete-case by prescription and no imputation is permitted
anywhere in this plan. The number is reported per centre in Stage 12's population audit entry, which
is what [§11] requires.

**Trigger.** A workbook in which the never-IVT centre's outcome completeness is materially worse, at
which point the transported population is small enough that the arm's interpretation changes.

**Status: open**, informational, and it is a denominator rather than a defect.

## Whether [§14a]'s hierarchical standardisation should be marginal rather than conditional (Stage 12 §12.6)

**What.** A random-intercept model offers two standardisations: integrate the prediction over
`b ~ N(0, σ²)`, or predict each patient at their own centre's conditional mode `b̂_c`. Stage 12
implements the second, because [§14a]'s *"patients at a never-IVT centre inform that centre's intercept
through their direct-EVT outcomes"* is only true under it. The first answers a different question —
the distribution for a patient at a *randomly drawn* centre — and a reviewer may well want it.

**Why it is not decided here.** It is a different estimand and [§14a] named one. Adding it is a [§14]
amendment, and it would need its own keys, its own interval and its own place in [§16]'s
different-populations statement.

**Trigger.** A reviewer asking what the arm says for a patient at an unobserved centre.

**Status: open**, and it is the PI's rather than a code change.

## Stage 12's spec is written against an unlanded Stage 11 (Stage 12 header, §21 item 11)

**What.** Stage 12's specification depends on Stage 11 in exactly three places, all named in its
header: `bootstrap.bucket` and `bootstrap.collect` being public (Stage 11 §5.4); the interval loop
being extractable as Stage 11 §16 describes; and `outcome.estimation_population` existing and
deliberately **not** being used (Stage 12 §4.1). Stage 11's spec is written and its code is not built.

**Why it is filed rather than resolved.** Stage 12's spec could not wait on Stage 11's implementation
without the roadmap's order becoming a serialisation of writing as well as of building, and the three
dependencies are small and named. But they are assumptions about code that does not exist.

**Trigger.** Stage 11 landing with `bucket`/`collect` private, or with the interval loop left
duplicated. Either makes Stage 12 §13.3 and §13.2 wrong rather than merely early.

**Status: open**, and it closes when Stage 11 lands.
