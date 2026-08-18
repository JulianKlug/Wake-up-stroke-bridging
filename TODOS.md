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
