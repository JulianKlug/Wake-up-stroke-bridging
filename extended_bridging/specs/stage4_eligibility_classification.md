# Stage 4 spec — eligibility classification

Implements roadmap Stage 4 [§3, §11]. Section references in brackets are to
`statistical_analysis_plan.md`. Numbers and decisions referenced as DECISION *n* are established in
Stage 0 and recorded in `../out/stage0_data_inventory.md`. Constants named in `SMALL_CAPS` are Stage
1's and live in `config.py`; `stage1_config_and_data_contract.md` is their specification. The frame
this stage receives is specified by `stage3_derived_variables.md` §11.

**Status.** Written 2026-08-10, against the frame and the workbook as verified in §18. **Revised
2026-08-12** after the engineering review recorded in §20, which found nine defects; all nine are
folded in above and every claim in §18 was re-run rather than re-read. It is the **sole source for
the Stage 4 implementation**: everything the implementer needs is here, and anything not here is not
to be invented.

**Two questions were put to the PI before this document was written**, because both change what the
code is rather than how it is written, and both are recorded where the code will read them:

1. **Whether [§14a]'s "all eligible patients" includes the indeterminate group.** It does. One
   predicate — *retained means not ineligible* — governs [§3] and [§14a] alike (§5). 43 of 126
   records turn on it.
2. **Whether the eligibility entry needs a new audit `kind`.** It does not; it reuses `derivation`
   (§9.1), so `data.py` is untouched by this stage and Stage 5 keeps a free hand for its cohort-flow
   table.

**AMENDED 2026-08-13 by DECISION 1a — eligibility is two-valued and the reason column is read by
nothing. §21 is the ledger and it supersedes this document wherever the two disagree.** The rule
blocks, the code blocks and the registry below are rewritten in place, because an implementer builds
from those. The surrounding *prose* about the three-class rule is left standing as the dated record of
how the module was specified and built, with a supersession pointer at each affected section — the
treatment [§3] and [§13] give their own amendments. Where a section is marked superseded, read §21.

**Nothing about Stage 4 is settled anywhere but here.** If a decision was made about this stage and
is not in this file, it is not a decision.

**Goal.** Every patient carries one of three eligibility classes, decided from the flag and the
presence of a reason under DECISION 1, with the one combination that would silently absorb a
contradiction raised on instead of resolved, and with the retained set declared once so that
`== "eligible"` and `!= "ineligible"` cannot mean different things in different stages.

**Not in scope.** Dropping any patient or any centre (Stage 5), the cohort-flow table (Stage 5), and
choosing a regime for an indeterminate patient under [§14b] (Stage 13). Stage 4 **adds one column**.
It drops no row, drops no column, and edits no value that Stage 2 or Stage 3 delivered.

---

## 0. Where Stage 4 sits

```
  data.load()  →  (df, audit)          25 analysis columns, 126 rows   [Stage 2 §11]
                       │
                       ▼
  derive(df, audit)                    31 columns                      [Stage 3 §11]
                       │
                       ▼
  ┌────────────────────────────────────────────────────────────────────────┐
  │  STAGE 4 — eligibility.py                reads config.py, data.Audit   │
  │                                                                        │
  │   classify(df, audit)                                                  │
  │     ├─ _assert_classifier_inputs(df)   E1…E5 → SchemaError             │
  │     ├─ _classify(df)                   the four cases, every mask      │
  │     │                                  filled, labels from config      │
  │     └─ _crosstab(df)                   centre x arm x class, every     │
  │                                        cell rendered; the cells must   │
  │                                        reconcile to len(df) or raise   │
  │                                                                        │
  │   retained(df)                         != ineligible, from             │
  │                                        ELIGIBILITY_RETAINED   [§5]     │
  └────────────────────────────────────────────────────────────────────────┘
                       │
                       ▼
              32 columns  →  Stage 5 restrictions  →  derive_cohort  →  33 columns
                       │                                                     │
                       │  the UNRESTRICTED frame keeps the column,            ▼
                       └─→ Stage 13 [§14b]                            Stages 6-11
                           the only analysis whose population
                           includes ineligible patients
```

`classify` appends one entry to the `Audit` that `load()` created. It writes no file.

### 0.1 Why a module of its own

Three reasons, in descending order of how much they would cost to get wrong.

- **Stage 13 needs this column on a frame Stage 5 has already thrown patients away from.** [§14b] is
  "the only §14 analysis whose population includes contraindicated patients", and its active regime
  assigns treatment *as a function of eligibility*. So the classification has to exist on the
  unrestricted 126-record frame and survive independently of cohort construction. A classifier living
  inside a Stage 5 cohort function would compute the one thing [§14b] needs and then discard the rows
  that need it.
- **Stage 3 has already declared it will not do this.** `derive.py`'s §14 states "Which patients are
  eligible. Stage 4. Stage 3 reads neither `ivt_contraindicated` nor `contraindication_reason`," and
  its shipped module docstring draws the arrow out to Stage 4. Folding the classifier into `derive.py`
  would make a shipped docstring false on the day it is edited.
- **The repository's unit is one module, one spec, one test file.** Stages 1, 2 and 3 each hold that
  shape and the acceptance criteria are numbered against it.

The module is small — one rule, five assertions, one table — and that is not an argument against it.
The rule decides who is in the analysis at all, and §18 records that a single mask in it moves 39
patients.

### 0.2 What Stage 4 does not touch

`classify` calls `df.copy()` before writing, as `derive` does, so a caller's frame is never mutated
underneath it. No existing column is edited, no row is dropped, no value is corrected. Unlike Stage 3
there is no private that writes in place, so the promise here is unconditional and needs no caller
contract.

**§12.14 asserts this promise directly rather than resting it on the `copy()` call**: the returned
column set is the input's plus `eligibility` and in that order, the shared columns compare equal
value for value, the row count is unchanged, and the caller's own frame is unchanged after the call.
Four one-line assertions for the sentence this whole section is about — a `copy()` that was moved,
dropped, or turned into a slice would otherwise be caught by nothing.

**Stage 4 has no dependency on Stage 3.** It reads `ivt_contraindicated`, `contraindication_reason`,
`ivt`, `center` and `case_id`, all of which Stage 2 delivers, and reads nothing `derive` produces. It
is placed after `derive` because Stage 3's shipped diagram places it there and one pipeline order is
better than two, not because it needs anything from there. §12.11 asserts it runs on a Stage 2 frame,
so the independence is a property rather than a claim — Stage 12 and Stage 13 both want it.

## 1. Deliverables

| Path | Contents |
|---|---|
| `extended_bridging/eligibility.py` | `classify`, `retained` |
| `extended_bridging/test_eligibility.py` | the acceptance tests in §12, and its **own** module-scoped `workbook` fixture (§10) — `test_data.py`'s does not cross files |
| `extended_bridging/config.py` | **amended** — `ELIGIBLE`, `INDETERMINATE`, `INELIGIBLE`, `ELIGIBILITY`, `ELIGIBILITY_RETAINED`; `ELIGIBILITY_ORDER` becomes a computed view; the `[§4]` citations corrected to `[§3]` (§10) |
| `extended_bridging/test_config.py` | **amended** — the new constants, 9.8's tuple list, `_RESOLVABLE`, and that `eligibility` reaches no covariate list (§12.12) |
| `extended_bridging/test_data.py` | **amended** — §12.14's two eligibility assertions **move** here (§10) |
| `extended_bridging/derive.py` | **amended** — two column counts in the docstrings, which Stage 4 makes wrong (§10) |
| `specs/stage3_derived_variables.md` | **amended** — the same two counts, in §0's diagram and §11's handover block |
| `extended_bridging/implementation_roadmap.md` | **amended, and already landed** on 2026-08-12 — Stage 4 gained its `**Spec:**` line, the retained-predicate reading, and an **Accept when** naming the every-cell-rendered and `!= ineligible` requirements. With this spec, not with the implementation; §17 T6 does not repeat it |

`data.py` is **not** amended. §9.1 gives the reason and §19 records it as a decision rather than an
omission.

`out/logs/audit_<label>.md` is an **output, not a deliverable**, as in Stage 2 §1 and Stage 3 §1:
gitignored, written only when a caller asks, and it names patients. Nothing under `specs/` may quote
a case identifier, and nothing in this file does — every patient below is a `HAND-N` fixture record
or a count.

## 2. Environment

Unchanged from Stage 2 §2.1. `uv`, Python 3.12, no new dependency — `eligibility.py` needs pandas and
nothing else. Commands run from `extended_bridging/`, flat module layout, `import config as C`.

```bash
cd extended_bridging && uv sync && uv run pytest -v
```

**What this costs to run.** One vectorised pass over at most 126 values for the classifier, five
boolean reductions for the assertions, 8 × 3 counts for the table, and one running sum over the same
eight cells for the reconciliation of §7.1. It is the cheapest stage in the pipeline. There is
nothing to cache and §15 says so as a standing instruction.

## 3. Module shape

```python
# eligibility.py
from __future__ import annotations

import pandas as pd

import config as C
from data import Audit
```

**No new types**, for Stage 3 §3's reason: everything this stage produces is a column or an audit
entry, and both already have types. The public surface is two functions:

```python
def classify(df: pd.DataFrame, audit: Audit) -> pd.DataFrame
def retained(df: pd.DataFrame) -> pd.Series
```

and the privates are `_assert_classifier_inputs`, `_classify`, `_crosstab`.

`retained` is public because Stages 5, 12 and 13 all call it, and because it is the single place the
[§3]-versus-[§14a] reading of "eligible" is decided (§5).

`eligibility.py` is **not** exempt from the Stage 1 §7 raw-name scan and must never become exempt.
`test_config.py` 9.4 walks `extended_bridging/**/*.py` minus a pinned three-file exemption set, so
this module is scanned from the moment it exists (§18). It names no raw header: `IVT_contraindicated_binary`
and `Contraindications_to_IVT` are read by their analysis names, and the free text of the second is
never read at all (§4.4).

### 3.1 The three pandas facts this stage turns on

Declared once here, verified by running them on pandas 2.3.3 as pinned by `uv.lock` (§18), because
this stage is one long chain of masks over nullable columns and the first two facts point in
**opposite** directions:

```
  (flag == 1)  on Int64 carrying <NA>          → <NA>                    propagates
  (flag == 0) & reason.notna()  with <NA>      → <NA> where reason is present
                                                 False where it is absent

  Series.mask(cond)      with <NA> in cond     → <NA> is treated as TRUE
                                                 the row IS masked, and takes the replacement
  df.loc[cond] / s[cond] with <NA> in cond     → <NA> is treated as FALSE
                                                 the row is NOT selected

  Series.isin(values)    with <NA> in the      → False, NEVER <NA>
                         series                  so ~isin SELECTS the missing row, and no fill
                                                 is needed — or wanted — on an isin mask
```

**The third fact is why `~…isin(…)` is the right shape for E4, E5 and `retained` and why none of
them carries a `.fillna(False)`.** `isin` answers a membership question, and a missing value is not
a member, so it returns `False` rather than propagating; negating it therefore *catches* the missing
row instead of silently passing it. A fill on an `isin` mask would be worse than redundant: on E5 it
would flip a missing exposure from caught to ignored, which is the opposite of what every other fill
in this module does. The two constructs look alike and behave oppositely, and §4.2 relies on that.
Verified, not reasoned (§18) — including the dtype asymmetry, which is real and harmless: `isin` on
the `string` `center` returns numpy `bool`, on the `Int64` exposure it returns pandas `boolean`, and
neither carries `<NA>`.

The asymmetry is the whole of §4.2 and §4.3, and it is worth being blunt about the direction of the
damage, because it is the opposite of the harmless one. Under the mask chain of §4.3 with the fills
removed, a record whose flag is missing is classified **`ineligible`** — verified, not reasoned
(§18) — and an ineligible record is one Stage 5 deletes. The failure is not a fabricated label in a
table. It is a patient silently removed from the analysis, with a full column, plausible counts and
every downstream denominator reconciling.

**One shipped comment states this wrongly and is corrected by §10.** `derive.py`'s
`_assert_onset_flags` says `.fillna(False)` "keeps `conflicting` free of `<NA>`, which `.loc` refuses
to mask with". On pandas 2.3.3 `.loc` does not refuse; it treats `<NA>` as False (§18). The sentence
is harmless there — `conflicting` carries no `<NA>` by construction — but it is the exact fact this
stage depends on, stated backwards, in the file an implementer will read for precedent.

## 4. The rule [§3, DECISION 1]

### 4.1 The four cases, and the one bit the reason column contributes

DECISION 1, restated once so the code has something to be checked against:

**Superseded by DECISION 1a (§21).** The rule is now two cases:

```
  ivt_contraindicated == 1                       →  ineligible     [§3] restriction 2
  ivt_contraindicated == 0                       →  eligible       including every patient whose
                                                                   reason was never documented
```

The three-class rule DECISION 1 specified, kept as the record of what was built and amended:

```
  treated                                        →  eligible       by revealed fact
  ivt_contraindicated == 1                       →  ineligible     [§3] restriction 2
  ivt_contraindicated == 0, a reason recorded    →  eligible       documented: no absolute
                                                                   contraindication
  ivt_contraindicated == 0, none recorded        →  indeterminate  [§3]: a blank is NEVER read
                                                                   as "no contraindication"
```

The **free text is never read**. `Contraindications_to_IVT` contributes exactly one bit — whether a
reason was recorded at all — and Stage 1's contract already says so in the column's own reason string.
Two consequences follow and both are load-bearing for what this stage tests:

- **There is no set of recognised reasons to validate a new value against**, so roadmap Stage 4's
  earlier "every observed reason string is classified" requirement has nothing to range over. What
  replaces it is §6.
- **The free-text normalisation problem does not exist.** The Stage 0 note records a `Clnician` typo
  and parenthetical annotations that would have made an exact-string classifier raise on four of five
  spellings of one reason. Nothing here can be broken by any of them.

**A note on the `[§4]` citations.** Roadmap Stage 4, Stage 0's DECISION 1 and three `config.py`
comments call `ivt_contraindicated` "the [§4] classifier". [§4] of the SAP is *Exposure* and is one
sentence long; it carries no classification requirement and no list of reasons. The rule this stage
implements is [§3]'s. There are **four** such citations in `config.py`, not three — lines 203, 412,
416 and 464 as the file stands (§18). Two of them are corrected to `[§3]` by §10, and the other two
disappear with §5.2's rewrite of the block they sit in, because both are inside it: line 412 is the
block's own header, `# --- eligibility [§3, §4, DECISION 1] ---`. The roadmap's two are left in place
and recorded in §13, because rewriting a stage entry the roadmap has already marked done is a larger
edit than the error justifies.

### 4.2 The assertion comes first

Before any label is assigned, and before the frame is even copied. Five checks, in three groups:
three decide whether the mask chain below can be trusted at all (E1, E2, E5), one is the
contradiction the *rule* would otherwise absorb (E3), and one is the contradiction the *table* would
otherwise absorb (E4).

```python
def _assert_classifier_inputs(df: pd.DataFrame) -> None:
    flag, treated = df["ivt_contraindicated"], df[C.TREATMENT] == 1
    bad: list[str] = []

    missing = df.loc[flag.isna(), "case_id"]
    if len(missing):
        bad.append(
            f"E1  ivt_contraindicated: {len(missing)} record(s) carry no classifier value: "
            f"{', '.join(sorted(missing))}. [§3] forbids reading an absent contraindication flag "
            "as 'no contraindication', and this stage will not invent a class for it.")

    outside = df.loc[flag.notna() & ~flag.isin((0, 1)), "case_id"]
    if len(outside):
        bad.append(
            f"E2  ivt_contraindicated: {len(outside)} record(s) outside {{0, 1}}: "
            f"{', '.join(sorted(outside))}. The classifier tests equality against both values; a "
            "third value would fall through to indeterminate and read as an undocumented reason.")

    conflicting = df.loc[(treated & (flag == 1)).fillna(False), "case_id"]
    if len(conflicting):
        bad.append(
            f"E3  {C.TREATMENT} and ivt_contraindicated: {len(conflicting)} treated record(s) are "
            f"flagged contraindicated: {', '.join(sorted(conflicting))}. One of the two is wrong "
            "and the revealed-fact rule would absorb the contradiction without a trace. Resolve it "
            "with the data owner; do not choose a class here.")

    off_centre = df.loc[~df["center"].isin(C.CENTER_ORDER), "case_id"]
    if len(off_centre):
        bad.append(
            f"E4  center: {len(off_centre)} record(s) carry a centre outside CENTER_ORDER: "
            f"{', '.join(sorted(off_centre))}. §7.1 renders one row per declared centre, so such a "
            "record would be classified and then vanish from the cross-tabulation, leaving a table "
            "whose empty cells no longer mean what [§3] restriction 1 reads them as meaning.")

    off_arm = df.loc[~df[C.TREATMENT].isin((0, 1)), "case_id"]
    if len(off_arm):
        bad.append(
            f"E5  {C.TREATMENT}: {len(off_arm)} record(s) whose exposure is missing or outside "
            f"{{0, 1}}: {', '.join(sorted(off_arm))}. The classifier reads the exposure twice — E3 "
            "and the revealed-fact mask — and both read an absent value as 'not treated', so a "
            "treated patient with no recorded exposure would be classified from the flag alone.")

    if bad:
        raise C.SchemaError(
            "\n  " + "\n  ".join(bad)
            + f"\n  {len(bad)} eligibility assertion(s) failed against {len(df)} records. [§3] "
            "classifies from time-zero information applied symmetrically across arms; none of "
            "these is a reason to relax that.")
```

**Collected rather than raised one at a time**, following Stage 2's `_assert_schema` and not Stage 3's
single check: five independent properties of four columns are worth reporting together, so a
corrected workbook is diagnosed in one run instead of five.

**No fill on E4 or E5, deliberately, and it is the one place in this module where a fill would be a
bug rather than belt-and-braces.** Both are `~…isin(…)`, and per §3.1 `isin` returns `False` for a
missing value rather than `<NA>`, so the negation already catches it. `.fillna(False)` on E5 would
convert "this patient's arm was never recorded" from *caught* into *ignored*. §4.3's masks are the
opposite case and carry the fill for the opposite reason. Two constructs, one line apart, that look
alike and must not be made uniform.

**E5 is checked after E3 although E3 reads the column E5 validates**, and the order is inert because
the messages are collected rather than raised: a frame that trips both reports both. What E5 changes
about E3 is worth stating anyway — E3's `.fillna(False)` reads a missing exposure as not-treated, so
on such a frame E3 *under-reports* rather than mis-reports, and E5 is what says so out loud.

**E3 is the assertion this stage exists to add.** It is the only one of the five that no runtime
check anywhere covers today — it lives solely as a data-gated test in `test_data.py` §12.14, which
runs on machines that have `data/` and nowhere else, and which cannot see a frame any other caller
builds. §10 moves it here, where it becomes a property of the classifier rather than a property of
one workbook.

**Every branch is unreachable on v7**, and that is recorded rather than treated as a reason to skip
one: §18 measures the flag as 0/1 on all 126 records, never missing, and never 1 on a treated
patient; the exposure likewise; and all 126 centres fall inside `CENTER_ORDER`, which is why §7.1's
table totals 126. Like Stage 3's `_assert_onset_flags`, these are written for the workbook that has
not arrived yet. §12.3 and §12.4 reach them on hand-built frames.

### 4.3 The build, with the labels declared rather than written

**Superseded by DECISION 1a (§21).** The shipped classifier is now one mask:

```python
def _classify(df: pd.DataFrame) -> pd.Series:
    flag = df["ivt_contraindicated"]
    out = pd.Series(C.ELIGIBLE, index=df.index, dtype="string")
    out = out.mask((flag == 1).fillna(False), C.INELIGIBLE)
    return out
```

Of the five properties below, four survive verbatim — the default is the class a record falls to when
no rule fires, no label is written in this module, `.fillna(False)` is load-bearing, and the column is
a `string` and total. The fifth, the revealed-fact mask, is deleted: see §21.

What DECISION 1 specified, kept as the record of what was built:

```python
def _classify(df: pd.DataFrame) -> pd.Series:
    flag, reason = df["ivt_contraindicated"], df["contraindication_reason"]
    out = pd.Series(C.INDETERMINATE, index=df.index, dtype="string")
    out = out.mask(((flag == 0) & reason.notna()).fillna(False), C.ELIGIBLE)
    out = out.mask((flag == 1).fillna(False), C.INELIGIBLE)
    out = out.mask((df[C.TREATMENT] == 1).fillna(False), C.ELIGIBLE)
    return out
```

Five properties, each load-bearing:

- **The default is `indeterminate`, and that is the [§3] requirement expressed structurally.** A
  record that matches no rule is one nothing is documented about, which is exactly what
  *indeterminate* means. Defaulting to `eligible` instead would read every blank as "no
  contraindication" — the one thing [§3] forbids — and no assertion could catch it, because the
  result would be a full column of plausible labels.
- **No label is written in this module.** `C.ELIGIBLE`, `C.INDETERMINATE` and `C.INELIGIBLE` are
  Stage 1's, and §5.2 declares them. Indexing into `ELIGIBILITY_ORDER` instead — `ORDER[0]`,
  `ORDER[2]` — would be worse than a literal: it couples the classification rule to a *display*
  order, so reordering the columns of a table would silently rewrite the rule. §12.6 scans for both
  mistakes.
- **`.fillna(False)` on every mask, and here it is not belt-and-braces.** Per §3.1 an unfilled
  condition carrying `<NA>` masks the row and takes the replacement, so a missing flag would be
  classified `ineligible` and deleted at Stage 5. The assertion in §4.2 makes that unreachable; the
  fills make it unreachable a second way, independently of whether E1 is ever weakened, and §12.4
  pins both cells of the square.
- **A `string`, not a `Categorical`**, for Stage 3 §4.4's reason about `onset_type`: Stage 5 removes
  a whole class from the frame, and a categorical that keeps a dead level makes every subsequent
  `groupby(observed=False)` resurrect `ineligible` as an all-missing row in the cohort tables.
- **The column is total.** Every record starts at `indeterminate` and every mask only overwrites, so
  no `<NA>` can appear. §6.3 states what that buys and §12.5 asserts it.

### 4.4 The order is not a precedence rule

The chain reads bottom-up as precedence — treated wins over ineligible, ineligible wins over
documented-eligible — and exactly one of those orderings can ever matter: treated **and** flagged, the
combination E3 has already refused. Written down rather than left in the reader's head, for the reason
Stage 3 §4.2 gives about its loop: if E3 were ever deleted, the last mask would win silently, which is
precisely why E3 is not optional.

The rest of the ordering is inert. Flagged-and-documented is not a conflict — §18 records that **all
19** flagged records carry a reason, so the second mask overwriting the first is the ordinary case,
not an edge one.

**And the last mask is the opposite of decoration.** §18 measures it: no treated patient carries a
recorded reason, so without the revealed-fact mask all **39** treated patients would classify
`indeterminate` rather than `eligible`. One line moves 39 of 126 records, and it moves them into the
class that [§3] retains anyway — so the primary cohort would be unchanged in *size* and wrong in
*composition*, with the indeterminate count rising from 43 to 82 and the limitation paragraph
[§3] demands describing a group twice its true size. §12.13 asserts the 39 directly.

### 4.5 `.notna()` is the whole of the reason column, and stripping it would be reading it

> **Superseded in full by DECISION 1a (§21).** The reason column contributes no bit at all: the
> classifier does not read it, so there is nothing to strip and A9 is no longer load-bearing here.
> The section is kept because its final paragraph — the A9 dependency — is the one a reader would
> otherwise still believe.

The bit DECISION 1 reads is presence. `.notna()` is presence. Nothing else in this stage touches
`contraindication_reason` — no strip, no case fold, no comparison, no membership test.

The temptation is to write `reason.str.strip().ne("")` for safety, and it must be resisted on two
grounds. The first is that Stage 2 already closed the hole: **A9 raises** on a whitespace-only cell,
with the message "a blank reads as 'a reason was recorded' and moves a patient from indeterminate to
eligible [DECISION 1, §3]" — so presence is already a clean bit on any frame that reached this stage.
The second is that Stage 2 §6 states the column "is not normalised at all: under DECISION 1 the free
text is never read, so there is nothing to normalise, and normalising it would create the appearance
that the text feeds a classifier." A strip here would create that appearance in the one module where
it is most misleading.

The dependency runs the other way from the usual one and is written down because of it: **A9 is
load-bearing for Stage 4, and it is the only Stage 2 assertion that is.** If A9 is ever relaxed, this
section is what says the classifier must be revisited in the same commit.

## 5. `eligible` is not the retained set

### 5.1 The predicate [§3, invariant 2]

Three classes, two of which are retained. Nothing in the SAP says "retained" in one word, so three
sections say it in three ways and the code must not:

- **[§3]** — "Where the reason a patient did not receive IVT was never recorded, eligibility is
  *indeterminate*; those patients are retained."
- **Roadmap invariant 2** — "No ineligible patient appears in any [§7] or [§14a] population."
  Phrased as an exclusion of `ineligible`, not as a restriction to `eligible`.
- **[§14a]** — "one proportional-odds model over all eligible patients from all centres", which
  taken as a bare word would exclude the indeterminate group.

**The PI's decision, taken before writing this spec: one predicate governs both — retained means not
ineligible.** [§3] and [§14a] draw from the same rule, and only ineligible patients are ever dropped.
So the code offers exactly one predicate and no stage gets to pick:

```python
def retained(df: pd.DataFrame) -> pd.Series:
    """[§3]'s retained set: every class but `ineligible`. Boolean, and total — never <NA>.

    The one place the difference between `== eligible` and `!= ineligible` is decided. It is 43 of
    126 records (§18) — 80% of the control arm in the primary cohort — so two stages resolving it
    differently would not look like a bug in either of them. Stage 5 restricts on this, [§14a]
    draws its population with it, and [§14b] uses the three-level column instead.

    Raises KeyError on a frame `classify` has not run on. That is the correct failure: a caller
    asking who is retained before anyone has been classified has a sequencing bug, and a predicate
    that answered anyway would answer about a different population.
    """
    return df[C.ELIGIBILITY].isin(C.ELIGIBILITY_RETAINED)
```

`isin` over the declared tuple, never `!= C.INELIGIBLE`. The two agree today and would diverge on the
one edit — a fourth class — that nothing else would catch, and the tuple is where a fourth class
would be declared. §12.7 asserts the equivalence *and* that a patched registry moves the predicate,
which is what proves it is registry-driven rather than a coincidence of the current three.

### 5.2 The registry

Replaces the current `ELIGIBILITY_ORDER` block, in the same file position — after
`outcome_model_covariates` and before `POST_TIME_ZERO`:

**Superseded by DECISION 1a (§21).** `INDETERMINATE` is removed, `ELIGIBILITY_ORDER` is
`(ELIGIBLE, INELIGIBLE)` and `ELIGIBILITY_RETAINED` is `(ELIGIBLE,)`; the shipped block is in
`config.py`. What DECISION 1 specified:

```python
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
```

`ELIGIBILITY_ORDER` was a literal tuple and becomes a computed view; it is read only by
`test_config.py` 9.8's immutability list, which asserts the type and not the contents (§18), so
nothing else moves.

## 6. What is asserted, and what is inherited

### 6.1 Five checks, and four of them have a Stage 2 counterpart

E1 (flag never missing) is Stage 2's A4b. E2 (flag in `{0, 1}`) is A5. E5 (exposure never missing,
and in `{0, 1}`) is A4 and A5 together. E4 (centre declared) is close to A3. Only E3 is new. They are
re-checked here deliberately, and the argument is not "in case A4b is wrong":

- **`classify` is called on frames Stage 2 never saw.** Stage 12 and Stage 13 take the unrestricted
  frame, every acceptance test below builds its own, and any future caller assembling a frame in a
  notebook gets no A4b at all. An assertion that only runs when the caller happened to come through
  `load()` is a guarantee about a code path, not about a function.
- **A5 alone permits missing.** Its own contract is "values in `{0, 1}` **or missing**"; A4 and A4b
  are what forbid the absence, and their coverage of these two columns is a two-element loop over a
  literal pair — `(("A4", C.TREATMENT), ("A4b", "ivt_contraindicated"))`. Nothing links that pair to
  the classifier. E1, E2 and E5 together make the classifier's precondition a property of the
  classifier.
- **E4 and A3 range over two different tuples.** A3 tests membership of `CENTER_RECODE.values()`;
  §7.1's table is built from `CENTER_ORDER`. The two are equal today — verified (§18) — but only by
  construction, and nothing asserts they stay equal. E4 checks the tuple the table actually uses,
  which is the only one that can make a row vanish.
- **The exposure is read twice by this module and validated by neither, before E5.** E3 and the
  revealed-fact mask both read it, and both read absence as "not treated". That is exactly the
  code-path-not-function gap the first bullet rejects, applied to the column that decides 39 of 126
  classifications (§4.4).
- **The cost is five boolean reductions over 126 rows**, once per run.

What is *not* re-checked here is A9's whitespace rule (§4.5). That one is left with Stage 2 because
checking it here would mean reading the text.

### 6.2 The classifier is total: no `<NA>`, ever

Stage 3's derived columns each inherit their source's missingness, and Stage 3 §8 extends the
missingness table to cover them. **Stage 4 does the opposite and it is worth saying why explicitly**,
because the obvious symmetry — a new column, therefore a new `absence_by_column` call — is wrong here.

The classification is total by construction: every record starts at `indeterminate` and the masks only
overwrite. A missing input does not produce a missing class; E1 raises on it instead. So the column
has no absence to classify, and a `missingness` entry over it would render one row reading `complete`
on every frame the pipeline can build — a table row that can never say anything.

What replaces it is an assertion: `eligibility` is never `<NA>` and every value is in
`ELIGIBILITY_ORDER` (§12.5). That is the stronger statement, and it is the one Stage 5 needs, because
a restriction written as `df[retained(df)]` over a column carrying `<NA>` would drop those patients
without counting them anywhere.

## 7. The cross-tabulation [roadmap Stage 4's acceptance criterion]

### 7.1 What it must show

Roadmap Stage 4: "a cross-tabulation of eligibility by centre and arm is produced". One table, in the
audit log, one row per centre × arm:

```python
def _crosstab(df: pd.DataFrame) -> tuple[tuple[str, ...], ...]:
    eligibility = df[C.ELIGIBILITY]
    header = ("centre", "arm", *C.ELIGIBILITY_ORDER, "n")
    rows: list[tuple[str, ...]] = []
    rendered = 0
    for centre in C.CENTER_ORDER:
        for code in sorted(C.TREATMENT_LABELS):
            cell = ((df["center"] == centre) & (df[C.TREATMENT] == code)).fillna(False)
            rendered += int(cell.sum())
            rows.append((
                centre, C.TREATMENT_LABELS[code],
                *(str(int((cell & (eligibility == cls)).sum())) for cls in C.ELIGIBILITY_ORDER),
                str(int(cell.sum()))))
    if rendered != len(df):
        raise C.SchemaError(
            f"\n  the cross-tabulation renders {rendered} of {len(df)} record(s). Every patient "
            "must fall in exactly one CENTER_ORDER x TREATMENT_LABELS cell. The absent ones would "
            "be classified and then missing from the table, and the `all / both` row would still "
            "reconcile against the frame, so nothing downstream would disagree.")
    rows.append((
        "all", "both",
        *(str(int((eligibility == cls).sum())) for cls in C.ELIGIBILITY_ORDER),
        str(len(df))))
    return (header, *rows)
```

**It takes the frame, not the frame and the Series.** `classify` assigns `eligibility` and then
tabulates from the column it is about to return, so the table is provably about the shipped column
rather than about a Series that happens to have been computed alongside it. One argument, one source
of truth, and the tabulated column and the returned column cannot be made to disagree by any future
edit that reorders the two statements.

Four rules, the first three of them the same rule `absence_by_column` follows for its per-centre
columns:

- **Every `CENTER_ORDER` × `TREATMENT_LABELS` cell is rendered, empty or not.** Never
  `pd.crosstab`, which drops absent combinations. The empty cells are the information: on v7 the
  table's `USZ / bridging` row is all zeros, and that row **is** the [§3] restriction 1 finding —
  the never-IVT centre, visible in Stage 4's own table, one stage before Stage 5 acts on it. A
  crosstab built from the data would omit the row and the finding with it.
- **The class columns come from `ELIGIBILITY_ORDER`** and the arm labels from `TREATMENT_LABELS`, so
  neither can be reordered or relabelled in one place only.
- **`sorted(C.TREATMENT_LABELS)`**, so the arm order is 0 then 1 — control before treated — in every
  interpreter. A dict iteration would be equally deterministic today and would stop being a
  *declared* order.

- **The cells must reconcile to `len(df)`, at run time, or the table is not produced.** This is the
  second of the two guards §4.2's E4 is the first of, and they are genuinely independent rather than
  belt-and-braces: E4 catches an *undeclared centre value*, this catches *any* reason a patient falls
  in no cell — a `TREATMENT_LABELS` that stopped covering the arm coding, a `CENTER_ORDER` shortened
  by an edit, a future third dimension. Verified on the raw `hand_frame()`, which is exactly this
  failure: three of its six records land in a cell, while E4 names only three of the six offending
  identifiers, because two of its raw centre codes happen to spell their own labels (§18). Neither
  guard alone sees the whole of it.

The `all / both` row lets a reader reconcile the table against the frame without adding nine numbers,
and §12.8 asserts that reconciliation too. Note what the `all / both` row cannot do on its own: it is
computed from the frame rather than from the cells, so a vanished patient leaves it *correct*. That
is precisely why the reconciliation is a check inside `_crosstab` and not merely a property a reader
could notice.

Verified rendering, on v7 under **DECISION 1a** (§18, re-run 2026-08-13). Every `n` and the
restriction-1 row are unchanged from the three-class table; the `indeterminate` column is gone and its
counts have moved into `eligible`, which is the table-level evidence that the amendment relabelled
records rather than moving them:

```
  | centre | arm       | eligible | ineligible | n   |
  | HUG    | EVT alone | 11       | 11         | 22  |
  | HUG    | bridging  | 30       | 0          | 30  |
  | CHUV   | EVT alone | 14       | 0          | 14  |
  | CHUV   | bridging  | 7        | 0          | 7   |
  | Lugano | EVT alone | 29       | 0          | 29  |
  | Lugano | bridging  | 2        | 0          | 2   |
  | USZ    | EVT alone | 14       | 8          | 22  |
  | USZ    | bridging  | 0        | 0          | 0   |   ← [§3] restriction 1
  | all    | both      | 107      | 19         | 126 |
```

Under DECISION 1 the same table read `64 / 43 / 19` across three class columns:

```
  | HUG EVT alone 11 0 11 22 · HUG bridging 30 0 0 30 · CHUV EVT alone 0 14 0 14
  | CHUV bridging 7 0 0 7 · Lugano EVT alone 0 29 0 29 · Lugano bridging 2 0 0 2
  | USZ EVT alone 14 0 8 22 · USZ bridging 0 0 0 0 · all both 64 43 19 126
```

### 7.2 What the table is not

It is **not** the per-centre completeness table. The pattern that produces the indeterminate group —
the reason field was collected at HUG and USZ and at neither other centre — is already in Stage 2's
`absence_by_column` row for `contraindication_reason`, which carries `kind = not recorded` and a count
per `CENTER_ORDER`. Rebuilding it here would give one fact two tables, and the one that drifts is
always the one nobody reads.

The two are complementary and the log reads in that order: Stage 2 says where the field was collected,
Stage 4 says what the classification made of it.

## 8. `eligibility` is not a derived name

`DERIVED_NAMES` is `("onset_type", *DERIVED_DICHOTOMIES, *SUBGROUPS)` and `eligibility` does **not**
join it, although it is a column no stage before this one can supply. The reason is a weakening that
Stage 3 §12.13 already had to repair once:

`test_config.py` 9.3's `test_every_ps_covariate_resolves` accepts any covariate that is
`in ANALYSIS_NAMES or in DERIVED_NAMES`. Widening `DERIVED_NAMES` therefore widens what may legally
appear in a covariate list — and [§14a] states plainly that "contraindication status is not a
covariate, the cohort being already restricted to eligible patients". A patient's eligibility must
never adjust anything; it decides who is in the population and nothing else.

So instead:

- `ELIGIBILITY` is its own constant (§5.2), like `TREATMENT`.
- `test_config.py`'s `_RESOLVABLE` gains `{config.ELIGIBILITY}`, so a future column-keyed constant
  naming it still resolves under 9.13.
- §12.12 asserts `ELIGIBILITY` appears in **no** covariate list — `PS_COVARIATES_FULL`,
  `BALANCE_ONLY`, `STANDARDISATION_COVARIATES`, or any outcome model's — exactly as Stage 3 asserts
  it of the subgroup keys.

It is not in `POST_TIME_ZERO` either, and that is deliberate rather than an omission: [§3] classifies
eligibility "using only information available at time zero", so it is a baseline attribute. Invariant
4's denylist is about variables that would be adjusted for; this one is barred from adjustment by the
assertion above instead, which is the stronger of the two because it names the property rather than
the timing.

## 9. The audit entry

### 9.1 The kind is `derivation`, and `data.py` is untouched

Decided with the PI. `KINDS` gains nothing, `_HEADINGS` gains nothing, `_MUST_NAME_CASES` is
unchanged.

The rule `KINDS` is organised on is *what the entry is about*, not which stage wrote it:
`derivation` means an entry describing a **column built from other columns**, which is what this is.
It names no patient — the one place Stage 4 names patients is §4.2's assertion, which raises, so such
a frame never reaches a log — so `_MUST_NAME_CASES` correctly does not cover it, for the same reason
it does not cover Stage 3's.

A new `cohort` kind was weighed and declined **for this stage**, not in general. Stage 5 has to log
rows *removed*, which no existing kind describes — `correction` is about values and `derivation` is
about columns — so Stage 5 will plausibly insert one. Inserting it here, one stage early, would put
Stage 4's classification under a heading whose only other content arrives later, and would spend a
`data.py` amendment on a decision that belongs to a spec not yet written. Stage 3 §9.2's repair
already made a later insertion a one-line change (the heading-adjacency test takes its neighbour from
`KINDS`), so nothing is lost by waiting.

### 9.2 The entry inventory

One entry. `step` is exact, and is `C.ELIGIBILITY` rather than the string, so the column and its log
entry cannot be renamed apart.

| kind | step | recorded by | `n` counts | table | `case_ids` |
|---|---|---|---|---|---|
| `derivation` | `eligibility` | `classify` | rows in the frame — every record is classified | 9 rows + header | — |

`detail`, with every count computed from the produced column and never from a second traversal:

> `[§3] eligibility under DECISION 1: treated → {ELIGIBLE}; ivt_contraindicated = 1 → {INELIGIBLE};
> flag = 0 with a reason recorded → {ELIGIBLE}; flag = 0 with none recorded → {INDETERMINATE}. The
> free text is never read — the reason column contributes one bit, whether a reason exists. {n_e}
> eligible, {n_i} indeterminate, {n_x} ineligible of {n} record(s). A blank is never read as "no
> contraindication" [§3]; the indeterminate group is retained and Stage 5 drops only the ineligible.`

`n` is the row count and not the eligible count, deliberately: this entry records that a
classification was made over the whole frame, and the three class counts are in the sentence and in
the table's last row. An `n` equal to one class's count would read as "63 records were affected",
which is not what happened.

The entry is unconditional — there is no frame on which it is skipped — so §12.9's inventory
assertion is a plain equality.

### 9.3 `classify`, written out

The public function is short and every line of it is load-bearing, so it is given here rather than
left to be reconstructed from §0's diagram and the table above.

```python
def classify(df: pd.DataFrame, audit: Audit) -> pd.DataFrame:
    """[§3] eligibility under DECISION 1. Appends one column and one audit entry; drops nothing.

    Read `ivt_contraindicated`, `contraindication_reason`, `ivt`, `center` and `case_id` — all of
    them Stage 2's — and nothing `derive` produces. §12.11 asserts that independence, which is what
    Stage 12 and Stage 13 rely on when they classify the unrestricted frame.
    """
    _assert_classifier_inputs(df)
    df = df.copy()
    df[C.ELIGIBILITY] = _classify(df)

    counts = df[C.ELIGIBILITY].value_counts()
    n_e, n_i, n_x = (int(counts.get(cls, 0))
                     for cls in (C.ELIGIBLE, C.INDETERMINATE, C.INELIGIBLE))
    audit.record(
        "derivation", C.ELIGIBILITY, len(df),
        f"[§3] eligibility under DECISION 1: treated → {C.ELIGIBLE}; ivt_contraindicated = 1 → "
        f"{C.INELIGIBLE}; flag = 0 with a reason recorded → {C.ELIGIBLE}; flag = 0 with none "
        f"recorded → {C.INDETERMINATE}. The free text is never read — the reason column contributes "
        f"one bit, whether a reason exists. {n_e} eligible, {n_i} indeterminate, {n_x} ineligible "
        f'of {len(df)} record(s). A blank is never read as "no contraindication" [§3]; the '
        "indeterminate group is retained and Stage 5 drops only the ineligible.",
        table=_crosstab(df))
    return df
```

Five things about those lines, each of which an implementer could reasonably get wrong:

- **The assertion runs before the copy**, so a frame that fails is never copied, and E1–E5 are
  reported about the caller's own frame. `_assert_classifier_inputs` only reads, so ordering it first
  is safe as well as cheaper.
- **The class names in `detail` are interpolated from `config`, never typed.** `{C.ELIGIBLE}` and its
  two siblings render as `eligible`, `indeterminate`, `ineligible` — the same words a literal would
  produce, and the reason §12.6's AST constant scan stays green. An implementer who writes the words
  into the f-string gets the right log and a failing test, which is the intended order of discovery.
- **The counts come from the produced column**, via one `value_counts()` over the column that was
  just assigned. Never from a second traversal of the flag and the reason, which would give the log
  and the frame two sources for one number.
- **`counts.get(cls, 0)`, not `counts[cls]`.** A frame in which no patient is ineligible — the
  fixture is one — has no `ineligible` entry in `value_counts()`, and indexing it would raise on the
  cleanest possible input.
- **`step` is `C.ELIGIBILITY`, not the string `"eligibility"`**, so the column and its log entry
  cannot be renamed apart. `kind` is the literal `"derivation"` because `data.KINDS` holds the
  literals and §9.1 is the decision that Stage 4 adds none.

## 10. What Stage 4 amends in Stages 1, 2 and 3

The full ledger, so that no amendment is discovered during implementation.

| File | Amendment | Why |
|---|---|---|
| `config.py` | add `ELIGIBLE`, `INDETERMINATE`, `INELIGIBLE`, `ELIGIBILITY_RETAINED`, `ELIGIBILITY`; `ELIGIBILITY_ORDER` becomes computed. All in the existing eligibility block, §5.2's text verbatim | the registries `eligibility.py` reads instead of literals |
| `config.py` | **four** comments say `[§4]`, at lines 203, 412, 416 and 464. Two change to `[§3]` by hand: the `IVT_contraindicated_binary` contract entry (203) and the `BINARY_COLUMNS` comment (464). The other two — the eligibility block's header (412) and its body (416) — go with §5.2's rewrite of that block | [§4] is *Exposure* and carries no classification rule (§4.1). Comment-only; no constant moves |
| `test_config.py` | 9.8's tuple list gains `ELIGIBILITY_RETAINED`; `_RESOLVABLE` gains `{config.ELIGIBILITY}`; the §12.12 assertions | 9.8's list and 9.13's union are explicit, so a new constant is outside both until added by hand (Stage 3 §12.13's warning, verbatim) |
| `test_data.py` | **move** `test_no_treated_patient_is_flagged_as_contraindicated` and `test_the_eligibility_classifier_is_never_missing` out of §12.14 and into `test_eligibility.py` §12.13 | both assert Stage 4's preconditions, and one of them is Stage 4's headline assertion (E3). Leaving them in `test_data.py` would mean the property is tested against one workbook and not against the classifier. Nothing else in §12.14 moves |
| `test_eligibility.py` | declares its **own** module-scoped `workbook` fixture — `load(WORKBOOK) → derive → classify`, under `DATA_GATED` | the two moved tests take a `workbook` fixture, and `test_data.py`'s (its line 628) is module-scoped and does not cross files. No `conftest.py` is added and nothing is imported from `test_data.py` for it: the repository has no conftest today, and a two-line fixture is not worth introducing one or making the scope act at a distance. §12.13 needs the classified frame anyway, which `test_data.py`'s fixture does not produce |
| `derive.py` | the module docstring's diagram (`32 columns → Stage 6`) and `derive_cohort`'s docstring (`a new frame of 32 columns`) | Stage 4 inserts a column between them: `derive` 31 → **`classify` 32** → `derive_cohort` **33**. Both counts are stale the moment this stage lands |
| `specs/stage3_derived_variables.md` §0, §11 | the same two counts, in the flow diagram and in the `derive_cohort` handover block | same |
| `implementation_roadmap.md` | Stage 4 gains its `**Spec:**` line; its **Accept when** gains the retained-predicate reading | precedent: Stage 3 §19 landed its spec line with the document, not with the code |

**Two things that look like they need amending and do not.**

- `test_derive.py`'s column-count assertion is `len(cohort.columns) == len(out.columns) + 1` —
  relative, so it survives the insertion unchanged (§18, read).
- `test_derive.py`'s Stage 3 audit-inventory test survives too, but for a **fragile** reason rather
  than a structural one: its `_full_pipeline()` helper does not call `classify`, so the tail slice it
  asserts on never spans this stage's entry. Verified by running it (§18). Do not "improve" it while
  here — §13 explains why, and `TODOS.md` carries the repair for Stage 5.
- `test_data.py`'s `_INVENTORY_ON_THE_FIXTURE` is the six entries `load()` appends and nothing else;
  Stage 4 appends to the `Audit` only when a caller invokes `classify`, which `data.load()` never
  does. Read, not assumed (§18).
- `data.py` — §9.1. No kind, no heading, no `_MUST_NAME_CASES` entry.

**One correction to a shipped comment, listed separately because it is a false statement rather than a
consequence.** `derive.py`'s `_assert_onset_flags` says `.loc` "refuses to mask with `<NA>`". It does
not (§3.1, §18). Correct the sentence to what pandas 2.3.3 does — `.loc` treats `<NA>` as False, which
is why the fill is belt-and-braces there and is **not** belt-and-braces here. The surrounding argument
for keeping the fill is unaffected, and no test asserts the sentence.

## 11. Data flow into Stages 5, 12, 13 and 14

```
  classify(df, audit) returns
   ├─ every column it was given, unchanged — no value edited, no row dropped
   ├─ eligibility       string, always one of ELIGIBILITY_ORDER, NEVER <NA>
   └─ 32 columns when called on derive()'s frame; 26 when called on Stage 2's

  retained(df) returns
   └─ bool, total; True for eligible and indeterminate, False for ineligible   [§5.1]

  audit
   └─ the same Audit object load() created. Stage 5 appends; nobody re-opens. No file
      written unless asked. The entries preceding Stage 4's one, counted by running the
      pipeline (§18) rather than carried:

        load()            6 on the fixture, 7 on the workbook
                          the 7th is `observation` / `onset_to_groin_999`, which fires
                          on no other source — so any positional assertion has to name
                          the source it counted
        derive()          4      onset_type, unknown_onset, dichotomies,
                                 absence_by_derived_column
        classify()        1  ←   Stage 4, here
        derive_cohort()   2      core_above_median, constant_covariates — these arrive
                                 at STAGE 5, after this entry, not before it

      Stage 3 contributes six entries in total, but only four of them exist when
      `classify` runs: `derive_cohort` is called by Stage 5, downstream of this stage.
      §12.9 asserts the position against the driver it prints, never against a
      remembered total.
```

What each later stage may rely on, and what each owes:

- **Stage 5** — restricts with `retained`, never with a comparison of its own, and calls
  `derive_cohort` *after* both [§3] restrictions. Its cohort-flow table must report the retained
  indeterminate count on its own line [roadmap Stage 5]; §18's number for it is **43** of the 54
  retained controls.
- **Stage 6 onwards** — never see an ineligible patient (invariant 2) and never adjust for
  `eligibility` (§8). The column is inert from here on for the primary analysis.
- **Stage 12 [§14a]** — draws its population with the same `retained` predicate, on the all-centre
  frame. §5.1 is the decision that makes that the same rule Stage 5 applies.
- **Stage 13 [§14b]** — takes the **three-level column**, not the predicate, because its active
  regime is a function of eligibility across a population that includes the ineligible. §14 records
  the one thing it must decide and Stage 4 does not.
- **Stage 14 [§16]** — every estimate carries its own denominator [§11], and the cross-tab in §7 is
  the source for the cohort-flow numbers that reconcile them.

## 12. Acceptance criteria

All tests live in `test_eligibility.py`, with section banners matching these numbers, as
`test_derive.py` does for Stage 3. Tests needing the private workbook reuse the `DATA_GATED` idiom;
everything else runs against hand-built frames. `test_eligibility.py` is not exempt from the raw-name
scan and builds its frames with analysis names.

The hand-built frame is `test_data.py`'s `hand_frame()`, with `corrupt()` for per-record variation —
imported, never re-declared, for Stage 3 §12's reason. It needs **no extending**: it already carries
both arms, both reason states and all four centres.

**But it carries `center` as a raw code, not as a label**, and every per-centre assertion in this
file turns on that. `CENTER_CODES = tuple(config.CENTER_RECODE)` is `("1", "Lausanne", "Lugano",
"USZ")`, and Stage 2's `centre_recode` is what maps them to `CENTER_ORDER`. Two of the four spell
their own labels and two do not, which is why the bare frame fails E4 on exactly three of six records
and renders exactly three of six in the table (§12.8). The classification below is therefore stated
for the frame **through `test_data.run()`**; §12.8 is the one test that deliberately uses the bare
frame, and it uses it as a fixture for the failure, not as a source of expected counts.

Its classification, verified rather than designed (§18):

```
  HAND-1  HUG     treated  flag 0  no reason        → eligible        revealed fact
  HAND-2  CHUV    treated  flag 0  no reason        → eligible        revealed fact
  HAND-3  Lugano  control  flag 0  "Anticoagul…"    → eligible        documented
  HAND-4  USZ     control  flag 0  "Clinician d…"   → eligible        documented
  HAND-5  HUG     treated  flag 0  no reason        → eligible        revealed fact
  HAND-6  Lugano  control  flag 0  no reason        → indeterminate   the [§3] blank

  5 eligible, 1 indeterminate, 0 ineligible — so every ineligible case below is made by
  corrupt("ivt_contraindicated", 1, where=…), which is also the only way to reach E3.
```

**One construction note that is not optional.** Four of the five assertions guard against inputs that
Stage 2 refuses to produce, so the frames that reach them cannot be built by `test_data.run()` and
must be corrupted *after* it:

```
  E1  ivt_contraindicated missing    A4b raises first
  E2  ivt_contraindicated = 2        A5  raises first
  E5  ivt missing                    A4  raises first
  E5  ivt = 2                        A5  raises first

  E3  treated and flagged            reachable through run() — corrupt() alone suffices
  E4  centre outside CENTER_ORDER    reachable, and the BARE hand_frame() already is one
```

§12.3 and §12.4 therefore corrupt the frame after Stage 2 and drive `_assert_classifier_inputs` and
the local variants directly. Verified — this is not a precaution. E4 is the one assertion whose
positive case falls out of an existing fixture at no cost, which is the second reason §12.8 keeps the
bare frame rather than building a new one.

1. **The four cases of the rule.** On `hand_frame()` through Stage 2 and `classify`: the two
   documented controls are `eligible`, the undocumented control is `indeterminate`, all three treated
   records are `eligible`. Asserted by identifier, never by row position. And with
   `corrupt("ivt_contraindicated", 1, where="HAND-3")` the documented control becomes `ineligible` —
   the flag overriding a recorded reason, which is the mask order of §4.4.
2. **A blank is never read as "no contraindication"** — roadmap Stage 4's acceptance criterion.
   HAND-6 has flag 0 and no reason and is `indeterminate`, never `eligible`. The complement is
   asserted in the same test: giving HAND-6 a reason moves it to `eligible`, so the test would fail
   if the classifier ignored the reason column entirely and defaulted everything to `indeterminate`.
   Both halves, because either alone passes against a different broken classifier.
3. **E2, E3, E4 and E5 fire, naming their cases.** `corrupt("ivt_contraindicated", 1,
   where="HAND-1")` — a treated record — raises `SchemaError` naming that identifier.
   `hand_frame(ivt_contraindicated=1)` raises naming **all three** treated identifiers, which is what
   proves the message enumerates rather than reports a count. A flagged *control* does not raise. A
   flag of `2` raises E2. A centre of `"Bern"` raises E4, and so does a **missing** centre — `isin`
   returns `False` for `<NA>`, so `~isin` catches it, and the test asserts that cell of §3.1's third
   fact rather than assuming it. An exposure of `2` raises E5, and so does a missing one. A frame
   tripping several at once reports all of them in one message, which is what the collected raise is
   for; assert the identifiers of two different checks in a single `SchemaError`.

   Two of these are unreachable through `test_data.run()` and must be corrupted after Stage 2, for
   the reason the construction note above gives: A5 raises on a flag of `2` and A4 raises on a
   missing exposure, both before `_classify` is reached.
4. **E1 fires, and the fabrication square.** A frame whose flag is `<NA>` on HAND-6 — corrupted after
   Stage 2, per the note above — raises `SchemaError` naming it. Then, against locally reimplemented
   two-line variants of the chain (never by monkey-patching `eligibility.py`, per Stage 3 §12.3):

   ```
     assertion + fills  (shipped)     → raises                         correct
     no assertion, no fillna          → "ineligible"   FABRICATED, and Stage 5 deletes the patient
     no assertion, with fillna        → "indeterminate" FABRICATED, but retained
   ```

   Both cells are pinned, because they fail differently and the shipped code guards both ways. The
   first is the one Definition of done 5 is watched against.
5. **The classifier is total, and its values are the declared ones.** `eligibility.isna().sum() == 0`
   and `set(eligibility) ⊆ set(ELIGIBILITY_ORDER)` on every frame in this file, including the
   corrupted ones that classify. The dtype is `string` and not `Categorical` (§4.3), asserted the way
   `test_derive.py` asserts it of `onset_type`.
6. **The labels come from the configuration.** An AST constant scan over `eligibility.py` for any
   string literal equal to a member of `ELIGIBILITY_ORDER` — Stage 1 9.4's technique, with a
   companion `test_the_label_scan_actually_fires` against a synthetic snippet that does contain one.
   **And a second scan, for the subtler mistake**: no `Subscript` of `C.ELIGIBILITY_ORDER` by an
   integer constant anywhere in the module, so the rule cannot be written against a display order
   (§4.3). Its companion test fires it too.
7. **`retained` is the declared set.** It equals `!= INELIGIBLE` on every frame here; it equals
   `isin(ELIGIBILITY_RETAINED)` under a `monkeypatch`ed registry whose retained tuple is
   `(ELIGIBLE,)` — where the two stop agreeing, which is what proves it is registry-driven. It is
   `bool` and total. It raises `KeyError` on a frame `classify` has not run on, and the test says so
   rather than leaving it to the docstring.
8. **The cross-tabulation.** Its header is `("centre", "arm", *ELIGIBILITY_ORDER, "n")`. It has
   `len(CENTER_ORDER) * len(TREATMENT_LABELS) + 1` rows plus the header — 9 + 1 — on **every** frame,
   including the two-record fixture, so an absent centre × arm combination renders as zeros rather
   than vanishing. Every row's class counts sum to its own `n`; every column's cell counts sum to the
   `all / both` row; that row's `n` is `len(df)`.

   **Asserted on `hand_frame()` put through `test_data.run()`** — never on the bare frame — where the
   totals are `5 / 1 / 0 / 6` and the nine rows are `HUG` 0/2, `CHUV` 0/1, `Lugano` 2/0, `USZ` 1/0
   by arm. The qualification is load-bearing and is verified rather than assumed (§18): `hand_frame()`
   carries **raw** centre codes, `CENTER_CODES = tuple(config.CENTER_RECODE)`, and Stage 2's
   `centre_recode` is what turns them into `CENTER_ORDER`'s labels. On the bare frame two of the four
   codes happen to spell their own labels and two do not, so the table would render three of six
   records and the column reconciliation would fail.

   **And the bare frame is kept, as a positive test of §7.1's guard.** `_crosstab` on `hand_frame()`
   raises `SchemaError` saying it renders 3 of 6, and `_assert_classifier_inputs` on it raises E4
   naming three of the six identifiers. Both numbers are asserted, because they differ — that
   difference is the evidence that E4 and the reconciliation are two guards and not one written
   twice, and it is the reason §7.1 keeps both.
9. **The audit inventory.** Exactly one new entry, `("derivation", "eligibility")`, with
   `n == len(df)` and the §9.2 table shape. **Its position is asserted against the driver below, not
   against a remembered total**: capture `len(audit.entries)` immediately before the `classify` call
   and assert the entry is the single one appended after that index. On the fixture that index is 10
   — six from `load()` and four from `derive()` — but the test must not name 10, because the count
   differs on the workbook (§11) and `derive_cohort`'s two entries are Stage 5's and arrive *after*
   this one. It renders under the existing
   `## Derivations` heading, and `data.KINDS` is unchanged — asserted directly, so a future
   implementer cannot quietly add a kind. Reproduction stays byte-identical within a process and
   across interpreters launched with `PYTHONHASHSEED=0` and `=1`. The driver is written out rather
   than sketched, for Stage 3 §12.12's reason:

   ```python
   script = ("import sys, data, derive, eligibility\n"
             "df, audit = data.load(data.FIXTURE)\n"
             "df = derive.derive(df, audit)\n"
             "df = eligibility.classify(df, audit)\n"
             "sys.stdout.write(audit.to_markdown())\n")
   ```

   `FIXTURE` and not `WORKBOOK`, so the test runs on a checkout with no `data/`. The fixture reaches
   the entry: 2 records, both `eligible` — one by revealed fact, one by a recorded reason — so it
   exercises two of the four rule branches and the full 9-row table with seven empty rows, which is
   precisely the shape §12.8 exists to protect.
10. **No bare `assert` in `eligibility.py`.** An AST scan for `ast.Assert`, as Stage 2 §12.10, with
    the companion test that proves the scan fires. And no raw workbook header appears in the module
    or in this test file — Stage 1 9.4 covers both automatically, and §12.10 asserts neither file has
    been added to `EXEMPT_FROM_RAW_NAME_SCAN`.
11. **Stage 4 does not depend on Stage 3.** `classify` runs on the Stage 2 frame directly and
    produces the identical column, asserted against the `derive`d frame's. This is what Stages 12 and
    13 rely on, and it is cheap to assert and expensive to rediscover.
12. **Stage 1 amendments**, in `test_config.py` rather than here: `ELIGIBILITY_ORDER ==
    (ELIGIBLE, INDETERMINATE, INELIGIBLE)` and has no duplicates; `ELIGIBILITY_RETAINED` is a subset
    of it and its complement is exactly `(INELIGIBLE,)`; `ELIGIBILITY` is disjoint from
    `ANALYSIS_NAMES` and from `DERIVED_NAMES`; **`ELIGIBILITY` appears in no covariate list** —
    `PS_COVARIATES_FULL`, `BALANCE_ONLY`, `STANDARDISATION_COVARIATES`, and
    `outcome_model_covariates(o)` for every registered outcome (§8); `ELIGIBILITY_RETAINED` is in
    9.8's tuple list; `_RESOLVABLE` contains `ELIGIBILITY`.
13. **Structural facts from the workbook.** Data-gated, against v7 and the §18 record:
    - the flag is 0/1 on all 126 records, never missing, and **never 1 on a treated patient** — the
      two assertions moved here from `test_data.py` §12.14 (§10);
    - the classification is **107 eligible, 19 ineligible** [DECISION 1a; was 64 / 43 / 19 under
      DECISION 1, and 107 = 64 + 43 — see §21.5];
    - the cross-tab reproduces §7.1's nine rows exactly, including the all-zero `USZ / bridging` row;
    - all **19** ineligible records carry a recorded reason, and **no treated patient** carries one;
    - and the one that pins §4.4: dropping the revealed-fact mask — a locally reimplemented chain,
      not a patch to the module — moves exactly **39** records, every one of them from `eligible` to
      `indeterminate`.

    These run against `test_eligibility.py`'s own module-scoped `workbook` fixture (§10) — `load` →
    `derive` → `classify`, `DATA_GATED` — not against `test_data.py`'s, which is module-scoped, does
    not cross files, and produces an unclassified frame in any case.

    The two E-assertion premises are asserted here as well, so that a workbook update re-checks what
    §4.2 says is unreachable rather than leaving it as prose: the exposure is 0/1 on all 126 records
    and never missing (E5's premise), and every `center` value is in `CENTER_ORDER` (E4's premise,
    and the reason §7.1's cells total 126).
14. **Stage 4 adds one column and changes nothing else** — §0.2's promise, asserted rather than
    rested on the `copy()` call. On the hand frame through Stage 2, and again on the fixture through
    `derive`:

    ```python
    before = df.copy()
    out = eligibility.classify(df, audit)

    assert list(out.columns) == [*df.columns, C.ELIGIBILITY]   # appended, and appended last
    assert len(out) == len(df)                                 # no row dropped
    assert out.drop(columns=[C.ELIGIBILITY]).equals(before)    # no value edited
    assert df.equals(before)                                   # the caller's frame is untouched
    ```

    The fourth line is the one that fails if `df.copy()` is ever moved below the assignment or turned
    into a slice, and nothing else in §12 would notice. The third is the one that fails if a future
    edit "tidies" a dtype on the way past. Both are one line each and the promise they check is the
    first sentence of this specification's scope.

### Coverage map

```
eligibility.py                                     test_eligibility.py
├── _assert_classifier_inputs(df)
│   ├── E1 flag <NA>          → SchemaError        ├── 12.4  corrupted AFTER Stage 2; A4b
│   │                                              │         raises first otherwise
│   ├── E2 flag outside {0,1} → SchemaError        ├── 12.3  unreachable through load(): A5
│   ├── E3 treated & flag = 1 → SchemaError        ├── 12.3  names ALL offending cases
│   ├── E4 centre ∉ CENTER_ORDER → SchemaError     ├── 12.3  incl. a MISSING centre
│   ├── E5 ivt missing or ∉ {0,1} → SchemaError    ├── 12.3  unreachable through load(): A4, A5
│   ├── several at once → one message              ├── 12.3  the collected raise
│   └── a flagged control     → no raise           └── 12.3
│
├── _classify(df)
│   ├── treated               → eligible           ├── 12.1  HAND-1/2/5
│   ├── flag 1                → ineligible         ├── 12.1  corrupt(HAND-3)
│   ├── flag 0 + reason       → eligible           ├── 12.1  HAND-3/4
│   ├── flag 0, no reason     → indeterminate      ├── 12.2  HAND-6 — THE criterion
│   ├── flag 0 → reason added → eligible           ├── 12.2  the complement half
│   ├── missing flag, no fillna → "ineligible"     ├── 12.4  the square, cell 2; DoD 5
│   ├── missing flag, fillna    → "indeterminate"  ├── 12.4  the square, cell 3
│   ├── total: never <NA>                          ├── 12.5
│   ├── string, not Categorical                    ├── 12.5
│   ├── labels from config                         ├── 12.6  AST literal scan + fires
│   └── never ELIGIBILITY_ORDER[i]                 └── 12.6  AST subscript scan + fires
│
├── _crosstab(df)                                            reads df[C.ELIGIBILITY]
│   ├── header from ELIGIBILITY_ORDER              ├── 12.8
│   ├── every centre x arm cell rendered           ├── 12.8  9 rows on a 2-record frame
│   ├── rows reconcile to their own n              ├── 12.8
│   ├── cells sum to len(df), else SchemaError     ├── 12.8  the RAW hand frame: 3 of 6
│   ├── columns reconcile to all/both              ├── 12.8
│   └── arm order from sorted(TREATMENT_LABELS)    └── 12.8
│
├── classify(df, audit)                                      written out in §9.3
│   ├── asserts BEFORE it copies                   ├── 12.14 caller's frame untouched
│   ├── column appended last, nothing else moved   ├── 12.14 §0.2's promise, asserted
│   ├── no row dropped, no value edited            ├── 12.14
│   ├── counts.get(cls, 0) — no ineligible on the  ├── 12.9  the fixture is such a frame
│   │   fixture and it must not raise              │
│   ├── one derivation entry, n = len(df)          ├── 12.9
│   ├── entry position vs len(entries) before      ├── 12.9  never a remembered total
│   ├── renders under ## Derivations               ├── 12.9  KINDS unchanged, asserted
│   ├── byte-identical, same process               ├── 12.9
│   ├── byte-identical, PYTHONHASHSEED 0 vs 1      ├── 12.9  driver written out
│   └── runs on a Stage 2 frame                    └── 12.11 what Stages 12/13 rely on
│
├── retained(df)
│   ├── != ineligible                              ├── 12.7
│   ├── moves under a patched registry             ├── 12.7  the registry-driven proof
│   ├── bool and total                             ├── 12.7
│   └── KeyError before classify                   └── 12.7  stated, not left to the docstring
│
└── no bare assert / no raw header / not exempt    └── 12.10 + both scans fire

config.py amendments                               test_config.py
├── ELIGIBILITY_ORDER == the three, no duplicates  ├── 12.12
├── RETAINED ⊂ ORDER, complement == (INELIGIBLE,)  ├── 12.12  the [§3] partition
├── ELIGIBILITY ∉ ANALYSIS_NAMES ∪ DERIVED_NAMES   ├── 12.12
├── ELIGIBILITY in no covariate list               ├── 12.12  [§14a]: not a covariate  [§8]
├── RETAINED in 9.8's immutability list            ├── 12.12  the list is explicit
└── _RESOLVABLE contains ELIGIBILITY               └── 12.12

Data-gated (skipif not DATA_XLSX.exists()):        └── 12.13  107/19 [1a], the 9-row table,
  own module-scoped `workbook` fixture                       19 reasons, 0 treated flagged,
  (load → derive → classify)                                 the 39-record mask, and E4's
                                                             and E5's premises on v7
Every branch above has a test. No branch is untested.
```

### Failure modes

| Failure | Test | Error handling | Visible or silent |
|---|---|---|---|
| A blank read as "no contraindication" | 12.2 | the default is `indeterminate`; nothing promotes a blank | visible |
| A treated patient flagged contraindicated, absorbed by the revealed-fact rule | 12.3 | E3, naming the cases | visible |
| A missing flag classified `ineligible` and deleted at Stage 5 | 12.4 | E1 raises; and `.fillna(False)` guards it a second, independent way | visible |
| A missing flag classified `indeterminate` and silently retained | 12.4 | E1. The fills alone would leave this one | visible |
| A third value in the flag falling through to `indeterminate` | 12.3 | E2 | visible |
| A whitespace-only reason read as "a reason was recorded" | — | **Stage 2's A9**, which Stage 4 depends on and deliberately does not re-check (§4.5) | visible, upstream |
| `== eligible` used where `!= ineligible` was meant — 43 records | 12.7, 12.12 | one `retained` predicate over one declared tuple (§5.1) | visible |
| The revealed-fact mask deleted — 39 records misclassified | 12.13 | the count is asserted against the workbook | visible |
| The rule rewritten by reordering a display tuple | 12.6 | AST subscript scan; the classes are named constants | visible |
| A label typed in `eligibility.py` | 12.6 | AST constant scan | visible |
| `eligibility` becoming a covariate | 12.12 | membership assertion over every covariate list — 9.3's resolution check would *not* catch it, which is why the column stays out of `DERIVED_NAMES` (§8) | visible |
| A centre outside `CENTER_ORDER` vanishing from the cross-tab | 12.3, 12.8 | **E4** at run time, naming the cases; and `_crosstab` raises when its cells do not sum to `len(df)` | visible |
| A patient falling in no cell for some *other* reason — a shortened `CENTER_ORDER`, an arm coding `TREATMENT_LABELS` stopped covering | 12.8 | the `_crosstab` reconciliation, which E4 alone would not catch | visible |
| A treated patient whose exposure was never recorded, classified from the flag alone | 12.3 | **E5**. Both readers of the exposure treat absence as "not treated", so nothing else would notice | visible |
| The tabulated column and the returned column disagreeing | — | **structurally impossible**: `_crosstab(df)` reads `df[C.ELIGIBILITY]` after assignment (§7.1) | n/a |
| `classify` mutating or reshaping the caller's frame | 12.14 | the four assertions of §12.14; `df.copy()` before the write | visible |
| A frame with no ineligible patient raising inside the audit `detail` | 12.9 | `counts.get(cls, 0)`, never `counts[cls]` (§9.3) — the fixture is such a frame | visible |
| A centre × arm combination absent from the data vanishing from the table | 12.8 | every cell rendered; `pd.crosstab` is not used (§7.1) | visible |
| An audit entry dropped, renamed or reordered | 12.9 | inventory comparison against §9.2 | visible |
| A new audit `kind` added without a spec | 12.9 | `data.KINDS` asserted unchanged | visible |
| The log made non-reproducible across processes | 12.9 | two-seed render comparison | visible |
| A check written as `assert`, run under `-O` | 12.10 | AST scan at test time | visible |
| `classify` silently requiring Stage 3 | 12.11 | asserted to run on a Stage 2 frame | visible |
| An indeterminate patient assigned a regime under [§14b] with no decision behind it | none | **deferred** — Stage 13's, recorded in §14 |

One row is deferred and one is upstream; every other failure is visible at the point it occurs.

## 13. Known gaps carried forward

- **The indeterminate group is centre-driven, not patient-driven, and no analysis brackets it.** §18
  records all 43 indeterminate records as controls at CHUV and Lugano — the two centres that never
  collected the reason field — which is 80% of the primary cohort's control arm. [§3] already states
  the consequence: retaining them "does assume they were eligible, and no analysis brackets that
  assumption — state it as a limitation". Stage 4 makes the group *countable*, in §7's table and in
  Stage 5's flow line. It does not make it smaller, and nothing in this stage should be read as
  having addressed it.
- **A patient's eligibility is a property of their centre's documentation practice, and the primary
  cohort keeps centre as a covariate.** So the indeterminate indicator and `center` are close to
  collinear by construction. This is not a Stage 4 problem — `eligibility` never enters a model
  (§8) — but it is why a future amendment proposing an indeterminate-versus-documented sensitivity
  analysis would be adjusting for something very like centre. Recorded so that proposal arrives with
  the objection already in hand.
- **The `[§4]` citation.** Corrected in `config.py` (§10); left standing in `implementation_roadmap.md`
  Stage 0 and Stage 4 and in `../out/stage0_data_inventory.md`, both of which are dated records of
  decisions taken. A reader following `[§4]` from any of those reaches the Exposure section and finds
  no classification rule; §4.1 is the pointer that resolves it.
- **E1 through E5 are all unreachable on v7**, so every one of them is tested only against hand-built
  frames. That is the same posture Stage 3 takes toward `_assert_onset_flags` and is recorded for the
  same reason: an implementation that omitted an assertion would be green on this workbook. The same
  holds for `_crosstab`'s reconciliation, which on v7 can only ever compute 126 = 126. §12.13 asserts
  each of their premises against the workbook, so a corrected workbook that broke one would say so.
- **`test_derive.py`'s Stage 3 audit-inventory test is a trap that Stage 4 must not spring and Stage 5
  will.** Its assertion is a *tail slice* — `audit.entries[-len(_STAGE_3_INVENTORY):]` at
  `test_derive.py:714` — which is correct only because its helper `_full_pipeline()` is
  `load → derive → derive_cohort` and never calls `classify`. The real pipeline puts this stage's
  entry between the two, so the slice will span it (§11). Stage 4 leaves the file untouched because
  T6's verification depends on that; the repair is recorded in `TODOS.md` and belongs to the commit
  that first builds the real order. Recorded here because the reason it is safe today is an accident
  of a test helper, and an accident that nobody wrote down is a defect waiting for its turn.
- **`retained` is not enforced.** Nothing prevents Stage 5, 12 or 13 from writing
  `df[df["eligibility"] == "eligible"]` by hand and quietly dropping 43 patients. §5.1 makes the
  right thing available and §12.12 keeps the constant honest; the invariant that would catch the
  wrong thing is Stage 5's to write, against its own cohort-flow counts.

## 14. What Stage 4 deliberately does not decide

Recorded so a later stage does not look here for an answer that was never placed here.

- **Which patients and centres are dropped.** Stage 5. Stage 4 classifies and counts; it removes
  nothing, and `retained` is a predicate rather than a filter for exactly that reason.
- **Which way an indeterminate patient goes under [§14b]'s active regime.** Stage 13. [§14b] says the
  active regime is "bridge if eligible, direct thrombectomy if contraindicated" — a binary decision
  over a three-class column, and the third class has no assignment. Both readings are defensible
  (bridge them, since [§3] retains them as assumed-eligible; or give them direct EVT, since nothing
  documents that IVT was available) and they answer different policy questions. Stage 4 supplies the
  three-level column precisely so that Stage 13 can decide in writing rather than by whichever
  comparison someone types.
- **Whether the cohort-flow table gets its own audit `kind`.** Stage 5 (§9.1).
- **How the limitation in §13's first bullet is worded in the manuscript.** Stage 14 [§16].

## 15. NOT in scope for Stage 4

| Considered | Why deferred |
|---|---|
| Classifying the free text of `Contraindications_to_IVT` | DECISION 1. It would also require naming a raw header, which Stage 1 §7 forbids, and would reopen the typo and annotation problems §4.1 records as retired |
| Normalising or stripping the reason column | Stage 2 §6 declines it and A9 already raises on the one case that matters. Stripping here would be reading the text (§4.5) |
| A three-category sensitivity analysis excluding the indeterminate group | Not in [§13], and [§13]'s sensitivity suite is deferred by decision (roadmap Stage 11). It is also close to a centre contrast (§13) |
| Making `eligibility` a covariate or an interaction term | [§14a] rules it out explicitly; §8 asserts against it |
| Making `eligibility` a `Categorical` | Stage 5 removes a whole class; Stage 3 §4.4's argument applies unchanged (§4.3) |
| Dropping ineligible patients here | Stage 5, which owns both [§3] restrictions and the flow table that reports them together |
| A CLI entry point | Stage 14 is the single entry point [§16]. `eligibility.py` is a library module |
| Caching the classification | 126 rows, one vectorised pass. A cache is a second source of truth |

## 16. What already exists, and what to lift

`stage0_data_inventory.py` is in the repository and is the only prior art; it predates the
configuration and reads raw headers by design, which is why it is one of the three files exempt from
the 9.4 scan.

| From Stage 0 | Status |
|---|---|
| `stage0_data_inventory.py:106`'s `eligibility(df)` — the four cases in this order, treated last | **lift the rule and the order**, and its docstring's "Order matters" sentence, which §4.4 expands |
| its `pd.Series("indeterminate", …, dtype=object)` | **do not lift.** `string`, per §4.3 and Stage 3 §4.4's precedent for `onset_type` |
| its `out[cond] = value` assignment form | **do not lift**, and know why: `[]` treats `<NA>` as False where `.mask` treats it as True (§3.1), so the two constructs fail *differently* on a missing flag. `.mask` with explicit fills is chosen because the fill is then visible at every site rather than implied by the construct |
| its raw header constants `ABS_FLAG`, `REASON`, `TREATMENT` | **do not lift.** Analysis names, from `config.py` |
| its assertion that no treated patient is flagged, with "resolve before proceeding" | **lift the check and the posture**; not the statement form (Stage 2 §8.3). E3 raises `SchemaError` and names the cases |
| its cross-tabulation by centre and arm | lift the **shape**, not `pd.crosstab`: §7.1 renders every cell, including the absent combinations `pd.crosstab` drops, and the all-zero `USZ / bridging` row is the [§3] restriction-1 finding |
| its cross-tabulation's *silence* on records that fall in no cell | **do not lift.** Stage 0 was an inventory and could afford a table that quietly covered a subset; this one is read as the [§3] restriction-1 finding and as Stage 5's flow source, so §7.1 raises rather than under-reporting. The `all / both` row cannot substitute: it is computed from the frame, so a vanished patient leaves it correct |
| its assertion posture generally — check the flag, resolve with the owner, do not choose | **lift and extend.** Stage 0 asserted one property of one column; §4.2 asserts five over four, because `classify` is called on frames Stage 0's script never was (§6.1) |

## 17. Implementation tasks

Ordered. Each independently verifiable. T1 touches Stage 1 files; T2–T5 are Stage 4's own; T6 is the
amendment sweep. The second lane depends on the first, so there is one lane in practice.

- [ ] **T1 (P1)** — `config.py`: the §5.2 block, in the file position it already occupies;
      `ELIGIBILITY_ORDER` becomes computed. `test_config.py`: the §12.12 assertions, **including**
      `ELIGIBILITY_RETAINED` added to 9.8's list and `ELIGIBILITY` added to `_RESOLVABLE`. Verify:
      acceptance 12.12 — and see the partition assertion fail first, by temporarily adding a fourth
      class to `ELIGIBILITY_RETAINED`, since nothing covered these constants before this task.
- [ ] **T2 (P1)** — `eligibility.py`: `_assert_classifier_inputs` with E1, E2, E3, E4, E5 and the
      collected raise. No `.fillna` on E4's or E5's `~isin` mask (§3.1, §4.2). Verify: acceptance
      12.3, 12.4's raise half.
- [ ] **T3 (P1)** — `_classify` and `retained`. Verify: acceptance 12.1, 12.2, 12.4's square, 12.5,
      12.6, 12.7.
- [ ] **T4 (P1)** — `_crosstab(df)` including its reconciliation raise, `classify` as written in
      §9.3, and the §9.2 entry. Verify: acceptance 12.8, 12.9, 12.11, 12.14.
- [ ] **T5 (P2)** — `test_eligibility.py` remainder: its own module-scoped `workbook` fixture, the
      AST scans, the reproduction driver, and the data-gated §18 facts including the 39-record
      assertion. Verify: `uv run pytest -v` green both with and without `data/`.
- [ ] **T6 (P1)** — the §10 sweep: move `test_data.py` §12.14's two assertions here **and give
      `test_eligibility.py` its own `workbook` fixture**, since the one they take does not cross
      files; the two `derive.py` docstring counts and the two in `specs/stage3_derived_variables.md`;
      the four `[§4]` citations at `config.py` 203, 412, 416 and 464 — two by hand, two with §5.2's
      block rewrite; and the `.loc`/`<NA>` sentence in `_assert_onset_flags`. Verify: `uv run pytest
      -v` green with **no test edited except the two moved**, which is where §10's claim that nothing
      else depends on those counts is checked rather than believed.

      **The roadmap is not in this sweep.** Its `**Spec:**` line and its rewritten **Accept when**
      landed with this document on 2026-08-12, which is what §1 and §19 always said should happen and
      what T6 previously contradicted (§20 defect 10). Nothing about the roadmap is left to do at
      implementation time; if a future edit needs one, it is a new decision and not this task.

### Diagrams that belong in the code, not only here

Two, and no more — a diagram nobody maintains is worse than none, because it is believed. Keeping them
true is part of any change that touches them, in the same commit.

- **`eligibility.py`'s module docstring** — §0's flow, and specifically the two arrows out: the
  restricted path into Stage 5, and the *unrestricted* frame into Stage 13. The second is the reason
  this module exists separately (§0.1) and is not visible from its functions.
- **Above `_classify`** — §4.1's four-case table, plus the `<NA>` row that §3.1 is about: a missing
  flag reaches no case, and without the fills it takes the last one. That is the line an implementer
  "tidying the masks" deletes.

One further comment, not a diagram, and it is the shortest of the three: **above E4 and E5, two lines
saying why they carry no `.fillna(False)` when every mask in `_classify` does** (§3.1's third fact).
The uniformity an implementer would impose here is the one that breaks E5.

### Definition of done

Stage 4 is complete when all of the following hold, and not before.

1. `uv run pytest -v` is green **with** `data/` present — every test, including the data-gated ones.
2. `uv run pytest -v` is green **with `data/` temporarily renamed** — the data-gated tests skip and
   nothing else fails or errors at collection.
3. `test_config.py`, `test_data.py` and `test_derive.py` are still green, and neither `eligibility.py`
   nor `test_eligibility.py` has been added to `EXEMPT_FROM_RAW_NAME_SCAN`.
4. The pipeline's audit log has been produced through all three stages, against the **workbook**, and
   reproduces byte-identically across hash seeds:

   ```bash
   cd extended_bridging
   S='import sys, data, derive, eligibility
   df, audit = data.load(data.WORKBOOK)
   df = derive.derive(df, audit)
   df = eligibility.classify(df, audit)
   sys.stdout.write(audit.to_markdown())'
   PYTHONHASHSEED=0 uv run python -c "$S" > /tmp/a.md
   PYTHONHASHSEED=1 uv run python -c "$S" > /tmp/b.md
   diff /tmp/a.md /tmp/b.md
   ```

   `diff` reports nothing, and the rendered eligibility table matches §7.1 cell for cell. The log
   names patients and must never be committed.
5. Acceptance 12.4 has been **seen to fail** with `.fillna(False)` removed from `_classify` *and* E1
   removed from the assertion — specifically its `ineligible` cell, which is the silent-exclusion
   failure §3.1 is about. Then seen to pass when both are restored. A green test that has never been
   watched failing is not evidence.
6. Acceptance 12.3 has been **seen to fail** with E3 removed: the treated flagged record classifies
   `eligible` and nothing anywhere says so. And 12.13's 39-record assertion seen to fail with the
   revealed-fact mask removed — it is the one data-gated test that would catch a plausible-looking
   rewrite of §4.4's order.
7. Acceptance 12.8 has been **seen to fail** with the `_crosstab` reconciliation removed: on the bare
   `hand_frame()` the table renders three of six records, every total agrees with itself, and nothing
   raises. Then seen to pass when it is restored. This is the silent-omission failure §7.1 is about
   and the reason the check is at run time rather than only in a test.
8. Acceptance 12.3 has been **seen to fail** with E4 removed, and separately with E5 removed. E4's
   half must be watched on the bare `hand_frame()`, where it and the reconciliation catch *different*
   subsets — if removing E4 changes nothing, the two guards have been collapsed into one and §7.1's
   argument for keeping both no longer holds.
9. A `.fillna(False)` pasted onto E5's mask has been **seen to make a test fail** — 12.3's
   missing-exposure half. This is the one uniformity an implementer is most likely to impose (§3.1),
   and it must not be silently absorbable.
10. Acceptance 12.14 has been **seen to fail** with `df.copy()` moved below the column assignment:
    the caller's frame gains the column. Nothing else in §12 notices, which is the point.
11. Acceptance 12.6's two scans have been **seen to fail** against a pasted label literal and a pasted
    `C.ELIGIBILITY_ORDER[0]` respectively, and 12.10's against a pasted `assert True`. All three pass
    trivially if they match nothing.
12. The audit log has been read end to end by a human, and the three class counts and the nine table
    rows agree with §18 and with `../out/stage0_data_inventory.md`.
13. `grep -n '\[§4\]' extended_bridging/config.py` returns **nothing**, and no `.py` file writes any
    member of `ELIGIBILITY_ORDER` as a literal outside `config.py`.

    The pattern is deliberately the **bracketed** citation and not a bare `§4`: `config.py:325` says
    "Stage 3 §4.1 asserts no record is…", which is a legitimate cross-reference to another spec's
    section and must survive. A bare `§4` grep can never return nothing and so can never be a
    definition-of-done step (§18).

## 18. Verification record

Checked on 2026-08-10, before this spec was finalised, by a read-only probe reporting aggregate counts
only — no outcome by arm, in the posture Stage 0 used. Recorded so a reader can tell which numbers
were verified rather than carried, and so the checks are re-runnable after a workbook update.
Acceptance 12.13 re-checks them in code.

The probe is the §12.9 driver with the classification chain of §4.3 pasted in, run against
`data.WORKBOOK`; it writes nothing.

| Claim | Where used | Verified |
|---|---|---|
| `ivt_contraindicated` is 0/1 on all 126 records and never missing | §4.2, §6.1, §12.13 | yes, run |
| **No treated patient carries the flag** — 0 of 39 | §4.2, §4.4, §12.13 | yes, run |
| Classification over 126 records: **64 eligible, 43 indeterminate, 19 ineligible** — DECISION 1's. Re-run under DECISION 1a as **107 / 19**, the 43 having moved to `eligible` (§21.5) | §7.1, §12.13 | yes, run, both rules |
| The nine-row cross-tab of §7.1, cell for cell, including the all-zero `USZ / bridging` row | §7.1, §12.8, §12.13 | yes, run |
| It reproduces `../out/stage0_data_inventory.md`'s eligibility table exactly | §7.1 | yes, compared |
| A reason is recorded on 44 of 126 records; **all 19** ineligible carry one; **no treated patient** does | §4.4, §12.13 | yes, run |
| Therefore, without the revealed-fact mask, **39** records move — every one from `eligible` to `indeterminate`, none in the other direction | §4.4, §12.13 | yes, run |
| All 43 indeterminate records are controls, at CHUV and Lugano only — under DECISION 1a the same 43 are the control-arm records with no documented reason (§21.4) | §7.2, §13 | yes, run, both rules |
| The 43 are 80% of the primary cohort's 54 retained controls | §5.1, §13 | carried from Stage 0's flow table |
| `classify` on `derive`'s frame takes it from 31 to 32 columns | §10, §11 | yes, run |
| `retained` returns numpy `bool`, total, and `df[retained(df)]` selects 107 of 126 on v7 | §5.1, §12.7 | yes, run |
| Every `center` value on v7 is in `CENTER_ORDER` — E4's premise, and why §7.1's cells total 126 | §4.2, §6.1, §12.13 | yes, run — re-run 2026-08-12 |
| `ivt` is 0/1 on all 126 records and never missing — E5's premise | §4.2, §6.1, §12.13 | yes, run — re-run 2026-08-12 |

**Re-verified on 2026-08-12**, during the engineering review of §20, by re-running the probe: every
number in the two tables above reproduced exactly — 0 missing flags, 0 outside `{0,1}`, 0 treated
flagged, 64/43/19, 44 reasons recorded, 19 of 19 ineligible carrying one, 0 treated carrying one, 39
records moved by the revealed-fact mask all from `indeterminate` to `eligible`, 107 retained, the
nine-row table cell for cell, and 31 → 32 columns.

**Claims added by the §20 review**, all counted or run rather than read off a previous sentence:

| Claim | Where used | Verified |
|---|---|---|
| `load()` appends **6** entries on the fixture and **7** on the workbook; the seventh is `observation` / `onset_to_groin_999`, which `test_data.py`'s `test_the_fixture_triggers_no_observation` pins as fixture-absent | §11, §12.9 | yes, run |
| `derive()` appends **4** and `derive_cohort()` a further **2**; `test_derive.py`'s `_STAGE_3_INVENTORY` is the six and `test_derives_four_entries_precede_derive_cohorts_two` pins the split | §11, §12.9 | yes, run and read |
| `derive_cohort` is called by **Stage 5**, so its two entries land *after* Stage 4's, not before — on the §12.9 driver the eligibility entry is at index 10 | §11, §12.9 | yes, run |
| `test_derive.py`'s `_full_pipeline()` does not call `classify`, so its inventory test is unaffected by this stage | §10 | yes, read |
| `Series.isin` returns `False` for `<NA>` and never `<NA>`; `~isin` therefore selects the missing row. Dtype differs by source and neither carries `<NA>`: numpy `bool` on the `string` `center`, pandas `boolean` on the `Int64` exposure | §3.1, §4.2 | yes, run |
| `hand_frame()` carries **raw** centre codes — `CENTER_CODES = tuple(config.CENTER_RECODE)` is `("1", "Lausanne", "Lugano", "USZ")` — so `_crosstab` on the bare frame renders **3 of 6** records while `all / both` still reads `5 1 0 6`, and E4 names **3 of 6** identifiers because two raw codes spell their own labels | §7.1, §12.8 | yes, run |
| `tuple(CENTER_RECODE.values()) == CENTER_ORDER` today, by construction and by nothing else; Stage 2's A3 tests membership of the former while §7.1's table is built from the latter | §6.1 | yes, run and read |
| `test_data.py`'s `workbook` fixture is module-scoped (its line 628) and both moved tests take it; module-scoped fixtures do not cross files and the repository has no `conftest.py` | §10, §12.13, T6 | yes, read |
| `config.py` carries **four** `[§4]` citations — 203, 412 (the eligibility block's own header), 416, 464 — and `grep -n '§4' config.py` additionally returns line 325's legitimate `Stage 3 §4.1` cross-reference | §4.1, §10, DoD 13 | yes, run |
| `Audit.record`'s signature is `(kind, step, n, detail, case_ids=(), table=None)`, and `AuditEntry.__post_init__` accepts `derivation` with no `case_ids` | §9.2, §9.3 | yes, read |
| `value_counts()` on the fixture's classified column has no `ineligible` key — both its records are `eligible` — so `counts[cls]` would raise there and `counts.get(cls, 0)` is required | §9.3, §12.9 | yes, run |

**pandas semantics, checked by running them** on pandas 2.3.3 as pinned by `uv.lock`. §3.1 is the
declaration; these are the observations behind it. Re-check on any pandas major bump — the first two
decide whether a missing flag is caught or silently excluded.

| Claim | Where used | Verified |
|---|---|---|
| `Series.mask(cond)` with `<NA>` in `cond` treats it as **True** and applies the replacement | §3.1, §4.3, §12.4 | yes, run |
| `df.loc[cond]` and `s[cond] = v` with `<NA>` in `cond` treat it as **False** — the opposite | §3.1, §16, §10 | yes, run |
| The §4.3 chain with the fills removed classifies a missing flag as **`ineligible`**; with the fills, as `indeterminate` | §3.1, §4.3, §12.4, DoD 5 | yes, run |
| `(flag == 0) & reason.notna()` on `Int64` × `string` yields `<NA>` where the flag is missing and a reason is present, `False` where neither is | §3.1 | yes, run |
| `Series.isin` on a `string` column returns numpy `bool`, never `<NA>` | §5.1, §12.7 | yes, run |
| `.sum()` on a `BooleanDtype` mask skips `<NA>` | §7.1 | yes, run |

**Claims checked against the committed repository**, so they are re-checkable with no `data/`:

| Claim | Where used | Verified |
|---|---|---|
| `hand_frame()` gives 5 eligible and 1 indeterminate, with the per-record classification of §12; `corrupt("ivt_contraindicated", 1, where="HAND-1")` raises E3 and `hand_frame(ivt_contraindicated=1)` names all three treated records | §12, §12.1, §12.3 | yes, run |
| A frame with a missing flag **cannot reach `_classify` through `test_data.run()`** — Stage 2's A4b raises first, so §12.4 must corrupt after Stage 2 | §12, §12.4 | yes, run |
| `tests/fixture_schema.xlsx` classifies as 2 eligible — one by revealed fact, one by a recorded reason — and renders the full 9-row table | §12.9 | yes, run |
| Stage 2's A4b asserts `ivt_contraindicated` is never missing, A5 that `BINARY_COLUMNS` are 0/1 **or missing**, and A9 that no reason cell is whitespace-only | §4.5, §6.1 | yes, read |
| Stage 2 §6 states `contraindication_reason` is deliberately not normalised | §4.5 | yes, read |
| `config.INFORMATIVE_ABSENCE` labels the column `not recorded` and `absence_by_column` already renders its per-centre counts | §7.2 | yes, read |
| `test_data.py` §12.14 holds `test_no_treated_patient_is_flagged_as_contraindicated` and `test_the_eligibility_classifier_is_never_missing`, both data-gated, and nothing else there touches eligibility | §10, §12.13 | yes, read |
| `test_config.py` 9.8's tuple list and 9.13's `parametrize` union both name their constants explicitly; `_RESOLVABLE` is `ANALYSIS_NAMES | DERIVED_NAMES | OUTCOMES` | §8, §12.12, T1 | yes, read |
| 9.3's `test_every_ps_covariate_resolves` accepts a covariate that is in `ANALYSIS_NAMES` **or** `DERIVED_NAMES`, so widening `DERIVED_NAMES` would widen what may be a covariate | §8 | yes, read |
| `ELIGIBILITY_ORDER` is read only by 9.8's immutability list, which asserts its type and not its contents — so making it computed moves nothing else | §5.2 | yes, read |
| `test_config.py` 9.4 walks `extended_bridging/**/*.py` minus a pinned three-file exemption set, so `eligibility.py` is scanned from the moment it exists | §3, §12.10, DoD 3 | yes, read |
| `test_derive.py`'s column-count assertion is relative (`len(cohort.columns) == len(out.columns) + 1`) and survives the insertion | §10 | yes, read |
| `derive.py`'s docstring diagram says `32 columns` after `derive_cohort`, and `derive_cohort`'s docstring says "a new frame of 32 columns" — both stale once this stage lands | §10 | yes, read |
| `derive.py`'s `_assert_onset_flags` comment says `.loc` "refuses to mask with `<NA>`", which is false on pandas 2.3.3 | §3.1, §10 | yes, read and run |
| SAP §4 is *Exposure*, one sentence, with no classification rule; `config.py` cites `[§4]` for the classifier at lines 203, 416 and 464, and the roadmap twice | §4.1, §10, §13 | yes, read |
| `stage0_data_inventory.py:106`'s classifier uses this rule order, `dtype=object`, `[]`-assignment and raw header constants | §16 | yes, read |
| `data.KINDS` already contains `derivation` and `_MUST_NAME_CASES` is `{correction, observation}`, so Stage 4 needs no `data.py` edit | §9.1, §12.9 | yes, read |

## 19. What this spec changed elsewhere

Four entries, all landing with this document rather than with the implementation, because a
specification that contradicts the roadmap is worse than no specification.

| # | Change | Files | The question it settles |
|---|---|---|---|
| 1 | One retained predicate: [§3] and [§14a] both retain the indeterminate group, and only `ineligible` is ever dropped | `implementation_roadmap.md` Stage 4's **Accept when**, and §5 here | Whether "all eligible patients" in [§14a] is `== eligible` or `!= ineligible`. It is the latter, by the PI's decision of 2026-08-10, and 43 records turn on it |
| 2 | Roadmap Stage 4 gains its `**Spec:**` line | `implementation_roadmap.md` | Stage 4 was the only stage without one |
| 3 | The eligibility entry reuses the `derivation` kind; no new audit kind lands with this stage | §9.1 here | Whether Stage 4 or Stage 5 owns a cohort-shaped audit heading. Stage 5, if it wants one — inserting it here would spend a `data.py` amendment on a spec not yet written |
| 4 | Roadmap Stage 4's **Accept when** gains the two requirements a reader of the old wording would not have known to check: every declared centre × arm cell rendered whether or not the data fills it, and the retained predicate asserted as `!= ineligible` including under a patched registry | `implementation_roadmap.md` | Whether "a cross-tabulation is produced" permits `pd.crosstab`. It does not (§7.1) — and the row it would drop is the [§3] restriction-1 finding |

Everything else this spec asks for is an amendment to Stage 1, 2 or 3 in service of Stage 4, and §10
is the ledger for those. Two of them are corrections of statements that are false rather than merely
stale — the `[§4]` citations and the `.loc`/`<NA>` sentence — and both are listed there separately
for that reason. None of them changes a statistical decision, a column's fate, or a number.

**What this spec has now had.** The engineering review §19's earlier draft asked for, on 2026-08-12.
It is recorded in §20. It did not have a cross-model second opinion: the Codex pass failed to
authenticate, and the alternative — a second Claude context — was declined as too weak a signal to be
worth the delay. That remains open, and §20 says where a second reader should start.

## 20. Engineering review, 2026-08-12

Ten defects, every one of them reproduced by running code against `tests/fixture_schema.xlsx`, the
hand frame, and v7, or by reading the named file — not by inspection. All ten are folded into the
sections above; this section is the record of what was wrong, so that a reader who knows the earlier
draft can see what moved.

| # | Defect | Where it was | Fixed in |
|---|---|---|---|
| 1 | The audit-entry counts were wrong twice over: Stage 2 was called seven when it is six on the fixture and seven only on the workbook, and Stage 4's entry was said to follow "Stage 3's six" when `derive_cohort`'s two arrive at Stage 5, *after* it | §11, §12.9 | §11's counted block, §12.9's position-against-the-driver rule, six §18 rows |
| 2 | A centre outside `CENTER_ORDER` vanished from the cross-tab at run time with nothing raising. The failure-mode table claimed §12.8's reconciliation covered it; that is a test-time assertion on hand frames, not a runtime guard | §7.1, failure modes | **E4** in §4.2 and the `_crosstab` reconciliation in §7.1 — two independent guards, and §12.8 asserts they catch different subsets |
| 3 | `_crosstab(df, eligibility)` tabulated a Series that need not be the column `classify` returns | §7.1 | `_crosstab(df)`, reading `df[C.ELIGIBILITY]` after assignment |
| 4 | `classify` — the public entry point, and the only function an implementer had to invent — was never written out. `detail`'s placeholders did not say which interpolate config constants | §3, §9.2 | §9.3, written out with the five things an implementer could reasonably get wrong |
| 5 | The classifier read `ivt` in two places and asserted nothing about it, while §6.1 argued at length for re-checking the *other* flag on reasoning that applies to this one verbatim | §4.2, §6.1 | **E5**, and §6.1 rewritten for five checks |
| 6 | The old Definition of done 9's `grep -rn '§4' config.py` could never return nothing — `config.py:325` legitimately says "Stage 3 §4.1" — and the citation inventory named three sites where there are four, missing the eligibility block's own header at 412 | DoD, §10, §18 | the bracketed `\[§4\]` pattern in **DoD 13**, and all four sites named in §4.1, §10 and §18 |
| 7 | §12.8's "asserted on the hand frame" was ambiguous and, read literally, wrong: `hand_frame()` carries raw centre codes, so the column reconciliation fails on it | §12.8 | the frame named, and the bare frame kept as a *positive* test of defect 2's guard |
| 8 | The two tests being moved out of `test_data.py` take a module-scoped `workbook` fixture that does not move with them, so T6's "no test edited except the two moved" could not hold | §10, §17 T6 | `test_eligibility.py` declares its own, which §12.13 needs anyway |
| 9 | §0.2's promise — one column added, no value edited, no row dropped — rested entirely on the `copy()` call and had no acceptance criterion | §0.2, §12 | **§12.14**, four one-line assertions, and DoD 10 |
| 10 | The document contradicted itself about when its own roadmap amendment lands: §1 and §19 said "with this spec", T6 put it in the implementation sweep — and it had landed nowhere. Stages 1, 2 and 3 each carry a `**Spec:**` line in `implementation_roadmap.md`; Stage 4 did not | §1, §19, §17 T6 | the roadmap edit **made on 2026-08-12**: Stage 4 gains its `**Spec:**` line, a retained-predicate paragraph, and an **Accept when** that names the every-cell-rendered and `!= ineligible` requirements. T6 no longer claims it |

**Four consequential edits the nine defects dragged behind them**, listed because each is a place a
reader of the earlier draft would otherwise carry a stale number: §0.1's "three assertions" and §2's
"three boolean reductions" became five; §12's construction note became a six-row table of which
assertions Stage 2 refuses to let a frame reach; §13 gained the `test_derive.py` tail-slice hazard
that Stage 5 will spring; and the Definition of done grew from nine items to thirteen, the four new
ones all of the *seen-to-fail* kind that the earlier list already used for E1, E3 and the AST scans.

**What the review did not change.** §4.4's ordering argument and §7.1's table shape — the two places
the earlier draft named as most likely to hide a plausible-looking simplification — were checked and
stand. The 39-record claim, the 43/64/19 split, the fabrication square and both pandas facts of §3.1
all reproduced exactly. The one substantive addition to §3.1 is a **third** fact, `isin`'s treatment
of `<NA>`, which the new E4 and E5 rest on and which points the opposite way from the two already
there.

**Where a second reader should start.** §9.3's `classify`, which is new and therefore the least
reviewed code in this document; and the E4/E5 fill discipline of §4.2, where the module now contains
two masks that must never carry `.fillna(False)` sitting four lines from three that must.

## 21. Amendment, 2026-08-13 — DECISION 1a: eligibility is two-valued

**The decision.** `ivt_contraindicated` is the classifier and `Contraindications_to_IVT` is read by
nothing — neither its text, which DECISION 1 already excluded, nor its presence, which DECISION 1 used
as one bit. Flag = 1 is `ineligible`; flag = 0 is `eligible`, including every patient whose
contraindication reason was never documented. The `indeterminate` class is withdrawn.

Recorded in `../out/stage0_data_inventory.md` as DECISION 1a and in [§3] as the amendment of the same
date. This section is the Stage 4 ledger and **supersedes every section above wherever the two
disagree**.

### 21.1 What did not change

**The population, and every count that describes it.** Verified by running both rules and comparing
identifier for identifier: the retained set is the same **107** records and the primary cohort the same
**93** patients — 39 bridging, 54 thrombectomy alone, HUG 30/11, CHUV 7/14, Lugano 2/29. The 43
formerly `indeterminate` records were already retained by [§3], so the amendment relabels them and
moves none. §7.1's `n` column and its all-zero `USZ / bridging` row are unchanged.

Also unchanged: all five preconditions E1–E5 and their messages; `_crosstab`'s every-cell-rendered rule
and its reconciliation raise; the `derivation` audit kind and the single entry; the `df.copy()` promise
of §0.2; the pandas facts of §3.1; and §8's rule that `eligibility` is not a derived name and reaches
no covariate list.

### 21.2 What changed, section by section

| Section | Change |
|---|---|
| §4.1 | Four cases become two. The rule blocks are rewritten in place with the superseded one kept beneath |
| §4.3 | One mask, not three. Rewritten in place. Four of its five properties survive verbatim |
| §4.4 | **Void.** The ordering argument, the precedence discussion and the 39-record claim all concern masks that no longer exist. There is one mask, so there is no order |
| §4.5 | **Void**, and it is the section most worth reading anyway: its last paragraph made A9 load-bearing for Stage 4, and that dependency is gone |
| §5.1 | The retained predicate survives; its *argument* does not. With two classes `== ELIGIBLE` and `!= INELIGIBLE` coincide, so the 43-record gap is closed. `retained()` stays because `ELIGIBILITY_RETAINED` is the seam a third class would return through |
| §5.2 | `INDETERMINATE` removed; `ELIGIBILITY_ORDER` is two-valued; `ELIGIBILITY_RETAINED` is `(ELIGIBLE,)` |
| §6.1 | Unchanged, except that **E3 moves from belt-and-braces to load-bearing** — see §21.3 |
| §7.1 | One fewer class column; the rendering is re-run and both versions are shown |
| §7.2 | Unchanged, and now more pointed: the per-centre completeness table in Stage 2's `absence_by_column` is the **only** place the reason column's collection pattern is recorded, since no class encodes it any more |
| §12.1 | "The four cases" becomes two. `test_the_flag_overrides_a_recorded_reason` survives unchanged — the flag still decides — but for a different reason: there is no competing rule left to override |
| §12.2 | Rewritten. "A blank is never read as no contraindication" is no longer the criterion; the criterion is that **blanking the reason column on every record leaves the classification unmoved**, which is the assertion that catches a partial revert |
| §12.4 | The fabrication square survives with one cell relabelled: without the fill a missing flag is still `ineligible` and still deleted at Stage 5; with the fill it is now `eligible` rather than `indeterminate`. Two cells, still failing differently |
| §12.6 | Both AST scans survive unchanged |
| §12.7 | The patched-registry test had to be inverted: it patched `ELIGIBILITY_RETAINED` to `(ELIGIBLE,)`, which is now the shipped value. It patches to `(INELIGIBLE,)` instead |
| §12.13 | 64/43/19 becomes 107/19. The 39-record revealed-fact test becomes its converse — that restoring the deleted mask changes nothing. The 43-record test becomes a count of retained **control-arm** records with no documented reason |
| §18 | Re-run; §21.5 records the new rows |

### 21.3 E3 is now load-bearing, and that is the one thing to carry forward

Under DECISION 1 the revealed-fact mask forced every treated patient to `eligible`, so a
treated-and-flagged record would have been classified `eligible` and E3's job was to stop that
contradiction being *absorbed silently* — a labelling error in a patient who stayed in the cohort.

Under DECISION 1a there is no such mask. A treated patient carrying the flag classifies **`ineligible`**
and Stage 5 deletes them. So E3 has stopped being a guard against a mislabelled retained patient and
become the guard against a **silently removed** one. It is the only check standing there, and it is
what makes deleting the mask safe rather than merely equivalent. §12.13 asserts the redundancy against
the workbook; E3 asserts it against every frame.

The corresponding sentence in §4.2 — "E3 is the assertion this stage exists to add" — was true when
written and is more true now.

### 21.4 What is quietly worse, and where it is recorded

Under DECISION 1 the 43 undocumented controls carried a label. Every table in the pipeline showed
them, and [§3] required a limitation about them. They now carry the same label as a patient who was
assessed and documented, and **nothing in the frame distinguishes them**.

The assumption did not go away; it moved. Eligibility now rests on `flag = 0` meaning the same thing at
the two centres that never collected a reason as at the two that collected one for every control. The
flag is 0/1 and never missing on all 126 records, which is equally consistent with its having been
assessed for everyone and with 0 being an uncompleted default; nothing in the data separates them. It
governs **43 of the 54** controls in the primary cohort.

Three things count the group, and they are the whole of the mitigation:

- the audit `detail` this stage writes, which names the arm-restricted count;
- Stage 5's cohort-flow table, whose `of which EVT alone, no reason on file` row exists for this;
- [§3]'s amended limitation paragraph.

If any of the three is lost in a later edit, a stated assumption becomes an unstated one. Over all arms
the same count is **82**, because no treated patient carries a reason either — so any count of this
group that is not arm-restricted is measuring nothing.

### 21.5 Verification record for the amendment

Run 2026-08-13, on v7, by driving both rules over the same frame and comparing.

| Claim | Verified |
|---|---|
| DECISION 1a classifies **107 eligible / 19 ineligible**; DECISION 1 classified 64 / 43 / 19; the 43 move to `eligible` and nothing else moves | yes, run |
| The retained set and the primary cohort are **identical** under both rules, compared identifier for identifier | yes, run |
| No treated patient carries the flag, so restoring the revealed-fact mask changes **no** record — the deletion is a no-op on this workbook | yes, run |
| §7.1's table under DECISION 1a, cell for cell, with every `n` unchanged | yes, run |
| **43** control-arm records carry no documented reason, at CHUV and Lugano only; over all arms it is **82** at 3 centres | yes, run |
| A reason is recorded on **44** records — exactly the HUG and USZ controls — and is now read by no classifier | yes, run |
| Restricting the comparator arm to documented controls leaves **11**, all at HUG, and a single-centre cohort of 41 | yes, run |
| The audit log reproduces byte-identically under `PYTHONHASHSEED=0` and `=1` after the amendment | yes, run |
| The suite is green with `data/` present (490 tests) and with `data/` absent (455 passed, 35 skipped, no collection error) | yes, run |
