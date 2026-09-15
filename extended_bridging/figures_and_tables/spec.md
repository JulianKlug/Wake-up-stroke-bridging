# Manuscript exhibits — specification

The four figures and two tables the paper carries, one notebook each, and the supplement of
model results they leave out. Section references in
brackets are to `../statistical_analysis_plan.md`; `SMALL_CAPS` constants live in `../config.py`; the
reporting layer this one sits on is specified by `../specs/stage14_outputs_and_guardrails.md`.

**Precondition.** Written against the landed Stages 1–14 and the analysis lock. This layer computes
no estimate, interval or p-value. It fits nothing, records nothing to the audit, and writes into a
directory `report.write` never touches.

**Status.** Written 2026-09-11, with the implementation. It takes no DECISION and changes no estimand.

---

## 1. Why the manuscript set is not the analysis set

`report.py` emits 22 tables and 5 figures in `OUTPUT_IDS` order, which is the order the *analysis*
produced them in. A manuscript argues, and the argument has an order:

```
Figure 1   who was compared          →  cohort flow
Figure 2   were they comparable      →  balance (Love plot) + propensity overlap by centre
Figure 3   what happened             →  weighted mRS distribution + the six cumulative RD_k
Figure 4   what else happened        →  secondary and safety risk differences
```

Two results make this ordering a requirement rather than a preference.

**T07 cannot be the main table unchanged.** It carries the right quantities — a common odds ratio and
six cumulative risk differences — but the differences do not share a sign, and [§16a] requires a reader
told that a single odds ratio summarising them is not a null result. A one-row table states the odds
ratio and buries the disagreement three rows down; a forest plot puts them in the same glance.

**The design diagnostics are not supplementary.** Treatment is nearly determined by centre. After
weighting, five covariates remain above [§9]'s threshold, and centre support is grossly uneven — one
centre contributes almost all the bridging arm, another almost all the control arm, and a third is
excluded for structural non-positivity. [§9]'s amendment already requires the exceedances named beside
the estimate; a figure showing *why* they are there belongs beside it too.

Everything else — T03, T04, T05, T06, T09, T13/T14, T16/T17, T19/T20, F03, F05 and the [§14]
transported analyses — stays exactly as `report.py` writes it and becomes supplementary material. The
[§14] analyses in particular are kept out of the main figures because they target different
populations and [§16] forbids reading them as like-for-like.

## 2. What this layer may do

`manuscript.py`'s neighbour below is `report.py`. It may call `report.*` and `config.*`, and it may
read any attribute of a result object it holds. It may **not** call an estimator module — including
`balance.levels`, which is why Table 1 routes through `report.baseline_rows`. `tests/test_manuscript.py`
asserts the import set.

`report.py` learns nothing about this layer: no upward import, and no manuscript id in `OUTPUT_IDS`,
whose comment now says what it is the one list *of*.

## 3. The cache

Seven notebooks, one run. `manuscript.load()` returns a `Bundle` — the Stage 1–11 result objects, which
is exactly what the seven exhibits read — pickled under a key that names everything capable of moving a
number: `DATA_SHA256`, a digest of the shipped modules, `SEED`, `N_BOOT`, `CI_LEVEL`,
`PERCENTILE_METHOD`, `BOOT_STRATUM`, the interpreter and the numpy and pandas versions
(`pyproject.toml` says why the last two are pinned), and a schema number.

Three properties of the record are load-bearing.

- **It is a `Bundle`, not a `report.Run`.** A `Run` also carries the [§14] objects, and its V4 guard
  exists so a partial one cannot reach `report.write` and print T22 without its comparators. Filling
  those fields with placeholders to satisfy the guard would build an object that lies; leaving them
  `None` is what the guard forbids. The field names match `Run`'s, which is what lets `report`'s row
  builders read a `Bundle` unchanged.
- **`__post_init__` is re-run on load.** `pickle` restores a frozen dataclass through `__dict__` and
  never calls it, so a cache written from a partial build would otherwise load silently.
- **Nothing unpickled re-enters an estimator.** `C.Specification` is compared by identity in four
  shipped call sites and `pickle` restores a copy, so a `ps` off the cache fed back to
  `sensitivity.full_covariate` would raise on a specification that is in fact correct. The bundle is
  input to formatting and drawing, and to nothing else.

## 4. The exhibits

| Exhibit | Panels / blocks | Reads | [§16] statements it must carry |
|---|---|---|---|
| **Figure 1** `cohort_flow` | Five boxes: source records → both [§3] restrictions → cohort → covariate-complete → outcome-present, with per-arm counts; exclusions annotated as computed differences. Also emitted as an **editable** slide diagram — see §7 | `audit.entry("cohort", "cohort_flow")`, `ps.in_model`, `primary.in_estimate` | own denominator per estimate [§11]; entry conditioning [§2] |
| **Figure 2** `design_diagnostics` | **A** Love plot, \|SMD\| before and after weighting, ordered by residual, threshold drawn, exceedances in red. **B** one propensity histogram per centre contributing a weighted contrast, arm counts in the title | `balance.covariates`, `balance.centres`, `balance.pooled`, `ps.e`, `ps.in_model`, `cohort` | `exceedance_clause`; the omitted centres named with their status [§9]; `ATO_DESCRIPTION` |
| **Figure 3** `primary_mrs` | **A** weighted mRS 0–6 by arm. **B** forest of RD₀…RD₅ with percentile intervals, null at 0, common OR annotated | `primary.cumulative`, `primary.rd`, `primary.odds_ratio`, `boot.intervals`, `balance` | `primary_clauses`: `CONSTANT_SHIFT`, `direction_clause`, `exceedance_clause`, `interval_note` |
| **Figure 4** `secondary_and_safety` | One forest; the two families separated by a rule; augmented points indented and marked | `secondary.by_family()`, `boot.intervals` | both families' `binary_rows` clauses — model-assisted, descriptive-only, which BH ran over how many |
| **Table 1** `baseline_and_weighting` | `report.baseline_summary` — `n (%)` for factor levels and `C.BINARY_COLUMNS` flags, `median (IQR)` otherwise — merged with the SMD columns off `balance.covariates`; ESS as the last row | `cohort`, `ps`, `balance.covariates` | `ATO_DESCRIPTION`; the full-confounder-set rule [§9]; `exceedance_clause` |
| **Table S1** `primary_and_sensitivity` | Per specification, a design block (n in model, ESS per arm, worst residual SMD) then an effect block (common OR, six RD_k). The [§13] arm carries no p-value | `primary`, `boot`, `ps`, `balance`, `arm` | `primary_clauses` + `arm_clauses` |
| **Table S2** `secondary_and_safety` | T10 and T11 in one table, keeping the columns through the augmented interval | `report.binary_rows` for both families | both families' `binary_rows` clauses |
| **Table 2** `outcomes_by_arm` | `report.outcome_summary` — the seven mRS classes then each binary outcome, `n (%)` as observed and the share in the weighted population. No effect, interval or p value | `cohort`, `primary.cumulative`, `secondary.estimates` | `DESCRIPTIVE_CRUDE` [§16b]; `ATO_DESCRIPTION`; own denominator per outcome [§11] |

### The main text carries no estimate

Tables 1 and 2 say who was compared and what happened to them. Every effect estimate — the primary
common odds ratio, its [§13] sensitivity arm, the secondary and safety outcomes, and the [§14]
analyses in other populations — is in the supplement, which `supplement_tables.ipynb` compiles from
the two `Table_S*` CSVs plus the locked `T22` and `T16`, and converts to Word with pandoc. That
notebook and its output are **untracked**, like the prose: the paper's documents are the authors'.

## 5. What lands on disk

`OUT/manuscript/`, which `.gitignore` already covers:

```
Figure_n_<slug>.png                    300 dpi
Figure_n_<slug>_source_data.csv        the plotted numbers, then the statements as `# ` lines
Table_n_<slug>.csv                     the table, then the statements as `# ` lines
Figure_1_cohort_flow.md                the slide-deck source for §7's editable diagram
Figure_1_cohort_flow.pptx              rendered from it, by the command in §7
figure_legends.md                      every figure legend, collected from the CSVs above
methods_and_results.md                 the two manuscript sections, numbers interpolated (untracked)
supplement_tables.md / .docx           the model-result supplement, compiled (untracked)
```

**No explanatory text is drawn into a PNG.** Panels are headed by a letter and nothing else; what each
panel shows, and any estimate quoted alongside it, is the legend's job. A PNG cannot carry the [§16]
statements anyway — they are required *beside the estimate* — so the source-data CSV is where the
caption and the statements live. It is rendered by `report.render_csv`, so both output families share
one convention: rows first, `# ` clauses beneath, caption first among them.

`manuscript.save_legends()` then collects those blocks into `figure_legends.md`, one section per
figure — the page a journal asks for separately from the figures. It is rebuilt from whatever is on
disk, so running one notebook leaves it complete rather than partial, and it cannot drift from the
figures because it is read back out of them. The provenance sentence is kept in the file as an HTML
comment: it says which run drew the figure, which the printed legend does not carry.

Figure appearance has none of the reproducibility guarantees the numbers have — `pyproject.toml` pins
numpy and pandas because they move results, and leaves matplotlib unbounded — so every exhibit's
clause block records the matplotlib version alongside the run's `DATA_SHA256`, seed and replicate
count.

## 6. What was promoted in `report.py`, and why

Eight names lost their underscore or were extracted, and `report.py`'s pinned public surface grew from
three names to eleven. The rule governing the list: **every promoted name is still called by
`report.py` itself**, so none can rot into a manuscript-only branch.

| name | still called by | needed by |
|---|---|---|
| `direction_clause` | `_t07` | Figure 3, Table 2 |
| `exceedance_clause` | `_t05`, `_t07`, `_t09` | Figure 2, Tables 1–2 |
| `interval_note` | `_t07`, `binary_rows`, `_t15`, `_t18` | Figures 3–4, Tables 2–3 |
| `pmf` | `_f01`, `_f05`, `_distribution_rows` | Figure 3 |
| `baseline_rows` | `_t02` | — (T02 keeps the means the SMDs summarise) |
| `baseline_summary` | — | Table 1's clinical `n (%)` / `median (IQR)` presentation |
| `outcome_summary` | — | Table 4's observed-and-weighted outcome rows |
| `binary_rows` | `_t10`, `_t11` | Figure 4, Table 3 |
| `arm_clauses` + `AGREEMENT_NOT_REASSURANCE` | `_t12` | Table 2 |
| `render_csv` | `write` | every exhibit |
| `fmt` (alias of `data._fmt`) | — | every exhibit; one formatter across both families |

`baseline_summary` and `outcome_summary` are the two exceptions to "still called by `report.py`".
For the first, deliberately: T02 keeps
the arm **means**, because those are what [§9]'s standardised mean differences summarise and a T02 of
medians would leave T05's numbers with no visible operand. Both walk `C.BALANCE_SET` through
`balance.levels` and `tests/test_manuscript.py` pins that they produce the same rows in the same
order, so the two presentations cannot come to disagree about what a row is. Which rows are counted is
declared — a [§6] factor, or a name in `C.BINARY_COLUMNS` — never inferred from the sample, and the
weighted medians are weighted quantiles under `C.PERCENTILE_METHOD`, the rule the [§10] intervals use.

`_populations_clause` stays private: only the [§14] exhibits need it, and those are supplementary.
The row builders, the figure builders and `_render_md` stay private, so the manuscript cannot reach
past the statements into the analysis outputs' own rendering.

## 7. Figure 1 twice: a picture, and shapes

A matplotlib figure is a picture of a diagram — a co-author cannot move a box in it, and a journal
cannot ask for the flow in its own house style without it being redrawn. The cohort flow is the one
exhibit here that is pure geometry: no data-derived shape, only counts and arrows. So the notebook
also writes a slide-deck source whose boxes, connectors and exclusion annotations render as **native**
PowerPoint shapes, editable in PowerPoint, Keynote and LibreOffice Impress.

The counts in it come from the same `steps` and `drops` the PNG is drawn from —
`tests/test_manuscript.py` pins that, because an editable box is also a box someone can retype.

Rendering is a separate tool and a separate command, which is why `manuscript.py` writes the spec and
does not shell out:

```sh
PYTHONPATH=~/.claude/skills/powerpoint python3 ~/.claude/skills/powerpoint/__main__.py \
    build out/manuscript/Figure_1_cohort_flow.md --out out/manuscript/Figure_1_cohort_flow.pptx
```

The deck is two slides: the brand's cover, then the diagram on a full-width slot with the layout's
chrome suppressed (`- notitle: true`). Delete the cover if only the figure is wanted. The renderer is
an external, machine-local skill and is deliberately not a dependency of this repository — without it
the `.md` is still a complete, readable specification of the diagram.

## 8. The prose

`manuscript_text.ipynb` writes the Methods and Results. Length and register follow the reference
research letter in `../../data/`: no subsections, past tense, roughly two hundred words and a hundred
and fifty.

**Neither the notebook nor its output is tracked.** `.gitignore` names the notebook and already covers
`out/`. The estimates, the exhibits and the code that makes them are this repository's; the paper's
text is the authors' and does not live under version control here. The notebook sits beside the
exhibit notebooks and runs exactly as they do — a checkout simply does not carry it, so `spec.md` is
the record that it exists and what it does.

**No number in it is typed.** Prose is the one manuscript surface a regeneration does not touch — a
figure is redrawn when the analysis moves, a sentence is left behind — so every count, estimate,
interval and p value is interpolated from the same bundle the exhibits are drawn from, and
`tests/test_manuscript.py` checks the quoted primary result and sensitivity arm against T07 and T12.
Covariate names are mapped to their clinical terms by a dictionary that asserts it covers
`C.PS_COVARIATES`, so a confounder added to [§6] and not named in the Methods fails rather than
silently vanishes.

The [§14a] and [§14b] analyses are described and quoted too, and they are the one place the prose
does not read the `Bundle`: those objects are Stage 12 and 13, the eighteen minutes of fitting no main
exhibit needs, so their numbers come from the locked `T22` and `T16` instead — quoted, still not
typed, and pinned by a test. [§16] requires the different-populations statement wherever a [§14]
estimate appears beside the primary, and a second test asserts the text carries it and never calls a
transported analysis a sensitivity analysis.

Two things the run cannot supply are left as `<!-- TO COMPLETE -->` comments and counted by a test:
the source cohort's inclusion criteria, and the ethics statement. Neither is in this repository and
neither is invented. The [§2] entry-conditioning limitation belongs in the Discussion, which this
notebook does not write.

## 9. Accept when

- Every exhibit's numbers match the locked `out/tables/*.csv` they derive from — T01, T05/T06, T07/T08,
  T10/T11, T12, T02 respectively.
- The seven notebooks execute headless from a cold cache and again from a warm one.
- `tests/test_manuscript.py` is green: the two families stay disjoint, the import set holds, the cache
  key moves when any of its eleven components moves, and every [§16] statement is present because it
  was computed from the run rather than typed into a caption.
- Figure 1's slide deck lints clean and its boxes carry the run's counts, not retyped ones.
- The Methods and Results quote the locked T07 and T12 numbers, and the two placeholders are still
  marked.
