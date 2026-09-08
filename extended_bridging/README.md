# IVT before thrombectomy versus thrombectomy alone in late-window, CTP-selected stroke

Analysis code for the extended-window bridging study: a four-centre cohort (HUG, CHUV, Lugano, USZ),
90-day mRS as the primary outcome, overlap weighting for the primary estimand and pooled ordinal
standardisation for the all-centre and feasible-policy analyses.

- `statistical_analysis_plan.md` — the SAP. Every `[§n]` citation in the code points here.
- `implementation_roadmap.md` — the fourteen stages, what each landed, and the analysis-lock milestone.
- `specs/` — one specification per stage, written before its implementation.
- `../TODOS.md` — deferred work, each item naming what blocks it and what it unblocks.

## Patient data availability

**The patient data cannot be shared. No version of it is in this repository, and none will
accompany publication or any archive release.** `data/` is gitignored, as is `out/`, which carries the
estimate. The workbook is verified by SHA-256 at load
(`config.DATA_SHA256`), so the analysis either runs against the exact file the results were produced
from or refuses to run.

### What this means for a checkout without the data

| | Runs | Notes |
|---|---|---|
| The test suite | Mostly | 130 of roughly 1300 test functions carry the data gate and skip; the rest run on frames, synthetic arrays and `tests/fixture_schema.xlsx`. |
| The pipeline end to end | **No** | `python -m report` needs the private workbook. The one committed fixture, `tests/fixture_schema.xlsx`, stops at Stage 5's both-arms postcondition — it reaches no cohort with a control arm at a treating centre, and neither does the hand-built frame the earlier stages test on. |
| Estimator internals | Yes | Every fitter is exercised on synthetic data with known answers, and four R oracles (`logistf`, `PSweight`, `ordinal::clm`, `ordinal::clmm`) plus statsmodels' `OrderedModel` check them independently where those packages are available. |

A committed both-arm fixture that would make the whole pipeline runnable from a checkout is filed in
`../TODOS.md` and has not been built. Because the data will not be released either, such a fixture is
the **only** way any reader outside the study team could ever execute this pipeline; the tests above
are otherwise the whole of what is externally runnable.

## Running it

```sh
uv run pytest                       # the suite; slow gates are opt-in per their skip messages
uv run python -m report             # Stages 1-13, then 22 tables, 5 figures and the audit log to out/
```

Reproducibility rests on a recorded seed (`config.SEED`) rather than on the environment: the
dependency bounds in `pyproject.toml` are pinned because numpy's NEP 50 promotion and pandas 3.0's
string-inference default both move results, and the notes there say where.

`notebooks/` explains how the estimates are derived, one notebook per family — the primary common odds
ratio, the binary secondary and safety outcomes, and the §14 standardisation and policy analyses. Each
calls the shipped modules rather than reimplementing them, so it needs the data to execute; the
committed copies are output-stripped by a clean filter (see `.gitattributes`).
