# IVT before thrombectomy versus thrombectomy alone in late-window

Analysis code for the extended-window bridging study: a four-centre cohort (HUG, CHUV, Lugano, USZ),
90-day mRS as the primary outcome, overlap weighting for the primary estimand and pooled ordinal
standardisation for the all-centre and feasible-policy analyses.

- `statistical_analysis_plan.md` — the SAP. Every `[§n]` citation in the code points here.
- `implementation_roadmap.md` — the fourteen stages, what each landed, and the analysis-lock milestone.
- `specs/` — one specification per stage, written before its implementation.


```sh
uv run pytest                       # the suite; slow gates are opt-in per their skip messages
uv run python -m report             # Stages 1-13, then 22 tables, 5 figures and the audit log to out/
```

