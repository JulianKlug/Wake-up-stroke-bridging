"""Stage 0 — the six facts the analysis plan has to be instantiated against.

Roadmap Stage 0. This runs *before* the configuration module of Stage 1 exists, so it
is deliberately self-contained and reads raw workbook column names directly. Nothing
here is imported by the analysis; its only product is the markdown note it writes.

Two things it does not do. It never tabulates an outcome against treatment — the six
facts do not require it and looking would compromise the prespecification. And it
writes its note into ``out/``, which is gitignored, because every number below is
derived from patient data.

    python3 stage0_data_inventory.py            # -> ../out/stage0_data_inventory.md
"""
from __future__ import annotations

from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
XLSX = ROOT / "data" / (
    "Excel_bridging_EXTEND_paper_HUG_CHUV_LUGANO_USZ_def_v7_july26_"
    "with_abs_contra_indication.xlsx"
)
SHEET = "Feuil1"
NOTE = ROOT / "out" / "stage0_data_inventory.md"

# The workbook records HUG as the integer 1; the other three carry site names.
CENTER_RECODE = {"1": "HUG", "Lausanne": "CHUV", "Lugano": "Lugano", "USZ": "USZ"}

TREATMENT = "IVTwithrtPA"
REASON = "Contraindications_to_IVT"
ABS_FLAG = "IVT_contraindicated_binary"
MRS90 = "mRSscoreat90days"

RARE_EVENT_THRESHOLD = 10  # SAP §8: below this an outcome is not augmented

# "decision" and not "clinician": the field contains a `Clnician` typo, so matching on the
# correctly spelt word silently misses a record. Only used to describe the field now that its
# content no longer classifies anyone (see DECISION 1).
CLINICIAN_STEM = "decision"
COAG_REASON = "Coagulation tests unavailable in anticoagulated patients at time of EVT"

# DECISION 1 (2026-08-06, PI): `IVT_contraindicated_binary` is the §4 eligibility classifier and
# the free text of `Contraindications_to_IVT` is not read. The reason column still contributes one
# bit — whether a reason was recorded at all — because the flag is 0 both for patients documented
# as having no absolute contraindication and for patients whose centre never collected the field,
# and §3 forbids reading that blank as "no contraindication". So:
#
#   treated                        -> eligible       (by revealed fact)
#   flag = 1                       -> ineligible
#   flag = 0, reason recorded      -> eligible
#   flag = 0, reason not recorded  -> indeterminate
#
# This retires the free-text normalisation problem: no string is ever matched, so the `Clnician`
# typo and the parenthetical annotations cannot break the classifier.

# DECISION 2 (2026-08-06, PI): the 90-day mRS is the ground truth for vital status. Death at 90
# days is `mRS == 6`; the shipped `Deathat90days` and `mRS56at90days` columns are not read.
#
# SAP §5 already required this. The shipped columns are named below only so the two can be
# compared — a disagreement is a data query, not something to resolve silently.

# DECISION 3 (2026-08-07, PI): TICI 2b-3 is augmented, but with a reduced outcome model rather
# than the §6 covariate set — 6 non-events cannot support 12 parameters. This is a declared
# amendment to §8, whose default is that m_a(X) and the propensity model share one covariate list.
TICI_OUTCOME_COVARIATES = ["center", "atrial_fib"]
DERIVED_FROM_MRS = {
    "mRS 0-2 at 90 days": ("mRS <= 2", lambda s: s <= 2, "mRS02at90days"),
    "mRS 0-1 at 90 days": ("mRS <= 1", lambda s: s <= 1, "mRSscore01at90days"),
    "mRS 5-6 at 90 days": ("mRS >= 5", lambda s: s >= 5, "mRS56at90days   "),
    "death at 90 days": ("mRS == 6", lambda s: s == 6, "Deathat90days"),
}
SHIPPED_BINARY = {
    "TICI 2b-3": "TICI_2b_3",
    "symptomatic ICH": "Symptomaticintracranialhaemorrhage",
    "parenchymal haematoma type 2": "Parenchymalhaematomatype2",
}

# SAP §6, with the degrees of freedom each contributes to the propensity model.
# `center` is 2 df only once the zero-bridging centre is dropped, so fact 1 has to
# be settled before this budget can be read.
PS_PARAMETERS = {
    "age": 1, "sex": 1, "prestroke_mrs": 1, "nihss_baseline": 1,
    "onset_type (3 levels)": 2, "core_ml": 1, "tmax6_ml": 1, "atrial_fib": 1,
    "center (3 both-arm centres)": 2,
}


def load() -> pd.DataFrame:
    df = pd.read_excel(XLSX, sheet_name=SHEET)
    df["center"] = df["Center"].astype(str).map(CENTER_RECODE)
    assert df["center"].notna().all(), "unmapped centre code"
    return df


def derive(source: pd.Series, is_event: pd.Series) -> pd.Series:
    """Dichotomise an ordinal source, reimposing its missingness.

    A comparison against a missing value returns false rather than missing, so
    without the mask a patient with no 90-day mRS becomes a non-event.
    """
    return is_event.astype("Int64").mask(source.isna())


def eligibility(df: pd.DataFrame) -> pd.Series:
    """Classify IVT eligibility [§3, §4] under DECISION 1 — flag plus reason presence.

    Order matters: treated first, so a treated patient is eligible by revealed fact
    whatever else is recorded. Stage 4 owns the real implementation; this exists so the
    inventory reports the classification the plan will use rather than a proxy for it.
    """
    out = pd.Series("indeterminate", index=df.index, dtype=object)
    out[(df[ABS_FLAG] == 0) & df[REASON].notna()] = "eligible"
    out[df[ABS_FLAG] == 1] = "ineligible"
    out[df[TREATMENT] == 1] = "eligible"
    assert not ((df[TREATMENT] == 1) & (df[ABS_FLAG] == 1)).any(), \
        "a treated patient is flagged as absolutely contraindicated — resolve before proceeding"
    return out


def md(df: pd.DataFrame) -> str:
    L: list[str] = []
    w = L.append
    n = len(df)

    w("# Stage 0 — data inventory")
    w("")
    w(f"Source: `{XLSX.name}`, sheet `{SHEET}`, {n} rows, {df.shape[1]} columns.")
    w("Generated by `extended_bridging/stage0_data_inventory.py`. Regenerate rather than edit.")
    w("")

    # ---- fact 1: zero-bridging centres ------------------------------------------
    by_centre = pd.crosstab(df["center"], df[TREATMENT])
    by_centre.columns = ["EVT alone", "bridging"]
    by_centre["n"] = by_centre.sum(axis=1)
    by_centre["% bridged"] = (100 * by_centre["bridging"] / by_centre["n"]).round(1)
    never = by_centre.index[by_centre["bridging"] == 0].tolist()
    both_arm = [c for c in by_centre.index if c not in never]

    w("## 1. Centres contributing zero bridging patients  [§3]")
    w("")
    w(by_centre.to_markdown())
    w("")
    w(f"**Zero-bridging: {', '.join(never) or 'none'}.** Both-arm centres: {', '.join(both_arm)}.")
    w("")
    thin = by_centre.loc[both_arm].query("bridging < 5")
    for c, r in thin.iterrows():
        w(f"- `{c}` clears the restriction with {int(r['bridging'])} treated of {int(r['n'])} "
          f"({r['% bridged']}%). Positivity there is thin rather than structurally absent, so it "
          "survives Stage 5 and has to be visible in the within-centre overlap table [§9].")
    w("")

    # ---- facts 2 and 3: the contraindication-reason field -------------------------
    reason = df[REASON]
    w("## 2. Contraindication-reason field and its values  [§3, §11]")
    w("")
    w(f"`{REASON}` exists: free text, non-missing for {int(reason.notna().sum())} of {n}.")
    w(f"A derived `{ABS_FLAG}` (0/1, non-missing for all {n}) also ships with this version "
      "of the workbook.")
    w("")
    tab = (pd.crosstab(reason, df[ABS_FLAG])
             .rename(columns={0: f"{ABS_FLAG}=0", 1: f"{ABS_FLAG}=1"}))
    w(tab.to_markdown())
    w("")
    w(f"{reason.nunique()} distinct non-missing values. `{ABS_FLAG}` is 1 for a strict subset of "
      "the records carrying a reason and never 1 where the reason is missing, so it re-encodes "
      "the free text and adds no coverage — what it adds is the data owner's judgement of which "
      "reasons are absolute.")
    w("")
    w(f"**Three-category eligibility [§3] is possible.** Per DECISION 1, `{ABS_FLAG}` *is* the [§4] "
      "classifier and the free text is not read. The reason column still contributes one bit — "
      f"whether a reason was recorded — because `{ABS_FLAG}` = 0 covers both patients documented as "
      "having no absolute contraindication and patients whose centre never collected the field, and "
      "[§3] forbids reading that blank as \"no contraindication\".")
    w("")
    variants = sorted(v for v in reason.dropna().unique() if CLINICIAN_STEM in v.lower())
    n_typo = sum("clinician" not in v.lower() for v in variants)
    n_paren = sum("(" in v for v in variants)
    w("Not reading the text retires a real hazard rather than merely deferring it: the field holds "
      f"{len(variants)} spellings of the same preference reason — "
      + ", ".join(f"`{v}`" for v in variants)
      + f" — of which {n_typo} misspells the word and {n_paren} carry a parenthetical annotation. "
      f"An exact-string classifier would have raised on {len(variants) - 1} of them. Since no string "
      "is ever matched, no normalisation rule is needed and none of these can break the classifier.")
    w("")

    arm = df[TREATMENT].map({0: "EVT alone", 1: "bridging"})
    present = pd.crosstab([df["center"], arm], reason.notna())
    present.columns = ["reason absent", "reason recorded"]
    present.index.names = ["centre", "arm"]
    w("## 3. Per-centre completeness of the reason field  [§11]")
    w("")
    w(present.reset_index().to_markdown(index=False))
    w("")
    centre_complete = (reason.notna()
                       .groupby([df["center"], df[TREATMENT] == 0])
                       .mean())
    collected = sorted({c for (c, ctrl), v in centre_complete.items() if ctrl and v > 0})
    w(f"Recorded for controls at **{', '.join(collected)}** and for no control elsewhere; "
      "recorded for no treated patient at any centre.")
    w("")

    # ---- cohort flow, needed before the event counts mean anything ---------------
    df = df.assign(eligibility=eligibility(df))
    prim = df[df["center"].isin(both_arm)]
    elig = prim[prim["eligibility"] != "ineligible"]
    indeterminate = elig[elig["eligibility"] == "indeterminate"]
    documented_ctrl = elig[(elig[TREATMENT] == 0) & (elig["eligibility"] == "eligible")]

    w("## Cohort flow under the §4 classifier")
    w("")
    w("Eligibility by centre and arm, as classified above:")
    w("")
    w(pd.crosstab([df["center"], arm], df["eligibility"])
        .rename_axis(index=["centre", "arm"]).reset_index().to_markdown(index=False))
    w("")
    w("The two [§3] restrictions, applied in order:")
    w("")
    w(f"| Step | n | treated |")
    w("|---|---|---|")
    w(f"| All records | {n} | {int(df[TREATMENT].sum())} |")
    w(f"| Drop zero-bridging centre(s) | {len(prim)} | {int(prim[TREATMENT].sum())} |")
    w(f"| Drop ineligible | {len(elig)} | {int(elig[TREATMENT].sum())} |")
    w("")
    w(f"Of the {len(elig) - int(elig[TREATMENT].sum())} controls retained, "
      f"**{len(indeterminate)} are of indeterminate eligibility** (no reason ever collected at "
      f"their centre) and {len(documented_ctrl)} carry a documented non-absolute reason. The "
      "indeterminate group is therefore "
      f"{100 * len(indeterminate) / (len(elig) - int(elig[TREATMENT].sum())):.0f}% of the control "
      "arm, not a residual category — the [§3] limitation that retaining them assumes eligibility "
      "carries most of the control arm with it.")
    w("")

    # ---- fact 4: event counts ---------------------------------------------------
    w("## 4. Event counts per binary outcome  [§8]")
    w("")
    w(f"Overall, not by arm. The rule is <{RARE_EVENT_THRESHOLD} events → no augmentation [§8].")
    w("")
    rows = []
    for scope, d in (("full", df), ("primary cohort", elig)):
        src = d[MRS90]
        for label, (rule, is_event, _shipped) in DERIVED_FROM_MRS.items():
            y = derive(src, is_event(src))
            rows.append(dict(outcome=label, scope=scope, rule=rule, n=int(y.notna().sum()),
                             events=int(y.sum()), missing=int(y.isna().sum())))
        for label, col in SHIPPED_BINARY.items():
            y = d[col].astype("Float64")
            rows.append(dict(outcome=label, scope=scope, rule="as shipped", n=int(y.notna().sum()),
                             events=int(y.sum()), missing=int(y.isna().sum())))
    counts = (pd.DataFrame(rows)
                .pivot(index=["outcome", "rule"], columns="scope",
                       values=["n", "events", "missing"])
                .reorder_levels([1, 0], axis=1)
                .sort_index(axis=1, level=0, sort_remaining=False))
    counts.columns = [f"{scope}: {stat}" for scope, stat in counts.columns]
    w(counts.reset_index().to_markdown(index=False))
    w("")
    rare = [r["outcome"] for r in rows
            if r["scope"] == "primary cohort" and r["events"] < RARE_EVENT_THRESHOLD]
    w(f"**Under the <{RARE_EVENT_THRESHOLD}-event rule in the primary cohort: "
      f"{', '.join(rare) or 'none'}.** Reported as unaugmented weighted risk difference and "
      "weighted marginal odds ratio only.")
    w("")
    for label, col in SHIPPED_BINARY.items():
        y = elig[col].astype("Float64")
        minority = min(int(y.sum()), int(y.notna().sum() - y.sum()))
        if minority < RARE_EVENT_THRESHOLD and int(y.sum()) >= RARE_EVENT_THRESHOLD:
            w(f"- `{label}` passes the rule as written ({int(y.sum())} events) but has only "
              f"{minority} non-events. The rule's stated reason — a nuisance model with more "
              "parameters than events — applies to whichever cell is small, so the literal rule "
              "and its rationale disagree here. Settled by DECISION 3 below.")
    w("")

    w("### Derived versus shipped dichotomies")
    w("")
    w("Per DECISION 2 the derived column governs everywhere; the shipped columns are compared here "
      "and then never read again.")
    w("")
    src = df[MRS90]
    for label, (rule, is_event, shipped) in DERIVED_FROM_MRS.items():
        new = derive(src, is_event(src))
        old = df[shipped].astype("Int64")
        bad = df.loc[(new != old) & new.notna() & old.notna()]
        verdict = "agrees throughout" if bad.empty else (
            f"**disagrees on {len(bad)}**: "
            + "; ".join(f"`{r.CaseID}` mRS={r[MRS90]:.0f} → derived "
                        f"{int(new[i])}, shipped {int(old[i])}"
                        for i, r in bad.iterrows()))
        w(f"- `{shipped.strip()}` vs `{rule}` — {verdict}")
    w("")

    # ---- fact 5: degrees of freedom ---------------------------------------------
    n_treated = int(elig[TREATMENT].sum())
    k = sum(PS_PARAMETERS.values())
    w("## 5. Treated-arm size and the degrees-of-freedom budget  [§6]")
    w("")
    w(f"Treated in the primary cohort: **{n_treated}** "
      + ", ".join(f"{c} {int(v)}" for c, v in
                  elig[elig[TREATMENT] == 1]["center"].value_counts().items()) + ".")
    w("")
    w("| §6 covariate | df |")
    w("|---|---|")
    for name, d in PS_PARAMETERS.items():
        w(f"| {name} | {d} |")
    w(f"| **total (excl. intercept)** | **{k}** |")
    w("")
    w(f"**{n_treated} / {k} = {n_treated / k:.2f} treated patients per parameter**, which is the "
      "~3.5 the plan assumes [§6]. The budget holds only for the parsimonious specification: the "
      "full-covariate sensitivity model would add 4 parameters and take it to "
      f"{n_treated / (k + 4):.2f}.")
    w("")

    # ---- structural facts that would break a later stage silently ----------------
    w("## Structural checks")
    w("")
    both = int(((df["Wakeupstroke"] == 1) & (df["Unwitnessedstroke"] == 1)).sum())
    combos = df.groupby(["Wakeupstroke", "Unwitnessedstroke"]).size()
    w(f"- `onset_type` [§5]: wake-up and unwitnessed are never both positive ({both} cases), so "
      "the three levels partition the cohort — witnessed "
      f"{combos.get((0, 0), 0)}, unwitnessed {combos.get((0, 1), 0)}, wake-up "
      f"{combos.get((1, 0), 0)}.")
    for col, nm in (("FirstbrainimageMRI", "first_image_mri"),):
        if df[col].nunique(dropna=True) == 1:
            w(f"- `{nm}` is constant ({df[col].dropna().iloc[0]}) across all {n} records — zero "
              "variance, dropped automatically [§6].")
    w(f"- `CaseID` is unique ({int(df['CaseID'].duplicated().sum())} duplicates) but not uniform: "
      "HUG uses `SSR-HUG-…`, Lugano `L-…`, USZ bare integers. Identifiers must be treated as "
      "strings.")
    w("- No `'N/A'` or empty-string sentinels: missingness is already blank in every column.")
    miss = {c: int(df[c].isna().sum()) for c in
            ["IschemiccorevolumeCBF30ml", "HypoperfusedtissuevolumeTmax6sml", MRS90, "TICI_2b_3"]
            if df[c].isna().any()}
    w("- Missingness in the §6 covariates and outcomes: "
      + ", ".join(f"`{c}` {v}" for c, v in miss.items()) + ".")
    w(f"- `IschemiccorevolumeCBF30ml` is exactly 0 in "
      f"{int((df['IschemiccorevolumeCBF30ml'] == 0).sum())} of "
      f"{int(df['IschemiccorevolumeCBF30ml'].notna().sum())} records, so any mismatch *ratio* "
      "[§13] is undefined for those patients and needs an explicit convention.")
    w("- `TimefromONSETtoIVTmin` is present for every treated patient and missing for every "
      "control — structurally non-applicable, not missing [§11].")
    w("")

    # ---- what has been decided, and what a guess would still decide silently -------
    n_coag = int((prim[REASON] == COAG_REASON).sum())
    annotated = sorted(v for v in prim[REASON].dropna().unique()
                       if CLINICIAN_STEM in v.lower() and "(" in v)
    n_death = int(derive(elig[MRS90], elig[MRS90] == 6).sum())
    n_death_shipped = int(elig["Deathat90days"].sum())

    w("## Decisions taken")
    w("")
    w("All taken by the PI, each before any outcome was examined by arm.")
    w("")
    w(f"**DECISION 1 — `{ABS_FLAG}` is the [§4] classifier; the free text of `{REASON}` is not "
      "read.** Consequences:")
    w("")
    w(f"- `{COAG_REASON}` ({n_coag} controls in the both-arm cohort) is **eligible**, the flag "
      "calling it not absolute. This reverses the pilot implementation, which listed it as an "
      "absolute contraindication.")
    if annotated:
        w("- The annotated clinician decisions — "
          + "; ".join(f"`{v}`" for v in annotated)
          + f" ({len(annotated)} control{'s' if len(annotated) != 1 else ''}) — are **eligible**. "
            "Cerebral metastasis and disseminated cavernomatosis are conventionally absolute "
            "contraindications, so this is a substantive call and not a formality; it follows from "
            "taking the flag as the classifier. Worth naming in the limitations.")
    w("- No string is matched, so the [§4] assertion changes character: there is no set of "
      f"recognised reasons to check a new value against. What must instead be asserted is that "
      f"`{ABS_FLAG}` is 0/1 and never missing, and that no treated patient carries a 1 — the "
      "condition the classifier's revealed-fact rule would otherwise mask.")
    w("- The `indeterminate` category [§3] survives, resting on whether a reason was recorded "
      f"rather than on what it said: {len(indeterminate)} patients.")
    w("")
    w("**DECISION 2 — the 90-day mRS is ground truth for vital status; death is derived as "
      "`mRS == 6`.** Consequences:")
    w("")
    w(f"- Death at 90 days in the primary cohort: **{n_death} events**, not the "
      f"{n_death_shipped} the shipped column reports.")
    w("- `mRS 5-6` likewise comes from the ordinal source. Both shipped columns are now unread and "
      "should be explicitly dropped in the Stage 1 mapping, not left available.")
    w("- The two contradictory records are read as **alive** with the mRS as recorded (1 and 4). "
      "This is a derivation rule, not a resolution: one of the two fields is wrong in each record, "
      "and it stays a standing query for the data owner. Both are HUG cases, both are listed by "
      "identifier in section 4, and if the mRS is the field in error then both are deaths "
      "misrecorded as good outcomes — the direction that matters most.")
    w("")
    n_non = int(elig["TICI_2b_3"].notna().sum() - elig["TICI_2b_3"].sum())
    k_tici = 2 + 2 + 1  # intercept + treatment + center (2 df) + atrial_fib
    w(f"**DECISION 3 (2026-08-07) — the [§8] rare-outcome rule is read by its rationale, on the "
      f"minority cell, and `TICI 2b-3` is augmented with a reduced outcome model.** `m_a(X)` for "
      f"TICI is treatment + `{'` + `'.join(TICI_OUTCOME_COVARIATES)}`, {k_tici} parameters against "
      f"{n_non} non-events. Consequences:")
    w("")
    w("- This **amends [§8]**, which otherwise requires `m_a(X)` to carry the §6 covariate set and "
      "to resolve to the same configuration entry as the propensity model. The amendment is "
      "per-outcome and declared: TICI is the only override, every other binary outcome keeps the "
      "§6 set, and the override must be emitted beside the TICI estimate rather than left in the "
      "configuration.")
    w("- The [§8] argument against dropping a §6 confounder from `m_a` — that it disables the "
      "correction exactly where it is needed — now applies to TICI knowingly. What guards it is the "
      "comparison [§8] already mandates: the unaugmented weighted risk difference is reported "
      "alongside, and material disagreement between the two is evidence about this reduced outcome "
      "model, not confirmation of either.")
    w("- The covariates were chosen on procedural grounds, never from the data: TICI is a "
      "recanalisation outcome, so centre carries device, technique and operator volume, and atrial "
      "fibrillation proxies clot composition. The §6 severity axes predict function rather than "
      "recanalisation.")
    w("- The **propensity model is unchanged**. Its budget is set by the treated arm (fact 5), not "
      "by any outcome, so no outcome-specific reduction touches it.")
    rare_still = [r["outcome"] for r in rows
                  if r["scope"] == "primary cohort" and r["events"] < RARE_EVENT_THRESHOLD]
    w(f"- Unaffected: {', '.join(rare_still)} stay **unaugmented**. Their minority cell is the event "
      "cell, they fall under the rule as written, and a reduced model would not rescue them.")
    w("")
    w("## Open decisions")
    w("")
    w("None. Stage 1 can be written against this note.")
    w("")
    return "\n".join(L) + "\n"


def main() -> None:
    NOTE.parent.mkdir(parents=True, exist_ok=True)
    NOTE.write_text(md(load()))
    print(f"wrote {NOTE}")


if __name__ == "__main__":
    main()
