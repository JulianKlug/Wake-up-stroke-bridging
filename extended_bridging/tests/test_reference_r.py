"""Stage 6 §16b.2 — the gated R oracle. The ONLY Python in this repository that knows R exists.

Imported by nothing, on no estimation path, and never producing a number that is reported. §12.7's
scan pins R's blast radius to this one file: `subprocess`, `Rscript`, `logistf` and `PSweight` appear
here and nowhere else under `extended_bridging/`.

**Why an oracle and not a dependency** [§16b.1]. R is where the methodologically authoritative
implementations live — `logistf` is the reference Firth implementation and `PSweight` is the reference
implementation of this SAP's exact estimand, written by the authors of the overlap-weight method — and
that is said plainly rather than worked around. It is declined as a *runtime* dependency because it
doubles the reproducibility surface, because the bridge is the fragile part (a serialisation boundary
crossed on every one of `N_BOOT` replicates), and because adopting `PSweight` would be adopting its
conventions for weights, variance and balance where [§7], [§8] and [§10] are already specified down to
the p-value's definition.

**The boundary is a CSV written to `tmp_path` and a CSV read back**, never `rpy2`: no embedded
interpreter, no shared-library build requirement, no R x Python x numpy compatibility matrix. The cost
is one process spawn per test, in a test that runs once — not 2000 times inside a bootstrap, which is
exactly the distinction §16b.1 turns on. Floats cross at `%.17g` and come back under
`options(digits = 17)`, so a disagreement is the estimator's and never the serialisation's, and
`Rscript --vanilla` means no user or site profile can change the result.

**The skip is loud on purpose.** Every AST scan in Stages 1-5 has a companion proving it fires,
precisely so a check matching nothing cannot pass as green. The same hazard applies here in a worse
form: this oracle skips on any machine without R and the two packages, which is *most* machines, so
the default state is "not run". `pytest -rs` therefore reads as an instruction rather than as noise.

**STATUS: the gate has been opened and passed, 2026-08-16 — DoD-16 is satisfied.** Measured against
`logistf` 1.26.1 and `PSweight` 2.1.2 on R 4.5.0; §18g carries the numbers.

**What it took to open it, because the obstacle was environmental and will recur.** R's startup runs
`system("uname -a", intern = TRUE)`. On this machine PATH carried a second toolchain prefix whose
`sh`/`uname` are linked against a different libc, so that call segfaulted with status 139 inside R's
forked child, `utils` never loaded, and `install.packages` did not exist — which looks exactly like a
broken R and is not one. `_r_environment` below is the fix, and it is why every subprocess in this
module gets an explicit environment rather than inheriting one. Two intermediate diagnoses were wrong
and are recorded so nobody repeats them: it is **not** the sandbox (reproduced with the sandbox
disabled and under `--vanilla --no-init-file`), and it is **not** a broken `r-base-core`.

Two things had to be true beyond that, both about libraries rather than about R: a glibc-matched
`cmake` on PATH for `nloptr` -> `lme4` -> `PSweight`, and a library path excluding
`/usr/local/lib/R/site-library`, which held 75 packages built under R 3.6.3 that abort dependent
builds with "installed before R 4.0.0". Hence `R_LIBS_SITE` and `R_LIBS_USER` being passed through.
"""
from __future__ import annotations

import os
import re
import shutil
import subprocess
from pathlib import Path
from typing import Final

import numpy as np
import pandas as pd
import pytest

import balance
import config
import model
import propensity

RSCRIPT = shutil.which("Rscript")
R_PACKAGES = ("logistf", "PSweight")
REFERENCE = Path(__file__).resolve().parent / "reference"
PROJECT = REFERENCE.parent.parent          # the flat module root, one level above tests/

# Stage 7 §16b. Where an oracle's packages plausibly live when `R_LIBS_USER` is unset — searched, not
# inherited, and that distinction is the whole finding.
_R_USER_LIB_CANDIDATES: Final[tuple[Path, ...]] = (Path.home() / "R" / "library",)


def _r_user_library() -> str | None:
    """A user library that EXISTS, when R_LIBS_USER is unset. §16b.

    Measured 2026-08-18: R's own default (~/R/x86_64-pc-linux-gnu-library/<ver>) does not exist on
    this machine and the oracle's packages are in ~/R/library, so **inheriting R's default finds
    nothing** — and "fall back to R's own default user library" is a no-op, because that default is
    exactly what R already uses. There is no ~/.Renviron either, so `--vanilla` suppresses nothing
    here and stays, because §16b.1 needs a run no user profile can change.

    This is the silently-skipping oracle Stage 6 §16b.2 was written to prevent, arriving by a route
    that section did not anticipate: the gate reported both packages missing while both were
    installed and loadable, on the machine they had been installed on.
    """
    return next((str(p) for p in _R_USER_LIB_CANDIDATES if p.is_dir()), None)


def _r_environment() -> dict[str, str]:
    """The environment every `Rscript` call below is given, and both edits to it are load-bearing.

    **System directories are put FIRST on PATH, not appended.** R's startup runs
    `system("uname -a", intern = TRUE)` through `sh`, so whichever `sh` and `uname` PATH resolves are
    the ones that execute inside R's forked child. On a machine carrying a second toolchain prefix —
    a Gentoo prefix, a conda environment, a cross-compilation sysroot — those binaries can be linked
    against a different libc, and running them under R segfaults: `utils` then fails to load, R
    continues without it, and `install.packages` and much else silently do not exist. Measured here,
    where the probe returned status 139 and the package check therefore reported the packages
    missing when they were installed and loadable.

    That failure mode is exactly the one §16b.2 says this gate must not have. A skip is the DEFAULT
    state of this module, so anything that makes the probe fail for an unrelated reason turns the
    oracle off permanently and silently, and the suite stays green while proving nothing. Preferring
    the system tools costs nothing when there is no second toolchain and is the whole difference when
    there is.

    **`R_LIBS_USER` and `R_LIBS_SITE` are passed through when set**, so an oracle installed into a
    user library is found, and so a site library full of packages built under an older R can be
    excluded by the caller without editing this file.

    **And when `R_LIBS_USER` is NOT set, a candidate library that exists is searched for** (§16b).
    Passing a variable through does nothing when the variable is unset, and R's own fallback points
    at a directory that does not exist here — so without this the gate stays shut on a machine where
    the packages are installed. An explicitly set variable still wins: a caller pointing at a
    different library is not overruled by a guess.
    """
    environment = {
        "PATH": "/usr/bin:/bin:" + os.environ.get("PATH", ""),
        "HOME": os.environ.get("HOME", ""),
        "LANG": "C",
    }
    for passthrough in ("R_LIBS_USER", "R_LIBS_SITE", "R_LIBS"):
        if passthrough in os.environ:
            environment[passthrough] = os.environ[passthrough]
    if "R_LIBS_USER" not in environment:
        found = _r_user_library()
        if found is not None:
            environment["R_LIBS_USER"] = found
    return environment


R_ENV = _r_environment()


def _r_has(packages: tuple[str, ...]) -> bool:
    """True only if every package LOADS. A package that installs but cannot load is not an oracle."""
    if RSCRIPT is None:
        return False
    probe = ";".join(f'if(!requireNamespace("{p}",quietly=TRUE)) quit(status=1)' for p in packages)
    try:
        return subprocess.run([RSCRIPT, "--vanilla", "-e", probe], env=R_ENV,
                              capture_output=True, timeout=180).returncode == 0
    except (OSError, subprocess.SubprocessError):
        return False


reference_r = pytest.mark.skipif(
    RSCRIPT is None or not _r_has(R_PACKAGES),
    reason=f"R oracle: needs Rscript and {', '.join(R_PACKAGES)}; see spec §16b.2. Searched for a "
           f"user library in {', '.join(str(p) for p in _R_USER_LIB_CANDIDATES)} and used "
           f"R_LIBS_USER={R_ENV.get('R_LIBS_USER', 'unset')}. Install into one of those and it is "
           "found without editing this file; if R reports that install.packages does not exist, "
           "that is the PATH collision _r_environment documents, not a broken R.")

DATA_GATED = pytest.mark.skipif(
    not config.DATA_XLSX.exists(),
    reason="the private workbook is gitignored and absent from this checkout")


def _run(script: str, source: Path, destination: Path, *extra: str) -> str:
    """Every R invocation goes through here, and every one gets `R_ENV` — including this one.

    The gate's probe and the scripts it gates must share an environment. An earlier version passed
    `env=R_ENV` in `_r_has` alone, so the probe found both packages, the tests un-skipped, and then
    every script died in R's startup with `could not find function "read.csv"` — `utils` unloadable
    for the PATH reason `_r_environment` documents. Loud rather than silent, but it is the same
    defect in the half of the pair that does the work.
    """
    done = subprocess.run(
        [RSCRIPT, "--vanilla", str(REFERENCE / script), str(source), str(destination), *extra],
        env=R_ENV, capture_output=True, text=True, timeout=600)
    assert done.returncode == 0, done.stderr
    return done.stdout


def _write(frame: pd.DataFrame, path: Path) -> Path:
    """Full precision across the boundary: 6 significant figures could not distinguish 1e-15
    agreement from 1e-7 agreement, and the first is the claim being made."""
    frame.to_csv(path, index=False, float_format="%.17g")
    return path


def cohort_design():
    """The workbook's design and response, built by `model.design` so R cannot disagree about it."""
    import cohort
    import data
    import derive
    import eligibility
    df, audit = data.load(data.WORKBOOK)
    df = cohort.build(eligibility.classify(derive.derive(df, audit), audit), audit)
    mask = model.complete_cases(df, config.PS_COVARIATES)
    X, _ = model.design(df.loc[mask], config.PS_COVARIATES)
    return X, df.loc[mask, config.TREATMENT].to_numpy(dtype=float)


# --- the two comparisons -------------------------------------------------------------------------

@reference_r
@DATA_GATED
def test_logistf_agrees_with_our_kernel_on_the_workbook_design(tmp_path):
    """DoD-16's second half: whether `logistf` and we agree on the coefficients.

    The tolerance is deliberately not machine precision, and the reason is a structural difference
    already established by reading the source (§18g): `logistf` stops when xconv AND gconv AND lconv
    are ALL met, ours stops when lconv OR gconv is met and never tests the coefficient step. On a
    penalised surface that §5.3 describes as genuinely flat, two different stopping rules stop in
    different places — so the honest comparison is on the penalised log-likelihood, with a loose
    coefficient bound beside it. The Python port of `logistf` differs from us by 3.66e-06 here for
    exactly this reason while agreeing on the objective to 2e-10, which is what sets these tolerances.
    """
    X, y = cohort_design()
    frame = X.copy()
    frame["y"] = y
    out = tmp_path / "coefficients.csv"
    version = _run("firth_logistf.R", _write(frame, tmp_path / "design.csv"), out,
                   str(config.FIRTH_MAX_ITER), str(config.FIRTH_MAX_HALVINGS),
                   str(config.FIRTH_MAX_STEP), str(config.FIRTH_TOL),
                   str(config.FIRTH_SCORE_TOL), str(config.FIRTH_TOL))
    assert version.strip(), "the R package version must be captured, not assumed"

    theirs = pd.read_csv(out)["coefficient"].to_numpy(dtype=float)
    ours = model.firth(X, y)
    Xc = np.column_stack([np.ones(len(y)), X.to_numpy(dtype=float)])
    assert abs(model._penalised_loglik(Xc, y, ours.beta)
               - model._penalised_loglik(Xc, y, theirs)) < 1e-6
    assert float(np.max(np.abs(ours.beta - theirs))) < 1e-3


@reference_r
@DATA_GATED
def test_psweight_agrees_on_e_and_w_before_the_ess_is_compared(tmp_path):
    """DoD-16's first half. **Compare `e` and `w` FIRST**, and the ESS separately.

    `PSweight` owns the whole ATO pipeline, so it may report an effective sample size that is not
    Kish's, or report one after its own normalisation of the weights. A convention difference in the
    third quantity must not be read as a weighting error in the first two, which is why this test
    asserts on `e` and `w` and merely RECORDS what the ESS came back as.

    `PSweight` fits its own unpenalised propensity model, so exact agreement with a Firth score is not
    expected and is not asserted: what is asserted is that the two scores describe the same population
    — rank correlation and the ATO weight identity — which is what a check on the *estimand* rather
    than on the *estimator* can establish.
    """
    X, y = cohort_design()
    frame = X.copy()
    frame["a"] = y
    out = tmp_path / "weights.csv"
    reported = _run("ato_psweight.R", _write(frame, tmp_path / "cohort.csv"), out)
    assert "PSweight" in reported

    theirs = pd.read_csv(out)
    ours = model.firth(X, y)
    assert np.all((theirs["e"].to_numpy() > 0.0) & (theirs["e"].to_numpy() < 1.0))
    # the ATO weight identity, which is the estimand and not the estimator
    expected = np.where(y == 1.0, 1.0 - theirs["e"].to_numpy(), theirs["e"].to_numpy())
    assert np.allclose(theirs["w"].to_numpy(), expected, atol=1e-12)
    # and the two scores order the cohort the same way, which is what a different fitter should preserve
    assert float(pd.Series(ours.p).corr(theirs["e"], method="spearman")) > 0.9

    # §16b.2's LAST open question, and it can only be answered by running the package: is the ESS it
    # REPORTS Kish's? Measured: yes, to every digit it prints. So there is no convention difference to
    # caveat, and the gap between its ESS and ours is entirely the propensity model — the estimator
    # [§7] prescribes — rather than the formula or a normalisation.
    #
    # Note WHY normalisation could never have produced one: Kish is scale-invariant, since
    # (Σcw)²/Σ(cw)² = (Σw)²/Σw². Only a different formula could differ, and it is not a different one.
    numbers = dict(zip(reported.split()[2::2], (float(v) for v in reported.split()[3::2])))
    for code, key in ((1, "reported_ess_treated"), (0, "reported_ess_control")):
        arm = theirs["w"].to_numpy()[y == float(code)]
        assert numbers[key] == pytest.approx(propensity.ess(arm), rel=1e-6)


# --- Stage 7 §16b — the balance oracle ---------------------------------------------------------------
#
# `PSweight` reports a per-covariate balance table on ANY supplied propensity score, which makes it an
# independent check on §5's arithmetic and on §4's covariate set at once. The score handed across is
# ours, so a disagreement is about the balance FORMULA and never about the propensity model.
#
# **The convention differs by a DEFAULT, not by a limitation.** `summary.SumStat` takes a
# `weighted.var` argument, and `weighted.var = FALSE` IS [§9]'s denominator:
#
#     PSweight, weighted.var = TRUE   |Δmean_w| / sqrt((SD_w,0² + SD_w,1²)/2)   WEIGHTED SDs
#     PSweight, weighted.var = FALSE  |Δmean_w| / sqrt((SD_0²   + SD_1²  )/2)   [§9]'s, and ours
#
# so the oracle is run BOTH ways: FALSE for a direct comparison against §5's weighted column, which is
# strictly stronger than reconstructing a numerator, and TRUE for the recorded convention difference.
# PSweight reports the ABSOLUTE value where §5 keeps the sign.


def imperfect_score(n: int = 80):
    """A synthetic design with a DELIBERATELY imperfect score, so the weighted numerator is non-zero.

    Synthetic rather than the workbook, so this comparison runs on a checkout with no `data/`: what
    is being checked is an arithmetic convention, and patient data adds nothing to that.
    """
    rng = np.random.default_rng(7)
    x1 = rng.normal(size=n)
    x2 = rng.binomial(1, 0.4, size=n).astype(float)
    a = rng.binomial(1, 1.0 / (1.0 + np.exp(-(0.8 * x1 + 0.4 * x2)))).astype(float)
    e = 1.0 / (1.0 + np.exp(-(0.4 * x1 + 0.2 * x2)))
    frame = pd.DataFrame({"x1": x1, "x2": x2, "a": a, "e": e})
    return frame, np.where(a == 1.0, 1.0 - e, e)


def psweight_balance(tmp_path):
    frame, w = imperfect_score()
    out = tmp_path / "balance.csv"
    reported = _run("balance_psweight.R", _write(frame, tmp_path / "cohort.csv"), out)
    assert "PSweight" in reported
    table = pd.read_csv(out).set_index(["convention", "covariate"])
    return frame, w, table


@reference_r
def test_the_unweighted_column_agrees_with_psweight_EXACTLY(tmp_path):
    """With w ≡ 1 the weighted SD is the unweighted SD, so the two conventions coincide — and the
    oracle therefore validates §5's arithmetic outright on half the table.
    """
    frame, _, table = psweight_balance(tmp_path)
    a = frame["a"].to_numpy(dtype=float)
    for covariate in ("x1", "x2"):
        ours = balance.smd(frame[covariate].to_numpy(dtype=float), a, np.ones(len(frame)))
        theirs = table.loc[("before_unweighted_var", covariate), "smd"]
        assert abs(ours) == pytest.approx(theirs, rel=1e-9)


@reference_r
def test_the_WEIGHTED_column_agrees_under_weighted_var_FALSE(tmp_path):
    """The strong assertion: [§9]'s denominator read straight out of the package, compared against
    §5's weighted column rather than reconstructed from a numerator."""
    frame, w, table = psweight_balance(tmp_path)
    a = frame["a"].to_numpy(dtype=float)
    for covariate in ("x1", "x2"):
        ours = balance.smd(frame[covariate].to_numpy(dtype=float), a, w)
        theirs = table.loc[("after_unweighted_var", covariate), "smd"]
        assert abs(ours) == pytest.approx(theirs, rel=1e-6)
        assert theirs >= 0.0                      # PSweight reports the absolute value


@reference_r
def test_psweights_ddof_is_ONE_and_a_ddof_zero_denominator_would_be_visible(tmp_path):
    """Asserted rather than assumed (§16b): a ddof = 0 convention differs by sqrt(n/(n-1)), which at
    n = 80 is 0.6% — visible at six significant figures and invisible to the eye, and enough to make
    an "agrees exactly" claim false in the fourth digit.
    """
    frame, w, table = psweight_balance(tmp_path)
    a = frame["a"].to_numpy(dtype=float)
    for covariate in ("x1", "x2"):
        x = frame[covariate].to_numpy(dtype=float)
        row = table.loc[("after_unweighted_var", covariate)]
        numerator = abs(row["mean_1"] - row["mean_0"])          # their own weighted means
        ddof_one = balance._pooled_sd(x, a)
        ddof_zero = float(np.sqrt((x[a == 1.0].var(ddof=0) + x[a == 0.0].var(ddof=0)) / 2.0))

        # their SMD is our numerator over OUR ddof = 1 denominator, and over no other
        assert numerator / ddof_one == pytest.approx(row["smd"], rel=1e-6)
        assert numerator / ddof_zero != pytest.approx(row["smd"], rel=1e-4)
        # the size of the difference the eye would miss: a percent or so — it is ~1/(2k) in the arm
        # size k, so 1.3% on these two arms of forty and ~0.5% on the workbook's ninety-two.
        assert 1e-3 < abs(ddof_zero - ddof_one) / ddof_one < 5e-2


@reference_r
def test_the_convention_difference_under_weighted_var_TRUE_is_RECONSTRUCTIBLE(tmp_path):
    """The denominators differ BY DECLARATION and the difference is recoverable, which is Stage 6
    §16b.2's "compare `e` and `w` first and the ESS separately" applied to a third quantity:
    `|ours| x our_pooled_unweighted_SD == theirs x their_pooled_weighted_SD`.

    [§9] prescribes the unweighted denominator and gives the reason — "so the yardstick does not
    move" — and PSweight under its own default recomputes the SD under each weighting scheme, so its
    before and after columns are divided by two different numbers. A legitimate convention, and not
    this SAP's. The oracle settles which convention each side uses, not which is better.
    """
    frame, w, table = psweight_balance(tmp_path)
    a = frame["a"].to_numpy(dtype=float)
    for covariate in ("x1", "x2"):
        x = frame[covariate].to_numpy(dtype=float)
        row = table.loc[("after_weighted_var", covariate)]
        theirs_sd = float(np.sqrt((row["sd_0"] ** 2 + row["sd_1"] ** 2) / 2.0))
        ours = balance.smd(x, a, w)
        assert abs(ours) * balance._pooled_sd(x, a) == pytest.approx(
            row["smd"] * theirs_sd, rel=1e-6)
        # and the ratios themselves do NOT agree, which is what makes this a convention difference
        assert abs(ours) != pytest.approx(row["smd"], rel=1e-4)


@reference_r
def test_the_gate_is_measured_OPEN_rather_than_assumed_open(monkeypatch):
    """§16b's closing assertion, and the state the oracle was silently skipping in.

    The parent shell's `R_LIBS_USER` is removed, `_r_environment` is rebuilt, and `_r_has` is asked
    whether both packages load under it. A gate that can only be opened by the caller exporting a
    variable is a gate that is shut by default, on the machine where the packages are installed.
    """
    monkeypatch.delenv("R_LIBS_USER", raising=False)
    environment = _r_environment()
    resolved = environment.get("R_LIBS_USER")
    assert resolved is not None and Path(resolved).is_dir()

    probe = ";".join(f'if(!requireNamespace("{p}",quietly=TRUE)) quit(status=1)'
                     for p in R_PACKAGES)
    done = subprocess.run([RSCRIPT, "--vanilla", "-e", probe], env=environment,
                          capture_output=True, timeout=180)
    assert done.returncode == 0, done.stderr.decode()


@reference_r
def test_the_r_package_versions_are_captured_and_non_empty(tmp_path):
    """An oracle that agreed at an unrecorded version is a claim with no date on it (§16b.2).

    Whoever opens this gate copies what this prints into the spec's §18, with the date it ran.
    """
    frame = pd.DataFrame({"x": np.arange(1.0, 11.0), "y": [0.0] * 5 + [1.0] * 5})
    version = _run("firth_logistf.R", _write(frame, tmp_path / "d.csv"), tmp_path / "c.csv",
                   str(config.FIRTH_MAX_ITER), str(config.FIRTH_MAX_HALVINGS),
                   str(config.FIRTH_MAX_STEP), str(config.FIRTH_TOL),
                   str(config.FIRTH_SCORE_TOL), str(config.FIRTH_TOL))
    assert "logistf" in version and any(ch.isdigit() for ch in version)


# --- the gate itself, which runs everywhere ----------------------------------------------------------

R_ONLY_NAMES = ("Rscript", "logistf", "PSweight")

# Every committed script, so a fourth one joins both hygiene checks below by being declared here.
R_SCRIPTS = ("firth_logistf.R", "ato_psweight.R", "balance_psweight.R", "polr_clm.R",
             "aug_psweight.R")


@pytest.mark.parametrize("name", R_ONLY_NAMES)
def test_R_is_confined_to_this_file(name):
    """§16b.2 and DoD-15: R's blast radius is this module and the two committed scripts.

    This is what keeps the oracle an oracle. R is where the methodologically authoritative
    implementations live, and it is declined as a RUNTIME dependency rather than as a reference — so
    the one thing that must stay true is that no estimation path can reach it, and no shipped module
    can so much as name it.

    The scan lives HERE rather than in `test_model.py` for a reason that is the rule itself: a scan
    must name the strings it searches for, so wherever it lives becomes a file that names them. Put it
    in the general test module and DoD-15's `grep` returns two files and a reader has to reason about
    which hit is real. Put it in the one file already licensed to know R exists, and the grep returns
    exactly that file.
    """
    offenders = [str(path.relative_to(PROJECT)) for path in _python_files()
                 if path.resolve() != Path(__file__).resolve()
                 and name in path.read_text(encoding="utf-8")]
    assert offenders == []


def _python_files() -> list[Path]:
    root = PROJECT
    found: list[Path] = []
    for dirpath, dirnames, filenames in os.walk(root):
        dirnames[:] = sorted(d for d in dirnames if d not in {".venv", "__pycache__"})
        found.extend(Path(dirpath) / n for n in sorted(filenames) if n.endswith(".py"))
    return found


def test_this_module_is_imported_by_nothing():
    offenders = [str(path.name) for path in _python_files()
                 if path.resolve() != Path(__file__).resolve()
                 and "test_reference_r" in path.read_text(encoding="utf-8")]
    assert offenders == []


def test_the_committed_r_scripts_exist_and_are_short_enough_to_read():
    # An oracle nobody can read is not evidence. Twenty-odd lines each: read CSV, fit, write
    # coefficients — no analysis logic and no covariate list, so the R side cannot disagree about the
    # SPECIFICATION while appearing to disagree about the ESTIMATOR.
    for script in R_SCRIPTS:
        path = REFERENCE / script
        assert path.exists(), path
        code = [line for line in path.read_text(encoding="utf-8").splitlines()
                if line.strip() and not line.strip().startswith("#")]
        assert len(code) <= 30, f"{script} has grown {len(code)} lines of code; keep it reviewable"


def test_the_scripts_name_no_covariate():
    # The design matrix arrives already built by `model.design`, which is what keeps invariant 4 and
    # the [§6] specification on the Python side of the boundary.
    for script in R_SCRIPTS:
        # Over the CODE only, and on whole words: the comments explain what the scripts are for, and
        # `packageVersion` contains "age".
        code = " ".join(line.split("#")[0]
                        for line in (REFERENCE / script).read_text(encoding="utf-8").splitlines())
        words = set(re.findall(r"[A-Za-z_][A-Za-z0-9_]*", code))
        assert words.isdisjoint(config.PS_COVARIATES), f"{script} names a covariate"


# --- Stage 8 §18b — `ordinal::clm` with weights, the ONE independent check on the weighted fit --------
#
# Stage 8's estimator carries five oracles. Four of them check us against ourselves — integer weights
# against row replication, `scipy` on our own objective, the analytic derivatives against central
# differences — or against an unweighted special case, `statsmodels.OrderedModel`, which §2 measured
# has no observation-weight support of any kind. **`clm` takes weights natively**, so it is the only
# INDEPENDENT implementation of the quantity Stage 8 computes, and it is therefore the one that
# matters most.
#
# It runs. An earlier draft of Stage 8 §18b recorded it as unrunnable on this machine, having measured
# `Rscript --vanilla` FROM A BARE SHELL — where R genuinely segfaults during its own startup,
# `system("uname -a")` returns status 139, `utils` and `stats` fail to load, and `requireNamespace`
# then returns FALSE for packages that are present. All of that is real and reproducible, and it is
# the PATH collision `_r_environment` above documents and repairs. The draft concluded that Stage 7's
# repair "does not help", which is the opposite of what is measured: through `R_ENV` it helps
# completely, and the fourteen R oracle tests of Stages 6 and 7 pass in this suite today. The evidence
# was one `pytest -q` away.
#
# THE GATE IS ITS OWN, and not `reference_r`. Stage 8 needs `ordinal`; Stages 6 and 7 need `logistf`
# and `PSweight`. Sharing one gate would skip this oracle on a machine that has `ordinal` and not
# `PSweight`, and Stage 6 §16b.2's whole argument is that a silently-skipping oracle is worse than no
# oracle.

ORDINAL_PACKAGES = ("ordinal",)


def _r_starts() -> bool:
    """True if R's own startup completes — which is a DIFFERENT question from whether a package is
    installed, and the skip reason below has to tell them apart.

    A skip reason that misattributes the cause sends the next reader to install a package that is
    already installed. The probe is `system()`, because that is what fails: R's startup runs
    `system("uname -a", intern = TRUE)` through `sh`, and on a machine carrying a second toolchain
    prefix those binaries can be linked against a different libc — measured status 139, a segfault
    inside R's forked child, after which `utils` never loads and `requireNamespace` returns FALSE for
    everything.
    """
    if RSCRIPT is None:
        return False
    try:
        done = subprocess.run(
            [RSCRIPT, "--vanilla", "-e", 'q(status = if (nzchar(R.version.string)) 0 else 1)'],
            env=R_ENV, capture_output=True, timeout=180)
    except (OSError, subprocess.SubprocessError):
        return False
    return done.returncode == 0


reference_r_ordinal = pytest.mark.skipif(
    RSCRIPT is None or not _r_has(ORDINAL_PACKAGES),
    reason=(
        "R oracle for Stage 8 §18b: needs Rscript and ordinal. "
        + ("Rscript is not on PATH." if RSCRIPT is None else
           ("`ordinal` did not load, and R itself started fine — install it into "
            f"R_LIBS_USER={R_ENV.get('R_LIBS_USER', 'unset')} and it is found without editing this "
            "file." if _r_starts() else
            "AND R ITSELF DID NOT START — this is not a missing package. R's startup runs "
            "system('uname -a') through sh; on a machine carrying a second toolchain prefix that "
            "call segfaults (status 139), utils never loads, and requireNamespace then returns "
            "FALSE for every package INCLUDING ones that are present. Do not install anything; see "
            "_r_environment above, which puts /usr/bin:/bin first on PATH for exactly this."))))


def _clm_weighted(y: np.ndarray, a: np.ndarray, w: np.ndarray, tmp_path: Path):
    """(thresholds, coefficient, version) from `ordinal::clm(ordered(y) ~ a, weights = w)`."""
    frame = pd.DataFrame({"y": y, "a": a, "w": w})
    out = tmp_path / "clm.csv"
    version = _run("polr_clm.R", _write(frame, tmp_path / "ordinal.csv"), out)
    table = pd.read_csv(out)
    thresholds = table.loc[table["kind"] == "threshold", "value"].to_numpy(dtype=float)
    coefficient = float(table.loc[table["kind"] == "coefficient", "value"].iloc[0])
    return thresholds, coefficient, version


@reference_r_ordinal
def test_clm_REPRODUCES_OUR_WEIGHTED_FIT_with_the_coefficient_NEGATED(tmp_path):
    """Stage 8 §18b's fifth oracle, and the strongest one this stage has.

    `clm` parametrises `logit P(Y <= k) = zeta_k - x'beta` where Stage 8 §5.1 uses [§14a]'s own
    `alpha_k + x'beta`, so **the coefficient comes back negated and the thresholds do not**. Both
    halves are asserted by name: a test asserting only that the two fits "agree" would pass on a
    comparison of the thresholds alone, which is exactly how the sign mistake survives a review — a
    reader sees six cutpoints agreeing to nine decimals and concludes the implementations agree.

    The third assertion is what makes the first two able to fail: `|ours - clm|` is asserted LARGE, so
    a frame on which both coefficients happened to be near zero cannot pass the sum test by having
    nothing to negate.

    Measured on `hand_ordinal()`: `|ours + clm|` 4.773e-12, `|ours - clm|` 1.739, and the thresholds
    matching to 7.177e-10 unnegated. This is the WEIGHTED fit validated against an independent
    implementation, which is the thing §18b says only this oracle can do.
    """
    from test_model import hand_ordinal
    y, a, w = hand_ordinal()
    thresholds, coefficient, version = _clm_weighted(y, a, w, tmp_path)
    assert version.strip(), "the R package version must be captured, not assumed"

    ours = model.polr(pd.DataFrame({config.TREATMENT: a}), y, w)
    assert len(thresholds) == len(ours.alpha) == 6

    assert abs(float(ours.beta[0]) + coefficient) < 1e-9        # NEGATED — measured 4.773e-12
    assert float(np.max(np.abs(ours.alpha - thresholds))) < 1e-8   # and NOT — measured 7.177e-10
    assert abs(float(ours.beta[0]) - coefficient) > 1e-2        # something to negate — measured 1.74


@reference_r_ordinal
def test_clm_agrees_on_a_SECOND_weighted_frame_so_the_agreement_is_not_one_fixtures(tmp_path):
    # `hand_ordinal()` has integer weights by construction (§14.3's replication oracle needs them).
    # This one has fractional weights, which is what the [§8] path actually passes — overlap weights
    # are strictly interior — so the oracle covers the weighting regime the stage uses.
    from test_model import replication_frame
    X, y, w = replication_frame(n=200, levels=4)
    a = X["x1"].to_numpy(dtype=float)
    fractional = w / 7.0
    thresholds, coefficient, _ = _clm_weighted(y, a, fractional, tmp_path)
    ours = model.polr(pd.DataFrame({"a": a}), y, fractional)
    assert abs(float(ours.beta[0]) + coefficient) < 1e-8
    assert float(np.max(np.abs(ours.alpha - thresholds))) < 1e-7
    assert abs(float(ours.beta[0])) > 1e-2


@reference_r_ordinal
def test_the_ordinal_gate_is_measured_OPEN_rather_than_assumed_open(monkeypatch):
    """Stage 6 §16b's closing assertion, for this stage's own package.

    The parent shell's `R_LIBS_USER` is removed, `_r_environment` is rebuilt, and `_r_has` is asked
    whether `ordinal` loads under it. A gate that can only be opened by the caller exporting a
    variable is a gate that is shut by default, on the machine where the package is installed — which
    is the state this oracle was recorded as being permanently in.
    """
    monkeypatch.delenv("R_LIBS_USER", raising=False)
    environment = _r_environment()
    probe = ";".join(f'if(!requireNamespace("{p}",quietly=TRUE)) quit(status=1)'
                     for p in ORDINAL_PACKAGES)
    done = subprocess.run([RSCRIPT, "--vanilla", "-e", probe], env=environment,
                          capture_output=True, timeout=180)
    assert done.returncode == 0, done.stderr.decode()


@reference_r_ordinal
def test_R_ITSELF_STARTS_under_R_ENV_which_is_the_claim_an_earlier_draft_got_backwards():
    # Recorded as a test rather than as a comment, because the failure it rules out is the one that
    # made a working oracle look unavailable for three days. From a BARE shell R does segfault here;
    # through `_r_environment` it does not.
    assert _r_starts()


def test_the_skip_reason_DISTINGUISHES_a_missing_package_from_an_R_THAT_DOES_NOT_START():
    """Stage 6 §16b.2's rule: a silently-skipping oracle is worse than no oracle, and a skip reason
    that misattributes the cause is a silently-skipping oracle with a misleading label.

    Asserted unconditionally — it is a property of the reason string, not of this machine.
    """
    reason = reference_r_ordinal.kwargs["reason"]
    assert "ordinal" in reason
    if RSCRIPT is not None and not _r_has(ORDINAL_PACKAGES) and not _r_starts():
        assert "R ITSELF DID NOT START" in reason
        assert "Do not install anything" in reason
    assert "Stage 8" in reason


# --- Stage 9 §19b — the augmented estimator, and why no package is its oracle --------------------------
#
# `PSweight` is the reference implementation of this SAP's exact WEIGHTING — overlap weights, written
# by the authors of the method — so it is the right oracle for the weighted marginal proportions and
# their difference, which is the same check Stage 6 §16b.2 made one function over.
#
# **It is NOT an oracle for the augmented estimate, and that is a fact about the ESTIMAND rather than
# about the package.** Every maintained augmented routine on offer targets the ATE or the ATT: it
# tilts by `1/e(1−e)`-style weights, or by `e/(1−e)`, and [§8] tilts by `h = e(1−e)`. A routine
# computing a different weighted average treatment effect is not a check on this one — it is a
# different quantity that would disagree for the right reason, which is the worst kind of oracle. So
# §8.1's formula is assembled TERM BY TERM in R and compared with the Python transcription.
#
# What that establishes and what it does not: agreement between two hand transcriptions in two
# languages shows the two TRANSCRIPTIONS agree. It does not show the formula is [§8]'s. The tilt tests
# in `test_outcome.py` are what check that, because only this repository knows which tilt [§7] means —
# and they are in Python for exactly that reason.


def augmented_frame():
    """The frame the oracle is run on: `golden_frame()` with `w`, `m1` and `m0` already computed.

    The GOLDEN frame and not the workbook, deliberately, and it is the one oracle in this file that
    is not data-gated. Stage 9 §4.3's rule is that no WEIGHTED quantity computed on patient data
    leaves the gitignored log — and an R oracle writes a CSV into a temporary directory and prints a
    version banner, which is a second place a weighted proportion would exist. The golden frame is
    synthetic, carries no patient data, and is the frame every other pin in Stage 9 is measured on.
    """
    import outcome
    from fixtures_stage9 import golden_arrays, golden_frame
    frame = golden_frame()
    y, a, w, e = golden_arrays(frame)
    X, _ = model.design(frame, config.outcome_model_covariates("tici_2b_3"))
    X = X.copy()
    X.insert(0, config.TREATMENT, a)
    fit = model.firth(X, y)
    m1, m0 = outcome._counterfactuals(fit, X)
    return pd.DataFrame({"y": y, "a": a, "e": e, "w": w, "m1": m1, "m0": m0}), (y, a, w, e, m1, m0)


@reference_r
def test_psweight_agrees_on_the_weighted_marginal_proportions_and_the_risk_difference(tmp_path):
    """The first row. The score is SUPPLIED, so a disagreement is about the weighting formula and
    never about the propensity model — which is the whole reason Stage 7's balance oracle supplies
    one too.

    Asserted tightly: both sides are computing the same weighted mean of the same 0/1 vector over the
    same weights, so there is no stopping rule and no optimiser between them, and anything looser than
    machine precision would be hiding a real difference.
    """
    import outcome
    frame, (y, a, w, e, m1, m0) = augmented_frame()
    out = tmp_path / "augmented.csv"
    version = _run("aug_psweight.R", _write(frame, tmp_path / "golden.csv"), out)
    assert "PSweight" in version

    theirs = pd.read_csv(out)
    share = outcome.weighted_proportion(y, a, w)
    assert theirs["p1"][0] == pytest.approx(share[1], abs=1e-12)
    assert theirs["p0"][0] == pytest.approx(share[0], abs=1e-12)
    assert theirs["rd"][0] == pytest.approx(outcome.weighted_rd(y, a, w), abs=1e-12)


@reference_r
def test_the_augmented_estimate_agrees_with_a_HAND_ASSEMBLED_R_COMPUTATION(tmp_path):
    """The second row, and the one that matters.

    Not against a package routine, for §19b's reason: every maintained AIPW implementation targets a
    different estimand. This is a cross-language check on the ARITHMETIC of §8.1 — three terms, two
    arm denominators and one tilt denominator — and explicitly not on the estimand.
    """
    import outcome
    frame, (y, a, w, e, m1, m0) = augmented_frame()
    out = tmp_path / "augmented.csv"
    _run("aug_psweight.R", _write(frame, tmp_path / "golden.csv"), out)

    theirs = float(pd.read_csv(out)["tau"][0])
    ours = outcome.augmented_rd(y, a, w, outcome._tilt(e), m1, m0)
    assert theirs == pytest.approx(ours, abs=1e-12)
    # and it is the PINNED value, so the oracle checks the transcription against the same number
    # every other Stage 9 test checks it against rather than against whatever we happen to compute
    assert ours == pytest.approx(0.2233300473, abs=1e-9)


def test_the_stage_9_oracle_does_NOT_ask_a_package_for_the_AUGMENTED_estimate():
    """By scan, because it is the one thing about this script that a passing comparison cannot show.

    If the augmented value came from a package routine, both tests above would still pass or fail
    together and nothing would record that the number came from an estimator targeting the ATE. The
    script's `tau` is built from `sum(...)` over the frame's own columns and from nothing else.
    """
    code = " ".join(line.split("#")[0]
                    for line in (REFERENCE / "aug_psweight.R").read_text(
                        encoding="utf-8").splitlines())
    tau_line = code[code.index("tau <-"):code.index("write.csv")]
    assert "::" not in tau_line, "the augmented term must not come from a package"
    assert tau_line.count("sum(") == 6
