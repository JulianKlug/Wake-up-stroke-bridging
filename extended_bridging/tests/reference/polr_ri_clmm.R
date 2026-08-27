# Stage 12 §12.7 — the ordinal::clmm oracle. Read a frame of (y, x1, x2, group), ask the package for
# a RANDOM-INTERCEPT proportional-odds fit, write the coefficients, the thresholds and the random-
# effect SD back.
#
# `clmm` is in the same package as `polr_clm.R`'s `clm`, so NO NEW R DEPENDENCY IS ADDED and the
# existing `ORDINAL_PACKAGES = ("ordinal",)` gate covers this file too.
#
# **THIS IS THE ONLY INDEPENDENT IMPLEMENTATION OF THE QUANTITY `model.polr_ri` COMPUTES.** Stage 12
# §12 measured that nothing in the Python environment can fit it: `statsmodels` has `MixedLM` for the
# linear case and Bayesian mixed GLMs for binomial and Poisson, and NO ordinal mixed model at any
# API; `scipy` is test-only by policy. So `polr_ri` is otherwise verified only against itself — its
# own analytic gradient, its node-count stability, and its collapse to `model.polr` at sigma = 0 —
# and this is the one comparison that can fail for the reason the estimator would be wrong.
#
# THE ASYMMETRY IS THREE-WAY AND THE CALLER ASSERTS ALL THREE PARTS.  `clmm` parametrises
# `logit P(Y <= k) = zeta_k - x'beta + b_c`, where Stage 12 §12.1 uses [§14a]'s own
# `alpha_k + x'beta + b_c`. So:
#
#     the COEFFICIENTS come back NEGATED
#     the THRESHOLDS do NOT
#     `sigma` does NOT either, because it is a SCALE and not a location
#
# A caller asserting that "the fits agree" would pass on the thresholds alone, which is exactly how
# the sign mistake survives a review — `polr_clm.R`'s banner makes the same point for two of the
# three, and the third is this file's own addition.
#
# `nAGQ` IS PASSED AND IS NOT LEFT TO DEFAULT.  `clmm`'s default is `nAGQ = 1`, the Laplace
# approximation, and Stage 12 §12.2's table shows that is a DIFFERENT OBJECTIVE — comparing against
# it would measure the approximation rather than the maximiser. It is passed as POLR_RI_NODES from
# the caller, so the two sides integrate the same likelihood.
#
# No covariate list and no analysis logic here, exactly as `polr_clm.R`, `firth_logistf.R`,
# `ato_psweight.R` and `balance_psweight.R` are shaped: the columns arrive already chosen, so this
# side cannot disagree about the SPECIFICATION while appearing to disagree about the ESTIMATOR.

args <- commandArgs(trailingOnly = TRUE)
options(digits = 17)

frame <- read.csv(args[1])
n_agq <- as.integer(args[3])

frame$group <- factor(frame$group)

fitted <- ordinal::clmm(
  ordered(y) ~ x1 + x2 + (1 | group),
  data = frame, link = "logit", nAGQ = n_agq)

thresholds <- fitted$alpha
coefficients <- fitted$beta

# THE RANDOM-EFFECT SD IS IN `$ST` AND **NOT** IN `$stDev`, WHICH IS NULL ON A `clmm` OBJECT.
# Measured: `names(fitted)` carries `ST` and no `stDev`, and `fitted$ST` is a LIST keyed by the
# grouping factor whose single entry is a 1x1 matrix — the `lme4`-style "ST" parametrisation, which
# for a lone random intercept IS the standard deviation. Verified against a frame generated at
# sigma = 0.8: `ST$group[1, 1]` comes back 0.826.
#
# An earlier version of this file read `fitted$stDev`, got `numeric(0)`, and made `data.frame()` fail
# on a length mismatch rather than on the comparison — so the oracle failed for a reason that had
# nothing to do with the estimator. It is recorded because `clm` and `clmm` differ here and the
# similarity of the two objects invites the assumption that they do not.
#
# It is written out under its own `kind` so the caller compares it at its own tolerance: a variance
# parameter at this cluster count is the least well determined thing in the fit, and pinning it to the
# coefficients' tolerance would be pinning the wrong number.
stopifnot(is.list(fitted$ST), length(fitted$ST) == 1L)
sigma <- as.numeric(fitted$ST[[1]][1, 1])

write.csv(
  data.frame(
    parameter = c(names(thresholds), names(coefficients), "sigma"),
    kind = c(rep("threshold", length(thresholds)),
             rep("coefficient", length(coefficients)),
             "sigma"),
    value = c(unname(thresholds), unname(coefficients), sigma),
    row.names = NULL),
  args[2], row.names = FALSE)

cat("ordinal", as.character(packageVersion("ordinal")), "nAGQ", n_agq, "\n")
