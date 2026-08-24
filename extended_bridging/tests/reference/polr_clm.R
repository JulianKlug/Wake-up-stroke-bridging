# Stage 8 §18b — the ordinal::clm oracle. Read a frame of (y, a, w), ask the package for a WEIGHTED
# proportional-odds fit, write the coefficient and the thresholds back.
#
# `clm` takes observation weights natively, which makes it the only INDEPENDENT implementation of the
# quantity Stage 8 computes: the other four oracles check us against ourselves (row replication, scipy
# on our own objective, central differences) or against an unweighted special case (statsmodels).
# There is no maintained WEIGHTED proportional-odds implementation for Python — §2 measured that
# `OrderedModel` has no weight support of any kind — so this is the one comparison that can fail for
# the reason the estimator would be wrong.
#
# `clm` parametrises `logit P(Y <= k) = zeta_k - x'beta`, where Stage 8 §5.1 uses [§14a]'s own
# `alpha_k + x'beta`. So THE COEFFICIENT COMES BACK NEGATED AND THE THRESHOLDS DO NOT, and the caller
# asserts that asymmetry by name rather than asserting the two fits "agree" — which would pass on a
# comparison of the thresholds alone, and is exactly how the sign mistake survives a review.
#
# No covariate list and no analysis logic here, exactly as firth_logistf.R, ato_psweight.R and
# balance_psweight.R are shaped: the columns arrive already chosen, so this side cannot disagree about
# the SPECIFICATION while appearing to disagree about the ESTIMATOR.

args <- commandArgs(trailingOnly = TRUE)
options(digits = 17)

frame <- read.csv(args[1])

fitted <- ordinal::clm(ordered(y) ~ a, weights = w, data = frame, link = "logit")

thresholds <- fitted$alpha
write.csv(
  data.frame(
    parameter = c(names(thresholds), "a"),
    kind = c(rep("threshold", length(thresholds)), "coefficient"),
    value = c(unname(thresholds), unname(fitted$beta)),
    row.names = NULL),
  args[2], row.names = FALSE)

cat("ordinal", as.character(packageVersion("ordinal")), "\n")
