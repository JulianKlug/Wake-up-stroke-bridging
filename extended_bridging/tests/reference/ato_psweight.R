# Stage 6 §16b.2 — the PSweight oracle. Read a frame, fit the ATO, write e, w and the ESS.
#
# PSweight is the reference implementation of this SAP's exact estimand — overlap weights, written by
# the authors of the method — and it is the only available check on the WHOLE [§7] pipeline rather
# than on its first step.
#
# No covariate list here either: the columns arrive already chosen, and the formula is built from
# whatever the frame carries so that this side cannot disagree about the specification. The response
# column is named `a` by the caller; every other column is a design column.
#
# TWO MISMATCHES TO EXPECT, neither of which is a bug in either implementation [§16b.2]:
#   * PSweight's propensity model is its own (glm, unpenalised) unless one is supplied, so `e` may
#     differ from a Firth score for reasons that are not an error. Compare e and w FIRST.
#   * PSweight's effective sample size may not be Kish's, or may be reported after its own
#     normalisation of the weights. A convention difference in the third quantity must not be read as
#     a weighting error, which is why the ESS is written separately below.

args <- commandArgs(trailingOnly = TRUE)
options(digits = 17)

frame <- read.csv(args[1])
covariates <- setdiff(names(frame), "a")
formula <- as.formula(paste("a ~", paste(covariates, collapse = " + ")))

fitted <- PSweight::SumStat(ps.formula = formula, data = frame, weight = "overlap")
e <- as.numeric(fitted$propensity[, 2])
w <- ifelse(frame$a == 1, 1 - e, e)

# `fitted$ess` is what the PACKAGE reports, and the row named "1" is the treated arm. Written out
# beside our own Kish sum on the same weights so the test can compare the two CONVENTIONS rather than
# only the two weight vectors — which is the question §16b.2 could not answer by reading source.
reported <- fitted$ess[, "overlap"]
write.csv(data.frame(e = e, w = w), args[2], row.names = FALSE)
cat("PSweight", as.character(packageVersion("PSweight")),
    "reported_ess_treated", reported[["1"]], "reported_ess_control", reported[["0"]], "\n")
