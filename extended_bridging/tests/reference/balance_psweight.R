# Stage 7 §16b — the PSweight balance oracle. Read a frame and a SUPPLIED propensity score, ask the
# package for its per-covariate balance table, write it back.
#
# PSweight is the reference implementation of this SAP's exact estimand — overlap weights, written by
# the authors of the method — so its balance table is an independent check on §5's arithmetic and on
# §4's covariate set at once.
#
# The score is supplied rather than fitted here, which is the whole point: with `ps.estimate` the
# package weights the population OUR estimator defined, so any disagreement is about the balance
# formula and never about the propensity model.
#
# The summary is taken BOTH WAYS, because `weighted.var` is a default and not a limitation:
# `weighted.var = FALSE` is [§9]'s unweighted pooled denominator and is directly comparable with §5's
# weighted column; `weighted.var = TRUE` is the package's own convention, whose difference from [§9]'s
# is recoverable rather than mysterious. No covariate list here either — the columns arrive already
# chosen, so this side cannot disagree about the SPECIFICATION while appearing to disagree about the
# ESTIMATOR.

args <- commandArgs(trailingOnly = TRUE)
options(digits = 17)

frame <- read.csv(args[1])
covariates <- setdiff(names(frame), c("a", "e"))
score <- cbind(1 - frame$e, frame$e)
colnames(score) <- c("0", "1")

fitted <- PSweight::SumStat(zname = "a", xname = covariates, data = frame,
                            ps.estimate = score, weight = "overlap")
plain <- summary(fitted, weighted.var = FALSE, metric = "ASD")
their_own <- summary(fitted, weighted.var = TRUE, metric = "ASD")

rows <- function(table, convention) data.frame(
  covariate = rownames(table), convention = convention,
  mean_0 = table[, "Mean 0"], mean_1 = table[, "Mean 1"],
  sd_0 = table[, "Weighted SD 0"], sd_1 = table[, "Weighted SD 1"],
  smd = table[, "SMD"], row.names = NULL)

write.csv(rbind(rows(plain$unweighted, "before_unweighted_var"),
                rows(plain$overlap, "after_unweighted_var"),
                rows(their_own$overlap, "after_weighted_var")), args[2], row.names = FALSE)
cat("PSweight", as.character(packageVersion("PSweight")), "\n")
