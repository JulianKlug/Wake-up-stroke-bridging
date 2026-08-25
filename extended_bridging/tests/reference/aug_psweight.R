# Stage 9 §19b — the augmented-estimator oracle. Read a frame, ask the package for the weighted
# marginal contrast, and assemble [§8]'s augmented formula BY HAND beside it.
#
# TWO ROWS ARE CHECKED AND THE SECOND IS THE ONE THAT MATTERS.
#
# 1. The weighted marginal proportions and their difference come from PSweight's overlap-weight
#    output on a SUPPLIED score, exactly as the Stage 7 balance oracle supplies one: with
#    `ps.estimate` the package weights the population OUR estimator defined, so a disagreement is
#    about the weighting formula and never about the propensity model. This is the same oracle
#    Stage 6 §16b.2 used one function over.
#
# 2. The augmented estimate is assembled TERM BY TERM here, and is NOT checked against an AIPW
#    routine from any package. That is a fact about the ESTIMAND and not about the package: every
#    maintained one targets the ATE or the ATT, so it tilts by 1/e(1-e)-style weights or by
#    e/(1-e) — not by h = e(1-e), which is what [§7] means by the overlap population. A routine
#    computing a different estimand is not an oracle for this one.
#
# WHAT THIS DOES NOT ESTABLISH, stated because the distinction is easy to lose: agreement with a
# hand-assembled formula in another language shows the two TRANSCRIPTIONS agree. It does not show the
# formula is [§8]'s. The Python suite's tilt tests are what check that, because only this repository
# knows which tilt [§7] means.
#
# No covariate list here either: `e`, `w`, `m1` and `m0` arrive already computed, so this side cannot
# disagree about the SPECIFICATION while appearing to disagree about the ESTIMATOR.

args <- commandArgs(trailingOnly = TRUE)
options(digits = 17)

frame <- read.csv(args[1])
score <- cbind(1 - frame$e, frame$e)
colnames(score) <- c("0", "1")

fitted <- PSweight::PSweight(ps.estimate = score, zname = "a", yname = "y",
                             data = frame, weight = "overlap")
p <- as.numeric(fitted$muhat)

# [§8]'s three terms, in the order [§8] writes them: two arm terms each normalised by its OWN arm's
# weight total, then one correction term normalised over BOTH arms — h is not arm-specific.
h <- frame$e * (1 - frame$e)
treated <- frame$a == 1
control <- frame$a == 0
tau <- sum(frame$w[treated] * (frame$y[treated] - frame$m1[treated])) / sum(frame$w[treated]) -
  sum(frame$w[control] * (frame$y[control] - frame$m0[control])) / sum(frame$w[control]) +
  sum(h * (frame$m1 - frame$m0)) / sum(h)

write.csv(data.frame(p0 = p[1], p1 = p[2], rd = p[2] - p[1], tau = tau), args[2], row.names = FALSE)
cat("PSweight", as.character(packageVersion("PSweight")), "\n")
