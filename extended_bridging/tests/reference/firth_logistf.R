# Stage 6 §16b.2 — the logistf oracle. Read a design, fit, write coefficients. Nothing else.
#
# No analysis logic and no covariate list: the design matrix arrives ALREADY BUILT by model.design,
# so this side cannot disagree about the specification while appearing to disagree about the
# estimator. The response is the last column; every other column is a design column.
#
# The control parameters are passed EXPLICITLY rather than left to logistf.control()'s defaults,
# because what makes the comparison like-for-like is the explicit pass and not a default two
# implementations happen to share [§16b.2]. Note that logistf's own defaults differ from ours in
# every one of maxit, maxhs, lconv and gconv, and that it additionally tests xconv, which our rule
# deliberately excludes — so a like-for-like comparison is only possible for the tolerances, not for
# the stopping RULE, which is a structural difference recorded in the spec's §18g.

args <- commandArgs(trailingOnly = TRUE)
options(digits = 17)                       # full precision back across the CSV boundary

frame <- read.csv(args[1])
y <- frame[[ncol(frame)]]
x <- as.matrix(frame[, -ncol(frame), drop = FALSE])

fit <- logistf::logistf(
  y ~ x, family = binomial(), pl = FALSE,
  control = logistf::logistf.control(
    maxit = as.integer(args[3]), maxhs = as.integer(args[4]), maxstep = as.numeric(args[5]),
    lconv = as.numeric(args[6]), gconv = as.numeric(args[7]), xconv = as.numeric(args[8])))

write.csv(data.frame(term = names(coef(fit)), coefficient = as.numeric(coef(fit))),
          args[2], row.names = FALSE)
cat("logistf", as.character(packageVersion("logistf")), "\n")
