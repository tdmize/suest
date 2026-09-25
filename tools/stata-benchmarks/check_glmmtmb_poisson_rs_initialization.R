# Run from tools/stata-benchmarks; uses the unchanged benchmark CSV.
suppressPackageStartupMessages(library(glmmTMB)) # Laplace mixed-model fits
dat <- read.csv("suest_r_glmmtmb_poisson_rs_crosslang_benchmark.csv")
dat$id <- factor(dat$id)
outcomes <- c("y1", "y2", "y1_partial", "y2_partial", "y1_left", "y2_right")
checks <- lapply(outcomes, function(outcome) {
  d <- dat[!is.na(dat[[outcome]]), ]
  formula <- reformulate(c("x", "z", "(1 + x | id)"), outcome)
  base <- glmmTMB(formula, data = d, family = poisson())
  # Fixed-only coefficients and unit random variances provide a distinct start.
  beta <- coef(glm(reformulate(c("x", "z"), outcome), data = d, family = poisson()))
  alternate <- glmmTMB(formula, data = d, family = poisson(),
    start = list(beta = unname(beta), theta = c(0, 0, 0)),
    control = glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS"),
      optCtrl = list(reltol = 1e-12, maxit = 1000)))
  S <- as.matrix(VarCorr(base)$cond$id)
  result <- data.frame(outcome = outcome, n = nobs(base),
    convergence = base$fit$convergence, pd_hessian = base$sdr$pdHess,
    alternate_convergence = alternate$fit$convergence,
    alternate_pd_hessian = alternate$sdr$pdHess,
    logLik = as.numeric(logLik(base)),
    abs_logLik_difference = abs(as.numeric(logLik(base)-logLik(alternate))),
    max_abs_parameter_difference = max(abs(base$fit$par-alternate$fit$par)),
    max_abs_gradient = max(abs(base$obj$gr(base$fit$par))),
    alternate_max_abs_gradient = max(abs(alternate$obj$gr(alternate$fit$par))),
    sd_intercept = sqrt(S[1, 1]), sd_slope = sqrt(S[2, 2]),
    correlation = cov2cor(S)[1, 2], min_covariance_eigenvalue = min(eigen(S)$values))
  stopifnot(result$convergence == 0, result$alternate_convergence == 0,
    result$pd_hessian, result$alternate_pd_hessian,
    result$abs_logLik_difference < 1e-5,
    result$max_abs_parameter_difference < 5e-4,
    result$max_abs_gradient < .005, result$alternate_max_abs_gradient < .005)
  result
})
checks <- do.call(rbind, checks)
print(checks, digits = 10, row.names = FALSE)
write.csv(checks, "glmmtmb_poisson_rs_initialization_checks.csv", row.names = FALSE)
cat("R_VERSION=", R.version.string, "\nGLMMTMB_VERSION=", as.character(packageVersion("glmmTMB")), "\n", sep = "")
cat("POISSON_RS_R_INITIALIZATION_CHECKS_COMPLETE=1\n")
