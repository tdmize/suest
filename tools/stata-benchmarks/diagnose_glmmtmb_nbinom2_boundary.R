# Run from tools/stata-benchmarks using the retained revision-1 CSV.
options(digits = 12)
dat <- read.csv("suest_r_glmmtmb_nbinom2_ri_crosslang_benchmark.csv")
dat$id <- factor(dat$id)
right <- droplevels(dat[as.integer(dat$id) > 40, ])
fit <- function(data, ...) glmmTMB::glmmTMB(y2 ~ x + z + (1 | id), data = data,
  family = glmmTMB::nbinom2("log"), ...)
m <- fit(right)
cat("ORIGINAL RIGHT-HALF FIT\n")
print(m$fit$par)
print(c(convergence = m$fit$convergence, pdHess = m$sdr$pdHess,
  variance = exp(2*m$fit$par["theta"]), logLik = as.numeric(logLik(m)),
  log_sd_se = sqrt(m$sdr$cov.fixed["theta", "theta"])))
cat("MULTIPLE STARTS\n")
for (sd in c(0.05, 0.3, 0.7, 1)) {
  f <- suppressWarnings(fit(right, start = list(theta = log(sd))))
  print(c(start_sd = sd, sd = exp(f$fit$par["theta"]),
    logLik = as.numeric(logLik(f)), convergence = f$fit$convergence, pdHess = f$sdr$pdHess))
}
cat("FIXED-VARIANCE PROFILE; OTHER PARAMETERS RE-ESTIMATED\n")
for (v in c(1e-12, 1e-8, 1e-6, 1e-4, .001, .01, .05, .1, .2)) {
  f <- fit(right, start = list(theta = .5*log(v)), map = list(theta = factor(NA)))
  print(c(variance = v, logLik = as.numeric(logLik(f))))
}
pooled <- MASS::glm.nb(y2 ~ x + z, data = right)
cat("NO RANDOM EFFECT NB2\n")
print(c(logLik = as.numeric(logLik(pooled)), log_phi = log(pooled$theta)))
cat("ALTERNATING-GROUP DISJOINT DESIGN\n")
for (outcome in c("y1", "y2")) {
  keep <- as.integer(dat$id) %% 2 == if (outcome == "y1") 1 else 0
  d <- droplevels(dat[keep, ])
  f <- glmmTMB::glmmTMB(reformulate(c("x", "z", "(1 | id)"), outcome), data = d,
    family = glmmTMB::nbinom2("log"))
  cat(outcome, "\n")
  print(c(variance = exp(2*f$fit$par["theta"]), logLik = as.numeric(logLik(f)),
    convergence = f$fit$convergence, pdHess = f$sdr$pdHess,
    log_sd_se = sqrt(f$sdr$cov.fixed["theta", "theta"])))
}
