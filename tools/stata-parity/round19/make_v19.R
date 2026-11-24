# R side of Stata round 19: median survival time predictions from survreg
# models (Weibull, exponential, lognormal, log-logistic), Stata's default
# prediction after streg. Uses batch4_fixture_v01.dta (from make_v18.R).
# Writes survreg_median_r_reference_v01.csv: for each model, the average
# median survival time and the average effect of x on it, with robust
# standard errors (Stata: the model's own vce(robust) fit, then margins), and
# the log likelihood.
# Usage: Rscript make_v19.R   (SUEST_DEV=<path> loads a development copy)
suppressPackageStartupMessages({
  if (nzchar(Sys.getenv("SUEST_DEV"))) pkgload::load_all(Sys.getenv("SUEST_DEV"), quiet = TRUE) else
    library(suest)
  library(haven); library(survival); library(marginaleffects)
})
d <- as.data.frame(read_dta("batch4_fixture_v01.dta"))
ctl <- survreg.control(rel.tolerance = 1e-13, maxiter = 100)
fit_one <- function(dist) survreg(Surv(tm, ev) ~ x + z, data = d, dist = dist, control = ctl)
models <- list(weibull = fit_one("weibull"), exponential = fit_one("exponential"),
  lognormal = fit_one("lognormal"), loglogistic = fit_one("loglogistic"))
out <- list()
row <- function(case, estimate, se) data.frame(case = case, estimate = estimate, se = se)
for (pair in list(c("weibull", "exponential"), c("lognormal", "loglogistic"))) {
  fit <- suest(models[[pair[1]]], models[[pair[2]]], model_names = pair, observation_id = "rid")
  p <- as.data.frame(avg_predictions(fit, newdata = d))
  s <- as.data.frame(avg_slopes(fit, variables = "x", newdata = d))
  for (m in pair) {
    out[[length(out) + 1]] <- row(paste0(m, "_median"), p$estimate[p$group == m], p$std.error[p$group == m])
    out[[length(out) + 1]] <- row(paste0(m, "_dydx_x"), s$estimate[s$group == m], s$std.error[s$group == m])
    out[[length(out) + 1]] <- row(paste0(m, "_ll"), models[[m]]$loglik[2], NA)
  }
}
ref <- do.call(rbind, out)
print(ref, digits = 12)
write.csv(ref, "survreg_median_r_reference_v01.csv", row.names = FALSE)
