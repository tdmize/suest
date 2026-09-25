# Run here after the existing six-model R reference has passed its preflight.
r <- readRDS("glmmtmb_nbinom2_rs_r_reference.rds")
outcomes <- c("y1_partial", "y2_partial")
stopifnot(all(outcomes %in% r$preflight$outcome),
  all(r$preflight$optimizer_loglik_gap < 1e-5), all(r$preflight$optimizer_parameter_gap < 1e-4))
seed <- lapply(r$components[outcomes], function(m) {
  b <- m$coefficients; S <- m$random_covariance
  unname(c(b[c("x", "z", "(Intercept)")], -b["log_phi"], S[2, 2], S[1, 1], S[1, 2]))
})
fmt <- function(x) paste(format(x, digits = 17, scientific = FALSE, trim = TRUE), collapse = ", ")
writeLines(c("* Generated from preflight-checked R Laplace fits; these are R starts, not Stata estimates.",
  paste0("matrix nb2rs_r_starts = (", fmt(seed[[1]]), " \\ ///"), paste0("    ", fmt(seed[[2]]), ")"),
  "matrix rownames nb2rs_r_starts = y1_partial y2_partial",
  "matrix colnames nb2rs_r_starts = x z _cons lnalpha var_slope var_intercept covariance",
  paste0("scalar nb2rs_r_", outcomes, "_ll = ", vapply(r$components[outcomes], function(m) fmt(m$logLik), character(1)))),
  "suest_r_glmmtmb_nbinom2_rs_partial_starts.do")
