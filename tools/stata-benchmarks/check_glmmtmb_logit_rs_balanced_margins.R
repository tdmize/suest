# Evaluate independent integrated margins at the returned Stata component points.
source("../../tests/testthat/helper-glmmtmb-logit-rs.R")
f <- readRDS("../../tests/testthat/fixtures/glmmtmb-logit-rs-stata-balanced-v1.rds")
result <- do.call(rbind, lapply(names(f$components), function(y) {
  s <- f$components[[y]]; a <- glmmtmb_logit_rs_native_margins(s$coefficients, s$vcov, f$data)
  returned <- vapply(s$native_margins, `[[`, numeric(1), "estimate")
  variance <- vapply(s$native_margins, `[[`, numeric(1), "variance")
  data.frame(outcome = y, quantity = names(a$estimates),
    returned_estimate = returned, independent_estimate = a$estimates,
    estimate_difference = returned-a$estimates,
    returned_variance = variance, independent_variance = diag(a$vcov),
    relative_se_difference = sqrt(variance/diag(a$vcov))-1)
}))
print(result, digits = 12, row.names = FALSE)
args <- commandArgs(trailingOnly = TRUE)
if (length(args)) write.csv(result, args[1L], row.names = FALSE)
stopifnot(max(abs(result$estimate_difference)) < 1e-8,
  max(abs(result$relative_se_difference)) < 1e-4)
cat("LOGIT_RS_BALANCED_NATIVE_MARGINS_COMPLETE=1\n")
