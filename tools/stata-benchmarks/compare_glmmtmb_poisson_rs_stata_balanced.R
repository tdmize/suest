# Run here; optionally pass a focused fixture and output CSV path.
args <- commandArgs(trailingOnly = TRUE)
input <- if (length(args)) args[1] else "../../tests/testthat/fixtures/glmmtmb-poisson-rs-stata-balanced-v3.rds"
s <- readRDS(input)
reference <- readRDS("glmmtmb_poisson_rs_r_reference.rds"); r <- reference$systems$balanced
stopifnot(identical(names(s$coefficients), names(r$coefficients)), s$n_clusters == r$n_clusters)
comparison <- data.frame(parameter = names(r$coefficients),
  coefficient_difference = s$coefficients-r$coefficients,
  stata_se = sqrt(diag(s$vcov)), r_se = sqrt(diag(r$vcov)),
  relative_se_difference = sqrt(diag(s$vcov)/diag(r$vcov))-1)
print(comparison, digits = 10, row.names = FALSE)
cat("MAX_ABS_COVARIANCE_DIFFERENCE=", max(abs(s$vcov-r$vcov)), "\n")
cat("MAX_ABS_CROSS_COVARIANCE_DIFFERENCE=", max(abs(s$vcov[1:6, 7:12]-r$vcov[1:6, 7:12])), "\n")
cat("ABS_LOGLIK_DIFFERENCE=", abs(s$logLik-sum(vapply(reference$components[1:2], `[[`, numeric(1), "logLik"))), "\n")
if (!is.null(s$modelbased)) {
  V <- matrix(0, 12, 12)
  V[1:6, 1:6] <- reference$components$y1$vcov
  V[7:12, 7:12] <- reference$components$y2$vcov
  cat("MODELBASED_MAX_ABS_COVARIANCE_DIFFERENCE=", max(abs(s$modelbased-V)), "\n")
  cat("MODELBASED_MAX_RELATIVE_SE_DIFFERENCE=", max(abs(sqrt(diag(s$modelbased)/diag(V))-1)), "\n")
  cat("MODELBASED_MAX_ABS_CROSS_COVARIANCE=", max(abs(s$modelbased[1:6, 7:12])), "\n")
}
if (length(args) >= 2L) write.csv(comparison, args[2], row.names = FALSE)
