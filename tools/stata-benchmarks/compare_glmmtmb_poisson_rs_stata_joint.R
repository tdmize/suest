# Run here; optional output CSV. Keep raw and cluster-aligned comparisons apart.
s <- c(list(balanced = readRDS("../../tests/testthat/fixtures/glmmtmb-poisson-rs-stata-balanced-v4.rds")),
  readRDS("../../tests/testthat/fixtures/glmmtmb-poisson-rs-stata-remaining-v5.rds"))
r <- readRDS("glmmtmb_poisson_rs_r_reference.rds")
outcomes <- list(balanced = c("y1", "y2"), partial_higher = c("y1_partial", "y2_partial"),
  disjoint = c("y1_left", "y2_right"))
comparison <- do.call(rbind, lapply(names(s), function(case) {
  a <- s[[case]]; b <- r$systems[[case]]
  stopifnot(identical(names(a$coefficients), names(b$coefficients)), a$n_clusters == b$n_clusters)
  do.call(rbind, lapply(c("raw", "aligned_cluster_correction"), function(scale) {
    factor <- if (case == "disjoint" && scale != "raw") 99/98 else 1
    V <- a$vcov*factor
    data.frame(case = case, scale = scale, covariance_factor = factor,
      max_abs_coefficient = max(abs(a$coefficients-b$coefficients)),
      max_abs_covariance = max(abs(V-b$vcov)),
      max_relative_se = max(abs(sqrt(diag(V)/diag(b$vcov))-1)),
      max_abs_cross_covariance_difference = max(abs(V[1:6, 7:12]-b$vcov[1:6, 7:12])),
      abs_loglik = abs(a$logLik-sum(vapply(r$components[outcomes[[case]]], `[[`, numeric(1), "logLik"))))
  }))
}))
print(comparison, digits = 10, row.names = FALSE)
args <- commandArgs(trailingOnly = TRUE)
if (length(args)) write.csv(comparison, args[1], row.names = FALSE)
