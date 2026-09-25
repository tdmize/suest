# Run after the v2 reference runner and fixture builder, from this directory.
options(digits = 12)
stata <- readRDS("../../tests/testthat/fixtures/glmmtmb-nbinom2-ri-stata.rds")
r <- readRDS("glmmtmb_nbinom2_ri_v2_r_reference.rds")
compare <- function(actual, expected, label) {
  V <- actual$vcov; W <- expected$vcov
  data.frame(case = label,
    max_abs_coefficient = max(abs(actual$coefficients - expected$coefficients)),
    max_abs_covariance = max(abs(V-W)),
    max_relative_se = max(abs(sqrt(diag(V)/diag(W))-1)),
    abs_loglik = if (!is.null(expected$logLik)) abs(actual$logLik-expected$logLik) else NA_real_)
}
results <- list()
outcomes <- names(stata$components)
for (i in seq_along(outcomes)) {
  case <- c("balanced", "partial_higher", "disjoint")[ceiling(i/2)]
  results[[outcomes[i]]] <- compare(stata$components[[i]],
    r[[case]][[paste0("model", if (i%%2) 1 else 2)]], outcomes[i])
}
for (case in names(stata$gsem)) {
  expected <- r[[case]]$system
  expected$logLik <- r[[case]]$model1$logLik + r[[case]]$model2$logLik
  results[[paste0(case, "_gsem")]] <- compare(stata$gsem[[case]], expected,
    paste0(case, "_gsem_raw"))
  corrected <- stata$gsem[[case]]
  corrected$vcov <- corrected$vcov * stata$metadata$suest_correction[case] /
    stata$metadata$gsem_correction[case]
  results[[paste0(case, "_gsem_corrected")]] <- compare(corrected, expected,
    paste0(case, "_gsem_aligned_cluster_correction"))
  results[[paste0(case, "_adaptive")]] <- compare(stata$adaptive[[case]], expected,
    paste0(case, "_adaptive_separate_approximation"))
}
results <- do.call(rbind, results)
print(results, row.names = FALSE)
args <- commandArgs(trailingOnly = TRUE)
if (length(args)) write.csv(results, args[1L], row.names = FALSE)
for (i in 1:2) {
  actual <- stata$components[[i]]$native_margins
  expected <- r$balanced[[paste0("model", i)]]$native_margins
  cat("\nNATIVE MARGINS:", outcomes[i], "\n")
  print(data.frame(estimand = names(expected$estimates),
    abs_estimate = abs(actual$estimates - expected$estimates),
    abs_se = abs(sqrt(actual$variances) - sqrt(diag(expected$vcov)))))
}
cat("\nDISJOINT MAX ABS CROSS-COVARIANCE\n")
print(c(R = max(abs(r$disjoint$system$vcov[1:5, 6:10])),
  gsem = max(abs(stata$gsem$disjoint$vcov[1:5, 6:10])),
  adaptive = max(abs(stata$adaptive$disjoint$vcov[1:5, 6:10]))))
