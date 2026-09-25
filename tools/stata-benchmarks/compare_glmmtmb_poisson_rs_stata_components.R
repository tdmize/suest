# Run here after building the component fixture and the existing R reference.
stata <- readRDS("../../tests/testthat/fixtures/glmmtmb-poisson-rs-stata-components.rds")
reference <- readRDS("glmmtmb_poisson_rs_r_reference.rds")
comparisons <- do.call(rbind, lapply(names(stata$components), function(y) {
  s <- stata$components[[y]]; r <- reference$components[[y]]
  stopifnot(identical(names(s$coefficients), names(r$coefficients)), s$nobs == r$nobs)
  data.frame(outcome = y, max_abs_coefficient = max(abs(s$coefficients-r$coefficients)),
    max_abs_vcov = max(abs(s$vcov-r$vcov)),
    max_relative_se = max(abs(sqrt(diag(s$vcov)/diag(r$vcov))-1)),
    abs_logLik = abs(s$logLik-r$logLik))
}))
print(comparisons, digits = 10, row.names = FALSE)
for (y in c("y1", "y2")) {
  s <- stata$components[[y]]$native_margins; r <- reference$components[[y]]$native_margins
  for (i in seq_along(s)) {
    cat(y, names(s)[i], "estimate difference=", s[[i]]["estimate"]-r$estimates[i],
      "relative SE difference=", sqrt(s[[i]]["variance"]/r$vcov[i, i])-1, "\n")
  }
}
write.csv(comparisons, "glmmtmb_poisson_rs_stata_component_comparison.csv", row.names = FALSE)
cat("BALANCED_STATA_COMPONENT_LOGLIK_SUM=", format(sum(vapply(stata$components[c("y1", "y2")],
  `[[`, numeric(1), "logLik")), digits = 16), "\n", sep = "")
