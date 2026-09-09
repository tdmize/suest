# Run from the package root after loading suest.
library(survey)
options(survey.lonely.psu = "fail", survey.adjust.domain.lonely = FALSE)
d <- read.csv("tools/stata-benchmarks/suest_r_survey_linear_benchmark.csv")
results <- lapply(1:6, function(case) {
  design <- if (case == 1) svydesign(~psu, strata = ~strata, weights = ~w, data = d) else
    svydesign(~psu, strata = ~strata, weights = ~w, fpc = ~pop, data = d)
  a_keep <- switch(case, rep(TRUE, nrow(d)), rep(TRUE, nrow(d)), d$unit <= 7,
    d$unit <= 5, d$psu %% 6 == 1, d$strata <= 2)
  b_keep <- switch(case, rep(TRUE, nrow(d)), rep(TRUE, nrow(d)), d$unit >= 4,
    d$unit > 5, d$psu %% 6 == 2, d$strata > 2)
  a <- svyglm(y ~ x, subset(design, a_keep))
  b <- svyglm(y2 ~ x + z, subset(design, b_keep))
  joint <- suest(a, b, survey_design = design, observation_id = "id")
  list(case = case, coefficients = coef(joint), vcov = vcov(joint),
    native_a = vcov(a), native_b = vcov(b), nobs = joint$nobs_models,
    overlap = joint$nobs_overlap, design_df = joint$survey$design_df,
    native_error = joint$survey$native_vcov_max_abs_error)
})
names(results) <- c("full_no_fpc", "full_fpc", "partial", "disjoint_shared_psu",
  "disjoint_different_psu", "disjoint_strata")
saveRDS(results, "tools/stata-benchmarks/survey_linear_r_reference.rds")
for (nm in names(results)) {
  cat("\nSURVEY CASE:", nm, "\n")
  print(results[[nm]], digits = 16)
}
