# R side of round 10 (Batch 2), the cases where R and Stata fit the same
# models: linear random-slope and three-level models (lmer; Stata mixed),
# with suest2, and the survey complementary log-log pair. From
# rs_reference_v01.R (November 7, 2026), reading the fixtures Stata used.
# Writes rs_r_reference_v01.csv.
suppressPackageStartupMessages({
  library(suest); library(haven); library(lme4); library(survey)
})
d <- as.data.frame(read_dta("rs_fixture_v01.dta"))
d$school <- factor(d$school); d$class <- factor(d$class)

out <- list()
add <- function(case, ..., survey_design = NULL) {
  fit <- suest(..., model_names = LETTERS[seq_len(...length())], survey_design = survey_design,
    observation_id = "id")
  b <- coef(fit); V <- vcov(fit)
  out[[case]] <<- rbind(
    data.frame(case = case, kind = "b", row = names(b), col = "", value = unname(b)),
    data.frame(case = case, kind = "V", row = rep(rownames(V), ncol(V)),
               col = rep(colnames(V), each = nrow(V)), value = as.vector(V)))
}
nlopt <- list(xtol_abs = 1e-12, ftol_abs = 1e-14, xtol_rel = 1e-12, ftol_rel = 1e-14, maxeval = 1e5)
tight_lmer <- lmerControl(optimizer = "nloptwrap", optCtrl = nlopt)

lmer_slopes <- lmer(yc ~ x1 + x2 + (1 + x1 + x2 | school), d, REML = FALSE, control = tight_lmer)
lmer_nested <- lmer(yc2 ~ x1 + x2 + (1 + x1 | school) + (1 | school:class), d, REML = FALSE,
  control = tight_lmer)
add("mixed_pair", lmer_slopes, lmer_nested)
add("mixed_logit", lmer_nested, glm(yb ~ x1 + x2, binomial, d))

s <- as.data.frame(read_dta("survey_models_fixture_v01.dta"))
design <- svydesign(~psu, strata = ~strata, weights = ~w, data = s)
glm_ctl <- list(epsilon = 1e-12, maxit = 100)
add("svy_cloglog_pair",
  svyglm(yb ~ x + z, design, family = quasibinomial("cloglog"), control = glm_ctl),
  svyglm(yb ~ x, design, family = quasibinomial("cloglog"), control = glm_ctl),
  survey_design = design)

write.csv(do.call(rbind, out), "rs_r_reference_v01.csv", row.names = FALSE)
