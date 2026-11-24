# R side of rounds 5-7, the cases where R and Stata fit the same models:
# linear random-intercept pairs (lmer, Stata mixed), survey systems, gologit2,
# and Heckman. Taken from families_reference_v01.R, survey_reference_v01.R,
# gologit_reference_v01.R, and heckman_reference_v01.R (October 2026), with
# the data read from the fixtures Stata used. Writes v05_v07_r_reference.csv.
# Round 5's multiple-imputation system is fitted in compare_v05_v07.R.
suppressPackageStartupMessages({
  library(suest); library(haven); library(lme4); library(survey); library(VGAM)
  library(sampleSelection)
})
options(survey.lonely.psu = "fail")
out <- list()
add <- function(case, fit) {
  b <- coef(fit); V <- vcov(fit)
  out[[case]] <<- rbind(
    data.frame(case = case, kind = "b", row = names(b), col = "", value = unname(b)),
    data.frame(case = case, kind = "V", row = rep(rownames(V), ncol(V)),
               col = rep(colnames(V), each = nrow(V)), value = as.vector(V)))
}

# Round 6: Gaussian random-intercept pair, and with a logit
d <- as.data.frame(read_dta("families_fixture_v01.dta"))
d$id <- factor(d$id)
add("lmer", suest(lmer(yc1 ~ x + z + (1 | id), data = d, REML = FALSE),
                  lmer(yc2 ~ x + (1 | id), data = d, REML = FALSE), model_names = c("A", "B")))
add("lmer_logit", suest(lmer(yc1 ~ x + z + (1 | id), data = d, REML = FALSE),
                        glm(yb2 ~ x, binomial, d), model_names = c("A", "B")))

# Round 7: survey systems (svy: with suest2)
d <- as.data.frame(read_dta("survey_fixture_v01.dta"))
design <- svydesign(~psu, strata = ~strata, weights = ~w, data = d)
ctl <- list(epsilon = 1e-12, maxit = 100)
lin <- svyglm(y ~ x + z, design)
logit <- svyglm(yb ~ x + z, design, family = quasibinomial("logit"), control = ctl)
pois <- svyglm(count ~ x, design, family = quasipoisson(), control = ctl)
add("regress_logit_poisson", suest(lin, logit, pois, survey_design = design,
    observation_id = "id", model_names = c("A", "B", "C")))
logit_sub <- svyglm(yb ~ x, subset(design, sub == 1), family = quasibinomial("logit"), control = ctl)
add("logit_sub_poisson", suest(logit_sub, pois, survey_design = design,
    observation_id = "id", model_names = c("A", "B")))

# Round 7: generalized ordered logit (gologit2)
d <- as.data.frame(read_dta("gologit_fixture_v01.dta"))
d$yo <- ordered(d$y)
vctl <- vglm.control(epsilon = 1e-12, maxit = 100)
full <- vglm(yo ~ x + z, cumulative(parallel = FALSE, reverse = TRUE), data = d, control = vctl)
partial <- vglm(yo ~ x + z, cumulative(parallel = FALSE ~ x, reverse = TRUE), data = d, control = vctl)
add("gologit_full_partial", suest(full, partial, model_names = c("A", "B")))
add("gologit_logit", suest(partial, glm(yb ~ x + z, binomial, d), model_names = c("A", "B")))

# Round 7: Heckman selection (maximum likelihood) and a linear model
d <- as.data.frame(read_dta("heckman_fixture_v01.dta"))
h <- selection(s ~ x + w, y ~ x + z, data = d, method = "ml",
               control = list(tol = 1e-12, reltol = 1e-14))
add("heckman_regress", suest(h, lm(yc ~ x + z, d), model_names = c("A", "B")))

write.csv(do.call(rbind, out), "v05_v07_r_reference.csv", row.names = FALSE)
