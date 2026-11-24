# R side of the check of the new survey models against Stata (svy: with
# suest2). Writes survey_models_r_reference_v01.csv. Fits use tight
# convergence so that differences from Stata reflect the estimators, not
# stopping rules. Replay copy: the data come from
# survey_models_fixture_v01.dta, the file Stata used.
suppressPackageStartupMessages({
  library(suest); library(haven); library(survey); library(svyVGAM); library(survival)
})
options(survey.lonely.psu = "fail")
set.seed(20261106)
n_strata <- 6; psu_per <- 8; m <- 12
d <- expand.grid(unit = seq_len(m), psu = seq_len(psu_per), strata = seq_len(n_strata))
d$psu <- (d$strata - 1) * psu_per + d$psu
d$id <- seq_len(nrow(d))
n <- nrow(d)
d$w <- 1 + (d$id %% 7) / 3
u <- rnorm(max(d$psu), sd = 0.5)[d$psu]
d$x <- rnorm(n) + u
d$z <- rnorm(n)
d$o <- as.integer(cut(0.8 * d$x - 0.5 * d$z + 0.3 * d$x * (d$x > 0) + u + rlogis(n),
  c(-Inf, -1, 0.5, 2, Inf)))
d$m <- as.integer(cut(0.7 * d$x + 0.4 * d$z + u + rnorm(n), c(-Inf, -0.4, 0.8, Inf)))
noise <- runif(n) < 0.3
d$m[noise] <- sample(1:3, sum(noise), replace = TRUE)
d$yb <- rbinom(n, 1, plogis(-0.3 + 0.7 * d$x + 0.4 * d$z + u))
d$cnt <- rnbinom(n, size = 1.5, mu = exp(0.4 + 0.4 * d$x - 0.3 * d$z + 0.5 * u))
ystar <- 0.5 + 0.8 * d$x - 0.4 * d$z + u + rnorm(n)
d$yc <- pmax(ystar, 0)
# Interval outcome: grouped into bands, open at both ends, a few exact values.
yint <- 2 + 0.6 * d$x + 0.5 * d$z + u + rnorm(n)
breaks <- c(-Inf, 0.5, 1.5, 2.5, 3.5, Inf)
band <- findInterval(yint, breaks[-c(1, length(breaks))]) + 1
d$lo <- breaks[band]; d$hi <- breaks[band + 1]
d$lo[!is.finite(d$lo)] <- NA; d$hi[!is.finite(d$hi)] <- NA
exact <- seq(1, n, by = 9)
d$lo[exact] <- d$hi[exact] <- yint[exact]
d$t <- rweibull(n, 1.4, exp(1 + 0.3 * d$x - 0.2 * d$z))
d$ev <- rbinom(n, 1, 0.8)
d$inst <- d$x + rnorm(n)
d$endo <- 0.6 * d$inst + 0.3 * d$z + u + rnorm(n)
d$yiv <- 1 + 0.5 * d$endo - 0.4 * d$z + u + rnorm(n)
d$dom <- as.integer(d$z > -0.8)
d <- use_fixture(d, "survey_models_fixture_v01.dta")

d$of <- ordered(d$o)
d$mf <- factor(d$m)
design <- svydesign(~psu, strata = ~strata, weights = ~w, data = d)
domain <- subset(design, dom == 1)
out <- list()
add <- function(case, ...) {
  fits <- list(...)
  fit <- suest(..., survey_design = design, observation_id = "id",
    model_names = LETTERS[seq_along(fits)])
  b <- coef(fit); V <- vcov(fit)
  out[[case]] <<- rbind(
    data.frame(case = case, kind = "b", row = names(b), col = "", value = unname(b)),
    data.frame(case = case, kind = "V", row = rep(rownames(V), ncol(V)),
               col = rep(colnames(V), each = nrow(V)), value = as.vector(V)))
}
olr <- function(f, des, method = "logistic")
  svyolr(f, des, method = method, control = list(reltol = 1e-15, maxit = 5000))
vctl <- VGAM::vglm.control(criterion = "coefficients", epsilon = 1e-9, maxit = 500)
sctl <- survreg.control(rel.tolerance = 1e-13, iter.max = 200)

add("ologit_pair", olr(of ~ x + z, design), olr(of ~ x, design))
add("oprobit_pair", olr(of ~ x + z, design, "probit"), olr(of ~ x, design, "probit"))
add("mlogit_pair",
  svy_vglm(mf ~ x + z, VGAM::multinomial(refLevel = 1), design = design, control = vctl),
  svy_vglm(mf ~ x, VGAM::multinomial(refLevel = 1), design = design, control = vctl))
add("nbreg_pair",
  svy_vglm(cnt ~ x + z, VGAM::negbinomial(), design = design, control = vctl),
  svy_vglm(cnt ~ x, VGAM::negbinomial(), design = design, control = vctl))
add("gologit_pair",
  svy_vglm(of ~ x + z, VGAM::cumulative(parallel = FALSE, reverse = TRUE), design = design,
    control = vctl),
  svy_vglm(of ~ x + z, VGAM::cumulative(parallel = FALSE ~ x, reverse = TRUE), design = design,
    control = vctl))
add("tobit_pair",
  svysurvreg(Surv(yc, yc > 0, type = "left") ~ x + z, design, dist = "gaussian", control = sctl),
  svysurvreg(Surv(yc, yc > 0, type = "left") ~ x, design, dist = "gaussian", control = sctl))
add("intreg_pair",
  svysurvreg(Surv(lo, hi, type = "interval2") ~ x + z, design, dist = "gaussian", control = sctl),
  svysurvreg(Surv(lo, hi, type = "interval2") ~ x, design, dist = "gaussian", control = sctl))
add("weibull_pair",
  svysurvreg(Surv(t, ev) ~ x + z, design, dist = "weibull", control = sctl),
  svysurvreg(Surv(t, ev) ~ x, design, dist = "weibull", control = sctl))
add("ivreg_pair", svyivreg(yiv ~ endo + z | inst + z, design), svyivreg(yiv ~ endo | inst, design))
add("mixed_six",
  svyglm(yiv ~ x + z, design),
  svyglm(yb ~ x + z, design, family = quasibinomial(), control = list(epsilon = 1e-12, maxit = 100)),
  olr(of ~ x + z, design),
  svy_vglm(mf ~ x + z, VGAM::multinomial(refLevel = 1), design = design, control = vctl),
  svy_vglm(cnt ~ x + z, VGAM::negbinomial(), design = design, control = vctl),
  svysurvreg(Surv(yc, yc > 0, type = "left") ~ x + z, design, dist = "gaussian", control = sctl))
add("domain_pair", olr(of ~ x + z, domain),
  svy_vglm(mf ~ x + z, VGAM::multinomial(refLevel = 1), design = domain, control = vctl))
write.csv(do.call(rbind, out), "survey_models_r_reference_v01.csv", row.names = FALSE)
