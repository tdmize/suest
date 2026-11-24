# R side of Stata round 20: predictions and marginal effects for each model
# type, to check that R's predictions are the same quantities as Stata's
# margins after suest2, and that the effects and their cross-model
# differences agree. Writes v20_r_reference.csv (R's results). Replay copy:
# the data come from v20_fixture.dta, the file Stata used.
# Usage: Rscript make_v20.R   (SUEST_DEV=<path> loads a development copy)
suppressPackageStartupMessages({
  if (nzchar(Sys.getenv("SUEST_DEV"))) pkgload::load_all(Sys.getenv("SUEST_DEV"), quiet = TRUE) else
    library(suest)
  library(marginaleffects); library(survival); library(pglm)
})
options(warn = 1)

# ---- fixture: 150 groups of 6 --------------------------------------------------
set.seed(20261124)
G <- 150L; m <- 6L; n <- G * m
id <- rep(seq_len(G), each = m); t <- rep(seq_len(m), G)
u <- rnorm(G, sd = 0.7)[id]; us <- rnorm(G, sd = 0.3)[id]
x <- rnorm(n); z <- 0.4 * x + rnorm(n); w <- rnorm(n)
f <- sample(1:3, n, TRUE, prob = c(0.4, 0.35, 0.25))
f2 <- as.numeric(f == 2); f3 <- as.numeric(f == 3)
eta <- -0.3 + 0.6 * x - 0.35 * z + 0.4 * f2 - 0.3 * f3
d <- data.frame(rid = seq_len(n), id, t, x, z, w, f2, f3)
d$y <- 2 + 0.6 * x - 0.35 * z + 0.4 * f2 - 0.3 * f3 + u + us * x + rnorm(n)
d$yb <- as.integer(eta + u + rlogis(n) > 0)
lat <- eta + 0.6 * u + rlogis(n)
d$yo <- as.integer(cut(lat, c(-Inf, -1.2, -0.1, 0.9, Inf), labels = FALSE))
d$ym <- as.integer(cut(eta + rlogis(n) + 0.8 * rnorm(n), c(-Inf, -1, 0, 1, Inf), labels = FALSE))
mu <- exp(0.4 + 0.3 * x - 0.2 * z + 0.2 * f2 - 0.15 * f3 + 0.5 * u)
d$yc <- rpois(n, mu)
d$ynb <- rnbinom(n, mu = exp(0.4 + 0.3 * x - 0.2 * z + 0.2 * f2 - 0.15 * f3), size = 1.4)
d$ypos <- rgamma(n, shape = 2, rate = 2 / exp(0.8 + 0.3 * x - 0.2 * z + 0.2 * f2))
d$yf <- plogis(eta + rnorm(n, sd = 0.7))
ystar <- 0.3 + 0.6 * x - 0.4 * z + 0.3 * f2 - 0.2 * f3 + rnorm(n)
d$ycens <- pmax(ystar, 0)
d$lo <- floor(ystar * 2) / 2; d$hi <- d$lo + 0.5
d$ytr <- ifelse(ystar > 0, ystar, NA)
d$ybeta <- rbeta(n, 5 * plogis(eta / 2), 5 * (1 - plogis(eta / 2)))
d$yh <- as.integer(0.2 + 0.6 * x - 0.3 * z + 0.3 * f2 + exp(0.35 * w) * rnorm(n) > 0)
e1 <- rnorm(n); e2 <- 0.45 * e1 + sqrt(1 - 0.45^2) * rnorm(n)
d$yb1 <- as.integer(0.2 + 0.6 * x - 0.3 * z + 0.3 * f2 + e1 > 0)
d$yb2 <- as.integer(-0.1 + 0.4 * x + 0.3 * z - 0.2 * f3 + e2 > 0)
d$instr1 <- rnorm(n); d$instr2 <- rnorm(n)
v <- rnorm(n)
d$endog <- 0.6 * d$instr1 + 0.4 * d$instr2 + 0.3 * x + v
d$yiv <- 1 + 0.5 * x - 0.3 * z + 0.6 * d$endog + 0.7 * v + rnorm(n)
d$yivb <- as.integer(0.1 + 0.5 * x - 0.3 * z + 0.4 * d$endog + 0.5 * v + rnorm(n, sd = sqrt(0.75)) > 0)
d$sel <- as.integer(0.3 + 0.5 * x + 0.8 * w + e1 > 0)
d$yheck <- ifelse(d$sel == 1, 1 + 0.6 * x - 0.4 * z + 0.3 * f2 + 1.1 * e2, NA)
d$wt <- round(runif(n, 0.5, 2.5), 3)
d$strata <- rep(1:10, length.out = G)[id]
d$psu <- id
d <- use_fixture(d, "v20_fixture.dta")

# R's view of the data: ordered and nominal outcomes as factors
r <- d
r$yo_f <- ordered(r$yo); r$ym_f <- factor(r$ym); r$id_f <- factor(r$id)
r$sel <- r$sel == 1

# ---- the systems --------------------------------------------------------------
fm <- function(y, rhs) stats::as.formula(paste(y, "~", rhs))
A <- "x + z + f2 + f3"; B <- "x + f2 + f3"
aghq <- 25L
tight <- glm.control(epsilon = 1e-12, maxit = 100)
S <- list()
add <- function(case, fit, link = FALSE, wts = FALSE) S[[case]] <<- list(fit = fit, link = link, wts = wts)

add("regress", function() list(lm(fm("y", A), r), lm(fm("y", B), r)))
add("logit", function() list(glm(fm("yb", A), binomial, r, control = tight), glm(fm("yb", B), binomial, r, control = tight)))
add("ologit", function() list(MASS::polr(fm("yo_f", A), r, Hess = TRUE, control = list(reltol = 1e-14)),
  MASS::polr(fm("yo_f", B), r, Hess = TRUE, control = list(reltol = 1e-14))))
add("oprobit", function() list(MASS::polr(fm("yo_f", A), r, method = "probit", Hess = TRUE, control = list(reltol = 1e-14)),
  MASS::polr(fm("yo_f", B), r, method = "probit", Hess = TRUE, control = list(reltol = 1e-14))))
add("mlogit", function() list(nnet::multinom(fm("ym_f", A), r, trace = FALSE, reltol = 1e-14, maxit = 500),
  nnet::multinom(fm("ym_f", B), r, trace = FALSE, reltol = 1e-14, maxit = 500)))
add("poisson", function() list(glm(fm("yc", A), poisson, r, control = tight), glm(fm("yc", B), poisson, r, control = tight)))
add("nbreg", function() list(MASS::glm.nb(fm("ynb", A), r, control = tight), MASS::glm.nb(fm("ynb", B), r, control = tight)))
add("nbreg_constant", function() list(glmmTMB::glmmTMB(fm("ynb", A), family = glmmTMB::nbinom1(), data = r),
  glmmTMB::glmmTMB(fm("ynb", B), family = glmmTMB::nbinom1(), data = r)))
add("glm_gamma", function() list(glm(fm("ypos", A), Gamma("log"), r, control = tight), glm(fm("ypos", B), Gamma("log"), r, control = tight)))
add("fracreg_logit", function() list(glm(fm("yf", A), quasibinomial("logit"), r, control = tight), glm(fm("yf", B), quasibinomial("logit"), r, control = tight)))
add("fracreg_probit", function() list(glm(fm("yf", A), quasibinomial("probit"), r, control = tight), glm(fm("yf", B), quasibinomial("probit"), r, control = tight)))
add("tobit", function() list(censReg::censReg(fm("ycens", A), left = 0, data = r), censReg::censReg(fm("ycens", B), left = 0, data = r)))
add("intreg", function() list(
  survreg(fm("Surv(lo, hi, type = 'interval2')", A), r, dist = "gaussian", control = survreg.control(rel.tolerance = 1e-13)),
  survreg(fm("Surv(lo, hi, type = 'interval2')", B), r, dist = "gaussian", control = survreg.control(rel.tolerance = 1e-13))))
add("truncreg", function() list(truncreg::truncreg(fm("ytr", A), data = r, point = 0, direction = "left"),
  truncreg::truncreg(fm("ytr", B), data = r, point = 0, direction = "left")))
add("betareg", function() list(betareg::betareg(fm("ybeta", A), r), betareg::betareg(fm("ybeta", B), r)))
add("hetprobit", function() list(Rchoice::hetprob(yh ~ x + z + f2 + f3 | w, data = r, link = "probit"),
  Rchoice::hetprob(yh ~ x + f2 + f3 | w, data = r, link = "probit")))
add("biprobit", function() list(
  mvProbit::mvProbit(cbind(yb1, yb2) ~ x + z + f2 + f3, data = r, algorithm = mvtnorm::TVPACK(), method = "BFGS",
    finalHessian = TRUE, iterlim = 300, reltol = 1e-12, printLevel = 0),
  mvProbit::mvProbit(cbind(yb1, yb2) ~ x + f2 + f3, data = r, algorithm = mvtnorm::TVPACK(), method = "BFGS",
    finalHessian = TRUE, iterlim = 300, reltol = 1e-12, printLevel = 0)))
add("ivregress", function() list(fixest::feols(yiv ~ x + z + f2 + f3 | endog ~ instr1 + instr2, data = r),
  fixest::feols(yiv ~ x + f2 + f3 | endog ~ instr1 + instr2, data = r)))
add("ivprobit", function() list(
  Rchoice::ivpml(yivb ~ x + z + f2 + f3 + endog | x + z + f2 + f3 + instr1 + instr2, data = r, messages = FALSE, printLevel = 0, reltol = 1e-12),
  Rchoice::ivpml(yivb ~ x + f2 + f3 + endog | x + f2 + f3 + instr1 + instr2, data = r, messages = FALSE, printLevel = 0, reltol = 1e-12)))
add("heckman", function() list(sampleSelection::selection(sel ~ x + w, fm("yheck", A), data = r, method = "ml"),
  sampleSelection::selection(sel ~ x + w, fm("yheck", B), data = r, method = "ml")))
add("gologit2", function() list(
  VGAM::vglm(fm("yo_f", A), VGAM::cumulative(reverse = TRUE, parallel = FALSE ~ x), data = r),
  VGAM::vglm(fm("yo_f", B), VGAM::cumulative(reverse = TRUE, parallel = TRUE), data = r)))
add("xtreg_fe", function() list(plm::plm(fm("y", A), r, index = c("id", "t"), model = "within"),
  plm::plm(fm("y", B), r, index = c("id", "t"), model = "within")))
add("xtreg_re", function() list(plm::plm(fm("y", A), r, index = c("id", "t"), model = "random"),
  plm::plm(fm("y", B), r, index = c("id", "t"), model = "random")))
add("xtreg_mle", function() list(nlme::lme(fm("y", A), random = ~ 1 | id, data = r, method = "ML", control = nlme::lmeControl(tolerance = 1e-12)),
  nlme::lme(fm("y", B), random = ~ 1 | id, data = r, method = "ML", control = nlme::lmeControl(tolerance = 1e-12))), link = TRUE)
assign("pd", r, envir = globalenv())
add("xtlogit_re", function() list(
  pglm::pglm(fm("yb", A), data = pd, index = c("id", "t"), family = binomial("logit"), model = "random", effect = "individual", R = 12, print.level = 0),
  pglm::pglm(fm("yb", B), data = pd, index = c("id", "t"), family = binomial("logit"), model = "random", effect = "individual", R = 12, print.level = 0)), link = TRUE)
add("xtprobit_re", function() list(
  pglm::pglm(fm("yb", A), data = pd, index = c("id", "t"), family = binomial("probit"), model = "random", effect = "individual", R = 12, print.level = 0),
  pglm::pglm(fm("yb", B), data = pd, index = c("id", "t"), family = binomial("probit"), model = "random", effect = "individual", R = 12, print.level = 0)), link = TRUE)
add("xtpoisson_re", function() list(
  pglm::pglm(fm("yc", A), data = pd, index = c("id", "t"), family = poisson, model = "random", effect = "individual", other = "sd", print.level = 0),
  pglm::pglm(fm("yc", B), data = pd, index = c("id", "t"), family = poisson, model = "random", effect = "individual", other = "sd", print.level = 0)), link = TRUE)
add("xtpoisson_fe", function() list(fixest::fepois(yc ~ x + z + f2 + f3 | id, data = r, notes = FALSE, glm.tol = 1e-12),
  fixest::fepois(yc ~ x + f2 + f3 | id, data = r, notes = FALSE, glm.tol = 1e-12)), link = TRUE)
add("xtgee_logit", function() list(geepack::geeglm(fm("yb", A), id = id, waves = t, data = r, family = binomial, corstr = "exchangeable"),
  geepack::geeglm(fm("yb", B), id = id, waves = t, data = r, family = binomial, corstr = "exchangeable")))
lmer_ctl <- lme4::lmerControl(optimizer = "bobyqa", optCtrl = list(rhoend = 1e-12))
add("mixed", function() list(lme4::lmer(fm("y", paste(A, "+ (1 + x | id)")), r, REML = FALSE, control = lmer_ctl),
  lme4::lmer(fm("y", paste(B, "+ (1 + x | id)")), r, REML = FALSE, control = lmer_ctl)))
glmer_ctl <- lme4::glmerControl(optimizer = "bobyqa", optCtrl = list(rhoend = 1e-12), tolPwrss = 1e-12)
add("melogit", function() list(lme4::glmer(fm("yb", paste(A, "+ (1 | id)")), r, binomial, nAGQ = aghq, control = glmer_ctl),
  lme4::glmer(fm("yb", paste(B, "+ (1 | id)")), r, binomial, nAGQ = aghq, control = glmer_ctl)))
add("meprobit", function() list(lme4::glmer(fm("yb", paste(A, "+ (1 | id)")), r, binomial("probit"), nAGQ = aghq, control = glmer_ctl),
  lme4::glmer(fm("yb", paste(B, "+ (1 | id)")), r, binomial("probit"), nAGQ = aghq, control = glmer_ctl)))
add("mepoisson", function() list(lme4::glmer(fm("yc", paste(A, "+ (1 | id)")), r, poisson, nAGQ = aghq, control = glmer_ctl),
  lme4::glmer(fm("yc", paste(B, "+ (1 | id)")), r, poisson, nAGQ = aghq, control = glmer_ctl)))
add("menbreg", function() list(glmmTMB::glmmTMB(fm("yc", paste(A, "+ (1 | id)")), family = glmmTMB::nbinom2(), data = r),
  glmmTMB::glmmTMB(fm("yc", paste(B, "+ (1 | id)")), family = glmmTMB::nbinom2(), data = r)))
add("meologit", function() list(ordinal::clmm(fm("yo_f", paste(A, "+ (1 | id_f)")), data = r, nAGQ = aghq),
  ordinal::clmm(fm("yo_f", paste(B, "+ (1 | id_f)")), data = r, nAGQ = aghq)))
design <- survey::svydesign(ids = ~psu, strata = ~strata, weights = ~wt, data = r, nest = TRUE)
add("svy_logit", function() list(survey::svyglm(fm("yb", A), design, family = quasibinomial()),
  survey::svyglm(fm("yb", B), design, family = quasibinomial())), wts = TRUE)
add("svy_ologit", function() list(survey::svyolr(fm("yo_f", A), design, control = list(reltol = 1e-14)),
  survey::svyolr(fm("yo_f", B), design, control = list(reltol = 1e-14))), wts = TRUE)
add("pweight_logit", function() list(suppressWarnings(glm(fm("yb", A), binomial, r, weights = wt, control = tight)),
  suppressWarnings(glm(fm("yb", B), binomial, r, weights = wt, control = tight))), wts = TRUE)

# ---- results ------------------------------------------------------------------
rows <- list()
keep <- function(case, stat, x) {
  x <- as.data.frame(x)
  g <- if ("group" %in% names(x)) as.character(x$group) else rep(NA_character_, nrow(x))
  model <- sub("::.*$", "", g); outcome <- ifelse(grepl("::", g), sub("^.*::", "", g), "")
  if ("hypothesis" %in% names(x)) { model <- "A-B"; outcome <- if (all(is.na(g))) "" else outcome }
  rows[[length(rows) + 1L]] <<- data.frame(case = case, stat = stat, model = model, outcome = outcome,
    estimate = x$estimate, se = x$std.error)
}
for (case in names(S)) {
  s <- S[[case]]
  t0 <- Sys.time()
  res <- tryCatch({
    ms <- s$fit()
    is_survey <- inherits(ms[[1]], "svyglm") || inherits(ms[[1]], "svyolr")
    fit <- if (is_survey) suest(ms[[1]], ms[[2]], model_names = c("A", "B"), survey_design = design, observation_id = "rid") else
      if (case == "pweight_logit") suest(ms[[1]], ms[[2]], model_names = c("A", "B"), weight_type = "pweight") else
      suest(ms[[1]], ms[[2]], model_names = c("A", "B"))
    wts <- if (s$wts) TRUE else FALSE
    keep(case, "pred", avg_predictions(fit, wts = wts))
    keep(case, "dydx", avg_slopes(fit, variables = "x", wts = wts))
    cmp <- avg_comparisons(fit, variables = "x", wts = wts)
    keep(case, "change", cmp)
    categorical <- any(fit$model_types %in% c("ologit", "oprobit", "multinom", "gologit", "re_ologit", "re_oprobit"))
    hyp <- if (categorical) function(x) {
      k <- nrow(x) / 2
      data.frame(term = as.character(x$group[seq_len(k)]), estimate = x$estimate[seq_len(k)] - x$estimate[k + seq_len(k)])
    } else difference ~ revpairwise
    dd <- hypotheses(cmp, hypothesis = hyp)
    if (categorical) dd$group <- paste0("A-B::", sub("^.*::", "", dd$term))
    keep(case, "change_diff", dd)
    if (s$link) {
      keep(case, "pred_link", avg_predictions(fit, type = "link"))
      keep(case, "dydx_link", avg_slopes(fit, variables = "x", type = "link"))
    }
    "ok"
  }, error = function(e) paste("ERROR:", conditionMessage(e)))
  cat(sprintf("%-16s %-6s %5.1fs\n", case, substr(res, 1, 60), as.numeric(difftime(Sys.time(), t0, units = "secs"))))
}
ref <- do.call(rbind, rows)
ref$model[ref$model == "A-B"] <- "A-B"
write.csv(ref, "v20_r_reference.csv", row.names = FALSE)
cat("Rows:", nrow(ref), "\n")
