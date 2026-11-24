# The same checks for every supported model: the joint covariance does not
# depend on the partner model, a model paired with itself has a zero
# difference, separate samples have zero cross-covariance, and the usual
# marginaleffects calls run. Survey models are checked in test-survey-models.R
# and two-part models in test-twopart.R.

every_model_data <- function(G = 40L, m = 6L) {
  set.seed(20261124)
  n <- G * m
  id <- rep(seq_len(G), each = m)
  time <- rep(seq_len(m), G)
  u <- rnorm(G, sd = 0.6)[id]
  us <- rnorm(G, sd = 0.8)[id]
  x <- rnorm(n); z <- 0.3 * x + rnorm(n); w <- rnorm(n)
  f <- factor(sample(c("a", "b", "c"), n, TRUE))
  fb <- (f == "b") * 0.4 - (f == "c") * 0.3
  eta <- -0.2 + 0.5 * x - 0.3 * z + fb
  d <- data.frame(id, time, x, z, w, f)
  d$y_lm <- 1 + eta + u + us * x + rnorm(n)
  d$y_bin <- as.integer(eta + u + rlogis(n) > 0)
  d$y_ord <- cut(eta + u + rlogis(n), c(-Inf, -1, 0, 1, Inf), labels = c("A", "B", "C", "D"),
    ordered_result = TRUE)
  d$y_nom <- factor(d$y_ord, ordered = FALSE)
  mu <- exp(0.3 + 0.3 * x - 0.2 * z + fb / 2 + u / 2)
  d$y_cnt <- rpois(n, mu)
  d$y_nb <- rnbinom(n, mu = mu, size = 1.5)
  d$y_zi <- ifelse(runif(n) < plogis(-1 + 0.5 * z), 0L, rpois(n, mu))
  d$y_frac <- plogis(eta + rnorm(n, sd = 0.5))
  d$y_beta <- pmin(pmax(rbeta(n, 4 * plogis(eta / 2), 4 * (1 - plogis(eta / 2))), 1e-4), 1 - 1e-4)
  d$y_pos <- rgamma(n, shape = 2, rate = 2 / exp(0.5 + 0.3 * x - 0.2 * z))
  t0 <- rweibull(n, shape = 1.3, scale = exp(1 + 0.4 * x - 0.3 * z))
  censor <- rexp(n, 0.1)
  d$time_s <- pmin(t0, censor); d$event <- as.integer(t0 <= censor)
  latent <- 0.2 + 0.6 * x - 0.4 * z + rnorm(n)
  d$y_cens <- pmax(latent, 0)
  d$y_lo <- floor(latent * 2) / 2; d$y_hi <- d$y_lo + 0.5
  d$instr1 <- rnorm(n); d$instr2 <- rnorm(n)
  v <- rnorm(n)
  d$endog <- 0.5 * d$instr1 + 0.4 * d$instr2 + 0.3 * x + v
  d$y_iv <- 1 + 0.5 * x + 0.7 * d$endog + 0.6 * v + rnorm(n)
  d$y_ivb <- as.integer(0.2 + 0.5 * x + 0.5 * d$endog + 0.5 * v + rnorm(n, sd = sqrt(0.75)) > 0)
  e1 <- rnorm(n); e2 <- 0.4 * e1 + sqrt(0.84) * rnorm(n)
  d$y_b1 <- as.integer(0.2 + 0.6 * x - 0.3 * z + e1 > 0)
  d$y_b2 <- as.integer(-0.1 + 0.4 * x + 0.3 * w + e2 > 0)
  d$y_het <- as.integer(0.2 + 0.6 * x - 0.3 * z + exp(0.3 * w) * rnorm(n) > 0)
  d$sel <- (0.3 + 0.5 * x + 0.8 * w + e1) > 0
  d$y_heck <- ifelse(d$sel, 1 + 0.6 * x - 0.4 * z + 1.2 * e2, NA)
  d$cls <- (time > 3) + 1L
  d$y_nest <- d$y_lm + rnorm(2 * G, sd = 0.8)[(id - 1L) * 2L + d$cls]
  d
}

# Each entry: packages, and a function of the data that returns three models
# (M1, M2, M3), all fit to the same data; factor = FALSE when f is not a predictor.
every_model_registry <- function() {
  three <- function(fit, rhs = c("x + z + f", "x + f", "x + w")) lapply(rhs, fit)
  list(
    lm = list(fit = function(d) three(function(r) lm(as.formula(paste("y_lm ~", r)), data = d))),
    logit = list(fit = function(d) three(function(r) glm(as.formula(paste("y_bin ~", r)), binomial, d))),
    probit = list(fit = function(d) three(function(r) glm(as.formula(paste("y_bin ~", r)), binomial("probit"), d))),
    cloglog = list(fit = function(d) three(function(r) glm(as.formula(paste("y_bin ~", r)), binomial("cloglog"), d))),
    poisson = list(fit = function(d) three(function(r) glm(as.formula(paste("y_cnt ~", r)), poisson, d))),
    negbin = list(fit = function(d) three(function(r) MASS::glm.nb(as.formula(paste("y_nb ~", r)), d))),
    gamma = list(fit = function(d) three(function(r) glm(as.formula(paste("y_pos ~", r)), Gamma("log"), d))),
    fractional = list(fit = function(d) three(function(r) glm(as.formula(paste("y_frac ~", r)), quasibinomial, d))),
    polr = list(fit = function(d) three(function(r) MASS::polr(as.formula(paste("y_ord ~", r)), d, Hess = TRUE))),
    clm = list(pkg = "ordinal", fit = function(d) three(function(r) ordinal::clm(as.formula(paste("y_ord ~", r)), data = d))),
    multinom = list(pkg = "nnet", fit = function(d) three(function(r)
      nnet::multinom(as.formula(paste("y_nom ~", r)), d, trace = FALSE))),
    nbinom1 = list(pkg = "glmmTMB", fit = function(d) three(function(r)
      glmmTMB::glmmTMB(as.formula(paste("y_nb ~", r)), family = glmmTMB::nbinom1(), data = d))),
    zeroinfl = list(pkg = "pscl", fit = function(d) three(function(r)
      pscl::zeroinfl(as.formula(paste("y_zi ~", r, "| z")), data = d))),
    censReg = list(pkg = "censReg", fit = function(d) three(function(r)
      censReg::censReg(as.formula(paste("y_cens ~", r)), data = d))),
    survreg = list(pkg = "survival", fit = function(d) three(function(r)
      survival::survreg(as.formula(paste("survival::Surv(time_s, event) ~", r)), data = d))),
    intreg = list(pkg = "survival", fit = function(d) three(function(r)
      survival::survreg(as.formula(paste("survival::Surv(y_lo, y_hi, type = 'interval2') ~", r)),
        data = d, dist = "gaussian"))),
    flexsurv = list(pkg = "flexsurv", fit = function(d) three(function(r)
      flexsurv::flexsurvreg(as.formula(paste("survival::Surv(time_s, event) ~", r)), data = d,
        dist = "gompertz"))),
    truncreg = list(pkg = "truncreg", fit = function(d) {
      dd <- d[d$y_cens > 0, ]
      three(function(r) truncreg::truncreg(as.formula(paste("y_cens ~", r)), data = dd, point = 0,
        direction = "left"))
    }),
    betareg = list(pkg = "betareg", fit = function(d) three(function(r)
      betareg::betareg(as.formula(paste("y_beta ~", r)), data = d))),
    hetprob = list(pkg = "Rchoice", fit = function(d) three(function(r)
      Rchoice::hetprob(as.formula(paste("y_het ~", r, "| w")), data = d, link = "probit"),
      c("x + z + f", "x + f", "x + z"))),
    gjrm = list(pkg = "GJRM", fit = function(d) list(
      GJRM::gjrm(list(y_b1 ~ x + z + f, y_b2 ~ x + w), data = d, margins = c("probit", "probit"), model = "B"),
      GJRM::gjrm(list(y_b1 ~ x + f, y_b2 ~ x + w), data = d, margins = c("probit", "probit"), model = "B"),
      GJRM::gjrm(list(y_b1 ~ x + f, y_b2 ~ x + z), data = d, margins = c("probit", "probit"), model = "B"))),
    feols = list(pkg = "fixest", fit = function(d) three(function(r)
      fixest::feols(as.formula(paste("y_iv ~", r, "| endog ~ instr1 + instr2")), data = d))),
    ivpml = list(pkg = "Rchoice", factor = FALSE, fit = function(d) list(
      Rchoice::ivpml(y_ivb ~ x + z + endog | x + z + instr1 + instr2, data = d, messages = FALSE, printLevel = 0),
      Rchoice::ivpml(y_ivb ~ x + endog | x + instr1 + instr2, data = d, messages = FALSE, printLevel = 0),
      Rchoice::ivpml(y_ivb ~ x + w + endog | x + w + instr1, data = d, messages = FALSE, printLevel = 0))),
    heckman = list(pkg = "sampleSelection", fit = function(d) three(function(r)
      sampleSelection::selection(sel ~ x + w, as.formula(paste("y_heck ~", r)), data = d, method = "ml"))),
    vglm = list(pkg = "VGAM", fit = function(d) list(
      VGAM::vglm(y_ord ~ x + z + f, VGAM::cumulative(reverse = TRUE, parallel = FALSE ~ x), data = d),
      VGAM::vglm(y_ord ~ x + f, VGAM::cumulative(reverse = TRUE, parallel = TRUE), data = d),
      VGAM::vglm(y_ord ~ x + w, VGAM::cumulative(reverse = TRUE, parallel = TRUE), data = d))),
    plm = list(pkg = "plm", factor = FALSE, fit = function(d) three(function(r)
      plm::plm(as.formula(paste("y_lm ~", r)), data = d, index = c("id", "time"), model = "within"),
      c("x + z", "x", "x + w"))),
    fepois = list(pkg = "fixest", fit = function(d) three(function(r)
      fixest::fepois(as.formula(paste("y_cnt ~", r, "| id")), data = d, notes = FALSE))),
    geeglm = list(pkg = "geepack", fit = function(d) three(function(r)
      geepack::geeglm(as.formula(paste("y_bin ~", r)), id = id, waves = time, data = d, family = binomial,
        corstr = "exchangeable"))),
    lme = list(pkg = "nlme", fit = function(d) three(function(r)
      nlme::lme(as.formula(paste("y_lm ~", r)), random = ~ 1 | id, data = d, method = "ML"))),
    glmmTMB = list(pkg = "glmmTMB", fit = function(d) three(function(r)
      glmmTMB::glmmTMB(as.formula(paste("y_cnt ~", r, "+ (1 | id)")), family = poisson, data = d))),
    lmer = list(pkg = "lme4", fit = function(d) three(function(r)
      lme4::lmer(as.formula(paste("y_lm ~", r, "+ (1 + x | id)")), data = d, REML = FALSE,
        control = lme4::lmerControl(optimizer = "bobyqa", optCtrl = list(rhoend = 1e-10))))),
    glmer = list(pkg = "lme4", fit = function(d) three(function(r)
      lme4::glmer(as.formula(paste("y_bin ~", r, "+ (1 | id)")), family = binomial, data = d, nAGQ = 7))),
    lmer_nested = list(pkg = "lme4", fit = function(d) three(function(r)
      lme4::lmer(as.formula(paste("y_nest ~", r, "+ (1 | id/cls)")), data = d, REML = FALSE))),
    clmm = list(pkg = "ordinal", fit = function(d) {
      d$id <- factor(d$id)
      three(function(r) ordinal::clmm(as.formula(paste("y_ord ~", r, "+ (1 | id)")), data = d))
    })
  )
}

every_model_block <- function(fit, i, j = i) fit$vcov[fit$index[[i]], fit$index[[j]], drop = FALSE]
every_model_relative <- function(a, b) max(abs(a - b)) / max(abs(b))

every_model_checks <- function(name, entry, d) {
  quiet <- function(expr) suppressWarnings(suppressMessages(expr))
  models <- quiet(entry$fit(d))
  fit <- quiet(suest(models[[1]], models[[2]], model_names = c("M1", "M2")))
  expect_s3_class(fit, "suest_model")
  expect_no_error(capture.output(print(fit), print(summary(fit))))

  # A model's block does not depend on its partner model
  other <- quiet(suest(models[[1]], models[[3]], model_names = c("M1", "M3")))
  expect_lt(every_model_relative(every_model_block(fit, 1), every_model_block(other, 1)), 1e-8)

  # A model paired with itself: cross-covariance equals its covariance
  self <- quiet(suest(models[[1]], models[[1]], model_names = c("A", "B")))
  expect_lt(every_model_relative(every_model_block(self, 1, 2), every_model_block(self, 1)), 1e-8)
  expect_lt(every_model_relative(every_model_block(self, 1), every_model_block(fit, 1)), 1e-8)

  # Separate samples: zero cross-covariance
  first <- d$id <= max(d$id) / 2
  part1 <- quiet(entry$fit(d[first, ]))
  part2 <- quiet(entry$fit(d[!first, ]))
  disjoint <- quiet(suest(part1[[1]], part2[[2]], model_names = c("A", "B")))
  expect_identical(max(abs(every_model_block(disjoint, 1, 2))), 0)

  # marginaleffects
  finite <- function(x) nrow(x) > 0 && all(is.finite(x$estimate))
  predictions <- quiet(marginaleffects::avg_predictions(fit))
  expect_true(finite(predictions))
  slopes <- quiet(marginaleffects::avg_slopes(fit, variables = "x"))
  expect_true(finite(slopes) && all(is.finite(slopes$std.error)))
  expect_true(finite(quiet(marginaleffects::avg_comparisons(fit))))
  if (!isFALSE(entry$factor)) {
    expect_true(finite(quiet(marginaleffects::avg_comparisons(fit, variables = "f"))))
    grid <- quiet(marginaleffects::predictions(fit,
      newdata = marginaleffects::datagrid(x = c(-1, 1), f = "b")))
    expect_true(finite(grid))
  }
  categorical <- any(fit$model_types %in% c("ologit", "oprobit", "multinom", "gologit", "re_ologit",
    "re_oprobit"))
  if (!categorical) {
    difference <- quiet(marginaleffects::avg_slopes(fit, variables = "x",
      hypothesis = difference ~ revpairwise))
    expect_true(finite(difference))
    zero <- quiet(marginaleffects::avg_predictions(self, hypothesis = difference ~ revpairwise))
    expect_lt(max(abs(zero$estimate)), 1e-8)
  }
}

registry <- every_model_registry()
for (name in names(registry)) local({
  entry <- registry[[name]]
  test_that(paste("every-model checks:", name), {
    skip_on_cran()
    for (pkg in entry$pkg) skip_if_not_installed(pkg)
    every_model_checks(name, entry, every_model_data())
  })
})
