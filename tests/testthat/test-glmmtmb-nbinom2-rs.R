glmmtmb_nbinom2_rs_data <- function(groups = 100L, periods = 10L) {
  set.seed(2432)
  id <- rep(seq_len(groups), each = periods)
  x <- runif(length(id), -1.5, 1.5); z <- rnorm(length(id))
  a <- rnorm(groups, sd = .6); b <- .35*a + rnorm(groups, sd = .45)
  a2 <- .6*a + rnorm(groups, sd = .45); b2 <- .6*b + rnorm(groups, sd = .4)
  woman <- rep(rep(0:1, length.out = groups), each = periods)
  data.frame(id = factor(id), time = rep(seq_len(periods), groups), higher = ceiling(id/5),
    x = x, z = z, woman = factor(woman),
    y1 = rnbinom(length(id), size = 2.5, mu = exp(.8 + .3*x - .2*z + .2*woman + a[id] + b[id]*x)),
    y2 = rnbinom(length(id), size = 3.5, mu = exp(.7 + .2*x + .15*z - .1*woman + a2[id] + b2[id]*x)))
}

glmmtmb_nbinom2_rs_fit <- function(formula, data, ...) {
  glmmTMB::glmmTMB(formula, data = data, family = glmmTMB::nbinom2("log"), ...)
}

# Independent bivariate Laplace likelihood with estimated NB2 dispersion.
glmmtmb_nbinom2_rs_loglik <- function(parameters, model, panels = NULL) {
  X <- model.matrix(model, component = "cond"); Z <- cbind(1, model$frame$x)
  eta <- drop(X %*% parameters[colnames(X)])
  y <- model.response(model.frame(model)); id <- as.character(model$frame$id)
  if (is.null(panels)) panels <- unique(id)
  phi <- exp(parameters["log_phi"])
  sd <- exp(parameters[c("log_sd_intercept", "log_sd_slope")]); rho <- tanh(parameters["atanh_rho"])
  S <- outer(sd, sd)*matrix(c(1, rho, rho, 1), 2); precision <- solve(S)
  sum(vapply(panels, function(panel) {
    rows <- id == panel; Zi <- Z[rows, , drop = FALSE]; yi <- y[rows]
    objective <- function(b) -sum(dnbinom(yi, size = phi, mu = exp(eta[rows] + drop(Zi %*% b)), log = TRUE)) +
      sum(b*(precision %*% b))/2
    gradient <- function(b) {
      mu <- exp(eta[rows] + drop(Zi %*% b))
      drop(crossprod(Zi, phi*(mu-yi)/(phi+mu)) + precision %*% b)
    }
    hessian <- function(b) {
      mu <- exp(eta[rows] + drop(Zi %*% b))
      crossprod(Zi, Zi*as.numeric(phi*mu*(phi+yi)/(phi+mu)^2)) + precision
    }
    mode <- optim(c(0, 0), objective, gradient, method = "BFGS",
      control = list(reltol = 1e-14, maxit = 500))$par
    for (iteration in 1:10) {
      g <- gradient(mode)
      if (max(abs(g)) < 1e-11) break
      mode <- mode-drop(solve(hessian(mode), g))
    }
    stopifnot(max(abs(gradient(mode))) < 1e-8)
    -objective(mode)-as.numeric(determinant(S, logarithm = TRUE)$modulus)/2-
      as.numeric(determinant(hessian(mode), logarithm = TRUE)$modulus)/2
  }, numeric(1)))
}

test_that("NB2 random-slope scores retain dispersion and covariance derivatives", {
  skip_if_not_installed("glmmTMB", minimum_version = "1.1.14")
  dat <- glmmtmb_nbinom2_rs_data()
  model <- glmmtmb_nbinom2_rs_fit(y1 ~ x + z + (1 + x | id), dat)
  adapter <- suest:::.suest_model_adapter(model)
  expect_identical(adapter$type, "glmm_nbinom2_rs")
  component <- suest:::.suest_model_components(model, adapter$type, adapter$engine)
  b <- component$parameters
  expect_identical(names(b), c("(Intercept)", "x", "z", "log_phi", "log_sd_intercept", "log_sd_slope", "atanh_rho"))
  expect_equal(unname(b["log_phi"]), log(sigma(model)), tolerance = 1e-12)
  S <- as.matrix(glmmTMB::VarCorr(model)$cond$id)
  expect_equal(unname(b[c("log_sd_intercept", "log_sd_slope")]), log(sqrt(diag(S))), ignore_attr = TRUE)
  expect_equal(unname(b["atanh_rho"]), atanh(cov2cor(S)[1, 2]), tolerance = 1e-10)
  expect_lt(abs(glmmtmb_nbinom2_rs_loglik(b, model)-as.numeric(logLik(model))), 1e-5)
  raw <- unname(tail(model$fit$par, 1)); J <- diag(c(rep(1, 6), 1/sqrt(1+raw^2)))
  expect_equal(unname(component$bread/nrow(dat)),
    unname(J %*% as.matrix(vcov(model, full = TRUE)) %*% J), tolerance = 1e-12)
  expect_equal(sum(rowSums(abs(component$score)) > 0), nlevels(dat$id))
  step <- .Machine$double.eps^(1/3)*pmax(1, abs(b))
  for (panel in c("1", "43", "100")) {
    numerical <- vapply(seq_along(b), function(j) {
      plus <- minus <- b; plus[j] <- plus[j]+step[j]; minus[j] <- minus[j]-step[j]
      (glmmtmb_nbinom2_rs_loglik(plus, model, panel)-glmmtmb_nbinom2_rs_loglik(minus, model, panel))/(2*step[j])
    }, numeric(1))
    expect_equal(unname(colSums(component$score[dat$id == panel, , drop = FALSE])), numerical, tolerance = 3e-5)
  }
})

test_that("NB2 random-slope means and effects include covariance but not a direct dispersion effect", {
  skip_if_not_installed("glmmTMB", minimum_version = "1.1.14")
  dat <- glmmtmb_nbinom2_rs_data()
  models <- lapply(c("y1", "y2"), function(y)
    glmmtmb_nbinom2_rs_fit(reformulate(c("x*woman", "z", "(1 + x | id)"), y), dat))
  fit <- suest(models[[1]], models[[2]], model_names = c("Y1", "Y2"), observation_id = c("id", "time"))
  means <- slopes <- comparisons <- numeric(2)
  G <- H <- D <- matrix(0, 2, length(coef(fit)))
  for (i in 1:2) {
    model <- models[[i]]; beta <- glmmTMB::fixef(model)$cond; k <- length(beta); p <- k+4L
    X <- model.matrix(model, component = "cond"); x <- dat$x; woman <- as.numeric(dat$woman == "1")
    S <- as.matrix(glmmTMB::VarCorr(model)$cond$id)
    A <- S[1, 1]; B <- S[2, 2]; C <- S[1, 2]; rho <- C/sqrt(A*B); E <- sqrt(A*B)*(1-rho^2)
    mu <- exp(drop(X %*% beta)+(A+2*C*x+B*x^2)/2)
    nuisance <- cbind(0, A+C*x, B*x^2+C*x, E*x)
    gradient <- mu*cbind(X, nuisance)
    multiplier <- beta["x"]+beta["x:woman1"]*woman+C+B*x
    direct <- matrix(0, nrow(dat), p)
    direct[, match(c("x", "x:woman1"), names(beta))] <- cbind(1, woman)
    direct[, k+2:4] <- cbind(C, C+2*B*x, E)
    index <- (i-1L)*p+seq_len(p)
    means[i] <- mean(mu); slopes[i] <- mean(mu*multiplier)
    G[i, index] <- colMeans(gradient); H[i, index] <- colMeans(gradient*multiplier+mu*direct)
    X0 <- X1 <- X; X0[, c("woman1", "x:woman1")] <- 0
    X1[, "woman1"] <- 1; X1[, "x:woman1"] <- x
    mu0 <- exp(drop(X0 %*% beta)+(A+2*C*x+B*x^2)/2)
    mu1 <- exp(drop(X1 %*% beta)+(A+2*C*x+B*x^2)/2)
    comparisons[i] <- mean(mu1-mu0)
    D[i, index] <- colMeans(cbind(mu1*X1-mu0*X0, (mu1-mu0)*nuisance))
  }
  pred <- marginaleffects::avg_predictions(fit, newdata = dat)
  slope <- marginaleffects::avg_slopes(fit, variables = "x", newdata = dat)
  comparison <- marginaleffects::avg_comparisons(fit, variables = "woman", newdata = dat)
  expect_equal(pred$estimate, means, tolerance = 1e-8)
  expect_equal(pred$std.error, sqrt(diag(G %*% vcov(fit) %*% t(G))), tolerance = 2e-6)
  expect_equal(slope$estimate, slopes, tolerance = 2e-6)
  expect_equal(slope$std.error, sqrt(diag(H %*% vcov(fit) %*% t(H))), tolerance = 2e-5)
  expect_equal(comparison$estimate, comparisons, tolerance = 1e-8)
  expect_equal(comparison$std.error, sqrt(diag(D %*% vcov(fit) %*% t(D))), tolerance = 2e-6)
  contrast <- suppressWarnings(marginaleffects::hypotheses(pred, "b1 - b2 = 0"))
  expect_equal(contrast$std.error, sqrt(drop((G[1, ]-G[2, ]) %*% vcov(fit) %*% (G[1, ]-G[2, ]))), tolerance = 2e-6)
  for (parameter in c("log_phi", "log_sd_intercept", "log_sd_slope", "atanh_rho")) {
    b <- coef(fit); b[paste0("Y1::", parameter)] <- b[paste0("Y1::", parameter)]+.2
    changed <- marginaleffects::set_coef(fit, b)
    result <- marginaleffects::avg_predictions(changed, newdata = dat, vcov = FALSE)
    if (parameter == "log_phi") expect_equal(result$estimate[1], means[1], tolerance = 1e-8)
    else expect_gt(abs(result$estimate[1]-means[1]), 1e-4)
    expect_equal(result$estimate[2], means[2], tolerance = 1e-8)
  }
  expect_equal(suest:::.suest_predict_values(fit$models[[1]], dat, "link", "glmmTMB::glmmTMB",
    "glmm_nbinom2_rs"), as.numeric(model.matrix(models[[1]], component = "cond") %*% glmmTMB::fixef(models[[1]])$cond))
  expect_equal(marginaleffects::avg_predictions(fit, newdata = suest_newdata(fit), vcov = FALSE)$estimate,
    means, tolerance = 1e-8)
})

test_that("NB2 random-slope systems align all seven parameter influences across samples", {
  skip_if_not_installed("glmmTMB", minimum_version = "1.1.14")
  dat <- glmmtmb_nbinom2_rs_data()
  fit_model <- function(y, d = dat) glmmtmb_nbinom2_rs_fit(reformulate(c("x", "z", "(1 + x | id)"), y), d)
  influence <- function(model) {
    raw <- unname(tail(model$fit$par, 1)); J <- diag(c(rep(1, 6), 1/sqrt(1+raw^2)))
    sandwich::estfun(model, full = TRUE) %*% as.matrix(vcov(model, full = TRUE)) %*% J
  }
  first <- fit_model("y1"); second <- fit_model("y2")
  U <- cbind(influence(first), influence(second))
  system <- suest(first, second, observation_id = c("id", "time"))
  expect_equal(unname(vcov(system)), unname(crossprod(U)*100/99), tolerance = 1e-9)
  higher <- suest(first, second, observation_id = c("id", "time"), cluster = "higher")
  grouped <- rowsum(U, ceiling(as.integer(rownames(U))/5))
  expect_equal(unname(vcov(higher)), unname(crossprod(grouped)*20/19), tolerance = 1e-9)
  expect_identical(higher$n_clusters, 20L)
  p1 <- dat[!(dat$time == 10 & as.integer(dat$id) %% 4 == 0), ]
  p2 <- dat[!(dat$time == 1 & as.integer(dat$id) %% 5 == 0), ]
  m1 <- fit_model("y1", p1); m2 <- fit_model("y2", p2)
  partial <- suest(m1, m2, observation_id = c("id", "time"), cluster = "higher")
  grouped <- rowsum(cbind(influence(m1), influence(m2)), ceiling(as.integer(rownames(U))/5))
  expect_equal(unname(vcov(partial)), unname(crossprod(grouped)*20/19), tolerance = 1e-9)
  expect_equal(partial$nobs_overlap, nrow(dat)-45L)
  expect_equal(partial$nobs_union, nrow(dat))
  left <- fit_model("y1", dat[as.integer(dat$id) %% 2 == 1, ])
  right <- fit_model("y2", dat[as.integer(dat$id) %% 2 == 0, ])
  disjoint <- suest(left, right, observation_id = c("id", "time"))
  expect_equal(unname(offdiag_vcov(disjoint)), matrix(0, 7, 7), tolerance = 1e-12)
  expect_equal(unname(vcov(disjoint)[1:7, 1:7]), unname(crossprod(influence(left))*50/49), tolerance = 1e-9)
  expect_equal(unname(vcov(disjoint)[8:14, 8:14]), unname(crossprod(influence(right))*50/49), tolerance = 1e-9)
  expect_error(suest(first, second, observation_id = c("id", "time"), cluster = "time"), "Panel IDs must be nested", fixed = TRUE)
  intercept <- glmmtmb_nbinom2_rs_fit(y1 ~ x + z + (1 | id), dat)
  poisson <- glmmTMB::glmmTMB(y1 ~ x + z + (1 + x | id), data = dat, family = poisson())
  expect_error(suest(first, intercept), "same type", fixed = TRUE)
  expect_error(suest(first, poisson), "same type", fixed = TRUE)
})

test_that("NB2 random slopes enforce the dispersion and random-effect support contract", {
  skip_if_not_installed("glmmTMB", minimum_version = "1.1.14")
  dat <- glmmtmb_nbinom2_rs_data(); dat$zero <- 0; dat$log_phi <- dat$x
  fit_model <- function(formula = y1 ~ x + z + (1 + x | id), ...)
    suppressWarnings(glmmtmb_nbinom2_rs_fit(formula, dat, ...))
  base <- fit_model()
  for (m in list(fit_model(dispformula = ~x), fit_model(dispformula = ~1 + offset(x)),
      fit_model(start = list(betadisp = log(2.5)), map = list(betadisp = factor(NA)))))
    expect_error(suest(m, base), "dispersion", ignore.case = TRUE)
  for (m in list(fit_model(map = list(beta = factor(c(3, 2, 1)))),
      fit_model(start = list(theta = c(-.5, -.8, .3)), map = list(theta = factor(c(1, 2, NA))))))
    expect_error(suest(m, base), "Mapped or constrained", fixed = TRUE)
  expect_error(suest(fit_model(y1 ~ log_phi + z + (1 + x | id)), base), "reserves coefficient names", fixed = TRUE)
  dat$atanh_rho <- dat$x
  expect_error(suest(fit_model(y1 ~ atanh_rho + z + (1 + x | id)), base), "reserve names", fixed = TRUE)
  expect_error(suest(fit_model(ziformula = ~1), base), "Zero-inflation", fixed = TRUE)
  expect_error(suest(fit_model(weights = rep(1, nrow(dat))), base), "Weighted", fixed = TRUE)
  expect_error(suest(fit_model(y1 ~ x + z + offset(zero) + (1 + x | id)), base), "Offsets", fixed = TRUE)
  expect_error(suest(fit_model(y1 ~ x + z + (1 + x || id)), base), "Structured", fixed = TRUE)
  for (formula in list(y1 ~ x + z + (1 + I(x^2) | id), y1 ~ x + z + (1 + woman | id)))
    expect_error(suest(fit_model(formula), base), "untransformed numeric", fixed = TRUE)
  for (formula in list(y1 ~ x + z + (1 + x + z | id), y1 ~ x + z + (0 + x | id),
      y1 ~ x + z + (1 + x | id) + (1 | higher)))
    expect_error(suest(fit_model(formula), base), "one conditional random intercept", fixed = TRUE)
  dat$constant_y <- 1
  expect_error(suest(fit_model(constant_y ~ x + z + (1 + x | id)), base),
    "boundaries|positive-definite Hessian|zero boundary")
  random_only <- fit_model(y1 ~ z + (1 + x | id))
  fit <- suest(random_only, random_only, observation_id = c("id", "time"))
  expect_error(marginaleffects::avg_predictions(fit, newdata = dat["z"]), "slope column")
  S <- as.matrix(glmmTMB::VarCorr(random_only)$cond$id)
  mu <- exp(drop(model.matrix(random_only, component = "cond") %*% glmmTMB::fixef(random_only)$cond)+
    (S[1, 1]+2*S[1, 2]*dat$x+S[2, 2]*dat$x^2)/2)
  expect_equal(marginaleffects::avg_predictions(fit, newdata = dat, vcov = FALSE)$estimate,
    rep(mean(mu), 2), tolerance = 1e-8)
})

test_that("NB2 random-slope units preserve inference with scale-aware native curvature", {
  skip_if_not_installed("glmmTMB", minimum_version = "1.1.14")
  dat <- glmmtmb_nbinom2_rs_data(); dat$scaled_x <- 1e4*dat$x
  base <- glmmtmb_nbinom2_rs_fit(y1 ~ x + z + (1 + x | id), dat,
    control = glmmTMB::glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS"),
      optCtrl = list(reltol = 1e-12, maxit = 500)))
  theta <- unname(tail(base$fit$par, 3)); theta[2] <- theta[2]-log(1e4)
  scaled <- glmmtmb_nbinom2_rs_fit(y1 ~ scaled_x + z + (1 + scaled_x | id), dat,
    start = list(beta = unname(glmmTMB::fixef(base)$cond*c(1, 1e-4, 1)),
      betadisp = unname(glmmTMB::fixef(base)$disp), theta = theta),
    control = glmmTMB::glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS"),
      optCtrl = list(parscale = c(1, 1e-4, 1, 1, 1, 1, 1), reltol = 1e-12, maxit = 500)))
  # TMB's default numerical Hessian uses an absolute .001 step, larger than
  # beta_scaled. Supply checked native curvature; suest must preserve it.
  coordinate_scale <- c(1, 1e-4, 1, 1, 1, 1, 1)
  curvature <- function(step) stats::optimHess(scaled$fit$par, scaled$obj$fn, scaled$obj$gr,
    control = list(ndeps = step*coordinate_scale))
  H <- curvature(.001); Hhalf <- curvature(.0005)
  J <- diag(1/coordinate_scale)
  expect_equal(J %*% solve(H) %*% J, J %*% solve(Hhalf) %*% J, tolerance = 1e-6)
  scaled$sdr <- TMB::sdreport(scaled$obj, par.fixed = scaled$fit$par, hessian.fixed = H)
  expect_equal(unname(vcov(base, full = TRUE)),
    unname(J %*% vcov(scaled, full = TRUE) %*% J), tolerance = 1e-6)
  first <- suest(base, base, observation_id = c("id", "time"))
  second <- suest(scaled, scaled, observation_id = c("id", "time"))
  a <- marginaleffects::avg_predictions(first, newdata = dat)
  b <- marginaleffects::avg_predictions(second, newdata = dat)
  expect_equal(a$estimate, b$estimate, tolerance = 1e-5)
  expect_equal(a$std.error, b$std.error, tolerance = 2e-5)
})
