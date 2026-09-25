glmmtmb_poisson_rs_data <- function(groups = 100L, periods = 8L) {
  set.seed(2431)
  id <- rep(seq_len(groups), each = periods)
  x <- runif(length(id), -1.5, 1.5); z <- rnorm(length(id))
  a <- rnorm(groups, sd = .5); b <- .3*a + rnorm(groups, sd = .4)
  a2 <- .6*a + rnorm(groups, sd = .4)
  b2 <- .5*b + rnorm(groups, sd = .35)
  woman <- rep(rep(0:1, length.out = groups), each = periods)
  data.frame(id = factor(id), time = rep(seq_len(periods), groups),
    higher = ceiling(id/5), x = x, z = z, woman = factor(woman),
    y1 = rpois(length(id), exp(.6 + .3*x - .2*z + .2*woman + a[id] + b[id]*x)),
    y2 = rpois(length(id), exp(.4 + .2*x + .15*z - .1*woman + a2[id] + b2[id]*x)))
}

glmmtmb_poisson_rs_fit <- function(formula, data, ...) {
  glmmTMB::glmmTMB(formula, data = data, family = poisson(), ...)
}

# Independent bivariate Laplace likelihood, in Fisher correlation coordinates.
glmmtmb_poisson_rs_loglik <- function(parameters, model, panels = NULL) {
  X <- model.matrix(model, component = "cond")
  Z <- cbind(1, model$frame$x)
  eta <- drop(X %*% parameters[colnames(X)])
  y <- model.response(model.frame(model)); id <- as.character(model$frame$id)
  if (is.null(panels)) panels <- unique(id)
  sd <- exp(parameters[c("log_sd_intercept", "log_sd_slope")])
  rho <- tanh(parameters["atanh_rho"])
  S <- outer(sd, sd)*matrix(c(1, rho, rho, 1), 2)
  precision <- solve(S)
  sum(vapply(panels, function(panel) {
    rows <- id == panel; Zi <- Z[rows, , drop = FALSE]
    objective <- function(b) -sum(dpois(y[rows], exp(eta[rows] + drop(Zi %*% b)), log = TRUE)) +
      sum(b*(precision %*% b))/2
    gradient <- function(b) drop(crossprod(Zi, exp(eta[rows] + drop(Zi %*% b)) - y[rows]) + precision %*% b)
    mode <- optim(c(0, 0), objective, gradient, method = "BFGS",
      control = list(reltol = 1e-14, maxit = 500))$par
    # Polish the conditional mode: finite differences of the Laplace Hessian
    # need a tighter mode than the optimizer's objective-change criterion.
    for (iteration in 1:10) {
      g <- gradient(mode)
      if (max(abs(g)) < 1e-11) break
      H <- crossprod(Zi, Zi*exp(eta[rows] + drop(Zi %*% mode))) + precision
      mode <- mode - drop(solve(H, g))
    }
    stopifnot(max(abs(gradient(mode))) < 1e-8)
    H <- crossprod(Zi, Zi*exp(eta[rows] + drop(Zi %*% mode))) + precision
    -objective(mode) - as.numeric(determinant(S, logarithm = TRUE)$modulus)/2 -
      as.numeric(determinant(H, logarithm = TRUE)$modulus)/2
  }, numeric(1)))
}

test_that("Poisson random-slope scores use the full Fisher correlation transformation", {
  skip_if_not_installed("glmmTMB", minimum_version = "1.1.14")
  dat <- glmmtmb_poisson_rs_data()
  model <- glmmtmb_poisson_rs_fit(y1 ~ x + z + (1 + x | id), dat)
  adapter <- suest:::.suest_model_adapter(model)
  expect_identical(adapter$type, "glmm_poisson_rs")
  component <- suest:::.suest_model_components(model, adapter$type, adapter$engine)
  b <- component$parameters
  S <- as.matrix(glmmTMB::VarCorr(model)$cond$id)
  expect_equal(unname(b[c("log_sd_intercept", "log_sd_slope")]), log(sqrt(diag(S))), ignore_attr = TRUE)
  expect_equal(unname(b["atanh_rho"]), atanh(S[1, 2]/sqrt(prod(diag(S)))), tolerance = 1e-10)
  expect_equal(glmmtmb_poisson_rs_loglik(b, model), as.numeric(logLik(model)), tolerance = 1e-7)
  raw <- unname(tail(model$fit$par, 1))
  J <- diag(c(rep(1, length(b)-1L), 1/sqrt(1+raw^2)))
  expect_equal(unname(component$bread/nrow(dat)),
    unname(J %*% as.matrix(vcov(model, full = TRUE)) %*% J), tolerance = 1e-12)
  step <- .Machine$double.eps^(1/3)*pmax(1, abs(b))
  for (panel in c("1", "43", "100")) {
    numerical <- vapply(seq_along(b), function(j) {
      plus <- minus <- b; plus[j] <- plus[j]+step[j]; minus[j] <- minus[j]-step[j]
      (glmmtmb_poisson_rs_loglik(plus, model, panel) -
        glmmtmb_poisson_rs_loglik(minus, model, panel))/(2*step[j])
    }, numeric(1))
    expect_equal(unname(colSums(component$score[dat$id == panel, , drop = FALSE])), numerical,
      tolerance = 3e-5)
  }
})

test_that("Poisson random slopes retain covariance uncertainty in means and effects", {
  skip_if_not_installed("glmmTMB", minimum_version = "1.1.14")
  dat <- glmmtmb_poisson_rs_data()
  models <- lapply(c("y1", "y2"), function(y)
    glmmtmb_poisson_rs_fit(reformulate(c("x", "z", "woman", "(1 + x | id)"), y), dat))
  fit <- suest(models[[1]], models[[2]], model_names = c("Y1", "Y2"), observation_id = c("id", "time"))
  means <- slopes <- comparisons <- numeric(2)
  G <- H <- D <- matrix(0, 2, length(coef(fit)))
  for (i in 1:2) {
    model <- models[[i]]; beta <- glmmTMB::fixef(model)$cond
    X <- model.matrix(model, component = "cond"); x <- dat$x
    S <- as.matrix(glmmTMB::VarCorr(model)$cond$id)
    A <- S[1, 1]; B <- S[2, 2]; C <- S[1, 2]
    rho <- C/sqrt(A*B); E <- sqrt(A*B)*(1-rho^2)
    mu <- exp(drop(X %*% beta) + (A + 2*C*x + B*x^2)/2)
    log_gradient <- cbind(X, A+C*x, B*x^2+C*x, E*x)
    gradient <- mu*log_gradient
    multiplier <- beta["x"] + C + B*x
    direct <- matrix(0, nrow(dat), ncol(gradient))
    direct[, 2] <- 1; direct[, 5:7] <- cbind(C, C+2*B*x, E)
    index <- ((i-1L)*7L+1L):(i*7L)
    means[i] <- mean(mu); slopes[i] <- mean(mu*multiplier)
    G[i, index] <- colMeans(gradient)
    H[i, index] <- colMeans(gradient*multiplier + mu*direct)
    X0 <- X1 <- X; X0[, "woman1"] <- 0; X1[, "woman1"] <- 1
    mu0 <- exp(drop(X0 %*% beta) + (A + 2*C*x + B*x^2)/2)
    mu1 <- exp(drop(X1 %*% beta) + (A + 2*C*x + B*x^2)/2)
    comparisons[i] <- mean(mu1-mu0)
    D[i, index] <- colMeans(cbind(mu1*X1-mu0*X0, (mu1-mu0)*log_gradient[, 5:7]))
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
  for (parameter in c("log_sd_intercept", "log_sd_slope", "atanh_rho")) {
    b <- coef(fit); b[paste0("Y1::", parameter)] <- b[paste0("Y1::", parameter)] + .2
    changed <- marginaleffects::set_coef(fit, b)
    result <- marginaleffects::avg_predictions(changed, newdata = dat, vcov = FALSE)
    expect_gt(abs(result$estimate[1]-means[1]), 1e-4)
    expect_equal(result$estimate[2], means[2], tolerance = 1e-8)
  }
  expect_equal(suest:::.suest_predict_values(fit$models[[1]], dat, "link", "glmmTMB::glmmTMB",
    "glmm_poisson_rs"), as.numeric(model.matrix(models[[1]], component = "cond") %*% glmmTMB::fixef(models[[1]])$cond))
  expect_equal(marginaleffects::avg_predictions(fit, newdata = suest_newdata(fit), vcov = FALSE)$estimate,
    means, tolerance = 1e-8)
})

test_that("Random-slope systems aggregate group influences across sample patterns", {
  skip_if_not_installed("glmmTMB", minimum_version = "1.1.14")
  dat <- glmmtmb_poisson_rs_data()
  fit_model <- function(y, d = dat)
    glmmtmb_poisson_rs_fit(reformulate(c("x", "z", "(1 + x | id)"), y), d)
  influence <- function(model) {
    raw <- tail(model$fit$par, 1)
    J <- diag(c(rep(1, 5), 1/sqrt(1+raw^2)))
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
  p1 <- dat[!(dat$time == 8 & as.integer(dat$id) %% 4 == 0), ]
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
  expect_equal(unname(offdiag_vcov(disjoint)), matrix(0, 6, 6), tolerance = 1e-12)
  expect_equal(unname(vcov(disjoint)[1:6, 1:6]), unname(crossprod(influence(left))*50/49), tolerance = 1e-9)
  expect_equal(unname(vcov(disjoint)[7:12, 7:12]), unname(crossprod(influence(right))*50/49), tolerance = 1e-9)
  expect_error(suest(first, second, observation_id = c("id", "time"), cluster = "time"),
    "Panel IDs must be nested", fixed = TRUE)
  intercept <- glmmtmb_poisson_rs_fit(y1 ~ x + z + (1 | id), dat)
  expect_error(suest(first, intercept), "same type", fixed = TRUE)
})

test_that("Random-slope support rejects unvalidated structures and boundary fits", {
  skip_if_not_installed("glmmTMB", minimum_version = "1.1.14")
  dat <- glmmtmb_poisson_rs_data()
  fit_model <- function(formula, ...) suppressWarnings(glmmtmb_poisson_rs_fit(formula, dat, ...))
  base <- fit_model(y1 ~ x + z + (1 + x | id))
  independent <- fit_model(y1 ~ x + z + (1 + x || id))
  transformed <- fit_model(y1 ~ x + z + (1 + I(x^2) | id))
  factor_slope <- fit_model(y1 ~ x + z + (1 + woman | id))
  multiple <- fit_model(y1 ~ x + z + (1 + x + z | id))
  crossed <- fit_model(y1 ~ x + z + (1 + x | id) + (1 | higher))
  no_intercept <- fit_model(y1 ~ x + z + (0 + x | id))
  mapped <- fit_model(y1 ~ x + z + (1 + x | id), map = list(beta = factor(c(3, 2, 1))))
  weighted <- fit_model(y1 ~ x + z + (1 + x | id), weights = rep(1, nrow(dat)))
  dat$zero <- 0
  offset <- fit_model(y1 ~ x + z + offset(zero) + (1 + x | id))
  expect_error(suest(independent, base), "Structured", fixed = TRUE)
  for (m in list(transformed, factor_slope))
    expect_error(suest(m, base), "untransformed numeric", fixed = TRUE)
  for (m in list(multiple, crossed, no_intercept))
    expect_error(suest(m, base), "one conditional random intercept", fixed = TRUE)
  expect_error(suest(mapped, base), "Mapped or constrained", fixed = TRUE)
  expect_error(suest(weighted, base), "Weighted", fixed = TRUE)
  expect_error(suest(offset, base), "Offsets", fixed = TRUE)
  dat$constant_y <- 1
  boundary <- fit_model(constant_y ~ x + z + (1 + x | id))
  expect_error(suest(boundary, base), "boundaries|positive-definite Hessian")
  # The random-slope variable may be absent from the fixed equation but must
  # still be supplied for population predictions.
  random_only <- fit_model(y1 ~ z + (1 + x | id))
  fit <- suest(random_only, random_only, observation_id = c("id", "time"))
  expect_error(marginaleffects::avg_predictions(fit, newdata = dat[c("z")]), "slope column")
  S <- as.matrix(glmmTMB::VarCorr(random_only)$cond$id)
  mu <- exp(drop(model.matrix(random_only, component = "cond") %*% glmmTMB::fixef(random_only)$cond) +
    (S[1, 1] + 2*S[1, 2]*dat$x + S[2, 2]*dat$x^2)/2)
  expect_equal(marginaleffects::avg_predictions(fit, newdata = dat, vcov = FALSE)$estimate,
    rep(mean(mu), 2), tolerance = 1e-8)
})

test_that("Random-slope acceptance and population means survive a change of units", {
  skip_if_not_installed("glmmTMB", minimum_version = "1.1.14")
  dat <- glmmtmb_poisson_rs_data(); dat$scaled_x <- 1e4*dat$x
  base <- glmmtmb_poisson_rs_fit(y1 ~ x + z + (1 + x | id), dat)
  # Start at the equivalent optimum to avoid upstream exploratory overflow.
  theta <- unname(tail(base$fit$par, 3)); theta[2] <- theta[2]-log(1e4)
  scaled <- glmmtmb_poisson_rs_fit(y1 ~ scaled_x + z + (1 + scaled_x | id), dat,
    start = list(beta = unname(glmmTMB::fixef(base)$cond*c(1, 1e-4, 1)), theta = theta),
    control = glmmTMB::glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS"),
      optCtrl = list(parscale = c(1, 1e-4, 1, 1, 1, 1), reltol = 1e-12)))
  original <- suest(base, base, observation_id = c("id", "time"))
  converted <- suest(scaled, scaled, observation_id = c("id", "time"))
  a <- marginaleffects::avg_predictions(original, newdata = dat)
  b <- marginaleffects::avg_predictions(converted, newdata = dat)
  expect_equal(a$estimate, b$estimate, tolerance = 1e-5)
  expect_equal(a$std.error, b$std.error, tolerance = 2e-5)
})
