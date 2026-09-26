glmmtmb_logit_rs_data <- function(groups = 100L, periods = 12L) {
  set.seed(2437)
  id <- rep(seq_len(groups), each = periods); x <- runif(length(id), -2, 2); z <- rnorm(length(id))
  a <- rnorm(groups, sd = .9); b <- .35*a + rnorm(groups, sd = .7)
  a2 <- .55*a + rnorm(groups, sd = .75); b2 <- .5*b + rnorm(groups, sd = .75)
  woman <- rep(rep(0:1, length.out = groups), each = periods)
  data.frame(id = factor(id), time = rep(seq_len(periods), groups), higher = ceiling(id/5),
    x = x, z = z, woman = factor(woman),
    y1 = rbinom(length(id), 1, plogis(-.3 + .55*x - .25*z + .2*woman + a[id] + b[id]*x)),
    y2 = rbinom(length(id), 1, plogis(-.15 + .4*x + .2*z - .15*woman + a2[id] + b2[id]*x)))
}

glmmtmb_logit_rs_fit <- function(formula, data, ...)
  glmmTMB::glmmTMB(formula, data = data, family = binomial("logit"), ...)

# Independent integration and analytic derivatives, including variance changes.
glmmtmb_logit_rs_moments <- function(model, data) {
  beta <- glmmTMB::fixef(model)$cond; X <- model.matrix(~x*woman+z, data)
  eta <- drop(X %*% beta); x <- data$x
  S <- as.matrix(glmmTMB::VarCorr(model)$cond$id)
  A <- S[1, 1]; B <- S[2, 2]; C <- S[1, 2]; E <- sqrt(A*B)*(1-C^2/(A*B))
  sigma <- sqrt(A+2*C*x+B*x^2)
  moments <- vapply(0:4, function(k) vapply(seq_along(x), function(i) integrate(function(u) {
    p <- plogis(eta[i]+sigma[i]*u); p1 <- p*(1-p)
    value <- switch(k+1L, p, p1, p1*(1-2*p), p1*(1-6*p+6*p^2),
      p1*(1-14*p+36*p^2-24*p^3))
    value*dnorm(u)
  }, -Inf, Inf, rel.tol = 1e-12, abs.tol = 1e-13)$value, numeric(1)), numeric(length(x)))
  n <- length(x); q <- length(beta); XX <- cbind(X, matrix(0, n, 3))
  Vg <- cbind(matrix(0, n, q), 2*(A+C*x), 2*(B*x^2+C*x), 2*E*x)
  effective <- beta["x"]+beta["x:woman1"]*(data$woman == "1"); D <- C+B*x
  direct_beta <- direct_D <- matrix(0, n, q+3)
  direct_beta[, match("x", names(beta))] <- 1
  direct_beta[, match("x:woman1", names(beta))] <- data$woman == "1"
  direct_D[, q+1:3] <- cbind(C, C+2*B*x, E)
  list(mean = moments[, 1], slope = effective*moments[, 2]+D*moments[, 3],
    G = moments[, 2]*XX+.5*moments[, 3]*Vg,
    H = effective*(moments[, 3]*XX+.5*moments[, 4]*Vg)+
      D*(moments[, 4]*XX+.5*moments[, 5]*Vg)+moments[, 2]*direct_beta+moments[, 3]*direct_D)
}

glmmtmb_logit_rs_loglik <- function(parameters, model, panels = NULL) {
  d <- model$frame; id <- as.character(d$id)
  if (is.null(panels)) panels <- unique(id)
  X <- model.matrix(model, component = "cond"); Z <- cbind(1, d$x)
  y <- model.response(d); if (is.factor(y)) y <- as.integer(y)-1L
  eta <- drop(X %*% parameters[colnames(X)])
  sd <- exp(parameters[c("log_sd_intercept", "log_sd_slope")]); rho <- tanh(parameters["atanh_rho"])
  S <- outer(sd, sd)*matrix(c(1, rho, rho, 1), 2); P <- solve(S)
  sum(vapply(panels, function(panel) {
    rows <- id == panel; zi <- Z[rows, , drop = FALSE]; yi <- y[rows]; ei <- eta[rows]
    objective <- function(u) {
      e <- ei+drop(zi %*% u)
      sum(pmax(e, 0)+log1p(exp(-abs(e)))-yi*e)+sum(u*(P %*% u))/2
    }
    gradient <- function(u) drop(crossprod(zi, plogis(ei+drop(zi %*% u))-yi)+P %*% u)
    hessian <- function(u) {
      p <- plogis(ei+drop(zi %*% u))
      crossprod(zi, zi*as.numeric(p*(1-p)))+P
    }
    u <- optim(c(0, 0), objective, gradient, method = "BFGS", control = list(reltol = 1e-14))$par
    for (j in 1:10) {
      g <- gradient(u); if (max(abs(g)) < 1e-11) break
      u <- u-drop(solve(hessian(u), g))
    }
    stopifnot(max(abs(gradient(u))) < 1e-8)
    -objective(u)-as.numeric(determinant(S, logarithm = TRUE)$modulus)/2-
      as.numeric(determinant(hessian(u), logarithm = TRUE)$modulus)/2
  }, numeric(1)))
}

test_that("logit random-slope scores match an independent bivariate Laplace likelihood", {
  skip_if_not_installed("glmmTMB", minimum_version = "1.1.14")
  d <- glmmtmb_logit_rs_data(); m <- glmmtmb_logit_rs_fit(y1 ~ x + z + (1 + x | id), d)
  a <- suest:::.suest_model_adapter(m)
  expect_identical(a$type, "glmm_logit_rs")
  c <- suest:::.suest_model_components(m, a$type, a$engine); b <- c$parameters
  expect_identical(names(b), c("(Intercept)", "x", "z", "log_sd_intercept", "log_sd_slope", "atanh_rho"))
  S <- as.matrix(glmmTMB::VarCorr(m)$cond$id)
  expect_equal(unname(b[4:5]), log(sqrt(diag(S))), ignore_attr = TRUE)
  expect_equal(unname(b[6]), atanh(cov2cor(S)[1, 2]), tolerance = 1e-10)
  expect_lt(abs(glmmtmb_logit_rs_loglik(b, m)-as.numeric(logLik(m))), 1e-7)
  J <- diag(c(rep(1, 5), 1/sqrt(1+unname(tail(m$fit$par, 1))^2)))
  expect_equal(unname(c$bread/nobs(m)), unname(J %*% vcov(m, full = TRUE) %*% J), tolerance = 1e-12)
  expect_equal(sum(rowSums(abs(c$score)) > 0), 100)
  h <- 1e-5*pmax(1, abs(b))
  for (panel in c("1", "43", "100")) {
    score <- vapply(seq_along(b), function(j) {
      plus <- minus <- b; plus[j] <- plus[j]+h[j]; minus[j] <- minus[j]-h[j]
      (glmmtmb_logit_rs_loglik(plus, m, panel)-glmmtmb_logit_rs_loglik(minus, m, panel))/(2*h[j])
    }, numeric(1))
    expect_equal(unname(colSums(c$score[d$id == panel, , drop = FALSE])), score, tolerance = 3e-5)
  }
})

test_that("logit random-slope population probabilities integrate the full random variance", {
  skip_if_not_installed("glmmTMB", minimum_version = "1.1.14")
  d <- glmmtmb_logit_rs_data(); m <- glmmtmb_logit_rs_fit(y1 ~ x + z + (1 + x | id), d)
  fit <- suest(m, m, model_names = c("Y1", "Y2"), observation_id = c("id", "time"))
  b <- coef(fit)[1:6]; names(b) <- sub("Y1::", "", names(b), fixed = TRUE)
  nd <- d[1:9, ]; nd$x <- seq(-8, 8, length.out = 9)
  eta <- drop(model.matrix(~x+z, nd) %*% b[1:3])
  S <- as.matrix(glmmTMB::VarCorr(m)$cond$id); sd <- sqrt(S[1, 1]+2*S[1, 2]*nd$x+S[2, 2]*nd$x^2)
  expected <- vapply(seq_along(eta), function(i) integrate(
    function(u) plogis(eta[i]+sd[i]*u)*dnorm(u), -Inf, Inf, rel.tol = 1e-12)$value, numeric(1))
  pred <- marginaleffects::predictions(fit, newdata = nd, vcov = FALSE)
  expect_equal(pred$estimate, rep(expected, 2), tolerance = 1e-10)
  expect_equal(suest:::.suest_predict_values(fit$models[[1]], nd, "link", "glmmTMB::glmmTMB", "glmm_logit_rs"), unname(eta))
  expect_equal(fit$comparison_scale, "predicted probabilities")
  expect_true(all(pred$estimate >= 0 & pred$estimate <= 1))
})

test_that("Logit random-slope systems aggregate group influences across sample patterns", {
  skip_if_not_installed("glmmTMB", minimum_version = "1.1.14")
  dat <- glmmtmb_logit_rs_data()
  fit_model <- function(y, d = dat)
    glmmtmb_logit_rs_fit(reformulate(c("x", "z", "(1 + x | id)"), y), d)
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
  p1 <- dat[!(dat$time == 12 & as.integer(dat$id) %% 4 == 0), ]
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
  intercept <- glmmtmb_logit_rs_fit(y1 ~ x + z + (1 | id), dat)
  expect_error(suest(first, intercept), "same type", fixed = TRUE)
})

test_that("logit random slopes propagate all parameters into marginal-effect uncertainty", {
  skip_if_not_installed("glmmTMB", minimum_version = "1.1.14")
  d <- glmmtmb_logit_rs_data()
  models <- lapply(c("y1", "y2"), function(y)
    glmmtmb_logit_rs_fit(as.formula(paste(y, "~x*woman+z+(1+x|id)")), d))
  fit <- suest(models[[1]], models[[2]], model_names = c("Y1", "Y2"), observation_id = c("id", "time"))
  nd <- d[seq(1, nrow(d), length.out = 24), ]
  G <- H <- D <- matrix(0, 2, length(coef(fit))); means <- slopes <- comparisons <- numeric(2)
  for (i in 1:2) {
    a <- glmmtmb_logit_rs_moments(models[[i]], nd); index <- (i-1)*8+1:8
    zero <- one <- nd; zero$woman <- factor(0, levels = 0:1); one$woman <- factor(1, levels = 0:1)
    a0 <- glmmtmb_logit_rs_moments(models[[i]], zero); a1 <- glmmtmb_logit_rs_moments(models[[i]], one)
    means[i] <- mean(a$mean); slopes[i] <- mean(a$slope); comparisons[i] <- mean(a1$mean-a0$mean)
    G[i, index] <- colMeans(a$G); H[i, index] <- colMeans(a$H); D[i, index] <- colMeans(a1$G-a0$G)
  }
  pred <- marginaleffects::avg_predictions(fit, newdata = nd)
  slope <- marginaleffects::avg_slopes(fit, variables = "x", newdata = nd,
    numderiv = list("fdcenter", eps = 1e-4))
  comparison <- marginaleffects::avg_comparisons(fit, variables = "woman", newdata = nd)
  expect_equal(pred$estimate, means, tolerance = 1e-10)
  expect_equal(slope$estimate, slopes, tolerance = 2e-7)
  expect_equal(comparison$estimate, comparisons, tolerance = 1e-10)
  for (pair in list(list(pred, G), list(slope, H), list(comparison, D)))
    expect_equal(pair[[1]]$std.error, sqrt(diag(pair[[2]] %*% vcov(fit) %*% t(pair[[2]]))), tolerance = 2e-5)
  # The documented centered-difference step check remains close to an
  # independent analytic gradient at both neighboring step sizes.
  for (step in c(5e-5, 2e-4)) {
    checked <- marginaleffects::avg_slopes(fit, variables = "x", newdata = nd,
      numderiv = list("fdcenter", eps = step))
    expect_equal(checked$std.error, sqrt(diag(H %*% vcov(fit) %*% t(H))),
      tolerance = 2e-5)
  }
  contrast <- suppressWarnings(marginaleffects::hypotheses(pred, "b1-b2=0"))
  delta <- G[1, ]-G[2, ]
  expect_equal(contrast$std.error, sqrt(drop(delta %*% vcov(fit) %*% delta)), tolerance = 2e-6)
  for (parameter in c("log_sd_intercept", "log_sd_slope", "atanh_rho")) {
    b <- coef(fit); b[paste0("Y1::", parameter)] <- b[paste0("Y1::", parameter)]+.2
    changed <- marginaleffects::avg_predictions(marginaleffects::set_coef(fit, b), newdata = nd, vcov = FALSE)
    expect_gt(abs(changed$estimate[1]-means[1]), 1e-5)
    expect_equal(changed$estimate[2], means[2], tolerance = 1e-10)
  }
  expect_equal(marginaleffects::avg_predictions(fit, newdata = suest_newdata(fit), vcov = FALSE)$estimate,
    marginaleffects::avg_predictions(fit, newdata = d, vcov = FALSE)$estimate, tolerance = 1e-10)
})

test_that("logit-normal integration is accurate at the quadrature boundary and in wide tails", {
  eta <- rep(c(-15, -3, 0, .2, 4, 20), 5)
  sigma <- rep(c(0, 1, 3-1e-7, 3+1e-7, 30), each = 6)
  reference <- vapply(seq_along(eta), function(i) integrate(function(u)
    plogis(eta[i]+sigma[i]*u)*dnorm(u), -Inf, Inf, rel.tol = 1e-12, abs.tol = 1e-13)$value, numeric(1))
  actual <- suest:::.suest_logit_normal_mean(eta, sigma)
  expect_lt(max(abs(actual-reference)), 1e-10)
  expect_equal(actual+suest:::.suest_logit_normal_mean(-eta, sigma), rep(1, length(eta)), tolerance = 1e-12)
  expect_equal(suest:::.suest_logit_normal_mean(c(-Inf, Inf, NA, 0), c(1, 1, 1, NA)), c(0, 1, NA, NA))
  expect_equal(suest:::.suest_logit_normal_mean(c(0, 1e5, -1e5), rep(1e6, 3)),
    pnorm(c(0, .1, -.1)), tolerance = 1e-11)
})

test_that("logit random-slope support stays restricted to regular Bernoulli fits", {
  skip_if_not_installed("glmmTMB", minimum_version = "1.1.14")
  d <- glmmtmb_logit_rs_data(); d$zero <- 0
  fit_model <- function(formula, ...) suppressWarnings(glmmtmb_logit_rs_fit(formula, d, ...))
  base <- fit_model(y1~x+z+(1+x|id))
  cases <- list(
    list(fit_model(y1~x+z+(1+x||id)), "Structured"),
    list(fit_model(y1~x+z+(1+I(x^2)|id)), "untransformed numeric"),
    list(fit_model(y1~x+z+(1+woman|id)), "untransformed numeric"),
    list(fit_model(y1~x+z+(1+x+z|id)), "one conditional random intercept"),
    list(fit_model(y1~x+z+(0+x|id)), "one conditional random intercept"),
    list(fit_model(y1~x+z+(1+x|id)+(1|higher)), "one conditional random intercept"),
    list(fit_model(y1~x+z+(1+x|id), map = list(theta = factor(c(1, 2, NA))), start = list(theta = c(0, 0, 0))), "Mapped or constrained"),
    list(fit_model(y1~x+z+(1+x|id), weights = rep(1, nrow(d))), "Weighted"),
    list(fit_model(y1~x+z+offset(zero)+(1+x|id)), "Offsets"),
    list(fit_model(cbind(y1, 2-y1)~x+z+(1+x|id)), "Grouped-binomial"),
    list(fit_model(y1~x+z+(1+x|id), ziformula = ~1), "Zero-inflation"))
  for (case in cases) expect_error(suest(case[[1]], base), case[[2]], fixed = TRUE)
  bad <- base; bad$fit$convergence <- 1L
  expect_error(suest(bad, base), "must converge")
  bad <- base; bad$fit$par[names(bad$fit$par) == "theta"] <- c(0, -30, 0)
  expect_error(suest(bad, base), "boundaries")
  random_only <- fit_model(y1~z+(1+x|id))
  fit <- suest(random_only, random_only, observation_id = c("id", "time"))
  expect_error(marginaleffects::avg_predictions(fit, newdata = d["z"]), "slope column")
  pred <- marginaleffects::avg_predictions(fit, newdata = d, vcov = FALSE)
  expect_true(all(is.finite(pred$estimate)))
  d$binary_factor <- factor(d$y1, levels = 0:1, labels = c("No", "Yes"))
  factor_model <- fit_model(binary_factor~x+z+(1+x|id))
  factor_system <- suest(factor_model, factor_model, observation_id = c("id", "time"))
  baseline <- suest(base, base, observation_id = c("id", "time"))
  expect_equal(coef(factor_system), coef(baseline), tolerance = 1e-8)
  expect_equal(vcov(factor_system), vcov(baseline), tolerance = 1e-8)
})

test_that("logit-normal integration stays smooth for mixed derivatives across its switch", {
  # At x=4 and log(sd_slope)=log(.5), eta=.2 and the integrated SD is exactly 3.
  # Equal variance perturbations in x and log SD expose a branch-level jump.
  x <- 4; theta <- log(.5); hx <- 4e-4; ht <- 1e-4
  fun <- function(x, t) suest:::.suest_logit_normal_mean(-1+.3*x, sqrt(5+exp(2*t)*x^2))
  actual <- (fun(x+hx, theta+ht)-fun(x-hx, theta+ht)-
    fun(x+hx, theta-ht)+fun(x-hx, theta-ht))/(4*hx*ht)
  moment <- function(k) integrate(function(u) {
    p <- plogis(.2+3*u); p1 <- p*(1-p)
    switch(k-1L, p1*(1-2*p), p1*(1-6*p+6*p^2),
      p1*(1-14*p+36*p^2-24*p^3))*dnorm(u)
  }, -Inf, Inf, rel.tol = 1e-12, abs.tol = 1e-13)$value
  expected <- .3*4*moment(3)+4*moment(4)+2*moment(2)
  expect_lt(abs(actual-expected), 2e-8)
})
