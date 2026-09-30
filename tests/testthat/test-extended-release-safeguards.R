test_that("beta bias-adjusted estimators are rejected explicitly", {
  skip_if_not_installed("betareg")
  set.seed(811)
  d <- data.frame(x = runif(350, -1, 1), z = runif(350, -1, 1))
  mu <- plogis(.3 + .5 * d$x)
  d$y <- rbeta(nrow(d), mu * 15, (1 - mu) * 15)
  ordinary <- betareg::betareg(y ~ x + z | z, data = d, type = "ML")
  expect_true(ordinary$converged)
  expect_s3_class(suest(ordinary, ordinary), "suest_model")
  for (kind in c("BR", "BC")) {
    adjusted <- betareg::betareg(y ~ x + z | z, data = d, type = kind)
    expect_true(adjusted$converged)
    expect_error(suest(ordinary, adjusted), "[Bb]ias|type = 'ML'")
  }
})

test_that("penalized survival fits receive an explicit estimator diagnostic", {
  skip_if_not_installed("survival")
  set.seed(813)
  d <- data.frame(x = runif(300, -1, 1), z = runif(300, -1, 1))
  d$time <- exp(.2 + .4 * d$x + rnorm(nrow(d), sd = .6))
  m <- survival::survreg(survival::Surv(time) ~ survival::ridge(x, z, theta = 2),
    data = d, dist = "lognormal", model = TRUE, x = TRUE, y = TRUE)
  expect_s3_class(m, "survreg.penal")
  expect_error(suest(m, m), "[Pp]enalized")
})

test_that("mixed interval survival scores use the log-scale likelihood derivative", {
  skip_if_not_installed("survival")
  set.seed(2491)
  d <- data.frame(x = runif(640, -1, 1), z = runif(640, -1, 1))
  latent <- .3 + .4 * d$x - .2 * d$z + rnorm(nrow(d), sd = .8)
  X <- model.matrix(~ x + z, d)
  censor <- seq_len(nrow(d)) %% 4L
  for (dist in c("gaussian", "lognormal")) for (fixed in c(FALSE, TRUE)) {
    y <- if (dist == "gaussian") latent else exp(latent)
    lo <- floor(latent / .4) * .4
    hi <- lo + .4
    if (dist == "lognormal") { lo <- exp(lo); hi <- exp(hi) }
    lo[censor == 0] <- hi[censor == 0] <- y[censor == 0]
    lo[censor == 1] <- NA_real_
    hi[censor == 2] <- NA_real_
    d$lo <- lo; d$hi <- hi
    m <- survival::survreg(survival::Surv(lo, hi, type = "interval2") ~ x + z,
      data = d, dist = dist, scale = if (fixed) .8 else 0,
      model = TRUE, x = TRUE, y = TRUE)
    par <- if (fixed) coef(m) else c(coef(m), `Log(scale)` = log(m$scale))
    ll <- function(b) {
      eta <- drop(X %*% b[1:3]); sig <- if (fixed) .8 else exp(b[4])
      lower <- if (dist == "gaussian") lo else log(lo)
      upper <- if (dist == "gaussian") hi else log(hi)
      lower[is.na(lower)] <- -Inf; upper[is.na(upper)] <- Inf
      out <- log(pnorm((upper - eta)/sig) - pnorm((lower - eta)/sig))
      exact <- censor == 0
      out[exact] <- if (dist == "gaussian") dnorm(y[exact], eta[exact], sig, log = TRUE) else
        dlnorm(y[exact], eta[exact], sig, log = TRUE)
      out
    }
    score <- vapply(seq_along(par), function(j) {
      plus <- minus <- par
      plus[j] <- plus[j] + 1e-6; minus[j] <- minus[j] - 1e-6
      (ll(plus) - ll(minus))/2e-6
    }, numeric(nrow(d)))
    components <- .suest_model_components(m, "survreg", "survival::survreg")
    expect_equal(unname(components$score), unname(score), tolerance = 2e-6)
    native_inverse_information <- m$var
    V <- native_inverse_information %*% crossprod(score) %*%
      native_inverse_information * nrow(d)/(nrow(d) - 1)
    joint <- suest(m, m)
    expect_equal(unname(vcov(joint)), rbind(cbind(V, V), cbind(V, V)), tolerance = 2e-6)
  }
})

test_that("survreg robust fitting variance is not used as sandwich bread", {
  skip_if_not_installed("survival")
  set.seed(815)
  d <- data.frame(x = runif(400, -1, 1), z = runif(400, -1, 1))
  d$time <- exp(.2 + .4 * d$x + rnorm(nrow(d), sd = .6))
  a <- survival::survreg(survival::Surv(time) ~ x, data = d, dist = "lognormal")
  b <- survival::survreg(survival::Surv(time) ~ x + z, data = d, dist = "lognormal")
  ar <- update(a, robust = TRUE)
  br <- update(b, robust = TRUE)
  ordinary <- suest(a, b, model_names = c("A", "B"))
  robust <- suest(ar, br, model_names = c("A", "B"))
  expect_equal(coef(robust), coef(ordinary), tolerance = 1e-12)
  expect_equal(vcov(robust), vcov(ordinary), tolerance = 1e-12)
})
