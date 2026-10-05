usability_data <- function() {
  set.seed(20261003)
  d <- data.frame(id = 1:400, x = rnorm(400), z = rnorm(400),
                  g = sample(1:3, 400, TRUE))
  d$y <- rbinom(400, 1, plogis(-.2 + .8 * d$x + .3 * d$g + .4 * d$z))
  d$pos <- exp(d$z)
  d
}

test_that("unsupported S4 models get the unsupported-model error", {
  methods::setClass("suestFakeS4", representation(x = "numeric"))
  d <- usability_data()
  m <- glm(y ~ x, binomial, d)
  expect_error(
    suest(methods::new("suestFakeS4", x = 1), m),
    "Unsupported model .* has class 'suestFakeS4'"
  )
})

test_that("na.exclude fits give the same system as na.omit fits", {
  d <- usability_data()
  d$z[c(3, 50, 77)] <- NA
  omit <- suest(glm(y ~ x + z, binomial, d), glm(y ~ x, binomial, d))
  excl <- suest(glm(y ~ x + z, binomial, d, na.action = na.exclude),
                glm(y ~ x, binomial, d, na.action = na.exclude))
  expect_equal(unname(vcov(excl)), unname(vcov(omit)), tolerance = 1e-12)
  expect_equal(unname(nobs(excl)), unname(nobs(omit)))
})

test_that("a data object changed between fits is not silently aligned", {
  d <- usability_data()
  m1 <- glm(y ~ x + z, binomial, d)
  d <- d[d$id %% 3 != 0, ]
  rownames(d) <- NULL
  m2 <- glm(y ~ x + z + g, binomial, d)
  expect_error(suest(m1, m2), "same data object")
})

test_that("overlapping subsets from different data objects are flagged", {
  d <- usability_data()
  dA <- d[d$id <= 300, ]
  dB <- d[d$id > 100, ]
  mA <- glm(y ~ x + z, binomial, dA)
  mB <- glm(y ~ x + z, binomial, dB)
  expect_warning(fit <- suest(mA, mB), "different data objects")
  expect_equal(fit$nobs_overlap, 0L)
  expect_no_warning(fit_id <- suest(mA, mB, observation_id = "id"))
  expect_equal(fit_id$nobs_overlap, 200L)

  dC <- d[d$id <= 200, ]
  dD <- d[d$id > 200, ]
  expect_no_warning(suest(glm(y ~ x, binomial, dC), glm(y ~ x, binomial, dD)))
})

test_that("suest_newdata keeps variables used inside formula transformations", {
  d <- usability_data()
  m1 <- glm(y ~ x + factor(g), binomial, d)
  m2 <- glm(y ~ x + factor(g) + log(pos), binomial, d)
  fit <- suest(m1, m2)
  nd <- suest_newdata(fit)
  expect_true(all(c("g", "pos") %in% names(nd)))
  own <- marginaleffects::avg_comparisons(fit, variables = "x", newdata = nd)
  common <- marginaleffects::avg_comparisons(fit, variables = "x", newdata = d)
  expect_equal(own$estimate, common$estimate, tolerance = 1e-10)
  expect_equal(own$std.error, common$std.error, tolerance = 1e-8)

  a <- glm(y ~ x + log(pos), binomial, d, subset = g == 1)
  b <- glm(y ~ x + log(pos), binomial, d, subset = g != 1)
  sub <- suest(a, b)
  joint <- marginaleffects::avg_slopes(sub, variables = "x", newdata = suest_newdata(sub))
  separate <- c(marginaleffects::avg_slopes(a, variables = "x")$estimate,
                marginaleffects::avg_slopes(b, variables = "x")$estimate)
  expect_equal(joint$estimate, separate, tolerance = 1e-6)
})

test_that("omitting newdata averages each model over its own sample", {
  d <- usability_data()
  d$z[1:60] <- NA
  m1 <- glm(y ~ x, binomial, d)
  m2 <- glm(y ~ x + z, binomial, d)
  fit <- suest(m1, m2)
  default <- marginaleffects::avg_comparisons(fit, variables = "x")
  own <- marginaleffects::avg_comparisons(fit, variables = "x", newdata = suest_newdata(fit))
  expect_equal(default$estimate, own$estimate, tolerance = 1e-12)
  expect_equal(default$std.error, own$std.error, tolerance = 1e-12)

  same <- suest(glm(y ~ x, binomial, d), glm(y ~ x + g, binomial, d))
  default <- marginaleffects::avg_comparisons(same, variables = "x")
  own <- marginaleffects::avg_comparisons(same, variables = "x",
    newdata = suest_newdata(same))
  expect_equal(default$estimate, own$estimate, tolerance = 1e-12)
  expect_equal(default$std.error, own$std.error, tolerance = 1e-12)
})

test_that("datagrid() builds grids from every model's variables", {
  d <- usability_data()
  fit <- suest(glm(y ~ x, binomial, d), glm(y ~ x + z, binomial, d),
    model_names = c("Base", "Full"))
  grid <- marginaleffects::predictions(fit,
    newdata = marginaleffects::datagrid(x = c(0, 1)))
  direct <- marginaleffects::predictions(fit,
    newdata = data.frame(x = c(0, 1), z = mean(d$z)))
  expect_equal(grid$estimate, direct$estimate, tolerance = 1e-12)
  expect_equal(grid$std.error, direct$std.error, tolerance = 1e-12)

  sub <- suest(glm(y ~ x, binomial, d, subset = g == 1),
    glm(y ~ x + z, binomial, d, subset = g != 1),
    model_names = c("One", "Other"))
  expect_error(
    marginaleffects::predictions(sub,
      newdata = marginaleffects::datagrid(x = c(0, 1))),
    'add .suest_model = c\\("One", "Other"\\)'
  )
  both <- marginaleffects::predictions(sub, newdata = marginaleffects::datagrid(
    x = c(0, 1), .suest_model = c("One", "Other")))
  expect_identical(nrow(both), 4L)
})

test_that("survreg models written with a bare Surv() share their data", {
  skip_if_not_installed("survival")
  d <- usability_data()
  d$time <- rexp(400, exp(-.3 * d$x)) + .1
  d$event <- rbinom(400, 1, .8)
  a <- survival::survreg(Surv(time, event) ~ x, data = d)
  b <- survival::survreg(Surv(time, event) ~ x + z, data = d)
  fit <- suest(a, b)
  expect_equal(fit$nobs_overlap, 400L)
  expect_equal(fit$nobs_union, 400L)
  a2 <- survival::survreg(survival::Surv(time, event) ~ x, data = d)
  b2 <- survival::survreg(survival::Surv(time, event) ~ x + z, data = d)
  expect_equal(unname(vcov(fit)), unname(vcov(suest(a2, b2))), tolerance = 1e-12)
})

test_that("zero-inflated negative binomial fits do not need x = TRUE", {
  skip_if_not_installed("pscl")
  d <- usability_data()
  d$count <- ifelse(runif(400) < plogis(-1 + .5 * d$z), 0L,
                    rnbinom(400, size = 2, mu = exp(.5 + .4 * d$x)))
  plain <- suest(pscl::zeroinfl(count ~ x | z, d, dist = "negbin"),
                 pscl::zeroinfl(count ~ x + z | z, d, dist = "negbin"))
  stored <- suest(pscl::zeroinfl(count ~ x | z, d, dist = "negbin", x = TRUE),
                  pscl::zeroinfl(count ~ x + z | z, d, dist = "negbin", x = TRUE))
  expect_equal(unname(vcov(plain)), unname(vcov(stored)), tolerance = 1e-12)
})

test_that("own-sample averaging works with matrix columns in the model frame", {
  skip_if_not_installed("survival")
  d <- usability_data()
  d$time <- rexp(400, exp(-.3 * d$x)) + .1
  d$event <- rbinom(400, 1, .8)
  sv <- suest(survival::survreg(survival::Surv(time, event) ~ x, data = d),
              survival::survreg(survival::Surv(time, event) ~ x + z, data = d))
  own <- marginaleffects::avg_slopes(sv, variables = "x", newdata = suest_newdata(sv))
  common <- marginaleffects::avg_slopes(sv, variables = "x", newdata = d)
  expect_equal(own$estimate, common$estimate, tolerance = 1e-10)

  pf <- suest(glm(y ~ poly(x, 2), binomial, d), glm(y ~ poly(x, 2) + z, binomial, d))
  own <- marginaleffects::avg_slopes(pf, variables = "x")
  common <- marginaleffects::avg_slopes(pf, variables = "x", newdata = d)
  expect_equal(own$estimate, common$estimate, tolerance = 1e-10)
})
