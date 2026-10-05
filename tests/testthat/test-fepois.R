fepois_test_data <- function() {
  set.seed(1009)
  groups <- 50L
  size <- 5L
  d <- data.frame(id = rep(seq_len(groups), each = size), time = rep(seq_len(size), groups),
    x = rnorm(groups*size), z = rnorm(groups*size))
  a <- rep(rnorm(groups), each = size)
  d$y <- rpois(nrow(d), exp(-0.5 + a + 0.4*d$x - 0.2*d$z))
  d$y2 <- rpois(nrow(d), exp(-0.2 + a + 0.3*d$x))
  d
}

test_that("fixed-effects Poisson blocks equal xtpoisson, fe vce(robust)", {
  skip_if_not_installed("fixest")
  d <- fepois_test_data()
  first <- suppressMessages(fixest::fepois(y ~ x + z | id, data = d))
  second <- suppressMessages(fixest::fepois(y2 ~ x | id, data = d))
  fit <- suest(first, second, model_names = c("A", "B"))
  expect_identical(unname(fit$model_types), rep("panel_poisson_fe", 2L))
  # Stata's native xtpoisson, fe vce(robust) has no G/(G-1) factor.
  clustered <- function(model) {
    out <- stats::vcov(model, cluster = ~id,
      ssc = fixest::ssc(adj = FALSE, cluster.adj = FALSE))
    matrix(out, nrow(out))
  }
  expect_equal(unname(vcov(fit)[1:2, 1:2]), unname(clustered(first)), tolerance = 1e-10)
  expect_equal(unname(vcov(fit)[3, 3]), unname(clustered(second))[1, 1], tolerance = 1e-10)

  # Cross block: inverse information around the unit-summed scores.
  bread <- function(model) solve(model$hessian)
  scores <- function(model) rowsum(sandwich::estfun(model), d$id[fixest::obs(model)])
  s1 <- scores(first)
  s2 <- scores(second)
  shared <- intersect(rownames(s1), rownames(s2))
  cross <- bread(first) %*%
    crossprod(s1[shared, , drop = FALSE], s2[shared, , drop = FALSE]) %*% bread(second)
  expect_equal(unname(vcov(fit)[1:2, 3, drop = FALSE]), unname(cross), tolerance = 1e-10)
})

test_that("fixed-effects Poisson predictions include each unit's fixed effect", {
  skip_if_not_installed("fixest")
  d <- fepois_test_data()
  model <- suppressMessages(fixest::fepois(y ~ x + z | id, data = d))
  fit <- suest(model, model, model_names = c("A", "B"))
  rows <- d[fixest::obs(model)[1:4], ]
  predicted <- marginaleffects::predictions(fit, newdata = rows)
  expect_equal(predicted$estimate[predicted$group == "A"],
    unname(stats::predict(model, newdata = rows)), tolerance = 1e-10)
  slopes <- marginaleffects::avg_slopes(fit, variables = "x")
  expect_equal(slopes$estimate[1], unname(coef(model)["x"])*mean(model$y), tolerance = 1e-6)
  link <- marginaleffects::avg_slopes(fit, variables = "x", type = "link")
  expect_equal(link$estimate[1], unname(coef(model)["x"]), tolerance = 1e-6)
  expect_equal(link$std.error[1], sqrt(vcov(fit)[1, 1]), tolerance = 1e-5)
})

test_that("unsupported fixed-effects Poisson models are refused", {
  skip_if_not_installed("fixest")
  d <- fepois_test_data()
  two <- suppressMessages(fixest::fepois(y ~ x | id + time, data = d))
  expect_error(suest(two, two), "exactly one plain fixed effect")
  d$w <- 1 + (d$id %% 3)
  weighted <- suppressMessages(fixest::fepois(y ~ x | id, data = d, weights = ~w))
  expect_error(suest(weighted, weighted), "Weighted")
  plain <- suppressMessages(fixest::fepois(y ~ x | id, data = d))
  expect_error(suest(plain, glm(y ~ x, poisson, d)), "same type")
})
