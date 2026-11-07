test_that("default slopes and comparisons skip the outcomes", {
  set.seed(1)
  sim <- data.frame(x = rnorm(300), z = rbinom(300, 1, 0.5))
  sim$m <- 0.5 * sim$x + rnorm(300)
  sim$y1 <- rbinom(300, 1, plogis(0.4 * sim$x + 0.6 * sim$m))
  sim$y2 <- rbinom(300, 1, plogis(0.3 * sim$x - 0.4 * sim$z))
  sim$amount <- ifelse(sim$y2 == 1, rgamma(300, 2, 1), 0)

  fit <- suest(glm(y1 ~ x + z, family = binomial(), data = sim),
               glm(y2 ~ x + z, family = binomial(), data = sim))
  expect_no_warning(s <- marginaleffects::avg_slopes(fit))
  expect_setequal(unique(s$term), c("x", "z"))
  expect_setequal(insight::find_response(fit), c("y1", "y2"))
  expect_setequal(insight::find_predictors(fit)$conditional, c("x", "z"))

  med <- suest(lm(m ~ x, data = sim), glm(y1 ~ x + m, family = binomial(), data = sim))
  expect_setequal(unique(marginaleffects::avg_slopes(med)$term), c("x", "m"))
  expect_equal(insight::find_response(med), "y1")

  two <- suest(glm(I(amount > 0) ~ x + z, family = binomial(), data = sim),
               glm(amount ~ x + z, family = quasipoisson(), data = sim, subset = amount > 0))
  expect_setequal(unique(marginaleffects::avg_slopes(two)$term), c("x", "z"))

  one <- suest(glm(y1 ~ 1, family = binomial(), data = sim),
               glm(y2 ~ x, family = binomial(), data = sim))
  expect_equal(unique(marginaleffects::avg_slopes(one)$term), "x")
})

test_that("default slopes skip the random-effect grouping variable", {
  skip_if_not_installed("lme4")
  dat <- mtcars
  dat$high <- as.integer(dat$mpg > 20)
  fit <- suest(lme4::glmer(high ~ wt + (1 | cyl), family = binomial(), data = dat),
               glm(am ~ wt, family = binomial(), data = dat))
  expect_equal(unique(marginaleffects::avg_slopes(fit)$term), "wt")
})

test_that("default slopes skip instruments, design variables, and fixed-effect ids", {
  skip_if_not_installed("fixest")
  skip_if_not_installed("survey")
  set.seed(2)
  n <- 400
  d <- data.frame(id = seq_len(n), unit = rep(1:40, each = 10), strata = rep(1:4, each = 100),
                  psu = rep(1:40, each = 10), z = rnorm(n), inst = rnorm(n))
  d$endo <- 0.6 * d$inst + 0.3 * d$z + rnorm(n)
  d$y <- 1 + 0.5 * d$endo - 0.4 * d$z + rnorm(n)
  d$cnt <- rpois(n, exp(0.2 + 0.3 * d$z))
  d$w <- 1 + (d$id %% 5) / 4
  iv <- suest(fixest::feols(y ~ z | endo ~ inst, data = d),
              fixest::feols(cnt ~ z | endo ~ inst, data = d))
  expect_setequal(unique(marginaleffects::avg_slopes(iv)$term), c("endo", "z"))
  fe <- suest(fixest::fepois(cnt ~ z | unit, data = d), fixest::fepois(cnt ~ z + inst | unit, data = d))
  expect_setequal(unique(marginaleffects::avg_slopes(fe)$term), c("z", "inst"))
  des <- survey::svydesign(~psu, strata = ~strata, weights = ~w, data = d)
  sv <- suest(survey::svyglm(y ~ z, design = des), survey::svyivreg(y ~ endo + z | inst + z, des),
              survey_design = des, observation_id = "id")
  expect_setequal(insight::find_predictors(sv)$conditional, c("z", "endo"))
})

test_that("tidy() and glance() give one row per parameter and the combined N", {
  dat <- mtcars
  m1 <- glm(am ~ wt + hp, family = binomial(), data = dat)
  m2 <- lm(mpg ~ wt, data = dat)
  fit <- suest(m1, m2, model_names = c("Logit", "Linear"))
  td <- generics::tidy(fit)
  expect_equal(td$term, c(names(coef(m1)), names(coef(m2)), "lnvar"))
  expect_equal(td$group, rep(c("Logit", "Linear"), c(3, 3)))
  expect_equal(td$estimate, unname(coef(fit)))
  expect_equal(td$std.error, unname(sqrt(diag(vcov(fit)))))
  expect_equal(td$p.value, 2 * stats::pnorm(-abs(td$estimate / td$std.error)))
  expect_equal(td$conf.high - td$estimate, stats::qnorm(0.975) * td$std.error)
  expect_false("conf.low" %in% names(generics::tidy(fit, conf.int = FALSE)))
  expect_equal(generics::glance(fit)$nobs, 32)

  skip_if_not_installed("modelsummary")
  tab <- modelsummary::modelsummary(fit, shape = term ~ group, output = "data.frame")
  cols <- c(grep("Logit$", names(tab)), grep("Linear$", names(tab)))
  expect_length(cols, 2)
  est <- tab[tab$term == "wt" & tab$statistic == "estimate", cols]
  expect_equal(as.numeric(unlist(est)),
               round(c(coef(m1)["wt"], coef(m2)["wt"]), 3), ignore_attr = TRUE)
})

test_that("tidy() for multiply imputed systems uses the pooled results", {
  skip_if_not_installed("mice")
  dat <- mtcars
  set.seed(1)
  dat$wt[sample(32, 5)] <- NA
  imp <- mice::mice(dat[c("am", "vs", "wt", "hp")], m = 3, printFlag = FALSE, seed = 1)
  fits <- lapply(seq_len(3), function(i) {
    d <- mice::complete(imp, i)
    suest(glm(am ~ wt, family = binomial(), data = d),
          glm(vs ~ wt, family = binomial(), data = d), model_names = c("A", "B"))
  })
  fit <- suest_mi(fits)
  td <- generics::tidy(fit)
  s <- summary(fit)
  expect_equal(td$estimate, s$estimate)
  expect_equal(td$std.error, s$std.error)
  expect_equal(td$group, rep(c("A", "B"), each = 2))
  expect_equal(td$term, rep(c("(Intercept)", "wt"), 2))
  expect_equal(generics::glance(fit)$nimp, 3)
})
