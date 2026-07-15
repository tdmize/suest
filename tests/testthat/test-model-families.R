test_that("negative binomial includes dispersion", {
  model1 <- MASS::glm.nb(breaks ~ wool, data = warpbreaks)
  model2 <- MASS::glm.nb(breaks ~ wool + tension, data = warpbreaks)

  fit <- suest(model1, model2)

  expect_true(any(grepl("ln_theta$", names(coef(fit)))))
  expect_true(all(is.finite(vcov(fit))))
})

test_that("ordered and multinomial models return category probabilities", {
  dat <- MASS::housing[
    rep(seq_len(nrow(MASS::housing)), MASS::housing$Freq),
    c("Sat", "Infl", "Type", "Cont")
  ]
  rownames(dat) <- NULL
  dat$Sat_ord <- ordered(dat$Sat, levels = c("Low", "Medium", "High"))
  dat$Sat_nom <- factor(dat$Sat, levels = c("Low", "Medium", "High"))

  ordered_model <- MASS::polr(
    Sat_ord ~ Cont + Infl + Type,
    data = dat,
    method = "logistic",
    Hess = TRUE
  )
  nominal_model <- nnet::multinom(
    Sat_nom ~ Cont + Infl + Type,
    data = dat,
    Hess = TRUE,
    trace = FALSE
  )

  fit <- suest(
    ordered_model,
    nominal_model,
    model_names = c("Ordered", "Multinomial")
  )

  predictions <- marginaleffects::avg_predictions(fit, newdata = dat)
  totals <- aggregate(
    predictions$estimate,
    list(model = sub(
      "::.*$",
      "",
      as.character(predictions$group)
    )),
    sum
  )

  expect_equal(totals$x, c(1, 1), tolerance = 1e-10)
})

test_that("allowed mixed scalar families work", {
  dat <- mtcars
  dat$am <- factor(dat$am)

  logit <- glm(am ~ wt, family = binomial("logit"), data = dat)
  probit <- glm(am ~ wt, family = binomial("probit"), data = dat)
  linear <- lm(as.numeric(as.character(am)) ~ wt, data = dat)

  expect_s3_class(suest(logit, probit), "suest_model")
  expect_s3_class(suest(logit, linear), "suest_model")
  expect_s3_class(suest(probit, linear), "suest_model")
})
