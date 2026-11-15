test_that("summary() shows each model's estimates as a panel", {
  dat <- mtcars
  model1 <- glm(am ~ wt, family = binomial(), data = dat)
  model2 <- glm(vs ~ wt + hp, family = binomial(), data = dat)
  fit <- suest(model1, model2)
  s <- summary(fit, conf.level = 0.90)
  td <- tidy.suest_model(fit, conf.level = 0.90)
  expect_s3_class(s, "summary.suest_model")
  expect_equal(s$coefficients, td)
  out <- capture.output(print(s))
  expect_match(out[2], "Simultaneous results for model1, model2 +Number of obs = 32")
  expect_true(all(c("model1        |", "model2        |") %in% out))
  expect_true(any(grepl("[90% conf. interval]", out, fixed = TRUE)))
  expect_equal(sum(grepl("^ +hp \\|", out)), 1L)
  expect_error(summary(fit, conf.level = 1), "strictly between")

  clustered <- suest(model1, model2, cluster = "cyl")
  expect_true(any(grepl("(Std. err. adjusted for 3 clusters)", capture.output(summary(clustered)), fixed = TRUE)))
})
