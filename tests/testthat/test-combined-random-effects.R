combined_re_data <- function() {
  set.seed(11)
  groups <- 60L
  size <- 6L
  d <- data.frame(id = rep(seq_len(groups), each = size),
    x = rnorm(groups*size), z = rnorm(groups*size))
  d$school <- (d$id - 1L) %/% 4L
  u <- rep(rnorm(groups, sd = 0.8), each = size)
  d$y1 <- rbinom(nrow(d), 1L, plogis(-0.2 + 0.6*d$x + 0.4*d$z + u))
  d$y2 <- rbinom(nrow(d), 1L, plogis(0.1 + 0.5*d$x + u))
  d$yc <- 1 + 0.5*d$x + 0.3*d$z + u + rnorm(nrow(d))
  d$count <- rpois(nrow(d), exp(0.2 + 0.3*d$x + 0.5*u))
  d$o <- cut(d$yc, c(-Inf, 0, 1.5, Inf), ordered_result = TRUE)
  d
}

block <- function(fit, model) {
  keep <- startsWith(names(coef(fit)), paste0(model, "::"))
  vcov(fit)[keep, keep, drop = FALSE]
}

test_that("an ordinary model combines with a glmmTMB model, clustered on its group", {
  skip_if_not_installed("glmmTMB")
  d <- combined_re_data()
  ordinary <- glm(y1 ~ x + z, binomial, d)
  multilevel <- glmmTMB::glmmTMB(y2 ~ x + (1 | id), family = binomial, data = d)
  fit <- suest(ordinary, multilevel, model_names = c("logit", "melogit"))

  expect_identical(fit$n_clusters, 60L)
  expect_equal(unname(block(fit, "logit")),
    unname(sandwich::vcovCL(ordinary, cluster = ~id, type = "HC0", cadjust = TRUE)),
    tolerance = 1e-10)
  same_type <- suest(multilevel,
    glmmTMB::glmmTMB(y2 ~ x + (1 | id), family = binomial, data = d),
    model_names = c("a", "b"))
  expect_equal(unname(block(fit, "melogit")), unname(block(same_type, "a")),
    tolerance = 1e-10)

  # Cross-model block: model-based covariances around the group-summed scores.
  parts <- suest:::.suest_model_components(multilevel, "glmm_logit_ri",
    "glmmTMB::glmmTMB", NULL)
  scores <- rowsum(parts$score, d$id)
  model_based <- parts$bread/nrow(parts$score)
  cross <- (60/59)*vcov(ordinary) %*%
    crossprod(rowsum(sandwich::estfun(ordinary), d$id), scores) %*% model_based
  keep1 <- startsWith(names(coef(fit)), "logit::")
  expect_equal(unname(vcov(fit)[keep1, !keep1]), unname(cross), tolerance = 1e-9)

  difference <- marginaleffects::avg_comparisons(fit, variables = "x",
    hypothesis = difference ~ revpairwise)
  expect_true(is.finite(difference$std.error) && difference$std.error > 0)
})

test_that("a linear model combines with nlme::lme using regress's cluster adjustment", {
  skip_if_not_installed("nlme")
  d <- combined_re_data()
  linear <- lm(yc ~ x + z, d)
  mixed <- nlme::lme(yc ~ x, random = ~ 1 | id, data = d, method = "ML")
  fit <- suest(linear, mixed, model_names = c("regress", "mixed"))
  beta <- block(fit, "regress")
  beta <- beta[!grepl("lnvar", rownames(beta)), !grepl("lnvar", colnames(beta))]
  expect_equal(unname(beta),
    unname(sandwich::vcovCL(linear, cluster = ~id, type = "HC1")), tolerance = 1e-10)
})

test_that("different random-effects families and ordinal models combine", {
  skip_if_not_installed("glmmTMB")
  d <- combined_re_data()
  logit <- glmmTMB::glmmTMB(y2 ~ x + (1 | id), family = binomial, data = d)
  count <- glmmTMB::glmmTMB(count ~ x + (1 | id), family = poisson, data = d)
  fit <- suest(logit, count, model_names = c("melogit", "mepoisson"))
  alone <- suest(count,
    glmmTMB::glmmTMB(count ~ x + (1 | id), family = poisson, data = d),
    model_names = c("a", "b"))
  expect_equal(unname(block(fit, "mepoisson")), unname(block(alone, "a")),
    tolerance = 1e-10)

  ordered <- MASS::polr(o ~ x, d, method = "logistic", Hess = TRUE)
  combined <- suest(ordered, logit, model_names = c("ologit", "melogit"))
  ordinary <- suest(ordered, glm(y1 ~ x, binomial, d), cluster = "id",
    model_names = c("ologit", "logit"))
  expect_equal(unname(block(combined, "ologit")), unname(block(ordinary, "ologit")),
    tolerance = 1e-10)
  expect_identical(nrow(marginaleffects::avg_comparisons(combined, variables = "x")), 4L)
})

test_that("combined random-effects systems check groups and model types", {
  skip_if_not_installed("glmmTMB")
  d <- combined_re_data()
  multilevel <- glmmTMB::glmmTMB(y2 ~ x + (1 | id), family = binomial, data = d)
  by_school <- glmmTMB::glmmTMB(y1 ~ x + (1 | school), family = binomial, data = d)
  expect_error(suest(multilevel, by_school), "same grouping variable")

  higher <- suest(glm(y1 ~ x, binomial, d), multilevel, cluster = "school")
  expect_identical(higher$n_clusters, length(unique(d$school)))
  expect_error(suest(glm(y1 ~ x, binomial, d), by_school, cluster = "id"),
    "Panel IDs must be nested")

  no_group <- d[, c("y1", "x")]
  expect_error(suppressWarnings(suest(glm(y1 ~ x, binomial, no_group), multilevel)),
    "Cluster-ID column(s) not found for model 'Model 1': id", fixed = TRUE)
  expect_error(suest(glm(count ~ x, Gamma(link = "log"), transform(d, count = count + 1)),
    multilevel), "Random-effects models")
})
