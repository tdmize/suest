suppressPackageStartupMessages(library(suest))
maxLik <- maxLik::maxLik

dat <- read.csv("suest_r_panel_probit_re_crosslang_benchmark.csv")
fit_re_probit <- function(formula) {
  start <- c(coef(glm(formula, data = dat, family = binomial("probit"))), sigma = 1)
  pglm::pglm(
    formula, data = dat, family = binomial("probit"), model = "random",
    effect = "individual", index = c("id", "time"), R = 12,
    method = "bfgs", print.level = 0, reltol = 1e-12, start = start
  )
}

transform_model <- function(model) {
  parameters <- coef(model)
  sigma_sign <- sign(parameters["sigma"])
  parameters["sigma"] <- abs(parameters["sigma"])
  covariance <- vcov(model)
  canonical <- diag(length(parameters))
  canonical[nrow(canonical), ncol(canonical)] <- sigma_sign
  covariance <- t(canonical) %*% covariance %*% canonical
  stata <- diag(length(parameters))
  stata[nrow(stata), ncol(stata)] <- 2/parameters["sigma"]
  parameters["sigma"] <- 2*log(parameters["sigma"])
  names(parameters)[length(parameters)] <- "lnsig2u"
  covariance <- stata %*% covariance %*% t(stata)
  dimnames(covariance) <- list(names(parameters), names(parameters))
  list(coefficients = parameters, vcov = covariance,
    logLik = as.numeric(logLik(model)))
}

transform_system <- function(fit) {
  parameters <- coef(fit)
  sigma <- grep("::sigma$", names(parameters))
  jacobian <- diag(length(parameters))
  jacobian[cbind(sigma, sigma)] <- 2/parameters[sigma]
  parameters[sigma] <- 2*log(parameters[sigma])
  names(parameters)[sigma] <- sub("sigma$", "lnsig2u", names(parameters)[sigma])
  covariance <- jacobian %*% vcov(fit) %*% t(jacobian)
  dimnames(covariance) <- list(names(parameters), names(parameters))
  list(
    coefficients = parameters,
    vcov = covariance,
    nobs_union = fit$nobs_union,
    nobs_overlap = fit$nobs_overlap,
    n_clusters = fit$n_clusters
  )
}

y1 <- fit_re_probit(y1 ~ x + z)
y2 <- fit_re_probit(y2 ~ x + z)
balanced <- suest(
  y1, y2, model_names = c("Y1", "Y2"),
  observation_id = c("id", "time")
)
balanced_predictions <- marginaleffects::avg_predictions(balanced, newdata = dat)
balanced_slopes <- marginaleffects::avg_slopes(
  balanced, variables = "x", newdata = dat
)

y1_partial <- fit_re_probit(y1_partial ~ x + z)
y2_partial <- fit_re_probit(y2_partial ~ x + z)
partial <- suest(
  y1_partial, y2_partial, model_names = c("Y1", "Y2"),
  observation_id = c("id", "time"), cluster = "higher"
)

y1_left <- fit_re_probit(y1_left ~ x + z)
y2_right <- fit_re_probit(y2_right ~ x + z)
disjoint <- suest(
  y1_left, y2_right, model_names = c("Y1", "Y2"),
  observation_id = c("id", "time")
)

reference <- list(
  metadata = list(
    R = R.version.string,
    suest = as.character(packageVersion("suest")),
    pglm = as.character(packageVersion("pglm")),
    marginaleffects = as.character(packageVersion("marginaleffects")),
    integration = "12-point nonadaptive Gauss-Hermite"
  ),
  balanced = list(
    model1 = transform_model(y1),
    model2 = transform_model(y2),
    system = transform_system(balanced),
    predictions = balanced_predictions,
    slopes = balanced_slopes
  ),
  partial_higher = list(
    model1 = transform_model(y1_partial),
    model2 = transform_model(y2_partial),
    system = transform_system(partial)
  ),
  disjoint = list(
    model1 = transform_model(y1_left),
    model2 = transform_model(y2_right),
    system = transform_system(disjoint)
  )
)
saveRDS(reference, "panel_probit_re_r_reference.rds")

cat("PANEL_PROBIT_RE_R_REFERENCE_COMPLETE=1\n")
print(reference$metadata)
for (case in names(reference)[-1L]) {
  cat("\n===", case, "===\n")
  print(reference[[case]]$system)
}
