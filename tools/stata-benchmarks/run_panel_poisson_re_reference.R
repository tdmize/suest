suppressPackageStartupMessages(library(suest))
maxLik <- maxLik::maxLik

dat <- read.csv("suest_r_panel_poisson_re_crosslang_benchmark.csv")
fit_re_poisson <- function(formula) {
  start <- c(coef(glm(formula, data = dat, family = poisson("log"))),
    alpha = 0.5)
  suppressWarnings(pglm::pglm(
    formula, data = dat, family = poisson("log"), model = "random",
    effect = "individual", other = "sd", index = c("id", "time"),
    method = "bfgs", print.level = 0, reltol = 1e-12, start = start
  ))
}

transform_model <- function(model) {
  parameters <- coef(model)
  names(parameters)[length(parameters)] <- "alpha"
  covariance <- vcov(model)
  dimnames(covariance) <- list(names(parameters), names(parameters))
  jacobian <- diag(length(parameters))
  jacobian[nrow(jacobian), ncol(jacobian)] <- 1/parameters["alpha"]
  parameters["alpha"] <- log(parameters["alpha"])
  names(parameters)[length(parameters)] <- "lnalpha"
  covariance <- jacobian %*% covariance %*% t(jacobian)
  dimnames(covariance) <- list(names(parameters), names(parameters))
  list(
    coefficients = parameters,
    vcov = covariance,
    logLik = as.numeric(logLik(model))
  )
}

transform_system <- function(fit) {
  parameters <- coef(fit)
  alpha <- grep("::alpha$", names(parameters))
  jacobian <- diag(length(parameters))
  jacobian[cbind(alpha, alpha)] <- 1/parameters[alpha]
  parameters[alpha] <- log(parameters[alpha])
  names(parameters)[alpha] <- sub("alpha$", "lnalpha", names(parameters)[alpha])
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

y1 <- fit_re_poisson(y1 ~ x + z)
y2 <- fit_re_poisson(y2 ~ x + z)
balanced <- suest(
  y1, y2, model_names = c("Y1", "Y2"),
  observation_id = c("id", "time")
)
balanced_predictions <- marginaleffects::avg_predictions(balanced, newdata = dat)
balanced_slopes <- marginaleffects::avg_slopes(
  balanced, variables = "x", newdata = dat
)

y1_partial <- fit_re_poisson(y1_partial ~ x + z)
y2_partial <- fit_re_poisson(y2_partial ~ x + z)
partial <- suest(
  y1_partial, y2_partial, model_names = c("Y1", "Y2"),
  observation_id = c("id", "time"), cluster = "higher"
)

y1_left <- fit_re_poisson(y1_left ~ x + z)
y2_right <- fit_re_poisson(y2_right ~ x + z)
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
    mixing = "multiplicative gamma random effect; natural-scale variance alpha"
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
saveRDS(reference, "panel_poisson_re_r_reference.rds")

cat("PANEL_POISSON_RE_R_REFERENCE_COMPLETE=1\n")
print(reference$metadata)
for (case in names(reference)[-1L]) {
  cat("\n===", case, "===\n")
  print(reference[[case]]$system)
}
