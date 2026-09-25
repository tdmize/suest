suppressPackageStartupMessages(library(suest)) # Joint score-based covariance and predictions

dat <- read.csv("suest_r_glmmtmb_nbinom2_ri_crosslang_benchmark.csv")
dat$id <- factor(dat$id)
fit_glmm_nbinom2 <- function(formula, data = dat) glmmTMB::glmmTMB(
  formula, data = data, family = glmmTMB::nbinom2("log"))

extract_model <- function(model) {
  parameters <- suest:::.suest_extract_parameters(
    model, "glmm_nbinom2_ri", "glmmTMB::glmmTMB")
  covariance <- as.matrix(vcov(model, full = TRUE))
  dimnames(covariance) <- list(names(parameters), names(parameters))
  X <- model.matrix(model, component = "cond")
  sigma2 <- exp(2*parameters["log_sigma"])
  mu <- exp(drop(X %*% parameters[colnames(X)]) + sigma2/2)
  mean_gradient <- c(colMeans(mu*X), 0, sigma2*mean(mu))
  slope_gradient <- mean_gradient*parameters["x"]
  slope_gradient[match("x", names(parameters))] <-
    slope_gradient[match("x", names(parameters))] + mean(mu)
  gradients <- rbind(mean = mean_gradient, slope_x = slope_gradient)
  list(coefficients = parameters, vcov = covariance,
    logLik = as.numeric(logLik(model)),
    native_margins = list(estimates = c(mean = mean(mu),
      slope_x = unname(parameters["x"])*mean(mu)),
      vcov = gradients %*% covariance %*% t(gradients)))
}

extract_system <- function(fit) list(
  coefficients = coef(fit), vcov = vcov(fit),
  nobs_union = fit$nobs_union, nobs_overlap = fit$nobs_overlap,
  n_clusters = fit$n_clusters)

y1 <- fit_glmm_nbinom2(y1 ~ x + z + (1 | id))
y2 <- fit_glmm_nbinom2(y2 ~ x + z + (1 | id))
balanced <- suest(y1, y2, model_names = c("Y1", "Y2"),
  observation_id = c("id", "time"))
balanced_predictions <- marginaleffects::avg_predictions(
  balanced, newdata = dat)
balanced_slopes <- marginaleffects::avg_slopes(
  balanced, variables = "x", newdata = dat)

y1_partial <- fit_glmm_nbinom2(y1_partial ~ x + z + (1 | id))
y2_partial <- fit_glmm_nbinom2(y2_partial ~ x + z + (1 | id))
partial <- suest(y1_partial, y2_partial, model_names = c("Y1", "Y2"),
  observation_id = c("id", "time"), cluster = "higher")

left_data <- dat[as.integer(dat$id) <= 40, ]
right_data <- dat[as.integer(dat$id) > 40, ]
y1_left <- fit_glmm_nbinom2(y1_left ~ x + z + (1 | id), left_data)
y2_right <- fit_glmm_nbinom2(y2_right ~ x + z + (1 | id), right_data)
disjoint <- suest(y1_left, y2_right, model_names = c("Y1", "Y2"),
  observation_id = c("id", "time"))

reference <- list(
  metadata = list(
    R = R.version.string, suest = as.character(packageVersion("suest")),
    glmmTMB = as.character(packageVersion("glmmTMB")),
    marginaleffects = as.character(packageVersion("marginaleffects")),
    fitting_integration = "Laplace",
    dispersion = "var(Y|b) = mu + mu^2/phi; Stata lnalpha = -log_phi",
    response_prediction = "exact Gaussian marginal mean"
  ),
  balanced = list(
    model1 = extract_model(y1), model2 = extract_model(y2),
    system = extract_system(balanced), predictions = balanced_predictions,
    slopes = balanced_slopes),
  partial_higher = list(
    model1 = extract_model(y1_partial), model2 = extract_model(y2_partial),
    system = extract_system(partial)),
  disjoint = list(
    model1 = extract_model(y1_left), model2 = extract_model(y2_right),
    system = extract_system(disjoint))
)
saveRDS(reference, "glmmtmb_nbinom2_ri_r_reference.rds")
cat("GLMMTMB_NBINOM2_RI_R_REFERENCE_COMPLETE=1\n")
print(reference$metadata)
for (case in names(reference)[-1L]) {
  cat("\n===", case, "===\n")
  print(reference[[case]])
}
