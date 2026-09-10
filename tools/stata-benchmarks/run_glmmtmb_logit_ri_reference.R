suppressPackageStartupMessages(library(suest))

dat <- read.csv("suest_r_glmmtmb_logit_ri_crosslang_benchmark.csv")
dat$id <- factor(dat$id)
fit_glmm_logit <- function(formula, data = dat) {
  glmmTMB::glmmTMB(
    formula, data = data, family = binomial("logit")
  )
}

extract_model <- function(model) {
  parameters <- suest:::.suest_extract_parameters(
    model, "glmm_logit_ri", "glmmTMB::glmmTMB")
  covariance <- as.matrix(vcov(model, full = TRUE))
  dimnames(covariance) <- list(names(parameters), names(parameters))
  list(
    coefficients = parameters,
    vcov = covariance,
    logLik = as.numeric(logLik(model))
  )
}

extract_system <- function(fit) {
  list(
    coefficients = coef(fit),
    vcov = vcov(fit),
    nobs_union = fit$nobs_union,
    nobs_overlap = fit$nobs_overlap,
    n_clusters = fit$n_clusters
  )
}

y1 <- fit_glmm_logit(y1 ~ x + z + (1 | id))
y2 <- fit_glmm_logit(y2 ~ x + z + (1 | id))
balanced <- suest(y1, y2, model_names = c("Y1", "Y2"),
  observation_id = c("id", "time"))
balanced_predictions <- marginaleffects::avg_predictions(
  balanced, newdata = dat)
balanced_slopes <- marginaleffects::avg_slopes(
  balanced, variables = "x", newdata = dat)

y1_partial <- fit_glmm_logit(y1_partial ~ x + z + (1 | id))
y2_partial <- fit_glmm_logit(y2_partial ~ x + z + (1 | id))
partial <- suest(y1_partial, y2_partial, model_names = c("Y1", "Y2"),
  observation_id = c("id", "time"), cluster = "higher")

left_data <- dat[as.integer(dat$id) <= 40, ]
right_data <- dat[as.integer(dat$id) > 40, ]
y1_left <- fit_glmm_logit(y1_left ~ x + z + (1 | id), left_data)
y2_right <- fit_glmm_logit(y2_right ~ x + z + (1 | id), right_data)
disjoint <- suest(y1_left, y2_right, model_names = c("Y1", "Y2"),
  observation_id = c("id", "time"))

reference <- list(
  metadata = list(
    R = R.version.string,
    suest = as.character(packageVersion("suest")),
    glmmTMB = as.character(packageVersion("glmmTMB")),
    marginaleffects = as.character(packageVersion("marginaleffects")),
    fitting_integration = "Laplace",
    response_prediction = "20-point Gaussian-Hermite marginal mean"
  ),
  balanced = list(
    model1 = extract_model(y1),
    model2 = extract_model(y2),
    system = extract_system(balanced),
    predictions = balanced_predictions,
    slopes = balanced_slopes
  ),
  partial_higher = list(
    model1 = extract_model(y1_partial),
    model2 = extract_model(y2_partial),
    system = extract_system(partial)
  ),
  disjoint = list(
    model1 = extract_model(y1_left),
    model2 = extract_model(y2_right),
    system = extract_system(disjoint)
  )
)
saveRDS(reference, "glmmtmb_logit_ri_r_reference.rds")

cat("GLMMTMB_LOGIT_RI_R_REFERENCE_COMPLETE=1\n")
print(reference$metadata)
for (case in names(reference)[-1L]) {
  cat("\n===", case, "===\n")
  print(reference[[case]]$system)
}
