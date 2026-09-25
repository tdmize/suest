suppressPackageStartupMessages(library(suest)) # Joint scores and integrated predictions
dat <- read.csv("suest_r_glmmtmb_poisson_rs_crosslang_benchmark.csv")
dat$id <- factor(dat$id)
fit_model <- function(outcome) {
  model <- glmmTMB::glmmTMB(reformulate(c("x", "z", "(1 + x | id)"), outcome),
    data = dat[!is.na(dat[[outcome]]), ], family = poisson())
  suest:::.suest_model_adapter(model)
  S <- as.matrix(glmmTMB::VarCorr(model)$cond$id)
  # Benchmark preflight is stricter than the package's numerical support check.
  stopifnot(min(eigen(S, symmetric = TRUE)$values) > .02,
    abs(cov2cor(S)[1, 2]) < .85,
    max(sqrt(diag(vcov(model, full = TRUE)))) < .8)
  model
}
extract_model <- function(model) {
  component <- suest:::.suest_model_components(model, "glmm_poisson_rs", "glmmTMB::glmmTMB")
  b <- component$parameters; V <- component$bread/nobs(model)
  X <- model.matrix(model, component = "cond"); x <- model$frame$x
  S <- as.matrix(glmmTMB::VarCorr(model)$cond$id)
  A <- S[1, 1]; B <- S[2, 2]; C <- S[1, 2]
  E <- sqrt(A*B)*(1-cov2cor(S)[1, 2]^2)
  mu <- exp(drop(X %*% b[colnames(X)]) + (A+2*C*x+B*x^2)/2)
  gradient <- mu*cbind(X, A+C*x, B*x^2+C*x, E*x)
  multiplier <- b["x"] + C+B*x
  direct <- cbind(0, 1, 0, C, C+2*B*x, E)
  gradients <- rbind(colMeans(gradient), colMeans(gradient*multiplier + mu*direct))
  list(coefficients = b, vcov = V, logLik = as.numeric(logLik(model)),
    random_covariance = S, nobs = nobs(model), n_groups = nlevels(droplevels(model$frame$id)),
    native_margins = list(estimates = c(mean = mean(mu), slope_x = mean(mu*multiplier)),
      vcov = gradients %*% V %*% t(gradients)))
}
outcomes <- c("y1", "y2", "y1_partial", "y2_partial", "y1_left", "y2_right")
models <- setNames(lapply(outcomes, fit_model), outcomes)
systems <- list(
  balanced = suest(models$y1, models$y2, model_names = c("Y1", "Y2"), observation_id = c("id", "time")),
  partial_higher = suest(models$y1_partial, models$y2_partial, model_names = c("Y1", "Y2"),
    observation_id = c("id", "time"), cluster = "higher"),
  disjoint = suest(models$y1_left, models$y2_right, model_names = c("Y1", "Y2"), observation_id = c("id", "time")))
reference <- list(metadata = list(R = R.version.string, suest = as.character(packageVersion("suest")),
  glmmTMB = as.character(packageVersion("glmmTMB")), marginaleffects = as.character(packageVersion("marginaleffects")),
  fitting_integration = "Laplace", random_structure = "correlated intercept and numeric x slope",
  nuisance_parameters = c("log_sd_intercept", "log_sd_slope", "atanh_rho")),
  components = lapply(models, extract_model),
  systems = lapply(systems, function(fit) list(coefficients = coef(fit), vcov = vcov(fit),
    nobs_union = fit$nobs_union, nobs_overlap = fit$nobs_overlap, n_clusters = fit$n_clusters)),
  predictions = marginaleffects::avg_predictions(systems$balanced, newdata = dat),
  slopes = marginaleffects::avg_slopes(systems$balanced, variables = "x", newdata = dat))
saveRDS(reference, "glmmtmb_poisson_rs_r_reference.rds")
print(reference)
cat("GLMMTMB_POISSON_RS_R_REFERENCE_COMPLETE=1\n")
