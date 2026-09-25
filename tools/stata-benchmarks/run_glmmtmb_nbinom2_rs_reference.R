# Run from this directory against the current development source.
pkgload::load_all("../..", quiet = TRUE)
dat <- read.csv("suest_r_glmmtmb_nbinom2_rs_crosslang_benchmark.csv"); dat$id <- factor(dat$id)
outcomes <- c("y1", "y2", "y1_partial", "y2_partial", "y1_left", "y2_right")
fit_model <- function(y) {
  d <- dat[!is.na(dat[[y]]), ]
  model <- glmmTMB::glmmTMB(reformulate(c("x", "z", "(1 + x | id)"), y), data = d,
    family = glmmTMB::nbinom2("log"))
  suest:::.suest_model_adapter(model)
  S <- as.matrix(glmmTMB::VarCorr(model)$cond$id)
  # Fixture-specific preflight, deliberately stronger than package restrictions.
  stopifnot(min(eigen(S, symmetric = TRUE)$values) > .02, abs(cov2cor(S)[1, 2]) < .85,
    sigma(model) > .5, sigma(model) < 20, max(sqrt(diag(vcov(model, full = TRUE)))) < .8)
  model
}
extract_model <- function(model) {
  part <- suest:::.suest_model_components(model, "glmm_nbinom2_rs", "glmmTMB::glmmTMB")
  b <- part$parameters; V <- part$bread/nobs(model)
  X <- model.matrix(model, component = "cond"); x <- model$frame$x
  S <- as.matrix(glmmTMB::VarCorr(model)$cond$id)
  A <- S[1, 1]; B <- S[2, 2]; C <- S[1, 2]; E <- sqrt(A*B)*(1-cov2cor(S)[1, 2]^2)
  mu <- exp(drop(X %*% b[colnames(X)])+(A+2*C*x+B*x^2)/2)
  gradient <- mu*cbind(X, 0, A+C*x, B*x^2+C*x, E*x)
  multiplier <- b["x"]+C+B*x
  direct <- cbind(0, 1, 0, 0, C, C+2*B*x, E)
  G <- rbind(colMeans(gradient), colMeans(gradient*multiplier+mu*direct))
  list(coefficients = b, vcov = V, logLik = as.numeric(logLik(model)),
    nobs = nobs(model), n_groups = nlevels(droplevels(model$frame$id)), random_covariance = S,
    native_margins = list(estimates = c(mean = mean(mu), slope_x = mean(mu*multiplier)), vcov = G %*% V %*% t(G)))
}
models <- setNames(lapply(outcomes, fit_model), outcomes)
# Independent optimizer check: target the same Laplace likelihood.
checks <- do.call(rbind, lapply(outcomes, function(y) {
  m <- models[[y]]
  second <- update(m, data = dat[!is.na(dat[[y]]), ], control = glmmTMB::glmmTMBControl(optimizer = optim,
    optArgs = list(method = "BFGS"), optCtrl = list(reltol = 1e-11, maxit = 500)))
  suest:::.suest_model_adapter(second)
  a <- extract_model(m); b <- extract_model(second)
  result <- data.frame(outcome = y, logLik = a$logLik,
    min_covariance_eigenvalue = min(eigen(a$random_covariance, symmetric = TRUE)$values),
    rho = cov2cor(a$random_covariance)[1, 2], phi = sigma(m),
    optimizer_loglik_gap = abs(a$logLik-b$logLik),
    optimizer_parameter_gap = max(abs(a$coefficients-b$coefficients)))
  stopifnot(result$optimizer_loglik_gap < 1e-5, result$optimizer_parameter_gap < 1e-4)
  result
}))
systems <- list(
  balanced = suest(models$y1, models$y2, model_names = c("Y1", "Y2"), observation_id = c("id", "time")),
  partial_higher = suest(models$y1_partial, models$y2_partial, model_names = c("Y1", "Y2"),
    observation_id = c("id", "time"), cluster = "higher"),
  disjoint = suest(models$y1_left, models$y2_right, model_names = c("Y1", "Y2"), observation_id = c("id", "time")))
r <- list(metadata = list(R = R.version.string, suest = as.character(packageVersion("suest")),
  glmmTMB = as.character(packageVersion("glmmTMB")), marginaleffects = as.character(packageVersion("marginaleffects")),
  integration = "Laplace", nuisance = c("log_phi", "log_sd_intercept", "log_sd_slope", "atanh_rho")),
  components = lapply(models, extract_model), preflight = checks,
  systems = lapply(systems, function(f) list(coefficients = coef(f), vcov = vcov(f),
    nobs_union = f$nobs_union, nobs_overlap = f$nobs_overlap, n_clusters = f$n_clusters)),
  predictions = marginaleffects::avg_predictions(systems$balanced, newdata = dat),
  slopes = marginaleffects::avg_slopes(systems$balanced, variables = "x", newdata = dat))
saveRDS(r, "glmmtmb_nbinom2_rs_r_reference.rds")
write.csv(checks, "glmmtmb_nbinom2_rs_preflight.csv", row.names = FALSE)
# R seeds are explicitly identified, then verified by zero-step Stata likelihoods.
# Natural Stata order: beta_x, beta_z, intercept, lnalpha, var_slope, var_intercept, covariance.
seed <- lapply(r$components[c("y1", "y2")], function(m) {
  b <- m$coefficients; S <- m$random_covariance
  unname(c(b[c("x", "z", "(Intercept)")], -b["log_phi"], S[2, 2], S[1, 1], S[1, 2]))
})
fmt <- function(x) paste(format(x, digits = 17, scientific = FALSE, trim = TRUE), collapse = ", ")
writeLines(c("* Generated from checked R Laplace fits; these are starting values, not Stata estimates.",
  paste0("matrix nb2rs_r_starts = (", fmt(seed[[1]]), " \\ ///"), paste0("    ", fmt(seed[[2]]), ")"),
  "matrix rownames nb2rs_r_starts = y1 y2",
  "matrix colnames nb2rs_r_starts = x z _cons lnalpha var_slope var_intercept covariance",
  paste0("scalar nb2rs_r_y1_ll = ", fmt(r$components$y1$logLik)),
  paste0("scalar nb2rs_r_y2_ll = ", fmt(r$components$y2$logLik))), "suest_r_glmmtmb_nbinom2_rs_start_values.do")
print(checks, digits = 10, row.names = FALSE)
cat("GLMMTMB_NB2_RS_REFERENCE_COMPLETE=1\n")
