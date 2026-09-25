# Run from this directory against the current development source.
pkgload::load_all("../..", quiet = TRUE)
dat <- read.csv("suest_r_glmmtmb_logit_rs_crosslang_benchmark.csv"); dat$id <- factor(dat$id)
outcomes <- c("y1", "y2", "y1_partial", "y2_partial", "y1_left", "y2_right")
fit_model <- function(y) {
  d <- dat[!is.na(dat[[y]]), ]
  model <- glmmTMB::glmmTMB(reformulate(c("x", "z", "(1 + x | id)"), y), data = d,
    family = binomial("logit"))
  suest:::.suest_model_adapter(model)
  S <- as.matrix(glmmTMB::VarCorr(model)$cond$id)
  J <- diag(c(rep(1, 5), 1/sqrt(1+unname(tail(model$fit$par, 1))^2)))
  # Assess precision in exported Fisher coordinates, not native scaled correlation.
  V <- J %*% vcov(model, full = TRUE) %*% J
  # Fixture-specific preflight, deliberately stronger than package restrictions.
  stopifnot(min(eigen(S, symmetric = TRUE)$values) > .02, abs(cov2cor(S)[1, 2]) < .85,
    max(sqrt(diag(V))) < .8)
  model
}
extract_model <- function(model) {
  part <- suest:::.suest_model_components(model, "glmm_logit_rs", "glmmTMB::glmmTMB")
  b <- part$parameters; V <- part$bread/nobs(model)
  X <- model.matrix(model, component = "cond"); x <- model$frame$x
  S <- as.matrix(glmmTMB::VarCorr(model)$cond$id)
  # Independent normal integration for component marginal estimates and gradients.
  margins <- function(p) {
    sd0 <- exp(p["log_sd_intercept"]); sd1 <- exp(p["log_sd_slope"]); rho <- tanh(p["atanh_rho"])
    eta <- drop(X %*% p[colnames(X)])
    sd <- sqrt(sd0^2+2*rho*sd0*sd1*x+sd1^2*x^2)
    moments <- vapply(0:2, function(k) vapply(seq_along(x), function(i) integrate(function(u) {
      q <- plogis(eta[i]+sd[i]*u)
      switch(k+1L, q, q*(1-q), q*(1-q)*(1-2*q))*dnorm(u)
    }, -Inf, Inf, rel.tol = 1e-11, abs.tol = 1e-12)$value, numeric(1)), numeric(length(x)))
    c(mean = mean(moments[, 1]), slope_x = mean(p["x"]*moments[, 2]+
      (rho*sd0*sd1+sd1^2*x)*moments[, 3]))
  }
  G <- numDeriv::jacobian(margins, b)
  list(coefficients = b, vcov = V, logLik = as.numeric(logLik(model)),
    nobs = nobs(model), n_groups = nlevels(droplevels(model$frame$id)), random_covariance = S,
    native_margins = list(estimates = margins(b), vcov = G %*% V %*% t(G)))
}
models <- setNames(lapply(outcomes, fit_model), outcomes)
# Independent optimizer check: target the same Laplace likelihood.
checks <- do.call(rbind, lapply(outcomes, function(y) {
  m <- models[[y]]
  second <- update(m, data = dat[!is.na(dat[[y]]), ], control = glmmTMB::glmmTMBControl(optimizer = optim,
    optArgs = list(method = "BFGS"), optCtrl = list(reltol = 1e-11, maxit = 500)))
  suest:::.suest_model_adapter(second)
  extract_preflight <- function(f) list(logLik = as.numeric(logLik(f)),
    coefficients = suest:::.suest_extract_parameters(f, "glmm_logit_rs", "glmmTMB::glmmTMB"),
    random_covariance = as.matrix(glmmTMB::VarCorr(f)$cond$id))
  a <- extract_preflight(m); b <- extract_preflight(second)
  result <- data.frame(outcome = y, logLik = a$logLik,
    min_covariance_eigenvalue = min(eigen(a$random_covariance, symmetric = TRUE)$values),
    rho = cov2cor(a$random_covariance)[1, 2],
    max_native_se = max(sqrt(diag(vcov(m, full = TRUE)))),
    max_exported_se = max(sqrt(diag(diag(c(rep(1, 5), 1/sqrt(1+unname(tail(m$fit$par, 1))^2))) %*%
      vcov(m, full = TRUE) %*% diag(c(rep(1, 5), 1/sqrt(1+unname(tail(m$fit$par, 1))^2)))))),
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
  integration = "Laplace", nuisance = c("log_sd_intercept", "log_sd_slope", "atanh_rho")),
  components = lapply(models, extract_model), preflight = checks,
  systems = lapply(systems, function(f) list(coefficients = coef(f), vcov = vcov(f),
    nobs_union = f$nobs_union, nobs_overlap = f$nobs_overlap, n_clusters = f$n_clusters)),
  predictions = marginaleffects::avg_predictions(systems$balanced, newdata = dat),
  slopes = marginaleffects::avg_slopes(systems$balanced, variables = "x", newdata = dat,
    numderiv = list("fdcenter", eps = 1e-4)))
saveRDS(r, "glmmtmb_logit_rs_r_reference.rds")
write.csv(checks, "glmmtmb_logit_rs_preflight.csv", row.names = FALSE)
# R seeds are explicitly identified, then verified by zero-step Stata likelihoods.
# Natural Stata order: beta_x, beta_z, intercept, var_slope, var_intercept, covariance.
seed <- lapply(r$components[c("y1", "y2")], function(m) {
  b <- m$coefficients; S <- m$random_covariance
  unname(c(b[c("x", "z", "(Intercept)")], S[2, 2], S[1, 1], S[1, 2]))
})
fmt <- function(x) paste(format(x, digits = 17, scientific = FALSE, trim = TRUE), collapse = ", ")
writeLines(c("* Generated from checked R Laplace fits; these are starting values, not Stata estimates.",
  paste0("matrix logitrs_r_starts = (", fmt(seed[[1]]), " \\ ///"), paste0("    ", fmt(seed[[2]]), ")"),
  "matrix rownames logitrs_r_starts = y1 y2",
  "matrix colnames logitrs_r_starts = x z _cons var_slope var_intercept covariance",
  paste0("scalar logitrs_r_y1_ll = ", fmt(r$components$y1$logLik)),
  paste0("scalar logitrs_r_y2_ll = ", fmt(r$components$y2$logLik))), "suest_r_glmmtmb_logit_rs_start_values.do")
print(checks, digits = 10, row.names = FALSE)
cat("GLMMTMB_LOGIT_RS_REFERENCE_COMPLETE=1\n")
