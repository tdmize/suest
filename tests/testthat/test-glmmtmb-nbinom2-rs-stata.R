for (case in c("balanced", "partial", "disjoint")) test_that(paste("NB2 random-slope", case,
  "Stata points and independently audited covariance agree"), {
  skip_if_not_installed("glmmTMB", minimum_version = "1.1.14")
  suffix <- switch(case, balanced = "balanced-v1", partial = "partial-v2", disjoint = "disjoint-v3")
  f <- readRDS(test_path("fixtures", paste0("glmmtmb-nbinom2-rs-stata-", suffix, ".rds")))
  audit <- readRDS(test_path("fixtures", paste0("glmmtmb-nbinom2-rs-", case, "-audit.rds")))
  dat <- f$data; dat$id <- factor(dat$id)
  models <- lapply(names(f$components), function(y) glmmTMB::glmmTMB(
    reformulate(c("x", "z", "(1 + x | id)"), y), data = dat[!is.na(dat[[y]]), ], family = glmmTMB::nbinom2("log")))
  expect_equal(vapply(models, nobs, numeric(1)),
    switch(case, balanced = c(1000, 1000), partial = c(975, 980), disjoint = c(500, 500)))
  for (i in 1:2) {
    model <- models[[i]]; s <- f$components[[i]]
    component <- suest:::.suest_model_components(model, "glmm_nbinom2_rs", "glmmTMB::glmmTMB")
    expect_identical(names(component$parameters), names(s$coefficients))
    delta <- component$parameters-s$coefficients
    if (case != "partial") expect_lt(max(abs(delta)), 5e-5) else {
      # Partial y2's log-size endpoint differs by 7.28e-5 (0.061% of its SE).
      expect_lt(max(abs(delta[-4L])), 5e-5)
      expect_lt(max(abs(delta)/sqrt(diag(component$bread/nobs(model)))), 1e-3)
    }
    expect_lt(abs(as.numeric(logLik(model))-s$logLik), 5e-6)
    # Component OIM from Stata differs numerically. Compare R's native
    # curvature with the separate Laplace finite-difference audit instead.
    B <- audit$R$components[[i]]$bread_native_delta
    V <- component$bread/nobs(model)
    expect_lt(max(abs(V-B)), 1e-6)
    expect_lt(max(abs(sqrt(diag(V)/diag(B))-1)), 1e-4)
    b <- s$coefficients; x <- model$frame$x; X <- model.matrix(model, component = "cond")
    A <- exp(2*b[5]); B <- exp(2*b[6]); rho <- tanh(b[7]); C <- sqrt(A*B)*rho
    E <- sqrt(A*B)*(1-rho^2)
    mu <- exp(drop(X %*% b[1:3])+(A+2*C*x+B*x^2)/2)
    G <- mu*cbind(X, 0, A+C*x, B*x^2+C*x, E*x)
    multiplier <- b[2]+C+B*x
    H <- G*multiplier+mu*cbind(0, 1, 0, 0, C, C+2*B*x, E)
    gradients <- rbind(colMeans(G), colMeans(H))
    estimates <- c(mean(mu), mean(mu*multiplier))
    variances <- diag(gradients %*% s$vcov %*% t(gradients))
    # Independently verify the printed native Stata margins and their delta
    # variances using Stata's own supplied full covariance, including nuisance.
    expect_equal(unname(vapply(s$native_margins, `[[`, numeric(1), "estimate")),
      unname(estimates), tolerance = 1e-8)
    native_variances <- vapply(s$native_margins, `[[`, numeric(1), "variance")
    if (case != "partial") expect_lt(max(abs(native_variances-variances)), 1e-7) else
      # Native Stata numerical derivatives differ by 1.25e-6 relative SE.
      expect_lt(max(abs(sqrt(native_variances/variances)-1)), 1e-5)
  }
  fit <- suest(models[[1]], models[[2]], model_names = c("Y1", "Y2"), observation_id = c("id", "time"),
    cluster = if (case == "partial") "higher" else NULL)
  expect_identical(names(coef(fit)), names(f$joint$coefficients))
  delta <- coef(fit)-f$joint$coefficients
  if (case != "partial") expect_lt(max(abs(delta)), 5e-5) else {
    expect_lt(max(abs(delta[-c(4L, 11L)])), 5e-5)
    expect_lt(max(abs(delta)/sqrt(diag(vcov(fit)))), 1e-3)
  }
  audited_vcov <- audit$R$native_coordinate_vcov
  expect_lt(max(abs(vcov(fit)-audited_vcov)), 1e-6)
  expect_lt(max(abs(sqrt(diag(vcov(fit))/diag(audited_vcov))-1)), 1e-4)
  expect_identical(fit$n_clusters, if (case == "partial") 20L else 100L)
  expect_equal(fit$nobs_union, 1000)
  expect_equal(fit$nobs_overlap, switch(case, balanced = 1000, partial = 955, disjoint = 0))
  if (case == "disjoint") {
    expect_equal(unname(audit$R$metadata$component_clusters), c(50L, 50L))
    expect_equal(unname(audit$R$metadata$normalization), (50/49)/(100/99), tolerance = 1e-14)
    expect_equal(max(abs(vcov(fit)[1:7, 8:14])), 0)
  }
  # Recovered score CROSS-PRODUCTS, not uniquely recovered score vectors.
  # This decomposition uses the independent audit at Stata's parameter point.
  Bs <- f$joint$modelbased; M <- audit$Stata$independent_meat
  G <- f$joint$n_clusters
  reproduced <- Bs %*% M %*% Bs*G/(G-1)
  expect_lt(max(abs(reproduced-f$joint$vcov)), 1e-6)
  expect_lt(max(abs(sqrt(diag(reproduced)/diag(f$joint$vcov))-1)), 1e-4)
})
