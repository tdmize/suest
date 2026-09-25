for (case in c("balanced", "partial", "disjoint")) test_that(paste("logit random-slope", case,
  "Stata points and audited covariance agree"), {
  skip_if_not_installed("glmmTMB", minimum_version = "1.1.14")
  f <- readRDS(test_path("fixtures", paste0("glmmtmb-logit-rs-stata-", case, "-v1.rds")))
  audit <- readRDS(test_path("fixtures", paste0("glmmtmb-logit-rs-", case, "-audit.rds")))
  d <- f$data; d$id <- factor(d$id)
  models <- lapply(names(f$components), function(y) glmmTMB::glmmTMB(
    reformulate(c("x", "z", "(1 + x | id)"), y), data = d[!is.na(d[[y]]), ], family = binomial("logit")))
  expect_equal(vapply(models, nobs, numeric(1)),
    switch(case, balanced = c(1200, 1200), partial = c(1175, 1180), disjoint = c(600, 600)))
  for (i in 1:2) {
    model <- models[[i]]; s <- f$components[[i]]
    part <- suest:::.suest_model_components(model, "glmm_logit_rs", "glmmTMB::glmmTMB")
    expect_identical(names(part$parameters), names(s$coefficients))
    expect_lt(max(abs(part$parameters-s$coefficients)), 5e-5)
    expect_lt(abs(as.numeric(logLik(model))-s$logLik), 5e-6)
    V <- part$bread/nobs(model)
    # The returned component covariance agrees directly, unlike raw joint GSEM.
    expect_lt(max(abs(V-s$vcov)), 1e-6)
    expect_lt(max(abs(sqrt(diag(V)/diag(s$vcov))-1)), 1e-4)
    B <- audit$R$components[[i]]$bread_native_delta
    expect_lt(max(abs(V-B)), 1e-6)
    expect_lt(max(abs(sqrt(diag(V)/diag(B))-1)), 1e-4)
    native <- glmmtmb_logit_rs_native_margins(s$coefficients, s$vcov, model$frame)
    expect_lt(max(abs(vapply(s$native_margins, `[[`, numeric(1), "estimate")-native$estimates)), 1e-8)
    returned_variance <- vapply(s$native_margins, `[[`, numeric(1), "variance")
    expect_lt(max(abs(sqrt(returned_variance/diag(native$vcov))-1)), 1e-4)
  }
  fit <- suest(models[[1]], models[[2]], model_names = c("Y1", "Y2"), observation_id = c("id", "time"),
    cluster = if (case == "partial") "higher" else NULL)
  expect_identical(names(coef(fit)), names(f$joint$coefficients))
  expect_lt(max(abs(coef(fit)-f$joint$coefficients)), 5e-5)
  # The raw GSEM covariance remains in f; this is a separate independent audit.
  independent <- audit$R$native_coordinate_vcov
  expect_lt(max(abs(vcov(fit)-independent)), 1e-6)
  expect_lt(max(abs(sqrt(diag(vcov(fit))/diag(independent))-1)), 1e-4)
  expect_identical(fit$n_clusters, if (case == "partial") 20L else 100L)
  expect_equal(fit$nobs_union, 1200)
  expect_equal(fit$nobs_overlap, switch(case, balanced = 1200, partial = 1155, disjoint = 0))
  if (case == "disjoint") {
    expect_equal(unname(audit$R$metadata$normalization), (50/49)/(100/99), tolerance = 1e-14)
    expect_equal(unname(audit$R$metadata$component_clusters), c(50L, 50L))
    expect_equal(unname(vcov(fit)[1:6, 7:12]), matrix(0, 6, 6), tolerance = 0)
    expect_equal(unname(independent[1:6, 7:12]), matrix(0, 6, 6), tolerance = 0)
  }
  # Reconstruct raw Stata V using its supplied bread and independent score
  # cross-products, evaluated at the actual Stata joint parameter point.
  Bs <- f$joint$modelbased; M <- audit$Stata$independent_meat
  G <- f$joint$n_clusters
  reconstructed <- Bs %*% M %*% Bs*G/(G-1)
  raw_gap <- max(abs(reconstructed-f$joint$vcov))
  if (case == "disjoint") {
    # This returned raw-Stata diagnostic remains outside the absolute bound.
    # Do not turn an independent-R pass into a claim of raw-Stata parity.
    expect_identical(audit$validation$raw_reconstruction_absolute_pass, raw_gap < 1e-6)
    expect_false(audit$validation$raw_reconstruction_absolute_pass)
    expect_true(audit$validation$independent_R_pass)
  } else expect_lt(raw_gap, 1e-6)
  expect_lt(max(abs(sqrt(diag(reconstructed)/diag(f$joint$vcov))-1)), 1e-4)
})
