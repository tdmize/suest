test_that("Poisson random-slope components reproduce the returned Stata Laplace fits", {
  skip_if_not_installed("glmmTMB", minimum_version = "1.1.14")
  fixture <- readRDS(test_path("fixtures", "glmmtmb-poisson-rs-stata-components.rds"))
  dat <- fixture$data; dat$id <- factor(dat$id)
  for (y in names(fixture$components)) {
    model <- glmmTMB::glmmTMB(reformulate(c("x", "z", "(1 + x | id)"), y),
      data = dat[!is.na(dat[[y]]), ], family = poisson())
    component <- suest:::.suest_model_components(model, "glmm_poisson_rs", "glmmTMB::glmmTMB")
    stata <- fixture$components[[y]]
    expect_identical(names(component$parameters), names(stata$coefficients))
    expect_lt(max(abs(component$parameters-stata$coefficients)), 5e-5)
    V <- component$bread/nobs(model)
    expect_lt(max(abs(V-stata$vcov)), 1e-5)
    expect_lt(max(abs(sqrt(diag(V)/diag(stata$vcov))-1)), 1e-4)
    expect_lt(abs(as.numeric(logLik(model))-stata$logLik), 1e-5)
  }
})

test_that("Poisson random-slope systems match Stata scores with validated component curvature", {
  skip_if_not_installed("glmmTMB", minimum_version = "1.1.14")
  f <- readRDS(test_path("fixtures", "glmmtmb-poisson-rs-stata-components.rds"))
  joint <- c(list(balanced = readRDS(test_path("fixtures", "glmmtmb-poisson-rs-stata-balanced-v4.rds"))),
    readRDS(test_path("fixtures", "glmmtmb-poisson-rs-stata-remaining-v5.rds")))
  outcomes <- list(balanced = c("y1", "y2"), partial_higher = c("y1_partial", "y2_partial"),
    disjoint = c("y1_left", "y2_right"))
  dat <- f$data; dat$id <- factor(dat$id)
  for (case in names(outcomes)) {
    ys <- outcomes[[case]]
    models <- lapply(ys, function(y) glmmTMB::glmmTMB(
      reformulate(c("x", "z", "(1 + x | id)"), y),
      data = dat[!is.na(dat[[y]]), ], family = poisson()))
    fit <- suest(models[[1]], models[[2]], model_names = c("Y1", "Y2"),
      observation_id = c("id", "time"), cluster = if (case == "partial_higher") "higher" else NULL)
    s <- joint[[case]]; component <- f$components[ys]
    # Joint GSEM's numerical curvature introduces small cross-equation terms.
    # Recover its score cross-product, then use independently validated Stata
    # component OIM curvature. Every reference quantity below comes from Stata;
    # raw joint V is retained unchanged in the fixture, not treated as exact.
    stata_correction <- s$n_clusters/(s$n_clusters-1)
    M <- solve(s$modelbased, s$vcov) %*% solve(s$modelbased)/stata_correction
    B <- matrix(0, 12, 12)
    B[1:6, 1:6] <- component[[1]]$vcov; B[7:12, 7:12] <- component[[2]]$vcov
    n_clusters <- if (case == "disjoint") 50 else s$n_clusters
    reference <- B %*% M %*% B*n_clusters/(n_clusters-1)
    expect_identical(names(coef(fit)), names(s$coefficients))
    expect_lt(max(abs(unname(coef(fit))-unname(unlist(lapply(component, `[[`, "coefficients"))))), 5e-5)
    expect_lt(max(abs(vcov(fit)-reference)), 1e-6)
    expect_lt(max(abs(sqrt(diag(vcov(fit))/diag(reference))-1)), 1e-4)
    expect_equal(fit$n_clusters, s$n_clusters)
    expect_equal(fit$nobs_union, 800)
    expect_equal(fit$nobs_overlap, c(balanced = 800, partial_higher = 755, disjoint = 0)[[case]])
    if (case == "disjoint") {
      expect_true(all(vcov(fit)[1:6, 7:12] == 0))
      expect_lt(max(abs(reference[1:6, 7:12])), 1e-12)
    }
  }
})
