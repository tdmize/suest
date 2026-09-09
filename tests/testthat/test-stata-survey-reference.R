# Full matrices from the user-returned Stata 19.5 revision-2 survey gate.
test_that("survey linear systems match returned Stata joint matrices", {
  skip_if_not_installed("survey")
  fixture <- readRDS(test_path("fixtures", "survey-stata.rds"))
  dat <- fixture$data
  old <- options(survey.lonely.psu = "fail",
    survey.adjust.domain.lonely = FALSE)
  on.exit(options(old), add = TRUE)

  for (case in seq_len(6L)) {
    design <- if (case == 1L) {
      survey::svydesign(~psu, strata = ~strata, weights = ~w, data = dat)
    } else {
      survey::svydesign(~psu, strata = ~strata, weights = ~w,
        fpc = ~pop, data = dat)
    }
    a_keep <- switch(case, rep(TRUE, nrow(dat)), rep(TRUE, nrow(dat)),
      dat$unit <= 7, dat$unit <= 5, dat$psu %% 6 == 1, dat$strata <= 2)
    b_keep <- switch(case, rep(TRUE, nrow(dat)), rep(TRUE, nrow(dat)),
      dat$unit >= 4, dat$unit > 5, dat$psu %% 6 == 2, dat$strata > 2)
    first <- survey::svyglm(y ~ x, subset(design, a_keep))
    second <- survey::svyglm(y2 ~ x + z, subset(design, b_keep))
    fit <- suest(first, second, model_names = c("A", "B"),
      survey_design = design, observation_id = "id")
    ref <- fixture$reference[[case]]

    expect_lt(max(abs(coef(fit) - ref$b)), 1e-12)
    expect_lt(max(abs(vcov(fit) - ref$V)), 1e-13)
    expect_equal(unname(vcov(fit)[fit$index[[1]], fit$index[[1]]]),
      ref$native_a, tolerance = 1e-13)
    expect_equal(unname(vcov(fit)[fit$index[[2]], fit$index[[2]]]),
      ref$native_b, tolerance = 1e-13)
  }
})

test_that("accepted suest2 survey systems equal the stacked Stata reference", {
  fixture <- readRDS(test_path("fixtures", "survey-stata.rds"))
  for (case in 1:2) {
    ref <- fixture$reference[[case]]
    expect_equal(ref$suest2_b, ref$b, tolerance = 1e-14)
    expect_equal(ref$suest2_V, ref$V, tolerance = 1e-14)
    expect_identical(ref$suest2_rc, 0L)
  }
  expect_identical(vapply(fixture$reference[3:6], `[[`, integer(1),
    "suest2_rc"), rep(322L, 4L))
})
