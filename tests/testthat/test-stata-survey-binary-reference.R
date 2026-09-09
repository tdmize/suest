# Returned Stata 19.5 survey-binary gate.
# Logit covariance is numerically comparable to survey::svyglm().
# The returned Stata probit covariance is reproduced by an observed-information
# bread, unlike survey::svyglm() (expected/Fisher information), so only coefficients are
# cross-language comparable for probit; native R covariance blocks remain the
# governing contract.
test_that("survey binary systems match returned Stata gate within the native-R contract", {
  skip_if_not_installed("survey")
  dat <- read.csv(test_path("fixtures", "survey-binary-stata-data.csv"))
  dat$f <- factor(dat$f)
  ref <- read.csv(test_path("fixtures", "survey-binary-stata-reference.csv"))
  old <- options(survey.lonely.psu = "fail",
    survey.adjust.domain.lonely = FALSE)
  on.exit(options(old), add = TRUE)

  get_b <- function(case) {
    z <- ref[ref$case == case & ref$kind == "b", ]
    z$value[order(z$col)]
  }
  get_V <- function(case) {
    z <- ref[ref$case == case & ref$kind == "V", ]
    n <- max(z$row)
    out <- matrix(NA_real_, n, n)
    out[cbind(z$row, z$col)] <- z$value
    out
  }

  for (case in seq_len(6L)) {
    link <- if (case %in% c(1L, 3L, 5L)) "logit" else "probit"
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

    if (case %in% c(3L, 5L)) {
      fa <- yb ~ x + f
      fb <- yb2 ~ x + z + f
    } else {
      fa <- yb ~ x
      fb <- yb2 ~ x + z
    }
    if (case %in% c(4L, 5L, 6L)) {
      fa <- update(fa, . ~ . + offset(off_a))
      fb <- update(fb, . ~ . + offset(off_b))
    }

    a <- survey::svyglm(fa, subset(design, a_keep),
      family = quasibinomial(link), influence = TRUE,
      control = list(epsilon = 1e-12, maxit = 100))
    b <- survey::svyglm(fb, subset(design, b_keep),
      family = quasibinomial(link), influence = TRUE,
      control = list(epsilon = 1e-12, maxit = 100))
    fit <- suest(a, b, model_names = c("A", "B"),
      survey_design = design, observation_id = "id")

    expect_lt(max(abs(unname(coef(fit)) - get_b(case))),
      if (link == "logit") 5e-8 else 2e-6)
    expect_equal(unname(vcov(fit)[fit$index[[1]], fit$index[[1]]]),
      unname(vcov(a)), tolerance = 2e-10)
    expect_equal(unname(vcov(fit)[fit$index[[2]], fit$index[[2]]]),
      unname(vcov(b)), tolerance = 2e-10)
    if (link == "logit")
      expect_lt(max(abs(unname(vcov(fit)) - get_V(case))), 5e-7)
  }

  rc <- ref$value[ref$kind == "suest2_rc"]
  expect_identical(as.integer(rc), c(0L, 0L, 322L, 322L, 322L, 322L))
})
