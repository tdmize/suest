stata_cluster_reference_data <- function(n = 1200L) {
  set.seed(6109)
  id <- seq_len(n)
  cluster4 <- ceiling(id / 4)
  cluster2 <- ceiling(id / 2)
  side <- id %% 2L
  cluster_effect <- rnorm(max(cluster4))[cluster4]
  x <- rnorm(n) + 0.25*cluster_effect
  z <- rnorm(n)
  z1 <- rnorm(n)
  z2 <- rnorm(n)
  disturbance <- rnorm(n)
  endogenous <- 0.75*z1 + 0.20*z2 + 0.35*disturbance + rnorm(n)
  y <- 1 + 0.45*x - 0.30*z + 0.60*cluster_effect + rnorm(n)
  y_iv <- 0.5 + 0.6*x + 1.1*endogenous + disturbance
  probability <- plogis(-0.25 + 0.50*x - 0.20*z + 0.45*cluster_effect)
  b <- rbinom(n, 1, probability)
  pw <- exp(0.20*z + 0.10*cluster_effect)

  data.frame(
    id, cluster4, cluster2, side, x, z, z1, z2, endogenous, y, y_iv, b, pw
  )
}

stata_cluster_order <- c(
  "LM::x", "LM::z", "LM::(Intercept)", "LM::lnvar",
  "Logit::x", "Logit::z", "Logit::(Intercept)"
)

stata_cluster_selected_vcov <- function(object) {
  V <- vcov(object)[stata_cluster_order, stata_cluster_order]
  c(V[1, 1], V[4, 4], V[5, 5], V[5, 1], V[5, 4])
}

expect_stata_reference <- function(actual, expected, tolerance) {
  expect_lt(max(abs(unname(actual) - unname(expected))), tolerance)
}

test_that("clustered lm-logit systems match Stata suest2 references", {
  dat <- stata_cluster_reference_data()

  identical_fit <- suest(
    lm(y ~ x + z, data = dat),
    glm(b ~ x + z, family = binomial(), data = dat),
    model_names = c("LM", "Logit"), cluster = "cluster4"
  )
  partial_fit <- suest(
    lm(y ~ x + z, data = dat, subset = id <= 900),
    glm(b ~ x + z, family = binomial(), data = dat, subset = id >= 301),
    model_names = c("LM", "Logit"), cluster = "cluster4"
  )
  disjoint_fit <- suest(
    lm(y ~ x + z, data = dat, subset = side == 0),
    glm(b ~ x + z, family = binomial(), data = dat, subset = side == 1),
    model_names = c("LM", "Logit"), cluster = "cluster2"
  )

  expect_stata_reference(
    coef(identical_fit)[stata_cluster_order],
    c(
      .631226331862798, -.339700060553379, .978243363734057,
      .310913460362295, .544528379126474, -.258508911766364,
      -.227017544979117
    ),
    3e-8
  )
  expect_stata_reference(
    stata_cluster_selected_vcov(identical_fit),
    c(
      .00122412589307078, .00169166167745193, .00338953672856865,
      .000514668379858386, -.000112381666637016
    ),
    2e-9
  )

  expect_stata_reference(
    coef(partial_fit)[stata_cluster_order],
    c(
      .652186360837345, -.349840539594271, .920742876848211,
      .312572766124102, .545174694175398, -.278365205818505,
      -.203990222017124
    ),
    3e-8
  )
  expect_stata_reference(
    stata_cluster_selected_vcov(partial_fit),
    c(
      .00150259022803089, .00239237619703571, .00424639864208483,
      .000365887010513948, .000238414858073437
    ),
    2e-9
  )

  expect_stata_reference(
    coef(disjoint_fit)[stata_cluster_order],
    c(
      .61613773223269, -.379203594323319, .971746605788846,
      .389388245365534, .467502458992567, -.259766884430857,
      -.155029473931919
    ),
    3e-8
  )
  expect_stata_reference(
    stata_cluster_selected_vcov(disjoint_fit),
    c(
      .00256111266609008, .00337744118473561, .00709796874434713,
      .000288137503705597, .000123832027348108
    ),
    2e-9
  )
})

test_that("clustered 2SLS matches the Stata suest2 raw-influence route", {
  skip_if_not_installed("fixest")
  dat <- stata_cluster_reference_data()
  base <- fixest::feols(y_iv ~ x | endogenous ~ z1, data = dat)
  adjusted <- fixest::feols(y_iv ~ x | endogenous ~ z1 + z2, data = dat)
  fit <- suest(
    base, adjusted, model_names = c("Base", "Adjusted"), cluster = "cluster4"
  )
  order <- c(
    "Base::fit_endogenous", "Base::x", "Base::(Intercept)",
    "Adjusted::fit_endogenous", "Adjusted::x", "Adjusted::(Intercept)"
  )
  V <- vcov(fit)[order, order]

  expect_stata_reference(
    coef(fit)[order],
    c(
      1.10443705314135, .628750235611305, .473666970462023,
      1.09793790847518, .628560118323013, .473842333313829
    ),
    6e-9
  )
  expect_stata_reference(
    c(V[1, 1], V[2, 2], V[4, 4], V[4, 1], V[5, 2]),
    c(
      .00101280586308529, .00069964815987473, .000986047406942866,
      .000977279081154815, .000700851297614898
    ),
    2e-11
  )
})

test_that("clustered pweighted lm-logit matches Stata including lnvar", {
  dat <- stata_cluster_reference_data()
  linear <- lm(y ~ x + z, weights = pw, data = dat)
  logit <- suppressWarnings(glm(
    b ~ x + z, family = binomial(), weights = pw, data = dat
  ))
  fit <- suest(
    linear, logit, model_names = c("LM", "Logit"),
    weight_type = "pweight", cluster = "cluster4"
  )

  expect_stata_reference(
    coef(fit)[stata_cluster_order],
    c(
      .630059198025592, -.333772639845326, 1.03045431639061,
      .31289210710126, .519877151687213, -.25324921044131,
      -.183776625501182
    ),
    2e-8
  )
  expect_stata_reference(
    stata_cluster_selected_vcov(fit),
    c(
      .00124157710090435, .00174938389455396, .00335029566961628,
      .000499593080480315, -.0000948762755362078
    ),
    2e-10
  )
})
