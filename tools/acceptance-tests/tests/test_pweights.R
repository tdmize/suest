cat("\n\nPWEIGHT LINEAR, BINARY, AND POISSON TESTS\n")

test_case("Pweights: load the Stata benchmark data", {
  benchmark_file <- file.path(
    .suest_test_state$root,
    "tools",
    "acceptance-tests",
    "data",
    "suest_pweight_benchmark.csv.gz"
  )
  pwdat <- utils::read.csv(gzfile(benchmark_file))
  pwdat$z <- factor(pwdat$z)
  expect_true(nrow(pwdat) == 2400, "Pweight benchmark should have N=2400.")
  expect_true("y_count" %in% names(pwdat), "Pweight benchmark should include y_count.")
  dim(pwdat)
})

test_case("Pweights: identical-sample lm matches Stata suest", {
  pw_base <- stats::lm(
    y_linear ~ x1 + I(x1^2) + x2 + z,
    weights = pw,
    data = pwdat
  )
  pw_adjusted <- stats::lm(
    y_linear ~ x1 + I(x1^2) + x2 + z + mediator,
    weights = pw,
    data = pwdat
  )
  pw_fit <- suest(
    pw_base,
    pw_adjusted,
    model_names = c("Base", "Adjusted"),
    weight_type = "pweight"
  )

  expected_base <- c(
    `(Intercept)` = 1.16916971854013,
    x1 = 1.08151251555599,
    `I(x1^2)` = -0.204305644358007,
    x2 = -0.596619415422722,
    z1 = 0.73940805231491
  )
  expected_adjusted <- c(
    `(Intercept)` = 1.13730398726956,
    x1 = 0.734206410164547,
    `I(x1^2)` = -0.187993075464129,
    x2 = -0.366353255805656,
    z1 = 0.534313853445243,
    mediator = 0.669675006345153
  )

  expect_near(
    stats::coef(pw_base)[names(expected_base)],
    expected_base,
    tolerance = 2e-12,
    label = "Stata pweight linear base coefficients"
  )
  expect_near(
    stats::coef(pw_adjusted)[names(expected_adjusted)],
    expected_adjusted,
    tolerance = 2e-12,
    label = "Stata pweight linear adjusted coefficients"
  )

  V <- stats::vcov(pw_fit)
  base <- pw_fit$index[[1]]
  adjusted <- pw_fit$index[[2]]
  base_names <- pw_fit$local_names[[1]]
  adjusted_names <- pw_fit$local_names[[2]]
  names(base) <- base_names
  names(adjusted) <- adjusted_names

  expect_near(
    c(
      V[base["x1"], base["x1"]],
      V[base["I(x1^2)"], base["I(x1^2)"]],
      V[base["x2"], base["x2"]],
      V[base["z1"], base["z1"]],
      V[base["(Intercept)"], base["(Intercept)"]]
    ),
    c(
      0.000908179802396165,
      0.000509531317363377,
      0.000624422782086817,
      0.00240468338485224,
      0.00154990201825923
    ),
    tolerance = 2e-12,
    label = "Stata suest pweight base covariance"
  )

  expect_near(
    c(
      V[adjusted["x1"], adjusted["x1"]],
      V[adjusted["I(x1^2)"], adjusted["I(x1^2)"]],
      V[adjusted["x2"], adjusted["x2"]],
      V[adjusted["z1"], adjusted["z1"]],
      V[adjusted["mediator"], adjusted["mediator"]],
      V[adjusted["(Intercept)"], adjusted["(Intercept)"]]
    ),
    c(
      0.000786616867552775,
      0.000367154123103818,
      0.000470817388774579,
      0.00165821472106434,
      0.000425565109255076,
      0.0010280124469469
    ),
    tolerance = 2e-12,
    label = "Stata suest pweight adjusted covariance"
  )

  expect_near(
    c(
      V[base["x1"], adjusted["x1"]],
      V[base["x1"], adjusted["mediator"]],
      V[base["(Intercept)"], adjusted["(Intercept)"]]
    ),
    c(
      0.000685594016223505,
      -0.0000301368906967427,
      0.00104734618268724
    ),
    tolerance = 2e-12,
    label = "Stata suest pweight cross-model covariance"
  )

  pw_fit
})

test_case("Pweights: partial overlap matches Stata suest", {
  partial_a <- stats::lm(
    y_linear ~ x1 + I(x1^2) + x2 + z,
    weights = pw,
    data = pwdat,
    subset = sample_a == 1
  )
  partial_b <- stats::lm(
    y_linear ~ x1 + I(x1^2) + x2 + z + mediator,
    weights = pw,
    data = pwdat,
    subset = sample_b == 1
  )
  partial_fit <- suest(
    partial_a,
    partial_b,
    model_names = c("Sample A", "Sample B"),
    weight_type = "pweight"
  )

  expect_true(
    partial_fit$nobs_overlap == 1200,
    "Partial pweight benchmark should have 1200 overlapping observations."
  )
  expect_true(
    partial_fit$nobs_union == 2400,
    "Partial pweight benchmark should have 2400 union observations."
  )

  V <- stats::vcov(partial_fit)
  first <- partial_fit$index[[1]]
  second <- partial_fit$index[[2]]
  names(first) <- partial_fit$local_names[[1]]
  names(second) <- partial_fit$local_names[[2]]

  expect_near(
    c(
      V[first["x1"], first["x1"]],
      V[second["x1"], second["x1"]],
      V[first["x1"], second["x1"]]
    ),
    c(
      0.00109981438547846,
      0.00108754243056852,
      0.000543627611255805
    ),
    tolerance = 2e-12,
    label = "Stata suest partial-overlap pweight covariance"
  )

  partial_fit
})


test_case("Pweights: identical-sample logit matches Stata suest", {
  logit_base <- suppressWarnings(stats::glm(
    y_binary ~ z + x1 + x2,
    family = stats::binomial(link = "logit"),
    weights = pw,
    data = pwdat
  ))
  logit_adjusted <- suppressWarnings(stats::glm(
    y_binary ~ z + x1 + x2 + mediator,
    family = stats::binomial(link = "logit"),
    weights = pw,
    data = pwdat
  ))
  logit_fit <- suest(
    logit_base,
    logit_adjusted,
    model_names = c("Base", "Adjusted"),
    weight_type = "pweight"
  )

  expect_near(
    stats::coef(logit_base),
    c(
      `(Intercept)` = -0.291464867252397,
      z1 = 0.610675630930217,
      x1 = 1.15478952979828,
      x2 = -0.426500374453094
    ),
    tolerance = 2e-7,
    label = "Stata pweight logit base coefficients"
  )
  expect_near(
    stats::coef(logit_adjusted),
    c(
      `(Intercept)` = -0.321675899521929,
      z1 = 0.504730841372048,
      x1 = 0.968760597401695,
      x2 = -0.285123420455741,
      mediator = 0.446442155427318
    ),
    tolerance = 2e-7,
    label = "Stata pweight logit adjusted coefficients"
  )

  V <- stats::vcov(logit_fit)
  first <- logit_fit$index[[1]]
  second <- logit_fit$index[[2]]
  names(first) <- logit_fit$local_names[[1]]
  names(second) <- logit_fit$local_names[[2]]

  expect_near(
    c(
      V[first["z1"], first["z1"]],
      V[first["x1"], first["x1"]],
      V[second["z1"], second["z1"]],
      V[second["x1"], second["x1"]],
      V[first["z1"], second["z1"]],
      V[first["x1"], second["x1"]],
      V[first["(Intercept)"], second["(Intercept)"]]
    ),
    c(
      0.0100789023522823,
      0.0038954110827882,
      0.01057275657174,
      0.00442270821969415,
      0.0101039758105563,
      0.00394534160074276,
      0.00532287607780682
    ),
    tolerance = 2e-8,
    label = "Stata suest pweight logit covariance"
  )

  logit_fit
})

test_case("Pweights: identical-sample probit matches Stata suest", {
  probit_base <- suppressWarnings(stats::glm(
    y_binary ~ z + x1 + x2,
    family = stats::binomial(link = "probit"),
    weights = pw,
    data = pwdat
  ))
  probit_adjusted <- suppressWarnings(stats::glm(
    y_binary ~ z + x1 + x2 + mediator,
    family = stats::binomial(link = "probit"),
    weights = pw,
    data = pwdat
  ))
  probit_fit <- suest(
    probit_base,
    probit_adjusted,
    model_names = c("Base", "Adjusted"),
    weight_type = "pweight"
  )

  expect_near(
    stats::coef(probit_base),
    c(
      `(Intercept)` = -0.173515606591847,
      z1 = 0.365133534393993,
      x1 = 0.688795935465237,
      x2 = -0.256385233932015
    ),
    tolerance = 2e-7,
    label = "Stata pweight probit base coefficients"
  )
  expect_near(
    stats::coef(probit_adjusted),
    c(
      `(Intercept)` = -0.190789270790028,
      z1 = 0.300753316313436,
      x1 = 0.577071454230468,
      x2 = -0.170400127258527,
      mediator = 0.266188627806507
    ),
    tolerance = 2e-7,
    label = "Stata pweight probit adjusted coefficients"
  )

  V <- stats::vcov(probit_fit)
  first <- probit_fit$index[[1]]
  second <- probit_fit$index[[2]]
  names(first) <- probit_fit$local_names[[1]]
  names(second) <- probit_fit$local_names[[2]]

  expect_near(
    c(
      V[first["z1"], first["z1"]],
      V[first["x1"], first["x1"]],
      V[second["z1"], second["z1"]],
      V[second["x1"], second["x1"]],
      V[first["z1"], second["z1"]],
      V[first["x1"], second["x1"]],
      V[first["(Intercept)"], second["(Intercept)"]]
    ),
    c(
      0.00357837816658024,
      0.0012263014865656,
      0.00371516410079853,
      0.00144110822686034,
      0.00356577542046508,
      0.0012528275318619,
      0.00185327262801978
    ),
    tolerance = 2e-8,
    label = "Stata suest pweight probit covariance"
  )

  probit_fit
})

test_case("Pweights: logit-probit pair matches Stata suest", {
  mixed_logit <- suppressWarnings(stats::glm(
    y_binary ~ z + x1 + x2 + mediator,
    family = stats::binomial(link = "logit"),
    weights = pw,
    data = pwdat
  ))
  mixed_probit <- suppressWarnings(stats::glm(
    y_binary ~ z + x1 + x2 + mediator,
    family = stats::binomial(link = "probit"),
    weights = pw,
    data = pwdat
  ))
  mixed_fit <- suest(
    mixed_logit,
    mixed_probit,
    model_names = c("Logit", "Probit"),
    weight_type = "pweight"
  )
  V <- stats::vcov(mixed_fit)
  first <- mixed_fit$index[[1]]
  second <- mixed_fit$index[[2]]
  names(first) <- mixed_fit$local_names[[1]]
  names(second) <- mixed_fit$local_names[[2]]

  expect_near(
    c(
      V[first["z1"], second["z1"]],
      V[first["x1"], second["x1"]],
      V[first["x2"], second["x2"]],
      V[first["mediator"], second["mediator"]],
      V[first["(Intercept)"], second["(Intercept)"]]
    ),
    c(
      0.00625228004734594,
      0.00251463935790408,
      0.00183464594417685,
      0.00175372869021398,
      0.00325361724835018
    ),
    tolerance = 2e-8,
    label = "Stata suest pweight logit-probit cross covariance"
  )

  mixed_fit
})

test_case("Pweights: partial logit and disjoint probit match Stata suest", {
  partial_a <- suppressWarnings(stats::glm(
    y_binary ~ z + x1 + x2,
    family = stats::binomial(),
    weights = pw,
    data = pwdat,
    subset = sample_a == 1
  ))
  partial_b <- suppressWarnings(stats::glm(
    y_binary ~ z + x1 + x2 + mediator,
    family = stats::binomial(),
    weights = pw,
    data = pwdat,
    subset = sample_b == 1
  ))
  partial_fit <- suest(
    partial_a,
    partial_b,
    model_names = c("Sample A", "Sample B"),
    weight_type = "pweight"
  )
  V <- stats::vcov(partial_fit)
  first <- partial_fit$index[[1]]
  second <- partial_fit$index[[2]]
  names(first) <- partial_fit$local_names[[1]]
  names(second) <- partial_fit$local_names[[2]]

  expect_near(
    c(
      V[first["z1"], first["z1"]],
      V[first["x1"], first["x1"]],
      V[second["z1"], second["z1"]],
      V[second["x1"], second["x1"]],
      V[first["z1"], second["z1"]],
      V[first["x1"], second["x1"]]
    ),
    c(
      0.0133500487082458,
      0.00502717903187792,
      0.0141740679744691,
      0.00569742344268775,
      0.00932733893068713,
      0.00339263853251683
    ),
    tolerance = 5e-8,
    label = "Stata suest partial-overlap pweight logit covariance"
  )

  disjoint_left <- suppressWarnings(stats::glm(
    y_binary ~ z + x1 + x2,
    family = stats::binomial(link = "probit"),
    weights = pw,
    data = pwdat,
    subset = sample_left == 1
  ))
  disjoint_right <- suppressWarnings(stats::glm(
    y_binary ~ z + x1 + x2 + mediator,
    family = stats::binomial(link = "probit"),
    weights = pw,
    data = pwdat,
    subset = sample_right == 1
  ))
  disjoint_fit <- suest(
    disjoint_left,
    disjoint_right,
    model_names = c("Left", "Right"),
    weight_type = "pweight"
  )
  V <- stats::vcov(disjoint_fit)
  first <- disjoint_fit$index[[1]]
  second <- disjoint_fit$index[[2]]
  names(first) <- disjoint_fit$local_names[[1]]
  names(second) <- disjoint_fit$local_names[[2]]

  expect_near(
    c(
      V[first["z1"], first["z1"]],
      V[first["x1"], first["x1"]],
      V[second["z1"], second["z1"]],
      V[second["x1"], second["x1"]]
    ),
    c(
      0.00719359651596096,
      0.00250224819667147,
      0.00750786567754631,
      0.00292495034624621
    ),
    tolerance = 5e-8,
    label = "Stata suest disjoint pweight probit covariance"
  )
  expect_true(
    all(V[first, second] == 0),
    "Disjoint pweighted probit models should have zero cross covariance."
  )

  list(partial_logit = partial_fit, disjoint_probit = disjoint_fit)
})

test_case("Pweights: linear-binary pairs match Stata suest", {
  linear <- stats::lm(
    y_linear ~ x1 + x2 + z,
    weights = pw,
    data = pwdat
  )
  logit <- suppressWarnings(stats::glm(
    y_binary ~ z + x1 + x2,
    family = stats::binomial(link = "logit"),
    weights = pw,
    data = pwdat
  ))
  probit <- suppressWarnings(stats::glm(
    y_binary ~ z + x1 + x2,
    family = stats::binomial(link = "probit"),
    weights = pw,
    data = pwdat
  ))

  logit_fit <- suest(
    linear,
    logit,
    model_names = c("Linear", "Logit"),
    weight_type = "pweight"
  )
  probit_fit <- suest(
    linear,
    probit,
    model_names = c("Linear", "Probit"),
    weight_type = "pweight"
  )

  extract_cross <- function(object) {
    V <- stats::vcov(object)
    first <- object$index[[1]]
    second <- object$index[[2]]
    names(first) <- object$local_names[[1]]
    names(second) <- object$local_names[[2]]
    c(
      V[first["z1"], second["z1"]],
      V[first["x1"], second["x1"]],
      V[first["x2"], second["x2"]],
      V[first["(Intercept)"], second["(Intercept)"]]
    )
  }

  expect_near(
    extract_cross(logit_fit),
    c(
      0.000773590563405227,
      0.00040340147395683,
      0.000294008560643943,
      0.000373285208307171
    ),
    tolerance = 2e-8,
    label = "Stata suest pweight linear-logit cross covariance"
  )
  expect_near(
    extract_cross(probit_fit),
    c(
      0.000480962742532784,
      0.000228765816666023,
      0.000188291778586675,
      0.000230416553599996
    ),
    tolerance = 2e-8,
    label = "Stata suest pweight linear-probit cross covariance"
  )

  list(linear_logit = logit_fit, linear_probit = probit_fit)
})

test_case("Pweights: identical-sample Poisson matches Stata suest", {
  poisson_base <- stats::glm(
    y_count ~ z + x1 + x2,
    family = stats::poisson(link = "log"),
    weights = pw,
    data = pwdat
  )
  poisson_adjusted <- stats::glm(
    y_count ~ z + x1 + x2 + mediator,
    family = stats::poisson(link = "log"),
    weights = pw,
    data = pwdat
  )
  poisson_fit <- suest(
    poisson_base,
    poisson_adjusted,
    model_names = c("Base", "Adjusted"),
    weight_type = "pweight"
  )

  expect_near(
    stats::coef(poisson_base),
    c(
      `(Intercept)` = 0.377833110250777,
      z1 = 0.399376687509099,
      x1 = 0.454534363465156,
      x2 = -0.259712042730786
    ),
    tolerance = 2e-8,
    label = "Stata pweight Poisson base coefficients"
  )
  expect_near(
    stats::coef(poisson_adjusted),
    c(
      `(Intercept)` = 0.341399960747014,
      z1 = 0.326438164018034,
      x1 = 0.334165570936692,
      x2 = -0.17898026429262,
      mediator = 0.244346960888699
    ),
    tolerance = 2e-8,
    label = "Stata pweight Poisson adjusted coefficients"
  )

  V <- stats::vcov(poisson_fit)
  first <- poisson_fit$index[[1]]
  second <- poisson_fit$index[[2]]
  names(first) <- poisson_fit$local_names[[1]]
  names(second) <- poisson_fit$local_names[[2]]

  expect_near(
    c(
      V[first["z1"], first["z1"]],
      V[first["x1"], first["x1"]],
      V[second["z1"], second["z1"]],
      V[second["x1"], second["x1"]],
      V[first["z1"], second["z1"]],
      V[first["x1"], second["x1"]],
      V[first["(Intercept)"], second["(Intercept)"]]
    ),
    c(
      0.0011274798055861,
      0.00026309329966446,
      0.0010505496238107,
      0.00030621991508818,
      0.00101250836849846,
      0.000237812124222109,
      0.000679336123875914
    ),
    tolerance = 2e-9,
    label = "Stata suest pweight Poisson covariance"
  )

  poisson_fit
})


test_case("Pweights: partial and disjoint Poisson match Stata suest", {
  partial_a <- stats::glm(
    y_count ~ z + x1 + x2,
    family = stats::poisson(),
    weights = pw,
    data = pwdat,
    subset = sample_a == 1
  )
  partial_b <- stats::glm(
    y_count ~ z + x1 + x2 + mediator,
    family = stats::poisson(),
    weights = pw,
    data = pwdat,
    subset = sample_b == 1
  )
  partial_fit <- suest(
    partial_a,
    partial_b,
    model_names = c("Sample A", "Sample B"),
    weight_type = "pweight"
  )
  V <- stats::vcov(partial_fit)
  first <- partial_fit$index[[1]]
  second <- partial_fit$index[[2]]
  names(first) <- partial_fit$local_names[[1]]
  names(second) <- partial_fit$local_names[[2]]

  expect_near(
    c(
      V[first["z1"], first["z1"]],
      V[first["x1"], first["x1"]],
      V[second["z1"], second["z1"]],
      V[second["x1"], second["x1"]],
      V[first["z1"], second["z1"]],
      V[first["x1"], second["x1"]]
    ),
    c(
      0.00144741477763976,
      0.000348425584142248,
      0.00143362495092192,
      0.000382413087492499,
      0.000902664263440914,
      0.000225455937067902
    ),
    tolerance = 2e-9,
    label = "Stata suest partial-overlap pweight Poisson covariance"
  )

  disjoint_left <- stats::glm(
    y_count ~ z + x1 + x2,
    family = stats::poisson(),
    weights = pw,
    data = pwdat,
    subset = sample_left == 1
  )
  disjoint_right <- stats::glm(
    y_count ~ z + x1 + x2 + mediator,
    family = stats::poisson(),
    weights = pw,
    data = pwdat,
    subset = sample_right == 1
  )
  disjoint_fit <- suest(
    disjoint_left,
    disjoint_right,
    model_names = c("Left", "Right"),
    weight_type = "pweight"
  )
  V <- stats::vcov(disjoint_fit)
  first <- disjoint_fit$index[[1]]
  second <- disjoint_fit$index[[2]]
  names(first) <- disjoint_fit$local_names[[1]]
  names(second) <- disjoint_fit$local_names[[2]]

  expect_near(
    c(
      V[first["z1"], first["z1"]],
      V[first["x1"], first["x1"]],
      V[second["z1"], second["z1"]],
      V[second["x1"], second["x1"]]
    ),
    c(
      0.00207804946963707,
      0.00051913342690592,
      0.00229213997404226,
      0.000527406314395024
    ),
    tolerance = 2e-9,
    label = "Stata suest disjoint pweight Poisson covariance"
  )
  expect_true(
    all(V[first, second] == 0),
    "Disjoint pweighted Poisson models should have zero cross covariance."
  )

  list(partial_poisson = partial_fit, disjoint_poisson = disjoint_fit)
})


test_case("Pweights: coherent extensions beyond Stata restrictions", {
  off_overlap_a <- stats::lm(
    y_linear ~ x1 + x2,
    weights = pw,
    data = pwdat,
    subset = sample_a == 1
  )
  off_overlap_b <- stats::lm(
    y_linear ~ x1 + x2 + mediator,
    weights = pw_partial_b,
    data = pwdat,
    subset = sample_b == 1
  )
  off_overlap_fit <- suest(
    off_overlap_a,
    off_overlap_b,
    weight_type = "pweight"
  )

  disjoint_left <- stats::lm(
    y_linear ~ x1 + x2,
    weights = pw_left,
    data = pwdat,
    subset = sample_left == 1
  )
  disjoint_right <- stats::lm(
    y_linear ~ x1 + x2 + mediator,
    weights = pw_right,
    data = pwdat,
    subset = sample_right == 1
  )
  disjoint_fit <- suest(
    disjoint_left,
    disjoint_right,
    weight_type = "pweight"
  )

  expect_true(
    off_overlap_fit$nobs_overlap == 1200,
    "Different off-overlap weights should retain the correct overlap."
  )
  expect_true(
    disjoint_fit$nobs_overlap == 0,
    "Different disjoint weights should have zero overlap."
  )
  expect_true(
    all(offdiag_vcov(disjoint_fit) == 0),
    "Disjoint pweighted models should have zero cross-model covariance."
  )

  disjoint_common_left <- stats::lm(
    y_linear ~ x1 + I(x1^2) + x2 + z,
    weights = pw,
    data = pwdat,
    subset = sample_left == 1
  )
  disjoint_common_right <- stats::lm(
    y_linear ~ x1 + I(x1^2) + x2 + z + mediator,
    weights = pw,
    data = pwdat,
    subset = sample_right == 1
  )
  disjoint_common_fit <- suest(
    disjoint_common_left,
    disjoint_common_right,
    weight_type = "pweight"
  )
  V <- stats::vcov(disjoint_common_fit)
  first <- disjoint_common_fit$index[[1]]
  second <- disjoint_common_fit$index[[2]]
  names(first) <- disjoint_common_fit$local_names[[1]]
  names(second) <- disjoint_common_fit$local_names[[2]]

  expect_near(
    c(
      V[first["x1"], first["x1"]],
      V[second["x1"], second["x1"]]
    ),
    c(
      0.0016114699738351,
      0.00174364808376926
    ),
    tolerance = 2e-12,
    label = "Stata suest disjoint-sample pweight covariance"
  )

  binary_off_a <- suppressWarnings(stats::glm(
    y_binary ~ x1 + x2,
    family = stats::binomial(),
    weights = pw,
    data = pwdat,
    subset = sample_a == 1
  ))
  binary_off_b <- suppressWarnings(stats::glm(
    y_binary ~ x1 + x2 + mediator,
    family = stats::binomial(),
    weights = pw_partial_b,
    data = pwdat,
    subset = sample_b == 1
  ))
  binary_off_fit <- suest(
    binary_off_a,
    binary_off_b,
    weight_type = "pweight"
  )

  binary_disjoint_left <- suppressWarnings(stats::glm(
    y_binary ~ x1 + x2,
    family = stats::binomial(link = "probit"),
    weights = pw_left,
    data = pwdat,
    subset = sample_left == 1
  ))
  binary_disjoint_right <- suppressWarnings(stats::glm(
    y_binary ~ x1 + x2 + mediator,
    family = stats::binomial(link = "probit"),
    weights = pw_right,
    data = pwdat,
    subset = sample_right == 1
  ))
  binary_disjoint_fit <- suest(
    binary_disjoint_left,
    binary_disjoint_right,
    weight_type = "pweight"
  )

  expect_true(
    binary_off_fit$nobs_overlap == 1200,
    "Binary weights may differ outside the overlapping observations."
  )
  expect_true(
    all(offdiag_vcov(binary_disjoint_fit) == 0),
    "Disjoint binary models with different weights should have zero cross covariance."
  )

  poisson_off_a <- stats::glm(
    y_count ~ x1 + x2,
    family = stats::poisson(),
    weights = pw,
    data = pwdat,
    subset = sample_a == 1
  )
  poisson_off_b <- stats::glm(
    y_count ~ x1 + x2 + mediator,
    family = stats::poisson(),
    weights = pw_partial_b,
    data = pwdat,
    subset = sample_b == 1
  )
  poisson_off_fit <- suest(
    poisson_off_a,
    poisson_off_b,
    weight_type = "pweight"
  )

  poisson_disjoint_left <- stats::glm(
    y_count ~ x1 + x2,
    family = stats::poisson(),
    weights = pw_left,
    data = pwdat,
    subset = sample_left == 1
  )
  poisson_disjoint_right <- stats::glm(
    y_count ~ x1 + x2 + mediator,
    family = stats::poisson(),
    weights = pw_right,
    data = pwdat,
    subset = sample_right == 1
  )
  poisson_disjoint_fit <- suest(
    poisson_disjoint_left,
    poisson_disjoint_right,
    weight_type = "pweight"
  )

  expect_true(
    poisson_off_fit$nobs_overlap == 1200,
    "Poisson weights may differ outside the overlapping observations."
  )
  expect_true(
    all(offdiag_vcov(poisson_disjoint_fit) == 0),
    "Disjoint Poisson models with different weights should have zero cross covariance."
  )

  list(
    partial_different_off_overlap = off_overlap_fit,
    disjoint_different_weights = disjoint_fit,
    binary_partial_different_off_overlap = binary_off_fit,
    binary_disjoint_different_weights = binary_disjoint_fit,
    poisson_partial_different_off_overlap = poisson_off_fit,
    poisson_disjoint_different_weights = poisson_disjoint_fit
  )
})

test_case("Pweights: invalid overlap and undeclared weights are rejected", {
  undeclared1 <- stats::lm(y_linear ~ x1, weights = pw, data = pwdat)
  undeclared2 <- stats::lm(
    y_linear ~ x1 + x2,
    weights = pw,
    data = pwdat
  )
  expect_error(
    suest(undeclared1, undeclared2),
    'weight_type = "pweight"'
  )

  scaled <- stats::lm(
    y_linear ~ x1 + x2,
    weights = pw10,
    data = pwdat
  )
  expect_error(
    suest(
      undeclared1,
      scaled,
      weight_type = "pweight"
    ),
    "must agree for observations included in both models"
  )

  TRUE
})
