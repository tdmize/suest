# Development notes

Version 0.1.4 is the certified survey-binary release candidate. Its
expanded model adapters, pweights, linear ancillary parameter, sample
alignment, cluster-robust covariance, linear-panel, random-intercept
Gaussian ML, bivariate-probit, IV-probit, GEE, multiple-imputation, and
narrow one-stage survey routes are covered by focused unit tests. The
established families are also covered by the full numerical acceptance
suite. Cross-language benchmarks live in `tools/stata-benchmarks/`.

Post-0.1.4 development prioritizes nonlinear panel families first and
general GLMM/multilevel models second. See
`tools/ROADMAP-POST-0.1.4-20260909.md` for the scoped sequence.

## First local setup

Open `suest.Rproj` and install the development packages:

``` r

install.packages(c(
  "brglm2", "callr", "devtools", "glm2", "ordinal", "pkgdown",
  "rlang", "roxygen2", "testthat"
))
```

Restart R after updating packages.

## Build and preview the website

From the package directory:

``` r

source("tools/build-site.R")
```

Then open `docs/articles/suest.html`.

The script documents the package and builds the site in separate clean R
processes. This avoids failures caused by the package or an older
dependency remaining loaded in the interactive session.

## Run all local checks

``` r

source("tools/preflight.R")
```

## Numerical acceptance tests

``` r

source("tools/acceptance-tests/run_pweight_tests.R")
source("tools/acceptance-tests/run_model_adapter_tests.R")
source("tools/acceptance-tests/run_all_tests.R")
```

The included GitHub Actions workflows run `R CMD check` on Windows,
macOS, and Linux, build and deploy the pkgdown site, and run the full
Mize, Doan, and Long (2019) numerical acceptance suite.
