# Development notes

Version 0.1.6 adds NB2 random-intercept support and correlated
random-intercept/ numeric-slope support for Poisson, NB2, and
Bernoulli-logit glmmTMB models. Independent validation covers balanced,
partial-overlap, and disjoint samples. Raw GSEM covariance differences
remain documented; strict raw-Stata parity is not claimed. Existing
survey, MI, panel, and single-level contracts are unchanged.

See the [release
notes](https://github.com/tdmize/suest/blob/main/tools/RELEASE-0.1.6.md)
and [current
handoff](https://github.com/tdmize/suest/blob/main/tools/SESSION-HANDOFF-0.1.6-20260925.md)
for scope, validation, and publication status. Cross-language benchmarks
live in `tools/stata-benchmarks/`. The post-0.1.4 roadmap and earlier
handoffs are historical records; choose subsequent model additions
explicitly.

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
