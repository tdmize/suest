# Development notes

## First local setup

Open `suest.Rproj` and install the development packages:

```r
install.packages(c(
  "brglm2", "callr", "devtools", "glm2", "ordinal", "pkgdown",
  "rlang", "roxygen2", "testthat"
))
```

Restart R after updating packages.

## Build and preview the website

From the package directory:

```r
source("tools/build-site.R")
```

Then open `docs/articles/suest.html`.

The script documents the package and builds the site in separate clean R
processes. This avoids failures caused by the package or an older dependency
remaining loaded in the interactive session.

## Run all local checks

```r
source("tools/preflight.R")
```

## Numerical acceptance tests

```r
source("tools/acceptance-tests/run_model_adapter_tests.R")
source("tools/acceptance-tests/run_all_tests.R")
```

The included GitHub Actions workflows run `R CMD check` on Windows, macOS, and
Linux, build and deploy the pkgdown site, and run the full Mize, Doan, and Long
(2019) numerical acceptance suite.
