# suest 0.1.3

* Corrected Example 6.5 to reproduce the exact published estimand: an age
  change from 20 to 30 with all other model-matrix columns held at their
  means.
* Added a supplied Stata-style `atmeans` replication helper and linked its
  output directly to the validated numerical acceptance test.
* Added explicit sample-size checks for all six paper examples so the website
  build fails instead of silently displaying results from a different sample.

* Matched the Example 6.4 sample-selection statement exactly to the
  `mecompare` Stata command example and added an explicit check and displayed
  result confirming `N = 5,062`.

* Simplified continuous-variable comparisons by using the native
  `marginaleffects` numeric forward-contrast syntax.
* Restored every example's common-sample restrictions to match the
  corresponding `mecompare` Stata example exactly, including variables used
  only to define the analysis sample.

* Revised the headings, setup descriptions, and interpretations for Examples
  6.1--6.6 to follow the language and organization of the `mecompare` Stata
  command website, adapted for the R workflow.

* Fixed package-level roxygen links and converted `NAMESPACE` and all manual
  pages to roxygen-managed files.
* Marked the package help topic as internal so it does not appear as a missing
  pkgdown reference topic.
* Added a clean-process website build script and updated the GitHub workflow
  to use the package already installed by `setup-r-dependencies`.
* Added explicit website dependencies and restored the exact validated model
  specifications in the consolidated Get Started page.
* Reorganized the pkgdown documentation around one comprehensive Get Started
  page, following the structure of the cleanplots package site.
* Moved all six paper replications and alternative-engine documentation onto
  the main page and removed the separate Articles menu.
* Enabled code evaluation during pkgdown builds so example output and results
  appear on the website, while keeping network-dependent chunks unevaluated
  during ordinary package checks.
* Condensed example code throughout the vignette, README, and function
  documentation.
* Clarified throughout that `mecompare` is a Stata command.

* Corrected the direct adapter tests to call the package's registered
  `suest_model` methods directly. The public `marginaleffects::get_predict()`
  helper normalizes `type = "response"` for native `clm` objects before S3
  dispatch and was therefore not an appropriate unit-test entry point.
* Retained the direct coefficient-perturbation check for `ordinal::clm()`
  while leaving end-to-end integration testing to `avg_comparisons()`.

# suest 0.1.2

* Fixed `ordinal::clm()` prediction and coefficient-replacement dispatch when
  model-engine vectors carry model-name attributes.
* Added direct tests that `clm()` probability predictions respond to
  coefficient perturbations.
* Corrected two focused acceptance-test reporting statements which treated
  returned error-message strings as condition objects.

# suest 0.1.1

* Made the acceptance-test runners remove stale global `suest()` and
  `suest_newdata()` functions before `devtools::load_all()`. This prevents a
  previously sourced standalone implementation from masking the package code.
* Added a test-runner guard and log entry confirming that tests use the
  package namespace and model-adapter implementation.

* Added a model-adapter layer that separates the statistical model family
  from the package or function used to fit it.
* Added tested support for `glm2::glm2()` binary logit, binary probit, and
  Poisson models.
* Added support for `ordinal::clm()` ordered logit and ordered probit models
  with flexible thresholds, proportional effects, and no scale model.
* Added defensive rejection of bias-reduced, adjusted-score, Firth, and
  penalized GLM fits such as `brglm2::brglmFit()`.
* Expanded the numerical acceptance suite to cover a full cross-model
  comparison-by-model-family matrix, including nested/mediator comparisons,
  alternative predictor operationalizations, different outcomes,
  sex-stratified samples, missing-data sample changes, and partially
  overlapping samples.
* Added same-sample and disjoint-sample tests for every supported cross-family
  comparison.
* Imported the `stats` generics used by the `suest_model` S3 methods so the
  namespace loads correctly during `R CMD check` and pkgdown builds.
* Updated GitHub Actions checkout steps to `actions/checkout@v6` and the
  GitHub Pages deployment action to version 4.8.0.
* Corrected vignette metadata so each article has a unique
  `VignetteIndexEntry`.
* Updated exact-zero covariance tests to ignore irrelevant matrix dimnames.
* Added a local package and pkgdown preflight script.

# suest 0.1.0

* Initial GitHub release.
* Combines exactly two supported cross-sectional models using a joint
  model-robust covariance matrix.
* Integrates with `marginaleffects` for predictions, comparisons, slopes,
  hypotheses, and plots.
* Supports identical, partially overlapping, and disjoint estimation samples.
* Includes replication vignettes for Examples 6.1--6.6 in Mize, Doan, and
  Long (2019).
