# suest 0.1.4

* Extended the narrow one-stage survey route to binary-response
  `survey::svyglm()` models fitted with `quasibinomial("logit")` or
  `quasibinomial("probit")`. Observation-level coefficient influences are
  reconstructed from the retained final GLM working residuals, working weights,
  and native expected-information bread, zero-padded to the
  full common design, and required to reproduce each native `svyglm()`
  covariance block before the joint covariance is formed. Focused coverage
  includes identical, overlapping, and disjoint domains, factors, formula
  offsets, and first-stage FPCs. Returned Stata 19.5 systems confirm logit
  coefficient/covariance parity. Probit coefficients agree, but the returned
  Stata covariance is reproduced by an observed-information bread whereas
  `svyglm()` uses expected/Fisher information; the R route deliberately
  preserves its native covariance convention rather than forcing cross-engine
  finite-sample VCE equality.
* Certified the survey binary-response increment under R 4.3.3 and `survey`
  4.5. The original certification used `marginaleffects` 0.31.0; release
  preflight also passed under `marginaleffects` 1.0.0 after correcting one test
  to ignore a new metadata attribute on a numerically unchanged standard
  error. The gates comprise 132 focused expectations, 147 complete-unit tests
  with 872 recorded outcomes, and 144 numerical acceptance cases; the source
  package builds with its vignette, installs and loads from a clean library,
  and completes `R CMD check --no-manual` with no errors, warnings, or notes.
* Added initial survey linear-model support through `survey_design =` and
  explicit observation-ID columns. Joint coefficient covariance uses one-stage
  PSU/stratum linearization, including FPCs and model-specific domains, and is
  checked against each native `svyglm()` covariance. This route excludes
  survey `lnvar`, nonlinear survey families, and complex design extensions.
* Survey predictions and effects use the joint coefficient covariance;
  `suest_newdata()` supplies sampling weights for explicit weighted averaging.
  Default inference remains asymptotic normal. The first returned Stata survey
  benchmark validates all 12 native model blocks. Six stacked Stata survey
  systems match the full R joint covariance to machine precision; the two
  same-sample systems also match `suest2`. `suest2` rejects different
  `subpop()` specifications, so the stacked system supplies those joint
  references.

* Added initial multiple-imputation pooling through `suest_mi()`. It applies
  Rubin's rules to compatible complete SUEST systems, retaining cross-model
  covariance in both the within- and between-imputation components. Pooled
  coefficients, covariance, summaries, and coefficient-level hypotheses are
  supported. Inputs may be ordinary lists, `mice::mira` objects from `with()`,
  analysis lists returned by `mice::getfit()`, or result lists from
  `mitools::with.imputationList()`. Predictions, slopes, and comparisons remain
  deferred until each nonlinear estimand can be pooled across imputations
  correctly.
* Corrected the pweighted linear ancillary bread to truncate the effective
  iweight count instead of rounding it. The returned weight-scaling diagnostic
  detects this distinction; regression mean parameters are unaffected.
* Returned Stata benchmarks now verify the full covariance for panel FE/BE,
  balanced RE, ten GEE cases, panel ML, bivariate probit, and IV probit. Added
  permanent cross-language regression fixtures and tests. The closest
  unbalanced Swamy-Arora RE comparison differs from Stata by at most 0.15% on
  a covariance diagonal; unbalanced RE with time indicators remains blocked by
  an upstream singular between design in plm.
* Recompute bivariate- and IV-probit information matrices by differentiating
  analytic score sums. This corrects inaccuracies in the fitting engines'
  Hessians while retaining their point estimates and observation scores.
* Added population-averaged GEE support through `geepack::geeglm()`:
  Gaussian identity, binary logit/probit/cloglog, and Poisson log, with
  independence or exchangeable working correlation. Reconstructed GEE scores
  and bread match native robust covariance before Stata's per-model cluster
  correction. Unequal samples, higher nested clusters, offsets, interactions,
  and `marginaleffects` are tested. Ten single-family Stata systems match the
  full R covariance. A mixed-family unequal-sample Stata case identifies a
  `suest2` score-extraction issue when excluded rows remain within panels.
* Corrected the provisional IV-probit response adapter to evaluate
  `pnorm(X beta)`, as documented, rather than the residual-conditioned
  probability selected by Rchoice's `asf = TRUE`. Independent prediction and
  derivative tests now protect this distinction.
* Retain full-precision random-intercept standard deviations from `nlme`
  rather than reading rounded `VarCorr()` display output, and reject
  random-slope-only fits from the random-intercept ML adapter.
* Added random-intercept Gaussian ML support through
  `nlme::lme(method = "ML")`. The adapter reproduces the panel likelihood's
  fixed-effect and variance-component scores, retains `sigma_u` and `sigma_e`,
  defaults to panel clustering, permits higher nested clusters, and works with
  `marginaleffects`. Returned Stata systems agree within 0.0081% on covariance
  diagonals; about half of that gap is already present in the native model-
  based covariance produced by the two fitting engines.
* Added provisional maximum-likelihood IV probit support through
  `Rchoice::ivpml()`, retaining the structural and reduced-form equations,
  `lnsigma`, and Stata-named `athrho`. Response-scale `marginaleffects` uses
  the average structural probability; link-scale predictions return the
  structural index.
* Added provisional bivariate-probit support through
  `mvProbit::mvProbit()` for models with a common regressor set across the two
  outcome equations. The adapter converts the natural correlation to Stata's
  `athrho` parameterization and predicts the joint-success probability.
* Added provisional linear-panel support through `plm::plm()` for individual
  fixed effects (`model = "within"`), between effects (`model = "between"`),
  and Swamy-Arora random effects (`model = "random"`). The joint covariance
  defaults to panel clustering, permits higher nested clusters, reproduces
  `suest2`'s route-specific finite-sample corrections, and works with
  `marginaleffects`. Cross-language Stata gates are included for final
  numerical verification.
* Added joint cluster-robust covariance through `cluster`, with system-level
  score aggregation, partially overlapping samples, cross-model clusters that
  span disjoint observations, direct or data-column cluster IDs, clustered
  2SLS influence functions, and `marginaleffects` integration. Five
  deterministic Stata 19.5 benchmarks verify identical, partially overlapping,
  disjoint, IV, and pweighted clustered systems.
* Added Stata-compatible `lnvar` ancillary parameters to unweighted and
  pweighted linear models. Unweighted fits use `log(RSS / df.residual)`;
  pweighted fits reproduce `suest2`'s iweight-reference point estimate, score,
  and bread. Ancillary parameters remain excluded from `marginaleffects`
  predictions.
* Matched Stata `suest2`'s specialized 2SLS covariance convention: unweighted
  `fixest::feols()` IV systems now use the native HC0 influence-function
  covariance rather than the ordinary-model `N/(N-1)` correction.
* Added heteroskedastic binary probit and logit support through
  `Rchoice::hetprob()`, including location and log-scale parameters. The
  adapter also repairs the package's prediction failure when the documented
  default probit link is omitted from the original call.
* Added unweighted instrumental-variable 2SLS support through
  `fixest::feols()` without absorbed fixed effects. Tests verify its
  coefficients, instrumented-regressor scores, IV bread, robust covariance,
  sample alignment, and `marginaleffects` integration.
* Added explicit fractional log-log support for quasi-binomial GLMs fitted
  with a user-supplied `link-glm` object. This is an R extension beyond
  `suest2`'s currently documented `fracreg logit` and `fracreg probit` route.
* Added tested direct Tobit support through `censReg::censReg()`, including
  left, right, and two-limit censoring, the log-scale parameter, latent-mean
  predictions for `marginaleffects`, and numerical agreement with the
  equivalent Gaussian `survival::survreg()` likelihood.
* Added tested truncated Gaussian regression support through
  `truncreg::truncreg()`, including the scale parameter and partially
  overlapping left- and right-truncated samples.
* Added zero-inflated Poisson and negative-binomial support through
  `pscl::zeroinfl()`. The negative-binomial adapter includes `log(theta)` and
  supplies its omitted dispersion score and observed-information blocks.
* Added standard beta-regression support through `betareg::betareg()` for
  logit, probit, complementary-log-log, and log-log mean links, including
  modeled precision parameters.
* Added parametric survival and censored-regression support through
  `survival::survreg()` with fixed or estimated common scales. Models with
  stratum-specific scales are rejected because their exposed score matrix does
  not align with the fitted parameter vector. Gaussian interval regression is
  tested explicitly as the R analogue of Stata's `intreg`.
* Added joint scalar--categorical systems and categorical models with different
  outcome levels. `marginaleffects` output uses model-qualified category labels
  and does not imply that different response scales are commensurate.
* Expanded `stats::glm()` and `glm2::glm2()` support to Gaussian, Gamma, and
  other families using supported response links; binary complementary-log-log;
  and fractional-response quasi-binomial logit, probit, and cloglog models.
* Added offsets for linear, binary, count, ordered, multinomial, and supported
  `ordinal::clm()` models, including pweighted custom-score calculations.
* Applied the observed-information binary score calculation to unweighted as
  well as pweighted probit models for closer Stata `suest` parity.
* Extended pweights to negative-binomial, ordered logit/probit, and multinomial
  logit models, completing the model-family coverage of Stata `suest2`'s
  ordinary pweight route.
* Extended `suest()` from exactly two models to two or more models, including
  heterogeneous scalar and categorical systems, pairwise overlap tracking,
  and `marginaleffects` predictions and comparisons for every component.
* Added `observation_id` to align overlapping observations explicitly across
  models fitted from different data objects. It accepts one or more ID column
  names or model-specific ID vectors and supports composite panel identifiers.
* Expanded valid heterogeneous combinations to every pair of supported
  scalar-response families and every pair of supported categorical-response
  families with matching categories. This adds, among others, ordered
  logit--ordered probit, ordered probit--multinomial, logit--Poisson, and
  linear--negative-binomial combinations.
* Clarified that comparisons spanning different scalar response scales report
  model-specific response values rather than a common probability or count
  scale.
* Extended `weight_type = "pweight"` to Poisson log-link models fitted with
  `stats::glm()` or `glm2::glm2()`.
* Matched Stata's identical, partially overlapping, and disjoint Poisson
  coefficient and joint-covariance benchmarks.
* Added Poisson pweight score, covariance, scaling, overlap, disjoint-sample,
  alternative-engine, and marginaleffects tests.
* Extended `weight_type = "pweight"` to binary logit and probit models fitted
  with `stats::glm()` or `glm2::glm2()`, including supported logit-probit,
  linear-logit, and linear-probit combinations.
* Matched Stata's identical, partially overlapping, disjoint, and cross-family
  binary-model joint covariance benchmarks. The probit implementation uses
  the observed information matrix required for Stata parity.
* Added weighted binary marginaleffects, scaling-invariance, overlap,
  disjoint-sample, and alternative-engine tests.
* Corrected the pweight scaling-invariance unit test to use fixed model
  labels, so the test compares numerical results rather than object names.
* Corrected the focused pweight test harness so its testthat expectations are
  not masked by the separate numerical-acceptance helper functions.
* Corrected the `suest_newdata()` marginaleffects test to compare results by
  model name rather than relying on output row order.
* Added the first stage of explicit sampling-weight support through
  `weight_type = "pweight"` for `stats::lm()` models.
* Matched Stata `suest` weighted linear-model coefficient and joint-covariance
  benchmarks, including identical and partially overlapping samples.
* Required evaluated pweights to agree on observations shared by both models,
  while allowing different weights on model-specific observations and in
  completely disjoint samples.
* Added scale-invariant weighted score and bread calculations, Stata's
  union-sample HC1 correction for different estimation samples, strict
  validation of finite positive pweights, and informative errors for
  undeclared or currently unsupported weighted models.
* Added `.suest_weight` to `suest_newdata()` for weighted averaging with the
  `wts` argument in `marginaleffects`.
* Added unit and numerical acceptance tests using the reproducible Stata
  pweight benchmark data.

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
