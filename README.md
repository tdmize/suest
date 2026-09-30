# suest for R

[![R-CMD-check](https://github.com/tdmize/suest/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/tdmize/suest/actions/workflows/R-CMD-check.yaml)

[Getting started, with worked examples and output](https://tdmize.github.io/suest/articles/suest.html)

`suest` lets you compare predictions and marginal effects across regression
models in R. It combines two or more separately fitted models into one object
with a joint robust covariance matrix. That object works directly with
[`marginaleffects`](https://marginaleffects.com/), so you can calculate
predictions and effects for each model and test whether they differ across
models.

The package implements the framework in:

> Mize, Trenton D., Long Doan, and J. Scott Long. 2019. “A General Framework for Comparing Predictions and Marginal Effects Across Models.” *Sociological Methodology* 49(1):152–189. [doi:10.1177/0081175019852763](https://doi.org/10.1177/0081175019852763)

The same framework is available in Stata through the
[`mecompare` command](https://www.trentonmize.com/software/mecompare).

## Installation

```r
# install.packages("remotes")
remotes::install_github("tdmize/suest")
```

## Basic usage

```r
library(suest)
library(marginaleffects)

dat <- mtcars
dat$am <- factor(dat$am)

base <- glm(am ~ wt, family = binomial("logit"), data = dat)
adjusted <- glm(am ~ wt + hp, family = binomial("logit"), data = dat)
combined <- suest(base, adjusted, model_names = c("Base", "Adjusted"))

effects <- avg_comparisons(combined, variables = "wt", newdata = dat)
effects
hypotheses(effects, hypothesis = difference ~ revpairwise)
```

`suest()` accounts for the covariance between estimates from different
models. This is essential when the models use the same or overlapping
observations.

## Supported models

- Linear regression: `lm()`
- Binary logit, probit, and complementary log-log: `glm()` or `glm2::glm2()`
- Poisson regression: `glm()` or `glm2::glm2()`
- Other GLMs with identity, log, logit, probit, complementary log-log, or
  log-log links
- Fractional response (quasi-binomial) with logit, probit, complementary
  log-log, or user-supplied log-log links
- Negative binomial regression: `MASS::glm.nb()`
- Parametric survival and censored regression with a common scale, including
  Gaussian interval regression: `survival::survreg()`
- Beta regression with logit, probit, complementary log-log, or log-log mean
  links: `betareg::betareg()`
- Zero-inflated Poisson and negative binomial regression: `pscl::zeroinfl()`
- Truncated Gaussian regression: `truncreg::truncreg()`
- Censored Gaussian (Tobit) regression with left, right, or two limits:
  `censReg::censReg()`. Response-scale predictions are the latent mean.
- Two-stage least squares, unweighted and without absorbed fixed effects:
  `fixest::feols()`. The original data object must remain available.
- Heteroskedastic binary probit and logit: `Rchoice::hetprob()`
- Maximum-likelihood instrumental-variable probit: `Rchoice::ivpml()`.
  Response-scale predictions use the average structural probability.
- Bivariate probit: `mvProbit::mvProbit()`, when both equations use the same
  regressors and the final Hessian is kept (`finalHessian = TRUE`)
- Linear panel models with individual fixed, between, or Swamy-Arora random
  effects (unweighted): `plm::plm()`. Comparisons of slopes and changes are
  supported. Confidence intervals and tests for average predicted levels at
  the estimation-sample means are not (their standard errors can be missing
  or near zero). Unbalanced random-effects panels can differ slightly from
  Stata.
- Random-intercept Gaussian panel models fitted by maximum likelihood
  (unweighted): `nlme::lme(method = "ML")`. Both variance components are
  included.
- Random-intercept binary logit or probit (unweighted, individual effects,
  exactly 12-point nonadaptive quadrature):
  `pglm::pglm(model = "random", effect = "individual")`. Predictions
  integrate over the random effect.
- Gamma random-effects panel Poisson with a log link (unweighted):
  `pglm::pglm(model = "random", effect = "individual", other = "sd")`. The
  gamma variance `alpha` is included.
- Mixed models with one Gaussian random intercept (unweighted):
  `glmmTMB::glmmTMB()` with binomial-logit, Poisson-log, or negative binomial
  (NB2) log with constant dispersion (`dispformula = ~1`). Predictions
  integrate over the random effect.
- Mixed models with a correlated random intercept and numeric slope,
  `(1 + x | id)` (unweighted): `glmmTMB::glmmTMB()` with binomial-logit,
  Poisson-log, or NB2-log. Predictions and marginal effects integrate over
  both random effects.
- Population-averaged GEE (unweighted): `geepack::geeglm()` with Gaussian
  identity, binary logit, probit, or complementary log-log, or Poisson log;
  independence or exchangeable correlation; numeric outcomes only; the
  default panel clusters or higher nested clusters
- Ordered logit and probit: `MASS::polr()`, or restricted `ordinal::clm()`
  specifications
- Multinomial logit: `nnet::multinom()`

Models can use the same, partially overlapping, or completely different
samples. Any supported models with a single outcome can be combined, and so
can models for categorical outcomes; the two kinds can also be mixed. Results
on different scales are labeled separately rather than treated as directly
comparable. Panel models can be combined only with panel models of the same
type. Offsets are supported.

As in Stata, linear models add an `lnvar` parameter (the log of the error
variance) to the combined coefficients and covariance matrix. Unweighted fits
use `log(RSS / df.residual)`; weighted fits follow `suest2`'s normalization.
It has no effect on predictions or marginal effects.

### Notes on the glmmTMB mixed models

- Weights, offsets, and zero inflation are not supported for these models.
  Other random-slope families and more complex random-effects structures are
  not supported.
- Center continuous predictors and keep them in well-scaled units. glmmTMB's
  covariance matrix can be inaccurate with extreme predictor units even when
  the model reports convergence, and `suest` uses that covariance as supplied.
  The refit and sensitivity workflow is in
  `tools/STABILITY-USABILITY-20260925.md`.
- For numerical slope standard errors in logit random-slope models, use
  `numderiv = list("fdcenter", eps = 1e-4)` in `avg_slopes()` and check
  that nearby step sizes give similar answers. This does not fix an inaccurate
  model covariance matrix.
- These models were checked against Stata. The validation reports document
  the comparisons (for the random-slope models, with identical, partially
  overlapping, and completely different samples) and the remaining numerical
  differences:
  [NB2 random intercept](https://github.com/tdmize/suest/blob/main/tools/GLMMTMB-NBINOM2-RI-VALIDATION-20260923.md),
  [Poisson random slope](https://github.com/tdmize/suest/blob/main/tools/GLMMTMB-POISSON-RS-VALIDATION-20260923.md),
  [NB2 random slope](https://github.com/tdmize/suest/blob/main/tools/GLMMTMB-NBINOM2-RS-VALIDATION-20260924.md), and
  [logit random slope](https://github.com/tdmize/suest/blob/main/tools/GLMMTMB-LOGIT-RS-VALIDATION-20260924.md).

### Fits that are not supported

- Bias-reduced, adjusted-score, Firth, and other penalized GLM fits, including
  `brglm2::brglmFit()`, are rejected rather than treated as ordinary
  maximum-likelihood GLMs.
- Ordinary `glm` and `glm2` fits must have converged.
- Penalized `nnet::multinom` fits (nonzero `decay`) are rejected.
- Beta regression must use ordinary maximum likelihood (`type = "ML"`);
  bias-reduced (`"BR"`) and bias-corrected (`"BC"`) fits are not supported.
- Penalized `survreg` fits are not supported. A regular `survreg` fit with
  `robust = TRUE` is supported: `suest` uses its model-based information when
  building the joint covariance matrix.

## Multiple imputation

To combine results across imputed datasets, fit the same `suest` system in
each one and pool them with `suest_mi()`:

```r
pooled <- suest_mi(list(fit_imp1, fit_imp2, fit_imp3, fit_imp4, fit_imp5))
summary(pooled)
```

Results from `mice::with()` can be pooled directly when each imputation fits
the same system:

```r
analyses <- with(imp, suest(
  lm(y ~ x),
  lm(y ~ x + z),
  model_names = c("Base", "Adjusted")
))
pooled <- suest_mi(analyses)
```

The result of `with(mitools::imputationList(...), suest(...))` can also be
passed directly to `suest_mi()`.

`suest_mi()` applies Rubin's rules to the full set of coefficients, keeping
the covariance between models both within and between imputations. Tests of
coefficients and their pooled standard errors work through
`marginaleffects::hypotheses()`. Predictions, slopes, and comparisons are not
yet available for pooled results, because those need to be calculated within
each imputation and then pooled.

## Cluster-robust covariance

Use `cluster` to get cluster-robust standard errors for the combined models:

```r
combined <- suest(base, adjusted, cluster = "person_id")
```

`cluster` takes one or more column names, or a list of cluster IDs matched to
each model's estimation sample. An observation that appears in more than one
model must have the same cluster ID in each. Observations in completely
different samples can share a cluster, so the covariance between models need
not be zero. The result works with `marginaleffects` as usual. Panel models
cluster on the panel identifier by default; a cluster you supply must contain
whole panels.

## Probability weights

`suest` supports probability weights. Fit each model with `weights=`, and tell
`suest()` how to interpret the weights:

```r
base <- lm(y ~ x + z, weights = pw, data = dat)
adjusted <- lm(y ~ x + z + mediator, weights = pw, data = dat)
combined <- suest(base, adjusted, weight_type = "pweight")
```

Weights must be positive and finite. An observation that appears in more than
one model must have the same weight in each; weights can differ for
observations that appear in only one model, including completely different
samples. To average over the sample with the weights in `marginaleffects`,
give the weight column to `wts`, such as `wts = "pw"`. When the models have
different estimation samples, use `suest_newdata(combined)` and
`wts = ".suest_weight"`.

Probability weights work with linear, binary logit and probit, Poisson,
negative binomial, ordered logit and probit, and multinomial logit models.
They are not yet supported for complementary log-log, other GLMs, or
fractional response models.

## Worked examples

The Getting started page has code and output for six R replications of
Examples 6.1–6.6 from the
[`mecompare` Stata command](https://www.trentonmize.com/software/mecompare):

1. Marginal effects to summarize curvilinear relationships and test mediation
2. Comparing marginal effects across nested logit models
3. Comparing marginal effects using alternative predictors
4. Comparing marginal effects across different outcomes
5. Comparing marginal effects across different model types (ordinal vs nominal)
6. Comparing marginal effects across different samples or groups

It also shows comparisons across fitting engines using `glm2::glm2()` and
`ordinal::clm()`.

## Survey-weighted models

`suest` supports Gaussian (identity) and binary `quasibinomial()` logit or
probit models fitted with `survey::svyglm()` from one common one-stage survey
design:

```r
library(survey)
design <- svydesign(~psu, strata = ~strata, weights = ~w, data = dat)
m1 <- svyglm(yb ~ x, design, family = quasibinomial("logit"))
m2 <- svyglm(yb ~ x + z, design, family = quasibinomial("logit"))
fit <- suest(m1, m2, survey_design = design, observation_id = "id")
avg_slopes(fit, variables = "x", newdata = suest_newdata(fit), wts = ".suest_weight")
```

If either model uses a subset (domain) of the data, keep `design` as the full
sampled design. The joint covariance accounts for strata, PSUs, and optional
one-stage finite population corrections.

Some limits:

- Survey systems include coefficients only, and Gaussian and binary models
  can't be mixed in the same system.
- Replicate weights, multistage designs, calibration, PPS designs, and other
  survey families are not yet supported.
- Lonely-PSU settings are limited to `fail`, `remove`, or `certainty`, with
  domain-lonely adjustments turned off.
- Predictions and effects use the joint coefficient covariance and treat the
  covariates being averaged over as fixed. Inference is asymptotic normal; the
  design degrees of freedom are stored in `fit$survey$design_df`.

The Gaussian and logit survey results match Stata to machine precision. For
survey probit, the coefficients match Stata, but the standard errors can
differ slightly: `survey::svyglm()` uses expected (Fisher) information, while
Stata uses observed information. `suest` keeps `svyglm()`'s covariance.
