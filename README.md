# suest for R

[![R-CMD-check](https://github.com/tdmize/suest/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/tdmize/suest/actions/workflows/R-CMD-check.yaml)

Full documentation, worked examples, and rendered output:
[tdmize.github.io/suest](https://tdmize.github.io/suest/articles/suest.html)

`suest` combines two or more separately fitted regression models into one object with a
joint model-robust covariance matrix. The combined object works with
[`marginaleffects`](https://marginaleffects.com/) to compare predictions and
marginal effects across models.

The package implements the framework developed in:

> Mize, Trenton D., Long Doan, and J. Scott Long. 2019. “A General Framework for Comparing Predictions and Marginal Effects Across Models.” *Sociological Methodology* 49(1):152–189. [doi:10.1177/0081175019852763](https://doi.org/10.1177/0081175019852763)

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

`suest()` accounts for the cross-model covariance between estimates. This is
essential when the models use the same or overlapping observations.

## Supported models

- linear regression using `lm()`
- binary logit, probit, and complementary-log-log using `glm()` or `glm2::glm2()`
- Poisson regression using `glm()` or `glm2::glm2()`
- other GLMs using identity, log, logit, probit, complementary-log-log, or
  log-log links
- fractional-response quasi-binomial logit, probit, complementary-log-log,
  and user-supplied log-log links
- negative-binomial regression using `MASS::glm.nb()`
- parametric survival and censored regression using `survival::survreg()` with
  a common scale, including Gaussian interval regression
- beta regression using `betareg::betareg()` with logit, probit,
  complementary-log-log, or log-log mean links
- Poisson and negative-binomial zero-inflated regression using
  `pscl::zeroinfl()`
- truncated Gaussian regression using `truncreg::truncreg()`
- left-, right-, and two-limit censored Gaussian regression using
  `censReg::censReg()`; response-scale predictions are the latent mean
- unweighted two-stage least squares using `fixest::feols()` without absorbed
  fixed effects; the original data object must remain available
- heteroskedastic binary probit and logit using `Rchoice::hetprob()`
- maximum-likelihood instrumental-variable probit using `Rchoice::ivpml()`;
  response-scale predictions use the average structural probability
- bivariate probit using `mvProbit::mvProbit()` when both outcome equations
  use the same regressors and an observed final Hessian is retained
- unweighted individual fixed-effects, between-effects, and Swamy-Arora
  random-effects linear panel models using `plm::plm()`; unbalanced random-
  effects panels can retain small engine-specific differences from Stata
- unweighted random-intercept Gaussian panel models using `nlme::lme()` with
  `method = "ML"`; both variance components are included in the joint system
- unweighted GEE using `geepack::geeglm()`: Gaussian identity,
  binary logit/probit/cloglog, and Poisson log, with independence or exchangeable
  correlation; numeric outcomes only, with default panel or higher nested clusters
- ordered logit and probit using `MASS::polr()`
- ordered logit and probit using restricted `ordinal::clm()` specifications
- multinomial logit using `nnet::multinom()`

Models may use identical, partially overlapping, or disjoint samples. All
combinations of supported scalar-response models are allowed. All combinations
of supported categorical-response models are allowed, and scalar and
categorical models may be combined. Outputs on different response scales are
labeled separately rather than treated as directly commensurate. Supported
panel systems can be combined only with panel models of the same type.
Offsets are supported. Linear models include a Stata-compatible `lnvar`
ancillary parameter in the combined coefficient and covariance matrices.
Unweighted fits use `log(RSS / df.residual)`; pweighted fits follow `suest2`'s
iweight-reference normalization. This nuisance parameter has no effect on
`marginaleffects` predictions.

## Multiple imputation

Pool compatible SUEST systems fitted in separate imputed datasets with
`suest_mi()`:

```r
pooled <- suest_mi(list(fit_imp1, fit_imp2, fit_imp3, fit_imp4, fit_imp5))
summary(pooled)
```

`mice::with()` results can be pooled directly when each imputation fits the
same SUEST system:

```r
analyses <- with(imp, suest(
  lm(y ~ x),
  lm(y ~ x + z),
  model_names = c("Base", "Adjusted")
))
pooled <- suest_mi(analyses)
```

The result of `with(mitools::imputationList(...), suest(...))` can likewise be
passed directly to `suest_mi()`.

The function applies Rubin's rules to the complete joint coefficient vector,
preserving cross-model covariance within and between imputations. Coefficient
hypotheses and their pooled standard errors work through
`marginaleffects::hypotheses()`. Predictions, slopes, and comparisons are not
yet exposed for pooled objects because nonlinear estimands must be estimated
within each imputation and then pooled.

## Cluster-robust covariance

Use `cluster` to aggregate the joint score contributions before constructing
the covariance matrix:

```r
combined <- suest(base, adjusted, cluster = "person_id")
```

The argument accepts one or more cluster-column names, or a list of cluster IDs
aligned to the estimation sample of each model. Shared observations must have
the same cluster ID in every model. Completely disjoint observations can share
a cluster, so their cross-model covariance need not be zero. The resulting
object remains directly compatible with `marginaleffects`. Supported panel
systems default to clustering on the panel identifier; a supplied higher-level
cluster must contain whole panels.

## Probability weights

Version 0.1.4 adds explicit probability-weight support. Fit each component model with
`weights=` and declare their interpretation when combining the models:

```r
base <- lm(y ~ x + z, weights = pw, data = dat)
adjusted <- lm(y ~ x + z + mediator, weights = pw, data = dat)
combined <- suest(base, adjusted, weight_type = "pweight")
```

Pweights must be finite and strictly positive. Their evaluated values must agree
for observations included in any pair of models. They may differ for observations
unique to one model, including completely disjoint samples. To calculate a
sampling-weighted average with `marginaleffects`, supply the weight column using
`wts`, such as `wts = "pw"`. With model-specific estimation samples, use
`suest_newdata(combined)` and `wts = ".suest_weight"`.

Pweights are supported for linear, binary logit/probit, Poisson,
negative-binomial, ordered logit/probit, and multinomial logit models. They are
not yet supported for complementary-log-log, other general GLMs, or fractional
response models.

## Worked examples

The single Get Started page contains code and rendered results for six R
replications corresponding to Examples 6.1–6.6 for the
[`mecompare` Stata command](https://www.trentonmize.com/software/mecompare):

1. Marginal effects to summarize curvilinear relationships and test mediation
2. Comparing marginal effects across nested logit models
3. Comparing marginal effects using alternative predictors
4. Comparing marginal effects across different outcomes
5. Comparing marginal effects across different model types (ordinal vs nominal)
6. Comparing marginal effects across different samples or groups

It also documents cross-engine comparisons using `glm2::glm2()` and
`ordinal::clm()`.

Bias-reduced, adjusted-score, Firth, and penalized GLM fits, including
`brglm2::brglmFit()`, are rejected rather than being treated incorrectly as
ordinary maximum-likelihood GLMs.


### Restricted survey-weighted models

Fit Gaussian identity or binary `quasibinomial()` logit/probit models from one
common one-stage survey design:

```r
library(survey)
design <- svydesign(~psu, strata = ~strata, weights = ~w, data = dat)
m1 <- svyglm(yb ~ x, design, family = quasibinomial("logit"))
m2 <- svyglm(yb ~ x + z, design, family = quasibinomial("logit"))
fit <- suest(m1, m2, survey_design = design, observation_id = "id")
avg_slopes(fit, variables = "x", newdata = suest_newdata(fit), wts = ".suest_weight")
```

Keep `design` as the full sampled design if either model uses a domain/subset.
The joint covariance incorporates strata, PSU totals, and optional one-stage
FPCs. Survey systems contain coefficients only; Gaussian and binary fits are
not mixed in the same survey system. Replicate weights, multistage designs,
calibration, PPS designs, and other survey families remain deferred.
Lonely-PSU settings are restricted to `fail`, `remove`, or `certainty`, with
domain-lonely adjustments disabled. Predictions/effects use the joint
coefficient covariance and treat averaging covariates as fixed. Default
inference remains asymptotic normal; the design DF is recorded in
`fit$survey$design_df`.

The certified Gaussian route matches six deterministic stacked Stata systems
to machine precision. The binary gate shows the same numerical parity for
survey logit. Survey probit coefficient estimates agree across R and Stata,
but their finite-sample sandwich covariances use different breads:
`survey::svyglm()` uses expected/Fisher information, while the returned Stata
covariance is reproduced by an observed-information bread. The R implementation
preserves the native `svyglm()` covariance.
