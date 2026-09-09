# Combine fitted models with seemingly unrelated estimation

`suest()` combines two or more separately fitted models and constructs a
joint model-robust covariance matrix from their observation-level score
contributions. The returned object can be passed directly to
[`marginaleffects::predictions()`](https://rdrr.io/pkg/marginaleffects/man/predictions.html),
[`marginaleffects::avg_comparisons()`](https://rdrr.io/pkg/marginaleffects/man/comparisons.html),
[`marginaleffects::avg_slopes()`](https://rdrr.io/pkg/marginaleffects/man/slopes.html),
and
[`marginaleffects::hypotheses()`](https://rdrr.io/pkg/marginaleffects/man/hypotheses.html).

## Usage

``` r
suest(
  ...,
  model_names = NULL,
  observation_id = NULL,
  cluster = NULL,
  weight_type = NULL,
  survey_design = NULL
)

# S3 method for class 'suest_model'
coef(object, ...)

# S3 method for class 'suest_model'
vcov(object, ...)

# S3 method for class 'suest_model'
nobs(object, ...)

# S3 method for class 'suest_model'
print(x, ...)
```

## Arguments

- ...:

  Two or more supported fitted model objects. For backward
  compatibility, a character vector supplied as the third unnamed
  argument is interpreted as `model_names` for a two-model system.

- model_names:

  Optional character vector containing one display name per model. By
  default, the object names supplied in the call are used.

- observation_id:

  Optional observation identifier used to align models fitted from
  different data objects. Supply one or more column names found in every
  model's original data, such as `"id"` or `c("id", "wave")`, or a list
  containing an ID vector, matrix, or data frame already aligned to each
  model's estimation sample. IDs must be complete and unique within each
  model. The default `NULL` uses data-source identity and model-frame
  row names.

- cluster:

  Optional cluster identifier for a joint cluster-robust covariance
  matrix. Supply one or more column names found in every model's
  original data, or a list containing a cluster vector, matrix, or data
  frame aligned to each model's estimation sample. Cluster IDs must be
  complete. Shared observations must have the same cluster ID in every
  model; disjoint observations may still share clusters across models.
  Supported panel systems default to the panel identifier when `cluster`
  is omitted. Supplied clusters for these models must contain whole
  panels.

- weight_type:

  Optional weight interpretation. The default `NULL` preserves the
  unweighted behavior and rejects nonunit estimation weights. Use
  `"pweight"` to treat model weights as sampling weights. In version
  0.1.4, pweights are supported for linear, binary logit/probit,
  Poisson, negative-binomial, ordered logit/probit, and multinomial
  logit models.

- survey_design:

  Full common one-stage
  [`survey::svydesign()`](https://rdrr.io/pkg/survey/man/svydesign.html)
  object for supported
  [`survey::svyglm()`](https://rdrr.io/pkg/survey/man/svyglm.html)
  models: Gaussian identity or binary
  [`quasibinomial()`](https://rdrr.io/r/stats/family.html) logit/probit.
  Requires observation-ID column names. Retain the full design before
  domain subsetting and specify weights and clusters in the design,
  without `weight_type` or `cluster`.

- object, x:

  A `"suest_model"` object.

## Value

An object of class `"suest_model"` containing the fitted models, their
joint coefficient vector, a joint model-robust covariance matrix, and
sample-alignment information.

## Details

Two or more models are supported. Models can use identical, partially
overlapping, or completely disjoint samples. When the model calls refer
to the same data source, model-frame row names identify overlapping
observations. Models fitted from different data objects are treated as
disjoint by default because shared observations cannot be inferred
safely. Use `observation_id` to identify their common observations
explicitly.

All combinations of the supported scalar-response models can be
combined, including models whose response variables or response scales
differ. All combinations of the supported categorical-response models
can be combined, and scalar and categorical models may appear in the
same system. Results on different response scales are labeled separately
for `marginaleffects`.

When `cluster` is supplied, observation-level score contributions are
summed within the system-level clusters before the joint covariance is
formed. Clusters can span observations from different, even disjoint,
model samples.

Ordinary linear models include an ancillary `lnvar` parameter, matching
Stata's `regress`/`suest` parameterization. Unweighted fits use
`log(RSS / df.residual)`; pweighted fits follow `suest2`'s reconstructed
iweight-reference normalization. Negative-binomial models include
`log(theta)` in the joint parameter vector. Ordered and multinomial
models use analytic score and observed-information calculations for
stable robust covariance estimation.

Aliased parameters are not supported. Nonunit weights are rejected
unless `weight_type = "pweight"`. Pweights must be finite and strictly
positive. For observations included in any pair of models, evaluated
weights must agree; weights may differ for observations unique to either
model, including completely disjoint samples. Pweight support is
available for linear, binary logit/probit, Poisson, negative-binomial,
ordered logit/probit, and multinomial logit models.

Bias-reduced, adjusted-score, Firth, and penalized GLM fits are rejected
because they do not use ordinary maximum-likelihood score equations.

## Survey models

Survey support combines coefficients from Gaussian identity-link and
binary [`quasibinomial()`](https://rdrr.io/r/stats/family.html)
logit/probit [`svyglm()`](https://rdrr.io/pkg/survey/man/svyglm.html)
fits under one common one-stage design. Gaussian and binary fits are not
mixed in the same survey system. Strata and optional first-stage
finite-population corrections are supported. Model-specific domains and
missing outcomes are aligned by observation IDs, with zero influence
outside each estimation sample. The full design retains PSUs outside all
model samples. Each native coefficient covariance must be reproduced
before the joint matrix is returned. Survey fits contain coefficients
only and do not add ancillary parameters.

Replicate-weight, multistage, two-phase, calibrated, raked,
post-stratified, and PPS designs are unsupported. Lonely-PSU options are
restricted to `fail`, `remove`, and `certainty`, with
`survey.adjust.domain.lonely = FALSE`. A singleton stratum with a
certainty FPC is supported under `fail`. Extra fitting weights are
unsupported; place offsets in the model formula.

Predictions and effects use the joint design-based coefficient
covariance. Averaging treats the supplied covariate distribution as
fixed and does not add design uncertainty from estimating that
distribution. Use
[`suest_newdata()`](https://tdmize.github.io/suest/reference/suest_newdata.md)
and `wts = ".suest_weight"` for model-specific weighted averages.
Default inference is asymptotic normal. The full design degrees of
freedom are recorded in `$survey$design_df`; no automatic survey t or F
adjustment is applied.

## Supported models

- [`stats::lm()`](https://rdrr.io/r/stats/lm.html)

- Restricted survey-weighted Gaussian identity and binary
  [`quasibinomial()`](https://rdrr.io/r/stats/family.html) logit/probit
  models from
  [`survey::svyglm()`](https://rdrr.io/pkg/survey/man/svyglm.html)

- binary logit, probit, and complementary-log-log models from
  [`stats::glm()`](https://rdrr.io/r/stats/glm.html) or
  [`glm2::glm2()`](https://rdrr.io/pkg/glm2/man/glm2.html)

- Poisson log-link models from
  [`stats::glm()`](https://rdrr.io/r/stats/glm.html) or
  [`glm2::glm2()`](https://rdrr.io/pkg/glm2/man/glm2.html)

- other GLMs using identity, log, logit, probit, complementary-log-log,
  or log-log links

- fractional-response GLMs using
  [`quasibinomial()`](https://rdrr.io/r/stats/family.html) with logit,
  probit, complementary-log-log, or a user-supplied log-log link

- negative-binomial log-link models from
  [`MASS::glm.nb()`](https://rdrr.io/pkg/MASS/man/glm.nb.html)

- parametric survival and censored-regression models from
  [`survival::survreg()`](https://rdrr.io/pkg/survival/man/survreg.html)
  with a common scale, including Gaussian interval regression

- beta regressions from
  [`betareg::betareg()`](https://rdrr.io/pkg/betareg/man/betareg.html)
  using logit, probit, complementary-log-log, or log-log mean links

- Poisson and negative-binomial zero-inflated models from
  [`pscl::zeroinfl()`](https://rdrr.io/pkg/pscl/man/zeroinfl.html)

- truncated Gaussian regressions from
  [`truncreg::truncreg()`](https://rdrr.io/pkg/truncreg/man/truncreg.html)

- left-, right-, and two-limit censored Gaussian regressions from
  [`censReg::censReg()`](https://rdrr.io/pkg/censReg/man/censReg.html);
  response-scale predictions are the latent mean

- unweighted two-stage least squares from
  [`fixest::feols()`](https://lrberge.github.io/fixest/reference/feols.html)
  without absorbed fixed effects; the original data object must remain
  available

- heteroskedastic binary probit and logit from
  [`Rchoice::hetprob()`](https://rdrr.io/pkg/Rchoice/man/hetprob.html)

- maximum-likelihood instrumental-variable probit from
  [`Rchoice::ivpml()`](https://rdrr.io/pkg/Rchoice/man/ivpml.html);
  response predictions use the average structural probability

- bivariate probit from
  [`mvProbit::mvProbit()`](https://rdrr.io/pkg/mvProbit/man/mvProbit.html)
  when both equations use the same regressors; fit with `intGrad = TRUE`
  and `finalHessian = TRUE`

- unweighted individual fixed-effects, between-effects, and Swamy-Arora
  random-effects linear panel models from
  [`plm::plm()`](https://rdrr.io/pkg/plm/man/plm.html); balanced random-
  effects panels have the closest Stata parity, while unbalanced panels
  can retain small engine-specific differences

- unweighted single-level random-intercept Gaussian panel models from
  [`nlme::lme()`](https://rdrr.io/pkg/nlme/man/lme.html) fitted with
  `method = "ML"`

- unweighted GEE from
  [`geepack::geeglm()`](https://rdrr.io/pkg/geepack/man/geeglm.html):
  Gaussian identity, binary logit/probit/cloglog, and Poisson log, with
  independence or exchangeable correlation and numeric outcomes

- ordered logit and probit models from
  [`MASS::polr()`](https://rdrr.io/pkg/MASS/man/polr.html)

- ordered logit and probit models from
  [`ordinal::clm()`](https://rdrr.io/pkg/ordinal/man/clm.html) with
  flexible thresholds, proportional effects, and no scale model

- multinomial logit models from
  [`nnet::multinom()`](https://rdrr.io/pkg/nnet/man/multinom.html)

## Examples

``` r
dat <- mtcars
dat$am <- factor(dat$am)

model1 <- glm(am ~ wt, family = binomial(), data = dat)
model2 <- glm(am ~ wt + hp, family = binomial(), data = dat)
fit <- suest(model1, model2, model_names = c("Base", "Adjusted"))
fit
#> Seemingly Unrelated Estimation
#> Models: Base + Adjusted 
#> Model types: Base=logit, Adjusted=logit 
#> Model engines: Base=stats::glm, Adjusted=stats::glm 
#> Comparison scale: predicted probabilities 
#> Observations: Base=32, Adjusted=32 
#> Overlapping observations: 32 
#> Union observations: 32 
#> Parameters: 5 

effects <- marginaleffects::avg_comparisons(fit, variables = "wt", newdata = dat)
marginaleffects::hypotheses(effects, hypothesis = difference ~ revpairwise)
#> 
#>           Hypothesis Estimate Std. Error    z Pr(>|z|)    S  2.5 % 97.5 %
#>  (Base) - (Adjusted)   0.0779     0.0141 5.54   <0.001 25.0 0.0503  0.105
#> 
#> 
```
