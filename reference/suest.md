# Combine fitted models with seemingly unrelated estimation

`suest()` combines two or more separately fitted models into one object
with a joint robust covariance matrix, so you can test whether
predictions and marginal effects differ across the models. Pass the
result directly to
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
  is omitted, and systems that combine random-effects models with
  ordinary models default to the random-effects grouping variable.
  Supplied clusters for these models must contain whole panels.

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
  models: Gaussian identity, binary
  [`quasibinomial()`](https://rdrr.io/r/stats/family.html)
  logit/probit/cloglog, or
  [`quasipoisson()`](https://rdrr.io/r/stats/family.html) log. Requires
  `observation_id` column names. Retain the design before model-specific
  domain subsetting. Specify weights and clusters in this design,
  without `weight_type` or `cluster`. Ordinary models leave this
  argument `NULL`.

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
explicitly. `suest()` warns when such models share row names with
identical values, and stops when models fit to the same data object
disagree on rows with the same row names (for example, after the data
were filtered and renumbered).

When `newdata` is omitted in `marginaleffects` functions, each model is
averaged over its own estimation sample, as with
`newdata = suest_newdata(fit)`.

All combinations of the supported scalar-response models can be
combined, including models whose response variables or response scales
differ. All combinations of the supported categorical-response models
can be combined, and scalar and categorical models may appear in the
same system. Results on different response scales are labeled separately
for `marginaleffects`.

Linear models include an ancillary `lnvar` parameter, the log of the
error variance (`log(RSS / df.residual)` for unweighted fits). It has no
effect on predictions or marginal effects. Negative-binomial models
include `log(theta)` in the joint parameter vector. Ordered and
multinomial models use analytic score and observed-information
calculations for stable robust covariance estimation. Probit and other
GLMs with non-canonical links, and beta regressions, use observed
information in the sandwich bread, as Stata does.

Aliased parameters are not supported. Nonunit weights are rejected
unless `weight_type = "pweight"`. Pweights must be finite and strictly
positive. For observations included in any pair of models, evaluated
weights must agree; weights may differ for observations unique to either
model, including completely disjoint samples. Pweight support is
available for linear, binary logit/probit, Poisson, negative-binomial,
ordered logit/probit, and multinomial logit models.

Bias-reduced, adjusted-score, Firth, and penalized GLM fits are rejected
because they do not use ordinary maximum-likelihood score equations.
Ordinary `glm` and `glm2` fits must have converged. Penalized
[`nnet::multinom`](https://rdrr.io/pkg/nnet/man/multinom.html) fits with
nonzero `decay` are also rejected. Beta regression requires
`type = "ML"`; bias-reduced and bias-corrected beta fits and penalized
`survreg` fits are not supported. Conventional `survreg` fits with
`robust = TRUE` use their model-based information for the joint sandwich
covariance.

## Survey models

Survey support combines coefficients from Gaussian identity-link, binary
[`quasibinomial()`](https://rdrr.io/r/stats/family.html)
logit/probit/cloglog, and
[`quasipoisson()`](https://rdrr.io/r/stats/family.html) log `svyglm()`
fits under one common one-stage design, in any combination. Strata and
first-stage finite-population corrections are supported. Model-specific
subsets and missing outcomes are aligned using observation IDs, with
zero influence outside each model's estimation sample. The full design
retains PSUs outside all model samples. Each native coefficient
covariance must be reproduced before the joint matrix is returned.
Survey fits contain coefficients only and do not add ancillary
parameters.

Replicate-weight, multistage, two-phase, calibrated, raked,
post-stratified, and PPS designs are unsupported. Lonely-PSU options are
restricted to `fail`, `remove`, or `certainty`, with
`survey.adjust.domain.lonely = FALSE`. A singleton stratum with a
certainty FPC is supported under `fail`. Extra fitting weights are
unsupported; place offsets in the model formula.

Predictions and effects use the joint design-based coefficient
covariance. Averaging treats the supplied covariate distribution as
fixed; it does not add design uncertainty from estimating that
distribution. Use
[`suest_newdata()`](https://tdmize.github.io/suest/reference/suest_newdata.md)
and `wts = ".suest_weight"` for model-specific weighted averages.
Inference uses the usual asymptotic normal default. The full design's
degrees of freedom are recorded in `$survey$design_df`; no automatic
survey t or F adjustment is applied.

## Supported models

### Single-level models

- [`stats::lm()`](https://rdrr.io/r/stats/lm.html)

- binary logit, probit, and complementary-log-log models from
  [`stats::glm()`](https://rdrr.io/r/stats/glm.html) or
  [`glm2::glm2()`](https://rdrr.io/pkg/glm2/man/glm2.html)

- ordered logit and probit models from
  [`MASS::polr()`](https://rdrr.io/pkg/MASS/man/polr.html)

- ordered logit and probit models from
  [`ordinal::clm()`](https://rdrr.io/pkg/ordinal/man/clm.html) with
  flexible thresholds, proportional effects, and no scale model

- multinomial logit models from
  [`nnet::multinom()`](https://rdrr.io/pkg/nnet/man/multinom.html)

- Poisson log-link models from
  [`stats::glm()`](https://rdrr.io/r/stats/glm.html) or
  [`glm2::glm2()`](https://rdrr.io/pkg/glm2/man/glm2.html)

- negative-binomial log-link models from
  [`MASS::glm.nb()`](https://rdrr.io/pkg/MASS/man/glm.nb.html)

- Poisson and negative-binomial zero-inflated models from
  [`pscl::zeroinfl()`](https://rdrr.io/pkg/pscl/man/zeroinfl.html)

- other GLMs using identity, log, logit, probit, complementary-log-log,
  or log-log links

- left-, right-, and two-limit censored Gaussian regressions from
  [`censReg::censReg()`](https://rdrr.io/pkg/censReg/man/censReg.html);
  response-scale predictions are the latent mean

- parametric survival and censored-regression models from
  [`survival::survreg()`](https://rdrr.io/pkg/survival/man/survreg.html)
  with a common scale, including Gaussian interval regression

- truncated Gaussian regressions from
  [`truncreg::truncreg()`](https://rdrr.io/pkg/truncreg/man/truncreg.html)

- fractional-response GLMs using
  [`quasibinomial()`](https://rdrr.io/r/stats/family.html) with logit,
  probit, complementary-log-log, or a user-supplied log-log link

- beta regressions from
  [`betareg::betareg()`](https://rdrr.io/pkg/betareg/man/betareg.html)
  using logit, probit, complementary-log-log, or log-log mean links

- heteroskedastic binary probit and logit from
  [`Rchoice::hetprob()`](https://rdrr.io/pkg/Rchoice/man/hetprob.html)

- bivariate probit from
  [`mvProbit::mvProbit()`](https://rdrr.io/pkg/mvProbit/man/mvProbit.html)
  when both equations use the same regressors; fit with `intGrad = TRUE`
  and `finalHessian = TRUE`

- unweighted two-stage least squares from
  [`fixest::feols()`](https://lrberge.github.io/fixest/reference/feols.html)
  without absorbed fixed effects; the original data object must remain
  available

- maximum-likelihood instrumental-variable probit from
  [`Rchoice::ivpml()`](https://rdrr.io/pkg/Rchoice/man/ivpml.html);
  response predictions use the average structural probability

### Selection models

- unweighted maximum-likelihood Heckman selection models from
  `sampleSelection::selection(method = "ml")`, Stata's `heckman`. The
  ancillary parameters use Stata's scale, `lnsigma` and `athrho`, and
  the equations are labeled `selection:` and `outcome:`. Predictions are
  the outcome equation's linear prediction, Stata's default

### Generalized ordered models

- unweighted [`VGAM::vglm()`](https://rdrr.io/pkg/VGAM/man/vglm.html)
  models with the `cumulative()` family and a logit, probit, or cloglog
  link: generalized ordered logit (Stata's `gologit2`), including the
  proportional (`parallel = TRUE`) and partial proportional odds forms
  (for example `parallel = FALSE ~ x`). Use `reverse = TRUE` for
  gologit2's `P(Y > j) = F(eta_j)` parameterization; both directions are
  supported. suest computes the scores analytically and the observed
  information (as Stata does) by differentiating them; VGAM's own
  covariance uses expected information. Predictions are category
  probabilities

### Panel models

- unweighted individual fixed-effects, between-effects, and Swamy-Arora
  random-effects linear panel models from
  [`plm::plm()`](https://rdrr.io/pkg/plm/man/plm.html). Within-model
  prediction uncertainty conditions on estimation-sample means and
  propagates slope uncertainty. At matching evaluation means, the
  covariance is structurally degenerate: average-level confidence
  intervals and hypothesis tests are unsupported, and numerical SEs can
  be missing or nearly zero. Slope and finite-change comparisons remain
  supported.

- unweighted individual random-intercept binary logit and probit models
  from [`pglm::pglm()`](https://rdrr.io/pkg/pglm/man/pglm.html) fitted
  with `model = "random"`, `effect = "individual"`, and `R = 12`;
  response predictions integrate over the random effect

- unweighted gamma random-effects Poisson log models from
  [`pglm::pglm()`](https://rdrr.io/pkg/pglm/man/pglm.html) fitted with
  `model = "random"`, `effect = "individual"`, and `other = "sd"`; the
  final parameter is exposed as gamma variance `alpha`

- unweighted fixed-effects Poisson models from
  [`fixest::fepois()`](https://lrberge.github.io/fixest/reference/feglm.html)
  with one absorbed fixed effect (Stata's `xtpoisson, fe`), clustered on
  that effect with no `G/(G-1)` factor, as in Stata's native
  `xtpoisson, fe vce(robust)`. Each unit's fixed effect is a function of
  the coefficients, `log(sum(y_i)) - log(sum(exp(x_i b)))`, so response
  predictions `exp(x b + alpha_i)` match
  [`predict()`](https://rdrr.io/r/stats/predict.html) for fixest and
  carry the coefficients' uncertainty; newdata must contain the
  fixed-effect variable. Link-scale effects equal the coefficients, as
  in Stata. These models combine only with other fixed-effects Poisson
  models

- unweighted GEE from
  [`geepack::geeglm()`](https://rdrr.io/pkg/geepack/man/geeglm.html):
  Gaussian identity, binary logit/probit/cloglog, and Poisson log, with
  independence or exchangeable correlation and numeric outcomes

### Multilevel models

- unweighted single-level random-intercept Gaussian models from
  [`nlme::lme()`](https://rdrr.io/pkg/nlme/man/lme.html) fitted with
  `method = "ML"`, matching Stata's `xtreg, mle` (full
  observed-information bread). For Stata's `mixed` layout, fit the same
  model with `lme4::lmer(..., REML = FALSE)`; the standard errors differ
  only slightly

- unweighted binomial (logit, probit, cloglog), Poisson-log,
  negative-binomial NB2 log, Gaussian-identity, and Gamma-log models
  from
  [`glmmTMB::glmmTMB()`](https://rdrr.io/pkg/glmmTMB/man/glmmTMB.html)
  with one grouping variable and one conditional random intercept; the
  final parameter is the log random-intercept standard deviation, and
  response predictions integrate over the Gaussian random effect (probit
  `pnorm(eta/sqrt(1 + v))`, log links `exp(eta + v/2)`, Gaussian `eta`,
  logit and cloglog by quadrature, where `v` is the random-effect
  variance). Gaussian models add `log_sigma_e`, the log residual
  standard deviation, and Gamma models `log_shape`, the log shape
  parameter (Stata's `/logs` equals `-log_shape/2`), before the
  random-effect parameters. Gaussian random-effects models from `lmer()`
  or `glmmTMB()` use the covariance layout of Stata's `mixed`: the fixed
  effects' bread is `(X'V^-1 X)^-1`, the variance parameters' bread is
  their block of the full observed-information inverse, and the two
  blocks are uncorrelated. NB2 models must use the default constant
  dispersion model (`dispformula = ~1`); its estimated log size
  parameter `log_phi` precedes `log_sigma`, with conditional variance
  `mu + mu^2/exp(log_phi)`. Weights, offsets, and zero inflation are
  unsupported; NB2 additionally excludes mapped or constrained
  parameters. Its Laplace log likelihood must exceed the
  zero-random-effect NB2 log likelihood at the same fixed effects and
  dispersion by more than
  `sqrt(.Machine$double.eps) * max(1, abs(logLik(model)))`. This
  numerical boundary check is not a significance test and does not
  guarantee an interior global maximum

- the same `glmmTMB` families with one correlated random intercept and
  numeric slope, `(1 + x | id)`, using an unstructured covariance
  matrix. The slope must be a single untransformed numeric column with a
  syntactically valid name. The nuisance parameters are
  `log_sd_intercept`, `log_sd_slope`, and `atanh_rho`. NB2 requires the
  default estimated constant dispersion (`dispformula = ~1`), includes
  `log_phi` before these three parameters, and uses the NB2 boundary
  check above; Gaussian and Gamma models include `log_sigma_e` or
  `log_shape` there. Response predictions use the random-intercept
  formulas with
  `v = var_intercept + 2*x*cov_intercept_slope + x^2*var_slope`;
  binomial responses must be Bernoulli. Link predictions are `X beta`.
  Near-singular random covariance is rejected in a centered,
  standardized predictor basis. Weights, offsets, zero inflation,
  constraints, diagonal covariance, multiple slopes, and slope-only
  terms are unsupported. The slope variable must be supplied for
  response predictions even when it is absent from the fixed-effects
  formula. Native model covariance is used as supplied. Poorly scaled
  predictors can produce inaccurate native numerical curvature even with
  convergence and a positive-definite Hessian; center and rescale
  continuous predictors before fitting, then compare predictions and
  effects after converting them to the same original units. For logit
  random-slope
  [`avg_slopes()`](https://rdrr.io/pkg/marginaleffects/man/slopes.html),
  use `numderiv = list("fdcenter", eps = 1e-4)` and compare results at
  nearby steps (for example, `5e-5` and `2e-4`). An unstable standard
  error requires further investigation; changing the finite-difference
  step does not repair an inaccurate native model covariance

- unweighted [`lme4::lmer()`](https://rdrr.io/pkg/lme4/man/lmer.html)
  models fit with `REML = FALSE`, and
  [`lme4::glmer()`](https://rdrr.io/pkg/lme4/man/glmer.html) binomial
  (logit, probit, cloglog) and Poisson-log models, with one grouping
  variable and a random intercept or one correlated numeric slope. They
  use the same parameters and predictions as the matching `glmmTMB`
  models. suest computes their group scores and observed information
  from its own per-group likelihood: exact for `lmer`, Laplace or
  adaptive Gauss-Hermite with the fit's `nAGQ` for `glmer`.
  `glmer(..., nAGQ = 7)` uses mode-curvature adaptive quadrature with 7
  points, Stata's `intmethod(mcaghermite) intpoints(7)`. lme4's default
  Laplace fits stop the inner mode search early; for closer agreement
  with other software, fit with
  `control = lme4::glmerControl(tolPwrss = 1e-12)`. Likewise, `lmer`'s
  default optimizer tolerance can leave random-slope variance parameters
  slightly short of the optimum;
  `control = lme4::lmerControl(optimizer = "bobyqa", optCtrl = list(rhoend = 1e-12))`
  gives closer agreement with Stata's `mixed`

- unweighted random-intercept ordered logit and probit models from
  [`ordinal::clmm()`](https://rdrr.io/pkg/ordinal/man/clmm.html) with
  one grouping variable and flexible thresholds (Stata's `meologit`,
  `meoprobit`, `xtologit`, `xtoprobit`), fit with Laplace (`nAGQ = 1`)
  or adaptive quadrature (`nAGQ > 1`; 7 points match Stata's
  `intmethod(mcaghermite) intpoints(7)`). Parameters are the thresholds,
  the coefficients, and `log_sigma`. suest computes scores and observed
  information from its own per-group likelihood, which reproduces
  clmm's. Predictions are category probabilities integrated over the
  random intercept

### Combining random-effects and ordinary models

Random-effects models from
[`nlme::lme()`](https://rdrr.io/pkg/nlme/man/lme.html),
[`pglm::pglm()`](https://rdrr.io/pkg/pglm/man/pglm.html),
[`glmmTMB::glmmTMB()`](https://rdrr.io/pkg/glmmTMB/man/glmmTMB.html),
lme4, and [`ordinal::clmm()`](https://rdrr.io/pkg/ordinal/man/clmm.html)
can be combined with each other, including different families, and with
linear, binary logit, probit, and cloglog, Poisson, negative binomial,
and ordered logit and probit models, as in Stata's suest2. Every
random-effects model must use the same grouping variable, and the
ordinary models' data must contain it. The system is clustered on that
group, or on `cluster`, which must contain whole groups. Each model's
block of the covariance then equals that model's own covariance
clustered on the group, including the `(N - 1)/(N - k)` adjustment for
linear models. `plm` and `geeglm` models combine only with models of the
same type.

### Survey models

- restricted survey-weighted Gaussian identity, binary
  [`quasibinomial()`](https://rdrr.io/r/stats/family.html)
  logit/probit/cloglog, and
  [`quasipoisson()`](https://rdrr.io/r/stats/family.html) log models
  from [`survey::svyglm()`](https://rdrr.io/pkg/survey/man/svyglm.html),
  in any combination

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
