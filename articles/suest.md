# Getting started with suest

[`suest()`](https://tdmize.github.io/suest/reference/suest.md) combines
two or more separately fitted models so you can test whether
predictions, marginal effects, or other quantities differ across the
models. It supplies the joint robust covariance matrix those tests need,
and leaves calculating and presenting the quantities themselves to
[`marginaleffects`](https://marginaleffects.com/).

## Functions at a glance

- [`suest()`](https://tdmize.github.io/suest/reference/suest.md)
  combines two or more fitted models and calculates their joint
  covariance matrix.
- [`suest_newdata()`](https://tdmize.github.io/suest/reference/suest_newdata.md)
  stacks the models’ estimation samples so each model’s effects can be
  averaged over its own sample.
- [`suest_mi()`](https://tdmize.github.io/suest/reference/suest_mi.md)
  pools [`suest()`](https://tdmize.github.io/suest/reference/suest.md)
  results across multiply imputed datasets.
- The usual `marginaleffects` functions, including
  [`predictions()`](https://rdrr.io/pkg/marginaleffects/man/predictions.html),
  [`avg_comparisons()`](https://rdrr.io/pkg/marginaleffects/man/comparisons.html),
  [`avg_slopes()`](https://rdrr.io/pkg/marginaleffects/man/slopes.html),
  [`hypotheses()`](https://rdrr.io/pkg/marginaleffects/man/hypotheses.html),
  and the plotting functions, work with the combined object.

## Installation and setup

``` r

remotes::install_github("tdmize/suest")
```

``` r

library(suest)
library(marginaleffects)
```

## Basic workflow

``` r

dat <- mtcars
dat$am <- factor(dat$am)

base <- glm(am ~ wt, family = binomial("logit"), data = dat)
adjusted <- glm(am ~ wt + hp, family = binomial("logit"), data = dat)
combined <- suest(base, adjusted, model_names = c("Base", "Adjusted"))
combined
#> Seemingly Unrelated Estimation
#> Models: Base + Adjusted 
#> Model types: Base=logit, Adjusted=logit 
#> Model engines: Base=stats::glm, Adjusted=stats::glm 
#> Comparison scale: predicted probabilities 
#> Observations: Base=32, Adjusted=32 
#> Overlapping observations: 32 
#> Union observations: 32 
#> Parameters: 5

effects <- avg_comparisons(combined, variables = "wt", newdata = dat)
effects
#> 
#>     Group Estimate Std. Error     z Pr(>|z|)     S  2.5 % 97.5 %
#>  Base       -0.277     0.0232 -11.9   <0.001 106.5 -0.322 -0.231
#>  Adjusted   -0.355     0.0184 -19.3   <0.001 273.3 -0.391 -0.319
#> 
#> Term: wt
#> Type: response
#> Comparison: +1
hypotheses(effects, hypothesis = difference ~ revpairwise)
#> 
#>           Hypothesis Estimate Std. Error    z Pr(>|z|)    S  2.5 % 97.5 %
#>  (Base) - (Adjusted)   0.0779     0.0141 5.54   <0.001 25.0 0.0503  0.105
```

The last command tests the difference between the two effects directly.
That is different from checking whether one effect is statistically
significant and the other is not.

## Supported models

### Single-level models

| Model | Fitted with |
|----|----|
| Linear regression | [`lm()`](https://rdrr.io/r/stats/lm.html) |
| Binary logit, probit, and complementary log-log | [`glm()`](https://rdrr.io/r/stats/glm.html), [`glm2::glm2()`](https://rdrr.io/pkg/glm2/man/glm2.html) |
| Ordered logit and probit | [`MASS::polr()`](https://rdrr.io/pkg/MASS/man/polr.html), [`ordinal::clm()`](https://rdrr.io/pkg/ordinal/man/clm.html) |
| Multinomial logit | [`nnet::multinom()`](https://rdrr.io/pkg/nnet/man/multinom.html) |
| Poisson regression | [`glm()`](https://rdrr.io/r/stats/glm.html), [`glm2::glm2()`](https://rdrr.io/pkg/glm2/man/glm2.html) |
| Negative binomial regression | [`MASS::glm.nb()`](https://rdrr.io/pkg/MASS/man/glm.nb.html) |
| Zero-inflated Poisson and negative binomial | [`pscl::zeroinfl()`](https://rdrr.io/pkg/pscl/man/zeroinfl.html) |
| Other GLMs with identity, log, logit, probit, complementary log-log, or log-log links | [`glm()`](https://rdrr.io/r/stats/glm.html), [`glm2::glm2()`](https://rdrr.io/pkg/glm2/man/glm2.html) |
| Tobit regression | [`censReg::censReg()`](https://rdrr.io/pkg/censReg/man/censReg.html); predictions are the latent mean |
| Interval, censored, and parametric survival regression | [`survival::survreg()`](https://rdrr.io/pkg/survival/man/survreg.html) with a common scale |
| Truncated regression | [`truncreg::truncreg()`](https://rdrr.io/pkg/truncreg/man/truncreg.html) |
| Fractional response | [`glm()`](https://rdrr.io/r/stats/glm.html) with [`quasibinomial()`](https://rdrr.io/r/stats/family.html): logit, probit, complementary log-log, or log-log link |
| Beta regression | [`betareg::betareg()`](https://rdrr.io/pkg/betareg/man/betareg.html) with a logit, probit, complementary log-log, or log-log link |
| Heteroskedastic probit and logit | [`Rchoice::hetprob()`](https://rdrr.io/pkg/Rchoice/man/hetprob.html) |
| Bivariate probit | [`mvProbit::mvProbit()`](https://rdrr.io/pkg/mvProbit/man/mvProbit.html) with the same predictors in both equations |
| Instrumental-variable regression (two-stage least squares) | [`fixest::feols()`](https://lrberge.github.io/fixest/reference/feols.html), unweighted and without absorbed fixed effects |
| Instrumental-variable probit | [`Rchoice::ivpml()`](https://rdrr.io/pkg/Rchoice/man/ivpml.html); predictions are average structural probabilities |
| Generalized ordered logit (gologit2), including partial proportional odds | [`VGAM::vglm()`](https://rdrr.io/pkg/VGAM/man/vglm.html) with `cumulative(reverse = TRUE)` |
| Heckman selection (maximum likelihood) | `sampleSelection::selection(method = "ml")` |

### Panel models

| Model | Fitted with |
|----|----|
| Linear fixed-, between-, and random-effects | [`plm::plm()`](https://rdrr.io/pkg/plm/man/plm.html) with individual effects |
| Random-effects logit and probit | `pglm::pglm(model = "random", effect = "individual")` |
| Random-effects Poisson | `pglm::pglm(model = "random", effect = "individual", other = "sd")` |
| Fixed-effects Poisson | `fixest::fepois(y ~ x | id)` |
| Population-averaged (GEE) linear, logit, probit, complementary log-log, and Poisson | [`geepack::geeglm()`](https://rdrr.io/pkg/geepack/man/geeglm.html) with independence or exchangeable correlation |

### Multilevel models

| Model | Fitted with |
|----|----|
| Linear random-intercept | `nlme::lme(method = "ML")` |
| Random-intercept logit, probit, cloglog, Poisson, negative binomial, linear, and Gamma (log link) | [`glmmTMB::glmmTMB()`](https://rdrr.io/pkg/glmmTMB/man/glmmTMB.html) |
| Random intercept and slope, `(1 + x | id)`, for the same families | [`glmmTMB::glmmTMB()`](https://rdrr.io/pkg/glmmTMB/man/glmmTMB.html) |
| Linear random intercept or slope (ML) | `lme4::lmer(REML = FALSE)` |
| Random-intercept or slope logit, probit, cloglog, and Poisson | [`lme4::glmer()`](https://rdrr.io/pkg/lme4/man/glmer.html), any `nAGQ` |
| Random-intercept ordered logit and probit | [`ordinal::clmm()`](https://rdrr.io/pkg/ordinal/man/clmm.html) |

### Survey models

| Model | Fitted with |
|----|----|
| Linear regression | [`survey::svyglm()`](https://rdrr.io/pkg/survey/man/svyglm.html) |
| Binary logit, probit, and cloglog | [`survey::svyglm()`](https://rdrr.io/pkg/survey/man/svyglm.html) with [`quasibinomial()`](https://rdrr.io/r/stats/family.html) |
| Poisson | [`survey::svyglm()`](https://rdrr.io/pkg/survey/man/svyglm.html) with [`quasipoisson()`](https://rdrr.io/r/stats/family.html) |

See [Survey data](#survey-data) below.

### Combining models

Models can use the same, partially overlapping, or completely different
samples. Any single-level models can be combined, including models for
different outcomes and models for categorical outcomes. Results on
different scales are labeled separately rather than treated as directly
comparable. Random-effects models
([`nlme::lme()`](https://rdrr.io/pkg/nlme/man/lme.html), `pglm()`,
`glmmTMB()`, lme4,
[`clmm()`](https://rdrr.io/pkg/ordinal/man/clmm.html)) can be combined
with each other, including different families, and with linear, binary,
count, and ordered models, such as a logit model with a random-intercept
logit model. They must share one grouping variable, which the system is
clustered on, and must be unweighted. Other panel models (`plm()`,
`geeglm()`) can be combined only with models of the same kind. Offsets
are supported, except in `pglm()` and `glmmTMB()` models.

A few models have settings worth knowing about; see [Notes on specific
models](https://tdmize.github.io/suest/articles/model-notes.md).

## Different samples or groups

When the models are fitted to different samples, such as separate
groups, use
[`suest_newdata()`](https://tdmize.github.io/suest/reference/suest_newdata.md)
so each model’s effects are averaged over its own sample:

``` r

four <- lm(mpg ~ wt + hp, data = mtcars, subset = cyl == 4)
other <- lm(mpg ~ wt + hp, data = mtcars, subset = cyl != 4)
fit_groups <- suest(four, other, model_names = c("Four cylinders", "Other"))

effects_groups <- avg_comparisons(fit_groups, variables = "wt",
                                  newdata = suest_newdata(fit_groups))
effects_groups
#> 
#>           Group Estimate Std. Error     z Pr(>|z|)    S 2.5 % 97.5 %
#>  Four cylinders    -5.12      1.083 -4.72   <0.001 18.7 -7.24  -2.99
#>  Other             -2.56      0.575 -4.45   <0.001 16.8 -3.69  -1.43
#> 
#> Term: wt
#> Type: response
#> Comparison: +1
hypotheses(effects_groups, hypothesis = difference ~ revpairwise)
#> 
#>                  Hypothesis Estimate Std. Error     z Pr(>|z|)   S 2.5 % 97.5 %
#>  (Four cylinders) - (Other)    -2.55       1.23 -2.08   0.0373 4.7 -4.96  -0.15
```

## Cluster-robust standard errors

Supply a cluster variable to get cluster-robust standard errors for the
combined models:

``` r

fit_clustered <- suest(base, adjusted, cluster = "person_id")
avg_comparisons(fit_clustered, variables = "x", newdata = dat)
```

`cluster` takes one or more column names, or a list with one vector,
matrix, or data frame matched to each model’s estimation sample. An
observation that appears in more than one model must have the same
cluster ID in each. Observations in completely different samples can
belong to the same cluster, in which case the covariance between models
can be nonzero.

Panel models cluster on the panel identifier by default. A cluster you
supply must contain whole panels.

## Probability weights

`suest` supports probability weights for linear, binary logit and
probit, Poisson, negative binomial, ordered logit and probit, and
multinomial logit models. Fit each model with `weights=`, and tell
[`suest()`](https://tdmize.github.io/suest/reference/suest.md) that they
are probability weights with `weight_type = "pweight"`:

``` r

set.seed(376)
weighted <- data.frame(x = rnorm(500), z = rnorm(500))
weighted$mediator <- 0.5 * weighted$x - 0.3 * weighted$z + rnorm(500)
weighted$pw <- exp(0.2 * weighted$z + 0.1 * rnorm(500))
weighted$y <- 1 + 0.7 * weighted$x - 0.4 * weighted$z +
  0.6 * weighted$mediator + rnorm(500)

base_w <- lm(y ~ x + z, weights = pw, data = weighted)
adjusted_w <- lm(y ~ x + z + mediator, weights = pw, data = weighted)
fit_w <- suest(base_w, adjusted_w, model_names = c("Base", "Adjusted"),
               weight_type = "pweight")
fit_w
#> Seemingly Unrelated Estimation
#> Models: Base + Adjusted 
#> Model types: Base=lm, Adjusted=lm 
#> Model engines: Base=stats::lm, Adjusted=stats::lm 
#> Comparison scale: fitted values 
#> Weight type: pweight 
#> Observations: Base=500, Adjusted=500 
#> Overlapping observations: 500 
#> Union observations: 500 
#> Parameters: 9

effects_w <- avg_comparisons(fit_w, variables = "x", newdata = weighted, wts = "pw")
effects_w
#> 
#>     Group Estimate Std. Error    z Pr(>|z|)     S 2.5 % 97.5 %
#>  Base        1.003     0.0537 18.7   <0.001 255.9 0.898  1.108
#>  Adjusted    0.722     0.0480 15.0   <0.001 167.3 0.628  0.816
#> 
#> Term: x
#> Type: response
#> Comparison: +1
hypotheses(effects_w, hypothesis = difference ~ revpairwise)
#> 
#>           Hypothesis Estimate Std. Error    z Pr(>|z|)    S 2.5 % 97.5 %
#>  (Base) - (Adjusted)    0.281     0.0322 8.74   <0.001 58.6 0.218  0.345
```

Weights must be positive and finite. An observation that appears in more
than one model must have the same weight in each; weights can differ for
observations in only one model, including completely different samples.
When the models have different estimation samples,
[`suest_newdata()`](https://tdmize.github.io/suest/reference/suest_newdata.md)
includes the weights as `.suest_weight`:

``` r

nd_w <- suest_newdata(fit_w)
avg_comparisons(fit_w, variables = "x", newdata = nd_w, wts = ".suest_weight")
```

## Survey data

`suest` supports linear, binary (logit, probit, cloglog), and Poisson
models fitted with
[`survey::svyglm()`](https://rdrr.io/pkg/survey/man/svyglm.html) from
one common one-stage survey design, in any combination:

``` r

library(survey)
design <- svydesign(~psu, strata = ~strata, weights = ~w, data = dat)
m1 <- svyglm(yb ~ x, design, family = quasibinomial("logit"))
m2 <- svyglm(yb ~ x + z, design, family = quasibinomial("logit"))
fit <- suest(m1, m2, survey_design = design, observation_id = "id")
avg_slopes(fit, variables = "x", newdata = suest_newdata(fit), wts = ".suest_weight")
```

If either model uses a subset (domain) of the data, keep `design` as the
full sampled design. The joint covariance accounts for strata, PSUs, and
optional one-stage finite population corrections. Replicate weights,
multistage designs, calibration, and PPS designs are not supported.

## Multiple imputation

To combine results across imputed datasets, fit the same `suest` system
in each one. With `mice`, [`with()`](https://rdrr.io/r/base/with.html)
does this, and `marginaleffects` functions then estimate within each
imputation and pool the results with Rubin’s rules, including
cross-model differences:

``` r

analyses <- with(imp, suest(
  glm(y ~ x, family = binomial),
  glm(y ~ x + z, family = binomial),
  model_names = c("Base", "Adjusted")
))
avg_comparisons(analyses, variables = "x",
  hypothesis = difference ~ revpairwise)
```

`mice` warns “Large sample assumed.” here. That is expected: SUEST
reports large-sample statistics, so the pooled degrees of freedom are
Rubin’s large-sample values.

[`suest_mi()`](https://tdmize.github.io/suest/reference/suest_mi.md)
pools the coefficients and their joint covariance:

``` r

pooled <- suest_mi(analyses)
summary(pooled)
```

`marginaleffects` functions give the same results on `pooled` as on
`analyses`. A list of
[`suest()`](https://tdmize.github.io/suest/reference/suest.md) results,
one per imputation, or the result of
`with(mitools::imputationList(...), suest(...))` can also be passed to
[`suest_mi()`](https://tdmize.github.io/suest/reference/suest_mi.md);
`marginaleffects` then warns that it could not recover the original data
from a `mids` object, which does not affect
[`suest()`](https://tdmize.github.io/suest/reference/suest.md) results.

## Examples

Six examples from Mize, Doan, and Long (2019), with code and output:

1.  [Curvilinear effects and
    mediation/attenuation](https://tdmize.github.io/suest/articles/mediation.md)
2.  [Nested logit models and
    mediation/attenuation](https://tdmize.github.io/suest/articles/nested-models.md)
3.  [Alternative
    predictors](https://tdmize.github.io/suest/articles/alternative-predictors.md)
4.  [Different
    outcomes](https://tdmize.github.io/suest/articles/different-outcomes.md)
5.  [Different model
    types](https://tdmize.github.io/suest/articles/model-types.md):
    ordinal versus nominal
6.  [Separate samples or
    groups](https://tdmize.github.io/suest/articles/groups.md)
