# Combine two fitted models with seemingly unrelated estimation

`suest()` combines two separately fitted models and constructs a joint
model-robust covariance matrix from their observation-level score
contributions. The returned object can be passed directly to
[`marginaleffects::predictions()`](https://rdrr.io/pkg/marginaleffects/man/predictions.html),
[`marginaleffects::avg_comparisons()`](https://rdrr.io/pkg/marginaleffects/man/comparisons.html),
[`marginaleffects::avg_slopes()`](https://rdrr.io/pkg/marginaleffects/man/slopes.html),
and
[`marginaleffects::hypotheses()`](https://rdrr.io/pkg/marginaleffects/man/hypotheses.html).

## Usage

``` r
suest(model1, model2, model_names = NULL)

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

- model1, model2:

  Two supported fitted model objects.

- model_names:

  Optional character vector containing two display names. By default,
  the object names supplied in the call are used.

- object, x:

  A `"suest_model"` object.

- ...:

  Additional arguments, currently ignored.

## Value

An object of class `"suest_model"` containing the two fitted models,
their joint coefficient vector, a joint model-robust covariance matrix,
and sample-alignment information.

## Details

Exactly two models are supported. Models can use identical, partially
overlapping, or completely disjoint samples. When the model calls refer
to the same data source, model-frame row names identify overlapping
observations. Models fitted from different data objects are treated as
disjoint because shared observations cannot be inferred safely without
an identifier.

Same-family comparisons are supported for all model types listed below.
The supported cross-family pairs are logit–probit, logit–linear,
probit–linear, Poisson–negative binomial, and ordered logit–multinomial
logit.

Negative-binomial models include `log(theta)` in the joint parameter
vector. Ordered and multinomial models use analytic score and
observed-information calculations for stable robust covariance
estimation.

Offsets, nonunit estimation weights, and aliased parameters are not
currently supported. Bias-reduced, adjusted-score, Firth, and penalized
GLM fits are rejected because they do not use ordinary
maximum-likelihood score equations.

## Supported models

- [`stats::lm()`](https://rdrr.io/r/stats/lm.html)

- binary logit and probit models from
  [`stats::glm()`](https://rdrr.io/r/stats/glm.html) or
  [`glm2::glm2()`](https://rdrr.io/pkg/glm2/man/glm2.html)

- Poisson log-link models from
  [`stats::glm()`](https://rdrr.io/r/stats/glm.html) or
  [`glm2::glm2()`](https://rdrr.io/pkg/glm2/man/glm2.html)

- negative-binomial log-link models from
  [`MASS::glm.nb()`](https://rdrr.io/pkg/MASS/man/glm.nb.html)

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
#>           Hypothesis Estimate Std. Error     z Pr(>|z|)    S  2.5 %  97.5 %
#>  (Adjusted) - (Base)  -0.0779     0.0141 -5.54   <0.001 25.0 -0.105 -0.0503
#> 
#> 
```
