# Two-part models

`suest_twopart()` combines the two parts of a two-part model (Stata's
`twopm`): a binary model for whether the outcome is positive, fit to
everyone, and a model for the outcome among those with a positive
outcome. Predictions, slopes, and comparisons from `marginaleffects` are
then on the overall expected outcome, `Pr(y > 0) * E(y | y > 0)`, with
standard errors from the joint covariance matrix of the two parts.

## Usage

``` r
suest_twopart(first, second, model_names = NULL, ...)

# S3 method for class 'suest_twopart'
print(x, ...)
```

## Arguments

- first:

  Binary logit, probit, or complementary log-log model for whether the
  outcome is positive.

- second:

  Model for the outcome among the observations with a positive outcome:
  a linear model or a GLM, for example a Gamma model or a Poisson
  pseudo-likelihood model
  ([`quasipoisson()`](https://rdrr.io/r/stats/family.html)) with a log
  link.

- model_names:

  Names of the two parts. By default, the object names supplied in the
  call.

- ...:

  Passed to
  [`suest()`](https://tdmize.github.io/suest/reference/suest.md), for
  example `cluster`, `weight_type`, or `survey_design` with
  `observation_id`.

- x:

  A `"suest_twopart"` object.

## Value

An object of class `c("suest_twopart", "suest_model")`.

## Details

Fit the second part to every observation in the first part's sample that
has a positive outcome, for example with `subset = y > 0`. Without
`newdata`, results are averaged over the first part's estimation sample.
`type = "response"` (the default) gives the overall expected outcome,
`type = "prob"` the probability of a positive outcome (the first part),
and `type = "conditional"` the expected outcome among those with a
positive outcome (the second part).

## Examples

``` r
dat <- mtcars
dat$y <- ifelse(dat$am == 1, dat$mpg, 0)
any <- glm(I(y > 0) ~ wt, family = binomial(), data = dat)
amount <- glm(y ~ wt, family = Gamma(link = "log"), data = dat, subset = y > 0)
fit <- suest_twopart(any, amount)
marginaleffects::avg_slopes(fit, variables = "wt")
#> 
#>  Estimate Std. Error     z Pr(>|z|)    S 2.5 % 97.5 %
#>     -10.7      0.953 -11.2   <0.001 94.5 -12.5  -8.82
#> 
#> Term: wt
#> Type: response
#> Comparison: dY/dX
#> 
```
