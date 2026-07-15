# Stack the estimation samples from a SUEST object

Creates a data frame containing each component model's own estimation
sample. This is useful when the two models were fitted on different
samples and marginal effects should be averaged separately within each
model's observed covariate distribution.

## Usage

``` r
suest_newdata(object)
```

## Arguments

- object:

  A `"suest_model"` returned by
  [`suest()`](https://tdmize.github.io/suest/reference/suest.md).

## Value

A data frame with the component model frames stacked vertically and two
internal columns, `.suest_model` and `.suest_rowid`, used to route rows
to the correct component model.

## Examples

``` r
dat <- mtcars
dat$am <- factor(dat$am)

model1 <- glm(am ~ wt + hp, family = binomial(), data = dat, subset = cyl == 4)
model2 <- glm(am ~ wt + hp, family = binomial(), data = dat, subset = cyl != 4)
#> Warning: glm.fit: algorithm did not converge
#> Warning: glm.fit: fitted probabilities numerically 0 or 1 occurred
fit <- suest(model1, model2, model_names = c("Four cylinders", "Other"))
nd <- suest_newdata(fit)
marginaleffects::avg_comparisons(fit, variables = "wt", newdata = nd)
#> 
#>           Group Estimate Std. Error         z Pr(>|z|)    S  2.5 % 97.5 %
#>  Four cylinders   -0.596   9.08e-02 -6.57e+00   <0.001 34.2 -0.774 -0.418
#>  Other            -0.238   2.90e-11 -8.20e+09   <0.001  Inf -0.238 -0.238
#> 
#> Term: wt
#> Type: response
#> Comparison: +1
#> 
```
