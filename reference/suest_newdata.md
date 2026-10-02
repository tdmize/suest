# Stack the estimation samples from a SUEST object

Creates a data frame that stacks each model's own estimation sample. Use
it when the models were fitted on different samples, so that each
model's marginal effects are averaged over its own sample.

## Usage

``` r
suest_newdata(object)
```

## Arguments

- object:

  A `"suest_model"` returned by
  [`suest()`](https://tdmize.github.io/suest/reference/suest.md).

## Value

A data frame with the component model frames stacked vertically and
internal columns `.suest_model` and `.suest_rowid`, used to route rows
to the correct component model. For pweighted and survey models,
`.suest_weight` contains each row's evaluated sampling weight and can be
supplied to the `wts` argument of `marginaleffects` averaging functions.

## Examples

``` r
dat <- mtcars

model1 <- lm(mpg ~ wt + hp, data = dat, subset = cyl == 4)
model2 <- lm(mpg ~ wt + hp, data = dat, subset = cyl != 4)
fit <- suest(model1, model2, model_names = c("Four cylinders", "Other"))
nd <- suest_newdata(fit)
marginaleffects::avg_comparisons(fit, variables = "wt", newdata = nd)
#> 
#>           Group Estimate Std. Error     z Pr(>|z|)    S 2.5 % 97.5 %
#>  Four cylinders    -5.12      1.083 -4.72   <0.001 18.7 -7.24  -2.99
#>  Other             -2.56      0.575 -4.45   <0.001 16.8 -3.69  -1.43
#> 
#> Term: wt
#> Type: response
#> Comparison: +1
#> 
```
