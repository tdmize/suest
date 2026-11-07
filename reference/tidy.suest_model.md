# Tidy and glance methods for regression tables

[`tidy()`](https://generics.r-lib.org/reference/tidy.html) returns one
row per parameter, with the model name in a `group` column, so
`modelsummary::modelsummary(fit, shape = term ~ group)` puts each model
in its own column. Standard errors come from the joint covariance
matrix. [`glance()`](https://generics.r-lib.org/reference/glance.html)
returns the number of observations in the combined sample.

## Usage

``` r
# S3 method for class 'suest_model'
tidy(x, conf.int = TRUE, conf.level = 0.95, ...)

# S3 method for class 'suest_mi'
tidy(x, conf.int = TRUE, conf.level = 0.95, ...)

# S3 method for class 'suest_model'
glance(x, ...)

# S3 method for class 'suest_mi'
glance(x, ...)
```

## Arguments

- x:

  A `"suest_model"` or `"suest_mi"` object.

- conf.int:

  Include confidence intervals.

- conf.level:

  Confidence level.

- ...:

  Not used.

## Value

A data frame.

## Examples

``` r
dat <- mtcars
model1 <- glm(am ~ wt, family = binomial(), data = dat)
model2 <- glm(vs ~ wt, family = binomial(), data = dat)
fit <- suest(model1, model2)
if (requireNamespace("generics", quietly = TRUE)) generics::tidy(fit)
#>          term  group  estimate std.error statistic     p.value  conf.low
#> 1 (Intercept) model1 12.040370 5.0352099  2.391235 0.016791806  2.171540
#> 2          wt model1 -4.023970 1.5826157 -2.542607 0.011002885 -7.125840
#> 3 (Intercept) model2  5.714682 2.0159983  2.834666 0.004587364  1.763397
#> 4          wt model2 -1.910548 0.6093721 -3.135273 0.001716945 -3.104895
#>    conf.high
#> 1 21.9091998
#> 2 -0.9221002
#> 3  9.6659656
#> 4 -0.7162003
```
