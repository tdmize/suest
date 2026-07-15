# Getting started with suest

[`suest()`](https://tdmize.github.io/suest/reference/suest.md) combines
two separately fitted models and supplies the joint robust covariance
needed to test whether predictions or effects differ across models. The
package leaves the calculation and presentation of substantive
quantities to `marginaleffects`.

## Installation and setup

``` r

library(suest)
library(marginaleffects)
```

## Fit two models

``` r

dat <- mtcars
dat$am <- factor(dat$am)

base <- glm(am ~ wt, family = binomial(), data = dat)
adjusted <- glm(am ~ wt + hp, family = binomial(), data = dat)
```

## Combine the models

``` r

combined <- suest(
  base,
  adjusted,
  model_names = c("Base", "Adjusted")
)

combined
```

The printed object reports the component model types, the comparison
scale, the component sample sizes, the number of overlapping
observations, and the number of parameters in the joint system.

## Compute and compare effects

``` r

effects <- avg_comparisons(
  combined,
  variables = "wt",
  newdata = dat
)

effects

hypotheses(
  effects,
  hypothesis = difference ~ revpairwise
)
```

The final command directly tests the difference between the effects. It
does not compare whether one effect is statistically significant while
the other is not.

## Different estimation samples

When models were fitted on different samples, use
[`suest_newdata()`](https://tdmize.github.io/suest/reference/suest_newdata.md)
to average each model’s effect over its own estimation sample.

``` r

model_4cyl <- glm(
  am ~ wt + hp,
  family = binomial(),
  data = dat,
  subset = cyl == 4
)

model_other <- glm(
  am ~ wt + hp,
  family = binomial(),
  data = dat,
  subset = cyl != 4
)

combined_samples <- suest(
  model_4cyl,
  model_other,
  model_names = c("Four cylinders", "Other")
)

nd <- suest_newdata(combined_samples)

avg_comparisons(
  combined_samples,
  variables = "wt",
  newdata = nd,
  hypothesis = difference ~ revpairwise
)
```

For disjoint samples, the cross-model covariance is zero, so the
variance of the difference is the sum of the two component variances.
