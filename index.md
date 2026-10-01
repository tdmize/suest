# suest for R

[![R-CMD-check](https://github.com/tdmize/suest/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/tdmize/suest/actions/workflows/R-CMD-check.yaml)

`suest` lets you compare predictions and marginal effects across
regression models in R. It combines two or more separately fitted models
into one object with a joint robust covariance matrix. That object works
directly with [`marginaleffects`](https://marginaleffects.com/), so you
can calculate predictions and effects for each model and test whether
they differ across models.

The package implements the framework in:

> Mize, Trenton D., Long Doan, and J. Scott Long. 2019. “A General
> Framework for Comparing Predictions and Marginal Effects Across
> Models.” *Sociological Methodology* 49(1):152–189.
> [doi:10.1177/0081175019852763](https://doi.org/10.1177/0081175019852763)

The same framework is available in Stata through the [`mecompare`
command](https://www.trentonmize.com/software/mecompare).

## Installation

``` r

# install.packages("remotes")
remotes::install_github("tdmize/suest")
```

## Basic usage

``` r

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

[`suest()`](https://tdmize.github.io/suest/reference/suest.md) accounts
for the covariance between estimates from different models. This is
essential when the models use the same or overlapping observations.

## Supported models

`suest` works with linear, binary, ordinal, multinomial, count,
censored, survival, and instrumental-variable models, as well as panel,
multilevel, and survey models. Models can use the same, partially
overlapping, or completely different samples. See the [full list of
supported
models](https://www.trentonmize.com/software/suest_r/articles/suest.html#supported-models).

## Learn more

- [Getting
  started](https://www.trentonmize.com/software/suest_r/articles/suest.html):
  the workflow, supported models, different samples, clustered standard
  errors, probability weights, survey data, and multiple imputation.
- Examples from Mize, Doan, and Long (2019), with code and output:
  [curvilinear effects and
  mediation](https://www.trentonmize.com/software/suest_r/articles/mediation.html),
  [nested logit
  models](https://www.trentonmize.com/software/suest_r/articles/nested-models.html),
  [alternative
  predictors](https://www.trentonmize.com/software/suest_r/articles/alternative-predictors.html),
  [different
  outcomes](https://www.trentonmize.com/software/suest_r/articles/different-outcomes.html),
  [different model
  types](https://www.trentonmize.com/software/suest_r/articles/model-types.html),
  and [separate samples or
  groups](https://www.trentonmize.com/software/suest_r/articles/groups.html).
- [Notes on specific
  models](https://www.trentonmize.com/software/suest_r/articles/model-notes.html).
