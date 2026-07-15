# suest for R

[![R-CMD-check](https://github.com/tdmize/suest/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/tdmize/suest/actions/workflows/R-CMD-check.yaml)

Full documentation and replication examples: [tdmize.github.io/suest](https://tdmize.github.io/suest/)

`suest` combines two separately fitted regression models into one object with a joint model-robust covariance matrix. The combined object works with [`marginaleffects`](https://marginaleffects.com/) to compare predictions and marginal effects across models.

The package implements the framework developed in:

> Mize, Trenton D., Long Doan, and J. Scott Long. 2019. “A General Framework for Comparing Predictions and Marginal Effects Across Models.” *Sociological Methodology* 49(1):152–189. [doi:10.1177/0081175019852763](https://doi.org/10.1177/0081175019852763)

## Installation

```r
# install.packages("remotes")
remotes::install_github("tdmize/suest")
```

## Basic usage

```r
library(suest)
library(marginaleffects)

dat <- mtcars
dat$am <- factor(dat$am)

base <- glm(am ~ wt, family = binomial(), data = dat)
adjusted <- glm(am ~ wt + hp, family = binomial(), data = dat)

combined <- suest(
  base,
  adjusted,
  model_names = c("Base", "Adjusted")
)

effects <- avg_comparisons(
  combined,
  variables = "wt",
  newdata = dat
)

effects
hypotheses(effects, hypothesis = difference ~ revpairwise)
```

`suest()` accounts for the cross-model covariance between estimates. This is essential when the models use the same or overlapping observations.

## Supported models

- linear regression using `lm()`
- binary logit and probit using `glm()`
- Poisson regression using `glm()`
- negative-binomial regression using `MASS::glm.nb()`
- ordered logit and probit using `MASS::polr()`
- multinomial logit using `nnet::multinom()`

Models may use identical, partially overlapping, or completely disjoint samples. The supported cross-family comparisons are logit–probit, logit–linear, probit–linear, Poisson–negative binomial, and ordered logit–multinomial logit.

## Replication vignettes

The package website includes six R replications corresponding to Examples 6.1–6.6 on the [`mecompare` page](https://www.trentonmize.com/software/mecompare):

1. Curvilinear effects and mediation
2. Nested logit models
3. Alternative predictor operationalizations
4. Different count outcomes
5. Ordered versus multinomial logit
6. Different samples
