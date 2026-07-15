# suest for R

[![R-CMD-check](https://github.com/tdmize/suest/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/tdmize/suest/actions/workflows/R-CMD-check.yaml)

Full documentation, worked examples, and rendered output:
[tdmize.github.io/suest](https://tdmize.github.io/suest/articles/suest.html)

`suest` combines two separately fitted regression models into one object with a
joint model-robust covariance matrix. The combined object works with
[`marginaleffects`](https://marginaleffects.com/) to compare predictions and
marginal effects across models.

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

base <- glm(am ~ wt, family = binomial("logit"), data = dat)
adjusted <- glm(am ~ wt + hp, family = binomial("logit"), data = dat)
combined <- suest(base, adjusted, model_names = c("Base", "Adjusted"))

effects <- avg_comparisons(combined, variables = "wt", newdata = dat)
effects
hypotheses(effects, hypothesis = difference ~ revpairwise)
```

`suest()` accounts for the cross-model covariance between estimates. This is
essential when the models use the same or overlapping observations.

## Supported models

- linear regression using `lm()`
- binary logit and probit using `glm()` or `glm2::glm2()`
- Poisson regression using `glm()` or `glm2::glm2()`
- negative-binomial regression using `MASS::glm.nb()`
- ordered logit and probit using `MASS::polr()`
- ordered logit and probit using restricted `ordinal::clm()` specifications
- multinomial logit using `nnet::multinom()`

Models may use identical, partially overlapping, or disjoint samples. Supported
cross-family comparisons are logit–probit, logit–linear, probit–linear,
Poisson–negative binomial, and ordered logit–multinomial logit.

## Worked examples

The single Get Started page contains code and rendered results for six R
replications corresponding to Examples 6.1–6.6 for the
[`mecompare` Stata command](https://www.trentonmize.com/software/mecompare):

1. Marginal effects to summarize curvilinear relationships and test mediation
2. Comparing marginal effects across nested logit models
3. Comparing marginal effects using alternative predictors
4. Comparing marginal effects across different outcomes
5. Comparing marginal effects across different model types (ordinal vs nominal)
6. Comparing marginal effects across different samples or groups

It also documents cross-engine comparisons using `glm2::glm2()` and
`ordinal::clm()`.

Bias-reduced, adjusted-score, Firth, and penalized GLM fits, including
`brglm2::brglmFit()`, are rejected rather than being treated incorrectly as
ordinary maximum-likelihood GLMs.
