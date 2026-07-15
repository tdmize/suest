# Example 6.1: Curvilinear effects and mediation

These data and specifications reproduce the corresponding example in
Mize, Doan, and Long (2019) and on the [`mecompare` documentation
page](https://www.trentonmize.com/software/mecompare). The code chunks
are not evaluated while building the package website because they
download the public replication data. Copy and run the code
interactively to reproduce the results.

## Load and prepare the data

``` r

library(haven)
library(marginaleffects)
library(suest)

dat <- read_dta(
  "https://tdmize.github.io/data/data/ah4_cme.dta"
)
dat <- zap_labels(dat)

vars <- c(
  "depsympB", "income", "age", "woman", "race", "jobsat"
)
dat <- dat[complete.cases(dat[vars]), ]
dat[c("woman", "race", "jobsat")] <-
  lapply(dat[c("woman", "race", "jobsat")], factor)
```

## Fit the base and mediator models

``` r

base <- lm(
  depsympB ~ income + I(income^2) + age + woman + race,
  data = dat
)

mediator <- lm(
  depsympB ~ income + I(income^2) + age + woman + race + jobsat,
  data = dat
)

combined <- suest(
  base,
  mediator,
  model_names = c("Base", "Mediator")
)
```

## Compare a one-standard-deviation increase in income

``` r

income_sd <- sd(dat$income)

forward_change <- function(amount) {
  function(x) data.frame(lo = x, hi = x + amount)
}

effects <- avg_comparisons(
  combined,
  variables = list(income = forward_change(income_sd)),
  newdata = dat
)

effects
hypotheses(effects, hypothesis = difference ~ revpairwise)
```

The expected effects are approximately -0.816 and -0.648
depressive-symptom units. Their difference is approximately -0.168 with
a robust standard error of 0.022.
