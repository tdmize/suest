# Example 6.6: Different samples

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

gss <- zap_labels(read_dta(
  "https://tdmize.github.io/data/data/gss_cme.dta"
))

vars <- c(
  "helpsickB", "conserv", "faminc", "employed", "woman",
  "age", "college", "married", "parent", "race", "year"
)
dat <- gss[complete.cases(gss[vars]), ]
dat[c("conserv", "employed", "woman", "college",
      "married", "parent", "race")] <-
  lapply(dat[c("conserv", "employed", "woman", "college",
               "married", "parent", "race")], factor)
```

## Fit the same model in two years

``` r

model_1986 <- glm(
  helpsickB ~ conserv + faminc + employed + woman + age +
    college + married + parent + race,
  family = binomial("logit"),
  data = dat,
  subset = year == 1986
)

model_2016 <- glm(
  helpsickB ~ conserv + faminc + employed + woman + age +
    college + married + parent + race,
  family = binomial("logit"),
  data = dat,
  subset = year == 2016
)

combined <- suest(
  model_1986,
  model_2016,
  model_names = c("1986", "2016")
)

combined
```

The samples are disjoint, so the printed object reports zero overlapping
observations and the cross-model covariance block is exactly zero.

## Average within each model’s estimation sample

``` r

nd <- suest_newdata(combined)

effects <- avg_comparisons(
  combined,
  variables = "conserv",
  newdata = nd
)

effects
hypotheses(effects, hypothesis = difference ~ revpairwise)
```

The expected average effects are approximately -0.092 in 1986 and -0.258
in 2016. Their difference is approximately 0.166 with a robust standard
error of 0.039.
