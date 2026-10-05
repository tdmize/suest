# Different model types

Example 6.5 of Mize, Doan, and Long (2019).

This example compares effects across two model types, an ordinal model
and a nominal model, with the same outcome and predictors. It also shows
a custom effect for a continuous predictor: a change in age from 20 to
30, with the other covariates held at their means.

``` r

library(suest)
library(marginaleffects)
library(haven)

gss <- zap_labels(read_dta("https://tdmize.github.io/data/data/gss_cme.dta"))

factorize <- function(data, variables) {
  data[variables] <- lapply(data[variables], factor)
  data
}
```

``` r

vars65 <- c("partyid5", "woman", "edyrs", "age", "parent", "married",
            "faminc", "employed", "region4", "year", "race")
dat65 <- subset(gss, year >= 2010)
dat65 <- dat65[complete.cases(dat65[vars65]), ]
nrow(dat65)
#> [1] 8179
dat65 <- factorize(dat65, c("woman", "parent", "married", "race",
                            "employed", "region4", "year"))

party_levels <- sort(unique(dat65$partyid5))
party_labels <- c("Strong Democrat", "Democrat", "Independent",
                  "Republican", "Strong Republican")
dat65$party_ord <- ordered(dat65$partyid5, levels = party_levels, labels = party_labels)
dat65$party_nom <- factor(dat65$partyid5, levels = party_levels, labels = party_labels)

ordered65 <- MASS::polr(party_ord ~ age + I(age^2) + woman + edyrs + parent +
                          married + race + faminc + employed + region4 + year,
                        data = dat65, method = "logistic", Hess = TRUE)
nominal65 <- nnet::multinom(party_nom ~ age + I(age^2) + woman + edyrs + parent +
                             married + race + faminc + employed + region4 + year,
                           data = dat65, Hess = TRUE, trace = FALSE)
fit65 <- suest(ordered65, nominal65, model_names = c("Ordered", "Multinomial"))
```

The effect of age compares age 20 with age 30, with every other column
of the model matrix held at its mean (including each category of the
factor predictors). This is not the same as averaging a ten-year
increase over everyone’s observed age and covariates, and the difference
matters a lot here because age is included as a quadratic.
`newdata = "mean"` in `marginaleffects` is different too: it uses means
for numeric variables but the most common category for factors. A helper
included with `suest` does the calculation. [View the helper
code](https://github.com/tdmize/suest/blob/main/inst/example-code/example-6-5-atmeans.R).

``` r

source(system.file("example-code", "example-6-5-atmeans.R", package = "suest"))
results65 <- example65_atmeans(fit65, age_lo = 20, age_hi = 30)
results65$effects
#>             category       model      estimate   std.error      p.value
#> 1    Strong Democrat     Ordered  0.0200043395 0.003873903 2.418938e-07
#> 2           Democrat     Ordered  0.0233309158 0.005309554 1.112128e-05
#> 3        Independent     Ordered -0.0030269099 0.000485863 4.665355e-10
#> 4         Republican     Ordered -0.0244327347 0.005230285 2.991588e-06
#> 5  Strong Republican     Ordered -0.0158756107 0.003786173 2.752300e-05
#> 6    Strong Democrat Multinomial  0.0324035507 0.003318210 1.585306e-22
#> 7           Democrat Multinomial -0.0102369435 0.011006328 3.523213e-01
#> 8        Independent Multinomial -0.0004415478 0.010034596 9.649024e-01
#> 9         Republican Multinomial -0.0314658066 0.010558017 2.879885e-03
#> 10 Strong Republican Multinomial  0.0097407472 0.003265056 2.851280e-03
results65$differences
#>            category     estimate   std.error      p.value
#> 1   Strong Democrat -0.012399211 0.003072737 5.454810e-05
#> 2          Democrat  0.033567859 0.009357678 3.342486e-04
#> 3       Independent -0.002585362 0.010024139 7.964736e-01
#> 4        Republican  0.007033072 0.009022795 4.356981e-01
#> 5 Strong Republican -0.025616358 0.003720332 5.758342e-12
```

**Example interpretation:** For someone who is 20 years old, the effect
of a ten-year increase in age differs significantly across the ordinal
and nominal models for three of the five outcome categories. For
example, aging increases the probability of identifying as a strong
Democrat in both models, but the increase is significantly larger in the
nominal model. More strikingly, the effects of age on identifying as a
strong Republican run in opposite directions across the two models.
