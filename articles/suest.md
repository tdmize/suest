# Getting started with suest

[`suest()`](https://tdmize.github.io/suest/reference/suest.md) combines
two separately fitted models and supplies the joint robust covariance
matrix needed to test whether predictions, marginal effects, or other
quantities differ across models. The package leaves calculation and
presentation of those quantities to
[`marginaleffects`](https://marginaleffects.com/).

The package implements the framework in Mize, Doan, and Long (2019). The
six replications below correspond to Examples 6.1–6.6 for the
[`mecompare` Stata
command](https://www.trentonmize.com/software/mecompare).

For numeric predictors, a single numeric value in `variables` requests a
forward change of that amount. Thus, `list(income = sd(dat61$income))`
calculates the change from each observed income value to that value plus
one standard deviation. This differs from `"sd"`, which requests a
centered one-standard-deviation contrast around the regressor mean.

## Functions at a glance

- [`suest()`](https://tdmize.github.io/suest/reference/suest.md)
  combines two fitted models and calculates their joint covariance
  matrix.
- [`suest_newdata()`](https://tdmize.github.io/suest/reference/suest_newdata.md)
  stacks the two estimation samples so each model’s effects can be
  averaged over its own sample.
- Standard `marginaleffects` functions—including
  [`predictions()`](https://rdrr.io/pkg/marginaleffects/man/predictions.html),
  [`avg_comparisons()`](https://rdrr.io/pkg/marginaleffects/man/comparisons.html),
  [`avg_slopes()`](https://rdrr.io/pkg/marginaleffects/man/slopes.html),
  [`hypotheses()`](https://rdrr.io/pkg/marginaleffects/man/hypotheses.html),
  and the plotting functions—work with the combined object.

## Installation and setup

``` r

remotes::install_github("tdmize/suest")
```

``` r

library(suest)
library(marginaleffects)
```

## Basic workflow

``` r

dat <- mtcars
dat$am <- factor(dat$am)

base <- glm(am ~ wt, family = binomial("logit"), data = dat)
adjusted <- glm(am ~ wt + hp, family = binomial("logit"), data = dat)
combined <- suest(base, adjusted, model_names = c("Base", "Adjusted"))
combined
#> Seemingly Unrelated Estimation
#> Models: Base + Adjusted 
#> Model types: Base=logit, Adjusted=logit 
#> Model engines: Base=stats::glm, Adjusted=stats::glm 
#> Comparison scale: predicted probabilities 
#> Observations: Base=32, Adjusted=32 
#> Overlapping observations: 32 
#> Union observations: 32 
#> Parameters: 5

effects <- avg_comparisons(combined, variables = "wt", newdata = dat)
effects
#> 
#>     Group Estimate Std. Error     z Pr(>|z|)     S  2.5 % 97.5 %
#>  Adjusted   -0.355     0.0184 -19.3   <0.001 273.3 -0.391 -0.319
#>  Base       -0.277     0.0232 -11.9   <0.001 106.5 -0.322 -0.231
#> 
#> Term: wt
#> Type: response
#> Comparison: +1
hypotheses(effects, hypothesis = difference ~ revpairwise)
#> 
#>           Hypothesis Estimate Std. Error     z Pr(>|z|)    S  2.5 %  97.5 %
#>  (Adjusted) - (Base)  -0.0779     0.0141 -5.54   <0.001 25.0 -0.105 -0.0503
```

The final command tests the difference between the two effects directly.
It does not compare whether one effect is statistically significant
while the other is not.

## Supported models

| Model | Supported functions |
|----|----|
| Linear regression | [`stats::lm()`](https://rdrr.io/r/stats/lm.html) |
| Binary logit and probit | [`stats::glm()`](https://rdrr.io/r/stats/glm.html), [`glm2::glm2()`](https://rdrr.io/pkg/glm2/man/glm2.html) |
| Poisson regression | [`stats::glm()`](https://rdrr.io/r/stats/glm.html), [`glm2::glm2()`](https://rdrr.io/pkg/glm2/man/glm2.html) |
| Negative binomial regression | [`MASS::glm.nb()`](https://rdrr.io/pkg/MASS/man/glm.nb.html) |
| Ordered logit and probit | [`MASS::polr()`](https://rdrr.io/pkg/MASS/man/polr.html), restricted [`ordinal::clm()`](https://rdrr.io/pkg/ordinal/man/clm.html) |
| Multinomial logit | [`nnet::multinom()`](https://rdrr.io/pkg/nnet/man/multinom.html) |

Models may use identical, partially overlapping, or disjoint samples.
Supported cross-family comparisons are logit–probit, logit–linear,
probit–linear, Poisson–negative binomial, and ordered logit–multinomial
logit.

Nonunit weights and offsets are not yet supported. Bias-reduced,
adjusted-score, Firth, and penalized GLM fits are rejected because their
estimating equations differ from ordinary maximum-likelihood GLM scores.

## Replicating Mize, Doan, and Long (2019)

The examples use the public replication data for the [`mecompare` Stata
command](https://www.trentonmize.com/software/mecompare). The data are
downloaded once and reused below.

``` r

library(haven)

ah <- zap_labels(read_dta("https://tdmize.github.io/data/data/ah4_cme.dta"))
gss <- zap_labels(read_dta("https://tdmize.github.io/data/data/gss_cme.dta"))

factorize <- function(data, variables) {
  data[variables] <- lapply(data[variables], factor)
  data
}
```

### Ex 6.1 - Marginal effects to summarize curvilinear relationships and test mediation

This example uses a linear regression model with a nonlinear effect
represented by a polynomial: both income and income squared are included
in the model. We test whether adding the mediator, job satisfaction,
reduces the average effect of income, which would suggest mediation. The
example also shows how to request a custom amount of change for a
continuous predictor—in this case, a one-standard-deviation increase in
income.

``` r

vars61 <- c("depsympB", "income", "inc10", "age", "woman", "race", "college", "jobsat")
dat61 <- ah[complete.cases(ah[vars61]), ]
stopifnot(nrow(dat61) == 4307)
nrow(dat61)
#> [1] 4307
dat61 <- factorize(dat61, c("woman", "race", "jobsat"))

base61 <- lm(depsympB ~ income + I(income^2) + age + woman + race, data = dat61)
mediator61 <- lm(depsympB ~ income + I(income^2) + age + woman + race + jobsat,
                 data = dat61)
fit61 <- suest(base61, mediator61, model_names = c("Base", "Mediator"))

effects61 <- avg_comparisons(fit61,
  variables = list(income = sd(dat61$income)), newdata = dat61)
effects61
#> 
#>     Group Estimate Std. Error      z Pr(>|z|)    S  2.5 % 97.5 %
#>  Base       -0.816     0.0756 -10.80   <0.001 87.8 -0.964 -0.668
#>  Mediator   -0.648     0.0739  -8.77   <0.001 59.0 -0.793 -0.503
#> 
#> Term: income
#> Type: response
#> Comparison: +27.2108562966562
hypotheses(effects61, hypothesis = difference ~ revpairwise)
#> 
#>           Hypothesis Estimate Std. Error     z Pr(>|z|)    S  2.5 % 97.5 %
#>  (Base) - (Mediator)   -0.168     0.0222 -7.59   <0.001 44.8 -0.212 -0.125
```

**Example interpretation for the effect of income:** In model 1, the
marginal effect is -0.816, indicating that a standard-deviation increase
in income is associated with about 0.8 fewer depressive symptoms. The
effect is reduced to -0.648 after accounting for job satisfaction. The
direct cross-model test shows that accounting for job satisfaction
decreases the average effect of income by 0.168, a statistically
significant reduction ($`p < .001`$) that is consistent with mediation.

### Ex 6.2 - Comparing marginal effects across nested logit models

This example uses binary logit models and examines how the effect of
college changes as additional variables are added to the model.

``` r

vars62 <- c("vhappy", "college", "wages", "occprest", "age", "married",
            "parent", "woman", "conserv", "reltrad", "year", "employed")
dat62 <- subset(gss, year >= 2000 & employed == 1)
dat62 <- dat62[complete.cases(dat62[vars62]), ]
stopifnot(nrow(dat62) == 9216)
nrow(dat62)
#> [1] 9216
dat62 <- factorize(dat62, c("college", "married", "parent", "woman", "conserv",
                            "reltrad", "year"))

model62a <- glm(vhappy ~ college, family = binomial("logit"), data = dat62)
model62b <- glm(vhappy ~ college + married + parent + woman + conserv + reltrad +
                  year + age + I(age^2), family = binomial("logit"), data = dat62)
fit62 <- suest(model62a, model62b, model_names = c("Model 1", "Model 2"))

effects62 <- avg_comparisons(fit62, variables = "college", newdata = dat62)
effects62
#> 
#>    Group Estimate Std. Error    z Pr(>|z|)    S  2.5 % 97.5 %
#>  Model 1   0.0718     0.0103 6.95   <0.001 38.0 0.0516 0.0921
#>  Model 2   0.0599     0.0105 5.69   <0.001 26.2 0.0393 0.0806
#> 
#> Term: college
#> Type: response
#> Comparison: 1 - 0
hypotheses(effects62, hypothesis = difference ~ revpairwise)
#> 
#>             Hypothesis Estimate Std. Error    z Pr(>|z|)   S 2.5 % 97.5 %
#>  (Model 1) - (Model 2)   0.0119    0.00402 2.95  0.00313 8.3 0.004 0.0198
```

**Example interpretation for the effect of college:** On average, the
probability of being very happy is 0.072 higher for people with a
college degree than for those without one ($`p < .001`$). Model 2 adds
demographic controls, reducing the average marginal effect of college to
0.060, which remains statistically significant. A direct test shows that
adding the controls decreases the effect of college by 0.012
($`p < .01`$).

### Ex 6.3 - Comparing marginal effects using alternative predictors

This example uses two binary logit models with alternative predictors.
The models use two different ways to measure sexuality: sexual behavior
(`sexbehav`) in one model and sexual identity (`sexident`) in the other.

``` r

vars63 <- c("samesexB", "sexident", "sexbehav", "college",
            "woman", "race", "age", "year")
dat63 <- gss[complete.cases(gss[vars63]), ]
stopifnot(nrow(dat63) == 4921)
nrow(dat63)
#> [1] 4921
dat63 <- factorize(dat63, c("sexident", "sexbehav", "college", "woman", "race", "year"))

behavior63 <- glm(samesexB ~ sexbehav + woman + college + age + race + year,
                  family = binomial("logit"), data = dat63)
identity63 <- glm(samesexB ~ sexident + woman + college + age + race + year,
                  family = binomial("logit"), data = dat63)
fit63 <- suest(behavior63, identity63, model_names = c("Behavior", "Identity"))

effects63 <- avg_comparisons(fit63,
  variables = list(sexbehav = "reference", sexident = "reference"), newdata = dat63)
effects63
#> 
#>      Term    Group Contrast Estimate Std. Error      z Pr(>|z|)     S  2.5 %  97.5 %
#>  sexbehav Behavior    2 - 1  -0.0972     0.0235  -4.14   <0.001  14.8 -0.143 -0.0511
#>  sexbehav Behavior    3 - 1  -0.3620     0.0384  -9.43   <0.001  67.7 -0.437 -0.2867
#>  sexbehav Identity    2 - 1   0.0000         NA     NA       NA    NA     NA      NA
#>  sexbehav Identity    3 - 1   0.0000         NA     NA       NA    NA     NA      NA
#>  sexident Behavior    2 - 1   0.0000         NA     NA       NA    NA     NA      NA
#>  sexident Behavior    3 - 1   0.0000         NA     NA       NA    NA     NA      NA
#>  sexident Identity    2 - 1  -0.2741     0.0415  -6.61   <0.001  34.6 -0.355 -0.1928
#>  sexident Identity    3 - 1  -0.4275     0.0337 -12.67   <0.001 119.8 -0.494 -0.3614
#> 
#> Type: response
```

Because each predictor appears in only one component model,
`marginaleffects` reports zero and `NA` for structurally absent effects.
The nonzero rows give the effects for each operationalization.

``` r

alternative_predictor_hypothesis <- function(x) {
  group <- as.character(x$group)
  term <- as.character(x$term)
  behavior <- x$estimate[term == "sexbehav" & group == "Behavior"]
  identity <- x$estimate[term == "sexident" & group == "Identity"]

  data.frame(
    term = c("Behavior bisexual - Identity bisexual",
             "Behavior gay - Identity gay",
             "(Behavior gay - bisexual) - (Identity gay - bisexual)"),
    estimate = c(behavior[1] - identity[1], behavior[2] - identity[2],
                 (behavior[2] - behavior[1]) - (identity[2] - identity[1])))
}

avg_comparisons(fit63,
  variables = list(sexbehav = "reference", sexident = "reference"),
  newdata = dat63, hypothesis = alternative_predictor_hypothesis)
#> 
#>                                                   Term Estimate Std. Error     z Pr(>|z|)    S
#>  Behavior bisexual - Identity bisexual                   0.1769     0.0425  4.16   <0.001 15.0
#>  Behavior gay - Identity gay                             0.0655     0.0406  1.61   0.1063  3.2
#>  (Behavior gay - bisexual) - (Identity gay - bisexual)  -0.1114     0.0587 -1.90   0.0578  4.1
#>    2.5 %  97.5 %
#>   0.0936 0.26023
#>  -0.0140 0.14507
#>  -0.2264 0.00368
#> 
#> Type: response
```

**Example interpretations:** The contrasts between heterosexual and
gay/lesbian respondents and between bisexual and gay/lesbian respondents
do not differ significantly across the two models. However, the
difference between heterosexual and bisexual respondents is
significantly larger when sexuality is measured using identity than when
it is measured using behavior.

### Ex 6.4 - Comparing marginal effects across different outcomes

This example compares effects across different dependent variables. Two
count outcomes—poor mental-health days and poor physical-health days—are
estimated using negative-binomial models, with effects calculated on the
predicted-rate scale. Effects and cross-model comparisons are calculated
for all predictors. To reproduce the `mecompare` Stata example exactly,
the common sample also excludes observations missing `reltrad`, even
though `reltrad` is not included in either displayed model.

``` r

vars64 <- c("mntlhlth", "physhlth", "woman", "married", "age",
            "faminc", "race", "college", "parent", "reltrad")
dat64 <- gss[complete.cases(gss[vars64]), ]
stopifnot(nrow(dat64) == 5062)
nrow(dat64)
#> [1] 5062

dat64 <- factorize(dat64, c("woman", "married", "parent", "college", "race", "year"))

mental64 <- MASS::glm.nb(mntlhlth ~ woman + married + parent + college + age +
                           faminc + race + year, data = dat64)
physical64 <- MASS::glm.nb(physhlth ~ woman + married + parent + college + age +
                             faminc + race + year, data = dat64)
fit64 <- suest(mental64, physical64, model_names = c("Mental", "Physical"))

variables64 <- list(woman = "reference", married = "reference", parent = "reference",
  college = "reference", age = sd(dat64$age), faminc = sd(dat64$faminc),
  race = "reference", year = "reference")
effects64 <- avg_comparisons(fit64, variables = variables64, newdata = dat64)
effects64
#> 
#>     Term    Group          Contrast Estimate Std. Error      z Pr(>|z|)    S  2.5 %   97.5 %
#>  age     Mental   +13.0755231484462  -0.4623      0.108 -4.294   <0.001 15.8 -0.673 -0.25129
#>  age     Physical +13.0755231484462   0.4912      0.118  4.179   <0.001 15.1  0.261  0.72155
#>  college Mental   1 - 0              -0.8794      0.229 -3.848   <0.001 13.0 -1.327 -0.43145
#>  college Physical 1 - 0              -0.5422      0.189 -2.870   0.0041  7.9 -0.912 -0.17198
#>  faminc  Mental   +35.054788229243   -0.4445      0.118 -3.769   <0.001 12.6 -0.676 -0.21334
#>  faminc  Physical +35.054788229243   -0.3808      0.102 -3.741   <0.001 12.4 -0.580 -0.18125
#>  married Mental   1 - 0              -1.0103      0.230 -4.397   <0.001 16.5 -1.461 -0.55992
#>  married Physical 1 - 0              -0.1591      0.192 -0.827   0.4082  1.3 -0.536  0.21795
#>  parent  Mental   1 - 0               0.2741      0.248  1.107   0.2682  1.9 -0.211  0.75944
#>  parent  Physical 1 - 0              -0.2692      0.219 -1.231   0.2182  2.2 -0.698  0.15931
#>  race    Mental   2 - 1              -1.0159      0.258 -3.941   <0.001 13.6 -1.521 -0.51070
#>  race    Mental   3 - 1              -0.4381      0.356 -1.230   0.2188  2.2 -1.136  0.26022
#>  race    Physical 2 - 1              -0.5308      0.217 -2.443   0.0146  6.1 -0.957 -0.10487
#>  race    Physical 3 - 1               0.1451      0.343  0.423   0.6724  0.6 -0.527  0.81749
#>  woman   Mental   1 - 0               0.9929      0.208  4.774   <0.001 19.1  0.585  1.40062
#>  woman   Physical 1 - 0               0.7734      0.174  4.446   <0.001 16.8  0.432  1.11430
#>  year    Mental   2006 - 2002        -0.9438      0.270 -3.492   <0.001 11.0 -1.474 -0.41408
#>  year    Mental   2010 - 2002        -0.0675      0.318 -0.213   0.8316  0.3 -0.690  0.55484
#>  year    Mental   2014 - 2002        -0.5819      0.302 -1.929   0.0537  4.2 -1.173  0.00922
#>  year    Physical 2006 - 2002        -0.2917      0.225 -1.299   0.1939  2.4 -0.732  0.14839
#>  year    Physical 2010 - 2002         0.3233      0.261  1.240   0.2148  2.2 -0.187  0.83400
#>  year    Physical 2014 - 2002        -0.2286      0.245 -0.933   0.3507  1.5 -0.709  0.25148
#> 
#> Type: response

age64 <- avg_comparisons(fit64, variables = list(age = sd(dat64$age)),
                         newdata = dat64)
hypotheses(age64, hypothesis = difference ~ revpairwise)
#> 
#>             Hypothesis Estimate Std. Error     z Pr(>|z|)    S 2.5 % 97.5 %
#>  (Mental) - (Physical)   -0.953      0.132 -7.22   <0.001 40.8 -1.21 -0.694
```

**Example interpretations:** Women report about 0.99 more poor
mental-health days and 0.77 more poor physical-health days per month
than men. Although the effect of gender is about 0.22 larger for mental
health, the cross-outcome difference is not statistically significant.

Being married significantly reduces poor mental-health days by about
1.01, whereas its effect on poor physical-health days is not
statistically significant. The direct comparison shows that the effect
of marriage is significantly larger for mental health than for physical
health (cross-model difference = -0.851, $`p < .01`$).

Similarly, the effect of age differs significantly across the outcomes:
aging is associated with fewer poor mental-health days but more poor
physical-health days.

### Ex 6.5 - Comparing marginal effects across different model types (ordinal vs nominal)

This example compares effects across two model types: an ordinal model
and a nominal model. The models are otherwise identical, using the same
outcome and predictors. It also illustrates a custom effect for a
continuous predictor: a change in age from 20 to 30, with the other
covariates held at their means.

``` r

vars65 <- c("partyid5", "woman", "edyrs", "age", "parent", "married",
            "faminc", "employed", "region4", "year", "race")
dat65 <- subset(gss, year >= 2010)
dat65 <- dat65[complete.cases(dat65[vars65]), ]
stopifnot(nrow(dat65) == 8179)
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

The estimand here is important. The `mecompare` Stata command compares
age 20 with age 30 and holds the other **model-matrix columns** at their
means. This is not the same as averaging a ten-year forward change over
observed ages and covariates. The distinction is especially
consequential because age enters quadratically.

It is also not reproduced exactly by `newdata = "mean"` in
`marginaleffects`, which uses means for numeric variables but modes for
categorical variables. The helper supplied with `suest` reproduces
Stata’s model-matrix `atmeans` calculation exactly. [View the helper
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

### Ex 6.6 - Comparing marginal effects across different samples or groups

This example compares effects from separate models estimated on
different samples: one sample from 1986 and another from 2016. The
models are otherwise identical, using binary logit with the same outcome
and predictors. The same approach can be used to compare models for
different groups, such as men and women. In R, the samples are selected
with `subset`, and
[`suest_newdata()`](https://tdmize.github.io/suest/reference/suest_newdata.md)
ensures that each marginal effect is averaged over the appropriate
estimation sample.

``` r

vars66 <- c("helpsickB", "polviews", "conserv", "faminc", "employed",
            "woman", "age", "college", "married", "parent", "race", "year")
dat66 <- gss[complete.cases(gss[vars66]), ]
dat66 <- factorize(dat66, c("conserv", "employed", "woman", "college",
                            "married", "parent", "race"))

model1986 <- glm(helpsickB ~ conserv + faminc + employed + woman + age +
                   college + married + parent + race, family = binomial("logit"),
                 data = dat66, subset = year == 1986)
model2016 <- glm(helpsickB ~ conserv + faminc + employed + woman + age +
                   college + married + parent + race, family = binomial("logit"),
                 data = dat66, subset = year == 2016)
stopifnot(nobs(model1986) == 1254, nobs(model2016) == 1670)
c(`1986` = nobs(model1986), `2016` = nobs(model2016))
#> 1986 2016 
#> 1254 1670
fit66 <- suest(model1986, model2016, model_names = c("1986", "2016"))
fit66
#> Seemingly Unrelated Estimation
#> Models: 1986 + 2016 
#> Model types: 1986=logit, 2016=logit 
#> Model engines: 1986=stats::glm, 2016=stats::glm 
#> Comparison scale: predicted probabilities 
#> Observations: 1986=1,254, 2016=1,670 
#> Overlapping observations: 0 
#> Union observations: 2,924 
#> Parameters: 22

effects66 <- avg_comparisons(fit66, variables = "conserv",
                             newdata = suest_newdata(fit66))
effects66
#> 
#>  Group Estimate Std. Error     z Pr(>|z|)    S  2.5 %  97.5 %
#>   1986  -0.0917     0.0296  -3.1  0.00193  9.0 -0.150 -0.0337
#>   2016  -0.2582     0.0253 -10.2  < 0.001 78.7 -0.308 -0.2085
#> 
#> Term: conserv
#> Type: response
#> Comparison: 1 - 0
hypotheses(effects66, hypothesis = difference ~ revpairwise)
#> 
#>       Hypothesis Estimate Std. Error    z Pr(>|z|)    S  2.5 % 97.5 %
#>  (1986) - (2016)    0.166     0.0389 4.28   <0.001 15.7 0.0902  0.243
```

**Example interpretation:** Conservatives had significantly lower
predicted probabilities of saying that government should be responsible
for providing health care in both 1986 (marginal effect = -0.092) and
2016 (marginal effect = -0.258). The gap between conservatives and
nonconservatives increased over time: the conservative effect is
significantly more negative in 2016 than in 1986 (cross-model difference
= 0.166, $`p < .001`$).

## Alternative fitting engines

The statistical model family and fitting engine are handled separately.
This allows comparisons across packages when likelihood scores,
parameter ordering, observation alignment, and prediction scales are
compatible.

### `glm2::glm2()`

``` r

library(glm2)

standard <- glm(am ~ wt + hp, family = binomial("logit"), data = dat)
alternative <- glm2(am ~ wt + hp, family = binomial("logit"), data = dat)
engine_fit <- suest(standard, alternative, model_names = c("glm", "glm2"))
engine_fit
#> Seemingly Unrelated Estimation
#> Models: glm + glm2 
#> Model types: glm=logit, glm2=logit 
#> Model engines: glm=stats::glm, glm2=glm2::glm2 
#> Comparison scale: predicted probabilities 
#> Observations: glm=32, glm2=32 
#> Overlapping observations: 32 
#> Union observations: 32 
#> Parameters: 6
avg_comparisons(engine_fit, variables = "wt", newdata = dat)
#> 
#>  Group Estimate Std. Error     z Pr(>|z|)     S  2.5 % 97.5 %
#>   glm    -0.355     0.0184 -19.3   <0.001 273.3 -0.391 -0.319
#>   glm2   -0.355     0.0184 -19.3   <0.001 273.3 -0.391 -0.319
#> 
#> Term: wt
#> Type: response
#> Comparison: +1
```

### `ordinal::clm()`

[`clm()`](https://rdrr.io/pkg/ordinal/man/clm.html) is supported for
proportional-effects ordered logit and probit models with flexible
thresholds and no scale or nominal-effects formula.

``` r

library(ordinal)

housing <- MASS::housing[rep(seq_len(nrow(MASS::housing)), MASS::housing$Freq),
                         c("Sat", "Infl", "Type", "Cont")]
rownames(housing) <- NULL
housing$Sat <- ordered(housing$Sat, levels = c("Low", "Medium", "High"))

clm_model <- clm(Sat ~ Cont + Infl + Type, data = housing, link = "logit")
polr_model <- MASS::polr(Sat ~ Cont + Infl + Type, data = housing,
                         method = "logistic", Hess = TRUE)
ordinal_fit <- suest(clm_model, polr_model, model_names = c("clm", "polr"))
avg_comparisons(ordinal_fit, variables = "Cont", newdata = housing)
#> 
#>         Group Estimate Std. Error     z Pr(>|z|)    S    2.5 %    97.5 %
#>  clm::Low     -0.07504    0.02010 -3.73   <0.001 12.4 -0.11443 -0.035653
#>  clm::Medium  -0.00374    0.00176 -2.13    0.033  4.9 -0.00718 -0.000303
#>  clm::High     0.07878    0.02071  3.80   <0.001 12.8  0.03819  0.119372
#>  polr::Low    -0.07504    0.02010 -3.73   <0.001 12.4 -0.11443 -0.035653
#>  polr::Medium -0.00374    0.00176 -2.13    0.033  4.9 -0.00718 -0.000303
#>  polr::High    0.07878    0.02071  3.80   <0.001 12.8  0.03819  0.119372
#> 
#> Term: Cont
#> Type: response
#> Comparison: High - Low
```

## For Stata users

The same framework is available in Stata through the [`mecompare` Stata
command](https://www.trentonmize.com/software/mecompare). The R package
uses `marginaleffects` for predictions and effects, while the
`mecompare` Stata command provides a Stata-native workflow for the same
class of cross-model comparisons.
