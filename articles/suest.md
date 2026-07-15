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
forward_change <- function(amount) function(x) data.frame(lo = x, hi = x + amount)
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

vars61 <- c("depsympB", "income", "age", "woman", "race", "jobsat")
dat61 <- ah[complete.cases(ah[vars61]), ]
dat61 <- factorize(dat61, c("woman", "race", "jobsat"))

base61 <- lm(depsympB ~ income + I(income^2) + age + woman + race, data = dat61)
mediator61 <- lm(depsympB ~ income + I(income^2) + age + woman + race + jobsat,
                 data = dat61)
fit61 <- suest(base61, mediator61, model_names = c("Base", "Mediator"))

effects61 <- avg_comparisons(fit61,
  variables = list(income = forward_change(sd(dat61$income))), newdata = dat61)
effects61
#> 
#>     Group Estimate Std. Error      z Pr(>|z|)    S  2.5 % 97.5 %
#>  Base       -0.816     0.0756 -10.80   <0.001 87.8 -0.964 -0.668
#>  Mediator   -0.648     0.0739  -8.77   <0.001 59.0 -0.793 -0.503
#> 
#> Term: income
#> Type: response
#> Comparison: custom
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

vars62 <- c("vhappy", "college", "age", "married", "parent", "woman",
            "conserv", "reltrad", "year", "employed")
dat62 <- subset(gss, year >= 2000 & employed == 1)
dat62 <- dat62[complete.cases(dat62[vars62]), ]
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
#>  Model 1   0.0709    0.00958 7.41   <0.001 42.8 0.0522 0.0897
#>  Model 2   0.0558    0.00974 5.73   <0.001 26.6 0.0367 0.0749
#> 
#> Term: college
#> Type: response
#> Comparison: 1 - 0
hypotheses(effects62, hypothesis = difference ~ revpairwise)
#> 
#>             Hypothesis Estimate Std. Error    z Pr(>|z|)    S   2.5 % 97.5 %
#>  (Model 1) - (Model 2)   0.0151    0.00376 4.02   <0.001 14.1 0.00777 0.0225
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
for all predictors.

``` r

vars64 <- c("mntlhlth", "physhlth", "woman", "married", "age",
            "faminc", "race", "college", "parent", "year")
dat64 <- gss[complete.cases(gss[vars64]), ]
dat64 <- factorize(dat64, c("woman", "married", "parent", "college", "race", "year"))

mental64 <- MASS::glm.nb(mntlhlth ~ woman + married + parent + college + age +
                           faminc + race + year, data = dat64)
physical64 <- MASS::glm.nb(physhlth ~ woman + married + parent + college + age +
                             faminc + race + year, data = dat64)
fit64 <- suest(mental64, physical64, model_names = c("Mental", "Physical"))

variables64 <- list(woman = "reference", married = "reference", parent = "reference",
  college = "reference", age = forward_change(sd(dat64$age)),
  faminc = forward_change(sd(dat64$faminc)), race = "reference", year = "reference")
effects64 <- avg_comparisons(fit64, variables = variables64, newdata = dat64)
effects64
#> 
#>     Term    Group    Contrast Estimate Std. Error      z Pr(>|z|)    S  2.5 % 97.5 %
#>  age     Mental   custom       -0.4946      0.105 -4.689  < 0.001 18.5 -0.701 -0.288
#>  age     Physical custom        0.4393      0.114  3.863  < 0.001 13.1  0.216  0.662
#>  college Mental   1 - 0        -0.9213      0.226 -4.077  < 0.001 14.4 -1.364 -0.478
#>  college Physical 1 - 0        -0.5483      0.188 -2.911  0.00361  8.1 -0.918 -0.179
#>  faminc  Mental   custom       -0.5254      0.115 -4.552  < 0.001 17.5 -0.752 -0.299
#>  faminc  Physical custom       -0.3991      0.101 -3.953  < 0.001 13.7 -0.597 -0.201
#>  married Mental   1 - 0        -1.0044      0.227 -4.422  < 0.001 16.6 -1.450 -0.559
#>  married Physical 1 - 0        -0.1781      0.191 -0.932  0.35113  1.5 -0.552  0.196
#>  parent  Mental   1 - 0         0.3691      0.241  1.532  0.12565  3.0 -0.103  0.841
#>  parent  Physical 1 - 0        -0.2889      0.214 -1.347  0.17783  2.5 -0.709  0.131
#>  race    Mental   2 - 1        -1.1590      0.249 -4.651  < 0.001 18.2 -1.647 -0.671
#>  race    Mental   3 - 1        -0.4042      0.359 -1.125  0.26044  1.9 -1.108  0.300
#>  race    Physical 2 - 1        -0.5833      0.215 -2.707  0.00678  7.2 -1.006 -0.161
#>  race    Physical 3 - 1         0.0776      0.337  0.231  0.81759  0.3 -0.582  0.738
#>  woman   Mental   1 - 0         0.9317      0.205  4.535  < 0.001 17.4  0.529  1.334
#>  woman   Physical 1 - 0         0.6997      0.173  4.054  < 0.001 14.3  0.361  1.038
#>  year    Mental   2006 - 2002  -0.8782      0.267 -3.287  0.00101  9.9 -1.402 -0.355
#>  year    Mental   2010 - 2002  -0.0951      0.309 -0.308  0.75803  0.4 -0.700  0.510
#>  year    Mental   2014 - 2002  -0.4676      0.300 -1.557  0.11948  3.1 -1.056  0.121
#>  year    Physical 2006 - 2002  -0.3091      0.222 -1.391  0.16408  2.6 -0.745  0.126
#>  year    Physical 2010 - 2002   0.3587      0.258  1.393  0.16368  2.6 -0.146  0.864
#>  year    Physical 2014 - 2002  -0.2117      0.243 -0.871  0.38365  1.4 -0.688  0.265
#> 
#> Type: response

age64 <- avg_comparisons(fit64, variables = list(age = forward_change(sd(dat64$age))),
                         newdata = dat64)
hypotheses(age64, hypothesis = difference ~ revpairwise)
#> 
#>             Hypothesis Estimate Std. Error     z Pr(>|z|)    S 2.5 % 97.5 %
#>  (Mental) - (Physical)   -0.934      0.128 -7.29   <0.001 41.5 -1.19 -0.683
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
dat65 <- factorize(dat65, c("woman", "parent", "married", "race",
                            "employed", "region4", "year"))

party_levels <- sort(unique(dat65$partyid5))
dat65$party_ord <- ordered(dat65$partyid5, levels = party_levels)
dat65$party_nom <- factor(dat65$partyid5, levels = party_levels)

ordered65 <- MASS::polr(party_ord ~ age + I(age^2) + woman + edyrs + parent +
                          married + race + faminc + employed + region4 + year,
                        data = dat65, method = "logistic", Hess = TRUE)
nominal65 <- nnet::multinom(party_nom ~ age + I(age^2) + woman + edyrs + parent +
                             married + race + faminc + employed + region4 + year,
                           data = dat65, Hess = TRUE, trace = FALSE)
fit65 <- suest(ordered65, nominal65, model_names = c("Ordered", "Multinomial"))
```

For an empirical-distribution average, `marginaleffects` calculates
category-specific ten-year age changes directly:

``` r

avg_comparisons(fit65, variables = list(age = forward_change(10)), newdata = dat65)
#> 
#>           Group Estimate Std. Error       z Pr(>|z|)     S     2.5 %    97.5 %
#>  Ordered::0      0.00527    0.00190   2.779  0.00546   7.5  0.001553  8.99e-03
#>  Ordered::1      0.00283    0.00163   1.736  0.08262   3.6 -0.000366  6.03e-03
#>  Ordered::2     -0.00100    0.00022  -4.557  < 0.001  17.6 -0.001433 -5.71e-04
#>  Ordered::3     -0.00428    0.00189  -2.270  0.02322   5.4 -0.007977 -5.84e-04
#>  Ordered::4     -0.00282    0.00146  -1.934  0.05313   4.2 -0.005674  3.81e-05
#>  Multinomial::0  0.02945    0.00304   9.695  < 0.001  71.4  0.023495  3.54e-02
#>  Multinomial::1 -0.01750    0.00335  -5.216  < 0.001  22.4 -0.024073 -1.09e-02
#>  Multinomial::2 -0.02667    0.00229 -11.663  < 0.001 102.0 -0.031148 -2.22e-02
#>  Multinomial::3 -0.00158    0.00338  -0.467  0.64076   0.6 -0.008191  5.04e-03
#>  Multinomial::4  0.01629    0.00282   5.781  < 0.001  27.0  0.010768  2.18e-02
#> 
#> Term: age
#> Type: response
#> Comparison: custom
```

The published comparison evaluates a change from age 20 to age 30 with
the other model-matrix columns held at their means. The full acceptance
suite reproduces the following benchmark:

| Category | Ordered effect | Ordered SE | Multinomial effect | Multinomial SE | Difference | Difference SE |
|----|---:|---:|---:|---:|---:|---:|
| Strong Democrat | 0.020 | 0.004 | 0.032 | 0.003 | -0.012 | 0.003 |
| Democrat | 0.023 | 0.005 | -0.010 | 0.011 | 0.034 | 0.009 |
| Independent | -0.003 | 0.000 | -0.000 | 0.010 | -0.003 | 0.010 |
| Republican | -0.024 | 0.005 | -0.031 | 0.011 | 0.007 | 0.009 |
| Strong Republican | -0.016 | 0.004 | 0.010 | 0.003 | -0.026 | 0.004 |

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

vars66 <- c("helpsickB", "conserv", "faminc", "employed", "woman",
            "age", "college", "married", "parent", "race", "year")
dat66 <- gss[complete.cases(gss[vars66]), ]
dat66 <- factorize(dat66, c("conserv", "employed", "woman", "college",
                            "married", "parent", "race"))

model1986 <- glm(helpsickB ~ conserv + faminc + employed + woman + age +
                   college + married + parent + race, family = binomial("logit"),
                 data = dat66, subset = year == 1986)
model2016 <- glm(helpsickB ~ conserv + faminc + employed + woman + age +
                   college + married + parent + race, family = binomial("logit"),
                 data = dat66, subset = year == 2016)
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
