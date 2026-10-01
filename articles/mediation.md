# Curvilinear effects and mediation/attenuation

Example 6.1 of Mize, Doan, and Long (2019). The [same example in
Stata](https://www.trentonmize.com/software/mecompare/examples/mediation.html)
uses `mecompare`.

This example uses a linear regression with a nonlinear effect of income,
included as income and income squared. We test whether adding a
mediator, job satisfaction, reduces the average effect of income, which
would suggest mediation. It also shows how to ask for a custom amount of
change in a continuous predictor: here, a one-standard-deviation
increase in income.

``` r

library(suest)
library(marginaleffects)
library(haven)

ah <- zap_labels(read_dta("https://tdmize.github.io/data/data/ah4_cme.dta"))

factorize <- function(data, variables) {
  data[variables] <- lapply(data[variables], factor)
  data
}
```

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

A single number in `variables` asks for a change of that amount from
each observed value. So `list(income = sd(dat61$income))` compares each
person’s observed income with that income plus one standard deviation.
This is different from `"sd"`, which asks for a one-standard-deviation
change centered on the mean.

**Example interpretation for the effect of income:** In model 1, the
marginal effect is -0.816, indicating that a standard-deviation increase
in income is associated with about 0.8 fewer depressive symptoms. The
effect is reduced to -0.648 after accounting for job satisfaction. The
direct cross-model test shows that accounting for job satisfaction
decreases the average effect of income by 0.168, a statistically
significant reduction ($`p < .001`$) that is consistent with mediation.
