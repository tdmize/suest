# Multilevel models

``` r

library(suest)
library(marginaleffects)
library(haven)
library(lme4)

ah <- zap_labels(read_dta("https://tdmize.github.io/data/data/ah_si_2026.dta"))
bangladesh <- zap_labels(read_dta("https://tdmize.github.io/data/data/bangladesh.dta"))
```

## Ordinary and multilevel logit

The coefficients of a multilevel logit are on a different scale than
those of an ordinary logit, so they cannot be compared directly; the
marginal effects can. The example uses wave 1 of Add Health: about 5,800
students in 132 schools. The outcome is being a current smoker. The
predictors are sex, age, race, and whether the student’s mother has a
college degree, plus a school-level variable: the share of students in
the school whose mother has a college degree. The same model is fit as
an ordinary logit, which ignores the schools, and as a logit with a
random intercept for the school.

``` r

ah <- subset(ah, wave == 1)
ah$school <- ah$cluster
ah <- ah[complete.cases(ah[c("smoker", "male", "age", "race", "momeduc")]), ]
ah$momcollege <- as.numeric(ah$momeduc >= 4)
ah$school_momcollege <- ave(ah$momcollege, ah$school)
ah[c("male", "race", "momcollege")] <- lapply(ah[c("male", "race", "momcollege")], factor)

ordinary <- glm(smoker ~ male + age + race + momcollege + school_momcollege,
                family = binomial("logit"), data = ah)
multilevel <- glmer(smoker ~ male + age + race + momcollege + school_momcollege +
                      (1 | school), family = binomial("logit"), data = ah)
fit_school <- suest(ordinary, multilevel, model_names = c("Ordinary", "Multilevel"))
```

[`suest()`](https://tdmize.github.io/suest/reference/suest.md) clusters
the standard errors on the school. The marginal effects average over the
random intercept, so they are changes in predicted probabilities, as
after the ordinary logit. Age and the school’s share of college-educated
mothers are increased by one standard deviation:

``` r

effects_school <- avg_comparisons(fit_school,
  variables = list(male = "reference", age = sd(ah$age), race = "reference",
                   momcollege = "reference",
                   school_momcollege = sd(ah$school_momcollege)),
  newdata = ah)
effects_school
#> 
#>               Term      Group           Contrast Estimate Std. Error         z Pr(>|z|)     S
#>  age               Ordinary   +1.7546893968729    0.07252    0.00754   9.62031   <0.001  70.4
#>  age               Multilevel +1.7546893968729    0.06549    0.00802   8.16757   <0.001  51.5
#>  male              Ordinary   1 - 0               0.01021    0.01167   0.87515    0.381   1.4
#>  male              Multilevel 1 - 0               0.01023    0.01141   0.89679    0.370   1.4
#>  momcollege        Ordinary   1 - 0              -0.06076    0.01001  -6.07128   <0.001  29.6
#>  momcollege        Multilevel 1 - 0              -0.06010    0.00995  -6.04172   <0.001  29.3
#>  race              Ordinary   2 - 1              -0.17637    0.01366 -12.91146   <0.001 124.3
#>  race              Ordinary   3 - 1              -0.11186    0.02200  -5.08417   <0.001  21.4
#>  race              Ordinary   4 - 1              -0.13385    0.03421  -3.91259   <0.001  13.4
#>  race              Ordinary   5 - 1              -0.00479    0.04078  -0.11739    0.907   0.1
#>  race              Multilevel 2 - 1              -0.16663    0.01380 -12.07166   <0.001 109.0
#>  race              Multilevel 3 - 1              -0.08536    0.02300  -3.71195   <0.001  12.2
#>  race              Multilevel 4 - 1              -0.12910    0.03356  -3.84698   <0.001  13.0
#>  race              Multilevel 5 - 1               0.00010    0.04468   0.00224    0.998   0.0
#>  school_momcollege Ordinary   +0.179735627486497 -0.01066    0.00936  -1.13910    0.255   2.0
#>  school_momcollege Multilevel +0.179735627486497 -0.00797    0.00947  -0.84165    0.400   1.3
#>    2.5 %   97.5 %
#>   0.0577  0.08729
#>   0.0498  0.08121
#>  -0.0127  0.03308
#>  -0.0121  0.03260
#>  -0.0804 -0.04114
#>  -0.0796 -0.04060
#>  -0.2031 -0.14959
#>  -0.1550 -0.06874
#>  -0.2009 -0.06680
#>  -0.0847  0.07515
#>  -0.1937 -0.13958
#>  -0.1304 -0.04029
#>  -0.1949 -0.06333
#>  -0.0875  0.08768
#>  -0.0290  0.00769
#>  -0.0265  0.01059
#> 
#> Type: response
hypotheses(effects_school, hypothesis = difference ~ revpairwise | term + contrast)
#> 
#>               Term           Contrast                Hypothesis   Estimate Std. Error       z
#>  age               +1.7546893968729   (Ordinary) - (Multilevel)  0.0070233   0.002856  2.4589
#>  male              1 - 0              (Ordinary) - (Multilevel) -0.0000231   0.000970 -0.0238
#>  momcollege        1 - 0              (Ordinary) - (Multilevel) -0.0006572   0.000727 -0.9034
#>  race              2 - 1              (Ordinary) - (Multilevel) -0.0097364   0.004796 -2.0302
#>  race              3 - 1              (Ordinary) - (Multilevel) -0.0265069   0.009107 -2.9107
#>  race              4 - 1              (Ordinary) - (Multilevel) -0.0047501   0.007296 -0.6511
#>  race              5 - 1              (Ordinary) - (Multilevel) -0.0048876   0.008318 -0.5876
#>  school_momcollege +0.179735627486497 (Ordinary) - (Multilevel) -0.0026925   0.001871 -1.4392
#>  Pr(>|z|)   S    2.5 %    97.5 %
#>   0.01394 6.2  0.00143  0.012622
#>   0.98101 0.0 -0.00192  0.001878
#>   0.36632 1.4 -0.00208  0.000769
#>   0.04233 4.6 -0.01914 -0.000337
#>   0.00361 8.1 -0.04436 -0.008658
#>   0.51499 1.0 -0.01905  0.009549
#>   0.55681 0.8 -0.02119  0.011416
#>   0.15010 2.7 -0.00636  0.000974
```

**Example interpretation:** Most effects are nearly the same in the two
models. The ordinary logit gives somewhat larger effects of age and
race: for example, its Latinx-White difference is 11.2 percentage
points, compared with 8.5 in the multilevel logit (the difference of 2.6
points has $`p < .01`$).

## Random-intercept logit and probit

The Bangladesh Fertility Survey: contraceptive use among about 1,900
women in 60 districts, modeled with a random-intercept logit and a
random-intercept probit.

``` r

bangladesh$urban <- factor(bangladesh$urban)

logit_re <- glmer(c_use ~ urban + age + (1 | district),
                  family = binomial("logit"), data = bangladesh, nAGQ = 5)
probit_re <- glmer(c_use ~ urban + age + (1 | district),
                   family = binomial("probit"), data = bangladesh, nAGQ = 5)
fit_district <- suest(logit_re, probit_re, model_names = c("Logit", "Probit"))

effects_district <- avg_comparisons(fit_district,
  variables = list(urban = "reference", age = 1), newdata = bangladesh)
effects_district
#> 
#>   Term  Group Contrast Estimate Std. Error    z Pr(>|z|)    S     2.5 %  97.5 %
#>  age   Logit     +1     0.00200    0.00129 1.56    0.120  3.1 -0.000519 0.00452
#>  age   Probit    +1     0.00201    0.00130 1.54    0.123  3.0 -0.000543 0.00457
#>  urban Logit     1 - 0  0.14988    0.03749 4.00   <0.001 13.9  0.076404 0.22335
#>  urban Probit    1 - 0  0.15051    0.03726 4.04   <0.001 14.2  0.077480 0.22354
#> 
#> Type: response
hypotheses(effects_district, hypothesis = difference ~ revpairwise | term + contrast)
#> 
#>   Term Contrast         Hypothesis   Estimate Std. Error      z Pr(>|z|)   S      2.5 %    97.5 %
#>  age      +1    (Logit) - (Probit) -0.0000113  0.0000331 -0.342    0.732 0.4 -0.0000763 0.0000536
#>  urban    1 - 0 (Logit) - (Probit) -0.0006333  0.0004510 -1.404    0.160 2.6 -0.0015172 0.0002507
```

**Example interpretation:** The logit and probit models give the same
marginal effects. Living in an urban area raises the probability of
using contraception by 15 percentage points in both models, and neither
effect differs across them.
