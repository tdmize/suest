# Different outcomes

Example 6.4 of Mize, Doan, and Long (2019).

This example compares effects across different outcomes. Two counts,
poor mental-health days and poor physical-health days, are modeled with
negative binomial regression, and effects are on the predicted-rate
scale. Effects and cross-model comparisons are calculated for all
predictors.

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

vars64 <- c("mntlhlth", "physhlth", "woman", "married", "age",
            "faminc", "race", "college", "parent", "reltrad")
dat64 <- gss[complete.cases(gss[vars64]), ]
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

`age = sd(dat64$age)` asks for a one-standard-deviation increase from
each person’s observed age.

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
