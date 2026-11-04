# Group comparisons

Section 4.2 of Mize and Han (forthcoming).

Separate models for each group let *every* effect differ across the
groups. These examples fit a model to each group and compare the effects
across the models.

``` r

library(suest)
library(marginaleffects)
library(haven)

gss <- zap_labels(read_dta("https://tdmize.github.io/data/data/gss_cme.dta"))
wvs <- zap_labels(read_dta("https://tdmize.github.io/data/data/lvm_wvs6.dta"))

factorize <- function(data, variables) {
  data[variables] <- lapply(data[variables], factor)
  data
}
```

## Separate models for men and women (Table 6)

The probability of being very happy among employed GSS respondents from
2000 on, modeled separately for men and for women.

``` r

vars6 <- c("vhappy", "college", "wages", "occprest", "age", "married",
           "parent", "woman", "conserv", "race")
dat6 <- subset(gss, year >= 2000 & employed == 1)
dat6 <- dat6[complete.cases(dat6[vars6]), ]
dat6 <- factorize(dat6, c("college", "married", "parent", "conserv", "race"))

men6 <- glm(vhappy ~ college + married + parent + conserv + race + age + I(age^2),
            family = binomial("logit"), data = dat6, subset = woman == 0)
women6 <- glm(vhappy ~ college + married + parent + conserv + race + age + I(age^2),
              family = binomial("logit"), data = dat6, subset = woman == 1)
fit6 <- suest(men6, women6, model_names = c("Men", "Women"))
```

The effect of every predictor, each averaged over its own group’s
sample, with all of the race contrasts and a one-standard-deviation
increase in age:

``` r

effects6 <- avg_comparisons(fit6,
  variables = list(college = "reference", married = "reference",
                   parent = "reference", conserv = "reference",
                   race = "pairwise", age = sd(dat6$age)),
  newdata = suest_newdata(fit6))
effects6
#> 
#>     Term Group          Contrast Estimate Std. Error      z Pr(>|z|)     S    2.5 %   97.5 %
#>  age     Men   +13.0646159309931  0.00671    0.00806  0.832  0.40533   1.3 -0.00909  0.02251
#>  age     Women +13.0646159309931  0.00905    0.00844  1.072  0.28381   1.8 -0.00750  0.02560
#>  college Men   1 - 0              0.04088    0.01406  2.907  0.00365   8.1  0.01332  0.06845
#>  college Women 1 - 0              0.06428    0.01437  4.474  < 0.001  17.0  0.03612  0.09243
#>  conserv Men   1 - 0              0.04570    0.01361  3.358  < 0.001  10.3  0.01903  0.07237
#>  conserv Women 1 - 0              0.04599    0.01431  3.214  0.00131   9.6  0.01795  0.07403
#>  married Men   1 - 0              0.19668    0.01443 13.631  < 0.001 138.1  0.16840  0.22496
#>  married Women 1 - 0              0.20439    0.01420 14.393  < 0.001 153.6  0.17656  0.23222
#>  parent  Men   1 - 0              0.02299    0.01608  1.429  0.15296   2.7 -0.00854  0.05451
#>  parent  Women 1 - 0             -0.04546    0.01663 -2.733  0.00627   7.3 -0.07806 -0.01286
#>  race    Men   2 - 1             -0.00238    0.02130 -0.112  0.91114   0.1 -0.04413  0.03938
#>  race    Men   3 - 1             -0.03481    0.02094 -1.662  0.09644   3.4 -0.07586  0.00623
#>  race    Men   3 - 2             -0.03243    0.02783 -1.165  0.24389   2.0 -0.08699  0.02212
#>  race    Women 2 - 1             -0.01969    0.01863 -1.057  0.29053   1.8 -0.05619  0.01682
#>  race    Women 3 - 1             -0.01473    0.02342 -0.629  0.52946   0.9 -0.06064  0.03118
#>  race    Women 3 - 2              0.00496    0.02774  0.179  0.85819   0.2 -0.04942  0.05934
#> 
#> Type: response
```

The difference between men and women in each effect:

``` r

hypotheses(effects6, hypothesis = difference ~ revpairwise | term + contrast)
#> 
#>     Term          Contrast      Hypothesis  Estimate Std. Error       z Pr(>|z|)   S   2.5 % 97.5 %
#>  age     +13.0646159309931 (Men) - (Women) -0.002341     0.0117 -0.2005  0.84107 0.2 -0.0252 0.0205
#>  college 1 - 0             (Men) - (Women) -0.023394     0.0201 -1.1636  0.24457 2.0 -0.0628 0.0160
#>  conserv 1 - 0             (Men) - (Women) -0.000289     0.0197 -0.0146  0.98832 0.0 -0.0390 0.0384
#>  married 1 - 0             (Men) - (Women) -0.007709     0.0202 -0.3808  0.70337 0.5 -0.0474 0.0320
#>  parent  1 - 0             (Men) - (Women)  0.068448     0.0231  2.9584  0.00309 8.3  0.0231 0.1138
#>  race    2 - 1             (Men) - (Women)  0.017310     0.0283  0.6117  0.54075 0.9 -0.0382 0.0728
#>  race    3 - 1             (Men) - (Women) -0.020082     0.0314 -0.6391  0.52273 0.9 -0.0817 0.0415
#>  race    3 - 2             (Men) - (Women) -0.037392     0.0393 -0.9515  0.34137 1.6 -0.1144 0.0396
```

**Example interpretation:** Married men and married women both report
being happier than their unmarried counterparts, and college degrees and
conservatism have similar positive effects for both. The one significant
gender difference is the effect of being a parent: it does not affect
the probability of being very happy for men, but women who are parents
are less likely to be very happy than women who are not.

## Separate models for Germany and Nigeria

Acceptance of same-sex relationships in the World Values Survey,
regressed on age (with a squared term) separately in Germany and
Nigeria. Each country has its own coefficients and its own observations,
so the cross-model difference reflects both.

``` r

wvs <- subset(wvs, age >= 18 & country %in% c(276, 566))

germany <- lm(Jgaylesbian ~ age + I(age^2), data = wvs, subset = country == 276)
nigeria <- lm(Jgaylesbian ~ age + I(age^2), data = wvs, subset = country == 566)
fit_countries <- suest(germany, nigeria, model_names = c("Germany", "Nigeria"))

age_effects <- avg_comparisons(fit_countries, variables = list(age = 10),
                               newdata = suest_newdata(fit_countries))
age_effects
#> 
#>    Group Estimate Std. Error      z Pr(>|z|)    S 2.5 %  97.5 %
#>  Germany  -0.4282     0.0417 -10.27   <0.001 79.8 -0.51 -0.3465
#>  Nigeria  -0.0902     0.0358  -2.52   0.0117  6.4 -0.16 -0.0201
#> 
#> Term: age
#> Type: response
#> Comparison: +10
hypotheses(age_effects, hypothesis = difference ~ revpairwise)
#> 
#>             Hypothesis Estimate Std. Error     z Pr(>|z|)    S  2.5 % 97.5 %
#>  (Germany) - (Nigeria)   -0.338     0.0549 -6.15   <0.001 30.3 -0.446  -0.23
```

**Example interpretation:** A ten-year increase in age reduces
acceptance of same-sex relationships in both countries, but much more in
Germany than in Nigeria.

## Reference

Mize, Trenton D. and Bing Han. Forthcoming. “Marginal effects: flexible
methods for interpretation across linear and nonlinear models.” In
*Understanding Data Modeling and Data Analysis*, edited by David
Weakliem. Edward Elgar Publishing.
