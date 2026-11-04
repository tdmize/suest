# Nested logit models and mediation/attenuation

Example 6.2 of Mize, Doan, and Long (2019).

This example uses binary logit models to see how the effect of college
changes as more variables are added: the bivariate model, then a full
set of controls, then wages, then occupational prestige.

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

vars62 <- c("vhappy", "college", "wages", "occprest", "age", "married",
            "parent", "woman", "conserv", "reltrad", "year", "employed")
dat62 <- subset(gss, year >= 2000 & employed == 1)
dat62 <- dat62[complete.cases(dat62[vars62]), ]
nrow(dat62)
#> [1] 9216
dat62 <- factorize(dat62, c("college", "married", "parent", "woman", "conserv",
                            "reltrad", "year"))

model62a <- glm(vhappy ~ college, family = binomial("logit"), data = dat62)
model62b <- glm(vhappy ~ college + married + parent + woman + conserv + reltrad +
                  year + age + I(age^2), family = binomial("logit"), data = dat62)
model62c <- glm(vhappy ~ college + wages + married + parent + woman + conserv +
                  reltrad + year + age + I(age^2), family = binomial("logit"),
                data = dat62)
model62d <- glm(vhappy ~ college + wages + occprest + married + parent + woman +
                  conserv + reltrad + year + age + I(age^2),
                family = binomial("logit"), data = dat62)
fit62 <- suest(model62a, model62b, model62c, model62d,
               model_names = c("Model 1", "Model 2", "Model 3", "Model 4"))

effects62 <- avg_comparisons(fit62, variables = "college", newdata = dat62)
effects62
#> 
#>    Group Estimate Std. Error    z Pr(>|z|)    S   2.5 % 97.5 %
#>  Model 1   0.0718     0.0103 6.95   <0.001 38.0  0.0516 0.0921
#>  Model 2   0.0599     0.0105 5.69   <0.001 26.2  0.0393 0.0806
#>  Model 3   0.0360     0.0111 3.24   0.0012  9.7  0.0142 0.0577
#>  Model 4   0.0192     0.0118 1.63   0.1034  3.3 -0.0039 0.0423
#> 
#> Term: college
#> Type: response
#> Comparison: 1 - 0
```

Whether the effect of college diminishes at each step:

``` r

steps62 <- hypotheses(effects62,
  hypothesis = c("b1 - b2 = 0", "b2 - b3 = 0", "b3 - b4 = 0"))
steps62
#> 
#>  Hypothesis Estimate Std. Error    z Pr(>|z|)    S   2.5 % 97.5 %
#>     b1-b2=0   0.0119    0.00402 2.95  0.00313  8.3 0.00400 0.0198
#>     b2-b3=0   0.0240    0.00403 5.94  < 0.001 28.4 0.01606 0.0319
#>     b3-b4=0   0.0168    0.00445 3.77  < 0.001 12.6 0.00805 0.0255
```

And whether the effect is the same in all four models:

``` r

hypotheses(steps62, joint = TRUE, joint_test = "chisq")
#> 
#> 
#> Joint hypothesis test:
#>   = 0
#>  
#>  ChiSq Pr(>|ChiSq|) Df
#>   60.1       <0.001  3
```

**Example interpretation for the effect of college:** On average, the
probability of being very happy is 0.072 higher for people with a
college degree than for those without one ($`p < .001`$). Each set of
added variables reduces the effect significantly: by 0.012 after adding
the demographic controls ($`p < .01`$), by a further 0.024 after adding
wages, and by a further 0.017 after adding occupational prestige (both
$`p < .001`$). With all of them in the model, the effect of college is
0.019 and no longer statistically significant.
