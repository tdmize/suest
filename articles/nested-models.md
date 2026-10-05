# Nested logit models and mediation/attenuation

Example 6.2 of Mize, Doan, and Long (2019).

This example uses binary logit models to see how the effect of college
changes when more variables are added to the model.

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
