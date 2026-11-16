# Alternative predictors

Example 6.3 of Mize, Doan, and Long (2019).

This example uses two binary logit models with alternative predictors:
two different ways to measure sexuality, sexual behavior (`sexbehav`) in
one model and sexual identity (`sexident`) in the other.

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

vars63 <- c("samesexB", "sexident", "sexbehav", "college",
            "woman", "race", "age", "year")
dat63 <- gss[complete.cases(gss[vars63]), ]
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

Because each predictor is in only one of the models, `marginaleffects`
reports zero and `NA` for the effects that don’t exist in a model. The
nonzero rows are the effects for each measure.

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

The article’s Table 4 also compares the predicted probabilities from the
two models. Setting both measures to the same category gives each
model’s prediction for that category:

``` r

orientation <- levels(dat63$sexident)
pred63 <- avg_predictions(fit63,
  variables = list(sexbehav = orientation, sexident = orientation),
  newdata = dat63)
pred63[pred63$sexbehav == pred63$sexident, ]
#> 
#>     Group sexbehav sexident Estimate Std. Error     z Pr(>|z|)     S  2.5 % 97.5 %
#>  Behavior        1        1    0.533    0.00710 75.07  < 0.001   Inf 0.5191  0.547
#>  Behavior        2        2    0.436    0.02236 19.49  < 0.001 278.6 0.3920  0.480
#>  Behavior        3        3    0.171    0.03772  4.53  < 0.001  17.4 0.0971  0.245
#>  Identity        1        1    0.531    0.00686 77.36  < 0.001   Inf 0.5176  0.545
#>  Identity        2        2    0.257    0.04085  6.29  < 0.001  31.6 0.1769  0.337
#>  Identity        3        3    0.104    0.03302  3.14  0.00171   9.2 0.0388  0.168
#> 
#> Type: response
```

The difference between the two models in each category:

``` r

orientation_differences <- function(x) {
  same <- x$sexbehav == x$sexident
  data.frame(
    term = c("Heterosexual", "Bisexual", "Gay/Lesbian"),
    estimate = x$estimate[same & x$group == "Behavior"] -
      x$estimate[same & x$group == "Identity"])
}

avg_predictions(fit63,
  variables = list(sexbehav = orientation, sexident = orientation),
  newdata = dat63, hypothesis = orientation_differences)
#> 
#>          Term Estimate Std. Error     z Pr(>|z|)    S    2.5 %  97.5 %
#>  Heterosexual  0.00192    0.00193 0.996   0.3194  1.6 -0.00186 0.00569
#>  Bisexual      0.17883    0.04136 4.324   <0.001 16.0  0.09778 0.25989
#>  Gay/Lesbian   0.06745    0.04013 1.681   0.0928  3.4 -0.01121 0.14611
#> 
#> Type: response
```

Only for people classified as bisexual do the two measures give a
significantly different predicted probability.
