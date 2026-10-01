# Notes on specific models

Most supported models need nothing special. The notes below cover the
few that do. Fits that `suest` can’t use, such as penalized fits or
models that did not converge, stop with an error that explains why.

## Linear models

Linear models add an `lnvar` parameter (the log of the error variance)
to the combined coefficients and covariance matrix. It has no effect on
predictions or marginal effects.

## Linear panel models

For [`plm::plm()`](https://rdrr.io/pkg/plm/man/plm.html) models,
comparisons of slopes and changes are supported. Confidence intervals
and tests for average predicted levels at the estimation-sample means
are not: their standard errors can be missing or near zero.

## Mixed models fitted with glmmTMB

- Predictions and marginal effects average over the random intercept,
  and over the random slope when there is one.
- Center continuous predictors and keep them in well-scaled units.
  glmmTMB’s covariance matrix can be inaccurate with extreme predictor
  units even when the model reports convergence, and `suest` uses that
  covariance as supplied.
- For slopes of continuous predictors in logit random-slope models, use
  `numderiv = list("fdcenter", eps = 1e-4)` in
  [`avg_slopes()`](https://rdrr.io/pkg/marginaleffects/man/slopes.html)
  for more stable standard errors, and check that nearby step sizes
  (such as `5e-5` and `2e-4`) give similar answers.

## Models fitted with different packages

Models fitted by different packages can be combined, as long as they are
supported models.

``` r

library(suest)
library(marginaleffects)
```

[`glm()`](https://rdrr.io/r/stats/glm.html) and
[`glm2::glm2()`](https://rdrr.io/pkg/glm2/man/glm2.html):

``` r

library(glm2)

dat <- mtcars
dat$am <- factor(dat$am)

standard <- glm(am ~ wt + hp, family = binomial("logit"), data = dat)
alternative <- glm2(am ~ wt + hp, family = binomial("logit"), data = dat)
engine_fit <- suest(standard, alternative, model_names = c("glm", "glm2"))
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

[`ordinal::clm()`](https://rdrr.io/pkg/ordinal/man/clm.html) and
[`MASS::polr()`](https://rdrr.io/pkg/MASS/man/polr.html).
[`clm()`](https://rdrr.io/pkg/ordinal/man/clm.html) models need flexible
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
