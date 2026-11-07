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

## Several random slopes and three-level models

- `lmer()`, `glmer()`, `glmmTMB()`, and
  [`clmm()`](https://rdrr.io/pkg/ordinal/man/clmm.html) models can have
  several random slopes for the grouping variable,
  `(1 + x1 + x2 | school)`, and one level nested in it,
  `(1 + x1 | school) + (1 | school:class)` or `(1 | school/class)`. In
  Stata, these are
  `|| school: x1 x2, covariance(unstructured) || class:`. Independent
  terms, `(1 + x1 || school)`, and glmmTMB’s
  [`diag()`](https://rdrr.io/r/base/diag.html) are Stata’s
  `covariance(independent)`. Crossed random effects and a third grouping
  level are not supported.
- Slopes must be numeric data columns. Make a numeric 0/1 version of a
  factor, and create transformed variables (such as squares) before
  fitting.
- The variance parameters are each term’s log standard deviations,
  `log_sd_school_x1`, and Fisher z correlations,
  `atanh_rho_school_intercept_x1`, as in Stata’s `lns` and `atr`. A
  nested level is named with a dot: `log_sd_school.class_intercept`.
- Nonlinear models with these structures use the Laplace approximation
  (lme4 allows adaptive quadrature only for a single random intercept).
  suest2 requires adaptive quadrature, so to check results in Stata, fit
  each model with `intmethod(laplace) vce(cluster school)`; its standard
  errors equal that model’s block of the
  [`suest()`](https://tdmize.github.io/suest/reference/suest.md)
  covariance. Linear models agree with `mixed` exactly.
- Fits near a boundary, with a variance near zero or correlations near
  -1 or 1, are refused even when lme4 does not flag them as singular.
  Simplify the random-effects structure.

## Survey models

- Binary probit and cloglog (`svyglm()`), ordered (`svyolr()`),
  multinomial, negative-binomial, and generalized ordered
  (`svy_vglm()`), and survival (`svysurvreg()`) models use exact
  observed information, as Stata does. Their standard errors can differ
  from the fitted model’s own: `svyglm()` and VGAM use expected
  information, `svyolr()` uses a numerical Hessian, and for
  interval-censored outcomes survival’s scale residuals are incorrect.
- These fits stop at fairly loose tolerances. For results that agree
  closely with Stata, tighten them: `control = list(reltol = 1e-12)` in
  `svyolr()` and
  `control = VGAM::vglm.control(criterion = "coefficients", epsilon = 1e-9)`
  in `svy_vglm()`.
- `svy_vglm()` multinomial models use the last category as the reference
  unless `multinomial(refLevel = )` says otherwise.
- [`suest()`](https://tdmize.github.io/suest/reference/suest.md) finds
  an `svyivreg()` model’s design from its call, so keep that design
  object available.

## Stata models without an R version

For these Stata models, `suest` supports a model that answers the same
question. The last column fits the same model in Stata.

| Stata model | Use in R | Same model in Stata |
|----|----|----|
| `xtlogit, fe` | Random-intercept logit with each time-varying predictor’s person mean added (`glmer()`, `glmmTMB()`, or `pglm()`); see [Within-person effects](https://tdmize.github.io/suest/articles/within-person.md) | `xtlogit, re` or `melogit` with the person means |
| `xtmlogit, fe` | [`nnet::multinom()`](https://rdrr.io/pkg/nnet/man/multinom.html) with the person means, clustered on person | `mlogit` with the person means, `vce(cluster id)` |
| `xtmlogit, re` | [`nnet::multinom()`](https://rdrr.io/pkg/nnet/man/multinom.html) clustered on person | `mlogit, vce(cluster id)` |
| `xtnbreg, fe` | Fixed-effects Poisson, [`fixest::fepois()`](https://lrberge.github.io/fixest/reference/feglm.html) | `xtpoisson, fe` |
| `xtnbreg, re` | Random-intercept negative binomial, `glmmTMB()` | `menbreg` |
| `xtnbreg, pa` | Population-averaged Poisson, [`geepack::geeglm()`](https://rdrr.io/pkg/geepack/man/geeglm.html) | `xtpoisson, pa vce(robust)` |
| `mestreg` | `survreg()` clustered on the group | `streg, vce(cluster group)` |
| `biprobit` with different predictors | Two probit models | Two `probit` models |

Fixed-effects Poisson with robust standard errors is valid for
overdispersed counts, and `xtnbreg, fe` is not a true fixed-effects
model (Allison and Waterman 2002). The multinomial and survival
stand-ins have no random effects, so they estimate effects averaged over
people rather than for a given person. Two probit models estimate the
same effects on each outcome as a bivariate probit, without the
correlation between their errors.

Allison, P. D., & Waterman, R. P. (2002). Fixed-effects negative
binomial regression models. *Sociological Methodology*, 32, 247-265.

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
