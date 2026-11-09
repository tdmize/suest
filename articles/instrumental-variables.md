# Instrumental-variable models

Two-stage least squares models fit with
[`fixest::feols()`](https://lrberge.github.io/fixest/reference/feols.html)
can be combined like any other model. Here median rent in the 50 states
(1980 Census) is regressed on median housing value, treated as
endogenous, and the percent of the population that is urban. The first
model instruments housing value with family income; the second adds the
census regions as instruments.

``` r

library(suest)
library(marginaleffects)
library(haven)
library(fixest)

hsng <- zap_labels(read_dta("https://tdmize.github.io/data/data/hsng2.dta"))

iv1 <- feols(rent ~ pcturban | hsngval ~ faminc, data = hsng)
iv2 <- feols(rent ~ pcturban | hsngval ~ faminc + factor(region), data = hsng)
fit_iv <- suest(iv1, iv2, model_names = c("Income", "Income and region"))
fit_iv
#> Seemingly Unrelated Estimation
#> Models: Income + Income and region 
#> Model types: Income=ivreg, Income and region=ivreg 
#> Model engines: Income=fixest::feols, Income and region=fixest::feols 
#> Comparison scale: model-specific response values 
#> Observations: Income=50, Income and region=50 
#> Overlapping observations: 50 
#> Union observations: 50 
#> Parameters: 6
```

The effect of percent urban in each model, and the difference:

``` r

urban <- avg_slopes(fit_iv, variables = "pcturban", newdata = hsng)
urban
#> 
#>              Group Estimate Std. Error      z Pr(>|z|)   S 2.5 % 97.5 %
#>  Income             -0.5064      0.543 -0.933    0.351 1.5 -1.57  0.557
#>  Income and region   0.0815      0.445  0.183    0.855 0.2 -0.79  0.953
#> 
#> Term: pcturban
#> Type: response
#> Comparison: dY/dX
hypotheses(urban, hypothesis = difference ~ revpairwise)
#> 
#>                      Hypothesis Estimate Std. Error    z Pr(>|z|)   S 2.5 %  97.5 %
#>  (Income) - (Income and region)   -0.588      0.255 -2.3   0.0214 5.5 -1.09 -0.0872
```

**Example interpretation:** Neither model shows a significant effect of
the percent urban on rent. But the effects differ significantly across
the two sets of instruments: adding the regions as instruments raises
the effect by 0.59 ($`p < .05`$).
