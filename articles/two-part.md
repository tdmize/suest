# Two-part models

A two-part model fits whether the outcome is positive (the first part)
and the outcome among those with a positive value (the second part).
[`suest_twopart()`](https://tdmize.github.io/suest/reference/suest_twopart.md)
combines the two parts, so predictions and marginal effects are on the
overall expected outcome, `Pr(y > 0) * E(y | y > 0)`, with standard
errors from the joint covariance matrix of the two parts.

The outcome is total medical expenditures of adults in the 2004 Medical
Expenditure Panel Survey; 18% had none. The first part is a logit for
having any expenditures; the second part is a Gamma model with a log
link for the amount among those with expenditures. Standard errors are
clustered on the dwelling unit, which family members share.

``` r

library(suest)
library(marginaleffects)
library(haven)

meps <- zap_labels(read_dta("https://tdmize.github.io/data/data/meps.dta"))
meps[c("female", "ins_unins", "anylim")] <-
  lapply(meps[c("female", "ins_unins", "anylim")], factor)

any <- glm(I(exp_tot > 0) ~ age + female + lninc + pcs12 + mcs12 +
             ins_unins + anylim, family = binomial("logit"), data = meps)
amount <- glm(exp_tot ~ age + female + lninc + pcs12 + mcs12 + ins_unins +
                anylim, family = Gamma(link = "log"), data = meps,
              subset = exp_tot > 0)
fit <- suest_twopart(any, amount, model_names = c("Any", "Amount"),
                     cluster = "duid")
fit
#> Two-part model: first part Any, second part Amount
#> Seemingly Unrelated Estimation
#> Models: Any + Amount 
#> Model types: Any=logit, Amount=glm 
#> Model engines: Any=stats::glm, Amount=stats::glm 
#> Comparison scale: model-specific response values 
#> Clusters: 10,940 
#> Observations: Any=19,386, Amount=15,946 
#> Overlapping observations: 15,946 
#> Union observations: 19,386 
#> Parameters: 16
```

For a Poisson pseudo-likelihood second part, use
`family = quasipoisson()` instead.

## Effects on expected expenditures

``` r

avg_slopes(fit, variables = c("age", "female", "ins_unins", "anylim"))
#> 
#>       Term Contrast Estimate Std. Error      z Pr(>|z|)     S   2.5 %  97.5 %
#>  age          dY/dX     62.5       3.92  15.97   <0.001 188.2    54.9    70.2
#>  anylim       1 - 0   1442.2     143.37  10.06   <0.001  76.7  1161.2  1723.1
#>  female       1 - 0    847.7     115.51   7.34   <0.001  42.1   621.3  1074.1
#>  ins_unins    1 - 0  -2613.2     111.69 -23.40   <0.001 399.8 -2832.1 -2394.3
#> 
#> Type: response
```

## Effects in each part

`type = "prob"` gives the effects on the probability of any
expenditures, and `type = "conditional"` the effects on expenditures
among those with any:

``` r

avg_slopes(fit, variables = "ins_unins", type = "prob")
#> 
#>  Estimate Std. Error     z Pr(>|z|)     S  2.5 % 97.5 %
#>    -0.242    0.00855 -28.4   <0.001 585.7 -0.259 -0.226
#> 
#> Term: ins_unins
#> Type: prob
#> Comparison: 1 - 0
avg_slopes(fit, variables = "ins_unins", type = "conditional")
#> 
#>  Estimate Std. Error     z Pr(>|z|)     S 2.5 % 97.5 %
#>     -2484        134 -18.6   <0.001 253.8 -2746  -2222
#> 
#> Term: ins_unins
#> Type: conditional
#> Comparison: 1 - 0
```

## Comparing groups

The effect of being uninsured for men (`female` 0) and women (`female`
1), and the difference between them:

``` r

by_sex <- avg_comparisons(fit, variables = "ins_unins", by = "female")
by_sex
#> 
#>  female Estimate Std. Error     z Pr(>|z|)     S 2.5 % 97.5 %
#>       0    -2170       95.6 -22.7   <0.001 376.3 -2357  -1982
#>       1    -2977      139.3 -21.4   <0.001 334.1 -3250  -2703
#> 
#> Term: ins_unins
#> Type: response
#> Comparison: 1 - 0
hypotheses(by_sex, hypothesis = difference ~ pairwise)
#> 
#>   Hypothesis Estimate Std. Error     z Pr(>|z|)    S 2.5 % 97.5 %
#>  (b2) - (b1)     -807       96.6 -8.35   <0.001 53.7  -996   -617
```
