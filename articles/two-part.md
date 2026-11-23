# Two-part models

A two-part model is used when an outcome is the product of a mix of two
processes, one for observed zeroes and non-zeroes, and one for an amount
among the non-zeroes. The expected outcome is the product of the two:

``` math
E(y) = \Pr(y > 0) \times E(y \mid y > 0)
```

With [`suest()`](https://tdmize.github.io/suest/reference/suest.md), fit
the two parts as separate models and combine them with `twopart`.
Predictions and marginal effects from `marginaleffects` are then on the
expected outcome, with standard errors from the joint covariance of the
two parts.

The examples use the medical expenditures of 19,386 adults in the 2004
Medical Expenditure Panel Survey. The first part is a logit of whether a
person had any medical expenditures. The second part is a Poisson model
of expenditures among the people who had any. With the robust standard
errors [`suest()`](https://tdmize.github.io/suest/reference/suest.md)
computes, this is Poisson pseudo-maximum likelihood: it fits a model of
the mean of any outcome that is zero or positive, not only counts, and
needs only the model of the mean to be right (Santos Silva and Tenreyro
2006). `family = quasipoisson` fits it without warnings about outcomes
that are not whole numbers.

## Fit the two parts and combine them

Both parts cluster on the dwelling unit (`duid`):

``` r

library(suest)
library(marginaleffects)
library(haven)

meps <- zap_labels(read_dta("https://tdmize.github.io/data/data/meps.dta"))
meps[c("female", "ins_unins", "anylim")] <-
  lapply(meps[c("female", "ins_unins", "anylim")], factor)
table(meps$exp_tot > 0)
#> 
#> FALSE  TRUE 
#>  3440 15946

any <- glm(I(exp_tot > 0) ~ age + female + lninc + pcs12 + mcs12 + ins_unins +
             anylim, family = binomial, data = meps)
amount <- glm(exp_tot ~ age + female + lninc + pcs12 + mcs12 + ins_unins +
                anylim, family = quasipoisson, data = meps, subset = exp_tot > 0)

fit <- suest(any, amount, model_names = c("Any", "Amount"), cluster = "duid",
             twopart = TRUE)
summary(fit)
#> 
#> Simultaneous results for Any, Amount                      Number of obs = 19,386
#>                                         (Std. err. adjusted for 10,940 clusters)
#> 
#> --------------------------------------------------------------------------------
#>               |                 Robust
#>               | Coefficient  std. err.        z   P>|z|     [95% conf. interval]
#> --------------+-----------------------------------------------------------------
#> Any           |
#>   (Intercept) |      -0.695      0.342   -2.033   0.042       -1.364      -0.025
#>           age |       0.030      0.002   19.358   0.000        0.027       0.033
#>       female1 |       0.957      0.044   21.744   0.000        0.871       1.044
#>         lninc |       0.349      0.026   13.273   0.000        0.298       0.401
#>         pcs12 |      -0.040      0.003  -12.130   0.000       -0.046      -0.033
#>         mcs12 |      -0.016      0.003   -6.430   0.000       -0.021      -0.011
#>    ins_unins1 |      -1.644      0.050  -32.946   0.000       -1.742      -1.546
#>       anylim1 |       0.799      0.073   10.991   0.000        0.657       0.942
#> --------------+-----------------------------------------------------------------
#> Amount        |
#>   (Intercept) |       9.183      0.284   32.360   0.000        8.627       9.739
#>           age |       0.011      0.001    8.146   0.000        0.008       0.014
#>       female1 |       0.046      0.039    1.180   0.238       -0.031       0.124
#>         lninc |       0.066      0.022    3.034   0.002        0.023       0.109
#>         pcs12 |      -0.036      0.002  -18.391   0.000       -0.040      -0.033
#>         mcs12 |      -0.012      0.002   -6.746   0.000       -0.015      -0.009
#>    ins_unins1 |      -0.910      0.064  -14.169   0.000       -1.036      -0.784
#>       anylim1 |       0.362      0.040    9.146   0.000        0.284       0.439
#> --------------------------------------------------------------------------------
```

With survey data, put the weights and clusters in a
[`svydesign()`](https://rdrr.io/pkg/survey/man/svydesign.html), fit the
first part to the full design and the second part to the people with any
expenditures, and give
[`suest()`](https://tdmize.github.io/suest/reference/suest.md) the full
design. `wts = TRUE` averages with the sampling weights:

``` r

design <- survey::svydesign(ids = ~duid, weights = ~wtdper, data = meps)
any_s <- survey::svyglm(I(exp_tot > 0) ~ age + female + ins_unins,
                        design = design, family = quasibinomial)
amount_s <- survey::svyglm(exp_tot ~ age + female + ins_unins,
                           design = subset(design, exp_tot > 0),
                           family = quasipoisson)
fit_s <- suest(any_s, amount_s, survey_design = design,
               observation_id = "dupersid", twopart = TRUE)
avg_slopes(fit_s, variables = "ins_unins", wts = TRUE)
#> 
#>  Estimate Std. Error     z Pr(>|z|)     S 2.5 % 97.5 %
#>     -2566        168 -15.3   <0.001 173.1 -2895  -2237
#> 
#> Term: ins_unins
#> Type: response
#> Comparison: 1 - 0
```

## Expected expenditures and effects on them

Predictions are the product of the two parts: the probability of any
expenditures from `Any` times the expected expenditures of someone who
has any from `Amount`:

``` r

avg_predictions(fit)
#> 
#>  Estimate Std. Error    z Pr(>|z|)   S 2.5 % 97.5 %
#>      3681         65 56.7   <0.001 Inf  3554   3809
#> 
#> Type: response
avg_slopes(fit, variables = c("age", "female", "ins_unins", "anylim"))
#> 
#>       Term Contrast Estimate Std. Error      z Pr(>|z|)     S   2.5 %  97.5 %
#>  age          dY/dX     48.6       4.96   9.80  < 0.001  72.9    38.9    58.3
#>  anylim       1 - 0   1532.8     151.64  10.11  < 0.001  77.4  1235.6  1830.0
#>  female       1 - 0    416.8     140.81   2.96  0.00308   8.3   140.8   692.7
#>  ins_unins    1 - 0  -2701.0     113.17 -23.87  < 0.001 415.8 -2922.8 -2479.2
#> 
#> Type: response
```

Expected expenditures average about \$3,680 a year. Being uninsured is
associated with \$2,700 less in expenditure; each year of age predicts a
\$49 increase; being female, a \$417 increase; and having any activity
limitation, a \$1,533 increase.

## Effects in each part

`type = "positive"` gives the effect on the probability of any
expenditures, and `type = "amount"` the effect on the expected
expenditures of someone who has any:

``` r

avg_comparisons(fit, variables = "ins_unins", type = "positive")
#> 
#>  Estimate Std. Error     z Pr(>|z|)     S  2.5 % 97.5 %
#>    -0.242    0.00855 -28.4   <0.001 585.7 -0.259 -0.226
#> 
#> Term: ins_unins
#> Type: positive
#> Comparison: 1 - 0
avg_comparisons(fit, variables = "ins_unins", type = "amount")
#> 
#>  Estimate Std. Error     z Pr(>|z|)     S 2.5 % 97.5 %
#>     -2603        136 -19.1   <0.001 267.9 -2870  -2336
#> 
#> Term: ins_unins
#> Type: amount
#> Comparison: 1 - 0
```

Being uninsured lowers the probability of any expenditures by 24
percentage points, and the expected expenditures of someone who has any
by about \$2,600. Both are averaged over the whole sample, so the two
effects are for the same people. `newdata = subset(exp_tot > 0)`
averages the amount effect over only the people with any expenditures,
the observations the second part was fit to:

``` r

avg_comparisons(fit, variables = "ins_unins", type = "amount",
                newdata = subset(exp_tot > 0))
#> 
#>  Estimate Std. Error     z Pr(>|z|)     S 2.5 % 97.5 %
#>     -2833        148 -19.2   <0.001 270.3 -3122  -2544
#> 
#> Term: ins_unins
#> Type: amount
#> Comparison: 1 - 0
```

Among the people with any expenditures, being uninsured is associated
with about \$2,830 less in expenditure. This describes the people who
spend, but it is averaged over different observations than the first
part’s effect, so the two can differ partly because the observations
used to calculate the marginal effect differ.

## How well the model predicts

The correlation between the outcome and its expected value summarizes
how well a two-part model predicts. It depends only on the expected
outcome, so it applies with any second part:

``` r

p <- predictions(fit)
cor(p$exp_tot, p$estimate)
#> [1] 0.3916091
```

The correlation is .392.

## Other second parts

The second part can also be a Gamma model with a log link
(`glm(family = Gamma(link = "log"))`), a linear model
([`lm()`](https://rdrr.io/r/stats/lm.html)), or a negative binomial
model ([`MASS::glm.nb()`](https://rdrr.io/pkg/MASS/man/glm.nb.html));
`twopart` works the same way with each. With a linear model of the log
of the outcome, add `retransform = "duan"` for Duan’s (1983) smearing
estimate or `retransform = "normal"` for the normal-theory
retransformation:

``` r

meps$ln_exp <- ifelse(meps$exp_tot > 0, log(meps$exp_tot), NA)
ln_amount <- lm(ln_exp ~ age + female + lninc + pcs12 + mcs12 + ins_unins +
                  anylim, data = meps, subset = exp_tot > 0)
fit_ln <- suest(any, ln_amount, model_names = c("Any", "LnAmount"),
                cluster = "duid", twopart = TRUE, retransform = "duan")
avg_slopes(fit_ln, variables = c("age", "ins_unins"))
#> 
#>       Term Contrast Estimate Std. Error     z Pr(>|z|)     S   2.5 % 97.5 %
#>  age          dY/dX     94.2       3.53  26.7   <0.001 517.9    87.2    101
#>  ins_unins    1 - 0  -3110.5      79.90 -38.9   <0.001   Inf -3267.1  -2954
#> 
#> Type: response
```

## Comparing two-part models

This example compares two second parts: a Gamma model with a log link,
and Poisson pseudo-maximum likelihood. Both share the first part, a
logit of whether a person had any expenditures. `twopart` names the two
pairs, each a first part and its second part; `Any` is in both:

``` r

any <- glm(I(exp_tot > 0) ~ age + female + lninc + ins_unins,
           family = binomial, data = meps)
gamma <- glm(exp_tot ~ age + female + lninc + ins_unins,
             family = Gamma(link = "log"), data = meps, subset = exp_tot > 0)
poisson <- glm(exp_tot ~ age + female + lninc + ins_unins,
               family = quasipoisson, data = meps, subset = exp_tot > 0)

fit2 <- suest(any, gamma, poisson, model_names = c("Any", "Gamma", "Poisson"),
              twopart = c("Any", "Gamma", "Any", "Poisson"))
effects2 <- avg_comparisons(fit2, variables = c("ins_unins", "age"))
effects2
#> 
#>       Term       Group Contrast Estimate Std. Error     z Pr(>|z|)     S   2.5 % 97.5 %
#>  age       Any+Gamma      +1         108       5.07  21.4   <0.001 335.1    98.5    118
#>  age       Any+Poisson    +1         100       4.47  22.4   <0.001 366.7    91.4    109
#>  ins_unins Any+Gamma      1 - 0    -2843     115.71 -24.6   <0.001 440.4 -3069.7  -2616
#>  ins_unins Any+Poisson    1 - 0    -2846     115.11 -24.7   <0.001 446.0 -3072.0  -2621
#> 
#> Type: response
hypotheses(effects2, hypothesis = difference ~ pairwise | term)
#> 
#>       Term                              Hypothesis Estimate Std. Error       z Pr(>|z|)    S 2.5 %
#>  age       (Any+Poisson +1) - (Any+Gamma +1)          -8.34       2.11 -3.9616   <0.001 13.7 -12.5
#>  ins_unins (Any+Poisson 1 - 0) - (Any+Gamma 1 - 0)    -3.52      36.81 -0.0956    0.924  0.1 -75.7
#>  97.5 %
#>   -4.21
#>   68.63
```

The two models give nearly the same effect of being uninsured: expected
expenditures about \$2,840 lower. The difference between them, about
\$4, is not significant. The effect of a year of age is about \$108 with
the Gamma second part and \$100 with the Poisson one, and that
difference, about \$8, is significant.

A two-part model can also be compared with a single model for the
outcome. Models not named in `twopart` keep their own prediction – here,
a Poisson model fit to everyone, zeros included:

``` r

one_part <- glm(exp_tot ~ age + female + lninc + ins_unins,
                family = quasipoisson, data = meps)

fit3 <- suest(any, poisson, one_part, model_names = c("Any", "Poisson", "OnePart"),
              twopart = c("Any", "Poisson"))
effects3 <- avg_comparisons(fit3, variables = c("ins_unins", "age"))
effects3
#> 
#>       Term       Group Contrast Estimate Std. Error     z Pr(>|z|)     S   2.5 % 97.5 %
#>  age       Any+Poisson    +1         100       4.47  22.4   <0.001 366.7    91.4    109
#>  age       OnePart        +1         101       4.60  22.0   <0.001 355.3    92.4    110
#>  ins_unins Any+Poisson    1 - 0    -2846     115.11 -24.7   <0.001 446.0 -3072.0  -2621
#>  ins_unins OnePart        1 - 0    -3028     109.03 -27.8   <0.001 561.6 -3241.9  -2815
#> 
#> Type: response
hypotheses(effects3, hypothesis = difference ~ pairwise | term)
#> 
#>       Term                            Hypothesis Estimate Std. Error      z Pr(>|z|)     S    2.5 %
#>  age       (OnePart +1) - (Any+Poisson +1)           1.29      0.305   4.22   <0.001  15.3    0.691
#>  ins_unins (OnePart 1 - 0) - (Any+Poisson 1 - 0)  -181.79     13.693 -13.28   <0.001 131.2 -208.626
#>   97.5 %
#>     1.89
#>  -154.95
```

The one-part model gives a larger effect of being uninsured, about
\$3,030, or about \$180 more than the two-part model, and the difference
is significant. The effects of age are close, about \$100 and \$101,
though the difference of about \$1 is significant.

## References

Belotti, F., Deb, P., Manning, W. G., & Norton, E. C. (2015). twopm:
Two-part models. *Stata Journal*, 15(1), 3-20.
<https://doi.org/10.1177/1536867X1501500102>

Duan, N. (1983). Smearing estimate: A nonparametric retransformation
method. *Journal of the American Statistical Association*, 78(383),
605-610. <https://doi.org/10.2307/2288126>

Santos Silva, J. M. C., & Tenreyro, S. (2006). The log of gravity.
*Review of Economics and Statistics*, 88(4), 641-658.
<https://doi.org/10.1162/rest.88.4.641>
