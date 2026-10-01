# Separate samples or groups

Example 6.6 of Mize, Doan, and Long (2019). The [same example in
Stata](https://www.trentonmize.com/software/mecompare/examples/groups.html)
uses `mecompare`.

This example compares effects from the same binary logit model fitted to
two different samples, one from 1986 and one from 2016. The same
approach works for comparing groups, such as men and women. The samples
are selected with `subset`, and
[`suest_newdata()`](https://tdmize.github.io/suest/reference/suest_newdata.md)
makes sure each marginal effect is averaged over its own model’s sample.

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

vars66 <- c("helpsickB", "polviews", "conserv", "faminc", "employed",
            "woman", "age", "college", "married", "parent", "race", "year")
dat66 <- gss[complete.cases(gss[vars66]), ]
dat66 <- factorize(dat66, c("conserv", "employed", "woman", "college",
                            "married", "parent", "race"))

model1986 <- glm(helpsickB ~ conserv + faminc + employed + woman + age +
                   college + married + parent + race, family = binomial("logit"),
                 data = dat66, subset = year == 1986)
model2016 <- glm(helpsickB ~ conserv + faminc + employed + woman + age +
                   college + married + parent + race, family = binomial("logit"),
                 data = dat66, subset = year == 2016)
stopifnot(nobs(model1986) == 1254, nobs(model2016) == 1670)
c(`1986` = nobs(model1986), `2016` = nobs(model2016))
#> 1986 2016 
#> 1254 1670
fit66 <- suest(model1986, model2016, model_names = c("1986", "2016"))
fit66
#> Seemingly Unrelated Estimation
#> Models: 1986 + 2016 
#> Model types: 1986=logit, 2016=logit 
#> Model engines: 1986=stats::glm, 2016=stats::glm 
#> Comparison scale: predicted probabilities 
#> Observations: 1986=1,254, 2016=1,670 
#> Overlapping observations: 0 
#> Union observations: 2,924 
#> Parameters: 22

effects66 <- avg_comparisons(fit66, variables = "conserv",
                             newdata = suest_newdata(fit66))
effects66
#> 
#>  Group Estimate Std. Error     z Pr(>|z|)    S  2.5 %  97.5 %
#>   1986  -0.0917     0.0296  -3.1  0.00193  9.0 -0.150 -0.0337
#>   2016  -0.2582     0.0253 -10.2  < 0.001 78.7 -0.308 -0.2085
#> 
#> Term: conserv
#> Type: response
#> Comparison: 1 - 0
hypotheses(effects66, hypothesis = difference ~ revpairwise)
#> 
#>       Hypothesis Estimate Std. Error    z Pr(>|z|)    S  2.5 % 97.5 %
#>  (1986) - (2016)    0.166     0.0389 4.28   <0.001 15.7 0.0902  0.243
```

**Example interpretation:** Conservatives had significantly lower
predicted probabilities of saying that government should be responsible
for providing health care in both 1986 (marginal effect = -0.092) and
2016 (marginal effect = -0.258). The gap between conservatives and
nonconservatives increased over time: the conservative effect is
significantly more negative in 2016 than in 1986 (cross-model difference
= 0.166, $`p < .001`$).
