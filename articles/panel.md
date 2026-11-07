# Panel models

## Random-effects logits

A 30% sample of households from the Health and Retirement Study, 2010 to
2020 (the households whose ID ends in 0, 1, or 2): about 31,500
interviews of 7,800 people aged 51 and older, interviewed every two
years. The two outcomes are depressive symptoms (four or more of the
eight CES-D items) and fair or poor self-rated health. Each is a
random-intercept logit fit with 25-point adaptive quadrature, with age
centered at 65.

``` r

library(suest)
library(marginaleffects)
library(haven)
library(lme4)

hrs <- zap_labels(read_dta("https://tdmize.github.io/data/data/hrs_si_2025.dta"))
hrs <- subset(hrs, year >= 2010 & age >= 51 & (hhidpn %/% 1000) %% 10 < 3)
hrs$depressed <- as.numeric(hrs$cesd >= 4)
hrs$poorhealth <- as.numeric(hrs$srhealthB == 0)
vars <- c("depressed", "poorhealth", "marstatB", "age", "woman", "race4cat", "collegeB")
hrs <- hrs[complete.cases(hrs[vars]), ]
hrs[c("marstatB", "woman", "race4cat", "collegeB")] <-
  lapply(hrs[c("marstatB", "woman", "race4cat", "collegeB")], factor)

dep <- glmer(depressed ~ marstatB + I(age - 65) + woman + race4cat + collegeB +
               (1 | hhidpn), family = binomial, data = hrs, nAGQ = 25,
             control = glmerControl(optimizer = "bobyqa"))
health <- glmer(poorhealth ~ marstatB + I(age - 65) + woman + race4cat + collegeB +
                  (1 | hhidpn), family = binomial, data = hrs, nAGQ = 25,
                control = glmerControl(optimizer = "bobyqa"))
fit_re <- suest(dep, health, model_names = c("Depressed", "Poor health"))
```

[`suest()`](https://tdmize.github.io/suest/reference/suest.md) clusters
the standard errors on the person. The marginal effects average over the
random intercept, so they are changes in the probability of each
outcome; age is changed by ten years:

``` r

effects_re <- avg_comparisons(fit_re,
  variables = list(marstatB = "reference", age = 10, woman = "reference",
                   race4cat = "reference", collegeB = "reference"),
  newdata = hrs)
effects_re
#> 
#>      Term       Group Contrast Estimate Std. Error      z Pr(>|z|)     S    2.5 %   97.5 %
#>  age      Depressed      +10    -0.0134    0.00266  -5.05  < 0.001  21.1 -0.01863 -0.00821
#>  age      Poor health    +10     0.0337    0.00365   9.25  < 0.001  65.2  0.02659  0.04090
#>  collegeB Depressed      1 - 0  -0.0629    0.00677  -9.29  < 0.001  65.9 -0.07617 -0.04964
#>  collegeB Poor health    1 - 0  -0.1636    0.00903 -18.11  < 0.001 241.1 -0.18132 -0.14591
#>  marstatB Depressed      1 - 0  -0.0956    0.00617 -15.50  < 0.001 177.5 -0.10773 -0.08354
#>  marstatB Poor health    1 - 0  -0.0689    0.00738  -9.35  < 0.001  66.6 -0.08340 -0.05448
#>  race4cat Depressed      1 - 0   0.0078    0.00779   1.00  0.31696   1.7 -0.00747  0.02307
#>  race4cat Depressed      2 - 0   0.0674    0.00990   6.81  < 0.001  36.6  0.04801  0.08680
#>  race4cat Depressed      3 - 0   0.0541    0.01697   3.19  0.00143   9.4  0.02084  0.08734
#>  race4cat Poor health    1 - 0   0.0935    0.01113   8.40  < 0.001  54.4  0.07169  0.11531
#>  race4cat Poor health    2 - 0   0.1979    0.01297  15.26  < 0.001 172.2  0.17248  0.22333
#>  race4cat Poor health    3 - 0   0.0742    0.02094   3.54  < 0.001  11.3  0.03313  0.11520
#>  woman    Depressed      1 - 0   0.0383    0.00624   6.14  < 0.001  30.1  0.02605  0.05051
#>  woman    Poor health    1 - 0  -0.0110    0.00839  -1.31  0.18897   2.4 -0.02747  0.00542
#> 
#> Type: response
hypotheses(effects_re, hypothesis = difference ~ revpairwise | term + contrast)
#> 
#>      Term Contrast                  Hypothesis Estimate Std. Error       z Pr(>|z|)     S   2.5 %
#>  age         +10   (Depressed) - (Poor health)  -0.0472    0.00374 -12.627   <0.001 119.0 -0.0545
#>  collegeB    1 - 0 (Depressed) - (Poor health)   0.1007    0.00858  11.734   <0.001 103.2  0.0839
#>  marstatB    1 - 0 (Depressed) - (Poor health)  -0.0267    0.00794  -3.361   <0.001  10.3 -0.0423
#>  race4cat    1 - 0 (Depressed) - (Poor health)  -0.0857    0.01066  -8.041   <0.001  50.0 -0.1066
#>  race4cat    2 - 0 (Depressed) - (Poor health)  -0.1305    0.01258 -10.371   <0.001  81.3 -0.1552
#>  race4cat    3 - 0 (Depressed) - (Poor health)  -0.0201    0.02124  -0.945    0.345   1.5 -0.0617
#>  woman       1 - 0 (Depressed) - (Poor health)   0.0493    0.00807   6.108   <0.001  29.9  0.0335
#>   97.5 %
#>  -0.0398
#>   0.1175
#>  -0.0111
#>  -0.0648
#>  -0.1058
#>   0.0216
#>   0.0651
```

**Example interpretation:** Being married is associated with a larger
reduction in the probability of depressive symptoms (9.6 percentage
points) than in the probability of fair or poor health (6.9 points); the
difference of 2.7 points has $`p < .001`$. Ten more years of age are
associated with a 1.3 point reduction in depressive symptoms but a 3.4
point increase in fair or poor health.

A random-effects model mixes differences between people with changes
within a person. The [within-person
effects](https://tdmize.github.io/suest/articles/within-person.md)
example separates them.

## Fixed-effects models

The National Longitudinal Survey of Young Women, 1968 to 1988: log wages
and hours worked, each regressed on total work experience, union
membership, their interaction, and year, with a fixed effect for each
woman.

``` r

library(plm)

nlswork <- zap_labels(read_dta("https://tdmize.github.io/data/data/nlswork.dta"))
nlswork[c("union", "year")] <- lapply(nlswork[c("union", "year")], factor)

wage_fe <- plm(ln_wage ~ ttl_exp * union + year, data = nlswork,
               index = c("idcode", "year"), model = "within")
hours_fe <- plm(hours ~ ttl_exp * union + year, data = nlswork,
                index = c("idcode", "year"), model = "within")
fit_fe <- suest(wage_fe, hours_fe, model_names = c("Log wage", "Hours"))
fit_fe
#> Seemingly Unrelated Estimation
#> Models: Log wage + Hours 
#> Model types: Log wage=panel_fe, Hours=panel_fe 
#> Model engines: Log wage=plm::plm, Hours=plm::plm 
#> Comparison scale: fitted values 
#> Clusters: 4,150 
#> Observations: Log wage=19,238, Hours=19,202 
#> Overlapping observations: 19,202 
#> Union observations: 19,238 
#> Parameters: 30
```

[`suest()`](https://tdmize.github.io/suest/reference/suest.md) clusters
the standard errors on the woman. The average effect of a year of work
experience on each outcome:

``` r

avg_slopes(fit_fe, variables = "ttl_exp")
#> 
#>     Group Estimate Std. Error    z Pr(>|z|)     S  2.5 % 97.5 %
#>  Log wage   0.0444    0.00236 18.8   <0.001 259.5 0.0397  0.049
#>  Hours      0.7115    0.07260  9.8   <0.001  72.9 0.5692  0.854
#> 
#> Term: ttl_exp
#> Type: response
#> Comparison: dY/dX
```

Each additional year of work experience is associated with 0.044 higher
log wages (about 4.5% higher wages) and 0.7 more hours of work per week.
