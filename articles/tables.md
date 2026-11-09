# Regression tables and coefficient plots

[`modelsummary`](https://modelsummary.com/) makes tables from a combined
object and from the `marginaleffects` results computed with it. Add
`shape = term ~ group` to put each model in its own column.
[`modelplot()`](https://modelsummary.com/man/modelplot.html) makes
coefficient plots and takes the same arguments.

The models are the two logits of union membership from Getting started:

``` r

library(suest)
library(marginaleffects)
library(modelsummary)
library(ggplot2)
library(haven)

nlsw <- zap_labels(read_dta("https://tdmize.github.io/data/data/nlsw88.dta"))
nlsw <- nlsw[!is.na(nlsw$union), ]
nlsw[c("collgrad", "south")] <- lapply(nlsw[c("collgrad", "south")], factor)

base <- glm(union ~ collgrad + age + south, family = binomial("logit"),
            data = nlsw)
adjusted <- glm(union ~ collgrad + age + south + wage + ttl_exp,
                family = binomial("logit"), data = nlsw)
combined <- suest(base, adjusted, model_names = c("Base", "Adjusted"))
```

## Coefficients

The standard errors are each model’s robust standard errors from the
joint covariance matrix.

``` r

modelsummary(combined, shape = term ~ group)
```

|             | \(1\)   |          |
|-------------|---------|----------|
|             | Base    | Adjusted |
| (Intercept) | -1.340  | -1.761   |
|             | (0.704) | (0.719)  |
| collgrad1   | 0.507   | 0.311    |
|             | (0.119) | (0.132)  |
| age         | 0.009   | 0.008    |
|             | (0.018) | (0.018)  |
| south1      | -0.747  | -0.677   |
|             | (0.116) | (0.118)  |
| wage        |         | 0.054    |
|             |         | (0.014)  |
| ttl_exp     |         | 0.006    |
|             |         | (0.013)  |
| Num.Obs.    | 1878    |          |

The same coefficients as a plot, without the intercepts:

``` r

modelplot(combined, shape = term ~ group, coef_omit = "Intercept") +
  scale_color_discrete(labels = c("Base", "Adjusted"))
```

![Logit coefficients and 95% confidence intervals for union membership
from the base and adjusted
models.](tables_files/figure-html/coefficients-plot-1.png)

## Marginal effects and their differences

``` r

effects <- avg_slopes(combined, variables = c("collgrad", "age", "south"))
differences <- hypotheses(effects, hypothesis = difference ~ revpairwise | term)
differences$group <- "Base - Adjusted"

modelsummary(
  list("Average marginal effects" = effects, "Difference" = differences),
  shape = term ~ group
)
```

|          | Average marginal effects |          | Difference      |
|----------|--------------------------|----------|-----------------|
|          | Base                     | Adjusted | Base - Adjusted |
| age      | 0.002                    | 0.001    | 0.000           |
|          | (0.003)                  | (0.003)  | (0.001)         |
| collgrad | 0.097                    | 0.057    | 0.039           |
|          | (0.024)                  | (0.025)  | (0.009)         |
| south    | -0.131                   | -0.117   | -0.013          |
|          | (0.019)                  | (0.020)  | (0.004)         |
| Num.Obs. | 1878                     |          | 1878            |

The same results as a plot:

``` r

modelplot(
  list("Average marginal effects" = effects, "Difference" = differences),
  shape = term ~ group,
  background = list(geom_vline(xintercept = 0, color = "grey50"))
) +
  labs(x = "Change in the probability of union membership")
```

![Average marginal effects on the probability of union membership from
the base and adjusted models, and their differences, with 95% confidence
intervals.](tables_files/figure-html/effects-plot-1.png)
