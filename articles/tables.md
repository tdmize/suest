# Regression tables

[`modelsummary`](https://modelsummary.com/) makes tables from a combined
object and from the `marginaleffects` results computed with it. Add
`shape = term ~ group` to put each model in its own column.

The models are the two logits of union membership from Getting started:

``` r

library(suest)
library(marginaleffects)
library(modelsummary)
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
