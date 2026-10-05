# Pool SUEST systems across multiple imputations

`suest_mi()` combines
[`suest()`](https://tdmize.github.io/suest/reference/suest.md) results
fitted in each of several imputed datasets, applying Rubin's rules to
the full set of coefficients and their joint covariance matrix. Supply a
list of compatible
[`suest()`](https://tdmize.github.io/suest/reference/suest.md) results,
the [`mice::mira`](https://amices.org/mice/reference/mira.html) object
returned by [`with()`](https://rdrr.io/r/base/with.html) when it fits a
SUEST system, or the result list from
[`mitools::with.imputationList()`](https://rdrr.io/pkg/mitools/man/with.imputationList.html).

## Usage

``` r
suest_mi(fits)

# S3 method for class 'suest_mi'
coef(object, ...)

# S3 method for class 'suest_mi'
vcov(object, ...)

# S3 method for class 'suest_mi'
nobs(object, ...)

# S3 method for class 'suest_mi'
summary(object, conf.level = 0.95, ...)

# S3 method for class 'suest_mi'
print(x, ...)
```

## Arguments

- fits:

  A list containing at least two compatible `"suest_model"` objects, one
  per imputation; a `"mira"` object containing those fits; or a result
  list returned by
  [`mitools::with.imputationList()`](https://rdrr.io/pkg/mitools/man/with.imputationList.html).

- object, x:

  A `"suest_mi"` object.

- ...:

  Additional arguments; currently ignored.

- conf.level:

  Confidence level for `summary.suest_mi()`.

## Value

A `"suest_mi"` object containing pooled coefficients, total,
within-imputation, and between-imputation covariance matrices, Rubin
large-sample degrees of freedom, and the original SUEST fits.

## Details

The pooled covariance is `Ubar + (1 + 1 / m) B`, where `Ubar` is the
mean within-imputation joint covariance, `B` is the between-imputation
covariance of the complete joint coefficient vector, and `m` is the
number of imputations. Because each within-imputation input is a
complete SUEST system, cross-model covariance is retained in both
components.

A `"mira"` object from
[`mice::with()`](https://amices.org/mice/reference/with.mids.html) is
read through its public `analyses` component. The list returned by
[`mice::getfit()`](https://amices.org/mice/reference/getfit.html) is
accepted as well. The `mice` package is not required when pooling an
ordinary list.

Results from
[`mitools::with.imputationList()`](https://rdrr.io/pkg/mitools/man/with.imputationList.html)
use the same list-based route, with their originating call retained in
the returned object. Legacy `"imputationResultList"` wrappers are
accepted as well. The `mitools` package is not required for an
already-created result list.

## Marginal effects

Predictions, slopes, comparisons, and hypothesis tests from
`marginaleffects` are estimated within each imputation, using that
imputation's joint SUEST covariance, and then pooled with Rubin's rules,
as Stata's `mimrgns` does. Cross-model differences, such as
`hypothesis = difference ~ revpairwise`, are pooled the same way. This
works on the `"mira"` object from `with(imp, suest(...))` and equally on
a `"suest_mi"` object, so the two give identical results:

    analyses <- with(imp, suest(glm(y ~ x, family = binomial),
                                glm(y ~ x + z, family = binomial)))
    marginaleffects::avg_comparisons(analyses, variables = "x",
      hypothesis = difference ~ revpairwise)
    marginaleffects::avg_comparisons(suest_mi(analyses), variables = "x",
      hypothesis = difference ~ revpairwise)

Without `newdata`, each imputation's completed data are used when the
imputations came from `mice`. For a plain list or `mitools` results,
each imputation's own estimation sample is used instead, and
`marginaleffects` warns that it could not recover the original data from
a `mids` object; that warning does not apply to
[`suest()`](https://tdmize.github.io/suest/reference/suest.md) results.
`mice` warns "Large sample assumed." because SUEST, like Stata's
`suest`, reports large-sample (z) statistics; the pooled degrees of
freedom are then Rubin's large-sample values, as in Stata's
`mi estimate`.
