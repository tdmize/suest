# Pool SUEST systems across multiple imputations

`suest_mi()` pools complete joint SUEST coefficient vectors and
covariance matrices using Rubin's rules. Supply either a list whose
elements are compatible
[`suest()`](https://tdmize.github.io/suest/reference/suest.md) results
fitted in separate imputed datasets, or the
[`mice::mira`](https://amices.org/mice/reference/mira.html) object
returned by [`with()`](https://rdrr.io/r/base/with.html) when that
expression fits a SUEST system. Result lists created by
[`mitools::with.imputationList()`](https://rdrr.io/pkg/mitools/man/with.imputationList.html)
are also accepted directly.

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

This initial interface supports coefficients, covariance matrices,
summaries, and coefficient-level hypotheses. Predictions, slopes, and
comparisons must be estimated within each imputation and pooled as
estimands; they are not yet implemented for `"suest_mi"` objects.
