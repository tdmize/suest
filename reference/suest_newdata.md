# Stack the Estimation Samples from a SUEST Object

Creates model-specific new data using each component model's own
estimation sample. This is useful when effects should be averaged
separately within different samples.

## Usage

``` r
suest_newdata(object)
```

## Arguments

- object:

  A `suest_model` returned by
  [`suest()`](https://tdmize.github.io/suest/reference/suest.md).

## Value

A data frame containing the component model frames and internal
model-routing columns.
