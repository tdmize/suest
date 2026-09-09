# Package index

## Core workflow

Combine models, pool imputed systems, and construct model-specific
evaluation data

- [`suest()`](https://tdmize.github.io/suest/reference/suest.md)
  [`coef(`*`<suest_model>`*`)`](https://tdmize.github.io/suest/reference/suest.md)
  [`vcov(`*`<suest_model>`*`)`](https://tdmize.github.io/suest/reference/suest.md)
  [`nobs(`*`<suest_model>`*`)`](https://tdmize.github.io/suest/reference/suest.md)
  [`print(`*`<suest_model>`*`)`](https://tdmize.github.io/suest/reference/suest.md)
  : Combine two fitted models with seemingly unrelated estimation
- [`suest_mi()`](https://tdmize.github.io/suest/reference/suest_mi.md)
  [`coef(`*`<suest_mi>`*`)`](https://tdmize.github.io/suest/reference/suest_mi.md)
  [`vcov(`*`<suest_mi>`*`)`](https://tdmize.github.io/suest/reference/suest_mi.md)
  [`nobs(`*`<suest_mi>`*`)`](https://tdmize.github.io/suest/reference/suest_mi.md)
  [`summary(`*`<suest_mi>`*`)`](https://tdmize.github.io/suest/reference/suest_mi.md)
  [`print(`*`<suest_mi>`*`)`](https://tdmize.github.io/suest/reference/suest_mi.md)
  : Pool SUEST systems across multiple imputations
- [`suest_newdata()`](https://tdmize.github.io/suest/reference/suest_newdata.md)
  : Stack the estimation samples from a SUEST object
