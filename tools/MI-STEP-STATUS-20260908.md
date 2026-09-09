# Initial multiple-imputation step

## Implemented

- `suest_mi()` accepts two or more compatible `suest_model` objects, one per
  imputation.
- It also accepts a `mice::mira` object returned by `with()` and a
  `mice::getfit()` analysis list directly.
- Result lists returned by `mitools::with.imputationList()` and legacy
  `imputationResultList` wrappers are accepted directly.
- Rubin pooling is applied to the complete joint coefficient vector.
- Within- and between-imputation cross-model covariance is preserved.
- `coef()`, `vcov()`, `nobs()`, `summary()`, and `print()` methods are included.
- `marginaleffects::hypotheses()` provides pooled coefficient contrasts and
  standard errors. `summary()` uses Rubin large-sample degrees of freedom for
  individual coefficients.
- Strict validation rejects incompatible model structures, coefficient layouts,
  sample counts, cluster counts, and nonfinite results.

## Deliberately deferred

Predictions, slopes, and comparisons are not evaluated at averaged
coefficients. They stop explicitly because correct nonlinear MI inference
requires estimating each estimand within every imputation and then pooling the
estimand and its covariance.

## Verification

- Focused MI tests reproduce independent Rubin calculations exactly.
- A real `mice::mice()` -> `with()` -> `mira` workflow and
  `mice::getfit()` extraction pass end to end. The always-run structural test
  follows the current `mira` public contract.
- A real `mitools::imputationList()` workflow matches
  `mitools::MIcombine()` exactly for the complete joint coefficient and
  covariance system.
- Full unit suite: 706 passed, 0 failed, 3 expected fitting warnings.
- Package build/check: all code, installation, namespace, documentation,
  examples, and tests pass. The two warnings are the existing unavailable
  vignette-builder warnings.
- GitHub remains unchanged.
