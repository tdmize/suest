# Survey binary-response final R validation — 2026-09-09

Package: `suest` 0.1.4 certified candidate  
Runtime: R 4.3.3 on x86_64 Ubuntu 24.04.3 LTS  
Key versions in the original certification: `survey` 4.5,
`marginaleffects` 0.31.0, `sandwich` 3.1-0, `MASS` 7.3-60.0.1

## Targeted correction

The original analytic binary influence formula was correct, but recomputing its
working quantities from the reported coefficient vector did not always match
the final `svyglm()` IRLS state at floating-point tolerance. The adapter now
reconstructs the same influence from the model's retained final working
residuals, working weights, and `naive.cov`. This reproduces the optional native
influence attribute exactly when present while remaining independent of that
optional attribute, and preserves exact native `svyglm()` covariance blocks.

The survey weight-rescaling tests now compare separately fitted models with
tight absolute numerical tolerances after requesting `epsilon = 1e-12` and
`maxit = 100`. This recognizes small IRLS stopping-state differences without
weakening the exact native-block checks.

## Results

| Gate | Result |
|---|---|
| `tests/testthat/test-survey.R` | 110 expectations passed |
| `tests/testthat/test-stata-survey-binary-reference.R` | 22 expectations passed |
| Focused survey total | 132 passed, 0 failed |
| Complete `tests/testthat` suite | 147 tests; 872 expectations; 0 failures, 0 errors, 0 skipped |
| Numerical acceptance suite | 144 passed, 0 failed |
| `R CMD build suest` | Passed; vignette created; `suest_0.1.4.tar.gz` built |
| Clean `R CMD INSTALL` and `library(suest)` | Passed; version 0.1.4 loaded from the fresh library |
| `R CMD check --no-manual suest_0.1.4.tar.gz` | `Status: OK`; 0 errors, 0 warnings, 0 notes |

The full unit suite emitted three expected warnings from deliberate fitting
stress cases: one non-convergence warning and two fitted-probability warnings.
They were test output, not test failures or `R CMD check` warnings.

## Dependency qualification corrected during release preflight

The earlier `marginaleffects` 1.0.0 `%||%` result was reproduced as an
incomplete-install problem. After completing its dependencies in a clean
library, version 1.0.0 passed a plain `lm()` prediction, all 132 focused survey
expectations, the complete 147-test unit suite, and all 144 numerical acceptance
cases. A single test assertion now strips the new `jacobian_source` attribute
from a scalar standard error before numeric comparison; the estimate and
standard error were unchanged. No `suest` calculation or dependency source was
patched.

## Cross-language interpretation retained

Survey logit retains full R/Stata numerical parity. Survey probit coefficients
agree, while finite-sample covariance differs because native `svyglm()` uses
expected/Fisher information and the returned Stata VCE is reproduced by
observed information. The package preserves exact native R covariance and does
not add a Stata-compatibility override.

No Stata source, Stata benchmark input, Git commit, Git push, or GitHub state was
changed during final validation.
