# GitHub release preparation — suest 0.1.4

## Release status

Version 0.1.4 is release-ready based on the local 2026-09-09 release-candidate
run. No GitHub commit, push, tag, or release was created during preparation.

Runtime and key dependencies:

- R 4.3.3 on Ubuntu 24.04.3 LTS;
- marginaleffects 1.0.0;
- survey 4.5;
- insight 1.5.4;
- data.table 1.18.6.1;
- sandwich 3.1.3.

Validation results:

| Gate | Result |
|---|---|
| Focused survey tests | 132 passed, 0 failed |
| Complete unit suite | 147 tests; 869 passed assertions plus 3 expected fitting warnings; 0 failures, 0 errors, 0 skipped |
| Numerical acceptance suite | 144 passed, 0 failed |
| `R CMD build` | Passed; vignette created; `suest_0.1.4.tar.gz` built |
| Clean installation/load | Passed; version 0.1.4 loaded from the fresh library |
| `R CMD check --no-manual` | `Status: OK`; 0 errors, 0 warnings, 0 notes |

The three unit-test warnings are deliberate fitting stress cases: one
non-convergence warning and two fitted-probability warnings. They are not
package-check warnings.

## Compatibility correction

The earlier report that marginaleffects 1.0.0 failed at `%||%` was caused by an
incomplete dependency installation. With a complete clean library, a plain
`lm()` prediction and every `suest` release gate above pass. One test assertion
now compares the numeric value of a contrast standard error while ignoring the
new `jacobian_source` metadata attribute; no package calculation was changed.

## Suggested release title

`suest 0.1.4 — survey binary models, MI pooling, GEE, and expanded adapters`

## Suggested release notes

Version 0.1.4 adds a narrow, explicitly scoped one-stage survey route for
linear, logit, and probit `survey::svyglm()` models, including domain overlap,
factors, offsets, and first-stage FPCs. Survey logit has full R/Stata numerical
parity. Survey probit coefficients agree across engines, while the R route
deliberately preserves native `svyglm()` expected-information covariance rather
than forcing Stata's observed-information finite-sample VCE convention.

The release also adds initial multiple-imputation pooling, expanded pweight and
cluster support, GEE models, panel FE/BE/RE and Gaussian random-intercept ML
routes, bivariate and IV probit, additional censored/truncated/count adapters,
cross-language fixtures, and a substantially expanded numerical test suite.
See `NEWS.md` and `tools/STATA-PARITY.md` for the complete contracts and
qualifications.

Release validation passed 132 focused survey assertions, 147 complete unit
tests, and all 144 numerical acceptance cases. The package built with its
vignette, installed and loaded from a clean library, and completed
`R CMD check --no-manual` with 0 errors, warnings, or notes.

## Post-release priorities

1. nonlinear panel model families, one independently benchmarked family at a
   time;
2. general GLMM/multilevel support, beginning with a narrow two-level
   random-intercept binary family;
3. MI estimand pooling and other extensions at lower priority.

The detailed scope is in `ROADMAP-POST-0.1.4-20260909.md`.

## Publication boundary

Before publication, review the local Git diff, commit the 0.1.4 source, push
`main`, allow the R-CMD-check and numerical-acceptance workflows to complete,
then create annotated tag `v0.1.4` and the GitHub release using the built source
tarball. Do not publish if either workflow is not green.
