# Survey binary-response step status — 2026-09-08

## Scope

This increment extends the existing one-stage common-design survey route to
binary-response `survey::svyglm()` fits using `quasibinomial("logit")` and
`quasibinomial("probit")`. It retains the established requirements:

- one ordinary one-stage `survey::svydesign()` supplied through
  `survey_design =`;
- explicit observation-ID columns;
- coefficient parameters only;
- full sampled design retained before model-specific domains/subsets;
- stratification and optional first-stage FPCs;
- zero padding of model influences outside each estimation sample;
- exact reproduction of every native `svyglm()` covariance block before the
  joint covariance is constructed;
- asymptotic-normal package inference by default, with survey design df retained
  only as metadata.

Gaussian and binary survey models are not mixed in the same system in this
increment. Replicate weights, multistage/two-phase designs, PPS,
calibration/raking/post-stratification, ancillary survey parameters, extra
fitting weights, and automatic survey t/F inference remain unsupported.

## Binary influence

For observation i, with model matrix row x_i, fitting weight w_i, mean mu_i,
link derivative mu'_i, and GLM variance V(mu_i), the coefficient influence is

`[x_i w_i (y_i-mu_i) mu'_i / V(mu_i)] A^{-1}`

where

`A = sum_i x_i x_i' w_i (mu'_i)^2 / V(mu_i)`.

For logit the score factor reduces to `w_i (y_i-mu_i)`. For probit it is
`w_i (y_i-mu_i) phi(eta_i) / [mu_i(1-mu_i)]`. A common constant rescaling of
prior weights cancels from the influence. Quasi-family dispersion is not
multiplied into the native design-based coefficient covariance.

The implementation reconstructs these influences directly rather than copying
the optional influence attribute returned by `survey`. To reproduce native
covariance exactly at the actual IRLS stopping state, it uses the retained final
working residuals, working weights, and `model$naive.cov`. These are the fitted-
state representation of the formula above; recomputing them at the printed
coefficient vector can differ slightly at the GLM convergence tolerance.

## Verification status

Focused tests cover native influence reproduction, exact native covariance
blocks, identical and partially overlapping samples, disjoint domains, FPCs,
factors/interactions, formula offsets, response/link predictions, and common
weight rescaling.

The returned Stata 19.5 gate confirms full survey-logit numerical parity.
Survey-probit point estimates agree, but the covariance has a documented
finite-sample engine difference: R `svyglm()` uses expected/Fisher information
for the sandwich bread, while the returned Stata covariance is reproduced by
an observed-information bread. The R route retains the native `svyglm()`
convention. See
`SURVEY-BINARY-BENCHMARK-RESULTS-20260908.md`.

## Final validation

Final validation completed under R 4.3.3 with `survey` 4.5 and
`marginaleffects` 0.31.0:

- focused survey files: 132 expectations passed;
- complete unit suite: 147 tests / 872 expectations passed, with no failures,
  errors, or skips and three expected model-fitting warnings from stress cases;
- numerical acceptance: 144 passed, 0 failed;
- standard source build with vignette creation: passed;
- clean installation and load: passed;
- `R CMD check --no-manual`: `OK` (0 errors, 0 warnings, 0 notes).

The increment is certified within the stated scope. See
`SURVEY-BINARY-FINAL-VALIDATION-RESULTS-20260909.md`.
