# Logit random slopes: disjoint validation and remaining raw-Stata diagnostic

Date: 2026-09-25. Development version: 0.1.5.9000.

Independent R validation passes for the final, disjoint sample pattern.
All three returned full fits converged after one iteration. Parameters, native
component covariance, and integrated probabilities/slopes agree closely.
The R joint sandwich meets the existing 1e-6 absolute covariance and 1e-4
relative-SE bounds. Strict raw GSEM covariance parity is not established.
A separate reconstruction of Stata's covariance narrowly misses its original
absolute bound and remains an open diagnostic. No production R code changed.

## Returned data and direct comparisons

Source: `suest_r_glmmtmb_logit_rs_disjoint_v1.log`, revision 1, with
`LOGIT_RS_DISJOINT_COMPLETE=1`. Each component uses 600 observations in 50
panels. The samples do not overlap; the union contains 1,200 observations and
100 panels. Stata's joint correction is 100/99 and R's per-component correction
is 50/49. Normalize Stata covariance by (50/49)/(100/99) = 99/98 for comparison.

| Fit | Maximum absolute coefficient difference | Maximum absolute covariance difference | Maximum relative SE difference | Absolute log-likelihood difference |
|---|---:|---:|---:|---:|
| y1_left | 1.211e-8 | 4.027e-7 | 7.552e-7 | 8.804e-8 |
| y2_right | 5.659e-6 | 3.115e-7 | 3.449e-6 | 3.339e-8 |
| Joint raw | 5.741e-6 | 4.462e-3 | 1.0557e-2 | 5.746e-8 |
| Joint with aligned correction | 5.741e-6 | 2.719e-3 | 1.5699e-2 | 5.746e-8 |

Coordinates are fixed effects, both log SDs, and Fisher correlation. The raw
joint SE difference after correction alignment is approximately 1.57%.

## Independent covariance gate: passed

The shared independent bivariate Laplace evaluator and coordinate transformations
are unchanged. At the R point, native-coordinate independent covariance differs
from R by at most 8.820e-7, with relative SE difference 1.553e-6 (0.0001553%).
Both native-coordinate and pure Fisher-coordinate sandwiches pass the original
bounds. Native component covariance audits differ by at most 4.001e-7 in
covariance and 7.503e-7 in relative SE. Independent score cross-products and
sandwich cross-equation blocks are exactly zero, as required for disjoint panels.
Finite-difference score sensitivity is below 5.330e-10 and curvature sensitivity
below 1.799e-8; independent likelihoods differ from the joint Stata result by
less than 5.839e-8.

## Raw Stata diagnostics: preserve the limitation

Stata's model-based covariance has a maximum cross-equation entry of 0.000840424,
although the independent likelihood curvature is block diagonal. Its raw robust
covariance has a maximum cross-block entry of 0.002691018. Recovered Stata score
cross-products have cross-block entries below 8.21e-12, versus exactly zero
independently. These observations locate the main raw discrepancy in numerical
curvature, rather than sample overlap or cluster correction.

Reconstructing raw covariance with Stata's supplied model-based covariance and
independent scores at Stata's endpoint gives maximum absolute error
**1.052009770e-6**, exceeding the existing **1e-6** bound. Only one entry exceeds
that bound: `Y1::atanh_rho` variance, 0.27954976440501 returned versus
0.27954871239524 reconstructed. Its scaled difference is 3.76323e-6; the largest
relative SE difference over all parameters is 4.531e-6, which meets the existing
1e-4 relative bound. Centering the scores changes the reconstruction by only
about 1.25e-13 and does not explain the discrepancy. The remaining difference
is consistent with numerical score derivatives in the returned Stata fit;
its exact cause is not certified without further Stata derivative output.

The first combined audit stopped on this raw-reconstruction assertion; its log
is preserved. The final audit explicitly separates the independent-R pass from
raw-reconstruction status, keeps the original thresholds, and records
`raw_reconstruction_absolute_pass = FALSE`. It does not relabel this diagnostic
as a pass. Regression tests retain the independent-R bounds and verify the
recorded raw limitation. Balanced and partial raw-reconstruction checks are
unchanged. This is a documented change in gate scope, not a claim that every
previous assertion now passes. A package-test pass therefore does not establish
strict raw-Stata parity.

## Native integrated probabilities and slopes

Independent adaptive integration and analytic gradients use each Stata
component's endpoint, its own 600-row estimation sample, and its full covariance
including nuisance parameters. Maximum probability difference is 8.229e-11;
maximum slope difference is 3.060e-10. Maximum relative SE difference is
8.310e-8. Existing probability and SE bounds pass without changes.

## Status and next work

Balanced, partial, and disjoint independent R validation is complete for the
restricted Bernoulli logit correlated random-intercept/numeric-slope route.
The raw Stata discrepancies and the reconstruction diagnostic above remain
recorded limitations. No repeat Stata run is requested for this R checkpoint.
If exact raw-reconstruction parity becomes a goal, collect scores and derivative
settings from the disjoint joint fit in a targeted diagnostic run; retain the
original log and fixture. Choose the next model increment separately.

The checkpoint includes returned logs, separate raw/audit fixtures, comparison
tables, reproduction scripts, source package, full Git history, and verification
logs. Source provenance is recorded in each raw fixture. Review confirmed the
sample geometry, correction alignment, independent gate, and the need to keep
the failed raw diagnostic explicit. Full package results are recorded below.

## Final package verification

- Focused balanced/partial/disjoint regression: 90 assertions passed, with no
  failures, warnings, errors, or skips.
- Full source-package suite: 1,626 passed, 0 failed, 0 skipped, and the same
  three existing test warnings.
- `R CMD build` succeeded; `R CMD check --no-manual` returned **Status: OK**,
  including examples and vignette rebuilding.
- Read-only review found no blocking issues in the final changes. It did not
  certify strict raw-Stata parity or the exact derivative-error mechanism.
- Environment: R 4.5.2, glmmTMB 1.1.14, TMB 1.9.25, marginaleffects 1.0.0.

The passing package suite checks the independent-R result and accurately records
the failed raw diagnostic; it does not resolve that diagnostic. Full logs and
the original combined-audit failure are retained in the checkpoint. Nothing
was pushed or released.
