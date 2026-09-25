# suest 0.1.6 release notes

Prepared 2026-09-25 UTC. This release collects the NB2 random-intercept and
Poisson, NB2, and Bernoulli-logit random-slope development increments after
0.1.5. Existing survey, MI, panel, and single-level model contracts are unchanged.

## Added models

| Route | Supported structure | Joint nuisance parameters |
|---|---|---|
| NB2 log random intercept | One Gaussian intercept, constant dispersion | Log size and log intercept SD |
| Poisson log random slope | One correlated Gaussian intercept/numeric slope | Both log SDs and Fisher correlation |
| NB2 log random slope | Same, with constant dispersion | Log size, both log SDs, Fisher correlation |
| Bernoulli logit random slope | Same; individual binary outcomes | Both log SDs and Fisher correlation |

These routes use unweighted `glmmTMB` fits, version 1.1.14 or later. They exclude
weights, offsets, zero inflation, parameter constraints, extra grouping factors,
and unsupported covariance structures. Population predictions integrate over
the random effects. `marginaleffects` inference includes nuisance covariance.
Use well-scaled predictors. For numerical logit slope SEs, use centered
coefficient differences and check sensitivity to the step size.

## Validation and limits

Balanced, partial-overlap, and disjoint samples have returned Stata benchmarks
and independent R checks. Cluster correction and likelihood approximation
choices are explicit. Raw Stata matrices are retained separately from
independent references. Strict raw GSEM covariance parity is not established.

Raw random-slope joint SE differences reach approximately 1% for Poisson,
1.16% for NB2, and 4.23% for logit. Independent audits explain the principal
differences through numerical curvature. The logit disjoint raw reconstruction
misses one original absolute bound, 1.052e-6 against 1e-6, while its independent
R covariance audit passes. This unresolved diagnostic remains documented; it
was not fixed by relaxing a numerical threshold. Package-test success does not
claim that every raw-Stata diagnostic passes.

Detailed evidence:

- [GLMMTMB-NBINOM2-RI-VALIDATION-20260923.md](https://github.com/tdmize/suest/blob/main/tools/GLMMTMB-NBINOM2-RI-VALIDATION-20260923.md)
- [GLMMTMB-POISSON-RS-VALIDATION-20260923.md](https://github.com/tdmize/suest/blob/main/tools/GLMMTMB-POISSON-RS-VALIDATION-20260923.md)
- [GLMMTMB-NBINOM2-RS-VALIDATION-20260924.md](https://github.com/tdmize/suest/blob/main/tools/GLMMTMB-NBINOM2-RS-VALIDATION-20260924.md)
- [GLMMTMB-LOGIT-RS-VALIDATION-20260924.md](https://github.com/tdmize/suest/blob/main/tools/GLMMTMB-LOGIT-RS-VALIDATION-20260924.md)
- [GLMMTMB-LOGIT-RS-DISJOINT-VALIDATION-20260925.md](https://github.com/tdmize/suest/blob/main/tools/GLMMTMB-LOGIT-RS-DISJOINT-VALIDATION-20260925.md)

The release changes metadata and documentation only relative to the validated
1e58447404b9d7e793525fc9b191392579ea541e development checkpoint. Final package,
website, and publication status are recorded in the current session handoff.
