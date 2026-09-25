# Poisson random intercept and slope: completed validation checkpoint

The restricted glmmTMB Poisson-log random-slope increment is implemented and
validated for balanced, partially overlapping/higher-clustered, and disjoint
samples. No further Stata run is required for this increment. This is still
unpublished development version 0.1.5.9000; nothing has been pushed to GitHub.

Validation uses returned Stata component information and recovered joint score cross-products,
plus an independent likelihood/curvature calculation. It does **not** assert
that raw joint GSEM covariance equals R covariance. Those raw discrepancies
are preserved and explained below.

## Supported model

Unweighted glmmTMB Poisson-log maximum likelihood with one group and one
correlated numeric random slope, `(1 + x | id)`. The system includes fixed
effects, both log standard deviations, and Fisher-transformed correlation.
Population predictions integrate the Gaussian random effects; marginal
slopes and uncertainty include covariate-dependent random variation.

Existing fixed-effect factors/interactions, default group clustering, higher
clusters containing whole groups, and unequal sample patterns are supported.
Weights, offsets, zero inflation, constraints, REML, multiple slopes/groups,
diagonal covariance, and other random-slope families remain excluded. The
full contract is in `GLMMTMB-POISSON-RS-DESIGN-20260923.md`.

## Returned Stata evidence

Six Laplace component fits and the four native balanced margins completed
in v2. Component coefficients differ from R by at most 6.76e-6, covariance
entries by 2.52e-7, and relative SEs by 5.56e-6. Native marginal means/slopes
differ by at most 6.25e-6, with relative SE differences at most 5.92e-6.

Cold joint optimization stalled. Supplying the successful Stata component
estimates resolved this: balanced v4 and both remaining v5 fits converged
after one iteration. Each zero-step likelihood reproduced the component sum,
and each final fit returned finite robust and model-based covariance.

| Joint case | Observations / clusters | Overlap | Log likelihood | Max coefficient difference from R |
|---|---:|---:|---:|---:|
| Balanced | 800 / 100 | 800 | -2812.42881825605 | 4.57e-5 |
| Partial, higher clusters | 800 / 20 | 755 | -2728.28808605051 | 3.56e-6 |
| Disjoint groups | 800 / 100 union | 0 | -1383.08810019762 | 2.28e-6 |

The largest likelihood difference from R is 7.56e-7. The zero-iteration
"convergence not achieved" messages are expected diagnostic output; final
acceptance separately requires convergence of the subsequent fitted model.

## Covariance agreement and the raw GSEM limitation

Maximum relative standard-error differences, expressed as percentages:

| Joint case | Raw GSEM, cluster correction aligned | Fully independent sandwich vs R | Stata component curvature + recovered Stata joint score cross-products vs R |
|---|---:|---:|---:|
| Balanced | 0.1796% | 0.00129% | 0.000832% |
| Partial, higher clusters | 0.9952% | 0.000341% | 0.000596% |
| Disjoint groups | 0.7049% | 0.000203% | 0.000417% |

For disjoint models, R uses 50/49 per component and GSEM uses 100/99 over the
union. Aligning GSEM covariance to R therefore multiplies it by **99/98**.
The unaligned disjoint SE difference is 0.8605%. Raw matrices are preserved;
normalization is applied only in labeled comparisons.

The independent calculation solves each group's bivariate Poisson conditional
mode, evaluates its Laplace likelihood including both determinant terms, and
uses central differences for scores and Richardson-extrapolated second
differences for curvature. Missing outcomes are omitted within each model;
wholly absent groups contribute zero. Higher-cluster scores are sums of whole
group scores. Halving numerical steps changes any score by at most 1.78e-9
and any inverse-information entry by at most 1.80e-9 across the six fits.

At the returned joint parameters, recover Stata's score cross-product M from
`V = c B M B'`, using returned robust V, model-based B, and cluster correction
c. This recovers the cross-product matrix, not individual score vectors.
Recombining independent scores with Stata B reproduces raw GSEM SEs within
0.000457%. The remaining differences are therefore chiefly in the returned
joint model-based covariance, not the group-score aggregation.

With cross-equation latent covariances fixed at zero, the likelihood
factorizes and its model-based curvature has zero cross-equation blocks.
Returned GSEM blocks instead reach 3.37e-5, 5.54e-5, and 1.36e-4 for the three
cases. The exact numerical mechanism causing this is not established.
A sensitivity calculation retains the nonstationary score term when changing
variance coordinates; it does not alter the conclusion.

The second reconstruction uses **only Stata-derived quantities**: validated
component OIM covariance as block-diagonal B and recovered joint M. It differs
from R covariance by at most 8.75e-7 overall. No R estimates enter this
reference. Independent curvature and scores differ from R covariance by at
most 6.71e-7. Both checks retain the existing independent-audit thresholds of
1e-6 absolute covariance difference and 1e-4 relative SE difference.

The disjoint cross-equation covariance is exactly zero in R and the fully
independent calculation, and at most 1.58e-14 in the Stata-only reconstruction.
Raw GSEM instead reports up to 2.48e-4 (2.51e-4 after correction alignment).
These raw numerical terms are not copied into the R implementation.

## Regression fixtures and reproduction

Permanent fixtures retain the Stata component estimates, v4 balanced output,
and v5 partial/disjoint output, including full raw 16-by-16 matrices and the
12-by-16 transformation Jacobian. Generic pending flags in parser metadata
mean that parsing alone does not certify raw covariance; this report and the
audit results record the completed validation assessment.

`test-glmmtmb-poisson-rs-stata.R` checks all six native component fits and all
three reconstructed joint sandwiches, including the full cross-equation
blocks, overlap/cluster counts, and exact zero R disjoint cross-covariance.
The reconstruction in the test is explicit; raw joint covariance is never
silently substituted or accepted with a widened tolerance.

From `tools/stata-benchmarks`, regenerate the R reference with
`run_glmmtmb_poisson_rs_reference.R`. Then use:

- `build_glmmtmb_poisson_rs_stata_components.R <v2-log>`
- `build_glmmtmb_poisson_rs_stata_balanced.R <v4-log>`
- `build_glmmtmb_poisson_rs_stata_remaining.R <v5-log>`
- `compare_glmmtmb_poisson_rs_stata_joint.R [output.csv]`
- `check_glmmtmb_poisson_rs_joint_audits.R [output.rds]`

The shared parser reproduces the earlier balanced fixture identically.
Independent review separately parsed both new raw covariance matrices,
checked the Jacobian, and reran all three audits successfully.

## Verification and handoff

On R 4.5.2, glmmTMB 1.1.14, and marginaleffects 1.0.0:

- The focused Stata regression file has 53 passing expectations and no warnings.
- The full suite has 1,334 passing expectations, zero failures/skips, and three
  warnings in test output (the same warning count as the earlier checkpoint).
- `R CMD check --no-manual` finishes with Status: OK: zero errors, warnings,
  or notes, including examples and vignette rebuilding.

The final source rebuild changes only NEWS/README wording and its packaging
timestamp; code, tests, fixtures, documentation, and generated vignettes are
byte-identical to the checked package. No package implementation changed
during the returned-log audit. GitHub's platform matrix has not run for this
unpublished checkpoint.

The checkpoint ZIP contains the complete Git history as a bundle, the built
R source package, these findings, returned Stata logs, numerical comparisons,
and verification logs. All archive entries are at its root. Keep this ZIP as
the current checkpoint; older benchmark rerun ZIPs are historical.

The next model increment can be one correlated numeric random slope for
NB2-log, followed by binomial-logit, using the same restricted structure and
staged validation. Those routes are not implemented by this checkpoint.
Do not push or merge until Trent authorizes publication.
