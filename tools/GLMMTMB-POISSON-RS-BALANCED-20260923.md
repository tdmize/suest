# Poisson random slopes: component validation and focused joint rerun

## What the returned revision-2 log establishes

All six default-start Laplace component attempts stalled and reached the
200-iteration limit with `e(converged)==0`, despite return code zero. The
script correctly rejected them. All six seven-point adaptive initializers
converged, and all six subsequent Laplace fits converged by iteration 2.
The four native margins calculations also completed successfully.

The balanced joint default-start fit stalled near -2881.223 and was
interrupted at iteration 116. The sum of the converged balanced component
log likelihoods is **-2812.42881830556**, about 68.79 units higher. No joint
fit has passed yet. Partial and disjoint joint models were not attempted.

This supplies direct evidence that different initialization resolves the
component failures. Reusing component fits is the next controlled experiment
for the joint model; its success has not yet been demonstrated in Stata.

## Completed comparison with R

The fixture retains original Stata estimates and full covariance matrices.
It transforms `(var_intercept, var_slope, covariance)` to
`(log_sd_intercept, log_sd_slope, atanh_rho)` with the full Jacobian, including
all three derivatives of the correlation transformation.

| Component | Max absolute coefficient difference | Max absolute covariance difference | Max relative SE difference |
| --- | ---: | ---: | ---: |
| y1 | 3.71e-6 | 1.39e-7 | 2.84e-6 |
| y2 | 6.76e-6 | 2.52e-7 | 5.56e-6 |
| y1_partial | 2.87e-6 | 1.11e-7 | 1.84e-6 |
| y2_partial | 3.56e-6 | 1.16e-7 | 2.87e-6 |
| y1_left | 5.19e-8 | 1.19e-8 | 1.86e-7 |
| y2_right | 2.28e-6 | 1.95e-7 | 1.48e-6 |

Maximum absolute log-likelihood difference: 1.06e-6. Across the two balanced
native marginal means and slopes, maximum absolute estimate difference is
6.25e-6 and maximum relative SE difference is 5.92e-6. These SE differences
are fractions, not percentages.

The new component regression file checks all six sets of parameters, full
covariance matrices, standard errors, and likelihoods against Stata. It has
30 expectations. The fixture explicitly labels joint validation as pending.

## Focused revision-3 run

Extract the flat ZIP into the benchmark working folder and run:

```stata
do suest_r_glmmtmb_poisson_rs_balanced_v3.do
```

Return **suest_r_glmmtmb_poisson_rs_balanced_v3.log**.

The do-file embeds the exact successful Stata `y1` and `y2` estimates from
the returned log. It does not refit the components. It obtains GSEM's native
parameter labels with `noestimate`, maps all 12 free parameters by name,
and retains the existing fixed loadings and covariance restrictions.

Before optimizing, a zero-iteration evaluation must reproduce the component
likelihood sum within 1e-4 and preserve the starting parameter vector. That
evaluation is diagnostic; it is never accepted as convergence. If either
check fails, the script stops and leaves the evidence in the log.

The only optimization is the balanced joint Laplace fit, initialized from
those component estimates and capped at **15 iterations**. Original
convergence tolerances and robust group covariance remain unchanged. There
are no automatic fallback fits or other sample patterns in this run.
The completion marker requires return code zero, `e(converged)==1`, 800
observations, 100 clusters, agreement with the component likelihood sum, and
a finite covariance matrix of the expected dimensions.
The returned full covariance matrix still needs comparison with R before
joint validation is considered complete.

## Verification and scope

The component fixture builder and comparison script run successfully under
R 4.5.2/glmmTMB 1.1.14. The new component regression tests pass. The Stata
do-file is reviewed against official documentation; Stata is unavailable in
this environment, so runtime verification requires the returned log.
No R implementation code changed, and the full package suite was not rerun
for this fixture/benchmark increment. No remote push was made.

Official Stata references:

- https://www.stata.com/manuals/semsemandgsemoptionfrom.pdf
- https://www.stata.com/manuals/semgsemestimationoptions.pdf

Next: resolve the balanced joint fit and compare its robust covariance.
Only then resume the partial and disjoint joint checks using component
starting values rather than repeating the failed default-start route.
