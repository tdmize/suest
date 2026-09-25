# Poisson random-slope Stata convergence diagnostics

## Returned evidence

The revision-1 log from Stata 19.5 stopped at `y1`, the first of six
`mepoisson` components. It printed "convergence not achieved" twice and
returned r(430). It contains no parameter estimates, likelihood, or iteration
history. The benchmark used `quietly`, which concealed the details needed to
distinguish poor starting values, an optimization problem, and other numerical
failures. No joint GSEM or marginal-effects comparison ran.

This does not establish that the model is too complicated or that it has a
variance boundary. Stata's failure cause is still unconfirmed.

## Fresh R evidence on the unchanged data

`check_glmmtmb_poisson_rs_initialization.R` fits each of the six components
twice: glmmTMB defaults and BFGS starting from fixed-only Poisson coefficients,
unit random variances, and zero correlation. All 12 fits report convergence
and positive-definite Hessians. Maximum absolute log-likelihood difference is
5.91e-9; maximum native-parameter difference is 5.65e-6. The maximum absolute
native gradient across these fits is 0.001971.

The default fits' random-covariance minimum eigenvalues range from 0.0772 to
0.1477, with correlations from 0.1115 to 0.5239. Thus these R fits show no
evidence of the near-zero random variance encountered in the earlier NB2
benchmark. This checks the R solutions; it does not prove why Stata failed or
validate Stata/R agreement. Runtime: R 4.5.2, glmmTMB 1.1.14, Linux.

## Revision-2 experiment

The CSV, formulas, samples, free within-model random covariances, fixed-zero
cross-model latent covariances, clustering, and convergence tolerances are
unchanged. No package implementation is modified.

For each component and joint system, the script:

1. Attempts Laplace with normal starting values and visible iteration and
   tolerance output, capped at 200 iterations.
2. If that fails, fits the same model by mean-variance adaptive quadrature to
   obtain starting values, capped at 100 iterations. Components use 7 points;
   joint systems use 3 to limit the cost of initialization. These fits are
   explicitly labeled as initializers only.
3. If the initializer converges, retries Laplace from its labeled `e(b)`,
   capped at 200 iterations, with the original convergence tolerances.
4. Accepts only return code zero AND `e(converged)==1` on a Laplace attempt.
   Clears estimation results between attempts, logs diagnostics for failures,
   and moves to the next case even if neither attempt succeeds.

No R estimates are used as starting values. The changed initialization is a
testable hypothesis, not a verified remedy. Increasing integration accuracy
is not the target; every accepted comparison still uses Laplace.

The script separately counts six components, three joint systems, and four
native marginal-effect calculations. Margins must return code zero, a finite
scalar estimate, and a finite nonnegative scalar variance. Only all 13 successful calculations
produce `GLMMTMB_POISSON_RS_CROSSLANG_BENCHMARK_COMPLETE=1`. The diagnostic
completion marker merely means the script reached its end. Individual
successes remain subject to numerical comparison with R.

## Run and return

Extract the flat ZIP into the same working folder. It contains the unchanged
CSV and a new do-file; no R installation or Git operation is needed.

```stata
do suest_r_glmmtmb_poisson_rs_crosslang_benchmark_v2.do
```

Return `suest_r_glmmtmb_poisson_rs_crosslang_benchmark_v2.log`, including if
some attempts fail. Let the script reach its final summary when possible.

The do-file has been reviewed against official Stata option documentation but
cannot be executed in this environment, which has no Stata installation.
Fresh R initialization checks pass. The earlier package tests and R CMD check
at b3f7414 still describe the unchanged implementation; they were not rerun
for this benchmark-only revision. Cross-language validation remains pending.

## Documentation checked

- https://www.stata.com/manuals/memepoisson.pdf: Laplace, `from()`, and
  maximization options.
- https://www.stata.com/manuals/memeglm.pdf: default fixed-only starting values
  and labeled starting-value vectors.
- https://www.stata.com/manuals/semgsemestimationoptions.pdf: integration,
  starting values, and maximization options for GSEM.
- https://www.stata.com/manuals/semintro12.pdf: convergence diagnostics and
  reusing estimates as starting values after changing integration methods.
