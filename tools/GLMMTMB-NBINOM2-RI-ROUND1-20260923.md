# NB2 returned Stata gate: boundary diagnosis and revision 2

Historical round-1 record. The revision-2 gate is now complete; see
`GLMMTMB-NBINOM2-RI-VALIDATION-20260923.md` for the final comparison.

Source: `suest_r_glmmtmb_nbinom2_ri_crosslang_benchmark.log`, Stata 19.5,
23 September 2026. Revision-1 development commit: `4b29641`.

## Finding

Five Laplace menbreg component fits completed. The sixth, `y2_right`, stopped
with r(430): numerical derivatives encountered a discontinuous region with
missing values. The log stops before joint-system comparisons.

The corresponding R fit is not a regular interior solution. Its estimated
random-intercept variance is `2.122698e-8`; the native standard error of
log(SD) is about `3098.8`. Nevertheless, glmmTMB returned convergence code zero
and a positive-definite Hessian flag. Checking only those flags and floating
point underflow was inadequate.

Refitting from four SD starting values produces the same limiting solution.
A profile with all other parameters re-estimated improves toward zero variance:

| Fixed variance | Laplace log likelihood |
|---:|---:|
| 1e-12 | -330.105983379 |
| 1e-8 | -330.105983391 |
| 1e-6 | -330.105984607 |
| .0001 | -330.106107462 |
| .001 | -330.107352612 |
| .01 | -330.131936180 |
| .05 | -330.453680895 |
| .10 | -331.163091323 |
| .20 | -332.986013870 |

Ordinary NB2 without a random intercept has log likelihood
`-330.105983376919`. This supports a zero-variance boundary diagnosis; it does
not justify increasing iterations or treating R's convergence flag as success.
The exact internal Stata optimizer trajectory cannot be reconstructed from
the quiet log, but the independent boundary evidence is decisive for excluding
this case from the regular inference benchmark.

## R correction

The NB2 adapter now evaluates the exact zero-random-effect likelihood at the
fitted fixed effects and dispersion. The fitted Laplace likelihood must improve
on that boundary point by more than
`sqrt(.Machine$double.eps) * max(1, abs(fitted_logLik))`. Likelihoods must be finite.

This is a numerical support restriction, not a significance test, a profiled
likelihood-ratio test, or a universal singularity diagnostic. It can miss
problematic fits because the other parameters are held fixed at the interior
fit's estimates. The implementation does not refit or mutate the model, and
does not claim to detect every weakly identified variance or the Poisson
dispersion limit. Existing logit/Poisson boundary checks are unchanged.

The original failing data is retained in a negative regression fixture.
The test failed before the guard was added and passed afterwards.

## Completed Stata comparisons

After converting lnalpha to log_phi and random-effect variance to log_sigma,
including both covariance Jacobians, the five completed component fits match R:

| Quantity | Largest absolute difference |
|---|---:|
| Coefficient | 3.721e-5 |
| Native covariance element | 6.741e-5 |
| Log likelihood | 1.791e-7 |
| Balanced marginal mean | 1.578e-5 |
| Balanced average slope | 9.129e-6 |
| Balanced margin standard error | 6.018e-5 |

The largest relative native standard-error difference is below 0.067%.
Five component fixtures are now tested. Joint NB2 parity remains unverified.

## Revision-2 positive benchmark

The base generated data is unchanged. Only the disjoint masks change from
first/last 40 groups to odd/even groups. This remains 240 observations and
40 nonoverlapping groups per model. This is a deliberately regular positive
test; the discarded allocation is retained as a boundary rejection test.

| Model | RI variance | Native SE of log(SD) |
|---|---:|---:|
| Y1, odd groups | .390917 | .201375 |
| Y2, even groups | .273061 | .215984 |

The new reference script checks every model with the adapter and requires
variance above .01 and SE of log(SD) below 1 for this positive fixture.
Those last two requirements are benchmark checks, not package acceptance rules.

Run `suest_r_glmmtmb_nbinom2_ri_crosslang_benchmark_v2.do` and return the v2 log.
No GitHub publication is included in this revision.

Method context: https://glmmtmb.github.io/glmmTMB/articles/troubleshooting.html
