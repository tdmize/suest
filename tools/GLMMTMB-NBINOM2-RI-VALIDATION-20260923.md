# NB2 random-intercept validation — 23 September 2026

The returned revision-2 Stata 19.5 benchmark completes the NB2 increment's
cross-language gate. All six Laplace `menbreg` components, all three joint
Laplace `gsem` systems, and all three optional 12-point adaptive `suest2`
systems completed. No further Stata execution is required for this increment.
The code remains unpublished development version 0.1.5.9000.

Source: `suest_r_glmmtmb_nbinom2_ri_crosslang_benchmark_v2.log`. The fixture
records its MD5 and retains both raw Stata matrices and transformed matrices.
R comparison: R 4.5.2, glmmTMB 1.1.14, marginaleffects 1.0.0.

## Scope and transformations

The supported model is unweighted `glmmTMB::nbinom2("log")` with one Gaussian
random intercept and one estimated constant dispersion. The parameter vector
retains fixed effects, `log_phi`, and `log_sigma`. The full nuisance covariance
enters the joint system and marginal-effect uncertainty calculation. Response
predictions integrate over the Gaussian random effect.

Stata's `lnalpha` becomes `log_phi = -lnalpha`; its random-intercept variance
becomes `log_sigma = log(variance)/2`. The covariance Jacobian uses -1 and
`1/(2*variance)`, respectively. GSEM's two fixed unit factor loadings and their
zero covariance rows/columns are checked and removed. No numerical estimates
from R are used to construct the Stata fixture.

## Component and marginal results

Maximum absolute differences across all six component models:

| Quantity | Maximum difference |
| --- | ---: |
| Coefficients, including both nuisance parameters | 3.721e-5 |
| Native covariance elements | 8.690e-5 |
| Log likelihood | 1.791e-7 |
| Relative native standard error | 0.188% |
| Balanced population means | 1.578e-5 |
| Balanced average slopes for x | 9.129e-6 |
| Native marginal mean/slope standard errors | 6.018e-5 |

The two balanced models have Stata `margins` output for means and slopes;
the other four component models do not include that output. Joint-system
prediction and contrast standard errors also have independent analytic
delta-method tests in R.

## Joint Laplace comparison

All systems have 480 union observations. Overlap counts are 480, 444, and 0;
cluster counts are 80, 16, and 80 for balanced, partial/higher, and disjoint.
Joint log likelihoods differ by at most 3.036e-7.

| System | Max coefficient difference | Max raw covariance difference | Max relative raw SE difference |
| --- | ---: | ---: | ---: |
| Balanced | 1.614e-5 | 5.748e-5 | 0.110% |
| Partial overlap, higher clusters | 3.883e-5 | 3.965e-5 | 0.073% |
| Disjoint groups | 2.347e-5 | 0.001272 | 1.017% |

The disjoint models each contain 40 groups. R's specialized panel SUEST route
uses each model's `40/39` correction; joint GSEM counts the union's 80 groups
and uses `80/79`. Multiplying the returned GSEM covariance by
`(40/39)/(80/79)` aligns these conventions. After that adjustment, the largest
disjoint covariance difference is 1.296e-4 and relative SE difference is
0.385%. The raw Stata covariance is preserved unchanged in the fixture.

GSEM also returns disjoint cross-covariance up to 6.426e-5, while R and
adaptive `suest2` return exactly zero. With disjoint groups and a separable
likelihood, the theoretical cross-covariance is zero. The small GSEM residual
is consistent with numerical derivative/information-matrix error; the log
does not isolate its internal source. Exact covariance equality is therefore
not claimed. R's block-diagonal information and group influences are checked
independently, including the `40/39` factor in both diagonal blocks.

The regression limits are 5e-5 for coefficients and 1e-6 for likelihoods.
After cluster-correction alignment, covariance limits are 1e-4 for balanced
and partial systems and 2e-4 for disjoint; relative SE limits are 0.2% and
0.5%, respectively. These bound the observed cross-engine numerical residuals
while the independent R sandwich checks use 1e-9 and exact-zero checks use
1e-12. They do not relax the package's statistical formula.

Stata's maximum-likelihood robust multiplier and cluster convention are
documented in [P] `_robust`, Options and Methods and formulas:
https://www.stata.com/manuals/p_robust.pdf. The returned `e(N_clust)` values
establish the cluster counts used in this benchmark.

## Adaptive comparison and boundary protection

`suest2` rejects Laplace fits with return code 322, requiring adaptive
quadrature. The optional 12-point adaptive comparisons all returned code 0.
Against R's Laplace fits, their maximum coefficient difference is 0.01432,
covariance difference is 0.002188, and relative SE difference is 1.55%.
These are distinct likelihood approximations, so adaptive coefficients and
covariance are retained as a secondary comparison, not an equality gate.

Revision 2 changes only the disjoint allocation to odd/even IDs. Original
outcomes and predictors remain unchanged. The original revision-1 right-half
Y2 fit is preserved as a negative regression fixture: it is numerically at
zero random-intercept variance despite nominal R convergence. The adapter's
zero-boundary likelihood guard remains in place. It is a numerical support
restriction, not a significance test or a comprehensive weak-identification
diagnostic. Random slopes, offsets, weights, zero inflation, dispersion
regressors, constraints, and more complex random-effects structures remain
outside this increment.

## Reproduction and status

From `tools/stata-benchmarks`, rebuild the fixture using
`Rscript build_glmmtmb_nbinom2_ri_stata_fixture.R <returned-v2-log>`.
Run `run_glmmtmb_nbinom2_ri_reference_v2.R` against the installed development
package, then `compare_glmmtmb_nbinom2_ri_stata.R [output.csv]` to reproduce
the numerical comparison. The archived Stata log and source package accompany
the checkpoint.

The built package passes `R CMD check --no-manual` on Linux/R 4.5.2 with
0 errors, 0 warnings, and 0 notes. Its test suite reports 1,233 passing
expectations, 0 failures, 0 skips, and three existing GLM warnings in test
output. The focused NB2 suite has 124 passing expectations with no warnings.
The fixture rebuilds byte for byte from the returned log. Independent review
found no blocking issues in the parsing, transformations, regression limits,
or documentation. No production-code change was needed after the v2 return.

The checkpoint handoff records the local Git commit. The published GitHub
release remains 0.1.5; this increment has not been pushed or tested on GitHub's
platform matrix. Random slopes are a separate next increment.
