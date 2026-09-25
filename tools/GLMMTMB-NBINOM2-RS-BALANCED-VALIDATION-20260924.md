# NB2 random slopes: balanced Stata validation and next gate

The returned `suest_r_glmmtmb_nbinom2_rs_balanced_v1.log` completed both
Laplace menbreg components and the balanced joint Laplace GSEM after one
full-model iteration each. Starting-value likelihood and parameter checks
passed. There are 1,000 observations and 100 groups. The implementation
remains development version 0.1.5.9000. No production R code changed during
this validation increment; no remote push was performed.

**The balanced point estimates agree, and independent audits support R's
joint covariance. Raw Stata covariance is numerically different and is not
treated as an exact reference. Partial-overlap and disjoint Stata validation
remain pending.**

## Returned comparisons

Full-precision output was parsed with strict matrix-dimension and convergence
checks. The component parameter order is `(intercept, x, z, log_phi,
log_sd_intercept, log_sd_slope, atanh_rho)`. Stata's `lnalpha` changes sign;
natural variance/covariance parameters transform with their complete Jacobian.
Four fixed latent loadings are excluded from the 18-column joint result.

| Comparison with R | y1 component | y2 component | Raw balanced joint |
| --- | ---: | ---: | ---: |
| Largest absolute parameter difference | 2.73e-5 | 2.13e-5 | 2.77e-5 |
| Absolute log-likelihood difference | 9.66e-7 | 6.25e-7 | 1.66e-6 |
| Largest absolute covariance difference | 6.59e-5 | 2.45e-5 | 2.68e-4 |
| Largest relative coefficient-SE difference | 0.2532% | 0.1103% | 1.1643% |

Native marginal means differ by at most 2.31e-8 and x slopes by 1.85e-6.
Native margin SEs differ by up to 0.2645%. Recalculating Stata's margin
variances analytically from its own transformed covariance agrees within
1e-7 absolute variance. Thus the margin-SE discrepancy follows the supplied
covariance, rather than a different definition of the integrated mean.

## Independent likelihood and covariance audit

The audit uses a separate bivariate Gaussian Laplace calculation with no
glmmTMB/TMB score or covariance code. Latent modes are found by Newton steps
with line search and verified gradients. The conditional NB2 gradient is
`phi*(mu-y)/(phi+mu)`, with curvature
`phi*mu*(phi+y)/(phi+mu)^2`. The Gaussian prior and its normalization are
included, together with the mode-Hessian determinant.

All seven parameter scores are finite-differenced for every group. Curvature
uses two Richardson extrapolations from steps .002, .001, and .0005. Score
step-halving changes are at most 1.43e-9; inverse-curvature changes are at
most 1.25e-9. Independent log likelihood agrees with the returned joint fit
within 1.70e-6.

Two evaluation points are retained separately:

- **R's fitted parameters:** the independent sandwich and suest covariance
  differ by at most 5.00e-7; relative SE differences are at most 7.85e-6
  (0.000785%). Independent component curvature agrees with R native
  covariance within 2.78e-7 and relative SE within 3.95e-6.
- **Stata's returned joint parameters:** using Stata's model-based covariance
  with independent score cross-products reproduces its raw robust covariance
  within 1.90e-7 and relative SE within 3.62e-6 (0.000362%).

This distinction matters. Evaluating the independent sandwich at Stata's
slightly different estimates and comparing it with R's original fitted
system gives an absolute covariance difference of 1.32e-6. That exceeded
the initial 1e-6 gate. Repeating the calculation at matching R estimates
passes the unchanged threshold; no comparison tolerance was widened.
Coordinate sensitivity away from an exact optimum was also recorded and
does not account for the 1.16% raw difference.

For the decomposition, `V = c B M B'` with `c = 100/99`. The reported
model-based covariance B and robust covariance V recover the score
**cross-product matrix** M. They do not uniquely identify individual score
vectors. Recovered Stata and independent M agree within 8.18e-6 on a
diagonal-standardized scale. Stata's joint model-based covariance has
cross-equation entries as large as 8.57e-5 even though this likelihood
factorizes and its exact cross-information is zero. These comparisons
locate the raw discrepancy chiefly in numerical likelihood curvature.

Unlike the preceding Poisson increment, substituting Stata component OIM
curvature does not produce a tight equality reference: reconstructed SEs
still differ by up to 0.7091%. This is retained as diagnostic evidence,
not accepted as exact Stata parity. No package covariance was altered to
force agreement.

## Verification and reproducibility

`build_glmmtmb_nbinom2_rs_stata_balanced.R` parses the returned log into the
raw fixture, retaining its source MD5. `compare_glmmtmb_nbinom2_rs_stata_balanced.R`
produces the comparison table. `check_glmmtmb_nbinom2_rs_balanced_audit.R`
reproduces both independent audits and checks the matching-point limits of
1e-6 absolute covariance and 1e-4 relative SE. Raw and independently audited
fixtures are stored separately. The new Stata regression test covers
point estimates, native curvature, joint covariance, margins, and the
Stata-curvature/independent-score decomposition.

The full 1,397-expectation suite passed at parent checkpoint `9e9ca45`.
This increment adds a focused 23-expectation regression test; runtime R
implementation files are unchanged. Fresh focused output and source-package
build/check records accompany this checkpoint; a package check with
`--no-tests` is explicitly separate from the already-run focused tests.

## Next bounded Stata run

Extract the supplied flat ZIP and run in its directory:

```stata
do suest_r_glmmtmb_nbinom2_rs_partial_v2.do
```

Return `suest_r_glmmtmb_nbinom2_rs_partial_v2.log`, including a failed log.
The CSV is unchanged. This run has only two partially observed components
and one joint GSEM: component sample sizes 975 and 980, overlap 955, union
1,000, 100 groups, and 20 higher clusters. Component seeds come from the
previously checked R fits. The joint starts from actual Stata component
estimates. All starts are checked by unchanged parameters and matching
zero-step likelihoods; every full mixed model has a 15-iteration cap.

The returned v1 log showed preliminary menbreg fixed-only fits despite
`startvalues(zero)`, contrary to the earlier script-review assumption.
They completed in four or five iterations. Revision 2 requests the documented
`startvalues(fixedonly, iterate(15))` control for component starting fits;
`from()` still supplies the full checked mixed-model vector. GSEM retains
`startvalues(zero)` and `adaptopts(tolerance(1e-12))`. The iteration limits
are not wall-clock timeouts. The script stops at the first failed check,
with no automatic restarts or quadrature comparisons. Stata is unavailable
in this runtime; revision 2 is statically reviewed, not executed here.

After the partial log, apply the same point/covariance audit before preparing
the odd/even disjoint-group gate. Do not describe all NB2 random-slope sample
patterns as Stata-validated yet.

Starting-value controls:
[Stata meglm manual, Maximization](https://www.stata.com/manuals/memeglm.pdf).
