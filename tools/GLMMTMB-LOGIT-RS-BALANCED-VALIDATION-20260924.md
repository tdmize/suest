# Binomial-logit random slopes: returned balanced validation

Date: 2026-09-24. Development version: 0.1.5.9000.
Returned source: `suest_r_glmmtmb_logit_rs_balanced_v1.log`, Stata 19.5.

**Balanced component results and the independent joint covariance audit pass.**
Raw joint GSEM covariance does not meet the strict agreement criterion: its
largest SE difference from R is 4.23%. The audit localizes this difference
to GSEM's numerical likelihood curvature. Both the raw and independently
computed results are retained. This is not a claim of raw GSEM covariance
parity. No production R covariance calculation was changed in this increment.

## Returned fits and direct comparisons

Both melogit components and the joint GSEM converged after one full-model
iteration, with 1,200 observations and 100 groups. Preliminary fixed-effects
fits converged in three iterations. All checked starts were preserved exactly.
Zero-step component log-likelihood gaps were -1.11e-7 and -1.21e-7; the joint
zero-step gap from the component sum was 7.26e-9. No restart was needed.

Comparisons transform Stata's natural variance/covariance parameters to the
package's log SDs and Fisher correlation, including the full Jacobian.

| Result | Maximum absolute coefficient difference | Absolute log-likelihood difference | Maximum absolute covariance difference | Maximum relative SE difference |
|---|---:|---:|---:|---:|
| y1 component | 3.30e-6 | 8.18e-8 | 1.54e-7 | 0.000195% |
| y2 component | 5.07e-6 | 1.19e-7 | 1.24e-7 | 0.000308% |
| Raw balanced GSEM | 5.20e-6 | 4.90e-7 | 0.00469208 | 4.22684% |

The 16-column raw GSEM result contains four fixed unit loadings and twelve
free parameters. The likelihood separates into two blocks because all four
cross-equation latent covariances are fixed at zero. Robust cross-equation
covariance is still permitted through correlated group scores.

## Independent covariance diagnosis

`logit_rs_joint_audit_helpers.R` independently computes each Bernoulli
conditional mode and bivariate Laplace group likelihood without package or
TMB score code. It evaluates at both the returned Stata joint parameters and
the R parameters. Centered group-score differences use steps 1e-5 and 5e-6;
Richardson-extrapolated likelihood Hessians use .002, .001, and .0005.
Maximum score-step change is 3.56e-10 and bread-step change is 1.43e-9.
No covariance or SE acceptance threshold was relaxed.

For the sandwich `V = c B M B'`, `B` is inverse observed information and `M`
is the group score cross-product matrix. Both systems use `c = 100/99` here.
The independent likelihood curvature is block diagonal. Returned GSEM
model-based covariance has a maximum cross-block entry of 0.00219223 where
the independent likelihood implies zero.

At Stata's parameter point, recover the score cross-products as
`M_stata = solve(B_stata, V_stata) %*% solve(B_stata) / c`.
This recovers a cross-product matrix, not unique score vectors. It agrees
closely with independent scores: maximum scaled discrepancy is 7.30e-6.
Using independent scores with Stata's supplied bread reconstructs raw Stata
covariance within 1.25e-7 absolute and 0.000348% relative SE. Replacing the
curvature with independently computed curvature removes the material gap.

At R's parameter point, the independent sandwich agrees with R covariance
within 5.49e-8 absolute and 0.000074% relative SE when evaluated in the same
native coordinate convention. The accepted bounds remain 1e-6 absolute
covariance and 1e-4 relative SE (0.01%).

Away from an exact stationary point, the Hessian must retain coordinate
chain-rule terms. With exported `z = atanh(rho)` and native glmmTMB
`t = sinh(z)`, the delta-transformed native inverse information is
`inverse(I_z + g_z*rho*E_66)`. The audit also retains the natural-variance
coordinate sensitivity calculation. These score-dependent terms are small
here and cannot explain the raw 4.23% GSEM SE gap. Raw GSEM outputs are never
overwritten with independently audited values.

## Native Stata probabilities and slopes

Independent adaptive normal integration, evaluated at each returned Stata
component parameter vector, reproduces its printed population probability
and continuous slope. Analytic gradients include the full variance function
`A + 2*C*x + B*x^2`, both log SDs, and the Fisher correlation. Variances are
then calculated with Stata's own returned full component covariance.

| Quantity | Stata estimate | Absolute independent estimate difference | Absolute relative SE difference |
|---|---:|---:|---:|
| y1 probability | 0.450037334369 | 1.28e-10 | 3.40e-11 |
| y1 slope of x | 0.095710934467 | 4.02e-11 | 3.27e-8 |
| y2 probability | 0.463437339379 | 5.07e-12 | 1.87e-9 |
| y2 slope of x | 0.076952223075 | 1.81e-11 | 8.87e-8 |

These checks separate numerical prediction/gradient agreement from small
optimizer endpoint differences. Comparing margins at R's slightly different
endpoints yields differences as large as 9.31e-7; that is not an integration
error. Stata used 40 quadrature points for prediction; its fitted likelihood
remained Laplace.

## Regression fixture and next gate

The returned-log parser requires revision/completion markers, all convergence
and sample checks, the three zero-step checks, finite full matrices, four
fixed unit loadings with zero covariance rows, and positive-definite free
covariance matrices. It stores raw matrices, full Jacobians, transformed
matrices, native margins, source filename, and source MD5. The audit is a
separate fixture. The added regression test checks points, native component
covariances/margins, independent joint covariance, and the raw-Stata
sandwich reconstruction at matching evaluation points.

The next gate is `suest_r_glmmtmb_logit_rs_partial_v1.do`. It runs only two
components and one joint fit: 1,175 and 1,180 observations, 1,155 overlapping,
1,200 in the union, 100 groups, and 20 higher clusters for the joint VCE.
R starts come from the already checked six-fit reference. Joint starts use
the actual returned Stata components. Full and preliminary fixed-effects
fits are capped at 15 iterations; zero-step checks precede full estimation,
and the first failed check stops the script. No automatic restarts.
Partial and disjoint Stata validation remain pending.

## Final verification

- New balanced regression: 28 assertions passed, with no failures, warnings,
  errors, or skips.
- Full source-package suite: 1,564 passed, 0 failed, 0 skipped, and the same
  three existing test warnings as the preceding checkpoint.
- R CMD build succeeded; R CMD check --no-manual returned **Status: OK**,
  including tests, examples, and vignette rebuilding.
- Read-only independent review verified fixture provenance and reparsing,
  covariance-coordinate corrections, native-margin gradients, and the bounded
  partial gate. Its analytic/numerical margin-gradient difference was 6.97e-12.
- Production R implementation is unchanged in this validation increment.

## Reproduction

Run in the repository's `tools/stata-benchmarks` directory:

```
Rscript build_glmmtmb_logit_rs_stata_balanced.R /path/to/suest_r_glmmtmb_logit_rs_balanced_v1.log
Rscript compare_glmmtmb_logit_rs_stata_balanced.R comparison.csv
Rscript check_glmmtmb_logit_rs_balanced_audit.R audit.rds
Rscript check_glmmtmb_logit_rs_balanced_margins.R margins.csv
Rscript make_glmmtmb_logit_rs_partial_starts.R
```

The independent margin helper is included under `tests/testthat`. The saved
R reference and preflight are included in the checkpoint alongside the raw
log, raw/audit fixtures, reports, and full verification output. The complete
Git bundle preserves the repository scripts excluded from R source builds.
