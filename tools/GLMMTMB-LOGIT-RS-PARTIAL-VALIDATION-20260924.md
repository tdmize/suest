# Logit random slopes: partial-overlap Stata validation

Date: 2026-09-24. Development version: 0.1.5.9000.

The returned partial-overlap gate passed. Both components and the joint
GSEM converged after one full-model iteration. Independent joint SEs agree
with R within 0.000078%. Raw GSEM joint SEs differ by up to 1.14%; those raw
results remain separate from the independent reference. No production R code
or acceptance tolerances changed in this increment. Disjoint validation is next.

## Returned fits and direct comparison

The log is `suest_r_glmmtmb_logit_rs_partial_v1.log`, revision 1, with
`LOGIT_RS_PARTIAL_HIGHER_COMPLETE=1`. Component samples contain 1,175 and
1,180 observations in 100 panels each; overlap is 1,155 and union is 1,200.
The joint VCE uses 20 higher clusters. All requested starting parameters are
preserved exactly in the zero-step checks; likelihood gaps are below 9.1e-8.

| Fit | Maximum absolute coefficient difference | Maximum absolute covariance difference | Maximum relative SE difference | Absolute log-likelihood difference |
|---|---:|---:|---:|---:|
| y1_partial | 2.367e-6 | 1.411e-7 | 3.627e-6 | 4.613e-8 |
| y2_partial | 6.355e-6 | 1.134e-7 | 1.988e-6 | 8.784e-8 |
| Joint, raw GSEM | 5.889e-6 | 9.548e-4 | 1.1445e-2 | 7.814e-8 |

Differences are in the package's fixed-effect, log-SD, and Fisher-correlation
coordinates. The joint raw covariance comparison does not meet the strict
parity threshold. The raw fixture preserves the original Stata matrices,
parameter stripes, transformations, sample counts, and source-log MD5.

## Independent covariance audit

The independent bivariate Laplace calculation polishes panel modes and
computes scores and curvature without using the package's score calculation.
It aggregates panel scores to the same 20 higher clusters. Both systems use
the 20/19 correction; no additional normalization is needed. Native-coordinate
curvature retains the score-dependent Hessian chain-rule term when changing
correlation coordinates away from an exactly stationary point.

At the R point, the native-coordinate independent sandwich differs from R
by at most 8.835e-8 in covariance and 7.730e-7 in relative SE (0.0000773%).
The pure Fisher-coordinate calculation has maximum covariance difference
8.555e-8 and relative SE difference 5.218e-7. Native component covariance
audits agree within 3.313e-8 absolute covariance and 3.847e-7 relative SE.
The existing acceptance limits remain 1e-6 absolute covariance and 1e-4
relative SE.

Using Stata's own model-based curvature with independently calculated scores
reproduces its raw sandwich within 2.836e-7 absolute covariance and 5.513e-6
relative SE. Recovered Stata score cross-products differ from the independent
ones by at most 7.751e-4 absolute and 1.104e-5 scaled. Stata's model-based
covariance has a maximum cross-equation block entry of 0.0004708864, although
the independent likelihood curvature is block diagonal under the imposed
cross-latent restrictions. This localizes the raw GSEM discrepancy to numerical
curvature. Independent likelihoods agree with Stata within 8.175e-8; score
step sensitivity is below 4.441e-10 and curvature step sensitivity below
2.971e-9. Raw and audited fixtures remain separate; production covariance
continues to use the supplied native covariance.

## Native probabilities and slopes

Independent adaptive integration and analytic gradients evaluate each Stata
component at its own endpoint and on its own estimation sample, retaining
all nuisance covariance. Mean probabilities differ by less than 1.435e-10;
mean slopes differ by less than 6.669e-11. Relative SE differences are below
4.166e-8. Comparisons at R's slightly different endpoint are recorded separately.

## Next bounded gate

Run `do suest_r_glmmtmb_logit_rs_disjoint_v1.do` and return
`suest_r_glmmtmb_logit_rs_disjoint_v1.log`, including an early-stopping log.
The gate uses 600 observations and 50 panels per component, with zero overlap
and 1,200 observations in 100 union panels. It retains checked starts,
zero-step checks, two components and one joint fit, 15-iteration caps on
preliminary and full fits, and no retries. The R preflight passes for both
components. On return, distinguish R's component correction 50/49 from Stata's
union correction 100/99 before comparing the disjoint sandwich matrices.

## Reproduction and verification

The partial parser, comparison, independent audit, and native-margin scripts
are in `tools/stata-benchmarks`. Regression tests cover both balanced and
partial fixtures, including higher-cluster aggregation, sample-specific
margins, strict covariance bounds, and separate raw-covariance reconstruction.
The raw fixture records source provenance; the audit fixture records the
independent computations. The checkpoint retains the original returned log.

Focused regression: 56 assertions passed, with no failures, warnings, errors,
or skips. Read-only review found no actionable issues in the partial scripts,
fixtures, tests, or disjoint gate. Full package verification is recorded below.

Full source-package verification: **1,592 passed, 0 failed, 0 skipped**, with
three existing test warnings. `R CMD build` succeeded and
`R CMD check --no-manual` returned **Status: OK**, including examples and
vignette rebuilding. Environment: R 4.5.2, glmmTMB 1.1.14, TMB 1.9.25,
marginaleffects 1.0.0. The checkpoint includes full verification logs and a
complete-history Git bundle. Nothing was pushed or released.
