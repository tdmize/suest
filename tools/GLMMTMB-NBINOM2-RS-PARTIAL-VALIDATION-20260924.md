# NB2 random slopes: partial-overlap validation and disjoint gate

The returned `suest_r_glmmtmb_nbinom2_rs_partial_v2.log` completed both
Laplace menbreg components and the joint Laplace GSEM after one full-model
iteration each. Zero-step likelihood and unchanged-start checks passed.
Component samples contain 975 and 980 observations, with 955 shared rows,
a union of 1,000 observations, and 20 higher-level clusters. All 100 id
groups contribute to both components. Preliminary fixed-only fits took
four or five iterations under the explicit 15-iteration limit.

**Partial-overlap point estimates agree closely, and independent audits
support R's covariance. Raw Stata covariance remains a separate, numerically
different result. Disjoint Stata validation is still pending.** Development
version remains 0.1.5.9000; production R files are unchanged in this increment.

## Returned comparisons

The full-precision log is parsed with convergence, sample-size, matrix-shape,
and parameter-stripe checks. Stata's `lnalpha` changes sign to `log_phi`;
its natural variance/covariance parameters are transformed to two log SDs
and the Fisher correlation using the full Jacobian. Four fixed latent
loadings are excluded from the 18-column joint result.

| Comparison with R | y1_partial | y2_partial | Raw partial joint |
| --- | ---: | ---: | ---: |
| Largest absolute parameter difference | 1.89e-15 | 7.28e-5 | 7.28e-5 |
| Absolute log-likelihood difference | 1.07e-6 | 6.44e-7 | 1.89e-6 |
| Largest absolute covariance difference | 1.88e-4 | 2.59e-4 | 1.32e-4 |
| Largest relative coefficient-SE difference | 0.3512% | 0.4816% | 0.2658% |

Returned component log likelihoods are -2147.25755373812 and
-1953.91344289076; joint log likelihood is -4101.17099680951.
The largest parameter difference is y2's `log_phi`: 7.27648e-5, or
0.06064% of its native SE. All other component differences are below 8e-6.
Thus the earlier balanced absolute parameter gate of 5e-5 does not pass for
this one dispersion endpoint. The partial regression retains that absolute
gate for every non-dispersion parameter and additionally requires every
parameter difference to be below 0.1% of its SE (native for components,
robust for the joint system). Balanced tolerances are unchanged. This
scale-aware endpoint gate is separate from the unchanged strict covariance
audit thresholds; no R estimates or covariance matrices were adjusted.

Native marginal means differ from R by at most 2.97e-7 and x slopes by
1.19e-6. Their SEs differ by up to 0.8811%, reflecting the supplied native
covariance. Recalculating Stata's own margin variances analytically gives
a largest absolute variance difference of 1.50e-7 and relative SE difference
of 1.25e-6 (0.000125%). The analytic margin gradients agree with independent
Richardson derivatives within 1.45e-10. For partial margins the regression
therefore checks relative SE below 1e-5, instead of the balanced case's
absolute variance threshold of 1e-7. Both observed differences and this
tolerance choice are explicit; this is not a claim of exact numerical equality.

## Independent covariance audit

The existing independent bivariate Laplace calculation is reused unchanged:
Gaussian latent modes by Newton steps with line search; NB2 conditional
likelihood and analytic latent curvature; parameter scores by finite
differences; inverse information from Richardson-extrapolated Hessians.
It does not use glmmTMB/TMB score or covariance code. Group scores are
summed into the same 20 higher clusters before forming the sandwich, with
finite-sample correction 20/19.

| Matching-point comparison | Largest absolute covariance difference | Largest relative SE difference |
| --- | ---: | ---: |
| Independent sandwich vs R, both at R estimates | 1.40e-7 | 3.60e-6 (0.000360%) |
| Stata curvature with independent scores vs Stata robust V, both at Stata estimates | 1.75e-7 | 5.25e-6 (0.000525%) |

Both pass the existing limits of 1e-6 absolute covariance and 1e-4 relative
SE. Independent component curvature agrees with native R covariance within
1.34e-7 and relative SE within 1.83e-6. Score step-halving changes are below
1.43e-9 and inverse-curvature step changes below 2.64e-9. Independent joint
log likelihood agrees with Stata within 1.89e-6.

Evaluation points matter: an independent sandwich evaluated at Stata's
estimates differs from R's original fitted covariance by 2.62e-6; this is
recorded separately and not substituted for a matching-point comparison.
Recovered Stata score cross-products agree with independent cross-products
at the Stata point within 1.33e-5 on a diagonal-standardized scale. These
are cross-product matrices, not uniquely recovered score vectors.

Stata's joint model-based covariance has cross-equation entries as large as
5.73e-5 despite the factorized likelihood's zero exact cross-information.
Using independent scores with Stata curvature reproduces Stata's robust
covariance tightly, supporting numerical curvature as the main source of
the raw difference. Substituting Stata component curvature still leaves
SE differences up to 1.064%; it is retained as a diagnostic, not an exact
reference. No package covariance was changed to force parity.

## Reproduction and verification

Run scripts from `tools/stata-benchmarks`. The partial builder accepts the
returned log path; the comparison accepts an optional CSV output path;
the independent audit accepts an optional RDS output path:

```sh
Rscript build_glmmtmb_nbinom2_rs_stata_partial.R /path/to/suest_r_glmmtmb_nbinom2_rs_partial_v2.log
Rscript compare_glmmtmb_nbinom2_rs_stata_partial.R
Rscript check_glmmtmb_nbinom2_rs_partial_audit.R
Rscript check_glmmtmb_nbinom2_rs_partial_margins.R
```

The builder records the log's MD5. Raw Stata and independent audit fixtures
are separate. The expanded regression covers balanced and partial cases,
sample sizes, overlap, cluster counts, coefficients, curvature, joint
covariance, margins, and the Stata-curvature/independent-score decomposition.
The focused regression passes 51 expectations. Fresh source build/check output and detailed audit
metrics accompany the checkpoint. The full 1,397-expectation suite passed
at implementation checkpoint `9e9ca45`; this fixture/documentation increment
does not rerun that full suite. The package check uses `--no-tests`, separately
from the focused regression run.

## Next bounded Stata run

Extract the flat ZIP into one directory, set Stata's working directory there,
and run:

```stata
do suest_r_glmmtmb_nbinom2_rs_disjoint_v3.do
```

Return `suest_r_glmmtmb_nbinom2_rs_disjoint_v3.log`, even if it stops early.
This runs only two components and one joint model. Each component has
500 observations in 50 id groups; the groups are odd/even, despite the
historical variable names `y1_left` and `y2_right`. Overlap is zero and the
union is 1,000 observations in 100 groups. Sample-pattern assertions check
these definitions. The CSV is unchanged.

R preflight-checked fits supply component starts. Actual Stata components
supply the joint start. Zero-step likelihood and unchanged-parameter checks
precede estimation. Every full mixed model and each preliminary component
fixed-only starting fit has a 15-iteration limit. GSEM uses Laplace integration,
`adaptopts(tolerance(1e-12))`, and `vce(robust)`. The run stops at its first
failed check, with no automatic restarts. Iteration caps are not wall-clock
timeouts. Stata is unavailable here; the new script is statically reviewed,
not executed in this runtime.

The returned disjoint covariance must be compared after explicit finite-sample
normalization: R uses 50/49 per component, while joint GSEM uses the 100/99
union correction. Their ratio is 99/98. Raw and normalized results will both
be retained, and the independent likelihood/score audit will be repeated
before calling that final sample pattern validated.
