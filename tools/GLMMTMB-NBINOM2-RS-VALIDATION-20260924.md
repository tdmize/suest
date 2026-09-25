# NB2 random slopes: completed Stata and independent validation

Development version 0.1.5.9000 supports the restricted glmmTMB NB2-log route
with one Gaussian random intercept and one correlated numeric random slope,
`(1 + x | id)`, and estimated constant dispersion. Returned Stata benchmarks
and independent likelihood/covariance audits now cover balanced samples,
partial overlap with 20 higher clusters, and disjoint odd/even groups.
No additional Stata run is required for this increment. No production R
implementation changed during this final validation step.

**Point estimates agree closely. Independent covariance agrees tightly with
R after matching parameter points, information coordinates, and cluster
corrections. Raw Stata covariance is numerically different and remains a
separate result; this is not a claim of exact raw GSEM covariance parity.**

## Final disjoint return

`suest_r_glmmtmb_nbinom2_rs_disjoint_v3.log` reports successful completion.
Both menbreg components and joint GSEM converged after one full-model
iteration each, under the 15-iteration cap. Preliminary fixed-only fits
took four or five iterations. Zero-step parameter changes were exactly zero;
component likelihood gaps were below 2.18e-7 and the joint gap was 4.55e-12.

Each component has 500 observations in 50 groups. The odd/even id groups
are disjoint despite the historical outcome names `y1_left` and `y2_right`.
The union has 1,000 observations and 100 groups; overlap is zero.

| Comparison with R | y1_left | y2_right | Joint, raw | Joint, correction aligned |
| --- | ---: | ---: | ---: | ---: |
| Largest absolute parameter difference | 2.994e-5 | 2.115e-5 | 2.994e-5 | 2.994e-5 |
| Absolute log-likelihood difference | 2.037e-7 | 1.761e-7 | 3.797e-7 | 3.797e-7 |
| Largest absolute covariance difference | 1.075e-4 | 3.858e-4 | 9.558e-4 | 3.544e-4 |
| Largest relative coefficient-SE difference | 0.1060% | 0.4223% | 0.9162% | 0.4119% |

Returned log likelihoods are -1053.32788051598, -1024.25346417252,
and -2077.5813446885, respectively. R applies 50/49 per component;
joint Stata GSEM applies 100/99 to the union. Multiplying Stata's raw joint
covariance by their ratio, 99/98, produces the aligned comparison. Both raw
and aligned matrices are retained; normalization does not remove the
remaining numerical-curvature difference.

Native marginal means differ by at most 5.09e-7 and x slopes by 1.84e-6.
Native margin SEs differ by up to 0.2839%. Using Stata's own covariance,
analytic margin variances reproduce the printed results within 1.84e-8,
or 4.11e-8 relative SE. Analytic gradients agree with independent Richardson
derivatives within 8.75e-11.

## Why the independent audit needed a coordinate refinement

The audit independently computes the bivariate Gaussian Laplace likelihood,
latent modes, group scores, and observed information without TMB likelihood,
score, or covariance code. Earlier audits inverted the observed information
directly in the exported Fisher-correlation coordinate `z = atanh(rho)`.
glmmTMB instead fits the native coordinate `t = sinh(z)` and transforms
the inverse information by the parameter Jacobian.

At an exact stationary point the two inverse-information calculations agree.
At a finite numerical optimizer endpoint, the fitted gradient is small but
nonzero, and the second-derivative chain rule must be retained. If `I_z` is
negative log-likelihood curvature and `g_z` is the log-likelihood gradient,
the native inverse information expressed in exported coordinates is

```text
B_native_delta = inverse(I_z + g_z[7] * rho * E77)
```

Here `E77` is zero except for a one in the correlation position. The sign
follows from `d2z/dt2 = -rho * (dz/dt)^2`. All other coordinates are shared.
This uses the independently computed gradient and information, not R's
native covariance as a replacement reference.

For disjoint samples, the original Fisher-coordinate sandwich differed from
R by 1.477e-6, exceeding the existing 1e-6 absolute covariance gate, although
its relative SE difference was only 1.353e-5. Including the native-coordinate
score term reduces that covariance difference to 3.872e-8 and relative SE
to 6.504e-7. **The tolerances were not widened.** Original Fisher-coordinate
and natural-variance-coordinate calculations remain recorded separately.

An additional diagnostic recomputes glmmTMB's native numerical Hessian at
steps .001, .0005, .00025, and .0001. The .001 result reproduces its supplied
covariance within 1.73e-14. The independent native-coordinate component
covariance agrees within 1.61e-8; at the two smallest numerical steps
(.00025 and .0001) the difference is below 3.49e-9. This separates the coordinate term from the
much smaller numerical differencing error. No package covariance was altered.

Balanced and partial audits were rerun with the same refinement. Their
original Fisher-coordinate comparisons still pass the earlier gates, and
their native-coordinate comparisons improve as expected.

## Independent covariance across all three cases

At R's fitted point, after matching its native information coordinates:

| Sample pattern | Largest absolute covariance difference | Largest relative SE difference |
| --- | ---: | ---: |
| Balanced, 100 groups | 2.256e-8 | 6.245e-7 (0.00006245%) |
| Partial overlap, 20 higher clusters | 1.895e-8 | 5.913e-7 (0.00005913%) |
| Disjoint, 50 groups per component | 3.872e-8 | 6.504e-7 (0.00006504%) |

At each Stata joint parameter point, using Stata's model-based curvature
with independently computed score cross-products reproduces its raw robust
covariance within 5.09e-7 across all cases, and relative SE within 5.82e-6.
All comparisons pass the unchanged 1e-6 absolute covariance and 1e-4 relative
SE gates. Independent disjoint log likelihood differs from Stata by less
than 3.81e-7. Score step-halving changes are below 1.07e-9; inverse-curvature
step changes are below 2.48e-9.

The raw disjoint Stata covariance has cross-equation entries as large as
6.785e-5; its model-based covariance has cross-entries up to 3.772e-5.
Recovering the score cross-product matrix from those two reported matrices
leaves cross-entries of at most 5.65e-11. Independent score cross-products,
independent covariance, and R's joint covariance have exactly zero cross-
blocks for disjoint groups. These are recovered cross-product matrices,
not uniquely recovered score vectors. The factorized likelihood and the
successful Stata-curvature/independent-score reconstruction locate the raw
discrepancy chiefly in numerical curvature. Substituting Stata component
curvature still leaves SE differences up to 1.034%, so that substitution
is diagnostic rather than an exact reference.

## Scope, limits, and verification

The joint parameters include fixed effects, `log_phi` (log NB2 size), both log
SDs, and the Fisher correlation. Population means and slopes integrate both
Gaussian random effects and retain all nuisance covariance. The route
remains restricted to ML, one grouping factor, one bare numeric random
slope, unstructured intercept/slope covariance, no weights or offsets,
no zero inflation, and estimated constant dispersion. Unsupported mappings,
priors, constrained covariance structures, and extra slopes remain rejected.
NB1 and binomial-logit random slopes remain unsupported. See the
implementation's support tests for exact guards.
Well-scaled predictors remain important: suest preserves native glmmTMB
covariance, including any upstream numerical-curvature limitations.

Earlier balanced and partial reports preserve their original results and
tolerance decisions. In particular, partial y2's log-size endpoint differs
by 7.28e-5 (0.061% of its native SE), and its printed margin derivatives
have a small documented numerical difference. Their existing point/margin
gates are retained. Disjoint coefficients use the balanced 5e-5 absolute
gate; disjoint native margin variances use its 1e-7 absolute gate.

The final regression covers all three sample patterns, six components,
parameter transformations, native curvature, integrated margins, joint
covariance, cluster normalization, and disjoint zero cross-blocks. All 78
focused expectations pass without warnings or skips. The full suite passes
1,475 expectations with no failures or skips and three test warnings, the
same warning count as the earlier implementation checkpoint. Source build
and full `R CMD check --no-manual` succeed with **Status: OK**. The handoff
and captured output record the environment and results. This checkpoint
remains development version 0.1.5.9000;
no release or remote push is implied.

## Reproduction

Run scripts from `tools/stata-benchmarks` in the repository recovered from
the complete-history bundle. The flat ZIP also includes the R reference,
returned logs, raw/audited fixtures, and calculation outputs. The source
package contains the installable package and regression fixtures.

```sh
Rscript build_glmmtmb_nbinom2_rs_stata_disjoint.R /path/to/suest_r_glmmtmb_nbinom2_rs_disjoint_v3.log
Rscript compare_glmmtmb_nbinom2_rs_stata_disjoint.R
Rscript check_glmmtmb_nbinom2_rs_disjoint_audit.R
Rscript check_glmmtmb_nbinom2_rs_disjoint_margins.R
Rscript check_glmmtmb_nbinom2_rs_disjoint_curvature.R
```

The balanced and partial audit scripts reproduce their updated independent
fixtures. `run_glmmtmb_nbinom2_rs_reference.R` rebuilds the six-model R
reference and optimizer preflight table. Raw fixtures preserve source-log
MD5 hashes and are never overwritten with reconstructed covariance.
