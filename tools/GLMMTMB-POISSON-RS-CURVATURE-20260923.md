# Poisson random slopes: balanced audit and v5 plan (historical)

V5 has returned and the remaining audits are complete. Current evidence and
verification are in `GLMMTMB-POISSON-RS-VALIDATION-20260923.md`; the balanced
audit and v5 plan below are retained for provenance.

## Returned balanced v4 result

The balanced joint model converged after one iteration with 800 observations
and 100 groups. Log likelihood is -2812.42881825605, within 7.56e-7 of the R
component sum. The maximum transformed coefficient difference is 4.57e-5.
The supplied component estimates again avoided the stalled cold optimization.

Tightening internal integration tolerance from 1e-8 to 1e-12 reduced the
maximum relative robust standard-error difference from 0.5665% to 0.1796%.
The maximum covariance difference is 5.01e-5, including 4.98e-5 in a
cross-equation block. Raw joint covariance therefore remains a numerical
diagnostic, not an exact-equality reference. No package code or existing
acceptance tolerance was changed to accommodate it.

## Where the remaining difference enters

The returned `e(V_modelbased)` permits a decomposition at the exact v4 joint
parameter values. Independently computed bivariate Laplace group likelihoods
supply group scores and a block-diagonal observed-information inverse. Scores
use central differences at 1e-5 and 5e-6. Curvature uses central second
differences with Richardson extrapolation at (.002, .001) and (.001, .0005).
Halving the steps changes any score by at most 1.78e-9 and any inverse-
information entry by at most 1.04e-9.

For a sandwich `V = c B M B'`, `c = 100/99`, the returned robust covariance
and model-based covariance B recover the Stata score cross-product M by
matrix inversion. This is diagnostic algebra; it does not replace Stata's
reported results or change package behavior.

| Calculation | Maximum relative SE difference |
|---|---:|
| Raw Stata joint vs R | 0.1796% |
| Independent curvature and scores vs R | 0.00129% |
| Stata curvature with independent scores vs raw Stata | 0.000363% |
| Independent curvature with recovered Stata scores vs fully independent | 0.000363% |

Thus the residual is located chiefly in the returned joint model-based
covariance. Its cross-equation entries reach 3.37e-5, whereas the likelihood
factorizes across equations and independent curvature has zero cross-blocks.
The recovered score cross-product agrees closely with independent scores.
This supports the R covariance and proceeding to the other sample patterns;
it does not establish a particular defect in Stata's numerical implementation.

The v4 point is close to, but not exactly at, a stationary point. Observed
information consequently depends slightly on parameter coordinates. The audit
also applies the total-score chain-rule correction for natural variance/
covariance coordinates. It changes independent SEs by at most 0.00214% and
leaves the zero cross-equation curvature unchanged. This sensitivity does not
alter the conclusion and makes no assumption about Stata's internal optimizer
coordinates.

The raw v4 fixture preserves the complete reported coefficients, robust
covariance, model-based covariance, gradient, and transformation Jacobian.
Its generic pending covariance flag prevents treating successful log parsing
as a covariance acceptance decision. This report records the subsequent audit.

## Next Stata run

Extract the flat ZIP into the Stata working directory and run:

```stata
do suest_r_glmmtmb_poisson_rs_remaining_v5.do
```

Return **suest_r_glmmtmb_poisson_rs_remaining_v5.log**, including if it stops.
The run contains only these two joint models:

| Case | Component outcomes | Union / overlap | Robust clusters |
|---|---|---|---|
| Partial overlap | y1_partial, y2_partial | 800 / 755 observations | 20 higher clusters |
| Disjoint | y1_left, y2_right | 800 / 0 observations | 100 union groups |

Both reuse the successful component estimates from the returned v2 log.
Each first evaluates the supplied parameters with zero optimization steps,
checks the joint likelihood against the component sum, then permits at most
15 joint optimization iterations. The zero-step diagnostic may print
"convergence not achieved"; acceptance requires the subsequent fitted model
to converge. The script stops at the first failed case. It uses the v4
integration tolerance and prints both covariance matrices. Completion markers
establish convergence and finite output, not cross-language agreement.

The data, likelihood, random-effect structure, and clustering specification
match the original benchmark. No component re-estimation, balanced rerun, or
automatic fallback fit is included. No R installation or Git operation is
needed for this check.

For the disjoint comparison, R applies 50/49 within each component; GSEM's
union correction is 100/99. The covariance-scale normalization to compare
GSEM to R is therefore 99/98. Keep raw and normalized results separately.
The R disjoint cross-equation covariance must remain zero.

## Reproduction and scope

From `tools/stata-benchmarks`:

1. Run `build_glmmtmb_poisson_rs_stata_balanced.R <returned-v4-log>`.
2. Run `compare_glmmtmb_poisson_rs_stata_balanced.R` with the v4 fixture path
   and an optional output CSV.
3. Run `check_glmmtmb_poisson_rs_balanced_curvature.R` (defaults to v4), or
   supply the fixture path and an output RDS path.

The parser, covariance comparison, and complete independent curvature audit
ran successfully under R 4.5.2. Independent review reproduced the raw matrix
parsing, Jacobian, and decomposition. The new Stata script was reviewed and
its four starting vectors and sample counts were checked against the saved
fixtures; it cannot be executed here because Stata is unavailable.

The R implementation and its existing tests are unchanged. The full package
suite was not repeated for this diagnostic increment. Partial and disjoint
joint validation remains pending; random slopes for NB2 and logit remain
out of scope until this increment is complete. No remote push was made.
