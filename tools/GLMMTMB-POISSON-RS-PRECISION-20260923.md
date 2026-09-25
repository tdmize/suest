# Poisson random slopes: revision-4 experiment (historical)

V4 has returned. The completed decomposition and next run are recorded in
`GLMMTMB-POISSON-RS-CURVATURE-20260923.md`; the experiment plan below is retained
for provenance.

## Revision-3 result

Using the returned Stata component estimates resolved the balanced joint
convergence problem. The fit converged after one iteration with 800
observations and 100 groups. Its log likelihood was -2812.42881818855,
only 1.17e-7 above the sum of the component log likelihoods. The zero-step
evaluation reproduced that sum within 4.74e-8 without changing any supplied
parameter values. The zero-step "convergence not achieved" warning was
expected because that diagnostic deliberately allowed no iterations.

The full 16-parameter Stata output (including four fixed loadings) is retained
in `glmmtmb-poisson-rs-stata-balanced-v3.rds`. A 12-by-16 Jacobian transforms
all covariance entries to the R coefficient/log-SD/Fisher-correlation scale.
The fixture marks covariance validation as pending; it is not a relaxed
covariance acceptance test.

Comparison with R:

- Maximum absolute coefficient difference: 7.26e-6.
- Maximum absolute covariance difference: 2.07405e-4.
- Maximum absolute cross-model covariance difference: 8.62211e-5.
- Maximum relative standard-error difference: 0.005665, or **0.5665%**,
  for the second model's intercept.

Both calculations use 100 groups and the 100/99 finite-sample correction.
The differences vary by parameter; a single correction factor does not explain
them. Coefficient and likelihood agreement alone does not validate the robust
covariance. No covariance tolerance was widened to accept this result.

## Independent check of the R covariance

`check_glmmtmb_poisson_rs_independent_sandwich.R` computes all 100 group
likelihoods per outcome using an independent bivariate Laplace calculation.
It solves each conditional mode by Newton iterations with backtracking,
includes the covariance and conditional-Hessian determinants, and obtains
all six group scores by central finite differences. Halving the step from
1e-5 to 5e-6 changes any score by at most 1.60e-9.

The calculation uses the **returned Stata component OIM covariance matrices**
as bread, combines both models' group influences, and applies 100/99. Its
covariance differs from R by at most 3.45e-7; its maximum relative SE
difference is 6.63e-6 (0.000663%). It retains approximately the same 0.566%
SE discrepancy with joint GSEM. This supports the R sandwich and motivates
checking the joint numerical calculation.

A further diagnostic: the zero-step joint output reported a gradient entry
of 31.234 for `var(S2[id])`. An independent finite difference at the supplied
component estimates gives approximately 6.85e-8 in that natural variance
coordinate. The final joint gradient is much smaller. A zero-step gradient
can reflect initialization or reporting behavior, so this discrepancy alone
does not establish an error in the final derivatives. It motivates further
diagnostics without identifying a particular Stata implementation defect.

## Revision-4 experiment

The only changed numerical setting is **`adaptopts(tolerance(1e-12))`**, from
the documented default 1e-8. Stata documents this control as applicable to
Laplace integration and documents numerical derivatives for GSEM's Laplace
route. Tighter internal integration is a hypothesis to test, not an already
verified fix.

The model, CSV, supplied component estimates, Laplace approximation, robust
VCE, outer convergence tolerances, zero-step checks, and 15-iteration maximum
remain unchanged. The new log also prints `e(V_modelbased)`, allowing the
joint bread to be compared with the block-diagonal component bread. This
does not require an additional model fit.

Extract the flat ZIP in the working directory and run:

```stata
do suest_r_glmmtmb_poisson_rs_balanced_v4.do
```

Return **suest_r_glmmtmb_poisson_rs_balanced_v4.log**. This is still only the
balanced model, with no component refits and no automatic fallback fits.
Partial and disjoint joint checks remain deferred until this covariance
question is resolved. No R installation or Git operation is needed.

## Reproduction and verification

From `tools/stata-benchmarks`, use
`build_glmmtmb_poisson_rs_stata_balanced.R <returned-log>` to record the
converged joint output, then `compare_glmmtmb_poisson_rs_stata_balanced.R`
with the chosen fixture path to compare it against the R reference. The
independent sandwich script also runs from that directory.

The v3 parser, full-matrix comparison, and independent audit ran successfully
under R 4.5.2. The v4 do-file was reviewed; it cannot be executed here because
Stata is unavailable. No R implementation changed, and the full package
suite was not repeated for this diagnostic increment. No remote push was made.

Primary documentation checked:

- https://www.stata.com/manuals/semgsemestimationoptions.pdf
- https://www.stata.com/manuals/semgsem.pdf
