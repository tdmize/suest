# NB2 random intercept and slope: development checkpoint

Latest checkpoint: all three Stata sample patterns have been audited.
See [the final validation report](GLMMTMB-NBINOM2-RS-VALIDATION-20260924.md)
for results, limits, and the native-coordinate curvature comparison. The
initial development record below preserves earlier scope and pending status.

Update after the returned balanced-v1 log: see the
[balanced validation report](GLMMTMB-NBINOM2-RS-BALANCED-VALIDATION-20260924.md).
Balanced point estimates and independent covariance audits are complete;
raw Stata numerical-curvature differences are documented. Partial/disjoint
Stata validation remains pending. The original checkpoint below records the
pre-return state.

This increment extends the completed Poisson random-slope route to NB2-log
models from glmmTMB 1.1.14+. R implementation and verification are available;
**Stata NB2 random-slope validation is pending**. This is development version
0.1.5.9000, not a new release. The parent checkpoint is
`524b08ce6db7de328b7351055328a38cbb32490f`; work is on
`dev/glmmtmb-nbinom2-rs`. No remote push is part of this checkpoint.

## Support contract

Unweighted maximum likelihood, one grouping variable, one Gaussian random
intercept and one correlated numeric slope, `(1 + x | id)`, with an
unstructured 2-by-2 covariance matrix. The slope must be one untransformed
numeric column with a syntactically valid name. Fixed-effect factors and
interactions are supported. The slope can be absent from the fixed formula,
but response-prediction data must contain it. NB2 dispersion must be estimated
and constant with `dispformula = ~1`; conditional variance is
`mu + mu^2/phi`.

The type is `glmm_nbinom2_rs`. Systems must contain models of that same type.
Default grouping-level or higher nested clustering, identical samples,
partially overlapping samples, and disjoint groups are covered. As in the
existing panel routes, equation-specific cluster corrections apply; disjoint
models with 50 groups use 50/49, not the union's 100/99.

Excluded: weights (including unit weights), offsets (including all-zero
offsets), zero inflation, dispersion regressors or fixed dispersion, REML,
priors, parameter maps/constraints, diagonal covariance, multiple slopes,
slope-only terms, transformed/factor slopes, multiple grouping variables,
NB1, and logit random slopes. Convergence, finite positive native curvature,
and interior random covariance are required. The covariance boundary guard
uses a centered, standardized slope basis. The NB2 likelihood must improve
numerically over the zero-random-effect likelihood at the same fixed effects
and dispersion. These are numerical support checks, not significance tests
or proof of a global optimum.

## Parameters and predictions

Joint parameter order: fixed effects, `log_phi`, `log_sd_intercept`,
`log_sd_slope`, `atanh_rho`. Native glmmTMB's final scaled-correlation parameter
`t` maps to Fisher correlation by `asinh(t)`; its covariance Jacobian is
`1/sqrt(1+t^2)`, and scores transform inversely. The full dispersion and
random-covariance blocks and all cross-covariances are retained.

Write the random covariance as `S = [[A,C],[C,B]]`. Population response means
are `exp(X beta + (A + 2*C*x + B*x^2)/2)`; link predictions are `X beta`.
The exact Gaussian mean has no direct `log_phi` derivative. Dispersion still
affects estimation and joint uncertainty through the full information and
score matrices. A continuous x effect includes `C + B*x` in addition to its
fixed-effect derivative. Marginaleffects predictions, slopes, factor
comparisons, and cross-model contrasts use the full joint covariance.

## Numerical checks and scaling limitation

Verification on R 4.5.2, glmmTMB 1.1.14, TMB 1.9.25, and marginaleffects 1.0.0:
63 focused expectations pass, with no warnings, failures, or skips. The full
installed-package suite passes 1,397 expectations with no failures or skips
and three existing warnings. Source build succeeds. The initial full
`R CMD check --no-manual` returns zero and records the complete passing test
output, but its main log is truncated. A final `R CMD check --no-manual
--no-tests` verifies all remaining gates without repeating tests: **Status:
OK**, zero errors, warnings, or notes. Both the full test output and complete
final check log accompany the checkpoint. Static review of the R changes and
bounded Stata gate found no remaining blockers; Stata execution is unverified.

The focused regression suite independently computes bivariate Laplace
likelihoods by optimizing the two latent effects and evaluating their
analytic Hessian. Finite differences verify all seven score coordinates on
three groups, including dispersion and random correlation. Further checks
cover the native covariance Jacobian, analytic mean/effect gradients,
coefficient replacement, sample alignment, clustering, unsupported
specifications, and slope data requirements.

An extreme-units test multiplies x by 10,000. glmmTMB 1.1.14/TMB 1.9.25's default
native numerical Hessian uses an absolute finite-difference step of .001,
which can exceed the scaled slope coefficient. In this fixture, both fits
report convergence and a positive-definite Hessian, but the mean's robust SE
changes from about .303 to 1.395 with the uncorrected supplied covariance.
TMB gradients and the likelihood remain consistent. A native `sdreport`
explicitly supplied with a scale-aware Hessian restores equivalent inference;
halving its steps verifies curvature stability. The regression test makes
that input correction explicit, compares native covariance under the unit
transformation, and retains its original mean/SE tolerances.

**suest does not silently replace native covariance or certify its numerical
accuracy. Fit predictors in well-scaled units and check sensitivity to
rescaling.** A general curvature-consistency diagnostic is outside this
increment. This limitation concerns upstream curvature, not a direct NB2
dispersion effect on the population mean.

## Benchmark and bounded Stata gate

The reproducible seed-2432 fixture has 1,000 observations, 100 groups of 10,
and 20 higher clusters. Partial samples have 975 and 980 observations with
955 shared rows. Disjoint samples use odd/even groups, 500 observations and
50 groups each. No data were changed after inspecting Stata results; none
have yet been returned for this increment.

Six R component fits passed a deliberately stronger fixture preflight:
smallest random-covariance eigenvalue above .02, absolute correlation below
.85, NB2 size between .5 and 20, and native parameter SE below .8. Observed
minimum eigenvalues are .0952–.1880, correlations .2291–.6402, and size
2.47–2.85. Separate nlminb and BFGS fits agree within 3.64e-8 in log
likelihood and 3.94e-5 in parameters. The saved reference contains six
components, all three joint R systems, native component margins, and balanced
joint predictions/slopes. These are R results, not Stata validation.

Run `suest_r_glmmtmb_nbinom2_rs_balanced_v1.do` with the CSV and generated
start-values do-file in the current Stata directory. It performs:

1. A no-estimation parameter-template lookup with `startvalues(zero)`, and a zero-iteration likelihood
   check for each balanced `menbreg` component, using explicitly labeled R
   Laplace seeds. Stata `lnalpha = -log_phi`; variance/covariance seeds use
   natural scales. A changed start or likelihood gap over 1e-4 stops the run.
   Correction from the returned log: `startvalues(zero)` did not suppress
   menbreg's preliminary fixed-only fits. They finished in four or five
   iterations. Revision 2 explicitly requests the documented fixed-only
   starting-fit limit through `startvalues(fixedonly, iterate(15))`.
2. Two component Laplace fits, each with at most 15 outer iterations.
   Successful estimates, native covariance, and native marginal mean/slope
   results are printed and checked for finite scalar results and nonnegative
   marginal variances.
3. One joint Laplace `gsem`, seeded by the successful **Stata** components.
   Each equation's intercept/slope covariance is free; all four latent
   cross-equation covariances are zero. The robust VCE retains dependence
   between observed outcomes. Another zero-step check precedes a maximum of
   15 outer iterations. Both robust and model-based covariance are printed
   so any discrepancy can be audited directly.

The joint GSEM uses `adaptopts(tolerance(1e-12))`, following the earlier Poisson
investigation. Components retain menbreg's default inner integration controls.
The outer iteration cap is not a wall-clock timeout. There
are no cold restarts, adaptive-quadrature comparisons, or partial/disjoint
Stata fits in this first gate. The script stops on the first failed check.
Zero-step r(430) is allowed only when a finite unchanged parameter vector
reproduces the expected likelihood; a full fit requires return code zero and
`e(converged) == 1`. Completion markers report fit success, not numerical
cross-language agreement. Stata is unavailable in the development runtime;
the script has been reviewed but not executed here.

Return `suest_r_glmmtmb_nbinom2_rs_balanced_v1.log`, including a stopped or
failed log. Compare component likelihoods, all seven transformed parameters,
native covariance, native margins, and joint robust/model-based covariance
before widening the Stata gate. Reuse the prior independent-Laplace and
curvature-audit approach if raw GSEM covariance differs. Do not describe NB2
random slopes as Stata-validated until that comparison is complete.

Reproduction scripts live in `tools/stata-benchmarks`:
`make_glmmtmb_nbinom2_rs_crosslang_data.R` and
`run_glmmtmb_nbinom2_rs_reference.R`. The latter writes the reference RDS,
preflight CSV, and R-seeded Stata start file. The delivered archive includes
the test/check output and a full-history local Git bundle.

Primary specification references:
[glmmTMB NB2 family](https://glmmtmb.github.io/glmmTMB/reference/nbinom2.html),
[Stata menbreg](https://www.stata.com/manuals/memenbreg.pdf),
and [Stata meglm estimation options](https://www.stata.com/manuals/memeglm.pdf).
