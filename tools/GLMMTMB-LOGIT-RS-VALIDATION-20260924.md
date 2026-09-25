# Binomial-logit random slopes: development validation

Date: 2026-09-24; validation updated 2026-09-25. Version: 0.1.5.9000.
**R development support independently validated in balanced, partial-overlap,
and disjoint samples. Strict raw GSEM covariance parity is not established.**
Raw joint SEs differ from R by up to 4.23%, 1.14%, and 1.57%, respectively,
after aligning cluster corrections. Numerical curvature accounts for the
main differences. One disjoint reconstruction diagnostic remains outside its
original absolute bound. See the
[balanced report](GLMMTMB-LOGIT-RS-BALANCED-VALIDATION-20260924.md),
[partial report](GLMMTMB-LOGIT-RS-PARTIAL-VALIDATION-20260924.md), and
[disjoint report](GLMMTMB-LOGIT-RS-DISJOINT-VALIDATION-20260925.md).
This increment adds Bernoulli binomial-logit models from glmmTMB 1.1.14 or later
with one correlated Gaussian random intercept and one bare numeric random
slope, `(1 + x | id)`. It follows the completed NB2 random-slope checkpoint.

## Parameter and prediction contract

The joint system retains fixed effects, `log_sd_intercept`, `log_sd_slope`,
and `atanh_rho`. If glmmTMB's native correlation coordinate is `t`, the exported
Fisher coordinate is `asinh(t)`, with covariance Jacobian `1/sqrt(1+t^2)`.
The full native covariance is retained under this transformation. Group
scores transform inversely; fixed/nuisance cross-covariances are retained.
Models in a system must share the supported type. Panel clusters are the
default; higher clusters must contain whole panels. Disjoint model blocks
use the existing component-cluster normalization and have zero cross-blocks.

Write the random covariance entries as `A = var(intercept)`,
`B = var(slope)`, `C = cov(intercept,slope)`. Response predictions are

```
eta = X beta
v(x) = A + 2 C x + B x^2
mu = E[plogis(eta + sqrt(v(x)) Z)],  Z ~ N(0,1).
```

The implementation evaluates `v` through a nonnegative sum of squares.
It uses a cached 320-point normalized normal Gauss-Hermite rule for SD <= 3. For larger
SDs it integrates the equivalent logistic-density convolution
`integral dlogis(t) * pnorm((eta-t)/SD) dt`, avoiding an increasingly sharp
normal-space integrand. Link predictions remain `X beta`. The existing
random-intercept-only predictor is unchanged.

If `m_k` is the expectation of the kth derivative of the logistic function,
the continuous slope is `beta_effective * m_1 + (C+B*x)*m_2`.
Its uncertainty includes changes in both the fixed predictor and the random
variance, and all covariance-parameter uncertainty. Factor comparisons
recompute the fixed design including interactions.

Bernoulli outcomes (including two-level factors) only. Grouped binomial,
weights, offsets, zero inflation, priors, REML, maps/constraints, diagonal
covariance, transformed/factor slopes, multiple slopes, slope-only terms,
and extra grouping factors are unsupported. Near-singular covariance is
rejected in a centered, standardized predictor basis. glmmTMB native
numerical covariance remains sensitive to poor scaling; use well-scaled data.

## Independent R evidence

`tests/testthat/test-glmmtmb-logit-rs.R` checks:

- A separately implemented bivariate Laplace likelihood with polished panel
  modes agrees with glmmTMB; finite-difference panel scores agree with its
  transformed full scores. Native covariance agrees after the full Jacobian.
- Integrated probabilities agree with independent adaptive normal integration,
  including SDs around the algorithm switch and very large SDs. A separate
  mixed-derivative test crosses the switch, since level accuracy alone does
  not establish slope-gradient accuracy.
- Independent analytic gradients verify population probability, continuous
  slope, factor comparison, and cross-equation contrast SEs. The fixed model
  includes a numeric-by-factor interaction. Perturbing each nuisance
  parameter changes only the intended equation's response predictions.
- Native coefficient influences reproduce balanced, higher-cluster,
  partially overlapping, and disjoint system covariance and corrections.
- Unsupported structures and boundary fits are rejected; factor outcomes and
  a slope appearing only in the random-effects formula are exercised.

Nested finite differences in marginaleffects' default slope calculation can
show small cancellation errors. During the initial interaction diagnostic
(with the 160-point rule), default slope SEs
differed from the analytic reference by 3.25e-6 and 1.88e-6 (about 0.021% and
0.011%). With `numderiv = list("fdcenter", eps = 1e-4)`, differences were below
3e-10; eps 1e-3 also agreed within 5e-10. Use centered coefficient differences
and check step-size sensitivity for precise slope inference. The package
does not silently change marginaleffects' global numerical settings.

Review identified that the initial 160-point quadrature differed from the
convolution by about 8.6e-12 at SD=3. Although acceptable for probability
levels, this produced a 5.38e-5 mixed-derivative error in a boundary-crossing
case. A regression test failed before the fix. The final 320-point symmetric, normalized
rule is tested against an independently integrated analytic mixed derivative
with a 2e-8 absolute bound. Enforcing the normal rule's exact node/weight
symmetry removes eigensolver asymmetry that otherwise amplified at the switch. The existing RI quadrature remains unchanged.

## Fixed benchmark and original balanced Stata gate

The seeded fixture has 1,200 rows in 100 groups, with 12 observations per
group and well-scaled predictors. Preflight covers both outcomes under all
three sample patterns, not just the first balanced gate. Required covariance
minimum eigenvalue > .02, absolute correlation < .85, and maximum exported
parameter SE < .8. Precision is assessed in the exported Fisher coordinate;
the native scaled-correlation coordinate can have a larger SE. Both values
are recorded. A second BFGS optimizer must agree within 1e-5 log likelihood
and 1e-4 exported parameter units. Data and R start values are explicitly
identified as generated references, not returned Stata estimates.

Checked R fits (native fits and both optimizer endpoints):

| Outcome | Minimum covariance eigenvalue | Correlation | Maximum exported SE | Optimizer parameter gap |
|---|---:|---:|---:|---:|
| y1 | 0.258236 | 0.590771 | 0.242011 | 4.77e-06 |
| y2 | 0.676053 | -0.047645 | 0.193203 | 5.09e-06 |
| y1_partial | 0.235359 | 0.609523 | 0.254455 | 4.37e-06 |
| y2_partial | 0.663910 | -0.049262 | 0.195578 | 6.47e-06 |
| y1_left | 0.111592 | 0.784290 | 0.516354 | 9.17e-05 |
| y2_right | 0.673331 | 0.265833 | 0.265652 | 6.85e-06 |

The largest two-optimizer log-likelihood difference is 1.95e-8.

The first do-file is `suest_r_glmmtmb_logit_rs_balanced_v1.do`:

1. Fit two balanced `melogit` components with Laplace integration and
   unstructured random intercept/slope covariance, starting from checked R
   parameters. Preliminary fixed-effects fits and full fits are capped at
   15 iterations. Zero-step likelihood checks must agree within 1e-4 and
   preserve the requested starting parameters within 1e-10.
2. Print native integrated probability and slope results with 40 quadrature
   points per random-effect dimension. These predictions integrate over the
   random effects; the fitted likelihood remains Laplace.
3. Fit one joint Bernoulli-logit GSEM from the actual Stata component
   estimates. Both within-equation random covariances are free; all four
   cross-equation latent covariances are fixed at zero. The robust VCE
   retains score dependence across outcomes. Check the joint zero-step
   likelihood against the component sum before the capped full fit.
4. Stop at the first failed check, with no automatic restarts. A completion
   marker proves successful fits only; it does not establish R/Stata parity.

No Stata runtime is available here. The user returned all three sample gates;
each component and joint fit converged after one full-model iteration.
Parameters and native integrated margins have been checked. Independent
covariance audits meet the unchanged R acceptance limits in every sample
pattern. Disjoint raw covariance reconstruction retains an explicit unresolved
absolute-bound diagnostic. No further Stata run is needed for independent R
validation; a targeted derivative export would be needed to resolve that
remaining Stata-specific diagnostic. Raw and audited results stay separate.

Syntax references: [melogit](https://www.stata.com/manuals/memelogit.pdf),
[melogit postestimation](https://www.stata.com/manuals/memelogitpostestimation.pdf),
[meglm](https://www.stata.com/manuals/memeglm.pdf),
and [gsem](https://www.stata.com/manuals/semgsem.pdf).
The Stata script follows the previously returned NB2 random-slope gate,
including full parameter-stripe matching and top-level Mata helper functions.

## Original implementation verification

- Focused new-family tests: 61 assertions passed in 7 tests; no failures,
  warnings, errors, or skips.
- Full source-package suite: 1,536 passed, 0 failed, 0 skipped, with the same
  three existing test warnings as the previous 1,475-pass NB2 checkpoint.
- R CMD build succeeded. R CMD check --no-manual returned **Status: OK**,
  including tests, examples, and vignette rebuilding.
- Independent read-only review found the integration-switch issue described
  above; the added regression test and exact-R diagnostics verify its fix.
- Environment: R 4.5.2, glmmTMB 1.1.14, TMB 1.9.25, marginaleffects 1.0.0.

## Reproduction

From the repository's `tools/stata-benchmarks` directory, with the development
source and suggested dependencies installed:

```
Rscript make_glmmtmb_logit_rs_crosslang_data.R
Rscript run_glmmtmb_logit_rs_reference.R
```

The generated reference contains component parameters/native covariances,
independently integrated native margins, system coefficients/covariances,
preflight diagnostics, and marginaleffects estimates. The checkpoint includes
full package-check and focused-test outputs. No remote push or release is
part of this increment.
