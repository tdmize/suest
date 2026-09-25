# Poisson random intercept and slope: development checkpoint

The next increment after validated NB2 random intercepts adds one numeric
random slope to the existing glmmTMB Poisson-log route. The R implementation
and validation against Stata component fits and joint score cross-products
are complete. Raw GSEM covariance discrepancies are documented in the
validation report. This remains development
version 0.1.5.9000, based on NB2 checkpoint `a2fbc1e`.

## Model and parameter contract

Supported: unweighted Poisson-log maximum likelihood using glmmTMB 1.1.14+
with `(1 + x | id)`, one grouping variable, and one unstructured 2-by-2
Gaussian random-effect covariance. `x` is a bare numeric data column with a
syntactically valid name. It may be absent from the fixed-effects equation.
Fixed effects may include factors and interactions under the existing adapter.

The model type is `glmm_poisson_rs`. Parameter order is conditional fixed
effects, `log_sd_intercept`, `log_sd_slope`, `atanh_rho`. These last three names
are reserved. The two models in a system must both use this route; different
fixed-effect formulas and sample patterns are allowed. Default clustering is
by random-effect group; higher clusters must contain whole groups.

Excluded: weights (including explicit unit weights), offsets (including
zero-valued offsets), zero inflation, dispersion submodels, REML, priors,
mapped/constrained parameters, diagonal covariance, multiple slopes, slope-only
terms, transformed/factor random slopes, additional grouping levels, and
random slopes for other families. Logit/NB2 random intercepts remain supported
under their existing contracts.

## Scores and covariance

For native glmmTMB correlation coordinate `t`, `rho = t/sqrt(1+t^2)`.
Fisher correlation is `atanh_rho = asinh(t)`. The transformation Jacobian is
identity except for its last diagonal, `1/sqrt(1+t^2)`. Native full covariance
is checked against retained TMB covariance and then transformed as `J V J'`.
The full group score is transformed by the inverse Jacobian. The sandwich
therefore preserves fixed-effect/variance/correlation uncertainty, including
cross-equation blocks. Group influence is aligned once per group; each model
uses its own cluster correction, as in the existing random-intercept routes.

Tests independently evaluate a bivariate Laplace likelihood. Conditional modes
are solved using the Poisson gradient and information, with Newton polishing
before finite differences. Group score derivatives verify all parameter
coordinates. Native group influences separately verify balanced, higher-cluster,
partial-overlap, and disjoint covariances. Disjoint cross-blocks must be zero.

## Predictions and effects

Write `A = Var(b0)`, `B = Var(b1)`, `C = Cov(b0,b1)`. Response predictions
integrate over a new Gaussian random effect:

`mu(x) = exp(X beta + (A + 2*C*x + B*x^2)/2)`.

Link predictions remain `X beta`. For a simple fixed linear term in x,
`d mu/dx = mu * (beta_x + C + B*x)`. Thus replacing a random intercept with a
random slope changes both the mean and its covariate derivative. Predictions
use a sum-of-squares expression for the variance to reduce cancellation.

Analytic tests check average means, slopes, factor comparisons, their standard
errors, and a cross-model mean contrast. Parameter replacement tests perturb
both log SDs and Fisher correlation. A slope-only predictor still must be
provided in `newdata` for response predictions.

## Numerical support and unit invariance

Fits must converge with a positive-definite Hessian. For the numerical
random-covariance check, transform `(b0,b1)` into the basis corresponding to
`(x-mean(x))/sd(x)` on the estimation sample. The smallest covariance
eigenvalue must exceed `sqrt(.Machine$double.eps)*max(1, largest_eigenvalue)`.
This detects numerical zero/singularity of random variation on the linear
predictor scale. It is not a significance test or a complete weak-identification
diagnostic.

Independent review caught a first guard that compared raw intercept and slope
variances, which have different units. The corrected guard uses the standardized
design basis. A regression changes x units by a factor of 10,000 and verifies
unchanged population means and standard errors. Its scaled fit uses a warm
start and scaled BFGS coordinates to avoid unrelated native optimizer overflow.

## Stata gate and preflight

The fixed CSV contains 800 rows in 100 groups, with eight observations each.
All six R components pass a deliberately stricter benchmark preflight:
minimum raw covariance eigenvalue above .02, absolute correlation below .85,
and native parameter SE below .8. These are fixture checks, not package limits.
Observed minimum eigenvalues range from .0772 to .1477 and correlations from
.1115 to .5239. Disjoint models use odd/even IDs, 50 groups each.

The do-file fits six Laplace `mepoisson` components with
`covariance(unstructured)` and three Laplace `gsem` joint systems. In each
GSEM equation, the intercept/slope covariance is free; all four covariances
between the two equations' latent effects are fixed to zero. Robust scores
retain cross-equation sampling covariance. The joint cases are balanced,
partial overlap with 20 higher clusters, and disjoint groups. Observation
overlaps are 800, 755, and 0. Balanced components also print native marginal
means and slopes and their covariance.

For the returned Stata variance coordinates `(A,B,C)`, the Fisher-correlation
Jacobian row is `(-rho/(2*A), -rho/(2*B), 1/sqrt(A*B))/(1-rho^2)`;
the log-SD Jacobian diagonals are `1/(2*A)` and `1/(2*B)`. Parameter names and
order must be checked against the actual returned log before constructing
fixtures. No Stata equality or convergence is claimed before that return.
This first gate uses Laplace only, avoiding a change of approximation.

## Reproduction and next action

Run `make_glmmtmb_poisson_rs_crosslang_data.R` and
`run_glmmtmb_poisson_rs_reference.R` from `tools/stata-benchmarks` with the
development package loaded/installed. The checkpoint includes package-check
results and full source history in one flat ZIP. GitHub has not been updated.

The returned v2 log contains all six converged components; v4 contains the
balanced joint fit; v5 contains partial and disjoint joint fits. All three
joint models converged after one iteration using successful component starts.
Raw joint GSEM covariance differs from R, chiefly through the returned joint
model-based covariance. Independent audits and a reconstruction using only
Stata component curvature and joint score cross-products validate the R covariance.

No additional Stata run is required for this increment. See
`GLMMTMB-POISSON-RS-VALIDATION-20260923.md` for final evidence, correction
conventions, current verification, and limitations. NB2/logit random slopes
remain future increments.

## Completed R validation

On R 4.5.2 (Linux), glmmTMB 1.1.14, and marginaleffects 1.0.0:

- 48 focused random-slope expectations pass.
- All 278 glmmTMB regression expectations pass with no warnings.
- The complete suite reports 1,281 passing expectations, no failures or skips,
  and the same three existing GLM warnings in test output.
- The built source passes `R CMD check --no-manual` with zero errors, warnings,
  or notes, including vignette building and examples.
- Independent review identified the unit-sensitive boundary guard, confirmed
  its correction, and found no remaining substantive issue in the implementation
  or Stata specification.

These are the original implementation checks. The later Stata component
results are recorded in `GLMMTMB-POISSON-RS-BALANCED-20260923.md`; joint
cross-language validation and its numerical limits are now recorded in the
validation report. No GitHub platform-matrix run or push is included.

## Primary implementation references

- https://glmmtmb.github.io/glmmTMB/reference/get_cor.html
- https://glmmtmb.github.io/glmmTMB/articles/covstruct.html
- https://www.stata.com/manuals/memepoisson.pdf
- https://www.stata.com/manuals/memepoissonpostestimation.pdf
- https://www.stata.com/manuals/semexample38g.pdf
