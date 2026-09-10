# `glmmTMB` random-intercept logit design

Design date: 2026-09-10
Release version: `0.1.5`
Selected engine: `glmmTMB::glmmTMB()` 1.1.14 or later
Selected likelihood: unweighted binomial logit with one Gaussian random
intercept, fitted by Laplace approximation

## Scope

The initial adapter accepts one conditional binary outcome coded 0/1 or as a
two-level factor, a logit link, and exactly one grouping variable with a random
intercept. It requires ordinary unweighted maximum likelihood, no offset, and
a converged fit with a positive-definite Hessian and a finite variance estimate
away from the zero boundary.

Grouped binomial responses, zero-inflation and dispersion submodels, priors,
random slopes, crossed or deeper random-effects structures, and non-logit
families remain outside this increment.

## Parameter and likelihood contract

The SUEST parameter vector contains the conditional fixed effects followed by
`log_sigma`, the native `glmmTMB` theta parameter. Thus
`sigma = exp(log_sigma)`. This is also the scale directly comparable to
Stata's random-intercept log-standard-deviation parameter.

For grouping level `i`, the Laplace log-likelihood is

\[
\ell_i^L = \ell_i(\widehat b_i)
  - \frac{\widehat b_i^2}{2\sigma^2}
  - \log(\sigma) - \frac{1}{2}\log(H_i),
\]

where `b-hat_i` is the conditional mode and

\[
H_i = \sum_{j \in i} p_{ij}(1-p_{ij}) + \sigma^{-2}.
\]

`sandwich::estfun(model, full = TRUE)` evaluates the corresponding TMB
Laplace objective after zeroing observations outside each grouping level. It
returns additive group scores for the fixed effects and `log_sigma`. The
adapter places each complete group score on one retained observation and zero
elsewhere, so ordinary observation alignment continues to work without
pretending that integrated-likelihood observations are independent.

`vcov(model, full = TRUE)` supplies the full native covariance. The adapter
checks it against `model$sdr$cov.fixed` and stores `N` times that covariance as
the sandwich bread.

## Clustering and overlap

The random-intercept grouping variable is the default system cluster. A
supplied higher cluster must contain complete grouping levels. The finite-
sample correction is `G/(G-1)`, where `G` is the number of system clusters.
Models may have identical, partially overlapping, or disjoint samples. Because
the meat is formed from complete group or higher-cluster totals, disjoint
groups produce an exactly zero cross-model covariance block.

## Predictions and marginal effects

`type = "link"` is the fixed linear predictor `X beta`. `type = "response"`
is the population probability integrated over a new Gaussian random effect,

\[
E_u\{\operatorname{logit}^{-1}(X\beta + \sigma u)\},
\qquad u \sim N(0,1),
\]

computed with 20-point Gauss-Hermite quadrature. This agrees with direct
numerical integration and matches the estimand requested by Stata
`predict, mu marginal`. It is deliberately not an empirical-Bayes conditional
prediction. Because the response prediction depends on `log_sigma`, its
uncertainty enters `marginaleffects` delta-method calculations.

## Validation status

- The focused file contains 54 passing assertions and no warnings.
- An independent implementation reproduces the fitted Laplace likelihood to
  `1e-8`.
- Complete group scores sum to zero at the optimum and selected group scores
  agree with numerical likelihood derivatives to `2e-6`.
- The full native covariance and stored bread agree exactly at the configured
  tolerance.
- Default grouping, higher clustering, partial overlap, disjoint groups,
  integrated predictions, slopes, and coefficient replacement all pass.
- The targeted regression gate and final `R CMD check --no-manual` complete
  with status OK.
- The returned Stata 19.5 revision-3 log completed every component, adaptive
  `suest2`, and Laplace `gsem` case.
- Laplace component coefficients, covariance, and log likelihood agree within
  `1.22e-5`, `6.61e-7`, and `9.09e-8`. Marginal mean predictions and average
  slopes agree within `2.37e-6` and `3.19e-6`, despite Stata reporting that its
  prediction used the default seven quadrature points.
- Against the joint Laplace `gsem` systems, coefficients agree within
  `2.14e-5`; all covariance elements agree within `0.00113`, fixed-effect
  elements within `0.000895`, and covariance diagonals within `3.08%`.
  Stata computes numerical derivatives for `gsem` under Laplace and uses its
  own joint robust-covariance conventions, whereas the R route uses the
  independently validated TMB group scores and native component bread.
- The R and adaptive `suest2` disjoint cross-blocks are exactly zero. Joint
  Laplace `gsem` leaves numerical cross-block entries no larger than `2.10e-4`.

## Stata benchmark contract

The component comparison fits Stata `melogit, intmethod(laplace)` and checks
coefficients, native full covariance, log likelihood, marginal mean
predictions, and average slopes. `suest2` 1.0.0 intentionally refuses those
Laplace fits, so the benchmark records that refusal and then uses explicit
12-point adaptive Gaussian-Hermite component fits for the balanced,
partial-overlap, higher-cluster, and disjoint joint systems.

Revision 2 confirmed exact Laplace component agreement and successful adaptive
joint systems. Relative to the R Laplace systems, the adaptive Stata fits
differ by up to `0.00592` in fixed coefficients and `5.48%` on covariance
diagonals. Revision 3 added joint Laplace `gsem, vce(robust)` systems, following
Stata's documented method for reproducing SUEST with random effects. The
returned revision-3 results are encoded by
`build_glmmtmb_logit_ri_stata_fixture.R`; no further Stata run is required for
this increment.

The R adapter remains Laplace-based; it does not silently substitute adaptive
quadrature merely to satisfy the current `suest2` command.
