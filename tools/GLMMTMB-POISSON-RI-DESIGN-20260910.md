# `glmmTMB` random-intercept Poisson design

Design date: 2026-09-10
Release version: `0.1.5`
Selected engine: `glmmTMB::glmmTMB()` 1.1.14 or later
Selected likelihood: unweighted Poisson log with one Gaussian random
intercept, fitted by Laplace approximation

## Scope and parameter contract

The adapter accepts a nonnegative integer outcome, log link, and exactly one
grouping variable with one conditional random intercept. It excludes weights,
offsets, zero-inflation and dispersion submodels, priors, random slopes,
crossed or deeper random effects, negative binomial families, boundary
variance estimates, and nonconverged fits.

The parameter vector contains the conditional fixed effects followed by
`log_sigma`, the native TMB log random-intercept standard deviation. Full
grouping-level scores come from `sandwich::estfun(model, full = TRUE)`, and the
matching bread is `N * vcov(model, full = TRUE)`, checked against the retained
TMB covariance.

## Likelihood and clustering

For group `i`, with conditional mean `mu_ij = exp(X_ij beta + b_i)`, the
Laplace log likelihood is

\[
\ell_i^L = \sum_j \log f(y_{ij}; \mu_{ij}(\widehat b_i))
  - \frac{\widehat b_i^2}{2\sigma^2} - \log(\sigma)
  - \frac{1}{2}\log\left\{\sum_j \mu_{ij}(\widehat b_i)
  + \sigma^{-2}\right\}.
\]

The grouping variable is the default cluster. Higher clusters must contain
complete grouping levels. The same partial-overlap and disjoint-group
alignment contract used by the logit adapter applies.

## Predictions

The link prediction is `X beta`. The population response mean integrates over
a new Gaussian random effect and has the exact closed form

\[
E_u\{\exp(X\beta + \sigma u)\} = \exp(X\beta + \sigma^2/2).
\]

It is not an empirical-Bayes conditional prediction. Because it depends on
`log_sigma`, variance-parameter uncertainty enters `marginaleffects`
predictions and slopes.

## Validation status

- The focused file contains 52 passing assertions and no warnings, including
  the returned Stata fixture.
- An independent implementation reproduces the fitted Laplace likelihood to
  `1e-8`.
- A selected group score agrees with numerical likelihood derivatives to
  `2e-6`; complete group scores and the full native covariance are validated.
- Default and higher clustering, partial overlap, disjoint groups, integrated
  predictions, slopes, coefficient replacement, and deliberate refusals pass.
- Returned Stata 19.5 Laplace component fits agree within `5.41e-6` for
  coefficients, `1.31e-7` for native covariance elements, and `1.83e-7` for
  log likelihoods. Marginal means and average slopes agree within `1.07e-5`
  and `4.02e-6`.
- The balanced, partial/higher-cluster, and disjoint R joint systems agree with
  the matching Laplace `gsem` references within `7.35e-6` for coefficients and
  `0.00110` for covariance elements. The largest relative covariance-diagonal
  difference is `4.98%` in the smaller disjoint case; the balanced and partial
  cases are within `2.38%` and `1.98%`.
- `suest2` returns code 322 for Laplace `mepoisson` fits because it requires
  adaptive Gaussian-Hermite quadrature. Its 12-point adaptive systems are
  retained as a separate approximation comparison: coefficients differ from
  R's Laplace fits by at most `0.00673` and covariance elements by at most
  `9.40e-5`. Both R and adaptive `suest2` give an exactly zero disjoint
  cross-model covariance block.
