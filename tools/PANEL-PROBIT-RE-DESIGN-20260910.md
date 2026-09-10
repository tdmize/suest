# Random-effects panel probit design

Design date: 2026-09-10
Release version: `0.1.5`
Selected engine: `pglm::pglm()`
Selected likelihood: individual random-intercept binary probit with 12-point
nonadaptive Gauss-Hermite quadrature

## Contract

The panel-probit route deliberately reuses the certified panel-logit
architecture. It admits only a numeric 0/1 or two-level-factor outcome, one
individual grouping variable, one Gaussian random intercept, ordinary
unweighted maximum likelihood, `model = "random"`, `effect = "individual"`,
and exactly `R = 12`. Offsets, weights, other panel effects/models,
nonconvergence, and a zero-boundary random-intercept SD are rejected.

The parameter vector contains every fixed effect and the natural-scale
random-intercept SD `sigma`. Scores are summed to one integrated-likelihood
contribution per panel. The symmetrized negative Hessian must reproduce the
canonicalized native covariance before the system is returned. Panel ID is the
default cluster; a supplied higher cluster must contain complete panels.

Response predictions integrate `pnorm(X beta + u)` over the retained 12-point
quadrature rule. Link predictions use `X beta`, corresponding to a zero random
effect. Because integrated response predictions depend on `sigma`, its
uncertainty enters `marginaleffects` delta-method standard errors.

## Validation status

- `pglm` analytic panel totals reproduce central numerical likelihood
  derivatives.
- Reconstructed information reproduces the native full covariance.
- Default panel clustering, higher clusters, partial overlap, disjoint panels,
  factor interactions, coefficient replacement, integrated predictions, and
  average slopes pass 20 focused assertions.
- The returned deterministic Stata 19.5 benchmark completed every case.
  Nonadaptive `xtprobit, re` component coefficients agree within `9.29e-7`,
  native covariance within `6.89e-8`, log likelihoods within `1.99e-7`, mean
  predictions within `2.99e-7`, and average slopes within `3.81e-8`.
- The adaptive Stata `suest2` systems differ from the nonadaptive R systems by
  at most `0.00113` in coefficients and `0.00117` across all covariance
  elements. Fixed-effect covariance elements differ by at most `1.33e-4` and
  covariance diagonals by at most `0.394%`. Disjoint cross-blocks are exactly
  zero in both implementations.
- The returned matrices are stored in a permanent regression fixture with
  explicit cross-integration tolerances. No further Stata run is required for
  this increment.
