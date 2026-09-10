# Gamma random-effects panel Poisson design

Design date: 2026-09-10
Release version: `0.1.5`
Selected engine: `pglm::pglm()`
Selected likelihood: individual gamma random-effects Poisson with log link

## Contract

The panel-Poisson route admits only a numeric nonnegative integer outcome, one
individual grouping variable, a multiplicative gamma random effect, ordinary
unweighted maximum likelihood, `model = "random"`, `effect = "individual"`,
and `other = "sd"`. In this `pglm` parameterization, the final ancillary
parameter is the natural-scale gamma variance `alpha`. Offsets, weights, other
panel effects/models, nonconvergence, and a zero-boundary `alpha` are rejected.

The parameter vector contains every fixed effect and `alpha`. Analytic scores
are summed to one integrated-likelihood contribution per panel. The symmetrized
negative Hessian must reproduce the native covariance before the system is
returned. Panel ID is the default cluster; a supplied higher cluster must
contain complete panels.

Response predictions are `exp(X beta)`, the marginal expected count because
the gamma effect has mean one. Link predictions use `X beta`. The ancillary
variance therefore remains in the joint covariance but does not affect
predictions or marginal effects when coefficients are replaced.

## Validation status

- `pglm` analytic panel totals reproduce central numerical likelihood
  derivatives.
- Reconstructed information reproduces the native full covariance.
- Default panel clustering, higher clusters, partial overlap, disjoint panels,
  factor interactions, coefficient replacement, expected-count predictions,
  and average slopes pass 24 focused assertions.
- The deterministic Stata benchmark compares `pglm`'s `alpha` after conversion
  to Stata's `lnalpha` scale. It covers both native component fits and joint
  `suest2` systems with balanced, partially overlapping, higher-clustered, and
  disjoint samples.
- The returned Stata 19.5 benchmark completed every case. Component
  coefficients, native covariance, log likelihoods, mean predictions, and
  average slopes agree within `1.44e-7`, `8.05e-9`, `5.80e-7`, `1.35e-7`, and
  `4.71e-8`, respectively. Joint coefficients agree within `1.46e-6`, and
  disjoint cross-blocks are exactly zero.
- Full joint covariance equality is deliberately not required. `suest2` 1.0.0
  repeats each full `lnalpha` cluster score on every observation. On this
  six-period benchmark that inflates the ancillary meat by roughly (6^2),
  with smaller spillover into fixed-effect covariance through the native
  bread. Its requested higher-cluster route also reconstructs the integrated
  likelihood using the higher cluster instead of the model's panel. Directly
  emulating those Stata rules from the returned log reproduces the balanced
  Stata covariance within `4.13e-9`, locating the discrepancy outside the R
  implementation.
- The returned matrices are stored in a permanent regression fixture. R keeps
  one exact integrated-likelihood score per panel and aggregates those panel
  scores when a valid higher cluster is requested.
