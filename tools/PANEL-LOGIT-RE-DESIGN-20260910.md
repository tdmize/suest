# Random-effects panel logit design — first nonlinear-panel increment

Design date: 2026-09-10
Release version: `0.1.5`
Selected engine: `pglm::pglm()`
Selected likelihood: individual random-intercept binary logit with 12-point
nonadaptive Gauss–Hermite quadrature

## Candidate audit

The design audit used R 4.5.2, `marginaleffects` 1.0.0, `lme4` 2.0-6,
`glmmTMB` 1.1.14, and `pglm` 0.2-4.

| Engine | Useful evidence | Blocking issue for the first adapter |
|---|---|---|
| `lme4::glmer()` | Native `marginaleffects` support; recoverable frame and grouping; selectable Laplace or adaptive quadrature | Ordinary `vcov()` contains fixed effects only; no direct panel-score matrix for the fixed effects plus random-effect SD |
| `glmmTMB::glmmTMB()` | Native `marginaleffects` support; recoverable frame; TMB objective and gradient internals | Initial binary route is Laplace-only; the ordinary coefficient view omits the variance parameter unless the full interface is requested |
| `pglm::pglm()` | Dedicated panel likelihood; retained frame/index/design; fixed effects plus natural-scale `sigma`; analytic score decomposition and full Hessian/covariance | The raw fit is not a native `marginaleffects` model and uses nonadaptive quadrature; the combined `suest_model` therefore supplies prediction and coefficient-replacement methods |

On the deterministic audit data, `pglm(..., R = 12)` and
`lme4::glmer(..., nAGQ = 12)` agreed to the displayed precision in log
likelihood and random-intercept SD, with fixed-effect differences below
`2e-5`. `pglm` was selected because it exposes the complete likelihood
parameter system and independently checkable analytic derivatives.

That selection was correct for the first `xtlogit, re` increment, but a later
audit of `glmmTMB` 1.1.14 found complete interfaces that were not used in the
initial comparison: `sandwich::estfun(model, full = TRUE)` returns one
Laplace-likelihood score per grouping level, including the random-effect
parameter, and `vcov(model, full = TRUE)` returns the matching covariance.
Those interfaces support a separate `melogit`-style adapter documented in
`GLMMTMB-LOGIT-RI-DESIGN-20260910.md`; they do not change this adapter's
nonadaptive-quadrature contract.

## Initial mathematical contract

1. **Model.** Binary logit, one individual grouping variable, one Gaussian
   random intercept, ordinary unweighted maximum likelihood, and no offset.
2. **Integration.** The fit must use `R = 12`. `pglm` uses nonadaptive
   Gauss–Hermite quadrature. The Stata benchmark must therefore use
   `intmethod(ghermite) intpoints(12)`, not Stata's default adaptive method.
3. **Parameters.** The SUEST vector contains all fixed effects and the
   natural-scale random-intercept standard deviation named `sigma`. `pglm`
   can return the equivalent negative root because its likelihood is symmetric
   in `sigma`; the adapter canonicalizes the estimate to `abs(sigma)` and
   transforms the score and covariance consistently.
4. **Scores.** `pglm$gradientObs` is an additive within-panel decomposition of
   the integrated likelihood score, not a collection of independent
   observation scores. The adapter sums it within panels and stores exactly
   one nonzero score row per panel. The system meat is then formed from panel
   totals.
5. **Bread.** The symmetrized negative fitted Hessian is inverted and must
   reproduce the transformed native `pglm` covariance within `1e-8` before a
   system is returned. The stored bread is `N` times this covariance, matching
   the package's sandwich convention.
6. **Clustering.** Panel ID is the default system cluster. The finite-sample
   correction is `G/(G-1)`, consistent with the existing integrated Gaussian
   panel and GEE routes. A supplied higher cluster must contain whole panels.
7. **Overlap.** Observation IDs align the retained rows, while cross-model
   covariance is formed from complete panel or higher-cluster score totals.
   Models may use identical, partially overlapping, or disjoint observation
   and panel samples.
8. **Predictions.** `type = "response"` is the probability integrated over
   the Gaussian random effect using the fitted 12-point rule. This matches the
   default `margins` estimand after Stata `xtlogit, re`. `type = "link"` is the
   fixed part `X beta`, equivalent to a zero random effect. Integrated response
   predictions depend on `sigma`, so its uncertainty enters the delta method.
9. **Marginal effects.** The original `pglm` object is not supported directly
   by `marginaleffects`; the combined `suest_model` is supported through the
   package's existing custom methods. Tests cover predictions, slopes, factor
   interactions, coefficient replacement, and standard errors.
10. **Refusals.** Probit, pooling/within/between models, other quadrature-point
    counts, time or two-way effects, weights, offsets, zero-boundary variance,
    and nonconverged fits remain outside this increment.

## Validation status

- Focused nonlinear-panel file: 25 passed assertions, 0 failures.
- Targeted regression gate: passed for core, GEE, marginaleffects, MI, linear
  panel, Gaussian panel ML, survey, and the new panel-logit route; two expected
  GLM fitting warnings from the pre-existing core stress test.
- The returned revision-2 Stata 19.5 log completed every benchmark case.
- The nonadaptive component fits agree on coefficients within `6.95e-8`, log
  likelihoods within `8.31e-8`, integrated mean predictions within `3.49e-9`,
  and average slopes within `3.39e-9`.
- The adaptive Stata `suest2` systems differ from the nonadaptive R systems by
  at most `1.06e-3` in coefficients and `3.10e-3` across all covariance
  elements. Fixed-effect covariance elements differ by at most `1.17e-4`; the
  maximum is an ancillary variance element in the disjoint case. Balanced and
  partial-system cross-block differences are below `1.77e-5`; disjoint
  cross-blocks are exactly zero in both implementations.
- The returned values are stored in a permanent regression fixture with
  explicit cross-integration tolerances. No further Stata run is required for
  this increment.

## Stata benchmark revision note

Benchmark revision 1 confirmed exact single-model agreement under 12-point
nonadaptive quadrature, including coefficients, log likelihoods, integrated
mean predictions, and average marginal effects. `suest2` then stopped with
`r(322)` because its random-effects panel-logit adapter requires adaptive
Gaussian-Hermite quadrature. Revision 2 therefore keeps the nonadaptive fits
for the exact pglm likelihood comparison and uses explicit 12-point adaptive
fits for the joint `suest2` covariance cases. The returned log must be used to
quantify the adaptive/nonadaptive difference before cross-language covariance
acceptance is certified.
