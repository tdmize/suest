# Post-0.1.4 development roadmap

This roadmap records the priority order selected after certification of
`suest` 0.1.4. Each model family remains a separate, benchmarked increment;
support is not inferred from a fitting class alone.

## Priority 1: nonlinear panel models

Add nonlinear panel families one at a time. Begin with a narrow binary-response
panel route whose likelihood scores, information convention, panel clustering,
sample alignment, native covariance blocks, predictions, and cross-model
covariance can all be independently verified. Choose the first fitting engine
only after a reproducible R/Stata benchmark establishes its parameterization
and ancillary-parameter contract.

Candidate sequence:

1. random-effects panel logit;
2. random-effects panel probit;
3. conditional fixed-effects logit;
4. panel count families after the binary routes are stable.

Do not combine these into a generic nonlinear-panel adapter before each family
has its own score, bread, covariance, overlap, and marginal-effects gates.

## Priority 2: general GLMM and multilevel models

Begin with a narrow two-level random-intercept GLMM family, then add random
slopes, additional families/links, and deeper nesting only after the first
route is numerically certified. The design phase must explicitly decide
whether ancillary covariance parameters participate in the joint system and
how likelihood approximations and integration settings affect the score and
bread.

Candidate sequence:

1. two-level random-intercept binary GLMM;
2. two-level random-intercept count GLMM;
3. random slopes and correlated random effects;
4. deeper nesting or crossed effects.

Engine candidates should be evaluated separately rather than treated as
interchangeable. Native covariance reproduction and stable observation- or
cluster-level score reconstruction are release gates.

## Lower priority

- MI estimand-level pooling for predictions, slopes, comparisons, and
  contrasts;
- `ivtobit` and Heckman selection when independently testable estimators are
  available;
- generalized ordered models after verified scores;
- remaining survival parameterizations;
- correlated-random-effects panel benchmarking;
- additional GEE correlations and count families;
- broader survey designs or inference extensions.

The one-stage survey common-design, coefficient-only, explicit-ID contract in
0.1.4 remains unchanged while these separate increments are developed.
