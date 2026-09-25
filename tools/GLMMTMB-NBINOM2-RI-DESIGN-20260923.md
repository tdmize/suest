# glmmTMB random-intercept NB2 development checkpoint

Date: 2026-09-23. Base release: 0.1.5 (`a46109f`). Development version:
0.1.5.9000. Engine: glmmTMB 1.1.14 or later. The returned revision-2 Stata gate
is complete; see `GLMMTMB-NBINOM2-RI-VALIDATION-20260923.md` for current evidence
and comparison limits.

The first returned Stata gate exposed a boundary fit missed by the original
validation. See `GLMMTMB-NBINOM2-RI-ROUND1-20260923.md` for the diagnosis,
new numerical boundary restriction, completed component comparisons, and
revision-2 odd/even benchmark. Validation counts below describe the initial
checkpoint; the validation report records the completed revision-2 comparison.

## Accepted model and parameter contract

Unweighted NB2 count outcomes with log link, one Gaussian random intercept,
Laplace maximum likelihood, and one estimated constant dispersion
(`dispformula = ~1`). The parameter order is conditional fixed effects,
`log_phi`, `log_sigma`. Conditional variance is `mu + mu^2/phi`; `phi` is NB2
size and `sigma` is the random-intercept SD. Stata's mean-dispersion parameter
is `alpha = 1/phi`, so `lnalpha = -log_phi`.

The adapter rejects dispersion regressors or offsets, zero inflation,
weights, conditional offsets, mapped parameters, priors, REML, random slopes,
multiple grouping variables, nonconvergence, and boundary random-effect
variance. Names `log_phi` and `log_sigma` are reserved. NB1 and other links
remain outside this increment. Systems contain NB2 models of the same type.

## Scores, bread, and alignment

The native full glmmTMB group scores retain all fixed effects and both nuisance
parameters. Full native covariance is checked against retained TMB covariance;
bread is `N * covariance`. The existing group alignment and finite-sample
correction apply: default clustering by random-intercept group, or higher
clusters containing whole groups. Partially overlapping observation samples
share group influence; fully disjoint groups have zero cross-model covariance.

An independent scalar Laplace check uses the conditional mode solving

`sum(phi*(y-mu)/(phi+mu)) - b/sigma^2 = 0`,

with `mu = exp(X beta + b)` and negative curvature

`sum(phi*mu*(phi+y)/(phi+mu)^2) + 1/sigma^2`.

The group log likelihood is the sum of NB2 log densities, minus
`b^2/(2*sigma^2) + log_sigma + log(curvature)/2`. Numerical derivatives of
this independent likelihood check the native score coordinates.

## Prediction contract

Response predictions are population means `exp(X beta + sigma^2/2)`;
link predictions are `X beta`. These integrate over a new Gaussian random
effect. The mean has zero direct derivative with respect to `log_phi`, but
dispersion remains in the joint covariance calculation. The derivative with
respect to `log_sigma` is `sigma^2 * mean`. The marginaleffects methods use
these predictions for means, slopes, discrete comparisons, and hypotheses.

## Initial validation (historical checkpoint)

The tests independently check likelihood, selected group score derivatives,
the full sandwich and higher-cluster aggregation, partial/disjoint samples,
analytic mean/slope/discrete-comparison gradients, cross-model contrasts,
parameter replacement, and unsupported specifications. Existing glmmTMB logit
and Poisson tests continue to cover shared adapter behavior.

The final glmmTMB regression run passes 150 assertions (44 NB2, 54 logit,
52 Poisson) without warnings. The 144-case numerical acceptance suite passes.
Independent code review identified and closed two validation loopholes:
permuted parameter maps and offsets with all-zero training values. Regression
tests demonstrated both failures before their fixes.

The built development package passes `R CMD check --no-manual` on R 4.5.2
(Linux) with 0 errors, 0 warnings, and 0 notes. Its full test suite reports
1,153 passing assertions, no failures or skips, and the same three pre-existing
GLM test warnings seen on v0.1.5. No new-platform CI run has been performed for
this unpublished development branch.

The flat Stata handoff provides fixed data, a do-file, and an R reference for
six component fits and three joint systems. Primary comparisons use Laplace
`menbreg` and `gsem`; optional adaptive `suest2` results are a distinct
approximation comparison. Component and joint nuisance parameters must be
transformed, including their covariance Jacobians, before comparison.

At this initial checkpoint, the Stata log had not yet been returned. It has
since been analyzed as recorded in the validation report above. Random slopes
remain a separate later increment. No new GitHub push is included in this work.

## Method references

- https://glmmtmb.github.io/glmmTMB/reference/nbinom2.html
- https://www.stata.com/manuals/memenbreg.pdf
- https://www.stata.com/manuals/memenbregpostestimation.pdf
- https://www.stata.com/manuals/semgsemfamily-and-linkoptions.pdf
