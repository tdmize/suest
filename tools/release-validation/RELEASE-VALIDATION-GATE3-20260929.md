# Panel, integrated-model, and survey contrast validation — 2026-09-29

Baseline 29ae851, suest 0.1.6.9000, with the local gate-1 and gate-2 corrections.
No GitHub push or new Stata run. R execution was performed by the assistant.

## Result and release implication

**100 of 100 cases pass the numerical gates**, covering 300 A-minus-B
contrasts: average predictions, finite changes, and average x slopes.
Ninety-six cases have no warnings. Four gamma–Poisson cases retain explicitly
reviewed native fitting warnings; they are not described as warning-free.
No new score, covariance, prediction, or gradient formula defect was found.

A supplemental boundary probe confirmed an important limitation for within/FE
average prediction levels at the training means: their covariance is
structurally degenerate under the current conditional-mean convention, and
numerical SEs can be missing or nearly zero. Those average-level confidence
intervals and tests are **unsupported**, not evidence of perfect precision.
The limitation is now explicit in README and help. Slope and finite-change
contrasts pass. No estimator or covariance formula was changed to hide it.

Confidence is stronger for the declared panel/GEE/survey contrast targets.
The remaining release work is a final targeted sweep of alternative-engine
and bivariate-probit contrast gaps, then current supported-platform checks of
the frozen candidate. Existing random-slope GSEM discrepancies remain open
convention/precision qualifications, not erased by these results.

## Coverage and independent references

| Route | Specifications × patterns | New evidence |
|---|---:|---|
| plm FE, BE, Swamy-Arora RE | 3 × 4 = 12 | Independent transformed-equation covariance plus analytic contrast gradients |
| Gaussian nlme ML | 1 × 4 = 4 | Independent mean/slope gradients, conditional on audited joint V |
| GEE | 5 means × 2 correlations × 4 = 40 | Independent mean-equation covariance plus analytic gradients |
| pglm RE | Logit/probit/gamma–Poisson × 4 = 12 | Independent integrated means/gradients, conditional on audited joint V |
| glmmTMB RI | Logit/Poisson × 4 = 8 | Independent integrated means/gradients, conditional on audited joint V |
| Survey | Gaussian/logit/probit × with/without FPC × 4 = 24 | Independent expected-information influence, full-design PSU covariance, analytic gradients |
| Total | 100 | 76 fresh covariance reconstructions; 24 conditional-propagation checks |

GEE means: Gaussian identity, binary logit/probit/cloglog, Poisson log;
working correlations: independence and exchangeable. Working scale/correlation
and RE quasi-demeaning parameters come from the native fitted engine. The
reference reconstructs their declared estimating equations independently;
it does not independently re-estimate those working parameters.

The four patterns are common observations; unequal, unbalanced partial
samples; disjoint panels sharing higher clusters; and disjoint higher
clusters. There are initially 160 panels × 5 visits, in 40 higher clusters.
The unequal samples retain 747 and 728 observations under fixed visit/panel
selection rules. Distinct outcomes are used for the two models.

The evaluation population is fixed: observed covariates shifted by x+0.25
and z+0.15. Survey averaging uses the fixed original weights; other routes use
equal weights. Finite changes set x from −0.5 to +0.5. Slopes are derivatives
with respect to x. No inference about estimating the averaging distribution
is claimed. All nuisance coordinates remain in V and the gradients.

Normal-mixture binary predictions use independently generated Golub-Welsch
quadrature matching the declared 12-point pglm or 20-point glmmTMB rule.
Separate adaptive integration for logit and the closed-form integrated probit
formula check the continuous target. Thus matching quadrature is not silently
presented as exact integration. Gamma–Poisson expected counts do not depend
on alpha; Gaussian-RE Poisson counts do depend on the random-effect variance.

The 24 conditional-propagation cases use the full suest covariance already
covered by earlier likelihood/external audits. They newly test independent
mean functions, gradients, and uncertainty propagation; they cannot detect
a covariance defect shared with those earlier audits. No suest prediction,
coefficient setter, or internal gradient supplies any reference target.

## Bounds and numerical results

No threshold was loosened to obtain passing comparisons.

| Quantity | Preset bound | Largest observed gap |
|---|---:|---:|
| Full covariance / marginal SE products (76 cases) | 2e-6 | 1.371e-7 |
| Analytic versus numerical target gradient, scaled | 2e-5 | 1.039e-9 |
| Contrast estimates, absolute | 2e-6 | 8.587e-12 |
| Prediction-contrast relative SE | 3e-4 | 1.730e-8 |
| Finite-change-contrast relative SE | 3e-4 | 5.029e-8 |
| Slope-contrast relative SE | 3e-4 | 1.567e-7 |
| Quadrature versus continuous estimate, absolute | 1e-5 | 7.416e-10 |
| Quadrature versus continuous relative SE | 1e-3 | 5.715e-9 |

These comparisons are sensitive to the terms they are meant to test.
Omitting random-effect variance derivatives changes SEs by as much as 1.87%
for pglm logit, 1.34% for pglm probit, 1.73% for glmmTMB logit, and 10.64%
for glmmTMB Poisson. Alpha and Gaussian-ML scale gradients are exactly zero
for their mean targets. Omitting cross-model covariance changes some SEs by
more than 100%. Ordinary disjoint-cluster cross-blocks are exactly zero.
All six corresponding survey cross-blocks remain nonzero because the full
stratum-centered design induces cross-domain covariance; zeroing them would
be wrong. The metrics file records exact values and case names.

## Native gamma–Poisson warnings: resolved and retained

pglm 0.2-4's BFGS optimizer visits negative trial values of its natural gamma
variance alpha. Its internal likelihood uses shape=1/alpha, causing log and
Hessian square-root warnings at those invalid trial points. Every retained
warning was traced while its native likelihood frame was live; every trial
alpha was negative, and every message was `NaNs produced` during fitting.
No prediction/reference warning was silently reclassified.

All eight final component fits have native convergence code zero. Independent
exact integrated gamma–Poisson log likelihoods agree within 1.14e-13;
information-scaled score sums are at most 9.78e-7; independent Richardson
information stability is at most 4.72e-8. Information is positive definite,
and re-evaluating the native objective/derivatives at the final coefficients
is finite and warning-free. Additional bounds (chosen before evaluating
these diagnostics) were 2e-6, 2e-3, and 5e-5, respectively.

The four statuses are `PASS_REVIEWED_FITTING_WARNING`. The original warnings
and trial values remain in the logs/RDS. This does not certify arbitrary
poorly converged fits or justify suppressing unexplained warnings.

## FE average-level boundary: explicit limitation

The reconstructed FE intercept is ybar − xbar beta. Its covariance propagates
slope uncertainty through the mapping M = rbind(−xbar, I). At an evaluation
mean xstar, the implied prediction variance is

    (xstar − xbar)' V_beta (xstar − xbar).

At xstar=xbar it is zero. A separate dummy-variable OLS calculation confirms
that nested-cluster scores annihilate this training-average direction:
its uncorrected SE is 3.67e-17, while the x-slope SE is 0.04486. Native
marginaleffects produced one missing average-prediction SE, another about
7.88e-12, and a contrast SE about 4.25e-12 with an enormous test statistic.
These are boundary/numerical manifestations, not valid average-level tests.
They do not establish zero sampling uncertainty in the observed response
mean or a population mean across fixed effects.

The initial exploratory gate encountered this zero-variance target before
performing the relative-SE comparison. The final gate uses the same fixed
shifted population for every route, and the original behavior is retained
in the supplemental probe. The final gate validates conditional slope-based
prediction propagation; it does not close the unsupported level-inference
case. README and the package help now state this boundary explicitly.

For background on the normalized FE intercept (not a claim about this
specific robust SE boundary), see StataCorp's primary explanation:
https://www.stata.com/support/faqs/statistics/intercept-in-fixed-effects-model/

## Review, execution history, and package verification

Independent read-only review found no blocking reference-formula issues.
It verified the analytic gradient and covariance conventions, sensitivity
metrics, disjoint-cluster geometry, and the conditional-propagation claims.
The subsequent warning tracing passed its own exact-likelihood checks.
The reviewer also confirmed the FE boundary algebra and scope interpretation.

The pilot initially omitted pglm Poisson's required `other="sd"`; the package
correctly rejected it. This harness error was corrected. An exploratory run
used balanced partial samples and unshifted covariates; final sampling was
made genuinely unbalanced and the evaluation population made nondegenerate.
These choices precede the final numerical conclusions. Pilot/exploratory
results are preserved and are not counted as clean final evidence.

The combined final RDS retains 96 cases from the final full run and replaces
the four warning-review entries with their focused replays. Their three
numerical contrast metric vectors were verified equal within 1e-12; only
warning diagnostics/classification were added. A fresh full default run of
the saved script reproduces the same 100-case design and warning checks.

Only README, NEWS, matching help/roxygen comments, and .Rbuildignore changed in this stage;
production executable R and package tests are unchanged from gate 2.
Gate 2's source/installed suites passed 1,667 assertions with no failures,
warnings, or skips; its acceptance suite passed 144 checks and its full
R CMD check was clean. Those are preserved prior results, not new suite runs.
The stage-3 documentation build/check result is recorded separately, with
package tests intentionally not rerun: `R CMD check --no-tests --no-manual`
finished with **Status: OK**, no errors, warnings, or notes. The initial check
reported the internal `.superpowers` folder; adding its build exclusion resolved
the note. Examples and vignette checks passed. The combined patch was applied
to a clean baseline and verified byte-for-byte against the candidate.

## Reproduction and remaining scope

From the candidate package source root:

    source("tools/release-validation/suest_release_gate3.R")

For a focused run:

    options(suest.gate3.define_only = TRUE)
    source("tools/release-validation/suest_release_gate3.R")
    run_suest_release_gate3(output="path/to/results", selected="pglm_log")

`GATE3_COMPLETED=1` means execution finished; inspect statuses and warnings.
The script does not equate a completed run with a passed gate.

This stage does not newly certify mixed-family systems, arbitrary links,
factor/offset designs, all integration variances, random-slope routes, or
population-level FE mean inference. Earlier evidence covers some separately;
this report does not promote it into a new unrestricted claim.
