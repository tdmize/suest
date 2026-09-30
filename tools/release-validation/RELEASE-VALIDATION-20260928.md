# suest release-validation record — 2026-09-28

**Goal:** establish model-specific evidence for widespread promotion of the R
package, including the standard errors of cross-model predictions and effects.

**Baseline:** GitHub main `29ae8519ed73fa4720b1138ef93e4aa8edc3d188`, version
0.1.6.9000. Published v0.1.6 is `a7a634a52db1e3de4f9eafa6a252753672c703f4`.
The evidence matrix describes the main baseline, not an assertion that a new
release candidate has passed. No remote push is authorized by this work.

**Approved scope:** the user accepted the prior recommendation to inventory
model evidence, fill substantive gaps, classify warnings/discrepancies, and
freeze/verify a release candidate. They specifically suggested official Stata
`suest`, Stata `gsem`, stacked datasets, and bootstrap comparisons.

## Plan and progress

- [x] Recover the exact source; inspect actual CI logs and numerical tests.
- [x] Inventory the supported routes and distinguish independent references
  from native-engine checks and finite-output checks (see evidence matrix).
- [x] Obtain a separate read-only source review of model admission and the
  less-common-family validation. The reviewer independently identified the
  convergence and penalized-multinomial guard gaps below.
- [x] Prepare minimal metadata corrections for the known check warning/note.
- [x] Prepare a first diagnostic gate using real invalid-fit examples,
  independent stacked GLMs, analytic contrast gradients, and a paired cluster
  bootstrap.
- [x] Execute gate 1 in R before and after the safeguards. All 40 stacked
  comparisons pass; the paired cluster bootstrap clears its diagnostic screen.
- [x] Add regression tests and correct confirmed model-admission defects.
  Three rejection expectations failed before the checks were added; all 14
  expectations in the new safeguard file pass afterward.
- [ ] Fill remaining family-specific full-covariance/contrast-SE gaps using
  the smallest suitable independent reference.
- [x] Run focused regression, existing acceptance suite, and local exact-source
  package checks on the current candidate; preserve output and package versions.
  Cross-platform CI awaits any later authorized push.
- [ ] Decide promotion readiness for each advertised route and publish only
  after the user's GitHub authorization.

## Execution environment and correction

The initial runtime search was too narrow: R was absent from the current PATH,
but R 4.5.2 and the installed package library remained in an earlier workspace.
After the user's correction, that runtime was located and restored. R checks
are executed by the assistant; the user only needs to execute Stata when a
specific additional reference is warranted. No new Stata run was performed.

The new runs use R 4.5.2 on Ubuntu 24.04, marginaleffects 1.0.0, sandwich 3.1-3,
nnet 7.3-20, pkgload 1.5.3, and testthat 3.3.2. Both gate logs record the actual
source fingerprints and session information. Historical CI evidence below
remains separate from these new executions.

The baseline's GitHub check run 36209160096 reports 1,630 passing expectations,
zero failures, and zero skips on Windows release, macOS release, Ubuntu release,
and Ubuntu devel. All four package checks have one warning (`maxLik` used in
tests but undeclared) and one note (nonstandard top-level files). Test-warning
counts are 3, 4, 3, and 21 respectively, mostly the existing separation examples,
one macOS unsupported-NB1 rejection fixture, and upstream R-devel deprecations.
Workflow success therefore does not mean a warning-free `R CMD check`.

## Findings requiring action

### RV-01: nonconverged ordinary GLMs can reach inference

Source evidence: `R/adapters.R`'s ordinary GLM branch does not inspect
`model$converged`. The subsequent `R/suest.R` checks inspect missing parameters,
score dimensions, and parameter ordering but do not establish convergence.
The baseline test `model-specific covariance blocks are preserved` itself
constructs a nonconverged GLM in the CI log and still passes.

Impact: finite coefficients/covariance can be mistaken for valid inference
from a fitted optimum. The first gate creates an actual one-iteration fit and
records whether it is accepted, with a converged control. It does not falsify a
model object's convergence flag. Complete/quasi-separation requires a separate
follow-up; a `converged=TRUE` flag alone is not a proof of a finite MLE.

Status: **reproduced and fixed locally**. The baseline accepted the actual
one-iteration fit with `converged=FALSE` and returned finite covariance. The
new check rejects it before constructing inference. Separate `glm` and `glm2`
regressions fail before the fix and pass afterward; converged controls pass.
The original core success fixture was replaced with a stable two-group example
while retaining its covariance assertions and adding convergence assertions.
No covariance formula was changed.

### RV-02: penalized multinomial fits enter ordinary-ML equations

Source evidence: the `multinom` branch in `R/adapters.R` accepts the class
without examining `decay`. `.suest_multinom_components` in `R/internals.R`
constructs the ordinary unpenalized scores and information, omitting penalty
terms. The existing adjusted-GLM safeguards do not cover this path.

Impact: a fit obtained with `nnet::multinom(decay > 0)` does not solve the
estimating equations used for its reported SUEST covariance. Rejecting this
unsupported estimator is the intended smallest correction. Gate 1 fits real
`decay=0` and `decay=2` controls, computes the unpenalized score sum independently,
and records acceptance/rejection. Do not claim a confirmed runtime fix until
the rejection regression is observed failing and then passing.

Status: **reproduced and fixed locally**. The baseline accepted `decay=2`
despite a maximum absolute ordinary score sum of 5.05289495. Nonzero decay now
produces an explicit unsupported-fit error. The regression fails before the
fix and passes afterward, while the `decay=0` control remains supported.

### RV-03: incomplete independent certification for extended families

Native single-model block preservation is valuable, but its reference often
calls the same `sandwich::estfun()` and `bread()` used by production. Finiteness
of a marginal-effect SE does not establish its magnitude. Generic/fractional
GLMs, beta, survival, truncated, heteroskedastic, and zero-inflated families
need the specific missing full-system/contrast checks listed in the matrix.

Gate 1 covers ten GLM family/link routes under four sample patterns. A separate
stacked fit supplies the complete covariance and independent analytic gradients
supply prediction-difference and finite-change-difference SEs. This independently
checks system assembly and postestimation; it still shares the native GLM and
sandwich implementation, so it is not an independent likelihood-engine audit.

Status: **gate 1 passes; remaining extended-family gaps stay open**. See the
numerical results below. These GLM checks do not certify beta, survival,
truncated, heteroskedastic, or zero-inflated routes.

### RV-04: estimator and specification boundaries need explicit audit

The generic GLM branch admits any family with one of six link names; that broad
admission is not numerical certification of arbitrary custom families/classes.
Beta bias-corrected/bias-reduced fits and penalized/robust survreg variants need
native-object probes before a support claim or a code change. These are
**candidates for investigation**, not demonstrated new numerical defects.

### RV-05: check metadata and warning classification

The local branch adds `maxLik` to Suggests and excludes `DEVELOPMENT.md`,
`UPDATE-INSTRUCTIONS.txt`, and `_pkgdown.yml` from the built package. The
complete candidate patch includes these four additions, the safeguards,
their regressions, improved success
fixtures, and user documentation. `R CMD build` succeeds with vignette creation.
The final package-check result is recorded in the execution summary.

## Gate 1 numerical results

The baseline run completed all calculations but exposed a cleanup error in
the diagnostic script's message-sink restoration, giving process status 1
after writing results. The cleanup was fixed; the final complete run exits 0.
This was a test-harness error, separate from the package safeguards.

The final gate has 42 PASS results (40 stacked comparisons and two valid-fit
controls), two PASS_REJECTION results, and one SCREEN_CLEAR bootstrap result.
The intentionally nonconverged input emits its expected native GLM warning;
no stacked comparison or bootstrap warning occurs.

| Quantity across the 40 stacked checks | Largest difference | Preset bound |
|---|---:|---:|
| Coefficient, absolute | 1.171448e-7 | 2e-6 |
| Joint covariance, absolute | 3.298793e-9 | Reported diagnostic |
| Joint covariance, standardized by marginal SE products | 3.732340e-7 | 1e-5 |
| Prediction-difference estimate, absolute | 6.234058e-9 | 2e-7 |
| Finite-change-difference estimate, absolute | 1.490822e-7 | 2e-7 |
| Prediction-difference SE, relative | 4.701311e-8 | 2e-4 |
| Finite-change-difference SE, relative | 1.083200e-7 | 2e-4 |

The paired bootstrap completes all 999 draws without failed fits. Its SE is
0.02841110409 versus analytic 0.02843940897, a ratio of 0.99900473 (about 0.10%
lower). Its approximate relative Monte Carlo SE is 2.24%; this is supportive
evidence for the tested conditional contrast, not a coverage guarantee.

## Execution summary

- New safeguard regression file: 3 expected failures before the fix; 14
  expectations pass afterward.
- Focused core, model-adapter, and safeguard files: 85 expectations pass,
  no warnings or failures.
- Full source suite: 1,649 expectations pass across 210 test blocks; zero
  failures, errors, warnings, or skips. All suggested packages are installed.
- Full numerical acceptance suite: 144 pass, none fail, including the paper
  replications. Previously cached public replication datasets were reused.
- The first package check identified a nonconverged logistic example in the
  `suest_newdata` documentation. It was changed to a stable linear example
  demonstrating the same separate-sample averaging workflow; the corrected
  example passes when executed directly. This initial check is not a clean
  release result.
- Final `R CMD build` succeeds, including vignette creation. Final
  `R CMD check --no-manual` exits 0 with **Status: OK**, no errors, warnings,
  or notes. Installed-package tests report 1,649 passing expectations, zero
  failures, warnings, or skips; examples and vignette rebuilding pass.
- A separate read-only source review found no material regression or scope
  problem in the safeguards, replacement fixtures, or gate.

These are local Ubuntu/R 4.5.2 results. The changed candidate has not been
pushed or rerun on Windows, macOS, or R-devel. The older CI results above
must not be represented as cross-platform verification of these changes.

## Gate 1 specification

Run the supplied `suest_release_gate1.R` from the existing package source root.
It uses `pkgload` to load local source and records source fingerprints, git HEAD,
working-tree state, dependency versions, all warnings/errors, and each result.
It does not install packages, apply patches, or alter package code.

The 40 stacked cases are Gaussian identity, binary logit/probit/cloglog,
Poisson log, Gamma log, and fractional logit/probit/cloglog/loglog, each with
identical samples, partial overlap, disjoint observations sharing clusters,
and disjoint clusters. Explicit clustering aligns the ordinary-model
`G/(G-1)` correction; the stacked reference uses HC0 with cluster adjustment.
Ordinary binary probit uses independently derived observed information from
`log Phi((2*y-1)*eta)`. Its expected-information result is logged separately;
fractional probit retains the native quasi-GLM expected-information convention.
This deliberately does not imply validation of every unclustered correction.

Predeclared numerical bounds: maximum coefficient gap < 2e-6; maximum covariance
gap standardized by reference marginal SE products < 1e-5; cross-model estimand
gap < 2e-7; relative contrast-SE gap < 2e-4. Failed bounds are recorded; they
must not be widened merely to obtain a pass.

The supplementary 999-replicate logit bootstrap resamples whole clusters,
jointly refits both models on the same draws, and evaluates both finite-change
effects on the original fixed covariate distribution. Failed fits are counted
and reported; no replacement draws hide failures. The bootstrap SE ratio is a
diagnostic, with an approximate Monte Carlo SE, not a deterministic equality
test or coverage guarantee. A `SCREEN_CLEAR` result does not certify the package.

Status vocabulary: `PASS` = the named deterministic bounds passed;
`PASS_REJECTION` = a targeted rejection occurred; `ACCEPTED_REQUIRES_REVIEW` = an
invalid-estimator probe was accepted; `FAIL` = a deterministic discrepancy;
`ERROR` = incomplete case; `REVIEW`/`REVIEW_WARNING`/`REVIEW_REJECTION` = inspect
the logged reason; `SCREEN_CLEAR` = bootstrap screening criterion met.
`GATE1_COMPLETED=1` means the diagnostic run finished, not that all cases passed.

## Choice of subsequent reference

| Situation | First independent reference | Interpretation requirement |
|---|---|---|
| Ordinary ML route supported by official Stata suest | Official `suest` | Match sample IDs, weights, nuisance coordinates, and finite-sample correction; identify official suest separately from suest2. |
| A distinct joint or stacked likelihood can represent both models | `gsem` or stacked fit | All equation-specific coefficients and nuisances must remain free; cross-equation latent covariance restrictions must reproduce separate fitting. |
| Disagreement caused by curvature or quadrature | Independent likelihood derivatives / matched approximation | Check point estimates, scores, information, then sandwich, then prediction gradients. |
| No convenient deterministic external route | Paired observation or cluster bootstrap | Refit both equations within every shared resample; hold evaluation covariates fixed for the current conditional estimand; quantify Monte Carlo error and failures. |
| Survey design | Existing design linearization or design-valid resampling | Never substitute an ordinary case bootstrap for stratification/PSU/FPC handling. |

Expected-versus-observed information must be matched separately for ordinary
binary probit (observed) and survey/fractional probit (native expected).
Stata-versus-R differences due to these choices, quadrature, or small-sample
corrections remain explicit
conventions. Raw GSEM discrepancies and the failed disjoint-logit reconstruction
remain open and separate from independent R passes. No additional broad Stata
run is requested until the first R results identify which comparison is needed.

## Next decision

RV-01 and RV-02 are fixed and their first validation gate is executed. Next,
prioritize the missing extended-family cross-covariance and contrast-SE
references. Continue R execution here; request only targeted Stata runs from
the user when needed. Do not promote all routes as equally certified, and do
not start new model support or MI functionality in this validation pass.

## Review decisions

- The first stacked-probit reference used native expected information. Fresh
  source review identified the mismatch with the ordinary-probit contract;
  it was replaced by independent observed-information algebra, retaining the
  expected-information result as a diagnostic. No production formula changed.
- Setup errors now print to the open log and return `GATE1_COMPLETED=0`.
  A failure before the log can be opened still requires the console message.
- The new admission regressions were run against unchanged production code
  and failed as intended, then passed after the minimal checks were added.
  No covariance or prediction formula was changed.
- A separate read-only review of the final safeguards, test fixtures, and
  diagnostic script identified no material regression or scope problem.
