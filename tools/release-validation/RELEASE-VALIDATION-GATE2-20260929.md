# Extended-family release validation — 2026-09-29

Baseline: `29ae8519ed73fa4720b1138ef93e4aa8edc3d188`, suest 0.1.6.9000.
This candidate includes the gate-1 safeguards and the gate-2 changes below.
Nothing has been pushed. No additional Stata run is needed to resolve the
findings in this stage. All R execution was performed by the assistant.

## Decision at this checkpoint

The extended-family gate passes **108 of 108 cases**, with no warnings.
This stage found and corrected a real interval-censoring covariance defect
that agreement with the native sandwich routine had missed. These results
increase confidence in the tested routes; they do not yet certify every
supported model, specification, dependency version, or operating system.
Panel/integrated-model and survey nonlinear contrast gaps in the evidence
matrix remain, as does testing the frozen candidate across supported platforms.
Previously documented random-slope GSEM differences remain unchanged.

## What was checked independently

Twenty-seven specifications each use four observation/cluster patterns:
identical availability, partial overlap, disjoint observations sharing clusters,
and disjoint clusters. Outcomes differ between the two components. Truncated
models additionally select on their own outcome, so even the first pattern can
have unequal realized samples. Each dataset has 1,600 initial observations
in 400 clusters; sample selection is explicit.

| Route | Specifications | Cases |
|---|---|---:|
| Beta ML | Four mean links with modeled log precision; constant identity precision; logit with observed information | 24 |
| Gaussian survreg | Left, right, interval censoring; fixed-scale left censoring | 16 |
| Parametric AFT | Weibull, lognormal, loglogistic; right and interval censoring | 24 |
| censReg Gaussian | Left, right, two-limit censoring | 12 |
| Truncated Gaussian | Left and right truncation | 8 |
| Heteroskedastic binary | Probit and logit | 8 |
| Zero-inflated counts | Poisson and negative binomial, each with logit/probit inflation | 16 |
| Total | 27 specifications × 4 patterns | 108 |

The reference writes each observation log likelihood using distribution
primitives, differentiates it independently, aggregates the complete score by
cluster ID, and combines it with independently reconstructed information.
Beta's default uses expected information derived from trigamma formulas;
`hessian = TRUE` uses observed information. All other routes use observed
likelihood curvature. Scale, precision, inflation, and dispersion parameters
remain in the complete covariance, including cross-model blocks. The common
finite-sample correction is G/(G−1).

Native fitted parameters locate the evaluation point. No suest or sandwich
score, bread, covariance, coefficient setter, or prediction helper supplies
the reference. Independent prediction gradients produce A-minus-B average
prediction and finite-change contrast SEs, compared with marginaleffects.
Predictions average over fixed evaluation covariates. The finite change sets
x from −0.5 to +0.5. Gaussian censored/truncated targets are latent means;
AFT response predictions are exp(linear predictor), not expected event time.
Zero-inflated response predictions include the zero-inflation probability.

## Preset thresholds and final maxima

No numerical threshold was loosened during this stage.

| Diagnostic | Required below | Maximum observed |
|---|---:|---:|
| Total log-likelihood absolute difference | 2e-6 | 4.55e-12 |
| Score sum divided by square root of information diagonal | 2e-3 | 8.70e-5 |
| Information stability, Richardson r=4 versus r=6 | 5e-5 | 1.89e-5 |
| Pointwise score stability, r=4 versus r=6 | 5e-5 | 1.41e-6 |
| Full covariance difference / marginal SE products | 3e-4 | 4.89e-5 |
| Prediction contrast estimate absolute difference | 2e-7 | 5.55e-17 |
| Prediction contrast relative SE difference | 3e-4 | 2.75e-7 |
| Finite-change contrast estimate absolute difference | 2e-7 | 1.39e-17 |
| Finite-change contrast relative SE difference | 3e-4 | 1.35e-7 |

Information must also be positive definite and all required outputs finite.
Fitting warnings and native optimizer fields are recorded separately for each
component. Unsupported/missing native convergence fields are not interpreted
as successful flags. The independent stationarity diagnostic applies to all
fits. `gate2-metrics.txt` provides exact maxima and the corresponding cases.

## Finding RV-03: interval-censored survival scale score

The initial 96-case run had four Gaussian interval full-covariance failures.
The largest standardized covariance discrepancy was 0.573622, well beyond
3e-4. Likelihood totals agreed to approximately 1e-12, coefficient score
columns agreed, and native bread agreed with inverse independent information
to 3.81e-13. The installed survival 3.8-3 residual routine supplied the
interval log-scale score with the opposite sign to the likelihood derivative.
The native sandwich method inherited that error. Aggregate scores near zero
did not reveal it, and latent-mean contrast SEs were much less sensitive than
the nuisance/cross-covariance entries.

For an interval with standardized limits a and b, the correct derivative is

    d log{F(b) − F(a)} / d log(scale)
      = {a f(a) − b f(b)} / {F(b) − F(a)}.

`suest` now computes this derivative directly for interval-censored rows with
estimated scale, using the native distribution's transformation and density.
It chooses a survivor-function difference in the upper tail. It does not
simply negate an upstream score, so a future upstream correction will not
reverse the result again. Exact, one-sided, and fixed-scale rows retain their
existing treatment; conventional robust fits retain model-based bread.

Regression tests failed before the correction and pass after it. They cover
mixed exact/left/right/interval observations, Gaussian and lognormal models,
and estimated/fixed scale, checking individual scores and public full
covariance. The standalone gate additionally covers Weibull and loglogistic
interval models. The prior test that used the affected native sandwich as its
interval reference was replaced; stronger independent derivative/covariance
assertions now supply that evidence.

## Estimator boundaries

Actual ML, bias-reduced (BR), and bias-corrected (BC) beta fits were probed.
Only ordinary ML is certified by these likelihood/information references.
BR uses adjusted native estimating equations and BC changes the estimator;
neither can be silently treated as the ordinary-ML contract. Both are now
rejected explicitly. This is a support restriction, not a claim that every
native BR score is an ordinary ML score or is wrong.

Penalized survreg already failed a generic parameter-order check. It now
receives an explicit unsupported-estimator error. Conventional `robust=TRUE`
survreg fits remain supported; a regression verifies equality to the ordinary
joint fit using model-based information. These changes add no model support
and no production dependency.

## Truncated-regression optimizer finding

The default BFGS route in truncreg 0.2-5 produced transient `log(sigma): NaNs
produced` fitting warnings in all eight cases. One right-truncation case
stopped with an information-scaled score sum of 0.002221968, exceeding 0.002,
despite a native successful-convergence message. Its full covariance still
agreed with the independent reconstruction to 4.58e-9 standardized units.
These results are retained as failures/warnings, not relabeled as clean passes.

Increasing iterations and supplying a tighter `control` did not cure this:
inspection showed that truncreg.fit did not forward that control to maxLik.
The supported `method="NR"` optimizer on exactly the same data and likelihood
converged with no warnings in all eight cases. Its largest stationarity
measure was 9.69e-7. The final 108-case run uses this optimizer, and the gate
accepts `truncreg_method="BFGS"` to reproduce the default-optimizer diagnostic.
No suest truncated-regression formula changed. These results do not certify
the convergence of arbitrary native truncated-regression fits or optimizers.

## Verification and review

- New extended-family safeguard file: 18 assertions pass; before-fix failures
  are retained in `gate2-regressions-before.log`.
- Focused beta, censored/survival, and both safeguard files: pass.
- Independent source review: no critical or important findings. Reviewer
  checked the interval derivative and additional native distributions/offsets.
- Final full source suite: **1,667 assertions pass, 0 failures/errors/warnings/skips**,
  across 214 test blocks (totals recovered from the saved testthat result object).
- Numerical acceptance suite: **144 pass, 0 fail**.
- Final source build and `R CMD check --no-manual`: **Status: OK**,
  no errors, warnings, or notes; installed-package tests also report
  **1,667 pass, 0 fail/warn/skip**.

Runtime: Linux x86_64, R 4.5.2; survival 3.8-3, sandwich 3.1-3,
betareg 3.2-6, censReg 0.5-38, truncreg 0.2-5, Rchoice 0.3-6,
pscl 1.5.9, marginaleffects 1.0.0, numDeriv 2016.8-1.1.
Full session information and source fingerprints are in the gate log.
No current cross-platform CI run is claimed.

## Reproduction and scope

From the candidate package source root in R:

    source("tools/release-validation/suest_release_gate2.R")

For focused execution or a different output directory:

    options(suest.gate2.define_only = TRUE)
    source("tools/release-validation/suest_release_gate2.R")
    run_suest_release_gate2(output = "path/to/results",
      selected = c("truncreg_left", "truncreg_right"),
      truncreg_method = "BFGS")

`GATE2_COMPLETED=1` only means execution finished. Inspect every status and
warning. The first pilot had a mismatched diagnostic label in the harness;
it was corrected before the 96-case baseline. Pilot statuses are not used as
package evidence. Score-stability recording and the 12 additional transformed
interval cases were added after review. The preserved baseline, post-fix
BFGS, focused NR, and final logs make these stages distinguishable.

This gate covers paired models within each specified route. It does not
exhaust mixed-family pairs, custom distributions, offsets, extreme-tail
intervals, weights, or arbitrary optimizers. Separate tests cover some of
those features; this gate does not promote them to newly certified results.
The next bounded stage is independent contrast-SE validation for the remaining
panel/integrated and survey routes identified in the evidence matrix.
