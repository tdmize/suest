# Remaining contrast and convergence validation — 2026-09-29

Candidate: suest 0.1.6.9000, baseline 29ae851 plus combined local corrections.
No push and no new Stata run. All R work was performed by the assistant.

## Scope and references

The gate covers 28 cases, four fixed sample patterns for each of seven routes:
matched glm/glm2 logit, probit and Poisson; matched polr/clm ordered logit and
probit; bivariate probit; and unweighted fixest 2SLS. Distinct outcomes are used
for models A and B. Patterns are common observations, unequal partial samples,
disjoint observations sharing clusters, and disjoint clusters. The generated
sample has 480 observations and 60 clusters. All patterns retain 60 clusters
in their union, while partial samples omit some union observations.

There are 132 A-minus-B contrasts: average predictions, x finite changes from
-0.5 to +0.5, and average x slopes. Ordinal targets cover each of three categories.
The evaluation population is 40 fixed covariate rows, with x shifted by +0.2;
no uncertainty in estimating this evaluation distribution is claimed.

Twenty cases compare matched engine systems, aligning coefficients and the
full covariance by name because polr and clm store thresholds in different
orders. These checks establish adapter equivalence, not a new independent
covariance certification. Earlier GLM/ordinal covariance references remain
necessary supporting evidence.

Four bivariate-probit cases independently differentiate the signed-outcome
bivariate normal log likelihood in beta/athrho coordinates. No native scores
or Hessians enter the reference. The bivariate CDF primitive is shared with
mvtnorm; derivative and covariance construction are separate. Observation
influences are aligned by explicit IDs, aggregated by cluster, and corrected
by the union-cluster G/(G-1). Prediction and finite-change parameter gradients
are analytic and checked against numerical derivatives. The x slope is
analytic, with its parameter gradient independently differentiated numerically.
The target is joint success, not a marginal outcome probability. Rho and
cross-model covariance sensitivity diagnostics are retained.

Four 2SLS cases independently project X onto Z, construct observation influences,
and aggregate with the existing uncorrected HC0 convention. All three mean
contrast gradients are analytic and verified numerically. This dataset is
just identified; existing overidentification tests remain separate evidence.

Preset limits: engine coefficient/estimate gap 2e-5, standardized covariance
gap 2e-4, relative contrast SE gap 5e-4. Independent bivariate/IV limits are
2e-5, 2e-6, and 3e-4 respectively. No limits were changed after results.

## Admission defects and corrections

Real low-iteration multinom and polr fits returned convergence code 1 but were
accepted. A clm fit returned code -1 (iteration limit / excessive gradient),
and a glm.nb fit returned `th.warn = "alternation limit reached"`, even though
its conditional GLM step reported convergence. Both were accepted too.

Five regression expectations failed before the fix; all 13 expectations pass
after, including successful controls. Admission now requires multinom/polr
code zero, clm diagnostic code zero, and glm.nb mean convergence without a
retained dispersion warning. clm positive codes denote unreliable information
or precision and are intentionally excluded as well. These guards do not
replace general numerical diagnostics or prove arbitrary fits well conditioned.
No covariance, likelihood, or prediction formula changed in this stage.

## Results and verification

All **28/28 cases meet the numerical gates**, covering **132 contrasts**.
Twenty-seven cases are warning-free; one bivariate partial-sample case is
`PASS_REVIEWED_FITTING_WARNING` after the supplemental review below.

| Quantity | Largest gap |
|---|---:|
| Engine-equivalence standardized full covariance | 3.756e-7 |
| Engine-equivalence contrast estimate | 1.007e-7 |
| Engine-equivalence relative contrast SE | 5.845e-6 |
| Independent bivariate/IV standardized full covariance | 2.475e-8 |
| Independent bivariate/IV contrast estimate | 3.803e-12 |
| Independent bivariate/IV relative contrast SE | 2.554e-7 |

Omitting rho uncertainty changes a bivariate prediction-contrast SE by 18.56%;
omitting cross-model covariance changes one by 22.04%. The new evidence is
sensitive to both contributions. The exact per-target maxima are in metrics.txt.

Five native warnings in model B's partial-sample fit were traced to
`log(pmvnormWrap(...))` at trial rho near -0.862385. Replaying each captured
trial gives a TVPACK probability between -1.016e-20 and -1.694e-21; independent
conditional-normal integration gives positive probabilities between 7.596e-24
and 3.546e-23. These are roundoff-negative tail probabilities at trial values,
not invalid final fitted correlations. Both final native optimizers return code
zero; native final log likelihoods reproduce the fitted values exactly,
likelihood/gradients are finite and warning-free, and maximum score sums are
below 2e-3. Independent full covariance and contrast checks pass unchanged.
The focused replay reproduces all three original metric vectors within 1e-12.
No warning was discarded from the raw logs or silently treated as warning-free.

The combined results retain the 24 engine/IV cases, three unchanged bivariate
cases, and the focused reviewed partial case. A default script run deliberately
still reports `REVIEW_WARNING` when native warnings occur; the supplemental
scripts and logs justify the final reviewed classification.

Full source suite: **1,680 passing assertions**, zero failures/errors/warnings/skips.
Numerical acceptance suite: **144 pass, zero fail** after correcting the
nonconverged fixture. Installed-package tests: **1,680 pass**, zero failures,
warnings or skips. `R CMD check --no-manual`: **Status: OK**, no errors,
warnings or notes, on Linux/R 4.5.2. The frozen built R/Rd/tests/README/NEWS
are verified against current source by the packet builder. Audit tools are
excluded from the built package; the revised acceptance script is separately
verified by its executed results and the combined patch.

## Review and execution history

Independent read-only review found no blocking reference or guard issue.
It checked coordinate ordering, signed likelihood, athrho chain rule, cluster
corrections, projected-X influences, and native convergence-code semantics.

The pilot detected ordinal parameter ordering, corrected by named alignment.
It also encountered a retained bivariate formula referencing the temporary
loop variable j; the script now stores the actual formula object in the model
call. These were harness errors, not adjusted statistical comparisons.
Original pilot and boundary-probe logs remain in the packet.

The broader acceptance suite exposed an old nonconverged iris multinomial
fixture. It was replaced with seeded overlapping categories while retaining
probability sums and effect sums, and now explicitly checks both fits converge.
The initial 143/144 result and subsequent 144/144 result are preserved.

## Release implication and remaining work

This closes the specified alternative-engine and remaining contrast checks
within the tested scope. It is not universal certification of arbitrary
links, datasets, covariance conditioning, or estimator options. Preserve the
FE average-level limitation and the raw random-slope GSEM discrepancies from
prior stages, including the failed absolute-tolerance diagnostic.

Next freeze the exact candidate and verify current supported-platform CI,
including Windows/macOS and available dependency versions. Historical baseline
CI is not a substitute. GitHub publication still requires the user's permission.
No new model support or nonlinear MI inference is included in this pass.

Reproduce from the package root with:

    source("tools/release-validation/suest_release_gate4.R")

`GATE4_COMPLETED=1` means execution finished; inspect every case status and
warning. Focused runs use `options(suest.gate4.define_only=TRUE)` before source,
then `run_suest_release_gate4(output="...", selected="biprobit")`.
