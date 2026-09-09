# SUEST R package development handoff — survey binary-response increment

Implementation date: 2026-09-08
Final R validation: 2026-09-09
Package version: 0.1.4 certified candidate
Authoritative source for the next session: `suest_r_0.1.4_survey_binary_certified_20260909.zip`

## Immediate status

The initial one-stage Gaussian identity `survey::svyglm()` increment remains
fully certified. This session implemented the next narrow survey increment:
binary-response `svyglm()` models using `quasibinomial("logit")` and
`quasibinomial("probit")`, while preserving the existing one-stage
common-design, coefficient-only, explicit-observation-ID contract.

The binary implementation, deterministic Stata gate, and complete R validation
sequence are finished. The increment is **finally certified** within the
documented narrow contract. The final R gate used R 4.3.3, `survey` 4.5, and
`marginaleffects` 0.31.0.

Read these files before changing the binary route:

- `tools/SURVEY-STEP-STATUS-20260908.md` — certified Gaussian survey contract;
- `tools/SURVEY-BINARY-STEP-STATUS-20260908.md` — binary implementation and
  influence formula;
- `tools/SURVEY-BINARY-BENCHMARK-RESULTS-20260908.md` — returned Stata gate and
  the logit/probit parity result;
- `tools/NEXT-SESSION-TODO-20260908.md` — exact next steps.

## Binary implementation completed this session

The production survey route now accepts systems in which all models are
`survey::svyglm()` fits from one ordinary one-stage `svydesign()` and are all
one of:

- Gaussian identity models under the previously certified linear route; or
- `quasibinomial("logit")` models; or
- `quasibinomial("probit")` models.

For this increment, Gaussian and binary survey models are deliberately not
mixed within one system, and logit/probit mixing is also kept outside the
narrow supported contract unless the code/tests explicitly permit it after
validation. Existing restrictions on replicate weights, multistage/two-phase
sampling, PPS, calibration/raking/post-stratification, survey ancillary
parameters, extra fitting weights, and automatic survey t/F inference remain.

Algebraically, the binary coefficient influence is as follows. For observation
`i`, with model-matrix row `x_i`, prior weight `w_i`, fitted mean `mu_i`, link
derivative `mu'_i`, and GLM variance `V(mu_i)`, the influence is

`[x_i w_i (y_i-mu_i) mu'_i / V(mu_i)] A^{-1}`

with

`A = sum_i x_i x_i' w_i (mu'_i)^2 / V(mu_i)`.

For logit, the score factor reduces to `w_i (y_i-mu_i)`. For probit it is
`w_i (y_i-mu_i) phi(eta_i) / [mu_i(1-mu_i)]`. Common prior-weight rescaling
cancels from the influence. Quasi dispersion is not multiplied into the native
design-based coefficient covariance. Operationally, the adapter evaluates this
identity from the final working residuals and working weights retained by
`svyglm()`, together with its retained `naive.cov` expected-information bread.
Recomputing the same quantities from the printed coefficients can differ at the
GLM stopping tolerance and fail exact native-covariance reproduction.

As in the certified linear route, observation influences are zero-padded onto
the full common survey design and passed through the common PSU/stratum/FPC
linearization. Critically, the system refuses to return unless each model's
own diagonal covariance block exactly reproduces its native `svyglm()`
covariance first.

## Focused tests added

The survey tests now include binary coverage for:

- direct/native influence reproduction;
- exact native covariance-block reproduction;
- identical samples;
- partial overlap;
- disjoint domains;
- first-stage FPCs;
- factors/interactions;
- formula offsets;
- response- and link-scale prediction behavior;
- common prior-weight rescaling;
- intended unsupported/refusal surfaces.

The returned Stata results are also retained as a permanent regression
reference in `tests/testthat/fixtures/survey-binary-stata-reference.csv`, with
`tests/testthat/test-stata-survey-binary-reference.R` exercising the reference.

## Returned Stata 19.5 benchmark

The user ran:

`do suest_r_survey_binary_benchmark.do`

and returned `suest_r_survey_binary_benchmark.log`. All six deterministic cases
completed. They cover logit and probit with identical samples, partial overlap,
disjoint domains, FPCs, factors, and offsets.

For the same-subpopulation cases, `suest2` returns code 0. For cases with
different `subpop()` expressions, `suest2` returns code 322, matching its
already documented interface restriction. The independently stacked Stata
survey GLM provides the complete joint reference in those cases.

### Logit

The logit systems have full numerical R/Stata parity under the independently
reconstructed R estimating equations and one-stage survey linearization:

- maximum absolute coefficient difference: about `1.79e-10`;
- maximum absolute covariance-entry difference: about `6.93e-11`;
- maximum relative covariance-diagonal difference: about `5.98e-10`.

### Probit

Probit point estimates agree to ordinary optimizer tolerance, but the finite-
sample covariance differs by fitting-engine convention. Native
`survey::svyglm()` uses the `stats::glm()` expected/Fisher-information bread,
whereas the returned Stata `svy: probit` covariance is reproduced by an
observed-information Hessian bread. Substituting the observed bread reproduces
the returned Stata covariance to roughly `1e-9`; using the native R/Fisher bread
produces finite-sample diagonal differences as large as about 7.9% in the gate.

This is **not treated as an R defect**. The package's governing contract is
exact reproduction of native `svyglm()` covariance blocks, so the R probit
route must retain the native R bread rather than force Stata covariance parity.
Document this as an estimator-convention qualification analogous to other
accepted cross-engine differences.

## Final R validation completed 2026-09-09

- Focused survey gate: 132 expectations passed, 0 failed.
- Complete unit suite: 147 tests / 872 expectations passed, 0 failures, 0
  errors, 0 skipped; three expected model-fitting warnings were emitted by
  deliberate convergence/separation stress cases.
- Numerical acceptance suite: 144 passed, 0 failed.
- Standard `R CMD build`: succeeded and created the vignette.
- Clean source-tarball installation and `library(suest)`: succeeded.
- `R CMD check --no-manual`: status `OK`, with 0 errors, 0 warnings, and 0 notes.

See `tools/SURVEY-BINARY-FINAL-VALIDATION-RESULTS-20260909.md` for the commands,
versions, and interpretation. No additional Stata run is required.

## Broader project state to preserve

The package's previously implemented model families, MI pooling contract,
explicit-ID alignment, cross-model score stacking, clustering rules, pweight
conventions, ancillary-parameter handling, and marginaleffects adapters remain
as documented in `tools/MODEL-ADAPTERS.md`, `tools/STATA-PARITY.md`, and the
prior handoff history. Do not reset or replace the development tree with the
old GitHub state.

No Stata `suest2` source was edited during this increment. No Git commit, push,
or other GitHub change was made or authorized.
