# SUEST R: initial survey linear-model increment

Date: 2026-09-08. Package: 0.1.4 development candidate.
This update supersedes the earlier handoff's statement that all survey support
is pending. All earlier implemented families and returned Stata fixtures are
retained. No Git commit, push, or Stata source edit was made.

## Implemented contract

```r
library(survey)
design <- svydesign(~psu, strata = ~strata, weights = ~w, data = dat)
m1 <- svyglm(y ~ x, design)
m2 <- svyglm(y ~ x + z, design)
fit <- suest(m1, m2, survey_design = design, observation_id = "id")
marginaleffects::avg_slopes(fit, variables = "x",
  newdata = suest_newdata(fit), wts = ".suest_weight")
```

- Two or more Gaussian identity-link `survey::svyglm()` models from one common,
  ordinary one-stage `svydesign` design. The full design and observation-ID
  column names are required; composite IDs are supported.
- Stratification, first-stage FPCs, overlapping and disjoint estimation samples,
  and domains are supported. Retain the full design before domain subsetting.
  Missing outcomes and reordered rows are aligned explicitly.
- Parameter vector contains regression coefficients only. Survey `lnvar` is
  deliberately deferred, rather than borrowing the ordinary pweight convention.
- Reconstruct each observation's coefficient influence as
  `(X * (w * residual)) %*% solve(crossprod(X, X * w))`.
  Pad with zeros on the full common design; stack columns across models; apply
  `survey::svyrecvar()` to obtain the complete joint covariance.
- Equivalently, within each stratum, sum influences within PSU, center PSU
  totals, and multiply their cross-product by `G_h/(G_h-1)` and the FPC
  multiplier `1 - G_h/N_h`. There is no extra ordinary HC1 adjustment.
- Before joining models, reproduce each entire native coefficient covariance
  block. All six deterministic benchmark cases agree within 2.43e-17 in maximum
  absolute covariance-entry difference.
- Nonoverlapping observations can have cross-model covariance through shared
  PSUs. Disjoint PSUs in the same stratum can also have cross-model covariance
  because stratum centering couples their totals. Different strata are independent.

## Deliberate restrictions and inference

No replicate-weight, multistage, two-phase, PPS, calibrated, raked, or
post-stratified designs; no mixed survey/nonsurvey systems or nonlinear survey
families. Extra `svyglm(weights=)` fitting weights are rejected. Formula offsets
are supported; the separate `offset=` argument is rejected. Aliases, duplicate
IDs, incompatible fitted-design metadata, and unreproduced native VCEs fail
explicitly. Do not supply `cluster` or `weight_type` on the survey route.

Lonely-PSU options: `fail`, `remove`, or `certainty`, with
`survey.adjust.domain.lonely = FALSE`. A certainty FPC permits a singleton
stratum under `fail`. `adjust` and `average` are deferred because their joint
and model-specific domain handling requires a separate contract.

Predictions and effects use the joint design-based coefficient covariance.
Weighted averaging is explicit through `.suest_weight`. The supplied covariate
distribution is treated as fixed; no additional sampling uncertainty for that
distribution is included. Default inference remains asymptotic normal. The full
design DF is retained in `fit$survey$design_df`, not automatically applied as
survey t/F or model residual DF.

## Verification in this increment

- Focused survey suite: **66 passed, 0 failures, 0 warnings, 0 skips**.
- Permanent returned-Stata survey suite: **31 passed, 0 failures, 0 warnings,
  0 skips**.
- Complete unit suite: **803 passed, 0 failures, 3 expected fitting warnings,
  0 skips**. This is the prior 706-test baseline plus 66 survey-contract checks
  and 31 returned-Stata regression checks.
- Numerical acceptance suite: **144 passed, 0 failures**.
- Tests include native diagonal blocks, independent explicit PSU covariance,
  duplicate/three-model systems, partial/disjoint domains, outside-domain PSUs,
  FPC and certainty strata, sampling-weight rescaling, nested PSU labels,
  composite IDs, missingness, reordering, factors/interactions, offsets,
  coefficient perturbation, weighted average slopes, and cross-model contrasts.
- Standard `R CMD build` with vignette creation: passed.
- Clean installation and loading from the built source tarball: passed.
- `R CMD check --no-manual`: **Status: OK**, with no errors, warnings, or notes;
  installation, namespace, code, documentation, examples, tests, and vignette
  execution/rebuilding all passed.
- Round 1 returned from Stata 19.5. All 12 native coefficient/VCE blocks agree
  closely with R: maximum absolute coefficient difference 2.14e-8, maximum
  absolute covariance-entry difference 8.92e-10, and maximum relative
  covariance-diagonal difference 1.09e-7. Revision 1 did not identify the joint
  covariance because `suest2` rejected the differently named `subpop()`
  expressions even when they evaluated identically. See
  `SURVEY-BENCHMARK-ROUND1-RESULTS-20260908.md`.
- Revision 2 used double-precision CSV import and returned all six independent
  stacked joint systems. The maximum absolute coefficient difference from R is
  4.44e-15 and the maximum covariance-entry difference is 5.03e-17. See
  `SURVEY-BENCHMARK-FINAL-RESULTS-20260908.md`.

Runtime: R 4.3.3, survey 4.5, marginaleffects 0.31.0, insight 1.5.4,
data.table 1.18.6.1. The initial marginaleffects 1.0.0 attempt failed with
`could not find function "%||%"` because its dependency installation was
incomplete. Release preflight later completed the clean library and passed the
full unit and numerical acceptance suites under version 1.0.0. No dependency
was patched. The prediction bridge was
verified using an unmodified compatible release, 0.31.0. Retest 1.0.0 on a newer
R runtime later; this is not evidence of a SUEST-specific defect.

A second native-interface observation: `predict.svyglm(newdata=)` in survey 4.5
omitted the formula offset in the tested case. The SUEST prediction adapter
explicitly evaluates `X beta + offset`; its offset test uses that analytic target.
Numerical marginaleffects slope-SE checks allow 2e-5 relative tolerance for
nested finite differences (observed absolute discrepancy under 7.4e-7).

## Completed revision-2 Stata gate

The returned Stata 19.5 revision-2 log contains the expected terminal marker
`SURVEY_LINEAR_BENCHMARK_V2_COMPLETE=1` and all six cases. Its independent
stacked `svy: regress` systems agree with R's complete joint coefficient vector
and covariance to essentially machine precision:

- maximum absolute coefficient difference: **4.45e-15**;
- maximum absolute covariance-entry difference: **5.04e-17**;
- maximum relative covariance-diagonal difference: **7.11e-15**.

In cases 1 and 2, which use the same estimation sample, `suest2` also accepts
the pair and agrees with the stacked Stata system (maximum covariance difference
**1.05e-17**). Cases 3 through 6 return code 322 because `suest2` refuses the
different `subpop()` specifications. This is a Stata command-interface
restriction, not a covariance discrepancy; the independent stacked Stata
systems identify and validate those four full joint matrices. Case 6's models
occupy different strata and its cross-model covariance block is exactly zero.

The returned matrices are stored as a permanent RDS fixture and exercised by
`tests/testthat/test-stata-survey-reference.R`. See
`SURVEY-BENCHMARK-FINAL-RESULTS-20260908.md` for the complete interpretation.

## Files and reproduction

Main addition: `R/survey.R`; new tests: `tests/testthat/test-survey.R`.
Small integration changes: `R/suest.R`, `R/adapters.R`, `R/suest-newdata.R`.
DESCRIPTION, NEWS, README, Rd files, adapter/parity trackers, handoff, and TODO
were updated. All prior source files remain in the complete source ZIP.

From the package root with compatible installed dependencies:

```r
pkgload::load_all(".", quiet = TRUE)
testthat::test_file("tests/testthat/test-survey.R", reporter = "summary")
testthat::test_dir("tests/testthat",
  filter = "^(survey|core|pweights|cluster|marginaleffects|multiple-models|mi)$",
  reporter = "summary")
source("tools/stata-benchmarks/make_survey_linear_data.R")
source("tools/stata-benchmarks/run_survey_linear_reference.R")
```

The gate ZIP includes run instructions, CSV/do-file, R generator/reference,
this status note, and current validation evidence. It is not a replacement for
the older session's ancillary archive, which was not supplied in this turn.

## Primary interfaces inspected

- Survey 4.5 source, `R/survey.R` (`svyglm`, `svy.varcoef`, influence extraction):
  https://github.com/cran/survey/blob/master/R/survey.R
- Survey 4.5 source, `R/multistage.R` (PSU totals, centering, FPCs, lonely PSUs):
  https://github.com/cran/survey/blob/master/R/multistage.R
- Author's variance documentation:
  https://r-survey.r-forge.r-project.org/pkgdown/docs/reference/svyrecvar.html
- Author's survey options documentation:
  https://r-survey.r-forge.r-project.org/pkgdown/docs/reference/surveyoptions.html
