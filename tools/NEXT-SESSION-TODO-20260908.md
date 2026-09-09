# SUEST R next-session TODO — survey binary validation completed

The validation work listed below was completed on 2026-09-09. Retain this file
as the audit checklist; do not rerun it merely to reopen the certified survey
binary-response increment.

## First task: validate the current candidate in an R-capable Work session

The final authority is
`suest_r_0.1.4_survey_binary_certified_20260909.zip`, superseding the validation
handoff ZIP.
Do not restart the survey binary implementation from the older linear closeout.
Read `SURVEY-BINARY-STEP-STATUS-20260908.md`,
`SURVEY-BINARY-BENCHMARK-RESULTS-20260908.md`, and
`SESSION-HANDOFF-20260908.md` first.

Completed in this order:

1. Extract/load the package and run the focused survey tests, especially
   `tests/testthat/test-survey.R` and
   `tests/testthat/test-stata-survey-binary-reference.R`.
2. Confirm that every binary model's diagonal block exactly reproduces its
   native `survey::svyglm()` covariance before any joint covariance is returned.
3. Confirm the logit/probit influence reconstruction, prior-weight rescaling,
   domain padding, FPC, factor, offset, and prediction tests actually pass in R.
4. If any focused test fails, make the smallest targeted correction and rerun
   the focused survey gate before broad testing.
5. Once focused tests are stable, run the complete unit suite.
6. Run the numerical acceptance suite.
7. Run a standard `R CMD build` with vignette creation.
8. Install the resulting source tarball into a clean library and load it.
9. Run `R CMD check --no-manual` and require no errors/warnings/notes unless an
   environmental qualification is clearly demonstrated and documented.
10. Reviewed the diff and updated the binary closeout/status documents with the
    actual R test counts, runtime/package versions, and final package-check
    result.

All ten items passed. Exact results are in
`SURVEY-BINARY-FINAL-VALIDATION-RESULTS-20260909.md`.

## Cross-language interpretation already established

Do not require identical Stata covariance for native-R survey probit as a final
acceptance condition. The returned Stata 19.5 gate establishes:

- survey logit: full numerical parity;
- survey probit: point-estimate parity, but finite-sample covariance differs
  because native `svyglm()` uses expected/Fisher-information bread and Stata's
  returned survey-probit VCE is reproduced by observed-information bread.

The R package must preserve native `svyglm()` covariance. Do not add a Stata-
compatibility bread override to the native survey-probit route.

No further Stata run is needed unless an R correction materially changes the
construction. The existing different-`subpop()` `suest2` code-322 behavior is
an interface restriction, not an R defect.

## Scope restrictions for this increment

Do not broaden the binary increment to:

- replicate-weight designs;
- multistage or two-phase designs;
- PPS designs;
- calibration, raking, or post-stratification;
- survey ancillary parameters;
- automatic survey t/F inference;
- explicit survey `lnvar`;
- broader lonely-PSU policies;
- uncertainty from averaging covariates;
- nonlinear-panel or multilevel models.

Discuss and benchmark those separately after this increment is closed.

## Dependency result corrected during release preflight

The earlier `marginaleffects` 1.0.0 `%||%` failure came from an incomplete
dependency installation, not from `suest` or R 4.3.3. With a complete clean
library, a plain `lm()` prediction succeeds and the `suest` focused survey,
complete unit, and 144-case numerical acceptance suites pass under
`marginaleffects` 1.0.0. One survey assertion was updated to compare the
numeric standard error while ignoring the new `jacobian_source` metadata
attribute. No package calculation or dependency code was patched.

## Post-release roadmap priority updated 2026-09-09

1. nonlinear panel families, one at a time;
2. general GLMM/multilevel support, beginning with a narrow two-level
   random-intercept family;
3. MI estimand-level pooling and all other extensions listed in
   `ROADMAP-POST-0.1.4-20260909.md` at lower priority.

## Separate Stata `suest2` issues

These remain documented findings and are not part of this R increment:

- help syntax shows exactly two models although the command accepts three;
- help/code disagreement about `xtreg, cre`;
- changelog version mismatch and `seemlessly` spelling error;
- pweighted linear `lnvar` changes under common weight rescaling;
- mixed-family unequal-sample PA score extraction fails until data are reduced
  to `e(sample)`;
- differing survey `subpop()` specifications are refused with code 322.

Do not edit Stata source unless explicitly authorized.

## Git/release boundary

Do not commit, push, or modify GitHub unless explicitly authorized. After final
R certification, prepare a clean closeout and wait for release instructions.
