# Prompt for the next development session

Continue development of the R version of `suest` from the authoritative
certified source ZIP:
`suest_r_0.1.4_survey_binary_certified_20260909.zip`.

Read `tools/SESSION-HANDOFF-20260908.md`,
`tools/NEXT-SESSION-TODO-20260908.md`,
`tools/SURVEY-BINARY-STEP-STATUS-20260908.md`, and
`tools/SURVEY-BINARY-BENCHMARK-RESULTS-20260908.md` before making changes.

The narrow one-stage survey binary-response logit/probit implementation, Stata
19.5 benchmark, focused R tests, complete unit and numerical acceptance suites,
source build with vignette, clean installation/load, and `R CMD check
--no-manual` are complete. Read
`tools/SURVEY-BINARY-FINAL-VALIDATION-RESULTS-20260909.md` before selecting a
new increment. Do not reopen or broaden the certified survey route without
explicit authorization.

Important established result: survey logit has full R/Stata numerical parity.
Survey probit coefficients agree, but native R and Stata finite-sample VCEs use
different bread conventions: `svyglm()` uses expected/Fisher information while
the returned Stata covariance is reproduced by observed information. Preserve
exact native `svyglm()` covariance; do not force Stata VCE parity for probit.

Retain the existing one-stage common-design, coefficient-only, explicit-ID
contract. Do not broaden this increment to replicate weights, multistage or
two-phase designs, PPS, calibration/raking/post-stratification, survey ancillary
parameters, or automatic survey t/F inference. Do not edit Stata `suest2`
source or make GitHub changes unless I explicitly authorize it.

Proceed in small stages and preserve the certified validation baseline.
