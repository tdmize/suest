# Survey linear-model cross-language certification

Date: 2026-09-08
R route: one-stage Gaussian identity `survey::svyglm()`
Stata: 19.5
`suest2`: 1.0.0, 30aug2026
Benchmark revision: 2
Completion marker: present

## Certification result

The initial survey-weighted linear-model addition passes its numerical
cross-language gate.

Across all six deterministic designs, the stacked Stata `svy: regress` system
and the R joint survey linearization agree to essentially machine precision:

- maximum absolute coefficient difference: 4.44e-15;
- maximum absolute covariance-entry difference: 5.03e-17;
- maximum relative covariance-diagonal difference: 7.11e-15.

The six cases cover identical samples with and without FPCs, partial overlap,
disjoint observations within shared PSUs, disjoint PSUs within shared strata,
and disjoint strata. The last case correctly returns an exactly zero
cross-model covariance block.

In identical-sample cases 1 and 2, `suest2` accepts the two stored survey models.
Its joint coefficient vector agrees with the stacked Stata reference within
9.99e-16, and its complete joint covariance agrees within 1.04e-17.

In cases 3-6, `suest2` returns code 322 because the models use different
`subpop()` specifications. This is an explicit Stata-interface restriction,
not a covariance discrepancy. The mathematically equivalent stacked survey
regression supplies the full Stata joint reference for these cases.

The `asdouble` revision also removes the small round-one import differences.
Every native Stata coefficient/VCE block now agrees with the corresponding
block of the independently reconstructed joint system at numerical precision.

## Scope of the conclusion

The certification applies only to the documented phase-one contract:

- ordinary one-stage `svydesign` objects;
- Gaussian identity-link `svyglm()` models;
- a common full design with explicit observation IDs;
- stratification and optional first-stage FPCs;
- overlapping, partial, and disjoint model samples;
- supported lonely-PSU policies `fail`, `remove`, and `certainty` under the
  documented restrictions;
- coefficient covariance only, without a survey `lnvar` ancillary parameter;
- asymptotic-normal inference by default.

It does not certify replicate-weight, multistage, two-phase, PPS, calibrated,
raked, or post-stratified designs, nonlinear survey families, survey ancillary
parameters, or automatic survey t/F inference.

In the disjoint-strata case, each native Stata model reports 10 residual design
degrees of freedom while the stacked full system reports 20. The R object stores
the full design's 20 degrees of freedom as metadata but does not apply them
automatically. This is consistent with the documented inference contract.

## Permanent evidence

- Returned log: `suest_r_survey_linear_benchmark_v2.log`.
- Deterministic data: `suest_r_survey_linear_benchmark.csv`.
- Stata gate: `suest_r_survey_linear_benchmark.do`.
- R reference generator: `run_survey_linear_reference.R`.
- Fixture builder: `build_survey_stata_fixture.R`.
- Permanent fixture: `tests/testthat/fixtures/survey-stata.rds`.
- Regression test: `tests/testthat/test-stata-survey-reference.R`.

No R implementation change was required after either returned Stata log. The
logs validated the implemented joint linearization and exposed only the
documented `suest2` different-subpopulation limitation.
