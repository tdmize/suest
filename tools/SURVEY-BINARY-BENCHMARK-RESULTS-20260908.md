# Survey binary-response benchmark results — 2026-09-08

Package: suest 0.1.4 development candidate  
R target: `survey::svyglm()` with `quasibinomial("logit")` or `quasibinomial("probit")`  
Stata gate: Stata 19.5, `svy: logit` / `svy: probit`, `suest2` 1.0.0

## Returned gate

The user-returned benchmark completed all six deterministic systems. Cases 1,
3, and 5 use logit; cases 2, 4, and 6 use probit. The cases cover identical
samples, partial overlap, disjoint domains, first-stage FPCs, factors, and
formula offsets. The full common one-stage design is retained in every case.

For cases 1 and 2, `suest2` accepts the same-subpopulation systems and returns
code 0. For cases 3–6, `suest2` returns code 322 because the model-specific
`subpop()` expressions differ, matching the already documented survey
interface restriction. The independently stacked survey GLM therefore supplies
the complete Stata joint reference for all six cases.

## Logit parity

The logit systems use the same coefficient estimating equation and Fisher bread
as the R `svyglm()` route. Independent reconstruction of the R estimating
functions and one-stage survey linearization matches the returned Stata stacked
systems extremely closely across cases 1, 3, and 5:

- maximum absolute coefficient difference: 1.79e-10;
- maximum absolute covariance-entry difference: 6.93e-11;
- maximum relative covariance-diagonal difference: 5.98e-10.

This is classified as full numerical parity for the narrow survey-logit route.

## Probit covariance convention

The probit coefficient estimating equation is also common across R and Stata,
and returned coefficients agree to ordinary optimizer tolerance (maximum
absolute difference about 4.45e-8). The finite-sample sandwich covariance is
not numerically identical because the fitting engines use different breads:

- `survey::svyglm()` delegates the fit to `stats::glm()`, stores
  `summary(glm)$cov.unscaled`, and uses that expected/Fisher-information inverse
  in its linearized influence and covariance;
- The returned Stata survey-probit covariance is reproduced by an observed-information Hessian bread.

Independent calculations reproduce the returned Stata probit covariance after
substituting the observed Hessian: maximum covariance-entry discrepancy is
1.22e-9 across cases 2, 4, and 6. Using the R/Fisher bread instead produces the
expected finite-sample differences, with a maximum returned benchmark diagonal
difference of about 7.89%.

This is an estimator-convention difference, not a score, domain-alignment, FPC,
factor, or offset error. The R package deliberately preserves the native
`svyglm()` covariance because the survey route requires exact reproduction of
each native R covariance block before a joint system is returned. No Stata-
compatibility override is added.

## Final R regression gate

The permanent benchmark regression file passed all 22 expectations under the
final R candidate. Its native-R covariance-block checks remain effectively
exact (maximum absolute errors at floating-point rounding scale). The logit
joint-covariance comparisons remain within sub-micro absolute tolerance; the
probit checks continue to require coefficient agreement only across engines,
as established above. No Stata source or benchmark input was changed.

## Files

- Deterministic data: `tools/stata-benchmarks/suest_r_survey_binary_benchmark.csv`
- Stata gate: `tools/stata-benchmarks/suest_r_survey_binary_benchmark.do`
- Returned log reviewed: `suest_r_survey_binary_benchmark.log`
- Machine-readable Stata reference:
  `tests/testthat/fixtures/survey-binary-stata-reference.csv`
- Regression test: `tests/testthat/test-stata-survey-binary-reference.R`
