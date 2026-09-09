# SUEST R survey binary-response session closeout

Implementation date: 2026-09-08
Final validation: 2026-09-09
Package: suest 0.1.4 certified candidate

## Outcome

This session implemented the narrow survey-weighted binary-response increment
for `survey::svyglm()` using `quasibinomial("logit")` and
`quasibinomial("probit")` under the already certified one-stage common-design,
coefficient-only, explicit-observation-ID contract.

The implementation, focused test additions, deterministic Stata benchmark,
returned-Stata regression fixture, and complete R validation sequence are
finished. The increment is **finally certified** within its stated scope.

## Statistical construction

Binary coefficient influences are reconstructed directly from each model's GLM
estimating equation and native R bread. Model influences are zero-padded onto
the full sampled design and combined through the common survey PSU/stratum/FPC
linearization. Every returned system is required to reproduce each component
model's native `svyglm()` coefficient covariance block exactly before forming
the cross-model covariance.

Focused tests were added for identical samples, partial overlap, disjoint
domains, FPCs, factors/interactions, offsets, influence reconstruction,
prior-weight rescaling, predictions, and refusal conditions.

## Returned Stata 19.5 gate

The user returned a complete six-case benchmark log.

- Logit: full numerical parity with the R estimating-equation construction and
  survey linearization (differences around `1e-10` or smaller).
- Probit: point estimates agree, but finite-sample covariance differs because
  native R `svyglm()` uses expected/Fisher-information bread while the returned
  Stata survey-probit VCE is reproduced by an observed-information bread.
- Same-subpopulation systems are accepted by `suest2`; systems with different
  `subpop()` expressions return the previously documented code 322.
- The disjoint-strata case has the expected exactly zero cross-model block.

The R probit route intentionally retains native `svyglm()` covariance rather
than forcing Stata VCE equality.

## Session boundary

Final R results were: 132/132 focused expectations; 147 tests and 872/872 full-
suite expectations; 144/144 numerical acceptance cases; successful source
build with vignette; successful clean install/load; and `R CMD check
--no-manual` status `OK` with no errors, warnings, or notes.

No Stata `suest2` source was edited. No GitHub change was made. No further Stata
run is presently required.

The next development session should start from the certified 2026-09-09 source
and treat broader survey features as separate, newly authorized increments.
