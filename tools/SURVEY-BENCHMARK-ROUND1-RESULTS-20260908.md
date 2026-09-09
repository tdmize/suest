# Survey benchmark round 1 results

Date: 2026-09-08  
Stata: 19.5  
`suest2`: 1.0.0, 30aug2026  
Benchmark revision: 1  
Completion marker: present

## Result

All 12 native Stata `svy: regress` coefficient and covariance blocks agree
closely with their R `survey::svyglm()` counterparts. Across the six designs
and two models per design:

- maximum absolute coefficient difference: 2.14e-8;
- maximum absolute covariance-entry difference: 8.92e-10;
- maximum relative covariance-diagonal difference: 1.09e-7.

These differences are consistent with revision 1 allowing Stata to select
storage types while importing the CSV. Revision 2 uses `asdouble`.

The native comparison therefore validates the R point estimates, each fitted
model's linearization, strata corrections, first-stage FPC handling, and domain
handling at the model-block level.

## Why no joint `suest2` result was returned

Every revision-1 case stored model A with `subpop(dom_a)` and model B with
`subpop(dom_b)`. Although both variables equal one in cases 1 and 2, `suest2`
compares the specifications rather than their evaluated values. It rejected all
six systems with return code 322:

> survey models use different subpopulation specifications

This is not evidence of disagreement in the R joint covariance. It means the
first gate did not test that covariance through `suest2`.

## Revision-2 response

The revised gate:

1. uses exactly `subpop(dom_a)` for both models in same-sample cases 1 and 2,
   allowing the supported `suest2` route to be tested;
2. retains genuinely different domains in cases 3-6, where return code 322 is
   currently expected from `suest2`;
3. adds an independent stacked `svy: regress` reference in all six cases. The
   stacked model has equation-specific intercept and slope columns, so its
   normal equations reproduce the separate coefficient estimates while its
   joint PSU totals and stratum centering identify the cross-model covariance;
4. imports floating-point CSV columns with `asdouble`.

Full survey parity remains pending revision 2. Do not yet turn the returned
values into permanent cross-language fixtures.

