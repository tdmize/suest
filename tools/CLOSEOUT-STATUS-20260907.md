# SUEST R 0.1.4 closeout status

## Release-candidate verification

- Unit suite: 706 passed, 0 failed, 0 errors, 3 expected fitting warnings.
- Numerical acceptance suite: 144 passed, 0 failed.
- Source build: passed.
- Package check: installation, namespace, code, documentation, examples, and
  tests passed. The two warnings are limited to unavailable `knitr` and
  `rmarkdown` and the intentionally skipped vignette build.
- `git diff --check`: passed.
- GitHub has not been modified.

## Cross-language results already settled

Returned Stata 19.5 benchmarks verify pweighted systems, 2SLS, linear `lnvar`,
joint clustering, panel FE/BE, balanced RE, ten single-family GEE systems,
bivariate probit, and IV probit. The fixed numerical references are included in
the R tests. Bivariate- and IV-probit curvature is recomputed from analytic
score derivatives because the upstream fitted Hessians lose precision.

## Final Stata gate completed

The returned `suest_r_followup_diagnostic.log` completed. It established:

1. the mixed-family unequal-sample GEE failure is a repairable `suest2` bridge
   issue: native score prediction fails with excluded within-panel rows still
   present and succeeds after keeping only `e(sample)`;
2. the small panel-ML covariance difference is substantially inherited from
   the native fitting-engine covariance and is closed as numerical convention/
   precision rather than a score defect;
3. unbalanced Swamy-Arora RE is close but not exact across engines, with a
   maximum 0.1504% relative covariance-diagonal difference.

## Subsequent MI increment

The closeout candidate now includes coefficient-level multiple-imputation
pooling through `suest_mi()`. Ordinary fit lists, `mice::mira`,
`mice::getfit()`, `mitools::with.imputationList()`, and legacy
`imputationResultList` inputs are supported. Rubin pooling retains the complete
within- and between-imputation cross-model covariance. Nonlinear predictions,
slopes, and comparisons remain deliberately deferred until they can be
estimated within each imputation and pooled correctly.

## Final closeout

All classified cross-language items have permanent regression coverage and the
final 2026-09-08 acceptance run passed 144/144. Broader survey,
nonlinear-panel, and general multilevel support remain future-version projects,
not blockers for this checkpoint. See `SESSION-HANDOFF-20260908.md` and
`NEXT-SESSION-TODO-20260908.md` for the clean restart point.
