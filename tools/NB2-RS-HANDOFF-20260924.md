# NB2 random-slope development checkpoint

## Current state

The restricted NB2-log glmmTMB random-intercept plus one correlated numeric
random-slope increment is complete in development version 0.1.5.9000.
Balanced, partial-overlap/higher-clustered, and disjoint Stata returns have
been compared and independently audited. No more Stata input is needed
for this increment. No production R code changed in its final validation
step. No remote push or release was performed.

The working branch is `dev/glmmtmb-nbinom2-rs`. Use the bundle ref and
`START_HERE.md` in the accompanying ZIP for the exact final commit. Preserve
the single working repository rather than creating extra checkout folders.

The current report is `GLMMTMB-NBINOM2-RS-VALIDATION-20260924.md`. Earlier
balanced, partial, and development reports are historical checkpoints.
Raw GSEM covariance differences are explicitly retained; do not replace
raw fixtures with independent covariance or claim exact raw Stata parity.

## Important audit detail

The disjoint Fisher-coordinate audit first missed the absolute covariance
gate (1.477e-6 versus 1e-6). Matching glmmTMB's native correlation coordinate
requires the nonstationary Hessian chain-rule term. With `z = atanh(rho)`,
`t = sinh(z)`, and log-likelihood gradient `g_z`, the transformed native
inverse information is `inverse(I_z + g_z[7] * rho * E77)`. The sandwich
formed from this bread and independent score cross-products agrees with
R joint covariance within 3.872e-8. No numeric tolerance or package
covariance was changed. Fisher, natural-variance, and native-coordinate
references are all separately retained.

The same refinement was independently rerun for balanced and partial cases.
Matching-coordinate SE differences are below 0.000066% for every case.
For disjoint samples, Stata's union correction is 100/99 and R's component
correction is 50/49; preserve the explicit ratio 99/98 in comparisons.

## Verification record

The focused Stata regression passed 78 expectations with no failures,
warnings, or skips. The full package suite passed 1,475 expectations with
zero failures and zero skips, plus three test warnings (the earlier
implementation checkpoint also recorded three). R CMD build succeeded;
full `R CMD check --no-manual`, including tests and vignettes, finished
with **Status: OK**. Captured outputs accompany the ZIP.

Environment: R 4.5.2 on Linux, glmmTMB 1.1.14, TMB 1.9.25,
marginaleffects 1.0.0, and testthat 3.3.2. GitHub CI has not been run for
this unpushed checkpoint.

## Next increment

The previously proposed next family is binomial-logit with the same single
correlated numeric-slope structure. It is not implemented by this checkpoint.
Keep it separate from NB2: establish the likelihood/parameter contract,
integrated probabilities and effects, and independent score/covariance
checks in R before requesting Stata comparisons. Use a small balanced gate
first, with checked warm starts and short visible iteration limits; then
partial and disjoint gates. Do not return to long cold-start joint fits.

Alternatively, this checkpoint can be prepared for a release before adding
another family. Versioning, remote integration, and GitHub publication remain
separate work; obtain Trent's explicit push authorization. Do not infer it
from his earlier push of version 0.1.5.

## Reproduction and working preferences

Key files are `test-glmmtmb-nbinom2-rs.R`, the three-case
`test-glmmtmb-nbinom2-rs-stata.R`, and separate raw/audited RDS fixtures.
Benchmark builders, comparison scripts, independent audit helpers, and
bounded Stata scripts live in `tools/stata-benchmarks`; run reproduction
scripts from that directory. The raw CSV and fitted R reference are included
in the flat ZIP; the full-history bundle preserves the entire repository.

Trent runs Stata 19.5 on Windows. Run R verification here. Deliver one flat
ZIP rather than individual downloads. For any future Stata gate give the
exact `do ...` command and expected log filename. Keep code compact and
put short purpose comments to the right of new `library(...)` calls.
