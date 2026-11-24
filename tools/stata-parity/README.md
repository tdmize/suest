# Stata parity replays

Each round folder holds one returned Stata check: the do-file Stata ran, its
log (gzipped), and the R scripts that fit the same models and compare them
with the log. `data/` holds the data files Stata used. `run_parity.R` reruns
the R side with the installed suest, compares it with the stored log, and
fails if any gap between R and Stata is more than twice its size when the
round was stored, plus 1e-6. No Stata run is needed.

```sh
Rscript tools/stata-parity/run_parity.R                   # every round, about 6 minutes
Rscript tools/stata-parity/run_parity.R round09 round18   # some rounds
Rscript tools/stata-parity/run_parity.R --update round09  # store the current gaps
```

Run it from the package folder after installing the package. The `twopart`
round downloads `meps.dta` from tdmize.github.io unless `SUEST_DATA` names a
folder that has it. The GitHub workflow `stata-parity.yaml` runs every round
when `R/` or this folder changes.

| Folder | Stata check | What it covers |
|---|---|---|
| `round01` | October 3, 2026 | Joint coefficients and covariances: complementary log-log, Gamma, fractional, zero-inflated, tobit, interval, survival (exponential, lognormal, log-logistic), truncated, beta, heteroskedastic probit, and different-sample pairs |
| `round05-07` | October 4-5 | Multiple imputation; linear random-intercept models (`mixed`); survey systems; `gologit2`; Heckman |
| `round08` | October 5 | Random-intercept binary, Poisson, and ordered models with 25-point adaptive quadrature; linear random slopes |
| `round09` | November 6 | Survey models through `svy:` (ordered, multinomial, negative binomial, generalized ordered, tobit, interval, Weibull, instrumental variables, a domain) |
| `round10` | November 7 | Linear models with several random slopes and a nested level; survey complementary log-log |
| `round14-17` | November 8 | Constant-dispersion negative binomial; population-averaged (GEE) systems; `gologit2` with cauchit and log-log links |
| `round18` | November 8 | Gompertz and generalized gamma survival models; bivariate probit with different predictors |
| `round19` | November 8 | Median survival time predictions after Weibull, exponential, lognormal, and log-logistic models |
| `zeroinfl` | November 8 | Prediction types of zero-inflated Poisson and negative binomial models |
| `twopart` | November 23 | Two-part models: first and second parts, log outcomes, weights, survey design |
| `round21` | November 24 | Predictions, marginal effects, and `mecompare` changes for 38 model types; `ivprobit` probabilities from round 22 |

Each `replay.R` says what its round leaves out, such as cases where R and
Stata fit different models (Laplace against adaptive quadrature, for
example) and differences that are Stata's (see `tools/STATA-PARITY.md`).
Earlier checks (before October 2026) are numerical tests in
`tests/testthat/` and `tools/acceptance-tests/`; their do-files are in
`tools/stata-benchmarks/`.

To add a round: make a folder with the do-file, its log gzipped, the R side,
and a `replay.R` that writes `gaps.csv` with `write_gaps()` (in `helpers.R`),
put its data in `data/`, and run `run_parity.R` on the new folder once to
store its gaps.
