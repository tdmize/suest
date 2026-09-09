# Stata pweight benchmarks

These do-files generate the probability-weight benchmark data and establish the
Stata reference results used by the R acceptance tests.

1. `suest_pweight_benchmark.do` generates
   `suest_pweight_benchmark_data.dta` and demonstrates that Stata `suest`
   rejects component models fitted directly with `pweight`.
2. `suest_iweight_pweight_reference.do` refits the same weighted estimating
   equations using `iweight`, as required by Stata `suest`, and logs the joint
   robust covariance matrices for linear models.
3. `suest_pweight_binary_benchmark.do` records logit, probit, logit-probit,
   linear-logit, linear-probit, partial-overlap, and disjoint-sample reference
   covariance matrices.

The compressed CSV in `tools/acceptance-tests/data/` contains the generated
benchmark observations used by the automated R tests. The acceptance tests
compare the R coefficients and covariance elements to the high-precision Stata
output.

4. `suest_pweight_poisson_benchmark.do` generates a reproducible count
   outcome and records identical, partial-overlap, disjoint, and weight-rule
   Poisson benchmarks.

5. `suest_r_pweight_extended_benchmark.do` records the negative-binomial
   two-model system and an ordered-logit, ordered-probit, multinomial
   three-model system. The 2026-09-06 Stata 19.5 log is encoded as numerical
   references in `test_pweights_extended.R`.
6. `suest_r_cluster_benchmark.do` compares official `suest` and `suest2` for
   identical, partially overlapping, and disjoint rows with shared clusters.
7. `suest_r_iv_ancillary_benchmark.do` covers 2SLS, the linear-model ancillary
   parameter, three-model syntax, and heteroskedastic probit. Store
   conventional IV and heteroskedastic-probit estimates; request robust VCE
   from `suest2`.
8. `suest_r_cluster_crosslang_benchmark.do` uses R-generated data to verify
   identical, partially overlapping, and disjoint clustered linear-logit
   systems, clustered 2SLS, and clustered pweighted linear-logit results.
9. `suest_r_pweight_lnvar_scale_diagnostic.do` isolates whether the
   iweight-reference `lnvar` changes when otherwise equivalent pweights are
   multiplied by a common constant.
10. `suest_r_panel_fe_crosslang_benchmark.do` checks the R `plm` fixed-effects
    adapter against balanced, partially overlapping unbalanced, and
    higher-level-clustered `xtreg, fe` systems.
11. `suest_r_panel_be_re_crosslang_benchmark.do` checks the R `plm` between and
    Swamy-Arora random-effects adapters against balanced, unbalanced, default
    panel-clustered, and higher-level-clustered `xtreg` systems.
12. `suest_r_panel_ml_crosslang_benchmark.do` checks the R `nlme` random-
    intercept ML adapter, including `sigma_u`, `sigma_e`, default panel
    clusters, unequal samples, and higher-level clusters.
13. `suest_r_endogenous_crosslang_benchmark.do` checks the R bivariate-probit
    and maximum-likelihood IV-probit adapters, including their ancillary
    correlation and scale equations.
14. `suest_r_gee_crosslang_benchmark.do` checks Gaussian, binary
    logit/probit/cloglog, and Poisson PA systems under independence/exchangeable
    correlation, plus a mixed-family unequal-sample higher-cluster system.
    Generate its CSV with `make_gee_crosslang_data.R`; the matching R fits are
    in `run_gee_crosslang_reference.R`.
15. `suest_r_followup_diagnostic.do` is the final focused diagnostic for the
    returned panel/GEE benchmarks. It isolates native GEE score generation on
    full versus restricted data, prints native `xtreg, mle` covariance, and
    tests unbalanced random effects without the time-indicator design that is
    singular in `plm`.

From the repository root, generate the two CSV inputs with:

```sh
Rscript tools/stata-benchmarks/make_extended_pweight_data.R \
  tools/acceptance-tests/data/suest_pweight_benchmark.csv.gz \
  suest_r_pweight_extended_benchmark.csv
Rscript tools/stata-benchmarks/make_iv_ancillary_data.R \
  suest_r_iv_ancillary_benchmark.csv
Rscript tools/stata-benchmarks/make_cluster_crosslang_data.R \
  suest_r_cluster_crosslang_benchmark.csv
Rscript tools/stata-benchmarks/make_panel_fe_crosslang_data.R \
  suest_r_panel_fe_crosslang_benchmark.csv
Rscript tools/stata-benchmarks/make_endogenous_crosslang_data.R \
  suest_r_endogenous_crosslang_benchmark.csv
```

Run each `.do` file with its generated CSV in Stata's current working
directory. A successful run ends with its `*_BENCHMARK_COMPLETE=1` marker.
