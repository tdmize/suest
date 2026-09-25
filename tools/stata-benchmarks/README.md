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
16. `suest_r_panel_logit_re_crosslang_benchmark.do` validates the restricted
    `pglm` random-effects panel-logit route. Nonadaptive 12-point component
    fits match the R likelihood exactly; the joint cases use 12-point adaptive
    quadrature because `suest2` rejects nonadaptive fits. The returned reference
    is encoded by `build_panel_logit_re_stata_fixture.R`.
17. `suest_r_panel_probit_re_crosslang_benchmark.do` applies the same two-track
    nonadaptive component/adaptive `suest2` design to the restricted random-
    effects panel-probit route. Generate its CSV with
    `make_panel_probit_re_crosslang_data.R`; the matching R systems are produced
    by `run_panel_probit_re_reference.R`.
18. `suest_r_panel_poisson_re_crosslang_benchmark.do` checks the restricted
    gamma random-effects panel-Poisson route, including native component fits,
    expected-count predictions and slopes, default/higher clusters, partial
    overlap, and disjoint panels. Generate its CSV with
    `make_panel_poisson_re_crosslang_data.R`; the matching R systems are
    produced by `run_panel_poisson_re_reference.R`. The returned Stata 19.5
    reference is encoded by `build_panel_poisson_re_stata_fixture.R`; see
    `PANEL-POISSON-RE-DESIGN-20260910.md` for the documented ancillary-score
    and higher-cluster differences in `suest2` 1.0.0.
19. `suest_r_glmmtmb_logit_ri_crosslang_benchmark.do` checks the restricted
    `glmmTMB` Laplace random-intercept logit route against Stata `melogit`.
    It compares Laplace component fits, marginal predictions, and slopes;
    records the expected `suest2` refusal of Laplace fits; runs balanced,
    partial-overlap, higher-cluster, and disjoint `suest2` systems with 12-point
    adaptive quadrature; and repeats the joint systems with Laplace `gsem` and
    robust covariance for the exact-integration comparison. Generate its CSV
    with `make_glmmtmb_logit_ri_crosslang_data.R`; the matching R systems are
    produced by `run_glmmtmb_logit_ri_reference.R`. The returned Stata 19.5
    reference is encoded by `build_glmmtmb_logit_ri_stata_fixture.R`; no further
    Stata run is required for this increment.
20. `suest_r_glmmtmb_poisson_ri_crosslang_benchmark.do` applies the same
    component, adaptive-`suest2`, and Laplace-`gsem` design to the restricted
    `glmmTMB` random-intercept Poisson-log route. Generate its CSV with
    `make_glmmtmb_poisson_ri_crosslang_data.R`; the matching R systems are
    produced by `run_glmmtmb_poisson_ri_reference.R`. The returned Stata 19.5
    reference is encoded by `build_glmmtmb_poisson_ri_stata_fixture.R`; no
    further Stata run is required for this increment.

The development Poisson random-slope increment has completed its Stata gate.
All six Laplace components and all three joint models converged. The v4
balanced and v5 partial/disjoint joint fits each converged in one iteration
from the returned component estimates. No further Stata run is required for
this increment; do not rerun the older cold-start scripts.

Raw GSEM covariance remains diagnostic evidence: cluster-aligned standard
errors differ from R by up to 0.9952%. Independent score/curvature audits agree
with R within 0.0013% in SE and localize the raw difference chiefly to the
returned joint model-based covariance. Combining recovered Stata joint score cross-products
with validated Stata component curvature agrees with R within 0.0009%.
Raw, cluster-aligned, and reconstructed results are kept separate. See
`../GLMMTMB-POISSON-RS-VALIDATION-20260923.md` for the complete interpretation.

From this directory, `make_glmmtmb_poisson_rs_crosslang_data.R` and
`run_glmmtmb_poisson_rs_reference.R` regenerate the data and R reference.
`build_glmmtmb_poisson_rs_stata_components.R <v2-log>` rebuilds the component
fixture; `build_glmmtmb_poisson_rs_stata_balanced.R <v4-log>` rebuilds balanced
results; `build_glmmtmb_poisson_rs_stata_remaining.R <v5-log>` rebuilds partial
and disjoint results. The two joint builders share `poisson_rs_joint_log_helpers.R`.
`compare_glmmtmb_poisson_rs_stata_joint.R [output.csv]` compares raw and
cluster-aligned covariance. `check_glmmtmb_poisson_rs_joint_audits.R [output.rds]`
runs all three independent audits, including missing-data handling, higher
cluster aggregation, disjoint correction, and coordinate sensitivity. Earlier
balanced-only audit scripts remain for historical reproduction.

The development NB2 increment uses
`suest_r_glmmtmb_nbinom2_ri_crosslang_benchmark_v2.do` and the matching CSV.
Run `make_glmmtmb_nbinom2_ri_crosslang_data_v2.R` and
`run_glmmtmb_nbinom2_ri_reference_v2.R` from this directory to regenerate the
data and R reference. Its six Laplace `menbreg` components and three joint
Laplace `gsem` systems are the primary gate; optional adaptive `suest2`
systems are a separate comparison. Revision 1's right-half sample had a
zero-variance boundary and is now a negative regression test. Revision 2 uses
odd/even disjoint groups with unchanged base outcomes. The returned Stata 19.5
revision-2 log completed all six components, all three joint systems, and all
three optional adaptive systems. No further Stata run is required for this
increment. Rebuild the full fixture with
`build_glmmtmb_nbinom2_ri_stata_fixture.R <returned-v2-log>` and compare it with
the R reference using `compare_glmmtmb_nbinom2_ri_stata.R [output.csv]`, both
from this directory. The [validation report](../GLMMTMB-NBINOM2-RI-VALIDATION-20260923.md)
records exact differences and the disjoint cluster-correction adjustment;
adaptive fits remain a separate approximation comparison.

The development NB2 random-slope route has also completed all three Stata
gates: `suest_r_glmmtmb_nbinom2_rs_balanced_v1.do`,
`suest_r_glmmtmb_nbinom2_rs_partial_v2.do`, and
`suest_r_glmmtmb_nbinom2_rs_disjoint_v3.do`. No new Stata run is required.
Run `make_glmmtmb_nbinom2_rs_crosslang_data.R` and
`run_glmmtmb_nbinom2_rs_reference.R` here to regenerate the data and R
reference. For each case (`balanced`, `partial`, `disjoint`), the corresponding
`build_glmmtmb_nbinom2_rs_stata_<case>.R <log>` parses raw output,
`compare_glmmtmb_nbinom2_rs_stata_<case>.R [output.csv]` compares it with R,
and `check_glmmtmb_nbinom2_rs_<case>_audit.R [output.rds]` independently
audits likelihood, curvature, and scores. Raw and audited fixtures remain
separate. The disjoint comparison explicitly aligns the cluster correction
by 99/98. The shared audit retains Fisher, natural-variance, and native
correlation-coordinate curvature, including the finite-gradient Hessian
chain-rule term. The disjoint `*_margins.R` and `*_curvature.R` diagnostics
check native Stata margins and the coordinate refinement. See the
[final validation report](../GLMMTMB-NBINOM2-RS-VALIDATION-20260924.md).

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
