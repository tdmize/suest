# Stata `suest2` parity tracker

Audit date: 2026-09-07. This is a development tracker, not a published support
claim. The Stata reference is the attached 1.0.0 `suest2.ado` and help file.
"Implemented" means the R adapter has unit tests, numerical covariance tests,
and an end-to-end `marginaleffects` test unless a narrower qualification is
shown.

## Ordinary and specialized single-level models

| Stata `suest2` family | R route | Status and qualification |
|---|---|---|
| `regress` | `stats::lm()` | Implemented for mean parameters, offsets, overlap, pweights, clustering, and `marginaleffects`. Unweighted systems include Stata's ancillary `lnvar = log(RSS / df.residual)`. Pweighted systems reproduce `suest2`'s iweight-reference `lnvar`, scores, bread, and cross-equation covariance. |
| `logit`, `logistic`, `probit`, `cloglog` | `stats::glm()`, `glm2::glm2()` | Implemented. Probit uses observed information for Stata parity. |
| `ologit`, `oprobit` | `MASS::polr()`, restricted `ordinal::clm()` | Implemented, including thresholds and pweights. `clm` scale/nominal formulas and non-flexible thresholds are rejected. |
| `mlogit` | `nnet::multinom()` | Implemented, including analytic scores, pweights, and category probabilities. |
| `poisson` | `stats::glm()`, `glm2::glm2()` | Implemented, including exposure/offset and pweights. |
| `nbreg` | `MASS::glm.nb()` | Implemented with `log(theta)`, including pweights. |
| `zip`, `zinb` | `pscl::zeroinfl()` | Implemented for logit/probit inflation links and response-scale expected counts. The ZINB adapter restores the dispersion score and observed-information block omitted by the package's default sandwich methods. |
| `glm` | `stats::glm()`, `glm2::glm2()` | Implemented for identity, log, logit, probit, cloglog, and user-supplied loglog links. |
| `fracreg logit`, `fracreg probit` | quasi-binomial `glm()`/`glm2()` | Implemented. R additionally supports fractional cloglog and user-supplied loglog links. |
| `tobit` | `censReg::censReg()` or Gaussian `survival::survreg()` | Implemented with scale and latent-mean predictions. The two R likelihoods agree in coefficients, scores, and bread. |
| `intreg` | Gaussian interval-censored `survival::survreg()` | Implemented and tested with `marginaleffects`. |
| parametric `streg` | `survival::survreg()` | Partially implemented. Weibull, exponential/fixed-scale, lognormal, and loglogistic AFT forms are feasible; Stata-only parameterizations/distributions still need mapping. |
| `betareg` | `betareg::betareg()` | Implemented for standard beta regression with logit, probit, cloglog, and loglog mean links, including modeled precision. |
| `truncreg` | `truncreg::truncreg()` | Implemented for left/right truncated Gaussian models with scale. |
| `ivregress 2sls` | unweighted IV `fixest::feols()` | Implemented without absorbed fixed effects. Coefficients, IV scores, bread, overlap, and `marginaleffects` are tested. The 2026-09-06 Stata benchmark confirms the identical- and partial-sample covariance, including Stata's route-specific HC0 convention. |
| maximum-likelihood `heckman` | likely `sampleSelection::selection()` | Blocked in this runtime: the current dependency chain requires `nloptr`, whose source build requires unavailable CMake. No support claim. |
| `gologit2` | candidate `VGAM::vglm(cumulative())` | Not admitted. VGAM exposes bread but no observation-level `estfun`; unrestricted fits can also produce crossing/negative category probabilities. A separately verified score adapter would be required. |
| `hetprobit` | `Rchoice::hetprob(link = "probit")` | Implemented with the full location/log-scale score system. R also supports the package's heteroskedastic logit model. |
| `biprobit` | `mvProbit::mvProbit()` | Verified for the tested bivariate models using the same regressors in both equations. Recomputed curvature in athrho coordinates resolves a fitting-engine Hessian precision issue; full joint covariance matches Stata to 8.67e-11. Joint-success prediction and marginaleffects are tested. |
| `ivprobit` | `Rchoice::ivpml()` | Verified for the tested ML models with one continuous endogenous regressor. Recomputed curvature resolves an inaccurate Rchoice Hessian. With reltol=1e-12, the full structural/reduced-form covariance, lnsigma, and athrho match Stata to 3.51e-9. Marginaleffects is tested. |
| `ivtobit` | candidate not yet selected | Pending. |

## Panel, multilevel, survey, and MI routes

| Stata route | R status | Reason/next requirement |
|---|---|---|
| `xtreg` FE/RE/BE | FE/BE and balanced RE verified through `plm::plm()` | Six returned cases verify full covariance, including unequal FE/BE samples and higher clusters. The simpler unbalanced RE comparison is close but not exact: maximum coefficient difference 1.39e-4, covariance difference 5.57e-6, and relative diagonal difference 0.1504%. Swamy-Arora is the closest tested `plm` method. The time-indicator benchmark fails inside plm because its unbalanced between design is rank deficient. |
| `xtreg, mle` | Verified within native engine precision through `nlme::lme(method = "ML")` | Single-level random-intercept Gaussian models include fixed effects, `sigma_u`, and `sigma_e`; default panel clusters, higher nested clusters, likelihood scores, observed bread, and `marginaleffects` pass R tests. Returned joint systems differ by at most 7.58e-7, or 0.00808% on a covariance diagonal. Stata's native model covariance already differs from R's independently verified likelihood covariance by up to 0.00412%, so this is classified as a small fitting-engine convention/precision difference. |
| `xtreg/xtlogit/xtprobit/xtcloglog/xtpoisson, pa` | Supported subset through `geepack::geeglm()` | Unweighted numeric outcomes, independence/exchangeable correlation. Ten single-family full covariance matrices match Stata to 4.61e-10. Unequal samples, mixed families, higher clusters, offsets, interactions, and prediction delta-method SEs pass R tests. Stata `suest2` 1.0.0 fails on a mixed-family unequal-sample case because native `predict, score` sees excluded rows in the panels; the same prediction succeeds after keeping only `e(sample)`. AR1, unstructured, and negative-binomial GEE remain unsupported. |
| correlated random effects | Provisional candidate through explicit panel means plus `plm::plm(model = "random")` | Uses the provisional Swamy-Arora random-effects route. A dedicated CRE comparison, including the treatment of unequal samples when constructing panel means, is still required. |
| `xtlogit`, `xtprobit`, `xtcloglog`, `xtologit`, `xtoprobit`, `xtmlogit`, `xtpoisson`, `xtnbreg` | Pending | Conditional likelihood, integrated random-effects scores, and panel-level alignment require family-specific verification. |
| `mixed`, `melogit`, `meprobit`, `mecloglog`, `mepoisson`, `menbreg`, `meologit`, `meoprobit`, `meglm`, `mestreg` | Pending | Integrated likelihood scores must be aggregated at the common highest-level group and matched to Stata's finite-sample correction. |
| linearized `svy:` | Gaussian and binary logit passed; binary probit has a documented bread-convention difference | One-stage common-design `svyglm`, explicit IDs, strata/FPCs, coefficient covariance only. The certified Gaussian route matches six stacked Stata systems to machine precision. The binary gate confirms full logit coefficient/covariance parity. Probit coefficients agree, but `survey::svyglm()` uses expected/Fisher information while the returned Stata survey-probit covariance is reproduced by an observed-information bread, producing finite-sample covariance differences up to about 7.9% on the returned benchmark diagonals. R preserves the exact native `svyglm()` covariance by contract. Different-subpopulation `suest2` systems return code 322. See `SURVEY-BENCHMARK-FINAL-RESULTS-20260908.md` and `SURVEY-BINARY-BENCHMARK-RESULTS-20260908.md`. |
| `mi estimate, post:` | Initial coefficient-level support through `suest_mi()` | Compatible per-imputation SUEST systems are pooled with Rubin's rules over the complete joint coefficient vector, preserving within- and between-imputation cross-model covariance. Direct `mice::mira`, `mice::getfit()`, and `mitools::with.imputationList()` inputs are supported. Nonlinear predictions, slopes, and comparisons remain pending. |
| pweighted ordinary models | Implemented subset | Linear, logit/probit, Poisson, NB, ordered logit/probit, and multinomial match the documented Stata pweight family list. The 2026-09-06 Stata 19.5 benchmark confirms the extended NB and three-model categorical systems numerically. |
| joint `cluster()` covariance | Implemented and cross-language verified | R accepts cluster columns or model-aligned cluster vectors, aggregates the aligned system scores by cluster, and supports partial overlap plus disjoint observations in shared clusters. Three ordinary systems confirm Stata's system-level `G/(G-1)` adjustment; clustered 2SLS confirms the specialized raw grouped-influence route. The pweighted clustered benchmark also verifies linear `lnvar` and mixed linear-logit blocks. |

## Confirmed R-side package issues or limitations

1. `sandwich::estfun.polr()` ignores model offsets. The R adapter therefore
   uses its own ordered-model score calculation.
2. `pscl`'s default ZINB `estfun()`/`vcov()` omit the dispersion parameter.
   The R adapter supplies the complete score and observed-information blocks.
3. `survival::survreg()` with stratum-specific scales exposes a score layout
   that does not match its fitted covariance parameters. These models are
   rejected rather than silently misaligned.
4. `censReg` has no native `predict()` method. The combined SUEST object
   supplies latent-mean prediction directly for `marginaleffects`.
5. `Rchoice::predict.hetprob()` fails when the documented default probit link
   is omitted from the fit call. The combined object restores the default on
   its internal copy before prediction.
6. The provisional IV-probit adapter previously selected Rchoice's
   residual-conditioned response probability despite documenting the
   structural probability. It now evaluates `pnorm(X beta)`; independent
   formula, instrument-invariance, and analytic-derivative tests cover this.
7. The provisional `nlme` ML adapter previously read a rounded standard
   deviation from `VarCorr()` display output. It now extracts the stored
   covariance at full precision. A single random slope is also explicitly
   rejected rather than mistaken for a random intercept.
8. `geeglm` inherits from `glm`, so generic GLM dispatch is not sufficient:
   GEE requires the fitted working-correlation adjustment. The explicit GEE
   adapter now intercepts these objects and refuses unsupported specifications.
9. Rchoice's IV-probit analytic Hessian disagrees with derivatives of its
   otherwise accurate scores for an overidentified fit. Re-evaluating its
   own Hessian function reproduces the discrepancy. The R adapter now computes
   observed information from central differences of analytic score sums.
   The pre-fix full covariance gap was 0.000671724 (largest relative diagonal
   gap 2.2146%); corrected curvature and tighter fitting reduce it to 3.51e-9.
10. The mvProbit final Hessian also has a smaller precision discrepancy.
    Independent bivariate-normal score derivatives reduce the full covariance
    gap from 2.55e-7 to 8.67e-11. The adapter now uses this curvature.
11. plm's unbalanced Swamy-Arora fit fails for the time-indicator benchmark:
    its between model has seven design columns but rank four. It tries to
    invert the singular cross-product before suest is called. Do not silently
    switch to another variance-component estimator to claim Stata parity.

## Current verification checkpoint

The GEE subset has 81 focused expectations against geepack 1.3.13, including
native naive/robust covariance, per-model finite-sample correction, identical
and partial samples, higher clusters, offsets/interactions, and independent
prediction delta-method standard errors. The deterministic R reference script
runs all 11 cases in the matching Stata benchmark. Ten returned single-family
cases verify cross-language equality. The eleventh fails in Stata's native
score extraction, while the corresponding R mixed-family system passes.

All returned follow-up logs have been reviewed. No further Stata run is needed
for the 0.1.4 closeout.

## Stata documentation discrepancies and items to verify

1. The 1.0.0 help syntax displays exactly two model names, but a 2026-09-06
   Stata 19.5 run confirms that `suest2` accepts three stored models. The
   runtime supports this case; the help syntax is stale.
2. `CHANGELOG-suest2.md` labels the main ado's "Current banner" as 0.1.93,
   while the attached `suest2.ado` and help file are version 1.0.0.
3. The help description spells "seamlessly" as "seemlessly".
4. The attached 1.0.0 ado actively scans, dispatches, and implements an
   `xtreg, cre` route, while the attached 1.0.0 help says that `xtreg, cre` is
   not supported. The code and documentation therefore disagree even though
   the same help recommends the explicit-Mundlak equivalent.
5. Confirmed ancillary scale sensitivity: the returned pweight diagnostic
   shows that multiplying weights by 10 changes suest2's lnvar by
   -0.00222088538830051 and its variance by 7.04424779679382e-6. Original
   pweight-fit lnvar and the mean coefficients/covariance remain invariant.
   The suest2 result agrees with log(weighted RSS/(sum(weights)-k)). Possible
   remedy: normalize weights in the ancillary reconstruction or use a
   scale-invariant pweight variance convention, with corresponding score and
   bread changes. The R implementation currently retains Stata parity; this
   behavior is documented, not silently changed.
6. Confirmed Stata execution failure: the mixed Gaussian/logit PA benchmark
   ends at suest2's native logit score extraction with r(134), after both
   component models fit successfully. Direct `predict, score` fails on the full
   data because 20 panels contain four rows while the stored fit used at most
   three; it succeeds after keeping only `e(sample)`. A possible `suest2` fix is
   to isolate the estimation rows around native PA score prediction and then
   merge the score back by the preserved observation identifier.
7. ML is closed as a fitting-engine precision difference: coefficients agree
   to 4.05e-9 and joint covariance differs by at most 7.58e-7 (0.00808% on a
   diagonal). Independent Gaussian likelihood checks support the R scores and
   information; Stata's native model covariance itself differs from R by up to
   0.00412% on a diagonal.

## Completed cross-language benchmarks

The 2026-09-06 Stata 19.5 extended-pweight log completed both systems. After
transforming Stata's `ln(alpha)` to R's `ln(theta) = -ln(alpha)`, the largest
absolute R/Stata differences were:

| System | Coefficients | Covariance elements |
|---|---:|---:|
| Negative binomial, two models | 5.62e-8 | 1.64e-9 |
| Ordered logit + ordered probit + multinomial | 5.19e-5 | 3.12e-7 |
| Instrumental-variables 2SLS, two models | 3.69e-9 | 1.83e-11 after removing the inappropriate HC1 correction |
| Linear ancillary `lnvar`, two models | 4.52e-9 | 2.54e-11 |
| Clustered ordinary models, three sample designs | 2.06e-8 | 1.45e-9 |
| Clustered 2SLS, two models | 4.73e-9 | 1.25e-11 |
| Clustered pweighted linear + logit, including `lnvar` | 1.73e-8 | 9.01e-11 |
| FE: balanced, unequal samples, higher clusters | 4.01e-9 | 6.07e-10 |
| BE: balanced and unbalanced/higher clusters | 8.57e-9 | 1.34e-9 |
| Balanced RE | 3.52e-9 | 2.16e-10 |
| Unbalanced RE, no time indicators | 1.39e-4 | 5.57e-6 (0.1504% maximum relative diagonal difference) |
| Random-intercept panel ML | 4.05e-9 | 7.58e-7 (0.00808% maximum relative diagonal difference) |
| GEE: ten single-family cases | 5.55e-8 | 4.61e-10 |
| Bivariate probit, after curvature correction | 1.08e-8 | 8.67e-11 |
| IV probit, after curvature correction and tighter R fit | 3.78e-7 | 3.51e-9 |

The categorical differences are consistent with optimizer stopping
tolerances. These Stata results are now fixed numerical references in the R
acceptance suite.
