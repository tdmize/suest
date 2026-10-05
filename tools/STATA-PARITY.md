# Stata `suest2` parity tracker

Audit date: 2026-09-10. This is a development tracker, not a published support
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
| maximum-likelihood `heckman` | `sampleSelection::selection(method = "ml")`, implemented October 5, 2026 | Native observation scores and analytic Hessian, transformed to Stata's `lnsigma` and `athrho`; scores and information match an independently coded likelihood (R tests). Predictions are the outcome equation's linear prediction (Stata's default). Stata check (`test_suest_r_stata_v07.do` Part 4, `heckman` + `regress` through `suest2`): joint coefficients within `1.1e-7`, standard errors within `1.8e-7` relative, covariance within `3.6e-7`. Passed. |
| `gologit2` | `VGAM::vglm(cumulative())`, implemented October 5, 2026 | suest's own analytic observation scores through VGAM's constraint matrices (unconstrained, proportional, and partial proportional odds; logit, probit, cloglog links) and observed information by differentiating them (VGAM's own covariance is expected information). Scores and information match an independently coded likelihood; a proportional-odds fit reproduces the `MASS::polr()` system. As in gologit2, unconstrained fits can give negative predicted probabilities away from the data. Stata check (`test_suest_r_stata_v07.do` Part 3, `suest2`): unconstrained + partial proportional odds, and partial + `logit`, joint coefficients within `9.0e-8`, standard errors within `6.0e-8` relative, covariance within `1.2e-7`. Passed. |
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
| `xtlogit, re` | Verified narrow route through `pglm::pglm()` | Binary logit with one individual random intercept, unweighted ML, and exactly 12-point nonadaptive quadrature. Fixed effects plus natural-scale `sigma`, panel-level scores, native-information bread, default/higher clustering, overlap/disjoint samples, integrated response predictions, and `marginaleffects` pass. Stata's nonadaptive coefficients/log likelihoods match below `9e-8`; because `suest2` requires adaptive quadrature, its joint covariance differs slightly as documented below. |
| `xtprobit, re` | Verified narrow route through `pglm::pglm()` | Same restricted contract as panel logit, with integrated normal-probability response predictions. Nonadaptive component coefficients and covariance match Stata within `9.29e-7` and `6.89e-8`. Because `suest2` requires adaptive quadrature, joint covariance differs slightly: fixed-effect elements by at most `1.33e-4` and all elements by at most `0.00117`. |
| `xtpoisson, re` | Component likelihood verified through `pglm::pglm()`; joint Stata covariance difference documented | Gamma random effects with individual panels, log link, unweighted ML, and `other = "sd"`. Fixed effects plus natural-scale gamma variance `alpha`, analytic panel scores, native-information bread, panel/higher clustering, overlap/disjoint samples, expected-count predictions, and `marginaleffects` pass. After converting `alpha` to `lnalpha`, component coefficients, native covariance, log likelihoods, predictions, and slopes match Stata within `1.44e-7`, `8.05e-9`, `5.80e-7`, `1.35e-7`, and `4.71e-8`. `suest2` 1.0.0 repeats the full ancillary cluster score on every observation, inflating `lnalpha` rows/columns and propagating through the bread; its higher-cluster route also integrates at the higher cluster rather than the panel. R preserves the exact panel likelihood instead. |
| `xtlogit, fe` | Deferred | `pglm(model = "within")` fails in its binomial likelihood implementation. `survival::clogit(method = "exact")` does not expose exact score residuals and its native `marginaleffects` support is Cox-style link/risk rather than the probability-scale contract needed here. No currently audited engine satisfies both the score and prediction requirements. |
| `xtpoisson, fe` | Implemented October 5, 2026 through `fixest::fepois()` with one absorbed fixed effect | Concentrated and conditional Poisson likelihoods give the same coefficients, unit scores, and information. Clustered on the fixed effect with no `G/(G-1)` factor, as native `xtpoisson, fe vce(robust)` (fixest's clustered covariance with no small-sample adjustment). Response predictions use `alpha_i(b) = log(sum y_i) - log(sum exp(x_i b))`. Stata check (`test_suest_r_stata_v07.do` Part 2): coefficients within `1.1e-9`; each model's block equals native `vce(robust)`. `suest2` 1.1.0 multiplies each block by the union-cluster `G/(G-1)` (here 116/115, variances 0.87% larger); reported to the suest2 project (`suest2_mcaghermite_xtpoisson_handoff_v01.md`). |
| `xtcloglog`, `xtologit`, `xtoprobit`, `xtmlogit`, `xtnbreg` | Pending | Conditional likelihood, integrated random-effects scores, and panel-level alignment require family-specific verification. |
| `melogit` | Verified narrow route through `glmmTMB::glmmTMB()` | Binary logit with one Gaussian random intercept, unweighted Laplace ML, fixed effects plus log SD, full group scores/covariance, default and higher clustering, overlap/disjoint samples, integrated population predictions, and `marginaleffects` pass. Laplace component coefficients and covariance match Stata within `1.22e-5` and `6.61e-7`. Joint Laplace `gsem` coefficients and covariance agree within `2.14e-5` and `0.00113`; the R and adaptive `suest2` disjoint cross-blocks are exactly zero. `suest2` rejects Laplace fits, so its 12-point adaptive results are retained as a separate approximation comparison. |
| `mepoisson` | Verified narrow route through `glmmTMB::glmmTMB()` | Poisson log with one Gaussian random intercept, unweighted Laplace ML, fixed effects plus log SD, full group scores/covariance, default and higher clustering, overlap/disjoint samples, exact integrated population means, and `marginaleffects` pass. Laplace component coefficients, covariance, and log likelihoods match Stata within `5.41e-6`, `1.31e-7`, and `1.83e-7`; joint Laplace `gsem` coefficients and covariance agree within `7.35e-6` and `0.00110`. `suest2` rejects Laplace fits, so its 12-point adaptive results are retained as a separate approximation comparison. |
| Ordinary + multilevel/panel systems (suest2 1.1.0) | Implemented October 5, 2026 for `lme`, `pglm`, and `glmmTMB` models with `lm`, logit, probit, cloglog, Poisson, NB, ordered logit/probit, and for different random-effects families together | Ordinary scores are summed within the shared group; each model's block uses its own `G_i/(G_i-1)` (and `(N-1)/(N-k)` for linear models), as suest2's heterogeneous route does. R tests: ordinary blocks equal `sandwich::vcovCL()`, multilevel blocks equal their same-type systems, cross blocks equal a by-hand assembly. Stata check `test_suest_r_stata_v05.do`: ordinary blocks (regress, logit, logit on a subsample) match exactly. `lme` + `regress` versus `mixed` + `regress` differs by up to 0.16% in standard errors because `lme` follows `xtreg, mle` (full observed-information bread) while `mixed` reports a block-diagonal covariance; `lmer` + `logit` versus `mixed` + `logit` (v06) agrees within `1.7e-7`. `logit` + random-intercept logit against `gsem` with Laplace: coefficients within `5.3e-6`, standard errors within 0.1%. `suest2` requires adaptive quadrature, so its 12-point results differ from R's Laplace fits by up to 8%. Exact check with 25-point adaptive quadrature in both programs (`test_suest_r_stata_v08.do` Part 3, `glmer(nAGQ = 25)`): `logit` + random-intercept logit and `regress` + `logit` + random-intercept logit match `suest2` within `1.1e-7` (coefficients and standard errors) and `2.2e-7` (covariance). Passed. |
| `meprobit`, `mecloglog`, `mixed`, `meglm, family(gamma) link(log)` | Implemented October 5, 2026 through `glmmTMB::glmmTMB()` (Laplace), random intercept or one correlated slope | Generic `glmmTMB` group scores and full covariance; each family's scores match an independently coded Laplace likelihood (R tests). Gaussian adds `log_sigma_e`; Gamma adds `log_shape` (Stata's `/logs` equals `-log_shape/2`). Gaussian models use `mixed`'s covariance layout: fixed-effect bread `(X'V^-1 X)^-1`, variance-parameter bread from the full observed-information inverse, no cross block (Stata's `mixed` `e(V)` is block-diagonal this way; reproduced to `3e-7`). Population-averaged predictions checked against numerical integration. Stata check `test_suest_r_stata_v06.do`: Gaussian pair against `suest2` with `mixed`, standard errors within `3.1e-5` (`glmmTMB`) and `2.9e-6` (`lmer`), passed. Probit and cloglog pairs against `gsem` with Laplace: coefficients within `2.2e-5`, but `gsem`'s robust standard errors differ by up to 0.55% (probit) and 3.6% (cloglog) from the exact Laplace calculation, whose scores match an independent likelihood; model-based diagnostic in `test_suest_r_stata_v08.do` Part 8. `suest2` with 12-point adaptive quadrature differs from R's Laplace fits by 1.3% to 8% (approximation, as expected). `test_suest_r_stata_v08.do`: with 25-point adaptive quadrature in both programs (`glmer(nAGQ = 25)`), logit, probit, cloglog, probit + random-intercept probit, and Poisson pairs match `suest2` within `9.2e-7` (coefficients), `2.3e-6` (standard errors), and `4.6e-6` (covariance); `glmmTMB` random slopes match `suest2` with `mixed, covariance(unstructured)` within `8.3e-6`. Passed. Gamma: Stata's `meglm, intmethod(laplace)` and `gsem` Laplace fits do not converge (v06, v08), so there is no Stata Laplace reference; an independent 25-point quadrature fit of the same pair in R matches `suest2` with 12-point quadrature within `1.1e-6` (standard errors) and `2.2e-6` (covariance), confirming the parameter mapping and joint covariance, and that the 1.3% gap is the Laplace approximation (`glmmTMB` offers only Laplace). `gsem`'s Laplace information differs from the exact Laplace information (model-based cloglog standard errors up to 1.8% apart, v08 Part 8), while native `melogit` and `mepoisson` Laplace fits match exactly; use the `me` commands, not `gsem`, for Laplace comparisons. |
| lme4 `lmer`, `glmer` (R-only engine; Stata counterparts `mixed`, `melogit`, `meprobit`, `mecloglog`, `mepoisson`) | Implemented October 5, 2026 | suest's own per-group likelihood: exact for `lmer`; Laplace or adaptive Gauss-Hermite with the fit's `nAGQ` for `glmer`. Reproduces lme4's log likelihood (with `nAGQ > 1`, lme4 reports it less the saturated log likelihood; for its default Laplace fits, the early-stopped inner mode search shifts it by about 1e-6 relative). `lmer` scores match merDeriv's analytic scores to 1e-6. `glmer(nAGQ = 7)` corresponds to `intmethod(mcaghermite) intpoints(7)`. Stata check `test_suest_r_stata_v06.do` Part 2: `lmer` pair against `mixed` within `2.9e-6` (standard errors), `lmer` + `logit` within `1.7e-7`, passed; `suest2` 1.1.0 refuses `mcaghermite` fits (reported to the suest2 project). `test_suest_r_stata_v08.do`: `glmer(nAGQ = 25)` against `intmethod(mvaghermite) intpoints(25)` matches within `2.3e-6` (standard errors) for logit, probit, cloglog, and Poisson pairs; `lmer` random slopes against `mixed, covariance(unstructured)` match within `1.6e-6` when `lmer` is fit with a tight optimizer tolerance (`lmerControl(optimizer = "bobyqa", optCtrl = list(rhoend = 1e-12))`); lme4's default tolerance left one model's variance parameters about `1e-4` from the optimum. Passed. |
| `meologit`, `meoprobit`, `xtologit`, `xtoprobit` | `ordinal::clmm()` random intercept, implemented October 5, 2026 | suest's own per-group likelihood (Laplace or adaptive quadrature as fitted) reproduces clmm's log likelihood (exactly with `nAGQ > 1`; clmm's Laplace stops its inner mode search at `gradTol = 1e-4`). A two-category fit reproduces the matching `glmer()` system. Stata check `test_suest_r_stata_v07.do` Part 5 did not run: `suest2` 1.1.0 refuses `intmethod(mcaghermite)` (reported to the suest2 project). Rerun with 25-point adaptive quadrature in both programs (`clmm(nAGQ = 25)`, `intmethod(mvaghermite) intpoints(25)`, `test_suest_r_stata_v08.do` Part 5): `meologit` + `meoprobit` and `meologit` + `logit` match `suest2` within `1.3e-5` (coefficients), `6.3e-6` (standard errors), and `1.3e-5` (covariance). Passed. |
| `menbreg, dispersion(constant)`, `mestreg`, multiple random slopes, three-level models | Pending | |
| linearized `svy:` | Gaussian and binary logit passed; binary probit has a documented bread-convention difference | One-stage common-design `svyglm`, explicit IDs, strata/FPCs, coefficient covariance only. The certified Gaussian route matches six stacked Stata systems to machine precision. The binary gate confirms full logit coefficient/covariance parity. Probit coefficients agree, but `survey::svyglm()` uses expected/Fisher information while the returned Stata survey-probit covariance is reproduced by an observed-information bread, producing finite-sample covariance differences up to about 7.9% on the returned benchmark diagonals. R preserves the exact native `svyglm()` covariance by contract. Different-subpopulation `suest2` systems return code 322. See `SURVEY-BENCHMARK-FINAL-RESULTS-20260908.md` and `SURVEY-BINARY-BENCHMARK-RESULTS-20260908.md`. Since October 5, 2026: linear, binary (logit, probit, cloglog), and `quasipoisson()` models combine in one system; each native covariance is reproduced and the joint matrix matches an independent PSU-total calculation (R tests). Stata check (`test_suest_r_stata_v07.do` Part 1, `suest2`): `svy:` regress + logit + Poisson joint standard errors within `2.6e-8` relative and covariance within `7.1e-8`; a logit on a domain (`subpop()`) matches its native block within `4.6e-9`. `svy:` cloglog coefficients agree (`6.6e-8`) but standard errors differ by up to 10%, the same expected-versus-observed information contract as probit. |
| `mi estimate, post:` | Initial coefficient-level support through `suest_mi()` | Compatible per-imputation SUEST systems are pooled with Rubin's rules over the complete joint coefficient vector, preserving within- and between-imputation cross-model covariance. Direct `mice::mira`, `mice::getfit()`, and `mitools::with.imputationList()` inputs are supported. Since October 5, 2026, `marginaleffects` predictions, slopes, comparisons, and cross-model differences are estimated within each imputation and pooled with Rubin's rules (as `mimrgns` does after `suest2`), on the `mira` object and on `suest_mi()` results; R tests match a by-hand pooling. Stata check `test_suest_r_stata_v05.do`: the pooled joint coefficients and covariance from `suest2` after `mi estimate, post:` match `suest_mi()` within `4.5e-8` (coefficients), `1.9e-6` (standard errors, relative), and `9.3e-6` (covariance); passed. The `mecompare` step stopped on an invalid `decimals(8)` (maximum 7). Rerun in `test_suest_r_stata_v08.do` Part 1: pooled marginal effects of `x` for both models and their difference match R (marginaleffects on the `mira` object) within `7.9e-9` (estimates) and `2.0e-6` (standard errors, relative). Passed. |
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
6. Stata `suest2` 1.0.0's specialized gamma `xtpoisson, re` route repeats each
   full ancillary cluster score on every observation. With six observations
   per panel, the benchmark's `lnalpha` meat is therefore inflated by roughly
   (6^2). When a higher cluster is requested, that route also reconstructs
   one gamma-integrated likelihood per higher cluster rather than retaining
   the original panel likelihood. The R adapter does neither: it uses one exact
   integrated-likelihood score per panel and then aggregates panels to valid
   higher clusters.
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

All returned 0.1.4 follow-up logs and the nonlinear-panel and multilevel
logit, probit, and Poisson benchmarks have been reviewed. No further Stata run
is needed for the 0.1.5 panel/multilevel increment.

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
| Random-effects panel logit, nonadaptive component fits | 6.95e-8 | 5.35e-9 native model covariance |
| Random-effects panel logit, adaptive Stata `suest2` comparison | 1.06e-3 | 3.10e-3; fixed-effect elements at most 1.17e-4 |
| Random-effects panel probit, nonadaptive component fits | 9.29e-7 | 6.89e-8 native model covariance |
| Random-effects panel probit, adaptive Stata `suest2` comparison | 1.13e-3 | 1.17e-3; fixed-effect elements at most 1.33e-4 |
| `glmmTMB`/`melogit`, Laplace component fits | 1.22e-5 | 6.61e-7 native model covariance |
| `glmmTMB`/`melogit`, joint Laplace `gsem` | 2.14e-5 | 1.13e-3 |
| `glmmTMB`/`mepoisson`, Laplace component fits | 5.41e-6 | 1.31e-7 native model covariance |
| `glmmTMB`/`mepoisson`, joint Laplace `gsem` | 7.35e-6 | 1.10e-3 |

The categorical differences are consistent with optimizer stopping
tolerances. These Stata results are now fixed numerical references in the R
acceptance suite.

## Returned parity run, October 3, 2026 (StataNow/MP 19.5, suest2 1.1.0)

`test_suest_r_stata_v01.do` stored suest2's joint e(b)/e(V) for 21 two-model
systems on a deterministic 800-row fixture (R side and comparison scripts:
`parity_reference_v01.R`, `compare_parity_v02.R`, kept in the project
records). Ancillary parameters were mapped by their transformation (sign
change for zinb `lnalpha` and Weibull `ln_p`, `exp(2x)` for tobit
`var(e.y)`, log for betareg `scale`). Gaps are relative to
sqrt(V_ii V_jj).

| System | Largest covariance gap |
|---|---:|
| fracreg logit, zip, zinb, tobit, intreg, truncreg, hetprobit | 1e-4 or less |
| streg exponential, lognormal, loglogistic (AFT) | 4e-7 or less |
| regress and logit, partial overlap and disjoint, unclustered | 3e-7 or less (regress 2e-14) |
| tobit, partial overlap | 2e-9 |
| cloglog, cloglog partial overlap, gamma log, fracreg probit, betareg | 6e-5 or less after switching these routes to observed information (before: 1.1% to 10%) |
| streg Weibull (AFT, `time`) | 0.89: a suest2 defect, not R. Native `streg, vce(robust)`, official `suest` on conventional fits, and a by-hand sandwich from `streg`'s own scores all give R's standard errors. suest2 uses `e(V_modelbased)` as the bread, and after `streg, distribution(weibull) [time] vce(robust)` that matrix is in a different parameterization (its non-constant entries equal the PH-metric model-based covariance); using it reproduces suest2's output exactly. Only Weibull is affected (PH and `time`); every other `streg` distribution stores a correct `e(V_modelbased)`. Handed to the suest2 project as `suest2_streg_weibull_handoff_v01.md` |

The partial-overlap and disjoint ordinary systems confirm the combined-sample
N/(N-1) correction. The 31-pair ME helper gate passed 30 of 31 after the R
references were aligned with the Stata estimands (centered continuous
change for Total ME of `z`; observed-information covariance for the single
plain probit); the remaining subgroup Total ME estimate differs by 4.4e-5,
mlogit's default convergence criterion: with tight convergence Stata's log likelihood and both subgroup estimates equal R's (October 4).
