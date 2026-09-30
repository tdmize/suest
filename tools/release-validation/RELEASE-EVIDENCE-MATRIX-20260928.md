# Model-specific release evidence — suest 0.1.6.9000

Baseline: `29ae8519ed73fa4720b1138ef93e4aa8edc3d188`. Source audit and new R
execution dated 2026-09-28; extended-family update dated 2026-09-29. Gate 1 passes all 40 stacked comparisons; its paired
cluster bootstrap clears the diagnostic screen. The gate-1 admission defects and the gate-2 interval-survival scale-score defect
are reproduced and fixed locally. Gate 2 passes all 108 independent likelihood,
full-covariance, and contrast-SE cases. Gate 3 passes 100 numerical panel/GEE/
survey cases (300 contrasts), with four explained fitting warnings and an
explicit FE average-level inference limitation. Gate 4 closes the specified
alternative-engine, bivariate-probit and IV contrast gaps and adds four
convergence safeguards. All 28 numerical cases pass, including one reviewed
fitting-warning case. No new Stata run was performed.
Historical baseline CI passes all 1,630 expectations; that result does not independently certify
every supported specification. This is a map of identified numerical evidence
and remaining questions, not a percentage confidence score or a new pass stamp.

## Reading the matrix

- **External:** a returned Stata matrix or numerical publication reference is
  asserted. Unless the cited source specifically says official `suest` or
  `gsem`, do not upgrade an existing `suest2` comparison to official-Stata evidence.
- **Independent formula/likelihood:** separately computed equations,
  derivatives, integration, or influence algebra check the named quantity.
- **Native:** agrees with the upstream engine or its sandwich methods. Useful
  for adapter fidelity, but does not rule out an upstream error shared by both sides.
- **Finite:** the workflow runs and returns finite values; magnitude remains
  uncertified by that particular test.
- **Gap:** an independent check was not identified for the named quantity in
  the reviewed tests/reports. This does not mean no tests exist or results are wrong.

Full covariance means both within-model and cross-model entries, including
nuisance parameters where applicable. An analytic contrast SE using `vcov(fit)`
checks differentiation/propagation conditional on that matrix; it must be read
together with the separate covariance evidence. New checks in gate 1 are
executed on R 4.5.2; see both release-validation records for bounds and results.

All test references below are in `tests/testthat/` unless prefixed otherwise.
Only the documented narrow contracts are considered; offsets, pweights,
clustering, and mixed-family combinations have separate restrictions.

## Ordinary and extended regression routes

| Route and engine | Component and covariance evidence | Prediction/effect evidence | Remaining release action | Primary references |
|---|---|---|---|---|
| Linear, `lm` | Returned full Stata ancillary/covariance matrices; clustered and weighted mixed systems | Paper examples assert effect differences and SEs | Preserve ancillary convention and test supported nonconverged inputs separately | `test-core.R`, `test-stata-cluster-reference.R`; `tools/acceptance-tests/tests/test_paper_examples.R` |
| Binary logit, `glm` | Native sandwich plus full clustered/pweight external references | Paper cross-model effects and SEs | RV-01 fixed; gate 1 stacked and bootstrap pass | `test-cluster.R`, `test-stata-cluster-reference.R`, paper tests |
| Binary probit, `glm` | Native components and pweight references | Existing effects workflow; not every combination has an external end-to-end reference | RV-01 fixed; gate 1 observed-information covariance and analytic contrasts pass | `test-model-families.R`, `test-pweights.R`, acceptance pweight tests |
| Binary cloglog, `glm` | Native diagonal covariance blocks | Finite effects | Gate 1 full stacked covariance and contrasts pass | `test-glm-families.R` |
| Poisson log, `glm` | Native and pweight numerical references | Existing effects workflow | Gate 1 full clustered covariance and contrasts pass | `test-pweights.R`, acceptance pweight tests |
| Other GLMs: Gaussian identity and Gamma log explicitly exercised | Native diagonal sandwich comparisons | Finite cross-model effects | Gate 1 passes for these two routes; enumerate further family/link combinations before a blanket claim | `test-glm-families.R`; generic admission in `R/adapters.R` |
| Fractional quasibinomial: logit, probit, cloglog, custom loglog | Native diagonal comparisons for all four links | Finite effects; loglog has a separate workflow case | Gate 1 passes for all four links, full joint covariance and analytic contrast SEs | `test-glm-families.R` |
| Alternative `glm2` engine: logit, probit, Poisson | Cross-engine adapter tests | Mainly finite-output checks | Gate 4 matched full-system covariance and all three contrast SEs pass across four sample patterns; engine equivalence, supported by prior independent GLM gates | `test-model-adapters.R`; `tools/acceptance-tests/tests/test_model_adapters.R` |
| Negative binomial log, `MASS::glm.nb` | Dispersion scores included; returned pweight full covariance | Paper different-outcome effect differences and SEs | Gate 4 rejects failed mean or dispersion convergence; preserve nuisance transformation | `test-model-families.R`; acceptance `test_invariants.R`, `test_pweights_extended.R`, paper tests |
| Ordered logit/probit, `MASS::polr` | Analytic scores; offset-specific tests; returned weighted categorical full matrix | Paper ordinal-versus-nominal effects/SEs; probability/effect-sum invariants | Gate 4 matches both links to clm including category-specific contrast SEs; failed fits rejected; retain offset regression | `test-offsets.R`, `test-model-families.R`; acceptance categorical and paper tests |
| Ordered logit/probit, restricted `ordinal::clm` | Native block covariance and comparison with `polr`; extensions rejected | Cross-engine effects often checked for finiteness | Gate 4 full-system named-parameter covariance and category-specific contrast SEs match polr; unreliable convergence diagnostics rejected | `test-model-adapters.R` |
| Multinomial logit, `nnet::multinom` | Analytic multinomial scores/information; weighted full-system reference | Paper ordinal-versus-nominal contrast SEs | RV-02 fixed; gate 4 also rejects native optimizer failure and verifies valid controls | `test-model-families.R`; paper and weighted categorical tests; `R/internals.R` |
| Gaussian censored/interval, `survival::survreg` | Gate 2 independent full covariance for left/right/interval and fixed-scale left censoring; mixed-censoring score regressions | Independent prediction and finite-change contrast SEs pass; latent-mean target | RV-03 interval log-scale score fixed; penalized fits explicitly rejected; conventional robust fits preserve model-based bread | `test-survreg.R`, `test-extended-release-safeguards.R`; gate 2 |
| Parametric AFT, `survival::survreg` | Gate 2 independent full covariance for Weibull/lognormal/loglogistic, right and interval censoring | Independent contrast SEs pass; response target is exp(linear predictor), not expected time | Preserve distribution/scale boundaries; no blanket claim for arbitrary custom distributions | `test-survreg.R`, extended safeguards; gate 2 |
| Beta: logit/probit/cloglog/loglog, `betareg` | Gate 2 independent expected information, observed-information variant, full covariance with precision parameters | Independent prediction and finite-change contrast SEs pass | Ordinary ML only; BR/BC now explicitly rejected; keep expected/observed convention explicit | `test-betareg.R`, extended safeguards; gate 2 |
| Zero-inflated Poisson, `pscl::zeroinfl` | Gate 2 independent mixture likelihood and full covariance, including count/inflation cross-blocks | Independent response prediction and finite-change contrast SEs pass | Preserve logit/probit inflation-link scope | `test-zero-inflated.R`; gate 2 |
| Zero-inflated NB, `pscl::zeroinfl` | Gate 2 independent mixture likelihood and full covariance, including log dispersion | Independent response prediction and finite-change contrast SEs pass | Preserve nuisance coordinates and supported inflation links | `test-zero-inflated.R`; gate 2 |
| Truncated Gaussian, `truncreg` | Gate 2 independent complete covariance for left/right truncation | Independent latent-mean prediction and finite-change contrast SEs pass | NR passes all 8 cases; default BFGS warnings and one stationarity failure are preserved separately; no formula change | `test-truncreg.R`; gate 2 optimizer finding |
| Censored Gaussian/Tobit, `censReg` | Gate 2 independent full covariance for left/right/two-limit censoring | Independent latent-mean prediction and finite-change contrast SEs pass | Preserve latent-mean target and fitted-sample definitions | `test-censreg.R`; gate 2 |
| Unweighted 2SLS, `fixest::feols` without absorbed FE | Independent IV equations and external full covariance, including clustered HC0 convention | Predictions/effects workflows | Gate 4 independent projected-X HC0 full covariance and analytic prediction/change/slope contrast SEs pass | `test-ivreg.R`, `test-cluster.R`, `test-stata-cluster-reference.R` |
| Heteroskedastic probit, `Rchoice::hetprob` | Gate 2 independently written likelihood, scores, curvature, and full covariance | Independent prediction and finite-change contrast SEs pass | Retain recorded numerical stability bounds and scale-equation coordinates | `test-hetprob.R`; gate 2 |
| Heteroskedastic logit, `Rchoice::hetprob` | Gate 2 independent full covariance for this R extension | Independent prediction and finite-change contrast SEs pass | Do not imply official Stata parity from this R likelihood check | `test-hetprob.R`; gate 2 |
| ML IV probit, `Rchoice::ivpml` | Returned full Stata system; independent likelihood scores and corrected curvature | Structural probability and analytic cross-model finite-change SE check | Preserve structural estimand, one-endogenous-variable contract, tight optimizer checks | `test-ivprobit.R`, `test-stata-endogenous-reference.R`, `test-returned-crosslang.R` |
| Bivariate probit, `mvProbit` | Returned full matrix after independently recomputed curvature; correlation transformation | Independent joint-success probability; most slope SE checks finite | Gate 4 independent likelihood full covariance and joint-success prediction/change/slope contrast SEs pass; one reviewed trial-CDF warning case; same-regressor restriction remains | `test-biprobit.R`, `test-returned-crosslang.R` |

## Panel and multilevel routes

| Route and engine | Component and covariance evidence | Prediction/effect evidence | Remaining release action | Primary references |
|---|---|---|---|---|
| Individual linear FE, `plm` | Gate 3 independent transformed covariance and slope/finite-change contrast SEs pass | Conditional prediction propagation passes away from training means; average-level inference at matching means is unsupported | Preserve/document structural zero-variance level boundary; do not interpret missing/tiny SEs as precision | `test-panel-fe.R`, `test-returned-crosslang.R`; gate 3 |
| Individual linear BE, `plm` | Gate 3 independent panel-mean covariance, including unbalanced samples | Independent prediction/slope/finite-change contrast SEs pass | Preserve panel-mean and finite-sample conventions | same; gate 3 |
| Swamy-Arora linear RE, `plm` | Gate 3 independent transformed covariance using native quasi-demeaning; prior external evidence retained | Independent prediction/slope/finite-change contrast SEs pass | Preserve native variance-estimator convention and documented unbalanced Stata differences | same; `tools/STATA-PARITY.md`; gate 3 |
| Gaussian random-intercept ML, `nlme::lme` | Independent likelihood scores/information; full Stata matrices within native precision | Gate 3 independent mean/slope/finite-change gradients and contrast SEs conditional on audited V | Preserve ML/common-variance restrictions; no new covariance certification implied | `test-panel-ml.R`, `test-returned-crosslang.R`; gate 3 |
| Random-intercept panel logit, `pglm`, 12-point nonadaptive | Independent likelihood/score checks; component Stata reference; adaptive joint comparison separately documented | Gate 3 independent quadrature gradients and all three contrast SEs pass; continuous integration difference measured separately | Preserve fitted 12-point nonadaptive approximation and prior Stata distinction | `test-panel-logit-re.R`; `tools/PANEL-LOGIT-RE-DESIGN-20260910.md`; gate 3 |
| Random-intercept panel probit, `pglm`, 12-point nonadaptive | Independent likelihood/information and external component evidence | Gate 3 independent quadrature/closed-form integrated-probit contrasts and SEs pass | Preserve approximation distinction; new propagation evidence conditions on prior-audited V | `test-panel-probit-re.R`; corresponding design report; gate 3 |
| Gamma random-effects Poisson, `pglm` | Exact panel-likelihood scores and external component agreement; known Stata suest2 score issue | Gate 3 independent expected-count/slope/change contrasts and SEs pass; alpha derivative correctly zero | Four fitting-warning cases traced to invalid native trial values; final exact-likelihood checks pass; do not match prior suest2 score bug | `test-panel-poisson-re.R`; `tools/STATA-PARITY.md`; gate 3 |
| GEE Gaussian identity, binary logit/probit/cloglog, Poisson log; independence/exchangeable | Ten prior returned full Stata matrices; gate 3 independent working-mean equations with unequal samples/higher clusters | Gate 3 independent prediction/slope/change contrast SEs pass across ten specifications | Preserve native fitted working correlation/scale; numeric/unweighted restrictions | `test-gee.R`, `test-returned-crosslang.R`; gate 3 |
| Logit Gaussian random intercept, `glmmTMB` | Independent Laplace score checks; component and joint Stata comparisons with approximation limits | Gate 3 independent normal-mixture gradients and all contrast SEs pass, including variance uncertainty | New propagation evidence conditions on prior-audited V; retain Laplace/GSEM qualifications | `test-glmmtmb-logit-ri.R`; corresponding design report; gate 3 |
| Poisson Gaussian random intercept, `glmmTMB` | Independent Laplace scores and returned component/joint Stata evidence | Gate 3 closed-form integrated mean/slope/change gradients and all contrast SEs pass | Random-effect variance derivatives materially affect SEs; retain full nuisance uncertainty | `test-glmmtmb-poisson-ri.R`; corresponding design report; gate 3 |
| NB2 Gaussian random intercept, `glmmTMB` | Independent likelihood/score checks, external components/joint matrices, boundary rejection | Analytic means, slopes, contrasts including variance uncertainty | Preserve stated curvature/correction limits; additional seeds only to address concrete risk | `test-glmmtmb-nbinom2-ri.R`; `tools/GLMMTMB-NBINOM2-RI-VALIDATION-20260923.md` |
| Poisson correlated intercept/numeric slope, `glmmTMB` | Independent Laplace and joint covariance audits for balanced/partial/disjoint cases | Analytic integrated means/effects and nuisance uncertainty | Raw GSEM curvature differences remain distinct; ordinary user-scale sensitivity warrants care | `test-glmmtmb-poisson-rs.R`, `test-glmmtmb-poisson-rs-stata.R`; validation report |
| NB2 correlated intercept/numeric slope, `glmmTMB` | Independent likelihood/score/curvature and full joint audits, including coordinate chain rule | Analytic effects and full nuisance uncertainty | Raw SE differences remain; the scale test supplies corrected native curvature and is not proof that an unmodified poorly scaled fit is safe | `test-glmmtmb-nbinom2-rs.R`, `test-glmmtmb-nbinom2-rs-stata.R`; validation report |
| Bernoulli logit correlated intercept/numeric slope, `glmmTMB` | Independent Laplace/score/full covariance audits for all three sample patterns | Independent integration and analytic prediction/slope/factor gradients and contrast SEs | Preserve up-to-4.23% raw GSEM SE difference and failed 1.052e-6 vs 1e-6 disjoint reconstruction; inspect scaling and differencing sensitivity | `test-glmmtmb-logit-rs.R`, `test-glmmtmb-logit-rs-stata.R`; balanced/partial/disjoint reports |

## Survey, pooling, and system-level features

| Route | Existing evidence | Remaining boundary/action | Primary references |
|---|---|---|---|
| One-stage survey Gaussian identity | Independent PSU algebra; six deterministic stacked Stata full matrices; analytic linear slope and contrast SE | Gate 3 independent full-design covariance and all three fixed-weight contrast SEs pass; retain fixed averaging interpretation | `test-survey.R`, `test-stata-survey-reference.R` |
| One-stage survey quasibinomial logit | Native design influence/block reproduction; independent PSU algebra; returned Stata full covariance | Gate 3 independent expected-information/PSU covariance and nonlinear contrast SEs pass with/without FPC; no ordinary case bootstrap | `test-survey.R`, `test-stata-survey-binary-reference.R` |
| One-stage survey quasibinomial probit | Native expected-information blocks and independent design algebra | Gate 3 independent expected-information covariance and nonlinear contrast SEs pass with/without FPC; preserve Stata information-convention distinction | same; `tools/SURVEY-BINARY-FINAL-VALIDATION-RESULTS-20260909.md` |
| MI coefficient pooling, `suest_mi` | Explicit Rubin within/between/full covariance algebra; `mitools::MIcombine`; coefficient-contrast SE | Nonlinear MI predictions/effects remain unsupported; no expansion in this release pass | `test-mi.R` |
| Pweights in admitted ordinary families | Returned linear/logit/probit/Poisson/NB/categorical evidence; invariance checks; full ancillary blocks | Keep supported family list and documented ancillary weight normalization; distinguish fitting weights from averaging weights | `test-pweights.R`, `test-offsets.R`, acceptance pweight tests |
| Cluster covariance | Full grouped-score algebra; external ordinary/2SLS/weighted systems | Gate 1 independently stacked shared/disjoint-cluster geometry passes | `test-cluster.R`, `test-stata-cluster-reference.R` |
| ID alignment, multiple models, heterogeneous outcomes | Explicit/composite IDs, pairwise overlap, three models, missingness, differing outcome scales | Broad combinations often prove execution only; choose distinct models and independently verify cross-block magnitudes | `test-core.R`, `test-multiple-models.R`, `test-comparison-universe.R` |
| Offsets | Scalar, weighted, ordered analytic-score regressions | Confirm each advertised route's own exclusions; generic “offsets supported” must not imply GLMM support | `test-offsets.R`; current README/model contracts |

## Prioritized next checks

1. Completed: gate 1, runtime reproductions, and tested safeguards for
   RV-01/RV-02. The 40 stacked cases and bootstrap pass again after the changes.
2. Completed for the declared scope: gate 2 independently checks beta, survival,
   censored/truncated, heteroskedastic, and ZIP/ZINB full covariance and
   prediction/effect contrast SEs (108 cases). RV-03 is fixed. Default truncated
   optimizer warnings/failure and the clean NR results remain distinguishable.
   See RELEASE-VALIDATION-GATE2-20260929.md for targets and exclusions.
3. Completed for the declared scope: gate 3, 100 cases and 300 analytic
   contrasts for panel/GEE/integrated and survey routes. FE average-level
   inference at training means remains explicitly unsupported. See the gate-3
   report for conditional covariance inputs and reviewed native warnings.
   Completed: gate 4, matched glm/glm2 and polr/clm contrasts, independent
   bivariate-probit and IV covariance/contrasts, and four convergence safeguards.
   All 28 cases/132 contrasts pass, with one reviewed native trial-CDF warning
   case. Source and installed suites pass 1,680; acceptance 144/144; Linux
   check clean. Freeze and verify the exact candidate on supported platforms.
4. Preserve the strong random-slope independent audits. If remaining numerical
   ambiguity affects the release decision, use a targeted independent gradient
   or paired cluster-bootstrap check, not repeated broad Stata convergence runs.
5. Check the exact frozen candidate on the supported platforms with dependencies
   available, no unexplained skips, and each warning classified. Historical
   passing results are not a substitute for checks of changed statistical code.

A baseline passing test demonstrates the assertion it contains. A gap in this
matrix is a request for additional evidence, not permission to relax a numerical
tolerance or rewrite an external reference to agree with the package.
