# Model-adapter architecture

`suest()` separates the statistical model family from the fitting engine.
This prevents package-specific object layouts from leaking into the joint
covariance and `marginaleffects` interfaces.

## Internal adapter functions

- `.suest_model_adapter(model)` identifies the statistical type and engine.
- `.suest_model_frame(model, engine)` recovers the estimation rows.
- `.suest_category_levels(model, engine)` returns ordered or nominal outcome
  categories.
- `.suest_extract_parameters(model, type, engine)` returns parameters in score
  order.
- `.suest_model_components(model, type, engine)` returns observation-level
  scores and the bread matrix.
- `.suest_set_parameters(model, parameters, type, engine)` updates a fitted
  object for numerical differentiation by `marginaleffects`.
- `.suest_predict_probabilities()` and `.suest_predict_values()` standardize
  predictions across engines.

## Supported engines

| Engine | Statistical types | Initial restrictions |
|---|---|---|
| `stats::lm()` | linear | ordinary least squares; unweighted and pweighted systems include Stata-compatible `lnvar` |
| `stats::glm()` | binary, count, fractional response, and other supported GLMs | ordinary maximum likelihood or quasi-likelihood |
| `glm2::glm2()` | supported GLMs | ordinary maximum likelihood or quasi-likelihood |
| `MASS::glm.nb()` | negative binomial | log link |
| `MASS::polr()` | ordered logit, ordered probit | standard threshold model |
| `ordinal::clm()` | ordered logit, ordered probit | flexible thresholds, no scale or nominal formula |
| `nnet::multinom()` | multinomial logit | `summ = 0` |
| `survival::survreg()` | parametric survival and censored regression | one common scale parameter |
| `betareg::betareg()` | beta regression | standard beta distribution and supported mean links |
| `pscl::zeroinfl()` | zero-inflated Poisson and negative binomial | response-scale predictions |
| `truncreg::truncreg()` | truncated Gaussian regression | left or right truncation |
| `censReg::censReg()` | censored Gaussian (Tobit) regression | latent-mean predictions |
| `fixest::feols()` | instrumental-variable 2SLS | unweighted, no absorbed fixed effects, original data available |
| `Rchoice::hetprob()` | heteroskedastic binary probit/logit | response/link predictions; no estimation weights |
| `Rchoice::ivpml()` | maximum-likelihood instrumental-variable probit | one continuous endogenous regressor; average-structural response prediction |
| `mvProbit::mvProbit()` | bivariate probit | same regressors in both equations; observed final Hessian required; joint-success prediction |
| `plm::plm()` | linear panel fixed, between, and random effects | individual effects; unweighted; Swamy-Arora random-effects method; same panel type within a system; unbalanced RE can differ slightly from Stata |
| `nlme::lme()` | Gaussian random-intercept panel ML | one grouping level and random intercept; ML; unweighted; independent homoskedastic residuals |
| `geepack::geeglm()` | population-averaged GEE | unweighted numeric Gaussian identity, binary logit/probit/cloglog, Poisson log; independence/exchangeable correlation; panel or higher nested clusters |

## Explicitly rejected GLM-like estimators

Objects fitted using bias reduction, adjusted score equations, Firth
corrections, or penalized likelihood are rejected. Inheriting from `"glm"` is
not enough: the observation-level estimating equations must match the ordinary
maximum-likelihood GLM scores used by the SUEST covariance.

## Candidate engines for later versions

The adapter layer is designed to support additional engines after separate
validation:

1. cross-sectional fixed-parameter `Rchoice::Rchoice()` models;
2. `fixest` models without absorbed fixed effects, followed by a separate
   investigation of absorbed effects;
3. ordinary maximum-likelihood `mclogit::mblogit()` models;
4. selected large-data engines when estimation rows and score contributions
   can be recovered safely.

Broader survey-weighted and general multilevel models require separate
covariance designs and are intentionally deferred to later development phases.
The admitted survey subset is documented below; the admitted `nlme::lme()` route
is the narrower random-intercept Gaussian model that maps directly to `xtreg, mle`.

## Joint clustering

Clustering is a system-level covariance operation, not a model adapter. After
scores are aligned to the union of the estimation samples, `suest()` sums the
joint score rows by `cluster` and forms the sandwich meat from those cluster
totals. This permits clusters to span partially overlapping or completely
disjoint model samples. Component prediction and coefficient-replacement
methods are unchanged, so clustered objects use the same `marginaleffects`
interface.


## Restricted one-stage survey route (2026-09-08)

`suest(..., survey_design = design, observation_id = "id")` dispatches before
ordinary adapters. `R/survey.R` reconstructs coefficient influences for Gaussian
identity and binary `quasibinomial()` logit/probit `svyglm()` fits, aligns them to
the full common design, verifies each native covariance block, and calls
`survey::svyrecvar()` once on the stacked influences. Gaussian and binary fits
are not mixed in one survey system. The ordinary adapter rejects `svyglm` to
prevent accidental use of ordinary GLM covariance. Survey systems retain
`suest_model` prediction integration with engine `survey::svyglm`; the parameter
vector contains coefficients only. See `SURVEY-STEP-STATUS-20260908.md` and
`SURVEY-BINARY-STEP-STATUS-20260908.md` for restrictions and validation.
