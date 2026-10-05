#' Combine fitted models with seemingly unrelated estimation
#'
#' `suest()` combines two or more separately fitted models into one object with
#' a joint robust covariance matrix, so you can test whether predictions and
#' marginal effects differ across the models. Pass the result directly to
#' [marginaleffects::predictions()], [marginaleffects::avg_comparisons()],
#' [marginaleffects::avg_slopes()], and [marginaleffects::hypotheses()].
#'
#' @param ... Two or more supported fitted model objects. For backward
#'   compatibility, a character vector supplied as the third unnamed argument
#'   is interpreted as `model_names` for a two-model system.
#' @param model_names Optional character vector containing one display name per model.
#'   By default, the object names supplied in the call are used.
#' @param observation_id Optional observation identifier used to align models
#'   fitted from different data objects. Supply one or more column names found
#'   in every model's original data, such as `"id"` or `c("id", "wave")`, or
#'   a list containing an ID vector, matrix, or data frame already aligned to
#'   each model's estimation sample. IDs must be complete and unique within
#'   each model. The default `NULL` uses data-source identity and model-frame
#'   row names.
#' @param cluster Optional cluster identifier for a joint cluster-robust
#'   covariance matrix. Supply one or more column names found in every model's
#'   original data, or a list containing a cluster vector, matrix, or data
#'   frame aligned to each model's estimation sample. Cluster IDs must be
#'   complete. Shared observations must have the same cluster ID in every
#'   model; disjoint observations may still share clusters across models.
#'   Supported panel systems default to the panel identifier when `cluster` is
#'   omitted, and systems that combine random-effects models with ordinary
#'   models default to the random-effects grouping variable. Supplied clusters
#'   for these models must contain whole panels.
#' @param weight_type Optional weight interpretation. The default `NULL`
#'   preserves the unweighted behavior and rejects nonunit estimation weights.
#'   Use `"pweight"` to treat model weights as sampling weights. In version
#'   0.1.4, pweights are supported for linear, binary logit/probit, Poisson,
#'   negative-binomial, ordered logit/probit, and multinomial logit models.
#' @param survey_design Full common one-stage `survey::svydesign()` object for
#'   supported `survey::svyglm()` models: Gaussian identity, binary
#'   `quasibinomial()` logit/probit/cloglog, or `quasipoisson()` log.
#'   Requires `observation_id` column names.
#'   Retain the design before model-specific domain subsetting. Specify weights
#'   and clusters in this design, without `weight_type` or `cluster`. Ordinary
#'   models leave this argument `NULL`.
#' @param object,x A `"suest_model"` object.
#'
#' @return An object of class `"suest_model"` containing the fitted models,
#'   their joint coefficient vector, a joint model-robust covariance matrix,
#'   and sample-alignment information.
#'
#' @details
#' Two or more models are supported. Models can use identical, partially
#' overlapping, or completely disjoint samples. When the model calls refer to
#' the same data source, model-frame row names identify overlapping
#' observations. Models fitted from different data objects are treated as
#' disjoint by default because shared observations cannot be inferred safely.
#' Use `observation_id` to identify their common observations explicitly.
#' `suest()` warns when such models share row names with identical values, and
#' stops when models fit to the same data object disagree on rows with the same
#' row names (for example, after the data were filtered and renumbered).
#'
#' When `newdata` is omitted in `marginaleffects` functions, each model is
#' averaged over its own estimation sample, as with
#' `newdata = suest_newdata(fit)`.
#'
#' All combinations of the supported scalar-response models can be combined,
#' including models whose response variables or response scales differ. All
#' combinations of the supported categorical-response models can be combined,
#' and scalar and categorical models may appear in the same system. Results on
#' different response scales are labeled separately for `marginaleffects`.
#'
#' Linear models include an ancillary `lnvar` parameter, the log of the error
#' variance (`log(RSS / df.residual)` for unweighted fits). It has no effect on
#' predictions or marginal effects.
#' Negative-binomial models include `log(theta)` in the joint parameter vector.
#' Ordered and multinomial models use analytic score and observed-information
#' calculations for stable robust covariance estimation.
#' Probit and other GLMs with non-canonical links, and beta regressions, use
#' observed information in the sandwich bread, as Stata does.
#'
#' Aliased parameters are not supported. Nonunit weights are rejected unless
#' `weight_type = "pweight"`. Pweights must be finite and
#' strictly positive. For observations included in any pair of models, evaluated
#' weights must agree; weights may differ for observations unique to either
#' model, including completely disjoint samples. Pweight support is available
#' for linear, binary logit/probit, Poisson, negative-binomial, ordered
#' logit/probit, and multinomial logit models.
#'
#' Bias-reduced, adjusted-score, Firth, and penalized
#' GLM fits are rejected because they do not use ordinary maximum-likelihood
#' score equations. Ordinary `glm` and `glm2` fits must have converged.
#' Penalized `nnet::multinom` fits with nonzero `decay` are also rejected.
#' Beta regression requires `type = "ML"`; bias-reduced and bias-corrected
#' beta fits and penalized `survreg` fits are not supported. Conventional
#' `survreg` fits with `robust = TRUE` use their model-based information for
#' the joint sandwich covariance.
#'
#' @section Survey models:
#' Survey support combines coefficients from Gaussian identity-link, binary
#' `quasibinomial()` logit/probit/cloglog, and `quasipoisson()` log
#' `svyglm()` fits under one common one-stage design, in any combination.
#' Strata and first-stage finite-population corrections are supported.
#' Model-specific subsets and missing outcomes are aligned using observation
#' IDs, with zero influence outside each model's estimation sample. The full
#' design retains PSUs outside all model samples. Each native coefficient
#' covariance must be reproduced before the joint matrix is returned. Survey
#' fits contain coefficients only and do not add ancillary parameters.
#'
#' Replicate-weight, multistage, two-phase, calibrated, raked, post-stratified,
#' and PPS designs are unsupported. Lonely-PSU options are restricted to
#' `fail`, `remove`, or `certainty`, with `survey.adjust.domain.lonely = FALSE`.
#' A singleton stratum with a certainty FPC is supported under `fail`.
#' Extra fitting weights are unsupported; place offsets in the model formula.
#'
#' Predictions and effects use the joint design-based coefficient covariance.
#' Averaging treats the supplied covariate distribution as fixed; it does not
#' add design uncertainty from estimating that distribution. Use
#' `suest_newdata()` and `wts = ".suest_weight"` for model-specific weighted
#' averages. Inference uses the usual asymptotic normal default. The full
#' design's degrees of freedom are recorded in `$survey$design_df`; no automatic
#' survey t or F adjustment is applied.
#'
#' @section Supported models:
#' ## Single-level models
#' * [stats::lm()]
#' * binary logit, probit, and complementary-log-log models from [stats::glm()]
#'   or `glm2::glm2()`
#' * ordered logit and probit models from [MASS::polr()]
#' * ordered logit and probit models from `ordinal::clm()` with flexible
#'   thresholds, proportional effects, and no scale model
#' * multinomial logit models from [nnet::multinom()]
#' * Poisson log-link models from [stats::glm()] or `glm2::glm2()`
#' * negative-binomial log-link models from [MASS::glm.nb()]
#' * Poisson and negative-binomial zero-inflated models from
#'   `pscl::zeroinfl()`
#' * other GLMs using identity, log, logit, probit, complementary-log-log, or
#'   log-log links
#' * left-, right-, and two-limit censored Gaussian regressions from
#'   `censReg::censReg()`; response-scale predictions are the latent mean
#' * parametric survival and censored-regression models from
#'   `survival::survreg()` with a common scale, including Gaussian interval
#'   regression
#' * truncated Gaussian regressions from `truncreg::truncreg()`
#' * fractional-response GLMs using `quasibinomial()` with logit, probit,
#'   complementary-log-log, or a user-supplied log-log link
#' * beta regressions from `betareg::betareg()` using logit, probit,
#'   complementary-log-log, or log-log mean links
#' * heteroskedastic binary probit and logit from `Rchoice::hetprob()`
#' * bivariate probit from `mvProbit::mvProbit()` when both equations use the
#'   same regressors; fit with `intGrad = TRUE` and `finalHessian = TRUE`
#' * unweighted two-stage least squares from `fixest::feols()` without
#'   absorbed fixed effects; the original data object must remain available
#' * maximum-likelihood instrumental-variable probit from `Rchoice::ivpml()`;
#'   response predictions use the average structural probability
#'
#' ## Selection models
#' * unweighted maximum-likelihood Heckman selection models from
#'   `sampleSelection::selection(method = "ml")`, Stata's `heckman`. The
#'   ancillary parameters use Stata's scale, `lnsigma` and `athrho`, and the
#'   equations are labeled `selection:` and `outcome:`. Predictions are the
#'   outcome equation's linear prediction, Stata's default
#'
#' ## Generalized ordered models
#' * unweighted `VGAM::vglm()` models with the `cumulative()` family and a
#'   logit, probit, or cloglog link: generalized ordered logit (Stata's
#'   `gologit2`), including the proportional (`parallel = TRUE`) and partial
#'   proportional odds forms (for example `parallel = FALSE ~ x`). Use
#'   `reverse = TRUE` for gologit2's `P(Y > j) = F(eta_j)` parameterization;
#'   both directions are supported. suest computes the scores analytically
#'   and the observed information (as Stata does) by differentiating them;
#'   VGAM's own covariance uses expected information. Predictions are
#'   category probabilities
#'
#' ## Panel models
#' * unweighted individual fixed-effects, between-effects, and Swamy-Arora
#'   random-effects linear panel models from `plm::plm()`. Within-model prediction
#'   uncertainty conditions on estimation-sample means and propagates slope
#'   uncertainty. At matching evaluation means, the covariance is structurally
#'   degenerate: average-level confidence intervals and hypothesis tests are
#'   unsupported, and numerical SEs can be missing or nearly zero. Slope and
#'   finite-change comparisons remain supported.
#' * unweighted individual random-intercept binary logit and probit models from
#'   `pglm::pglm()` fitted with `model = "random"`, `effect = "individual"`,
#'   and `R = 12`; response predictions integrate over the random effect
#' * unweighted gamma random-effects Poisson log models from `pglm::pglm()`
#'   fitted with `model = "random"`, `effect = "individual"`, and
#'   `other = "sd"`; the final parameter is exposed as gamma variance `alpha`
#' * unweighted fixed-effects Poisson models from `fixest::fepois()` with one
#'   absorbed fixed effect (Stata's `xtpoisson, fe`), clustered on that
#'   effect with no `G/(G-1)` factor, as in Stata's native
#'   `xtpoisson, fe vce(robust)`. Each unit's fixed effect is a function of the coefficients,
#'   `log(sum(y_i)) - log(sum(exp(x_i b)))`, so response predictions
#'   `exp(x b + alpha_i)` match `predict()` for fixest and carry the
#'   coefficients' uncertainty; newdata must contain the fixed-effect
#'   variable. Link-scale effects equal the coefficients, as in Stata. These
#'   models combine only with other fixed-effects Poisson models
#' * unweighted GEE from `geepack::geeglm()`: Gaussian identity,
#'   binary logit/probit/cloglog, and Poisson log, with independence or
#'   exchangeable correlation and numeric outcomes
#'
#' ## Multilevel models
#' * unweighted single-level random-intercept Gaussian models from
#'   `nlme::lme()` fitted with `method = "ML"`, matching Stata's
#'   `xtreg, mle` (full observed-information bread). For Stata's `mixed`
#'   layout, fit the same model with `lme4::lmer(..., REML = FALSE)`; the
#'   standard errors differ only slightly
#' * unweighted binomial (logit, probit, cloglog), Poisson-log, negative-binomial
#'   NB2 log, Gaussian-identity, and Gamma-log models from
#'   `glmmTMB::glmmTMB()` with one grouping variable and one conditional random
#'   intercept; the final parameter is the log random-intercept standard
#'   deviation, and response predictions integrate over the Gaussian random
#'   effect (probit `pnorm(eta/sqrt(1 + v))`, log links `exp(eta + v/2)`,
#'   Gaussian `eta`, logit and cloglog by quadrature, where `v` is the
#'   random-effect variance). Gaussian models add `log_sigma_e`, the log
#'   residual standard deviation, and Gamma models `log_shape`, the log shape
#'   parameter (Stata's `/logs` equals `-log_shape/2`), before the
#'   random-effect parameters. Gaussian random-effects models from `lmer()` or
#'   `glmmTMB()` use the covariance layout of Stata's `mixed`: the fixed
#'   effects' bread is `(X'V^-1 X)^-1`, the variance parameters' bread is
#'   their block of the full observed-information inverse, and the two blocks
#'   are uncorrelated. NB2 models must use the default constant dispersion model
#'   (`dispformula = ~1`); its estimated log size parameter `log_phi` precedes
#'   `log_sigma`, with conditional variance `mu + mu^2/exp(log_phi)`.
#'   Weights, offsets, and zero inflation are unsupported;
#'   NB2 additionally excludes mapped or constrained parameters. Its Laplace
#'   log likelihood must exceed the zero-random-effect NB2 log likelihood at
#'   the same fixed effects and dispersion by more than
#'   `sqrt(.Machine$double.eps) * max(1, abs(logLik(model)))`.
#'   This numerical boundary check is not a significance test and does not
#'   guarantee an interior global maximum
#' * the same `glmmTMB` families with one correlated random
#'   intercept and numeric slope, `(1 + x | id)`, using an unstructured
#'   covariance matrix. The slope must be a single untransformed numeric
#'   column with a syntactically valid name. The nuisance parameters are
#'   `log_sd_intercept`, `log_sd_slope`, and `atanh_rho`. NB2 requires the
#'   default estimated constant dispersion (`dispformula = ~1`), includes
#'   `log_phi` before these three parameters, and uses the NB2 boundary check
#'   above; Gaussian and Gamma models include `log_sigma_e` or `log_shape`
#'   there. Response predictions use the random-intercept formulas with
#'   `v = var_intercept + 2*x*cov_intercept_slope + x^2*var_slope`; binomial
#'   responses must be Bernoulli.
#'   Link predictions are `X beta`. Near-singular random
#'   covariance is rejected in a centered, standardized predictor basis.
#'   Weights, offsets, zero inflation, constraints, diagonal covariance,
#'   multiple slopes, and slope-only terms are unsupported. The slope variable must be supplied for response predictions even when
#'   it is absent from the fixed-effects formula. Native model covariance is
#'   used as supplied. Poorly scaled predictors can produce inaccurate native
#'   numerical curvature even with convergence and a positive-definite Hessian;
#'   center and rescale continuous predictors before fitting, then compare
#'   predictions and effects after converting them to the same original units.
#'   For logit random-slope `avg_slopes()`, use
#'   `numderiv = list("fdcenter", eps = 1e-4)` and compare results at nearby
#'   steps (for example, `5e-5` and `2e-4`). An unstable standard error
#'   requires further investigation; changing the finite-difference step
#'   does not repair an inaccurate native model covariance
#'
#' * unweighted `lme4::lmer()` models fit with `REML = FALSE`, and
#'   `lme4::glmer()` binomial (logit, probit, cloglog) and Poisson-log models,
#'   with one grouping variable and a random intercept or one correlated
#'   numeric slope. They use the same parameters and predictions as the
#'   matching `glmmTMB` models. suest computes their group scores and observed
#'   information from its own per-group likelihood: exact for `lmer`, Laplace
#'   or adaptive Gauss-Hermite with the fit's `nAGQ` for `glmer`.
#'   `glmer(..., nAGQ = 7)` uses mode-curvature adaptive quadrature with 7
#'   points, Stata's `intmethod(mcaghermite) intpoints(7)`. lme4's default
#'   Laplace fits stop the inner mode search early; for closer agreement with
#'   other software, fit with `control = lme4::glmerControl(tolPwrss = 1e-12)`
#'
#' * unweighted random-intercept ordered logit and probit models from
#'   `ordinal::clmm()` with one grouping variable and flexible thresholds
#'   (Stata's `meologit`, `meoprobit`, `xtologit`, `xtoprobit`), fit with
#'   Laplace (`nAGQ = 1`) or adaptive quadrature (`nAGQ > 1`; 7 points match
#'   Stata's `intmethod(mcaghermite) intpoints(7)`). Parameters are the
#'   thresholds, the coefficients, and `log_sigma`. suest computes scores and
#'   observed information from its own per-group likelihood, which
#'   reproduces clmm's. Predictions are category probabilities integrated
#'   over the random intercept
#'
#' ## Combining random-effects and ordinary models
#' Random-effects models from `nlme::lme()`, `pglm::pglm()`,
#' `glmmTMB::glmmTMB()`, lme4, and `ordinal::clmm()` can be combined with each other, including different
#' families, and with linear, binary logit, probit, and cloglog, Poisson,
#' negative binomial, and ordered logit and probit models, as in Stata's
#' suest2. Every random-effects model must use the same grouping variable,
#' and the ordinary models' data must contain it. The system is clustered on
#' that group, or on `cluster`, which must contain whole groups. Each model's
#' block of the covariance then equals that model's own covariance clustered
#' on the group, including the `(N - 1)/(N - k)` adjustment for linear
#' models. `plm` and `geeglm` models combine only with models of the same
#' type.
#'
#' ## Survey models
#' * restricted survey-weighted Gaussian identity, binary `quasibinomial()`
#'   logit/probit/cloglog, and `quasipoisson()` log models from
#'   `survey::svyglm()`, in any combination
#'
#' @examples
#' dat <- mtcars
#' dat$am <- factor(dat$am)
#'
#' model1 <- glm(am ~ wt, family = binomial(), data = dat)
#' model2 <- glm(am ~ wt + hp, family = binomial(), data = dat)
#'
#' fit <- suest(model1, model2, model_names = c("Base", "Adjusted"))
#' fit
#'
#' effects <- marginaleffects::avg_comparisons(fit, variables = "wt", newdata = dat)
#' marginaleffects::hypotheses(effects, hypothesis = difference ~ revpairwise)
#'
#' @export
suest <- function(
    ...,
    model_names = NULL,
    observation_id = NULL,
    cluster = NULL,
    weight_type = NULL,
    survey_design = NULL) {
  if (!requireNamespace("sandwich", quietly = TRUE))
    stop("Package 'sandwich' is required.", call. = FALSE)
  if (!requireNamespace("marginaleffects", quietly = TRUE))
    stop("Package 'marginaleffects' is required.", call. = FALSE)

  weight_type <- .suest_validate_weight_type(weight_type)
  models <- list(...)
  model_calls <- as.list(substitute(list(...)))[-1L]

  legacy_names <- is.null(model_names) && length(models) == 3L &&
    is.character(models[[3L]]) && length(models[[3L]]) == 2L
  if (legacy_names) {
    model_names <- models[[3L]]
    models <- models[1:2]
    model_calls <- model_calls[1:2]
  }

  n_models <- length(models)
  if (n_models < 2L)
    stop("Supply at least two fitted models.", call. = FALSE)
  models <- lapply(models, function(model)
    .suest_prepare_heckman(.suest_wrap_vglm(.suest_wrap_lme4(.suest_unexclude(model)))))

  if (is.null(model_names)) {
    model_names <- vapply(
      model_calls,
      function(x) paste(deparse(x, width.cutoff = 500L), collapse = ""),
      character(1)
    )
    if (any(nchar(model_names) > 30) || anyDuplicated(model_names))
      model_names <- paste("Model", seq_len(n_models))
  }
  if (length(model_names) != n_models || anyNA(model_names) ||
      any(model_names == ""))
    stop(
      "'model_names' must contain one nonempty name per model.",
      call. = FALSE
    )

  model_names <- make.unique(as.character(model_names))
  names(models) <- model_names

  survey_models <- vapply(models, inherits, logical(1), "svyglm")
  if (any(survey_models) || !is.null(survey_design))
    return(.suest_survey_system(models, model_names, observation_id,
      survey_design, cluster, weight_type, match.call()))

  adapters <- lapply(models, .suest_model_adapter)
  model_types <- vapply(adapters, `[[`, character(1), "type")
  model_engines <- vapply(adapters, `[[`, character(1), "engine")

  unsupported <- which(model_types == "unsupported")
  if (length(unsupported))
    stop(
      paste0(
        "Each model must be a supported stats::lm, stats::glm, glm2::glm2, ",
        "MASS::glm.nb, survival::survreg, betareg::betareg, pscl::zeroinfl, ",
        "truncreg::truncreg, censReg::censReg, MASS::polr, ordinal::clm, ",
        "nnet::multinom, mvProbit::mvProbit, Rchoice::hetprob/ivpml, a supported ",
        "plm::plm linear ",
        "panel model, an nlme::lme random-intercept ML model, a pglm::pglm ",
        "random-effects binary or Poisson model, a supported ",
        "glmmTMB::glmmTMB or lme4 random-effects model, a VGAM::vglm generalized ordered model, a sampleSelection::selection ML Heckman model, an ordinal::clmm random-intercept ordered model, a supported ",
        "geepack::geeglm model, ",
        "a fixest::feols IV model, or a fixest::fepois fixed-effects Poisson model. Unsupported model '",
        model_names[unsupported[1L]], "' has class '",
        class(models[[unsupported[1L]]])[1L], "'."
      ),
      call. = FALSE
    )

  plm_panel_types <- c("panel_fe", "panel_be", "panel_re")
  glmm_types <- unique(model_types[.suest_is_glmm(model_types)])
  panel_types <- c(
    plm_panel_types, "panel_ml", "panel_gee", "panel_logit_re",
    "panel_probit_re", "panel_poisson_re", "panel_poisson_fe", "re_ologit",
    "re_oprobit", glmm_types
  )
  # Random-effects models can be combined with each other and with these
  # ordinary models, all clustered on the shared highest-level group (as in
  # Stata suest2's multilevel route).
  multilevel_types <- c(
    "panel_ml", "panel_logit_re", "panel_probit_re", "panel_poisson_re",
    "re_ologit", "re_oprobit", glmm_types
  )
  partner_types <- c(
    "lm", "logit", "probit", "cloglog", "poisson", "negbin", "ologit", "oprobit"
  )
  is_panel <- model_types %in% panel_types
  heterogeneous <- any(is_panel) && length(unique(model_types)) > 1L
  if (heterogeneous && (!all(model_types[is_panel] %in% multilevel_types) ||
      !all(model_types[!is_panel] %in% partner_types)))
    stop(
      paste0(
        "Random-effects models (nlme::lme, pglm, glmmTMB) can be combined with ",
        "each other and with linear, logit, probit, cloglog, Poisson, negative ",
        "binomial, ordered logit, and ordered probit models. Other panel ",
        "models (plm, geeglm) can be combined only with models of the same type."
      ),
      call. = FALSE
    )

  pair <- .suest_pair_info(model_types)

  if (any(model_types == "multinom") &&
      !requireNamespace("nnet", quietly = TRUE))
    stop(
      "Package 'nnet' is required for multinomial models.",
      call. = FALSE
    )

  if (any(model_types == "survreg") &&
      !requireNamespace("survival", quietly = TRUE))
    stop(
      "Package 'survival' is required for survreg models.",
      call. = FALSE
    )

  if (any(model_types == "betareg") &&
      !requireNamespace("betareg", quietly = TRUE))
    stop(
      "Package 'betareg' is required for beta-regression models.",
      call. = FALSE
    )

  if (any(model_types %in% c("zip", "zinb")) &&
      !requireNamespace("pscl", quietly = TRUE))
    stop(
      "Package 'pscl' is required for zero-inflated models.",
      call. = FALSE
    )

  if (any(model_types == "truncreg") &&
      !requireNamespace("truncreg", quietly = TRUE))
    stop(
      "Package 'truncreg' is required for truncated-regression models.",
      call. = FALSE
    )

  if (any(model_types == "censreg") &&
      !requireNamespace("censReg", quietly = TRUE))
    stop(
      "Package 'censReg' is required for censored-regression models.",
      call. = FALSE
    )

  if (any(model_types == "heckman") &&
      !requireNamespace("sampleSelection", quietly = TRUE))
    stop(
      "Package 'sampleSelection' is required for Heckman selection models.",
      call. = FALSE
    )

  if (any(model_types %in% c("ivreg", "panel_poisson_fe")) &&
      !requireNamespace("fixest", quietly = TRUE))
    stop(
      "Package 'fixest' is required for instrumental-variable models.",
      call. = FALSE
    )

  if (any(model_types %in% c("hetprobit", "hetlogit", "ivprobit")) &&
      !requireNamespace("Rchoice", quietly = TRUE))
    stop(
      "Package 'Rchoice' is required for heteroskedastic or IV probit models.",
      call. = FALSE
    )

  if (any(model_types == "biprobit") &&
      (!requireNamespace("mvProbit", quietly = TRUE) ||
       !requireNamespace("mvtnorm", quietly = TRUE)))
    stop(
      "Packages 'mvProbit' and 'mvtnorm' are required for bivariate probit models.",
      call. = FALSE
    )

  if (any(model_types %in% plm_panel_types) &&
      !requireNamespace("plm", quietly = TRUE))
    stop(
      "Package 'plm' is required for linear panel models.",
      call. = FALSE
    )

  if (any(model_types == "panel_ml") &&
      !requireNamespace("nlme", quietly = TRUE))
    stop(
      "Package 'nlme' is required for random-intercept ML panel models.",
      call. = FALSE
    )

  if (any(model_types %in% c(
        "panel_logit_re", "panel_probit_re", "panel_poisson_re"
      )) &&
      !requireNamespace("pglm", quietly = TRUE))
    stop(
      "Package 'pglm' is required for random-effects panel models.",
      call. = FALSE
    )

  if (any(.suest_is_glmm(model_types)) &&
      !requireNamespace("glmmTMB", quietly = TRUE))
    stop(
      "Package 'glmmTMB' is required for GLMM random-effects models.",
      call. = FALSE
    )

  categorical <- model_types %in% c("ologit", "oprobit", "multinom", "gologit", "re_ologit", "re_oprobit")

  category_levels <- Map(
    .suest_category_levels,
    models,
    model_engines
  )

  model_frames <- Map(
    .suest_model_frame,
    models,
    model_engines
  )

  for (i in seq_along(models)) {
    model <- models[[i]]
    type <- model_types[i]
    engine <- .suest_engine_name(model_engines[[i]])
    mf <- model_frames[[i]]

    parameters <- .suest_extract_parameters(model, type, engine)
    if (anyNA(parameters))
      stop(
        sprintf(
          "Model '%s' contains aliased or missing parameters.",
          model_names[i]
        ),
        call. = FALSE
      )

    weights <- .suest_model_weights(mf, model, engine)

    if (is.null(weight_type)) {
      if (any(weights != 1))
        stop(
          paste0(
            "Nonunit estimation weights require ",
            "weight_type = \"pweight\"."
          ),
          call. = FALSE
        )
    } else {
      if (!type %in% c(
        "lm", "logit", "probit", "poisson", "negbin",
        "ologit", "oprobit", "multinom"
      ))
        stop(
          paste0(
            "Pweight support is not available for this model type. See ",
            "the supported-model documentation for current coverage."
          ),
          call. = FALSE
        )
      .suest_validate_pweights(weights, model_names[i])
    }

    y <- stats::model.response(mf)

    if (type %in% c(
          "lm", "panel_fe", "panel_be", "panel_re", "panel_ml", "panel_gee"
        ) &&
        (!is.numeric(y) || is.matrix(y)))
      stop("Linear-model outcomes must be numeric vectors.",
           call. = FALSE)

    if (type == "ivreg" && (!is.numeric(y) || is.matrix(y)))
      stop(
        "Instrumental-variable outcomes must be numeric vectors.",
        call. = FALSE
      )

    if ((type %in% c("panel_poisson_re", "panel_poisson_fe") || .suest_glmm_count(type)) &&
        (!is.numeric(y) || is.matrix(y) ||
         any(y < 0 | abs(y - round(y)) > 1e-8)))
      stop(
        "Random-effects count outcomes must be nonnegative integer counts.",
        call. = FALSE
      )

    if (type %in% c(
      "logit", "probit", "cloglog", "ivprobit", "hetprobit", "hetlogit",
      "panel_logit_re", "panel_probit_re"
    ) || .suest_glmm_binary(type)) {
      binary <- if (is.factor(y)) {
        nlevels(y) == 2L
      } else if (is.matrix(y)) {
        FALSE
      } else {
        all(stats::na.omit(unique(y)) %in% c(0, 1))
      }

      if (!binary)
        stop(
          paste0(
            "Binary-response outcomes must ",
            "be two-level factors or numeric 0/1."
          ),
          call. = FALSE
        )
    }

    if (type == "biprobit") {
      binary_pair <- is.matrix(y) && ncol(y) == 2L &&
        all(stats::na.omit(unique(as.numeric(y))) %in% c(0, 1))
      if (!binary_pair)
        stop(
          "Bivariate-probit outcomes must be two numeric 0/1 columns.",
          call. = FALSE
        )
    }

    if (type %in% c(
      "fraclogit", "fracprobit", "fraccloglog", "fracloglog"
    )) {
      fractional <- is.numeric(y) && !is.matrix(y) && all(is.finite(y)) &&
        all(y >= 0 & y <= 1)
      if (!fractional)
        stop(
          "Fractional-response outcomes must be numeric values from zero to one.",
          call. = FALSE
        )
    }

    if (type == "betareg") {
      beta_outcome <- is.numeric(y) && !is.matrix(y) && all(is.finite(y)) &&
        all(y > 0 & y < 1)
      if (!beta_outcome)
        stop(
          "Beta-regression outcomes must be numeric values strictly between zero and one.",
          call. = FALSE
        )
    }

    if (type %in% c("poisson", "negbin", "zip", "zinb")) {
      count <- is.numeric(y) && !is.matrix(y) && all(is.finite(y)) &&
        all(y >= 0) &&
        all(abs(y - round(y)) < sqrt(.Machine$double.eps))

      if (!count)
        stop(
          paste0(
            "Poisson, negative-binomial, and zero-inflated outcomes must be ",
            "nonnegative integer counts."
          ),
          call. = FALSE
        )
    }

    if (type %in% c("ologit", "oprobit") &&
        (!is.ordered(y) || nlevels(y) < 3L))
      stop(
        "Ordered-model outcomes must be ordered factors with at least three levels.",
        call. = FALSE
      )

    if (type == "multinom") {
      cf <- stats::coef(model)
      levels <- category_levels[[i]]
      valid_multinom <- !is.null(levels) &&
        length(levels) >= 3L &&
        is.matrix(cf) &&
        nrow(cf) == length(levels) - 1L

      if (!valid_multinom)
        stop(
          "Multinomial models must have at least three outcome categories.",
          call. = FALSE
        )
    }
  }

  model_weights <- Map(
    .suest_model_weights,
    model_frames,
    models,
    model_engines
  )

  components <- Map(
    function(model, type, engine)
      .suest_model_components(model, type, engine, weight_type),
    models,
    model_types,
    model_engines
  )
  scores <- lapply(components, `[[`, "score")
  breads <- lapply(components, `[[`, "bread")
  parameters <- lapply(components, `[[`, "parameters")
  local_names <- lapply(parameters, names)

  for (i in which(model_types %in% c(
    "panel_fe", "panel_logit_re", "panel_probit_re", "panel_poisson_re"
  ) | .suest_is_glmm(model_types)))
    models[[i]] <- .suest_set_parameters(
      models[[i]],
      parameters[[i]],
      model_types[i],
      model_engines[[i]]
    )

  for (i in seq_along(models)) {
    U <- scores[[i]]
    B <- breads[[i]]
    mf <- model_frames[[i]]

    if (nrow(U) != nrow(mf))
      stop(
        paste0(
          "Each model must provide one score contribution per observation. ",
          "For nnet::multinom(), use the default summ = 0."
        ),
        call. = FALSE
      )

    if (nrow(B) != ncol(U) || ncol(B) != ncol(U))
      stop(
        sprintf(
          "The bread and score dimensions do not match for model '%s'.",
          model_names[i]
        ),
        call. = FALSE
      )

    if (!identical(colnames(U), local_names[[i]]))
      stop(
        sprintf(
          "The score parameter order does not match model '%s'.",
          model_names[i]
        ),
        call. = FALSE
      )

    if (nrow(U) <= 1L)
      stop("Each model must contain at least two observations.",
           call. = FALSE)
  }

  sources <- vapply(models, .suest_data_source, character(1))
  sample_rows <- lapply(model_frames, rownames)

  for (i in seq_along(sample_rows)) {
    if (is.null(sample_rows[[i]]))
      sample_rows[[i]] <- as.character(seq_len(nrow(model_frames[[i]])))
    if (is.null(observation_id) && anyDuplicated(sample_rows[[i]]))
      stop(
        sprintf(
          "Model '%s' has duplicated model-frame row names and cannot be aligned automatically.",
          model_names[i]
        ),
        call. = FALSE
      )
  }

  if (is.null(observation_id))
    .suest_check_row_alignment(model_frames, sources, sample_rows, model_names)

  sample_ids <- if (is.null(observation_id)) {
    NULL
  } else {
    .suest_observation_ids(
      models,
      model_frames,
      observation_id,
      model_names
    )
  }

  sample_keys <- if (is.null(sample_ids)) {
    Map(
      function(source, rows) paste(source, rows, sep = "\r"),
      sources,
      sample_rows
    )
  } else {
    lapply(sample_ids, `[[`, "keys")
  }

  union_keys <- unique(unlist(sample_keys, use.names = FALSE))
  overlap_keys <- Reduce(intersect, sample_keys)

  group_name <- if (is.null(cluster) && any(is_panel) &&
      all(model_types[is_panel] %in% multilevel_types))
    .suest_shared_group_name(models[is_panel], model_names[is_panel],
      required = !all(is_panel))
  cluster_info <- if (is.null(cluster) && any(is_panel)) {
    .suest_group_clusters(models, model_frames, is_panel, group_name, model_names)
  } else if (is.null(cluster)) {
    NULL
  } else {
    .suest_clusters(models, model_frames, cluster, model_names)
  }
  if (any(is_panel))
    .suest_validate_panel_cluster_nesting(
      models[is_panel], cluster_info[is_panel], model_names[is_panel]
    )
  union_clusters <- if (is.null(cluster_info)) {
    NULL
  } else {
    .suest_union_clusters(
      cluster_info,
      sample_keys,
      union_keys,
      model_names
    )
  }

  if (identical(weight_type, "pweight"))
    .suest_validate_overlap_weights(
      model_weights,
      sample_keys,
      model_names
    )

  n_model <- vapply(scores, nrow, integer(1))
  n_union <- length(union_keys)
  n_overlap <- length(overlap_keys)

  n_clusters <- if (is.null(union_clusters)) {
    NULL
  } else {
    length(unique(union_clusters))
  }
  if (!is.null(n_clusters) && n_clusters < 2L)
    stop("The joint system must contain at least two clusters.", call. = FALSE)

  corrected_scores <- if (!is.null(union_clusters)) {
    lapply(seq_along(scores), function(i) {
      # Stata's specialized 2SLS route stacks raw coefficient influences.
      # Ordinary suest models use the system-level G/(G-1) correction.
      # Stata's native xtpoisson, fe vce(robust) applies no G/(G-1) factor.
      correction <- if (model_types[i] %in% c("ivreg", "panel_poisson_fe")) {
        1
      } else if (model_types[i] %in% c(
                   "panel_ml", "panel_gee", "panel_logit_re",
                   "panel_probit_re", "panel_poisson_re",
                   "re_ologit", "re_oprobit"
                 ) || .suest_is_glmm(model_types[i])) {
        model_clusters <- length(unique(cluster_info[[i]]$keys))
        if (model_clusters < 2L)
          stop(
            sprintf(
              "Panel model '%s' must contain at least two clusters.",
              model_names[i]
            ),
            call. = FALSE
          )
        model_clusters / (model_clusters - 1)
      } else if (model_types[i] %in% panel_types) {
        model_clusters <- length(unique(cluster_info[[i]]$keys))
        if (model_clusters < 2L)
          stop(
            sprintf(
              "Panel model '%s' must contain at least two clusters.",
              model_names[i]
            ),
            call. = FALSE
          )
        k <- components[[i]]$cluster_df_k
        df_n <- components[[i]]$cluster_df_n
        if (is.null(df_n))
          df_n <- n_model[i]
        if (df_n <= k)
          stop(
            sprintf(
              "Panel model '%s' has insufficient residual degrees of freedom.",
              model_names[i]
            ),
            call. = FALSE
          )
        (model_clusters / (model_clusters - 1)) *
          ((df_n - 1) / (df_n - k))
      } else if (heterogeneous) {
        # Each ordinary model's block equals its own fit with
        # vce(cluster group), including regress's (N-1)/(N-k).
        model_clusters <- length(unique(cluster_info[[i]]$keys))
        if (model_clusters < 2L)
          stop(
            sprintf(
              "Model '%s' must contain at least two clusters.",
              model_names[i]
            ),
            call. = FALSE
          )
        small_sample <- if (model_types[i] == "lm") {
          (n_model[i] - 1) / (n_model[i] - models[[i]]$rank)
        } else {
          1
        }
        (model_clusters / (model_clusters - 1)) * small_sample
      } else {
        n_clusters / (n_clusters - 1)
      }
      scores[[i]] * sqrt(correction)
    })
  } else if (identical(weight_type, "pweight")) {
    lapply(
      scores,
      function(U) U * sqrt(n_union / (n_union - 1))
    )
  } else {
    Map(
      function(U, type) {
        # Stata's native robust 2SLS covariance is HC0. Its specialized
        # suest2 IV route stacks coefficient influence functions without the
        # N/(N-1) correction used by ordinary suest models.
        correction <- if (type == "ivreg") 1 else n_union / (n_union - 1)
        U * sqrt(correction)
      },
      scores,
      model_types
    )
  }

  aligned_scores <- Map(
    .suest_align_scores,
    corrected_scores,
    sample_keys,
    MoreArgs = list(union_keys = union_keys)
  )

  joint_breads <- Map(
    function(B, n) B * n_union / n,
    breads,
    n_model
  )

  U <- do.call(cbind, aligned_scores)
  B <- .suest_block_diag(joint_breads)
  meat_scores <- if (is.null(union_clusters)) {
    U
  } else {
    rowsum(U, group = union_clusters, reorder = FALSE)
  }
  meat <- crossprod(meat_scores) / n_union
  V <- (B %*% meat %*% t(B)) / n_union
  V <- (V + t(V)) / 2

  joint_names <- unlist(Map(
    function(model, term) paste0(model, "::", term),
    rep(model_names, lengths(local_names)),
    unlist(local_names, use.names = FALSE)
  ), use.names = FALSE)

  b <- unlist(parameters, use.names = FALSE)
  names(b) <- joint_names
  dimnames(V) <- list(joint_names, joint_names)

  p <- lengths(local_names)
  starts <- cumsum(c(1L, p[-length(p)]))
  index <- Map(
    function(start, size) start - 1L + seq_len(size),
    starts,
    p
  )
  names(index) <- model_names

  theta <- stats::setNames(rep(NA_real_, n_models), model_names)
  for (i in seq_along(models)) {
    if (model_types[i] == "negbin")
      theta[i] <- as.numeric(models[[i]]$theta)
  }
  if (all(is.na(theta)))
    theta <- NULL

  out <- list(
    models = models,
    model_frames = model_frames,
    raw_data = Map(.suest_raw_data, models, model_frames),
    model_weights = model_weights,
    weight_type = weight_type,
    coefficients = b,
    vcov = V,
    model_names = model_names,
    local_names = local_names,
    index = index,
    model_types = model_types,
    model_engines = model_engines,
    category_levels = category_levels,
    pair_key = pair$key,
    mixed_models = pair$mixed,
    comparison_scale = pair$scale,
    theta = theta,
    sample_sources = sources,
    sample_rows = sample_rows,
    sample_ids = if (is.null(sample_ids)) {
      NULL
    } else {
      lapply(sample_ids, `[[`, "values")
    },
    cluster_ids = if (is.null(cluster_info)) {
      NULL
    } else {
      lapply(cluster_info, `[[`, "values")
    },
    cluster_keys = union_clusters,
    n_clusters = n_clusters,
    sample_keys = sample_keys,
    nobs_models = stats::setNames(n_model, model_names),
    nobs_union = n_union,
    nobs_overlap = n_overlap,
    nobs_overlap_pairwise = .suest_pairwise_overlap(sample_keys, model_names),
    call = match.call()
  )
  class(out) <- "suest_model"
  out
}

#' @rdname suest
#' @export
coef.suest_model <- function(object, ...) object$coefficients

#' @rdname suest
#' @export
vcov.suest_model <- function(object, ...) object$vcov

#' @rdname suest
#' @export
nobs.suest_model <- function(object, ...) object$nobs_union

#' @rdname suest
#' @export
print.suest_model <- function(x, ...) {
  cat("Seemingly Unrelated Estimation\n")
  cat("Models:", paste(x$model_names, collapse = " + "), "\n")
  cat(
    "Model types:",
    paste(
      paste0(x$model_names, "=", unname(x$model_types)),
      collapse = ", "
    ),
    "\n"
  )
  if (!is.null(x$model_engines)) {
    cat(
      "Model engines:",
      paste(
        paste0(x$model_names, "=", unname(x$model_engines)),
        collapse = ", "
      ),
      "\n"
    )
  }
  cat("Comparison scale:", x$comparison_scale, "\n")
  if (!is.null(x$weight_type))
    cat("Weight type:", x$weight_type, "\n")
  if (!is.null(x$n_clusters))
    cat("Clusters:", format(x$n_clusters, big.mark = ","), "\n")

  if (!is.null(x$theta)) {
    keep <- !is.na(x$theta)
    cat(
      "Theta:",
      paste(
        paste0(
          names(x$theta)[keep],
          "=",
          format(x$theta[keep], digits = 5)
        ),
        collapse = ", "
      ),
      "\n"
    )
  }

  cat(
    "Observations:",
    paste(
      paste0(names(x$nobs_models), "=", format(x$nobs_models, big.mark = ",")),
      collapse = ", "
    ),
    "\n"
  )
  overlap_label <- if (length(x$models) == 2L) {
    "Overlapping observations:"
  } else {
    "Observations common to all models:"
  }
  cat(overlap_label, format(x$nobs_overlap, big.mark = ","), "\n")
  cat("Union observations:", format(x$nobs_union, big.mark = ","), "\n")
  cat("Parameters:", length(x$coefficients), "\n")
  invisible(x)
}

# Create model-specific newdata using each model's own estimation sample.
# This is useful when models were estimated on different samples and effects
# should be averaged separately within those samples.
