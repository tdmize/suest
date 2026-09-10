# Internal model-adapter layer ---------------------------------------------
#
# Each adapter identifies the statistical model type separately from the
# package/function used to fit it. This lets multiple fitting engines share
# the same SUEST covariance and marginaleffects workflow.

.suest_call_head <- function(model) {
  call <- model$call
  if (is.null(call) || length(call) == 0L)
    return("")

  paste(deparse(call[[1L]], width.cutoff = 500L), collapse = "")
}

.suest_engine_name <- function(engine) {
  engine <- unname(as.character(engine))

  if (length(engine) != 1L || is.na(engine) || engine == "")
    stop("Internal error: invalid model-engine label.", call. = FALSE)

  engine
}

.suest_reject_adjusted_glm <- function(model) {
  adjusted_classes <- c(
    "brglmFit",
    "brglm",
    "brnb",
    "brmultinom",
    "bracl",
    "mdyplFit"
  )

  method_text <- if (!is.null(model$call$method)) {
    paste(
      deparse(model$call$method, width.cutoff = 500L),
      collapse = ""
    )
  } else {
    ""
  }

  adjusted <- any(vapply(
    adjusted_classes,
    function(class_name) inherits(model, class_name),
    logical(1)
  )) || grepl(
    "brglm|mdypl|bias.?reduc|firth",
    method_text,
    ignore.case = TRUE
  )

  if (adjusted)
    stop(
      paste0(
        "Bias-reduced, adjusted-score, Firth, and penalized GLM fits are ",
        "not supported. Their estimating equations differ from ordinary ",
        "maximum-likelihood GLM scores."
      ),
      call. = FALSE
    )

  invisible(TRUE)
}

.suest_hetprob_link <- function(model) {
  link_call <- model$call$link
  if (is.null(link_call))
    return("probit")

  if (!is.character(link_call) || length(link_call) != 1L)
    stop(
      paste0(
        "Refit Rchoice::hetprob() with link written explicitly as ",
        "'probit' or 'logit'."
      ),
      call. = FALSE
    )

  tolower(link_call)
}

.suest_is_pglm <- function(model) {
  inherits(model, "maxLik") &&
    .suest_call_head(model) %in% c("pglm", "pglm::pglm") &&
    is.function(model$objectiveFn) &&
    is.environment(environment(model$objectiveFn))
}

.suest_pglm_environment <- function(model) {
  if (!.suest_is_pglm(model))
    stop("Internal error: invalid pglm model object.", call. = FALSE)

  fit_environment <- environment(model$objectiveFn)
  required <- c(
    "data", "effect", "family", "formula", "id", "link", "model",
    "X", "y"
  )
  if (!all(vapply(
    required,
    exists,
    logical(1),
    envir = fit_environment,
    inherits = FALSE
  )))
    stop(
      "The pglm fit does not retain the data needed to reconstruct its likelihood.",
      call. = FALSE
    )

  fit_environment
}

.suest_pglm_formula <- function(model) {
  get("formula", .suest_pglm_environment(model), inherits = FALSE)
}

.suest_pglm_xlevels <- function(model) {
  frame <- get(
    "data",
    .suest_pglm_environment(model),
    inherits = FALSE
  )
  factors <- vapply(frame, is.factor, logical(1))
  lapply(frame[factors], levels)
}

.suest_glmmtmb_xlevels <- function(model) {
  frame <- stats::model.frame(model)
  variables <- all.vars(stats::delete.response(stats::terms(model)))
  factors <- names(frame) %in% variables & vapply(frame, is.factor, logical(1))
  lapply(frame[factors], levels)
}

.suest_normal_quadrature <- function(points) {
  points <- as.integer(points)
  if (length(points) != 1L || is.na(points) || points < 1L)
    stop("Internal error: invalid normal-quadrature order.", call. = FALSE)

  jacobi <- matrix(0, nrow = points, ncol = points)
  if (points > 1L) {
    off_diagonal <- sqrt(seq_len(points - 1L))
    jacobi[cbind(seq_len(points - 1L), 2:points)] <- off_diagonal
    jacobi[cbind(2:points, seq_len(points - 1L))] <- off_diagonal
  }
  decomposition <- eigen(jacobi, symmetric = TRUE)
  order <- order(decomposition$values)
  list(
    nodes = decomposition$values[order],
    weights = decomposition$vectors[1L, order]^2
  )
}

.suest_glmmtmb_random_info <- function(model) {
  random <- model$modelInfo$reTrms$cond
  if (!is.list(random) || length(random$cnms) != 1L ||
      length(random$flist) != 1L ||
      !identical(unname(random$cnms[[1L]]), "(Intercept)"))
    stop(
      paste0(
        "glmmTMB support is currently limited to one grouping variable ",
        "with one conditional random intercept."
      ),
      call. = FALSE
    )

  structure <- model$modelInfo$reStruc$condReStruc
  if (length(structure) != 1L || structure[[1L]]$blockSize != 1L ||
      structure[[1L]]$blockNumTheta != 1L)
    stop(
      "Structured glmmTMB random effects are not supported.",
      call. = FALSE
    )

  list(
    name = names(random$flist)[1L],
    id = droplevels(random$flist[[1L]])
  )
}

.suest_model_adapter <- function(model) {
  if (inherits(model, "svyglm"))
    stop("Survey models require the dedicated survey_design route in suest().", call. = FALSE)
  .suest_reject_adjusted_glm(model)

  if (inherits(model, "glmmTMB")) {
    family <- model$modelInfo$family
    model_type <- if (identical(family$family, "binomial") &&
                      identical(family$link, "logit")) {
      "glmm_logit_ri"
    } else if (identical(family$family, "poisson") &&
               identical(family$link, "log")) {
      "glmm_poisson_ri"
    } else {
      stop(
        paste0(
          "glmmTMB support is currently limited to binomial-logit and ",
          "Poisson-log models with one conditional random intercept."
        ),
        call. = FALSE
      )
    }

    random <- .suest_glmmtmb_random_info(model)
    if (length(model$modelInfo$reStruc$ziReStruc) ||
        length(model$modelInfo$reStruc$dispReStruc) ||
        length(glmmTMB::fixef(model)$zi) ||
        length(glmmTMB::fixef(model)$disp))
      stop(
        "Zero-inflation and dispersion submodels are not supported for glmmTMB models.",
        call. = FALSE
      )
    if (isTRUE(model$modelInfo$REML))
      stop("glmmTMB models must use maximum likelihood.", call. = FALSE)
    if (!is.null(model$modelInfo$priors) && nrow(model$modelInfo$priors))
      stop("glmmTMB models with priors are not supported.", call. = FALSE)

    fit_data <- model$obj$env$data
    if (!is.null(model$call$weights) ||
        length(fit_data$weights) != nrow(model$frame) ||
        any(!is.finite(fit_data$weights)) || any(fit_data$weights != 1))
      stop("Weighted glmmTMB models are not supported.", call. = FALSE)
    if (length(fit_data$offset) != nrow(model$frame) ||
        any(!is.finite(fit_data$offset)) || any(fit_data$offset != 0))
      stop("Offsets are not supported for glmmTMB models.", call. = FALSE)
    if (model_type == "glmm_logit_ri" &&
        (length(fit_data$size) != nrow(model$frame) ||
         any(fit_data$size != 1)))
      stop(
        "Grouped-binomial glmmTMB responses are not supported; use Bernoulli 0/1 outcomes.",
        call. = FALSE
      )
    if (length(random$id) != nrow(model$frame) || anyNA(random$id))
      stop(
        "The glmmTMB grouping variable does not align with its estimation sample.",
        call. = FALSE
      )
    if (!identical(model$fit$convergence, 0L) || !isTRUE(model$sdr$pdHess))
      stop(
        "The glmmTMB model must converge with a positive-definite Hessian.",
        call. = FALSE
      )

    theta <- model$fit$par[names(model$fit$par) == "theta"]
    if (length(theta) != 1L || !is.finite(theta) ||
        exp(theta) <= sqrt(.Machine$double.eps))
      stop(
        paste0(
          "The glmmTMB random-intercept standard deviation must be finite ",
          "and strictly away from its zero boundary."
        ),
        call. = FALSE
      )

    if (utils::packageVersion("glmmTMB") < "1.1.14")
      stop(
        paste0(
          "This glmmTMB version does not expose the full cluster scores and ",
          "covariance required by suest; use glmmTMB 1.1.14 or later."
        ),
        call. = FALSE
      )

    return(list(engine = "glmmTMB::glmmTMB", type = model_type))
  }

  if (inherits(model, "geeglm")) {
    family <- model$family$family
    link <- model$family$link
    valid <- (family == "gaussian" && link == "identity") ||
      (family == "binomial" && link %in% c("logit", "probit", "cloglog")) ||
      (family == "poisson" && link == "log")
    if (!valid || !model$corstr %in% c("independence", "exchangeable"))
      stop("GEE currently supports Gaussian identity, binary logit/probit/cloglog, and Poisson log with independence or exchangeable correlation.", call. = FALSE)
    if (!identical(model$geese$error, 0L))
      stop("The GEE fit did not converge.", call. = FALSE)
    if (any(model$prior.weights != 1) || any(model$geese$weights != 1))
      stop("Weighted GEE models are not supported.", call. = FALSE)
    if (length(model$geese$gamma) != 1L)
      stop("GEE requires a common scale parameter.", call. = FALSE)
    id <- model$id
    if (length(id) != nrow(model$model) || anyNA(id) ||
        anyDuplicated(rle(as.character(id))$values))
      stop("GEE IDs must identify contiguous panels in the retained model sample.", call. = FALSE)
    y <- stats::model.response(model$model)
    if (!is.numeric(y) || is.matrix(y) || any(!is.finite(y)) ||
        (family == "binomial" && !all(y %in% c(0, 1))) ||
        (family == "poisson" && any(y < 0 | abs(y - round(y)) > 1e-8)))
      stop("GEE outcomes must be numeric: binary 0/1, integer counts, or Gaussian responses.", call. = FALSE)
    return(list(engine = "geepack::geeglm", type = "panel_gee"))
  }

  if (.suest_is_pglm(model)) {
    fit_environment <- .suest_pglm_environment(model)
    family <- get("family", fit_environment, inherits = FALSE)
    link <- get("link", fit_environment, inherits = FALSE)
    panel_model <- get("model", fit_environment, inherits = FALSE)
    effect <- get("effect", fit_environment, inherits = FALSE)
    terms <- attr(get("data", fit_environment, inherits = FALSE), "terms")
    parameters <- stats::coef(model)

    binary <- identical(family, "binomial") && link %in% c("logit", "probit")
    poisson <- identical(family, "poisson") && identical(link, "log")
    if ((!binary && !poisson) || !identical(panel_model, "random") ||
        !identical(effect, "individual"))
      stop(
        paste0(
          "pglm support is currently limited to individual random-effects ",
          "binary logit/probit or Poisson log models."
        ),
        call. = FALSE
      )
    if (binary) {
      quadrature <- get("rn", fit_environment, inherits = FALSE)
      if (!is.list(quadrature) || length(quadrature$nodes) != 12L ||
          length(quadrature$weights) != 12L)
        stop(
          "Refit the pglm random-effects binary model with R = 12.",
          call. = FALSE
        )
    }
    if (poisson) {
      other <- get("other", fit_environment, inherits = FALSE)
      if (!identical(other, "sd"))
        stop(
          paste0(
            "Refit the pglm random-effects Poisson model with other = \"sd\" ",
            "so its ancillary parameter is the natural-scale gamma variance."
          ),
          call. = FALSE
        )
    }
    if (!is.null(terms) && length(attr(terms, "offset")))
      stop(
        "Offsets are not supported for pglm random-effects models.",
        call. = FALSE
      )
    if (!is.null(model$call$weights))
      stop(
        "Weighted pglm random-effects models are not supported.",
        call. = FALSE
      )
    if (!identical(model$code, 0L))
      stop(
        "The pglm random-effects model did not converge.",
        call. = FALSE
      )
    if (binary) {
      if (!"sigma" %in% names(parameters) ||
          !is.finite(unname(parameters["sigma"])) ||
          abs(unname(parameters["sigma"])) <= sqrt(.Machine$double.eps))
        stop(
          paste0(
            "The pglm random-intercept standard deviation must be finite and ",
            "strictly away from its zero boundary."
          ),
          call. = FALSE
        )
    } else if (!is.finite(parameters[length(parameters)]) ||
               parameters[length(parameters)] <= sqrt(.Machine$double.eps)) {
      stop(
        paste0(
          "The pglm gamma random-effect variance must be finite and strictly ",
          "away from its zero boundary."
        ),
        call. = FALSE
      )
    }
    if (!is.matrix(model$gradientObs) || !is.matrix(model$hessian))
      stop(
        "The pglm fit must retain its analytic scores and Hessian.",
        call. = FALSE
      )

    type <- if (identical(link, "logit")) {
      "panel_logit_re"
    } else if (identical(link, "probit")) {
      "panel_probit_re"
    } else {
      "panel_poisson_re"
    }
    return(list(engine = "pglm::pglm", type = type))
  }

  if (inherits(model, "mvProbit")) {
    if (!identical(model$nDep, 2L))
      stop(
        "Only bivariate mvProbit::mvProbit models are supported.",
        call. = FALSE
      )
    if (any(model$fixed))
      stop(
        "mvProbit models with fixed parameters are not supported.",
        call. = FALSE
      )
    observed_hessian <- identical(model$call$finalHessian, TRUE)
    if (!is.matrix(model$gradientObs) || !is.matrix(model$hessian) ||
        !observed_hessian)
      stop(
        paste0(
          "The mvProbit model must retain observation scores and a final ",
          "Hessian; refit with intGrad = TRUE and finalHessian = TRUE."
        ),
        call. = FALSE
      )

    return(list(engine = "mvProbit::mvProbit", type = "biprobit"))
  }

  if (inherits(model, "clm")) {
    if (!model$link %in% c("logit", "probit"))
      return(list(engine = "ordinal::clm", type = "unsupported"))

    if (!identical(model$threshold, "flexible"))
      stop(
        "ordinal::clm models currently require threshold = 'flexible'.",
        call. = FALSE
      )

    if (!is.null(model$S.terms) && length(model$S.terms) > 0L)
      stop(
        "ordinal::clm scale models are not supported yet.",
        call. = FALSE
      )

    if (!is.null(model$nom.terms) && length(model$nom.terms) > 0L)
      stop(
        "ordinal::clm nominal-effects models are not supported yet.",
        call. = FALSE
      )

    if (is.null(model$model))
      stop(
        "Refit ordinal::clm models with model = TRUE.",
        call. = FALSE
      )

    return(list(
      engine = "ordinal::clm",
      type = if (identical(model$link, "logit")) "ologit" else "oprobit"
    ))
  }

  if (inherits(model, "multinom")) {
    return(list(engine = "nnet::multinom", type = "multinom"))
  }

  if (inherits(model, "polr")) {
    type <- if (identical(model$method, "logistic")) {
      "ologit"
    } else if (identical(model$method, "probit")) {
      "oprobit"
    } else {
      "unsupported"
    }
    return(list(engine = "MASS::polr", type = type))
  }

  if (inherits(model, "negbin")) {
    type <- if (identical(unname(model$family$link), "log")) {
      "negbin"
    } else {
      "unsupported"
    }
    return(list(engine = "MASS::glm.nb", type = type))
  }

  if (inherits(model, "survreg")) {
    if (length(model$scale) != 1L)
      stop(
        "survival::survreg models with stratum-specific scales are not supported yet.",
        call. = FALSE
      )

    return(list(engine = "survival::survreg", type = "survreg"))
  }

  if (inherits(model, "betareg")) {
    distribution <- if (is.null(model$dist)) "beta" else model$dist
    mean_link <- model$link[[1L]]$name

    if (!identical(distribution, "beta"))
      stop(
        "Extended-support beta regressions are not supported yet.",
        call. = FALSE
      )
    if (!mean_link %in% c("logit", "probit", "cloglog", "loglog"))
      return(list(engine = "betareg::betareg", type = "unsupported"))

    return(list(engine = "betareg::betareg", type = "betareg"))
  }

  if (inherits(model, "zeroinfl")) {
    type <- switch(
      model$dist,
      "poisson" = "zip",
      "negbin" = "zinb",
      "unsupported"
    )
    return(list(engine = "pscl::zeroinfl", type = type))
  }

  if (inherits(model, "truncreg"))
    return(list(engine = "truncreg::truncreg", type = "truncreg"))

  if (inherits(model, "censReg"))
    return(list(engine = "censReg::censReg", type = "censreg"))

  if (inherits(model, "lme")) {
    random <- model$modelStruct$reStruct
    random_intercept <- length(random) == 1L &&
      identical(colnames(nlme::pdMatrix(random[[1L]])), "(Intercept)")

    if (!identical(model$method, "ML"))
      stop(
        "nlme::lme models must be fitted with method = \"ML\".",
        call. = FALSE
      )
    if (!random_intercept || length(names(model$groups)) != 1L)
      stop(
        "Only single-level random-intercept nlme::lme models are supported.",
        call. = FALSE
      )
    if (!is.null(model$modelStruct$corStruct) ||
        !is.null(model$modelStruct$varStruct))
      stop(
        "Correlated or heteroskedastic nlme::lme residuals are not supported.",
        call. = FALSE
      )

    return(list(engine = "nlme::lme", type = "panel_ml"))
  }

  if (inherits(model, "plm")) {
    panel_type <- if (identical(model$args$model, "within") &&
                      identical(model$args$effect, "individual")) {
      "panel_fe"
    } else if (identical(model$args$model, "between") &&
               identical(model$args$effect, "individual")) {
      "panel_be"
    } else if (identical(model$args$model, "random") &&
               identical(model$args$effect, "individual") &&
               (is.null(model$args$random.method) ||
                identical(model$args$random.method, "swar"))) {
      "panel_re"
    } else {
      "unsupported"
    }
    if (identical(panel_type, "unsupported"))
      return(list(engine = "plm::plm", type = "unsupported"))

    if (!is.null(model$call$weights))
      stop("Weighted plm fixed-effects models are not supported yet.", call. = FALSE)
    if (attr(stats::terms(model), "intercept") != 1L)
      stop(
        "plm fixed-effects models currently require an intercept in the formula.",
        call. = FALSE
      )

    return(list(engine = "plm::plm", type = panel_type))
  }

  if (inherits(model, "fixest")) {
    valid_iv <- identical(model$method, "feols") && isTRUE(model$is_iv) &&
      identical(as.integer(model$iv_stage), 2L)

    if (valid_iv && length(model$fixef_vars))
      stop(
        "Instrumental-variable fixest models with absorbed fixed effects are not supported yet.",
        call. = FALSE
      )

    return(list(
      engine = "fixest::feols",
      type = if (valid_iv) "ivreg" else "unsupported"
    ))
  }

  if (inherits(model, "ivpml"))
    return(list(engine = "Rchoice::ivpml", type = "ivprobit"))

  if (inherits(model, "hetprob")) {
    link <- .suest_hetprob_link(model)
    type <- switch(
      link,
      "probit" = "hetprobit",
      "logit" = "hetlogit",
      "unsupported"
    )
    return(list(engine = "Rchoice::hetprob", type = type))
  }

  if (inherits(model, "glm")) {
    family <- tolower(model$family$family)
    link <- tolower(model$family$link)
    key <- paste(family, link, sep = ":")
    type <- switch(
      key,
      "binomial:logit" = "logit",
      "binomial:probit" = "probit",
      "binomial:cloglog" = "cloglog",
      "quasibinomial:logit" = "fraclogit",
      "quasibinomial:probit" = "fracprobit",
      "quasibinomial:cloglog" = "fraccloglog",
      "quasibinomial:loglog" = "fracloglog",
      "poisson:log" = "poisson",
      "unsupported"
    )
    if (identical(type, "unsupported") &&
        link %in% c("identity", "log", "logit", "probit", "cloglog", "loglog"))
      type <- "glm"

    head <- .suest_call_head(model)
    engine <- if (head %in% c("glm2", "glm2::glm2")) {
      "glm2::glm2"
    } else {
      "stats::glm"
    }

    return(list(engine = engine, type = type))
  }

  if (inherits(model, "lm"))
    return(list(engine = "stats::lm", type = "lm"))

  list(engine = class(model)[1L], type = "unsupported")
}

.suest_model_type <- function(model) {
  .suest_model_adapter(model)$type
}

.suest_model_engine <- function(model) {
  .suest_model_adapter(model)$engine
}

.suest_model_frame <- function(model, engine) {
  engine <- .suest_engine_name(engine)
  if (identical(engine, "Rchoice::hetprob"))
    return(model$mf)
  if (identical(engine, "Rchoice::ivpml"))
    return(model$mf)

  if (identical(engine, "glmmTMB::glmmTMB")) {
    frame <- as.data.frame(stats::model.frame(model))
    formula <- stats::formula(model)
    original <- try(
      eval(model$call$data, envir = environment(formula)),
      silent = TRUE
    )
    if (!inherits(original, "try-error") && is.data.frame(original)) {
      index <- match(rownames(frame), rownames(original))
      if (!anyNA(index)) {
        extra <- setdiff(names(original), names(frame))
        for (variable in extra)
          frame[[variable]] <- original[[variable]][index]
      }
    }
    return(frame)
  }

  if (identical(engine, "mvProbit::mvProbit")) {
    fit_environment <- environment(model$objectiveFn)
    formula <- if (is.environment(fit_environment) &&
                   exists("formula", fit_environment, inherits = FALSE)) {
      get("formula", fit_environment, inherits = FALSE)
    } else {
      stats::formula(model)
    }
    data <- if (is.environment(fit_environment) &&
                exists("data", fit_environment, inherits = FALSE)) {
      get("data", fit_environment, inherits = FALSE)
    } else {
      try(
        eval(model$call$data, envir = environment(formula)),
        silent = TRUE
      )
    }
    if (inherits(data, "try-error") || !is.data.frame(data))
      stop(
        paste0(
          "Could not recover the original data for the mvProbit model. ",
          "Keep the data object used in mvProbit::mvProbit() available."
        ),
        call. = FALSE
      )

    frame <- if (is.environment(fit_environment) &&
                 exists("mf", fit_environment, inherits = FALSE)) {
      get("mf", fit_environment, inherits = FALSE)
    } else {
      stats::model.frame(
        formula,
        data = data,
        na.action = stats::na.pass
      )
    }
    if (nrow(frame) != model$nObs)
      stop(
        "Could not align the mvProbit estimation sample to its data.",
        call. = FALSE
      )
    index <- match(rownames(frame), rownames(data))
    if (anyNA(index))
      stop(
        "Could not align the mvProbit estimation sample to its data.",
        call. = FALSE
      )
    extra <- setdiff(names(data), names(frame))
    for (variable in extra)
      frame[[variable]] <- data[[variable]][index]
    return(frame)
  }

  if (identical(engine, "pglm::pglm")) {
    fit_environment <- .suest_pglm_environment(model)
    retained <- get("data", fit_environment, inherits = FALSE)
    frame <- as.data.frame(retained)
    attr(frame, "terms") <- attr(retained, "terms")
    rownames(frame) <- rownames(retained)

    index <- attr(retained, "index")
    if (is.data.frame(index) && nrow(index) == nrow(frame)) {
      for (variable in names(index))
        if (!variable %in% names(frame))
          frame[[variable]] <- index[[variable]]
    }

    formula <- .suest_pglm_formula(model)
    original <- try(
      eval(model$call$data, envir = environment(formula)),
      silent = TRUE
    )
    if (!inherits(original, "try-error") && is.data.frame(original)) {
      original_index <- match(rownames(frame), rownames(original))
      if (!anyNA(original_index)) {
        extra <- setdiff(names(original), names(frame))
        for (variable in extra)
          frame[[variable]] <- original[[variable]][original_index]
      }
    }

    return(frame)
  }

  if (identical(engine, "nlme::lme")) {
    data <- try(nlme::getData(model), silent = TRUE)
    used_rows <- rownames(model$groups)
    if (inherits(data, "try-error") || !is.data.frame(data) ||
        is.null(used_rows))
      stop(
        paste0(
          "Could not recover the estimation data for the nlme::lme model. ",
          "Refit with keep.data = TRUE."
        ),
        call. = FALSE
      )

    index <- match(used_rows, rownames(data))
    if (anyNA(index))
      stop(
        "Could not align the nlme::lme estimation sample to its data.",
        call. = FALSE
      )

    estimation_data <- data[index, , drop = FALSE]
    frame <- stats::model.frame(
      stats::formula(model),
      data = estimation_data,
      na.action = stats::na.pass
    )
    extra <- setdiff(names(estimation_data), names(frame))
    for (variable in extra)
      frame[[variable]] <- estimation_data[[variable]]
    rownames(frame) <- used_rows
    return(frame)
  }

  if (identical(engine, "fixest::feols")) {
    data_call <- model$call$data
    data_environment <- if (is.environment(model$call_env)) {
      model$call_env
    } else {
      environment(stats::formula(model))
    }
    data <- try(eval(data_call, envir = data_environment), silent = TRUE)

    if (inherits(data, "try-error") || !is.data.frame(data))
      stop(
        paste0(
          "Could not recover the original data for the fixest IV model. ",
          "Keep the data object used in fixest::feols() available."
        ),
        call. = FALSE
      )

    rows <- fixest::obs(model)
    if (length(rows) != stats::nobs(model) || anyNA(rows) ||
        any(rows < 1L | rows > nrow(data)))
      stop(
        "Could not align the fixest IV estimation sample to its original data.",
        call. = FALSE
      )

    estimation_data <- data[rows, , drop = FALSE]
    frame <- stats::model.frame(
      model$fml_all$linear,
      data = estimation_data,
      na.action = stats::na.pass
    )
    extra <- setdiff(names(estimation_data), names(frame))
    for (variable in extra)
      frame[[variable]] <- estimation_data[[variable]]
    rownames(frame) <- rownames(estimation_data)
    return(frame)
  }

  stats::model.frame(model)
}

.suest_category_levels <- function(model, engine) {
  engine <- .suest_engine_name(engine)
  if (identical(engine, "ordinal::clm")) {
    model$y.levels
  } else if (engine %in% c("MASS::polr", "nnet::multinom")) {
    model$lev
  } else {
    NULL
  }
}

.suest_set_parameters <- function(model, parameters, type, engine) {
  engine <- .suest_engine_name(engine)
  if (type %in% c("glmm_logit_ri", "glmm_poisson_ri")) {
    model$suest_parameters <- parameters
    return(model)
  }
  if (identical(engine, "ordinal::clm")) {
    n_alpha <- length(model$alpha)
    n_beta <- length(model$beta)

    model$coefficients <- parameters
    model$alpha <- parameters[seq_len(n_alpha)]

    if (n_beta > 0L) {
      model$beta <- parameters[n_alpha + seq_len(n_beta)]
    } else {
      model$beta <- numeric(0)
    }

    return(model)
  }

  if (type %in% c("ologit", "oprobit")) {
    k <- length(model$coefficients)
    q <- length(model$zeta)
    model$coefficients <- parameters[seq_len(k)]
    model$zeta <- parameters[k + seq_len(q)]
    return(model)
  }

  if (type == "multinom") {
    r <- length(model$vcoefnames)
    q <- length(model$lev) - 1L

    cf <- matrix(
      parameters,
      nrow = q,
      ncol = r,
      byrow = TRUE,
      dimnames = list(model$lev[-1L], model$vcoefnames)
    )

    W <- matrix(model$wts, nrow = model$n[3L], byrow = TRUE)
    W[-1L, 1L + seq_len(r)] <- cf
    model$wts <- as.vector(t(W))
    return(model)
  }

  if (type == "negbin") {
    beta_names <- setdiff(names(parameters), "ln_theta")
    model$coefficients <- parameters[beta_names]
    model$theta <- exp(unname(parameters["ln_theta"]))
    model$family <- MASS::negative.binomial(
      theta = model$theta,
      link = model$family$link
    )
    return(model)
  }

  if (type == "survreg") {
    scale_name <- "Log(scale)"
    beta_names <- setdiff(names(parameters), scale_name)
    model$coefficients <- parameters[beta_names]

    if (scale_name %in% names(parameters)) {
      model$scale <- exp(unname(parameters[scale_name]))
      if (!is.null(model$icoef) && scale_name %in% names(model$icoef))
        model$icoef[scale_name] <- unname(parameters[scale_name])
    }

    return(model)
  }

  if (type == "betareg") {
    component_lengths <- lengths(model$coefficients)
    starts <- cumsum(c(1L, component_lengths[-length(component_lengths)]))

    for (i in seq_along(model$coefficients)) {
      index <- starts[i] - 1L + seq_len(component_lengths[i])
      model$coefficients[[i]][] <- parameters[index]
    }

    return(model)
  }

  if (type %in% c("zip", "zinb")) {
    component_lengths <- lengths(model$coefficients)
    starts <- cumsum(c(1L, component_lengths[-length(component_lengths)]))

    for (i in seq_along(model$coefficients)) {
      index <- starts[i] - 1L + seq_len(component_lengths[i])
      model$coefficients[[i]][] <- parameters[index]
    }

    if (type == "zinb")
      model$theta <- exp(unname(parameters["ln_theta"]))

    return(model)
  }

  if (type == "censreg") {
    model$estimate[] <- parameters
    return(model)
  }

  if (type %in% c("hetprobit", "hetlogit")) {
    model$estimate[] <- parameters
    return(model)
  }

  if (type == "ivprobit") {
    names(parameters)[names(parameters) == "athrho"] <- "atanhrho"
    model$estimate[] <- parameters[names(model$estimate)]
    return(model)
  }

  if (type == "lm") {
    beta_names <- setdiff(names(parameters), "lnvar")
    model$coefficients <- parameters[beta_names]
    return(model)
  }

  if (type == "panel_fe") {
    beta_names <- setdiff(names(parameters), "(Intercept)")
    model$coefficients <- parameters[beta_names]
    model$suest_intercept <- unname(parameters["(Intercept)"])
    return(model)
  }

  if (type == "panel_ml") {
    beta_names <- names(model$coefficients$fixed)
    model$coefficients$fixed[] <- parameters[beta_names]
    return(model)
  }

  if (type %in% c("panel_logit_re", "panel_probit_re")) {
    model$estimate[] <- parameters[names(model$estimate)]
    return(model)
  }

  if (type == "panel_poisson_re") {
    beta_names <- names(model$estimate)[-length(model$estimate)]
    model$estimate[beta_names] <- parameters[beta_names]
    model$estimate[length(model$estimate)] <- parameters["alpha"]
    return(model)
  }

  if (type == "biprobit") {
    raw_names <- names(model$estimate)
    beta_names <- raw_names[-length(raw_names)]
    model$estimate[beta_names] <- parameters[beta_names]
    model$estimate[length(raw_names)] <- tanh(unname(parameters["athrho"]))
    return(model)
  }

  model$coefficients <- parameters
  model
}

.suest_predict_probabilities <- function(model, newdata, engine) {
  engine <- .suest_engine_name(engine)
  if (identical(engine, "ordinal::clm")) {
    response <- all.vars(stats::formula(model))[1L]
    prediction_data <- newdata
    prediction_data[[response]] <- NULL

    prediction <- stats::predict(
      model,
      newdata = prediction_data,
      type = "prob"
    )
    probability <- prediction$fit
  } else {
    probability <- stats::predict(
      model,
      newdata = newdata,
      type = "probs"
    )
  }

  if (is.null(dim(probability))) {
    category <- names(probability)
    probability <- matrix(
      as.numeric(probability),
      nrow = 1L,
      dimnames = list(NULL, category)
    )
  }

  probability
}

.suest_predict_values <- function(model, newdata, type, engine, model_type) {
  if (identical(engine, "survey::svyglm")) {
    terms <- stats::delete.response(stats::terms(model))
    mf <- stats::model.frame(terms, newdata, na.action = stats::na.pass,
      xlev = model$xlevels)
    X <- stats::model.matrix(terms, mf, contrasts.arg = model$contrasts)
    offset <- stats::model.offset(mf)
    if (is.null(offset)) offset <- rep(0, nrow(X))
    eta <- as.numeric(X %*% stats::coef(model)) + offset
    if (identical(type, "link") || identical(model_type, "survey_lm")) return(eta)
    return(model$family$linkinv(eta))
  }
  engine <- .suest_engine_name(engine)

  if (model_type %in% c("glmm_logit_ri", "glmm_poisson_ri")) {
    parameters <- if (!is.null(model$suest_parameters)) {
      model$suest_parameters
    } else {
      .suest_extract_parameters(model, model_type, engine)
    }
    beta_names <- setdiff(names(parameters), "log_sigma")
    terms <- stats::delete.response(stats::terms(model))
    prediction_frame <- stats::model.frame(
      terms, newdata, na.action = stats::na.pass,
      xlev = .suest_glmmtmb_xlevels(model)
    )
    X <- stats::model.matrix(
      terms, prediction_frame,
      contrasts.arg = model$modelInfo$contrasts$cond[
        intersect(names(model$modelInfo$contrasts$cond), names(prediction_frame))
      ]
    )
    missing <- setdiff(beta_names, colnames(X))
    if (length(missing))
      stop(
        "Unable to construct the glmmTMB random-intercept prediction matrix.",
        call. = FALSE
      )
    eta <- as.numeric(X[, beta_names, drop = FALSE] %*% parameters[beta_names])
    if (identical(type, "link"))
      return(eta)

    sigma <- exp(unname(parameters["log_sigma"]))
    if (model_type == "glmm_poisson_ri")
      return(exp(eta + sigma^2/2))

    quadrature <- .suest_normal_quadrature(20L)
    probabilities <- vapply(
      seq_along(quadrature$nodes),
      function(i) quadrature$weights[i] * stats::plogis(
        eta + sigma * quadrature$nodes[i]
      ),
      numeric(length(eta))
    )
    return(rowSums(probabilities))
  }

  if (model_type == "censreg") {
    beta <- stats::coef(model)
    beta <- beta[setdiff(names(beta), "logSigma")]
    terms <- stats::delete.response(stats::terms(model))
    X <- stats::model.matrix(terms, newdata)

    missing <- setdiff(names(beta), colnames(X))
    if (length(missing))
      stop(
        "Unable to construct the censored-regression prediction matrix.",
        call. = FALSE
      )

    return(as.numeric(X[, names(beta), drop = FALSE] %*% beta))
  }

  if (model_type %in% c("hetprobit", "hetlogit")) {
    model$call$link <- if (model_type == "hetprobit") "probit" else "logit"
    prediction_type <- if (identical(type, "link")) "xb" else "pr"
    return(stats::predict(
      model,
      newdata = newdata,
      type = prediction_type
    ))
  }

  if (model_type == "ivprobit") {
    prediction_type <- if (identical(type, "link")) "xb" else "pr"
    return(stats::predict(
      model,
      newdata = newdata,
      type = prediction_type,
      # Rchoice's asf = TRUE evaluates the residual-conditioned probability.
      # The structural response probability integrates out that disturbance.
      asf = FALSE
    ))
  }

  if (model_type == "panel_fe") {
    beta <- c(`(Intercept)` = model$suest_intercept, stats::coef(model))
    terms <- stats::delete.response(stats::terms(model))
    X <- stats::model.matrix(terms, newdata)
    missing <- setdiff(names(beta), colnames(X))
    if (length(missing))
      stop(
        "Unable to construct the panel fixed-effects prediction matrix.",
        call. = FALSE
      )
    return(as.numeric(X[, names(beta), drop = FALSE] %*% beta))
  }

  if (model_type == "panel_ml") {
    beta <- model$coefficients$fixed
    terms <- stats::delete.response(stats::terms(model))
    X <- stats::model.matrix(
      terms,
      newdata,
      contrasts.arg = model$contrasts
    )
    missing <- setdiff(names(beta), colnames(X))
    if (length(missing))
      stop(
        "Unable to construct the random-intercept ML prediction matrix.",
        call. = FALSE
      )
    return(as.numeric(X[, names(beta), drop = FALSE] %*% beta))
  }

  if (model_type %in% c("panel_logit_re", "panel_probit_re")) {
    fit_environment <- .suest_pglm_environment(model)
    training <- get("data", fit_environment, inherits = FALSE)
    training_matrix <- get("X", fit_environment, inherits = FALSE)
    terms <- stats::delete.response(attr(training, "terms"))
    prediction_frame <- stats::model.frame(
      terms,
      newdata,
      na.action = stats::na.pass,
      xlev = .suest_pglm_xlevels(model)
    )
    X <- stats::model.matrix(
      terms,
      prediction_frame,
      contrasts.arg = attr(training_matrix, "contrasts")
    )
    parameters <- stats::coef(model)
    beta_names <- setdiff(names(parameters), "sigma")
    missing <- setdiff(beta_names, colnames(X))
    if (length(missing))
      stop(
        "Unable to construct the pglm random-effects binary prediction matrix.",
        call. = FALSE
      )
    eta <- as.numeric(X[, beta_names, drop = FALSE] %*% parameters[beta_names])
    if (identical(type, "link"))
      return(eta)

    sigma <- unname(parameters["sigma"])
    quadrature <- get("rn", fit_environment, inherits = FALSE)
    inverse_link <- if (model_type == "panel_logit_re") {
      stats::plogis
    } else {
      stats::pnorm
    }
    probabilities <- vapply(
      seq_along(quadrature$nodes),
      function(i) quadrature$weights[i] * inverse_link(
        eta + sqrt(2) * sigma * quadrature$nodes[i]
      ),
      numeric(length(eta))
    )
    return(rowSums(probabilities) / sqrt(pi))
  }

  if (model_type == "panel_poisson_re") {
    fit_environment <- .suest_pglm_environment(model)
    training <- get("data", fit_environment, inherits = FALSE)
    training_matrix <- get("X", fit_environment, inherits = FALSE)
    terms <- stats::delete.response(attr(training, "terms"))
    prediction_frame <- stats::model.frame(
      terms, newdata, na.action = stats::na.pass,
      xlev = .suest_pglm_xlevels(model)
    )
    X <- stats::model.matrix(
      terms, prediction_frame,
      contrasts.arg = attr(training_matrix, "contrasts")
    )
    parameters <- stats::coef(model)
    beta_names <- names(parameters)[-length(parameters)]
    missing <- setdiff(beta_names, colnames(X))
    if (length(missing))
      stop(
        "Unable to construct the pglm random-effects Poisson prediction matrix.",
        call. = FALSE
      )
    eta <- as.numeric(X[, beta_names, drop = FALSE] %*% parameters[beta_names])
    return(if (identical(type, "link")) eta else exp(eta))
  }

  if (model_type == "biprobit") {
    if (identical(type, "link"))
      stop(
        "Bivariate probit models are supported on the response scale only.",
        call. = FALSE
      )
    if (!requireNamespace("mvtnorm", quietly = TRUE))
      stop(
        "Package 'mvtnorm' is required for bivariate-probit predictions.",
        call. = FALSE
      )

    formula <- stats::formula(model)
    terms <- stats::delete.response(stats::terms(formula))
    training <- .suest_model_frame(model, "mvProbit::mvProbit")
    training_matrix <- stats::model.matrix(terms, training)
    X <- stats::model.matrix(
      terms,
      newdata,
      contrasts.arg = attr(training_matrix, "contrasts")
    )
    k <- model$nReg
    beta1 <- model$estimate[seq_len(k)]
    beta2 <- model$estimate[k + seq_len(k)]
    if (ncol(X) != k)
      stop(
        "Unable to construct the bivariate-probit prediction matrix.",
        call. = FALSE
      )
    eta1 <- as.numeric(X %*% beta1)
    eta2 <- as.numeric(X %*% beta2)
    rho <- unname(model$estimate[length(model$estimate)])
    sigma <- matrix(c(1, rho, rho, 1), nrow = 2L)
    return(vapply(seq_along(eta1), function(i) {
      as.numeric(mvtnorm::pmvnorm(
        lower = c(-Inf, -Inf),
        upper = c(eta1[i], eta2[i]),
        mean = c(0, 0),
        sigma = sigma,
        algorithm = mvtnorm::TVPACK()
      ))
    }, numeric(1)))
  }

  prediction_type <- if (model_type == "lm") {
    "response"
  } else if (model_type == "survreg" && identical(type, "link")) {
    "linear"
  } else if (model_type %in% c("zip", "zinb") && identical(type, "link")) {
    stop(
      "Zero-inflated models are currently supported on the response scale only.",
      call. = FALSE
    )
  } else {
    type
  }

  stats::predict(
    model,
    newdata = newdata,
    type = prediction_type
  )
}
