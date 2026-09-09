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

.suest_model_adapter <- function(model) {
  if (inherits(model, "svyglm"))
    stop("Survey models require the dedicated survey_design route in suest().", call. = FALSE)
  .suest_reject_adjusted_glm(model)

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
