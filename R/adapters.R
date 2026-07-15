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

.suest_model_adapter <- function(model) {
  .suest_reject_adjusted_glm(model)

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

    if (length(attr(model$terms, "offset")))
      stop("Offsets are not supported in this draft.", call. = FALSE)

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

  if (inherits(model, "glm")) {
    key <- paste(model$family$family, model$family$link, sep = ":")
    type <- switch(
      key,
      "binomial:logit" = "logit",
      "binomial:probit" = "probit",
      "poisson:log" = "poisson",
      "unsupported"
    )

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
  engine <- .suest_engine_name(engine)
  prediction_type <- if (model_type == "lm") "response" else type

  stats::predict(
    model,
    newdata = newdata,
    type = prediction_type
  )
}
