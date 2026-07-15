# Shared helpers for the SUEST acceptance-test suite

.suest_test_state <- new.env(parent = emptyenv())
.suest_test_state$results <- list()
.suest_test_state$root <- NULL

test_initialize <- function(root) {
  .suest_test_state$root <- normalizePath(root, mustWork = TRUE)
  .suest_test_state$results <- list()
  invisible(TRUE)
}

test_record <- function(name, status, message = "") {
  .suest_test_state$results[[length(.suest_test_state$results) + 1L]] <-
    data.frame(
      test = name,
      status = status,
      message = message,
      stringsAsFactors = FALSE
    )
  invisible(TRUE)
}

test_case <- function(name, expr) {
  cat("\n", paste(rep("=", 78), collapse = ""), "\n", sep = "")
  cat(name, "\n")
  cat(paste(rep("=", 78), collapse = ""), "\n", sep = "")

  value <- tryCatch(
    force(expr),
    error = function(e) {
      test_record(name, "FAIL", conditionMessage(e))
      cat("FAIL:", conditionMessage(e), "\n")
      return(structure(
        list(message = conditionMessage(e)),
        class = "suest_test_failure"
      ))
    }
  )

  if (!inherits(value, "suest_test_failure")) {
    test_record(name, "PASS", "")
    cat("PASS\n")
    if (!is.null(value))
      print(value)
  }

  invisible(value)
}

expect_true <- function(x, message = "Expected condition to be TRUE.") {
  if (length(x) != 1L || is.na(x) || !isTRUE(x))
    stop(message, call. = FALSE)
  invisible(TRUE)
}

expect_equal <- function(actual, expected, tolerance = 1e-8,
                         message = NULL) {
  ok <- isTRUE(all.equal(
    actual,
    expected,
    tolerance = tolerance,
    check.attributes = FALSE
  ))

  if (!ok) {
    if (is.null(message)) {
      message <- paste0(
        "Objects are not equal within tolerance ", tolerance, "."
      )
    }
    stop(message, call. = FALSE)
  }

  invisible(TRUE)
}

expect_near <- function(actual, expected, tolerance,
                        label = "value") {
  actual <- as.numeric(actual)
  expected <- as.numeric(expected)

  if (length(actual) != length(expected))
    stop(
      sprintf(
        "%s length mismatch: actual=%s, expected=%s.",
        label,
        length(actual),
        length(expected)
      ),
      call. = FALSE
    )

  error <- abs(actual - expected)
  if (any(!is.finite(error)) || any(error > tolerance)) {
    bad <- which(!is.finite(error) | error > tolerance)
    details <- paste(
      sprintf(
        "%s[%s]: actual=%0.8f, expected=%0.8f, error=%0.8f",
        label,
        bad,
        actual[bad],
        expected[bad],
        error[bad]
      ),
      collapse = "; "
    )
    stop(
      paste0(
        "Benchmark tolerance exceeded (", tolerance, "): ",
        details
      ),
      call. = FALSE
    )
  }

  invisible(TRUE)
}

expect_error <- function(expr, pattern = NULL) {
  error <- tryCatch(
    {
      force(expr)
      NULL
    },
    error = function(e) e
  )

  if (is.null(error))
    stop("Expected an error, but none was raised.", call. = FALSE)

  if (!is.null(pattern) &&
      !grepl(pattern, conditionMessage(error), ignore.case = TRUE))
    stop(
      paste0(
        "Error did not match '", pattern, "': ",
        conditionMessage(error)
      ),
      call. = FALSE
    )

  invisible(conditionMessage(error))
}

test_summary <- function() {
  if (length(.suest_test_state$results) == 0L) {
    out <- data.frame(
      test = character(),
      status = character(),
      message = character()
    )
  } else {
    out <- do.call(rbind, .suest_test_state$results)
  }

  cat("\n\n", paste(rep("=", 78), collapse = ""), "\n", sep = "")
  cat("TEST SUMMARY\n")
  cat(paste(rep("=", 78), collapse = ""), "\n", sep = "")
  print(out, row.names = FALSE)

  cat("\nPassed:", sum(out$status == "PASS"), "\n")
  cat("Failed:", sum(out$status == "FAIL"), "\n")

  invisible(out)
}

require_test_packages <- function() {
  packages <- c(
    "brglm2",
    "glm2",
    "haven",
    "marginaleffects",
    "ordinal",
    "sandwich",
    "MASS",
    "nnet"
  )

  missing <- packages[
    !vapply(packages, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))
  ]

  if (length(missing))
    stop(
      "Install required packages: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )

  invisible(TRUE)
}

test_cache_dir <- function() {
  path <- file.path(.suest_test_state$root, "data-cache")
  dir.create(path, showWarnings = FALSE, recursive = TRUE)
  path
}

load_stata_data <- function(filename) {
  url <- paste0(
    "https://tdmize.github.io/data/data/",
    filename
  )
  destination <- file.path(test_cache_dir(), filename)

  if (!file.exists(destination)) {
    cat("Downloading ", filename, "...\n", sep = "")
    utils::download.file(
      url,
      destination,
      mode = "wb",
      quiet = TRUE
    )
  }

  haven::zap_labels(haven::read_dta(destination))
}

complete_data <- function(data, variables) {
  data[stats::complete.cases(data[variables]), , drop = FALSE]
}

factor_numeric <- function(x, ordered = FALSE) {
  values <- sort(unique(x[!is.na(x)]))
  factor(
    x,
    levels = values,
    ordered = ordered
  )
}

factorize <- function(data, variables, ordered = FALSE) {
  for (variable in variables)
    data[[variable]] <- factor_numeric(data[[variable]], ordered = ordered)
  data
}

forward_change <- function(amount) {
  force(amount)
  function(x) {
    data.frame(
      lo = x,
      hi = x + amount
    )
  }
}

vcov_block <- function(object, model = 1L) {
  index <- object$index[[model]]
  out <- object$vcov[index, index, drop = FALSE]
  names <- object$local_names[[model]]
  dimnames(out) <- list(names, names)
  out
}

robust_vcov_direct <- function(model) {
  U <- sandwich::estfun(model)
  B <- sandwich::bread(model)
  n <- nrow(U)
  B %*% crossprod(U) %*% B / n^2 * n / (n - 1)
}

offdiag_vcov <- function(object) {
  object$vcov[
    object$index[[1]],
    object$index[[2]],
    drop = FALSE
  ]
}

select_effects <- function(x, term, groups = NULL) {
  keep <- as.character(x$term) == term
  if (!is.null(groups))
    keep <- keep & as.character(x$group) %in% groups
  x[keep, , drop = FALSE]
}

pair_by_term_hypothesis <- function(model1, model2) {
  force(model1)
  force(model2)

  function(x) {
    term <- as.character(x$term)
    contrast <- if ("contrast" %in% names(x)) {
      as.character(x$contrast)
    } else {
      rep("", nrow(x))
    }
    group <- as.character(x$group)
    keys <- unique(paste(term, contrast, sep = "\r"))

    output <- lapply(keys, function(key) {
      rows <- paste(term, contrast, sep = "\r") == key
      i1 <- which(rows & group == model1)
      i2 <- which(rows & group == model2)

      if (length(i1) != 1L || length(i2) != 1L)
        return(NULL)

      data.frame(
        term = paste(
          term[i1],
          contrast[i1],
          paste0(model1, " - ", model2),
          sep = ": "
        ),
        estimate = x$estimate[i1] - x$estimate[i2]
      )
    })

    output <- Filter(Negate(is.null), output)
    if (!length(output))
      stop("No matched effects were found for the two models.")

    do.call(rbind, output)
  }
}

alternative_predictor_hypothesis <- function(x) {
  group <- as.character(x$group)
  term <- as.character(x$term)

  behavior <- x$estimate[
    term == "sexbehav" & group == "Behavior"
  ]
  identity <- x$estimate[
    term == "sexident" & group == "Identity"
  ]

  if (length(behavior) != 2L || length(identity) != 2L)
    stop("Expected two nonreference effects for each sexuality measure.")

  data.frame(
    term = c(
      "Behavior bisexual - Identity bisexual",
      "Behavior gay - Identity gay",
      "(Behavior gay - bisexual) - (Identity gay - bisexual)"
    ),
    estimate = c(
      behavior[1] - identity[1],
      behavior[2] - identity[2],
      (behavior[2] - behavior[1]) -
        (identity[2] - identity[1])
    )
  )
}

category_difference_hypothesis <- function(model1, model2, category) {
  force(model1)
  force(model2)
  force(category)

  function(x) {
    group <- as.character(x$group)
    key1 <- paste0(model1, "::", category)
    key2 <- paste0(model2, "::", category)
    i1 <- match(key1, group)
    i2 <- match(key2, group)

    if (is.na(i1) || is.na(i2))
      stop("Requested model-category combination was not found.")

    data.frame(
      term = paste0(key1, " - ", key2),
      estimate = x$estimate[i1] - x$estimate[i2]
    )
  }
}

model_matrix_atmeans <- function(model, age) {
  mf <- stats::model.frame(model)
  tt <- stats::delete.response(stats::terms(model))
  X <- stats::model.matrix(
    tt,
    data = mf,
    contrasts.arg = model$contrasts
  )
  x <- colMeans(X)

  if ("age" %in% names(x))
    x["age"] <- age
  if ("I(age^2)" %in% names(x))
    x["I(age^2)"] <- age^2

  x
}

polr_prob_atmeans <- function(model, age) {
  x <- model_matrix_atmeans(model, age)
  beta <- model$coefficients
  x <- x[names(beta)]
  eta <- sum(x * beta)
  cumulative <- stats::plogis(model$zeta - eta)
  c(
    cumulative[1],
    diff(cumulative),
    1 - cumulative[length(cumulative)]
  )
}

multinom_prob_atmeans <- function(model, age) {
  x <- model_matrix_atmeans(model, age)
  cf <- stats::coef(model)
  x <- x[colnames(cf)]
  eta <- c(0, as.numeric(cf %*% x))
  eta <- eta - max(eta)
  probability <- exp(eta) / sum(exp(eta))
  names(probability) <- model$lev
  probability
}

polr_probability_jacobian_atmeans <- function(model, age) {
  x <- model_matrix_atmeans(model, age)
  beta <- model$coefficients
  zeta <- model$zeta
  x <- x[names(beta)]

  eta <- sum(x * beta)
  cumulative <- stats::plogis(zeta - eta)
  density <- cumulative * (1 - cumulative)
  categories <- length(zeta) + 1L

  probability <- c(
    cumulative[1],
    diff(cumulative),
    1 - cumulative[length(cumulative)]
  )

  derivative_eta <- c(
    -density[1],
    density[-length(density)] - density[-1L],
    density[length(density)]
  )
  J_beta <- outer(derivative_eta, x)

  J_zeta <- matrix(
    0,
    nrow = categories,
    ncol = length(zeta)
  )
  J_zeta[1L, 1L] <- density[1L]

  if (categories > 2L) {
    for (category in 2L:(categories - 1L)) {
      J_zeta[category, category - 1L] <- -density[category - 1L]
      J_zeta[category, category] <- density[category]
    }
  }

  J_zeta[categories, length(zeta)] <- -density[length(zeta)]
  J <- cbind(J_beta, J_zeta)
  colnames(J) <- c(names(beta), names(zeta))

  list(
    probability = probability,
    jacobian = J
  )
}

multinom_probability_jacobian_atmeans <- function(model, age) {
  x <- model_matrix_atmeans(model, age)
  coefficients <- stats::coef(model)
  x <- x[colnames(coefficients)]

  eta <- c(0, as.numeric(coefficients %*% x))
  eta <- eta - max(eta)
  probability <- exp(eta) / sum(exp(eta))

  categories <- length(probability)
  regressors <- length(x)
  J <- matrix(
    0,
    nrow = categories,
    ncol = (categories - 1L) * regressors
  )

  for (outcome in 2L:categories) {
    derivative <- probability *
      ((seq_len(categories) == outcome) - probability[outcome])
    columns <- (outcome - 2L) * regressors + seq_len(regressors)
    J[, columns] <- outer(derivative, x)
  }

  colnames(J) <- as.vector(t(outer(
    model$lev[-1L],
    colnames(coefficients),
    paste,
    sep = ":"
  )))

  list(
    probability = probability,
    jacobian = J
  )
}

paper65_results <- function(object, age_lo = 20, age_hi = 30) {
  ordered_lo <- polr_probability_jacobian_atmeans(
    object$models[[1]],
    age_lo
  )
  ordered_hi <- polr_probability_jacobian_atmeans(
    object$models[[1]],
    age_hi
  )
  nominal_lo <- multinom_probability_jacobian_atmeans(
    object$models[[2]],
    age_lo
  )
  nominal_hi <- multinom_probability_jacobian_atmeans(
    object$models[[2]],
    age_hi
  )

  ordered_effect <- ordered_hi$probability - ordered_lo$probability
  nominal_effect <- nominal_hi$probability - nominal_lo$probability
  ordered_J <- ordered_hi$jacobian - ordered_lo$jacobian
  nominal_J <- nominal_hi$jacobian - nominal_lo$jacobian

  categories <- object$models[[1]]$lev
  K <- length(categories)
  P <- length(stats::coef(object))

  J <- matrix(0, nrow = 2L * K, ncol = P)
  J[seq_len(K), object$index[[1]]] <- ordered_J
  J[K + seq_len(K), object$index[[2]]] <- nominal_J

  estimate <- c(ordered_effect, nominal_effect)
  V <- J %*% stats::vcov(object) %*% t(J)
  V <- (V + t(V)) / 2

  labels <- c(
    paste0("Ordered::", categories),
    paste0("Multinomial::", categories)
  )
  names(estimate) <- labels
  dimnames(V) <- list(labels, labels)

  C <- cbind(diag(K), -diag(K))
  difference <- as.numeric(C %*% estimate)
  difference_vcov <- C %*% V %*% t(C)
  difference_vcov <- (difference_vcov + t(difference_vcov)) / 2

  list(
    effects = data.frame(
      group = labels,
      estimate = estimate,
      std.error = sqrt(diag(V)),
      row.names = NULL
    ),
    vcov = V,
    differences = data.frame(
      category = categories,
      estimate = difference,
      std.error = sqrt(diag(difference_vcov)),
      row.names = NULL
    ),
    difference_vcov = difference_vcov,
    jacobian = J
  )
}

check_basic_invariants <- function(object) {
  V <- stats::vcov(object)
  b <- stats::coef(object)

  expect_true(
    length(b) == nrow(V) && nrow(V) == ncol(V),
    "Coefficient and covariance dimensions do not agree."
  )
  expect_near(
    V,
    t(V),
    tolerance = 1e-8,
    label = "covariance symmetry"
  )
  expect_true(
    all(is.finite(V)),
    "The joint covariance contains nonfinite values."
  )

  restored <- suest:::set_coef.suest_model(object, suest:::get_coef.suest_model(object))
  expect_equal(
    suest:::get_coef.suest_model(restored),
    suest:::get_coef.suest_model(object),
    tolerance = 0
  )

  invisible(TRUE)
}
