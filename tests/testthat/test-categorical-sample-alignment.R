alignment_data <- function() {
  set.seed(20261001)
  n <- 240L
  d <- data.frame(id = seq_len(n), f = factor(rep(LETTERS[1:4], length.out = n)),
    b = factor(rep(c("No", "Yes"), each = 3, length.out = n)), z = rnorm(n),
    s = factor(rep(c("Low", "High"), each = 5, length.out = n)))
  eta <- c(-.4, .3, .8, -.2)[d$f] + .4 * d$z
  p <- exp(cbind(0, eta, -.3 - eta, .2 + .3 * eta)); p <- p / rowSums(p)
  d$y <- factor(vapply(seq_len(n), function(i) sample(letters[1:4], 1, prob = p[i, ]), character(1)))
  d$yo <- ordered(cut(eta + rlogis(n), c(-Inf, -.5, .5, 1.5, Inf), labels = letters[1:4]))
  d
}

compare_categorical_samples <- function(fit, variables, by = NULL) {
  nd <- suest_newdata(fit)
  nd <- nd[rev(seq_len(nrow(nd))), , drop = FALSE]
  args <- list(variables = variables, vcov = FALSE)
  if (!is.null(by)) args$by <- c("group", "term", "contrast", by)
  joint <- do.call(marginaleffects::avg_comparisons, c(list(model = fit, newdata = nd), args))
  for (i in seq_along(fit$models)) {
    separate <- do.call(marginaleffects::avg_comparisons,
      c(list(model = fit$models[[i]], newdata = fit$model_frames[[i]]), args))
    key <- function(x, group) paste(x$term, group, x$contrast,
      if (is.null(by)) "" else x[[by]], sep = "|")
    group <- if (is.null(separate$group)) rep(fit$model_names[i], nrow(separate)) else
      paste0(fit$model_names[i], "::", separate$group)
    expected <- key(separate, group)
    idx <- match(expected, key(joint, joint$group))
    expect_false(anyNA(idx))
    expect_equal(sum(joint$group %in% unique(group)), nrow(separate))
    expect_equal(joint$estimate[idx], separate$estimate, tolerance = 1e-8)
  }
}

test_that("stacked multinomial samples retain every outcome-by-contrast estimate", {
  skip_if_not_installed("nnet")
  d <- alignment_data()
  for (samples in list(list(1:120, 121:240), list(1:160, 81:240), list(1:80, 81:160, 161:240))) {
    mods <- lapply(samples, function(ix) nnet::multinom(y ~ f + b + z + s, d[ix, ], trace = FALSE))
    fit <- do.call(suest, c(mods, list(model_names = paste0("Sample", seq_along(mods)))))
    compare_categorical_samples(fit, list(f = "pairwise", b = "reference", z = .5))
    compare_categorical_samples(fit, list(f = "pairwise"), by = "s")
  }
})

test_that("stacked scalar and categorical models retain their own prediction rows", {
  skip_if_not_installed("nnet")
  d <- alignment_data()
  m1 <- glm(I(y == "a") ~ f + z, d[1:120, ], family = binomial())
  m2 <- nnet::multinom(y ~ f + z, d[121:240, ], trace = FALSE)
  fit <- suest(m1, m2, model_names = c("Binary", "Nominal"))
  compare_categorical_samples(fit, list(f = "pairwise", z = .5))
})

test_that("stacked ordered and multinomial models retain different outcome sets", {
  skip_if_not_installed("nnet")
  d <- alignment_data()
  m1 <- MASS::polr(yo ~ f + z, d[1:120, ], Hess = TRUE)
  reduced <- d[121:240, ]; reduced$y <- droplevels(factor(ifelse(reduced$y == "d", "c", as.character(reduced$y))))
  m2 <- nnet::multinom(y ~ f + z, reduced, trace = FALSE)
  fit <- suest(m1, m2, model_names = c("Ordered", "Nominal"))
  compare_categorical_samples(fit, list(f = "pairwise", z = .5))
})

test_that("stacked categorical models retain their own predictor level sets", {
  skip_if_not_installed("nnet")
  d <- alignment_data()
  a <- d[1:120, ]
  b <- droplevels(d[d$id > 120 & d$f != "D", ])
  m1 <- nnet::multinom(y ~ f + z, a, trace = FALSE)
  m2 <- nnet::multinom(y ~ f + z, b, trace = FALSE)
  fit <- suest(m1, m2, model_names = c("Four levels", "Three levels"))
  compare_categorical_samples(fit, list(f = "pairwise", z = .5))
})
