#' Seemingly Unrelated Estimation for Model Comparisons
#'
#' The `suest` package combines two separately fitted cross-sectional
#' regression models into a single object with a joint model-robust covariance
#' matrix. The combined object works with the
#' [marginaleffects][marginaleffects::marginaleffects-package] package to compare
#' predictions, slopes, and average comparisons across models.
#'
#' @section Supported models:
#' `lm()`, binary logit and probit models fitted by `glm()`, Poisson models
#' fitted by `glm()`, negative-binomial models fitted by
#' [MASS::glm.nb()], ordered logit and probit models fitted by
#' [MASS::polr()], and multinomial logit models fitted by
#' [nnet::multinom()].
#'
#' @section Reference:
#' Mize, Trenton D., Long Doan, and J. Scott Long. 2019.
#' "A General Framework for Comparing Predictions and Marginal Effects Across
#' Models." *Sociological Methodology* 49(1):152--189.
#' \doi{10.1177/0081175019852763}
#'
#' @keywords internal
"_PACKAGE"
