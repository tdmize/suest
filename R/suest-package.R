#' Seemingly Unrelated Estimation for Model Comparisons
#'
#' The `suest` package lets you compare predictions and marginal effects across
#' regression models. It combines two or more separately fitted models into one
#' object with a joint robust covariance matrix. The combined object works with
#' the [marginaleffects package](https://marginaleffects.com/) to compare
#' predictions, slopes, and average comparisons across models.
#'
#' @section Supported models:
#' Linear, binary, ordinal, multinomial, count, censored, survival, and
#' instrumental-variable models, as well as panel, multilevel, and survey
#' models. See [suest()] for the full list.
#'
#' Probability weights (`weight_type = "pweight"`) are supported for linear,
#' binary logit and probit, Poisson, negative-binomial, ordered logit and
#' probit, and multinomial logit models. Joint cluster-robust covariance is
#' available through the `cluster` argument to [suest()]. Systems fitted in
#' multiply imputed datasets can be pooled with [suest_mi()] using Rubin's
#' rules.
#'
#' @section Reference:
#' Mize, Trenton D., Long Doan, and J. Scott Long. 2019.
#' "A General Framework for Comparing Predictions and Marginal Effects Across
#' Models." *Sociological Methodology* 49(1):152--189.
#' \doi{10.1177/0081175019852763}
#'
#' @aliases suest-package NULL
#' @importFrom marginaleffects predictions
#' @importFrom stats coef nobs vcov
#' @keywords internal
"_PACKAGE"
