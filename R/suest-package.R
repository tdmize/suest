#' Seemingly Unrelated Estimation for Model Comparisons
#'
#' The `suest` package combines two or more separately fitted regression
#' models into a single object with a joint model-robust covariance
#' matrix. The combined object works with the
#' [marginaleffects package](https://marginaleffects.com/) to compare
#' predictions, slopes, and average comparisons across models.
#'
#' @section Supported models:
#' `lm()`, binary logit, probit, and complementary-log-log models fitted by
#' `glm()` or `glm2::glm2()`, Poisson and other supported GLMs,
#' negative-binomial models fitted by [MASS::glm.nb()], ordered logit and
#' probit models fitted by [MASS::polr()] or restricted `ordinal::clm()`
#' specifications, multinomial logit models fitted by [nnet::multinom()], and
#' parametric survival or censored-regression models fitted by
#' `survival::survreg()`, and beta regressions fitted by `betareg::betareg()`.
#' Poisson and negative-binomial zero-inflated models fitted by
#' `pscl::zeroinfl()` are also supported.
#' Truncated Gaussian regressions fitted by `truncreg::truncreg()` are
#' supported as well, as are direct Tobit models fitted by
#' `censReg::censReg()` and unweighted two-stage least squares fitted by
#' `fixest::feols()` without absorbed fixed effects. Heteroskedastic binary
#' probit and logit fitted by `Rchoice::hetprob()` are supported.
#' Maximum-likelihood IV probit fitted by `Rchoice::ivpml()` and bivariate
#' probit fitted by `mvProbit::mvProbit()` are supported with their documented
#' adapter restrictions.
#' Unweighted individual fixed-effects, between-effects, and Swamy-Arora
#' random-effects linear panel models fitted by `plm::plm()` are supported;
#' unbalanced random-effects fits can retain small engine-specific differences.
#' Single-level random-intercept Gaussian models fitted by `nlme::lme()` with
#' `method = "ML"` are supported, including both variance components.
#' Unweighted GEE fitted by `geepack::geeglm()` is supported for Gaussian
#' identity, binary logit/probit/cloglog, and Poisson log models with
#' independence or exchangeable correlation.
#' Explicit `weight_type = "pweight"` support is available for the model
#' families corresponding to Stata `suest2`'s ordinary pweight route.
#' Joint cluster-robust covariance is available through the `cluster` argument
#' to [suest()].
#' Compatible systems fitted in multiple imputed datasets can be pooled with
#' [suest_mi()] using Rubin's rules. The initial MI interface supports joint
#' coefficients, covariance, summaries, and coefficient-level hypotheses.
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
