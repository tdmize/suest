# Seemingly Unrelated Estimation for Model Comparisons

The `suest` package combines two separately fitted cross-sectional
regression models into a single object with a joint model-robust
covariance matrix. The combined object works with the [marginaleffects
package](https://marginaleffects.com/) to compare predictions, slopes,
and average comparisons across models.

## Supported models

[`lm()`](https://rdrr.io/r/stats/lm.html), binary logit and probit
models fitted by [`glm()`](https://rdrr.io/r/stats/glm.html) or
[`glm2::glm2()`](https://rdrr.io/pkg/glm2/man/glm2.html), Poisson models
fitted by either GLM engine, negative-binomial models fitted by
[`MASS::glm.nb()`](https://rdrr.io/pkg/MASS/man/glm.nb.html), ordered
logit and probit models fitted by
[`MASS::polr()`](https://rdrr.io/pkg/MASS/man/polr.html) or restricted
[`ordinal::clm()`](https://rdrr.io/pkg/ordinal/man/clm.html)
specifications, and multinomial logit models fitted by
[`nnet::multinom()`](https://rdrr.io/pkg/nnet/man/multinom.html).

## Reference

Mize, Trenton D., Long Doan, and J. Scott Long. 2019. "A General
Framework for Comparing Predictions and Marginal Effects Across Models."
*Sociological Methodology* 49(1):152–189.
[doi:10.1177/0081175019852763](https://doi.org/10.1177/0081175019852763)
