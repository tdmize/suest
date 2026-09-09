# Seemingly Unrelated Estimation for Model Comparisons

The `suest` package combines two or more separately fitted regression
models into a single object with a joint model-robust covariance matrix.
The combined object works with the [marginaleffects
package](https://marginaleffects.com/) to compare predictions, slopes,
and average comparisons across models.

## Supported models

[`lm()`](https://rdrr.io/r/stats/lm.html), binary logit, probit, and
complementary-log-log models fitted by
[`glm()`](https://rdrr.io/r/stats/glm.html) or
[`glm2::glm2()`](https://rdrr.io/pkg/glm2/man/glm2.html), Poisson and
other supported GLMs, negative-binomial models fitted by
[`MASS::glm.nb()`](https://rdrr.io/pkg/MASS/man/glm.nb.html), ordered
logit and probit models fitted by
[`MASS::polr()`](https://rdrr.io/pkg/MASS/man/polr.html) or restricted
[`ordinal::clm()`](https://rdrr.io/pkg/ordinal/man/clm.html)
specifications, multinomial logit models fitted by
[`nnet::multinom()`](https://rdrr.io/pkg/nnet/man/multinom.html),
parametric survival and censored regression fitted by
[`survival::survreg()`](https://rdrr.io/pkg/survival/man/survreg.html),
beta regression fitted by
[`betareg::betareg()`](https://rdrr.io/pkg/betareg/man/betareg.html),
zero-inflated Poisson and negative-binomial models fitted by
[`pscl::zeroinfl()`](https://rdrr.io/pkg/pscl/man/zeroinfl.html),
truncated Gaussian regression fitted by
[`truncreg::truncreg()`](https://rdrr.io/pkg/truncreg/man/truncreg.html),
and direct Tobit models fitted by
[`censReg::censReg()`](https://rdrr.io/pkg/censReg/man/censReg.html),
and unweighted two-stage least squares fitted by
[`fixest::feols()`](https://lrberge.github.io/fixest/reference/feols.html)
without absorbed fixed effects. Heteroskedastic binary probit and logit
fitted by
[`Rchoice::hetprob()`](https://rdrr.io/pkg/Rchoice/man/hetprob.html) are
supported. Maximum-likelihood IV probit fitted by
[`Rchoice::ivpml()`](https://rdrr.io/pkg/Rchoice/man/ivpml.html) and
bivariate probit fitted by
[`mvProbit::mvProbit()`](https://rdrr.io/pkg/mvProbit/man/mvProbit.html)
are supported with their documented adapter restrictions. Unweighted
individual fixed-effects, between-effects, and Swamy-Arora
random-effects linear panel models fitted by
[`plm::plm()`](https://rdrr.io/pkg/plm/man/plm.html) are supported;
unbalanced random-effects fits can retain small engine-specific
differences. Single-level random-intercept Gaussian models fitted by
[`nlme::lme()`](https://rdrr.io/pkg/nlme/man/lme.html) with
`method = "ML"` are supported, including both variance components.
Unweighted GEE fitted by
[`geepack::geeglm()`](https://rdrr.io/pkg/geepack/man/geeglm.html) is
supported for Gaussian identity, binary logit/probit/cloglog, and
Poisson log models with independence or exchangeable correlation.
Explicit `weight_type = "pweight"` support is available for the model
families corresponding to Stata `suest2`'s ordinary pweight route. Joint
cluster-robust covariance is available through the `cluster` argument to
[`suest()`](https://tdmize.github.io/suest/reference/suest.md).
Compatible systems fitted in multiple imputed datasets can be pooled
with
[`suest_mi()`](https://tdmize.github.io/suest/reference/suest_mi.md)
using Rubin's rules. The initial MI interface supports joint
coefficients, covariance, summaries, and coefficient-level hypotheses.

## Reference

Mize, Trenton D., Long Doan, and J. Scott Long. 2019. "A General
Framework for Comparing Predictions and Marginal Effects Across Models."
*Sociological Methodology* 49(1):152–189.
[doi:10.1177/0081175019852763](https://doi.org/10.1177/0081175019852763)
