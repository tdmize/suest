# Probability-weight implementation

Version 0.1.4 introduces explicit `weight_type = "pweight"` support for
`stats::lm()`, binary logit/probit and Poisson models from `stats::glm()` or
`glm2::glm2()`, `MASS::glm.nb()`, ordered logit/probit from `MASS::polr()`,
and `nnet::multinom()`.

For model \(m\), let \(X_m\) be the model matrix, \(e_m\) the residual vector,
and \(w_m\) the evaluated sampling weights. The estimating-function
contribution for observation \(i\) is

\[
u_{mi} = w_{mi} x_{mi} e_{mi},
\]

and the coefficient bread is based on

\[
A_m = X_m^\mathsf{T} W_m X_m.
\]

The existing SUEST alignment machinery places score contributions on the union
of the two estimation samples. Shared observations generate the cross-model
meat; model-specific observations contribute zero to the other model. For any
within- or cross-model block, Stata `suest` applies the same union-sample
finite-sample correction:

\[
\frac{n_{\mathrm{union}}}{n_{\mathrm{union}}-1}
A_j^{-1}
\left(\sum_i u_{ji}u_{ki}^\mathsf{T}\right)
A_k^{-1}.
\]

For a within-model block, \(j=k\). For a cross-model block, the score product
is nonzero only for observations included in both members of a model pair.

## Linear-model ancillary parameter

Like Stata `suest2`, pweighted linear systems include a Gaussian ancillary
parameter \(\lambda=\log(\sigma^2)\). The iweight-reference construction used
by `suest2` implies

\[
\sigma^2 = \frac{\sum_i w_i e_i^2}{\sum_i w_i-k},
\]

where \(k\) is the number of mean parameters. Its score contributions are

\[
u_{\beta i}=\frac{w_i x_i e_i}{\sigma^2}, \qquad
u_{\lambda i}=\frac{w_i}{2}\left(\frac{e_i^2}{\sigma^2}-1\right).
\]

In the package's `sandwich::bread()` convention, the block-diagonal bread is

\[
B_\beta=n\sigma^2(X^\mathsf{T}WX)^{-1}, \qquad
B_\lambda=\frac{2n}{N_w},
\]

where \(N_w=\lfloor\sum_iw_i+0.5\rfloor\) reproduces Stata's effective
integer iweight count. A deterministic Stata 19.5 clustered benchmark agrees
for the full coefficient and covariance matrices, including all mean--ancillary
and cross-model blocks.


## Binary logit and probit models

For binary models, let \(p_i\) be the fitted probability, \(\eta_i\) the
linear predictor, and \(f_i\) the standard-normal density for a probit model.
The weighted score contribution is

\[
u_i = w_i x_i s_i,
\]

where \(s_i = y_i-p_i\) for logit and

\[
s_i = \frac{(y_i-p_i)f_i}{p_i(1-p_i)}
\]

for probit. The logit information contribution is
\(p_i(1-p_i)\). For probit, the implementation uses the observed-information
contribution, not the Fisher-scoring working weight used by `glm()`. This
choice is necessary to reproduce Stata's pweighted `suest` covariance matrix.
The same union-sample correction and overlap alignment used for linear models
then apply to logit-logit, probit-probit, logit-probit, linear-logit, and
linear-probit pairs.

## Poisson log-link models

For Poisson models with conditional mean \(\mu_i=\exp(x_i^\mathsf{T}\beta)\),
the weighted score contribution is

\[
u_i = w_i x_i(y_i-\mu_i),
\]

and the information matrix is based on

\[
A = X^\mathsf{T} W\,\mathrm{diag}(\mu)X.
\]

Because the log link is canonical for the Poisson model, the observed and
expected coefficient information coincide. The same union-sample correction
and score alignment apply to identical, partially overlapping, and disjoint
Poisson model pairs.

## Weight compatibility

- Weights are compared by evaluated value, not variable name.
- Weights must agree, within numerical tolerance, on shared observations.
- Weights may differ on observations unique to one model.
- Disjoint samples may use entirely different weight vectors.
- Multiplying all weights in every model by the same constant does not change
  mean-parameter coefficients, their covariance blocks, or
  `marginaleffects` results.
- The linear-model `lnvar` nuisance parameter follows Stata's iweight-reference
  normalization and is therefore not generally invariant to arbitrary weight
  rescaling. A dedicated Stata diagnostic records this behavior separately
  from the scale-invariant substantive results.
- Zero, negative, missing, and infinite pweights are rejected.

Stata imposes a stricter implementation restriction and requires the same
weighting expression across component models. The R implementation preserves
Stata results in supported cases but permits different weights where no
observation receives conflicting weights.

## Finite-sample correction

The union-sample correction differs from applying a separate
`n_model / (n_model - 1)` factor to each component model when their estimation
samples differ. The Stata partial-overlap and disjoint-sample benchmarks test
this distinction directly.

With `cluster`, ordinary weighted scores are instead aggregated at the system
cluster level and receive Stata's \(G/(G-1)\) adjustment, where \(G\) is the
number of union-sample clusters.
