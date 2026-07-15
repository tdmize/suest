# Model-adapter architecture

`suest()` separates the statistical model family from the fitting engine.
This prevents package-specific object layouts from leaking into the joint
covariance and `marginaleffects` interfaces.

## Internal adapter functions

- `.suest_model_adapter(model)` identifies the statistical type and engine.
- `.suest_model_frame(model, engine)` recovers the estimation rows.
- `.suest_category_levels(model, engine)` returns ordered or nominal outcome
  categories.
- `.suest_extract_parameters(model, type, engine)` returns parameters in score
  order.
- `.suest_model_components(model, type, engine)` returns observation-level
  scores and the bread matrix.
- `.suest_set_parameters(model, parameters, type, engine)` updates a fitted
  object for numerical differentiation by `marginaleffects`.
- `.suest_predict_probabilities()` and `.suest_predict_values()` standardize
  predictions across engines.

## Engines supported in version 0.1.1

| Engine | Statistical types | Initial restrictions |
|---|---|---|
| `stats::lm()` | linear | unweighted, no offset |
| `stats::glm()` | logit, probit, Poisson | ordinary maximum likelihood |
| `glm2::glm2()` | logit, probit, Poisson | ordinary maximum likelihood |
| `MASS::glm.nb()` | negative binomial | log link |
| `MASS::polr()` | ordered logit, ordered probit | standard threshold model |
| `ordinal::clm()` | ordered logit, ordered probit | flexible thresholds, no scale or nominal formula |
| `nnet::multinom()` | multinomial logit | `summ = 0` |

## Explicitly rejected GLM-like estimators

Objects fitted using bias reduction, adjusted score equations, Firth
corrections, or penalized likelihood are rejected. Inheriting from `"glm"` is
not enough: the observation-level estimating equations must match the ordinary
maximum-likelihood GLM scores used by the SUEST covariance.

## Candidate engines for later versions

The adapter layer is designed to support additional engines after separate
validation:

1. cross-sectional fixed-parameter `Rchoice::Rchoice()` models;
2. `fixest` models without absorbed fixed effects, followed by a separate
   investigation of absorbed effects;
3. ordinary maximum-likelihood `mclogit::mblogit()` models;
4. selected large-data engines when estimation rows and score contributions
   can be recovered safely.

Survey-weighted and random-intercept models require separate covariance
designs and are intentionally deferred to later development phases.
