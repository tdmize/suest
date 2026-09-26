# Post-0.1.6 stability and usability check

This development pass adds clearer unsupported-model diagnostics and a practical
check for numerical sensitivity. It does not alter model estimates, the native
`glmmTMB` covariance, or the covariance conventions of `suest`.

## Scale predictors before fitting random slopes

When a continuous predictor has very large or very small units, or is far from
zero, center and scale it **in the data before fitting both component models**.
For example, when the random slope is `x` and the data contain `y1`, `y2`, `z`,
and panel `id`:

```r
sx <- sd(d$x)
mx <- mean(d$x)
d$x_scaled <- (d$x - mx) / sx

m1 <- glmmTMB::glmmTMB(y1 ~ x_scaled + z + (1 + x_scaled | id),
  data = d, family = binomial("logit"))
m2 <- glmmTMB::glmmTMB(y2 ~ x_scaled + z + (1 + x_scaled | id),
  data = d, family = binomial("logit"))
fit <- suest::suest(m1, m2, observation_id = "id")
```

Use the actual observation ID columns required to distinguish rows in your
data; for repeated observations within a panel, supply both panel and occasion
IDs, for example `observation_id = c("id", "time")`. Refit the original-unit
models for a sensitivity comparison when they converge. Predictions evaluated
on corresponding rows should agree after the change of units. Compare effect
estimates and standard errors *on the same units*: a slope per original unit is
the slope for `x_scaled` divided by `sx`. This is a sensitivity check, with no
universal numerical cutoff. If the native covariance changes materially under
an equivalent rescaling, investigate the fits and prefer a well-scaled fit.
Reported convergence and a positive-definite Hessian alone do not certify
numerical curvature.

## Check numerical logit slope standard errors

For response-scale effects from Bernoulli logit random-slope models, use centered
coefficient differences, then compare a few nearby step sizes:

```r
steps <- c(5e-5, 1e-4, 2e-4)
slopes <- lapply(steps, function(step)
  marginaleffects::avg_slopes(fit, variables = "x_scaled", newdata = d,
    numderiv = list("fdcenter", eps = step)))

lapply(slopes, function(tab)
  data.frame(group = tab$group, estimate = tab$estimate / sx,
    std.error = tab$std.error / sx))
```

The example converts the reported average slopes and standard errors to effects
per original unit of `x`. A material change across steps calls for numerical
investigation. A stable finite-difference step does not establish that the
native `glmmTMB` covariance is accurate; the rescaling check addresses a
different source of error. The package does not change the global
`marginaleffects` numerical-derivative settings.

## Unsupported specifications

For a model family outside the supported routes, `suest()` now reports the
family and link for `glmmTMB`; the generic unsupported-model error identifies
the model name and class. Existing specific errors still cover unsupported
random-effects structures, weights, offsets, constraints, dispersion models,
and convergence failures. See `?suest` for the supported model contract.

## Validation boundary

The release's independent R likelihood, score, covariance, and marginal-effect
checks remain distinct from raw Stata GSEM comparisons. The unresolved raw
GSEM covariance differences and the disjoint logit reconstruction error
(1.052e-6 against the original 1e-6 bound) remain documented in
`RELEASE-0.1.6.md` and the linked validation reports. This usability pass does
not reclassify that diagnostic as passing.
