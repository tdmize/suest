# Cross-model comparison test universe

The expanded acceptance suite crosses each supported model family with six
common comparison designs.

| Comparison design | lm | logit | probit | Poisson | negative binomial | ordered logit | ordered probit | multinomial |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| Nested model / mediator adjustment | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |
| Alternative predictor operationalizations | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |
| Different outcomes | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |
| Men versus women using `subset=` | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |
| Different samples from mediator missingness | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |
| Partially overlapping row subsets | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |

For every cell, the suite verifies:

1. The two models combine successfully.
2. The reported component, overlap, and union sample sizes are correct.
3. The cross-model covariance is nonzero for overlapping samples.
4. The cross-model covariance is exactly zero for disjoint samples.
5. `marginaleffects::avg_comparisons()` produces finite estimates and standard
   errors.
6. A direct cross-model hypothesis produces finite estimates and standard
   errors.
7. Category-specific effects are matched correctly for ordered and
   multinomial outcomes.

The suite also tests every allowed cross-family pair on both the same sample
and disjoint men/women samples:

- logit versus probit;
- logit versus linear probability model;
- probit versus linear probability model;
- Poisson versus negative binomial;
- ordered logit versus multinomial logit.

Finally, a separate test confirms the documented behavior that models fitted
from separately filtered data objects are treated as disjoint when common
observations cannot be inferred safely.


## Alternative engines

A separate adapter suite tests that the comparison universe is not tied to a
single fitting function:

- `glm2::glm2()` is tested for binary logit, binary probit, and Poisson models;
- `ordinal::clm()` is tested for ordered logit and ordered probit models;
- `stats::glm()` and `glm2::glm2()` are compared directly;
- `ordinal::clm()` and `MASS::polr()` are compared directly;
- `clm()` is tested with identical, disjoint, and partially overlapping
  samples;
- unsupported `clm()` extensions and adjusted-score GLMs are rejected.
