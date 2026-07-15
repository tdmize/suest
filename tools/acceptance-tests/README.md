# Full numerical acceptance tests

From the repository root, run:

```r
source("tools/acceptance-tests/run_all_tests.R")
```

The suite downloads the two public replication datasets on first use and
checks all six Mize, Doan, and Long (2019) examples, supported model families,
sample-overlap cases, and validation errors. The GitHub Actions
`paper-replications` workflow runs this suite automatically.


## Cross-model comparison universe

`tests/test_comparison_universe.R` systematically crosses all supported model
families with nested/mediator, alternative-predictor, different-outcome,
sex-stratified, missing-data, and partially overlapping-sample comparisons.
It also tests all allowed cross-family pairs on identical and disjoint samples.
See `COMPARISON-UNIVERSE.md` for the full matrix.
