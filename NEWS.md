# suest 0.1.0

* Initial GitHub release.
* Combines exactly two supported cross-sectional models using a joint
  model-robust covariance matrix.
* Integrates with `marginaleffects` for predictions, comparisons, slopes,
  hypotheses, and plots.
* Supports identical, partially overlapping, and disjoint estimation samples.
* Includes replication vignettes for Examples 6.1--6.6 in Mize, Doan, and
  Long (2019).


* Expanded the numerical acceptance suite to cover a full cross-model
  comparison-by-model-family matrix, including nested/mediator comparisons,
  alternative predictor operationalizations, different outcomes,
  sex-stratified samples, missing-data sample changes, and partially
  overlapping samples.
* Added same-sample and disjoint-sample tests for every supported cross-family
  comparison.

* Imported the `stats` generics used by the `suest_model` S3 methods so the
  namespace loads correctly during `R CMD check` and pkgdown builds.
* Updated GitHub Actions checkout steps to `actions/checkout@v6` and the
  GitHub Pages deployment action to version 4.8.0.
