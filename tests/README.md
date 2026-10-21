# Test-suite map

Run tests from the repository root with `devtools::test()`. The R package
suite lives directly in `tests/testthat/`; `tests/ci/` contains separate
Python checks for repository and workflow policy.

## Finding a test

Filenames describe the behavior under test. Use the matching `R/` module as
the first place to look, then the relevant integration or numerical suite.

| Area | Files in `tests/testthat/` |
| --- | --- |
| Saved geometry and shared inputs | `test-fitted-geometry-contract.R`, `test-input-route-contract.R`, `test-prediction-selectors.R` |
| Weight control and scalar/grid consumers | `test-weight-control-contract.R`, `test-tuning-consumer-contract.R` |
| Shared grid callbacks, failure rows, duplicate keys and cached target/inference reuse | `test-tuning-grid-engine.R` |
| Local origins, offsets, translated SD and covariance | `test-centering-translation.R` |
| Neighbor extras and CATE nuisance contracts | `test-result-neighbors.R`, `test-cate-nuisance-contract.R`, `test-score-cate-accumulation.R` |
| Exact-table represented parity and near-full support | `test-exact-table-portability.R` |
| Weighted CDF boundary indices, physical radius and outcome SD ranges | `test-weighted-distribution-range.R`, `fixtures/weighted-cdf.csv` |
| Component contracts | `test-api-*`, `test-distance-*`, `test-order-*`, `test-score-*`, `test-utils-validate.R` |
| Compensated weighted sums, native/R parity, retained deletion corrections and buffer boundaries | `test-compensated-accumulation.R` |
| Numerical references and optimized paths | `test-rank-weight-*`, `test-jackknife-rank-*`, `test-dnn-scalar-reference.R`, `test-dnn-jackknife-prefix.R`, `test-tdnn-*`, `test-compiled-fallback-parity.R`, `test-errors-fastpaths.R` |
| Native/R distance arithmetic, large-offset and accumulation-sensitive ranks, scalar/batch agreement and FMA-sensitive ties | `test-distance-scaled-arithmetic.R` |
| Huge/tiny finite L2 norms, shared exceptional boundaries, exact ranks and serialized prediction | `test-distance-range-stability.R` |
| Close-order TDNN coefficients and tiny Beta interval masses against independent references | `test-weight-numerical-stability.R` |
| Support-aware jackknife centering across estimators, compiled/R and fast/naive paths | `test-jackknife-support-centering.R` |
| Shared NP/CATE/paired-ratio prediction boundary | `test-prediction-targets-common.R` |
| CATE and nonparametric regression | `test-user-cate-*`, `test-user-np-dnn*`, `test-user-np-methods.R`, `test-user-tdnn-np*`, `test-demo-data-npr.R` |
| Conditional LATE and paired ratios | `test-clate-*`, `test-user-clate-*`, `test-user-affine-ratio-*`, `test-weights-single-row.R` |
| Tuning and stored candidate views | `test-api-tuning.R`, `test-user-dnn-tuning.R`, `test-tuning-*`, `test-candidate-fit.R`, `test-user-*-tuning-grid.R` |
| Diagnostics and plots | `test-user-dnn-diagnostics.R`, `test-plot-*`, `test-tuning-grid-diagnostics-plots.R` |
| CATE simulations | `test-sim-cate-data.R`, `test-sim-cate-one-rep.R`, `test-sim-cate-run.R`, `test-sim-cate-summary.R`, `test-sim-cate-tuning*.R` |
| Experimental PLQR and GMM/IV | `test-user-plqr-*`, `test-user-dnn-gmm*`, `test-user-gmm-methods.R`, `test-gmm-models.R`, `test-score-plqr.R`, `test-demo-data-plqr.R`, `test-demo-data-iv.R`, `test-sim-iv-*` |
| Integration, replication, and documentation | `test-integration-*`, `test-paper-*`, `test-result-*`, `test-*-documentation-*`, `test-notation-*`, `test-readme-*`, `test-citation-*`, `test-cate-paper-diagnostics.R`, `test-s-default-paper-rule.R` |

The larger wrapper families separate estimator identities, data interfaces,
validation, and result methods where those contracts warrant separate files.
The generic score derivative API is in `test-api-scores-diff.R`; the public
PLQR score matrix is in `test-user-plqr-score-matrix.R`.

## Focused runs

`filter` is a regular expression matched against filenames after removing
`test-` and `.R`. Anchor both ends to select one file:

```r
devtools::test(filter = "^sim-cate-summary$")
devtools::test(filter = "^user-dnn-gmm-validation$")
devtools::test(filter = "^jackknife-support-centering$")
```

Select a complete family when a change crosses contracts:

```r
devtools::test(filter = "^user-cate")
devtools::test(filter = "^user-(np-dnn|tdnn-np)")
devtools::test(filter = "^(api-tuning|tuning-|candidate-fit|user-dnn-tuning|user-.*-tuning-grid|sim-cate-tuning)")
devtools::test(filter = "^(plot-|user-dnn-diagnostics|tuning-grid-diagnostics-plots)")
devtools::test(filter = "^sim-cate")
devtools::test(filter = "^(clate-|user-clate|user-affine-ratio|weights-single-row)")
devtools::test(filter = "^(user-plqr|user-dnn-gmm|user-gmm|gmm-models|score-plqr|demo-data-(plqr|iv)|sim-iv)")
```

Use the full suite for changes shared across families. The release check
also runs these tests against the installed package; see
[RELEASING.md](../RELEASING.md).

## References, engines and optional checks

Domain-specific `helper-*.R` files define reusable fixtures without fitting
models eagerly. Single-suite fixtures stay beside their tests.
`helper-dnn-reference.R` supplies explicit scalar deletions,
`helper-sim-cate.R` independent AIPW calculations, `helper-user-dnn-gmm.R`
custom models, and `helper-clate-reference.R` complete-subset enumeration.
CLATE tests retain independent dense ridge solves, one-/two-deleted fits and
confidence-set geometry references. Small simulations check execution, not
coverage or manuscript claims.

Force and establish each engine; test automatic fallback separately. Use
`options(RdistNN.use_compiled = FALSE)` for the R path and require actual native
symbols for native acceptance. Exact ranks, support, statuses and missingness
accompany quantity-specific numeric tolerances. The
[testing policy](../markdown_notes/design/09_testing_and_regression_policy.md)
and [compiled-code policy](../markdown_notes/design/10_rcpp_and_compiled_code_policy.md)
own these acceptance rules. Shared primitive changes also require their
estimator, candidate, covariance and saved-prediction consumers.

Tests must run independently through public interfaces. Exercise S3 methods
through their generics on fitted objects; missing unexported aliases must not
skip installed-package coverage. Helper loading follows testthat's
[special-files guide](https://testthat.r-lib.org/articles/special-files.html).
Optional packages and compiled capabilities have explicit guards; native
requirements in native checks are failures rather than silent skips.
The Linux R-release job runs installed tests in an isolated library containing
Rcpp, testthat and their strong dependencies. Reproduce it with
`Rscript tools/ci/check_minimal_dependencies.R path/to/RdistNN_VERSION.tar.gz`.

Visual snapshots require a deliberate opt-in:

```sh
RDISTNN_VDIFFR=true Rscript -e 'devtools::test(filter = "^plot-vdiffr$")'
```

Inspect changed SVGs in `tests/testthat/_snaps/plot-vdiffr/`. Ordinary runs
retain numerical plot-data checks without platform-dependent visual snapshots.

## Package-maintenance checks

[CONTRIBUTING.md](../CONTRIBUTING.md#editing-and-checking) lists the commands.
`tests/ci/` uses isolated positive/negative fixtures for repository/workflow
policy, generated-file drift, candidate hygiene and function ownership.
Release checks run the installed-package suite; see
[RELEASING.md](../RELEASING.md).
