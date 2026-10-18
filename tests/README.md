# Test-suite map

Run tests from the repository root with `devtools::test()`. The R package
suite lives directly in `tests/testthat/`; `tests/ci/` contains separate
Python checks for repository and workflow policy.

## Finding a test

Filenames describe the behavior under test. Use the matching `R/` module as
the first place to look, then the relevant integration or numerical suite.

| Area | Files in `tests/testthat/` |
| --- | --- |
| Component contracts | `test-api-*`, `test-distance-*`, `test-order-*`, `test-score-*`, `test-utils-validate.R` |
| Numerical references and optimized paths | `test-rank-weight-*`, `test-jackknife-rank-*`, `test-dnn-scalar-reference.R`, `test-dnn-jackknife-prefix.R`, `test-tdnn-*`, `test-compiled-fallback-parity.R`, `test-errors-fastpaths.R` |
| Support-aware jackknife centering across estimators, compiled/R and fast/naive paths | `test-jackknife-support-centering.R` |
| CATE and nonparametric regression | `test-user-cate-*`, `test-user-np-dnn*`, `test-user-tdnn-np*`, `test-demo-data-npr.R` |
| Experimental conditional LATE | `test-clate-*`, `test-user-clate-*`, `test-user-affine-ratio-*`, `test-weights-single-row.R` |
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

## Helpers and optional tests

`test-jackknife-support-centering.R` covers opposite-sign zero-weight outliers,
extreme unused values, real local variation, target-specific deletion support,
the boundary rank that gains weight after deletion, signed TDNN weights and
asymptotic tails with both normalization settings. Its analytical variance-9
fixture is also enforced in the Stata companion. Both packages retain their
explicit deletion algorithms for parity comparisons and fallback use.

Shared fixtures and reference calculations live in domain-specific
`helper-*.R` files. Fixtures used only by one suite stay beside its tests.
`helper-dnn-reference.R` preserves an explicit scalar deletion reference;
`helper-sim-cate.R` holds independent AIPW calculations;
`helper-user-dnn-gmm.R` supplies custom GMM model fixtures. Helper files define
support functions rather than fitting models eagerly.

`helper-clate-reference.R` enumerates complete subsamples independently of the
package's weight/deletion engines. CLATE suites cover paired covariance,
signed and zero denominators, deterministic ties, R/compiled/explicit-loop
parity, methods and serialization. The bounded ridge suite checks independent
dense normal equations and all required one-/two-deleted fits. Confidence-set
tests cover every quadratic geometry, membership, singular covariance and
extreme scales. Small generated examples test execution, not coverage.

The CLATE alignment suites additionally check literal asymptotic weights and
paired deletions, exhaustive exact general-k subset averages, perfect-compliance
reduction to CATE, common summary/tidy schemas and every candidate-dependent
covariance/extras field. Plot tests cover instrument versus treatment overlap,
paired score contributions, unavailable bounds, serialization, optional packages,
and real ggplot builds. The shared rank-weight and plot-builder changes also
require the CATE/NP suites. The implementation/extension map is in
[result architecture](../markdown_notes/design/13_estimator_and_result_architecture.md)
and [paired-ratio contract](../markdown_notes/design/17_clate_ratio_inference.md).

Test files should run independently. Testthat's helper loading and file
conventions are documented in its
[special-files guide](https://testthat.r-lib.org/articles/special-files.html).
The package's seed, numerical, validation, and fallback requirements remain
in [the testing policy](../markdown_notes/design/09_testing_and_regression_policy.md),
routed from [AGENTS.md](../AGENTS.md).

Optional plotting packages and compiled capabilities have explicit guards.
Visual snapshots additionally require a deliberate local opt-in:

```bash
RDISTNN_VDIFFR=true Rscript -e 'devtools::test(filter = "^plot-vdiffr$")'
```

Review any changed SVGs in `tests/testthat/_snaps/plot-vdiffr/` before accepting
them. A normal run retains numerical plot-data checks without enabling
platform-dependent visual snapshots.
