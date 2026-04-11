# RdistNN news

## RdistNN (development)

### Diagnostics and reporting

- Added `diagnostics()` methods for `np_dnn()`, `cate_dnn()`, and `plqr_dnn()`
  fitted objects, including generic DNN weight/tuning summaries plus
  estimator-specific NP outcome, CATE overlap, and PLQR Gram conditioning
  diagnostics.
- Added a diagnostics/reporting vignette and updated README, getting-started,
  conditional Z-estimation, kernel-order tuning, and pkgdown references to show
  the diagnostics workflow.

### CI and package engineering

- Added `tools/ci/check_tests_declared_deps.R` and a CI guard step to fail when
  `tests/` uses `pkg::` or `pkg:::` for undeclared packages.
- Added `tools/ci/check_pkgdown_reference_coverage.R` and a CI guard step to
  fail when public `man/*.Rd` topics are missing from `_pkgdown.yml` reference
  contents.
- Added `tools/ci/check_generated_file_drift.R` and a CI guard step to fail
  when generated files (`NAMESPACE`, `R/RcppExports.R`,
  `src/RcppExports.cpp`) are stale.
- Added workflow `.github/workflows/alignment-guardrails.yaml` to run these
  guardrail checks on push and pull requests.

### Dependency and docs alignment

- Added `tibble` to `Suggests` for test-time namespace declarations.
- Filled missing `_pkgdown.yml` reference entries for TDNN and Mahalanobis
  topics that were failing site checks.
