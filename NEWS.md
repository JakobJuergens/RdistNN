# RdistNN news

## RdistNN (development)

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
