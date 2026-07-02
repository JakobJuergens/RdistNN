# RdistNN news

## RdistNN (development)

### Breaking change: default kernel order `s`

- `np_dnn()`, `cate_dnn()`, `plqr_dnn()`, and `gmm_dnn()` no longer default
  to `s = 1L`, which produced exactly uniform DNN weights and therefore no
  localization (every target point silently received the global weighted
  mean). When both `s` and `tuning` are omitted, `s` is now resolved from
  the paper undersmoothing rule `tuning_dnn_paper()`,
  `s = floor(n^((kx + 1)/(kx + 2)))`, matching the default already used by
  the Stata companion (`cate_dnn.ado`). The resolved rule is recorded in
  `meta$tuning` with `source = "paper_default"`. Calls that passed `s`
  explicitly are unaffected.
- All DNN and TDNN user wrappers now emit a warning when the computed
  weights are (near-)uniform (max effective sample size at least `0.9 * n`,
  for example `s = 1` or `k = s`), instead of silently returning an
  unlocalized fit.
- Added `tests/testthat/test-s-default-paper-rule.R`, the `s`-analogue of
  the `k = 1` documentation guardrail: it pins the default-`s` resolution
  for all four wrappers, the manual-`s` path, the degenerate-weight
  warning, and R/Stata default-formula agreement.
- `tdnn_np()` still uses explicit `s1`/`s2` defaults; a two-scale default
  rule is a separate follow-up.

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
