# RdistNN news

## RdistNN (development)

### Plotting layer (Phase 3)

- `autoplot()`/`plot()` methods for `plqr_result` (zero-referenced moment
  facets with a beta/tau subtitle) and `gmm_result` (coefficient facets
  with the model name); zero reference lines are now drawn for PLQR
  moments as well as CATE effects.
- `plot_conditioning()` for PLQR fits: smallest Gram eigenvalue and log10
  condition number across targets with caution/ill-conditioned guides at
  `cond_limits`.
- `plot_calibration()` for PLQR fits: the weighted below-threshold share
  against the nominal `tau` with a normal reference band scaled by the
  local effective sample size; the indicator is stored on new fits and
  can be supplied as `residuals` otherwise.
- `plot_influence()` for CATE fits: per-observation weighted score
  contributions against rank or distance, with the largest contributions
  labelled by observation index.
- `plot_tdnn_scales()` overlays a TDNN fit with internal single-scale
  refits at `s1` and `s2`.
- `plot_tuning()` for `np_result` and `cate_result`: refits across the
  tuning grid's candidate kernel orders (always including the fitted `s`)
  and plots estimate stability with the fitted order marked.
- Fitted objects store the additional components these panels need:
  training data on `np_dnn()`/`np_tdnn()` fits, and the below-threshold
  indicator on `plqr_dnn()` fits.
- Optional visual regression layer: vdiffr snapshots for the main panels,
  opt-in via the `RDISTNN_VDIFFR=true` environment variable so platform
  font differences cannot break CI.

### Plotting layer (Phase 2)

- `plot_overlap()` for CATE fits: `type = "mass"` draws the weighted
  treated/control mass per target with dotted share guides and shaded
  warning bands outside `share_limits` (violations reported in the
  caption); `type = "propensity"` draws weighted propensity densities at
  selected targets when observation-level propensities are available and
  a weighted q10-q90 summary otherwise. Inputs resolve through the
  documented precedence chain (stored fit components, then
  `treatment`/`propensity` arguments, then an actionable error).
- `autoplot()`/`plot()` on `diagnostics()` frames: effective sample size
  (with the degenerate-weights guide at `0.9 n`), maximum weight, and
  top-5 weight share faceted across targets.
- `plot_dashboard()` for `cate_result` and `np_result`: patchwork
  composition with a metadata annotation (`n`, `m`, `k`, `s` with
  provenance, distance, jackknife type); unavailable panels are
  substituted (diagnostics facets replace a missing overlap panel) or
  dropped, and without patchwork the panels are returned as a named list
  with a message. patchwork joins Suggests.
- Fitted `np_dnn()`, `cate_dnn()`, `plqr_dnn()`, and `np_tdnn()` objects
  record the jackknife type in `meta$jackknife` (`"none"` when `se =
  FALSE`).
- New vignette `visualising-fits` walks from a CATE fit through the
  estimate, weight, overlap, and diagnostics panels to the dashboard and
  panel customisation.
- The agent instruction files are excluded from the build, removing the
  long-standing top-level-files NOTE in R CMD check.

### Plotting layer (Phase 1)

- Added ggplot2-based plotting (ggplot2 is a suggested dependency; all
  entry points fail with an actionable message when it is missing):
  `plot_estimates()` draws point estimates with normal confidence bounds
  as a profile band (targets varying along one stored coordinate), a
  forest layout (small non-profile target sets, labelled from the stored
  coordinates), or pointranges, with an honesty caption carrying the
  median effective sample size and the kernel-order provenance on every
  panel; CATE panels add a dashed zero reference line and an opt-in
  `shade_significant` rug that marks intervals, never point estimates.
- `plot_weights()` draws weight decay by rank (support-edge guide at
  `n - s + 1`, optional log scale), cumulative weight shares (0.5/0.9
  guides and top-share markers), and the weight-vs-distance "effective
  kernel" view (rank-`s` radius guide); two-scale TDNN panels add a zero
  line for signed weights, and large target sets are subset to
  first/median-`n_eff`/last by default.
- `autoplot()` (registered on ggplot2) and `plot()` methods for
  `np_result` and `cate_result` delegate to `plot_estimates()`.
- README gains a visualisation section with the CATE profile figure.

### Plotting data layer (Phase 0)

- Added the tidy data extractors behind the upcoming ggplot2-based plotting
  interface, with no plotting dependency: `plot_data_estimates()` (S3
  generic with methods for `np_result`, `cate_result`, `plqr_result`, and
  `gmm_result`; TDNN fits dispatch through `np_result`) returns one row
  per target and coefficient with estimates, stored jackknife standard
  errors, and normal confidence bounds at a chosen `level`;
  `plot_data_weights()` returns per-(target, rank) weight profiles with
  cumulative shares and distances when stored; `plot_data_diagnostics()`
  returns the per-target diagnostics frame with an integer `target_id`
  aligned with the other extractors.
- Design and phase plan: `markdown_notes/plotting_layer_implementation_plan.md`.

### Fitted objects store evaluation points

- `np_dnn()`, `cate_dnn()`, `plqr_dnn()`, `gmm_dnn()`, and `np_tdnn()`
  result objects now include `x0`, the `m x p` matrix of evaluation
  points, with covariate column names propagated from the training data
  when the targets carry none. This is Phase 0 of the plotting layer
  (`markdown_notes/plotting_layer_implementation_plan.md`): profile
  plots need the target coordinates, which were previously not stored.

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
- `np_tdnn()` still uses explicit `s1`/`s2` defaults; a two-scale default
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
