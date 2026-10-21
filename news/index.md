# Changelog

## RdistNN 0.1.1 (unreleased)

- Resolve exact distance ties lexicographically by the original
  covariates across R and native implementations. Identical covariate
  rows retain input order. Rank weights and effective sample size are
  preserved.
- Preserve PLQR fitted distance specifications and resolved covariate
  names during saved prediction; align scalar tuning-plot uncertainty
  with the requested policy. Permit full-order analytic PLQR inference
  when no deletion calculation is needed.
- Correct compensated PLQR score reduction and invalid
  coefficient/subset validation, signed rank-table mass validation, and
  range-sensitive weighted SD and Mahalanobis origin selection. Cover
  one-covariate bump demos.
- Share adaptive PLQR comparisons across score constructors, native
  kernels and retained diagnostics. Exact accumulation resolves ordinary
  cancellation, intermediate overflow, tiny residuals and inclusive
  threshold equality.
- Share finite-operand weighted-contribution recovery between CATE
  diagnostics and CATE/CLATE influence plots; preserve representable
  dominant observations and concentration when intermediate residuals or
  totals exceed the range.
- Detect undeclared ordinary namespace references in the test dependency
  check.
- Refresh pinned paper-replication displays with embedded PDF fonts and
  isolate nested renderer dependencies from caller bindings, including
  on errors.
- Remove obsolete development-fit compatibility readers and their
  dedicated fixtures. Refit objects created by older development
  versions.
- Name simulation outputs `localization_drift` and
  `standardized_localization_drift` and validate the current
  simulation-row schema rather than filling missing fields as successful
  results.
- Separate estimator fitting, grid adapters, retained-state support and
  public methods; share identical grid warning/metadata policy.
  Reorganize test ownership, locally scope fixture RNG state and
  complete the documented focused-suite selections.
- Record CATE signal provenance (`aipw` or `supplied`) and the original
  fit’s loaded package version. Summaries and candidates preserve that
  original provenance and version; current retained state is required
  for reuse.
- Correct demo subset helpers to return integer `K` by `s` matrices even
  for singleton subsets, preserving draws and caller RNG behavior.
- Clarify current CATE inputs and simulation reproducibility: reproduce
  simulations under the same RNG kind and backend configuration.
- Credit Zakareya Chehime for quality assurance, testing and validation
  in contributor metadata and documentation.

## RdistNN 0.1.0.9005 (development snapshot)

The maintenance and documentation changes below ship together as one
flattened (squash-merged) production commit. Version 0.1.0.9004 records
the previous estimator update; tags and assets of the existing v0.1.0
release are unchanged.

### Maintenance

- Limit native OpenMP teams to two threads by default, with
  package-specific opt-in controls that respect lower runtime limits.
- Restore caller saved RNG state after explicitly seeded demos and
  simulations, preserving seeded data and learner draws; `NULL` uses the
  current stream. Make stochastic CATE tuning callbacks reproducible
  under the replicate seed.
- Preserve the generated vignette index in source packages.
- Add [`nobs()`](https://rdrr.io/r/stats/nobs.html) for NP, TDNN and
  CATE results, reporting the training sample size.
- Validate demo targets before touching RNG state; preserve valid-call
  draws.
- Guard plotting-only tests when ggplot2 is absent and assert requested
  CATE outputs instead of skipping their checks. Add an isolated
  dependency check to Linux CI. Keep CI-only lintr in the lint workflow
  rather than package suggestions.

### Documentation

- Shorten READMEs and introductory guides around runnable user
  workflows. Point to the paper for theory and keep detailed syntax and
  engine contracts in their reference documents.
- Refresh merged-source installation pins, website copy and the R/Stata
  documentation crosswalk.
- Reconcile function help with current prediction, diagnostics, score
  evaluation, PLQR uncertainty and simulation interfaces; correct
  outdated examples and labels.
- Focus vignette workflows, execute lightweight examples during
  rendering and seed random reference examples. Keep theory and full
  replication protocols in their primary references.
- Review all 13 indexed vignettes; update sample-count, seed and thread
  guidance, simplify extension examples and remove repeated navigation
  prose.
- Simplify replication and contributor/test documentation, consolidate
  repeated design guidance and retain implementation contracts in their
  owning references.

## RdistNN 0.1.0.9004 (merged development snapshot)

The estimator update was published as one reviewed, flattened
(squash-merged) production commit in
[PR](https://github.com/JakobJuergens/RdistNN/pull/16)
[\#16](https://github.com/JakobJuergens/RdistNN/issues/16).

- Add on-demand joint CATE covariance from common deletions, including
  retained geometry, candidate views and serialized fits without
  optional weights.

- Draw finite estimates independently of unavailable intervals and
  preserve profile gaps. Support ggplot2 \>= 3.4.1 with focused
  older-version checks.

- Report order fractions and available causal support/contribution
  summaries without new warning thresholds. Clarify CATE migration,
  sampling assumptions and row-order ties; remove six unused suggested
  dependencies.

- Remove the normalized Robinson CATE pseudo-outcome route and its demo
  signals; it is outside the accompanying paper’s theoretical results.
  Use AIPW nuisance inputs or an externally justified supplied signal.

- Clarify accepted-evidence rendering, historical dependency
  reproduction and current-engine validation. Bind release candidates
  and dependency comparisons to source/archive/installed-file
  identities; correct development chronology and local documentation
  news links.

## RdistNN 0.1.0.9003 (merged development snapshot)

The numerical changes below first entered the public production source
in `f2c7bc36e8da040831b668281e50e18b2f1a4b3c`. Version 0.1.0.9003 is a
merged development snapshot; candidate verification, public tags and
downloadable release assets are separate. The existing v0.1.0 release is
unchanged.

- Use ordered compensated Mahalanobis projections so scalar and batch
  targets retain the same cancellation residuals across compiled and
  fallback engines.

- Recover cancellation-sensitive covariance cross terms with compensated
  products and mean-subtraction residuals; test weight mass using the
  portable reduction.

- Preserve scalar tuning sensitivity sweeps across optional return
  fields for NP, CATE, paired-ratio and CLATE fits with a manually
  chosen kernel order.

- Retain custom tuning in experimental moment plots and clarify the
  original-unit Mahalanobis rank cutoff, with independent extreme-scale
  regression references.

- Split shared input resolution, tuning schemas, grid evaluation and
  jackknife orchestration into dedicated modules while preserving public
  interfaces and numerical calculations. Remove unused duplicate
  helpers.

- Strengthen cached-grid target and uncertainty alignment tests.
  Reorganize regression suites by behavior and improve fixture-state
  cleanup.

- Check unique package-level helper ownership, candidate-file hygiene,
  generated Rd/namespace/Rcpp drift, and mirrored native headers. Keep
  the Stata companion on the same development version and refresh the
  documentation entry points.

- Stabilize mixed-sign, large-offset jackknife signals and ordinary
  SD/covariance centering without dropping represented tails. Share
  portable exact rank-table recurrences across C++, R and Mata,
  including subnormal-tail safeguards.

- Preserve data-frame selectors in NP/TDNN/CATE prediction, compute
  requested neighbors independently of SEs, and honor scalar
  tuning/uncertainty controls. Make NP/CATE tidy extraction grid-aware
  and name CATE coefficients by target.

- Reject unrepresentable clipping intervals and partial AIPW nuisance
  inputs; support one-dimensional Mahalanobis geometry and stable
  signed-zero ties. Clarify estimator-specific covariance/LOO contracts.

- Retain resolved fitted geometry across NP/TDNN, CATE and
  paired-ratio/CLATE prediction. Preserve legacy geometry-refitting LOO
  inference. Compensate pure subset-score sums consistently in scalar,
  native and R evaluation.

- Share an explicit compensated weighted-CDF boundary rule with Stata;
  keep selected observations exact across engines. Recover finite
  physical radii and weighted outcome SDs when intermediate squares
  overflow or underflow.

- Reject fractional weight controls and conflicting core input routes;
  mark orders above n infeasible even for point-only grids. Honor
  normalize=FALSE in asymptotic DNN wrappers and correct TDNN
  zero-weight documentation.

- Share compensated rank-weight normalization and weighted signal means
  across native and R paths. Retain correction state through
  prefix/suffix deletion merges, reuse the existing working-buffer
  budget, and align estimator, candidate and saved-prediction
  calculations with the Stata companion.

- Fit Mahalanobis whitening factors directly from the centered training
  design with an explicit singular-value rank cutoff shared with Stata.
  Retain resolvable near-collinear directions, reuse C++/R L2 norms
  after projection, and preserve legacy saved inverse contexts.

- Centralize covariance centering, binary range normalization and
  symmetric assembly across jackknife, NP/TDNN, paired-ratio and
  moment-covariance paths. Recover representable extreme variances;
  report unrepresentable nonzero variances as rescaling errors, or
  existing paired-ratio failure statuses.

## RdistNN 0.1.0.9002 (earlier development snapshot)

- Recover finite fitted SD scales and scaled L1 differences after
  intermediate overflow or underflow; reject nonfinite fitted scales
  before clamping. NP/TDNN prediction standard errors now compute
  marginal variances without allocating a full covariance matrix for the
  prediction grid.

- Stabilize close-order TDNN coefficients with `log1p`/`expm1`, retain
  small asymptotic Beta-increment weights with direct tail
  probabilities, and recover representable extreme L2 norms through
  scaled sum-of-squares evaluation. Share the calculations across
  affected native and fallback paths and retain ordinary-input tie
  ordering and fitted geometry.

- Align R L1/L2 distance fallbacks with compiled subtraction, scaling
  and coordinate accumulation order, preserving stable ranks near
  numerical ties and small differences at large coordinate offsets.
  Scaled evaluation reuses fitted scales with the supplied training
  matrix in both engine modes.

- Add [`summary()`](https://rdrr.io/r/base/summary.html),
  [`coef()`](https://rdrr.io/r/stats/coef.html),
  [`confint()`](https://rdrr.io/r/stats/confint.html),
  [`vcov()`](https://rdrr.io/r/stats/vcov.html) and saved-fit
  [`predict()`](https://rdrr.io/r/stats/predict.html) methods for
  single- and two-scale NP regression, using retained training state and
  full common-deletion target covariance. NP, CATE and paired-ratio
  prediction share target validation; tabular covariate selectors retain
  resolved names so reordered prediction columns remain aligned.

- Promote supplied-prediction conditional LATE and ordinary
  paired-signal ratios to supported workflows, including their shared
  reporting, diagnostics, plots and frozen-model reuse. Pointwise
  inference conditions are unchanged; ratio confidence-set inversion,
  PLQR and local IV-GMM remain experimental.

## RdistNN 0.1.0.9001 (merged development snapshot)

- Advance the shared development version and document local development
  branches, reviewed production PRs and R/Stata parity expectations.

## RdistNN 0.1.0.9000 (merged development snapshot)

- Stabilized frozen-score jackknife calculations with target-specific,
  support-aware numerical centering across compiled, R, fast and
  explicit deletion paths. Zero-weight outliers no longer manufacture
  variance in locally constant signals. Exact and asymptotic weight
  rules are unchanged.

- Reworked the README for package users, with current installation
  guidance, executable CATE, regression and CLATE examples, and shared
  reporting and plotting documentation. Updated the website reference
  descriptions and added contributor navigation.

- Aligned CLATE and generic affine ratios with exact/asymptotic
  single-scale DNN, general `k`, normalization, tuning grids and
  candidate views. Prediction and diagnostics retain the chosen weights
  and frozen paired signals.

- Standardized ratio summaries on the CATE coefficient-table schema,
  retaining component uncertainty and target statuses in separate
  fields. Added CATE-style estimate, weight, instrument-overlap,
  paired-score influence, tuning and dashboard plots, including
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) and
  [`ggplot2::autoplot()`](https://ggplot2.tidyverse.org/reference/autoplot.html)
  dispatch. Shared builders and rank-weight construction support both
  methods; TDNN remains separate NP-regression functionality.

- Added experimental conditional LATE from supplied nuisance predictions
  with
  [`clate_dnn()`](https://JakobJuergens.github.io/RdistNN/reference/clate_dnn.md),
  [`iv_augmented_signals()`](https://JakobJuergens.github.io/RdistNN/reference/iv_augmented_signals.md)
  and the generic
  [`affine_ratio_dnn()`](https://JakobJuergens.github.io/RdistNN/reference/affine_ratio_dnn.md).
  Both input formats retain paired components, covariance scales and
  explicit zero/negative-denominator statuses. Reporting, prediction and
  diagnostics reuse frozen signals; nuisance learning remains external.

- Added
  [`ratio_confidence_set()`](https://JakobJuergens.github.io/RdistNN/reference/ratio_confidence_set.md)
  with bounded, singleton, ray, half-line, all-real, empty and
  unresolved numerical results, plus a bounded
  [`clate_demo_data()`](https://JakobJuergens.github.io/RdistNN/reference/clate_demo_data.md)
  benchmark and internal ridge-series reference adapter. New guides
  explain the estimand, covariance and pointwise inference scope.

- Added independent subset/deletion references, signal identities,
  interface and validation tests, R/compiled parity, ridge refits and
  confidence-set geometry/rescaling tests. Hardened single-observation
  rank shapes, empty CLATE/ratio inputs and near-one confidence
  quantiles.

- Replaced the superseded CATE and IV paper launchers with the current
  affine simulation and SIPP workflows. The bundle pins unchanged paper
  scientific source, all four accepted production configurations, and
  compressed evidence.

- Added portable evidence rendering and validation runs, Linux measured
  production launches with resume checks, and SIPP full-fit/OOB/held-out
  replication with frozen-score and localization diagnostics.

- Updated the replication article and package guidance to distinguish
  the current study from general simulation APIs. Accepted evidence
  regenerates all 21 manuscript tables and six PDF figures without
  refitting.

## RdistNN 0.1.0

### Initial public release scope

- Refocused onboarding, package metadata, and the documentation site on
  CATE and nonparametric regression, with executable examples of both
  input formats, tuning, diagnostics, and the methods each result class
  supports.
- Marked all PLQR functionality as experimental in help pages and
  articles, and grouped its estimator, scores, demos, inference methods,
  and specialized plots in the experimental reference section. PLQR
  continues to evaluate moments at supplied coefficients without solving
  for local coefficients.
- Added contributor guidance, a source-package release checklist, and
  website announcement copy. This polish does not change estimator
  behavior or APIs.

### Paper tuning defaults

- [`tuning_dnn_paper()`](https://JakobJuergens.github.io/RdistNN/reference/tuning_dnn_paper.md)
  now uses the paper’s finite-sample default `lambda = 0.4` and
  sensitivity grid `c(0.3, 0.4, 0.5, 0.75)`, while preserving
  `s(lambda) = floor(lambda * n^((kx + 1)/(kx + 2)))` and all explicit
  tuning overrides.

### Pre-release documentation and metadata

- Added a notation and package-conventions vignette that maps public
  arguments and fitted-object fields to the companion papers,
  distinguishes effective local scale from physical support, and
  separates PLQR moment diagnostics from coefficient-level sandwich
  inference.
- Corrected the
  [`np_dnn()`](https://JakobJuergens.github.io/RdistNN/reference/np_dnn.md),
  [`cate_dnn()`](https://JakobJuergens.github.io/RdistNN/reference/cate_dnn.md),
  [`plqr_dnn()`](https://JakobJuergens.github.io/RdistNN/reference/plqr_dnn.md),
  and
  [`gmm_dnn()`](https://JakobJuergens.github.io/RdistNN/reference/gmm_dnn.md)
  help pages so their automatic-order formula includes the finite-sample
  `lambda = 0.4` multiplier used by the implementation.
- Added the CATE inference contract to the README and command help:
  orthogonality controls first-order population nuisance drift, while
  same-sample reuse separately requires the paper’s nuisance-rate and
  full-state stability conditions.
- Added the current paper link, repository and issue-tracker metadata,
  and explicit build exclusions for local check, archive, and
  rendered-README outputs.
- Expanded release CI to run on default-branch pushes and to cover both
  Apple Silicon and Intel macOS alongside Windows and Linux. Clarified
  that OpenMP is optional at runtime but source installation still
  requires C++17.
- Added `options(RdistNN.use_compiled = FALSE)` as a single switch for
  forcing supported R distance, score, and jackknife fallbacks, with
  compiled/fallback parity tests in the cross-platform matrix.
- Removed local agent settings and paper PDF snapshots from version
  control, eliminated the dormant CI log-commit path, tightened default
  workflow permissions, and added a repository-hygiene guard for
  generated artifacts and finite artifact retention.

### Theory-facing CATE simulation diagnostics

- User-facing `diagnostics(cate_result)` now reuses the paper-simulation
  definitions for the ordinary-Euclidean 90-percent weight-mass radius,
  normalized treated/control arm ESS, and top-five absolute weighted
  score-contribution share, with explicit availability fields.
- The CATE simulation now decomposes DNN error into realized
  localized-branch drift and estimator remainder, and records the
  corresponding root-effective- sample-size standardized drift.
- Weight-based simulation rows now report a physical 90-percent-mass
  radius, treated and control effective sample sizes, score-contribution
  concentration, and the studentized error.
  [`summarise_sim_results()`](https://JakobJuergens.github.io/RdistNN/reference/summarise_sim_results.md)
  adds finite-sample studentized quantiles, skewness, excess kurtosis,
  and support summaries.
- Integrated tuning-grid rows compute these diagnostics separately for
  every candidate while retaining shared-draw and shared-fitted-state
  parity.

### PLQR jackknife-meat inference, supplied-theta semantics, and plot defaults

- [`plqr_dnn()`](https://JakobJuergens.github.io/RdistNN/reference/plqr_dnn.md)
  gains `theta_meat = c("analytic", "jackknife")`: with `"jackknife"`,
  the full delete-1 jackknife covariance of the localized vector moment
  (`moment_vcov`) is used as the coefficient-sandwich meat with the
  density-backed bread. `meta$theta_meat_method` records the route;
  design note 15 now carries the full route table.
- Supplied-theta semantics are explicit across the PLQR S3 methods:
  [`summary()`](https://rdrr.io/r/base/summary.html) prints
  `theta_source` and suppresses `statistic`/`p.value` for supplied
  evaluation values; [`coef()`](https://rdrr.io/r/stats/coef.html),
  [`confint()`](https://rdrr.io/r/stats/confint.html), and
  [`predict()`](https://rdrr.io/r/stats/predict.html) attach a
  `theta_source` attribute and document that current fits repeat the
  supplied `theta_eval`.
- PLQR estimate plots default to moment coordinates:
  [`plot_data_estimates()`](https://JakobJuergens.github.io/RdistNN/reference/plot_data_estimates.md)/[`plot_estimates()`](https://JakobJuergens.github.io/RdistNN/reference/plot_estimates.md)
  accept `quantity = c("moment", "theta")`, zero reference lines are
  drawn for moment plots only, and candidate-grid plots stay in moment
  coordinates.
- `.plqr_direct_sandwich()` and the jackknife-meat route internally
  require normalized weight columns, failing clearly on rescaled or
  corrupted weight matrices.
- PLQR Jacobian diagnostics rename `density_min_weighted` /
  `density_max_weighted` to `density_min_support` /
  `density_max_support` (they are support extrema) and add the weighted
  mean `density_mean_weighted`.
- [`plqr_demo_data()`](https://JakobJuergens.github.io/RdistNN/reference/plqr_demo_data.md)
  now ships a true oracle threshold density: `data$density_oracle` and
  `truth$density_fn`, with the ALD error generator made
  density-continuous at the quantile (`f(0) = 4 tau (1 - tau) / sigma`;
  `tau = 0.5` draws are unchanged). README and vignette examples use it
  instead of a hand-set constant.

### Pre-public breaking cleanup (paper-notation alignment)

This pass deliberately breaks pre-release APIs so the first public
version matches the conditional Z-estimation paper’s notation. No
aliases for the old names are kept.

- [`plqr_dnn()`](https://JakobJuergens.github.io/RdistNN/reference/plqr_dnn.md):
  renamed the shape-nuisance inputs `q`/`q_col` to `g`/`g_col` (the
  paper’s (g_0(U))) and the coefficient input `beta` to `theta_eval`.
  The result stores `theta_eval` (the supplied vector) instead of
  `beta`, and `meta$theta_source = "supplied"` records that no root
  solving was performed.
- PLQR standard errors are no longer overloaded: coefficient-level
  direct-sandwich results live in `theta_se`/`theta_vcov`/`theta_omega`
  (previously `se`/`vcov`/`omega`), while the jackknife moment
  diagnostic lives in `moment_se`. `meta$se_route` was replaced by
  `meta$theta_se_method` and `meta$moment_se_method`.
- Renamed the PLQR return extra `"g"` to `"score"` (stored as
  `out$score`), reserving `g` for the shape nuisance.
- [`plqr_demo_data()`](https://JakobJuergens.github.io/RdistNN/reference/plqr_demo_data.md)
  now returns `g_oracle`, `theta_true`, `g_fn`, and `g_x0` (previously
  `q_oracle`, `beta_true`, `q_fn`, `q_x0`).
- Replaced the demo-bound `plqr_moment()` with the general
  `plqr_score_matrix(data, y_col, D_cols, g_col, mD_cols, theta, tau)`.
- [`score_plqr_component()`](https://JakobJuergens.github.io/RdistNN/reference/score_plqr_component.md)
  gained an explicit `y_col` argument (the outcome column was hardcoded
  as `"Y"`), dropped the ignored compatibility-only `k` argument, and
  now uses `g_col`/`theta` notation throughout.
- Renamed the rank-weight registry functions: `register_rank()` is now
  [`register_rank_weight()`](https://JakobJuergens.github.io/RdistNN/reference/register_rank_weight.md),
  `get_rank()` is
  [`get_rank_weight()`](https://JakobJuergens.github.io/RdistNN/reference/get_rank_weight.md),
  and `rank_names()` is
  [`rank_weight_names()`](https://JakobJuergens.github.io/RdistNN/reference/rank_weight_names.md).
- Direct PLQR coefficient inference now requires `normalize = TRUE`;
  fits combining `density`/`density_col`, `se = TRUE`, and
  `normalize = FALSE` error instead of returning scale-dependent
  standard errors.
- Package `Title`/`Description` now lead with distributional
  nearest-neighbor inference for local conditional moment models rather
  than generalized U-statistics; the README labels the local IV-GMM
  layer as experimental companion functionality.
- Renamed the paper-rule tuning constant from `kappa` to `lambda` to
  avoid a collision with the paper’s own use of kappa:
  `tuning_dnn_paper(lambda =, lambda_grid =)`, candidate-table column
  `lambda`, candidate ids `lambda_*`, diagnostics column
  `tuning_lambda_default`, and the `lambda` column in `sim_cate_*`
  result and summary tables. No `kappa` alias is kept.
- Removed the ignored compatibility-only `k` argument from
  [`score_cate_rho()`](https://JakobJuergens.github.io/RdistNN/reference/score_cate_rho.md)
  and
  [`score_np_y()`](https://JakobJuergens.github.io/RdistNN/reference/score_np_y.md),
  matching the same cleanup in
  [`score_plqr_component()`](https://JakobJuergens.github.io/RdistNN/reference/score_plqr_component.md).
- Tuning-grid diagnostics payloads no longer carry the legacy duplicate
  `top_share` column; `top5_weight_share` is the single concentration
  field, and
  [`diagnostics()`](https://JakobJuergens.github.io/RdistNN/reference/diagnostics.md)
  no longer accepts legacy payloads without it.
- The experimental IV layer no longer aliases its target as tau:
  [`sim_iv_data()`](https://JakobJuergens.github.io/RdistNN/reference/sim_iv_data.md)/[`iv_demo_data()`](https://JakobJuergens.github.io/RdistNN/reference/iv_demo_data.md)
  drop `tau_x0`, `truth$tau_fn`, and the duplicate `tau_oracle` data
  column, keeping `theta_x0`/`theta_fn`/ `theta_oracle` only. `tau`
  remains reserved for the CATE target and the PLQR quantile level.
- Exported-function documentation no longer names internal compiled
  `rustat_*` symbols; compiled fast paths are described generically.
- Pruned the public namespace: the raw distance calculators and order
  helpers (`dist_l1()`, `dist_l2()`, `dist_l1_scaled()`,
  `dist_l2_scaled()`, `order_by_l1()`, `order_by_l2()`,
  `order_by_l1_scaled()`, `order_by_l2_scaled()`), the rank-to-weight
  plumbing (`dnn_weights_*_full_from_ranks()`,
  `tdnn_rank_table_exact()`, `tdnn_rank_table_asymp()`), and
  `ord_to_ranks()` are now internal. Use the `distance_*()` spec
  constructors and the `rank_weight_*()` / `weights_*()` layers instead;
  the internal implementations remain test-covered.

### PLQR inference, prediction, and tuning-grid mode

- Added the direct PLQR coefficient sandwich (`density`/`density_col`):
  [`plqr_dnn()`](https://JakobJuergens.github.io/RdistNN/reference/plqr_dnn.md)
  now reports coefficient-level `theta_se`, `theta_vcov`, `theta_omega`,
  `jacobian`, `score_covariance`, and `neff` when a local conditional
  density is supplied, matching the direct `J, S` sandwich this
  estimator’s source theory specifies (see
  `markdown_notes/design/15_plqr_variance_estimation_rationale.md`).
  [`confint.plqr_result()`](https://JakobJuergens.github.io/RdistNN/reference/confint.plqr_result.md)
  and
  [`summary.plqr_result()`](https://JakobJuergens.github.io/RdistNN/reference/summary.plqr_result.md)
  report coefficient inference only through this route.
- Added
  [`predict.plqr_result()`](https://JakobJuergens.github.io/RdistNN/reference/predict.plqr_result.md),
  including coefficient-level `se.fit`/ `interval = "confidence"`
  support gated on the same direct-sandwich requirement.
- Added `tuning_mode = c("default", "grid")` to
  [`cate_dnn()`](https://JakobJuergens.github.io/RdistNN/reference/cate_dnn.md),
  [`np_dnn()`](https://JakobJuergens.github.io/RdistNN/reference/np_dnn.md),
  and
  [`plqr_dnn()`](https://JakobJuergens.github.io/RdistNN/reference/plqr_dnn.md):
  grid mode evaluates a tuning object’s feasible kernel- order
  candidates and attaches a tidy `tuning_grid` payload (`candidates`,
  `estimates`, `diagnostics`, `meta`), reusing cached ranks and
  rank-weight tables rather than refitting per candidate. See
  `markdown_notes/design/14_integrated_tuning_grid_reuse.md`.
- Recorded, but deliberately did not build, a PLQR analogue of the
  `sim-cate-*.R` paper-simulation harness: the Cond_Z_Estimation paper
  has no PLQR numerical study to mirror, so a simulation harness would
  require inventing a DGP and comparator methods from scratch. Left as
  an explicit, flagged future-work item rather than undertaken
  speculatively.

### Plotting layer (Phase 3)

- `autoplot()`/[`plot()`](https://rdrr.io/r/graphics/plot.default.html)
  methods for `plqr_result` (zero-referenced moment facets with a
  theta_eval/tau subtitle) and `gmm_result` (coefficient facets with the
  model name); zero reference lines are now drawn for PLQR moments as
  well as CATE effects.
- [`plot_conditioning()`](https://JakobJuergens.github.io/RdistNN/reference/plot_conditioning.md)
  for PLQR fits: smallest Gram eigenvalue and log10 condition number
  across targets with caution/ill-conditioned guides at `cond_limits`.
- [`plot_calibration()`](https://JakobJuergens.github.io/RdistNN/reference/plot_calibration.md)
  for PLQR fits: the weighted below-threshold share against the nominal
  `tau` with a normal reference band scaled by the local effective
  sample size; the indicator is stored on new fits and can be supplied
  as `residuals` otherwise.
- [`plot_influence()`](https://JakobJuergens.github.io/RdistNN/reference/plot_influence.md)
  for CATE fits: per-observation weighted score contributions against
  rank or distance, with the largest contributions labelled by
  observation index.
- [`plot_tdnn_scales()`](https://JakobJuergens.github.io/RdistNN/reference/plot_tdnn_scales.md)
  overlays a TDNN fit with internal single-scale refits at `s1` and
  `s2`.
- [`plot_tuning()`](https://JakobJuergens.github.io/RdistNN/reference/plot_tuning.md)
  for `np_result` and `cate_result`: refits across the tuning grid’s
  candidate kernel orders (always including the fitted `s`) and plots
  estimate stability with the fitted order marked.
- Fitted objects store the additional components these panels need:
  training data on
  [`np_dnn()`](https://JakobJuergens.github.io/RdistNN/reference/np_dnn.md)/[`np_tdnn()`](https://JakobJuergens.github.io/RdistNN/reference/np_tdnn.md)
  fits, and the below-threshold indicator on
  [`plqr_dnn()`](https://JakobJuergens.github.io/RdistNN/reference/plqr_dnn.md)
  fits.
- Optional visual regression layer: vdiffr snapshots for the main
  panels, opt-in via the `RDISTNN_VDIFFR=true` environment variable so
  platform font differences cannot break CI.

### Plotting layer (Phase 2)

- [`plot_overlap()`](https://JakobJuergens.github.io/RdistNN/reference/plot_overlap.md)
  for CATE fits: `type = "mass"` draws the weighted treated/control mass
  per target with dotted share guides and shaded warning bands outside
  `share_limits` (violations reported in the caption);
  `type = "propensity"` draws weighted propensity densities at selected
  targets when observation-level propensities are available and a
  weighted q10-q90 summary otherwise. Inputs resolve through the
  documented precedence chain (stored fit components, then
  `treatment`/`propensity` arguments, then an actionable error).
- `autoplot()`/[`plot()`](https://rdrr.io/r/graphics/plot.default.html)
  on
  [`diagnostics()`](https://JakobJuergens.github.io/RdistNN/reference/diagnostics.md)
  frames: effective sample size (with the degenerate-weights guide at
  `0.9 n`), maximum weight, and top-5 weight share faceted across
  targets.
- [`plot_dashboard()`](https://JakobJuergens.github.io/RdistNN/reference/plot_dashboard.md)
  for `cate_result` and `np_result`: patchwork composition with a
  metadata annotation (`n`, `m`, `k`, `s` with provenance, distance,
  jackknife type); unavailable panels are substituted (diagnostics
  facets replace a missing overlap panel) or dropped, and without
  patchwork the panels are returned as a named list with a message.
  patchwork joins Suggests.
- Fitted
  [`np_dnn()`](https://JakobJuergens.github.io/RdistNN/reference/np_dnn.md),
  [`cate_dnn()`](https://JakobJuergens.github.io/RdistNN/reference/cate_dnn.md),
  [`plqr_dnn()`](https://JakobJuergens.github.io/RdistNN/reference/plqr_dnn.md),
  and
  [`np_tdnn()`](https://JakobJuergens.github.io/RdistNN/reference/np_tdnn.md)
  objects record the jackknife type in `meta$jackknife` (`"none"` when
  `se = FALSE`).
- New vignette `visualising-fits` walks from a CATE fit through the
  estimate, weight, overlap, and diagnostics panels to the dashboard and
  panel customisation.
- The agent instruction files are excluded from the build, removing the
  long-standing top-level-files NOTE in R CMD check.

### Plotting layer (Phase 1)

- Added ggplot2-based plotting (ggplot2 is a suggested dependency; all
  entry points fail with an actionable message when it is missing):
  [`plot_estimates()`](https://JakobJuergens.github.io/RdistNN/reference/plot_estimates.md)
  draws point estimates with normal confidence bounds as a profile band
  (targets varying along one stored coordinate), a forest layout (small
  non-profile target sets, labelled from the stored coordinates), or
  pointranges, with an honesty caption carrying the median effective
  sample size and the kernel-order provenance on every panel; CATE
  panels add a dashed zero reference line and an opt-in
  `shade_significant` rug that marks intervals, never point estimates.
- [`plot_weights()`](https://JakobJuergens.github.io/RdistNN/reference/plot_weights.md)
  draws weight decay by rank (support-edge guide at `n - s + 1`,
  optional log scale), cumulative weight shares (0.5/0.9 guides and
  top-share markers), and the weight-vs-distance “effective kernel” view
  (rank-`s` radius guide); two-scale TDNN panels add a zero line for
  signed weights, and large target sets are subset to
  first/median-`n_eff`/last by default.
- `autoplot()` (registered on ggplot2) and
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) methods for
  `np_result` and `cate_result` delegate to
  [`plot_estimates()`](https://JakobJuergens.github.io/RdistNN/reference/plot_estimates.md).
- README gains a visualisation section with the CATE profile figure.

### Plotting data layer (Phase 0)

- Added the tidy data extractors behind the upcoming ggplot2-based
  plotting interface, with no plotting dependency:
  [`plot_data_estimates()`](https://JakobJuergens.github.io/RdistNN/reference/plot_data_estimates.md)
  (S3 generic with methods for `np_result`, `cate_result`,
  `plqr_result`, and `gmm_result`; TDNN fits dispatch through
  `np_result`) returns one row per target and coefficient with
  estimates, stored jackknife standard errors, and normal confidence
  bounds at a chosen `level`;
  [`plot_data_weights()`](https://JakobJuergens.github.io/RdistNN/reference/plot_data_weights.md)
  returns per-(target, rank) weight profiles with cumulative shares and
  distances when stored;
  [`plot_data_diagnostics()`](https://JakobJuergens.github.io/RdistNN/reference/plot_data_diagnostics.md)
  returns the per-target diagnostics frame with an integer `target_id`
  aligned with the other extractors.
- Maintainer architecture:
  `markdown_notes/design/13_estimator_and_result_architecture.md`.

### Fitted objects store evaluation points

- [`np_dnn()`](https://JakobJuergens.github.io/RdistNN/reference/np_dnn.md),
  [`cate_dnn()`](https://JakobJuergens.github.io/RdistNN/reference/cate_dnn.md),
  [`plqr_dnn()`](https://JakobJuergens.github.io/RdistNN/reference/plqr_dnn.md),
  [`gmm_dnn()`](https://JakobJuergens.github.io/RdistNN/reference/gmm_dnn.md),
  and
  [`np_tdnn()`](https://JakobJuergens.github.io/RdistNN/reference/np_tdnn.md)
  result objects now include `x0`, the `m x p` matrix of evaluation
  points, with covariate column names propagated from the training data
  when the targets carry none. Profile plots need the target
  coordinates, which were previously not stored.

### Breaking change: default kernel order `s`

- [`np_dnn()`](https://JakobJuergens.github.io/RdistNN/reference/np_dnn.md),
  [`cate_dnn()`](https://JakobJuergens.github.io/RdistNN/reference/cate_dnn.md),
  [`plqr_dnn()`](https://JakobJuergens.github.io/RdistNN/reference/plqr_dnn.md),
  and
  [`gmm_dnn()`](https://JakobJuergens.github.io/RdistNN/reference/gmm_dnn.md)
  no longer default to `s = 1L`, which produced exactly uniform DNN
  weights and therefore no localization (every target point silently
  received the global weighted mean). When both `s` and `tuning` are
  omitted, `s` is now resolved from the paper undersmoothing rule
  [`tuning_dnn_paper()`](https://JakobJuergens.github.io/RdistNN/reference/tuning_dnn_paper.md),
  `s = floor(0.4 * n^((kx + 1)/(kx + 2)))`, matching the current default
  used by the Stata companion (`cate_dnn.ado`). The resolved rule is
  recorded in `meta$tuning` with `source = "paper_default"`. Calls that
  passed `s` explicitly are unaffected.
- All DNN and TDNN user wrappers now emit a warning when the computed
  weights are (near-)uniform (max effective sample size at least
  `0.9 * n`, for example `s = 1` or `k = s`), instead of silently
  returning an unlocalized fit.
- Added `tests/testthat/test-s-default-paper-rule.R`, the `s`-analogue
  of the `k = 1` documentation guardrail: it pins the default-`s`
  resolution for all four wrappers, the manual-`s` path, the
  degenerate-weight warning, and R/Stata default-formula agreement.
- [`np_tdnn()`](https://JakobJuergens.github.io/RdistNN/reference/np_tdnn.md)
  still uses explicit `s1`/`s2` defaults; a two-scale default rule is a
  separate follow-up.

### Diagnostics and reporting

- Added
  [`diagnostics()`](https://JakobJuergens.github.io/RdistNN/reference/diagnostics.md)
  methods for
  [`np_dnn()`](https://JakobJuergens.github.io/RdistNN/reference/np_dnn.md),
  [`cate_dnn()`](https://JakobJuergens.github.io/RdistNN/reference/cate_dnn.md),
  and
  [`plqr_dnn()`](https://JakobJuergens.github.io/RdistNN/reference/plqr_dnn.md)
  fitted objects, including generic DNN weight/tuning summaries plus
  estimator-specific NP outcome, CATE overlap, and PLQR Gram
  conditioning diagnostics.
- Added a diagnostics/reporting vignette and updated README,
  getting-started, conditional Z-estimation, kernel-order tuning, and
  pkgdown references to show the diagnostics workflow.

### CI and package engineering

- Added `tools/ci/check_tests_declared_deps.R` and a CI guard step to
  fail when `tests/` uses `pkg::` or `pkg:::` for undeclared packages.
- Added `tools/ci/check_pkgdown_reference_coverage.R` and a CI guard
  step to fail when public `man/*.Rd` topics are missing from
  `_pkgdown.yml` reference contents.
- Added `tools/ci/check_generated_file_drift.R` and a CI guard step to
  fail when generated files (`NAMESPACE`, `R/RcppExports.R`,
  `src/RcppExports.cpp`) are stale.
- Added workflow `.github/workflows/alignment-guardrails.yaml` to run
  these guardrail checks on push and pull requests.

### Dependency and docs alignment

- Added `tibble` to `Suggests` for test-time namespace declarations.
- Filled missing `_pkgdown.yml` reference entries for TDNN and
  Mahalanobis topics that were failing site checks.
