# Plotting Layer Implementation Plan

Created: 2026-07-02.
Status: planning note; no plotting code exists yet.  This is the working
document for the implementation and should be updated as phases land.

## 1. Purpose and Motivation

The package currently has no `plot()`/`autoplot()` methods.  For the
package's target audience the product *is* the picture: a CATE profile with
a confidence band, a weight-decay profile with `n_eff` annotated, and an
overlap diagnostic are what a seminar, referee, or applied user looks at
first.  The plotting layer closes the largest remaining gap between
"researcher's lab bench" and "applied user's instrument" (see the
package-review findings of 2026-07-02), and it produces the demo assets for
the job-market paper's talks and README.

Technology decision (author preference, confirmed): **ggplot2** for all
plots, **patchwork** for composition.  No base-graphics fallback.

## 2. Design Decisions

### 2.1 Dependency policy

Per `design/06_namespace_and_dependency_policy.md`, heavy dependencies stay
out of `Imports`:

- `ggplot2`, `patchwork` → `Suggests`.
- `vdiffr` → `Suggests` (visual regression tests only).
- No `scales`, `broom`, or `ggtext` dependencies in phase 1–3; use
  ggplot2-native scales and plain labels.  Revisit `broom` integration
  (§10) only after the layer stabilizes.

Every plotting entry point guards with `requireNamespace("ggplot2",
quietly = TRUE)` and fails through the standardized error helpers
(`design/05_error_messages_and_warnings.md`) with an actionable message
("install.packages('ggplot2')").  Dashboard functions degrade gracefully
when `patchwork` is missing: return a named list of ggplots with a
one-time message instead of a composed figure.

### 2.2 Generic surface

- Implement **`autoplot()` methods** (ggplot2's generic) for the fitted
  classes, and **`plot()` methods that delegate to `autoplot()`** and
  return the object invisibly after printing.  This is the established
  ggplot2-ecosystem convention and keeps both `plot(fit)` (casual user)
  and `autoplot(fit) + theme_bw()` (power user) working.
- Implement a small **`plot_*()` family** of explicitly named builders for
  the individual panels, exported per the export discipline in design/06
  (each has a stable user story and public-API tests):
  - `plot_estimates()` — estimates with CIs across target points;
  - `plot_weights()` — weight profiles;
  - `plot_overlap()` — CATE treated/control weighted mass;
  - `plot_tuning()` — estimate stability across the `s` candidate grid;
  - `plot_dashboard()` — patchwork composition of the above.
- All builders return ggplot (or patchwork) objects; they never print.

### 2.3 Data layer first

Each plot is a thin ggplot over an exported, documented **data extractor**
returning a plain `data.frame`:

- `plot_data_estimates(fit)` → one row per target × coefficient:
  `target_id`, `coef` (for vector-valued fits), `estimate`, `se`,
  `ci_lower`, `ci_upper` (given `level`), plus target coordinates
  `x0_1..x0_p` when available.
- `plot_data_weights(fit)` → one row per (rank, target): `target_id`,
  `rank`, `weight`, `cum_share`, and `distance` when `dx` was returned.
- `plot_data_diagnostics(fit)` → the per-target frame already computed by
  `.dnn_weight_diagnostics()` (`n_eff`, `max_weight`, `top1/5/10`,
  `n_nonzero`, rank radii, weighted distances) in tidy form.

Rationale: users who dislike our defaults can build their own plots from
the same frames; tests can assert on data frames without rendering; and
the Stata companion can later mirror the *data* contract even though it
cannot mirror ggplot.

### 2.4 Prerequisite: store target coordinates in fitted objects

Profile plots need the evaluation points.  Currently `*_result` objects
carry `meta` (`n`, `m`, `p`, `k`, `s`, `method`, `distance`) but **not**
`X0`.  Phase 0 therefore adds `x0` (the `m x p` target matrix) to the
result objects of all wrappers (`np_dnn`, `cate_dnn`, `plqr_dnn`,
`gmm_dnn`, `tdnn_np`).  Memory cost is trivial (`m x p` doubles).  This is
a non-breaking additive change; update the `@return` docs and the
structure tests accordingly.

### 2.5 Axis semantics for multi-dimensional targets

Rule ("profile detection"): if the target matrix varies in **exactly one
column**, that coordinate is the x-axis (labelled with the column name
when available, else `x[j]`).  Otherwise default to `target_id` (ordered
index) with a caption noting that targets vary in several coordinates,
and accept an explicit `x_var =` argument (column index/name) to override.
Never silently pick "the first varying coordinate" when several vary.

## 3. Plot Inventory by Class

### 3.1 `np_result` / `tdnn_np_result`

- **Estimates panel**: `mu` vs. profile axis; `geom_ribbon` CI band when
  `se` present and targets form a profile, `geom_pointrange` otherwise.
- **Weights panel**: weight vs. rank (line, one facet or colour per
  target when `m` small; summarize to selected targets when `m` large),
  with `n_eff` annotated per target and a vertical guide at rank
  `n - s + 1` (support edge).
- TDNN: overlay the two scales' weight profiles; negative weights make
  `geom_hline(yintercept = 0)` mandatory.

### 3.2 `cate_result`

- **Estimates panel**: `tau` with CI band/pointrange (the flagship plot;
  this is the README image).
- **Overlap panel**: per-target weighted treated mass `sum(w * T)` and
  control mass side by side (stacked or dodged bars for few targets, two
  lines for profiles), with a reference band marking severe imbalance;
  reuses the quantities from `diagnostics.cate_result()`.
- **Weights panel**: as in 3.1.
- **Dashboard**: estimates / overlap / weights composed via patchwork
  with `plot_annotation()` carrying `n`, `s` (and its provenance:
  `paper_default` vs manual), `k`, distance, and jackknife type.

### 3.3 `plqr_result`

- **Estimates panel**: coefficient-wise pointrange, faceted by
  coefficient across targets (profile axis when detected).
- **Conditioning panel**: smallest eigenvalue / condition number of the
  weighted Gram matrix per target (already in the PLQR diagnostics).

### 3.4 `gmm_result`

- **Estimates panel**: coefficient-wise pointrange (as PLQR).
- Overidentification visuals (local `J` profiles / `J_n(x)` maps) are
  **out of scope** here; they belong to the Cond_GMM paper timeline and
  get their own plan when that paper's inference layer lands.

### 3.5 `dnn_diagnostics`

- `autoplot()` on the diagnostics object itself: `n_eff`, `max_weight`,
  and top-share panels across targets — cheap, since the frame exists.

### 3.6 Tuning sensitivity (`plot_tuning`)

For a fitted object (or `tuning_eval()` grid), refit at each candidate `s`
in `kappa_grid` and plot estimate ± CI against `s`, marking the selected
default.  This directly answers the "how sensitive is this to `s`?"
seminar question.  Refitting is the expensive part: accept a fitted
object plus data, or precomputed fits; document cost; cap default grid at
the four built-in kappas.

## 4. Visual Style

- All panels use a package-level `.rdistnn_theme()` built on
  `theme_minimal(base_size = 11)`; do **not** export a theme, and apply it
  only inside our builders so user themes can override by `+`.
- Colour: single-hue defaults; treated/control uses a colour-blind-safe
  pair (Okabe–Ito palette hard-coded, no dependency).
- CI level default `0.95`, argument `level` everywhere; label axes with
  units-neutral names (`"Estimate"`, `"DNN weight"`, `"Rank"`), titles off
  by default (dashboards set them via patchwork annotation).
- No plot ever drops the `n_eff` information silently: estimate panels
  carry a caption with `median n_eff` (and range) across shown targets.
  This keeps the package's honesty-first stance visible in every figure.

## 5. Errors, Edge Cases, Guards

- Missing ggplot2 → standardized `.abort_missing_pkg()`-style error; the
  message names the function and the package to install.
- `se = FALSE` fits: estimate panels draw points/lines without bands and
  message once ("no standard errors stored; rerun with se = TRUE").
- Degenerate weights (the `s = 1`/`k = s` warning path): weight panels
  must still render (flat profile) — the picture *is* the diagnostic.
- Single target (`m = 1`): estimates panel becomes a one-point
  pointrange; dashboard remains valid.
- Large `m` (say > 50): weight panel defaults to a target subsample
  (first/median/last plus `targets =` override) to avoid overplotting.
- `W` not stored (`return` without `"W"`): `plot_weights()` errors with
  the exact `return = c("W", "dx")` remedy; `plot_data_weights()` same.

## 6. Testing Policy

Per `design/09_testing_and_regression_policy.md` and the CI guard on
declared test dependencies:

- **Structural tests (always run)**: extractor frames have the promised
  columns/types/rows on all five classes; builders return objects
  inheriting `"ggplot"` (or `"patchwork"`); guard errors fire without
  ggplot2 (simulated via mocked `requireNamespace` where practical);
  profile detection picks the right axis; `k = 1` in every example
  (the k-1 documentation guardrail scans plot docs too).
- **Visual regression (optional layer)**: `vdiffr::expect_doppelganger`
  for one canonical figure per builder, `skip_if_not_installed`, and
  skipped on CRAN.  Baselines regenerate only via an explicit local run.
- Cross-field consistency: dashboard annotation values must equal
  `meta` fields (test compares rendered labels' data, not pixels).

## 7. Documentation and Promotion

- New vignette `visualising-fits.Rmd`: one narrative from `cate_dnn()`
  fit → dashboard → customizing panels (`autoplot(fit) + theme_bw()`),
  guarded by `requireNamespace` chunks per design/11.
- Extend `diagnostics-reporting.Rmd` with the `autoplot(diagnostics(fit))`
  hook.
- README gains the CATE profile figure near the top (pre-rendered PNG via
  README.Rmd so CRAN/GitHub render without Suggests installed).
- pkgdown: new "Visualisation" reference section (`_pkgdown.yml` must list
  every exported topic — the coverage guard enforces this).
- NEWS entry per phase.

## 8. Phased Rollout (acceptance gates per phase)

- **Phase 0 — data plumbing.**  Store `x0` in all result objects; export
  `plot_data_estimates()`, `plot_data_weights()`, `plot_data_diagnostics()`
  with tests.  Gate: structural tests green; no ggplot2 dependency yet.
  *Status 2026-07-02: Phase 0 complete.  (a) `x0` storage — all five
  wrappers (`np_dnn`, `cate_dnn`, `plqr_dnn`, `gmm_dnn`, `np_tdnn`) store
  the `m x p` target matrix via `.dnn_x0_with_names()` with covariate
  names propagated from training data (user-supplied `X0` names win);
  tests in `test-result-x0-target-matrix.R`.  (b) Extractors in
  `R/plot-data.R`: `plot_data_estimates()` (S3 generic, methods for the
  four base classes; TDNN dispatches via `np_result`; `level` argument;
  `x0_1..x0_p` columns with `attr(., "x0_names")` and `attr(., "level")`),
  `plot_data_weights()` (distance-based ranks when `dx` stored, weight-
  ordered otherwise; actionable error when `W` absent), and
  `plot_data_diagnostics()` (integer `target_id`, original label kept in
  `target_label`).  Tests in `test-plot-data.R`; pkgdown section
  "Plotting data layer" added; full suite 2295 passing and all three CI
  guardrails green.  Note: the TDNN wrapper is `np_tdnn()`, not
  `tdnn_np()` as earlier drafts of this note said.  Phase 1 (autoplot/
  plot + estimates/weights panels, ggplot2 in Suggests) is next.*
- **Phase 1 — core panels.**  `autoplot()`/`plot()` for `np_result` and
  `cate_result` (estimates + weights), `plot_estimates()`,
  `plot_weights()`; Suggests updated; guard errors tested.  Gate: full
  suite green; README figure generated.
- **Phase 2 — diagnostics and composition.**  `plot_overlap()`,
  `autoplot.dnn_diagnostics()`, `plot_dashboard()` with patchwork
  degradation path; dashboard annotation contract test.  Gate: vignette
  draft renders.
- **Phase 3 — remaining classes + tuning.**  `plqr_result`, `gmm_result`,
  TDNN overlays, `plot_tuning()`.  Gate: pkgdown section complete;
  vdiffr baselines recorded.
- **Phase 4 — polish.**  README/vignette finalized; NEWS consolidated;
  Stata parity decision recorded (§9).

Each phase is a separate commit series with the standard verification
gate (full testthat suite; pkgdown coverage check when exports change).

## 9. Out of Scope (recorded deliberately)

- Stata graphical parity: at most a `twoway`-based estimates plot much
  later; the *data contract* (§2.3) is the cross-language surface for now.
- Local `J_n(x)` overidentification maps (Cond_GMM paper timeline).
- Simulation-layer figures (`sim_cate_*` coverage/bias plots): these are
  paper-replication assets, not user API; they can consume the same
  extractors from `inst/paper-replication/` scripts without new exports.
- Interactive (plotly/shiny) output.

## 10. Open Decisions (resolve before or during Phase 1)

1. **`plot(fit)` default**: single estimates panel (proposed) vs. full
   dashboard.  Proposal: estimates panel; dashboards are explicit
   (`plot_dashboard()`), so casual `plot()` stays fast and legible.
2. **CI type**: plain Wald bands from stored `se` only, or also expose
   `level` recomputation from jackknife internals?  Proposal: stored `se`
   with `level` argument; nothing recomputed.
3. **broom methods** (`tidy.cate_result` etc.): defer; `plot_data_*`
   covers the need without a new dependency surface.  Revisit if users
   ask.
4. **Naming**: `plot_dashboard()` vs `plot_summary()`; proposal
   `plot_dashboard()` to avoid clashing with `summary()` semantics.
5. Whether `plot_tuning()` belongs in Phase 3 or ships earlier as the
   public answer to the `s`-sensitivity question (it is the most
   seminar-valuable panel after the CATE profile).
