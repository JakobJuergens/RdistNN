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
  (each has a stable user story and public-API tests).  Where one builder
  hosts several closely related views, a `type =` argument selects the
  view instead of multiplying exports:
  - `plot_estimates()` — estimates with CIs across target points
    (profile / pointrange / forest variants, §3.0);
  - `plot_weights(type = c("profile", "cumulative", "kernel"))` — weight
    decay by rank, cumulative concentration, and weight-vs-distance
    ("effective kernel") views;
  - `plot_overlap(type = c("mass", "propensity"))` — CATE treated/control
    weighted mass and local weighted propensity summaries;
  - `plot_influence()` — per-observation weighted score contributions
    (CATE first; NP analogue later);
  - `plot_calibration()` — PLQR local quantile calibration;
  - `plot_conditioning()` — PLQR weighted-Gram conditioning;
  - `plot_tuning()` — estimate stability across the `s` candidate grid;
  - `plot_dashboard()` — patchwork composition per class (§3.1–3.3).
- All builders return ggplot (or patchwork) objects; they never print.
- **Sample-level inputs precedence.**  Panels that need observation-level
  data beyond the stored fit (raw outcome scatter, propensity values,
  score contributions) resolve inputs in a fixed order: (1) components
  stored on the fit (e.g. `model$t_vec`, `model$e` on data-interface CATE
  fits), then (2) stored `return` extras (`psi`, `nuisance`, `g`, `dx`),
  then (3) explicit user-supplied arguments (`X`, `Y`, `T`, `e`), and
  otherwise (4) an actionable error naming the `return =` or argument
  remedy.  No panel silently recomputes first stages.

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

## 3. Plot Inventory: Precise Panel Specifications

Panels are identified as `NP-*` (nonparametric regression), `C-*` (CATE),
`P-*` (PLQR), and `G-*` (GMM).  Each spec states the question the panel
answers, its data source in extractor terms, the exact geometry, and edge
cases.  Data sources refer to the Phase-0 extractors
(`plot_data_estimates()`, `plot_data_weights()`,
`plot_data_diagnostics()`) and the confirmed diagnostics columns
(`treated_mass`, `control_mass`, `treated_share`,
`weighted_mean_propensity`, `weighted_propensity_q10/q90`,
`gram_min_eig`, `gram_max_eig`, `gram_condition_number`, plus the generic
`n_eff`, `max_weight`, `top1/5/10_weight_share`, `radius_rank_k/s`,
`weighted_mean/median_distance`, `tuning_*` fields).

### 3.0 Shared panel grammar

Rules applied by every builder, so individual specs below only state
deviations:

- **Estimate display mode** (from `plot_data_estimates()`): if targets
  form a profile (§2.5) **and** `m >= 4`: `geom_line` + `geom_ribbon`
  (CI band, `alpha = 0.25`) over the profile axis.  If `m < 4` or targets
  are not a profile: `geom_pointrange` over `target_id` (or forest
  layout, §3.2 C-2).  If `m = 1`: single `geom_pointrange`.
- **CI handling**: bands/ranges from `ci_lower`/`ci_upper`; when `se` is
  all `NA`, draw the point/line layer only and add the caption line
  "no standard errors stored (se = FALSE)".
- **Honesty caption**: every estimate-bearing panel appends
  "median n_eff = X (min–max across targets); s = S (paper_default|manual)"
  built from `plot_data_diagnostics()` and `meta$tuning$source`.
- **Multi-coefficient fits** (PLQR, GMM): `facet_wrap(~ coef,
  scales = "free_y")`, coefficient order as returned by the extractor.
- **Reference lines**: zero lines (`geom_hline(yintercept = 0,
  linetype = "dashed")`) are mandatory wherever zero is meaningful
  (CATE effects, PLQR moments, TDNN weights) and omitted elsewhere
  (NP conditional means).
- **Weight-domain scales**: rank-domain panels use a linear y by default;
  `log_weights = TRUE` switches to `scale_y_log10` after dropping exact
  zeros with a caption noting how many were dropped.  Distance-domain
  panels always use linear axes.
- **Many targets**: weight-domain panels default to a target subset
  (first, median-`n_eff`, last; `targets =` overrides) once `m > 8`;
  colour maps `factor(target_id)` with an `x0`-based legend label when a
  profile axis exists.
- **Palette**: Okabe–Ito, hard-coded; treated = "#0072B2" (blue),
  control = "#D55E00" (vermillion); single-series panels use near-black.

### 3.1 Nonparametric regression (`np_result`, `tdnn_np_result`)

- **NP-1 — Fit profile** (flagship; `plot_estimates()` /
  `autoplot.np_result`).  *Question: what does the fitted conditional
  mean look like, and how uncertain is it?*  Data: estimates frame
  (`coef == "mu"`).  Geometry per §3.0.  Optional raw-data underlay:
  when the caller supplies `X`/`Y` (precedence chain §2.2) and a profile
  axis exists, add `geom_point(alpha = 0.25, colour = "grey55")` of `Y`
  against the matching covariate *before* the fit layer — the classic
  nonparametric-regression picture.  Off-profile covariates of the
  sample are *not* filtered; the caption notes that the scatter ignores
  the non-varying coordinates (pointwise honesty, not a curve-fitting
  claim).
- **NP-2 — Effective kernel** (`plot_weights(type = "kernel")`).
  *Question: what kernel did the DNN weights effectively place at this
  target?*  Data: weights frame with `distance` (requires `dx`).
  Geometry: `geom_step(aes(distance, weight))` per target; vertical
  dashed guide at `radius_rank_s` (the s-th neighbour distance, from the
  diagnostics frame) labelled "rank-s radius"; `n_eff` in the subtitle.
  For TDNN fits the zero hline is mandatory (§3.0) and the caption warns
  that two-scale weights are signed.
- **NP-3 — Weight decay and concentration**
  (`plot_weights(type = "profile")` and `type = "cumulative"`).
  *Questions: how fast do weights decay in rank; how concentrated is the
  estimate?*  Profile view: `geom_step(aes(rank, weight))`, vertical
  guide at rank `n - s + 1` (support edge, exact zeros beyond).
  Cumulative view: `geom_line(aes(rank, cum_share))`, horizontal guides
  at 0.5 and 0.9 with the crossing ranks annotated ("half the weight
  from the nearest r* points"), and points marking `top1/5/10` shares
  from the diagnostics frame for cross-checking.
- **NP-4 — n_eff profile** (`autoplot.dnn_diagnostics`, first facet).
  *Question: where in covariate space does the fit have effective
  data?*  Data: diagnostics frame.  Geometry: `geom_line`/`geom_point`
  of `n_eff` over the profile axis (or `target_id`); companion facets
  for `max_weight` and `top5_weight_share`.  A horizontal guide at
  `0.9 * n` marks the degenerate-weights warning threshold so the plot
  and the runtime warning tell one story.
- **NP-5 — TDNN scale decomposition** (`plot_tdnn_scales()`; Phase 3).
  *Question: what does the two-scale combination do relative to its
  single-scale components?*  Implementation decision: the builder refits
  `np_dnn` at `s1` and `s2` internally (cheap; same `X`, `Y` via the
  precedence chain) and overlays three NP-1 layers: DNN(s1), DNN(s2)
  (thin, muted), TDNN (bold), with a legend keyed by estimator.  No new
  storage on the fit object.
- **NP-6 — s-sensitivity** (`plot_tuning()`).  *Question: how much does
  the estimate move across admissible kernel orders?*  Data: refits at
  each `s` in the tuning grid (default: the four built-in kappas; accept
  precomputed fits to avoid recomputation).  Geometry: pointrange of
  estimate ± CI against `s` (x continuous, breaks at candidate values),
  the selected default marked by a filled point and vertical dotted
  guide, `n_eff` per candidate as secondary text labels.  One target per
  panel; multiple targets facet.

### 3.2 CATE (`cate_result`)

- **C-1 — CATE profile** (flagship; `plot_estimates()` /
  `autoplot.cate_result`; the README image).  *Question: how does the
  treatment effect vary along the profile, and where is it
  distinguishable from zero?*  Data: estimates frame (`coef == "tau"`).
  Geometry per §3.0 plus the mandatory zero hline; optional
  `shade_significant = TRUE` (default `FALSE`) underlines profile
  segments whose CI excludes zero with a subtle `geom_rug`-style bar
  rather than colour fills (no traffic-light styling of point
  estimates).  Caption carries the §3.0 honesty line.
- **C-2 — Forest view** (`plot_estimates(style = "forest")`; automatic
  when `m <= 12` and targets are not a profile).  *Question: how do
  effects compare across a set of named covariate profiles?*  Rows are
  targets ordered by `target_id`, labelled by compact `x0` descriptions
  built from `x0_*` columns and `attr(., "x0_names")` (e.g.
  "X1 = 0, X2 = 1"); `geom_pointrange` horizontal, zero vline.
- **C-3 — Overlap mass** (`plot_overlap(type = "mass")`).  *Question: do
  both treatment arms carry weight near the target?*  Data: diagnostics
  frame columns `treated_mass`, `control_mass`, `treated_share`.
  Geometry: for profiles, two lines (treated/control, Okabe–Ito pair)
  of arm mass over the profile axis; otherwise dodged bars per target.
  A shaded warning band where `treated_share` leaves `[0.1, 0.9]`
  (configurable `share_limits =`), matching the overlap gate language in
  the paper.  Unavailable overlap columns (matrix-interface fit without
  stored treatment) produce the precedence-chain error naming
  `t_col`/`T` remedies.
- **C-4 — Local propensity summary** (`plot_overlap(type =
  "propensity")`).  *Question: what does the weighted propensity
  distribution look like at the target?*  Data: diagnostics columns
  `weighted_mean_propensity`, `weighted_propensity_q10/q90` (interval
  plot over targets: pointrange mean with q10–q90 whiskers, guides at
  the clipping bounds); when observation-level `e` is available
  (stored `model$e` or user-supplied), add a per-target weighted density
  overlay (`geom_density(aes(weight = w))`, one facet per selected
  target) — the richer view supersedes the summary whiskers in that
  facet.
- **C-5 — Influence / contribution** (`plot_influence()`; Phase 2–3).
  *Question: is the estimate driven by a handful of observations?*
  Data: weights frame joined with the stored pseudo-outcome (`psi`
  extra; else user-supplied): contribution
  `w_i * (psi_i - tau_hat)` against rank (or distance when available),
  `geom_point` with the top three |contribution| observations labelled
  by `obs_id`.  Caption ties to `max_weight`/`top1_weight_share`
  diagnostics.  This is the plot-level counterpart of the delete-one
  jackknife's sensitivity logic.
- **C-6 — CATE dashboard** (`plot_dashboard.cate_result`).  patchwork
  layout `(C-1 | C-3) / (C-4 | NP-3 cumulative)` with
  `plot_annotation()` carrying `n`, `m`, `k`, `s` + provenance,
  distance, jackknife type; panels degrade per availability (a missing
  overlap panel is replaced by NP-4 rather than an empty slot).

### 3.3 PLQR (`plqr_result`)

The current wrapper evaluates localized PLQR moments at a supplied
`beta`; panels are phrased for that object and will transfer unchanged
to the future estimation wrapper (whose coefficient profiles then reuse
the §3.0 estimate grammar like GMM).

- **P-1 — Moment profile** (`plot_estimates()` /
  `autoplot.plqr_result`).  *Question: how far are the localized
  orthogonal moments from zero at this beta?*  Data: estimates frame
  (one `coef` level per moment coordinate).  Geometry per §3.0 with
  coefficient facets and the mandatory zero hline; the subtitle names
  `beta` and `tau`.  Reading: a specification/centering check — moments
  within their CIs of zero are consistent with `beta` at these targets.
- **P-2 — Gram conditioning** (`plot_conditioning()`).  *Question: is
  the local residualized design well-conditioned enough to trust the
  moment evaluation (and later, inversion)?*  Data: diagnostics columns
  `gram_min_eig`, `gram_condition_number`.  Geometry: two facets over
  targets — `gram_min_eig` (linear y, guide at 0) and
  `gram_condition_number` (`scale_y_log10`, dashed guides at 1e3/1e6
  labelled "caution"/"ill-conditioned", `cond_limits =` configurable).
  This panel is the software face of the paper's local-design
  conditioning diagnostics.
- **P-3 — Local quantile calibration** (`plot_calibration()`; Phase 3).
  *Question: does the weighted share of below-threshold residuals match
  tau locally?*  Data: requires the stored score matrix (`g` extra) or
  user-supplied residual indicators; per target compute the weighted
  indicator share and plot as points over targets with a horizontal
  line at `tau` and a normal-approximation reference band
  `tau ± z * sqrt(tau(1-tau)/n_eff)`.  Deviations flag local
  misspecification of `beta`/first stages — the applied reading of the
  quantile moment condition.
- **P-4 — Threshold band mass** (optional, Phase 3+, theory-facing).
  *Question: how much weighted mass sits in a shrinking band around the
  quantile threshold?*  Weighted mass of `|residual| <= h` against a
  small grid of `h`, one line per target.  This is the fitted-band
  weighted-square object from the paper's nonsmooth ScoreCov route made
  visible; it is deliberately labelled as a research diagnostic and can
  live in the paper-replication material first if the API surface should
  stay lean.
- **P-5 — PLQR dashboard** (`plot_dashboard.plqr_result`): layout
  `P-1 / (P-2 | NP-3 cumulative)` with the standard annotation block.

### 3.4 GMM (`gmm_result`)

- **G-1 — Coefficient estimates**: §3.0 estimate grammar with
  coefficient facets from `meta$coef_names` (works today via
  `plot_data_estimates.gmm_result`).  Overidentification visuals
  (local `J` profiles, `J_n(x)` maps) remain out of scope pending the
  Cond_GMM paper timeline.

### 3.5 Data requirements matrix

| Panel | needs `se` | needs `W` | needs `dx` | needs extras/args |
|---|---|---|---|---|
| NP-1 / C-1 / C-2 / P-1 / G-1 | optional | no | no | `X`,`Y` for NP-1 underlay |
| NP-2 kernel | no | yes | yes | — |
| NP-3 profile/cumulative | no | yes | optional | — |
| NP-4 diagnostics facets | no | via diagnostics | optional | — |
| NP-5 TDNN scales | optional | no | no | `X`,`Y` (refits) |
| NP-6 / tuning | optional | no | no | data for refits |
| C-3 mass / C-4 propensity | no | via diagnostics | no | stored `t_vec`/`e` or args |
| C-5 influence | no | yes | optional | `psi` extra or arg |
| P-2 conditioning | no | via diagnostics | no | — |
| P-3 calibration | no | yes | no | `g` extra or residuals |
| P-4 band mass | no | yes | no | residuals + `h` grid |

Every "needs" failure routes through the §2.2 precedence chain and ends
in an actionable error; no panel fabricates inputs.

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
- **Phase 1 — core panels.**  Panels NP-1, C-1, C-2, NP-3
  (profile + cumulative): `autoplot()`/`plot()` for `np_result` and
  `cate_result`, `plot_estimates()` (incl. forest style),
  `plot_weights()` (all three types; kernel view NP-2 included since it
  shares the weights-frame plumbing); Suggests updated; guard errors
  tested.  Gate: full suite green; README figure (C-1) generated.
  *Status 2026-07-03: complete.  `R/plot-common.R` (guard via mockable
  `.plot_pkg_installed()`, Okabe–Ito constants, theme, profile-axis
  detection with `x_var` override, honesty caption, mode selection,
  forest labels), `R/plot-estimates.R` (builder + `autoplot`/`plot`
  methods; delayed `S3method(ggplot2::autoplot, ...)` registration so
  ggplot2 stays in Suggests; `.data` pronoun via the tidy-eval mask with
  a `globalVariables` declaration), `R/plot-weights.R` (three `type`
  views, target subsetting, `log_weights`, TDNN zero line).  46
  structural tests in `test-plot-builders.R`; full suite 2422 passing;
  all CI guardrails green.  README gained the Visualisation section and
  the C-1 figure; regenerating README.md exposed a pre-existing
  source/generated drift (README.Rmd used `k = 2L` where the committed
  README.md had been hand-fixed to `k = 1L`) — fixed in the source.
  Phase 2 (overlap, diagnostics autoplot, dashboard) is next.*
- **Phase 2 — diagnostics and composition.**  Panels C-3, C-4, NP-4,
  C-6: `plot_overlap()` (both types), `autoplot.dnn_diagnostics()`,
  `plot_dashboard()` for cate/np with the patchwork degradation path and
  panel-availability substitution; dashboard annotation contract test.
  Gate: vignette draft renders.
- **Phase 3 — remaining classes, influence, tuning.**  Panels P-1, P-2,
  P-3, P-5, G-1, C-5, NP-5, NP-6: `autoplot.plqr_result` /
  `autoplot.gmm_result`, `plot_conditioning()`, `plot_calibration()`,
  `plot_influence()`, `plot_tdnn_scales()`, `plot_tuning()`.  P-4 lands
  in paper-replication material first (§10.7).  Gate: pkgdown section
  complete; vdiffr baselines recorded.
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

## 10. Decisions

Resolved during the 2026-07-02 detailing pass (§3):

1. **Sample-level inputs**: no training data is stored on fit objects for
   plotting; panels use the fixed precedence chain in §2.2 (stored model
   components → `return` extras → user arguments → actionable error).
2. **TDNN scale decomposition (NP-5)**: implemented by internal refits at
   `s1`/`s2`, not by new storage on the fit.
3. **Weight-panel scales**: linear y default; `log_weights = TRUE` opt-in
   with zero-dropping caption; distance-domain panels always linear.
4. **View multiplexing**: `type =` arguments on `plot_weights()` and
   `plot_overlap()` instead of one export per view.

Still open (resolve before or during the phase that ships them):

5. **`plot(fit)` default**: single estimates panel (proposed) vs. full
   dashboard.  Proposal: estimates panel; dashboards are explicit
   (`plot_dashboard()`), so casual `plot()` stays fast and legible.
6. **CI type**: plain Wald bands from stored `se` only, or also expose
   `level` recomputation from jackknife internals?  Proposal: stored `se`
   with `level` argument; nothing recomputed.  (Phase 0 implemented the
   proposal in `plot_data_estimates()`.)
7. **P-4 band-mass panel**: research diagnostic — ship in
   paper-replication material first (proposed) or in the public API.
8. **Default guide values**: `share_limits = c(0.1, 0.9)` for C-3 and
   `cond_limits = c(1e3, 1e6)` for P-2 are provisional; confirm against
   the paper's diagnostics prose before Phase 2/3 ship them.
9. **`shade_significant`** on C-1: keep default `FALSE` (proposed) so the
   flagship plot never invites naive significance reading; revisit only
   with a design that shades intervals, not point estimates.
10. **broom methods** (`tidy.cate_result` etc.): defer; `plot_data_*`
    covers the need without a new dependency surface.  Revisit if users
    ask.
11. **Naming**: `plot_dashboard()` vs `plot_summary()`; proposal
    `plot_dashboard()` to avoid clashing with `summary()` semantics.
12. Whether `plot_tuning()` belongs in Phase 3 or ships earlier as the
    public answer to the `s`-sensitivity question (it is the most
    seminar-valuable panel after the CATE profile).
