> **Status (2026-07-09): completed and archived.** The diagnostics layer
> landed 2026-04-26, and the remaining tuning-sensitivity follow-up this
> note carries was closed by the integrated tuning-grid work (design note
> 14; `tuning_mode = "grid"`, `candidate_fit()`, `plot_tuning()`).
> Code identifiers reflect pre-2026-07-08 API names (e.g. `kappa`, now
> `lambda`).

# Diagnostics and Reporting Implementation Plan

*Created: 2026-04-10. Updated: 2026-04-26 after the public diagnostics layer landed. Status: historical implementation plan plus the remaining tuning-sensitivity follow-up.*

## Purpose

This note turns `markdown_notes/RdistNN_implementation_note.md` into a concrete
implementation plan for the next package-facing work.  The immediate goal is to
make the existing estimator workflows easier to evaluate, report, and defend in
paper applications without prematurely expanding the package into a fully
general conditional Z-estimation solver.

The main conclusion is now:

1. The stable diagnostics/reporting layer for existing DNN wrappers has landed.
2. The main remaining simulation-facing package task is tuning-grid sensitivity on top of that diagnostics layer.
3. Compact diagnostics exist through `diagnostics()`, but summary/print integration remains partial.
4. The general conditional Z-estimation front end should still be deferred until the existing workflows are diagnosable and well documented.

## Status Update

As of 2026-04-26, the following pieces described in this plan are now live in
the repo:

- `diagnostics()` plus methods for `np_result`, `cate_result`, and `plqr_result`;
- universal weight/distance/tuning diagnostics;
- CATE overlap diagnostics;
- PLQR residual-Gram conditioning diagnostics;
- a `dnn_diagnostics` print method;
- README and vignette coverage, including `diagnostics-reporting` and
  `paper-cate-simulation`.

What remains live from this note is narrower:

- there is still no exported `tuning_sensitivity()` helper;
- estimator `print()` / `summary()` methods do not yet surface diagnostics as
  fully as they could; and
- the general conditional Z-estimation front end remains a later project.

Read the detailed phase material below accordingly:

- Phases 1-4 and the documentation pass are now historical route records.
- The tuning-sensitivity phase is still the main open simulation-facing follow-up.

## Design-Note Alignment

This plan follows the repository design notes in `markdown_notes/design/`.

### Layering

The diagnostics feature should preserve the existing package layering:

- Core layer: low-level diagnostics helpers that operate on canonical numeric
  matrices and vectors, especially `W`, optional `dx`, optional ranks/orders,
  and metadata.
- User layer: S3 diagnostics methods for `np_result`, `cate_result`, and
  `plqr_result`.
- Documentation / replication layer: examples, vignettes, and simulation
  reporting that consume the public diagnostics API rather than reaching into
  internal fields.

This avoids mixing estimator computation, reporting, and paper-specific
replication code in the same helper.

### Interface Discipline

No new user-facing estimator input convention should be introduced in the first
diagnostics pass.  The public wrappers already resolve matrix and data-frame
inputs to canonical internal forms.  Diagnostics should consume fitted objects
or canonical internal objects, not reimplement data-frame parsing.

If later tuning sensitivity needs to re-run an estimator, the replay contract
should be explicit.  Do not hide refitting behind fragile assumptions about
which fields happened to be retained in an object.

### Validation And Messages

All new exported behavior should validate early and use the existing message
style:

- no `stopifnot()` on user-facing paths;
- clear argument names in backticks;
- `call. = FALSE` for user-facing errors and warnings;
- separate missing, shape, non-finite, and semantic failures.

Diagnostics methods should fail clearly when a fit does not contain enough
information.  For example, if `W` is absent and the object does not have
replay-safe model state, the method should say that diagnostics require
retained weights or a refit with `return = "W"`.

### Return-Object Stability

Diagnostics should be additive.  They should not change the meaning of existing
fit objects or estimator return values.  A successful diagnostics method should
return a stable classed object or data frame with predictable columns.

The first implementation should not make diagnostics mandatory fields on all
fit objects.  Instead, expose diagnostics through methods and optionally cache a
compact diagnostics object later if there is a clear performance reason.

### Namespace Discipline

Export only the stable public surface:

- a `diagnostics()` S3 generic;
- methods for package fit classes;
- later, a tuning sensitivity helper if its user story is stable.

Keep low-level helpers internal.  Do not hand-edit `NAMESPACE`; use roxygen and
regenerate generated files in the implementation pass.

### Testing Standard

Every new user-visible behavior should have tests.  The first pass needs tests
for numeric identities, return shapes, unavailable-diagnostic behavior, summary
integration, and tuning-grid sensitivity.

## Current Implementation Facts

The relevant source files are:

- `R/user-dnn-core.R`: shared DNN helpers, common extras, jackknife support.
- `R/user-dnn-npr.R`: `np_dnn()` and `np_result` construction.
- `R/user-dnn-cate.R`: `cate_dnn()`, CATE result construction, hidden
  prediction model state.
- `R/user-dnn-cate-methods.R`: CATE summary, confidence intervals, prediction.
- `R/user-dnn-plqr.R`: `plqr_dnn()` and PLQR moment construction.
- `R/user-dnn-np-methods.R`: `print.np_result()`.
- `R/user-dnn-plqr-methods.R`: `print.plqr_result()`.
- `R/api-tuning.R`: tuning object API and grid/default evaluation.

Important object-state constraints:

- `np_dnn()` still does not retain replay-safe model state for arbitrary refits.
- `cate_dnn()` now retains hidden `.model` state that includes `x_mat`, `psi`,
  `t_vec`, `e`, distance settings, and tuning parameters.
- `plqr_dnn()` now retains hidden `.model$d_resid_mat` for Gram diagnostics,
  but it still computes localized moments at a supplied `beta`; it does not
  estimate a local root or a Jacobian/bread object.
- All three wrappers usually return `W` by default, but users can request other
  `return` values and may create objects without distance or rank extras.

These facts now shape the **remaining** tuning-sensitivity helper.  Do not
design that helper as if every fit object can be replayed.

## Core Diagnostic Concepts

### Universal Weight Diagnostics

For each target point, report:

- `target_id`
- `n`, `p`, `k`, `s`, `method`, `distance`
- `n_eff = 1 / sum_i w_i^2`
- `max_weight`
- `top1_weight_share`
- `top5_weight_share`
- `top10_weight_share`
- `n_nonzero_weight`
- `n_above_1pct`, optional
- `n_above_5pct`, optional

The top-share calculations should use sorted absolute or signed weights only
after an explicit design decision.  For ordinary DNN weights, weights are
nonnegative and columns sum to one.  Because TDNN weights can be signed, the
first DNN diagnostics pass should either be DNN-only or clearly define signed
weight behavior before extending to TDNN.

### Distance And Neighborhood Diagnostics

If distance information is available, report:

- `radius_rank_k`: distance to the `k`-th nearest observation;
- `radius_rank_s`: distance to the `s`-th nearest observation;
- `weighted_mean_distance`;
- `weighted_median_distance`;
- optional weighted distance quantiles.

Terminology matters.  Exact DNN weights are not an `s`-nearest-neighbor
truncation.  Therefore `radius_rank_s` should be presented as a realized
rank-scale diagnostic, not as the support radius of the DNN estimator.

If `dx` is absent, distance fields should be `NA_real_` and a diagnostics
metadata flag should record that distance diagnostics were unavailable.

### Tuning Diagnostics

If tuning metadata is available through `object$meta$tuning` or an explicitly
returned `object$tuning`, report:

- `tuning_source`;
- `tuning_name`;
- `tuning_default_s`;
- `tuning_grid_available`;
- `tuning_grid_n_candidates`;
- optional default `kappa`, when present.

Diagnostics should not duplicate the full tuning grid in every target row.
Store the full grid as an attribute or a separate component on the diagnostics
object if needed.

## CATE-Specific Diagnostics

CATE diagnostics should focus on localized overlap.

For each target, report:

- `overlap_available`;
- `treated_mass = sum_i w_i D_i`;
- `control_mass = sum_i w_i (1 - D_i)`;
- `treated_share = treated_mass / (treated_mass + control_mass)`.

If propensity scores are available, also report:

- `propensity_available`;
- `weighted_mean_propensity`;
- `weighted_propensity_q10`;
- `weighted_propensity_q90`.

Implementation implication:

- When `cate_dnn()` constructs pseudo-outcomes from nuisance columns, retain
  `t_vec` and `e` in hidden model state for diagnostics.
- When `cate_dnn()` is called with direct `pseudo_outcome`, overlap diagnostics
  should be unavailable unless the API later adds an explicit treatment vector.

Do not infer treatment status from the pseudo-outcome.

## PLQR-Specific Diagnostics

The first PLQR diagnostics pass should not claim to estimate the bread unless a
true local Jacobian/bread is computed.

Current `plqr_dnn()` computes the localized moment

```text
sum_i w_i(x0) * g_i(beta)
```

at a user-supplied `beta`.  A useful first identification diagnostic is the
weighted residual-regressor Gram matrix

```text
G(x0) = sum_i w_i(x0) * (D_i - m_D(X_i)) (D_i - m_D(X_i))'
```

For each target, report:

- `gram_available`;
- `gram_min_eig`;
- `gram_max_eig`;
- `gram_condition_number`;
- `gram_rank_tol`, if a rank tolerance is used.

Implementation implication:

- Retain `d_mat` and `mD_mat`, or at least `d_resid_mat = d_mat - mD_mat`, in
  hidden model state for PLQR diagnostics.
- Keep this model state internal.  Do not add it to the public return payload
  unless the user explicitly requests a documented optional field in a later
  pass.

Later, a general Z-estimation front end can add finite-difference or analytic
bread diagnostics using the score API.

## Proposed API

### Public Generic

Add:

```r
diagnostics <- function(object, ...) UseMethod("diagnostics")
```

Methods:

```r
diagnostics.np_result <- function(object, ...)
diagnostics.cate_result <- function(object, ...)
diagnostics.plqr_result <- function(object, ...)
```

Return a classed data-frame-like object, for example:

```r
class(out) <- c("dnn_diagnostics", "data.frame")
```

This keeps the first public surface small and R-native.

### Internal Helpers

Add internal helpers such as:

```r
.dnn_weight_diagnostics()
.dnn_tuning_diagnostics()
.dnn_weighted_quantile()
.cate_overlap_diagnostics()
.plqr_gram_diagnostics()
```

These should operate on canonical numeric matrices and vectors only.

### Print Method

Add:

```r
print.dnn_diagnostics <- function(x, ...)
```

The print method should be compact:

- number of targets;
- `n_eff` range;
- maximum top-5 share;
- whether distance diagnostics are available;
- object-specific headline diagnostics if present.

## Implementation Phases

### Phase 1: Generic Weight Diagnostics

Add the internal weight diagnostics helper and the `diagnostics()` generic.

Minimum deliverable:

- `diagnostics(np_result)` works when `W` is present.
- `diagnostics(cate_result)` and `diagnostics(plqr_result)` return universal
  weight diagnostics when `W` is present.
- Distance fields are populated only when `dx` is present.
- Tuning metadata is summarized when available.

Tests:

- hand-built `W` gives exact `n_eff`;
- top-share calculations are correct;
- multi-target shape is correct;
- missing `dx` produces `NA` distance diagnostics without error;
- missing `W` errors clearly.

### Phase 2: CATE Overlap Diagnostics

Extend `cate_dnn()` hidden `.model` state with treatment and propensity fields
when they are available.

Minimum deliverable:

- `diagnostics(cate_result)` adds overlap diagnostics for nuisance-column
  calls.
- pseudo-outcome-only calls return `overlap_available = FALSE`.
- no public CATE estimator arguments change.

Tests:

- known `D` and `W` produce exact treated/control masses;
- pseudo-outcome-only path has unavailable overlap diagnostics;
- data-frame and matrix/pseudo-outcome paths remain backward compatible.

### Phase 3: PLQR Gram Diagnostics

Extend `plqr_dnn()` hidden model state with `d_resid_mat`.

Minimum deliverable:

- `diagnostics(plqr_result)` adds Gram conditioning diagnostics when residual
  regressors are available.
- If old objects lack model state, diagnostics still return universal weight
  diagnostics when possible and mark Gram diagnostics unavailable.

Tests:

- simple residual-regressor design gives known Gram matrix eigenvalues;
- singular local design produces infinite or very large condition number in a
  documented, stable way;
- existing PLQR return extras and print method are not broken.

### Phase 4: Summary And Print Integration

Update print methods to show compact diagnostics when available from retained
fields.  Do not make printing expensive.

Recommended behavior:

- If `W` is present, print `n_eff` range and max top-5 share.
- If `W` is absent, leave current print behavior unchanged.
- For CATE, show local treated/control mass only when overlap diagnostics are
  available.
- For PLQR, show Gram condition number only when Gram diagnostics are available.

Update `summary.cate_result()` to include a `diagnostics` component.  Decide
later whether `np_result` and `plqr_result` need exported summary methods.

Tests:

- `print()` still runs for all three fit classes;
- `summary.cate_result()` keeps existing coefficient columns;
- diagnostics component appears only when enough information is present.

### Phase 5: Tuning-Grid Sensitivity

Build tuning sensitivity after diagnostics are stable.
As of 2026-04-26, this is the first phase in this note that is still genuinely open.

Do not start with `sensitivity(fit)` as the only API because current fit
objects are not all replay-safe.  Use an explicit estimator-call API first:

```r
tuning_sensitivity(
  estimator,
  args,
  tuning = NULL,
  level = 0.95,
  se = TRUE
)
```

where:

- `estimator` is `np_dnn`, `cate_dnn`, or `plqr_dnn`;
- `args` is a named list of estimator arguments excluding `s` and `tuning`;
- `tuning` supplies the grid, or the helper uses tuning metadata from `args`.

Minimum output fields:

- `candidate_id`;
- `kappa`, if present;
- `s`;
- `feasible`;
- `reason`;
- target/component identifiers;
- estimate or moment;
- standard error, if available;
- confidence interval length, where scalar inference is available;
- core diagnostics;
- CATE overlap or PLQR Gram diagnostics when available.

Tests:

- one row per candidate-target for NP and CATE;
- PLQR output has component and target identifiers;
- infeasible candidates are retained with `NA` estimates and a reason;
- helper rejects missing/invalid `args` with clear messages.

### Phase 6: Documentation And Site Integration

Update:

- `README.Rmd`;
- `vignettes/getting-started.Rmd`;
- `vignettes/conditional-z-estimation.Rmd`;
- `vignettes/kernel-order-tuning.Rmd`;
- `_pkgdown.yml`;
- `NEWS.md`.

Docs should show the preferred public API:

```r
fit <- cate_dnn(...)
diagnostics(fit)
```

and, after Phase 5:

```r
tuning_sensitivity(cate_dnn, args = list(...), tuning = tuning_dnn_paper())
```

Do not make diagnostics examples depend on optional packages.

## Non-Goals For The First Diagnostics Project

Do not build a full general conditional Z-estimation front end in this pass.

Do not add automatic nuisance fitting in this pass.

Do not expose every internal helper as public API.

Do not call PLQR Gram diagnostics a bread estimate unless the implementation
actually estimates the local Jacobian.

Do not rely on raw return extras as the only diagnostics interface once the
public `diagnostics()` method exists.

## Longer-Run Framework Task

After diagnostics and tuning sensitivity are stable, build a general localized
conditional Z-estimation front end on top of the score API.  That future task
should be treated as a separate design project.

A future front end should likely:

- accept a score or moment specification;
- accept nuisance inputs and derivative information;
- localize with DNN weights;
- solve the local estimating equation;
- compute or approximate bread/Jacobian information;
- return estimate, standard error, and diagnostics through the same diagnostics
  interface designed here.

## Practical Implementation Order

Recommended commit-sized sequence:

1. Land generic weight diagnostics, `diagnostics()` generic, and
   `diagnostics.np_result()`. Landed.
2. Extend generic diagnostics to CATE and PLQR universal fields. Landed.
3. Add CATE overlap state and diagnostics. Landed.
4. Add PLQR residual-regressor state and Gram diagnostics. Landed.
5. Integrate compact diagnostics into print/summary methods. Partially landed.
6. Add explicit-call tuning sensitivity helper. Still open.
7. Update docs, pkgdown reference groups, and NEWS. Largely landed.
8. Regenerate roxygen outputs and run tests. Historical implementation step.

Verification commands for the final implementation pass:

```powershell
Rscript -e "devtools::document()"
Rscript -e "devtools::test()"
Rscript -e "devtools::check(document = FALSE)"
```

If local dependency setup makes `devtools::check()` impractical, at minimum run
the relevant `testthat` files plus any repository CI guardrail scripts that
check generated-file drift and pkgdown reference coverage.
