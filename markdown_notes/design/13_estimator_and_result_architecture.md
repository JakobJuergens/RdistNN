# Estimator and result architecture

This document describes the shared implementation and presentation contracts.
Read the [interface policy](01_interface_standardization_decisions.md),
[shape contract](04_api_shapes_and_input_formats.md) and
[validation policy](03_input_validation_contract.md) for boundary rules.

## Computation layers

High-level estimators validate options and the selected input interface, resolve
columns and targets, and convert inputs to numeric matrices and vectors. Core
computation then builds weights, evaluates signals or moments, calculates the
estimate and uncertainty, and assembles a classed result. Optional `return`
fields are additive. Lower-level primitives keep explicit numeric interfaces.

Reuse the distance/ranking, rank-weight, tuning, validation and jackknife layers.
Compiled routines implement the same numerical contract as their R fallbacks;
they do not define another API. Keep user-facing validation at the R boundary.

Jackknife covariance paths subtract a common numerical reference per target
and signal before summing deletions. Choose that reference only from the union
of observations with nonzero weight in any deletion. A nonzero entry at rank
`r` in the `n-1` table can use original ranks `r` and `r+1`; support must include
both. Ignore unused observations before subtraction, including extreme values
that could overflow despite their zero weights. Signed TDNN tables use nonzero
support. Asymptotic tables retain all represented nonzero weights without an
extra cutoff or normalization. Covariance calculations keep the shifted values;
requested raw deletion estimates restore the common weighted offset. Explicit
deletion algorithms remain available and use the same centering rule. When a
reference path recomputes geometry, determine its union from the actual deletion
weights rather than imposing frozen full-sample ranks.

## Distinct result classes and common interfaces

| Result | Main quantity | Uncertainty and interpretation |
| --- | --- | --- |
| `np_result` (including TDNN) | `mu` | Conditional mean, `se`, joint target covariance |
| `cate_result` | `tau` | Conditional effect, `se` and on-demand joint target covariance |
| `affine_ratio_result` / `clate_result` | `theta`, paired `components` | Component covariance, ratio variance and per-target computational status |
| `plqr_result` | `moment`, supplied `theta_eval` | Moment diagnostics; coefficient inference requires the documented sandwich |
| `gmm_result` | Solved `theta` | GMM covariance, convergence and identification diagnostics |

CLATE inherits from `affine_ratio_result`, not `cate_result`. A shared display
does not make distinct estimands or covariance calculations interchangeable.
PLQR evaluated moments must remain distinguishable from solved coefficients.

Existing generics select the relevant adapter from the result class:
`summary()`, `coef()`, `vcov()`, `confint()`, `predict()`, `diagnostics()`,
`plot_data_estimates()`, `candidate_fit()` and `candidate_fits()`, where each
class documents support. Builders consume common tables rather than guessing
the meaning of similarly named raw fields. Do not change a fit's class merely
to reuse another estimator's methods.

CATE and ratio summaries use a coefficient table with `estimate`, `std.error`,
`statistic`, `p.value`, `conf.low` and `conf.high`. Ratios additionally retain
component uncertainty and target-status tables. Undefined or zero-score
variance targets can have unavailable Wald bounds even when covariance was
requested; retain these rows and report their status.

## Shared prediction boundary and retained state

NP, CATE and paired-ratio/CLATE methods validate their own fitted scientific
state, then resolve `newdata`, `X0` or `x0` through
`.dnn_resolve_prediction_targets()` in `R/user-dnn-prediction.R`. The shared
boundary owns interface exclusivity, tabular column extraction and finite
canonical target matrices. Explicit target resolution delegates to the same
`.dnn_resolve_targets()` used by estimator wrappers. Tabular models retain
resolved covariate names, so numeric fit-time selectors do not become positional
selectors in reordered prediction data.

Retained state follows the same structural grammar: canonical `x_mat`, resolved
`x_cols`, weight controls (`method`, `k`, `s`, `normalize`) and geometry
(`distance`, `ctx`, cached ranks where available). Signal fields remain
explicit about their meaning: NP uses `y_vec`, CATE uses `psi`, and paired
ratios use a two-column `signal`. TDNN adds `s1`, `s2` and `dim`. Common
candidate adapters update the retained order; estimator adapters own the
resulting means, effects or paired covariance/status quantities. Sharing the
target boundary does not change these inference contracts.

`R/user-dnn-geometry.R` owns `.dnn_fitted_geometry()` for NP/TDNN, CATE and
paired-ratio/CLATE fits. It resolves the specification and training context once,
retains a setup-free specification for prediction, and separately retains the
original resolved specification/context for legacy deletion refitting. Registry
replacement therefore cannot change a saved fit. NP/TDNN and CATE `loo` keep
their established deletion-specific geometry; paired-ratio `loo` remains the
frozen-rank reference. Custom closures with mutable external state retain their
own reproducibility obligations.

## Nonparametric regression reporting and prediction

CATE also exposes `vcov()` from retained frozen signals, without requiring
optional weights. Common deletion rows supply cross-target covariance;
`fast` freezes fitted geometry while `loo` uses the retained original metric
and inference context for deletion-specific refitting. A missing original
context field differs from an explicitly retained `NULL`. Covariance requires
an active fitted jackknife, preserves zero/singular matrices and follows
coefficient target names. This calculation does not establish joint inference
assumptions. Printing reports order fractions descriptively and causal
diagnostics summarize existing quantities without warning thresholds or schema changes.

DNN and TDNN regression share `np_result` methods for `summary()`, `coef()`,
`confint()`, `vcov()` and `predict()`. The adapter retains canonical training
outcomes/covariates, resolved weight controls, column names and fitted distance
specification/context independently of optional extras. Resolving a registered
metric at fitting time protects saved fits from subsequent registry changes.
Prediction passes the full retained training matrix and fitted context to the
metric, never the new targets as a training sample. Custom metrics own the
contract to respect their supplied context; arbitrary callable side effects
cannot be frozen by the package.

`vcov()` constructs the full joint estimator covariance across targets from
common deletion rows, preserving off-diagonal dependence. It calculates this
on request instead of retaining an additional quadratic target-size matrix in
every fit. Candidate views update the retained order, so covariance and
prediction reconstruct the selected candidate. The fitted covariance matches
the stored SE calculation: `fast` uses frozen ranks, while the legacy `loo`
reference recomputes deletion geometry using the original context. Prediction
uncertainty always uses frozen training geometry and fast rank deletions.

NP summaries follow the CATE coefficient-table grammar, labeling local means.
Absent uncertainty stays missing; covariance/interval extraction requires an
active fitted jackknife. Prediction can independently request uncertainty from
a fit with `se = FALSE` when deletion orders are feasible. Exact zero variance
is retained as zero SE/covariance. The shared computational Wald grammar
reports collapsed bounds; nonzero means give infinite statistics and zero
p-values, while zero-over-zero statistics and p-values are undefined. These
outputs do not establish negligible bias or valid coverage. This is pointwise
reporting and does not supply simultaneous bands.

## Implementation ownership

| Responsibility | Shared module or helper | Estimator adapter owns |
| --- | --- | --- |
| Covariance arithmetic | `R/utils-covariance.R` | Statistical normalization, result classes and per-target failure statuses |
| Fit input routes and canonical matrices | `R/user-dnn-inputs.R` | Allowed signal/selector combinations and estimator-specific validation |
| Mahalanobis geometry | `R/dist_mahal.R` | Retained fitted context; existing L2 engines evaluate projected norms |
| Prediction target interface | `.dnn_resolve_prediction_targets()` in `R/user-dnn-prediction.R` | Fitted-state validation and method-specific signal/uncertainty |
| NP retained state and joint covariance | `R/user-dnn-np-core.R` | Local means and DNN/TDNN controls; public methods in `R/user-dnn-np-methods.R` |
| DNN rank weights | `.dnn_rank_weight_spec()` in `R/user-dnn-core.R` | Exact/asymptotic rule, `k`, `s`, normalization |
| Order and candidate identities | `.dnn_resolve_s()` in `R/user-dnn-tuning.R` | Statistical feasibility and tuning specification |
| Grid evaluation and cache reuse | `.dnn_grid_eval()` in `R/user-dnn-grid.R` | Estimate/SE callbacks and target-specific statuses |
| Shared deletion orchestration | `R/user-dnn-jackknife.R` | Signal meaning, geometry refitting and statistical normalization |
| Paired ratios | `.affine_ratio_evaluate()` | CLATE signal construction and IV metadata |
| Candidate views | `R/candidate-fit.R` | Every candidate-dependent estimate, covariance and optional extra |
| Tidy estimate frames | `R/plot-data.R`, `R/plot-common.R` | Estimate scale, coefficient names and unavailable bounds |
| Estimate and weight panels | `R/plot-estimates.R`, `R/plot-weights.R` | Tidy extractor and meaningful reference lines |
| Overlap and influence | `R/plot-overlap.R`, `R/plot-influence.R` | Treatment versus instrument arms; single versus paired score residuals |
| Tuning and dashboards | `R/plot-tuning.R`, `R/plot-dashboard.R` | Frozen-state grid adapter and appropriate diagnostics |

The [tuning contract](14_integrated_tuning_grid_reuse.md) specifies grid
schemas and reuse. The [ratio contract](17_clate_ratio_inference.md) specifies
paired uncertainty. The [PLQR rationale](15_plqr_variance_estimation_rationale.md)
explains why moment uncertainty does not suffice for coefficient inference.

## CATE and CLATE presentation

Both methods use `plot_estimates()`, `plot()` and `ggplot2::autoplot()` with
profile, forest and pointrange styles; `plot_weights()` with profile,
cumulative and kernel views; and the common overlap, influence, tuning and
dashboard calls. Weight panels require the documented retained `W`/`dx`
extras. Optional ggplot2 and patchwork behavior follows the same policy.

CATE overlap uses treatment arms and treatment propensity. CLATE uses
instrument arms and instrument propensity. CATE influence panels show
weighted pseudo-outcome residuals; CLATE uses
`w * (phi_y - theta * phi_d)`. These score contributions are not individual
causal effects or deletion pseudo-values. Ratio plots keep finite estimates
visible when their interval is unavailable, and explain unavailable output
in captions. The tidy tables retain all target rows.

Exact and asymptotic single-scale DNN are the CATE/CLATE weighting choices.
Two-scale TDNN remains NP-regression functionality. Supporting a computational
choice does not establish its statistical assumptions or extend a theorem.

## Extending an estimator

1. Define the estimate, covariance scale, class and retained numerical state.
   Preserve the shared input/validation grammar and additive optional outputs.
2. Implement the scientifically meaningful existing generic methods and tidy
   extractors. Reuse validators and frame builders before adding a new API.
3. Supply grid callbacks using frozen signals and geometry. Candidate views
   must rebuild every dependent field, including covariance and deletion
   arrays, without training nuisances or recomputing a retained ranker.
4. Reuse plot builders, keeping score, uncertainty and label differences in
   explicit adapters. Preserve candidate selection and dependency guards.
5. Test independent analytic references, special-case reductions, supported
   weighting controls, candidate/scalar equality, mixed target statuses,
   tabular inputs, serialization, prediction and R/compiled parity. Check
   actual plot rendering as well as tidy data.

See [the test-suite map](../../tests/README.md) for focused runs and
[CONTRIBUTING.md](../../CONTRIBUTING.md) for documentation generation and checks.
