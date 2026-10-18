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

## Distinct result classes and common interfaces

| Result | Main quantity | Uncertainty and interpretation |
| --- | --- | --- |
| `np_result` | `mu` | Conditional mean and its `se` |
| `cate_result` | `tau` | Conditional effect and its `se` |
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

## Implementation ownership

| Responsibility | Shared module or helper | Estimator adapter owns |
| --- | --- | --- |
| DNN rank weights | `.dnn_rank_weight_spec()` in `R/user-dnn-core.R` | Exact/asymptotic rule, `k`, `s`, normalization |
| Order and grid evaluation | `.dnn_resolve_s()`, `.dnn_grid_eval()` | Estimate/SE callbacks and target-specific statuses |
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
