# API shapes and input formats

The [interface decisions](01_interface_standardization_decisions.md) define
the dual public interface and canonical internal boundary.

## Inputs

| Input | Shape/meaning |
| --- | --- |
| `X` | n by p training covariates |
| `X0` | m by p reporting targets |
| `x0` | Length-p single target; exclusive with `X0` |
| Outcome/signal/nuisance vectors | Length n, in training-row order |
| `data`, `x_cols` | Data frame (including tibbles) and explicit covariate selectors |
| `y_col`, `t_col`, `m0_col`, `m1_col`, `e_col` | Explicit outcome, treatment and nuisance column selectors |

Use documented matrix-like inputs where supported. Data routes use explicit
selectors rather than NSE or guessed columns, validate extracted values and
convert them to numeric matrices/vectors. Internal single targets may use a
1 by p matrix. Batch documentation must identify which axes index observations,
targets, signal components or dimensions.

Core NP/TDNN/CATE boundaries reject ignored direct inputs and selectors without
`data`. NP/TDNN allow `data` plus direct `y` instead of `y_col`, but not
matrix-route `X` or `Y`. CATE allows `data` plus `pseudo_outcome`, but not `X`.
`.dnn_validate_input_route()` checks conflicts before canonicalization. Paired
ratios/CLATE retain strict column/direct signal exclusivity.

PLQR uses either `data` with selectors or direct `X`, `Y`, `D`, `g` and `mD`.
Its direct `density` vector is allowed with either route; `density_col`
requires `data`, and those density inputs are mutually exclusive.

Shared batch evaluators preserve `n x m` distance/rank-weight matrices when
either dimension is one. Custom batch callbacks must return actual matrices;
scalar adapters assemble validated columns. Distance setup runs once per batch
when context is absent. Scalar scores return one finite numeric value per
subset. Analytic batch gradients have exactly two dimensions, `K x p`; a
length-`K` vector is accepted only for `p = 1` and canonicalized to a matrix.

## Results

Successful high-level calls return classed objects with stable estimates,
metadata and documented methods, rather than unnamed lists. Optional `return`
extras are additive and unknown entries warn and are ignored under the current
interface. Main fields must not change meaning when an extra is requested.

Use the common argument names `X`, `X0`, `x0`, `x_cols`, `k`, `s`, `method`,
`distance`, `jackknife` and `return`; introduce alternatives only for an
explicitly documented current-interface requirement. New estimators follow the same sequence: resolve/validate,
canonicalize, compute weights/signals, estimate, calculate requested uncertainty,
and assemble the result. See [architecture](13_estimator_and_result_architecture.md)
for class-specific fields and shared methods.
