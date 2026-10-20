# API Shapes and Input Formats

This note standardizes how inputs and outputs should be structured across
`RdistNN`.

## 1. Dual-interface policy

Where appropriate, user-facing estimators may support both:

- matrix interface; and
- `data.frame` + column-name interface.

This is a feature, not a burden, but only if both interfaces remain predictable.

## 2. Matrix interface

The matrix interface is the lower-level, more explicit contract.

Standard expectations:

- `X` is an `n x p` matrix-like object;
- `X0` is an `m x p` matrix-like object when multiple targets are allowed;
- `x0` is a length-`p` vector for a single target;
- signal or pseudo-outcome vectors are length `n`.

Whenever both `X0` and `x0` are supported:

- exactly one should be supplied;
- the resolver should convert them into one internal canonical target format.

## 3. Data-frame interface

The data-frame interface is the more convenient user path.

Standard expectations:

- `data` must be a `data.frame`;
- tibbles are allowed as a subtype of `data.frame`;
- `x_cols` identifies the covariate columns;
- signal / outcome / treatment / nuisance columns are passed explicitly by name;
- extracted columns are validated immediately after resolution.

The data-frame interface should never rely on NSE or implicit column guessing.
Core NP/TDNN/CATE boundaries reject ignored direct inputs and selectors without
`data`. NP/TDNN allow `data` plus a direct `y` vector in place of `y_col`, but
not matrix-route `X` or `Y`. CATE allows `data` plus `pseudo_outcome`, but not
`X`. `.dnn_validate_input_route()` checks these route conflicts before numeric
canonicalization. Paired ratio and CLATE retain their strict column/direct
signal exclusivity.

## 4. Canonical internal forms

After validation, code should move to canonical internal forms.

Typical canonical forms:

- covariates as numeric matrix;
- single-target input as `1 x p` matrix where convenient;
- signals / pseudo-outcomes as numeric vector;
- optional metadata collected in named lists.

This reduces branching deeper in the code.

## 5. Return-value policy

Returned objects should be predictable and stable.

### Core rule
Requested optional returns should be additive.
They should not change the meaning of the main result.

### Standard behavior
- main fields should always exist when the call succeeds;
- optional fields should only be added when requested;
- unknown requested return fields should either:
  - warn and be ignored; or
  - error, if strict mode is later introduced.

If warnings are used, the message should list valid options.

## 6. Object classes

High-level estimators should return classed objects with:

- a stable core payload;
- clearly named metadata;
- methods for `print()`, `summary()`, `confint()`, `predict()` where relevant.

Do not return ad hoc unnamed lists for user-facing estimation results.

## 7. Naming conventions for arguments

Prefer consistent names across the package.

### Covariates / targets
- `X`, `X0`, `x0`
- `x_cols`

### Outcomes / treatment / nuisances
- `y_col`, `t_col`, `m0_col`, `m1_col`, `e_col`
- analogous matrix/vector names in lower-level APIs

### DNN controls
- `k`, `s`, `method`, `distance`, `jackknife`

### Optional return control
- `return`

Avoid introducing synonyms unless there is a compelling compatibility reason.

## 8. Shape conventions for batch operations

Batch APIs should use the convention:

- rows index observations / training sample;
- columns index targets or dimensions, depending on the context;
- documentation must say which.

Any function returning batched distances, ranks, or weights must document the
expected shape explicitly.

## 9. Interface discipline for future estimators

A new user-level estimator should follow the same overall grammar as the
existing ones:

- validate and resolve interface;
- convert to canonical matrix/vector form;
- compute weights or score contributions;
- compute estimate;
- optionally compute jackknife SE;
- assemble classed result.

This consistency is more important than local cleverness.
