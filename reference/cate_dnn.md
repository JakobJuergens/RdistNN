# CATE via DNN rank weights (orthogonalized local averaging)

Estimate conditional average treatment effects by combining a
pseudo-outcome with DNN-based local weights.

## Usage

``` r
cate_dnn(
  data = NULL,
  x_cols = NULL,
  y_col = NULL,
  t_col = NULL,
  m0_col = NULL,
  m1_col = NULL,
  m_col = NULL,
  e_col = NULL,
  pseudo_outcome = NULL,
  X = NULL,
  X0 = NULL,
  x0 = NULL,
  clip = 1e-08,
  distance = "l2",
  method = c("exact", "asymp"),
  k = 1L,
  s = NULL,
  tuning = NULL,
  tuning_mode = c("default", "grid"),
  normalize = TRUE,
  se = TRUE,
  jackknife = c("fast", "loo", "none"),
  ctx = NULL,
  return = c("W")
)
```

## Arguments

- data:

  Data frame with covariates and optional nuisance columns. Column names
  must be non-missing, non-empty, and unique.

- x_cols:

  Names or indices of covariate columns in `data`.

- y_col:

  Name or index of the outcome column in `data`.

- t_col:

  Name or index of the treatment column in `data` (must be 0/1).

- m0_col:

  Name or index for \\E\[Y(0)\mid X\]\\ (full AIPW interface).

- m1_col:

  Name or index for \\E\[Y(1)\mid X\]\\ (full AIPW interface).

- m_col:

  Name or index for \\E\[Y\mid X\]\\ (Robinson interface).

- e_col:

  Name or index for propensity \\P(T=1\mid X)\\.

- pseudo_outcome:

  Numeric vector of length `n`. If supplied, nuisance columns are not
  required.

- X:

  Numeric matrix (`n x p`) for the matrix interface.

- X0:

  Numeric matrix (`m x p`) of target points; mutually exclusive with
  `x0`.

- x0:

  Numeric vector of length `p` for a single target point.

- clip:

  Finite propensity clipping constant in `(0, 0.5)`, default `1e-8`. Raw
  nuisance propensities must be finite and lie in `[0, 1]`;
  nuisance-based AIPW and Robinson scores use
  `pmin(pmax(e, clip), 1 - clip)`. This argument does not alter a
  supplied `pseudo_outcome`.

- distance:

  Distance name or specification (for example `"l2"`).

- method:

  One of `"exact"` or `"asymp"` for DNN weights.

- k:

  Integer DNN neighbour count inside each size-`s` subsample. The paper
  theory uses `k = 1` (the default).

- s:

  Integer kernel order (undersmoothing parameter), or `NULL`. When both
  `s` and `tuning` are `NULL` (the default), `s` is resolved from the
  paper undersmoothing rule
  [`tuning_dnn_paper()`](https://JakobJuergens.github.io/RdistNN/reference/tuning_dnn_paper.md),
  \\s = \lfloor 0.4 n^{(k_x+1)/(k_x+2)} \rfloor\\. Note that small fixed
  values such as `s = 1` (or `k = s`) yield uniform weights and hence no
  localization; a warning is emitted in that case.

- tuning:

  Optional kernel-order tuning specification. When supplied, `s` must be
  omitted and is resolved from `tuning`.

- tuning_mode:

  One of `"default"` or `"grid"`. `"default"` evaluates the default `s`
  only. `"grid"` evaluates the feasible candidates from `tuning`, or
  from
  [`tuning_dnn_paper()`](https://JakobJuergens.github.io/RdistNN/reference/tuning_dnn_paper.md)
  when `tuning` is omitted, and stores candidate-level estimates,
  standard errors, diagnostics, and metadata in `tuning_grid`. Grid mode
  cannot be combined with explicit `s`.

- normalize:

  Logical; renormalize columns of `W` to sum to 1.

- se:

  Logical; compute delete-1 jackknife standard errors.

- jackknife:

  One of `"fast"`, `"loo"`, or `"none"`. `"fast"` (default) uses cached
  full-sample ranks with delete-1 rank updates, `"loo"` recomputes
  leave-one-out weights directly, and `"none"` skips SEs.

- ctx:

  Optional context passed to distance and ranking components.

- return:

  Optional extras in addition to the always-returned `tau`, `se`, and
  `meta`. Allowed values:
  `c("W", "ranks", "ord", "dx", "psi", "nuisance", "ess", "neighbors", "tuning")`.
  The `neighbors` extra returns the first `k` nearest observations used
  inside each size-`s` subsample, not the full DNN weight support.

## Value

An object of class `cate_result` with components:

- `tau`: estimated CATE values at target points;

- `se`: optional jackknife standard errors;

- `x0`: the `m x p` matrix of evaluation points, with covariate column
  names propagated from the training data when available;

- `meta`: metadata (`n`, `m`, `p`, `k`, `s`, `method`, `distance`);

- `tuning_grid`: in grid mode, a list with `candidates`, `estimates`,
  `diagnostics`, and `meta` tables. Candidate diagnostics report top-5
  weight concentration separately from the estimator's `k`;

- optional extras requested via `return`.

## Details

Two interfaces are supported.

1.  Data-frame interface (recommended): pass `data` and `x_cols`, then
    provide either `pseudo_outcome` directly or nuisance columns used to
    construct it. For nuisance columns, use either:

    - full AIPW form: `y_col`, `t_col`, `e_col`, `m0_col`, `m1_col`;

    - Robinson form: `y_col`, `t_col`, `e_col`, `m_col`.

2.  Matrix interface: pass `X`, one of `X0` or `x0`, and
    `pseudo_outcome`.

As in
[`weights_dnn()`](https://JakobJuergens.github.io/RdistNN/reference/weights_dnn.md),
exact DNN localization uses the full rank permutation at each target
rather than an `s`-nearest-neighbor truncation. The `"fast"` jackknife
option reuses the full-sample rank matrix and updates ranks under
deletion. For exact DNN fits, the finite-population support edge is rank
`n - s + k`, not rank `s`.

`tuning_mode = "default"` evaluates only the default kernel order and
preserves the ordinary scalar-fit return shape. `tuning_mode = "grid"`
evaluates the feasible candidates in `tuning`, or the paper default
tuning grid when `tuning` is omitted, and attaches a tidy `tuning_grid`
payload. The pseudo-outcome and nuisance columns are resolved once
before grid evaluation. The primary `tau`, `se`, `meta`, and optional
extras such as `W`, `ranks`, `ess`, `psi`, and `nuisance` remain
anchored to the default candidate. Full per-candidate weight matrices
are not stored. Candidate diagnostics include `ess`, `max_weight`,
`top_n`, `top5_weight_share`, and `degenerate`; `top_n` is fixed at 5
and is independent of the estimator's `k`.

Neyman orthogonality controls first-order population nuisance drift, but
it does not by itself validate arbitrary nuisance predictions fitted on
the same observations used for localization and inference. Valid
full-sample reuse additionally requires the nuisance-rate and full-state
stability conditions in the accompanying paper. This function computes
the estimator and its requested jackknife procedure; it cannot diagnose
whether an arbitrary learner satisfies those conditions. Supply oracle
quantities, cross-fitted or out-of-fold predictions, or full-sample
predictions from a learner for which the paper's rate and stability
conditions are credible.

## See also

[`weights_dnn()`](https://JakobJuergens.github.io/RdistNN/reference/weights_dnn.md),
[`weights_matrix()`](https://JakobJuergens.github.io/RdistNN/reference/weights_matrix.md),
[`score_cate_rho()`](https://JakobJuergens.github.io/RdistNN/reference/score_cate_rho.md),
[`cate_demo_data()`](https://JakobJuergens.github.io/RdistNN/reference/cate_demo_data.md).

## Examples

``` r
set.seed(1)
n <- 50
x_mat <- matrix(rnorm(n * 2), n, 2)
x0_mat <- matrix(rnorm(6), 3, 2)
psi <- rnorm(n)
out <- cate_dnn(
  X = x_mat, X0 = x0_mat, pseudo_outcome = psi,
  distance = "l2", method = "exact", k = 1L, s = 3L,
  se = FALSE, return = c("W", "ess")
)
out$tau
#> [1] -0.02293282 -0.04260484 -0.24413829
```
