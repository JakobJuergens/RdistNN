# Nonparametric regression via DNN rank weights

Estimate a localized nonparametric mean by combining observed outcomes
with DNN-based local weights.

## Usage

``` r
np_dnn(
  data = NULL,
  x_cols = NULL,
  y_col = NULL,
  y = NULL,
  X = NULL,
  Y = NULL,
  X0 = NULL,
  x0 = NULL,
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

  Optional data frame containing covariates and outcome columns. Column
  names must be non-missing, non-empty, and unique.

- x_cols:

  Names or indices of covariate columns in `data`.

- y_col:

  Name or index of the outcome column in `data`.

- y:

  Optional numeric outcome vector of length `n`. In the data-frame
  interface, use this instead of `y_col`.

- X:

  Numeric matrix (`n x p`) for the matrix interface.

- Y:

  Numeric outcome vector of length `n` for the matrix interface.

- X0:

  Numeric matrix (`m x p`) of target points; mutually exclusive with
  `x0`.

- x0:

  Numeric vector of length `p` for a single target point.

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

  Optional extras in addition to the always-returned `mu`, `se`, and
  `meta`. Allowed values:
  `c("W", "ranks", "ord", "dx", "y", "ess", "neighbors", "tuning")`.
  `neighbors` is a `k x m` matrix of nearest training-observation
  indices, one column per target. It is a rank summary, not the full DNN
  weight support.

## Value

An object of class `rdistnn_np_result` with components:

- `mu`: estimated local mean values at target points;

- `se`: optional jackknife standard errors;

- `x0`: the `m x p` matrix of evaluation points, with covariate column
  names propagated from the training data when available;

- `meta`: metadata (`n`, `m`, `p`, `k`, `s`, `method`, `distance`);

- `tuning_grid`: in grid mode, a list with `candidates`, `estimates`,
  `diagnostics`, and `meta` tables. Candidate diagnostics report top-5
  weight concentration separately from the estimator's `k`;

- `.model`: retained canonical training state and fitted geometry used
  by prediction and covariance methods, independently of optional
  extras;

- optional extras requested via `return`.

## Details

The documented sampling framework uses independent observations; cluster
and survey-design inference are not supplied. Observations with exactly
equal computed distance to a target share their block's average rank
weight, even when their covariates differ. Groups are target-specific
and use exact represented equality, without a tolerance. Ordering uses
distance then input row index; outcomes and nuisances do not define
groups. Every jackknife deletion averages the surviving blocks with the
sample-size-n-1 rule. A warning is emitted once per public call with
distance ties: the accompanying theory assumes a continuous distance
distribution and does not establish inference under persistent ties.
Refit saved models after package updates during experimental
development.

Two interfaces are supported.

1.  Data-frame interface (recommended): pass `data`, `x_cols`, and
    either `y_col` or `y` (exactly one).

2.  Matrix interface: pass `X`, one of `X0` or `x0`, and `Y`.

With `data`, omit `X` and `Y`; use exactly one of `y_col` or `y`. Column
selectors require `data`. Conflicting routes are rejected.

For `method = "exact"`, DNN weights are computed from the full rank
permutation at each target. The unaveraged rank-table support extends
through rank `n - s + k` (`n - s + 1` in the paper-faithful `k = 1`
case); this is not an `s`-nearest-neighbor truncation. The `"fast"`
jackknife option reuses the full-sample rank matrix and updates ranks
under deletion.

Averaging can extend nonzero observation weights through the end of a
distance block that straddles the nominal rank-table edge. Effective
sample sizes and weight diagnostics use the resulting observation
weights.

`tuning_mode = "default"` evaluates only the default kernel order and
preserves the ordinary scalar-fit return shape. `tuning_mode = "grid"`
evaluates the feasible candidates in `tuning`, or the paper default
tuning grid when `tuning` is omitted, and attaches a tidy `tuning_grid`
payload. The primary `mu`, `se`, `meta`, and optional extras such as
`W`, `ranks`, `ess`, and `neighbors` remain anchored to the default
candidate. Full per-candidate weight matrices are not stored. Candidate
diagnostics include `ess`, `max_weight`, `top_n`, `top5_weight_share`,
and `degenerate`; `top_n` is fixed at 5 and is independent of the
estimator's `k`.

## See also

[`weights_dnn()`](https://JakobJuergens.github.io/RdistNN/reference/weights_dnn.md),
[`np_demo_data()`](https://JakobJuergens.github.io/RdistNN/reference/np_demo_data.md),
[rdistnn_np_result](https://JakobJuergens.github.io/RdistNN/reference/rdistnn_np_result.md)
for summary, coefficient, joint covariance, confidence interval and
prediction methods.

## Examples

``` r
set.seed(1)
n <- 60
x_mat <- matrix(rnorm(n * 2), n, 2)
y_vec <- x_mat[, 1] - 0.5 * x_mat[, 2] + rnorm(n, sd = 0.3)
out <- np_dnn(
  X = x_mat, x0 = c(0, 0), Y = y_vec,
  distance = "l2", method = "exact", k = 1L, s = 3L,
  se = FALSE, return = c("W", "ess")
)
out$mu
#> [1] 0.05695052
```
