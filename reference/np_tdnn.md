# Nonparametric regression via TDNN rank weights

Estimate a localized nonparametric mean by combining observed outcomes
with Demirkaya two-scale DNN (TDNN) local weights. TDNN is a signed
linear combination of two DNN tables at orders `s1 < s2`; the
combination cancels the leading bias term and achieves a faster MSE rate
than single-scale DNN under sufficient smoothness.

## Usage

``` r
np_tdnn(
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
  s1 = 1L,
  s2 = 2L,
  dim = NULL,
  normalize = TRUE,
  se = TRUE,
  jackknife = c("fast", "loo", "none"),
  ctx = NULL,
  return = c("W")
)
```

## Arguments

- data:

  Optional data frame (or tibble) containing covariates and the outcome
  column. Column names must be non-missing, non-empty, and unique.

- x_cols:

  Names or integer indices of covariate columns in `data`.

- y_col:

  Name or integer index of the outcome column in `data`.

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

  One of `"exact"` or `"asymp"` for TDNN rank weights.

- s1:

  Positive integer; smaller U-statistic order. Must satisfy `s1 < s2`.

- s2:

  Positive integer; larger U-statistic order. Must satisfy
  `s1 < s2 <= n`.

- dim:

  Positive integer; localization dimension. Defaults to `ncol(X)` when
  `NULL`.

- normalize:

  Logical; renormalize the asymptotic rank table to sum 1 (only used
  when `method = "asymp"`; default `TRUE`).

- se:

  Logical; compute delete-1 jackknife standard errors.

- jackknife:

  One of `"fast"`, `"loo"`, or `"none"`. `"fast"` (default) uses the
  prefix/suffix rank-order engine; `"loo"` recomputes leave-one-out
  weights from scratch; `"none"` skips SEs entirely.

- ctx:

  Optional context passed through to distance and rank-weight
  evaluators.

- return:

  Optional extras in addition to the always-returned `mu`, `se`, and
  `meta`. Allowed values:
  `c("W", "ranks", "ord", "dx", "y", "ess", "neighbors")`. The
  `neighbors` extra returns the nearest observation for each target; it
  is not the TDNN signed weight support.

## Value

An object of class `tdnn_np_result` with components:

- `mu`: estimated local mean values at target points;

- `se`: optional jackknife standard errors;

- `x0`: the `m x p` matrix of evaluation points, with covariate column
  names propagated from the training data when available;

- `meta`: metadata (`n`, `m`, `p`, `s1`, `s2`, `dim`, `method`,
  `distance`);

- `.model`: retained canonical training state and fitted geometry used
  by the inherited prediction and covariance methods;

- optional extras requested via `return`.

## Details

Two interfaces are supported.

1.  Data-frame interface (recommended): pass `data`, `x_cols`, and
    exactly one of `y_col` or `y`.

2.  Matrix interface: pass `X`, one of `X0` or `x0`, and `Y`.

With `data`, omit `X` and `Y`; use exactly one of `y_col` or `y`. Column
selectors require `data`. Conflicting routes are rejected.

**Signed weights.** The columns of the TDNN weight matrix sum to 1, but
individual entries can be negative. This is expected and reflects the
bias-correction structure; it does not indicate numerical problems.

**Localization dimension.** `dim` governs the rate at which TDNN order
weights decay with rank. It defaults to `ncol(X)` (the covariate
dimension) and should be overridden only when a custom distance has a
different effective dimension.

**Jackknife.** `"fast"` (default) uses the prefix/suffix rank-order
engine
([`jackknife_rank_full()`](https://JakobJuergens.github.io/RdistNN/reference/jackknife_rank_full.md)),
which precomputes the TDNN rank tables at `n` and `n-1` once and runs an
\\O(n \cdot m)\\ loop without re-evaluating the distance. `"loo"`
recomputes weights from scratch for each deletion and is useful for
validation.

## See also

[`weights_tdnn()`](https://JakobJuergens.github.io/RdistNN/reference/weights_tdnn.md),
[`np_dnn()`](https://JakobJuergens.github.io/RdistNN/reference/np_dnn.md),
[`np_demo_data()`](https://JakobJuergens.github.io/RdistNN/reference/np_demo_data.md),
[np_result](https://JakobJuergens.github.io/RdistNN/reference/np_result.md)
for inherited summary, covariance, confidence interval and prediction
methods.

## Examples

``` r
set.seed(1)
n <- 80
x_mat <- matrix(rnorm(n * 2), n, 2)
y_vec <- x_mat[, 1] - 0.5 * x_mat[, 2] + rnorm(n, sd = 0.3)
out <- np_tdnn(
  X = x_mat, x0 = c(0, 0), Y = y_vec,
  distance = "l2", method = "exact", s1 = 2L, s2 = 4L,
  se = FALSE, return = c("W", "ess")
)
out$mu
#> [1] 0.04952335
```
