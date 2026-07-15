# TDNN weights (asymptotic, k = 1), distance-based ranking

Computes Demirkaya two-scale DNN (TDNN) weights using the asymptotic
rank-weight table (`tdnn_rank_table_asymp()`). Full rank support: all
`n` observations receive non-zero weight (though some are negative).

## Usage

``` r
weights_tdnn_asymp(
  X,
  X0 = NULL,
  x0 = NULL,
  distance = "l2",
  s1,
  s2,
  dim = NULL,
  normalize = TRUE,
  ctx = NULL,
  return = c("W")
)
```

## Arguments

- X:

  Numeric matrix (n x p).

- X0:

  Numeric matrix (m x p); mutually exclusive with `x0`.

- x0:

  Numeric length-p single target (convenience for m = 1).

- distance:

  Distance spec or registered name (e.g., `"l2"`).

- s1:

  Positive integer; smaller U-statistic order (`s1 < s2 <= n`).

- s2:

  Positive integer; larger U-statistic order.

- dim:

  Positive integer; localization dimension. Defaults to `ncol(X)` when
  `NULL`.

- normalize:

  Logical; renormalize rank table columns to sum 1 (default `TRUE`).

- ctx:

  Optional context passed through to distance and rank-weight
  evaluators.

- return:

  Character vector of extras to include; any of
  `c("W", "ord", "ranks", "dx")`.

## Value

A `weights_result` object (same structure as
[`weights_dnn_asymp()`](https://JakobJuergens.github.io/RdistNN/reference/weights_dnn_asymp.md)).

## See also

[`weights_tdnn_exact()`](https://JakobJuergens.github.io/RdistNN/reference/weights_tdnn_exact.md),
[`weights_tdnn()`](https://JakobJuergens.github.io/RdistNN/reference/weights_tdnn.md),
[`rank_weight_tdnn_asymp()`](https://JakobJuergens.github.io/RdistNN/reference/rank_weight_tdnn_asymp.md)

## Examples

``` r
set.seed(1)
X  <- matrix(rnorm(30), 15, 2)
x0 <- c(0, 0)
out <- weights_tdnn_asymp(X, x0 = x0, s1 = 3L, s2 = 6L)
colSums(out$W)  # 1
#> [1] 1
```
