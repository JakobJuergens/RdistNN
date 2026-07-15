# Rank-weight spec: TDNN asymptotic (Demirkaya two-scale DNN, k = 1)

Returns a `rank_weight` spec whose `eval_batch` uses the asymptotic TDNN
rank table from `tdnn_rank_table_asymp()`. Structurally identical to
[`rank_weight_tdnn_exact()`](https://JakobJuergens.github.io/RdistNN/reference/rank_weight_tdnn_exact.md)
but uses the Beta order-statistic approximation, which gives full rank
support (all n entries are non-zero).

## Usage

``` r
rank_weight_tdnn_asymp(s1, s2, dim, normalize = TRUE)
```

## Arguments

- s1:

  Positive integer; smaller U-statistic order (`s1 < s2`).

- s2:

  Positive integer; larger U-statistic order.

- dim:

  Positive integer; localization dimension.

- normalize:

  Logical; renormalize the rank table to sum 1 (default `TRUE`; guards
  floating-point drift).

## Value

A `rank_weight` spec (list with class `"rank_weight"`).

## See also

[`rank_weight_tdnn_exact()`](https://JakobJuergens.github.io/RdistNN/reference/rank_weight_tdnn_exact.md),
`tdnn_rank_table_asymp()`,
[`weights_tdnn_asymp()`](https://JakobJuergens.github.io/RdistNN/reference/weights_tdnn_asymp.md)

## Examples

``` r
rw <- rank_weight_tdnn_asymp(s1 = 3L, s2 = 6L, dim = 2L)
R <- matrix(c(2L, 1L, 3L, 4L, 5L, 6L,
              6L, 2L, 1L, 3L, 4L, 5L), nrow = 6)
W <- rank_weight_eval_batch(rw, R)
colSums(W)  # c(1, 1)
#> [1] 1 1
```
