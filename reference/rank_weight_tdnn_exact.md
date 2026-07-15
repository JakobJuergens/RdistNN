# Rank-weight spec: TDNN exact (Demirkaya two-scale DNN, k = 1)

Returns a `rank_weight` spec whose `eval_batch` uses the exact TDNN rank
table from `tdnn_rank_table_exact()`. The rank table is computed once
per sample size `n` and cached inside the spec closure, so calling
[`rank_weight_eval_batch()`](https://JakobJuergens.github.io/RdistNN/reference/rank_weight_eval_batch.md)
repeatedly on rank matrices of the same size does not recompute the
table.

## Usage

``` r
rank_weight_tdnn_exact(s1, s2, dim)
```

## Arguments

- s1:

  Positive integer; smaller U-statistic order (`s1 < s2`).

- s2:

  Positive integer; larger U-statistic order.

- dim:

  Positive integer; localization dimension. For standard Euclidean
  ranking this is `ncol(X)`; for custom distances it should be set
  explicitly.

## Value

A `rank_weight` spec (list with class `"rank_weight"`).

## Details

**Signed weights.** TDNN weights are a signed linear combination of two
DNN tables and can be negative. They sum to 1 across all `n`
observations but the negative entries are part of the bias-correction
mechanism and should not be interpreted as a standard probability mass.

## See also

[`rank_weight_tdnn_asymp()`](https://JakobJuergens.github.io/RdistNN/reference/rank_weight_tdnn_asymp.md),
`tdnn_rank_table_exact()`,
[`weights_tdnn_exact()`](https://JakobJuergens.github.io/RdistNN/reference/weights_tdnn_exact.md)

## Examples

``` r
rw <- rank_weight_tdnn_exact(s1 = 3L, s2 = 6L, dim = 2L)
R <- matrix(c(2L, 1L, 3L, 4L, 5L, 6L,
              6L, 2L, 1L, 3L, 4L, 5L), nrow = 6)
W <- rank_weight_eval_batch(rw, R)
colSums(W)  # c(1, 1)
#> [1] 1 1
```
