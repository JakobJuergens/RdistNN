# Asymptotic TDNN weights from ranks (Demirkaya two-scale DNN, k = 1)

Given a rank permutation `ranks`, computes the TDNN weight for each
observation using the asymptotic rank-weight table from
[`tdnn_rank_table_asymp()`](https://JakobJuergens.github.io/RdistNN/reference/tdnn_rank_table_asymp.md).

## Usage

``` r
dnn_weights_tdnn_asymp_full_from_ranks(ranks, s1, s2, dim, normalize = TRUE)
```

## Arguments

- ranks:

  Integer permutation of `1:n`; `ranks[i]` is the rank (1 = nearest) of
  observation `i` relative to the target.

- s1:

  Positive integer; smaller U-statistic order.

- s2:

  Positive integer; larger U-statistic order.

- dim:

  Positive integer; localization dimension.

- normalize:

  Logical; renormalize to sum 1 (default `TRUE`).

## Value

A list with

- `w`: numeric length-`n` weights in data order (sum to 1, some entries
  may be negative).

- `w_rank`: numeric length-`n` weights in rank order (1..n).

- `meta`: list with `n`, `s1`, `s2`, `dim`, `normalize`, and `scalars`.

## See also

[`tdnn_rank_table_asymp()`](https://JakobJuergens.github.io/RdistNN/reference/tdnn_rank_table_asymp.md),
[`rank_weight_tdnn_asymp()`](https://JakobJuergens.github.io/RdistNN/reference/rank_weight_tdnn_asymp.md),
[`dnn_weights_tdnn_exact_full_from_ranks()`](https://JakobJuergens.github.io/RdistNN/reference/dnn_weights_tdnn_exact_full_from_ranks.md)

## Examples

``` r
set.seed(1)
n <- 15L; ranks <- sample.int(n)
out <- dnn_weights_tdnn_asymp_full_from_ranks(ranks, s1 = 3L, s2 = 6L, dim = 2L)
sum(out$w)        # 1
#> [1] 1
any(out$w < 0)    # TRUE
#> [1] TRUE
```
