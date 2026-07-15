# Asymptotic TDNN rank-weight table (Demirkaya two-scale DNN, k = 1)

Computes the rank-order weight sequence for the two-scale DNN (TDNN)
estimator using the Beta order-statistic approximation for the
underlying single-scale tables. Structurally identical to
[`tdnn_rank_table_exact()`](https://JakobJuergens.github.io/RdistNN/reference/tdnn_rank_table_exact.md)
but uses the asymptotic approximation, which gives full rank support
(all `n` entries are non-zero).

## Usage

``` r
tdnn_rank_table_asymp(n, s1, s2, dim, normalize = TRUE)
```

## Arguments

- n:

  Positive integer; sample size.

- s1:

  Positive integer; smaller U-statistic order (`s1 < s2 <= n`).

- s2:

  Positive integer; larger U-statistic order (`s1 < s2 <= n`).

- dim:

  Positive integer; localization dimension (typically `ncol(X)`).

- normalize:

  Logical; renormalize so entries sum exactly to 1 (default `TRUE`;
  guards floating-point drift).

## Value

A numeric vector of length `n` giving the rank-order TDNN asymptotic
weights `w_rank[1:n]`. Entries sum to 1 but some can be negative.

## Details

\$\$w^{\text{TDNN,asymp}}\_j = w_1^\* \\ \Delta_j\text{Beta}(1, s_1) +
w_2^\* \\ \Delta_j\text{Beta}(1, s_2)\$\$

where \\\Delta_j\text{Beta}(1, s) = F\_{\text{Beta}(1,s)}(j/n) -
F\_{\text{Beta}(1,s)}((j-1)/n)\\ and the bias-correction scalars
\\w_1^\*, w_2^\*\\ are the same as in
[`tdnn_rank_table_exact()`](https://JakobJuergens.github.io/RdistNN/reference/tdnn_rank_table_exact.md).

## See also

[`tdnn_rank_table_exact()`](https://JakobJuergens.github.io/RdistNN/reference/tdnn_rank_table_exact.md),
[`dnn_weights_tdnn_asymp_full_from_ranks()`](https://JakobJuergens.github.io/RdistNN/reference/dnn_weights_tdnn_asymp_full_from_ranks.md),
[`rank_weight_tdnn_asymp()`](https://JakobJuergens.github.io/RdistNN/reference/rank_weight_tdnn_asymp.md)

## Examples

``` r
n <- 20; s1 <- 3L; s2 <- 6L; d <- 2L
w <- tdnn_rank_table_asymp(n, s1, s2, d)
sum(w)      # 1
#> [1] 1
any(w < 0)  # TRUE
#> [1] TRUE
```
