# Exact TDNN rank-weight table (Demirkaya two-scale DNN, k = 1)

Computes the rank-order weight sequence for the two-scale DNN (TDNN)
estimator of Demirkaya et al. as a signed linear combination of two
exact single-scale DNN rank-weight tables evaluated at `(k=1, s1)` and
`(k=1, s2)`:

## Usage

``` r
tdnn_rank_table_exact(n, s1, s2, dim)
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

## Value

A numeric vector of length `n` giving the rank-order TDNN weights
`w_rank[1:n]` (rank 1 = nearest). Entries sum to 1 but some can be
negative.

## Details

\$\$w^{\text{TDNN}}\_j = w_1^\* \\ w^{\text{DNN}}\_j(s_1) + w_2^\* \\
w^{\text{DNN}}\_j(s_2)\$\$

with bias-correction scalars \$\$w_1^\* = \frac{1}{1 -
(s_1/s_2)^{-2/d}}, \quad w_2^\* = 1 - w_1^\*,\$\$ where \\d\\ is the
localization dimension.

**Note:** TDNN weights are signed. They sum to 1 but individual entries
can be negative (specifically, ranks near 1 receive negative weight from
the \\s_1\\ term). This is by design: the signed combination cancels the
leading bias term in the DNN MSE expansion.

## See also

[`tdnn_rank_table_asymp()`](https://JakobJuergens.github.io/RdistNN/reference/tdnn_rank_table_asymp.md),
[`dnn_weights_tdnn_exact_full_from_ranks()`](https://JakobJuergens.github.io/RdistNN/reference/dnn_weights_tdnn_exact_full_from_ranks.md),
[`rank_weight_tdnn_exact()`](https://JakobJuergens.github.io/RdistNN/reference/rank_weight_tdnn_exact.md)

## Examples

``` r
n <- 20; s1 <- 3L; s2 <- 6L; d <- 2L
w <- tdnn_rank_table_exact(n, s1, s2, d)
sum(w)          # 1
#> [1] 1
any(w < 0)      # TRUE (bias-correction weights)
#> [1] TRUE
```
