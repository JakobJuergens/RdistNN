# Asymptotic DNN weights from ranks (full U-statistic, WOR)

Compute **asymptotic** DNN (k-NN) weights via the L-statistic
representation over the full U-statistic of order `s`, using Beta
order-statistic increments as an approximation to the exact
finite-population law.

## Usage

``` r
dnn_weights_asymp_full_from_ranks(ranks, k, s, normalize = TRUE)
```

## Arguments

- ranks:

  Integer permutation of `1:n`; `ranks[i]` is the global rank of obs
  `i`.

- k:

  Integer, number of nearest neighbors used inside each s-subset; must
  satisfy `1 <= k <= s`.

- s:

  Integer, U-statistic order (subset size); must satisfy `1 <= s <= n`.

- normalize:

  Logical; if TRUE (default) renormalizes weights to sum exactly 1
  (guards roundoff).

## Value

A list with

- `w`: numeric length-`n` weights in data order (sum to 1).

- `w_rank`: numeric length-`n` weights in rank order (1..n).

- `meta`: list with `n`, `k`, `s`, `approx = "beta(OrderStat)"`.

## Details

Unlike the exact finite-population rule, this asymptotic approximation
is evaluated on the full rank grid `1:n` and is therefore typically
treated as a full-support rank-weight sequence.

## Examples

``` r
set.seed(1)
n <- 50
s <- 10
k <- 3
ranks <- sample.int(n)
out <- dnn_weights_asymp_full_from_ranks(ranks, k = k, s = s)
sum(out$w) # 1
#> [1] 1
sum(out$w_rank) # 1
#> [1] 1
```
