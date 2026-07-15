# Exact DNN weights from ranks (full U-statistic, WOR)

Compute exact DNN (k-NN) weights for a target point using the
L-statistic representation over the **full** U-statistic of order `s`
(all \\\binom{n}{s}\\ subsets sampled *without* replacement). The
function is always given the permutation of ranks `ranks`, where
`ranks[i]` is the global rank (1 = nearest) of observation `i` relative
to the target.

## Usage

``` r
dnn_weights_exact_full_from_ranks(ranks, k, s)
```

## Arguments

- ranks:

  integer permutation of `1:n`; `ranks[i]` is the rank of obs `i`.

- k:

  integer, number of nearest neighbors used inside each s-subset; must
  satisfy `1 <= k <= s`.

- s:

  integer, U-statistic order (subset size); must satisfy `1 <= s <= n`.

## Value

A list with

- `w`: numeric length-`n` weights in data order (sum to 1).

- `w_rank`: numeric length-`n` weights in rank order (1..n).

- `meta`: list with `n`, `k`, `s`.

## Details

For the exact rule, weights are computed from the full rank permutation,
not by truncating to the first `s` observations. The resulting exact
support in rank order is:

- `j <= n - s + k` in the general `k`-NN-within-subset case;

- in particular, `j <= n - s + 1` when `k = 1`.

## Examples

``` r
set.seed(1)
n <- 10
s <- 4
k <- 2
ranks <- sample.int(n) # ranks[i] in 1..n
out <- dnn_weights_exact_full_from_ranks(ranks, k = k, s = s)
sum(out$w) # 1
#> [1] 1
sum(out$w_rank) # 1
#> [1] 1
```
