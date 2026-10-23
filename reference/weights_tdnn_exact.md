# TDNN weights (exact, k = 1), distance-based ranking

Computes Demirkaya two-scale DNN (TDNN) weights using the exact
rank-weight table (`tdnn_rank_table_exact()`). TDNN is a signed linear
combination of two single-scale exact DNN tables at orders `s1` and
`s2`; the combination coefficients are chosen to cancel the leading bias
term in the DNN MSE expansion. The unaveraged table's nominal support
ends at `n - s1 + 1`; averaging can extend observation weights through a
distance block that straddles this edge.

## Usage

``` r
weights_tdnn_exact(
  X,
  X0 = NULL,
  x0 = NULL,
  distance = "l2",
  s1,
  s2,
  dim = NULL,
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

- ctx:

  Optional context passed through to distance and rank-weight
  evaluators.

- return:

  Character vector of extras to include; any of
  `c("W", "ord", "ranks", "dx")`.

## Value

A `weights_result` object (same structure as
[`weights_dnn_exact()`](https://JakobJuergens.github.io/RdistNN/reference/weights_dnn_exact.md)).

## Details

**Signed weights.** The returned weight matrix `W` can contain negative
entries; columns still sum to 1. This is by design and reflects the
bias-correction structure of TDNN.

## Exact distance ties

For each target, observations with exactly equal computed distance share
the average of their assigned rank weights, including distinct covariate
vectors. Equality is exact represented equality, without an epsilon or
rounding. Groups may differ by target. Block mass is preserved for
signed and unnormalized rules. Equal distances are ordered by original
row index. A bare supplied rank or order matrix defines a strict order
and is not averaged. A list containing ranks or orders together with
`dx` preserves that authoritative permutation and averages by the
supplied distances. No outcomes, treatments or nuisance values define
membership.

For a block at ranks `a:b`, each member receives
`sum(w_n[a:b]) / (b-a+1)`. Estimator jackknife deletions remove the
actual observation, update ranks, and average remaining blocks with
`w_(n-1)`. If a member is deleted, each survivor receives
`sum(w_(n-1)[a:(b-1)]) / (b-a)`; a deleted singleton contributes zero.
Exact support edges describe the unaveraged rank table; averaging can
extend nonzero observation weights through a block crossing that edge.
Lower-level rank-only APIs retain their rank-primitive semantics.

A warning is emitted once per public call with distance ties. The
accompanying theory assumes a continuous distance distribution and does
not establish inference under persistent ties, including for `k = 1`.

## See also

[`weights_tdnn_asymp()`](https://JakobJuergens.github.io/RdistNN/reference/weights_tdnn_asymp.md),
[`weights_tdnn()`](https://JakobJuergens.github.io/RdistNN/reference/weights_tdnn.md),
[`rank_weight_tdnn_exact()`](https://JakobJuergens.github.io/RdistNN/reference/rank_weight_tdnn_exact.md)

## Examples

``` r
set.seed(1)
X  <- matrix(rnorm(30), 15, 2)
x0 <- c(0, 0)
out <- weights_tdnn_exact(X, x0 = x0, s1 = 3L, s2 = 6L)
colSums(out$W)  # 1
#> [1] 1
any(out$W < 0)  # TRUE
#> [1] TRUE
```
