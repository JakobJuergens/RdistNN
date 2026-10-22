# Compute weight matrix W via an arbitrary ranking + weight mapping

This is the general, ranker-agnostic top-layer. It consumes a *ranking
provider* (distance, custom function, or precomputed
ranks/orders/distances) and a weighting rule (`rank_weight`), then
returns the weight matrix `W` (n x m).

## Usage

``` r
weights_matrix(
  X,
  X0 = NULL,
  x0 = NULL,
  ranking,
  weight = NULL,
  ctx = NULL,
  return = c("W")
)
```

## Arguments

- X:

  numeric matrix (n x p), reference data.

- X0:

  numeric matrix (m x p), target points (rows). If missing, use `x0`.

- x0:

  optional numeric length-p single target (convenience for m=1).

- ranking:

  one of:

  - a `dnn_distance` spec, a registered distance name (e.g., "l2","l1"),

  - a function `(X, X0, ctx) ->` list with `$ranks`/`$ord`/`$dx`, or an
    `n x m` ranks matrix,

  - precomputed object: list with `$ranks`, `$ord`, or `$dx`, or a ranks
    matrix.

- weight:

  a `rank_weight` spec (or zero-arg factory returning one).

- ctx:

  optional shared context passed to ranking and weight evaluators. On
  the ranking side, if `ctx` is `NULL` and a distance spec defines
  `setup(X)`, that setup path is used automatically. It is evaluated
  once for the target batch, including scalar-only distance
  specifications; singleton samples retain an `n x m` distance matrix.

- return:

  which extras to include besides `W`: any of c("ord","ranks","dx").

## Value

An object of class `weights_result` with components:

- `W` : numeric matrix n x m (columns typically sum to 1),

- optional: `ord`, `ranks`, `dx`,

- `meta`: list with `n, m, p, weight_kind, weight_name`.

## Duplicate-covariate blocks

Rows identical on every supplied raw covariate receive the average of
their assigned rank weights, preserving each block's total mass. This
includes caller-supplied ranks and signed or unnormalized rank rules.
Distinct covariate vectors at equal distance remain separate. A warning
identifies this permutation-invariant modification as outside the
paper's theory, including for `k = 1`. Estimator jackknife deletions
repeat the average over each remaining duplicate block using the
selected rule's sample-size-n-1 table. For a block at consecutive ranks
`a:b`, the full weight is `sum(w_n[a:b]) / (b-a+1)` per member; after
deleting a member, each survivor receives
`sum(w_(n-1)[a:(b-1)]) / (b-a)`. A deleted singleton contributes zero.
The smaller block size must be used in the deletion. An exact support
edge describes the unaveraged rank table; a straddling duplicate block
can extend nonzero observation weights beyond it. Lower-level rank-only
APIs retain their rank-primitive semantics.

## See also

[`weights_dnn()`](https://JakobJuergens.github.io/RdistNN/reference/weights_dnn.md),
[`rank_weight_eval_batch()`](https://JakobJuergens.github.io/RdistNN/reference/rank_weight_eval_batch.md),
[`distance_eval_batch()`](https://JakobJuergens.github.io/RdistNN/reference/distance_eval_batch.md)

## Examples

``` r
set.seed(1)
x_mat <- matrix(rnorm(20), 10, 2)
x0_mat <- matrix(rnorm(6), 3, 2)
rw <- as_rank_weight(function(r) {
  w <- 1 / as.numeric(r)
  w / sum(w)
}, name = "inv_rank")
out <- weights_matrix(x_mat, x0_mat, ranking = "l2", weight = rw, return = c("W", "ranks"))
dim(out$W)
#> [1] 10  3
```
