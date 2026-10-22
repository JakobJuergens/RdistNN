# Evaluate a distance spec for many targets (batch)

Scalar-only specifications are supported through the shared fallback.
The result retains both matrix dimensions for singleton samples/targets.

## Usage

``` r
distance_eval_batch(spec, X, X0, ctx = NULL)
```

## Arguments

- spec:

  Distance spec (see
  [`is_distance()`](https://JakobJuergens.github.io/RdistNN/reference/is_distance.md)).

- X:

  Numeric matrix \\n x p\\.

- X0:

  Numeric matrix \\m x p\\ (targets as rows).

- ctx:

  Optional context; if `NULL` and `spec$setup` exists, it will be called
  as `spec$setup(X)` to create a context.

## Value

Numeric matrix \\n x m\\.

## Examples

``` r
set.seed(1)
x_mat <- matrix(rnorm(12), 4, 3)
x0_mat <- matrix(rnorm(6), 2, 3)
distance_eval_batch(distance_l2(), x_mat, x0_mat)
#>          [,1]     [,2]
#> [1,] 0.991542 1.672781
#> [2,] 2.125100 2.813195
#> [3,] 1.669451 1.583606
#> [4,] 2.286328 3.928913
scalar_l1 <- as_distance(function(X, x0) {
  rowSums(abs(sweep(X, 2L, x0, "-")))
})
dim(distance_eval_batch(scalar_l1, x_mat[1L, , drop = FALSE], x0_mat))
#> [1] 1 2
```
