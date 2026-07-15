# Evaluate a distance spec for many targets (batch)

Evaluate a distance spec for many targets (batch)

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
x_mat <- matrix(rnorm(12), 4, 3)
x0_mat <- matrix(rnorm(6), 2, 3)
distance_eval_batch(distance_l2(), x_mat, x0_mat)
#>          [,1]     [,2]
#> [1,] 1.808042 1.291032
#> [2,] 1.795455 1.221198
#> [3,] 1.977835 2.050604
#> [4,] 4.127344 3.380991
```
