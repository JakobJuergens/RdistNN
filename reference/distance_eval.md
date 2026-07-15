# Evaluate a distance spec at a single target

Evaluate a distance spec at a single target

## Usage

``` r
distance_eval(spec, X, x0, ctx = NULL)
```

## Arguments

- spec:

  Distance spec (see
  [`is_distance()`](https://JakobJuergens.github.io/RdistNN/reference/is_distance.md)).

- X:

  Numeric matrix \\n x p\\.

- x0:

  Numeric vector length \\p\\.

- ctx:

  Optional context; if `NULL` and `spec$setup` exists, it will be called
  as `spec$setup(X)` to create a context.

## Value

Numeric vector length `n`.

## Examples

``` r
x_mat <- matrix(rnorm(12), 4, 3)
x0_vec <- c(0, 0, 0)
distance_eval(distance_l2(), x_mat, x0_vec)
#> [1] 1.664026 1.179849 0.938389 2.523089
```
