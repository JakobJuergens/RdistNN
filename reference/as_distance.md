# Wrap a function as a distance spec

Wrap a function as a distance spec

## Usage

``` r
as_distance(fun, name = "user_distance")
```

## Arguments

- fun:

  Function `(X, x0, ctx) -> numeric(n)`.

- name:

  Character name for the spec.

## Value

A distance spec with class `rdistnn_distance`.

## Examples

``` r
spec <- as_distance(
  function(X, x0, ctx = NULL) {
    rowSums(abs(X - matrix(x0, nrow(X), ncol(X), TRUE)))
  }
)
x_mat <- matrix(1:9, nrow = 3)
distance_eval(spec, x_mat, c(0, 0, 0))
#> [1] 12 15 18
```
