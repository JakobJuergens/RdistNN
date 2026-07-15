# Is this a distance spec?

A distance spec is a list with `eval(X, x0, ctx)` and/or
`eval_batch(X, X0, ctx)`. It may optionally define `setup(X)` to
precompute a context (e.g., scaling).

## Usage

``` r
is_distance(x)
```

## Arguments

- x:

  Any R object.

## Value

Logical.

## Examples

``` r
is_distance(distance_l2())
#> [1] TRUE
is_distance(list(eval = function(X, x0, ctx) rep(0, nrow(X))))
#> [1] TRUE
```
