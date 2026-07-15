# Construct a distance spec

Construct a distance spec

## Usage

``` r
new_distance(name, eval = NULL, eval_batch = NULL, setup = NULL)
```

## Arguments

- name:

  Character scalar.

- eval:

  Optional function `(X, x0, ctx) -> numeric(n)`.

- eval_batch:

  Optional function `(X, X0, ctx) -> numeric(n x m)`.

- setup:

  Optional function `(X) -> ctx`.

## Value

A validated `dnn_distance` spec.

## Examples

``` r
spec <- new_distance(
  name = "sq_l2",
  eval = function(X, x0, ctx = NULL) rowSums((X - matrix(x0, nrow(X), ncol(X), TRUE))^2)
)
is_distance(spec)
#> [1] TRUE
```
