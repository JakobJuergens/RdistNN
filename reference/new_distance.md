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

A validated `rdistnn_distance` spec.

## Details

A batch evaluator must return an actual `n x m` numeric matrix,
including when either dimension is one. With only `eval`, the shared
batch evaluator assembles those dimensions for you and validates each
result. When `ctx` is `NULL`, `setup(X)` runs once per evaluation call,
not once per target in a batch.
[`weights_matrix()`](https://JakobJuergens.github.io/RdistNN/reference/weights_matrix.md)
uses this same evaluator.

## Examples

``` r
spec <- new_distance(
  name = "sq_l2",
  eval = function(X, x0, ctx = NULL) rowSums((X - matrix(x0, nrow(X), ncol(X), TRUE))^2)
)
is_distance(spec)
#> [1] TRUE
```
