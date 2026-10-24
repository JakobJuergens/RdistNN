# Construct a rank-weight spec

Construct a rank-weight spec

## Usage

``` r
new_rank_weight(name, eval = NULL, eval_batch = NULL, setup = NULL)
```

## Arguments

- name:

  Character scalar.

- eval:

  Optional function `(ranks, ctx) -> numeric(n)`.

- eval_batch:

  Optional function `(rank_matrix, ctx) -> numeric(n x m)`.

- setup:

  Optional function to build context. At the top level,
  [`weights_matrix()`](https://JakobJuergens.github.io/RdistNN/reference/weights_matrix.md)
  passes any supplied `ctx` through to the rank-weight evaluator, but it
  does not automatically call `setup()` for you.

## Value

A validated `rdistnn_rank_weight` specification.

## Examples

``` r
rw <- new_rank_weight(
  name = "inv_rank",
  eval = function(ranks, ctx = NULL) {
    w <- 1 / as.numeric(ranks)
    w / sum(w)
  }
)
is_rank_weight(rw)
#> [1] TRUE
```
