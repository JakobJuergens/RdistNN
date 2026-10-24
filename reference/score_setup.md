# Set up score context (precomputation)

Set up score context (precomputation)

## Usage

``` r
score_setup(score, data, ...)
```

## Arguments

- score:

  rdistnn_ustat_score

- data:

  Data object understood by the score, such as a vector, matrix or data
  frame.

- ...:

  passed to score\$setup

## Value

context object or NULL

## Examples

``` r
sc <- new_score(
  name = "ctx_sum",
  order = 1L,
  setup = function(data, ...) list(offset = 1),
  eval = function(idxs, data, ctx, rng) data[idxs] + ctx$offset
)
ctx <- score_setup(sc, 1:3)
score_eval(sc, 2L, 1:3, ctx)
#> [1] 3
```
