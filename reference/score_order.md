# Resolve kernel order `s`

Resolve kernel order `s`

## Usage

``` r
score_order(score, n, data, ctx = NULL)
```

## Arguments

- score:

  rdistnn_ustat_score

- n:

  integer sample size

- data:

  data used by the score (for dynamic orders)

- ctx:

  optional context

## Value

integer order `s`

## Examples

``` r
sc <- new_score(
  name = "dyn",
  order = function(n, data, ctx) if (n > 5) 3L else 2L,
  eval = function(idxs, data, ctx, rng) sum(data[idxs])
)
score_order(sc, 4L, 1:4, NULL)
#> [1] 2
```
