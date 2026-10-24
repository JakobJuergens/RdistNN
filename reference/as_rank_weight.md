# Wrap a function as a rdistnn_rank_weight spec

Wrap a function as a rdistnn_rank_weight spec

## Usage

``` r
as_rank_weight(fun, name = "user_rank_weight")
```

## Arguments

- fun:

  Function `(ranks, ctx) -> numeric(n)` OR just `(ranks) -> numeric(n)`.

- name:

  Character name.

## Value

A rank-weight spec (list) with fields `name`, `eval`.

## Examples

``` r
rw <- as_rank_weight(function(ranks) {
  w <- 1 / as.numeric(ranks)
  w / sum(w)
})
rank_weight_eval(rw, c(3L, 1L, 2L))
#> [1] 0.1818182 0.5454545 0.2727273
```
