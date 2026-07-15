# Register a rank_weight under a name

If `spec` is a function, it is wrapped with
[`as_rank_weight()`](https://JakobJuergens.github.io/RdistNN/reference/as_rank_weight.md).

## Usage

``` r
register_rank(name, spec)
```

## Arguments

- name:

  Character scalar registry key.

- spec:

  A rank_weight spec or function `(ranks, ctx) -> numeric(n)`.

## Value

Invisibly, the registered spec.

## Examples

``` r
nm <- "tmp_rank_demo"
register_rank(nm, function(ranks) rep(1 / length(ranks), length(ranks)))
nm %in% rank_names()
#> [1] TRUE
```
