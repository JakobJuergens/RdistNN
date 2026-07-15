# Retrieve a registered rank_weight

Retrieve a registered rank_weight

## Usage

``` r
get_rank(name)
```

## Arguments

- name:

  Character scalar.

## Value

A rank_weight spec.

## Examples

``` r
register_rank("tmp_rank_get", function(ranks) rep(1 / length(ranks), length(ranks)))
is_rank_weight(get_rank("tmp_rank_get"))
#> [1] TRUE
```
