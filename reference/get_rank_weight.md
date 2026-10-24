# Retrieve a registered rdistnn_rank_weight

Retrieve a registered rdistnn_rank_weight

## Usage

``` r
get_rank_weight(name)
```

## Arguments

- name:

  Character scalar.

## Value

A rdistnn_rank_weight spec.

## Examples

``` r
register_rank_weight("tmp_rank_get", function(ranks) rep(1 / length(ranks), length(ranks)))
is_rank_weight(get_rank_weight("tmp_rank_get"))
#> [1] TRUE
```
