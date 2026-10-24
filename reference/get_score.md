# Retrieve a registered score

Retrieve a registered score

## Usage

``` r
get_score(name)
```

## Arguments

- name:

  character scalar

## Value

a `rdistnn_ustat_score`

## Examples

``` r
sc0 <- as_score(function(idxs, data) sum(data[idxs]), s = 2L, name = "tmp_score_get")
register_score("tmp_score_get", sc0)
sc <- get_score("tmp_score_get")
is_score(sc)
#> [1] TRUE
```
