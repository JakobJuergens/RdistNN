# Register a score

Register a score

## Usage

``` r
register_score(name, score)
```

## Arguments

- name:

  Character scalar registry key.

- score:

  A `rdistnn_ustat_score` object.

## Value

Invisibly, the registered name.

## Examples

``` r
sc <- as_score(function(idxs, data) sum(data[idxs]), s = 2L, name = "tmp_score_demo")
register_score("tmp_score_demo", sc)
"tmp_score_demo" %in% list_scores()
#> [1] TRUE
```
