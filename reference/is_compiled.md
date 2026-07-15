# Does a score have a compiled implementation?

Does a score have a compiled implementation?

## Usage

``` r
is_compiled(score)
```

## Arguments

- score:

  a `ustat_score`

## Value

logical

## Examples

``` r
is_compiled(new_score("r_only", order = 1L, eval = function(idxs, data, ctx, rng) data[idxs]))
#> [1] FALSE
```
