# Is object a `ustat_score`?

Is object a `ustat_score`?

## Usage

``` r
is_score(x)
```

## Arguments

- x:

  any object

## Value

logical

## Examples

``` r
sc <- as_score(function(idxs, data) sum(data[idxs]), s = 2L)
is_score(sc)
#> [1] TRUE
```
