# Print a `ustat_score` object

Print a `ustat_score` object

## Usage

``` r
# S3 method for class 'ustat_score'
print(x, ...)
```

## Arguments

- x:

  A `ustat_score` object.

- ...:

  Not used.

## Value

`x`, invisibly.

## Examples

``` r
sc <- as_score(function(idxs, data) sum(data[idxs]), s = 2L, name = "sum2")
print(sc)
#> <ustat_score: sum2>
#>   order: 2
#>   impl : R eval 
#>   diff : none
```
