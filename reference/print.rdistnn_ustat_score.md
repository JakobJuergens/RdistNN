# Print a `rdistnn_ustat_score` object

Print a `rdistnn_ustat_score` object

## Usage

``` r
# S3 method for class 'rdistnn_ustat_score'
print(x, ...)
```

## Arguments

- x:

  A `rdistnn_ustat_score` object.

- ...:

  Not used.

## Value

`x`, invisibly.

## Examples

``` r
sc <- as_score(function(idxs, data) sum(data[idxs]), s = 2L, name = "sum2")
print(sc)
#> <rdistnn_ustat_score: sum2>
#>   order: 2
#>   impl : R eval 
#>   diff : none
```
