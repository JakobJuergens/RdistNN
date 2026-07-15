# Convert an index order to ranks (inverse permutation)

Convert an index order to ranks (inverse permutation)

## Usage

``` r
ord_to_ranks(ord)
```

## Arguments

- ord:

  Integer vector (values `1..n`) giving increasing order.

## Value

Integer ranks vector `r` where `r[i]` is the rank of item `i`.

## Examples

``` r
ord_to_ranks(c(3L, 1L, 2L))
#> [1] 2 3 1
```
