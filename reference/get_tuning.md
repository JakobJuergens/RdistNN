# Retrieve a registered tuning spec

Retrieve a registered tuning spec

## Usage

``` r
get_tuning(name)
```

## Arguments

- name:

  Character scalar.

## Value

A tuning spec.

## Examples

``` r
get_tuning("paper_undersmooth")
#> $name
#> [1] "paper_undersmooth"
#> 
#> $eval
#> function (n, kx, k, se, jackknife, ctx = NULL) 
#> {
#>     s <- .compute_s(n = n, kx = kx, kap = lambda)
#>     list(s = s, meta = list(lambda = as.numeric(lambda), formula = "lambda * n^((kx + 1)/(kx + 2))", 
#>         rounding = rounding))
#> }
#> <bytecode: 0x55c9d22e0110>
#> <environment: 0x55c9d225a420>
#> 
#> $eval_grid
#> function (n, kx, k, se, jackknife, ctx = NULL) 
#> {
#>     s_vec <- vapply(lambda_grid, function(kap) .compute_s(n = n, 
#>         kx = kx, kap = kap), integer(1L))
#>     data.frame(lambda = as.numeric(lambda_grid), s = as.integer(s_vec), 
#>         stringsAsFactors = FALSE)
#> }
#> <bytecode: 0x55c9d2259268>
#> <environment: 0x55c9d225a420>
#> 
#> $setup
#> NULL
#> 
#> attr(,"class")
#> [1] "rdistnn_tuning" "list"          
```
