# Wrap a function as a tuning spec

Wrap a function as a tuning spec

## Usage

``` r
as_tuning(fun, name = "user_tuning")
```

## Arguments

- fun:

  Function that returns either a scalar `s` or a list with field `s`.

- name:

  Character name for the spec.

## Value

A `dnn_tuning` spec.

## Examples

``` r
tun <- as_tuning(function(n, kx, ...) floor(n^((kx + 1)/(kx + 2))))
tuning_eval(tun, n = 200, kx = 2, k = 1)
#> $s
#> [1] 53
#> 
#> $meta
#> $meta$name
#> [1] "user_tuning"
#> 
#> $meta$n
#> [1] 200
#> 
#> $meta$kx
#> [1] 2
#> 
#> 
#> $feasibility
#> $feasibility$k_le_s
#> [1] TRUE
#> 
#> $feasibility$s_lt_n
#> [1] TRUE
#> 
#> 
```
