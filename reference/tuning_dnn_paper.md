# Paper-default undersmoothing rule for kernel order

Builds a tuning spec implementing \$\$s(\lambda) = \lfloor \lambda
n^{(k_x+1)/(k_x+2)} \rfloor\$\$ with finite-sample default \\\lambda =
0.4\\ and optional candidate grid.

## Usage

``` r
tuning_dnn_paper(
  lambda = 0.4,
  lambda_grid = c(0.3, 0.4, 0.5, 0.75),
  rounding = "floor",
  kx_override = NULL
)
```

## Arguments

- lambda:

  Numeric scalar tuning constant.

- lambda_grid:

  Numeric vector of candidate constants for grid mode.

- rounding:

  One of `"floor"`, `"ceiling"`, or `"round"`.

- kx_override:

  Optional positive integer to override runtime `kx`.

## Value

A `rdistnn_tuning` spec.

## Examples

``` r
tun <- tuning_dnn_paper()
tuning_eval(tun, n = 1000, kx = 2, k = 1)
#> $s
#> [1] 71
#> 
#> $meta
#> $meta$lambda
#> [1] 0.4
#> 
#> $meta$formula
#> [1] "lambda * n^((kx + 1)/(kx + 2))"
#> 
#> $meta$rounding
#> [1] "floor"
#> 
#> $meta$name
#> [1] "paper_undersmooth"
#> 
#> $meta$n
#> [1] 1000
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
