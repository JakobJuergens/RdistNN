# Quick helper: sample a batch of subsets (singleton by default)

Quick helper: sample a batch of subsets (singleton by default)

## Usage

``` r
np_demo_subsets(n, K = 1000L, s = 1L, seed = NULL)
```

## Arguments

- n:

  Integer; number of observations.

- K:

  Integer; number of subsets to sample.

- s:

  Integer; subset size (order of the kernel).

- seed:

  Optional integer seed for reproducibility.

## Examples

``` r
np_demo_subsets(n = 10, K = 4, s = 1, seed = 1)
#>      [,1] [,2] [,3] [,4]
#> [1,]    9    4    7    1
```
