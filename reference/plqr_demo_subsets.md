# Quick helper: sample a batch of subsets (for U-stat / localization experiments) (experimental)

Quick helper: sample a batch of subsets (for U-stat / localization
experiments) (experimental)

## Usage

``` r
plqr_demo_subsets(n, K = 1000L, s = 8L, seed = NULL)
```

## Arguments

- n:

  Integer; number of observations.

- K:

  Integer; number of subsets to sample.

- s:

  Integer; subset size (order of the kernel).

- seed:

  Optional integer seed, used locally without changing the caller's RNG
  state. `NULL` uses the current stream.

## Value

For `s > 1`, an integer `K x s` matrix with one subset per row. For
`s = 1`, the current helper returns a `1 x K` matrix; reshape it with
`matrix(result, ncol = 1L)` before passing it to
[`score_eval_batch()`](https://JakobJuergens.github.io/RdistNN/reference/score_eval_batch.md).

## Lifecycle

**Experimental.** The PLQR API and inference interfaces may change.

## Examples

``` r
plqr_demo_subsets(n = 10, K = 4, s = 3, seed = 1)
#>      [,1] [,2] [,3]
#> [1,]    4    7    9
#> [2,]    1    2    5
#> [3,]    2    3    7
#> [4,]    1    3    5
```
