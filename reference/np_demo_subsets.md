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

  Optional integer seed, used locally without changing the caller's RNG
  state. `NULL` uses the current stream.

## Value

For `s > 1`, an integer `K x s` matrix with one subset per row. For
`s = 1`, the current helper returns a `1 x K` matrix; reshape it with
`matrix(result, ncol = 1L)` before passing it to
[`score_eval_batch()`](https://JakobJuergens.github.io/RdistNN/reference/score_eval_batch.md).

## Examples

``` r
draws <- np_demo_subsets(n = 10, K = 4, s = 1, seed = 1)
matrix(draws, ncol = 1L)  # one singleton subset per row
#>      [,1]
#> [1,]    9
#> [2,]    4
#> [3,]    7
#> [4,]    1
```
