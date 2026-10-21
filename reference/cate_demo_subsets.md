# Quick helper: sample a batch of subsets (for U-stat experiments)

Draws \\K\\ subsets of size \\s\\ uniformly at random from
\\\\1,\ldots,n\\\\. Within each subset, elements are sampled **without**
replacement (all indices are distinct). Across draws, subsets are
sampled independently (the same index may reappear in different rows).

## Usage

``` r
cate_demo_subsets(n, K = 1000L, s = 8L, seed = NULL)
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

An integer `K x s` matrix with one sorted subset per row, including when
`K = 1` or `s = 1`. The result can be passed to
[`score_eval_batch()`](https://JakobJuergens.github.io/RdistNN/reference/score_eval_batch.md).

## Examples

``` r
subs <- cate_demo_subsets(n = 200, K = 1000, s = 8, seed = 1)
dim(subs)
#> [1] 1000    8
```
