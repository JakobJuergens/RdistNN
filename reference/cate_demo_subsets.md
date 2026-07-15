# Quick helper: sample a batch of subsets (for U-stat experiments)

Draws \\K\\ subsets of size \\s\\ uniformly at random from
\\\\1,\ldots,n\\\\. Within each subset, elements are sampled **without**
replacement (all indices are distinct). Across rows, subsets are sampled
independently (the same index may reappear in different rows).

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

  Optional integer seed for reproducibility.

## Value

An integer matrix of size \\K \times s\\. Each row contains \\s\\
distinct integers in \\\\1,\ldots,n\\\\ representing one subset.

## Examples

``` r
subs <- cate_demo_subsets(n = 200, K = 1000, s = 8, seed = 1)
dim(subs)
#> [1] 1000    8
```
