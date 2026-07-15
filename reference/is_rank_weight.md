# Is this a rank_weight spec?

A rank-weight spec is a list with `eval(ranks, ctx)` and/or
`eval_batch(rank_matrix, ctx)`, where `rank_matrix` is an n x m rank
matrix (each column a permutation of 1..n). Specs may also define
`table(n)` as an optional rank-to-weight lookup table and `setup(data)`
for manual or wrapper-level context construction.

## Usage

``` r
is_rank_weight(x)
```

## Arguments

- x:

  Any object.

## Value

Logical.

## Examples

``` r
is_rank_weight(rank_weight_dnn_exact(k = 1, s = 1))
#> [1] TRUE
```
