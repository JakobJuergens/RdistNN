# Evaluate rank-weight for many targets (batch)

A scalar-only rule is evaluated and validated separately for each
target. The assembled result remains an `n x m` matrix when `n = 1` or
`m = 1`. A supplied `eval_batch` must itself return that matrix shape.

## Usage

``` r
rank_weight_eval_batch(spec, R, ctx = NULL)
```

## Arguments

- spec:

  Rank-weight spec (see
  [`is_rank_weight()`](https://JakobJuergens.github.io/RdistNN/reference/is_rank_weight.md)).

- R:

  Integer matrix (n x m), each column a permutation of `1..n`.

- ctx:

  Optional context (if the spec uses it).

## Value

Numeric matrix (n x m).

## Examples

``` r
rw <- as_rank_weight(function(r) rep(1 / length(r), length(r)))
rank_matrix <- cbind(c(1L, 2L, 3L), c(2L, 1L, 3L))
rank_weight_eval_batch(rw, rank_matrix)
#>           [,1]      [,2]
#> [1,] 0.3333333 0.3333333
#> [2,] 0.3333333 0.3333333
#> [3,] 0.3333333 0.3333333
dim(rank_weight_eval_batch(rw, matrix(1L, nrow = 1L, ncol = 2L)))
#> [1] 1 2
```
