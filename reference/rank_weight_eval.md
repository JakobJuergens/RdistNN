# Evaluate rank-weight for a single target

Evaluate rank-weight for a single target

## Usage

``` r
rank_weight_eval(spec, ranks, ctx = NULL)
```

## Arguments

- spec:

  Rank-weight spec (see
  [`is_rank_weight()`](https://JakobJuergens.github.io/RdistNN/reference/is_rank_weight.md)).

- ranks:

  Integer permutation of `1..n` (length `n`).

- ctx:

  Optional context (if the spec uses it).

## Value

Numeric vector of length `n`.

## Examples

``` r
rw <- as_rank_weight(function(r) rep(1 / length(r), length(r)))
rank_weight_eval(rw, c(2L, 1L, 3L))
#> [1] 0.3333333 0.3333333 0.3333333
```
