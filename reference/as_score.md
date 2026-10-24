# Coerce a plain function to a score (simple cases)

Coerce a plain function to a score (simple cases)

## Usage

``` r
as_score(f, s, name = deparse(substitute(f)))
```

## Arguments

- f:

  function(idxs, data, ctx, rng) -\> scalar

- s:

  integer order

- name:

  character scalar used as the score name (optional)

## Value

A `rdistnn_ustat_score`

## See also

[`new_score()`](https://JakobJuergens.github.io/RdistNN/reference/new_score.md),
[`score_eval()`](https://JakobJuergens.github.io/RdistNN/reference/score_eval.md)

## Examples

``` r
sc <- as_score(function(idxs, data) prod(data[idxs]), s = 2L, name = "prod2")
score_eval(sc, c(1L, 3L), c(2, 4, 6), NULL)
#> [1] 12
```
