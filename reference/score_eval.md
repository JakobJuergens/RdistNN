# Evaluate kernel on one subset

For maximum robustness, this prefers the scalar
[`eval()`](https://rdrr.io/r/base/eval.html) path when available and
otherwise falls back to `eval_batch()` for a single-row batch.

## Usage

``` r
score_eval(score, idxs, data, ctx = NULL, rng = NULL)
```

## Arguments

- score:

  ustat_score

- idxs:

  integer vector of length `s`

- data:

  data.frame/matrix

- ctx:

  optional context

- rng:

  optional RNG key/seed

## Value

numeric(1)

## Examples

``` r
sc <- as_score(function(idxs, data) sum(data[idxs]), s = 2L)
score_eval(sc, c(1L, 2L), c(3, 4, 5), NULL)
#> [1] 7
```
