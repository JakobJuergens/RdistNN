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

  rdistnn_ustat_score

- idxs:

  integer vector of length `s`

- data:

  Data object understood by the score, such as a vector, matrix or data
  frame.

- ctx:

  optional context

- rng:

  Optional object passed unchanged to the evaluator. The caller or
  evaluator manages random-number seeding.

## Value

numeric(1)

## Examples

``` r
sc <- as_score(function(idxs, data) sum(data[idxs]), s = 2L)
score_eval(sc, c(1L, 2L), c(3, 4, 5), NULL)
#> [1] 7
```
