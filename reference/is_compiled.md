# Does a score carry compiled-score pointer metadata?

Checks whether `ptr` is non-NULL. This does not validate the pointer or
indicate which evaluator will run;
[`score_eval()`](https://JakobJuergens.github.io/RdistNN/reference/score_eval.md)
and
[`score_eval_batch()`](https://JakobJuergens.github.io/RdistNN/reference/score_eval_batch.md)
use the R functions stored on the score.

## Usage

``` r
is_compiled(score)
```

## Arguments

- score:

  a `ustat_score`

## Value

logical

## Examples

``` r
is_compiled(new_score("r_only", order = 1L, eval = function(idxs, data, ctx, rng) data[idxs]))
#> [1] FALSE
```
