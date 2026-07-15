# Evaluate kernel on many subsets (R fast-path)

Evaluate kernel on many subsets (R fast-path)

## Usage

``` r
score_eval_batch(score, subs, data, ctx = NULL, rng = NULL)
```

## Arguments

- score:

  ustat_score

- subs:

  integer matrix with nrow=K, ncol=s

- data:

  data.frame/matrix

- ctx:

  optional context

- rng:

  optional RNG key/seed

## Value

numeric vector of length `K`

## Examples

``` r
sc <- as_score(function(idxs, data) sum(data[idxs]), s = 2L)
subs <- rbind(c(1L, 2L), c(2L, 3L))
score_eval_batch(sc, subs, c(3, 4, 5), NULL)
#> [1] 7 9
```
