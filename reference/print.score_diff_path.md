# Print a `score_diff_path` object

Print a `score_diff_path` object

## Usage

``` r
# S3 method for class 'score_diff_path'
print(x, ...)
```

## Arguments

- x:

  A `score_diff_path` object created by
  [`score_diff_path()`](https://JakobJuergens.github.io/RdistNN/reference/score_diff_path.md).

- ...:

  Not used.

## Value

`x`, invisibly.

## Examples

``` r
sc <- new_score(
  name = "order1",
  order = 1L,
  eval_param = function(idxs, data, theta, ctx, rng) theta[1] * data[idxs]
)
print(score_diff_path(sc, data = c(1, 2, 3)))
#> <score_diff_path>
#>   score : order1
#>   grad  : auto
#>   fd    : step=1e-06, scheme=central
```
