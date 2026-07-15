# L2 distance spec (single and batched evaluation)

Returns a spec with [`eval()`](https://rdrr.io/r/base/eval.html) and
`eval_batch()` closures for L2 distances. Uses compiled paths when
present.

## Usage

``` r
distance_l2(squared = FALSE)
```

## Arguments

- squared:

  Logical; if `TRUE`, return squared distances.

## Value

A list with:

- `name`: `"l2"` or `"l2_sq"`

- `eval(X, x0, ctx=NULL)`

- `eval_batch(X, X0, ctx=NULL)`

## Batch evaluation

`eval_batch(X, X0)` expects an \\m \times p\\ matrix of queries `X0`.

## Examples

``` r
X <- matrix(rnorm(12), 4, 3)
X0 <- matrix(rnorm(6), 2, 3)
dspec <- distance_l2()
dspec$eval_batch(X, X0)
#>          [,1]     [,2]
#> [1,] 1.758800 3.425654
#> [2,] 2.780884 1.186027
#> [3,] 1.184679 2.962925
#> [4,] 2.149800 3.330523
```
