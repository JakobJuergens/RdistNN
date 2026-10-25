# L1 distance spec (single and batched evaluation)

Returns a spec with [`eval()`](https://rdrr.io/r/base/eval.html) and
`eval_batch()` closures for L1 distances. Compiled fast paths are used
when available; otherwise R fallbacks are used.

## Usage

``` r
distance_l1()
```

## Value

A list with:

- `name`: `"l1"`

- `eval(X, x0, ctx=NULL)`: numeric vector of length \\n\\

- `eval_batch(X, X0, ctx=NULL)`: numeric matrix \\n \times m\\

## Batch evaluation

`eval_batch(X, X0)` expects an \\m \times p\\ matrix of queries `X0`.

## See also

[`distance_l1_scaled()`](https://JakobJuergens.github.io/RdistNN/reference/distance_l1_scaled.md),
[`distance_eval()`](https://JakobJuergens.github.io/RdistNN/reference/distance_eval.md),
[`weights_matrix()`](https://JakobJuergens.github.io/RdistNN/reference/weights_matrix.md)

## Examples

``` r
set.seed(3)
X <- matrix(rnorm(12), 4, 3)
X0 <- matrix(rnorm(6), 2, 3)
dspec <- distance_l1()
dspec$eval(X, X0[1, ])
#> [1] 0.5551521 2.7661406 1.2500104 1.5785391
dspec$eval_batch(X, X0)
#>           [,1]      [,2]
#> [1,] 0.5551521 2.2886396
#> [2,] 2.7661406 2.7985700
#> [3,] 1.2500104 0.4957488
#> [4,] 1.5785391 3.3120267
```
