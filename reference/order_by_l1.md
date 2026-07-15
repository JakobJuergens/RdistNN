# Order rows by L1 distance to a point (stable ties)

Returns distances and a stable ordering. If `rustat_order_by_l1()`
exists, it is used; otherwise `order(dx, seq_along(dx))` provides a
deterministic tie-break.

## Usage

``` r
order_by_l1(X, x0)
```

## Arguments

- X:

  Numeric matrix \\n \times p\\. Each row is an observation.

- x0:

  Numeric vector of length \\p\\. The query point.

## Value

A list with:

- `dx`: numeric vector (length \\n\\)

- `ord`: integer vector of 1-based row indices

## See also

[`order_by_l1_scaled()`](https://JakobJuergens.github.io/RdistNN/reference/order_by_l1_scaled.md),
[`dist_l1()`](https://JakobJuergens.github.io/RdistNN/reference/dist_l1.md)

## Examples

``` r
set.seed(2)
X <- matrix(rnorm(10), 5, 2)
x0 <- c(0, 0)
order_by_l1(X, x0)
#> $dx
#> [1] 1.0293348 0.8928039 1.8275434 3.1148496 0.2190388
#> 
#> $ord
#> [1] 5 2 1 3 4
#> 
```
