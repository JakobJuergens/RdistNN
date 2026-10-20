# Joint jackknife covariance across CATE reporting profiles

Joint jackknife covariance across CATE reporting profiles

## Usage

``` r
# S3 method for class 'cate_result'
vcov(object, ...)
```

## Arguments

- object:

  A fitted `cate_result` from
  [`cate_dnn()`](https://JakobJuergens.github.io/RdistNN/reference/cate_dnn.md).

- ...:

  Not used.

## Value

An `m` by `m` covariance matrix, with axes in `coef(object)` order.

## Details

Uses common deleted observations with frozen supplied signals. The
diagonal agrees with the fit's squared standard errors. The `fast` route
freezes fitted geometry; `loo` retains deletion-specific geometry
refitting. Covariance is calculated on request from retained model
state, independently of optional `return` extras, including after
serialization. Fitted uncertainty and an active jackknife are required.
Singular matrices and exact zero covariance are retained. Joint
covariance supplies the arithmetic for prespecified contrasts; valid
contrast inference requires the corresponding joint assumptions and does
not supply simultaneous bands.

## Examples

``` r
fit <- cate_dnn(X = matrix(1:7, ncol = 1),
  X0 = matrix(c(2, 5), ncol = 1),
  pseudo_outcome = c(2, 1, 4, 3, 6, 5, 7), s = 3,
  return = character())
vcov(fit)
#>           target_1  target_2
#> target_1 0.9459184 0.1111224
#> target_2 0.1111224 1.0193878
```
