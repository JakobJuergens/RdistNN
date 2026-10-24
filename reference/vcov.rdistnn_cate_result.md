# Joint jackknife covariance across CATE reporting profiles

Joint jackknife covariance across CATE reporting profiles

## Usage

``` r
# S3 method for class 'rdistnn_cate_result'
vcov(object, ...)
```

## Arguments

- object:

  A fitted `rdistnn_cate_result` from
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
#> Warning: cate_dnn(): Exact distance ties detected; weights are averaged within distance blocks, including jackknife deletions. The accompanying theory assumes a continuous distance distribution and does not establish inference under persistent ties.
vcov(fit)
#> Warning: vcov.rdistnn_cate_result(): Exact distance ties detected; weights are averaged within distance blocks, including jackknife deletions. The accompanying theory assumes a continuous distance distribution and does not establish inference under persistent ties.
#>           target_1  target_2
#> target_1 1.0322449 0.1215306
#> target_2 0.1215306 0.7526020
```
