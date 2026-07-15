# Mahalanobis distance spec (with one-time covariance setup)

Returns a distance spec that estimates the sample covariance from the
training matrix X in `setup(X)`, then evaluates \$\$d(X_i, x_0) =
\sqrt{(X_i - x_0)^\top \hat\Sigma^{-1}(X_i - x_0)}\$\$ in
[`eval()`](https://rdrr.io/r/base/eval.html) / `eval_batch()`. Requires
at least two covariates.

## Usage

``` r
distance_mahalanobis()
```

## Value

A distance spec (list with `name`, `setup`, `eval`, `eval_batch`) for
use with
[`cate_dnn()`](https://JakobJuergens.github.io/RdistNN/reference/cate_dnn.md),
[`weights_dnn()`](https://JakobJuergens.github.io/RdistNN/reference/weights_dnn.md),
or
[`register_distance()`](https://JakobJuergens.github.io/RdistNN/reference/register_distance.md).
The string alias `"mahalanobis"` is registered at package load.

## Details

The covariance inverse is computed once in `setup()` and reused for all
evaluations via the `ctx` argument, following the same convention as
[`distance_l2_scaled()`](https://JakobJuergens.github.io/RdistNN/reference/distance_l2_scaled.md).
When the sample covariance is (near-)singular a pseudoinverse via SVD is
used automatically; no external packages are required.

Particularly useful for DGPs with correlated covariates (e.g.\\
`"dgp5_correlated"`) where L2 or diagonal-scaled L2 distances fail to
adapt to the off-diagonal covariance structure.

## See also

[`distance_l2()`](https://JakobJuergens.github.io/RdistNN/reference/distance_l2.md),
[`distance_l2_scaled()`](https://JakobJuergens.github.io/RdistNN/reference/distance_l2_scaled.md),
[`cate_dnn()`](https://JakobJuergens.github.io/RdistNN/reference/cate_dnn.md)

## Examples

``` r
set.seed(1)
X  <- matrix(rnorm(100), 50, 2)
x0 <- c(0, 0)
dspec <- distance_mahalanobis()
ctx   <- dspec$setup(X)
head(dspec$eval(X, x0, ctx))
#> [1] 0.8446870 0.6615278 1.0527165 2.2075438 1.5473831 2.2365731
```
