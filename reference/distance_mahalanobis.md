# Mahalanobis distance spec (with one-time factorization)

Returns a distance spec that fits the sample covariance geometry from
the training matrix X in `setup(X)`, then evaluates \$\$d(X_i, x_0) =
\sqrt{(X_i - x_0)^\top \hat\Sigma^{-1}(X_i - x_0)}\$\$ in
[`eval()`](https://rdrr.io/r/base/eval.html) / `eval_batch()`. Requires
at least two observations and one covariate. The inverse in this formula
is the Moore-Penrose inverse when the centered design is numerically
rank deficient.

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

A whitening factor is computed once by SVD of the binary-normalized,
centered training design, avoiding an explicit covariance inverse and
squaring the design's condition number. Singular values greater than
`max(nrow(X), ncol(X)) * .Machine$double.eps * max(singular_values)` are
retained. Null directions contribute zero distance; an all-constant
training design gives zero distances. The cutoff uses variation in the
original covariate units; changing individual column units can therefore
change the retained numerical rank. The factor is reused through `ctx`,
including for subsets and new targets, following
[`distance_l2_scaled()`](https://JakobJuergens.github.io/RdistNN/reference/distance_l2_scaled.md).
Euclidean norms after projection use the selected C++ or R distance
engine. The explicit rank policy can change neighbors relative to older
fits that truncated covariance singular values after an inverse failed.
Different linear-algebra libraries can still differ at a numerical rank
boundary.

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
