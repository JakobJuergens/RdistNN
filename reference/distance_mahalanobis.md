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
engine. SVD fitting and projection rounding can change represented
distance equality or ordering across packages or platforms, even for
well-conditioned full-rank designs away from the cutoff. Ranks,
exact-distance groups, weights and estimates may then differ. Group
averaging uses each fit's computed distances and cannot reconcile these
differences. See the [numerical geometry
contract](https://github.com/JakobJuergens/RdistNN/blob/master/markdown_notes/design/18_geometry_and_covariance_arithmetic.md).

Particularly useful for DGPs with correlated covariates (e.g.
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
