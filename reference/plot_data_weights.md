# Extract a tidy weight-profile frame from a fitted RdistNN object

Returns one row per (target point, rank) with the DNN weight, cumulative
weight share, the originating observation index, and — when distances
were stored — the distance to the target. Ranks are assigned by
ascending stored distance when `dx` is available and by descending
weight otherwise; the two orderings coincide for monotone single-scale
DNN weights, while two-scale (TDNN) fits should request
`return = c("W", "dx")` for distance-based ranks. For exact DNN fits the
unaveraged rank-table support edge is rank `n - s + k`, not rank `s`; a
duplicate block can extend nonzero observation weights beyond that edge
after averaging. `s` is the subsample/localization order.

## Usage

``` r
plot_data_weights(object)
```

## Arguments

- object:

  A fitted result object holding the weight matrix `W` (stored by
  default via `return = c("W")`).

## Value

A `data.frame` with columns `target_id` (integer), `rank` (integer,
`1..n`), `obs_id` (integer row index into the training data), `weight`,
`cum_share` (cumulative weight share in rank order), and `distance`
(`NA` when distances were not stored).

## Examples

``` r
set.seed(1)
n <- 60
x_mat <- matrix(rnorm(n * 2), n, 2)
y_vec <- x_mat[, 1] + rnorm(n, sd = 0.3)
fit <- np_dnn(
  X = x_mat, x0 = c(0, 0), Y = y_vec,
  k = 1L, s = 10L, se = FALSE, return = c("W", "dx")
)
head(plot_data_weights(fit))
#>   target_id rank obs_id     weight cum_share  distance
#> 1         1    1      2 0.16666667 0.1666667 0.1877888
#> 2         1    2     60 0.14124294 0.3079096 0.2229032
#> 3         1    3     42 0.11932593 0.4272355 0.2568383
#> 4         1    4     16 0.10048499 0.5277205 0.2948897
#> 5         1    5     26 0.08433562 0.6120561 0.3376483
#> 6         1    6     49 0.07053525 0.6825914 0.4002750
```
