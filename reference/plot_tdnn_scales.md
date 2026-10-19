# Plot the TDNN two-scale decomposition

Overlays the two-scale TDNN estimate (bold) with the single-scale DNN
estimates at `s1` and `s2` (muted) obtained by internal refits on the
stored training data. The comparison shows what the two-scale
bias-cancelling combination does relative to its components.

## Usage

``` r
plot_tdnn_scales(object, X = NULL, Y = NULL, ...)
```

## Arguments

- object:

  A `tdnn_np_result` from
  [`np_tdnn()`](https://JakobJuergens.github.io/RdistNN/reference/np_tdnn.md).

- X, Y:

  Optional training data; used when the fit does not store it.

- ...:

  Ignored; reserved for future options.

## Value

A ggplot object.

## Examples

``` r
if (requireNamespace("ggplot2", quietly = TRUE)) {
  set.seed(1)
  n <- 100
  x_mat <- matrix(rnorm(n * 2), n, 2,
                  dimnames = list(NULL, c("X1", "X2")))
  y_vec <- x_mat[, 1] + rnorm(n, sd = 0.3)
  fit <- np_tdnn(
    X = x_mat, X0 = cbind(seq(-1, 1, length.out = 6), 0), Y = y_vec,
    s1 = 5L, s2 = 10L, se = FALSE
  )
  plot_tdnn_scales(fit)
}
```
