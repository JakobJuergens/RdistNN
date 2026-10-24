# Plot local quantile calibration for PLQR fits (experimental)

For each target, the weighted share of observations below the quantile
threshold, `sum_i w_ij 1{Y_i <= D_i' theta + g_i}`, is plotted against
the nominal level `tau` with a reference band
`tau +/- z * sqrt(tau (1 - tau) / n_eff)`. This is a descriptive
calibration diagnostic for the supplied threshold; the band does not
test the full coefficient moment condition or account for
fitted-nuisance uncertainty.

## Usage

``` r
plot_calibration(object, level = 0.95)
```

## Arguments

- object:

  A `rdistnn_plqr_result` with the weight matrix `W` stored (the
  default).

- level:

  Confidence level of the reference band (default `0.95`).

## Value

A ggplot object.

## Details

The below-threshold indicator is taken from the fit's retained state.
Refit if this indicator is missing or invalid.

## Lifecycle

**Experimental.** The PLQR API and inference interfaces may change.

## Examples

``` r
if (requireNamespace("ggplot2", quietly = TRUE)) {
  demo <- plqr_demo_data(n = 90, p = 2, p_d = 2, tau = 0.5, seed = 1)
  df <- demo$data
  fit <- plqr_dnn(
    X = as.matrix(df[paste0("X", 1:2)]),
    D = as.matrix(df[paste0("D", 1:2)]),
    g = as.numeric(df$g_oracle),
    mD = as.matrix(df[paste0("mD", 1:2)]),
    Y = as.numeric(df$Y),
    theta_eval = demo$truth$theta_true, tau = 0.5,
    X0 = rbind(c(0, 0), c(0.5, 0.5)), k = 1L, s = 8L, se = FALSE
  )
  plot_calibration(fit)
}
```
