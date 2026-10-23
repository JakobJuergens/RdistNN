# Experimental PLQR: moments and coefficient uncertainty

[`plqr_dnn()`](https://JakobJuergens.github.io/RdistNN/reference/plqr_dnn.md)
evaluates localized quantile-score moments at a **supplied coefficient
vector**. It does not solve for coefficients. This guide shows how
moment uncertainty differs from density-backed coefficient uncertainty;
the entire PLQR interface remains experimental.

PLQR uses the shared localization weights: observations at exactly equal
represented selected-metric distance to a target receive the mean of
their assigned rank weights. Each deletion averages its `n-1` table over
surviving groups. Regressors, outcomes and nuisance predictions do not
define the groups. The package warns that inference with persistent
distance ties is not established; see the [weight
guide](https://JakobJuergens.github.io/RdistNN/articles/dnn-weights.html#exact-distance-ties).

For statistical background, see the
[paper](https://jakobjuergens.com/research/#orthogonal-inference-conditional-z).
For the supported CATE, regression and CLATE workflows, use [Getting
started](https://JakobJuergens.github.io/RdistNN/articles/getting-started.md)
and [Conditional
LATE](https://JakobJuergens.github.io/RdistNN/articles/conditional-late.md).

## Supply the predictions and coefficient

The demo provides oracle nuisance values. In an application, supply
aligned external predictions: `g` is the model’s shape component, `mD`
contains the regressor means, and `density` is the conditional outcome
density at the threshold `D' theta_eval + g`. Nuisance fitting is
external.

This example requests both the full moment covariance and a
density-backed coefficient sandwich with `theta_meat = "jackknife"`:

``` r

demo <- plqr_demo_data(n = 160, p = 2, p_d = 2, tau = 0.5, seed = 1)
fit <- plqr_dnn(
  data = demo$data, x_cols = c("X1", "X2"),
  D_cols = c("D1", "D2"), g_col = "g_oracle",
  mD_cols = c("mD1", "mD2"), y_col = "Y",
  density_col = "density_oracle", theta_eval = demo$truth$theta_true,
  tau = 0.5, x0 = c(0, 0), k = 1L, s = 10L,
  theta_meat = "jackknife", se = TRUE
)

data.frame(
  supplied_coefficient = as.numeric(coef(fit)),
  moment = as.numeric(fit$moment),
  moment_se = as.numeric(fit$moment_se),
  coefficient_se = as.numeric(fit$theta_se)
)
#>   supplied_coefficient      moment  moment_se coefficient_se
#> 1                  0.5 -0.12784942 0.05934042     0.06939776
#> 2                  0.5 -0.03771978 0.05977131     0.07589505
```

Use either `data` with column selectors or the direct `X`, `Y`, `D`, `g`
and `mD` inputs; mixing those routes is rejected. A direct `density`
vector is allowed with either route. Use `density_col` only with `data`,
and supply `density` or `density_col`, never both.

The two SE columns have different units. `moment_se` describes variation
in the localized score; `theta_se` additionally uses the density-backed
Jacobian. A small moment is not an estimated coefficient. `coef(fit)`
repeats `theta_eval`, and intervals are centered at that supplied value:

``` r

confint(fit)
#>     target term  conf.low conf.high
#> 1 target_1   D1 0.3639829 0.6360171
#> 2 target_1   D2 0.3512484 0.6487516
```

Without a density input, leave `theta_meat` at its default and use the
fit for moment diagnostics only; coefficient intervals are unavailable.
With a density input, the default `theta_meat = "analytic"` uses
analytic score covariance. See
[`?plqr_dnn`](https://JakobJuergens.github.io/RdistNN/reference/plqr_dnn.md)
for the two sandwich routes and returned arrays.

## Inspect conditioning and new profiles

A well-conditioned regressor Gram matrix does not guarantee a
well-conditioned density-weighted Jacobian. Inspect both:

``` r

d <- as.data.frame(diagnostics(fit))
d[, c("target_id", "gram_condition_number", "jacobian_condition_number",
      "theta_vcov_available")]
#>          target_id gram_condition_number jacobian_condition_number
#> target_1  target_1              1.228876                  1.263684
#>          theta_vcov_available
#> target_1                 TRUE
```

Prediction relocalizes the retained sample. It repeats the supplied
coefficient and computes new uncertainty; it does not fit a new PLQR
root. Prediction uncertainty uses the analytic sandwich even when the
original fit used jackknife meat.

``` r

pred <- predict(fit, X0 = rbind(c(0, 0), c(0.25, -0.25)), se.fit = TRUE)
pred$se.fit
#>      target_1   target_2
#> D1 0.07219991 0.08500969
#> D2 0.07391976 0.08562597
```

For array shapes and scales, see [Notation and
conventions](https://JakobJuergens.github.io/RdistNN/articles/notation-and-conventions.md).
For the optional conditioning and calibration plots, see
[`?plot_conditioning`](https://JakobJuergens.github.io/RdistNN/reference/plot_conditioning.md)
and
[`?plot_calibration`](https://JakobJuergens.github.io/RdistNN/reference/plot_calibration.md).
Intervals are pointwise; numerical conditioning does not verify nuisance
accuracy, localization bias or the paper’s assumptions.
