# Observation-level orthogonal PLQR scores at a supplied coefficient

Computes, for each observation \\i\\, the orthogonal PLQR score \$\$
\psi\_{i} \\=\\ (\tau - \mathbf{1}\\\\Y\_{i} \le D\_{i}^{\top}\theta +
g\_{0}(U\_{i})\\\\)\\ (D\_{i} - m\_{D}(U\_{i})) \\ , \$\$ from supplied
nuisance columns. No weights are applied and no localization is
performed;
[`plqr_dnn()`](https://JakobJuergens.github.io/RdistNN/reference/plqr_dnn.md)
aggregates these scores with DNN rank weights.

## Usage

``` r
plqr_score_matrix(data, y_col, D_cols, g_col, mD_cols, theta, tau)
```

## Arguments

- data:

  A `data.frame` containing the outcome, regressor, and nuisance
  columns.

- y_col:

  Name or index of the outcome column in `data`.

- D_cols:

  Names or indices of partially linear regressor columns in `data`.

- g_col:

  Name or index of the shape-nuisance column in `data` holding the
  values of \\g_0(U)\\.

- mD_cols:

  Names or indices of nuisance columns for `m_D(U)`; must have the same
  length and order as `D_cols`.

- theta:

  Numeric coefficient vector with length `length(D_cols)`.

- tau:

  Quantile level in `(0,1)`.

## Value

A numeric matrix of size \\n \times p_d\\ of observation-level scores.

## See also

[`plqr_dnn()`](https://JakobJuergens.github.io/RdistNN/reference/plqr_dnn.md),
[`plqr_demo_data()`](https://JakobJuergens.github.io/RdistNN/reference/plqr_demo_data.md).

## Examples

``` r
demo <- plqr_demo_data(n = 100, p = 2, p_d = 2, tau = 0.5, seed = 1)
scores <- plqr_score_matrix(
  data = demo$data,
  y_col = "Y",
  D_cols = c("D1", "D2"),
  g_col = "g_oracle",
  mD_cols = c("mD1", "mD2"),
  theta = demo$truth$theta_true,
  tau = 0.5
)
head(scores)
#>               D1         D2
#> [1,]  0.12282055  0.2681021
#> [2,] -0.50666199  0.3141894
#> [3,]  0.47597653  0.5914012
#> [4,] -0.09927234 -0.1150896
#> [5,]  0.68557066 -0.4962436
#> [6,]  0.74929848  0.4536638
```
