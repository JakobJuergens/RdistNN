# PLQR diagnostic score (component \\j\\); pure score (no weights, no fitting)

Constructs a **pure-score** kernel for PLQR diagnostics that stores the
per-observation score contributions \$\$ \psi\_{i} \\=\\ \Bigl(\tau -
\mathbf{1}\\\\Y\_{i} \le D\_{i}^{\top}\theta +
g\_{0}(U\_{i})\\\\\Bigr)\\ \bigl(D\_{i,j} - m\_{D,ij}(U\_{i})\bigr)\\,
\$\$ and returns subset sums \\\sum\_{i \in S} \psi_i\\ in
`eval`/`eval_batch`. This score performs **no** nearest-neighbor
ordering or L-stat weighting; any aggregation/weights are handled
outside the score layer.

## Usage

``` r
score_plqr_component(
  x0,
  x_cols,
  D_cols,
  g_col,
  mD_cols,
  y_col,
  theta,
  j = 1L,
  tau = 0.5,
  s
)
```

## Arguments

- x0:

  Numeric length-\\p\\ target point. Used only for basic dimension
  checks against `x_cols`; not used in the kernel.

- x_cols:

  Integer or character indices of the \\X\\ columns in `data`; used for
  dimension checks with `x0`.

- D_cols:

  Integer or character indices of the \\D\\ columns in `data` (length
  \\p_d\\).

- g_col:

  Character; name of the column in `data` holding the shape-nuisance
  values \\g_0(U)\\.

- mD_cols:

  Integer or character indices of the \\m\_{D}(U)\\ columns in `data`;
  must have the same length and order as `D_cols`. The \\j\\-th column
  represents \\m\_{D,ij}(U)\\.

- y_col:

  Name or index of the outcome column in `data`.

- theta:

  Numeric vector of length `length(D_cols)` giving the fixed coefficient
  vector \\\theta\\.

- j:

  Integer component index in \\1,\ldots,\mathrm{length}(D\\cols)\\.

- tau:

  Quantile level in \\(0,1)\\.

- s:

  Integer kernel order (subset size). Recorded in the score metadata;
  not used internally by this pure score.

## Value

An `ustat_score` object whose `$setup` stores the vector `score` of
contributions in the context, and whose `$eval`/`$eval_batch` return
subset sums of `score` (no weights applied).

## Examples

``` r
demo <- plqr_demo_data(n = 200, p = 3, p_d = 2, tau = 0.5, seed = 1)
sc <- score_plqr_component(
  x0      = rep(0, 3),
  x_cols  = paste0("X", 1:3),
  D_cols  = paste0("D", 1:2),
  g_col   = "g_oracle",
  mD_cols = paste0("mD", 1:2),
  y_col   = "Y",
  theta   = demo$truth$theta_true,
  j       = 1L,
  tau     = 0.5,
  s       = 5
)
ctx <- score_setup(sc, demo$data)
subs <- cate_demo_subsets(n = nrow(demo$data), K = 3, s = 5, seed = 1)
score_eval_batch(sc, subs, demo$data, ctx)
#> [1]  0.6484292 -0.1609615 -0.7863173
```
