# CATE score from a precomputed pseudo-outcome \\\rho\\ (pure score; no weights)

Constructs a **pure-score** kernel for CATE experiments that uses only
the precomputed pseudo-outcome \\\rho\\. Per observation \\i\\, the
stored contribution is \$\$g_i \equiv \rho_i\\,\$\$ and for a subset \\S
\subset \\1,\ldots,n\\\\ of size \\s\\ the kernel evaluates to the **sum
of contributions** \$\$h(S) \\=\\ \sum\_{i \in S} g_i \\=\\ \sum\_{i \in
S} \rho_i\\.\$\$

## Usage

``` r
score_cate_rho(x0, x_cols, s, rho_col = NULL, rho = NULL)
```

## Arguments

- x0:

  Numeric length-\\p\\ target point. Used only to check that
  `length(x0)` matches `length(x_cols)`; it is not used in the kernel
  itself.

- x_cols:

  Integer or character indices of the \\X\\ columns in `data`, used for
  dimension checks with `x0`.

- s:

  Integer subset size (kernel order). Recorded in the score metadata;
  the pure kernel internally computes \\\sum\_{i \in S} \rho_i\\ for any
  subset \\S\\ of size \\s\\.

- rho_col:

  Optional name of the column in `data` that contains the pseudo-outcome
  \\\rho\\. If `NULL`, supply `rho`.

- rho:

  Optional numeric vector (length \\n\\) of \\\rho\\ values to use when
  `rho_col` is `NULL`.

## Value

An `ustat_score` object with:

- `$name`: a descriptive name,

- `$order`: the kernel order \\s\\,

- `$setup`: stores `g <- rho` in the context,

- `$eval`/`$eval_batch`: return \\\sum\_{i \in S} g_i\\ for each subset
  \\S\\ (no weights applied).

## Details

This score does **not** perform any nearest-neighbor ordering, kernel
weighting, or L-stat aggregation. Any weighting/averaging of subset
evaluations is handled outside the score layer (e.g., by the
rank-/distance- weight modules). The arguments `x0`, `x_cols`, and `s`
are accepted for API consistency and basic dimension checks; only
\\\rho\\ is used to construct the score.

## Examples

``` r
demo <- cate_demo_data(n = 200, p = 2, design = "hetero", seed = 1)
sc <- score_cate_rho(
  x0 = rep(0, 2), x_cols = paste0("X", 1:2),
  s = 5, rho_col = "aipw_oracle"
)
ctx <- score_setup(sc, demo$data)
subs <- cate_demo_subsets(n = nrow(demo$data), K = 3, s = 5, seed = 1)
score_eval_batch(sc, subs, demo$data, ctx) # three subset evaluations
#> [1] 3.187566 2.527508 3.965492
```
