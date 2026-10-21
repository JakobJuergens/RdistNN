# Delete-one jackknife SE for TDNN via prefix/suffix rank-order engine

Computes the jackknife standard error for a TDNN rank-weighted estimator
using the shared centered prefix/suffix engine. It precomputes the TDNN
rank tables at sizes \\n\\ and \\n-1\\ and passes them to
`.dnn_loo_centered()`, as does the `fast` route of
[`.tdnn_jackknife_se_loo_matrix()`](https://JakobJuergens.github.io/RdistNN/reference/dot-tdnn_jackknife_se_loo_matrix.md).

## Usage

``` r
.tdnn_jackknife_se_loo_full(
  x_mat,
  x0_mat,
  signal,
  distance,
  method,
  s1,
  s2,
  dim,
  normalize,
  ctx,
  signal_arg = "signal",
  precomputed_ranks = NULL
)
```

## Arguments

- x_mat:

  Numeric matrix (n x p).

- x0_mat:

  Numeric matrix (m x p) of target points.

- signal:

  Numeric vector or matrix (n x q) of per-observation signals.

- distance:

  Distance spec or registered name.

- method:

  `"exact"` or `"asymp"`.

- s1:

  Integer; smaller U-statistic order.

- s2:

  Integer; larger U-statistic order.

- dim:

  Integer; localization dimension.

- normalize:

  Logical; passed to the asymptotic rank-table constructor.

- ctx:

  Optional context.

- signal_arg:

  Name of the signal argument used in error messages.

- precomputed_ranks:

  Optional precomputed n x m integer rank matrix. When `NULL` (default)
  ranks are computed internally via
  [`weights_tdnn()`](https://JakobJuergens.github.io/RdistNN/reference/weights_tdnn.md).

## Value

Numeric matrix (q x m) of jackknife standard errors.

## Details

The C++ engine is used automatically when the package has been compiled;
otherwise the same centered calculation uses its R fallback.
