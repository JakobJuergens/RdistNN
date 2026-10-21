# Delete-one jackknife SE for TDNN rank-weighted estimators (internal)

Mirrors `.dnn_jackknife_se_loo_matrix()` with the DNN rank-weight spec
replaced by the TDNN rank-weight spec.

## Usage

``` r
.tdnn_jackknife_se_loo_matrix(
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
  jackknife = c("fast", "loo"),
  precomputed_ranks = NULL
)
```

## Arguments

- x_mat:

  Numeric matrix (n x p).

- x0_mat:

  Numeric matrix (m x p) of target points.

- signal:

  Numeric vector/matrix of per-observation signal values (n x q).

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

  Logical; passed to the asymptotic rank-weight constructor.

- ctx:

  Optional context.

- signal_arg:

  Name of signal argument for error messages.

- jackknife:

  One of `"fast"` (default) or `"loo"`.

- precomputed_ranks:

  Optional precomputed n x m integer rank matrix (only used by the fast
  path; if `NULL` ranks are computed internally).

## Value

Numeric matrix (q x m) of jackknife standard errors.

## Details

**Fast path** (`jackknife = "fast"`): computes full-sample ranks and the
deletion weight table once, then uses the shared centered prefix/suffix
engine. No distances are recomputed.

**LOO fallback** (`jackknife = "loo"`): for each `i`, recomputes weights
from scratch on the `(n-1)` reduced dataset. Strictly slower; useful for
validation.
