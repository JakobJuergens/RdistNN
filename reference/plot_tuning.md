# Plot estimate sensitivity to the kernel order

Evaluates or reuses the candidate kernel orders of a tuning grid
(default: the paper rule's lambda grid via
[`tuning_dnn_paper()`](https://JakobJuergens.github.io/RdistNN/reference/tuning_dnn_paper.md),
always including the fitted `s`) and plots the estimates with confidence
bounds against `s`, marking the fitted kernel order with a filled point
and a dotted guide. For grid fits, the stored `tuning_grid` payload is
consumed directly. For scalar fits, the internal grid engine is used on
the stored model state; public estimators are not refit.

## Usage

``` r
plot_tuning(
  object,
  tuning = tuning_dnn_paper(),
  targets = NULL,
  level = 0.95,
  se = TRUE,
  candidates = "all",
  ...
)

# Default S3 method
plot_tuning(
  object,
  tuning = tuning_dnn_paper(),
  targets = NULL,
  level = 0.95,
  se = TRUE,
  candidates = "all",
  ...
)

# S3 method for class 'np_result'
plot_tuning(
  object,
  tuning = tuning_dnn_paper(),
  targets = NULL,
  level = 0.95,
  se = TRUE,
  candidates = "all",
  ...
)

# S3 method for class 'cate_result'
plot_tuning(
  object,
  tuning = tuning_dnn_paper(),
  targets = NULL,
  level = 0.95,
  se = TRUE,
  candidates = "all",
  ...
)

# S3 method for class 'affine_ratio_result'
plot_tuning(
  object,
  tuning = tuning_dnn_paper(),
  targets = NULL,
  level = 0.95,
  se = TRUE,
  candidates = "all",
  ...
)

# S3 method for class 'plqr_result'
plot_tuning(
  object,
  tuning = tuning_dnn_paper(),
  targets = NULL,
  level = 0.95,
  se = TRUE,
  candidates = "all",
  ...
)
```

## Arguments

- object:

  A fitted `np_result`, `cate_result`, `affine_ratio_result` (including
  `clate_result`), or `plqr_result`.

- tuning:

  A tuning specification for the candidate grid (default
  [`tuning_dnn_paper()`](https://JakobJuergens.github.io/RdistNN/reference/tuning_dnn_paper.md)).

- targets:

  Optional integer vector of target ids (defaults to at most three:
  first, middle, last).

- level:

  Confidence level for the bounds (default `0.95`).

- se:

  Logical; compute jackknife bounds for each candidate when `object` is
  a scalar fit (default `TRUE`). Stored grid fits reuse the standard
  errors already present in `object$tuning_grid`.

- candidates:

  Candidate rows to display. The default `"all"` shows the tuning sweep;
  a character vector of candidate ids or `"default"` can be used to
  restrict the display.

- ...:

  Ignored; reserved for future options.

## Value

A ggplot object.

## Details

Scalar grid evaluation uses the training data stored on the fit; fits
from package versions predating that storage must be re-estimated first.

## Experimental methods

Methods for `plqr_result` objects are experimental, as is all PLQR
functionality in this release; see
[`plqr_dnn()`](https://JakobJuergens.github.io/RdistNN/reference/plqr_dnn.md).

## Examples

``` r
# \donttest{
if (requireNamespace("ggplot2", quietly = TRUE)) {
  set.seed(1)
  n <- 100
  x_mat <- matrix(rnorm(n * 2), n, 2)
  y_vec <- x_mat[, 1] + rnorm(n, sd = 0.3)
  fit <- np_dnn(X = x_mat, x0 = c(0, 0), Y = y_vec, k = 1L)
  plot_tuning(fit)
}

# }
```
