# Extract one candidate view from a grid fit

`candidate_fit()` turns one row of a grid-mode fit into a scalar-fit
view. The returned object has the same class as `object`, with primary
estimates and standard errors taken from `object$tuning_grid`,
candidate-dependent extras such as `W`, `ess`, and `neighbors` rebuilt
from cached ranks when the parent carried them, and
candidate-independent state carried over. The rebuilt `neighbors` extra
is the first `k` nearest observations, not the full DNN weight support.

## Usage

``` r
candidate_fit(object, candidate, ...)

# Default S3 method
candidate_fit(object, candidate, ...)

# S3 method for class 'np_result'
candidate_fit(object, candidate, ...)

# S3 method for class 'cate_result'
candidate_fit(object, candidate, ...)

# S3 method for class 'plqr_result'
candidate_fit(object, candidate, ...)

candidate_fits(object, candidates = "all", ...)

# Default S3 method
candidate_fits(object, candidates = "all", ...)

# S3 method for class 'np_result'
candidate_fits(object, candidates = "all", ...)

# S3 method for class 'cate_result'
candidate_fits(object, candidates = "all", ...)

# S3 method for class 'plqr_result'
candidate_fits(object, candidates = "all", ...)
```

## Arguments

- object:

  A grid fit returned by
  [`np_dnn()`](https://JakobJuergens.github.io/RdistNN/reference/np_dnn.md),
  [`cate_dnn()`](https://JakobJuergens.github.io/RdistNN/reference/cate_dnn.md),
  or
  [`plqr_dnn()`](https://JakobJuergens.github.io/RdistNN/reference/plqr_dnn.md)
  with `tuning_mode = "grid"`.

- candidate:

  A single `candidate_id` string or integer row index into
  `object$tuning_grid$candidates`.

- ...:

  Reserved for future methods.

- candidates:

  `"all"` (default), a character vector of candidate ids, or an integer
  vector of candidate row indices.

## Value

A scalar-style `np_result`, `cate_result`, or `plqr_result` view.

A named list of scalar-style candidate views.

## Examples

``` r
set.seed(1)
x_mat <- matrix(rnorm(40 * 2), 40, 2)
y_vec <- x_mat[, 1] + rnorm(40, sd = 0.2)
fit <- np_dnn(
  X = x_mat, x0 = c(0, 0), Y = y_vec,
  tuning_mode = "grid", se = FALSE
)
candidate_fit(fit, fit$tuning_grid$candidates$candidate_id[[1L]])
#> <np_result>
#>   data  : n=40, p=2, m=1
#>   DNN   : method=exact, k=1, s=4, distance=l2
#>   est   : mu[1:min(5,m)] = 0.0386  (SEs not requested)
```
