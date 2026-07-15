# Slow delete-one jackknife for rank-weighted L-stat rules (fixed ranks)

A reference implementation that recomputes *aggregates* per deletion
using the provided neighbor ranks. This **does not** recompute ranks
after deletion (useful as a correctness oracle and for testing).

## Usage

``` r
jackknife_slow_rank(ranks, a, phi, normalize = TRUE)
```

## Arguments

- ranks:

  Integer matrix (T x (s+1)): for each target `t`, `ranks[t, k]` gives
  the positive 1-based index of the k-th neighbor, for
  `k = 1, ..., s+1`. (Must exclude `t` itself if targets are a subset of
  reference.)

- a:

  Finite numeric length-`s` vector of L-stat coefficients
  \\(a_1,\dots,a_s)\\.

- phi:

  Either a finite numeric matrix (n_ref x q) with row `j` equal to
  \\\phi(j)\\, or a function `phi(j, t)` returning a finite length-`q`
  numeric vector.

- normalize:

  Scalar logical; if `TRUE` (default) divide aggregates by `sum(a)`.

## Value

A list with:

- theta0:

  numeric length-`q` baseline aggregate

- theta_loo:

  numeric matrix (n_ref x q) of delete-one aggregates

- info:

  list with fields `s`, `T`, `q`, `normalize`

## Details

Like
[`jackknife_fast_rank()`](https://JakobJuergens.github.io/RdistNN/reference/jackknife_fast_rank.md),
this helper targets compact-support finite-rank L-stat rules. It is not
the internal engine used by the user-level exact DNN wrappers, which
instead update the full rank matrix under deletion and reevaluate the
DNN rank-weight rule.

Works with:

- `phi` as a numeric **matrix** (n_ref x q): row `j` equals \\\phi(j)\\.

- `phi` as a **function** `function(j, t)` that returns a length-`q`
  numeric, allowing target-dependent contributions.

## Examples

``` r
set.seed(1)
n <- 20
p <- 3
s <- 4
X <- matrix(rnorm(n * p), n, p)
# toy ranks: nearest-neighbor indices precomputed elsewhere (T x (s+1))
# here we fake something deterministic just for the example
ranks <- cbind(matrix(rep(1:s, each = n), n, s), s + 1L)
y <- rnorm(n)
a <- rep(1, s)
out <- jackknife_slow_rank(ranks, a, cbind(y))
str(out)
#> List of 3
#>  $ theta0   : num 0.77
#>  $ theta_loo: num [1:20, 1] -0.0162 0.594 0.4118 0.5772 0.77 ...
#>  $ info     :List of 4
#>   ..$ s        : int 4
#>   ..$ T        : int 20
#>   ..$ q        : int 1
#>   ..$ normalize: logi TRUE
```
