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
  `k = 1, ..., s+1`. Use distinct indices within each row for
  observation deletion. A target may coincide with a reference
  observation; self-exclusion is not required. Repeated indices follow
  the additive convention in
  [`jackknife_fast_rank()`](https://JakobJuergens.github.io/RdistNN/reference/jackknife_fast_rank.md).

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

  numeric length-`q` baseline aggregate, averaged across targets

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
  numeric, allowing target-dependent contributions. With matrix `phi`,
  `n_ref = nrow(phi)`. With function-valued `phi`, it is inferred as
  `max(ranks)`.

## Examples

``` r
set.seed(1)
n <- 20
s <- 3L
X <- matrix(rnorm(n * 2), n, 2)
wt <- weights_dnn(X, X0 = rbind(c(0, 0), c(1, 0)), k = 1L, s = s, return = "ord")
ranks <- t(wt$ord[seq_len(s + 1L), , drop = FALSE])
y <- rnorm(n)
a <- rep(1, s)
out <- jackknife_slow_rank(ranks, a, cbind(y))
str(out)
#> List of 3
#>  $ theta0   : num -0.0484
#>  $ theta_loo: num [1:20, 1] -0.0484 -0.0484 -0.0484 -0.0484 -0.0484 ...
#>  $ info     :List of 4
#>   ..$ s        : int 3
#>   ..$ T        : int 2
#>   ..$ q        : int 1
#>   ..$ normalize: logi TRUE
```
