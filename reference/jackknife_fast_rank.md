# Fast delete-one jackknife for rank-weighted functionals (auto C++ when available, R fallback otherwise)

Computes the baseline aggregate `theta0` and all delete-one aggregates
`theta_loo` for generic rank-weighted L-stat rules with finite support
`s`.

## Usage

``` r
jackknife_fast_rank(ranks, a, phi, normalize = TRUE)
```

## Arguments

- ranks:

  Integer matrix (T x (s+1)) of positive neighbor indices (1-based).
  Each row t lists the first s+1 neighbors of target t.

- a:

  Finite numeric length-s vector of L-stat coefficients a_1..a_s.

- phi:

  either:

  - finite numeric matrix (n_ref x q): row j is the q-vector phi(j), or

  - function `function(j, t)` returning a finite numeric length-q vector
    (allows target dependence). For a scalar \\\phi\\, use `cbind(y)` as
    the matrix case.

- normalize:

  Scalar logical; if TRUE (default) divide by sum(a).

## Value

list with:

- `theta0`: numeric length-q baseline aggregate,

- `theta_loo`: numeric matrix (n_ref x q) of delete-one aggregates,

- `info`: list(s, T, q, normalize).

## Details

This utility is for rules that are fully characterized by:

- the first `s + 1` ranked observations for each target, and

- a finite coefficient vector `a_1, ..., a_s`.

It is therefore appropriate for compact-support L-statistics. The
user-level DNN wrappers
([`np_dnn()`](https://JakobJuergens.github.io/RdistNN/reference/np_dnn.md),
[`cate_dnn()`](https://JakobJuergens.github.io/RdistNN/reference/cate_dnn.md),
[`plqr_dnn()`](https://JakobJuergens.github.io/RdistNN/reference/plqr_dnn.md))
do **not** call this helper directly for exact DNN weights. Those
wrappers keep the full rank matrix and reevaluate the DNN rank-weight
rule after delete-one updates.

- If `phi` is a numeric **matrix** (n_ref x q) with rows
  `phi[j, ] = phi(j)` and the package was compiled with C++ support,
  this calls the **C++ fast path**.

- Otherwise, it uses a general **R fallback** that also supports
  target-dependent per-neighbor contributions via a function
  `phi(j, t)`.

## Examples

``` r
set.seed(1)
n <- 20
p <- 2
X <- matrix(rnorm(n * p), n, p)
ranks <- cbind(matrix(sample.int(n, n * 3, replace = TRUE), n, 3), sample.int(n, n, replace = TRUE))
a <- c(1, 1, 1)
phi <- cbind(rnorm(n))
out <- jackknife_fast_rank(ranks, a, phi, normalize = TRUE)
str(out$theta0)
#>  num -0.311
```
