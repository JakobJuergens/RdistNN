# Generic rank-order delete-one jackknife (prefix/suffix engine)

Computes all \\n\\ leave-one-out estimates for a rank-weighted linear
estimator \\\hat\theta_t = \sum_i w_n\[r\_{it}\]\\\phi_i\\ at each of
\\m\\ target points, using an \\O(n \cdot m \cdot q)\\ prefix/suffix
algorithm.

## Usage

``` r
jackknife_rank_full(ranks, PHI, wn, wnm1)
```

## Arguments

- ranks:

  Integer matrix (\\n \times m\\); column \\t\\ is a permutation of
  `1:n` giving the rank of each observation for target \\t\\ (rank 1 =
  nearest).

- PHI:

  Finite numeric matrix (\\n \times q\\); row \\i\\ is the
  \\q\\-dimensional signal contribution of observation \\i\\. For a
  scalar signal, pass `cbind(y)`.

- wn:

  Finite numeric vector of length \\n\\; rank-weight table at full
  sample size \\n\\. Entry `wn[r]` is the weight assigned to rank \\r\\.

- wnm1:

  Finite numeric vector of length \\n-1\\; rank-weight table at
  leave-one-out sample size \\n-1\\.

## Value

A list with:

- `theta0`:

  Numeric matrix (\\q \times m\\); baseline estimates at full sample
  size.

- `theta_loo`:

  Array with `dim = c(n, q, m)`; element `[i, qq, t]` is the
  leave-one-out estimate after deleting observation \\i\\ for signal
  dimension \\qq\\ and target \\t\\.

- `info`:

  List with fields `n`, `m`, `q`.

## Details

The rank-weight rule is supplied via precomputed tables `wn` and `wnm1`
rather than a spec object, so the function is completely agnostic to the
weight construction. Callers are responsible for providing consistent
tables:

- `wn` should be the rank-weight sequence evaluated at sample size \\n\\
  (i.e. for the full-sample estimator).

- `wnm1` should be the **same rule** evaluated at sample size \\n-1\\
  (i.e. for each leave-one-out sub-sample).

When the package was compiled with C++ support, the C++ engine is used.
Otherwise the function falls back silently to a pure-R implementation of
the same algorithm.

## See also

[`jackknife_fast_rank()`](https://JakobJuergens.github.io/RdistNN/reference/jackknife_fast_rank.md)
for the finite-support L-stat engine;
[`jackknife_slow_rank()`](https://JakobJuergens.github.io/RdistNN/reference/jackknife_slow_rank.md)
for the reference implementation.

## Examples

``` r
set.seed(1)
n <- 20L; m <- 3L; q <- 2L
ranks <- apply(matrix(runif(n * m), n, m), 2L, rank)
storage.mode(ranks) <- "integer"
PHI <- matrix(rnorm(n * q), n, q)
# toy uniform rank weights
wn   <- rep(1 / n, n)
wnm1 <- rep(1 / (n - 1L), n - 1L)
out <- jackknife_rank_full(ranks, PHI, wn, wnm1)
dim(out$theta_loo)   # c(n, q, m)
#> [1] 20  2  3
dim(out$theta0)      # c(q, m)
#> [1] 2 3
```
