### Name: jackknife_slow_rank
### Title: Slow delete-one jackknife for rank-weighted L-stat rules (fixed
###   ranks)
### Aliases: jackknife_slow_rank

### ** Examples

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



