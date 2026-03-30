### Name: jackknife_fast_rank
### Title: Fast delete-one jackknife for rank-weighted functionals (auto
###   C++ when available, R fallback otherwise)
### Aliases: jackknife_fast_rank

### ** Examples

set.seed(1)
n <- 20
p <- 2
X <- matrix(rnorm(n * p), n, p)
ranks <- cbind(matrix(sample.int(n, n * 3, replace = TRUE), n, 3), sample.int(n, n, replace = TRUE))
a <- c(1, 1, 1)
phi <- cbind(rnorm(n))
out <- jackknife_fast_rank(ranks, a, phi, normalize = TRUE)
str(out$theta0)



