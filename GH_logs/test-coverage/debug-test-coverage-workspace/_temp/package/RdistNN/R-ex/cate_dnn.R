### Name: cate_dnn
### Title: CATE via DNN rank weights (orthogonalized local averaging)
### Aliases: cate_dnn

### ** Examples

set.seed(1)
n <- 50
x_mat <- matrix(rnorm(n * 2), n, 2)
x0_mat <- matrix(rnorm(6), 3, 2)
psi <- rnorm(n)
out <- cate_dnn(
  X = x_mat, X0 = x0_mat, pseudo_outcome = psi,
  distance = "l2", method = "exact", k = 2L, s = 3L,
  se = FALSE, return = c("W", "ess")
)
out$tau



