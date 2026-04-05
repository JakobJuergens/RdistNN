### Name: np_dnn
### Title: Nonparametric regression via DNN rank weights
### Aliases: np_dnn

### ** Examples

set.seed(1)
n <- 60
x_mat <- matrix(rnorm(n * 2), n, 2)
y_vec <- x_mat[, 1] - 0.5 * x_mat[, 2] + rnorm(n, sd = 0.3)
out <- np_dnn(
  X = x_mat, x0 = c(0, 0), Y = y_vec,
  distance = "l2", method = "exact", k = 2L, s = 3L,
  se = FALSE, return = c("W", "ess")
)
out$mu



