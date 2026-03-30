### Name: weights_dnn_exact
### Title: DNN weights (exact full-U), distance-based ranking
### Aliases: weights_dnn_exact

### ** Examples

set.seed(1)
x_mat <- matrix(rnorm(20), 10, 2)
x0_mat <- matrix(rnorm(6), 3, 2)
out <- weights_dnn_exact(x_mat, x0_mat, distance = "l2", k = 2L, s = 3L, return = "W")
dim(out$W)



