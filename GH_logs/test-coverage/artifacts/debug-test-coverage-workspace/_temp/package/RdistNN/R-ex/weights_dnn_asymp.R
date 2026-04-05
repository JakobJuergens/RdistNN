### Name: weights_dnn_asymp
### Title: DNN weights (asymptotic full-U), distance-based ranking
### Aliases: weights_dnn_asymp

### ** Examples

set.seed(1)
x_mat <- matrix(rnorm(20), 10, 2)
x0_mat <- matrix(rnorm(4), 2, 2)
out <- weights_dnn_asymp(x_mat, x0_mat, k = 2L, s = 3L, normalize = TRUE)
dim(out$W)



