### Name: weights_dnn
### Title: DNN weights (choose method), distance-based ranking
### Aliases: weights_dnn

### ** Examples

set.seed(1)
x_mat <- matrix(rnorm(20), 10, 2)
x0 <- c(0, 0)
out <- weights_dnn(x_mat, x0 = x0, k = 2L, s = 3L, method = "exact")
dim(out$W)



