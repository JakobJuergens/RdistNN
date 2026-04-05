### Name: print.weights_result
### Title: Print a 'weights_result' object
### Aliases: print.weights_result

### ** Examples

set.seed(1)
x_mat <- matrix(rnorm(20), 10, 2)
out <- weights_dnn_exact(x_mat, x0 = c(0, 0), k = 1L, s = 1L)
print(out)



