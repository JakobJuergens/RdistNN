### Name: summary.cate_result
### Title: Summarize a 'cate_result' object
### Aliases: summary.cate_result

### ** Examples

set.seed(1)
n <- 30
x_mat <- matrix(rnorm(n * 2), n, 2)
out <- cate_dnn(X = x_mat, x0 = c(0, 0), pseudo_outcome = rnorm(n))
summary(out)



