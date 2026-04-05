### Name: confint.cate_result
### Title: Confidence intervals for 'cate_result'
### Aliases: confint.cate_result

### ** Examples

set.seed(2)
n <- 25
x_mat <- matrix(rnorm(n * 2), n, 2)
out <- cate_dnn(X = x_mat, x0 = c(0, 0), pseudo_outcome = rnorm(n), se = TRUE)
confint(out)



