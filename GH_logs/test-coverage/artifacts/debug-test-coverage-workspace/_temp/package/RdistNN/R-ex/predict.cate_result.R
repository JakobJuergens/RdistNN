### Name: predict.cate_result
### Title: Predict from a 'cate_result' object at new targets
### Aliases: predict.cate_result

### ** Examples

set.seed(3)
n <- 40
x_mat <- matrix(rnorm(n * 2), n, 2)
fit <- cate_dnn(X = x_mat, x0 = c(0, 0), pseudo_outcome = rnorm(n), se = TRUE)
predict(fit, x0 = c(0.2, -0.1))



