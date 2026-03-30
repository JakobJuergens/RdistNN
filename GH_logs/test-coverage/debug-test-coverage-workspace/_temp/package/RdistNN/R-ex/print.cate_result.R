### Name: print.cate_result
### Title: Print a 'cate_result' object
### Aliases: print.cate_result

### ** Examples

set.seed(1)
n <- 20
x_mat <- matrix(rnorm(n * 2), n, 2)
out <- cate_dnn(
  X = x_mat,
  x0 = c(0, 0),
  pseudo_outcome = rnorm(n),
  se = FALSE
)
print(out)



