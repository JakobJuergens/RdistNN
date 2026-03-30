### Name: order_by_l1_scaled
### Title: Order rows by scaled L1 (stable ties)
### Aliases: order_by_l1_scaled

### ** Examples

set.seed(6)
X <- matrix(rnorm(10), 5, 2)
x0 <- c(0, 0)
order_by_l1_scaled(X, x0, rule = "sd")




