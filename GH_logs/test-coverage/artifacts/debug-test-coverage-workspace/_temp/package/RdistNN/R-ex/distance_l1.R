### Name: distance_l1
### Title: L1 distance spec (single and batched evaluation)
### Aliases: distance_l1

### ** Examples

set.seed(3)
X <- matrix(rnorm(12), 4, 3)
X0 <- matrix(rnorm(6), 2, 3)
dspec <- distance_l1()
dspec$eval(X, X0[1, ])
dspec$eval_batch(X, X0)




