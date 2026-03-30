### Name: distance_l2
### Title: L2 distance spec (single and batched evaluation)
### Aliases: distance_l2

### ** Examples

X <- matrix(rnorm(12), 4, 3)
X0 <- matrix(rnorm(6), 2, 3)
dspec <- distance_l2()
dspec$eval_batch(X, X0)



