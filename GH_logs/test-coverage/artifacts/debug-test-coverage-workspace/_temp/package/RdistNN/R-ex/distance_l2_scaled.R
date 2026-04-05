### Name: distance_l2_scaled
### Title: Scaled L2 distance spec (one-time setup + eval/batch)
### Aliases: distance_l2_scaled

### ** Examples

X <- matrix(rnorm(12), 4, 3)
X0 <- matrix(rnorm(6), 2, 3)
dspec <- distance_l2_scaled(rule = "sd")
ctx <- dspec$setup(X)
dspec$eval_batch(X, X0, ctx)



