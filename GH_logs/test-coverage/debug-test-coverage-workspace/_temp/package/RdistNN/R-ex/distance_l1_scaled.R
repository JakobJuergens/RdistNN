### Name: distance_l1_scaled
### Title: Scaled L1 distance spec (one-time setup + eval/batch)
### Aliases: distance_l1_scaled

### ** Examples

X <- matrix(rnorm(12), 4, 3)
X0 <- matrix(rnorm(6), 2, 3)
dspec <- distance_l1_scaled(rule = "sd")
ctx <- dspec$setup(X)
dspec$eval_batch(X, X0, ctx)



