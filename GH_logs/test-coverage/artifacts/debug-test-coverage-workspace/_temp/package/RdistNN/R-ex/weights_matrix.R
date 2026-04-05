### Name: weights_matrix
### Title: Compute weight matrix W via an arbitrary ranking + weight
###   mapping
### Aliases: weights_matrix

### ** Examples

set.seed(1)
x_mat <- matrix(rnorm(20), 10, 2)
x0_mat <- matrix(rnorm(6), 3, 2)
rw <- as_rank_weight(function(r) {
  w <- 1 / as.numeric(r)
  w / sum(w)
}, name = "inv_rank")
out <- weights_matrix(x_mat, x0_mat, ranking = "l2", weight = rw, return = c("W", "ranks"))
dim(out$W)



