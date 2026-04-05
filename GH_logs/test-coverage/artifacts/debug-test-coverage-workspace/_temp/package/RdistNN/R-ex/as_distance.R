### Name: as_distance
### Title: Wrap a function as a distance spec
### Aliases: as_distance

### ** Examples

spec <- as_distance(
  function(X, x0, ctx = NULL) {
    rowSums(abs(X - matrix(x0, nrow(X), ncol(X), TRUE)))
  }
)
x_mat <- matrix(1:9, nrow = 3)
distance_eval(spec, x_mat, c(0, 0, 0))



