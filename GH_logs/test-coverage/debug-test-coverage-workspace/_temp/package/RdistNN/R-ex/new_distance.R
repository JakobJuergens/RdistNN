### Name: new_distance
### Title: Construct a distance spec
### Aliases: new_distance

### ** Examples

spec <- new_distance(
  name = "sq_l2",
  eval = function(X, x0, ctx = NULL) rowSums((X - matrix(x0, nrow(X), ncol(X), TRUE))^2)
)
is_distance(spec)



