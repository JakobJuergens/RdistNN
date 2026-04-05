### Name: register_distance
### Title: Register a distance spec under a name
### Aliases: register_distance

### ** Examples

nm <- "tmp_l1_demo"
register_distance(
  nm,
  function(X, x0, ctx = NULL) {
    rowSums(abs(X - matrix(x0, nrow(X), ncol(X), TRUE)))
  }
)
nm %in% distance_names()



