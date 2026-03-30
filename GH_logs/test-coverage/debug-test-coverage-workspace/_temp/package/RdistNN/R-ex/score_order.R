### Name: score_order
### Title: Resolve kernel order 's'
### Aliases: score_order

### ** Examples

sc <- new_score(
  name = "dyn",
  order = function(n, data, ctx) if (n > 5) 3L else 2L,
  eval = function(idxs, data, ctx, rng) sum(data[idxs])
)
score_order(sc, 4L, 1:4, NULL)



