### Name: score_setup
### Title: Set up score context (precomputation)
### Aliases: score_setup

### ** Examples

sc <- new_score(
  name = "ctx_sum",
  order = 1L,
  setup = function(data, ...) list(offset = 1),
  eval = function(idxs, data, ctx, rng) data[idxs] + ctx$offset
)
ctx <- score_setup(sc, 1:3)
score_eval(sc, 2L, 1:3, ctx)



