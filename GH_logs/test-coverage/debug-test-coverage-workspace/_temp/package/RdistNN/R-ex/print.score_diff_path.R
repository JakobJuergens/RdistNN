### Name: print.score_diff_path
### Title: Print a 'score_diff_path' object
### Aliases: print.score_diff_path

### ** Examples

sc <- new_score(
  name = "order1",
  order = 1L,
  eval_param = function(idxs, data, theta, ctx, rng) theta[1] * data[idxs]
)
print(score_diff_path(sc, data = c(1, 2, 3)))



