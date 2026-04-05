### Name: score_diff_path
### Title: Build a modular differentiable path for a score
### Aliases: score_diff_path

### ** Examples

sc <- new_score(
  name = "affine_order1",
  order = 1L,
  eval_param = function(idxs, data, theta, ctx, rng) {
    theta[1] * data[idxs] + theta[2]
  }
)
path <- score_diff_path(sc, data = c(1, 2, 4), method = "fd")
path$value(2L, theta = c(0.5, 1))
path$grad(2L, theta = c(0.5, 1))



