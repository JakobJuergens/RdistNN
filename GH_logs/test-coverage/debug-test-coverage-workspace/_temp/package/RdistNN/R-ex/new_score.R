### Name: new_score
### Title: Score specification (S3) for U-stat kernels
### Aliases: new_score

### ** Examples

sc <- new_score(
  name = "sum2",
  order = 2L,
  eval = function(idxs, data, ctx, rng) sum(data[idxs])
)
is_score(sc)



