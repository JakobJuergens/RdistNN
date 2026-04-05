### Name: new_rank_weight
### Title: Construct a rank-weight spec
### Aliases: new_rank_weight

### ** Examples

rw <- new_rank_weight(
  name = "inv_rank",
  eval = function(ranks, ctx = NULL) {
    w <- 1 / as.numeric(ranks)
    w / sum(w)
  }
)
is_rank_weight(rw)



