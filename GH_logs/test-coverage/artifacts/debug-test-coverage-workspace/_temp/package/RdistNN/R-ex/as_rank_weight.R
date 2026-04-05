### Name: as_rank_weight
### Title: Wrap a function as a rank_weight spec
### Aliases: as_rank_weight

### ** Examples

rw <- as_rank_weight(function(ranks) {
  w <- 1 / as.numeric(ranks)
  w / sum(w)
})
rank_weight_eval(rw, c(3L, 1L, 2L))



