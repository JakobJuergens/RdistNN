### Name: rank_weight_eval
### Title: Evaluate rank-weight for a single target
### Aliases: rank_weight_eval

### ** Examples

rw <- as_rank_weight(function(r) rep(1 / length(r), length(r)))
rank_weight_eval(rw, c(2L, 1L, 3L))



