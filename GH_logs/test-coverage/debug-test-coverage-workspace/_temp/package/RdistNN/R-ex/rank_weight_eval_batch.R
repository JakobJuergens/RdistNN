### Name: rank_weight_eval_batch
### Title: Evaluate rank-weight for many targets (batch)
### Aliases: rank_weight_eval_batch

### ** Examples

rw <- as_rank_weight(function(r) rep(1 / length(r), length(r)))
rank_matrix <- cbind(c(1L, 2L, 3L), c(2L, 1L, 3L))
rank_weight_eval_batch(rw, rank_matrix)



