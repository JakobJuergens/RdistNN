### Name: score_eval_batch
### Title: Evaluate kernel on many subsets (R fast-path)
### Aliases: score_eval_batch

### ** Examples

sc <- as_score(function(idxs, data) sum(data[idxs]), s = 2L)
subs <- rbind(c(1L, 2L), c(2L, 3L))
score_eval_batch(sc, subs, c(3, 4, 5), NULL)



