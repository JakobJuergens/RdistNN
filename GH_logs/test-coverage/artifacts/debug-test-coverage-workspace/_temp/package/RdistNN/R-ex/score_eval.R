### Name: score_eval
### Title: Evaluate kernel on one subset
### Aliases: score_eval

### ** Examples

sc <- as_score(function(idxs, data) sum(data[idxs]), s = 2L)
score_eval(sc, c(1L, 2L), c(3, 4, 5), NULL)



