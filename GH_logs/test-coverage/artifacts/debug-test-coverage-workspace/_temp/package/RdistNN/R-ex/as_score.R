### Name: as_score
### Title: Coerce a plain function to a score (simple cases)
### Aliases: as_score

### ** Examples

sc <- as_score(function(idxs, data) prod(data[idxs]), s = 2L, name = "prod2")
score_eval(sc, c(1L, 3L), c(2, 4, 6), NULL)



