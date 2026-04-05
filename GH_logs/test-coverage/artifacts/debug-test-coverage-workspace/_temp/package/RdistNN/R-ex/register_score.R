### Name: register_score
### Title: Register a score
### Aliases: register_score

### ** Examples

sc <- as_score(function(idxs, data) sum(data[idxs]), s = 2L, name = "tmp_score_demo")
register_score(sc)
"tmp_score_demo" %in% list_scores()



