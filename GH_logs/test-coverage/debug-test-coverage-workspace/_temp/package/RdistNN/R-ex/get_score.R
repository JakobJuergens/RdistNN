### Name: get_score
### Title: Retrieve a registered score
### Aliases: get_score

### ** Examples

sc0 <- as_score(function(idxs, data) sum(data[idxs]), s = 2L, name = "tmp_score_get")
register_score(sc0)
sc <- get_score("tmp_score_get")
is_score(sc)



