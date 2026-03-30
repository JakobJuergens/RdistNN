### Name: register_rank
### Title: Register a rank_weight under a name
### Aliases: register_rank

### ** Examples

nm <- "tmp_rank_demo"
register_rank(nm, function(ranks) rep(1 / length(ranks), length(ranks)))
nm %in% rank_names()



