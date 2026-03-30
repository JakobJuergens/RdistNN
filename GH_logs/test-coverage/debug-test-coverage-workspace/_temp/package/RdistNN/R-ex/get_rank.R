### Name: get_rank
### Title: Retrieve a registered rank_weight
### Aliases: get_rank

### ** Examples

register_rank("tmp_rank_get", function(ranks) rep(1 / length(ranks), length(ranks)))
is_rank_weight(get_rank("tmp_rank_get"))



