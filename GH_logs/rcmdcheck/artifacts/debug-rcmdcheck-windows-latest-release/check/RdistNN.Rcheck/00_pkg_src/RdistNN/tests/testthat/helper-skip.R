# Helper predicates for conditional skips used in tests
has_cate_cpp <- function() {
  exists("rustat_cate_ctx_create", mode = "function") &&
    exists("rustat_cate_eval_batch", mode = "function")
}
