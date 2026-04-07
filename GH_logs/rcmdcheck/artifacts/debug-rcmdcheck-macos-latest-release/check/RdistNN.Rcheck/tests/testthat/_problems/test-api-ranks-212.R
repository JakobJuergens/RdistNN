# Extracted from test-api-ranks.R:212

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "RdistNN", path = "..")
attach(test_env, warn.conflicts = FALSE)

# prequel ----------------------------------------------------------------------
.rand_ranks_mat <- function(n, m) {
  R <- matrix(NA_integer_, n, m)
  set.seed(123)
  for (j in seq_len(m)) R[, j] <- sample.int(n, n, replace = FALSE)
  R
}

# test -------------------------------------------------------------------------
skip_if_not(exists("rank_weight_eval"), "rank_weight_eval() not found")
skip_if_not(exists("rank_weight_eval_batch"), "rank_weight_eval_batch() not found")
skip_if_not(exists("as_rank_weight"), "as_rank_weight() not found")
rw <- as_rank_weight(function(r, ctx = NULL) {
    w <- 1 / as.numeric(r)
    w / sum(w)
  }, name = "inv_rank")
expect_error(rank_weight_eval(rw, c(1L, 1L, 2L)), "permutation", ignore.case = TRUE)
expect_error(rank_weight_eval(rw, c(1L, NA_integer_, 3L)), "cannot contain NA", ignore.case = TRUE)
expect_error(rank_weight_eval(rw, integer(0)), "positive length", ignore.case = TRUE)
bad_spec <- list(name = "bad")
bad_spec$.__kind <- "rank_weight"
class(bad_spec) <- c("rank_weight", class(bad_spec))
expect_error(rank_weight_eval(bad_spec, c(1L, 2L)), "neither eval nor eval_batch")
bad_len <- as_rank_weight(function(r, ctx = NULL) 1, name = "bad_len")
expect_error(rank_weight_eval(bad_len, c(1L, 2L, 3L)), "length\\(ranks\\)")
bad_inf <- as_rank_weight(function(r, ctx = NULL) rep(Inf, length(r)), name = "bad_inf")
expect_error(rank_weight_eval(bad_inf, c(1L, 2L, 3L)), "finite", ignore.case = TRUE)
rw_bad_batch_shape <- list(
    name = "bad_batch_shape",
    eval_batch = function(R, ctx = NULL) matrix(0, 1, 1)
  )
rw_bad_batch_shape$.__kind <- "rank_weight"
class(rw_bad_batch_shape) <- c("rank_weight", class(rw_bad_batch_shape))
R <- .rand_ranks_mat(5L, 2L)
expect_error(rank_weight_eval_batch(rw_bad_batch_shape, R), "nrow\\(R\\) x ncol\\(R\\)")
