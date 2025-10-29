# -------------------------------------------------------------------
# Rank registry — register/get/list/overwrite/errors
# Tags:
#   [registry] Basic ops & overwrite
#   [errors]   Missing/invalid inputs
# -------------------------------------------------------------------

test_that("[registry] register/get works and overwrite replaces the object", {
  # minimal valid rank_weights differing in behavior
  rw1 <- new_rank_weight(
    name = "tmp_rank_registry_A",
    eval = function(ranks, ctx = NULL) {
      # uniform weights
      rep(1 / length(ranks), length(ranks))
    }
  )
  nm <- register_rank(rw1)
  # cleanup so we don't leak into other tests
  on.exit(rm(list = nm, envir = RdistNN:::.rank_registry), add = TRUE)

  expect_identical(nm, "tmp_rank_registry_A")

  got <- get_rank("tmp_rank_registry_A")
  expect_s3_class(got, "rank_weight")
  w1 <- rank_weight_eval(got, ranks = 1:5)
  expect_equal(w1, rep(1/5, 5))

  # Overwrite with same name but different rule: all mass on top-1 rank
  rw2 <- new_rank_weight(
    name = "tmp_rank_registry_A",
    eval = function(ranks, ctx = NULL) {
      w <- numeric(length(ranks))
      w[ranks == 1L] <- 1
      w
    }
  )
  register_rank(rw2)

  got2 <- get_rank("tmp_rank_registry_A")
  expect_s3_class(got2, "rank_weight")
  w2 <- rank_weight_eval(got2, ranks = 1:5)
  expect_equal(w2, c(1, 0, 0, 0, 0))
})

test_that("[registry] list_ranks includes the newly registered name and is sorted", {
  rw <- new_rank_weight(
    name = "tmp_rank_registry_B",
    eval = function(ranks, ctx = NULL) rep(1/length(ranks), length(ranks))
  )
  nm <- register_rank(rw)
  on.exit(rm(list = nm, envir = RdistNN:::.rank_registry), add = TRUE)

  lst <- list_ranks()
  expect_true(is.character(lst))
  expect_true("tmp_rank_registry_B" %in% lst)
  # confirm alphabetical sorting
  expect_identical(lst, sort(lst))
})

test_that("[registry] errors on missing name and non-rank inputs", {
  expect_error(get_rank("i_do_not_exist"),
               "No rank_weight registered under 'i_do_not_exist'")
  expect_error(register_rank(123))  # any error is fine
})
