# -------------------------------------------------------------------
# Score registry — register/get/list/overwrite/errors
# Tags:
#   [registry] Basic ops & overwrite
#   [errors]   Missing/invalid inputs
#   [order]    Overwrite changes are effective
# -------------------------------------------------------------------

test_that("[registry] register/get works and overwrite replaces the object", {
  # minimal valid scores differing only in order
  sc1 <- new_score(
    name  = "tmp_sc_registry_A",
    order = 2L,
    eval  = function(idxs, data, ctx = NULL, rng = NULL) 1
  )
  nm <- register_score(sc1)
  # cleanup so we don't leak into other tests
  on.exit(rm(list = nm, envir = RdistNN:::.score_registry), add = TRUE)

  expect_identical(nm, "tmp_sc_registry_A")

  got <- get_score("tmp_sc_registry_A")
  expect_s3_class(got, "ustat_score")
  expect_identical(score_order(got, n = 10L, data = NULL), 2L)

  # Overwrite with same name but different order
  sc2 <- new_score(
    name  = "tmp_sc_registry_A",
    order = 3L,
    eval  = function(idxs, data, ctx = NULL, rng = NULL) 2
  )
  register_score(sc2)

  got2 <- get_score("tmp_sc_registry_A")
  expect_s3_class(got2, "ustat_score")
  expect_identical(score_order(got2, n = 10L, data = NULL), 3L)
})

test_that("[registry] list_scores includes the newly registered name and is sorted", {
  sc <- new_score(
    name  = "tmp_sc_registry_B",
    order = 2L,
    eval  = function(idxs, data, ctx = NULL, rng = NULL) 0
  )
  nm <- register_score(sc)
  on.exit(rm(list = nm, envir = RdistNN:::.score_registry), add = TRUE)

  lst <- list_scores()
  expect_true(is.character(lst))
  expect_true("tmp_sc_registry_B" %in% lst)
  # confirm it's alphabetically sorted
  expect_identical(lst, sort(lst))
})

test_that("[registry] errors on missing name and non-score inputs", {
  expect_error(get_score("i_do_not_exist"), "No score registered under 'i_do_not_exist'")
  expect_error(register_score(123)) # any error is fine
})
