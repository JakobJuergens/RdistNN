# tests/testthat/test-wdist-registry.R
# -------------------------------------------------------------------
# Registry for distance-based weights (R/wdist-registry.R)
# Tags:
#   [smoke]   Round-trip register/get
#   [list]    Listing & sorted behavior
#   [overwrite] Re-registering same name overwrites
#   [errors]  Robust error messages
# -------------------------------------------------------------------

# helper constructors
.make_uniform_dw <- function(name) {
  new_dist_weight(
    name = name,
    eval_dx = function(dx, ctx = NULL) {
      n <- length(dx)
      rep(1 / n, n)
    }
  )
}
.make_pointmass_dw <- function(name, pos = 1L) {
  new_dist_weight(
    name = name,
    eval_dx = function(dx, ctx = NULL) {
      n <- length(dx); w <- numeric(n); w[pos] <- 1; w
    }
  )
}

test_that("[smoke] register_wdist() and get_wdist() round-trip works", {
  nm <- "wdist_reg_a"
  dw <- .make_uniform_dw(nm)

  expect_invisible(register_wdist(dw))
  got <- get_wdist(nm)

  expect_s3_class(got, "dist_weight")
  expect_identical(got$name, nm)

  # behavior sanity: returns proper length and unit-sum weights
  dx <- c(3, 2, 5)
  w  <- dist_weight_eval(got, dx)
  expect_length(w, length(dx))
  expect_equal(sum(w), 1, tolerance = 1e-12)
})

test_that("[list] list_wdists() includes registered names and is sorted", {
  nm1 <- "wdist_reg_b"
  nm2 <- "wdist_reg_c"
  register_wdist(.make_uniform_dw(nm1))
  register_wdist(.make_uniform_dw(nm2))

  lst <- list_wdists()
  expect_true(all(c(nm1, nm2) %in% lst))
  expect_identical(lst, sort(lst))  # list is sorted
})

test_that("[overwrite] re-registering the same name overwrites the object", {
  nm <- "wdist_reg_same"
  register_wdist(.make_uniform_dw(nm))           # old
  register_wdist(.make_pointmass_dw(nm, pos=2))  # new with same name

  got <- get_wdist(nm)
  dx  <- c(10, 0.1, 5)
  w   <- dist_weight_eval(got, dx)
  expect_equal(w, c(0, 1, 0))  # reflects the new definition
})

test_that("[errors] bad inputs raise informative errors", {
  # not a dist_weight
  expect_error(register_wdist(123), "is_dist_weight\\(dw\\) is not TRUE")

  # missing name
  expect_error(get_wdist("<<<does_not_exist>>>"),
               "No dist_weight registered under '<<<does_not_exist>>>'.",
               fixed = TRUE)
})
