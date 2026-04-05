# Shared helpers for subset-based score tests
# -------------------------------------------------------------------

`%||%` <- function(x, y) if (is.null(x)) y else x

# Data generators
gen_y <- function(n = 200L, seed = 1L) {
  set.seed(seed)
  rnorm(n)
}
gen_off <- function(n = 200L, seed = 2L) {
  set.seed(seed)
  rnorm(n)
}
gen_d <- function(n = 200L, p = 0.5, seed = 3L) {
  set.seed(seed)
  rbinom(n, 1L, p)
}

# Subset matrix: K rows, size s each (1-based indices)
make_subs <- function(n, s, K = 128L, seed = 123) {
  set.seed(seed)
  matrix(sample.int(n, K * s, replace = TRUE), nrow = K, ncol = s)
}

# Pull g from ctx if stashed there
ctx_g <- function(ctx) {
  if (is.list(ctx) && !is.null(ctx$g)) {
    return(as.numeric(ctx$g))
  }
  NULL
}

# Reference rowwise subset sum in R
subset_row_sums_R <- function(subs, g) {
  rowSums(matrix(g[subs], nrow = nrow(subs), ncol = ncol(subs)))
}

# Tolerant numeric equality
expect_equal_num <- function(x, y, tol = 1e-12) {
  expect_true(isTRUE(all.equal(as.numeric(x), as.numeric(y), tolerance = tol)))
}

# Optional C++ subset-sum availability
has_subset_sum_cpp <- function() exists("rustat_subset_row_sums_g", mode = "function")

# Standard battery for a subset-based score
#   score_factory: function() -> spec with $name, $setup(data), $eval(sub_row, data, ctx, rng),
#                  $eval_batch(subs, data, ctx, rng)
#   mk_data:       function() -> data.frame/list with required named columns
#   s:             subset size
#   corrupt:       function(dat) -> dat with appropriate NA injected (targets correct field)
run_subset_score_battery <- function(score_name, score_factory, mk_data, s, corrupt) {
  test_that(paste0(score_name, ": setup/eval/eval_batch basic contract"), {
    sp <- score_factory()
    dat <- mk_data()
    ctx <- sp$setup(dat)
    n <- nrow(dat) %||% length(dat$y) %||% length(dat$Y)

    subs <- make_subs(n, s, K = 16L)
    s1 <- sp$eval(subs[1L, ], dat, ctx, NULL)
    S <- sp$eval_batch(subs, dat, ctx, NULL)
    expect_equal_num(s1, S[1L])
  })

  test_that(paste0(score_name, ": batched equals loop over eval()"), {
    sp <- score_factory()
    dat <- mk_data()
    ctx <- sp$setup(dat)
    n <- nrow(dat) %||% length(dat$y) %||% length(dat$Y)

    subs <- make_subs(n, s, K = 64L)
    S1 <- sp$eval_batch(subs, dat, ctx, NULL)
    S2 <- vapply(seq_len(nrow(subs)), function(i) sp$eval(subs[i, ], dat, ctx, NULL), numeric(1))
    expect_equal_num(S1, S2)
  })

  test_that(paste0(score_name, ": NA policy errors (targeted)"), {
    sp <- score_factory()
    dat <- mk_data()
    dat_bad <- corrupt(dat)
    expect_error(sp$setup(dat_bad))
  })
}
