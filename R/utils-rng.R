# Capture before seeding; the public caller registers this closure with on.exit.
# An unseeded call must keep consuming its caller's stream instead.
.rdistnn_rng_restore <- function() {
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  seed <- get0(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  kind <- RNGkind()
  function() {
    if (!identical(RNGkind(), kind)) {
      suppressWarnings(do.call(RNGkind, as.list(kind)))
    }
    if (had_seed) {
      assign(".Random.seed", seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
    invisible(NULL)
  }
}

# Shared subset-drawing contract; public wrappers retain their own defaults.
.rdistnn_demo_subsets <- function(n, K, s, seed) {
  n <- .validate_scalar_int(n, arg = "n", min = 1L)
  s <- .validate_scalar_int(s, arg = "s", min = 1L)
  K <- .validate_scalar_int(K, arg = "K", min = 1L)
  if (s > n) {
    .abort_bad_arg("s", sprintf("must satisfy s <= n (got s=%d, n=%d)", s, n))
  }
  if (!is.null(seed)) {
    restore_rng <- .rdistnn_rng_restore()
    on.exit(restore_rng(), add = TRUE)
    set.seed(seed)
  }
  t(matrix(replicate(K, sort(sample.int(n, s))), nrow = s, ncol = K))
}
