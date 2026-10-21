# Test-only state scope, independent of production RNG restoration.
# Base R does not expose the Box-Muller normal cache.
.test_rng_scope <- function(code, seed = NULL) {
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  saved_seed <- get0(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  saved_kind <- RNGkind()
  on.exit({
    suppressWarnings(do.call(RNGkind, as.list(saved_kind)))
    if (had_seed) {
      assign(".Random.seed", saved_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  })
  if (!is.null(seed)) set.seed(seed)
  force(code)
}
