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
