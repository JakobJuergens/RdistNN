# inst/paper-replication/iv/grid.R
#
# Shared grid helpers for the IV-GMM simulation replication scripts.
# These helpers are intentionally script-local infrastructure, not package API.

iv_paper_grid <- function(smoke = identical(Sys.getenv("RDISTNN_IV_SMOKE"), "1")) {
  if (isTRUE(smoke)) {
    return(list(
      designs = c("dgp1_iv_justid_linear", "dgp2_iv_homog_overid"),
      n_vals = 80L,
      kx_vals = 2L,
      q_vals = 2L,
      n_reps = 2L,
      master_seed = 42L,
      k = 1L,
      s_fn = function(n, kx) floor(n^((kx + 1L) / (kx + 2L))),
      methods = "iv_oracle",
      alpha = 0.05,
      sigma = 0.5,
      clip = 1e-3,
      n_folds = 2L
    ))
  }

  list(
    designs = c(
      "dgp1_iv_justid_linear",
      "dgp2_iv_homog_overid",
      "dgp3_iv_heterog_overid",
      "dgp5_iv_correlated_instruments",
      "dgp6_iv_nonlinear_firststage"
    ),
    n_vals = c(250L, 500L, 1000L, 2000L),
    kx_vals = 2L,
    q_vals = c(2L, 3L),
    n_reps = 10000L,
    master_seed = 42L,
    k = 1L,
    s_fn = function(n, kx) floor(n^((kx + 1L) / (kx + 2L))),
    methods = c(
      "dnn_iv_fullsample",
      "dnn_iv_crossfit",
      "local_2sls",
      "2sls_interact"
    ),
    alpha = 0.05,
    sigma = 0.5,
    clip = 1e-3,
    n_folds = 5L
  )
}

iv_diagnostic_grid <- function(smoke = identical(Sys.getenv("RDISTNN_IV_SMOKE"), "1")) {
  grid <- iv_paper_grid(smoke = smoke)
  grid$designs <- "dgp4_iv_weak_local"
  grid
}

iv_paper_cells <- function(grid = iv_paper_grid()) {
  rows <- list()
  idx <- 1L

  for (design in grid$designs) {
    q_for_design <- if (identical(design, "dgp1_iv_justid_linear")) {
      1L
    } else {
      unique(as.integer(grid$q_vals))
    }

    for (n in grid$n_vals) {
      for (kx in grid$kx_vals) {
        for (q in q_for_design) {
          rows[[idx]] <- data.frame(
            cell = idx,
            design = design,
            n = as.integer(n),
            kx = as.integer(kx),
            q = as.integer(q),
            stringsAsFactors = FALSE
          )
          idx <- idx + 1L
        }
      }
    }
  }

  if (length(rows) == 0L) {
    return(data.frame(
      cell = integer(),
      design = character(),
      n = integer(),
      kx = integer(),
      q = integer(),
      stringsAsFactors = FALSE
    ))
  }

  do.call(rbind, rows)
}

iv_expected_cells <- function(grid = iv_paper_grid()) {
  nrow(iv_paper_cells(grid))
}

IV_PAPER_GRID <- iv_paper_grid(smoke = FALSE)
IV_PAPER_CELL_COUNT <- iv_expected_cells(IV_PAPER_GRID)
