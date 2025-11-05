# tests/testthat/helper-fastpaths.R
.has_row_sums_sym <- function() {
  exists("rustat_subset_row_sums_g", mode = "function") &&
    isTRUE(is.loaded("_RdistNN_rustat_subset_row_sums_g", PACKAGE = "RdistNN"))
}

.has_plqr_ctx_syms <- function() {
  exists("rustat_plqr_ctx_new", mode = "function") &&
    exists("rustat_plqr_g_all", mode = "function") &&
    isTRUE(is.loaded("_RdistNN_rustat_plqr_ctx_new",  PACKAGE = "RdistNN")) &&
    isTRUE(is.loaded("_RdistNN_rustat_plqr_g_all",    PACKAGE = "RdistNN"))
}

skip_if_no_row_sums <- function() testthat::skip_if_not(.has_row_sums_sym(),
                                                        "Compiled helper not available at runtime; skipping parity check."
)
skip_if_no_plqr <- function() testthat::skip_if_not(.has_plqr_ctx_syms(),
                                                    "Compiled PLQR ctx not available at runtime; skipping parity check."
)
