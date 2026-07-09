# Minimal IV-GMM demo data with oracle nuisances and test injection hooks.

#' Minimal IV-GMM demo data
#'
#' Generates a small synthetic local IV-GMM dataset for experimenting with
#' oracle doubly-robust IV moments. The output is a thin, test-friendly wrapper
#' around [sim_iv_data()] with shorter design names.
#'
#' @param n Integer sample size.
#' @param p Integer feature dimension, at least 2.
#' @param q Integer number of instruments, between 1 and 3.
#' @param design Character; one of `"homog"`, `"heterog"`, or `"weak"`.
#' @param sigma Noise standard deviation in `Y`.
#' @param clip Clipping bound in `(0, 0.5)` for binary instrument probabilities.
#' @param seed Optional integer seed for reproducibility.
#' @param x0 Optional numeric length-`p` target point.
#' @param inject Optional, one of `NULL`, `"na"`, or `"inf"`; if provided,
#'   injects a single non-finite value into `inject_col` at `idx_inject`.
#' @param inject_col Character injection target.
#' @param idx_inject Integer row to corrupt.
#'
#' @return A list with the same structure as [sim_iv_data()].
#'
#' @examples
#' demo <- iv_demo_data(n = 120, p = 2, q = 2, design = "homog", seed = 1)
#' str(demo$data)
#' demo$theta_x0
#'
#' @export
iv_demo_data <- function(n = 200,
                         p = 2,
                         q = 2,
                         design = c("homog", "heterog", "weak"),
                         sigma = 0.5,
                         clip = 1e-3,
                         seed = NULL,
                         x0 = NULL,
                         inject = NULL,
                         inject_col = "Y",
                         idx_inject = 3L) {
  design <- match.arg(design)
  dgp_name <- switch(design,
    homog = "dgp2_iv_homog_overid",
    heterog = "dgp3_iv_heterog_overid",
    weak = "dgp4_iv_weak_local"
  )

  out <- sim_iv_data(
    design = dgp_name,
    n = n,
    kx = p,
    q = q,
    sigma = sigma,
    clip = clip,
    seed = seed,
    x0 = x0
  )

  if (!is.null(inject)) {
    if (!inject %in% c("na", "inf")) {
      .abort_bad_arg("inject", "must be one of NULL, 'na', or 'inf'", quoted = FALSE)
    }
    idx_inject <- .validate_scalar_int(idx_inject, arg = "idx_inject", min = 1L)
    if (idx_inject > nrow(out$data)) {
      .abort_bad_arg(
        "idx_inject",
        sprintf("must be <= n (got idx_inject=%d, n=%d)", idx_inject, nrow(out$data))
      )
    }

    allowed <- c(
      "Y", "D", paste0("Z", seq_len(out$q)),
      paste0("piZ", seq_len(out$q), "_oracle"),
      "muY_oracle", "muD_oracle", "theta_oracle"
    )
    if (!inject_col %in% allowed) {
      stop(sprintf(
        "Column '%s' (specified by `inject_col`) is not a valid injection target.\nValid options: %s",
        inject_col, paste(allowed, collapse = ", ")
      ), call. = FALSE)
    }
    if (identical(inject, "na")) out$data[[inject_col]][idx_inject] <- NA_real_
    if (identical(inject, "inf")) out$data[[inject_col]][idx_inject] <- Inf
  }

  out
}

#' Ready-made multi-scenario generator for IV-GMM tests
#'
#' @param n Integer number of rows.
#' @param p Integer number of X columns.
#' @param q Integer number of instruments.
#' @return Named list of deterministic scenarios.
#' @keywords internal
#' @noRd
iv_demo_sets <- function(n = 80L, p = 2L, q = 2L) {
  n <- .validate_scalar_int(n, arg = "n", min = 6L)
  p <- .validate_scalar_int(p, arg = "p", min = 2L)
  q <- .validate_scalar_int(q, arg = "q", min = 1L)

  cols <- list(
    X = paste0("X", seq_len(p)),
    Y = "Y",
    D = "D",
    Z = paste0("Z", seq_len(q)),
    pi = paste0("piZ", seq_len(q), "_oracle"),
    muY = "muY_oracle",
    muD = "muD_oracle",
    theta = "theta_oracle"
  )

  make_case <- function(dd) {
    list(
      data = dd$data,
      x0 = dd$x0,
      q = dd$q,
      theta_x0 = dd$theta_x0,
      x_cols = cols$X,
      z_cols = cols$Z,
      cols = cols
    )
  }

  ok_homog <- iv_demo_data(n = n, p = p, q = q, design = "homog", seed = 1)
  ok_heterog <- iv_demo_data(n = n, p = p, q = q, design = "heterog", seed = 2)
  ok_weak <- iv_demo_data(n = n, p = p, q = q, design = "weak", seed = 3)
  nf_y_na <- iv_demo_data(n = n, p = p, q = q, design = "homog", seed = 4,
                          inject = "na", inject_col = "Y")
  nf_z_inf <- iv_demo_data(n = n, p = p, q = q, design = "homog", seed = 5,
                           inject = "inf", inject_col = "Z1")

  list(
    ok_homog = make_case(ok_homog),
    ok_heterog = make_case(ok_heterog),
    ok_weak = make_case(ok_weak),
    nonfinite_y_na = make_case(nf_y_na),
    nonfinite_z_inf = make_case(nf_z_inf)
  )
}
