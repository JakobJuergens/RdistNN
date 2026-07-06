# Dashboard composition (plan panels C-6 and the np analogue): patchwork
# layouts of the individual panels with a metadata annotation. patchwork is
# a suggested dependency; when it is missing the dashboard degrades to a
# named list of ggplots with a message.

.plot_dashboard_panels <- function(specs) {
  panels <- list()
  for (nm in names(specs)) {
    p <- tryCatch(specs[[nm]](), error = function(e) NULL)
    if (!is.null(p)) panels[[nm]] <- p
  }
  panels
}

.plot_dashboard_compose <- function(panels, object) {
  if (length(panels) == 0L) {
    .abort_bad_state(
      "no dashboard panel could be built for this fit; see the individual plot_* builders for the missing inputs."
    )
  }
  if (!.plot_pkg_installed("patchwork")) {
    message(
      "patchwork is not installed; returning the dashboard panels as a ",
      "named list of ggplot objects. install.packages(\"patchwork\") to ",
      "compose them."
    )
    return(panels)
  }
  patchwork::wrap_plots(panels, ncol = 2L) +
    patchwork::plot_annotation(
      caption = .plot_dashboard_annotation(object)
    )
}

#' Compose a diagnostics dashboard for a fitted object
#'
#' A patchwork composition of the estimate panel with the weight and
#' overlap diagnostics, annotated with the fit's metadata (`n`, `m`, `k`,
#' `s` with provenance, distance, and jackknife type). CATE dashboards use
#' estimates, overlap mass, local propensity, and cumulative weight share;
#' nonparametric-regression dashboards use estimates, the diagnostics
#' facets, and the weight decay and concentration views. Panels whose
#' inputs are unavailable are substituted (the diagnostics facets take the
#' place of a missing overlap panel) or dropped rather than failing the
#' dashboard.
#'
#' When patchwork is not installed the panels are returned as a named list
#' of ggplot objects with a message.
#'
#' @param object A fitted result object (`cate_result` or `np_result`).
#' @param ... Passed on to [plot_estimates()] for the estimates panel.
#'
#' @return A patchwork object (or a named list of ggplots when patchwork
#'   is unavailable).
#'
#' @examples
#' if (requireNamespace("ggplot2", quietly = TRUE) &&
#'     requireNamespace("patchwork", quietly = TRUE)) {
#'   set.seed(1)
#'   n <- 100
#'   x_mat <- matrix(rnorm(n * 2), n, 2)
#'   y_vec <- x_mat[, 1] + rnorm(n, sd = 0.3)
#'   fit <- np_dnn(
#'     X = x_mat, X0 = cbind(seq(-1, 1, length.out = 6), 0), Y = y_vec,
#'     k = 1L, s = 12L, return = c("W", "dx")
#'   )
#'   plot_dashboard(fit)
#' }
#' @export
plot_dashboard <- function(object, ...) {
  UseMethod("plot_dashboard")
}

#' @rdname plot_dashboard
#' @export
plot_dashboard.default <- function(object, ...) {
  .abort_bad_arg(
    "object",
    sprintf(
      "plot_dashboard() is not implemented for class '%s'",
      class(object)[1L]
    ),
    quoted = FALSE
  )
}

#' @rdname plot_dashboard
#' @export
plot_dashboard.cate_result <- function(object, ...) {
  .plot_require_ggplot2("plot_dashboard")
  specs <- list(
    estimates = function() plot_estimates(object, ...),
    overlap = function() plot_overlap(object, type = "mass"),
    propensity = function() plot_overlap(object, type = "propensity"),
    weights = function() plot_weights(object, type = "cumulative")
  )
  if (!is.null(object$tuning_grid)) {
    specs$tuning <- function() plot_tuning(object)
  }
  panels <- .plot_dashboard_panels(specs)
  # Availability substitution (plan C-6): a missing overlap panel is
  # replaced by the diagnostics facets rather than an empty slot.
  if (!"overlap" %in% names(panels)) {
    diag_panel <- tryCatch(
      ggplot2::autoplot(diagnostics(object)),
      error = function(e) NULL
    )
    if (!is.null(diag_panel)) {
      panels <- append(panels, list(diagnostics = diag_panel), after = 1L)
    }
  }
  .plot_dashboard_compose(panels, object)
}

#' @rdname plot_dashboard
#' @export
plot_dashboard.np_result <- function(object, ...) {
  .plot_require_ggplot2("plot_dashboard")
  specs <- list(
    estimates = function() plot_estimates(object, ...),
    diagnostics = function() ggplot2::autoplot(diagnostics(object)),
    weights = function() plot_weights(object, type = "profile"),
    concentration = function() plot_weights(object, type = "cumulative")
  )
  if (!is.null(object$tuning_grid)) {
    specs$tuning <- function() plot_tuning(object)
  }
  panels <- .plot_dashboard_panels(specs)
  .plot_dashboard_compose(panels, object)
}
