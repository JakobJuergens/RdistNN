# CATE overlap panels (plan panels C-3 and C-4): treated/control weighted
# mass and the local weighted propensity view. Inputs resolve through the
# plan's precedence chain: stored fit components, then user-supplied
# vectors, then an actionable error.

#' Plot local overlap for CATE and CLATE fits
#'
#' Two views of treatment overlap inside the DNN neighbourhoods of a
#' [cate_dnn()] fit. `type = "mass"` shows the weighted treated and control
#' mass per target (for simplex weights the two masses partition one, so
#' the treated mass is the local treated share) with dotted guides and
#' shaded warning bands outside `share_limits`. `type = "propensity"` shows
#' the local weighted propensity distribution: weighted densities per
#' selected target when observation-level propensities are available, and
#' the weighted mean with a q10–q90 whisker summary otherwise.
#'
#' Inputs resolve in a fixed order: components stored on the fit
#' (`t_vec`/`e` from the data-frame interface), then the `treatment`/
#' `propensity` arguments, then an actionable error. Nothing is
#' recomputed from first stages.
#'
#' For [clate_dnn()] fits these same panels describe instrument overlap:
#' the arms are A=1 and A=0 and the propensity is P(A=1|V). Both are retained
#' by the CLATE fit. Treatment uptake D is not used as the instrument arm.
#'
#' @param object A `cate_result` or `clate_result` with the weight matrix `W` stored (the
#'   default).
#' @param type One of `"mass"` or `"propensity"`.
#' @param share_limits Numeric length-2 vector (default `c(0.1, 0.9)`);
#'   treated shares outside this range are flagged as overlap warnings.
#' @param treatment Optional 0/1 vector of length `n`; used when the fit
#'   does not store treatment (matrix interface). For CLATE use `instrument`.
#' @param propensity Optional numeric vector of length `n` in (0, 1); used
#'   when the fit does not store propensities.
#' @param instrument Optional 0/1 instrument vector for a CLATE fit whose
#'   instrument record is unavailable. Stored instrument data take precedence.
#' @param targets Optional integer vector of target ids for the density
#'   view (defaults to the first/median-`n_eff`/last subset above eight
#'   targets).
#' @param ... Ignored; reserved for future options.
#'
#' @return A ggplot object.
#'
#' @examples
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   set.seed(1)
#'   n <- 120
#'   df <- data.frame(X1 = rnorm(n), X2 = rnorm(n))
#'   df$e <- plogis(0.5 * df$X1)
#'   df$T <- as.integer(runif(n) < df$e)
#'   df$m0 <- df$X1
#'   df$m1 <- df$X1 + 1
#'   df$Y <- ifelse(df$T == 1, df$m1, df$m0) + rnorm(n, sd = 0.3)
#'   fit <- cate_dnn(
#'     data = df, x_cols = c("X1", "X2"), y_col = "Y", t_col = "T",
#'     e_col = "e", m0_col = "m0", m1_col = "m1",
#'     X0 = cbind(seq(-1, 1, length.out = 5), 0),
#'     k = 1L, s = 15L, se = FALSE
#'   )
#'   plot_overlap(fit, type = "mass")
#' }
#' @export
plot_overlap <- function(object,
                         type = c("mass", "propensity"),
                         share_limits = c(0.1, 0.9),
                         treatment = NULL,
                         propensity = NULL,
                         targets = NULL,
                         instrument = NULL,
                         ...) {
  .dnn_reuse_scope(object)
  .plot_require_ggplot2("plot_overlap")
  type <- match.arg(type)
  clate <- inherits(object, "clate_result")
  if (!inherits(object, "cate_result") && !clate) {
    .abort_bad_arg(
      "object",
      sprintf(
        "plot_overlap() supports cate_result or clate_result fits; got class '%s'",
        class(object)[1L]
      ),
      quoted = FALSE
    )
  }
  if (clate && !is.null(treatment)) {
    .abort_bad_arg("treatment", "is a CATE argument; use instrument for a CLATE overlap plot")
  }
  if (!is.numeric(share_limits) || length(share_limits) != 2L ||
      any(!is.finite(share_limits)) || share_limits[1L] >= share_limits[2L]) {
    .abort_bad_arg("share_limits", "must be two increasing finite numbers")
  }

  W <- object$W
  if (is.null(W)) {
    .abort_bad_arg(
      "object",
      "must store the weight matrix `W` for overlap panels; refit with `return = c(\"W\")`",
      quoted = FALSE
    )
  }
  W <- as.matrix(W)
  n <- nrow(W)
  m <- ncol(W)
  model <- object$.model

  dd <- .plot_diag_frame(object)
  axis <- .plot_profile_axis(dd)
  honesty <- .plot_honesty_caption(object)

  if (type == "mass") {
    arm_arg <- if (clate) "instrument" else "treatment"
    t_use <- if (clate) {
      if (!is.null(model$iv$A)) model$iv$A else instrument
    } else if (is.list(model) && !is.null(model$t_vec)) {
      model$t_vec
    } else {
      treatment
    }
    if (is.null(t_use)) {
      .abort_bad_arg(
        arm_arg,
        paste0(
          "is required for type = \"mass\": fit through the data ",
          if (clate) "interface with `a_col`, or supply `instrument` (0/1, length n)" else
            "interface with `t_col`, or supply `treatment` (0/1, length n)"
        ),
        quoted = FALSE
      )
    }
    t_num <- .validate_binary01(t_use, arg = arm_arg)
    if (length(t_num) != n) {
      .abort_bad_arg(
        arm_arg,
        sprintf("must have length n = %d", n),
        quoted = FALSE
      )
    }

    treated_mass <- as.numeric(crossprod(t_num, W))
    control_mass <- as.numeric(crossprod(1 - t_num, W))
    total <- treated_mass + control_mass
    treated_share <- ifelse(is.finite(total) & total > 0,
                            treated_mass / total, NA_real_)
    df <- data.frame(
      target_id = rep(seq_len(m), times = 2L),
      arm = rep(if (clate) c("instrument 1", "instrument 0") else c("treated", "control"), each = m),
      mass = c(treated_mass, control_mass),
      stringsAsFactors = FALSE
    )
    if (!is.null(axis$var)) {
      df[[axis$var]] <- rep(dd[[axis$var]], times = 2L)
    }
    n_violations <- sum(
      !is.na(treated_share) &
        (treated_share < share_limits[1L] | treated_share > share_limits[2L])
    )

    xcol <- if (!is.null(axis$var)) axis$var else "target_id"
    p <- ggplot2::ggplot(
      df,
      ggplot2::aes(x = .data[[xcol]], y = .data$mass, colour = .data$arm,
                   fill = .data$arm)
    ) +
      ggplot2::annotate(
        "rect", xmin = -Inf, xmax = Inf,
        ymin = -Inf, ymax = share_limits[1L],
        alpha = 0.08, fill = "grey40"
      ) +
      ggplot2::annotate(
        "rect", xmin = -Inf, xmax = Inf,
        ymin = share_limits[2L], ymax = Inf,
        alpha = 0.08, fill = "grey40"
      ) +
      ggplot2::geom_hline(yintercept = share_limits, linetype = "dotted")
    if (!is.null(axis$var) && m >= 4L) {
      p <- p + ggplot2::geom_line()
    } else {
      p <- p + ggplot2::geom_col(position = "dodge", width = 0.6)
    }
    p <- p +
      ggplot2::scale_colour_manual(
        values = stats::setNames(c(.plot_col_treated, .plot_col_control),
                                if (clate) c("instrument 1", "instrument 0") else c("treated", "control")),
        aesthetics = c("colour", "fill")
      ) +
      ggplot2::labs(
        x = axis$label, y = "Weighted arm mass",
        colour = "arm", fill = "arm"
      )
    captions <- c(
      if (n_violations > 0L) {
        sprintf(
          "%d target(s) with %s share outside [%.2f, %.2f]",
          n_violations, if (clate) "instrument 1" else "treated", share_limits[1L], share_limits[2L]
        )
      },
      honesty
    )
    if (length(captions) > 0L) {
      p <- p + ggplot2::labs(caption = paste(captions, collapse = "\n"))
    }
    return(p + .rdistnn_theme())
  }

  # type == "propensity"
  e_use <- if (clate) {
    if (!is.null(model$iv$propensity$used)) model$iv$propensity$used else propensity
  } else if (is.list(model) && !is.null(model$e)) {
    model$e
  } else {
    propensity
  }
  if (!is.null(e_use)) {
    e_num <- as.numeric(e_use)
    if (length(e_num) != n || any(!is.finite(e_num)) ||
        any(e_num <= 0 | e_num >= 1)) {
      .abort_bad_arg(
        "propensity",
        sprintf("must be a numeric vector of length n = %d with values in (0, 1)", n),
        quoted = FALSE
      )
    }
    neff_by_target <- 1 / colSums(W^2)
    shown <- .plot_select_targets(seq_len(m), neff_by_target, targets)
    dens <- do.call(rbind, lapply(shown, function(j) {
      data.frame(
        target_id = j,
        e = e_num,
        w = W[, j] / sum(W[, j]),
        stringsAsFactors = FALSE
      )
    }))
    dens$target <- factor(
      sprintf("target %d", dens$target_id),
      levels = sprintf("target %d", shown)
    )
    p <- ggplot2::ggplot(dens, ggplot2::aes(x = .data$e)) +
      ggplot2::geom_density(
        ggplot2::aes(weight = .data$w),
        fill = .plot_col_treated, colour = .plot_col_treated, alpha = 0.3
      ) +
      ggplot2::facet_wrap(ggplot2::vars(.data$target), scales = "free_y") +
      ggplot2::coord_cartesian(xlim = c(0, 1)) +
      ggplot2::labs(x = if (clate) "Instrument propensity P(A = 1 | V)" else "Propensity", y = "Weighted density")
  } else if ("propensity_available" %in% names(dd) &&
             any(dd$propensity_available)) {
    xcol <- if (!is.null(axis$var)) axis$var else "target_id"
    p <- ggplot2::ggplot(
      dd,
      ggplot2::aes(x = .data[[xcol]], y = .data$weighted_mean_propensity)
    ) +
      ggplot2::geom_pointrange(
        ggplot2::aes(
          ymin = .data$weighted_propensity_q10,
          ymax = .data$weighted_propensity_q90
        ),
        colour = .plot_col_treated
      ) +
      ggplot2::coord_cartesian(ylim = c(0, 1)) +
      ggplot2::labs(x = axis$label, y = "Weighted propensity (q10-q90)")
  } else {
    .abort_bad_arg(
      "propensity",
      paste0(
        "is required for type = \"propensity\": fit through the data ",
        "interface with `e_col`, or supply `propensity` ",
        "(length n, values in (0, 1))"
      ),
      quoted = FALSE
    )
  }
  if (nzchar(honesty)) {
    p <- p + ggplot2::labs(caption = honesty)
  }
  p + .rdistnn_theme()
}
