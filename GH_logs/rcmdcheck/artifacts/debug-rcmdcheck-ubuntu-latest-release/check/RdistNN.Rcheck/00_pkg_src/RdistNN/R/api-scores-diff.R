# Optional derivative paths for score objects.

.score_check_theta <- function(theta, arg = "theta") {
  if (!is.numeric(theta) || length(theta) < 1L) {
    .abort_bad_arg(arg, "must be a non-empty numeric vector")
  }
  theta <- as.numeric(theta)
  if (any(!is.finite(theta))) {
    .abort_bad_arg(arg, "must be finite")
  }
  theta
}

.score_fd_step_vec <- function(fd_step, theta, arg = "fd_step") {
  step <- as.numeric(fd_step)
  if (length(step) == 1L) {
    if (!is.finite(step) || step <= 0) {
      .abort_bad_arg(arg, "must be positive and finite")
    }
    return(step * pmax(1, abs(theta)))
  }
  if (length(step) != length(theta) || any(!is.finite(step)) || any(step <= 0)) {
    .abort_bad_state("`fd_step` must be a positive finite scalar or length(theta) vector.")
  }
  step
}

.score_validate_scalar_numeric <- function(x, label) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x)) {
    stop(label, " must return a finite scalar numeric.", call. = FALSE)
  }
  as.numeric(x)
}

.score_validate_batch_numeric <- function(x, K, label) {
  if (!is.numeric(x) || length(x) != K || any(!is.finite(x))) {
    stop(label, " must return numeric(K) with finite entries.", call. = FALSE)
  }
  as.numeric(x)
}

.score_validate_grad_scalar <- function(g, p, label) {
  if (!is.numeric(g) || length(g) != p || any(!is.finite(g))) {
    stop(label, " must return finite numeric length(theta).", call. = FALSE)
  }
  as.numeric(g)
}

.score_validate_grad_batch <- function(G, K, p, label) {
  if (is.null(dim(G))) {
    if (p == 1L && is.numeric(G) && length(G) == K && all(is.finite(G))) {
      G <- matrix(as.numeric(G), nrow = K, ncol = 1L)
    } else {
      stop(label, " must return a K x length(theta) finite matrix.", call. = FALSE)
    }
  }
  if (!is.numeric(G) || nrow(G) != K || ncol(G) != p || any(!is.finite(G))) {
    stop(label, " must return a K x length(theta) finite matrix.", call. = FALSE)
  }
  matrix(as.numeric(G), nrow = K, ncol = p)
}

.score_fun_supports_theta <- function(fun) {
  fmls <- tryCatch(names(formals(fun)), error = function(e) character())
  any(c("theta", "param", "...") %in% fmls)
}

.score_assert_theta_support <- function(fun, label) {
  if (!.score_fun_supports_theta(fun)) {
    stop(
      label, " does not declare `theta`/`param` (or `...`); ",
      "cannot use parameter-aware path.",
      call. = FALSE
    )
  }
}

# Internal low-level helper; user-facing entry is score_diff_path().
score_eval_param <- function(score, idxs, data, theta, ctx = NULL, rng = NULL) {
  if (!is_score(score)) .abort_bad_arg("score", "must be a `ustat_score`")
  n <- .score_n(data)
  s <- score_order(score, n, data, ctx)
  idxs <- .score_check_idxs(idxs, n, s)
  theta <- .score_check_theta(theta)

  args_scalar <- list(idxs = idxs, data = data, theta = theta, param = theta, ctx = ctx, rng = rng)
  args_batch1 <- list(
    subs = matrix(as.integer(idxs), nrow = 1L),
    data = data, theta = theta, param = theta, ctx = ctx, rng = rng
  )

  if (is.function(score$eval_param)) {
    .score_assert_theta_support(score$eval_param, "score$eval_param")
    out <- .call_with_supported_args(score$eval_param, args_scalar)
    return(.score_validate_scalar_numeric(out, "score$eval_param"))
  }
  if (is.function(score$eval_param_batch)) {
    .score_assert_theta_support(score$eval_param_batch, "score$eval_param_batch")
    out <- .call_with_supported_args(score$eval_param_batch, args_batch1)
    return(.score_validate_scalar_numeric(out, "score$eval_param_batch"))
  }
  if (is.function(score$eval)) {
    .score_assert_theta_support(score$eval, "score$eval")
    out <- .call_with_supported_args(score$eval, args_scalar)
    return(.score_validate_scalar_numeric(out, "score$eval"))
  }
  if (is.function(score$eval_batch)) {
    .score_assert_theta_support(score$eval_batch, "score$eval_batch")
    out <- .call_with_supported_args(score$eval_batch, args_batch1)
    return(.score_validate_scalar_numeric(out, "score$eval_batch"))
  }

  stop(
    "No parameter-aware evaluator available. Provide `eval_param`/`eval_param_batch` ",
    "or an `eval`/`eval_batch` that accepts `theta`.",
    call. = FALSE
  )
}

# Internal low-level helper; user-facing entry is score_diff_path().
score_eval_param_batch <- function(score, subs, data, theta, ctx = NULL, rng = NULL) {
  if (!is_score(score)) .abort_bad_arg("score", "must be a `ustat_score`")
  n <- .score_n(data)
  s <- score_order(score, n, data, ctx)
  subs <- .score_check_subs(subs, n, s)
  K <- nrow(subs)
  theta <- .score_check_theta(theta)

  args_batch <- list(subs = subs, data = data, theta = theta, param = theta, ctx = ctx, rng = rng)

  if (is.function(score$eval_param_batch)) {
    .score_assert_theta_support(score$eval_param_batch, "score$eval_param_batch")
    out <- .call_with_supported_args(score$eval_param_batch, args_batch)
    return(.score_validate_batch_numeric(out, K, "score$eval_param_batch"))
  }
  if (is.function(score$eval_param)) {
    .score_assert_theta_support(score$eval_param, "score$eval_param")
    out <- numeric(K)
    for (k in seq_len(K)) {
      out[k] <- .call_with_supported_args(
        score$eval_param,
        list(idxs = subs[k, ], data = data, theta = theta, param = theta, ctx = ctx, rng = rng)
      )
    }
    return(.score_validate_batch_numeric(out, K, "score$eval_param"))
  }
  if (is.function(score$eval_batch)) {
    .score_assert_theta_support(score$eval_batch, "score$eval_batch")
    out <- .call_with_supported_args(score$eval_batch, args_batch)
    return(.score_validate_batch_numeric(out, K, "score$eval_batch"))
  }
  if (is.function(score$eval)) {
    .score_assert_theta_support(score$eval, "score$eval")
    out <- numeric(K)
    for (k in seq_len(K)) {
      out[k] <- .call_with_supported_args(
        score$eval,
        list(idxs = subs[k, ], data = data, theta = theta, param = theta, ctx = ctx, rng = rng)
      )
    }
    return(.score_validate_batch_numeric(out, K, "score$eval"))
  }

  stop(
    "No parameter-aware evaluator available. Provide `eval_param`/`eval_param_batch` ",
    "or an `eval`/`eval_batch` that accepts `theta`.",
    call. = FALSE
  )
}

# Internal low-level helper; user-facing entry is score_diff_path().
score_grad <- function(
  score, idxs, data, theta, ctx = NULL, rng = NULL,
  method = c("auto", "analytic", "fd"),
  fd_step = NULL, fd_scheme = NULL
) {
  if (!is_score(score)) .abort_bad_arg("score", "must be a `ustat_score`")
  n <- .score_n(data)
  s <- score_order(score, n, data, ctx)
  idxs <- .score_check_idxs(idxs, n, s)
  theta <- .score_check_theta(theta)
  p <- length(theta)
  method <- match.arg(method)

  args_scalar <- list(idxs = idxs, data = data, theta = theta, param = theta, ctx = ctx, rng = rng)
  args_batch1 <- list(
    subs = matrix(as.integer(idxs), nrow = 1L),
    data = data, theta = theta, param = theta, ctx = ctx, rng = rng
  )

  if (method %in% c("auto", "analytic")) {
    if (is.function(score$grad)) {
      g <- .call_with_supported_args(score$grad, args_scalar)
      return(.score_validate_grad_scalar(g, p, "score$grad"))
    }
    if (is.function(score$grad_batch)) {
      g <- .call_with_supported_args(score$grad_batch, args_batch1)
      G <- .score_validate_grad_batch(g, K = 1L, p = p, label = "score$grad_batch")
      return(as.numeric(G[1L, ]))
    }
    if (method == "analytic") {
      stop("No analytic derivative provided (`grad` or `grad_batch`).", call. = FALSE)
    }
  }

  scheme <- if (is.null(fd_scheme)) score$fd_scheme else match.arg(fd_scheme, c("central", "forward"))
  base_step <- if (is.null(fd_step)) score$fd_step else fd_step
  h <- .score_fd_step_vec(base_step, theta)

  g <- numeric(p)
  if (scheme == "forward") {
    f0 <- score_eval_param(score, idxs, data, theta, ctx = ctx, rng = rng)
  }
  for (j in seq_len(p)) {
    th_plus <- theta
    th_plus[j] <- th_plus[j] + h[j]
    f_plus <- score_eval_param(score, idxs, data, th_plus, ctx = ctx, rng = rng)
    if (scheme == "central") {
      th_minus <- theta
      th_minus[j] <- th_minus[j] - h[j]
      f_minus <- score_eval_param(score, idxs, data, th_minus, ctx = ctx, rng = rng)
      g[j] <- (f_plus - f_minus) / (2 * h[j])
    } else {
      g[j] <- (f_plus - f0) / h[j]
    }
  }

  .score_validate_grad_scalar(g, p, "Finite-difference gradient")
}

# Internal low-level helper; user-facing entry is score_diff_path().
score_grad_batch <- function(
  score, subs, data, theta, ctx = NULL, rng = NULL,
  method = c("auto", "analytic", "fd"),
  fd_step = NULL, fd_scheme = NULL
) {
  if (!is_score(score)) .abort_bad_arg("score", "must be a `ustat_score`")
  n <- .score_n(data)
  s <- score_order(score, n, data, ctx)
  subs <- .score_check_subs(subs, n, s)
  K <- nrow(subs)
  theta <- .score_check_theta(theta)
  p <- length(theta)
  method <- match.arg(method)

  args_batch <- list(subs = subs, data = data, theta = theta, param = theta, ctx = ctx, rng = rng)

  if (method %in% c("auto", "analytic")) {
    if (is.function(score$grad_batch)) {
      G <- .call_with_supported_args(score$grad_batch, args_batch)
      return(.score_validate_grad_batch(G, K, p, "score$grad_batch"))
    }
    if (is.function(score$grad)) {
      G <- matrix(NA_real_, nrow = K, ncol = p)
      for (k in seq_len(K)) {
        gk <- .call_with_supported_args(
          score$grad,
          list(idxs = subs[k, ], data = data, theta = theta, param = theta, ctx = ctx, rng = rng)
        )
        G[k, ] <- .score_validate_grad_scalar(gk, p, "score$grad")
      }
      return(G)
    }
    if (method == "analytic") {
      stop("No analytic derivative provided (`grad` or `grad_batch`).", call. = FALSE)
    }
  }

  scheme <- if (is.null(fd_scheme)) score$fd_scheme else match.arg(fd_scheme, c("central", "forward"))
  base_step <- if (is.null(fd_step)) score$fd_step else fd_step
  h <- .score_fd_step_vec(base_step, theta)

  G <- matrix(NA_real_, nrow = K, ncol = p)
  if (scheme == "forward") {
    f0 <- score_eval_param_batch(score, subs, data, theta, ctx = ctx, rng = rng)
  }
  for (j in seq_len(p)) {
    th_plus <- theta
    th_plus[j] <- th_plus[j] + h[j]
    f_plus <- score_eval_param_batch(score, subs, data, th_plus, ctx = ctx, rng = rng)
    if (scheme == "central") {
      th_minus <- theta
      th_minus[j] <- th_minus[j] - h[j]
      f_minus <- score_eval_param_batch(score, subs, data, th_minus, ctx = ctx, rng = rng)
      G[, j] <- (f_plus - f_minus) / (2 * h[j])
    } else {
      G[, j] <- (f_plus - f0) / h[j]
    }
  }

  .score_validate_grad_batch(G, K, p, "Finite-difference gradient")
}

#' Build a modular differentiable path for a score
#'
#' Returns closures for value and gradient evaluation of a score at any
#' parameter vector `theta`, using analytic or finite-difference derivatives.
#'
#' @param score A `ustat_score`.
#' @param data Data object passed to score evaluators.
#' @param ctx Optional context from [score_setup()]. If `NULL`, one is
#'   created via [score_setup()].
#' @param method One of `c("auto", "analytic", "fd")` for gradients.
#' @param fd_step Optional finite-difference step override.
#' @param fd_scheme Optional finite-difference scheme override.
#'
#' @return A list with closures:
#'   - `value(idxs, theta, rng = NULL)` for one subset;
#'   - `value_batch(subs, theta, rng = NULL)` for many subsets;
#'   - `grad(idxs, theta, rng = NULL)` for one-subset gradients;
#'   - `grad_batch(subs, theta, rng = NULL)` for batch gradients;
#'   - `meta` with derivative settings.
#'
#' @seealso [new_score()], [score_setup()], [score_eval()], [score_eval_batch()].
#'
#' @examples
#' sc <- new_score(
#'   name = "affine_order1",
#'   order = 1L,
#'   eval_param = function(idxs, data, theta, ctx, rng) {
#'     theta[1] * data[idxs] + theta[2]
#'   }
#' )
#' path <- score_diff_path(sc, data = c(1, 2, 4), method = "fd")
#' path$value(2L, theta = c(0.5, 1))
#' path$grad(2L, theta = c(0.5, 1))
#' @export
score_diff_path <- function(
  score, data, ctx = NULL,
  method = c("auto", "analytic", "fd"),
  fd_step = NULL, fd_scheme = NULL
) {
  if (!is_score(score)) {
    .abort_bad_arg("score", "must be a `ustat_score`")
  }
  method <- match.arg(method)
  if (is.null(ctx)) {
    ctx <- score_setup(score, data)
  }

  out <- list(
    value = function(idxs, theta, rng = NULL) {
      score_eval_param(score, idxs, data, theta, ctx = ctx, rng = rng)
    },
    value_batch = function(subs, theta, rng = NULL) {
      score_eval_param_batch(score, subs, data, theta, ctx = ctx, rng = rng)
    },
    grad = function(idxs, theta, rng = NULL) {
      score_grad(
        score, idxs, data, theta, ctx = ctx, rng = rng,
        method = method, fd_step = fd_step, fd_scheme = fd_scheme
      )
    },
    grad_batch = function(subs, theta, rng = NULL) {
      score_grad_batch(
        score, subs, data, theta, ctx = ctx, rng = rng,
        method = method, fd_step = fd_step, fd_scheme = fd_scheme
      )
    },
    meta = list(
      score_name = score$name,
      method = method,
      fd_step = if (is.null(fd_step)) score$fd_step else fd_step,
      fd_scheme = if (is.null(fd_scheme)) score$fd_scheme else fd_scheme
    )
  )
  class(out) <- "score_diff_path"
  out
}

#' Print a `score_diff_path` object
#'
#' @param x A `score_diff_path` object created by [score_diff_path()].
#' @param ... Not used.
#' @return `x`, invisibly.
#' @examples
#' sc <- new_score(
#'   name = "order1",
#'   order = 1L,
#'   eval_param = function(idxs, data, theta, ctx, rng) theta[1] * data[idxs]
#' )
#' print(score_diff_path(sc, data = c(1, 2, 3)))
#' @export
print.score_diff_path <- function(x, ...) {
  cat("<score_diff_path>\n", sep = "")
  cat("  score : ", x$meta$score_name, "\n", sep = "")
  cat("  grad  : ", x$meta$method, "\n", sep = "")
  cat("  fd    : step=", paste0(x$meta$fd_step, collapse = ","), ", scheme=", x$meta$fd_scheme, "\n", sep = "")
  invisible(x)
}
