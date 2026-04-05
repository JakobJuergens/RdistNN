# Lightweight score (kernel) API mirroring the rank/distance style,
# while preserving the S3 class "ustat_score" and your current engine.

# --- tiny helpers -------------------------------------------------------------

# minimal internal helpers (no extra deps)
.is_scalar_int <- function(x) is.numeric(x) && length(x) == 1L && is.finite(x) && as.integer(x) == x
.is_fun <- function(f, n = 0L) is.function(f) && (length(formals(f)) >= n)
.call_with_supported_args <- function(fun, args) {
  fmls <- tryCatch(names(formals(fun)), error = function(e) character())
  if ("..." %in% fmls) {
    return(do.call(fun, args))
  }
  do.call(fun, args[intersect(names(args), fmls)])
}

# derive n from data (matrix/data.frame/vector)
.score_n <- function(data) {
  if (is.matrix(data) || is.data.frame(data)) nrow(data) else length(data)
}

.score_check_idxs <- function(idxs, n, s) {
  .validate_score_idxs(idxs, n, s)
}
.score_check_subs <- function(subs, n, s) {
  .validate_score_subs(subs, n, s)
}

# --- constructors & predicates ------------------------------------------------

#' Score specification (S3) for U-stat kernels
#'
#' A score (kernel) is a symmetric function evaluated on s-sized subsets.
#' Provide either `ptr` (compiled fast path), or an R `eval_batch`, or an R `eval`.
#'
#' If your score depends on a parameter vector (for example Z-estimation), use
#' `eval_param`/`eval_param_batch` and optionally `grad`/`grad_batch` to enable
#' derivative-aware workflows via [score_diff_path()].
#'
#' @param name Character scalar.
#' @param order Integer order s, or function(n, data, ctx) -> s.
#' @param setup Optional function(data, ...) -> ctx object (list or externalptr).
#' @param eval R fallback: function(idxs, data, ctx, rng) -> numeric(1).
#' @param eval_batch R fast path: function(subs, data, ctx, rng) -> numeric(nrow(subs)).
#' @param eval_param Optional parameter-aware scalar evaluator:
#'   function(idxs, data, theta, ctx, rng) -> numeric(1).
#' @param eval_param_batch Optional parameter-aware batch evaluator:
#'   function(subs, data, theta, ctx, rng) -> numeric(nrow(subs)).
#' @param grad Optional analytic scalar gradient:
#'   function(idxs, data, theta, ctx, rng) -> numeric(length(theta)).
#' @param grad_batch Optional analytic batch gradient:
#'   function(subs, data, theta, ctx, rng) -> matrix(nrow(subs), length(theta)).
#' @param fd_step Finite-difference step used by [score_diff_path()] when
#'   `method = "fd"` and no explicit step is supplied.
#' @param fd_scheme Finite-difference scheme used by [score_diff_path()] when
#'   `method = "fd"` and no explicit scheme is supplied. One of
#'   `c("central", "forward")`.
#' @param ptr External pointer to compiled score (see vignette); if present, used by engine.
#' @param randomized Logical: does the kernel draw randomness?
#' @param symmetric Logical: permutation symmetry (engine assumes TRUE).
#' @return An object of class `ustat_score`.
#' @seealso [as_score()], [score_eval()], [score_eval_batch()], [score_diff_path()]
#' @examples
#' sc <- new_score(
#'   name = "sum2",
#'   order = 2L,
#'   eval = function(idxs, data, ctx, rng) sum(data[idxs])
#' )
#' is_score(sc)
#' @export
new_score <- function(
  name,
  order,
  setup = NULL,
  eval = NULL,
  eval_batch = NULL,
  eval_param = NULL,
  eval_param_batch = NULL,
  grad = NULL,
  grad_batch = NULL,
  fd_step = 1e-6,
  fd_scheme = "central",
  ptr = NULL,
  randomized = FALSE,
  symmetric = TRUE
) {
  fd_scheme <- match.arg(fd_scheme, c("central", "forward"))
  sc <- list(
    name = as.character(name),
    order = order,
    setup = setup,
    eval = eval,
    eval_batch = eval_batch,
    eval_param = eval_param,
    eval_param_batch = eval_param_batch,
    grad = grad,
    grad_batch = grad_batch,
    fd_step = as.numeric(fd_step),
    fd_scheme = fd_scheme,
    ptr = ptr,
    randomized = isTRUE(randomized),
    symmetric = isTRUE(symmetric)
  )
  class(sc) <- "ustat_score"
  validate_score(sc)
}

#' Coerce a plain function to a score (simple cases)
#'
#' @param f function(idxs, data, ctx, rng) -> scalar
#' @param s integer order
#' @param name character scalar used as the score name (optional)
#' @return A \code{ustat_score}
#' @seealso [new_score()], [score_eval()]
#' @examples
#' sc <- as_score(function(idxs, data) prod(data[idxs]), s = 2L, name = "prod2")
#' score_eval(sc, c(1L, 3L), c(2, 4, 6), NULL)
#' @export
as_score <- function(f, s, name = deparse(substitute(f))) {
  if (!.is_fun(f, 2L)) {
    .abort_bad_arg("f", "must be a function with at least 2 arguments")
  }
  if (!.is_scalar_int(s)) {
    .abort_bad_arg("s", "must be an integer scalar")
  }
  new_score(
    name = name,
    order = as.integer(s),
    eval = function(idxs, data, ctx, rng) {
      .call_with_supported_args(
        f,
        list(idxs = idxs, data = data, ctx = ctx, rng = rng)
      )
    }
  )
}

#' Is object a \code{ustat_score}?
#' @param x any object
#' @return logical
#' @examples
#' sc <- as_score(function(idxs, data) sum(data[idxs]), s = 2L)
#' is_score(sc)
#' @export
is_score <- function(x) inherits(x, "ustat_score")

#' Does a score have a compiled implementation?
#' @param score a \code{ustat_score}
#' @return logical
#' @examples
#' is_compiled(new_score("r_only", order = 1L, eval = function(idxs, data, ctx, rng) data[idxs]))
#' @export
is_compiled <- function(score) is_score(score) && !is.null(score$ptr)

#' @noRd
validate_score <- function(sc) {
  if (!inherits(sc, "ustat_score")) {
    .abort_bad_arg("sc", "must inherit from class 'ustat_score'")
  }
  if (!nzchar(sc$name)) stop("score name must be non-empty", call. = FALSE)
  if (!(.is_scalar_int(sc$order) || .is_fun(sc$order, 1L))) {
    stop("`order` must be integer or function(n, data, ctx).", call. = FALSE)
  }
  if (!is.null(sc$eval_param) && !is.function(sc$eval_param)) {
    stop("`eval_param` must be NULL or a function.", call. = FALSE)
  }
  if (!is.null(sc$eval_param_batch) && !is.function(sc$eval_param_batch)) {
    stop("`eval_param_batch` must be NULL or a function.", call. = FALSE)
  }
  if (!is.null(sc$grad) && !is.function(sc$grad)) {
    stop("`grad` must be NULL or a function.", call. = FALSE)
  }
  if (!is.null(sc$grad_batch) && !is.function(sc$grad_batch)) {
    stop("`grad_batch` must be NULL or a function.", call. = FALSE)
  }
  if (!is.numeric(sc$fd_step) || length(sc$fd_step) != 1L ||
    !is.finite(sc$fd_step) || sc$fd_step <= 0) {
    stop("`fd_step` must be a positive finite scalar.", call. = FALSE)
  }
  if (!is.character(sc$fd_scheme) || length(sc$fd_scheme) != 1L ||
    !(sc$fd_scheme %in% c("central", "forward"))) {
    stop("`fd_scheme` must be one of: 'central', 'forward'.", call. = FALSE)
  }
  if (is.null(sc$ptr) &&
    is.null(sc$eval_batch) && is.null(sc$eval) &&
    is.null(sc$eval_param_batch) && is.null(sc$eval_param)) {
    stop("Provide at least one of: ptr, eval_batch, eval, eval_param_batch, eval_param", call. = FALSE)
  }
  if (!isTRUE(sc$symmetric)) {
    warning("Score marked non-symmetric; engine assumes symmetry.", call. = FALSE)
  }
  sc
}

# --- pretty printing ----------------------------------------------------------

#' Print a `ustat_score` object
#'
#' @param x A `ustat_score` object.
#' @param ... Not used.
#' @return `x`, invisibly.
#' @examples
#' sc <- as_score(function(idxs, data) sum(data[idxs]), s = 2L, name = "sum2")
#' print(sc)
#' @export
print.ustat_score <- function(x, ...) {
  cat("<ustat_score: ", x$name, ">\n", sep = "")
  cat("  order: ", if (is.function(x$order)) "function(n, data, ctx)" else x$order, "\n", sep = "")
  cat("  impl : ",
    if (!is.null(x$ptr)) {
      "compiled(ptr) "
    } else if (!is.null(x$eval_batch)) {
      "R eval_batch "
    } else if (!is.null(x$eval)) {
      "R eval "
    } else {
      "NONE "
    },
    if (isTRUE(x$randomized)) "(randomized)" else "",
    "\n",
    sep = ""
  )
  has_analytic <- is.function(x$grad) || is.function(x$grad_batch)
  has_param_eval <- is.function(x$eval_param) || is.function(x$eval_param_batch)
  diff_mode <- if (has_analytic) {
    "analytic"
  } else if (has_param_eval) {
    paste0("fd-ready (", x$fd_scheme, ")")
  } else {
    "none"
  }
  cat("  diff : ", diff_mode, "\n", sep = "")
  invisible(x)
}

# --- setup & order ------------------------------------------------------------

#' Set up score context (precomputation)
#' @param score ustat_score
#' @param data data.frame/matrix
#' @param ... passed to score$setup
#' @return context object or NULL
#' @examples
#' sc <- new_score(
#'   name = "ctx_sum",
#'   order = 1L,
#'   setup = function(data, ...) list(offset = 1),
#'   eval = function(idxs, data, ctx, rng) data[idxs] + ctx$offset
#' )
#' ctx <- score_setup(sc, 1:3)
#' score_eval(sc, 2L, 1:3, ctx)
#' @export
score_setup <- function(score, data, ...) UseMethod("score_setup")
#' @export
score_setup.ustat_score <- function(score, data, ...) {
  if (is.null(score$setup)) {
    return(NULL)
  }
  score$setup(data, ...)
}

#' Resolve kernel order \code{s}
#' @param score ustat_score
#' @param n integer sample size
#' @param data data used by the score (for dynamic orders)
#' @param ctx optional context
#' @return integer order \code{s}
#' @examples
#' sc <- new_score(
#'   name = "dyn",
#'   order = function(n, data, ctx) if (n > 5) 3L else 2L,
#'   eval = function(idxs, data, ctx, rng) sum(data[idxs])
#' )
#' score_order(sc, 4L, 1:4, NULL)
#' @export
score_order <- function(score, n, data, ctx = NULL) UseMethod("score_order")
#' @export
score_order.ustat_score <- function(score, n, data, ctx = NULL) {
  if (.is_scalar_int(score$order)) {
    return(as.integer(score$order))
  }
  as.integer(score$order(n, data, ctx))
}

# --- evaluators (robust fallbacks) -------------------------------------------

#' Evaluate kernel on one subset
#'
#' For maximum robustness, this prefers the scalar `eval()` path when available
#' and otherwise falls back to `eval_batch()` for a single-row batch.
#'
#' @param score ustat_score
#' @param idxs integer vector of length \code{s}
#' @param data data.frame/matrix
#' @param ctx optional context
#' @param rng optional RNG key/seed
#' @return numeric(1)
#' @examples
#' sc <- as_score(function(idxs, data) sum(data[idxs]), s = 2L)
#' score_eval(sc, c(1L, 2L), c(3, 4, 5), NULL)
#' @export
score_eval <- function(score, idxs, data, ctx = NULL, rng = NULL) UseMethod("score_eval")
#' @export
score_eval.ustat_score <- function(score, idxs, data, ctx = NULL, rng = NULL) {
  if (is.function(score$eval)) {
    out <- score$eval(idxs, data, ctx, rng)
    if (!is.numeric(out) || length(out) != 1L || !is.finite(out)) {
      stop("score$eval must return a finite scalar numeric.", call. = FALSE)
    }
    return(out)
  }

  if (is.function(score$eval_batch)) {
    subs1 <- if (is.matrix(idxs)) idxs else matrix(as.integer(idxs), nrow = 1L)
    out <- score$eval_batch(subs1, data, ctx, rng)
    if (!is.numeric(out) || length(out) != 1L || !is.finite(out)) {
      stop("score$eval_batch must return a finite scalar for a single subset.", call. = FALSE)
    }
    return(as.numeric(out))
  }

  # Match test: scalar path complains "eval or eval_batch"
  stop("No R `eval` or `eval_batch` provided", call. = FALSE)
}

#' Evaluate kernel on many subsets (R fast-path)
#'
#' @param score ustat_score
#' @param subs integer matrix with nrow=K, ncol=s
#' @param data data.frame/matrix
#' @param ctx optional context
#' @param rng optional RNG key/seed
#' @return numeric vector of length \code{K}
#' @examples
#' sc <- as_score(function(idxs, data) sum(data[idxs]), s = 2L)
#' subs <- rbind(c(1L, 2L), c(2L, 3L))
#' score_eval_batch(sc, subs, c(3, 4, 5), NULL)
#' @export
score_eval_batch <- function(score, subs, data, ctx = NULL, rng = NULL) UseMethod("score_eval_batch")
#' @export
score_eval_batch.ustat_score <- function(score, subs, data, ctx = NULL, rng = NULL) {
  # Fail fast with batch-first wording to match test
  if (!is.function(score$eval_batch) && !is.function(score$eval)) {
    stop("No R `eval_batch` or `eval` provided", call. = FALSE)
  }

  if (is.function(score$eval_batch)) {
    out <- score$eval_batch(subs, data, ctx, rng)
    n_subsets <- if (is.matrix(subs)) nrow(subs) else length(subs)
    if (!is.numeric(out) || length(out) != n_subsets || any(!is.finite(out))) {
      stop("score$eval_batch must return numeric(K) with finite entries.", call. = FALSE)
    }
    return(as.numeric(out))
  }

  # Fallback: validate inputs and loop through scalar eval
  n <- .score_n(data)
  s <- score_order(score, n, data, ctx)
  subsets_matrix <- .score_check_subs(subs, n, s)
  n_subsets <- nrow(subsets_matrix)
  ans <- numeric(n_subsets)
  for (k in seq_len(n_subsets)) {
    ans[k] <- score_eval(score, subsets_matrix[k, ], data, ctx, rng)
  }
  ans
}

# --- registry (dual-form registration, like ranks) ----------------------------

.score_registry <- new.env(parent = emptyenv())

#' Register a score
#'
#' Accepts **either** `(score)` where `score$name` is used as key,
#' or `(name, score)` where `name` overrides/sets `score$name`.
#'
#' @param name Character scalar **or** a score object (legacy single-arg form).
#' @param score Optional `ustat_score` when using the two-arg form.
#' @return Invisibly, the registered name.
#' @examples
#' sc <- as_score(function(idxs, data) sum(data[idxs]), s = 2L, name = "tmp_score_demo")
#' register_score(sc)
#' "tmp_score_demo" %in% list_scores()
#' @export
register_score <- function(name, score = NULL) {
  if (missing(score)) {
    # Legacy: register_score(score)
    score <- name
    if (!is_score(score) || !nzchar(score$name)) {
      .abort_bad_arg("score", "must be a score object with a non-empty `$name`")
    }
    nm <- score$name
  } else {
    # New: register_score(name, score)
    if (!is.character(name) || length(name) != 1L) {
      .abort_bad_arg("name", "must be a character scalar")
    }
    if (!is_score(score)) {
      .abort_bad_arg("score", "must be a `ustat_score` object")
    }
    if (is.null(score$name) || !nzchar(score$name)) score$name <- name
    nm <- name
  }
  assign(nm, score, envir = .score_registry)
  invisible(nm)
}

#' Retrieve a registered score
#' @param name character scalar
#' @return a \code{ustat_score}
#' @examples
#' sc0 <- as_score(function(idxs, data) sum(data[idxs]), s = 2L, name = "tmp_score_get")
#' register_score(sc0)
#' sc <- get_score("tmp_score_get")
#' is_score(sc)
#' @export
get_score <- function(name) {
  if (!exists(name, envir = .score_registry, inherits = FALSE)) {
    stop("No score registered under '", name, "'.", call. = FALSE)
  }
  get(name, envir = .score_registry, inherits = FALSE)
}

#' List registered scores
#' @return character vector of names
#' @examples
#' head(list_scores(), 3)
#' @export
list_scores <- function() sort(ls(envir = .score_registry, all.names = FALSE))

#' Back-compat alias for [list_scores()]
#' @return character vector of names
#' @examples
#' identical(score_names(), list_scores())
#' @export
score_names <- function() list_scores()
