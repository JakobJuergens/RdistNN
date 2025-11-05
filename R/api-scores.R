# R/score-api.R
# -----------------------------------------------------------------------------
# Lightweight score (kernel) API mirroring the rank/distance style,
# while preserving the S3 class "ustat_score" and your current engine.

# --- tiny helpers -------------------------------------------------------------

# minimal internal helpers (no extra deps)
.is_scalar_int <- function(x) is.numeric(x) && length(x) == 1L && is.finite(x) && as.integer(x) == x
.is_fun <- function(f, n = 0L) is.function(f) && (length(formals(f)) >= n)

# derive n from data (matrix/data.frame/vector)
.score_n <- function(data) {
  if (is.matrix(data) || is.data.frame(data)) nrow(data) else length(data)
}

.score_check_idxs <- function(idxs, n, s) {
  idxs <- as.integer(idxs)
  if (length(idxs) != s) stop(sprintf("`idxs` must have %d elements (order s).", s), call. = FALSE)
  if (anyNA(idxs))       stop("`idxs` cannot contain NA.", call. = FALSE)
  if (any(idxs < 1L | idxs > n)) stop("`idxs` must lie in 1..n.", call. = FALSE)  # <- normalized
  idxs
}
.score_check_subs <- function(subs, n, s) {
  if (!is.matrix(subs)) subs <- as.matrix(subs)
  if (ncol(subs) != s)  stop(sprintf("`subs` must have %d columns (order s).", s), call. = FALSE)
  if (nrow(subs) == 0L) stop("`subs` must have at least one row.", call. = FALSE)
  if (!is.integer(subs)) storage.mode(subs) <- "integer"
  if (any(!is.finite(subs))) stop("`subs` contains non-finite entries.", call. = FALSE)
  if (any(subs < 1L | subs > n)) stop("`subs` indices must lie in 1..n.", call. = FALSE)  # <- normalized
  subs
}

# --- constructors & predicates ------------------------------------------------

#' Score specification (S3) for U-stat kernels
#'
#' A score (kernel) is a symmetric function evaluated on s-sized subsets.
#' Provide either `ptr` (compiled fast path), or an R `eval_batch`, or an R `eval`.
#'
#' @param name Character scalar.
#' @param order Integer order s, or function(n, data, ctx) -> s.
#' @param setup Optional function(data, ...) -> ctx object (list or externalptr).
#' @param eval R fallback: function(idxs, data, ctx, rng) -> numeric(1).
#' @param eval_batch R fast path: function(subs, data, ctx, rng) -> numeric(nrow(subs)).
#' @param ptr External pointer to compiled score (see vignette); if present, used by engine.
#' @param randomized Logical: does the kernel draw randomness?
#' @param symmetric Logical: permutation symmetry (engine assumes TRUE).
#' @return An object of class `ustat_score`.
#' @export
new_score <- function(
    name,
    order,
    setup = NULL,
    eval = NULL,
    eval_batch = NULL,
    ptr = NULL,
    randomized = FALSE,
    symmetric = TRUE
) {
  sc <- list(
    name = as.character(name),
    order = order,
    setup = setup,
    eval = eval,
    eval_batch = eval_batch,
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
#' @export
as_score <- function(f, s, name = deparse(substitute(f))) {
  stopifnot(.is_fun(f, 2L), .is_scalar_int(s))
  new_score(name = name, order = as.integer(s), eval = f)
}

#' Is object a \code{ustat_score}?
#' @param x any object
#' @return logical
#' @export
is_score <- function(x) inherits(x, "ustat_score")

#' Does a score have a compiled implementation?
#' @param score a \code{ustat_score}
#' @return logical
#' @export
is_compiled <- function(score) is_score(score) && !is.null(score$ptr)

#' @noRd
validate_score <- function(sc) {
  stopifnot(inherits(sc, "ustat_score"))
  if (!nzchar(sc$name)) stop("score name must be non-empty", call. = FALSE)
  if (!(.is_scalar_int(sc$order) || .is_fun(sc$order, 1L))) {
    stop("`order` must be integer or function(n, data, ctx).", call. = FALSE)
  }
  if (is.null(sc$ptr) && is.null(sc$eval_batch) && is.null(sc$eval)) {
    stop("Provide at least one of: ptr, eval_batch, eval", call. = FALSE)
  }
  if (!isTRUE(sc$symmetric)) {
    warning("Score marked non-symmetric; engine assumes symmetry.", call. = FALSE)
  }
  sc
}

# Optional parity helper (export if you like)
score_names <- function() list_scores()

# --- pretty printing ----------------------------------------------------------

#' @export
print.ustat_score <- function(x, ...) {
  cat("<ustat_score: ", x$name, ">\n", sep = "")
  cat("  order: ", if (is.function(x$order)) "function(n, data, ctx)" else x$order, "\n", sep = "")
  cat("  impl : ",
      if (!is.null(x$ptr))       "compiled(ptr) "
      else if (!is.null(x$eval_batch)) "R eval_batch "
      else if (!is.null(x$eval)) "R eval "
      else                       "NONE ",
      if (isTRUE(x$randomized)) "(randomized)" else "",
      "\n", sep = "")
  invisible(x)
}

# --- setup & order ------------------------------------------------------------

#' Set up score context (precomputation)
#' @param score ustat_score
#' @param data data.frame/matrix
#' @param ... passed to score$setup
#' @return context object or NULL
#' @export
score_setup <- function(score, data, ...) UseMethod("score_setup")
#' @export
score_setup.ustat_score <- function(score, data, ...) {
  if (is.null(score$setup)) return(NULL)
  score$setup(data, ...)
}

#' Resolve kernel order \code{s}
#' @param score ustat_score
#' @param n integer sample size
#' @param data data used by the score (for dynamic orders)
#' @param ctx optional context
#' @return integer order \code{s}
#' @export
score_order <- function(score, n, data, ctx = NULL) UseMethod("score_order")
#' @export
score_order.ustat_score <- function(score, n, data, ctx = NULL) {
  if (.is_scalar_int(score$order)) return(as.integer(score$order))
  as.integer(score$order(n, data, ctx))
}

# --- evaluators (robust fallbacks) -------------------------------------------

#' Evaluate kernel on one subset
#'
#' For maximum robustness, this prefers the batch path when available and
#' otherwise falls back to the scalar path.
#'
#' @param score ustat_score
#' @param idxs integer vector of length \code{s}
#' @param data data.frame/matrix
#' @param ctx optional context
#' @param rng optional RNG key/seed
#' @return numeric(1)
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
    K <- if (is.matrix(subs)) nrow(subs) else length(subs)
    if (!is.numeric(out) || length(out) != K || any(!is.finite(out))) {
      stop("score$eval_batch must return numeric(K) with finite entries.", call. = FALSE)
    }
    return(as.numeric(out))
  }

  # Fallback: validate inputs and loop through scalar eval
  n <- .score_n(data)
  s <- score_order(score, n, data, ctx)
  SS <- .score_check_subs(subs, n, s)
  K  <- nrow(SS)
  ans <- numeric(K)
  for (k in seq_len(K)) {
    ans[k] <- score_eval(score, SS[k, ], data, ctx, rng)
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
#' @export
register_score <- function(name, score = NULL) {
  if (is.null(score)) {
    # Legacy: register_score(score)
    score <- name
    stopifnot(is_score(score), nzchar(score$name))
    nm <- score$name
  } else {
    # New: register_score(name, score)
    stopifnot(is.character(name), length(name) == 1L, is_score(score))
    if (is.null(score$name) || !nzchar(score$name)) score$name <- name
    nm <- name
  }
  assign(nm, score, envir = .score_registry)
  invisible(nm)
}

#' Retrieve a registered score
#' @param name character scalar
#' @return a \code{ustat_score}
#' @export
get_score <- function(name) {
  if (!exists(name, envir = .score_registry, inherits = FALSE)) {
    stop("No score registered under '", name, "'.", call. = FALSE)
  }
  get(name, envir = .score_registry, inherits = FALSE)
}

#' List registered scores
#' @return character vector of names
#' @export
list_scores <- function() sort(ls(envir = .score_registry, all.names = FALSE))

#' Back-compat alias for [list_scores()]
#' @return character vector of names
#' @export
score_names <- function() list_scores()
