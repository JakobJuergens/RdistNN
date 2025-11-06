# R/rank-api.R
# -----------------------------------------------------------------------------
# Lightweight rank-weight API + tiny registry (non-S3), consistent with distance API

.rank_registry <- new.env(parent = emptyenv())

# ---- predicates & adapters ---------------------------------------------------

#' Is this a rank_weight spec?
#'
#' A rank-weight spec is a list with `eval(ranks, ctx)` and/or
#' `eval_batch(R, ctx)`, where `R` is an n x m rank matrix (each column a
#' permutation of 1..n). Optionally it may define `setup(data)`.
#'
#' @param x Any object.
#' @return Logical.
#' @export
is_rank_weight <- function(x) {
  inherits(x, "rank_weight") ||
    (is.list(x) && identical(x$.__kind, "rank_weight"))
}

# -----------------------------------------------------------------------------
.as_fun_accepts_ctx <- function(fun) {
  fmls <- formals(fun)
  if (is.null(fmls)) {
    return(FALSE)
  }
  # Accept either explicitly named 'ctx' or any >= 2-arg signature
  "ctx" %in% names(fmls) || length(fmls) >= 2L
}

#' Wrap a function as a rank_weight spec
#'
#' @param fun Function `(ranks, ctx) -> numeric(n)` OR just `(ranks) -> numeric(n)`.
#' @param name Character name.
#' @return A rank-weight spec (list) with fields `name`, `eval`.
#' @export
as_rank_weight <- function(fun, name = "user_rank_weight") {
  stopifnot(is.function(fun))
  takes_ctx <- .as_fun_accepts_ctx(fun)

  spec <- list(
    name = name,
    eval = function(ranks, ctx = NULL) {
      r <- .rank_check_vec(ranks)
      w <- if (takes_ctx) fun(r, ctx) else fun(r)
      if (!is.numeric(w) || length(w) != length(r) || any(!is.finite(w))) {
        stop("Rank-weight function must return finite numeric length(ranks).", call. = FALSE)
      }
      as.numeric(w)
    }
  )
  spec$.__kind <- "rank_weight"
  class(spec) <- c("rank_weight", class(spec))
  spec
}

# ---- registry ---------------------------------------------------------------

#' Register a rank_weight under a name
#'
#' Accepts **either** `(name, spec)` or legacy single-arg form `(spec)` where
#' `spec$name` is used as the key. If `spec` is a function, it is wrapped with
#' [as_rank_weight()].
#'
#' @param name Character scalar or a spec (legacy).
#' @param spec A rank_weight spec or function `(ranks, ctx) -> numeric(n)`.
#' @return Invisibly, the registered spec.
#' @export
register_rank <- function(name, spec = NULL) {
  if (is.null(spec)) {
    # legacy form: register_rank(spec)
    spec <- name
    if (is.function(spec)) {
      nm <- "user_rank_weight"
      spec <- as_rank_weight(spec, name = nm)
      name <- nm
    } else if (is.list(spec) && is.character(spec$name) && length(spec$name) == 1L) {
      name <- spec$name
    } else {
      stop("Legacy register_rank(spec): spec must have a non-empty $name.", call. = FALSE)
    }
  } else {
    # new form: register_rank(name, spec)
    stopifnot(is.character(name), length(name) == 1L)
    if (is.function(spec)) spec <- as_rank_weight(spec, name = name)
    if (is.null(spec$name)) spec$name <- name
  }
  if (!is_rank_weight(spec)) stop("`spec` is not a rank_weight.", call. = FALSE)
  assign(name, spec, envir = .rank_registry)
  invisible(spec)
}

#' Retrieve a registered rank_weight
#' @param name Character scalar.
#' @return A rank_weight spec.
#' @export
get_rank <- function(name) {
  stopifnot(is.character(name), length(name) == 1L)
  if (!exists(name, envir = .rank_registry, inherits = FALSE)) {
    known <- ls(.rank_registry, all.names = FALSE)
    stop(sprintf("Unknown rank_weight '%s'. Known: %s", name, paste(known, collapse = ", ")), call. = FALSE)
  }
  get(name, envir = .rank_registry, inherits = FALSE)
}

#' List registered rank_weight names
#' @return Character vector.
#' @export
rank_names <- function() {
  ls(.rank_registry, all.names = FALSE)
}

#' Back-compat alias for [rank_names()]
#' @export
list_ranks <- function() rank_names()

# ---- evaluators -------------------------------------------------------------

# internal validators used by the evaluators
.rank_check_vec <- function(ranks) {
  r <- as.integer(ranks)
  n <- length(r)
  if (n == 0L) stop("`ranks` must have positive length.", call. = FALSE)
  if (anyNA(r)) stop("`ranks` cannot contain NA.", call. = FALSE)
  if (!setequal(r, seq_len(n))) stop("`ranks` must be a permutation of 1..n.", call. = FALSE)
  r
}
.rank_check_mat <- function(R) {
  R <- as.matrix(R)
  if (!is.integer(R)) storage.mode(R) <- "integer"
  if (nrow(R) == 0L || ncol(R) == 0L) stop("`R` must have positive nrow and ncol.", call. = FALSE)
  for (j in seq_len(ncol(R))) .rank_check_vec(R[, j])
  R
}

#' Evaluate rank-weight for a single target
#'
#' @param spec Rank-weight spec (see [is_rank_weight()]).
#' @param ranks Integer permutation of `1..n` (length `n`).
#' @param ctx Optional context (if the spec uses it).
#' @return Numeric vector of length `n`.
#' @export
rank_weight_eval <- function(spec, ranks, ctx = NULL) {
  if (!is_rank_weight(spec)) stop("`spec` must be a rank_weight.", call. = FALSE)
  r <- .rank_check_vec(ranks)

  if (is.function(spec$eval_batch)) {
    # Prefer batch path whenever available (robust to accidental/legacy evals)
    R1 <- matrix(r, nrow = length(r), ncol = 1L)
    W1 <- rank_weight_eval_batch(spec, R1, ctx)
    w <- W1[, 1L]
  } else if (is.function(spec$eval)) {
    w <- spec$eval(r, ctx)
  } else {
    stop("rank_weight spec has neither eval nor eval_batch.", call. = FALSE)
  }

  if (!is.numeric(w) || length(w) != length(r) || any(!is.finite(w))) {
    stop("rank weight must return finite numeric of length(ranks).", call. = FALSE)
  }
  as.numeric(w)
}

#' Evaluate rank-weight for many targets (batch)
#'
#' @param spec Rank-weight spec (see [is_rank_weight()]).
#' @param R Integer matrix (n x m), each column a permutation of `1..n`.
#' @param ctx Optional context (if the spec uses it).
#' @return Numeric matrix (n x m).
#' @export
rank_weight_eval_batch <- function(spec, R, ctx = NULL) {
  if (!is_rank_weight(spec)) stop("`spec` must be a rank_weight.", call. = FALSE)
  RR <- .rank_check_mat(R)

  if (is.function(spec$eval_batch)) {
    W <- spec$eval_batch(RR, ctx)
  } else if (is.function(spec$eval)) {
    # Fallback: apply eval to each column; coerce to matrix if m == 1
    W <- apply(RR, 2L, function(r) spec$eval(r, ctx))
    if (!is.matrix(W)) W <- matrix(W, nrow = nrow(RR), ncol = 1L)
  } else {
    stop("rank_weight spec has neither eval nor eval_batch.", call. = FALSE)
  }

  if (!is.matrix(W) || any(!is.finite(W)) ||
    nrow(W) != nrow(RR) || ncol(W) != ncol(RR)) {
    stop("rank_weight$eval_batch must return an nrow(R) x ncol(R) finite matrix.", call. = FALSE)
  }
  W
}
