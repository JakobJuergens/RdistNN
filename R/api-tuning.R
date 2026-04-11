# Unified API + registry for kernel-order tuning specifications

# lightweight registry (name -> spec)
.tuning_registry <- new.env(parent = emptyenv())

# --- type + adapters ----------------------------------------------------------
.tuning_has_eval <- function(spec) !is.null(spec$eval) && is.function(spec$eval)
.tuning_has_eval_grid <- function(spec) !is.null(spec$eval_grid) && is.function(spec$eval_grid)

.tuning_call_user <- function(fun, args) {
  fmls <- tryCatch(names(formals(fun)), error = function(e) character())
  if ("..." %in% fmls) {
    return(do.call(fun, args))
  }
  do.call(fun, args[intersect(names(args), fmls)])
}

.tuning_validate_spec <- function(spec) {
  if (!is.list(spec)) .abort_bad_arg("spec", "must be a list")
  if (!(.tuning_has_eval(spec) || .tuning_has_eval_grid(spec))) {
    .abort_bad_arg("spec", "must provide at least one of `eval` or `eval_grid`")
  }
  if (is.null(spec$name) || !nzchar(as.character(spec$name))) {
    spec$name <- "unnamed_tuning"
  }
  class(spec) <- union("dnn_tuning", class(spec))
  spec
}

.tuning_needs_jackknife <- function(se, jackknife) {
  isTRUE(se) && !identical(jackknife, "none")
}

.tuning_round_int <- function(x, rounding = c("floor", "ceiling", "round")) {
  rounding <- match.arg(rounding)
  if (identical(rounding, "floor")) return(as.integer(floor(x)))
  if (identical(rounding, "ceiling")) return(as.integer(ceiling(x)))
  as.integer(round(x))
}

.tuning_validate_s_scalar_output <- function(x, arg) {
  if (
    !is.numeric(x) ||
      length(x) != 1L ||
      !is.finite(x) ||
      x != floor(x) ||
      x < 1L ||
      x > .Machine$integer.max
  ) {
    .abort_bad_arg(
      arg,
      "must return a positive finite integer-like numeric scalar",
      quoted = FALSE
    )
  }
  as.integer(x)
}

.tuning_validate_s_grid_output <- function(x, arg) {
  if (
    !is.numeric(x) ||
      length(x) < 1L ||
      any(!is.finite(x)) ||
      any(x != floor(x)) ||
      any(x < 1L) ||
      any(x > .Machine$integer.max)
  ) {
    .abort_bad_arg(
      arg,
      "column `s` must be positive finite integer-like numeric values",
      quoted = FALSE
    )
  }
  as.integer(x)
}

.tuning_validate_default_output <- function(
  out, n, kx, k, se, jackknife, name
) {
  if (is.numeric(out) && length(out) == 1L) {
    out <- list(s = .tuning_validate_s_scalar_output(out, arg = "tuning$eval()"))
  }
  if (!is.list(out) || is.null(out$s)) {
    .abort_bad_arg(
      "tuning$eval()",
      "must return a list with field `s` or a numeric scalar",
      quoted = FALSE
    )
  }
  s <- .validate_scalar_int(out$s, arg = "s")
  needs_jk <- .tuning_needs_jackknife(se = se, jackknife = jackknife)
  feasible <- list(
    k_le_s = as.integer(k) <= as.integer(s),
    s_lt_n = if (needs_jk) as.integer(s) < as.integer(n) else TRUE
  )
  meta <- out$meta
  if (is.null(meta) || !is.list(meta)) meta <- list()
  meta_name <- meta[["name", exact = TRUE]]
  meta_n <- meta[["n", exact = TRUE]]
  meta_kx <- meta[["kx", exact = TRUE]]
  meta$name <- as.character(if (is.null(meta_name)) name else meta_name)
  meta_n <- if (is.null(meta_n)) n else meta_n
  meta_kx <- if (is.null(meta_kx)) kx else meta_kx
  meta$n <- as.integer(meta_n)
  meta$kx <- as.integer(meta_kx)
  list(
    s = as.integer(s),
    meta = meta,
    feasibility = feasible
  )
}

.tuning_validate_grid_output <- function(
  out, n, kx, k, se, jackknife, name
) {
  if (is.numeric(out) && length(out) == 0L) {
    .abort_bad_arg("tuning$eval_grid()", "must return at least one candidate row", quoted = FALSE)
  }
  if (is.numeric(out) && length(out) > 0L) {
    out <- data.frame(
      kappa = NA_real_,
      s = .tuning_validate_s_grid_output(out, arg = "tuning$eval_grid()"),
      stringsAsFactors = FALSE
    )
  }
  if (!is.data.frame(out)) {
    .abort_bad_arg(
      "tuning$eval_grid()",
      "must return a data.frame with at least column `s`",
      quoted = FALSE
    )
  }
  if (!("s" %in% names(out))) {
    .abort_bad_arg(
      "tuning$eval_grid()",
      "must return a data.frame with column `s`",
      quoted = FALSE
    )
  }
  if (nrow(out) < 1L) {
    .abort_bad_arg("tuning$eval_grid()", "must return at least one candidate row", quoted = FALSE)
  }
  out$s <- .tuning_validate_s_grid_output(out$s, arg = "tuning$eval_grid()")
  needs_jk <- .tuning_needs_jackknife(se = se, jackknife = jackknife)
  out$feasible <- (as.integer(k) <= out$s) & (!needs_jk | (out$s < as.integer(n)))
  out$reason <- ifelse(
    out$feasible,
    "ok",
    ifelse(as.integer(k) > out$s, "k > s", "jackknife requires s < n")
  )
  list(
    candidates = out,
    meta = list(
      name = as.character(name),
      n = as.integer(n),
      kx = as.integer(kx)
    )
  )
}

.tuning_eval_default <- function(tuning, n, kx, k, se, jackknife, ctx) {
  args <- list(
    n = as.integer(n),
    kx = as.integer(kx),
    k = as.integer(k),
    se = isTRUE(se),
    jackknife = jackknife,
    ctx = ctx
  )

  if (!.tuning_has_eval(tuning)) {
    if (!.tuning_has_eval_grid(tuning)) {
      .abort_bad_arg("tuning", "has no `eval` or `eval_grid` method")
    }
    grid <- .tuning_call_user(tuning$eval_grid, args)
    grid_ok <- .tuning_validate_grid_output(
      grid,
      n = n,
      kx = kx,
      k = k,
      se = se,
      jackknife = jackknife,
      name = tuning$name
    )
    if (nrow(grid_ok$candidates) < 1L) {
      .abort_bad_arg("tuning$eval_grid()", "must return at least one candidate row", quoted = FALSE)
    }
    return(.tuning_validate_default_output(
      list(
        s = grid_ok$candidates$s[1L],
        meta = grid_ok$meta
      ),
      n = n,
      kx = kx,
      k = k,
      se = se,
      jackknife = jackknife,
      name = tuning$name
    ))
  }

  out <- .tuning_call_user(tuning$eval, args)
  .tuning_validate_default_output(
    out,
    n = n,
    kx = kx,
    k = k,
    se = se,
    jackknife = jackknife,
    name = tuning$name
  )
}

.tuning_eval_grid <- function(tuning, n, kx, k, se, jackknife, ctx) {
  args <- list(
    n = as.integer(n),
    kx = as.integer(kx),
    k = as.integer(k),
    se = isTRUE(se),
    jackknife = jackknife,
    ctx = ctx
  )
  if (.tuning_has_eval_grid(tuning)) {
    out <- .tuning_call_user(tuning$eval_grid, args)
    return(.tuning_validate_grid_output(
      out,
      n = n,
      kx = kx,
      k = k,
      se = se,
      jackknife = jackknife,
      name = tuning$name
    ))
  }

  # Fall back to one-row grid from eval()
  def <- .tuning_eval_default(
    tuning = tuning,
    n = n,
    kx = kx,
    k = k,
    se = se,
    jackknife = jackknife,
    ctx = ctx
  )
  .tuning_validate_grid_output(
    data.frame(kappa = NA_real_, s = def$s, stringsAsFactors = FALSE),
    n = n,
    kx = kx,
    k = k,
    se = se,
    jackknife = jackknife,
    name = tuning$name
  )
}

#' Is this a tuning spec?
#'
#' A tuning spec is a list with an `eval(n, kx, k, se, jackknife, ctx)` function
#' and optional `eval_grid(...)` for candidate generation.
#'
#' @param x Any R object.
#' @return Logical.
#' @examples
#' is_tuning(tuning_dnn_paper())
#' @export
is_tuning <- function(x) {
  inherits(x, "dnn_tuning") ||
    (is.list(x) && (is.function(x$eval) || is.function(x$eval_grid)))
}

#' Construct a tuning spec
#'
#' @param name Character scalar.
#' @param eval Optional function `(n, kx, k, se, jackknife, ctx) -> list|numeric`.
#' @param eval_grid Optional function `(n, kx, k, se, jackknife, ctx) -> data.frame`.
#' @param setup Optional setup function for future use.
#' @return A validated `dnn_tuning` spec.
#' @examples
#' tun <- new_tuning(
#'   name = "fixed_s",
#'   eval = function(n, kx, k, se, jackknife, ctx = NULL) list(s = 5L)
#' )
#' is_tuning(tun)
#' @export
new_tuning <- function(name, eval = NULL, eval_grid = NULL, setup = NULL) {
  .tuning_validate_spec(list(
    name = as.character(name),
    eval = eval,
    eval_grid = eval_grid,
    setup = setup
  ))
}

#' Wrap a function as a tuning spec
#'
#' @param fun Function that returns either a scalar `s` or a list with field `s`.
#' @param name Character name for the spec.
#' @return A `dnn_tuning` spec.
#' @examples
#' tun <- as_tuning(function(n, kx, ...) floor(n^((kx + 1)/(kx + 2))))
#' tuning_eval(tun, n = 200, kx = 2, k = 5)
#' @export
as_tuning <- function(fun, name = "user_tuning") {
  if (!is.function(fun)) {
    .abort_bad_arg("fun", "must be a function")
  }
  new_tuning(
    name = name,
    eval = function(n, kx, k, se, jackknife, ctx = NULL) {
      .tuning_call_user(
        fun,
        list(n = n, kx = kx, k = k, se = se, jackknife = jackknife, ctx = ctx)
      )
    }
  )
}

#' Register a tuning spec under a name
#'
#' @param name Character scalar registry key.
#' @param tuning A tuning spec or function coercible via [as_tuning()].
#' @return (Invisibly) the registered spec.
#' @examples
#' register_tuning("tmp_fixed_4", function(...) 4L)
#' "tmp_fixed_4" %in% tuning_names()
#' @export
register_tuning <- function(name, tuning) {
  if (!is.character(name) || length(name) != 1L || !nzchar(name)) {
    .abort_bad_arg("name", "must be a non-empty character scalar")
  }
  if (is.function(tuning)) tuning <- as_tuning(tuning, name = name)
  if (is.list(tuning) && (is.null(tuning$name) || !nzchar(as.character(tuning$name)))) {
    tuning$name <- name
  }
  tuning <- .tuning_validate_spec(tuning)
  assign(name, tuning, envir = .tuning_registry)
  invisible(tuning)
}

#' Retrieve a registered tuning spec
#'
#' @param name Character scalar.
#' @return A tuning spec.
#' @examples
#' get_tuning("paper_undersmooth")
#' @export
get_tuning <- function(name) {
  if (!is.character(name) || length(name) != 1L) {
    .abort_bad_arg("name", "must be a character scalar")
  }
  if (!exists(name, envir = .tuning_registry, inherits = FALSE)) {
    known <- ls(.tuning_registry, all.names = FALSE)
    .abort_bad_arg(
      "name",
      "Unknown tuning",
      expected = if (length(known) > 0L) paste(known, collapse = ", ") else "<none registered>",
      actual = as.character(name)
    )
  }
  get(name, envir = .tuning_registry, inherits = FALSE)
}

#' List registered tuning names
#' @return Character vector of names.
#' @examples
#' tuning_names()
#' @export
tuning_names <- function() {
  sort(ls(.tuning_registry, all.names = FALSE))
}

#' Evaluate a tuning spec
#'
#' @param tuning Tuning spec or registered name.
#' @param n Integer sample size.
#' @param kx Integer localization covariate dimension.
#' @param k Integer DNN neighbor count.
#' @param se Logical; whether jackknife SEs are requested.
#' @param jackknife One of `"fast"`, `"loo"`, or `"none"`.
#' @param ctx Optional context passed to tuning evaluators.
#' @param mode One of `"default"` or `"grid"`.
#'
#' @return
#' - `mode = "default"`: list with `s`, `meta`, and `feasibility`.
#' - `mode = "grid"`: list with `candidates` data.frame and `meta`.
#'
#' @examples
#' tuning_eval("paper_undersmooth", n = 500, kx = 2, k = 5)
#' tuning_eval("paper_undersmooth", n = 500, kx = 2, k = 5, mode = "grid")
#' @export
tuning_eval <- function(
  tuning,
  n,
  kx,
  k,
  se = FALSE,
  jackknife = "none",
  ctx = NULL,
  mode = c("default", "grid")
) {
  mode <- match.arg(mode)
  n <- .validate_scalar_int(n, arg = "n")
  kx <- .validate_scalar_int(kx, arg = "kx")
  k <- .validate_scalar_int(k, arg = "k")
  if (!is.logical(se) || length(se) != 1L || is.na(se)) {
    .abort_bad_arg("se", "must be TRUE or FALSE")
  }
  if (!is.character(jackknife) || length(jackknife) != 1L) {
    .abort_bad_arg("jackknife", "must be a character scalar")
  }
  jackknife <- match.arg(jackknife, choices = c("fast", "loo", "none"))

  tun <- if (is.character(tuning)) {
    get_tuning(tuning)
  } else if (is.function(tuning)) {
    as_tuning(tuning)
  } else {
    tuning
  }
  if (!is_tuning(tun)) {
    .abort_bad_arg("tuning", "must be a tuning spec, registered name, or function")
  }

  if (identical(mode, "default")) {
    return(.tuning_eval_default(
      tuning = tun,
      n = n,
      kx = kx,
      k = k,
      se = se,
      jackknife = jackknife,
      ctx = ctx
    ))
  }
  .tuning_eval_grid(
    tuning = tun,
    n = n,
    kx = kx,
    k = k,
    se = se,
    jackknife = jackknife,
    ctx = ctx
  )
}

#' Paper-default undersmoothing rule for kernel order
#'
#' Builds a tuning spec implementing
#' \deqn{s(\kappa) = \lfloor \kappa n^{(k_x+1)/(k_x+2)} \rfloor}
#' with default \eqn{\kappa = 1} and optional candidate grid.
#'
#' @param kappa Numeric scalar tuning constant.
#' @param kappa_grid Numeric vector of candidate constants for grid mode.
#' @param rounding One of `"floor"`, `"ceiling"`, or `"round"`.
#' @param kx_override Optional positive integer to override runtime `kx`.
#' @return A `dnn_tuning` spec.
#' @examples
#' tun <- tuning_dnn_paper()
#' tuning_eval(tun, n = 1000, kx = 2, k = 5)
#' @export
tuning_dnn_paper <- function(
  kappa = 1,
  kappa_grid = c(0.75, 1, 1.25, 1.5),
  rounding = "floor",
  kx_override = NULL
) {
  if (!is.numeric(kappa) || length(kappa) != 1L || !is.finite(kappa) || kappa <= 0) {
    .abort_bad_arg("kappa", "must be a positive finite numeric scalar")
  }
  if (!is.numeric(kappa_grid) || length(kappa_grid) < 1L || any(!is.finite(kappa_grid)) || any(kappa_grid <= 0)) {
    .abort_bad_arg("kappa_grid", "must be a positive finite numeric vector")
  }
  rounding <- match.arg(rounding, choices = c("floor", "ceiling", "round"))
  if (!is.null(kx_override)) {
    kx_override <- .validate_scalar_int(kx_override, arg = "kx_override")
  }

  .compute_s <- function(n, kx, kap) {
    kx_eff <- if (is.null(kx_override)) as.integer(kx) else as.integer(kx_override)
    exponent <- (kx_eff + 1) / (kx_eff + 2)
    raw <- as.numeric(kap) * as.numeric(n)^exponent
    s <- .tuning_round_int(raw, rounding = rounding)
    # no silent clipping; retain explicit feasibility checks downstream
    s
  }

  new_tuning(
    name = "paper_undersmooth",
    eval = function(n, kx, k, se, jackknife, ctx = NULL) {
      s <- .compute_s(n = n, kx = kx, kap = kappa)
      list(
        s = s,
        meta = list(
          kappa = as.numeric(kappa),
          formula = "kappa * n^((kx + 1)/(kx + 2))",
          rounding = rounding
        )
      )
    },
    eval_grid = function(n, kx, k, se, jackknife, ctx = NULL) {
      s_vec <- vapply(kappa_grid, function(kap) .compute_s(n = n, kx = kx, kap = kap), integer(1L))
      data.frame(
        kappa = as.numeric(kappa_grid),
        s = as.integer(s_vec),
        stringsAsFactors = FALSE
      )
    }
  )
}

#' @noRd
.register_builtin_tuning <- function() {
  register_tuning("paper_undersmooth", tuning_dnn_paper())
}
