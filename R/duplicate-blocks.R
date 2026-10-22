# Exact raw-covariate groups. Sorting compares represented numeric values;
# string keys, rounding, distances and fitted signals never define membership.
.dnn_tie_policy <- function() "distance_then_covariates_with_duplicate_block_averaging"

.dnn_duplicate_groups <- function(X) {
  n <- nrow(X)
  if (n < 2L) return(seq_len(n))
  ord <- order(.covariate_tie_key(X))
  different <- rowSums(X[ord[-1L], , drop = FALSE] !=
                        X[ord[-n], , drop = FALSE]) != 0L
  group <- integer(n)
  group[ord] <- cumsum(c(TRUE, different))
  group
}

.dnn_average_duplicate_weights <- function(W, groups) {
  if (is.null(groups) || !anyDuplicated(groups)) return(W)
  for (ids in split(seq_along(groups), groups)) {
    if (length(ids) < 2L) next
    for (j in seq_len(ncol(W))) {
      W[ids, j] <- .dnn_duplicate_block_mean(W[ids, j])
    }
  }
  W
}

# Sum first to retain subnormal contributions; divide first only when needed
# to keep an otherwise finite block mean inside the floating-point range.
.dnn_duplicate_block_mean <- function(values) {
  total <- .compensated_sum(values)
  if (is.finite(total)) return(total / length(values))
  .compensated_sum(values / length(values))
}

# Sum the deleted table over this group's remaining ranks. The usual geometry
# gives a contiguous block and one common mass; supplied ranks may interleave
# identical raw rows, so keep a linear compensated prefix/suffix reference.
.dnn_duplicate_deleted_mass <- function(rr, table) {
  h <- length(rr)
  if (all(diff(rr) == 1L)) {
    return(rep(.compensated_sum(table[seq.int(rr[1L], rr[h] - 1L)]), h))
  }
  before <- after <- before_low <- after_low <- numeric(h)
  high <- low <- 0
  for (i in seq_len(h - 1L)) {
    state <- .compensated_add_R(high, low, table[rr[i]])
    high <- state$high
    low <- state$low
    before[i + 1L] <- high
    before_low[i + 1L] <- low
  }
  high <- low <- 0
  for (i in seq.int(h, 2L)) {
    state <- .compensated_add_R(high, low, table[rr[i] - 1L])
    high <- state$high
    low <- state$low
    after[i - 1L] <- high
    after_low[i - 1L] <- low
  }
  .compensated_merge_R(before, after, before_low, after_low, 0)
}

# State immediately before each represented product; optional block resets
# also supply within-group survivor states. Reverse the inputs for suffixes.
.dnn_block_states_R <- function(values, reset = NULL) {
  high <- low <- numeric(length(values))
  leading <- correction <- 0
  for (i in seq_along(values)) {
    if (!is.null(reset) && reset[i]) leading <- correction <- 0
    high[i] <- leading
    low[i] <- correction
    value <- values[i]
    next_high <- leading + value
    if (is.finite(next_high)) {
      correction <- correction + if (abs(leading) >= abs(value)) {
        (leading - next_high) + value
      } else {
        (value - next_high) + leading
      }
    } else {
      correction <- 0
    }
    leading <- next_high
  }
  list(high = high, low = low)
}

# Explicit direct regrouping is the independent algorithm for a requested
# reference path, and the O(n^2*q) boundary for noncanonical supplied ranks.
.dnn_block_loo_direct_R <- function(signal, ranks, table, groups) {
  n <- nrow(signal)
  out <- matrix(0, n, ncol(signal))
  ord <- order(ranks)
  members <- split(seq_len(n), groups)
  for (deleted in seq_len(n)) {
    kept <- ord[ord != deleted]
    weights <- numeric(n)
    for (ids in members) {
      ids <- ids[ids != deleted]
      if (!length(ids)) next
      remaining <- ranks[ids] - (ranks[ids] > ranks[deleted])
      weights[ids] <- .dnn_duplicate_block_mean(table[remaining])
    }
    out[deleted, ] <- .compensated_crossprod_R(signal[kept, , drop = FALSE], matrix(weights[kept], ncol = 1L))
  }
  out
}

.dnn_block_loo_R <- function(signal, ranks, table, groups, direct = FALSE) {
  n <- nrow(signal)
  ord <- order(ranks)
  ranked_groups <- groups[ord]
  runs <- rle(ranked_groups)
  ordered_members <- !any(diff(ord)[ranked_groups[-n] == ranked_groups[-1L]] < 0L)
  if (direct || anyDuplicated(runs$values) || !ordered_members) {
    return(.dnn_block_loo_direct_R(signal, ranks, table, groups))
  }
  sizes <- runs$lengths
  ends <- cumsum(sizes)
  starts <- c(1L, utils::head(ends, -1L) + 1L)
  padded <- c(0, table, 0)
  before <- after <- survivor <- numeric(length(sizes))
  for (g in seq_along(sizes)) {
    a <- starts[g]
    b <- ends[g]
    before[g] <- .dnn_duplicate_block_mean(padded[seq.int(a + 1L, b + 1L)])
    after[g] <- .dnn_duplicate_block_mean(padded[seq.int(a, b)])
    if (sizes[g] > 1L) survivor[g] <- .dnn_duplicate_block_mean(table[seq.int(a, b - 1L)])
  }
  before <- rep(before, sizes)
  after <- rep(after, sizes)
  survivor <- rep(survivor, sizes)
  reset_left <- reset_right <- rep(FALSE, n)
  reset_left[starts] <- TRUE
  reset_right[ends] <- TRUE
  out <- matrix(0, n, ncol(signal))
  for (component in seq_len(ncol(signal))) {
    z <- signal[ord, component]
    prefix <- .dnn_block_states_R(z * before)
    suffix <- lapply(.dnn_block_states_R(rev(z * after)), rev)
    local_left <- .dnn_block_states_R(z * survivor, reset_left)
    local_right <- lapply(.dnn_block_states_R(rev(z * survivor), rev(reset_right)), rev)
    # The four states stay split until the final merge. Rounding each state
    # first would again lose a finite survivor after cancellation.
    states <- list(rep(prefix$high[starts], sizes), rep(suffix$high[ends], sizes),
      local_left$high, local_right$high,
      rep(prefix$low[starts], sizes), rep(suffix$low[ends], sizes),
      local_left$low, local_right$low)
    for (a in seq.int(1L, n, by = 2048L)) {
      ids <- seq.int(a, min(n, a + 2047L))
      high <- low <- numeric(length(ids))
      for (state in states) {
        updated <- .compensated_add_R(high, low, state[ids])
        high <- updated$high
        low <- updated$low
      }
      out[ord[ids], component] <- high + low
    }
  }
  out
}

.dnn_block_loo <- function(signal, ranks, table, groups, direct = FALSE) {
  groups <- match(groups, unique(groups))
  cpp <- .rdistnn_compiled_fun("rustat_block_rank_loo")
  if (!is.null(cpp)) {
    return(matrix(cpp(matrix(as.integer(ranks), ncol = 1L), signal, table,
                      as.integer(groups), direct), nrow(signal), ncol(signal)))
  }
  .dnn_block_loo_R(signal, ranks, table, groups, direct)
}

# A dynamic scope lives only for the active public call. Nested weight builds,
# candidate/deletion loops and CLATE delegation share its single warning.
.dnn_duplicate_scope <- function(X, groups = .dnn_duplicate_groups(X)) {
  scoped <- any(vapply(sys.frames(), function(frame) {
    exists("..dnn_duplicate_warning_scope", envir = frame, inherits = FALSE)
  }, logical(1L)))
  assign("..dnn_duplicate_warning_scope", TRUE, envir = parent.frame())
  if (!scoped && anyDuplicated(groups)) warning(paste0(
    "Identical covariate rows detected; weights are averaged within duplicate-covariate blocks, ",
    "including jackknife deletions. This permutation-invariant modification is not covered ",
    "by the paper's theory."), call. = FALSE)
  invisible(groups)
}

.dnn_check_tie_model <- function(model) {
  if (is.list(model) && is.matrix(model$x_mat) &&
      !identical(model$tie_policy, .dnn_tie_policy()) &&
      anyDuplicated(.dnn_duplicate_groups(model$x_mat))) {
    stop(paste0("This retained fit predates duplicate-covariate block averaging. ",
                "Refit with the current implementation before reporting or reusing it."),
         call. = FALSE)
  }
  invisible(NULL)
}

.dnn_reuse_scope <- function(object) {
  if (!is.list(object)) return(invisible(NULL))
  .dnn_check_tie_model(object$.model)
  if (!is.null(object$.model$x_mat)) {
    # Keep the marker in the calling public method, not this helper's frame.
    scoped <- any(vapply(sys.frames(), function(frame) {
      exists("..dnn_duplicate_warning_scope", envir = frame, inherits = FALSE)
    }, logical(1L)))
    if (!scoped) .dnn_duplicate_scope(object$.model$x_mat)
    assign("..dnn_duplicate_warning_scope", TRUE, envir = parent.frame())
  }
  invisible(NULL)
}
