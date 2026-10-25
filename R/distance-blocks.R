# Exact represented-distance groups, separately for every reporting target.
# The policy token versions retained scientific state; old fits require refitting.
.dnn_tie_policy <- function() "exact_distance_block_averaging_v1"

.dnn_prepare_distance_groups <- function(dx) {
  if (is.null(dx)) return(list(groups = NULL, tied_targets = logical()))
  dx <- as.matrix(dx)
  tied <- logical(ncol(dx))
  groups <- matrix(vapply(seq_len(ncol(dx)), function(j) {
    values <- unique(dx[, j])
    tied[j] <<- length(values) < nrow(dx)
    match(dx[, j], values)
  }, integer(nrow(dx))), nrow(dx), ncol(dx))
  list(groups = groups, tied_targets = tied)
}

.dnn_remember_group_facts <- function(groups, tied_targets, scope = .dnn_reporting_scope()) {
  facts <- list(tied_targets = tied_targets)
  if (!is.null(scope) && !is.null(groups)) {
    entry <- list(groups = groups, facts = facts)
    primary <- scope$group_facts
    rows <- function(x) if (is.matrix(x)) nrow(x) else length(x)
    # Deletion refits must not evict full-sample preparation. Retain at most
    # two geometry references, independent of the number of deletions.
    if (is.null(primary) || rows(groups) >= rows(primary$groups)) {
      scope$group_scratch <- primary
      scope$group_facts <- entry
    } else {
      scope$group_scratch <- entry
    }
  }
  facts
}


# Invert validated rank permutations in linear work, preserving target shape.
.dnn_orders_from_ranks <- function(ranks) {
  n <- nrow(ranks)
  matrix(vapply(seq_len(ncol(ranks)), function(j) {
    ord <- integer(n)
    ord[ranks[, j]] <- seq_len(n)
    ord
  }, integer(n)), n, ncol(ranks))
}

.dnn_group_column <- function(groups, j) {
  if (is.null(groups)) NULL else if (is.matrix(groups)) groups[, j] else groups
}

.dnn_group_subset <- function(groups, keep) {
  if (is.null(groups)) NULL else if (is.matrix(groups)) groups[keep, , drop = FALSE] else groups[keep]
}

.dnn_compute_group_facts <- function(groups) {
  list(tied_targets = if (is.null(groups)) logical() else if (!is.matrix(groups)) {
    anyDuplicated(groups) > 0L
  } else {
    vapply(seq_len(ncol(groups)), function(j) anyDuplicated(groups[, j]) > 0L, logical(1L))
  })
}

.dnn_group_facts <- function(groups) {
  if (is.null(groups)) return(list(tied_targets = logical()))
  scope <- .dnn_reporting_scope()
  cached <- if (!is.null(scope)) scope$group_facts else NULL
  if (!is.null(cached) && identical(groups, cached$groups)) return(cached$facts)
  scratch <- if (!is.null(scope)) scope$group_scratch else NULL
  if (!is.null(scratch) && identical(groups, scratch$groups)) return(scratch$facts)
  facts <- .dnn_compute_group_facts(groups)
  .dnn_remember_group_facts(groups, facts$tied_targets, scope)
}

.dnn_has_ties <- function(groups, facts = NULL) {
  if (is.null(facts)) facts <- .dnn_group_facts(groups)
  any(facts$tied_targets)
}

.dnn_average_distance_weights <- function(W, groups, facts = NULL) {
  if (is.null(facts)) facts <- .dnn_group_facts(groups)
  tied <- facts$tied_targets
  if (!any(tied)) return(W)
  for (j in which(rep_len(tied, ncol(W)))) {
    g <- .dnn_group_column(groups, j)
    repeated <- duplicated(g) | duplicated(g, fromLast = TRUE)
    for (ids in split(which(repeated), g[repeated])) {
      W[ids, j] <- .dnn_distance_block_mean(W[ids, j])
    }
  }
  W
}

# Sum first to retain subnormal contributions; divide first only when needed
# to keep an otherwise finite block mean inside the floating-point range.
.dnn_distance_block_mean <- function(values) {
  total <- .compensated_sum(values)
  if (is.finite(total)) return(total / length(values))
  .compensated_sum(values / length(values))
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
      weights[ids] <- .dnn_distance_block_mean(table[remaining])
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
    before[g] <- .dnn_distance_block_mean(padded[seq.int(a + 1L, b + 1L)])
    after[g] <- .dnn_distance_block_mean(padded[seq.int(a, b)])
    if (sizes[g] > 1L) survivor[g] <- .dnn_distance_block_mean(table[seq.int(a, b - 1L)])
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

# One mutable warning scope per public call. An early geometry-free scope or
# untied batch must not mark a later tied batch as already warned.
.dnn_distance_scope <- function(groups = NULL, facts = NULL) {
  scope <- .dnn_reporting_scope()
  if (is.null(scope)) {
    scope <- new.env(parent = emptyenv())
    scope$emitted <- character()
  }
  assign("..dnn_distance_warning_scope", scope, envir = parent.frame())
  if (!"distance_ties" %in% scope$emitted && .dnn_has_ties(groups, facts)) {
    .dnn_warn("distance_ties", paste0(
      "Exact distance ties detected; weights are averaged within distance blocks, ",
      "including jackknife deletions. The accompanying theory assumes a continuous distance distribution ",
      "and does not establish inference under persistent ties."))
  }
  invisible(groups)
}

.dnn_check_tie_policy <- function(policy) {
  if (!identical(policy, .dnn_tie_policy())) {
    .dnn_abort("state_invalid", paste0("This result lacks the current exact-distance localization policy. ",
                "Refit with the current implementation before reporting or reusing it."))
  }
  invisible(NULL)
}

.dnn_check_retained_distance <- function(spec, field) {
  if (!is.list(spec) || !inherits(spec, "rdistnn_distance")) {
    .dnn_abort("state_invalid", sprintf(
      "Current fitted `%s` lacks current geometry identity. Refit with the current estimator.", field),
      argument = field)
  }
  slots <- c("name", "eval", "eval_batch", "setup")
  counts <- vapply(slots, function(slot) sum(names(spec) == slot, na.rm = TRUE), integer(1L))
  if (counts[["name"]] != 1L || any(counts > 1L)) {
    .dnn_abort("state_invalid", sprintf(
      "Current fitted `%s` has missing or duplicated geometry fields. Refit with the current estimator.", field),
      argument = field)
  }
  # Reject a missing exact slot only when downstream `$` would resolve an
  # unknown alias. Known slots and metadata beside exact fields stay valid.
  matched <- pmatch(slots[counts == 0L], names(spec), duplicates.ok = TRUE)
  if (any(!is.na(matched) & !(names(spec)[matched] %in% slots))) {
    .dnn_abort("state_invalid", sprintf(
      "Current fitted `%s` has noncanonical geometry field aliases. Refit with the current estimator.", field),
      argument = field)
  }
  # Validate an exact-field view through the existing structural owner. Its
  # canonical return is discarded; retained specifications are never repaired.
  exact <- stats::setNames(lapply(slots, function(slot) spec[[slot, exact = TRUE]]), slots)
  .dnn_validate_retained({
    if (is.null(exact$name)) {
      .abort_bad_arg(paste0(field, "$name"), "must retain a non-empty geometry name")
    }
    .distance_validate_spec(exact)
  })
  invisible(NULL)
}

.dnn_check_tie_model <- function(model) {
  scope <- .dnn_reporting_scope()
  if (!is.null(scope) && !is.null(scope$validated_model) &&
      identical(model, scope$validated_model)) {
    .dnn_remember_group_facts(model$groups, scope$validated_model_facts$tied_targets, scope)
    return(invisible(NULL))
  }
  required <- c("x_mat", "tie_policy", "distance", "ctx", "inference_distance",
    "inference_ctx", "groups", "dx", "ranks")
  absent <- if (is.list(model)) setdiff(required, names(model)) else required
  if (length(absent)) {
    .dnn_abort("state_invalid", sprintf(
      "Current fitted `object$.model$%s` is missing. Refit with the current estimator.", absent[[1L]]),
      argument = paste0("object$.model$", absent[[1L]]))
  }
  if (!is.list(model) || !is.matrix(model$x_mat) || !is.numeric(model$x_mat) ||
      is.complex(model$x_mat) || any(!is.finite(model$x_mat)) ||
      nrow(model$x_mat) < 1L || ncol(model$x_mat) < 1L) {
    .dnn_abort("state_invalid", paste0("This retained fit lacks current localization model state. ",
                "Refit with the current implementation."))
  }
  .dnn_check_tie_policy(model$tie_policy)
  for (field in c("distance", "inference_distance")) {
    .dnn_check_retained_distance(model[[field, exact = TRUE]], paste0("object$.model$", field))
  }
  n <- nrow(model$x_mat)
  groups <- model$groups
  dx <- model$dx
  ranks <- model$ranks
  if (!is.matrix(groups) || !is.numeric(groups) || is.complex(groups) || anyNA(groups) ||
      any(!is.finite(groups)) || any(groups < 1L | groups > n) ||
      any(groups != as.integer(groups)) || nrow(groups) != n || ncol(groups) < 1L) {
    .dnn_abort("state_invalid", paste0("This retained fit lacks valid target-specific distance groups. ",
                "Refit with the current implementation."))
  }
  if (!is.matrix(dx) || !is.numeric(dx) || is.complex(dx) || any(!is.finite(dx)) ||
      !identical(dim(dx), dim(groups))) {
    .dnn_abort("state_invalid", paste0("This retained fit has inconsistent target-specific distance groups. ",
                "Refit with the current implementation."))
  }
  if (!is.matrix(ranks) || !is.numeric(ranks) || is.complex(ranks) ||
      !identical(dim(ranks), dim(dx)) || any(!is.finite(ranks)) ||
      any(ranks < 1L | ranks > n) || any(ranks != as.integer(ranks)) ||
      any(vapply(seq_len(ncol(ranks)), function(j) anyDuplicated(ranks[, j]) > 0L, logical(1L)))) {
    .dnn_abort("state_invalid", paste0("This retained fit lacks valid cached rank permutations. ",
                "Refit with the current implementation."))
  }
  ord <- .dnn_orders_from_ranks(ranks)
  tied <- logical(ncol(dx))
  for (j in seq_len(ncol(dx))) {
    left <- ord[-n, j]
    right <- ord[-1L, j]
    if (any(dx[left, j] > dx[right, j]) ||
        any(dx[left, j] == dx[right, j] & left > right)) {
      .dnn_abort("state_invalid",
        paste0("This retained fit has cached ranks inconsistent with distance-then-index order. ",
                  "Refit with the current implementation."))
    }
    # Validate equality against the already ordered distances. This supplies
    # tie facts without rebuilding/hashing distance groups for every panel.
    ordered_groups <- groups[ord[, j], j]
    starts <- c(TRUE, dx[right, j] != dx[left, j])
    group_starts <- c(TRUE, ordered_groups[-1L] != ordered_groups[-n])
    count <- sum(starts)
    first <- match(seq_len(count), groups[, j])
    if (!identical(unname(starts), unname(group_starts)) || anyDuplicated(ordered_groups[starts]) ||
        any(groups[, j] > count) || anyNA(first) || any(diff(first) <= 0L)) {
      .dnn_abort("state_invalid", paste0("This retained fit has inconsistent target-specific distance groups. ",
        "Refit with the current implementation."))
    }
    tied[j] <- count < n
  }
  .dnn_remember_group_facts(groups, tied, scope)
  if (!is.null(scope)) {
    scope$validated_model <- model
    scope$validated_model_facts <- list(tied_targets = tied)
  }
  invisible(NULL)
}

.dnn_reuse_scope <- function(object, groups = object$.model$groups) {
  if (!is.list(object)) return(invisible(NULL))
  .dnn_distance_scope()
  .dnn_check_report_state(object)
  .dnn_distance_scope(groups = groups)
  # Install the same environment in the public caller after the helper exits.
  assign("..dnn_distance_warning_scope", get("..dnn_distance_warning_scope", inherits = FALSE),
         envir = parent.frame())
  invisible(NULL)
}
