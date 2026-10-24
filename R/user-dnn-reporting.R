# Shared pointwise normal critical value and descriptive localization labels.
.dnn_normal_critical <- function(level) {
  level <- .validate_scalar_prob(level, "level")
  if (level < 1e-4) {
    # sqrt(2)*inverf(level), through its cubic term. Below this threshold
    # the omitted relative remainder is <1.5e-17 (its leading coefficient
    # is 7*pi^2/480); unlike 1-level, this preserves tiny central levels.
    return((sqrt(pi / 2) * level) * (1 + (pi / 12) * level^2))
  }
  stats::qnorm((1 - level) / 2, lower.tail = FALSE)
}

.dnn_normal_critical_alpha <- function(alpha) {
  alpha <- .validate_scalar_prob(alpha, "alpha")
  if (alpha > 1 - 1e-4) return(.dnn_normal_critical(1 - alpha))
  stats::qnorm(log(alpha) - log(2), lower.tail = FALSE, log.p = TRUE)
}

.dnn_localization_label <- function(meta) {
  orders <- if (!is.null(meta$s1) && !is.null(meta$s2)) {
    c(s1 = meta$s1, s2 = meta$s2)
  } else if (!is.null(meta$s)) {
    c(s = meta$s)
  } else {
    return("")
  }
  n <- meta$n
  valid_n <- is.numeric(n) && length(n) == 1L && is.finite(n) && n > 0
  paste(vapply(names(orders), function(name) {
    s <- orders[[name]]
    fraction <- if (valid_n) format(signif(s / n, 3L), trim = TRUE) else "unavailable"
    sprintf("%s = %s, %s/n = %s", name, s, name, fraction)
  }, character(1L)), collapse = "; ")
}

.dnn_print_localization <- function(meta) {
  label <- .dnn_localization_label(meta)
  if (nzchar(label)) cat("  localization: ", label, "\n", sep = "")
  invisible(NULL)
}

.dnn_print_diagnostic_range <- function(x, field, label) {
  if (!field %in% names(x)) return(invisible(NULL))
  values <- x[[field]][is.finite(x[[field]])]
  text <- if (length(values)) {
    paste(format(signif(range(values), 3L), trim = TRUE), collapse = " to ")
  } else {
    "unavailable"
  }
  cat("  ", label, ": ", text, "\n", sep = "")
  invisible(NULL)
}

# Reports and reconstruction must refer to the same fitted computation. Check
# duplicated controls directly, without rebuilding weights or consulting the
# live registry, which may have changed since the fit was created.
.dnn_check_fitted_controls <- function(object, invalid, require_fields) {
  meta <- object$meta
  model <- object$.model
  tdnn <- inherits(object, "rdistnn_tdnn_np_result")
  paired <- inherits(object, "rdistnn_affine_ratio_result")
  orders <- if (tdnn) c("s1", "s2", "dim") else c("k", "s")
  require_fields(c("method", "distance", "tie_rule", orders), meta, "meta$")
  require_fields(c("method", "normalize", "k", orders), model, ".model$")
  for (location in c("meta", ".model")) {
    value <- if (location == "meta") meta$method else model$method
    if (!is.character(value) || !is.null(dim(value)) || length(value) != 1L ||
        is.na(value) || !value %in% c("exact", "asymp")) {
      invalid(paste0(location, "$method"), "must retain 'exact' or 'asymp'")
    }
  }
  if (!identical(meta$method, model$method)) {
    invalid("meta$method/.model$method", "contains inconsistent retained weight methods")
  }
  positive_integer <- function(value, field) {
    if (!is.numeric(value) || is.complex(value) || !is.null(dim(value)) ||
        length(value) != 1L || !is.finite(value) || value < 1 ||
        value > .Machine$integer.max || value != trunc(value)) {
      invalid(field, "must retain a positive integer scalar")
    }
    value
  }
  for (field in orders) {
    public <- positive_integer(meta[[field]], paste0("meta$", field))
    retained <- positive_integer(model[[field]], paste0(".model$", field))
    if (public != retained) {
      invalid(paste0("meta$", field, "/.model$", field), "contains inconsistent fitted controls")
    }
  }
  if (tdnn && positive_integer(model$k, ".model$k") != 1L) {
    invalid(".model$k", "must retain TDNN's fixed k = 1 control")
  }
  logical_scalar <- function(value, field) {
    if (!is.logical(value) || !is.null(dim(value)) || length(value) != 1L || is.na(value)) {
      invalid(field, "must retain TRUE or FALSE")
    }
    value
  }
  normalize <- logical_scalar(model$normalize, ".model$normalize")
  # Only paired ratios duplicate normalization in their current metadata.
  if (paired) {
    require_fields("normalize", meta, "meta$")
    if (logical_scalar(meta$normalize, "meta$normalize") != normalize) {
      invalid("meta$normalize/.model$normalize", "contains inconsistent fitted normalization")
    }
  }
  label <- function(value, field) {
    if (!is.character(value) || !is.null(dim(value)) || length(value) != 1L ||
        is.na(value) || !nzchar(value)) invalid(field, "must retain one nonempty geometry label")
    value
  }
  label(meta$distance, "meta$distance")
  frozen <- label(model$distance$name, ".model$distance$name")
  original <- label(model$inference_distance$name, ".model$inference_distance$name")
  if (!identical(frozen, original)) {
    invalid(".model$distance$name/.model$inference_distance$name", "contains inconsistent fitted geometry names")
  }
  # meta$distance is the caller's selector, which can be a registry alias;
  # the two retained specs share a name even when that selector differs.
  if (!identical(meta$tie_rule, model$tie_policy)) {
    invalid("meta$tie_rule/.model$tie_policy", "contains inconsistent localization policies")
  }
  invisible(NULL)
}

# Validate cached public arrays before reports or reuse consume them. This
# checks the retained schema, without recomputing estimates or metric values.
.dnn_check_report_state <- function(object) {
  # Prediction validators reuse this operation's validated geometry rather
  # than rechecking cached ranks and distance groups independently.
  .dnn_distance_scope()
  scope <- .dnn_reporting_scope()
  if (!is.null(scope) && !is.null(scope$validated_report) &&
      identical(object, scope$validated_report)) return(invisible(NULL))
  .dnn_check_tie_model(object[[".model", exact = TRUE]])
  if (!inherits(object, c("rdistnn_np_result", "rdistnn_cate_result", "rdistnn_affine_ratio_result",
    "rdistnn_plqr_result", "rdistnn_gmm_result"))) {
    return(invisible(NULL))
  }
  invalid <- function(field, reason) {
    .dnn_abort(
      "state_invalid",
      sprintf("Current fitted `object$%s` %s. Refit with the current estimator.", field, reason),
      argument = paste0("object$", field)
    )
  }
  require_fields <- function(fields, value = object, prefix = "") {
    absent <- setdiff(fields, names(value))
    if (length(absent)) invalid(paste0(prefix, absent[[1L]]), "is missing")
  }
  check_array <- function(value, field, shape, allow_na = FALSE, nonnegative = FALSE) {
    shaped <- if (length(shape) == 1L) is.null(dim(value)) && length(value) == shape else
      is.matrix(value) && identical(dim(value), as.integer(shape))
    if (!is.numeric(value) || is.complex(value) || !shaped) {
      invalid(field, paste0("must be a real numeric array with dimensions ", paste(shape, collapse = " by ")))
    }
    admissible <- is.finite(value) | (allow_na & is.na(value) & !is.nan(value))
    if (!all(admissible) || (nonnegative && any(value[!is.na(value)] < 0))) {
      invalid(field, "contains unavailable, nonfinite or negative values outside its reporting contract")
    }
    invisible(NULL)
  }
  check_covariance <- function(value, field, components) {
    if (!is.list(value) || length(value) != m) invalid(field, "must retain one covariance matrix per target")
    for (j in seq_len(m)) {
      check_array(value[[j]], field, c(components, components), allow_na = TRUE)
      diagonal <- diag(value[[j]])
      if (any(diagonal[!is.na(diagonal)] < 0)) invalid(field, "contains a negative variance")
    }
    invisible(NULL)
  }
  require_fields(c("meta", "x0"))
  model <- object$.model
  n <- nrow(model$x_mat)
  p <- ncol(model$x_mat)
  m <- ncol(model$groups)
  if (!is.list(object$meta)) invalid("meta", "must retain the current metadata")
  for (field in c("n", "p", "m")) {
    value <- object$meta[[field]]
    if (!is.numeric(value) || is.complex(value) || length(value) != 1L ||
        !is.finite(value) || value < 1 || value > .Machine$integer.max || value != trunc(value)) {
      invalid(paste0("meta$", field), "must be a positive integer scalar")
    }
  }
  if (object$meta$n != n || object$meta$p != p || object$meta$m != m) {
    invalid("meta$n/p/m", "is inconsistent with the retained observation and target dimensions")
  }
  require_fields(c("target_notices", "n_tied_targets", "n_weak_targets"), object$meta, "meta$")
  notices <- object$meta$target_notices
  if (!is.data.frame(notices) || nrow(notices) != m ||
      !all(c("target_id", "distance_ties", "weight_concentration", "weak_localization") %in% names(notices)) ||
      !identical(notices$target_id, seq_len(m))) {
    invalid("meta$target_notices", "must retain aligned target diagnostic facts")
  }
  for (field in c("distance_ties", "weak_localization")) {
    if (!is.logical(notices[[field]]) || !is.null(dim(notices[[field]])) ||
        length(notices[[field]]) != m || anyNA(notices[[field]])) {
      invalid(paste0("meta$target_notices$", field), "must retain one logical flag per target")
    }
  }
  concentration <- notices$weight_concentration
  if (!is.numeric(concentration) || is.complex(concentration) || !is.null(dim(concentration)) ||
      length(concentration) != m || anyNA(concentration) || any(concentration < 0)) {
    invalid("meta$target_notices$weight_concentration", "must retain nonnegative concentration values")
  }
  if (!identical(notices$distance_ties, unname(.dnn_group_facts(model$groups)$tied_targets)) ||
      !identical(notices$weak_localization, unname(n >= 2L & concentration >= .9 * n))) {
    invalid("meta$target_notices", "is inconsistent with the validated distance groups or concentration")
  }
  for (kind in c("tied", "weak")) {
    field <- paste0("n_", kind, "_targets")
    flag <- if (kind == "tied") "distance_ties" else "weak_localization"
    value <- object$meta[[field]]
    if (!is.numeric(value) || length(value) != 1L || !is.finite(value) ||
        value != sum(notices[[flag]])) invalid(paste0("meta$", field), "does not match the target flags")
  }
  if (!is.matrix(object$x0) || !identical(dim(object$x0), as.integer(c(m, p)))) {
    invalid("x0", "must match `object$meta$m` target rows and `object$meta$p` covariate columns")
  }
  check_array(object$x0, "x0", c(m, p))
  require_fields("x_cols", model, ".model$")
  x_cols <- model[["x_cols", exact = TRUE]]
  # Repeated selectors and matrix-valued columns can expand differently
  # from canonical matrix columns; preserve the original selector grammar.
  if (!is.null(x_cols) && (!is.character(x_cols) || !is.null(dim(x_cols)) ||
      length(x_cols) < 1L || anyNA(x_cols) || any(!nzchar(x_cols)))) {
    invalid(".model$x_cols", "must retain NULL or nonempty character column selectors")
  }
  .dnn_check_fitted_controls(object, invalid, require_fields)
  has_grid <- !is.null(object[["tuning_grid", exact = TRUE]])
  has_tuning <- "tuning" %in% names(object$meta)
  if (has_tuning || has_grid) {
    if (sum(names(object$meta) == "tuning", na.rm = TRUE) != 1L) {
      invalid("meta$tuning", "must be retained exactly once")
    }
    .dnn_validate_retained_tuning(object$meta[["tuning", exact = TRUE]], object$meta$s)
  }
  if (has_grid) {
    grid <- .candidate_grid_payload(object)
    if (!identical(grid$meta$jackknife, object$meta$jackknife)) {
      invalid("tuning_grid$meta$jackknife", "is inconsistent with the fitted uncertainty policy")
    }
    anchor <- which(grid$candidates$is_default)
    if (!isTRUE(grid$candidates$feasible[anchor]) ||
        !isTRUE(grid$candidates$status[anchor] == "ok")) {
      invalid("tuning_grid$candidates$is_default",
        "must identify a feasible, successfully evaluated retained default")
    }
  }
  if (inherits(object, "rdistnn_gmm_result")) {
    .gmm_check_report_arrays(object, n, m, invalid, require_fields, check_array)
    if (!is.null(scope)) scope$validated_report <- object
    return(invisible(NULL))
  }
  if (inherits(object, "rdistnn_np_result")) {
    require_fields(c("y_vec", "kind"), model, ".model$")
    .np_prediction_state(object)
    check_array(model[["y_vec", exact = TRUE]], ".model$y_vec", n)
    if (inherits(object, "rdistnn_tdnn_np_result") &&
        (!identical(model$kind, "tdnn") || model$k != 1L)) {
      invalid(".model$kind/k", "must retain the fitted TDNN rule with k = 1")
    }
  }
  if (inherits(object, "rdistnn_cate_result")) {
    .cate_result_meta(object$meta)
    require_fields(c("psi", "t_vec", "e"), model, ".model$")
    check_array(model[["psi", exact = TRUE]], ".model$psi", n)
    .cate_validate_prediction_model(model)
    aipw <- identical(object$meta[["signal_source", exact = TRUE]], "aipw")
    for (field in c("t_vec", "e")) {
      value <- model[[field, exact = TRUE]]
      if (is.null(value)) {
        if (aipw) invalid(paste0(".model$", field), "must retain the fitted AIPW overlap inputs")
      } else {
        check_array(value, paste0(".model$", field), n)
        if ((field == "t_vec" && any(!value %in% c(0, 1))) ||
            (field == "e" && any(value <= 0 | value >= 1))) {
          invalid(paste0(".model$", field), "must retain aligned binary treatment or interior propensity values")
        }
      }
    }
  }
  if (inherits(object, "rdistnn_clate_result")) {
    require_fields("iv", model, ".model$")
    iv <- model[["iv", exact = TRUE]]
    if (!is.list(iv)) invalid(".model$iv", "must retain the fitted instrument record")
    require_fields(c("A", "D", "one_sided", "propensity", "training"), iv, ".model$iv$")
    arm <- iv[["A", exact = TRUE]]
    check_array(arm, ".model$iv$A", n)
    if (any(!arm %in% c(0, 1))) invalid(".model$iv$A", "must retain aligned binary instrument values")
    check_array(iv$D, ".model$iv$D", n)
    if (any(!iv$D %in% c(0, 1))) invalid(".model$iv$D", "must retain aligned binary treatment values")
    if (!is.logical(iv$one_sided) || length(iv$one_sided) != 1L || is.na(iv$one_sided) ||
        (iv$one_sided && any(iv$D > arm))) {
      invalid(".model$iv$one_sided", "must retain the fitted one-sided declaration and compatible treatment values")
    }
    propensity <- iv[["propensity", exact = TRUE]]
    if (!is.list(propensity)) invalid(".model$iv$propensity", "must retain fitted propensity inputs")
    require_fields(c("raw", "used", "clipped"), propensity, ".model$iv$propensity$")
    check_array(propensity$raw, ".model$iv$propensity$raw", n)
    if (any(propensity$raw < 0 | propensity$raw > 1)) {
      invalid(".model$iv$propensity$raw", "must retain propensity values in [0, 1]")
    }
    used <- propensity[["used", exact = TRUE]]
    check_array(used, ".model$iv$propensity$used", n)
    if (any(used <= 0 | used >= 1)) invalid(".model$iv$propensity$used", "must retain interior propensity values")
    clipped <- propensity[["clipped", exact = TRUE]]
    if (!is.logical(clipped) || !is.null(dim(clipped)) || length(clipped) != n || anyNA(clipped)) {
      invalid(".model$iv$propensity$clipped", "must retain aligned logical clipping flags")
    }
    if (!identical(unname(clipped), unname(propensity$raw != used))) {
      invalid(".model$iv$propensity$clipped", "is inconsistent with the raw and used propensity values")
    }
    training <- iv[["training", exact = TRUE]]
    if (!is.character(training) || length(training) != 1L || is.na(training) ||
        !training %in% c("supplied", "oracle", "full_sample", "oob", "cross_fit")) {
      invalid(".model$iv$training", "must retain the fitted training declaration")
    }
  }
  if (inherits(object, "rdistnn_plqr_result")) {
    require_fields(c("theta", "moment", "theta_eval", "theta_se", "moment_se"))
    q <- ncol(model$d_mat)
    if (is.null(q) || q < 1L) invalid(".model$d_mat", "must retain the coefficient dimension")
    check_array(object$theta_eval, "theta_eval", q)
    check_array(object$theta, "theta", c(q, m))
    check_array(object$moment, "moment", c(q, m))
    for (field in c("theta_se", "moment_se")) {
      if (!is.null(object[[field]])) check_array(object[[field]], field, c(q, m), allow_na = TRUE, nonnegative = TRUE)
    }
    if (!is.null(object$theta_se)) {
      require_fields("theta_vcov")
      check_covariance(object$theta_vcov, "theta_vcov", q)
    }
    .plqr_check_report_model(object, n, m, invalid, require_fields, check_array, check_covariance)
    if (!is.null(scope)) scope$validated_report <- object
    return(invisible(NULL))
  }
  paired <- inherits(object, "rdistnn_affine_ratio_result")
  estimate <- if (paired) "theta" else if (inherits(object, "rdistnn_cate_result")) "tau" else "mu"
  require_fields(c(estimate, "se"))
  check_array(object[[estimate]], estimate, m, allow_na = paired)
  jackknife <- object$meta[["jackknife", exact = TRUE]]
  if (!is.character(jackknife) || length(jackknife) != 1L || is.na(jackknife) ||
      !jackknife %in% c("none", "fast", "loo")) {
    invalid("meta$jackknife", "must retain 'none', 'fast' or 'loo'")
  }
  if ((jackknife == "none") != is.null(object$se)) {
    invalid("se/meta$jackknife", "contains inconsistent uncertainty availability")
  }
  if (!is.null(object$se)) {
    check_array(object$se, "se", m, allow_na = paired, nonnegative = TRUE)
  }
  if (paired) {
    require_fields(c("signal", "row_ids", "denominator_tol"), model, ".model$")
    require_fields("row_ids")
    check_array(model$signal, ".model$signal", c(n, 2L))
    check_array(model$denominator_tol, ".model$denominator_tol", 1L, nonnegative = TRUE)
    .dnn_validate_retained({
      .affine_prediction_state(object)
      k <- .validate_scalar_int(model$k, "object$.model$k")
      s <- .validate_scalar_int(model$s, "object$.model$s")
      if (k > s || s > n) .abort_bad_arg("object$.model$k/s", "must satisfy k <= s <= n")
      if (is.null(model$row_ids) || is.null(object$row_ids)) {
        invalid(".model$row_ids/row_ids", "must retain the fitted observation identities")
      }
      .affine_row_ids(model$row_ids, n)
      .affine_row_ids(object$row_ids, n)
    })
    if (!identical(object$row_ids, model$row_ids)) {
      invalid(".model$row_ids/row_ids", "must retain the same observation identities in training-row order")
    }
    require_fields(c("components", "status", "component_vcov", "theta_vcov", "theta_omega"))
    check_array(object$components, "components", c(m, 2L), allow_na = TRUE)
    covariance <- object$meta[["covariance", exact = TRUE]]
    if (!is.character(covariance) || length(covariance) != 1L || is.na(covariance) ||
        !covariance %in% c("none", "paired_frozen_score")) {
      invalid("meta$covariance", "must retain 'none' or 'paired_frozen_score'")
    }
    uncertainty <- c("se", "component_vcov", "theta_vcov", "theta_omega")
    disabled <- vapply(object[uncertainty], is.null, logical(1L))
    if ((covariance == "none" && !all(disabled)) ||
        (covariance == "paired_frozen_score" && any(disabled))) {
      invalid("se/covariance", "is inconsistent with the retained covariance mode")
    }
    for (field in c("component_vcov", "theta_vcov", "theta_omega")) {
      if (!is.null(object[[field]])) check_covariance(object[[field]], field,
        if (field == "component_vcov") 2L else 1L)
    }
    status <- object$status
    if (!is.data.frame(status) || nrow(status) != m ||
        !"target_id" %in% names(status) ||
        !is.numeric(status$target_id) || anyNA(status$target_id) ||
        any(status$target_id != seq_len(m))) {
      invalid("status", "must contain aligned target identities and availability flags")
    }
    .affine_validate_retained_status(status, object$theta, object$se,
      covariance == "paired_frozen_score", invalid)
  }
  if (!is.null(scope)) scope$validated_report <- object
  invisible(NULL)
}
