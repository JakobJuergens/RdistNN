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

# Validate cached public arrays before reports or reuse consume them. This
# checks the retained schema, without recomputing estimates or metric values.
.dnn_check_report_state <- function(object) {
  scope <- .dnn_reporting_scope()
  if (!is.null(scope) && !is.null(scope$validated_report) &&
      identical(object, scope$validated_report)) return(invisible(NULL))
  .dnn_check_tie_model(object$.model)
  if (!inherits(object, c("np_result", "cate_result", "affine_ratio_result", "plqr_result"))) {
    return(invisible(NULL))
  }
  invalid <- function(field, reason) {
    .abort_bad_state(sprintf("Current fitted `object$%s` %s. Refit with the current estimator.", field, reason),
      code = "state_invalid")
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
  if (!is.matrix(object$x0) || !identical(dim(object$x0), as.integer(c(m, p)))) {
    invalid("x0", "must match `object$meta$m` target rows and `object$meta$p` covariate columns")
  }
  check_array(object$x0, "x0", c(m, p))
  if (!inherits(object, "plqr_result")) {
    require_fields("x_cols", model, ".model$")
    x_cols <- model[["x_cols", exact = TRUE]]
    # Repeated selectors and matrix-valued columns can expand differently
    # from canonical matrix columns; preserve the original selector grammar.
    if (!is.null(x_cols) && (!is.character(x_cols) || !is.null(dim(x_cols)) ||
        length(x_cols) < 1L || anyNA(x_cols) || any(!nzchar(x_cols)))) {
      invalid(".model$x_cols", "must retain NULL or nonempty character column selectors")
    }
  }
  if (inherits(object, "cate_result")) {
    require_fields(c("psi", "t_vec", "e"), model, ".model$")
    check_array(model[["psi", exact = TRUE]], ".model$psi", n)
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
  if (inherits(object, "clate_result")) {
    require_fields("iv", model, ".model$")
    iv <- model[["iv", exact = TRUE]]
    if (!is.list(iv)) invalid(".model$iv", "must retain the fitted instrument record")
    require_fields(c("A", "propensity", "training"), iv, ".model$iv$")
    arm <- iv[["A", exact = TRUE]]
    check_array(arm, ".model$iv$A", n)
    if (any(!arm %in% c(0, 1))) invalid(".model$iv$A", "must retain aligned binary instrument values")
    propensity <- iv[["propensity", exact = TRUE]]
    if (!is.list(propensity)) invalid(".model$iv$propensity", "must retain fitted propensity inputs")
    require_fields(c("used", "clipped"), propensity, ".model$iv$propensity$")
    used <- propensity[["used", exact = TRUE]]
    check_array(used, ".model$iv$propensity$used", n)
    if (any(used <= 0 | used >= 1)) invalid(".model$iv$propensity$used", "must retain interior propensity values")
    clipped <- propensity[["clipped", exact = TRUE]]
    if (!is.logical(clipped) || !is.null(dim(clipped)) || length(clipped) != n || anyNA(clipped)) {
      invalid(".model$iv$propensity$clipped", "must retain aligned logical clipping flags")
    }
    training <- iv[["training", exact = TRUE]]
    if (!is.character(training) || length(training) != 1L || is.na(training) ||
        !training %in% c("supplied", "oracle", "full_sample", "oob", "cross_fit")) {
      invalid(".model$iv$training", "must retain the fitted training declaration")
    }
  }
  if (inherits(object, "plqr_result")) {
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
    return(invisible(NULL))
  }
  paired <- inherits(object, "affine_ratio_result")
  estimate <- if (paired) "theta" else if (inherits(object, "cate_result")) "tau" else "mu"
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
  if (paired) {
    model_method <- model[["method", exact = TRUE]]
    meta_method <- object$meta[["method", exact = TRUE]]
    for (field in c(".model$method", "meta$method")) {
      method <- if (field == ".model$method") model_method else meta_method
      if (!is.character(method) || length(method) != 1L || is.na(method) ||
          !method %in% c("exact", "asymp")) {
        invalid(field, "must retain 'exact' or 'asymp'")
      }
    }
    if (model_method != meta_method) {
      invalid(".model$method/meta$method", "contains inconsistent retained weight methods")
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
        !all(c("target_id", "ratio_status", "wald_available") %in% names(status)) ||
        !is.numeric(status$target_id) || anyNA(status$target_id) ||
        any(status$target_id != seq_len(m)) || !is.character(status$ratio_status) ||
        anyNA(status$ratio_status) || !is.logical(status$wald_available) || anyNA(status$wald_available)) {
      invalid("status", "must contain aligned target identities and availability flags")
    }
    if (any(is.na(object$theta) & (status$ratio_status == "ok" | status$wald_available))) {
      invalid("theta/status", "contains inconsistent unavailable ratios")
    }
    if (is.null(object$se) && any(status$wald_available)) {
      invalid("se/status", "contains inconsistent Wald availability")
    }
  }
  if (!is.null(object$se)) {
    check_array(object$se, "se", m, allow_na = paired, nonnegative = TRUE)
    if (paired && any(object$status$wald_available & (is.na(object$se) | object$se <= 0))) {
      invalid("se/status", "contains inconsistent Wald availability")
    }
  }
  if (!is.null(scope)) scope$validated_report <- object
  invisible(NULL)
}
