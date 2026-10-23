#' Confidence sets for a ratio of paired component estimates
#'
#' Inverts the pointwise normal-reference moment inequality without dividing
#' by an estimated denominator. This experimental helper accepts the components
#' and paired estimator-scale covariance returned by [affine_ratio_dnn()].
#'
#' @param components Finite numeric m-by-2 matrix, with numerator in the first
#'   column and denominator in the second; m must be positive.
#' @param component_vcov List of m finite symmetric positive semidefinite 2-by-2
#'   covariance matrices on the estimator scale. Singular matrices are allowed.
#' @param level Confidence level strictly between zero and one.
#' @param tol Finite nonnegative relative numerical tolerance. It controls
#'   covariance symmetry/PSD checks and near-zero discriminant diagnostics.
#' @param x A ratio_confidence_sets object.
#' @param digits Printing precision.
#' @param ... Not used.
#'
#' @details For candidate theta, invert
#'   \deqn{(\hat\Delta_Y-\theta\hat\Delta_D)^2
#'     \le z^2(V_{YY}-2\theta V_{YD}+\theta^2V_{DD}).}
#'   Inputs must retain numerator-denominator pairing. Do not divide the supplied
#'   covariance by effective sample size again.
#'
#'   Each set can be bounded, a singleton, two rays, a half-line, all real
#'   numbers or empty. Finite endpoints are included; infinities are open in
#'   printed interval notation. An unresolved numerical result has no intervals
#'   and is distinguished from a mathematically empty set by its type/status.
#'   Near-boundary calculations are reported rather than rounding coefficients
#'   or discriminants to zero. Materially indefinite covariance is rejected;
#'   a negative eigenvalue within tolerance gives an unresolved result.
#'   Deterministic-denominator intervals are computed directly from the
#'   components to preserve small positive uncertainty. Repeated roots of a
#'   rounded polynomial alone do not certify zero uncertainty; such cases are
#'   unresolved unless a component-level zero-variance identity establishes it.
#'
#'   Power-of-two coordinate scaling avoids intermediate overflow. The numerical
#'   record retains original quadratic coefficients where representable, plus
#'   scaled coefficients for t = theta * 2^(-ratio_scale_exponent). Multiplying
#'   that polynomial by 2^quadratic_scale_exponent gives the polynomial of the
#'   coordinate-scaled components. This metadata permits reconstruction when
#'   the original coefficients exceed floating-point range.
#'
#'   At a fixed true candidate, inclusion follows from validity of that
#'   candidate's local score test under its design, learner, bias and positive
#'   variance conditions. This calculation establishes no uniform coverage over
#'   increasingly weak instruments. Algebraic zero-variance branches do not
#'   automatically satisfy a positive-variance test theorem.
#'
#' @return A ratio_confidence_sets list with one record per target. Records
#'   contain target_id, type, intervals (r-by-2 lower/upper matrix), quadratic,
#'   discriminant, numerical_status, level and numerical scaling metadata.
#'   Printing returns its input invisibly.
#' @examples
#' f <- affine_ratio_dnn(
#'   X = matrix(1:4, ncol = 1), x0 = 0,
#'   numerator = c(2, 1, 4, 3), denominator = c(1, -1, 2, 1), s = 2
#' )
#' ratio_confidence_set(f$components, f$component_vcov)
#' # An exactly zero denominator can give two unbounded rays.
#' ratio_confidence_set(matrix(c(3, 0), 1), list(diag(2)))
#' @export
ratio_confidence_set <- function(
  components, component_vcov, level = .95, tol = 64 * .Machine$double.eps
) {
  components <- .validate_matrix(components, "components",
    min_rows = 1L,
    min_cols = 2L, numeric = TRUE, finite = TRUE
  )
  if (ncol(components) != 2L) .abort_bad_arg("components", "must have exactly two columns")
  m <- nrow(components)
  if (!is.list(component_vcov) || length(component_vcov) != m) {
    .abort_bad_arg("component_vcov", sprintf("must be a list of %d paired 2-by-2 matrices", m))
  }
  level <- .validate_scalar_prob(level, "level")
  tol <- .affine_nonnegative_scalar(tol, "tol")
  critical <- .dnn_normal_critical(level)
  out <- lapply(seq_len(m), function(j) {
    checked <- .ratio_validate_covariance(component_vcov[[j]], j, tol)
    V <- checked$covariance
    delta <- components[j, ]
    exponents <- vapply(1:2, function(k) {
      .ratio_scale_exponent(max(abs(delta[k]), sqrt(V[k, k])))
    }, numeric(1L))
    normalized <- .ratio_pow2(delta, -exponents)
    v <- matrix(.ratio_pow2(as.numeric(V), -as.numeric(outer(exponents, exponents, "+"))), 2L)
    abc <- c(
      a = normalized[2L]^2 - critical^2 * v[2L, 2L],
      b = -2 * normalized[1L] * normalized[2L] + 2 * critical^2 * v[1L, 2L],
      c = normalized[1L]^2 - critical^2 * v[1L, 1L]
    )
    names(abc) <- c("a", "b", "c")
    answer <- do.call(.ratio_quadratic_set, c(as.list(abc), list(tol = tol)))
    reason <- NULL
    if (any(V != 0) && answer$scaled_discriminant == 0 &&
      !(delta[1L] == 0 && V[1L, 1L] == 0)) {
      reason <- "coefficient_rounding"
    }
    if (delta[2L] != 0 && V[2L, 2L] == 0 && V[1L, 1L] > 0 && critical > 0) {
      # Expansion can absorb z^2 V_YY into Delta_Y^2 even when the interval
      # width is representable. Compute that width before forming the square.
      center <- normalized[1L] / normalized[2L]
      width <- critical * sqrt(v[1L, 1L]) / abs(normalized[2L])
      bounds <- center + c(-1, 1) * width
      answer$type <- "bounded"
      answer$intervals <- .ratio_intervals(bounds)
      answer$numerical_status <- "ok"
      reason <- if (width == 0 || bounds[1L] == bounds[2L]) "endpoint_precision" else NULL
      direct_discriminant <- 4 * critical^2 * normalized[2L]^2 * v[1L, 1L]
      answer$scaled_discriminant <- .ratio_pow2(
        direct_discriminant, -2 * answer$quadratic_scale_exponent
      )
      if (answer$scaled_discriminant == 0) answer$scaled_discriminant <- NA_real_
    }
    answer$target_id <- j
    answer$level <- level
    answer$covariance_status <- checked$status
    answer$ratio_scale_exponent <- exponents[1L] - exponents[2L]
    original_exponents <- c(2 * exponents[2L], sum(exponents), 2 * exponents[1L])
    original <- .ratio_pow2(abc, original_exponents)
    representable <- is.finite(original) & (original != 0 | abc == 0)
    original[!representable] <- NA_real_
    answer$quadratic <- stats::setNames(original, c("a", "b", "c"))
    answer$quadratic_representable <- all(representable)
    answer$discriminant <- .ratio_pow2(
      answer$scaled_discriminant,
      2 * answer$quadratic_scale_exponent + 2 * sum(exponents)
    )
    if (!is.finite(answer$discriminant) ||
      (answer$discriminant == 0 && answer$scaled_discriminant != 0)) {
      answer$discriminant <- NA_real_
    }
    if (checked$status == "near_psd_boundary") reason <- "near_psd_boundary"
    if (any(V != 0 & v == 0)) reason <- "scaling_underflow"
    if (critical == 0 && any(V != 0)) reason <- "confidence_level_precision"
    if (all(V == 0)) {
      if (delta[2L] != 0) {
        root <- delta[1L] / delta[2L]
        answer$type <- "singleton"
        answer$intervals <- .ratio_intervals(c(root, root))
        answer$numerical_status <- "ok"
        if (!is.finite(root) || (root == 0 && delta[1L] != 0)) reason <- "endpoint_range"
      }
    } else if (nrow(answer$intervals)) {
      endpoints <- answer$intervals
      finite <- is.finite(endpoints)
      endpoints[finite] <- .ratio_pow2(endpoints[finite], answer$ratio_scale_exponent)
      if (any(!is.finite(endpoints[finite])) ||
        any(endpoints[finite] == 0 & answer$intervals[finite] != 0)) {
        reason <- "endpoint_range"
      } else {
        answer$intervals <- endpoints
      }
    }
    if (!is.null(reason)) {
      answer$type <- "unresolved"
      answer$intervals <- .ratio_intervals()
      answer$numerical_status <- reason
    }
    answer
  })
  structure(out, class = c("ratio_confidence_sets", "list"), level = level)
}

# Apply powers in representable chunks; exponent may exceed the range of 2^x.
.ratio_pow2 <- function(x, exponent) {
  exponent <- rep_len(exponent, length(x))
  repeat {
    step <- pmax(-1022, pmin(1023, exponent))
    x <- x * 2^step
    exponent <- exponent - step
    if (!any(exponent != 0)) {
      return(x)
    }
  }
}

.ratio_scale_exponent <- function(magnitude) {
  if (magnitude == 0) {
    return(0)
  }
  max(-1074, min(1023, floor(log2(magnitude))))
}

.ratio_validate_covariance <- function(V, target, tol) {
  arg <- sprintf("component_vcov[[%d]]", target)
  V <- .validate_matrix(V, arg, numeric = TRUE, finite = TRUE)
  if (!identical(dim(V), c(2L, 2L))) .abort_bad_arg(arg, "must be a paired 2-by-2 matrix")
  if (any(diag(V) < 0)) .abort_bad_arg(arg, "must be positive semidefinite")
  difference <- abs(V[1L, 2L] - V[2L, 1L])
  if (difference > tol * max(abs(V[1L, 2L]), abs(V[2L, 1L]))) {
    .abort_bad_arg(arg, "must be symmetric within the relative tolerance")
  }
  status <- if (difference == 0) "ok" else "symmetrized"
  if (difference != 0) V[1L, 2L] <- V[2L, 1L] <- V[1L, 2L] / 2 + V[2L, 1L] / 2
  if (any(diag(V) == 0)) {
    if (V[1L, 2L] != 0) .abort_bad_arg(arg, "must be positive semidefinite")
  } else {
    powers <- vapply(sqrt(diag(V)), .ratio_scale_exponent, numeric(1L))
    diagonal <- .ratio_pow2(diag(V), -2 * powers)
    off <- .ratio_pow2(V[1L, 2L], -sum(powers))
    rho <- off / sqrt(diagonal[1L]) / sqrt(diagonal[2L])
    if (!is.finite(rho)) .abort_bad_arg(arg, "must be positive semidefinite")
    # Normalization bounds the diagonal entries near one. Rounded square
    # roots/divisions can hide a negative determinant or put an exact singular
    # covariance just above rho=1. Certify the represented products inside
    # that rounding envelope, without discarding either product's residual.
    boundary <- abs(abs(rho) - 1) <= 8 * .Machine$double.eps
    if (!boundary && abs(rho) > 1 + tol) .abort_bad_arg(arg, "must be positive semidefinite")
    negative <- if (boundary) {
      !.exact_sum_products_nonnegative(c(diagonal[1L], off), c(diagonal[2L], -off))
    } else {
      abs(rho) > 1
    }
    if (negative) {
      if (tol == 0 || abs(rho) > 1 + tol) .abort_bad_arg(arg, "must be positive semidefinite")
      status <- "near_psd_boundary"
    }
  }
  list(covariance = V, status = status)
}

.ratio_intervals <- function(...) {
  values <- as.numeric(unlist(list(...), use.names = FALSE))
  matrix(values, ncol = 2L, byrow = TRUE, dimnames = list(NULL, c("lower", "upper")))
}

.ratio_quadratic_set <- function(a, b, c, tol = 64 * .Machine$double.eps) {
  abc <- c(a = a, b = b, c = c)
  exponent <- .ratio_scale_exponent(max(abs(abc)))
  scaled <- .ratio_pow2(abc, -exponent)
  a <- scaled[1L]
  b <- scaled[2L]
  c <- scaled[3L]
  discriminant <- b * b - 4 * a * c
  out <- list(
    type = "unresolved", intervals = .ratio_intervals(), quadratic = abc,
    scaled_quadratic = stats::setNames(scaled, c("a", "b", "c")),
    quadratic_scale_exponent = exponent, scaled_discriminant = unname(discriminant),
    discriminant = .ratio_pow2(discriminant, 2 * exponent), numerical_status = "ok"
  )
  finish <- function(type, intervals = .ratio_intervals(), status = "ok") {
    out$type <- type
    out$intervals <- intervals
    out$numerical_status <- status
    out
  }
  unresolved <- function(reason) finish("unresolved", status = reason)
  if (any(abc != 0 & scaled == 0)) {
    return(unresolved("scaling_underflow"))
  }
  if (a == 0) {
    if (b == 0) {
      if (c <= 0) {
        return(finish("all_real", .ratio_intervals(c(-Inf, Inf))))
      }
      return(finish("empty"))
    }
    root <- -c / b
    if (!is.finite(root) || (root == 0 && c != 0)) {
      return(unresolved("endpoint_range"))
    }
    bounds <- if (b > 0) c(-Inf, root) else c(root, Inf)
    return(finish("half_line", .ratio_intervals(bounds)))
  }
  # Exact factorization avoids squaring a tiny linear coefficient or treating
  # theta^2 times an arbitrary nonzero constant as an uncertain repeated root.
  if (c == 0) {
    if (b == 0) {
      if (a < 0) {
        return(finish("all_real", .ratio_intervals(c(-Inf, Inf))))
      }
      return(finish("singleton", .ratio_intervals(c(0, 0))))
    }
    root <- -b / a
    if (!is.finite(root) || root == 0) {
      return(unresolved("endpoint_range"))
    }
    roots <- sort(c(0, root))
    if (a > 0) {
      return(finish("bounded", .ratio_intervals(roots)))
    }
    return(finish("two_rays", .ratio_intervals(c(-Inf, roots[1L]), c(roots[2L], Inf))))
  }
  # Certify exact repeated roots on a dyadic grid whose integer products fit
  # exactly in a double. Other near-zero discriminants remain unresolved.
  integers <- scaled * 2^22
  certified <- all(integers == round(integers)) && all(abs(integers) <= 2^24)
  if (certified) {
    discriminant <- (integers[2L]^2 - 4 * integers[1L] * integers[3L]) / 2^44
    out$scaled_discriminant <- unname(discriminant)
    out$discriminant <- .ratio_pow2(discriminant, 2 * exponent)
  }
  size <- b * b + abs(4 * a * c)
  if (abs(discriminant) <= tol * size && !(discriminant == 0 && certified)) {
    return(unresolved("near_boundary"))
  }
  if (discriminant < 0 || (discriminant == 0 && a < 0)) {
    if (a < 0) {
      return(finish("all_real", .ratio_intervals(c(-Inf, Inf))))
    }
    return(finish("empty"))
  }
  if (discriminant == 0) {
    root <- -b / (2 * a)
    if (!is.finite(root)) {
      return(unresolved("endpoint_range"))
    }
    return(finish("singleton", .ratio_intervals(c(root, root))))
  }
  q <- -.5 * (b + if (b >= 0) sqrt(discriminant) else -sqrt(discriminant))
  roots <- sort(c(q / a, c / q))
  if (any(!is.finite(roots)) || (c != 0 && any(roots == 0)) || roots[1L] == roots[2L]) {
    return(unresolved("endpoint_range"))
  }
  if (a > 0) {
    return(finish("bounded", .ratio_intervals(roots)))
  }
  finish("two_rays", .ratio_intervals(c(-Inf, roots[1L]), c(roots[2L], Inf)))
}

#' @rdname ratio_confidence_set
#' @export
print.ratio_confidence_sets <- function(x, digits = max(3L, getOption("digits") - 2L), ...) {
  cat("Paired ratio confidence sets\n")
  for (j in seq_along(x)) {
    record <- x[[j]]
    shown <- vapply(seq_len(nrow(record$intervals)), function(k) {
      bounds <- record$intervals[k, ]
      paste0(
        if (is.infinite(bounds[1L])) "(" else "[",
        format(bounds[1L], digits = digits), ", ", format(bounds[2L], digits = digits),
        if (is.infinite(bounds[2L])) ")" else "]"
      )
    }, character(1L))
    cat(sprintf(
      "  target %d: %s %s (%s)\n", j, record$type,
      paste(shown, collapse = " U "), record$numerical_status
    ))
  }
  invisible(x)
}
