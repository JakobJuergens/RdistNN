#!/usr/bin/env Rscript
# Release validation only. Never changes the paper or promotes numerical output.
# See --help for the frozen-packet and installation-identity contracts.

file_sha256 <- function(path) digest::digest(file = path, algo = "sha256")

checked_hash <- function(path, expected) {
  if (!file.exists(path) || !identical(file_sha256(path), expected))
    stop("Hash mismatch: ", path, call. = FALSE)
  invisible(path)
}

new_output <- function(path) {
  if (file.exists(path)) stop("Output directory already exists: ", path, call. = FALSE)
  if (!dir.create(path, recursive = TRUE)) stop("Cannot create output directory.", call. = FALSE)
  normalizePath(path, winslash = "/", mustWork = TRUE)
}

compare_discrete <- function(a, b) list(pass = identical(a, b), bitwise = identical(a, b))

compare_numeric <- function(a, b, budget) {
  if (!identical(dim(a), dim(b)) || length(a) != length(b) ||
      !identical(is.na(a), is.na(b)) || !identical(is.finite(a), is.finite(b)))
    return(list(pass = FALSE, bitwise = FALSE, max_absolute = NULL, max_scaled = NULL))
  use <- is.finite(a) & is.finite(b)
  delta <- abs(a[use] - b[use])
  expand <- function(value) {
    block <- if (is.null(budget$block)) 1L else budget$block
    rep_len(rep(as.numeric(unlist(value)), each = block), length(a))[use]
  }
  allowance <- expand(budget$atol) + expand(budget$rtol) *
    pmax(abs(a[use]), abs(b[use]), expand(budget$scale))
  list(pass = all(delta <= allowance) && identical(a[!use], b[!use]),
    bitwise = identical(serialize(a, NULL, version = 3), serialize(b, NULL, version = 3)),
    max_absolute = max(c(0, delta)),
    max_scaled = max(c(0, delta / allowance)))
}

comparison_args <- function(args) {
  allowed <- c("paper", "historical-library", "candidate-library", "historical-identity",
    "candidate-identity", "cases", "output-dir", "child", "library", "identity", "help")
  out <- list()
  for (arg in args) {
    key <- sub("=.*$", "", sub("^--", "", arg))
    if (!startsWith(arg, "--") || !key %in% allowed)
      stop("Unknown argument: ", arg, call. = FALSE)
    if (key %in% names(out)) stop("Repeated argument: --", key, call. = FALSE)
    if (key %in% c("help", "child")) out[[key]] <- TRUE else {
      if (!grepl("=.+$", arg)) stop("Supply --", key, "=PATH", call. = FALSE)
      out[[key]] <- sub("^[^=]*=", "", arg)
    }
  }
  out
}

read_packets <- function(manifest_path, paper) {
  manifest <- jsonlite::read_json(manifest_path, simplifyVector = FALSE)
  if (!identical(manifest$schema_version, 1L) || !length(manifest$cases))
    stop("Expected a schema_version 1 manifest with cases.", call. = FALSE)
  for (name in names(manifest$paper_files))
    checked_hash(file.path(paper, name), manifest$paper_files[[name]])
  ids <- vapply(manifest$cases, `[[`, character(1), "id")
  if (anyDuplicated(ids)) stop("Case identifiers must be unique.", call. = FALSE)
  packets <- lapply(manifest$cases, function(case) {
    path <- file.path(dirname(manifest_path), case$file)
    checked_hash(path, case$sha256)
    x <- readRDS(path)
    n <- nrow(x$X)
    if (!identical(x$id, case$id) || n != nrow(x$PHI) ||
        any(!is.finite(x$PHI)) || any(!is.finite(x$X)) ||
        n != length(x$row_ids) || anyDuplicated(x$row_ids) ||
        x$s < 1 || x$s >= n || x$k < 1 || x$k > x$s ||
        !identical(dim(x$ranks), c(n, nrow(x$X0))) ||
        length(x$wn) != n || length(x$wnm1) != n - 1L ||
        any(!is.finite(c(x$wn, x$wnm1))) ||
        !x$geometry %in% c("l2", "mahalanobis") ||
        any(!is.finite(x$ranks)) || any(x$ranks != floor(x$ranks)) ||
        !all(apply(x$ranks, 2L, function(r) identical(sort(as.integer(r)), seq_len(n)))))
      stop("Invalid frozen packet: ", case$id, call. = FALSE)
    required <- c("full_weights", "deletion_weights", "full_mass", "deletion_mass", "ess",
      "components", "deleted", "deleted_averages", "centered", "pseudo", "component_covariance",
      "cate_target_covariance", "cate_estimate", "cate_se", "estimate", "se", "lower", "upper",
      "local_limit_variance", "distances")
    if (!all(required %in% names(case$budgets)))
      stop("Missing prespecified budgets: ", case$id, call. = FALSE)
    for (budget in case$budgets) for (field in c("atol", "rtol", "scale")) {
      value <- unlist(budget[[field]])
      if (!length(value) || !is.numeric(value) || any(!is.finite(value)) || any(value < 0))
        stop("Invalid prespecified budget: ", case$id, call. = FALSE)
    }
    x
  })
  list(manifest = manifest, packets = packets)
}

load_identified_package <- function(lib, identity_path) {
  if ("RdistNN" %in% loadedNamespaces())
    stop("RdistNN was already loaded; use a fresh Rscript --vanilla process.", call. = FALSE)
  identity <- jsonlite::read_json(identity_path, simplifyVector = TRUE)
  expected <- normalizePath(file.path(lib, "RdistNN"), winslash = "/", mustWork = TRUE)
  # Explicit lib.loc prevents the user's normal library from selecting RdistNN.
  ns <- loadNamespace("RdistNN", lib.loc = normalizePath(lib, mustWork = TRUE))
  actual <- normalizePath(getNamespaceInfo(ns, "path"), winslash = "/", mustWork = TRUE)
  if (!identical(actual, expected)) stop("Wrong RdistNN installation loaded.", call. = FALSE)
  version <- as.character(utils::packageVersion("RdistNN", lib.loc = lib))
  dll <- getLoadedDLLs()[["RdistNN"]][["path"]]
  dll <- normalizePath(dll, winslash = "/", mustWork = TRUE)
  if (!startsWith(dll, paste0(actual, "/"))) stop("Native library outside selected package.", call. = FALSE)
  if (version != identity$version) stop("Package version identity mismatch.", call. = FALSE)
  files <- names(identity$files)
  if (is.null(identity$source_identity) || !length(files))
    stop("Identity must bind source_identity and installed files including the DLL.", call. = FALSE)
  for (name in files) {
    if (grepl("^/|^[A-Za-z]:|\\\\", name) ||
        any(strsplit(name, "/", fixed = TRUE)[[1L]] %in% c("", ".", "..")) ||
        endsWith(name, "/"))
      stop("Identity files require canonical relative paths: ", name, call. = FALSE)
    path <- normalizePath(file.path(actual, name), winslash = "/", mustWork = TRUE)
    if (!startsWith(path, paste0(actual, "/")))
      stop("Identity file resolves outside the selected package.", call. = FALSE)
    checked_hash(path, identity$files[[name]])
  }
  dll_relative <- substring(dll, nchar(actual) + 2L)
  if (!dll_relative %in% files)
    stop("Loaded native library is not bound by the expected identity hashes.", call. = FALSE)
  list(namespace = ns, record = list(package_path = actual, version = version,
    source_identity = identity$source_identity, identity_sha256 = file_sha256(identity_path),
    installed_files = as.list(identity$files), native_library = dll,
    native_sha256 = file_sha256(dll), pid = Sys.getpid(), R = R.version.string,
    platform = R.version$platform, runtime = capture.output(sessionInfo()),
    external_libraries = as.list(extSoftVersion()),
    thread_environment = as.list(Sys.getenv(c("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"))),
    public_engine_sha256 = digest::digest(body(get("jackknife_rank_full", ns)), algo = "sha256")))
}

evaluate_tables <- function(x, ranks, wn, wd, ns) {
  n <- nrow(ranks)
  jk <- get("jackknife_rank_full", ns)(ranks, x$PHI, wn, wd)
  q <- ncol(x$PHI); m <- ncol(ranks)
  centered <- pseudo <- jk$theta_loo
  component_cov <- array(0, c(q, q, m))
  target_cov <- matrix(0, m, m)
  estimate <- se <- numeric(m)
  status <- character(m)
  for (j in seq_len(m)) {
    deleted <- matrix(jk$theta_loo[, , j], n, q)
    centered[, , j] <- sweep(deleted, 2L, colMeans(deleted))
    pseudo[, , j] <- n * matrix(jk$theta0[, j], n, q, byrow = TRUE) - (n - 1) * deleted
    component_cov[, , j] <- (n - 1) / n * crossprod(centered[, , j])
    estimate[j] <- -jk$theta0[3L, j] / jk$theta0[2L, j]
    score <- pseudo[, 2L, j] * estimate[j] + pseudo[, 3L, j]
    score <- score - mean(score)
    se[j] <- sqrt(sum(score^2) / (n * (n - 1))) / abs(jk$theta0[2L, j])
    status[j] <- if (!is.finite(estimate[j]) || !is.finite(se[j])) "unavailable" else if (se[j] == 0) "zero" else "ok"
  }
  cate_centered <- matrix(centered[, 1L, ], n, m)
  target_cov[,] <- (n - 1) / n * crossprod(cate_centered)
  list(row_ids = x$row_ids, target_ids = x$target_ids, ranks = ranks,
    support = which(wn != 0), deletion_support = which(wd != 0),
    full_weights = wn, deletion_weights = wd,
    full_mass = sum(wn), deletion_mass = sum(wd), ess = 1 / sum(wn^2),
    components = jk$theta0, deleted = jk$theta_loo,
    deleted_averages = apply(jk$theta_loo, c(2L, 3L), mean), centered = centered,
    pseudo = pseudo, component_covariance = component_cov,
    cate_target_covariance = target_cov, cate_estimate = as.numeric(jk$theta0[1L, ]),
    cate_se = sqrt(diag(target_cov)), estimate = estimate, se = se,
    lower = estimate - qnorm(.975) * se, upper = estimate + qnorm(.975) * se,
    local_limit_variance = (se^2) / sum(wn^2), status = status,
    includes_zero = (estimate - qnorm(.975) * se <= 0) & (estimate + qnorm(.975) * se >= 0),
    displayed_estimate = round(estimate, x$display_digits),
    displayed_lower = round(estimate - qnorm(.975) * se, x$display_digits),
    displayed_upper = round(estimate + qnorm(.975) * se, x$display_digits))
}

evaluate_case <- function(x, ns) {
  # A fixes the represented historical/paper tables. B changes only their construction.
  arithmetic <- evaluate_tables(x, x$ranks, x$wn, x$wnm1, ns)
  weight_spec <- get("rank_weight_dnn_exact", ns)(x$k, x$s)
  wn <- weight_spec$table(nrow(x$X))
  wd <- weight_spec$table(nrow(x$X) - 1L)
  weights <- evaluate_tables(x, x$ranks, wn, wd, ns)
  # The actual user wrapper is a separate consumer: SIPP used this route,
  # whereas the simulation called the full-rank deletion engine directly.
  fit <- get("cate_dnn", ns)(X = x$X, X0 = x$X0, pseudo_outcome = x$PHI[, 1L],
    k = x$k, s = x$s, se = TRUE, return = c("ranks", "W", "dx"),
    distance = get(if (x$geometry == "mahalanobis") "distance_mahalanobis" else "distance_l2", ns)())
  wrapper <- list(row_ids = x$row_ids, target_ids = x$target_ids,
    ranks = fit$ranks, support = fit$W != 0, full_weights = fit$W, distances = fit$dx,
    cate_estimate = as.numeric(fit$tau), cate_se = as.numeric(fit$se))
  wrapper$lower <- wrapper$cate_estimate - qnorm(.975) * wrapper$cate_se
  wrapper$upper <- wrapper$cate_estimate + qnorm(.975) * wrapper$cate_se
  wrapper$displayed_estimate <- round(wrapper$cate_estimate, x$display_digits)
  wrapper$displayed_lower <- round(wrapper$lower, x$display_digits)
  wrapper$displayed_upper <- round(wrapper$upper, x$display_digits)
  wrapper$includes_zero <- wrapper$lower <= 0 & wrapper$upper >= 0
  # C explicitly uses raw coordinates; it never changes the paper's Euclidean specification.
  metric <- get(if (x$geometry == "mahalanobis") "distance_mahalanobis" else "distance_l2", ns)()
  ctx <- if (is.function(metric$setup)) metric$setup(x$X) else NULL
  distances <- metric$eval_batch(x$X, x$X0, ctx)
  ranks <- apply(distances, 2L, function(d) { r <- integer(length(d)); r[order(d, seq_along(d))] <- seq_along(d); r })
  ranks <- matrix(as.integer(ranks), nrow(x$X), nrow(x$X0))
  geometry <- evaluate_tables(x, ranks, wn, wd, ns)
  geometry$distances <- distances
  geometry$numerical_rank <- if (!is.null(ctx$rank)) ctx$rank else if (!is.null(ctx$Si)) {
    covariance <- stats::cov(x$X)
    if (!inherits(try(solve(covariance), silent = TRUE), "try-error")) ncol(x$X) else {
      singular <- svd(covariance, nu = 0, nv = 0)$d
      as.integer(sum(singular > 1e-10 * max(singular, 1e-300)))
    }
  } else ncol(x$X)
  # Exercise backward-compatible inverse snapshots independently of new fitted geometry.
  if (!is.null(x$legacy_context))
    geometry$legacy_distances <- metric$eval_batch(x$X, x$X0, x$legacy_context)
  list(A_arithmetic = arithmetic, B_exact_weights = weights, C_geometry = geometry,
    D_wrapper = wrapper)
}

comparison_child <- function(opts) {
  loaded <- load_identified_package(opts$library, opts$identity)
  frozen <- read_packets(opts$cases, opts$paper)
  # Trace the native wrapper to establish which route actually ran, independently
  # of the option value. This instrumentation does not alter its calculation.
  assign(".paper_native_calls", 0L, envir = .GlobalEnv)
  trace("rustat_jackknife_rank_full", where = loaded$namespace,
    tracer = quote(.GlobalEnv$.paper_native_calls <- .GlobalEnv$.paper_native_calls + 1L),
    print = FALSE)
  on.exit(untrace("rustat_jackknife_rank_full", where = loaded$namespace), add = TRUE)
  modes <- list()
  engine_calls <- list()
  for (mode in c("compiled", "forced_R")) {
    assign(".paper_native_calls", 0L, envir = .GlobalEnv)
    options(RdistNN.use_compiled = identical(mode, "compiled"))
    if (mode == "compiled" && !is.loaded("_RdistNN_rustat_jackknife_rank_full", PACKAGE = "RdistNN"))
      stop("Required native deletion symbol is unavailable.", call. = FALSE)
    modes[[mode]] <- setNames(lapply(frozen$packets, evaluate_case, ns = loaded$namespace),
      vapply(frozen$packets, `[[`, character(1), "id"))
    engine_calls[[mode]] <- get(".paper_native_calls", .GlobalEnv)
    if ((mode == "compiled" && engine_calls[[mode]] == 0L) ||
        (mode == "forced_R" && engine_calls[[mode]] != 0L))
      stop("Observed engine route disagrees with requested mode: ", mode, call. = FALSE)
  }
  # Recheck all installation hashes after calculation as well as before it.
  for (name in names(loaded$record$installed_files))
    checked_hash(file.path(loaded$record$package_path, name), loaded$record$installed_files[[name]])
  loaded$record$engine_calls <- engine_calls
  saveRDS(list(identity = loaded$record, modes = modes), file.path(opts[["output-dir"]], "result.rds"), version = 3)
}

summarize_comparison <- function(a, b, manifest) {
  discrete <- c("row_ids", "target_ids", "ranks", "support", "deletion_support", "status", "includes_zero", "numerical_rank")
  discrete <- c(discrete, "wrapper_includes_zero")
  display <- c("displayed_estimate", "displayed_lower", "displayed_upper",
    "wrapper_displayed_estimate", "wrapper_displayed_lower", "wrapper_displayed_upper")
  results <- list()
  for (mode in names(a$modes)) for (case in manifest$cases) for (part in names(a$modes[[mode]][[case$id]])) {
    old <- a$modes[[mode]][[case$id]][[part]]; new <- b$modes[[mode]][[case$id]][[part]]
    fields <- lapply(names(old), function(field) {
      if (field %in% c(discrete, display)) r <- compare_discrete(old[[field]], new[[field]]) else {
        budget <- case$budgets[[field]]
        if (is.null(budget)) stop("Missing prespecified budget for ", case$id, "/", field, call. = FALSE)
        r <- compare_numeric(old[[field]], new[[field]], budget)
        r$budget <- budget
      }
      r$quantity <- field
      r$classification <- if (r$bitwise) "BITWISE_IDENTICAL" else if (field %in% display)
        if (r$pass) "DISPLAY_IDENTICAL" else "DISPLAY_CHANGE_REVIEW" else if (field %in% discrete)
        if (r$pass) "EXACT_DISCRETE_VALUES" else "DISCRETE_CHANGE_REVIEW" else if (r$pass)
        "WITHIN_BUDGET" else "UNEXPLAINED_CHANGE_REVIEW"
      r
    })
    results[[length(results) + 1L]] <- list(case = case$id, mode = mode, comparison = part,
      scope = case$scope, input_sha256 = case$sha256, fields = fields,
      pass = all(vapply(fields, `[[`, logical(1), "pass")))
  }
  results
}

comparison_main <- function(args = commandArgs(TRUE)) {
  opts <- comparison_args(args)
  if (isTRUE(opts$help)) {
    cat("Usage: Rscript --vanilla compare_paper_dependency.R --paper=READ_ONLY_CHECKOUT\n",
      " --historical-library=LIB --historical-identity=IDENTITY.json\n",
      " --candidate-library=LIB --candidate-identity=IDENTITY.json\n",
      " --cases=FROZEN_MANIFEST.json --output-dir=NEW_DIRECTORY\n\n",
      "Identity: version, source_identity, files {relative_installed_path:sha256}, including libs/.\n",
      "Cases: schema_version 1, paper_files hashes, cases [{id,file,sha256,scope,budgets}].\n",
      "Each private RDS packet: id,X,X0,PHI (phi,B,b),ranks,wn,wnm1,k,s,row_ids,target_ids,\n",
      "geometry (l2|mahalanobis),display_digits; optional legacy_context. Budgets map each\n",
      "continuous output to {atol,rtol,scale} declared before execution. Absolute units\n",
      "must match the quantity; component matrices can use conformable scale matrices.\n",
      "Runs both compiled and forced-R modes in separate fresh installation processes.\n",
      "Outputs JSON/Markdown summaries and PRIVATE observation-level RDS results.\n",
      "Exit 2 means review is required; no historical evidence is replaced.\n", sep = "")
    return(invisible(NULL))
  }
  required <- if (isTRUE(opts$child)) c("library", "identity", "paper", "cases", "output-dir") else
    c("paper", "historical-library", "candidate-library", "historical-identity", "candidate-identity", "cases", "output-dir")
  if (any(!required %in% names(opts))) stop("Missing arguments; use --help.", call. = FALSE)
  if (isTRUE(opts$child)) return(comparison_child(opts))
  # Inspect hashes and budgets before starting either calculation.
  frozen <- read_packets(opts$cases, opts$paper)
  output <- new_output(opts[["output-dir"]])
  script <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)))
  for (side in c("historical", "candidate")) {
    path <- file.path(output, side); dir.create(path)
    child_args <- c("--vanilla", script, "--child", paste0("--library=", opts[[paste0(side, "-library")]]),
      paste0("--identity=", opts[[paste0(side, "-identity")]]), paste0("--paper=", opts$paper),
      paste0("--cases=", normalizePath(opts$cases)), paste0("--output-dir=", path))
    code <- system2(file.path(R.home("bin"), "Rscript"), shQuote(child_args),
      stdout = file.path(path, "process.log"), stderr = file.path(path, "process.log"),
      env = c("SIM_RDISTNN_SOURCE_DIR=", "OMP_NUM_THREADS=1", "OPENBLAS_NUM_THREADS=1", "MKL_NUM_THREADS=1"))
    if (code != 0L) stop(side, " process failed; see ", path, "/process.log", call. = FALSE)
  }
  a <- readRDS(file.path(output, "historical/result.rds")); b <- readRDS(file.path(output, "candidate/result.rds"))
  results <- summarize_comparison(a, b, frozen$manifest)
  pass <- all(vapply(results, `[[`, logical(1), "pass"))
  summary <- list(schema_version = 1L, created_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
    classification = if (pass) "WITHIN_DECLARED_BUDGET" else "REVIEW_REQUIRED",
    manifest_sha256 = file_sha256(opts$cases), manifest = frozen$manifest,
    historical = a$identity, candidate = b$identity, comparisons = results,
    reproduction_reference = "Historical dependency remains authoritative; candidate output is validation only.")
  jsonlite::write_json(summary, file.path(output, "summary.json"), pretty = TRUE, auto_unbox = TRUE, digits = NA)
  lines <- c("# Conditional Z-estimation dependency comparison", "", summary$classification, "",
    paste("Historical package:", a$identity$package_path), paste("Candidate package:", b$identity$package_path), "",
    "A fixes ranks and represented tables; B constructs exact tables at fixed ranks; C fits raw-coordinate geometry.",
    "D separately evaluates the end-to-end CATE wrapper and records its ranks, support, weights and distances.",
    "Mahalanobis policy comparisons do not describe the paper's fixed Euclidean estimator.",
    "All covariance fields are on estimator scale except the explicitly named local_limit_variance.",
    "Observation-level results remain private. The historical dependency remains the reproduction reference.", "",
    "| Case | Mode | Comparison | Result | Nonidentical quantities (maximum absolute error) |",
    "|---|---|---|---|---|")
  for (r in results) {
    changed <- Filter(function(f) !f$bitwise, r$fields)
    detail <- vapply(changed, function(f) paste0(f$quantity, ": ", f$classification,
      if (!is.null(f$max_absolute)) paste0(" (", format(f$max_absolute, digits = 4), ")") else ""), character(1))
    lines <- c(lines, paste0("| ", r$case, " | ", r$mode, " | ", r$comparison, " | ",
      if (r$pass) "PASS" else "REVIEW", " | ", paste(detail, collapse = "; "), " |"))
  }
  if (length(frozen$manifest$untested)) lines <- c(lines, "", "Not checked:", unlist(frozen$manifest$untested))
  writeLines(lines, file.path(output, "report.md"))
  cat(summary$classification, ": ", output, "\n", sep = "")
  if (!pass) quit(status = 2L)
  invisible(summary)
}

if (sys.nframe() == 0L) comparison_main()
