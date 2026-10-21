#!/usr/bin/env Rscript
# Manuscript displays from accepted portable evidence only. No fits or RDS reads.
# From the repository root:
#   Rscript Codes/Simulations/affine_scores/campaign_tools/render_manuscript_simulations.R
#   Rscript Codes/Simulations/affine_scores/campaign_tools/render_manuscript_simulations.R --check
#   Rscript Codes/Simulations/affine_scores/campaign_tools/render_manuscript_simulations.R --tables-only
# --tables-only preserves verified figures when only tabular formatting changes.
# Captions and table environments belong to the manuscript, not these fragments.
args <- commandArgs(TRUE)
if (length(args) > 1L || (length(args) && !args %in% c("--check", "--tables-only")))
  stop("Usage: render_manuscript_simulations.R [--check|--tables-only]")
check_only <- identical(args, "--check")
tables_only <- identical(args, "--tables-only")
script <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)))
repo <- normalizePath(file.path(dirname(script), "../../../.."))
relative <- function(x) sub(paste0("^", repo, "/"), "", normalizePath(x, mustWork = FALSE))
accepted <- file.path(repo, "markdown_notes/affine_scores/simulation_program/full_10000_20261012")
analysis <- file.path(repo, "markdown_notes/affine_scores/simulation_program/production_analysis/data")
out <- file.path(repo, "Artifacts/affine_scores/manuscript_simulations_20261012")
tables <- file.path(repo, "Tables/Simulations/affine_production")
previous <- if (tables_only) jsonlite::read_json(file.path(out, "manifest.json"), simplifyVector = TRUE) else NULL
for (pkg in c("ggplot2", "patchwork", "jsonlite", "digest"))
  if (!requireNamespace(pkg, quietly = TRUE)) stop("Required package: ", pkg)
sha <- function(f) digest::digest(file = f, algo = "sha256")
inputs <- list()
validate_hash <- function(f, expected = NULL) {
  actual <- sha(f)
  if (!is.null(expected) && !identical(actual, expected)) stop("Input hash mismatch: ", relative(f))
  inputs[[relative(f)]] <<- actual
  invisible(actual)
}

# Validate the accepted packet before reading numerical tables.
checksum_file <- file.path(accepted, "SHA256SUMS")
validate_hash(checksum_file)
checksums <- readLines(checksum_file)
for (line in checksums) {
  expected <- substr(line, 1L, 64L)
  f <- sub("^[0-9a-f]{64}  ", "", line)
  if (!grepl("^[0-9a-f]{64}$", expected) || f == line) stop("Malformed SHA256SUMS entry")
  validate_hash(file.path(accepted, f), expected)
}
am_path <- file.path(analysis, "analysis_manifest.json")
validate_hash(am_path)
am <- jsonlite::read_json(am_path, simplifyVector = TRUE)
audit <- jsonlite::read_json(file.path(accepted, "audit.json"), simplifyVector = TRUE)
stopifnot(audit$status == "completed_and_audited_production", audit$total_rows == 2516400,
          audit$common_independent_draw_families == 10000,
          am$status == "analysis_of_audited_saved_production", !am$estimation_performed,
          am$source_audit_sha256 == sha(file.path(accepted, "audit.json")))
for (i in seq_len(nrow(am$data_hashes))) {
  h <- am$data_hashes[i, ]
  validate_hash(file.path(analysis, h$file), h$sha256)
}
proposal <- jsonlite::read_json(file.path(accepted, "provenance/proposal.json"), simplifyVector = TRUE)
for (f in names(proposal$source_hashes))
  validate_hash(file.path(dirname(dirname(script)), f), proposal$source_hashes[[f]])
read <- function(root, name) read.csv(file.path(root, paste0(name, ".csv")), stringsAsFactors = FALSE)
R <- read(analysis, "raw_metrics"); P <- read(accepted, "performance")
D <- read(analysis, "diagnostic_moments"); Q <- read(analysis, "qq")
C <- read(analysis, "paired_contrasts"); J <- read(analysis, "projection_references")
G <- read(analysis, "gram_conditional")
A <- read(accepted, "paired"); S <- read(accepted, "spatial_decomposition")
V <- read(accepted, "covariance_ratios"); W <- read(accepted, "variance_pairs")
T <- read(accepted, "training_pairs"); L <- read(accepted, "localizer_pairs")
E <- read(accepted, "gram_events")
key <- function(x, fields) do.call(paste, c(x[fields], sep = "|"))
row_keys <- c("profile", "cell_id", "application", "training_method", "localization_method", "variance_method")
cell_keys <- c("profile", "cell_id", "application", "training_method", "localization_method")
stopifnot(nrow(R) == 918L, identical(key(R, row_keys), key(P, row_keys)),
          !anyDuplicated(key(R, row_keys)))
overlap <- unique(R$cell_id[R$profile == "dnn"])
canonical <- function(x) {
  rule <- !(x$profile == "kernel_comparison" & x$localization_method == "dnn" & x$cell_id %in% overlap)
  if ("canonical" %in% names(x)) stopifnot(identical(rule, x$canonical))
  x[rule, , drop = FALSE]
}
P$canonical <- R$canonical
P$s <- R$s; P$lambda <- R$lambda
P$coverage_empirical_center <- R$coverage_empirical_center
P$coverage_realized_weighted_truth <- R$coverage_realized_weighted_truth
P$point_t_excess_kurtosis <- R$point_t_excess_kurtosis
P$mc_z_excess_kurtosis <- R$mc_z_excess_kurtosis
R <- canonical(R); P <- canonical(P); D <- canonical(D); Q <- canonical(Q)
A <- canonical(A); S <- canonical(S); V <- canonical(V); W <- canonical(W)
T <- canonical(T)
stopifnot(nrow(P[P$design_id == "C0" & P$n == 8192 & P$localization_method == "dnn" &
                P$training_method == "full" & P$variance_method == "jackknife", ]) == 2L)
main_methods <- c("full", "pooled_five_fold", "oracle")
stopifnot(all(P$completed[P$training_method %in% main_methods] == 10000),
          all(P$calibration_pairs[P$training_method %in% main_methods] == 10000),
          all(P$failures[P$training_method %in% main_methods] == 0),
          all(C$paired == 10000), all(C$excluded == 0), all(C$excluded_points == 0))
for (f in c("decomposition_scaled_mean", "decomposition_scaled_centered", "grid_prediction_radius_squared")) {
  m3 <- D[D$profile == "M3" & D$field == f, ]
  stopifnot(nrow(m3) > 0, all(m3$observed == 0), all(is.na(m3$mean)))
}

one <- function(x, ...) {
  filters <- list(...)
  keep <- rep(TRUE, nrow(x))
  for (f in names(filters)) keep <- keep & !is.na(x[[f]]) & x[[f]] %in% filters[[f]]
  z <- x[keep, , drop = FALSE]
  if (nrow(z) != 1L) stop("Expected one row; found ", nrow(z), " for ", paste(names(filters), filters, collapse = ", "))
  z
}
perf <- function(design, n, app, order = 1, training = "full", localizer = "dnn", route = "jackknife")
  one(P, design_id = design, n = n, application = app, order_mult = order,
      training_method = training, localization_method = localizer, variance_method = route)
diag <- function(design, n, app, field, order = 1, training = "full", localizer = "dnn")
  one(D, design_id = design, n = n, application = app, order_mult = order,
      training_method = training, localization_method = localizer, field = field)
omega <- function(app) ifelse(app == "cate", 1 / (4 * sqrt(15)), 2 / 31)
app_name <- function(app) ifelse(app == "cate", "CATE", "Slope")
decorate <- function(z) {
  z$app <- factor(app_name(z$application), c("CATE", "Slope"))
  if ("training_method" %in% names(z))
    z$method <- factor(ifelse(z$training_method == "full", "Full fit", "True nuisance"), c("Full fit", "True nuisance"))
  z
}
selections <- list()
record <- function(name, x, source, fields, extra = list()) {
  fields <- intersect(fields, names(x))
  selections[[name]] <<- c(list(source = source, rows = nrow(x), selected = unique(x[fields])), extra)
  invisible(x)
}
case_fields <- c(row_keys, "design_id", "n", "M", "K", "s", "order_mult")

library(ggplot2)
library(patchwork)
theme_set(theme_bw(base_size = 9.5, base_family = "sans") +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_line(colour = "grey92", linewidth = .2),
        strip.background = element_rect(fill = "grey95", colour = "grey70", linewidth = .25),
        strip.text = element_text(size = 9.5), axis.text = element_text(size = 9),
        axis.title = element_text(size = 9.5), legend.text = element_text(size = 9),
        legend.position = "bottom", legend.title = element_blank(),
        legend.key.width = grid::unit(1.3, "lines"), legend.key.height = grid::unit(.8, "lines"),
        panel.spacing = grid::unit(.6, "lines"), plot.margin = margin(5, 5, 5, 5)))
method_colours <- c("Full fit" = "#1F4B63", "True nuisance" = "#8A4F13")
method_styles <- list(scale_colour_manual(values = method_colours),
                      scale_linetype_manual(values = c("Full fit" = "solid", "True nuisance" = "dashed")),
                      scale_shape_manual(values = c("Full fit" = 16, "True nuisance" = 1)))
figures <- list()
save_plot <- function(name, p, height, receipt) {
  figures[[name]] <<- c(list(pdf = paste0(name, ".pdf"), png = paste0(name, ".png"),
                             width_inches = 5.4, height_inches = height, base_font_points = 9.5), receipt)
  if (tables_only) {
    stopifnot(isTRUE(all.equal(figures[[name]], previous$figures[[name]], tolerance = 0)))
    for (extension in c(".pdf", ".png")) {
      f <- relative(file.path(out, paste0(name, extension)))
      expected <- previous$output_hashes$sha256[previous$output_hashes$file == f]
      stopifnot(length(expected) == 1L, sha(file.path(repo, f)) == expected)
    }
    for (i in seq_len(nrow(previous$input_hashes))) {
      h <- previous$input_hashes[i, ]
      stopifnot(!is.null(inputs[[h$file]]), inputs[[h$file]] == h$sha256)
    }
  } else if (!check_only) {
    ggsave(file.path(out, paste0(name, ".pdf")), p, device = cairo_pdf, width = 5.4, height = height)
    ggsave(file.path(out, paste0(name, ".png")), p, device = "png", width = 5.4, height = height, dpi = 220)
  }
}
if (!check_only) {
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  dir.create(tables, recursive = TRUE, showWarnings = FALSE)
}

# Main performance: full and true nuisances, primary DNN jackknife, default order.
performance <- P[P$design_id %in% c("C0", "C1", "C2") & P$order_mult == 1 &
                   P$localization_method == "dnn" & P$variance_method == "jackknife" &
                   P$training_method %in% c("full", "oracle"), ]
performance <- decorate(performance[order(performance$application, performance$design_id, performance$mean_N), ])
stopifnot(nrow(performance) == 52L, all(performance$N_reference == performance$mean_N),
          all(performance$coverage >= .855 & performance$coverage <= .968))
record("performance", performance, "accepted/performance.csv + analysis/raw_metrics.csv", case_fields,
       list(metric = "point-target coverage", replication_count = 10000, x = "exact deterministic DNN N"))
nominal <- unique(performance[c("app", "design_id")]); nominal$x <- 9; nominal$y <- .954
transition <- performance[performance$application == "cate" & performance$design_id == "C0" &
                            performance$training_method == "full" & performance$n %in% c(512, 4096, 531441), ]
transition$label <- c("512\nM1", "4,096\nM2", "531,441\nM3")[match(transition$n, c(512, 4096, 531441))]
transition$label_y <- c(.878, .878, .878)[match(transition$n, c(512, 4096, 531441))]
transition$label_x <- ifelse(transition$n == 531441, 40, ifelse(transition$n == 512, 9.5, transition$mean_N))
p <- ggplot(performance, aes(mean_N, coverage, colour = method, linetype = method, shape = method)) +
  geom_hline(yintercept = .95, colour = "grey35", linetype = "dotted", linewidth = .45) +
  geom_line(linewidth = .55) + geom_point(size = 1.6, stroke = .55) + method_styles +
  geom_text(data = nominal, aes(x, y, label = "0.95"), inherit.aes = FALSE, size = 3.2, hjust = 0, colour = "grey30") +
  geom_text(data = transition, aes(label_x, label_y, label = label), inherit.aes = FALSE,
            size = 3.2, lineheight = .95, colour = "grey30") +
  facet_grid(app ~ design_id) +
  scale_x_log10(breaks = c(8.5, 15, 30, 53), labels = c("8.5", "15", "30", "53"),
                limits = c(7.5, 60), expand = c(0, 0)) +
  scale_y_continuous(breaks = c(.86, .90, .95), limits = c(.855, .968), expand = c(0, 0)) +
  labs(x = "Exact DNN local information N", y = "Point-target coverage") +
  guides(colour = guide_legend(nrow = 1), linetype = guide_legend(nrow = 1), shape = guide_legend(nrow = 1))
save_plot("performance", p, 4.15,
          list(panels = 6, coverage_reference = .95, coverage_limits = c(.855, .968),
               annotations = "C0 CATE: n=512/M1, n=4096/M2, n=531441/M3; dictionary fixed within each M"))

# Main mechanisms: paired q95 gaps and exact realized-truth decomposition.
gaps <- A[A$design_id == "C0" & A$order_mult == 1 & A$localization_method == "dnn", ]
gaps$app <- factor(app_name(gaps$application), c("CATE", "Slope"))
gaps$N <- vapply(seq_len(nrow(gaps)), function(i) perf("C0", gaps$n[i], gaps$application[i])$mean_N, numeric(1))
gaps$gap <- factor(ifelse(gaps$comparison == "full_vs_CF", "Full minus five-fold", "Full minus true nuisance"),
                   c("Full minus five-fold", "Full minus true nuisance"))
gaps <- gaps[order(gaps$app, gaps$gap, gaps$N), ]
stopifnot(nrow(gaps) == 28L, all(gaps$paired == 10000), all(gaps$excluded == 0),
          all(gaps$q95_absolute_Omega_units >= 0 & gaps$q95_absolute_Omega_units <= 1.8),
          max(abs(gaps$q95_absolute_Omega_units - gaps$q95_absolute / sqrt(omega(gaps$application)))) < 1e-12)
record("mechanisms_gaps", gaps, "accepted/paired.csv + accepted/performance.csv",
       c("profile", "cell_id", "application", "localization_method", "comparison", "n", "M", "order_mult"),
       list(metric = "type-8 q95(abs(sqrt(N)*(full-comparator)))/sqrt(Omega_point)", paired_count = 10000))
gap_plot <- ggplot(gaps, aes(N, q95_absolute_Omega_units, colour = gap, linetype = gap, shape = gap)) +
  geom_line(linewidth = .6) + geom_point(size = 1.65) + facet_grid(app ~ .) +
  scale_colour_manual(values = c("Full minus five-fold" = "#1F4B63", "Full minus true nuisance" = "#8A4F13")) +
  scale_linetype_manual(values = c("Full minus five-fold" = "solid", "Full minus true nuisance" = "dashed")) +
  scale_shape_manual(values = c("Full minus five-fold" = 16, "Full minus true nuisance" = 1)) +
  scale_x_log10(breaks = c(8.5, 15, 30, 53), labels = c("8.5", "15", "30", "53")) +
  scale_y_continuous(limits = c(0, 1.8), breaks = c(0, .5, 1, 1.5)) +
  labs(x = "Exact DNN N (C0)", y = "Paired absolute-gap q95\n(asymptotic SD units)") +
  guides(colour = guide_legend(nrow = 2, byrow = TRUE), linetype = guide_legend(nrow = 2, byrow = TRUE),
         shape = guide_legend(nrow = 2, byrow = TRUE))
centers <- S[S$localization_method == "dnn" & S$training_method == "full" &
               ((S$design_id == "C1" & S$n == 2048 & S$order_mult == 1) |
                (S$design_id == "C2" & S$n == 512 & S$order_mult == .4)), ]
centers <- decorate(centers)
centers$case <- factor(ifelse(centers$design_id == "C1", "C1\nn=2,048", "C2\nn=512\norder multiplier 0.4"),
                       c("C1\nn=2,048", "C2\nn=512\norder multiplier 0.4"))
centers$N <- vapply(seq_len(nrow(centers)), function(i)
  perf(centers$design_id[i], centers$n[i], centers$application[i], centers$order_mult[i])$mean_N, numeric(1))
centers$normalizer <- sqrt(centers$N) / sqrt(omega(centers$application))
centers$total_bias <- (centers$mean_drift + centers$mean_remainder) * centers$normalizer
stopifnot(nrow(centers) == 4L, all(centers$paired == 10000), all(centers$excluded == 0))
for (i in seq_len(nrow(centers))) {
  r <- centers[i, ]; stopifnot(abs((r$mean_drift + r$mean_remainder) -
    perf(r$design_id, r$n, r$application, r$order_mult)$bias) < 1e-12)
}
components <- do.call(rbind, lapply(c("mean_drift", "mean_remainder"), function(field) {
  z <- centers; z$value <- z[[field]] * z$normalizer
  z$component <- if (field == "mean_drift") "Spatial drift" else "Remainder"
  z
}))
components$component <- factor(components$component, c("Spatial drift", "Remainder"))
record("mechanisms_centering", centers, "accepted/spatial_decomposition.csv + accepted/performance.csv",
       c(cell_keys, "design_id", "n", "M", "order_mult"),
       list(normalization = "sqrt(exact N)/sqrt(Omega_point)", paired_count = 10000,
            remainder = "estimate minus application-specific realized weighted truth; includes oracle sampling noise"))
center_plot <- ggplot(components, aes(case, value, fill = component)) +
  geom_hline(yintercept = 0, colour = "grey55", linewidth = .3) +
  geom_col(position = position_dodge(width = .65), width = .28, colour = "grey20", linewidth = .35) +
  geom_point(data = centers, aes(case, total_bias, shape = "Total bias"), inherit.aes = FALSE, size = 2.5) +
  facet_grid(app ~ .) + scale_fill_manual(values = c("Spatial drift" = "grey65", "Remainder" = "white")) +
  scale_shape_manual(values = c("Total bias" = 18)) +
  scale_y_continuous(limits = c(-.04, 1.8), breaks = c(0, .5, 1, 1.5)) +
  labs(x = "Design and sample size", y = "Mean component\n(asymptotic SD units)") +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE), shape = guide_legend(nrow = 1)) +
  theme(legend.box = "vertical", legend.spacing.y = grid::unit(0, "pt"))
save_plot("mechanisms", gap_plot + center_plot + plot_layout(widths = c(1, 1.05)), 4.65,
          list(panels = 4, gap_reference = "Omega_CATE=1/(4 sqrt(15)); Omega_slope=2/31",
               normalization = "sqrt(N)/sqrt(Omega_point)", paired_count = 10000))

# Appendix QQ: central view retains .025/.975, with full tails in the old atlas.
qq_selected <- Q[Q$variance_method == "jackknife" & Q$localization_method == "dnn" &
                   Q$training_method %in% c("full", "oracle") & Q$definition == "point_t" &
                   Q$probability >= .025 & Q$probability <= .975 &
                   ((Q$design_id == "C0" & Q$n %in% c(512, 531441) & Q$order_mult == 1) |
                    (Q$design_id == "C2" & Q$n == 512 & Q$order_mult == .4)), ]
qq_selected <- decorate(qq_selected)
qq_selected$N <- vapply(seq_len(nrow(qq_selected)), function(i)
  perf(qq_selected$design_id[i], qq_selected$n[i], qq_selected$application[i], qq_selected$order_mult[i])$mean_N, numeric(1))
qq_selected$case <- paste0(qq_selected$design_id, ", n=",
  format(qq_selected$n, scientific = FALSE, big.mark = ",", trim = TRUE),
  ifelse(qq_selected$design_id == "C2", "\norder multiplier 0.4", ""),
  "\nN=", sprintf("%.1f", qq_selected$N), "; M=", qq_selected$M)
levels_qq <- unique(qq_selected$case[order(qq_selected$N)])
qq_selected$case <- factor(qq_selected$case, levels_qq)
qq_selected <- qq_selected[order(qq_selected$app, qq_selected$case, qq_selected$method, qq_selected$probability), ]
stopifnot(length(levels_qq) == 3L, nrow(qq_selected) == 1116L, all(qq_selected$observed == 10000),
          all(qq_selected$empirical_quantile >= -3 & qq_selected$empirical_quantile <= 4.7))
record("studentization", qq_selected, "analysis/qq.csv + accepted/performance.csv",
       c(row_keys, "definition", "design_id", "n", "M", "order_mult"),
       list(curves = 12, probabilities = sort(unique(qq_selected$probability)), observed_per_curve = 10000,
            order = "increasing exact DNN N", central_view = c(.025, .975), quantile_type = 8,
            full_tails = "Artifacts/affine_scores/production_simulations_20261012/12_CATE_point_student_QQ.pdf and 14_Slope_point_student_QQ.pdf"))
p <- ggplot(qq_selected, aes(normal_quantile, empirical_quantile, colour = method, linetype = method)) +
  geom_abline(slope = 1, intercept = 0, colour = "grey50", linetype = "dotted", linewidth = .45) +
  geom_line(linewidth = .65) + method_styles[1:2] + facet_grid(app ~ case) +
  scale_x_continuous(breaks = c(-2, 0, 2), limits = c(-2.05, 2.05)) +
  scale_y_continuous(breaks = c(-2, 0, 2, 4), limits = c(-3, 4.7)) +
  labs(x = "Standard normal quantile: central probabilities [0.025, 0.975]",
       y = "(Estimate - point target) / fitted SE")
save_plot("studentization", p, 4.5,
          list(panels = 6, definition = "point-centered jackknife-studentized statistic",
               central_probabilities = c(.025, .975), observed_per_curve = 10000, quantile_type = 8))

# Tables are tabular fragments, with no floats, captions, labels, or local macros.
fragments <- list(); table_receipts <- list()
fmt <- function(x, digits = 4) ifelse(is.na(x), "---", formatC(x, format = "f", digits = digits))
signed <- function(x, digits = 4) ifelse(is.na(x), "---", sprintf(paste0("%+.", digits, "f"), x))
# prettyNum drops the backslash in big.mark="\\,"; insert TeX spacing after formatting.
integer_tex <- function(x) vapply(format(x, scientific = FALSE, big.mark = ",", trim = TRUE),
  function(z) paste(strsplit(z, ",", fixed = TRUE)[[1]], collapse = "\\,"), character(1))
stopifnot(integer_tex(531441) == "531\\,441")
case_tex <- function(design, n, order = NULL)
  paste0(design, " $(", integer_tex(n), if (!is.null(order)) paste0(",", format(order, trim = TRUE)) else "", ")$")
row <- function(...) paste0(paste(c(...), collapse = " & "), " \\\\")
tabular <- function(name, columns, header, rows, receipt) {
  text <- c("% Generated by render_manuscript_simulations.R; edit the renderer, not this fragment.",
            paste0("\\begin{tabular}{", columns, "}"), "\\hline", header, "\\hline", rows,
            "\\hline", "\\end{tabular}")
  fragments[[paste0(name, ".tex")]] <<- text
  table_receipts[[name]] <<- c(list(rows = length(rows), columns = columns), receipt)
  if (!check_only) writeLines(text, file.path(tables, paste0(name, ".tex")), useBytes = TRUE)
}
apps <- c("cate", "slope")
main_cases <- data.frame(design = c("C0", "C2", "C0H", "C2S"), n = c(2048, 512, 8192, 2048))
main_rows <- main_records <- list()
for (i in seq_len(nrow(main_cases))) for (app in apps) {
  design <- main_cases$design[i]; n <- main_cases$n[i]
  pp <- rbind(perf(design, n, app), perf(design, n, app, localizer = "kernel_gaussian"),
              perf(design, n, app, localizer = "kernel_gaussian", route = "direct_squared"))
  main_records[[length(main_records) + 1L]] <- pp
  pairs <- paste0(fmt(pp$coverage), " (", fmt(pp$mean_length, 3), ")")
  main_rows[[length(main_rows) + 1L]] <- row(case_tex(design, n), app_name(app), pairs[1], pairs[2], pairs[3])
}
main_records <- do.call(rbind, main_records)
record("main_comparison", main_records, "accepted/performance.csv", case_fields, list(replications = 10000))
tabular("main_comparison", "@{}llrrr@{}",
        row("Case $(n)$", "App.", "\\shortstack{DNN\\\\jackknife}", "\\shortstack{Kernel\\\\jackknife}",
            "\\shortstack{Kernel\\\\conventional}"), unlist(main_rows),
        list(entries = "point coverage (mean interval length)", order_mult = 1, training = "full", replications = 10000))

schedule <- do.call(rbind, lapply(c(512, 2048, 4096, 8192, 16384, 65536, 531441), function(n) perf("C0", n, "cate")))
record("design_schedule", schedule, "accepted/performance.csv + analysis/raw_metrics.csv", case_fields)
schedule_rows <- vapply(seq_len(nrow(schedule)), function(i) {
  z <- schedule[i, ]; row(integer_tex(z$n), z$M, z$K, integer_tex(z$s), fmt(z$N_reference, 2))
}, character(1))
tabular("design_schedule", "@{}rrrrr@{}", row("$n$", "$M$", "$K$", "$s$", "$N$"), schedule_rows,
        list(design = "C0", order_mult = 1, N = "exact deterministic rank-weight effective information"))

appendix_cases <- data.frame(design = c("C0", "C0", "C1", "C2", "C2", "C0"),
                             n = c(512, 65536, 2048, 512, 512, 531441), order = c(1, 1, 1, .4, 1, 1))
arows <- crows <- list(); appendix_records <- list()
for (i in seq_len(nrow(appendix_cases))) for (app in apps) {
  design <- appendix_cases$design[i]; n <- appendix_cases$n[i]; order <- appendix_cases$order[i]
  pp <- do.call(rbind, lapply(main_methods, function(training) perf(design, n, app, order, training)))
  appendix_records[[length(appendix_records) + 1L]] <- pp
  z <- pp[1, ]
  if (i %in% c(1, 3, 4))
    arows[[length(arows) + 1L]] <- row(case_tex(design, n, order), app_name(app), fmt(z$bias), fmt(z$rmse),
                                      fmt(pp$coverage[1]), fmt(pp$coverage[2]), fmt(pp$coverage[3]))
  if (i %in% c(1, 4, 6))
    crows[[length(crows) + 1L]] <- row(case_tex(design, n, order), app_name(app), fmt(z$mean_SE_over_SD, 3),
                                      fmt(z$rms_SE_over_SD, 3),
                                      fmt(z$coverage_empirical_center), fmt(z$coverage_realized_weighted_truth))
}
appendix_records <- do.call(rbind, appendix_records)
performance_records <- appendix_records[!(appendix_records$design_id == "C0" & appendix_records$n %in% c(65536, 531441)) &
  !(appendix_records$design_id == "C2" & appendix_records$order_mult == 1), ]
calibration_records <- appendix_records[appendix_records$training_method == "full" &
  ((appendix_records$design_id == "C0" & appendix_records$n %in% c(512, 531441)) |
   (appendix_records$design_id == "C2" & appendix_records$order_mult == .4)), ]
stopifnot(length(arows) == 6L, length(crows) == 6L,
          nrow(performance_records) == 18L, nrow(calibration_records) == 6L)
record("appendix_performance", performance_records, "accepted/performance.csv", case_fields, list(replications = 10000))
record("appendix_calibration", calibration_records,
       "accepted/performance.csv + analysis/raw_metrics.csv", case_fields,
       list(calibration_pairs = 10000, diagnostic_coverages = "original interval with empirical-mean or draw-specific realized-truth center",
            empirical_center_uncertainty = "no binomial uncertainty attached"))
tabular("appendix_performance", "@{}llrrrrr@{}",
        c(row("", "", "\\multicolumn{2}{c}{Full fitting}", "\\multicolumn{3}{c}{Coverage}"),
          "\\cline{3-4}\\cline{5-7}",
          row("Case $(n,\\lambda_s)$", "App.", "Bias", "RMSE", "Full", "Five-fold", "True")), unlist(arows),
        list(metrics = "full-fit bias/RMSE; DNN-jackknife point coverage under three nuisance training routes", replications = 10000))
tabular("appendix_calibration", "@{}llrrrr@{}",
        row("Case $(n,\\lambda_s)$", "App.", "\\shortstack{Mean\\\\SE/SD}", "\\shortstack{RMS\\\\SE/SD}",
            "\\shortstack{Empirical\\\\center}", "\\shortstack{Weighted\\\\truth}"), unlist(crows),
        list(training = "full", localizer = "dnn", variance_route = "jackknife", matched_calibration_pairs = 10000,
             mean_SE_SD = "mean(SE on finite positive interval mask)/SD(estimates on the same mask)",
             rms_SE_SD = "sqrt(mean(SE^2 on that mask))/SD(estimates on the same mask)",
             diagnostic_center_uncertainty = "no binomial interval added"))

pair_cases <- data.frame(design = c("C0", "C2", "C2S"), n = c(2048, 512, 2048), order = c(1, .4, 1))
pair_rows <- pair_records <- list()
for (i in seq_len(nrow(pair_cases))) for (app in apps) {
  design <- pair_cases$design[i]; n <- pair_cases$n[i]; order <- pair_cases$order[i]
  within <- one(W, design_id = design, n = n, order_mult = order, application = app,
                training_method = "full", localization_method = "kernel_gaussian", comparison = "jackknife_vs_direct_squared")
  loss <- one(C, design_id = design, n = n, order_mult = order, application = app,
              training_method = "full", comparison = "dnn_jk_minus_kernel_jk")
  stopifnot(within$paired == 10000, within$excluded == 0)
  pair_records[[length(pair_records) + 1L]] <- data.frame(design_id = design, n = n, order_mult = order, application = app,
    kernel_JK_minus_conventional_coverage = within$coverage_difference, coverage_mcse = within$coverage_difference_mcse,
    kernel_JK_minus_conventional_length = within$length_difference, length_mcse = within$length_difference_mcse,
    DNN_minus_kernel_MSE = loss$MSE_difference, MSE_mcse = loss$MSE_difference_mcse, paired = 10000)
  pair_rows[[length(pair_rows) + 1L]] <- row(case_tex(design, n, order), app_name(app),
    paste0(signed(within$coverage_difference), " (", fmt(within$coverage_difference_mcse), ")"),
    signed(within$length_difference, 3),
    paste0(signed(1000 * loss$MSE_difference, 3), " (", fmt(1000 * loss$MSE_difference_mcse, 3), ")"))
}
pair_records <- do.call(rbind, pair_records)
stopifnot(all(pair_records$DNN_minus_kernel_MSE[pair_records$design_id == "C2"] > 0))
record("appendix_pairs", pair_records, "accepted/variance_pairs.csv + analysis/paired_contrasts.csv",
       c("design_id", "n", "order_mult", "application"), list(sign = "left minus right; MCSE from directly paired draws", paired_count = 10000))
tabular("appendix_pairs", "@{}llrrr@{}",
        row("Case $(n,\\lambda_s)$", "App.", "\\shortstack{Kernel JK$-$conv.\\\\coverage (MCSE)}",
            "\\shortstack{Length\\\\cost}", "\\shortstack{DNN$-$kernel MSE\\\\$10^{-3}$ (MCSE)}"), unlist(pair_rows),
        list(sign = "kernel jackknife minus kernel conventional; point MSE is DNN minus kernel, common within each localizer's covariance routes", units_MSE = "1e-3", paired_count = 10000))

# Within-DNN covariance comparisons: point estimator and nuisance fit held fixed.
dnn_cov_rows <- dnn_pair_rows <- dnn_cov_records <- dnn_pair_records <- list()
for (n in c(512, 65536)) for (app in apps) {
  pp <- do.call(rbind, lapply(c("jackknife", "direct_moment", "direct_squared"),
    function(route) perf("C0", n, app, route = route)))
  stopifnot(all(pp$completed == 10000), max(abs(pp$bias - pp$bias[1])) < 1e-12,
            max(abs(pp$sd - pp$sd[1])) < 1e-12)
  dnn_cov_records[[length(dnn_cov_records) + 1L]] <- pp
  if (app == "cate") dnn_cov_rows[[length(dnn_cov_rows) + 1L]] <- row(integer_tex(n), app_name(app),
    paste0(fmt(pp$coverage), " (", fmt(pp$mean_length, 3), ")"))
  for (route in c("direct_moment", "direct_squared")) {
    pair <- one(W, design_id = "C0", n = n, order_mult = 1, application = app,
      training_method = "full", localization_method = "dnn", comparison = paste0("jackknife_vs_", route))
    other <- pp[pp$variance_method == route, ]
    stopifnot(pair$paired == 10000, pair$excluded == 0,
      abs(pair$coverage_difference - (pp$coverage[1] - other$coverage)) < 1e-12,
      abs(pair$length_difference - (pp$mean_length[1] - other$mean_length)) < 1e-12)
    dnn_pair_records[[length(dnn_pair_records) + 1L]] <- pair
    if (app == "cate") dnn_pair_rows[[length(dnn_pair_rows) + 1L]] <- row(integer_tex(n), app_name(app),
      ifelse(route == "direct_moment", "Moment", "Squared"),
      paste0(signed(pair$coverage_difference), " (", fmt(pair$coverage_difference_mcse), ")"),
      signed(pair$length_difference, 3))
  }
}
dnn_cov_records <- do.call(rbind, dnn_cov_records)
dnn_pair_records <- do.call(rbind, dnn_pair_records)
record("appendix_dnn_covariance", dnn_cov_records, "accepted/performance.csv", case_fields,
  list(design = "C0", training = "full", common_point_estimator = TRUE, replications = 10000))
record("appendix_dnn_covariance_pairs", dnn_pair_records, "accepted/variance_pairs.csv",
  c(cell_keys, "design_id", "n", "order_mult", "comparison"),
  list(sign = "jackknife minus alternative", paired_count = 10000, excluded = 0))
tabular("appendix_dnn_covariance", "@{}llrrr@{}",
  row("$n$", "App.", "Jackknife", "Direct moment", "Direct squared"),
  c(unlist(dnn_cov_rows), "\\hline", row("$n$", "App.", "Alternative", "\\shortstack{Coverage gain\\\\(paired MCSE)}", "Length cost"),
    "\\hline", unlist(dnn_pair_rows)),
  list(entries = "upper panel coverage (mean interval length); lower panel jackknife-minus-alternative coverage difference (paired MCSE) and mean length difference",
       design = "C0", n = c(512, 65536), printed_application = "cate", supporting_applications = apps,
       training = "full", paired_count = 10000))

endpoint_cases <- expand.grid(application = apps, n = c(65536, 531441), stringsAsFactors = FALSE)
endpoint_records <- lapply(seq_len(nrow(endpoint_cases)), function(i) {
  app <- endpoint_cases$application[i]; n <- endpoint_cases$n[i]; pp <- perf("C0", n, app)
  cf <- one(A, design_id = "C0", n = n, application = app, localization_method = "dnn", comparison = "full_vs_CF")
  oracle <- one(A, design_id = "C0", n = n, application = app, localization_method = "dnn", comparison = "full_vs_oracle")
  risk_field <- if (app == "cate") "point_cate_prediction_MSE_sum" else "point_slope_prediction_MSE_sum"
  fresh <- diag("C0", n, app, risk_field)
  projection <- one(J, design_id = "C0", M = pp$M)
  estimated <- diag("C0", n, app, "ratio_target_asymptotic")
  finite <- diag("C0", n, app, "ratio_target_point_oracle")
  control <- diag("C0", n, app, risk_field, training = "population_projection")
  intensive <- diag("C0", n, app, "decomposition_scaled_mean")
  stopifnot(fresh$observed == 10000, estimated$observed == 10000,
            control$observed == ifelse(n == 531441, 100, 1000),
            intensive$observed == ifelse(n == 531441, 0, 200))
  data.frame(application = app, n = n, M = pp$M, K = pp$K, s = pp$s, N = pp$N_reference,
    bias = pp$bias, rmse = pp$rmse,
    point_coverage = pp$coverage, coverage_mcse = pp$coverage_mcse,
    coverage_lower = pp$coverage_lower, coverage_upper = pp$coverage_upper,
    mean_length = pp$mean_length, fresh_point_risk_sum = fresh$mean,
    analytic_projection_risk = if (app == "cate") projection$cate_MSE_sum else projection$slope_MSE_sum,
    q95_full_fivefold_point_SD = cf$q95_absolute_Omega_units,
    q95_full_oracle_point_SD = oracle$q95_absolute_Omega_units,
    q95_absolute_full_oracle = oracle$q95_absolute / sqrt(pp$N_reference),
    mean_estimated_normalized_variance_over_point_Omega = estimated$mean,
    actual_normalized_MC_variance_over_point_Omega = pp$N_reference * pp$sd^2 / omega(app),
    mean_estimated_normalized_variance_over_finite_oracle_array = finite$mean,
    RMS_SE_over_MC_SD = pp$rms_SE_over_SD, main_observed = 10000, fresh_observed = fresh$observed,
    projection_control_observed = control$observed, intensive_observed = intensive$observed,
    intensive_mean = intensive$mean)
})
endpoint_records <- do.call(rbind, endpoint_records)
record("appendix_endpoint", endpoint_records, "accepted/performance.csv + accepted/paired.csv + analysis/diagnostic_moments.csv + analysis/projection_references.csv",
       c("application", "n", "M", "K", "s"), list(reference = "asymptotic point Omega, distinct from finite point-oracle array",
       absolute_gap = "paired q95 of abs(sqrt(N)*(full-oracle)) divided by sqrt(exact N); N constant within these DNN cells"))
endpoint_metrics <- c("bias", "rmse", "point_coverage", "coverage_lower", "coverage_upper", "mean_length", "fresh_point_risk_sum", "analytic_projection_risk",
                      "q95_full_fivefold_point_SD", "q95_full_oracle_point_SD", "q95_absolute_full_oracle",
                      "mean_estimated_normalized_variance_over_point_Omega",
                      "actual_normalized_MC_variance_over_point_Omega", "RMS_SE_over_MC_SD")
endpoint_labels <- c("Bias", "RMSE", "Point coverage", "95\\% coverage lower", "95\\% coverage upper", "Mean interval length", "Fresh-point risk sum", "Analytic projection risk",
                     "\\shortstack[l]{Absolute full/five-fold q95\\\\(asymptotic SD units)}",
                     "\\shortstack[l]{Absolute full/true-nuisance q95\\\\(asymptotic SD units)}",
                     "\\shortstack[l]{Absolute full/true-nuisance q95\\\\(estimate units)}",
                     "$\\Ex{\\widehat\\Omega}/\\Omega_x$", "$N\\Varbi{\\rm MC}{\\widehat a}/\\Omega_x$", "RMS SE/MC SD")
endpoint_digits <- c(4, 4, 4, 4, 4, 3, 5, 5, 3, 3, 4, 3, 3, 3)
endpoint_rows <- vapply(seq_along(endpoint_metrics), function(i)
  row(endpoint_labels[i], fmt(endpoint_records[[endpoint_metrics[i]]], endpoint_digits[i])), character(1))
tabular("appendix_endpoint", "@{}lrrrr@{}",
        row("Metric", "\\shortstack{M2\\\\CATE}", "\\shortstack{M2\\\\Slope}", "\\shortstack{M3\\\\CATE}", "\\shortstack{M3\\\\Slope}"),
        endpoint_rows, list(M2_n = 65536, M3_n = 531441, normalized_gap_units = "asymptotic SD",
             estimated_variance_reference = "asymptotic point constant; not finite oracle-array ratio",
             main_and_risk_count = 10000, projection_control_count = c(M2 = 1000, M3 = 100),
             intensive_count = c(M2 = 200, M3 = 0)))

covariance_records <- do.call(rbind, lapply(c(2048, 65536, 531441), function(n) {
  dd <- do.call(rbind, lapply(c("component", "target", "cross", "transfer", "oracle_cross"), function(field)
    diag("C0", n, "slope", paste0("covariance_", field, "_ratio_point_oracle"))))
  stopifnot(all(dd$observed == 10000), abs(dd$mean[4] - sum(dd$mean[1:3])) < 1e-12)
  data.frame(n = n, M = dd$M[1], application = "slope", component = dd$mean[1], target = dd$mean[2],
             cross = dd$mean[3], transfer = dd$mean[4], oracle_cross = dd$mean[5], observed = 10000)
}))
record("appendix_covariance", covariance_records, "analysis/diagnostic_moments.csv", c("n", "M", "application"),
       list(reference = "mean of per-draw ratios to finite point-oracle score covariance O",
            identity = "T = component + target + signed cross; oracle cross = F-O-T by subtraction"))
covariance_rows <- vapply(seq_len(nrow(covariance_records)), function(i) {
  z <- covariance_records[i, ]; row(integer_tex(z$n), fmt(z$component, 4), fmt(z$target, 5),
                                   signed(z$cross, 5), fmt(z$transfer, 4), signed(z$oracle_cross, 5))
}, character(1))
tabular("appendix_covariance", "@{}rrrrrr@{}",
        row("$n$", "Component", "Target", "Cross", "$T$", "\\shortstack{Oracle cross\\\\$F-O-T$}"), covariance_rows,
        list(application = "slope", design = "C0", training = "full", localizer = "dnn", observed = 10000,
             entries = "means of draw-specific ratios to point-oracle score covariance; cross terms signed"))

# Numerical prose receipt: unrounded selected cells and direct paired contrasts.
anchor_cases <- unique(rbind(appendix_cases, data.frame(design = c("C0", "C2", "C0H", "C2S", "C0H", "C2S"),
  n = c(2048, 2048, 8192, 2048, 2048, 8192), order = 1)))
anchor_perf <- do.call(rbind, lapply(seq_len(nrow(anchor_cases)), function(i) {
  z <- anchor_cases[i, ]; P[P$design_id == z$design & P$n == z$n & P$order_mult == z$order &
    P$training_method %in% main_methods & P$variance_method == "jackknife" & P$localization_method == "dnn", ]
}))
anchor_perf <- anchor_perf[, c(case_fields, "planned", "completed", "intervals", "studentized", "calibration_pairs", "calibration_SD",
  "bias", "bias_mcse", "sd", "rmse", "mean_N", "mean_weighted_radius", "coverage", "coverage_mcse", "coverage_lower", "coverage_upper",
  "mean_length", "mean_SE_over_SD", "rms_SE_over_SD", "student_q025", "student_q975", "student_lower_tail", "student_upper_tail",
  "coverage_empirical_center", "coverage_realized_weighted_truth", "point_t_excess_kurtosis", "mc_z_excess_kurtosis")]
anchor_pairs <- C[C$training_method == "full" & C$localization_method == "dnn" &
  C$comparison %in% c("C2S_minus_baseline", "C0H_minus_baseline", "C1_minus_C0", "M3_endpoint_minus_M2_endpoint",
                      "order_one_minus_point_four", "dnn_jk_minus_kernel_sandwich", "dnn_jk_minus_kernel_jk", "full_minus_pooled_five_fold") &
  ((C$n %in% c(512, 2048, 65536, 531441)) | (C$design_id == "C0H" & C$n == 8192)), ]
anchor_fields <- c("point_cate_prediction_MSE_sum", "point_slope_prediction_MSE_sum", "mean_score_at_point",
  "ratio_target_asymptotic", "ratio_target_point_oracle", "decomposition_scaled_mean", "decomposition_scaled_centered",
  "grid_prediction_radius_squared")
anchor_diagnostics <- D[D$localization_method == "dnn" & D$training_method %in% c("full", "population_projection") &
  D$field %in% anchor_fields & ((D$design_id == "C0" & D$n %in% c(65536, 531441)) |
  (D$design_id == "C1" & D$n %in% c(2048, 4096))), ]
record("prose_performance", anchor_perf, "accepted/performance.csv + analysis/raw_metrics.csv", case_fields)
record("prose_pairs", anchor_pairs, "analysis/paired_contrasts.csv",
       c("design_id", "n", "order_mult", "application", "training_method", "localization_method", "comparison", "left_cell", "right_cell", "left_route", "right_route"))
record("prose_diagnostics", anchor_diagnostics, "analysis/diagnostic_moments.csv",
       c(cell_keys, "design_id", "n", "M", "order_mult", "field", "observed"))
anchors <- list(status = "unrounded_manuscript_prose_anchors", independent_draw_families = 10000,
  sign_convention = "left minus right; normalized paired gaps use absolute values before q95",
  point_Omega = c(cate = omega("cate"), slope = omega("slope")),
  performance = anchor_perf,
  main_comparison = main_records[, c(case_fields, "coverage", "coverage_mcse", "mean_length", "mean_SE_over_SD", "rms_SE_over_SD")],
  paired_contrasts = anchor_pairs[, c("design_id", "n", "order_mult", "application", "comparison", "left_cell", "right_cell", "left_route", "right_route",
    "paired", "excluded", "coverage_difference", "coverage_difference_mcse", "length_difference", "length_difference_mcse", "MSE_difference", "MSE_difference_mcse",
    "mean_gap", "mean_gap_mcse")],
  kernel_pairs = pair_records, dnn_covariance = dnn_cov_records,
  dnn_covariance_pairs = dnn_pair_records, normalized_training_gaps = gaps[, c("cell_id", "n", "M", "application", "comparison", "paired", "N",
    "q95_absolute", "q95_absolute_Omega_units", "q95_absolute_mcse")],
  centering = centers[, c("cell_id", "design_id", "n", "order_mult", "application", "paired", "mean_drift", "mean_drift_mcse", "mean_remainder", "N", "total_bias")],
  endpoint = endpoint_records, covariance = covariance_records,
  projection_references = J[(J$design_id == "C0" & J$M %in% c(2, 3)) | (J$design_id == "C1" & J$M %in% c(1, 2)), ],
  diagnostics = anchor_diagnostics[, c("cell_id", "design_id", "n", "M", "application", "training_method", "field", "observed", "mean", "mcse")],
  gram_events = E[E$design_id %in% c("C0", "C1") & E$n == 4096, ],
  gram_cohorts = G[G$design_id %in% c("C0", "C1") & G$n == 4096, ],
  unavailable_M3_intensive = "observed=0; mean=NA; JSON null means unavailable, never zero")
anchor_path <- file.path(out, "prose_anchors.json")
anchor_text <- jsonlite::toJSON(anchors, pretty = TRUE, auto_unbox = TRUE, digits = NA, na = "null", dataframe = "rows")
if (!check_only) writeLines(anchor_text, anchor_path)

if (check_only) {
  stored <- jsonlite::read_json(file.path(out, "manifest.json"), simplifyVector = TRUE)
  stopifnot(stored$renderer_sha256 == sha(script), !stored$simulation_performed,
            !stored$estimation_performed, !stored$native_RDS_read,
            nrow(stored$output_hashes) == 15L,
            identical(names(stored$output_hashes), c("file", "sha256")))
  for (i in seq_len(nrow(stored$input_hashes))) {
    h <- stored$input_hashes[i, ]; stopifnot(sha(file.path(repo, h$file)) == h$sha256)
  }
  for (i in seq_len(nrow(stored$output_hashes))) {
    h <- stored$output_hashes[i, ]; stopifnot(sha(file.path(repo, h$file)) == h$sha256)
  }
  for (name in names(fragments)) stopifnot(identical(readLines(file.path(tables, name)), fragments[[name]]))
  stopifnot(identical(readLines(anchor_path), strsplit(anchor_text, "\n", fixed = TRUE)[[1]]))
  if (nzchar(Sys.which("qpdf"))) for (name in names(figures)) {
    status <- system2("qpdf", c("--check", shQuote(file.path(out, paste0(name, ".pdf")))), stdout = TRUE, stderr = TRUE)
    if (!is.null(attr(status, "status"))) stop("PDF validation failed: ", name)
  }
  cat("PASS: accepted and analysis hashes, frozen scientific source, canonical rows, matched masks, paired signs, requested selections, eight tabular fragments, prose anchors, 15 output hashes and three PDFs.\n")
} else {
  output_files <- c(file.path(out, paste0(rep(names(figures), each = 2), c(".pdf", ".png"))),
                    file.path(tables, names(fragments)), anchor_path)
  stopifnot(length(output_files) == 15L, all(file.exists(output_files)))
  output_hashes <- data.frame(file = vapply(output_files, relative, character(1)),
                             sha256 = vapply(output_files, sha, character(1)), row.names = NULL)
  input_hashes <- data.frame(file = names(inputs), sha256 = unlist(inputs), row.names = NULL)
  manifest <- list(status = "manuscript_displays_from_accepted_portable_evidence", simulation_performed = FALSE,
    estimation_performed = FALSE, native_RDS_read = FALSE, renderer = relative(script), renderer_sha256 = sha(script),
    figure_renderer_sha256 = if (tables_only) {
      if (!is.null(previous$figure_renderer_sha256)) previous$figure_renderer_sha256 else previous$renderer_sha256
    } else sha(script),
    figure_refresh = if (tables_only) "verified figure bytes preserved after tabular-only formatting change" else "rendered from accepted portable evidence",
    source_commit = system2("git", c("-C", shQuote(repo), "rev-parse", "HEAD"), stdout = TRUE),
    scientific_source_fingerprint = proposal$source_fingerprint, RdistNN_production_commit = proposal$RdistNN_source_commit,
    runtime = list(R = R.version.string, packages = setNames(lapply(c("ggplot2", "patchwork", "jsonlite", "digest"),
      function(pkg) as.character(utils::packageVersion(pkg))), c("ggplot2", "patchwork", "jsonlite", "digest"))),
    input_hashes = input_hashes, output_hashes = output_hashes,
    canonical_rule = "exclude kernel_comparison DNN rows only where cell_id exists in dnn profile; retain unique C0 n8192 DNN",
    join_keys = row_keys, diagnostic_keys = c(cell_keys, "field"),
    pairing = "replication_id in accepted summaries and analysis; direct left-minus-right contrast and paired MCSE; 10000 families, zero exclusions",
    calibration_mask = "successful finite estimates with finite positive SE and available finite interval; SD and SE moments on identical mask; calibration_pairs=10000 in selected main cells",
    diagnostic_centers = "Monte Carlo group-mean and draw-specific realized weighted-truth centers are infeasible diagnostic centers; no independent-binomial uncertainty is added for empirical centering.",
    transformations = list(gaps = "q95(abs(sqrt(N)*(full-comparator)))/sqrt(Omega_point), type 8",
      centering = "sqrt(N)*(mean spatial drift or remainder)/sqrt(Omega_point); remainder includes oracle noise",
      absolute_oracle_gap = "accepted q95_absolute/sqrt(exact DNN N)",
      estimated_variance = "mean ratio_target_asymptotic = mean(N*Vhat_JK)/Omega_point; distinct finite oracle-array reference",
      actual_MC_variance = "exact N * accepted matched-mask Monte Carlo SD^2 / Omega_point",
      covariance = "signed means of draw-specific component/O, target/O, cross/O, T/O, (F-O-T)/O; not ratios of means",
      MSE_units = "1000*(DNN minus kernel MSE) and 1000*paired MCSE"),
    replication_accounting = list(main_per_cell = 10000, independent_draw_families = 10000,
      duplicated_raw_DNN_rows = audit$duplicate_DNN_rows,
      diagnostic_field_count_distribution = as.data.frame(table(observed = D$observed)),
      endpoint_counts = endpoint_records[, c("n", "application", "main_observed", "fresh_observed", "projection_control_observed", "intensive_observed", "intensive_mean")]),
    selections = selections, figures = figures, tables = table_receipts)
  jsonlite::write_json(manifest, file.path(out, "manifest.json"), pretty = TRUE, auto_unbox = TRUE, digits = NA, na = "null", dataframe = "rows")
  cat("Rendered three PDF/PNG figures, eight tabular fragments, prose_anchors.json and manifest.json.\n")
}
