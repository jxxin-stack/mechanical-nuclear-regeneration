#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
script_dir <- dirname(normalizePath(script_file, winslash = "/"))
repo <- normalizePath(file.path(script_dir, ".."), winslash = "/")
qc_file <- file.path(repo, "data", "atac", "library_qc.tsv")
out_dir <- file.path(repo, "results", "atac_peak_metrics")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

if (!file.exists(qc_file)) stop("Missing input: ", qc_file)
qc <- fread(qc_file)
required <- c("sample_id", "state", "donor", "peak_count", "frip")
if (!all(required %in% names(qc))) stop("library_qc.tsv must contain: ", paste(required, collapse = ", "))
qc <- qc[state %in% c("Fresh", "NRC")]
qc[, state := factor(state, levels = c("Fresh", "NRC"))]
if (any(table(qc$state) < 2L)) stop("At least two libraries per state are required")

state_colors <- c(Fresh = "#3B82B8", NRC = "#D6604D")
format_p <- function(p) {
  if (p < 0.001) formatC(p, format = "e", digits = 1) else formatC(p, format = "f", digits = 3)
}

plot_metric <- function(column, y_label, output_name, scale = 1) {
  d <- copy(qc)
  d[, value := get(column) / scale]
  p_value <- t.test(value ~ state, data = d, paired = FALSE)$p.value
  ymax <- max(d$value, na.rm = TRUE)
  p <- ggplot(d, aes(state, value, fill = state)) +
    geom_boxplot(width = 0.55, outlier.shape = NA, linewidth = 0.35, alpha = 0.25) +
    geom_point(position = position_jitter(width = 0.09), size = 1.2, alpha = 0.8) +
    annotate("text", x = 1.5, y = ymax * 1.08,
      label = paste0("Welch P = ", format_p(p_value)), size = 2.2) +
    scale_fill_manual(values = state_colors) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.20))) +
    labs(x = NULL, y = y_label) +
    theme_classic(base_size = 7) +
    theme(legend.position = "none")
  ggsave(file.path(out_dir, paste0(output_name, ".pdf")), p,
    width = 38 / 25.4, height = 32 / 25.4, device = cairo_pdf)
  data.table(metric = column, p_value = p_value)
}

statistics <- rbindlist(list(
  plot_metric("peak_count", "Accessible regions (thousands)", "Fig2c_peak_count", 1000),
  plot_metric("frip", "FRiP", "Fig2c_FRiP")
))
fwrite(qc, file.path(out_dir, "Fig2c_source_data.tsv"), sep = "\t")
fwrite(statistics, file.path(out_dir, "Fig2c_statistics.tsv"), sep = "\t")
capture.output(sessionInfo(), file = file.path(out_dir, "sessionInfo.txt"))
message("ATAC peak metrics complete: ", out_dir)

