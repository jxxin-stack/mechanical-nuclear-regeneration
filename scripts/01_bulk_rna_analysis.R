#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
script_dir <- dirname(normalizePath(script_file, winslash = "/"))
repo <- normalizePath(file.path(script_dir, ".."), winslash = "/")
data_dir <- file.path(repo, "data", "bulk")
config_dir <- file.path(repo, "config")
gene_set_file <- file.path(config_dir, "bulk_gene_sets.gmt")
out_dir <- file.path(repo, "results", "bulk_rna")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(gridExtra)
})

matrix_file <- file.path(data_dir, "merged_fpkm_qn.tsv")
metadata_file <- file.path(data_dir, "sample_metadata.tsv")
required <- c(matrix_file, metadata_file,
  file.path(config_dir, "gene_panels.tsv"),
  file.path(config_dir, "fig3f_gene_sets.tsv"))
missing <- required[!file.exists(required)]
if (length(missing)) stop("Missing inputs:\n", paste(missing, collapse = "\n"))

write_tsv <- function(x, name) {
  fwrite(as.data.table(x), file.path(out_dir, name), sep = "\t")
}

read_gmt <- function(file) {
  if (!file.exists(file)) return(list())
  lines <- readLines(file, warn = FALSE)
  sets <- lapply(strsplit(lines, "\t", fixed = TRUE), function(x) unique(toupper(x[-c(1, 2)])))
  names(sets) <- vapply(strsplit(lines, "\t", fixed = TRUE), `[`, character(1), 1)
  sets
}

collapse_symbols <- function(symbols, x) {
  symbols <- toupper(trimws(symbols))
  keep <- !is.na(symbols) & nzchar(symbols)
  symbols <- symbols[keep]
  x <- x[keep, , drop = FALSE]
  ord <- order(rowMeans(x, na.rm = TRUE), decreasing = TRUE)
  ord <- ord[!duplicated(symbols[ord])]
  x <- x[ord, , drop = FALSE]
  rownames(x) <- symbols[ord]
  x
}

row_zscore <- function(x) {
  z <- t(scale(t(as.matrix(x))))
  z[!is.finite(z)] <- 0
  pmax(-2.5, pmin(2.5, z))
}

gsea_curve <- function(ranked_scores, gene_set) {
  hit <- names(ranked_scores) %in% gene_set
  n_hit <- sum(hit)
  n <- length(ranked_scores)
  if (n_hit < 5L || n_hit >= n) return(NULL)
  increment <- rep(-1 / (n - n_hit), n)
  weights <- abs(ranked_scores[hit])
  if (!sum(weights)) weights[] <- 1
  increment[hit] <- weights / sum(weights)
  running <- cumsum(increment)
  if (max(running) >= abs(min(running))) {
    es <- max(running)
    peak <- which.max(running)
    leading <- names(ranked_scores)[hit & seq_len(n) <= peak]
  } else {
    es <- min(running)
    peak <- which.min(running)
    leading <- names(ranked_scores)[hit & seq_len(n) >= peak]
  }
  list(ES = es, size = n_hit, leading = leading)
}

run_gsea <- function(scores, gene_sets, nperm = 500L, seed = 20260830L,
                     positive = "NRC", negative = "Fresh") {
  scores <- sort(scores[is.finite(scores)], decreasing = TRUE)
  scores <- scores[!duplicated(names(scores))]
  set.seed(seed)
  rows <- lapply(names(gene_sets), function(set_name) {
    genes <- intersect(unique(toupper(gene_sets[[set_name]])), names(scores))
    observed <- gsea_curve(scores, genes)
    if (is.null(observed)) return(NULL)
    null <- replicate(nperm, {
      random_genes <- names(scores)[sample.int(length(scores), observed$size)]
      gsea_curve(scores, random_genes)$ES
    })
    if (observed$ES >= 0) {
      same_sign <- null[null >= 0]
      p <- (1 + sum(same_sign >= observed$ES)) / (1 + length(same_sign))
      direction <- positive
    } else {
      same_sign <- null[null < 0]
      p <- (1 + sum(same_sign <= observed$ES)) / (1 + length(same_sign))
      direction <- negative
    }
    data.table(
      gene_set = set_name,
      size = observed$size,
      ES = observed$ES,
      NES = observed$ES / mean(abs(same_sign)),
      nominal_p = p,
      direction = direction,
      leading_edge = paste(observed$leading, collapse = ";")
    )
  })
  ans <- rbindlist(rows, fill = TRUE)
  ans[, FDR := p.adjust(nominal_p, method = "BH")]
  ans[order(-NES)]
}

raw <- fread(matrix_file, check.names = FALSE)
meta <- fread(metadata_file)
required_meta <- c("sample_id", "state", "donor", "cell_type", "reference_panel")
if (!all(required_meta %in% names(meta))) {
  stop("sample_metadata.tsv must contain: ", paste(required_meta, collapse = ", "))
}
if (!"Symbol" %in% names(raw)) stop("merged_fpkm_qn.tsv must contain a Symbol column")
sample_ids <- intersect(meta$sample_id, names(raw))
if (!length(sample_ids)) stop("No metadata sample IDs match the expression matrix")
meta <- meta[match(sample_ids, sample_id)]
expr <- collapse_symbols(raw$Symbol, as.matrix(raw[, ..sample_ids]))
storage.mode(expr) <- "numeric"
expr[!is.finite(expr)] <- 0

state_indices <- lapply(c("Fresh", "NRC"), function(s) which(meta$state == s))
names(state_indices) <- c("Fresh", "NRC")
if (any(lengths(state_indices) == 0L)) stop("Fresh and NRC samples are required")

donor_mean <- function(state) {
  idx <- state_indices[[state]]
  donors <- unique(meta$donor[idx])
  ans <- vapply(donors, function(d) {
    rowMeans(expr[, idx[meta$donor[idx] == d], drop = FALSE])
  }, numeric(nrow(expr)))
  if (is.null(dim(ans))) ans <- matrix(ans, ncol = 1)
  dimnames(ans) <- list(rownames(expr), donors)
  ans
}

fresh_donor <- donor_mean("Fresh")
nrc_donor <- donor_mean("NRC")
fresh_mean <- rowMeans(fresh_donor)
nrc_mean <- rowMeans(nrc_donor)
log2fc <- log2((nrc_mean + 1) / (fresh_mean + 1))
rank_score <- log2fc * log2(pmax(nrc_mean, fresh_mean) + 1)
transition <- data.table(
  gene = rownames(expr), Fresh_mean = fresh_mean, NRC_mean = nrc_mean,
  log2FC_NRC_vs_Fresh = log2fc, rank_score = rank_score,
  NRC_higher_in_all_donors = apply(nrc_donor, 1, min) > apply(fresh_donor, 1, max),
  Fresh_higher_in_all_donors = apply(fresh_donor, 1, min) > apply(nrc_donor, 1, max)
)
setorder(transition, -abs(rank_score))
write_tsv(transition, "Fresh_to_NRC_ranked_genes.tsv")

gene_sets <- read_gmt(gene_set_file)
if (length(gene_sets)) {
  scores <- setNames(rank_score, rownames(expr))
  scores <- scores[pmax(fresh_mean, nrc_mean) > 1]
  gsea <- run_gsea(scores, gene_sets)
  write_tsv(gsea, "Fresh_to_NRC_GSEA.tsv")
  plot_data <- gsea[order(NES)]
  plot_data[, gene_set := factor(gene_set, levels = gene_set)]
  p_gsea <- ggplot(plot_data, aes(NES, gene_set, fill = direction)) +
    geom_col(width = 0.72) +
    geom_vline(xintercept = 0, linewidth = 0.25) +
    scale_fill_manual(values = c(Fresh = "#4A8BC2", NRC = "#D86A61")) +
    labs(x = "Normalized enrichment score", y = NULL, title = "GSEA: Fresh to NRC") +
    theme_classic(base_size = 7) +
    theme(legend.position = "none", plot.title = element_text(hjust = 0.5))
  ggsave(file.path(out_dir, "Fresh_to_NRC_GSEA.pdf"), p_gsea,
    width = 64 / 25.4, height = 48 / 25.4, device = cairo_pdf)
}

# Gene-expression panels used across the main figures.
panel_config <- fread(file.path(config_dir, "gene_panels.tsv"))
aliases <- list(OCT4 = c("POU5F1", "OCT4"), H3F3A = c("H3F3A", "H3-3A"))
cell_order <- unique(meta$cell_type)
plot_gene <- function(gene, show_y = FALSE) {
  candidates <- if (gene %in% names(aliases)) aliases[[gene]] else gene
  query <- candidates[candidates %in% rownames(expr)][1]
  if (is.na(query)) query <- gene
  if (!query %in% rownames(expr)) {
    warning("Gene not found: ", gene)
    return(NULL)
  }
  d <- data.table(FPKM = expr[query, ], cell_type = factor(meta$cell_type, levels = cell_order))
  ggplot(d, aes(cell_type, FPKM, fill = cell_type)) +
    geom_boxplot(outlier.shape = NA, linewidth = 0.25) +
    geom_point(position = position_jitter(width = 0.15), size = 0.65, alpha = 0.45) +
    labs(x = NULL, y = if (show_y) "FPKM" else NULL,
      title = bquote(italic(.(gene)))) +
    theme_classic(base_size = 6) +
    theme(
      legend.position = "none", plot.title = element_text(hjust = 0.5, size = 7),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 4.5),
      axis.text.y = element_text(size = 4.5), plot.margin = margin(1, 1, 1, 1, "mm")
    )
}

for (fig in unique(panel_config$figure)) {
  genes <- panel_config[figure == fig, gene]
  plots <- Filter(Negate(is.null), lapply(seq_along(genes), function(i) plot_gene(genes[i], i == 1L)))
  if (!length(plots)) next
  ncol_panel <- if (length(plots) > 7L) 6L else length(plots)
  nrow_panel <- ceiling(length(plots) / ncol_panel)
  width_mm <- if (length(plots) >= 7L) 170 else 170 / 6 * length(plots)
  pdf(file.path(out_dir, paste0(fig, "_gene_panel.pdf")), width = width_mm / 25.4,
    height = (28 * nrow_panel) / 25.4, useDingbats = FALSE)
  do.call(grid.arrange, c(plots, ncol = ncol_panel))
  dev.off()
}

# Fig. 3f gene-level heatmaps.
fig3f <- fread(file.path(config_dir, "fig3f_gene_sets.tsv"))
study_idx <- which(meta$state %in% c("Fresh", "NRC"))
study_expr <- log2(expr[, study_idx, drop = FALSE] + 1)
study_labels <- paste(meta$state[study_idx], meta$donor[study_idx], sep = "-")
for (panel_name in unique(fig3f$panel)) {
  cfg <- fig3f[panel == panel_name]
  genes <- unique(cfg$gene[cfg$gene %in% rownames(study_expr)])
  if (!length(genes)) next
  z <- row_zscore(study_expr[genes, , drop = FALSE])
  d <- as.data.table(as.table(z))
  setnames(d, c("gene", "sample", "z"))
  d[, sample := factor(sample, levels = colnames(z), labels = study_labels)]
  d[, gene := factor(gene, levels = rev(genes))]
  p_heat <- ggplot(d, aes(sample, gene, fill = z)) +
    geom_tile() +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
      midpoint = 0, limits = c(-2.5, 2.5), name = "z-score") +
    labs(x = NULL, y = NULL) +
    theme_classic(base_size = 7) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
      axis.text.y = element_text(face = "italic"))
  ggsave(file.path(out_dir, paste0("Fig3f_", panel_name, ".pdf")), p_heat,
    width = 80 / 25.4, height = max(55, 1.8 * length(genes)) / 25.4,
    device = cairo_pdf)
  source_table <- as.data.table(z, keep.rownames = "gene")
  write_tsv(source_table, paste0("Fig3f_", panel_name, "_source_data.tsv"))
}

# Correlations with each bulk reference panel are calculated separately.
study_mean <- cbind(Fresh = fresh_mean, NRC = nrc_mean)
for (panel in setdiff(unique(meta$reference_panel), c("study", "Study", NA, ""))) {
  idx <- which(meta$reference_panel == panel & !is.na(meta$cell_type))
  if (!length(idx)) next
  types <- unique(meta$cell_type[idx])
  reference_mean <- vapply(types, function(type) {
    rowMeans(expr[, idx[meta$cell_type[idx] == type], drop = FALSE])
  }, numeric(nrow(expr)))
  colnames(reference_mean) <- types
  keep <- rowMeans(cbind(study_mean, reference_mean)) > 1 &
    !grepl("^(MT-|RPL|RPS)", rownames(expr))
  correlation <- cor(log2(reference_mean[keep, , drop = FALSE] + 1),
    log2(study_mean[keep, , drop = FALSE] + 1), method = "spearman")
  panel_id <- gsub("[^A-Za-z0-9]+", "_", panel)
  write_tsv(data.table(cell_type = rownames(correlation), correlation),
    paste0("bulk_correlation_", panel_id, ".tsv"))

  # Blood-reference GSEA uses an equal-weight mean across reference cell types.
  if (length(gene_sets)) {
    blood_mean <- rowMeans(reference_mean)
    comparable <- blood_mean > 0 & pmax(fresh_mean, nrc_mean, blood_mean) > 1 &
      !grepl("^MT-", rownames(expr))
    for (state in c("Fresh", "NRC")) {
      target <- if (state == "Fresh") fresh_mean else nrc_mean
      blood_score <- log2((target[comparable] + 1) / (blood_mean[comparable] + 1))
      names(blood_score) <- rownames(expr)[comparable]
      blood_gsea <- run_gsea(
        blood_score, gene_sets, nperm = 1000L,
        seed = if (state == "Fresh") 20260831L else 20260832L,
        positive = state, negative = "Blood"
      )
      write_tsv(blood_gsea, paste0(state, "_vs_", panel_id, "_GSEA.tsv"))
    }
  }
}

capture.output(sessionInfo(), file = file.path(out_dir, "sessionInfo.txt"))
message("Bulk RNA analysis complete: ", out_dir)
