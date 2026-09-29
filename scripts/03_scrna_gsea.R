#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
script_dir <- dirname(normalizePath(script_file, winslash = "/"))
repo <- normalizePath(file.path(script_dir, ".."), winslash = "/")
data_dir <- file.path(repo, "data", "scrna")
gene_set_dir <- file.path(repo, "data", "gene_sets")
config_dir <- file.path(repo, "config")
out_dir <- file.path(repo, "results", "scrna_gsea")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(Matrix)
  library(data.table)
  library(ggplot2)
  library(fgsea)
})

required <- c(
  file.path(data_dir, "nrc_counts.rds"),
  file.path(data_dir, "nrc_metadata.tsv"),
  file.path(data_dir, "pbmc_v31_counts.rds"),
  file.path(data_dir, "pbmc_v31_metadata.tsv"),
  file.path(gene_set_dir, "ReactomePathways.gmt"),
  file.path(config_dir, "macrophage_signature.tsv"),
  file.path(config_dir, "fig5b_pathways.txt")
)
missing <- required[!file.exists(required)]
if (length(missing)) stop("Missing inputs:\n", paste(missing, collapse = "\n"))

as_counts <- function(x) {
  if (!inherits(x, "dgCMatrix")) x <- as(x, "dgCMatrix")
  genes <- toupper(trimws(rownames(x)))
  keep <- !is.na(genes) & nzchar(genes)
  x <- x[keep, , drop = FALSE]
  genes <- genes[keep]
  levels <- unique(genes)
  map <- sparseMatrix(i = match(genes, levels), j = seq_along(genes), x = 1,
    dims = c(length(levels), length(genes)))
  ans <- as(map %*% x, "dgCMatrix")
  dimnames(ans) <- list(levels, colnames(x))
  ans
}

read_dataset <- function(count_name, meta_name) {
  x <- as_counts(readRDS(file.path(data_dir, count_name)))
  meta <- fread(file.path(data_dir, meta_name))
  idx <- match(colnames(x), meta$cell_id)
  if (anyNA(idx)) stop("Metadata does not cover every cell in ", count_name)
  meta <- meta[idx]
  keep <- as.logical(meta$qc_pass)
  list(counts = x[, keep, drop = FALSE], meta = meta[keep])
}

normalize_cp10k <- function(x) {
  ans <- x %*% Diagonal(x = 1e4 / colSums(x))
  ans@x <- log2(ans@x + 1)
  dimnames(ans) <- dimnames(x)
  ans
}

nrc <- read_dataset("nrc_counts.rds", "nrc_metadata.tsv")
pbmc <- read_dataset("pbmc_v31_counts.rds", "pbmc_v31_metadata.tsv")
if (!"cluster" %in% names(nrc$meta)) stop("NRC metadata requires cluster")
if (!"cell_type" %in% names(pbmc$meta)) stop("PBMC metadata requires cell_type")

common <- sort(intersect(rownames(nrc$counts), rownames(pbmc$counts)))
nrc$counts <- nrc$counts[common, , drop = FALSE]
pbmc$counts <- pbmc$counts[common, , drop = FALSE]
nrc_log <- normalize_cp10k(nrc$counts)
pbmc_log <- normalize_cp10k(pbmc$counts)

monocyte_types <- intersect(
  c("Classical monocytes", "Non-classical monocytes"),
  unique(pbmc$meta$cell_type)
)
if (!length(monocyte_types)) stop("No PBMC monocyte annotations were found")

monocyte_mean_by_type <- vapply(monocyte_types, function(type) {
  rowMeans(pbmc_log[, pbmc$meta$cell_type == type, drop = FALSE])
}, numeric(length(common)))
monocyte_detect_by_type <- vapply(monocyte_types, function(type) {
  rowMeans(pbmc$counts[, pbmc$meta$cell_type == type, drop = FALSE] > 0)
}, numeric(length(common)))
monocyte_mean <- rowMeans(monocyte_mean_by_type)
monocyte_detect <- rowMeans(monocyte_detect_by_type)

reactome <- fgsea::gmtPathways(file.path(gene_set_dir, "ReactomePathways.gmt"))
names(reactome) <- make.unique(paste0("Reactome: ", names(reactome)))
signature <- fread(file.path(config_dir, "macrophage_signature.tsv"))
custom <- split(toupper(signature$gene), signature$programme)
pathways <- c(custom, reactome)

clusters <- paste0("C", 0:4)
all_results <- list()
all_ranks <- list()
for (cluster in clusters) {
  idx <- nrc$meta$cluster == cluster
  if (!any(idx)) next
  nrc_mean <- rowMeans(nrc_log[, idx, drop = FALSE])
  nrc_detect <- rowMeans(nrc$counts[, idx, drop = FALSE] > 0)
  keep <- !grepl("^(MT-|RPL|RPS)", common) &
    (nrc_detect >= 0.05 | monocyte_detect >= 0.05)
  ranking <- nrc_mean[keep] - monocyte_mean[keep]
  names(ranking) <- common[keep]
  ranking <- sort(ranking, decreasing = TRUE)
  gsea <- suppressWarnings(fgseaMultilevel(
    pathways = pathways,
    stats = ranking,
    minSize = 10,
    maxSize = 500,
    eps = 0
  ))
  gsea[, cluster := cluster]
  all_results[[cluster]] <- gsea
  all_ranks[[cluster]] <- data.table(cluster = cluster, gene = names(ranking), score = ranking)
}

results <- rbindlist(all_results, fill = TRUE)
ranks <- rbindlist(all_ranks)
fwrite(results, file.path(out_dir, "Fig5b_all_GSEA_results.tsv.gz"), sep = "\t")
fwrite(ranks, file.path(out_dir, "Fig5b_gene_rankings.tsv.gz"), sep = "\t")

display_pathways <- readLines(file.path(config_dir, "fig5b_pathways.txt"), warn = FALSE)
display <- results[pathway %in% display_pathways]
missing_pathways <- setdiff(display_pathways, display$pathway)
if (length(missing_pathways)) warning("Displayed pathways absent from GSEA output: ",
  paste(missing_pathways, collapse = ", "))
display[, cluster := factor(cluster, levels = clusters)]
display[, short_pathway := sub("^Reactome: ", "", pathway)]
display[pathway == "Macrophage core", short_pathway := "Macrophage-associated"]
display[, short_pathway := vapply(short_pathway, function(x) {
  paste(strwrap(x, width = 31), collapse = "\n")
}, character(1))]
path_order <- display_pathways[display_pathways %in% display$pathway]
labels <- unique(display[, .(pathway, short_pathway)])
labels[, order := match(pathway, path_order)]
setorder(labels, order)
display[, short_pathway := factor(short_pathway, levels = rev(labels$short_pathway))]
display[, significance := fifelse(
  pval < 0.001, "***",
  fifelse(pval < 0.01, "**", fifelse(pval < 0.05, "*", ""))
)]
limit <- max(2, min(3, max(abs(display$NES), na.rm = TRUE)))

p <- ggplot(display, aes(cluster, short_pathway,
    fill = pmax(-limit, pmin(limit, NES)))) +
  geom_tile(color = "white", linewidth = 0.3) +
  geom_text(aes(label = significance), size = 2.3) +
  scale_fill_gradient2(
    low = "#83ADD0", mid = "#FAFAFA", high = "#D97872",
    midpoint = 0, limits = c(-limit, limit), name = "NES"
  ) +
  labs(x = NULL, y = NULL, title = "Enrichment relative to PBMC monocytes") +
  theme_classic(base_size = 6.2) +
  theme(
    axis.text.x = element_text(size = 5.8),
    axis.text.y = element_text(size = 5.0, lineheight = 0.88),
    plot.title = element_text(size = 7.2, hjust = 0.5),
    legend.title = element_text(size = 5.6),
    legend.text = element_text(size = 5.2),
    legend.key.height = grid::unit(4.5, "mm"),
    plot.margin = margin(2, 2, 2, 2)
  )
ggsave(file.path(out_dir, "Fig5b.pdf"), p,
  width = 65 / 25.4, height = 48 / 25.4, device = cairo_pdf)
ggsave(file.path(out_dir, "Fig5b.png"), p,
  width = 65 / 25.4, height = 48 / 25.4, dpi = 600, bg = "white")
fwrite(display, file.path(out_dir, "Fig5b_source_data.tsv"), sep = "\t")

capture.output(sessionInfo(), file = file.path(out_dir, "sessionInfo.txt"))
message("Single-cell GSEA complete: ", out_dir)
