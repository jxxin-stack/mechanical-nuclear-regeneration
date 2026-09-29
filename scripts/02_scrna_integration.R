#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
script_dir <- dirname(normalizePath(script_file, winslash = "/"))
repo <- normalizePath(file.path(script_dir, ".."), winslash = "/")
data_dir <- file.path(repo, "data", "scrna")
config_dir <- file.path(repo, "config")
out_dir <- file.path(repo, "results", "scrna_integration")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(Matrix)
  library(data.table)
  library(ggplot2)
})
stopifnot(
  requireNamespace("scran", quietly = TRUE),
  requireNamespace("irlba", quietly = TRUE),
  requireNamespace("harmony", quietly = TRUE),
  requireNamespace("uwot", quietly = TRUE),
  requireNamespace("pheatmap", quietly = TRUE),
  requireNamespace("ggrepel", quietly = TRUE)
)

seed <- 20260904L
set.seed(seed)

input_files <- c(
  nrc_counts = file.path(data_dir, "nrc_counts.rds"),
  nrc_meta = file.path(data_dir, "nrc_metadata.tsv"),
  pbmc_counts = file.path(data_dir, "pbmc_v31_counts.rds"),
  pbmc_meta = file.path(data_dir, "pbmc_v31_metadata.tsv"),
  markers = file.path(config_dir, "scrna_marker_sets.tsv")
)
missing <- input_files[!file.exists(input_files)]
if (length(missing)) stop("Missing inputs:\n", paste(missing, collapse = "\n"))

as_counts <- function(x) {
  if (!inherits(x, "dgCMatrix")) x <- as(x, "dgCMatrix")
  stopifnot(!is.null(rownames(x)), !is.null(colnames(x)), all(x@x >= 0))
  x
}

collapse_duplicate_genes <- function(x) {
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

read_dataset <- function(count_file, meta_file, type) {
  counts <- collapse_duplicate_genes(as_counts(readRDS(count_file)))
  meta <- fread(meta_file)
  if (!all(c("cell_id", "qc_pass") %in% names(meta))) {
    stop(basename(meta_file), " must contain cell_id and qc_pass")
  }
  idx <- match(colnames(counts), meta$cell_id)
  if (anyNA(idx)) stop("Metadata missing cell IDs for ", type)
  meta <- meta[idx]
  keep <- as.logical(meta$qc_pass)
  counts <- counts[, keep, drop = FALSE]
  meta <- meta[keep]
  meta[, dataset := type]
  list(counts = counts, meta = meta)
}

normalize_cp10k <- function(x) {
  size <- colSums(x)
  if (any(size <= 0)) stop("Zero-library cells remain after QC")
  ans <- x %*% Diagonal(x = 1e4 / size)
  ans@x <- log2(ans@x + 1)
  dimnames(ans) <- dimnames(x)
  ans
}

group_means <- function(x, group) {
  levels <- unique(group)
  ans <- vapply(levels, function(z) rowMeans(x[, group == z, drop = FALSE]), numeric(nrow(x)))
  if (is.null(dim(ans))) ans <- matrix(ans, ncol = 1)
  dimnames(ans) <- list(rownames(x), levels)
  ans
}

nrc <- read_dataset(input_files[["nrc_counts"]], input_files[["nrc_meta"]], "NRC")
pbmc <- read_dataset(input_files[["pbmc_counts"]], input_files[["pbmc_meta"]], "PBMC")
if (!"cluster" %in% names(nrc$meta)) stop("nrc_metadata.tsv must contain cluster")
if (!"cell_type" %in% names(pbmc$meta)) stop("pbmc_v31_metadata.tsv must contain cell_type")

common <- sort(intersect(rownames(nrc$counts), rownames(pbmc$counts)))
nrc_log <- normalize_cp10k(nrc$counts[common, , drop = FALSE])
pbmc_log <- normalize_cp10k(pbmc$counts[common, , drop = FALSE])
non_mt_rp <- !grepl("^(MT-|RPL|RPS)", common)
eligible <- non_mt_rp & rowSums(nrc$counts[common, , drop = FALSE] > 0) >= 10 &
  rowSums(pbmc$counts[common, , drop = FALSE] > 0) >= 10

# Equal weighting prevents the larger reference dataset from determining all HVGs.
var_nrc <- scran::modelGeneVar(nrc_log[eligible, , drop = FALSE])$bio
var_pbmc <- scran::modelGeneVar(pbmc_log[eligible, , drop = FALSE])$bio
genes_eligible <- common[eligible]
rank_score <- rowMeans(cbind(
  rank(-var_nrc, ties.method = "average") / length(var_nrc),
  rank(-var_pbmc, ties.method = "average") / length(var_pbmc)
))
hvg <- genes_eligible[order(rank_score, genes_eligible)][seq_len(min(5000L, length(rank_score)))]
fwrite(data.table(gene = genes_eligible, equal_weight_variance_rank = rank_score),
  file.path(out_dir, "shared_HVG_statistics.tsv"), sep = "\t")

# Harmony is used for joint visualization only.
X <- cbind(nrc_log[hvg, ], pbmc_log[hvg, ])
source <- c(rep("NRC", ncol(nrc_log)), rep("PBMC", ncol(pbmc_log)))
cell_id <- c(paste0("NRC_", nrc$meta$cell_id), paste0("PBMC_", pbmc$meta$cell_id))
set.seed(seed)
pca <- irlba::prcomp_irlba(t(X), n = 50, center = TRUE, scale. = TRUE, tol = 1e-5)
set.seed(seed)
harmony_coordinates <- harmony::RunHarmony(
  data_mat = pca$x,
  meta_data = data.frame(dataset = source),
  vars_use = "dataset",
  theta = 2,
  lambda = 1,
  max_iter = 20,
  ncores = 2,
  verbose = FALSE
)
set.seed(seed)
umap <- uwot::umap(
  harmony_coordinates,
  n_neighbors = 30,
  min_dist = 0.3,
  n_epochs = 500,
  metric = "euclidean",
  seed = seed,
  n_threads = 2,
  n_sgd_threads = 1,
  verbose = FALSE
)

cell_type <- c(rep("Regenerated cells", nrow(nrc$meta)), pbmc$meta$cell_type)
cluster <- c(as.character(nrc$meta$cluster), rep(NA_character_, nrow(pbmc$meta)))
coordinates <- data.table(
  cell_id = cell_id, dataset = source, cell_type = cell_type,
  NRC_cluster = cluster, UMAP1 = umap[, 1], UMAP2 = umap[, 2]
)
fwrite(coordinates, file.path(out_dir, "Fig4c_UMAP_coordinates.tsv"), sep = "\t")

key_types <- c(
  "Regenerated cells", "Classical monocytes", "Intermediate monocytes",
  "Non-classical monocytes", "Myeloid DC", "Plasmacytoid DC",
  "Naive B", "Memory B", "Naive/central-memory CD4 T",
  "Effector/memory CD4 T", "Naive/central-memory CD8 T",
  "Effector/memory CD8 T", "Regulatory T", "MAIT T", "NK",
  "Platelet-enriched / mixed", "Unresolved"
)
palette <- setNames(grDevices::hcl.colors(length(key_types), "Dynamic"), key_types)
palette["Regenerated cells"] <- "#D73027"
palette["Classical monocytes"] <- "#1B9E77"
palette["Intermediate monocytes"] <- "#66A61E"
palette["Non-classical monocytes"] <- "#A6D854"
unknown <- setdiff(unique(coordinates$cell_type), names(palette))
if (length(unknown)) palette <- c(palette, setNames(grDevices::hcl.colors(length(unknown), "Set 2"), unknown))

anchors <- coordinates[, {
  center_x <- median(UMAP1)
  center_y <- median(UMAP2)
  j <- which.min((UMAP1 - center_x)^2 + (UMAP2 - center_y)^2)
  .(UMAP1 = UMAP1[j], UMAP2 = UMAP2[j], cells = .N)
}, by = cell_type]
set.seed(seed)
p_umap <- ggplot(coordinates[sample(.N)], aes(UMAP1, UMAP2, color = cell_type)) +
  geom_point(size = 0.16, alpha = 0.75, stroke = 0) +
  ggrepel::geom_text_repel(
    data = anchors, aes(label = cell_type), color = "black", size = 2.2,
    box.padding = 0.18, point.padding = 0.08, segment.size = 0.2,
    min.segment.length = 0, max.overlaps = Inf, seed = seed
  ) +
  scale_color_manual(values = palette) +
  coord_equal() +
  labs(x = "UMAP 1", y = "UMAP 2") +
  theme_classic(base_size = 7.5) +
  theme(legend.position = "none")
ggsave(file.path(out_dir, "Fig4c_NRC_PBMC_v31_UMAP.pdf"), p_umap,
  width = 80 / 25.4, height = 72 / 25.4, device = cairo_pdf)

# Cell-type correlations use uncorrected log-normalized expression.
clusters <- paste0("C", 0:4)
nrc_mean <- group_means(nrc_log, factor(nrc$meta$cluster, levels = clusters))[, clusters, drop = FALSE]
valid_types <- !is.na(pbmc$meta$cell_type) &
  !(pbmc$meta$cell_type %in% c("Unresolved", "Platelet-enriched / mixed"))
pbmc_mean <- group_means(pbmc_log[, valid_types, drop = FALSE], pbmc$meta$cell_type[valid_types])
correlation <- cor(pbmc_mean[non_mt_rp, , drop = FALSE], nrc_mean[non_mt_rp, , drop = FALSE],
  method = "spearman")
fwrite(as.data.table(correlation, keep.rownames = "cell_type"),
  file.path(out_dir, "Fig4d_Spearman_correlations.tsv"), sep = "\t")

pdf(file.path(out_dir, "Fig4d_Spearman_heatmap.pdf"), width = 80 / 25.4,
  height = 65 / 25.4, useDingbats = FALSE)
pheatmap::pheatmap(
  correlation,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  color = colorRampPalette(c("#D6E7F4", "#FFFFFF", "#E58A83"))(101),
  breaks = seq(min(correlation), max(correlation), length.out = 102),
  border_color = "white",
  display_numbers = TRUE,
  number_format = "%.2f",
  fontsize = 7,
  fontsize_number = 5.5
)
dev.off()

# A priori monocyte and dendritic-cell marker comparison.
marker_config <- fread(input_files[["markers"]])
marker_genes <- unique(marker_config$gene)
comparison_types <- c("Classical monocytes", "Non-classical monocytes", "Myeloid DC")
comparison_types <- comparison_types[comparison_types %in% unique(pbmc$meta$cell_type)]
marker_rows <- list()
for (cl in c("C0", "C1", "C2")) {
  idx <- nrc$meta$cluster == cl
  present <- intersect(marker_genes, rownames(nrc_log))
  marker_rows[[length(marker_rows) + 1L]] <- data.table(
    group = cl, gene = present,
    mean_expression = rowMeans(nrc_log[present, idx, drop = FALSE]),
    percent_expressing = 100 * rowMeans(nrc$counts[present, idx, drop = FALSE] > 0)
  )
}
for (type in comparison_types) {
  idx <- pbmc$meta$cell_type == type
  present <- intersect(marker_genes, rownames(pbmc_log))
  marker_rows[[length(marker_rows) + 1L]] <- data.table(
    group = type, gene = present,
    mean_expression = rowMeans(pbmc_log[present, idx, drop = FALSE]),
    percent_expressing = 100 * rowMeans(pbmc$counts[present, idx, drop = FALSE] > 0)
  )
}
markers <- rbindlist(marker_rows)
markers[, scaled_mean := as.numeric(scale(mean_expression)), by = gene]
markers[, panel := marker_config$panel[match(gene, marker_config$gene)]]
fwrite(markers, file.path(out_dir, "EDFig7_marker_dotplot_source_data.tsv"), sep = "\t")

p_dot <- ggplot(markers, aes(gene, group, size = percent_expressing, color = scaled_mean)) +
  geom_point() +
  facet_grid(. ~ panel, scales = "free_x", space = "free_x") +
  scale_color_gradient2(low = "#4C78A8", mid = "white", high = "#D95F59", midpoint = 0) +
  scale_size(range = c(0.2, 4)) +
  labs(x = NULL, y = NULL, color = "Scaled mean", size = "% expressing") +
  theme_classic(base_size = 7) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "italic"),
    strip.background = element_blank(), strip.text = element_text(size = 6.5))
ggsave(file.path(out_dir, "EDFig7_monocyte_DC_marker_dotplot.pdf"), p_dot,
  width = 145 / 25.4, height = 62 / 25.4, device = cairo_pdf)

capture.output(sessionInfo(), file = file.path(out_dir, "sessionInfo.txt"))
message("Single-cell integration complete: ", out_dir)

