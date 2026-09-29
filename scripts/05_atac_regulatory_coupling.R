#!/usr/bin/env Rscript

# Core Fresh-NRC ATAC/RNA analysis for Extended Data Fig. 3.
# Existing Duren-style openness scores are used without distribution-forcing
# normalization. Technical libraries are averaged within independent donors.

options(stringsAsFactors = FALSE)

script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
script_dir <- dirname(normalizePath(script_file, winslash = "/"))
repo <- normalizePath(file.path(script_dir, ".."), winslash = "/")
atac_dir <- file.path(repo, "data", "atac")
out_dir <- file.path(repo, "results", "atac_regulatory_coupling")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

atac_file <- file.path(atac_dir, "union_peak_openness.tsv")
sample_file <- file.path(atac_dir, "sample_metadata.tsv")
blacklist_file <- file.path(atac_dir, "blacklisted_peaks.txt")
fresh_network_file <- file.path(atac_dir, "fresh_network.tsv")
nrc_network_file <- file.path(atac_dir, "nrc_network.tsv")
rna_file <- file.path(atac_dir, "rna_fpkm_qn.tsv")

required <- c(
  atac_file, sample_file, blacklist_file,
  fresh_network_file, nrc_network_file, rna_file
)
missing <- required[!file.exists(required)]
if (length(missing)) stop("Missing inputs:\n", paste(missing, collapse = "\n"))

fresh_colour <- "#3B82B8"
nrc_colour <- "#D6604D"
stable_colour <- "#A6A6A6"
shared_colour <- "#7B6FA8"
dark_grey <- "#555555"
seed <- 20260917L
n_perm <- 10000L
atac_cutoff <- 1
minimum_open_openness <- 2
stable_cutoff <- 0.25
rna_cutoff <- 1
mean_fpkm_cutoff <- 2
minimum_re_coverage <- 0.80
minimum_linked_peaks <- 2L
nrc_opening_link_fraction <- 0.75

format_n <- function(x) format(x, big.mark = ",", scientific = FALSE, trim = TRUE)

save_plot <- function(stem, width_mm, height_mm, draw_fun) {
  pdf(
    file.path(out_dir, paste0(stem, ".pdf")),
    width = width_mm / 25.4, height = height_mm / 25.4,
    family = "Helvetica", useDingbats = FALSE
  )
  draw_fun()
  dev.off()

  png(
    file.path(out_dir, paste0(stem, ".png")),
    width = width_mm, height = height_mm, units = "mm", res = 600,
    type = if (capabilities("cairo")) "cairo" else "windows",
    bg = "white"
  )
  draw_fun()
  dev.off()

  jpeg(
    file.path(out_dir, paste0(stem, ".jpg")),
    width = width_mm, height = height_mm, units = "mm", res = 600,
    quality = 100, bg = "white"
  )
  draw_fun()
  dev.off()
}

parse_regions <- function(ids) {
  data.frame(
    ID = ids,
    Chr = sub("_([0-9]+)_([0-9]+)$", "", ids),
    Start = suppressWarnings(as.numeric(sub("^.*_([0-9]+)_([0-9]+)$", "\\1", ids))),
    End = suppressWarnings(as.numeric(sub("^.*_([0-9]+)_([0-9]+)$", "\\2", ids))),
    stringsAsFactors = FALSE
  )
}

build_peak_index <- function(ids) {
  regions <- parse_regions(ids)
  valid <- is.finite(regions$Start) & is.finite(regions$End) &
    regions$End > regions$Start
  regions <- regions[valid, , drop = FALSE]
  regions$Peak_row <- which(valid)
  by_chr <- split(seq_len(nrow(regions)), regions$Chr)
  by_chr <- lapply(by_chr, function(i) i[order(regions$Start[i])])
  list(regions = regions, by_chr = by_chr)
}

map_re_to_peak <- function(re_ids, peak_index) {
  re <- parse_regions(re_ids)
  peak_row <- rep(NA_integer_, nrow(re))
  coverage <- rep(0, nrow(re))
  peaks <- peak_index$regions

  for (chr in intersect(unique(re$Chr), names(peak_index$by_chr))) {
    ri <- which(re$Chr == chr & is.finite(re$Start) & is.finite(re$End))
    pi <- peak_index$by_chr[[chr]]
    ps <- peaks$Start[pi]
    pe <- peaks$End[pi]
    for (j in ri) {
      left <- findInterval(re$Start[j], ps)
      candidates <- unique(c(left - 1L, left, left + 1L, left + 2L))
      candidates <- candidates[candidates >= 1L & candidates <= length(pi)]
      if (!length(candidates)) next
      overlap <- pmax(
        0,
        pmin(re$End[j], pe[candidates]) - pmax(re$Start[j], ps[candidates])
      )
      best <- which.max(overlap)
      if (overlap[best] > 0) {
        coverage[j] <- overlap[best] / (re$End[j] - re$Start[j])
        peak_row[j] <- peaks$Peak_row[pi[candidates[best]]]
      }
    }
  }
  data.frame(RE = re_ids, Peak_row = peak_row, RE_coverage = coverage)
}

read_network_pairs <- function(path, state) {
  net <- read.delim(
    path, header = TRUE, sep = "\t", quote = "", comment.char = "",
    check.names = FALSE
  )
  if (!all(c("TG", "FDR", "REs") %in% names(net))) {
    stop(basename(path), " must contain TG, FDR and REs columns")
  }
  net <- net[
    is.finite(net$FDR) & net$FDR < 0.05 &
      !is.na(net$TG) & nzchar(net$TG) & !is.na(net$REs),
    c("TG", "REs"), drop = FALSE
  ]
  net$TG <- toupper(trimws(net$TG))
  net <- unique(net)
  re_lists <- strsplit(net$REs, ";", fixed = TRUE)
  gene <- rep(net$TG, lengths(re_lists))
  re <- trimws(unlist(re_lists, use.names = FALSE))
  keep <- nzchar(gene) & nzchar(re)
  pairs <- unique(data.frame(
    Gene = gene[keep], RE = re[keep], State = state,
    stringsAsFactors = FALSE
  ))
  pairs
}

# -----------------------------------------------------------------------------
# Common Fresh-NRC peak universe and state-level openness.
# -----------------------------------------------------------------------------

message("Reading the common Fresh-NRC openness matrix without quantile normalization...")
atac_table <- read.delim(atac_file, header = TRUE, sep = "\t", quote = "",
  comment.char = "", check.names = FALSE)
sample_meta <- read.delim(sample_file, header = TRUE, sep = "\t",
  stringsAsFactors = FALSE, check.names = FALSE)
if (!all(c("sample_id", "state", "donor") %in% names(sample_meta))) {
  stop("sample_metadata.tsv must contain sample_id, state and donor")
}
sample_names <- intersect(sample_meta$sample_id, colnames(atac_table))
if (!length(sample_names)) stop("No ATAC sample IDs match the openness matrix")
sample_meta <- sample_meta[match(sample_names, sample_meta$sample_id), , drop = FALSE]
peak_ids <- trimws(atac_table[[1]])
atac <- as.matrix(atac_table[, sample_names, drop = FALSE])
storage.mode(atac) <- "double"
colnames(atac) <- sample_names
rm(atac_table)

blacklist <- unique(trimws(readLines(blacklist_file, warn = FALSE)))
keep <- grepl("^chr[^_]+_[0-9]+_[0-9]+$", peak_ids) &
  !peak_ids %in% blacklist & !grepl("^chrM_", peak_ids, ignore.case = TRUE)
peak_ids <- peak_ids[keep]
atac <- atac[keep, , drop = FALSE]
rm(keep)

fresh_donors <- unique(sample_meta$donor[sample_meta$state == "Fresh"])
nrc_donors <- unique(sample_meta$donor[sample_meta$state == "NRC"])
if (length(fresh_donors) != 2L || length(nrc_donors) != 2L) {
  stop("This analysis expects two independent donors per state")
}
fresh_d1_cols <- which(sample_meta$state == "Fresh" & sample_meta$donor == fresh_donors[1])
fresh_d2_cols <- which(sample_meta$state == "Fresh" & sample_meta$donor == fresh_donors[2])
nrc_d1_cols <- which(sample_meta$state == "NRC" & sample_meta$donor == nrc_donors[1])
nrc_d2_cols <- which(sample_meta$state == "NRC" & sample_meta$donor == nrc_donors[2])
if (any(lengths(list(fresh_d1_cols, fresh_d2_cols, nrc_d1_cols, nrc_d2_cols)) == 0L)) {
  stop("Every state-by-donor group must contain at least one ATAC library")
}

# The input values are local-background-normalized openness scores. Work on the
# log scale and average technical libraries within each biological donor.
atac_log <- log2(atac + 1)
fresh_d1_atac <- rowMeans(atac_log[, fresh_d1_cols, drop = FALSE])
fresh_d2_atac <- rowMeans(atac_log[, fresh_d2_cols, drop = FALSE])
nrc_d1_atac <- rowMeans(atac_log[, nrc_d1_cols, drop = FALSE])
nrc_d2_atac <- rowMeans(atac_log[, nrc_d2_cols, drop = FALSE])
fresh_atac <- (fresh_d1_atac + fresh_d2_atac) / 2
nrc_atac <- (nrc_d1_atac + nrc_d2_atac) / 2
atac_change <- nrc_atac - fresh_atac
for (object_name in c(
  "fresh_d1_atac", "fresh_d2_atac", "nrc_d1_atac", "nrc_d2_atac",
  "fresh_atac", "nrc_atac", "atac_change"
)) {
  assign(object_name, setNames(get(object_name), peak_ids))
}

# With only two independent donors per state, call robust changes by effect size
# plus complete separation of the two donor-level state distributions.
fresh_separated <- pmin(fresh_d1_atac, fresh_d2_atac) >
  pmax(nrc_d1_atac, nrc_d2_atac)
nrc_separated <- pmin(nrc_d1_atac, nrc_d2_atac) >
  pmax(fresh_d1_atac, fresh_d2_atac)
minimum_open_log <- log2(minimum_open_openness + 1)
fresh_absolutely_open <- fresh_d1_atac > minimum_open_log &
  fresh_d2_atac > minimum_open_log
nrc_absolutely_open <- nrc_d1_atac > minimum_open_log &
  nrc_d2_atac > minimum_open_log
fresh_open <- atac_change < -atac_cutoff & fresh_separated &
  fresh_absolutely_open
nrc_open <- atac_change > atac_cutoff & nrc_separated &
  nrc_absolutely_open

peak_table <- data.frame(
  Peak = peak_ids,
  Fresh_donor1_mean_log2_openness = fresh_d1_atac,
  Fresh_donor2_mean_log2_openness = fresh_d2_atac,
  NRC_donor1_mean_log2_openness = nrc_d1_atac,
  NRC_donor2_mean_log2_openness = nrc_d2_atac,
  Fresh_mean_log2_openness = fresh_atac,
  NRC_mean_log2_openness = nrc_atac,
  log2_NRC_over_Fresh = atac_change,
  Complete_donor_separation = fresh_separated | nrc_separated,
  Fresh_openness_above_2_in_both_donors = fresh_absolutely_open,
  NRC_openness_above_2_in_both_donors = nrc_absolutely_open,
  Direction = ifelse(
    fresh_open, "Fresh-opening",
    ifelse(nrc_open, "NRC-opening", "Other")
  )
)
write.table(
  peak_table, file.path(out_dir, "01_common_peak_state_changes.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

draw_peak_comparison <- function() {
  par(
    family = "sans", mar = c(4.0, 4.4, 0.7, 0.5),
    mgp = c(2.05, 0.52, 0), tcl = -0.22, las = 1,
    cex.axis = 0.68, cex.lab = 0.64, bty = "l"
  )
  set.seed(seed)
  other <- which(!fresh_open & !nrc_open)
  show_other <- sample(other, min(length(other), 45000L))
  show_fresh <- sample(which(fresh_open), min(sum(fresh_open), 12000L))
  show_nrc <- sample(which(nrc_open), min(sum(nrc_open), 12000L))
  lim <- range(
    quantile(c(fresh_atac, nrc_atac), c(0.002, 0.998), na.rm = TRUE),
    finite = TRUE
  )
  plot(
    fresh_atac[show_other], nrc_atac[show_other],
    pch = 16, cex = 0.16, col = adjustcolor(stable_colour, alpha.f = 0.18),
    xlim = lim, ylim = lim,
    xlab = expression("Fresh mean log"[2] * "(openness + 1)"),
    ylab = expression("NRC mean log"[2] * "(openness + 1)")
  )
  abline(a = 0, b = 1, col = "#777777", lty = 2, lwd = 0.7)
  points(
    fresh_atac[show_fresh], nrc_atac[show_fresh],
    pch = 16, cex = 0.20, col = adjustcolor(fresh_colour, alpha.f = 0.35)
  )
  points(
    fresh_atac[show_nrc], nrc_atac[show_nrc],
    pch = 16, cex = 0.20, col = adjustcolor(nrc_colour, alpha.f = 0.35)
  )
  legend(
    "topleft",
    legend = c(
      paste0("Fresh-opening: ", format_n(sum(fresh_open))),
      paste0("NRC-opening: ", format_n(sum(nrc_open))),
      paste0("Other: ", format_n(sum(!fresh_open & !nrc_open)))
    ),
    text.col = c(fresh_colour, nrc_colour, dark_grey),
    bty = "n", cex = 0.53
  )
}
save_plot("EDFig3a_Fresh_NRC_common_peak_comparison", 76, 64, draw_peak_comparison)

# -----------------------------------------------------------------------------
# Map both independently inferred networks to the same common peak universe.
# -----------------------------------------------------------------------------

message("Reading Fresh and NRC networks...")
fresh_pairs <- read_network_pairs(fresh_network_file, "Fresh")
nrc_pairs <- read_network_pairs(nrc_network_file, "NRC")
all_re <- union(fresh_pairs$RE, nrc_pairs$RE)
peak_index <- build_peak_index(peak_ids)
re_map <- map_re_to_peak(all_re, peak_index)
re_map <- re_map[
  is.finite(re_map$Peak_row) & re_map$RE_coverage >= minimum_re_coverage,
  , drop = FALSE
]

attach_mapping <- function(pairs) {
  pairs$Peak_row <- re_map$Peak_row[match(pairs$RE, re_map$RE)]
  pairs <- pairs[is.finite(pairs$Peak_row), , drop = FALSE]
  pairs$Peak <- peak_ids[pairs$Peak_row]
  pairs <- unique(pairs[, c("Gene", "Peak", "State")])
  pairs$ATAC_log2FC <- atac_change[pairs$Peak]
  pairs
}
fresh_mapped <- attach_mapping(fresh_pairs)
nrc_mapped <- attach_mapping(nrc_pairs)
mapped_pairs <- unique(rbind(fresh_mapped, nrc_mapped))
write.table(
  mapped_pairs, file.path(out_dir, "mapped_Fresh_NRC_network_RE_TG_pairs.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

fresh_network_peaks <- unique(fresh_mapped$Peak)
nrc_network_peaks <- unique(nrc_mapped$Peak)
network_peak_category <- data.frame(
  Peak = union(fresh_network_peaks, nrc_network_peaks),
  stringsAsFactors = FALSE
)
network_peak_category$Category <- ifelse(
  network_peak_category$Peak %in% fresh_network_peaks &
    network_peak_category$Peak %in% nrc_network_peaks,
  "Shared",
  ifelse(
    network_peak_category$Peak %in% fresh_network_peaks,
    "Fresh network only", "NRC network only"
  )
)
network_peak_category$ATAC_log2FC <- atac_change[network_peak_category$Peak]
network_peak_category$Category <- factor(
  network_peak_category$Category,
  levels = c("Fresh network only", "Shared", "NRC network only")
)
write.table(
  network_peak_category, file.path(out_dir, "02_network_peak_categories.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

network_category_summary <- do.call(rbind, lapply(
  levels(network_peak_category$Category), function(category) {
    x <- network_peak_category$ATAC_log2FC[
      network_peak_category$Category == category
    ]
    data.frame(
      Category = category,
      N_peaks = length(x),
      Median_state_difference = median(x, na.rm = TRUE),
      Q1 = unname(quantile(x, 0.25, na.rm = TRUE)),
      Q3 = unname(quantile(x, 0.75, na.rm = TRUE)),
      Fraction_NRC_higher = mean(x > 0, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }
))
write.table(
  network_category_summary,
  file.path(out_dir, "02_network_peak_category_summary.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

draw_network_peak_comparison <- function() {
  par(
    family = "sans", mar = c(6.8, 3.8, 0.7, 0.5),
    mgp = c(2.0, 0.52, 0), tcl = -0.22, las = 1,
    cex.axis = 0.62, cex.lab = 0.70, bty = "l"
  )
  values <- split(
    network_peak_category$ATAC_log2FC,
    network_peak_category$Category, drop = TRUE
  )
  ylim <- as.numeric(quantile(unlist(values), c(0.01, 0.99), na.rm = TRUE))
  labels <- paste0(
    c("Fresh network\nonly", "Shared", "NRC network\nonly"),
    "\n(n = ", format_n(lengths(values)), ")"
  )
  boxplot(
    values, outline = FALSE, ylim = ylim,
    col = c(
      adjustcolor(fresh_colour, alpha.f = 0.68),
      adjustcolor(shared_colour, alpha.f = 0.68),
      adjustcolor(nrc_colour, alpha.f = 0.68)
    ),
    border = c(fresh_colour, shared_colour, nrc_colour),
    names = rep("", 3), boxwex = 0.56, xaxt = "n",
    ylab = expression(Delta * " log"[2] * "(openness + 1), NRC - Fresh")
  )
  axis(1, at = 1:3, labels = FALSE, tick = FALSE)
  mtext(labels, side = 1, at = 1:3, line = 1.35, cex = 0.55)
  abline(h = 0, col = "#777777", lty = 2, lwd = 0.7)
}
# This category comparison is retained in the result tables for auditing but is
# not included in the focused Extended Data figure.

# -----------------------------------------------------------------------------
# Symmetric target-gene analysis using the union of Fresh and NRC network links.
# -----------------------------------------------------------------------------

message("Reading quantile-normalized FPKM...")
rna <- read.delim(rna_file, check.names = FALSE, quote = "", comment.char = "")
symbols <- toupper(trimws(rna$Symbol))
value_cols <- setdiff(names(rna), c("EnsemblID", "Symbol", "Alias"))
expr <- as.matrix(rna[, value_cols, drop = FALSE])
storage.mode(expr) <- "double"
valid <- !is.na(symbols) & nzchar(symbols)
expr <- expr[valid, , drop = FALSE]
symbols <- symbols[valid]
if (anyDuplicated(symbols)) {
  sums <- rowsum(expr, symbols, reorder = FALSE, na.rm = TRUE)
  counts <- as.numeric(table(factor(symbols, levels = rownames(sums))))
  expr <- sums / counts
} else {
  rownames(expr) <- symbols
}

mean_pattern <- function(patterns) {
  pattern <- paste(patterns, collapse = "|")
  cols <- grep(pattern, colnames(expr), value = TRUE)
  if (!length(cols)) stop("No RNA columns matched: ", pattern)
  rowMeans(expr[, cols, drop = FALSE], na.rm = TRUE)
}

fresh_d1_rna <- mean_pattern("^Fresh_1_")
fresh_d2_rna <- mean_pattern("^Fresh_2_")
nrc_d1_rna <- mean_pattern(c("^NRC_1_", "^Act_1_"))
nrc_d2_rna <- mean_pattern(c("^NRC_2_", "^Act_2_"))
fresh_rna <- (fresh_d1_rna + fresh_d2_rna) / 2
nrc_rna <- (nrc_d1_rna + nrc_d2_rna) / 2

peak_direction <- setNames(
  ifelse(fresh_open, "Fresh-opening",
         ifelse(nrc_open, "NRC-opening", "Other")),
  peak_ids
)
union_pairs <- unique(mapped_pairs[, c("Gene", "Peak")])
gene_rows <- split(seq_len(nrow(union_pairs)), union_pairs$Gene)
gene_table <- do.call(rbind, lapply(names(gene_rows), function(gene) {
  i <- gene_rows[[gene]]
  linked_peaks <- unique(union_pairs$Peak[i])
  changes <- atac_change[linked_peaks]
  directions <- peak_direction[linked_peaks]
  data.frame(
    Gene = gene,
    Linked_peaks = length(linked_peaks),
    NRC_opening_linked_peaks = sum(directions == "NRC-opening", na.rm = TRUE),
    Fresh_opening_linked_peaks = sum(directions == "Fresh-opening", na.rm = TRUE),
    Linked_peak_ATAC_median = median(changes, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}))
common_genes <- intersect(gene_table$Gene, rownames(expr))
gene_table <- gene_table[match(common_genes, gene_table$Gene), , drop = FALSE]
gene_table$Fresh_mean_FPKM <- fresh_rna[common_genes]
gene_table$NRC_mean_FPKM <- nrc_rna[common_genes]
gene_table$Baseline_FPKM <- (
  gene_table$Fresh_mean_FPKM + gene_table$NRC_mean_FPKM
) / 2
gene_table$RNA_log2FC <- log2(gene_table$NRC_mean_FPKM + 1) -
  log2(gene_table$Fresh_mean_FPKM + 1)
gene_table <- gene_table[
  is.finite(gene_table$Linked_peak_ATAC_median) &
    is.finite(gene_table$RNA_log2FC) & gene_table$Baseline_FPKM > 0.1,
  , drop = FALSE
]
gene_table$Target_class <- ifelse(
  gene_table$Linked_peaks >= minimum_linked_peaks &
    gene_table$NRC_opening_linked_peaks / gene_table$Linked_peaks >=
      nrc_opening_link_fraction,
  "NRC-opening-dominant targets", "Other network-linked targets"
)
tested_levels <- c("Other network-linked targets", "NRC-opening-dominant targets")
gene_table$Target_class <- factor(gene_table$Target_class, levels = tested_levels)

gene_table$Expression_bin <- pmin(
  10L,
  ceiling(rank(gene_table$Baseline_FPKM, ties.method = "average") /
            nrow(gene_table) * 10L)
)
gene_table$Degree_bin <- cut(
  gene_table$Linked_peaks,
  breaks = c(0, 1, 2, 4, 9, Inf),
  labels = c("1", "2", "3-4", "5-9", "10+"), right = TRUE
)
gene_table$Stratum <- interaction(
  gene_table$Expression_bin, gene_table$Degree_bin,
  drop = TRUE, lex.order = TRUE
)

matched_test <- function() {
  target_idx <- which(gene_table$Target_class == "NRC-opening-dominant targets")
  control_idx <- which(gene_table$Target_class == "Other network-linked targets")
  observed <- median(gene_table$RNA_log2FC[target_idx])
  if (length(target_idx) < 10L) {
    return(list(
      observed = observed, null = numeric(), p = NA_real_,
      n = length(target_idx), reason = "fewer than 10 target genes"
    ))
  }
  target_counts <- table(gene_table$Stratum[target_idx])
  controls <- split(control_idx, gene_table$Stratum[control_idx])
  sample_control <- function() {
    unlist(lapply(names(target_counts), function(s) {
      pool <- controls[[s]]
      n <- as.integer(target_counts[[s]])
      if (is.null(pool) || !length(pool)) {
        expression_bin <- strsplit(s, "\\.", fixed = FALSE)[[1]][1]
        pool <- control_idx[gene_table$Expression_bin[control_idx] ==
                              as.integer(expression_bin)]
      }
      if (!length(pool)) pool <- control_idx
      sample(pool, n, replace = length(pool) < n)
    }), use.names = FALSE)
  }
  set.seed(seed + 1L)
  null <- replicate(n_perm, {
    idx <- sample_control()
    median(gene_table$RNA_log2FC[idx])
  })
  p <- (1 + sum(null >= observed)) / (length(null) + 1)
  list(observed = observed, null = null, p = p, n = length(target_idx), reason = "tested")
}

nrc_test <- matched_test()
write.table(
  gene_table, file.path(out_dir, "EDFig3b_network_gene_level_data.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

target_class_summary <- do.call(rbind, lapply(tested_levels, function(category) {
  x <- gene_table$RNA_log2FC[gene_table$Target_class == category]
  data.frame(
    Target_class = category,
    N_genes = length(x),
    Median_RNA_log2FC = median(x, na.rm = TRUE),
    Q1 = unname(quantile(x, 0.25, na.rm = TRUE)),
    Q3 = unname(quantile(x, 0.75, na.rm = TRUE)),
    Fraction_NRC_higher = mean(x > 0, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}))
write.table(
  target_class_summary,
  file.path(out_dir, "EDFig3b_target_class_RNA_summary.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

format_p <- function(p) {
  if (!is.finite(p)) {
    "not tested (n < 10)"
  } else if (p < 1e-4) {
    "P < 1 x 10^-4"
  } else {
    paste0("P = ", formatC(p, digits = 2, format = "g"))
  }
}

draw_target_expression <- function() {
  par(
    family = "sans", mar = c(6.8, 3.8, 0.7, 0.5),
    mgp = c(2.0, 0.52, 0), tcl = -0.22, las = 1,
    cex.axis = 0.61, cex.lab = 0.70, bty = "l"
  )
  values <- split(gene_table$RNA_log2FC, gene_table$Target_class, drop = TRUE)
  ylim <- as.numeric(quantile(unlist(values), c(0.01, 0.99), na.rm = TRUE))
  labels <- paste0(
    c("Other network-linked\ntargets", "NRC-opening-dominant\ntargets"),
    "\n(n = ", format_n(lengths(values)), ")"
  )
  boxplot(
    values, outline = FALSE, ylim = ylim,
    col = c(
      adjustcolor(stable_colour, alpha.f = 0.68),
      adjustcolor(nrc_colour, alpha.f = 0.68)
    ),
    border = c("#777777", nrc_colour),
    names = rep("", 2), boxwex = 0.55, xaxt = "n",
    ylab = expression("RNA log"[2] * "(NRC / Fresh)")
  )
  axis(1, at = 1:2, labels = FALSE, tick = FALSE)
  mtext(labels, side = 1, at = 1:2, line = 1.35, cex = 0.54)
  abline(h = 0, col = "#777777", lty = 2, lwd = 0.7)
  legend(
    "topleft",
    legend = paste0("Expression- and degree-matched ", format_p(nrc_test$p)),
    text.col = nrc_colour, bty = "n", cex = 0.50
  )
}
save_plot("EDFig3b_NRC_opening_dominant_target_expression", 78, 62, draw_target_expression)

# -----------------------------------------------------------------------------
# Symmetrically selected state-associated target-gene heatmaps.
# -----------------------------------------------------------------------------

fresh_donor_matrix <- cbind(
  `Fresh 1` = fresh_d1_rna,
  `Fresh 2` = fresh_d2_rna
)
nrc_donor_matrix <- cbind(
  `NRC 1` = nrc_d1_rna,
  `NRC 2` = nrc_d2_rna
)
rownames(fresh_donor_matrix) <- rownames(expr)
rownames(nrc_donor_matrix) <- rownames(expr)

profiles <- cbind(
  fresh_donor_matrix,
  nrc_donor_matrix,
  `H1-hESC` = mean_pattern("^H1_hESC_"),
  `iPSC` = mean_pattern("^iPSC_"),
  `MSC` = mean_pattern("^mesenchymalStemCell_ENC"),
  `MSC of BM` = mean_pattern("^mesenchymalStemCellOfTheBoneMarrow_"),
  `HMPC` = mean_pattern("^hematopoieticMultipotentProgenitorCell_"),
  `Monocyte` = mean_pattern("^CD14_positiveMonocyte_"),
  `PBMC` = mean_pattern("^peripheralBloodMononuclearCell_")
)

gene_table$Fresh_donor_min <- apply(
  fresh_donor_matrix[gene_table$Gene, , drop = FALSE], 1, min
)
gene_table$Fresh_donor_max <- apply(
  fresh_donor_matrix[gene_table$Gene, , drop = FALSE], 1, max
)
gene_table$NRC_donor_min <- apply(
  nrc_donor_matrix[gene_table$Gene, , drop = FALSE], 1, min
)
gene_table$NRC_donor_max <- apply(
  nrc_donor_matrix[gene_table$Gene, , drop = FALSE], 1, max
)

nrc_selected <- gene_table[
  gene_table$Target_class == "NRC-opening-dominant targets" &
    gene_table$RNA_log2FC > rna_cutoff &
    gene_table$NRC_mean_FPKM >= mean_fpkm_cutoff &
    gene_table$NRC_donor_min > gene_table$Fresh_donor_max,
  , drop = FALSE
]
nrc_selected$Rank_score <- nrc_selected$RNA_log2FC *
  log2(nrc_selected$NRC_mean_FPKM + 1) *
  abs(nrc_selected$Linked_peak_ATAC_median)
nrc_selected <- nrc_selected[order(nrc_selected$Rank_score, decreasing = TRUE), ]

write.table(
  nrc_selected, file.path(out_dir, "EDFig3c_NRC_opening_dominant_target_genes.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

heat_col <- colorRampPalette(c("#3B82B8", "#F7F7F7", "#D6604D"))(101)
column_groups <- c(
  "Fresh", "Fresh", "NRC", "NRC", "Stem", "Stem", "Stem", "Stem",
  "Blood", "Blood", "Blood"
)
annotation_colours <- c(
  Fresh = fresh_colour, NRC = nrc_colour,
  Stem = "#65A765", Blood = "#B07AA1"
)

draw_heatmap <- function(selected, state, max_genes = 60L) {
  if (!nrow(selected)) stop("No selected genes for ", state)
  shown <- head(selected, max_genes)
  m <- log2(profiles[shown$Gene, , drop = FALSE] + 1)
  z <- t(scale(t(m)))
  z[!is.finite(z)] <- 0
  z[z > 2.5] <- 2.5
  z[z < -2.5] <- -2.5
  row_order <- if (nrow(z) > 1L) hclust(dist(z), method = "complete")$order else 1L
  z <- z[row_order, , drop = FALSE]
  gene_labels <- rownames(z)
  n_gene <- nrow(z)

  function() {
    layout(matrix(c(1, 2), nrow = 1), widths = c(14, 3.0))
    par(
      family = "sans", mar = c(7.0, 1.0, 2.8, 6.2),
      mgp = c(1.8, 0.45, 0), tcl = -0.2, xpd = NA
    )
    image(
      x = seq_len(ncol(z)), y = seq_len(n_gene),
      z = t(z[n_gene:1, , drop = FALSE]),
      col = heat_col, zlim = c(-2.5, 2.5), axes = FALSE,
      xlab = "", ylab = "", useRaster = TRUE
    )
    abline(v = seq(1.5, ncol(z) - 0.5, by = 1), col = "white", lwd = 0.35)
    if (n_gene > 1L) {
      abline(h = seq(1.5, n_gene - 0.5, by = 1), col = "white", lwd = 0.25)
    }
    box(col = "#777777", lwd = 0.45)
    axis(
      1, at = seq_len(ncol(z)), labels = colnames(z),
      las = 2, tick = FALSE, line = -0.2, cex.axis = 0.68
    )
    italic_labels <- as.expression(lapply(
      rev(gene_labels), function(g) bquote(italic(.(g)))
    ))
    axis(
      4, at = seq_len(n_gene), labels = italic_labels,
      las = 1, tick = FALSE, line = -0.2,
      cex.axis = if (n_gene > 45) 0.46 else if (n_gene > 30) 0.53 else 0.62
    )
    usr <- par("usr")
    y0 <- usr[4] + 0.15
    rect(
      seq_len(ncol(z)) - 0.48, y0,
      seq_len(ncol(z)) + 0.48, y0 + 0.36,
      col = annotation_colours[column_groups], border = NA
    )
    mtext(state, side = 3, line = 1.25, adj = 0.5, cex = 0.86)
    par(
      family = "sans", mar = c(7.0, 0.2, 3.7, 2.3),
      mgp = c(1.4, 0.3, 0), tcl = -0.18, xpd = FALSE
    )
    key_values <- seq(-2.5, 2.5, length.out = 101)
    image(
      x = 1, y = key_values, z = matrix(key_values, nrow = 1),
      col = heat_col, zlim = c(-2.5, 2.5), axes = FALSE,
      xlab = "", ylab = "", useRaster = TRUE
    )
    axis(4, at = c(-2.5, 0, 2.5), labels = c("-2.5", "0", "2.5"),
         las = 1, cex.axis = 0.58)
    mtext("Row\nz-score", side = 3, line = 0.2, cex = 0.58)
    box(col = "#777777", lwd = 0.4)
    layout(1)
  }
}

nrc_n_show <- min(60L, nrow(nrc_selected))
nrc_height <- max(80, 34 + 3.1 * nrc_n_show)

save_plot(
  "EDFig3c_NRC_opening_dominant_targets_heatmap", 140, nrc_height,
  draw_heatmap(nrc_selected, "NRC-opening-dominant target genes", max_genes = 60L)
)

summary_lines <- c(
  "Fresh-NRC Extended Data Fig. 3 core analysis using raw Duren-style openness scores",
  "Fresh and NRC donors were independent; technical libraries were averaged within donor.",
  "No quantile normalization was applied; input openness scores were already normalized to local background.",
  paste0(
    "Opening calls required openness > ", minimum_open_openness,
    " in both donors of the higher state, |state difference| > 1, and complete donor separation."
  ),
  paste0("Retained common peaks: ", format_n(length(peak_ids))),
  paste0("Fresh-opening common peaks: ", format_n(sum(fresh_open))),
  paste0("NRC-opening common peaks: ", format_n(sum(nrc_open))),
  paste0("Median common-peak log2(NRC/Fresh): ",
         formatC(median(atac_change), digits = 4, format = "f")),
  paste0("Mapped Fresh network peaks: ", format_n(length(fresh_network_peaks))),
  paste0("Mapped NRC network peaks: ", format_n(length(nrc_network_peaks))),
  paste0("Mapped shared network peaks: ",
         format_n(length(intersect(fresh_network_peaks, nrc_network_peaks)))),
  paste0(
    "Network category median state differences: ",
    paste(
      paste0(network_category_summary$Category, " = ",
             formatC(network_category_summary$Median_state_difference,
                     digits = 4, format = "f")),
      collapse = "; "
    )
  ),
  paste0(
    "NRC-opening-dominant target genes (at least ", minimum_linked_peaks,
    " linked peaks; fraction NRC-opening >= ", nrc_opening_link_fraction,
    "): ", nrc_test$n,
         "; matched empirical result = ", format_p(nrc_test$p)),
  paste0("Other network-linked target genes: ",
         sum(gene_table$Target_class == "Other network-linked targets")),
  paste0("NRC-opening-dominant, transcriptionally activated heatmap genes: ",
         nrow(nrc_selected)),
  "The regulatory networks integrate chromatin and RNA information; results show integrative concordance, not causality or independent validation."
)
writeLines(summary_lines, file.path(out_dir, "analysis_summary.txt"))
message("Finished: ", out_dir)
