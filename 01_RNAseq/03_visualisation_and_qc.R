#!/usr/bin/env Rscript

# Generate RNA-seq PCA, sample-distance, MA, volcano and heatmap figures.
# Inputs are the fitted DESeq2 object, result list and sample metadata.

suppressPackageStartupMessages({
  library(DESeq2)
  library(dplyr)
  library(readr)
  library(tibble)
  library(ggplot2)
  library(ggrepel)
  library(pheatmap)
  library(RColorBrewer)
  library(ComplexHeatmap)
  library(circlize)
  library(msigdbr)
})

base_output <- Sys.getenv("RNASEQ_OUTPUT_DIR", "code_availability/RNAseq/output")
deseq_dir <- Sys.getenv("RNASEQ_DESEQ2_DIR", file.path(base_output, "deseq2"))
dds_file <- Sys.getenv("DDS_FILE", file.path(deseq_dir, "deseq2_dds.rds"))
results_rds <- Sys.getenv("DESEQ_RESULTS_RDS", file.path(deseq_dir, "deseq2_results_list.rds"))
metadata_file <- Sys.getenv("METADATA_FILE", "")
outdir <- Sys.getenv("RNASEQ_PLOT_DIR", file.path(base_output, "plots"))
if (!file.exists(dds_file) || !file.exists(results_rds)) stop("Missing DESeq2 outputs.")
if (!nzchar(metadata_file) || !file.exists(metadata_file)) stop("Set METADATA_FILE.")
sample_id_col <- "sample"
color_by <- "condition"
shape_by <- NA_character_
top_variable_genes <- 1000
padj_threshold <- 0.05
label_top_n <- 20

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

dds <- readRDS(dds_file)
res_list <- readRDS(results_rds)
metadata <- readr::read_csv(metadata_file, show_col_types = FALSE) %>% as.data.frame()
rownames(metadata) <- metadata[[sample_id_col]]

# Reuse DESeq2 transformations so all plots share the same normalization.
normalized_counts <- counts(dds, normalized = TRUE)
log2_normalized_counts <- log2(normalized_counts + 1)
rld <- rlog(dds, blind = FALSE)
rld_mat <- assay(rld)

write.csv(normalized_counts, file.path(outdir, "normalized_counts_from_dds.csv"))
write.csv(rld_mat, file.path(outdir, "rlog_counts_from_dds.csv"))

###############################################################################
# PCA on rlog-transformed counts
###############################################################################
row_vars <- matrixStats::rowVars(rld_mat)
select <- order(row_vars, decreasing = TRUE)[seq_len(min(top_variable_genes, length(row_vars)))]
pca <- prcomp(t(rld_mat[select, , drop = FALSE]), center = TRUE, scale. = TRUE)
metadata_for_join <- metadata
metadata_for_join[[sample_id_col]] <- rownames(metadata_for_join)
pca_df <- as.data.frame(pca$x) %>%
  tibble::rownames_to_column(sample_id_col) %>%
  dplyr::left_join(metadata_for_join, by = sample_id_col)
percent_var <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)

pca_mapping <- if (!is.na(color_by) && color_by %in% colnames(pca_df)) {
  aes(x = PC1, y = PC2, label = .data[[sample_id_col]], colour = .data[[color_by]])
} else {
  aes(x = PC1, y = PC2, label = .data[[sample_id_col]])
}

p <- ggplot(pca_df, pca_mapping) +
  geom_point(size = 3) +
  geom_text_repel(size = 3, max.overlaps = Inf) +
  labs(x = paste0("PC1 (", percent_var[1], "%)"),
       y = paste0("PC2 (", percent_var[2], "%)")) +
  theme_bw(base_size = 11)
ggsave(file.path(outdir, "PCA_rlog_PC1_PC2.pdf"), p, width = 6, height = 5)

###############################################################################
# Sample-distance heatmap
###############################################################################
sample_dist <- dist(t(rld_mat))
sample_dist_mat <- as.matrix(sample_dist)
annotation_col <- metadata[colnames(sample_dist_mat), , drop = FALSE]
pdf(file.path(outdir, "sample_distance_heatmap.pdf"), width = 7, height = 6)
pheatmap::pheatmap(sample_dist_mat,
                  clustering_distance_rows = sample_dist,
                  clustering_distance_cols = sample_dist,
                  annotation_col = annotation_col,
                  main = "Sample distances from rlog counts")
dev.off()

###############################################################################
# Count-distribution plot
###############################################################################
plot_df <- as.data.frame(log2_normalized_counts) %>%
  tibble::rownames_to_column("gene") %>%
  tidyr::pivot_longer(-gene, names_to = "sample", values_to = "log2_normalized_count")
p <- ggplot(plot_df, aes(x = sample, y = log2_normalized_count)) +
  geom_boxplot(outlier.size = 0.2) +
  coord_flip() +
  theme_bw(base_size = 10) +
  labs(x = NULL, y = "log2(normalized counts + 1)")
ggsave(file.path(outdir, "count_distribution_boxplots.pdf"), p, width = 7, height = 6)

###############################################################################
# MA and volcano plots
###############################################################################
for (nm in names(res_list)) {
  res <- res_list[[nm]]
  logfc_col <- if ("logFC_MMSE" %in% colnames(res)) "logFC_MMSE" else if ("logFC_MLE" %in% colnames(res)) "logFC_MLE" else "log2FoldChange"
  gene_col <- if ("SYMBOL" %in% colnames(res)) "SYMBOL" else if ("genesymb" %in% colnames(res)) "genesymb" else "gene_id"
  res <- res %>%
    mutate(significant = !is.na(padj) & padj < padj_threshold,
           neg_log10_padj = -log10(padj))

  p_ma <- ggplot(res, aes(x = baseMean, y = .data[[logfc_col]])) +
    geom_point(data = dplyr::filter(res, !significant), color = "dodgerblue2", size = 0.5, alpha = 0.08) +
    geom_point(data = dplyr::filter(res, significant), color = "darkorange3", size = 0.6, alpha = 0.4) +
    scale_x_log10() +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    labs(x = "Mean expression (baseMean, log10 scale)",
         y = "log2 fold change",
         title = nm) +
    theme_bw(base_size = 11)
  ggsave(file.path(outdir, paste0(nm, "_MA_plot.pdf")), p_ma, width = 6, height = 5)

  label_df <- res %>%
    dplyr::filter(significant, !is.na(.data[[gene_col]])) %>%
    dplyr::arrange(padj) %>%
    dplyr::slice_head(n = label_top_n)
  p_volcano <- ggplot(res, aes(x = .data[[logfc_col]], y = neg_log10_padj)) +
    geom_point(data = dplyr::filter(res, !significant), color = "dodgerblue2", size = 0.5, alpha = 0.08) +
    geom_point(data = dplyr::filter(res, significant), color = "darkorange3", size = 0.7, alpha = 0.45) +
    geom_hline(yintercept = -log10(padj_threshold), linetype = "dashed", color = "grey50") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    geom_text_repel(data = label_df, aes(label = .data[[gene_col]]), size = 2.5, max.overlaps = Inf) +
    labs(x = "log2 fold change", y = "-log10(adjusted P value)", title = nm) +
    theme_bw(base_size = 11)
  ggsave(file.path(outdir, paste0(nm, "_volcano_plot.pdf")), p_volcano, width = 6, height = 5)
}

###############################################################################
# Heatmap of significant DE genes from the first comparison
###############################################################################
first_name <- names(res_list)[1]
first_res <- res_list[[first_name]]
gene_col <- if ("SYMBOL" %in% colnames(first_res)) "SYMBOL" else if ("genesymb" %in% colnames(first_res)) "genesymb" else "gene_id"
sig_genes <- first_res %>%
  dplyr::filter(!is.na(padj), padj < padj_threshold) %>%
  dplyr::arrange(padj) %>%
  dplyr::slice_head(n = 50) %>%
  dplyr::pull(.data[[gene_col]]) %>%
  unique()
id_lookup <- first_res %>%
  transmute(gene_id, display_gene = .data[[gene_col]]) %>%
  filter(!is.na(display_gene), !duplicated(display_gene))
sig_ids <- id_lookup$gene_id[match(sig_genes, id_lookup$display_gene)]
sig_genes <- intersect(na.omit(sig_ids), rownames(rld_mat))

if (length(sig_genes) >= 2) {
  heatmap_mat <- rld_mat[sig_genes, , drop = FALSE]
  heatmap_mat <- heatmap_mat[rowSums(is.na(heatmap_mat)) < ncol(heatmap_mat), , drop = FALSE]
  scaled_mat <- t(scale(t(heatmap_mat)))
  pdf(file.path(outdir, paste0(first_name, "_top_DEG_heatmap.pdf")), width = 7, height = 8)
  pheatmap::pheatmap(scaled_mat,
                    annotation_col = annotation_col[colnames(scaled_mat), , drop = FALSE],
                    cluster_cols = TRUE,
                    cluster_rows = TRUE,
                    show_rownames = TRUE,
                    main = paste(first_name, "top DE genes"))
  dev.off()
}

###############################################################################
# Optional Hallmark gene-set heatmaps
###############################################################################
h_gene_sets <- msigdbr::msigdbr(species = "Homo sapiens", category = "H") %>%
  dplyr::select(gs_name, gene_symbol)
hallmarks_to_plot <- c("HALLMARK_E2F_TARGETS", "HALLMARK_ESTROGEN_RESPONSE_EARLY")

for (gs in hallmarks_to_plot) {
  genes <- h_gene_sets %>% dplyr::filter(gs_name == gs) %>% dplyr::pull(gene_symbol)
  symbol_lookup <- first_res %>%
    transmute(gene_id, symbol = .data[[gene_col]]) %>%
    filter(!is.na(symbol), !duplicated(symbol))
  genes <- symbol_lookup$gene_id[match(genes, symbol_lookup$symbol)]
  genes <- intersect(na.omit(genes), rownames(log2_normalized_counts))
  if (length(genes) < 2) next
  mat <- log2_normalized_counts[genes, , drop = FALSE]
  mat <- mat[rowSums(mat == 0) < ncol(mat), , drop = FALSE]
  scaled_mat <- t(scale(t(mat)))
  pdf(file.path(outdir, paste0(gs, "_heatmap.pdf")), width = 7, height = 8)
  pheatmap::pheatmap(scaled_mat,
                    annotation_col = annotation_col[colnames(scaled_mat), , drop = FALSE],
                    show_rownames = FALSE,
                    main = gs)
  dev.off()
}

session_info <- utils::capture.output(sessionInfo())
writeLines(session_info, file.path(outdir, "sessionInfo_visualisation_qc.txt"))
