# Additional ER-expression figure helpers
#
# Clean public helpers for expression and GSEA panels used across ER-response
# analyses: patient CDK4/6 inhibitor barcode plots, ER target heatmaps,
# CCND1/single-gene expression panels, and multi-cell-line expression summaries.
# Runtime input paths are optional so panels can be generated independently.

suppressPackageStartupMessages({
  library(ComplexHeatmap)
  library(circlize)
  library(clusterProfiler)
  library(enrichplot)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
  library(qs)
  library(tidyverse)
})

project_dir <- normalizePath(Sys.getenv("PROJECT_DIR", "."), winslash = "/", mustWork = FALSE)
output_dir <- file.path(Sys.getenv("FIGURES_OUTPUT_DIR", file.path(project_dir, "code_availability", "Figures", "output")), "05_ER_expression")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

read_object <- function(file) {
  if (!file.exists(file)) {
    warning("Missing input object: ", file, call. = FALSE)
    return(NULL)
  }
  if (grepl("\\.qs$", file, ignore.case = TRUE)) qs::qread(file) else readRDS(file)
}

# Mark gene-set members along a precomputed ranked expression vector.
plot_gsea_barcode_from_rank <- function(rank_df, genes, comparison_label, title = comparison_label) {
  stopifnot(all(c("gene", "rank_metric") %in% colnames(rank_df)))
  barcode_data <- rank_df %>%
    arrange(desc(rank_metric)) %>%
    mutate(rank = row_number(), hit = gene %in% genes) %>%
    filter(hit)

  ggplot(barcode_data, aes(x = rank, y = 1)) +
    geom_linerange(aes(ymin = 0, ymax = 1), color = "black", linewidth = 0.25) +
    labs(title = title, x = "Ranked genes", y = NULL) +
    theme_classic() +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
}

# Subset, optionally row-scale and draw an expression heatmap.
plot_expression_heatmap <- function(expr_mat, genes, column_annotation = NULL, z_score = TRUE, title = "") {
  genes <- intersect(genes, rownames(expr_mat))
  mat <- expr_mat[genes, , drop = FALSE]
  mat <- mat[rowSums(is.na(mat)) < ncol(mat), , drop = FALSE]
  if (z_score) mat <- t(scale(t(mat)))
  mat[is.na(mat)] <- 0

  top_anno <- NULL
  if (!is.null(column_annotation)) {
    top_anno <- HeatmapAnnotation(df = column_annotation)
  }

  Heatmap(
    mat,
    name = if (z_score) "z-score" else "expression",
    top_annotation = top_anno,
    col = circlize::colorRamp2(c(-2, 0, 2), c("dodgerblue3", "white", "orange")),
    cluster_rows = TRUE,
    cluster_columns = FALSE,
    show_row_names = TRUE,
    row_names_gp = grid::gpar(fontsize = 6),
    column_title = title
  )
}

plot_single_gene_expression <- function(expr_df, gene, group_col = "group", value_col = "expression") {
  expr_df %>%
    filter(.data$gene == !!gene) %>%
    ggplot(aes(x = .data[[group_col]], y = .data[[value_col]], fill = .data[[group_col]])) +
    geom_boxplot(outlier.shape = NA, alpha = 0.75) +
    geom_jitter(width = 0.15, size = 1.5, alpha = 0.8) +
    labs(x = NULL, y = gene) +
    theme_classic() +
    theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))
}

# Optional runtime inputs. These are intentionally generic so local/private
# object names are not embedded in the public script.
rank_file <- Sys.getenv("ER_RANK_TABLE", unset = "")
geneset_file <- Sys.getenv("ER_GENESET_FILE", unset = "")
expression_matrix_file <- Sys.getenv("ER_EXPRESSION_MATRIX", unset = "")
single_gene_file <- Sys.getenv("ER_SINGLE_GENE_TABLE", unset = "")

if (nzchar(rank_file) && nzchar(geneset_file)) {
  rank_df <- readr::read_csv(rank_file, show_col_types = FALSE)
  genes <- readr::read_lines(geneset_file)
  p <- plot_gsea_barcode_from_rank(rank_df, genes, comparison_label = "ER response")
  ggsave(file.path(output_dir, "ER_response_barcodeplot.pdf"), plot = p, width = 5, height = 2.2)
}

if (nzchar(expression_matrix_file) && nzchar(geneset_file)) {
  expr_mat <- as.matrix(read.csv(expression_matrix_file, row.names = 1, check.names = FALSE))
  genes <- readr::read_lines(geneset_file)
  pdf(file.path(output_dir, "ER_target_expression_heatmap.pdf"), width = 5, height = 6)
  draw(plot_expression_heatmap(expr_mat, genes, title = "ER target genes"))
  dev.off()
}

if (nzchar(single_gene_file)) {
  expr_df <- readr::read_csv(single_gene_file, show_col_types = FALSE)
  if (all(c("gene", "group", "expression") %in% names(expr_df))) {
    p <- plot_single_gene_expression(expr_df, gene = Sys.getenv("ER_SINGLE_GENE", unset = "CCND1"))
    ggsave(file.path(output_dir, "single_gene_expression.pdf"), plot = p, width = 4, height = 3.5)
  }
}

message("Saved ER-expression figure outputs to: ", output_dir)
capture.output(sessionInfo(), file = file.path(output_dir, "sessionInfo.txt"))
