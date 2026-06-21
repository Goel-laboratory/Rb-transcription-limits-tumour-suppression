# KDM5A/KDM5B expression and mechanism figure helpers
#
# Clean public helpers for KDM5A/KDM5B mechanism panels beyond the mass-spec
# volcano: RNA-seq heatmaps, selected gene-expression changes, and overlap
# summaries.
# Optional runtime inputs allow each panel family to run independently.

suppressPackageStartupMessages({
  library(ComplexHeatmap)
  library(circlize)
  library(qs)
  library(tidyverse)
})

project_dir <- normalizePath(Sys.getenv("PROJECT_DIR", "."), winslash = "/", mustWork = FALSE)
output_dir <- file.path(Sys.getenv("FIGURES_OUTPUT_DIR", file.path(project_dir, "code_availability", "Figures", "output")), "07_KDM5A_KDM5B_expression")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

read_object <- function(file) {
  if (!file.exists(file)) {
    warning("Missing input object: ", file, call. = FALSE)
    return(NULL)
  }
  if (grepl("\\.qs$", file, ignore.case = TRUE)) qs::qread(file) else readRDS(file)
}

# Express selected genes relative to configured baseline samples.
make_relative_heatmap <- function(expr_mat, genes, baseline_cols, display_cols = colnames(expr_mat), z_score = TRUE, title = "") {
  genes <- intersect(genes, rownames(expr_mat))
  mat <- expr_mat[genes, display_cols, drop = FALSE]
  baseline <- rowMeans(expr_mat[genes, baseline_cols, drop = FALSE], na.rm = TRUE)
  mat <- sweep(mat, 1, baseline, FUN = "-")
  mat <- mat[rowSums(mat != 0, na.rm = TRUE) > 0, , drop = FALSE]
  if (z_score) mat <- t(scale(t(mat)))
  mat[is.na(mat)] <- 0

  Heatmap(
    mat,
    name = if (z_score) "z-score" else "relative expression",
    col = circlize::colorRamp2(c(-2, 0, 2), c("dodgerblue3", "white", "orange")),
    cluster_rows = TRUE,
    cluster_columns = FALSE,
    row_names_gp = grid::gpar(fontsize = 6),
    column_title = title
  )
}

# Extract fold changes and significance for a defined mechanism gene list.
summarise_gene_changes <- function(dge_df, genes, symbol_col = "SYMBOL", logfc_col = "log2FoldChange", padj_col = "padj", padj_cutoff = 0.05) {
  dge_df %>%
    as_tibble() %>%
    filter(.data[[symbol_col]] %in% genes) %>%
    transmute(
      gene = .data[[symbol_col]],
      log2FC = .data[[logfc_col]],
      padj = .data[[padj_col]],
      significant = !is.na(padj) & padj < padj_cutoff
    ) %>%
    arrange(desc(log2FC))
}

plot_gene_change_bars <- function(gene_change_df) {
  gene_change_df %>%
    mutate(gene = factor(gene, levels = gene)) %>%
    ggplot(aes(x = gene, y = log2FC, fill = significant)) +
    geom_col(width = 0.75) +
    geom_hline(yintercept = 0, linewidth = 0.3) +
    coord_flip() +
    scale_fill_manual(values = c("TRUE" = "orange", "FALSE" = "gray70")) +
    labs(x = NULL, y = "log2 fold change") +
    theme_classic()
}

expr_mat_file <- Sys.getenv("KDM5_EXPRESSION_MATRIX", unset = "")
genes_file <- Sys.getenv("KDM5_GENE_LIST", unset = "")
baseline_cols <- strsplit(Sys.getenv("KDM5_BASELINE_COLUMNS", unset = ""), ",", fixed = TRUE)[[1]]
dge_file <- Sys.getenv("KDM5_DGE_TABLE", unset = "")

if (nzchar(expr_mat_file) && nzchar(genes_file) && length(baseline_cols) > 0 && nzchar(baseline_cols[[1]])) {
  expr_mat <- as.matrix(read.csv(expr_mat_file, row.names = 1, check.names = FALSE))
  genes <- readr::read_lines(genes_file)
  pdf(file.path(output_dir, "KDM5_relative_expression_heatmap.pdf"), width = 5, height = 6)
  draw(make_relative_heatmap(expr_mat, genes, baseline_cols = baseline_cols, title = "KDM5 mechanism genes"))
  dev.off()
}

if (nzchar(dge_file) && nzchar(genes_file)) {
  dge_df <- read_object(dge_file)
  genes <- readr::read_lines(genes_file)
  gene_changes <- summarise_gene_changes(as.data.frame(dge_df), genes)
  write.csv(gene_changes, file.path(output_dir, "KDM5_selected_gene_changes.csv"), row.names = FALSE)
  p <- plot_gene_change_bars(gene_changes)
  ggsave(file.path(output_dir, "KDM5_selected_gene_changes.pdf"), plot = p, width = 4, height = 5)
}

message("Saved KDM5A/KDM5B expression figure outputs to: ", output_dir)
capture.output(sessionInfo(), file = file.path(output_dir, "sessionInfo.txt"))
