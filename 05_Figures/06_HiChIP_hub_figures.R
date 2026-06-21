# HiChIP and enhancer-hub figure helpers
#
# Clean public helpers for RB non-promoter/HiChIP panels, including
# hypergeometric/ORA dot plots, hub count plots, hub gene-expression heatmaps,
# and selected loop/hub gene visualisations.
#
# Method note: hypergeometric/ORA result objects are expected to contain
# `GeneRatio`, `FoldEnrichment`, and BH-adjusted `p.adjust`. Figure panels use
# `p.adjust < 0.05` as the display/significance cutoff.
# Inputs are saved hub summaries/ORA objects; this script does not rebuild hubs.

suppressPackageStartupMessages({
  library(ComplexHeatmap)
  library(circlize)
  library(qs)
  library(tidyverse)
})

project_dir <- normalizePath(Sys.getenv("PROJECT_DIR", "."), winslash = "/", mustWork = FALSE)
output_dir <- file.path(Sys.getenv("FIGURES_OUTPUT_DIR", file.path(project_dir, "code_availability", "Figures", "output")), "06_HiChIP_hubs")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

read_object <- function(file) {
  if (!file.exists(file)) {
    warning("Missing input object: ", file, call. = FALSE)
    return(NULL)
  }
  if (grepl("\\.qs$", file, ignore.case = TRUE)) qs::qread(file) else readRDS(file)
}

# Convert clusterProfiler ratio strings such as 12/140 to numeric values.
gene_ratio_number <- function(x) {
  vapply(strsplit(as.character(x), "/", fixed = TRUE), function(parts) {
    as.numeric(parts[[1]]) / as.numeric(parts[[2]])
  }, numeric(1))
}

# Combine selected hub-gene ORA comparisons into one dot plot.
plot_hypergeo_dotplot <- function(hypergeo_list, selected_groups = names(hypergeo_list), fdr_cutoff = 0.05, top_n = 8) {
  hypergeo_df <- hypergeo_list[selected_groups] %>%
    lapply(as_tibble)

  top_genesets <- hypergeo_df %>%
    lapply(function(x) x %>% arrange(p.adjust) %>% pull(ID) %>% head(top_n)) %>%
    unlist() %>%
    unique()

  plot_df <- bind_rows(lapply(names(hypergeo_df), function(nm) {
    hypergeo_df[[nm]] %>%
      filter(ID %in% top_genesets) %>%
      mutate(group = nm)
  }))

  if ("GeneRatio" %in% names(plot_df)) {
    plot_df$GeneRatio_number <- gene_ratio_number(plot_df$GeneRatio)
  } else {
    plot_df$GeneRatio_number <- NA_real_
  }

  plot_df <- plot_df %>%
    mutate(
      ID = factor(ID, levels = rev(top_genesets)),
      group = factor(group, levels = selected_groups),
      p.adjust = if_else(p.adjust > fdr_cutoff, NA_real_, p.adjust)
    )

  ggplot(plot_df, aes(x = group, y = ID, color = FoldEnrichment, size = GeneRatio_number)) +
    geom_point(alpha = 0.85) +
    scale_color_gradientn(colors = c("orange", "red", "darkred"), na.value = "gray85") +
    labs(x = NULL, y = NULL, color = "Fold enrichment", size = "Gene ratio") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

plot_hub_counts <- function(count_df, group_col = "group", count_col = "count") {
  count_df %>%
    ggplot(aes(x = .data[[group_col]], y = .data[[count_col]], fill = .data[[group_col]])) +
    geom_col(width = 0.7) +
    labs(x = NULL, y = "Count") +
    theme_classic() +
    theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))
}

plot_hub_gene_heatmap <- function(expr_mat, genes, title = "Hub genes") {
  genes <- intersect(genes, rownames(expr_mat))
  mat <- expr_mat[genes, , drop = FALSE]
  mat <- mat[rowSums(mat != 0, na.rm = TRUE) > 0, , drop = FALSE]
  mat <- t(scale(t(mat)))
  mat[is.na(mat)] <- 0

  Heatmap(
    mat,
    name = "z-score",
    col = circlize::colorRamp2(c(-2, 0, 2), c("dodgerblue3", "white", "orange")),
    cluster_rows = TRUE,
    cluster_columns = FALSE,
    show_row_names = TRUE,
    row_names_gp = grid::gpar(fontsize = 6),
    column_title = title
  )
}

hypergeo_file <- Sys.getenv("HICHIP_HUB_HYPERGEO_QS", unset = "")
hub_counts_file <- Sys.getenv("HICHIP_HUB_COUNTS_CSV", unset = "")
hub_expr_file <- Sys.getenv("HICHIP_HUB_EXPRESSION_MATRIX", unset = "")
hub_genes_file <- Sys.getenv("HICHIP_HUB_GENES", unset = "")

if (nzchar(hypergeo_file)) {
  hypergeo_output <- read_object(hypergeo_file)
  if (!is.null(hypergeo_output)) {
    p <- plot_hypergeo_dotplot(hypergeo_output)
    ggsave(file.path(output_dir, "HiChIP_hub_hypergeometric_dotplot.pdf"), plot = p, width = 6, height = 4)
  }
}

if (nzchar(hub_counts_file)) {
  count_df <- readr::read_csv(hub_counts_file, show_col_types = FALSE)
  p <- plot_hub_counts(count_df)
  ggsave(file.path(output_dir, "HiChIP_hub_counts.pdf"), plot = p, width = 4, height = 3)
}

if (nzchar(hub_expr_file) && nzchar(hub_genes_file)) {
  expr_mat <- as.matrix(read.csv(hub_expr_file, row.names = 1, check.names = FALSE))
  genes <- readr::read_lines(hub_genes_file)
  pdf(file.path(output_dir, "HiChIP_hub_gene_expression_heatmap.pdf"), width = 5, height = 6)
  draw(plot_hub_gene_heatmap(expr_mat, genes))
  dev.off()
}

message("Saved HiChIP/hub figure outputs to: ", output_dir)
capture.output(sessionInfo(), file = file.path(output_dir, "sessionInfo.txt"))
