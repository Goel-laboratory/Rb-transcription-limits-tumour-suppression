# ER-expression paper figures from RNA-seq/GSEA
#
# Cleaned from:
#   analysis_paper/03_ER_exprs/01_RNASeq_GSEA_BubblePlot/MCF7M_Par_sgRB1_ZR_GSEA_bubbleplot.Rmd
#
# This script reruns the lightweight GSEA and figure assembly from saved DESeq2
# result tables. It does not rerun alignment, featureCounts, or DESeq2.
# Outputs are Hallmark GSEA tables and ER-response bubble plots.

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(msigdbr)
  library(tidyverse)
})

project_dir <- normalizePath(Sys.getenv("PROJECT_DIR", "."), winslash = "/", mustWork = FALSE)
output_dir <- file.path(Sys.getenv("FIGURES_OUTPUT_DIR", file.path(project_dir, "code_availability", "Figures", "output")), "03_ER_exprs")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

path <- function(...) file.path(project_dir, ...)

load_rds_or_warn <- function(file) {
  if (!file.exists(file)) {
    warning("Missing input object: ", file, call. = FALSE)
    return(NULL)
  }
  readRDS(file)
}

# Prefer the DESeq2 Wald statistic, with documented project fallbacks.
get_ranks <- function(res_input) {
  symbol_col <- if ("SYMBOL" %in% names(res_input)) "SYMBOL" else "genesymb"
  rank_col <- if ("stat" %in% names(res_input) && !all(is.na(res_input$stat))) {
    "stat"
  } else if ("logFC_MMSE" %in% names(res_input)) {
    "logFC_MMSE"
  } else {
    "log2FoldChange"
  }

  res_input %>%
    dplyr::select(SYMBOL = all_of(symbol_col), rank_value = all_of(rank_col)) %>%
    na.omit() %>%
    group_by(SYMBOL) %>%
    summarise(rank_value = mean(rank_value, na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(rank_value)) %>%
    deframe()
}

# Run Hallmark GSEA independently for each saved comparison.
run_hallmark_gsea <- function(dge_list) {
  h_gene_sets <- msigdbr(species = "Homo sapiens", category = "H") %>%
    dplyr::select(gs_name, gene_symbol)

  ranks <- lapply(dge_list, get_ranks)
  lapply(ranks, function(x) {
    clusterProfiler::GSEA(sort(x, decreasing = TRUE), TERM2GENE = h_gene_sets, pvalueCutoff = 1.0, eps = 0)
  })
}

prepare_bubble_df <- function(gsea_output, fdr_cutoff = 0.05, top_n = 25) {
  gsea_df <- lapply(gsea_output, as_tibble)
  names(gsea_df) <- names(gsea_output)

  combined_genesets_up <- gsea_df %>%
    lapply(function(x) x %>% filter(p.adjust < fdr_cutoff, NES > 0) %>% arrange(desc(NES)) %>% pull(ID) %>% head(top_n)) %>%
    unlist() %>%
    unique()

  combined_genesets_down <- gsea_df %>%
    lapply(function(x) x %>% filter(p.adjust < fdr_cutoff, NES < 0) %>% arrange(NES) %>% pull(ID) %>% head(top_n)) %>%
    unlist() %>%
    unique()

  combined_genesets <- unique(c(combined_genesets_up, combined_genesets_down))

  plot_df <- bind_rows(lapply(names(gsea_df), function(nm) {
    gsea_df[[nm]] %>%
      mutate(group = nm) %>%
      dplyr::select(-any_of(c("leading_edge", "core_enrichment")))
  })) %>%
    filter(ID %in% combined_genesets) %>%
    mutate(group = factor(group, levels = names(gsea_df)))

  shared_parent_genesets <- plot_df %>%
    filter(grepl("Par", group), p.adjust < fdr_cutoff) %>%
    filter(duplicated(ID) | duplicated(ID, fromLast = TRUE)) %>%
    group_by(ID) %>%
    filter(all(NES > 0) | all(NES < 0)) %>%
    ungroup() %>%
    pull(ID) %>%
    unique()

  plot_df <- plot_df %>%
    filter(ID %in% shared_parent_genesets)

  ordered_ids <- plot_df %>%
    filter(group == "MCF7_Par_LY_5days_vs_DMSO_5days") %>%
    arrange(NES) %>%
    pull(ID) %>%
    unique()

  plot_df %>%
    mutate(
      ID = factor(ID, levels = ordered_ids),
      p.adjust = if_else(p.adjust > fdr_cutoff, NA_real_, p.adjust)
    )
}

plot_gsea_bubble <- function(plot_df) {
  ggplot(plot_df, aes(x = group, y = ID, col = NES)) +
    scale_color_gradientn(colours = c("dodgerblue3", "white", "orange")) +
    geom_point(aes(size = -log10(p.adjust))) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 6),
      axis.text.y = element_text(size = 6),
      axis.title = element_blank()
    )
}

MCF7_dge <- load_rds_or_warn(path("analysis", "RNAseq", "240108_MCF7M_Parent_sgRB1_DMSO_LY_RNAseq", "R_analysis", "3.DGE", "data", "rds", "DESeq2_list.rds"))
ZR751_dge <- load_rds_or_warn(path("analysis", "RNAseq", "20240910_ZR751_DMSO_LY_RNAseq", "R_analysis", "2.DE", "data", "DESeq2.rds"))

dge_list <- list()
if (!is.null(MCF7_dge)) {
  names(MCF7_dge) <- paste0("MCF7_", names(MCF7_dge))
  dge_list <- c(dge_list, MCF7_dge)
}
if (!is.null(ZR751_dge)) {
  ZR751_dge$SYMBOL <- ZR751_dge$genesymb
  dge_list[["ZR_Par_LY_4days_vs_DMSO_4days"]] <- ZR751_dge
}

selected_comparisons <- c(
  "MCF7_Par_LY_5days_vs_DMSO_5days",
  "MCF7_sgRB_LY_5days_vs_DMSO_5days",
  "ZR_Par_LY_4days_vs_DMSO_4days"
)
dge_list <- dge_list[intersect(selected_comparisons, names(dge_list))]

if (length(dge_list) > 0) {
  gsea_output <- run_hallmark_gsea(dge_list)
  names(gsea_output) <- names(dge_list)
  plot_df <- prepare_bubble_df(gsea_output)
  write.csv(plot_df, file.path(output_dir, "ER_expression_GSEA_bubble_data.csv"), row.names = FALSE)
  p <- plot_gsea_bubble(plot_df)
  ggsave(file.path(output_dir, "ER_expression_GSEA_bubble_plot.pdf"), plot = p, width = 8, height = 8)
}

message("Saved ER-expression GSEA outputs to: ", output_dir)
capture.output(sessionInfo(), file = file.path(output_dir, "sessionInfo.txt"))
