# KDM5A/KDM5B mechanism paper figures
#
# This script starts from saved limma/Spectronaut outputs. It does not rerun
# raw mass-spectrometry processing or limma.
# Outputs are mass-spec volcano and selected-protein abundance panels.

suppressPackageStartupMessages({
  library(ggrepel)
  library(qs)
  library(tidyverse)
})

project_dir <- normalizePath(Sys.getenv("PROJECT_DIR", "."), winslash = "/", mustWork = FALSE)
output_dir <- file.path(Sys.getenv("FIGURES_OUTPUT_DIR", file.path(project_dir, "code_availability", "Figures", "output")), "04_KDM5A_KDM5B_mechanism")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

path <- function(...) file.path(project_dir, ...)

load_rds_or_warn <- function(file) {
  if (!file.exists(file)) {
    warning("Missing input object: ", file, call. = FALSE)
    return(NULL)
  }
  readRDS(file)
}

round_up_to <- function(x, to = 2) ceiling(x / to) * to

# Label proteins central to the proposed mechanism on the complete comparison.
simple_volcano <- function(comparison, labeled_genes = c("KMT2C", "KDM5B", "RB1", "E2F1", "E2F3", "E2F2", "TFDP1", "KDM5A", "ESR1")) {
  comparison$adj.P.Val.plot <- pmax(comparison$adj.P.Val, .Machine$double.xmin)
  filtered_data <- comparison %>% filter(rownames(comparison) %in% labeled_genes)
  x_axis_limits <- round_up_to(max(abs(comparison$logFC), na.rm = TRUE), 2)

  comparison %>%
    ggplot(aes(x = logFC, y = -log10(adj.P.Val.plot))) +
    geom_point(alpha = 0.2, size = 1.5, color = "dodgerblue2") +
    geom_point(data = filtered_data, aes(x = logFC, y = -log10(adj.P.Val)), color = "darkorange", size = 2) +
    geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "darkred") +
    geom_hline(yintercept = -log10(0.2), linetype = "dashed", color = "darkred") +
    scale_x_continuous("log2 Fold Change", limits = c(-x_axis_limits, x_axis_limits), breaks = seq(-x_axis_limits, x_axis_limits, by = 4)) +
    ylab("-log10(padj)") +
    theme_bw(base_size = 12) +
    geom_text_repel(
      data = filtered_data,
      aes(label = rownames(filtered_data)),
      box.padding = 0.35,
      point.padding = 0.5,
      segment.color = "lightgrey",
      size = 3,
      max.overlaps = 50
    )
}

# Highlight a supplied protein list without changing the limma statistics.
volcano_highlight <- function(comparison, protlist, labeled_genes) {
  comparison$adj.P.Val.plot <- pmax(comparison$adj.P.Val, .Machine$double.xmin)
  highlighted <- comparison %>% filter(rownames(comparison) %in% protlist)
  labeled <- comparison %>% filter(rownames(comparison) %in% labeled_genes)
  x_axis_limits <- round_up_to(max(abs(comparison$logFC), na.rm = TRUE), 2)

  comparison %>%
    ggplot(aes(x = logFC, y = -log10(adj.P.Val.plot))) +
    geom_point(alpha = 0.2, size = 1.5, color = "dodgerblue2") +
    geom_point(data = highlighted, color = "orange2", size = 1.5, shape = 10) +
    geom_point(data = labeled, color = "red", size = 2) +
    geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "red3") +
    geom_hline(yintercept = -log10(0.2), linetype = "dashed", color = "red3") +
    scale_x_continuous("log2 Fold Change", limits = c(-x_axis_limits, x_axis_limits), breaks = seq(-x_axis_limits, x_axis_limits, by = 4)) +
    ylab("-log10(padj)") +
    theme_bw(base_size = 12) +
    geom_text_repel(
      data = labeled,
      aes(label = rownames(labeled)),
      box.padding = 0.35,
      point.padding = 0.5,
      segment.color = "lightgrey",
      size = 3,
      max.overlaps = 50
    )
}

massspec_result_file <- path("analysis", "MassSpec", "20250109_RB_ER_IP_MS", "10.Spectronaut_data_analysis_outlier_excluded", "rds", "differential_analysis_result.rds")
diff_result <- load_rds_or_warn(massspec_result_file)

if (!is.null(diff_result)) {
  # The paper used Abema RB vs Abema IgG for the main mass-spec volcano figure.
  main_comparison <- diff_result$Abema_RBvAbema_IgG
  p <- simple_volcano(main_comparison)
  ggsave(file.path(output_dir, "Fig4a_Abema_RB_vs_Abema_IgG_volcano.pdf"), plot = p, width = 6, height = 5)

  sigprot1_table <- diff_result$Abema_RBvDMSO_RB %>%
    filter(logFC > 0.5, adj.P.Val < 0.2)
  sigprot2_table <- diff_result$Abema_RBvAbema_IgG %>%
    filter(logFC > 0.5, adj.P.Val < 0.2)
  sigprot_overlap <- intersect(rownames(sigprot1_table), rownames(sigprot2_table))

  write.csv(sigprot1_table, file.path(output_dir, "Abema_RB_vs_DMSO_RB_logFC_gt_0.5_padj_lt_0.2.csv"))
  write.csv(sigprot2_table, file.path(output_dir, "Abema_RB_vs_Abema_IgG_logFC_gt_0.5_padj_lt_0.2.csv"))
  write.csv(data.frame(Protein = sigprot_overlap), file.path(output_dir, "Abema_RB_enriched_overlap_proteins.csv"), row.names = FALSE)

  p_highlight <- volcano_highlight(
    diff_result$Abema_RBvAbema_IgG,
    protlist = sigprot_overlap,
    labeled_genes = c("KDM5B", "KMT2C", "KDM5A", "RB1", "ESR1")
  )
  ggsave(file.path(output_dir, "Abema_RB_vs_Abema_IgG_overlap_highlight_volcano.pdf"), plot = p_highlight, width = 6, height = 5)
}

message("Saved KDM5A/KDM5B mechanism outputs to: ", output_dir)
capture.output(sessionInfo(), file = file.path(output_dir, "sessionInfo.txt"))
