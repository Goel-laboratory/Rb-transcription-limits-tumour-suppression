#!/usr/bin/env Rscript

# Mass Spectrometry visualization and export.
# Inputs are Spectronaut objects and limma results; outputs include PCA and
# volcano plots.

library(tidyverse)
library(ggplot2)
library(ggrepel)
library(qs)

project_dir <- normalizePath(Sys.getenv("PROJECT_DIR", "."), winslash = "/", mustWork = FALSE)
if (!dir.exists(project_dir)) {
  stop("PROJECT_DIR does not exist: ", project_dir, ". Set PROJECT_DIR or run from the project root.")
}
output_dir <- normalizePath(
  Sys.getenv("MASSSPEC_OUTPUT_DIR", file.path(project_dir, "code_availability/MassSpec/output")),
  winslash = "/", mustWork = FALSE
)
input_dir <- output_dir
fs::dir_create(output_dir)

objects_path <- file.path(input_dir, "massspec_data_objects.rds")
diff_path <- file.path(input_dir, "massspec_diff_results.rds")

if (!file.exists(objects_path)) stop("Missing imported objects: ", objects_path)
if (!file.exists(diff_path)) stop("Missing differential results: ", diff_path)

imported <- readRDS(objects_path)
Spectro <- imported$Spectro
sampleInfo <- imported$sampleInfo
Spectro.im.log2 <- imported$Spectro.im.log2

diff_result <- readRDS(diff_path)

# Perform PCA on sample-wise protein abundances for one comparison.
plot_pca <- function(comparison_name) {
  data <- Spectro[[comparison_name]]
  data <- data %>% distinct(Protein, .keep_all = TRUE)
  data <- data %>% filter(Protein != "Nas", Protein != "NA")
  data_matrix <- data %>% column_to_rownames("Protein") %>% as.matrix()
  data_for_pca <- t(data_matrix)
  data_for_pca[is.nan(data_for_pca) | is.infinite(data_for_pca)] <- NA_real_
  keep <- apply(data_for_pca, 2, function(x) sum(!is.na(x)) >= 2)
  data_for_pca <- data_for_pca[, keep, drop = FALSE]
  data_for_pca <- apply(data_for_pca, 2, function(x) {
    x[is.na(x)] <- median(x, na.rm = TRUE)
    x
  })
  pca_result <- prcomp(data_for_pca, center = TRUE, scale. = TRUE)
  sample_info <- sampleInfo[match(rownames(data_for_pca), sampleInfo$SampleID), ]
  var_explained <- round(summary(pca_result)$importance[2, ] * 100, 1)
  pca_scores <- data.frame(Sample = rownames(data_for_pca), PC1 = pca_result$x[,1], PC2 = pca_result$x[,2], group = sample_info$group)

  ggplot(pca_scores, aes(x = PC1, y = PC2, color = group)) +
    geom_point(size = 4, alpha = 0.8) +
    labs(title = paste0("PCA: ", comparison_name), x = paste0("PC1 (", var_explained[1], "% )"), y = paste0("PC2 (", var_explained[2], "% )")) +
    theme_bw() +
    theme(legend.title = element_blank())
}

for (comparison in names(Spectro)) {
  p <- plot_pca(comparison)
  ggsave(file.path(output_dir, paste0("PCA_", comparison, ".png")), plot = p, width = 8, height = 6)
}

# Plot the Abemaciclib RB versus Abemaciclib IgG comparison as used in the
# original MassSpec analysis.
volcano_fig <- function(results_df, protlist = "RB1", protlist2 = "RB1"){
  # comparisons
  # diff_result$Abema_RBvDMSO_RB
  # protlist: up proteins
  # protlist2: c("your specific genes of interest")
  # Add labels to proteins of interest
  filtered_data <- results_df %>%
    dplyr::filter(rownames(results_df) %in% protlist2)
  # Set x-axis limits
  x_axis_limits <- DescTools::RoundTo(max(abs(results_df$logFC)), 5, ceiling)
  # Plot
  results_df %>%
    ggplot(aes(x = logFC, y = -log10(adj.P.Val))) +
    geom_point(alpha = 0.08, size = 0.4, color = "dodgerblue2") + 
    geom_point(data = results_df %>%
                 filter(rownames(.) %in% protlist), 
               aes(x = logFC, y = -log(adj.P.Val, 10)), 
               color = "orange", size = 0.3, alpha = 0.2) +
    # Highlight labeled points in orange
    geom_point(data = filtered_data, 
               aes(x = logFC, y = -log10(adj.P.Val)), 
               color = "darkorange3", size = 1.4) +
    geom_vline(xintercept = c(-1, 1), linetype = 'dashed', color = 'red3', linewidth = 0.2, alpha = 0.6) +
    geom_hline(yintercept = -log10(0.05), linetype = 'dashed', color = 'red3', linewidth = 0.2, alpha = 0.6) +
    scale_x_continuous('log2 Fold Change', 
                       limits = c(-16, 16),
                       #limits = c(-x_axis_limits, x_axis_limits), 
                       breaks = seq(-x_axis_limits, x_axis_limits, by = 4)) +
    ylab("-log10 adjusted P-value") +
    theme_bw(base_size = 6) +
    # Add labels
    geom_text_repel(data = filtered_data, aes(label = rownames(filtered_data)),
                    box.padding = 0.35, point.padding = 0.5, 
                    segment.color = 'black', size = 2, max.overlaps = 50) +
    theme(
      panel.grid.major = element_line(color = "grey", linetype = "dotted", linewidth = 0.15),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", size = 0.25, fill = NA),
      axis.title.x = element_text(size = 6, color = "black"),
      axis.title.y = element_text(size = 6, color = "black"),
      axis.text.x = element_text(size = 5, color = "black"),
      axis.text.y = element_text(size = 5, color = "black")
    )
}

protlist2 <- c("RB1", "E2F3", "E2F4", "TFDP1", "KDM5A", "ESR1")
Abema_RB_vs_Abema_IgG_volcano <- volcano_fig(diff_result$Abema_RBvAbema_IgG, protlist2 = protlist2)
ggsave(
  file.path(output_dir, "volcano_Abema_RBvAbema_IgG.png"),
  plot = Abema_RB_vs_Abema_IgG_volcano,
  width = 8,
  height = 6
)

message("Saved MassSpec visualization outputs to ", output_dir)
capture.output(sessionInfo(), file = file.path(output_dir, "sessionInfo_04_plots.txt"))
