# RB non-promoter function paper figures
#
# Cleaned from:
#   analysis_paper/02_RB_nonprom_function/01_RBup_ChiPEnrich_ORA/02_ChIPEnrich_ORA_nonpromoter.Rmd
#   analysis_paper/02_RB_nonprom_function/02_RBup_geneexprsup_hypergeometric/01_RBup_geneexprsup_hypergeometric.Rmd
#   analysis_paper/02_RB_nonprom_function/04_RBup_HiChIP_H3K27ac_hypergeometric/01_RBup_HiChIP_H3K27ac_hypergeometric.Rmd
#
# This script consumes saved ChIPEnrich and hypergeometric/ORA result objects
# and recreates the final dot-plot style summaries. It does not rerun
# ChIPEnrich, peak overlaps, or HiChIP processing.
#
# Method notes:
# - ChIPEnrich panels used Hallmark gene sets with `chipenrich`,
#   `method = "chipenrich"`, `locusdef = "nearest_tss"`, and genome `hg38`.
# - The hypergeometric/ORA outputs consumed here contain `GeneRatio`,
#   `FoldEnrichment`, and BH-adjusted `p.adjust`; figure panels used
#   `p.adjust < 0.05` as the display/significance cutoff.
# Outputs are publication-style ChIPEnrich and ORA summary plots.

suppressPackageStartupMessages({
  library(qs)
  library(tidyverse)
})

project_dir <- normalizePath(Sys.getenv("PROJECT_DIR", "."), winslash = "/", mustWork = FALSE)
output_dir <- file.path(Sys.getenv("FIGURES_OUTPUT_DIR", file.path(project_dir, "code_availability", "Figures", "output")), "02_RB_nonprom_function")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

path <- function(...) file.path(project_dir, ...)

load_qs_or_rds <- function(file) {
  if (!file.exists(file)) {
    warning("Missing input object: ", file, call. = FALSE)
    return(NULL)
  }
  if (grepl("\\.qs$", file, ignore.case = TRUE)) {
    qs::qread(file)
  } else {
    readRDS(file)
  }
}

select_chipenrich_result <- function(obj, key, new_name) {
  if (is.null(obj) || !key %in% names(obj)) return(NULL)
  out <- obj[[key]]
  out$source <- new_name
  out
}

# Standardize source-specific columns before combining comparisons.
standardise_chipenrich <- function(df, cell_line, treatment, test = "CE") {
  df %>%
    as_tibble() %>%
    dplyr::select(Geneset.ID, FDR, Odds.Ratio, gene_symbol) %>%
    mutate(Cell.line = cell_line, Treatment = treatment, Test = test)
}

plot_nonprom_chipenrich <- function(ce_np) {
  geneset_order <- c(
    "HALLMARK_ESTROGEN_RESPONSE_EARLY",
    "HALLMARK_MYOGENESIS",
    "HALLMARK_TGF_BETA_SIGNALING",
    "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
    "HALLMARK_APICAL_JUNCTION"
  )

  plot_df <- ce_np %>%
    filter(Geneset.ID %in% geneset_order) %>%
    mutate(
      Cell.line = factor(Cell.line, levels = c("MCF7", "ZR751", "PDX3837")),
      Geneset.ID = factor(Geneset.ID, levels = rev(geneset_order)),
      FDR_plot = pmax(FDR, .Machine$double.xmin),
      Sig = if_else(!is.na(FDR) & FDR < 0.05, "sig", "ns")
    )

  ggplot(plot_df, aes(x = Cell.line, y = Geneset.ID)) +
    facet_grid(~Treatment, scales = "free_x", space = "free_x") +
    geom_point(
      aes(size = -log10(FDR_plot), color = Sig, fill = if_else(Sig == "sig", Odds.Ratio, NA_real_)),
      alpha = 0.75,
      shape = 21,
      stroke = 0.5
    ) +
    scale_fill_gradientn(colors = c("orange", "red", "darkred"), na.value = "gray85") +
    scale_color_manual(values = c("sig" = "darkred", "ns" = "gray70")) +
    scale_size(range = c(1, 4.5)) +
    theme_bw() +
    theme(
      axis.text.y = element_text(color = "black", size = 6),
      axis.text.x.bottom = element_text(color = "black", size = 6),
      axis.title = element_blank(),
      legend.position = "right",
      legend.text = element_text(size = 5),
      legend.title = element_text(size = 5)
    )
}

gene_ratio_number <- function(x) {
  vapply(strsplit(as.character(x), "/", fixed = TRUE), function(parts) {
    as.numeric(parts[[1]]) / as.numeric(parts[[2]])
  }, numeric(1))
}

standardise_hypergeo <- function(df, cell_line, treatment, test = "Hypergeometric") {
  df %>%
    as_tibble() %>%
    mutate(
      GeneRatio_number = if ("GeneRatio" %in% names(.)) gene_ratio_number(GeneRatio) else NA_real_,
      Cell.line = cell_line,
      Treatment = treatment,
      Test = test
    ) %>%
    dplyr::select(ID, p.adjust, GeneRatio_number, FoldEnrichment, geneID, Cell.line, Treatment, Test, any_of("group"))
}

plot_nonprom_hypergeo <- function(hypergeo_df, fdr_cutoff = 0.05) {
  geneset_order <- c(
    "HALLMARK_ESTROGEN_RESPONSE_EARLY",
    "HALLMARK_MYOGENESIS",
    "HALLMARK_TGF_BETA_SIGNALING",
    "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
    "HALLMARK_APICAL_JUNCTION"
  )

  plot_df <- hypergeo_df %>%
    filter(ID %in% geneset_order) %>%
    mutate(
      ID = factor(ID, levels = rev(geneset_order)),
      Sig = if_else(!is.na(p.adjust) & p.adjust < fdr_cutoff, "sig", "ns")
    )

  ggplot(plot_df, aes(x = Cell.line, y = ID)) +
    facet_grid(~Treatment, scales = "free_x", space = "free_x") +
    geom_point(
      aes(size = -log10(p.adjust), color = Sig, fill = if_else(Sig == "sig", FoldEnrichment, NA_real_)),
      alpha = 0.75,
      shape = 21,
      stroke = 0.5
    ) +
    scale_fill_gradientn(colors = c("orange", "red", "darkred"), na.value = "gray85") +
    scale_color_manual(values = c("sig" = "darkred", "ns" = "gray70")) +
    scale_size(range = c(1, 4.5)) +
    theme_bw() +
    theme(axis.title = element_blank(), axis.text = element_text(size = 6), legend.position = "right")
}

chipenrich_inputs <- list(
  MCF7_abema = path("analysis", "CutRun", "240108_MCF7M_Parent_sgRB1_DMSO_LY_12targets", "R_analysis", "4.annotation", "chipenrich", "data", "rds", "db_up_ce_HM_df.rds"),
  ZR751_abema = path("analysis", "CutRun", "240416_MCF7M_ZR751_Parent_sgRB1_DMSO_LY_6targets_CnR", "R_analysis", "3.ZR_enrichment", "4.Chipenrich", "data", "qs", "ZR751_LY_prom_nonprom_chipenrich.qs"),
  PDX3837_abema = Sys.getenv("PDX3837_CHIPENRICH_QS", unset = ""),
  palbo = path("analysis", "CutRun", "241217_MCF7_ZR751_Parent_sgRB_palbo_RB_CnR", "R_analysis", "4.annotation", "chipenrich", "data", "qs", "db_up_ce_HM_df.qs")
)

MCF7_ce <- load_qs_or_rds(chipenrich_inputs$MCF7_abema)
ZR751_ce <- load_qs_or_rds(chipenrich_inputs$ZR751_abema)
PDX3837_ce <- if (nzchar(chipenrich_inputs$PDX3837_abema)) load_qs_or_rds(chipenrich_inputs$PDX3837_abema) else NULL
palbo_ce <- load_qs_or_rds(chipenrich_inputs$palbo)

ce_np <- list(
  if (!is.null(MCF7_ce)) standardise_chipenrich(MCF7_ce$Par_RB.other, "MCF7", "Abema"),
  if (!is.null(ZR751_ce)) standardise_chipenrich(ZR751_ce$ZR751_nonprom, "ZR751", "Abema"),
  if (!is.null(PDX3837_ce)) standardise_chipenrich(PDX3837_ce$PDX3837_nonprom, "PDX3837", "Abema"),
  if (!is.null(palbo_ce)) standardise_chipenrich(palbo_ce$MCF7M_parent_RB.other, "MCF7", "Palbo"),
  if (!is.null(palbo_ce)) standardise_chipenrich(palbo_ce$ZR751_parent_RB.other, "ZR751", "Palbo")
) %>%
  compact() %>%
  bind_rows()

if (nrow(ce_np) > 0) {
  write.csv(ce_np, file.path(output_dir, "RBup_nonpromoter_ChIPEnrich_summary.csv"), row.names = FALSE)
  p <- plot_nonprom_chipenrich(ce_np)
  ggsave(file.path(output_dir, "RBup_ChIPEnrich_nonprom.pdf"), plot = p, width = 5.2, height = 1.4)
}

message("Saved RB non-promoter function outputs to: ", output_dir)
capture.output(sessionInfo(), file = file.path(output_dir, "sessionInfo.txt"))
