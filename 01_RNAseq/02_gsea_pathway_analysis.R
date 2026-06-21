#!/usr/bin/env Rscript

# Run GSEA from saved DESeq2 result tables and export pathway summaries.
#
# Source analyses used MSigDB Hallmark gene sets and either clusterProfiler::GSEA
# or fgsea. This template defaults to clusterProfiler::GSEA, matching the most
# recent current-analysis Rmds; set method <- "fgsea" to use fgseaMultilevel.

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tibble)
  library(ggplot2)
  library(clusterProfiler)
  library(msigdbr)
  library(enrichplot)
})

deseq_results_rds <- Sys.getenv(
  "DESEQ_RESULTS_RDS",
  file.path(Sys.getenv("RNASEQ_DESEQ2_DIR", "code_availability/RNAseq/output/deseq2"), "deseq2_results_list.rds")
)
outdir <- Sys.getenv(
  "RNASEQ_GSEA_DIR",
  file.path(Sys.getenv("RNASEQ_OUTPUT_DIR", "code_availability/RNAseq/output"), "gsea")
)
method <- Sys.getenv("GSEA_METHOD", "clusterProfiler")
if (!file.exists(deseq_results_rds)) stop("Missing DESeq2 results: ", deseq_results_rds)
species <- "Homo sapiens"
msigdb_category <- "H"
custom_term2gene_file <- NA_character_

# Recent source GSEA Rmds ranked by DESeq2 Wald statistic. Older notes mention
# using shrunken log2 fold change when many statistic values were zero.
rank_column <- "stat"
fallback_rank_column <- "logFC_MMSE"
gene_symbol_col <- "SYMBOL"
padj_threshold <- 0.05
gsea_padj_threshold <- 0.05
min_gs_size <- 10
max_gs_size <- 500
set.seed(4)

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

result_list <- readRDS(deseq_results_rds)

# Use a supplied TERM2GENE file when configured, otherwise retrieve MSigDB.
load_term2gene <- function() {
  if (!is.na(custom_term2gene_file)) {
    x <- readr::read_csv(custom_term2gene_file, show_col_types = FALSE)
    stopifnot(all(c("gs_name", "gene_symbol") %in% colnames(x)))
    return(x %>% dplyr::select(gs_name, gene_symbol))
  }
  msigdbr::msigdbr(species = species, category = msigdb_category) %>%
    dplyr::select(gs_name, gene_symbol)
}

term2gene <- load_term2gene()
term2gene <- term2gene %>% dplyr::filter(!is.na(gs_name), !is.na(gene_symbol))
pathways <- split(term2gene$gene_symbol, term2gene$gs_name)

# Build decreasing rank vectors with deliberate duplicate-symbol handling.
make_rank_vector <- function(df) {
  rank_col <- rank_column
  if (!rank_col %in% colnames(df) || all(is.na(df[[rank_col]]))) {
    rank_col <- fallback_rank_column
  }
  if (!rank_col %in% colnames(df)) {
    stop("No usable rank column found. Tried: ", rank_column, " and ", fallback_rank_column)
  }
  if (!gene_symbol_col %in% colnames(df)) {
    if ("genesymb" %in% colnames(df)) {
      df[[gene_symbol_col]] <- df$genesymb
    } else {
      stop("Gene symbol column not found: ", gene_symbol_col)
    }
  }
  df %>%
    dplyr::filter(!is.na(.data[[gene_symbol_col]]), !is.na(.data[[rank_col]])) %>%
    dplyr::group_by(.data[[gene_symbol_col]]) %>%
    dplyr::summarise(rank_metric = mean(.data[[rank_col]], na.rm = TRUE), .groups = "drop") %>%
    dplyr::arrange(desc(rank_metric)) %>%
    tibble::deframe()
}

run_clusterprofiler <- function(ranks) {
  clusterProfiler::GSEA(
    geneList = sort(ranks, decreasing = TRUE),
    TERM2GENE = term2gene,
    pvalueCutoff = 1.0,
    exponent = 1,
    minGSSize = min_gs_size,
    maxGSSize = max_gs_size,
    eps = 0,
    seed = TRUE
  )
}

run_fgsea <- function(ranks) {
  if (!requireNamespace("fgsea", quietly = TRUE)) {
    stop("Package 'fgsea' is required when method <- 'fgsea'.")
  }
  fgsea::fgseaMultilevel(
    pathways = pathways,
    stats = sort(ranks, decreasing = TRUE),
    minSize = min_gs_size,
    maxSize = max_gs_size,
    eps = 0
  )
}

ranks <- lapply(result_list, make_rank_vector)
saveRDS(ranks, file.path(outdir, "gsea_rank_vectors.rds"))

# Run the selected GSEA implementation with common gene-set size limits.
gsea_output <- lapply(ranks, function(x) {
  if (method == "clusterProfiler") {
    run_clusterprofiler(x)
  } else if (method == "fgsea") {
    run_fgsea(x)
  } else {
    stop("method must be 'clusterProfiler' or 'fgsea'")
  }
})
saveRDS(gsea_output, file.path(outdir, paste0("gsea_output_", method, ".rds")))

as_gsea_tbl <- function(x) {
  if (inherits(x, "gseaResult")) {
    return(as_tibble(x@result))
  }
  as_tibble(x)
}

gsea_tables <- lapply(gsea_output, as_gsea_tbl)
for (nm in names(gsea_tables)) {
  tbl <- gsea_tables[[nm]]
  write.csv(tbl, file.path(outdir, paste0(nm, "_GSEA_all_results.csv")), row.names = FALSE)
  sig_col <- if ("p.adjust" %in% colnames(tbl)) "p.adjust" else "padj"
  sig <- tbl %>% dplyr::filter(!is.na(.data[[sig_col]]), .data[[sig_col]] < gsea_padj_threshold)
  write.csv(sig, file.path(outdir, paste0(nm, "_GSEA_significant.csv")), row.names = FALSE)
}

plot_df <- bind_rows(lapply(names(gsea_tables), function(nm) {
  gsea_tables[[nm]] %>% mutate(Comparison = nm)
}))

if (nrow(plot_df) > 0 && all(c("NES", "ID", "p.adjust") %in% colnames(plot_df))) {
  top_terms <- plot_df %>%
    dplyr::filter(!is.na(p.adjust)) %>%
    dplyr::group_by(Comparison) %>%
    dplyr::arrange(p.adjust, .by_group = TRUE) %>%
    dplyr::slice_head(n = 10) %>%
    dplyr::ungroup() %>%
    dplyr::pull(ID) %>%
    unique()
  p <- plot_df %>%
    dplyr::filter(ID %in% top_terms) %>%
    dplyr::mutate(log10padj = -log10(p.adjust),
                  Significance = ifelse(p.adjust < padj_threshold, "sig", "ns")) %>%
    ggplot(aes(x = Comparison, y = reorder(ID, NES))) +
    geom_point(aes(size = log10padj, fill = NES, color = Significance),
               shape = 21, alpha = 0.8) +
    scale_fill_gradient2(low = "dodgerblue4", mid = "white", high = "darkorange3") +
    scale_color_manual(values = c(sig = "black", ns = "grey75"), guide = "none") +
    labs(x = NULL, y = NULL, size = expression("-log"[10] * " adjusted P value")) +
    theme_bw(base_size = 10) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave(file.path(outdir, "GSEA_top_pathways_dotplot.pdf"), p, width = 8, height = 7)
}

if (method == "clusterProfiler") {
  for (nm in names(gsea_output)) {
    result_tbl <- as_tibble(gsea_output[[nm]]@result)
    top_ids <- result_tbl %>%
      dplyr::filter(!is.na(p.adjust)) %>%
      dplyr::arrange(p.adjust) %>%
      dplyr::slice_head(n = 2) %>%
      dplyr::pull(ID)
    if (length(top_ids) > 0) {
      p <- enrichplot::gseaplot2(gsea_output[[nm]], geneSetID = top_ids, title = nm, pvalue_table = TRUE)
      ggsave(file.path(outdir, paste0(nm, "_top_enrichment_plots.pdf")), p, width = 8, height = 5)
    }
  }
}

session_info <- utils::capture.output(sessionInfo())
writeLines(session_info, file.path(outdir, "sessionInfo_gsea.txt"))
