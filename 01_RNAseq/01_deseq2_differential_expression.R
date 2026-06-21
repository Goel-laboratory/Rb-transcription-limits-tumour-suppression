#!/usr/bin/env Rscript

# Run DESeq2 differential expression from a count matrix and sample metadata.
# Outputs include normalized matrices, the fitted DESeqDataSet and contrast tables.
#
# Edit the configuration block, then run:
#   Rscript 02_deseq2_differential_expression.R

suppressPackageStartupMessages({
  library(DESeq2)
  library(readr)
  library(dplyr)
  library(tibble)
  library(stringr)
  library(ggplot2)
  library(ashr)
})

counts_file <- Sys.getenv("COUNTS_FILE", "")
metadata_file <- Sys.getenv("METADATA_FILE", "")
annotation_file <- Sys.getenv("ANNOTATION_FILE", NA_character_)
outdir <- Sys.getenv(
  "RNASEQ_DESEQ2_DIR",
  file.path(Sys.getenv("RNASEQ_OUTPUT_DIR", "code_availability/RNAseq/output"), "deseq2")
)
if (!nzchar(counts_file) || !file.exists(counts_file)) stop("Set COUNTS_FILE.")
if (!nzchar(metadata_file) || !file.exists(metadata_file)) stop("Set METADATA_FILE.")

sample_id_col <- "sample"
design_formula <- ~ condition
reference_levels <- list(condition = "control")

# Use either DESeq2 contrast triples, DESeq2 coefficient names, or numeric contrasts.
contrasts <- list(
  condition_treated_vs_control = c("condition", "treated", "control")
)
coef_contrasts <- list()
numeric_contrasts <- list()

# Low-count filtering thresholds used across source Rmd files commonly included:
#   rowSums(counts > 5) >= 2, >= 3, or >= 5 depending on sample number/design.
# The representative publication template uses at least 5 counts in at least
# min_samples samples. Set min_samples to the smallest biological group size.
min_count <- 5
min_samples <- 2

alpha <- 0.05
lfc_threshold <- 0
shrink_lfc <- TRUE
shrink_type <- "ashr"
write_rlog <- TRUE

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# Accept either a featureCounts table or a CSV count matrix.
read_featurecounts <- function(path) {
  fc <- readr::read_tsv(path, comment = "#", show_col_types = FALSE)
  annotation <- fc %>% dplyr::select(Geneid:Length)
  count_data <- fc %>%
    dplyr::select(-c(Chr:Length)) %>%
    as.data.frame()
  rownames(count_data) <- count_data$Geneid
  count_data$Geneid <- NULL
  count_data <- round(as.matrix(count_data))
  storage.mode(count_data) <- "integer"
  list(counts = count_data, feature_annotation = annotation)
}

load_counts <- function(path) {
  ext <- tools::file_ext(path)
  if (ext %in% c("csv")) {
    mat <- read.csv(path, row.names = 1, check.names = FALSE)
    mat <- round(as.matrix(mat))
    storage.mode(mat) <- "integer"
    return(mat)
  }
  read_featurecounts(path)
}

counts_input <- load_counts(counts_file)
if (is.list(counts_input)) {
  count_mat <- counts_input$counts
  feature_annotation <- counts_input$feature_annotation
} else {
  count_mat <- counts_input
  feature_annotation <- NULL
}

metadata <- readr::read_csv(metadata_file, show_col_types = FALSE) %>%
  as.data.frame()
stopifnot(sample_id_col %in% colnames(metadata))
rownames(metadata) <- metadata[[sample_id_col]]

# Remove common BAM/path suffixes so count columns match sample IDs.
clean_count_colnames <- function(x) {
  x %>%
    basename() %>%
    str_remove("\\.Aligned\\.sortedByCoord\\.out\\.bam$") %>%
    str_remove("Aligned\\.sortedByCoord\\.out\\.bam$") %>%
    str_remove("\\.bam$")
}

colnames(count_mat) <- clean_count_colnames(colnames(count_mat))

missing_in_counts <- setdiff(rownames(metadata), colnames(count_mat))
missing_in_metadata <- setdiff(colnames(count_mat), rownames(metadata))
if (length(missing_in_counts) > 0 || length(missing_in_metadata) > 0) {
  stop(
    "Sample names do not match.\n",
    "Missing in count matrix: ", paste(missing_in_counts, collapse = ", "), "\n",
    "Missing in metadata: ", paste(missing_in_metadata, collapse = ", ")
  )
}

count_mat <- count_mat[, rownames(metadata), drop = FALSE]
stopifnot(identical(colnames(count_mat), rownames(metadata)))

keep <- rowSums(count_mat > min_count) >= min_samples
filtered_counts <- count_mat[keep, , drop = FALSE]
message("Retained ", nrow(filtered_counts), " of ", nrow(count_mat), " genes after filtering.")

for (var in names(reference_levels)) {
  if (var %in% colnames(metadata)) {
    metadata[[var]] <- relevel(factor(metadata[[var]]), ref = reference_levels[[var]])
  }
}

dds <- DESeqDataSetFromMatrix(
  countData = filtered_counts,
  colData = metadata,
  design = design_formula
)
dds <- DESeq(dds)

# Save reusable normalized transformations before extracting contrasts.
normalized_counts <- counts(dds, normalized = TRUE)
write.csv(normalized_counts, file.path(outdir, "normalized_counts.csv"))
write.csv(log2(normalized_counts + 1), file.path(outdir, "log2_normalized_counts.csv"))

vst_counts <- assay(vst(dds, blind = FALSE))
write.csv(vst_counts, file.path(outdir, "vst_counts.csv"))

if (write_rlog) {
  rlog_counts <- assay(rlog(dds, blind = FALSE))
  write.csv(rlog_counts, file.path(outdir, "rlog_counts.csv"))
}

if (!is.na(annotation_file)) {
  gene_annotation <- readr::read_csv(annotation_file, show_col_types = FALSE)
} else {
  gene_annotation <- NULL
}

annotate_result <- function(res_mle, res_shrunk = NULL) {
  out <- as.data.frame(res_mle)
  if (!is.null(res_shrunk)) {
    out$logFC_MMSE <- res_shrunk$log2FoldChange
    out$lfcSE_MMSE <- res_shrunk$lfcSE
    if ("svalue" %in% colnames(as.data.frame(res_shrunk))) {
      out$svalue <- res_shrunk$svalue
    }
    names(out)[names(out) == "log2FoldChange"] <- "logFC_MLE"
    names(out)[names(out) == "lfcSE"] <- "lfcSE_MLE"
  }
  out <- tibble::rownames_to_column(out, "gene_id")
  out$gene_id_no_version <- sub("\\..*$", "", out$gene_id)
  if (!is.null(gene_annotation)) {
    join_column <- if ("gene_id" %in% names(gene_annotation)) "gene_id" else
      if ("gene_id_no_version" %in% names(gene_annotation)) "gene_id_no_version" else
        stop("Annotation requires gene_id or gene_id_no_version.")
    out <- dplyr::left_join(out, gene_annotation, by = join_column)
  }
  out %>% dplyr::arrange(padj)
}

# Run one requested contrast and optionally apply ashr LFC shrinkage.
run_one_result <- function(label, spec, type) {
  if (type == "contrast") {
    res_mle <- results(dds, contrast = spec, alpha = alpha, lfcThreshold = lfc_threshold)
    res_shrunk <- if (shrink_lfc) lfcShrink(dds, contrast = spec, type = shrink_type, res = res_mle, lfcThreshold = lfc_threshold) else NULL
  } else if (type == "coef") {
    res_mle <- results(dds, name = spec, alpha = alpha, lfcThreshold = lfc_threshold)
    res_shrunk <- if (shrink_lfc) lfcShrink(dds, coef = spec, type = shrink_type, res = res_mle, lfcThreshold = lfc_threshold) else NULL
  } else {
    res_mle <- results(dds, contrast = spec, alpha = alpha, lfcThreshold = lfc_threshold)
    res_shrunk <- if (shrink_lfc) lfcShrink(dds, contrast = spec, type = shrink_type, res = res_mle, lfcThreshold = lfc_threshold) else NULL
  }
  res <- annotate_result(res_mle, res_shrunk)
  write.csv(res, file.path(outdir, paste0(label, "_DESeq2_all_results.csv")), row.names = FALSE)
  sig <- res %>% dplyr::filter(!is.na(padj), padj < alpha)
  write.csv(sig, file.path(outdir, paste0(label, "_DESeq2_padj_lt_", alpha, ".csv")), row.names = FALSE)
  res
}

result_list <- list()
for (nm in names(contrasts)) {
  result_list[[nm]] <- run_one_result(nm, contrasts[[nm]], "contrast")
}
for (nm in names(coef_contrasts)) {
  result_list[[nm]] <- run_one_result(nm, coef_contrasts[[nm]], "coef")
}
for (nm in names(numeric_contrasts)) {
  result_list[[nm]] <- run_one_result(nm, numeric_contrasts[[nm]], "numeric")
}

saveRDS(dds, file.path(outdir, "deseq2_dds.rds"))
saveRDS(result_list, file.path(outdir, "deseq2_results_list.rds"))
writeLines(resultsNames(dds), file.path(outdir, "deseq2_results_names.txt"))

session_info <- utils::capture.output(sessionInfo())
writeLines(session_info, file.path(outdir, "sessionInfo_deseq2.txt"))
