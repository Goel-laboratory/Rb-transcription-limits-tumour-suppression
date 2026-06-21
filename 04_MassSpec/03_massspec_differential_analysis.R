#!/usr/bin/env Rscript

# Run limma differential-abundance analysis from cleaned protein matrices.
# Contrasts and sample groups are supplied in MASSSPEC_CONTRASTS_CSV.

suppressPackageStartupMessages({
  library(limma)
  library(readr)
})

output_dir <- normalizePath(
  Sys.getenv("MASSSPEC_OUTPUT_DIR", file.path(Sys.getenv("PROJECT_DIR", "."), "code_availability/MassSpec/output")),
  winslash = "/", mustWork = FALSE
)
matrix_path <- file.path(output_dir, "Spectro.im.log2.clean.rds")
contrast_path <- Sys.getenv("MASSSPEC_CONTRASTS_CSV", "")
if (!file.exists(matrix_path)) stop("Run 02 first; missing ", matrix_path)
if (!nzchar(contrast_path) || !file.exists(contrast_path)) {
  stop("Set MASSSPEC_CONTRASTS_CSV with columns matrix, sample, group, and contrast.")
}

matrices <- readRDS(matrix_path)
config <- read_csv(contrast_path, show_col_types = FALSE)
required <- c("matrix", "sample", "group", "contrast")
if (length(setdiff(required, names(config)))) {
  stop("Contrast CSV requires: ", paste(required, collapse = ", "))
}

results <- list()
# Build and fit one no-intercept design for each configured matrix.
for (matrix_name in unique(config$matrix)) {
  if (!matrix_name %in% names(matrices)) stop("Unknown matrix: ", matrix_name)
  current <- config[config$matrix == matrix_name, ]
  mat <- matrices[[matrix_name]]
  if (!all(current$sample %in% colnames(mat))) stop("Samples missing from matrix ", matrix_name)
  mat <- mat[, current$sample, drop = FALSE]

  group <- factor(current$group)
  design <- model.matrix(~ 0 + group)
  colnames(design) <- levels(group)
  fit <- lmFit(mat, design, na.action = na.exclude)

# Test every requested contrast using empirical-Bayes moderation.
  for (contrast_string in unique(current$contrast)) {
    contrast_matrix <- makeContrasts(contrasts = contrast_string, levels = design)
    fit2 <- eBayes(contrasts.fit(fit, contrast_matrix))
    label <- paste(matrix_name, gsub("[^A-Za-z0-9._-]+", "_", contrast_string), sep = "__")
    results[[label]] <- topTable(fit2, adjust.method = "BH", number = Inf, sort.by = "P")
  }
}

saveRDS(results, file.path(output_dir, "massspec_diff_results.rds"))
for (name in names(results)) {
  write.csv(results[[name]], file.path(output_dir, paste0(name, ".csv")), row.names = TRUE)
}
capture.output(sessionInfo(), file = file.path(output_dir, "sessionInfo_03_limma.txt"))
