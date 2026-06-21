#!/usr/bin/env Rscript

# Clean log2 Spectronaut abundance tables and apply the configured missing-value rule.
# Output is a list of protein-by-sample matrices for limma.

suppressPackageStartupMessages(library(tidyverse))

output_dir <- normalizePath(
  Sys.getenv("MASSSPEC_OUTPUT_DIR", file.path(Sys.getenv("PROJECT_DIR", "."), "code_availability/MassSpec/output")),
  winslash = "/", mustWork = FALSE
)
input_path <- file.path(output_dir, "massspec_data_objects.rds")
if (!file.exists(input_path)) stop("Run 01 first; missing ", input_path)

objects <- readRDS(input_path)
min_observed <- as.integer(Sys.getenv("MIN_OBSERVED_SAMPLES", "3"))
missing_strategy <- Sys.getenv("MISSING_STRATEGY", "retain")

# Keep one row per protein and require observations in enough samples.
clean_matrix <- function(df) {
  if (!"Protein" %in% names(df)) stop("Each abundance table must contain a Protein column.")
  df <- df %>% filter(!is.na(Protein), !duplicated(Protein))
  proteins <- df$Protein
  mat <- as.matrix(select(df, -Protein))
  storage.mode(mat) <- "numeric"
  mat[is.nan(mat) | is.infinite(mat)] <- NA_real_
  keep <- rowSums(!is.na(mat)) >= min_observed
  mat <- mat[keep, , drop = FALSE]
  rownames(mat) <- proteins[keep]

  if (missing_strategy == "zero") {
    mat[is.na(mat)] <- 0
  } else if (missing_strategy != "retain") {
    stop("MISSING_STRATEGY must be 'retain' or 'zero'.")
  }
  mat
}

# Apply identical cleaning to every Spectronaut comparison.
clean_matrices <- lapply(objects$Spectro.im.log2, clean_matrix)
saveRDS(clean_matrices, file.path(output_dir, "Spectro.im.log2.clean.rds"))
capture.output(sessionInfo(), file = file.path(output_dir, "sessionInfo_02_cleaning.txt"))
