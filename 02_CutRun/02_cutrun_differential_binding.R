#!/usr/bin/env Rscript

# Configure DiffBind contrasts, run DESeq2-backed differential analysis and
# export complete, increased and reduced peak sets.

suppressPackageStartupMessages({
  library(DiffBind)
  library(tidyverse)
})

env_path <- function(name, default = "") {
  value <- Sys.getenv(name, default)
  if (!nzchar(value)) "" else normalizePath(value, winslash = "/", mustWork = FALSE)
}

output_dir <- env_path("CUTRUN_OUTPUT_DIR", file.path(Sys.getenv("PROJECT_DIR", "."), "code_availability/CutRun/output"))
contrasts_path <- env_path("CONTRASTS_CSV")
fdr_threshold <- as.numeric(Sys.getenv("FDR_THRESHOLD", "0.10"))
normalized_path <- file.path(output_dir, "dbObj_normalized.rds")

if (!file.exists(normalized_path)) stop("Run 01 first; missing ", normalized_path)
if (!nzchar(contrasts_path) || !file.exists(contrasts_path)) {
  stop("Set CONTRASTS_CSV to a CSV with object, factor, numerator, denominator, and min_members columns.")
}

db_normalized <- readRDS(normalized_path)
contrasts <- readr::read_csv(contrasts_path, show_col_types = FALSE)
required_columns <- c("object", "factor", "numerator", "denominator", "min_members")
missing_columns <- setdiff(required_columns, names(contrasts))
if (length(missing_columns)) stop("Missing contrast columns: ", paste(missing_columns, collapse = ", "))

unknown_objects <- setdiff(unique(contrasts$object), names(db_normalized))
if (length(unknown_objects)) stop("Contrasts reference unknown objects: ", paste(unknown_objects, collapse = ", "))

db_contrast <- db_normalized
# Add all CSV-defined contrasts to their corresponding DBA objects.
for (object_name in unique(contrasts$object)) {
  object_contrasts <- dplyr::filter(contrasts, object == object_name)
  for (i in seq_len(nrow(object_contrasts))) {
    current <- object_contrasts[i, ]
    db_contrast[[object_name]] <- dba.contrast(
      db_contrast[[object_name]],
      minMembers = current$min_members,
      contrast = c(current$factor, current$numerator, current$denominator)
    )
  }
}

# Analyze each contrasted object and retain standard chromosomes.
db_analyze <- lapply(db_contrast, dba.analyze, method = DBA_DESEQ2)
standard_chromosomes <- c(paste0("chr", 1:22), "chrX", "chrY")

db_report <- lapply(db_analyze, function(object) {
  report <- dba.report(
    object,
    contrast = 1,
    method = DBA_DESEQ2,
    th = 1,
    bCalled = TRUE,
    bCounts = TRUE,
    bNormalized = TRUE,
    fold = 0
  )
  report[as.character(seqnames(report)) %in% standard_chromosomes]
})

# Split significant peaks by the sign of the DiffBind fold change.
db_up <- lapply(db_report, function(x) x[!is.na(x$FDR) & x$FDR < fdr_threshold & x$Fold > 0])
db_down <- lapply(db_report, function(x) x[!is.na(x$FDR) & x$FDR < fdr_threshold & x$Fold < 0])

saveRDS(db_contrast, file.path(output_dir, "dbObj_contrast.rds"))
saveRDS(db_analyze, file.path(output_dir, "dbObj_analyze.rds"))
saveRDS(db_report, file.path(output_dir, "db_report.rds"))
saveRDS(db_up, file.path(output_dir, "db_up.rds"))
saveRDS(db_down, file.path(output_dir, "db_down.rds"))

peak_counts <- tibble(
  object = names(db_report),
  tested = vapply(db_report, length, integer(1)),
  up = vapply(db_up, length, integer(1)),
  down = vapply(db_down, length, integer(1)),
  fdr_threshold = fdr_threshold
)
readr::write_csv(peak_counts, file.path(output_dir, "differential_peak_counts.csv"))
capture.output(sessionInfo(), file = file.path(output_dir, "sessionInfo_02_differential.txt"))

message("Saved differential binding objects to ", output_dir)
