#!/usr/bin/env Rscript

# Import DiffBind objects, count reads over consensus peaks and normalize them.
# Inputs are a saved DBA object/list or sample sheet.

suppressPackageStartupMessages({
  library(DiffBind)
  library(tidyverse)
})

env_path <- function(name, default = "") {
  value <- Sys.getenv(name, default)
  if (!nzchar(value)) "" else normalizePath(value, winslash = "/", mustWork = FALSE)
}

output_dir <- env_path("CUTRUN_OUTPUT_DIR", file.path(Sys.getenv("PROJECT_DIR", "."), "code_availability/CutRun/output"))
objects_rds <- env_path("DBA_OBJECTS_RDS")
sample_sheet_path <- env_path("SAMPLE_SHEET")
group_column <- Sys.getenv("GROUP_COLUMN", "Genotype_Antibody")
histone_pattern <- Sys.getenv("HISTONE_PATTERN", "H3K4me2|H3K27ac")
tf_summit <- as.integer(Sys.getenv("TF_SUMMIT", "200"))
histone_summit <- as.integer(Sys.getenv("HISTONE_SUMMIT", "500"))

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Standardize a single DBA object or list to one validated named list.
as_dba_list <- function(x) {
  if (inherits(x, "DBA")) {
    x <- list(all_samples = x)
  }
  if (!is.list(x) || !all(vapply(x, inherits, logical(1), what = "DBA"))) {
    stop("Input must be a DBA object or a list containing only DBA objects.")
  }
  if (is.null(names(x)) || any(!nzchar(names(x)))) {
    names(x) <- paste0("object_", seq_along(x))
  }
  x
}

if (nzchar(objects_rds)) {
  if (!file.exists(objects_rds)) stop("DBA_OBJECTS_RDS not found: ", objects_rds)
  db_objects <- as_dba_list(readRDS(objects_rds))
} else if (nzchar(sample_sheet_path)) {
  if (!file.exists(sample_sheet_path)) stop("SAMPLE_SHEET not found: ", sample_sheet_path)
  sample_sheet <- readr::read_csv(sample_sheet_path, show_col_types = FALSE)

  if (group_column %in% names(sample_sheet)) {
    grouped_sheets <- split(sample_sheet, sample_sheet[[group_column]], drop = TRUE)
    grouped_sheets <- grouped_sheets[vapply(grouped_sheets, nrow, integer(1)) > 0]
    db_objects <- lapply(grouped_sheets, function(sheet) {
      dba(sampleSheet = as.data.frame(sheet), minOverlap = 1)
    })
  } else {
    db_objects <- list(all_samples = dba(sampleSheet = as.data.frame(sample_sheet), minOverlap = 1))
  }
} else {
  stop("Set either DBA_OBJECTS_RDS or SAMPLE_SHEET.")
}

summit_widths <- ifelse(
  grepl(histone_pattern, names(db_objects), ignore.case = TRUE),
  histone_summit,
  tf_summit
)
names(summit_widths) <- names(db_objects)

message("Counting peaks with summit widths: ",
        paste(names(summit_widths), summit_widths, sep = "=", collapse = ", "))

# Use source-supported summit widths for TF and histone-mark objects.
db_count <- Map(function(object, summit) {
  dba.count(
    object,
    minOverlap = 1,
    summit = summit,
    bUseSummarizeOverlaps = TRUE,
    filter = 0
  )
}, db_objects, summit_widths)

# Apply the DESeq2/RLE background-library normalization used downstream.
db_normalized <- lapply(db_count, function(object) {
  dba.normalize(
    object,
    method = DBA_DESEQ2,
    normalize = DBA_NORM_RLE,
    library = DBA_LIBSIZE_BACKGROUND,
    background = TRUE
  )
})

saveRDS(db_count, file.path(output_dir, "dbObj_count.rds"))
saveRDS(db_normalized, file.path(output_dir, "dbObj_normalized.rds"))
write.csv(
  data.frame(object = names(summit_widths), summit = unname(summit_widths)),
  file.path(output_dir, "summit_widths.csv"),
  row.names = FALSE
)
capture.output(sessionInfo(), file = file.path(output_dir, "sessionInfo_01_setup.txt"))

message("Saved counted and normalized DiffBind objects to ", output_dir)
