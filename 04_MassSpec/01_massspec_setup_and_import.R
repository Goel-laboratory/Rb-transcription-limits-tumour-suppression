#!/usr/bin/env Rscript

# Import saved Spectronaut abundance and sample-information objects.
# Output is one bundled RDS used by downstream MassSpec scripts.

suppressPackageStartupMessages({
  library(qs)
  library(tidyverse)
})

env_path <- function(name, default = "") {
  value <- Sys.getenv(name, default)
  if (!nzchar(value)) "" else normalizePath(value, winslash = "/", mustWork = FALSE)
}

output_dir <- env_path("MASSSPEC_OUTPUT_DIR", file.path(Sys.getenv("PROJECT_DIR", "."), "code_availability/MassSpec/output"))
spectro_path <- env_path("SPECTRO_QS")
sample_info_path <- env_path("SAMPLE_INFO_QS")
log2_path <- env_path("SPECTRO_LOG2_QS")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

required <- c(SPECTRO_QS = spectro_path, SAMPLE_INFO_QS = sample_info_path, SPECTRO_LOG2_QS = log2_path)
# Validate all external objects before attempting deserialization.
missing <- names(required)[!nzchar(required) | !file.exists(required)]
if (length(missing)) stop("Set valid input paths for: ", paste(missing, collapse = ", "))

objects <- list(
  Spectro = qread(spectro_path),
  sampleInfo = qread(sample_info_path),
  Spectro.im.log2 = qread(log2_path)
)

# Create a consistent treatment/IP group for downstream figures.
if ("SampleID" %in% names(objects$sampleInfo) &&
    all(c("Treatment", "IP") %in% names(objects$sampleInfo))) {
  objects$sampleInfo <- objects$sampleInfo %>%
    mutate(group = paste(Treatment, IP, sep = "_"))
}

saveRDS(objects, file.path(output_dir, "massspec_data_objects.rds"))
capture.output(sessionInfo(), file = file.path(output_dir, "sessionInfo_01_import.txt"))
