#!/usr/bin/env Rscript

# Run HOMER motif discovery and optional ChIPEnrich on exported CUT&RUN BED sets.
# Outputs are written per peak set beneath the motif_enrichment directory.

suppressPackageStartupMessages(library(tidyverse))

env_path <- function(name, default = "") {
  value <- Sys.getenv(name, default)
  if (!nzchar(value)) "" else normalizePath(value, winslash = "/", mustWork = FALSE)
}

output_dir <- env_path("CUTRUN_OUTPUT_DIR", file.path(Sys.getenv("PROJECT_DIR", "."), "code_availability/CutRun/output"))
peak_dir <- env_path("PEAK_DIR", file.path(output_dir, "peaksets"))
analysis_dir <- file.path(output_dir, "motif_enrichment")
homer_genome <- Sys.getenv("HOMER_GENOME", "hg38")
homer_threads <- as.integer(Sys.getenv("HOMER_THREADS", "14"))
chipenrich_cores <- as.integer(Sys.getenv("CHIPENRICH_CORES", "5"))
run_chipenrich <- tolower(Sys.getenv("RUN_CHIPENRICH", "true")) %in% c("1", "true", "yes")

dir.create(analysis_dir, recursive = TRUE, showWarnings = FALSE)
peak_files <- list.files(peak_dir, pattern = "\\.bed$", full.names = TRUE)
peak_files <- peak_files[file.info(peak_files)$size > 0]
if (!length(peak_files)) stop("No non-empty BED files found in ", peak_dir)

homer_bin <- Sys.which("findMotifsGenome.pl")
if (!nzchar(homer_bin)) stop("findMotifsGenome.pl was not found on PATH.")

# Run HOMER independently for every non-empty peak set.
for (peak_path in peak_files) {
  peak_name <- tools::file_path_sans_ext(basename(peak_path))
  homer_out <- file.path(analysis_dir, "HOMER", peak_name)
  preparsed_dir <- file.path(homer_out, "preparsed")
  dir.create(preparsed_dir, recursive = TRUE, showWarnings = FALSE)

  status <- system2(
    homer_bin,
    args = c(
      shQuote(peak_path),
      shQuote(homer_genome),
      shQuote(homer_out),
      "-size", "200",
      "-preparsedDir", shQuote(preparsed_dir),
      "-p", as.character(homer_threads),
      "-len", "6,10,15,20"
    )
  )
  if (!identical(status, 0L)) stop("HOMER failed for ", peak_path, " with status ", status)
}

# Run ChIPEnrich only when requested and available.
if (run_chipenrich) {
  if (!requireNamespace("chipenrich", quietly = TRUE)) {
    stop("RUN_CHIPENRICH is true but the chipenrich package is not installed.")
  }

  chip_dir <- file.path(analysis_dir, "ChipEnrich")
  dir.create(chip_dir, recursive = TRUE, showWarnings = FALSE)
  selected_files <- peak_files[grepl("_promoter\\.bed$|_nonpromoter\\.bed$", peak_files)]

  for (peak_path in selected_files) {
    peak_name <- tools::file_path_sans_ext(basename(peak_path))
    locus_definition <- if (grepl("_promoter\\.bed$", peak_path)) "1kb" else "nearest_tss"
    result <- chipenrich::chipenrich(
      peaks = peak_path,
      genome = "hg38",
      method = "chipenrich",
      locusdef = locus_definition,
      genesets = "hallmark",
      n_cores = chipenrich_cores,
      out_name = NULL,
      qc = TRUE
    )
    saveRDS(result, file.path(chip_dir, paste0(peak_name, "_chipenrich.rds")))
    write_csv(as_tibble(result$results), file.path(chip_dir, paste0(peak_name, "_chipenrich.csv")))
  }
}

capture.output(sessionInfo(), file = file.path(output_dir, "sessionInfo_05_motif.txt"))
message("Saved motif and enrichment results to ", analysis_dir)
