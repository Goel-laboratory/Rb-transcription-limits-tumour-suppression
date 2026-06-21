#!/usr/bin/env Rscript

# Add ER promoter/non-promoter labels to the RB-annotated loops.
# Output: HC_integrated_overlap_df.rds for graph and paired-anchor analyses.

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(rtracklayer)
  library(tidyverse)
})

output_dir <- normalizePath(
  Sys.getenv("HICHIP_OUTPUT_DIR", file.path(Sys.getenv("PROJECT_DIR", "."), "code_availability/HiChIP_H3K27ac/output")),
  winslash = "/", mustWork = FALSE
)
input_path <- file.path(output_dir, "HC_RB_overlap_df.rds")
if (!file.exists(input_path)) stop("Run 02 first; missing ", input_path)
interactions <- readRDS(input_path)

er_promoter_path <- Sys.getenv("ER_PROMOTER_BED", "")
er_nonpromoter_path <- Sys.getenv("ER_NONPROMOTER_BED", "")
if (!nzchar(er_promoter_path) || !nzchar(er_nonpromoter_path) ||
    !file.exists(er_promoter_path) || !file.exists(er_nonpromoter_path)) {
  warning("ER BED inputs not supplied; copying RB-only object.")
  saveRDS(interactions, file.path(output_dir, "HC_integrated_overlap_df.rds"))
  quit(save = "no")
}

er_promoter <- import(er_promoter_path)
er_nonpromoter <- import(er_nonpromoter_path)

# As in the source analysis, promoter overlap takes priority.
classify_er <- function(anchor) case_when(
  overlapsAny(anchor, er_promoter) ~ "ERprom",
  overlapsAny(anchor, er_nonpromoter) ~ "ERnonprom",
  TRUE ~ "none"
)

# Use the same orientation-independent pairing scheme as the RB labels.
canonical_pair <- function(a, b) {
  order <- c(none = 1L, ERnonprom = 2L, ERprom = 3L)
  ifelse(order[a] >= order[b], paste(a, b, sep = "_"), paste(b, a, sep = "_"))
}

# Annotate both anchors while retaining all existing loop metadata.
integrated <- lapply(interactions, function(df) {
  anchor1 <- GRanges(df$chr1, IRanges(df$s1, df$e1))
  anchor2 <- GRanges(df$chr2, IRanges(df$s2, df$e2))
  df$peak1_ER <- classify_er(anchor1)
  df$peak2_ER <- classify_er(anchor2)
  df$interaction_ER_simple <- canonical_pair(df$peak1_ER, df$peak2_ER)
  df
})

saveRDS(integrated, file.path(output_dir, "HC_integrated_overlap_df.rds"))
capture.output(sessionInfo(), file = file.path(output_dir, "sessionInfo_03_ER_integration.txt"))
