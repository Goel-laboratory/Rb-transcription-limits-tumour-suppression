#!/usr/bin/env Rscript

# Label both anchors by RB promoter/non-promoter peak overlap.
# Promoter overlap has priority when a 5,001-bp anchor overlaps both sets.
# Output: HC_RB_overlap_df.rds plus loop-class counts.

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(rtracklayer)
  library(tidyverse)
})

output_dir <- normalizePath(
  Sys.getenv("HICHIP_OUTPUT_DIR", file.path(Sys.getenv("PROJECT_DIR", "."), "code_availability/HiChIP_H3K27ac/output")),
  winslash = "/", mustWork = FALSE
)
promoter_bed <- Sys.getenv("RB_PROMOTER_BED", "")
nonpromoter_bed <- Sys.getenv("RB_NONPROMOTER_BED", "")
for (path in c(promoter_bed, nonpromoter_bed)) {
  if (!nzchar(path) || !file.exists(path)) stop("Set valid RB_PROMOTER_BED and RB_NONPROMOTER_BED.")
}

input_path <- file.path(output_dir, "HC_hg38_df.rds")
if (!file.exists(input_path)) stop("Run 01 first; missing ", input_path)
interactions <- readRDS(input_path)
rb_promoter <- import(promoter_bed)
rb_nonpromoter <- import(nonpromoter_bed)

# Classify the full HiChIP anchor window, not the narrower RB peak itself.
classify_anchor <- function(anchor) {
  case_when(
    overlapsAny(anchor, rb_promoter) ~ "RBprom",
    overlapsAny(anchor, rb_nonpromoter) ~ "RBnonprom",
    TRUE ~ "none"
  )
}

# Store an orientation-independent label for the pair of anchor classes.
canonical_pair <- function(a, b) {
  order <- c(none = 1L, RBnonprom = 2L, RBprom = 3L)
  ifelse(order[a] >= order[b], paste(a, b, sep = "_"), paste(b, a, sep = "_"))
}

# Keep the original loop rows and add an RB label to each anchor.
annotated <- lapply(interactions, function(df) {
  anchor1 <- GRanges(df$chr1, IRanges(df$s1, df$e1))
  anchor2 <- GRanges(df$chr2, IRanges(df$s2, df$e2))
  df$peak1_RB <- classify_anchor(anchor1)
  df$peak2_RB <- classify_anchor(anchor2)
  df$interaction_simple <- canonical_pair(df$peak1_RB, df$peak2_RB)
  df
})

saveRDS(annotated, file.path(output_dir, "HC_RB_overlap_df.rds"))
imap_dfr(annotated, ~ count(.x, interaction_simple) %>% mutate(condition = .y)) %>%
  write_csv(file.path(output_dir, "HC_RB_overlap_counts.csv"))
capture.output(sessionInfo(), file = file.path(output_dir, "sessionInfo_02_RB_annotation.txt"))
