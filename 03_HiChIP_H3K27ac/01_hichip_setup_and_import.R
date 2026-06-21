#!/usr/bin/env Rscript

# Import paired HiChIP loops and place all samples in hg38 coordinates.
# Inputs: HICHIP_MANIFEST and, for hg19 data, LIFTOVER_CHAIN.
# Output: HC_hg38_df.rds.

suppressPackageStartupMessages({
  library(data.table)
  library(GenomicRanges)
  library(rtracklayer)
  library(tidyverse)
})

output_dir <- normalizePath(
  Sys.getenv("HICHIP_OUTPUT_DIR", file.path(Sys.getenv("PROJECT_DIR", "."), "code_availability/HiChIP_H3K27ac/output")),
  winslash = "/", mustWork = FALSE
)
manifest_path <- Sys.getenv("HICHIP_MANIFEST", "")
chain_path <- Sys.getenv("LIFTOVER_CHAIN", "")
expected_width <- as.integer(Sys.getenv("EXPECTED_ANCHOR_WIDTH", "5001"))
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

if (!nzchar(manifest_path) || !file.exists(manifest_path)) {
  stop("Set HICHIP_MANIFEST to a CSV with condition, file, and genome_build columns.")
}
manifest <- readr::read_csv(manifest_path, show_col_types = FALSE)
if (length(setdiff(c("condition", "file", "genome_build"), names(manifest)))) {
  stop("HICHIP_MANIFEST requires condition, file, and genome_build.")
}
if (any(!file.exists(manifest$file))) stop("Manifest contains missing interaction files.")

interaction_columns <- c(
  "chr1", "s1", "e1", "chr2", "s2", "e2", "cc", "Coverage1", "isPeak1",
  "Bias1", "Mapp1", "GCContent1", "RESites1", "p", "Q.Value_Bias"
)

# Import the loop table and assign a stable row-level loop ID.
read_interactions <- function(path) {
  x <- fread(path, sep = "\t", header = FALSE, data.table = FALSE)
  if (ncol(x) < 6) stop("Interaction file has fewer than six columns: ", path)
  names(x)[seq_len(min(length(interaction_columns), ncol(x)))] <-
    interaction_columns[seq_len(min(length(interaction_columns), ncol(x)))]
  x %>% mutate(loop_id = paste0("loop_", row_number()))
}

# Lift both anchors together and retain only loops with a mapped pair.
lift_pairs <- function(df, chain) {
  anchor1 <- GRanges(df$chr1, IRanges(df$s1, df$e1), loop_id = df$loop_id)
  anchor2 <- GRanges(df$chr2, IRanges(df$s2, df$e2), loop_id = df$loop_id)
  lifted1 <- unlist(liftOver(anchor1, chain))
  lifted2 <- unlist(liftOver(anchor2, chain))

  if (!is.na(expected_width) && expected_width > 0) {
    lifted1 <- lifted1[width(lifted1) == expected_width]
    lifted2 <- lifted2[width(lifted2) == expected_width]
  }
  ids <- intersect(lifted1$loop_id, lifted2$loop_id)
  lifted1 <- lifted1[match(ids, lifted1$loop_id)]
  lifted2 <- lifted2[match(ids, lifted2$loop_id)]
  metadata <- df[match(ids, df$loop_id), , drop = FALSE]

  metadata %>% mutate(
    chr1 = as.character(seqnames(lifted1)),
    s1 = start(lifted1),
    e1 = end(lifted1),
    chr2 = as.character(seqnames(lifted2)),
    s2 = start(lifted2),
    e2 = end(lifted2)
  )
}

needs_liftover <- any(tolower(manifest$genome_build) == "hg19")
chain <- NULL
if (needs_liftover) {
  if (!nzchar(chain_path) || !file.exists(chain_path)) stop("hg19 input requires LIFTOVER_CHAIN.")
  chain <- import.chain(chain_path)
}

interactions <- setNames(vector("list", nrow(manifest)), manifest$condition)
# Process each condition independently while preserving manifest names.
for (i in seq_len(nrow(manifest))) {
  current <- read_interactions(manifest$file[i])
  build <- tolower(manifest$genome_build[i])
  interactions[[manifest$condition[i]]] <- if (build == "hg19") {
    lift_pairs(current, chain)
  } else if (build == "hg38") {
    current
  } else {
    stop("Unsupported genome_build: ", manifest$genome_build[i])
  }
}

saveRDS(interactions, file.path(output_dir, "HC_hg38_df.rds"))
capture.output(sessionInfo(), file = file.path(output_dir, "sessionInfo_01_import.txt"))
