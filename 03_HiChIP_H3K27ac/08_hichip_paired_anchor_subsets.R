#!/usr/bin/env Rscript

# Export promoter/non-promoter classes at two coordinate scales:
# the 5,001-bp HiChIP windows and the narrower RB peaks inside them.
# CSVs retain peak-window multiplicity; unique regions are also BED-exported.

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- if (length(script_arg)) dirname(normalizePath(sub("^--file=", "", script_arg[[1]]))) else "."
source(file.path(script_dir, "00_hichip_config.R"))

suppressPackageStartupMessages({
  library(dplyr)
  library(GenomicRanges)
  library(purrr)
  library(readr)
  library(rtracklayer)
})

input_path <- hichip_output("HC_integrated_overlap_df.rds")
if (!file.exists(input_path)) stop("Run steps 01-03 first; missing ", input_path)
loops <- readRDS(input_path)
rb_promoter_path <- hichip_env_file("RB_PROMOTER_BED")
rb_nonpromoter_path <- hichip_env_file("RB_NONPROMOTER_BED")
rb_peaks <- list(
  RBprom = rtracklayer::import(rb_promoter_path),
  RBnonprom = rtracklayer::import(rb_nonpromoter_path)
)

paired_dir <- hichip_output("paired_anchor_subsets")
dir.create(paired_dir, recursive = TRUE, showWarnings = FALSE)
window_dir <- file.path(paired_dir, "anchor_windows_5001bp")
peak_dir <- file.path(paired_dir, "RB_peaks_within_anchor_windows")
dir.create(window_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(peak_dir, recursive = TRUE, showWarnings = FALSE)

# Export complete loop rows for every orientation-independent class.
write_classes <- function(df, condition, prefix, col1, col2) {
  hichip_require_columns(df, c(col1, col2), paste(condition, prefix, "loop table"))
  classes <- sort(unique(c(df[[col1]], df[[col2]])))
  pairs <- expand.grid(class1 = classes, class2 = classes, stringsAsFactors = FALSE) %>%
    filter(class1 <= class2)

  map_dfr(seq_len(nrow(pairs)), function(i) {
    a <- pairs$class1[[i]]
    b <- pairs$class2[[i]]
    keep <- (df[[col1]] == a & df[[col2]] == b) | (df[[col1]] == b & df[[col2]] == a)
    subset <- df[keep, , drop = FALSE]
    filename <- paste(condition, prefix, a, "TO", b, sep = "_")
    write_csv(subset, file.path(paired_dir, paste0(filename, ".csv")))
    data.frame(condition = condition, factor = prefix, class1 = a, class2 = b, number_loops = nrow(subset))
  })
}

summary <- imap_dfr(loops, function(df, condition) {
  rb <- write_classes(df, condition, "RB", "peak1_RB", "peak2_RB")
  er <- if (all(c("peak1_ER", "peak2_ER") %in% names(df))) {
    write_classes(df, condition, "ER", "peak1_ER", "peak2_ER")
  } else {
    data.frame()
  }
  bind_rows(rb, er)
})

write_csv(summary, hichip_output("HiChIP_paired_anchor_subset_counts.csv"))

# Reproduce the source's two-coordinate-scale analysis:
#   1. retain the 5001-bp HiChIP anchor/window for every loop occurrence;
#   2. extract the narrower RB peak(s) that overlap each window;
#   3. copy the loop/window metadata onto every peak-window overlap.
#
# Multiplicity is intentionally retained. The same RB peak can occur more than
# once when its window participates in multiple loops, and one peak may overlap
# two adjacent HiChIP windows. This matches the active source analysis.
make_anchor_table <- function(df, condition) {
  hichip_require_columns(
    df,
    c("chr1", "s1", "e1", "chr2", "s2", "e2", "peak1_RB", "peak2_RB"),
    paste(condition, "loop table")
  )
  if (!"loop_id" %in% names(df)) df$loop_id <- paste0("loop_", seq_len(nrow(df)))

  loop_class <- mapply(
    function(a, b) paste(sort(c(a, b)), collapse = "_TO_"),
    df$peak1_RB,
    df$peak2_RB,
    USE.NAMES = FALSE
  )
  loop_metadata <- df
  names(loop_metadata) <- paste0("loop_", names(loop_metadata))

  bind_rows(
    bind_cols(
      tibble(
        condition,
        loop_id = df$loop_id,
        anchor_number = 1L,
        anchor_id = paste0(df$chr1, ":", df$s1, "-", df$e1),
        anchor_chr = df$chr1,
        anchor_start = df$s1,
        anchor_end = df$e1,
        anchor_RB_class = df$peak1_RB,
        partner_RB_class = df$peak2_RB,
        loop_RB_class = loop_class
      ),
      loop_metadata
    ),
    bind_cols(
      tibble(
        condition,
        loop_id = df$loop_id,
        anchor_number = 2L,
        anchor_id = paste0(df$chr2, ":", df$s2, "-", df$e2),
        anchor_chr = df$chr2,
        anchor_start = df$s2,
        anchor_end = df$e2,
        anchor_RB_class = df$peak2_RB,
        partner_RB_class = df$peak1_RB,
        loop_RB_class = loop_class
      ),
      loop_metadata
    )
  )
}

# Return the actual RB peaks with metadata from every overlapping window.
extract_peaks_from_windows <- function(anchor_table) {
  windows <- makeGRangesFromDataFrame(
    anchor_table,
    keep.extra.columns = TRUE,
    seqnames.field = "anchor_chr",
    start.field = "anchor_start",
    end.field = "anchor_end"
  )

  map_dfr(names(rb_peaks), function(rb_class) {
    class_windows <- windows[windows$anchor_RB_class == rb_class]
    if (!length(class_windows)) return(data.frame())

    hits <- findOverlaps(rb_peaks[[rb_class]], class_windows)
    if (!length(hits)) return(data.frame())

    peak_df <- as.data.frame(rb_peaks[[rb_class]][queryHits(hits)])
    window_df <- as.data.frame(mcols(class_windows[subjectHits(hits)]))
    names(peak_df)[1:5] <- c("peak_chr", "peak_start", "peak_end", "peak_width", "peak_strand")

    bind_cols(peak_df, window_df) %>%
      mutate(
        peak_RB_class = rb_class,
        peak_id = paste0(peak_chr, ":", peak_start, "-", peak_end),
        .before = 1
      )
  })
}

anchor_windows <- imap_dfr(loops, make_anchor_table)
rb_peaks_in_windows <- extract_peaks_from_windows(anchor_windows)

write_csv(anchor_windows, hichip_output("HiChIP_anchor_windows_5001bp.csv"))
write_csv(rb_peaks_in_windows, hichip_output("HiChIP_RB_peaks_within_anchor_windows.csv"))

window_counts <- anchor_windows %>%
  count(condition, loop_RB_class, anchor_RB_class, name = "number_anchor_windows")
peak_counts <- if (nrow(rb_peaks_in_windows)) {
  rb_peaks_in_windows %>%
    group_by(condition, loop_RB_class, peak_RB_class) %>%
    summarise(
      number_peak_window_overlaps = n(),
      number_unique_RB_peaks = n_distinct(peak_id),
      .groups = "drop"
    )
} else {
  tibble(
    condition = character(),
    loop_RB_class = character(),
    peak_RB_class = character(),
    number_peak_window_overlaps = integer(),
    number_unique_RB_peaks = integer()
  )
}
write_csv(window_counts, hichip_output("HiChIP_anchor_window_counts.csv"))
write_csv(peak_counts, hichip_output("HiChIP_RB_peak_within_window_counts.csv"))

# Write class-specific HiChIP-window tables and unique BED regions.
walk(
  split(anchor_windows, interaction(anchor_windows$condition, anchor_windows$loop_RB_class, anchor_windows$anchor_RB_class, drop = TRUE)),
  function(x) {
    label <- paste(x$condition[[1]], x$loop_RB_class[[1]], x$anchor_RB_class[[1]], sep = "_")
    write_csv(x, file.path(window_dir, paste0(label, "_windows5001bp.csv")))
    window_gr <- makeGRangesFromDataFrame(
      x,
      keep.extra.columns = FALSE,
      seqnames.field = "anchor_chr",
      start.field = "anchor_start",
      end.field = "anchor_end"
    )
    rtracklayer::export(unique(window_gr), file.path(window_dir, paste0(label, "_windows5001bp.bed")))
  }
)

# Write the corresponding narrow RB-peak tables and unique BED regions.
if (nrow(rb_peaks_in_windows)) {
  walk(
    split(
      rb_peaks_in_windows,
      interaction(
        rb_peaks_in_windows$condition,
        rb_peaks_in_windows$loop_RB_class,
        rb_peaks_in_windows$peak_RB_class,
        drop = TRUE
      )
    ),
    function(x) {
      label <- paste(x$condition[[1]], x$loop_RB_class[[1]], x$peak_RB_class[[1]], sep = "_")
      write_csv(x, file.path(peak_dir, paste0(label, "_RBpeaks.csv")))
      peak_gr <- makeGRangesFromDataFrame(
        x,
        keep.extra.columns = FALSE,
        seqnames.field = "peak_chr",
        start.field = "peak_start",
        end.field = "peak_end"
      )
      rtracklayer::export(unique(peak_gr), file.path(peak_dir, paste0(label, "_RBpeaks.bed")))
    }
  )
}

capture.output(sessionInfo(), file = hichip_output("sessionInfo_08_paired_anchor_subsets.txt"))
