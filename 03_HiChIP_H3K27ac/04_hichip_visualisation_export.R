#!/usr/bin/env Rscript

# Export RB loop-class summaries, a diagnostic PDF and WashU loop tracks.
# This reporting step does not alter the loop or hub definitions.

library(tidyverse)
library(GenomicRanges)

project_dir <- normalizePath(Sys.getenv("PROJECT_DIR", "."), winslash = "/", mustWork = FALSE)
if (!dir.exists(project_dir)) {
  stop("PROJECT_DIR does not exist: ", project_dir, ". Set PROJECT_DIR or run from the project root.")
}
output_dir <- normalizePath(
  Sys.getenv("HICHIP_OUTPUT_DIR", file.path(project_dir, "code_availability/HiChIP_H3K27ac/output")),
  winslash = "/", mustWork = FALSE
)
fs::dir_create(output_dir)

report_path <- file.path(output_dir, "HC_RB_overlap_df.rds")
if (!file.exists(report_path)) stop("Missing RB overlap HiChIP data at: ", report_path)

HC_RB_overlap <- readRDS(report_path)

# Save one interaction-class count table per condition.
save_summary <- function(df, condition) {
  summary <- df %>% count(interaction_simple) %>% arrange(desc(n))
  write_csv(summary, file.path(output_dir, paste0(condition, "_RB_interaction_counts.csv")))
  summary
}

summaries <- imap(HC_RB_overlap, save_summary)

plot_counts <- function(counts, condition) {
  ggplot(counts, aes(x = interaction_simple, y = n, fill = interaction_simple)) +
    geom_col() +
    labs(title = paste(condition, "RB HiChIP interaction anchor counts"), x = "RB overlap class", y = "Number of loops") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")
}

plots <- imap(summaries, plot_counts)

pdf(file.path(output_dir, "HC_RB_overlap_counts.pdf"), width = 8, height = 5)
walk(plots, print)
dev.off()

# Convert each paired loop directly to the WashU interaction format.
write_washu <- function(df, condition) {
  washudf <- df %>%
    mutate(
      s1 = pmax(s1 - 1L, 0L),
      second_cord = paste0(chr2, ":", s2, "-", e2),
      value1 = if ("cc" %in% names(.)) cc else 1,
      value2 = if ("Q.Value_Bias" %in% names(.)) Q.Value_Bias else 0
    ) %>%
    select(chr1, s1, e1, second_cord, value1, value2)
  write_tsv(washudf, file.path(output_dir, paste0(condition, "_WashU_hg38.bed")), col_names = FALSE)
}

imap(HC_RB_overlap, write_washu)

message("Export complete. Saved summaries and plots to ", output_dir)
capture.output(sessionInfo(), file = file.path(output_dir, "sessionInfo_04_export.txt"))
