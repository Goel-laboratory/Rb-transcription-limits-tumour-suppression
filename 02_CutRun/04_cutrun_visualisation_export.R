#!/usr/bin/env Rscript

# Generate CUT&RUN peak-count and volcano plots plus coordinate-preserving tables.
# Input is the complete DiffBind report list from step 02.

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(tidyverse)
})

env_path <- function(name, default = "") {
  value <- Sys.getenv(name, default)
  if (!nzchar(value)) "" else normalizePath(value, winslash = "/", mustWork = FALSE)
}

output_dir <- env_path("CUTRUN_OUTPUT_DIR", file.path(Sys.getenv("PROJECT_DIR", "."), "code_availability/CutRun/output"))
plot_dir <- file.path(output_dir, "plots")
table_dir <- file.path(output_dir, "tables")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

report_path <- file.path(output_dir, "db_report.rds")
if (!file.exists(report_path)) stop("Run 02 first; missing ", report_path)
db_report <- readRDS(report_path)

fdr_threshold <- as.numeric(Sys.getenv("FDR_THRESHOLD", "0.10"))
safe_name <- function(x) gsub("[^A-Za-z0-9._-]+", "_", x)

# Convert genomic reports to plotting tables with significance classes.
report_to_table <- function(report) {
  as.data.frame(report) %>%
    as_tibble() %>%
    mutate(
      FDR_plot = pmax(FDR, .Machine$double.xmin),
      direction = case_when(
        !is.na(FDR) & FDR < fdr_threshold & Fold > 0 ~ "increased",
        !is.na(FDR) & FDR < fdr_threshold & Fold < 0 ~ "reduced",
        TRUE ~ "not_significant"
      )
    )
}

report_tables <- lapply(db_report, report_to_table)

# Summarize increased, reduced and non-significant peaks per comparison.
peak_counts <- imap_dfr(report_tables, function(x, object_name) {
  count(x, direction, name = "count") %>% mutate(object = object_name)
}) %>%
  complete(
    object = names(report_tables),
    direction = c("increased", "reduced", "not_significant"),
    fill = list(count = 0)
  )

write_csv(peak_counts, file.path(table_dir, "peak_counts_summary.csv"))

counts_plot <- peak_counts %>%
  filter(direction != "not_significant") %>%
  ggplot(aes(x = object, y = count, fill = direction)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  geom_text(
    aes(label = count),
    position = position_dodge(width = 0.8),
    vjust = -0.3,
    size = 3
  ) +
  scale_fill_manual(values = c(increased = "#D55E00", reduced = "#0072B2")) +
  labs(x = NULL, y = "Differential peak count", fill = NULL) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(
  file.path(plot_dir, "differential_peak_counts.pdf"),
  counts_plot,
  width = max(7, 0.65 * length(report_tables)),
  height = 5
)

# Export each report and its corresponding volcano plot.
for (object_name in names(report_tables)) {
  result <- report_tables[[object_name]]
  file_stub <- safe_name(object_name)
  write_csv(
    select(result, -FDR_plot),
    file.path(table_dir, paste0(file_stub, "_results.csv"))
  )

  volcano <- ggplot(result, aes(x = Fold, y = -log10(FDR_plot), color = direction)) +
    geom_point(size = 0.7, alpha = 0.8) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey50") +
    geom_hline(yintercept = -log10(fdr_threshold), linetype = "dashed", color = "grey50") +
    scale_color_manual(values = c(
      increased = "#D55E00",
      reduced = "#0072B2",
      not_significant = "grey70"
    )) +
    labs(
      title = object_name,
      x = "log2 fold change",
      y = "-log10 FDR",
      color = NULL
    ) +
    theme_bw() +
    theme(legend.position = "top")

  ggsave(file.path(plot_dir, paste0("volcano_", file_stub, ".pdf")), volcano, width = 7, height = 5)
}

capture.output(sessionInfo(), file = file.path(output_dir, "sessionInfo_04_visualisation.txt"))
message("Saved coordinate-preserving tables and plots to ", output_dir)
