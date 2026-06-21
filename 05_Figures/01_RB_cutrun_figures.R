# RB CUT&RUN paper figures
#
# This script reproduces the figure-generation logic from the paper Rmd using
# saved DiffBind report objects. It does not rerun DiffBind.
# Outputs are combined report objects and MA-style PDFs.

suppressPackageStartupMessages({
  library(DiffBind)
  library(GenomicRanges)
  library(qs)
  library(tidyverse)
})

project_dir <- normalizePath(Sys.getenv("PROJECT_DIR", "."), winslash = "/", mustWork = FALSE)
output_dir <- file.path(Sys.getenv("FIGURES_OUTPUT_DIR", file.path(project_dir, "code_availability", "Figures", "output")), "01_RB_cutrun")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

path <- function(...) file.path(project_dir, ...)

# Read either serialization format used by the source analyses.
load_qs_or_rds <- function(file) {
  if (!file.exists(file)) {
    warning("Missing input object: ", file, call. = FALSE)
    return(NULL)
  }
  if (grepl("\\.qs$", file, ignore.case = TRUE)) {
    qs::qread(file)
  } else {
    readRDS(file)
  }
}

load_named_reports <- function(input_map) {
  reports <- list()
  for (set_name in names(input_map)) {
    obj <- load_qs_or_rds(input_map[[set_name]]$file)
    if (is.null(obj)) next

    if (!is.null(input_map[[set_name]]$select)) {
      obj <- obj[input_map[[set_name]]$select]
    }
    if (!is.null(input_map[[set_name]]$rename)) {
      names(obj) <- input_map[[set_name]]$rename
    }
    reports <- c(reports, obj)
  }
  reports
}

# Standard MA plot used for replicated DiffBind comparisons.
make_MAplot <- function(x, plot_title, FDR_val = 0.1) {
  ma_df <- x %>% mcols() %>% as_tibble()
  up_number <- ma_df %>% filter(FDR < FDR_val, Fold > 0) %>% nrow()
  down_number <- ma_df %>% filter(FDR < FDR_val, Fold < 0) %>% nrow()
  y_max <- max(abs(ma_df$Fold), na.rm = TRUE) + 0.5

  ma_df %>%
    ggplot(aes(x = Conc, y = Fold)) +
    geom_point(data = filter(ma_df, FDR >= FDR_val), size = 1, colour = "dodgerblue3", alpha = 0.08) +
    geom_point(data = filter(ma_df, FDR < FDR_val), size = 1, colour = "orange", alpha = 0.32) +
    geom_hline(yintercept = 0, linetype = 2, alpha = 0.5, color = "red", linewidth = 1.2) +
    labs(x = "Mean expression (log2)", y = "log2 Fold Change") +
    ggtitle(plot_title) +
    theme_classic(base_line_size = 0.5, base_rect_size = 0.5) +
    theme(axis.title = element_text(size = 6), axis.text = element_text(size = 6)) +
    ylim(c(-y_max, y_max)) +
    annotate(
      "label",
      x = max(ma_df$Conc, na.rm = TRUE) * 0.8,
      y = c(-(y_max - 0.5), y_max - 0.5),
      label = c(paste0("Downregulated: ", down_number), paste0("Upregulated: ", up_number)),
      col = "red",
      size = 4
    )
}

make_PDX_MAplot <- function(x, plot_title, lfc_cutoff = 0.5) {
  ma_df <- x %>% mcols() %>% as_tibble()
  ma_df$Conc[ma_df$Conc > 400] <- 400
  up_number <- ma_df %>% filter(Fold > lfc_cutoff) %>% nrow()
  down_number <- ma_df %>% filter(Fold < -lfc_cutoff) %>% nrow()
  y_max <- max(abs(ma_df$Fold), na.rm = TRUE) + 0.5

  ma_df %>%
    ggplot(aes(x = Conc, y = Fold)) +
    geom_point(data = filter(ma_df, abs(Fold) <= lfc_cutoff), size = 1, colour = "dodgerblue3", alpha = 0.08) +
    geom_point(data = filter(ma_df, abs(Fold) > lfc_cutoff), size = 1, colour = "orange", alpha = 0.32) +
    geom_point(data = filter(ma_df, Conc > 399), shape = 17, size = 1, colour = "red", alpha = 1) +
    geom_hline(yintercept = 0, linetype = 2, alpha = 0.5, color = "red", linewidth = 1.2) +
    labs(x = "Mean expression (log2)", y = "log2 Fold Change") +
    ggtitle(plot_title) +
    theme_classic(base_line_size = 0.5, base_rect_size = 0.5) +
    theme(axis.title = element_text(size = 6), axis.text = element_text(size = 6)) +
    ylim(c(-y_max, y_max)) +
    xlim(0, 400) +
    annotate(
      "label",
      x = max(ma_df$Conc, na.rm = TRUE) * 0.8,
      y = c(-4, 4),
      label = c(paste0("Downregulated: ", down_number), paste0("Upregulated: ", up_number)),
      col = "red",
      size = 4
    )
}

# Saved report objects used by the paper Rmd. Paths are relative to the project
# root except for the external PDX3837 DynaTag object, which must be copied into
# the project or supplied by overriding `external_inputs` below.
report_inputs <- list(
  MCF7_LY = list(
    file = path("analysis", "CutRun", "240108_MCF7M_Parent_sgRB1_DMSO_LY_12targets", "R_analysis", "3.diffbind", "data", "rds", "db_report.rds"),
    select = c("Par_RB", "sgRB_RB"),
    rename = c("MCF7_Par_LY", "MCF7_sgRB_LY")
  ),
  MCF7_ZR751_PB = list(
    file = path("analysis", "CutRun", "241217_MCF7_ZR751_Parent_sgRB_palbo_RB_CnR", "R_analysis", "3.differential_analysis", "data", "qs", "db_report_separate.qs"),
    select = c("MCF7M_parent_RB", "ZR751_parent_RB"),
    rename = c("MCF7M_Par_PB", "ZR751_Par_PB")
  ),
  ZR751_LY = list(
    file = path("analysis_paper", "01_RB_cutrun", "01_MAplots", "data", "qs", "dbareport_ZR751_LY.qs"),
    rename = "ZR751_Par_LY"
  ),
  OVCAR_Kuramochi_INX = list(
    file = path("analysis", "CutRun", "241127_OVCAR3_Kuramochi_INX_DMSO_RB_CnR", "R_analysis", "antonioahn", "1.differential_binding_analysis", "data", "qs", "db_report.qs"),
    rename = c("Kuramochi_Par_INX", "OVCAR3_Par_INX")
  ),
  PDX4433_LY = list(
    file = path("analysis", "CutRun", "250406d_PDX4433-2_Veh_abema_03_RB_CnR", "R_analysis", "3.differential_analysis", "data", "qs", "db_report_separate_1x1.qs"),
    select = "RB",
    rename = "PDX4433_LY"
  )
)

# Optional collaborator-exported object used for the PDX3837 figure panel.
# The public object should use a neutral list entry name. If a different entry
# name is used, set PDX3837_DB_REPORT_NAME in the runtime environment.
pdx3837_report <- Sys.getenv("PDX3837_DB_REPORT_QS", unset = "")
if (nzchar(pdx3837_report)) {
  report_inputs$PDX3837_LY <- list(
    file = pdx3837_report,
    select = Sys.getenv("PDX3837_DB_REPORT_NAME", unset = "PDX3837_LY"),
    rename = "PDX3837_LY"
  )
}

dbreport_RB_comb <- load_named_reports(report_inputs)

# The original Rmd manually set PDX4433 concentration and FDR because this was a
# one-replicate comparison plotted by fold change rather than DiffBind FDR.
if ("PDX4433_LY" %in% names(dbreport_RB_comb)) {
  mcols(dbreport_RB_comb$PDX4433_LY)$Conc <-
    (mcols(dbreport_RB_comb$PDX4433_LY)$Abema_1_RB + mcols(dbreport_RB_comb$PDX4433_LY)$Vehicle_1_RB) / 2
  mcols(dbreport_RB_comb$PDX4433_LY)$FDR <- 0
}

saveRDS(dbreport_RB_comb, file.path(output_dir, "combined_RB_DiffBind_reports.rds"))

for (nm in names(dbreport_RB_comb)) {
  p <- if (grepl("^PDX", nm)) {
    make_PDX_MAplot(dbreport_RB_comb[[nm]], nm)
  } else {
    make_MAplot(dbreport_RB_comb[[nm]], nm, FDR_val = 0.1)
  }
  ggsave(file.path(output_dir, paste0(nm, "_MAplot.pdf")), plot = p, width = 6, height = 4)
}

message("Saved RB CUT&RUN MA plots to: ", output_dir)
capture.output(sessionInfo(), file = file.path(output_dir, "sessionInfo.txt"))
