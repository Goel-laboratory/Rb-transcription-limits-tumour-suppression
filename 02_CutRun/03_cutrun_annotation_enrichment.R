#!/usr/bin/env Rscript

# Annotate CUT&RUN peaks with genomic class, promoter status and gene symbols.
# Outputs include annotated objects and promoter/non-promoter BED peak sets.

suppressPackageStartupMessages({
  library(ChIPseeker)
  library(GenomicFeatures)
  library(GenomicRanges)
  library(rtracklayer)
  library(tidyverse)
})

env_path <- function(name, default = "") {
  value <- Sys.getenv(name, default)
  if (!nzchar(value)) "" else normalizePath(value, winslash = "/", mustWork = FALSE)
}

output_dir <- env_path("CUTRUN_OUTPUT_DIR", file.path(Sys.getenv("PROJECT_DIR", "."), "code_availability/CutRun/output"))
annotation_dir <- env_path("ANNOTATION_DIR", file.path(Sys.getenv("PROJECT_DIR", "."), "resources/annotation"))
gtf_file <- env_path("GTF_FILE", file.path(annotation_dir, "gencode.v35.basic.annotation.gtf"))
gene_map_file <- env_path("GENE_MAP_FILE", file.path(annotation_dir, "mart_aug2020.RData"))
bed_dir <- file.path(output_dir, "peaksets")
dir.create(bed_dir, recursive = TRUE, showWarnings = FALSE)

# Load required outputs from the differential-binding stage.
read_required <- function(filename) {
  path <- file.path(output_dir, filename)
  if (!file.exists(path)) stop("Run 02 first; missing ", path)
  readRDS(path)
}

db_up <- read_required("db_up.rds")
db_down <- read_required("db_down.rds")
db_report <- read_required("db_report.rds")

if (file.exists(gtf_file)) {
  txdb_hg38 <- makeTxDbFromGFF(gtf_file, format = "gtf")
} else if (requireNamespace("TxDb.Hsapiens.UCSC.hg38.knownGene", quietly = TRUE)) {
  txdb_hg38 <- getExportedValue(
    "TxDb.Hsapiens.UCSC.hg38.knownGene",
    "TxDb.Hsapiens.UCSC.hg38.knownGene"
  )
} else {
  stop("Set GTF_FILE or install TxDb.Hsapiens.UCSC.hg38.knownGene.")
}

gene_map <- NULL
if (file.exists(gene_map_file)) {
  loaded_names <- load(gene_map_file)
  if ("t2g_aug2020" %in% loaded_names) {
    gene_map <- get("t2g_aug2020")
  } else {
    warning("GENE_MAP_FILE did not contain t2g_aug2020; SYMBOL annotation skipped.")
  }
}

# Collapse detailed ChIPseeker labels into the classes used in figures.
simplify_annotation <- function(annotation) {
  case_when(
    grepl("Exon", annotation) ~ "Exon",
    grepl("Intron", annotation) ~ "Intron",
    grepl("Promoter", annotation) ~ "Promoter",
    grepl("Downstream|Distal Intergenic", annotation) ~ "Intergenic",
    TRUE ~ annotation
  )
}

expand_intergenic <- function(annotation_simple, distance_to_tss) {
  case_when(
    annotation_simple != "Intergenic" ~ annotation_simple,
    abs(distance_to_tss) < 10000 ~ "Intergenic_Proximal",
    abs(distance_to_tss) < 100000 ~ "Intergenic_Distal",
    TRUE ~ "Intergenic_Desert"
  )
}

# Add promoter status, intergenic-distance groups and optional symbols.
annotate_one <- function(peaks) {
  annotation <- annotatePeak(
    peaks,
    tssRegion = c(-2000, 500),
    TxDb = txdb_hg38,
    genomicAnnotationPriority = c("Promoter", "5UTR", "3UTR", "Exon", "Intron", "Intergenic")
  )

  annotation@anno$annotation_simple <- simplify_annotation(annotation@anno$annotation)
  annotation@anno$annotation_expand <- expand_intergenic(
    annotation@anno$annotation_simple,
    annotation@anno$distanceToTSS
  )
  annotation@anno$annotation_promoter <- ifelse(
    annotation@anno$annotation_simple == "Promoter",
    "Promoter",
    "Other"
  )

  if (!is.null(gene_map) && all(c("ext_gene", "ens_gene_ver") %in% names(gene_map))) {
    annotation@anno$SYMBOL <- gene_map$ext_gene[
      match(annotation@anno$geneId, gene_map$ens_gene_ver)
    ]
  }
  annotation
}

annotate_list <- function(peaks) {
  peaks <- peaks[vapply(peaks, length, integer(1)) > 0]
  lapply(peaks, annotate_one)
}

db_up_peakAnno <- annotate_list(db_up)
db_down_peakAnno <- annotate_list(db_down)
db_report_peakAnno <- annotate_list(db_report)

saveRDS(db_up_peakAnno, file.path(output_dir, "db_up_peakAnno.rds"))
saveRDS(db_down_peakAnno, file.path(output_dir, "db_down_peakAnno.rds"))
saveRDS(db_report_peakAnno, file.path(output_dir, "db_report_peakAnno.rds"))

safe_name <- function(x) gsub("[^A-Za-z0-9._-]+", "_", x)

# Export complete, promoter and non-promoter BED files for each peak set.
export_peak_sets <- function(annotation_list, direction) {
  for (object_name in names(annotation_list)) {
    peaks <- annotation_list[[object_name]]@anno
    prefix <- paste(safe_name(object_name), direction, sep = "_")
    export(peaks, file.path(bed_dir, paste0(prefix, ".bed")), format = "BED")
    export(
      peaks[peaks$annotation_promoter == "Promoter"],
      file.path(bed_dir, paste0(prefix, "_promoter.bed")),
      format = "BED"
    )
    export(
      peaks[peaks$annotation_promoter == "Other"],
      file.path(bed_dir, paste0(prefix, "_nonpromoter.bed")),
      format = "BED"
    )
  }
}

export_peak_sets(db_up_peakAnno, "up")
export_peak_sets(db_down_peakAnno, "down")
capture.output(sessionInfo(), file = file.path(output_dir, "sessionInfo_03_annotation.txt"))

message("Saved annotations and BED peak sets to ", output_dir)
