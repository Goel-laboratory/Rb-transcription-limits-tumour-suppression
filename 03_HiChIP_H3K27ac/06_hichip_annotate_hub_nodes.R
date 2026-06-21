#!/usr/bin/env Rscript

# Annotate hub nodes with genes, promoter status and RB/ER overlaps.
# Optional expression columns are joined by geneId.
# Outputs: annotated node CSV and GRanges RDS.

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- if (length(script_arg)) dirname(normalizePath(sub("^--file=", "", script_arg[[1]]))) else "."
source(file.path(script_dir, "00_hichip_config.R"))

suppressPackageStartupMessages({
  library(ChIPseeker)
  library(dplyr)
  library(GenomicFeatures)
  library(GenomicRanges)
  library(purrr)
  library(readr)
  library(rtracklayer)
  library(stringr)
  library(tibble)
})

membership_path <- hichip_output("HiChIP_hub_node_membership.csv")
if (!file.exists(membership_path)) stop("Run step 05 first; missing ", membership_path)

gtf_path <- hichip_env_file("HICHIP_GENCODE_GTF")
gene_map_path <- hichip_env_file("HICHIP_GENE_MAP")
peak_manifest_path <- hichip_env_file("HICHIP_PEAK_MANIFEST")
expression_path <- hichip_env_file("HICHIP_EXPRESSION_TABLE", required = FALSE)

nodes <- read_csv(membership_path, show_col_types = FALSE)
parts <- str_match(nodes$node_id, "^([^:]+):(\\d+)-(\\d+)$")
if (anyNA(parts[, 2:4])) stop("One or more graph node IDs are not chr:start-end coordinates.")

# Reconstruct genomic ranges from the coordinate strings used as graph nodes.
node_gr <- GRanges(
  seqnames = parts[, 2],
  ranges = IRanges(as.integer(parts[, 3]), as.integer(parts[, 4])),
  node_id = nodes$node_id,
  condition = nodes$condition,
  hub_id = nodes$hub_id
)

# Original annotation: GENCODE v35 basic, promoter -2000/+500, promoter-first
# annotation priority, then simplified genomic classes.
txdb <- GenomicFeatures::makeTxDbFromGFF(gtf_path, format = "gtf")
peak_anno <- ChIPseeker::annotatePeak(
  node_gr,
  tssRegion = c(-2000, 500),
  TxDb = txdb,
  genomicAnnotationPriority = c("Promoter", "5UTR", "3UTR", "Exon", "Intron", "Intergenic")
)
annotated <- as.data.frame(peak_anno)

annotate_simple <- function(x) {
  ifelse(
    grepl("Exon", x), "Exon",
    ifelse(
      grepl("Intron", x), "Intron",
      ifelse(
        grepl("Promoter", x), "Promoter",
        ifelse(grepl("Downstream|Distal Intergenic", x), "Intergenic", x)
      )
    )
  )
}

annotated$annotation_simple <- annotate_simple(annotated$annotation)
annotated$annotation_expand <- with(
  annotated,
  ifelse(
    annotation_simple == "Intergenic" & distanceToTSS > -10000 & distanceToTSS < 10000,
    "Intergenic_Proximal",
    ifelse(
      annotation_simple == "Intergenic" & abs(distanceToTSS) >= 10000 & abs(distanceToTSS) < 100000,
      "Intergenic_Distal",
      ifelse(annotation_simple == "Intergenic" & abs(distanceToTSS) >= 100000, "Intergenic_Desert", annotation_simple)
    )
  )
)
annotated$annotation_promoter <- ifelse(annotated$annotation_simple == "Promoter", "Promoter", "Other")

# Add public gene symbols and biotypes in place of the private biomart object.
gene_map <- read_csv(gene_map_path, show_col_types = FALSE)
hichip_require_columns(gene_map, c("geneId", "SYMBOL", "gene_biotype"), "HICHIP_GENE_MAP")
annotated <- annotated %>%
  left_join(gene_map %>% distinct(geneId, .keep_all = TRUE), by = "geneId")

peak_manifest <- read_csv(peak_manifest_path, show_col_types = FALSE)
hichip_require_columns(peak_manifest, c("condition", "factor", "class", "file"), "HICHIP_PEAK_MANIFEST")
if (any(!file.exists(peak_manifest$file))) stop("HICHIP_PEAK_MANIFEST contains missing BED files.")
if (any(!peak_manifest$class %in% c("promoter", "nonpromoter"))) {
  stop("Peak-manifest class values must be promoter or nonpromoter.")
}
peak_sets <- lapply(peak_manifest$file, rtracklayer::import)

node_gr_annotated <- makeGRangesFromDataFrame(
  annotated,
  keep.extra.columns = TRUE,
  seqnames.field = "seqnames",
  start.field = "start",
  end.field = "end"
)

# Apply condition-matched peak sets; ALL denotes a shared differential set.
label_factor <- function(gr, condition, factor_name) {
  rows <- which(
    peak_manifest$factor == factor_name &
      peak_manifest$condition %in% c(condition, "ALL")
  )
  promoter_rows <- rows[peak_manifest$class[rows] == "promoter"]
  nonpromoter_rows <- rows[peak_manifest$class[rows] == "nonpromoter"]
  promoter <- if (length(promoter_rows)) do.call(c, peak_sets[promoter_rows]) else GRanges()
  nonpromoter <- if (length(nonpromoter_rows)) do.call(c, peak_sets[nonpromoter_rows]) else GRanges()
  prefix <- if (factor_name == "RB_up") "RB" else factor_name
  ifelse(
    overlapsAny(gr, promoter), paste0(prefix, "prom"),
    ifelse(overlapsAny(gr, nonpromoter), paste0(prefix, "nonprom"), "none")
  )
}

# Label RB, differential RB-up and ER occupancy for every condition.
for (condition in unique(node_gr_annotated$condition)) {
  idx <- node_gr_annotated$condition == condition
  current <- node_gr_annotated[idx]
  node_gr_annotated$RB_overlap[idx] <- label_factor(current, condition, "RB")
  node_gr_annotated$RB_up_overlap[idx] <- label_factor(current, condition, "RB_up")
  node_gr_annotated$ER_overlap[idx] <- label_factor(current, condition, "ER")
}

annotated <- as.data.frame(node_gr_annotated)
# Expression input is optional because hub construction does not require it.
if (nzchar(expression_path)) {
  expression <- read_csv(expression_path, show_col_types = FALSE)
  hichip_require_columns(expression, "geneId", "HICHIP_EXPRESSION_TABLE")
  annotated <- annotated %>% left_join(expression, by = "geneId", suffix = c("", "_expression"))
}

write_csv(annotated, hichip_output("HiChIP_hub_nodes_annotated.csv"))
saveRDS(node_gr_annotated, hichip_output("HiChIP_hub_nodes_annotated.rds"))
capture.output(sessionInfo(), file = hichip_output("sessionInfo_06_hub_annotation.txt"))
