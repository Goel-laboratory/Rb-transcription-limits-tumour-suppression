#!/usr/bin/env Rscript

# Summarise each hub and optionally run Hallmark-style over-representation analysis.
# Inputs are annotated nodes, graph edge counts and optional TERM2GENE.
# Outputs include hub summaries, category counts and enrichment tables.

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- if (length(script_arg)) dirname(normalizePath(sub("^--file=", "", script_arg[[1]]))) else "."
source(file.path(script_dir, "00_hichip_config.R"))

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(dplyr)
  library(purrr)
  library(readr)
  library(stringr)
  library(tidyr)
})

node_path <- hichip_output("HiChIP_hub_nodes_annotated.csv")
edge_summary_path <- hichip_output("HiChIP_hub_graph_summary.csv")
if (!file.exists(node_path) || !file.exists(edge_summary_path)) {
  stop("Run steps 05 and 06 before hub summarisation.")
}

nodes <- read_csv(node_path, show_col_types = FALSE)
edge_summary <- read_csv(edge_summary_path, show_col_types = FALSE)

# Summarise RB/ER occupancy and protein-coding promoter genes per hub.
hub_summary <- nodes %>%
  group_by(condition, hub_id) %>%
  summarise(
    number_nodes = n(),
    n_edges_legacy_source_field = n(),
    total_RBprom = sum(RB_overlap == "RBprom", na.rm = TRUE),
    total_RBnonprom = sum(RB_overlap == "RBnonprom", na.rm = TRUE),
    total_RBprom_up = sum(RB_up_overlap == "RBprom", na.rm = TRUE),
    total_RBnonprom_up = sum(RB_up_overlap == "RBnonprom", na.rm = TRUE),
    total_ERprom = sum(ER_overlap == "ERprom", na.rm = TRUE),
    total_ERnonprom = sum(ER_overlap == "ERnonprom", na.rm = TRUE),
    promoter_protein_coding_genes = paste(
      unique(na.omit(SYMBOL[annotation_promoter == "Promoter" & gene_biotype == "protein_coding"])),
      collapse = "|"
    ),
    .groups = "drop"
  ) %>%
  left_join(
    edge_summary %>% select(condition, hub_id, number_intra_hub_edges),
    by = c("condition", "hub_id")
  )

write_csv(hub_summary, hichip_output("HiChIP_hub_biological_summary.csv"))

# Reproduce the key RB-associated and RB/ER-associated hub counts.
hub_counts <- hub_summary %>%
  group_by(condition) %>%
  summarise(
    total_hubs = n(),
    hubs_RBnonprom = sum(total_RBnonprom > 0),
    hubs_RBnonprom_up = sum(total_RBnonprom_up > 0),
    hubs_RBnonprom_and_ER = sum(total_RBnonprom > 0 & (total_ERprom > 0 | total_ERnonprom > 0)),
    hubs_RBnonprom_up_and_ER = sum(total_RBnonprom_up > 0 & (total_ERprom > 0 | total_ERnonprom > 0)),
    .groups = "drop"
  )
write_csv(hub_counts, hichip_output("HiChIP_hub_category_counts.csv"))

gene_sets_path <- hichip_env_file("HICHIP_TERM2GENE", required = FALSE)
# Enrichment is optional so hub summaries can be produced without MSigDB data.
if (nzchar(gene_sets_path)) {
  term2gene <- read_csv(gene_sets_path, show_col_types = FALSE)
  if (ncol(term2gene) < 2) stop("HICHIP_TERM2GENE must contain at least two columns.")
  term2gene <- term2gene[, 1:2]
  names(term2gene) <- c("term", "gene")

  extract_genes <- function(x) {
    genes <- str_trim(unlist(str_split(replace_na(x, ""), "\\|")))
    unique(genes[!is.na(genes) & nzchar(genes)])
  }

  gene_groups <- hub_summary %>%
    group_by(condition) %>%
    group_split() %>%
    set_names(map_chr(., ~ unique(.x$condition))) %>%
    map(function(x) {
      list(
        all_hub_genes = extract_genes(x$promoter_protein_coding_genes),
        RBnonprom_up_hub_genes = extract_genes(x$promoter_protein_coding_genes[x$total_RBnonprom_up > 0]),
        no_RBnonprom_up_hub_genes = extract_genes(x$promoter_protein_coding_genes[x$total_RBnonprom_up == 0]),
        RBnonprom_ER_hub_genes = extract_genes(
          x$promoter_protein_coding_genes[x$total_RBnonprom > 0 & (x$total_ERprom > 0 | x$total_ERnonprom > 0)]
        ),
        # This preserves the original source expression, which used OR here.
        RBnonprom_noER_source_rule_genes = extract_genes(
          x$promoter_protein_coding_genes[x$total_RBnonprom > 0 & (x$total_ERprom == 0 | x$total_ERnonprom == 0)]
        )
      )
    })

# Preserve source cutoffs and export all tested terms for transparent filtering.
  enrichment <- imap(gene_groups, function(groups, condition) {
    imap(groups, function(genes, group) {
      set.seed(1)
      result <- clusterProfiler::enricher(
        genes,
        TERM2GENE = term2gene,
        pvalueCutoff = 1,
        qvalueCutoff = 1
      )
      if (is.null(result)) return(data.frame())
      as.data.frame(result) %>% mutate(condition = condition, gene_group = group, .before = 1)
    })
  }) %>%
    flatten() %>%
    bind_rows()

  write_csv(enrichment, hichip_output("HiChIP_hub_enrichment_results.csv"))
  saveRDS(gene_groups, hichip_output("HiChIP_hub_gene_groups.rds"))
}

capture.output(sessionInfo(), file = hichip_output("sessionInfo_07_hub_summary_enrichment.txt"))
