#!/usr/bin/env Rscript

# Convert HiChIP loops to condition-specific graphs and identify Louvain hubs.
# Nodes are unique anchor coordinates; loop rows are undirected, equal-weight
# edges. Outputs include graphs, memberships and connectivity summaries.

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- if (length(script_arg)) dirname(normalizePath(sub("^--file=", "", script_arg[[1]]))) else "."
source(file.path(script_dir, "00_hichip_config.R"))

suppressPackageStartupMessages({
  library(dplyr)
  library(igraph)
  library(purrr)
  library(readr)
  library(tibble)
})

input_path <- hichip_output("HC_integrated_overlap_df.rds")
if (!file.exists(input_path)) {
  input_path <- hichip_output("HC_RB_overlap_df.rds")
}
if (!file.exists(input_path)) stop("Run steps 01-03 first; no annotated loop table was found.")

loops <- readRDS(input_path)
seed <- as.integer(Sys.getenv("HICHIP_LOUVAIN_SEED", "1"))
if (is.na(seed)) stop("HICHIP_LOUVAIN_SEED must be an integer.")

# The first two columns are deliberately from/to for graph_from_data_frame().
make_edges <- function(x) {
  hichip_require_columns(x, c("chr1", "s1", "e1", "chr2", "s2", "e2"), "HiChIP loop table")
  x %>%
    mutate(
      from = paste0(chr1, ":", s1, "-", e1),
      to = paste0(chr2, ":", s2, "-", e2)
    ) %>%
    select(from, to, any_of("cc"), everything())
}

edge_tables <- lapply(loops, make_edges)

# This deliberately reproduces the source call. The graph is undirected and is
# not simplified: parallel edges and self-loops, if present, remain in place.
graphs <- lapply(
  edge_tables,
  function(x) igraph::graph_from_data_frame(d = x, directed = FALSE)
)

# The source used lapply(igraph_l, cluster_louvain) without a seed. igraph >=1.3
# processes vertices in random order, so the reviewed workflow sets one seed
# immediately before the otherwise identical condition-wise calls.
set.seed(seed)
communities <- lapply(graphs, igraph::cluster_louvain)
memberships <- lapply(communities, igraph::membership)

# Export a simple node-to-hub map for downstream genomic annotation.
node_tables <- imap_dfr(memberships, function(member, condition) {
  tibble(
    condition = condition,
    node_id = names(member),
    hub_id = as.integer(unname(member))
  )
})

# Count nodes and true intra-community graph edges for each hub.
summarise_graph <- function(graph, community, condition) {
  member <- igraph::membership(community)
  edges <- igraph::as_data_frame(graph, what = "edges")
  edges$hub_from <- unname(member[edges$from])
  edges$hub_to <- unname(member[edges$to])
  intra <- edges[edges$hub_from == edges$hub_to, , drop = FALSE]

  node_summary <- tibble(
    hub_id = as.integer(member),
    node_id = names(member)
  ) %>%
    count(hub_id, name = "number_nodes")

  edge_summary <- as.data.frame(table(intra$hub_from), stringsAsFactors = FALSE)
  names(edge_summary) <- c("hub_id", "number_intra_hub_edges")
  edge_summary$hub_id <- as.integer(edge_summary$hub_id)
  edge_summary$number_intra_hub_edges <- as.integer(edge_summary$number_intra_hub_edges)

  left_join(node_summary, edge_summary, by = "hub_id") %>%
    mutate(
      number_intra_hub_edges = coalesce(number_intra_hub_edges, 0L),
      condition = condition,
      .before = 1
    )
}

hub_graph_summary <- imap_dfr(
  seq_along(graphs),
  function(i, condition) summarise_graph(graphs[[i]], communities[[i]], names(graphs)[[i]])
)

# Record graph properties that were implicit in the historical source.
graph_diagnostics <- imap_dfr(graphs, function(graph, condition) {
  tibble(
    condition = condition,
    number_nodes = igraph::vcount(graph),
    number_edges = igraph::ecount(graph),
    number_communities = length(communities[[condition]]),
    modularity = igraph::modularity(communities[[condition]]),
    directed = igraph::is_directed(graph),
    simple_graph = igraph::is_simple(graph),
    self_loops = sum(igraph::which_loop(graph)),
    edges_in_parallel_sets = sum(igraph::which_multiple(graph)),
    weighted_louvain = FALSE,
    louvain_seed = seed
  )
})

saveRDS(edge_tables, hichip_output("HiChIP_edge_tables.rds"))
saveRDS(graphs, hichip_output("HiChIP_igraph_objects.rds"))
saveRDS(communities, hichip_output("HiChIP_Louvain_communities.rds"))
saveRDS(memberships, hichip_output("HiChIP_Louvain_memberships.rds"))
write_csv(node_tables, hichip_output("HiChIP_hub_node_membership.csv"))
write_csv(hub_graph_summary, hichip_output("HiChIP_hub_graph_summary.csv"))
write_csv(graph_diagnostics, hichip_output("HiChIP_graph_diagnostics.csv"))
capture.output(sessionInfo(), file = hichip_output("sessionInfo_05_Louvain_hubs.txt"))
