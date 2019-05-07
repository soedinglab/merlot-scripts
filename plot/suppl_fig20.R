rm(list = ls()) # clear the environment 
suppressMessages(library(merlot))
library(Matrix)

source("~/Documents/repos/merlot-scripts/plot/get_Monocle2_Paul.R")
source("~/Documents/repos/merlot-scripts/plot/plot_monocle_pie_chart.R")

set.seed(42)
N <- 800
scaf <- CalculateScaffoldTree(CellCoordinates, python_location = "~/miniconda3/envs/py36/bin/python3",
                              reduced=N, BranchMinLengthSensitive = sqrt(N)/2, random_seed=42)

elastic <- CalculateElasticTree(ScaffoldTree=scaf, N_yk=50, NBranchScaffoldNodes = 0)
elastic <- inflate_elastic_tree(elastic, CellCoordinates)

graph_yk <- get_node_graph(elastic)
# make a tree layout
l <- igraph::layout_as_tree(graph_yk, root = 8)
plot_flattened_tree(elastic, cell_annot = broad_type, node_size="topology", legend_position = "topleft", layout=l)

# for more plots we need the embedded tree
FullCoords <- as.matrix(t(exprs(valid_subset_GSE72857_cds2)))
embedded <- GenesSpaceEmbedding(FullCoords, elastic)
pt <- CalculatePseudotimes(embedded, T0=8)
plot_flattened_tree_gene_expression(embedded, "Vwf")
plot_flattened_tree_gene_expression(embedded, "Mpo")
plot_flattened_tree_gene_expression(embedded, "Gata2")

plot_pseudotime_expression_gene(GeneName = "Gata2",
                                EmbeddedTree = embedded,
                                Pseudotimes = pt,
                                addlegend = F, range_y = "tree")
plot_pseudotime_expression_gene(GeneName = "Vwf",
                                EmbeddedTree = embedded,
                                Pseudotimes = pt,
                                addlegend = F, range_y = "tree")
plot_pseudotime_expression_gene(GeneName = "Klf1",
                                EmbeddedTree = embedded,
                                Pseudotimes = pt,
                                addlegend = F, range_y = "cells")
