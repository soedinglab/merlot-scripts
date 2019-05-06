
# Plotting of Monocle tree
cell_annot <- broad_type
mondata <- valid_subset_GSE72857_cds2

graph_yk <- mondata@minSpanningTree
cells2nodes <- mondata@auxOrderingData[["DDRTree"]]$pr_graph_cell_proj_closest_vertex
cell_annot=broad_type
layout=NULL
legend_position="topleft"
node_size="topology"
node_radius=12
cols=NULL

NumberOfNodes <- vcount(graph_yk)
N <- length(which(degree(graph_yk) == 1))
B <- length(which(degree(graph_yk) > 2))

if (is.null(cols)) {
  selected_colors <- c("forestgreen", "firebrick3", "dodgerblue3", "darkorchid",
                       "darkorange3", "orange", "blue", "aquamarine", "magenta",
                       "brown", "gray", "wheat1", "azure4", "lightsalmon4",
                       "navy", "sienna1", "gold4", "red4", "violetred", "black")
} else {
  selected_colors <- cols
}

# first decide what plot we are making
# if cell_annot is empty, we are making the normal plot
# if cell_annot has values we are making pie charts
# if we are making pie charts, their size is either uniform
# or it depends on how many cells belong to each pie
annot_levels <- levels(as.factor(cell_annot))
if (length(annot_levels) <= length(selected_colors)) {
  pie_cols <- c("white", selected_colors[1:length(annot_levels)])
} else {
  print(paste("Supplied colors not enough; creating LSD palette with", length(annot_levels), "colors"))
  tmp <- LSD::distinctcolors(length(annot_levels), show = FALSE, bw = TRUE)
  pie_cols <- c("white", tmp)
}
annot_levels <- c("empty", annot_levels)
pie_values <- lapply(1:NumberOfNodes, function(x) table(annot_levels) - 2)
pie_cols <- pie_cols[order(annot_levels)]
annot_levels <- annot_levels[order(annot_levels)]

for (n in 1:NumberOfNodes) {
  cells_n <- (cells2nodes[,1] == n)
  for (real in annot_levels) {
    pie_values[[n]][real] <- sum(cell_annot[cells_n] == real)
  }
  if (sum(pie_values[[n]]) == 0) {
    pie_values[[n]]["empty"] <- 1
  }
}

nodes_labels=c(seq(1, N + B), rep(NA, NumberOfNodes - N - B))

avg_cells_per_node <- length(cell_annot) / NumberOfNodes
nodes_sizes = rep(5, NumberOfNodes)
cells <- rep(0, NumberOfNodes)
for (n in 1:NumberOfNodes) {
  cells[n] <- sum((cells2nodes[,1] == n))
  nodes_sizes[n] <- nodes_sizes[n] * (cells[n] / avg_cells_per_node)
}
sizes <- seq(0, max(cells), 10)[-1]
size_length <- length(sizes)
keep <- c(1,
          as.integer(size_length/4),
          as.integer(size_length/2),
          as.integer(size_length) * 0.75,
          size_length)

l <- igraph::layout_as_tree(graph_yk, root = 79)
plot(graph_yk, layout=l,
     vertex.label=nodes_labels,
     vertex.size=nodes_sizes,
     vertex.shape="pie",
     vertex.pie=pie_values,
     vertex.pie.color=list(pie_cols),
     vertex.label.cex=0.6)
l <- legend(x=legend_position,
            legend=annot_levels,
            col="black",
            pt.bg=pie_cols,
            pch=21)
a <- legend(x = l$rect$left, y=l$rect$top - l$rect$h,
            legend=sizes[keep],
            pt.cex= keep,
            col='black')
x <- (a$text$x + a$rect$left) / 2
y <- a$text$y
graphics::symbols(x,y,circles=sizes[keep] * 5 / (200 * avg_cells_per_node) ,inches=FALSE,add=TRUE,bg='orange')

