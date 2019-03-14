library(igraph)
source("~/Documents/repos/merlot-scripts/scripts/clusters_to_branches.R")

benchmark <- "benchmark5"
sim <- "sim2"
job <- paste("~/Documents/data/examine", benchmark, sim, sim, sep="/")

dm <- readRDS(paste(job, "destiny_log_k", sep="_"))
coords <- dm@eigenvectors[, 1:6]
cell_params <- read.table(file = paste(job, "cellparams.txt", sep = "_"), sep = "\t", header = T, row.names = 1)
time <- cell_params$pseudotime - min(cell_params$pseudotime) + 1
start <- min(which(time == min(time)))
labels <- cell_params$branches + 1

mclusters <- mclustBIC(CellCoordinates, G = 5:50)
mod1 <- Mclust(CellCoordinates, x = mclusters)

# use clusters to get lineages
# specify the cluster of cell 1 as start to get best pseudotime results
sds <- getLineages(CellCoordinates, mod1$classification, start.clus = mod1$classification[start])

# get connectivity of clusters and re-create the MST
nodes_order <- order(as.numeric(colnames(sds@adjacency)))
adj_matrix <- sds@adjacency[nodes_order, nodes_order]
mstree <- graph_from_adjacency_matrix(adj_matrix, mode = "undirected")

# map cells to clusters, get cluster centroids and collapse co-linear clusters
# to form branch assignments for each cell
cells2nodes <- mod1$classification
centroids <- t(mod1$parameters$mean)
clusterid <- as.numeric(cells2nodes)
sling_branches <- assign_branches(mstree, clusterid, centroids, CellCoordinates)

# predict pseudotime per trajectory for each cell
sds <- getCurves(sds)


which_sling <- "slingshot_destiny_log_k"
sling <- readRDS(paste(job, which_sling, sep="_"))
interpreted <- read.table(paste(job, which_sling, "branches.csv", sep="_"), sep=",", header = T)

predicted_labels <- as.factor(apply(sling@clusterLabels, 1, function(x) names(which(x == 1))))

par(mfrow=c(1,2))
plot(reducedDims(sling), col = predicted_labels, pch=16, asp = 1)
lines(sling, lwd=2, type = 'lineages', col="black")

plot(reducedDims(sling), col = interpreted$x, pch=16, asp = 1)
lines(sling, lwd=2, type = 'lineages', col="black")
par(mfrow=c(1,1))
