library(igraph)
library(mclust, quietly = TRUE)
library(slingshot)

various <- paste(mscripts, "/scripts/various.R", sep = "")
evaluat <- paste(mscripts, "/scripts/evaluate_method.R", sep = "")

suppressPackageStartupMessages(source("~/Documents/repos/merlot-scripts/scripts/evaluate_method.R"))
suppressPackageStartupMessages(source("~/Documents/repos/merlot-scripts/scripts/various.R"))
source("~/Documents/repos/merlot-scripts/scripts/clusters_to_branches.R")

benchmark <- "benchmark5"
sim <- "sim68"
job <- paste("~/Documents/data/examine", benchmark, sim, sim, sep="/")

dm <- readRDS(paste(job, "destiny_log_k", sep="_"))
CellCoordinates <- dm@eigenvectors[, 1:6]
cell_params <- read.table(file = paste(job, "cellparams.txt", sep = "_"), sep = "\t", header = T, row.names = 1)
time <- cell_params$pseudotime - min(cell_params$pseudotime) + 1
start <- min(which(time == min(time)))
labels <- cell_params$branches + 1
par_loc <- paste(job, "params.txt", sep = "_")
N <- dim(CellCoordinates)[1]
K <- dim(CellCoordinates)[2]

mclusters <- mclustBIC(CellCoordinates, G = 5:50)
mod1 <- Mclust(CellCoordinates, x = mclusters)

# use clusters to get lineages
# specify the cluster of cell 1 as start to get best pseudotime results
sds <- tryCatch( {
  getLineages(CellCoordinates, mod1$classification, start.clus = mod1$classification[start])
  }, 
  error = function(cond) {
    print("System is computationally singular. Rerunning with some random noise.")
    print(cond)
    # LOG_MESSAGE <- paste(LOG_MESSAGE, cond, "\n")
    # assign("LOG_MESSAGE", LOG_MESSAGE, envir = .GlobalEnv)
    span <- range(CellCoordinates)[2] - range(CellCoordinates)[1]
    noise <- rnorm(N*K, mean = 0, sd = span/100)
    dim(noise) <- c(N, K)
    getLineages(CellCoordinates+noise, mod1$classification, start.clus = mod1$classification[start])
  })

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
adjusted_mi(labels, sling_branches, "bla")

# predict pseudotime per trajectory for each cell
sds <- getCurves(sds)
sling <- sds

# which_sling <- "slingshot_destiny_log_k"
# sling <- readRDS(paste(job, which_sling, sep="_"))
# interpreted <- read.table(paste(job, which_sling, "branches.csv", sep="_"), sep=",", header = T)

predicted_labels <- as.factor(apply(sling@clusterLabels, 1, function(x) names(which(x == 1))))

par(mfrow=c(1,2))
plot(reducedDims(sling), col = predicted_labels, pch=16, asp = 1)
lines(sling, lwd=2, type = 'lineages', col="black")

plot(reducedDims(sling), col = interpreted$x, pch=16, asp = 1)
lines(sling, lwd=2, type = 'lineages', col="black")
par(mfrow=c(1,1))
