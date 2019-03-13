#!/usr/bin/Rscript

suppressPackageStartupMessages(library(slingshot))
suppressPackageStartupMessages(library(mclust))
suppressPackageStartupMessages(library("optparse"))

option_list <- list(
  make_option(c("-t", "--mscripts"), type = "character", default = "~/Documents/repos/mscripts", 
              help = "Location of the package [default= %default]", metavar = "/path/to/mscripts"),
  make_option(c("-o", "--out"), type = "character", 
              help = "Location where the output is stored", metavar = "/output/folder"),
  make_option(c("-j", "--job"), type = "character", 
              help = "job name", metavar = "jobname"),
  make_option(c("-d", "--dimensions"), type = "integer", default = 2,
              help = "how many dimensions to use for the embedding [default = %default]", metavar = "INT"),
  make_option(c("-i", "--input"), type = "character", default = "destiny",
              help="The input file suffix. [default = %default]")
);

opt_parser <- OptionParser(option_list = option_list);
opt <- parse_args(opt_parser);

mscripts <- opt$mscripts
JobFolder <- opt$out
JobName <- opt$job
dimensions <- opt$dimensions
filein <- opt$input

# mscripts <- "~/Documents/repos/mscripts"
# JobFolder <- "/data/niko/test/benchmark1/"
# JobName <- "splat15"
# JobFolder <- paste(JobFolder, JobName, "/", sep="")
# filein <- paste(JobFolder, JobName, "_monocl2", sep="")
# dimensions <- 2

job <- paste(JobFolder, JobName, sep = "")

CellCoordinates <- NULL
split_name <- strsplit(filein, JobName)[[1]]
suffix <- split_name[length(split_name)]
input_data <- readRDS(filein)

if (grepl("destiny", filein)) {
  input_coordinates <- input_data@eigenvectors[,1:dimensions]
  dimensions <- min(c(dimensions, dim(input_coordinates)[2]))
  CellCoordinates <- input_coordinates[,1:dimensions]
} else if (grepl("monocl", filein)) {
  input_coordinates <- t(input_data@reducedDimS)[,1:dimensions]
  dimensions <- min(c(dimensions, dim(input_coordinates)[2]))
  CellCoordinates <- input_coordinates[,1:dimensions]
} else if (grepl("TSCAN", filein)) {
  input_coordinates <- input_data$pcareduceres[,1:dimensions]
  dimensions <- min(c(dimensions, dim(input_coordinates)[2]))
  CellCoordinates <- input_coordinates[,1:dimensions]
}

# load functions
various <- paste(mscripts, "/scripts/various.R", sep = "")
evaluat <- paste(mscripts, "/scripts/evaluate_method.R", sep = "")
clust2branch <- paste(mscripts, "/scripts/clusters_to_branches.R", sep = "")
suppressPackageStartupMessages(source(various))
suppressPackageStartupMessages(source(evaluat))
suppressPackageStartupMessages(source(clust2branch))

# read the data
par_loc <- paste(job, "params.txt", sep = "_")
methname <- paste("slingshot", suffix, sep = "")

# read the simulation parameters
cell_params <- read.table(file = paste(job, "cellparams.txt", sep = "_"), sep = "\t", header = T, row.names = 1)
time <- cell_params$pseudotime - min(cell_params$pseudotime) + 1
start <- min(which(time == min(time)))
labels <- cell_params$branches + 1

# find number of clusters and perform clustering with mclust
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
clusterid <- as.numeric(sds@clusterLabels)
sling_branches <- assign_branches(mstree, clusterid, centroids, CellCoordinates)

# predict pseudotime per trajectory for each cell
sds <- getCurves(sds)
traj_pseudotime <- slingPseudotime(sds)
# average over the different trajectories for global pt. This is ok to do, because:
# 1. all trajectories start at the same cell, so we are averaging geodesic distances from the same point
# 2. if cell c doesn't belong to trajectory t it doesn't get a pseudotime for it
# 3. cells on overlapping parts of the trajectory have very similar pt
sling_times <- apply(X = traj_pseudotime, MARGIN = 1, FUN = function(x) mean(na.omit(x)))

where <- paste(JobFolder, JobName, "_", methname, sep = "")
write.csv(sling_branches, paste(where, "branches.csv", sep = "_"))
write.csv(sling_times, paste(where, "times.csv", sep = "_"))

saveRDS(object = sds, file = paste(where, sep = "_"))

res <- evaluate_method(methname, sling_branches, sling_times, cell_params, par_loc)
where <- paste(JobFolder, JobName, sep = "")
read_write_output(res, paste(where, "_eval.txt", sep = ""))


