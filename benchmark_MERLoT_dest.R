#!/usr/bin/Rscript

suppressPackageStartupMessages(library(merlot))
suppressPackageStartupMessages(library(destiny))
library("optparse")

LOG_MESSAGE <- ""

option_list <- list(
  make_option(c("-t", "--mscripts"), type = "character", default = "~/Documents/repos/merlot-scripts", 
              help = "Location of the package [default= %default]", metavar = "/path/to/merlot-scripts"),
  make_option(c("-o", "--out"), type = "character", 
              help = "Location where the output is stored", metavar = "/output/folder"),
  make_option(c("-j", "--job"), type = "character", 
              help = "job name", metavar = "jobname"),
  make_option(c("-d", "--dimensions"), type = "integer", default = 2,
              help = "how many dimensions to use for the embedding [default = %default]", metavar = "INT"),
  make_option(c("-k", "--n_yk"), type = "integer", default = 100,
              help = "k for elastic tree [default = %default]", metavar = "INT"),
  make_option(c("-m", "--mu0"), type = "double", default = 0.00625,
              help = "mu0 for elastic tree [default = %default]", metavar = "double"),
  make_option(c("-l", "--lambda0"), type = "double", default = 2.03e-09,
              help = "lambda0 for elastic tree [default = %default]", metavar = "double"),
  make_option(c("-s", "--showparams"), action = "store_true", default = FALSE,
              help = "Whether to print the values of k, l0, mu0 in the result [default = %default]"),
  make_option(c("-g", "--grid"), action = "store_true", default = FALSE,
              help = "Toggle to change the output file name when performing a grid search [default = %default]"),
  make_option(c("-b", "--embedded"), action = "store_true", default = FALSE,
              help = "Toggle to use the embedded instead of the elastic tree. [default = %default]"),
  make_option(c("-f", "--fixed"), action = "store_true", default = FALSE,
              help = "Toggle if the number of branches was known while constructing the tree. [default = %default]"),
  make_option(c("--sens"), action = "store_true", default = FALSE,
              help = "Toggle in order to use the number of cells instead of the number of scaffold nodes for branch detection. [default = %default]"),
  make_option(c("-n", "--knn"), action = "store_true", default = FALSE,
              help = "Toggle to use a simple heuristic for choosing a k for diffusion maps. [default = %default]"),
  make_option("--log", action = "store_true", default = FALSE,
              help = "Toggle to transform data to log(data+1) before analysis. [default = %default]"),
  make_option(c("-i", "--interpolate"), type = "character", action = "store", default = "none",
              help = "Whether and where to add interpolated nodes to the tree before evaluation. One of 'elastic', 'emb', 'elemb', 'none'. [default = %default]")
); 

opt_parser <- OptionParser(option_list = option_list);
opt <- parse_args(opt_parser);

mscripts <- opt$mscripts
JobFolder <- opt$out
JobName <- opt$job
dimensions <- opt$dimensions
N_yk <- opt$n_yk
mu0 <- opt$mu0
lambda0 <- opt$lambda0
showparams <- opt$showparams
isgrid <- opt$grid
embed <- opt$embedded
fixed <- opt$fixed
sens <- opt$sens
knn <- opt$knn
iflog <- opt$log
interp <- opt$interpolate

choices <- c("elastic", "emb", "elemb", "none")
if(!interp %in% choices) {
  sprintf("Interpolation option %s not set correctly, setting to 'none'", interp)
  interp <- "none"
}

if (fixed) {
  sens <- FALSE
}

# mscripts <- "~/Documents/repos/mscripts"
# JobFolder <- "/data/niko/final/benchmark2/test0/"
# JobName <- "test0"
# dimensions <- 2
# N_yk <- 100
# mu0 <- 0.00625
# lambda0 <- 2.03e-09
# showparams <- FALSE
# isgrid <- FALSE
# embed <- TRUE
# fixed <- FALSE
# sens <- TRUE
# knn <- FALSE
# iflog <- TRUE
# interp <- "none"


if (interp == "elastic" && embed == TRUE) {
  stop("Interpolation and evaluation of elastic tree not compatible with \\
  '--embed TRUE' option.")
  warning("Interpolation and evaluation of elastic tree not compatible with \\
  '--embed TRUE' option.")
}

OgName <- JobName

job <- paste(JobFolder, JobName, sep = "")

prefix <- ""
if (iflog) prefix <- paste(prefix, "_log", sep = "")
if (knn) prefix <- paste(prefix, "_k", sep = "")
diffmap <- paste("destiny", prefix, sep = "")
# print(paste("using ", job, diffmap, sep = "_"))
dif <- readRDS(file = paste(job, diffmap, sep = "_"))
CellCoordinates <- dif@eigenvectors[,1:dimensions]
LOG_MESSAGE <- paste(LOG_MESSAGE, "dimensions:", dim(CellCoordinates), "\n")
assign("LOG_MESSAGE", LOG_MESSAGE, envir = .GlobalEnv)
rm(dif)

# if (knn) {
#   prefix = paste(prefix, "k", sep = "")
#   diffmap <- paste(JobFolder, JobName, "_destiny_", prefix, sep = "")
#   CellCoordinates = read.table(diffmap)
#   prefix = paste(prefix, "_", sep = "")
# } else {
#   diffmap <- paste(JobFolder, JobName, "_destiny", sep = "")
#   CellCoordinates = read.table(diffmap)
# }

if (fixed) {
  prefix <- paste(prefix, "fixed", sep = "_")
} else {
  prefix <- paste(prefix, "auto", sep = "_")
}

if (embed) {
  res_prefix <- paste(prefix, "emb", sep = "_")
} else {
  res_prefix <- paste(prefix, "el", sep = "_")
}
# print(prefix)

various <- paste(mscripts, "/various.R", sep = "")
evaluat <- paste(mscripts, "/evaluate_method.R", sep = "")

suppressPackageStartupMessages(source(various))
suppressPackageStartupMessages(source(evaluat))

Dataset <- ReadDataset(paste(job, "simulation.txt", sep = "_"))

sum1 <- apply(Dataset$ExpressionMatrix, 1, sum)
scalings <- sum1/mean(sum1)
Dataset$ExpressionMatrix = (1/scalings)*Dataset$ExpressionMatrix
if (iflog) Dataset$ExpressionMatrix = log(Dataset$ExpressionMatrix+1)

cells <- rownames(raw)
genes <- colnames(raw)

cell_params <- read.table(file = paste(job, "cellparams.txt", sep = "_"), sep = "\t", header = T, row.names = 1)
time <- cell_params$pseudotime - min(cell_params$pseudotime) + 1
labels <- cell_params$branches + 1

if(embed) {
  methname <- paste("_MERLoT_dest", prefix, sep = "")

  if (interp == "elemb") {
    methname <- paste(methname, "_el_il", sep = "")
    res_prefix <- paste(res_prefix, "ilm", sep = "_")
  } else {
    methname <- paste(methname, "_el", sep = "")
  }

  if (sens) {
    res_prefix <- paste(res_prefix, "sens", sep = "_")
    methname <- paste(methname, "sens", sep = "_")
  }
  etree <- readRDS(paste(job, methname, ".elastic", sep = ""))
  tree <- GenesSpaceEmbedding(ExpressionMatrix = Dataset$ExpressionMatrix, ElasticTree = etree)

  if (interp == "emb") {
    tree <- DuplicateTreeNodes(tree)
    methname <- paste(methname, "_emb_im", sep = "")
    res_prefix <- paste(res_prefix, "im", sep = "_")
  }

} else {
  methname <- paste("_MERLoT_dest", res_prefix, sep = "")
  # since we run many versions of the elastic tree maybe we have already calculated the scaffold tree
  if (interp == "elastic") {
    scaffold <- readRDS(paste(job, methname, ".scaf", sep = ""))
    tree <- readRDS(paste(job, methname, ".elastic", sep = ""))

    tree <- DuplicateTreeNodes(tree)
    res_prefix <- paste(res_prefix, "il", sep = "_")
    saveRDS(object = tree, file = paste(job, methname, "_il.elastic", sep = ""))
  } else {

    if(fixed) {
      scaffold <- CalculateScaffoldTree(CellCoordinates, NEndpoints = dimensions + 1)
    } else {
      if (sens) {
        N <- length(cells)
        res_prefix <- paste(res_prefix, "sens", sep = "_")
        scaffold <- CalculateScaffoldTree(CellCoordinates, BranchMinLengthSensitive = sqrt(N))
      } else {
        scaffold <- CalculateScaffoldTree(CellCoordinates)
      }
    }
    tree <- CalculateElasticTree(scaffold, N_yk, FixEndpoints = F)
    methname <- paste("_MERLoT_dest", res_prefix, sep = "")
    saveRDS(object = tree, file = paste(job, methname, ".elastic", sep = ""))
    saveRDS(object = scaffold, file = paste(job, methname, ".scaf", sep = ""))
  }
}

nodes = tree$Cells2TreeNodes[, 2]
Pseudotimes <- CalculatePseudotimes(tree, C0 = 1)
# C0 = min(which(nodes %in% tree$Topology$Endpoints))
# if(is.infinite(C0)) {
#   Pseudotimes = CalculatePseudotimes(tree, C0 = 1)
# } else {
#   Pseudotimes = CalculatePseudotimes(tree, C0 = C0)
# }


hhbranches <- tree$Cells2Branches
hhtimes <- Pseudotimes$Times_cells
hhcorr <- Pseudotimes$Proyected_Times_Cells

# sometimes if a topology is predicted linear we don't get
# a branch prediction for every cell
if(length(hhbranches) == 1) {
  hhbranches <- replicate(length(hhtimes), 1)
  hhcorr <- hhtimes
}


# print(hhtimes)
par_loc <- paste(job, "params.txt", sep = "_")
methname <- paste("MERLoT_dest", res_prefix, sep = "")
res <- evaluate_method(methname, hhbranches, hhtimes, cell_params, par_loc)
# print(hhres)
if (showparams) {
  res <- cbind(res, N_yk, mu0, lambda0)
}

# save coordinates, branch assignments and pseudotime predictions
# coords <- CellCoordinates
branches <- hhbranches
times <- hhtimes
where <- paste(JobFolder, OgName, sep = "")
read_write_output(res, paste(where, "_eval.txt", sep = ""))

# methname <- paste("chTree", res_prefix, sep = "_")
# res = evaluate_method(methname, hhbranches, hhcorr, labels, time)
# read_write_output(res, paste(where, "_eval.txt", sep = ""))

if (isgrid) {
  where <- paste(where, N_yk, mu0, lambda0, sep = "_")
}

# mark everything with exact method name in order to help evaluation
where <- paste(JobFolder, OgName, "_", methname, sep = "")
write.csv(branches, paste(where, "branches.csv", sep = "_"))
write.csv(times, paste(where, "times.csv", sep = "_"))
# write.csv(hhcorr, paste(where, "ctimes.csv", sep = "_"))
write(x = LOG_MESSAGE, file = paste(where, ".log", sep = ""))