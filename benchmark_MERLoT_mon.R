#!/usr/bin/Rscript

suppressPackageStartupMessages(library(merlot))
suppressPackageStartupMessages(library(monocle))
suppressPackageStartupMessages(library("optparse"))

LOG_MESSAGE <- ""

option_list = list(
  make_option(c("-t", "--mscripts"), type="character", default="~/Documents/repos/merlot-scripts", 
              help="Location of the package [default= %default]", metavar="/path/to/merlot-scripts"),
  make_option(c("-o", "--out"), type="character", 
              help="Location where the output is stored", metavar="/output/folder"),
  make_option(c("-j", "--job"), type="character", 
              help="job name", metavar="jobname"),
  make_option(c("-d", "--dimensions"), type="integer", default=2,
              help="how many dimensions to use for the embedding [default = %default]", metavar="INT"),
  make_option(c("-k", "--n_yk"), type="integer", default=100,
              help="k for elastic tree [default = %default]", metavar="INT"),
  make_option(c("-m", "--mu0"), type="double", default=0.00625,
              help="mu0 for elastic tree [default = %default]", metavar="double"),
  make_option(c("-l", "--lambda0"), type="double", default=2.03e-09,
              help="lambda0 for elastic tree [default = %default]", metavar="double"),
  make_option(c("-s", "--showparams"), action="store_true", default=FALSE,
              help="Whether to print the values of k, l0, mu0 in the result [default = %default]"),
  make_option(c("-b", "--embedded"), action="store_true", default=FALSE,
              help="Use the embedded instead of the elastic tree. [default = %default]"),
  make_option(c("--sens"), action="store_true", default=FALSE,
              help="Toggle in order to use the number of cells instead of the number of scaffold nodes for branch detection. [default = %default]"),
  make_option(c("-f", "--fixed"), action="store_true", default=FALSE,
              help="Toggle if the number of branches was known while constructing the tree. [default = %default]"),
  make_option(c("-u", "--unconstrained"), action="store_true", default=FALSE,
              help="Toggle to use the unconstrained Monocle2 run (file ending in _unc). [default = %default]")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

mscripts <- opt$mscripts
JobFolder <- opt$out
JobName <- opt$job
dimensions <- opt$dimensions
N_yk <- opt$n_yk
mu0 <- opt$mu0
lambda0 <- opt$lambda0
showparams <- opt$showparams
embed <- opt$embedded
fixed <- opt$fixed
sens <- opt$sens
unc <- opt$unconstrained

if (fixed) {
  sens <- FALSE
}

OgName <- JobName
prefix <- ""

if (unc) {
  prefix <- "unc"
  if (fixed) {
    prefix <- paste(prefix, "fixed", sep = "_")
  } else {
    prefix <- paste(prefix, "auto", sep = "_")
  }
} else {
  if (fixed) {
    prefix <- "fixed"
  } else {
    prefix <- "auto"
  }
}

if (embed) {
  res_prefix <- paste(prefix, "emb", sep = "_")
} else {
  res_prefix <- paste(prefix, "el", sep = "_")
}

various <- paste(mscripts, "/various.R", sep="")
evaluat <- paste(mscripts, "/evaluate_method.R", sep="")

suppressPackageStartupMessages(source(various))
suppressPackageStartupMessages(source(evaluat))

job <- paste(JobFolder, JobName, sep = "")
Dataset <- ReadDataset(paste(job, "simulation.txt", sep = "_"))

where <- paste(JobFolder, JobName, "_monocl2", sep = "")
if (unc) {
  where <- paste(where, "unc", sep = "_")
}
data <- readRDS(where)
CellCoordinates <- t(reducedDimS(data))
LOG_MESSAGE <- paste(LOG_MESSAGE, "dimensions", dim(CellCoordinates), "\n")
assign("LOG_MESSAGE", LOG_MESSAGE, envir = .GlobalEnv)

cell_params <- read.table(file = paste(job, "cellparams.txt", sep = "_"),
                         sep = "\t", header = TRUE, row.names = 1)
time <- cell_params$pseudotime - min(cell_params$pseudotime) + 1
labels <- cell_params$branches + 1
N <- length(labels)

# now calculate tree:
# prepare data for TreeTopology.py
if (embed) {
  methname <- paste("_MERLoT_mon_", prefix, "_el", sep = "")
  if (sens) {
    res_prefix <- paste(res_prefix, "sens", sep = "_")
    methname <- paste(methname, "sens", sep = "_")
  }
  etree <- readRDS(paste(job, methname, ".elastic", sep = ""))
  normed_exp <- 1 / sizeFactors(data) * Dataset$ExpressionMatrix
  tree <- GenesSpaceEmbedding(ExpressionMatrix = normed_exp,
            ElasticTree = etree)
} else {
  if(fixed) {
    scaffold <- CalculateScaffoldTree(CellCoordinates,
                  NEndpoints = dimensions + 1)
  } else if (sens) {
    scaffold <- CalculateScaffoldTree(CellCoordinates,
                  BranchMinLengthSensitive = sqrt(N))
    res_prefix <- paste(res_prefix, "sens", sep = "_")
  } else {
    scaffold <- CalculateScaffoldTree(CellCoordinates)
  }

  tree <- CalculateElasticTree(scaffold, N_yk)
  methname <- paste("_MERLoT_mon_", res_prefix, sep = "")
  saveRDS(object = tree, file = paste(job, methname, ".elastic", sep = ""))
  saveRDS(object = scaffold, file = paste(job, methname, ".scaf", sep = ""))
}


nodes <- tree$Cells2TreeNodes[,2]
Pseudotimes <- CalculatePseudotimes(tree, C0=1)
# C0 = min(which(nodes %in% tree$Topology$Endpoints))
# if (is.infinite(C0)) {
#   Pseudotimes = CalculatePseudotimes(tree, C0=1)
# } else {
#   Pseudotimes = CalculatePseudotimes(tree, C0=C0)
# }

hhbranches <- tree$Cells2Branches
hhtimes <- Pseudotimes$Times_cells
hhcorr <- Pseudotimes$Proyected_Times_Cells

# sometimes if a topology is predicted linear we don't get
# a branch prediction for every cell
if (length(hhbranches) == 1) {
  hhbranches <- rep(1, length(hhtimes))
}

par_loc <- paste(job, "params.txt", sep = "_")
methname <- paste("MERLoT_mon", res_prefix, sep = "_")
where <- paste(JobFolder, OgName, "_MERLoT_mon", sep = "")

res <- evaluate_method(methname, hhbranches, hhtimes, cell_params, par_loc)
# print(hhres)
if (showparams) {
  res <- cbind(res, N_yk, mu0, lambda0)
}

# save coordinates, branch assignments and pseudotime predictions
branches <- hhbranches
times <- hhtimes
where <- paste(JobFolder, OgName, sep = "")
read_write_output(res, paste(where, "_eval.txt", sep = ""))

# methname <- paste("conTree", res_prefix, sep="_")
# res <- evaluate_method(methname, hhbranches, hhcorr, labels, time)
# read_write_output(res, paste(where, "_eval.txt", sep=""))

where <- paste(JobFolder, OgName, "_MERLoT_mon", sep = "")
if (unc) {
  where <- paste(where, "_unc", sep = "")
}
# the coordinates are the monocle2 coordinates
# write.csv(CellCoordinates, paste(where, "coordinates.csv", sep="_"))

# but the branches and pseudotime depend on the parameters, so the filenames should reflect that
where <- paste(JobFolder, OgName, "_", methname, sep = "")
write.csv(branches, paste(where, "branches.csv", sep = "_"))
write.csv(times, paste(where, "times.csv", sep = "_"))
# write.csv(hhcorr, paste(where, "ctimes.csv", sep="_"))
write(x = LOG_MESSAGE, file = paste(where, ".log", sep = ""))
