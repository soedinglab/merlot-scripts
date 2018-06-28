#!/usr/bin/Rscript

suppressPackageStartupMessages(library(splatter))
library(magrittr)
library("optparse")

option_list <- list(
  make_option(c("-t", "--mscripts"), type="character", default="~/Documents/repos/mscripts", 
              help="Location of the package [default= %default]", metavar="/path/to/mscripts"),
  make_option(c("-o", "--out"), type = "character", 
              help = "Location where the output is stored", metavar = "/output/folder"),
  make_option(c("-j", "--job"), type = "character", 
              help = "Job ID (prepended to all generated files)", metavar = "jobname"),
  make_option(c("-n", "--num_brpoints"), type = "numeric",
              help = "How many branching points the simulation contains.")
); 

opt_parser <- OptionParser(option_list = option_list);
opt <- parse_args(opt_parser);

mscripts <- opt$mscripts
JobFolder <- opt$out
JobName <- opt$job
num_bif <- opt$num_brpoints

splat_funcs <- paste(mscripts, "/scripts/splat_funcs.R", sep="")
suppressPackageStartupMessages(source(splat_funcs))

# JobFolder <- "/data/niko/test/splat0/"
# JobName <- "splat0"
# num_bif <- 2

# set random seed for reproducibility
seed <- sample(0:999, 1)
set.seed(seed)
# sample number of genes
G <- sample(1000:10000, 1)
# set branch length
branch_length <- 50

num_branches <- 2 * num_bif + 1
group_probs <- rep(1/num_branches, num_branches)
total_cells <- branch_length * num_branches
lineage_tree <- gen_random_topology(num_bif)
  
params <- newSplatParams(seed=seed,
                         nGenes=G,
                         batchCells=total_cells,
                         path.length=branch_length,
                         group.prob = group_probs,
                         path.from=lineage_tree)
sim.paths <- splatSimulate(method = "paths", verbose = FALSE, params=params)

# count matrix
X <- t(counts(sim.paths))

# branch identity of each cell
tmp <- as.factor(sim.paths@colData$Group)
branches <- unclass(tmp) %>% as.numeric

# library size of each cell
splat_scalings <- sim.paths@colData$ExpLibSize
scalings <- splat_scalings / mean(splat_scalings)

# pseudotime of each cell
pseudotime <- sim.paths@colData$Step

path.order <- get_path_order(lineage_tree)
path.offsets <- unlist(lapply(lineage_tree, FUN = trace_back, path.from=lineage_tree)) * branch_length
pseudotime <- pseudotime + path.offsets[branches]

## save the same way as PROSSTT
job <- paste(JobFolder, "/", JobName, sep = "")
# the simulation
write.table(x = X, file = paste(job, "simulation.txt", sep="_"), sep = "\t")

# cell parameters/annotation
cellparams <- data.frame(pseudotime, branches, scalings)
colnames(cellparams) <- c("pseudotime", "branches", "scalings")
rownames(cellparams) <- rownames(X)
write.table(x = cellparams, file = paste(job, "cellparams.txt", sep="_"), sep = "\t")

# simulation parameters
topology <- branch_connectivity(lineage_tree)
save_splat_parameters(paste(job, "params.txt", sep="_"), G, branch_length, topology, seed)