#!/usr/bin/Rscript

library(TSCAN)
library("optparse")

option_list <- list(
  make_option(c("-t", "--mscripts"), type = "character", default = "~/Documents/repos/mscripts", 
              help = "Location of the package [default= %default]", metavar = "/path/to/mscripts"),
  make_option(c("-o", "--out"), type = "character", 
              help = "Location where the output is stored", metavar = "/output/folder"),
  make_option(c("-j", "--job"), type = "character", 
              help = "job name", metavar = "jobname")#,
  # make_option(c("-c", "--clusters"), type = "int", default = 1:50,
  #             help = "number of clusters to search for", metavar = "clusters")
);

opt_parser <- OptionParser(option_list = option_list);
opt <- parse_args(opt_parser);

mscripts <- opt$mscripts
JobFolder <- opt$out
JobName <- opt$job

# mscripts <- "~/Documents/repos/mscripts"
# JobFolder <- "/data/niko/repeat/benchmark1/test2/"
# JobName <- "test2"

# load functions
various <- paste(mscripts, "/scripts/various.R", sep = "")
evaluat <- paste(mscripts, "/scripts/evaluate_method.R", sep = "")
splat_funcs <- paste(mscripts, "/scripts/clusters_to_branches.R", sep = "")
suppressPackageStartupMessages(source(various))
suppressPackageStartupMessages(source(evaluat))
suppressPackageStartupMessages(source(splat_funcs))


# read the data
job <- paste(JobFolder, JobName, sep = "")
raw <- read.table(file = paste(job, "simulation.txt", sep = "_"), sep = "\t", header = T, row.names = 1, stringsAsFactors = T)
X <- t(as.matrix(raw))
N <- dim(X)[2]
par_loc <- paste(job, "params.txt", sep = "_")
methname <- "TSCAN"

# read the simulation parameters
cell_params <- read.table(file = paste(job, "cellparams.txt", sep = "_"), sep = "\t", header = T, row.names = 1)
time <- cell_params$pseudotime - min(cell_params$pseudotime) + 1
labels <- cell_params$branches + 1

procdata <- preprocess(X, minexpr_value = 0.5)
tscanned <- run_with_maximal_clusters(procdata, 50)
# plotmclust(tscanned, show_tree = T)

# get "correct" branch assignments by collapsing co-linear branches
tscan_branches <- sling_branches <- assign_branches(tscanned$MSTtree,
                                                    tscanned$clusterid,
                                                    tscanned$clucenter,
                                                    tscanned$pcareduceres)

# re-run with 2 clusters so that we get global pseudotime
tscanned_time <- exprmclust(procdata, clusternum = 2)
ordered_cells <- TSCANorder(tscanned_time)
scan_times <- matrix(1:N)
cells <- colnames(X)
tscan_times <- scan_times[match(cells, ordered_cells)]

where <- paste(JobFolder, JobName, "_", methname, sep = "")
write.csv(tscan_branches, paste(where, "branches.csv", sep = "_"))
write.csv(tscan_times, paste(where, "times.csv", sep = "_"))

saveRDS(object = tscanned, file = paste(where, sep = "_"))

res <- evaluate_method(methname, tscan_branches, tscan_times, cell_params, par_loc)
where <- paste(JobFolder, JobName, sep = "")
read_write_output(res, paste(where, "_eval.txt", sep = ""))
