#!/usr/bin/Rscript

suppressPackageStartupMessages(library(SLICER))
library("optparse")

LOG_MESSAGE <- ""

option_list <- list(
  make_option(c("-t", "--hhtree"), type="character", default="~/Documents/repos/hhtree", 
              help="Location of the package [default= %default]", metavar="/path/to/hhtree"),
  make_option(c("-o", "--out"), type="character", 
              help="Location where the output is stored", metavar="/output/folder"),
  make_option(c("-j", "--job"), type="character", 
              help="job name", metavar="jobname"),
  make_option(c("-d", "--dimensions"), type="integer", default = 2,
              help="how many dimensions to use for the embedding [default = %default]", metavar="INT"),
  make_option(c("-l", "--log"), action="store_true", default = FALSE,
              help="Toggle to log the expression matrix [default = %default]")
); 

opt_parser <- OptionParser(option_list = option_list);
opt <- parse_args(opt_parser);

hhtree <- opt$hhtree
JobFolder <- opt$out
JobName <- opt$job
dimensions <- opt$dimensions
iflog <- opt$log

various <- paste(hhtree, "/scripts/various.R", sep="")
evaluat <- paste(hhtree, "/scripts/evaluate_method.R", sep="")

suppressPackageStartupMessages(source(various))
suppressPackageStartupMessages(source(evaluat))

job <- paste(JobFolder, JobName, sep="")
raw <- read.table(file = paste(job, "simulation.txt", sep="_"), sep="\t", header = T, row.names = 1, stringsAsFactors = T)
data <- as.matrix(raw)
sum1 <- apply(data, 1, sum)
scalings <- sum1/mean(sum1)
data <- (1/scalings)*data

if(iflog) {
  data <- log(data)
}

cell_params <- read.table(file = paste(job, "cellparams.txt", sep="_"), sep="\t", header = T, row.names = 1)
time <- cell_params$pseudotime - min(cell_params$pseudotime) + 1
labels <- cell_params$branches + 1

# SLICER ----
genes <- select_genes(data)
k <- select_k(data[,genes], kmin = 5)
data_lle <- lle(data[,genes], m = dimensions, k)$Y
LOG_MESSAGE <- paste(LOG_MESSAGE, "dimensions:", dim(data_lle), "\n")
assign("LOG_MESSAGE", LOG_MESSAGE, envir = .GlobalEnv)
data_graph <- conn_knn_graph(data_lle,k)
ends <- find_extreme_cells(data_graph, data_lle)
start <- 1
cells_ordered <- cell_order(data_graph, start)
stimes <- cells_ordered

sbranches <- tryCatch({assign_branches(data_graph, start)},
                      error = function(cond) {
                        print("SLICER failed with message:")
                        print(cond)
                        LOG_MESSAGE <- paste(LOG_MESSAGE, cond, "\n")
                        assign("LOG_MESSAGE", LOG_MESSAGE, envir = .GlobalEnv)
                        return(NA)
                      })

par_loc <- paste(job, "params.txt", sep="_")
mname <- "SLICER"
if (iflog) {
  mname <- paste(mname, "log", sep="_")
}
res <- evaluate_method(mname, sbranches, stimes, cell_params, par_loc)
read_write_output(res, paste(JobFolder, JobName, "_eval.txt", sep=""))

# save coordinates, branch assignments and pseudotime predictions
coords <- data_lle
branches <- sbranches
times <- stimes
where <- paste(JobFolder, JobName, "_", mname, sep="")

write.csv(coords, paste(where, "coordinates.csv", sep="_"))
write.csv(branches, paste(where, "branches.csv", sep="_"))
write.csv(times, paste(where, "times.csv", sep="_"))
write(x = LOG_MESSAGE, file = paste(where, ".log", sep=""))