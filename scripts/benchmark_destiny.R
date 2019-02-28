#!/usr/bin/Rscript

library(destiny)
library("optparse")

LOG_MESSAGE <- ""

option_list <- list(
  make_option(c("-t", "--mscripts"), type = "character", default = "~/Documents/repos/mscripts",
              help = "Location of the package [default= %default]", metavar = "/path/to/mscripts"),
  make_option(c("-o", "--out"), type = "character",
              help = "Location where the output is stored", metavar = "/output/folder"),
  make_option(c("-j", "--job"), type = "character",
              help = "job name", metavar = "jobname"),
  make_option(c("-d", "--dimensions"), type = "integer", default = 2,
              help = "how many dimensions to use for the embedding [default = %default]", metavar = "INT"),
  make_option(c("-l", "--log"), action = "store_true", default = FALSE,
              help = "Toggle to log-transform the expression matrix [default = %default]"),
  make_option(c("-n", "--knn"), action = "store_true", default = FALSE,
              help = "Toggle to use a simple heuristic for choosing a k for diffusion maps. [default = %default]"),
  make_option(c("-e", "--select"), type="character", default = "none",
              help = "The <name> (merlot|monocle) of a file with a subset of the dataset with only informative genes.")
);

opt_parser <- OptionParser(option_list = option_list);
opt <- parse_args(opt_parser);

mscripts <- opt$mscripts
JobFolder <- opt$out
JobName <- opt$job
dimensions <- opt$dimensions
iflog <- opt$log
knn <- opt$knn
select <- opt$select

# mscripts <- "~/Documents/repos/mscripts"
# JobFolder <- "/data/niko/final/benchmark2/test0/"
# JobName <- "test0"
# dimensions <- 3
# iflog <- TRUE
# knn <- FALSE
# select <- "none"

various <- paste(mscripts, "/scripts/various.R", sep = "")
evaluat <- paste(mscripts, "/scripts/evaluate_method.R", sep = "")

suppressPackageStartupMessages(source(various))
suppressPackageStartupMessages(source(evaluat))

job <- paste(JobFolder, JobName, sep = "")

cell_params <- read.table(file = paste(job, "cellparams.txt", sep = "_"), sep = "\t", header = T, row.names = 1)
time <- cell_params$pseudotime - min(cell_params$pseudotime) + 1
start <- min(which(time == min(time)))
labels <- cell_params$branches + 1

# destiny ----
diffmap <- ""
if (iflog) {diffmap <- paste(diffmap, "_log", sep = "")}
if (knn) {diffmap <- paste(diffmap, "_k", sep = "")}
if (select == "merlot") {
  diffmap <- paste(diffmap, "_sel_merlot", sep = "")
} else if (select == "monocle") {
  diffmap <- paste(diffmap, "_sel_monocle", sep = "")
} else if (select == "seurat") {
  diffmap <- paste(diffmap, "_sel_seurat", sep = "")
}
diffmap <- paste("destiny", diffmap, sep = "")

# check if present:
res_file <- paste(JobFolder, JobName, "_eval.txt", sep = "")
eval <- t(as.matrix(read.table(res_file, check.names=FALSE, stringsAsFactors = FALSE)))
if (diffmap %in% colnames(eval)) {
  stop(paste(diffmap, "already evaluated!"))
}

dif <- readRDS(file = paste(job, diffmap, sep = "_"))

LOG_MESSAGE <- paste(LOG_MESSAGE, "dimensions:", dim(dif@eigenvectors), "\n")
assign("LOG_MESSAGE", LOG_MESSAGE, envir = .GlobalEnv)

dpt <- tryCatch( {DPT(dif)}, 
                error = function(cond) {
                  print("branching assignment failed")
                  print(cond)
                  LOG_MESSAGE <- paste(LOG_MESSAGE, cond, "\n")
                  assign("LOG_MESSAGE", LOG_MESSAGE, envir = .GlobalEnv)
                  return(1)
                 })

if (!is.atomic(dpt)) {
  dbranches <- dpt@branch[,1]
  dtimes <- dpt[, start]
  dbranches[is.na(dbranches)] <- 0
  dbranches[dbranches == 0] <- max(dbranches) + 1
} else {
  dbranches <- NA
  dtimes <- NA
}

par_loc <- paste(job, "params.txt", sep = "_")
res <- evaluate_method(diffmap, dbranches, dtimes, cell_params, par_loc)
read_write_output(res, paste(JobFolder, JobName, "_eval.txt", sep = ""))

# save coordinates, branch assignments and pseudotime predictions
coords <- eigenvectors(dif)[, 1:dimensions]
branches <- dbranches
times <- dtimes
# where <- paste(job, diffmap, sep = "_")

# write.csv(coords, paste(where, "coordinates.csv", sep = "_"))
write.csv(branches, paste(job, diffmap, "branches.csv", sep = "_"))
write.csv(times, paste(job, diffmap, "times.csv", sep = "_"))
write(x = LOG_MESSAGE, file = paste(job, "_", diffmap, ".log", sep = ""))