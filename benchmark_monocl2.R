#!/usr/bin/Rscript

suppressPackageStartupMessages(library(monocle))
library("optparse")

LOG_MESSAGE <- ""

option_list <- list(
  make_option(c("-t", "--hhtree"), type = "character", default = "~/Documents/repos/hhtree", 
              help = "Location of the package [default= %default]", metavar = "/path/to/hhtree"),
  make_option(c("-o", "--out"), type = "character", 
              help = "Location where the output is stored", metavar = "/output/folder"),
  make_option(c("-j", "--job"), type = "character", 
              help = "job name", metavar = "jobname"),
  make_option(c("--unconstrained"), action="store_true", default=FALSE,
              help="Run monocle2 on default, with only 2 dimensions. [default = %default]"),
  make_option(c("-d", "--dimensions"), type = "integer", default = 2,
              help = "how many dimensions to use for the embedding [default = %default]", metavar = "INT")
); 

opt_parser <- OptionParser(option_list = option_list);
opt <- parse_args(opt_parser);

# hhtree <- "~/Documents/repos/hhtree"
# JobFolder <- "/data/niko/final/benchmark10/test12/"
# JobName <- "test12"
# dimensions <- 11

hhtree <- opt$hhtree
JobFolder <- opt$out
JobName <- opt$job
unconstrained <- opt$unconstrained
dimensions <- opt$dimensions

methname <- "monocl2"
if (unconstrained) {
  methname <- paste(methname, "unc", sep = "_")
  dimensions <- 2
}
various <- paste(hhtree, "/scripts/various.R", sep = "")
evaluat <- paste(hhtree, "/scripts/evaluate_method.R", sep = "")

suppressPackageStartupMessages(source(various))
suppressPackageStartupMessages(source(evaluat))

job <- paste(JobFolder, JobName, sep = "")

cell_params <- read.table(file = paste(job, "cellparams.txt", sep = "_"), sep = "\t", header = T, row.names = 1)
time <- cell_params$pseudotime - min(cell_params$pseudotime) + 1
labels <- cell_params$branches + 1

# set up monocle ----
# days <- data.frame(day = assign_days(length(cells), 5), row.names = cells)
# pd <- new("AnnotatedDataFrame", data = days)

where <- paste(JobFolder, JobName, "_", methname, sep = "")
data <- readRDS(where)

LOG_MESSAGE <- paste(LOG_MESSAGE, "dimensions:", dim(data@reducedDimS), "\n")
assign("LOG_MESSAGE", LOG_MESSAGE, envir = .GlobalEnv)

mbranches <- NA
mtimes <- NA
if (!is.na(data)) {  
  mbranches <- pData(data)$State
  mtimes <- pData(data)$Pseudotime

  # help monocle out with determining pseudotime:
  # ideal: make it start pseudotime at cell 1
  # but monocle sometimes will put cell_1 in an inner branch.
  new_root <- pData(data)$State[1]
  tryCatch( {
    data2 <- orderCells(data, root_state = new_root)
    mbranches <- pData(data2)$State
    mtimes <- pData(data2)$Pseudotime}, 
    error = function(cond) {
      print("heuristic found inner branch D:")
      print(cond)
      LOG_MESSAGE <- paste(LOG_MESSAGE, cond, "\n")
      assign("LOG_MESSAGE", LOG_MESSAGE, envir = .GlobalEnv)
      })
}

par_loc <- paste(job, "params.txt", sep = "_")
res <- evaluate_method(methname, mbranches, mtimes, cell_params, par_loc)
read_write_output(res, paste(JobFolder, JobName, "_eval.txt", sep = ""))

# save coordinates, branch assignments and pseudotime predictions
coords <- t(reducedDimS(data))
branches <- mbranches
times <- mtimes
where <- paste(JobFolder, JobName, "_", methname, sep = "")

write.csv(coords, paste(where, "coordinates.csv", sep = "_"))
write.csv(branches, paste(where, "branches.csv", sep = "_"))
write.csv(times, paste(where, "times.csv", sep = "_"))
write(x = LOG_MESSAGE, file = paste(where, ".log", sep = ""))