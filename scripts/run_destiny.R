#!/usr/bin/Rscript

library("optparse")
suppressPackageStartupMessages(library(FNN))
suppressPackageStartupMessages(library(destiny))

LOG_MESSAGE <- ""

optimize_k <- function(data, d) {
  opt_k <- -10
  k_gap <- 0

  pb <- txtProgressBar(min = 5, max = 100, style = 3)

  print("optimizing k...")
  for (local_k in seq(from = 5, to = 100, by = 5)) {
    setTxtProgressBar(pb, local_k)
    dm <- tryCatch({
      DiffusionMap(data, density_norm = T, sigma = "local", k = local_k)
      }
      , error = function(cond) {
        LOG_MESSAGE <- paste(LOG_MESSAGE, cond, "\n")
        assign("LOG_MESSAGE", LOG_MESSAGE, envir = .GlobalEnv)
        return(NULL)
        })
    if (is.null(dm)) {
      next
    }
    eigv <- eigenvalues(dm)
    # plot(eigv)
    local_k_gap <- abs((eigv[d+1]-eigv[d+2]) - (eigv[d]-eigv[d+1]))
    # print(local_k_gap)
    if(local_k_gap > k_gap)
    {
      opt_k <- local_k
      k_gap <- local_k_gap
    }
  }
  close(pb)

  # if (opt_k < 0) {
  #   opt_k <- find_dm_k(dim(data)[1] - 1)
  #   msg = paste("could not find optimal k, using", opt_k)
  #   LOG_MESSAGE <- paste(LOG_MESSAGE, msg, "\n")
  # }
  print(paste("found", opt_k))
  return(opt_k)
}


option_list <- list(
  make_option(c("-o", "--out"), type = "character", help = "Location where the output is stored", metavar = "/output/folder"),
  make_option(c("-j", "--job"), type = "character", help = "job name", metavar = "jobname"),
  make_option(c("-d", "--dimensions"), type = "integer", default = 2, help = "how many dimensions to use for the embedding [default = %default]", metavar = "INT"),
  make_option(c("-l", "--log"), action = "store_true", default = FALSE, help = "Toggle to log the expression matrix [default = %default]"),
  make_option(c("-n", "--knn"), action = "store_true", default = FALSE, help = "Toggle to use a simple heuristic for choosing a k for diffusion maps. [default = %default]"),
  make_option(c("-s", "--select"), type="character", default = "none", help = "The <name> of a file with a subset of the dataset with only informative genes.")
); 

opt_parser <- OptionParser(option_list = option_list);
opt <- parse_args(opt_parser);

JobFolder <- opt$out
JobName <- opt$job
dimensions <- opt$dimensions
knn <- opt$knn
iflog <- opt$log
select <- opt$select

# JobFolder <- "/data/niko/test/benchmark1/splat0/"
# JobName <- "splat0"
# dimensions <- 2
# knn <- TRUE
# iflog <- TRUE
# select <- "merlot"

job <- paste(JobFolder, JobName, sep = "")

preselected <- FALSE
sim_name <- ""
if (select == "none") {
  sim_name <- paste(job, "simulation.txt", sep = "_")
} else {
  preselected <- TRUE
  sim_name <- paste(job, "_simulation_sel_", select, ".txt", sep = "")
}
# print(sim_name)

raw <- read.table(file = sim_name, sep = "\t", header = T, row.names = 1, stringsAsFactors = T)
data <- as.matrix(raw)

if (!preselected) {
  sum1 <- apply(data, 1, sum)
  scalings <- sum1/mean(sum1)
  data = (1/scalings)*data
}

if (iflog) {
  if (select != "seurat") {data <- log(data+1)}
} else {
  if (select == "seurat") {data <- exp(data)-1}
}

if (knn) {
  k <- optimize_k(data, dimensions)
  dif <- DiffusionMap(data, density_norm = T, verbose = T, sigma = "local", k = k)
} else {
  dif <- DiffusionMap(data)
}

JobName <- paste(JobName, "destiny", sep = "_")

if (iflog) {JobName <- paste(JobName, "_log", sep = "")}
if (knn) {JobName <- paste(JobName, "_k", sep = "")}

if (preselected) {
  JobName <- paste(JobName, "sel", select, sep = "_")
  saveRDS(dif, file = paste(JobFolder, JobName, sep = ""))
} else {
  saveRDS(dif, file = paste(JobFolder, JobName, sep = ""))
}
write(x = LOG_MESSAGE, file = paste(JobFolder, JobName, ".log", sep = ""))