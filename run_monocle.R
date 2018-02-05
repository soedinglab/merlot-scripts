#!/usr/bin/Rscript

library("optparse")
suppressPackageStartupMessages(library(monocle))

LOG_MESSAGE <- ""

monocle_dimensions <- function(data, dimensions, auto_params=TRUE) {
  if (dimensions < 2) {
    return (NA)
  }
  tryCatch( {
      tmp = reduceDimension(data, max_components = dimensions,
        auto_param_selection = auto_params)

      LOG_MESSAGE <- paste(LOG_MESSAGE, "ran with", dimensions,
        "dimensions, auto_param_selection =", auto_params, "\n")
      assign("LOG_MESSAGE", LOG_MESSAGE, envir = .GlobalEnv)
      return (tmp)
    },
    error = function(cond) {
      LOG_MESSAGE <- paste(LOG_MESSAGE, cond, "\n")
      assign("LOG_MESSAGE", LOG_MESSAGE, envir = .GlobalEnv)

      if (!auto_params) {
        LOG_MESSAGE <- paste(LOG_MESSAGE, dimensions,
          "dimensions are too many, trying with", dimensions-1, "\n")
        assign("LOG_MESSAGE", LOG_MESSAGE, envir = .GlobalEnv)

        return (monocle_dimensions(data, dimensions-1))
      } else {
        LOG_MESSAGE <- paste(LOG_MESSAGE, dimensions,
          "trying with auto_params=F\n")
        assign("LOG_MESSAGE", LOG_MESSAGE, envir = .GlobalEnv)

        return (monocle_dimensions(data, dimensions, auto_params = FALSE))
      }
    })
}

option_list <- list(
  make_option(c("-o", "--out"), type="character", 
              help="Location where the output is stored", metavar="/output/folder"),
  make_option(c("-j", "--job"), type="character", 
              help="job name", metavar="jobname"),
  make_option(c("-d", "--dimensions"), type="integer", default=2,
              help="how many dimensions to use for the embedding [default = %default]", metavar="INT"),
  make_option("--unconstrained", action = "store_true", default = FALSE,
              help = "Toggle to let monocle do unconstrained dimensionality reduction [default = %default]")
);

opt_parser <- OptionParser(option_list=option_list);
opt <- parse_args(opt_parser);

JobFolder <- opt$out
JobName <- opt$job
unconstr <- opt$unconstrained

if (unconstr) {
  dimensions <- 2
} else {
  dimensions <- opt$dimensions
}

# JobFolder = "/data/niko/final/benchmark10/test82/"
# JobName = "test82"
# dimensions = 11
# unconstr = FALSE

job <- paste(JobFolder, JobName, sep="")
raw <- read.table(file=paste(job, "simulation.txt", sep="_"), sep="\t", header=T, row.names=1, stringsAsFactors=T)
data <- as.matrix(raw)
exprs <- t(as.matrix(raw))
# monocle does the normalization itself, so no need to
# sum1 = apply(exprs, 1, sum)
# scalings = sum1/mean(sum1)
# exprs = (1/scalings)*exprs

data <- newCellDataSet(exprs, expressionFamily = negbinomial(), lowerDetectionLimit = 1)

# run monocle ----
data <- estimateSizeFactors(data)
data <- tryCatch( {estimateDispersions(data)}, 
                  error = function(cond) {
                    print(cond)
                    LOG_MESSAGE <- paste(LOG_MESSAGE, cond, "\n")
                    assign("LOG_MESSAGE", LOG_MESSAGE, envir = .GlobalEnv)
                    return(NA)
                  })

if (!is.na(data)) {
  data <- detectGenes(data, min_expr=0.1)
  
  disp_table <- dispersionTable(data)
  ordering_genes <- subset(disp_table, mean_expression >= 0.5 & dispersion_empirical >= 2*dispersion_fit)$gene_id
  
  data <- setOrderingFilter(data, ordering_genes)
  if (unconstr) {
    data <- monocle_dimensions(data, dimensions = 2)
    LOG_MESSAGE <- paste(LOG_MESSAGE, "ran with 2 dimensions\n")
    assign("LOG_MESSAGE", LOG_MESSAGE, envir = .GlobalEnv)
  } else {
    data <- monocle_dimensions(data, dimensions)
  }
  data <- orderCells(data, reverse=FALSE)
}

# coords <- t(reducedDimS(data))
# now calculate tree:
# prepare data for TreeTopology.py
JobName <- paste(JobName, "monocl2", sep = "_")
if (unconstr) {
  JobName <- paste(JobName, "unc", sep = "_")
}

saveRDS(object=data, file = paste(JobFolder, JobName, sep=""))
# run TreeTopology.py on JobFolder_JobName_CellCoord
write(x = LOG_MESSAGE, file = paste(JobFolder, JobName, ".log", sep=""))