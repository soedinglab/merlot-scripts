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

library(monocle)

benchmark <- "benchmark8"
dimensions=9
sim <- "sim18"
job <- paste("~/Documents/data/examine", benchmark, sim, sim, sep="/")

cell_params <- read.table(file = paste(job, "cellparams.txt", sep = "_"), sep = "\t", header = T, row.names = 1)
time <- cell_params$pseudotime - min(cell_params$pseudotime) + 1
start <- min(which(time == min(time)))
labels <- cell_params$branches + 1

full <- read.table(file=paste(job, "simulation.txt", sep="_"), sep="\t", header=T, row.names=1, stringsAsFactors=T)

exprs <- t(as.matrix(full))
data <- newCellDataSet(exprs, expressionFamily = negbinomial(), lowerDetectionLimit = 1)

# run monocle ----
data <- estimateSizeFactors(data)
data <- estimateDispersions(data)
data <- detectGenes(data, min_expr=0.1)
disp_table <- dispersionTable(data)
ordering_genes <- subset(disp_table, mean_expression >= 0.5 & dispersion_empirical >= 2*dispersion_fit)$gene_id
data <- setOrderingFilter(data, ordering_genes)
data <- monocle_dimensions(data, dimensions)
data <- orderCells(data, reverse=FALSE)

# where <- paste(job, "_monocl2", sep = "")
# data <- readRDS(where)

mbranches <- pData(data)$State
mtimes <- pData(data)$Pseudotime

coords <- t(reducedDimS(data))
plot(coords, pch=16, col=mbranches)
plot(coords, pch=16, col=labels)








