#!/usr/bin/Rscript
library(igraph)

evaluate_method <- function(method, mbranches, mtimes, cell_params, par_loc, returnNA=FALSE) {

  funclist <- list()
  funclist[[1]] <- randInd_manual
  funclist[[2]] <- matthews_cor
  funclist[[3]] <- f_measure
  funclist[[4]] <- jaccard
  funclist[[5]] <- fowkles_mallows
  funclist[[6]] <- adjusted_mi

  timefuncs <- list()
  timefuncs[[1]] <- goodman_kruskal_index
  timefuncs[[2]] <- goodman_kruskal_index
  timefuncs[[3]] <- kendall_index
  timefuncs[[4]] <- kendall_index

  funcnames <- c("rand index", "matthews corr. coef.",
                 "F1 measure", "jaccard index", "fowkles-mallows index", "adjusted MI")

  timenames <- c("goodman-kruskall (unweighted)", "goodman-kruskall (weighted)",
                 "kendall index (unweighted)", "kendall index (weighted)")

  res <- data.frame(matrix(NA, nrow = length(funcnames)+length(timenames)+2, ncol = 1))
  rownames(res) <- c(funcnames, timenames, "LPGK", "#branches")

  colnames(res) <- method
  # print(res)

  if (returnNA) { return(t(res)) }

  mstatus <- assign_status(mbranches, cell_params$branches+1)

  for (i in 1:length(funcnames)) {
    res[i, 1] <- funclist[[i]](mbranches, cell_params$branches+1, mstatus)
  }

  ltimes <- cell_params$pseudotime - min(cell_params$pseudotime) + 1
  for (i in 1:length(timenames)) {
    weighted <- (i %% 2 == 0)
    res[length(funcnames)+i, 1] <- timefuncs[[i]](mtimes, ltimes, weighted=weighted)
  }

  # longest path goodman kruskal:
  lp_id = get_lpgk_indices(par_loc, cell_params)
  lp_cells = rownames(cell_params)[lp_id]
  ltimes = cell_params[lp_cells,]$pseudotime
  mtimes = mtimes[lp_id]
  res[length(funcnames) + length(timenames) + 1,] <- goodman_kruskal_index(mtimes, ltimes, weighted=FALSE)

  # add number of predicted branches
  num_branches <- length(levels(as.factor(mbranches)))
  res[length(funcnames) + length(timenames) + 2,] <- num_branches

  return(t(res))
}

get_lpgk_indices <- function(par_loc, cell_params) {
  cell_params$branches = cell_params$branches + 1
  time = cell_params$pseudotime - min(cell_params$pseudotime) + 1
  
  StartCell = rownames(cell_params)[which(time==min(time))][1]
  EndCell = rownames(cell_params)[which(time==max(time))][1]

  StartBranch = cell_params[StartCell,]$branches
  EndBranch = cell_params[EndCell,]$branches
  # loc = "~/Desktop/modern_sims/sim3/sim3_params.txt"
  results = system(sprintf("awk '{if($0~/topology/) print}' %s", par_loc), intern = T)
  digits = gregexpr("[0-9]+", results)
  Branches = as.numeric(regmatches(results, digits)[[1]])
  # if the simulation contains a linear topology then the line will read
  # topology: []
  # of course this means that digits will not contain any numbers,
  # so the regex will not match and we have to account for that
  if (length(Branches) == 0) {
    Branches <- 1
  } else {
    Branches = Branches + 1
  }

  # Create adjacency matrix and calculate shortest path in between start and end branches
  BranchConnections = matrix(0, max(Branches), max(Branches))
  for (i in seq(from=1, to=length(Branches), by=2)) {
    BranchConnections[Branches[i], Branches[i+1]]=1
    BranchConnections[Branches[i+1], Branches[i]]=1
  }

  Graph_Branches=graph_from_adjacency_matrix(BranchConnections)
  LongestPath=unlist(shortest_paths(Graph_Branches, from=StartBranch, to=EndBranch)$vpath)

  LongestPathCellsID=which(cell_params$branches %in% LongestPath)
  return(LongestPathCellsID)
}

read_write_output <- function(res, out) {
  if (!file.exists(out)) {
    write.table(res, out)
  } else {
    old <- read.table(out, check.names=FALSE)
    res <- rbind(old, res)
  
    write.table(res, out, append=F)
  }
}


#' Finds all paths from the root to an endpoint for MERLoT
#'
#' Finds all paths that lead from the lineage tree root to an endpoint and the
#' cells on them
#' 
#' @param tree_name suffix of the binary file that contains the predicted tree
#' @param job location and prefix of the file that contains the predicted tree
#' @param embed flag for MERLoT: embed tree or no
#'
#' @return a list of cell IDs, one for each terminal path
#'
#' @export
#'
#' @importFrom igraph graph_from_edgelist shortest_paths
#' @importFrom utils read.table
#' @importFrom merlot GenesSpaceEmbedding
merlot_trajectories <- function(tree_name, job,  embed = TRUE) {
  # if (!requireNamespace("merlot", quietly = TRUE)) {
  #   print("The package `merlot` is not installed but is required for this function.")
  #   return(NA)
  # }
  cell_params <- read.table(file = paste(job, "cellparams.txt", sep = "_"), sep = "\t", header = T, row.names = 1)
  merlot <- readRDS(paste(job, tree_name, sep = ""))

  if (embed == TRUE) {
    raw_data <- read.table(file = paste(job, "simulation.txt", sep = "_"),
                           sep = "\t", header = T, row.names = 1, stringsAsFactors = T)
    X <- as.matrix(raw_data)
    sum1 <- apply(X, 1, sum)
    scalings <- sum1/mean(sum1)
    X <- (1/scalings)*X
    X <- log(X + 1)
    merlot <- merlot::GenesSpaceEmbedding(ExpressionMatrix = X, ElasticTree = merlot)
  }

  g <- igraph::graph_from_edgelist(merlot$Edges, directed = FALSE)
  endpoints <- merlot$Topology$Endpoints

  res <- list()
  n <- length(endpoints)

  for (i in 1:(n-1) ) {
    end_from <- endpoints[i]
    for (j in (i+1):n) {
      end_to <- endpoints[j]
      paths <- igraph::shortest_paths(g, end_from, end_to)
      nodes_in_path <- paths$vpath[[1]]
      cells_in_path <- merlot$Cells2TreeNodes[,2] %in% nodes_in_path
      k <- (n*(n-1)/2) - (n-(i-1))*((n-(i-1))-1)/2 + j - i
      # print(c(i, j, k))
      res[[k]] <- rownames(cell_params)[cells_in_path]
    }
  }
  res
}

#' Finds all paths from the root to an endpoint for TSCAN
#'
#' Finds all paths that lead from the lineage tree root to an endpoint and the
#' cells on them
#' 
#' @param tree_name suffix of the binary file that contains the predicted tree
#' @param job location and prefix of the file that contains the predicted tree
#' @param embed flag for MERLoT: embed tree or no. Needed here for batch runs
#'
#' @return a list of cell IDs, one for each terminal path
#'
#' @export
#'
#' @importFrom igraph shortest_paths degree
#' @importFrom utils read.table
tscan_trajectories <- function(tree_name, job, embed=FALSE) {
  cell_params <- read.table(file = paste(job, "cellparams.txt", sep = "_"), sep = "\t", header = T, row.names = 1)
  tscanned <- readRDS(paste(job, tree_name, sep = ""))

  g <- tscanned$MSTtree
  degrees <- unname(igraph::degree(g))
  endpoints <- which(degrees == 1)

  res <- list()
  n <- length(endpoints)

  for (i in 1:(n-1) ) {
    end_from <- endpoints[i]
    for (j in (i+1):n) {
      end_to <- endpoints[j]
      paths <- igraph::shortest_paths(g, end_from, end_to)
      nodes_in_path <- paths$vpath[[1]]
      cells_in_path <- tscanned$clusterid %in% nodes_in_path
      k <- (n*(n-1)/2) - (n-(i-1))*((n-(i-1))-1)/2 + j - i
      # print(c(i, j, k))
      res[[k]] <- rownames(cell_params)[cells_in_path]
    }
  }
  res
}

#' Finds all paths from the root to an endpoint for Monocle
#'
#' Finds all paths that lead from the lineage tree root to an endpoint and the
#' cells on them
#' 
#' @param tree_name suffix of the binary file that contains the predicted tree
#' @param job location and prefix of the file that contains the predicted tree
#' @param embed flag for MERLoT: embed tree or no. Needed here for batch runs
#'
#' @return a list of cell IDs, one for each terminal path
#'
#' @export
#'
#' @importFrom igraph shortest_paths degree
#' @importFrom utils read.table
monocle_trajectories <- function(tree_name, job, embed=FALSE) {
  cell_params <- read.table(file = paste(job, "cellparams.txt", sep = "_"), sep = "\t", header = T, row.names = 1)
  monocle2 <- readRDS(paste(job, tree_name, sep = ""))

  g <- monocle2@minSpanningTree
  degrees <- unname(igraph::degree(g))
  endpoints <- which(degrees == 1)

  res <- list()
  n <- length(endpoints)

  for (i in 1:(n-1) ) {
    end_from <- endpoints[i]
    for (j in (i+1):n) {
      end_to <- endpoints[j]
      paths <- igraph::shortest_paths(g, end_from, end_to)
      nodes_in_path <- paths$vpath[[1]]
      cells2nodes <- monocle2@auxOrderingData$DDRTree$pr_graph_cell_proj_closest_vertex
      cells_in_path <- cells2nodes %in% nodes_in_path
      k <- (n*(n-1)/2) - (n-(i-1))*((n-(i-1))-1)/2 + j - i
      # print(c(i, j, k))
      res[[k]] <- rownames(cell_params)[cells_in_path]
    }
  }
  res
}

#' Calculates truth table for cells in longest overlapping paths
#'
#' Finds the predicted trajectory with the largest overlap with the simulated
#' longest path and then calculates a truth table for this prediction.
#'
#' @param tree_name suffix of the binary file that contains the predicted tree
#' @param job location and prefix of the file that contains the predicted tree
#' @param tool_function function to locate longest trajectory in prediction
#' @param embed flag for MERLoT: embed tree or no. Needed here for batch runs
#'
#' @return (flat) truth table
#'
#' @export
#'
#' @importFrom utils read.table
on_longest_path <- function(tree_name, job, tool_function, embed = FALSE) {
  par_loc <- paste(job, "params.txt", sep = "_")

  # read the simulation parameters
  cell_params <- read.table(file = paste(job, "cellparams.txt", sep = "_"), sep = "\t", header = T, row.names = 1)
  labels <- cell_params$branches + 1

  lp_id <- get_lpgk_indices(par_loc, cell_params)
  lp_cells <- rownames(cell_params)[lp_id]
  other_id <- which(!seq(length(labels)) %in% lp_id)
  other_cells <- rownames(cell_params)[other_id]

  res <- tool_function(tree_name, job, embed)

  best_traj <- 0
  max_overlap <- 0
  for (i in seq_along(res)) {
    trajectory <- res[[i]]
    overlap <- sum(lp_cells %in% trajectory)
    if (overlap > max_overlap) {
      best_traj <- i
      max_overlap <- overlap
    }
  }

  trajectory <- res[[best_traj]]
  outside_id <- (!rownames(cell_params) %in% res[[best_traj]])
  outside <- rownames(cell_params)[outside_id]

  TP <- sum(trajectory %in% lp_cells)
  FP <- sum(trajectory %in% other_cells)
  FN <- sum(outside %in% lp_cells)
  TN <- sum(outside %in% other_cells)
  c(TP, TN, FN, FP)
}

