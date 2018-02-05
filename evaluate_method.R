#!/usr/bin/Rscript



evaluate_method <- function(method, mbranches, mtimes, cells_params, par_loc, returnNA=FALSE) {

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