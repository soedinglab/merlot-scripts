maximal_nonbranching <- function(itree) {
  degrees <- unname(degree(itree))
  maximal.nonbranching <- list()
  
  for (i in seq_along(V(itree))) {
    v <- V(itree)[i]
    all_neighbors <- neighbors(itree, v)
    if (degrees[i] != 2) {
      for (j in all_neighbors) {
        w <- V(itree)[j]
        thispath <- c(v, w)
        prev <- v
        while (degree(itree, w) == 2) {
          u <- neighbors(itree, w)
          keep <- which(u != prev)
          thispath <- c(thispath, u[keep])
          prev <- w
          w <- u[keep]
        }
        maximal.nonbranching <- push(maximal.nonbranching, thispath)
      }
    }
  }
  
  temp <- lapply(maximal.nonbranching, function(x) sort(as.numeric(x)))
  keep <- !duplicated(temp)
  branches <- maximal.nonbranching[keep]
  branches
}

push <- function(foo, bar){
  foo[[length(foo)+1]] <- bar
  foo
}

assign_branches_legacy <- function(tscan_output) {
  # how to reassign stuff:
  # find all branch points and their nearest neighbors.
  # then reassign all cells of the branch point cluster to
  # the closest neighbor centroid
  itree <- tscan_output$MSTtree
  degrees <- unname(degree(itree))
  scan_branches <- tscan_output$clusterid
  branches <- rep(0, length(scan_branches))
  topology <- maximal_nonbranching(tscan_output$MSTtree)
  
  branch_points <- which(degrees > 2)
  for (i in seq_along(topology)) {
    branch <- topology[[i]]
    # endpoints and intermediate branches can be
    # converted immediately
    convert_at_once <- which(degrees[branch] < 3)
    # translate to cells
    cells <- (scan_branches %in% branch[convert_at_once])
    branches[cells] <- -i
    
    branch_points <- branch[which(degrees[branch] > 2)]
    for (b in branch_points) {
      b_neighbors <- as.numeric(neighbors(itree, b))
      b_centers <- tscan_output$clucenter[b_neighbors, ]
      b_cells <- (scan_branches == b) & (branches==0)
      if (sum(b_cells) == 0) {
        next
      }
      # we have to format cells as a matrix because if there is exactly one 
      # unassigned cell it will be treated as a vector D:
      ndim <- dim(tscan_output$pcareduceres)[2]
      cells <- matrix(tscan_output$pcareduceres[b_cells, ], ncol = ndim)
      
      # lots of help by https://stackoverflow.com/questions/16283347
      distances <- outer(
        1:nrow(cells), 
        1:nrow(b_centers), 
        Vectorize( function(i,j) { 
          sum( (cells[i,] - b_centers[j,])^2 )
        } )
      )
      clusters <- apply(distances, 1, which.min)
      branch_representative <- which(b_neighbors %in% branch)
      branches[b_cells] <- clusters
      branches[branches == branch_representative] <- -i
      branches[branches > 0] <- 0
    }
  }
  branches <- branches * (-1)
  branches
}


assign_branches <- function(g, clusterid, clucenter, coords) {
  # how to reassign stuff:
  # find all branch points and their nearest neighbors.
  # then reassign all cells of the branch point cluster to
  # the closest neighbor centroid
  branches <- rep(0, dim(coords)[1])
  topology <- maximal_nonbranching(g)
  
  # branch_points <- V(g)[which(degrees > 2)]
  for (i in seq_along(topology)) {
    branch <- topology[[i]]
    # endpoints and intermediate branches can be
    # converted immediately
    convert_at_once <- branch[which(degree(g, branch) < 3)]
    # translate to cells
    cells <- (clusterid %in% as.numeric(names(convert_at_once)))
    branches[cells] <- -i
    
    branch_points <- branch[which(degree(g, branch) > 2)]
    for (j in seq_along(branch_points)) {
      b <- branch_points[j]
      b_neighbors <- neighbors(g, b)
      b_centers <- clucenter[b_neighbors, ]
      b_cells <- (clusterid == as.numeric(names(b))) & (branches==0)
      if (sum(b_cells) == 0) {
        next
      }
      # we have to format cells as a matrix because if there is exactly one 
      # unassigned cell it will be treated as a vector D:
      ndim <- dim(coords)[2]
      if (ndim == 1) {
        cells <- matrix(coords[b_cells, ], ncol = ndim)
      } else {
        cells <- as.matrix(coords[b_cells, ])
      }
      
      # lots of help by https://stackoverflow.com/questions/16283347
      distances <- outer(
        1:nrow(cells), 
        1:nrow(b_centers), 
        Vectorize( function(i,j) { 
          sum( (cells[i,] - b_centers[j,])^2 )
        } )
      )
      clusters <- apply(distances, 1, which.min)
      branch_representative <- which(b_neighbors %in% branch)
      branches[b_cells] <- clusters
      branches[branches == branch_representative] <- -i
      branches[branches > 0] <- 0
    }
  }
  branches <- branches * (-1)
  branches
}

run_with_maximal_clusters <- function(procdata, no_clusters) {
  res <- NULL
  res <- tryCatch({
    exprmclust(procdata, clusternum = no_clusters)
  }, error = function(cond) {
    return(run_with_maximal_clusters(procdata, no_clusters - 1))
  })
  return(res)
}
