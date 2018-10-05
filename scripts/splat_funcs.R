#!/usr/bin/Rscript

map2color<-function(x,pal,limits=NULL){
  if(is.null(limits)) limits=range(x)
  pal[findInterval(x,seq(limits[1],limits[2],length.out=length(pal)+1), all.inside=TRUE)]
}

save_splat_parameters <- function(output, genes, pseudotime, topology, random_seed, modules = 0) {
  destination <- file(output)
  
  gene_string <- paste("Genes:", genes, sep = " ")
  
  top <- "["
  for (i in 1:dim(topology)[1]) {
    branch <- topology[i, ]
    top <- paste(top, "[", branch[1], ", ", branch[2], "]", sep = "")
    if (i < dim(topology)[1]) {
      top <- paste(top, ", ", sep = "")
    }
  }
  top <- paste(top, "]", sep = "")
  top_string <- paste("topology:", top, sep = " ")
  
  pt_string <- "pseudotimes: ["
  for (i in 1:(dim(topology)[1] + 1)) {
    pt_string <- paste(pt_string, pseudotime, sep = "")
    if (i < (dim(topology)[1] + 1)) {
      pt_string <- paste(pt_string, ", ", sep = "")
    }
  }
  pt_string <- paste(pt_string, "]", sep = "")
  
  mod_string <- paste("#modules:", modules, sep = " ")
  seed_string <- paste("random seed:", random_seed, sep = " ")
  
  param_text <- c(gene_string, pt_string, top_string, mod_string, seed_string)
  writeLines(param_text, destination)
  close(destination)
}

gen_random_topology <- function(branch_points) {
  total_branches <- 2 * branch_points + 1
  res <- rep(-1, total_branches)
  res[1] <- 0
  
  available_branchpoints <- c(1)
  available_destinations <- 2:total_branches
  
  for (i in 1:branch_points) {
    # select new branch point
    new_branch_point <- sample(available_branchpoints, 1)
    # select destinations
    goes_to <- sample(available_destinations, 2, replace = FALSE)
    # map from-to
    res[goes_to] <- new_branch_point
    
    # update available branchpoints
    available_branchpoints <- append(available_branchpoints, goes_to)
    without_new_branchpoint <- (available_branchpoints != new_branch_point)
    available_branchpoints <- available_branchpoints[without_new_branchpoint]
    
    # update available destinations
    without_new_destinations <- (!(available_destinations %in% goes_to))
    available_destinations <- available_destinations[without_new_destinations]
  }
  return(res)
}

get_path_order <- function(path.from) {
  
  # Transform the vector into a list of (from, to) pairs
  path.pairs <- list()
  for (idx in seq_along(path.from)) {
    path.pairs[[idx]] <- c(path.from[idx], idx)
  }
  
  # Determine the processing order
  # If a path is in the "done" vector any path originating here can be
  # completed
  done <- 0
  while (length(path.pairs) > 0) {
    path.pair <- path.pairs[[1]]
    path.pairs <- path.pairs[-1]
    from <- path.pair[1]
    to <- path.pair[2]
    if (from %in% done) {
      done <- c(done, to)
    } else {
      path.pairs <- c(path.pairs, list(path.pair))
    }
  }
  
  # Remove the origin from the vector
  done <- done[-1]
  
  return(done)
}

trace_back <- function(path.from, n) {
  if (n == 0) {
    return(0)
  } else {
    return(trace_back(path.from, path.from[n]) + 1)
  }
}

branch_connectivity <- function(lineage_tree) {
  total_branches <- length(lineage_tree)
  pairings <- matrix(-1, ncol=2, nrow=total_branches-1)
  branches <- data.frame(from=lineage_tree, to=1:total_branches)
  counter <- 0
  
  for (b in 1:total_branches) {
    connections <- which(branches$from == b)
    if (length(connections) > 0) {
      for (i in seq_along(connections)) {
        pairings[counter + i, ] <- c(b, connections[i])
      }
      counter <- counter + length(connections)
    }
  }
  pairings
}
