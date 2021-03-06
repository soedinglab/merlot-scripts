library(pdist)
library(hashmap)
library(igraph)
library(ggplot2)
library(reshape2)


branches_not_adjacent <- function(branches, topology) {
  res <- array(FALSE, dim = c(length(branches), length(branches)))
  for (i in seq_along(branches)) {
    for (j in seq_along(branches)) {
      y <- branches[i]
      z <- branches[j]
      res[i, j] <- any(apply(topology, 1, FUN = function(x) { all((x == c(y, z)) | (x == c(z, y))) | y == z }))
    }
  }
  
  !res
}

branch_eucl_length <- function(coords, cellparams) {
  branches <- unique(cellparams$branches)
  branch_length <- hashmap(keys = branches, values = rep(0, length(branches)))
  for (branch in branches) {
    min_time <- min(cellparams$pseudotime[cellparams$branches == branch])
    max_time <- max(cellparams$pseudotime[cellparams$branches == branch])
    minima <- which(cellparams$branches == branch & cellparams$pseudotime == min_time)
    maxima <- which(cellparams$branches == branch & cellparams$pseudotime == max_time)
    dists <- pdist(coords, indices.A = minima, indices.B = maxima)
    branch_length[[branch]] <- mean(dists@dist)
  }
  branch_length
}

# euclidean distance from the beginning of b1 to the end of b2
direct_dist <- function(coords, cellparams, b1, b2) {
  min_time <- min(cellparams$pseudotime[cellparams$branches == b1])
  max_time <- max(cellparams$pseudotime[cellparams$branches == b2])
  from <- which(cellparams$branches == b1 & cellparams$pseudotime == min_time)
  to <- which(cellparams$branches == b2 & cellparams$pseudotime == max_time)
  direct <- pdist(coords, indices.A = from, indices.B = to)
  direct@dist
}




splatter_from_origin <- list()
for (l in 1:11) {
  splatter_from_origin[[l]] <- list()
}

prosstt_from_origin <- list()
for (l in 1:11) {
  prosstt_from_origin[[l]] <- list()
}

pb <- txtProgressBar(min = 0, max = 1000, style = 3)
for (b in 1:10) {
  benchmark_dir <- paste("/path/to/splatter/benchmark", b, "/", sep="")
  dim_keep <- b + 2
  for (i in 0:99) {
    JobName <- paste("sim", i, sep="")
    JobFolder <- paste(benchmark_dir, JobName, "/", sep="")
    job <- paste(JobFolder, JobName, sep="")
    og_cellparams <- read.table(file = paste(job, "cellparams.txt", sep = "_"), sep = "\t", header = T, row.names = 1)
    topology <- read_prosstt_topology(job)
    
    dm_name <- paste(job, "destiny_log", sep = "_")
    dm <- readRDS(dm_name)
    data <- eigenvectors(dm)[, 1:dim_keep]
    
    N <- dim(data)[1]
    sample_down <- sample(N, N/2)
    data <- data[sample_down, ]
    cellparams <- data.frame("pseudotime"= og_cellparams$pseudotime[sample_down], "branches"=og_cellparams$branches[sample_down])
    
    branch_lengths <- branch_eucl_length(data, cellparams)
    
    g <- igraph::graph_from_edgelist(topology)
    root <- which(sapply(sapply(V(g), function(x) neighbors(g, x, mode="in")), length) == 0)
    
    paths <- igraph::all_simple_paths(g, root, V(g))
    for (k in seq_along(paths)) {
      p <- paths[[k]]
      direct_distance <- mean(direct_dist(data, cellparams, p[1], tail(p, n = 1)))
      tree_distance <- mean(branch_lengths[[p]])
      splatter_from_origin[[length(p)]] <- append(splatter_from_origin[[length(p)]], direct_distance / tree_distance)
    }
    
    setTxtProgressBar(pb, (b-1) * 100 + (i+1))
  }
}

pb <- txtProgressBar(min = 0, max = 1000, style = 3)
for (b in 1:10) {
  benchmark_dir <- paste("/path/to/prosstt/benchmark", b, "/", sep="")
  dim_keep <- b + 1
  for (i in 0:99) {
    JobName <- paste("sim", i, sep="")
    JobFolder <- paste(benchmark_dir, JobName, "/", sep="")
    job <- paste(JobFolder, JobName, sep="")
    cellparams <- read.table(file = paste(job, "cellparams.txt", sep = "_"), sep = "\t", header = T, row.names = 1)
    cellparams$branches <- cellparams$branches + 1
    cellparams$pseudotime <- cellparams$pseudotime + 1
    topology <- read_prosstt_topology(job) + 1
    
    dm_name <- paste(job, "destiny_log", sep = "_")
    dm <- readRDS(dm_name)
    data <- eigenvectors(dm)[, 1:dim_keep]
    
    branch_lengths <- branch_eucl_length(data, cellparams)
    
    g <- igraph::graph_from_edgelist(topology)
    root <- which(sapply(sapply(V(g), function(x) neighbors(g, x, mode="in")), length) == 0)
    
    paths <- igraph::all_simple_paths(g, root, V(g))
    for (k in seq_along(paths)) {
      p <- paths[[k]]
      direct_distance <- min(direct_dist(data, cellparams, p[1], tail(p, n = 1)))
      tree_distance <- mean(branch_lengths[[p]])
      prosstt_from_origin[[length(p)]] <- append(prosstt_from_origin[[length(p)]], direct_distance/tree_distance)
    }
    setTxtProgressBar(pb, (b-1) * 100 + (i+1))
  }
}

prosstt_dists <- lapply(prosstt_from_origin, unlist)
prosstt_howmany <- unlist(lapply(prosstt_dists, length))
prosstt_dists[[1]] <- c(1)
prosstt_where <- unlist(lapply(prosstt_dists, max))

splatter_dists <- lapply(splatter_from_origin, unlist)
splatter_howmany <- unlist(lapply(splatter_dists, length))
splatter_dists[[1]] <- c(1)
splatter_where <- unlist(lapply(splatter_dists, max))


svg("/save/path/prosstt_dm_dists.svg", width = 7, height = 5)
boxplot(prosstt_dists,
        main="Gene expression space",
        xlab="Time from origin (pseudotime units)",
        ylab="Normalized on-tree distance (average branch lengths)",
        xaxt="n", ylim=c(0., 5.0), xlim=c(0.5, 9.5))
box(which ="plot", lwd=2)
keep <- prosstt_where > -Inf
text(x=1:sum(keep), y=prosstt_where[keep] + 0.2, labels=prosstt_howmany[keep])
axis(1, at=1:9, labels = (1:9)*50)
abline(h = 1, lty=2)
dev.off()

svg("/save/path/splatter_dm_dists.svg", width = 7, height = 5)
boxplot(splatter_dists,
        main="Diffusion maps",
        xlab="Time from origin (pseudotime units)",
        ylab="Normalized on-tree distance (average branch lengths)",
        xaxt="n", ylim=c(0, 5.), xlim=c(0.5, 9.5))
keep <- splatter_where > -Inf
text(x=1:sum(keep), y=splatter_where[keep] + 0.2, labels=splatter_howmany[keep])
axis(1, at=1:9, labels = (1:9)*50)
abline(h = 1, lty=2)
box(which ="plot", lwd=2)
dev.off()