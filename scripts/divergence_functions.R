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

read_prosstt_topology <- function(job) {
  params <- suppressMessages(read_delim(paste(job, "_params.txt", sep=""), ":", col_names = FALSE, trim_ws = TRUE))
  top_index <- which(params$X1 == "topology")
  group <- strsplit(params$X2[top_index], "],")[[1]]
  processed <- gsub("\\[|\\]|\\s", "", group, perl = TRUE)
  aslist <- strsplit(processed, ",")
  topology <- matrix(0, ncol = 2, nrow = length(aslist))
  for (i in seq_along(aslist)) {
    topology[i, ] <- as.numeric(aslist[[i]])
  }
  topology
}