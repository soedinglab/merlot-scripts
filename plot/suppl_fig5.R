# single bifurcation
this <- "test5"
params_loc <- paste("path/to/benchmark1/", this, "/", this, "_cellparams.txt", sep="")
slicer_loc <- paste("path/to/benchmark1/", this, "/", this, "_SLICER_branches.csv", sep="")
monocl_loc <- paste("path/to/benchmark1/", this, "/", this, "_monocl2_branches.csv", sep="")
hhtree_loc <- paste("path/to/benchmark1/", this, "/", this, "_hhTree_log_k_fixed_emb_branches.csv", sep="")
hhfree_loc <- paste("path/to/benchmark1/", this, "/", this, "_hhTree_log_k_free_emb_sens_branches.csv", sep="")
montree_loc <- paste("path/to/benchmark1/", this, "/", this, "_monTree_fixed_emb_branches.csv", sep="")

labels <- list()
slicer <- read.csv(slicer_loc, row.names = 1)$x
hhtree <- read.csv(hhtree_loc, row.names = 1)$x
monocl <- read.csv(monocl_loc, row.names = 1)$x
hhfree <- read.csv(hhfree_loc, row.names = 1)$x
montree <- read.csv(montree_loc, row.names = 1)$x

labels[[1]] <- as.factor(slicer)
labels[[2]] <- as.factor(monocl)
labels[[3]] <- as.factor(hhtree)
labels[[4]] <- as.factor(hhfree)
labels[[5]] <- as.factor(montree)

names <- c("slicer", "monocl", "hhtree", "hhfree", "montree")

for (i in 1:5) {
svg(paste("path/out/single_", names[i], ".svg", sep=""), width = 8, height = 8)
  br_length = 50
  x1 <- seq(1-sqrt(2), 1, length.out = br_length)
  x2 <- seq(1, 2, length.out = br_length)
  y1 <- rep(1, br_length)
  y2 <- seq(1, 2, length.out = br_length)
  y3 <- seq(1, 0, length.out = br_length)
  x <- c(x1, x2, x2)
  y <- c(y1, y2, y3)
  plot(x, y, cex=7.5, pch="|", col=labels[[i]], axes=F, ann=F)
  dev.off()
}


# triple bifurcation
this <- "test30"
params_loc <- paste("path/to/benchmark3/", this, "/", this, "_cellparams.txt", sep="")
slicer_loc <- paste("path/to/benchmark3/", this, "/", this, "_SLICER_branches.csv", sep="")
monocl_loc <- paste("path/to/benchmark3/", this, "/", this, "_monocl2_branches.csv", sep="")
hhtree_loc <- paste("path/to/benchmark3/", this, "/", this, "_hhTree_log_k_fixed_emb_branches.csv", sep="")
hhfree_loc <- paste("path/to/benchmark3/", this, "/", this, "_hhTree_log_k_free_emb_sens_branches.csv", sep="")
montree_loc <- paste("path/to/benchmark3/", this, "/", this, "_monTree_fixed_emb_branches.csv", sep="")

labels <- list()
slicer <- read.csv(slicer_loc, row.names = 1)$x
hhtree <- read.csv(hhtree_loc, row.names = 1)$x
monocl <- read.csv(monocl_loc, row.names = 1)$x
hhfree <- read.csv(hhfree_loc, row.names = 1)$x
montree <- read.csv(montree_loc, row.names = 1)$x

labels[[1]] <- as.factor(slicer)
labels[[2]] <- as.factor(monocl)
labels[[3]] <- as.factor(hhtree)
labels[[4]] <- as.factor(hhfree)
labels[[5]] <- as.factor(montree)

names <- c("slicer", "monocl", "hhtree", "hhfree", "montree")

for (i in 1:5) {
  svg(paste("path/out/triple_", names[i], ".svg", sep=""), width = 8, height = 8)
  br_length = 50
  # topology: [[0, 1], [0, 2], [2, 3], [2, 4], [4, 5], [4, 6]]
  x1 <- seq(1-sqrt(2), 1, length.out = br_length)
  x2 <- seq(1, 2, length.out = br_length)
  x3 <- seq(2, 3, length.out = br_length)
  x4 <- seq(3, 4, length.out = br_length)
  y1 <- rep(1, br_length) # branch 0
  y2 <- seq(1, 2, length.out = br_length) # branch 1
  y3 <- seq(1, 0, length.out = br_length) # branch 2
  # [0, 2], [2, 3], [2, 4]
  y4 <- seq(0, -1, length.out = br_length) # branch 3
  y5 <- seq(0, 1, length.out = br_length) # branch 4
  # [2, 4], [4, 5], [4, 6]
  y6 <- seq(1, 2, length.out = br_length)
  y7 <- seq(1, 0, length.out = br_length)
  x <- c(x1, x2, x2, x3, x3, x4, x4)
  y <- c(y1, y2, y3, y4, y5, y6, y7)
  plot(x, y, cex=4, pch="|", col=labels[[i]], axes=F, ann=F)
  dev.off()
}


# quintuple bifurcation
this <- "test44"
params_loc <- paste("path/to/benchmark5/", this, "/", this, "_cellparams.txt", sep="")
slicer_loc <- paste("path/to/benchmark5/", this, "/", this, "_SLICER_branches.csv", sep="")
monocl_loc <- paste("path/to/benchmark5/", this, "/", this, "_monocl2_branches.csv", sep="")
hhtree_loc <- paste("path/to/benchmark5/", this, "/", this, "_hhTree_log_k_fixed_emb_branches.csv", sep="")
hhfree_loc <- paste("path/to/benchmark5/", this, "/", this, "_hhTree_log_k_free_emb_sens_branches.csv", sep="")
montree_loc <- paste("path/to/benchmark5/", this, "/", this, "_monTree_fixed_emb_branches.csv", sep="")

labels <- list()
slicer <- read.csv(slicer_loc, row.names = 1)$x
hhtree <- read.csv(hhtree_loc, row.names = 1)$x
monocl <- read.csv(monocl_loc, row.names = 1)$x
hhfree <- read.csv(hhfree_loc, row.names = 1)$x
montree <- read.csv(montree_loc, row.names = 1)$x

labels[[1]] <- as.factor(slicer)
labels[[2]] <- as.factor(monocl)
labels[[3]] <- as.factor(hhtree)
labels[[4]] <- as.factor(hhfree)
labels[[5]] <- as.factor(montree)

names <- c("slicer", "monocl", "hhtree", "hhfree", "montree")

for (i in 1:5) {
  svg(paste("path/out/quint_", names[i], ".svg", sep=""), width = 8, height = 8)
  br_length = 50
  # topology: [[0, 1], [0, 2], [2, 3], [2, 4], [3, 5], [3, 6], [6, 7], [6, 8], [7, 9], [7, 10]]
  x0 <- seq(1-sqrt(2), 1, length.out = br_length)
  x1 <- seq(1, 2, length.out = br_length)
  x2 <- seq(2, 3, length.out = br_length)
  r2 <- seq(3, 2, length.out = br_length)
  x3 <- seq(3, 4, length.out = br_length)
  x4 <- seq(4, 5, length.out = br_length)
  x5 <- seq(5, 6, length.out = br_length)
  # [0]
  y0 <- rep(1, br_length) # branch 0
  # [0, 1], [0, 2]
  y1 <- seq(1, 2, length.out = br_length) # branch 1
  y2 <- seq(1, 0, length.out = br_length) # branch 2
  # [1, 3], [1, 4]
  y3 <- seq(2, 3, length.out = br_length) # branch 3
  y4 <- seq(2, 1, length.out = br_length) # branch 4
  # [4, 5], [4, 6]
  y5 <- seq(1, 1.7, length.out = br_length)
  y6 <- seq(1, 0, length.out = br_length)
  # [3, 7], [3, 8]
  y7 <- seq(3, 4, length.out = br_length)
  y8 <- seq(3, 2.3, length.out = br_length)
  # [6, 9], [6, 10]
  y9 <- seq(0, 1, length.out = br_length)
  y10 <- seq(0, -1, length.out = br_length)
  x <- c(x0, x1, x1, x2, x2, x3, x3, x3, x3, x4, x4)
  y <- c(y0, y1, y2, y3, y4, y5, y6, y7, y8, y9, y10)
  plot(x, y, cex=4, pch="|", col=labels[[i]], axes=F, ann=F)
  dev.off()
}
