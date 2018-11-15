library(Rtsne)
library(viridis)
library(LSD)
source("/path/to/merlot-scipts/scripts/splat_funcs.R")

set.seed(42)

i <- 9
j <- 28

main_dir <- "path/to/benchmark"
out_dir <- "/output/path"
benchmark_dir = paste(main_dir, i, "/", sep="")
this <- paste("sim", j, sep="")
job <- paste(benchmark_dir, this, "/", this, sep="")
sim_name <- paste(job, "simulation.txt", sep="_")

X <- as.matrix(read.table(file = sim_name, sep = "\t", header = T, row.names = 1, stringsAsFactors = T))
sum1 <- apply(X, 1, sum)
scalings <- sum1/mean(sum1)
X = (1/scalings)*X
subsample <- sample(x = 1:dim(X)[1], size=dim(X)[1]/2, replace = FALSE)
# subsample <- 1:dim(X)[1]
X <- X[subsample, ]

params <- read.table(file = paste(job, "cellparams.txt", sep = "_"), sep = "\t", header = T, row.names = 1)
labels <- params$branches + 1
time_cols <- viridis::viridis(100)
pseudotime <- map2color(params$pseudotime, viridis(100))

num_col <- length(unique(labels))
branch_cols <- LSD::distinctcolors(nrcol=num_col, method = "RGB", bw=TRUE)

pca_test <- prcomp(log(X+1))
# plot(pca_test$x[, 1:2], pch=16, col=as.factor(params$branches[subsample]), main=paste("splatter", this))
tsne_test <- Rtsne(pca_test$x[, 1:20], initial_dims=20, perplexity = 15)
# plot(tsne_test$Y, pch=16, col=branch_cols[as.factor(labels[subsample])], main=paste("prosstt", "bif", i, "sim", j))

svg(paste(out_dir, this, ".svg", sep=""))
plot(tsne_test$Y, pch=16, col=branch_cols[as.factor(labels[subsample])], main=paste("splatter", "bif", i, "sim", j))
dev.off()

svg(paste(out_dir, this, "pseudotime.svg", sep=""))
plot(tsne_test$Y, pch=16, col=pseudotime[subsample], main=paste("splatter", "bif", i, "sim", j))
dev.off()
