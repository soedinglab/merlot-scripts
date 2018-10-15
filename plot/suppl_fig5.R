library(prosstt)

pal <- palette(c("black", "red", "green3", "blue", "cyan", "magenta", "yellow", "gray", "darkorange", "saddlebrown",
                 "chartreuse", "aquamarine4", "purple3", "gray", "lightpink", "steelblue1"))

mscripts <- "path/to/merlot-scripts/"
source(paste(mscripts, "plot/auxiliary.R", sep=""))

selected <- prosstt_selected

# plot a single bifurcation
this <- "sim22"
benchmark_dir <- "path/to/splatter/benchmark1/"
params_loc <- paste(benchmark_dir, this, "/", this, "_params.txt", sep="")
scores_loc <- paste(benchmark_dir, this, "/", this, "_eval.txt", sep="")
job <- paste(benchmark_dir, this, "/", this, sep="")
params <- read.table(file = paste(job, "cellparams.txt", sep = "_"), sep = "\t", header = T, row.names = 1)

blarb <- flat_simulation(params, params_loc, mode = "splatter")
scores <- read.table(scores_loc, check.names = FALSE)

# svg("/save/location/bif1.svg", width = 10, height = 8)
par(mfrow=c(3,4))
plot_flat_tree(params, blarb, params$branches+1, plot_title = "truth", pcex = 3, col_pal = pal)
plot_flat_tree(params, blarb, params$pseudotime, plot_title = "truth", pcex = 3)
for (m in methnames[selected]) {
  tryCatch(expr = {
    pred_loc <- paste(job, m, "branches.csv", sep="_")
    pred <- read.csv(pred_loc, row.names = 1)
    plot_flat_tree(params, blarb, pred$x, plot_title = m, pcex = 3, col_pal = pal)
    legend("topleft", legend=round(scores[m, ]$`adjusted MI`, 3), bty="n")
  }, error = function(cond) {
    plot(1, type="n", axes=F, xlab="", ylab="", main = m)
  })
}
# dev.off()

# plot a triple bifurcation
this <- "sim27"
benchmark_dir <- "path/to/slatter/benchmark3/"
params_loc <- paste(benchmark_dir, this, "/", this, "_params.txt", sep="")
scores_loc <- paste(benchmark_dir, this, "/", this, "_eval.txt", sep="")
job <- paste(benchmark_dir, this, "/", this, sep="")
params <- read.table(file = paste(job, "cellparams.txt", sep = "_"), sep = "\t", header = T, row.names = 1)

blarb <- flat_simulation(params, params_loc, mode = "splatter")
scores <- read.table(scores_loc, check.names = FALSE)

# we need more colors here because the TSCAN prediction creates 33 clusters
cols <- c("black", "red", "green3", "blue", "cyan", "magenta", "yellow", "gray", "darkorange", "saddlebrown",
          "chartreuse", "aquamarine4", "purple3", "gray", "lightpink", "steelblue1")
fill <- LSD::distinctcolors(nrcol = 20, bw=TRUE, show=FALSE)
pal <- palette(c(cols, fill))

# svg("/save/locations/bif3.svg", width = 10, height = 8)
par(mfrow=c(3,4))
plot_flat_tree(params, blarb, params$branches+1, plot_title = "truth", pcex = 2, col_pal = pal)
plot_flat_tree(params, blarb, params$pseudotime, plot_title = "truth", pcex = 2)
for (m in methnames[selected]) {
  tryCatch(expr = {
    pred_loc <- paste(job, m, "branches.csv", sep="_")
    pred <- read.csv(pred_loc, row.names = 1)
    plot_flat_tree(params, blarb, pred$x, plot_title = m, pcex = 2, col_pal = pal)
    legend("topleft", legend=round(scores[m, ]$`adjusted MI`, 3), bty="n")
  }, error = function(cond) {
    plot(1, type="n", axes=F, xlab="", ylab="", main = m)
  })
}
# dev.off()

# plot a five-fold bifurcation
this <- "sim31"
benchmark_dir <- "path/to/splatter/benchmark5/"
params_loc <- paste(benchmark_dir, this, "/", this, "_params.txt", sep="")
scores_loc <- paste(benchmark_dir, this, "/", this, "_eval.txt", sep="")
job <- paste(benchmark_dir, this, "/", this, sep="")
params <- read.table(file = paste(job, "cellparams.txt", sep = "_"), sep = "\t", header = T, row.names = 1)

blarb <- flat_simulation(params, params_loc, mode = "splatter")
scores <- read.table(scores_loc, check.names = FALSE)
# revert to original color palette
pal <- palette(cols)

# svg("/save/location/bif5.svg", width = 10, height = 8)
par(mfrow=c(3,4))
plot_flat_tree(params, blarb, params$branches+1, plot_title = "truth", pcex = 2, col_pal = pal)
plot_flat_tree(params, blarb, params$pseudotime, plot_title = "truth", pcex = 2)
for (m in methnames[selected]) {
  tryCatch(expr = {
    pred_loc <- paste(job, m, "branches.csv", sep="_")
    pred <- read.csv(pred_loc, row.names = 1)
    plot_flat_tree(params, blarb, pred$x, plot_title = m, pcex = 2, col_pal = pal)
    legend("topleft", legend=round(scores[m, ]$`adjusted MI`, 3), bty="n")
  }, error = function(cond) {
    plot(1, type="n", axes=F, xlab="", ylab="", main = m)
  })
}
# dev.off()