library(RColorBrewer)
# pal <- palette(c("black", "red", "green3", "blue", "cyan", "magenta", "yellow", "gray"))

mscripts <- "path/to/merlot-scripts/"
source(paste(mscripts, "scripts/flat_trees.R", sep=""))

pal <- brewer.pal("Set3", n=12)
methnames <- c("destiny_log_k",
               "destiny_log_k_sel_seurat",
               "destiny_log",
               "destiny_log_sel_seurat",
               "SLICER",
               "monocl2",
               "monocl2_unc",
               "monGraph_free_el_sens",
               "monGraph_unc_free_el_sens",
               "monGraph_free_emb_sens",
               "monGraph_unc_free_emb_sens",
               "monGraph_fixed_el",
               "monGraph_unc_fixed_el",
               "monGraph_fixed_emb",
               "monGraph_unc_fixed_emb",
               "LPGraph_log_k_sel_seurat_free_el_sens",
               "LPGraph_log_k_free_el_sens",
               "LPGraph_log_k_sel_seurat_free_emb_sens",
               "LPGraph_log_k_free_emb_sens",
               "LPGraph_log_k_sel_seurat_fixed_el",
               "LPGraph_log_k_fixed_el",
               "LPGraph_log_k_sel_seurat_fixed_emb",
               "LPGraph_log_k_fixed_emb",
               "TSCAN",
               "slingshot_destiny_log",
               "slingshot_monocl2",
               "slingshot_destiny_log_k",
               "slingshot_destiny_log_k_sel_seurat",
               "slingshot_destiny_log_sel_seurat")

selected <- c(5, 6, 22, 14, 18, 10, 29, 26, 24)

# plot a single bifurcation
this <- "sim36"
benchmark_dir <- "path/to/splatter/benchmark1/"
params_loc <- paste(benchmark_dir, this, "/", this, "_params.txt", sep="")
scores_loc <- paste(benchmark_dir, this, "/", this, "_eval.txt", sep="")
job <- paste(benchmark_dir, this, "/", this, sep="")
params <- read.table(file = paste(job, "cellparams.txt", sep = "_"), sep = "\t", header = T, row.names = 1)

blarb <- flat_simulation(params, params_loc, mode = "splatter")
scores <- read.table(scores_loc, check.names = FALSE)

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

# plot a triple bifurcation
this <- "sim52"
benchmark_dir <- "path/to/splatter/benchmark3/"
params_loc <- paste(benchmark_dir, this, "/", this, "_params.txt", sep="")
scores_loc <- paste(benchmark_dir, this, "/", this, "_eval.txt", sep="")
job <- paste(benchmark_dir, this, "/", this, sep="")
params <- read.table(file = paste(job, "cellparams.txt", sep = "_"), sep = "\t", header = T, row.names = 1)

blarb <- flat_simulation(params, params_loc, mode = "splatter")
scores <- read.table(scores_loc, check.names = FALSE)

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

# plot a five-fold bifurcation
this <- "sim97"
benchmark_dir <- "path/to/splatter/benchmark5/"
params_loc <- paste(benchmark_dir, this, "/", this, "_params.txt", sep="")
scores_loc <- paste(benchmark_dir, this, "/", this, "_eval.txt", sep="")
job <- paste(benchmark_dir, this, "/", this, sep="")
params <- read.table(file = paste(job, "cellparams.txt", sep = "_"), sep = "\t", header = T, row.names = 1)

blarb <- flat_simulation(params, params_loc, mode = "splatter")
scores <- read.table(scores_loc, check.names = FALSE)

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
