suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(scales))
suppressPackageStartupMessages(library(plotly))

.pardefault <- par(no.readonly = T)

# find the longest path in a tree and score it according to a method
all_longest_paths <- function(name, method, mtimes, cell_params, par_loc,
                              returnNA = FALSE) {
  timefuncs <- list()
  timefuncs[[1]] <- goodman_kruskal_index
  timefuncs[[2]] <- goodman_kruskal_index
  timefuncs[[3]] <- kendall_index
  timefuncs[[4]] <- kendall_index
  
  timenames <- c("goodman-kruskall (unweighted)", "goodman-kruskall (weighted)",
                 "kendall index (unweighted)", "kendall index (weighted)")
  
  res <- data.frame(matrix(NA, nrow = length(timenames), ncol = 1))
  rownames(res) <- timenames
  colnames(res) <- name
  
  if (returnNA){
    return(t(res))
  }
  
  ltimes <- cell_params$pseudotime - min(cell_params$pseudotime) + 1
  # longest path goodman kruskal:
  lp_id <- get_lpgk_indices(par_loc, cell_params)
  lp_cells <- rownames(cell_params)[lp_id]
  ltimes <- cell_params[lp_cells, ]$pseudotime
  mtimes <- mtimes[lp_id]
  
  for (i in 1:length(timenames)) {
    weighted <- (i %% 2 == 0)
    res[i, 1] <- timefuncs[[i]](mtimes, ltimes, weighted = weighted)
  }
  return(t(res))
}

hhtree <- "~/Documents/repos/hhtree"
various <- paste(hhtree, "/scripts/various.R", sep="")
evaluat <- paste(hhtree, "/scripts/evaluate_method.R", sep="")

suppressPackageStartupMessages(source(various))
suppressPackageStartupMessages(source(evaluat))

bif_names = c("single", "double", "triple", "quadruple", "quintuple", "sextuple", "septuple", "octuple", "nonuple", "decuple")

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


# structures that will aggregate the results
# branch evaluation results
bif_num <- 1:10
res = data.frame(matrix(data=0, ncol=length(methnames)+1, nrow=length(bif_num)))
colnames(res) = c("X", methnames)
res$X = bif_num

# variance of branch evaluation
err = data.frame(matrix(data=0, ncol=length(methnames)+1, nrow=length(bif_num)))
colnames(err) = c("X", methnames)
err$X = bif_num

# functions for branch detection evaluation, names to be used in legends, and short names for generated files
funcnames <- c("rand index", "matthews corr. coef.",
               "F1 measure", "jaccard index", "fowkles-mallows index", "adjusted MI", "#branches")
leg_names <- c("Rand Index", "Matthews Corr. Coef.", "F1 Measure", "Jaccard index", "Fowkles-Mallows index", "Normalized MI")
shortnames <- c("rand", "mcc", "f1", "jaccard", "fowklesm", "adj_mi")

timenames <- c("goodman-kruskall (unweighted)", "goodman-kruskall (weighted)",
               "kendall index (unweighted)", "kendall index (weighted)")

ylims <- list(c(0,1), c(-1., 0.75), c(0, 0.8), c(0, 0.75), c(0, 0.8), c(0, 0.75))

# plot MI accross bif number
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

# selected methods to plot
# in this order: SLICER, monocle2, merlot+dest fixed, merlot+monocl fixed, merlot+dest free, merlot+mon fixed, slingshot+destiny, slingshot+monocle, TSCAN
splatter_selected <- c(5, 6, 22, 14, 18, 10, 29, 26, 24)
prosstt_selected <-  c(5, 6, 23, 12, 19, 8, 27, 26, 24)
# colors for the paper plots
cols <- c("dodgerblue", "brown2", "darkgreen", "forestgreen", "green3", "lawngreen",
          "mediumorchid1", "mediumorchid4", "darkorange")
pchs <- c(1, 2, 18, 0, 18, 0, 6, 6, 4)
ltys <- c(1, 1, 1, 1, 2, 2, 1, 1, 1)

# legend colors in case we are doing the paper plots
legend_text <- c("SLICER",
                 "Monocle2",
                 "MERLoT + destiny fixed",
                 "MERLoT + DDRTree fixed",
                 "MERLoT + destiny auto",
                 "MERLoT + DDRTree auto",
                 "slingshot + destiny",
                 "slingshot + DDRTree",
                 "TSCAN + PCA")