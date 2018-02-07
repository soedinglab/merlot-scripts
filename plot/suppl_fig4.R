library(ggplot2)
library(reshape2)
library(scales)
library(plotly)

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

# import auxiliary functions for I/O and parsing
mscripts <- "~/Documents/repos/merlot-scripts"

various <- paste(mscripts, "/various.R", sep="")
evaluat <- paste(mscripts, "/evaluate_method.R", sep="")

suppressPackageStartupMessages(source(various))
suppressPackageStartupMessages(source(evaluat))

# benchmark directory
main_dir <- "/data/niko/repeat/benchmark"
bif_num = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)

# the methods that were run in the benchmark
methnames <- c("destiny_log",
               "destiny_log_k",
               "SLICER",
               "monocl2",
               "monocl2_unc",
               "hhTree_log_free_el_sens",
               "hhTree_log_k_free_el_sens",
               "hhTree_log_k_free_emb_sens",
               "hhTree_log_k_fixed_el",
               "hhTree_log_k_fixed_emb",
               "monTree_free_el_sens", "monTree_unc_free_el_sens",
               "monTree_free_emb_sens", "monTree_unc_free_emb_sens",
               "monTree_fixed_el", "monTree_unc_fixed_el",
               "monTree_fixed_emb", "monTree_unc_fixed_emb")

# the methods that are plotted here
selected <- c(3,4,10,17,8,14)

# structures that will aggregate the results
# branch evaluation results
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

ylims <- list(c(0,1), c(-1., 0.6), c(0, 0.8), c(0, 0.6), c(0, 0.8), c(0, 0.75))

# colors to plot all methods
cols = gg_color_hue(length(methnames))
pchs <- rep(20, length(methnames))
ltys <- rep(1, length(methnames))

# colors for the paper plots
paper_cols <- c("dodgerblue", "brown2", "darkgreen", "forestgreen", "green3", "lawngreen")
paper_pch <- c(1, 2, 18, 0, 18, 0)
paper_lty <- c(1, 1, 1, 1, 2, 2)
cols[1] <- "darkorchid" # color for destiny
pchs[1] <- 18
cols[selected] <- paper_cols
pchs[selected] <- paper_pch
ltys[selected] <- paper_lty

# legend colors in case we are doing the paper plots
legend_text <- c("SLICER", "Monocle2", "MERLoT + destiny fixed", "MERLoT + DDRTree fixed", "MERLoT + destiny free", "MERLoT + DDRTree free")
legend_col <- cols[selected]
legend_pch <- pchs[selected]
legend_lty <- ltys[selected]
permutation <- c(4, 6, 3, 5, 1, 2) # different order than the plotting

for (f in 1:6) {
  chosen <- funcnames[f]
  for (i in 1:length(bif_num)) {
    benchmark_dir = paste(main_dir, bif_num[i], "/", sep="")
    
    # remove missing files
    subfolders <- dir(benchmark_dir)
    remove <- c()
    for (s in subfolders) {
      eval <- paste(benchmark_dir, s, "/", s, "_eval.txt", sep="")
      if (!file.exists(eval)) {
        remove <- cbind(remove, s)
      }
    }
    
    subfolders <- subfolders[!is.element(subfolders, remove)]
    
    # parse the benchmark files and convert everything to numerics. Choose the relevant method.
    br_all <- parse_benchmark(subfolders, benchmark_dir, funcnames, methnames, "branch", replNA = TRUE)
    br_all <- as.data.frame(br_all)
    br_res <- br_all[,1:length(methnames)]
    br_res <- apply(br_res, 2, as.numeric)
    br_res = br_res[br_all$measure == chosen,]
    
    # calculate confidence intervals for error bars
    q = qnorm(0.975)
    k <- which(res$X == bif_num[i])
    for (m in methnames) {
      res[k, m] = mean(br_res[,m])
      err[k, m] = q * sqrt(var(br_res[,m])) / 10
    }
  }
  
  # (symmetric) error bars
  ups = res + err
  downs = res - err
  
  # svg(paste("some/directory/", shortnames[f], ".svg", sep=""), height=5, width=8)
  plot(bif_num, res[bif_num, methnames[4]], ylim=ylims[[f]], col=cols[4], type="n", lwd=3, pch=pchs[4], cex=1.5,
       ylab=leg_names[f], xlab="#cell fates", axes=FALSE, main="Branch assignment quality")
  
  box(which = "plot", lwd=2)
  axis(1, at=bif_num, labels = bif_num+2)
  axis(2, at=c(-1, -0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75, 1.), labels = c(-1, -0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75, 1))
  
  for (i in selected) {
    m = methnames[i]
    mc = cols[i]
    points(bif_num, res[bif_num, m], pch=pchs[i], col=mc, type="b", lwd=3, cex=1.5, lty=ltys[i])
    arrows(x0=bif_num, y0 = downs[bif_num, m], y1=ups[bif_num, m],
           code=3, angle=90, length=0.1, col=mc)
  }
  legend("bottom", legend=legend_text[permutation], lty=legend_lty[permutation], lwd=3,
         col=legend_col[permutation], ncol=3, pt.cex=1.5, pch=legend_pch[permutation])
  # dev.off()
}

# how to produce the legend used in the paper
# svg("some/directory", height=5, width=8)
plot.new()
selected <- c(1, 3, 4, 10, 17, 8, 14)
permutation <- c(4, 6, 1, 2, 5, 7, 3, 8)
legend_text <- c("Destiny", "SLICER", "Monocle2", "MERLoT + destiny fixed", "MERLoT + DDRTree fixed",
                 "MERLoT + destiny auto", "MERLoT + DDRTree auto")
legend_col <- cols[selected]
legend_pch <- pchs[selected]
legend_lty <- ltys[selected]
legend("bottom", legend=legend_text[permutation], lty=legend_lty[permutation], lwd=3,
       col=legend_col[permutation], ncol=2, pt.cex=1.5, pch=legend_pch[permutation])
# dev.off()
