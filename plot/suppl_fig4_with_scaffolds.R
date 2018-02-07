#!/usr/bin/Rscript

library(ggplot2)
library(reshape2)
library(scales)

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

# import auxiliary functions for I/O and parsing
mscripts <- "path/to/merlot-scripts"

various <- paste(mscripts, "/various.R", sep="")
evaluat <- paste(mscripts, "/evaluate_method.R", sep="")

suppressPackageStartupMessages(source(various))
suppressPackageStartupMessages(source(evaluat))

# read branch assignments from scaffold tree so we can
# easily evaluate performance
branchify <- function(scaffold) {
  N <- dim(scaffold$CellCoordinates)[1]
  Cells2Branches <- rep(0, N)
  for(i in 1:length(scaffold$Cells2BranchesAssignments)) {
    Cells2Branches[scaffold$Cells2BranchesAssignments[[i]]]=i
  }
  return(Cells2Branches)
}

# retrieve all information for the scaffold trees
all_scaffolds <- function(name, method, mbranches, cell_params, par_loc,
                              returnNA = FALSE) {
  res <- data.frame(matrix(NA, nrow = length(funcnames), ncol = 1))
  rownames(res) <- funcnames
  colnames(res) <- name
  
  if (returnNA){
    return(t(res))
  }
  
  lbranches <- cell_params$branches + 1
  status <- assign_status(mbranches, lbranches)
  
  for (i in 1:length(funcnames)) {
    # weighted <- (i %% 2 == 0)
    res[i, 1] <- funclist[[i]](mbranches, lbranches, status)
  }
  return(t(res))
}

main_dir <- "benchmark/directory"
bif_nums <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)
bif_names <- c("single", "double", "triple", "quadruple", "quintuple",
               "sextuple", "septuple", "octuple", "nonuple", "decuple")

# scaffold methnames
methnames <- c("hhTree_log_k_fixed_el",
               "hhTree_log_k_free_el_sens",
               "monTree_fixed_el",
               "monTree_free_el_sens",
               "monTree_unc_fixed_el",
               "monTree_unc_free_el_sens")

# names of all the methods that were run in the benchmark
all_methnames <- c("destiny_log",
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
selected <- c(3, 4, 10, 17, 8, 14)
scaf_sel <- c(0, 0, 1, 3, 2, 6)

all_funcnames <- c("rand index", "matthews corr. coef.",
               "F1 measure", "jaccard index", "fowkles-mallows index", "adjusted MI", "#branches")
ylims <- list(c(0,1), c(-1., 0.75), c(0, 0.8), c(0, 0.75), c(0, 0.8), c(0, 0.75))

funclist <- list()
funclist[[1]] <- randInd_manual
funclist[[2]] <- matthews_cor
funclist[[3]] <- f_measure
funclist[[4]] <- jaccard
funclist[[5]] <- fowkles_mallows
funclist[[6]] <- adjusted_mi

funcnames <- c("rand index", "matthews corr. coef.",
               "F1 measure", "jaccard index", "fowkles-mallows index", "adjusted MI")
leg_names <- c("Rand Index", "Matthews Corr. Coef.", "F1 Measure", "Jaccard index",
               "Fowkles-Mallows index", "Normalized MI")
shortnames <- c("rand", "mcc", "f1", "jaccard", "fowklesm", "adj_mi")

res <- array(0, dim = c(length(selected), length(bif_nums), length(funcnames)))
dimnames(res) <- list(all_methnames[selected], bif_names, funcnames)
err <- array(0, dim = c(length(selected), length(bif_nums), length(funcnames)))

scaf_res <- array(0, dim = c(length(methnames), length(bif_nums), length(funcnames)))
dimnames(scaf_res) <- list(methnames, bif_names, funcnames)
scaf_err <- array(0, dim = c(length(methnames), length(bif_nums), length(funcnames)))


for (b in 1:10) {
  tmp <- array(0, dim = c(length(methnames), length(funcnames), 100))
  
  print(paste(b, "/", length(bif_nums)))
  benchmark_dir <- paste(main_dir, b, "/", sep = "")
  bif_num <- bif_names[b]
  subfolders <- dir(benchmark_dir)
  remove <- c()
  
  for (s in subfolders) {
    eval <- paste(benchmark_dir, s, "/", s, "_eval.txt", sep = "")
    
    if (!file.exists(eval)) {
      remove <- cbind(remove, s)
    }
  }
  
  subfolders <- subfolders[!is.element(subfolders, remove)]
  print(paste(length(subfolders), " / 100", sep = ""))
  
  dimnames(tmp) <- list(methnames, funcnames, subfolders)
  
  for (i in 1:100) {
    s <- subfolders[i]
    job <- paste(benchmark_dir, s, "/", s, sep = "")
    par_loc <- paste(job, "params.txt", sep = "_")
    cell_params <- read.table(file = paste(job, "cellparams.txt", sep = "_"),
                              sep = "\t", header = T, row.names = 1)
    for (j in 1:length(methnames)) {
      m <- methnames[j]
      tryCatch({
        scaf_name <- paste(job, "_", m, ".scaf", sep="")
        scaf <- readRDS(scaf_name)
        mbranches <- branchify(scaf)
        row_res <- all_scaffolds(s, m, mbranches, cell_params, par_loc)
      },
      error = function(e) {
        print(paste("couldn't find", m, "for", s, "in benchmark set", b))
        row_res <- c(0, -1, 0, 0, 0, 0)
      },
      finally = {
        tmp[m, , s] <- as.numeric(row_res)
      })
    }
  }
  storage.mode(tmp) <- "numeric"
  tmp[which(is.nan(tmp))] <- 0
  # saveRDS(object = tmp, file = paste(benchmark_dir, "all_scaffolds.csv",
  #                                    sep = ""))
  scaf_res[, b, ] <- apply(tmp, c(1, 2), mean)
  
  
  br_all <- parse_benchmark(subfolders, benchmark_dir, all_funcnames, all_methnames, "branch", replNA = TRUE)
  br_all <- as.data.frame(br_all)
  br_num <- br_all[,1:length(all_methnames)]
  br_num <- apply(br_num, 2, as.numeric)
  q <- qnorm(0.975)
  for (f in 1:6) {
    chosen <- funcnames[f]
    br_res = br_num[br_all$measure == chosen,]
    res[, b, chosen] <- apply(br_res[, selected], 2, mean)
    err[, b, f] = q * sqrt(apply(br_res[, selected], 2, var)) / 10
  }
  scaf_err[, b, ] <- q * sqrt(apply(tmp, c(1, 2), var)) / 10
}

scaf_ups <- scaf_res + scaf_err
scaf_downs <- scaf_res - scaf_err

ups <- res + err
downs <- res - err

# plotting
# colors to plot all methods
cols = gg_color_hue(length(all_methnames))
pchs <- rep(20, length(all_methnames))
ltys <- rep(1, length(all_methnames))
paper_cols <- c("dodgerblue", "brown2", "darkgreen", "forestgreen", "green3", "lawngreen")
paper_pch <- c(1, 2, 18, 0, 18, 0)
paper_lty <- c(1, 1, 1, 1, 5, 5)

scaf_cols <- rep("red", length(methnames))
scaf_cols[scaf_sel] <- c("darkgreen", "forestgreen", "green3", "lawngreen")
scaf_pch <- rep(1, length(methnames))
scaf_pch[scaf_sel] <- c(18, 0, 18, 0)

cols[1] <- "darkorchid" # color for destiny
pchs[1] <- 18
cols[selected] <- paper_cols
pchs[selected] <- paper_pch
ltys[selected] <- paper_lty

# legend colors in case we are doing the paper plots

for (f in 1:6) {
  svg(paste("~/Documents/collaborations/gonzalo/all_branch/scaf_", shortnames[f], ".svg", sep=""), height=5, width=8)
  plot(bif_nums, res[all_methnames[4], bif_nums, f], ylim=ylims[[f]], col=cols[4], type="n", lwd=3, pch=pchs[4], cex=1.5,
       ylab=leg_names[f], xlab="#cell fates", axes=FALSE, main="Branch assignment quality")
  
  box(which = "plot", lwd=2)
  axis(1, at=bif_nums, labels = bif_nums+2)
  axis(2, at=c(-1, -0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75, 1.), labels = c(-1, -0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75, 1))
  
  for (i in selected) {
    m = all_methnames[i]
    mc = cols[i]
    points(bif_nums, res[m, bif_nums, f], pch=pchs[i], col=mc, type="b", lwd=3, cex=1.5, lty=ltys[i])
    arrows(x0=bif_nums, y0 = downs[m, bif_nums, f], y1=ups[m, bif_nums, f],
           code=3, angle=90, length=0.1, col=mc)
  }
  
  for (i in scaf_sel[scaf_sel!=0]) {
    m = methnames[i]
    mc = scaf_cols[i]
    points(bif_num, scaf_res[m, bif_num, f], pch=scaf_pch[i], col=mc, type="b", lwd=3, cex=1.5, lty=3)
    arrows(x0=bif_nums, y0 = scaf_downs[m, bif_nums, f], y1=scaf_ups[m, bif_nums, f],
           code=3, angle=90, length=0.1, col=mc)
  }
  dev.off()
}

# legend separately to make connecting the SVGs easier
# svg("path/out/scaf_legend.svg", height=3, width=12)
plot.new()

permutation <- c(3, 5, 4, 6, 1, 2, 7)
legend_text <- c("SLICER", "Monocle2", "MERLoT + destiny fixed",
                 "MERLoT + DDRTree fixed", "MERLoT + destiny auto", "MERLoT + DDRTree auto", "scaffold")
legend_col <- c(paper_cols, "black")
legend_lty <- c(1, 1, 1, 1, 5, 5, 3)
legend_pch <- c(pchs[selected], NA)
legend("bottom", legend=legend_text[permutation], lty=legend_lty[permutation], lwd=3,
       col=legend_col[permutation], ncol=4, pt.cex=1.5, pch=legend_pch[permutation])
dev.off()
