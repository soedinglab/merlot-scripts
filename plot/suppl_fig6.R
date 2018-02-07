#!/usr/bin/Rscript

#!/usr/bin/Rscript

library(ggplot2)
library(reshape2)
library(scales)

# import auxiliary functions for I/O and parsing
mscripts <- "path/to/merlot-scripts"

various <- paste(mscripts, "/various.R", sep="")
evaluat <- paste(mscripts, "/evaluate_method.R", sep="")

suppressPackageStartupMessages(source(various))
suppressPackageStartupMessages(source(evaluat))

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

main_dir <- "path/to/benchmark"
folder_name <- "" # name of folder with each set. For published benchmarks: "benchmark"
main_dir <- paste(main_dir, folder_name, sep="")
bif_nums <- 1:10
bif_names <- c("single", "double", "triple", "quadruple", "quintuple",
              "sextuple", "septuple", "octuple", "nonuple", "decuple")

methnames <- c("destiny_log",
               "SLICER",
               "monocl2",
               "hhTree_log_k_fixed_emb",
               "monTree_fixed_emb",
               "hhTree_log_k_free_emb_sens",
               "monTree_unc_free_emb_sens")

timenames <- c("goodman-kruskall (unweighted)", "goodman-kruskall (weighted)",
               "kendall index (unweighted)", "kendall index (weighted)")
shortnames <- c("gku", "gkw", "kiu", "kiw")
leg_names <- c("Goodman-Kruskall (unweighted)", "Goodman-Kruskall (weighted)",
               "Kendall index (unweighted)", "Kendall index (weighted)")
res <- array(0, dim = c(length(methnames), length(bif_nums), length(timenames)))
dimnames(res) <- list(methnames, bif_names, timenames)
err <- array(0, dim = c(length(methnames), length(bif_nums), length(timenames)))

for (b in 1:10) {
  tmp <- array(0, dim = c(length(methnames), length(timenames), 100))

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

  dimnames(tmp) <- list(methnames, timenames, subfolders)

  for (i in 1:100) {
    s <- subfolders[i]
    job <- paste(benchmark_dir, s, "/", s, sep = "")
    par_loc <- paste(job, "params.txt", sep = "_")
    cell_params <- read.table(file = paste(job, "cellparams.txt", sep = "_"),
                              sep = "\t", header = T, row.names = 1)
    for (j in 1:length(methnames)) {
      m <- methnames[j]
      tryCatch({
        mtimes <- read.csv(paste(job, m, "times.csv", sep = "_"))$x
        row_res <- all_longest_paths(s, m, mtimes, cell_params, par_loc)
      },
      error = function(e) {
        print(paste("couldn't find", m, "for", s, "in benchmark set", b))
        row_res <- c(0, 0, 0, 0)
      },
      finally = {
        tmp[m, , s] <- as.numeric(row_res)
      })
    }
  }
  storage.mode(tmp) <- "numeric"
  tmp[which(is.nan(tmp))] <- 0
  # saveRDS(object = tmp, file = paste(benchmark_dir, "longest_time.csv",
          # sep = ""))
  res[, b, ] <- apply(tmp, c(1, 2), mean)

  q <- qnorm(0.975)
  err[, b, ] <- q * sqrt(apply(tmp, c(1, 2), var)) / 10
}

tups <- res + err
tdowns <- res - err

paper_cols <- c("darkorchid", "dodgerblue", "brown2", "darkgreen", "forestgreen", "green3", "lawngreen")
paper_pch <- c(18, 1, 2, 18, 0, 18, 0)
paper_lty <- c(1, 1, 1, 1, 1, 2, 2)

for (j in 1:length(timenames)) {
  svg(paste("~/Documents/collaborations/gonzalo/PaperTree/Figuras/Supplementary/pseudotime_benchmark/all_pseudotimes/lp_", shortnames[j], ".svg",
  sep = ""), height = 5, width = 8)
  plot(bif_nums, res[methnames[1],,timenames[j]],
       ylim=c(0, 0.8),
       pch=paper_pch[1],
       col=paper_cols[1],
       type="b",
       lwd=3,
       lty=paper_lty[j],
       cex=1.5,
       ylab=paste("Longest Tree", timenames[j]),
       xlab="#cell fates",
       axes=FALSE,
       main="Pseudotime assignment quality")
  arrows(x0=bif_nums,
         y0 = tdowns[methnames[1],,timenames[j]],
         y1=tups[methnames[1],,timenames[j]],
         code=3,
         angle=90,
         length=0.1,
         col=cols[1])
  
  box(which = "plot", lwd=2)
  axis(1, at=bif_nums, labels = bif_nums+2)
  axis(2, at=c(-1, -0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75, 1.), labels = c(-1, -0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75, 1))
  
  for (i in 2:length(methnames)) {
    m = methnames[i]
    mc = cols[i]
    points(bif_nums, res[methnames[i],,timenames[j]],
         ylim=c(0, 0.75),
         pch=paper_pch[i],
         col=paper_cols[i],
         type="b",
         lwd=3,
         lty=paper_lty[i],
         cex=1.5)
    arrows(x0=bif_nums,
           y0 = tdowns[methnames[i],,timenames[j]],
           y1=tups[methnames[i],,timenames[j]],
           code=3,
           angle=90,
           length=0.1,
           col=cols[i])
  }
  dev.off()
}

svg("~/Documents/collaborations/gonzalo/PaperTree/Figuras/Supplementary/pseudotime_benchmark/legend.svg", height=5, width=12)
plot.new()
permutation <- c(4, 6, 1, 2, 5, 7, 3)
legend_text <- c("Destiny", "SLICER", "Monocle2", "MERLoT + destiny fixed", "MERLoT + DDRTree fixed", "MERLoT + destiny free", "MERLoT + DDRTree free")
legend("bottom", legend=legend_text[permutation], lty=paper_lty[permutation], lwd=3,
       col=paper_cols[permutation], ncol=4, pt.cex=1.5, pch=paper_pch[permutation])
dev.off()