#!/usr/bin/Rscript
source("/path/to/merlot-scripts/scripts/auxiliary.R")

benchmarks <- c("/path/to/prosstt/benchmark",
                "/path/to/splatter/benchmark")
benchmark_names <- c("PROSSTT", "splatter")
selecteds <- list(prosstt_selected,
               prosstt_selected)

par(mfcol=c(4,2))
for (bs in seq_along(benchmarks)){
  bench_dir <- benchmarks[bs]
  selected <- selecteds[[bs]]
  for (f in c(6, 4, 3, 2)) {
    chosen <- funcnames[f]
    
    for (b in 1:length(bif_num)) {
      benchmark_dir <- paste(bench_dir, b, "/", sep="")
      
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
      # sim_names <- br_all$experiment[br_all$measure == chosen]
      # suspicious <- sim_names[which(br_res[,"destiny_log_k_sel_seurat"] == 0)]
      
      # calculate confidence intervals for error bars
      q = qnorm(0.975)
      k <- which(res$X == bif_num[b])
      for (m in methnames) {
        res[k, m] = mean(br_res[,m])
        err[k, m] = q * sqrt(var(br_res[,m])) / 10
      }
    }
    
    # (symmetric) error bars
    ups = res + err
    downs = res - err
    if(bs == 1) {
      par(mar=c(1,3,1,1))
    }
    
    if (f == 2) {
      par(mar=c(3,1,1,1))
      if (bs == 1) {
        par(mar=c(3,3,1,1))
      }
    }
    
    # plot branch assignment quality
    plot(bif_num, res[bif_num, methnames[4]], ylim=ylims[[f]], axes=FALSE, type="n")
    
    if(bs == 1) {
      mtext(leg_names[f], side=2, line=2)
      axis(2, at=c(-1, -0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75, 1.),
           labels = c(-1, -0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75, 1))
    }
    
    if(f == 6) {
      title(main=benchmark_names[bs])
      axis(1, at=bif_num, labels=FALSE)
      axis(2, at=c(-1, -0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75, 1.), labels = FALSE)
    } else if (f == 2) {
      axis(1, at=bif_num, labels = bif_num+2)
      mtext("#cell fates", side=1, line=2)
    } else {
      axis(1, at=bif_num, labels=FALSE)
      axis(2, at=c(-1, -0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75, 1.), labels = FALSE)
    }
    box(which = "plot", lwd=2)
    
    
    for (i in seq_along(selected)) {
      s <- selected[i]
      m = methnames[s]
      mc = cols[i]
      points(bif_num, res[, m], pch=pchs[i], col=mc, type="b", lwd=3, cex=1.5, lty=ltys[i])
      arrows(x0=bif_num, y0 = downs[bif_num, m], y1=ups[bif_num, m],
             code=3, angle=90, length=0.1, col=mc)
    }
    par(mar = c(1,1,1,1))
  }
}

# legend
par(mfrow=c(1,1))
plot(1, type="n", axes=F, ylab="", xlab="")
permutation <- c(1, 2, 9, 5, 3, 7, 6, 4, 8)
legend("bottom", lwd=3, ncol=3, pt.cex=1.5,
       legend=legend_text[permutation],
       lty=ltys[permutation],
       col=cols[permutation],
       pch=pchs[permutation])
