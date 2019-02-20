source("/path/to/merlot-scipts/scripts/auxiliary.R")

bench_dir <- "/path/to/benchmark"
selected <- prosstt_selected

for (i in 1:length(bif_num)) {
  benchmark_dir <- paste(bench_dir, i, "/", sep="")
  
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
  br_res = br_res[br_all$measure == "#branches",]
  
  # calculate confidence intervals for error bars
  q = qnorm(0.975)
  k <- which(res$X == bif_num[i])
  for (m in methnames) {
    res[k, m] = mean(br_res[,m][br_res[,m] > 0])
    err[k, m] = q * sqrt(var(br_res[,m])) / 10
  }
}

# (symmetric) error bars
ups = res + err
downs = res - err

par(mfrow=c(1,1))
svg("/save/location", height=5, width=8)
plot(bif_num, res[bif_num, methnames[4]], ylim=c(0, 26), col=cols[4], type="n", lwd=3, pch=pchs[4], cex=1.5,
     ylab="#predicted branches", xlab="#cell fates", axes=FALSE, main="Predicted number of branches")
box(which = "plot", lwd=2)
axis(1, at=bif_num, labels = bif_num+2)
axis(2, at=seq(0, 28, 5), labels = seq(0, 28, 5))

points(bif_num,
       bif_num*2 + 1,
       type="b", lwd = 3, col="gray")

for (i in seq_along(selected)) {
  s <- selected[i]
  m = methnames[s]
  mc = cols[i]
  points(bif_num, res[bif_num, m], pch=pchs[i], col=mc, type="b", lwd=2, cex=1.5, lty=ltys[i])
  arrows(x0=bif_num, y0 = downs[bif_num, m], y1=ups[bif_num, m],
         code=3, angle=90, length=0.1, col=mc)
}
dev.off()

# legend
svg("/save/location/", height=5, width=12)
plot.new()
permutation <- c(11, 4, 6, 1, 2, 5, 7, 3, 8, 9, 10)
legend_text <- c("Destiny", "SLICER", "Monocle2", "MERLoT + destiny fixed", "MERLoT + DDRTree fixed",
                 "MERLoT + destiny auto", "MERLoT + DDRTree auto", "slingshot + destiny", "slingshot + monocle",
                 "TSCAN", "Expected branch number")
paper_cols <- c("royalblue4", cols, "gray")
paper_pch <- c(18, pchs, 1)
paper_lty <- c(1, ltys, 1)
legend("bottom", legend=legend_text[permutation], lty=paper_lty[permutation], lwd=3,
       col=paper_cols[permutation], ncol=2, pt.cex=1.5, pch=paper_pch[permutation])
dev.off()
