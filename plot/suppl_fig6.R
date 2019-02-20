#!/usr/bin/Rscript
source("/path/to/merlot-scripts/scripts/auxiliary.R")

benchmarks <- c("/path/to/prosstt/benchmark",
                "/path/to/splatter/benchmark")
benchmark_names <- c("PROSSTT", "splatter")
selecteds <- list(c(2, prosstt_selected),
                  c(1, prosstt_selected))

leg_names <- c("Goodman-Kruskall (unweighted)", "Goodman-Kruskall (weighted)",
               "Kendall index (unweighted)", "Kendall index (weighted)")

res_unit <- array(0, dim = c(10, length(bif_num), length(timenames)))
res <- list(res_unit, res_unit, res_unit)
err <- list(res_unit, res_unit, res_unit)
tups <- list(res_unit, res_unit, res_unit)
tdowns <- list(res_unit, res_unit, res_unit)

for (bs in seq_along(benchmarks)) {
  selected <- selecteds[[bs]]
  main_dir <- benchmarks[bs]
  local_methods <- methnames[selected]
  dimnames(res[[bs]]) <- list(local_methods, bif_names[bif_num], timenames)
  for (b in 1:10) {
    benchmark_dir <- paste(main_dir, b, "/", sep = "")
    
    tmp <- readRDS(file = paste(benchmark_dir, "longest_time.csv", sep = ""))
    res[[bs]][, b, ] <- apply(tmp, c(1, 2), mean)
    
    q <- qnorm(0.975)
    err[[bs]][, b, ] <- q * sqrt(apply(tmp, c(1, 2), var)) / 10
  }
  tups[[bs]] <- res[[bs]] + err[[bs]]
  tdowns[[bs]] <- res[[bs]] - err[[bs]]
}

paper_cols <- c("royalblue4", cols)
paper_pch <- c(18, pchs)
paper_lty <- c(1, ltys)

par(mfcol=c(4,2))
for (bs in seq_along(benchmarks)) {
  selected <- selecteds[[bs]]
  main_dir <- benchmarks[bs]
  local_methnames <- methnames[selected]
  for (j in 1:length(timenames)) {
    
    if(bs == 1) {
      par(mar=c(1,3,1,1))
    }
    
    if (j == length(timenames)) {
      par(mar=c(3,1,1,1))
      if (bs == 1) {
        par(mar=c(3,3,1,1))
      }
    }
    
    plot(bif_num, res[[bs]][local_methnames[1],,timenames[j]], ylim=c(0, 1), axes=FALSE, type="n")
    
    if(bs == 1) {
      mtext(leg_names[j], side=2, line=2)
      axis(2, at=c(-1, -0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75, 1.),
           labels = c(-1, -0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75, 1))
    }
    
    if(j == 1) {
      title(main=benchmark_names[bs])
      axis(1, at=bif_num, labels=FALSE)
      axis(2, at=c(-1, -0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75, 1.), labels = FALSE)
    } else if (j == length(timenames)) {
      axis(1, at=bif_num, labels = bif_num+2)
      mtext("#cell fates", side=1, line=2)
    } else {
      axis(1, at=bif_num, labels=FALSE)
      axis(2, at=c(-1, -0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75, 1.), labels = FALSE)
    }
    
    box(which = "plot", lwd=2)
    
    for (i in 1:length(local_methnames)) {
      m = methnames[i]
      mc = cols[i]
      points(bif_num, res[[bs]][local_methnames[i],,timenames[j]],
             ylim=c(0, 0.75),
             pch=paper_pch[i],
             col=paper_cols[i],
             type="b",
             lwd=3,
             lty=paper_lty[i],
             cex=1.5)
      arrows(x0=bif_num,
             y0 = tdowns[[bs]][local_methnames[i],,timenames[j]],
             y1=tups[[bs]][local_methnames[i],,timenames[j]],
             code=3,
             angle=90,
             length=0.1,
             col=paper_cols[i])
    }
  }
}
