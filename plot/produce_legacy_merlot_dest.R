#!/usr/bin/Rscript
source("~/Documents/repos/merlot-scripts/plot/auxiliary.R")

tres = data.frame(matrix(data=0, ncol=length(methnames)+1, nrow=length(bif_num)))
colnames(tres) = c("X", methnames)
tres$X = bif_num

terr = data.frame(matrix(data=0, ncol=length(methnames)+1, nrow=length(bif_num)))
colnames(terr) = c("X", methnames)
terr$X = bif_num

timenames <- c("goodman-kruskall (unweighted)", "goodman-kruskall (weighted)",
               "kendall index (unweighted)", "kendall index (weighted)", "LPGK")

bench_dir <- "~/Documents/data/prosstt/benchmark"
selected <- prosstt_selected[c(3, 5)]
  
f <- 6
chosen <- funcnames[f]
tchosen = timenames[5]

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
  ti_all <- parse_benchmark(subfolders, benchmark_dir, timenames, methnames, "time", replNA = TRUE)
  ti_all <- as.data.frame(ti_all)
  ti_res <- ti_all[,1:length(methnames)]
  ti_res <- apply(ti_res, 2, as.numeric)
  ti_res = ti_res[ti_all$measure == tchosen, ]
  
  # calculate confidence intervals for error bars
  q = qnorm(0.975)
  k <- which(res$X == bif_num[b])
  for (m in methnames) {
    res[k, m] = mean(br_res[,m])
    err[k, m] = q * sqrt(var(br_res[,m])) / 10
    
    tres[k, m] = mean(ti_res[, m])
    terr[k, m] = q * sqrt(var(ti_res[,m])) / 10
  }
}

# (symmetric) error bars
ups = res + err
downs = res - err

tups = tres + terr
tdowns = tres - terr

leg_res <- res[methnames[selected]]
leg_ups <- ups[methnames[selected]]
leg_downs <- downs[methnames[selected]]

leg_tres <- tres[methnames[selected]]
leg_tups <- tups[methnames[selected]]
leg_tdowns <- tdowns[methnames[selected]]
