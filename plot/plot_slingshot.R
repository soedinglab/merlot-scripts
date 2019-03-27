library(ggplot2)
library(reshape2)
library(scales)
library(plotly)

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

hhtree <- "~/Documents/repos/merlot-scripts"
various <- paste(hhtree, "/scripts/various.R", sep="")
evaluat <- paste(hhtree, "/scripts/evaluate_method.R", sep="")
suppressPackageStartupMessages(source(various))
suppressPackageStartupMessages(source(evaluat))

nosing <- "~/Documents/data/no_singularity/benchmark"
sing <- "~/Documents/data/singularity_corrected/benchmark"
bif_num = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)

methnames <- c("destiny_log_k",
               "MERLoT_log_k_free_knn_elpi_emb_sens",
               "MERLoT_log_k_fixed_knn_elpi_emb",
               "slingshot_destiny_log_k",
               "slingshot_destiny_log",
               "monocl2")

res = data.frame(matrix(data=0, ncol=length(methnames)+3, nrow=length(bif_num)))
colnames(res) = c("X", methnames, "sing_slingshot_destiny_log", "sing_slingshot_destiny_log_k")
res$X = bif_num

tres = data.frame(matrix(data=0, ncol=length(methnames)+3, nrow=length(bif_num)))
colnames(tres) = c("X", methnames, "sing_slingshot_destiny_log", "sing_slingshot_destiny_log_k")
tres$X = bif_num

err = data.frame(matrix(data=0, ncol=length(methnames)+3, nrow=length(bif_num)))
colnames(err) = c("X", methnames, "sing_slingshot_destiny_log", "sing_slingshot_destiny_log_k")
err$X = bif_num

terr = data.frame(matrix(data=0, ncol=length(methnames)+3, nrow=length(bif_num)))
colnames(terr) = c("X", methnames, "sing_slingshot_destiny_log", "sing_slingshot_destiny_log_k")
terr$X = bif_num

funcnames <- c("rand index", "matthews corr. coef.",
               "F1 measure", "jaccard index", "fowkles-mallows index", "adjusted MI", "#branches")
leg_names <- c("Rand Index", "Matthews Corr. Coef.", "F1 Measure", "Jaccard index", "Fowkles-Mallows index", "Normalized MI")
shortnames <- c("rand", "mcc", "f1", "jaccard", "fowklesm", "adj_mi")
# funcnames <- c("matthews corr. coef.", "adjusted MI")
timenames <- c("goodman-kruskall (unweighted)", "goodman-kruskall (weighted)",
               "kendall index (unweighted)", "kendall index (weighted)", "LPGK")
ylims <- list(c(0,1), c(-1., 0.6), c(0, 0.8), c(0, 0.6), c(0, 0.8), c(0, 0.75))


tchosen = timenames[5]
f <- 6
chosen <- funcnames[f]
for (i in 1:length(bif_num)) {
  benchmark_dir = paste(nosing, bif_num[i], "/", sep="")
  subfolders <- dir(benchmark_dir)
  
  br_all <- parse_cluster_output(subfolders, benchmark_dir, funcnames, methnames, "branch", replNA = TRUE)
  br_all <- as.data.frame(br_all)
  br_res <- br_all[,1:length(methnames)]
  br_res <- apply(br_res, 2, as.numeric)
  br_res = br_res[br_all$measure == chosen,]
  rownames(br_res) <- br_all$experiment[br_all$measure == chosen]
  
  ti_all <- parse_cluster_output(subfolders, benchmark_dir, timenames, methnames, "time", replNA = TRUE)
  ti_all <- as.data.frame(ti_all)
  ti_res <- ti_all[,1:length(methnames)]
  ti_res <- apply(ti_res, 2, as.numeric)
  ti_res = ti_res[ti_all$measure == tchosen, ]
  
  q = qnorm(0.975)
  k <- which(res$X == bif_num[i])
  for (m in methnames) {
    res[k, m] = mean(br_res[,m])
    err[k, m] = q * sqrt(var(br_res[,m])) / 10
    
    tres[k, m] = mean(ti_res[, m])
    terr[k, m] = q * sqrt(var(ti_res[,m])) / 10
  }
  
  # now from the non-singularity corrected
  benchmark_dir = paste(sing, bif_num[i], "/", sep="")
  subfolders <- dir(benchmark_dir)
  
  br_all <- parse_cluster_output(subfolders, benchmark_dir, funcnames, methnames, "branch", replNA = TRUE)
  br_all <- as.data.frame(br_all)
  br_res <- br_all[,1:length(methnames)]
  br_res <- apply(br_res, 2, as.numeric)
  br_res = br_res[br_all$measure == chosen,]
  rownames(br_res) <- br_all$experiment[br_all$measure == chosen]
  
  ti_all <- parse_cluster_output(subfolders, benchmark_dir, timenames, methnames, "time", replNA = TRUE)
  ti_all <- as.data.frame(ti_all)
  ti_res <- ti_all[,1:length(methnames)]
  ti_res <- apply(ti_res, 2, as.numeric)
  ti_res = ti_res[ti_all$measure == tchosen, ]
  
  res[k, 8] = mean(br_res[, "slingshot_destiny_log"])
  err[k, 8] = q * sqrt(var(br_res[, "slingshot_destiny_log"])) / 10
  tres[k, 8] = mean(ti_res[, "slingshot_destiny_log"])
  terr[k, 8] = q * sqrt(var(ti_res[, "slingshot_destiny_log"])) / 10
  res[k, 9] = mean(br_res[, "slingshot_destiny_log_k"])
  err[k, 9] = q * sqrt(var(br_res[, "slingshot_destiny_log_k"])) / 10
  tres[k, 9] = mean(ti_res[, "slingshot_destiny_log_k"])
  terr[k, 9] = q * sqrt(var(ti_res[, "slingshot_destiny_log_k"])) / 10
}



ups = res + err
downs = res - err

tups = tres + terr
tdowns = tres - terr

cols = gg_color_hue(length(methnames))
cols <- c(cols, cols[c(5, 4)])
ltys <- c(1,1,1,1,1,1,3,3)
legend_ltys <- c(NA, NA, NA, NA, NA, NA, 3)
pchs <- c(21, 21, 21, 21, 21, 21, NA)

svg("/home/npapado/Documents/presentations/2019-03_benchmark_knn/branch.svg", height=6, width=8)
# plot fixed methods
plot(1, type="n", xlim=c(min(bif_num), max(bif_num)), ylim=c(0.25, 0.85),
     ylab=leg_names[f], xlab="#cell fates", axes=FALSE, main="Branch assignment quality")
box(which = "plot", lwd=2)
axis(1, at=bif_num, labels = bif_num+2)
axis(2, at=c(-1, -0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75, 1.), labels = c(-1, -0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75, 1))

for (i in seq_along(colnames(res[2:9]))) {
  m = colnames(res[2:9])[i]
  mc = cols[i]
  points(bif_num, res[, m], pch=20, col=mc, type="b", lwd=3, cex=1.5, lty=ltys[i])
  arrows(x0 = bif_num, y0 = downs[, m], y1=ups[, m],
         code=3, angle=90, length=0.1, col=mc)
}
legend("bottomleft", legend=c(methnames, "singularity corrected"), ncol=1, lwd=c(1,1,1,1,1,1,3),
       pch=pchs, pt.bg = c(cols[1:6], NA), lty=legend_ltys, pt.cex=1.5)
dev.off()

# plot fixed pseudotime
svg("/home/npapado/Documents/presentations/2019-03_benchmark_knn/pseudotime.svg", height=6, width=8)
plot(1, type="n", xlim=c(min(bif_num), max(bif_num)), ylim=c(0.4, 1.0),
     ylab=leg_names[f], xlab="#cell fates", axes=FALSE, main="Pseudotime")
box(which = "plot", lwd=2)
axis(1, at=bif_num, labels = bif_num+2)
axis(2, at=c(-1, -0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75, 1.), labels = c(-1, -0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75, 1))

for (i in seq_along(colnames(res[2:9]))) {
  m = colnames(res[2:9])[i]
  mc = cols[i]
  points(bif_num, tres[, m], pch=20, col=mc, type="b", lwd=3, cex=1.5, lty=ltys[i])
  arrows(x0 = bif_num, y0 = tdowns[, m], y1=tups[, m],
         code=3, angle=90, length=0.1, col=mc)
}
legend("bottomright", legend=c(methnames, "singularity corrected"), ncol=1, lwd=c(1,1,1,1,1,1,3),
       pch=pchs, pt.bg = c(cols[1:6], NA), lty=legend_ltys, pt.cex=1.5)
dev.off()


# # error_bars <- ups - downs
# # error_bars$X <- 1:10
# # melted_error <- melt(error_bars, id.vars = 1)
# melted <- melt(res, id.vars = 1)
# # melted[, "error"] <- melted_error$value
# 
# p <- plot_ly(melted, x = ~X, y = ~value, #error_y = ~list(array = error),
#              type = "scatter",
#              mode="markers+lines",
#              colors = cols,
#              color = ~variable) %>%
#   layout(title = "Average branch assignment performance",
#          xaxis = list(title = "#bifurcations"),
#          yaxis = list (title = "adjusted MI", range=c(0, 1)))
# htmlwidgets::saveWidget(as_widget(p), "/home/npapado/Documents/presentations/2019-03_benchmark_knn/branches_nosingularity.html")
# 
# 
# melted <- melt(tres, id.vars = 1)
# p <- plot_ly(melted, x = ~X, y = ~value,
#              type = "scatter",
#              mode="markers+lines",
#              colors = cols,
#              color = ~variable) %>%
#   layout(title = "Average pseudotime prediction performance",
#          xaxis = list(title = "#bifurcations"),
#          yaxis = list (title = "Longest Path GK index", range=c(0., 1)))
# htmlwidgets::saveWidget(as_widget(p), "/home/npapado/Documents/presentations/2019-03_benchmark_knn/pt_nosingularity.html")
# 
# 
