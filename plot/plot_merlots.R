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

methnames <- c("MERLoT_log_free_knn_merlot_el_sens",
               "MERLoT_log_free_knn_elpi_el_sens",
               "MERLoT_log_k_free_knn_merlot_el_sens",
               "MERLoT_log_k_free_knn_elpi_el_sens",
               "MERLoT_log_free_knn_merlot_emb_sens",
               "MERLoT_log_free_knn_elpi_emb_sens",
               "MERLoT_log_k_free_knn_merlot_emb_sens",
               "MERLoT_log_k_free_knn_elpi_emb_sens",
               "MERLoT_log_fixed_knn_merlot_el",
               "MERLoT_log_fixed_knn_elpi_el",
               "MERLoT_log_k_fixed_knn_merlot_el",
               "MERLoT_log_k_fixed_knn_elpi_el",
               "MERLoT_log_fixed_knn_merlot_emb",
               "MERLoT_log_fixed_knn_elpi_emb",
               "MERLoT_log_k_fixed_knn_merlot_emb",
               "MERLoT_log_k_fixed_knn_elpi_emb")


res = data.frame(matrix(data=0, ncol=length(methnames)+1, nrow=length(bif_num)))
colnames(res) = c("X", methnames)
res$X = bif_num

tres = data.frame(matrix(data=0, ncol=length(methnames)+1, nrow=length(bif_num)))
colnames(tres) = c("X", methnames)
tres$X = bif_num

err = data.frame(matrix(data=0, ncol=length(methnames)+1, nrow=length(bif_num)))
colnames(err) = c("X", methnames)
err$X = bif_num

terr = data.frame(matrix(data=0, ncol=length(methnames)+1, nrow=length(bif_num)))
colnames(terr) = c("X", methnames)
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
}

ups = res + err
downs = res - err

tups = tres + terr
tdowns = tres - terr

cols = gg_color_hue(4)
selected <- c(3,4,7,8,11,12,14,15)
ltys <- c(1,2,1,2,1,2,2,1)

svg("/home/npapado/Documents/presentations/2019-03_benchmark_knn/merlots_branch.svg", height=5, width=6.5)
# plot fixed methods
plot(1, type="n", xlim=c(min(bif_num), max(bif_num)), ylim=c(0.5, 0.85),
     ylab=leg_names[f], xlab="#cell fates", axes=FALSE, main="Branch assignment quality")
box(which = "plot", lwd=2)
axis(1, at=bif_num, labels = bif_num+2)
axis(2, at=c(0.5, 0.6, 0.7, 0.8, 0.9, 1.), labels = c(0.5, 0.6, 0.7, 0.8, 0.9, 1.))

for (i in seq_along(methnames[selected])) {
  m = methnames[selected][i]
  mc = cols[as.integer(i/2) + i%%2]
  points(bif_num, res[, m], pch=20, col=mc, type="b", lwd=3, cex=1.5, lty=ltys[i])
  arrows(x0 = bif_num, y0 = downs[, m], y1=ups[, m],
         code=3, angle=90, length=0.1, col=mc)
}
dev.off()

# plot fixed pseudotime
svg("/home/npapado/Documents/presentations/2019-03_benchmark_knn/merlots_pt.svg", height=5, width=6.5)
plot(1, type="n", xlim=c(min(bif_num), max(bif_num)), ylim=c(0.5, 0.95),
     ylab=leg_names[f], xlab="#cell fates", axes=FALSE, main="Pseudotime assignment quality")
box(which = "plot", lwd=2)
axis(1, at=bif_num, labels = bif_num+2)
axis(2, at=c(0.5, 0.6, 0.7, 0.8, 0.9, 1.), labels = c(0.5, 0.6, 0.7, 0.8, 0.9, 1.))

for (i in seq_along(methnames[selected])) {
  m = methnames[selected][i]
  mc = cols[as.integer(i/2) + i%%2]
  points(bif_num, tres[, m], pch=20, col=mc, type="b", lwd=3, cex=1.5, lty=ltys[i])
  arrows(x0 = bif_num, y0 = tdowns[, m], y1=tups[, m],
         code=3, angle=90, length=0.1, col=mc)
}
dev.off()

svg("/home/npapado/Documents/presentations/2019-03_benchmark_knn/legend_merlots.svg", height=6, width=8)
plot(1, type="n", xlim=c(min(bif_num), max(bif_num)), ylim=c(0., 1.0), axes=FALSE)
legend_names <- c("auto, elastic",
                  "auto, embedded",
                  "fixed, elastic",
                  "fixed, embedded",
                  "MERLoT hyperparameters",
                  "ElPiGraph.R hyperparameters")
legend("bottomleft", legend=legend_names, ncol=1, lwd=c(1,1,1,1,3,3),
       pch=c(21, 21, 21, 21, NA, NA), pt.bg = c(cols, NA, NA), lty=c(NA, NA, NA, NA, 1, 2), pt.cex=1.5)
dev.off()

# # error_bars <- ups - downs
# # error_bars$X <- 1:10
# # melted_error <- melt(error_bars, id.vars = 1)
melted <- melt(res, id.vars = 1)
# melted[, "error"] <- melted_error$value

p <- plot_ly(melted, x = ~X, y = ~value, #error_y = ~list(array = error),
             type = "scatter",
             mode="markers+lines",
             colors = cols,
             color = ~variable) %>%
  layout(title = "Average branch assignment performance",
         xaxis = list(title = "#bifurcations"),
         yaxis = list (title = "adjusted MI", range=c(0, 1)))
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
