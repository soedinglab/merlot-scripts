library(ggplot2)
library(reshape2)
library(scales)
library(plotly)

source("~/Documents/repos/merlot-scripts/plot/produce_legacy_merlot_dest.R")

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

hhtree <- "~/Documents/repos/merlot-scripts"
various <- paste(hhtree, "/scripts/various.R", sep="")
evaluat <- paste(hhtree, "/scripts/evaluate_method.R", sep="")
suppressPackageStartupMessages(source(various))
suppressPackageStartupMessages(source(evaluat))

main_dir <- "~/Documents/data/only_knn/benchmark"
bif_num = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)

methnames <- c("destiny_log_k",
               "MERLoT_log_free_knn_merlot_el_sens",
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

methtitles <- c("destiny_log_k",
                "free_merlot_el_sens",
                "free_elpi_el_sens",
                "k_free_merlot_el_sens",
                "k_free_elpi_el_sens",
                "free_merlot_emb_sens",
                "free_elpi_emb_sens",
                "k_free_merlot_emb_sens",
                "k_free_elpi_emb_sens",
                "fixed_merlot_el",
                "fixed_elpi_el",
                "k_fixed_merlot_el",
                "k_fixed_elpi_el",
                "fixed_merlot_emb",
                "fixed_elpi_emb",
                "k_fixed_merlot_emb",
                "k_fixed_elpi_emb")

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
  benchmark_dir = paste(main_dir, bif_num[i], "/", sep="")
  
  # remove missing files
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
  
  bplot_list <- apply(br_res, 2, function(x) x[x>0])
  boxplot(bplot_list, main=bif_num[i], xaxt = "n",  xlab = "")
  abline(h = leg_res[b, "LPGraph_log_k_free_emb_sens"], col="red")
  axis(1, labels = FALSE, at = seq_along(methnames))
  text(x = seq_along(methnames), y = par("usr")[3] + 0.5, srt = 90, adj = 1,
       labels = methnames, xpd = TRUE)
  
  
  
  q = qnorm(0.975)
  k <- which(res$X == bif_num[i])
  for (m in methnames) {
    res[k, m] = mean(br_res[,m][br_res[,m]>0])
    err[k, m] = q * sqrt(var(br_res[,m][br_res[,m]>0])) / 10
    
    tres[k, m] = mean(ti_res[, m][br_res[,m]>0])
    terr[k, m] = q * sqrt(var(ti_res[,m][br_res[,m]>0])) / 10
  }
  
  # q = qnorm(0.975)
  # k <- which(res$X == bif_num[i])
  # for (m in methnames) {
  #   res[k, m] = mean(br_res[,m])
  #   err[k, m] = q * sqrt(var(br_res[,m])) / 10
  #   
  #   tres[k, m] = mean(ti_res[, m])
  #   terr[k, m] = q * sqrt(var(ti_res[,m])) / 10
  # }
}



ups = res + err
downs = res - err

tups = tres + terr
tdowns = tres - terr

cols = gg_color_hue(length(methnames))
pchs <- rep(20, length(methnames))
ltys <- rep(1, length(methnames))

# plot fixed methods
plot(1, type="n", xlim=c(min(bif_num), max(bif_num)), ylim=c(0.25, 0.9),
     ylab=leg_names[f], xlab="#cell fates", axes=FALSE, main="Branch assignment quality")
box(which = "plot", lwd=2)
axis(1, at=bif_num, labels = bif_num+2)
axis(2, at=c(-1, -0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75, 1.), labels = c(-1, -0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75, 1))

m <- "LPGraph_log_k_fixed_emb"
points(bif_num, leg_res[, m], pch=20, col="black", type="b", lwd=3, cex=1.5)
arrows(x0 = bif_num, y0 = leg_downs[, m], y1=leg_ups[, m], code=3, angle=90, length=0.1, col="black")
for (i in seq_along(methnames)) {
  m = methnames[i]
  if (!grepl("fixed", m)) {
    next
  }
  mc = cols[i]
  points(bif_num, res[, m], pch=pchs[i], col=mc, type="b", lwd=3, cex=1.5, lty=ltys[i])
  arrows(x0 = bif_num, y0 = downs[, m], y1=ups[, m],
         code=3, angle=90, length=0.1, col=mc)
}
fixed <- grepl("fixed", methnames)
legend("bottomright", legend=c("legacy", methtitles[fixed]), pch=21,
       pt.bg = c("black", cols[fixed]), pt.cex=1.5, ncol=2)


# plot free methods
plot(1, type="n", xlim=c(min(bif_num), max(bif_num)), ylim=c(0.25, 0.9),
     ylab=leg_names[f], xlab="#cell fates", axes=FALSE, main="Branch assignment quality")
box(which = "plot", lwd=2)
axis(1, at=bif_num, labels = bif_num+2)
axis(2, at=c(-1, -0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75, 1.), labels = c(-1, -0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75, 1))

m <- "LPGraph_log_k_free_emb_sens"
points(bif_num, leg_res[, m], pch=20, col="black", type="b", lwd=3, cex=1.5)
arrows(x0 = bif_num, y0 = leg_downs[, m], y1=leg_ups[, m], code=3, angle=90, length=0.1, col="black")
for (i in seq_along(methnames)) {
  m = methnames[i]
  if (!grepl("free", m)) {
    next
  }
  mc = cols[i]
  points(bif_num, res[, m], pch=pchs[i], col=mc, type="b", lwd=3, cex=1.5, lty=ltys[i])
  arrows(x0 = bif_num, y0 = downs[, m], y1=ups[, m],
         code=3, angle=90, length=0.1, col=mc)
}
free <- grepl("free", methnames)
legend("topright", legend=c("legacy", methtitles[free]), pch=21,
       pt.bg = c("black", cols[free]), pt.cex=1.5, ncol = 3)




# plot fixed pseudotime
plot(1, type="n", xlim=c(min(bif_num), max(bif_num)), ylim=c(0.4, 1.0),
     ylab=leg_names[f], xlab="#cell fates", axes=FALSE, main="Pseudotime")
box(which = "plot", lwd=2)
axis(1, at=bif_num, labels = bif_num+2)
axis(2, at=c(-1, -0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75, 1.), labels = c(-1, -0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75, 1))

m <- "LPGraph_log_k_fixed_emb"
points(bif_num, leg_tres[, m], pch=20, col="black", type="b", lwd=3, cex=1.5)
arrows(x0 = bif_num, y0 = leg_tdowns[, m], y1=leg_tups[, m], code=3, angle=90, length=0.1, col="black")
for (i in seq_along(methnames)) {
  m = methnames[i]
  if (!grepl("fixed", m)) {
    next
  }
  mc = cols[i]
  points(bif_num, tres[, m], pch=pchs[i], col=mc, type="b", lwd=3, cex=1.5, lty=ltys[i])
  arrows(x0 = bif_num, y0 = tdowns[, m], y1=tups[, m],
         code=3, angle=90, length=0.1, col=mc)
}
fixed <- grepl("fixed", methnames)
legend("bottom", legend=c("legacy", methtitles[fixed]), pch=21,
       pt.bg = c("black", cols[fixed]), pt.cex=1.5)



# plot free pseudotime
plot(1, type="n", xlim=c(min(bif_num), max(bif_num)), ylim=c(0.4, 1.0),
     ylab=leg_names[f], xlab="#cell fates", axes=FALSE, main="Pseudotime")
box(which = "plot", lwd=2)
axis(1, at=bif_num, labels = bif_num+2)
axis(2, at=c(-1, -0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75, 1.), labels = c(-1, -0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75, 1))

m <- "LPGraph_log_k_free_emb_sens"
points(bif_num, leg_tres[, m], pch=20, col="black", type="b", lwd=3, cex=1.5)
arrows(x0 = bif_num, y0 = leg_tdowns[, m], y1=leg_tups[, m], code=3, angle=90, length=0.1, col="black")
for (i in seq_along(methnames)) {
  m = methnames[i]
  if (!grepl("free", m)) {
    next
  }
  mc = cols[i]
  points(bif_num, tres[, m], pch=pchs[i], col=mc, type="b", lwd=3, cex=1.5, lty=ltys[i])
  arrows(x0 = bif_num, y0 = tdowns[, m], y1=tups[, m],
         code=3, angle=90, length=0.1, col=mc)
}
fixed <- grepl("free", methnames)
legend("bottom", legend=c("legacy", methtitles[fixed]), pch=21,
       pt.bg = c("black", cols[fixed]), pt.cex=1.5)











p <- plot_ly(type = "scatter", mode="markers+lines") %>%
  layout(title = "Average branch assignment performance",
         xaxis = list(title = "#bifurcations"),
         yaxis = list (title = "adjusted MI", range=c(0, 1)))
for (i in seq_along(methnames)) {
  p <- add_trace(p, x=1:10, y=res[,methnames[i]], name=methtitles[i])
}
p
