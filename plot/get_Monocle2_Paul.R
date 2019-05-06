rm(list = ls()) # clear the environment 
# load all the necessary libraries 
options(warn=-1) # turn off warning message globally 
suppressMessages(library(monocle))
suppressMessages(library(plyr))
suppressMessages(library(dplyr))

setwd("~/Documents/data/paul/")
# this RData is from Maren Büttner (https://github.com/theislab/scAnalysisTutorial)
load('./Paul_Cell_MARSseq_GSE72857.RData')
# the following code is used to select feature genes used by Maren 
gene.names <-sapply(strsplit(rownames(data.debatched), ";"), "[", 1)
is.informative <- gene.names %in% info.genes[order(info.genes)]
data.info.genes <- data.debatched[is.informative,]
rownames(data.info.genes) <- gene.names[is.informative]

CellTypes= paste(find.package("merlot"), "/example/PaulCellsMapping.csv", sep="")
MAP_cells_clusters <- read.table(file =CellTypes, header=F, stringsAsFactors = F, sep=",")
row.names(MAP_cells_clusters) <- MAP_cells_clusters$V1
# this .txt is from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE72857
design_mat <- read.table('./GSE72857_experimental_design.txt', header = T, row.names = 1, skip = 19, sep = '\t')
design_mat$cluster <- MAP_cells_clusters[row.names(design_mat), 'V2']
valid_design_mat <- subset(design_mat, !is.na(cluster))

#filtering cells to include only the ones which were assigned a cluster id: 
# this .txt is from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE72857
valid_subset_GSE72857_exprs <- read.table('./GSE72857_umitab.txt', header = T, row.names = 1)

# Get the intersect gene used by Maren Büttner and the genes we have 
common_genes <- rownames(valid_subset_GSE72857_exprs)[rownames(valid_subset_GSE72857_exprs) %in% info.genes]
fd <- new("AnnotatedDataFrame", data = data.frame(gene_short_name = common_genes, row.names = common_genes))
pd <- new("AnnotatedDataFrame", data = valid_design_mat)

# create a CDS with data.info.genes 
valid_subset_GSE72857_cds <- newCellDataSet(as(as.matrix(data.info.genes[common_genes, ]), 'sparseMatrix'), 
                                            phenoData = pd, 
                                            featureData = fd,
                                            lowerDetectionLimit=1,
                                            expressionFamily=negbinomial.size())
valid_subset_GSE72857_cds <- estimateSizeFactors(valid_subset_GSE72857_cds)
valid_subset_GSE72857_cds <- estimateDispersions(valid_subset_GSE72857_cds)

pData(valid_subset_GSE72857_cds)$cell_type <- revalue(as.character(pData(valid_subset_GSE72857_cds)$cluster), 
                                                      c("1" = 'erythroid',
                                                        "2" = 'erythroid',
                                                        "3" = 'erythroid',
                                                        "4" = 'erythroid',
                                                        "5" = 'erythroid',
                                                        "6" = 'erythroid', 
                                                        "7" = 'CMP',
                                                        "8" = 'CMP',
                                                        "9" = 'CMP',
                                                        "10" = 'CMP',
                                                        "11" = 'DC', 
                                                        "12" = 'GMP',
                                                        "13" = 'GMP',
                                                        "14" = 'GMP',
                                                        "15" = 'GMP',
                                                        "16" = 'GMP',
                                                        "17" = 'GMP',
                                                        "18" = 'GMP', 
                                                        "19" = 'lymphoid'))

#remove all lymphoid cells as it is not belong to myeloid lineage 
valid_subset_GSE72857_cds <- valid_subset_GSE72857_cds[, pData(valid_subset_GSE72857_cds)$cell_type != 'lymphoid']

options(repr.plot.width=4, repr.plot.height=3)
plot_pc_variance_explained(valid_subset_GSE72857_cds) + geom_vline(xintercept = 10)
valid_subset_GSE72857_cds2 <- reduceDimension(valid_subset_GSE72857_cds, norm_method = 'log', verbose = F, max_components = 10) 
valid_subset_GSE72857_cds2 <- orderCells(valid_subset_GSE72857_cds2, reverse = T)

# make multiple branch plot
detailed_cell_type_color <- c("B" = "forestgreen",
                              "DC" = "firebrick3",
                              "E" = "dodgerblue3",
                              "Ery" = "darkorchid",
                              "M" = "orange",
                              "MP/EP" = "aquamarine",
                              "GMP" = "darkorange3",
                              "MK" = "blue",
                              "N" = "magenta")

pData(valid_subset_GSE72857_cds2)$cell_type2 <- revalue(as.character(pData(valid_subset_GSE72857_cds2)$cluster), 
                                                        c("1" = 'Ery', "2" = 'Ery', "3" = 'Ery', "4" = 'Ery', "5" = 'Ery', "6" = 'Ery', 
                                                          "7" = 'MP/EP', "8" = 'MK', "9" = 'GMP', "10" = 'GMP',
                                                          "11" = 'DC', 
                                                          "12" = 'B', "13" = 'B', "14" = 'M', "15" = 'M', "16" = 'N', "17" = 'N', "18" = 'E', 
                                                          "19" = 'lymphoid'))

options(repr.plot.width=4, repr.plot.height=3)
# svg("~/Documents/presentations/2019-05_Paul/monocle_plot.svg", width = 7, height = 7)
plot_complex_cell_trajectory(valid_subset_GSE72857_cds2, color_by = "cell_type2", show_branch_points = T, 
                             cell_size = 1, cell_link_size = 0.3, root_states = c(10)) + scale_size(range = c(0.2, 0.2)) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
  theme (legend.position="left", legend.title=element_blank()) + scale_color_manual(values = detailed_cell_type_color)
# dev.off()

# MERLoT scaffold+elastic tree
broad_type <- (valid_subset_GSE72857_cds2@phenoData$cell_type2)
CellCoordinates <- t(valid_subset_GSE72857_cds2@reducedDimS)

