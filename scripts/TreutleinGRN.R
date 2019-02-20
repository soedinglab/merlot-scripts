source('/home/gonzalo/Desktop/Johannes/MERLoT/GRNCode/gene_network_analysis.R')
library(merlot)
library(data.table)

# Read the example Guo dataset that is distributed with the package
DataFile= paste(find.package("merlot"), "/example/Treutlein2014.txt", sep="")
Dataset=ReadDataset(DataFile)

# Embed Cells into their manifold, in this case we use Diffusion Maps as calculated by Destiny
library(destiny)
DatasetDM <- DiffusionMap(Dataset$ExpressionMatrix)
# we take the first two diffusion coordinates
CellCoordinates=DatasetDM@eigenvectors[,1:2]
# We calculate the scaffold tree for MERLoT
ScaffoldTree=CalculateScaffoldTree(CellCoordinates = CellCoordinates)

# This is the number of nodes to be used for the Principal Elastic Tree
NumberOfNodes=200

# We calculate the elastic principal tree using the scaffold tree for its initialization.
ElasticTree= CalculateElasticTree(ScaffoldTree = ScaffoldTree, N_yk = NumberOfNodes)
# We embed the elastic tree into the expression space.
#We modified the elasticity parameters in order to produce more stretched pseudotemporal gene expression profiles
EmbeddedTree=GenesSpaceEmbedding(Dataset$ExpressionMatrix, ElasticTree, lambda_0 = 2.03e-09,
                                    mu_0 = 0.00625, increaseFactor_mu = 100, increaseFactor_lambda = 100)

Pseudotimes=CalculatePseudotimes(EmbeddedTree, T0=2)
plot_pseudotimes(CellCoordinates, Pseudotimes)

# We calculate differentially expressed genes for the different branches in the tree
FbrDE<-branch_differential_expression(3, EmbeddedTree = EmbeddedTree)
NbrDE<-branch_differential_expression(2, EmbeddedTree = EmbeddedTree)
MbrDE<-branch_differential_expression(1, EmbeddedTree = EmbeddedTree)

# Make summary of differentially expressed genes
FDEgenes<-gatherDEgenes(FbrDE)
NDEgenes<-gatherDEgenes(NbrDE)
MDEgenes<-gatherDEgenes(MbrDE)

# Plot the reconstructed lineage tree
# svg("/home/gonzalo/Desktop/Johannes/MERLoT/GRNs/Tree.svg", width = 6, height = 5)
plot_elastic_tree(ElasticTree)
# dev.off()

# We calculate the GRN. Modes can be "tree" for using tree imputed gene expression values or "cells" for raw cell gene expression values
# using imputed values and an edge pearson's correlation threshold of 0.9
treutgrns<-list()
treutgrns$imputed0.9 = makeGRN(emb_tree = EmbeddedTree, globalmode = T, cutoff = .9, mode = 'tree' )
treutgrns$nonimputed0.9 = makeGRN(emb_tree = EmbeddedTree, globalmode = T, cutoff = .9, mode = 'cells')

# Density plots whole tree
density_imputed=density(treutgrns$imputed0.9$full_cors)
density_nonimputed=density(treutgrns$nonimputed0.9$full_cors)

# svg("/home/gonzalo/Desktop/Johannes/MERLoT/GRNs/correlations.svg", width = 6, height = 5)
plot(density_imputed$x, density_imputed$y, xlim=c(-1.1,1.1), ylim=c(0,1.4), type="l", col="darkgreen", lwd=2, xlab="Pearson's correlation coeff", ylab="Density")
lines(density_nonimputed$x, density_nonimputed$y, col="darkred", lwd=2)
abline(v=0.9)
box(lwd=2)
# dev.off()

# svg("/home/gonzalo/Desktop/Johannes/MERLoT/GRNs/gene_neu.svg", width = 6, height = 5)
plot_pseudotime_expression_gene(GeneName = "Ank2" , EmbeddedTree = EmbeddedTree, Pseudotimes = Pseudotimes, addlegend = F, range_y = "cells")
# dev.off()

# svg("/home/gonzalo/Desktop/Johannes/MERLoT/GRNs/gene_myo.svg", width = 6, height = 5)
plot_pseudotime_expression_gene(GeneName = MDEgenes$GeneName[1] , EmbeddedTree = EmbeddedTree, Pseudotimes = Pseudotimes, addlegend = F, range_y = "cells")
# dev.off()

# svg("/home/gonzalo/Desktop/GRNNames.svg", width = 12, height = 10)
plotGRN(treutgrns$imputed0.9, DEdifference = T, DEgenes= NDEgenes, sizeNodes = 3, title='Imputed neuron', plotnamesofgenes = '')
# dev.off()

########################
# Neuronal figures
########################
###############################################
# We create different GRNs using different
# thresholds for the non-imputed values.
###############################################
# svg("/home/gonzalo/Desktop/Johannes/MERLoT/GRNs/non_imputed0.9.svg", width = 5, height = 5)
plotGRN(treutgrns$nonimputed0.9, DEdifference = T, DEgenes= NDEgenes, sizeNodes = 3, title='', plotnamesofgenes = '')
# dev.off()

# We try different lower thresholds for reconstructing the GRN with non imputed values
# svg("/home/gonzalo/Desktop/Johannes/MERLoT/GRNs/non_imputed0.8.svg", width = 5, height = 5)
treutgrns$nonimputed0.8 = makeGRN(emb_tree = EmbeddedTree, globalmode = T, cutoff = .8, mode = 'cells')
plotGRN(treutgrns$nonimputed0.8, DEdifference = T, DEgenes = NDEgenes, sizeNodes = 3, title='', plotnamesofgenes = '')
# dev.off()

# svg("/home/gonzalo/Desktop/Johannes/MERLoT/GRNs/non_imputed0.7.svg", width = 5, height = 5)
treutgrns$nonimputed0.7 = makeGRN(emb_tree = EmbeddedTree, globalmode = T, cutoff = .7, mode = 'cells')
plotGRN(treutgrns$nonimputed0.7, DEdifference = T, DEgenes = NDEgenes, sizeNodes = 3, title='', plotnamesofgenes = '')
# dev.off()

# svg("/home/gonzalo/Desktop/Johannes/MERLoT/GRNs/non_imputed0.6.svg", width = 5, height = 5)
treutgrns$nonimputed0.6 = makeGRN(emb_tree = EmbeddedTree, globalmode = T, cutoff = .6, mode = 'cells')
plotGRN(treutgrns$nonimputed0.6, DEdifference = T, DEgenes = NDEgenes, sizeNodes = 3, title='', plotnamesofgenes = '')
# dev.off()

# svg("/home/gonzalo/Desktop/Johannes/MERLoT/GRNs/non_imputed0.5.svg", width = 5, height = 5)
treutgrns$nonimputed0.5 = makeGRN(emb_tree = EmbeddedTree, globalmode = T, cutoff = .5, mode = 'cells')
plotGRN(treutgrns$nonimputed0.5, DEdifference = T, DEgenes = NDEgenes, sizeNodes = 3, title='', plotnamesofgenes = '')
# dev.off()

# svg("/home/gonzalo/Desktop/Johannes/MERLoT/GRNs/non_imputed0.4.svg", width = 5, height = 5)
treutgrns$nonimputed0.4 = makeGRN(emb_tree = EmbeddedTree, globalmode = T, cutoff = .4, mode = 'cells')
plotGRN(treutgrns$nonimputed0.4, DEdifference = T, DEgenes = NDEgenes, sizeNodes = 3, title='', plotnamesofgenes = '')
# dev.off()

# Clustering
# svg("/home/gonzalo/Desktop/Johannes/MERLoT/GRNs/Mnon_imputed0.4.svg", width = 5, height = 5)
treutclust<-lapply(treutgrns, clusterGRN, method='label_propagation', min.clust=3)
plotGRN(treutgrns$imputed0.9, clustering = treutclust$imputed0.9, sizeNodes = 2)
# dev.off()

## Annotations
annot<-makeGOAnnotation(mapping='org.Mm.eg.db',ID='alias')
GObrWhBP <- GRN_GO_termenrichment(annotation=annot,clustering= treutclust$nonimputed0.4,ontology='BP', n = 10)
# GObrWhMF <- GRN_GO_termenrichment(annotation=annot,clustering= treutclust$wholetree,ontology='MF', n = 10)
# GObrWhCC <- GRN_GO_termenrichment(annotation=annot,clustering= treutclust$wholetree,ontology='CC', n = 10)
#
Whole_MF_GRN_GO_termenrichment <-GRN_GO_termenrichment(annotation=annot,clustering= treutclust$nonimputed0.4,ontology='BP', n = 10)
# F_MF_GRN_GO_termenrichment <-GRN_GO_termenrichment(annotation=annot,clustering= treutclust$brF,ontology='BP', n = 10)
# M_MF_GRN_GO_termenrichment <-GRN_GO_termenrichment(annotation=annot,clustering= treutclust$brM,ontology='BP', n = 10)
# N_MF_GRN_GO_termenrichment <-GRN_GO_termenrichment(annotation=annot,clustering= treutclust$brN,ontology='BP', n = 10)
#
write.csv(x = rbindlist(GObrWhBP, idcol="ID"), file = '/home/gonzalo/Desktop/Johannes/MERLoT/GRNs/GObrWhBP_cells_50.csv')
# write.csv(x = rbindlist(GObrWhMF, idcol="ID"), file = 'GObrWhMF_cells_50.csv')
# write.csv(x = rbindlist(GObrWhCC, idcol="ID"), file = 'GObrWhCC_cells_50.csv')


## Annotations for imputed values GRN

GO_imp_Wh_MF <- GRN_GO_termenrichment(annotation=annot,clustering= treutclust$imputed0.9,ontology='MF', n = 10)
GO_imp_Wh_CC <- GRN_GO_termenrichment(annotation=annot,clustering= treutclust$imputed0.9,ontology='CC', n = 10)
GO_imp_Wh_BP <- GRN_GO_termenrichment(annotation=annot,clustering= treutclust$imputed0.9,ontology='BP', n = 10)

write.csv(x = rbindlist(GO_imp_Wh_MF, idcol="ID"), file = '/home/gonzalo/Desktop/Johannes/MERLoT/GRNs/GO_imp_Whole_MF.csv')
write.csv(x = rbindlist(GO_imp_Wh_CC, idcol="ID"), file = '/home/gonzalo/Desktop/Johannes/MERLoT/GRNs/GO_imp_Whole_CC.csv')
write.csv(x = rbindlist(GO_imp_Wh_BP, idcol="ID"), file = '/home/gonzalo/Desktop/Johannes/MERLoT/GRNs/GO_imp_Whole_BP.csv')

