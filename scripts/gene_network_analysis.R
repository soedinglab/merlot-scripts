# setwd("/home/gonzalo/Desktop/Johannes/MERLoT/GRNs")
library(igraph)
library(viridis)
library(ggplot2)
library(ggnetwork) # to convert between different graph objects - necessary to make igraph work with ggplot
library(network)
library(intergraph)
library(ggrepel)
source("http://bioconductor.org/biocLite.R")
library(AnnotationDbi)
library(org.Mm.eg.db)
library(org.Hs.eg.db)
library(topGO)


# Auxialliary function to define a rescaling function given original range and desired range. Output is a function!
def.rescale<-function(x,y,xprime=0,yprime=1){
    FUN<-function(z){
      rescalingA = (yprime-xprime)/(y-x)
      rescalingB = xprime - (rescalingA * x)
      rescaledValue<-z*rescalingA+rescalingB
      return(rescaledValue)
    }
    return(FUN)
}

# auxilliary and strictly internal function for converting a list of interesting gene sets to strings (which can then be coerced to colors)
assigncolor<-function(gene,geneSetsToColor)
{
  for(list in 1:length(geneSetsToColor))
  {
    if(gene %in% geneSetsToColor[[list]])
    {
      return(names(geneSetsToColor)[list])
    }
  }
  return('')
}

######## how traverse_branches was supposed to be all along
# Given two endpoints,  returns branches traversed between them. Also specifies whether the branch is in the correct orientation or reversed.
traverse_branches_by_endpoints<-function(x,y,embeddedTree){
  connectivity<-matrix(nrow = length(embeddedTree$Branches),ncol=2)
  for(i in 1:length(embeddedTree$Branches)) connectivity[i,]<-c(embeddedTree$Branches[[i]][1],embeddedTree$Branches[[i]][length(embeddedTree$Branches[[i]])])
  adj.matrix<-matrix(0,nrow=max(connectivity),ncol=max(connectivity))
  for(row in 1:nrow(connectivity)) adj.matrix[connectivity[row,1],connectivity[row,2]]<-1
  adj.net<-graph_from_adjacency_matrix(adjmatrix=adj.matrix,mode = 'undirected')
  shortestpath<-as.vector(shortest_paths(adj.net,x,y)$vpath[[1]])
  branches<-c()
  reversed<-c()
  # HORRIBLE piece of logic: tests for all of the 2-mers of shortestpath if they are the same as the edges in connectivity (or the inversion of them). 
  # isTRUE(all.equal.. is necessary because identical() cares about some aspect of 
  for(i in 1:(length(shortestpath)-1)){
    for(j in 1:nrow(connectivity)){
      if(isTRUE(all.equal(connectivity[j,],shortestpath[c(i,i+1)]))) 
      {
        branches<-c(branches,j)
        reversed<-c(reversed,F)
      }
      if(isTRUE(all.equal(connectivity[j,],shortestpath[c(i+1,i)]))) 
      {
        branches<-c(branches,j)
        reversed<-c(reversed,T)
      }
    }
  }
  
  return(list(Branches=branches,Reversed=reversed))
}

# Given two endpoints, calls traverse_branches to calculate the branches in between. Then computes the correlation (or MI) between the genes in those branches, and plots this.
makeGRN<- function(endpoint1,
                   endpoint2,
                   emb_tree,
                   cutoff=0.70,
                   quantilecutoff=NULL, # alternative way of specifying cutoff - as a quantile (e.g. only show the top 10% of edges)
                   mode=c('tree','cells'),
                   globalmode=F, # globalmode means that all nodes/cells are included in the analysis
                   minconnections=cutoff,
                   cellsubset=NULL, # input a vector of coordinates, relative to the trajectory under analysis to be included in net. mainly for sliding window net generation
                   genesubset=NULL) # a vector of genes to make a network based on
{
  
  if('tree'%in%mode) mode<-'tree'
  if(globalmode){
    Branches<-1:length(emb_tree$Branches)
    Reversed<-rep(F,length(emb_tree$Branches))
    tree<-data.frame(Branches=Branches,Reversed=Reversed)
    endpoints<-'globalmode'  
    }  
    else{ 
    tree<-traverse_branches_by_endpoints(endpoint1,endpoint2,emb_tree)
    endpoints<-c(endpoint1,endpoint2)
    }
  cellnums<-c()
  if(mode=='tree') 
  { 
    for(i in 1:length(tree$Branches))
    {
      cellstemp<-emb_tree$Branches[[tree$Branches[i]]]
      if(tree$Reversed[i]) cellstemp<-cellstemp[length(cellstemp):1]
      cellnums<-c(cellnums,cellstemp)}
  }
  cells<-emb_tree$Nodes[cellnums,]
  if(mode=='cells')
  { 
    for(i in 1:length(tree$Branches)) cellnums<-c(cellnums,which(emb_tree$Cells2Branches==tree$Branches[i]))
    cells<-emb_tree$CellCoords[cellnums,]
  }
  cells<-unique(cells)
  
  if(!is.null(cellsubset)) cells<-cells[cellsubset,]
  if(!is.null(genesubset)) cells<-cells[,genesubset[!is.na(genesubset)]]
  
  cors<-cor(cells)
  full_cors=cors
  names<-colnames(cells)
  rownames(cors)<-names
  colnames(cors)<-names
  if(!is.null(quantilecutoff)) cutoff<-quantile(cors,1-2*(1-quantilecutoff))
  cors[which(cors<cutoff)]<-0
  diag(cors)<-0                                                     
  connectedcors<-cors[,colSums(abs(cors))>=minconnections]           # removes unconnected genes
  connectedcors<-connectedcors[rowSums(abs(cors))>=minconnections,]  # removes unconnected genes

  # net<-as.network.matrix(connectedcors,matrix.type ='adjacency' ,directed = FALSE,ignore.eval = FALSE) # without intergraph - doesn't work with clustering
  net<-graph_from_adjacency_matrix(connectedcors,mode='undirected',weighted=T,diag=F)
  return(list(net=net,connectedcors=connectedcors,cutoff=cutoff,mode=mode,endpoints=endpoints,exprData=cells,cell.indexes=cellnums,full_cors=full_cors))
}  

# calculates the clustering of the network according to a method of choice, then 'declusters' any clusters below a set size
clusterGRN<-function(GRN,method='fastgreedy',min.cluster=2){
  if(method=='fastgreedy') genes2cluster<-membership(fastgreedy.community(GRN$net))               # seems good
  if(method=='edge_betweenness') genes2cluster<-membership(edge.betweenness.community(GRN$net))   # seems good
  if(method=='label_propagation') genes2cluster<-membership(label.propagation.community(GRN$net)) # seems good
  if(method=='leading_eigen') genes2cluster<-membership(leading.eigenvector.community(GRN$net))   # seems bad
  if(method=='louvain') genes2cluster<-membership(cluster_louvain(GRN$net))                       # seems bad
  if(method=='spinglass') genes2cluster<-membership(spinglass.community(GRN$net))                 # seems bad
  if(method=='walktrap') genes2cluster<-membership(walktrap.community(GRN$net))                   # seems okay
  
  clusters<-list(unclustered=NULL)
  for(i in 1:max(genes2cluster))
  {
    if(sum(genes2cluster==i)<min.cluster)
    {
      clusters[['unclustered']]<-c(clusters[['unclustered']],names(genes2cluster[which(genes2cluster==i)]))
      genes2cluster[genes2cluster==i] <- 0
    }
    else clusters[[i+1]]<-names(genes2cluster[which(genes2cluster==i)])
  }
  names(clusters) <- c('unclustered',1:(length(clusters)-1))
  clusters = clusters[-which(sapply(clusters, is.null))]
#  clusters = clusters[length(clusters):1]
  return(list(genes2cluster=genes2cluster,clusters=clusters))
}

# calculates the centrality of the network according to a method of choice
centralityGRN<-function(GRN,method='betweenness'){
  if(method=='betweenness') return(sna::betweenness(GRN$net,cmode='undirected',ignore.eval = FALSE))    # good
  if(method=='closeness') return(sna::closeness(GRN$net,cmode='undirected',ignore.eval = FALSE))        # not good
  
}
  

  
# plots the network. colorscheme takes the same inputs as the viridis 'option' parameter. GeneSetsToColor takes a list of gene sets, passes it through
# assign color to color them differently. Meant to color e.g. DE genes or interesting genes and their connected genes 
plotGRN<-function(GRN,
                  clustering=NULL,
                  centrality=NULL,
                  title=NULL,
                  colorscheme='magma',
                  plotnamesofgenes=NULL,
                  geneSetsToColor=NULL, 
                  DEgenes = NULL, # takes the output of gatherDEgenes, and colors the differentially up and downregulated on a gradient scale
                  DEdifference=F, # default for DEgenes is to color according to e-value - this argument can instead take 'difference' and then the function colors according to mean difference of expression value
                  sizeNodes =4,sizeEdges =3 ,sizeLabels = .2 ,sizeTitle = 10) { # size variables - optimal values depend on the final size of the plot
   set.seed(42)
  if(is.null(title)){
  if(GRN$endpoints=='globalmode') title<-'Whole tree, '
  else title<-paste0('Endpoints ',GRN$endpoints[1],' and ',GRN$endpoints[2],', ')
  title<-paste0(title,'mode: ',GRN$mode,', cutoff = ',GRN$cutoff)
  }
  # sets a default value to color (for 'uncolored' plots)
  color<-'1'
  
  if(!is.null(clustering)){
    color<-clustering$genes2cluster
  }
  if(!is.null(centrality)){
    centralitymean<-mean(centrality)
    centralitysd<-sd(centrality)
    centrality<-sapply(centrality,FUN = function(x){
      if(x>centralitymean+2*centralitysd) x<-centralitymean+2*centralitysd
      return(x)    
      })
    color<-centrality
    
    color.rescale <- def.rescale(min(color),max(color),1,100)
    color = ceiling(color.rescale(color))
    names(color)<-row.names(GRN$connectedcors)
  }
  
  rescaling.cors<-def.rescale(GRN$cutoff,1,.05,.35)

  plottingobject<-ggnetwork(GRN$net,weighted=T)
  colnames(plottingobject)[7]<-'weight'
  plottingobject['color']<-NA
  
  for(name in names(color)){plottingobject$color[which(plottingobject$vertex.names==name)]<-color[name]}

  if(!is.null(plotnamesofgenes)){
      plottingobject['interesting']<-''    
      for(gene in plotnamesofgenes){
        plottingobject$interesting[which(plottingobject$vertex.names==gene)]<-gene
      }
  }
  else plottingobject['interesting']<- plottingobject$vertex.names
  
  if(!is.null(geneSetsToColor))
  {
    plottingobject$color<-sapply(plottingobject$vertex.names,FUN=assigncolor,geneSetsToColor=geneSetsToColor) 
  }

  if(!is.null(DEgenes))
  {
    for(i in 1:nrow(DEgenes))
    {
      if(DEdifference==T) plottingobject$color[which(plottingobject$vertex.names==as.character(DEgenes$GeneName[i]))]<-as.numeric(DEgenes$MeanDifference[i])
      else plottingobject$color[which(plottingobject$vertex.names==as.character(DEgenes$GeneName[i]))]<-as.numeric(DEgenes$MinusLogE[i])
    }
  }
  if(is.null(DEgenes)) plottingobject$color<-as.factor(plottingobject$color)
  if(all(is.na(plottingobject$color))) plottingobject$color <- '1'

  plot<-ggplot(plottingobject, aes(x = x, y = y, xend = xend, yend = yend)) +
    geom_edges(alpha=0.2,aes(size=rescaling.cors(weight))) +
    geom_nodetext_repel(aes(label = interesting,size=sizeLabels), fontface = "bold", box.padding = unit(.25, "lines"),point.padding = unit(0.005,'cm'),segment.size = .2) +
    geom_nodes(aes(fill = color),shape=21, size = sizeNodes,stroke=0) +
    scale_color_manual(values=c('grey30','red'),guide='legend')+
    scale_size_area(max_size=sizeEdges)+
    theme_blank(plot.title = element_text(size=sizeTitle,hjust = 0),legend.position='bottom',legend.text=element_text(size=6),legend.title=element_text(size=7))+
    ggtitle(title)
  if(!is.null(DEgenes))   plot<-plot + scale_fill_gradient2(high='red',mid='black',low='deepskyblue',na.value = 'black')
  if(!is.null(clustering)) plot<- plot + scale_fill_discrete()
    
  if(!is.null(DEgenes)) plot<-plot + guides(fill=guide_colorbar(title='Rank',ticks=F),size=F)
  if(!is.null(centrality)) plot<-plot + guides(fill=guide_colorbar(title='Centrality',ticks=F),size=F)
  if(!is.null(clustering)|!is.null(geneSetsToColor)) plot<-plot+guides(fill=guide_legend(title='Clusters'),size=F)
  if(is.null(centrality)&is.null(clustering)&is.null(geneSetsToColor)&is.null(DEgenes)) plot<-plot + guides(fill=F,size=F)
  return(plot)
  }

# To do GO term enrichment first you need to load an annotation. mapping is essentially a list of genes, with corresponding GO
# terms (for example for a given organism). ID is the required type of gene ID, e.g. 'alias' or 'entrez'.

makeGOAnnotation<-function(mapping,ID){
  annot<-list()
  annot$BP<-inverseList(annFUN.org(whichOnto='BP',mapping=mapping,ID=ID))  
  annot$CC<-inverseList(annFUN.org(whichOnto='CC',mapping=mapping,ID=ID))
  annot$MF<-inverseList(annFUN.org(whichOnto='MF',mapping=mapping,ID=ID))
  return(annot)
}
# given an annotation, a clustering of a GRN it gives the top n GO terms of the ontology declared in 'ontology' (BP, MF or CC)
# for each cluster

GRN_GO_termenrichment<-function(annotation, clustering,ontology,n=10){
  result=list()
  test.stat <- new("classicCount", testStatistic = GOFisherTest, name = "Fisher test")
  
  for(i in 1:length(clustering$clusters))
  {
    geneNames<-names(annotation[[ontology]])
    universe<-geneNames[geneNames %in% names(clustering$genes2cluster)]
    geneList <- factor(as.integer(universe %in% clustering$clusters[[i]]))
    if(length(levels(geneList))!=2) next
    names(geneList) <- universe
    GOdata <- new("topGOdata", ontology = ontology, allGenes = geneList,annot = annFUN.gene2GO, gene2GO = annotation[[ontology]])
    resultFisher <- getSigGroups(GOdata, test.stat)
    result[[i]] <- GenTable(GOdata, classicFisher = resultFisher, ranksOf = "classicFisher", topNodes = 10)
  }
  return(result)
}



# function for the production of plots with ranked DE genes colored according to their rank. produces a data.frame that contains the genes
# and their ranks. takes output of branch_differential_expression or subpopulation_differential_expression
gatherDEgenes <-function(differential_expression,cutoff=10^-3)
{
  output<-data.frame(GeneName = differential_expression$GeneName, 
                     MinusLogE = -log(differential_expression$evals),
                     MeanDifference = differential_expression$MeanDifference)
  output$MinusLogE[output$MinusLogE<0]<-0
  
 # MeanDiffRescalingMax<-max(abs(range(output)))
  output$MeanDifference[which(output$MinusLogE< (-log(cutoff)))]<-0
  rescalingneg<-def.rescale(min(output$MinusLogE[which(differential_expression$MeanDifference<0)]),
                            max(output$MinusLogE[which(differential_expression$MeanDifference<0)]),
                            0,1)
  rescalingpos<-def.rescale(min(output$MinusLogE[which(differential_expression$MeanDifference>0)]),
                            max(output$MinusLogE[which(differential_expression$MeanDifference>0)]),
                            0,1)
  output$MinusLogE[which(differential_expression$MeanDifference<0)]<- -rescalingneg(output$MinusLogE[which(differential_expression$MeanDifference<0)])
  output$MinusLogE[which(differential_expression$MeanDifference>0)]<- rescalingpos(output$MinusLogE[which(differential_expression$MeanDifference>0)])
  output$Rank<-1:nrow(output)
  return(output)
}











