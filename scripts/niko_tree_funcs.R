nikoElastic <- function(ScaffoldTree,
                                 N_yk = 100,
                                 lambda = 0.01,
                                 mu = 0.1,
                                 FixEndpoints = F,
                                 plot = F,
                                 NBranchScaffoldNodes = 1,
                                 NCores = 1)
{
  # Testing
  # Default parameters taken from adjustment in real datasets with 100 N_yks
  # N_yk=100
  # lambda_0=2.03e-09
  # mu_0=0.00625
  #  End testing
  
  Coords=c()
  Edges=c()
  
  # scaling according to N_yk
  # mu=(N_yk-1)*mu_0
  # lambda=((N_yk-2)**3)*lambda_0
  
  # I apply unique in case there are repeated nodes, which is a consequence of having trifurcations or higher order connections
  TopologyNodes=unique(c(ScaffoldTree$Endpoints, ScaffoldTree$Branchpoints))
  # Calculate connectivity in between branches
  TopologyCoords=ScaffoldTree$CellCoordinates[TopologyNodes,]
  # Given initial topology, calculate connectivity in between branches
  TopologyNodes=c(ScaffoldTree$Endpoints, ScaffoldTree$Branchpoints)
  TopologyEdges=c()
  
  if(length(ScaffoldTree$Endpoints)==2)
  {
    TopologyEdges=rbind(TopologyEdges, (c(which(TopologyNodes==ScaffoldTree$Branches[1])[1], which(TopologyNodes==ScaffoldTree$Branches[2])[1])))
  }else
  {
    for(i in (1:dim(ScaffoldTree$Branches)[1]))
    {
      TopologyEdges=rbind(TopologyEdges, (c(which(TopologyNodes==ScaffoldTree$Branches[i,1])[1], which(TopologyNodes==ScaffoldTree$Branches[i,2])[1])))
    }
  }
  
  
  TopologyCoordsAux=TopologyCoords
  TopologyEdgesAux=TopologyEdges
  TopologyNodesAux=TopologyNodes
  
  # ---------------------------------------
  # -------Add branch middle scaffold nodes
  # ---------------------------------------
  if(NBranchScaffoldNodes)
  {
    TopologyEdgesAux=c()
    
    for(i in 1:dim(ScaffoldTree$Branches)[1])
    {
      # Distance between to branch extreme cells
      extremepointsDistance=ScaffoldTree$DijkstraDistances[ScaffoldTree$Branches[i,1], ScaffoldTree$Branches[i,2]]
      # Cells between those two branches
      nodesBranchi=ScaffoldTree$SkeletonNodes[which(ScaffoldTree$SkeletonNodes %in% ScaffoldTree$Cells2BranchesAssignments[[i]])]
      # take the endpoints out of the array
      nodesBranchi=nodesBranchi[which(!(nodesBranchi %in% ScaffoldTree$Branches[i,]))]
      
      minDistanceCenter=min((ScaffoldTree$DijkstraDistances[ScaffoldTree$Branches[i,1],nodesBranchi])-(extremepointsDistance/2))
      Distances2Nodes=ScaffoldTree$DijkstraDistances[ScaffoldTree$Branches[i,1],nodesBranchi]
      
      # New ScaffoldNode
      NewScaffoldNode=nodesBranchi[which(abs(Distances2Nodes-(extremepointsDistance/2)) == min(abs(Distances2Nodes-extremepointsDistance/2)))]
      TopologyCoordsAux=rbind(TopologyCoordsAux, ScaffoldTree$CellCoordinates[NewScaffoldNode,])
      NewScaffoldNodePosition=length(TopologyNodesAux)+1
      TopologyNodesAux=c(TopologyNodesAux, NewScaffoldNode)
      # NedEdges
      Edge1=c(which(TopologyNodesAux==ScaffoldTree$Branches[i,1]), NewScaffoldNodePosition)
      Edge2=c(which(TopologyNodesAux==ScaffoldTree$Branches[i,2]), NewScaffoldNodePosition)
      TopologyEdgesAux=rbind(TopologyEdgesAux, Edge1)
      TopologyEdgesAux=rbind(TopologyEdgesAux, Edge2)
      rownames(TopologyEdgesAux)<-NULL
    }
  }
  
  ElasticTree=ElPiGraph.R::computeElasticPrincipalCurve(X = ScaffoldTree$CellCoordinates,
                                                        NumNodes = N_yk,
                                                        InitNodePositions = TopologyCoordsAux,
                                                        InitEdges = TopologyEdgesAux,
                                                        Lambda = lambda,
                                                        Mu = mu,
                                                        Do_PCA = F,
                                                        verbose = F,
                                                        drawAccuracyComplexity = F,
                                                        drawPCAView = F,
                                                        drawEnergy = F,
                                                        Mode = 1,
                                                        n.cores = NCores)
  
  # Unlist the ElasticTree structure
  ElasticTree=ElasticTree[[1]]
  names(ElasticTree)[1]="Nodes"
  ElasticTree$Edges=ElasticTree$Edges$Edges
  
  # Add the topology in the tree structure
  ETreeTopology= list(Endpoints=seq(from=1, to=length(ScaffoldTree$Endpoints), by=1), Branchpoints=seq(from=(length(ScaffoldTree$Endpoints) +1), to= (length(ScaffoldTree$Endpoints)) + length(ScaffoldTree$Branchpoints), by=1))
  ElasticTree=c(ElasticTree, Topology=1)
  ElasticTree$Topology <- ETreeTopology
  
  # Coords for the topology elements
  ElasticTree=c(ElasticTree, EndpointsCoords=1)
  ElasticTree$EndpointsCoords=ElasticTree$Nodes[ElasticTree$Topology$Endpoints,]
  ElasticTree=c(ElasticTree, BranchpointsCoords=1)
  ElasticTree$BranchpointsCoords=ElasticTree$Nodes[ElasticTree$Topology$Branchpoints,]
  
  # Add Input Coordinates to the tree structure
  ElasticTree=c(ElasticTree, CellCoords=1)
  ElasticTree$CellCoords <- ScaffoldTree$CellCoordinates
  
  # Add Tree Connectivity
  ElasticTree=c(ElasticTree, Connectivity=1)
  ElasticTree$Connectivity  <- Edges
  
  # Assign nodes from the embedded tree to the different branches and save the structure in the tree structure
  BranchesNodes=list()
  EdgesTree=matrix(0, N_yk, N_yk)
  
  for(i in (1: dim(ElasticTree$Edges)[1]))
  {
    EdgesTree[ElasticTree$Edges[i, 1], ElasticTree$Edges[i, 2]]=1
    EdgesTree[ElasticTree$Edges[i, 2], ElasticTree$Edges[i, 1]]=1
  }
  Graph_yk=igraph::graph_from_adjacency_matrix(EdgesTree)
  for(i in 1:dim(TopologyEdges)[1])
  {
    path_brach_i=igraph::get.shortest.paths(Graph_yk,from = TopologyEdges[i,1], to = TopologyEdges[i,2])
    BranchesNodes[[i]]=path_brach_i$vpath[[1]]
  }
  
  cell2yk=c()
  
  if (length(ScaffoldTree$Endpoints)==2)
  {
    cell_i=matrix(ScaffoldTree$CellCoordinates[1,], nrow=1)
    dist_cell_i=as.matrix(stats::dist(rbind(cell_i, ElasticTree$Nodes), method = "euclidean", diag = FALSE, upper = TRUE, p = 2))
    #find the closest yk index. Decrease the index in 1, since the 1 element is the element itself
    closest_yk=sort(dist_cell_i[,1], index.return=T)$ix[2]-1
    cell2yk=rbind(cell2yk, c(1, closest_yk))
  } else {
    for (i in 1:dim(ScaffoldTree$CellCoordinates)[1])
    {
      cell_i=matrix(ScaffoldTree$CellCoordinates[i,], nrow=1)
      dist_cell_i=as.matrix(stats::dist(rbind(cell_i, ElasticTree$Nodes), method = "euclidean", diag = FALSE, upper = TRUE, p = 2))
      #find the closest yk index. Decrease the index in 1, since the 1 element is the element itself
      closest_yk=sort(dist_cell_i[,1], index.return=T)$ix[2]-1
      cell2yk=rbind(cell2yk, c(i, closest_yk))
    }
  }
  
  
  # Add cells 2 yks mapping
  ElasticTree=c(ElasticTree, Cells2TreeNodes=1)
  ElasticTree$Cells2TreeNodes <- cell2yk
  
  # assign cells to the branches
  allnodes=c()
  for(i in 1:length(BranchesNodes))
  {
    allnodes=c(allnodes, BranchesNodes[[i]])
  }
  allnodes=sort(unique(allnodes))
  cells_branchs_assigments=c()
  for (i in 1:length(BranchesNodes))
  {
    cells_branchs_assigments[which(cell2yk[,2] %in% BranchesNodes[[i]])]=i
  }
  
  ElasticTree=c(ElasticTree, Branches=1)
  ElasticTree$Branches <- BranchesNodes
  
  ElasticTree=c(ElasticTree, Cells2Branches=1)
  ElasticTree$Cells2Branches <- cells_branchs_assigments
  
  # Fixing endpoints coordinates
  if(FixEndpoints==T)
  {
    ElasticTree$Nodes[1:length(ScaffoldTree$Endpoints),]=ScaffoldTree$CellCoordinates[ScaffoldTree$Endpoints,]
  }
  
  return (ElasticTree)
}





nikoEmbed <- function(ExpressionMatrix, ElasticTree,  lambda=0.01, mu=0.1, NCores=1)
{
  # The number of nodes for the embedding tree is the same as the ones for the input low dimensional one
  N_yk=dim(ElasticTree$Nodes)[1]
  
  # scaling according to N_yk and multiply the values by the increase factor
  # mu=(N_yk-1) * mu_0 * increaseFactor_mu
  # lambda=((N_yk-2)**3) * lambda_0 * increaseFactor_lambda
  
  # Map each xn (cells) to the closest y_k (tree node)
  cell2yk=c()
  for (i in 1:dim(ExpressionMatrix)[1])
  {
    cell_i=matrix(ElasticTree$CellCoords[i,], nrow = 1)
    dist_cell_i=as.matrix(stats::dist(rbind(cell_i, ElasticTree$Nodes), method = "euclidean", diag = FALSE, upper = TRUE, p = 2))
    #find the closest yk index. Decrease the index in 1, since the 1 element is the element itself
    closest_yk=sort(dist_cell_i[,1], index.return=T)$ix[2]-1
    cell2yk=rbind(cell2yk, c(i, closest_yk))
  }
  
  # count number of xn per y_k
  yk_counts <- rep(0, N_yk)
  yk_profiles <- matrix(data = 0, ncol = dim(ExpressionMatrix)[2], nrow = N_yk)
  colnames(yk_profiles) <- colnames(ExpressionMatrix)
  # calculate the transcriptional profile for y_ks
  for (i in 1:dim(ElasticTree$Nodes)[1])
  {
    # count cells associated to y_ki
    yk_counts[i] <- length(cell2yk[which(cell2yk[,2]==i),1])
    if (yk_counts[i] > 1) {
      yk_profiles[i,] <- colMeans(ExpressionMatrix[which(cell2yk[,2]==i),])
    } else if (yk_counts[i] == 1) {
      yk_profiles[i,] <- unlist(ExpressionMatrix[which(cell2yk[,2]==i),])
    }
  }
  
  # set NA values for y_k with no mapped x_n to 0
  yk_profiles[yk_counts == 0, ] <- 0
  
  CellCoordinates=data.matrix(ExpressionMatrix)
  InitialNodesCoordinates=yk_profiles
  InitialEdges= ElasticTree$Edges
  
  EmbeddedTree=ElPiGraph.R::computeElasticPrincipalCurve(
    X = CellCoordinates, NumNodes = N_yk,
    InitNodePositions = InitialNodesCoordinates, InitEdges = InitialEdges,
    Do_PCA = F, verbose = T, drawAccuracyComplexity = F,
    drawPCAView = F, drawEnergy = F, Lambda = lambda, Mu = mu, Mode = 1, n.cores = NCores)
  
  # Unlist EmbeddedTree structure
  EmbeddedTree=EmbeddedTree[[1]]
  names(EmbeddedTree)[1]="Nodes"
  EmbeddedTree$Edges=EmbeddedTree$Edges$Edges
  
  # reassigning cells to nodes in the full dimensional space
  cell2yk_post=c()
  for (i in 1:dim(ExpressionMatrix)[1])
  {
    cell_i= matrix(ExpressionMatrix[i,], nrow=1)
    dist_cell_i=as.matrix(stats::dist(rbind(cell_i, EmbeddedTree$Nodes), method = "euclidean", diag = FALSE, upper = TRUE, p = 2))
    #find the closest yk index. Decrease the index in 1, since the 1 element is the element itself
    closest_yk=sort(dist_cell_i[,1], index.return=T)$ix[2]-1
    cell2yk_post=rbind(cell2yk_post, c(i, closest_yk))
  }
  
  allnodes=1:N_yk
  cells_branchs_assigments=c()
  for (i in 1:length(ElasticTree$Branches))
  {
    cells_branchs_assigments[which(cell2yk_post[,2] %in% ElasticTree$Branches[[i]])]=i
  }
  
  # Add tree topology to the tree structure
  EmbeddedTree=c(EmbeddedTree, Topology=1)
  EmbeddedTree$Topology <- ElasticTree$Topology
  
  # Add tree connectivity to the tree structure
  EmbeddedTree=c(EmbeddedTree, Connectivity=1)
  EmbeddedTree$Connectivity <- ElasticTree$Connectivity
  
  # Add gene names to columns
  colnames(EmbeddedTree$Nodes) <-colnames(ExpressionMatrix)
  
  # Add Input Coordinates to the tree structure
  EmbeddedTree=c(EmbeddedTree, CellCoords=1)
  EmbeddedTree$CellCoords <- ExpressionMatrix
  
  # Add cells 2 yks mapping
  EmbeddedTree=c(EmbeddedTree, Cells2TreeNodes=1)
  EmbeddedTree$Cells2TreeNodes <- cell2yk_post
  
  # Add cells 2 branches mapping
  EmbeddedTree=c(EmbeddedTree, Cells2Branches=1)
  EmbeddedTree$Cells2Branches <- cells_branchs_assigments
  
  # Assigning branches from input tree to this one
  EmbeddedTree=c(EmbeddedTree, Branches=1)
  EmbeddedTree$Branches <- ElasticTree$Branches
  
  # Calculate diffusion map for resulting tree taking out vectors equal to 0
  yk_profiles_nozeros=yk_profiles[which(!(seq(1,N_yk,1) %in% which(rowMeans((yk_profiles))==0))),]
  
  # testing
  EmbeddedTree=c(EmbeddedTree, AveragedNodes=1)
  EmbeddedTree$AveragedNodes <- yk_profiles
  
  return (EmbeddedTree)
}

rescale <- function(x, ref) {
  scale_factor <- (max(x) - min(x)) / (max(ref) - min(ref))
  x <- x - min(x)
  array(scale(x, scale_factor, center=FALSE)) + min(ref)
}

match_clusters <- function(truth, pred, cutoff=0.7) {
  id_truth <- sort(unique(truth))
  id_pred <- sort(unique(pred))
  
  test <- sapply(id_truth, FUN = function(x) {
    cluster_truth <- which(truth == x)
    sapply(id_pred, FUN = function(y) {
      cluster_pred <- which(pred == y)
      return(length(intersect(cluster_truth, cluster_pred)) / length(cluster_truth))
    })
  })
  
  best_match <- apply(test, 2, function(x) which(x > cutoff))
  corresponding <- lapply(best_match, function(x) {
    if (length(x) == 0) {
      return(NA)
    } else {
      return(id_pred[x])
    }
  })
  unlist(corresponding)
}

inflate_scaffold_tree <- function(ElasticTree, FullCoordinates) {
  ElasticTree$CellCoordinates <- FullCoordinates
  nodelist <- unlist(ElasticTree$Nodes)
  NodeMatrix <- ElasticTree$CellCoordinates[nodelist,]
  # map the cells to nodes
  cell2yk=c()
  for (i in 1:dim(ElasticTree$CellCoordinates)[1]) {
    cell_i=matrix(ElasticTree$CellCoordinates[i,], nrow=1)
    dist_cell_i=as.matrix(stats::dist(rbind(cell_i, NodeMatrix), method = "euclidean", diag = FALSE, upper = TRUE, p = 2))
    #find the closest yk index. Decrease the index in 1, since the 1 element is the element itself
    closest_yk=sort(dist_cell_i[,1], index.return=T)$ix[2]-1
    cell2yk=rbind(cell2yk, c(i, nodelist[closest_yk]))
  }
  ElasticTree$Cells2TreeNodes <- cell2yk
  
  # now map cells to branches based on their node assignment
  newbranches <- rep(0, dim(ElasticTree$CellCoordinates)[1])
  
  for (i in seq_along(ElasticTree$Nodes)) {
    # take branch i and all the nodes that belong to it
    nodes_i <- unlist(ElasticTree$Nodes[[i]])
    # take all cells that map to any of these nodes
    cells_i <- ElasticTree$Cells2TreeNodes[,2] %in% nodes_i
    newbranches[cells_i] <- i
  }
  ElasticTree$Cells2Branches <- newbranches
  
  return(ElasticTree)
}

run_niko_tree <- function(job, reduced_coords, full_coords, start, funcname) {
  cell_params <- read.table(file = paste(job, "resampled.txt", sep = "_"), sep = "\t", header = T, row.names = 1)
  par_loc <- paste(job, "params.txt", sep = "_")
  N <- dim(reduced_coords)[1]
  scaf <- CalculateScaffoldTree(reduced_coords, BranchMinLengthSensitive = sqrt(N), python_location="~/miniconda3/envs/py36/bin/python")
  if (length(scaf$Endpoints) == 2) {
      result <- evaluate_method(funcname, 0, 0, cell_params, par_loc)
      return(result)
  }
  elastic <- nikoElastic(scaf)
  elastic <- inflate_elastic_tree(elastic, full_coords)
  pseudotime <- CalculatePseudotimes(elastic, C0=start)

  pred_pseudotime <- pseudotime$Times_cells
  pred_branches <- elastic$Cells2Branches

  result <- evaluate_method(funcname, pred_branches, pred_pseudotime, cell_params, par_loc)
  return(result)
}

run_normal_tree <- function(job, reduced_coords, full_coords, start, funcname) {
  cell_params <- read.table(file = paste(job, "resampled.txt", sep = "_"), sep = "\t", header = T, row.names = 1)
  par_loc <- paste(job, "params.txt", sep = "_")
  N <- dim(reduced_coords)[1]
  scaf <- CalculateScaffoldTree(reduced_coords, BranchMinLengthSensitive = sqrt(N), python_location="~/miniconda3/envs/py36/bin/python")
  if (length(scaf$Endpoints) == 2) {
      result <- evaluate_method(funcname, 0, 0, cell_params, par_loc)
      return(result)
  }
  elastic <- CalculateElasticTree(scaf)
  elastic <- inflate_elastic_tree(elastic, full_coords)
  pseudotime <- CalculatePseudotimes(elastic, C0=start)

  pred_pseudotime <- pseudotime$Times_cells
  pred_branches <- elastic$Cells2Branches

  result <- evaluate_method(funcname, pred_branches, pred_pseudotime, cell_params, par_loc)
  return(result)
}

run_merlot <- function(job, full_coords, start, funcname) {
  cell_params <- read.table(file = paste(job, "resampled.txt", sep = "_"), sep = "\t", header = T, row.names = 1)
  par_loc <- paste(job, "params.txt", sep = "_")
  N <- dim(full_coords)[1]
  scaf <- CalculateScaffoldTree(full_coords, BranchMinLengthSensitive = sqrt(N), python_location="~/miniconda3/envs/py36/bin/python")
  if (length(scaf$Endpoints) == 2) {
      result <- evaluate_method(funcname, 0, 0, cell_params, par_loc)
      return(result)
  }
  elastic <- CalculateElasticTree(scaf)
  # elastic <- inflate_elastic_tree(elastic, full_coords)
  pseudotime <- CalculatePseudotimes(elastic, C0=start)

  pred_pseudotime <- pseudotime$Times_cells
  pred_branches <- elastic$Cells2Branches

  result <- evaluate_method(funcname, pred_branches, pred_pseudotime, cell_params, par_loc)
  return(result)
}