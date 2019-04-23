library(infotheo)
library(igraph)
library(entropy)

map2color<-function(x,pal,limits=NULL){
  if(is.null(limits)) limits=range(x)
  pal[findInterval(x,seq(limits[1],limits[2],length.out=length(pal)+1), all.inside=TRUE)]
}

norm_lib_size <- function(X) {
  sum1 <- apply(X, 1, sum)
  scalings <- sum1 / mean(sum1)
  X <- (1 / scalings) * X
  return (X)
}

# from http://stackoverflow.com/a/38010916/7195775
expected_mi <- function(s1, s2, l1, l2, n){    # expected mutual information
  s_emi <- 0
  for(i in 1:l1){
    for (j in 1:l2){
      min_nij <- max(1,s1[i]+s2[j]-n)
      max_nij <- min(s1[i],s2[j])
      n.ij <- seq(min_nij, max_nij)   # sequence of consecutive numbers
      t1<- (n.ij / n) * log((n.ij * n) / (s1[i]*s2[j]))
      t2 <- exp(lfactorial(s1[i]) + lfactorial(s2[j]) + lfactorial(n - s1[i]) + lfactorial(n - s2[j]) - lfactorial(n) - lfactorial(n.ij) - lfactorial(s1[i] - n.ij) - lfactorial(s2[j] - n.ij) - lfactorial(n - s1[i] - s2[j] + n.ij))
      emi <- sum(t1*t2)
      s_emi <- s_emi + emi
    }
  }
  return(s_emi)
}

# from http://stackoverflow.com/a/38010916/7195775
adjusted_mi <- function(U, V, status){
  if (all(status == 0)) {
    return (NaN)
  }
  
  U <- as.numeric(U)
  V <- as.numeric(V)
  
  # apparently NA values in U are just ignored D:
  # so we need to compensate
  U[is.na(U)] <- 0
  U[U==0] <- max(U) + 1

  U = as.factor(U)
  V = as.factor(V)
  
  s1 <- tabulate(U)
  s2 <- tabulate(V)
  l1 <- length(s1)
  l2 <- length(s2)
  N <- length(U)
  tij <- table(U, V)
  mi <- mutinformation(U, V)
  
  # avoid log(0)
  ls1 <- s1/N
  ls1[ls1==0] <- 1
  ls1 <- log(ls1)
  
  ls2 <- s2/N
  ls2[ls2==0] <- 1
  ls2 <- log(ls2)
  
  h1 <- -sum(s1*ls1)/N
  h2 <- -sum(s2*ls2)/N
  
  emi <- expected_mi(s1,s2,l1,l2,N) # EMI Expected MI
  ami <- (mi-emi)/(max(h1,h2) - emi)  #AMI Adjusted MI
  return(ami)
}


assign_status <- function(U, V) {
  # V is truth!
  # TP: pair of points in same cluster in U, V
  # TN: pair of points in diff cluster in U, V
  # FN: pair of points in diff cluster in U, same in V
  # FP: pair of points in same cluster in U, diff in V
  # when there are NA values in the predictions:
  # u=NA : FN++
  
  if (!all(is.na(U))) { # this will happen when a method is unable to run
    bu <- outer(U, U, "==")
    bv <- outer(V, V, "==")
    TP <- length(intersect(which(bu), which(bv))) - length(U)
    TN <- length(intersect(which(!bu), which(!bv)))
    FP <- length(intersect(which(bu), which(!bv)))
    FN <- length(intersect(which(!bu), which(bv)))
    FN <- FN + length(which(is.na(bu)))
    return(c(TP, TN, FN, FP)/2)
  }
  return(c(0,0,0,0))
}




randInd_manual <- function(U, V, status) {
  TP <- status[1]
  TN <- status[2]
  FN <- status[3]
  FP <- status[4]
  
  res <- (TP+TN) / (TP+FP+TN+FN)
  return(res)
}

# adjusted_randInd <- function(status) {
#   TP <- status[1]
#   TN <- status[2]
#   FN <- status[3]
#   FP <- status[4]

#   res <- 2*(TN*TP − FP*FN) / ((TN+FP)*(FP+TP) + (TN+FN)*(FN+TP))

#   return(res)
# }

matthews_cor <- function(U, V, status) {
  TP <- status[1]
  TN <- status[2]
  FN <- status[3]
  FP <- status[4]
  
  res <- (TP*TN - FP*FN) / sqrt((TP+FP) * (TP+FN) * (TN+FP) * (TN+FN))
  return(res)
}

f_measure <- function(U, V, status, b=1) {
  TP <- status[1]
  TN <- status[2]
  FN <- status[3]
  FP <- status[4]
  
  P <- TP / (TP+FP)
  R <- TP / (TP+FN)
  res <- (b^2 + 1) * P * R / (b^2 * P + R)
  return(res)
}


jaccard <- function(U, V, status) {
  TP <- status[1]
  TN <- status[2]
  FN <- status[3]
  FP <- status[4]
  
  res <- TP / (TP+FP+FN)
}


fowkles_mallows <- function(U, V, status) {
  TP <- status[1]
  TN <- status[2]
  FN <- status[3]
  FP <- status[4]
  
  P <- TP / (TP+FP)
  R <- TP / (TP+FN)
  res <- sqrt(P*R)
}


assign_days <- function(N, d) {
  day <- rep(0, N)
  for (i in 1:d) {
    from <- (i-1)*N/d + 1
    to <- i*N/d
    day[from:to] = i
  }
  return(day)
}


goodman_kruskal_weights <- function(a, b) {
  l <- length(a)
  
  a_width <- max(a) - min(a)
  b_width <- max(b) - min(b)

  wA <- matrix(0, ncol = l, nrow = l)
  wB <- matrix(0, ncol = l, nrow = l)

  if (a_width>0) {
    wA <- outer(a, a, '-') / a_width
  }
  if (b_width>0) {
    wB <- outer(b, b, '-') / b_width
  }

  W <- array(0, dim=c(l, l, 2))
  sA <- sign(wA)
  sB <- sign(wB)

  same <- (sA*sB == 1)
  opp <- (sA*sB == -1)
  null <- (sA == 0 & sB == 0)

  W[,,1] = wA / wB
  W[,,2] = wB / wA

  res <- array(0, dim=c(l,l))

  res[same] = pmin(W[,,1], W[,,2])[same]
  res[opp] = pmax(W[,,1], W[,,2])[opp]
  res[null] = 1

  return(res)
}


goodman_kruskal_signs <- function(a, b) {
  l <- length(a)
  
  wA <- matrix(0, ncol = l, nrow = l)

  wA <- sign(outer(a, a, '-'))
  wB <- sign(outer(b, b, '-'))

  uneq <- (wB != 0)
  null <- (wA == 0 & wB == 0)

  res <- array(0, dim=c(l,l))

  res[uneq] = (wA/wB)[uneq]
  res[null] = 1

  return(res)
}



goodman_kruskal_index <- function(a, b, weighted=TRUE) {
  if (length(a) == 1) {
    return(NaN)
  }

  if (weighted) {
    W <- goodman_kruskal_weights(a, b)
  } else {
    W <- goodman_kruskal_signs(a, b)
  }

  tri <- upper.tri(W, diag=FALSE)

  up <- sum(W[tri])
  down <- sum(abs(W)[tri])
  return(up/down)
}


kendall_index <- function(a, b, weighted=TRUE) {
  if (length(a) == 1) {
    return(NaN)
  }

  if (weighted) {
    W <- goodman_kruskal_weights(a, b)
  } else {
    W <- goodman_kruskal_signs(a, b)
  }

  tri <- upper.tri(W, diag=FALSE)

  up <- sum(W[tri])
  n <- dim(W)[1]
  down <- n*(n-1)/2

  return(up/down)
}


hor_colorbar <- function(lut, min, max) {
  scale = (length(lut)-1)/(max-min)
  
  plot(c(min,max), c(0,10), type='n', bty='n', xaxt='n', xlab='', yaxt='n', ylab='')
  title(main="pseudotime", cex=1.5)
  
  axis(1, at=c(0, N), labels=c("0", "max(pt)"), lwd=0, lwd.ticks=1, cex.axis=1.5)
  for (i in 1:(length(lut)-1)) {
    y = (i-1)/scale + min
    rect(y,0,y+1/scale,10, col=lut[i], border=NA)
  }
  rect(min,0,max,10, border="black", lwd = 2)
}


parse_benchmark <- function(subfolders, benchmark_dir, functions, methods, type, replNA = TRUE) {
  fac <- length(functions)
  met <- length(methods)
  res <- matrix(0, ncol=length(methods)+2,
                nrow=length(subfolders)*fac)
  
  colnames(res) <- c(methods, "measure", "experiment")
  
  if(type == "branch") {
    def <- c(0,-1,0,0,0,0,0)
  } else {
    def <- c(-1, -1, -1, -1, -1)
  }
  
  default <- matrix(rep(def, met), ncol=met)
  
  for (i in 1:length(subfolders)) {
    r_eval = matrix(rep(def, met+2), ncol=met+2)
    colnames(r_eval) <- c(methods, "measure", "experiment")

    # print(i)
    s <- subfolders[i]
    eval <- paste(benchmark_dir, s, "/", s, "_eval.txt", sep="")
    s_eval <- t(as.matrix(read.table(eval, check.names=FALSE, stringsAsFactors = FALSE)))
    keep <- rownames(s_eval) %in% functions
    s_eval <- s_eval[keep,]
    
    keep <- colnames(s_eval) %in% methods
    s_eval <- s_eval[,keep]
    if(replNA) {
      default <- matrix(rep(def, ncol(s_eval)), ncol=ncol(s_eval))
      repl <- is.na(s_eval)
      s_eval[repl] <- default[repl]
    }
    
    for (m in colnames(s_eval)) {
      r_eval[,m] = s_eval[,m]
    }
    r_eval[,"measure"] = rownames(s_eval)
    r_eval[,"experiment"] <- subfolders[i]
    # s_eval <- cbind(s_eval, names, exp)
    start <- (i-1)* fac + 1
    end <- start + fac - 1
    
    # print(dim(s_eval))
    # print(dim(res[start:end,]))
    # print("======")
    res[start:end, 1:dim(r_eval)[2]] <- r_eval
  }
  return(res)
}


parse_cluster_output <- function(subfolders, benchmark_dir, functions, methods, type, replNA = TRUE) {
  fac <- length(functions)
  met <- length(methods)
  res <- matrix(0, ncol=length(methods)+2,
                nrow=length(subfolders)*fac)
  
  colnames(res) <- c(methods, "measure", "experiment")
  
  if(type == "branch") {
    def <- c(0,-1,0,0,0,0,0)
  } else {
    def <- c(-1, -1, -1, -1, 0)
  }
  
  default <- matrix(rep(def, met), ncol=met)
  
  for (i in 1:length(subfolders)) {
    r_eval = matrix(rep(def, met+2), ncol=met+2)
    colnames(r_eval) <- c(methods, "measure", "experiment")

    # print(i)
    s <- subfolders[i]
    eval <- paste(benchmark_dir, s, sep="")
    s_eval <- t(as.matrix(read.table(eval, check.names=FALSE, stringsAsFactors = FALSE)))
    keep <- rownames(s_eval) %in% functions
    s_eval <- s_eval[keep,]
    
    keep <- colnames(s_eval) %in% methods
    s_eval <- s_eval[,keep]
    if(replNA) {
      default <- matrix(rep(def, ncol(s_eval)), ncol=ncol(s_eval))
      repl <- is.na(s_eval)
      s_eval[repl] <- default[repl]
    }
    
    for (m in colnames(s_eval)) {
      r_eval[,m] = s_eval[,m]
    }
    r_eval[,"measure"] = rownames(s_eval)
    r_eval[,"experiment"] <- subfolders[i]
    # s_eval <- cbind(s_eval, names, exp)
    start <- (i-1)* fac + 1
    end <- start + fac - 1
    
    # print(dim(s_eval))
    # print(dim(res[start:end,]))
    # print("======")
    res[start:end, 1:dim(r_eval)[2]] <- r_eval
  }
  return(res)
}


parse_grid <- function(subfolders, benchmark_dir, funcnames, timenames, params) {
  fac <- length(funcnames) + length(timenames)
  length_m <- prod(unlist(lapply(params, length)))
  
  def <- c(0,-1,0,0,0,0,0,-1,-1,-1,-1, 0,0,0, "test")
  res <- t(matrix(rep(def, length(subfolders)*length_m), nrow=fac+length(params)+1))
  
  colnames(res) <- c(funcnames, timenames, "N_yk", "mu0", "lambda0", "experiment")
  
  for (i in 1:length(subfolders)) {
    s <- subfolders[i]
    print(s)
    files <- list.files(path=paste(benchmark_dir, s, sep=""), pattern = "*eval")
    s_eval <- t(matrix(rep(def[1:14], length_m), nrow=fac+length(params)))
    
    for (j in seq_along(files)) {
      eval <- paste(benchmark_dir, s, "/", files[j], sep="")
      # print(eval)
      s_eval[j,] <- as.matrix(read.table(eval, check.names=FALSE))
    }
    
    exp <- rep(subfolders[i], length_m)
    s_eval <- cbind(s_eval, exp)
    start <- (i-1)* length_m + 1
    end <- start + length_m - 1
    
    # print(dim(s_eval))
    # print(dim(res[start:end,]))
    # print("======")
    res[start:end,] <- s_eval
  }
  res <- unique(res)
  return(res)
}


read_parameters <- function(subfolders, benchmark_dir) {
  res <- matrix(0, nrow = length(subfolders), ncol = 6)
  # colnames(res) = c("alpha", "beta")

  for (i in 1:length(subfolders)) {
    s <- subfolders[i]
    # print(s)
    # paramfile <- paste(benchmark_dir, "/", s, "/", s, "_geneparams.txt", sep = "")
    # params <- read.csv(paramfile, sep = "\t")
    # res[i, "alpha"] = mean(params$alpha)
    # res[i, "beta"] = mean(params$beta)
    sim <- paste(benchmark_dir, "/", s, "/", s, "_simulation.txt", sep = "")
    Dataset <- read.csv(sim, sep = "\t")
    X <- data.matrix(Dataset[2:dim(Dataset)[2]])

    gene_entropy = apply(X, 1, function(x) entropy(x, numBins=100))
    res[i, 1] = sum(gene_entropy)
    dX <- discretize(X, numBins=100)
    res[i, 2] <- entropy(dX)

    X <- norm_lib_size(X)
    gene_entropy = apply(X, 1, memb_ent)
    res[i, 3] = sum(gene_entropy)
    dX <- discretize(X, numBins=100)
    res[i, 4] <- entropy(dX)

    logX <- log(X + 1)
    gene_entropy = apply(logX, 1, memb_ent)
    res[i, 5] = sum(gene_entropy)
    dX <- discretize(logX, numBins=100)
    res[i, 6] <- entropy(dX)
  }
  return (res)
}

memb_ent <- function(x) {
  d = discretize(x, numBins=100)
  entropy(d)
}