getDistances <- function(sgf) {
  Adj <- getAdj(sgf)
  # Split the nodes in type a and b
  zeroMatrix <- Adj * 0;
  # Assign distances
  lengthexons <- width(sgf[type(sgf)=="E"])
  names(lengthexons) <- featureID(sgf[type(sgf)=="E"])
  eye <- diag(c(1,lengthexons+1,1))
  BigAdjacency <- rbind(cbind(zeroMatrix, eye),cbind(Adj, zeroMatrix))
  colnames(BigAdjacency) <- rownames(BigAdjacency) <- c(paste(colnames(Adj),"a"), paste(colnames(Adj),"b"))
  g <- graph_from_adjacency_matrix(BigAdjacency, weighted = TRUE)
  E(g)$weight <- E(g)$weight -1
  # Assign distances
  minDistances <- distances(g,weights = E(g)$weight, mode = "out", algorithm = "johnson")
  maxDistances <- -distances(g,weights = -E(g)$weight, mode = "out", algorithm = "johnson")
  Dejar <- grep("a", colnames(minDistances))
  minDistances <- minDistances[Dejar, Dejar]
  maxDistances <- maxDistances[Dejar, Dejar]
  colnames(minDistances) <-  rownames(minDistances) <- gsub(" a", "", colnames(minDistances))
  colnames(maxDistances) <-  rownames(maxDistances) <- gsub(" a", "", colnames(maxDistances))
  maxDistances[maxDistances==-Inf] <-  Inf
  # minDistances <- pmin(minDistances, t(minDistances))
  # minDistances[lower.tri(minDistances)] <- -minDistances[lower.tri(minDistances)] 
  # maxDistances <- pmin(maxDistances, t(maxDistances))
  # maxDistances[lower.tri(maxDistances)] <- -maxDistances[lower.tri(maxDistances)] 
  return(list(minDistances=minDistances,maxDistances=maxDistances,lengthexons=lengthexons,Adj=Adj))
}
  
  
