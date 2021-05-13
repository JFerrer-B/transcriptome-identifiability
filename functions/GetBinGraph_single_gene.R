GetBinGraph_single_gene <- function(miSG,milr){
  X <- miSG
  Ntrans <- length(unique(unlist(txName(X$sgf))))
  adjbin <- getAdjBinGraph(X$sgf, lr = milr)
  adjbin$NumRealTrans <- Ntrans
  adjbin$NumRealTransInPaths <- length(adjbin$Transcripts)
  adjbin$Adjacency["S","E"]<-0
  nPaths <- solve(diag(ncol(adjbin$Adjacency)) - adjbin$Adjacency)["S","E"]
  adjbin <- contractGraph(adjbin = adjbin)
  adjbin$nPaths_BG <- nPaths
  adjbin$nPaths_SG <- X$nPaths
  
  return(adjbin)
  
}
  
