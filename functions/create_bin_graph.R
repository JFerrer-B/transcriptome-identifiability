##create_bin_graph


create_bin_graph <- function(SG_list,milr,ncores=1){

  Bin_Graph <- mclapply(SG_list, function(X){
    Ntrans <- length(unique(unlist(txName(X$sgf))))
    # if(X$nPaths > 100000 | X$nPaths == 0 | Ntrans > X$nPaths | Ntrans == 1){
    if(X$nPaths > 100000 | X$nPaths == 0 | Ntrans > X$nPaths){
      return(NULL)
    }else{
      adjbin <- getAdjBinGraph(X$sgf, lr = milr)
      if(nrow(adjbin$Bindata)==0){
        return(NULL)
      }
      # plotNice(adjbin$Adjacency)
      adjbin$NumRealTrans <- Ntrans
      ## Delete NULL trancripts to avoid errors
      if(sum(sapply(adjbin$TranscriptBins,length)==0)>0){
        NullTranscript<-which(sapply(adjbin$TranscriptBins,length)==0)
        adjbin$TranscriptBins<-adjbin$TranscriptBins[-NullTranscript]
        adjbin$Transcripts<-adjbin$Transcripts[-NullTranscript]
      }
      adjbin$NumRealTransInPaths <- length(adjbin$Transcripts)
      adjbin$Adjacency["S","E"]<-0
      nPaths <- solve(diag(ncol(adjbin$Adjacency)) - adjbin$Adjacency)["S","E"]
      
      adjbin <- contractGraph(adjbin = adjbin)
      adjbin$nPaths_BG <- nPaths
      adjbin$nPaths_SG <- X$nPaths
      return(adjbin)
      # plotNice(adjbin$Adjacency)
      
    }
  },mc.cores = ncores)
  
  Bin_Graph <- Filter(Negate(is.null), Bin_Graph)
  
  return(Bin_Graph)
  
}