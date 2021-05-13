GetConvolvedMatrixRank <- function(AA,FlData,lr) {
  
  ## Function that with a given bin graph and a set of parameters of fl and lr
  ## returns the rank of the convoluted matrix for single fl and fl combinations
  if( (AA$nPaths_BG < 1000) & (nrow(AA$Bindata)^2 < 100000)){
    adjbin <- AA
    
    
    
    # FlData <- FlData_2
    # lr <- milr
    
    ##### 1. Initialization ######################################################## 
    ## Initial data
    
    Adjacency<-adjbin$Adjacency 
    Bindata<-adjbin$Bindata
    rownames(Bindata)<-adjbin$Bindata$names
    Pathnum<-as.numeric(adjbin$nPaths_BG)
    PairedM2<-vector("list",length =  round(Pathnum))## Organized as vector
    length(PairedM2)
    
    g <-  graph_from_adjacency_matrix(Adjacency)
    possibleTranscripts <- all_simple_paths(g,"S","E")
    ## Variable initialization
    
    Allbins<-vector() ## All possible bin pairs will be saved here
    TBins<-vector("list", Pathnum) ## Save specific bin pairs for transcript
    Tlength<-vector("list", Pathnum) ## Transcript lengths
    posBin<-vector("list",Pathnum) ## List with transcript bin positions
    # RankV<-vector("numeric",ncol(FlM))
    # AllConvolved<-vector("list",Flsize[2])
    # MissedTranscripts<-vector("numeric",Flsize[2])
    
    # names(RankV)<-colnames(FlM) #esto no tiene sentido... esta mal
    
    
    
    
    ##### 2. Bin pair length Calculation ###################################################### 
    
    minicounter <- 0
    # for(FlID in 1:ncol(FlData$fl)){ ## loop through all the given fl lengths
    FlID <- 1
    
    for (FragmentID in 1:nrow(FlData$fl)){ ## loop through the distribution of the given fl length
      # FragmentID <- 6
      Flvalue<-FlData$fl[FragmentID,FlID] ## Fragment length
      Dens<-FlData$weight[FragmentID,FlID] ## Normalized density function value
      
      for (TranscriptID in 1:Pathnum) { ## loop through all the paths for a given gene and fl
        
        # TranscriptID <- 1
        ##### 2.1 Parameters and initialization ########################################
        ## Number of bins
        Pathlength<-length(names(possibleTranscripts[[TranscriptID]])) 
        Pathind<-2:(Pathlength-1) ## Indices used to eliminate the S and the E
        ## Possible bins in the transcript
        TranscriptBins<-names(possibleTranscripts[[TranscriptID]])[c(Pathind)] 
        l <- length(TranscriptBins)
        posF<-c(0,cumsum(l:2)) ## Variable to store the values correctly on v
        ## indexes to save in FinalM
        # indF<-cumsum(c(0,l:2))
        ## Temporal vector to save obtained bin values
        v<-vector(mode="numeric", length=(l*(l+1))/2) #Vector of bin lengths = combinacion de n elementos tomados de 2 en 2 con rep
        
        if(minicounter<Pathnum){ # esto valores no dependen de mu_fl ni de fl ni nada... son solo del Bin GRAPH.
          ind<-match(TranscriptBins,Bindata[[1]]) #Indices of bins in this transcript
          Tlength[[TranscriptID]]<-sum(Bindata[ind,2]) #Length of transcript
          posBin[[TranscriptID]]<-cumsum(Bindata[TranscriptBins,2])
          minicounter <- minicounter + 1
        }
        
        ## Position of the paired reads end and beginning
        posP<-c(1,Flvalue-lr+1)
        
        if(Flvalue<Tlength[[TranscriptID]]+lr-1){ ## In case fl is not longer than the whole potential transcript
          
          posBinV<-posBin[[TranscriptID]] ## Position of Bins for this TrancriptID
          
          ## Decided by now checking like a scanner through the transcript
          while(posP[2]<Tlength[[TranscriptID]]){
            
            posN1<-minpositive(posBinV-posP[1]) 
            posN2<-minpositive(posBinV-posP[2])
            
            pos1<-posN1+posP[1]
            pos2<-posN2+posP[2]
            
            r<-.Internal(which(posBinV==pos1)) ## The minimum positive will be the next bin pair
            Col<-.Internal(which(posBinV==pos2))
            
            #match? Faster? Try it
            
            
            if(posN1<posN2){
              
              posP<-posP+posN1+1 ## Updating the read positions
              v[posF[r]+Col-r+1]<-posN1[1]+1 ## Saving the obtained value for the actual bin pair
              
            }else{
              
              posP<-posP+posN2+1 ## Updating the read positions
              v[posF[r]+Col-r+1]<-posN2[1]+1 ## Saving the obtained value for the actual bin pair
              
            }
          }
          
          if(FragmentID==1){
            ## Obtain the bin combination names only for the first fragment length
            Index1 <-  rep(1:l,l:1)
            di <- diff(Index1)
            cs <- cumsum(di)
            Index2 <- 1:(l*(l+1)/2) + c(0,cumsum(di*cs))-c(0,cs)*l
            UpDiag <- paste0(TranscriptBins[Index1]," & ",TranscriptBins[Index2])
            TBins[[TranscriptID]]<-UpDiag
            names(v)<-TBins[[TranscriptID]]
            
            PairedM2[[TranscriptID]]<-v*Dens
            Allbins<-c(Allbins,UpDiag)
            Allbins<-unique(Allbins)
          }else{
            PairedM2[[TranscriptID]]<-PairedM2[[TranscriptID]]+v*Dens
          }
          
          
        }else{
          # MissedTranscripts[FlID]<-MissedTranscripts[FlID]+1
          # cat("missed ",TranscriptID)
          if(FragmentID==1){
            ## Obtain the bin combination names only for the first fragment length
            Index1 <-  rep(1:l,l:1)
            di <- diff(Index1)
            cs <- cumsum(di)
            Index2 <- 1:(l*(l+1)/2) + c(0,cumsum(di*cs))-c(0,cs)*l
            UpDiag <- paste0(TranscriptBins[Index1],"-",TranscriptBins[Index2])
            TBins[[TranscriptID]]<-UpDiag
            names(v)<-TBins[[TranscriptID]]
            
            Allbins<-c(Allbins,UpDiag)
            Allbins<-unique(Allbins)
          }
        }
      }
    }
    
    ## We only want the bin pairs once
    AllbinsU<-unique(Allbins)
    
    ## Initialization of the FinalMatrix
    FinalM<-matrix(0,nrow=length(AllbinsU),ncol=length(PairedM2))
    dim(FinalM)
    rownames(FinalM)<-AllbinsU
    
    ## Save the results in the Final M
    Rows<-match(names(unlist(PairedM2)),rownames(FinalM))
    Columns<-rep(1:Pathnum,lapply(PairedM2, length))
    FinalM[cbind(Rows,Columns)]<-unlist(PairedM2)
    # AllConvolved[[FlID]]<-FinalM
    
    ## Matrix Ranks
    # RankV[FlID]<-rankmatrix(round(FinalM),"chol")
    RankV<-rankmatrix(round(FinalM),"chol")
    # }
    
    Data<-list(FinalM=FinalM,possibleTranscripts=Pathnum,Rank=RankV)
    return(Data)
    
    ##### 3. Generation of combined fl libraries and returning data ########################################  
    # RankC<-combn(1:Flsize[2],2)
    # for(RankID in 1:dim(RankC)[2]){
    #   RankV[Flsize[2]+RankID]<-rankmatrix(round(rbind(AllConvolved[[RankC[1,RankID]]],AllConvolved[[RankC[2,RankID]]])),"chol")
    # }
    # ## Final list with all fl ranks
    # Data<-list(possibleTranscripts=Pathnum,Rank=RankV,Missed=MissedTranscripts,FinalM=FinalM)
    # return(Data)
  }else{
    return("to big to build the M matrix")
  }
  
}