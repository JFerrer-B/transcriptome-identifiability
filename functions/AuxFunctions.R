#AuxFunctions

# Axuliar function called by the main functions.

processFeatures2 <- function(features, coerce = FALSE, 
                             merge = FALSE) {
  junctions <- granges(features)[type(features) == 
                                   "J"]
  junctions_D <- flank(junctions, -1, TRUE)
  junctions_A <- flank(junctions, -1, FALSE)
  mcols(junctions)$type <- rep("J", length(junctions))
  if (is(features, "TxFeatures")) {
    exons <- features[type(features) %in% 
                        c("I", "F", "L", "U")]
    exons_D <- flank(features[type(features) %in% 
                                c("I", "F")], -1, FALSE)
    exons_A <- flank(features[type(features) %in% 
                                c("I", "L")], -1, TRUE)
  } else if (is(features, "SGFeatures")) {
    exons <- features[type(features) == 
                        "E"]
    exons_D <- flank(features[splice3p(features)], 
                     -1, FALSE)
    exons_A <- flank(features[splice5p(features)], 
                     -1, TRUE)
  }
  exons <- granges(exons)
  exons_D <- granges(exons_D)
  exons_A <- granges(exons_A)
  D <- unique(c(junctions_D, exons_D))
  mcols(D)$type <- rep("D", length(D))
  A <- unique(c(junctions_A, exons_A))
  mcols(A)$type <- rep("A", length(A))
  splicesites <- c(D, A)
  other <- c(junctions, splicesites)
  exons <- disjoin(exons)
  exons_start <- flank(exons, -1, TRUE)
  exons_end <- flank(exons, -1, FALSE)
  i_q <- which(!exons_end %over% splicesites)
  i_s <- which(!exons_start %over% splicesites)
  ol <- findOverlaps(suppressWarnings(flank(exons[i_q], 
                                            1, FALSE)), exons_start[i_s])
  if ((length(ol) > 0)) {
    qH <- i_q[queryHits(ol)]
    sH <- i_s[subjectHits(ol)]
    i_to_be_merged <- union(qH, sH)
    d <- data.frame(from = qH, to = sH)
    v <- data.frame(name = i_to_be_merged)
    g <- graph.data.frame(d = d, directed = TRUE, 
                          vertices = v)
    k <- clusters(g)$membership
    exons_to_be_merged <- split(exons[i_to_be_merged], 
                                k)
    exons_merged <- unlist(reduce(exons_to_be_merged))
    if (length(exons_to_be_merged) != 
        length(exons_merged)) {
      stop("cannot merge non-adjacent exons")
    }
    if (merge) {
      exons <- c(exons[-i_to_be_merged], 
                 exons_merged)
    }
  }
  exons_start <- flank(exons, -1, TRUE)
  exons_end <- flank(exons, -1, FALSE)
  splice5p <- rep(FALSE, length(exons))
  i_spliced <- unique(queryHits(findOverlaps(exons_start, 
                                             A)))
  i_adjacent <- unique(queryHits(findOverlaps(suppressWarnings(flank(exons, 
                                                                     1, TRUE)), exons)))
  splice5p[setdiff(i_spliced, i_adjacent)] <- TRUE
  splice3p <- rep(FALSE, length(exons))
  i_spliced <- unique(queryHits(findOverlaps(exons_end, 
                                             D)))
  i_adjacent <- unique(queryHits(findOverlaps(suppressWarnings(flank(exons, 
                                                                     1, FALSE)), exons)))
  splice3p[setdiff(i_spliced, i_adjacent)] <- TRUE
  mcols(exons)$type <- rep("E", length(exons))
  mcols(exons)$splice5p <- splice5p
  mcols(exons)$splice3p <- splice3p
  mcols(other)$splice5p <- rep(NA, length(other))
  mcols(other)$splice3p <- rep(NA, length(other))
  features <- setNames(c(exons, other), 
                       NULL)
  features <- sort(features)
  return(features)
}

getAdj <- function(sgf) {
  nodes <- sgf[type(sgf)=="E"]
  # edges <- sgf[type(sgf)!="E"]
  edges <- sgf[type(sgf)=="J"]
  
  nodenames <- c("S",featureID(nodes),"E")
  Adj <-  matrix(0, ncol = length(nodenames), nrow = length(nodenames))
  colnames(Adj) <-  rownames(Adj) <- nodenames
  
  # Join subexons
  Adj[cbind(match(start(nodes), end(nodes)+1)[-1]+1, match(end(nodes)+1, start(nodes))[-length(nodes)]+1)] <- 1
  # Join junctions
  Adj[cbind(match(start(edges), end(nodes))+1, match(end(edges), start(nodes))+1)] <- 1
  # Include Start nodes
  Trans <- unique(unlist(txName(nodes)))
  for (n in 1:length(Trans)) {
    Adj["S",1+min(which(unlist(lapply(match(txName(nodes), Trans[n]), any))))] <- 1
    Adj[1+max(which(unlist(lapply(match(txName(nodes), Trans[n]), any)))),"E"] <- 1
  }
  return(Adj)
}


getAdjBinGraph <- function(sgf, lr = 200) {
  # Transcripts
  Transcripts <- unique(unlist(as.list(txName(sgf))))
  Adjacency <- NULL
  Bindata <- NULL
  TranscriptBins <- vector(mode = "list",length = length(Transcripts))
  names(TranscriptBins) <- Transcripts
  
  for (TranscriptID in 1:length(Transcripts)) {
    # print(TranscriptID)
    # if(TranscriptID == 8) browser()
    subexons <- intersect(grep(Transcripts[TranscriptID], as.list(txName(sgf))), which(type(sgf)=="E"))
    sgftxiki <- sgf[subexons,]
    endtrans <- cumsum(width(sgftxiki))
    exonlengths <-  width(sgftxiki)
    begintrans <- c(1,endtrans[-length(endtrans)])
    
    A1 <- outer(endtrans, begintrans, FUN ="-")
    B1 <- outer(begintrans, endtrans, FUN ="-")
    A <- A1>lr
    B <- B1<lr
    if (!all(dim(A)==dim(B))) {
      # warning("Null transcript!!??")
      # Transcripts <- Transcripts[-TranscriptID]
      next
    }
    Aux <- A*B
    colnames(Aux) <- rownames(Aux) <- featureID(sgftxiki)
    
    Subexons <- which(Aux>0, arr.ind=TRUE)
    Bins <- data.frame(names = rep("",nrow(Subexons)), sizes = rep(NA,nrow(Subexons)),stringsAsFactors = FALSE)
    for (n in 1:nrow(Subexons)) {
      if(nrow(Subexons)==0) {
        # warning("Read length larger than some transcripts!")
        break;
      }
      Bins[n,1] <- paste(rownames(Aux)[Subexons[n,2]:Subexons[n,1]], collapse="|",sep="|")
      # Include the lengths of the SubexonBins
      # If only one: le -2(lr-1); later on correct for start and end nodes
      # if several exons: 
      #   two exons: 
      #   More exons: 2*(lr-1) + sum(lengths_middle_exons)
      if (Subexons[n,2] == Subexons[n,1]) { # One exon
        Bins[n,2] <- exonlengths[Subexons[n,1]]+1-lr
      } else {
        lm <- ifelse((Subexons[n,1] - Subexons[n,2] ==1), 
                     0, sum(exonlengths[(1+Subexons[n,2]):(Subexons[n,1]-1)]))
        lt <- sum(exonlengths[Subexons[n,2]:Subexons[n,1]])
        Bins[n,2] <- min(lt - lr, exonlengths[Subexons[n,2]] -1) -
          max(0,exonlengths[Subexons[n,2]]+lm-(lr-1))+1
      }
    }
    if(is.null(Adjacency)) {
      Adjacency <- matrix(0, ncol = 2 + length(Bins$names), nrow = 2+length(Bins$names))
      Bindata <- Bins
      colnames(Adjacency) <- rownames(Adjacency) <- c("S",Bins$names, "E")
    }
    TranscriptBins[[TranscriptID]] <- Bins$names
    moreBins <- setdiff(Bins$names, colnames(Adjacency))
    if(length(moreBins)>0) {
      Bindata <- unique(rbind(Bins,Bindata))
      newnames <- c(colnames(Adjacency),moreBins)
      Adjacency <- cbind(Adjacency,matrix(0, ncol= length(moreBins), nrow=nrow(Adjacency)))
      Adjacency <- rbind(Adjacency,matrix(0, ncol = ncol(Adjacency), nrow=length(moreBins)))
      colnames(Adjacency) <- rownames(Adjacency) <- newnames
    }
    Indices <- match(c("S",Bins$names,"E"), colnames(Adjacency))
    Aux <- Adjacency
    Aux[cbind(Indices[1:(length(Indices)-1)], Indices[2:length(Indices)])] <- 1
    Adjacency <- Adjacency + Aux
  }
  Adjacency[Adjacency>1] <- 1
  return(list(Adjacency=Adjacency, Bindata = Bindata, 
              TranscriptBins=TranscriptBins, Transcripts=Transcripts))
}


contractGraph <- function(adjbin){
  
  Adj <- adjbin$Adjacency
  
  Inc <- matrix(0, nrow = ncol(Adj), ncol = sum(Adj > 0))
  Inc[cbind(which(Adj > 0, arr.ind = 1)[, 2], seq_len(ncol(Inc)))] <- 1
  Inc[cbind(which(Adj > 0, arr.ind = 1)[, 1], seq_len(ncol(Inc)))] <- -1
  
  rownames(Inc) <- rownames(Adj)
  colnames(Inc) <- seq_len(ncol(Inc))
  
  mis_inc <-names(which(rowSums(Inc==1)>1))
  
  Adj_2 <- Adj
  Adj_2["S",mis_inc] <- 1
  
  Inc_2 <- matrix(0, nrow = ncol(Adj_2), ncol = sum(Adj_2 > 0))
  Inc_2[cbind(which(Adj_2 > 0, arr.ind = 1)[, 2], seq_len(ncol(Inc_2)))] <- 1
  Inc_2[cbind(which(Adj_2 > 0, arr.ind = 1)[, 1], seq_len(ncol(Inc_2)))] <- -1
  
  rownames(Inc_2) <- rownames(Adj_2)
  colnames(Inc_2) <- seq_len(ncol(Inc_2))
  
  randSol <- getRandomFlow(Inc_2,10)
  randSol <- (Inc_2==1) %*% randSol
  X <- as.matrix(dist(randSol))
  Inc_2 <- (X < 1e-08)
  g <- graph_from_adjacency_matrix(Inc_2)
  Groups <- clusters(g)
  misgrupos <- Groups$membership
  
  startnode <- names(which(Adj["S",]==1))
  endnode <- names(which(Adj[,"E"]==1))
  misgrupos <- misgrupos[-which(names(misgrupos) %in% c(startnode,endnode,"S","E"))]
  
  finalgruups <- lapply(unique(misgrupos),function(X){
    nnx <- which(misgrupos==X)
    if(length(nnx)==1){
      NULL
    }else{
      names(nnx)
    }
  })
  
  finalgruups <- Filter(Negate(is.null), finalgruups)
  newbins <- sapply(finalgruups,function(X) paste(X,collapse = "-"))
  names(finalgruups) <- newbins
  
  
  if(length(finalgruups)>0){
    Adj_2 <- adjbin$Adjacency
    Bindata_2 <- adjbin$Bindata
    TranscriptBins_2 <- adjbin$TranscriptBins
    Transcripts_2 <- adjbin$Transcripts
    
    for(jjx in 1:length(finalgruups)){
      # jjx <- 1
      rrx <- match(finalgruups[[jjx]],rownames(Adj_2))
      newline <- colSums2(Adj_2[rrx,])
      Adj_2 <- Adj_2[-rrx,]
      Adj_2 <- rbind(Adj_2,newline)
      rownames(Adj_2)[nrow(Adj_2)] <- names(finalgruups)[jjx]
      
      rrx <- match(finalgruups[[jjx]],colnames(Adj_2))
      newline <- rowSums2(Adj_2[,rrx])
      Adj_2 <- Adj_2[,-rrx]
      Adj_2 <- cbind(Adj_2,newline)
      colnames(Adj_2)[ncol(Adj_2)] <- names(finalgruups)[jjx]
      
      rrx <- match(finalgruups[[jjx]],Bindata_2$names)
      toad <- data.frame(names=names(finalgruups)[jjx],sizes=sum(Bindata_2$sizes[rrx]),stringsAsFactors = F)
      Bindata_2 <- Bindata_2[-rrx,]
      Bindata_2 <- rbind(Bindata_2,toad)
      
      TranscriptBins_2 <- lapply(TranscriptBins_2,function(X){
        rrx <- match(finalgruups[[jjx]],X)
        rrx <- rrx[!is.na(rrx)]
        if(length(rrx)>1){
          X[min(rrx)] <- names(finalgruups)[jjx]
          X <- X[-max(rrx)]
          X  
        }else{
          X
        }
        
      })
      
      
    }
    diag(Adj_2) <- 0
    Adj_2
    Bindata_2
    TranscriptBins_2
    
    adjbin2 <- adjbin
    adjbin2$Adjacency <- Adj_2
    adjbin2$Bindata <- Bindata_2
    adjbin2$TranscriptBins <- TranscriptBins_2
    return(adjbin2)
  }else{
    return(adjbin)
  }
  
  
  
}

getRandomFlow <- function(Incidence, ncol = 1) {
  # With the incidence matrix, it is
  # possible to get its null-space and
  # generate an arbitrary flow on it. Using
  # the flow it is possible to get the
  # triplets of events.
  
  # The seed is set to ensure the order of
  # events remains the same
  # set.seed("0xABBA")
  
  # Solve the Null Space for the Incidence
  # Matrix
  solh <- Null(t(Incidence))
  
  # Condition to ensure that everything
  # that exits the Start node (-1), exits
  # at the End Node (1)
  myvec <- vector(mode = "numeric",length = nrow(Incidence))
  myvec[which(rownames(Incidence)=="S")] <- -1
  myvec[which(rownames(Incidence)=="E")] <- 1
  solp <- ginv(Incidence) %*% myvec
  
  # Matrix of fluxes, with as many columns
  # as specified by the user
  v <- matrix(runif(ncol(solh) * ncol), 
              ncol = ncol)
  randSol <- as.vector(solp) + solh %*% 
    v
  
  return(randSol)
}


minpositive<- function(x) {
  # Small function to calculate the minimum positive value on a vector
  min(x[x > 0])
}

rankmatrix <- function(A, method = "chol") {
  ## Calculation of the rank of the matrix (crossproduct with transposed for faster computation)
  dims <- dim(A)
  if (dims[1] > dims[2]) {
    A <- crossprod(A)
  } else {A <- tcrossprod(A)}
  if (method == "chol") {
    
    rango <- suppressWarnings(attr(chol(A, pivot = T), 'rank'))
    
    
  } else { rango <- qr(A)$rank }
  return(rango)
}



data_summary <- function(ConvMatrix){
  
  sum_index <- t(sapply(ConvMatrix,function(X){
    if(class(X)=="list"){
      c(X$possibleTranscripts,X$Rank)
    }else{
      c(0,0)
    }
  }))
  
  isid <- sum_index[,2]/sum_index[,1]
  a <- table(isid==1)
  return(a[2]/(a[1]+a[2])*100)
}

