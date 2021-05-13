plotNice <- function(Adj,SS=5,h=1.5,w=3) {
  library(Rgraphviz)
  Adj2 <- Adj

  colnames(Adj2) <- rownames(Adj2) <- gsub("[|]",":",colnames(Adj))
  startnodes  <- which(Adj2["S",]>0)
  endnodes <- which(Adj2[,"E"]>0)
  fakes <- match(c("S","E"), colnames(Adj2))
  Adj3 <- Adj2[-fakes, -fakes, drop = FALSE] # Remove fake nodes
  Salida  <- as(Adj3, "graphNEL") # Adj3: without start and end nodes

  coloring <- rep("light grey",ncol(Adj3))
  names(coloring) <- colnames(Adj3)
  coloring[names(startnodes)] <- "green"
  coloring[names(endnodes)] <- "red"

  nAttrs  <- makeNodeAttrs(Salida, fillcolor = coloring[match(nodes(Salida), names(coloring))], 
                             shape = "box",height=h,width=w,
                             fontsize=max(nchar(colnames(Adj3)))^.6+SS,
                             fixedsize=FALSE)
  
  plot(Salida,nodeAttrs = nAttrs)
  
  # TODO: Put a different color if a node is a start node and an end node.
}