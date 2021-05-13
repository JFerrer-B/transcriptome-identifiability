#gencode24_chr22

#this code apply the algorithm over the chr22 of the gencode24 transcriptome.

######### load required libraries #######
library(SGSeq)
library(igraph)
library(Matrix)
library(graph)
library(MASS)
library(pracma)
library(PoissonBinomial)
library(addreg)
library(doParallel)
library(foreach)
library(iterators)
library(numDeriv)
library(ggplot2)
library(hrbrthemes)
library(SparseM)
library(Rgraphviz)
library(IRanges)
library(S4Vectors)
library(turboEM)

########## load functions and internal functions ############

source("./functions/convertToSGFeatures2.R")
source("./functions/AuxFunctions.R")
source("./functions/getDistances.R")
source("./functions/plotNice.R")
source("./functions/GetBinGraph_single_gene.R")
source("./functions/Get_Fragment_Distribution.R")
source("./functions/GetConvolvedMatrixRank.R")
source("./functions/create_bin_graph.R")

#########  LOAD TRANSCRIPTOME DATA ##########


Transcripts <- importTranscripts(file = "input_data/genecode.24.chr22.gtf")
txf_ucsc <- convertToTxFeatures(Transcripts) 
sgf_ucsc <- convertToSGFeatures(txf_ucsc) 

allgenes <- unique(unlist(geneName(sgf_ucsc)))
listgeneNames <- as.list(geneName(sgf_ucsc))

# Build a huge indicial matrix
A <-  Matrix(0,nrow = length(allgenes), 
             ncol = length(listgeneNames), sparse = T)
rownames(A) <- allgenes

A[cbind(match(unlist(listgeneNames),rownames(A)),rep(1:ncol(A), sapply(listgeneNames, length)))] <- 1


########### built the splicing graph for each gene ##############
# as it doesn't depend on the lr value, we can built all the SG:


SG_list <- mclapply(allgenes, function(X){
  selectedLoci <- which(A[X,]==1)
  sgf <- sgf_ucsc[selectedLoci]
  featureID(sgf) <- featureID(sgf)-min(featureID(sgf)) +1L
  Distances <- getDistances(sgf)
  Adjtxiki <- Distances$Adj
  nPaths <- solve(diag(ncol(Adjtxiki)) - Adjtxiki)["S","E"]
  return(list(Adjtxiki=Adjtxiki,sgf=sgf,nPaths=nPaths))
},mc.cores = 16)

names(SG_list) <- allgenes
SG_list <- Filter(Negate(is.null), SG_list)

save(SG_list,file="./output/gencode24_chr22/encode24_SplicingGraph_list.RData")



###### build the bing graph for each gene and for each lr #######


### we used 75, 100, 150, 200 and 300 read length
# load("./output/gencode24_chr22/encode24_SplicingGraph_list.RData")
all_posible_lr <- c(75,100,150,200,300)
# nnx <- 1

for (nnx in 1:length(all_posible_lr)) {
  milr <- all_posible_lr[nnx]
  command <- paste0("Bin_Graph_",milr," <- create_bin_graph(SG_list,milr,ncores = 16)")
  eval(parse(text = command))
  
  command_2 <- paste0('save(Bin_Graph_',milr,',file = "./output/gencode24_chr22/Bin_Graph_',milr,'.RData")')
  eval(parse(text = command_2))
  
}


########### FOR EACH LR (AND SEVERAL FLDATA) WE BUILD THE Convolved MATRIX FOR EACH GENE ###############
# load("./output/gencode24_chr22/Bin_Graph_75.RData")
# load("./output/gencode24_chr22/Bin_Graph_100.RData")
# load("./output/gencode24_chr22/Bin_Graph_200.RData")
# load("./output/gencode24_chr22/Bin_Graph_150.RData")
# load("./output/gencode24_chr22/Bin_Graph_300.RData")

STD = 50
num <- 11
los_lr <- c(75,75,rep(100,6),150,150,300,300,300,200,200,200)
los_fl <- c(150,400,100,200,300,400,600,800,300,600,600,900,500,350,400,500)

mis_lr_fl_STD_AMOUNT <- cbind(los_lr,los_fl,rep(STD,length(los_fl)),rep(num,length(los_fl)))
colnames(mis_lr_fl_STD_AMOUNT) <- c("read_length","fragment_length","STD","AMOUNT")
lr <- mis_lr_fl_STD_AMOUNT[,1]
FlData <- Get_Fragment_Distribution(mis_lr_fl_STD_AMOUNT)


all_posible_lr <- unique(lr)
for(nnx in 1:length(all_posible_lr)){
  milr <- all_posible_lr[nnx]
  indexfld <- which(as.numeric(gsub("-.*","",colnames(FlData$fl)))==milr)
  for(jjx in 1:length(indexfld)){
    FlData_2 <- FlData
    FlData_2$fl <- FlData_2$fl[,indexfld[jjx],drop=F]
    FlData_2$weight <- FlData_2$weight[,indexfld[jjx],drop=F]
    
    possiblefile <- paste0('./output/gencode24_chr22/ConvolvedMatrix_',gsub("-","_",colnames(FlData_2$fl)),'.RData')
    
    if(file.exists(possiblefile)) next
    
    
    command_3 <- paste0("ConvolvedMatrix_",gsub("-","_",colnames(FlData_2$fl))," <- mclapply(Bin_Graph_",milr,",function(X){
      GetConvolvedMatrixRank(AA = X,FlData = FlData_2,lr = milr)
    },mc.cores = 16)")
    cat(command_3,sep = "\n")
    eval(parse(text = command_3))
    
        command_4 <- paste0('save(ConvolvedMatrix_',gsub("-","_",colnames(FlData_2$fl)),',file="./output/gencode24_chr22/ConvolvedMatrix_',gsub("-","_",colnames(FlData_2$fl)),'.RData")')
    cat(command_4,sep = "\n")
    eval(parse(text = command_4))
    
    
  }
  
}



######### summary of the results #############
## function to get the % of identifiable genes

data_summary(ConvolvedMatrix_75_150_50_11)
data_summary(ConvolvedMatrix_75_400_50_11)
data_summary(ConvolvedMatrix_100_100_50_11)
data_summary(ConvolvedMatrix_100_200_50_11)
data_summary(ConvolvedMatrix_100_300_50_11)
data_summary(ConvolvedMatrix_100_400_50_11)
data_summary(ConvolvedMatrix_100_600_50_11)
data_summary(ConvolvedMatrix_100_800_50_11)
data_summary(ConvolvedMatrix_150_300_50_11)
data_summary(ConvolvedMatrix_150_600_50_11)
data_summary(ConvolvedMatrix_200_350_50_11)
data_summary(ConvolvedMatrix_200_400_50_11)
data_summary(ConvolvedMatrix_200_500_50_11)
data_summary(ConvolvedMatrix_300_500_50_11)
data_summary(ConvolvedMatrix_300_600_50_11)
data_summary(ConvolvedMatrix_300_900_50_11)


######### merge different fragment length #######

Conv_100_200_800 <- lapply(1:length(ConvolvedMatrix_100_200_50_11),function(X){
  if(class(ConvolvedMatrix_100_200_50_11[[X]])=="list"){
    CombM <- rbind(ConvolvedMatrix_100_200_50_11[[X]]$FinalM,ConvolvedMatrix_100_800_50_11[[X]]$FinalM)
    therank_comb <- rankmatrix(CombM)
    return(list(FinalM=CombM,possibleTranscripts=ConvolvedMatrix_100_200_50_11[[X]]$possibleTranscripts,Rank=therank_comb))
  }
})

data_summary(ConvolvedMatrix_100_200_50_11)
data_summary(ConvolvedMatrix_100_800_50_11)
data_summary(Conv_100_200_800)
















