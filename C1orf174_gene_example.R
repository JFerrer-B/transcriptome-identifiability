#Code to run the algorithm over the gene C1orf174

######### load required libraries #######
library(SGSeq)
library(igraph)
library(Matrix)
library(graph)
library(MASS)



########## load functions and internal functions ############

source("./functions/convertToSGFeatures2.R")
source("./functions/AuxFunctions.R")
source("./functions/getDistances.R")
source("./functions/plotNice.R")
source("./functions/GetBinGraph_single_gene.R")
source("./functions/Get_Fragment_Distribution.R")
source("./functions/GetConvolvedMatrixRank.R")




########## load gene gtf file #########
Transcripts <- importTranscripts(file = "./input_data/C1orf174_gene.gtf")
txf_ucsc <- convertToTxFeatures(Transcripts) 
sgf_ucsc <- convertToSGFeatures2(txf_ucsc) 

allgenes <- unique(unlist(geneName(sgf_ucsc)))
listgeneNames <- as.list(geneName(sgf_ucsc))

# Build a huge indicial matrix
A <-  Matrix(0,nrow = length(allgenes), 
             ncol = length(listgeneNames), sparse = T)
rownames(A) <- allgenes

A[cbind(match(unlist(listgeneNames),rownames(A)),rep(1:ncol(A), sapply(listgeneNames, length)))] <- 1


######### get the splicing graph ##########

# the splicing graph does not depend on the read and fragment length.

selectedLoci <- which(A[allgenes[1],]==1)
sgf <- sgf_ucsc[selectedLoci]
featureID(sgf) <- featureID(sgf)-min(featureID(sgf)) +1L

Distances <- getDistances(sgf)
Adjtxiki <- Distances$Adj
nPaths <- solve(diag(ncol(Adjtxiki)) - Adjtxiki)["S","E"]

# plot of the splicnig graph
C1orf174_splicing_graph <- plotNice(Adjtxiki)

jpeg(filename = "./output/C1orf174_gene_example/C1orf174_splicing_graph.jpeg")
plot(C1orf174_splicing_graph)
dev.off()

SplicingGraph_info <-list(Adjtxiki=Adjtxiki,sgf=sgf,nPaths=nPaths)


######  get BIN graph for each read length (75 and 114) ########

C1orf174_BinGraph_75nt <- GetBinGraph_single_gene(miSG = SplicingGraph_info,milr = 75)
C1orf174_BinGraph_114nt <- GetBinGraph_single_gene(miSG = SplicingGraph_info,milr = 114)

jpeg(filename = "./output/C1orf174_gene_example/C1orf174_BinGraph_75nt.jpeg")
plotNice(C1orf174_BinGraph_75nt$Adjacency)
dev.off()

jpeg(filename = "./output/C1orf174_gene_example/C1orf174_BinGraph_114nt.jpeg")
plotNice(C1orf174_BinGraph_114nt$Adjacency)
dev.off()

###### get RANK convolved Matrix ########

#get convolved matrix for read length = 75 and fragmengt length of 150

STD = 5 #select de standard deviation of the fragment length
num <- 11 #number of elements to approach the convolution
my_lr <- c(75) # the  read length
my_fl <- c(150) # the fragment length 
mis_lr_fl_STD_AMOUNT <- cbind(my_lr,my_fl,rep(STD,length(my_fl)),rep(num,length(my_fl)))
colnames(mis_lr_fl_STD_AMOUNT) <- c("read_length","fragment_length","STD","AMOUNT")
FlData <- Get_Fragment_Distribution(mis_lr_fl_STD_AMOUNT)

C1orf174_ConvMatrix_75_150 <- GetConvolvedMatrixRank(AA = C1orf174_BinGraph_75nt,FlData = FlData,lr = my_lr)
ncol(C1orf174_ConvMatrix_75_150$FinalM)
C1orf174_ConvMatrix_75_150$Rank


#get convolved matrix for read length = 114 and fragmengt length of 250

STD = 5 #select de standard deviation of the fragment length
num <- 11 #number of elements to approach the convolution
my_lr <- c(114) # the  read length
my_fl <- c(250) # the fragment length 
mis_lr_fl_STD_AMOUNT <- cbind(my_lr,my_fl,rep(STD,length(my_fl)),rep(num,length(my_fl)))
colnames(mis_lr_fl_STD_AMOUNT) <- c("read_length","fragment_length","STD","AMOUNT")
FlData <- Get_Fragment_Distribution(mis_lr_fl_STD_AMOUNT)

C1orf174_ConvMatrix_114_250 <- GetConvolvedMatrixRank(AA = C1orf174_BinGraph_114nt,FlData = FlData,lr = my_lr)
ncol(C1orf174_ConvMatrix_114_250$FinalM)
C1orf174_ConvMatrix_114_250$Rank
























































































































