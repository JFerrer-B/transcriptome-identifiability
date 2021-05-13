#Code to run the algorithm over the a simulated gene to test the advantages of using two libraries

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


######### load gene data ########
Transcripts <- importTranscripts(file = "./input_data/simulated_gene.gtf")
txf_ucsc <- convertToTxFeatures(Transcripts) 
sgf_ucsc <- convertToSGFeatures(txf_ucsc) 
sgf_ucsc@txName[10][[1]] <- sgf_ucsc@txName[1][[1]]
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
jpeg(filename = "./output/twolibraries_simulation/SimulatedGene_splicing_graph.jpeg")
plotNice(Adjtxiki)
dev.off()

miSG <-list(Adjtxiki=Adjtxiki,sgf=sgf,nPaths=nPaths)

######  get BIN graph for read length of 100 nt ########

Simulated_BinGraph_100nt <- GetBinGraph_single_gene(miSG = miSG,milr = 100)

jpeg(filename = "./output/twolibraries_simulation/Simulated_BinGraph_100nt.jpeg")
plotNice(Simulated_BinGraph_100nt$Adjacency)
dev.off()


####### get Convolved Matrix for fragment length of 200 and fragment length of 800 #######
STD = 5 #select de standard deviation of the fragment length
num <- 11 #number of elements to approach the convolution
my_lr <- c(100) # the  read length
my_fl <- c(200) # the fragment length 
mis_lr_fl_STD_AMOUNT <- cbind(my_lr,my_fl,rep(STD,length(my_fl)),rep(num,length(my_fl)))
colnames(mis_lr_fl_STD_AMOUNT) <- c("read_length","fragment_length","STD","AMOUNT")
FlData <- Get_Fragment_Distribution(mis_lr_fl_STD_AMOUNT)

Mmatrices_200 <- GetConvolvedMatrixRank(AA = Simulated_BinGraph_100nt,FlData = FlData,lr = my_lr)


STD = 5 #select de standard deviation of the fragment length
num <- 11 #number of elements to approach the convolution
my_lr <- c(100) # the  read length
my_fl <- c(800) # the fragment length 
mis_lr_fl_STD_AMOUNT <- cbind(my_lr,my_fl,rep(STD,length(my_fl)),rep(num,length(my_fl)))
colnames(mis_lr_fl_STD_AMOUNT) <- c("read_length","fragment_length","STD","AMOUNT")
FlData <- Get_Fragment_Distribution(mis_lr_fl_STD_AMOUNT)

Mmatrices_800 <- GetConvolvedMatrixRank(AA = Simulated_BinGraph_100nt,FlData = FlData,lr = my_lr)


######### bootstrap simulacion ########

Mmatrices_200$FinalM <- Mmatrices_200$FinalM + 1e-6
Mmatrices_800$FinalM <- Mmatrices_800$FinalM + 1e-6

expresion <- c(1,1,1,1,0,0)
y1 <- Mmatrices_200$FinalM %*% expresion
y1s <- rpois(length(y1),y1)
names(y1s) <- 1:length(y1s)

y2 <- 0.1*Mmatrices_800$FinalM %*% expresion
y2s <- rpois(length(y2),y2)
names(y2s) <- 1:length(y2s)


boots_fl_200 <- rbind()
boots_fl_800 <- rbind()
boots_fl_both <- rbind()

for(n in 1:500){
  y1sbb <- y1s * 0
  A <- rep(names(y1s),y1s)
  dummy <- table(sample(A,replace = TRUE,size = sum(y1s)))
  y1sbb[as.numeric(names(dummy))] <- dummy 
  
  # y1s
  # y1sbb
  # sum(y1s)
  # sum(y1sbb)
  
  y2sbb <- y2s * 0
  A <- rep(names(y2s),y2s)
  dummy <- table(sample(A,replace = TRUE,size = sum(y2s)))
  y2sbb[as.numeric(names(dummy))] <- dummy 
  
  # y2s
  # y2sbb
  # sum(y2s)
  # sum(y2sbb)
  
  
  expest1 <- coef(nnpois(x=Mmatrices_200$FinalM,y = y1sbb, standard = rep(1,length(y1sbb)), offset = 0, start = runif(length(expresion))))
  expest2 <- coef(nnpois(x=Mmatrices_800$FinalM,y = y2sbb, standard = rep(1,length(y2sbb)), offset = 0, start = runif(length(expresion))))
  
  lambda <- sum(y2sbb)/sum(y1sbb)
  
  for(k in 1:100){
    M2 <- Mmatrices_800$FinalM*lambda
    A <- rbind(rbind(Mmatrices_200$FinalM,M2))
    # nnls::nnls(A=A,b = c(y1s,y2s))
    expest12 <- coef(nnpois(x=A,y = c(y1sbb,y2sbb), standard = rep(1,length(c(y1s,y2s))), offset = 0, start = runif(length(expresion))))
    lambda_1 <- sum(y2sbb)/sum(Mmatrices_800$FinalM %*% expest12)
    mydiff <- abs(lambda_1-lambda)
    if(mydiff < 1e-6){
      break
    }else{
      lambda <- lambda_1
    }
  }
  
  boots_fl_200 <- rbind(boots_fl_200,expest1)
  boots_fl_800 <- rbind(boots_fl_800,expest2)
  boots_fl_both <- rbind(boots_fl_both,expest12)
}




####### tran1 vs tran 2
Transcript_1 <- c(boots_fl_200[,1],boots_fl_800[,1],boots_fl_both[,1])
Transcript_2 <- c(boots_fl_200[,2],boots_fl_800[,2],boots_fl_both[,2])
fragmenet_length <- rep(c("200nt","800nt","200|800nt"),c(nrow(boots_fl_200),nrow(boots_fl_800),nrow(boots_fl_both)))

df <- data.frame(Transcript_1 = Transcript_1,Transcript_2 = Transcript_2,fragmen_length = fragmenet_length)
df <- rbind(df,data.frame(Transcript_1=1,Transcript_2=1,fragmen_length = "Expected Value"))
head(df)
colnames(df) <- c("Transcript 1","Transcript 2","Fragment Length")
df[nrow(df),]

df$`Fragment Length`
title_axis <- 14
text_size <- 12
mishape <- rep(16,nrow(df))
mishape[length(mishape)] <- 23
misize <- rep(1,nrow(df))
misize[length(misize)] <- 5

# nn <- 0
# nn <- nn+1
pp <- ggplot(df, aes(x=`Transcript 1`, y=`Transcript 2`, group=`Fragment Length`)) + 
  geom_point(aes(shape=`Fragment Length`,color=`Fragment Length`,size=`Fragment Length`,fill=`Fragment Length`)) + 
  scale_size_manual(values = c(2,2,2,4)) +
  scale_shape_manual(values = c(16,16,16,24)) +
  scale_color_manual(values = c("#00798c","#d1495b","#11773390","#000000")) +
  scale_fill_manual(values = c("#00798c","#d1495b","#7CAE0050","white")) +
  theme_bw() + ylim(0,1.25) + xlim(0,2.5)
pp

pp <- pp + theme(axis.title.x =element_text(size=title_axis),axis.title.y =element_text(size=title_axis),
           axis.text = element_text(size=text_size),legend.text = element_text(size = text_size),legend.title = element_blank())
pp

jpeg(filename = "./output/twolibraries_simulation/tran1_vs_tran2_simulacion.jpeg")
pp
dev.off()



####### tran1 vs tran 3
Transcript_1 <- c(boots_fl_200[,1],boots_fl_800[,1],boots_fl_both[,1])
Transcript_2 <- c(boots_fl_200[,3],boots_fl_800[,3],boots_fl_both[,3])
fragmenet_length <- rep(c("200nt","800nt","200|800nt"),c(nrow(boots_fl_200),nrow(boots_fl_800),nrow(boots_fl_both)))

df <- data.frame(Transcript_1 = Transcript_1,Transcript_2 = Transcript_2,fragmen_length = fragmenet_length)
df <- rbind(df,data.frame(Transcript_1=1,Transcript_2=1,fragmen_length = "Expected Value"))
head(df)
colnames(df) <- c("Transcript 1","Transcript 3","Fragment Length")
df[nrow(df),]


df$`Fragment Length`
title_axis <- 14
text_size <- 12
mishape <- rep(16,nrow(df))
mishape[length(mishape)] <- 23
misize <- rep(1,nrow(df))
misize[length(misize)] <- 5

pp <- ggplot(df, aes(x=`Transcript 1`, y=`Transcript 3`, group=`Fragment Length`)) + 
  geom_point(aes(shape=`Fragment Length`,color=`Fragment Length`,size=`Fragment Length`,fill=`Fragment Length`)) + 
  scale_size_manual(values = c(2,2,2,4)) +
  scale_shape_manual(values = c(16,16,16,24)) +
  scale_color_manual(values = c("#00798c","#d1495b","#11773390","#000000")) +
  scale_fill_manual(values = c("#00798c","#d1495b","#7CAE0050","white")) +
  theme_bw() + ylim(0,1.25) + xlim(0,2.5)
pp

pp <- pp + theme(axis.title.x =element_text(size=title_axis),axis.title.y =element_text(size=title_axis),
                 axis.text = element_text(size=text_size),legend.text = element_text(size = text_size),legend.title = element_blank())
pp

jpeg(filename = "./output/twolibraries_simulation/tran1_vs_tran3_simulacion.jpeg")
pp
dev.off()

####### tran1 vs tran 4
Transcript_1 <- c(boots_fl_200[,1],boots_fl_800[,1],boots_fl_both[,1])
Transcript_2 <- c(boots_fl_200[,4],boots_fl_800[,4],boots_fl_both[,4])
fragmenet_length <- rep(c("200nt","800nt","200|800nt"),c(nrow(boots_fl_200),nrow(boots_fl_800),nrow(boots_fl_both)))

df <- data.frame(Transcript_1 = Transcript_1,Transcript_2 = Transcript_2,fragmen_length = fragmenet_length)
df <- rbind(df,data.frame(Transcript_1=1,Transcript_2=1,fragmen_length = "Expected Value"))
head(df)
colnames(df) <- c("Transcript 1","Transcript 4","Fragment Length")
df[nrow(df),]


df$`Fragment Length`
title_axis <- 14
text_size <- 12
mishape <- rep(16,nrow(df))
mishape[length(mishape)] <- 23
misize <- rep(1,nrow(df))
misize[length(misize)] <- 5

pp <- ggplot(df, aes(x=`Transcript 1`, y=`Transcript 4`, group=`Fragment Length`)) + 
  geom_point(aes(shape=`Fragment Length`,color=`Fragment Length`,size=`Fragment Length`,fill=`Fragment Length`)) + 
  scale_size_manual(values = c(2,2,2,4)) +
  scale_shape_manual(values = c(16,16,16,24)) +
  scale_color_manual(values = c("#00798c","#d1495b","#11773390","#000000")) +
  scale_fill_manual(values = c("#00798c","#d1495b","#7CAE0050","white")) +
  theme_bw() + ylim(0,1.25) + xlim(0,2.5)
pp

pp <- pp + theme(axis.title.x =element_text(size=title_axis),axis.title.y =element_text(size=title_axis),
                 axis.text = element_text(size=text_size),legend.text = element_text(size = text_size),legend.title = element_blank())
pp

jpeg(filename = "./output/twolibraries_simulation/tran1_vs_tran4_simulacion.jpeg")
pp
dev.off()


####### tran2 vs tran 4
Transcript_1 <- c(boots_fl_200[,2],boots_fl_800[,2],boots_fl_both[,2])
Transcript_2 <- c(boots_fl_200[,4],boots_fl_800[,4],boots_fl_both[,4])
fragmenet_length <- rep(c("200nt","800nt","200|800nt"),c(nrow(boots_fl_200),nrow(boots_fl_800),nrow(boots_fl_both)))

df <- data.frame(Transcript_1 = Transcript_1,Transcript_2 = Transcript_2,fragmen_length = fragmenet_length)
df <- rbind(df,data.frame(Transcript_1=1,Transcript_2=1,fragmen_length = "Expected Value"))
head(df)
colnames(df) <- c("Transcript 2","Transcript 4","Fragment Length")
df[nrow(df),]


df$`Fragment Length`
title_axis <- 14
text_size <- 12
mishape <- rep(16,nrow(df))
mishape[length(mishape)] <- 23
misize <- rep(1,nrow(df))
misize[length(misize)] <- 5

pp <- ggplot(df, aes(x=`Transcript 2`, y=`Transcript 4`, group=`Fragment Length`)) + 
  geom_point(aes(shape=`Fragment Length`,color=`Fragment Length`,size=`Fragment Length`,fill=`Fragment Length`)) + 
  scale_size_manual(values = c(2,2,2,4)) +
  scale_shape_manual(values = c(16,16,16,24)) +
  scale_color_manual(values = c("#00798c","#d1495b","#11773390","#000000")) +
  scale_fill_manual(values = c("#00798c","#d1495b","#7CAE0050","white")) +
  theme_bw() + ylim(0,1.25) + xlim(0,1.5)
pp

pp <- pp + theme(axis.title.x =element_text(size=title_axis),axis.title.y =element_text(size=title_axis),
                 axis.text = element_text(size=text_size),legend.text = element_text(size = text_size),legend.title = element_blank())
pp

jpeg(filename = "./output/twolibraries_simulation/tran2_vs_tran4_simulacion.jpeg")
pp
dev.off()


####### tran3 vs tran 4
Transcript_1 <- c(boots_fl_200[,3],boots_fl_800[,3],boots_fl_both[,3])
Transcript_2 <- c(boots_fl_200[,4],boots_fl_800[,4],boots_fl_both[,4])
fragmenet_length <- rep(c("200nt","800nt","200|800nt"),c(nrow(boots_fl_200),nrow(boots_fl_800),nrow(boots_fl_both)))

df <- data.frame(Transcript_1 = Transcript_1,Transcript_2 = Transcript_2,fragmen_length = fragmenet_length)
df <- rbind(df,data.frame(Transcript_1=1,Transcript_2=1,fragmen_length = "Expected Value"))
head(df)
colnames(df) <- c("Transcript 3","Transcript 4","Fragment Length")
df[nrow(df),]


df$`Fragment Length`
title_axis <- 14
text_size <- 12
mishape <- rep(16,nrow(df))
mishape[length(mishape)] <- 23
misize <- rep(1,nrow(df))
misize[length(misize)] <- 5

pp <- ggplot(df, aes(x=`Transcript 3`, y=`Transcript 4`, group=`Fragment Length`)) + 
  geom_point(aes(shape=`Fragment Length`,color=`Fragment Length`,size=`Fragment Length`,fill=`Fragment Length`)) + 
  scale_size_manual(values = c(2,2,2,4)) +
  scale_shape_manual(values = c(16,16,16,24)) +
  scale_color_manual(values = c("#00798c","#d1495b","#11773390","#000000")) +
  scale_fill_manual(values = c("#00798c","#d1495b","#7CAE0050","white")) +
  theme_bw() + ylim(0,1.25) + xlim(0,1.5)
pp

pp <- pp + theme(axis.title.x =element_text(size=title_axis),axis.title.y =element_text(size=title_axis),
                 axis.text = element_text(size=text_size),legend.text = element_text(size = text_size),legend.title = element_blank())
pp

jpeg(filename = "./output/twolibraries_simulation/tran3_vs_tran4_simulacion.jpeg")
pp
dev.off()




####### tran2 vs tran 3
Transcript_1 <- c(boots_fl_200[,2],boots_fl_800[,2],boots_fl_both[,2])
Transcript_2 <- c(boots_fl_200[,3],boots_fl_800[,3],boots_fl_both[,3])
fragmenet_length <- rep(c("200nt","800nt","200|800nt"),c(nrow(boots_fl_200),nrow(boots_fl_800),nrow(boots_fl_both)))

df <- data.frame(Transcript_1 = Transcript_1,Transcript_2 = Transcript_2,fragmen_length = fragmenet_length)
df <- rbind(df,data.frame(Transcript_1=1,Transcript_2=1,fragmen_length = "Expected Value"))
head(df)
colnames(df) <- c("Transcript 2","Transcript 3","Fragment Length")
df[nrow(df),]


df$`Fragment Length`
title_axis <- 14
text_size <- 12
mishape <- rep(16,nrow(df))
mishape[length(mishape)] <- 23
misize <- rep(1,nrow(df))
misize[length(misize)] <- 5

pp <- ggplot(df, aes(x=`Transcript 2`, y=`Transcript 3`, group=`Fragment Length`)) + 
  geom_point(aes(shape=`Fragment Length`,color=`Fragment Length`,size=`Fragment Length`,fill=`Fragment Length`)) + 
  scale_size_manual(values = c(2,2,2,4)) +
  scale_shape_manual(values = c(16,16,16,24)) +
  scale_color_manual(values = c("#00798c","#d1495b","#11773390","#000000")) +
  scale_fill_manual(values = c("#00798c","#d1495b","#7CAE0050","white")) +
  theme_bw() + ylim(0,1.25) + xlim(0,1.5)
pp

pp <- pp + theme(axis.title.x =element_text(size=title_axis),axis.title.y =element_text(size=title_axis),
                 axis.text = element_text(size=text_size),legend.text = element_text(size = text_size),legend.title = element_blank())
pp

jpeg(filename = "./output/twolibraries_simulation/tran2_vs_tran3_simulacion.jpeg")
pp
dev.off()















































































































































































