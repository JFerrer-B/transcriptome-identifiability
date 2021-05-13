Get_Fragment_Distribution<-function(mis_lr_fl_STD_AMOUNT){
  # This function creates a set of fragment length distributions following a Gaussian with a set of inputs
  # amount = Total number of final fragment lengths from distribution
  # fl = mean length of the distribution to create
  # std = standard deviation of the distribution
  # MaxFl= Maximum value of the selected fragment length from the distribution
  # MinFl= Minimum value of the selected fragment length from the distribution
  
  amount <- mis_lr_fl_STD_AMOUNT[,4]
  fl <- mis_lr_fl_STD_AMOUNT[,c(1,2),drop=FALSE]
  std <- mis_lr_fl_STD_AMOUNT[,3]
  maxfl=4*std
  minfl=-4*std
  
  
  # D<-dim(fl)
  # minfl<-rep(MinFl,D[1])
  # maxfl<-rep(MaxFl,D[1])
  
  #To avoid negative values of fragment length(fl)
  minfl[(fl[,2]+minfl)<fl[,1]]<-fl[(fl[,2]+minfl)<fl[,1],1]+1-fl[(fl[,2]+minfl)<fl[,1],2]
  
  #Contains maximum and minimum values for each lr fl pair
  fl2<-cbind(fl[,2]+minfl,fl[,2]+maxfl) 
  
  #Generate distribution
  seq2 <- Vectorize(seq.default, vectorize.args = c("from", "to"))
  
  Sequence <- sapply(1:nrow(fl2),function(i){
    cbind(seq2(from = fl2[i,1], to = fl2[i,2],length.out = amount[i]))
  })
  # Sequence<-cbind(seq2(from = fl2[,1], to = fl2[,2],length.out = amount))
  
  #Create the weights according to a normal distr
  Density<-t(dnorm(t(Sequence),mean=fl[,2],sd=std)) 
  Normalized<-t(Density)/apply(Density,2,sum) 
  #And normalize them, so that they sum up to 1 for the posterior convolution
  Normalized<-t(Normalized) 
  
  #Name the col corresponding to lr and fl
  # colnames(Sequence)<-paste(fl[,1],fl[,2],sep="-") 
  # colnames(Normalized)<-paste(fl[,1],fl[,2],sep="-")
  colnames(Sequence) <- paste(mis_lr_fl_STD_AMOUNT[,1],mis_lr_fl_STD_AMOUNT[,2],mis_lr_fl_STD_AMOUNT[,3],mis_lr_fl_STD_AMOUNT[,4],sep="-") 
  colnames(Normalized) <- paste(mis_lr_fl_STD_AMOUNT[,1],mis_lr_fl_STD_AMOUNT[,2],mis_lr_fl_STD_AMOUNT[,3],mis_lr_fl_STD_AMOUNT[,4],sep="-") 
  
  #Save data
  FragmentData<-list(fl=Sequence,weight=Normalized)
  return(FragmentData)
  
}
