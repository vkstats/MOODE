
### Obtaining the  design matrix, no lables, transformed to [-1,1]

extract.design<-function(X1)
{
  labels<-as.vector(X1[,1])                        # extract labels
  cand<-as.matrix(expand.grid(Levels))
  candl<-cbind(label(cand),cand)                 # labelling
  index<-rep(0,Nruns)
  for (i in 1:Nruns)
  {
    index[i]<-which(candl[,1]==labels[i])
  }
  design<-apply(cand[index,],2,transform)
}
