
extend<-function(candlt, K, Parameters)            # treatment matrix (with labels) -> extended design matrix (with labels)
{
  treat<-as.matrix(candlt[,-1])
  X<-matrix(1,nrow=nrow(candlt), ncol=1)  # intercept
  X<-cbind(X,treat[, Parameters[1:K]==1])
  count<-ncol(X)-1
  for (i in 1:K)
  {
    count<-count+1
    if (Parameters[count]>0) X<-cbind(X, treat[,i]^2)
  }

  if (K>=2)
  {
    for (i in 1:(K-1))
    {
      for (j in (i+1):K)
      {
        count<-count+1
        if (Parameters[count]>0) X<-cbind(X, treat[,i]*treat[,j])
      }
    }
  }

  X<-cbind(candlt[,1],X)
  return (X)
}
