
### Generating initial design

initial<-function(cand,Nruns=Nruns)          # cand - matrix with labels, primary terms
{
  eps<-10^(-6)
  det<-0
  while (det<eps)
  {
    index<-sample(1:nrow(cand),size=Nruns,replace=TRUE)
    X<-cand[index,]
    det<-round(prod(eigen(t(X[,-1])%*%X[,-1],symmetric=TRUE,only.values=TRUE)$values),6) # det of the information matrix
  }
  list(X=X, det=det)
}
