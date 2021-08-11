

initial.full<-function(cand.full)          # cand.full - extended model matrix with labels, returns initial primary and potential matrices
{
  eps<-10^(-6)
  det<-0
  while (det<eps)
  {
    index<-sample(1:nrow(cand.full),size=Nruns,replace=TRUE)
    X1<-cand.full[index,1:(P+1)]
    X2<-cand.full[index,c(1,(P+2):(P+Q+1))]
    det<-round(prod(eigen(t(X1[,-1])%*%X1[,-1],symmetric=TRUE,only.values=TRUE)$values),6)  # det of the information matrix
  }
  list(X1=X1, X2=X2, det=det)  # X1, X2 -- matrices of primary and potential terms, both with labels
}