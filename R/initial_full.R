

initial.full<-function(cand.full, Nruns,
                       primary.terms, potential.terms)          # cand.full - extended model matrix with labels, returns initial primary and potential matrices
{
  eps<-10^(-6)
  det<-0
  while (det<eps)
  {
    index<-sample(1:nrow(cand.full), size=Nruns, replace=TRUE)
    X1<-cand.full[index, c("label", primary.terms), drop = F]
#    cat(index, "\n")
#    cat(potential.terms, "\n")
#    print(cand.full)
    X2<-cand.full[index, c("label", potential.terms), drop = F]
    det<-round(prod(eigen(t(X1[,-1])%*%X1[,-1], symmetric=TRUE, only.values=TRUE)$values), 6)  # det of the information matrix
  }
  list(X1=X1, X2=X2, det=det)  # X1, X2 -- matrices of primary and potential terms, both with labels
}
