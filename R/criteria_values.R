
## Functions calculating different criteria values
### Obtain simple criteria values: Ds, DP, L, LP (all non-orthogonal here)

criteria.values<-function(X1)
{

  ##if X1.orth extra functions
  if(orth=='Y')
  {
    cand<-candidate_set(Levels)
    labels<-as.vector(X1[,1])       # labels
    index<-rep(0,Nruns)
    for (i in 1:Nruns){
      index[i]<-which(cand[,1]==labels[i])
    }
    X1<-cand[index,]                     # extended matrix, not orthonormalized, with labels
  }###


  M<-crossprod(X1[,-1])
  D<-prod(round(eigen(M,symmetric=TRUE,only.values=TRUE)$values,6))
  Ds<-(D/Nruns)^(1.0/(P-1))
  Minv<-solve(M)
  Ls<-1.0/(W%*%(diag(Minv)[-1]))        # tr(W(X'X)^-1) - weighted A
  DF<-nlevels(as.factor(X1[,1]))        # df = "pure error" degrees of freedom
  df<-Nruns-DF
  if (df>0){
    DP<-Ds/qf(1-alpha.DP,P-1,df)
    LP<-Ls/qf(1-alpha.LP,1,df)
    LoF<-ifelse ((Nruns-P-df)>0,1.0/qf(1-alpha.LoF,Nruns-P-df,df),0)
  } else {DP<-0; LP<-0; LoF<-0;}
  compound<-(Ds^kappa.Ds)*(Ls^kappa.Ls)*(DP^kappa.DP)*(LP^kappa.LP)*(DF^kappa.DF)*(LoF^kappa.LoF)

  list (Ds=Ds, Ls=Ls, DP=DP, LP=LP, df=df, LoF=LoF, compound=compound)
}
