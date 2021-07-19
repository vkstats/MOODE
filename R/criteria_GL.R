
### Functions used to find GL-optimal designs (nearly optimal)

### Minimising Components: Ls, LoF(L), bias(L)

criteria.GL<-function(X1,X2,eps=10^-23)      # X1, X2 -- matrices of primary and potential terms, both with labels
{
  Ls<-0; LoF<-0; bias<-0;
  DF<-nlevels(as.factor(X1[,1]))
  df<-Nruns-DF                                # df - pure error degrees of freedom

  M<-crossprod(X1[,-1])                       # information matrix of primary terms
  D<-prod(round(eigen(M,symmetric=TRUE,only.values=TRUE)$values,8))/Nruns
  if (D>eps)
  {
    Minv<-solve(M)
  } else {return (list (Ls=0, LP=0, LoF=0, bias=0, df=df, compound=10^6));}

  if (kappa.Ls>0)
  {
    Ls<-W%*%(diag(Minv)[-1])
  }

  if ((kappa.LoF>0)||(kappa.bias>0))           # check for A calculation
  {
    M12<-crossprod(X1[,-1],X2[,-1])
    A<-Minv%*%M12
  }
  if (kappa.LoF>0)
  {
    L0<-crossprod(X2[,-1])-t(M12)%*%A+diag(1./tau2,nrow=Q)     # dispersion matrix + Iq/tau2
    LoF<-1./(sum(diag(L0))/Q)                                  # inverse of the (averaged) trace of the matrix
  }
  if (kappa.bias>0)
  {
    A0<-crossprod(A)+diag(1,nrow=Q)
    bias<-sum(diag(A0))/Q                                        # (averaged) trace of the A'A+Iq
  }

  compound<-Ls^kappa.Ls*LoF^kappa.LoF*bias^kappa.bias
  list (Ls=Ls, LoF=LoF, bias=bias, df=df, compound=compound)
}
