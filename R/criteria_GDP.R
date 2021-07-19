
### Functions used to find GDP-optimal designs (nearly optimal)
### Minimising. Components: Ds, DPs, LoF(DP), bias(D)
criteria.GDP<-function(X1,X2,eps=10^-23)      # X1, X2 -- matrices of primary and potential terms, both with labels
{
  Ds<-0; DP<-0; LoF<-0; bias<-0;
  DF<-nlevels(as.factor(X1[,1]))
  df<-Nruns-DF                                # df - pure error degrees of freedom

  M<-crossprod(X1[,-1])                       # information matrix of primary terms
  D<-prod(round(eigen(M,symmetric=TRUE,only.values=TRUE)$values,8))/Nruns
  if (D>eps)
  {
    Minv<-solve(M)
  } else {return (list (eval=0, Ds=0, DP=0, LoF=0, bias=0, df=df, compound=10^6));} # if M is comp. singular

  if ((kappa.Ds>0) || (kappa.DP>0))
  {
    if (D^(1.0/(P-1))>0)
    {
      Ds<-D^(-1.0/(P-1))
    } else {return (list (eval=0, Ds=0, DP=0, LoF=0, bias=0, df=df, compound=10^6));} # Ds
  }
  if (kappa.DP>0)
  {
    if (df>0)
    {
      DP<-Ds*qf(1-alpha.DP,P-1,df)
    } else {return (list (eval=0, Ds=0, DP=0, LoF=0, bias=0, df=df, compound=10^6));} # df=0
  }

  if (((kappa.LoF>0) && (df>0))||(kappa.bias>0))  # check for A calculation
  {
    M12<-crossprod(X1[,-1],X2[,-1])
    A<-Minv%*%M12
  }
  if (kappa.LoF>0)
  {
    if (df>0)
    {
      L0<-crossprod(X2[,-1])-t(M12)%*%A+diag(1./tau2,nrow=Q)
      LoF<-(det(L0))^(-1.0/Q)*qf(1-alpha.LoF,Q,df)
    } else {return (list (eval=0, Ds=0, DP=0, LoF=0, bias=0, df=df, compound=10^6));} # df=0
  }
  if (kappa.bias>0)
  {
    A0<-crossprod(A)+diag(1,nrow=Q)
    bias<-prod(round(eigen(A0,only.values=TRUE)$values,10))^(1.0/Q)
  }
  compound<-Ds^kappa.Ds*DP^kappa.DP*LoF^kappa.LoF*bias^kappa.bias
  list (eval=1,Ds=Ds, DP=DP, LoF=LoF, bias=bias, df=df, compound=compound)
}
