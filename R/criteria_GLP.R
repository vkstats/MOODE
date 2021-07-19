
### Functions used to find GLP-optimal designs (nearly optimal)

### Minimising Components: Ls, LP, LoF(LP), bias(L)

criteria.GLP<-function(X1,X2,eps=10^-23)      # X1, X2 -- matrices of primary and potential terms, both with labels
{
  Ls<-0; LP<-0; LoF<-0; bias<-0;
  DF<-nlevels(as.factor(X1[,1]))
  df<-Nruns-DF                                # df - pure error degrees of freedom

  M<-crossprod(X1[,-1])                       # information matrix of primary terms
  D<-prod(round(eigen(M,symmetric=TRUE,only.values=TRUE)$values,8))/Nruns
  if (D>eps)
  {
    Minv<-solve(M)
  } else {return (list (Ls=0, LP=0, LoF=0, bias=0, df=df, compound=10^6));}

  if ((kappa.Ls>0)||(kappa.LP>0))
  {
    Ls<-W%*%(diag(Minv)[-1])
  }
  if (kappa.LP>0)
  {
    if (df>0)
    {
      LP<-Ls*qf(1-alpha.LP,1,df)
    } else {return (list (Ls=0, LP=0, LoF=0, bias=0, df=df, compound=10^6));}
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
      L0<-crossprod(X2[,-1])-t(M12)%*%A+diag(1./tau2,nrow=Q)         # dispersion matrix + Iq/tau2
      L0.inv.trace<-Re(sum(1./eigen(L0,only.values=TRUE)$values))    # trace of the inverse matrix
      LoF<-L0.inv.trace*qf(1-alpha.LoFL,1,df)/Q
    } else {return (list (Ls=0, LP=0, LoF=0, bias=0, df=df, compound=10^6));}
  }
  if (kappa.bias>0)
  {
    A0<-crossprod(A)+diag(1,nrow=Q)
    bias<-sum(diag(A0))/Q                                        # (averaged) trace of the A'A+Iq
  }

  compound<-Ls^kappa.Ls*LP^kappa.LP*LoF^kappa.LoF*bias^kappa.bias
  list (Ls=Ls, LP=LP, LoF=LoF, bias=bias, df=df, compound=compound)
}

