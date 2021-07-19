
### Functions used to find unblocked MSE(D)-optimal designs (nearly optimal)
###           using the point prior for the mse(D)-component estimation

### Minimising. Components: DPs, LoF(DP), mse(D) (point estimation)
criteria.mseP<-function(X1,X2,eps=10^-23)      # X1, X2 -- matrices of primary and potential terms, with labels
{
  Ds<-0; DP<-0; LoF<-0; bias<-0;mse<-0;
  DF<-nlevels(as.factor(X1[,1]))
  df<-Nruns-DF                                # df - pure error degrees of freedom

  M<-crossprod(X1[,-1])                       # information matrix of primary terms
  D<-prod(round(eigen(M,symmetric=TRUE,only.values=TRUE)$values,8))
  if (D>eps)
  {
    Minv<-solve(M)
  } else {return (list (eval=0, Ds=0, DP=0, LoF=0, bias=0, df=df, compound=10^6));} # if M is computationally singular

  if ((kappa.Ds>0) || (kappa.DP>0))
  {
    if (D^(1.0/(P-1))>0)
    {
      Ds<-(D/Nruns)^(-1.0/(P-1))
    } else {return (list (eval=0, Ds=0, DP=0, LoF=0, bias=0, df=df, compound=10^6));} # Ds
  }
  if (kappa.DP>0)
  {
    if (df>0)
    {
      DP<-Ds*qf(1-alpha.DP,P-1,df)
    } else {return (list (eval=0, Ds=0, DP=0, LoF=0, bias=0, df=df, compound=10^6));} # if df=0
  }

  if (((kappa.LoF>0) && (df>0))||(kappa.mse>0))  # check for A calculation
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
    } else {return (list (eval=0, Ds=0, DP=0, LoF=0, bias=0, df=df, compound=10^6));} # if df=0
  }
  if (kappa.mse>0)
  {
    Tvalue<-0
    M12<-t(X1[,-c(1,2)])%*%Z0%*%X2[,-1]
    MM<-t(M12)%*%Minv[-1,-1]%*%M12
    beta2<-rep(tau,Q)                                        # prior point estimate = rep(tau,Q)
    Tvalue<-(1+t(beta2)%*%MM%*%beta2)
    mse<-(Tvalue*Nruns/D)^(1./(P-1))                         # MSE(D)_s point estimate
  }
  compound<-DP^kappa.DP*LoF^kappa.LoF*mse^kappa.mse
  list (eval=1,DP=DP, LoF=LoF, mse=mse, df=df, compound=compound)
}
