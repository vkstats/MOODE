
##### GD, GDP GL, GLP criteria evaluation

criteria.values.G<-function(X1,X2,eps=10^-23)             # X1 - orthonormalised matrix
{
  DP<-0; LoFDP<-0; LP<-0; LoFLP<-0;
  DF<-nlevels(as.factor(X1[,1]))              # d.f. = N-number of unique design points
  df<-Nruns-DF                                # df - pure error degrees of freedom

  M<-crossprod(X1[,-1])                       # information matrix of primary terms
  D<-prod(round(eigen(M,symmetric=TRUE,only.values=TRUE)$values,8))/Nruns
  if (D>eps)
  {
    Minv<-solve(M)
  } else {return (eval=0)}
  if (D^(1.0/(P-1))>0)
  {
    Ds<-D^(-1.0/(P-1))                        # Ds component
  } else {return (eval=0)}
  Ls<-W%*%(diag(Minv)[-1])                    # Ls component

  M12<-crossprod(X1[,-1],X2[,-1])
  A<-Minv%*%M12
  L0<-crossprod(X2[,-1])-t(M12)%*%A+diag(1./tau2,nrow=Q)
  LoFD<-(det(L0))^(-1.0/Q)                                          # LoF (D)
  LoFL<-1./(sum(diag(L0))/Q)                                        # LoF (L)

  A0<-crossprod(A)+diag(1,nrow=Q)
  biasD<-prod(round(eigen(A0,only.values=TRUE)$values,10))^(1.0/Q)  # bias (D)
  biasL<-sum(diag(A0))/Q                                            # bias (L)

  if (df>0)
  {
    DP<-Ds*qf(1-alpha.DP,P-1,df)                               # DP component
    LP<-Ls*qf(1-alpha.LP,1,df)                                 # LP component
    LoFDP<-LoFD*qf(1-alpha.LoF,Q,df)                           # LoF (DP)
    L0.inv.trace<-sum(1./eigen(L0,only.values=TRUE)$values)
    LoFLP<-L0.inv.trace*qf(1-alpha.LoFL,1,df)/Q                 # LoF (LP)
  } else {DP<-0;LoFDP<0; LoFLP<-0}

  list (Ds=Ds,DP=DP,LoFD=LoFD,LoFDP=LoFDP,biasD=biasD, Ls=Ls,LP=LP,LoFL=LoFL,LoFLP=LoFLP,biasL=biasL)
}
