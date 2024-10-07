#' Evaluating individual criteria of the designs, from the Generalized compound criteria (Goos et al., 2005), (Egorova, 2017)
#' @description Calculating values of determinant- and trace-based components of 
#' Generalized D-, DP-, L- and LP- criteria for an output of a search object, 
#' with model and control parameters set in a mood object.
#' 
#' @param search.obj Output of the `Search' function
#' @param mood.obj Output of the `mood' function
#' @param eps Computational tolerance, default 10^-20
#' 
#' @return List of the calculated values:
#' \itemize{
#' \item `df` pure error degrees of freedom
#' \item `Ds` Ds-criterion value, intercept excluded
#' \item `DP` DPs-criterion value, intercept excluded
#' \item `LoFD` LoF(D)-criterion value from the GD-criterion
#' \item `LoFDP` LoF(DP)-criterion value from the GDP-criterion
#' \item `biasD` bias(D)-criterion value from the GD-criterion
#' \item `Ls` Ls-criterion value, intercept excluded
#' \item `LP` LPs-criterion value, intercept excluded
#' \item `LoFL` LoF(L)-criterion value from the GL-criterion
#' \item `LoFLP` LoF(LP)-criterion value from the GLP-criterion
#' \item `biasL` bias(L)-criterion value from the GL-criterion
#' }

criteria.values.G<-function(search.obj, mood.obj, eps=10^-23)             
{
  DP<-0; LoFDP<-0; LP<-0; LoFLP<-0;
  X1 <- search.obj$X1
  X2 <- search.obj$X2
  df<- search.obj$df                     # df - pure error degrees of freedom
  
  P <- mood.obj$P; Q <- mood.obj$Q
  tau <- mood.obj$tau; tau2 <- mood.obj$tau2
  Nruns <- mood.obj$Nruns
  
  alpha.DP <- mood.obj$alpha.DP
  alpha.LP <- mood.obj$alpha.LP
  alpha.LoF <- mood.obj$alpha.LoF
  alpha.LoFL <- mood.obj$alpha.LoFL
  Z0 <- mood.obj$Z0
  W <- mood.obj$W

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
  
  Ls<-W[,-1]%*%(diag(Minv)[-1])               # Ls component

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
    DP<-Ds*qf(1-alpha.DP, P-1, df)                               # DP component
    LP<-Ls*qf(1-alpha.LP, 1, df)                                 # LP component
    LoFDP<-LoFD*qf(1-alpha.LoF, Q, df)                           # LoF (DP)
    L0.inv.trace<-sum(1./eigen(L0, only.values=TRUE)$values)
    LoFLP<-L0.inv.trace*qf(1-alpha.LoFL, 1, df)/Q                 # LoF (LP)
  } else {DP<-0;LoFDP<0; LoFLP<-0}

  list (df = df, 
        Ds=Ds, DP=DP, LoFD=LoFD, LoFDP=LoFDP, biasD=biasD, 
        Ls=Ls, LP=LP, LoFL=LoFL, LoFLP=LoFLP, biasL=biasL)
}
