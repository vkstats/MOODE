#' Evaluating individual criteria of the designs
#' @description Calculating values of determinant- and trace-based components of 
#' MSE(D)- and MSE(L)- criteria for an output of a search object, 
#' with model and control parameters set in a mood object.
#' 
#' @param search.obj Output of the `Search' function
#' @param mood.obj Output of the `mood' function
#' @param eps Computational tolerance, default 10^-20
#' @param Biter MC sample size for evaluating the mse(D)-component
#' 
#' @return List of the calculated values:
#' \itemize{
#' \item `df` pure error degrees of freedom
#' \item `Ds` Ds-criterion value, intercept excluded
#' \item `DP` DPs-criterion value, intercept excluded
#' \item `LoFDP` LoF(DP)-criterion value
#' \item `mseD` mse(D)-criterion value, obtained via MC sampling
#' \item `mseP` mse(D)-criterion value, obtained using point prior 
#' \item `Ls` Ls-criterion value, intercept excluded
#' \item `LP` LPs-criterion value, intercept excluded
#' \item `LoFLP` LoF(LP)-criterion value
#' \item `mseL` mse(L)-criterion value
#' }
#' @export

criteria.values.mse<-function(search.obj, mood.obj, eps=10^-20, Biter=1000)      # X1,X2 - not orthonormalised matrices
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

  M<-crossprod(X1[,-1])                  # information matrix of primary terms
  D<-prod(round(eigen(M, symmetric=TRUE, only.values=TRUE)$values,8))
  if (D>eps)
  {
    Minv<-solve(M)
  } else {return (eval=0)}
  if (D^(1.0/(P-1))>0)
  {
    Ds<-(D/Nruns)^(-1.0/(P-1))                        # Ds component
  } else {return (eval=0)}

  Ls<-W[,-1]%*%(diag(Minv)[-1])                            # Ls component

  M12<-crossprod(X1[,-1],X2[,-1])
  A<-Minv%*%M12
  L0<-crossprod(X2[,-1])-t(M12)%*%A+diag(1./tau2,nrow=Q)
  LoFD<-(det(L0))^(-1.0/Q)                                          # LoF (D)
  LoFL<-1./(sum(diag(L0))/Q)                                        # LoF (L)

  if (df>0)
  {
    DP<-Ds * stats::qf(1-alpha.DP,P-1,df)                               # DPs component
    LP<-Ls * stats::qf(1-alpha.LP, 1, df)                                 # LPs component
    LoFDP<-LoFD * stats::qf(1-alpha.LoF, Q, df)                           # LoF (DP)
    L0.inv.trace<-sum(1./eigen(L0, only.values=TRUE)$values)
    LoFLP<-L0.inv.trace * stats::qf(1-alpha.LoFL, 1, df)/Q                # LoF (LP)
  } else {DP<-0;LoFDP<0; LoFLP<-0}

  mseL<-sum(diag(Minv+tau2*tcrossprod(A))[-1])/(P-1)                   # MSE(Ls)

  Tvalue<-0
  M120<-t(X1[,-c(1,2)])%*%Z0%*%X2[,-1]
  for (j in 1:Biter)                              # random beta-s, MC for the expectation
  {
    beta2 <- stats::rnorm(Q,mean=0,sd=tau)
    M12b<-M120%*%beta2
    Evalue<-t(M12b)%*%(Minv[-1,-1])%*%M12b
    Tvalue<-Tvalue+log(1+Evalue)
  }
  MC<-Tvalue/Biter
  mseD<-(exp(MC)*Nruns/D)^(1./(P-1))                           # MSE(Ds)
  
  MM <- t(M120)%*%(Minv[-1,-1])%*%M120
  betaP <- rep(tau, Q)                                        # prior point estimate = rep(tau,Q)
  TvalueP <-(1+t(betaP)%*%MM%*%betaP)
  mseP <- (Tvalue*Nruns/D)^(1./(P-1))                         # MSE(D)_s_point

  list (df = df, Ds=Ds, DP=DP, LoFDP=LoFDP, mseD=mseD, mseP=mseP, 
        Ls=Ls, LP=LP, LoFLP=LoFLP,mseL = mseL)
}
