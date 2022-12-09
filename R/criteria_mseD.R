#' Calculates the values of the MSE DPs-criterion and its components
#'
#' This function evaluates the MSE DPs-criterion for given primary and potential model matrices. Candidate full model matrices do not have to be orthonormalised.
#' Components: DPs-, LoF(DP)- and MSE(D)-optimality.
#' The weights kappa.DP, kappa.LoF and kappa.mse, and other parameters (tau2, alpha-s) are taken from the global environment.
#'
#' @param X1 The primary model matrix, with the first column containing the labels of treatments, and the second -- the intercept term.
#' @param X2 The matrix of potential terms, with the first column containing the labels of treatments.
#' @param eps Computational tolerance, the default value is 10^-23
#'
#' @return A list of values: indicator of whether the evaluation was successful ("eval"), DPs-criterion value -- intercept excluded ("DP"),
#' Lack-of-fit(DP) criterion value ("LoF"), the MSE(D) component value ("mse"), 
#' the number of pure error degrees of freedom ("df") and the value of the compound criterion ("compound").
#' @export
#' @examples 
#'# Experiment: one 5-level factor, primary model -- full quadratic, X^3 and X^4 potential terms.
#' K<-1; P<-3; Q<-2;
#'Levels <- list(1:5)   
#'# Generating candidate sets: primary and full ones
#'cand.primary <- candidate_set(Levels);
#'cand.full <- cbind(cand.primary, cand.primary[,3]^3, cand.primary[,4]^2) #X^3 and X^4 potential terms
#'# Choosing a design
#'index <- c(rep(1,2),3,rep(4,2),rep(5,3)); Nruns<- length(index)
#'Z0<-diag(1,Nruns)-matrix(1/Nruns,nrow=Nruns,ncol=Nruns)
#'X.primary <- cand.full[index, 1:(P+1)]
#'X.potential <- cand.full[index, (c(1,(P+2):(P+Q+1)))]
#'# Evaluating a compound MSE(DP)-criterion
#'kappa.DP = kappa.LoF = kappa.mse = 1./3; tau2 <-1; tau <- sqrt(tau2)
#'alpha.DP = alpha.LoF = 0.05; Biter <- 1000; 
#'criteria.mseD(X1 = X.primary, X2 = X.potential)
#' # Output: eval = 1, DP = 2.682, LoF = 6.455, mse = .8854, df = 4, compound = 2.4842

criteria.mseD<-function(X1, X2, P, Q, Nruns, eps=10^-23)      # X1, X2 -- matrices of primary and potential terms, both with labels
{
  Ds<-0; DP<-0; LoF<-0; bias<-0; mse<-0;
  DF<-nlevels(as.factor(X1[,1]))
  df<-Nruns-DF                                # df - pure error degrees of freedom

  M<-crossprod(X1[,-1])                       # information matrix of primary terms
  D<-prod(round(eigen(M,symmetric=TRUE,only.values=TRUE)$values,8)) 
  if (D>eps)
  {
    Minv<-solve(M)
  } else {return (list (eval=0, DP=0, LoF=0, bias=0, df=df, compound=10^6));} # if M is computationally singular
  if (kappa.Ds>0 || (kappa.DP>0))  
  {
    if (D^(1.0/(P-1))>0)
    {
      Ds<-(D/Nruns)^(-1.0/(P-1))
    } else {return (list (eval=0, DP=0, LoF=0, bias=0, df=df, compound=10^6));} # Ds
  }
  if (kappa.DP>0)
  {
    if (df>0)
    {
      DP<-Ds*qf(1-alpha.DP,P-1,df)             # DPs
    } else {return (list (eval=0, DP=0, LoF=0, bias=0, df=df, compound=10^6));} # if df=0
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
      LoF<-(det(L0))^(-1.0/Q)*qf(1-alpha.LoF,Q,df) #Ds
    } else {return (list (eval=0, DP=0, LoF=0, bias=0, df=df, compound=10^6));} # if df=0
  }
  if (kappa.mse>0)
  {
    Tvalue<-0
    M120<-t(X1[,-c(1,2)])%*%Z0%*%X2[,-1]
    for (j in 1:Biter)                                         # MC estimation of the second part of the MSE(D)-component
    {
      beta2<-rnorm(Q,mean=0,sd=tau)
      M12b<-M120%*%beta2
      Evalue<-t(M12b)%*%(Minv[-1,-1])%*%M12b
      Tvalue<-Tvalue+log(1+Evalue)
    }
    MC<-Tvalue/Biter
    mse<-(exp(MC)*Nruns/D)^(1./(P-1))                                 # MSE(Ds)
  }
  compound<-DP^kappa.DP*LoF^kappa.LoF*mse^kappa.mse
  list (eval=1, DP=DP, LoF=LoF, mse=mse, df=df, compound=compound)
}
