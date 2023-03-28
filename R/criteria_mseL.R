#' Calculates the values of the MSE LPs-criterion and its components
#'
#' This function evaluates the MSE LPs-criterion for given primary and potential model matrices. Candidate full model matrices do not have to be orthonormalised.
#' Components: LPs-, LoF(LP)- and MSE(L)-optimality.
#' The weights kappa.LP, kappa.LoF and kappa.mse, and other parameters (tau2, alpha-s) are taken from the global environment.
#'
#' @param X1 The primary model matrix, with the first column containing the labels of treatments, and the second -- the intercept term.
#' @param X2 The matrix of potential terms, with the first column containing the labels of treatments.
#' @param eps Computational tolerance, the default value is 10^-23
#'
#' @return A list of values: indicator of whether the evaluation was successful ("eval"), LPs-criterion value -- intercept excluded ("LP"),
#' Lack-of-fit(LP) criterion value ("LoF"), the MSE(L) component value ("mse"), 
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
#'# Evaluating a compound MSE(LP)-criterion
#'kappa.LP = kappa.LoF = kappa.mse = 1./3; tau2 <-1; tau <- sqrt(tau2)
#'alpha.LP = alpha.LoFL = 0.05;  
#'criteria.mseL(X1 = X.primary, X2 = X.potential)
#' # Output: eval = 1, LP = 2.3863, LoF = 7.1846, mse = 1.5994, df = 4, compound = 2.8047
#' 
criteria.mseL<-function(X1, X2, search.object, eps=10^-23)      # X1, X2 -- matrices of primary and potential terms, both with labels
{
  Ls<-0; LP<-0; LoF<-0; mse<-0;
  DF<-nlevels(as.factor(X1[,1]))  
  
  Nruns<-search.object$Nruns
  P<-search.object$P; Q<-search.object$Q
  tau2<-search.object$tau2
  W<-search.object$W
  
  kappa.Ls<-search.object$kappa.Ls
  kappa.LP<-search.object$kappa.LP
  kappa.LoF<-search.object$kappa.LoF
  kappa.mse<-search.object$kappa.mse
  
  alpha.LP<-search.object$alpha.LP
  alpha.LoFL<-search.object$alpha.LoFL
  
  df<-Nruns-DF                                # df - pure error degrees of freedom

  M<-crossprod(X1[,-1])                       # information matrix of primary terms
  D<-prod(round(eigen(M,symmetric=TRUE,only.values=TRUE)$values,8))
  if (D>eps)
  {
    Minv<-solve(M)
  } else {return (list (eval=0, LP=0, LoF=0, mse=0, df=df, compound=10^6));}

  if ((kappa.Ls>0)||(kappa.LP>0))
  {
    Ls<-W%*%(diag(Minv)[-1])                   # Ls
  }
  if (kappa.LP>0)
  {
    if (df>0)
    {
      LP<-Ls*qf(1-alpha.LP,1,df)              # LPs
    } else {return (list (eval=0, LP=0, LoF=0, mse=0, df=df, compound=10^6));}
  }
  if (((kappa.LoF>0) && (df>0))||(kappa.mse>0))         # check for A calculation
  {
    M12<-crossprod(X1[,-1],X2[,-1])
    A<-Minv%*%M12
  }
  if (kappa.LoF>0)
  {
    if (df>0)
    {
      L0<-crossprod(X2[,-1])-t(M12)%*%A+diag(1./tau2,nrow=Q)            # dispersion matrix + Iq/tau2
      L0.inv.trace<-Re(sum(1./eigen(L0,only.values=TRUE)$values))       # trace of the inverse matrix
      LoF<-L0.inv.trace*qf(1-alpha.LoFL,1,df)/Q
    } else {return (list (eval=0, LP=0, LoF=0, mse=0, df=df, compound=10^6));}
  }
  if (kappa.mse>0)
  {
    A0<-Minv+tau2*tcrossprod(A)
    mse<-sum(diag(A0)[-1])/(P-1)                                        # averaged trace of the A'A+Iq
  }
  compound<-LP^kappa.LP*LoF^kappa.LoF*mse^kappa.mse
  list (eval=1, LP=LP, LoF=LoF, mse=mse, df=df, compound=compound)
}
