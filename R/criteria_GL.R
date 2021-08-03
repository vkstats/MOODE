#' Calculates the values of the Generalised Ls-criterion and its components
#'
#' This function evaluates the Generalised Ls-criterion (Goos et al., 2005) for given primary and potential model matrices.
#' The weights kappa.Ls, kappa.LoF and kappa.bias are taken from the global environment.
#'
#' @param X1 The primary model matrix, with the first column containing the labels of treatments, and the second -- the intercept term.
#' @param X2 The matrix of potential terms, with the first column containing the labels of treatments.
#' @param eps Computational tolerance, the default value is 10^-23
#'
#' @return A list of values: indicator of whether the evaluation was successful ("eval"), Ls-criterion value -- intercept excluded ("Ls"),
#' Lack-of-fit criterion value ("LoF"), the bias component value ("bias"), the number of pure error degrees of freedom ("df") 
#' and the value of the compound criterion ("compound").
#' @export
#' @examples 
#' kappa.Ls = kappa.LoF = kappa.bias = 1./3;
#' cand.primary <- candidate_set(rep(list(1:3), 2));
#' index <- c(1,2,4,6,7); P <- 6; Q <- 2;
#' X.primary <- cand.primary[index,1:(P+1)]
#' X.potential <- cbind(X.primary[,1], X.primary[,2]^3, X.primary[,3]^3) # cubic terms
#' criteria.GL(X1 = X.primary, X2 = X.potential)

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
