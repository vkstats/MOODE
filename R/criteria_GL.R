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
#' #'#Experiment: one 5-level factor, primary model -- full quadratic, one potential (cubic) term
#'K <-1; P<-3; Q<-1; Levels <- list(1:5)
#' # Generating candidate sets: primary and full orthonormalised ones
#'cand.primary <- candidate_set(Levels);
#'cand.not.orth <-cbind(cand.primary[,-1], cand.primary[,3]*cand.primary[,4])
#'cand.full.orth <- cbind(cand.primary[,1], far::orthonormalization(cand.not.orth,basis=FALSE))
#' # Choosing a design
#'index <- c(rep(1,2),3,4, rep(5,3)); Nruns<- length(index)
#'X.primary <- cand.full.orth[index, 1:(P+1)]
#'X.potential <- cand.full.orth[index, (c(1,(P+2):(P+Q+1)))]
#' # Evaluating a compound GL-criterion
#'kappa.Ls = kappa.LoF = kappa.bias = 1./3; tau2 <-1;
#'W = matrix(c(0.8, 0.2), nrow = 1)
#'criteria.GL(X1 = X.primary, X2 = X.potential)
#' # Output: eval = 1, Ls = .5613, LoF = .7213, bias = 1.4331, df = 3, compound = .8418


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
  } else {return (list (eval=0, Ls=0, LP=0, LoF=0, bias=0, df=df, compound=10^6));}

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
  list (eval=1, Ls=Ls, LoF=LoF, bias=bias, df=df, compound=compound)
}
