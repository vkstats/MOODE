#' Calculates the values of the MSE LP-criterion and its components
#'
#' This function evaluates the MSE LP-criterion for given primary and potential model matrices. Candidate full model matrices do not have to be orthonormalised.
#' Components: LP-, LoF(LP)- and MSE(L)-optimality.
#'
#' @param X1 The primary model matrix, with the first column containing the labels of treatments, and the second -- the intercept term.
#' @param X2 The matrix of potential terms, with the first column containing the labels of treatments.
#' @param search.object Object of class [mood()] specifying experiment parameters.
#' @param eps Computational tolerance, the default value is 10^-23
#'
#' @return A list of values: indicator of whether the evaluation was successful ("eval"), LP-criterion value -- intercept excluded ("LP"),
#' Lack-of-fit(LP) criterion value ("LoF"), the MSE(L) component value ("mse"), 
#' the number of pure error degrees of freedom ("df") and the value of the compound criterion ("compound").
#' @export
#' @examples 
#' #'# Experiment: one 5-level factor, primary model -- full quadratic, X^3 and X^4 potential terms.
#' ex.mood <- mood(K = 1, Levels = 5, Nruns = 8, criterion.choice = "MSE.L", 
#'                kappa = list(kappa.LP = 1./3, kappa.LoF = 1./3, kappa.mse = 1./3), 
#'                model_terms = list(primary.model = "second_order", potential.terms = "x14"))
#' # Generating candidate sets: primary and full orthonormalised ones
#' K <- ex.mood$K
#' Levels <- ex.mood$Levels 
#' cand.not.orth <- candidate_set_full(candidate_trt_set(Levels, K), K)
#' cand.full.orth <- candidate_set_orth(cand.not.orth, ex.mood$primary.terms, ex.mood$potential.terms)
#' # Choosing a design
#' index <- c(rep(1,2),3,rep(4,2),rep(5,3))
#' X.primary <- cand.full.orth[index, c(1, match(ex.mood$primary.terms, colnames(cand.full.orth)))]
#' X.potential <- cand.full.orth[index, 
#' (c(1, match(ex.mood$potential.terms, colnames(cand.full.orth))))]
#' # Evaluating a compound GDP-criterion
#' criteria.mseL(X.primary, X.potential, ex.mood)
#' # Output: eval = 1, LP = 4.584705, LoF = 3.895182, mse = 0.3926842, df = 4, compound = 1.914084

criteria.mseL<-function(X1, X2, search.object, eps=10^-23)      # X1, X2 -- matrices of primary and potential terms, both with labels
{
  Ls<-0; LP<-0; LoF<-0; mse<-0;
  DF<-nlevels(as.factor(X1[,1]))  
  
  Nruns<-search.object$Nruns
  P<-search.object$P; Q<-search.object$Q
  tau2<-search.object$tau2
  W<-search.object$W
  
  kappa.L<-search.object$kappa.L
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

  if ((kappa.L>0)||(kappa.LP>0))
  {
    Ls<-W[, -1]%*%(diag(Minv)[-1])                   # Ls
  }
  if (kappa.LP>0)
  {
    if (df>0)
    {
      LP<-Ls * stats::qf(1-alpha.LP,1,df)              # LPs
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
      LoF<-L0.inv.trace * stats::qf(1-alpha.LoFL,1,df)/Q
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
