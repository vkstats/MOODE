#' Calculates the values of the Generalised Ls-criterion and its components
#'
#' This function evaluates the Generalised Ls-criterion (Goos et al., 2005) for given primary and potential model matrices.
#'
#' @param X1 The primary model matrix, with the first column containing the labels of treatments, and the second -- the intercept term.
#' @param X2 The matrix of potential terms, with the first column containing the labels of treatments.
#' @param search.object Object of class [mood()] specifying experiment parameters.
#' @param eps Computational tolerance, the default value is 10^-23
#'
#' @return A list of values: indicator of whether the evaluation was successful ("eval"), Ls-criterion value -- intercept excluded ("Ls"),
#' Lack-of-fit criterion value ("LoF"), the bias component value ("bias"), the number of pure error degrees of freedom ("df") 
#' and the value of the compound criterion ("compound").
#' @export
#' @examples 
#'#Experiment: one 5-level factor, primary model -- full quadratic, one potential (cubic) term
#'# setting up the example
#'ex.mood <- mood(K = 1, Levels = 5, Nruns = 7, criterion.choice = "GL", 
#'                kappa = list(kappa.Ls = 1./3, kappa.LoF = 1./3, kappa.bias = 1./3), 
#'                model_terms = list(primary.model = "second_order", potential.model = "cubic_terms"))
#'# Generating candidate set: orthonormalised
#'K <- ex.mood$K
#'Levels <- ex.mood$Levels 
#'cand.not.orth <- candidate_set_full(candidate_trt_set(Levels, K), K)
#'cand.full.orth <- candidate_set_orth(cand.not.orth, ex.mood$primary.terms, ex.mood$potential.terms)
#'# Choosing a design
#'index <- c(rep(1, 2), 3, 4, rep(5, 3))
#'X.primary <- cand.full.orth[index, c(1, match(ex.mood$primary.terms, colnames(cand.full.orth)))]
#'X.potential <- cand.full.orth[index, 
#'(c(1, match(ex.mood$potential.terms, colnames(cand.full.orth))))]
#'# Evaluating a compound GD-criterion
#'criteria.GL(X1 = X.primary, X2 = X.potential, ex.mood)
#' # Output: eval = 1, Ls = 0.3118626, LoF = 0.7212544, bias = 1.473138, df = 3, compound = 0.6919878
#'

criteria.GL<-function(X1, X2, search.object, eps=10^-23)      # X1, X2 -- matrices of primary and potential terms, both with labels
{
  Ls<-0; LoF<-0; bias<-0;
  DF<-nlevels(as.factor(X1[,1]))
  
  Nruns<-search.object$Nruns
  P<-search.object$P; Q<-search.object$Q
  tau2<-search.object$tau2
  W<-search.object$W
  
  kappa.Ls<-search.object$kappa.Ls
  kappa.LoF<-search.object$kappa.LoF
  kappa.bias<-search.object$kappa.bias
  
  df<-Nruns-DF                                # df - pure error degrees of freedom

  # kept this form, using the *full* info matrix and dividing by the sum of squared intercept column (which could be orth.)
  M<-crossprod(X1[,-1])                       # information matrix of primary terms
  D<-prod(round(eigen(M, symmetric=TRUE, only.values=TRUE)$values, 8)) / sum(X1[, 2]^2)
  if (D>eps)
  {
    Minv<-solve(M)
  } else {return (list (eval=0, Ls=0, LP=0, LoF=0, bias=0, df=df, compound=10^6));}

  if (kappa.Ls>0)
  {
    # to remove intercept as nuisance parameter  
    M0 <- search.object$Z0 %*% X1[, -c(1, 2)]
    M0 <- crossprod(X1[, -c(1, 2)], M0)
    M0inv <- solve(M0)
    Ls<-W[, -1] %*% diag(M0inv)
  }

  if ((kappa.LoF>0)||(kappa.bias>0))           # check for A calculation
  {
    M12<-crossprod(X1[,-1], X2[,-1])
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
