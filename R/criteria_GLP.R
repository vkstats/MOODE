#' Calculates the values of the Generalised LPs-criterion and its components
#'
#' This function evaluates the Generalised LPs-criterion for given primary and potential model matrices. 
#' Components: Ls-, LPs-, LoF(LP)- and Bias(L)-optimality.
#'
#' @param X1 The primary model matrix, with the first column containing the labels of treatments, and the second -- the intercept term.
#' @param X2 The matrix of potential terms, with the first column containing the labels of treatments.
#' @param search.object Object of class [mood()] specifying experiment parameters.
#' @param eps Computational tolerance, the default value is 10^-23
#'
#' @return A list of values: indicator of whether the evaluation was successful ("eval"), Ls-criterion value -- intercept excluded ("Ls"),
#' LPs-criterion value -- intercept excluded ("LPs"), Lack-of-fit(LP) criterion value ("LoF"), the bias component value ("bias"), 
#' the number of pure error degrees of freedom ("df") and the value of the compound criterion ("compound").
#' @export
#' @examples 
#'#' # Experiment: one 5-level factor, primary model -- full quadratic, X^3 and X^4 potential terms.
#' ex.mood <- mood(K = 1, Levels = 5, Nruns = 8, criterion.choice = "GLP", 
#'                kappa = list(kappa.Ls = .25, kappa.LoF = .25, kappa.bias = .25, kappa.LP = .25), 
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
#' criteria.GLP(X1 = X.primary, X2 = X.potential, ex.mood)
#'# Output: eval = 1, Ls = 0.2952603, LP = 4.584705, LoF = 3.895182, 
#'# bias = 1.03807, df = 4, compound = 1.529564
#'

criteria.GLP<-function(X1, X2, search.object, eps = 1e-23)      # X1, X2 -- matrices of primary and potential terms, both with labels
{
  Ls<-0; LP<-0; LoF<-0; bias<-0;
  DF<-nlevels(as.factor(X1[,1]))
  
  Nruns<-search.object$Nruns
  P<-search.object$P; Q<-search.object$Q
  tau2<-search.object$tau2
  W<-search.object$W
  
  kappa.Ls<-search.object$kappa.Ls
  kappa.LP<-search.object$kappa.LP
  kappa.LoF<-search.object$kappa.LoF
  kappa.bias<-search.object$kappa.bias
  
  alpha.LP<-search.object$alpha.LP
  alpha.LoFL<-search.object$alpha.LoFL
  
  df<-Nruns-DF                                # df - pure error degrees of freedom

  # kept this form, using the *full* info matrix and dividing by the sum of squared intercept column (which could be orth.)
  M<-crossprod(X1[,-1])                       # information matrix of primary terms
  D<-prod(round(eigen(M, symmetric=TRUE, only.values=TRUE)$values,8)) / sum(X1[, 2]^2)
  if (D>eps)
  {
    Minv<-solve(M)
  } else {return (list (eval=0,Ls=0, LP=0, LoF=0, bias=0, df=df, compound=10^6));}

  if ((kappa.Ls>0)||(kappa.LP>0))
  {
    # to remove intercept as nuisance parameter  
    M0 <- search.object$Z0 %*% X1[, -c(1, 2)]
    M0 <- crossprod(X1[, -c(1, 2)], M0)
    M0inv <- solve(M0)
    Ls<-W[, -1] %*% diag(M0inv)
  }
  if (kappa.LP>0)
  {
    if (df>0)
    {
      LP <- Ls * stats::qf(1-alpha.LP,1,df)
    } else {return (list (eval=0, Ls=Ls, LP=0, LoF=0, bias=0, df=df, compound=10^6));}
  }
  if (((kappa.LoF>0) && (df>0))||(kappa.bias>0))  # check for A calculation
  {
    M12<-crossprod(X1[,-1],X2[,-1])
    A<-Minv%*%M12
  }
  if (kappa.LoF>0)
  {
    if (df>0)
    {
      L0<-crossprod(X2[,-1])-t(M12)%*%A+diag(1./tau2,nrow=Q)         # dispersion matrix + Iq/tau2
      L0.inv.trace<-Re(sum(1./eigen(L0,only.values=TRUE)$values))    # trace of the inverse matrix
      LoF<-L0.inv.trace * stats::qf(1-alpha.LoFL,1,df)/Q
    } else {return (list (Ls=0, LP=0, LoF=0, bias=0, df=df, compound=10^6));}
  }
  if (kappa.bias>0)
  {
    A0<-crossprod(A)+diag(1,nrow=Q)
    bias<-sum(diag(A0))/Q                                        # (averaged) trace of the A'A+Iq
  }

  compound<-Ls^kappa.Ls*LP^kappa.LP*LoF^kappa.LoF*bias^kappa.bias
  list (eval=1, Ls=Ls, LP=LP, LoF=LoF, bias=bias, df=df, compound=compound)
}

