#' Calculates the values of the Generalised DPs-criterion and its components
#'
#' This function evaluates the Generalised DPs-criterion for given primary and potential model matrices. 
#' Components: Ds-, DPs-, LoF(DP)- and Bias(D)-optimality.
#'
#' @param X1 The primary model matrix, with the first column containing the labels of treatments, and the second -- the intercept term.
#' @param X2 The matrix of potential terms, with the first column containing the labels of treatments.
#' @param search.object Object of class [mood()] specifying experiment parameters.
#' @param eps Computational tolerance, the default value is 10^-23
#'
#' @return A list of values: indicator of whether the evaluation was successful ("eval"), Ds-criterion value -- intercept excluded ("Ds"),
#' DPs-criterion value -- intercept excluded ("DPs"), Lack-of-fit(DP) criterion value ("LoF"), the bias component value ("bias"), 
#' the number of pure error degrees of freedom ("df") and the value of the compound criterion ("compound").
#' @export
#' @examples 
#' # Experiment: one 5-level factor, primary model -- full quadratic, X^3 and X^4 potential terms.
#' ex.mood <- mood(K = 1, Levels = 5, Nruns = 8, criterion.choice = "GDP", 
#'                kappa = list(kappa.Ds = .25, kappa.LoF = .25, kappa.bias = .25, kappa.DP = .25), 
#'                model_terms = list(primary.model = "second_order", potential.terms = "x14"))
#' # Generating candidate sets: primary and full orthonormalised ones
#' K <- ex.mood$K
#' Levels <- ex.mood$Levels 
#' cand.not.orth <- candidate_set_full(candidate_trt_set(Levels, K), K)
#' cand.full.orth <- candidate_set_orth(cand.not.orth, ex.mood$primary.terms, ex.mood$potential.terms)
#' # Choosing a design
#' index <- c(rep(1,2),3,rep(4,2),rep(5,3))
#' X.primary <- cand.full.orth[index, c(1, match(ex.mood$primary.terms, colnames(cand.full.orth)))]
#' X.potential <- cand.full.orth[index, (c(1, match(ex.mood$potential.terms, colnames(cand.full.orth))))]
#' # Evaluating a compound GDP-criterion
#' criteria.GDP(X1 = X.primary, X2 = X.potential, ex.mood)
#' 
#' # Output: eval = 1, Ds = 0.6884783, DP = 4.4538023, LoF = 3.895182, bias = 1.03807, df = 4, compound = 2.465318
#'

criteria.GDP<-function(X1, X2, search.object, eps=10^-23)      # X1, X2 -- matrices of primary and potential terms, both with labels
{
  Ds<-0; DP<-0; LoF<-0; bias<-0;
  DF<-nlevels(as.factor(X1[,1]))
  
  Nruns<-search.object$Nruns
  P<-search.object$P; Q<-search.object$Q
  tau2<-search.object$tau2
  
  kappa.Ds<-search.object$kappa.Ds
  kappa.DP<-search.object$kappa.DP
  kappa.LoF<-search.object$kappa.LoF
  kappa.bias<-search.object$kappa.bias
  
  alpha.DP<-search.object$alpha.DP
  alpha.LoF<-search.object$alpha.LoF
  
  df<-Nruns-DF # df - pure error degrees of freedom

  M<-crossprod(X1[,-1])                       # information matrix of primary terms
  
  # correction to remove intercept and do Ds?
  # M <- t(X1[, -1]) %*% search.object$Z0 %*% X1[, -1] 
  
  # kept this form, using the *full* info matrix and dividing by the sum of squared intercept column (which could be orth.) 
  # to make the other criteria here work
  D<-prod(round(eigen(M, symmetric=TRUE, only.values=TRUE)$values,8)) / sum(X1[, 2]^2) 
  if (D>eps)
  {
    Minv<-solve(M)
  } else {return (list (eval=0, Ds=0, DP=0, LoF=0, bias=0, df=df, compound=10^6));} # if M is comp. singular

  if ((kappa.Ds>0) || (kappa.DP>0))
  {
    if (D^(1.0/(P-1))>0)
    {
      Ds<-D^(-1.0/(P-1))
    } else {return (list (eval=0, Ds=0, DP=0, LoF=0, bias=0, df=df, compound=10^6));} # Ds
  }
  if (kappa.DP>0)
  {
    if (df>0)
    {
      DP<-Ds*qf(1-alpha.DP,P-1,df)
    } else {return (list (eval=0, Ds=0, DP=0, LoF=0, bias=0, df=df, compound=10^6));} # df=0
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
      L0<-crossprod(X2[,-1])-t(M12)%*%A+diag(1./tau2,nrow=Q)
      LoF<-(det(L0))^(-1.0/Q)*qf(1-alpha.LoF,Q,df)
    } else {return (list (eval=0, Ds=0, DP=0, LoF=0, bias=0, df=df, compound=10^6));} # df=0
  }
  if (kappa.bias>0)
  {
    A0<-crossprod(A)+diag(1,nrow=Q)
    bias<-prod(round(eigen(A0,only.values=TRUE)$values,10))^(1.0/Q)
  }
  compound<-Ds^kappa.Ds*DP^kappa.DP*LoF^kappa.LoF*bias^kappa.bias
  list (eval=1,Ds=Ds, DP=DP, LoF=LoF, bias=bias, df=df, compound=compound)
}

