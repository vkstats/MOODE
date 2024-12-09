#' Calculates the values of the MSE DPs-criterion using the point prior for the MSE(D)-component estimation
#'
#' This function evaluates the MSE DPs-criterion for given primary and potential model matrices, using point MSE(D)-component estimation. Candidate full model matrices do not have to be orthonormalised.
#' Components: DPs-, LoF(DP)- and MSE(D)-optimality.
#'
#' @param X1 The primary model matrix, with the first column containing the labels of treatments, and the second -- the intercept term.
#' @param X2 The matrix of potential terms, with the first column containing the labels of treatments.
#' @param search.object Object of class [mood()] specifying experiment parameters.
#' @param eps Computational tolerance, the default value is 10^-23
#'
#' @return A list of values: indicator of whether the evaluation was successful ("eval"), DPs-criterion value -- intercept excluded ("DP"),
#' Lack-of-fit(DP) criterion value ("LoF"), the MSE(D) component value ("mse"), 
#' the number of pure error degrees of freedom ("df") and the value of the compound criterion ("compound").
#' @export
#' @examples 
#' # Experiment: one 5-level factor, primary model -- full quadratic, X^3 and X^4 potential terms.
#' ex.mood <- mood(K = 1, Levels = 5, Nruns = 8, criterion.choice = "MSE.P", 
#'                kappa = list(kappa.DP = 1./3, kappa.LoF = 1./3, kappa.mse = 1./3), 
#'                model_terms = list(primary.model = "second_order", potential.terms = "x14"))
#' # Generating candidate sets: primary and full orthonormalised
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
#' criteria.mseP(X.primary, X.potential, ex.mood)
#' # Output: eval = 1, DP = 4.538023, LoF = 3.895182, mse = 0.6992699, df = 4, compound = 2.312135

criteria.mseP<-function(X1, X2, search.object, eps=10^-23)      # X1, X2 -- matrices of primary and potential terms, both with labels
{
  Ds<-0; DP<-0; LoF<-0; mse<-0;
  DF<-nlevels(as.factor(X1[,1]))
  
  Nruns<-search.object$Nruns
  P<-search.object$P; Q<-search.object$Q
  tau2<-search.object$tau2
  tau<-search.object$tau
  Z0<-search.object$Z0
  Biter<-search.object$Biter
  
  kappa.Ds<-search.object$kappa.Ds
  kappa.DP<-search.object$kappa.DP
  kappa.LoF<-search.object$kappa.LoF
  kappa.mse<-search.object$kappa.mse
  
  alpha.DP<-search.object$alpha.DP
  alpha.LoF<-search.object$alpha.LoF
  
  df<-Nruns-DF                                # df - pure error degrees of freedom

  M<-crossprod(X1[,-1])                       # information matrix of primary terms
  D<-prod(round(eigen(M,symmetric=TRUE,only.values=TRUE)$values,8)) / sum(X1[, 2]^2)
  if (D>eps)
  {
    Minv<-solve(M)
  } else {return (list (eval=0, DP=0, LoF=0, mse=0, df=df, compound=10^6));} # if M is computationally singular

  if ((kappa.Ds>0) || (kappa.DP>0)) # kappa.Ds needed?
  {
    if (D^(1.0/(P-1))>0)
    {
      Ds<-(D)^(-1.0/(P-1))
    } else {return (list (eval=0, DP=0, LoF=0, mse=0, df=df, compound=10^6));}
  }
  if (kappa.DP>0)
  {
    if (df>0)
    {
      DP <- Ds * stats::qf(1-alpha.DP,P-1,df)
    } else {return (list (eval=0, DP=0, LoF=0, mse=0, df=df, compound=10^6));} # if df=0
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
       
      LoF<-(det(L0))^(-1.0/Q) * stats::qf(1-alpha.LoF,Q,df)
    } else {return (list (eval=0, DP=0, LoF=0, mse=0, df=df, compound=10^6));} # if df=0
  }
  if (kappa.mse>0)
  {
    Tvalue<-0
    M12<-t(X1[,-c(1,2)])%*%Z0%*%X2[,-1]
    MM<-t(M12)%*%Minv[-1,-1]%*%M12
    beta2<-rep(tau,Q)                                        # prior point estimate = rep(tau,Q)
    Tvalue<-(1+t(beta2)%*%MM%*%beta2)
    mse<-(Tvalue / D)^(1./(P-1))                         # MSE(D)_s point estimate
  }
  compound<-DP^kappa.DP*LoF^kappa.LoF*mse^kappa.mse
  list (eval=1, DP=DP, LoF=LoF, mse=mse, df=df, compound=compound)
}
