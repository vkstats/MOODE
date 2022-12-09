#' Calculates the values of the Generalised LPs-criterion and its components
#'
#' This function evaluates the Generalised LPs-criterion for given primary and potential model matrices. 
#' Components: Ls-, LPs-, LoF(LP)- and Bias(L)-optimality.
#' The weights kappa.Ls, kappa.LP, kappa.LoF and kappa.bias are taken from the global environment.
#'
#' @param X1 The primary model matrix, with the first column containing the labels of treatments, and the second -- the intercept term.
#' @param X2 The matrix of potential terms, with the first column containing the labels of treatments.
#' @param eps Computational tolerance, the default value is 10^-23
#'
#' @return A list of values: indicator of whether the evaluation was successful ("eval"), Ls-criterion value -- intercept excluded ("Ls"),
#' LPs-criterion value -- intercept excluded ("LPs"), Lack-of-fit(LP) criterion value ("LoF"), the bias component value ("bias"), 
#' the number of pure error degrees of freedom ("df") and the value of the compound criterion ("compound").
#' @export
#' @examples 
#'#Experiment: one 5-level factor, primary model -- full quadratic, X^3 and X^4 potential terms.
#'K <-1; P<-3; Q<-2; Levels <- list(1:5)
#' # Generating candidate sets: primary and full orthonormalised ones
#'cand.primary <- candidate_set(Levels);
#'cand.not.orth <-cbind(cand.primary[,-1], cand.primary[,3]^3, cand.primary[,4]^2)
#'cand.full.orth <- cbind(cand.primary[,1], far::orthonormalization(cand.not.orth,basis=FALSE))
#' # Choosing a design
#'index <- c(rep(1,2),3,rep(4,2),rep(5,3)); Nruns<- length(index)
#'X.primary <- cand.full.orth[index, 1:(P+1)]
#'X.potential <- cand.full.orth[index, (c(1,(P+2):(P+Q+1)))]
#' # Evaluating a compound GLP-criterion
#'kappa.Ls = kappa.LP = kappa.LoF = kappa.bias = 0.25; tau2 <-1;
#'alpha.LP = alpha.LoFL = 0.05;
#'criteria.GLP(X1 = X.primary, X2 = X.potential)
#' 
#'Output: eval = 1, Ls = .5315, LP = 4.0969, LoF = 5.3727, bias = 1.4013, df = 4, compound = 2.0122
#'
criteria.GLP<-function(X1, X2, P, Q, Nruns, eps=10^-23)      # X1, X2 -- matrices of primary and potential terms, both with labels
{
  Ls<-0; LP<-0; LoF<-0; bias<-0;
  DF<-nlevels(as.factor(X1[,1]))
  df<-Nruns-DF                                # df - pure error degrees of freedom

  M<-crossprod(X1[,-1])                       # information matrix of primary terms
  D<-prod(round(eigen(M, symmetric=TRUE, only.values=TRUE)$values,8))/Nruns
  if (D>eps)
  {
    Minv<-solve(M)
  } else {return (list (eval=0,Ls=0, LP=0, LoF=0, bias=0, df=df, compound=10^6));}

  if ((kappa.Ls>0)||(kappa.LP>0))
  {
    Ls<-W%*%(diag(Minv)[-1])
  }
  if (kappa.LP>0)
  {
    if (df>0)
    {
      LP<-Ls*qf(1-alpha.LP,1,df)
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
      LoF<-L0.inv.trace*qf(1-alpha.LoFL,1,df)/Q
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

