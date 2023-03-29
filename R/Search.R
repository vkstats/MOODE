Search <- function(mood.object, pointex = T)
{
  K <- mood.object$K
  Levels <- mood.object$Levels
  Klev <- mood.object$Klev
  
  P <- mood.object$P
  Q <- mood.object$Q
  primary.terms <- mood.object$primary.terms
  potential.terms <- mood.object$potential.terms
  
  Nruns <- mood.object$Nruns
  Nstarts <- mood.object$Nstarts
  Biter <- mood.object$Biter
  
  W <- mood.object$W; Z0 <- mood.object$Z0
  tau2 <- mood.object$tau2; tau <- mood.object$tau
  
  #cand <- mood.object$cand
  #cand.full <- mood.object$cand.full
  criterion.choice <- mood.object$criterion.choice
  #Parameters <- mood.object$Parameters
  Cubic <- mood.object$Cubic
  orth<- mood.object$orth
  
  kappa.Ls <- mood.object$kappa.Ls
  kappa.LP <- mood.object$kappa.LP
  kappa.Ds <- mood.object$kappa.Ds
  kappa.DP <- mood.object$kappa.DP
  kappa.LoF <- mood.object$kappa.LoF
  kappa.bias <- mood.object$kappa.bias
  kappa.mse <- mood.object$kappa.mse
  
  #kappa.all<-mood.object$kappa.all
  #alpha.all<-mood.object$alpha.all
  
  alpha.DP <- mood.object$alpha.DP
  alpha.LP <- mood.object$alpha.LP
  alpha.LoF <- mood.object$alpha.LoF
  alpha.LoFL <- mood.object$alpha.LoFL
  
  search.object <- list("Nruns" = Nruns, "P" = P, "Q" = Q,
                        "criterion.choice" = criterion.choice,
                        "primary.terms" = primary.terms,
                        "potential.terms" = potential.terms,
                        "Biter" = Biter, "tau2" = tau2,"tau" = tau, 
                        "W" = W, "Z0" = Z0,
                        "alpha.DP" = alpha.DP,"alpha.LP" = alpha.LP,
                        "alpha.LoF" = alpha.LoF, "alpha.LoFL" = alpha.LoFL,
                        "kappa.Ds" = kappa.Ds, "kappa.DP" = kappa.DP,
                        "kappa.Ls" = kappa.Ls, "kappa.LP" = kappa.LP,
                        "kappa.LoF" = kappa.LoF, "kappa.bias" = kappa.bias, 
                        "kappa.mse" = kappa.mse)
  
  start_time <- Sys.time()
  
  cand.trt <- candidate_trt_set(Levels, K, Cubic)   # form the candidate set of treatments
  cand.full <- candidate_set_full(cand.trt, K)      # build candidate set, with potential terms
  
  if (orth){
    cand.full <- candidate_set_orth(cand.full,
                                    primary.terms, potential.terms)
    print("Point exchange algorithm will be  used for othonormalised candidate sets")
  }

  crit.values<-matrix(0,ncol=1, nrow=Nstarts)
  for (k in 1:Nstarts){
    
    initial <- initial.full(cand.full, Nruns, primary.terms, potential.terms)
    
    X1 <- initial$X1
    X2 <- initial$X2
    if (k==1){
      
      crit.opt <- objfun(X1, X2, search.object)$compound
      X1.opt <- X1; X2.opt <- X2
    }
    s <- 1
    while (s==1)
    {
      if(pointex) {
        Xs <- point.swap(X1, X2, cand.full, search.object)
      } else {
        Xs <- coord.swap(X1, X2, K, Levels, search.object)
      }
      X1 <- Xs$X1; X2 <- Xs$X2;
      s <- Xs$search
    }
    crit.values[k] <- Xs$compound # track the change of criterion values
    if (crit.opt>Xs$compound)
    {
      X1.opt <- Xs$X1; X2.opt<-Xs$X2
      crit.opt <- Xs$compound
    }
  }
  finish_time <- Sys.time()
  time <- finish_time-start_time

  criteria.opt <- objfun(X1, X2, search.object)
  list (time=time, X1=X1.opt, X2=X2.opt,  df=criteria.opt$df, compound=criteria.opt$compound,
        path=crit.values,criteria.opt=criteria.opt)
}
