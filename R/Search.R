Search <- function(mood.object, pointex = T)
{
  K <- mood.object$K
  Levels <- mood.object$Levels
  Klev <- mood.object$Klev
  
  P <- mood.object$P
  Q <- mood.object$Q
  
  Nruns <- mood.object$Nruns
  Nstarts <- mood.object$Nstarts
  
  cand <- mood.object$cand
  cand.full <- mood.object$cand.full
  criterion.choice <- mood.object$criterion.choice
  Parameters <- mood.object$Parameters
  Cubic <- mood.object$Cubic
  orth<- mood.object$orth
  
  kappa.Ls <- mood.object$kappa.Ls
  kappa.LP <- mood.object$kappa.LP
  kappa.Ds <- mood.object$kappa.Ds
  kappa.DP <- mood.object$kappa.DP
  kappa.LoF <- mood.object$kappa.LoF
  kappa.bias <- mood.object$kappa.bias
  
  alpha.DP <- mood.object$alpha.DP
  alpha.LP <- mood.object$alpha.LP
  alpha.LoF <- mood.object$alpha.LoF
  alpha.LoFL <- mood.object$alpha.LoFL
  
  start_time <- Sys.time()
  cand <- candidate_set(Levels, K, Cubic)   # form the candidate set of treatments, primary terms
  
  # build candidate set, with potential term
  if(orth=='Y')
    {cand.full <- candidate_set_orth(cand)}
  if(orth == 'N')
    {cand.full <- candidate_set_full(cand, K)}

  crit.values<-matrix(0,ncol=1, nrow=Nstarts)
  for (k in 1:Nstarts)
  {
    initial <- initial.full(cand.full, P, Q, Nruns)
    X1 <- initial$X1
    X2 <- initial$X2
    if (k==1)
    {
      crit.opt <- objfun(X1, X2, P, Q, 
                       kappa.Ls, kappa.LP, kappa.Ds, kappa.DP, 
                       kappa.LoF, kappa.bias, 
                       Nruns, criterion.choice)$compound
      X1.opt <- X1; X2.opt <- X2
    }
    s <- 1
    while (s==1)
    {
      if(pointex) {
        Xs <- point.swap(X1, X2, P, Q, Nruns, cand.full, kappa.Ls, kappa.LP, kappa.Ds, kappa.DP, 
                       kappa.LoF, kappa.bias, criterion.choice)
      } else {
        Xs <- coord.swap(X1, X2, K, P, Q, Nruns, Levels, kappa.Ls, kappa.LP, kappa.Ds, kappa.DP, 
                         kappa.LoF, kappa.bias, criterion.choice)
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

  criteria.opt <- objfun(X1, X2, P, Q, 
                       kappa.Ls, kappa.LP, kappa.Ds, kappa.DP, 
                       kappa.LoF, kappa.bias, 
                       Nruns, criterion.choice)
  list (time=time, X1=X1.opt, X2=X2.opt,  df=criteria.opt$df, compound=criteria.opt$compound,
        path=crit.values,criteria.opt=criteria.opt)
}
