Search<-function(object.settings)
{
  
  Levels  <- object.settings$Levels
  X1 <- object.settings$X1
  X2 <- object.settings$X2
  P <- object.settings$P
  Q <- object.settings$Q
  Nruns <- object.settings$Nruns
  Nstarts <- object.settings$Nstarts
  cand <- object.settings$cand
  cand.full <- object.settings$cand.full
  criterion.choice <- object.settings$criterion.choice
  K <- object.settings$K
  Parameters <- object.settings$Parameters
  Cubic <- object.settings$Cubic
  orth<- object.settings$orth
    
  start_time<-Sys.time()
  cand<-candidate_set(Levels, K, Parameters, Cubic)   # form the candidate set of treatments, primary terms

  # build candidate set, with potential term
  if(orth=='Y')
  { cand.full<-candidate_set_orth(cand)}
  if(orth=='N')
  {cand.full<-candidate_set_full(cand, K)}

  crit.values<-matrix(0,ncol=1, nrow=Nstarts)
  for (k in 1:Nstarts)
  {
    initial<-initial.full(cand.full, P, Q)
    X1<-initial$X1
    X2<-initial$X2
    if (k==1)
    {
      crit.opt<-criteria(X1, X2, P, Q, Nruns, criterion.choice)$compound # choose the ***opt. criterion*** to be used
      X1.opt<-X1; X2.opt<-X2
    }
    s<-1
    while (s==1)
    {
      Xs<-swap(X1, X2, P, Q, Nruns, cand.full, criterion.choice)
      X1<-Xs$X1; X2<-Xs$X2;
      s<-Xs$search
    }
    crit.values[k]<-Xs$compound # track the change of criterion values
    if (crit.opt>Xs$compound)
    {
      X1.opt<-Xs$X1; X2.opt<-Xs$X2
      crit.opt<-Xs$compound
    }
  }
  finish_time<-Sys.time()
  time<-finish_time-start_time

  criteria.opt<-criteria(X1.opt, X2.opt, P, Q, 
                         kappa.Ls, kappa.LP, kappa.Ds, kappa.DP, 
                         kappa.LoF, kappa.bias, 
                         Nruns, criterion.choice) #***opt. criterion***
  list (time=time, X1=X1.opt, X2=X2.opt,  df=criteria.opt$df, compound=criteria.opt$compound,
        path=crit.values,criteria.opt=criteria.opt)
}
