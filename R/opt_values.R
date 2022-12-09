
opt.values<-function(design)
{ #MSE
  Nruns=nrow(design) # Number of runs of the experiment
  Z0<-diag(1, Nruns)-matrix(1/Nruns, nrow=Nruns, ncol=Nruns)   # for excluding the intercept

  #Obtain the model matrix
  # treatment matrix (with no labels) -> extended design matrix (with labels)
  if(dim(table(design[,1]))<=Klev){design=cbind(rep(1, nrow(design)), design)}
  X1=extend(design)
  #Obtain the potential matrix (no labels)
  X2=potential.matrix(X1)

  #Optimal values
  optimal.values=design_evalation(X1, X2, P, Q, Nruns)

  return(optimal.values)
}
