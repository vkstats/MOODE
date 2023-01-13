#' Setting up the parameters of the experiment, criteria and search procedure
#'
#' This function 
#' @param K 
#' @return The object
#' @export
#' @examples
#' 
#'  
#' 
#' 



settings <- function(K, Klev, Nruns = 40,
            criterion.choice="MSE.P", kappa.DP=1, kappa.LoF=0, kappa.mse=0,
            kappa.Ds=0, kappa.Ls=0,
            Nstarts = 10,     # Number of random starts of the search
            Biter=50,         # Number of MC iterations for the MSE criteria
            Cubic = "Y",      # 'Y' - cubic, 'N' - spheric coordinates
            tau2 = 1, 
            MC = "Y",           # 'Y' - Yes, 'N' - No
            prob.LP=0.95,  prob.LoFL=0.95,
            alpha.LoF=0.05,  alpha.DP=0.05,
            orth='N'){
  
  Levels<-rep(list(1:Klev), K)         # Levels of each factor
  primary.terms <- c("x1", "x2", "x3", "x12", "x22", "x32", "x1x2", "x1x3", "x2x3")       # primary model terms
  potential.terms <- c("x12x2", "x1x22", "x12x3", "x1x32", "x22x3", "x2x32", "x1x2x3")    # potential terms
  
  P <- length(primary.terms)+1    # number of the primary parameters, including intercept
  Q <- length(potential.terms)    # number of the potential terms 
  Z0<-diag(1, Nruns)-matrix(1/Nruns, nrow=Nruns, ncol=Nruns)   # for excluding the intercept
  tau<-tau2^0.5 

  # Weights on parameters for Ls/LP-criteria. 
  # Default: quadratic terms are assigned 1/4 of the weights on the rest of the parameters  
  W<-matrix(1, nrow=1, ncol = P-1)  
  colnames(W) <- primary.terms
  W[c("x12", "x22", "x32")] = 1./4
  W <- matrix(W/sum(W), nrow=1)
  
  # Confidence levels
  prob.LP<-ifelse(MC=='Y', prob.LP^(1./(P-1)), prob.LP)            # Correction for multiple comparisons             
  prob.LoFL<-ifelse(MC=='Y', prob.LoFL^(1./Q), prob.LoFL)
  alpha.LP<-1-prob.LP; alpha.LoFL<-1-prob.LoFL; 
  # MSE(L/D) weights (kappa.(L/D)P, kappa.LoF, kappa.mse), sum of the weights = 1
  kappa.Ls<-0; kappa.Ds<-0  # weights set to zero
  
  # create some output
  out <- list("K" = K, "Klev" = Klev, "Nruns" = Nruns, "criterion.choice" = criterion.choice,
              "Nstarts" = Nstarts, "Biter" = Biter, "tau2" = tau2,"tau" = tau, "Levels" = Levels, "Cubic" = Cubic,
              "MC" = MC,"prob.LP" = prob.LP,"prob.LoFL" = prob.LoFL,
              "alpha.LP" = alpha.LP,"alpha.LoF" = alpha.LoF, "alpha.DP" = alpha.DP,
               "orth" = orth, "Z0" = Z0,"W" = W, "primary terms" = primary.terms, 
              "potential terms" = potential.terms,
              "kappa.DP" = kappa.DP,"kappa.LoF" = kappa.LoF,"kappa.mse" = kappa.mse, #weights
              "P" = P, "Q" = Q, "kappa.Ls" = kappa.Ls, "kappa.Ds" = kappa.Ds)
  class(out) <- append(class(out), "settings")
#  return(attach(out, warn.conflicts = F)) ## remember, this also basically does "print(out)"
 return(out) 
}

 

