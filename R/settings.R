#' Setting up the parameters of the experiment, criteria and search procedure.
#'
#' This function 
#' @param K Number of factors.
#' @param Levels List of length K of the vectors containing values of the factors.
#' @param Klev If all factors have the same number of levels: number of levels of each factor.
#' @return The object
#' @export
#' @examples
#' 
#'  
#' 
#' 



settings <- function(K, Levels, Klev,
                     Nruns = 40, 
                     criterion.choice="MSE.P", 
                     kappa.Ds = 0.0, kappa.DP = 1.0, kappa.Ls = 0.0, kappa.LP = 0.0,
                     kappa.LoF = 0.0, kappa.mse = 0.0,  
                     Nstarts = 10,     # Number of random starts of the search
                     Cubic = "Y",      # 'Y' - cubic, 'N' - spheric coordinates
                     tau2 = 1, 
                     Biter=50,         # Number of MC iterations for the MSE criteria
                     MC = "Y",         # Whether MC is used in MSE 'Y' - Yes, 'N' - No
                     
                     prob.DP = 0.95, prob.LP = 0.95,  prob.LoF = 0.95,
                     alpha.DP = 0.05, alpha.LP = 0.05, alpha.LoF= 0.05,  
                     
                     primary.model = "first order",
                     potential.model = "second order",
                     
                     primary.terms = NA, potential.terms = NA,
                     
                     orth='N'){
  
  
  
  if (length(Klev) == 1) {
    Levels<-rep(list(1:Klev), K)   # Levels of each factor if all factors have the same number of levels
  } 
  
  if (is.na(primary.terms)) {  # if primary terms are not specified explicitly
    primary.terms = c()
    if (primary.model %in% c("first order", "second order")) {
      for (k in 1:K){
        primary.terms = c(primary.terms, paste("x", as.character(k), sep = ""))
      }
    }
    
    if (primary.model == c("second order")){
      for (k in 1:K){      # quadratic terms
        primary.terms = c(primary.terms, paste("x", as.character(k), as.character(2), sep = ""))
      }
      for (k in 1:(K-1)){  # interaction terms
        for (j in (k+1):K) {
          primary.terms = c(primary.terms, 
                            paste("x", as.character(k), "x", as.character(j), sep = ""))
        }
      }
    }
  }
  
  if (is.na(potential.terms)) {
    potential.terms = c()
    if ("linear interactions" %in% potential.model) {   #
      for (k in 1:(K-1)){  # interaction terms
        for (j in (k+1):K) {
          primary.terms = c(primary.terms, 
                            paste("x", as.character(k), "x", as.character(j), sep = ""))
        }
      }
    }
    if ("second order" %in% potential.model){
      for (k in 1:K){      # quadratic terms
        potetial.terms = c(primary.terms, paste("x", as.character(k), as.character(2), sep = ""))
      }
    }
    if ("third order" %in% potential.model){
      for (k in 1:(K-1)){
        for (j in (k+1):K) {  # QxL  and LxQ interaction terms: "x12x2", "x1x22", etc.
          potential.terms = c(potential.terms, 
                              paste("x", as.character(k), as.character(2), "x", as.character(j), sep = ""),
                              paste("x", as.character(k), "x", as.character(j), as.character(2), sep = "")) 
        }
      }
      if (K > 2) {
        for (k in 1:(K-2)){
          for (j in (k+1):(K-1)) {
            for (i in (j+1):K){        # LxLxL terms: "x1x2x3", "x1x3x4", etc.
              potential.terms = c(potential.terms, 
                                  paste("x", as.character(k), "x", as.character(j), "x", as.character(i), sep = "")   
            }
          }
        }
      }
      for (k in 1:K){     # cubic terms: "x13", "x33", etc.
        potential.terms = c(potential.terms, 
                            paste("x", as.character(k) , as.character(3), sep = "")   
      }
      
    }
  }
  
  #primary.terms <- c("x1", "x2", "x3", "x12", "x22", "x32", "x1x2", "x1x3", "x2x3")       # primary model terms
  #potential.terms <- c("x12x2", "x1x22", "x12x3", "x1x32", "x22x3", "x2x32", "x1x2x3")    # potential terms
  
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

 

