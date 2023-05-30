#' Setting up the parameters of a factorial experiment to search for multi-objective optimal completely randomised design.
#' @description Creates an object containing the parameters of the experiment, compound optimality criterion with the
#'  weights and parameters of the search.
#' 
#' @param K Number of factors.
#' @param Levels Either (a) a common number of levels for each factor or (b) a list of length K of the vectors containing levels of each factor.
#' @param Nruns Number of runs of the experiment.
#' @param criterion.choice Compound criterion to be used for the optimal design search or evaluation. 
#' Possible values are: 
#' \itemize{
#' \item "GL", "GD" for Generalised D- and L-optimality (Goos et al., 2005) 
#' \item "GDP" and "GLP" for Generalised DP- and LP-optimality
#' \item "MSE.D", "MSE.L" and "MSE.P" for compound criteria with MSE-component: determinant-based, trace-based and determinant-based but with point estimates for the MSE(D)-component
#'}
#' @param kappa.Ds Weight of the Ds-criterion; value must be between 0 and 1.
#' @param kappa.DP Weight of the DP-criterion; value must be between 0 and 1.
#' @param kappa.Ls Weight of the Ls-criterion; value must be between 0 and 1.
#' @param kappa.LP Weight of the LP-criterion; value must be between 0 and 1.
#' @param kappa.LoF Weight of the lack-of-fit component criterion; value must be between 0 and 1.
#' @param kappa.bias Weight of the bias component criterion; value must be between 0 and 1.
#' @param kappa.mse Weight of the MSE component criterion; value must be between 0 and 1.
#' @param Nstarts The number of randomly generated start designs of the search algorithm.
#' @param Cubic Indicator of whether the experimental region is cubic (`TRUE`) or spherical (`FALSE`).
#' @param tau2 The variance scaling parameter for the prior distribution of the potential terms.
#' @param Biter Number of samples for evaluating the MSE determinant-based component criterion.
#' @param MC Indicator of the multiple comparison (Bonferroni) correction for trace-based criteria.
#' @param prob.DP Confidence level for the DP-criterion.  
#' @param prob.LP Confidence level for the LP-criterion; pre-Bonferroni correction.
#' @param prob.LoF Confidence level for the Lack-of-fit criterion.
#' @param primary.model The order of the fitted polynomial model. Alternatively polynomial terms can be directly specified through the `primary.terms` parameter, see Details.
#' @param potential.model The order of the potential/extra polynomial terms. Alternatively can be specified through the `potential.terms` parameter, see Details.
#' @param primary.terms Vector of the names of the primary terms, see Details.
#' @param potential.terms Vector of the names of the potential terms, see Details.
#' @param orth Indicator of whether to orthonormalise the potential and primary terms (`TRUE') or not (`FALSE').
#' 
#' @export
#' @details The function provides different ways of specifying the levels of the factors and the models.
#' 
#' Specifying the factors and levels
#' 
#' If all `K` factors have the same number of levels, `Levels` parameter is used to input that number.
#' Otherwise, `Levels` is set to be a list of vectors containing the values of the factors, e.g.
#' `list(1:3, 1:2, 1:4)` for 3 factors with equally spaced \eqn{3, 2} and \eqn{4} levels respectively.
#' 
#' Specifying the fitted model and the potential terms
#' 
#' There are two ways to describe the primary and potential sets of model terms.
#' `primary.model` and `potential.model` arguments can be used to specify the fitted model and the potential terms in a string from. 
#' They are used to generate the sets of `primary.terms` and `potential.terms` which can be input directly. 
#' Possible values of the `primary.model` argument are: 
#' \itemize{
#' \item `main_effects` -- main effects for all the factors
#' \item `first_order` -- main effects and linear interactions
#' \item `second_order` -- full second order polynomial 
#' \item `third_order` -- full second order polynomial model and all interactions of degree 3
#' \item `cubic` -- third order polynomial model with cubic terms 
#' }
#' The intercept is always included as a primary term. 
#'
#' Possible elements of the `potential.model` vector argument:
#' \itemize{ 
#' \item `linear_interactions` -- linear interactions among the factors
#' \item `quadratic_terms` -- quadratic terms for all the factors
#' \item `third_order_terms` --  all interactions of degree 3: linear-by-linear-by-linear and quadratic terms
#' \item `cubic_terms` -- cubic terms for all the factors 
#' \item `fourth_order_terms` -- all interactions of degree 4, similar to `third_order_terms`
#' }
#' `primary.terms` and `potential.terms` arguments are used to specify the fitted model and the potential terms explicitly -- up to the total power of 4.
#' \itemize{
#' \item Single factor powers,  are coded as a string starting with with "x" and followed by 
#' the index of the factor and the power: `"x32"`. 
#' For example, \eqn{$x_3^2$} is coded as `"x32"`; `"x22"` stands for \eqn{x_2^2}, and `"x4"` stands for the linear term \eqn{x_4}.
#' \item Interaction of factors' powers are coded by merging the individual factors' 
#' powers, ordered by the factors' indexes. For example, \eqn{x_2^2\times x_1} is coded as `"x1x22"`, 
#' \eqn{x_3x_12x_4} -- as `"x12x3x4"`.
#' }
#' 
#' 
#' @return List of parameters of the experiment, compound criterion of choice, and primary and potential model terms.
#' \itemize{
#' \item `K` Number of factors.
#' \item `Klev` Number of levels of each factor, if all factors have the same number of levels.
#' \item `Levels` List of length K of the vectors containing values of the factors.
#' \item `Nruns` Number of runs of the experiment.
#' \item `criterion.choice` Compound criterion to be used for the optimal design search or evaluation.
#' \item `Nstarts` The number of randomly generated start designs of the search algorithm.
#' \item `Biter` Number of samples for evaluating the MSE determinant-based component criterion.
#' \item `tau2` The variance scaling parameter for the prior distribution of the potential terms.
#' \item `tau`  The square root of `tau2`
#' \item `Cubic` Whether the experimental region is cubic (`TRUE`) or spherical (`FALSE`).
#' \item `MC` Indicator of the multiple comparison (Bonferroni) correction for trace-based criteria.
#' \item `prob.DP, prob.LP, prob.LoF, prob.LoFL` Confidence levels for the DP-, LP-, lack of fit determinant- and trace-based criteria.
#' \item `alpha.DP, alpha.LP, alpha.LoF, alpha.LoFL` Significance levels for the DP-, LP-, lack of fit determinant- and trace-based criteria.
#' \item `orth` Whether the candidate sets are orthonormalised (`TRUE`) or not (`FALSE`).
#' \item `Z0` Z0 matrix.
#' \item `W` Weight matrix for Ls criterion.
#' \item `primary.terms` Fitted (primary) model terms.
#' \item `potential.terms` Potential terms.
#' \item `P` The number of terms in the fitted model (including intercept).
#' \item `Q` The number of potential terms.
#' \item `kappa.Ds, kappa.DP, kappa.Ls, kappa.LP, 
#' kappa.LoF, kappa.bias, kappa.mse` Compound criterion weights.
#' \item `warning.msg` Warning messages.
#' }
#' @examples
#' 
#'example1 <- mood(K = 5, Levels = 3, Nruns = 40, criterion.choice = "GDP", 
#'kappa.Ds = 1./3, kappa.DP = 1./3, kappa.LoF = 1./3, 
#'Nstarts = 50, tau2 = 0.1, primary.model = "second_order",
#' potential.model = NA, potential.terms = c("x12x2", "x22x3", "x32x4", "x42x5"))
#'example1
#'
#'example2 <- mood(K = 3, Levels = list(1:3, 1:2, 1:2), criterion.choice = "MSE.L",
#'kappa.LP = 1./2, kappa.LoF = 1./4, kappa.mse = 1./4,
#'Nstarts = 50, tau2 = 1, primary.terms = "first_order",
#' potential.model = NA, potential.terms = c("x12", "x12x2", "x12x3"))
#' example2
#' 
mood <- function(K, 
                 Levels, 
                 Nruns = 15, 
                 criterion.choice = "MSE.P", 
                 kappa.Ds = 0.0, kappa.DP = 1.0, 
                 kappa.Ls = 0.0, kappa.LP = 0.0,
                 kappa.LoF = 0.0, kappa.bias = 0.0, kappa.mse = 0.0,  
                 Nstarts = 10,     # Number of random starts of the search
                 Cubic = TRUE,     # Cubic experimental region (T) or spherical (F)
                 tau2 = 1, 
                 Biter=50,         # Number of MC iterations for the MSE criteria
                 MC = TRUE,        # Bonferroni correction for multiple comparisons
                 
                 prob.DP = 0.95, prob.LP = 0.95, prob.LoF = 0.95, prob.LoFL = 0.95,
                 
                 primary.model = "first_order",
                 potential.model = NA, 
                 
                 primary.terms = NA, potential.terms = NA,
                 
                 orth=FALSE){      # Orthonormalised terms 
  
  warning.msg <- c()
  
# need some checks on inputs
  Klev = NA
  if (identical(length(Levels), as.integer(1))) {
    Klev = Levels
    Levels <- rep(list(1:Klev), K)   # Levels of each factor if all factors have the same number of levels
  } 
  
  if (any(is.na(primary.terms))){   # if primary terms are not specified explicitly
    primary.terms <- c()
    if (any(primary.model %in% c("main_effects","first_order", "second_order",
                                 "third order", "cubic"))) {
      for (k in 1:K){
        primary.terms <- c(primary.terms, paste("x", as.character(k), sep = ""))
      }
    }
    
    if (any(primary.model %in% c("second_order", "third_order", "cubic"))){
      for (k in 1:K){      # quadratic terms
        primary.terms <- c(primary.terms, paste("x", as.character(k), as.character(2), sep = ""))
      }
    }
    
    if (any(primary.model %in% c("first_order", "second_order", 
                                 "third_order", "cubic"))){
      if(K > 1) { 
        for (k in 1:(K-1)){  # interaction terms
          for (j in (k+1):K) {
            primary.terms <- c(primary.terms, 
                              paste("x", as.character(k), "x", as.character(j), sep = ""))
          }
        }
      }
    }
    if (any(primary.model %in% c("third_order", "cubic"))){
      if (K >1){
        for (k in 1:(K-1)){
          for (j in (k+1):K) {  # QxL  and LxQ interaction terms: "x12x2", "x1x22", etc.
            primary.terms <- c(primary.terms, 
                              paste("x", as.character(k), as.character(2), "x", as.character(j), sep = ""),
                              paste("x", as.character(k), "x", as.character(j), as.character(2), sep = "")) 
          }
        }
      } 
      
      if (K > 2) {
        for (k in 1:(K-2)){
          for (j in (k+1):(K-1)) {
            for (i in (j+1):K){        # LxLxL terms: "x1x2x3", "x1x3x4", etc.
              primary.terms <- c(primary.terms, 
                                paste("x", as.character(k), "x", as.character(j), "x", as.character(i), sep = ""))   
            }
          }
        }
      }
    }
    if ("cubic" %in% primary.model){
      for (k in 1:K){     # cubic terms: "x13", "x33", etc.
        primary.terms <- c(primary.terms, 
                          paste("x", as.character(k) , as.character(3), sep = ""))   
      }
    }
  }
    
  # primary terms always contains the intercept
  primary.terms<- c("intercept", primary.terms)
  
  if (any(is.na(potential.terms))) {  # if potential terms are not specified explicitly
    potential.terms <- c()
    if ("linear_interactions" %in% potential.model) {
      if (K>1){
        for (k in 1:(K-1)){  # interaction terms
          for (j in (k+1):K) {
            potential.terms <- c(potential.terms, 
                                paste("x", as.character(k), "x", as.character(j), sep = ""))
          }
        }
      } else {
        warning.msg <- append(warning.msg, 
                              "Warning: no linear interactions can be added, too few factors.")
      }
      
    }
    if ("quadratic_terms" %in% potential.model){
      for (k in 1:K){      # quadratic terms
        potential.terms <- c(potential.terms, paste("x", as.character(k), as.character(2), sep = ""))
      }
    }
    
    if ("third_order_terms" %in% potential.model){
      if (K>1){
        for (k in 1:(K-1)){
          for (j in (k+1):K) {  # QxL  and LxQ interaction terms: "x12x2", "x1x22", etc.
            potential.terms <- c(potential.terms, 
                                paste("x", as.character(k), as.character(2), "x", as.character(j), sep = ""),
                                paste("x", as.character(k), "x", as.character(j), as.character(2), sep = "")) 
          }
        }
      } else {
        warning.msg <- append(warning.msg, 
                              "Warning: no third order terms can be added, too few factors.")
        }
      
      if (K > 2) {
        for (k in 1:(K-2)){
          for (j in (k+1):(K-1)) {
            for (i in (j+1):K){        # LxLxL terms: "x1x2x3", "x1x3x4", etc.
              potential.terms <- c(potential.terms, 
                                  paste("x", as.character(k), "x", as.character(j), "x", as.character(i), sep = ""))  
            }
          }
        }
      }
    }
    if ("cubic_terms" %in% potential.model){
      for (k in 1:K){     # cubic terms: "x13", "x33", etc.
        potential.terms <- c(potential.terms, 
                            paste("x", as.character(k) , as.character(3), sep = ""))   
      }
    }
    if ("fourth_order_terms" %in% potential.model){
      if (K>1){
        for (k in 1:(K-1)){
          for (j in (k+1):K) { 
            # Quadratic x Quadratic
            potential.terms <- c(potential.terms, 
                                paste("x", as.character(k), as.character(2), 
                                      "x", as.character(j), as.character(2), sep = ""))
            # Cubic x Linear
            potential.terms <- c(potential.terms, 
                                paste("x", as.character(k), as.character(3), "x", as.character(j), sep = ""),
                                paste("x", as.character(k), "x", as.character(j), as.character(3), sep = "")) 
          }
        }
      } else {
        warning.msg <- append(warning.msg,
                              "Warning: no fourth order terms can be added, too few factors.")
        }
      
      if (K > 2) {
        for (k in 1:(K-2)){
          for (j in (k+1):(K-1)) {
            for (i in (j+1):K){        # QxLxL, LxQxL and LxLxQ terms: "x12x2x3", "x1x22x4", etc.
              potential.terms <- c(potential.terms, 
                                  paste("x", as.character(k), as.character(2), "x", as.character(j), "x", as.character(i), sep = ""),
                                  paste("x", as.character(k), "x", as.character(j), as.character(2), "x", as.character(i), sep = ""),
                                  paste("x", as.character(k), "x", as.character(j), "x", as.character(i), as.character(2), sep = ""))
            }
          }
        }
      }
      if (K > 3) {
        for (k in 1:(K-3)){
          for (j in (k+1):(K-2)) {
            for (i in (j+1):(K-1)){
              for (l in (i+1):K){   # LxLxLxL terms: "x1x2x3x4", etc.
                potential.terms <- c(potential.terms, 
                                    paste("x", as.character(k), "x", as.character(j),
                                          "x", as.character(i), "x", as.character(l),sep = ""))
              }
            }
          }
        }
      }
    }
  }
  
  P <- length(primary.terms) + 1    # number of the primary parameters, including intercept
  Q <- length(potential.terms)      # number of the potential terms 
  
  Z0 <- diag(1, Nruns)-matrix(1/Nruns, nrow=Nruns, ncol=Nruns)   # for excluding the intercept
  tau<-tau2^0.5 
  
  # Weights on parameters for Ls/LP-criteria. 
  # Default: quadratic terms are assigned 1/4 of the weights on the rest of the parameters  
  W<-matrix(1, nrow=1, ncol = P-1)  
  colnames(W) <- primary.terms
  quadratic.names <- c()
  for (k in 1:K){      # quadratic terms
    quadratic.names <- c(quadratic.names, paste("x", as.character(k), as.character(2), sep = ""))
  }
  W[quadratic.names] <- 1./4
  W <- matrix(W/sum(W), nrow=1)
  
  # Confidence levels
  prob.LP<-ifelse(MC, prob.LP^(1./(P-1)), prob.LP)            # Corrections for multiple comparisons             
  prob.LoFL<-ifelse(MC, prob.LoFL^(1./Q), prob.LoFL)
  
  alpha.LP<-1-prob.LP; alpha.LoFL<-1-prob.LoFL
  alpha.DP <- 1-prob.DP; alpha.LoF<-1-prob.LoF
  
  # Criterion choice check
  
 
  
  if (!(criterion.choice %in% c("GD", "GL", "GDP", "GLP", "MSE.D", "MSE.L", "MSE.P"))) {
    print ("Error: invalid criterion choice. Please choose any of the following: 
           GD, GL, GDP, GLP, MSE.D, MSE.L, MSE.P")
    return()
  }
 
 
   
  # Kappa-s check
  if ((any(is.na(potential.terms)) || is.null(potential.terms) || (length(potential.terms) == 0)) && 
      (any(is.na(potential.model)) || is.null(potential.model) || (length(potential.model) == 0))){
    kappa.LoF <- 0; kappa.mse <- 0; kappa.bias <- 0
    warning.msg <- append(warning.msg,  "No potential terms have been specified. 
                                        Corresponding weights have been set to 0. 
                                        Please check the other weights.")
  }
  
  kappa.all <- c(kappa.Ds, kappa.Ls, kappa.DP, kappa.LP, 
                 kappa.LoF, kappa.bias, kappa.mse)
  if ((any(kappa.all < 0)) || (any(kappa.all>1))) {
    print("Error: criteria weights should be between 0 and 1")
    return()
  }
  
  if ((criterion.choice == "GD") && ((kappa.Ds + kappa.LoF + kappa.bias) != 1.0)) {
    print("Error: the sum of criteria weights kappa.Ds, kappa.LoF and kappa.bias should be equal to 1")
    return()
  }
  if ((criterion.choice == "GL") && ((kappa.Ls + kappa.LoF + kappa.bias) != 1.0)) {
    print("Error: the sum of criteria weights kappa.Ls, kappa.LoF and kappa.bias should be equal to 1")
    return()
  }
  if ((criterion.choice == "GDP") && ((kappa.Ds + kappa.DP +kappa.LoF + kappa.bias) != 1.0)) {
    print("Error: the sum of criteria weights kappa.Ds, kappa.DP, kappa.LoF and kappa.bias should be equal to 1")
    return()
  }
  if ((criterion.choice == "GLP") && ((kappa.Ls + kappa.LP +kappa.LoF + kappa.bias) != 1.0)) {
    print("Error: the sum of criteria weights kappa.Ls, kappa.LP, kappa.LoF and kappa.bias should be equal to 1")
    return()
  }
  if ((criterion.choice %in% c("MSE.D", "MSE.P")) && ((kappa.DP + kappa.LoF + kappa.mse) != 1.0)) {
    print("Error: the sum of criteria weights kappa.DP, kappa.LoF and kappa.mse should be equal to 1")
    return()
  }
  
  
  if ((criterion.choice == "MSE.L") && ((kappa.LP + kappa.LoF + kappa.mse) != 1.0)) {
    print("Error: the sum of criteria weights kappa.LP, kappa.LoF and kappa.mse should be equal to 1")
    return()
  }

  
    
  cat(warning.msg, sep = "\n") # print out warning messages
  
  
  # create the output
  out <- list("K" = K,
              "Klev" = Klev,
              "Levels" = Levels, "Nruns" = Nruns, 
              "criterion.choice" = criterion.choice,
              "Nstarts" = Nstarts, "Biter" = Biter, "tau2" = tau2,
              "tau" = tau, "Levels" = Levels, "Cubic" = Cubic,
              "MC" = MC,
              "prob.DP"=prob.DP, "prob.LP" = prob.LP,
              "prob.LoF" = prob.LoF, "prob.LoFL" = prob.LoFL,
              "alpha.DP" = alpha.DP,"alpha.LP" = alpha.LP,
              "alpha.LoF" = alpha.LoF, "alpha.LoFL" = alpha.LoFL,
              "orth" = orth, "Z0" = Z0, "W" = W, 
              "primary.terms" = primary.terms, "potential.terms" = potential.terms,
              "P" = P, "Q" = Q,
              "kappa.Ds" = kappa.Ds, "kappa.DP" = kappa.DP,
              "kappa.Ls" = kappa.Ls, "kappa.LP" = kappa.LP,
              "kappa.LoF" = kappa.LoF, "kappa.bias" = kappa.bias, "kappa.mse" = kappa.mse,
              "warning messages" = warning.msg)
  
  class(out) <- append(class(out), "settings")
  #  return(attach(out, warn.conflicts = F)) ## remember, this also basically does "print(out)"
  return(out) 
}

