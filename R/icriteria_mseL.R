#' Internal function to find criteria values for vignette example
#' 
#' @param X1 model matrix for primary terms
#' @param X2 model matrix for the potential terms
#' @param search.object object of class mood specifying the experiment details and model
#' 
#' @noRd

icriteria.mseL <- function(X1, X2, search.object) {
  crit_values <- NULL
  ## LP
  search.object$kappa.LP <- 1; search.object$kappa.LoF <- 0; search.object$kappa.mse <- 0
  crit_values[1] <- criteria.mseL(X1, X2, search.object)$LP
  ## LoF
  search.object$kappa.LP <- 0; search.object$kappa.LoF <- 1; search.object$kappa.mse <- 0
  crit_values[2] <- criteria.mseL(X1, X2, search.object)$LoF
  ## MSE
  search.object$kappa.LP <- 0; search.object$kappa.LoF <- 0; search.object$kappa.mse <- 1
  crit_values[3] <- criteria.mseL(X1, X2, search.object)$mse
  ## DF
  search.object$kappa.LP <- 1; search.object$kappa.LoF <- 1; search.object$kappa.mse <- 1
  crit_values[4] <- criteria.mseL(X1, X2, search.object)$df
  ## Ls
  search.object$kappa.LP <- 0; search.object$kappa.LoF <- 0; search.object$kappa.mse <- 0; search.object$kappa.Ls <- 1
  search.object$criterion.choice <- "GL"
  crit_values[5] <- criteria.GL(X1, X2, search.object)$Ls
  list(LP = crit_values[1], LoF = crit_values[2], mse = crit_values[3], df = crit_values[4], Ls = crit_values[5])
}