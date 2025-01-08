#' S3 print method
#' @param x mood object
#' @description
#' Prints a summary of the mood object, including parameters that define the 
#' experiment and the (compound) criterion under which the design 
#' will be sought. 
#' @param ... further arguments passed to or from other methods
#' @return No return value, prints summary of object to output
#' @export
print.settings <- function(x, ...) {
  cat("Number of factors:", x$K, "\n", "\n")
  cat("Number of levels:", x$Klev, "\n", "\n")
  cat("Number of runs of the experiment:", x$Nruns, "\n", "\n")
  cat("Number of random starts of the search:", x$Nstarts , "\n", "\n")
  cat("Number of MC iterations for the MSE criteria:", x$Biter, "\n", "\n")
  cat("Criterion choice: GL or GLP or GD or GDP or MSE.L or MSE.D or MSE.P:", x$criterion.choice, "\n", "\n")
  cat("Scaling parameter of the prior variance, squared:", x$tau2, "\n", "\n")
  cat("Weight for DP:", x$kappa.DP, "\n", "\n")
  cat("Weight for loF:", x$kappa.LoF, "\n", "\n")
  cat("Weight for MSE:", x$kappa.mse, "\n", "\n")
  }
