print.settings <- function(object) {
  cat("Number of factors:", object$K, "\n", "\n")
  cat("Number of levels:", object$Klev, "\n", "\n")
  cat("Number of runs of the experiment:", object$Nruns, "\n", "\n")
  cat("Number of random starts of the search:", object$Nstarts , "\n", "\n")
  cat("Number of MC iterations for the MSE criteria:", object$Biter, "\n", "\n")
  cat("Criterion choice: MSE.L or MSE.D or MSE.P:", object$criterion.choice, "\n", "\n")
  cat("Scaling parameter of the prior variance, squared:", object$tau2, "\n", "\n")
  cat("Weight for DP:", object$kappa.DP, "\n", "\n")
  cat("Weight for loF:", object$kappa.LoF, "\n", "\n")
  cat("Weight for MSE:", object$kappa.mse, "\n", "\n")
  }
