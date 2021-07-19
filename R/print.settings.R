print.settings <- function(object) {
  cat("Number of factors:", K, "\n", "\n")
  cat("Number of levels:", Klev, "\n", "\n")
  cat("Number of runs of the experiment:", Nruns, "\n", "\n")
  cat("Number of random starts of the search:", Nstarts , "\n", "\n")
  cat("Number of MC iterations for the MSE criteria:", Biter, "\n", "\n")
  cat("Criterion choice: MSE.L or MSE.D or MSE.P:", criterion.choice, "\n", "\n")
  cat("Scaling parameter of the prior variance, squared:", tau2, "\n", "\n")
  cat("Weight for DP:", kappa.DP, "\n", "\n")
  cat("Weight for loF:", kappa.LoF, "\n", "\n")
  cat("Weight for MSE:", kappa.mse, "\n", "\n")
  }
