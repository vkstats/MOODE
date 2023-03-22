 # Co-ordinate exchange algorithm 

coord.swap <- function(X1, X2, K, P, Q, Nruns, Levels, kappa.Ls, kappa.LP, kappa.Ds, kappa.DP, 
                       kappa.LoF, kappa.bias, criterion.choice) { #orth = T
  
  Xcrit <- objfun(X1, X2, P, Q, kappa.Ls, kappa.LP, kappa.Ds, kappa.DP, 
                  kappa.LoF, kappa.bias, Nruns, criterion.choice) 
  Xcomp <- Xcrit$compound
  search <- 0
  
  levels_scaled <- sapply(Levels, Transform) # columns are factors, rows are levels
  
  d <- X1[, 1:(K + 1)]
  dc <- d
  Xc1 <- X1
  Xc2 <- X2
  
  for(i in 1:Nruns) {
    if(i==1 || ((i>1) && (X1[i, 1] != X1[(i-1), 1]))) {
      for(j in 2:(K+1)) {
        dc <- d
        Xc1 <- X1
        Xc2 <- X2
      
        for(k in levels_scaled[!levels_scaled[, j - 1] == dc[i, j], j - 1]) {
          dc[i, j] <- k
          candij <- candidate_set_full(matrix(dc[i, ], nrow = 1, dimnames = list(1, colnames(d))), K)
          
          # if(orth){ 
          #   candij <- candidate_set_orth(matrix(dc[i, ], nrow = 1, dimnames = list(1, colnames(d))))
          # } else{          candij <- candidate_set_full(matrix(dc[i, ], nrow = 1, dimnames = list(1, colnames(d))), K)
          # }
          Xc1[i, ] <- candij[, 1:(P + 1)]
          Xc2[i, ] <- candij[, c(1, (P + 2):(P + Q + 1))]
          Ccrit <- objfun(X1 = Xc1, X2 = Xc2, P, Q, kappa.Ls, kappa.LP, kappa.Ds, kappa.DP, 
                        kappa.LoF, kappa.bias,  Nruns, criterion.choice) 
          Ccomp <- Ccrit$compound
          if (Xcomp > Ccomp) {   # if the new design is better (minimising)
            d[i, j] <- dc[i, j] 
            X1[i, ] <- Xc1[i, ]
            X2[i, ] <- Xc2[i, ]
            Xcomp <- Ccomp
            search <- 1
          }
        }
      }
    }
  }
 list(X1 = X1, X2 = X2, compound=Xcomp, search=search, 
     crit=objfun(X1, X2, P, Q, kappa.Ls, kappa.LP, kappa.Ds, kappa.DP, kappa.LoF, kappa.bias, Nruns, criterion.choice))
}
  
  