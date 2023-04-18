# Co-ordinate exchange algorithm, only for all same-levels factors

coord.swap.same <- function(X1, X2, K, Levels, search.object) { 
  
  Xcrit <- objfun(X1, X2, search.object) 
  Xcomp <- Xcrit$compound
  search <- 0
  
  # columns are factors, rows are levels
  levels_scaled <- sapply(Levels, Transform) 
  
  # list of factors' levels scaled to [-1,1]
  #levels_scaled <- lapply(Levels, Transform) 
  
  Nruns <- search.object$Nruns
  primary.terms <- search.object$primary.terms
  potential.terms <- search.object$potential.terms
  
  all.terms = c(primary.terms, potential.terms)
  # all linear terms, i.e. factor levels to be swapped
  treatments = all.terms[sapply(all.terms, function(x) (nchar(x) == 2))]   
  
  d <- X1[, c("label", treatments)] # swapping just the linear terms
  dc <- d
  Xc1 <- X1
  Xc2 <- X2
  
  for (i in 1:Nruns) {
    if(i==1 || ((i>1) && (X1[i, "label"] != X1[(i-1), "label"]))) {
      for(j in 2:(K+1)) {  
        dc <- d
        Xc1 <- X1
        Xc2 <- X2
        
        for (k in levels_scaled[!levels_scaled[, j - 1] == dc[i, j], j - 1]) {
        #for (k in levels_scaled[[j-1]][!levels_scaled[[j - 1]] == dc[i, j]]) {
          dc[i, j] <- k
          
          candij <- candidate_set_full(matrix(dc[i, ], nrow = 1, 
                                              dimnames = list(1, colnames(d))), K)
          
          Xc1[i, ] <- candij[, c("label", primary.terms)]
          Xc2[i, ] <- candij[, c("label", potential.terms)]
          Ccrit <- objfun(X1 = Xc1, X2 = Xc2, search.object) 
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
       crit=objfun(X1, X2, search.object))
}

