 # Co-ordinate exchange algorithm 

coord.swap <- function(X1, X2, K, Levels, search.object) { 
  
  Xcrit <- objfun(X1, X2, search.object) 
  Xcomp <- Xcrit$compound
  search <- 0
  
  Nruns <- search.object$Nruns
  primary.terms <- search.object$primary.terms
  potential.terms <- search.object$potential.terms
  
  levels_scaled <- search.object$levels_scaled # list of levels, scaled to [-1,1]
  levels_steps <- search.object$levels_steps   # difference in labels when levels moved to the next
  #steps <- search.object$steps # length of one "step" for each factor
  levels_lengths <- search.object$levels_lengths # number of levels of each factor
  
  all.terms = c(primary.terms, potential.terms)
  # all linear terms, i.e. factor levels to be swapped
  treatments = all.terms[sapply(all.terms, function(x) (nchar(x) < 4))]   
  
  d <- X1[, c("label", treatments)] # swapping just the linear terms
  dc <- d
  Xc1 <- X1
  Xc2 <- X2
  
  for (i in 1:Nruns) {
    if(i==1 || ((i>1) && (X1[i, "label"] != X1[(i-1), "label"]))) {
      for(j in 2:(K+1)) {  
        dc <- d
        dc.ind <- which(levels_scaled[[j-1]] == dc[i,j]) # index of the current level value 
        Xc1 <- X1
        Xc2 <- X2
        
        exchange.ind <- (1:levels_lengths[[j-1]])[!(1:levels_lengths[[j-1]]) %in% dc.ind]
        
        for (k in exchange.ind){ # working through the indexes
          
          dc.ind <- which(levels_scaled[[j-1]] == dc[i,j])
          dc[i,1] <- dc[i, 1] + (k - dc.ind)*levels_steps[[j-1]] # change in the index * label change per one "step"
          dc[i,j] <- levels_scaled[[j-1]][k]

          candij <- candidate_set_full(matrix(dc[i, ], nrow = 1, 
                                              dimnames = list(1, colnames(d))), K)

          Xc1[i, ] <- candij[, c("label", primary.terms)]
          Xc2[i, ] <- candij[, c("label", potential.terms)]
          
          Ccrit <- objfun(X1 = Xc1, X2 = Xc2, search.object) 
          Ccomp <- Ccrit$compound
          if (Xcomp > Ccomp) {   # if the new design is better (minimising)
            d[i, j] <- dc[i, j] 
            d[i, 1] <- dc[i, 1]  # label 
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
  
  