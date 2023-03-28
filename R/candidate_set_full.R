#' Forms the full candidate set of treatments polynomial terms
#'
#' This function forms the full extended candidate set with all polynomial terms, with labels, not orthonormalised.
#' @param cand Candidate set of treatments, the first column contains treatment labels. 
#' Usually obtained as output from the "candidate_set" function.
#' @return The full extended candidate set: the column of treatment labels, then named columns with polynomial terms up to the order of 3. 
#' For example, "x12" stands for $x_1^2$, and "x1x2 stands for $x_1x_2$, and "x23x4" -- for $x_2^3x_4$.
#' @export
#' @examples
#' 
#' # Full extended candidate set for two 3-level factors
#' 
#' K<-2; Levels <- rep(list(1:3),K);
#' candidate_set_full(candidate_set(Levels, K), K)

candidate_set_full = function(cand, K) {
  
  cand.terms = cand[, -1, drop = F]  #  linear terms

    ### Naming linear and generating and naming quadratic terms
  for (k in 1:K) {
    colnames(cand.terms)[k] = paste("x", as.character(k), sep = "")  # named linear terms, e.g. "x3"
    cand.terms = cbind(cand.terms, cand.terms[,k]^2)
    colnames(cand.terms)[K+k] = paste("x", as.character(k), as.character(2), sep = "") # quadratic terms: "x12", "x22", "x32", etc.
  }
  
  ### Interaction terms: linear by linear
  if (K > 1) {
    int.count = 0
    for (k in 1:(K-1)){
      for (j in (k+1):K) {
        int.count = int.count + 1                  # count the number of interaction terms
        cand.terms = cbind(cand.terms, cand.terms[,k]* cand.terms[,j])
        colnames(cand.terms)[2*K + int.count] = 
          paste("x", as.character(k), "x", as.character(j), sep = "")   # interaction terms: "x1x2", "x2x3", "x3x4", etc.
      }
    }
    
    ### Interaction terms: quadratic by linear 
    int2.count = 0
    nterms = ncol(cand.terms)
    for (k in 1:(K-1)){
      for (j in (k+1):K) {
        int2.count = int2.count + 2                # count the number of QxL interaction terms
        cand.terms = cbind(cand.terms, cand.terms[,k]^2 * cand.terms[,j], cand.terms[,k] * cand.terms[,j]^2)
        colnames(cand.terms)[nterms + int2.count - c(1,0)] = 
          c(paste("x", as.character(k), as.character(2), "x", as.character(j), sep = ""),
            paste("x", as.character(k), "x", as.character(j), as.character(2), sep = "")) # QxL  and LxQ interaction terms: "x12x2", "x1x22", etc.
      }
    }
  }
  
  ### Linear-by-linear-by-linear terms
  if (K > 2) {
    int3.count = 0
    nterms = ncol(cand.terms)
    for (k in 1:(K-2)){
      for (j in (k+1):(K-1)) {
        for (i in (j+1):K){
          int3.count = int3.count + 1                     # count the number of LxLxL interaction terms
          cand.terms = cbind(cand.terms, cand.terms[,k] * cand.terms[,j]* cand.terms[,i])
          colnames(cand.terms)[nterms + int3.count] = 
            paste("x", as.character(k), "x", as.character(j), "x", as.character(i), sep = "")   # LxLxL terms: "x1x2x3", "x1x3x4", etc.
        }
        
      }
    }
  }
  
  ### Cubic terms
  nterms = ncol(cand.terms)
  for (k in 1:K){
    cand.terms = cbind(cand.terms, cand.terms[,k]^3)
    colnames(cand.terms)[nterms + k] = paste("x", as.character(k) , as.character(3), sep = "")   # cubic terms: "x13", "x33", etc.
  }
  
  cand.terms = cbind(cand[,1], cand.terms)    # append the column of labels
  colnames(cand.terms)[1] <- "label"
  
  return (cand.terms)
}
