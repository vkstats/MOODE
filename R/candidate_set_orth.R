#' Forms the orthonormalised full candidate set
#'
#' This function forms the full extended orthonormalised candidate set of primary and potential terms, with labels.
#' @param cand Candidate set containing primary terms, with labels in the first column.
#' @return The orthonormalised full candidate set containing primary and potential terms, with labels.
#' @export
#' @examples
#' 
#' # Full extended orthonormalised candidate set for two 3-level factors, 
#' # full quadratic polynomial model as primary model and all three-order terms as potential.
#' 
#' K<-2; Levels <- rep(list(1:3),K);
#' Parameters <- c(1, rep(1,K), rep(1,K), K*(K-1)/2) 
#' candidate_set_orth(candidate_set(Levels)) 

candidate_set_orth<-function(cand.full, primary.terms, potential.terms)
{
  cand.not.orth<-cand.full[, c(primary.terms, potential.terms)]    # extended model matrix, no labels
  cand.full.orth<-cbind(cand.full[,"label"],
                        orthonormalization(cand.full.not.orth, basis=FALSE))   # orthonormalisation, adding labels
  return (cand.full.orth)
}

