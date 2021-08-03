#' Forms the candidate set
#'
#' This function forms the full extended orthonormalised candidate set, with labels.
#' @param cand Candidate set containing primary terms, with labels in the first column.
#' @return The orthonormalised full candidate set containing primary and potential terms, with labels.
#' @export
#' @examples
#' candidate_set_orth(candidate_set(rep(list(1:3),5))) # Candidate set of 5 factors each of 3 levels

candidate_set_orth<-function(cand)
{
  potential<-potential.matrix(cand)           # potential matrix for extended model, labelled
  cand.full.not.orth<-cbind(cand[,-1],potential[,-1])                                 # extended model matrix, no labels
  cand.full<-cbind(cand[,1],orthonormalization(cand.full.not.orth,basis=FALSE))        # orthonormalisation, adding labels
  return (cand.full)
}

