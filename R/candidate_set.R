#' Forms the candidate set of primary terms
#'
#' This function forms the candidate set of primary terms from the factors' levels, adds labels, with optional transforming the coordinates.
#' @param Levels Levels of each factor.
#' @return The extended and labelled candidate set.
#' @export
#' @examples
#' candidate_set(rep(list(1:3),5)) # Candidate set of 5 factors each of 3 levels

candidate_set<-function(Levels)
{
  cand<-as.matrix(expand.grid(Levels))
  candl<-cbind(label(cand),cand)                           # labelling
  candlt<-cbind(candl[,1],apply(candl[,-1],2,transform))   # rescaling
  if (Cubic=='N') {candlt<-as.matrix(spheric(candlt))}     # spheric coords
  candset<-extend(candlt)                                  # extended matrix
  return (candset)
}
