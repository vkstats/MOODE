#' Forms the candidate set of primary terms
#'
#' This function forms the candidate set of primary terms from the factors' levels, adds labels, with optional transforming the coordinates.
#' @param Levels Levels of each factor.
#' @return The extended and labelled candidate set.
#' @export
#' @examples
#' 
#' # Candidate set for five 3-level factors, full quadratic polynomial model
#' 
#' K<-5; Levels <- rep(list(1:3),K);
#' Parameters <- c(1, rep(1,K), rep(1,K), K*(K-1)/2)
#' candidate_set(Levels) 

candidate_set<-function(Levels)
{
  cand<-as.matrix(expand.grid(Levels))
  candl<-cbind(label(cand),cand)                           # labelling
  candlt<-cbind(candl[,1],apply(as.matrix(candl[,-1]),2,transform))   # rescaling
  if (Cubic=='N') {candlt<-as.matrix(spheric(candlt))}     # spheric coords
  candset<-extend(candlt)                                  # extended matrix
  return (candset)
}
