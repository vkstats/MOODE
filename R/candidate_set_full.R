#' Forms the full candidate set
#'
#' This function forms the full extended candidate set by joining primary and potential terms, with labels, not orthonormalised.
#' @param cand Candidate set of primary terms, the first column contains treatment labels. 
#' Usually obtained as output from the "candidate_set" function.
#' @return The full extended candidate set: the column of treatment labels, then P primary terms, then Q potential terms.
#' @export
#' @examples
#' 
#' # Full extended candidate set for two 3-level factors, full quadratic polynomial model
#' 
#' K<-2; Levels <- rep(list(1:3),K);
#' Parameters <- c(1, rep(1,K), rep(1,K), K*(K-1)/2)
#' candidate_set_full(candidate_set(Levels))

candidate_set_full<-function(cand)
{
  cand.full<-cbind(cand, potential.matrix(cand)[,-1])
  return (cand.full)
}
