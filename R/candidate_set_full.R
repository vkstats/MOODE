#' Forms the full extended candidate set
#'
#' This function forms the full extended candidate set, primary terms, with labels, not orthonormalised.
#' @param cand Candidate set.
#' @return The full extended candidate set.
#' @export
#' @examples
#' candidate_set_full(candidate_set(rep(list(1:3),5))) # Full extended candidate set of 5 factors each of 3 levels

#cand -  candidate set, primary terms (full quadratic model),with labels, not orthonormalised
candidate_set_full<-function(cand)
{
  cand.full<-cbind(cand,potential.matrix(cand)[,-1])
  return (cand.full)
}
