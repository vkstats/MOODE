#' Forms the labelled candidate set of treatments
#'
#' This function forms the candidate set of treatments from the factors' levels, 
#' adds labels, with optional spherical transformation of the coordinates.
#' @param Levels Levels of each factor.
#' @param K Number of factors.
#' @param Cubic Indicator whether the experimental region is cubic ('Y') or spherical ('N').
#' @return Matrix of candidate set of treatments, with treatment labels contained in the first column.
#' @export
#' @examples
#'
#' # Candidate treatment set for five 3-level factors
#' 
#' K<-5; Levels <- rep(list(1:3),K);
#' candidate_set(Levels, K) 

candidate_set <- function(Levels, K, Cubic = TRUE)
{
  cand <- as.matrix(expand.grid(Levels))
  candl <- cbind(label(cand, Levels, K), cand)                         # creating a column of treatment labels
  cand.trt <- cbind(candl[,1], apply(as.matrix(candl[,-1]), 2, Transform))   # rescaling the factors' values to [-1, 1]
      if (!Cubic) {cand.trt<-as.matrix(spheric(cand.trt, K))}     # spheric coords
  return (cand.trt)
}
