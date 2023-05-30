#' Forms the orthonormalised full candidate set
#'
#' This function forms the full extended orthonormalised candidate set of primary and potential terms, with an intercept column and labels.
#' @param cand.full Candidate set containing terms up to 4th order, with labels in the first column.
#' @param primary.terms Character vector identifying primary model terms.
#' @param potential.terms Character vector identifying potential model terms.
#' @return The orthonormalised full candidate set containing primary and potential terms, with labels.
#' @export
#' @examples
#' 
#' # Full extended orthonormalised candidate set for two 4-level factors, 
#' # full quadratic polynomial model as primary model and all three-order terms as potential.
#' 
#' K<-2; Levels <- rep(list(1:4),K)
#' cand.trt <- candidate_trt_set(Levels, K)
#' cand.full <- candidate_set_full(cand.trt, K)
#' prime.terms <- colnames(cand.full)[2:7]
#' poten.terms <- colnames(cand.full)[8:11]
#' Parameters <- c(1, rep(1,K), rep(1,K), K*(K-1)/2) 
#' candidate_set_orth(cand.full, prime.terms, poten.terms) 

candidate_set_orth<-function(cand.full, primary.terms, potential.terms)
{
  cand.not.orth<-cand.full[, c(primary.terms, potential.terms)]    # extended model matrix, no labels
  cand.full.orth<-cbind(label = cand.full[,"label"],
                        far::orthonormalization(cand.not.orth, basis=FALSE))   # orthonormalisation, adding labels
  # cand.full.orth[, "intercept"] <- 1 # exclude the intercept from the orthonormalisation 
  return (cand.full.orth)
}

