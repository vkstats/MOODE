#' Evaluates the design based on a specific optimality criterion
#'
#' This function evaluates the design based on a specific optimality criterion.
#' @param criterion.choice Criterion choice can be "GL", "GLP", "GD", "GDP", "MSE.L", "MSE.D", "MSE.P".
#' @return A list...
#' @export


objfun<-function(X1, X2, search.object){
  
  criterion.choice <- search.object$criterion.choice
  if(criterion.choice=="GL"){Xcrit<-criteria.GL(X1, X2, search.object)}
  if(criterion.choice=="GLP"){Xcrit<-criteria.GLP(X1, X2, search.object)}
  if(criterion.choice=="GD") {Xcrit<-criteria.GD(X1, X2, search.object)}
  if(criterion.choice=="GDP") {Xcrit<-criteria.GDP(X1, X2, search.object)}
  if(criterion.choice=="MSE.L") {Xcrit<-criteria.mseL(X1, X2, search.object)}
  if(criterion.choice=="MSE.D") {Xcrit<-criteria.mseD(X1, X2, search.object)}
  if(criterion.choice=="MSE.P") {Xcrit<-criteria.mseP(X1, X2, search.object)}
  
  return(Xcrit)
}
