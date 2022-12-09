
design_evalation<-function(X1, X2, P, Q, Nruns)
{
  matrix=matrix(criteria.values(X1, P, Nruns),nrow=1,byrow=T)
  colnames(matrix)<-c("Ds","Ls","DP","LP","df","LoF","compound")
  matrix.mse=matrix(criteria.values.mse(X1, X2, P, Q, Nruns),nrow=1,byrow=T)
  colnames(matrix.mse)<-c("Ds","DP","LoFDP","mseD","Ls","LP","LoFLP","mseL")#,"compoundD","compoundL")
  matrix.G=matrix(criteria.values.G(X1, X2, P, Q, Nruns),nrow=1,byrow=T)
  colnames(matrix.G)<-c("Ds","DP","LoFD","LoFDP","biasD","Ls","LP","LoFL","LoFLP","biasL")
  list(matrix,matrix.mse)
}