
### Output function
output<-function(S.criterion)
{
  list1<-c(S.criterion$df, S.criterion$compound, S.criterion$time)
  list(out=list1, values=criteria.values(S.criterion$X1, P, Nruns))
}
