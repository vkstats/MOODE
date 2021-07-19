
output.mse<-function(S.criterion)
{
  list1<-c(S.criterion$df, S.criterion$compound, S.criterion$mse, S.criterion$time)
  list(out=list1, values=criteria.values.mse(S.criterion$X1, S.criterion$X2), mse.point=MSE.point(S.criterion))
}
