
### Obtaining a point estimate of the MSE(Ds)-component

MSE.point<-function(S.mse, P, Q, Nruns)
{
  M12<-t(S.mse$X1[, -c(1,2)])%*%Z0%*%(S.mse$X2[,-1])
  M<-t(S.mse$X1[,-1])%*%S.mse$X1[,-1]
  D<-prod(round(eigen(M, symmetric=TRUE, only.values=TRUE)$values, 8))
  Minv<-solve(M)
  MM<-t(M12)%*%Minv[-1,-1]%*%M12
  beta2<-rep(tau, Q)                                        # prior point estimate = rep(tau,Q)
  Tvalue<-(1+t(beta2)%*%MM%*%beta2)
  mse<-(Tvalue*Nruns/D)^(1./(P-1))                         # MSE(D)_s_point
  return (mse)
}
