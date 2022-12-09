
spheric<-function(candlt, K)          # cand - treatment matrix, with labels
{
  treat<-candlt[,-1]
  cand<-treat[apply(treat^2,1,sum)==0,]
  for (i in 1:K)
  {
    cand<-rbind(cand, sqrt(K/i)*treat[apply(treat^2,1,sum)==i,])
  }
  return (cbind(candlt[,1],cand))
}
