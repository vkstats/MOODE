
spheric<-function(candlt, K)          # cand - treatment matrix, with labels
{
  treat <- candlt[, -1]
  rs <- rowSums(treat^2)
  cand <- matrix(NA, nrow = K, ncol = ncol(treat))
  cand[1, ] <- treat[rs == 0, ]
  for (i in 1:K)
  {
    cand[i + 1, ] <- rbind(cand, sqrt(K / i) * treat[rs == i, ])
  }
  return (cbind(candlt[, 1], cand))
}
