#' Forms the matrix of the potential terms
#'
#' This function forms the matrix of the potential (polynomial) terms. 
#' 
#' @param cand Candidate set of primary terms, the first column contains treatment labels. 
#' Usually obtained as output from the "candidate_set" function.
#' @param pmat Potential terms can be reduced cubic, cubic, quartic
#' @return 
#' @export
#' @examples

potential.matrix<-function(cand, K)   # cand -- matrix with labels
{
  cand.X2<-matrix(0,nrow=nrow(cand),ncol=1)
  cand.X3<-matrix(0,nrow=nrow(cand),ncol=1)
  cand.X4<-matrix(0,nrow=nrow(cand),ncol=1)

  #Model:full quadratic
  #Note that the order is: 1 (labels), 2 (intercept),
  # 3:K+2 (linear terms), K+3:2K+2 (quadratic terms)
  #

  nL=K+2
  nQ1=K+3
  nQ2=2*K+2

  n_2fi=2*K+3
  n_fi=0
  N_2fi=length(n_2fi:ncol(cand))

  Xl= cand[,3:nL]                                     #linear terms #head(Xl)
  Xq= cand[,nQ1:nQ2]                                  #quadratic terms (Xq)

  for (i in 1:K){                                   # all L*Q terms
    for (j in 1:K){
      if(j!=i){
        cand.X2<-cbind(cand.X2,Xl[,i]*Xq[,j])
      }
    }
  }

  for (i in 1:K){                                   # all L*L*L terms
    X_2fi=cand[,n_2fi:ncol(cand)]#2fiead(Xq)
    n_fi=n_fi+K-i
    if(n_fi<N_2fi){
      X_2fi=as.matrix(X_2fi[,-c(1:n_fi)])
      for (j in 1:ncol(X_2fi)){
        cand.X3<-cbind(cand.X3,Xl[,i]*X_2fi[,j])
      }
    }
  }

  X2<-cbind(cand[,1],cand.X2[,-1],cand.X3[,-1])        # adding column with labels


  for (i in 1:K)                               # All cubic terms
  {
    cand.X4<-cbind(cand.X4,Xl[,i]^3)
  }

  X3<-cbind(cand[,1],cand.X2[,-1],cand.X3[,-1],cand.X4[,-1])        # adding column with labels

  return (X3) # return X2 is the reduced cubic, return X3 is the full cubic

  # if(pmat=="cubic") {
  #   return(X3)
  # } else {
  #   return(X2)
  # }

}

#quartic
# we need more:
# all L*C terms        -- K(K-1)
# all QxLxL terms      -- 3K(K-1)(K-2)/6
# all quartic terms    -- K
# in total 70 terms

# cand=candidate_set(Levels)
#
# cand.X5<-matrix(0,nrow=nrow(cand),ncol=1)
# for (i in 1:K){                                    # all L*C terms
#   for (j in 1:K){
#     if(j!=i){
#       cand.X5<-cbind(cand.X5,Xl[,i]*Xl[,j]^3)
#     }
#   }
# }
# dim(cand.X5)
# cand.X6<-matrix(0,nrow=nrow(cand),ncol=1)
# for (i in 1:K){                                # all L*L*Q terms
#   for (j in 1:K){
#       if(j>i){
#          for (k in 1:K){
#            if(k!=i & k!=j){
#                cand.X6<-cbind(cand.X6,Xl[,i]*Xl[,j]*Xl[,k]^2)
#         }}}}}
# dim(cand.X6)
#
# cand.X7<-matrix(0,nrow=nrow(cand),ncol=1)
# for (i in 1:K)                               # All cubic terms
# {
#   cand.X7<-cbind(cand.X7,Xl[,i]^4)
# }
#
# dim(cand.X7)
