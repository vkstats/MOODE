### Swapping treatments
swap<-function(X1,X2,cand.full) {
  Xcrit<-criteria(X1,X2,criterion.choice)  #***opt. criterion***
  Xcomp<-Xcrit$compound
  search<-0
  n<-nrow(cand.full)
  for (l in 1:Nruns)
  {
    move<-ifelse(l==1||((l>1)&&(X1[l,1]!=X1[(l-1),1])),1,0)
    if (move == 1){
      Xc1<-X1
      Xc2<-X2
      for (i in 1:n)
      {
        if (X1[l,1]!=cand.full[i,1])  # look at labels
        {
          Xc1[l,]<-cand.full[i,1:(P+1)]
          Xc2[l,]<-cand.full[i,c(1,(P+2):(P+Q+1))]
          Ccrit<-criteria(X1=Xc1,X2=Xc2,criterion.choice) #***opt. criterion***
          Ccomp<-Ccrit$compound
          if (Xcomp>Ccomp)    # if the new design is better (minimising)
          {
            X1<-Xc1; X2<-Xc2
            Xcomp<-Ccomp
            search<-1
          }
        }
      }
    }
  }
  list (X1=X1, X2=X2, compound=Xcomp, search=search, crit=criteria(X1,X2,criterion.choice))
}
