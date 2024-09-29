### Swapping treatments
#' Swapping points between the current design and candidate set 
#' @description Performing point-exchange algorithm, extensive swap of points procedure between the current design 
#' and candidate set.
#' 
#' @param X1 Current fitted (primary) model matrix
#' @param X2 Current potential terms matrix
#' @param cand.full Full candidate matrix
#' @param search.object Object for the search
#' 
#' @details \code{point.swap} is called within the \code{Search} function
#' 
#' 
#' @return 
#' @examples
#' 
point.swap<-function(X1, X2, cand.full, search.object) {
  
  Xcrit<-objfun(X1, X2, search.object)  
  Xcomp<-Xcrit$compound
  search<-0
  n<-nrow(cand.full)
  
  Nruns<-search.object$Nruns
  P<-search.object$P
  Q<-search.object$Q
  primary.terms<-search.object$primary.terms
  potential.terms<-search.object$potential.terms
  
  for (l in 1:Nruns)
  {
    move<-ifelse(l==1||((l>1)&&(X1[l,"label"]!=X1[(l-1),"label"])),1,0)
    if (move == 1){
      Xc1<-X1
      Xc2<-X2
      for (i in 1:n)
      {
        if (X1[l,1]!=cand.full[i,"label"])  # look at labels
        {
          Xc1[l,]<-cand.full[i, c("label", primary.terms)]
          Xc2[l,]<-cand.full[i, c("label", potential.terms)]
          
          Ccrit<-objfun(X1=Xc1, X2=Xc2, search.object)
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
  list (X1=X1, X2=X2, compound=Xcomp, search=search, 
        crit=objfun(X1, X2, search.object))
}
